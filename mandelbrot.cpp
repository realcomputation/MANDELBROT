#include <iRRAM.h>
#include "iRRAMx/polynomial.hpp"

#include <fstream>
#include <set>

using namespace iRRAM;

constexpr int multiplier = 1024;

class painter {
    int real_resolution, imag_resolution;
    REAL real_lo, real_hi;
    REAL imag_lo, imag_hi;
    REAL distance2;
    bool computed;

    std::set<std::pair<int, int>> unknown;
    std::set<std::pair<std::pair<int, int>, INTEGER>> exterior;
    std::set<std::pair<std::pair<int, int>, INTEGER>> interior;
    std::set<std::pair<int, int>> boundary;

    COMPLEX int_pair_to_COMPLEX(std::pair<int, int> const &p) const {
        return { real_lo + (real_hi - real_lo) / (real_resolution << 1) + (real_hi - real_lo) / real_resolution * p.first,
                 imag_lo + (imag_hi - imag_lo) / (imag_resolution << 1) + (imag_hi - imag_lo) / imag_resolution * p.second };
    }

    COMPLEX INTEGER_pair_to_COMPLEX(std::pair<INTEGER, INTEGER> const &p, INTEGER const &resolution) const {
        return { real_lo + (real_hi - real_lo) / (resolution << 1) + (real_hi - real_lo) / resolution * p.first,
                 imag_lo + (imag_hi - imag_lo) / (resolution << 1) + (imag_hi - imag_lo) / resolution * p.second };
    }

public:
    painter(int real_resolution, int imag_resolution, REAL real_lo, REAL real_hi, REAL imag_lo, REAL imag_hi)
    : real_resolution(real_resolution), imag_resolution(imag_resolution),
      real_lo(std::move(real_lo)), real_hi(std::move(real_hi)),
      imag_lo(std::move(imag_lo)), imag_hi(std::move(imag_hi)),
      distance2(maximum((real_hi - real_lo) / real_resolution, (imag_hi - imag_lo) / imag_resolution)),
      computed(false) {
        distance2 = distance2 * distance2;
        for(int i = 0; i < real_resolution; i++) {
            for(int j = 0; j < imag_resolution; j++) {
                unknown.insert({ i, j });
            }
        }
    }

    void compute() {
        INTEGER iteration = 0;
        INTEGER iteration_power2 = 1;
        while(!unknown.empty()) {
            cout << "iteration #" << iteration << "\n";

            // EXTERIOR
            // std::cout << "EXTERIOR\n";
            for(auto it = unknown.begin(); it != unknown.end();) {
                auto const c = int_pair_to_COMPLEX(*it);
                auto x = c;
                bool flag = false;
                for(INTEGER iter = 0; iter < iteration * multiplier; iter = iter + 1) {
                    x = x * x + c;
                }
                for(int iter = 0; iter < multiplier; iter++) {
                    auto const d = real(x) * real(x) + imag(x) * imag(x);
                    if(choose(d > 4, d < 5) == 1) {
                        exterior.insert({ std::move(*it), -(iteration * multiplier + iter + 1) });
                        it = unknown.erase(it);
                        flag = true;
                        break;
                    }
                    x = x * x + c;
                }
                if(!flag) ++it;
            }

            // INTERIOR
            // std::cout << "INTERIOR\n";
            for(auto it = unknown.begin(); it != unknown.end();) {
                auto const c = int_pair_to_COMPLEX(*it);
                POLYNOMIAL const p(2, { c, REAL(0), REAL(1) });
                POLYNOMIAL P = p;
                //std::cout << "Working on: " << it->first << ' ' << it->second << '\n';
                //cout << x << "\n";
                for(INTEGER i = 0; i < iteration; i = i + 1) {
                    P = composite(P, p);
                }

                P.coef[1] = REAL(-1);
                auto const v = roots(P);
                P.coef[1] = REAL(0);
                auto const D = deriv(P, 1);
                bool flag = false;

                if(std::any_of(v.cbegin(), v.cend(), [&D](const COMPLEX& c) {
                    auto const val = abs(D(c));
                    return choose(val < 1, val > 0.5) == 1;
                })) {
                    interior.insert({ std::move(*it), iteration + 1 });
                    it = unknown.erase(it);
                    flag = true;
                    break;
                }
                if(!flag) ++it;
            }

            // HELPER GRID
            // std::cout << "HELPER GRID\n";
            std::set<std::pair<INTEGER, INTEGER>> grid_exterior;
            std::set<std::pair<INTEGER, INTEGER>> grid_interior;
            for(INTEGER i = 1; i < iteration_power2; i = i + 2) {
                for(INTEGER j = 1; j < iteration_power2; j = j + 2) {
                    auto const c = INTEGER_pair_to_COMPLEX({ i, j }, iteration_power2);
                    auto x = c;
                    bool flag = false;

                    // EXTERIOR
                    for(INTEGER i = 0; i < iteration * multiplier; i = i + 1) {
                        auto const d = real(x) * real(x) + imag(x) * imag(x);
                        if(choose(d > 4, d < 5) == 1) {
                            grid_exterior.insert({ i, j });
                            flag = true;
                            break;
                        }
                        x = x * x + c;
                    }
                    if(flag) break;

                    // INTERIOR
                    POLYNOMIAL const p(2, { c, REAL(0), REAL(1) });
                    POLYNOMIAL P = p;
                    for(INTEGER i = 0; i < iteration; i = i + 1) {
                        P.coef[1] = REAL(-1);
                        auto const v = roots(P);
                        P.coef[1] = REAL(0);
                        auto const D = deriv(P, 1);

                        if(std::any_of(v.cbegin(), v.cend(), [&D](const COMPLEX& c) {
                            auto const val = abs(D(c));
                            return choose(val < 1, val > 0.5) == 1;
                        })) {
                            grid_interior.insert({ i, j });
                            break;
                        }

                        if(i != iteration - 1) P = composite(P, p);
                    }
                }
            }

            // BOUNDARY
            // std::cout << "BOUNDARY\n";
            for(auto it = unknown.begin(); it != unknown.end();) {
                bool flag = false;
                auto const point = int_pair_to_COMPLEX(*it);
                for(auto extr = grid_exterior.begin(); extr != grid_exterior.end(); ++extr) {
                    auto const ext_point = INTEGER_pair_to_COMPLEX(*extr, iteration_power2);
                    auto const ext_diff = ext_point - point;
                    auto const ext_dist2 = real(ext_diff) * real(ext_diff) + imag(ext_diff) * imag(ext_diff);
                    if(choose(ext_dist2 < distance2, ext_dist2 > distance2 / 2) != 1) continue;

                    for(auto intr = grid_interior.begin(); intr != grid_interior.end(); ++intr) {
                        auto const int_point = INTEGER_pair_to_COMPLEX(*intr, iteration_power2);
                        auto const int_diff = int_point - point;
                        auto const int_dist2 = real(int_diff) * real(int_diff) + imag(int_diff) * imag(int_diff);
                        if(choose(int_dist2 < distance2, int_dist2 > distance2 / 2) != 1) continue;

                        boundary.insert(std::move(*it));
                        it = unknown.erase(it);
                        flag = true;
                        break;
                    }

                    if(flag) break;
                }
                if(!flag) ++it;
            }

            iteration = iteration + 1;
            iteration_power2 = (iteration_power2 << 1);
        }
        computed = true;
    }

    void draw(std::ostream &os) {
        if(!computed) compute();
        auto canvas = new INTEGER*[imag_resolution];
        for(int j = 0; j < imag_resolution; j++) {
            canvas[j] = new INTEGER[real_resolution];
            for(int i = 0; i < real_resolution; i++) {
                canvas[j][i] = 0;
            }
        }
        for(auto const &now : exterior) {
            canvas[now.first.second][now.first.first] = now.second;
        }
        for(auto const &now : interior) {
            canvas[now.first.second][now.first.first] = now.second;
        }
        for(auto const &now : boundary) {
            canvas[now.second][now.first] = 0;
        }

        for(int i = 0; i < imag_resolution; i++) {
            for(int j = 0; j < real_resolution; j++) {
                os << std::setfill(' ') << std::setw(8) << (int)canvas[i][j];
            } os << '\n';
        }

        for(int i = 0; i < imag_resolution; i++) {
            delete[] canvas[i];
        } delete[] canvas;
    }
};

constexpr int resolution = 17;

void compute() {
    painter p(resolution, resolution, -2, 2, -2, 2);
    std::ofstream ofs("result.txt");
    p.draw(ofs);
    ofs.close();
}