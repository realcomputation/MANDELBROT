#include <iRRAM/lib.h>
#include <iRRAMx/polynomial.hpp>

#include <ctime>
#include <memory>
#include <set>

void print_current(int resolution, std::unique_ptr<std::unique_ptr<std::unique_ptr<iRRAM::INTEGER>[]>[]> const &canvas) {
    for(int j = 0; j < resolution; j++) {
        for(int i = 0; i < resolution; i++) {
            if(!canvas[i][j]) { printf("?"); continue; }
            if(*canvas[i][j] < 0) { printf(" "); continue; }
            if(*canvas[i][j] > 0) { printf("#"); continue; }
            printf("~");
        } printf("\n");
    }
}

inline iRRAM::REAL iterate_function(iRRAM::INTEGER const &n) {
    auto const p = power(2, size(n) - 1);
    return ((n - p) * 2 + 1) / (p * 2);
}

iRRAM::COMPLEX int_pair_to_COMPLEX(int resolution, std::pair<int, int> const &p, iRRAM::INTEGER const &iteration) {
    iRRAM::INTEGER l = 1, r = 1;
    while(r * (r + 1) / 2 <= iteration) r = (r << 1);
    while(l < r) {
        auto const m = (l + r) / 2;
        if(m * (m + 1) / 2 >= iteration) r = m;
        else l = m + 1;
    }
    l = 1;
    auto const d = r * (r + 1) / 2 - iteration;
    l += d, r -= d;

    return {
            -2 + iRRAM::REAL(4) * p.first / resolution + iRRAM::REAL(4) / resolution * iterate_function(l),
        -2 + iRRAM::REAL(4) * p.second / resolution + iRRAM::REAL(4) / resolution * iterate_function(r)
    };
}
iRRAM::COMPLEX INTEGER_pair_to_COMPLEX(int resolution, std::pair<INTEGER, INTEGER> const &p) {
    return {
            -2 + iRRAM::REAL(2) / resolution + iRRAM::REAL(4) * p.first / resolution,
            -2 + iRRAM::REAL(2) / resolution + iRRAM::REAL(4) * p.second / resolution
    };
}

uint64_t reiteration_count;

void compute(int const &resolution, std::unique_ptr<std::unique_ptr<std::unique_ptr<iRRAM::INTEGER>[]>[]> const &canvas) {
    for(int i = 0; i < resolution; i++) {
        for(int j = 0; j < resolution; j++) {
            canvas[i][j] = nullptr;
        }
    }
    bool loop = true;
    iRRAM::REAL const distance(REAL(4) / resolution);
    iRRAM::REAL const distance2(distance * distance);
    iRRAM::INTEGER iteration(0), iteration_power2(2);
    iRRAM::INTEGER iteration_exterior_next(2048), iteration_interior_next(2);
    std::unique_ptr<std::unique_ptr<std::unique_ptr<iRRAM::COMPLEX>[]>[]> storage_exterior;
    std::unique_ptr<std::unique_ptr<std::unique_ptr<iRRAM::POLYNOMIAL>[]>[]> storage_interior;

    /// INITIALIZATION
    storage_exterior = std::make_unique<std::unique_ptr<std::unique_ptr<iRRAM::COMPLEX>[]>[]>(resolution);
    storage_interior = std::make_unique<std::unique_ptr<std::unique_ptr<iRRAM::POLYNOMIAL>[]>[]>(resolution);
    for(int i = 0; i < resolution; i++) {
        storage_exterior[i] = std::make_unique<std::unique_ptr<iRRAM::COMPLEX>[]>(resolution);
        storage_interior[i] = std::make_unique<std::unique_ptr<iRRAM::POLYNOMIAL>[]>(resolution);
    }

//printf("BEGIN: reiteration = %lu\n", reiteration_count++);
//    clock_t start, end;

    while(loop) {
        loop = false;
//printf("iteration_iteration_next = %s\niteration_interior = %s\n",
//       swrite(iteration_exterior_next).c_str(), swrite(iteration_interior_next).c_str());

        for(int i = 0; i < resolution; i++) {
            for(int j = 0; j < resolution; j++) {
                if(canvas[i][j]) continue;
//printf("working on (%d, %d)...\n", i, j);
                loop = true;
                auto const c = int_pair_to_COMPLEX(resolution, { i, j }, iteration);
                bool flag = false;

//printf("### EXTERIOR ###\n");
                /// EXTERIOR
//                start = clock();
                auto z = iRRAM::COMPLEX(iRRAM::REAL(0));
                for(iRRAM::INTEGER iter = 1; iter < iteration_exterior_next; iter = iter + 1) {
                    z = z * z + c;
                    auto const d2 = real(z) * real(z) + imag(z) * imag(z);
                    if(choose(d2 > 4, d2 < 5) == 1) {
                        canvas[i][j] = std::make_unique<INTEGER>(-iter);
                        storage_exterior[i][j] = nullptr;
                        flag = true;
                        break;
                    }
                }
//                end = clock();
//printf("Elapsed time: %Lf\n", ((long double)(end - start)) / CLOCKS_PER_SEC);
//                print_current(resolution, canvas);
                if(flag) continue;

//printf("### INTERIOR ###\n");
                /// INTERIOR
//                start = clock();
                auto const P = iRRAM::POLYNOMIAL(2, { c, REAL(0), REAL(1) });
                auto p = iRRAM::POLYNOMIAL(1, { REAL(0), REAL(1) });
                for(auto iter = 1; iter < iteration_interior_next; iter = iter + 1) {
                    p = iRRAM::composite(p, P);
                    p.coef[1] = REAL(-1);
                    auto const rootv = iRRAM::roots(p);
                    p.coef[1] = REAL(0);
                    auto const D = iRRAM::deriv(p, 1);
                    for(auto const &r : rootv) {
                        auto const x = D(r);
                        auto const d2 = real(x) * real(x) + imag(x) * imag(x);
                        if(choose(d2 < 1, d2 > 0.5) == 1) {
                            canvas[i][j] = std::make_unique<INTEGER>(iter);
                            storage_exterior[i][j] = nullptr;
                            flag = true;
                            break;
                        }
                    }
                }
//                end = clock();
//printf("Elapsed time: %Lf\n", ((long double)(end - start)) / CLOCKS_PER_SEC);
//                print_current(resolution, canvas);
                if(flag) continue;

//printf("### BOUNDARY ###\n");
                /// HELPER GRID & BOUNDARY
//                start = clock();
                bool grid_exterior = false;
                bool grid_interior = false;
                for(iRRAM::INTEGER x = 1; x < iteration_power2; x = x + 2) {
                    for(iRRAM::INTEGER y = 1; y < iteration_power2; y = y + 2) {
                        auto cur = INTEGER_pair_to_COMPLEX(resolution, { x, y });
                        cur = cur - c;
                        iRRAM::REAL d2 = real(cur) * real(cur) + imag(cur) * imag(cur);

                        /// Check if point on grid is in exterior
                        if(choose(d2 < distance2, d2 > distance2 / 2) == 1) {
                            bool inner_flag = false;
                            z = REAL(0);
                            for(iRRAM::INTEGER iter = 1; iter < iteration_exterior_next; iter = iter + 1) {
                                z = z * z + c;
                                d2 = real(z) * real(z) + imag(z) * imag(z);
                                if(choose(d2 > 4, d2 < 5) == 1) {
                                    grid_exterior = true;
                                    inner_flag = true;
                                    break;
                                }
                            }
                            if(inner_flag) {
                                if(grid_interior) {
                                    canvas[i][j] = std::make_unique<iRRAM::INTEGER>(0);
                                    storage_exterior[i][j] = nullptr;
                                    storage_interior[i][j] = nullptr;
                                    goto END;
                                }
                                else continue;
                            }

                            /// Check if point on grid is in interior
                            auto const P = iRRAM::POLYNOMIAL(2, { cur, REAL(0), REAL(1) });
                            p = iRRAM::POLYNOMIAL(1, { REAL(0), REAL(1) });
                            for(iRRAM::INTEGER iter = 1; iter < iteration_interior_next; iter = iter + 1) {
                                p = composite(p, P);
                                p.coef[1] = REAL(-1);
                                auto const rootv = iRRAM::roots(p);
                                p.coef[1] = REAL(0);
                                auto const D = iRRAM::deriv(p, 1);
                                for(auto const &r : rootv) {
                                    auto const tmp = D(r);
                                    auto const tmp2 = real(tmp) * real(tmp) + imag(tmp) * imag(tmp);
                                    if(choose(tmp2 < 1, tmp2 > 0.5) == 1) {
                                        grid_interior = true;
                                        inner_flag = true;
                                        break;
                                    }
                                }
                            }
                            if(inner_flag) {
                                if(grid_exterior) {
                                    canvas[i][j] = std::make_unique<iRRAM::INTEGER>(0);
                                    storage_exterior[i][j] = nullptr;
                                    storage_interior[i][j] = nullptr;
                                    goto END;
                                }
                                else continue;
                            }
                        }
                    }
                }
//                end = clock();
//printf("Elapsed time: %Lf\n", ((long double)(end - start)) / CLOCKS_PER_SEC);
//                print_current(resolution, canvas);
            END:;
            }
        }

        iteration = iteration + 1;
        //iteration_exterior = iteration_exterior_next;
        iteration_exterior_next = (iteration_exterior_next << 1);
        //iteration_interior = iteration_interior_next;
        iteration_interior_next = (iteration_interior_next + 1);
        iteration_power2 = (iteration_power2 << 1);
    }
}

int main(int argc, char *argv[]) {
    int resolution;
    reiteration_count = 0;
    iRRAM_initialize(argc, argv);

    std::cout << "Input resolution: ";
    std::cin >> resolution;

    clock_t start, end;
    auto canvas = std::make_unique<std::unique_ptr<std::unique_ptr<iRRAM::INTEGER>[]>[]>(resolution);
    for(int i = 0; i < resolution; i++) {
        canvas[i] = std::make_unique<std::unique_ptr<iRRAM::INTEGER>[]>(resolution);
    }

    start = clock();
    iRRAM::exec(compute, resolution, canvas);
    end = clock();
    for(int j = 0; j < resolution; j++) {
		for(int i = 0; i < resolution; i++) {
			std::cout << swrite(*canvas[i][j]) << ' ';
		} std::cout << '\n';
	}
    std::cout << "Elapsed time: " << ((long double)(end - start)) / CLOCKS_PER_SEC << '\n';
}

