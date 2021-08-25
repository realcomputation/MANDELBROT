#include <iRRAM/lib.h>
#include <iRRAMx/polynomial.hpp>

#include <memory>
#include <set>
#include <ctime>

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

iRRAM::COMPLEX int_pair_to_COMPLEX(int resolution, std::pair<int, int> const &p) {
    return {
        -2 + iRRAM::REAL(2) / resolution + iRRAM::REAL(4) * p.first / resolution,
        -2 + iRRAM::REAL(2) / resolution + iRRAM::REAL(4) * p.second / resolution
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
    iRRAM::INTEGER iteration_power2(2);
    iRRAM::INTEGER iteration_exterior(1), iteration_exterior_next(2048), iteration_interior(1), iteration_interior_next(2);
    std::unique_ptr<std::unique_ptr<std::unique_ptr<iRRAM::COMPLEX>[]>[]> storage_exterior;
    std::unique_ptr<std::unique_ptr<std::unique_ptr<iRRAM::POLYNOMIAL>[]>[]> storage_interior;

    /// INITIALIZATION
    storage_exterior = std::make_unique<std::unique_ptr<std::unique_ptr<iRRAM::COMPLEX>[]>[]>(resolution);
    storage_interior = std::make_unique<std::unique_ptr<std::unique_ptr<iRRAM::POLYNOMIAL>[]>[]>(resolution);
    for(int i = 0; i < resolution; i++) {
        storage_exterior[i] = std::make_unique<std::unique_ptr<iRRAM::COMPLEX>[]>(resolution);
        storage_interior[i] = std::make_unique<std::unique_ptr<iRRAM::POLYNOMIAL>[]>(resolution);
    }

printf("BEGIN: reiteration = %lu\n", reiteration_count++);
    clock_t start, end;

    while(loop) {
        loop = false;
        for(int i = 0; i < resolution; i++) {
            for(int j = 0; j < resolution; j++) {
                if(canvas[i][j]) continue;
printf("working on (%d, %d)...\n", i, j);
                loop = true;
                auto const c = int_pair_to_COMPLEX(resolution, { i, j });
                bool flag = false;

printf("### EXTERIOR ###\n");
                /// EXTERIOR
                start = clock();
                auto z = (storage_exterior[i][j] ? *storage_exterior[i][j] : COMPLEX(REAL(0)));
                for(auto iter = iteration_exterior; iter < iteration_exterior_next; iter = iter + 1) {
                    z = z * z + c;
                    auto const d2 = real(z) * real(z) + imag(z) * imag(z);
                    if(choose(d2 > 4, d2 < 5) == 1) {
                        canvas[i][j] = std::make_unique<INTEGER>(-iter);
                        storage_exterior[i][j] = nullptr;
                        flag = true;
                        break;
                    }
                }
                end = clock();
printf("Elapsed time: %Lf\n", ((long double)(end - start)) / CLOCKS_PER_SEC);
                print_current(resolution, canvas);
                if(flag) continue;

printf("### INTERIOR ###\n");
                /// INTERIOR
                start = clock();
                auto const P = iRRAM::POLYNOMIAL(2, { c, REAL(0), REAL(1) });
                auto p = (storage_interior[i][j] ? *storage_interior[i][j] : iRRAM::POLYNOMIAL(1, { REAL(0), REAL(1) }));
                for(auto iter = iteration_interior; iter < iteration_interior_next; iter = iter + 1) {
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
                end = clock();
printf("Elapsed time: %Lf\n", ((long double)(end - start)) / CLOCKS_PER_SEC);
                print_current(resolution, canvas);
                if(flag) continue;

printf("### BOUNDARY ###\n");
                /// HELPER GRID & BOUNDARY
                start = clock();
                bool grid_exterior = false;
                bool grid_interior = false;
                for(iRRAM::INTEGER x = 1; x < iteration_power2; x = x + 2) {
                    for(iRRAM::INTEGER y = 1; y < iteration_power2; y = y + 2) {
                        auto const cur = INTEGER_pair_to_COMPLEX(resolution, { x, y });
                        auto const d = cur - c;
                        auto const d2 = real(d) * real(d) + imag(d) * imag(d);
                        if(choose(d2 < distance2, d2 > distance2 / 2) == 1) {
                            bool inner_flag = false;
                            z = REAL(0);
                            for(iRRAM::INTEGER iter = 1; iter < iteration_exterior_next; iter = iter + 1) {
                                z = z * z + c;
                                auto const d2 = real(z) * real(z) + imag(z) * imag(z);
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
                end = clock();
printf("Elapsed time: %Lf\n", ((long double)(end - start)) / CLOCKS_PER_SEC);
                print_current(resolution, canvas);
            END:;
            }
        }

        iteration_exterior = iteration_exterior_next;
        iteration_exterior_next = (iteration_exterior_next << 1);
        iteration_interior = iteration_interior_next;
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

    auto canvas = std::make_unique<std::unique_ptr<std::unique_ptr<iRRAM::INTEGER>[]>[]>(resolution);
    for(int i = 0; i < resolution; i++) {
        canvas[i] = std::make_unique<std::unique_ptr<iRRAM::INTEGER>[]>(resolution);
    }

    iRRAM::exec(compute, resolution, canvas);

    for(int j = 0; j < resolution; j++) {
		for(int i = 0; i < resolution; i++) {
			std::cout << swrite(*canvas[i][j]) << ' ';
		} std::cout << '\n';
	}
}
