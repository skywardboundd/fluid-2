#include <bits/stdc++.h>

using namespace std;

constexpr size_t N = 36, M = 84;
// constexpr size_t N = 14, M = 5;
constexpr size_t T = 1'000;
constexpr std::array<pair<int, int>, 4> deltas{{{-1, 0}, {1, 0}, {0, -1}, {0, 1}}};

constexpr size_t s_N = 32, s_K = 16;

size_t THREADS = 4;
// char field[N][M + 1] = {
//     "#####",
//     "#.  #",
//     "#.# #",
//     "#.# #",
//     "#.# #",
//     "#.# #",
//     "#.# #",
//     "#.# #",
//     "#...#",
//     "#####",
//     "#   #",
//     "#   #",
//     "#   #",
//     "#####",
// };

char t_field[N][M + 1] = {
    "####################################################################################",
    "#                                                                                  #",
    "#                                                                                  #",
    "#                                                                                  #",
    "#                                                                                  #",
    "#                                                                                  #",
    "#                                       .........                                  #",
    "#..............#            #           .........                                  #",
    "#..............#            #           .........                                  #",
    "#..............#            #           .........                                  #",
    "#..............#            #                                                      #",
    "#..............#            #                                                      #",
    "#..............#            #                                                      #",
    "#..............#            #                                                      #",
    "#..............#............#                                                      #",
    "#..............#............#                                                      #",
    "#..............#............#                                                      #",
    "#..............#............#                                                      #",
    "#..............#............#                                                      #",
    "#..............#............#                                                      #",
    "#..............#............#                                                      #",
    "#..............#............#                                                      #",
    "#..............#............################                     #                 #",
    "#...........................#....................................#                 #",
    "#...........................#....................................#                 #",
    "#...........................#....................................#                 #",
    "##################################################################                 #",
    "#                                                                                  #",
    "#                                                                                  #",
    "#                                                                                  #",
    "#                                                                                  #",
    "#                                                                                  #",
    "#                                                                                  #",
    "#                                                                                  #",
    "#                                                                                  #",
    "####################################################################################",
};

template <int N, int K>
class Fixed {
public:
    constexpr Fixed() : value(0) {}
    constexpr Fixed(double d) : value(static_cast<int64_t>(d * (1LL << K))) {}

    static constexpr Fixed from_raw(int64_t raw) {
        Fixed f;
        f.value = raw;
        return f;
    }

    double to_double() const {
        return static_cast<double>(value) / (1LL << K);
    }

    constexpr Fixed operator+(Fixed other) const {
        return Fixed::from_raw(value + other.value);
    }

    constexpr Fixed operator-(Fixed other) const {
        return Fixed::from_raw(value - other.value);
    }

    constexpr Fixed operator*(Fixed other) const {
        return Fixed::from_raw((value * other.value) >> K);
    }

    constexpr Fixed operator/(Fixed other) const {
        return Fixed::from_raw((value << K) / other.value);
    }

    Fixed& operator+=(Fixed other) {
        value += other.value;
        return *this;
    }

    Fixed& operator-=(Fixed other) {
        value -= other.value;
        return *this;
    }

    Fixed& operator*=(Fixed other) {
        value = (value * other.value) >> K;
        return *this;
    }

    Fixed& operator/=(Fixed other) {
        value = (value << K) / other.value;
        return *this;
    }


    constexpr bool operator==(const Fixed& other) const = default;

    constexpr auto operator<=>(const Fixed& other) const = default;

private:
    int64_t value; 
};


#ifdef FLOAT
#error "FLOAT is already defined"
#endif
#ifdef DOUBLE
#error "DOUBLE is already defined"
#endif
#ifdef FIXED
#error "FIXED is already defined"
#endif
#ifdef FAST_FIXED
#error "FAST_FIXED is already defined"
#endif
#ifdef S
#error "S is defined"
#endif

#define FLOAT            float
#define DOUBLE           double
#define FAST_FIXED(N, K) types::FastFixed<N, K>
#define FIXED(N, K)      types::Fixed<N, K>



template <int N, int K>
std::ostream& operator<<(std::ostream& os, const Fixed<N, K>& f) {
    return os << f.to_double();
}   

template <int N, int K>
constexpr Fixed<N, K> inf = Fixed<N, K>::from_raw(std::numeric_limits<int64_t>::max());

template <int N, int K>
constexpr Fixed<N, K> eps = Fixed<N, K>::from_raw(1);

template <int N, int K>
Fixed<N, K> abs(Fixed<N, K> x) {
    return x.value < 0 ? Fixed<N, K>::from_raw(-x.value) : x;
}

template <int N, int K>
ostream &operator<<(ostream &out, Fixed<N, K> x) {
    return out << x.v / (double) (1 << 16);
}

template <typename T, size_t N, size_t M>
struct VectorField {
    array<T, deltas.size()> v[N][M];
    T &add(int x, int y, int dx, int dy, T dv) {
        return get(x, y, dx, dy) += dv;
    }

    T &get(int x, int y, int dx, int dy) {
        size_t i = 0;
        if(dx == 1 && dy == 0) i = 1;
        else if (dx == 0 && dy == -1) i = 2;
        else if (dx == 0 && dy == 1) i = 3;
        return v[x][y][i];
    }
};



template <typename T, size_t N, size_t M>
class FluidSimulator {
public:
    FluidSimulator();
    void run_simulation(size_t steps);

private:
    char field[N][M + 1];
    int dirs[N][M]{};
    T p[N][M]{}, old_p[N][M];
    VectorField<T, N, M> velocity{}, velocity_flow{};
    T rho[256];
    int last_use[N][M]{};
    int UT = 0;

    tuple<T, bool, pair<int, int>> propagate_flow(int x, int y, T lim);
    T random01();
    bool propagate_move(int x, int y, bool is_first=false);
    T move_prob(int x, int y);
    void propagate_stop(int x, int y, bool force=false);
};


mt19937 rnd(1337);

template <typename T, size_t N, size_t M>
FluidSimulator<T, N, M>::FluidSimulator() {
    rho[' '] = T(0.01);
    rho['.'] = T(1000);

    memcpy(field, t_field, sizeof(t_field));

    memset(p, 0, sizeof(p));
    memset(old_p, 0, sizeof(old_p));
}


template <typename T, size_t N, size_t M>
tuple<T, bool, pair<int, int>> FluidSimulator<T, N, M>::propagate_flow(int x, int y, T lim) {
    last_use[x][y] = UT - 1;
    T ret = 0;
    for (auto [dx, dy] : deltas) {
        int nx = x + dx, ny = y + dy;
        if (field[nx][ny] != '#' && last_use[nx][ny] < UT) {
            auto cap = velocity.get(x, y, dx, dy);
            auto flow = velocity_flow.get(x, y, dx, dy);
            if (flow == cap) {
                continue;
            }
            auto vp = min(lim, cap - flow);
            if (last_use[nx][ny] == UT - 1) {
                velocity_flow.add(x, y, dx, dy, vp);
                last_use[x][y] = UT;
                return {vp, 1, {nx, ny}};
            }
            auto [t, prop, end] = propagate_flow(nx, ny, vp);
            ret += t;
            if (prop) {
                velocity_flow.add(x, y, dx, dy, t);
                last_use[x][y] = UT;
                return {t, prop && end != pair(x, y), end};
            }
        }
    }
    last_use[x][y] = UT;
    return {ret, 0, {0, 0}};
}

template <typename T, size_t N, size_t M>
T FluidSimulator<T, N, M>::random01() {
    return T::from_raw((rnd() & ((1 << 16) - 1)));
}

template <typename T, size_t N, size_t M>
void FluidSimulator<T, N, M>::propagate_stop(int x, int y, bool force) {
    if (!force) {
        bool stop = true;
        for (auto [dx, dy] : deltas) {
            int nx = x + dx, ny = y + dy;
            if (field[nx][ny] != '#' && last_use[nx][ny] < UT - 1 && velocity.get(x, y, dx, dy) > 0) {
                stop = false;
                break;
            }
        }
        if (!stop) {
            return;
        }
    }
    last_use[x][y] = UT;
    for (auto [dx, dy] : deltas) {
        int nx = x + dx, ny = y + dy;
        if (field[nx][ny] == '#' || last_use[nx][ny] == UT || velocity.get(x, y, dx, dy) > 0) {
            continue;
        }
        propagate_stop(nx, ny);
    }
}

template <typename T, size_t N, size_t M>
T FluidSimulator<T, N, M>::move_prob(int x, int y) {
    T sum = 0;
    for (size_t i = 0; i < deltas.size(); ++i) {
        auto [dx, dy] = deltas[i];
        int nx = x + dx, ny = y + dy;
        if (field[nx][ny] == '#' || last_use[nx][ny] == UT) {
            continue;
        }
        auto v = velocity.get(x, y, dx, dy);
        if (v < 0) {
            continue;
        }
        sum += v;
    }
    return sum;
}

template <typename T, size_t N, size_t M>
struct ParticleParams {
    char type;
    T cur_p;
    array<T, deltas.size()> v;

    void swap_with(int x, int y, char field[N][M + 1], T p[N][M], VectorField<T, N, M> &velocity) {
        swap(field[x][y], type);
        swap(p[x][y], cur_p);
        swap(velocity.v[x][y], v);
    }
};

template <typename T, size_t N, size_t M>
bool FluidSimulator<T, N, M>::propagate_move(int x, int y, bool is_first) {
    last_use[x][y] = UT - is_first;
    bool ret = false;
    int nx = -1, ny = -1;
    do {
        std::array<T, deltas.size()> tres;
        T sum = 0;
        for (size_t i = 0; i < deltas.size(); ++i) {
            auto [dx, dy] = deltas[i];
            int nx = x + dx, ny = y + dy;
            if (field[nx][ny] == '#' || last_use[nx][ny] == UT) {
                tres[i] = sum;
                continue;
            }
            auto v = velocity.get(x, y, dx, dy);
            if (v < 0) {
                tres[i] = sum;
                continue;
            }
            sum += v;
            tres[i] = sum;
        }

        if (sum == 0) {
            break;
        }

        T p = random01() * sum;
        size_t d = std::ranges::upper_bound(tres, p) - tres.begin();

        auto [dx, dy] = deltas[d];
        nx = x + dx;
        ny = y + dy;
        assert(velocity.get(x, y, dx, dy) > 0 && field[nx][ny] != '#' && last_use[nx][ny] < UT);

        ret = (last_use[nx][ny] == UT - 1 || propagate_move(nx, ny, false));
    } while (!ret);
    last_use[x][y] = UT;
    for (size_t i = 0; i < deltas.size(); ++i) {
        auto [dx, dy] = deltas[i];
        int nx = x + dx, ny = y + dy;
        if (field[nx][ny] != '#' && last_use[nx][ny] < UT - 1 && velocity.get(x, y, dx, dy) < 0) {
            propagate_stop(nx, ny);
        }
    }
    if (ret) {
        if (!is_first) {
            ParticleParams<T, N, M> pp{};
            pp.swap_with(x, y, this->field, this->p, this->velocity);
            pp.swap_with(nx, ny, this->field, this->p, this->velocity);
            pp.swap_with(x, y, this->field, this->p, this->velocity);

        }
    }
    return ret;
}

#include<thread>

template <typename T, size_t N, size_t M>
void FluidSimulator<T, N, M>::run_simulation(size_t steps) {
    rho[' '] = 0.01;
    rho['.'] = 1000;
    T g = 9.8;

    auto f_dirs = [&](size_t n1, size_t n2) {
        for (size_t x = n1; x < n2; ++x) {
            for (size_t y = 0; y < M; ++y) {
                if (field[x][y] == '#')
                    continue;
                for (auto [dx, dy] : deltas) {
                    dirs[x][y] += (field[x + dx][y + dy] != '#');
                }
            }
        }
        
    };

    std::vector<std::thread> threads;
    for (size_t i = 0; i < THREADS; ++i) {
        size_t j = N / THREADS + 1;
        size_t n1 = i * j;
        size_t n2 = n1 + j;
        threads.emplace_back(f_dirs, n1, n2);
    }

    for (auto &thread : threads) {
        thread.join();
    }

    // for (size_t x = 0; x < N; ++x) {
    //     for (size_t y = 0; y < M; ++y) {
    //         if (field[x][y] == '#')
    //             continue;
    //         for (auto [dx, dy] : deltas) {
    //             dirs[x][y] += (field[x + dx][y + dy] != '#');
    //         }
    //     }
    // }

    

    for (size_t i = 0; i < steps; ++i) {
        
        T total_delta_p = 0;
        // Apply external forces
        for (size_t x = 0; x < N; ++x) {
            for (size_t y = 0; y < M; ++y) {
                if (field[x][y] == '#')
                    continue;
                if (field[x + 1][y] != '#')
                    velocity.add(x, y, 1, 0, g);
            }
        }

        // Apply forces from p
        memcpy(old_p, p, sizeof(p));
        for (size_t x = 0; x < N; ++x) {
            for (size_t y = 0; y < M; ++y) {
                if (field[x][y] == '#')
                    continue;
                for (auto [dx, dy] : deltas) {
                    int nx = x + dx, ny = y + dy;
                    if (field[nx][ny] != '#' && old_p[nx][ny] < old_p[x][y]) {
                        auto delta_p = old_p[x][y] - old_p[nx][ny];
                        auto force = delta_p;
                        auto &contr = velocity.get(nx, ny, -dx, -dy);
                        if (contr * rho[(int) field[nx][ny]] >= force) {
                            contr -= force / rho[(int) field[nx][ny]];
                            continue;
                        }
                        force -= contr * rho[(int) field[nx][ny]];
                        contr = 0;
                        velocity.add(x, y, dx, dy, force / rho[(int) field[x][y]]);
                        p[x][y] -= force / dirs[x][y];
                        total_delta_p -= force / dirs[x][y];
                    }
                }
            }
        }

        // Make flow from velocities
        velocity_flow = {};
        bool prop = false;
        do {
            UT += 2;
            prop = 0;
            for (size_t x = 0; x < N; ++x) {
                for (size_t y = 0; y < M; ++y) {
                    if (field[x][y] != '#' && last_use[x][y] != UT) {
                        auto [t, local_prop, _] = propagate_flow(x, y, 1);
                        if (t > 0) {
                            prop = 1;
                        }
                    }
                }
            }
        } while (prop);

        // Recalculate p with kinetic energy
        for (size_t x = 0; x < N; ++x) {
            for (size_t y = 0; y < M; ++y) {
                if (field[x][y] == '#')
                    continue;
                for (auto [dx, dy] : deltas) {
                    auto old_v = velocity.get(x, y, dx, dy);
                    auto new_v = velocity_flow.get(x, y, dx, dy);
                    if (old_v > 0) {
                        assert(new_v <= old_v);
                        velocity.get(x, y, dx, dy) = new_v;
                        auto force = (old_v - new_v) * rho[(int) field[x][y]];
                        if (field[x][y] == '.')
                            force *= 0.8;
                        if (field[x + dx][y + dy] == '#') {
                            p[x][y] += force / dirs[x][y];
                            total_delta_p += force / dirs[x][y];
                        } else {
                            p[x + dx][y + dy] += force / dirs[x + dx][y + dy];
                            total_delta_p += force / dirs[x + dx][y + dy];
                        }
                    }
                }
            }
        }

        UT += 2;
        prop = false;
        for (size_t x = 0; x < N; ++x) {
            for (size_t y = 0; y < M; ++y) {
                if (field[x][y] != '#' && last_use[x][y] != UT) {
                    if (random01() < move_prob(x, y)) {
                        prop = true;
                        propagate_move(x, y, true);
                    } else {
                        propagate_stop(x, y, true);
                    }
                }
            }
        }

        if (prop) {
            cout << "Tick " << i << ":\n";
            for (size_t x = 0; x < N; ++x) {
                cout << field[x] << "\n";
            }
        }
    }
}

#include <chrono>
#include <thread>

int main() {

    auto start_time = std::chrono::high_resolution_clock::now();

    FluidSimulator<Fixed<32, 16>, N, M> simulator;
    simulator.run_simulation(T);

    auto end_time = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed = end_time - start_time;
    std::cout << "Simulation completed in: " << elapsed.count() << " seconds\n";
    
    return 0;
}
