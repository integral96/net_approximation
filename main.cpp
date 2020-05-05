#include <solution_equation.hpp>

#include <thread>
#include <mutex>
#include <numeric>

template<size_t N, size_t M>
using  cranc_type = cranc_nicolson1<ubls::fixed_matrix<double, N + 1, M + 1>, N, double, double, container_algebra, default_operation<N>>;

template<size_t N, size_t M>
using  state_type = ubls::fixed_matrix<double, N + 1, M + 1>;

static auto f_f([](auto x) {
    return x*x*(5 - x);
});
static auto p_0([](auto t) {
    return t*t;
});
static auto p_1([](auto t) {
    return 70 + t;
});

template<size_t N, size_t M>
struct system_equation
{
    const double sigma, R, b;
    system_equation(const double sigma, const double R, const double b):sigma(sigma), R(R), b(b) {}
    void operator() (const state_type<N, M>& F, state_type<N, M>& U, double t, double x)
    {
        for(size_t i = 0; i < N; ++i)
            U(i, 0) = sigma*F(i, 0);
        for(size_t i = 0; i < N; ++i)
            for(size_t j = 0; j < N; ++j)
                U(i, j) = .5*b*R*(U(i - 1, j) - 2*U(i, j) + U(i + 1, j));
        for(size_t j = 0; j < M - 1; ++j)
            U(0, j + 1) = U(M, j + 1) = 0;
    }
};

int main() {

    cranc_nicolson<double, 10, 10, decltype (f_f), decltype (p_0), decltype (p_1)> solv(.1, .1, 1.);
    solv.boundary_conditions(f_f, p_0, p_1);
    solv.print();
    cranc_type<10, 10> crnc;
    system_equation<10, 10> system(10., 28., 8./3.);
    ubls::fixed_matrix<double, 11, 11> U;
    for(size_t i = 0; i < 10; ++i)
        for(size_t j = 0; j < 10; ++j)
            crnc.do_step(system, U, j*.1, 1., i*.1, 1);

    std::cout << "BINGO\n";
	return 0;
}
