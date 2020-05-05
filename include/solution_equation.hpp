#pragma once

#include <iostream>
#include <vector>
#include <tuple>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/spirit/include/karma.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/numeric/ublas/triangular.hpp>
 #include <boost/numeric/ublas/lu.hpp>
 #include <boost/numeric/ublas/io.hpp>


namespace ubls = boost::numeric::ublas;
namespace krm  = boost::spirit::karma;


template<typename state_type>
void resize(const state_type& in, state_type& out) {
    out.resize(std::size(in));
}

template <typename T, size_t N, size_t M>
void resize(const ubls::fixed_matrix<T, N, M> & in, ubls::fixed_matrix<T, N, M> & out) {

}

template<typename System, typename T>
void euler_step(System system, ubls::matrix<T>& u, const T t, const T dt, const T dx) {
    ubls::matrix<T> k(u.size1(), u.size2());
    system(u, k, t);
    for(size_t i = 0; i < u.size1(); ++i) {
        for(size_t j = 0; j < u.size2(); ++j) {
            u(i, j) += dt*k(i, j) + dx*k(i, j);
        }
    }
}

struct container_algebra
{
    template<typename Op, typename ...Args>
    void for_eachN(Op&& op, Args&& ... args) {
        for(size_t i = 0; i < sizeof... (args); i++)
        std::invoke(std::forward<Op>(op), std::forward<Args>(args)...);
    }
};

template <typename ...Args> struct sum_type {};
template <typename T> struct sum_type<T> {
    using type = T;
};
template <typename T, typename ...Args> struct sum_type<T, Args...> {
    using type = decltype (T() + typename sum_type<Args...>::type());
};
template <typename ...Args>
using sum_type_t = typename sum_type<Args...>::type;

template <typename T, typename ...Args>
sum_type_t<T, Args...> sum_m(T t, Args&&... args) {
    return t + sum_m(args...);
}

template<size_t N = 10>
struct default_operation
{
    template<typename ...ArgsF>
    struct scale_sumN {
        const std::tuple<ArgsF...> alpha;
        using result_type = void;
        scale_sumN(ArgsF&& ... args) : alpha(std::forward<ArgsF>(args)...) {}
        template<typename T0, typename ...Args>
        void operator() (T0 t0, Args&& ... args) {
            t0 = sum_m(std::forward<Args>(std::get<N>(alpha)*args)...);
        }
    };
};
template <typename state_type, size_t N,
          typename value_type = double,
          typename time_type = value_type,
          typename algebra = container_algebra,
          typename operations = default_operation<N>>
class cranc_nicolson1
{
public:
    cranc_nicolson1() {}
    template<typename System>
    void do_step(System& system, state_type& U, time_type t, time_type dt, time_type x, time_type dx) {
        adjust_size(U);
        const value_type one = 1;
        const time_type dt2 = dt/2, dt3 = dt/3, dt6 = dt/6;
        const time_type dx2 = dx/2, dx3 = dx/3, dx6 = dx/6;
        using  scale_sum2 = typename operations::template scale_sumN<value_type, time_type>;
        using  scale_sum5 = typename operations::template scale_sumN<value_type, time_type, time_type, time_type, time_type>;
        system(U, K1, x, t);
        m_algebra.for_eachN(U_tmp, U, K1/*, scale_sum2(one, dt2*/);
        system(U_tmp, K2, x + dx2, t+dt2);
        m_algebra.for_eachN(U_tmp, U, K2/*, scale_sum2(one, dt2)*/);
        system(U_tmp, K3, x + dx3, t+dt2);
        m_algebra.for_eachN(U_tmp, U, K3/*, scale_sum2(one, dt)*/);
        system(U_tmp, K4, x + dx6, t+dt);
        m_algebra.for_eachN(U, U, K1, K2, K3, K4/*, scale_sum5(one, dt6, dt3, dt3, dt6)*/);
    }
private:
    state_type U_tmp, K1, K2, K3, K4;
    algebra m_algebra;
    void adjust_size(const state_type& U) {
        resize(U, U_tmp);
        resize(U, K1);
        resize(U, K2);
        resize(U, K3);
        resize(U, K4);
    }
};
/////
template <typename T, size_t N, size_t M, typename F_x, typename P0_t, typename P1_t>
class cranc_nicolson
{
private:
    ubls::fixed_matrix<T, N + 1, M + 1> U;
    ubls::fixed_matrix<T, N + 1, M + 1> A;
    ubls::fixed_matrix<T, N + 1, M + 1> B;
    using pmatrix = ubls::permutation_matrix<size_t>;
    using cranc_types = boost::mpl::vector4<T, F_x, P0_t, P1_t>;
    T q_fix, h;
    T D;
public:
    cranc_nicolson(T s, T h, T d) : h(h), q_fix(s/(h*h)), D(d) {
        A(0, 0) = 1;
        A(N, M) = 1;
        for(size_t i = 0; i < N; ++i) {
            for(size_t j = 1; j < M; ++j) {
               if(i == j) {
                   A(i, j) = 1 + q_fix*D;
                   A(i, j + 1) = - q_fix*D/2;
                   A(i, j - 1) = - q_fix*D/2;
               }
            }
        }
    }
    void boundary_conditions(const F_x& f_x, const P0_t& p0_t, const P1_t& p1_t) {
        for(size_t j = 0; j < M; ++j) {
            U(0, j) = B(0, j) = p0_t(j);
            U(N, j) = B(N, j) = p1_t(j);
        }
        for(size_t i = 0; i < N; ++i) {
            U(i, 0) = B(i, 0) = f_x(i*h);
            U(i, M) = B(i, M) = f_x(i*h);
        }
        for(size_t j = 1; j <= M - 1; ++j) {
            for(size_t i = 1; i <= N - 1; ++i) {
                B(i, j) = U(i, j) + .5*q_fix*D*(U(i - 1, j) - 2*U(i, j) + U(i + 1, j));
            }
        }
        U = ubls::prod(invert_matrix(A), B);
    }
    ubls::fixed_matrix<T, N + 1, M + 1> invert_matrix(const ubls::fixed_matrix<T, N + 1, M + 1>& A_in) {
        ubls::fixed_matrix<T, N + 1, M + 1> A_inv;
        ubls::fixed_matrix<T, N + 1, M + 1> A_i(A_in);
        pmatrix pm(A_i.size1());
        auto res = ubls::lu_factorize(A_i, pm);
        A_inv.assign(ubls::identity_matrix<T>(A_i.size1()));
        ubls::lu_substitute(A_i, pm, A_inv);
        return  A_inv;
    }
    void print() {
        std::cout << krm::format_delimited(krm::columns(A.size2()) [krm::auto_], '\t', A.data()) << std::endl;
        std::cout << krm::format_delimited(krm::columns(this->invert_matrix(A).size2()) [krm::auto_], '\t', this->invert_matrix(A).data()) << std::endl;
        std::cout << krm::format_delimited(krm::columns(U.size2()) [krm::auto_], '\t', U.data()) << std::endl;
    }
};
