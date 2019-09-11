#include <iostream>
#include <cmath>
#include <boost/numeric/mtl/mtl.hpp>

using namespace mtl;

const double X1L {.0};
const double X1D {.5};
const double X1R {1.};
const double X2L {.0};
const double X2D {.5};
const double X2R {1.};
const double KAPPA {100.};
const double BI {10.};

const size_t N1 {11};
const size_t N2 {11};
const double EPS {.0001};

//===============Lambds
	static auto H1([](const double& X1R, const double& X1L)
	{
		return (X1R - X1L)/(N1 - 1);
	});
	static auto H2([](const double& X2R, const double& X2L){
		return (X2R - X2L)/(N2 - 1);
	});
//=================
template<template<typename> typename Cont, typename T>
void FDS01(Cont<T>& A0, Cont<T>& A1, Cont<T>& A2, Cont<T>& F) 
{
	T X1;
    T X2;
	T H12 = H1(X1R, X1L)/H2(X2R, X2L);
    T H21 = H2(X2R, X2L)/H1(X1R, X1L);
	//=========================
	static auto K([Kz {1.}](const T& X1, const T& X2) mutable
	{
		if (X1 < X1D && X2 < X2D)
		{
			return KAPPA;
		}
		return Kz;
	});
	
	auto K1([](const T& X1, const T& X2) {
		return .5*(K(X1, X2 - .5*H2(X2R, X2L)) + K(X1, X2 + .5*H2(X2R, X2L)));
	});
	auto K2([](const T& X1, const T& X2) {
		return .5*(K(X1 - .5*H1(X1R, X1L), X2) + K(X1 + .5*H1(X1R, X1L), X2));
	});
//=======================
	for (size_t j = 1; j < N2; j++)
    {
        X2 = X2L + (j - 1)*H2(X2R, X2L);
        for (size_t i = 1; i < N1; i++)
        {
            X1 = X1L + (i - 1)*H1(X1R, X1L);
            A1[i-1][j] = H21*K1(X1 - .5*H1(X1R, X1L), X2);
            A1[i][j] = H21*K1(X1 + .5*H1(X1R, X1L), X2);
            A2[i][j-1] = H12*K2(X1, X2 - .5*H2(X2R, X2L));
            A2[i][j] = H12*K2(X1, X2 + .5*H2(X2R, X2L));
            A0[i][j] = A1[i-1][j] + A1[i][j] + A2[i][j-1] + A2[i][j];
            F[i][j] = .0;
        }
        A1[0][j] = .0;
        A2[0][j] = .0;
        A0[0][j] = 1.0;
        F[0][j] = .0;
    }
        for (size_t i = 1; i < N1; i++)
        {
            X1 = X1L + (i - 1)*H1(X1R, X1L);
            A1[i-1][0] = H21*K1(X1 - .5*H1(X1R, X1L), X2L);
            A1[i][0] = H21*K1(X1 + .5*H1(X1R, X1L), X2L);
            A0[i][0] = A1[i-1][0] + A1[i][0] + A2[i][0];
            F[i][0] = .0;
            A1[i-1][N2-1] = H21*K1(X1 - .5*H1(X1R, X1L), X2R);
            A1[i][N2-1] = H21*K1(X1 + .5*H1(X1R, X1L), X2R);
            A2[i][N2-1] = .0;
            A0[i][N2-1] = A1[i-1][N2-1] + A1[i][N2-1] + A2[i][N2-1];
            F[i][N2-1] = .0;
        }
    for (size_t j = 1; j < N2; j++)
    {
        X2 = X2L + (j - 1)*H2(X2R, X2L);
        A1[N1-1][j] = .0;
        A2[N1-1][j-1] = .5*H12*K2(X1R, X2 - .5*H2(X2R, X2L));
        A2[N1-1][j] = .5*H12*K2(X1R, X2 + .5*H2(X2R, X2L));
        A0[N1-1][j] = A1[N1-2][j] + A2[N1-1][j-1] + A2[N1-1][j] + H2(X2R, X2L)*BI;
        F[N1-1][j] = H2(X2R, X2L)*BI;
    }
    A1[0][0] = .0;
    A2[0][0] = .0;
    A0[0][0] = 1.;
    F[0][0] = .0;

    A1[0][N2-1] = .0;
    A2[0][N2-1] = .0;
    A0[0][N2-1] = 1.;
    F[0][N2-1] = .0;

    A1[N1-1][0] = .0;
    A0[N1-1][0] = A1[N1-2][0] + A2[N1-1][0] + .5*H2(X2R, X2L)*BI;
    F[N1-1][0] = .5*H2(X2R, X2L)*BI;

    A1[N1-1][N2-1] = .0;
    A2[N1-1][N2-1] = .0;
    A0[N1-1][N2-1] = A1[N1-2][N2-1] + A2[N1-1][N2-2] + .5*H2(X2R, X2L)*BI;
    F[N1-1][N2-1] = .5*H2(X2R, X2L)*BI;
}
template<template<typename> typename Cont, typename T>
void SOLVE1(const size_t& N, size_t L, Cont<T>& A, Cont<T>& Y, Cont<T>& F, const double& EPS)
{
	size_t I1R = N + 1;
	size_t I2R = 2*N + L;
	size_t IG = 3*N + L;
	size_t IT = IG + N;
	size_t ITL = IT + L;
	size_t IQ = IT + N + L;

	size_t NIT;
	T RR, RRI, BK, AK, EPSNIT, TQ, RR0;
	Cont<T> G(N, N);
	try
	{
		for (size_t j = 0; j < N; j++)
		{
			G[j][j] = A[j][j] - A[I1R+j-1][I1R+j-1]*A[IT+j-1][IT+j-1]*A[I2R+j-1][I2R+j-1] + A[I2R+j-1][I2R+j-1] 
			- A[I2R+j-L][I2R+j-L]*A[IT+j-L][IT+j-L]*A[I1R+j-L][I1R+j-L] + A[I2R+j-L][I2R+j-L];
			A[IT+j][IT+j] = 1./G[j][j];
			A[IG+j][IG+j] = std::sqrt(A[IT+j][IT+j]);
		}
	}	
	catch(const std::exception& e)
	{
		std::cerr << e.what() <<  " 0 HELP!!!" <<'\n';
	}
	try 
	{
		for (size_t j = 1; j <= N; j++)
		{
			A[j][j] *= A[IG+j][IG+j]*A[IG+j][IG+j] - 2.;
			A[I1R+j][I1R+j] = A[IG+j][IG+j]*A[I1R+j][I1R+j]*A[IG+j+1][IG+j+1];
			A[I2R+j][I2R+j] = A[IG+j][IG+j]*A[I2R+j][I2R+j]*A[IG+j+L][IG+j+L];
		} 
	}
	catch(const std::exception& e)
	{
		std::cerr << e.what() <<  " 1 HELP!!!" <<'\n';
	}
	try 
	{
		for (size_t j = N; j > 0; j--)
		{
			A[ITL+j][ITL+j] = Y[j][j]/A[IG+j][IG+j];
			Y[j][j] = A[ITL+j][ITL+j] - A[I1R+j][I1R+j]*A[ITL+j+1][ITL+j+1] - A[I2R+j][I2R+j]*A[ITL+j+L][ITL+j+L];
		}
	}
	catch(const std::exception& e)
	{
		std::cerr << e.what() <<  " 2 HELP!!!" <<'\n';
	}
	RR0 = .0;
	
	try 
	{
		
		for (size_t j = 0; j < N; j++)
		{
			A[IT+j][IT+j] = F[j][j]*A[IG+j][IG+j] - A[j][j]*A[ITL+j][ITL+j] - Y[j][j]
				+ A[I1R+j-1][I1R+j-1]*A[IT+j-1][IT+j-1] + A[I2R+j-L][I2R+j-L]*A[IT+j-L][IT+j-L];
			F[j][j] = A[IT+j][IT+j] - A[ITL+j][ITL+j];
			RR0 += F[j][j]*F[j][j];
		}
		std::cout << "RR0 = " << RR0 << std::endl;
	}
	catch(const std::exception& e)
	{
		std::cerr << e.what() <<  " 3 HELP!!!" <<'\n';
	}
	NIT = 0;
	RR = RR0;
	RRI  = .0;
	do
	{
		
		BK = RR*RRI;
		RRI = 1./RR;
		try 
		{
			for (size_t j = N; j > 0; j--)
			{
				A[IQ+j][IQ+j] = F[j][j] + BK*A[IQ+j][IQ+j];
				A[ITL+j][ITL+j] = A[IQ+j][IQ+j] + A[I1R+j][I1R+j]*A[ITL+j+1][ITL+j+1] + A[I2R+j][I2R+j]*A[ITL+j+L][ITL+j+L];
			}
		}
		catch(const std::exception& e)
		{
			std::cerr << e.what() <<  " 4 HELP!!!" <<'\n';
		}
		TQ = .0;
		try 
		{
			for (size_t j = 0; j < N; j++)
			{
				G[j][j] = A[ITL+j][ITL+j];
				A[ITL+j][ITL+j] = A[IQ+j][IQ+j] + A[j][j]*G[j][j] + A[I1R+j-1][I1R+j-1]*A[ITL+j-1][ITL+j-1]
					+ A[I2R+j-L][ITL+j-L];
				A[IT+j][IT+j] = G[j][j] + A[ITL+j][ITL+j];
				TQ += A[IT+j][IT+j]*A[IQ+j][IQ+j];
			}
		}
		catch(const std::exception& e)
		{
			std::cerr << e.what() <<  " 5 HELP!!!" <<'\n';
		}
		AK = RR/TQ;
		RR = .0;
			for (size_t j = 1; j <= N; j++)
			{
				Y[j][j] += AK*A[IQ+j][IQ+j];
				F[j][j] -= AK*A[IT+j][IT+j];
				RR += F[j][j]*F[j][j];
			}
		NIT++;
		EPSNIT = std::sqrt(RR/RR0);
		std::cout << "NIT & EPSNIT" << NIT << " " << EPSNIT << std::endl;
	} while (EPSNIT > EPS);
		for (size_t j = N; j > 0; j--)
		{
			A[ITL+j][ITL+j] = Y[j][j] + A[I1R+j][I1R+j]*A[ITL+j+1][ITL+j+1] + A[I2R+j][I2R+j]*A[ITL+j+L][ITL+j+L];
			Y[j][j] = A[IG+j][IG+j]*A[ITL+j][ITL+j] ;
		}
	//std::cout << "NIT & EPSNIT" << NIT << " " << EPSNIT << std::endl;
}
int main(int, char**)
{
	size_t N {N1*N2};
    // dense2D<double> A0(8*N+1, 8*N+1), A1(8*N+1, 8*N+1), A2(8*N+1, 8*N+1), F(8*N+1, 8*N+1);
	dense2D<double> A0(N1, N2), A1(N1, N2), A2(N1, N2), F(N1, N2);
	//dense2D<double> A(N, N), Y(7*N+1, 7*N+1), F1(8*N+1, 8*N+1);
    FDS01(A0, A1, A2, F);
	//SOLVE1(N, N1, A0, A1, F, EPS);
    //A = 0;
	std::cout << "A0 is \n" << A0 << "\n";
    std::cout << "A2 is \n" << A1 << "\n";
	std::cout << "A2 is \n" << A2 << "\n";
    std::cout << "F is \n" << F << "\n";
	
    return 0;
}