/* Robust Estimation of Multiple Motions
 *
 * This motion estimation based on
 * M.J.Black and P.Anandan, "The Robust Estimation of Multiple Motions: Parametric and Piecewise-Smooth Flow Fields," Computer Vision and Image Understanding Vol. 63, No. 1, 1996, pp. 75-104.
 */
#ifndef LIB_ImgClass_AffineParametric
#define LIB_ImgClass_AffineParametric

#include <string>
#include <vector>

#include "Vector.h"
#include "BlockMatching.h"


	const double sigma = 0.2 / sqrt(2.0); //4.0 / sqrt(2.0);

template <class T>
class AffineParametric {
	private:
		int _width;
		int _height;
		double sigma;
		int IterMax;
		double Error_Min_Threshold;

		ImgVector<T> interest;
		ImgVector<T> reference;
		std::vector<std::vector<VECTOR_2D<int> > > regions;

		ImgVector<VECTOR_2D<double> > grad;
		ImgVector<double> dt;
		ImgVector<std::vector<double> > u_affine;

	public:
		std::vector<ImgVector<VECTOR_2D<double> > > AffineParametric(const int IterMax);

	protected:
		void IRLS_AffineParametric_region(std::vector<double>* u_affine, const std::vector<VECTOR_2D<int> >& region);
		std::vector<double> Error_a_region(const std::vector<double>* u_affine, const std::vector<VECTOR_2D<int> >& region);
		std::vector<double> sup_Error_aa_region(const std::vector<VECTOR_2D<int> >& region);
		double Error_region(const std::vector<double>* u_affine, const std::vector<VECTOR_2D<int> >& region);
}


namespace LIB_ImgClass_AffineParametric {
	template <class T>
	inline
	double
	Geman_McClure_rho(const T& x, const T& sigma)
	{
		return POW2(x) / (sigma + POW2(x));
	}

	template <class T>
	inline
	double
	Geman_McClure_psi(const T& x, const T& sigma)
	{
		return 2.0 * x * sigma / POW2(sigma + POW2(x));
	}

	template <class T>
	inline
	double
	Lorentzian_rho(const T& x, const T& sigma)
	{
		return log1p(0.5 * POW2(x / sigma));
	}

	template <class T>
	inline
	double
	Lorentzian_psi(const T& x, const T& sigma)
	{
		return 2.0 * x / (2.0 * POW2(sigma) + POW2(x));
	}
}

#include "AffineParametric_private.h"

#endif

