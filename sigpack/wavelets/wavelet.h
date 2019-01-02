// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// @author Juan Pablo Cordova Echeverria
// @email jpcordovae@gmail.com

// Wavelet Version
// 0.1.0 JPCE 2017-06-14 Initial codes, virtual class and 

#ifndef SP_WAVELET_H
#define SP_WAVELET_H

#include "sigpack.h"

namespace sp {
	///
	/// @defgroup wavelets Wavelets
	/// \brief base class for wavelets and derivates for wavelets signal analisys
	/// The shifting values are keeped in to the class to avoid call ovehead 
	/// 
	/// @{

	///
	/// \brief Wavelet class.
	///
	/// Implements Wavelets basic functions
	///
	template<class T>
	class Wavelet
	{
	public:
		Wavelet() {
			// default values
			m_t_shift = 0.0;
			m_scale = 1.0;
		}
		~Wavelet() {}

		////////////////////////////////////////////////////////////////////////////////////////////
		/// \brief continuous value of the mother wavelet.
		///
		/// The calue is calculated the continuoe mother walet for a group of time values
		/// @returns A T class vector
		/// @param _in time input vector
		////////////////////////////////////////////////////////////////////////////////////////////
		arma::Col<T> mother(arma::Col<T> &_in) 
		{
			arma::Col<T> ret(_in);
			for (auto &it : ret) {
				it = mother(it);
			}
			return ret;
		}

		////////////////////////////////////////////////////////////////////////////////////////////
		/// \brief continuous value of the child wavelet.
		///
		/// The calue is calculated the continuoe child walet for a group of time values
		/// @returns A T class vector
		/// @param _in time input vector
		/// @param _t_shift shift time for child wavelet
		/// @param _scale scale parameter forchild wavelet
		////////////////////////////////////////////////////////////////////////////////////////////
		arma::Col<T> child(arma::Col<T> &_in, double _t_shit, double _scale) {
			arma::Col<T> ret(_in);
			set_child(_t_shit, _scale);
			for (auto &it : ret) {
				it = child(it);
			}
			return ret;
		}

		////////////////////////////////////////////////////////////////////////////////////////////
		/// \brief continuous value of the child wavelet.
		///
		/// The calue is calculated the continuoe child walet for a group of time values
		/// @returns A T class vector
		/// @param _in time input vector
		////////////////////////////////////////////////////////////////////////////////////////////
		arma::Col<T> child(arma::Col<T> &_in) {
			arma::Col<T> ret(_in);
			for (auto &it : ret) {
				it = child(it);
			}
			return ret;
		}

		////////////////////////////////////////////////////////////////////////////////////////////
		/// \brief general description of a mother wavelet
		///
		/// calculate the value of the continue mother wavelet at time t
		/// @returns T class value
		/// @param _t time 
		////////////////////////////////////////////////////////////////////////////////////////////
		arma_inline T mother(double _t) {
			return psi(_t);
		}

		////////////////////////////////////////////////////////////////////////////////////////////
		/// \brief general description of a child wavelet
		///
		/// calculate the value of the continue child wavelet at time t
		/// @returns A T class value
		/// @param _t time 
		/// @param _t_shift shift time for child wavelet
		/// @param _scale scale parameter forchild wavelet
		////////////////////////////////////////////////////////////////////////////////////////////
		arma_inline T child(double _t, double _t_shit, double _scale) {
			set_child(_t_shit, _scale);
			return child(_t);
		}

		////////////////////////////////////////////////////////////////////////////////////////////
		/// \brief general description of a child wavelet using setted scale and time shift parameters
		///
		/// calculate the value of the continue child wavelet at time t
		/// @returns A T class value
		/// @param _t time 
		///////////////////////////////////////////////////////////////////////////////////////////
		arma_inline T child(double _t) {
			return (1 / std::sqrt(m_scale)) * mother((_t - m_t_shift) / m_scale);
		}

		////////////////////////////////////////////////////////////////////////////////////////////
		/// \brief setup of shit time and scale parameters
		///
		/// setup of parameters to avoid useless calculations over and over and over....
		/// @returns A T class value
		/// @param _t time 
		///////////////////////////////////////////////////////////////////////////////////////////
		void set_child(double _t_shift, double _scale) {
			if (_scale == double(0.0)) {
				throw std::invalid_argument("cannot set a scale parameter to 0");
			}
			m_scale = _scale;
			m_t_shift = _t_shift;
		}

	private:
		double m_t_shift; ///! time shift for child wavelet
		double m_scale; ///! scale value for child wavelet

		virtual T psi(double _arg) = 0;

	};

	class wavelet_haar : public Wavelet<double> 
	{
	private:
		double psi(double _t) 
		{
			if (_t >= double(0.0) && _t < double(0.5)) return double(1.0);
			if (_t >= double(0.5) && _t < double(1.0)) return double(-1.0);
			return double(0.0);
		}
	};

	class wavelet_shanon : public Wavelet<double> 
	{
	private:
		double psi(double _t) 
		{
			return 2.0 * sp::sinc( 2.0 * _t - 1 ) - sp::sinc( _t );
		}
	};

	class wavelet_mexican_hat : public Wavelet<double>
	{
	protected:
		double psi(double _t)
		{
			return (2/std::sqrt(3)*pow(sp::PI,-1/4))*(1 - pow(_t, 2))*arma::trunc_exp(-pow(_t, 2) / 2);
		}
	};

	class wavelet_morlet : public Wavelet<double>
	{
	private:
		double psi(double _t)
		{
			arma::trunc_exp(-pow(_t, 2) / 2)*std::cos(5 * _t);
		}
	};


	/// @}
}// end namespace

#endif