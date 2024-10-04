#pragma once
#include "matrix_vecor_mapper.hpp"
#include <future>

#define complex_mult_real(a, b) a.real() * b.real() - a.imag() * b.imag()
#define complex_mult_imag(a, b) a.real() * b.imag() + a.imag() * b.real()
#define complex_conj_mult_real(a, b) a.real() * b.real() + a.imag() * b.imag() //b is conjugated
#define complex_conj_mult_imag(a, b) a.imag() * b.real() - a.real() * b.imag() //b is conjugated

#define NUM_DFT_THREADS 8
#define PREALLOC_SIZE 8192 // 512 complex_t

namespace YADRO_TEST {
	class FFT {
	public:
		enum class dft_direction {
			forward,
			inverse
		};
	private:

#define DFT_METHOD template<dft_direction direction> static void
#define DFT_METHOD_DEF template<FFT::dft_direction direction> static void

		static void fill_vector_map(const complex_t* src_ptr, complex_t* dst_ptr, const matrix_vector_mapper& src_map, matrix_vector_mapper& dst_map) {
			for (int i = 0; i < src_map.count(); i++)
				for (int j = 0; j < src_map.length(); j++)
					dst_map.get(dst_ptr, i, j) = src_map.get(src_ptr, i, j);
		}

		//takes maps pointers because std::async doenst allow references
		DFT_METHOD fill_dft_separate(const complex_t* src_ptr, complex_t* dst_ptr, const matrix_vector_mapper* src_map, matrix_vector_mapper* dst_map, uint32_t rowInd) {
			
			const number_t m2pi_over_N = -d_pi / src_map->length();
			for (size_t f = 0; f < src_map->length(); f++) {
				complex_t& elem = dst_map->get(dst_ptr, rowInd, f);
				for (size_t n = 0; n < src_map->length(); n++) {
					const number_t w = n * f * m2pi_over_N;
					const complex_t e = { cos_f(w), sin_f(w) };

					if constexpr (direction == dft_direction::forward) {
						elem._Val[0] += complex_mult_real(src_map->get(src_ptr, rowInd, n), e);
						elem._Val[1] += complex_mult_imag(src_map->get(src_ptr, rowInd, n), e);
					}
					else {
						elem._Val[0] += complex_conj_mult_real(src_map->get(src_ptr, rowInd, n), e);
						elem._Val[1] += complex_conj_mult_imag(src_map->get(src_ptr, rowInd, n), e);
					}
				}
				if constexpr (direction == dft_direction::inverse) {
					elem._Val[0] = elem._Val[0] / src_map->length();
					elem._Val[1] = -elem._Val[1] / src_map->length();
				}
			}
		}

		//src map and dst map must "live" to the end of the function
		DFT_METHOD mapped_dft(const complex_t* src_ptr, complex_t* dst_ptr, const matrix_vector_mapper& src_map, matrix_vector_mapper& dst_map) {
			std::future<void> futures[NUM_DFT_THREADS];
			int32_t current_row = src_map.count() - 1;

			while (current_row >= 0) {
				//all threads should take about same time, since they doing the same job for same amount of data
				//then no need to use every vacant thread upon completion just wait while all finish
				uint8_t used_threads = 0;
				for (int i = 0; i < std::min(current_row + 1, NUM_DFT_THREADS); i++) {
					futures[i] = std::async(std::launch::async, fill_dft_separate<direction>, src_ptr, dst_ptr, &src_map, &dst_map, current_row);
					current_row -= 1;
					used_threads++;
				}

				for (int i = 0; i < used_threads; i++) {
					futures[i].wait();
				}
			}
		}
	public:
		//It is assumed that memory for the arrays is preallocated
		DFT_METHOD dft(const complex_t* const src_ptr, complex_t* const dst_ptr, const size_t N);
		DFT_METHOD fft_subdivided(const complex_t* const src_ptr, complex_t* const dst_ptr, const size_t N, const size_t N1);
		static void fft_recursive(const complex_t* const src_ptr, complex_t* const dst_ptr, const size_t N);
		static void ifft_recursive(const complex_t* const src_ptr, complex_t* const dst_ptr, const size_t N);

		DFT_METHOD dft(const std::vector<complex_t>& points, std::vector<complex_t>& spectrum);
		DFT_METHOD fft_subdivided(const std::vector<complex_t>& points, std::vector<complex_t>& spectrum, const size_t N1);
		static void fft_recursive(const std::vector<complex_t>& points, std::vector<complex_t>& spectrum);
		static void ifft_recursive(const std::vector<complex_t>& points, std::vector<complex_t>& spectrum);

	};


	DFT_METHOD_DEF FFT::dft(const complex_t* const src_ptr, complex_t* const dst_ptr, const size_t N) {
		if (src_ptr == nullptr || dst_ptr == nullptr)
			throw std::exception("null pointer received");

		const number_t m_pi_over_n = -d_pi / N;

		for (uint32_t f = 0; f < N; f++) {
			complex_t& elem = dst_ptr[f];
			for (uint32_t n = 0; n < N; n++) {
				const number_t w = m_pi_over_n * n * f;
				const complex_t e = { cos_f(w), sin_f(w) };

				if constexpr (direction == dft_direction::inverse) {
					elem._Val[0] += complex_conj_mult_real(src_ptr[n], e);
					elem._Val[1] += complex_conj_mult_imag(src_ptr[n], e);
				}
				else {
					elem._Val[0] += complex_mult_real(src_ptr[n], e);
					elem._Val[1] += complex_mult_imag(src_ptr[n], e);
				}
			}
			if constexpr (direction == dft_direction::inverse) {
				elem._Val[0] = dst_ptr[f]._Val[0] / N;
				elem._Val[1] = dst_ptr[f]._Val[1] / N;
			}
		}
	}
	DFT_METHOD_DEF FFT::dft(const std::vector<complex_t>& points, std::vector<complex_t>& spectrum) {
		spectrum.resize(points.size());
		dft<direction>(points.data(), spectrum.data(), points.size());
	}

	DFT_METHOD_DEF FFT::fft_subdivided(const complex_t* const points, complex_t* const dst_ptr, const size_t N, const size_t N1) {
		if (points == nullptr || dst_ptr == nullptr)
			throw std::exception("null pointer received");
		if (N % N1 != 0)
			throw std::exception("number of points must be divisible by N1");

		const uint32_t N2 = N / N1;

		//complex_t* temp = new complex_t[N];
		complex_t* temp = FFT::alloc_temp(N);

		const matrix_vector_mapper pointsMap(N, N1, N2, vector_map_direction::crosswise);
		matrix_vector_mapper tempMap(N, N1, N2, vector_map_direction::lengthwise);
		matrix_vector_mapper resultMap(N, N2, N1, vector_map_direction::crosswise);

		mapped_dft<direction>(points, temp, pointsMap, tempMap);

		tempMap.swap_dimensions();
		tempMap.change_direction();

		const number_t m2pi_over_N = -d_pi / N;
		complex_t temp_val;

		for (uint32_t el = 0; el < N1; el++) {
			for (uint32_t seq = 0; seq < N2; seq++) {
				const number_t w = m2pi_over_N * seq * el;
				const complex_t e = { cos_f(w), sin_f(w) };

				temp_val = tempMap.get(temp, seq, el);

				tempMap.get(temp, seq, el)._Val[0] = complex_mult_real(temp_val, e);
				tempMap.get(temp, seq, el)._Val[1] = complex_mult_imag(temp_val, e);
			}
		}

		mapped_dft<FFT::dft_direction::forward>(temp, dst_ptr, tempMap, resultMap);

		if constexpr (direction == dft_direction::inverse) {
			const number_t k = 1.0 / N * N2;
			for (size_t p = 0; p < N; p++) {
				dst_ptr[p]._Val[0] = dst_ptr[p]._Val[0] * k;
				dst_ptr[p]._Val[1] = -dst_ptr[p]._Val[1] * k;
			}

		}
		//delete[] temp;
		FFT::reset_prealloc();
	}
	DFT_METHOD_DEF FFT::fft_subdivided(const std::vector<complex_t>& points, std::vector<complex_t>& spectrum, const size_t N1) {
		if (points.size() % N1 != 0)
			throw std::exception("number of points must be divisible by N1");
		spectrum.resize(points.size());

		fft_subdivided<direction>(points.data(), spectrum.data(), points.size(), N1);
	}

	void FFT::fft_recursive(const complex_t* const src_ptr, complex_t* const dst_ptr, const size_t N) {
		if (src_ptr == nullptr || dst_ptr == nullptr)
			throw std::exception("null pointer received");

		size_t N1 = 1;

		if (N % 2 == 0) {
			N1 = 2;
		}
		else if (N % 3 == 0) {
			N1 = 3;
		}
		else if (N % 4 == 0) {
			N1 = 4;
		}
		else if (N % 5 == 0) {
			N1 = 5;
		}
		else if (N % 7 == 0) {
			N1 = 7;
		}
		else if (N % 8 == 0) {
			N1 = 8;
		}
		else if (N % 16 == 0) {
			N1 = 16;
		}
		else if (N % 32 == 0) {
			N1 = 32;
		}
		

		if (N1 == 1) {//that means N is not divisible by 2|3|5|...
			dft<dft_direction::forward>(src_ptr, dst_ptr, N);
		}
		else {
			const size_t N2 = N / N1;

			if (N2 > 1) {//that means N = 2|3|5|... * something

				std::unique_ptr<complex_t[]> temp = std::make_unique<complex_t[]>(N);
				
				const matrix_vector_mapper pointsMap(N, N2, N1, vector_map_direction::crosswise);
				matrix_vector_mapper tempMap(N, N2, N1, vector_map_direction::lengthwise);
				matrix_vector_mapper resultMap(N, N1, N2, vector_map_direction::crosswise);

				FFT::fill_vector_map(src_ptr, temp.get(), pointsMap, tempMap);

				for (int i = 0; i < N2; i++) {
					fft_recursive(temp.get() + N1 * i, dst_ptr + N1 * i, N1);
				}

				const number_t m2pi_over_N = -d_pi / N;
				complex_t temp_val;

				for (uint32_t el = 0; el < N2; el++) {
					for (uint32_t seq = 0; seq < N1; seq++) {
						const number_t w = m2pi_over_N * seq * el;
						const complex_t e = { cos_f(w), sin_f(w) };

						temp_val = resultMap.get(dst_ptr, seq, el);

						resultMap.get(dst_ptr, seq, el)._Val[0] = complex_mult_real(temp_val, e);
						resultMap.get(dst_ptr, seq, el)._Val[1] = complex_mult_imag(temp_val, e);
					}
				}

				temp.reset(new complex_t[N]);
				mapped_dft<dft_direction::forward>(dst_ptr, temp.get(), resultMap, resultMap);
				memcpy(dst_ptr, temp.get(), sizeof(complex_t) * N);

			}
			else {//that means N1 = N = 2|3|5|...
				dft<dft_direction::forward>(src_ptr, dst_ptr, N);
			}
		}
	}
	void FFT::fft_recursive(const std::vector<complex_t>& points, std::vector<complex_t>& spectrum) {
		spectrum.resize(points.size());
		fft_recursive(points.data(), spectrum.data(), points.size());
	}

	void FFT::ifft_recursive(const complex_t* const src_ptr, complex_t* const dst_ptr, const size_t N) {
		const number_t norm = 1.0 / N;

		complex_t* tempa = new complex_t[N];

		for (uint32_t i = 0; i < N; i++) {
			tempa[i] = std::conj(src_ptr[i]);
		}

		FFT::fft_recursive(tempa, dst_ptr, N);

		for (uint32_t i = 0; i < N; i++) {
			dst_ptr[i] = std::conj(dst_ptr[i]) * norm;
		}

		delete[] tempa;
	}
	void FFT::ifft_recursive(const std::vector<complex_t>& points, std::vector<complex_t>& spectrum) {
		spectrum.resize(points.size());
		ifft_recursive(points.data(), spectrum.data(), points.size());
	}
}
	

#undef DFT_METHOD
#undef DFT_METHOD_DEF
#undef complex_mult_real
#undef complex_mult_imag
#undef complex_conj_mult_real
#undef complex_conj_mult_imag