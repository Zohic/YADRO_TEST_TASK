#pragma once
#include "matrix_vecor_mapper.hpp"

#define complex_mult_real(a, b) a.real() * b.real() - a.imag() * b.imag()
#define complex_mult_imag(a, b) a.real() * b.imag() + a.imag() * b.real()
#define complex_conj_mult_real(a, b) a.real() * b.real() + a.imag() * b.imag() //b is conjugated
#define complex_conj_mult_imag(a, b) a.imag() * b.real() - a.real() * b.imag() //b is conjugated

inline constexpr uint8_t NUM_DFT_THREADS = 4;
inline constexpr uint32_t PREALLOC_SIZE = 1 << 15; // of complex<double>


#define ARRAY_DFT_ASYNC 0
#define ARRAY_DFT_THREADS 1

#define ARRAY_DFT_METHOD ARRAY_DFT_THREADS

#if ARRAY_DFT_METHOD == ARRAY_DFT_ASYNC
#include <future>
#elif ARRAY_DFT_METHOD == ARRAY_DFT_THREADS
#include <mutex>
#include <atomic>
#include <thread>
#elif
static_assert(false, "ARRAY_DFT_METHOD must be set to one of options");
#endif


namespace YADRO_TEST {
	class FFT {
	public:
		enum class dft_direction {
			forward,
			inverse
		};

	private:

#define DFT_METHOD template<FFT::dft_direction direction> static void

#if ARRAY_DFT_METHOD == ARRAY_DFT_ASYNC
		//takes maps as pointers because std::async doenst allow references
		DFT_METHOD fill_dft_separate(const complex_t* src_ptr, complex_t* dst_ptr, const matrix_vector_mapper* src_map, matrix_vector_mapper* dst_map, uint32_t row_ind);
#elif ARRAY_DFT_METHOD == ARRAY_DFT_THREADS
		template<uint32_t num_threads> struct thread_queue {
			std::thread _threads[num_threads];
			volatile std::atomic_uint8_t _job_counter = num_threads;
			bool _alive = true;

			std::condition_variable _cv;
			std::mutex _m;

			const complex_t* _src_ptr;
			complex_t* _dst_ptr;
			const matrix_vector_mapper* _src_map;
			matrix_vector_mapper* _dst_map;

			uint32_t _per_thread_count;
			uint32_t _leftover;

			FFT::dft_direction _dir;
			
			template<FFT::dft_direction direction>
			void compute_parallel(const complex_t* src_ptr, complex_t* dst_ptr, const matrix_vector_mapper* src_map, matrix_vector_mapper* dst_map) {
				_src_ptr = src_ptr;
				_dst_ptr = dst_ptr;
				_src_map = src_map;
				_dst_map = dst_map;
				_dir = direction;

				_per_thread_count = src_map->count() / num_threads;
				_leftover = src_map->count() % num_threads;

				if (_job_counter.load() == num_threads || _job_counter.load() == num_threads * 2) {// num_threads * 2 for initiliazation
					_job_counter.store(0);
					_cv.notify_all();
					//printf("began computing...\n");
				}
				else
					throw std::exception("not all threads finished last time?");
				
				while (_job_counter.load() < num_threads) {
					//printf("now done %i\n", _job_counter.operator unsigned char());
				}
				//printf("finished computing!\n");
			}

			template<uint32_t ind = 0> void set_thread_func() {
				_threads[ind] = std::thread([this]() mutable {
					try {
						std::unique_lock u_lock(this->_m);
						while (this->_alive) {
							//printf("\tthread %zu waiting...\n", std::this_thread::get_id());
							this->_job_counter++;
							this->_cv.wait(u_lock);
							//printf("\tthread %zu working...\n", std::this_thread::get_id());
							
							if (!this->_alive)
								break;

							uint32_t leftover = 0;
							if constexpr (ind == num_threads - 1)
								leftover = this->_leftover;
							if (this->_dir == FFT::dft_direction::forward)
								fill_dft_some<FFT::dft_direction::forward>(
									this->_src_ptr, this->_dst_ptr,
									this->_src_map, this->_dst_map,
									this->_per_thread_count * ind, this->_per_thread_count * (ind + 1) + leftover);
							
							else
								fill_dft_some<FFT::dft_direction::inverse>(
									this->_src_ptr, this->_dst_ptr,
									this->_src_map, this->_dst_map,
									this->_per_thread_count * ind, this->_per_thread_count * (ind + 1) + leftover);
							//printf("\tthread %zu done!\n", std::this_thread::get_id());
						}
					}
					catch (std::exception& e) {
						printf("THREAD ERROR: %s\n", e.what());
					}
					});

				if constexpr (num_threads > 1 && ind < num_threads - 1)
					set_thread_func<ind + 1>();
			}

			thread_queue() {
				set_thread_func();
				std::this_thread::sleep_for(std::chrono::milliseconds(50));
			}
			~thread_queue() {
				_alive = false;
				_cv.notify_all();
				for (int i = 0; i < num_threads; i++)
					_threads[i].join();
			}
		};
		DFT_METHOD fill_dft_some(const complex_t* src_ptr, complex_t* dst_ptr, const matrix_vector_mapper* src_map, matrix_vector_mapper* dst_map, uint32_t from_row, uint32_t to_row);
		static thread_queue<NUM_DFT_THREADS> _thread_queue;
#endif // ARRAY_DFT_METHOD == ARRAY_DFT_THREADS

		struct preallocated {
			complex_t* _prealloc;
			size_t _prealloc_count;

			preallocated() {
				_prealloc = new complex_t[PREALLOC_SIZE];
				_prealloc_count = 0;
				memset(_prealloc, 0, sizeof(complex_t) * PREALLOC_SIZE);
			}

			complex_t* get(size_t N) {
				if (_prealloc_count + N >= PREALLOC_SIZE)
					throw std::exception("out of preallocated memory");

				//printf("alloceted %zu complex\n", N);

				size_t prev = _prealloc_count;
				_prealloc_count += N;
				return _prealloc + prev;
			}

			void reset() {
				//printf("reset, %zu was allocated\n", _prealloc_count);
				memset(_prealloc, 0, _prealloc_count * sizeof(complex_t));
				_prealloc_count = 0;
			}

		};
		static preallocated _prealloc;

		static void fill_vector_map(const complex_t* src_ptr, complex_t* dst_ptr, const matrix_vector_mapper& src_map, matrix_vector_mapper& dst_map);
		DFT_METHOD mapped_dft(const complex_t* src_ptr, complex_t* dst_ptr, const matrix_vector_mapper& src_map, matrix_vector_mapper& dst_map);
		static void fft_recursive_base(complex_t* src_ptr, complex_t* dst_ptr, const size_t N) {
			if (src_ptr == nullptr || dst_ptr == nullptr)
				throw std::exception("null pointer received");

			size_t N1 = 1;

			if (N % 32 == 0) {
				N1 = 32;
			}
			else if (N % 16 == 0) {
				N1 = 16;
			}
			else if (N % 8 == 0) {
				N1 = 8;
			}
			else if (N % 7 == 0) {
				N1 = 7;
			}
			else if (N % 5 == 0) {
				N1 = 5;
			}
			else if (N % 4 == 0) {
				N1 = 4;
			}
			else if (N % 3 == 0) {
				N1 = 3;
			}
			else if (N % 2 == 0) {
				N1 = 2;
			}


			if (N1 == 1) {//that means N is not divisible by 2|3|5|...
				dft<dft_direction::forward>(src_ptr, dst_ptr, N);
			}
			else {
				const size_t N2 = N / N1;

				if (N2 > 1) {//that means N = 2|3|5|... * something
					const matrix_vector_mapper pointsMap(N, N1, N2, vector_map_direction::crosswise);
					matrix_vector_mapper tempMap(N, N1, N2, vector_map_direction::lengthwise);
					matrix_vector_mapper resultMap(N, N2, N1, vector_map_direction::crosswise);


					//memset(dst_ptr, 0, sizeof(complex_t) * N);
					memcpy(dst_ptr, src_ptr, sizeof(complex_t) * N);
					fill_vector_map(dst_ptr, src_ptr, pointsMap, tempMap);
					memset(dst_ptr, 0, sizeof(complex_t) * N);
					for (int i = 0; i < N1; i++)
					{
						fft_recursive_base(src_ptr + i * N2, dst_ptr + i * N2, N2);
					}

					const number_t m2pi_over_N = -d_pi / N;
					complex_t temp_val;
					for (uint32_t el = 0; el < N1; el++) {
						for (uint32_t seq = 0; seq < N2; seq++) {
							const number_t w = m2pi_over_N * seq * el;
							const complex_t e = { cos_f(w), sin_f(w) };

							temp_val = resultMap.get(dst_ptr, seq, el);

							resultMap.get(dst_ptr, seq, el)._Val[0] = complex_mult_real(temp_val, e);
							resultMap.get(dst_ptr, seq, el)._Val[1] = complex_mult_imag(temp_val, e);
						}
					}

					memset(src_ptr, 0, sizeof(complex_t) * N);
					mapped_dft<dft_direction::forward>(dst_ptr, src_ptr, resultMap, resultMap);
					memcpy(dst_ptr, src_ptr, sizeof(complex_t) * N);
				}
				else {//that means N1 = N = 2|3|5|...
					dft<dft_direction::forward>(src_ptr, dst_ptr, N);
				}
			}
		}
	public:
		//It is assumed that memory for the arrays is preallocated and dst is zero initialized
		DFT_METHOD dft(const complex_t* const src_ptr, complex_t* const dst_ptr, const size_t N);
		//It is assumed that memory for the arrays is preallocated and dst is zero initialized
		DFT_METHOD fft_subdivided(const complex_t* const src_ptr, complex_t* const dst_ptr, const size_t N, const size_t N1);
		//It is assumed that memory for the arrays is preallocated and dst is zero initialized
		DFT_METHOD fft_recursive(complex_t* src_ptr, complex_t* dst_ptr, const size_t N);

		DFT_METHOD dft(const std::vector<complex_t>& points, std::vector<complex_t>& spectrum);
		DFT_METHOD fft_subdivided(const std::vector<complex_t>& points, std::vector<complex_t>& spectrum, const size_t N1);
		DFT_METHOD fft_recursive(std::vector<complex_t>& points, std::vector<complex_t>& spectrum);
	};

#if ARRAY_DFT_METHOD == ARRAY_DFT_ASYNC
	DFT_METHOD FFT::fill_dft_separate(const complex_t* src_ptr, complex_t* dst_ptr, const matrix_vector_mapper* src_map, matrix_vector_mapper* dst_map, uint32_t row_ind) {
		const number_t m2pi_over_N = -d_pi / src_map->length();
		for (size_t f = 0; f < src_map->length(); f++) {
			complex_t& elem = dst_map->get(dst_ptr, row_ind, f);
			for (size_t n = 0; n < src_map->length(); n++) {
				const number_t w = n * f * m2pi_over_N;
				const complex_t e = { cos_f(w), sin_f(w) };

				if constexpr (direction == FFT::dft_direction::forward) {
					elem._Val[0] += complex_mult_real(src_map->get(src_ptr, row_ind, n), e);
					elem._Val[1] += complex_mult_imag(src_map->get(src_ptr, row_ind, n), e);
				}
				else {
					elem._Val[0] += complex_conj_mult_real(src_map->get(src_ptr, row_ind, n), e);
					elem._Val[1] += complex_conj_mult_imag(src_map->get(src_ptr, row_ind, n), e);
				}
			}
			if constexpr (direction == FFT::dft_direction::inverse) {
				elem._Val[0] = elem._Val[0] / src_map->length();
				elem._Val[1] = -elem._Val[1] / src_map->length();
				}
			}
		}

	DFT_METHOD FFT::mapped_dft(const complex_t* src_ptr, complex_t* dst_ptr, const matrix_vector_mapper& src_map, matrix_vector_mapper& dst_map) {
		std::future<void> futures[NUM_DFT_THREADS];
		int32_t current_row = src_map.count() - 1;

		while (current_row >= 0) {
			//all threads should take about same time, since they doing the same job for same amount of data
			//then no need to use every vacant thread upon completion just wait while all finish
			uint8_t used_threads = 0;
			uint32_t a = current_row + 1;
			uint32_t b = (uint32_t)NUM_DFT_THREADS;
			uint32_t minab = std::min(a, b);
			for (int i = 0; i < minab; i++) {
				futures[i] = std::async(std::launch::async, FFT::fill_dft_separate<direction>, src_ptr, dst_ptr, &src_map, &dst_map, current_row);
				current_row -= 1;
				used_threads++;
			}

			for (int i = 0; i < used_threads; i++) {
				futures[i].wait();
			}
		}
	}
#elif ARRAY_DFT_METHOD == ARRAY_DFT_THREADS
	DFT_METHOD FFT::fill_dft_some(const complex_t* src_ptr, complex_t* dst_ptr, const matrix_vector_mapper* src_map, matrix_vector_mapper* dst_map, uint32_t from_row, uint32_t to_row) {
		const number_t m2pi_over_N = -d_pi / src_map->length();

		for (uint32_t c = from_row; c < to_row; c++)
			for (size_t f = 0; f < src_map->length(); f++) {

				complex_t& elem = dst_map->get(dst_ptr, c, f);

				for (size_t n = 0; n < src_map->length(); n++) {
					const number_t w = n * f * m2pi_over_N;
					const complex_t e = { cos_f(w), sin_f(w) };

					if constexpr (direction == FFT::dft_direction::forward) {
						elem._Val[0] += complex_mult_real(src_map->get(src_ptr, c, n), e);
						elem._Val[1] += complex_mult_imag(src_map->get(src_ptr, c, n), e);
					}
					else {
						elem._Val[0] += complex_conj_mult_real(src_map->get(src_ptr, c, n), e);
						elem._Val[1] += complex_conj_mult_imag(src_map->get(src_ptr, c, n), e);
					}
				}
				if constexpr (direction == FFT::dft_direction::inverse) {
					elem._Val[0] = elem._Val[0] / src_map->length();
					elem._Val[1] = -elem._Val[1] / src_map->length();
				}
			}
	}
	DFT_METHOD FFT::mapped_dft(const complex_t * src_ptr, complex_t * dst_ptr, const matrix_vector_mapper & src_map, matrix_vector_mapper & dst_map) {
		//synchronous operation. will stop current flow until finished
		FFT::_thread_queue.compute_parallel<direction>(src_ptr, dst_ptr, &src_map, &dst_map);
	}
#endif // ARRAY_DFT_METHOD
	
	void FFT::fill_vector_map(const complex_t* src_ptr, complex_t* dst_ptr, const matrix_vector_mapper& src_map, matrix_vector_mapper& dst_map) {
		for (int i = 0; i < src_map.count(); i++)
			for (int j = 0; j < src_map.length(); j++)
				dst_map.get(dst_ptr, i, j) = src_map.get(src_ptr, i, j);
	}

	DFT_METHOD FFT::dft(const complex_t* const src_ptr, complex_t* const dst_ptr, const size_t N) {
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
	DFT_METHOD FFT::dft(const std::vector<complex_t>& points, std::vector<complex_t>& spectrum) {
		spectrum.resize(points.size());
		dft<direction>(points.data(), spectrum.data(), points.size());
	}

	DFT_METHOD FFT::fft_subdivided(const complex_t* const points, complex_t* const dst_ptr, const size_t N, const size_t N1) {
		if (points == nullptr || dst_ptr == nullptr)
			throw std::exception("null pointer received");
		if (N % N1 != 0)
			throw std::exception("number of points must be divisible by N1");

		const uint32_t N2 = N / N1;

		//complex_t* temp = new complex_t[N];
		complex_t* temp = FFT::_prealloc.get(N);

		const matrix_vector_mapper pointsMap(N, N1, N2, vector_map_direction::crosswise);
		matrix_vector_mapper tempMap(        N, N1, N2, vector_map_direction::lengthwise);
		matrix_vector_mapper resultMap(      N, N2, N1, vector_map_direction::crosswise);

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
		FFT::_prealloc.reset();
	}
	DFT_METHOD FFT::fft_subdivided(const std::vector<complex_t>& points, std::vector<complex_t>& spectrum, const size_t N1) {
		if (points.size() % N1 != 0)
			throw std::exception("number of points must be divisible by N1");
		spectrum.resize(points.size());

		fft_subdivided<direction>(points.data(), spectrum.data(), points.size(), N1);
	}

	DFT_METHOD FFT::fft_recursive(complex_t* src_ptr, complex_t* dst_ptr, const size_t N) {
		if (src_ptr == nullptr || dst_ptr == nullptr)
			throw std::exception("null pointer received");

		if constexpr (direction == dft_direction::inverse) {
			complex_t* temp = FFT::_prealloc.get(N);

			for (int i = 0; i < N; i++) {
				temp[i] = std::conj(src_ptr[i]);
			}

			fft_recursive_base(temp, dst_ptr, N);

			const number_t inv_N = 1.0 / N;

			for (int i = 0; i < N; i++) {
				dst_ptr[i] = inv_N * std::conj(dst_ptr[i]);
			}

			return;
		}
		else {
			complex_t* temp = FFT::_prealloc.get(N);

			memcpy(temp, src_ptr, sizeof(complex_t) * N);
			fft_recursive_base(temp, dst_ptr, N);
		}
		FFT::_prealloc.reset();

	}
	DFT_METHOD FFT::fft_recursive(std::vector<complex_t>& points, std::vector<complex_t>& spectrum) {
		spectrum.resize(points.size());
		fft_recursive<direction>(points.data(), spectrum.data(), points.size());
	}
}
	


YADRO_TEST::FFT::preallocated YADRO_TEST::FFT::_prealloc = YADRO_TEST::FFT::preallocated();

#if ARRAY_DFT_METHOD == ARRAY_DFT_THREADS
	YADRO_TEST::FFT::thread_queue<NUM_DFT_THREADS> YADRO_TEST::FFT::_thread_queue = YADRO_TEST::FFT::thread_queue<NUM_DFT_THREADS>();
#endif //ARRAY_DFT_METHOD == ARRAY_DFT_THREADS

#undef DFT_METHOD
#undef complex_mult_real
#undef complex_mult_imag
#undef complex_conj_mult_real
#undef complex_conj_mult_imag