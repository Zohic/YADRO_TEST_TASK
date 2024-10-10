#include "fft.hpp"
#include <chrono>

#define COUNT_TIME_BEGIN(TIMER) TIMER = std::chrono::system_clock::now();
#define COUNT_TIME_END(TIMER, RES) \
RES = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - TIMER).count();\

#define TEST_TIME(TIMER, RES, FUNC)\
COUNT_TIME_BEGIN(TIMER)\
FUNC;\
COUNT_TIME_END(TIMER, RES)

#define TEST_TIME_PRINT(TEST_NAME, TIMER, RES, FUNC)\
COUNT_TIME_BEGIN(TIMER)\
FUNC;\
COUNT_TIME_END(TIMER, RES)\
printf(TEST_NAME " took %zu milliseconds\n", RES);



int main(int argc, char** argv) {
	using namespace YADRO_TEST;
	using dft_dir = FFT::dft_direction;
	
	constexpr size_t pts_cnt = 3 * 5 * 7 * 13;
	//constexpr int pts_cnt = 3 * 5 * 7 * 13 * 10;
	std::vector<complex_t> points(pts_cnt);

	for (int i = 0; i < pts_cnt; i++) {
		points[i]._Val[0] = -5 + (number_t)rand() / RAND_MAX * 10;
		points[i]._Val[1] = -5 + (number_t)rand() / RAND_MAX * 10;
	}

	std::vector<complex_t> DFTspectrum(points.size());
	std::vector<complex_t> IDFTspectrum(points.size());
	std::vector<complex_t> FFTSspectrum(points.size());
	std::vector<complex_t> IFFTSspectrum(points.size());
	std::vector<complex_t> FFTspectrum(points.size());
	std::vector<complex_t> IFFTspectrum(points.size());

	std::chrono::time_point t1 = std::chrono::system_clock::now();
	size_t passed_ms[6] = { 0 };
	size_t total_time[6] = { 0 };

	const uint32_t tests = 10;

	//TEST_TIME_PRINT("fft", t1, passed_ms[4], FFT::fft_recursive<dft_dir::forward>(points, FFTspectrum));
	//TEST_TIME_PRINT("ifft", t1, passed_ms[5], FFT::fft_recursive<dft_dir::inverse>(FFTspectrum, IFFTspectrum));

	for (int t = 0; t < tests; t++) {
		TEST_TIME(t1, passed_ms[0], FFT::dft<dft_dir::forward>(points, DFTspectrum));
		TEST_TIME(t1, passed_ms[1], FFT::dft<dft_dir::inverse>(DFTspectrum, IDFTspectrum));
		TEST_TIME(t1, passed_ms[2], FFT::fft_subdivided<dft_dir::forward>(points, FFTSspectrum, 13));
		TEST_TIME(t1, passed_ms[3], FFT::fft_subdivided<dft_dir::inverse>(FFTSspectrum, IFFTSspectrum, 13));
		TEST_TIME(t1, passed_ms[4], FFT::fft_recursive<dft_dir::forward>(points, FFTspectrum));
		TEST_TIME(t1, passed_ms[5], FFT::fft_recursive<dft_dir::inverse>(FFTspectrum, IFFTspectrum));
		printf("done test %i\n", t);
		for (int i = 0; i < 6; i++)
			total_time[i] += passed_ms[i];
	}

	for (int i = 0; i < 6; i++)
		printf("avertage time %i = %10.20f ms\n", i, (double)total_time[i] / tests);
}
