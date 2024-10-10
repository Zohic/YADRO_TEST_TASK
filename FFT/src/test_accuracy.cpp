#include "fft.hpp"
#include <iostream>
#include <algorithm>

int main(int argc, char** argv) {


	/*const number_t fmax = 100;
	const number_t fs = fmax * 2;
	const number_t Ts = 1.0 / fs;

	const number_t df = 4;
	const size_t pts_cnt = std::ceil(1.0 / (df * Ts));

	const number_t dw = 2 * pi * df;

	const number_t f = df * 3;
	const number_t w = 2 * pi * f;*/


	using namespace YADRO_TEST;
	using dft_dir = FFT::dft_direction;

	constexpr size_t pts_cnt = 3 * 5 * 7 * 13;
	constexpr int a = 3 * 5 * 7;
	std::vector<complex_t> points(pts_cnt);
	srand(2);
	for (int i = 0; i < pts_cnt; i++) {
		points[i]._Val[0] = -5 + (number_t)rand() / RAND_MAX * 10;
		points[i]._Val[1] = -5 + (number_t)rand() / RAND_MAX * 10;
	}

	std::vector<complex_t> DFTspectrum(pts_cnt);
	std::vector<complex_t> IDFTspectrum(pts_cnt);
	std::vector<complex_t> FFTSspectrum(pts_cnt);
	std::vector<complex_t> IFFTSspectrum(pts_cnt);
	std::vector<complex_t> FFTspectrum(pts_cnt);
	std::vector<complex_t> IFFTspectrum(pts_cnt);

	FFT::dft<dft_dir::forward>(points, DFTspectrum);
	FFT::dft<dft_dir::inverse>(DFTspectrum, IDFTspectrum);
	FFT::fft_subdivided<dft_dir::forward>(points, FFTSspectrum, 13);
	FFT::fft_subdivided<dft_dir::inverse>(FFTSspectrum, IFFTSspectrum, 13);
	FFT::fft_recursive<dft_dir::forward>(points, FFTspectrum);
	FFT::fft_recursive<dft_dir::inverse>(FFTspectrum, IFFTspectrum);

	std::vector<number_t> difference[3];
	difference[0].resize(pts_cnt);
	difference[1].resize(pts_cnt);
	difference[2].resize(pts_cnt);
	
	number_t max_err[3] = {0.0f, 0.0f , 0.0f };

	std::generate(difference[0].begin(), difference[0].end(), [&max_err, a = points.begin(), b = IDFTspectrum.begin()]() mutable {
		number_t val = std::abs(*a - *b);
		if (val > max_err[0])
			max_err[0] = val;
		a++;
		b++;
		return val;
		});
	std::generate(difference[1].begin(), difference[1].end(), [&max_err, a = points.begin(), b = IFFTSspectrum.begin()]() mutable {
		number_t val = std::abs(*a - *b);
		if (val > max_err[1])
			max_err[1] = val;
		a++;
		b++;
		return val;
	});
	std::generate(difference[2].begin(), difference[2].end(), [&max_err, a = points.begin(), b = IFFTspectrum.begin()]() mutable {
		number_t val = std::abs(*a - *b);
		if (val > max_err[2])
			max_err[2] = val;
		a++;
		b++;
		return val;
	});

	printf("max dft err: %10.25f\n", max_err[0]);
	printf("max subd fft err: %10.25f\n", max_err[1]);
	printf("max fft err: %10.25f\n", max_err[2]);
	
	number_t max_amp = std::abs(*std::max_element(points.begin(), points.end(), [](complex_t& a, complex_t& b) {return std::abs(a) < std::abs(b); }));

	printf("max abs of point = %10.10f\n", max_amp);

	number_t rel_err[3] = { max_err[0] / max_amp * 100, max_err[1] / max_amp * 100, max_err[2] / max_amp * 100 };

	printf("dft error: %5.20f%%\nfft_subdivide error: %5.20f%%\nfft error: %5.20f%%\n", rel_err[0], rel_err[1], rel_err[2]);
	printf("dft error: %5.20f\nfft_subdivide error: %5.20f\nfft error: %5.20f\n", 20*log10(rel_err[0]), 20 * log10(rel_err[1]), 20 * log10(rel_err[2]));
	
	for (int i = points.size() * 0; i < points.size(); i++) {
		printf("point: (%f, %f) | dft: (%f, %f) | fft: (%f, %f) | idft(%f, %f)",
			points[i].real(), points[i].imag(),
			DFTspectrum[i].real(), DFTspectrum[i].imag(),
			FFTspectrum[i].real(), FFTspectrum[i].imag(),
			IFFTspectrum[i].real(), IFFTspectrum[i].imag());

		std::cout << std::endl;
	}
		
}
