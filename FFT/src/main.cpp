#include <iostream>
#include <fstream>
#include "fft.hpp"
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

	size_t pts_cnt = 100;
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

	FFT::dft<dft_dir::forward>(points, DFTspectrum);
	FFT::dft<dft_dir::inverse>(DFTspectrum, IDFTspectrum);

	FFT::dft<dft_dir::forward>(points, FFTSspectrum);
	FFT::dft<dft_dir::inverse>(FFTSspectrum, IFFTSspectrum);

	FFT::fft_recursive(points, FFTspectrum);
	FFT::ifft_recursive(FFTspectrum, IFFTspectrum);

	
	number_t max_amp = 0.0;
	number_t max_err[3]{ 0.0, 0.0, 0.0 };
	for (int i = 0; i < pts_cnt; i++) {

		number_t ampl = std::abs(points[i]);
		number_t dif[3] = { std::abs(points[i] - IDFTspectrum[i]) , std::abs(points[i] - IFFTSspectrum[i]) , std::abs(points[i] - IFFTspectrum[i]) };
		
		if (ampl > max_amp)
			max_amp = ampl;

		if (dif[0] > max_err[0])
			max_err[0] = dif[0];

		if (dif[1] > max_err[1])
			max_err[1] = dif[1];

		if (dif[2] > max_err[2])
			max_err[2] = dif[2];
	}

	number_t rel_err[3] = { max_err[0] / max_amp * 100, max_err[1] / max_amp * 100, max_err[2] / max_amp * 100 };

	

	printf("dft error: %5.20f%%, fft_subdivide error: %5.20f%%, fft error: %5.20f%%", rel_err[0], rel_err[1], rel_err[2]);

	/*for (int i = 0; i < points.size(); i++) {
		printf("point: (%f, %f) | dft: (%f, %f) | idft(%f, %f) | diff(%10.20f)",
			points[i].real(), points[i].imag(),
			std::abs(DFTspectrum[i]), std::arg(DFTspectrum[i]),
			inverseDFT[i].real(), inverseDFT[i].imag(), difference[i]);

		std::cout << std::endl;
	}*/
		
}
