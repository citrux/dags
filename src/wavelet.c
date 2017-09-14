#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdbool.h>
#include "wavelet.h"

#define wavelet_f0      6
#define wavelet_dofmin  2
#define wavelet_cdelta  0.776
#define wavelet_gamma   2.32
#define wavelet_deltaj0 0.6
#define wavelet_flambda ((4 * M_PI) / (wavelet_f0 + sqrt(2.0 + wavelet_f0 * wavelet_f0)))

double wavelet_psi_ft(double f) {
    return pow(M_PI, -0.25) * exp(-0.5 * pow(f - wavelet_f0, 2));
}

/********** FFT ********/
void fft_bit_rev(double complex const *in,
                 double complex *out,
                 size_t n,
                 bool inverse) {

    size_t bit_reverse(size_t i, size_t n) {
        size_t result = 0;
        while (n) {
            result <<= 1;
            result |= i & 1;
            i >>= 1;
            n >>= 1;
        }
        return result;
    }

    for (size_t i = 0; i < n; ++i) {
        out[i] = in[bit_reverse(i, n-1)];
    }
    for (size_t m = 1; m < n; m *= 2) {
        double complex em = inverse ? cexp(M_PI / m * I) : cexp(-M_PI / m * I);
        for (size_t k = 0; k < n; k += 2 * m) {
            double complex e = 1;
            for (size_t j = 0; j < m; ++j) {
                double complex u = out[k + j];
                double complex v = e * out[k + j + m];
                out[k + j] = u + v;
                out[k + j + m] = u - v;
                e *= em;
            }
        }
    }
}


int fft(double complex const *in, double complex *out, size_t n) {
    if (n & (n - 1)) {
        return 1;
    }
    fft_bit_rev(in, out, n, false);
    return 0;
}

int ifft(double complex const *in, double complex *out, size_t n) {
    if (n & (n - 1)) {
        return 1;
    }
    fft_bit_rev(in, out, n, true);
    for (size_t i = 0; i < n; ++i) {
        out[i] /= n;
    }
    return 0;
}

/********* FFT ********/

// Estimate of the lag-one autocorrelation.
double ar1(const double *signal, int n) {
    // Estimates the lag zero and one covariance
    double c0 = 0, c1 = 0;
    for (int i = 0; i < n - 1; ++i)
    {
        c0 += signal[i] * signal[i];
        c1 += signal[i] * signal[i+1];
    }
    c0 += signal[n-1] * signal[n-1];
    c0 /= n;
    c1 /= (n - 1);

    // According to A. Grinsteds' substitutions
    double A = c0 * n * n;
    double B = -c1 * n - c0 * n * n - 2 * c0 + 2 * c1 - c1 * n * n + c0 * n;
    double C = n * (c0 + c1 * n - c1);
    double D = B * B - 4 * A * C;
    if (D < 0) { return 0.9; }

    return (-B - sqrt(D)) / (2 * A);
}

void remove_trend(double *signal, size_t n) {
    double sumx = (n - 1) * n / 2,
          sumx2 = n * (n - 1) * (2 * n - 1) / 6,
          sumy=0, sumxy=0;
    for (size_t i = 0; i < n; ++i) {
        sumy += crealf(signal[i]);
        sumxy += i * crealf(signal[i]);
    }

    double a = (n * sumxy - sumx * sumy) / (n * sumx2 - sumx * sumx);
    double b = -(sumx * sumxy - sumx2 * sumy) / (n * sumx2 - sumx * sumx);

    for (size_t i = 0; i < n; ++i) {
        signal[i] -= a * i + b;
    }
}

double chi2_ppf_95(double dof) {
    if (dof < 20) {
        return 2.10308 + 1.98911 * dof - 0.0464057 * powf(dof, 2)
             + 0.00102171 * powf(dof, 3);
    }
    return 9.27863 + 1.157 * dof;
}

void wavelet_power(double *signal, int n, double dt, double *power,
                   double *frequencies, double *signif, int m) {
    /* scale step */
    double dj = 1.0 / 12;
    /* Minimal wavelet scale */
    double s0 = 2 * dt / wavelet_flambda;

    remove_trend(signal, n);
    /* Find the nearest power of 2, greater than n */
    int N = 1;
    while (N < n) { N *= 2; }

    double complex image[N];
    double complex image_conv[N];
    double complex preimage[N];

    /* convert double to double complex */
    for (size_t i = 0; i < n; ++i) {
        preimage[i] = signal[i];
    }
    for (size_t i = n; i < N; ++i) {
        preimage[i] = 0;
    }

    /* perform FT */
    fft(preimage, image, N);

    double ftfrequencies[N];
    for (int i = 0; i < N; ++i) {
        if (2 * i < N) {
            ftfrequencies[i] = 2 * M_PI * i / dt / N;
        } else {
            ftfrequencies[i] = 2 * M_PI * (i - N) / dt / N;
        }
    }

    double scales[m];
    for (size_t i = 0; i < m; ++i) {
        scales[i] = s0 * pow(2.0, i * dj);
        frequencies[i] = 1 / (wavelet_flambda * scales[i]);
    }
    for (size_t i = 0; i < m; ++i) {
        for (size_t k = 0; k < N; ++k) {
            double psi_ft_bar = sqrt(scales[i] * ftfrequencies[1] * N) *
                                  wavelet_psi_ft(scales[i] * ftfrequencies[k]);
            image_conv[k] = image[k] * psi_ft_bar;
        }

        ifft(image_conv, preimage, N);

        power[i] = 0;
        for (size_t k = 0; k < n; ++k) {
            power[i] += pow(cabs(preimage[k]), 2);
        }
        power[i] /= n;
    }

    /* calculate autocorrelation */
    double alpha = ar1(signal, n);

    /* now calculate significance */
    double mean = 0, m2 = 0;
    for (int i = 0; i < n; ++i) {
        double delta = signal[i] - mean;
        mean += delta / (i + 1);
        double delta2 = signal[i] - mean;
        m2 += delta * delta2;
    }
    double variance = m2 / n;

    double dofmin = wavelet_dofmin;     // Degrees of freedom with no smoothing
    double Cdelta = wavelet_cdelta;     // Reconstruction factor
    double gamma_fac = wavelet_gamma;   // Time-decorrelation factor
    double dj0 = wavelet_deltaj0;       // Scale-decorrelation factor

    /* Time-averaged significance*/
    for (int i = 0; i < m; ++i) {
        double dof = n; //- scales[i];
        if (dof < 1) { dof = 1; }
        dof = dofmin * sqrt(1 + pow(dof * dt / gamma_fac / scales[i], 2));
        if (dof < dofmin) { dof = dofmin; }
        double chisquare = chi2_ppf_95(dof) / dof;
        double fft_theor = variance * (1 - pow(alpha, 2)) /
            (1 + pow(alpha, 2) -
             2 * alpha * cos(2 * M_PI * frequencies[i] * dt / n));
        signif[i] = fft_theor * chisquare;
    }
}
