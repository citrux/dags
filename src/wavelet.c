#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdbool.h>
#include "wavelet.h"

#define MAX_BUFFER_SIZE 256 // must be 2^n
#define N_FREQS 85

#define wavelet_f0      6
#define wavelet_dofmin  2
#define wavelet_cdelta  0.776
#define wavelet_gamma   2.32
#define wavelet_deltaj0 0.6
#define wavelet_flambda ((4 * M_PI) / (wavelet_f0 + sqrtf(2.0 + wavelet_f0 * wavelet_f0)))

float wavelet_psi_ft(float f) {
    return powf(M_PI, -0.25) * expf(-0.5 * powf(f - wavelet_f0, 2));
}

float complex BUFFER[MAX_BUFFER_SIZE];
float complex FFT[MAX_BUFFER_SIZE];
float complex FFT_CONV[MAX_BUFFER_SIZE];

size_t BUFFER_SIZE = 0;

/********** FFT ********/
void fft_bit_rev(float complex const *in,
                 float complex *out,
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
        float complex em = inverse ? cexpf(M_PI / m * I) : cexpf(-M_PI / m * I);
        for (size_t k = 0; k < n; k += 2 * m) {
            float complex e = 1;
            for (size_t j = 0; j < m; ++j) {
                float complex u = out[k + j];
                float complex v = e * out[k + j + m];
                out[k + j] = u + v;
                out[k + j + m] = u - v;
                e *= em;
            }
        }
    }
}


int fft(float complex const *in, float complex *out, size_t n) {
    if (n & (n - 1)) {
        return 1;
    }
    fft_bit_rev(in, out, n, false);
    return 0;
}

int ifft(float complex const *in, float complex *out, size_t n) {
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
float ar1(const float complex *signal, size_t n) {
    // Estimates the lag zero and one covariance
    float c0 = 0, c1 = 0;
    for (size_t i = 0; i < n - 1; ++i)
    {
        c0 += signal[i] * signal[i];
        c1 += signal[i] * signal[i+1];
    }
    c0 += signal[n-1] * signal[n-1];
    c0 /= n;
    c1 /= (n - 1);

    // According to A. Grinsteds' substitutions
    float A = c0 * n * n;
    float B = -c1 * n - c0 * n * n - 2 * c0 + 2 * c1 - c1 * n * n + c0 * n;
    float C = n * (c0 + c1 * n - c1);
    float D = B * B - 4 * A * C;

    return (-B - sqrtf(D)) / (2 * A);
}


// A chi-squared percent point function for 95%
float chi2_ppf_95(float dof) {
    if (dof < 20) {
        return 2.10308 + 1.98911 * dof - 0.0464057 * powf(dof, 2)
             + 0.00102171 * powf(dof, 3);
    }
    return 9.27863 + 1.157 * dof;
}

float ftfreq(size_t i, size_t n) {
    if (2 * i < n) {
        return 2 * M_PI * i / n;
    }
    return -2 * M_PI * (n - i) / n;
}

void remove_trend(float complex *signal, size_t n) {
    float sumx = (n - 1) * n / 2,
          sumx2 = n * (n - 1) * (2 * n - 1) / 6,
          sumy=0, sumxy=0;
    for (size_t i = 0; i < n; ++i) {
        printf("%f %f\n", crealf(signal[i]), sumy);
        sumy += crealf(signal[i]);
        sumxy += i * crealf(signal[i]);
    }

    float a = (n * sumxy - sumx * sumy) / (n * sumx2 - sumx * sumx);
    float b = -(sumx * sumxy - sumx2 * sumy) / (n * sumx2 - sumx * sumx);

    printf("%f %f %f %f %f %f\n", a, b, sumx, sumx2, sumy, sumxy);

    for (size_t i = 0; i < n; ++i) {
        signal[i] -= a * i + b;
    }
}

uint16_t steps() {
    remove_trend(BUFFER, BUFFER_SIZE);
    /* calculate autocorrelation */
    float alpha = ar1(BUFFER, BUFFER_SIZE);

    /* calculate variance */
    float mean = 0, m2 = 0;
    for (size_t i = 0; i < BUFFER_SIZE; ++i) {
        float delta = crealf(BUFFER[i]) - mean;
        mean += delta / (i + 1);
        float delta2 = crealf(BUFFER[i]) - mean;
        m2 += delta * delta2;
    }
    float variance = m2 / BUFFER_SIZE;
    printf("%f %f\n", variance, m2);


    /* scale step */
    float dj = 1.0 / 12;
    /* Minimal wavelet scale */
    float s0 = 2 / wavelet_flambda;

    for (size_t i = BUFFER_SIZE; i < MAX_BUFFER_SIZE; ++i) {
        BUFFER[i] = 0;
    }

    /* perform FT */
    fft(BUFFER, FFT, MAX_BUFFER_SIZE);

    float max_power = 0;
    float max_freq = 0;
    float max_signif = 0;
    float dofmin = wavelet_dofmin;     // Degrees of freedom with no smoothing
    float gamma_fac = wavelet_gamma;   // Time-decorrelation factor
    
    for (size_t i = 0; i < N_FREQS; ++i) {
        float scale = s0 * powf(2.0, i * dj);
        float freq = 1 / (wavelet_flambda * scale);
        for (size_t k = 0; k < MAX_BUFFER_SIZE; ++k) {
            float psi_ft_bar = sqrtf(scale
                                   * ftfreq(1, MAX_BUFFER_SIZE)
                                   * MAX_BUFFER_SIZE)
                             * wavelet_psi_ft(scale
                                            * ftfreq(k, MAX_BUFFER_SIZE));
            FFT_CONV[k] = FFT[k] * psi_ft_bar;
        }

        ifft(FFT_CONV, BUFFER, MAX_BUFFER_SIZE);

        float power = 0;
        for (size_t k = 0; k < BUFFER_SIZE; ++k) {
            power += powf(cabsf(BUFFER[k]), 2);
        }

        if (power > max_power) {
            max_power = power;
            max_freq = freq;

            float dof = BUFFER_SIZE - scale;
            if (dof < 1) { dof = 1; }
            dof = dofmin * sqrtf(1 + powf(dof / gamma_fac / scale, 2));
            if (dof < dofmin) { dof = dofmin; }
            float chisquare = chi2_ppf_95(dof) / dof;
            float fft_theor = variance * (1 - powf(alpha, 2)) /
                (1 + powf(alpha, 2) -
                 2 * alpha * cos(2 * M_PI * freq / BUFFER_SIZE));
            max_signif = fft_theor * chisquare;
        }
    }

    if (max_power > max_signif) {
        return max_freq * BUFFER_SIZE;
    }

    return 0;
}

void waveletProcessNewData(uint16_t x, uint16_t y, uint16_t z, uint16_t time) {
    if (BUFFER_SIZE < MAX_BUFFER_SIZE) {
        BUFFER[BUFFER_SIZE++] = sqrtf(x * x + y * y + z * z);
    }
}

uint16_t waveletGetStepsCount(void) {
    uint16_t count = steps();
    BUFFER_SIZE = 0;
    return count;
}