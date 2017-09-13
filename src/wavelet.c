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
#define wavelet_flambda ((4 * M_PI) / (wavelet_f0 + sqrtf(2.0 + wavelet_f0 * wavelet_f0)))

float wavelet_psi_ft(float f) {
    return powf(M_PI, -0.25) * expf(-0.5 * powf(f - wavelet_f0, 2));
}

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
float ar1(const float *signal, int n) {
    // Estimates the lag zero and one covariance
    float c0 = 0, c1 = 0;
    for (int i = 0; i < n - 1; ++i)
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

void wavelet_power(const float *signal, size_t n, float dt, float *power,
                   float *frequencies, float *signif, size_t m) {
    /* scale step */
    float dj = 1.0 / 12;
    /* Minimal wavelet scale */
    float s0 = 2 * dt / wavelet_flambda;

    /* Find the nearest power of 2, greater than n */
    size_t N = 1;
    while (N < n) { N *= 2; }

    float complex image[N];
    float complex image_conv[N];
    float complex preimage[N];

    /* convert float to float complex */
    for (size_t i = 0; i < n; ++i) {
        preimage[i] = signal[i];
    }
    for (size_t i = n; i < N; ++i) {
        preimage[i] = 0;
    }

    /* perform FT */
    fft(preimage, image, N);

    float ftfrequencies[N];
    for (size_t i = 0; i < N; ++i) {
        if (2 * i < N) {
            ftfrequencies[i] = 2 * M_PI * i / dt / N;
        } else {
            ftfrequencies[i] = 2 * M_PI * (i - N) / dt / N;
        }
    }

    float scales[m];
    for (size_t i = 0; i < m; ++i) {
        scales[i] = s0 * powf(2.0, i * dj);
        frequencies[i] = 1 / (wavelet_flambda * scales[i]);
    }
    for (size_t i = 0; i < m; ++i) {
        for (size_t k = 0; k < N; ++k) {
            float psi_ft_bar = sqrtf(scales[i] * ftfrequencies[1] * N) *
                                  wavelet_psi_ft(scales[i] * ftfrequencies[k]);
            image_conv[k] = image[k] * psi_ft_bar;
        }

        ifft(image_conv, preimage, N);

        power[i] = 0;
        for (size_t k = 0; k < n; ++k) {
            power[i] += powf(cabs(preimage[k]), 2);
        }
        power[i] /= n;
    }

    /* calculate autocorrelation */
    float alpha = ar1(signal, n);

    /* now calculate significance */
    float mean = 0, m2 = 0;
    for (size_t i = 0; i < n; ++i) {
        float delta = signal[i] - mean;
        mean += delta / (i + 1);
        float delta2 = signal[i] - mean;
        m2 += delta * delta2;
    }
    float variance = m2 / n;

    float dofmin = wavelet_dofmin;     // Degrees of freedom with no smoothing
    // float Cdelta = wavelet_cdelta;     // Reconstruction factor
    float gamma_fac = wavelet_gamma;   // Time-decorrelation factor
    // float dj0 = wavelet_deltaj0;       // Scale-decorrelation factor

    /* Time-averaged significance*/
    for (size_t i = 0; i < m; ++i) {
        float dof = n - scales[i];
        if (dof < 1) { dof = 1; }
        dof = dofmin * sqrtf(1 + powf(dof * dt / gamma_fac / scales[i], 2));
        if (dof < dofmin) { dof = dofmin; }
        float chisquare = chi2_ppf_95(dof) / dof;
        float fft_theor = variance * (1 - powf(alpha, 2)) /
            (1 + powf(alpha, 2) -
             2 * alpha * cos(2 * M_PI * frequencies[i] * dt / n));
        signif[i] = fft_theor * chisquare;
    }
}

// remove trend before use it
float steps(const float *signal, int n, float dt) {
    size_t m = 85;
    float power[m], frequencies[m], signif[m];
    wavelet_power(signal, n, dt, power, frequencies, signif, m);

    int max_ind = 0;
    for (size_t i = 1; i < m; ++i) {
        if (power[i] > power[max_ind]) {
            max_ind = i;
        }
    }
    
    if (signif[max_ind] > power[max_ind]) {
        return 0;
    }

    return frequencies[max_ind] * n * dt;
}


int main() {
    FILE *f;
    f = fopen("input.dat", "r");

    size_t n;
    fscanf(f, "%lu", &n);

    float signal[n];
    for (size_t i = 0; i < n; ++i) {
        fscanf(f, "%f", signal + i);
    }

    fclose(f);

    size_t m = 85;
    float power[m], frequencies[m], signif[m];
    wavelet_power(signal, n, 1, power, frequencies, signif, m);

    f = fopen("output.dat", "w");
    for (size_t i = 0; i < m; ++i) {
        fprintf(f, "%f %f %f\n", frequencies[i], power[i], signif[i]);
    }
    fclose(f);
}
