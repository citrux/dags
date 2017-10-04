#include <stdbool.h>
#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <math.h>

#define MAX_BUFFER_SIZE 256

int16_t SIN[MAX_BUFFER_SIZE * 3 / 4] = {0, 804, 1608, 2410, 3212, 4011, 4808, 5602, 6393, 7179, 7962, 8739, 9512, 10278, 11039, 11793, 12539, 13279, 14010, 14732, 15446, 16151, 16846, 17530, 18204, 18868, 19519, 20159, 20787, 21403, 22005, 22594, 23170, 23731, 24279, 24811, 25329, 25832, 26319, 26790, 27245, 27683, 28105, 28510, 28898, 29268, 29621, 29956, 30273, 30571, 30852, 31113, 31356, 31580, 31785, 31971, 32137, 32285, 32412, 32521, 32609, 32678, 32728, 32757, 32767, 32757, 32728, 32678, 32609, 32521, 32412, 32285, 32137, 31971, 31785, 31580, 31356, 31113, 30852, 30571, 30273, 29956, 29621, 29268, 28898, 28510, 28105, 27683, 27245, 26790, 26319, 25832, 25329, 24811, 24279, 23731, 23170, 22594, 22005, 21403, 20787, 20159, 19519, 18868, 18204, 17530, 16846, 16151, 15446, 14732, 14010, 13279, 12539, 11793, 11039, 10278, 9512, 8739, 7962, 7179, 6393, 5602, 4808, 4011, 3212, 2410, 1608, 804, 0, -804, -1608, -2410, -3212, -4011, -4808, -5602, -6393, -7179, -7962, -8739, -9512, -10278, -11039, -11793, -12539, -13279, -14010, -14732, -15446, -16151, -16846, -17530, -18204, -18868, -19519, -20159, -20787, -21403, -22005, -22594, -23170, -23731, -24279, -24811, -25329, -25832, -26319, -26790, -27245, -27683, -28105, -28510, -28898, -29268, -29621, -29956, -30273, -30571, -30852, -31113, -31356, -31580, -31785, -31971, -32137, -32285, -32412, -32521, -32609, -32678, -32728, -32757};

typedef struct {int32_t re; int32_t im;} complex;

complex cadd(complex u, complex v) {
    complex w = {u.re + v.re,  u.im + v.im};
    return w;
}

complex csub(complex u, complex v) {
    complex w = {u.re - v.re,  u.im - v.im};
    return w;
}

complex cmul(complex u, complex v) {
    complex w = {u.re * v.re - u.im * v.im,
                 u.re * v.im + u.im * v.re};
    return w;
}

void fft_bit_rev(complex const *in,
                 complex *out,
                 bool inverse) {

    uint16_t bit_reverse(uint16_t a, uint16_t b) {
        uint16_t result = 0;
        while (b) {
            result <<= 1;
            result |= a & 1;
            a >>= 1;
            b >>= 1;
        }
        return result;
    }

    for (uint16_t i = 0; i < MAX_BUFFER_SIZE; ++i) {
        out[i] = in[bit_reverse(i, MAX_BUFFER_SIZE - 1)];
    }
    for (uint16_t m = 1; m < MAX_BUFFER_SIZE; m *= 2) {
        for (uint16_t k = 0; k < MAX_BUFFER_SIZE; k += 2 * m) {
            for (uint16_t j = 0; j < m; ++j) {
                complex e = {SIN[j * MAX_BUFFER_SIZE / 2 / m + MAX_BUFFER_SIZE / 4],
                             -SIN[j * MAX_BUFFER_SIZE / 2 / m]};
                if (inverse) { e.im = -e.im; }
                complex u = out[k + j];
                complex v = cmul(e, out[k + j + m]);
                v.re >>= 15;
                v.im >>= 15;
                out[k + j] = cadd(u, v);
                out[k + j + m] = csub(u, v);
            }
        }
    }
}


int fft(complex const *in, complex *out) {
    fft_bit_rev(in, out, false);
    return 0;
}

int ifft(complex const *in, complex *out) {
    fft_bit_rev(in, out, true);
    for (uint16_t i = 0; i < MAX_BUFFER_SIZE; ++i) {
        out[i].re /= MAX_BUFFER_SIZE;
        out[i].im /= MAX_BUFFER_SIZE;
    }
    return 0;
}

int main(int argc, char const *argv[])
{
    if (argc < 1) { return 1; }

    complex signal[MAX_BUFFER_SIZE];
    complex spectrum[MAX_BUFFER_SIZE];
    complex signal_rev[MAX_BUFFER_SIZE];
    
    FILE *fd;
    fd = fopen(argv[1], "r");
    int16_t x, y, z, a, b;
    uint16_t t, i = 0;
    while(fscanf(fd, "%"SCNd16" %"SCNd16" %"SCNd16" %"SCNu16" %"SCNd16" %"SCNd16, &x, &y, &z, &t, &a, &b) != EOF) {
        float xf = x;
        float yf = y;
        float zf = z;
        signal[i].im = 0;
        signal[i++].re = sqrtf(xf * xf + yf * yf + zf * zf);
    }
    fclose(fd);
    fft(signal, spectrum);
    fd = fopen("spectrum.dat", "w");
    for (int i = 0; i < MAX_BUFFER_SIZE; ++i)
    {
        fprintf(fd, "%d %d\n", spectrum[i].re, spectrum[i].im);
    }
    fclose(fd);

    ifft(spectrum, signal_rev);
    fd = fopen("signal_rev.dat", "w");
    for (int i = 0; i < MAX_BUFFER_SIZE; ++i)
    {
        fprintf(fd, "%d %d\n", signal_rev[i].re, signal_rev[i].im);
    }
    fclose(fd);
    return 0;
}