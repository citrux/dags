#include <stdio.h>
#include "../src/wavelet.h"

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