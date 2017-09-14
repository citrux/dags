#include <stdio.h>
#include <math.h>
#include "../src/wavelet.h"

int main(int argc, char const *argv[])
{
    double signal[200];
    FILE *in = fopen("input.dat", "r");
    int i = 0, j = 0;
    int x, y, z, t, a, b;
    while(fscanf(in, "%d %d %d %d %d %d", &x, &y, &z, &t, &a, &b) != EOF) {
        if (i < 200) {
            signal[i] = sqrt(x * x + y * y + z * z);
            ++i;
        } else {
            double dt = 1e-2;
            size_t m = 85;
            double power[m], frequencies[m], signif[m];
            wavelet_power(signal, 200, dt, power, frequencies, signif, m);
            
            char outname[64];
            sprintf(outname, "output_%02d.dat", j);
            FILE *out = fopen(outname, "w");
            for (int k = 0; k < m; ++k)
            {
                fprintf(out, "%f %f %f\n", frequencies[k], power[k], signif[k]);
            }
            fclose(out);
            i = 0;
            ++j;
        }
    }
    fclose(in);
    return 0;
}