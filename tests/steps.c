#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include "../src/wavelet.h"

int main(int argc, char const *argv[])
{
    if (argc < 2) {
        exit(1);
    }

    FILE *f;
    f = fopen(argv[1], "r");
    int16_t x, y, z, a, b, steps = 0;
    uint16_t t;
    size_t i = 0;
    while(fscanf(f, "%"SCNd16" %"SCNd16" %"SCNd16" %"SCNu16" %"SCNd16" %"SCNd16, &x, &y, &z, &t, &a, &b) != EOF) {
        if (i < 200) {
            waveletProcessNewData(x, y, z, t);
            ++i;
        } else {
            steps += waveletGetStepsCount();
            i = 0;
        }
    }
    fclose(f);
    printf("%d\n", steps);
    return 0;
}