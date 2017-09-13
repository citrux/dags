#include <stdint.h>

void waveletProcessNewData(uint16_t x, uint16_t y, uint16_t z, uint16_t time);
uint16_t waveletGetStepsCount(void);

// for testing purposes
void wavelet_power(const float *signal, size_t n, float dt, float *power,
                   float *frequencies, float *signif, size_t m);