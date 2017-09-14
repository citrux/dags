#!/usr/bin/env python3

import numpy as np
import pycwt as wavelet
from matplotlib import pyplot as plt


def steps(a, dt=1e-2, minFreq=0.5, maxFreq=10):
    """
    Функция считает количество шагов на интервале,
    используя непрерывное вейвлет-преобразование.

    Изначально считается, что временной ряд однороден.
    Исполльзуется вейвлет Морле.

    Parameters
    ----------
    a : np.array
        массив модулей ускорений
    dt : float
        период дискретизации временного ряда
    minFreq : float
        минимальная детектируемая частота шагов
    maxFreq : float
        максимальная детектируемая частота шагов
    Returns
    -------
    int
        Количество шагов

    """
    n = a.size
    t = np.arange(n) * dt

    # убираем тренд
    p = np.polyfit(t, a, 1)
    aNoTrend = a - np.polyval(p, t)

    # нормируем
    std = aNoTrend.std()
    var = std ** 2
    aNorm = aNoTrend / std

    # Lag-1 autocorrelation for red noise
    try:
        alpha, _, _ = wavelet.ar1(aNorm)
    except:
        alpha = 0.9  # костыль, потом надо пофиксить

    # The next step is to define some parameters of our wavelet analysis. We
    # select the mother wavelet, in this case the Morlet wavelet with
    # :math:`\omega_0=6`.
    mother = wavelet.Morlet(6)
    s0 = 2 * dt  # Starting scale
    dj = 1 / 12  # Twelve sub-octaves per octaves
    J = 7 / dj  # Seven powers of two with dj sub-octaves
    # вейвлет преобразование
    wave, scales, freqs, _, _, _ = wavelet.cwt(aNorm, dt, dj, s0, J,
                                               mother)
    power = np.abs(wave) ** 2

    # Then, we calculate the global wavelet spectrum and determine its
    # significance level.
    glblPower = var * power.mean(axis=1)
    dof = n - scales  # Correction for padding at edges
    glblSignif, _ = wavelet.significance(var, dt, scales, 1, alpha,
                                         significance_level=0.95, dof=dof,
                                         wavelet=mother)
    index = glblPower.argmax()
    m = n * freqs[index] * dt

    plt.plot(freqs, glblPower)
    plt.plot(freqs, glblSignif)
    plt.show()

    # если спектр "хороший", то выводим посчитаное число шагов, если нет -- 0
    if minFreq < freqs[index] < maxFreq and \
                                  2 * glblSignif[index] < glblPower[index]:
        return m

    return 0

if __name__ == '__main__':
    data = np.loadtxt("input.dat")
    for i in range(0, len(data) - 200, 200):
        a = np.sqrt(np.sum(data[i:i+200,:3] ** 2, axis=1))
        steps(a)