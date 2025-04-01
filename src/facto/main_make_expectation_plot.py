import fractions
import math

import numpy as np
import sympy
from matplotlib import pyplot as plt
import multiprocessing


def _count_samples_for_target_number(n: int):
    if sympy.isprime(n) or n % 2 == 0:
        return n, 1, 1, 0, 0

    m_range = 4 ** n.bit_length()
    mps = []
    for m in range(m_range):
        r = fractions.Fraction(m, m_range)
        p = r.limit_denominator(n).denominator
        mps.append((m, p))

    attempts = 0
    hits_noiseless = 0
    hits_shared = 0
    hits_random = 0
    for g in range(2, n - 2):
        if math.gcd(n, g) != 1:
            hits_shared += m_range
            attempts += m_range
        else:
            # Brute force the period.
            period = 1
            acc = g
            while acc != 1:
                period += 1
                acc *= g
                acc %= n

            # Make a signal with that period, and sample its frequency spectrum.
            signal = np.zeros(m_range, dtype=np.complex128)
            signal[::period] = 1 / math.sqrt(math.ceil(m_range / period))
            signal = np.fft.fft(signal, norm='ortho')
            signal = (signal*np.conj(signal)).real

            attempts += m_range
            for m, p in mps:
                f = math.gcd(n, pow(g, p >> 1, n) + 1)
                if 1 < f < n and n % f == 0:
                    hits_noiseless += m_range * signal[m]
                    hits_random += 1
    return n, attempts, hits_shared, hits_noiseless, hits_random


def stats_to_expected_samples(attempts: int, hits_sample: int, hits_shared: int) -> int:
    e = 0
    p_early = hits_shared / attempts
    p_sample = hits_sample / attempts
    for k in range(50):
        q = (1 - p_early - p_sample) ** k
        e += q * p_early * k
        e += q * p_sample * (k + 1)
    return e


def main():
    vals = []
    for v in range(2, 256):
        vals.append((v, str(v)))

    fig: plt.Figure
    ax: plt.Axes
    fig, ax = plt.subplots()
    pool = multiprocessing.Pool()
    tots = pool.map(_count_samples_for_target_number, [n for n, _ in vals])
    stats = []
    for n, attempts, hits_shared, hits_noiseless, hits_random in tots:
        e1 = stats_to_expected_samples(attempts, hits_sample=hits_noiseless, hits_shared=hits_shared)
        e2 = stats_to_expected_samples(attempts, hits_sample=hits_random, hits_shared=hits_shared)
        if e1 > 0:
            stats.append((n, e1, e2))

    stats = sorted(stats, key=lambda e: e[2])
    ax.bar(range(len(stats)), [e[2] for e in stats], color='red', label='random samples', width=0.9)
    ax.bar(range(len(stats)), [e[1] for e in stats], color='blue', label='noiseless quantum samples', width=0.8)
    ax.set_xticks(range(len(stats)), labels=[str(e[0]) for e in stats], rotation=90)
    ax.set_ylabel("Expected Samples to Factor Number", fontdict={'size': 20})
    ax.set_yticks(range(8), [str(e) for e in range(8)], fontdict={'size': 16})
    ax.legend(fontsize=12)
    fig.set_dpi(200)
    fig.set_size_inches(12, 8)
    fig.tight_layout()
    fig.savefig('expectation.png')
    plt.show()


if __name__ == '__main__':
    main()
