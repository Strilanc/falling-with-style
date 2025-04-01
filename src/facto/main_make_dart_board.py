import fractions
import itertools
import math
import random

import sympy
from matplotlib import pyplot as plt

lines = []

lines.append("""<svg viewBox="-100 0 2580 2680"  version="1.1" xmlns="http://www.w3.org/2000/svg">""")
lines.append(f"""<text x="0" y="{2580/2}" font-size="100" text-anchor="end" dominant-baseline="middle">g</text>""")
lines.append(f"""<text x="{2580/2}" y="{2580}" font-size="100" text-anchor="middle" dominant-baseline="hanging">m</text>""")

n = 7*13
m_range = 4**n.bit_length()
mps = []
for m in range(m_range):
    r = fractions.Fraction(m, m_range)
    p = r.limit_denominator(n).denominator
    mps.append((m, p))

h = 2580 / (n - 4)
w = 2580 / m_range
y = 0
hits = 0
attempts = 0
hits2 = 0
for g in range(2, n - 2):
    if math.gcd(n, g) != 1:
        hits2 += m_range
        attempts += m_range
        lines.append(f"""<rect x="{0}" y="{y}" width="{2580}" height="{h}" stroke="none" fill="green"/>""")
    else:
        x = 0
        for m, p in mps:
            f = math.gcd(n, pow(g, p >> 1, n) + 1)
            if 1 < f < n and n % f == 0:
                hits += 1
                lines.append(f"""<rect x="{x}" y="{y}" width="{w}" height="{h}" stroke="none" fill="blue"/>""")
            attempts += 1
            x += w
    y += h
lines.append(f"""<rect x="0" y="0" width="{2580}" height="{2580}" stroke="black" fill="none"/>""")
lines.append("""</svg>""")

with open('dartboard.svg', 'w') as f:
    print('\n'.join(lines), file=f)

