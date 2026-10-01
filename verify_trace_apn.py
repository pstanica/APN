"""Reproduce the finite checks in Trace_Bivariate_APN.tex.

Uses exact binary arithmetic; numpy is used only for integer table counting.
Run: python verify_trace_apn.py
The default run covers all four parameter orbits in dimension 12 and the
cubic checks in dimensions 10 and 14. No external source files are needed.
"""
from collections import Counter
import json
from pathlib import Path
import numpy as np


class Field:
    def __init__(self, m, modulus):
        self.m, self.size, self.modulus = m, 1 << m, modulus
        self.product = [[self.multiply(a, b) for b in range(self.size)]
                        for a in range(self.size)]
        self.trace = [self.tr(a) for a in range(self.size)]
        assert set(self.trace) == {0, 1}
        # This exhaustive identity also detects a reducible modulus here.
        assert all(self.power(a, self.size-1) == 1 for a in range(1, self.size))

    def multiply(self, a, b):
        r = 0
        while b:
            if b & 1:
                r ^= a
            b >>= 1
            a <<= 1
            if a & self.size:
                a ^= self.modulus
        return r

    def power(self, a, e):
        r = 1
        while e:
            if e & 1:
                r = self.product[r][a]
            a = self.product[a][a]
            e >>= 1
        return r

    def tr(self, a):
        r = 0
        for _ in range(self.m):
            r ^= a
            a = self.product[a][a]
        return r


def tables(k, c):
    mul, pw, tr = k.product, k.power, k.trace
    inv, c2 = pw(c, k.size-2), pw(c ^ 1, 2)
    seed, switch = [], []
    for z in range(k.size*k.size):
        x, y = z & (k.size-1), z >> k.m
        t = mul[x][y]
        p = pw(x, 3) ^ mul[x][pw(y, 2)] ^ mul[c][pw(y, 3)] ^ t
        q = (pw(x, 5) ^ mul[c][mul[pw(x, 4)][y]]
             ^ mul[c2][mul[x][pw(y, 4)]] ^ mul[c][pw(y, 5)]
             ^ t ^ mul[c][pw(t, 2)])
        seed.append(p | (q << k.m))
        switch.append((p ^ mul[inv][tr[mul[c][t]]]) | (q << k.m))
    return seed, switch


def echelon(rows):
    pivots = {}
    for r in rows:
        while r:
            i = r.bit_length()-1
            if i in pivots:
                r ^= pivots[i]
            else:
                pivots[i] = r
                break
    return pivots


def ortho(vals, n):
    result = [0]
    for a in range(1, len(vals)):
        rows = [vals[a ^ (1 << i)] ^ vals[a] ^ vals[1 << i] ^ vals[0]
                for i in range(n)]
        pivots = echelon(rows)
        assert len(pivots) == n-1, ('APN rank failure', a)
        free = next(i for i in range(n) if i not in pivots)
        normal = 1 << free
        for i in sorted(pivots):
            if (normal & pivots[i]).bit_count() & 1:
                normal ^= 1 << i
        assert normal and all((normal & r).bit_count() % 2 == 0 for r in rows)
        result.append(normal)
    return result


def differential_histogram(vals):
    size = len(vals)
    vv, xx = np.array(vals, dtype=np.int64), np.arange(size)
    hist = np.zeros(size+1, dtype=np.int64)
    for a in range(1, size):
        counts = np.bincount(vv ^ vv[xx ^ a], minlength=size)
        assert int(counts.sum()) == size
        hist += np.bincount(counts, minlength=size+1)
    assert int(hist.sum()) == (size-1)*size
    assert sum(i*int(v) for i, v in enumerate(hist)) == (size-1)*size
    return [[i, int(v)] for i, v in enumerate(hist) if v]


EXPECTED_SEED = [[0,10172943],[2,5089491],[4,1269000],[6,211509],
                 [8,26829],[10,3006],[12,288],[16,18],[18,36]]
EXPECTED_SWITCH = [[0,10173780],[2,5090166],[4,1265283],[6,213372],
                   [8,27441],[10,2727],[12,234],[14,63],[16,18],[18,36]]


def comparison12():
    k = Field(6, 0x43)
    root_image = {k.power(z, 3) ^ z for z in range(k.size)}
    alphas = [a for a in range(1, k.size) if a not in root_image]
    orbits, seen = [], set()
    for a in alphas:
        if a in seen:
            continue
        orbit, z = [], a
        while z not in orbit:
            orbit.append(z)
            z = k.power(z, 2)
        seen.update(orbit)
        orbits.append(orbit)
    assert orbits == [[6,20,28,31,26,11],[14,23,25],
                      [34,52,44,35,53,45],[40,51,57,62,43,54]]
    report = []
    for orbit in orbits:
        seed, switch = tables(k, orbit[0])
        sf = differential_histogram(ortho(seed, 12))
        sg = differential_histogram(ortho(switch, 12))
        assert sf == EXPECTED_SEED and sg == EXPECTED_SWITCH
        # Directly check Frobenius covariance for both maps as well.
        fs, gs = tables(k, k.power(orbit[0], 2))
        for z in range(4096):
            z2 = k.power(z & 63, 2) | (k.power(z >> 6, 2) << 6)
            for original, squared in ((seed, fs), (switch, gs)):
                w = original[z]
                assert squared[z2] == (k.power(w & 63, 2)
                                       | (k.power(w >> 6, 2) << 6))
        item = dict(orbit=orbit, seed_image=len(set(seed)),
                    switch_image=len(set(switch)), seed_histogram=sf,
                    switch_histogram=sg)
        report.append(item)
        print('dimension 12 orbit', orbit, ': APN ranks, spectra and covariance pass',
              flush=True)
    return report


def anf_degree(vals, n):
    coefficients = vals[:]
    for i in range(n):
        for z in range(len(vals)):
            if z & (1 << i):
                coefficients[z] ^= coefficients[z ^ (1 << i)]
    return max((z.bit_count() for z, c in enumerate(coefficients) if c), default=0)


def cubic_check(m, modulus):
    k = Field(m, modulus)
    _, vals = tables(k, 1)
    mul, pw, tr = k.product, k.power, k.trace
    h, cubic = [], []
    for z in range(k.size*k.size):
        x, y = z & (k.size-1), z >> m
        p = pw(x, 3) ^ mul[x][pw(y, 2)] ^ pw(y, 3)
        q = pw(x, 5) ^ mul[pw(x, 4)][y] ^ pw(y, 5)
        h.append(tr[p ^ q])
        assert h[-1] == tr[(vals[z] & (k.size-1)) ^ (vals[z] >> m)]
    phi = [z ^ h[z] for z in range(len(vals))]
    assert len(set(phi)) == len(vals)
    assert all(phi[phi[z]] == z and h[z ^ 1] == h[z] for z in range(len(vals)))
    for z in range(len(vals)):
        x, y = z & (k.size-1), z >> m
        d1 = pw(x, 2) ^ x ^ 1 ^ pw(y, 2) ^ y ^ tr[y]
        d2 = pw(x, 4) ^ x ^ 1 ^ pw(y, 2)
        cubic.append(vals[z] ^ ((d1 | (d2 << m)) if h[z] else 0))
        assert cubic[-1] == vals[phi[z]]
    radical = [a for a in range(len(vals))
               if all(h[a ^ (1 << i)] ^ h[a] ^ h[1 << i] ^ h[0] == 0
                      for i in range(2*m))]
    assert radical == [0, 1, 1 << m, (1 << m) | 1]
    assert anf_degree(h, 2*m) == 2
    assert anf_degree(cubic, 2*m) == 3
    ortho(vals, 2*m)  # full quadratic polar-rank verification
    assert Counter(vals) == Counter(cubic)
    if m == 5:
        # Full derivative check of the nonquadratic representative.
        hist = differential_histogram(cubic)
        assert {j for j, _ in hist} == {0, 2}
    print('dimension', 2*m, ': cubic formula, radical, involution, exact degree pass',
          flush=True)
    return dict(dimension=2*m, modulus=modulus, radical=radical,
                cubic_degree=3, image_size=len(set(cubic)))


if __name__ == '__main__':
    report = dict(comparison12=comparison12(),
                  cubic=[cubic_check(5,0x25), cubic_check(7,0x83)])
    path = Path(__file__).with_name('verification_results.json')
    path.write_text(json.dumps(report, indent=2)+'\n')
    print('All finite checks passed. Results:', path.name)
