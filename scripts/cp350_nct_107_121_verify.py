#!/usr/bin/env python3
"""Verify the d=107 and d=121 NCT closures used in Paper B v3.3.

Checks:
  * given factors multiply to the Pell primitive part Psi_{4d};
  * all listed factors are prime (SymPy plus PARI/GP isprime);
  * split factors r == 1 mod 8 have composite ord_r(2)/2;
  * inert factors are irrelevant to the split compatible-triple search.
"""

from __future__ import annotations

from fractions import Fraction
import re
import subprocess

import sympy as sp


DATA = {
    107: [
        161_522_909_729,
        1_262_880_768_427,
        17_449_564_219_651,
        2_906_272_939_091_099,
        132_168_496_320_829_356_461_361_353_683,
    ],
    121: [
        1_451,
        568_699,
        1_197_595_081,
        1_614_230_507,
        135_855_589_627,
        7_499_151_961_531,
        249_507_118_668_323,
        4_005_463_110_187_137_017,
    ],
}


def pell_u(n: int) -> int:
    if n == 0:
        return 0
    a, b = 0, 1
    for _ in range(1, n):
        a, b = b, 2 * b + a
    return b


def psi(n: int) -> int:
    """Möbius primitive part Psi_n = prod_{e|n} U_e^{mu(n/e)}."""
    value = Fraction(1, 1)
    for e in sp.divisors(n):
        mu = sp.mobius(n // e)
        if mu == 1:
            value *= pell_u(e)
        elif mu == -1:
            value /= pell_u(e)
    assert value.denominator == 1
    return int(value)


def pari_isprime(n: int) -> bool:
    proc = subprocess.run(
        ["gp", "-q"],
        input=f"isprime({n})\n\\q\n",
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True,
    )
    cleaned = re.sub(r"\x1b\[[0-9;]*m", "", proc.stdout).strip()
    return cleaned == "1"


def main() -> None:
    for d, factors in DATA.items():
        n = 4 * d
        primitive_part = psi(n)
        product = sp.prod(factors)
        print(f"d={d}, n={n}")
        print(f"  digits(Psi_{{4d}})={len(str(primitive_part))}")
        print(f"  product matches Psi_{{4d}}: {product == primitive_part}")
        assert product == primitive_part

        for r in factors:
            sympy_prime = sp.isprime(r)
            pari_prime = pari_isprime(r)
            assert sympy_prime and pari_prime
            residue = r % 8
            print(f"  r={r} mod8={residue} prime=yes")
            if residue == 1:
                order = sp.n_order(2, r)
                half = order // 2
                half_factors = sp.factorint(half)
                print(f"    ord_r(2)={order}")
                print(f"    ord_r(2)/2={half} factors={half_factors}")
                assert order % 2 == 0
                assert not sp.isprime(half)

    print("OK: d=107 and d=121 split primitive factors have composite half-orders.")


if __name__ == "__main__":
    main()
