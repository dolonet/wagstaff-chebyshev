"""cp354: FactorDB sweep for V_{2d}, d in {59,81,83,99,107,121} — the
remaining k<=500 diagonal d's (same-d-exclusion closures need complete
V_{2d} factorizations). Prints the numbers; query via curl below."""
def V(n):
    a, b = 2, 2
    for _ in range(n - 1):
        a, b = b, 2 * b + a
    return b
for d in [59, 81, 83, 99, 107, 121]:
    v = V(2 * d)
    print(f"d={d} V_{2*d} ({len(str(v))} digits): {v}")
