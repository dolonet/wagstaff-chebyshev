#!/usr/bin/env python3
"""
Расширение Theorem 4.1: BLS N+1 для Вагстаффа с СОСТАВНЫМ W_{p-2}.
Для каждого известного простого/PRP Вагстаффа:
  1. Факторизовать W_{p-2} (или частично)
  2. Проверить F = 4 × (факторизованная часть) > √W_p
  3. Верифицировать BLS условия с a=3
  4. Если всё проходит → ДОКАЗАННАЯ простота
"""

from sympy import isprime, factorint, Integer
import time
import math

def wagstaff(p):
    return (2**p + 1) // 3

def omega3_power_mod(n, W):
    """ω₃^n mod W, ω₃ = 3+2√2. Элемент (T, U) = T + U·√8 (disc=8)."""
    if n == 0:
        return (1, 0)
    d = 8  # a²-1 = 9-1 = 8
    def mul(x, y):
        return ((x[0]*y[0] + d*x[1]*y[1]) % W, (x[0]*y[1] + x[1]*y[0]) % W)
    result = (1, 0)
    base = (3 % W, 1)  # ω₃ = 3 + 1·√8
    nn = abs(n)
    while nn > 0:
        if nn % 2 == 1:
            result = mul(result, base)
        base = mul(base, base)
        nn //= 2
    return result

# Все известные простые Вагстаффа (доказанные через ECPP)
proven_wagstaff_p = [
    3, 5, 7, 11, 13, 17, 19, 23, 31, 43, 61, 79, 101, 127,
    167, 191, 199, 313, 347, 701, 1709, 2617, 3539, 5807,
    10501, 10691, 11279, 12391, 14479, 42737, 83339, 95369,
    117239, 127031, 138937, 141079
]

# Wagstaff PRP (не доказаны)
wagstaff_prp_p = [267017, 269987, 374321, 986191, 4031399, 13347311, 13372531, 15135397]

print("=" * 80)
print("РАСШИРЕННАЯ ТЕОРЕМА: BLS N+1 для Вагстаффа с составным W_{p-2}")
print("=" * 80)
print()

results = []

for p in proven_wagstaff_p:
    if p < 5:
        continue

    W = wagstaff(p)
    W_pm2_val = wagstaff(p - 2)
    Wp1 = W + 1  # = 4 * W_{p-2}

    sqrt_W = int(math.isqrt(W))

    # Пытаемся факторизовать W_{p-2}
    t0 = time.time()
    digits = len(str(W_pm2_val))

    if digits <= 70:
        # sympy может справиться
        try:
            factors = factorint(W_pm2_val, limit=10**8)  # с ограничением для скорости
            factored_part = 1
            unfactored = W_pm2_val
            for q, e in factors.items():
                factored_part *= q**e
                unfactored //= q**e

            if unfactored == 1:
                status = "ПОЛНАЯ"
            else:
                status = f"ЧАСТИЧНАЯ (остаток {len(str(unfactored))} цифр)"
        except Exception as e:
            factors = {}
            factored_part = 1
            unfactored = W_pm2_val
            status = f"ОШИБКА: {e}"
    else:
        # Слишком большое — пробуем только пробное деление
        factors = {}
        factored_part = 1
        unfactored = W_pm2_val
        # Делители W_{p-2} ≡ 1 mod 2(p-2)
        step = 2 * (p - 2)
        q = step + 1
        while q < min(10**9, unfactored) and q * q <= unfactored:
            if unfactored % q == 0:
                e = 0
                while unfactored % q == 0:
                    unfactored //= q
                    e += 1
                factors[q] = e
                factored_part *= q**e
            q += step

        if unfactored > 1 and unfactored != W_pm2_val:
            if isprime(unfactored):
                factors[unfactored] = 1
                factored_part *= unfactored
                unfactored = 1

        if unfactored == 1:
            status = "ПОЛНАЯ"
        elif factored_part > 1:
            status = f"ЧАСТИЧНАЯ ({len(str(unfactored))} цифр осталось)"
        else:
            status = "НЕ УДАЛОСЬ"

    dt = time.time() - t0

    F = 4 * factored_part
    bls_size_ok = F > sqrt_W

    result = {
        'p': p,
        'digits_W': len(str(W)),
        'digits_Wpm2': digits,
        'factors': factors,
        'factored_part': factored_part,
        'unfactored': unfactored,
        'F': F,
        'sqrt_W': sqrt_W,
        'bls_size_ok': bls_size_ok,
        'status': status,
        'time': dt,
    }

    # Если размер ОК — проверяем BLS условия
    if bls_size_ok:
        # Условие 1: ω₃^{W+1} = 1
        r1 = omega3_power_mod(Wp1, W)
        cond1 = (r1 == (1, 0))

        # Условие 2: ω₃^{(W+1)/2} = -1
        r2 = omega3_power_mod(Wp1 // 2, W)
        cond2 = (r2 == (W - 1, 0))

        # Условие 3: для каждого нечётного простого q | factored_part: ω₃^{(W+1)/q} ≠ 1
        cond3 = True
        failed_q = None
        for q in factors:
            if q == 2:
                continue  # обработано в cond2
            exp = Wp1 // q
            r3 = omega3_power_mod(exp, W)
            if r3 == (1, 0):
                cond3 = False
                failed_q = q
                break

        result['cond1'] = cond1
        result['cond2'] = cond2
        result['cond3'] = cond3
        result['proven'] = cond1 and cond2 and cond3
    else:
        result['proven'] = False

    results.append(result)

    # Вывод
    tag = "✓✓✓ ДОКАЗАНО" if result.get('proven', False) else ("F < √W" if not bls_size_ok else "условия не выполнены")
    print(f"p={p:6d} ({result['digits_W']:3d} цифр): W_{{p-2}} {status:30s} "
          f"F/√W = {float(F)/float(sqrt_W):.3f}  [{dt:.1f}s]  {tag}")

    if result.get('proven', False):
        prime_factors_str = " × ".join(f"{q}^{e}" if e > 1 else str(q) for q, e in sorted(factors.items()))
        print(f"         W_{{p-2}} = {prime_factors_str}"
              f"{(' × ' + str(result['unfactored'])) if result['unfactored'] > 1 else ''}")

print()
print("=" * 80)
print("ИТОГ")
print("=" * 80)

proven_count = sum(1 for r in results if r.get('proven', False))
total = len(results)
print(f"Доказано через BLS: {proven_count} из {total}")
print()
print("Список доказанных:")
for r in results:
    if r.get('proven', False):
        print(f"  p = {r['p']:6d}  ({r['digits_W']:3d} цифр)  W_{{p-2}} факторизовано: {r['status']}")
