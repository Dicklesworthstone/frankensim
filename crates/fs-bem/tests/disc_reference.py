"""Independent development oracle, never a production dependency.

Mellow & Karkkainen, JASA 118 (2005) 1311, doi:10.1121/1.2000828,
eqs. 47-49, 51-54, rigid unbaffled limit b=a. In those equations the
prefactor is sqrt(pi), not pi. Half-integer gamma ratios below are evaluated
as exact recurrences with Decimal arithmetic. The paper uses exp(+iwt);
the final impedance is conjugated for FrankenSim's exp(-iwt).
Run: python3 -I crates/fs-bem/tests/disc_reference.py
It prints truncation/precision ladders, then Rust reference values. No fit to
the BEM result or to a plot enters the oracle. Truncation M includes 0..M.
"""
from decimal import Decimal as D, localcontext
from math import factorial

PI = D("3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148086513282306647093844609550582231725359408128")


def add(a, b):
    return a[0] + b[0], a[1] + b[1]


def mul(a, b):
    return a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0]


def div(a, b):
    den = b[0]*b[0]+b[1]*b[1]
    return (a[0]*b[0]+a[1]*b[1])/den, (a[1]*b[0]-a[0]*b[1])/den


def reference(ka, truncation, precision):
    with localcontext() as ctx:
        ctx.prec = precision
        n = truncation+1
        x = D(ka)/2
        h = [D(1)]
        for j in range(1, 2*n+4):
            h.append(h[-1]*(D(j)-D("0.5")))
        f = [D(factorial(j)) for j in range(3*n+4)]
        matrix = []
        for q in range(n):
            row = []
            for m in range(n):
                b, s = D(0), D(0)
                for r in range(n):
                    b += D((-1)**(q+r))*h[m+2]*f[q+r]*x**(2*(q+r)+3)/(f[r]*f[q]**2*h[r+m+2]*h[q+r+2])
                    if q+r >= m:
                        poch = D(1)
                        for j in range(q):
                            poch *= D(r-m+j)-D("0.5")
                        s += PI*h[m+2]*D((-1)**(q+r+m))*poch*x**(2*(q+r-m))/(f[r]*f[q]**2*f[q+r-m])
                row.append((b, -s))
            matrix.append(row)
        rhs = [(D(-1) if i == 0 else D(0), D(0)) for i in range(n)]
        for j in range(n):
            pivot = max(range(j, n), key=lambda i: sum(v*v for v in matrix[i][j]))
            matrix[j], matrix[pivot] = matrix[pivot], matrix[j]
            rhs[j], rhs[pivot] = rhs[pivot], rhs[j]
            for i in range(j+1, n):
                ratio = div(matrix[i][j], matrix[j][j])
                for k in range(j+1, n):
                    product = mul(ratio, matrix[j][k])
                    matrix[i][k] = add(matrix[i][k], (-product[0], -product[1]))
                product = mul(ratio, rhs[j])
                rhs[i] = add(rhs[i], (-product[0], -product[1]))
        tau = [(D(0), D(0)) for _ in range(n)]
        for i in reversed(range(n)):
            value = rhs[i]
            for j in range(i+1, n):
                product = mul(matrix[i][j], tau[j])
                value = add(value, (-product[0], -product[1]))
            tau[i] = div(value, matrix[i][i])
        scale = -2*PI*D(ka)
        return scale*sum(t[0] for t in tau), -scale*sum(t[1] for t in tau)


if __name__ == "__main__":
    for ka in ["0.01", "0.5", "1", "2", "5"]:
        previous = None
        for truncation, precision in [(12, 80), (18, 100), (24, 120), (30, 140)]:
            value = reference(ka, truncation, precision)
            discrepancy = None if previous is None else max(abs(value[i]-previous[i]) for i in range(2))
            print(f"ka={ka} M={truncation} digits={precision} Z/(rho*c*a^2)={value} delta={discrepancy}")
            previous = value
        check = reference(ka, 30, 120)
        assert max(abs(check[i]-value[i]) for i in range(2)) < D("1e-20")
        # The final truncation pair must resolve the benchmark tolerance amply.
        assert discrepancy < D("1e-12")
        if ka == "0.01":
            assert abs(-value[1]/D(ka)-D(8)/3) < D("0.001")
        print(f"RUST ({ka}, {float(value[0]):.17e}, {float(value[1]):.17e}),")
