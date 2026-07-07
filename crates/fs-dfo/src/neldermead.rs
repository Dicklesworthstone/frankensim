//! Nelder–Mead simplex (the cheap deterministic polish baseline).
//! Standard reflection/expansion/contraction/shrink coefficients,
//! deterministic vertex ordering (`total_cmp`, lowest index on ties).

/// Minimize `f` from `x0` with initial simplex scale `h`. Returns
/// (x_best, f_best, evals). Deterministic (no randomness at all).
#[must_use]
pub fn nelder_mead<F: FnMut(&[f64]) -> f64>(
    f: &mut F,
    x0: &[f64],
    h: f64,
    max_evals: usize,
    f_target: f64,
) -> (Vec<f64>, f64, usize) {
    let n = x0.len();
    let (alpha, gamma, rho, sigma) = (1.0, 2.0, 0.5, 0.5);
    // Initial simplex: x0 plus axis steps.
    let mut xs: Vec<Vec<f64>> = Vec::with_capacity(n + 1);
    xs.push(x0.to_vec());
    for i in 0..n {
        let mut v = x0.to_vec();
        v[i] += h;
        xs.push(v);
    }
    let mut fs: Vec<f64> = xs.iter_mut().map(|x| f(x)).collect();
    let mut evals = n + 1;
    while evals < max_evals {
        // Order vertices (deterministic).
        let mut idx: Vec<usize> = (0..=n).collect();
        idx.sort_by(|&a, &b| fs[a].total_cmp(&fs[b]).then(a.cmp(&b)));
        let (best, worst, second_worst) = (idx[0], idx[n], idx[n - 1]);
        if fs[best] <= f_target {
            break;
        }
        // Centroid of all but worst.
        let mut c = vec![0.0f64; n];
        for &k in idx.iter().take(n) {
            for i in 0..n {
                c[i] += xs[k][i] / n as f64;
            }
        }
        let blend = |a: &[f64], b: &[f64], t: f64| -> Vec<f64> {
            a.iter().zip(b).map(|(p, q)| t.mul_add(q - p, *p)).collect()
        };
        // Reflect.
        let xr = blend(&c, &xs[worst], -alpha);
        let fr = f(&xr);
        evals += 1;
        if fr < fs[best] {
            // Expand.
            let xe = blend(&c, &xs[worst], -gamma);
            let fe = f(&xe);
            evals += 1;
            if fe < fr {
                xs[worst] = xe;
                fs[worst] = fe;
            } else {
                xs[worst] = xr;
                fs[worst] = fr;
            }
        } else if fr < fs[second_worst] {
            xs[worst] = xr;
            fs[worst] = fr;
        } else {
            // Contract (outside/inside).
            let (xc, fc) = if fr < fs[worst] {
                let x = blend(&c, &xs[worst], -rho);
                let v = f(&x);
                (x, v)
            } else {
                let x = blend(&c, &xs[worst], rho);
                let v = f(&x);
                (x, v)
            };
            evals += 1;
            if fc < fs[worst].min(fr) {
                xs[worst] = xc;
                fs[worst] = fc;
            } else {
                // Shrink toward best.
                let xb = xs[best].clone();
                for k in 0..=n {
                    if k == best {
                        continue;
                    }
                    xs[k] = blend(&xb, &xs[k], sigma);
                    fs[k] = f(&xs[k]);
                    evals += 1;
                }
            }
        }
    }
    let mut bi = 0usize;
    for k in 1..=n {
        if fs[k] < fs[bi] {
            bi = k;
        }
    }
    (xs[bi].clone(), fs[bi], evals)
}
