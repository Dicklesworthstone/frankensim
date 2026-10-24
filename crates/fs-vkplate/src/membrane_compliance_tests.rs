use super::*;

fn plate(h: f64, a: f64, nu: f64) -> Body {
    Body::plate(
        "membrane compliance",
        &PlateSpec {
            radius: 0.65,
            hole: 0.0,
            grid: h,
            edge: Edge::Free,
            thickness: &|r| 0.01 / (1.0 + a * r * r),
            rest_height: &|_| 0.0,
            young: 100.0,
            poisson: nu,
            density: 1000.0,
            tension: 0.0,
            bending: true,
            nonlinear: Nonlinear::VonKarman,
            t60: [(100.0, f64::INFINITY), (1000.0, f64::INFINITY)],
            support: None,
            muffle: None,
        },
    )
    .unwrap()
}

#[test]
fn uniform_compliance_retains_the_laplacian_squared_operator() {
    for nu in [-0.2, 0.0, 0.33, 0.49] {
        let b = plate(0.08, 0.0, nu);
        let vk = b.vk.as_ref().unwrap();
        let nc = vk.nb.len();
        let f: Vec<_> = (0..nc).map(|i| (i as f64 * 1.7).sin()).collect();
        let mut full = vec![0.0; b.n];
        for (stencil, &f) in vk.nb.iter().zip(&f) {
            for (node, weight) in stencil[..5].iter().zip([-4.0, 1.0, 1.0, 1.0, 1.0]) {
                if let Some(node) = node {
                    full[*node] += weight * f / (b.h * b.h);
                }
            }
        }
        let expected: f64 = full.iter().map(|x| x * x).sum();
        let mut af = vec![0.0; nc];
        vk.matrix.mul_add(&f, 1.0, &mut af);
        let actual: f64 = f.iter().zip(af).map(|(f, a)| f * a).sum();
        assert!(
            (actual - expected).abs() < 2e-13 * expected,
            "nu={nu}: {actual} != {expected}"
        );
    }
}

#[test]
fn local_compliance_matches_manufactured_plane_stress_equilibrium() {
    // E*t = 1/(1+a*r^2). The plane-stress complementary energy gives
    // A(F) = dxx[c(Fxx-nu Fyy)] + dyy[c(Fyy-nu Fxx)]
    //        + 2(1+nu) dxy[c Fxy], c=1+a*r^2.
    // Analytic independent polynomial forcings:
    // F=x^4+y^4 -> 48+24*a*(7-nu)*r^2;
    // F=x^2*y^2 -> 8+(28-4*nu)*a*r^2. The normal terms contribute
    // -8*nu+4*a*(1-7*nu)*r^2; the shear term contributes
    // 8*(1+nu)*(1+3*a*r^2). Retain the Poisson term for variable c.
    for nu in [-0.2, 0.33, 0.49] {
        for mixed in [false, true] {
            let mut previous = f64::INFINITY;
            for h in [0.08, 0.04, 0.02] {
                let a = 4.0;
                let b = plate(h, a, nu);
                let vk = b.vk.as_ref().unwrap();
                let f: Vec<_> = vk
                    .nb
                    .iter()
                    .map(|s| {
                        let [x, y] = b.xy[s[0].unwrap()];
                        if mixed {
                            x * x * y * y
                        } else {
                            x.powi(4) + y.powi(4)
                        }
                    })
                    .collect();
                let mut af = vec![0.0; f.len()];
                vk.matrix.mul_add(&f, 1.0, &mut af);
                let mut error = 0.0_f64;
                for (s, computed) in vk.nb.iter().zip(af) {
                    let [x, y] = b.xy[s[0].unwrap()];
                    if x.abs() > 0.12 || y.abs() > 0.12 {
                        continue;
                    }
                    let expected = if mixed {
                        8.0 + (28.0 - 4.0 * nu) * a * (x * x + y * y)
                    } else {
                        48.0 + 24.0 * a * (7.0 - nu) * (x * x + y * y)
                    };
                    error = error.max((computed / vk.eh - expected).abs());
                }
                eprintln!("compliance nu={nu} mixed={mixed} h={h} error={error:e}");
                assert!(
                    error < 0.27 * previous,
                    "local compliance must converge at second order"
                );
                previous = error;
            }
            assert!(
                previous < 0.15,
                "manufactured force did not converge: {previous}"
            );
        }
    }
}

#[test]
fn variable_compliance_force_is_the_potential_gradient() {
    let mut b = plate(0.08, 4.0, 0.33);
    for (w, [x, y]) in b.w.iter_mut().zip(&b.xy) {
        *w = 0.01 * (3.0 * x + 2.0 * y).sin();
    }
    let v = b.nonlinear_gradient();
    assert!(v > 0.0);
    for i in (0..b.n).step_by(11) {
        let original = b.w[i];
        b.nonlinear_gradient();
        let expected = b.g[i];
        // V is quartic in w. Richardson cancels the central difference's
        // exact quadratic truncation term, including for tiny edge forces.
        // The remaining bound accounts explicitly for energy subtraction.
        let mut derivatives = [0.0; 2];
        let mut roundoff = 0.0;
        for (j, delta) in [1e-5, 5e-6].into_iter().enumerate() {
            b.w[i] = original + delta;
            let upper = b.nonlinear_gradient();
            b.w[i] = original - delta;
            let lower = b.nonlinear_gradient();
            derivatives[j] = (upper - lower) / (2.0 * delta);
            roundoff += 128.0 * f64::EPSILON * (upper.abs() + lower.abs()) / delta;
        }
        b.w[i] = original;
        let observed = (4.0 * derivatives[1] - derivatives[0]) / 3.0;
        assert!(
            (observed - expected).abs() <= 2e-6 * expected.abs() + roundoff,
            "force is not the derivative of stored membrane energy: {observed} != {expected}"
        );
    }
}
