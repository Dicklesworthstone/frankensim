use super::*;

fn limits(n: usize) -> ComplexOptions {
    ComplexOptions {
        max_panels: n,
        max_bytes: 64 * n * n,
        max_work: 512 * (n * n) as u64,
    }
}

fn box_surface() -> SpherePanels {
    let points = [
        [-0.2, -0.15, -0.001],
        [0.2, -0.15, -0.001],
        [0.2, 0.15, -0.001],
        [-0.2, 0.15, -0.001],
        [-0.2, -0.15, 0.001],
        [0.2, -0.15, 0.001],
        [0.2, 0.15, 0.001],
        [-0.2, 0.15, 0.001],
    ];
    let faces = [
        [0, 2, 1],
        [0, 3, 2],
        [4, 5, 6],
        [4, 6, 7],
        [0, 1, 5],
        [0, 5, 4],
        [1, 2, 6],
        [1, 6, 5],
        [2, 3, 7],
        [2, 7, 6],
        [3, 0, 4],
        [3, 4, 7],
    ];
    SpherePanels::from_triangles(faces.iter().map(|f| f.map(|i| points[i])).collect()).unwrap()
}

fn sphere(level: usize, radius: f64) -> SpherePanels {
    let p = [
        [1., 0., 0.],
        [-1., 0., 0.],
        [0., 1., 0.],
        [0., -1., 0.],
        [0., 0., 1.],
        [0., 0., -1.],
    ];
    let faces = [
        [0, 2, 4],
        [2, 1, 4],
        [1, 3, 4],
        [3, 0, 4],
        [2, 0, 5],
        [1, 2, 5],
        [3, 1, 5],
        [0, 3, 5],
    ];
    let mut triangles: Vec<_> = faces.iter().map(|f| f.map(|i| p[i])).collect();
    let mid = |a: [f64; 3], b: [f64; 3]| {
        let s: [f64; 3] = core::array::from_fn(|i| a[i] + b[i]);
        let l = norm(s);
        s.map(|v| v / l)
    };
    for _ in 0..level {
        triangles = triangles
            .iter()
            .flat_map(|&[a, b, c]| {
                let ab = mid(a, b);
                let bc = mid(b, c);
                let ca = mid(c, a);
                [[a, ab, ca], [ab, b, bc], [ca, bc, c], [ab, bc, ca]]
            })
            .collect();
    }
    SpherePanels::from_triangles(
        triangles
            .into_iter()
            .map(|t| t.map(|v| v.map(|x| x * radius)))
            .collect(),
    )
    .unwrap()
}

fn relative_error(actual: &[C64], expected: &[C64]) -> f64 {
    let scale = expected
        .iter()
        .map(|v| v.abs())
        .fold(0.0_f64, f64::max)
        .max(1e-30);
    actual
        .iter()
        .zip(expected)
        .map(|(a, b)| (*a - *b).abs() / scale)
        .fold(0.0, f64::max)
}

#[test]
fn g0_positive_real_axis_is_original_operator_bitwise() {
    let surface = box_surface();
    let options = limits(12);
    for formulation in [
        Formulation::PlainCbie,
        Formulation::BurtonMiller,
        Formulation::BurtonMillerWrongAlphaSign,
    ] {
        for k in [0.03, 0.7, 2.0] {
            let old = assemble_dense(&surface, k, alpha_for(formulation, k));
            let new =
                assemble_with_cancel(&surface, C64::from_re(k), formulation, &options, || false)
                    .unwrap();
            for (a, b) in old.0.iter().chain(&old.1).zip(new.a.iter().chain(&new.b)) {
                assert_eq!(
                    (a.re.to_bits(), a.im.to_bits()),
                    (b.re.to_bits(), b.im.to_bits())
                );
            }
        }
    }
}

#[test]
fn g1_complex_derivative_matches_two_independent_frequency_directions() {
    let surface = box_surface();
    let options = limits(12);
    for formulation in [Formulation::PlainCbie, Formulation::BurtonMiller] {
        for k in [
            C64::new(1.2, -0.3),
            C64::new(0.04, 0.013),
            C64::new(1.2, 0.3),
        ] {
            let op = assemble_with_cancel(&surface, k, formulation, &options, || false).unwrap();
            for direction in [C64::ONE, C64::new(0., 1.)] {
                // An eighth-order centered oracle resolves weak frequency
                // slopes beside O(1) static jumps on close opposite faces.
                // A tiny two-point step loses those slopes to cancellation.
                // This metre-scale fixture needs a wavenumber reference of
                // 1/m to resolve its weak low-frequency CBIE slopes. The BM
                // coupling has a pole at zero, so its stencil stays local to k.
                let step = if formulation == Formulation::PlainCbie {
                    2e-3 * k.abs().max(1.0)
                } else {
                    2e-3 * k.abs()
                };
                let h = direction.scale(step);
                let mut da = vec![C64::ZERO; op.a.len()];
                let mut db = vec![C64::ZERO; op.b.len()];
                for (j, coefficient) in [0.8, -0.2, 4.0 / 105.0, -1.0 / 280.0]
                    .into_iter()
                    .enumerate()
                {
                    let delta = h.scale((j + 1) as f64);
                    let plus =
                        assemble_with_cancel(&surface, k + delta, formulation, &options, || false)
                            .unwrap();
                    let minus =
                        assemble_with_cancel(&surface, k - delta, formulation, &options, || false)
                            .unwrap();
                    for index in 0..da.len() {
                        da[index] = da[index] + (plus.a[index] - minus.a[index]).scale(coefficient);
                        db[index] = db[index] + (plus.b[index] - minus.b[index]).scale(coefficient);
                    }
                }
                for value in da.iter_mut().chain(&mut db) {
                    *value = *value * h.recip();
                }
                assert!(
                    relative_error(&da, &op.da_dk) < 3e-8,
                    "A {:?} {:?} {}",
                    formulation,
                    k,
                    relative_error(&da, &op.da_dk)
                );
                assert!(
                    relative_error(&db, &op.db_dk) < 3e-8,
                    "B {:?} {:?} {}",
                    formulation,
                    k,
                    relative_error(&db, &op.db_dk)
                );
            }
        }
    }
}

#[test]
fn g1_complex_pulsating_sphere_converges_to_outgoing_analytic_impedance() {
    let radius = 0.2;
    let k = C64::new(2.0, -0.35);
    let medium = Medium {
        density: 1.2,
        sound_speed: 343.0,
    };
    let ka = k.scale(radius);
    let exact = (ka * (ka + C64::new(0., 1.)).recip()).scale(medium.density * medium.sound_speed);
    let mut errors = Vec::new();
    for level in 0..3 {
        let surface = sphere(level, radius);
        let n = surface.areas().len();
        let op = assemble_with_cancel(&surface, k, Formulation::BurtonMiller, &limits(n), || false)
            .unwrap();
        let mut pressure = vec![C64::ZERO; n];
        let q = (C64::new(0., 1.) * k).scale(medium.density * medium.sound_speed);
        for i in 0..n {
            for j in 0..n {
                pressure[i] = pressure[i] + op.b[i * n + j] * q;
            }
        }
        lu_complex(&op.a, n).unwrap().solve(&mut pressure);
        let area: f64 = surface.areas().iter().sum();
        let mean = pressure
            .iter()
            .zip(surface.areas())
            .fold(C64::ZERO, |s, (p, a)| s + p.scale(*a / area));
        let error = (mean - exact).abs() / exact.abs();
        errors.push(error);
        println!("complex sphere panels={n} relative_impedance_error={error:.8e}");
    }
    assert!(
        errors[1] < errors[0] * 0.7 && errors[2] < errors[1] * 0.7,
        "{errors:?}"
    );
    assert!(errors[2] < 0.06, "{errors:?}");
}

#[test]
fn g0_complex_operator_caps_and_cancellation_are_transactional() {
    let surface = box_surface();
    let mut options = limits(12);
    let k = C64::new(1., -0.1);
    assemble_with_cancel(&surface, k, Formulation::PlainCbie, &options, || false).unwrap();
    options.max_bytes -= 1;
    assert!(matches!(
        assemble_with_cancel(&surface, k, Formulation::PlainCbie, &options, || false),
        Err(ComplexError::Budget { .. })
    ));
    options = limits(12);
    options.max_work -= 1;
    assert!(matches!(
        assemble_with_cancel(&surface, k, Formulation::PlainCbie, &options, || false),
        Err(ComplexError::Budget { .. })
    ));
    options = limits(12);
    options.max_panels = 11;
    assert!(matches!(
        assemble_with_cancel(&surface, k, Formulation::PlainCbie, &options, || false),
        Err(ComplexError::PanelBudget { .. })
    ));
    options = limits(12);
    let mut calls = 0;
    assert_eq!(
        assemble_with_cancel(&surface, k, Formulation::PlainCbie, &options, || {
            calls += 1;
            calls == 4
        })
        .unwrap_err(),
        ComplexError::Cancelled
    );
    assert_eq!(calls, 4);
    for invalid in [
        C64::ZERO,
        C64::new(f64::NAN, 0.0),
        C64::new(0.0, f64::INFINITY),
    ] {
        assert!(
            assemble_with_cancel(
                &surface,
                invalid,
                Formulation::BurtonMiller,
                &options,
                || false
            )
            .is_err()
        );
    }
    let first =
        assemble_with_cancel(&surface, k, Formulation::PlainCbie, &options, || false).unwrap();
    let second =
        assemble_with_cancel(&surface, k, Formulation::PlainCbie, &options, || false).unwrap();
    for (a, b) in first
        .a
        .iter()
        .chain(&first.b)
        .chain(&first.da_dk)
        .chain(&first.db_dk)
        .zip(
            second
                .a
                .iter()
                .chain(&second.b)
                .chain(&second.da_dk)
                .chain(&second.db_dk),
        )
    {
        assert_eq!(
            (a.re.to_bits(), a.im.to_bits()),
            (b.re.to_bits(), b.im.to_bits())
        );
    }
}
