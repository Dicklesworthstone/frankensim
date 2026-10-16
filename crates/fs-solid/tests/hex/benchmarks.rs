use super::{support::*, tet_support};
use fs_math::det;
const PI: f64 = std::f64::consts::PI;

fn overlap(a: &TetElasticAssembly, phi: &[f64], reference: &[f64]) -> f64 {
    let mut mr = vec![0.0; reference.len()];
    a.mass.spmv(reference, &mut mr);
    let norm = reference.iter().zip(&mr).map(|(x, y)| x * y).sum::<f64>();
    phi.iter().zip(mr).map(|(x, y)| x * y).sum::<f64>().abs() / det::sqrt(norm)
}
pub fn roller() {
    let started = std::time::Instant::now();
    with_cx(|cx| {
        for family in FAMILIES {
            let mut errors = [Vec::new(), Vec::new()];
            let mut dofs = Vec::new();
            let hs = vec![1.0 / 8.0, 1.0 / 16.0, 1.0 / 32.0];
            for nx in [8, 16, 32] {
                let size = [1.0, 0.7, 0.4];
                let mesh = mesh([nx, 1, 1], size, family, cx);
                let fixed = mesh
                    .nodes_m
                    .iter()
                    .enumerate()
                    .flat_map(|(i, p)| {
                        (0..3).filter_map(move |k| {
                            if p[k].abs() < 1e-13 || (p[k] - size[k]).abs() < 1e-13 {
                                Some(3 * i + k)
                            } else {
                                None
                            }
                        })
                    })
                    .collect::<Vec<_>>();
                let a = assembly(&mesh, &material(1.0, 0.3), &fixed, cx);
                dofs.push(a.free_dofs.len());
                for mode in 1..=2 {
                    let mu = 1.0 / 2.6;
                    let lambda = 0.3 / (1.3 * 0.4);
                    let exact = (mode as f64 * PI).powi(2) * (lambda + 2.0 * mu);
                    let reference = a
                        .free_dofs
                        .iter()
                        .map(|&i| {
                            if i % 3 == 0 {
                                det::sin(mode as f64 * PI * mesh.nodes_m[i / 3][0])
                            } else {
                                0.0
                            }
                        })
                        .collect::<Vec<_>>();
                    let found = modes(&a, (0.95 * exact, 1.05 * exact));
                    let selected = found
                        .iter()
                        .max_by(|a1, b1| {
                            overlap(&a, &a1.phi, &reference)
                                .total_cmp(&overlap(&a, &b1.phi, &reference))
                        })
                        .unwrap();
                    assert!(overlap(&a, &selected.phi, &reference) > 0.99);
                    let error = (selected.lambda / exact - 1.0).abs();
                    assert!(selected.residual < 0.02 * (selected.lambda - exact).abs());
                    errors[mode - 1].push(error);
                    log(
                        "g1_roller_block_exact_p_modes_order",
                        "G1",
                        family_name(family),
                        &[1.0 / nx as f64],
                        &[a.free_dofs.len()],
                        &[error],
                        &[],
                        &[],
                        0.3,
                        started.elapsed().as_secs_f64(),
                    );
                }
            }
            for (mode, error) in errors.iter().enumerate() {
                let ladder = fs_mms::RefinementLadder::new(hs.clone(), error.clone()).unwrap();
                let verdict = fs_mms::OrderGate { theoretical: 4.0 }
                    .check(
                        &format!("{}-roller-P{}", family_name(family), mode + 1),
                        fs_mms::LadderSide::Primal,
                        &ladder,
                    )
                    .unwrap();
                assert!((verdict.fit.observed - 4.0).abs() <= 0.3);
                println!("{}", verdict.json_line(true));
            }
        }
    });
}

fn hyperbolic(x: f64) -> (f64, f64) {
    let e = det::exp(x);
    ((e + 1.0 / e) / 2.0, (e - 1.0 / e) / 2.0)
}
/// Free-free Timoshenko characteristic from zero end shear and moment.
/// Khasawneh & Segalman (2019), Table 7, arXiv:1811.03222, equivalent
/// cosh-scaled form. r=I/(AL²), s=EI/(kGA L²), q=rho A L⁴ omega²/(EI).
fn timo_parameters(beta: f64, slenderness: f64, nu: f64) -> (f64, f64, f64) {
    let r = 1.0 / (12.0 * slenderness * slenderness);
    let kappa = 10.0 * (1.0 + nu) / (12.0 + 11.0 * nu);
    let s = r * 2.0 * (1.0 + nu) / kappa;
    let q = beta.powi(4);
    let discriminant = det::sqrt(q * q * (r - s) * (r - s) + 4.0 * q);
    let b2 = (discriminant + q * (r + s)) / 2.0;
    let a2 = 2.0 * q * (1.0 - q * r * s) / (discriminant + q * (r + s));
    let eta = (a2 + q * s) / (b2 - q * s);
    (det::sqrt(a2), det::sqrt(b2), eta)
}
fn timo_characteristic(beta: f64, slenderness: f64, nu: f64) -> f64 {
    let (a, b, eta) = timo_parameters(beta, slenderness, nu);
    let (ch, sh) = hyperbolic(a);
    let ratio = eta * a / b;
    det::cos(b) - 1.0 / ch + 0.5 * (ratio - 1.0 / ratio) * sh / ch * det::sin(b)
}
fn roots(characteristic: impl Fn(f64) -> f64) -> [f64; 3] {
    core::array::from_fn(|i| {
        let mut beta = (2.0 * (i + 1) as f64 + 1.0) * PI / 2.0;
        for _ in 0..32 {
            let h = 1e-5;
            let f = characteristic(beta);
            let derivative = (characteristic(beta + h) - characteristic(beta - h)) / (2.0 * h);
            assert!(derivative.is_finite() && derivative.abs() > 0.1);
            let step = f / derivative;
            beta -= step;
            if step.abs() < 2e-14 {
                break;
            }
        }
        assert!(characteristic(beta).abs() < 1e-12);
        beta
    })
}
fn timo_oracle() -> ([f64; 3], [f64; 3]) {
    let beta = roots(|b| timo_characteristic(b, 50.0, 0.3));
    let eb = roots(|b| {
        let (ch, _) = hyperbolic(b);
        det::cos(b) - 1.0 / ch
    });
    for &root in &beta {
        // Independent free-end shear/moment residuals of the four spatial
        // functions, rather than self-verifying only the Newton scalar.
        let (a, b, eta) = timo_parameters(root, 50.0, 0.3);
        let (ch, sh) = hyperbolic(a);
        let coefficient = -(sh + eta * a / b * det::sin(b)) / (ch - det::cos(b));
        let c = eta;
        let d = b / a * coefficient;
        let shear = sh / a + coefficient * ch / a + c * det::sin(b) / b - d * det::cos(b) / b;
        let moment = ch + coefficient * sh - c / eta * det::cos(b) - d / eta * det::sin(b);
        assert!(shear.abs() < 1e-10 * (sh / a).abs().max(1.0));
        assert!(moment.abs() < 1e-10 * ch);
        println!(
            "{{\"suite\":\"irp-s1\",\"stage\":\"timoshenko-root\",\"beta\":{root},\"characteristic\":{},\"end_shear_residual\":{shear},\"end_moment_residual\":{moment},\"cowper_kappa\":{}}}",
            timo_characteristic(root, 50.0, 0.3),
            10.0 * 1.3 / (12.0 + 11.0 * 0.3)
        );
    }
    let scale = 0.02 / det::sqrt(12.0);
    (beta.map(|b| b * b * scale), eb.map(|b| b * b * scale))
}
fn bar_modes(a: &TetElasticAssembly, nodes: &[[f64; 3]]) -> ([f64; 7], [f64; 7]) {
    let flex = modes(a, (1e-5, 1.0));
    assert!(flex.len() >= 6);
    let eb_roots = roots(|b| {
        let (ch, _) = hyperbolic(b);
        det::cos(b) - 1.0 / ch
    });
    for (pair, &root) in flex[..6].chunks_exact(2).zip(&eb_roots) {
        let (ch, sh) = hyperbolic(root);
        let sigma = (ch - det::cos(root)) / (sh - det::sin(root));
        let mut projection_trace = 0.0;
        for component in [1, 2] {
            let reference = a
                .free_dofs
                .iter()
                .map(|&i| {
                    if i % 3 == component {
                        let x = root * nodes[i / 3][0];
                        let (c, s) = hyperbolic(x);
                        c + det::cos(x) - sigma * (s + det::sin(x))
                    } else {
                        0.0
                    }
                })
                .collect::<Vec<_>>();
            for mode in pair {
                projection_trace += overlap(a, &mode.phi, &reference).powi(2);
            }
        }
        assert!(
            projection_trace > 1.8,
            "both independent bending planes must be present: {projection_trace}"
        );
    }
    let mut values = [0.0; 7];
    let mut bounds = [0.0; 7];
    for i in 0..6 {
        values[i] = flex[i].lambda;
        bounds[i] = flex[i].residual;
    }
    // First longitudinal mode is selected by mass overlap, not by its index
    // among interleaved higher flexural/torsional modes.
    let target = PI * PI;
    let reference = a
        .free_dofs
        .iter()
        .map(|&i| {
            if i % 3 == 0 {
                det::cos(PI * nodes[i / 3][0])
            } else {
                0.0
            }
        })
        .collect::<Vec<_>>();
    let window = modes(a, (0.95 * target, 1.05 * target));
    let axial = window
        .iter()
        .max_by(|a1, b1| {
            overlap(a, &a1.phi, &reference).total_cmp(&overlap(a, &b1.phi, &reference))
        })
        .unwrap();
    assert!(overlap(a, &axial.phi, &reference) > 0.99);
    values[6] = axial.lambda;
    bounds[6] = axial.residual;
    (values, bounds)
}
pub fn slender_bar() {
    let started = std::time::Instant::now();
    let (timo, eb) = timo_oracle();
    with_cx(|cx| {
        for family in FAMILIES {
            let mut values = Vec::new();
            let mut residuals = Vec::new();
            let mut dofs = Vec::new();
            let mut hs = Vec::new();
            for r in [1, 2, 4] {
                let mesh = mesh([16 * r, r, r], [1.0, 0.02, 0.02], family, cx);
                let a = assembly(&mesh, &material(1.0, 0.3), &[], cx);
                let (row, bounds) = bar_modes(&a, &mesh.nodes_m);
                values.push(row);
                residuals.push(bounds);
                dofs.push(a.free_dofs.len());
                hs.push(1.0 / (16 * r) as f64);
                let errors = core::array::from_fn::<_, 7, _>(|i| {
                    let exact = if i < 6 { timo[i / 2] } else { PI };
                    (det::sqrt(row[i]) / exact - 1.0).abs()
                });
                let eb_errors = core::array::from_fn::<_, 6, _>(|i| {
                    (det::sqrt(row[i]) / eb[i / 2] - 1.0).abs()
                });
                log(
                    "g2_slender_square_bar_free_free_vs_timoshenko",
                    "G2",
                    family_name(family),
                    &[1.0 / (16 * r) as f64],
                    &[a.free_dofs.len()],
                    &errors,
                    &eb_errors,
                    &[],
                    0.001,
                    started.elapsed().as_secs_f64(),
                );
            }
            let mut errors = Vec::new();
            let mut orders = Vec::new();
            let mut estimates = Vec::new();
            for i in 0..7 {
                let row = core::array::from_fn(|r| values[r][i]);
                let bounds = core::array::from_fn(|r| residuals[r][i]);
                let (_limit, estimate, order) = richardson(row, bounds);
                let finest = values[2][i];
                assert!(
                    estimate.is_finite() && estimate >= 0.0 && estimate < finest,
                    "Richardson estimate must resolve a positive eigenvalue"
                );
                let relative_frequency_estimate = estimate / (2.0 * (finest - estimate));
                assert!(
                    relative_frequency_estimate <= 0.0002,
                    "{} mode {i}: Richardson frequency estimate {relative_frequency_estimate:e}, order {order}",
                    family_name(family)
                );
                let exact = if i < 6 { timo[i / 2] } else { PI };
                let error = (det::sqrt(finest) / exact - 1.0).abs();
                assert!(
                    error <= 0.001,
                    "{} mode {i}: {error:e}",
                    family_name(family)
                );
                errors.push(error);
                orders.push(order);
                estimates.push(relative_frequency_estimate);
            }
            println!(
                "{{\"suite\":\"irp-s1\",\"stage\":\"bar-discretization-estimates\",\"family\":\"{}\",\"relative_frequency_estimates\":{estimates:?},\"tolerance\":0.0002}}",
                family_name(family)
            );
            log(
                "g2_slender_square_bar_free_free_vs_timoshenko",
                "G2",
                family_name(family),
                &hs,
                &dofs,
                &errors,
                &[],
                &orders,
                0.001,
                started.elapsed().as_secs_f64(),
            );
        }
    });
}
pub fn cross_family() {
    let started = std::time::Instant::now();
    with_cx(|cx| {
        let mut hex_values = Vec::new();
        let mut tet_values = Vec::new();
        let mut hex_bounds = Vec::new();
        let mut tet_bounds = Vec::new();
        let mut dofs = Vec::new();
        for r in [1, 2, 4] {
            let extent = [1.0, 0.1, 0.07];
            let h = mesh([6 * r, r, r], extent, HexFamily::Lagrange27, cx);
            let ha = assembly(&h, &material(1.0, 0.25), &[], cx);
            let hm = modes(&ha, (1e-5, 12.0));
            assert!(hm.len() >= 6);
            hex_values.push(core::array::from_fn::<_, 6, _>(|i| hm[i].lambda));
            hex_bounds.push(core::array::from_fn::<_, 6, _>(|i| hm[i].residual));
            let linear = tet_support::block([6 * r, r, r], extent);
            let (nodes, elements) = tet_support::promote(&linear, cx);
            let ta = tet_support::assembly(&nodes, &elements, &material(1.0, 0.25), &[], cx);
            let tm = modes(&ta, (1e-5, 12.0));
            assert!(tm.len() >= 6);
            tet_values.push(core::array::from_fn::<_, 6, _>(|i| tm[i].lambda));
            tet_bounds.push(core::array::from_fn::<_, 6, _>(|i| tm[i].residual));
            dofs.push(ha.free_dofs.len());
        }
        let mut differences = Vec::new();
        let mut orders = Vec::new();
        let mut relative_allowance = 0.0f64;
        for i in 0..6 {
            let (_, he, hp) = richardson(
                core::array::from_fn(|r| hex_values[r][i]),
                core::array::from_fn(|r| hex_bounds[r][i]),
            );
            let (_, te, tp) = richardson(
                core::array::from_fn(|r| tet_values[r][i]),
                core::array::from_fn(|r| tet_bounds[r][i]),
            );
            let difference = (hex_values[2][i] - tet_values[2][i]).abs();
            assert!(
                difference <= he + te,
                "mode {i}: {difference:e} > HEX {he:e} + TET {te:e}"
            );
            let scale = hex_values[2][i];
            differences.push(difference / scale);
            relative_allowance = relative_allowance.max((he + te) / scale);
            println!(
                "{{\"suite\":\"irp-s1\",\"stage\":\"hex27-tet10-comparison\",\"mode\":{i},\"absolute_difference\":{difference},\"hex_discretization_estimate\":{he},\"tet_discretization_estimate\":{te},\"difference_over_estimate_sum\":{}}}",
                difference / (he + te)
            );
            orders.extend([hp, tp]);
        }
        log(
            "g3_hex27_vs_tet10_same_bar",
            "G3",
            "hex27-vs-tet10",
            &[1.0, 0.5, 0.25],
            &dofs,
            &differences,
            &[],
            &orders,
            relative_allowance,
            started.elapsed().as_secs_f64(),
        );
    });
}
