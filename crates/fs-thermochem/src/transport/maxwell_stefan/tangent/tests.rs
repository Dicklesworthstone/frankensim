use super::super::{MolarMass, SpeciesId};
use super::*;

fn owner(masses: &[f64]) -> MaxwellStefanV1 {
    MaxwellStefanV1::new(
        masses
            .iter()
            .enumerate()
            .map(|(i, &mass)| {
                (
                    SpeciesId::new(&format!("T{i:03}")).unwrap(),
                    MolarMass::new(mass),
                )
            })
            .collect(),
    )
    .unwrap()
}
fn pairs(n: usize) -> Vec<f64> {
    let mut d = vec![0.0; n * n];
    for i in 0..n {
        for j in i + 1..n {
            let value = 1e-5 * (1 + (3 * i + 7 * j) % 17) as f64;
            d[i * n + j] = value;
            d[j * n + i] = value;
        }
    }
    d
}
fn direction<'a>(
    dx: &'a [f64],
    dd: &'a [f64],
    dc: f64,
    d_pair: &'a [f64],
) -> MaxwellStefanDirectionV1<'a> {
    MaxwellStefanDirectionV1 {
        mole_fraction_derivatives: dx,
        driving_force_derivatives_per_m: dd,
        concentration_derivative_mol_m3: dc,
        binary_diffusivity_derivatives_m2_s: d_pair,
    }
}
fn near(actual: f64, expected: f64, relative: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs() + absolute,
        "actual={actual:e}, expected={expected:e}"
    );
}
fn run(
    owner: &MaxwellStefanV1,
    x: &[f64],
    d: &[f64],
    c: f64,
    pairs: &[f64],
    derivative: MaxwellStefanDirectionV1<'_>,
    mass: bool,
    work: &mut MaxwellStefanTangentWorkspaceV1,
) -> (Vec<f64>, Vec<f64>) {
    if mass {
        let r = owner
            .evaluate_driving_forces_tangent(
                x,
                d,
                AmountConcentration::new(c),
                pairs,
                derivative,
                MaxwellStefanSolveV1::default(),
                work,
                || false,
            )
            .unwrap();
        assert!(r.maximum_gradient_derivative_residual_per_m().is_finite());
        assert!(r.mass_frame_derivative_residual_kg_m2_s().abs() < 1e-12);
        (
            r.value().mass_fluxes().to_vec(),
            r.mass_flux_derivatives().to_vec(),
        )
    } else {
        let r = owner
            .evaluate_molar_driving_forces_tangent(
                x,
                d,
                AmountConcentration::new(c),
                pairs,
                derivative,
                MaxwellStefanMolarSolveV1::default(),
                work,
                || false,
            )
            .unwrap();
        assert!(r.maximum_gradient_derivative_residual_per_m().is_finite());
        assert!(r.molar_frame_derivative_residual_mol_m2_s().abs() < 1e-12);
        (
            r.value().molar_fluxes().to_vec(),
            r.molar_flux_derivatives().to_vec(),
        )
    }
}

#[test]
fn g1_binary_nonideal_direction_matches_analytic_thermodynamic_response_in_both_frames() {
    let masses = [0.018, 0.046];
    let law = owner(&masses);
    let (x0, dx0, chi, dchi, gradient, dgradient) = (0.25, 0.03, 1.5, 0.2, 0.2, -0.07);
    let (c, dc, diffusion, ddiffusion) = (18000.0, 500.0, 8e-10, 1e-10);
    let factor = 1.0 - 2.0 * chi * x0 * (1.0 - x0);
    let dfactor = -2.0 * dchi * x0 * (1.0 - x0) - 2.0 * chi * dx0 * (1.0 - 2.0 * x0);
    let force = factor * gradient;
    let dforce = dfactor * gradient + factor * dgradient;
    let x = [x0, 1.0 - x0];
    let d = [force, -force];
    let pair = [0.0, diffusion, diffusion, 0.0];
    let d_pair = [0.0, ddiffusion, ddiffusion, 0.0];
    let dx = [dx0, -dx0];
    let dd = [dforce, -dforce];
    let derivative = direction(&dx, &dd, dc, &d_pair);
    let mean = x0 * masses[0] + (1.0 - x0) * masses[1];
    let dmean = dx0 * (masses[0] - masses[1]);
    let molar = -c * diffusion * force;
    let dmolar = -(dc * diffusion * force + c * ddiffusion * force + c * diffusion * dforce);
    for mass in [false, true] {
        let (value, tangent) = run(
            &law,
            &x,
            &d,
            c,
            &pair,
            derivative,
            mass,
            &mut law.tangent_workspace().unwrap(),
        );
        let (expected, response) = if mass {
            let product = masses[0] * masses[1];
            (
                molar * product / mean,
                product * (dmolar / mean - molar * dmean / (mean * mean)),
            )
        } else {
            (molar, dmolar)
        };
        near(value[0], expected, 4e-14, 0.0);
        near(value[1], -expected, 4e-14, 0.0);
        near(tangent[0], response, 8e-14, 0.0);
        near(tangent[1], -response, 8e-14, 0.0);
    }
}

#[test]
fn g1_original_differentiated_friction_recovers_dense_absent_and_trace_fluxes() {
    for x in [
        vec![0.2, 0.3, 0.5],
        vec![0.5, 0.0, 1e-250, 0.125, 0.125, 0.0, 0.0, 0.25],
        vec![1.0 / 32.0; 32],
        vec![1.0 / 128.0; 128],
    ] {
        let n = x.len();
        let masses: Vec<_> = (0..n).map(|i| 0.005 + 0.002 * (i % 13) as f64).collect();
        let law = owner(&masses);
        let pair = pairs(n);
        let pair_direction: Vec<_> = pair
            .iter()
            .enumerate()
            .map(|(q, &v)| v * (0.03 * ((q / n + q % n) % 5) as f64 - 0.05))
            .collect();
        let mut dx: Vec<_> = (0..n)
            .map(|i| {
                if x[i] == 0.0 {
                    0.01
                } else {
                    0.003 * (i % 5) as f64
                }
            })
            .collect();
        dx[0] = -dx.iter().skip(1).sum::<f64>();
        let (c, dc) = (2.5, -0.7);
        for mass in [false, true] {
            let mut primary: Vec<_> = (0..n).map(|i| ((i % 7) as f64 - 3.0) * 1e-6).collect();
            let mut primary_dot: Vec<_> = (0..n).map(|i| ((i % 9) as f64 - 4.0) * 2e-7).collect();
            primary[0] = -primary.iter().skip(1).sum::<f64>();
            primary_dot[0] = -primary_dot.iter().skip(1).sum::<f64>();
            let z: Vec<_> = primary
                .iter()
                .zip(&masses)
                .map(|(&v, &m)| if mass { v / m } else { v })
                .collect();
            let dz: Vec<_> = primary_dot
                .iter()
                .zip(&masses)
                .map(|(&v, &m)| if mass { v / m } else { v })
                .collect();
            let mut force = vec![0.0; n];
            let mut dforce = vec![0.0; n];
            for i in 0..n {
                // Independent differentiation of the ORIGINAL dimensional rows.
                let mut friction = 0.0;
                let mut response = 0.0;
                for j in 0..n {
                    if i != j {
                        let numerator = x[j] * z[i] - x[i] * z[j];
                        friction += numerator / pair[i * n + j];
                        response += (dx[j] * z[i] + x[j] * dz[i] - dx[i] * z[j] - x[i] * dz[j])
                            / pair[i * n + j]
                            - numerator * pair_direction[i * n + j] / pair[i * n + j].powi(2);
                    }
                }
                force[i] = -friction / c;
                dforce[i] = -response / c + friction * dc / (c * c);
            }
            let derivative = direction(&dx, &dforce, dc, &pair_direction);
            let mut workspace = law.tangent_workspace().unwrap();
            let (actual, actual_dot) =
                run(&law, &x, &force, c, &pair, derivative, mass, &mut workspace);
            for i in 0..n {
                near(actual[i], primary[i], 2e-10, 2e-18);
                near(actual_dot[i], primary_dot[i], 3e-9, 5e-17);
            }
            let repeat = run(&law, &x, &force, c, &pair, derivative, mass, &mut workspace);
            assert_eq!(actual, repeat.0);
            assert_eq!(actual_dot, repeat.1);
        }
    }
}

#[test]
fn g3_interior_parameter_family_central_difference_converges_to_analytic_direction() {
    let law = owner(&[0.005, 0.02, 0.04, 0.08, 0.12]);
    let x = [0.25, 0.125, 0.125, 0.25, 0.25];
    let dx = [-0.07, 0.02, -0.01, 0.04, 0.02];
    let force = [0.2, -0.1, 0.05, -0.07, -0.08];
    let dforce = [-0.04, 0.03, 0.02, -0.01, 0.0];
    let d = pairs(5);
    let dd: Vec<_> = d
        .iter()
        .enumerate()
        .map(|(q, &v)| v * (0.11 * ((q / 5 + q % 5) % 3) as f64 - 0.09))
        .collect();
    let (c, dc) = (3.0, 0.6);
    for mass in [false, true] {
        let (_, response) = run(
            &law,
            &x,
            &force,
            c,
            &d,
            direction(&dx, &dforce, dc, &dd),
            mass,
            &mut law.tangent_workspace().unwrap(),
        );
        let value = |s: f64| {
            let x: Vec<_> = x.iter().zip(dx).map(|(x, dx)| x + s * dx).collect();
            let f: Vec<_> = force.iter().zip(dforce).map(|(f, df)| f + s * df).collect();
            let d: Vec<_> = d.iter().zip(&dd).map(|(d, dd)| d + s * dd).collect();
            let mut work = law.workspace().unwrap();
            if mass {
                law.evaluate_driving_forces(
                    &x,
                    &f,
                    AmountConcentration::new(c + s * dc),
                    &d,
                    MaxwellStefanSolveV1::default(),
                    &mut work,
                    || false,
                )
                .unwrap()
                .mass_fluxes()
                .to_vec()
            } else {
                law.evaluate_molar_driving_forces(
                    &x,
                    &f,
                    AmountConcentration::new(c + s * dc),
                    &d,
                    MaxwellStefanMolarSolveV1::default(),
                    &mut work,
                    || false,
                )
                .unwrap()
                .molar_fluxes()
                .to_vec()
            }
        };
        let errors = [0.01, 0.005, 0.0025].map(|h| {
            let plus = value(h);
            let minus = value(-h);
            plus.iter()
                .zip(minus)
                .zip(&response)
                .map(|((&p, m), r)| ((p - m) / (2.0 * h) - r).abs())
                .fold(0.0, f64::max)
        });
        for ratio in [errors[0] / errors[1], errors[1] / errors[2]] {
            assert!((3.9..4.1).contains(&ratio), "{errors:?}");
        }
        let scale = response.iter().map(|v| v.abs()).fold(0.0, f64::max);
        assert!(errors[2] < 2e-6 * scale);
    }
}

#[test]
fn g1_zero_forces_zero_direction_and_trace_concentration_scaling_keep_physical_limits() {
    let law = owner(&[1.0; 3]);
    let x = [0.2, 0.3, 0.5];
    let d = pairs(3);
    let zeros = [0.0; 3];
    let zero_pairs = [0.0; 9];
    let forcing_direction = direction(&zeros, &[0.2, -0.1, -0.1], 0.0, &zero_pairs);
    for mass in [true, false] {
        let (primal, response) = run(
            &law,
            &x,
            &zeros,
            1.0,
            &d,
            forcing_direction,
            mass,
            &mut law.tangent_workspace().unwrap(),
        );
        assert_eq!(primal, zeros);
        assert!(response.iter().any(|v| *v != 0.0));
    }
    for trace in [1e-100, 1e-200, 1e-300] {
        let x = [1.0, trace, trace];
        let f = [0.0, trace, -trace];
        let d = [0.0, 1e-10, 1e-10, 1e-10, 0.0, 1e-10, 1e-10, 1e-10, 0.0];
        for mass in [true, false] {
            let (value, response) = run(
                &law,
                &x,
                &f,
                1.0,
                &d,
                direction(&zeros, &zeros, 1.0, &zero_pairs),
                mass,
                &mut law.tangent_workspace().unwrap(),
            );
            assert_ne!(response[1], 0.0);
            assert_ne!(response[2], 0.0);
            assert_eq!(response, value);
            let (_, zero) = run(
                &law,
                &x,
                &f,
                1.0,
                &d,
                direction(&zeros, &zeros, 0.0, &zero_pairs),
                mass,
                &mut law.tangent_workspace().unwrap(),
            );
            assert_eq!(zero, zeros);
        }
    }
}

#[test]
fn g0_derivative_sources_are_admitted_before_any_tangent_is_published() {
    let law = owner(&[0.02, 0.03, 0.04]);
    let x = [0.2, 0.3, 0.5];
    let f = [0.2, -0.1, -0.1];
    let d = pairs(3);
    let zero = [0.0; 3];
    let pair_zero = [0.0; 9];
    let mut work = law.tangent_workspace().unwrap();
    let mut diagonal = pair_zero;
    diagonal[0] = 1.0;
    let mut asymmetric = pair_zero;
    asymmetric[1] = 1e-5;
    let mut nan_pair = pair_zero;
    nan_pair[1] = f64::NAN;
    nan_pair[3] = f64::NAN;
    for direction in [
        direction(&zero[..2], &zero, 0.0, &pair_zero),
        direction(&zero, &zero[..2], 0.0, &pair_zero),
        direction(&zero, &zero, 0.0, &pair_zero[..8]),
        direction(&[0.1, 0.0, 0.0], &zero, 0.0, &pair_zero),
        direction(&zero, &[0.1, 0.0, 0.0], 0.0, &pair_zero),
        direction(&[f64::NAN, 0.0, 0.0], &zero, 0.0, &pair_zero),
        direction(&zero, &[f64::INFINITY, 0.0, 0.0], 0.0, &pair_zero),
        direction(&zero, &zero, f64::INFINITY, &pair_zero),
        direction(&zero, &zero, 0.0, &diagonal),
        direction(&zero, &zero, 0.0, &asymmetric),
        direction(&zero, &zero, 0.0, &nan_pair),
    ] {
        assert!(matches!(
            law.evaluate_driving_forces_tangent(
                &x,
                &f,
                AmountConcentration::new(1.0),
                &d,
                direction,
                MaxwellStefanSolveV1::default(),
                &mut work,
                || false,
            ),
            Err(MaxwellStefanErrorV1::Invalid(_))
        ));
    }
    let mut invalid_source = d.clone();
    invalid_source[1] = 0.0;
    assert!(
        law.evaluate_driving_forces_tangent(
            &x,
            &f,
            AmountConcentration::new(1.0),
            &invalid_source,
            direction(&zero, &zero, 0.0, &pair_zero),
            MaxwellStefanSolveV1::default(),
            &mut work,
            || false,
        )
        .is_err()
    );
    run(
        &law,
        &x,
        &f,
        1.0,
        &d,
        direction(&zero, &zero, 0.0, &pair_zero),
        true,
        &mut work,
    );
}

#[test]
fn g3_g4_direction_scaling_cancellation_and_retained_storage_preserve_original_value() {
    let law = owner(&[0.005, 0.02, 0.04, 0.08, 0.12]);
    let x = [0.25, 0.125, 0.125, 0.25, 0.25];
    let dx = [-0.07, 0.02, -0.01, 0.04, 0.02];
    let f = [0.2, -0.1, 0.05, -0.07, -0.08];
    let df = [-0.04, 0.03, 0.02, -0.01, 0.0];
    let d = pairs(5);
    let dd: Vec<_> = d.iter().map(|v| 0.125 * v).collect();
    let original = direction(&dx, &df, 0.6, &dd);
    let mut work = law.tangent_workspace().unwrap();
    let expected_bytes =
        law.workspace().unwrap().workspace_bytes() + 2 * 5 * core::mem::size_of::<f64>();
    assert_eq!(work.workspace_bytes(), expected_bytes);
    let pointers = (
        work.value_mass.as_ptr(),
        work.value_molar.as_ptr(),
        work.solver.matrix.as_ptr(),
    );
    let base = run(&law, &x, &f, 3.0, &d, original, true, &mut work);
    for scale in [-3.0, 0.0, 0.5, 2.0] {
        let dx = dx.map(|v| scale * v);
        let df = df.map(|v| scale * v);
        let dd: Vec<_> = dd.iter().map(|v| scale * v).collect();
        let result = run(
            &law,
            &x,
            &f,
            3.0,
            &d,
            direction(&dx, &df, scale * 0.6, &dd),
            true,
            &mut work,
        );
        assert_eq!(result.0, base.0);
        for (&a, &b) in result.1.iter().zip(&base.1) {
            near(a, scale * b, 3e-12, 1e-22);
        }
        assert_eq!(work.workspace_bytes(), expected_bytes);
        assert_eq!(
            pointers,
            (
                work.value_mass.as_ptr(),
                work.value_molar.as_ptr(),
                work.solver.matrix.as_ptr()
            )
        );
    }
    let mut total = 0;
    law.evaluate_driving_forces_tangent(
        &x,
        &f,
        AmountConcentration::new(3.0),
        &d,
        original,
        MaxwellStefanSolveV1::default(),
        &mut work,
        || {
            total += 1;
            false
        },
    )
    .unwrap();
    for stop in [1, total / 3, 2 * total / 3, total] {
        let mut count = 0;
        assert!(matches!(
            law.evaluate_driving_forces_tangent(
                &x,
                &f,
                AmountConcentration::new(3.0),
                &d,
                original,
                MaxwellStefanSolveV1::default(),
                &mut work,
                || {
                    count += 1;
                    count == stop
                },
            ),
            Err(MaxwellStefanErrorV1::Cancelled)
        ));
        assert_eq!(run(&law, &x, &f, 3.0, &d, original, true, &mut work), base);
    }
}
