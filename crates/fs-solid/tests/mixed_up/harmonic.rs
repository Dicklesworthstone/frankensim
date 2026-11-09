//! G1/G3: full mixed forced equations, reactive pressure and support work.
use super::*;
use fs_math::c64::C64;
use fs_modal::harmonic::{HarmonicResponseControls, HarmonicResponseError};
use fs_solid::hyper3::mixed::vibration::harmonic::MixedHarmonicError;

fn controls() -> HarmonicResponseControls {
    HarmonicResponseControls {
        maximum_coordinates: 128,
        maximum_payload_bytes: 8 * 1024 * 1024,
        maximum_relative_residual: 2e-8,
    }
}
fn complex_near(a: C64, b: C64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:?} != {b:?}");
}
fn mesh(cx: &Cx<'_>) -> (Vec<[f64; 3]>, Vec<[usize; 10]>) {
    promote_p1_to_p2(
        &[[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        &[[0, 1, 2, 3]],
        None,
        cx,
    )
    .unwrap()
}

#[test]
fn finite_bulk_forcing_retains_loaded_tangent_pressure_and_cross_support_inertia() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(cx);
        let (mu, bulk, rho, stretch) = (2.0, 20.0, 7.0, 1.15);
        let material = law(mu, bulk);
        let initial: Vec<_> = nodes
            .iter()
            .map(|x| [(stretch - 1.0) * x[0], 0.0, 0.0])
            .collect();
        let supports: Vec<_> = (0..3 * nodes.len())
            .filter(|&d| d != 3)
            .map(|d| (d, initial[d / 3][d % 3]))
            .collect();
        let mut p = problem(
            &nodes,
            &elements,
            &material,
            &supports,
            &[],
            VolumetricTreatment::OriginalMaterial,
        );
        p.settings.load_steps = 1;
        let spectrum = p
            .solve_vibration_from(
                &initial,
                &[bulk * (stretch - 1.0); 4],
                &[rho],
                modal_settings(),
                cx,
            )
            .unwrap();
        let f = [stretch, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        // Exact affine P2 integrals: integral (dN1/dX1)^2 = 1/10,
        // integral N1² = 1/420. Full material tangent includes current prestress.
        let k = material.tangent(&f).unwrap()[0][0] / 10.0;
        let m = rho / 420.0;
        let mut unit = vec![[0.0; 3]; nodes.len()];
        unit[1][0] = 1.0;
        let mut mass_column = vec![0.0; 3 * nodes.len()];
        for (d, v) in spectrum.reference_acceleration_load_n(&unit, cx).unwrap() {
            mass_column[d] = v;
        }
        assert!(
            mass_column
                .iter()
                .enumerate()
                .any(|(d, &v)| d != 3 && v != 0.0)
        );
        let mut force = vec![C64::ZERO; 3 * nodes.len()];
        force[3] = C64::new(0.002, -0.001);
        force[0] = C64::new(0.03, 0.02);
        for ratio in [0.0_f64, 0.37, 1.37] {
            let omega = ratio * (k / m).sqrt();
            let prepared = spectrum
                .prepare_harmonic_response(omega, controls(), cx)
                .unwrap();
            let response = prepared.solve(&force, cx).unwrap();
            let expected = force[3].scale(1.0 / (k - omega * omega * m));
            complex_near(response.displacement_m[3], expected, 2e-9 * expected.abs());
            complex_near(
                response.velocity_m_s[3],
                C64::new(omega * expected.im, -omega * expected.re),
                2e-9 * (1.0 + omega) * expected.abs(),
            );
            assert!(response.momentum_relative_residual <= controls().maximum_relative_residual);
            assert!(
                response.pressure_relative_residual <= modal_settings().pressure_residual_tolerance
            );
            for (j, &pressure) in response.pressure_pa.iter().enumerate() {
                complex_near(
                    pressure,
                    expected.scale(spectrum.pressure_response_pa_per_m()[j]),
                    2e-9 * (1.0 + pressure.abs()),
                );
            }
            for &(d, reaction) in &response.support_reaction_n {
                let mut expected_reaction = expected
                    .scale(spectrum.mixed_tangent().get(d, 3) - omega * omega * mass_column[d])
                    - force[d];
                for (j, &pressure) in response.pressure_pa.iter().enumerate() {
                    expected_reaction = expected_reaction
                        + pressure.scale(spectrum.mixed_tangent().get(d, 3 * nodes.len() + j));
                }
                complex_near(
                    reaction,
                    expected_reaction,
                    2e-10 * (1.0 + expected_reaction.abs()),
                );
                assert_eq!(response.displacement_m[d], C64::ZERO);
            }
            assert!(response.average_input_power_w.abs() < 1e-14);
            let repeated = prepared.solve(&force, cx).unwrap();
            assert_eq!(response.displacement_m, repeated.displacement_m);
            assert_eq!(response.pressure_pa, repeated.pressure_pa);
        }
    });
}

#[test]
fn exact_hydrostatic_forcing_has_pressure_without_admissible_displacement() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(cx);
        let material = law(2.0, 20.0);
        // Full face clamp retains all four independent pressure coordinates.
        let supports: Vec<_> = nodes
            .iter()
            .enumerate()
            .filter(|(_, x)| x[2] == 0.0)
            .flat_map(|(node, _)| (0..3).map(move |a| (3 * node + a, 0.0)))
            .collect();
        let p = problem(
            &nodes,
            &elements,
            &material,
            &supports,
            &[],
            VolumetricTreatment::ExactIncompressible,
        );
        let spectrum = p
            .solve_incompressible_vibration(&[7.0], modal_settings(), cx)
            .unwrap();
        let nu = 3 * nodes.len();
        let amplitude = C64::new(0.04, -0.07);
        let mut force = vec![C64::ZERO; nu];
        // Uniform incremental pressure is a full hydrostatic traction. Its
        // free force lies in range(B^T), outside the admissible nullspace.
        for (d, force) in force.iter_mut().enumerate() {
            for j in 0..4 {
                *force = *force + amplitude.scale(spectrum.mixed_tangent().get(d, nu + j));
            }
        }
        let omega = 0.31 * spectrum.modes()[0].coordinates.lambda.sqrt();
        let prepared = spectrum
            .prepare_harmonic_response(omega, controls(), cx)
            .unwrap();
        let response = prepared.solve(&force, cx).unwrap();
        assert!(response.displacement_m.iter().all(|u| u.abs() < 1e-12));
        for &pressure in &response.pressure_pa {
            complex_near(pressure, amplitude, 2e-11);
        }
        assert!(
            response
                .support_reaction_n
                .iter()
                .all(|(_, r)| r.abs() < 1e-11)
        );
        assert!(response.average_input_power_w.abs() < 1e-15);
        // A free-space modal-pressure sum alone would incorrectly return zero.
        assert!(response.pressure_pa.iter().all(|p| p.abs() > 0.05));
    });
}

#[test]
fn exact_forced_modes_preserve_gauge_and_all_original_support_reactions() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(cx);
        let nu = 3 * nodes.len();
        let material = law(2.0, 20.0);
        let free = [0, 3, 7, 11];
        let supports: Vec<_> = (0..nu)
            .filter(|d| !free.contains(d))
            .map(|d| (d, 0.0))
            .collect();
        let mut p = problem(
            &nodes,
            &elements,
            &material,
            &supports,
            &[],
            VolumetricTreatment::ExactIncompressible,
        );
        p.pressure_gauge = Some((0, 0.0));
        p.settings.load_steps = 1;
        let spectrum = p
            .solve_incompressible_vibration_with_reference_density(
                &[[3.0, 15.0, 3.0, 3.0]],
                modal_settings(),
                cx,
            )
            .unwrap();
        let omega = 0.43 * spectrum.modes()[0].coordinates.lambda.sqrt();
        let displacement = [
            C64::new(0.001, -0.002),
            C64::new(-0.001, 0.002),
            C64::new(-0.001, 0.002),
            C64::new(-0.001, 0.002),
        ];
        let pressure = [
            C64::ZERO,
            C64::new(0.01, 0.03),
            C64::new(-0.04, 0.01),
            C64::new(0.02, -0.03),
        ];
        let mut real = vec![[0.0; 3]; nodes.len()];
        let mut imaginary = real.clone();
        for (&d, &u) in free.iter().zip(&displacement) {
            real[d / 3][d % 3] = u.re;
            imaginary[d / 3][d % 3] = u.im;
        }
        let mut inertia = vec![C64::ZERO; nu];
        for (d, v) in spectrum.reference_acceleration_load_n(&real, cx).unwrap() {
            inertia[d].re = v;
        }
        for (d, v) in spectrum
            .reference_acceleration_load_n(&imaginary, cx)
            .unwrap()
        {
            inertia[d].im = v;
        }
        let mut force = vec![C64::ZERO; nu];
        let mut expected_reactions = Vec::new();
        for d in 0..nu {
            let mut action = -inertia[d].scale(omega * omega);
            for (&j, &u) in free.iter().zip(&displacement) {
                action = action + u.scale(spectrum.mixed_tangent().get(d, j));
            }
            for (j, &pressure) in pressure.iter().enumerate() {
                action = action + pressure.scale(spectrum.mixed_tangent().get(d, nu + j));
            }
            if free.contains(&d) {
                force[d] = action;
            } else {
                force[d] = C64::new(0.0001 * (d + 1) as f64, -0.0002);
                expected_reactions.push((d, action - force[d]));
            }
        }
        let prepared = spectrum
            .prepare_harmonic_response(omega, controls(), cx)
            .unwrap();
        let response = prepared.solve(&force, cx).unwrap();
        for (&d, &u) in free.iter().zip(&displacement) {
            complex_near(response.displacement_m[d], u, 2e-11);
        }
        for (&actual, &expected) in response.pressure_pa.iter().zip(&pressure) {
            complex_near(actual, expected, 2e-10);
        }
        assert_eq!(response.pressure_pa[0], C64::ZERO);
        for ((d, a), (e, b)) in response.support_reaction_n.iter().zip(expected_reactions) {
            assert_eq!(*d, e);
            complex_near(*a, b, 2e-10);
        }
        assert!(
            response.pressure_relative_residual <= modal_settings().pressure_residual_tolerance
        );
        assert!(response.momentum_relative_residual <= controls().maximum_relative_residual);
        assert!(response.average_input_power_w.abs() < 1e-14);
    });
}

#[test]
fn fully_constrained_pressure_response_and_admission_refusals() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(cx);
        let nu = 3 * nodes.len();
        let material = law(2.0, 20.0);
        // These four shape divergences span P1, leaving no admissible motion.
        let free = [3, 7, 11, 13];
        let supports: Vec<_> = (0..nu)
            .filter(|d| !free.contains(d))
            .map(|d| (d, 0.0))
            .collect();
        let p = problem(
            &nodes,
            &elements,
            &material,
            &supports,
            &[],
            VolumetricTreatment::ExactIncompressible,
        );
        let spectrum = p
            .solve_incompressible_vibration(&[7.0], modal_settings(), cx)
            .unwrap();
        assert_eq!(spectrum.dimension(), 0);
        let prepared = spectrum
            .prepare_harmonic_response(3.0, controls(), cx)
            .unwrap();
        let mut force = vec![C64::ZERO; nu];
        force[3] = C64::new(0.01, 0.02);
        let response = prepared.solve(&force, cx).unwrap();
        assert!(response.projected_response.is_none());
        assert!(response.displacement_m.iter().all(|u| *u == C64::ZERO));
        assert!(response.pressure_pa.iter().any(|p| p.abs() > 1e-3));
        let small = HarmonicResponseControls {
            maximum_coordinates: 3,
            ..controls()
        };
        assert!(matches!(
            spectrum.prepare_harmonic_response(3.0, small, cx),
            Err(MixedHarmonicError::Harmonic(
                HarmonicResponseError::Budget { .. }
            ))
        ));
        let small = HarmonicResponseControls {
            maximum_payload_bytes: 0,
            ..controls()
        };
        assert!(matches!(
            spectrum.prepare_harmonic_response(3.0, small, cx),
            Err(MixedHarmonicError::Harmonic(
                HarmonicResponseError::Budget { .. }
            ))
        ));
        assert!(prepared.solve(&force[..nu - 1], cx).is_err());
        force[3].re = f64::NAN;
        assert!(prepared.solve(&force, cx).is_err());
        force[3].re = 0.01;
        assert!(
            spectrum
                .prepare_harmonic_response(f64::INFINITY, controls(), cx)
                .is_err()
        );
        let gate = CancelGate::new();
        gate.request();
        with_gate(&gate, |cancelled| {
            assert!(prepared.solve(&force, cancelled).is_err());
            assert!(
                spectrum
                    .prepare_harmonic_response(3.0, controls(), cancelled)
                    .is_err()
            );
        });
        assert!(
            prepared.solve(&force, cx).is_ok(),
            "cancelled use must preserve prepared factor"
        );
    });
}
