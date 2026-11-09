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

fn full_mass(nu: usize, action: impl Fn(&[[f64; 3]]) -> Vec<(usize, f64)>) -> Vec<f64> {
    let mut matrix = vec![0.0; nu * nu];
    for j in 0..nu {
        let mut unit = vec![[0.0; 3]; nu / 3];
        unit[j / 3][j % 3] = 1.0;
        for (i, value) in action(&unit) {
            matrix[i * nu + j] = value;
        }
    }
    matrix
}

// Independent original mixed KKT solve. It does not use pressure condensation,
// a constraint nullspace, modal forces or the adapter's affine lift.
#[allow(clippy::too_many_arguments)]
fn original_kkt(
    tangent: &fs_sparse::Csr,
    mass: &[f64],
    free: &[usize],
    pressures: &[usize],
    omega: f64,
    force: &[C64],
    supports: &[(usize, C64)],
) -> (Vec<C64>, Vec<C64>, Vec<(usize, C64)>) {
    let nu = force.len();
    let np = tangent.nrows() - nu;
    let indices: Vec<_> = free
        .iter()
        .copied()
        .chain(pressures.iter().map(|&p| nu + p))
        .collect();
    let n = indices.len();
    let entry = |i, j| {
        tangent.get(i, j)
            - if i < nu && j < nu {
                omega * omega * mass[i * nu + j]
            } else {
                0.0
            }
    };
    let mut matrix = vec![0.0; n * n];
    let mut real = vec![0.0; n];
    let mut imaginary = real.clone();
    for (i, &row) in indices.iter().enumerate() {
        let mut rhs = if row < nu { force[row] } else { C64::ZERO };
        for &(column, value) in supports {
            rhs = rhs - value.scale(entry(row, column));
        }
        real[i] = rhs.re;
        imaginary[i] = rhs.im;
        for (j, &col) in indices.iter().enumerate() {
            matrix[i * n + j] = entry(row, col);
        }
    }
    let factor = fs_la::factor::lu(&matrix, n).unwrap();
    factor.solve(&mut real);
    factor.solve(&mut imaginary);
    let mut displacement = vec![C64::ZERO; nu];
    let mut pressure = vec![C64::ZERO; np];
    for &(d, value) in supports {
        displacement[d] = value;
    }
    for (i, &d) in free.iter().enumerate() {
        displacement[d] = C64::new(real[i], imaginary[i]);
    }
    for (i, &p) in pressures.iter().enumerate() {
        pressure[p] = C64::new(real[free.len() + i], imaginary[free.len() + i]);
    }
    let reactions = (0..nu)
        .filter(|d| !free.contains(d))
        .map(|d| {
            let mut reaction = -force[d];
            for (j, &value) in displacement.iter().enumerate() {
                reaction = reaction + value.scale(entry(d, j));
            }
            for (p, &value) in pressure.iter().enumerate() {
                reaction = reaction + value.scale(tangent.get(d, nu + p));
            }
            (d, reaction)
        })
        .collect();
    (displacement, pressure, reactions)
}

fn check_kkt(
    response: &fs_solid::hyper3::mixed::vibration::harmonic::MixedHarmonicResponse,
    expected: (Vec<C64>, Vec<C64>, Vec<(usize, C64)>),
) {
    for (&a, &b) in response.displacement_m.iter().zip(&expected.0) {
        complex_near(a, b, 2e-10 * (1.0 + b.abs()));
    }
    for (&a, &b) in response.pressure_pa.iter().zip(&expected.1) {
        complex_near(a, b, 2e-9 * (1.0 + b.abs()));
    }
    for (&(i, a), &(j, b)) in response.support_reaction_n.iter().zip(&expected.2) {
        assert_eq!(i, j);
        complex_near(a, b, 2e-9 * (1.0 + b.abs()));
    }
    assert!(response.power_balance_residual_w.abs() < 2e-11);
}

#[test]
fn finite_bulk_base_excitation_retains_pressure_lift_and_original_cross_mass() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(cx);
        let nu = 3 * nodes.len();
        let material = law(2.0, 20.0);
        for free in [&[3][..], &[3, 7, 11, 13][..]] {
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
                VolumetricTreatment::OriginalMaterial,
            );
            let spectrum = p
                .solve_vibration_with_reference_density(
                    &[[3.0, 15.0, 3.0, 3.0]],
                    modal_settings(),
                    cx,
                )
                .unwrap();
            let mass = full_mass(nu, |a| {
                spectrum.reference_acceleration_load_n(a, cx).unwrap()
            });
            let motion = [
                (0, C64::new(0.001, -0.002)),
                (12, C64::new(-0.003, 0.002)),
                (21, C64::new(0.004, 0.001)),
            ];
            let mut force = vec![C64::ZERO; nu];
            force[3] = C64::new(0.003, 0.004);
            force[0] = C64::new(-0.003, 0.002);
            for ratio in [0.0_f64, 0.4, 1.4] {
                let omega = ratio * spectrum.modes()[0].displacement.lambda.sqrt();
                let prepared = spectrum
                    .prepare_harmonic_response(omega, controls(), cx)
                    .unwrap();
                let response = prepared
                    .solve_with_support_motion(&force, &motion, cx)
                    .unwrap();
                check_kkt(
                    &response,
                    original_kkt(
                        spectrum.mixed_tangent(),
                        &mass,
                        spectrum.free_dofs(),
                        &[0, 1, 2, 3],
                        omega,
                        &force,
                        &motion,
                    ),
                );
                assert!(response.pressure_pa.iter().any(|p| p.abs() > 1e-3));
                for &(d, value) in &motion {
                    assert_eq!(response.displacement_m[d], value);
                }
                if omega > 0.0 {
                    assert!(response.average_support_power_w.abs() > 1e-7);
                    assert!(response.average_input_power_w.abs() > 1e-7);
                }
                let again = prepared
                    .solve_with_support_motion(&force, &motion, cx)
                    .unwrap();
                assert_eq!(response.displacement_m, again.displacement_m);
                assert_eq!(response.pressure_pa, again.pressure_pa);
                assert_eq!(
                    prepared.solve(&force, cx).unwrap().displacement_m,
                    prepared
                        .solve_with_support_motion(&force, &[], cx)
                        .unwrap()
                        .displacement_m
                );
            }
        }
    });
}

#[test]
fn exact_base_motion_lifts_weak_volume_including_zero_dimensional_motion_space() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(cx);
        let nu = 3 * nodes.len();
        let material = law(2.0, 20.0);
        for only_four in [false, true] {
            let supports: Vec<_> = (0..nu)
                .filter(|&d| {
                    if only_four {
                        ![3, 7, 11, 13].contains(&d)
                    } else {
                        nodes[d / 3][2] == 0.0
                    }
                })
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
                .solve_incompressible_vibration_with_reference_density(
                    &[[3.0, 15.0, 3.0, 3.0]],
                    modal_settings(),
                    cx,
                )
                .unwrap();
            let mass = full_mass(nu, |a| {
                spectrum.reference_acceleration_load_n(a, cx).unwrap()
            });
            // A moving quadratic edge shape has nonzero volume work, so the
            // free particular displacement cannot be dropped by Z^T.
            let motion = if only_four {
                vec![(14, C64::new(0.003, -0.004))]
            } else {
                vec![(13, C64::new(0.003, -0.004))]
            };
            let mut force = vec![C64::ZERO; nu];
            force[spectrum.free_dofs()[0]] = C64::new(0.002, 0.001);
            force[motion[0].0] = C64::new(0.001, 0.002);
            for omega in [0.0, 0.7, 2.3] {
                let prepared = spectrum
                    .prepare_harmonic_response(omega, controls(), cx)
                    .unwrap();
                let response = prepared
                    .solve_with_support_motion(&force, &motion, cx)
                    .unwrap();
                check_kkt(
                    &response,
                    original_kkt(
                        spectrum.mixed_tangent(),
                        &mass,
                        spectrum.free_dofs(),
                        spectrum.independent_pressure_indices(),
                        omega,
                        &force,
                        &motion,
                    ),
                );
                let mut maximum_constraint_action = 0.0_f64;
                for p in 0..4 {
                    let action = spectrum.free_dofs().iter().fold(C64::ZERO, |s, &d| {
                        s + response.displacement_m[d]
                            .scale(spectrum.mixed_tangent().get(nu + p, d))
                    });
                    maximum_constraint_action = maximum_constraint_action.max(action.abs());
                }
                assert!(maximum_constraint_action > 1e-5);
                assert_eq!(response.projected_response.is_none(), only_four);
                if only_four {
                    assert!(
                        spectrum
                            .free_dofs()
                            .iter()
                            .any(|&d| response.displacement_m[d].abs() > 1e-3)
                    );
                }
            }
        }
    });
}

#[test]
fn gauge_base_motion_preserves_all_constraints_and_rejects_incompatible_volume() {
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
        let spectrum = p
            .solve_incompressible_vibration(&[7.0], modal_settings(), cx)
            .unwrap();
        let mass = full_mass(nu, |a| {
            spectrum.reference_acceleration_load_n(a, cx).unwrap()
        });
        let omega = 0.37 * spectrum.modes()[0].coordinates.lambda.sqrt();
        let prepared = spectrum
            .prepare_harmonic_response(omega, controls(), cx)
            .unwrap();
        let translation = [
            C64::new(0.002, -0.003),
            C64::new(-0.001, 0.001),
            C64::new(0.003, 0.002),
        ];
        let expected: Vec<_> = (0..nu).map(|d| translation[d % 3]).collect();
        let motion: Vec<_> = supports.iter().map(|&(d, _)| (d, expected[d])).collect();
        let pressure = [
            C64::ZERO,
            C64::new(0.04, -0.01),
            C64::new(-0.03, 0.02),
            C64::new(0.02, 0.01),
        ];
        let force: Vec<_> = (0..nu)
            .map(|d| {
                let mut f = C64::ZERO;
                for j in 0..nu {
                    f = f + expected[j].scale(
                        spectrum.mixed_tangent().get(d, j) - omega * omega * mass[d * nu + j],
                    );
                }
                for (j, &q) in pressure.iter().enumerate() {
                    f = f + q.scale(spectrum.mixed_tangent().get(d, nu + j));
                }
                f
            })
            .collect();
        let response = prepared
            .solve_with_support_motion(&force, &motion, cx)
            .unwrap();
        for (&actual, &expected) in response.displacement_m.iter().zip(&expected) {
            complex_near(actual, expected, 2e-11);
        }
        for (&actual, &expected) in response.pressure_pa.iter().zip(&pressure) {
            complex_near(actual, expected, 2e-10);
        }
        assert_eq!(response.pressure_pa[0], C64::ZERO);
        assert!(
            response
                .support_reaction_n
                .iter()
                .all(|(_, f)| f.abs() < 1e-10)
        );
        assert!(matches!(
            prepared.solve_with_support_motion(
                &vec![C64::ZERO; nu],
                &[(13, C64::new(0.01, 0.0))],
                cx
            ),
            Err(MixedHarmonicError::Residual {
                equation: "all original pressure equations",
                ..
            })
        ));
        for invalid in [
            vec![(free[0], C64::ZERO)],
            vec![(nu, C64::ZERO)],
            vec![(1, C64::ZERO), (1, C64::ZERO)],
            vec![(2, C64::ZERO), (1, C64::ZERO)],
            vec![(1, C64::new(f64::NAN, 0.0))],
        ] {
            assert!(
                prepared
                    .solve_with_support_motion(&force, &invalid, cx)
                    .is_err()
            );
        }
        let gate = CancelGate::new();
        gate.request();
        with_gate(&gate, |cancelled| {
            assert!(
                prepared
                    .solve_with_support_motion(&force, &motion, cancelled)
                    .is_err()
            )
        });
        assert_eq!(
            response.displacement_m,
            prepared
                .solve_with_support_motion(&force, &motion, cx)
                .unwrap()
                .displacement_m
        );
    });
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
