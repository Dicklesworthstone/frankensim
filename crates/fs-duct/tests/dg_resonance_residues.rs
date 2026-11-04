//! G1/G3: residues of original DG port transfer, including spatial participation.
use fs_duct::dg::continuous::poles::{PoleOptions, PoleResponseOptions};
use fs_duct::dg::{DgBoundary, DgDuct, DgMesh, DgSidePort, PortId};
use fs_duct::{Duct, Segment};
use fs_material::gas::{GasSpec, GasState};
use fs_math::{c64::C64, det};

const WORK: usize = 1_000_000_000;
fn options(omega: f64) -> PoleOptions {
    PoleOptions {
        spectral_scale: omega,
        step_tolerance: 1e-11,
        residual_tolerance: 1e-10,
        max_iterations: 24,
        max_dimension: 1024,
        max_dense_bytes: 128 * 1024 * 1024,
        max_work: 2_000_000_000,
    }
}
fn response_options(omega: f64) -> PoleResponseOptions {
    PoleResponseOptions {
        spectral_scale: omega,
        pole_step_tolerance: 1e-10,
        pole_residual_tolerance: 1e-10,
        linear_residual_tolerance: 1e-11,
        max_dimension: 1024,
        max_vectors: 2,
        max_dense_bytes: 128 * 1024 * 1024,
        max_work: 2_000_000_000,
    }
}

#[test]
fn original_cylinder_multiport_residues_converge_to_standing_wave_participation() {
    let air = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap();
    let (radius, length) = (0.01, 0.5);
    let fundamental = core::f64::consts::PI * air.sound_speed / (2.0 * length);
    let residue_scale = air.density * air.sound_speed * air.sound_speed
        / (core::f64::consts::PI * radius * radius * length);
    for harmonic in [1.0, 3.0] {
        let omega = harmonic * fundamental;
        // Z_ij(s) has residue rho*c²/(S*L)*cos(k*x_i)*cos(k*x_j).
        // At the second resonance the side port at L/3 is a pressure node.
        let side = det::cos(harmonic * core::f64::consts::PI / 6.0);
        let exact = [1.0, side, side, side * side];
        let mut errors = Vec::new();
        for elements in [3, 6, 12] {
            let duct = DgDuct::new(
                &Duct {
                    segments: vec![Segment::Cylinder { radius, length }],
                },
                &air,
                DgMesh { elements, order: 2 },
                DgBoundary::Port(PortId(7)),
                DgBoundary::IdealOpen,
            )
            .unwrap()
            .with_side_ports(&[DgSidePort {
                id: PortId(11),
                downstream_element: elements / 3,
            }])
            .unwrap();
            let seed = duct
                .harmonic_response(
                    0.98 * omega,
                    &[(PortId(7), C64::ONE)],
                    WORK,
                    1e-10,
                    &mut || false,
                )
                .unwrap()
                .energy_coordinates;
            let prepared = duct
                .prepare_resonance_poles(omega, WORK as u64, || false)
                .unwrap();
            let pole = prepared
                .resonance_pole(C64::new(0.0, -omega), &seed, &options(omega), || false)
                .unwrap();
            let result = prepared
                .port_residue(&pole, &response_options(omega), || false)
                .unwrap();
            assert_eq!(result.ports, vec![PortId(7), PortId(11)]);
            let error = result
                .pressure_per_volume_residue_pa_m3
                .iter()
                .zip(exact)
                .map(|(&got, expected)| {
                    (got.scale(1.0 / residue_scale) - C64::from_re(expected)).abs()
                })
                .fold(0.0_f64, f64::max);
            errors.push(error);
            assert!(result.original_relative_residual < 1e-10);
            assert!(result.forcing_relative_residuals.iter().all(|r| *r < 1e-11));
            assert_eq!(
                prepared
                    .port_residue(&pole, &response_options(omega), || false)
                    .unwrap(),
                result
            );
        }
        assert!(errors[1] < errors[0] && errors[2] < errors[1], "{errors:?}");
        assert!(errors[2] < 3e-3, "{errors:?}");
    }
}

// Independent telegraph law for the retained original wall branches. Spatial
// resonances satisfy p(s) q(s)+(c*k)^2=0. Differentiating the original
// Z(s)=rho*c/S*sqrt(q/p)*tanh(L*sqrt(p*q)/c) gives the complex port residue.
fn telegraph_coefficients(model: &fs_phs::ZkDiffusive, s: C64) -> (C64, C64, C64) {
    let mut p = s;
    let mut q = s + C64::from_re(model.series_drag_per_s());
    let (mut dp, mut dq) = (C64::ONE, C64::ONE);
    for (terms, value, derivative) in [
        (model.shunt_terms(), &mut p, &mut dp),
        (model.series_terms(), &mut q, &mut dq),
    ] {
        for term in terms {
            let denominator = s + C64::from_re(term.rate_per_s);
            *value = *value + (s / denominator).scale(term.gain_per_s);
            *derivative = *derivative
                + C64::from_re(term.gain_per_s * term.rate_per_s) / (denominator * denominator);
        }
    }
    (p, q, dp * q + p * dq)
}

#[test]
fn original_viscothermal_residue_converges_with_all_material_memories() {
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap();
    let (radius, length) = (0.001, 0.3);
    let omega = core::f64::consts::PI * gas.sound_speed / (2.0 * length);
    let residue_scale = gas.density * gas.sound_speed * gas.sound_speed
        / (core::f64::consts::PI * radius * radius * length);
    let mut errors = Vec::new();
    for elements in [1, 2, 4] {
        let duct = DgDuct::new(
            &Duct {
                segments: vec![Segment::Cylinder { radius, length }],
            },
            &gas,
            DgMesh { elements, order: 2 },
            DgBoundary::Port(PortId(7)),
            DgBoundary::IdealOpen,
        )
        .unwrap()
        .with_viscothermal(4, [80.0, 1200.0], 0.2)
        .unwrap();
        let model = duct.viscothermal_model().unwrap();
        let mut expected_pole = C64::new(-10.0, -omega);
        for _ in 0..40 {
            let (p, q, derivative) = telegraph_coefficients(model, expected_pole);
            let correction = (p * q + C64::from_re(omega * omega)) / derivative;
            expected_pole = expected_pole - correction;
            if correction.abs() < 1e-13 * omega {
                break;
            }
        }
        let (p, q, derivative) = telegraph_coefficients(model, expected_pole);
        assert!((p * q + C64::from_re(omega * omega)).abs() < 1e-10 * omega * omega);
        let expected_residue = (q / derivative).scale(2.0 * residue_scale);
        assert!(expected_residue.im.abs() > 1e-4 * residue_scale);
        let seed = duct
            .harmonic_response(
                0.98 * expected_pole.im.abs(),
                &[(PortId(7), C64::ONE)],
                WORK,
                1e-10,
                &mut || false,
            )
            .unwrap()
            .energy_coordinates;
        let prepared = duct
            .prepare_resonance_poles(omega, WORK as u64, || false)
            .unwrap();
        let pole = prepared
            .resonance_pole(expected_pole, &seed, &options(omega), || false)
            .unwrap();
        assert!(pole.energy_coordinates.len() > 2 * duct.field_dofs());
        assert!(pole.auxiliary_stored_energy_j > 0.0 && pole.wall_dissipation_w > 0.0);
        let residue = prepared
            .port_residue(&pole, &response_options(omega), || false)
            .unwrap();
        errors.push(
            (residue.pressure_per_volume_residue_pa_m3[0] - expected_residue).abs() / residue_scale,
        );
    }
    assert!(errors[1] < errors[0] && errors[2] < errors[1], "{errors:?}");
    assert!(errors[2] < 3e-3, "{errors:?}");
}

#[test]
fn original_mode_limits_and_cancellation_are_rechecked_before_publication() {
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap();
    let omega = core::f64::consts::PI * gas.sound_speed;
    let duct = DgDuct::new(
        &Duct {
            segments: vec![Segment::Cylinder {
                radius: 0.01,
                length: 0.5,
            }],
        },
        &gas,
        DgMesh {
            elements: 1,
            order: 2,
        },
        DgBoundary::Port(PortId(7)),
        DgBoundary::IdealOpen,
    )
    .unwrap();
    let seed = duct
        .harmonic_response(
            0.98 * omega,
            &[(PortId(7), C64::ONE)],
            WORK,
            1e-10,
            &mut || false,
        )
        .unwrap()
        .energy_coordinates;
    let prepared = duct
        .prepare_resonance_poles(omega, WORK as u64, || false)
        .unwrap();
    let pole = prepared
        .resonance_pole(C64::new(0.0, -omega), &seed, &options(omega), || false)
        .unwrap();
    let limits = response_options(omega);
    let before = prepared.port_residue(&pole, &limits, || false).unwrap();
    for reference in [omega, omega * 1e100] {
        let rescaled = duct
            .prepare_resonance_poles(reference, WORK as u64, || false)
            .unwrap();
        for gain in [C64::new(1e200, -3e200), C64::new(-1e-200, 2e-200)] {
            let mut changed = pole.clone();
            for value in &mut changed.energy_coordinates {
                *value = *value * gain;
            }
            let actual = rescaled.port_residue(&changed, &limits, || false).unwrap();
            for (&got, &expected) in actual
                .pressure_per_volume_residue_pa_m3
                .iter()
                .zip(&before.pressure_per_volume_residue_pa_m3)
            {
                assert!((got - expected).abs() < 1e-10 * expected.abs());
            }
        }
    }
    let mut forged = pole.clone();
    forged.laplace_pole_per_s = forged.laplace_pole_per_s.scale(0.9);
    forged.relative_residual = 0.0;
    assert!(prepared.port_residue(&forged, &limits, || false).is_err());
    forged = pole.clone();
    forged.energy_coordinates.pop();
    assert!(prepared.port_residue(&forged, &limits, || false).is_err());
    for limited in [
        PoleResponseOptions {
            max_dimension: 1,
            ..limits
        },
        PoleResponseOptions {
            max_vectors: 0,
            ..limits
        },
        PoleResponseOptions {
            max_dense_bytes: 1,
            ..limits
        },
        PoleResponseOptions {
            max_work: 1,
            ..limits
        },
    ] {
        assert!(prepared.port_residue(&pole, &limited, || false).is_err());
    }
    for boundary in [0, 10, 25] {
        let mut calls = 0;
        assert!(
            prepared
                .port_residue(&pole, &limits, || {
                    calls += 1;
                    calls > boundary
                })
                .is_err()
        );
    }
    assert_eq!(
        prepared.port_residue(&pole, &limits, || false).unwrap(),
        before
    );
}
