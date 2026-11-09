//! G1/G3: actual material stiffness, coupled sound and separately positive heat.
use super::*;
use crate::vibroacoustic::{
    CavityModes, ModalRadiationLoading, StructuralModes, VibroacousticModel, radiation_sample,
};

fn arms() -> Vec<RelaxationBranch> {
    vec![
        RelaxationBranch {
            projection: vec![1.0, -0.5],
            stiffness: 4.0,
            relaxation_time_s: 0.4,
        },
        RelaxationBranch {
            projection: vec![0.3, 1.0],
            stiffness: 1.0,
            relaxation_time_s: 0.9,
        },
    ]
}

fn material(branches: &[RelaxationBranch], band: (f64, f64)) -> ModalStructuralRelaxation {
    ModalStructuralRelaxation::try_new(
        2,
        branches,
        band,
        branches.len(),
        &CancelGate::new_clock_free(),
    )
    .unwrap()
}

fn model(coupling: Vec<f64>, radiation: Option<ModalRadiationLoading>) -> VibroacousticModel {
    VibroacousticModel::try_new(
        &StructuralModes {
            omegas: vec![3.0, 5.0],
            shapes: vec![vec![1.0], vec![1.0]],
            loss_factor: 0.0,
        },
        &CavityModes {
            omegas: vec![4.0],
            lambdas: vec![0.7],
            interface: vec![vec![1.0]],
            loss_factor: 0.04,
            rho0: 1.2,
            c0: 3.0,
        },
        coupling,
        radiation,
    )
    .unwrap()
}

fn near(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 3e-12 * a.abs().max(b.abs()).max(1e-25),
        "{a:e} != {b:e}"
    );
}

fn near_complex(a: C64, b: C64) {
    assert!(
        (a - b).abs() <= 3e-12 * a.abs().max(b.abs()).max(1e-25),
        "{a:?} != {b:?}"
    );
}

// Independent scalar spring/dashpot equation: e = -i*w*tau/(1-i*w*tau) * c.q.
// This oracle does not call the production spectral-weight owner.
fn spring_fraction(omega: f64, tau: f64) -> C64 {
    let s = C64::new(0.0, -omega * tau);
    s / (C64::ONE + s)
}

#[test]
fn hereditary_material_changes_cavity_and_radiation_with_independent_heat_closure() {
    let branches = arms();
    let coupling = [0.03, -0.06];
    let force = [C64::new(0.8, 0.3), C64::new(-0.2, 0.1)];
    let flow = C64::new(0.007, -0.002);
    let z = [
        C64::new(0.2, -0.1),
        C64::new(0.04, 0.02),
        C64::new(0.04, 0.02),
        C64::new(0.3, -0.15),
    ];
    for omega in [0.7, 1.9, 3.4, 5.8] {
        let original = model(
            coupling.to_vec(),
            Some(radiation_sample(omega, 2, z.to_vec())),
        );
        let dry = original.frf_with_sources(omega, &force, &[flow]).unwrap();
        let loaded = original
            .with_structural_relaxation(material(&branches, (0.1, 10.0)))
            .unwrap();
        let response = loaded.frf_with_sources(omega, &force, &[flow]).unwrap();
        let da = C64::new(16.0 - omega * omega, -0.04 * 16.0);
        let gamma = 1.2 * 3.0 * 3.0 / 0.7 * omega * omega;
        let source = C64::new(0.0, -omega * (1.2 * 3.0 * 3.0 / 0.7)) * flow;
        let mut d = [C64::ZERO; 4];
        for r in 0..2 {
            for s in 0..2 {
                d[2 * r + s] = C64::from_re(if r == s {
                    [9.0, 25.0][r] - omega * omega
                } else {
                    0.0
                }) - C64::new(0.0, omega) * z[2 * r + s];
                for arm in &branches {
                    d[2 * r + s] = d[2 * r + s]
                        + spring_fraction(omega, arm.relaxation_time_s)
                            .scale(arm.stiffness * arm.projection[r] * arm.projection[s]);
                }
                d[2 * r + s] = d[2 * r + s] - C64::from_re(gamma * coupling[r] * coupling[s]) / da;
            }
        }
        let rhs: [C64; 2] = core::array::from_fn(|r| force[r] + source.scale(coupling[r]) / da);
        let determinant = d[0] * d[3] - d[1] * d[2];
        let expected = [
            (d[3] * rhs[0] - d[1] * rhs[1]) / determinant,
            (d[0] * rhs[1] - d[2] * rhs[0]) / determinant,
        ];
        let pressure = (source
            + (expected[0].scale(coupling[0]) + expected[1].scale(coupling[1])).scale(gamma))
            / da;
        for r in 0..2 {
            near_complex(response.b[r], expected[r]);
        }
        near_complex(response.a[0], pressure);
        assert!((response.a[0] - dry.a[0]).abs() > 1e-6 * response.a[0].abs());
        let mut heat = 0.0;
        let mut energy = 0.0;
        for arm in &branches {
            let strain =
                expected[0].scale(arm.projection[0]) + expected[1].scale(arm.projection[1]);
            let elastic = spring_fraction(omega, arm.relaxation_time_s) * strain;
            heat += 0.5 * arm.stiffness / arm.relaxation_time_s * elastic.norm_sq();
            energy += 0.25 * arm.stiffness * elastic.norm_sq();
        }
        near(response.power.structural_relaxation, heat);
        near(response.power.structural, heat);
        near(response.structural_relaxation_stored_energy_j, energy);
        assert!(heat > 0.0 && response.power.radiated > 0.0 && response.power.cavity > 0.0);
        assert!(response.power.balance_residual.abs() < 2e-12);
        assert!(response.power.interface_residual.abs() < 2e-12);

        // The same full source rows supply the retained principal block.
        let truncated = loaded
            .frf_truncated_with_sources(omega, &force, &[flow], 1, 1)
            .unwrap();
        let x = rhs[0] / d[0];
        near_complex(truncated.b[0], x);
        near_complex(truncated.a[0], (source + x.scale(gamma * coupling[0])) / da);
        assert!(truncated.power.balance_residual.abs() < 2e-12);
    }
}

#[test]
fn original_material_channels_keep_relaxed_and_unrelaxed_stiffness_and_loss_nullspace() {
    let branches = [RelaxationBranch {
        projection: vec![1.0, -1.0],
        stiffness: 2.0,
        relaxation_time_s: 0.5,
    }];
    let material = material(&branches, (0.0, 1e12));
    for omega in [0.0, 2.0, 1e10] {
        let sample = material.at_frequency(omega).unwrap();
        let mut matrix = [C64::ZERO; 4];
        sample.add_stiffness(&mut matrix, 2, 2).unwrap();
        let expected = spring_fraction(omega, 0.5).scale(2.0);
        for (actual, sign) in matrix.into_iter().zip([1.0, -1.0, -1.0, 1.0]) {
            near_complex(actual, expected.scale(sign));
        }
        let (stored, heat) = sample.energy_and_heat(&[C64::new(0.7, -0.3); 2]).unwrap();
        assert_eq!((stored, heat), (0.0, 0.0));
    }
    let mut loaded = model(vec![0.0; 2], None);
    loaded.omegas_s = vec![3.0; 2];
    let loaded = loaded.with_structural_relaxation(material).unwrap();
    let common = loaded.frf(2.0, &[C64::ONE; 2]).unwrap();
    for value in common.b {
        near_complex(value, C64::from_re(0.2));
    }
    assert!(common.power.structural_relaxation.abs() < 1e-28);
    let differential = loaded.frf(2.0, &[C64::ONE, -C64::ONE]).unwrap();
    let x = (C64::from_re(5.0) + spring_fraction(2.0, 0.5).scale(4.0)).recip();
    near_complex(differential.b[0], x);
    near_complex(differential.b[1], -x);
    assert!(differential.power.structural_relaxation > 0.0);
    assert!(differential.power.balance_residual.abs() < 2e-12);
}

#[test]
fn original_source_band_caps_cancellation_and_scalar_loss_are_respected() {
    let branches = arms();
    let gate = CancelGate::new_clock_free();
    assert!(ModalStructuralRelaxation::try_new(2, &branches, (1.0, 4.0), 1, &gate).is_err());
    assert!(ModalStructuralRelaxation::try_new(1, &branches, (1.0, 4.0), 2, &gate).is_err());
    assert!(ModalStructuralRelaxation::try_new(2, &branches, (4.0, 1.0), 2, &gate).is_err());
    let mut invalid = branches.clone();
    invalid[0].relaxation_time_s = 0.0;
    assert!(ModalStructuralRelaxation::try_new(2, &invalid, (1.0, 4.0), 2, &gate).is_err());
    let mut scalar_loss = model(vec![0.0; 2], None);
    scalar_loss.eta_s = 0.01;
    assert!(
        scalar_loss
            .with_structural_relaxation(material(&branches, (1.0, 4.0)))
            .is_err()
    );
    let coupled = model(vec![0.03, -0.06], None)
        .with_structural_relaxation(material(&branches, (1.0, 4.0)))
        .unwrap();
    for omega in [0.9, 4.1] {
        assert!(coupled.frf(omega, &[C64::ONE; 2]).is_err());
        assert!(coupled.frf_truncated(omega, &[C64::ONE; 2], 1, 1).is_err());
        assert!(
            coupled
                .frf_with_convergence(omega, &[C64::ONE; 2], 1.0)
                .is_err()
        );
    }
    for omega in [1.0, 4.0] {
        assert!(coupled.frf(omega, &[C64::ONE; 2]).is_ok());
    }
    assert!(
        coupled
            .with_structural_relaxation(material(&branches, (1.0, 4.0)))
            .is_err()
    );
    gate.request();
    assert!(matches!(
        ModalStructuralRelaxation::try_new(2, &branches, (1.0, 4.0), 2, &gate),
        Err(VibroError::Cancelled),
    ));
}

#[test]
fn empty_original_material_preserves_the_existing_coupled_response() {
    let original = model(vec![0.03, -0.06], None);
    let force = [C64::new(0.5, 0.1), C64::new(-0.3, 0.2)];
    let before = original.frf(2.3, &force).unwrap();
    let after = original
        .with_structural_relaxation(material(&[], (0.0, 10.0)))
        .unwrap()
        .frf(2.3, &force)
        .unwrap();
    assert_eq!(before.b, after.b);
    assert_eq!(before.a, after.a);
    assert_eq!(before.power.input, after.power.input);
    assert_eq!(after.power.structural_relaxation, 0.0);
    assert_eq!(after.structural_relaxation_stored_energy_j, 0.0);
}
