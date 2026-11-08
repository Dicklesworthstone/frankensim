//! Physical modal sensitivity when active work cannot measure the response.

use super::*;

fn conservative_model(
    structural: &[f64],
    acoustic: &[f64],
    coupling: Vec<f64>,
) -> VibroacousticModel {
    VibroacousticModel::try_new(
        &StructuralModes {
            omegas: structural.to_vec(),
            shapes: vec![vec![1.0]; structural.len()],
            loss_factor: 0.0,
        },
        &CavityModes {
            omegas: acoustic.to_vec(),
            lambdas: vec![1.0; acoustic.len()],
            interface: vec![vec![1.0]; acoustic.len()],
            loss_factor: 0.0,
            rho0: 1.0,
            c0: 1.0,
        },
        coupling,
        None,
    )
    .unwrap()
}

#[test]
fn conservative_structural_resonance_cannot_converge_from_zero_work() {
    let model = conservative_model(&[1.0, 2.0], &[0.0], vec![0.0; 2]);
    let force = [C64::ONE; 2];
    let full = model.frf(1.99, &force).unwrap();
    let half = model.frf_truncated(1.99, &force, 1, 1).unwrap();
    assert_eq!(full.power.input, 0.0);
    assert_eq!(half.power.input, 0.0);
    assert!(full.b[1].abs() > 70.0 * full.b[0].abs());
    let err = model.frf_with_convergence(1.99, &force, 0.05).unwrap_err();
    let VibroError::TruncationNotConverged { relative_delta, .. } = err else {
        panic!("wrong refusal: {err:?}");
    };
    assert!(relative_delta > 0.99);
}

#[test]
fn conservative_acoustic_resonance_survives_small_source_amplitudes() {
    let model = conservative_model(&[3.0], &[0.0, 2.0], vec![0.0; 2]);
    let forces = [C64::ZERO];
    let mut deltas = Vec::new();
    // The smaller drive has representable pressure but underflowed pressure
    // squares. The relative physical response must be independent of level.
    for level in [1.0, 1e-180] {
        let flows = [C64::from_re(0.1 * level), C64::from_re(level)];
        let checked = model
            .frf_with_sources_and_convergence(1.99, &forces, &flows, 2.0)
            .unwrap();
        assert_eq!(checked.response.power.input, 0.0);
        assert_eq!(checked.input_power_delta, 0.0);
        assert!(checked.energy_norm_delta > 0.999);
        assert_eq!(checked.truncation_delta, checked.energy_norm_delta);
        deltas.push(checked.energy_norm_delta);
        assert!(matches!(
            model.frf_with_sources_and_convergence(1.99, &forces, &flows, 0.05),
            Err(VibroError::TruncationNotConverged { .. }),
        ));
    }
    assert!((deltas[0] - deltas[1]).abs() < 2e-14);
}

#[test]
fn coherent_sources_cannot_hide_an_omitted_mode_with_cancelling_work() {
    let model = conservative_model(&[1.0, 1.01], &[0.0, 2.0], vec![1.0, 0.0, 0.0, 0.0]);
    let forces = [C64::ONE; 2];
    let flows = [C64::ONE, C64::ZERO];
    let full = model.frf_with_sources(1.0, &forces, &flows).unwrap();
    let half = model
        .frf_truncated_with_sources(1.0, &forces, &flows, 1, 1)
        .unwrap();
    // Exact first-block solution: b_0 = 1+i, p_0 = -1. The mechanical
    // source supplies 1/2 W and the acoustic source absorbs the same work;
    // the omitted near-resonant second structure mode stores reactive energy.
    for response in [&full, &half] {
        assert_eq!(response.power.structural_input, 0.5);
        assert_eq!(response.power.acoustic_input, -0.5);
        assert_eq!(response.power.input, 0.0);
    }
    let err = model
        .frf_with_sources_and_convergence(1.0, &forces, &flows, 0.05)
        .unwrap_err();
    let VibroError::TruncationNotConverged { relative_delta, .. } = err else {
        panic!("wrong refusal: {err:?}");
    };
    assert!(relative_delta > 0.99);
}

fn scaled_cavity_model(scales: [f64; 2]) -> VibroacousticModel {
    VibroacousticModel::try_new(
        &StructuralModes {
            omegas: vec![2.0, 3.5],
            shapes: vec![vec![1.0]; 2],
            loss_factor: 0.03,
        },
        &CavityModes {
            // Unit-length, unit-area cavity: psi_0 = s_0,
            // psi_1 = s_1 cos(pi x), so Lambda = [s_0^2, s_1^2/2].
            omegas: vec![0.0, core::f64::consts::PI],
            lambdas: vec![scales[0] * scales[0], 0.5 * scales[1] * scales[1]],
            interface: scales.iter().map(|&s| vec![s]).collect(),
            loss_factor: 0.04,
            rho0: 1.0,
            c0: 1.0,
        },
        vec![
            0.3 * scales[0],
            0.25 * scales[1],
            -0.2 * scales[0],
            0.1 * scales[1],
        ],
        None,
    )
    .unwrap()
}

fn physical_energy(omega: f64, b: &[C64], a: &[C64], scales: [f64; 2]) -> f64 {
    let structural = 0.25
        * b.iter()
            .zip([2.0_f64, 3.5])
            .map(|(&b, natural)| (omega * omega + natural * natural) * b.norm_sq())
            .sum::<f64>();
    // Independently reconstruct physical pressure and gradient-derived
    // velocity, then integrate their cycle-averaged energy in the unit cavity.
    let mut acoustic = 0.0;
    for i in 0..64 {
        let phase = core::f64::consts::PI * (i as f64 + 0.5) / 64.0;
        let pressure =
            a[0].scale(scales[0]) + a[1].scale(scales[1] * fs_math::det::cos(phase));
        let gradient =
            a[1].scale(-scales[1] * core::f64::consts::PI * fs_math::det::sin(phase));
        let velocity = C64::new(0.0, -1.0 / omega) * gradient;
        acoustic += 0.25 * (pressure.norm_sq() + velocity.norm_sq()) / 64.0;
    }
    structural + acoustic
}

#[test]
fn complex_response_norm_matches_physical_energy_and_cavity_basis_scaling() {
    let omega = 2.6;
    let forces = [C64::new(1.0, 0.4), C64::new(-0.3, 0.2)];
    let original_flows = [C64::new(0.6, -0.2), C64::new(0.2, 0.3)];
    let mut reference: Option<f64> = None;
    for scales in [[1.0, 1.0], [-7.0, 1e-4], [1e4, -1e4]] {
        let model = scaled_cavity_model(scales);
        let flows = [
            original_flows[0].scale(scales[0]),
            original_flows[1].scale(scales[1]),
        ];
        let checked = model
            .frf_with_sources_and_convergence(omega, &forces, &flows, 10.0)
            .unwrap();
        let full = &checked.response;
        let half = model
            .frf_truncated_with_sources(omega, &forces, &flows, 1, 1)
            .unwrap();
        let db = [full.b[0] - half.b[0], full.b[1]];
        let da = [full.a[0] - half.a[0], full.a[1]];
        let expected = (physical_energy(omega, &db, &da, scales)
            / physical_energy(omega, &full.b, &full.a, scales))
        .sqrt();
        assert!((checked.energy_norm_delta - expected).abs() < 2e-13);
        assert_eq!(
            checked.truncation_delta,
            checked.input_power_delta.max(checked.energy_norm_delta),
        );
        if let Some(delta) = reference {
            assert!((checked.energy_norm_delta - delta).abs() < 2e-12);
        } else {
            reference = Some(checked.energy_norm_delta);
        }
        assert!(matches!(
            model.frf_with_sources_and_convergence(
                omega,
                &forces,
                &flows,
                0.5 * checked.energy_norm_delta,
            ),
            Err(VibroError::TruncationNotConverged { .. }),
        ));
    }
}

#[test]
fn zero_fields_have_zero_relative_error() {
    let model = scaled_cavity_model([1.0, 1.0]);
    let checked = model
        .frf_with_convergence(2.6, &[C64::ZERO; 2], 1e-12)
        .unwrap();
    assert_eq!(checked.input_power_delta, 0.0);
    assert_eq!(checked.energy_norm_delta, 0.0);
    assert_eq!(checked.truncation_delta, 0.0);
}
