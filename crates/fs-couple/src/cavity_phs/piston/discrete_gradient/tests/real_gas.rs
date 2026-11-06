//! Original PR calorics and entropy drive the same wall, not an ideal proxy.

use super::*;
use crate::cavity_phs::piston::tests::real_gas::{real_controls, real_parcel};
use fs_material::fluid_acoustics::FluidAcoustics;
use fs_material::fluid_transport::FluidTransport;
use fs_material::thermochemical_gas::nonideal::energy::parcel::IsochoricRealGasParcel;

fn input(heat_j: f64) -> PistonInput {
    PistonInput {
        duration_s: 0.00005,
        heat_j,
        exterior_pressure_pa: 4.5e6,
        external_force_n: -100.0,
    }
}

fn entropy(initial: &IsochoricRealGasParcel, endpoint: &IsochoricRealGasParcel) -> f64 {
    initial.mass_kg()
        * initial
            .model()
            .specific_entropy_change_between_densities(
                initial.density_kg_m3(),
                initial.state().temperature_k(),
                endpoint.density_kg_m3(),
                endpoint.state().temperature_k(),
            )
            .unwrap()
}

#[test]
fn dense_gas_motion_preserves_original_entropy_and_total_trajectory_energy() {
    for velocity in [-0.4, 0.4] {
        let mut piston = GasPiston::new(real_parcel(300.0, 5e6), 2.0, 0.01, velocity).unwrap();
        let initial = piston.clone();
        let mut external = 0.0;
        for _ in 0..24 {
            let report = piston
                .step_discrete_gradient(input(0.0), real_controls(), || false)
                .unwrap();
            external += report.external_force_work_j - report.exterior_pressure_work_j;
            assert!(report.momentum_residual_n_s.abs() <= real_controls().momentum_tolerance_n_s);
            assert!(report.total_energy_residual_j.abs() <= real_controls().energy_tolerance_j);
            let mean = &report.thermodynamic_mean;
            let transport = mean.transport_state().unwrap();
            assert_eq!(transport.density_kg_m3(), mean.density_kg_m3());
            assert_eq!(
                transport.specific_heat_j_kg_k(),
                mean.specific_heat_cp_j_kg_k()
            );
            assert_eq!(
                mean.acoustic_state().unwrap().sound_speed_m_s(),
                mean.sound_speed_m_s()
            );
            assert!(
                (report.gas.work_by_gas_j
                    - report.gas.work_pressure_pa
                        * (report.gas.final_volume_m3 - report.gas.initial_volume_m3))
                    .abs()
                    < 1e-10
            );
        }
        let energy = initial.gas().mass_kg()
            * initial
                .gas()
                .model()
                .specific_internal_energy_change_between_densities(
                    initial.gas().density_kg_m3(),
                    initial.gas().state().temperature_k(),
                    piston.gas().density_kg_m3(),
                    piston.gas().state().temperature_k(),
                )
                .unwrap();
        let kinetic = 0.5
            * piston.mass_kg()
            * (piston.velocity_m_s().powi(2) - initial.velocity_m_s().powi(2));
        assert!((energy + kinetic - external).abs() < 24.0 * real_controls().energy_tolerance_j);
        let bound =
            8.0 * 24.0 * initial.gas().mass_kg() * real_controls().gas.energy_tolerance_j_kg
                / real_controls().gas.minimum_temperature_k;
        assert!(entropy(initial.gas(), piston.gas()).abs() <= bound);
        assert_eq!(initial.gas().mass_kg(), piston.gas().mass_kg());
        assert_eq!(initial.gas().model(), piston.gas().model());
        assert_eq!(piston.gas().point().axes().get("specimen-lot"), Some(&42.0));
        assert!(piston.audit_last_window().is_ok());
    }
}

#[test]
fn fixed_dense_gas_heat_uses_distinct_original_sampling_state_and_work_effort() {
    let initial = real_parcel(300.0, 5e6);
    let mut reference = initial.clone();
    reference
        .add_heat(3000.0, real_controls().gas, || false)
        .unwrap();
    let mut piston = GasPiston::new(initial.clone(), 2.0, 0.01, 0.0)
        .unwrap()
        .with_fixed_support()
        .unwrap();
    let config = PistonSolve {
        volume_bracket_m3: [initial.volume_m3(); 2],
        ..real_controls()
    };
    let report = piston
        .step_discrete_gradient(input(3000.0), config, || false)
        .unwrap();
    let mean = &report.thermodynamic_mean;
    // Independent PR/NASA coordinate-gradient evaluation of this synthetic
    // source, rather than a pressure sampled from the ideal reference.
    assert!((piston.gas().state().temperature_k() - 364.77082151484).abs() < 1e-7);
    assert!((mean.temperature_k() - 331.42100073968).abs() < 1e-7);
    assert!((report.gas.work_pressure_pa - 5645146.27826).abs() < 0.01);
    assert!((mean.pressure_pa() - 5643512.96843).abs() < 0.01);
    assert!(
        (piston.gas().state().temperature_k() - reference.state().temperature_k()).abs() < 1e-8
    );
    assert!((entropy(&initial, piston.gas()) - 3000.0 / mean.temperature_k()).abs() < 1e-9);
    assert_eq!(report.gas.work_by_gas_j, 0.0);
    assert_eq!(report.support_work_j, 0.0);
    assert_eq!(piston.gas().density_kg_m3(), initial.density_kg_m3());
    // A discrete conjugate effort is generally not the PR pressure at T_bar.
    // Keep the original material's density, EOS and transport coherent.
    assert!((report.gas.work_pressure_pa - mean.pressure_pa()).abs() > 1.0);
    assert_eq!(mean.density_kg_m3(), initial.density_kg_m3());
    let context = initial
        .point()
        .clone()
        .with("temperature", mean.temperature_k())
        .unwrap();
    let resolved = initial
        .model()
        .evaluate_density_temperature(&context, mean.density_kg_m3())
        .unwrap();
    assert_eq!(*mean, resolved);
    let expected_support = -input(3000.0).duration_s
        * (piston.area_m2() * (report.gas.work_pressure_pa - input(3000.0).exterior_pressure_pa)
            + input(3000.0).external_force_n);
    assert!((report.support_impulse_n_s - expected_support).abs() < 1e-12);
    assert!(piston.audit_last_window().is_ok());
}

#[test]
fn dense_gas_impulse_tangent_and_heat_entropy_use_the_same_source_gradient() {
    let piston = GasPiston::new(real_parcel(300.0, 5e6), 2.0, 0.01, 0.4).unwrap();
    let volume = 0.001001;
    let delta = 1e-8;
    let drive = input(100.0);
    let center = piston
        .discrete_gradient_trial(volume, drive, real_controls(), &mut || false)
        .unwrap();
    let plus = piston
        .discrete_gradient_trial(volume + delta, drive, real_controls(), &mut || false)
        .unwrap();
    let minus = piston
        .discrete_gradient_trial(volume - delta, drive, real_controls(), &mut || false)
        .unwrap();
    let derivative = (plus.residual - minus.residual) / (2.0 * delta);
    assert!((center.tangent / derivative - 1.0).abs() < 2e-6);
    for heat in [-50.0, 100.0] {
        let mut next = piston.clone();
        let report = next
            .step_discrete_gradient(input(heat), real_controls(), || false)
            .unwrap();
        assert!(
            (entropy(piston.gas(), next.gas()) - heat / report.thermodynamic_mean.temperature_k())
                .abs()
                < 1e-9
        );
        assert!(report.total_energy_residual_j.abs() <= real_controls().energy_tolerance_j);
    }
}

#[test]
fn nonideal_source_refusal_and_late_cancellation_preserve_gas_motion_and_ports() {
    let mut piston = GasPiston::new(real_parcel(300.0, 5e6), 2.0, 0.01, 0.4).unwrap();
    piston
        .step_discrete_gradient(input(0.0), real_controls(), || false)
        .unwrap();
    let before = piston.clone();
    let transfers = piston.last_window_transfers().map(<[_]>::to_vec);
    let duration = piston.last_window_duration_s();
    assert!(
        piston
            .step_discrete_gradient(input(1e9), real_controls(), || false)
            .is_err()
    );
    let mut probe = before.clone();
    let mut polls = 0;
    probe
        .step_discrete_gradient(input(100.0), real_controls(), || {
            polls += 1;
            false
        })
        .unwrap();
    let mut calls = 0;
    assert!(matches!(
        piston.step_discrete_gradient(input(100.0), real_controls(), || {
            calls += 1;
            calls == polls
        }),
        Err(PistonError::Cancelled)
    ));
    assert_eq!(piston, before);
    assert_eq!(piston.last_window_transfers().map(<[_]>::to_vec), transfers);
    assert_eq!(piston.last_window_duration_s(), duration);
}
