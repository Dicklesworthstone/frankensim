//! G0/G1/G3/G4: local vortex memory, represented storage and native wave ports.
//! Coefficients below are authored numerical fixtures, not helium calibration.

use std::cell::Cell;

use fs_duct::superfluid::{
    BodyForce, BoundaryVelocities, ChannelError, HeatStepError, LinearizationLimits,
    MutualFriction, TwoFluidChannel,
};
use fs_material::superfluid::vinen::{VinenParameters, VinenStepControls, VinenTangle};
use fs_material::superfluid::{NormalViscousSpectrum, TwoFluidProperties, TwoFluidState};

const FRICTION: MutualFriction = MutualFriction::EvolvingVinenTangle;
const LINE_ENERGY: f64 = 1e-15;
const DECAY: f64 = 1e-7;

fn material() -> TwoFluidState {
    TwoFluidState::try_new(TwoFluidProperties {
        temperature_k: 2.0,
        density_kg_m3: 4.0,
        superfluid_fraction: 0.75,
        entropy_j_kg_k: 4.0,
        specific_heat_cp_j_kg_k: 16.0,
        adiabatic_sound_speed_m_s: 10.0,
        thermal_expansion_per_k: -0.2,
        normal_viscosity_pa_s: 1e-4,
    })
    .unwrap()
}

fn model(source: TwoFluidState) -> VinenTangle {
    VinenTangle::try_new(
        source,
        VinenParameters {
            circulation_quantum_m2_s: 1e-7,
            dissipative_b: 1.2,
            production_coefficient: 0.8,
            decay_coefficient_m2_s: DECAY,
            line_energy_j_m: LINE_ENERGY,
        },
        VinenStepControls {
            maximum_fractional_step: 0.04,
            maximum_substeps: 512,
        },
    )
    .unwrap()
}

fn channel(n: usize, profile: bool) -> TwoFluidChannel {
    let source = material();
    let limits = LinearizationLimits {
        max_pressure_perturbation_pa: 1.0,
        max_temperature_perturbation_k: 0.02,
        max_component_speed_m_s: 0.1,
    };
    if profile {
        TwoFluidChannel::try_new_with_normal_profile(
            source,
            1.0,
            0.01,
            n,
            limits,
            NormalViscousSpectrum::circular_capillary(source, 0.01, 3).unwrap(),
        )
    } else {
        TwoFluidChannel::try_new(source, 1.0, 0.01, n, limits)
    }
    .unwrap()
}

fn seed(duct: &mut TwoFluidChannel, densities: &[f64]) {
    duct.configure_vinen_tangle(model(*duct.material()), densities, || false)
        .unwrap();
}

fn pulse(duct: &mut TwoFluidChannel) {
    let n = duct.cells();
    let shape: Vec<_> = (0..=n)
        .map(|i| (std::f64::consts::PI * i as f64 / n as f64).sin())
        .collect();
    duct.initialize(
        &vec![0.0; n],
        &vec![0.0; n],
        &shape.iter().map(|x| 0.003 * x).collect::<Vec<_>>(),
        &shape.iter().map(|x| -0.001 * x).collect::<Vec<_>>(),
        || false,
    )
    .unwrap();
}

fn close(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() <= 3e-11 * scale,
        "actual={actual:e}, expected={expected:e}, scale={scale:e}"
    );
}

fn snapshot(duct: &TwoFluidChannel) -> Vec<u64> {
    let mut values = vec![
        duct.time_s(),
        duct.wave_energy_j().unwrap(),
        duct.vortex_energy_j().unwrap(),
        duct.represented_energy_j().unwrap(),
        duct.cumulative_mutual_friction_heat_j(),
        duct.cumulative_normal_viscous_heat_j(),
    ];
    for i in 0..duct.cells() {
        let cell = duct.cell(i).unwrap();
        values.extend([
            cell.pressure_perturbation_pa,
            cell.temperature_perturbation_k,
            cell.density_perturbation_kg_m3,
            cell.specific_entropy_perturbation_j_kg_k,
        ]);
    }
    for i in 0..=duct.cells() {
        let face = duct.face(i).unwrap();
        values.extend([face.normal_velocity_m_s, face.superfluid_velocity_m_s]);
        if let Some(density) = duct.vortex_line_density_per_m2(i) {
            values.push(density);
        }
        if let Some(spectrum) = duct.normal_viscous_spectrum() {
            for mode in 0..spectrum.modes().len() {
                if let Some(velocity) = duct.normal_profile_mode_velocity_m_s(i, mode) {
                    values.push(velocity);
                }
            }
        }
    }
    values.extend(duct.last_dissipated_heat_j());
    values.into_iter().map(f64::to_bits).collect()
}

#[test]
fn zero_flow_decays_each_seed_analytically_and_routes_its_line_energy_once() {
    let mut duct = channel(5, false);
    let initial = [0.0, 2e7, 1e8, 3e8];
    seed(&mut duct, &initial);
    let volume = duct.area_m2() * duct.length_m() / duct.cells() as f64;
    let initial_energy = LINE_ENERGY * initial.iter().sum::<f64>() * volume;
    close(
        duct.vortex_energy_j().unwrap(),
        initial_energy,
        initial_energy,
    );
    assert_eq!(duct.vortex_line_density_per_m2(0), None);
    assert_eq!(duct.vortex_line_density_per_m2(5), None);
    let dt = duct.maximum_step_s() * 0.7;
    let mut total_heat = 0.0;
    for _ in 0..8 {
        let before = duct.vortex_energy_j().unwrap();
        let (report, accepted_heat) = duct
            .step_with_heat(
                dt,
                BoundaryVelocities::CLOSED,
                FRICTION,
                || false,
                |_, heat| Ok::<_, ()>(heat.to_vec()),
            )
            .unwrap();
        assert_eq!(report.wave_energy_before_j, 0.0);
        assert_eq!(report.wave_energy_after_j, 0.0);
        assert_eq!(report.vortex_energy_before_j, before);
        assert_eq!(
            report.vortex_energy_after_j,
            duct.vortex_energy_j().unwrap()
        );
        close(
            accepted_heat.iter().sum(),
            before - report.vortex_energy_after_j,
            initial_energy,
        );
        assert_eq!(accepted_heat[0], 0.0);
        total_heat += report.mutual_friction_heat_j;
        for (i, &density) in initial.iter().enumerate() {
            let expected = density / (1.0 + DECAY * density * duct.time_s());
            close(
                duct.vortex_line_density_per_m2(i + 1).unwrap(),
                expected,
                density.max(1.0),
            );
        }
    }
    close(
        duct.represented_energy_j().unwrap() + total_heat,
        initial_energy,
        initial_energy,
    );
    close(
        duct.cumulative_dissipated_heat_j(),
        total_heat,
        initial_energy,
    );
}

#[test]
fn a_zero_seed_stays_empty_and_preserves_the_native_no_friction_wave() {
    for profile in [false, true] {
        let mut bare = channel(12, profile);
        let mut seeded = channel(12, profile);
        pulse(&mut bare);
        pulse(&mut seeded);
        seed(&mut seeded, &[0.0; 11]);
        let dt = 0.6 * bare.maximum_step_s();
        for _ in 0..5 {
            bare.step(dt, BoundaryVelocities::CLOSED, MutualFriction::None, || {
                false
            })
            .unwrap();
            let report = seeded
                .step(dt, BoundaryVelocities::CLOSED, FRICTION, || false)
                .unwrap();
            assert_eq!(report.mutual_friction_heat_j, 0.0);
            assert_eq!(report.vortex_energy_after_j, 0.0);
            for i in 0..bare.cells() {
                assert_eq!(bare.cell(i), seeded.cell(i));
            }
            for i in 0..=bare.cells() {
                assert_eq!(bare.face(i), seeded.face(i));
            }
        }
    }
}

#[test]
fn wave_tangle_and_viscous_profile_share_one_storage_and_heat_account() {
    for profile in [false, true] {
        let mut duct = channel(16, profile);
        pulse(&mut duct);
        let densities: Vec<_> = (1..16).map(|i| if i < 8 { 1e7 } else { 1e8 }).collect();
        seed(&mut duct, &densities);
        let initial = duct.represented_energy_j().unwrap();
        let mut heat = 0.0;
        let mut grew = false;
        for _ in 0..80 {
            let report = duct
                .step(
                    0.5 * duct.maximum_step_s(),
                    BoundaryVelocities::CLOSED,
                    FRICTION,
                    || false,
                )
                .unwrap();
            grew |= report.vortex_energy_after_j > report.vortex_energy_before_j;
            assert!(report.mutual_friction_heat_j >= 0.0);
            assert!(report.normal_viscous_heat_j >= 0.0);
            heat += report.dissipated_heat_j();
            close(report.wave_balance_defect_j, 0.0, initial);
            close(
                duct.represented_energy_j().unwrap() + heat,
                initial,
                initial,
            );
        }
        assert!(
            grew,
            "the pulse must actually create represented vortex length"
        );
        assert!(heat > 0.0);
        assert!(duct.vortex_line_density_per_m2(5).unwrap() > 0.0);
        assert_ne!(
            duct.vortex_line_density_per_m2(5),
            duct.vortex_line_density_per_m2(10)
        );
        if profile {
            assert!(duct.cumulative_normal_viscous_heat_j() > 0.0);
        }
    }
}

#[test]
fn tangle_history_survives_probes_refusals_cancellation_and_retry() {
    let mut duct = channel(9, true);
    pulse(&mut duct);
    seed(&mut duct, &[1e8; 8]);
    let dt = 0.6 * duct.maximum_step_s();
    duct.step(dt, BoundaryVelocities::CLOSED, FRICTION, || false)
        .unwrap();
    let before = snapshot(&duct);
    let old_temperature = [
        duct.cell(0).unwrap().temperature_perturbation_k,
        duct.cell(8).unwrap().temperature_perturbation_k,
    ];
    let response = duct
        .counterflow_ends_response(dt, [0.0; 2], BodyForce::None, FRICTION, || false)
        .unwrap();
    assert_eq!(snapshot(&duct), before);
    duct.impermeable_ends_response(dt, [0.0; 2], BodyForce::None, FRICTION, || false)
        .unwrap();
    assert_eq!(snapshot(&duct), before);

    let refused = duct.step_with_heat(
        dt,
        BoundaryVelocities::CLOSED,
        FRICTION,
        || false,
        |_, _| Err::<(), _>("caloric owner refuses"),
    );
    assert_eq!(refused, Err(HeatStepError::Heat("caloric owner refuses")));
    assert_eq!(snapshot(&duct), before);
    let prepared = Cell::new(false);
    let cancelled = duct.step_with_heat(
        dt,
        BoundaryVelocities::CLOSED,
        FRICTION,
        || prepared.get(),
        |_, _| {
            prepared.set(true);
            Ok::<_, ()>(())
        },
    );
    assert_eq!(
        cancelled,
        Err(HeatStepError::Channel(ChannelError::Cancelled))
    );
    assert_eq!(snapshot(&duct), before);

    let polls = Cell::new(0);
    duct.counterflow_ends_response(dt, [0.0; 2], BodyForce::None, FRICTION, || {
        polls.set(polls.get() + 1);
        false
    })
    .unwrap();
    for stop in [1, polls.get() / 2] {
        let seen = Cell::new(0);
        let result = duct.step(dt, BoundaryVelocities::CLOSED, FRICTION, || {
            seen.set(seen.get() + 1);
            seen.get() == stop
        });
        assert_eq!(result, Err(ChannelError::Cancelled));
        assert_eq!(snapshot(&duct), before);
    }
    let heat = [-0.0001, 0.00015];
    let boundaries = duct.counterflow_boundary_velocities(heat).unwrap();
    duct.step(dt, boundaries, FRICTION, || false).unwrap();
    for (i, cell) in [0, 8].into_iter().enumerate() {
        let expected = response.reference_temperature_perturbation_k[i]
            - response.temperature_heat_impedance_k_w[i][0] * heat[0]
            - response.temperature_heat_impedance_k_w[i][1] * heat[1];
        let actual =
            0.5 * (old_temperature[i] + duct.cell(cell).unwrap().temperature_perturbation_k);
        close(actual, expected, 0.001);
    }
    assert_ne!(snapshot(&duct), before);
}

#[test]
fn configuration_and_exclusive_selector_refuse_without_changing_accepted_state() {
    let mut duct = channel(5, false);
    pulse(&mut duct);
    let dt = 0.5 * duct.maximum_step_s();
    assert!(
        duct.step(dt, BoundaryVelocities::CLOSED, FRICTION, || false)
            .is_err()
    );
    seed(&mut duct, &[1e8; 4]);
    let before = snapshot(&duct);
    let different_material = TwoFluidState::try_new(TwoFluidProperties {
        temperature_k: 1.9,
        density_kg_m3: 4.0,
        superfluid_fraction: 0.75,
        entropy_j_kg_k: 4.0,
        specific_heat_cp_j_kg_k: 16.0,
        adiabatic_sound_speed_m_s: 10.0,
        thermal_expansion_per_k: -0.2,
        normal_viscosity_pa_s: 1e-4,
    })
    .unwrap();
    assert!(
        duct.configure_vinen_tangle(model(different_material), &[1e8; 4], || false)
            .is_err()
    );
    assert_eq!(snapshot(&duct), before);
    for densities in [
        &[1e8; 3][..],
        &[1e8, f64::NAN, 1e8, 1e8],
        &[1e8, -1.0, 1e8, 1e8],
    ] {
        assert!(
            duct.configure_vinen_tangle(model(material()), densities, || false)
                .is_err()
        );
        assert_eq!(snapshot(&duct), before);
    }
    assert_eq!(
        duct.configure_vinen_tangle(model(material()), &[1e7; 4], || true),
        Err(ChannelError::Cancelled)
    );
    assert_eq!(snapshot(&duct), before);
    for friction in [
        MutualFriction::None,
        MutualFriction::Linear { gamma_kg_m3_s: 1.0 },
        MutualFriction::GorterMellink { a_m_s_kg: 1.0 },
        MutualFriction::IsotropicVortexTangle {
            circulation_quantum_m2_s: 1e-7,
            line_density_per_m2: 1e8,
            dissipative_b: 1.2,
        },
    ] {
        assert!(
            duct.step(dt, BoundaryVelocities::CLOSED, friction, || false)
                .is_err()
        );
        assert_eq!(snapshot(&duct), before);
    }
    duct.step(dt, BoundaryVelocities::CLOSED, FRICTION, || false)
        .unwrap();
    let stepped = snapshot(&duct);
    assert!(
        duct.configure_vinen_tangle(model(material()), &[1e8; 4], || false)
            .is_err()
    );
    assert_eq!(snapshot(&duct), stepped);
}

#[test]
fn exhausted_local_source_budget_refuses_the_whole_wave_and_can_retry() {
    let mut duct = channel(9, false);
    pulse(&mut duct);
    let base = model(material());
    let limited = VinenTangle::try_new(
        material(),
        base.parameters(),
        VinenStepControls {
            maximum_fractional_step: 1e-4,
            maximum_substeps: 1,
        },
    )
    .unwrap();
    duct.configure_vinen_tangle(limited, &[0.0, 1e8, 1e8, 1e8, 1e8, 1e8, 1e8, 1e8], || false)
        .unwrap();
    let before = snapshot(&duct);
    let dt = 0.5 * duct.maximum_step_s();
    assert!(matches!(
        duct.step(dt, BoundaryVelocities::CLOSED, FRICTION, || false),
        Err(ChannelError::Material(_))
    ));
    assert_eq!(snapshot(&duct), before);
    duct.step(dt / 4096.0, BoundaryVelocities::CLOSED, FRICTION, || false)
        .unwrap();
    assert!(duct.time_s() > 0.0);
    assert_eq!(duct.vortex_line_density_per_m2(1), Some(0.0));
}
