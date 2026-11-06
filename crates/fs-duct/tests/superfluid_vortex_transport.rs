//! Native conservative vortex transport: ports, positivity, accuracy and rollback.
//! All coefficients are authored numerical fixtures, not helium calibration.

use fs_duct::superfluid::{
    BodyForce, BoundaryVelocities, ChannelError, HeatStepError, LinearizationLimits,
    MutualFriction, TwoFluidChannel, VortexTransportBoundaries, VortexTransportBoundary,
    VortexTransportControls,
};
use fs_material::superfluid::vinen::{VinenParameters, VinenStepControls, VinenTangle};
use fs_material::superfluid::{NormalViscousSpectrum, TwoFluidProperties, TwoFluidState};
use std::cell::Cell;

const FRICTION: MutualFriction = MutualFriction::EvolvingVinenTangle;
const LINE_ENERGY: f64 = 1e-15;
const CLOSED: VortexTransportBoundaries = VortexTransportBoundaries::CLOSED;
const CONTROLS: VortexTransportControls = VortexTransportControls {
    maximum_outgoing_courant: 0.45,
    maximum_substeps: 128,
};

fn channel(n: usize, density: &[f64], active_source: bool, profile: bool) -> TwoFluidChannel {
    let material = TwoFluidState::try_new(TwoFluidProperties {
        temperature_k: 2.0,
        density_kg_m3: 4.0,
        superfluid_fraction: 0.75,
        entropy_j_kg_k: 4.0,
        specific_heat_cp_j_kg_k: 16.0,
        adiabatic_sound_speed_m_s: 10.0,
        thermal_expansion_per_k: -0.2,
        normal_viscosity_pa_s: 1e-4,
    })
    .unwrap();
    let limits = LinearizationLimits {
        max_pressure_perturbation_pa: 1.0,
        max_temperature_perturbation_k: 0.02,
        max_component_speed_m_s: 0.1,
    };
    let mut duct = if profile {
        TwoFluidChannel::try_new_with_normal_profile(
            material,
            1.0,
            0.01,
            n,
            limits,
            NormalViscousSpectrum::circular_capillary(material, 0.01, 3).unwrap(),
        )
    } else {
        TwoFluidChannel::try_new(material, 1.0, 0.01, n, limits)
    }
    .unwrap();
    let model = VinenTangle::try_new(
        material,
        VinenParameters {
            circulation_quantum_m2_s: 1e-7,
            dissipative_b: if active_source { 1.2 } else { 0.0 },
            production_coefficient: if active_source { 0.8 } else { 0.0 },
            // Keep a positive admitted decay coefficient; pure-transport fixtures
            // make its physical change far below the transport accuracy measured.
            decay_coefficient_m2_s: if active_source { 1e-7 } else { 1e-30 },
            line_energy_j_m: LINE_ENERGY,
        },
        VinenStepControls {
            maximum_fractional_step: 0.04,
            maximum_substeps: 512,
        },
    )
    .unwrap();
    duct.configure_vinen_tangle(model, density, || false)
        .unwrap();
    duct
}

fn pulse(duct: &mut TwoFluidChannel) {
    let n = duct.cells();
    let shape: Vec<_> = (0..=n)
        .map(|i| (std::f64::consts::PI * i as f64 / n as f64).sin())
        .collect();
    duct.initialize(
        &vec![0.0; n],
        &vec![0.0; n],
        &shape.iter().map(|v| 0.003 * v).collect::<Vec<_>>(),
        &shape.iter().map(|v| -0.001 * v).collect::<Vec<_>>(),
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
    let mut data = vec![
        duct.time_s(),
        duct.wave_energy_j().unwrap(),
        duct.vortex_energy_j().unwrap(),
        duct.cumulative_mutual_friction_heat_j(),
        duct.cumulative_normal_viscous_heat_j(),
    ];
    for cell in 0..duct.cells() {
        let sample = duct.cell(cell).unwrap();
        data.extend([
            sample.pressure_perturbation_pa,
            sample.temperature_perturbation_k,
        ]);
    }
    for face in 0..=duct.cells() {
        let sample = duct.face(face).unwrap();
        data.extend([sample.normal_velocity_m_s, sample.superfluid_velocity_m_s]);
        data.extend(duct.vortex_line_density_per_m2(face));
        if let Some(spectrum) = duct.normal_viscous_spectrum() {
            for mode in 0..spectrum.modes().len() {
                data.extend(duct.normal_profile_mode_velocity_m_s(face, mode));
            }
        }
    }
    data.extend(duct.last_dissipated_heat_j());
    data.into_iter().map(f64::to_bits).collect()
}

fn reservoir(density: f64) -> VortexTransportBoundary {
    VortexTransportBoundary::Reservoir {
        incoming_line_density_per_m2: density,
    }
}

#[test]
fn closed_transport_moves_sharp_fronts_without_negative_density_or_lost_line_energy() {
    let n = 33;
    let initial: Vec<_> = (0..n - 1)
        .map(|i| {
            if (4..10).contains(&i) || (20..26).contains(&i) {
                1e8
            } else {
                0.0
            }
        })
        .collect();
    let mut duct = channel(n, &initial, false, false);
    let drift: Vec<_> = (0..n)
        .map(|edge| {
            if edge == 0 || edge == n - 1 {
                0.0
            } else if edge < n / 2 {
                30.0
            } else {
                -30.0
            }
        })
        .collect();
    duct.configure_vortex_line_transport(&drift, CLOSED, CONTROLS, || false)
        .unwrap();
    let before = duct.represented_energy_j().unwrap();
    let mut heat = 0.0;
    for _ in 0..30 {
        let report = duct
            .step(
                0.8 * duct.maximum_step_s(),
                BoundaryVelocities::CLOSED,
                FRICTION,
                || false,
            )
            .unwrap();
        assert!(report.vortex_transport_substeps > 2);
        assert_eq!(report.vortex_transport_boundary_line_m, [0.0; 2]);
        assert_eq!(report.vortex_transport_boundary_energy_j, [0.0; 2]);
        close(
            report.vortex_transport_line_balance_defect_m * LINE_ENERGY,
            0.0,
            before,
        );
        close(report.wave_balance_defect_j, 0.0, before);
        heat += report.dissipated_heat_j();
        for face in 1..n {
            let density = duct.vortex_line_density_per_m2(face).unwrap();
            assert!(density.is_finite() && density >= 0.0);
        }
    }
    close(duct.represented_energy_j().unwrap() + heat, before, before);
    assert!(duct.vortex_line_density_per_m2(n / 2).unwrap() > 1e8);
    assert!(duct.vortex_line_density_per_m2(5).unwrap() < initial[4]);
}

#[test]
fn reservoir_inflow_without_outgoing_velocity_and_outflow_have_paired_energy_ports() {
    let mut duct = channel(2, &[0.0], false, false);
    let boundaries = VortexTransportBoundaries {
        left: reservoir(1e8),
        right: VortexTransportBoundary::Reflecting,
    };
    duct.configure_vortex_line_transport(&[0.8, 0.0], boundaries, CONTROLS, || false)
        .unwrap();
    let dt = 0.6 * duct.maximum_step_s();
    let report = duct
        .step(dt, BoundaryVelocities::CLOSED, FRICTION, || false)
        .unwrap();
    let line = 0.8 * 1e8 * dt * duct.area_m2();
    close(report.vortex_transport_boundary_line_m[0], line, line);
    assert_eq!(report.vortex_transport_boundary_line_m[1], 0.0);
    assert_eq!(
        report.vortex_transport_boundary_energy_j[0],
        LINE_ENERGY * report.vortex_transport_boundary_line_m[0]
    );
    close(
        duct.vortex_energy_j().unwrap() + report.dissipated_heat_j(),
        LINE_ENERGY * line,
        LINE_ENERGY * line,
    );
    close(report.wave_balance_defect_j, 0.0, LINE_ENERGY * line);

    let mut outgoing = channel(2, &[1e8], false, false);
    outgoing
        .configure_vortex_line_transport(
            &[0.0, 1.0],
            VortexTransportBoundaries {
                left: VortexTransportBoundary::Reflecting,
                right: reservoir(0.0),
            },
            VortexTransportControls {
                maximum_outgoing_courant: 0.01,
                ..CONTROLS
            },
            || false,
        )
        .unwrap();
    let initial = outgoing.vortex_energy_j().unwrap();
    let report = outgoing
        .step(dt, BoundaryVelocities::CLOSED, FRICTION, || false)
        .unwrap();
    let exact = 1e8 * (-dt / 0.5).exp(); // One uniform dual volume with open outflow.
    assert!((outgoing.vortex_line_density_per_m2(1).unwrap() - exact).abs() < 3e-6 * exact);
    assert!(report.vortex_transport_boundary_line_m[1] < 0.0);
    assert!(report.dissipated_heat_j() < 1e-20 * initial);
    close(
        outgoing.vortex_energy_j().unwrap() + report.dissipated_heat_j()
            - report
                .vortex_transport_boundary_energy_j
                .iter()
                .sum::<f64>(),
        initial,
        initial,
    );
}

#[test]
fn throughflow_and_reversal_preserve_a_uniform_tangle_and_the_two_port_signs() {
    let mut duct = channel(8, &[1e8; 7], false, false);
    let initial = duct.vortex_energy_j().unwrap();
    for speed in [2.0, -3.0] {
        duct.configure_vortex_line_transport(
            &[speed; 8],
            VortexTransportBoundaries {
                left: reservoir(1e8),
                right: reservoir(1e8),
            },
            CONTROLS,
            || false,
        )
        .unwrap();
        let dt = 0.7 * duct.maximum_step_s();
        let report = duct
            .step(dt, BoundaryVelocities::CLOSED, FRICTION, || false)
            .unwrap();
        let left = speed * 1e8 * dt * duct.area_m2();
        close(report.vortex_transport_boundary_line_m[0], left, left.abs());
        close(
            report.vortex_transport_boundary_line_m[1],
            -left,
            left.abs(),
        );
        for face in 1..8 {
            close(duct.vortex_line_density_per_m2(face).unwrap(), 1e8, 1e8);
        }
        close(duct.vortex_energy_j().unwrap(), initial, initial);
        close(report.wave_balance_defect_j, 0.0, initial);
    }
}

fn advection_error(n: usize) -> f64 {
    let shape = |x: f64| {
        let offset = (x - 0.3) / 0.14;
        if offset.abs() < 1.0 {
            0.25 * (1.0 + (std::f64::consts::PI * offset).cos()).powi(2)
        } else {
            0.0
        }
    };
    let initial: Vec<_> = (1..n)
        .map(|face| 1e8 * shape(face as f64 / n as f64))
        .collect();
    let mut duct = channel(n, &initial, false, false);
    let mut drift = vec![4.0; n];
    drift[0] = 0.0;
    drift[n - 1] = 0.0;
    duct.configure_vortex_line_transport(&drift, CLOSED, CONTROLS, || false)
        .unwrap();
    let duration = 0.025;
    let steps = (duration / (0.7 * duct.maximum_step_s())).ceil() as usize;
    for _ in 0..steps {
        duct.step(
            duration / steps as f64,
            BoundaryVelocities::CLOSED,
            FRICTION,
            || false,
        )
        .unwrap();
    }
    (1..n)
        .map(|face| {
            (duct.vortex_line_density_per_m2(face).unwrap() / 1e8
                - shape(face as f64 / n as f64 - 4.0 * duration))
            .abs()
        })
        .sum::<f64>()
        / n as f64
}

#[test]
fn smooth_translation_converges_under_spatial_and_temporal_refinement() {
    let coarse = advection_error(48);
    let fine = advection_error(96);
    assert!(
        fine < 0.45 * coarse,
        "coarse L1={coarse:e}, fine L1={fine:e}"
    );
}

#[test]
fn advected_line_density_changes_native_drag_and_keeps_the_profile_heat_balance() {
    let initial: Vec<_> = (0..15).map(|i| if i < 5 { 1e8 } else { 0.0 }).collect();
    let mut moving = channel(16, &initial, true, true);
    let mut local = channel(16, &initial, true, true);
    pulse(&mut moving);
    pulse(&mut local);
    let mut drift = [30.0; 16];
    drift[0] = 0.0;
    drift[15] = 0.0;
    moving
        .configure_vortex_line_transport(&drift, CLOSED, CONTROLS, || false)
        .unwrap();
    let before = moving.represented_energy_j().unwrap();
    let mut heat = 0.0;
    for _ in 0..8 {
        let dt = 0.6 * moving.maximum_step_s();
        let report = moving
            .step(dt, BoundaryVelocities::CLOSED, FRICTION, || false)
            .unwrap();
        local
            .step(dt, BoundaryVelocities::CLOSED, FRICTION, || false)
            .unwrap();
        heat += report.dissipated_heat_j();
        close(report.wave_balance_defect_j, 0.0, before);
        close(
            moving.represented_energy_j().unwrap() + heat,
            before,
            before,
        );
    }
    assert_eq!(local.vortex_line_density_per_m2(10), Some(0.0));
    assert!(moving.vortex_line_density_per_m2(10).unwrap() > 0.0);
    assert_ne!(moving.face(10), local.face(10));
    assert!(moving.cumulative_normal_viscous_heat_j() > 0.0);
}

#[test]
fn queries_heat_refusal_cancellation_and_retry_preserve_the_accepted_transport_history() {
    let make = || {
        let mut duct = channel(9, &[1e8, 1e8, 1e8, 0.0, 0.0, 0.0, 0.0, 0.0], true, true);
        pulse(&mut duct);
        duct.configure_vortex_line_transport(
            &[10.0; 9],
            VortexTransportBoundaries {
                left: reservoir(1e7),
                right: reservoir(0.0),
            },
            CONTROLS,
            || false,
        )
        .unwrap();
        duct.step(
            0.5 * duct.maximum_step_s(),
            BoundaryVelocities::CLOSED,
            FRICTION,
            || false,
        )
        .unwrap();
        duct
    };
    let mut duct = make();
    let mut direct = make();
    let before = snapshot(&duct);
    let dt = 0.5 * duct.maximum_step_s();
    let old_theta = [
        duct.cell(0).unwrap().temperature_perturbation_k,
        duct.cell(8).unwrap().temperature_perturbation_k,
    ];
    let response = duct
        .counterflow_ends_response(dt, [0.0; 2], BodyForce::None, FRICTION, || false)
        .unwrap();
    duct.impermeable_ends_response(dt, [0.0; 2], BodyForce::None, FRICTION, || false)
        .unwrap();
    assert_eq!(snapshot(&duct), before);
    let polls = Cell::new(0);
    let refused = duct.step_with_heat(
        dt,
        BoundaryVelocities::CLOSED,
        FRICTION,
        || {
            polls.set(polls.get() + 1);
            false
        },
        |_, _| Err::<(), _>("heat owner refusal"),
    );
    assert_eq!(refused, Err(HeatStepError::Heat("heat owner refusal")));
    assert_eq!(snapshot(&duct), before);
    for stop in [1, polls.get() / 3, 2 * polls.get() / 3] {
        let count = Cell::new(0);
        let result = duct.step(dt, BoundaryVelocities::CLOSED, FRICTION, || {
            count.set(count.get() + 1);
            count.get() == stop
        });
        assert_eq!(result, Err(ChannelError::Cancelled));
        assert_eq!(snapshot(&duct), before);
    }
    let prepared = Cell::new(false);
    let result = duct.step_with_heat(
        dt,
        BoundaryVelocities::CLOSED,
        FRICTION,
        || prepared.get(),
        |_, _| {
            prepared.set(true);
            Ok::<_, ()>(())
        },
    );
    assert_eq!(result, Err(HeatStepError::Channel(ChannelError::Cancelled)));
    assert_eq!(snapshot(&duct), before);
    let heat = [-0.0001, 0.00015];
    let boundary = duct.counterflow_boundary_velocities(heat).unwrap();
    let retry = duct.step(dt, boundary, FRICTION, || false).unwrap();
    let fresh = direct.step(dt, boundary, FRICTION, || false).unwrap();
    assert_eq!(retry, fresh);
    assert_eq!(snapshot(&duct), snapshot(&direct));
    for (end, cell) in [0, 8].into_iter().enumerate() {
        let predicted = response.reference_temperature_perturbation_k[end]
            - response.temperature_heat_impedance_k_w[end][0] * heat[0]
            - response.temperature_heat_impedance_k_w[end][1] * heat[1];
        close(
            0.5 * (old_theta[end] + duct.cell(cell).unwrap().temperature_perturbation_k),
            predicted,
            0.001,
        );
    }
}

#[test]
fn configuration_budget_and_zero_drift_preserve_the_existing_source_path() {
    let initial = [1e8; 3];
    let mut duct = channel(4, &initial, true, false);
    let mut local = channel(4, &initial, true, false);
    pulse(&mut duct);
    pulse(&mut local);
    duct.configure_vortex_line_transport(&[0.0; 4], CLOSED, CONTROLS, || false)
        .unwrap();
    let model = *duct.vinen_tangle().unwrap();
    duct.configure_vinen_tangle(model, &initial, || false)
        .unwrap();
    assert_eq!(duct.vortex_transport_velocity_m_s(1), Some(0.0));
    let dt = 0.7 * duct.maximum_step_s();
    assert_eq!(
        duct.step(dt, BoundaryVelocities::CLOSED, FRICTION, || false)
            .unwrap(),
        local
            .step(dt, BoundaryVelocities::CLOSED, FRICTION, || false)
            .unwrap()
    );
    assert_eq!(snapshot(&duct), snapshot(&local));
    let before = snapshot(&duct);
    assert!(
        duct.configure_vortex_line_transport(&[0.0; 3], CLOSED, CONTROLS, || false)
            .is_err()
    );
    assert!(
        duct.configure_vortex_line_transport(&[1.0; 4], CLOSED, CONTROLS, || false)
            .is_err()
    );
    assert!(
        duct.configure_vortex_line_transport(&[0.0, f64::NAN, 0.0, 0.0], CLOSED, CONTROLS, || {
            false
        })
        .is_err()
    );
    assert!(
        duct.configure_vortex_line_transport(
            &[0.0; 4],
            VortexTransportBoundaries {
                left: reservoir(-1.0),
                right: reservoir(0.0)
            },
            CONTROLS,
            || false
        )
        .is_err()
    );
    assert!(
        duct.configure_vortex_line_transport(
            &[0.0; 4],
            CLOSED,
            VortexTransportControls {
                maximum_outgoing_courant: 0.51,
                ..CONTROLS
            },
            || false
        )
        .is_err()
    );
    assert_eq!(
        duct.configure_vortex_line_transport(&[0.0, 1.0, 1.0, 0.0], CLOSED, CONTROLS, || true),
        Err(ChannelError::Cancelled)
    );
    assert_eq!(duct.vortex_transport_velocity_m_s(1), Some(0.0));
    assert_eq!(snapshot(&duct), before);
    let drift = [0.0, 1000.0, 1000.0, 0.0];
    duct.configure_vortex_line_transport(
        &drift,
        CLOSED,
        VortexTransportControls {
            maximum_substeps: 1,
            ..CONTROLS
        },
        || false,
    )
    .unwrap();
    assert!(
        duct.step(dt, BoundaryVelocities::CLOSED, FRICTION, || false)
            .is_err()
    );
    assert_eq!(snapshot(&duct), before);
    duct.configure_vortex_line_transport(&drift, CLOSED, CONTROLS, || false)
        .unwrap();
    duct.step(dt, BoundaryVelocities::CLOSED, FRICTION, || false)
        .unwrap();
    assert!(duct.time_s() > f64::from_bits(before[0]));
}
