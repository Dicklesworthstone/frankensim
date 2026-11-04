use super::*;
use crate::dg_port::{DgBranch, DgTerminalChannelInput};
use crate::modal_acoustic_time::{ModalAcousticMode, ModalAcousticTimeBudget};
use crate::render::plate::impact::radiation::{Model, Pole};
use fs_duct::dg::{DgMesh, DgSidePort, MAX_DG_PROPAGATOR_DOFS};
use fs_duct::{Duct, Segment};
use fs_material::gas::{GasSpec, GasState};
use fs_vfit::mimo::{MimoModel, load::MimoLoadLimits};
use fs_vfit::model::PoleTerm;

const INLET: PortId = PortId(7);
const OUTLET: PortId = PortId(3);
const SIDE: PortId = PortId(19);
const WORK: usize = 100_000_000;
const DT: f64 = 1.0 / 48_000.0;
const STEP: f64 = 0.0002;

fn gas() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap()
}

pub(in super::super) fn oscillator() -> CoupledModalSystem {
    let gate = CancelGate::new_clock_free();
    let model = ModalAcousticTimeModel::try_new(
        48_000,
        vec![ModalAcousticMode {
            angular_frequency_rad_s: 1800.0,
            damping_ratio: 0.003,
            pressure_per_modal_velocity: C64::ZERO,
        }],
        ModalAcousticTimeBudget::audible_reference(),
    )
    .unwrap();
    CoupledModalSystem::new_with_columns(
        vec![model],
        vec![ModalColumnConnection {
            column: vec![0.6],
            stiffness_n_m: 20_000.0,
            damping_n_s_m: 4.0,
            rest_extension_m: 0.0,
        }],
        ModalCouplingConfig {
            max_modes: 4,
            max_connections: 4,
            max_setup_terms: WORK,
            nyquist_guard_fraction: 0.9,
            maximum_total_energy_j: 100.0,
            maximum_abs_pressure_pa: 1e6,
            maximum_abs_connection_force_n: 1e6,
            solve_relative_tolerance: 1e-11,
            energy_absolute_tolerance_j: 1e-12,
            energy_relative_tolerance: 1e-9,
        },
        &gate,
    )
    .unwrap()
    .with_material_relaxation(
        vec![fs_phs::RelaxationBranch {
            projection: vec![1.0, 0.0],
            stiffness: 30_000.0,
            relaxation_time_s: 0.002,
        }],
        Some(vec![0.004]),
        &gate,
    )
    .unwrap()
    .with_radiation_load(
        &Model {
            ports: 1,
            poles: vec![Pole {
                omega: 2200.0,
                zeta: 0.02,
                coupling: vec![100.0],
            }],
        },
        &[0],
        1,
        &gate,
    )
    .unwrap()
}

fn bare() -> DgDuct {
    DgDuct::new(
        &Duct {
            segments: vec![Segment::Cylinder {
                radius: 0.007,
                length: 0.16,
            }],
        },
        &gas(),
        DgMesh {
            elements: 1,
            order: 1,
        },
        DgBoundary::Port(INLET),
        DgBoundary::IdealOpen,
    )
    .unwrap()
}

pub(in super::super) fn loaded(circular: bool, reverse: bool) -> DgPortLoad {
    let bore = DgDuct::new(
        &Duct {
            segments: vec![Segment::Cylinder {
                radius: 0.007,
                length: 0.16,
            }],
        },
        &gas(),
        DgMesh {
            elements: 2,
            order: 1,
        },
        DgBoundary::Port(INLET),
        DgBoundary::Port(OUTLET),
    )
    .unwrap()
    .with_side_ports(&[DgSidePort {
        id: SIDE,
        downstream_element: 1,
    }])
    .unwrap();
    let branch = DgDuct::new(
        &Duct {
            segments: vec![Segment::Cylinder {
                radius: 0.004,
                length: 0.35,
            }],
        },
        &gas(),
        DgMesh {
            elements: 1,
            order: 1,
        },
        DgBoundary::Port(INLET),
        DgBoundary::Port(OUTLET),
    )
    .unwrap()
    .with_viscothermal(16, [50.0, 2000.0], 0.02)
    .unwrap();
    let zc = gas().density * gas().sound_speed / (core::f64::consts::PI * 0.007_f64.powi(2));
    let mut residues: Vec<_> = [0.25, 0.06, 0.06, 0.2]
        .into_iter()
        .map(|value| C64::from_re(value * zc * 600.0))
        .collect();
    let ports = if reverse {
        residues.reverse();
        vec![OUTLET, SIDE]
    } else {
        vec![SIDE, OUTLET]
    };
    // Strictly proper, reciprocal and mutually coupled: no D inverse exists.
    let terminal = PreparedMimoLoad::new(
        &MimoModel {
            ports: 2,
            poles: vec![PoleTerm::Real {
                pole: -600.0,
                residue: 0.0,
            }],
            residues: vec![residues],
            d: vec![0.0; 4],
            e: vec![0.0; 4],
            reciprocal: true,
        },
        DT,
        MimoLoadLimits {
            band_rad_s: [100.0, core::f64::consts::PI / DT],
            sampled_passivity_points: 33,
            max_states: 8,
            max_memory_bytes: 1 << 20,
        },
    )
    .unwrap();
    DgPortLoad::new_with_terminal_branches(
        bore,
        DT,
        1e-11,
        MAX_DG_PROPAGATOR_DOFS,
        ports,
        terminal,
        vec![DgBranch {
            port: SIDE,
            duct: branch,
        }],
    )
    .unwrap()
    .with_terminal_channels(
        vec![
            DgTerminalChannelInput {
                port: SIDE,
                kind: if circular {
                    DgTerminalChannelKind::Circular {
                        radius_m: 0.0008,
                        length_m: 0.02,
                        radial_cells: 8,
                        maximum_mach: 0.1,
                    }
                } else {
                    DgTerminalChannelKind::Slit {
                        width_m: 0.01,
                        height_m: 0.0005,
                        length_m: 0.02,
                        transverse_cells: 8,
                        maximum_mach: 0.1,
                    }
                },
            },
            DgTerminalChannelInput {
                port: OUTLET,
                kind: DgTerminalChannelKind::Closed,
            },
        ],
        8,
    )
    .unwrap()
}

fn close(actual: C64, expected: C64) {
    assert!(
        (actual - expected).abs() <= 3e-6 * actual.abs().max(expected.abs()) + 1e-25,
        "{actual:?} != {expected:?}"
    );
}
fn scalar(actual: f64, lo: f64, hi: f64) {
    let expected = (hi - lo) * 0.5 / STEP;
    assert!(
        (actual - expected).abs() <= 3e-6 * actual.abs().max(expected.abs()) + 1e-25,
        "{actual:e} != {expected:e}"
    );
}
fn vector(actual: &[C64], lo: &[C64], hi: &[C64]) {
    assert_eq!(actual.len(), lo.len());
    assert_eq!(actual.len(), hi.len());
    for ((&value, &lo), &hi) in actual.iter().zip(lo).zip(hi) {
        close(value, (hi - lo).scale(0.5 / STEP));
    }
}
fn difference_with_roundoff(actual: C64, lo: C64, hi: C64, primal_scale: f64) {
    let expected = (hi - lo).scale(0.5 / STEP);
    // A rigid closure has exactly zero flow and zero directly driven source
    // states. Differencing the joint solves amplifies roundoff by 1/STEP. Use
    // the uncancelled primal norm of the same physical quantity and units.
    let roundoff = 128.0 * f64::EPSILON * primal_scale / STEP;
    assert!(
        (actual - expected).abs() <= 3e-6 * actual.abs().max(expected.abs()) + roundoff + 1e-25,
        "{actual:?} != {expected:?} (roundoff {roundoff:e})"
    );
}
fn duct_derivatives(
    actual: &DgHarmonicDerivative,
    lo: &DgHarmonicResponse,
    hi: &DgHarmonicResponse,
) {
    vector(
        &actual.energy_coordinates,
        &lo.energy_coordinates,
        &hi.energy_coordinates,
    );
    vector(&actual.pressure_pa, &lo.pressure_pa, &hi.pressure_pa);
    vector(
        &actual.volume_flow_m3_s,
        &lo.volume_flow_m3_s,
        &hi.volume_flow_m3_s,
    );
    for (actual, lo, hi) in [
        (
            actual.mean_stored_energy_j,
            lo.mean_stored_energy_j,
            hi.mean_stored_energy_j,
        ),
        (
            actual.mean_wall_energy_j,
            lo.mean_wall_energy_j,
            hi.mean_wall_energy_j,
        ),
        (
            actual.wall_dissipation_w,
            lo.wall_dissipation_w,
            hi.wall_dissipation_w,
        ),
        (
            actual.molecular_dissipation_w,
            lo.molecular_dissipation_w,
            hi.molecular_dissipation_w,
        ),
        (
            actual.matched_absorption_w,
            lo.matched_absorption_w,
            hi.matched_absorption_w,
        ),
        (actual.input_power_w, lo.input_power_w, hi.input_power_w),
    ] {
        scalar(actual, lo, hi);
    }
    let flow_scale = lo
        .ports
        .iter()
        .chain(&hi.ports)
        .map(|port| port.inflow_m3_s.abs())
        .fold(0.0_f64, f64::max);
    for ((actual, lo), hi) in actual.ports.iter().zip(&lo.ports).zip(&hi.ports) {
        assert_eq!(actual.id, lo.id);
        assert_eq!(actual.id, hi.id);
        close(
            actual.pressure_pa,
            (hi.pressure_pa - lo.pressure_pa).scale(0.5 / STEP),
        );
        difference_with_roundoff(
            actual.inflow_m3_s,
            lo.inflow_m3_s,
            hi.inflow_m3_s,
            flow_scale,
        );
        scalar(
            actual.mean_input_power_w,
            lo.mean_input_power_w,
            hi.mean_input_power_w,
        );
    }
}

#[test]
fn full_loaded_tangent_retains_mutual_terminal_branch_memory_and_finite_gap_physics() {
    let gate = CancelGate::new_clock_free();
    let system = oscillator();
    let force = C64::new(1e-5, -4e-6);
    let slope = C64::new(2e-8, -1e-8);
    let column = [0.003];
    let memories = system.material_memory().to_vec();
    let air = system.radiation_state().to_vec();
    for circular in [false, true] {
        let mut prior = None;
        for reverse in [false, true] {
            let load = loaded(circular, reverse);
            let state = load.state().clone();
            let branch_state = load.branches()[0].state().clone();
            let terminal_state = load.terminal().unwrap().pole_states().to_vec();
            let channel_state = load.terminal_channels()[0]
                .channel()
                .unwrap()
                .velocities_m_s()
                .to_vec();
            let omega = 1300.0;
            let result = system
                .harmonic_response_with_load_frequency_derivative(
                    omega,
                    &[force],
                    Some(&[slope]),
                    &load,
                    &column,
                    WORK,
                    &gate,
                )
                .unwrap();
            assert_eq!(
                result.response,
                system
                    .harmonic_response_with_load(omega, &[force], &load, &column, WORK, &gate)
                    .unwrap()
            );
            let lo = system
                .harmonic_response_with_load(
                    omega - STEP,
                    &[force - slope.scale(STEP)],
                    &load,
                    &column,
                    WORK,
                    &gate,
                )
                .unwrap();
            let hi = system
                .harmonic_response_with_load(
                    omega + STEP,
                    &[force + slope.scale(STEP)],
                    &load,
                    &column,
                    WORK,
                    &gate,
                )
                .unwrap();
            assert_eq!(result.modal.response, result.response.modal);
            vector(
                &result.modal.displacement_derivative,
                &lo.modal.displacement_m_sqrt_kg,
                &hi.modal.displacement_m_sqrt_kg,
            );
            vector(
                &result.modal.velocity_derivative,
                &lo.modal.velocity_m_sqrt_kg_per_s,
                &hi.modal.velocity_m_sqrt_kg_per_s,
            );
            vector(
                &result.modal.material_memory_derivative,
                &lo.modal.material_memory,
                &hi.modal.material_memory,
            );
            vector(
                &result.modal.radiation_state_derivative,
                &lo.modal.radiation_state,
                &hi.modal.radiation_state,
            );
            vector(
                &result.modal.generalized_connection_forces_derivative,
                &lo.modal.generalized_connection_forces_n_per_sqrt_kg,
                &hi.modal.generalized_connection_forces_n_per_sqrt_kg,
            );
            scalar(
                result.modal.external_power_derivative,
                lo.modal.external_power_w,
                hi.modal.external_power_w,
            );
            duct_derivatives(&result.duct, &lo.duct, &hi.duct);
            for ((actual, lo), hi) in result.branches.iter().zip(&lo.branches).zip(&hi.branches) {
                duct_derivatives(&actual.duct, &lo.duct, &hi.duct);
                close(
                    actual.junction_flow_m3_s,
                    (hi.junction_flow_m3_s - lo.junction_flow_m3_s).scale(0.5 / STEP),
                );
                close(
                    actual.bore_pressure_pa,
                    (hi.bore_pressure_pa - lo.bore_pressure_pa).scale(0.5 / STEP),
                );
                close(
                    actual.inlet_pressure_pa,
                    (hi.inlet_pressure_pa - lo.inlet_pressure_pa).scale(0.5 / STEP),
                );
                scalar(actual.bore_power_w, lo.bore_power_w, hi.bore_power_w);
                scalar(actual.branch_power_w, lo.branch_power_w, hi.branch_power_w);
            }
            let terminal = result.terminal.as_ref().unwrap();
            let lo_terminal = lo.terminal.as_ref().unwrap();
            let hi_terminal = hi.terminal.as_ref().unwrap();
            let flow_scale = lo_terminal
                .outward_flow_m3_s
                .iter()
                .chain(&hi_terminal.outward_flow_m3_s)
                .map(|flow| flow.abs())
                .fold(
                    lo.inlet_flow_m3_s.abs().max(hi.inlet_flow_m3_s.abs()),
                    f64::max,
                );
            let state_scale = lo_terminal
                .state
                .iter()
                .chain(&hi_terminal.state)
                .map(|state| state.abs())
                .fold(0.0_f64, f64::max);
            for ((&actual, &lo), &hi) in terminal
                .state
                .iter()
                .zip(&lo_terminal.state)
                .zip(&hi_terminal.state)
            {
                difference_with_roundoff(actual, lo, hi, state_scale);
            }
            for ((&actual, &lo), &hi) in terminal
                .outward_flow_m3_s
                .iter()
                .zip(&lo_terminal.outward_flow_m3_s)
                .zip(&hi_terminal.outward_flow_m3_s)
            {
                difference_with_roundoff(actual, lo, hi, flow_scale);
            }
            vector(
                &terminal.pressure_pa,
                &lo_terminal.pressure_pa,
                &hi_terminal.pressure_pa,
            );
            for ((&actual, &lo), &hi) in terminal
                .mean_input_power_w
                .iter()
                .zip(&lo_terminal.mean_input_power_w)
                .zip(&hi_terminal.mean_input_power_w)
            {
                scalar(actual, lo, hi);
            }
            for ((actual, lo), hi) in result.channels.iter().zip(&lo.channels).zip(&hi.channels) {
                close(
                    actual.duct_pressure_pa,
                    (hi.duct_pressure_pa - lo.duct_pressure_pa).scale(0.5 / STEP),
                );
                close(
                    actual.exterior_pressure_pa,
                    (hi.exterior_pressure_pa - lo.exterior_pressure_pa).scale(0.5 / STEP),
                );
                close(
                    actual.pressure_drop_pa,
                    (hi.pressure_drop_pa - lo.pressure_drop_pa).scale(0.5 / STEP),
                );
                difference_with_roundoff(
                    actual.outward_flow_m3_s,
                    lo.outward_flow_m3_s,
                    hi.outward_flow_m3_s,
                    flow_scale,
                );
                scalar(actual.input_power_w, lo.input_power_w, hi.input_power_w);
                if let Some(field) = &actual.channel {
                    let (lo, hi) = (lo.channel.as_ref().unwrap(), hi.channel.as_ref().unwrap());
                    vector(
                        &field.velocities_m_s,
                        &lo.velocities_m_s,
                        &hi.velocities_m_s,
                    );
                    scalar(
                        field.mean_stored_energy_j,
                        lo.mean_stored_energy_j,
                        hi.mean_stored_energy_j,
                    );
                    scalar(
                        field.viscous_dissipation_w,
                        lo.viscous_dissipation_w,
                        hi.viscous_dissipation_w,
                    );
                    scalar(field.input_power_w, lo.input_power_w, hi.input_power_w);
                }
            }
            scalar(
                result.mean_stored_energy_j,
                lo.mean_stored_energy_j,
                hi.mean_stored_energy_j,
            );
            scalar(
                result.authored_external_power_w,
                lo.authored_external_power_w,
                hi.authored_external_power_w,
            );
            scalar(
                result.terminal_power_w,
                lo.terminal_power_w,
                hi.terminal_power_w,
            );
            scalar(
                result.total_dissipation_derivative_w(),
                lo.total_dissipation_w(),
                hi.total_dissipation_w(),
            );
            let scale = result.authored_external_power_w.abs()
                + result.terminal_power_w.abs()
                + result.total_dissipation_derivative_w().abs();
            assert!(result.power_residual_w.abs() <= 3e-8 * scale + 1e-27);
            assert!(result.interface_power_residual_w.abs() <= 3e-8 * scale + 1e-27);
            let closure = result
                .channels
                .iter()
                .find(|channel| channel.port == OUTLET)
                .unwrap();
            assert!(closure.channel.is_none());
            assert!(closure.exterior_pressure_pa.abs() > 0.0);
            assert!(closure.outward_flow_m3_s.abs() < 1e-9 * result.inlet_flow_m3_s.abs() + 1e-25);
            if let Some(old) = prior {
                close(result.inlet_pressure_pa, old);
            }
            prior = Some(result.inlet_pressure_pa);
            assert_eq!(load.state(), &state);
            assert_eq!(load.branches()[0].state(), &branch_state);
            assert_eq!(load.terminal().unwrap().pole_states(), terminal_state);
            assert_eq!(
                load.terminal_channels()[0]
                    .channel()
                    .unwrap()
                    .velocities_m_s(),
                channel_state
            );
        }
    }
    assert_eq!(system.material_memory(), memories);
    assert_eq!(system.radiation_state(), air);
    assert_eq!(system.samples_rendered(), 0);
}

#[test]
fn bare_duct_resonance_keeps_a_finite_joint_frequency_tangent() {
    let gate = CancelGate::new_clock_free();
    let duct = bare();
    let view = duct.continuous_operator();
    let l = view.generator_per_s();
    assert_eq!(view.dimension(), 4);
    let g = [l[2], l[3], l[6], l[7]];
    let (a, d, c) = (
        g[0] * g[0] + g[1] * g[1],
        g[2] * g[2] + g[3] * g[3],
        g[0] * g[2] + g[1] * g[3],
    );
    let omega = det::sqrt(0.5 * (a + d - det::sqrt((a - d) * (a - d) + 4.0 * c * c)));
    let system = oscillator();
    let force = C64::new(1e-5, -4e-6);
    let result = system
        .harmonic_response_with_duct_frequency_derivative(
            omega,
            &[force],
            None,
            &duct,
            INLET,
            &[0.02],
            WORK,
            &gate,
        )
        .unwrap();
    let lo = system
        .harmonic_response_with_duct(omega - STEP, &[force], &duct, INLET, &[0.02], WORK, &gate)
        .unwrap();
    let hi = system
        .harmonic_response_with_duct(omega + STEP, &[force], &duct, INLET, &[0.02], WORK, &gate)
        .unwrap();
    vector(
        &result.modal.displacement_derivative,
        &lo.modal.displacement_m_sqrt_kg,
        &hi.modal.displacement_m_sqrt_kg,
    );
    vector(
        &result.duct.energy_coordinates,
        &lo.duct.energy_coordinates,
        &hi.duct.energy_coordinates,
    );
    assert!(result.modal.displacement_derivative[0].abs() > 0.0);
    assert!(result.solve_relative_residual <= 1e-11);
}

#[test]
fn loaded_zero_force_slope_and_refusals_keep_all_original_histories() {
    let gate = CancelGate::new_clock_free();
    let system = oscillator();
    let load = loaded(false, false);
    let original = load.state().clone();
    let slope = C64::new(1e-5, -4e-6);
    let result = system
        .harmonic_response_with_load_frequency_derivative(
            1300.0,
            &[C64::ZERO],
            Some(&[slope]),
            &load,
            &[0.003],
            WORK,
            &gate,
        )
        .unwrap();
    let linear = system
        .harmonic_response_with_load(1300.0, &[slope], &load, &[0.003], WORK, &gate)
        .unwrap();
    assert_eq!(
        result.modal.displacement_derivative,
        linear.modal.displacement_m_sqrt_kg
    );
    assert_eq!(
        result.duct.energy_coordinates,
        linear.duct.energy_coordinates
    );
    for value in [
        result.mean_stored_energy_j,
        result.authored_external_power_w,
        result.terminal_power_w,
        result.total_dissipation_derivative_w(),
    ] {
        assert_eq!(value, 0.0);
    }
    assert!(
        system
            .harmonic_response_with_load_frequency_derivative(
                1300.0,
                &[slope],
                None,
                &load,
                &[0.003],
                1,
                &gate
            )
            .is_err()
    );
    assert!(
        system
            .harmonic_response_with_load_frequency_derivative(
                1300.0,
                &[slope],
                Some(&[]),
                &load,
                &[0.003],
                WORK,
                &gate
            )
            .is_err()
    );
    assert!(
        system
            .harmonic_response_with_load_frequency_derivative(
                1300.0,
                &[slope],
                Some(&[C64::new(f64::NAN, 0.0)]),
                &load,
                &[0.003],
                WORK,
                &gate
            )
            .is_err()
    );
    assert!(
        system
            .harmonic_response_with_load_frequency_derivative(
                1e9,
                &[slope],
                None,
                &load,
                &[0.003],
                WORK,
                &gate
            )
            .is_err()
    );
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        system.harmonic_response_with_load_frequency_derivative(
            1300.0,
            &[slope],
            None,
            &load,
            &[0.003],
            WORK,
            &cancelled
        ),
        Err(DgHarmonicCouplingError::Modal(
            ModalCouplingError::Cancelled
        ))
    ));
    assert_eq!(load.state(), &original);
    assert_eq!(system.samples_rendered(), 0);
}

fn resonance_controls() -> super::super::resonance::ResonanceControls {
    super::super::resonance::ResonanceControls {
        absolute_frequency_tolerance_rad_s: 1e-4,
        relative_frequency_tolerance: 1e-10,
        maximum_evaluations: 256,
        measure_half_power_bandwidth: true,
    }
}

#[test]
fn loaded_driven_peaks_keep_source_memory_gaps_and_mutual_terminal_port_order() {
    use super::super::resonance::{DgHarmonicObservation as Observe, ModalHarmonicObservation};
    let gate = CancelGate::new_clock_free();
    let system = oscillator();
    let force = [C64::new(1e-5, -4e-6)];
    let column = [0.003];
    let memory = system.material_memory().to_vec();
    let air = system.radiation_state().to_vec();
    let interval = (1720.0, 1950.0);
    for circular in [false, true] {
        let mut original: Option<(f64, f64)> = None;
        for reverse in [false, true] {
            let load = loaded(circular, reverse);
            let state = load.state().clone();
            let branch_state = load.branches()[0].state().clone();
            let terminal_state = load.terminal().unwrap().pole_states().to_vec();
            let channel_state = load.terminal_channels()[0]
                .channel()
                .unwrap()
                .velocities_m_s()
                .to_vec();
            let motion = system
                .driven_resonance_with_load(
                    &force,
                    &load,
                    &column,
                    Observe::Modal(ModalHarmonicObservation::Displacement(&[1.0])),
                    interval,
                    resonance_controls(),
                    WORK,
                    &gate,
                )
                .unwrap();
            let pressure = system
                .driven_resonance_with_load(
                    &force,
                    &load,
                    &column,
                    Observe::TerminalPressure(SIDE),
                    interval,
                    resonance_controls(),
                    WORK,
                    &gate,
                )
                .unwrap();
            // A local acoustic observation includes its dispersive transfer.
            // Its true driven peak need not coincide with the body's motion.
            assert!(
                (motion.peak.angular_frequency_rad_s - pressure.peak.angular_frequency_rad_s).abs()
                    > 1e-4
            );
            for (result, pressure_observation) in [(&motion, false), (&pressure, true)] {
                let field = |response: &DgLoadedModalHarmonicResponse| {
                    if pressure_observation {
                        let terminal = response.terminal.as_ref().unwrap();
                        terminal.pressure_pa[terminal
                            .ports
                            .iter()
                            .position(|port| *port == SIDE)
                            .unwrap()]
                    } else {
                        response.modal.displacement_m_sqrt_kg[0]
                    }
                };
                let response = system
                    .harmonic_response_with_load(
                        result.peak.angular_frequency_rad_s,
                        &force,
                        &load,
                        &column,
                        WORK,
                        &gate,
                    )
                    .unwrap();
                close(result.observation, field(&response));
                assert!(
                    (result.observation_derivative / result.observation)
                        .re
                        .abs()
                        < 2e-5
                );
                let band = result.half_power.unwrap();
                for edge in [band.lower, band.upper] {
                    let response = system
                        .harmonic_response_with_load(
                            edge.angular_frequency_rad_s,
                            &force,
                            &load,
                            &column,
                            WORK,
                            &gate,
                        )
                        .unwrap();
                    let ratio = field(&response).norm_sq() / result.observation.norm_sq();
                    assert!((ratio - 0.5).abs() < 2e-5, "half-height ratio {ratio}");
                }
                assert!(band.bandwidth_rad_s > 0.0);
                assert!(result.maximum_relative_residual < 1e-9);
            }
            if let Some((old_motion, old_pressure)) = original {
                assert!((motion.peak.angular_frequency_rad_s - old_motion).abs() < 2e-4);
                assert!((pressure.peak.angular_frequency_rad_s - old_pressure).abs() < 2e-4);
            }
            original = Some((
                motion.peak.angular_frequency_rad_s,
                pressure.peak.angular_frequency_rad_s,
            ));
            assert_eq!(load.state(), &state);
            assert_eq!(load.branches()[0].state(), &branch_state);
            assert_eq!(load.terminal().unwrap().pole_states(), terminal_state);
            assert_eq!(
                load.terminal_channels()[0]
                    .channel()
                    .unwrap()
                    .velocities_m_s(),
                channel_state
            );
        }
    }
    assert_eq!(system.material_memory(), memory);
    assert_eq!(system.radiation_state(), air);
    assert_eq!(system.samples_rendered(), 0);
}

#[test]
fn local_loaded_observations_use_original_pressure_planes_and_flow_signs() {
    use super::super::resonance::{DgHarmonicObservation as Observe, ModalHarmonicObservation};
    let gate = CancelGate::new_clock_free();
    let system = oscillator();
    let load = loaded(false, true);
    let force = [C64::new(1e-5, -4e-6)];
    let column = [0.003];
    let mut controls = resonance_controls();
    controls.measure_half_power_bandwidth = false;
    for observation in [
        Observe::Modal(ModalHarmonicObservation::Velocity(&[1.0])),
        Observe::Modal(ModalHarmonicObservation::Acceleration(&[1.0])),
        Observe::BorePressure(INLET),
        Observe::BoreInflow(INLET),
        Observe::BranchPressure {
            branch: SIDE,
            port: OUTLET,
        },
        Observe::BranchInflow {
            branch: SIDE,
            port: OUTLET,
        },
        Observe::TerminalOutwardFlow(SIDE),
        Observe::ChannelPressureDrop(SIDE),
        Observe::ChannelFlow(SIDE),
        Observe::ChannelPressureDrop(OUTLET),
    ] {
        let peak = system
            .driven_resonance_with_load(
                &force,
                &load,
                &column,
                observation,
                (1720.0, 1950.0),
                controls,
                WORK,
                &gate,
            )
            .unwrap();
        assert!(peak.half_power.is_none());
        let response = system
            .harmonic_response_with_load(
                peak.peak.angular_frequency_rad_s,
                &force,
                &load,
                &column,
                WORK,
                &gate,
            )
            .unwrap();
        let expected = match observation {
            Observe::Modal(ModalHarmonicObservation::Velocity(_)) => {
                response.modal.velocity_m_sqrt_kg_per_s[0]
            }
            Observe::Modal(ModalHarmonicObservation::Acceleration(_)) => {
                response.modal.displacement_m_sqrt_kg[0]
                    .scale(-response.modal.omega_rad_s * response.modal.omega_rad_s)
            }
            Observe::BorePressure(_) => response.inlet_pressure_pa,
            Observe::BoreInflow(_) => response.inlet_flow_m3_s,
            Observe::BranchPressure { .. } => {
                response
                    .channels
                    .iter()
                    .find(|channel| channel.port == SIDE)
                    .unwrap()
                    .duct_pressure_pa
            }
            Observe::BranchInflow { .. } => {
                -response
                    .channels
                    .iter()
                    .find(|channel| channel.port == SIDE)
                    .unwrap()
                    .outward_flow_m3_s
            }
            Observe::TerminalOutwardFlow(_) | Observe::ChannelFlow(_) => {
                response
                    .channels
                    .iter()
                    .find(|channel| channel.port == SIDE)
                    .unwrap()
                    .outward_flow_m3_s
            }
            Observe::ChannelPressureDrop(port) => {
                let channel = response
                    .channels
                    .iter()
                    .find(|channel| channel.port == port)
                    .unwrap();
                channel.duct_pressure_pa - channel.exterior_pressure_pa
            }
            _ => unreachable!(),
        };
        close(peak.observation, expected);
        assert!((peak.observation_derivative / peak.observation).re.abs() < 2e-5);
    }
    assert_eq!(load.state().step_count(), 0);
    assert_eq!(system.samples_rendered(), 0);
}

#[test]
fn bare_loaded_peak_and_native_source_refusals_preserve_all_histories() {
    use super::super::resonance::{
        DgHarmonicObservation as Observe, ModalHarmonicObservation, ResonanceSearchError,
    };
    let gate = CancelGate::new_clock_free();
    let system = oscillator();
    let duct = bare();
    let force = [C64::new(1e-5, -4e-6)];
    let column = [0.003];
    let observation = Observe::Modal(ModalHarmonicObservation::Velocity(&[1.0]));
    let result = system
        .driven_resonance_with_duct(
            &force,
            &duct,
            INLET,
            &column,
            observation,
            (1600.0, 1950.0),
            resonance_controls(),
            WORK,
            &gate,
        )
        .unwrap();
    assert!(result.half_power.is_some());
    let original = system
        .harmonic_response_with_duct(
            result.peak.angular_frequency_rad_s,
            &force,
            &duct,
            INLET,
            &column,
            WORK,
            &gate,
        )
        .unwrap();
    close(
        result.observation,
        original.modal.velocity_m_sqrt_kg_per_s[0],
    );
    let load = loaded(false, false);
    let state = load.state().clone();
    let memory = system.material_memory().to_vec();
    let air = system.radiation_state().to_vec();
    // Terminal/branch band violations and unknown physical locations retain
    // native errors. A known exact closure flow must not become a noise peak.
    let outside_branch_band =
        1.01 * load.branches()[0].duct().viscothermal_band_rad_s().unwrap()[1];
    for (observation, interval, work) in [
        (observation, (80.0, 1950.0), WORK),
        (observation, (1720.0, outside_branch_band), WORK),
        (observation, (1720.0, 1950.0), 1),
        (Observe::BorePressure(PortId(999)), (1720.0, 1950.0), WORK),
        (
            Observe::TerminalPressure(PortId(999)),
            (1720.0, 1950.0),
            WORK,
        ),
        (
            Observe::BranchPressure {
                branch: SIDE,
                port: PortId(999),
            },
            (1720.0, 1950.0),
            WORK,
        ),
        (Observe::BoreInflow(OUTLET), (1720.0, 1950.0), WORK),
        (Observe::TerminalOutwardFlow(OUTLET), (1720.0, 1950.0), WORK),
        (Observe::ChannelFlow(OUTLET), (1720.0, 1950.0), WORK),
    ] {
        let result = system.driven_resonance_with_load(
            &force,
            &load,
            &column,
            observation,
            interval,
            resonance_controls(),
            work,
            &gate,
        );
        assert!(
            matches!(&result, Err(ResonanceSearchError::Source(_))),
            "{observation:?} {interval:?} work={work}: {result:?}"
        );
    }
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        system.driven_resonance_with_load(
            &force,
            &load,
            &column,
            observation,
            (1720.0, 1950.0),
            resonance_controls(),
            WORK,
            &cancelled,
        ),
        Err(ResonanceSearchError::Source(
            DgHarmonicCouplingError::Modal(ModalCouplingError::Cancelled)
        )) | Err(ResonanceSearchError::Cancelled)
    ));
    assert_eq!(load.state(), &state);
    assert_eq!(system.material_memory(), memory);
    assert_eq!(system.radiation_state(), air);
    assert_eq!(system.samples_rendered(), 0);
}
