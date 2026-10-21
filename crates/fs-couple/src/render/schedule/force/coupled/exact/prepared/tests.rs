use super::*;
use crate::modal_acoustic_time::{ModalAcousticMode, ModalAcousticTimeBudget};
use fs_math::c64::C64;
use fs_ports::PortedModule;

fn config(modes: usize, links: usize) -> ModalCouplingConfig {
    ModalCouplingConfig {
        max_modes: modes, max_connections: links,
        max_setup_terms: modes * (links + 1) * (links + 1),
        nyquist_guard_fraction: 0.9, maximum_total_energy_j: 100.0,
        maximum_abs_pressure_pa: 1e6, maximum_abs_connection_force_n: 1e8,
        solve_relative_tolerance: 1e-11, energy_absolute_tolerance_j: 1e-12,
        energy_relative_tolerance: 1e-10,
    }
}

fn budget(modes: usize) -> PreparedModalBudget {
    PreparedModalBudget {
        max_modes: modes, max_dense_bytes: 128 * 1024 * 1024,
        max_multiply_terms: 100_000_000_000,
    }
}

fn fixture(rate: u32, prepared: bool) -> CoupledModalSystem {
    // Noncommuting spring/damper columns, intrinsic loss, a free coordinate,
    // nonzero prestress, and independent held input on all original coordinates.
    let mut model = ModalAcousticTimeModel::try_new_with_free_coordinates(rate,
        [(500.0, 0.01), (900.0, 0.03), (0.0, 0.0)].into_iter()
            .map(|(w, z)| ModalAcousticMode { angular_frequency_rad_s: w,
                damping_ratio: z, pressure_per_modal_velocity: C64::ZERO }).collect(),
        ModalAcousticTimeBudget::audible_reference()).unwrap();
    model.restore_states(&[
        ModalAcousticState { displacement_m_sqrt_kg: 1e-4, velocity_m_sqrt_kg_per_s: 0.02 },
        ModalAcousticState { displacement_m_sqrt_kg: -2e-4, velocity_m_sqrt_kg_per_s: -0.01 },
        ModalAcousticState { displacement_m_sqrt_kg: 3e-5, velocity_m_sqrt_kg_per_s: 0.005 },
    ]).unwrap();
    let links = vec![
        ModalColumnConnection { column: vec![1.0, -0.7, 0.2],
            stiffness_n_m: 2e4, damping_n_s_m: 40.0, rest_extension_m: 1e-5 },
        ModalColumnConnection { column: vec![-0.2, 0.5, 1.0],
            stiffness_n_m: 1e4, damping_n_s_m: 25.0, rest_extension_m: -1e-5 },
    ];
    let gate = CancelGate::new_clock_free();
    if prepared {
        CoupledModalSystem::new_with_columns_prepared(vec![model], links, config(3, 2), budget(3), &gate).unwrap()
    } else {
        CoupledModalSystem::new_with_columns(vec![model], links, config(3, 2), &gate).unwrap()
    }
}

fn close(a: f64, b: f64) {
    assert!((a - b).abs() <= 2e-10 * a.abs().max(b.abs()) + 2e-15, "{a:e} != {b:e}");
}

#[test]
fn prepared_keeps_full_noncommuting_operator_rest_load_and_each_port() {
    let mut original = fixture(48_000, false);
    let mut prepared = fixture(48_000, true);
    let before_switch = flatten(&original.models);
    let mut switched = fixture(48_000, false)
        .with_prepared_linear_propagation(budget(3), &CancelGate::new_clock_free()).unwrap();
    assert_eq!(flatten(&switched.models), before_switch);
    assert_eq!(original.prepared_budget(), None);
    assert_eq!(prepared.prepared_budget(), Some(budget(3)));
    assert_eq!(switched.prepared_budget(), Some(budget(3)));
    let column = [0.3, -0.7, 0.1];
    for (a, b) in original.linear_response(&column, &CancelGate::new_clock_free()).unwrap()
        .iter().zip(prepared.linear_response(&column, &CancelGate::new_clock_free()).unwrap()) {
        close(*a, b);
    }
    for tick in 0..12 {
        let external = [0.003, if tick % 2 == 0 { -0.002 } else { 0.004 }, 0.001];
        let a = original.step(&external).unwrap().clone();
        let b = prepared.step(&external).unwrap().clone();
        let switched_frame = switched.step(&external).unwrap().clone();
        assert_eq!(switched_frame, b);
        assert_eq!(flatten(&switched.models), flatten(&prepared.models));
        for (left, right) in original.components()[0].states().iter().zip(prepared.components()[0].states()) {
            close(left.displacement_m_sqrt_kg, right.displacement_m_sqrt_kg);
            close(left.velocity_m_sqrt_kg_per_s, right.velocity_m_sqrt_kg_per_s);
        }
        close(a.external_work_j, b.external_work_j);
        close(a.component_dissipation_j, b.component_dissipation_j);
        close(a.connection_dissipation_j, b.connection_dissipation_j);
        assert!(b.component_dissipation_j > 0.0 && b.connection_dissipation_j > 0.0);
        for (left, right) in original.last_window_transfers().unwrap().iter()
            .zip(prepared.last_window_transfers().unwrap()) {
            close(left.energy_into().value(), right.energy_into().value());
        }
        assert!(prepared.last_window_transfer_error().is_none());
        let fs_ports::ModuleWindowAudit::Available(audit) = prepared.audit_last_window().unwrap() else {
            panic!("prepared exact physical account unavailable")
        };
        assert!(audit.violations().is_empty());
    }
}

#[test]
fn exact_loss_and_connection_work_obey_held_force_time_partition() {
    let mut full = fixture(48_000, true);
    let mut halves = fixture(96_000, true);
    let external = [0.003, -0.002, 0.001];
    let reference = full.step(&external).unwrap().clone();
    let first = halves.step(&external).unwrap().clone();
    let ports_first: Vec<_> = halves.last_window_transfers().unwrap().iter()
        .map(|v| v.energy_into().value()).collect();
    let second = halves.step(&external).unwrap().clone();
    for (a, b) in full.components()[0].states().iter().zip(halves.components()[0].states()) {
        close(a.displacement_m_sqrt_kg, b.displacement_m_sqrt_kg);
        close(a.velocity_m_sqrt_kg_per_s, b.velocity_m_sqrt_kg_per_s);
    }
    close(reference.external_work_j, first.external_work_j + second.external_work_j);
    close(reference.component_dissipation_j, first.component_dissipation_j + second.component_dissipation_j);
    close(reference.connection_dissipation_j, first.connection_dissipation_j + second.connection_dissipation_j);
    for ((a, b), c) in full.last_window_transfers().unwrap().iter()
        .zip(halves.last_window_transfers().unwrap()).zip(ports_first) {
        close(a.energy_into().value(), b.energy_into().value() + c);
    }
}

#[test]
fn larger_connected_network_retains_all_states_and_matches_analytic_common_motion() {
    // Forty-three identical masses share ONE physical viscous coordinate.
    // Its 43-component signed column creates one genuine 43-mode K/C island,
    // exceeding the ordinary 42-mode cap; no state or connection is discarded.
    let n = 43;
    let omega = 200.0_f64;
    let q0 = 1e-4_f64;
    let v0 = 0.02_f64;
    let mut model = ModalAcousticTimeModel::try_new(48_000,
        vec![ModalAcousticMode { angular_frequency_rad_s: omega,
            damping_ratio: 0.0, pressure_per_modal_velocity: C64::ZERO }; n],
        ModalAcousticTimeBudget::audible_reference()).unwrap();
    model.restore_states(&vec![ModalAcousticState {
        displacement_m_sqrt_kg: q0, velocity_m_sqrt_kg_per_s: v0 }; n]).unwrap();
    let links = vec![ModalColumnConnection { column: vec![1.0; n],
        stiffness_n_m: 0.0, damping_n_s_m: 2.0, rest_extension_m: 0.0 }];
    let gate = CancelGate::new_clock_free();
    assert!(CoupledModalSystem::new_with_columns(vec![model.clone()], links.clone(), config(n, 1), &gate).is_err());
    let mut system = CoupledModalSystem::new_with_columns_prepared(vec![model], links, config(n, 1), budget(n), &gate).unwrap();
    assert_eq!(system.mode_count(), n);
    let initial = 0.5 * n as f64 * (omega.powi(2) * q0.powi(2) + v0.powi(2));
    let mut loss = 0.0;
    let gamma = n as f64;
    let wd = (omega * omega - gamma * gamma).sqrt();
    for step in 1..=4 {
        let frame = system.step(&vec![0.0; n]).unwrap().clone();
        let t = f64::from(step) / 48_000.0;
        let q = (-gamma * t).exp() * (q0 * (wd * t).cos() + (v0 + gamma * q0) / wd * (wd * t).sin());
        let v = (-gamma * t).exp() * (v0 * (wd * t).cos() - (gamma * v0 + omega * omega * q0) / wd * (wd * t).sin());
        for state in system.components()[0].states() {
            close(state.displacement_m_sqrt_kg, q);
            close(state.velocity_m_sqrt_kg_per_s, v);
        }
        assert_eq!(system.last_window_transfers().unwrap().len(), 1);
        close(-system.last_window_transfers().unwrap()[0].energy_into().value(), frame.connection_dissipation_j);
        assert_eq!(frame.component_dissipation_j, 0.0);
        loss += frame.connection_dissipation_j;
    }
    close(initial, system.total_energy_j().unwrap() + loss);
}

#[test]
fn stiff_free_coordinate_memory_uses_short_loss_integrals_without_inverse_overflow() {
    // The full-interval inverse Van Loan block would contain exp(1000).
    // Positive damping is still integrated over all of the physical timestep.
    let model = ModalAcousticTimeModel::try_new_with_free_coordinates(48_000,
        vec![ModalAcousticMode { angular_frequency_rad_s: 0.0,
            damping_ratio: 0.0, pressure_per_modal_velocity: C64::ZERO }],
        ModalAcousticTimeBudget::audible_reference()).unwrap();
    let rate = 48_000_000.0_f64;
    let mut model = model;
    model.restore_states(&[ModalAcousticState {
        displacement_m_sqrt_kg: 0.0, velocity_m_sqrt_kg_per_s: 0.01 }]).unwrap();
    let mut system = CoupledModalSystem::new_with_columns_prepared(vec![model],
        vec![ModalColumnConnection { column: vec![1.0], stiffness_n_m: 0.0,
            damping_n_s_m: rate, rest_extension_m: 0.0 }],
        config(1, 1), budget(1), &CancelGate::new_clock_free()).unwrap();
    let frame = system.step(&[0.0]).unwrap().clone();
    let state = system.components()[0].states()[0];
    close(state.displacement_m_sqrt_kg, 0.01 / rate);
    assert!(state.velocity_m_sqrt_kg_per_s.abs() < 1e-100);
    close(frame.connection_dissipation_j, 0.5 * 0.01_f64.powi(2));
    close(system.total_energy_j().unwrap() + frame.connection_dissipation_j, 0.5 * 0.01_f64.powi(2));
}

#[test]
fn explicit_preparation_refusals_and_negative_loss_gate_do_not_invent_accuracy() {
    let seed = fixture(48_000, false);
    let columns = vec![ModalColumnConnection { column: vec![1.0, -1.0, 0.0],
        stiffness_n_m: 100.0, damping_n_s_m: 1.0, rest_extension_m: 0.0 }];
    let gate = CancelGate::new_clock_free();
    for refused in [
        PreparedModalBudget { max_modes: 2, ..budget(3) },
        PreparedModalBudget { max_dense_bytes: 1, ..budget(3) },
        PreparedModalBudget { max_multiply_terms: 1, ..budget(3) },
        PreparedModalBudget { max_modes: MAX_PREPARED_MODAL_MODES + 1, ..budget(3) },
    ] {
        assert!(CoupledModalSystem::new_with_columns_prepared(seed.models.clone(),
            columns.clone(), config(3, 1), refused, &gate).is_err());
    }
    assert!(CoupledModalSystem::new_with_columns_prepared(seed.models.clone(),
        columns.clone(), ModalCouplingConfig { max_connections: 0, ..config(3, 1) }, budget(3), &gate).is_err());
    assert!(quadratic(&[-1.0], &[1.0]).is_err());
    let mut system = fixture(48_000, true);
    let before = flatten(&system.models);
    assert!(system.step(&[f64::NAN, 0.0, 0.0]).is_err());
    assert_eq!(flatten(&system.models), before);
    assert_eq!(system.samples_rendered(), 0);
}

#[test]
fn one_stage_costs_one_loss_sweep_and_the_shared_values_equal_a_direct_evaluation_bit_for_bit() {
    let mut system = fixture(48_000, true);
    let sweeps = |s: &CoupledModalSystem| s.exact.as_ref().unwrap().prepared.as_ref().unwrap().sweeps();
    let before = sweeps(&system);
    let steps = 6u64;
    for tick in 0..steps {
        let external = [0.003, if tick % 2 == 0 { -0.002 } else { 0.004 }, 0.001];
        system.step(&external).unwrap();
    }
    // `advance`, `losses` and `connection_work` of a stage evaluate the same vector: one sweep per step, not three.
    assert_eq!(sweeps(&system) - before, steps);
    let prepared = system.exact.as_ref().unwrap().prepared.as_ref().unwrap();
    let x = flatten(&system.models);
    for u in [[0.003, -0.002, 0.001], [0.0, 0.0, 0.0], [0.003, -0.002, 0.001]] {
        let z = prepared.state_input(&x, &u, prepared.dt).unwrap();
        let direct_component = quadratic(prepared.component_loss.as_ref().unwrap(), &z).unwrap();
        let mut direct_connection = 0.0;
        for w in prepared.connection_loss.iter().flatten() {
            direct_connection = finite(direct_connection + quadratic(w, &z).unwrap()).unwrap();
        }
        let (component, connection) = prepared.loss_sums(&z).unwrap();
        assert_eq!(component.to_bits(), direct_component.to_bits());
        assert_eq!(connection.to_bits(), direct_connection.to_bits());
        assert!(component > 0.0 && connection > 0.0);
        // A second request for the same vector is served from the sweep, bit for bit.
        let count = prepared.sweeps();
        let (again_component, again_connection) = prepared.loss_sums(&z).unwrap();
        assert_eq!(prepared.sweeps(), count);
        assert_eq!((again_component.to_bits(), again_connection.to_bits()), (component.to_bits(), connection.to_bits()));
    }
}
