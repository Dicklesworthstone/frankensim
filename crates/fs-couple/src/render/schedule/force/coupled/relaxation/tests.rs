use super::*;
use crate::modal_acoustic_time::{ModalAcousticMode, ModalAcousticTimeBudget};
use super::super::contact::friction::{
    FrictionDrive, FrictionModalSystem, ModalFrictionConfig, ModalFrictionPort,
};
use fs_math::c64::C64;
use fs_ports::PortedModule;
use fs_tribo::FrictionLaw;

mod tensor;

fn config() -> ModalCouplingConfig {
    ModalCouplingConfig {
        max_modes: 16, max_connections: 8, max_setup_terms: 1_000_000,
        nyquist_guard_fraction: 0.9, maximum_total_energy_j: 100.0,
        maximum_abs_pressure_pa: 1e6, maximum_abs_connection_force_n: 1e6,
        solve_relative_tolerance: 1e-11, energy_absolute_tolerance_j: 1e-12,
        energy_relative_tolerance: 1e-9,
    }
}

fn model(rate: u32, omega: f64, damping: f64, q: f64, v: f64) -> ModalAcousticTimeModel {
    let mut model = ModalAcousticTimeModel::try_new_with_free_coordinates(rate,
        vec![ModalAcousticMode { angular_frequency_rad_s: omega, damping_ratio: damping,
            pressure_per_modal_velocity: C64::ZERO }], ModalAcousticTimeBudget::audible_reference()).unwrap();
    model.restore_states(&[ModalAcousticState {
        displacement_m_sqrt_kg: q, velocity_m_sqrt_kg_per_s: v,
    }]).unwrap();
    model
}

fn arm(projection: Vec<f64>, stiffness: f64, tau: f64) -> RelaxationBranch {
    RelaxationBranch { projection, stiffness, relaxation_time_s: tau }
}

fn mechanical(rate: u32, omega: f64, q: f64) -> CoupledModalSystem {
    CoupledModalSystem::new(vec![model(rate, omega, 0.0, q, 0.0)], vec![], config(),
        &CancelGate::new_clock_free()).unwrap()
}

fn close(a: f64, b: f64) {
    assert!((a-b).abs() <= 2e-9 * a.abs().max(b.abs()) + 2e-13, "{a:e} != {b:e}");
}

fn states(system: &CoupledModalSystem) -> Vec<ModalAcousticState> {
    system.components().iter().flat_map(|m| m.states().iter().copied()).collect()
}

#[test]
fn first_order_memory_joins_disconnected_masses_and_matches_analytic_maxwell_motion() {
    let gate = CancelGate::new_clock_free();
    let (c, tau, q0) = (1000.0_f64, 0.001_f64, 1e-4_f64);
    let mut system = CoupledModalSystem::new(vec![
        model(48_000, 0.0, 0.0, q0/2.0, 0.0),
        model(48_000, 0.0, 0.0, -q0/2.0, 0.0),
    ], vec![], config(), &gate).unwrap().with_material_relaxation(
        vec![arm(vec![1.0, 0.0, -1.0, 0.0], c*c, tau)], Some(vec![0.0]), &gate,
    ).unwrap();
    assert_eq!(system.mode_count(), 2); // z is not an invented third mass.
    assert_eq!(system.material_memory().len(), 1);
    let response = system.linear_response(&[1.0, 0.0], &gate).unwrap();
    assert!(response[1] > 0.0, "shared material must react onto the other mass");
    let initial = 0.5 * (c*q0).powi(2);
    let gamma = 0.5 / tau;
    let omega = (2.0*c*c - gamma*gamma).sqrt();
    let mut loss = 0.0;
    for tick in 1..=16 {
        let frame = system.step(&[0.0, 0.0]).unwrap().clone();
        let t = f64::from(tick) / 48_000.0;
        let relative_q = q0 * (-gamma*t).exp()
            * ((omega*t).cos() + gamma/omega * (omega*t).sin());
        let relative_v = -q0 * 2.0*c*c/omega * (-gamma*t).exp() * (omega*t).sin();
        let z = -relative_v/(2.0*c*tau);
        let pair = states(&system);
        close(pair[0].displacement_m_sqrt_kg, relative_q/2.0);
        close(pair[1].displacement_m_sqrt_kg, -relative_q/2.0);
        close(pair[0].velocity_m_sqrt_kg_per_s, relative_v/2.0);
        close(pair[1].velocity_m_sqrt_kg_per_s, -relative_v/2.0);
        close(system.material_memory()[0], z);
        close(frame.material_energy_j, 0.5*(c*relative_q-z).powi(2));
        let work = system.material_relaxation_observation_work(0..1).unwrap();
        let observed = system.material_relaxation_observation(0..1, work, &gate).unwrap();
        assert_eq!(observed.branch_count, 1);
        assert_eq!(observed.stored_energy_j, system.material_energy_j().unwrap());
        close(observed.stored_energy_j, frame.material_energy_j);
        close(observed.dissipated_power_w, (c*relative_q-z).powi(2)/tau);
        assert_eq!(frame.component_dissipation_j, 0.0);
        assert_eq!(frame.connection_dissipation_j, 0.0);
        assert!(frame.material_dissipation_j > 0.0);
        loss += frame.material_dissipation_j;
        close(system.total_energy_j().unwrap() + loss, initial);
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
    }
}

#[test]
fn stiff_coupled_memory_matches_overdamped_motion_and_integrated_heat_in_both_exact_owners() {
    let gate = CancelGate::new_clock_free();
    let (c, tau, q0) = (100.0_f64, 1e-7_f64, 1e-3_f64);
    let rate = 48_000;
    assert!(1.0 / (f64::from(rate) * tau) > 200.0);
    let gamma = 0.5 / tau;
    let root = (gamma * gamma - 2.0 * c * c).sqrt();
    // Product-preserving slow root avoids cancelling gamma - root. The
    // relative coordinate solves q'' + q'/tau + 2*c^2*q = 0, with z(0)=0.
    let fast = -gamma - root;
    let slow = 2.0 * c * c / fast;
    let separation = slow - fast;
    let initial = 0.5 * (c * q0).powi(2);
    for prepared in [false, true] {
        let network = CoupledModalSystem::new(vec![
            model(rate, 0.0, 0.0, q0 / 2.0, 0.0),
            model(rate, 0.0, 0.0, -q0 / 2.0, 0.0),
        ], vec![], config(), &gate).unwrap();
        let network = if prepared {
            network.with_prepared_linear_propagation(PreparedModalBudget {
                max_modes: 2, max_dense_bytes: 8 * 1024 * 1024,
                max_multiply_terms: 1_000_000_000,
            }, &gate).unwrap()
        } else { network };
        let mut system = network.with_material_relaxation(vec![
            arm(vec![1.0, 0.0, -1.0, 0.0], c * c, tau),
        ], Some(vec![0.0]), &gate).unwrap();
        let mut heat = 0.0;
        for tick in 1..=8 {
            let frame = system.step(&[0.0, 0.0]).unwrap().clone();
            let t = f64::from(tick) / f64::from(rate);
            let es = (slow * t).exp();
            let ef = (fast * t).exp();
            let q = q0 * (-fast * es + slow * ef) / separation;
            let v = q0 * 2.0 * c * c * (ef - es) / separation;
            let z = -v / (2.0 * c * tau);
            let pair = states(&system);
            close(pair[0].displacement_m_sqrt_kg, q / 2.0);
            close(pair[1].displacement_m_sqrt_kg, -q / 2.0);
            close(pair[0].velocity_m_sqrt_kg_per_s, v / 2.0);
            close(pair[1].velocity_m_sqrt_kg_per_s, -v / 2.0);
            close(system.material_memory()[0], z);
            // Integrate the independent analytic mismatch squared / tau,
            // rather than infer heat from the owner's energy residual.
            let a = c * q0 * slow / separation;
            let b = -c * q0 * fast / separation;
            let integral = |pole: f64| (pole * t).exp_m1() / pole;
            let expected_heat = (a * a * integral(2.0 * slow)
                + 2.0 * a * b * integral(slow + fast)
                + b * b * integral(2.0 * fast)) / tau;
            heat += frame.material_dissipation_j;
            close(heat, expected_heat);
            close(system.total_energy_j().unwrap() + heat, initial);
            assert!(frame.material_dissipation_j >= 0.0);
            assert_eq!(frame.component_dissipation_j, 0.0);
            assert_eq!(frame.connection_dissipation_j, 0.0);
            assert_eq!(frame.external_work_j, 0.0);
            assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        }
    }
}

fn coupled(rate: u32) -> CoupledModalSystem {
    let gate = CancelGate::new_clock_free();
    CoupledModalSystem::new_with_columns(vec![
        model(rate, 500.0, 0.01, 1e-4, 0.02),
        model(rate, 900.0, 0.04, -2e-4, -0.01),
    ], vec![ModalColumnConnection { column: vec![1.0, -0.3], stiffness_n_m: 4e4,
        damping_n_s_m: 20.0, rest_extension_m: 1e-5 }], config(), &gate).unwrap()
        .with_material_relaxation(vec![
            arm(vec![1.0, 0.0, -0.4, 0.0], 3e4, 0.002),
            arm(vec![0.6, 0.0, 1.0, 0.0], 2e4, 0.005),
        ], Some(vec![0.007, -0.004]), &gate).unwrap()
}

#[test]
fn material_observation_preserves_coherent_projection_and_disjoint_source_ranges() {
    let gate = CancelGate::new_clock_free();
    let q = 1.0 / 1024.0;
    let system = CoupledModalSystem::new(vec![
        model(48_000, 100.0, 0.0, q, 0.0),
        model(48_000, 100.0, 0.0, q, 0.0),
    ], vec![], config(), &gate).unwrap().with_material_relaxation(vec![
        arm(vec![1.0, 0.0, -1.0, 0.0], 4.0, 0.5),
        arm(vec![1.0, 0.0, 1.0, 0.0], 9.0, 0.25),
    ], Some(vec![0.0, q]), &gate).unwrap();
    let observe = |range: core::ops::Range<usize>| {
        let work = system.material_relaxation_observation_work(range.clone()).unwrap();
        system.material_relaxation_observation(range, work, &gate).unwrap()
    };
    let cancelled_projection = observe(0..1);
    assert_eq!(cancelled_projection.stored_energy_j, 0.0);
    assert_eq!(cancelled_projection.dissipated_power_w, 0.0);
    let constructive = observe(1..2);
    // The signed sum is formed before squaring, and supplied memory is retained.
    let mismatch = 5.0*q;
    assert_eq!(constructive.stored_energy_j, 0.5*mismatch*mismatch);
    assert_eq!(constructive.dissipated_power_w, mismatch*mismatch/0.25);
    let all = observe(0..2);
    assert_eq!(all.branch_count, 2);
    assert_eq!(all.stored_energy_j,
        cancelled_projection.stored_energy_j + constructive.stored_energy_j);
    assert_eq!(all.dissipated_power_w,
        cancelled_projection.dissipated_power_w + constructive.dissipated_power_w);
    assert_eq!(all.stored_energy_j, system.material_energy_j().unwrap());
    assert_eq!(observe(2..2), MaterialRelaxationObservation::default());
    assert_eq!(system.material_relaxation_observation_work(2..2).unwrap(), 0);
}

#[test]
fn material_observation_distinguishes_relaxed_and_unrelaxed_initial_history() {
    let gate = CancelGate::new_clock_free();
    let branch = arm(vec![1.0, 0.0], 4e4, 0.002);
    let relaxed = mechanical(48_000, 300.0, 1e-4)
        .with_material_relaxation(vec![branch.clone()], None, &gate).unwrap();
    let unrelaxed = mechanical(48_000, 300.0, 1e-4)
        .with_material_relaxation(vec![branch], Some(vec![0.0]), &gate).unwrap();
    let work = relaxed.material_relaxation_observation_work(0..1).unwrap();
    let settled = relaxed.material_relaxation_observation(0..1, work, &gate).unwrap();
    assert_eq!(settled.stored_energy_j, 0.0);
    assert_eq!(settled.dissipated_power_w, 0.0);
    let initial = unrelaxed.material_relaxation_observation(0..1, work, &gate).unwrap();
    close(initial.stored_energy_j, 0.0002);
    close(initial.dissipated_power_w, 0.2);
    assert_eq!(states(&relaxed), states(&unrelaxed));
    assert_eq!(relaxed.samples_rendered(), 0);
    assert_eq!(unrelaxed.samples_rendered(), 0);
}

#[test]
fn bounded_or_cancelled_observation_never_changes_the_accepted_state_or_next_step() {
    let gate = CancelGate::new_clock_free();
    let stopped = CancelGate::new_clock_free();
    stopped.request();
    let (mut observed, mut control) = (coupled(48_000), coupled(48_000));
    observed.step(&[0.02, -0.01]).unwrap();
    control.step(&[0.02, -0.01]).unwrap();
    let before = states(&observed);
    let memory = observed.material_memory().to_vec();
    let frame = observed.last_frame().unwrap().clone();
    let work = observed.material_relaxation_observation_work(0..2).unwrap();
    assert!(work > 0);
    let good = observed.material_relaxation_observation(0..2, work, &gate).unwrap();
    assert!(observed.material_relaxation_observation(0..2, work-1, &gate).is_err());
    assert!(matches!(observed.material_relaxation_observation(0..2, work, &stopped),
        Err(ModalCouplingError::Cancelled)));
    assert!(observed.material_relaxation_observation_work(
        core::ops::Range { start: 2, end: 1 }).is_err());
    assert!(observed.material_relaxation_observation(0..3, usize::MAX, &gate).is_err());
    assert_eq!(states(&observed), before);
    assert_eq!(observed.material_memory(), memory);
    assert_eq!(observed.last_frame(), Some(&frame));
    assert_eq!(good.stored_energy_j, observed.material_energy_j().unwrap());
    assert_eq!(observed.samples_rendered(), 1);
    let after = observed.step(&[0.02, -0.01]).unwrap().clone();
    assert_eq!(after, control.step(&[0.02, -0.01]).unwrap().clone());
    assert_eq!(states(&observed), states(&control));
    assert_eq!(observed.material_memory(), control.material_memory());
}

#[test]
fn history_work_and_three_independent_loss_channels_obey_time_partition() {
    let mut whole = coupled(48_000);
    let mut half = coupled(96_000);
    let initial = whole.total_energy_j().unwrap();
    let full = whole.step(&[0.02, -0.01]).unwrap().clone();
    let first = half.step(&[0.02, -0.01]).unwrap().clone();
    let first_work = half.last_window_transfers().unwrap()[0].energy_into().value();
    let second = half.step(&[0.02, -0.01]).unwrap().clone();
    for (a, b) in states(&whole).iter().zip(states(&half)) {
        close(a.displacement_m_sqrt_kg, b.displacement_m_sqrt_kg);
        close(a.velocity_m_sqrt_kg_per_s, b.velocity_m_sqrt_kg_per_s);
    }
    for (a, b) in whole.material_memory().iter().zip(half.material_memory()) { close(*a, *b); }
    close(full.external_work_j, first.external_work_j + second.external_work_j);
    close(full.component_dissipation_j, first.component_dissipation_j + second.component_dissipation_j);
    close(full.connection_dissipation_j, first.connection_dissipation_j + second.connection_dissipation_j);
    close(full.material_dissipation_j, first.material_dissipation_j + second.material_dissipation_j);
    assert!(full.component_dissipation_j > 0.0 && full.connection_dissipation_j > 0.0
        && full.material_dissipation_j > 0.0);
    close(whole.last_window_transfers().unwrap()[0].energy_into().value(),
        first_work + half.last_window_transfers().unwrap()[0].energy_into().value());
    close(full.total_energy_j() + full.total_dissipation_j() - full.external_work_j, initial);
    let fs_ports::ModuleWindowAudit::Available(audit) = whole.audit_last_window().unwrap()
        else { panic!("material account was not published"); };
    assert!(audit.violations().is_empty());
}

#[test]
fn relaxed_preload_preserves_equilibrium_but_supplied_history_changes_the_motion() {
    let gate = CancelGate::new_clock_free();
    let branch = arm(vec![1.0, 0.0], 4e4, 0.002);
    let mut settled = mechanical(48_000, 300.0, 0.0)
        .with_material_relaxation(vec![branch.clone()], None, &gate).unwrap();
    let energy = settled.initialize_static_equilibrium(&[9.0], &gate).unwrap();
    close(settled.components()[0].states()[0].displacement_m_sqrt_kg, 1e-4);
    close(settled.material_memory()[0], 0.02);
    close(energy, settled.total_energy_j().unwrap());
    assert_eq!(settled.material_energy_j().unwrap(), 0.0);
    let mut unrelaxed = mechanical(48_000, 300.0, 1e-4)
        .with_material_relaxation(vec![branch], Some(vec![0.0]), &gate).unwrap();
    for _ in 0..4 {
        settled.step(&[9.0]).unwrap();
        unrelaxed.step(&[9.0]).unwrap();
    }
    close(settled.components()[0].states()[0].displacement_m_sqrt_kg, 1e-4);
    assert!(unrelaxed.components()[0].states()[0].velocity_m_sqrt_kg_per_s < -1e-4);
    assert!(unrelaxed.material_memory()[0] > 0.0);
}

fn friction_config() -> ModalFrictionConfig {
    ModalFrictionConfig { max_iterations: 96, maximum_force_n: 1e4,
        maximum_slip_speed_m_s: 100.0, force_absolute_tolerance_n: 1e-10,
        force_relative_tolerance: 1e-9, sticking_tolerance_m: 1e-12 }
}

fn friction(network: CoupledModalSystem) -> FrictionModalSystem {
    FrictionModalSystem::new(network, ModalFrictionPort {
        left: ModalAttachment { component: 0, shapes: vec![1.0] }, right: None,
        law: FrictionLaw::Coulomb { static_mu: 0.8, kinetic_mu: 0.4 },
    }, friction_config(), &CancelGate::new_clock_free()).unwrap()
}

#[test]
fn friction_and_acceleration_pressure_see_the_same_material_reaction() {
    let gate = CancelGate::new_clock_free();
    for prepared in [false, true] {
        let network = mechanical(48_000, 100.0, 1e-4);
        let network = if prepared { network.with_prepared_linear_propagation(PreparedModalBudget {
            max_modes: 1, max_dense_bytes: 8 * 1024 * 1024, max_multiply_terms: 1_000_000_000,
        }, &gate).unwrap() } else { network };
        let mut owner = friction(network)
            .with_acceleration_pressure(vec![2.0]).unwrap()
            .with_material_relaxation(vec![arm(vec![1.0, 0.0], 4e6, 0.001)],
                Some(vec![0.0]), &gate).unwrap();
        for _ in 0..8 {
            let frame = *owner.step(&[0.03], FrictionDrive { speed_m_s: 0.04, normal_force_n: 0.2 }).unwrap();
            let state = owner.components()[0].states()[0];
            let material_force = -2000.0*(2000.0*state.displacement_m_sqrt_kg-owner.material_memory()[0]);
            let acceleration = 0.03 + frame.traction_n + material_force - 1e4*state.displacement_m_sqrt_kg;
            close(frame.observer_pressure_pa, 2.0*acceleration);
            close(frame.network_energy_j, owner.total_energy_j().unwrap());
            assert!(owner.coupled_frame().unwrap().material_dissipation_j > 0.0);
            assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        }
    }
}

#[test]
fn refused_pressure_and_cancelled_trials_retain_memory_and_retry_identically() {
    for prepared in [false, true] {
        let make = || {
            let mut network = mechanical(48_000, 100.0, 0.0);
            network.config.maximum_abs_pressure_pa = 0.001;
            let network = if prepared { network.with_prepared_linear_propagation(PreparedModalBudget {
                max_modes: 1, max_dense_bytes: 8 * 1024 * 1024, max_multiply_terms: 1_000_000_000,
            }, &CancelGate::new_clock_free()).unwrap() } else { network };
            friction(network).with_material_relaxation(vec![arm(vec![1.0, 0.0], 4e4, 0.002)],
                None, &CancelGate::new_clock_free()).unwrap()
                .with_acceleration_pressure(vec![2.0]).unwrap()
        };
        let (mut owner, mut clean) = (make(), make());
        let released = FrictionDrive { speed_m_s: 0.0, normal_force_n: 0.0 };
        owner.step(&[0.0], released).unwrap();
        clean.step(&[0.0], released).unwrap();
        let before = owner.components()[0].states().to_vec();
        let memory = owner.material_memory().to_vec();
        let frame = *owner.last_frame().unwrap();
        assert!(owner.step(&[10.0], released).is_err()); // fails after the z candidate exists
        assert_eq!(owner.components()[0].states(), before);
        assert_eq!(owner.material_memory(), memory);
        assert_eq!(owner.last_frame(), Some(&frame));
        assert_eq!(owner.samples_rendered(), 1);
        let gate = CancelGate::new_clock_free();
        gate.request();
        assert!(owner.step_under_gate(&[1e-6], released, &gate).is_err());
        assert_eq!(owner.material_memory(), memory);
        let retry = *owner.step(&[1e-6], released).unwrap();
        let expected = *clean.step(&[1e-6], released).unwrap();
        assert_eq!(retry, expected);
        assert_eq!(owner.material_memory(), clean.material_memory());
        assert_eq!(owner.components()[0].states(), clean.components()[0].states());
    }
}

#[test]
fn component_extension_and_nonlinear_reactions_keep_original_material_history() {
    use super::super::contact::friction::nonlinear::{
        MappedModalPotentialBinding, ModalCoordinate, ModalPotential, NonlinearModalConfig,
    };
    struct Quartic;
    impl fs_phs::Storage for Quartic {
        fn hamiltonian(&self, q: &[f64]) -> f64 { 0.25e11*q[0].powi(4) }
        fn gradient(&self, q: &[f64], out: &mut [f64]) { out[0] = 1e11*q[0].powi(3); }
    }
    impl ModalPotential for Quartic {
        fn coordinate_count(&self) -> usize { 1 }
        fn validate_coordinates(&self, q: &[f64]) -> Result<(), ModalCouplingError> {
            if q.len() != 1 || !q[0].is_finite() { return Err(invalid("quartic test coordinate")); }
            Ok(())
        }
    }
    let gate = CancelGate::new_clock_free();
    let owner = friction(mechanical(48_000, 100.0, 1e-4))
        .with_material_relaxation(vec![arm(vec![1.0, 0.0], 4e4, 0.002)],
            Some(vec![0.03]), &gate).unwrap()
        .with_internal_components(vec![model(48_000, 200.0, 0.0, 0.0, 0.0)], &gate).unwrap()
        .with_mapped_nonlinear_potentials(vec![MappedModalPotentialBinding {
            coordinates: vec![ModalCoordinate { component: 0, mode: 0 }],
            replaces_linear: vec![0], potential: Box::new(Quartic),
        }], NonlinearModalConfig { max_iterations: 64, maximum_generalized_force: 1e4,
            force_absolute_tolerance: 1e-9, force_relative_tolerance: 1e-8,
        }, &gate).unwrap();
    assert_eq!(owner.mode_count(), 2);
    assert_eq!(owner.material_memory(), [0.03]);
    assert_eq!(owner.material_relaxation_branches()[0].projection, [1.0, 0.0, 0.0, 0.0]);
    let mut owner = owner.with_additional_ports(vec![(ModalFrictionPort {
        left: ModalAttachment { component: 1, shapes: vec![1.0] }, right: None,
        law: FrictionLaw::Coulomb { static_mu: 0.8, kinetic_mu: 0.4 },
    }, friction_config())],
        super::super::contact::friction::multiple::MultiFrictionConfig {
            max_ports: 2, max_sweeps: 64, max_setup_terms: 10_000,
        }, &gate).unwrap();
    for _ in 0..4 {
        let frame = owner.step(&[0.0, 0.0], &[
            FrictionDrive { speed_m_s: 0.02, normal_force_n: 0.1 },
            FrictionDrive { speed_m_s: -0.01, normal_force_n: 0.05 },
        ]).unwrap();
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        assert!(frame.nonlinear.is_some());
    }
    assert_ne!(owner.material_memory(), [0.03]);
}

#[test]
fn no_memory_is_unchanged_and_augmented_capacity_setup_limits_remain_explicit() {
    let gate = CancelGate::new_clock_free();
    let branch = || arm(vec![1.0, 0.0], 100.0, 0.002);
    let mut plain = mechanical(48_000, 100.0, 1e-4);
    let mut empty = mechanical(48_000, 100.0, 1e-4)
        .with_material_relaxation(vec![], None, &gate).unwrap();
    assert!(plain.exact.is_none() && empty.exact.is_none());
    for _ in 0..4 {
        assert_eq!(plain.step(&[0.01]).unwrap(), empty.step(&[0.01]).unwrap());
        assert_eq!(states(&plain), states(&empty));
    }
    assert!(plain.with_material_relaxation(vec![branch()], None, &gate).is_err());
    let mut tight = mechanical(48_000, 100.0, 0.0);
    tight.config.max_setup_terms = 11; // 1 + (2+1)^2 + 2*1 requires twelve.
    assert!(tight.with_material_relaxation(vec![branch()], None, &gate).is_err());
    let oversized = vec![branch(); 125]; // 2*(2 mechanics + 125 memory + 1 input + 1)>256.
    assert!(matches!(mechanical(48_000, 100.0, 0.0)
        .with_material_relaxation(oversized, None, &gate),
        Err(ModalCouplingError::Propagation(fs_phs::PhsError::LinearSplitBudget { .. }))));
    let budget = PreparedModalBudget { max_modes: 1, max_dense_bytes: 8*1024*1024,
        max_multiply_terms: 1_000_000_000 };
    let material = mechanical(48_000, 100.0, 0.0)
        .with_material_relaxation(vec![branch()], None, &gate).unwrap();
    let material = material.with_prepared_linear_propagation(budget, &gate).unwrap();
    assert_eq!(material.material_memory().len(), 1);
    let prepared = mechanical(48_000, 100.0, 0.0)
        .with_prepared_linear_propagation(budget, &gate).unwrap();
    let prepared = prepared.with_material_relaxation(vec![branch()], None, &gate).unwrap();
    assert_eq!(prepared.material_memory(), material.material_memory());
    for bad in [arm(vec![1.0, 1.0], 1.0, 0.001),
        arm(vec![1.0, 0.0], 1.0, 0.0), arm(vec![1.0, 0.0], 1e12, 0.001)] {
        assert!(mechanical(48_000, 100.0, 0.0).with_material_relaxation(vec![bad], None, &gate).is_err());
    }
    assert!(mechanical(48_000, 100.0, 0.0)
        .with_material_relaxation(vec![branch()], Some(vec![f64::NAN]), &gate).is_err());
}
