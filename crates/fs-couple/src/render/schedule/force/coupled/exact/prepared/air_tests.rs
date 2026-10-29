//! Two-way air on the prepared owner (`air`): the held exchange against the dense exact augmented owner, its energy
//! identity, the air-loaded affine response the contacts use, and its composition with constant-loss links.
use super::*;
use crate::modal_acoustic_time::{ModalAcousticMode, ModalAcousticTimeBudget};
use crate::render::plate::impact::radiation::{Model, Pole};
use crate::render::schedule::force::coupled::ConstantLossConnection;
use fs_math::c64::C64;
use fs_ports::PortedModule;

const TAU: f64 = core::f64::consts::TAU;

fn config() -> ModalCouplingConfig {
    ModalCouplingConfig {
        max_modes: 8, max_connections: 4, max_setup_terms: 1_000_000,
        nyquist_guard_fraction: 0.9, maximum_total_energy_j: 100.0,
        maximum_abs_pressure_pa: 1e6, maximum_abs_connection_force_n: 1e8,
        solve_relative_tolerance: 1e-11, energy_absolute_tolerance_j: 1e-12,
        energy_relative_tolerance: 1e-10,
    }
}

fn budget(modes: usize) -> PreparedModalBudget {
    PreparedModalBudget { max_modes: modes, max_dense_bytes: 128 * 1024 * 1024, max_multiply_terms: 100_000_000_000 }
}

fn model(rate: u32, modes: &[(f64, f64)], states: &[(f64, f64)]) -> ModalAcousticTimeModel {
    let mut model = ModalAcousticTimeModel::try_new(rate,
        modes.iter().map(|&(w, z)| ModalAcousticMode { angular_frequency_rad_s: w, damping_ratio: z,
            pressure_per_modal_velocity: C64::ZERO }).collect(),
        ModalAcousticTimeBudget::audible_reference()).unwrap();
    model.restore_states(&states.iter().map(|&(q, v)| ModalAcousticState {
        displacement_m_sqrt_kg: q, velocity_m_sqrt_kg_per_s: v }).collect::<Vec<_>>()).unwrap();
    model
}

/// One mode and one pole near it; or two modes under a signed link and three poles with signed cross terms on both,
/// one pole far above the modes.
struct Fixture {
    modes: Vec<(f64, f64)>,
    states: Vec<(f64, f64)>,
    links: Vec<ModalColumnConnection>,
    air: Model,
}

fn fixture(mimo: bool) -> Fixture {
    if !mimo {
        return Fixture { modes: vec![(3000.0, 0.005)], states: vec![(1e-4, 0.0)], links: vec![],
            air: Model { ports: 1, poles: vec![Pole { omega: 3500.0, zeta: 0.2, coupling: vec![400.0] }] } };
    }
    Fixture {
        modes: vec![(TAU * 480.0, 0.004), (TAU * 1150.0, 0.006)],
        states: vec![(2e-5, 0.01), (-1e-5, 0.02)],
        links: vec![ModalColumnConnection { column: vec![1.0, -0.6], stiffness_n_m: 2e6, damping_n_s_m: 15.0,
            rest_extension_m: 0.0 }],
        air: Model { ports: 2, poles: vec![
            Pole { omega: TAU * 520.0, zeta: 0.3, coupling: vec![350.0, -120.0] },
            Pole { omega: TAU * 1300.0, zeta: 0.15, coupling: vec![-90.0, 260.0] },
            Pole { omega: TAU * 4000.0, zeta: 0.5, coupling: vec![60.0, 80.0] },
        ] },
    }
}

fn mechanics(rate: u32, mimo: bool, prepared: bool) -> CoupledModalSystem {
    let f = fixture(mimo);
    let gate = CancelGate::new_clock_free();
    let models = vec![model(rate, &f.modes, &f.states)];
    if prepared {
        CoupledModalSystem::new_with_columns_prepared(models, f.links, config(), budget(f.modes.len()), &gate).unwrap()
    } else {
        CoupledModalSystem::new_with_columns(models, f.links, config(), &gate).unwrap()
    }
}

fn loaded(rate: u32, mimo: bool, prepared: bool) -> CoupledModalSystem {
    let f = fixture(mimo);
    let sources: Vec<_> = (0..f.air.ports).collect();
    mechanics(rate, mimo, prepared)
        .with_radiation_load(&f.air, &sources, f.air.poles.len(), &CancelGate::new_clock_free()).unwrap()
}

/// The energy norm sqrt(sum (omega q)^2 + v^2) over the modes and the poles, with each one's own frequency.
fn energy_norm(system: &CoupledModalSystem, mechanical: &[f64], air: &[f64]) -> f64 {
    let modes = system.models[0].modes();
    let poles = &system.radiation_model().unwrap().poles;
    let mut sum = 0.0;
    for (mode, x) in modes.iter().zip(mechanical.chunks_exact(2)) {
        sum += (mode.angular_frequency_rad_s * x[0]).powi(2) + x[1] * x[1];
    }
    for (pole, x) in poles.iter().zip(air.chunks_exact(2)) {
        sum += (pole.omega * x[0]).powi(2) + x[1] * x[1];
    }
    det::sqrt(sum)
}

// The held exchange is the implicit midpoint rule on the interconnection: against the dense exact augmented owner its
// error over 4 ms falls by 4 per halving of dt, for one mode and one pole and for two linked modes under three poles
// with signed cross terms. An independent numpy model of the same step (scratch pred_order.py) predicted 1.22e-4 and
// 1.88e-4 at 48 kHz, ratio 4.000 to the third digit.
#[test]
fn the_held_air_exchange_converges_at_second_order_to_the_exact_augmented_owner() {
    for mimo in [false, true] {
        let mut previous: Option<f64> = None;
        for rate in [48_000_u32, 96_000, 192_000, 384_000] {
            let mut exact = loaded(rate, mimo, false);
            let mut held = loaded(rate, mimo, true);
            assert!(exact.prepared_budget().is_none() && held.prepared_budget().is_some());
            let external = vec![0.0; exact.mode_count()];
            for _ in 0..rate / 250 {
                exact.step(&external).unwrap();
                held.step(&external).unwrap();
            }
            let (a, b) = (flatten(&exact.models), flatten(&held.models));
            let mechanical: Vec<_> = a.iter().zip(&b).map(|(x, y)| y - x).collect();
            let air: Vec<_> = exact.radiation_state().iter().zip(held.radiation_state()).map(|(x, y)| y - x).collect();
            let error = energy_norm(&exact, &mechanical, &air) / energy_norm(&exact, &a, exact.radiation_state());
            let ratio = previous.map_or(f64::NAN, |p| p / error);
            eprintln!("{} {rate} Hz: relative error {error:.4e}, ratio {ratio:.4}", if mimo { "mimo" } else { "siso" });
            if let Some(p) = previous {
                assert!((p / error - 4.0).abs() < 0.05, "{rate} Hz: ratio {}", p / error);
            }
            previous = Some(error);
        }
        assert!(previous.unwrap() < 5e-6);
    }
}

// Every step closes mechanics, links and air together: storage change plus every loss equals the external work, to
// 1e-11 of the storage. Released, the network never gains energy, and the air carries part of the loss.
#[test]
fn every_air_loaded_step_closes_its_energy_and_a_free_decay_never_gains() {
    let mut system = loaded(48_000, true, true);
    let mut energy = system.total_energy_j().unwrap();
    for tick in 0..480 {
        let external = [0.002 * (tick % 3) as f64, if tick % 2 == 0 { -0.003 } else { 0.001 }];
        let frame = system.step(&external).unwrap().clone();
        let residual = frame.total_energy_j() + frame.total_dissipation_j() - frame.external_work_j - energy;
        assert!(residual.abs() <= 1e-11 * energy.max(frame.total_energy_j()), "tick {tick}: residual {residual:e}");
        assert!(frame.radiation_dissipation_j > 0.0);
        close(frame.total_energy_j(), system.total_energy_j().unwrap());
        energy = frame.total_energy_j();
    }
    let (start, mut radiated, mut lost) = (energy, 0.0, 0.0);
    for tick in 0..9600 {
        let frame = system.step(&[0.0, 0.0]).unwrap().clone();
        assert!(frame.total_energy_j() <= energy, "tick {tick}: {:e} > {energy:e}", frame.total_energy_j());
        let residual = frame.total_energy_j() + frame.total_dissipation_j() - energy;
        assert!(residual.abs() <= 1e-11 * energy, "tick {tick}: residual {residual:e}");
        radiated += frame.radiation_dissipation_j;
        lost += frame.total_dissipation_j();
        energy = frame.total_energy_j();
    }
    eprintln!("free decay over 0.2 s: {start:e} -> {energy:e} J, radiated share {:.4}", radiated / lost);
    assert!(energy < start && radiated > 0.0 && radiated < lost);
}

fn close(a: f64, b: f64) {
    assert!((a - b).abs() <= 1e-12 * a.abs().max(b.abs()) + 1e-18, "{a:e} != {b:e}");
}

// A zero coupling binds nothing, and the two attachment orders (prepare, then air; air on the dense owner, then
// prepare) build the same prepared owner: the same steps bit for bit.
#[test]
fn zero_air_changes_nothing_and_both_attachment_orders_step_identically() {
    let gate = CancelGate::new_clock_free();
    let mut bare = mechanics(48_000, true, true);
    let mut zero_air = Model { ports: 2, poles: fixture(true).air.poles };
    for pole in &mut zero_air.poles { pole.coupling = vec![0.0, 0.0]; }
    let mut zero = mechanics(48_000, true, true).with_radiation_load(&zero_air, &[0, 1], 3, &gate).unwrap();
    assert!(zero.radiation_model().is_none());
    let mut first = loaded(48_000, true, true);
    let mut later = loaded(48_000, true, false).with_prepared_linear_propagation(budget(2), &gate).unwrap();
    for tick in 0..64 {
        let external = [0.001 * tick as f64, -0.002];
        assert_eq!(bare.step(&external).unwrap().clone(), zero.step(&external).unwrap().clone());
        assert_eq!(first.step(&external).unwrap().clone(), later.step(&external).unwrap().clone());
        assert_eq!(flatten(&first.models), flatten(&later.models));
        assert_eq!(first.radiation_state(), later.radiation_state());
    }
}

#[test]
fn prepared_air_material_and_anchored_joint_preserve_history_in_every_cold_order() {
    let gate = CancelGate::new_clock_free();
    let branch = fs_phs::RelaxationBranch { projection: vec![1.0, 0.0, -0.6, 0.0],
        stiffness: 8e5, relaxation_time_s: 0.001 };
    let anchors = [ConstantLossConnection { connection: 0, reference_omega_rad_s: TAU * 200.0 }];
    let attach = |system: CoupledModalSystem| system.with_material_relaxation(
        vec![branch.clone()], Some(vec![0.012]), &gate).unwrap();
    let mut first = attach(loaded(48_000, true, true))
        .with_constant_loss_connections(&anchors, &gate).unwrap();
    let f = fixture(true);
    let mut last = attach(mechanics(48_000, true, false))
        .with_prepared_linear_propagation(budget(2), &gate).unwrap()
        .with_constant_loss_connections(&anchors, &gate).unwrap()
        .with_radiation_load(&f.air, &[0, 1], 3, &gate).unwrap();
    assert_eq!(first.material_memory(), [0.012]);
    assert_eq!(first.material_memory(), last.material_memory());
    let mut energy = first.total_energy_j().unwrap();
    for tick in 0..32 {
        let force = [0.001 * tick as f64, -0.002];
        let a = first.step(&force).unwrap().clone();
        let b = last.step(&force).unwrap().clone();
        assert_eq!(a, b);
        assert_eq!(first.material_memory(), last.material_memory());
        assert_eq!(first.radiation_state(), last.radiation_state());
        assert!(a.material_dissipation_j > 0.0 && a.radiation_dissipation_j > 0.0);
        assert!(a.component_dissipation_j > 0.0 && a.connection_dissipation_j > 0.0);
        close(a.total_energy_j() + a.total_dissipation_j() - a.external_work_j, energy);
        assert_eq!(first.last_window_transfers().unwrap(), last.last_window_transfers().unwrap());
        assert!(first.last_window_transfer_error().is_none());
        energy = a.total_energy_j();
    }
}

// The contacts' affine model: the fast free prediction is the full advance's end state bit for bit, and a step under an
// extra force f along a column ends at that prediction plus f times the air-loaded response, which differs from the
// response of the same mechanics without air.
#[test]
fn the_air_loaded_prediction_and_response_compose_to_the_step() {
    let gate = CancelGate::new_clock_free();
    let mut system = loaded(48_000, true, true);
    for _ in 0..40 {
        system.step(&[0.002, -0.001]).unwrap();
    }
    assert!(system.radiation_state().iter().all(|&x| x != 0.0));
    let external = [0.003, -0.001];
    let exact = system.exact.as_ref().unwrap();
    let (fast, x) = exact.predict(&system.models, &system.material_memory, &system.radiation_state, &external,
        system.dt, false).unwrap();
    let (full, start) = exact.predict(&system.models, &system.material_memory, &system.radiation_state, &external,
        system.dt, true).unwrap();
    assert_eq!(fast[..], full[..4]);
    assert_eq!(x[..], start[..4]);
    assert_eq!(start[4..], *system.radiation_state());
    let column = [0.7, -0.4];
    let response = system.linear_response(&column, &gate).unwrap();
    let unloaded = mechanics(48_000, true, true).linear_response(&column, &gate).unwrap();
    assert!(response.iter().zip(&unloaded).all(|(a, b)| a != b));
    let force = 0.05;
    let forced: Vec<_> = external.iter().zip(&column).map(|(e, c)| e + force * c).collect();
    system.step(&forced).unwrap();
    let after = flatten(&system.models);
    for k in 0..2 {
        let (moved, free, extra) = (after[2 * k] - x[2 * k], fast[2 * k] - x[2 * k], force * response[k]);
        let scale = moved.abs().max(free.abs()).max(extra.abs());
        assert!((moved - free - extra).abs() <= 1e-12 * scale, "mode {k}: {moved:e} != {free:e} + {extra:e}");
    }
}

// A constant-loss link and the air together: every step's energy closes, the per-port ledger sums to the connection
// loss with a clean audit, and anchoring before or after the air steps identically.
#[test]
fn a_constant_loss_link_and_the_air_keep_both_ledgers() {
    let gate = CancelGate::new_clock_free();
    let reference = TAU * 200.0;
    let modes = [(TAU * 180.0, 0.002), (TAU * 950.0, 0.004), (TAU * 3100.0, 0.003)];
    let states = [(2e-5, 0.01), (-1e-5, -0.02), (4e-6, 0.03)];
    let links = || vec![
        ModalColumnConnection { column: vec![1.0, -0.6, 0.3], stiffness_n_m: 4e6,
            damping_n_s_m: 0.0162 * 4e6 / reference, rest_extension_m: 1e-7 },
        ModalColumnConnection { column: vec![-0.2, 0.7, 0.9], stiffness_n_m: 2e7,
            damping_n_s_m: 30.0, rest_extension_m: -1e-7 },
    ];
    let air = Model { ports: 2, poles: vec![
        Pole { omega: TAU * 200.0, zeta: 0.4, coupling: vec![300.0, -150.0] },
        Pole { omega: TAU * 2900.0, zeta: 0.2, coupling: vec![120.0, 200.0] },
    ] };
    let anchors = [ConstantLossConnection { connection: 0, reference_omega_rad_s: reference }];
    let build = |anchored: bool| CoupledModalSystem::new_with_columns_inner(vec![model(48_000, &modes, &states)], links(),
        config(), Some(budget(3)), if anchored { &anchors[..] } else { &[] }, &gate).unwrap();
    let mut system = build(true).with_radiation_load(&air, &[0, 2], 2, &gate).unwrap();
    let mut later = build(false).with_radiation_load(&air, &[0, 2], 2, &gate).unwrap()
        .with_constant_loss_connections(&anchors, &gate).unwrap();
    let mut energy = system.total_energy_j().unwrap();
    for tick in 0..240 {
        let external = [0.002, if tick % 2 == 0 { -0.003 } else { 0.001 }, 0.0005];
        let frame = system.step(&external).unwrap().clone();
        assert_eq!(frame, *later.step(&external).unwrap());
        let residual = frame.total_energy_j() + frame.total_dissipation_j() - frame.external_work_j - energy;
        assert!(residual.abs() <= 1e-11 * energy.max(frame.total_energy_j()), "tick {tick}: residual {residual:e}");
        assert!(frame.radiation_dissipation_j > 0.0 && frame.connection_dissipation_j > 0.0);
        energy = frame.total_energy_j();
        let transfers = system.last_window_transfers().unwrap();
        assert_eq!(transfers.len(), 2);
        assert!(system.last_window_transfer_error().is_none());
        let fs_ports::ModuleWindowAudit::Available(audit) = system.audit_last_window().unwrap() else {
            panic!("air-loaded prepared physical account unavailable")
        };
        assert!(audit.violations().is_empty());
    }
}

// The pole ceiling: the dense owner keeps its 256, the prepared owner admits a bank beyond it and still closes every
// step.
#[test]
fn the_prepared_owner_admits_an_air_bank_beyond_the_dense_ceiling() {
    let gate = CancelGate::new_clock_free();
    let poles: Vec<_> = (0..300).map(|i| {
        let x = f64::from(i) / 299.0;
        Pole { omega: TAU * (300.0 + 4700.0 * x), zeta: 0.05 + 0.4 * x,
            coupling: vec![12.0 * (1.0 - x) - 4.0, 9.0 * x - 3.0] }
    }).collect();
    let air = Model { ports: 2, poles };
    assert!(mechanics(48_000, true, false).with_radiation_load(&air, &[0, 1], 300, &gate).is_err());
    let mut system = mechanics(48_000, true, true).with_radiation_load(&air, &[0, 1], 300, &gate).unwrap();
    assert_eq!(system.radiation_state().len(), 600);
    let mut energy = system.total_energy_j().unwrap();
    for tick in 0..96 {
        let frame = system.step(&[0.001, if tick % 2 == 0 { -0.002 } else { 0.0 }]).unwrap().clone();
        let residual = frame.total_energy_j() + frame.total_dissipation_j() - frame.external_work_j - energy;
        assert!(residual.abs() <= 1e-11 * energy.max(frame.total_energy_j()), "tick {tick}: residual {residual:e}");
        energy = frame.total_energy_j();
    }
}

// A mode resting on its link's preload, released 1e-12 of its static deflection off rest, as a body mode rests on the
// strings' preload: its first increments are a few units in the last place of the deflection, so the two held port
// forces agree only to the rounding of those differences, not to 1e-9 of the force. Every step is admitted and closes
// its energy, and the mode stays at its rest.
#[test]
fn a_mode_resting_on_a_preload_steps_under_the_air_from_its_first_increments() {
    let gate = CancelGate::new_clock_free();
    let (omega, stiffness, extension) = (3000.0_f64, 2e6, 1e-3);
    let rest = stiffness * extension / (omega * omega + stiffness);
    let links = vec![ModalColumnConnection { column: vec![1.0], stiffness_n_m: stiffness, damping_n_s_m: 15.0,
        rest_extension_m: extension }];
    let models = vec![model(48_000, &[(omega, 0.005)], &[(rest * (1.0 + 1e-12), 0.0)])];
    let air = fixture(false).air;
    let mut system = CoupledModalSystem::new_with_columns_prepared(models, links, config(), budget(1), &gate).unwrap()
        .with_radiation_load(&air, &[0], 1, &gate).unwrap();
    let mut energy = system.total_energy_j().unwrap();
    for tick in 0..480 {
        let frame = system.step(&[0.0]).unwrap().clone();
        let residual = frame.total_energy_j() + frame.total_dissipation_j() - frame.external_work_j - energy;
        assert!(residual.abs() <= 1e-11 * energy.max(frame.total_energy_j()), "tick {tick}: residual {residual:e}");
        energy = frame.total_energy_j();
    }
    let q = flatten(&system.models)[0];
    assert!((q / rest - 1.0).abs() < 1e-11, "{q:e} left the rest {rest:e}");
}
