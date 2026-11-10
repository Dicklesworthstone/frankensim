use super::*;
use crate::modal_acoustic_time::{ModalAcousticMode, ModalAcousticTimeBudget};
use fs_math::c64::C64;
use fs_math::det;
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

/// At the mode ceiling the admission hands the exponential its Van Loan dimension, 6 x modes, and the budgeted
/// exponential accepts it; one more mode refuses.
#[test]
fn the_prepared_mode_ceiling_admits_its_van_loan_dimension_and_refuses_one_more_mode() {
    let wide = PreparedModalBudget {
        max_modes: MAX_PREPARED_MODAL_MODES,
        max_dense_bytes: usize::MAX,
        max_multiply_terms: u64::MAX,
    };
    let exponential = wide.admit(MAX_PREPARED_MODAL_MODES, 1, 1, 0, None).unwrap();
    assert_eq!(exponential.max_dimension, 6 * MAX_PREPARED_MODAL_MODES);
    let identity = fs_phs::linear_exact::exponential_with_budget(&[0.0], 1, exponential).unwrap();
    assert_eq!(identity, vec![1.0]);
    let over = PreparedModalBudget {
        max_modes: MAX_PREPARED_MODAL_MODES + 1,
        ..wide
    };
    assert!(over.validate().is_err());
    assert!(wide.admit(MAX_PREPARED_MODAL_MODES + 1, 1, 1, 0, None).is_err());
}

/// The wound E2 string at 96 modes per polarization with the four-interval string-air row, on the 128-mode body with its
/// cavity, is 128 + 1 + 96 x 2 x (1 + 4) = 1089 coordinates: past the former 1024 ceiling, inside this one. The admission
/// hands the exponential its Van Loan dimension, 6 x coordinates, and one coordinate past the ceiling refuses.
#[test]
fn the_prepared_admission_takes_the_wound_e2_string_with_the_string_air_row() {
    const E2_WITH_STRING_AIR: usize = 128 + 1 + 96 * 2 * (1 + 4);
    assert_eq!(E2_WITH_STRING_AIR, 1089);
    assert!(E2_WITH_STRING_AIR > 1024 && E2_WITH_STRING_AIR <= MAX_PREPARED_MODAL_MODES);
    let wide = PreparedModalBudget {
        max_modes: MAX_PREPARED_MODAL_MODES,
        max_dense_bytes: usize::MAX,
        max_multiply_terms: u64::MAX,
    };
    let exponential = wide.admit(E2_WITH_STRING_AIR, 1, 1, 0, None).unwrap();
    assert_eq!(exponential.max_dimension, 6 * E2_WITH_STRING_AIR);
    assert!(exponential.max_dimension <= fs_phs::linear_exact::MAX_PREPARED_EXPONENTIAL_DIM);
    assert!(wide.admit(MAX_PREPARED_MODAL_MODES + 1, 1, 1, 0, None).is_err());
}

/// A refused admission prints the bytes (or multiply terms) it asked for and the limit it was held to, not only that
/// the limit was exceeded: the card's author needs the figure to set the limit.
#[test]
fn a_refused_prepared_admission_prints_the_figure_it_asked_for_and_the_limit() {
    // count 3, no memory, 2 losses, 1 link, no constant-loss anchor, no air bank: d = 9, extent 81, 130 matrices,
    // (81 * 130) f64 entries and one link table entry.
    let asked = (81 * 130 * core::mem::size_of::<f64>() + core::mem::size_of::<Option<Vec<f64>>>()) as f64;
    let tight = PreparedModalBudget { max_modes: 3, max_dense_bytes: 84_000, max_multiply_terms: u64::MAX };
    let Err(refused) = tight.admit(3, 2, 1, 0, None) else { panic!("84000 bytes must not admit a 3-mode network") };
    match refused {
        ModalCouplingError::Budget { what, value, limit } => {
            assert_eq!(what, "prepared exact dense bytes (prepared-linear max_dense_bytes)");
            assert_eq!(value, asked);
            assert_eq!(limit, 84_000.0);
        }
        other => panic!("expected the budget refusal with its figures, got {other}"),
    }
    let text = tight.admit(3, 2, 1, 0, None).err().unwrap().to_string();
    assert_eq!(text, format!("coupled prepared exact dense bytes (prepared-linear max_dense_bytes): {asked:e} exceeds {:e}", 84_000.0_f64));
    assert!(text.contains("8.4264e4 exceeds 8.4e4"), "{text}");

    let roomy = PreparedModalBudget { max_modes: 3, max_dense_bytes: usize::MAX, max_multiply_terms: 1_000 };
    let Err(ModalCouplingError::Budget { what, value, limit }) = roomy.admit(3, 2, 1, 0, None) else {
        panic!("1000 multiply terms must not admit a 3-mode network")
    };
    assert_eq!(what, "prepared exact multiply terms (prepared-linear max_multiply_terms)");
    assert_eq!(limit, 1_000.0);
    assert!(value > limit && value.is_finite(), "{value}");
    // within both limits it still admits
    assert!(PreparedModalBudget { max_modes: 3, max_dense_bytes: 84_264, max_multiply_terms: u64::MAX }.admit(3, 2, 1, 0, None).is_ok());
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
    // `advance` and `losses` of a stage evaluate the same vector: one sweep (two dense forms) per step, and no per-link form.
    assert_eq!(sweeps(&system) - before, steps);
    let prepared = system.exact.as_ref().unwrap().prepared.as_ref().unwrap();
    assert_eq!(prepared.ledger_evaluations(), 0);
    // Neither construction nor stepping computed a link's own integral: only the component and summed-connection ones.
    assert!(prepared.ledger.gramians.get().is_none());
    let x = flatten(&system.models);
    for u in [[0.003, -0.002, 0.001], [0.0, 0.0, 0.0], [0.003, -0.002, 0.001]] {
        let z = prepared.state_input(&x, &u, prepared.dt).unwrap();
        let direct_component = quadratic(prepared.component_loss.as_ref().unwrap(), &z).unwrap();
        let direct_connection = quadratic(prepared.connection_total.as_ref().unwrap(), &z).unwrap();
        let (component, connection, material) = prepared.loss_sums(&z).unwrap();
        assert_eq!(material, 0.0);
        assert_eq!(component.to_bits(), direct_component.to_bits());
        assert_eq!(connection.to_bits(), direct_connection.to_bits());
        assert!(component > 0.0 && connection > 0.0);
        // The one integral of the summed resistance is the sum of the links' own integrals up to roundoff, inside the
        // ledger's 1e-11.
        let mut by_links = 0.0;
        for w in prepared.ledger.gramians(None).unwrap().iter().flatten() {
            by_links = finite(by_links + quadratic(w, &z).unwrap()).unwrap();
        }
        assert!((connection - by_links).abs() <= 1e-12 * connection.max(by_links), "{connection:e} vs {by_links:e}");
        // A second request for the same vector is served from the sweep, bit for bit.
        let count = prepared.sweeps();
        let (again_component, again_connection, _) = prepared.loss_sums(&z).unwrap();
        assert_eq!(prepared.sweeps(), count);
        assert_eq!((again_component.to_bits(), again_connection.to_bits()), (component.to_bits(), connection.to_bits()));
    }
}

#[test]
fn the_per_link_ledger_is_evaluated_only_when_asked_once_and_agrees_with_the_pre_summed_loss() {
    let mut system = fixture(48_000, true);
    let evaluations = |s: &CoupledModalSystem| s.exact.as_ref().unwrap().prepared.as_ref().unwrap().ledger_evaluations();
    let (mut asked, mut worst) = (0u64, 0.0_f64);
    for tick in 0..200u64 {
        let external = [0.003 * (tick as f64 * 0.37).sin(), if tick % 2 == 0 { -0.002 } else { 0.004 }, 0.001 * (tick % 7) as f64];
        // The sum of the links' own forms against the one pre-summed form, on the very vector this step uses.
        let relative = {
            let prepared = system.exact.as_ref().unwrap().prepared.as_ref().unwrap();
            let u: Vec<f64> = external.iter().zip(&prepared.rest_force).map(|(f, r)| f + r).collect();
            let z = prepared.state_input(&flatten(&system.models), &u, prepared.dt).unwrap();
            let total = quadratic(prepared.connection_total.as_ref().unwrap(), &z).unwrap();
            let by_links: f64 = prepared.ledger.gramians(None).unwrap().iter().flatten().map(|w| quadratic(w, &z).unwrap()).sum();
            (total - by_links).abs() / total.max(by_links).max(f64::MIN_POSITIVE)
        };
        worst = worst.max(relative);
        system.step(&external).unwrap();
        // Stepping evaluated no ledger; the first request evaluates it once and every later request shares it.
        assert_eq!(evaluations(&system), asked, "tick {tick}");
        if tick % 10 == 3 {
            let first: Vec<f64> = system.last_window_transfers().unwrap().iter().map(|t| t.energy_into().value()).collect();
            asked += 1;
            assert_eq!(evaluations(&system), asked);
            let again: Vec<f64> = system.last_window_transfers().unwrap().iter().map(|t| t.energy_into().value()).collect();
            assert_eq!(first, again);
            assert!(system.last_window_transfer_error().is_none(), "tick {tick}");
            assert_eq!(evaluations(&system), asked);
        }
    }
    eprintln!("summed-resistance vs per-link connection loss: worst relative difference {worst:e} over 200 steps (the ledger check allows 1e-11)");
    assert!(worst < 1e-12, "{worst:e}");
}

#[test]
fn free_prediction_keeps_its_energy_check_on_every_64th_sample_and_equals_the_full_advance_bit_for_bit() {
    let mut system = fixture(48_000, true);
    let sweeps = |s: &CoupledModalSystem| s.exact.as_ref().unwrap().prepared.as_ref().unwrap().sweeps();
    let gate = CancelGate::new_clock_free();
    for sample in 0..130u64 {
        let external = [0.003, if sample % 2 == 0 { -0.002 } else { 0.004 }, 0.001];
        assert_eq!(system.frame.sample, sample);
        // The free prediction the friction owner makes before the sample.
        let before = sweeps(&system);
        system.prepare_prediction(&external, Some(&gate)).unwrap();
        // The sentinel: a loss sweep (and with it the trajectory energy-residual check) only on samples 0, 64 and 128.
        assert_eq!(sweeps(&system) - before, u64::from(sample % 64 == 0), "sample {sample}");
        // Skipped or not, the prediction is exactly the full advance's end state, bit for bit.
        let (record, x, _) = system.exact.as_ref().unwrap().advance(&system.models, &[], &[], &external, system.dt).unwrap();
        for k in 0..system.mode_count() {
            assert_eq!(system.old_q[k].to_bits(), x[2 * k].to_bits());
            assert_eq!(system.free_delta[k].to_bits(), (record.x[2 * k] - x[2 * k]).to_bits());
        }
        system.step(&external).unwrap();
    }
}

/// Deterministic sign-mixed entries over forty binades, with exact zeros of both signs, for the bit-identity checks.
fn entries(count: usize, seed: &mut u64) -> Vec<f64> {
    (0..count).map(|i| {
        *seed = seed.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1_442_695_040_888_963_407);
        match i % 11 {
            3 => 0.0,
            7 => -0.0,
            _ => ((*seed >> 11) as f64 / (1u64 << 53) as f64 - 0.5) * 2f64.powi((*seed % 40) as i32 - 20),
        }
    }).collect()
}

#[test]
fn the_blocked_product_and_interleaved_apply_equal_the_textbook_loops_bit_for_bit() {
    let gate = CancelGate::new_clock_free();
    let mut seed = 0x9e37_79b9_7f4a_7c15_u64;
    // Sizes below, at and between multiples of the four-row block, and one the size of a small network.
    for n in [1, 2, 3, 4, 5, 7, 8, 9, 13, 64, 99] {
        let (a, b, z) = (entries(n * n, &mut seed), entries(n * n, &mut seed), entries(n, &mut seed));
        let mut textbook = vec![0.0; n * n];
        for i in 0..n {
            for k in 0..n {
                for j in 0..n {
                    textbook[i * n + j] += a[i * n + k] * b[k * n + j];
                }
            }
        }
        let blocked = product(&a, &b, n, &gate).unwrap();
        assert!(blocked.iter().zip(&textbook).all(|(x, y)| x.to_bits() == y.to_bits()), "product, n = {n}");
        let rows: Vec<f64> = a.chunks_exact(n).map(|row| row.iter().zip(&z).fold(0.0, |s, (x, y)| s + x * y)).collect();
        let applied = apply(&a, &z).unwrap();
        assert!(applied.iter().zip(&rows).all(|(x, y)| x.to_bits() == y.to_bits()), "apply, n = {n}");
    }
    // A row that overflows refuses, inside an eight-row block and in the remainder rows alike.
    for (n, row) in [(9, 1), (9, 8)] {
        let mut a = vec![1.0; n * n];
        a[row * n] = f64::MAX;
        assert!(apply(&a, &vec![2.0; n]).is_err(), "row {row}");
    }
}

/// The previous loss quadratic, a check of every term and partial sum, kept here as the reference the fused sweep must equal.
fn textbook_quadratic(matrix: &[f64], z: &[f64]) -> Result<f64, ModalCouplingError> {
    let (mut sum, mut magnitude) = (0.0, 0.0);
    for (i, row) in matrix.chunks_exact(z.len()).enumerate() {
        for (j, value) in row.iter().enumerate() {
            let term = finite(z[i] * value * z[j])?;
            sum = finite(sum + term)?;
            magnitude = finite(magnitude + term.abs())?;
        }
    }
    let roundoff = 4096.0 * f64::EPSILON * magnitude.max(f64::MIN_POSITIVE);
    if sum < -roundoff { return Err(invalid("prepared resistance integral is negative beyond roundoff")); }
    Ok(sum.max(0.0))
}

#[test]
fn the_fused_loss_pair_equals_two_textbook_quadratics_bit_for_bit_refusals_included() {
    let outcome = |r: Result<f64, ModalCouplingError>| r.map(f64::to_bits).map_err(|e| e.to_string());
    let mut seed = 0x2545_f491_4f6c_dd1d_u64;
    for n in [1, 2, 3, 5, 8, 13, 64, 99] {
        // Gramian-like PSD matrices B^T B, and the sign-mixed matrices themselves (some quadratics negative beyond roundoff).
        let (b1, b2, z) = (entries(n * n, &mut seed), entries(n * n, &mut seed), entries(n, &mut seed));
        let gram = |b: &[f64]| -> Vec<f64> {
            (0..n * n).map(|ij| (0..n).fold(0.0, |s, k| s + b[k * n + ij / n] * b[k * n + ij % n])).collect()
        };
        for (first, second) in [(gram(&b1), gram(&b2)), (b1.clone(), b2.clone()), (gram(&b1), b2.clone())] {
            let expected = (outcome(textbook_quadratic(&first, &z)), outcome(textbook_quadratic(&second, &z)));
            assert_eq!(outcome(quadratic(&first, &z)), expected.0, "n = {n}");
            let fused = quadratic_pair(&first, &second, &z);
            match &expected {
                (Ok(a), Ok(b)) => assert_eq!(fused.map(|(x, y)| (x.to_bits(), y.to_bits())).unwrap(), (*a, *b), "n = {n}"),
                (Err(e), _) | (Ok(_), Err(e)) => assert_eq!(fused.unwrap_err().to_string(), *e, "n = {n}"),
            }
        }
    }
    // A term that overflows refuses in either matrix, as the per-term check did.
    let mut huge = vec![1.0; 9];
    huge[4] = f64::MAX;
    for (first, second) in [(huge.clone(), vec![1.0; 9]), (vec![1.0; 9], huge.clone())] {
        assert!(textbook_quadratic(&huge, &[2.0; 3]).is_err());
        assert!(quadratic_pair(&first, &second, &[2.0; 3]).is_err());
    }
}

// ---- wohg3: the prepared owner against a closed form (the secular KAT) ----
//
// A 60 Hz harmonic string ladder (omega_i = 2 pi 60 i) on ONE bridge spring k = omega_1^2 at station 0.8137, with a
// dashpot. With proportional loss C = alpha K (zeta_i = alpha omega_i / 2, c = alpha k) the step has a closed form: the
// eigenpairs of K = W^2 + k phi phi^T from the secular equation, then one damped 2x2 exponential and one held-force
// integral per coupled mode. The prepared transition, the budgeted exponential and a 1000-step trajectory are held to the
// stated bound of fs-phs CONTRACT.md ("Stated accuracy") at 40, 64 and 200 modes (dimensions 120, 192 and 600), and the
// bound must reject two planted wrong answers. The non-proportional dashpot has no closed form; there the ordinary and
// prepared owners are compared with each other.

const KAT_DT: f64 = 1.0 / 48_000.0;

/// The ladder's modes, the bridge column, spring and dashpot, and alpha of the proportional loss.
struct Ladder {
    omega: Vec<f64>,
    zeta: Vec<f64>,
    phi: Vec<f64>,
    k: f64,
    c: f64,
    alpha: f64,
}

/// `k_factor` scales the bridge spring (0 drops the bridge); proportional loss keeps c = alpha k, the non-proportional
/// dashpot is zeta = 1e-3 on every mode and c = 300 N s/m.
fn ladder(modes: usize, k_factor: f64, proportional: bool) -> Ladder {
    let w1 = 2.0 * core::f64::consts::PI * 60.0;
    let alpha = 2e-3 / w1;
    let omega: Vec<f64> = (1..=modes).map(|i| w1 * i as f64).collect();
    let phi = (1..=modes).map(|i| core::f64::consts::SQRT_2 * det::sin(i as f64 * core::f64::consts::PI * 0.8137)).collect();
    let k = w1 * w1 * k_factor;
    let (zeta, c) = if proportional {
        (omega.iter().map(|w| 0.5 * alpha * w).collect(), alpha * k)
    } else {
        (vec![1e-3; modes], 300.0)
    };
    Ladder { omega, zeta, phi, k, c, alpha }
}

fn ladder_start(l: &Ladder) -> Vec<ModalAcousticState> {
    l.omega.iter().enumerate().map(|(i, w)| ModalAcousticState {
        displacement_m_sqrt_kg: 0.02 / w * det::sin(1.3 * i as f64 + 0.2),
        velocity_m_sqrt_kg_per_s: 0.02 * det::cos(0.7 * i as f64 + 0.1),
    }).collect()
}

fn ladder_force(modes: usize) -> Vec<f64> {
    (0..modes).map(|i| 0.5 * det::sin(2.1 * i as f64 + 0.3)).collect()
}

fn ladder_system(l: &Ladder, prepared: bool) -> Result<CoupledModalSystem, ModalCouplingError> {
    let modes = l.omega.len();
    let mut model = ModalAcousticTimeModel::try_new(48_000, l.omega.iter().zip(&l.zeta).map(|(&w, &z)| ModalAcousticMode {
        angular_frequency_rad_s: w, damping_ratio: z, pressure_per_modal_velocity: C64::ZERO }).collect(),
        ModalAcousticTimeBudget::audible_reference()).unwrap();
    model.restore_states(&ladder_start(l)).unwrap();
    let links = vec![ModalColumnConnection { column: l.phi.clone(), stiffness_n_m: l.k, damping_n_s_m: l.c, rest_extension_m: 0.0 }];
    let config = ModalCouplingConfig { max_setup_terms: 1_000_000_000, ..config(modes, 1) };
    let gate = CancelGate::new_clock_free();
    if prepared {
        let budget = PreparedModalBudget { max_modes: modes, max_dense_bytes: 1 << 30, max_multiply_terms: 2_000_000_000_000 };
        CoupledModalSystem::new_with_columns_prepared(vec![model], links, config, budget, &gate)
    } else {
        CoupledModalSystem::new_with_columns(vec![model], links, config, &gate)
    }
}

/// Eigenpairs of K = diag(d) + k z z^T (d increasing, every z_i nonzero, k > 0). Each root mu_j is bisected to adjacent
/// floats in a variable shifted to its nearer pole and kept as (pole, offset), so mu_j - d_i keeps its relative
/// accuracy; the vectors follow Gu and Eisenstat (1994): the z for which these roots are exact, then (D - mu_j)^-1 z.
struct Secular {
    d: Vec<f64>,
    roots: Vec<(usize, f64)>,
    /// Row-major; column j belongs to mu_j.
    v: Vec<f64>,
}

impl Secular {
    fn new(l: &Ladder) -> Self {
        let m = l.omega.len();
        let d: Vec<f64> = l.omega.iter().map(|w| w * w).collect();
        let mut s = Self { roots: (0..m).map(|i| (i, 0.0)).collect(), v: vec![0.0; m * m], d };
        if l.k == 0.0 {
            (0..m).for_each(|i| s.v[i * m + i] = 1.0);
            return s;
        }
        let (d, z, k) = (&s.d, &l.phi, l.k);
        // f rises from -inf to +inf on each (d_j, d_j+1); the last root lies below d_m + k |z|^2.
        let f = |pole: usize, x: f64| 1.0 + k * (0..m).fold(0.0, |sum, i| sum + z[i] * z[i] / ((d[i] - d[pole]) - x));
        let reach = k * z.iter().fold(0.0, |sum, v| sum + v * v);
        s.roots = (0..m).map(|j| {
            let upper = if j + 1 < m { d[j + 1] } else { d[j] + reach };
            let pole = if j + 1 < m && f(j, 0.5 * (upper - d[j])) <= 0.0 { j + 1 } else { j };
            let (mut lo, mut hi) = (d[j] - d[pole], upper - d[pole]);
            loop {
                let mid = 0.5 * (lo + hi);
                if mid <= lo || mid >= hi { break (pole, mid); }
                if f(pole, mid) > 0.0 { hi = mid } else { lo = mid }
            }
        }).collect();
        let zhat: Vec<f64> = (0..m).map(|i| {
            let mut p = s.gap(i, m - 1) / k;
            for j in 0..i { p *= s.gap(i, j) / (s.d[j] - s.d[i]); }
            for j in i..m - 1 { p *= s.gap(i, j) / (s.d[j + 1] - s.d[i]); }
            assert!(p > 0.0, "secular interlacing failed at {i}");
            det::sqrt(p).copysign(z[i])
        }).collect();
        for j in 0..m {
            let column: Vec<f64> = (0..m).map(|i| -zhat[i] / s.gap(i, j)).collect();
            let norm = det::sqrt(column.iter().fold(0.0, |sum, x| sum + x * x));
            (0..m).for_each(|i| s.v[i * m + j] = column[i] / norm);
        }
        s
    }

    /// mu_j - d_i.
    fn gap(&self, i: usize, j: usize) -> f64 { (self.d[self.roots[j].0] - self.d[i]) + self.roots[j].1 }

    fn mu(&self, j: usize) -> f64 { self.d[self.roots[j].0] + self.roots[j].1 }

    /// (orthogonality defect ||V^T V - I||_inf, largest relative residual ||K v_j - mu_j v_j||_2 / mu_j).
    fn defects(&self, l: &Ladder) -> (f64, f64) {
        let m = self.d.len();
        let (mut eta, mut rho) = (0.0_f64, 0.0_f64);
        for j in 0..m {
            let col = |i: usize| self.v[i * m + j];
            eta = eta.max((0..m).fold(0.0, |sum, c| {
                let dot = (0..m).fold(0.0, |acc, i| acc + col(i) * self.v[i * m + c]);
                sum + (dot - if c == j { 1.0 } else { 0.0 }).abs()
            }));
            let zv = (0..m).fold(0.0, |acc, i| acc + l.phi[i] * col(i));
            let r = (0..m).fold(0.0, |acc, i| {
                let x = self.d[i] * col(i) + l.k * l.phi[i] * zv - self.mu(j) * col(i);
                acc + x * x
            });
            rho = rho.max(det::sqrt(r) / self.mu(j));
        }
        (eta, rho)
    }
}

/// One coupled mode's exact step over t under C = alpha K: the (q, v) transition [E11, E12, E21, E22] and the held unit
/// force's displacement and velocity [G1, G2].
fn modal_step(mu: f64, alpha: f64, t: f64) -> [f64; 6] {
    let sigma = 0.5 * alpha * mu;
    assert!(mu > sigma * sigma, "the KAT's coupled modes are underdamped");
    let wd = det::sqrt(mu - sigma * sigma);
    let e = det::exp(-sigma * t);
    let (s, c) = det::sin_cos(wd * t);
    let half = det::sin(0.5 * wd * t);
    [e * (c + sigma * s / wd), e * s / wd, -e * mu * s / wd, e * (c - sigma * s / wd),
        (-det::expm1(-sigma * t) * c + 2.0 * half * half - sigma / wd * e * s) / mu, e * s / wd]
}

/// The exact held-force transition over t in the prepared owner's coordinates z = [w_i q_i, v_i, ..., u_i] (3m x 3m).
fn closed_form(l: &Ladder, s: &Secular, w: &[f64], t: f64) -> Vec<f64> {
    let (m, n) = (w.len(), 3 * w.len());
    let blocks: Vec<[f64; 6]> = (0..m).map(|j| modal_step(s.mu(j), l.alpha, t)).collect();
    let mut f = vec![0.0; n * n];
    for i in 0..m {
        for c in 0..m {
            let mut sum = [0.0; 6];
            for (j, block) in blocks.iter().enumerate() {
                let vv = s.v[i * m + j] * s.v[c * m + j];
                sum.iter_mut().zip(block).for_each(|(acc, b)| *acc += vv * b);
            }
            f[2 * i * n + 2 * c] = w[i] * sum[0] / w[c];
            f[2 * i * n + 2 * c + 1] = w[i] * sum[1];
            f[(2 * i + 1) * n + 2 * c] = sum[2] / w[c];
            f[(2 * i + 1) * n + 2 * c + 1] = sum[3];
            f[2 * i * n + 2 * m + c] = w[i] * sum[4];
            f[(2 * i + 1) * n + 2 * m + c] = sum[5];
        }
        f[(2 * m + i) * n + 2 * m + i] = 1.0;
    }
    f
}

/// The exact z state (2m entries) at time t from the fixture's start under its constant held force.
fn closed_path<'a>(l: &'a Ladder, s: &'a Secular, w: &'a [f64]) -> impl Fn(f64) -> Vec<f64> + 'a {
    let m = w.len();
    let (start, force) = (ladder_start(l), ladder_force(m));
    let project = |x: &dyn Fn(usize) -> f64| -> Vec<f64> {
        (0..m).map(|j| (0..m).fold(0.0, |acc, i| acc + s.v[i * m + j] * x(i))).collect()
    };
    let (q0, v0, f0) = (project(&|i| start[i].displacement_m_sqrt_kg), project(&|i| start[i].velocity_m_sqrt_kg_per_s),
        project(&|i| force[i]));
    move |t| {
        let mut z = vec![0.0; 2 * m];
        for j in 0..m {
            let b = modal_step(s.mu(j), l.alpha, t);
            let (q, v) = (b[0] * q0[j] + b[1] * v0[j] + b[4] * f0[j], b[2] * q0[j] + b[3] * v0[j] + b[5] * f0[j]);
            for i in 0..m {
                z[2 * i] += w[i] * s.v[i * m + j] * q;
                z[2 * i + 1] += s.v[i * m + j] * v;
            }
        }
        z
    }
}

fn gamma(count: usize) -> f64 {
    let x = count as f64 * 0.5 * f64::EPSILON;
    x / (1.0 - x)
}

/// ||a||_inf, each row summed left to right as fs-phs sums the operator norm that chooses its squarings.
fn norm_inf(a: &[f64], n: usize) -> f64 {
    a.chunks_exact(n).fold(0.0, |m, row| m.max(row.iter().fold(0.0, |s, v| s + v.abs())))
}

fn distance(a: &[f64], b: &[f64], n: usize) -> f64 {
    a.chunks_exact(n).zip(b.chunks_exact(n)).fold(0.0, |m, (x, y)| m.max(x.iter().zip(y).fold(0.0, |s, (p, q)| s + (p - q).abs())))
}

fn largest_gap(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).fold(0.0, |m, (x, y)| m.max((x - y).abs()))
}

/// The squarings exponential_with_budget selects for an operator of computed norm `norm`, and the Taylor-stage norm
/// theta, inflated by gamma_n for the norm's own rounding.
fn taylor_stage(norm: f64, n: usize) -> (usize, f64) {
    let (mut scale, mut squarings) = (1.0, 0);
    while norm * scale > 0.5 {
        scale *= 0.5;
        squarings += 1;
    }
    (squarings, norm * scale * (1.0 + gamma(n)))
}

/// fs-phs CONTRACT.md "Stated accuracy": the bound on ||computed - exp(X)||_inf after the Taylor stage at norm theta and
/// `squarings` squarings. `level(i, e_i)` returns G_i >= ||exp(2^i X')||_inf; `input` is the relative rounding of the
/// operator's entries against the exact one.
fn stated_bound(theta: f64, n: usize, squarings: usize, input: f64, mut level: impl FnMut(usize, f64) -> f64) -> f64 {
    let (g, sums) = (gamma(n + 1), gamma(32));
    let terms = det::exp(theta) * det::expm1(theta * det::expm1(g));
    let mut e = terms + sums * (1.0 + sums) * det::exp(theta * det::exp(g))
        + (n as f64 * 0.25 * f64::EPSILON + terms) * theta / (2.0 - theta)
        + input * theta * det::exp(theta * (1.0 + input));
    for i in 0..squarings {
        let big = level(i, e);
        e = (2.0 * big + e) * e + gamma(n) * (big + e) * (big + e);
    }
    e
}

/// What the prepared transition's bound needs. Its exponential received the generator times the network's short
/// interval h and chose `inner` squarings of its own; the network then squared `levels.len() - 1` times. Level i of
/// the whole chain is exp(2^i X') with X' that operator over `tau`.
struct PreparedStage {
    weights: Vec<f64>,
    transition: Vec<f64>,
    generator: Vec<f64>,
    levels: Vec<Vec<f64>>,
    inner: usize,
    theta: f64,
    tau: f64,
}

impl PreparedStage {
    fn of(system: &CoupledModalSystem) -> Self {
        let prepared = system.exact.as_ref().unwrap().prepared.as_ref().unwrap();
        let integral = prepared.ledger.integral.as_ref().unwrap();
        let n = prepared.weights.len() * 3;
        let scaled: Vec<f64> = integral.generator.iter().map(|v| v * integral.h).collect();
        let (inner, theta) = taylor_stage(norm_inf(&scaled, n), n);
        Self { weights: prepared.weights.clone(), transition: prepared.transition.clone(),
            generator: integral.generator.clone(), levels: integral.levels.clone(), inner, theta,
            tau: integral.h * 0.5_f64.powi(inner as i32) }
    }

    fn squarings(&self) -> usize { self.inner + self.levels.len() - 1 }

    fn bound(&self, level: impl FnMut(usize, f64) -> f64) -> f64 {
        stated_bound(self.theta, 3 * self.weights.len(), self.squarings(), gamma(8), level)
    }

    /// Level norms without a closed form: exp(2^i theta) inside the exponential, the computed level plus its own bound
    /// after it.
    fn computed_level(&self, i: usize, e: f64) -> f64 {
        if i < self.inner {
            det::exp(self.theta * 2.0_f64.powi(i as i32) * (1.0 + gamma(8)))
        } else {
            norm_inf(&self.levels[i - self.inner], 3 * self.weights.len()) + e
        }
    }
}

/// The owner's z trajectory over `steps` held-force steps (z_0 first), and the worst energy balance of a step, from
/// the frames: |dH + loss - work| over the step's largest term.
fn walk(system: &mut CoupledModalSystem, w: &[f64], steps: usize) -> (Vec<Vec<f64>>, f64) {
    let force = ladder_force(w.len());
    let z = |s: &CoupledModalSystem| -> Vec<f64> {
        s.components()[0].states().iter().zip(w)
            .flat_map(|(x, w)| [w * x.displacement_m_sqrt_kg, x.velocity_m_sqrt_kg_per_s]).collect()
    };
    let (mut path, mut worst) = (vec![z(system)], 0.0_f64);
    for _ in 0..steps {
        let before = system.total_energy_j().unwrap();
        let frame = system.step(&force).unwrap().clone();
        let after = system.total_energy_j().unwrap();
        let loss = frame.total_dissipation_j();
        let scale = before.abs().max(after.abs()).max(frame.external_work_j.abs()).max(loss);
        worst = worst.max(((after - before) + loss - frame.external_work_j).abs() / scale);
        path.push(z(system));
    }
    (path, worst)
}

/// e_N = P sum_{k<N} (B + gamma_{n+4} ||T||_inf) ||[z_k, u]||_inf for N = 0..=steps. P = sqrt(2m kappa), kappa the
/// condition of the balanced storage blockdiag(W^-1 K W^-1, I): the exact passive step never increases that energy norm,
/// so P bounds every power of it in the inf-norm. The local term covers the transition's error, the apply and the
/// weight scaling's roundoff.
fn trajectory_bound(l: &Ladder, w: &[f64], path: &[Vec<f64>], per_step: f64) -> Vec<f64> {
    let hi = 1.0 + l.k * l.phi.iter().zip(w).fold(0.0, |s, (p, w)| s + (p / w) * (p / w));
    let lo = l.omega.iter().zip(w).fold(1.0_f64, |m, (o, w)| m.min(o * o / (w * w)));
    let p = det::sqrt(2.0 * w.len() as f64 * hi / lo);
    let held = ladder_force(w.len()).iter().fold(0.0_f64, |m, f| m.max(f.abs()));
    let mut sum = 0.0;
    let mut bound = vec![0.0];
    for z in &path[..path.len() - 1] {
        sum += per_step * z.iter().fold(held, |m, x| m.max(x.abs()));
        bound.push(p * sum);
    }
    bound
}

/// The owner's one-step map in z coordinates, column by column from unit single steps (a unit held force for the input
/// columns); the input rows are the identity.
fn stepped_transition(system: &mut CoupledModalSystem, w: &[f64]) -> Vec<f64> {
    let (m, n) = (w.len(), 3 * w.len());
    let mut t = vec![0.0; n * n];
    for col in 0..n {
        let mut states = vec![ModalAcousticState::default(); m];
        let mut external = vec![0.0; m];
        match col {
            c if c >= 2 * m => external[c - 2 * m] = 1.0,
            c if c % 2 == 0 => states[c / 2].displacement_m_sqrt_kg = 1.0 / w[c / 2],
            c => states[c / 2].velocity_m_sqrt_kg_per_s = 1.0,
        }
        system.models[0].restore_states(&states).unwrap();
        system.step(&external).unwrap();
        for (i, (x, w)) in system.components()[0].states().iter().zip(w).enumerate() {
            t[2 * i * n + col] = w * x.displacement_m_sqrt_kg;
            t[(2 * i + 1) * n + col] = x.velocity_m_sqrt_kg_per_s;
        }
    }
    (0..m).for_each(|i| t[(2 * m + i) * n + 2 * m + i] = 1.0);
    t
}

#[test]
fn prepared_owner_meets_the_stated_bound_against_the_secular_closed_form_at_40_64_and_200_modes() {
    let steps = 1000;
    for modes in [40usize, 64, 200] {
        let (l, n) = (ladder(modes, 1.0, true), 3 * modes);
        let eig = Secular::new(&l);
        let (eta, rho) = eig.defects(&l);
        // Above 42 modes the ordinary owner refuses on its 256 block; nothing is dropped.
        let ordinary = ladder_system(&l, false);
        if modes > 42 {
            let Err(ModalCouplingError::Propagation(fs_phs::PhsError::LinearSplitBudget { what })) = &ordinary else {
                panic!("{modes} modes: the ordinary owner did not refuse on its block capacity")
            };
            assert!(what.contains("dense exponential block capacity"), "{what}");
        }
        let started = std::time::Instant::now();
        let mut system = ladder_system(&l, true).unwrap();
        let built_s = started.elapsed().as_secs_f64();
        let st = PreparedStage::of(&system);
        let w = st.weights.clone();
        let bound = st.bound(|i, _| {
            norm_inf(&closed_form(&l, &eig, &w, st.tau * 2.0_f64.powi(i as i32)), n) * (1.0 + 1e-12)
        });
        let (transition, stage) = (&st.transition, format!("{}+{} squarings, theta {:.4}", st.inner, st.levels.len() - 1, st.theta));
        let exact = closed_form(&l, &eig, &w, KAT_DT);
        let error = distance(transition, &exact, n);
        // The same operator through fs-phs alone: X = generator * dt, the exponential choosing its own squarings; the
        // ordinary entry runs it bit for bit up to dimension 256 and refuses above.
        let operator: Vec<f64> = st.generator.iter().map(|v| v * KAT_DT).collect();
        let raw = exponential_with_budget(&operator, n, DenseExponentialBudget {
            max_dimension: n, max_dense_bytes: usize::MAX, max_multiply_terms: u64::MAX }).unwrap();
        let (squarings, theta) = taylor_stage(norm_inf(&operator, n), n);
        let raw_bound = stated_bound(theta, n, squarings, gamma(8), |i, _| {
            norm_inf(&closed_form(&l, &eig, &w, KAT_DT * 0.5_f64.powi((squarings - i) as i32)), n) * (1.0 + 1e-12)
        });
        let raw_error = distance(&raw, &exact, n);
        match fs_phs::linear_exact::exponential(&operator, n) {
            Ok(plain) => assert!(n <= 256 && plain.iter().zip(&raw).all(|(a, b)| a.to_bits() == b.to_bits())),
            Err(fs_phs::PhsError::LinearSplitBudget { what }) => {
                assert!(n > 256 && what == "invalid or oversized exponential operator", "{what}");
            }
            Err(other) => panic!("{other}"),
        }
        // A 1000-step trajectory under a constant held force, against the exact modal solution.
        let started = std::time::Instant::now();
        let (path, energy) = walk(&mut system, &w, steps);
        let step_ms = 1e3 * started.elapsed().as_secs_f64() / steps as f64;
        let bounds = trajectory_bound(&l, &w, &path, bound + gamma(n + 4) * norm_inf(transition, n));
        let exact_path = closed_path(&l, &eig, &w);
        let errors: Vec<f64> = path.iter().enumerate().map(|(k, z)| largest_gap(z, &exact_path(k as f64 / 48_000.0))).collect();
        let worst = errors.iter().zip(&bounds).skip(1).fold(0.0_f64, |m, (e, b)| m.max(e / b));
        println!("{{\"test\":\"wohg3_secular_kat\",\"modes\":{modes},\"n\":{n},\"transition_error\":{error:.3e},\"bound\":{bound:.3e},\"ratio\":{:.1},\"stage\":\"{stage}\",\"raw_error\":{raw_error:.3e},\"raw_bound\":{raw_bound:.3e},\"raw_squarings\":{squarings},\"oracle_eta\":{eta:.1e},\"oracle_rho\":{rho:.1e},\"path_error\":{:.3e},\"path_bound\":{:.3e},\"worst_path_ratio\":{worst:.2e},\"energy_residual\":{energy:.1e},\"build_s\":{built_s:.2},\"step_ms\":{step_ms:.3}}}",
            bound / error, errors[steps], bounds[steps]);
        assert!(error <= bound, "{modes} modes: transition error {error:e} over the stated bound {bound:e}");
        assert!(raw_error <= raw_bound, "{modes} modes: exponential error {raw_error:e} over the stated bound {raw_bound:e}");
        assert!(worst <= 1.0, "{modes} modes: trajectory error over its bound (worst ratio {worst:e})");
        assert!(energy <= 1e-11, "{modes} modes: energy residual {energy:e}");
        // The planted wrong answers: each must fail the bound that the true answer met.
        for (name, wrong) in [("bridge dropped", ladder(modes, 0.0, true)), ("bridge k x (1 + 1e-6)", ladder(modes, 1.0 + 1e-6, true))] {
            let planted_eig = Secular::new(&wrong);
            let planted = closed_form(&wrong, &planted_eig, &w, KAT_DT);
            let (against, raw_against) = (distance(transition, &planted, n), distance(&raw, &planted, n));
            let path_against = largest_gap(&path[steps], &closed_path(&wrong, &planted_eig, &w)(steps as f64 / 48_000.0));
            println!("{{\"test\":\"wohg3_secular_kat\",\"modes\":{modes},\"planted\":\"{name}\",\"transition_error\":{against:.3e},\"over_bound\":{:.1e},\"raw_over_bound\":{:.1e},\"path_over_bound\":{:.1e}}}",
                against / bound, raw_against / raw_bound, path_against / bounds[steps]);
            assert!(against > bound && raw_against > raw_bound, "{modes} modes: planted '{name}' passed the bound");
            if name == "bridge dropped" {
                assert!(path_against > bounds[steps], "{modes} modes: the dropped bridge passed the trajectory bound");
            }
        }
        // At 40 modes the ordinary (eigen-split) owner is held to the same bounds: its one-step map and trajectory.
        if let Ok(mut ordinary) = ordinary {
            let stepped = stepped_transition(&mut ladder_system(&l, false).unwrap(), &w);
            let ordinary_error = distance(&stepped, &exact, n);
            let (ordinary_path, ordinary_energy) = walk(&mut ordinary, &w, steps);
            let ordinary_worst = ordinary_path.iter().enumerate().skip(1)
                .fold(0.0_f64, |m, (k, z)| m.max(largest_gap(z, &exact_path(k as f64 / 48_000.0)) / bounds[k]));
            let local = bound + gamma(n + 4) * norm_inf(transition, n);
            println!("{{\"test\":\"wohg3_secular_kat\",\"modes\":{modes},\"owner\":\"ordinary\",\"stepped_transition_error\":{ordinary_error:.3e},\"local_bound\":{local:.3e},\"worst_path_ratio\":{ordinary_worst:.2e},\"energy_residual\":{ordinary_energy:.1e}}}");
            assert!(ordinary_error <= local && ordinary_worst <= 1.0 && ordinary_energy <= 1e-11);
        }
    }
}

#[test]
fn ordinary_and_prepared_owners_agree_within_the_stated_bound_on_a_non_proportional_bridge_dashpot() {
    // zeta = 1e-3 on every mode plus a 300 N s/m dashpot at the bridge: C is not a multiple of K, so no closed form. The
    // ordinary (Jacobi eigen-split) and prepared (Taylor) owners are compared with each other. The cross-check bound
    // is twice the prepared bound: it ASSUMES the ordinary owner obeys the prepared owner's stated per-step bound, which
    // the proportional test checks at this size but nothing proves.
    let (steps, modes) = (1000, 40);
    let (l, n) = (ladder(modes, 1.0, false), 3 * modes);
    let (mut prepared, mut ordinary) = (ladder_system(&l, true).unwrap(), ladder_system(&l, false).unwrap());
    let st = PreparedStage::of(&prepared);
    let bound = st.bound(|i, e| st.computed_level(i, e));
    let local = bound + gamma(n + 4) * norm_inf(&st.transition, n);
    let w = &st.weights;
    let cross = distance(&stepped_transition(&mut ladder_system(&l, false).unwrap(), w),
        &stepped_transition(&mut ladder_system(&l, true).unwrap(), w), n);
    let (prepared_path, prepared_energy) = walk(&mut prepared, w, steps);
    let (ordinary_path, ordinary_energy) = walk(&mut ordinary, w, steps);
    let bounds = trajectory_bound(&l, w, &prepared_path, local);
    let worst = prepared_path.iter().zip(&ordinary_path).zip(&bounds).skip(1)
        .fold(0.0_f64, |m, ((a, b), e)| m.max(largest_gap(a, b) / (2.0 * e)));
    let last = largest_gap(&prepared_path[steps], &ordinary_path[steps]);
    println!("{{\"test\":\"wohg3_cross_check\",\"modes\":{modes},\"n\":{n},\"transition_difference\":{cross:.3e},\"bound\":{:.3e},\"ratio\":{:.1},\"path_difference\":{last:.3e},\"path_bound\":{:.3e},\"worst_path_ratio\":{worst:.2e},\"energy_residuals\":[{prepared_energy:.1e},{ordinary_energy:.1e}]}}",
        2.0 * local, 2.0 * local / cross, 2.0 * bounds[steps]);
    assert!(cross <= 2.0 * local, "one-step maps differ by {cross:e}, over twice the stated bound {local:e}");
    assert!(worst <= 1.0, "trajectories differ beyond twice their bound (worst ratio {worst:e})");
    assert!(prepared_energy <= 1e-11 && ordinary_energy <= 1e-11, "{prepared_energy:e} {ordinary_energy:e}");
}
