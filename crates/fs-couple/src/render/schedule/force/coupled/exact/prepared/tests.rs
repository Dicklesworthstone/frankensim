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
