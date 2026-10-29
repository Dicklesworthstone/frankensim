//! The existing prepared pencil against the original port-Hamiltonian owner,
//! a closed-form Maxwell oscillator, and the physical four-channel ledger.
use super::*;
use crate::modal_acoustic_time::{ModalAcousticMode, ModalAcousticTimeBudget};
use fs_math::c64::C64;
use fs_phs::RelaxationBranch;

fn budget(count: usize) -> PreparedModalBudget {
    PreparedModalBudget { max_modes: count, max_dense_bytes: 128 * 1024 * 1024,
        max_multiply_terms: 100_000_000_000 }
}

fn model(rate: u32, omega: &[f64], q: &[f64], v: &[f64]) -> ModalAcousticTimeModel {
    let mut result = ModalAcousticTimeModel::try_new_with_free_coordinates(rate,
        omega.iter().map(|&w| ModalAcousticMode { angular_frequency_rad_s: w,
            damping_ratio: if w == 0.0 { 0.0 } else { 0.01 },
            pressure_per_modal_velocity: C64::ZERO }).collect(),
        ModalAcousticTimeBudget::audible_reference()).unwrap();
    result.restore_states(&q.iter().zip(v).map(|(&q, &v)| ModalAcousticState {
        displacement_m_sqrt_kg: q, velocity_m_sqrt_kg_per_s: v,
    }).collect::<Vec<_>>()).unwrap();
    result
}

fn close(a: f64, b: f64) {
    assert!((a - b).abs() <= 2e-9 * a.abs().max(b.abs()) + 2e-14, "{a:e} != {b:e}");
}

fn branches() -> Vec<RelaxationBranch> {
    vec![
        RelaxationBranch { projection: vec![1.0, 0.0, -0.4, 0.0],
            stiffness: 3e4, relaxation_time_s: 0.002 },
        RelaxationBranch { projection: vec![0.6, 0.0, 1.0, 0.0],
            stiffness: 2e4, relaxation_time_s: 0.005 },
    ]
}

fn prepared(rate: u32) -> (PreparedNetwork, Vec<ModalAcousticTimeModel>, Vec<f64>) {
    let models = vec![model(rate, &[500.0, 900.0], &[1e-4, -2e-4], &[0.02, -0.01])];
    let links = vec![ConnectionLaw { stiffness_n_m: 4e4, damping_n_s_m: 20.0,
        rest_extension_m: 0.0, constant_loss_reference_omega_rad_s: None }];
    let owner = PreparedNetwork::new(&models, &[vec![1.0, -0.3]], &links, 2,
        1.0 / f64::from(rate), budget(2), &branches(), None,
        &CancelGate::new_clock_free()).unwrap();
    (owner, models, vec![0.007, -0.004])
}

fn publish(models: &mut [ModalAcousticTimeModel], memory: &mut [f64], endpoint: &[f64]) {
    let count = models[0].modes().len();
    models[0].restore_states(&(0..count).map(|i| ModalAcousticState {
        displacement_m_sqrt_kg: endpoint[2 * i], velocity_m_sqrt_kg_per_s: endpoint[2 * i + 1],
    }).collect::<Vec<_>>()).unwrap();
    memory.copy_from_slice(&endpoint[2 * count..2 * count + memory.len()]);
}

#[test]
fn prepared_maxwell_retains_full_noncommuting_storage_and_independent_loss_channels() {
    // Assemble the ORIGINAL canonical port-Hamiltonian equations independently
    // of the prepared generator, including signed inter-mode spring/dashpot work.
    let (owner, mut models, mut memory) = prepared(48_000);
    let (n, count) = (4, 2);
    let mut q = vec![0.0; n * n];
    let mut j = vec![0.0; n * n];
    let mut r = vec![0.0; n * n];
    let mut g = vec![0.0; n * count];
    for (i, mode) in models[0].modes().iter().enumerate() {
        q[(2 * i) * n + 2 * i] = mode.angular_frequency_rad_s.powi(2);
        q[(2 * i + 1) * n + 2 * i + 1] = 1.0;
        j[(2 * i) * n + 2 * i + 1] = 1.0;
        j[(2 * i + 1) * n + 2 * i] = -1.0;
        r[(2 * i + 1) * n + 2 * i + 1] = 2.0 * mode.damping_ratio * mode.angular_frequency_rad_s;
        g[(2 * i + 1) * count + i] = 1.0;
    }
    let b = [1.0, -0.3];
    for i in 0..count { for k in 0..count {
        q[(2 * i) * n + 2 * k] += 4e4 * b[i] * b[k];
        r[(2 * i + 1) * n + 2 * k + 1] += 20.0 * b[i] * b[k];
    }}
    let system = PortHamiltonian::new(n, count, j, r, g,
        Box::new(QuadraticStorage::new(q, n).unwrap())).unwrap()
        .with_relaxation_branches(branches()).unwrap();
    let split = split_quadratic(&system,
        LinearSplitBudget { max_modes: 6, max_omega_dt: core::f64::consts::PI }).unwrap();
    let mut truth = flatten(&models);
    truth.extend_from_slice(&memory);
    for tick in 0..12 {
        let force = [0.02, if tick % 2 == 0 { -0.01 } else { 0.03 }];
        let reference = step_render(&system, &split, &truth, &force, owner.dt).unwrap();
        let (actual, start, input) = owner.advance(&models, &memory, &[], &force, owner.dt).unwrap();
        for (&a, &b) in actual.x.iter().zip(&reference.x) { close(a, b); }
        close(actual.delta_h, reference.delta_h);
        close(actual.supplied, reference.supplied);
        close(actual.dissipated, reference.dissipated);
        let (component, connection, material, air) = owner.losses(&start, &input, owner.dt, actual.dissipated).unwrap();
        assert!(component > 0.0 && connection > 0.0 && material > 0.0);
        assert_eq!(air, 0.0);
        let mut resistance = vec![0.0; 36];
        resistance[4 * 6 + 4] = 1.0 / 0.002;
        resistance[5 * 6 + 5] = 1.0 / 0.005;
        close(material, integrated_resistance_loss(&system, &split, &truth, &force,
            owner.dt, &resistance).unwrap());
        let transfers = owner.ledger.work(&start, &input, &actual.x, owner.dt, connection,
            &owner.ledger.columns, &owner.ledger.springs, None).unwrap();
        let extension = |x: &[f64]| x[0] - 0.3 * x[2];
        close(transfers[0].energy_j,
            -0.5 * 4e4 * (extension(&actual.x).powi(2) - extension(&start).powi(2)) - connection);
        publish(&mut models, &mut memory, &actual.x);
        truth = reference.x;
    }
}

#[test]
fn prepared_memory_obeys_exact_time_partition_for_state_work_and_all_material_heat() {
    let (full, models, memory) = prepared(48_000);
    let (half, mut fine_models, mut fine_memory) = prepared(96_000);
    let force = [0.02, -0.01];
    let (a, start, input) = full.advance(&models, &memory, &[], &force, full.dt).unwrap();
    let (b, s1, u1) = half.advance(&fine_models, &fine_memory, &[], &force, half.dt).unwrap();
    let lb = half.losses(&s1, &u1, half.dt, b.dissipated).unwrap();
    publish(&mut fine_models, &mut fine_memory, &b.x);
    let (c, s2, u2) = half.advance(&fine_models, &fine_memory, &[], &force, half.dt).unwrap();
    let la = full.losses(&start, &input, full.dt, a.dissipated).unwrap();
    let lc = half.losses(&s2, &u2, half.dt, c.dissipated).unwrap();
    for (&x, &y) in a.x.iter().zip(&c.x) { close(x, y); }
    close(a.supplied, b.supplied + c.supplied);
    close(a.dissipated, b.dissipated + c.dissipated);
    close(la.0, lb.0 + lc.0); close(la.1, lb.1 + lc.1); close(la.2, lb.2 + lc.2);
}

#[test]
fn prepared_retains_more_material_states_than_the_small_exponential_capacity() {
    let (count, rate, q0) = (125_usize, 48_000, 1e-4_f64);
    let mut models = vec![model(rate, &[0.0], &[q0], &[0.0])];
    let arms = vec![RelaxationBranch { projection: vec![1.0, 0.0], stiffness: 100.0,
        relaxation_time_s: 0.002 }; count];
    let mut memory = vec![0.0; count];
    assert!(2 * (3 + count + 1) > fs_phs::linear_exact::MAX_LINEAR_EXPONENTIAL_BLOCK_DIM);
    let owner = PreparedNetwork::new(&models, &[], &[], 1, 1.0 / f64::from(rate),
        budget(1), &arms, None, &CancelGate::new_clock_free()).unwrap();
    let k = count as f64 * 100.0;
    let root = (250.0_f64.powi(2) - k).sqrt();
    let (slow, fast) = (-250.0 + root, -250.0 - root);
    let initial = 0.5 * k * q0 * q0;
    let mut heat = 0.0;
    for tick in 1..=4 {
        let (frame, start, input) = owner.advance(&models, &memory, &[], &[0.0], owner.dt).unwrap();
        let time = f64::from(tick) * owner.dt;
        let q = q0 * (-fast * (slow * time).exp() + slow * (fast * time).exp()) / (slow - fast);
        let v = q0 * slow * fast * ((fast * time).exp() - (slow * time).exp()) / (slow - fast);
        close(frame.x[0], q); close(frame.x[1], v);
        for &z in &frame.x[2..] { close(z, -v / (0.002 * count as f64 * 10.0)); }
        let (component, connection, material, air) = owner.losses(&start, &input, owner.dt, frame.dissipated).unwrap();
        assert_eq!((component, connection, air), (0.0, 0.0, 0.0));
        assert!(material > 0.0);
        heat += material;
        close(owner.hamiltonian(&frame.x) + heat, initial);
        publish(&mut models, &mut memory, &frame.x);
    }
    let wide = PreparedModalBudget { max_modes: 1, max_dense_bytes: usize::MAX, max_multiply_terms: u64::MAX };
    assert!(wide.admit_with_memory(1, fs_phs::linear_exact::MAX_PREPARED_EXPONENTIAL_DIM / 2,
        1, 0, 0, None).is_err());
    assert!(PreparedModalBudget { max_dense_bytes: 1024, ..wide }
        .admit_with_memory(1, count, 1, 0, 0, None).is_err());
    assert!(PreparedModalBudget { max_multiply_terms: 1024, ..wide }
        .admit_with_memory(1, count, 1, 0, 0, None).is_err());
}

#[test]
fn stiff_unloaded_material_history_releases_exact_heat_without_inverse_overflow() {
    let models = vec![model(48_000, &[0.0], &[0.0], &[0.0])];
    let dt = 1.0 / 48_000.0;
    let arms = [RelaxationBranch { projection: vec![0.0, 0.0], stiffness: 1.0,
        relaxation_time_s: dt / 1000.0 }];
    let owner = PreparedNetwork::new(&models, &[], &[], 1, dt, budget(1), &arms,
        None, &CancelGate::new_clock_free()).unwrap();
    let (step, start, input) = owner.advance(&models, &[0.02], &[], &[0.0], dt).unwrap();
    assert_eq!(&step.x[..2], &[0.0, 0.0]);
    assert!(step.x[2].abs() < 1e-100);
    let loss = owner.losses(&start, &input, dt, step.dissipated).unwrap();
    assert_eq!((loss.0, loss.1, loss.3), (0.0, 0.0, 0.0));
    close(loss.2, 0.5 * 0.02_f64.powi(2));
    close(step.delta_h + loss.2, 0.0);
    assert_eq!(step.supplied, 0.0);
}
