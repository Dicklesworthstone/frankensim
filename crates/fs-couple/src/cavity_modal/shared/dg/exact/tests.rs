use super::super::tests::{cavity, load};
use super::*;

fn budget() -> ExactDgBudget {
    ExactDgBudget {
        max_dimension: 256,
        max_dense_bytes: 64 * 1024 * 1024,
        max_multiply_terms: 20_000_000_000,
    }
}
fn create(rate: u32, memory: bool) -> SharedCavityDg {
    cavity(rate, memory, 1e6)
        .with_exact_downstream_duct(
            load(rate, memory, 0.001),
            budget(),
            &CancelGate::new_clock_free(),
        )
        .unwrap()
}
fn close(a: f64, b: f64, absolute: f64) {
    assert!(
        (a - b).abs() <= absolute + 2e-8 * a.abs().max(b.abs()),
        "{a:e} != {b:e}"
    );
}
fn same_state(a: &SharedCavityDg, b: &SharedCavityDg) {
    for (ma, mb) in a
        .cavity()
        .network()
        .components()
        .iter()
        .zip(b.cavity().network().components())
    {
        for (sa, sb) in ma.states().iter().zip(mb.states()) {
            close(sa.displacement_m_sqrt_kg, sb.displacement_m_sqrt_kg, 2e-16);
            close(
                sa.velocity_m_sqrt_kg_per_s,
                sb.velocity_m_sqrt_kg_per_s,
                2e-12,
            );
        }
    }
    for (x, y) in a
        .cavity()
        .network()
        .material_memory()
        .iter()
        .zip(b.cavity().network().material_memory())
    {
        close(*x, *y, 2e-12);
    }
    for (x, y) in a
        .load()
        .state()
        .pressure()
        .iter()
        .zip(b.load().state().pressure())
    {
        close(*x, *y, 2e-8);
    }
    for (x, y) in a
        .load()
        .state()
        .volume_flow()
        .iter()
        .zip(b.load().state().volume_flow())
    {
        close(*x, *y, 2e-16);
    }
    close(
        a.stored_energy_j().unwrap(),
        b.stored_energy_j().unwrap(),
        2e-15,
    );
}

#[test]
fn complete_joint_trajectory_is_invariant_under_three_clock_partitions() {
    // Cooperative cancellation reaches a later cold form, with no cache or
    // physical state published and no modification to the original kernel.
    let source = fs_phs::modal_bank_ports(&[100.0], &[0.0], &[&[1.0]]).unwrap();
    let split = fs_phs::linear_exact::split_quadratic(
        &source,
        fs_phs::linear_exact::LinearSplitBudget {
            max_modes: 2,
            max_omega_dt: 3.0,
        },
    )
    .unwrap();
    let mut polls = 0;
    let cancelled = fs_phs::linear_exact::prepared::PreparedQuadraticFlow::new_controlled(
        &source,
        &split,
        1.0 / 48000.0,
        &[vec![0.0, 0.0, 0.0, 1.0], vec![0.0, 0.0, 0.0, 1.0]],
        &[],
        budget(),
        || {
            polls += 1;
            polls == 15
        },
    );
    assert!(matches!(
        cancelled,
        Err(fs_phs::PreparedStepError::Cancelled)
    ));
    let gate = CancelGate::new_clock_free();
    let mut whole = create(24000, false);
    let mut halves = create(48000, false);
    let mut quarters = create(96000, false);
    let mut covariance = 0.0_f64;
    for _ in 0..32 {
        let frame = whole.step(&[0.001], 0.1, &gate).unwrap().clone();
        let mut half = [0.0; 3];
        let mut quarter = [0.0; 3];
        for (system, count, totals) in [
            (&mut halves, 2, &mut half),
            (&mut quarters, 4, &mut quarter),
        ] {
            for _ in 0..count {
                let f = system.step(&[0.001], 0.1, &gate).unwrap();
                totals[0] += f.authored_external_work_j;
                totals[1] += f.inferred_total_dissipation_j;
                totals[2] += f.inlet_work_j;
            }
        }
        same_state(&whole, &halves);
        same_state(&whole, &quarters);
        for totals in [half, quarter] {
            close(frame.authored_external_work_j, totals[0], 2e-15);
            close(frame.inferred_total_dissipation_j, totals[1], 2e-15);
            close(frame.inlet_work_j, totals[2], 2e-15);
        }
        assert!(frame.exact_joint_evolution && frame.duct.ledger.balance_checked);
        covariance = covariance.max(
            (frame.inlet_work_j - frame.mean_pressure_pa * frame.mean_flow_m3_s / 24000.0).abs(),
        );
        close(
            frame.modal.external_work_j,
            frame.authored_external_work_j - frame.inlet_work_j,
            1e-22,
        );
        assert!(frame.balance_residual_j.abs() <= frame.balance_tolerance_j);
    }
    assert!(
        covariance > 1e-12,
        "true continuous work must differ from a product of means: {covariance:e}"
    );
}

fn endpoint(system: &SharedCavityDg) -> Vec<f64> {
    let mut result = Vec::new();
    for s in system
        .cavity()
        .network()
        .components()
        .iter()
        .flat_map(|m| m.states())
    {
        result.extend([
            1400.0 * s.displacement_m_sqrt_kg,
            s.velocity_m_sqrt_kg_per_s,
        ]);
    }
    let gas = system.cavity().source().gas();
    let scale = det::sqrt(
        core::f64::consts::PI * 0.001_f64.powi(2) * 0.05 / (gas.density * gas.sound_speed.powi(2)),
    );
    result.extend(system.load().state().pressure().iter().map(|p| p * scale));
    result
}
#[test]
fn original_mean_connection_refines_toward_the_continuous_joint_operator() {
    let gate = CancelGate::new_clock_free();
    let mut exact = create(48000, false);
    for _ in 0..128 {
        exact.step(&[0.0], 0.0, &gate).unwrap();
    }
    let reference = endpoint(&exact);
    let errors: Vec<_> = [96000, 192000, 384000]
        .into_iter()
        .map(|rate| {
            let mut mean = cavity(rate, false, 1e6)
                .with_downstream_duct(load(rate, false, 0.001), &gate)
                .unwrap();
            for _ in 0..rate / 375 {
                mean.step(&[0.0], 0.0, &gate).unwrap();
            }
            det::sqrt(
                endpoint(&mean)
                    .iter()
                    .zip(&reference)
                    .map(|(x, y)| (x - y).powi(2))
                    .sum(),
            )
        })
        .collect();
    assert!(
        errors[1] < 0.4 * errors[0] && errors[2] < 0.4 * errors[1],
        "mean-port errors: {errors:?}"
    );
}

#[test]
fn actual_thermal_wall_histories_independent_loss_and_atomic_refusal_survive_joint_evolution() {
    let gate = CancelGate::new_clock_free();
    let mut whole = create(48000, true);
    let mut half = create(96000, true);
    let initial = whole.stored_energy_j().unwrap();
    let (mut work, mut loss) = (0.0, 0.0);
    for _ in 0..8 {
        let f = whole.step(&[0.001], 0.03, &gate).unwrap().clone();
        half.step(&[0.001], 0.03, &gate).unwrap();
        half.step(&[0.001], 0.03, &gate).unwrap();
        same_state(&whole, &half);
        work += f.authored_external_work_j;
        loss += f.inferred_total_dissipation_j;
        assert!(f.duct.ledger.balance_checked && f.modal.material_dissipation_j > 0.0);
    }
    assert!(whole.load().state().wall_stored_energy_j() > 0.0);
    assert!(whole.load().state().cumulative_dissipated_j() > 0.0);
    close(
        whole.stored_energy_j().unwrap() + loss,
        initial + work,
        2e-15,
    );
    let frame = whole.last_frame().unwrap().clone();
    let duct = whole.load().state().clone();
    let memory = whole.cavity().network().material_memory().to_vec();
    let state = endpoint(&whole);
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        whole.step(&[0.0], 0.0, &cancelled),
        Err(SharedCavityDgError::Modal(ModalCouplingError::Cancelled))
    ));
    assert!(whole.step(&[f64::NAN], 0.0, &gate).is_err());
    assert!(whole.step(&[1e20], 0.0, &gate).is_err());
    assert_eq!(whole.last_frame().unwrap(), &frame);
    assert_eq!(whole.load().state(), &duct);
    assert_eq!(whole.cavity().network().material_memory(), memory);
    assert_eq!(endpoint(&whole), state);
    assert!(
        cavity(48000, false, 1e6)
            .with_exact_downstream_duct(
                load(48000, false, 0.001),
                ExactDgBudget {
                    max_dense_bytes: 1,
                    ..budget()
                },
                &gate
            )
            .is_err()
    );
    let mut limited = cavity(48000, false, 1e-12)
        .with_exact_downstream_duct(load(48000, false, 0.001), budget(), &gate)
        .unwrap();
    let before = limited.load().state().clone();
    assert!(matches!(
        limited.step(&[0.0], 0.0, &gate),
        Err(SharedCavityDgError::Modal(
            ModalCouplingError::Budget { .. }
        ))
    ));
    assert_eq!(limited.load().state(), &before);
    assert_eq!(limited.cavity().network().samples_rendered(), 0);
}

mod terminal;
