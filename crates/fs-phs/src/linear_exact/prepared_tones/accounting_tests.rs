use super::*;
use crate::QuadraticStorage;

fn fixture() -> (PortHamiltonian, [ResistanceWitness; 2]) {
    let system = PortHamiltonian::new(
        1,
        2,
        vec![0.0],
        vec![2.0],
        vec![1.0, -1.0],
        Box::new(QuadraticStorage::new(vec![1.0], 1).unwrap()),
    )
    .unwrap()
    .with_resistance_witness(
        ResistanceWitness::try_from_source_rows(1, vec![2.0], vec![det::sqrt(2.0)]).unwrap(),
    )
    .unwrap();
    let channels = [0.5, 1.5]
        .map(|r| ResistanceWitness::try_from_source_rows(1, vec![r], vec![det::sqrt(r)]).unwrap());
    (system, channels)
}

fn budget() -> DenseExponentialBudget {
    DenseExponentialBudget {
        max_dimension: 64,
        max_dense_bytes: 16 * 1024 * 1024,
        max_multiply_terms: 100_000_000,
    }
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 2e-11,
        "{actual:e} != {expected:e}"
    );
}

// G1: independent scalar integrating-factor solution of x'=-2x+u0-u1,
// followed by composite Simpson integration of EACH physical power. This
// reference does not use the native exponential, Gramian or source-work cache.
fn independent_accounts() -> [f64; 5] {
    let dt = 0.7;
    let panels = 4096;
    let h = dt / panels as f64;
    let mut integrals = [0.0; 4];
    let mut endpoint = 0.0;
    for i in 0..=panels {
        let t = i as f64 * h;
        let (sn, cs) = (3.0 * t).sin_cos();
        let decay = (-2.0 * t).exp();
        let cosine = (2.0 * cs + 3.0 * sn - 2.0 * decay) / 13.0;
        let sine = (2.0 * sn - 3.0 * cs + 3.0 * decay) / 13.0;
        let x = 0.3 * decay + 0.75 * (1.0 - decay) + 1.15 * cosine - 1.1 * sine;
        let u0 = 2.0 + 0.75 * cs - 0.2 * sn;
        let u1 = 0.5 - 0.4 * cs + 0.9 * sn;
        let weight = if i == 0 || i == panels {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        for (sum, value) in integrals.iter_mut().zip([u0 * x, -u1 * x, x * x, x]) {
            *sum += weight * value * h / 3.0;
        }
        endpoint = x;
    }
    [
        endpoint,
        integrals[0],
        integrals[1],
        integrals[2],
        integrals[3] / dt,
    ]
}

#[test]
fn borrowed_channels_and_port_work_match_independent_physical_integrals() {
    let (system, channels) = fixture();
    let mut prepared = PreparedExactTones::new_with_resistance_channels_controlled(
        &system,
        &[0.0, 3.0],
        0.7,
        LinearSplitBudget {
            max_modes: 1,
            max_omega_dt: 1.0,
        },
        budget(),
        &channels,
        || false,
    )
    .unwrap();
    assert!(std::ptr::eq(prepared.system(), &system));
    let mut next = [0.0];
    let mut output = [0.0; 2];
    let mut heat = [0.0; 2];
    let mut work = [0.0; 2];
    let record = prepared
        .step_accounted_into(
            &[0.3],
            &[2.0, 0.75, 0.5, -0.4],
            &[0.0, -0.2, 0.0, 0.9],
            &mut next,
            &mut output,
            &mut heat,
            &mut work,
        )
        .unwrap();
    let reference = independent_accounts();
    near(next[0], reference[0]);
    near(work[0], reference[1]);
    near(work[1], reference[2]);
    near(heat[0], 0.5 * reference[3]);
    near(heat[1], 1.5 * reference[3]);
    near(output[0], reference[4]);
    near(output[1], -reference[4]);
    near(work[0] + work[1], record.supplied);
    near(heat[0] + heat[1], record.dissipated);
    near(record.delta_h + record.dissipated, record.supplied);
}

#[test]
fn opposing_continuous_sources_retain_signed_exchange_when_net_supply_is_zero() {
    let (system, _) = fixture();
    let mut prepared = PreparedExactTones::new(
        &system,
        &[0.0, 3.0],
        0.7,
        LinearSplitBudget {
            max_modes: 1,
            max_omega_dt: 1.0,
        },
        budget(),
    )
    .unwrap();
    let mut next = [0.0];
    let mut output = [0.0; 2];
    let mut work = [0.0; 2];
    let record = prepared
        .step_accounted_into(
            &[0.3],
            &[2.0, 0.75, 2.0, 0.75],
            &[0.0, -0.2, 0.0, -0.2],
            &mut next,
            &mut output,
            &mut [],
            &mut work,
        )
        .unwrap();
    near(next[0], 0.3 * (-1.4_f64).exp());
    assert!(work[0] > 0.01 && work[1] < -0.01);
    near(work[0], -work[1]);
    near(record.supplied, 0.0);
    assert!(record.dissipated > 0.0);
}

#[test]
fn work_accounting_refusals_publish_no_outputs_and_preserve_retry() {
    let (system, channels) = fixture();
    let mut prepared = PreparedExactTones::new_with_resistance_channels_controlled(
        &system,
        &[0.0, 3.0],
        0.7,
        LinearSplitBudget {
            max_modes: 1,
            max_omega_dt: 1.0,
        },
        budget(),
        &channels,
        || false,
    )
    .unwrap();
    let cosine = [2.0, 0.75, 0.5, -0.4];
    let sine = [0.0, -0.2, 0.0, 0.9];
    let mut next = [91.0];
    let mut output = [92.0; 2];
    let mut heat = [93.0; 2];
    let mut work = [94.0; 2];
    assert!(
        prepared
            .step_accounted_into(
                &[0.3],
                &cosine,
                &sine,
                &mut next,
                &mut output,
                &mut heat,
                &mut work[..1],
            )
            .is_err()
    );
    let mut checks = 0;
    let refused = prepared.step_accounted_into_controlled(
        &[0.3],
        &cosine,
        &sine,
        &mut next,
        &mut output,
        &mut heat,
        &mut work,
        || {
            checks += 1;
            checks == 4
        },
    );
    assert!(matches!(refused, Err(PreparedStepError::Cancelled)));
    assert_eq!(next, [91.0]);
    assert_eq!(output, [92.0; 2]);
    assert_eq!(heat, [93.0; 2]);
    assert_eq!(work, [94.0; 2]);
    prepared
        .step_accounted_into(
            &[0.3],
            &cosine,
            &sine,
            &mut next,
            &mut output,
            &mut heat,
            &mut work,
        )
        .unwrap();
    let reference = independent_accounts();
    near(next[0], reference[0]);
    near(work[0], reference[1]);
    near(work[1], reference[2]);
}
