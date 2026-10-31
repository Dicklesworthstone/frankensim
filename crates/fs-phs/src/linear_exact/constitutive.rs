//! Preserve the original finite-step law under a declared quadratic change.
//! Dense rows mirror QuadraticStorage's ordered sparse mean/action arithmetic;
//! no fitted law, extra state or callback allocation is introduced.

pub(super) fn add_quadratic_mean(
    q: &[f64],
    a: &[f64],
    b: &[f64],
    mid: &[f64],
    sign: f64,
    out: &mut [f64],
) {
    let n = a.len();
    for i in 0..n {
        let (mut gradient, mut correction) = (0.0, 0.0);
        for j in 0..n {
            let qij = q[i * n + j];
            let qji = q[j * n + i];
            if qij == 0.0 && qji == 0.0 {
                continue;
            }
            gradient += qij * mid[j];
            let offset = (a[j] - mid[j]) + (b[j] - mid[j]);
            correction += (0.25 * qij + 0.25 * qji) * offset + (0.5 * qji - 0.5 * qij) * mid[j];
        }
        // Keep the SAME two sums and final addition as the native quadratic
        // owner so subtracting an original block leaves exactly zero effort.
        out[i] += sign * (gradient + correction);
    }
}

pub(super) fn add_quadratic_end_action(q: &[f64], direction: &[f64], sign: f64, out: &mut [f64]) {
    let n = direction.len();
    for i in 0..n {
        let mut action = 0.0;
        for j in 0..n {
            let qij = q[i * n + j];
            let qji = q[j * n + i];
            if qij != 0.0 || qji != 0.0 {
                action += (0.25 * qij + 0.25 * qji) * direction[j];
            }
        }
        out[i] += sign * action;
    }
}

#[cfg(test)]
pub(super) fn source_fixture(
    scales: [f64; 2],
) -> (
    crate::circuit::CircuitDae,
    [crate::circuit::capacitance::CapacitorChargeCurve; 2],
) {
    use crate::circuit::capacitance::{
        CapacitorChargeCurve, CapacitorCurve, assemble_circuit_with_capacitor_curves,
    };
    use crate::circuit::{Branch, CircuitGraph, MutualInductance};
    let knots = [
        [
            (-4.0, -7.0),
            (-1.0, -1.0),
            (0.0, 0.0),
            (1.0, 1.0),
            (4.0, 7.0),
        ],
        [
            (-4.0, -5.0),
            (-0.5, -1.5),
            (0.0, 0.0),
            (0.5, 1.5),
            (4.0, 5.0),
        ],
    ];
    let laws = std::array::from_fn(|i| {
        CapacitorChargeCurve::try_new(knots[i].iter().map(|&(q, v)| (q * scales[i], v)).collect())
            .unwrap()
    });
    let dae = assemble_circuit_with_capacitor_curves(
        &CircuitGraph {
            node_count: 2,
            branches: vec![
                (1, 0, Branch::Capacitor { farads: 1.0 }),
                (1, 0, Branch::Capacitor { farads: 1.0 }),
                (1, 0, Branch::Inductor { henries: 2.0 }),
                (1, 0, Branch::Inductor { henries: 3.0 }),
                (1, 0, Branch::Capacitor { farads: 1.0 }),
            ],
            transformers: vec![],
        },
        &[
            CapacitorCurve {
                branch: 0,
                curve: laws[0].clone(),
            },
            CapacitorCurve {
                branch: 1,
                curve: laws[1].clone(),
            },
        ],
        &[MutualInductance {
            branches: (2, 3),
            henries: 0.4,
        }],
    )
    .unwrap();
    (dae, laws)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::linear_exact::{
        LinearPropagator, LinearSplitBudget, MatchedStorage, split_quadratic, step_split,
    };
    use crate::{Storage, discrete_gradient};

    #[test]
    fn matched_wrapper_keeps_original_component_work_and_actual_end_derivative() {
        let (dae, laws) = source_fixture([1.0; 2]);
        let n = dae.system.state_dim();
        let original = dae.system.storage();
        let mut correction = original.quadratic_part().unwrap().matrix;
        for value in &mut correction {
            *value *= 0.25;
        }
        let storage = MatchedStorage {
            original,
            correction: &correction,
        };
        let (mut a, mut b) = (vec![0.0; n], vec![0.0; n]);
        a[dae.charge_index[0]] = 2.0;
        b[dae.charge_index[0]] = 0.5;
        a[dae.charge_index[1]] = 0.25;
        b[dae.charge_index[1]] = 1.5;
        a[dae.charge_index[2]] = 0.7;
        b[dae.charge_index[2]] = -0.4;
        a[dae.flux_index[0]] = 0.2;
        b[dae.flux_index[0]] = 0.3;
        a[dae.flux_index[1]] = -0.1;
        b[dae.flux_index[1]] = 0.2;
        let effort = discrete_gradient(&storage, &a, &b);
        for (slot, law) in dae.charge_index.iter().take(2).zip(&laws) {
            assert_eq!(
                effort[*slot],
                law.discrete_voltage_v(a[*slot], b[*slot]).unwrap()
            );
            assert!(
                (effort[*slot] * (b[*slot] - a[*slot])
                    - law.energy_change_j(a[*slot], b[*slot]).unwrap())
                .abs()
                    < 1e-14
            );
        }
        let work: f64 = (0..n).map(|i| effort[i] * (b[i] - a[i])).sum();
        assert!((work - (storage.hamiltonian(&b) - storage.hamiltonian(&a))).abs() < 1e-13);
        let direction: Vec<_> = (0..n).map(|i| 0.03 * (i as f64 - 2.0)).collect();
        let mut action = vec![0.0; n];
        assert!(
            storage
                .discrete_gradient_rule()
                .unwrap()
                .end_jacobian_action(&a, &b, &direction, &mut action)
        );
        let (mut hi, mut lo) = (b.clone(), b.clone());
        for i in 0..n {
            hi[i] += 1e-5 * direction[i];
            lo[i] -= 1e-5 * direction[i];
        }
        let gp = discrete_gradient(&storage, &a, &hi);
        let gm = discrete_gradient(&storage, &a, &lo);
        for i in 0..n {
            assert!((action[i] - (gp[i] - gm[i]) / 2e-5).abs() < 2e-10);
        }
        let mut outside = b;
        outside[dae.charge_index[0]] = 5.0;
        assert!(
            discrete_gradient(&storage, &a, &outside)
                .iter()
                .any(|v| !v.is_finite())
        );
        assert!(
            !storage
                .discrete_gradient_rule()
                .unwrap()
                .end_jacobian_action(&a, &outside, &direction, &mut action)
        );
        let quadratic = crate::QuadraticStorage::new(vec![1.0], 1).unwrap();
        assert!(
            MatchedStorage {
                original: &quadratic,
                correction: &[0.2]
            }
            .discrete_gradient_rule()
            .is_none()
        );
    }

    #[test]
    fn matched_mode_resolves_both_small_charge_units_beside_a_large_held_store() {
        let scales = [1e-9, 1e-12];
        let (dae, laws) = source_fixture(scales);
        let n = dae.system.state_dim();
        let mut ports = vec![0.0; 2 * n];
        ports[2 * dae.charge_index[0]] = 1.0;
        ports[2 * dae.charge_index[1] + 1] = 1.0;
        let system = crate::PortHamiltonian::new(
            n,
            2,
            vec![0.0; n * n],
            vec![0.0; n * n],
            ports,
            dae.system.storage().shared_clone().unwrap(),
        )
        .unwrap();
        let split = split_quadratic(
            &system,
            LinearSplitBudget {
                max_modes: n,
                max_omega_dt: 1.0,
            },
        )
        .unwrap();
        let mut initial = vec![0.0; n];
        initial[dae.charge_index[0]] = 2.0 * scales[0];
        initial[dae.charge_index[1]] = 0.25 * scales[1];
        initial[dae.charge_index[2]] = 1e8;
        let input = [-1.5 * scales[0], 1.25 * scales[1]];
        let record = step_split(
            &system,
            &split,
            LinearPropagator::FrequencyMatched,
            &initial,
            &input,
            1.0,
        )
        .unwrap();
        assert!(record.newton_iters > 0);
        let mut total = 0.0;
        for i in 0..2 {
            let slot = dae.charge_index[i];
            let expected = initial[slot] + input[i];
            assert!((record.x[slot] - expected).abs() < 1e-13 * scales[i]);
            let work = laws[i].energy_change_j(initial[slot], expected).unwrap();
            assert!((record.y[i] * input[i] - work).abs() < 3e-13 * scales[i]);
            total += work;
        }
        assert_eq!(record.x[dae.charge_index[2]], 1e8);
        assert_eq!(system.hamiltonian(&initial), system.hamiltonian(&record.x));
        assert!((record.delta_h - total).abs() < 1e-22);
        assert!((record.supplied - total).abs() < 1e-22);
        assert_eq!(record.dissipated, 0.0);
    }
}
