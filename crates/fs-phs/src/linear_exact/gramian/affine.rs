//! Exact affine support of the existing symmetric Gramian recurrence.
//!
//! A generator with a zero final row and column preserves forms supported only
//! on that coordinate's row/column. Keep their symmetric half-row, with the same
//! ordered products, stopping decisions and ordinary/double-double arithmetic.
//! Forms share Taylor terms of E, but each convergence group is doubled from
//! its OWN stopping-order E immediately. No later E replaces an earlier one.

use super::{
    Compensated, dd_add, dd_div, dd_mul, dd_mul_plain, doubling_rounding_bound, finite, identity,
    largest, product, two_prod,
};

pub(crate) struct AffineFlowIntegral {
    pub(crate) transition: Vec<f64>,
    pub(crate) half_rows: Vec<Vec<f64>>,
    #[cfg(test)]
    pub(crate) convergence_orders: Vec<u32>,
}

/// The caller admits the same dense reference byte/work budget before entry.
/// Every supplied half-row has a zero final entry. None means cancellation;
/// no partial transition or observation is returned on cancellation or refusal.
pub(crate) fn flow_integrals(
    generator: &[f64],
    half_rows: Vec<Vec<f64>>,
    n: usize,
    dt: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Option<AffineFlowIntegral>, &'static str> {
    if cancelled() {
        return Ok(None);
    }
    if n == 0
        || n.checked_mul(n) != Some(generator.len())
        || half_rows.is_empty()
        || !finite(generator)
        || !dt.is_finite()
        || dt < 0.0
        || half_rows.iter().any(|row| row.len() != n || !finite(row))
    {
        return Err("invalid affine flow integral");
    }
    let last = n - 1;
    if (0..n).any(|i| generator[i * n + last] != 0.0 || generator[last * n + i] != 0.0)
        || half_rows.iter().any(|row| row[last] != 0.0)
    {
        return Err("affine flow requires an unchanged constant coordinate");
    }
    if dt == 0.0 {
        let mut rows = half_rows;
        for row in &mut rows {
            if cancelled() {
                return Ok(None);
            }
            row.fill(0.0);
        }
        if cancelled() {
            return Ok(None);
        }
        return Ok(Some(AffineFlowIntegral {
            transition: identity(n),
            #[cfg(test)]
            convergence_orders: vec![0; rows.len()],
            half_rows: rows,
        }));
    }
    // At d<=2, per-form group bookkeeping would consume the old dense byte
    // allowance. Reuse the original tiny full kernel sequentially, retaining
    // compact results in place and the unchanged successful budget admission.
    if n <= 2 {
        let mut rows = half_rows;
        let mut transition = Vec::new();
        for row in &mut rows {
            if cancelled() {
                return Ok(None);
            }
            let mut form = vec![0.0; n * n];
            for i in 0..n {
                form[i * n + last] = row[i];
                form[last * n + i] = row[i];
            }
            let prepared = super::flow_integral(generator, &form, n, dt)?;
            if cancelled() {
                return Ok(None);
            }
            if transition.is_empty() {
                transition = prepared.transition;
            }
            row.copy_from_slice(&prepared.loss_integral[last * n..]);
        }
        if cancelled() {
            return Ok(None);
        }
        return Ok(Some(AffineFlowIntegral {
            transition,
            #[cfg(test)]
            convergence_orders: vec![0; rows.len()],
            half_rows: rows,
        }));
    }
    // Identical scaling and precision selection to the full Gramian owner.
    let mut norm = 0.0_f64;
    for i in 0..n {
        let row: f64 = (0..n).map(|j| (dt * generator[i * n + j]).abs()).sum();
        let column: f64 = (0..n).map(|j| (dt * generator[j * n + i]).abs()).sum();
        norm = norm.max(row).max(column);
    }
    if !norm.is_finite() {
        return Err("nonfinite quadratic flow operator");
    }
    let mut scale = 1.0;
    let mut doublings = 0_u32;
    while norm * scale > 0.25 {
        if doublings == 48 {
            return Err("quadratic flow scaling exceeds 48 doublings");
        }
        scale *= 0.5;
        doublings += 1;
    }
    let h = dt * scale;
    if doubling_rounding_bound(scale, n) > super::super::EXACT_BALANCE_RELATIVE {
        compensated_chain(generator, half_rows, n, h, doublings, cancelled)
    } else {
        plain_chain(generator, half_rows, n, h, doublings, cancelled)
    }
}

fn row_product(out: &mut [f64], row: &[f64], matrix: &[f64], n: usize) {
    out.fill(0.0);
    for k in 0..n {
        for j in 0..n {
            out[j] += row[k] * matrix[k * n + j];
        }
    }
}

fn plain_chain(
    generator: &[f64],
    mut rows: Vec<Vec<f64>>,
    n: usize,
    h: f64,
    doublings: u32,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Option<AffineFlowIntegral>, &'static str> {
    let last = n - 1;
    let a: Vec<_> = generator.iter().map(|value| h * value).collect();
    let mut transition = identity(n);
    let mut e_term = transition.clone();
    for row in &mut rows {
        if cancelled() {
            return Ok(None);
        }
        for value in row {
            *value = h * f64::midpoint(*value, *value);
        }
    }
    let mut terms = Vec::with_capacity(rows.len());
    for row in &rows {
        if cancelled() {
            return Ok(None);
        }
        terms.push(row.clone());
    }
    let mut orders = vec![0_u32; rows.len()];
    let mut first_transition = Vec::new();
    let mut newly_converged = Vec::with_capacity(rows.len());
    let mut completed = 0;
    let mut we_row = vec![0.0; n];
    let mut we_column = vec![0.0; n];
    for order in 1..=32 {
        if cancelled() {
            return Ok(None);
        }
        e_term = product(&e_term, &a, n);
        for (sum, term) in transition.iter_mut().zip(&mut e_term) {
            *term /= f64::from(order);
            *sum += *term;
        }
        if !finite(&transition) {
            return Err("quadratic flow series is unrepresentable");
        }
        let e_converged = largest(&e_term) <= 0.25 * f64::EPSILON * largest(&transition);
        newly_converged.clear();
        for index in 0..rows.len() {
            if cancelled() {
                return Ok(None);
            }
            if orders[index] != 0 {
                continue;
            }
            row_product(&mut we_row, &terms[index], &a, n);
            for j in 0..n {
                // The full symmetric recurrence uses WA[i,last]+WA[last,i].
                let other = if j == last { we_row[j] } else { 0.0 };
                let value = (other + we_row[j]) / f64::from(order + 1);
                terms[index][j] = value;
                rows[index][j] += value;
            }
            if !finite(&rows[index]) {
                return Err("quadratic flow series is unrepresentable");
            }
            if e_converged && largest(&terms[index]) <= 0.25 * f64::EPSILON * largest(&rows[index])
            {
                orders[index] = order;
                newly_converged.push(index);
            }
        }
        if newly_converged.is_empty() {
            continue;
        }
        // Keeping this snapshot only until its group is finished bounds live
        // dense storage independently of the number of convergence groups.
        let mut group_transition = transition.clone();
        for _ in 0..doublings {
            if cancelled() {
                return Ok(None);
            }
            for &index in &newly_converged {
                if cancelled() {
                    return Ok(None);
                }
                row_product(&mut we_row, &rows[index], &group_transition, n);
                for k in 0..last {
                    we_column[k] = 0.0 + rows[index][k] * group_transition[last * n + last];
                }
                we_column[last] = we_row[last];
                for i in 0..n {
                    let mut upper = 0.0;
                    let mut lower = 0.0;
                    for k in 0..n {
                        upper += group_transition[k * n + i] * we_column[k];
                        let entry = if k == last {
                            we_row[i]
                        } else if i == last {
                            we_column[k]
                        } else {
                            0.0
                        };
                        lower += group_transition[k * n + last] * entry;
                    }
                    rows[index][i] += f64::midpoint(upper, lower);
                }
            }
            group_transition = product(&group_transition, &group_transition, n);
            if !finite(&group_transition) {
                return Err("quadratic flow doubling is unrepresentable");
            }
            for &index in &newly_converged {
                if cancelled() {
                    return Ok(None);
                }
                if !finite(&rows[index]) {
                    return Err("quadratic flow doubling is unrepresentable");
                }
            }
        }
        if orders[0] == order {
            first_transition = group_transition;
        }
        completed += newly_converged.len();
        if completed == rows.len() {
            if cancelled() {
                return Ok(None);
            }
            return Ok(Some(AffineFlowIntegral {
                transition: first_transition,
                half_rows: rows,
                #[cfg(test)]
                convergence_orders: orders,
            }));
        }
    }
    Err("quadratic flow series did not converge in 32 terms")
}

fn compensated_row_times_plain(out: &mut Compensated, row: &Compensated, matrix: &[f64], n: usize) {
    out.0.fill((0.0, 0.0));
    for k in 0..n {
        for j in 0..n {
            out.0[j] = dd_add(out.0[j], dd_mul_plain(row.0[k], matrix[k * n + j]));
        }
    }
}

fn compensated_row_product(
    out: &mut Compensated,
    row: &Compensated,
    matrix: &Compensated,
    n: usize,
) {
    out.0.fill((0.0, 0.0));
    for k in 0..n {
        for j in 0..n {
            out.0[j] = dd_add(out.0[j], dd_mul(row.0[k], matrix.0[k * n + j]));
        }
    }
}

fn compensated_chain(
    generator: &[f64],
    half_rows: Vec<Vec<f64>>,
    n: usize,
    h: f64,
    doublings: u32,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Option<AffineFlowIntegral>, &'static str> {
    let last = n - 1;
    let a: Vec<_> = generator.iter().map(|value| h * value).collect();
    let mut transition = Compensated::identity(n);
    let mut e_term = transition.clone();
    let mut rows = Vec::with_capacity(half_rows.len());
    for row in half_rows {
        if cancelled() {
            return Ok(None);
        }
        rows.push(Compensated(
            row.into_iter()
                .map(|value| two_prod(h, f64::midpoint(value, value)))
                .collect(),
        ));
    }
    let mut terms = Vec::with_capacity(rows.len());
    for row in &rows {
        if cancelled() {
            return Ok(None);
        }
        terms.push(row.clone());
    }
    let mut orders = vec![0_u32; rows.len()];
    let mut first_transition = Vec::new();
    let mut newly_converged = Vec::with_capacity(rows.len());
    let mut completed = 0;
    let mut we_row = Compensated(vec![(0.0, 0.0); n]);
    let mut we_column = vec![(0.0, 0.0); n];
    let tolerance = 0.25 * f64::EPSILON * f64::EPSILON;
    for order in 1..=32 {
        if cancelled() {
            return Ok(None);
        }
        e_term = e_term.times_plain(&a, n);
        for (sum, term) in transition.0.iter_mut().zip(&mut e_term.0) {
            *term = dd_div(*term, f64::from(order));
            *sum = dd_add(*sum, *term);
        }
        if !transition.finite() {
            return Err("quadratic flow series is unrepresentable");
        }
        let e_converged = e_term.largest() <= tolerance * transition.largest();
        newly_converged.clear();
        for index in 0..rows.len() {
            if cancelled() {
                return Ok(None);
            }
            if orders[index] != 0 {
                continue;
            }
            compensated_row_times_plain(&mut we_row, &terms[index], &a, n);
            for j in 0..n {
                let other = if j == last { we_row.0[j] } else { (0.0, 0.0) };
                let value = dd_div(dd_add(other, we_row.0[j]), f64::from(order + 1));
                terms[index].0[j] = value;
                rows[index].0[j] = dd_add(rows[index].0[j], value);
            }
            if !rows[index].finite() {
                return Err("quadratic flow series is unrepresentable");
            }
            if e_converged && terms[index].largest() <= tolerance * rows[index].largest() {
                orders[index] = order;
                newly_converged.push(index);
            }
        }
        if newly_converged.is_empty() {
            continue;
        }
        let mut group_transition = transition.clone();
        for _ in 0..doublings {
            if cancelled() {
                return Ok(None);
            }
            for &index in &newly_converged {
                if cancelled() {
                    return Ok(None);
                }
                compensated_row_product(&mut we_row, &rows[index], &group_transition, n);
                for k in 0..last {
                    we_column[k] = dd_add(
                        (0.0, 0.0),
                        dd_mul(rows[index].0[k], group_transition.0[last * n + last]),
                    );
                }
                we_column[last] = we_row.0[last];
                for i in 0..n {
                    let mut upper = (0.0, 0.0);
                    let mut lower = (0.0, 0.0);
                    for k in 0..n {
                        upper = dd_add(upper, dd_mul(group_transition.0[k * n + i], we_column[k]));
                        let entry = if k == last {
                            we_row.0[i]
                        } else if i == last {
                            we_column[k]
                        } else {
                            (0.0, 0.0)
                        };
                        lower = dd_add(lower, dd_mul(group_transition.0[k * n + last], entry));
                    }
                    // Preserve both differently ordered DD dot products and
                    // the original midpoint, including its low-part rounding.
                    let (sh, sl) = dd_add(upper, lower);
                    rows[index].0[i] = dd_add(rows[index].0[i], (0.5 * sh, 0.5 * sl));
                }
            }
            group_transition = group_transition.product(&group_transition, n);
            if !group_transition.finite() {
                return Err("quadratic flow doubling is unrepresentable");
            }
            for &index in &newly_converged {
                if cancelled() {
                    return Ok(None);
                }
                if !rows[index].finite() {
                    return Err("quadratic flow doubling is unrepresentable");
                }
            }
        }
        if orders[0] == order {
            first_transition = group_transition.high();
        }
        completed += newly_converged.len();
        if completed == rows.len() {
            if cancelled() {
                return Ok(None);
            }
            let mut half_rows = Vec::with_capacity(rows.len());
            for row in rows {
                if cancelled() {
                    return Ok(None);
                }
                half_rows.push(row.high());
            }
            if cancelled() {
                return Ok(None);
            }
            return Ok(Some(AffineFlowIntegral {
                transition: first_transition,
                half_rows,
                #[cfg(test)]
                convergence_orders: orders,
            }));
        }
    }
    Err("quadratic flow series did not converge in 32 terms")
}

#[cfg(test)]
mod tests {
    use super::*;

    fn bits(values: &[f64]) -> Vec<u64> {
        values.iter().map(|value| value.to_bits()).collect()
    }

    fn dense(row: &[f64]) -> Vec<f64> {
        let n = row.len();
        let last = n - 1;
        let mut form = vec![0.0; n * n];
        for i in 0..n {
            form[i * n + last] = row[i];
            form[last * n + i] = row[i];
        }
        form
    }

    fn parity(generator: &[f64], rows: &[Vec<f64>], n: usize, dt: f64) -> AffineFlowIntegral {
        let batch = flow_integrals(generator, rows.to_vec(), n, dt, &mut || false)
            .unwrap()
            .unwrap();
        for (index, row) in rows.iter().enumerate() {
            let reference = super::super::flow_integral(generator, &dense(row), n, dt).unwrap();
            if index == 0 {
                assert_eq!(bits(&batch.transition), bits(&reference.transition));
            }
            assert_eq!(
                bits(&batch.half_rows[index]),
                bits(&reference.loss_integral[(n - 1) * n..]),
                "row {index}, dt={dt}",
            );
            for i in 0..n {
                assert_eq!(
                    batch.half_rows[index][i].to_bits(),
                    reference.loss_integral[i * n + n - 1].to_bits(),
                );
            }
        }
        assert_eq!(
            batch.half_rows.iter().map(Vec::len).sum::<usize>(),
            rows.len() * n
        );
        batch
    }

    #[test]
    fn ordinary_affine_batch_matches_every_dense_bit_and_keeps_compact_payload() {
        let generator = [
            -0.4, 3.0, 0.2, 1.0, 0.0, -3.0, -0.1, 0.5, -0.4, 0.0, 0.0, 0.0, -80.0, 0.3, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        ];
        let tiny = f64::from_bits(1);
        let rows = vec![
            vec![0.3, -0.2, 0.0, 0.0, 0.0],
            vec![-0.7, 0.4, 1.25, -0.1, 0.0],
            vec![0.0; 5],
            vec![-0.0, 0.0, -0.0, 0.0, 0.0],
            vec![tiny, -tiny, tiny, 0.0, 0.0],
            vec![1e-200, -2e-200, 3e-200, 0.0, 0.0],
            vec![1e200, -2e200, 3e200, 0.0, 0.0],
        ];
        for dt in [0.0, 0.001, 0.013] {
            parity(&generator, &rows, 5, dt);
        }
    }

    #[test]
    fn compensated_affine_batch_preserves_both_dd_dot_orders_and_final_rounding() {
        let omega = 1_048_576.0;
        let generator = [
            -1.0, omega, 0.25, 0.0, -omega, -1.0, -0.7, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        ];
        let rows = vec![
            vec![0.3, -0.2, 0.0, 0.0],
            vec![-0.7, 0.4, 0.1, 0.0],
            vec![0.0; 4],
            vec![-0.0, 0.0, -0.0, 0.0],
            vec![1e-200, -2e-200, 0.0, 0.0],
            vec![1e200, -2e200, 0.0, 0.0],
        ];
        assert!(
            doubling_rounding_bound(0.5_f64.powi(23), 4)
                > super::super::super::EXACT_BALANCE_RELATIVE
        );
        parity(&generator, &rows, 4, 1.0);
    }

    #[test]
    fn each_convergence_group_uses_its_own_transition_in_plain_and_dd_chains() {
        // For a dense rank-one decay, individual E entries have a 1/32 factor
        // that the all-ones observation cancels. Its W series therefore needs
        // one more term than the zero/coordinate observations. Put that slower
        // form first to check that the earliest completed group's E cannot win.
        let physical = 32;
        let n = physical + 1;
        let mut rows = vec![vec![0.0; n]; 3];
        rows[0][..physical].fill(1.0);
        rows[2][0] = 1.0;
        for rate in [0.21875, 262_144.0] {
            let mut generator = vec![0.0; n * n];
            for i in 0..physical {
                for j in 0..physical {
                    generator[i * n + j] = -rate / physical as f64;
                }
            }
            let batch = parity(&generator, &rows, n, 1.0);
            assert!(
                batch.convergence_orders[0] > batch.convergence_orders[1],
                "rate={rate}, orders={:?}",
                batch.convergence_orders,
            );
            if rate == 0.21875 {
                let earlier =
                    super::super::flow_integral(&generator, &dense(&rows[1]), n, 1.0).unwrap();
                assert_ne!(bits(&batch.transition), bits(&earlier.transition));
            }
        }
    }

    #[test]
    fn tiny_dimension_fallback_keeps_many_outputs_with_no_group_bookkeeping() {
        let generator = [-2.0, 0.0, 0.0, 0.0];
        let rows: Vec<_> = (0..128)
            .map(|i| vec![(i as f64 - 63.0) / 17.0, 0.0])
            .collect();
        let batch = parity(&generator, &rows, 2, 0.1);
        assert!(batch.convergence_orders.iter().all(|order| *order == 0));
    }

    #[test]
    fn cancellation_publishes_no_partial_batch_and_does_not_change_replay() {
        let generator = [-0.4, 3.0, 0.0, -3.0, -0.1, 0.0, 0.0, 0.0, 0.0];
        let rows = vec![vec![0.3, -0.2, 0.0], vec![-0.7, 0.4, 0.0], vec![0.0; 3]];
        let mut polls = 0;
        let baseline = flow_integrals(&generator, rows.clone(), 3, 0.8, &mut || {
            polls += 1;
            false
        })
        .unwrap()
        .unwrap();
        for stop in [1, polls / 3, polls - 1] {
            let mut seen = 0;
            assert!(
                flow_integrals(&generator, rows.clone(), 3, 0.8, &mut || {
                    seen += 1;
                    seen >= stop
                })
                .unwrap()
                .is_none()
            );
        }
        let replay = flow_integrals(&generator, rows, 3, 0.8, &mut || false)
            .unwrap()
            .unwrap();
        assert_eq!(bits(&baseline.transition), bits(&replay.transition));
        for (a, b) in baseline.half_rows.iter().zip(&replay.half_rows) {
            assert_eq!(bits(a), bits(b));
        }
    }

    #[test]
    fn affine_support_is_exact_and_invalid_operators_still_refuse() {
        let mut generator = vec![0.0; 9];
        generator[0] = -1.0;
        let rows = vec![vec![0.5, 0.0, 0.0]];
        assert!(flow_integrals(&generator, rows.clone(), 3, -1.0, &mut || false).is_err());
        assert!(flow_integrals(&generator, Vec::new(), 3, 1.0, &mut || false).is_err());
        generator[2] = f64::from_bits(1);
        assert!(flow_integrals(&generator, rows.clone(), 3, 1.0, &mut || false).is_err());
        generator[2] = 0.0;
        generator[6] = f64::from_bits(1);
        assert!(flow_integrals(&generator, rows.clone(), 3, 1.0, &mut || false).is_err());
        generator[6] = 0.0;
        generator[0] = f64::NAN;
        assert!(flow_integrals(&generator, rows, 3, 1.0, &mut || false).is_err());
        assert!(
            flow_integrals(
                &[0.0; 9],
                vec![vec![0.0, 0.0, f64::from_bits(1)]],
                3,
                1.0,
                &mut || false
            )
            .is_err()
        );
    }
}
