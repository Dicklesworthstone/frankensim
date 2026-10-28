//! Joint propagation and integrated quadratic loss without an unstable inverse
//! propagator. For E(t) = exp(A t) and W(t) = integral E(s)^T C E(s) ds,
//! W(2t) = W(t) + E(t)^T W(t) E(t). Only the physical propagator is squared;
//! a rapidly decaying mode never requires forming exp(-A t).

pub(crate) fn identity(n: usize) -> Vec<f64> {
    let mut a = vec![0.0; n * n];
    for i in 0..n {
        a[i * n + i] = 1.0;
    }
    a
}

pub(crate) fn product(a: &[f64], b: &[f64], n: usize) -> Vec<f64> {
    let mut out = vec![0.0; n * n];
    accumulate(&mut out, a, b, n);
    out
}

// Preserve the existing i-k-j addition order, sharing each row of b between
// four output rows. Both the ordinary exponential and integral use this kernel.
pub(crate) fn accumulate(out: &mut [f64], a: &[f64], b: &[f64], n: usize) {
    if n == 0 {
        return;
    }
    for (a_rows, out_rows) in a.chunks(4 * n).zip(out.chunks_mut(4 * n)) {
        for (k, b_row) in b.chunks_exact(n).enumerate() {
            for (a_row, out_row) in a_rows.chunks_exact(n).zip(out_rows.chunks_exact_mut(n)) {
                let scale = a_row[k];
                for (value, &right) in out_row.iter_mut().zip(b_row) {
                    *value += scale * right;
                }
            }
        }
    }
}

pub(crate) struct FlowIntegral {
    pub(crate) transition: Vec<f64>,
    pub(crate) loss_integral: Vec<f64>,
}

fn largest(a: &[f64]) -> f64 {
    a.iter().fold(0.0_f64, |norm, value| norm.max(value.abs()))
}

fn finite(a: &[f64]) -> bool {
    a.iter().all(|value| value.is_finite())
}

/// The caller admits the dense dimension/work budget before entry. The
/// original 32-term/48-doubling bounds remain. C is a quadratic loss form:
/// its symmetric part has the identical scalar value. No eigenvalue clipping,
/// inferred loss, state removal or change to the requested timestep is used.
pub(crate) fn flow_integral(
    generator: &[f64],
    loss: &[f64],
    n: usize,
    dt: f64,
) -> Result<FlowIntegral, &'static str> {
    if n == 0
        || n.checked_mul(n) != Some(generator.len())
        || loss.len() != generator.len()
        || !finite(generator)
        || !finite(loss)
        || !dt.is_finite()
        || dt < 0.0
    {
        return Err("invalid quadratic flow integral");
    }
    if dt == 0.0 {
        return Ok(FlowIntegral {
            transition: identity(n),
            loss_integral: vec![0.0; n * n],
        });
    }
    // Bound A and A^T so L(X)=A^T X+X A has norm <= 1/2 on the short
    // interval. C's magnitude does not force unnecessary squarings: its
    // series is linear in C and convergence is measured relatively.
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
    let mut doublings = 0;
    while norm * scale > 0.25 {
        if doublings == 48 {
            return Err("quadratic flow scaling exceeds 48 doublings");
        }
        scale *= 0.5;
        doublings += 1;
    }
    let h = dt * scale;
    let a: Vec<_> = generator.iter().map(|value| h * value).collect();
    let mut transition = identity(n);
    let mut e_term = transition.clone();
    let mut loss_integral = vec![0.0; n * n];
    for i in 0..n {
        for j in i..n {
            let value = h * f64::midpoint(loss[i * n + j], loss[j * n + i]);
            loss_integral[i * n + j] = value;
            loss_integral[j * n + i] = value;
        }
    }
    let mut w_term = loss_integral.clone();
    let mut converged = false;
    for order in 1..=32 {
        e_term = product(&e_term, &a, n);
        for (sum, term) in transition.iter_mut().zip(&mut e_term) {
            *term /= f64::from(order);
            *sum += *term;
        }
        // W = h C + L(h C)/2! + L^2(h C)/3! + ... . All terms are
        // symmetric, so one product gives both halves of L(term).
        let wa = product(&w_term, &a, n);
        for i in 0..n {
            for j in i..n {
                let value = (wa[i * n + j] + wa[j * n + i]) / f64::from(order + 1);
                w_term[i * n + j] = value;
                w_term[j * n + i] = value;
                loss_integral[i * n + j] += value;
                if i != j {
                    loss_integral[j * n + i] += value;
                }
            }
        }
        if !finite(&transition) || !finite(&loss_integral) {
            return Err("quadratic flow series is unrepresentable");
        }
        let tolerance = 0.25 * f64::EPSILON;
        if largest(&e_term) <= tolerance * largest(&transition)
            && largest(&w_term) <= tolerance * largest(&loss_integral)
        {
            converged = true;
            break;
        }
    }
    if !converged {
        return Err("quadratic flow series did not converge in 32 terms");
    }
    for _ in 0..doublings {
        let we = product(&loss_integral, &transition, n);
        for i in 0..n {
            for j in i..n {
                let mut upper = 0.0;
                let mut lower = 0.0;
                for k in 0..n {
                    upper += transition[k * n + i] * we[k * n + j];
                    lower += transition[k * n + j] * we[k * n + i];
                }
                // Averaging preserves the quadratic form, not a PSD repair.
                loss_integral[i * n + j] += f64::midpoint(upper, lower);
                loss_integral[j * n + i] = loss_integral[i * n + j];
            }
        }
        transition = product(&transition, &transition, n);
        if !finite(&transition) || !finite(&loss_integral) {
            return Err("quadratic flow doubling is unrepresentable");
        }
    }
    Ok(FlowIntegral {
        transition,
        loss_integral,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(actual: f64, expected: f64, scale: f64) {
        assert!(
            (actual - expected).abs() <= 8.0e-13 * scale,
            "actual={actual:.17e}, expected={expected:.17e}"
        );
    }

    #[test]
    fn passive_decay_resolves_heat_after_a_million_relaxation_times() {
        for rate in [1.0e-8_f64, 0.25, 10.0, 1.0e3, 1.0e6] {
            let flow = flow_integral(&[-rate], &[rate], 1, 1.0).unwrap();
            let decay = (-rate).exp();
            let heat = -0.5 * (-2.0 * rate).exp_m1();
            close(flow.transition[0], decay, decay.max(f64::MIN_POSITIVE));
            close(flow.loss_integral[0], heat, heat);
            close(
                0.5 * flow.transition[0].powi(2) + flow.loss_integral[0],
                0.5,
                0.5,
            );
        }
    }

    #[test]
    fn stiff_driven_storage_closes_independent_source_work_and_heat() {
        // x'=-a*x+u, H=x^2/2. Hold the source and integrate the conjugate
        // output x as additional states: the real exact-ZOH augmentation.
        let rate = 10_000.0;
        let input = 2_500.0;
        let flow = flow_integral(
            &[-rate, input, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
            &[rate, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            3,
            1.0,
        )
        .unwrap();
        let initial = [1.0, 1.0, 0.0];
        let end: Vec<f64> = flow
            .transition
            .chunks_exact(3)
            .map(|row| row.iter().zip(initial).map(|(a, b)| a * b).sum())
            .collect();
        let heat: f64 = (0..3)
            .map(|i| {
                (0..3)
                    .map(|j| initial[i] * flow.loss_integral[3 * i + j] * initial[j])
                    .sum::<f64>()
            })
            .sum();
        close(end[0], 0.25, 1.0);
        close(end[2], 0.25 + 0.75 / rate, 1.0);
        close(heat, 625.65625, 625.65625);
        let source_work = input * end[2];
        close(
            0.5 * (end[0].powi(2) - 1.0) + heat,
            source_work,
            source_work,
        );
    }

    #[test]
    fn nilpotent_transport_retains_zero_modes_and_exact_polynomial_work() {
        for t in [0.0, 0.001, 0.5, 25.0] {
            let flow = flow_integral(&[0.0, 1.0, 0.0, 0.0], &[1.0, 0.0, 0.0, 1.0], 2, t).unwrap();
            for (actual, expected) in flow.transition.iter().zip([1.0, t, 0.0, 1.0]) {
                close(*actual, expected, expected.abs().max(1.0));
            }
            for (actual, expected) in
                flow.loss_integral
                    .iter()
                    .zip([t, t * t / 2.0, t * t / 2.0, t + t * t * t / 3.0])
            {
                close(*actual, expected, expected.abs().max(1.0e-30));
            }
        }
    }

    #[test]
    fn noncommuting_damping_and_rotation_match_the_analytic_energy_integral() {
        let (r1, r2, omega, t) = (2.0_f64, 7.0_f64, 3.0_f64, 0.8_f64);
        let generator = [-r1, omega, -omega, -r2];
        let flow = flow_integral(&generator, &[r1, 0.0, 0.0, r2], 2, t).unwrap();
        let first = flow_integral(&generator, &[r1, 0.0, 0.0, 0.0], 2, t).unwrap();
        let second = flow_integral(&generator, &[0.0, 0.0, 0.0, r2], 2, t).unwrap();
        let mu = 0.5 * (r1 + r2);
        let delta = 0.5 * (r2 - r1);
        let nu = (omega * omega - delta * delta).sqrt();
        let decay = (-mu * t).exp();
        let sine = (nu * t).sin() / nu;
        let cosine = (nu * t).cos();
        let e = [
            decay * (cosine + delta * sine),
            decay * omega * sine,
            -decay * omega * sine,
            decay * (cosine - delta * sine),
        ];
        for i in 0..2 {
            for j in 0..2 {
                close(flow.transition[i * 2 + j], e[i * 2 + j], 1.0);
                // H=x^T x/2 implies W=(I-E^T E)/2; independent of the
                // production Taylor recurrence and its semigroup doubling.
                let identity = if i == j { 1.0 } else { 0.0 };
                let expected = 0.5 * (identity - e[i] * e[j] - e[2 + i] * e[2 + j]);
                close(flow.loss_integral[i * 2 + j], expected, 1.0);
                close(
                    first.loss_integral[i * 2 + j] + second.loss_integral[i * 2 + j],
                    expected,
                    1.0,
                );
            }
        }
        assert_eq!(flow.loss_integral[1], flow.loss_integral[2]);
    }

    #[test]
    fn malformed_or_unrepresentable_integrals_refuse() {
        assert!(flow_integral(&[], &[], 0, 1.0).is_err());
        assert!(flow_integral(&[0.0], &[0.0; 2], 1, 1.0).is_err());
        assert!(flow_integral(&[f64::NAN], &[1.0], 1, 1.0).is_err());
        assert!(flow_integral(&[-1.0], &[f64::INFINITY], 1, 1.0).is_err());
        assert!(flow_integral(&[-1.0], &[1.0], 1, -1.0).is_err());
        assert!(flow_integral(&[-1.0e18], &[1.0], 1, 1.0).is_err());
    }
}
