//! Exact discrete RK4 polynomial oracles, independent of stage integration.

use super::*;

struct Linear {
    dimension: usize,
    matrix: Vec<f64>,
}

impl OrbitProblem for Linear {
    fn dim(&self) -> usize {
        self.dimension
    }
    fn island(&self, _: f64, state: &[f64], out: &mut [f64]) {
        for (i, target) in out.iter_mut().enumerate() {
            *target = state
                .iter()
                .enumerate()
                .map(|(j, x)| self.matrix[i * self.dimension + j] * x)
                .sum();
        }
    }
    fn autonomous(&self) -> bool {
        true
    }
    fn time_domain_rhs(&self) -> Option<TimeDomainRhs<Self>> {
        Some(Self::island)
    }
    fn time_domain_state_jacobian(&self) -> Option<TimeDomainStateJacobian<Self>> {
        Some(|problem, _, _, out| {
            out.copy_from_slice(&problem.matrix);
            true
        })
    }
}

// One closed polynomial evaluation, with no RK stages or numerical probes.
fn stability(z: f64) -> (f64, f64) {
    let square = z * z;
    let cube = square * z;
    (
        1.0 + z + square / 2.0 + cube / 6.0 + square * square / 24.0,
        1.0 + z + square / 2.0 + cube / 6.0,
    )
}

fn period_map(rate: f64, period: f64, steps: usize) -> (f64, f64) {
    let (r, derivative) = stability(rate * period / steps as f64);
    // d R(a*T/N)^N / dT: the factor N cancels d(a*T/N)/dT.
    (
        r.powi(steps as i32),
        rate * derivative * r.powi(steps as i32 - 1),
    )
}

fn close(actual: f64, expected: f64, steps: usize) {
    let tolerance = 64.0 * f64::EPSILON * (steps + 1) as f64;
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1.0),
        "{steps} steps: {actual:.17e} != {expected:.17e}"
    );
}

#[allow(clippy::too_many_arguments)] // Complete original physical chart and discrete period map.
fn integrate(
    problem: &Linear,
    initial: &[f64],
    period: f64,
    steps: usize,
    scales: &[f64],
    period_scale: f64,
) -> Image {
    let mut workspace = Workspace::new(problem.dimension, problem.dimension + 1);
    let mut counts = Counts {
        maps: 0,
        maximum: 1,
        rhs: 0,
        jacobian: 0,
    };
    let image = workspace
        .integrate(
            problem,
            problem.time_domain_rhs().unwrap(),
            problem.time_domain_state_jacobian(),
            initial,
            period,
            steps,
            scales,
            period_scale,
            &mut counts,
            &mut || false,
        )
        .unwrap();
    assert_eq!(
        (counts.maps, counts.rhs, counts.jacobian),
        (1, 4 * steps, 4 * steps)
    );

    // The variational path must differentiate exactly the primal-only map.
    let mut primal = Workspace::new(problem.dimension, 0);
    let mut counts = Counts {
        maps: 0,
        maximum: 1,
        rhs: 0,
        jacobian: 0,
    };
    let primal = primal
        .integrate(
            problem,
            problem.time_domain_rhs().unwrap(),
            None,
            initial,
            period,
            steps,
            scales,
            period_scale,
            &mut counts,
            &mut || false,
        )
        .unwrap();
    assert_eq!(primal.end, image.end);
    assert_eq!(primal.excursion, image.excursion);
    assert!(primal.tangent.is_empty());
    assert_eq!(
        (counts.maps, counts.rhs, counts.jacobian),
        (1, 4 * steps, 0)
    );
    image
}

#[test]
fn scalar_state_and_period_columns_are_the_discrete_stability_polynomial_derivatives() {
    for rate in [-1.0, 0.0, 0.75] {
        let problem = Linear {
            dimension: 1,
            matrix: vec![rate],
        };
        for initial in [-0.8, 1.7] {
            for steps in [1, 2, 7, 32, 128] {
                let (period, period_scale) = (1.4, 0.6);
                let image = integrate(&problem, &[initial], period, steps, &[1.0], period_scale);
                let (map, derivative) = period_map(rate, period, steps);
                close(image.end[0], map * initial, steps);
                close(image.tangent[0], map, steps);
                close(image.tangent[1], derivative * initial * period_scale, steps);
                if rate == 0.0 {
                    assert_eq!(image.end[0], initial);
                    assert_eq!(image.tangent, vec![1.0, 0.0]);
                }
            }
        }
    }
    // A coarse RK4 period derivative is observably NOT f(Phi_T(x)).
    // R(-1)=3/8 and R'(-1)=1/3, whereas exp(-1) has equal value/derivative.
    let problem = Linear {
        dimension: 1,
        matrix: vec![-1.0],
    };
    let image = integrate(&problem, &[1.0], 1.0, 1, &[1.0], 1.0);
    close(image.end[0], 3.0 / 8.0, 1);
    close(image.tangent[1], -1.0 / 3.0, 1);
    assert!((image.tangent[1] - (-image.end[0])).abs() > 0.04);
}

#[test]
fn scalar_period_column_preserves_mixed_physical_state_and_time_units() {
    let (rate, period, period_scale, initial) = (-0.7, 1.3, 0.6, 1.4);
    for state_exponent in [-300, 0, 300] {
        let state_unit = 2.0_f64.powi(state_exponent);
        for time_exponent in [-400, 0, 400] {
            let time_unit = 2.0_f64.powi(time_exponent);
            let problem = Linear {
                dimension: 1,
                matrix: vec![rate / time_unit],
            };
            for steps in [1, 9, 128] {
                let image = integrate(
                    &problem,
                    &[initial],
                    period * time_unit,
                    steps,
                    &[state_unit],
                    period_scale * time_unit,
                );
                let (map, derivative) = period_map(rate, period, steps);
                close(image.end[0], map * initial, steps);
                close(image.tangent[0], map, steps);
                close(image.tangent[1], derivative * initial * period_scale, steps);
                // Convert the normalized period column back to physical dx/dT,
                // then recover the same base-unit result without squaring units.
                let physical = (image.tangent[1] / (period_scale * time_unit)) * state_unit;
                close(
                    (physical / state_unit) * time_unit,
                    derivative * initial,
                    steps,
                );
            }
        }
    }
}

#[test]
fn nonnormal_triangular_state_and_period_columns_retain_all_physical_coupling() {
    let (a, d, coupling, period, period_scale) = (-0.5, 0.3, 32.0, 1.1, 0.6);
    let initial = [0.4, -0.7];
    for state_units in [[1.0, 1.0], [2.0_f64.powi(120), 2.0_f64.powi(-120)]] {
        for time_exponent in [-200, 0, 200] {
            let time_unit = 2.0_f64.powi(time_exponent);
            let problem = Linear {
                dimension: 2,
                matrix: vec![
                    a / time_unit,
                    ((coupling * state_units[0]) / state_units[1]) / time_unit,
                    0.0,
                    d / time_unit,
                ],
            };
            for steps in [1, 3, 32] {
                let image = integrate(
                    &problem,
                    &initial,
                    period * time_unit,
                    steps,
                    &state_units,
                    period_scale * time_unit,
                );
                // For distinct triangular eigenvalues, every analytic matrix
                // function has off-diagonal b*(g(a)-g(d))/(a-d). Applying this
                // to g=R(h*z)^N is independent of variational stage integration.
                let (ga, dga) = period_map(a, period, steps);
                let (gd, dgd) = period_map(d, period, steps);
                let off_diagonal = coupling * (ga - gd) / (a - d);
                let off_diagonal_derivative = coupling * (dga - dgd) / (a - d);
                let expected = [
                    ga,
                    off_diagonal,
                    period_scale * (dga * initial[0] + off_diagonal_derivative * initial[1]),
                    0.0,
                    gd,
                    period_scale * dgd * initial[1],
                ];
                close(
                    image.end[0],
                    ga * initial[0] + off_diagonal * initial[1],
                    steps,
                );
                close(image.end[1], gd * initial[1], steps);
                for (actual, expected) in image.tangent.iter().zip(expected) {
                    close(*actual, expected, steps);
                }
                assert_eq!(image.tangent[3], 0.0);
            }
        }
    }
}
