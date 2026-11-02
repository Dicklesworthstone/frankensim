//! Stamp original structural and load equations; the shared pH owner propagates.
use super::*;
use fs_math::det;
use fs_phs::linear_exact::prepared::PreparedQuadraticFlow;
use fs_phs::linear_exact::{MAX_LINEAR_EXPONENTIAL_BLOCK_DIM, split_quadratic};
use fs_phs::{PortHamiltonian, PreparedStepError};

pub(super) fn invalid(message: &str) -> AssemblyMimoTransientError {
    AssemblyMimoTransientError::Structure(message.into())
}
pub(super) fn check(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(AssemblyMimoTransientError::Cancelled)
    } else {
        Ok(())
    }
}

pub(super) fn preflight(
    n: usize,
    p: usize,
    source: &PreparedMimoLoad,
    physical: u128,
    components: usize,
    budget: AssemblyMimoTransientBudget,
) -> Result<()> {
    let states = 2 * n as u128 + source.continuous_operator().dimension() as u128;
    let d = states + n as u128 + 1;
    if n == 0
        || p == 0
        || p != source.ports()
        || p > budget.maximum_ports
        || states > budget.maximum_real_states as u128
        || d > MAX_LINEAR_EXPONENTIAL_BLOCK_DIM as u128
        || d > budget.exact.max_dimension as u128
    {
        return Err(invalid(
            "complete joint source or physical port capacity refused",
        ));
    }
    let p128 = p as u128;
    let bytes = 8
        * (5 * states * states + 2 * p128 * states + 24 * states + 24 * p128 + 8 * physical)
        + 240 * components as u128
        + 2 * source.limits().max_memory_bytes as u128;
    let terms = 16 * states * states * states
        + 16 * p128 * states * states
        + 16 * p128 * p128 * states
        + 16 * physical * (p128 + n as u128 + 1);
    if bytes > budget.maximum_connection_bytes as u128
        || terms > budget.maximum_connection_terms as u128
        || !budget.maximum_relative_balance_error.is_finite()
        || budget.maximum_relative_balance_error <= 0.
        || budget.maximum_relative_balance_error > 1e-8
    {
        return Err(invalid(
            "complete joint source, physical ports or connection budget refused",
        ));
    }
    Ok(())
}

pub(super) struct PreparedJointFlow {
    pub(super) flow: PreparedQuadraticFlow,
    flow_rows: Vec<f64>,
    effort_rows: Vec<f64>,
    force_map: Vec<f64>,
    n: usize,
    p: usize,
}
impl PreparedJointFlow {
    #[allow(clippy::too_many_arguments)]
    pub(super) fn new(
        system: &PortHamiltonian,
        stiffness: &[f64],
        inverse_mass: &[f64],
        damping: &[f64],
        force_map: &[f64],
        source: &PreparedMimoLoad,
        budget: AssemblyMimoTransientBudget,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        check(cancelled)?;
        let n = system.port_dim();
        let p = source.ports();
        preflight(n, p, source, 0, 0, budget)?;
        if system.state_dim() != 2 * n
            || stiffness.len() != n * n
            || inverse_mass.len() != n * n
            || damping.len() != n * n
            || force_map.len() != n * p
            || stiffness
                .iter()
                .chain(inverse_mass)
                .chain(damping)
                .chain(force_map)
                .any(|v| !v.is_finite())
        {
            return Err(invalid(
                "original canonical source and physical load map dimensions disagree",
            ));
        }
        let continuous = source.continuous_operator();
        let t = continuous.dimension();
        let states = 2 * n + t;
        let mut interrupted = false;
        let original = continuous.state_space(budget.maximum_connection_bytes, &mut || {
            let stop = cancelled();
            interrupted |= stop;
            stop
        });
        if interrupted {
            return Err(AssemblyMimoTransientError::Cancelled);
        }
        let original = original.map_err(AssemblyMimoTransientError::Source)?;
        check(cancelled)?;
        let mut velocities = vec![vec![0.; states]; n];
        for i in 0..n {
            velocities[i][n..2 * n].copy_from_slice(&inverse_mass[i * n..(i + 1) * n]);
        }
        let mut flow_rows = vec![0.; p * states];
        for a in 0..p {
            check(cancelled)?;
            for j in 0..n {
                flow_rows[a * states + n + j] = (0..n)
                    .map(|i| force_map[i * p + a] * inverse_mass[i * n + j])
                    .sum();
            }
        }
        let mut effort_rows = vec![0.; p * states];
        for a in 0..p {
            check(cancelled)?;
            for j in 0..states {
                effort_rows[a * states + j] = (0..p)
                    .map(|b| original.d[a * p + b] * flow_rows[b * states + j])
                    .sum();
            }
            effort_rows[a * states + 2 * n..(a + 1) * states]
                .copy_from_slice(&original.c[a * t..(a + 1) * t]);
        }
        let mut generator = vec![0.; states * states];
        let mut damped_velocity = vec![0.; n * states];
        for i in 0..n {
            check(cancelled)?;
            generator[i * states..(i + 1) * states].copy_from_slice(&velocities[i]);
            for j in 0..states {
                damped_velocity[i * states + j] =
                    (0..n).map(|k| damping[i * n + k] * velocities[k][j]).sum();
                generator[(n + i) * states + j] = -damped_velocity[i * states + j]
                    - (0..p)
                        .map(|a| force_map[i * p + a] * effort_rows[a * states + j])
                        .sum::<f64>();
            }
            for j in 0..n {
                generator[(n + i) * states + j] -= stiffness[i * n + j];
            }
        }
        for i in 0..t {
            check(cancelled)?;
            for j in 0..states {
                generator[(2 * n + i) * states + j] = (0..p)
                    .map(|a| original.b[i * p + a] * flow_rows[a * states + j])
                    .sum();
            }
            generator[(2 * n + i) * states + 2 * n..(2 * n + i + 1) * states]
                .copy_from_slice(&original.a[i * t..(i + 1) * t]);
        }
        // Symmetric observation coefficients preserve all signed/mutual work.
        // This symmetry changes no x^T Q x and never modifies the generator.
        let mut glue = vec![0.; states * states];
        let mut work = vec![0.; states * states];
        for i in 0..states {
            check(cancelled)?;
            for j in 0..=i {
                let c = 0.5
                    * (0..n)
                        .map(|k| {
                            velocities[k][i] * damped_velocity[k * states + j]
                                + velocities[k][j] * damped_velocity[k * states + i]
                        })
                        .sum::<f64>();
                let w = 0.5
                    * (0..p)
                        .map(|a| {
                            flow_rows[a * states + i] * effort_rows[a * states + j]
                                + flow_rows[a * states + j] * effort_rows[a * states + i]
                        })
                        .sum::<f64>();
                glue[i * states + j] = c;
                glue[j * states + i] = c;
                work[i * states + j] = w;
                work[j * states + i] = w;
            }
        }
        let mut extra_scales = Vec::with_capacity(t);
        for j in 0..t {
            let magnitude = (0..p)
                .map(|a| original.c[a * t + j].abs())
                .fold(0., f64::max);
            // Numerical conditioning only; no load storage or heat is inferred.
            extra_scales.push(if magnitude > 0. {
                det::sqrt(magnitude)
            } else {
                1.
            });
        }
        let mut linear = velocities;
        linear.extend(effort_rows.chunks_exact(states).map(<[f64]>::to_vec));
        let split = split_quadratic(system, budget.structure.linear)
            .map_err(|e| AssemblyMimoTransientError::Flow(e.to_string()))?;
        check(cancelled)?;
        let flow = PreparedQuadraticFlow::new_extended_controlled(
            system,
            &split,
            &generator,
            &extra_scales,
            source.time_step_s(),
            &[glue, work],
            &linear,
            budget.exact,
            cancelled,
        )
        .map_err(|e| match e {
            PreparedStepError::Cancelled => AssemblyMimoTransientError::Cancelled,
            PreparedStepError::Solver(e) => AssemblyMimoTransientError::Flow(e.to_string()),
        })?;
        Ok(Self {
            flow,
            flow_rows,
            effort_rows,
            force_map: force_map.to_vec(),
            n,
            p,
        })
    }

    pub(super) fn endpoint_flows(&self, state: &[f64]) -> Result<Vec<f64>> {
        multiply(&self.flow_rows, self.p, state)
    }
    pub(super) fn endpoint_efforts(&self, state: &[f64]) -> Result<Vec<f64>> {
        multiply(&self.effort_rows, self.p, state)
    }
    pub(super) fn project_velocities(&self, velocities: &[f64]) -> Result<Vec<f64>> {
        if velocities.len() != self.n {
            return Err(invalid("complete velocity observation required"));
        }
        let mut result = vec![0.; self.p];
        for (a, value) in result.iter_mut().enumerate() {
            *value = (0..self.n)
                .map(|i| self.force_map[i * self.p + a] * velocities[i])
                .sum();
        }
        if result.iter().any(|v| !v.is_finite()) {
            return Err(invalid("port flow integral overflow"));
        }
        Ok(result)
    }
}

pub(super) fn dot(a: &[f64], b: &[f64]) -> Result<f64> {
    if a.len() != b.len() {
        return Err(invalid("incomplete physical contraction"));
    }
    let sum: f64 = a.iter().zip(b).map(|(a, b)| a * b).sum();
    if !sum.is_finite() {
        return Err(invalid("physical contraction overflow"));
    }
    Ok(sum)
}
pub(super) fn multiply(matrix: &[f64], rows: usize, vector: &[f64]) -> Result<Vec<f64>> {
    if vector.is_empty() || matrix.len() != rows * vector.len() {
        return Err(invalid("complete original matrix/state shape required"));
    }
    matrix
        .chunks_exact(vector.len())
        .map(|row| dot(row, vector))
        .collect()
}

pub(super) fn energy_change(k: &[f64], inverse_mass: &[f64], a: &[f64], b: &[f64]) -> Result<f64> {
    let n = a.len() / 2;
    let mut value = 0.;
    for (matrix, start) in [(k, 0), (inverse_mass, n)] {
        for i in 0..n {
            let effort: f64 = (0..n)
                .map(|j| matrix[i * n + j] * (0.5 * a[start + j] + 0.5 * b[start + j]))
                .sum();
            value += (b[start + i] - a[start + i]) * effort;
        }
    }
    if !value.is_finite() {
        return Err(invalid("structural energy change overflow"));
    }
    Ok(value)
}
