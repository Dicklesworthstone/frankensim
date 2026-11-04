//! One RK4 stage owner for primal and consistent discrete variational flows.

use super::{Counts, Result, VariationalShootError, invalid, poll};
use crate::{OrbitProblem, TimeDomainRhs, TimeDomainStateJacobian};
use fs_la::gemm::gemm_scalar_checked;

#[cfg(test)]
mod tests;

pub(super) struct Image {
    pub end: Vec<f64>,
    /// Row-major d by (d + autonomous-period-column), in the scaled chart.
    pub tangent: Vec<f64>,
    pub excursion: f64,
}

/// All stages reuse this storage across period maps and line searches.
pub(super) struct Workspace {
    n: usize,
    columns: usize,
    state: Vec<f64>,
    stage: Vec<f64>,
    physical: Vec<f64>,
    slopes: [Vec<f64>; 4],
    tangent: Vec<f64>,
    stage_tangent: Vec<f64>,
    tangent_slopes: [Vec<f64>; 4],
    jacobian: Vec<f64>,
}

impl Workspace {
    pub fn new(n: usize, columns: usize) -> Self {
        Self {
            n,
            columns,
            state: vec![0.0; n],
            stage: vec![0.0; n],
            physical: vec![0.0; n],
            slopes: core::array::from_fn(|_| vec![0.0; n]),
            tangent: vec![0.0; n * columns],
            stage_tangent: vec![0.0; n * columns],
            tangent_slopes: core::array::from_fn(|_| vec![0.0; n * columns]),
            jacobian: vec![0.0; if columns == 0 { 0 } else { n * n }],
        }
    }

    #[allow(clippy::too_many_arguments)]
    pub fn integrate<P: OrbitProblem>(
        &mut self,
        problem: &P,
        rhs: TimeDomainRhs<P>,
        jacobian: Option<TimeDomainStateJacobian<P>>,
        initial: &[f64],
        period: f64,
        steps: usize,
        scales: &[f64],
        period_scale: f64,
        counts: &mut Counts,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Image> {
        poll(cancelled)?;
        counts.begin()?;
        let (n, columns) = (self.n, self.columns);
        let h = period / steps as f64;
        if !h.is_finite() || h <= 0.0 {
            return Err(invalid("positive representable RK4 step required"));
        }
        let with_tangent = jacobian.is_some();
        let period_column = with_tangent && columns > n;
        let dh = period_scale / steps as f64;
        self.state.copy_from_slice(initial);
        if with_tangent {
            self.tangent.fill(0.0);
            for i in 0..n {
                self.tangent[i * columns + i] = 1.0;
            }
        }
        let mut excursion = 0.0_f64;
        for step in 0..steps {
            poll(cancelled)?;
            let time = h * step as f64;
            for stage in 0..4 {
                let fraction = [0.0_f64, 0.5, 0.5, 1.0][stage];
                for i in 0..n {
                    self.stage[i] = if stage == 0 {
                        self.state[i]
                    } else if stage == 3 {
                        h.mul_add(self.slopes[stage - 1][i], self.state[i])
                    } else {
                        fraction.mul_add(h * self.slopes[stage - 1][i], self.state[i])
                    };
                    excursion = excursion.max((self.stage[i] - initial[i]).abs());
                    self.physical[i] = self.stage[i] * scales[i];
                }
                if self.physical.iter().any(|v| !v.is_finite()) {
                    return Err(invalid("physical RK4 stage is not representable"));
                }
                let at = time + fraction * h;
                self.slopes[stage].fill(f64::NAN);
                rhs(problem, at, &self.physical, &mut self.slopes[stage]);
                counts.rhs += 1;
                poll(cancelled)?;
                if self.slopes[stage].iter().any(|v| !v.is_finite()) {
                    return Err(VariationalShootError::SourceRefused {
                        time: at,
                        tangent: false,
                    });
                }
                for i in 0..n {
                    self.slopes[stage][i] /= scales[i];
                }
                if let Some(jacobian) = jacobian {
                    for i in 0..n {
                        poll(cancelled)?;
                        for j in 0..columns {
                            let index = i * columns + j;
                            self.stage_tangent[index] = if stage == 0 {
                                self.tangent[index]
                            } else {
                                fraction.mul_add(
                                    h * self.tangent_slopes[stage - 1][index],
                                    self.tangent[index],
                                )
                            };
                            // Differentiate the DISCRETE stage's step size too.
                            // The period variable is T / period_scale. Its time
                            // argument contributes zero only for autonomous f.
                            if period_column && j == n && stage != 0 {
                                self.stage_tangent[index] +=
                                    fraction * dh * self.slopes[stage - 1][i];
                            }
                        }
                    }
                    self.jacobian.fill(f64::NAN);
                    let admitted = jacobian(problem, at, &self.physical, &mut self.jacobian);
                    counts.jacobian += 1;
                    poll(cancelled)?;
                    if !admitted || self.jacobian.iter().any(|v| !v.is_finite()) {
                        return Err(VariationalShootError::SourceRefused {
                            time: at,
                            tangent: true,
                        });
                    }
                    for i in 0..n {
                        poll(cancelled)?;
                        for j in 0..n {
                            self.jacobian[i * n + j] =
                                (self.jacobian[i * n + j] * scales[j]) / scales[i];
                        }
                    }
                    if self
                        .jacobian
                        .iter()
                        .chain(&self.stage_tangent)
                        .any(|v| !v.is_finite())
                    {
                        return Err(invalid("scaled time-domain tangent is not representable"));
                    }
                    // Existing deterministic matrix-product owner; no separate
                    // variational matrix multiplication or solver is introduced.
                    // The allocation-free shared scalar kernel avoids packing
                    // buffers at every RK4 stage; row tiles bound cancellation.
                    for row in 0..n {
                        poll(cancelled)?;
                        gemm_scalar_checked(
                            1,
                            columns,
                            n,
                            1.0,
                            &self.jacobian[row * n..(row + 1) * n],
                            &self.stage_tangent,
                            0.0,
                            &mut self.tangent_slopes[stage][row * columns..(row + 1) * columns],
                        )
                        .map_err(|_| invalid("complete variational matrix product required"))?;
                    }
                    poll(cancelled)?;
                }
                if self.slopes[stage].iter().any(|v| !v.is_finite()) {
                    return Err(invalid("scaled time-domain flow is not representable"));
                }
            }
            for i in 0..n {
                let sum = self.slopes[3][i]
                    + self.slopes[0][i]
                    + 2.0 * (self.slopes[1][i] + self.slopes[2][i]);
                self.state[i] += h / 6.0 * sum;
                excursion = excursion.max((self.state[i] - initial[i]).abs());
                if with_tangent {
                    for j in 0..columns {
                        let index = i * columns + j;
                        self.tangent[index] += h / 6.0
                            * (self.tangent_slopes[3][index]
                                + self.tangent_slopes[0][index]
                                + 2.0
                                    * (self.tangent_slopes[1][index]
                                        + self.tangent_slopes[2][index]));
                        if period_column && j == n {
                            self.tangent[index] += dh / 6.0 * sum;
                        }
                    }
                }
            }
            if self.state.iter().any(|v| !v.is_finite())
                || (with_tangent && self.tangent.iter().any(|v| !v.is_finite()))
                || !excursion.is_finite()
            {
                return Err(invalid(
                    "discrete period flow or tangent is not representable",
                ));
            }
        }
        poll(cancelled)?;
        Ok(Image {
            end: self.state.clone(),
            tangent: if with_tangent {
                self.tangent.clone()
            } else {
                Vec::new()
            },
            excursion,
        })
    }
}
