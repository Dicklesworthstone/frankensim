//! Symmetric composition of exact anisotropic liquid-cell/body source maps.
//!
//! K maps mass-whitened body velocity to one cell's rigid-point velocity.
//! With cell mass m, relative mobility is G = I + m K K^T. A principal
//! rate factor Z gives Gamma = rate_scale Z Z^T. The symmetric metric
//! Z^T G Z, rather than the generally nonsymmetric G Gamma, supplies the
//! exact local exponential. Zero principal rates are never inverted.
//!
//! Only source increments are accumulated. No native state, clock, population
//! or thermal owner is touched by this compact, cancellable numerical kernel.

use fs_math::dd::Dd;

#[derive(Clone, Debug)]
pub(super) struct Cell {
    pub mass: f64,
    /// Initial liquid velocity minus the initial rigid-point velocity [m/s].
    pub slip: [f64; 3],
    /// Columns of the original body's kinetic-metric point-velocity map.
    pub body_columns: [[f64; 3]; 6],
    pub principal_rates_per_s: [f64; 3],
    /// Original proper principal axes resolved in the common world frame.
    pub principal_axes: [[f64; 3]; 3],
}

#[derive(Clone, Debug, PartialEq)]
pub(super) struct Plan {
    pub impulses: Vec<[f64; 3]>,
    pub body_increment: [f64; 6],
    pub dissipation_j: f64,
    pub pair_exchanges: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub(super) enum Error {
    Invalid(&'static str),
    Cancelled,
}

struct Response {
    // These columns have unit G metric, avoiding division by small rates.
    modes: [[f64; 3]; 3],
    half_impulse: [f64; 3],
    full_impulse: [f64; 3],
    half_loss: [f64; 3],
    full_loss: [f64; 3],
    active: bool,
}

fn dd(x: f64) -> Dd {
    Dd::from_f64(x)
}
fn dot(a: [f64; 3], b: [f64; 3]) -> Dd {
    (0..3).fold(Dd::ZERO, |sum, i| sum + dd(a[i]) * dd(b[i]))
}
fn invalid(reason: &'static str) -> Error {
    Error::Invalid(reason)
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}

/// Worst-case count, before allocating any selected-cell state.
pub(super) fn work_bound(cells: usize, substeps: usize, maximum: usize) -> Result<usize, Error> {
    if cells == 0 || substeps == 0 {
        return Err(invalid(
            "heterogeneous drag needs nonempty cells and positive substeps",
        ));
    }
    let work = cells
        .checked_mul(2)
        .and_then(|n| n.checked_sub(1))
        .and_then(|n| n.checked_mul(substeps))
        .ok_or_else(|| invalid("heterogeneous drag pair count overflows"))?;
    if work > maximum {
        return Err(invalid(
            "heterogeneous drag exceeds the declared pair-exchange budget",
        ));
    }
    Ok(work)
}

/// Prefer a normal intermediate product: a huge rate and duration can still
/// have a small exposure in a weak principal direction. Positive overflow is
/// the exponential-zero limit only after the third factor has been included.
fn exposure(a: f64, b: f64, c: f64) -> f64 {
    if a == 0.0 || b == 0.0 || c == 0.0 {
        return 0.0;
    }
    for (x, y, z) in [(a, b, c), (a, c, b), (b, c, a)] {
        let product = x * y;
        if product.is_normal() {
            return product * z;
        }
    }
    a * b * c
}

impl Response {
    fn new(cell: &Cell, half_step: f64) -> Result<Self, Error> {
        if !cell.mass.is_finite()
            || cell.mass <= 0.0
            || cell.slip.iter().any(|v| !v.is_finite())
            || cell.body_columns.iter().flatten().any(|v| !v.is_finite())
            || cell.principal_axes.iter().flatten().any(|v| !v.is_finite())
            || cell
                .principal_rates_per_s
                .iter()
                .any(|v| !v.is_finite() || *v < 0.0)
        {
            return Err(invalid(
                "heterogeneous drag needs finite inertia, slip and nonnegative rates",
            ));
        }
        for i in 0..3 {
            for j in i..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                if (dot(cell.principal_axes[i], cell.principal_axes[j]).to_f64() - expected).abs()
                    > 256.0 * f64::EPSILON
                {
                    return Err(invalid(
                        "drag principal axes must be an original orthonormal frame",
                    ));
                }
            }
        }
        let rate = cell
            .principal_rates_per_s
            .into_iter()
            .fold(0.0_f64, f64::max);
        let mut response = Self {
            modes: [[0.0; 3]; 3],
            half_impulse: [0.0; 3],
            full_impulse: [0.0; 3],
            half_loss: [0.0; 3],
            full_loss: [0.0; 3],
            active: rate > 0.0 && half_step > 0.0,
        };
        if !response.active {
            return Ok(response);
        }
        let mut root = [[0.0; 3]; 3];
        for i in 0..3 {
            let ratio = cell.principal_rates_per_s[i] / rate;
            if cell.principal_rates_per_s[i] > 0.0 && ratio == 0.0 {
                return Err(invalid("principal drag rate ratio is unrepresentable"));
            }
            root[i] = cell.principal_axes[i].map(|v| v * ratio.sqrt());
        }
        let mapped: [[f64; 6]; 3] = core::array::from_fn(|i| {
            core::array::from_fn(|j| dot(root[i], cell.body_columns[j]).to_f64())
        });
        let mut matrix = [0.0; 9];
        for i in 0..3 {
            for j in i..3 {
                let value = (0..6)
                    .fold(dot(root[i], root[j]), |sum, k| {
                        sum + dd(cell.mass) * dd(mapped[i][k]) * dd(mapped[j][k])
                    })
                    .to_f64();
                if !value.is_finite() {
                    return Err(invalid("heterogeneous drag mobility is unrepresentable"));
                }
                matrix[3 * i + j] = value;
                matrix[3 * j + i] = value;
            }
        }
        let scale = matrix.into_iter().fold(0.0_f64, |a, x| a.max(x.abs()));
        if !scale.is_finite() || scale <= 0.0 {
            return Err(invalid(
                "active directional drag has no representable kinetic metric",
            ));
        }
        let normalized = matrix.map(|v| v / scale);
        let (values, vectors) = fs_la::eigen::jacobi_eigh(&normalized, 3);
        let rank = cell
            .principal_rates_per_s
            .iter()
            .filter(|&&v| v > 0.0)
            .count();
        for k in 0..3 {
            let value = values[k];
            if !value.is_finite() || value < -256.0 * f64::EPSILON {
                return Err(invalid(
                    "heterogeneous drag positive metric failed eigensolve",
                ));
            }
            // Exactly null source directions have zero rows/columns in the
            // principal basis. Never turn an unresolved positive rate into a
            // free-slip direction or divide by a zero eigenvalue.
            if k < 3 - rank {
                if value != 0.0 {
                    return Err(invalid("heterogeneous drag null mode was not retained"));
                }
                continue;
            }
            let lambda = value * scale;
            if !lambda.is_finite() || lambda <= 0.0 {
                return Err(invalid("positive directional drag mode is unresolved"));
            }
            let mode = core::array::from_fn(|axis| {
                ((0..3).fold(Dd::ZERO, |sum, j| {
                    sum + dd(root[j][axis]) * dd(vectors[3 * j + k])
                }) / dd(lambda.sqrt()))
                .to_f64()
            });
            if mode.iter().any(|v| !v.is_finite()) {
                return Err(invalid("directional drag response is unrepresentable"));
            }
            let norm = (0..6)
                .fold(dot(mode, mode), |sum, j| {
                    let projection = dot(mode, cell.body_columns[j]);
                    sum + dd(cell.mass) * projection * projection
                })
                .to_f64();
            if !norm.is_finite() || (norm - 1.0).abs() > 4096.0 * f64::EPSILON {
                return Err(invalid(
                    "directional drag eigenmode has unresolved kinetic norm",
                ));
            }
            response.modes[k] = mode;
            let theta = exposure(rate, half_step, lambda);
            if theta == 0.0 {
                return Err(invalid(
                    "positive directional drag exposure is unrepresentable",
                ));
            }
            let half = -(-theta).exp_m1();
            let full = -(-2.0 * theta).exp_m1();
            let full_loss = -(-4.0 * theta).exp_m1();
            if [half, full, full_loss]
                .iter()
                .any(|v| !v.is_finite() || *v < 0.0 || *v > 1.0)
            {
                return Err(invalid("directional drag exponential is unrepresentable"));
            }
            response.half_impulse[k] = half;
            response.full_impulse[k] = full;
            response.half_loss[k] = full;
            response.full_loss[k] = full_loss;
        }
        // Independent metric orthogonality admission catches an unresolved
        // eigensystem before its local physical dissipation is trusted.
        for i in 0..3 {
            for j in i + 1..3 {
                let a = response.modes[i];
                let b = response.modes[j];
                let overlap = (0..6)
                    .fold(dot(a, b), |sum, k| {
                        sum + dd(cell.mass)
                            * dot(a, cell.body_columns[k])
                            * dot(b, cell.body_columns[k])
                    })
                    .to_f64();
                if !overlap.is_finite() || overlap.abs() > 4096.0 * f64::EPSILON {
                    return Err(invalid(
                        "directional drag eigenmodes are not kinetically orthogonal",
                    ));
                }
            }
        }
        Ok(response)
    }
}

struct Accumulated {
    impulses: Vec<[Dd; 3]>,
    body: [Dd; 6],
    loss: Dd,
    pair_exchanges: usize,
}

impl Accumulated {
    fn exchange(
        &mut self,
        index: usize,
        cell: &Cell,
        response: &Response,
        full: bool,
    ) -> Result<(), Error> {
        if !response.active {
            return Ok(());
        }
        let slip: [f64; 3] = core::array::from_fn(|axis| {
            let body_change = (0..6).fold(Dd::ZERO, |sum, k| {
                sum + dd(cell.body_columns[k][axis]) * self.body[k]
            });
            (dd(cell.slip[axis]) + self.impulses[index][axis] / dd(cell.mass) - body_change)
                .to_f64()
        });
        if slip.iter().any(|v| !v.is_finite()) {
            return Err(invalid("composed liquid/body slip is unrepresentable"));
        }
        let mut impulse = [Dd::ZERO; 3];
        let mut loss = Dd::ZERO;
        for mode in 0..3 {
            let amplitude = dot(response.modes[mode], slip);
            let weight = if full {
                response.full_impulse[mode]
            } else {
                response.half_impulse[mode]
            };
            let loss_weight = if full {
                response.full_loss[mode]
            } else {
                response.half_loss[mode]
            };
            for axis in 0..3 {
                impulse[axis] = impulse[axis]
                    - dd(cell.mass) * dd(response.modes[mode][axis]) * dd(weight) * amplitude;
            }
            loss = loss + dd(0.5) * dd(cell.mass) * dd(loss_weight) * amplitude * amplitude;
        }
        if impulse
            .iter()
            .any(|j| !j.hi.is_finite() || !j.lo.is_finite())
            || !loss.hi.is_finite()
            || !loss.lo.is_finite()
            || loss.hi < 0.0
        {
            return Err(invalid(
                "exact local drag impulse or dissipation is unrepresentable",
            ));
        }
        for axis in 0..3 {
            self.impulses[index][axis] = self.impulses[index][axis] + impulse[axis];
        }
        for k in 0..6 {
            let recoil = (0..3).fold(Dd::ZERO, |sum, axis| {
                sum + dd(cell.body_columns[k][axis]) * impulse[axis]
            });
            self.body[k] = self.body[k] - recoil;
        }
        self.loss = self.loss + loss;
        self.pair_exchanges += 1;
        Ok(())
    }
}

pub(super) fn solve(
    cells: &[Cell],
    duration_s: f64,
    substeps: usize,
    maximum_pair_exchanges: usize,
    mut cancelled: impl FnMut() -> bool,
) -> Result<Plan, Error> {
    poll(&mut cancelled)?;
    work_bound(cells.len(), substeps, maximum_pair_exchanges)?;
    if !duration_s.is_finite() || duration_s < 0.0 {
        return Err(invalid(
            "heterogeneous drag duration must be finite and nonnegative",
        ));
    }
    let half_step = (duration_s / substeps as f64) * 0.5;
    if duration_s > 0.0 && half_step == 0.0 {
        return Err(invalid("heterogeneous drag subdivision is unrepresentable"));
    }
    let mut responses = Vec::new();
    responses
        .try_reserve_exact(cells.len())
        .map_err(|_| invalid("drag response allocation refused"))?;
    let mut impulses = Vec::new();
    impulses
        .try_reserve_exact(cells.len())
        .map_err(|_| invalid("drag increment allocation refused"))?;
    for cell in cells {
        poll(&mut cancelled)?;
        responses.push(Response::new(cell, half_step)?);
        impulses.push([Dd::ZERO; 3]);
    }
    let mut state = Accumulated {
        impulses,
        body: [Dd::ZERO; 6],
        loss: Dd::ZERO,
        pair_exchanges: 0,
    };
    // Merge the adjacent half exchanges of the last cell into one exact full
    // exchange. The remaining order is palindromic, so the frozen source is
    // second order under subdivision, without a false global-exactness claim.
    if duration_s != 0.0 {
        for _ in 0..substeps {
            for index in 0..cells.len() {
                poll(&mut cancelled)?;
                state.exchange(
                    index,
                    &cells[index],
                    &responses[index],
                    index + 1 == cells.len(),
                )?;
            }
            for index in (0..cells.len() - 1).rev() {
                poll(&mut cancelled)?;
                state.exchange(index, &cells[index], &responses[index], false)?;
            }
        }
    }
    poll(&mut cancelled)?;
    let body_increment = state.body.map(Dd::to_f64);
    let dissipation_j = state.loss.to_f64();
    let mut impulses = Vec::new();
    impulses
        .try_reserve_exact(cells.len())
        .map_err(|_| invalid("drag result allocation refused"))?;
    for impulse in state.impulses {
        impulses.push(impulse.map(Dd::to_f64));
    }
    if body_increment
        .iter()
        .chain(impulses.iter().flatten())
        .chain([&dissipation_j])
        .any(|v| !v.is_finite())
        || dissipation_j < 0.0
    {
        return Err(invalid("composed drag source is unrepresentable"));
    }
    Ok(Plan {
        impulses,
        body_increment,
        dissipation_j,
        pair_exchanges: state.pair_exchanges,
    })
}

#[cfg(test)]
mod tests;
