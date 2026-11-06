//! Lagrangian averaging of the existing Favre--Germano stress identity.
//!
//! The two least-squares contractions travel with the resolved flow:
//! `D I_LM / Dt = (L:M - I_LM)/T`,
//! `D I_MM / Dt = (M:M - I_MM)/T`.
//! Backtracking uses the current physical velocity and trilinear interpolation
//! of the accepted history, followed by exact exponential relaxation for a
//! frozen current source. This is first-order trajectory/source integration;
//! its departure-distance limit is explicit. Numerators and denominators are
//! transported separately, and negative fitted coefficients produce zero
//! viscosity while their signed history is retained.
//!
//! This follows the pathline-averaging construction of Meneveau, Lund & Cabot
//! (1996), doi:10.1017/S0022112096007379. Here the positive memory time is
//! caller-declared: the paper's adaptive timescale and calibration are NOT
//! implemented. Initial history uses the existing local 27-cell average.
//! The grid/test filters remain one/three cells on a uniform periodic box.
//! No wall, inflow, backscatter-energy, scalar-flux or DNS claim is added.
//!
//! A prepared trial borrows its owner. Dropping it rejects the history;
//! committing it after the consuming flow step succeeds publishes history
//! without allocation or a fallible operation. Failed or rejected flow trials
//! therefore cannot age the turbulence memory.

use super::*;

/// Explicit controls; lengths, velocities, sample times and memory time must
/// use one consistent unit system. Native lattice adapters use lattice units.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LagrangianSgsConfig {
    /// Positive relaxation time T in the same units as the sample time.
    pub memory_time: f64,
    /// Maximum backtracked distance per axis in cells, in (0, 1].
    /// This bounds the first-order trajectory approximation; it is not a CFL
    /// stability theorem for the flow solver.
    pub maximum_departure_cells: f64,
    /// Upper admission bound on the fitted C_s squared.
    pub maximum_coefficient_squared: f64,
    /// Enforced by the consuming molecular transport owner, never a clamp.
    pub maximum_eddy_to_molecular: f64,
    pub usage: SgsUse,
}

impl LagrangianSgsConfig {
    pub fn validate(self) -> Result<(), DynamicSgsError> {
        DynamicSgsConfig {
            averaging: DynamicAveraging::LocalBox,
            maximum_coefficient_squared: self.maximum_coefficient_squared,
            maximum_eddy_to_molecular: self.maximum_eddy_to_molecular,
            usage: self.usage,
        }
        .validate()?;
        if !self.memory_time.is_finite() || self.memory_time <= 0.0 {
            return Err(DynamicSgsError::Invalid(
                "positive finite Lagrangian memory time required",
            ));
        }
        if !self.maximum_departure_cells.is_finite()
            || self.maximum_departure_cells <= 0.0
            || self.maximum_departure_cells > 1.0
        {
            return Err(DynamicSgsError::Invalid(
                "Lagrangian departure limit must be in (0, 1]",
            ));
        }
        Ok(())
    }
}

/// Reusable periodic filter plus two history banks (168 numerical bytes/cell).
/// A workspace belongs to one evolving trajectory; a new trajectory requires
/// a new owner. Preparation performs no allocation.
pub struct LagrangianDynamicSgs {
    flow: PeriodicDynamicSgs,
    history: Vec<[f64; 2]>,
    candidate: Vec<[f64; 2]>,
    accepted_time: Option<f64>,
    // Physical contractions scale as rho_scale^2 * velocity_scale^4.
    // A shared logarithmic scale avoids overflowing those fourth powers and
    // keeps histories comparable when current density/velocity scales change.
    history_log_scale: Option<f64>,
}

/// A fully admitted candidate. Its viscosities are valid until it is dropped
/// or committed. Commit only after the consuming flow transaction succeeds.
#[must_use = "commit after a successful flow step, or drop to reject the history"]
pub struct LagrangianSgsTrial<'a> {
    owner: &'a mut LagrangianDynamicSgs,
    time: f64,
    log_scale: Option<f64>,
    report: DynamicSgsReport,
}

impl LagrangianSgsTrial<'_> {
    pub fn viscosities(&self) -> &[f64] {
        self.owner.flow.viscosities()
    }

    pub const fn report(&self) -> DynamicSgsReport {
        self.report
    }

    /// Infallible, allocation-free publication of the successfully used fit.
    pub fn commit(self) {
        core::mem::swap(&mut self.owner.history, &mut self.owner.candidate);
        self.owner.accepted_time = Some(self.time);
        self.owner.history_log_scale = self.log_scale;
    }
}

impl LagrangianDynamicSgs {
    pub fn required_bytes(cells: usize) -> Result<usize, DynamicSgsError> {
        cells
            .checked_mul(21 * core::mem::size_of::<f64>())
            .ok_or(DynamicSgsError::Invalid(
                "Lagrangian workspace size overflow",
            ))
    }

    pub fn new(
        shape: [usize; 3],
        spacing: [f64; 3],
        maximum_bytes: usize,
    ) -> Result<Self, DynamicSgsError> {
        FilterGeometry::Periodic.validate_shape(shape)?;
        let cells = shape
            .into_iter()
            .try_fold(1usize, |n, d| n.checked_mul(d))
            .ok_or(DynamicSgsError::Invalid("Lagrangian cell count overflow"))?;
        let required = Self::required_bytes(cells)?;
        if required > maximum_bytes {
            return Err(DynamicSgsError::WorkspaceBytes {
                required,
                maximum: maximum_bytes,
            });
        }
        let result = Self {
            flow: PeriodicDynamicSgs::new_cartesian(shape, spacing, maximum_bytes)?,
            history: allocated(cells, [0.0; 2])?,
            candidate: allocated(cells, [0.0; 2])?,
            accepted_time: None,
            history_log_scale: None,
        };
        let required = result.workspace_bytes();
        if required > maximum_bytes {
            return Err(DynamicSgsError::WorkspaceBytes {
                required,
                maximum: maximum_bytes,
            });
        }
        Ok(result)
    }

    pub fn workspace_bytes(&self) -> usize {
        self.flow.workspace_bytes()
            + (self.history.capacity() + self.candidate.capacity())
                * core::mem::size_of::<[f64; 2]>()
    }

    pub const fn shape(&self) -> [usize; 3] {
        self.flow.shape()
    }

    /// Time of the last committed fit, not the consuming solver's endpoint.
    pub const fn accepted_time(&self) -> Option<f64> {
        self.accepted_time
    }

    /// Fit from current same-time physical fields. The first sample starts the
    /// history with a local-box fit. Later times must strictly increase; the
    /// elapsed interval is inferred from the last COMMITTED sample, so a
    /// rejected trial never shortens the retry's memory interval.
    pub fn prepare(
        &mut self,
        config: LagrangianSgsConfig,
        time: f64,
        density: &[f64],
        velocity: &[[f64; 3]],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<LagrangianSgsTrial<'_>, DynamicSgsError> {
        poll(&mut cancelled, 0)?;
        config.validate()?;
        if !time.is_finite()
            || time < 0.0
            || self.accepted_time.is_some_and(|previous| time <= previous)
        {
            return Err(DynamicSgsError::Invalid(
                "finite increasing Lagrangian sample times required",
            ));
        }
        let elapsed = self.accepted_time.map_or(0.0, |previous| time - previous);
        let (scales, _) =
            self.flow
                .prepare_contractions(density, velocity, [[0.0; 3]; 2], &mut cancelled)?;
        let current_log = (scales.velocity > 0.0)
            .then(|| 2.0 * det::ln(scales.density) + 4.0 * det::ln(scales.velocity));
        let common_log = match (self.history_log_scale, current_log) {
            (Some(a), Some(b)) => a.max(b),
            (Some(a), None) | (None, Some(a)) => a,
            (None, None) => 0.0,
        };
        let old_scale = self
            .history_log_scale
            .map_or(0.0, |value| det::exp(value - common_log));
        let new_scale = current_log.map_or(0.0, |value| det::exp(value - common_log));
        let relaxation = elapsed / config.memory_time;
        let fresh = -det::expm1(-relaxation);
        // Subtracting fresh from one erases representable, physically relevant
        // history once elapsed/T exceeds about 37, even while exp(-elapsed/T)
        // remains nonzero. This matters when strong eddies decay into weak flow.
        let old_weight = det::exp(-relaxation);
        if self.accepted_time.is_some() && fresh == 0.0 {
            return Err(DynamicSgsError::Invalid(
                "unrepresentable Lagrangian relaxation interval",
            ));
        }
        let mut report = DynamicSgsReport {
            minimum_raw_coefficient_squared: f64::INFINITY,
            ..DynamicSgsReport::default()
        };
        let mut amplitude: f64 = 0.0;
        for cell in 0..self.history.len() {
            poll(&mut cancelled, cell)?;
            let pair = if self.accepted_time.is_none() {
                let mut pair = [0.0; 2];
                let mut determined = false;
                each_box(self.flow.shape, cell, |j| {
                    determined |= self.flow.contraction[j][1] != 0.0;
                    for k in 0..2 {
                        pair[k] += self.flow.contraction[j][k] / 27.0;
                    }
                });
                if determined && pair[1] == 0.0 {
                    return Err(DynamicSgsError::Cell {
                        cell,
                        reason: "underflow of initial Lagrangian denominator",
                    });
                }
                pair
            } else {
                let departure = core::array::from_fn(|axis| {
                    velocity[cell][axis] * (elapsed / self.flow.spacing[axis])
                });
                if departure
                    .iter()
                    .any(|d: &f64| !d.is_finite() || d.abs() > config.maximum_departure_cells)
                {
                    return Err(DynamicSgsError::Cell {
                        cell,
                        reason: "Lagrangian departure exceeds declared cell limit",
                    });
                }
                let advected = sample_history(self.flow.shape, cell, departure, &self.history);
                let current = self.flow.contraction[cell];
                if (old_weight > 0.0 && old_scale == 0.0 && advected[1] != 0.0)
                    || (fresh > 0.0 && new_scale == 0.0 && current[1] != 0.0)
                {
                    return Err(DynamicSgsError::Cell {
                        cell,
                        reason: "unrepresentable Lagrangian history scale separation",
                    });
                }
                core::array::from_fn(|k| {
                    old_weight * (old_scale * advected[k]) + fresh * (new_scale * current[k])
                })
            };
            if pair.iter().any(|v| !v.is_finite())
                || pair[1] < 0.0
                || (pair[1] == 0.0 && pair[0] != 0.0)
            {
                return Err(DynamicSgsError::Cell {
                    cell,
                    reason: "unrepresentable Lagrangian contraction",
                });
            }
            let raw = if pair[1] == 0.0 {
                report.undetermined_cells += 1;
                0.0
            } else {
                pair[0] / pair[1]
            };
            if !raw.is_finite() {
                return Err(DynamicSgsError::Cell {
                    cell,
                    reason: "unrepresentable Lagrangian coefficient",
                });
            }
            report.minimum_raw_coefficient_squared =
                report.minimum_raw_coefficient_squared.min(raw);
            report.backscatter_limited_cells += usize::from(raw < 0.0);
            let coefficient = raw.max(0.0);
            if coefficient > config.maximum_coefficient_squared {
                return Err(DynamicSgsError::Coefficient {
                    cell,
                    observed: coefficient,
                    maximum: config.maximum_coefficient_squared,
                });
            }
            let viscosity = coefficient
                * self.flow.strain_rate[cell]
                * scales.velocity
                * self.flow.filter_width;
            if !viscosity.is_finite() {
                return Err(DynamicSgsError::Cell {
                    cell,
                    reason: "unrepresentable Lagrangian viscosity",
                });
            }
            self.flow.viscosity[cell] = viscosity;
            self.flow.basis[cell][0] = coefficient;
            self.candidate[cell] = pair;
            amplitude = amplitude.max(pair[0].abs()).max(pair[1]);
            report.maximum_coefficient_squared =
                report.maximum_coefficient_squared.max(coefficient);
            report.maximum_eddy_viscosity = report.maximum_eddy_viscosity.max(viscosity);
        }
        // Renormalize after the fit. Otherwise a formerly strong flow would
        // eventually underflow its decaying memory while weaker new eddies form.
        let log_scale = if amplitude == 0.0 {
            None
        } else {
            for (cell, pair) in self.candidate.iter_mut().enumerate() {
                poll(&mut cancelled, cell)?;
                for value in pair {
                    let normalized = *value / amplitude;
                    if normalized == 0.0 && *value != 0.0 {
                        return Err(DynamicSgsError::Cell {
                            cell,
                            reason: "Lagrangian history normalization underflow",
                        });
                    }
                    *value = normalized;
                }
            }
            Some(common_log + det::ln(amplitude))
        };
        poll(&mut cancelled, 0)?;
        Ok(LagrangianSgsTrial {
            owner: self,
            time,
            log_scale,
            report,
        })
    }
}

// Positive trilinear weights preserve the sign of the denominator. The
// physical velocity, including uniform advection, defines the departure;
// subtracting the covariance's reference velocity here would pin eddies.
fn sample_history(
    shape: [usize; 3],
    cell: usize,
    departure: [f64; 3],
    history: &[[f64; 2]],
) -> [f64; 2] {
    let p = coordinates(shape, cell);
    let coordinate = core::array::from_fn::<_, 3, _>(|a| {
        (p[a] as f64 - departure[a]).rem_euclid(shape[a] as f64)
    });
    let lower = coordinate.map(|v| v.floor() as usize);
    let fraction = core::array::from_fn::<_, 3, _>(|a| coordinate[a] - lower[a] as f64);
    let mut result = [0.0; 2];
    for z in 0..2 {
        for y in 0..2 {
            for x in 0..2 {
                let offset = [x, y, z];
                let q = core::array::from_fn(|a| (lower[a] + offset[a]) % shape[a]);
                let weight = (0..3)
                    .map(|a| {
                        if offset[a] == 0 {
                            1.0 - fraction[a]
                        } else {
                            fraction[a]
                        }
                    })
                    .product::<f64>();
                for k in 0..2 {
                    result[k] += weight * history[index(shape, q)][k];
                }
            }
        }
    }
    result
}

#[cfg(test)]
#[path = "lagrangian/tests.rs"]
mod tests;
