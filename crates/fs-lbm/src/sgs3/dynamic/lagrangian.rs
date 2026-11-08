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
//! (1996), doi:10.1017/S0022112096007379. The default positive memory time is
//! caller-declared. [`AdaptiveLagrangianMemory`] optionally derives a bounded
//! local timescale from the advected accepted stress contractions, with an
//! explicit density normalization for the Favre extension. It supplies no
//! empirical calibration. Periodic history starts from the local 27-cell average.
//! A resolved plane channel instead uses the existing tangential test filter
//! and a local nine-cell initial average, not homogeneous-plane averaging.
//! Its actual no-slip wall velocities enter the original Favre contractions.
//!
//! Channel backtracking wraps only the two tangential axes. Normal departures
//! must remain between the physical walls, half a cell beyond the outer
//! centers. Positive interpolation extends the nearest center through each
//! outer half cell; this declared first-order boundary extrapolation never
//! joins opposite walls or invents inflow history. A departure beyond a wall
//! refuses the trial. No wall law, inflow, backscatter-energy or DNS claim is
//! added. Independent scalar-flux histories use the opt-in [`scalar`] owner.
//!
//! A prepared trial borrows its owner. Dropping it rejects the history;
//! committing it after the consuming flow step succeeds publishes history
//! without allocation or a fallible operation. Failed or rejected flow trials
//! therefore cannot age the turbulence memory.

use super::*;

#[path = "lagrangian/memory.rs"]
mod memory;
pub use memory::AdaptiveLagrangianMemory;
use memory::AdaptiveUpdate;

/// Joint momentum/scalar pathline histories with independent scalar transport.
#[path = "lagrangian/scalar.rs"]
pub mod scalar;

/// Explicit controls; lengths, velocities, sample times and memory time must
/// use one consistent unit system. Native lattice adapters use lattice units.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LagrangianSgsConfig {
    /// Positive relaxation time T in the same units as the sample time.
    /// With an explicitly adaptive owner this is the upper bound on local T.
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
    adaptive_memory: Option<AdaptiveLagrangianMemory>,
    accepted_memory_time_range: Option<[f64; 2]>,
}

/// A fully admitted candidate. Its viscosities are valid until it is dropped
/// or committed. Commit only after the consuming flow transaction succeeds.
#[must_use = "commit after a successful flow step, or drop to reject the history"]
pub struct LagrangianSgsTrial<'a> {
    owner: &'a mut LagrangianDynamicSgs,
    time: f64,
    log_scale: Option<f64>,
    report: DynamicSgsReport,
    memory_time_range: [f64; 2],
}

impl LagrangianSgsTrial<'_> {
    pub fn viscosities(&self) -> &[f64] {
        self.owner.flow.viscosities()
    }

    /// Nonnegative fitted C_s squared, independently of centered cell strain.
    /// Conservative face consumers must combine it with their own face strain.
    /// This iterator is available only on a fully admitted history trial.
    pub fn coefficients_squared(&self) -> impl ExactSizeIterator<Item = f64> + '_ {
        self.owner.flow.coefficients_squared()
    }

    /// Fitted coefficient for one cell; out-of-domain indices return None.
    pub fn coefficient_squared(&self, cell: usize) -> Option<f64> {
        self.owner.flow.coefficient_squared(cell)
    }

    pub const fn report(&self) -> DynamicSgsReport {
        self.report
    }

    /// Minimum/maximum relaxation time used by this admitted candidate.
    /// The initial local-box sample reports the authored upper time bound.
    pub const fn memory_time_range(&self) -> [f64; 2] {
        self.memory_time_range
    }

    /// Infallible, allocation-free publication of the successfully used fit.
    pub fn commit(self) {
        core::mem::swap(&mut self.owner.history, &mut self.owner.candidate);
        self.owner.accepted_time = Some(self.time);
        self.owner.history_log_scale = self.log_scale;
        self.owner.accepted_memory_time_range = Some(self.memory_time_range);
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
        Self::new_geometry(shape, spacing, FilterGeometry::Periodic, maximum_bytes)
    }

    /// Resolved no-slip channel with two periodic axes and one fixed wall normal.
    /// Wall velocities are supplied to each `prepare_with_walls` call. Storage
    /// and allocation bounds are identical to the periodic owner.
    pub fn new_channel(
        shape: [usize; 3],
        spacing: [f64; 3],
        normal_axis: usize,
        maximum_bytes: usize,
    ) -> Result<Self, DynamicSgsError> {
        Self::new_geometry(
            shape,
            spacing,
            FilterGeometry::Channel(normal_axis),
            maximum_bytes,
        )
    }

    fn new_geometry(
        shape: [usize; 3],
        spacing: [f64; 3],
        geometry: FilterGeometry,
        maximum_bytes: usize,
    ) -> Result<Self, DynamicSgsError> {
        geometry.validate_shape(shape)?;
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
            flow: PeriodicDynamicSgs::new_cartesian_geometry(
                shape,
                spacing,
                geometry,
                maximum_bytes,
            )?,
            history: allocated(cells, [0.0; 2])?,
            candidate: allocated(cells, [0.0; 2])?,
            accepted_time: None,
            history_log_scale: None,
            adaptive_memory: None,
            accepted_memory_time_range: None,
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

    /// Fixed impermeable wall normal, or None for the periodic owner.
    pub const fn normal_axis(&self) -> Option<usize> {
        match self.flow.geometry {
            FilterGeometry::Periodic => None,
            FilterGeometry::Channel(axis) => Some(axis),
        }
    }

    /// Physical grid-filter width used by the Germano model.
    pub const fn filter_width(&self) -> f64 {
        self.flow.filter_width()
    }

    /// Time of the last committed fit, not the consuming solver's endpoint.
    pub const fn accepted_time(&self) -> Option<f64> {
        self.accepted_time
    }

    /// Select a bounded invariant timescale before this trajectory starts.
    /// Momentum history remains signed; heat/component memories are separate.
    /// No field arrays or additional numerical workspace are allocated.
    pub fn with_adaptive_memory(
        mut self,
        config: AdaptiveLagrangianMemory,
    ) -> Result<Self, DynamicSgsError> {
        config.validate()?;
        if self.accepted_time.is_some() {
            return Err(DynamicSgsError::Invalid(
                "select adaptive Lagrangian memory before the first accepted fit",
            ));
        }
        self.adaptive_memory = Some(config);
        Ok(self)
    }

    pub const fn adaptive_memory(&self) -> Option<AdaptiveLagrangianMemory> {
        self.adaptive_memory
    }

    /// Timescale range of the last committed fit, never a rejected trial.
    pub const fn accepted_memory_time_range(&self) -> Option<[f64; 2]> {
        self.accepted_memory_time_range
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
        cancelled: impl FnMut() -> bool,
    ) -> Result<LagrangianSgsTrial<'_>, DynamicSgsError> {
        self.prepare_scaled(config, time, density, velocity, cancelled)
            .map(|(trial, _)| trial)
    }

    /// Prepare channel history from current physical fields and both no-slip
    /// wall velocities (lower, upper). Normal wall speeds must be zero.
    /// This entry requires a channel owner; periodic preparation requires `prepare`.
    /// The geometry and accepted history remain fixed while wall speeds may change.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare_with_walls(
        &mut self,
        config: LagrangianSgsConfig,
        time: f64,
        density: &[f64],
        velocity: &[[f64; 3]],
        wall_velocity: [[f64; 3]; 2],
        cancelled: impl FnMut() -> bool,
    ) -> Result<LagrangianSgsTrial<'_>, DynamicSgsError> {
        self.prepare_scaled_with_walls(
            config,
            time,
            density,
            velocity,
            Some(wall_velocity),
            cancelled,
        )
        .map(|(trial, _)| trial)
    }

    // Joint scalar preparation reuses these exact current Favre fields and
    // scales. The standalone momentum path retains its original arithmetic.
    fn prepare_scaled(
        &mut self,
        config: LagrangianSgsConfig,
        time: f64,
        density: &[f64],
        velocity: &[[f64; 3]],
        cancelled: impl FnMut() -> bool,
    ) -> Result<(LagrangianSgsTrial<'_>, InputScales), DynamicSgsError> {
        self.prepare_scaled_with_walls(config, time, density, velocity, None, cancelled)
    }

    #[allow(clippy::too_many_arguments)]
    fn prepare_scaled_with_walls(
        &mut self,
        config: LagrangianSgsConfig,
        time: f64,
        density: &[f64],
        velocity: &[[f64; 3]],
        walls: Option<[[f64; 3]; 2]>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<(LagrangianSgsTrial<'_>, InputScales), DynamicSgsError> {
        self.prepare_scaled_at(config, time, density, velocity, walls, false, cancelled)
    }

    // Native owners may refresh a constitutive field after a same-time wall
    // or force change. With zero elapsed time the accepted contractions do
    // not advect or age; only their emission uses the current physical state.
    // Public trajectory preparation continues to require increasing times.
    #[allow(clippy::too_many_arguments)]
    fn prepare_scaled_at(
        &mut self,
        config: LagrangianSgsConfig,
        time: f64,
        density: &[f64],
        velocity: &[[f64; 3]],
        walls: Option<[[f64; 3]; 2]>,
        allow_refresh: bool,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(LagrangianSgsTrial<'_>, InputScales), DynamicSgsError> {
        poll(&mut cancelled, 0)?;
        config.validate()?;
        if let Some(memory) = self.adaptive_memory {
            memory.validate_upper_bound(config.memory_time)?;
        }
        if !time.is_finite()
            || time < 0.0
            || self.accepted_time.is_some_and(|previous| {
                time < previous || (!allow_refresh && time == previous)
            })
        {
            return Err(DynamicSgsError::Invalid(
                "finite increasing Lagrangian sample times required",
            ));
        }
        let walls = match (self.normal_axis(), walls) {
            (None, None) => [[0.0; 3]; 2],
            (Some(_), Some(walls)) => walls,
            _ => {
                return Err(DynamicSgsError::Invalid(
                    "Lagrangian wall data must match workspace geometry",
                ));
            }
        };
        let elapsed = self.accepted_time.map_or(0.0, |previous| time - previous);
        let (scales, _) =
            self.flow
                .prepare_contractions(density, velocity, walls, &mut cancelled)?;
        let current_log = (scales.velocity > 0.0
            && (self.accepted_time.is_none() || elapsed > 0.0))
            .then(|| 2.0 * det::ln(scales.density) + 4.0 * det::ln(scales.velocity));
        let update = HistoryUpdate::new(
            &self.flow,
            self.accepted_time.is_some(),
            self.history_log_scale,
            current_log,
            elapsed,
            // The adaptive path computes its own per-cell weights. Admit its
            // fastest possible relaxation here: an enormous upper time bound
            // must not reject a representable local adaptive update.
            self.adaptive_memory.map_or(config.memory_time, |memory| memory.minimum_time),
            config.maximum_departure_cells,
        )?;
        let adaptive = self.adaptive_memory.map(|memory| {
            AdaptiveUpdate::new(
                memory,
                config.memory_time,
                self.flow.filter_width,
                self.history_log_scale,
                current_log,
                &update,
            )
        });
        let mut report = DynamicSgsReport {
            minimum_raw_coefficient_squared: f64::INFINITY,
            ..DynamicSgsReport::default()
        };
        let mut amplitude: f64 = 0.0;
        let mut memory_time_range = [f64::INFINITY, 0.0_f64];
        for cell in 0..self.history.len() {
            poll(&mut cancelled, cell)?;
            let (pair, memory_time) = if let Some(adaptive) = &adaptive {
                adaptive.pair(
                    &update,
                    cell,
                    velocity[cell],
                    density[cell],
                    &self.history,
                    &self.flow.contraction,
                )?
            } else {
                (
                    update.pair(cell, velocity[cell], &self.history, &self.flow.contraction)?,
                    config.memory_time,
                )
            };
            memory_time_range[0] = memory_time_range[0].min(memory_time);
            memory_time_range[1] = memory_time_range[1].max(memory_time);
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
        let common_log = adaptive.as_ref().map_or(update.common_log, |a| a.common_log());
        let update = HistoryUpdate { common_log, ..update };
        let log_scale = update.normalize(&mut self.candidate, amplitude, &mut cancelled)?;
        poll(&mut cancelled, 0)?;
        Ok((
            LagrangianSgsTrial {
                owner: self,
                time,
                log_scale,
                report,
                memory_time_range,
            },
            scales,
        ))
    }
}

// Both stress and scalar contractions use the same pathline transport and
// logarithmic history arithmetic. Only their physical normalization differs.
struct HistoryUpdate {
    shape: [usize; 3],
    geometry: FilterGeometry,
    spacing: [f64; 3],
    initialized: bool,
    elapsed: f64,
    maximum_departure_cells: f64,
    common_log: f64,
    old_scale: f64,
    new_scale: f64,
    old_weight: f64,
    fresh: f64,
}

impl HistoryUpdate {
    fn new(
        flow: &PeriodicDynamicSgs,
        initialized: bool,
        previous_log: Option<f64>,
        current_log: Option<f64>,
        elapsed: f64,
        memory_time: f64,
        maximum_departure_cells: f64,
    ) -> Result<Self, DynamicSgsError> {
        // A same-time constitutive refresh has no current-source weight.
        // Ignore that source's scale instead of risking loss of accepted
        // history against an arbitrarily larger, but unused, current source.
        let current_log = if initialized && elapsed == 0.0 { None } else { current_log };
        let common_log = match (previous_log, current_log) {
            (Some(a), Some(b)) => a.max(b),
            (Some(a), None) | (None, Some(a)) => a,
            (None, None) => 0.0,
        };
        let old_scale = previous_log.map_or(0.0, |value| det::exp(value - common_log));
        let new_scale = current_log.map_or(0.0, |value| det::exp(value - common_log));
        let relaxation = elapsed / memory_time;
        let fresh = -det::expm1(-relaxation);
        // Subtracting fresh from one loses representable old history when
        // elapsed/T exceeds about 37. Evaluate the exponential independently.
        let old_weight = det::exp(-relaxation);
        if initialized && elapsed > 0.0 && fresh == 0.0 {
            return Err(DynamicSgsError::Invalid(
                "unrepresentable Lagrangian relaxation interval",
            ));
        }
        Ok(Self {
            shape: flow.shape,
            geometry: flow.geometry,
            spacing: flow.spacing,
            initialized,
            elapsed,
            maximum_departure_cells,
            common_log,
            old_scale,
            new_scale,
            old_weight,
            fresh,
        })
    }

    fn pair(
        &self,
        cell: usize,
        velocity: [f64; 3],
        history: &[[f64; 2]],
        contraction: &[[f64; 2]],
    ) -> Result<[f64; 2], DynamicSgsError> {
        let pair = if !self.initialized {
            let mut pair = [0.0; 2];
            let mut determined = false;
            let samples = self.geometry.samples();
            self.geometry.each_box(self.shape, cell, |j| {
                determined |= contraction[j][1] != 0.0;
                for k in 0..2 {
                    pair[k] += contraction[j][k] / samples;
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
            let advected = self.advected(cell, velocity, history)?;
            let current = contraction[cell];
            if (self.old_weight > 0.0 && self.old_scale == 0.0 && advected[1] != 0.0)
                || (self.fresh > 0.0 && self.new_scale == 0.0 && current[1] != 0.0)
            {
                return Err(DynamicSgsError::Cell {
                    cell,
                    reason: "unrepresentable Lagrangian history scale separation",
                });
            }
            core::array::from_fn(|k| {
                self.old_weight * (self.old_scale * advected[k])
                    + self.fresh * (self.new_scale * current[k])
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
        Ok(pair)
    }

    fn advected(
        &self,
        cell: usize,
        velocity: [f64; 3],
        history: &[[f64; 2]],
    ) -> Result<[f64; 2], DynamicSgsError> {
        let departure =
            core::array::from_fn(|axis| velocity[axis] * (self.elapsed / self.spacing[axis]));
        if departure
            .iter()
            .any(|d: &f64| !d.is_finite() || d.abs() > self.maximum_departure_cells)
        {
            return Err(DynamicSgsError::Cell {
                cell,
                reason: "Lagrangian departure exceeds declared cell limit",
            });
        }
        sample_history(self.geometry, self.shape, cell, departure, history)
    }

    fn normalize(
        &self,
        candidate: &mut [[f64; 2]],
        amplitude: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Option<f64>, DynamicSgsError> {
        // A formerly strong flow must not underflow its decaying memory while
        // weaker new eddies form. Normalize only the candidate, never history.
        if amplitude == 0.0 {
            return Ok(None);
        }
        for (cell, pair) in candidate.iter_mut().enumerate() {
            poll(cancelled, cell)?;
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
        Ok(Some(self.common_log + det::ln(amplitude)))
    }
}

// Positive trilinear weights preserve the sign of the denominator. The
// physical velocity, including uniform advection, defines the departure;
// subtracting the covariance's reference velocity here would pin eddies.
fn sample_history(
    geometry: FilterGeometry,
    shape: [usize; 3],
    cell: usize,
    departure: [f64; 3],
    history: &[[f64; 2]],
) -> Result<[f64; 2], DynamicSgsError> {
    let p = coordinates(shape, cell);
    let normal = match geometry {
        FilterGeometry::Periodic => None,
        FilterGeometry::Channel(axis) => {
            let at = p[axis] as f64 - departure[axis];
            if !at.is_finite() || at < -0.5 || at > shape[axis] as f64 - 0.5 {
                return Err(DynamicSgsError::Cell {
                    cell,
                    reason: "Lagrangian departure crosses impermeable channel wall",
                });
            }
            Some(axis)
        }
    };
    let coordinate = core::array::from_fn::<_, 3, _>(|a| {
        let at = p[a] as f64 - departure[a];
        if normal == Some(a) {
            at.clamp(0.0, (shape[a] - 1) as f64)
        } else {
            at.rem_euclid(shape[a] as f64)
        }
    });
    let lower = coordinate.map(|v| v.floor() as usize);
    let upper = core::array::from_fn::<_, 3, _>(|a| {
        if normal == Some(a) {
            (lower[a] + 1).min(shape[a] - 1)
        } else {
            (lower[a] + 1) % shape[a]
        }
    });
    let fraction = core::array::from_fn::<_, 3, _>(|a| coordinate[a] - lower[a] as f64);
    // rem_euclid can round a tiny negative periodic coordinate to the extent.
    // Retain the original wrapped index even in that endpoint case.
    let lower = core::array::from_fn::<_, 3, _>(|a| lower[a] % shape[a]);
    let mut result = [0.0; 2];
    for z in 0..2 {
        for y in 0..2 {
            for x in 0..2 {
                let offset = [x, y, z];
                let q = core::array::from_fn(|a| if offset[a] == 0 { lower[a] } else { upper[a] });
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
    Ok(result)
}

#[cfg(test)]
#[path = "lagrangian/tests.rs"]
mod tests;
