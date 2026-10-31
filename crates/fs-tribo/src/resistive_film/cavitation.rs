//! Transient Elrod--Adams liquid content on the original Reynolds-film graph.
//!
//! Each cell retains `m = rho * area * gap * theta`. Absolute pressure and
//! liquid content satisfy `p >= p_c`, `0 <= theta <= 1` and
//! `(p - p_c) * (1 - theta) = 0`. Backward Euler conserves this mass through
//! shared Poiseuille fluxes and upwinded, caller-supplied Couette entrainment.
//! Rupture and reformation therefore change liquid content, not just pressure.
//!
//! This is the isothermal, incompressible-liquid, negligible-void-mass thin-film
//! approximation. `theta` is retained lubricant content, not vapor mass. Bubble
//! inertia, finite vapor density, latent heat, dissolved-gas transport, elastic
//! deformation and thermal feedback need their own owners. In particular, a
//! fully filled sealed pocket cannot be compressed by this model.
//!
//! Reference: H. G. Elrod, "A Cavitation Algorithm", Journal of Lubrication
//! Technology 103 (1981), 350--354, DOI 10.1115/1.3251669. This implementation
//! uses pressure/content complementarity on the existing finite-volume graph;
//! it does not infer a physical cavitation pressure from a fluid name.

use super::{FilmError, MAX_CELLS, MAX_CHANNELS, MAX_PORTS, ResistiveFilm, dot, fail};
use std::sync::Arc;

/// Explicit frozen liquid and reservoir data, all in SI units.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CavitationLiquid {
    /// Constant liquid density [kg m^-3]. Void mass is neglected explicitly.
    pub density_kg_m3: f64,
    /// Absolute rupture/reformation pressure [Pa], supplied by a material law.
    pub cavitation_pressure_pa: f64,
    /// Absolute pressure [Pa] at every `FilmChannel::to == None` reservoir.
    pub ambient_pressure_pa: f64,
    /// Reservoir liquid content for incoming Couette transport. Must be one
    /// above cavitation pressure; partial filling is admitted only at `p_c`.
    pub ambient_liquid_fraction: f64,
}

/// Work cap and tolerances for the actual cell and global liquid-mass balances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CavitationSolveConfig {
    /// Maximum deterministic active-set solves; must lie in 1..=512.
    pub maximum_iterations: usize,
    /// Absolute tolerance [kg] on each cell and the complete graph.
    pub absolute_mass_tolerance_kg: f64,
    /// Relative tolerance on stored/transported liquid mass, in [0, 1).
    pub relative_mass_tolerance: f64,
}

/// Accepted configuration, content and absolute pressure. Creating a new
/// candidate never changes this state; the surrounding assembly owns commit.
#[derive(Clone, Debug)]
pub struct CavitatingFilmState {
    owner: Arc<()>,
    configuration: Vec<f64>,
    liquid_fraction: Vec<f64>,
    pressure_pa: Vec<f64>,
}
impl PartialEq for CavitatingFilmState {
    fn eq(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.owner, &other.owner)
            && self.configuration == other.configuration
            && self.liquid_fraction == other.liquid_fraction
            && self.pressure_pa == other.pressure_pa
    }
}
impl CavitatingFilmState {
    /// Mechanical configuration used to compute the retained cell volumes.
    #[must_use]
    pub fn configuration(&self) -> &[f64] {
        &self.configuration
    }
    /// Liquid volume/content fractions; a void has zero retained liquid mass.
    #[must_use]
    pub fn liquid_fraction(&self) -> &[f64] {
        &self.liquid_fraction
    }
    /// Absolute cell pressure [Pa].
    #[must_use]
    pub fn pressure_pa(&self) -> &[f64] {
        &self.pressure_pa
    }
}

/// Direct conservation diagnostics. This is not a thermal or work receipt.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CavitationReport {
    /// Retained liquid before the proposed step [kg].
    pub liquid_mass_before_kg: f64,
    /// Retained liquid in the candidate state [kg].
    pub liquid_mass_after_kg: f64,
    /// Signed outward reservoir mass over this step [kg].
    pub boundary_mass_out_kg: f64,
    /// `after - before + boundary_out`, independently rebuilt from fluxes [kg].
    pub mass_residual_kg: f64,
    /// Largest absolute independently reconstructed cell mass residual [kg].
    pub maximum_cell_mass_residual_kg: f64,
    /// Smallest candidate cell liquid content, in [0, 1].
    pub minimum_liquid_fraction: f64,
    /// Largest candidate absolute cell pressure [Pa].
    pub maximum_pressure_pa: f64,
    /// Number of active-set linear solves used by this candidate.
    pub iterations: usize,
}

/// An immutable candidate. The same face flux appears with opposite signs in
/// adjacent cells. Forces are endpoint gauge-pressure projections, not average
/// forces or an assertion of time-integrated mechanical passivity.
#[derive(Clone, Debug, PartialEq)]
pub struct CavitationStep {
    /// Candidate state; the caller commits it after its complete step accepts.
    pub state: CavitatingFilmState,
    /// Signed total liquid-volume flux from `channel.from` towards `to` [m^3/s].
    pub channel_volume_flow_m3_s: Vec<f64>,
    /// Positive resisting force in the original work-conjugate port basis.
    pub generalized_force: Vec<f64>,
    /// Direct liquid-mass accounting for this candidate.
    pub report: CavitationReport,
}

/// Cavitating liquid storage and transport on the existing geometric owner.
#[derive(Clone, Debug)]
pub struct CavitatingFilm {
    owner: Arc<()>,
    geometry: ResistiveFilm,
    liquid: CavitationLiquid,
    solve: CavitationSolveConfig,
}

impl CavitatingFilm {
    /// Reuse the original cells, channels, gap map, viscosity and validity
    /// bounds. No pressure floor, gap floor or hidden leakage is introduced.
    pub fn new(
        geometry: ResistiveFilm,
        liquid: CavitationLiquid,
        solve: CavitationSolveConfig,
    ) -> Result<Self, FilmError> {
        if !liquid.density_kg_m3.is_finite()
            || liquid.density_kg_m3 <= 0.0
            || !liquid.cavitation_pressure_pa.is_finite()
            || liquid.cavitation_pressure_pa <= 0.0
            || !liquid.ambient_pressure_pa.is_finite()
            || liquid.ambient_pressure_pa < liquid.cavitation_pressure_pa
            || !liquid.ambient_liquid_fraction.is_finite()
            || !(0.0..=1.0).contains(&liquid.ambient_liquid_fraction)
            || (liquid.ambient_pressure_pa > liquid.cavitation_pressure_pa
                && liquid.ambient_liquid_fraction != 1.0)
            || liquid.ambient_pressure_pa - liquid.cavitation_pressure_pa
                > geometry.limits.maximum_pressure_pa
            || !(1..=512).contains(&solve.maximum_iterations)
            || !solve.absolute_mass_tolerance_kg.is_finite()
            || solve.absolute_mass_tolerance_kg < 0.0
            || !solve.relative_mass_tolerance.is_finite()
            || !(0.0..1.0).contains(&solve.relative_mass_tolerance)
            || (solve.absolute_mass_tolerance_kg == 0.0 && solve.relative_mass_tolerance == 0.0)
        {
            return fail(
                "invalid cavitating liquid, absolute reservoir pressure, filling or solve budget",
            );
        }
        Ok(Self {
            owner: Arc::new(()),
            geometry,
            liquid,
            solve,
        })
    }

    /// The original immutable pressure graph and viscosity coefficient.
    #[must_use]
    pub const fn geometry(&self) -> &ResistiveFilm {
        &self.geometry
    }
    /// Frozen density, cavitation pressure and explicit reservoir conditions.
    #[must_use]
    pub const fn liquid(&self) -> CavitationLiquid {
        self.liquid
    }
    /// Number of retained liquid-content cells.
    #[must_use]
    pub fn cell_count(&self) -> usize {
        self.geometry.cell_count()
    }
    /// Number of transport faces, including exterior reservoir faces.
    #[must_use]
    pub fn channel_count(&self) -> usize {
        self.geometry.channels.len()
    }

    /// Full liquid at the declared reservoir pressure. An initial state is an
    /// initial condition, not a solution of steady transport or wall motion.
    pub fn initial_state(&self, q: &[f64]) -> Result<CavitatingFilmState, FilmError> {
        self.initialize(q, &[1.0; MAX_CELLS][..self.cell_count()])
    }

    /// Explicit initial starvation. Partial cells start at `p_c`; full cells
    /// start at ambient pressure. No liquid mass is inferred for missing volume.
    pub fn initialize(
        &self,
        q: &[f64],
        liquid_fraction: &[f64],
    ) -> Result<CavitatingFilmState, FilmError> {
        self.geometry.validate_configuration(q)?;
        if liquid_fraction.len() != self.cell_count()
            || liquid_fraction
                .iter()
                .any(|x| !x.is_finite() || !(0.0..=1.0).contains(x))
        {
            return fail("invalid initial cavitating-film liquid fractions");
        }
        let state = CavitatingFilmState {
            owner: Arc::clone(&self.owner),
            configuration: q.to_vec(),
            liquid_fraction: liquid_fraction.to_vec(),
            pressure_pa: liquid_fraction
                .iter()
                .map(|x| {
                    if *x < 1.0 {
                        self.liquid.cavitation_pressure_pa
                    } else {
                        self.liquid.ambient_pressure_pa
                    }
                })
                .collect(),
        };
        self.validate_state(&state)?;
        Ok(state)
    }

    fn validate_state(&self, state: &CavitatingFilmState) -> Result<(), FilmError> {
        if !Arc::ptr_eq(&self.owner, &state.owner) {
            return fail("cavitating-film state belongs to another geometry or liquid owner");
        }
        self.geometry.validate_configuration(&state.configuration)?;
        if state.liquid_fraction.len() != self.cell_count()
            || state.pressure_pa.len() != self.cell_count()
        {
            return fail("cavitating-film state belongs to another cell layout");
        }
        for (&theta, &p) in state.liquid_fraction.iter().zip(&state.pressure_pa) {
            if !theta.is_finite()
                || !(0.0..=1.0).contains(&theta)
                || !p.is_finite()
                || p < self.liquid.cavitation_pressure_pa
                || (p - self.liquid.ambient_pressure_pa).abs()
                    > self.geometry.limits.maximum_pressure_pa
                || (theta < 1.0 && p != self.liquid.cavitation_pressure_pa)
            {
                return fail("invalid cavitating-film pressure/content complementarity");
            }
        }
        Ok(())
    }

    /// Advance prescribed geometry and transport through one backward-Euler
    /// step. `entrainment_m3_s[k]` is the signed FULL-film Couette flux on the
    /// original channel; the solver multiplies it by the upwind liquid content.
    /// The caller owns gap/velocity quadrature and the no-slip thin-film regime.
    ///
    /// The cancellation callback is polled at assembly, active-set and dense
    /// elimination boundaries. Every refusal leaves the input state untouched.
    /// All work is bounded by the existing graph limits and the explicit solve
    /// cap; iterations allocate nothing. Retained outputs allocate once.
    pub fn advance(
        &self,
        state: &CavitatingFilmState,
        q_next: &[f64],
        entrainment_m3_s: &[f64],
        dt_s: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CavitationStep, FilmError> {
        self.advance_kernel(state, q_next, entrainment_m3_s, dt_s, &mut cancelled)
    }

    // Keep one numerical kernel for all caller cancellation closures. Dynamic
    // dispatch occurs only at bounded solver boundaries, never per matrix entry.
    #[allow(clippy::too_many_lines, clippy::large_stack_arrays)]
    fn advance_kernel(
        &self,
        state: &CavitatingFilmState,
        q_next: &[f64],
        entrainment_m3_s: &[f64],
        dt_s: f64,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<CavitationStep, FilmError> {
        poll(cancelled)?;
        self.validate_state(state)?;
        self.geometry.validate_configuration(q_next)?;
        if !dt_s.is_finite()
            || dt_s <= 0.0
            || entrainment_m3_s.len() != self.channel_count()
            || entrainment_m3_s.iter().any(|x| !x.is_finite())
        {
            return fail("invalid cavitating-film duration or channel entrainment");
        }
        let n = self.cell_count();
        let mut a = [0.0; MAX_CELLS * MAX_CELLS];
        let mut b = [0.0; MAX_CELLS * MAX_CELLS];
        let mut rhs = [0.0; MAX_CELLS];
        let mut volume = [0.0; MAX_CELLS];
        let mut old_liquid_volume = [0.0; MAX_CELLS];
        let mut conductance = [0.0; MAX_CHANNELS];
        let mut full = [false; MAX_CELLS];
        for (i, cell) in self.geometry.cells.iter().enumerate() {
            let h0 = cell.gap.reference_m - dot(&cell.gap.closure, state.configuration());
            let h1 = cell.gap.reference_m - dot(&cell.gap.closure, q_next);
            volume[i] = cell.area_m2 * h1;
            old_liquid_volume[i] = cell.area_m2 * h0 * state.liquid_fraction[i];
            if !volume[i].is_finite() || volume[i] <= 0.0 || !old_liquid_volume[i].is_finite() {
                return fail("unrepresentable cavitating-film volume");
            }
            b[i * n + i] = volume[i];
            rhs[i] = old_liquid_volume[i];
            full[i] = old_liquid_volume[i] > volume[i]
                || (old_liquid_volume[i] == volume[i]
                    && state.pressure_pa[i] > self.liquid.cavitation_pressure_pa);
        }
        let boundary_excess = self.liquid.ambient_pressure_pa - self.liquid.cavitation_pressure_pa;
        for (k, (edge, factor)) in self
            .geometry
            .channels
            .iter()
            .zip(&self.geometry.factors)
            .enumerate()
        {
            let h = edge.gap.reference_m - dot(&edge.gap.closure, q_next);
            let g = if h > 0.0 { factor * h * h * h } else { 0.0 };
            if !g.is_finite() || (h > 0.0 && g == 0.0) || (h <= 0.0 && entrainment_m3_s[k] != 0.0) {
                return fail(
                    "invalid cavitating-film conductance or entrainment through a closed channel",
                );
            }
            conductance[k] = g;
            let ag = dt_s * g;
            let e = dt_s * entrainment_m3_s[k];
            let i = edge.from;
            a[i * n + i] += ag;
            if let Some(j) = edge.to {
                a[j * n + j] += ag;
                a[i * n + j] -= ag;
                a[j * n + i] -= ag;
                let donor = if e >= 0.0 { i } else { j };
                b[i * n + donor] += e;
                b[j * n + donor] -= e;
            } else {
                rhs[i] += ag * boundary_excess;
                if e >= 0.0 {
                    b[i * n + i] += e;
                } else {
                    rhs[i] -= e * self.liquid.ambient_liquid_fraction;
                }
            }
        }
        if a[..n * n]
            .iter()
            .chain(&b[..n * n])
            .chain(&rhs[..n])
            .any(|x| !x.is_finite())
        {
            return fail("cavitating-film transport assembly overflow");
        }
        let mut pressure_excess = [0.0; MAX_CELLS];
        let mut theta = [0.0; MAX_CELLS];
        for iteration in 1..=self.solve.maximum_iterations {
            poll(cancelled)?;
            let mut matrix = [0.0; MAX_CELLS * MAX_CELLS];
            let mut x = rhs;
            let mut pressure_scale = [0.0; MAX_CELLS];
            for j in 0..n {
                if full[j] {
                    if a[j * n + j] <= 0.0 {
                        return fail("fully filled undrained film requires liquid compressibility");
                    }
                    pressure_scale[j] = b[j * n + j] / a[j * n + j];
                    if !pressure_scale[j].is_finite() || pressure_scale[j] <= 0.0 {
                        return fail("unrepresentable cavitating-film pressure scale");
                    }
                }
                for i in 0..n {
                    if full[j] {
                        matrix[i * n + j] = (a[i * n + j] / a[j * n + j]) * b[j * n + j];
                        x[i] -= b[i * n + j];
                    } else {
                        matrix[i * n + j] = b[i * n + j];
                    }
                }
            }
            solve_local(&mut matrix, &mut x, n, cancelled)?;
            let mut changed = false;
            let mut negative_mass = false;
            // Normalize only machine-roundoff at a complementarity boundary;
            // independently rebuilt mass balances still have to meet the caller
            // tolerance. This does not clamp an unconstrained pressure solve.
            let roundoff = 32.0 * n as f64 * f64::EPSILON;
            for i in 0..n {
                if full[i] {
                    if x[i] < -roundoff {
                        full[i] = false;
                        changed = true;
                    }
                    pressure_excess[i] = x[i].max(0.0) * pressure_scale[i];
                    theta[i] = 1.0;
                } else {
                    if x[i] > 1.0 + roundoff {
                        full[i] = true;
                        changed = true;
                    }
                    negative_mass |= x[i] < -roundoff;
                    theta[i] = x[i].clamp(0.0, 1.0);
                    pressure_excess[i] = 0.0;
                }
            }
            if changed {
                continue;
            }
            if negative_mass {
                return fail("cavitating-film solve produced negative liquid mass");
            }
            return self.finish(
                q_next,
                entrainment_m3_s,
                dt_s,
                iteration,
                &volume,
                &old_liquid_volume,
                &conductance,
                &theta,
                &pressure_excess,
                cancelled,
            );
        }
        fail("cavitating-film active-set iteration budget exhausted")
    }

    #[allow(clippy::too_many_arguments)]
    fn finish(
        &self,
        q_next: &[f64],
        entrainment: &[f64],
        dt: f64,
        iterations: usize,
        volume: &[f64; MAX_CELLS],
        old_volume: &[f64; MAX_CELLS],
        conductance: &[f64; MAX_CHANNELS],
        theta: &[f64; MAX_CELLS],
        excess: &[f64; MAX_CELLS],
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<CavitationStep, FilmError> {
        let n = self.cell_count();
        let rho = self.liquid.density_kg_m3;
        let mut residual = [0.0; MAX_CELLS];
        let mut scales = [0.0; MAX_CELLS];
        let mut forces = [0.0; MAX_PORTS];
        let mut flows = [0.0; MAX_CHANNELS];
        let mut boundary = 0.0;
        let mut boundary_scale = 0.0;
        let mut before = 0.0;
        let mut after = 0.0;
        for i in 0..n {
            residual[i] = volume[i] * theta[i] - old_volume[i];
            scales[i] = (volume[i] * theta[i]).abs().max(old_volume[i].abs());
            before += old_volume[i];
            after += volume[i] * theta[i];
        }
        let boundary_excess = self.liquid.ambient_pressure_pa - self.liquid.cavitation_pressure_pa;
        for (k, edge) in self.geometry.channels.iter().enumerate() {
            let donor = if entrainment[k] >= 0.0 {
                theta[edge.from]
            } else {
                edge.to
                    .map_or(self.liquid.ambient_liquid_fraction, |j| theta[j])
            };
            let flux = conductance[k]
                * (excess[edge.from] - edge.to.map_or(boundary_excess, |j| excess[j]))
                + entrainment[k] * donor;
            flows[k] = flux;
            let amount = dt * flux;
            residual[edge.from] += amount;
            scales[edge.from] += amount.abs();
            if let Some(j) = edge.to {
                residual[j] -= amount;
                scales[j] += amount.abs();
            } else {
                boundary += amount;
                boundary_scale += amount.abs();
            }
        }
        let mut max_residual = 0.0_f64;
        for i in 0..n {
            let r = rho * residual[i];
            let tolerance = self.solve.absolute_mass_tolerance_kg
                + self.solve.relative_mass_tolerance * rho * scales[i];
            if !r.is_finite() || !tolerance.is_finite() || r.abs() > tolerance {
                return fail("cavitating-film cell liquid-mass balance exceeds tolerance");
            }
            max_residual = max_residual.max(r.abs());
            let gauge = excess[i] - boundary_excess;
            if !gauge.is_finite() || gauge.abs() > self.geometry.limits.maximum_pressure_pa {
                return fail("cavitating-film pressure exceeds declared validity domain");
            }
            for (j, weight) in self.geometry.cells[i].gap.closure.iter().enumerate() {
                forces[j] += self.geometry.cells[i].area_m2 * gauge * weight;
            }
        }
        let mass_residual = rho * (after - before + boundary);
        let tolerance = self.solve.absolute_mass_tolerance_kg
            + self.solve.relative_mass_tolerance
                * rho
                * (before.abs() + after.abs() + boundary_scale);
        if !mass_residual.is_finite()
            || !tolerance.is_finite()
            || mass_residual.abs() > tolerance
            || forces[..self.geometry.port_count()]
                .iter()
                .any(|x| !x.is_finite())
            || flows[..self.channel_count()].iter().any(|x| !x.is_finite())
            || !((before + after + boundary.abs()) * rho).is_finite()
        {
            return fail("cavitating-film global liquid-mass or force closure failed");
        }
        poll(cancelled)?;
        let candidate = CavitatingFilmState {
            owner: Arc::clone(&self.owner),
            configuration: q_next.to_vec(),
            liquid_fraction: theta[..n].to_vec(),
            pressure_pa: excess[..n]
                .iter()
                .map(|p| p + self.liquid.cavitation_pressure_pa)
                .collect(),
        };
        self.validate_state(&candidate)?;
        let report = CavitationReport {
            liquid_mass_before_kg: rho * before,
            liquid_mass_after_kg: rho * after,
            boundary_mass_out_kg: rho * boundary,
            mass_residual_kg: mass_residual,
            maximum_cell_mass_residual_kg: max_residual,
            minimum_liquid_fraction: theta[..n].iter().copied().fold(1.0, f64::min),
            maximum_pressure_pa: candidate.pressure_pa.iter().copied().fold(0.0, f64::max),
            iterations,
        };
        Ok(CavitationStep {
            state: candidate,
            channel_volume_flow_m3_s: flows[..self.channel_count()].to_vec(),
            generalized_force: forces[..self.geometry.port_count()].to_vec(),
            report,
        })
    }
}

fn poll(cancelled: &mut dyn FnMut() -> bool) -> Result<(), FilmError> {
    if cancelled() {
        fail("cavitating-film solve cancelled")
    } else {
        Ok(())
    }
}

/// Bounded pressure-island elimination. Upwind entrainment makes the selected
/// columns nonsymmetric, so the parent full-film Cholesky cannot be reused.
/// Row equilibration and deterministic partial pivoting introduce no jitter.
#[allow(clippy::needless_range_loop)]
fn solve_local(
    matrix: &mut [f64; MAX_CELLS * MAX_CELLS],
    rhs: &mut [f64; MAX_CELLS],
    n: usize,
    cancelled: &mut dyn FnMut() -> bool,
) -> Result<(), FilmError> {
    for i in 0..n {
        let scale = matrix[i * n..(i + 1) * n]
            .iter()
            .map(|x| x.abs())
            .fold(0.0, f64::max);
        if !scale.is_finite() || scale <= 0.0 {
            return fail("singular undrained cavitating-film pressure island");
        }
        for j in 0..n {
            matrix[i * n + j] /= scale;
        }
        rhs[i] /= scale;
    }
    for k in 0..n {
        poll(cancelled)?;
        let mut pivot = k;
        for i in k + 1..n {
            if matrix[i * n + k].abs() > matrix[pivot * n + k].abs() {
                pivot = i;
            }
        }
        if !matrix[pivot * n + k].is_finite() || matrix[pivot * n + k].abs() <= 64.0 * f64::EPSILON
        {
            return fail("singular or unresolved saturated film island requires compressibility");
        }
        if pivot != k {
            for j in 0..n {
                matrix.swap(k * n + j, pivot * n + j);
            }
            rhs.swap(k, pivot);
        }
        for i in k + 1..n {
            let factor = matrix[i * n + k] / matrix[k * n + k];
            matrix[i * n + k] = 0.0;
            for j in k + 1..n {
                matrix[i * n + j] -= factor * matrix[k * n + j];
            }
            rhs[i] -= factor * rhs[k];
        }
    }
    for i in (0..n).rev() {
        for j in i + 1..n {
            rhs[i] -= matrix[i * n + j] * rhs[j];
        }
        rhs[i] /= matrix[i * n + i];
    }
    if rhs[..n].iter().any(|x| !x.is_finite()) {
        return fail("nonfinite cavitating-film pressure solution");
    }
    Ok(())
}

#[cfg(test)]
mod tests;
