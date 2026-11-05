//! Barotropic liquid storage on the original Elrod--Adams film graph.
//!
//! The authored constant bulk modulus B satisfies rho*dp/drho=B. Density is
//! referenced to the native ambient pressure, so rho_c=rho_ref*exp((pc-pref)/B).
//! The retained content u=rho(p)*theta/rho_c is fractional filling below one
//! and compressed full liquid above one: p=pc+B*ln(u), theta=min(u,1).
//! See Brewe and Jacobson, NASA TM-88826 (1986), background theory eqs. 1--5:
//! https://ntrs.nasa.gov/api/citations/19870002529/downloads/19870002529.pdf.
//!
//! Storage and upwind entrainment are linear in u. Pressure mass transport
//! integrates rho/mu along each original constant-gap passage, including its
//! native common Barus law. With x=u-1 >=0 in full liquid,
//! Phi(x)=B*expm1((1-alpha*B)*ln(1+x))/(1-alpha*B). Its differential is
//! B*exp(-alpha*(p-pc)); for alpha=0 the pressure operator is exactly linear
//! in x. Native conductance already includes exp(alpha*(pref-pc)).
//!
//! Compression x and partial filling theta are retained separately during the
//! solve: small physical compressibility is never rounded away by forming 1+x.
//! One bounded semismooth Newton system reuses the native dense elimination,
//! compensated incidence sums, reservoirs, state identity and force projection.
//! This is an isothermal Newtonian thin-film reduction with negligible void
//! mass, not a caloric/phase-energy solver or a time-integrated work receipt.

use std::sync::Arc;

use fs_math::det;

use super::{
    CavitatingFilm, CavitatingFilmState, CavitationLiquid, CavitationReport, CavitationReservoir,
    FilmError, FlowBalance, MAX_CELLS, MAX_CHANNELS, MAX_PORTS, ResistiveFilm, dot, fail, poll,
    solve_local,
};

/// The native film with an explicitly authored positive isothermal bulk modulus.
/// The native liquid density is the density at its absolute ambient pressure.
#[derive(Clone, Debug)]
pub struct CompressibleCavitatingFilm {
    inventory: CavitatingFilm,
    bulk_modulus_pa: f64,
    density_at_cavitation_kg_m3: f64,
    exponent: f64,
    supplies: Arc<[Point]>,
}

/// An immutable mass-conserving candidate; accept its state only after the
/// surrounding assembly accepts. A compressible passage has one MASS flux;
/// dividing it by an arbitrary density would not define a unique volume flux.
#[derive(Clone, Debug, PartialEq)]
pub struct CompressibleCavitationStep {
    pub state: CavitatingFilmState,
    /// Signed liquid mass flow from each original `from` to `to` [kg/s].
    pub channel_mass_flow_kg_s: Vec<f64>,
    /// Endpoint normal gauge-pressure effort in the original port basis.
    pub generalized_force: Vec<f64>,
    /// Direct mass audit, including physical compressive storage.
    pub report: CavitationReport,
}

/// Endpoint material owner and its matching conservative inventory. Accept both
/// together after the surrounding thermal/mechanical proposal succeeds.
#[derive(Clone, Debug)]
pub struct CompressibleCavitationContinuation {
    pub film: CompressibleCavitatingFilm,
    pub step: CompressibleCavitationStep,
}

#[derive(Clone, Copy, Debug, Default)]
struct Point {
    theta: f64,
    compression: f64,
    pressure_excess: f64,
    differential: f64,
    fluidity: f64,
    full: bool,
}

impl CompressibleCavitatingFilm {
    /// Reuse the geometry, viscosities, Barus law, supplies and solve budget.
    /// A distinct state identity prevents an incompressible owner from treating
    /// this pressure-dependent inventory as constant-density liquid. Initialize
    /// through this owner; cloning or viscosity/supply refresh retains identity.
    pub fn new(inventory: CavitatingFilm, bulk_modulus_pa: f64) -> Result<Self, FilmError> {
        let mut film = Self::compile(inventory, bulk_modulus_pa)?;
        film.inventory.owner = Arc::new(());
        Ok(film)
    }

    // The native material-transition owner already assigns the correct new or
    // retained inventory identity. Compile coefficients without resetting it.
    fn compile(inventory: CavitatingFilm, bulk_modulus_pa: f64) -> Result<Self, FilmError> {
        if !bulk_modulus_pa.is_finite() || bulk_modulus_pa <= 0.0 {
            return fail("compressible cavitating film needs a positive finite bulk modulus");
        }
        let exponent = 1.0 - inventory.pressure_law.alpha * bulk_modulus_pa;
        let liquid = inventory.liquid;
        let density_at_cavitation_kg_m3 = liquid.density_kg_m3
            * det::exp(
                (liquid.cavitation_pressure_pa - liquid.ambient_pressure_pa) / bulk_modulus_pa,
            );
        if !exponent.is_finite()
            || !density_at_cavitation_kg_m3.is_finite()
            || density_at_cavitation_kg_m3 <= 0.0
        {
            return fail("unrepresentable compressible cavitating-film reference density or law");
        }
        let mut film = Self {
            inventory,
            bulk_modulus_pa,
            density_at_cavitation_kg_m3,
            exponent,
            supplies: Arc::from([]),
        };
        film.supplies = film.compile_supplies()?;
        Ok(film)
    }

    #[must_use]
    pub const fn bulk_modulus_pa(&self) -> f64 {
        self.bulk_modulus_pa
    }

    #[must_use]
    pub const fn density_at_cavitation_kg_m3(&self) -> f64 {
        self.density_at_cavitation_kg_m3
    }

    #[must_use]
    pub const fn geometry(&self) -> &ResistiveFilm {
        self.inventory.geometry()
    }

    /// Density is referenced to `ambient_pressure_pa`, not to every cell.
    #[must_use]
    pub const fn liquid(&self) -> CavitationLiquid {
        self.inventory.liquid()
    }

    pub fn reservoir_conditions(
        &self,
    ) -> impl ExactSizeIterator<Item = Option<CavitationReservoir>> + '_ {
        self.inventory.reservoir_conditions()
    }

    pub fn initial_state(&self, q: &[f64]) -> Result<CavitatingFilmState, FilmError> {
        self.inventory.initial_state(q)
    }

    pub fn initialize(
        &self,
        q: &[f64],
        fractions: &[f64],
    ) -> Result<CavitatingFilmState, FilmError> {
        self.inventory.initialize(q, fractions)
    }

    /// Explicit compressed initial inventory. Partial cells must be exactly at
    /// rupture; pressure and filling are independent data only for full liquid.
    /// This supplies an initial condition, not a steady-flow solution.
    pub fn initialize_with_pressure(
        &self,
        q: &[f64],
        fractions: &[f64],
        pressure_pa: &[f64],
    ) -> Result<CavitatingFilmState, FilmError> {
        let mut state = self.inventory.initialize(q, fractions)?;
        state.pressure_pa = pressure_pa.to_vec();
        self.inventory.validate_state(&state)?;
        for (&pressure, &theta) in state.pressure_pa.iter().zip(&state.liquid_fraction) {
            self.from_pressure(pressure, theta)?;
        }
        Ok(state)
    }

    /// Continue the same compressed inventory with updated native viscosity.
    /// Density and bulk modulus remain fixed; no heat or temperature is inferred.
    pub fn with_viscosity(&self, viscosity_pa_s: f64) -> Result<Self, FilmError> {
        let mut film = self.clone();
        film.inventory = self.inventory.with_viscosity(viscosity_pa_s)?;
        Ok(film)
    }

    pub fn with_channel_viscosities(&self, viscosities_pa_s: &[f64]) -> Result<Self, FilmError> {
        let mut film = self.clone();
        film.inventory = self.inventory.with_channel_viscosities(viscosities_pa_s)?;
        Ok(film)
    }

    /// Prescribe absolute supplies without resetting pressure or stored mass.
    pub fn with_reservoirs(
        &self,
        supplies: &[Option<CavitationReservoir>],
    ) -> Result<Self, FilmError> {
        let mut film = self.clone();
        film.inventory = self.inventory.with_reservoirs(supplies)?;
        film.supplies = film.compile_supplies()?;
        Ok(film)
    }

    fn compile_supplies(&self) -> Result<Arc<[Point]>, FilmError> {
        self.inventory
            .reservoir_conditions()
            .map(|supply| match supply {
                Some(s) => self.from_pressure(s.pressure_pa, s.liquid_fraction),
                None => Ok(Point::default()),
            })
            .collect::<Result<Vec<_>, _>>()
            .map(Arc::from)
    }

    fn from_pressure(&self, pressure: f64, theta: f64) -> Result<Point, FilmError> {
        let compression = if theta == 1.0 {
            det::expm1((pressure - self.liquid().cavitation_pressure_pa) / self.bulk_modulus_pa)
        } else {
            0.0
        };
        self.point(theta, compression, theta == 1.0)
    }

    fn point(&self, theta: f64, compression: f64, full: bool) -> Result<Point, FilmError> {
        let pressure_excess = self.bulk_modulus_pa * log1p(compression);
        let liquid = self.liquid();
        let gauge = pressure_excess - (liquid.ambient_pressure_pa - liquid.cavitation_pressure_pa);
        let attenuation = det::exp(-self.inventory.pressure_law.alpha * pressure_excess);
        let differential = self.bulk_modulus_pa * attenuation;
        let fluidity = self.inventory.pressure_law.fluidity_scale * attenuation;
        if !theta.is_finite()
            || !(0.0..=1.0).contains(&theta)
            || !compression.is_finite()
            || compression < 0.0
            || (theta < 1.0 && compression != 0.0)
            || !pressure_excess.is_finite()
            || !gauge.is_finite()
            || gauge.abs() > self.geometry().limits.maximum_pressure_pa
            || !differential.is_finite()
            || differential <= 0.0
            || !fluidity.is_finite()
            || fluidity <= 0.0
            || !(self.density_at_cavitation_kg_m3 * (1.0 + compression)).is_finite()
        {
            return fail("compressible cavitating-film state exceeds its physical domain");
        }
        Ok(Point {
            theta,
            compression,
            pressure_excess,
            differential,
            fluidity,
            full,
        })
    }

    /// Backward-Euler storage and conservative upwind mean-wall transport.
    /// `entrainment_m3_s` is the FULL-film wall-volume transport on each original
    /// oriented passage, multiplied by the upwind density and liquid fraction.
    /// Physical bulk modulus is never softened to make the solve easier.
    /// The native iteration cap bounds Newton/phase corrections; each correction
    /// has at most 32 domain-preserving backtracks. Iterations allocate nothing.
    pub fn advance(
        &self,
        state: &CavitatingFilmState,
        q_next: &[f64],
        entrainment_m3_s: &[f64],
        dt_s: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CompressibleCavitationStep, FilmError> {
        self.advance_kernel(self, state, q_next, entrainment_m3_s, dt_s, &mut cancelled)
    }

    /// Continue the ORIGINAL pressure-dependent mass with endpoint liquid data.
    ///
    /// The old density is reconstructed at the OLD reference and rupture
    /// pressures; new density, saturation pressure and channel viscosities enter
    /// the endpoint balance. New density is referenced to the NEW ambient
    /// pressure. B remains the independently authored isothermal bulk modulus.
    /// A decrease in reference density can therefore compress sealed full
    /// liquid instead of manufacturing or rejecting its original mass.
    ///
    /// Geometry, the Barus coefficient and explicit absolute supplies stay
    /// fixed. Native material admission rebuilds their pressure transforms and
    /// refuses supplies incompatible with the new rupture/reference domain.
    /// Changed liquid data receive a new inventory identity; unchanged liquid
    /// and viscosity-only changes retain compatible state identity. The new
    /// film and its candidate must be accepted together.
    ///
    /// This is a first-order prescribed endpoint constitutive split. The caller
    /// owns temperature, source validity and temporal refinement. No energy,
    /// thermal expansion law, phase heat or temperature-dependent B is inferred.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_liquid(
        &self,
        state: &CavitatingFilmState,
        q_next: &[f64],
        entrainment_m3_s: &[f64],
        dt_s: f64,
        liquid: CavitationLiquid,
        channel_viscosities_pa_s: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CompressibleCavitationContinuation, FilmError> {
        poll(&mut cancelled)?;
        self.inventory.validate_state(state)?;
        let inventory = self
            .inventory
            .prepare_liquid(liquid, channel_viscosities_pa_s)?;
        let film = Self::compile(inventory, self.bulk_modulus_pa)?;
        let step =
            film.advance_kernel(self, state, q_next, entrainment_m3_s, dt_s, &mut cancelled)?;
        Ok(CompressibleCavitationContinuation { film, step })
    }

    #[allow(
        clippy::large_stack_arrays,
        clippy::too_many_lines,
        clippy::too_many_arguments
    )]
    fn advance_kernel(
        &self,
        previous: &Self,
        state: &CavitatingFilmState,
        q_next: &[f64],
        entrainment_m3_s: &[f64],
        dt_s: f64,
        mut cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<CompressibleCavitationStep, FilmError> {
        poll(&mut cancelled)?;
        previous.inventory.validate_state(state)?;
        self.geometry().validate_configuration(q_next)?;
        if !dt_s.is_finite()
            || dt_s <= 0.0
            || entrainment_m3_s.len() != self.inventory.channel_count()
            || entrainment_m3_s.iter().any(|x| !x.is_finite())
        {
            return fail("invalid compressible cavitating-film duration or entrainment");
        }
        let n = self.inventory.cell_count();
        let mut step = Step {
            film: self,
            q: q_next,
            entrainment: entrainment_m3_s,
            dt: dt_s,
            volume: [0.0; MAX_CELLS],
            old_volume: [0.0; MAX_CELLS],
            old: [Point::default(); MAX_CELLS],
            previous_density: previous.density_at_cavitation_kg_m3,
            old_content: [[0.0; 2]; MAX_CELLS],
            conductance: [0.0; MAX_CHANNELS],
        };
        let mut candidate = Candidate::new();
        for (i, cell) in self.geometry().cells.iter().enumerate() {
            step.volume[i] = cell.area_m2 * (cell.gap.reference_m - dot(&cell.gap.closure, q_next));
            step.old_volume[i] = cell.area_m2
                * (cell.gap.reference_m - dot(&cell.gap.closure, state.configuration()));
            step.old[i] = previous.from_pressure(state.pressure_pa[i], state.liquid_fraction[i])?;
            candidate.points[i] = if self.liquid() == previous.liquid() {
                step.old[i]
            } else {
                // Warm-start in the new material's admitted pressure domain.
                // The original state below supplies mass, not new pressure.
                self.from_pressure(self.liquid().ambient_pressure_pa, state.liquid_fraction[i])?
            };
            if !step.volume[i].is_finite()
                || step.volume[i] <= 0.0
                || !step.old_volume[i].is_finite()
                || step.old_volume[i] <= 0.0
            {
                return fail("unrepresentable compressible cavitating-film volume");
            }
            for (part, value) in [step.old[i].theta, step.old[i].compression]
                .into_iter()
                .enumerate()
            {
                let content = step.old_volume[i] * value;
                let normalized = if step.previous_density == self.density_at_cavitation_kg_m3 {
                    content
                } else {
                    (step.previous_density * content) / self.density_at_cavitation_kg_m3
                };
                if !normalized.is_finite() || (content > 0.0 && normalized == 0.0) {
                    return fail("unrepresentable compressible-film material mass transfer");
                }
                step.old_content[i][part] = normalized;
            }
        }
        for (k, edge) in self.geometry().channels.iter().enumerate() {
            let h = edge.gap.reference_m - dot(&edge.gap.closure, q_next);
            let g = if h > 0.0 {
                self.geometry().factors[k] * h * h * h * self.inventory.pressure_law.fluidity_scale
            } else {
                0.0
            };
            if !g.is_finite() || (h > 0.0 && g <= 0.0) || (h <= 0.0 && entrainment_m3_s[k] != 0.0) {
                return fail("invalid compressible film conductance or closed-channel entrainment");
            }
            step.conductance[k] = g;
        }
        step.evaluate(&mut candidate, &mut cancelled)?;
        let mut trial = Candidate::new();
        for iteration in 0..=self.inventory.solve.maximum_iterations {
            poll(&mut cancelled)?;
            if candidate.admitted {
                return step.finish(&candidate, iteration, &mut cancelled);
            }
            if iteration == self.inventory.solve.maximum_iterations {
                break;
            }
            let mut direction = candidate.residual.map(|r| -r);
            solve_local(&mut candidate.jacobian, &mut direction, n, &mut cancelled)?;
            let mut changed = false;
            for i in 0..n {
                let p = &mut candidate.points[i];
                if p.theta == 1.0 && p.compression == 0.0 && direction[i] != 0.0 {
                    let full = direction[i] > 0.0;
                    changed |= full != p.full;
                    p.full = full;
                }
            }
            if changed {
                step.evaluate(&mut candidate, &mut cancelled)?;
                continue;
            }
            let mut fraction = 1.0_f64;
            for (i, &d) in direction[..n].iter().enumerate() {
                if d < 0.0 {
                    let p = candidate.points[i];
                    fraction = fraction.min(-(p.theta + p.compression) / d);
                }
            }
            if fraction < 1.0 {
                fraction *= 0.99;
            }
            let mut accepted = false;
            for _ in 0..32 {
                poll(&mut cancelled)?;
                let mut interrupted = false;
                let result = (|| {
                    for (i, &d) in direction[..n].iter().enumerate() {
                        let p = candidate.points[i];
                        let increment = fraction * d;
                        trial.points[i] = if p.full {
                            let x = p.compression + increment;
                            if x >= 0.0 {
                                self.point(1.0, x, true)?
                            } else {
                                self.point(1.0 + x, 0.0, false)?
                            }
                        } else {
                            let theta = p.theta + increment;
                            if theta > 1.0 {
                                self.point(1.0, theta - 1.0, true)?
                            } else {
                                self.point(theta, 0.0, false)?
                            }
                        };
                    }
                    step.evaluate(&mut trial, &mut || {
                        let stop = cancelled();
                        interrupted |= stop;
                        stop
                    })
                })();
                // A rejected constitutive trial can be shortened; cancellation
                // remains an immediate refusal, never a reason to try again.
                if interrupted {
                    return fail("cavitating-film solve cancelled");
                }
                if result.is_ok()
                    && (trial.admitted
                        || trial.maximum_residual
                            <= candidate.maximum_residual * (1.0 - 1e-4 * fraction))
                {
                    core::mem::swap(&mut candidate, &mut trial);
                    accepted = true;
                    break;
                }
                fraction *= 0.5;
            }
            if !accepted {
                return fail("compressible cavitating-film line-search budget exhausted");
            }
        }
        fail("compressible cavitating-film iteration budget exhausted")
    }
}

// Compensated log1p using the native deterministic logarithm. Keeping the
// authored compression x avoids losing p/B when 1+x rounds to one.
fn log1p(x: f64) -> f64 {
    let u = 1.0 + x;
    if u == 1.0 {
        x
    } else {
        det::ln(u) * (x / (u - 1.0))
    }
}

struct Candidate {
    points: [Point; MAX_CELLS],
    jacobian: [f64; MAX_CELLS * MAX_CELLS],
    residual: [f64; MAX_CELLS],
    maximum_residual: f64,
    admitted: bool,
}
impl Candidate {
    fn new() -> Self {
        Self {
            points: [Point::default(); MAX_CELLS],
            jacobian: [0.0; MAX_CELLS * MAX_CELLS],
            residual: [0.0; MAX_CELLS],
            maximum_residual: 0.0,
            admitted: false,
        }
    }
}

struct Step<'a> {
    film: &'a CompressibleCavitatingFilm,
    q: &'a [f64],
    entrainment: &'a [f64],
    dt: f64,
    volume: [f64; MAX_CELLS],
    old_volume: [f64; MAX_CELLS],
    old: [Point; MAX_CELLS],
    previous_density: f64,
    old_content: [[f64; 2]; MAX_CELLS],
    conductance: [f64; MAX_CHANNELS],
}

impl Step<'_> {
    // Integrate the positive mobility between the ACTUAL endpoints, avoiding
    // subtraction of nearly equal transformed potentials at high pressure.
    fn pressure_flow(&self, k: usize, a: Point, b: Point) -> Result<f64, FilmError> {
        let mu = self.film.geometry().channel_viscosities_pa_s()[k];
        for p in [a, b] {
            let viscosity = mu / p.fluidity;
            if !viscosity.is_finite() || viscosity <= 0.0 {
                return fail("unrepresentable compressible film pressure-dependent viscosity");
            }
        }
        let delta = a.compression - b.compression;
        if delta == 0.0 || self.conductance[k] == 0.0 {
            return Ok(0.0);
        }
        if self.film.inventory.pressure_law.alpha == 0.0 {
            return Ok(self.conductance[k] * self.film.bulk_modulus_pa * delta);
        }
        let low = if delta > 0.0 { b } else { a };
        let ratio = delta.abs() / (1.0 + low.compression);
        let logarithm = log1p(ratio);
        let z = self.film.exponent * logarithm;
        let log_ratio = if ratio == 0.0 { 1.0 } else { logarithm / ratio };
        let exp_ratio = if z == 0.0 { 1.0 } else { det::expm1(z) / z };
        Ok(self.conductance[k] * low.differential * delta * log_ratio * exp_ratio)
    }

    fn evaluate(
        &self,
        state: &mut Candidate,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<(), FilmError> {
        let n = self.film.inventory.cell_count();
        state.jacobian[..n * n].fill(0.0);
        state.maximum_residual = 0.0;
        state.admitted = true;
        let mut pressure = [FlowBalance::default(); MAX_CELLS];
        let mut advected = [FlowBalance::default(); MAX_CELLS];
        let mut pressure_scale = [0.0; MAX_CELLS];
        for (k, edge) in self.film.geometry().channels.iter().enumerate() {
            if k.is_multiple_of(32) {
                poll(cancelled)?;
            }
            let i = edge.from;
            let a = state.points[i];
            let b = edge.to.map_or(self.film.supplies[k], |j| state.points[j]);
            let amount = self.dt * self.pressure_flow(k, a, b)?;
            let e = self.dt * self.entrainment[k];
            let donor = if e >= 0.0 { a } else { b };
            let ga = if a.full {
                self.dt * self.conductance[k] * a.differential
            } else {
                0.0
            };
            let gb = if b.full {
                self.dt * self.conductance[k] * b.differential
            } else {
                0.0
            };
            pressure[i].add(amount);
            pressure_scale[i] += amount.abs();
            state.jacobian[i * n + i] += ga;
            if let Some(j) = edge.to {
                pressure[j].add(-amount);
                pressure_scale[j] += amount.abs();
                state.jacobian[j * n + i] -= ga;
                state.jacobian[i * n + j] -= gb;
                state.jacobian[j * n + j] += gb;
            }
            for part in [e * donor.theta, e * donor.compression] {
                advected[i].add(part);
                if let Some(j) = edge.to {
                    advected[j].add(-part);
                }
            }
            let donor_index = if e >= 0.0 { Some(i) } else { edge.to };
            if let Some(j) = donor_index {
                state.jacobian[i * n + j] += e;
                if let Some(to) = edge.to {
                    state.jacobian[to * n + j] -= e;
                }
            }
        }
        for i in 0..n {
            let p = state.points[i];
            let stored = self.volume[i] * (p.theta + p.compression);
            let old = if self.previous_density == self.film.density_at_cavitation_kg_m3 {
                self.old_volume[i] * (self.old[i].theta + self.old[i].compression)
            } else {
                self.old_content[i][0] + self.old_content[i][1]
            };
            let mut residual = pressure[i];
            residual.add(self.volume[i] * p.theta - self.old_content[i][0]);
            residual.add(self.volume[i] * p.compression - self.old_content[i][1]);
            residual.add(advected[i].value());
            state.residual[i] = residual.value();
            state.maximum_residual = state.maximum_residual.max(state.residual[i].abs());
            state.jacobian[i * n + i] += self.volume[i];
            let scale = stored.max(old) + advected[i].value().abs() + pressure_scale[i];
            let tolerance = self.film.inventory.solve.absolute_mass_tolerance_kg
                + self.film.inventory.solve.relative_mass_tolerance
                    * self.film.density_at_cavitation_kg_m3
                    * scale;
            let defect = self.film.density_at_cavitation_kg_m3 * state.residual[i];
            if !scale.is_finite() || !tolerance.is_finite() || !defect.is_finite() {
                return fail("unrepresentable compressible cavitating-film mass balance");
            }
            state.admitted &= defect.abs() <= tolerance;
        }
        if state.jacobian[..n * n].iter().any(|x| !x.is_finite()) {
            return fail("unrepresentable compressible cavitating-film Jacobian");
        }
        Ok(())
    }

    fn finish(
        &self,
        state: &Candidate,
        iterations: usize,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<CompressibleCavitationStep, FilmError> {
        let n = self.film.inventory.cell_count();
        let rho = self.film.density_at_cavitation_kg_m3;
        let mut cell_balance = [FlowBalance::default(); MAX_CELLS];
        let mut before = FlowBalance::default();
        let mut after = FlowBalance::default();
        let mut boundary = FlowBalance::default();
        let mut boundary_scale = 0.0;
        let mut forces = [0.0; MAX_PORTS];
        let mut flows = [0.0; MAX_CHANNELS];
        for i in 0..n {
            for (a, b) in [
                (state.points[i].theta, self.old[i].theta),
                (state.points[i].compression, self.old[i].compression),
            ] {
                let old = self.previous_density * (self.old_volume[i] * b);
                let new = rho * (self.volume[i] * a);
                before.add(old);
                after.add(new);
                cell_balance[i].add(new - old);
            }
            let gauge = state.points[i].pressure_excess
                - (self.film.liquid().ambient_pressure_pa
                    - self.film.liquid().cavitation_pressure_pa);
            let cell = &self.film.geometry().cells[i];
            for (j, &weight) in cell.gap.closure.iter().enumerate() {
                forces[j] += cell.area_m2 * gauge * weight;
            }
        }
        for (k, edge) in self.film.geometry().channels.iter().enumerate() {
            if k.is_multiple_of(32) {
                poll(cancelled)?;
            }
            let a = state.points[edge.from];
            let b = edge.to.map_or(self.film.supplies[k], |j| state.points[j]);
            let donor = if self.entrainment[k] >= 0.0 { a } else { b };
            let mut flow = FlowBalance::default();
            for component in [
                self.pressure_flow(k, a, b)?,
                self.entrainment[k] * donor.theta,
                self.entrainment[k] * donor.compression,
            ] {
                let mass_flow = rho * component;
                flow.add(mass_flow);
                let amount = self.dt * mass_flow;
                cell_balance[edge.from].add(amount);
                if let Some(j) = edge.to {
                    cell_balance[j].add(-amount);
                } else {
                    boundary.add(amount);
                    boundary_scale += amount.abs();
                }
            }
            flows[k] = flow.value();
        }
        let before = before.value();
        let after = after.value();
        let boundary = boundary.value();
        let mut total = FlowBalance::default();
        let mut maximum_cell = 0.0_f64;
        for (i, balance) in cell_balance[..n].iter().enumerate() {
            let residual = balance.value();
            let stored =
                rho * self.volume[i] * (state.points[i].theta + state.points[i].compression);
            let old = self.previous_density
                * self.old_volume[i]
                * (self.old[i].theta + self.old[i].compression);
            let allowance = self.film.inventory.solve.absolute_mass_tolerance_kg
                + self.film.inventory.solve.relative_mass_tolerance * (stored + old);
            if !residual.is_finite() || residual.abs() > allowance {
                return fail("compressible cavitating-film cell mass audit exceeds tolerance");
            }
            total.add(residual);
            maximum_cell = maximum_cell.max(residual.abs());
        }
        let tolerance = self.film.inventory.solve.absolute_mass_tolerance_kg
            + self.film.inventory.solve.relative_mass_tolerance * (before + after + boundary_scale);
        let total_residual = after - before + boundary;
        if !tolerance.is_finite()
            || !total_residual.is_finite()
            || total_residual.abs() > tolerance
            || total.value().abs() > tolerance
            || forces[..self.film.geometry().port_count()]
                .iter()
                .any(|x| !x.is_finite())
            || flows[..self.film.inventory.channel_count()]
                .iter()
                .any(|x| !x.is_finite())
        {
            return fail("compressible cavitating-film global mass or force audit failed");
        }
        poll(cancelled)?;
        let candidate = CavitatingFilmState {
            owner: Arc::clone(&self.film.inventory.owner),
            configuration: self.q.to_vec(),
            liquid_fraction: state.points[..n].iter().map(|p| p.theta).collect(),
            pressure_pa: state.points[..n]
                .iter()
                .map(|p| p.pressure_excess + self.film.liquid().cavitation_pressure_pa)
                .collect(),
        };
        self.film.inventory.validate_state(&candidate)?;
        let report = CavitationReport {
            liquid_mass_before_kg: before,
            liquid_mass_after_kg: after,
            boundary_mass_out_kg: boundary,
            mass_residual_kg: total_residual,
            maximum_cell_mass_residual_kg: maximum_cell,
            minimum_liquid_fraction: candidate
                .liquid_fraction
                .iter()
                .copied()
                .fold(1.0, f64::min),
            maximum_pressure_pa: candidate.pressure_pa.iter().copied().fold(0.0, f64::max),
            iterations,
        };
        Ok(CompressibleCavitationStep {
            state: candidate,
            channel_mass_flow_kg_s: flows[..self.film.inventory.channel_count()].to_vec(),
            generalized_force: forces[..self.film.geometry().port_count()].to_vec(),
            report,
        })
    }
}
