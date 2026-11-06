//! Mass-conserving Elrod--Adams transport of a generalized-Newtonian liquid.
//!
//! The original cavitation owner retains volume, inventory, absolute-pressure
//! supplies, complementarity and the independent mass audit. Each face uses
//! `Q = Qp(delta_p, Uupper-Ulower) + Qentrained*theta_upwind`; Qp and its
//! positive pressure differential come from the existing source-owned nonlinear
//! Couette--Poiseuille profile. There is no separately superposed shear law.
//!
//! One unknown per cell represents either positive excess pressure or missing
//! liquid content. Its selected pressure/content Jacobian uses the native
//! bounded dense elimination and upwind storage operator. Damped corrections
//! preserve pressure/content complementarity and the original validity domain.
//! See the conservative p-theta formulation in Profito et al. (2015),
//! https://doi.org/10.1007/s11249-015-0588-0, sections 2.2--2.3; extension to
//! generalized-Newtonian lubricants requires the actual constitutive face flux.
//!
//! Only endpoint NORMAL pressure force and liquid inventory are returned.
//! Cavitated tangential traction, viscous heat partition, latent heat, thermal
//! evolution and compressible trapped liquid are not supplied by this model.
//! The private borrowed transport below also serves original Carreau passages;
//! both families retain this same pressure/content iteration and mass audit.

use super::{
    CavitatingFilm, CavitatingFilmState, CavitationLiquid, CavitationReservoir,
    CavitationSolveConfig, CavitationStep, FilmError, FlowBalance, MAX_CELLS, MAX_CHANNELS,
    ResistiveFilm, dot, fail, poll, solve_local,
};
use crate::resistive_film::carreau::CarreauFilm;
use crate::resistive_film::power_law::{PowerLawControls, PowerLawFilm, RegularizedPowerLaw};

/// The original liquid-content owner with explicit SI power-law rheology.
///
/// Density, rupture pressure and geometry remain fixed during each step.
/// Relative wall speed and full-film mean entrainment are separately authored
/// in original passage orientation. The original active-set budget bounds
/// pressure/content corrections and corner phase pivots; the constitutive
/// owner's iteration and backtracking ceilings bound local stress work and
/// candidate trials.
/// Mass admission uses `CavitationSolveConfig`, not a full-film power receipt.
#[derive(Clone, Debug)]
pub struct CavitatingPowerLawFilm {
    inventory: CavitatingFilm,
    rheology: PowerLawFilm,
}

/// Endpoint constitutive owner and the matching conservative liquid candidate.
/// Accept both together after the surrounding thermal/mechanical step succeeds.
#[derive(Clone, Debug)]
pub struct CavitatingPowerLawContinuation {
    pub film: CavitatingPowerLawFilm,
    pub step: CavitationStep,
}

impl CavitatingFilm {
    /// Bind explicit nonlinear laws to the SAME accepted liquid inventory.
    ///
    /// Existing states, geometry, supplies and mass tolerances remain valid.
    /// Only the transport law changes; source validity and thermal/rheological
    /// history remain with the caller. A native Barus law is not compatible
    /// with this pressure-independent generalized-Newtonian reduction.
    pub fn with_power_law(
        &self,
        laws: Vec<RegularizedPowerLaw>,
        controls: PowerLawControls,
    ) -> Result<CavitatingPowerLawFilm, FilmError> {
        let rheology = PowerLawFilm::new(self.geometry(), laws, controls)?;
        let inventory =
            self.with_channel_viscosities(rheology.pressure_operator().channel_viscosities_pa_s())?;
        Ok(CavitatingPowerLawFilm {
            inventory,
            rheology,
        })
    }
}

impl CavitatingPowerLawFilm {
    /// Bind the existing nonlinear liquid law and native conservative inventory.
    /// The supplied law already admits positive K, n and its physical core.
    /// Pressure-dependent viscosity is excluded by that original law owner.
    pub fn new(
        rheology: PowerLawFilm,
        liquid: CavitationLiquid,
        solve: CavitationSolveConfig,
    ) -> Result<Self, FilmError> {
        let inventory = CavitatingFilm::new(rheology.pressure_operator().clone(), liquid, solve)?;
        Ok(Self {
            inventory,
            rheology,
        })
    }

    /// Original geometry. Its core-viscosity pressure evaluator is not a
    /// substitute for this nonlinear, mass-conserving film step.
    #[must_use]
    pub fn geometry(&self) -> &ResistiveFilm {
        self.inventory.geometry()
    }

    /// Explicit homogeneous liquid density, absolute rupture and reference data.
    #[must_use]
    pub fn liquid(&self) -> CavitationLiquid {
        self.inventory.liquid()
    }

    /// Original channel-ordered generalized-Newtonian material laws.
    #[must_use]
    pub fn laws(&self) -> &[RegularizedPowerLaw] {
        self.rheology.laws()
    }

    /// Prescribed exterior supplies; internal passages carry `None`.
    pub fn reservoir_conditions(
        &self,
    ) -> impl ExactSizeIterator<Item = Option<CavitationReservoir>> + '_ {
        self.inventory.reservoir_conditions()
    }

    /// Rebind original absolute-pressure and incoming-content supplies without
    /// resetting existing liquid content, pressure or state identity.
    pub fn with_reservoirs(
        &self,
        reservoirs: &[Option<CavitationReservoir>],
    ) -> Result<Self, FilmError> {
        Ok(Self {
            inventory: self.inventory.with_reservoirs(reservoirs)?,
            rheology: self.rheology.clone(),
        })
    }

    /// Refresh viscosity laws while continuing exactly the accepted inventory.
    /// Density/rupture pressure do not change; this is not a thermal time step.
    pub fn with_laws(&self, laws: Vec<RegularizedPowerLaw>) -> Result<Self, FilmError> {
        if laws == self.laws() {
            return Ok(self.clone());
        }
        let rheology = PowerLawFilm::new(self.geometry(), laws, self.rheology.controls())?;
        let inventory = self
            .inventory
            .with_channel_viscosities(rheology.pressure_operator().channel_viscosities_pa_s())?;
        Ok(Self {
            inventory,
            rheology,
        })
    }

    /// Native full-liquid initial condition at the original reference pressure.
    pub fn initial_state(&self, q: &[f64]) -> Result<CavitatingFilmState, FilmError> {
        self.inventory.initial_state(q)
    }

    /// Explicit starvation; partial cells start at the absolute rupture pressure.
    pub fn initialize(
        &self,
        q: &[f64],
        liquid_fraction: &[f64],
    ) -> Result<CavitatingFilmState, FilmError> {
        self.inventory.initialize(q, liquid_fraction)
    }

    /// Advance original liquid mass with a simultaneous nonlinear face profile.
    ///
    /// `entrainment_m3_s` is the signed FULL-film mean-wall volume transport,
    /// multiplied here by the upwind liquid content. `relative_wall_speed_m_s`
    /// is upper minus lower wall speed along the same passage. For homogeneous
    /// no-slip walls the former is width*gap*(Ulower+Uupper)/2. The caller owns
    /// geometry/wall-time quadrature and the thin-gap creeping-flow regime.
    ///
    /// All n=1 laws dispatch directly to the existing native cavitation solve.
    /// Failed mass/pressure admission, solver budgets or cancellation publish
    /// no candidate and never change the accepted state or rheology.
    #[allow(clippy::too_many_arguments)]
    pub fn advance(
        &self,
        state: &CavitatingFilmState,
        q_next: &[f64],
        entrainment_m3_s: &[f64],
        relative_wall_speed_m_s: &[f64],
        dt_s: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CavitationStep, FilmError> {
        poll(&mut cancelled)?;
        self.validate_relative_motion(q_next, relative_wall_speed_m_s)?;
        self.advance_kernel(
            &self.inventory,
            state,
            q_next,
            entrainment_m3_s,
            relative_wall_speed_m_s,
            dt_s,
            &mut cancelled,
        )
    }

    /// Continue the ORIGINAL mass using explicit endpoint liquid and rheology.
    ///
    /// The native balance remains
    /// `rho_next*(V_next*theta_next + dt*div Q_next) = rho_old*V_old*theta_old`.
    /// Endpoint channel laws drive the simultaneous pressure/sliding profiles;
    /// a new density never reprices the old liquid volume. Pressure/content is
    /// solved at the endpoint rupture pressure. Explicit reservoirs retain
    /// their original pressure/content and must satisfy that new domain.
    ///
    /// This is a first-order endpoint constitutive split. The caller supplies
    /// the source/temperature proposal and temporal refinement; no caloric
    /// energy, phase heat or rheological evolution is inferred. A changed
    /// liquid has a distinct inventory identity, while a rheology-only refresh
    /// keeps the original state compatibility. If every ENDPOINT law is n=1,
    /// the existing Newtonian mass-transfer kernel is used exactly.
    ///
    /// Any material, domain, budget or cancellation refusal leaves the accepted
    /// state and both original owners intact. The returned film and candidate
    /// are a pair, to be accepted together by the surrounding time owner.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_liquid(
        &self,
        state: &CavitatingFilmState,
        q_next: &[f64],
        entrainment_m3_s: &[f64],
        relative_wall_speed_m_s: &[f64],
        dt_s: f64,
        liquid: CavitationLiquid,
        laws: Vec<RegularizedPowerLaw>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CavitatingPowerLawContinuation, FilmError> {
        poll(&mut cancelled)?;
        self.inventory.validate_state(state)?;
        self.validate_relative_motion(q_next, relative_wall_speed_m_s)?;
        let rheology = PowerLawFilm::new(self.geometry(), laws, self.rheology.controls())?;
        let inventory = self.inventory.prepare_liquid(
            liquid,
            rheology.pressure_operator().channel_viscosities_pa_s(),
        )?;
        let film = Self {
            inventory,
            rheology,
        };
        let step = film.advance_kernel(
            &self.inventory,
            state,
            q_next,
            entrainment_m3_s,
            relative_wall_speed_m_s,
            dt_s,
            &mut cancelled,
        )?;
        Ok(CavitatingPowerLawContinuation { film, step })
    }

    fn validate_relative_motion(
        &self,
        q_next: &[f64],
        relative_wall_speed_m_s: &[f64],
    ) -> Result<(), FilmError> {
        if relative_wall_speed_m_s.len() != self.inventory.channel_count()
            || relative_wall_speed_m_s.iter().any(|v| !v.is_finite())
        {
            return fail("invalid cavitating power-law relative wall speeds");
        }
        self.geometry().validate_configuration(q_next)?;
        for (edge, &relative) in self.geometry().channels.iter().zip(relative_wall_speed_m_s) {
            if edge.gap.reference_m - dot(&edge.gap.closure, q_next) <= 0.0 && relative != 0.0 {
                return fail("relative sliding through a closed cavitating passage");
            }
        }
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    fn advance_kernel(
        &self,
        previous: &CavitatingFilm,
        state: &CavitatingFilmState,
        q_next: &[f64],
        entrainment_m3_s: &[f64],
        relative_wall_speed_m_s: &[f64],
        dt_s: f64,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<CavitationStep, FilmError> {
        let transport = if self.rheology.is_newtonian() {
            None
        } else {
            Some(Transport::PowerLaw {
                rheology: &self.rheology,
                relative: relative_wall_speed_m_s,
            })
        };
        self.inventory.advance_kernel(
            previous,
            state,
            q_next,
            entrainment_m3_s,
            dt_s,
            transport,
            cancelled,
        )
    }
}

// Both original material owners use the same conservative pressure/content
// iteration. The bounded borrowed choice adds no per-face heap allocation or
// independent constitutive approximation.
pub(super) enum Transport<'a> {
    PowerLaw {
        rheology: &'a PowerLawFilm,
        relative: &'a [f64],
    },
    Carreau {
        rheology: &'a CarreauFilm,
        relative: &'a [f64],
    },
}

impl Transport<'_> {
    fn pressure_response(
        &self,
        channel: usize,
        gap: f64,
        core_conductance: f64,
        dp: f64,
        check: &mut impl FnMut() -> Result<(), FilmError>,
    ) -> Result<(f64, f64), FilmError> {
        match self {
            Self::PowerLaw { rheology, relative } => rheology.passage_pressure_response(
                channel,
                gap,
                core_conductance,
                dp,
                relative[channel],
                check,
            ),
            Self::Carreau { rheology, relative } => rheology.passage_pressure_response(
                channel,
                gap,
                core_conductance,
                dp,
                relative[channel],
                check,
            ),
        }
    }

    fn maximum_line_search_steps(&self) -> usize {
        match self {
            Self::PowerLaw { rheology, .. } => rheology.controls().maximum_line_search_steps,
            Self::Carreau { rheology, .. } => {
                rheology.controls().pressure.maximum_line_search_steps
            }
        }
    }
}

struct Candidate {
    pressure: [f64; MAX_CELLS],
    theta: [f64; MAX_CELLS],
    full: [bool; MAX_CELLS],
    jacobian: [f64; MAX_CELLS * MAX_CELLS],
    residual: [f64; MAX_CELLS],
    maximum_residual: f64,
    admitted: bool,
}

impl Candidate {
    fn new() -> Self {
        Self {
            pressure: [0.0; MAX_CELLS],
            theta: [0.0; MAX_CELLS],
            full: [false; MAX_CELLS],
            jacobian: [0.0; MAX_CELLS * MAX_CELLS],
            residual: [0.0; MAX_CELLS],
            maximum_residual: 0.0,
            admitted: false,
        }
    }
}

// All fields are borrowed from the SAME native step assembly and original
// state. This constitutive leaf does not rebuild topology, source reservoirs,
// inventory identity, density epochs or force projection.
struct Step<'a> {
    film: &'a CavitatingFilm,
    transport: Transport<'a>,
    entrainment: &'a [f64],
    dt: f64,
    gaps: [f64; MAX_CHANNELS],
    volume: &'a [f64; MAX_CELLS],
    old_volume: &'a [f64; MAX_CELLS],
    core_conductance: &'a [f64; MAX_CHANNELS],
    previous_density: f64,
}

impl Step<'_> {
    fn original_content(&self, cell: usize) -> f64 {
        if self.previous_density == self.film.liquid.density_kg_m3 {
            self.old_volume[cell]
        } else {
            (self.previous_density * self.old_volume[cell]) / self.film.liquid.density_kg_m3
        }
    }

    fn pressure_flow(
        &self,
        channel: usize,
        dp: f64,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<(f64, f64), FilmError> {
        self.transport.pressure_response(
            channel,
            self.gaps[channel],
            self.core_conductance[channel],
            dp,
            &mut || poll(cancelled),
        )
    }

    fn evaluate(
        &self,
        state: &mut Candidate,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<(), FilmError> {
        let n = self.film.cell_count();
        state.jacobian[..n * n].fill(0.0);
        state.maximum_residual = 0.0;
        state.admitted = true;
        let mut transported = [FlowBalance::default(); MAX_CELLS];
        let mut pressure_transfer = [FlowBalance::default(); MAX_CELLS];
        let mut pressure_scale = [0.0; MAX_CELLS];
        for (k, edge) in self.film.geometry.channels.iter().enumerate() {
            if k.is_multiple_of(32) {
                poll(cancelled)?;
            }
            let dp = state.pressure[edge.from]
                - edge
                    .to
                    .map_or(self.film.boundary_potential(k), |j| state.pressure[j]);
            let (flow, differential) = self.pressure_flow(k, dp, cancelled)?;
            let pressure_amount = self.dt * flow;
            let g = self.dt * differential;
            let donor = if self.entrainment[k] >= 0.0 {
                state.theta[edge.from]
            } else {
                edge.to
                    .map_or(self.film.boundary_liquid_fraction(k), |j| state.theta[j])
            };
            let entrained_amount = (self.dt * self.entrainment[k]) * donor;
            if [pressure_amount, g, entrained_amount]
                .iter()
                .any(|x| !x.is_finite())
                || g < 0.0
            {
                return fail("unrepresentable nonlinear cavitating face response");
            }
            let i = edge.from;
            pressure_transfer[i].add(pressure_amount);
            transported[i].add(entrained_amount);
            pressure_scale[i] += pressure_amount.abs();
            state.jacobian[i * n + i] += g;
            if let Some(j) = edge.to {
                pressure_transfer[j].add(-pressure_amount);
                transported[j].add(-entrained_amount);
                pressure_scale[j] += pressure_amount.abs();
                state.jacobian[j * n + j] += g;
                state.jacobian[i * n + j] -= g;
                state.jacobian[j * n + i] -= g;
            }
        }
        let rho = self.film.liquid.density_kg_m3;
        for i in 0..n {
            let stored = self.volume[i] * state.theta[i];
            let old = self.original_content(i);
            let advection = transported[i].value();
            let mut residual = pressure_transfer[i];
            residual.add(stored - old);
            residual.add(advection);
            state.residual[i] = residual.value();
            state.maximum_residual = state.maximum_residual.max(state.residual[i].abs());
            // Gross divergence-free wall transport cannot hide a small
            // pressure defect. Only resolved inventory/net transport and
            // pressure-flow scales enter this nonlinear acceptance check.
            let scale = stored.abs().max(old.abs()) + advection.abs() + pressure_scale[i];
            let tolerance = self.film.solve.absolute_mass_tolerance_kg
                + self.film.solve.relative_mass_tolerance * rho * scale;
            let defect = rho * state.residual[i];
            if !scale.is_finite() || !tolerance.is_finite() || !defect.is_finite() {
                return fail("unrepresentable nonlinear cavitating mass residual");
            }
            state.admitted &= defect.abs() <= tolerance;
        }
        Ok(())
    }
}

#[allow(clippy::too_many_arguments, clippy::large_stack_arrays)]
pub(super) fn advance(
    film: &CavitatingFilm,
    previous: &CavitatingFilmState,
    q: &[f64],
    entrainment: &[f64],
    dt: f64,
    previous_density: f64,
    volume: &[f64; MAX_CELLS],
    old_volume: &[f64; MAX_CELLS],
    core_conductance: &[f64; MAX_CHANNELS],
    storage: &[f64; MAX_CELLS * MAX_CELLS],
    initial_full: &[bool; MAX_CELLS],
    transport: Transport<'_>,
    cancelled: &mut dyn FnMut() -> bool,
) -> Result<CavitationStep, FilmError> {
    let n = film.cell_count();
    let mut step = Step {
        film,
        transport,
        entrainment,
        dt,
        gaps: [0.0; MAX_CHANNELS],
        volume,
        old_volume,
        core_conductance,
        previous_density,
    };
    for (k, edge) in film.geometry.channels.iter().enumerate() {
        step.gaps[k] = edge.gap.reference_m - dot(&edge.gap.closure, q);
    }
    let maximum_pressure = (film.liquid.ambient_pressure_pa - film.liquid.cavitation_pressure_pa)
        + film.geometry.limits.maximum_pressure_pa;
    if !maximum_pressure.is_finite() || maximum_pressure <= 0.0 {
        return fail("unrepresentable nonlinear cavitating pressure domain");
    }
    let mut state = Candidate::new();
    state.full = *initial_full;
    for i in 0..n {
        if state.full[i] {
            state.pressure[i] =
                (previous.pressure_pa[i] - film.liquid.cavitation_pressure_pa).max(0.0);
            state.theta[i] = 1.0;
        } else {
            state.theta[i] = (step.original_content(i) / volume[i]).min(1.0);
        }
    }
    step.evaluate(&mut state, cancelled)?;
    let mut trial = Candidate::new();
    for iteration in 0..=film.solve.maximum_iterations {
        poll(cancelled)?;
        if iteration == film.solve.maximum_iterations && !state.admitted {
            return fail("nonlinear cavitating pressure/content iteration budget exhausted");
        }
        let mut matrix = [0.0; MAX_CELLS * MAX_CELLS];
        let mut direction = state.residual.map(|r| -r);
        let mut scale = [0.0; MAX_CELLS];
        let mut coordinate = [0.0; MAX_CELLS];
        for j in 0..n {
            let a = state.jacobian[j * n + j];
            let b = storage[j * n + j];
            if a > 0.0 {
                scale[j] = b / a;
                if !scale[j].is_finite() || scale[j] <= 0.0 {
                    return fail("unrepresentable nonlinear cavitating pressure scale");
                }
            } else if state.full[j] {
                return fail("fully filled undrained film requires liquid compressibility");
            }
            coordinate[j] = if state.full[j] {
                state.pressure[j] / scale[j]
            } else {
                state.theta[j] - 1.0
            };
            for i in 0..n {
                matrix[i * n + j] = if state.full[j] {
                    (state.jacobian[i * n + j] / a) * b
                } else {
                    storage[i * n + j]
                };
            }
        }
        // Connectivity/admission still applies at zero residual; no pressure
        // datum or compressibility is invented for a sealed saturated island.
        solve_local(&mut matrix, &mut direction, n, cancelled)?;
        if state.admitted {
            let mut pressure = state.pressure;
            return film.finish::<true>(
                q,
                entrainment,
                dt,
                iteration.max(1),
                previous_density,
                volume,
                old_volume,
                |k, dp, cancelled| {
                    step.pressure_flow(k, dp, cancelled)
                        .map(|response| response.0)
                },
                &state.theta,
                &mut pressure,
                cancelled,
            );
        }
        // At a full/rupture corner the opposite phase has a different
        // Jacobian. A content direction mapped into pressure need not descend
        // even for a linear constitutive law on a branching graph. Select its
        // entered phase and relinearize at the SAME physical state first.
        // These pivots consume the original bounded active-set budget.
        let mut phase_changed = false;
        for i in 0..n {
            if state.pressure[i] == 0.0 && state.theta[i] == 1.0 && direction[i] != 0.0 {
                let full = direction[i] > 0.0;
                phase_changed |= state.full[i] != full;
                state.full[i] = full;
            }
        }
        if phase_changed {
            continue;
        }
        let mut fraction = 1.0_f64;
        for i in 0..n {
            let bound = if direction[i] > 0.0 {
                let top = if scale[i] > 0.0 {
                    maximum_pressure / scale[i]
                } else {
                    0.0
                };
                (top - coordinate[i]) / direction[i]
            } else if direction[i] < 0.0 {
                (-1.0 - coordinate[i]) / direction[i]
            } else {
                1.0
            };
            fraction = fraction.min(bound);
        }
        if fraction < 1.0 {
            fraction *= 0.99;
        }
        if !fraction.is_finite() || fraction <= 0.0 {
            return fail("nonlinear cavitating correction cannot remain in its physical domain");
        }
        let mut accepted = false;
        for _ in 0..step.transport.maximum_line_search_steps() {
            poll(cancelled)?;
            for i in 0..n {
                let next = coordinate[i] + fraction * direction[i];
                trial.full[i] = if next == 0.0 {
                    state.full[i]
                } else {
                    next > 0.0
                };
                if next >= 0.0 {
                    trial.pressure[i] = next * scale[i];
                    trial.theta[i] = 1.0;
                } else {
                    trial.pressure[i] = 0.0;
                    trial.theta[i] = 1.0 + next;
                }
            }
            let mut callback_refused = false;
            let result = step.evaluate(&mut trial, &mut || {
                let refused = cancelled();
                callback_refused |= refused;
                refused
            });
            if callback_refused {
                return result.and_then(|()| fail("cavitating-film step cancelled"));
            }
            poll(cancelled)?;
            if result.is_ok()
                && (trial.admitted
                    || trial.maximum_residual <= state.maximum_residual * (1.0 - 1e-4 * fraction))
            {
                core::mem::swap(&mut state, &mut trial);
                accepted = true;
                break;
            }
            fraction *= 0.5;
        }
        if !accepted {
            return fail("nonlinear cavitating line-search budget exhausted");
        }
    }
    fail("nonlinear cavitating pressure/content iteration budget exhausted")
}
