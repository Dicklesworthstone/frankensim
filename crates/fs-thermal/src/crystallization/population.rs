//! Common-temperature material containing distinct original kinetic histories.
//!
//! Instant thermal equilibration is an explicit closure: incoming masses share
//! one temperature without advancing their kinetic clocks. Each mass keeps its
//! own Avrami exposure. Averaging exposure, or replacing the population with an
//! equivalent age inferred from its current crystallinity, generally changes
//! future conversion when n != 1 and is deliberately unavailable here.
//!
//! The admitted cohorts have the same material, caloric reference, pressure,
//! composition and kinetic law. Their original storage is linear in mass
//! fraction, so h(T) = h_a(T) - mean(X)*(h_a(T)-h_c(T)). At shared temperature
//! every nonsaturated cohort receives the same integral of K(T). This scalar
//! increment is the only evolving history coordinate; immutable cohort data
//! are shared across trial states. Each fraction query is O(cohorts), while
//! time steps allocate no cohort buffers. No histories are binned or truncated.
//!
//! This is a parallel-cohort constitutive model, not a spatial transport or
//! molecular mixing solver. It supplies no inter-cohort nucleation/impingement,
//! melting, pressure change, mixing enthalpy, mechanical work or flow geometry.
//! Input masses must come from the caller's conservative transport owner.
//! The explicit midpoint heat-path step is second order on smooth trajectories;
//! its timestep and convergence are the transport owner's responsibility.
//! Optional bounded step doubling resolves every retained cohort's conversion
//! and shared age increment on that same prescribed enthalpy path.

use std::sync::Arc;

use super::{
    CrystallizationError, CrystallizationLaw, CrystallizationSolve, CrystallizationState,
    check_cancelled, check_duration, check_solve, finite, positive,
};

const MAXIMUM_COHORTS: usize = 1_000_000;
const CHECKPOINT_STRIDE: usize = 32;

#[path = "population/adaptive.rs"]
mod adaptive;
pub use adaptive::AdaptiveCrystallizationPopulationStep;

#[derive(Clone, Copy)]
struct Cohort {
    mass_kg: f64,
    initial_exposure: f64,
    initially_saturated: bool,
}

/// Accepted thermally equilibrated masses with distinct, retained crystal ages.
#[derive(Clone)]
pub struct CrystallizationPopulation<'a> {
    law: CrystallizationLaw<'a>,
    cohorts: Arc<Vec<Cohort>>,
    mass_kg: f64,
    specific_enthalpy_j_kg: f64,
    temperature_k: f64,
    crystalline_mass_fraction: f64,
    exposure_increment: f64,
    fully_crystallized: bool,
    caloric: CrystallizationSolve,
}

/// One accepted population interval; latent heat remains in original enthalpy.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CrystallizationPopulationStep {
    pub duration_s: f64,
    pub specific_heat_j_kg: f64,
    pub crystallized_mass_kg: f64,
    pub caloric_residual_j_kg: f64,
    pub midpoint_rate_s_inv: f64,
}

impl<'a> CrystallizationPopulation<'a> {
    /// Equilibrate incoming masses at unchanged phase ages and total enthalpy.
    ///
    /// The complete incoming state is validated against the same original law.
    /// Cohorts remain in input order for deterministic reductions. The supplied
    /// count budget is checked before allocating, and an empty/zero-mass cohort
    /// is refused rather than removed. A mixture requiring melting refuses.
    pub fn new(
        law: CrystallizationLaw<'a>,
        incoming: &[(f64, CrystallizationState)],
        maximum_cohorts: usize,
        solve: CrystallizationSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, CrystallizationError> {
        check_cancelled(&mut cancelled)?;
        check_solve(solve)?;
        if maximum_cohorts == 0
            || maximum_cohorts > MAXIMUM_COHORTS
            || incoming.is_empty()
            || incoming.len() > maximum_cohorts
        {
            return Err(invalid(
                "crystallization population exceeds its cohort budget",
            ));
        }
        let mut cohorts = Vec::new();
        cohorts
            .try_reserve_exact(incoming.len())
            .map_err(|_| invalid("crystallization cohort allocation refused"))?;
        let mut mass_kg = 0.0;
        let mut enthalpy = 0.0;
        let mut fraction = 0.0;
        let mut fully_crystallized = true;
        for (index, &(mass, state)) in incoming.iter().enumerate() {
            if index % CHECKPOINT_STRIDE == 0 {
                check_cancelled(&mut cancelled)?;
            }
            if !positive(mass) {
                return Err(invalid(
                    "crystallization cohorts require positive finite mass",
                ));
            }
            law.validate_state(state, solve)?;
            let next_mass = add_mass(mass_kg, mass)?;
            let weight = mass / next_mass;
            enthalpy = mean(enthalpy, state.specific_enthalpy_j_kg(), weight)?;
            fraction = mean(fraction, state.crystalline_mass_fraction(), weight)?;
            let saturated =
                state.crystalline_mass_fraction() == law.config().maximum_crystalline_mass_fraction;
            fully_crystallized &= saturated;
            cohorts.push(Cohort {
                mass_kg: mass,
                initial_exposure: state.exposure(),
                initially_saturated: saturated,
            });
            mass_kg = next_mass;
        }
        let temperature_k =
            law.temperature_at_enthalpy_fraction(enthalpy, fraction, solve, &mut cancelled)?;
        check_cancelled(&mut cancelled)?;
        Ok(Self {
            law,
            cohorts: Arc::new(cohorts),
            mass_kg,
            specific_enthalpy_j_kg: enthalpy,
            temperature_k,
            crystalline_mass_fraction: fraction,
            exposure_increment: 0.0,
            fully_crystallized,
            caloric: solve,
        })
    }

    pub const fn law(&self) -> &CrystallizationLaw<'a> {
        &self.law
    }
    pub const fn mass_kg(&self) -> f64 {
        self.mass_kg
    }
    pub const fn specific_enthalpy_j_kg(&self) -> f64 {
        self.specific_enthalpy_j_kg
    }
    pub const fn temperature_k(&self) -> f64 {
        self.temperature_k
    }
    pub const fn crystalline_mass_fraction(&self) -> f64 {
        self.crystalline_mass_fraction
    }
    pub fn cohort_count(&self) -> usize {
        self.cohorts.len()
    }
    pub fn cohort_mass_kg(&self, index: usize) -> Option<f64> {
        self.cohorts.get(index).map(|cohort| cohort.mass_kg)
    }

    /// Recover an individual original-law state for a later conservative
    /// transfer. Its enthalpy is its own phase mixture at the common temperature;
    /// no fictitious single-age state represents the whole population.
    pub fn cohort_state(
        &self,
        index: usize,
    ) -> Result<Option<CrystallizationState>, CrystallizationError> {
        self.cohorts
            .get(index)
            .map(|cohort| {
                self.law.state_at_temperature(
                    self.temperature_k,
                    exposure(cohort, self.exposure_increment)?,
                )
            })
            .transpose()
    }

    /// Total H on the original source datum. Use step heat for small increments
    /// that a large arbitrary formation datum would hide in total subtraction.
    pub fn enthalpy_j(&self) -> Result<f64, CrystallizationError> {
        finite(self.mass_kg * self.specific_enthalpy_j_kg)
    }

    /// Recheck original calorics and every retained kinetic history without
    /// changing any state. Cross-module consumers can use the same local budget.
    pub fn validate(
        &self,
        solve: CrystallizationSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), CrystallizationError> {
        check_solve(solve)?;
        check_cancelled(&mut cancelled)?;
        let (fraction, saturated) = self.fraction_at(self.exposure_increment, &mut cancelled)?;
        if fraction.to_bits() != self.crystalline_mass_fraction.to_bits()
            || saturated != self.fully_crystallized
        {
            return Err(CrystallizationError::StateMismatch);
        }
        let reconstructed = self
            .law
            .specific_enthalpy_j_kg(self.temperature_k, fraction)?;
        if (reconstructed - self.specific_enthalpy_j_kg).abs()
            > solve.specific_enthalpy_tolerance_j_kg
        {
            return Err(CrystallizationError::StateMismatch);
        }
        check_cancelled(&mut cancelled)
    }

    /// Exact isothermal advancement of all original cohort exposures. Negative
    /// specific heat is the energy removed by the common-temperature thermostat.
    pub fn advance_isothermal(
        &self,
        duration_s: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(Self, CrystallizationPopulationStep), CrystallizationError> {
        check_duration(duration_s)?;
        self.validate(self.caloric, &mut cancelled)?;
        let rate = self.law.rate_s_inv(self.temperature_k)?;
        let increment = self.next_exposure(rate, duration_s)?;
        let (fraction, saturated) = self.fraction_at(increment, &mut cancelled)?;
        let h = self
            .law
            .specific_enthalpy_j_kg(self.temperature_k, fraction)?;
        let mut next = self.clone();
        next.exposure_increment = increment;
        next.crystalline_mass_fraction = fraction;
        next.fully_crystallized = saturated;
        next.specific_enthalpy_j_kg = h;
        let step = self.report(&next, duration_s, rate, 0.0)?;
        check_cancelled(&mut cancelled)?;
        Ok((next, step))
    }

    /// Advance the common thermal state along a prescribed linear enthalpy path.
    /// Both midpoint and endpoint use the original mixture inversion with the
    /// actual mass-weighted cohort fractions. Latent release therefore changes
    /// the shared temperature and subsequent rate, without an added heat source.
    /// A refused stage requires a smaller step or a different admitted model;
    /// it never clips age, crystallinity, temperature or energy.
    pub fn advance_specific_enthalpy(
        &self,
        final_h_j_kg: f64,
        duration_s: f64,
        solve: CrystallizationSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(Self, CrystallizationPopulationStep), CrystallizationError> {
        check_duration(duration_s)?;
        self.validate(solve, &mut cancelled)?;
        let heat = finite(final_h_j_kg - self.specific_enthalpy_j_kg)?;
        let rate = self.law.rate_s_inv(self.temperature_k)?;
        if heat == 0.0 && (rate == 0.0 || self.fully_crystallized) {
            // No physical derivative is active. Re-inverting an already
            // accepted caloric residual can introduce artificial undercooling
            // at Tm and restart growth during a stationary hold.
            let residual = finite(
                self.law
                    .specific_enthalpy_j_kg(self.temperature_k, self.crystalline_mass_fraction)?
                    - self.specific_enthalpy_j_kg,
            )?;
            let mut next = self.clone();
            next.caloric = solve;
            let step = self.report(&next, duration_s, rate, residual)?;
            check_cancelled(&mut cancelled)?;
            return Ok((next, step));
        }
        let midpoint_increment = self.next_exposure(rate, 0.5 * duration_s)?;
        let (midpoint_fraction, _) = self.fraction_at(midpoint_increment, &mut cancelled)?;
        let midpoint_t = self.law.temperature_at_enthalpy_fraction(
            finite(self.specific_enthalpy_j_kg + 0.5 * heat)?,
            midpoint_fraction,
            solve,
            &mut cancelled,
        )?;
        let midpoint_rate = self.law.rate_s_inv(midpoint_t)?;
        let increment = self.next_exposure(midpoint_rate, duration_s)?;
        let (fraction, saturated) = self.fraction_at(increment, &mut cancelled)?;
        let temperature_k = self.law.temperature_at_enthalpy_fraction(
            final_h_j_kg,
            fraction,
            solve,
            &mut cancelled,
        )?;
        let residual =
            finite(self.law.specific_enthalpy_j_kg(temperature_k, fraction)? - final_h_j_kg)?;
        let mut next = self.clone();
        next.exposure_increment = increment;
        next.crystalline_mass_fraction = fraction;
        next.fully_crystallized = saturated;
        next.specific_enthalpy_j_kg = final_h_j_kg;
        next.temperature_k = temperature_k;
        next.caloric = solve;
        let step = self.report(&next, duration_s, midpoint_rate, residual)?;
        check_cancelled(&mut cancelled)?;
        Ok((next, step))
    }

    fn next_exposure(&self, rate: f64, duration_s: f64) -> Result<f64, CrystallizationError> {
        if self.fully_crystallized {
            Ok(self.exposure_increment)
        } else {
            finite(self.exposure_increment + finite(rate * duration_s)?)
        }
    }

    fn fraction_at(
        &self,
        increment: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(f64, bool), CrystallizationError> {
        let mut total_mass = 0.0;
        let mut fraction = 0.0;
        let mut saturated = true;
        for (index, cohort) in self.cohorts.iter().enumerate() {
            if index % CHECKPOINT_STRIDE == 0 {
                check_cancelled(cancelled)?;
            }
            let x = self
                .law
                .crystalline_mass_fraction(exposure(cohort, increment)?)?;
            let next_mass = add_mass(total_mass, cohort.mass_kg)?;
            fraction = mean(fraction, x, cohort.mass_kg / next_mass)?;
            saturated &= x == self.law.config().maximum_crystalline_mass_fraction;
            total_mass = next_mass;
        }
        Ok((fraction, saturated))
    }

    fn report(
        &self,
        next: &Self,
        duration_s: f64,
        rate: f64,
        residual: f64,
    ) -> Result<CrystallizationPopulationStep, CrystallizationError> {
        Ok(CrystallizationPopulationStep {
            duration_s,
            specific_heat_j_kg: finite(next.specific_enthalpy_j_kg - self.specific_enthalpy_j_kg)?,
            crystallized_mass_kg: finite(
                self.mass_kg * (next.crystalline_mass_fraction - self.crystalline_mass_fraction),
            )?,
            caloric_residual_j_kg: residual,
            midpoint_rate_s_inv: rate,
        })
    }
}

fn exposure(cohort: &Cohort, increment: f64) -> Result<f64, CrystallizationError> {
    if cohort.initially_saturated {
        Ok(cohort.initial_exposure)
    } else {
        finite(cohort.initial_exposure + increment)
    }
}

fn add_mass(previous: f64, mass: f64) -> Result<f64, CrystallizationError> {
    let total = finite(previous + mass)?;
    if total <= previous {
        Err(invalid(
            "crystallization cohort mass increment is unrepresentable",
        ))
    } else {
        Ok(total)
    }
}

fn mean(previous: f64, value: f64, weight: f64) -> Result<f64, CrystallizationError> {
    // Difference form does not multiply a large arbitrary enthalpy datum by
    // independently rounded weights. Exact endpoints retain saturated fractions.
    if weight == 1.0 {
        Ok(value)
    } else if weight <= 0.5 {
        finite(previous + weight * finite(value - previous)?)
    } else {
        finite(value - (1.0 - weight) * finite(value - previous)?)
    }
}

fn invalid(message: &'static str) -> CrystallizationError {
    CrystallizationError::Invalid(message)
}

#[cfg(test)]
#[path = "population/tests.rs"]
mod tests;
