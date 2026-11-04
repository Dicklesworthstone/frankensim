//! Finite-rate crystallization with conservative thermal feedback.
//!
//! The Nakamura/JMAK exposure `z` obeys `dz/dt = K(T)`; relative conversion is
//! `alpha = 1 - exp(-z^n)`. Absolute crystalline MASS fraction is `Xmax*alpha`,
//! not an assumed unit final crystallinity. Integrating exposure starts from
//! zero without an artificial nucleus or the zero-rate trap of an alpha ODE.
//!
//! Original amorphous and crystalline caloric branches define the ONE storage
//! account `h(T,X) = (1-X)*h_a(T) + X*h_c(T)`. A change in X therefore changes
//! temperature at fixed enthalpy. Distinct heat capacities and a temperature-
//! dependent latent difference remain in the original branches; no separate
//! latent heater is added. Both branches must share a formation reference.
//!
//! A thermal step uses explicit midpoint for exposure and inverts this same
//! enthalpy at the midpoint and endpoint. Its input enthalpy is linear in time
//! over that step. This is second order on smooth trajectories, not an exact
//! nonisothermal integration or a temperature-error certificate. The transport
//! owner must refine time for changing rates and heat flux. A rejected stage
//! or cancellation returns no candidate. Scalar work is bounded and allocation-
//! free; the same deterministic math owner supplies all transcendentals.
//!
//! This irreversible isokinetic model is an explicitly admitted approximation.
//! It supplies no stochastic nucleation, crystal morphology, segregation,
//! melting/reset, glass transition, stress locking or experimentally calibrated
//! material constants. An explicitly admitted superheated source domain permits
//! completely amorphous material to cool into crystallization. The supplied rate
//! must be zero at and above melting; any nonzero crystal fraction above that
//! boundary refuses instead of inventing a melt law or resetting phase history.
//! [`EquivalentAgeRemelting`] is a separate opt-in approximation: instantaneous
//! equilibrium melting at the original melting temperature, with the surviving
//! fraction defining an equivalent Avrami age. The default growth APIs remain
//! irreversible; no universal polymer melt-memory kinetics are inferred.
//!
//! Governing model: Nakamura, Katayama & Amano (1973),
//! <https://doi.org/10.1002/app.1973.070170404>. The heat balance and integral
//! exposure appear together in Bhatla & Yao (2008), equations 2–3,
//! <https://www.columbia.edu/~yly1/PDFs3/NAMRC36-.pdf>.

use core::fmt;

#[path = "crystallization/adaptive.rs"]
mod adaptive;
pub use adaptive::{AdaptiveCrystallizationStep, CrystallizationAdaptivity};

#[path = "crystallization/remelting.rs"]
mod remelting;
pub use remelting::{CrystallizationCaloricTransition, EquivalentAgeRemelting};

/// A refusal publishes no physical state or history increment.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum CrystallizationError {
    Invalid(&'static str),
    TemperatureOutOfRange {
        temperature_k: f64,
        minimum_k: f64,
        maximum_k: f64,
    },
    EnthalpyOutOfRange {
        specific_enthalpy_j_kg: f64,
        minimum_j_kg: f64,
        maximum_j_kg: f64,
    },
    NonConverged,
    StateMismatch,
    Cancelled,
}
impl fmt::Display for CrystallizationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "crystallization: {self:?}")
    }
}
impl std::error::Error for CrystallizationError {}

/// Original immutable, strictly increasing caloric branch on the whole domain.
///
/// Amorphous and crystalline charts share the same enthalpy reference. Both
/// need positive heat capacity; the amorphous enthalpy must exceed crystalline
/// enthalpy everywhere admitted. These are constitutive contracts, not facts
/// inferred from material names or endpoint samples. No extrapolation is allowed.
pub trait CrystallizationCaloric {
    fn specific_enthalpy_j_kg(&self, temperature_k: f64) -> Result<f64, CrystallizationError>;
    /// Identity of this immutable ORIGINAL caloric chart, including its energy
    /// datum, complete state domain and constitutive data. Thermal adapters use
    /// this same identity to bind their registered law to this actual branch;
    /// agreement at one temperature is not a source-identity check. A material
    /// name or the combined kinetic-law identity is not a substitute. Generic
    /// closures remain unidentified and cannot claim such a thermal binding.
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        None
    }
    /// Optional positive lower bound on Cp throughout the ENTIRE admitted
    /// domain, J/(kg K). The original source owner proves this bound, never
    /// this solver's endpoint samples. It permits fast inverse interpolation;
    /// absence retains the bounded temperature-bracket fallback.
    fn minimum_heat_capacity_j_kg_k(&self) -> Option<f64> {
        None
    }
}
impl<F: Fn(f64) -> Result<f64, CrystallizationError>> CrystallizationCaloric for F {
    fn specific_enthalpy_j_kg(&self, temperature_k: f64) -> Result<f64, CrystallizationError> {
        self(temperature_k)
    }
}

/// Nonnegative Nakamura rate K(T), in 1/s, over the entire admitted domain.
///
/// This is the nth root of the isothermal Avrami coefficient k(T), whose units
/// are s^-n. The rate is immutable and continuous; the source owner admits its
/// material/process applicability. It must vanish at and above the melting
/// temperature; the solver never replaces an inadmissible supplied rate by zero.
pub trait CrystallizationRate {
    fn rate_s_inv(&self, temperature_k: f64) -> Result<f64, CrystallizationError>;
}
impl<F: Fn(f64) -> Result<f64, CrystallizationError>> CrystallizationRate for F {
    fn rate_s_inv(&self, temperature_k: f64) -> Result<f64, CrystallizationError> {
        self(temperature_k)
    }
}

/// Explicit model parameters. No material identity or source authority is inferred.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CrystallizationConfig {
    /// Caller-owned binding of BOTH caloric sources, the kinetic source, their
    /// common enthalpy reference, applicability and every parameter below.
    pub identity: [u8; 32],
    pub avrami_exponent: f64,
    /// Absolute mass crystallinity at complete relative conversion, in (0,1].
    pub maximum_crystalline_mass_fraction: f64,
    /// Inclusive overlap of both original caloric branches and the kinetic law.
    /// Above melting, only completely amorphous states are admitted. Both
    /// original branches must still cover this full source domain.
    pub temperature_range_k: [f64; 2],
    /// Thermodynamic growth boundary, K. A source domain may extend above it
    /// for amorphous feed, without admitting remelting of existing crystals.
    pub melting_temperature_k: f64,
}

/// Local caloric inversion budget, separate from time-integration accuracy.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CrystallizationSolve {
    pub maximum_iterations: usize,
    /// Absolute reconstruction error in the original specific enthalpy, J/kg.
    /// No relative tolerance is multiplied by an arbitrary formation datum.
    pub specific_enthalpy_tolerance_j_kg: f64,
    /// Final enclosing temperature interval, K.
    pub temperature_tolerance_k: f64,
}
impl Default for CrystallizationSolve {
    fn default() -> Self {
        Self {
            maximum_iterations: 96,
            specific_enthalpy_tolerance_j_kg: 1e-8,
            temperature_tolerance_k: 1e-9,
        }
    }
}

impl CrystallizationSolve {
    /// Validate the original bounded caloric inversion controls without
    /// evaluating a material source or changing a physical state.
    pub fn validate(self) -> Result<(), CrystallizationError> {
        check_solve(self)
    }
}

/// Accepted history bound to one immutable caloric/kinetic law.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CrystallizationState {
    identity: [u8; 32],
    specific_enthalpy_j_kg: f64,
    temperature_k: f64,
    exposure: f64,
    crystalline_mass_fraction: f64,
}
impl CrystallizationState {
    pub const fn identity(self) -> [u8; 32] {
        self.identity
    }
    pub const fn specific_enthalpy_j_kg(self) -> f64 {
        self.specific_enthalpy_j_kg
    }
    pub const fn temperature_k(self) -> f64 {
        self.temperature_k
    }
    pub const fn exposure(self) -> f64 {
        self.exposure
    }
    pub const fn crystalline_mass_fraction(self) -> f64 {
        self.crystalline_mass_fraction
    }
}

/// One proposal; both states and the entire original energy account are retained.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CrystallizationStep {
    pub initial: CrystallizationState,
    pub state: CrystallizationState,
    pub duration_s: f64,
    /// Heat INTO unit mass, h_final-h_initial, J/kg. Includes crystallization once.
    pub specific_heat_j_kg: f64,
    /// Reconstructed branch-mixture enthalpy minus stored final enthalpy, J/kg.
    pub caloric_residual_j_kg: f64,
    pub midpoint_rate_s_inv: f64,
}

/// Reusable immutable constitutive closure; scalar queries allocate nothing.
#[derive(Clone, Copy)]
pub struct CrystallizationLaw<'a> {
    amorphous: &'a dyn CrystallizationCaloric,
    crystalline: &'a dyn CrystallizationCaloric,
    rate: &'a dyn CrystallizationRate,
    config: CrystallizationConfig,
    minimum_heat_capacities_j_kg_k: Option<[f64; 2]>,
}
impl<'a> CrystallizationLaw<'a> {
    pub fn new(
        amorphous: &'a dyn CrystallizationCaloric,
        crystalline: &'a dyn CrystallizationCaloric,
        rate: &'a dyn CrystallizationRate,
        config: CrystallizationConfig,
    ) -> Result<Self, CrystallizationError> {
        let [lower, upper] = config.temperature_range_k;
        if !positive(config.avrami_exponent)
            || !positive(config.maximum_crystalline_mass_fraction)
            || config.maximum_crystalline_mass_fraction > 1.0
            || !positive(lower)
            || !positive(upper)
            || lower >= upper
            || !positive(config.melting_temperature_k)
        {
            return Err(CrystallizationError::Invalid(
                "invalid crystallization parameters or domain",
            ));
        }
        let minimums = [
            amorphous.minimum_heat_capacity_j_kg_k(),
            crystalline.minimum_heat_capacity_j_kg_k(),
        ];
        if minimums
            .into_iter()
            .flatten()
            .any(|minimum| !positive(minimum))
        {
            return Err(CrystallizationError::Invalid(
                "invalid original heat-capacity lower bound",
            ));
        }
        let minimum_heat_capacities_j_kg_k = minimums[0].zip(minimums[1]).map(|(a, c)| [a, c]);
        let law = Self {
            amorphous,
            crystalline,
            rate,
            config,
            minimum_heat_capacities_j_kg_k,
        };
        // Original source owners prove whole-domain contracts. Endpoint checks
        // reject common mismatches without pretending to certify arbitrary closures.
        let at_lower = law.branches(lower)?;
        let at_upper = law.branches(upper)?;
        if at_upper.0 <= at_lower.0 || at_upper.1 <= at_lower.1 {
            return Err(CrystallizationError::Invalid(
                "caloric branches must increase with temperature",
            ));
        }
        law.rate_s_inv(lower)?;
        law.rate_s_inv(upper)?;
        if lower < config.melting_temperature_k && config.melting_temperature_k < upper {
            law.rate_s_inv(config.melting_temperature_k)?;
        }
        Ok(law)
    }

    pub const fn config(&self) -> CrystallizationConfig {
        self.config
    }
    pub const fn identity(&self) -> [u8; 32] {
        self.config.identity
    }

    /// Identity of the original borrowed parent caloric branch, independently
    /// of the combined kinetic-law identity. `None` forbids source-bound thermal
    /// attachment; it does not prevent standalone authored constitutive queries.
    pub fn parent_caloric_identity(&self) -> Option<[u8; 32]> {
        self.amorphous.caloric_identity()
    }

    /// Identity of the original crystalline branch, for adapters that must bind
    /// both caloric sources rather than infer a phase entropy or mechanical law.
    pub fn crystalline_caloric_identity(&self) -> Option<[u8; 32]> {
        self.crystalline.caloric_identity()
    }

    /// Recheck a state against the original law and an explicit caloric budget.
    /// This is useful at cross-module boundaries; states themselves have private
    /// fields and can only be constructed through admitted original-law queries.
    pub fn validate_state(
        &self,
        state: CrystallizationState,
        solve: CrystallizationSolve,
    ) -> Result<(), CrystallizationError> {
        check_solve(solve)?;
        self.check_state(state)?;
        if self.crystalline_mass_fraction(state.exposure)?.to_bits()
            != state.crystalline_mass_fraction.to_bits()
        {
            return Err(CrystallizationError::StateMismatch);
        }
        let reconstructed =
            self.specific_enthalpy_j_kg(state.temperature_k, state.crystalline_mass_fraction)?;
        if (reconstructed - state.specific_enthalpy_j_kg).abs()
            > solve.specific_enthalpy_tolerance_j_kg
        {
            return Err(CrystallizationError::StateMismatch);
        }
        Ok(())
    }

    /// Absolute crystalline mass fraction, retaining small conversions with expm1.
    pub fn crystalline_mass_fraction(&self, exposure: f64) -> Result<f64, CrystallizationError> {
        if !nonnegative(exposure) {
            return Err(CrystallizationError::Invalid(
                "exposure must be finite and nonnegative",
            ));
        }
        let extended = fs_math::det::pow(exposure, self.config.avrami_exponent);
        let conversion = if extended.is_infinite() {
            1.0
        } else {
            -fs_math::det::expm1(-extended)
        };
        let fraction = self.config.maximum_crystalline_mass_fraction * conversion;
        if !nonnegative(fraction) || fraction > self.config.maximum_crystalline_mass_fraction {
            return Err(CrystallizationError::Invalid(
                "crystallinity is not representable",
            ));
        }
        Ok(fraction)
    }

    /// Current original-branch storage, J/kg. Fractions are MASS fractions.
    pub fn specific_enthalpy_j_kg(
        &self,
        temperature_k: f64,
        crystalline_mass_fraction: f64,
    ) -> Result<f64, CrystallizationError> {
        if !nonnegative(crystalline_mass_fraction)
            || crystalline_mass_fraction > self.config.maximum_crystalline_mass_fraction
        {
            return Err(CrystallizationError::Invalid(
                "crystalline mass fraction outside model",
            ));
        }
        self.check_phase_temperature(temperature_k, crystalline_mass_fraction)?;
        let (amorphous, crystalline) = self.branches(temperature_k)?;
        if crystalline_mass_fraction == 0.0 {
            return Ok(amorphous);
        }
        if crystalline_mass_fraction == 1.0 {
            return Ok(crystalline);
        }
        // Difference form avoids separately multiplying two large formation datums.
        finite(amorphous - crystalline_mass_fraction * (amorphous - crystalline))
    }

    pub fn rate_s_inv(&self, temperature_k: f64) -> Result<f64, CrystallizationError> {
        self.check_temperature(temperature_k)?;
        let rate = self.rate.rate_s_inv(temperature_k)?;
        if !nonnegative(rate) {
            return Err(CrystallizationError::Invalid(
                "Nakamura rate must be finite and nonnegative",
            ));
        }
        if temperature_k >= self.config.melting_temperature_k && rate != 0.0 {
            return Err(CrystallizationError::Invalid(
                "crystallization rate must vanish at and above melting",
            ));
        }
        Ok(rate)
    }

    /// Author an initial state by temperature and dimensionless kinetic exposure.
    pub fn state_at_temperature(
        &self,
        temperature_k: f64,
        exposure: f64,
    ) -> Result<CrystallizationState, CrystallizationError> {
        let fraction = self.crystalline_mass_fraction(exposure)?;
        let h = self.specific_enthalpy_j_kg(temperature_k, fraction)?;
        Ok(self.state(h, temperature_k, exposure, fraction))
    }

    /// Invert original mixture enthalpy while retaining the exact supplied account.
    pub fn state_at_enthalpy(
        &self,
        specific_enthalpy_j_kg: f64,
        exposure: f64,
        solve: CrystallizationSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CrystallizationState, CrystallizationError> {
        check_solve(solve)?;
        check_cancelled(&mut cancelled)?;
        let h = finite(specific_enthalpy_j_kg)?;
        let fraction = self.crystalline_mass_fraction(exposure)?;
        let minimum_cp = self
            .minimum_heat_capacities_j_kg_k
            .map(|[a, c]| (1.0 - fraction) * a + fraction * c);
        if minimum_cp.is_some_and(|cp| !positive(cp)) {
            return Err(CrystallizationError::Invalid(
                "mixture heat-capacity bound is not representable",
            ));
        }
        let [mut lower, mut upper] = self.phase_temperature_range_k(fraction)?;
        let mut h_lower = self.specific_enthalpy_j_kg(lower, fraction)?;
        let mut h_upper = self.specific_enthalpy_j_kg(upper, fraction)?;
        if h_lower > h_upper || (lower < upper && h_lower == h_upper) {
            return Err(CrystallizationError::Invalid(
                "mixture caloric branch is not increasing",
            ));
        }
        if h < h_lower || h > h_upper {
            return Err(CrystallizationError::EnthalpyOutOfRange {
                specific_enthalpy_j_kg: h,
                minimum_j_kg: h_lower,
                maximum_j_kg: h_upper,
            });
        }
        if h == h_lower || h == h_upper {
            check_cancelled(&mut cancelled)?;
            return Ok(self.state(
                h,
                if h == h_lower { lower } else { upper },
                exposure,
                fraction,
            ));
        }
        for iteration in 0..solve.maximum_iterations {
            check_cancelled(&mut cancelled)?;
            // False-position interpolation is exact for affine source segments.
            // Bisect every second trial to guarantee bracket progress even for
            // strongly curved calorics. Also bisect at an endpoint or overflow.
            // Missing Cp bounds still require the final temperature bracket.
            let weight = (h - h_lower) / (h_upper - h_lower);
            let interpolated = lower + weight * (upper - lower);
            let temperature = if iteration % 2 == 0
                && interpolated.is_finite()
                && interpolated > lower
                && interpolated < upper
            {
                interpolated
            } else {
                lower + 0.5 * (upper - lower)
            };
            let reconstructed = self.specific_enthalpy_j_kg(temperature, fraction)?;
            if reconstructed < h_lower || reconstructed > h_upper {
                return Err(CrystallizationError::Invalid(
                    "nonmonotone original caloric branch",
                ));
            }
            let residual = (reconstructed - h).abs();
            // This roundoff allowance only guards the shortcut; it NEVER
            // inflates the requested absolute enthalpy tolerance. A large
            // common formation datum therefore cannot purchase cheap accuracy.
            let arithmetic_scale = h
                .abs()
                .max(reconstructed.abs())
                .max(h_lower.abs())
                .max(h_upper.abs());
            let shortcut_temperature_error =
                minimum_cp.map(|cp| (residual + 8.0 * f64::EPSILON * arithmetic_scale) / cp);
            if residual <= solve.specific_enthalpy_tolerance_j_kg
                && (upper - lower <= solve.temperature_tolerance_k
                    || shortcut_temperature_error
                        .is_some_and(|error| error <= solve.temperature_tolerance_k))
            {
                check_cancelled(&mut cancelled)?;
                return Ok(self.state(h, temperature, exposure, fraction));
            }
            if temperature == lower || temperature == upper {
                break;
            }
            if reconstructed < h {
                lower = temperature;
                h_lower = reconstructed;
            } else {
                upper = temperature;
                h_upper = reconstructed;
            }
        }
        Err(CrystallizationError::NonConverged)
    }

    /// Exact frozen-temperature JMAK step, including the thermostat heat required.
    /// Negative returned heat is crystallization energy removed by that thermostat.
    pub fn advance_isothermal(
        &self,
        current: CrystallizationState,
        duration_s: f64,
    ) -> Result<CrystallizationStep, CrystallizationError> {
        self.check_state(current)?;
        check_duration(duration_s)?;
        let rate = self.rate_s_inv(current.temperature_k)?;
        let exposure = self.advance_exposure(current, rate, duration_s)?;
        let state = self.state_at_temperature(current.temperature_k, exposure)?;
        Ok(CrystallizationStep {
            initial: current,
            state,
            duration_s,
            specific_heat_j_kg: finite(
                state.specific_enthalpy_j_kg - current.specific_enthalpy_j_kg,
            )?,
            caloric_residual_j_kg: 0.0,
            midpoint_rate_s_inv: rate,
        })
    }

    /// One midpoint kinetic step with prescribed LINEAR enthalpy evolution.
    ///
    /// Transport supplies the final enthalpy; zero change is adiabatic. Both
    /// intermediate and final temperatures are solved on the actual evolving
    /// mixture, so released latent energy feeds back into the rate. Insufficient
    /// time resolution may produce a refused out-of-domain stage; refine the
    /// timestep instead of clipping temperature, crystallinity or energy.
    pub fn advance_specific_enthalpy(
        &self,
        current: CrystallizationState,
        final_h_j_kg: f64,
        duration_s: f64,
        solve: CrystallizationSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CrystallizationStep, CrystallizationError> {
        self.check_state(current)?;
        check_duration(duration_s)?;
        check_solve(solve)?;
        check_cancelled(&mut cancelled)?;
        let heat = finite(final_h_j_kg - current.specific_enthalpy_j_kg)?;
        let initial_rate = self.rate_s_inv(current.temperature_k)?;
        let midpoint_exposure = self.advance_exposure(current, initial_rate, 0.5 * duration_s)?;
        let midpoint = self.state_at_enthalpy(
            finite(current.specific_enthalpy_j_kg + 0.5 * heat)?,
            midpoint_exposure,
            solve,
            &mut cancelled,
        )?;
        let midpoint_rate = self.rate_s_inv(midpoint.temperature_k)?;
        let exposure = self.advance_exposure(current, midpoint_rate, duration_s)?;
        let state = self.state_at_enthalpy(final_h_j_kg, exposure, solve, &mut cancelled)?;
        let residual = finite(
            self.specific_enthalpy_j_kg(state.temperature_k, state.crystalline_mass_fraction)?
                - final_h_j_kg,
        )?;
        check_cancelled(&mut cancelled)?;
        Ok(CrystallizationStep {
            initial: current,
            state,
            duration_s,
            specific_heat_j_kg: heat,
            caloric_residual_j_kg: residual,
            midpoint_rate_s_inv: midpoint_rate,
        })
    }

    fn advance_exposure(
        &self,
        current: CrystallizationState,
        rate: f64,
        duration_s: f64,
    ) -> Result<f64, CrystallizationError> {
        // Once conversion rounds to its exact upper endpoint, further exposure
        // cannot change any represented material state. Do not overflow an
        // otherwise physically constant saturated state by accumulating its age.
        if current.crystalline_mass_fraction == self.config.maximum_crystalline_mass_fraction {
            return Ok(current.exposure);
        }
        finite(current.exposure + finite(rate * duration_s)?)
    }

    fn state(
        &self,
        h: f64,
        temperature: f64,
        exposure: f64,
        fraction: f64,
    ) -> CrystallizationState {
        CrystallizationState {
            identity: self.config.identity,
            specific_enthalpy_j_kg: h,
            temperature_k: temperature,
            exposure,
            crystalline_mass_fraction: fraction,
        }
    }

    fn check_state(&self, state: CrystallizationState) -> Result<(), CrystallizationError> {
        if state.identity != self.config.identity {
            return Err(CrystallizationError::StateMismatch);
        }
        self.check_phase_temperature(state.temperature_k, state.crystalline_mass_fraction)
    }

    fn check_phase_temperature(
        &self,
        temperature_k: f64,
        crystalline_mass_fraction: f64,
    ) -> Result<(), CrystallizationError> {
        let [minimum_k, maximum_k] = self.phase_temperature_range_k(crystalline_mass_fraction)?;
        if !temperature_k.is_finite() || temperature_k < minimum_k || temperature_k > maximum_k {
            return Err(CrystallizationError::TemperatureOutOfRange {
                temperature_k,
                minimum_k,
                maximum_k,
            });
        }
        Ok(())
    }

    fn phase_temperature_range_k(
        &self,
        crystalline_mass_fraction: f64,
    ) -> Result<[f64; 2], CrystallizationError> {
        let [lower, mut upper] = self.config.temperature_range_k;
        if crystalline_mass_fraction > 0.0 {
            upper = upper.min(self.config.melting_temperature_k);
        }
        if lower > upper {
            return Err(CrystallizationError::Invalid(
                "source domain admits no crystalline temperature at or below melting",
            ));
        }
        Ok([lower, upper])
    }

    fn check_temperature(&self, temperature_k: f64) -> Result<(), CrystallizationError> {
        let [minimum_k, maximum_k] = self.config.temperature_range_k;
        if !temperature_k.is_finite() || temperature_k < minimum_k || temperature_k > maximum_k {
            return Err(CrystallizationError::TemperatureOutOfRange {
                temperature_k,
                minimum_k,
                maximum_k,
            });
        }
        Ok(())
    }

    fn branches(&self, temperature: f64) -> Result<(f64, f64), CrystallizationError> {
        self.check_temperature(temperature)?;
        let amorphous = finite(self.amorphous.specific_enthalpy_j_kg(temperature)?)?;
        let crystalline = finite(self.crystalline.specific_enthalpy_j_kg(temperature)?)?;
        if !positive(amorphous - crystalline) {
            return Err(CrystallizationError::Invalid(
                "positive original crystallization enthalpy difference required",
            ));
        }
        Ok((amorphous, crystalline))
    }
}

/// One source-authored Nakamura rate point, not an Avrami coefficient k in s^-n.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CrystallizationRateKnot {
    pub temperature_k: f64,
    pub rate_s_inv: f64,
}

/// Explicit piecewise-linear rate with strict refusal outside its source domain.
#[derive(Clone, Debug, PartialEq)]
pub struct TabulatedCrystallizationRate {
    knots: Vec<CrystallizationRateKnot>,
}
impl TabulatedCrystallizationRate {
    pub fn new(knots: Vec<CrystallizationRateKnot>) -> Result<Self, CrystallizationError> {
        if knots.len() < 2
            || knots.len() > 65_536
            || knots
                .iter()
                .any(|point| !positive(point.temperature_k) || !nonnegative(point.rate_s_inv))
            || knots
                .windows(2)
                .any(|pair| pair[0].temperature_k >= pair[1].temperature_k)
        {
            return Err(CrystallizationError::Invalid(
                "invalid bounded crystallization-rate table",
            ));
        }
        Ok(Self { knots })
    }
    pub fn knots(&self) -> &[CrystallizationRateKnot] {
        &self.knots
    }
}
impl CrystallizationRate for TabulatedCrystallizationRate {
    fn rate_s_inv(&self, temperature_k: f64) -> Result<f64, CrystallizationError> {
        let minimum_k = self.knots[0].temperature_k;
        let maximum_k = self.knots[self.knots.len() - 1].temperature_k;
        if !temperature_k.is_finite() || temperature_k < minimum_k || temperature_k > maximum_k {
            return Err(CrystallizationError::TemperatureOutOfRange {
                temperature_k,
                minimum_k,
                maximum_k,
            });
        }
        let high = self
            .knots
            .partition_point(|point| point.temperature_k < temperature_k);
        if self.knots[high].temperature_k == temperature_k {
            return Ok(self.knots[high].rate_s_inv);
        }
        let low = self.knots[high - 1];
        let high = self.knots[high];
        let fraction =
            (temperature_k - low.temperature_k) / (high.temperature_k - low.temperature_k);
        finite((1.0 - fraction) * low.rate_s_inv + fraction * high.rate_s_inv)
    }
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn nonnegative(value: f64) -> bool {
    value.is_finite() && value >= 0.0
}
fn finite(value: f64) -> Result<f64, CrystallizationError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(CrystallizationError::Invalid(
            "nonfinite crystallization arithmetic",
        ))
    }
}
fn check_duration(duration_s: f64) -> Result<(), CrystallizationError> {
    if positive(duration_s) {
        Ok(())
    } else {
        Err(CrystallizationError::Invalid(
            "duration must be positive and finite",
        ))
    }
}
fn check_cancelled(cancelled: &mut impl FnMut() -> bool) -> Result<(), CrystallizationError> {
    if cancelled() {
        Err(CrystallizationError::Cancelled)
    } else {
        Ok(())
    }
}
fn check_solve(solve: CrystallizationSolve) -> Result<(), CrystallizationError> {
    if !(1..=256).contains(&solve.maximum_iterations)
        || !positive(solve.specific_enthalpy_tolerance_j_kg)
        || !positive(solve.temperature_tolerance_k)
    {
        Err(CrystallizationError::Invalid(
            "invalid bounded caloric solve controls",
        ))
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn config() -> CrystallizationConfig {
        CrystallizationConfig {
            identity: [0x43; 32],
            avrami_exponent: 2.0,
            maximum_crystalline_mass_fraction: 0.6,
            temperature_range_k: [250.0, 499.0],
            melting_temperature_k: 500.0,
        }
    }
    fn amorphous(t: f64) -> Result<f64, CrystallizationError> {
        Ok(2000.0 * (t - 300.0) + 100_000.0)
    }
    fn crystalline(t: f64) -> Result<f64, CrystallizationError> {
        Ok(1000.0 * (t - 300.0))
    }
    fn near(actual: f64, expected: f64, tolerance: f64) {
        assert!(
            (actual - expected).abs() <= tolerance,
            "{actual} != {expected}"
        );
    }

    #[test]
    fn g1_isothermal_jmak_starts_from_zero_and_accounts_for_thermostat_heat() {
        let rate = |_: f64| Ok(0.1);
        let law = CrystallizationLaw::new(&amorphous, &crystalline, &rate, config()).unwrap();
        let initial = law.state_at_temperature(350.0, 0.0).unwrap();
        let step = law.advance_isothermal(initial, 10.0).unwrap();
        let expected = 0.6 * (1.0 - (-1.0_f64).exp());
        near(step.state.crystalline_mass_fraction(), expected, 2e-15);
        near(step.specific_heat_j_kg, -150_000.0 * expected, 2e-10);
        assert_eq!(step.state.temperature_k(), initial.temperature_k());
        let first = law.advance_isothermal(initial, 4.0).unwrap();
        let second = law.advance_isothermal(first.state, 6.0).unwrap();
        near(second.state.crystalline_mass_fraction(), expected, 2e-15);
        near(
            first.specific_heat_j_kg + second.specific_heat_j_kg,
            step.specific_heat_j_kg,
            2e-10,
        );
    }

    #[test]
    fn g1_adiabatic_recalescence_retains_distinct_heat_capacities() {
        let rate = |_: f64| Ok(0.1);
        let law = CrystallizationLaw::new(&amorphous, &crystalline, &rate, config()).unwrap();
        let initial = law.state_at_temperature(300.0, 0.0).unwrap();
        let step = law
            .advance_specific_enthalpy(
                initial,
                100_000.0,
                10.0,
                CrystallizationSolve::default(),
                || false,
            )
            .unwrap();
        let fraction = 0.6 * (1.0 - (-1.0_f64).exp());
        let expected_t = 300.0 + 100_000.0 * fraction / (2000.0 - 1000.0 * fraction);
        near(step.state.temperature_k(), expected_t, 2e-10);
        assert_eq!(step.state.specific_enthalpy_j_kg(), 100_000.0);
        assert_eq!(step.specific_heat_j_kg, 0.0);
        assert!(step.caloric_residual_j_kg.abs() <= 1e-8);
        assert!(step.state.temperature_k() > initial.temperature_k());
    }

    #[test]
    fn g1_superheated_amorphous_feed_cools_into_growth_on_one_enthalpy_account() {
        let rate = |t: f64| Ok(0.001 * (500.0 - t).max(0.0));
        let cfg = CrystallizationConfig {
            temperature_range_k: [250.0, 650.0],
            ..config()
        };
        let law = CrystallizationLaw::new(&amorphous, &crystalline, &rate, cfg).unwrap();
        let initial = law.state_at_temperature(600.0, 0.0).unwrap();
        let dwell = law.advance_isothermal(initial, 10.0).unwrap();
        assert_eq!(dwell.state, initial);
        assert_eq!(dwell.specific_heat_j_kg, 0.0);
        assert_eq!(law.rate_s_inv(500.0).unwrap(), 0.0);
        assert_eq!(law.rate_s_inv(600.0).unwrap(), 0.0);

        let mut current = initial;
        let mut heat = 0.0;
        for temperature in [550.0, 500.0] {
            let target_h = amorphous(temperature).unwrap();
            let step = law
                .advance_specific_enthalpy(
                    current,
                    target_h,
                    1.0,
                    CrystallizationSolve::default(),
                    || false,
                )
                .unwrap();
            assert_eq!(step.state.exposure(), 0.0);
            assert_eq!(step.state.crystalline_mass_fraction(), 0.0);
            assert_eq!(step.state.specific_enthalpy_j_kg(), target_h);
            near(step.state.temperature_k(), temperature, 1e-9);
            heat += step.specific_heat_j_kg;
            current = step.state;
        }

        // h falls from 500 to 100 kJ/kg. The zero initial rate leaves the
        // midpoint amorphous at 400 K, giving z=0.1 and this independent
        // unequal-Cp mixture inversion. Latent release reduces the cooling;
        // it is already in h, with no second thermal credit.
        let step = law
            .advance_specific_enthalpy(
                current,
                100_000.0,
                1.0,
                CrystallizationSolve::default(),
                || false,
            )
            .unwrap();
        let fraction = 0.6 * (1.0 - (-0.01_f64).exp());
        let expected_t = 300.0 + 100_000.0 * fraction / (2000.0 - 1000.0 * fraction);
        near(step.state.exposure(), 0.1, 1e-11);
        near(step.state.crystalline_mass_fraction(), fraction, 1e-12);
        near(step.state.temperature_k(), expected_t, 1e-9);
        assert!(step.state.temperature_k() > 300.0);
        assert_eq!(step.state.specific_enthalpy_j_kg(), 100_000.0);
        assert!(step.caloric_residual_j_kg.abs() <= 1e-8);
        heat += step.specific_heat_j_kg;
        assert_eq!(heat, -600_000.0);
        assert_eq!(initial.specific_enthalpy_j_kg() + heat, 100_000.0);
    }

    #[test]
    fn g0_superheated_rates_and_source_domains_are_never_silently_extended() {
        let cfg = CrystallizationConfig {
            temperature_range_k: [250.0, 650.0],
            ..config()
        };
        let positive_above_melt = |_: f64| Ok(0.1);
        assert!(
            CrystallizationLaw::new(&amorphous, &crystalline, &positive_above_melt, cfg,).is_err()
        );
        let nonzero_at_melt = |t: f64| Ok(0.001 * (650.0 - t));
        assert!(CrystallizationLaw::new(&amorphous, &crystalline, &nonzero_at_melt, cfg,).is_err());

        // Endpoints cannot certify an arbitrary source. The actual interior
        // query must also refuse a positive superheated rate, not clamp it.
        let bad_interior = |t: f64| {
            Ok(if t <= 500.0 {
                0.001 * (500.0 - t)
            } else {
                1e-6 * (t - 500.0) * (650.0 - t)
            })
        };
        let law = CrystallizationLaw::new(&amorphous, &crystalline, &bad_interior, cfg).unwrap();
        assert!(law.rate_s_inv(550.0).is_err());
        let hot = law.state_at_temperature(550.0, 0.0).unwrap();
        assert!(law.advance_isothermal(hot, 1.0).is_err());

        let rate = |t: f64| Ok(0.001 * (500.0 - t).max(0.0));
        let cold_amorphous = |t: f64| {
            if t > 500.0 {
                Err(CrystallizationError::Invalid("parent source limit"))
            } else {
                amorphous(t)
            }
        };
        let cold_crystalline = |t: f64| {
            if t > 500.0 {
                Err(CrystallizationError::Invalid("crystal source limit"))
            } else {
                crystalline(t)
            }
        };
        assert!(CrystallizationLaw::new(&cold_amorphous, &crystalline, &rate, cfg).is_err());
        assert!(CrystallizationLaw::new(&amorphous, &cold_crystalline, &rate, cfg).is_err());
        let cold_rate = |t: f64| {
            if t > 500.0 {
                Err(CrystallizationError::Invalid("kinetic source limit"))
            } else {
                rate(t)
            }
        };
        assert!(CrystallizationLaw::new(&amorphous, &crystalline, &cold_rate, cfg).is_err());
    }

    #[test]
    fn g0_existing_crystals_cannot_cross_melting_in_an_extended_source_domain() {
        let rate = |t: f64| Ok(0.001 * (500.0 - t).max(0.0));
        let cfg = CrystallizationConfig {
            temperature_range_k: [250.0, 650.0],
            ..config()
        };
        let law = CrystallizationLaw::new(&amorphous, &crystalline, &rate, cfg).unwrap();
        let initial = law.state_at_temperature(450.0, 0.6).unwrap();
        let fraction = initial.crystalline_mass_fraction();
        let hot_h =
            (1.0 - fraction) * amorphous(550.0).unwrap() + fraction * crystalline(550.0).unwrap();
        assert!(law.specific_enthalpy_j_kg(550.0, fraction).is_err());
        assert!(law.state_at_temperature(550.0, initial.exposure()).is_err());
        assert!(matches!(
            law.state_at_enthalpy(
                hot_h,
                initial.exposure(),
                CrystallizationSolve::default(),
                || false
            ),
            Err(CrystallizationError::EnthalpyOutOfRange { .. })
        ));
        assert!(
            law.advance_specific_enthalpy(
                initial,
                hot_h,
                0.01,
                CrystallizationSolve::default(),
                || false,
            )
            .is_err()
        );

        let melt_h = law.specific_enthalpy_j_kg(500.0, fraction).unwrap();
        let at_melt = law
            .state_at_enthalpy(
                melt_h,
                initial.exposure(),
                CrystallizationSolve::default(),
                || false,
            )
            .unwrap();
        assert_eq!(at_melt.temperature_k(), 500.0);
        assert_eq!(at_melt.crystalline_mass_fraction(), fraction);
        assert_eq!(law.advance_isothermal(at_melt, 1.0).unwrap().state, at_melt);
        assert_eq!(initial.temperature_k(), 450.0);
        assert_eq!(initial.exposure(), 0.6);
    }

    #[test]
    fn g1_nonisothermal_latent_feedback_converges_at_second_order() {
        // Independent closed form: equal Cp gives T=300+50*X; n=1 and
        // K(T)=0.1*(350-T)/50 give dz/dt=0.1*exp(-z), hence
        // z=ln(1+0.1*t), X=0.1*t/(1+0.1*t), T=300+50*X.
        let crystal = |t: f64| Ok(2000.0 * (t - 300.0));
        let rate = |t: f64| Ok(0.1 * (350.0 - t) / 50.0);
        let law = CrystallizationLaw::new(
            &amorphous,
            &crystal,
            &rate,
            CrystallizationConfig {
                avrami_exponent: 1.0,
                maximum_crystalline_mass_fraction: 1.0,
                temperature_range_k: [250.0, 350.0],
                melting_temperature_k: 350.0,
                ..config()
            },
        )
        .unwrap();
        let initial = law.state_at_temperature(300.0, 0.0).unwrap();
        let run = |steps: usize| {
            let mut state = initial;
            for _ in 0..steps {
                state = law
                    .advance_specific_enthalpy(
                        state,
                        100_000.0,
                        12.0 / steps as f64,
                        CrystallizationSolve::default(),
                        || false,
                    )
                    .unwrap()
                    .state;
            }
            state
        };
        let expected_exposure = 2.2_f64.ln();
        let expected_fraction = 1.2 / 2.2;
        let coarse = run(8);
        let fine = run(16);
        let coarse_error = (coarse.exposure() - expected_exposure).abs();
        let fine_error = (fine.exposure() - expected_exposure).abs();
        assert!(
            coarse_error > 3.5 * fine_error,
            "{coarse_error} / {fine_error}"
        );
        near(fine.crystalline_mass_fraction(), expected_fraction, 3e-4);
        near(
            fine.temperature_k(),
            300.0 + 50.0 * expected_fraction,
            0.015,
        );
        assert!(
            fine.exposure() < 1.2,
            "recalescence must slow this declared rate"
        );
        near(
            fine.specific_enthalpy_j_kg(),
            initial.specific_enthalpy_j_kg(),
            0.0,
        );
        assert_eq!(run(16), fine);
    }

    #[test]
    fn g0_small_and_saturated_conversion_preserve_bounds() {
        let rate = |_: f64| Ok(0.1);
        let law = CrystallizationLaw::new(&amorphous, &crystalline, &rate, config()).unwrap();
        near(
            law.crystalline_mass_fraction(1e-10).unwrap(),
            0.6e-20,
            1e-34,
        );
        let saturated = law.state_at_temperature(350.0, f64::MAX).unwrap();
        assert_eq!(saturated.crystalline_mass_fraction(), 0.6);
        assert_eq!(
            law.advance_isothermal(saturated, 1.0).unwrap().state,
            saturated
        );
        assert!(law.crystalline_mass_fraction(-1.0).is_err());
    }

    #[test]
    fn g0_tables_domains_controls_and_identity_refuse_without_mutation() {
        let table = TabulatedCrystallizationRate::new(vec![
            CrystallizationRateKnot {
                temperature_k: 250.0,
                rate_s_inv: 0.0,
            },
            CrystallizationRateKnot {
                temperature_k: 350.0,
                rate_s_inv: 0.2,
            },
            CrystallizationRateKnot {
                temperature_k: 500.0,
                rate_s_inv: 0.0,
            },
        ])
        .unwrap();
        near(table.rate_s_inv(300.0).unwrap(), 0.1, 0.0);
        assert!(table.rate_s_inv(501.0).is_err());
        let mut cfg = config();
        cfg.temperature_range_k[1] = 500.0;
        let law = CrystallizationLaw::new(&amorphous, &crystalline, &table, cfg).unwrap();
        let initial = law.state_at_temperature(300.0, 0.0).unwrap();
        assert!(law.state_at_temperature(501.0, 0.0).is_err());
        assert!(matches!(
            law.state_at_enthalpy(1e9, 0.0, CrystallizationSolve::default(), || false),
            Err(CrystallizationError::EnthalpyOutOfRange { .. })
        ));
        assert!(matches!(
            law.state_at_enthalpy(
                100_000.0,
                0.0,
                CrystallizationSolve {
                    maximum_iterations: 1,
                    ..CrystallizationSolve::default()
                },
                || false
            ),
            Err(CrystallizationError::NonConverged)
        ));
        cfg.identity = [0x44; 32];
        let other = CrystallizationLaw::new(&amorphous, &crystalline, &table, cfg).unwrap();
        assert_eq!(
            other.advance_isothermal(initial, 1.0),
            Err(CrystallizationError::StateMismatch)
        );
        let nonzero_at_melt = |_: f64| Ok(0.1);
        assert!(CrystallizationLaw::new(&amorphous, &crystalline, &nonzero_at_melt, cfg).is_err());
        assert_eq!(initial, law.state_at_temperature(300.0, 0.0).unwrap());
    }

    #[test]
    fn g4_cancellation_interrupts_a_local_solve_without_publishing_a_step() {
        let rate = |_: f64| Ok(0.1);
        let law = CrystallizationLaw::new(&amorphous, &crystalline, &rate, config()).unwrap();
        let initial = law.state_at_temperature(300.0, 0.0).unwrap();
        let mut calls = 0;
        let result = law.advance_specific_enthalpy(
            initial,
            100_000.0,
            1.0,
            CrystallizationSolve::default(),
            || {
                calls += 1;
                calls == 7
            },
        );
        assert_eq!(result, Err(CrystallizationError::Cancelled));
        assert_eq!(calls, 7);
        assert_eq!(initial.exposure(), 0.0);
    }

    struct PolynomialCaloric {
        cp_at_300: f64,
        quadratic: f64,
        datum: f64,
        queries: std::cell::Cell<usize>,
    }
    impl CrystallizationCaloric for PolynomialCaloric {
        fn specific_enthalpy_j_kg(&self, t: f64) -> Result<f64, CrystallizationError> {
            self.queries.set(self.queries.get() + 1);
            let delta = t - 300.0;
            Ok(self.datum + delta * (self.cp_at_300 + self.quadratic * delta))
        }
        fn minimum_heat_capacity_j_kg_k(&self) -> Option<f64> {
            // Exact derivative minimum on the test law's whole [250,499] K.
            Some(self.cp_at_300 - 100.0 * self.quadratic)
        }
    }

    #[test]
    fn g1_source_cp_bounds_accelerate_affine_and_nonlinear_original_calorics() {
        let a = PolynomialCaloric {
            cp_at_300: 2000.0,
            quadratic: 0.0,
            datum: 100_000.0,
            queries: std::cell::Cell::new(0),
        };
        let c = PolynomialCaloric {
            cp_at_300: 1000.0,
            quadratic: 0.0,
            datum: 0.0,
            queries: std::cell::Cell::new(0),
        };
        let rate = |_: f64| Ok(0.1);
        let law = CrystallizationLaw::new(&a, &c, &rate, config()).unwrap();
        a.queries.set(0);
        c.queries.set(0);
        let state = law
            .state_at_enthalpy(100_000.0, 0.8, CrystallizationSolve::default(), || false)
            .unwrap();
        let x = 0.6 * (1.0 - (-0.64_f64).exp());
        near(
            state.temperature_k(),
            300.0 + 100_000.0 * x / (2000.0 - 1000.0 * x),
            2e-10,
        );
        assert_eq!(
            a.queries.get(),
            3,
            "two endpoints and one affine inverse query"
        );
        assert_eq!(c.queries.get(), 3);
        law.validate_state(state, CrystallizationSolve::default())
            .unwrap();

        let a = PolynomialCaloric {
            quadratic: 1.5,
            ..a
        };
        let c = PolynomialCaloric {
            cp_at_300: 1200.0,
            quadratic: 0.4,
            ..c
        };
        let law = CrystallizationLaw::new(&a, &c, &rate, config()).unwrap();
        let solve = CrystallizationSolve {
            maximum_iterations: 32,
            ..CrystallizationSolve::default()
        };
        let state = law
            .state_at_enthalpy(160_000.0, 0.8, solve, || false)
            .unwrap();
        let linear = (1.0 - x) * 2000.0 + x * 1200.0;
        let quadratic = (1.0 - x) * 1.5 + x * 0.4;
        let energy = 160_000.0 - (1.0 - x) * 100_000.0;
        let delta = 2.0 * energy / (linear + (linear * linear + 4.0 * quadratic * energy).sqrt());
        near(state.temperature_k(), 300.0 + delta, 2e-10);
        law.validate_state(state, solve).unwrap();
    }

    #[test]
    fn g0_datum_shift_does_not_inflate_absolute_energy_accuracy() {
        let datum = 1e12;
        let a = PolynomialCaloric {
            cp_at_300: 2000.0,
            quadratic: 0.0,
            datum: 100_000.0 + datum,
            queries: std::cell::Cell::new(0),
        };
        let c = PolynomialCaloric {
            cp_at_300: 1000.0,
            quadratic: 0.0,
            datum,
            queries: std::cell::Cell::new(0),
        };
        let rate = |_: f64| Ok(0.1);
        let law = CrystallizationLaw::new(&a, &c, &rate, config()).unwrap();
        let initial = law.state_at_temperature(300.0, 0.0).unwrap();
        let step = law
            .advance_specific_enthalpy(
                initial,
                initial.specific_enthalpy_j_kg(),
                1.0,
                CrystallizationSolve::default(),
                || false,
            )
            .unwrap();
        assert_eq!(step.specific_heat_j_kg, 0.0);
        assert_eq!(step.caloric_residual_j_kg, 0.0);
        assert!(
            a.queries.get() > 6,
            "large common datum must disable the optimistic Cp shortcut"
        );
        assert_eq!(
            step.state.specific_enthalpy_j_kg(),
            initial.specific_enthalpy_j_kg()
        );
    }
}
