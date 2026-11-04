//! Explicit equilibrium remelting and a declared equivalent-age history rule.
//!
//! At the original Tm, h = h_a(Tm) - X * (h_a(Tm) - h_c(Tm)). Heating can
//! decrease X on this plateau; complete melting returns the original amorphous
//! branch and zero effective exposure. After partial melting, surviving X is
//! represented by the same forward Avrami relation used by the growth owner.
//! This is a chosen reduced history rule, not a claim about retained nuclei,
//! lamellar populations, dissolution rates, molecular memory or stress locking.
//!
//! Growth still uses the original explicit midpoint step and caloric inverse.
//! Midpoint and endpoint states are projected onto the admitted melting plateau
//! without changing stored enthalpy. Smooth unmelted paths retain their original
//! second-order method; transition steps need timestep refinement and carry no
//! second-order event-accuracy claim. An oversized projected midpoint that
//! would erase real growth below Tm refuses rather than accepting false arrest.

use super::{
    CrystallizationConfig, CrystallizationError, CrystallizationLaw, CrystallizationSolve,
    CrystallizationState, CrystallizationStep, check_cancelled, check_duration, check_solve,
    finite,
};

/// Admitted instantaneous melting with equivalent Avrami age after partial melt.
///
/// Construct only through `CrystallizationLaw::equivalent_age_remelting`.
/// The original immutable source/domain binding and latent enthalpy are retained.
/// A full melt erases effective age; this does not model universal melt memory.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EquivalentAgeRemelting {
    config: CrystallizationConfig,
    amorphous_h_at_melting_j_kg: f64,
    crystalline_h_at_melting_j_kg: f64,
    latent_enthalpy_j_kg: f64,
}

impl EquivalentAgeRemelting {
    pub const fn identity(self) -> [u8; 32] {
        self.config.identity
    }

    fn enthalpy_at_melting(self, fraction: f64) -> Result<f64, CrystallizationError> {
        // Same endpoint handling and difference form as the original law.
        if fraction == 0.0 {
            Ok(self.amorphous_h_at_melting_j_kg)
        } else if fraction == 1.0 {
            Ok(self.crystalline_h_at_melting_j_kg)
        } else {
            finite(self.amorphous_h_at_melting_j_kg - fraction * self.latent_enthalpy_j_kg)
        }
    }
}

impl CrystallizationLaw<'_> {
    /// Explicitly select instantaneous equilibrium remelting and equivalent age.
    ///
    /// Both original branches and the kinetic law must admit Tm. Their existing
    /// common-datum, positive-Cp and positive-latent-gap contracts remain intact.
    /// The supplied rate must vanish at Tm; it is never replaced by zero here.
    pub fn equivalent_age_remelting(&self) -> Result<EquivalentAgeRemelting, CrystallizationError> {
        let temperature = self.config.melting_temperature_k;
        self.rate_s_inv(temperature)?;
        let (amorphous, crystalline) = self.branches(temperature)?;
        Ok(EquivalentAgeRemelting {
            config: self.config,
            amorphous_h_at_melting_j_kg: amorphous,
            crystalline_h_at_melting_j_kg: crystalline,
            latent_enthalpy_j_kg: finite(amorphous - crystalline)?,
        })
    }

    /// Prescribed linear enthalpy evolution with explicitly admitted remelting.
    ///
    /// Heat remains `h_final - h_initial`, with no separate latent heater. A
    /// partial melt retains an equivalent age for its surviving crystal MASS
    /// fraction; complete melting sets both fraction and effective age to zero.
    /// Subsequent cooling resumes the original kinetic growth law from that age.
    ///
    /// `NonConverged` also requests a smaller timestep when midpoint projection
    /// would erase growth below Tm. Neither caloric convergence nor a projection
    /// certifies time accuracy through a melting/growth transition. Cancellation,
    /// source/domain loss and a failed represented-age solve publish no step.
    pub fn advance_specific_enthalpy_with_remelting(
        &self,
        current: CrystallizationState,
        final_h_j_kg: f64,
        duration_s: f64,
        solve: CrystallizationSolve,
        remelting: EquivalentAgeRemelting,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CrystallizationStep, CrystallizationError> {
        check_solve(solve)?;
        check_duration(duration_s)?;
        check_cancelled(&mut cancelled)?;
        if remelting.config != self.config {
            return Err(CrystallizationError::StateMismatch);
        }
        self.validate_state(current, solve)?;
        let heat = finite(final_h_j_kg - current.specific_enthalpy_j_kg)?;
        let initial_rate = self.rate_s_inv(current.temperature_k)?;
        if heat == 0.0 && initial_rate == 0.0 {
            // An accepted plateau state may carry a small admitted caloric
            // residual. Re-inverting it can create a spurious undercooling and
            // restart growth during an exactly stationary hold. Preserve the
            // accepted history when both physical derivatives vanish.
            let residual = finite(
                self.specific_enthalpy_j_kg(
                    current.temperature_k,
                    current.crystalline_mass_fraction,
                )? - current.specific_enthalpy_j_kg,
            )?;
            check_cancelled(&mut cancelled)?;
            return Ok(CrystallizationStep {
                initial: current,
                state: current,
                duration_s,
                specific_heat_j_kg: 0.0,
                caloric_residual_j_kg: residual,
                midpoint_rate_s_inv: 0.0,
            });
        }
        let midpoint_exposure = self.advance_exposure(current, initial_rate, 0.5 * duration_s)?;
        let midpoint_h = finite(current.specific_enthalpy_j_kg + 0.5 * heat)?;
        let (midpoint, midpoint_projected) = self.remelting_state_at_enthalpy(
            midpoint_h,
            midpoint_exposure,
            solve,
            remelting,
            &mut cancelled,
        )?;
        let midpoint_rate = self.rate_s_inv(midpoint.temperature_k)?;
        let exposure = self.advance_exposure(current, midpoint_rate, duration_s)?;
        let (state, _) = self.remelting_state_at_enthalpy(
            final_h_j_kg,
            exposure,
            solve,
            remelting,
            &mut cancelled,
        )?;
        if midpoint_projected
            && initial_rate > 0.0
            && midpoint_rate == 0.0
            && state.crystalline_mass_fraction <= current.crystalline_mass_fraction
            && state.temperature_k < self.config.melting_temperature_k
        {
            // In an adiabatic recalescence overshoot, clipping the midpoint to
            // Tm makes Kmid=0. Accepting the resulting old exposure would freeze
            // the same still-undercooled state on every repeated timestep.
            return Err(CrystallizationError::NonConverged);
        }
        let residual = finite(
            self.specific_enthalpy_j_kg(state.temperature_k, state.crystalline_mass_fraction)?
                - final_h_j_kg,
        )?;
        if residual.abs() > solve.specific_enthalpy_tolerance_j_kg {
            return Err(CrystallizationError::NonConverged);
        }
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

    fn remelting_state_at_enthalpy(
        &self,
        h: f64,
        exposure: f64,
        solve: CrystallizationSolve,
        remelting: EquivalentAgeRemelting,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(CrystallizationState, bool), CrystallizationError> {
        check_cancelled(cancelled)?;
        let h = finite(h)?;
        let fraction = self.crystalline_mass_fraction(exposure)?;
        if h >= remelting.amorphous_h_at_melting_j_kg {
            return self
                .state_at_enthalpy(h, 0.0, solve, cancelled)
                .map(|state| (state, exposure != 0.0));
        }
        if fraction == 0.0 || h <= remelting.enthalpy_at_melting(fraction)? {
            return self
                .state_at_enthalpy(h, exposure, solve, cancelled)
                .map(|state| (state, false));
        }

        // Reconstruct effective age, then ALWAYS obtain the stored fraction
        // from the original forward relation. An independently rounded X would
        // violate state validation and change history on the next growth step.
        let target_fraction =
            finite(remelting.amorphous_h_at_melting_j_kg - h)? / remelting.latent_enthalpy_j_kg;
        let maximum = self.config.maximum_crystalline_mass_fraction;
        let relative = target_fraction / maximum;
        let extended = if relative < 0.5 {
            let remainder = 1.0 - relative;
            if remainder == 1.0 {
                relative
            } else {
                // Compensate the rounded subtraction before the shared ln.
                -fs_math::det::ln(remainder) * (relative / (1.0 - remainder))
            }
        } else {
            -fs_math::det::ln((maximum - target_fraction) / maximum)
        };
        let guess = fs_math::det::pow(extended, self.config.avrami_exponent.recip());
        // Positive finite f64 bit patterns are ordered by value. Bisection in
        // this ordering reaches adjacent representable ages in at most 63
        // iterations, even if the original saturated exposure is f64::MAX.
        // Ordinary affine latent closure normally accepts the analytic guess.
        let mut lower = 0_u64;
        let mut upper = exposure.to_bits();
        for iteration in 0..solve.maximum_iterations {
            check_cancelled(cancelled)?;
            let age = if iteration == 0 && guess.is_finite() && guess > 0.0 && guess < exposure {
                guess
            } else {
                f64::from_bits(lower + (upper - lower) / 2)
            };
            let surviving = self.crystalline_mass_fraction(age)?;
            let reconstructed = remelting.enthalpy_at_melting(surviving)?;
            if surviving <= fraction
                && (reconstructed - h).abs() <= solve.specific_enthalpy_tolerance_j_kg
            {
                check_cancelled(cancelled)?;
                return Ok((
                    self.state(h, self.config.melting_temperature_k, age, surviving),
                    true,
                ));
            }
            let bits = age.to_bits();
            if bits <= lower || bits >= upper {
                break;
            }
            if reconstructed > h {
                lower = bits;
            } else {
                upper = bits;
            }
        }
        Err(CrystallizationError::NonConverged)
    }
}

#[cfg(test)]
#[path = "remelting/tests.rs"]
mod tests;
