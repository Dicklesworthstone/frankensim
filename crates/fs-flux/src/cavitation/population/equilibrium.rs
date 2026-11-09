//! Pressure continuation of the original, fixed-number bubble distribution.
//!
//! Each occupied bin follows its stable static branch without replacing its
//! gas inventory. The same radius supplies volume, acoustic coefficients and
//! nonlinear initial state. This is a sequence of equilibrium queries with
//! frozen liquid properties and an authored polytropic exponent, not a dynamic
//! compression, heat-transfer or bubble-transport process.

use super::super::inception::{
    BubbleEquilibrium, EquilibriumBranch, EquilibriumControls, equilibrium_at_pressure_controlled,
};
use super::super::{BubblePressureSample, BubbleState, CavitationError, poll};
use super::{
    BubblePopulation, BubblePopulationClass, BubblePopulationError, EvolvingBubblePopulation,
    PopulationEvolutionError, PopulationEvolutionLimits, Sum,
};

/// Numerical root evidence in original class order. Empty bins have no root.
#[derive(Clone, Debug, PartialEq)]
pub struct PopulationEquilibriumReport {
    pub absolute_pressure_pa: f64,
    pub equilibria: Vec<Option<BubbleEquilibrium>>,
    pub total_iterations: usize,
    pub maximum_absolute_pressure_residual_pa: f64,
}

/// Complete distribution at one requested pressure; the original is unchanged.
#[derive(Clone, Debug, PartialEq)]
pub struct PopulationEquilibrium {
    pub population: BubblePopulation,
    pub report: PopulationEquilibriumReport,
}

impl BubblePopulationClass {
    /// Follow this class's stable branch, preserving its ORIGINAL model and n.
    ///
    /// This allocates no storage. An empty class stays exactly empty without
    /// requiring its unused model to possess a stable equilibrium. An occupied
    /// class must have strictly positive restoring stiffness; the Blake fold
    /// itself cannot supply a stable acoustic linearization.
    pub fn at_equilibrium_pressure(
        &self,
        absolute_pressure_pa: f64,
        controls: &EquilibriumControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Self, Option<BubbleEquilibrium>), super::super::inception::InceptionError> {
        poll(cancelled)?;
        controls.validate()?;
        if !absolute_pressure_pa.is_finite() {
            return Err(CavitationError::InvalidInput("population equilibrium pressure").into());
        }
        if self.number_density_per_m3 == 0.0 {
            return Ok((*self, None));
        }
        // The radial equation checks the instantaneous wall/far-field jump.
        // A static equilibrium makes that jump zero even after an arbitrarily
        // large ambient shift. Also bound the shift from the ORIGINAL liquid
        // operating point so repeated continuation cannot walk frozen rho/c
        // outside their weak-compressibility approximation.
        let reference = self.bubble.linearization()?;
        let liquid = self.bubble.liquid();
        let maximum_pressure_ratio = self.bubble.compressibility_limits().maximum_pressure_ratio;
        let admit_pressure = |pressure: f64| -> Result<(), CavitationError> {
            let pressure_ratio = (pressure - reference.equilibrium_absolute_pressure_pa).abs()
                / liquid.liquid_density_kg_m3()
                / liquid.liquid_sound_speed_m_s()
                / liquid.liquid_sound_speed_m_s();
            if !pressure_ratio.is_finite() || pressure_ratio > maximum_pressure_ratio {
                return Err(CavitationError::PressureRatioExceeded {
                    pressure_ratio,
                    maximum_pressure_ratio,
                });
            }
            Ok(())
        };
        admit_pressure(absolute_pressure_pa)?;
        let equilibrium = equilibrium_at_pressure_controlled(
            &self.bubble,
            absolute_pressure_pa,
            EquilibriumBranch::Stable,
            controls,
            cancelled,
        )?;
        let linearization = self
            .bubble
            .linearization_at_radius(equilibrium.response.radius_m)?;
        // Root tolerances are numerical controls, not permission to exceed a
        // physical bound. Also admit the actual returned equilibrium, even
        // when an explicitly loose solve leaves a large pressure residual.
        admit_pressure(linearization.equilibrium_absolute_pressure_pa)?;
        let radius = linearization.equilibrium_radius_m;
        let volume = (4.0 * core::f64::consts::PI / 3.0) * radius * radius * radius;
        let void_fraction = self.number_density_per_m3 * volume;
        if !volume.is_finite()
            || volume <= 0.0
            || !void_fraction.is_finite()
            || void_fraction <= 0.0
        {
            return Err(CavitationError::NonFiniteDerived("equilibrated class volume").into());
        }
        poll(cancelled)?;
        Ok((
            Self {
                bubble: self.bubble,
                number_density_per_m3: self.number_density_per_m3,
                void_fraction,
                linearization: Some(linearization),
            },
            Some(equilibrium),
        ))
    }
}

impl BubblePopulation {
    /// Prepare all original nuclei at one absolute pressure, transactionally.
    ///
    /// Each original gas law and number density remain bit-identical; no mean
    /// radius, fitted exponent, newly equilibrated gas inventory or renormalized
    /// number distribution is introduced. The common pressure is satisfied to
    /// the explicit root tolerance, recorded separately for every occupied bin.
    /// Subsequent wave queries reuse these cached native coefficients.
    ///
    /// Work is bounded by the original class budget times the per-class root
    /// budget, with cancellation in every root iteration and before publication.
    /// Below-threshold, out-of-radius, unstable or overfilled candidates fail
    /// without changing any source class. Consumers still enforce their own
    /// stricter dilute-mixture and acoustic applicability limits. The absolute
    /// pressure shift from each ORIGINAL reference equilibrium is bounded by
    /// its weak-compressibility pressure/bulk-modulus ratio; successive queries
    /// cannot accumulate unbounded shifts of the same frozen liquid state.
    pub fn at_equilibrium_pressure(
        &self,
        absolute_pressure_pa: f64,
        controls: &EquilibriumControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PopulationEquilibrium, BubblePopulationError> {
        poll(cancelled)?;
        controls.validate()?;
        if !absolute_pressure_pa.is_finite()
            || self
                .occupied_classes
                .checked_mul(controls.maximum_iterations)
                .is_none()
        {
            return Err(BubblePopulationError::InvalidInput(
                "finite equilibrium pressure and representable total root budget required",
            ));
        }
        let mut classes = Vec::with_capacity(self.classes.len());
        let mut equilibria = Vec::with_capacity(self.classes.len());
        let mut volume = Sum::default();
        let mut total_iterations = 0;
        let mut maximum_absolute_pressure_residual_pa: f64 = 0.0;
        for (class_index, class) in self.classes.iter().enumerate() {
            let (next, equilibrium) = class
                .at_equilibrium_pressure(absolute_pressure_pa, controls, cancelled)
                .map_err(|cause| BubblePopulationError::Equilibrium { class_index, cause })?;
            if let Some(root) = equilibrium {
                total_iterations += root.iterations;
                maximum_absolute_pressure_residual_pa =
                    maximum_absolute_pressure_residual_pa.max(root.pressure_residual_pa.abs());
            }
            volume.add(next.void_fraction);
            classes.push(next);
            equilibria.push(equilibrium);
        }
        let void_fraction = volume.value();
        if !void_fraction.is_finite()
            || void_fraction >= 1.0
            || (self.occupied_classes > 0 && void_fraction <= 0.0)
        {
            return Err(BubblePopulationError::InvalidInput(
                "equilibrated occupied volume must be representable and below one",
            ));
        }
        poll(cancelled)?;
        Ok(PopulationEquilibrium {
            population: Self {
                liquid: self.liquid,
                classes,
                maximum_classes: self.maximum_classes,
                number_density_per_m3: self.number_density_per_m3,
                void_fraction,
                equilibrium_absolute_pressure_pa: (self.occupied_classes > 0)
                    .then_some(absolute_pressure_pa),
                occupied_classes: self.occupied_classes,
            },
            report: PopulationEquilibriumReport {
                absolute_pressure_pa,
                equilibria,
                total_iterations,
                maximum_absolute_pressure_residual_pa,
            },
        })
    }
}

impl PopulationEquilibrium {
    /// Start the ORIGINAL nonlinear models at the prepared equilibrium radii.
    ///
    /// Every slot retains its original gas reference; only physical radius and
    /// wall velocity are initialized. Empty slots retain their unused reference
    /// radius. This is initialization, not an instantaneous reset of an existing
    /// trajectory. The native evolution owner rechecks its entire dilution
    /// envelope and all original model domains before accepting the population.
    pub fn into_evolving(
        self,
        time_s: f64,
        limits: PopulationEvolutionLimits,
    ) -> Result<EvolvingBubblePopulation, PopulationEvolutionError> {
        let states = self
            .population
            .classes()
            .iter()
            .map(|class| {
                let radius = class
                    .linearization()
                    .map_or(class.bubble().gas().reference_radius_m(), |linear| {
                        linear.equilibrium_radius_m
                    });
                BubbleState::try_new(time_s, radius, 0.0)
            })
            .collect::<Result<Vec<_>, _>>()?;
        EvolvingBubblePopulation::try_new(
            self.population,
            time_s,
            states,
            limits,
            BubblePressureSample {
                absolute_pressure_pa: self.report.absolute_pressure_pa,
                time_derivative_pa_s: 0.0,
            },
        )
    }
}

#[cfg(test)]
#[path = "equilibrium/tests.rs"]
mod tests;
