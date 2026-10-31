//! Frozen number-density distribution of original spherical bubble owners.
//!
//! Each class retains its Keller--Miksis gas inventory, interface coefficients,
//! reference radius and native linearization. Number density is an absolute
//! number per mixture volume, not a probability or a volume-fraction weight.
//! The occupied volume is sum_j (4 pi/3) n_j R_j^3. Inputs are never normalized,
//! merged into a mean radius, or assigned an inferred gas law.
//!
//! This represents a specified equilibrium distribution, not a population-
//! balance model for nucleation, breakup, coalescence, diffusion or transport.
//! Independent consumers must still admit dilution and bubble interactions.
//! Commander and Prosperetti (1989), equations 13--14, use this same fixed
//! number distribution and a separate radial equation for each equilibrium size:
//! <https://pages.jh.edu/aprospe1/publications/PapersPublished/BubblyLiquids/Commander.pdf>.

use super::{BubbleLinearization, CavitationError, CavitationState, KellerMiksisBubble};

#[path = "population/evolution.rs"]
mod evolution;
pub use evolution::{
    EvolvingBubblePopulation, PopulationAdvance, PopulationEvolutionError,
    PopulationEvolutionLimits, PopulationEvolutionNumerics, PopulationEvolutionStep,
    PopulationObservation,
};

/// Invalid population data or refusal by an original radial model.
#[derive(Clone, Debug, PartialEq)]
pub enum BubblePopulationError {
    Bubble(CavitationError),
    InvalidInput(&'static str),
    ClassBudgetExceeded {
        requested: usize,
        maximum: usize,
    },
    MaterialMismatch {
        class_index: usize,
    },
    EquilibriumMismatch {
        class_index: usize,
        reference_pressure_pa: f64,
        class_pressure_pa: f64,
    },
}
impl core::fmt::Display for BubblePopulationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "bubble population refused: {self:?}")
    }
}
impl std::error::Error for BubblePopulationError {}
impl From<CavitationError> for BubblePopulationError {
    fn from(error: CavitationError) -> Self {
        Self::Bubble(error)
    }
}

/// One original bubble model and its absolute equilibrium number density.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BubblePopulationClass {
    bubble: KellerMiksisBubble,
    number_density_per_m3: f64,
    void_fraction: f64,
    linearization: Option<BubbleLinearization>,
}
impl BubblePopulationClass {
    pub fn try_new(
        bubble: KellerMiksisBubble,
        number_density_per_m3: f64,
    ) -> Result<Self, BubblePopulationError> {
        if !number_density_per_m3.is_finite() || number_density_per_m3 < 0.0 {
            return Err(BubblePopulationError::InvalidInput(
                "number density must be finite and nonnegative",
            ));
        }
        let (void_fraction, linearization) = if number_density_per_m3 == 0.0 {
            // No occupied bubble exists, so no fictitious resonance or volume
            // arithmetic restricts the exact empty-population endpoint.
            (0.0, None)
        } else {
            let linear = bubble.linearization()?;
            let radius = linear.equilibrium_radius_m;
            let volume = (4.0 * core::f64::consts::PI / 3.0) * radius * radius * radius;
            let fraction = number_density_per_m3 * volume;
            if !volume.is_finite() || volume <= 0.0 || !fraction.is_finite() || fraction <= 0.0 {
                return Err(BubblePopulationError::InvalidInput(
                    "occupied bubble volume and void fraction must be representable",
                ));
            }
            (fraction, Some(linear))
        };
        Ok(Self {
            bubble,
            number_density_per_m3,
            void_fraction,
            linearization,
        })
    }

    pub const fn bubble(&self) -> &KellerMiksisBubble {
        &self.bubble
    }
    pub const fn number_density_per_m3(&self) -> f64 {
        self.number_density_per_m3
    }
    pub const fn void_fraction(&self) -> f64 {
        self.void_fraction
    }
    pub const fn linearization(&self) -> Option<BubbleLinearization> {
        self.linearization
    }
}

/// An admitted, immutable distribution in one original liquid and ambient state.
/// Preparation consumes the supplied classes and performs no sorting or fitting.
/// The caller's explicit class budget bounds all subsequent traversals.
#[derive(Clone, Debug, PartialEq)]
pub struct BubblePopulation {
    liquid: CavitationState,
    classes: Vec<BubblePopulationClass>,
    maximum_classes: usize,
    number_density_per_m3: f64,
    void_fraction: f64,
    equilibrium_absolute_pressure_pa: Option<f64>,
    occupied_classes: usize,
}
impl BubblePopulation {
    pub fn try_new(
        liquid: CavitationState,
        classes: Vec<BubblePopulationClass>,
        maximum_classes: usize,
    ) -> Result<Self, BubblePopulationError> {
        if classes.len() > maximum_classes {
            return Err(BubblePopulationError::ClassBudgetExceeded {
                requested: classes.len(),
                maximum: maximum_classes,
            });
        }
        let mut number = Sum::default();
        let mut volume = Sum::default();
        let mut reference: Option<(f64, f64)> = None;
        let mut occupied_classes = 0;
        for (class_index, class) in classes.iter().enumerate() {
            if class.bubble.liquid() != liquid {
                return Err(BubblePopulationError::MaterialMismatch { class_index });
            }
            let Some(linear) = class.linearization else {
                continue;
            };
            let scale = equilibrium_scale(&class.bubble, linear);
            if let Some((pressure, reference_scale)) = reference {
                // Only roundoff in reconstructing pg0 + pv - 2 sigma/R0 is
                // tolerated. This is not a physical pressure-spread fit.
                let tolerance = 64.0 * f64::EPSILON * scale.max(reference_scale);
                if (linear.equilibrium_absolute_pressure_pa - pressure).abs() > tolerance {
                    return Err(BubblePopulationError::EquilibriumMismatch {
                        class_index,
                        reference_pressure_pa: pressure,
                        class_pressure_pa: linear.equilibrium_absolute_pressure_pa,
                    });
                }
            } else {
                reference = Some((linear.equilibrium_absolute_pressure_pa, scale));
            }
            occupied_classes += 1;
            number.add(class.number_density_per_m3);
            volume.add(class.void_fraction);
        }
        let number_density_per_m3 = number.value();
        let void_fraction = volume.value();
        if !number_density_per_m3.is_finite()
            || !void_fraction.is_finite()
            || void_fraction >= 1.0
            || (occupied_classes > 0 && (number_density_per_m3 <= 0.0 || void_fraction <= 0.0))
        {
            return Err(BubblePopulationError::InvalidInput(
                "population number/volume totals must be representable with void fraction below one",
            ));
        }
        Ok(Self {
            liquid,
            classes,
            maximum_classes,
            number_density_per_m3,
            void_fraction,
            equilibrium_absolute_pressure_pa: reference.map(|(pressure, _)| pressure),
            occupied_classes,
        })
    }

    pub const fn liquid(&self) -> CavitationState {
        self.liquid
    }
    pub fn classes(&self) -> &[BubblePopulationClass] {
        &self.classes
    }
    pub const fn maximum_classes(&self) -> usize {
        self.maximum_classes
    }
    pub const fn occupied_classes(&self) -> usize {
        self.occupied_classes
    }
    pub const fn number_density_per_m3(&self) -> f64 {
        self.number_density_per_m3
    }
    pub const fn void_fraction(&self) -> f64 {
        self.void_fraction
    }
    pub const fn equilibrium_absolute_pressure_pa(&self) -> Option<f64> {
        self.equilibrium_absolute_pressure_pa
    }
}

fn equilibrium_scale(bubble: &KellerMiksisBubble, linear: BubbleLinearization) -> f64 {
    bubble
        .gas()
        .reference_pressure_pa()
        .abs()
        .max(bubble.liquid().saturation_pressure_pa())
        .max(2.0 * bubble.liquid().surface_tension_n_m() / linear.equilibrium_radius_m)
        .max(linear.equilibrium_absolute_pressure_pa.abs())
        .max(1.0)
}

// Neumaier summation retains small classes alongside much larger populations.
// The supplied order is fixed; no nondeterministic or parallel reduction occurs.
#[derive(Default)]
struct Sum {
    sum: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) {
        let next = self.sum + value;
        self.correction += if self.sum.abs() >= value.abs() {
            (self.sum - next) + value
        } else {
            (value - next) + self.sum
        };
        self.sum = next;
    }
    fn value(&self) -> f64 {
        self.sum + self.correction
    }
}

#[cfg(test)]
#[path = "population/tests.rs"]
mod tests;
