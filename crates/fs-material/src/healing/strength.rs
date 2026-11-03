//! An explicit measured-or-declared relation from mobility dose to strength.
//! Data supplied here are not automatically qualified as calibrated evidence.

use super::{ContactHistory, HealingError, nonnegative, positive, validate_identity};
use crate::cohesive::BilinearCohesive;
use crate::cohesive::mixed_mode::MixedModeParameters;

/// Independent physical strength/fracture measurements or declared model inputs.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StrengthKnot {
    pub dose: f64,
    pub peak_traction_pa: f64,
    pub fracture_energy_j_m2: f64,
}

/// Independently supplied opening and shear observations at one mobility dose.
/// Both shear directions share the declared shear law; no opening/shear ratio
/// is inferred from mobility, temperature, contact conductance or bulk strength.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MixedModeStrengthKnot {
    pub dose: f64,
    pub normal_strength_pa: f64,
    pub tangential_strength_pa: f64,
    pub mode_i_fracture_energy_j_m2: f64,
    pub mode_ii_fracture_energy_j_m2: f64,
}

/// Dose-independent parts of an explicitly selected mixed-mode strength map.
/// The original cohesive owner supplies quadratic initiation and frozen-onset
/// Benzeggagh--Kenane toughness, with its original validity checks.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MixedModeStrengthConfig {
    pub normal_stiffness_pa_m: f64,
    pub tangential_stiffness_pa_m: f64,
    /// At least one, as required by the original mixed-mode cohesive law.
    pub bk_exponent: f64,
    pub compression_stiffness_pa_m: f64,
    pub maximum_separation_m: f64,
}

#[derive(Debug, Clone)]
enum StrengthEnvelope {
    Common {
        knots: Vec<StrengthKnot>,
        stiffness_pa_m: f64,
    },
    MixedMode {
        knots: Vec<MixedModeStrengthKnot>,
        config: MixedModeStrengthConfig,
    },
}

/// A frozen interpretation of exposure. Linear interpolation only inside the
/// declared dose range. Above-range dose refuses; there is no silent saturation.
/// K and the compression law are explicitly supplied, not inferred from heat flow.
#[derive(Debug, Clone)]
pub struct StrengthMap {
    identity: String,
    mobility_identity: String,
    envelope: StrengthEnvelope,
    compression_pa_m: f64,
    separation_limit_m: f64,
}

/// Fresh, zero-damage cohesive law for a separate mechanical stage, together
/// with the exposure/mapping identities. Never use this to reset an existing
/// damaged cohesive history or to claim simultaneous healing and fracture.
#[derive(Debug, Clone)]
pub struct FrozenBond {
    pub map_identity: String,
    pub mobility_identity: String,
    pub dose: f64,
    /// None denotes no developed mechanical bond under the supplied map.
    pub law: Option<BilinearCohesive>,
}

impl StrengthMap {
    pub fn new(
        identity: String,
        mobility_identity: String,
        knots: Vec<StrengthKnot>,
        stiffness_pa_m: f64,
        compression_pa_m: f64,
        separation_limit_m: f64,
    ) -> Result<Self, HealingError> {
        validate_identity(&identity)?;
        validate_identity(&mobility_identity)?;
        positive(stiffness_pa_m)?;
        positive(compression_pa_m)?;
        positive(separation_limit_m)?;
        if !(2..=1024).contains(&knots.len())
            || knots[0].dose != 0.0
            || knots[0].peak_traction_pa != 0.0
            || knots[0].fracture_energy_j_m2 != 0.0
        {
            return Err(HealingError::Invalid(
                "strength map starts with (dose,strength,energy)=(0,0,0)",
            ));
        }
        for (i, knot) in knots.iter().enumerate() {
            nonnegative(knot.dose)?;
            nonnegative(knot.peak_traction_pa)?;
            nonnegative(knot.fracture_energy_j_m2)?;
            if i > 0
                && (knots[i - 1].dose >= knot.dose
                    || knots[i - 1].peak_traction_pa > knot.peak_traction_pa
                    || knots[i - 1].fracture_energy_j_m2 > knot.fracture_energy_j_m2)
            {
                return Err(HealingError::Invalid(
                    "dose increases and strength/energy cannot decrease",
                ));
            }
            if knot.peak_traction_pa == 0.0 {
                if knot.fracture_energy_j_m2 != 0.0 {
                    return Err(HealingError::Invalid(
                        "zero-strength knot must have zero fracture energy",
                    ));
                }
            } else {
                BilinearCohesive::new(
                    stiffness_pa_m,
                    knot.peak_traction_pa,
                    knot.fracture_energy_j_m2,
                    compression_pa_m,
                    separation_limit_m,
                )
                .map_err(|_| {
                    HealingError::Invalid("strength knot is not an admissible cohesive law")
                })?;
            }
        }
        Ok(Self {
            identity,
            mobility_identity,
            envelope: StrengthEnvelope::Common {
                knots,
                stiffness_pa_m,
            },
            compression_pa_m,
            separation_limit_m,
        })
    }

    /// Freeze distinct opening/shear strengths and fracture work from the SAME
    /// accepted contact exposure. Each column interpolates linearly only within
    /// the supplied dose domain. Zero dose has no bond; later all-zero knots may
    /// explicitly retain an induction interval. A developed knot requires all
    /// four observations and an admissible original mixed-mode envelope.
    ///
    /// This creates a law for a new mechanical stage. It does not evolve an
    /// already damaged bond, establish a contact area, or infer calibration.
    pub fn new_mixed_mode(
        identity: String,
        mobility_identity: String,
        knots: Vec<MixedModeStrengthKnot>,
        config: MixedModeStrengthConfig,
    ) -> Result<Self, HealingError> {
        validate_identity(&identity)?;
        validate_identity(&mobility_identity)?;
        for value in [
            config.normal_stiffness_pa_m,
            config.tangential_stiffness_pa_m,
            config.bk_exponent,
            config.compression_stiffness_pa_m,
            config.maximum_separation_m,
        ] {
            positive(value)?;
        }
        if config.bk_exponent < 1.0 {
            return Err(HealingError::Invalid("BK exponent must be at least one"));
        }
        if !(2..=1024).contains(&knots.len())
            || knots[0].dose != 0.0
            || knot_values(knots[0]) != [0.0; 4]
        {
            return Err(HealingError::Invalid(
                "mixed-mode strength map starts with zero dose, strengths and energies",
            ));
        }
        for (i, &knot) in knots.iter().enumerate() {
            nonnegative(knot.dose)?;
            for value in knot_values(knot) {
                nonnegative(value)?;
            }
            if i > 0
                && (knots[i - 1].dose >= knot.dose
                    || knot_values(knots[i - 1])
                        .iter()
                        .zip(knot_values(knot))
                        .any(|(&previous, current)| previous > current))
            {
                return Err(HealingError::Invalid(
                    "dose increases and mode strengths/energies cannot decrease",
                ));
            }
            mixed_law(knot, config)?;
        }
        Ok(Self {
            identity,
            mobility_identity,
            envelope: StrengthEnvelope::MixedMode { knots, config },
            compression_pa_m: config.compression_stiffness_pa_m,
            separation_limit_m: config.maximum_separation_m,
        })
    }

    pub fn freeze(&self, history: &ContactHistory<'_>) -> Result<FrozenBond, HealingError> {
        if history.law_identity() != self.mobility_identity {
            return Err(HealingError::Invalid(
                "strength map belongs to a different mobility law",
            ));
        }
        let dose = history.dose();
        let law = match &self.envelope {
            StrengthEnvelope::Common {
                knots,
                stiffness_pa_m,
            } => self.freeze_common(knots, *stiffness_pa_m, dose)?,
            StrengthEnvelope::MixedMode { knots, config } => {
                if dose > knots[knots.len() - 1].dose {
                    return Err(HealingError::Dose(dose));
                }
                let upper = knots.partition_point(|k| k.dose < dose);
                let values = if upper == 0 {
                    [0.0; 4]
                } else {
                    let a = knots[upper - 1];
                    let b = knots[upper];
                    let w = (dose - a.dose) / (b.dose - a.dose);
                    let from = knot_values(a);
                    let to = knot_values(b);
                    core::array::from_fn(|i| (1.0 - w) * from[i] + w * to[i])
                };
                mixed_law(
                    MixedModeStrengthKnot {
                        dose,
                        normal_strength_pa: values[0],
                        tangential_strength_pa: values[1],
                        mode_i_fracture_energy_j_m2: values[2],
                        mode_ii_fracture_energy_j_m2: values[3],
                    },
                    *config,
                )?
            }
        };
        Ok(FrozenBond {
            map_identity: self.identity.clone(),
            mobility_identity: self.mobility_identity.clone(),
            dose,
            law,
        })
    }

    fn freeze_common(
        &self,
        knots: &[StrengthKnot],
        stiffness_pa_m: f64,
        dose: f64,
    ) -> Result<Option<BilinearCohesive>, HealingError> {
        if dose > knots[knots.len() - 1].dose {
            return Err(HealingError::Dose(dose));
        }
        let upper = knots.partition_point(|k| k.dose < dose);
        let (strength, energy) = if upper == 0 {
            (0.0, 0.0)
        } else {
            let a = knots[upper - 1];
            let b = knots[upper];
            let w = (dose - a.dose) / (b.dose - a.dose);
            (
                (1.0 - w) * a.peak_traction_pa + w * b.peak_traction_pa,
                (1.0 - w) * a.fracture_energy_j_m2 + w * b.fracture_energy_j_m2,
            )
        };
        if strength == 0.0 {
            if energy != 0.0 {
                return Err(HealingError::Invalid(
                    "unrepresentable zero-strength interpolation",
                ));
            }
            Ok(None)
        } else {
            Ok(Some(
                BilinearCohesive::new(
                    stiffness_pa_m,
                    strength,
                    energy,
                    self.compression_pa_m,
                    self.separation_limit_m,
                )
                .map_err(|_| {
                    HealingError::Invalid("interpolated strength is not an admissible cohesive law")
                })?,
            ))
        }
    }
}

fn knot_values(knot: MixedModeStrengthKnot) -> [f64; 4] {
    [
        knot.normal_strength_pa,
        knot.tangential_strength_pa,
        knot.mode_i_fracture_energy_j_m2,
        knot.mode_ii_fracture_energy_j_m2,
    ]
}

fn mixed_law(
    knot: MixedModeStrengthKnot,
    config: MixedModeStrengthConfig,
) -> Result<Option<BilinearCohesive>, HealingError> {
    if knot_values(knot) == [0.0; 4] {
        return Ok(None);
    }
    let parameters = MixedModeParameters {
        normal_stiffness_pa_m: config.normal_stiffness_pa_m,
        tangential_stiffness_pa_m: config.tangential_stiffness_pa_m,
        normal_strength_pa: knot.normal_strength_pa,
        tangential_strength_pa: knot.tangential_strength_pa,
        mode_i_fracture_energy_j_m2: knot.mode_i_fracture_energy_j_m2,
        mode_ii_fracture_energy_j_m2: knot.mode_ii_fracture_energy_j_m2,
        bk_exponent: config.bk_exponent,
    };
    BilinearCohesive::new_mixed_mode(
        parameters,
        config.compression_stiffness_pa_m,
        config.maximum_separation_m,
    )
    .map(Some)
    .map_err(|_| HealingError::Invalid("mixed-mode strength is not an admissible cohesive law"))
}

#[cfg(test)]
mod tests;
