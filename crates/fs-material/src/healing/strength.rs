//! An explicit measured-or-declared relation from mobility dose to strength.
//! Data supplied here are not automatically qualified as calibrated evidence.

use super::{ContactHistory, HealingError, nonnegative, positive, validate_identity};
use crate::cohesive::BilinearCohesive;

/// Independent physical strength/fracture measurements or declared model inputs.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StrengthKnot {
    pub dose: f64,
    pub peak_traction_pa: f64,
    pub fracture_energy_j_m2: f64,
}

/// A frozen interpretation of exposure. Linear interpolation only inside the
/// declared dose range. Above-range dose refuses; there is no silent saturation.
/// K and the compression law are explicitly supplied, not inferred from heat flow.
#[derive(Debug, Clone)]
pub struct StrengthMap {
    identity: String,
    mobility_identity: String,
    knots: Vec<StrengthKnot>,
    stiffness_pa_m: f64,
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
    pub fn new(identity: String, mobility_identity: String, knots: Vec<StrengthKnot>,
        stiffness_pa_m: f64, compression_pa_m: f64, separation_limit_m: f64,
    ) -> Result<Self, HealingError> {
        validate_identity(&identity)?;
        validate_identity(&mobility_identity)?;
        positive(stiffness_pa_m)?;
        positive(compression_pa_m)?;
        positive(separation_limit_m)?;
        if !(2..=1024).contains(&knots.len()) || knots[0].dose != 0.0
            || knots[0].peak_traction_pa != 0.0 || knots[0].fracture_energy_j_m2 != 0.0 {
            return Err(HealingError::Invalid("strength map starts with (dose,strength,energy)=(0,0,0)"));
        }
        for (i,knot) in knots.iter().enumerate() {
            nonnegative(knot.dose)?;
            nonnegative(knot.peak_traction_pa)?;
            nonnegative(knot.fracture_energy_j_m2)?;
            if i > 0 && (knots[i-1].dose >= knot.dose
                || knots[i-1].peak_traction_pa > knot.peak_traction_pa
                || knots[i-1].fracture_energy_j_m2 > knot.fracture_energy_j_m2) {
                return Err(HealingError::Invalid("dose increases and strength/energy cannot decrease"));
            }
            if knot.peak_traction_pa == 0.0 {
                if knot.fracture_energy_j_m2 != 0.0 {
                    return Err(HealingError::Invalid("zero-strength knot must have zero fracture energy"));
                }
            } else {
                BilinearCohesive::new(stiffness_pa_m,knot.peak_traction_pa,
                    knot.fracture_energy_j_m2,compression_pa_m,separation_limit_m)
                    .map_err(|_| HealingError::Invalid("strength knot is not an admissible cohesive law"))?;
            }
        }
        Ok(Self { identity, mobility_identity, knots, stiffness_pa_m,
            compression_pa_m, separation_limit_m })
    }

    pub fn freeze(&self, history: &ContactHistory<'_>) -> Result<FrozenBond, HealingError> {
        if history.law_identity() != self.mobility_identity {
            return Err(HealingError::Invalid("strength map belongs to a different mobility law"));
        }
        let dose = history.dose();
        if dose > self.knots[self.knots.len()-1].dose { return Err(HealingError::Dose(dose)); }
        let upper = self.knots.partition_point(|k| k.dose < dose);
        let (strength,energy) = if upper == 0 { (0.0,0.0) } else {
            let a = self.knots[upper-1];
            let b = self.knots[upper];
            let w = (dose-a.dose)/(b.dose-a.dose);
            ((1.0-w)*a.peak_traction_pa+w*b.peak_traction_pa,
             (1.0-w)*a.fracture_energy_j_m2+w*b.fracture_energy_j_m2)
        };
        let law = if strength == 0.0 {
            if energy != 0.0 { return Err(HealingError::Invalid("unrepresentable zero-strength interpolation")); }
            None
        } else {
            Some(BilinearCohesive::new(self.stiffness_pa_m,strength,energy,
                self.compression_pa_m,self.separation_limit_m)
                .map_err(|_| HealingError::Invalid("interpolated strength is not an admissible cohesive law"))?)
        };
        Ok(FrozenBond { map_identity: self.identity.clone(),
            mobility_identity: self.mobility_identity.clone(), dose, law })
    }
}
