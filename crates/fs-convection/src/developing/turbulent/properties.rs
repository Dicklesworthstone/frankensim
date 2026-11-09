//! Explicit wall-property multipliers of the original Gnielinski tube model.
//!
//! Gnielinski (2013), Eqs. (3) and (4), doi:10.1016/j.ijheatmasstransfer.2013.04.015.
//! These correct the bulk-property result, not its Reynolds number or entrance.

use super::{
    Ambition, CARD, CONSTANT_PROPERTY_ASSUMPTION, CorrelationError, Evidence,
    HeatTransferCoefficient, ModelCard, ProvenanceHash, SensitivitySummary,
    SourceProvenance, Temperature, ValidityDomain, det,
};
use crate::DomainViolation;

/// Explicit fluid-dependent wall correction for a developing turbulent tube.
///
/// All base Re, Pr and conductivity queries still use bulk properties. The
/// liquid correction needs the SAME fluid's Prandtl number at the wall; it is
/// not a viscosity-ratio substitution. The air exponent must not be applied
/// to arbitrary gases. No caller-supplied exponent or empirical blend exists.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum TurbulentTubePropertyCorrection {
    /// Eq. (3): `K = (Pr_bulk / Pr_wall)^0.11`.
    Liquid { wall_prandtl: f64 },
    /// Eq. (4), air specialization: `K = (T_bulk / T_wall)^0.45`.
    Air {
        bulk_temperature: Temperature,
        wall_temperature: Temperature,
    },
}

impl TurbulentTubePropertyCorrection {
    /// Formula authority, separate from the uncorrected mean-card source.
    pub const fn source(self) -> SourceProvenance {
        SourceProvenance {
            citation: match self {
                Self::Liquid { .. } => "V. Gnielinski, On heat transfer in tubes, International Journal of Heat and Mass Transfer 63 (2013), 134-140, Eq. (3), liquid property correction",
                Self::Air { .. } => "V. Gnielinski, On heat transfer in tubes, International Journal of Heat and Mass Transfer 63 (2013), 134-140, Eq. (4), air exponent n=0.45",
            },
            identifier: "doi:10.1016/j.ijheatmasstransfer.2013.04.015",
        }
    }

    /// The correction's explicit admission and no-claim boundaries.
    ///
    /// The ratio boxes below are implementation policies, not validation
    /// ranges established in the 2013 paper. They never widen the original
    /// tube card's independent Re, Pr or length support. No numerical value
    /// for the corrected model discrepancy is inferred from the K=1 card.
    pub fn model_card(self) -> ModelCard {
        let (name, law, admission) = match self {
            Self::Liquid { .. } => (
                "convection.gnielinski-liquid-wall-properties",
                "Gnielinski (2013) Eq. (3): K=(Pr_bulk/Pr_wall)^0.11; wall Pr comes from the same liquid at its actual wall temperature",
                "Pr_bulk/Pr_wall in [0.05,20] is an implementation admission, not a source-published jointly validated range",
            ),
            Self::Air { .. } => (
                "convection.gnielinski-air-wall-properties",
                "Gnielinski (2013) Eq. (4): K=(T_bulk/T_wall)^0.45 for air; both temperatures are absolute kelvin",
                "T_bulk/T_wall in [0.5,1.5] is an implementation admission, not a source-published jointly validated range",
            ),
        };
        let (axis, low, high) = self.ratio_domain();
        ModelCard::new(
            name, "1.0.0", Ambition::Solid,
            [
                law, admission,
                "the retained baseline means use K=1; exactly one explicit wall-property factor multiplies their actual section integral",
                "Re, Pr and conductivity retain the original bulk-fluid state; a held correction does not resolve streamwise property or wall-temperature variation",
                "formula authority: doi:10.1016/j.ijheatmasstransfer.2013.04.015; no independent corrected-model discrepancy bound is assigned",
            ].into_iter().map(str::to_string).collect(),
            ValidityDomain::unconstrained().with(axis, low, high),
            [
                "substituting wall viscosity alone for wall Prandtl number, or using the air exponent for another gas",
                "transition, roughness, phase change, buoyancy or a property ratio outside the declared admission",
                "compounding corrections or claiming a validated local coefficient from a cumulative mean",
            ].into_iter().map(str::to_string).collect(),
            f64::INFINITY,
        )
    }

    const fn ratio_domain(self) -> (&'static str, f64, f64) {
        match self {
            Self::Liquid { .. } => ("Pr_over_Pr_wall", 0.05, 20.0),
            Self::Air { .. } => ("T_bulk_over_T_wall", 0.5, 1.5),
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub(super) struct PreparedCorrection {
    pub(super) input: TurbulentTubePropertyCorrection,
    bulk_prandtl: f64,
    ratio: f64,
    factor: f64,
}

impl PreparedCorrection {
    pub(super) fn new(
        input: TurbulentTubePropertyCorrection,
        bulk_prandtl: f64,
    ) -> Result<Self, CorrelationError> {
        let (ratio, exponent) = match input {
            TurbulentTubePropertyCorrection::Liquid { wall_prandtl } => {
                if !wall_prandtl.is_finite() || wall_prandtl <= 0.0 {
                    return Err(CorrelationError::InvalidGroup {
                        axis: "Pr_wall", value_bits: wall_prandtl.to_bits(),
                    });
                }
                (bulk_prandtl / wall_prandtl, 0.11)
            }
            TurbulentTubePropertyCorrection::Air { bulk_temperature, wall_temperature } => {
                let positive_kelvin = |field, temperature: Temperature| {
                    let value = temperature.value();
                    if value.is_finite() && value > 0.0 {
                        Ok(value)
                    } else {
                        Err(CorrelationError::InvalidDimensionalInput {
                            field, value_bits: value.to_bits(),
                            requirement: "supply a finite strictly positive absolute temperature in kelvin",
                        })
                    }
                };
                let bulk = positive_kelvin("bulk air temperature", bulk_temperature)?;
                let wall = positive_kelvin("air wall temperature", wall_temperature)?;
                (bulk / wall, 0.45)
            }
        };
        let (axis, low, high) = input.ratio_domain();
        if !ratio.is_finite() || ratio < low || ratio > high {
            return Err(CorrelationError::OutOfDomain {
                correlation: CARD,
                violations: vec![DomainViolation {
                    axis: axis.into(), value: Some(ratio), low, high,
                }],
            });
        }
        Ok(Self { input, bulk_prandtl, ratio, factor: det::pow(ratio, exponent) })
    }

    pub(super) const fn factor(self) -> f64 { self.factor }

    pub(super) const fn bulk_prandtl_elasticity(self) -> f64 {
        match self.input {
            TurbulentTubePropertyCorrection::Liquid { .. } => 0.11,
            TurbulentTubePropertyCorrection::Air { .. } => 0.0,
        }
    }

    /// Composition concerns the final coefficient; the original source means
    /// remain unmodified and separately inspectable through the section.
    pub(super) fn attach_evidence(self, coefficient: &mut Evidence<HeatTransferCoefficient>) {
        let card = self.input.model_card();
        coefficient.model.cards.push(card.name);
        coefficient.model.cards.sort_unstable();
        coefficient.model.cards.dedup();
        coefficient.model.assumptions.retain(|assumption| assumption != CONSTANT_PROPERTY_ASSUMPTION);
        coefficient.model.assumptions.extend(card.assumptions);
        coefficient.model.assumptions.sort_unstable();
        coefficient.model.assumptions.dedup();
        coefficient.model.validity = coefficient.model.validity.intersect(&card.validity);
        coefficient.model.discrepancy_rel = f64::INFINITY;
        // Hash actual coordinates, not just K: equal ratios at different
        // temperatures have different absolute temperature derivatives.
        let mut values = [0_u8; 25];
        let (tag, first, second) = match self.input {
            TurbulentTubePropertyCorrection::Liquid { wall_prandtl } => (0, wall_prandtl, 0.0),
            TurbulentTubePropertyCorrection::Air { bulk_temperature, wall_temperature } => {
                (1, bulk_temperature.value(), wall_temperature.value())
            }
        };
        values[0] = tag;
        for (bytes, value) in values[1..].chunks_exact_mut(8).zip([
            self.bulk_prandtl, first, second,
        ]) {
            bytes.copy_from_slice(&value.to_bits().to_le_bytes());
        }
        coefficient.provenance = ProvenanceHash::chain(
            "convection-gnielinski-wall-property-correction-v1",
            &[coefficient.provenance, ProvenanceHash::of_bytes(&values)],
        );
    }

    pub(super) fn attach_sensitivities(
        self,
        sensitivity: &mut SensitivitySummary,
        qoi: f64,
    ) -> Result<(), CorrelationError> {
        let (_, low, high) = self.input.ratio_domain();
        if self.ratio <= low || self.ratio >= high {
            if matches!(self.input, TurbulentTubePropertyCorrection::Liquid { .. }) {
                sensitivity.d_qoi.remove("Pr");
            }
            return Ok(());
        }
        let mut insert = |axis: &str, value: f64| {
            if value.is_finite() {
                sensitivity.d_qoi.insert(axis.into(), value);
                Ok(())
            } else {
                Err(CorrelationError::NonFiniteResult {
                    stage: "Gnielinski wall-property input derivative", value_bits: value.to_bits(),
                })
            }
        };
        match self.input {
            TurbulentTubePropertyCorrection::Liquid { wall_prandtl } => {
                insert("Pr_wall", -0.11 * qoi / wall_prandtl)?;
            }
            TurbulentTubePropertyCorrection::Air { bulk_temperature, wall_temperature } => {
                insert("T_bulk", 0.45 * qoi / bulk_temperature.value())?;
                insert("T_wall", -0.45 * qoi / wall_temperature.value())?;
            }
        }
        Ok(())
    }
}
