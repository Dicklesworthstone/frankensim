//! Gnielinski's finite-length smooth-tube relation and its integrated sections.
//!
//! V. Gnielinski, IJHMT 63 (2013), 134–140, Eqs. (2) and (6), with K=1.
//! The 2015 corrigendum, IJHMT 81, 638, restores the missing -2 exponent
//! in the friction-factor equations. This is the paper's 1.8 log10(Re)-1.5
//! convention, not the older friction fit of the separate `Gnielinski` card.

use super::{
    CorrelationError, CorrelationId, CorrelationInputs, CorrelationRobinBoundary, Evidence,
    HeatTransferCoefficient, Length, NumericalCertificate, NusseltEvaluation, ProvenanceHash,
    SensitivitySummary, Temperature, ThermalConductivity, attach_section_sensitivities, det, evaluate,
};
use crate::{
    Ambition, CorrelationCard, DiscrepancyBasis, ModelCard, SourceProvenance, ValidityDomain,
};

const CARD: CorrelationId = CorrelationId::GnielinskiDeveloping;

fn friction_root_inverse(reynolds: f64) -> f64 {
    1.8_f64.mul_add(det::ln(reynolds) / core::f64::consts::LN_10, -1.5)
}

fn developed_nusselt(reynolds: f64, prandtl: f64) -> f64 {
    let z = friction_root_inverse(reynolds);
    let friction = 1.0 / (z * z);
    let denominator = 1.0 + 12.7 * det::sqrt(friction / 8.0) * (det::pow(prandtl, 2.0 / 3.0) - 1.0);
    (friction / 8.0) * (reynolds - 1_000.0) * prandtl / denominator
}

pub(crate) fn mean_nusselt(reynolds: f64, prandtl: f64, length_over_diameter: f64) -> f64 {
    developed_nusselt(reynolds, prandtl) * (1.0 + det::pow(length_over_diameter, -2.0 / 3.0))
}

// The geometry multiplier is independent of Re and Pr, for both inlet means
// and finite sections. Compute these two transport partials together.
fn transport_elasticities(reynolds: f64, prandtl: f64) -> [f64; 2] {
    let z = friction_root_inverse(reynolds);
    let coefficient = 12.7 / (det::sqrt(8.0) * z);
    let pr_power = det::pow(prandtl, 2.0 / 3.0);
    let correction = coefficient * (pr_power - 1.0);
    [
        reynolds / (reynolds - 1_000.0)
            - (3.6 / core::f64::consts::LN_10) / z
                * (1.0 - 0.5 * correction / (1.0 + correction)),
        1.0 - (2.0 / 3.0) * coefficient * pr_power / (1.0 + correction),
    ]
}

/// Logarithmic partials of the same mean formula at admitted interior points.
pub(crate) fn elasticity(
    reynolds: f64,
    prandtl: f64,
    length_over_diameter: f64,
    axis: &str,
) -> f64 {
    match axis {
        "Re" => transport_elasticities(reynolds, prandtl)[0],
        "Pr" => transport_elasticities(reynolds, prandtl)[1],
        "L_over_Dh" => {
            let entrance = det::pow(length_over_diameter, -2.0 / 3.0);
            -(2.0 / 3.0) * entrance / (1.0 + entrance)
        }
        _ => 0.0,
    }
}

/// Positive endpoint interval of the same primitive, without subtraction.
fn section_nusselt(reynolds: f64, prandtl: f64, start: f64, end: f64) -> f64 {
    let a = det::pow(start, 1.0 / 3.0);
    let b = det::pow(end, 1.0 / 3.0);
    developed_nusselt(reynolds, prandtl) * (1.0 + 1.0 / (a * a + a * b + b * b))
}

pub(crate) fn card() -> CorrelationCard {
    CorrelationCard {
        id: CARD,
        source: SourceProvenance {
            citation: "V. Gnielinski, On heat transfer in tubes, International Journal of Heat and Mass Transfer 63 (2013), 134-140, Eqs. (2) and (6), K=1; corrected friction exponent in Corrigendum, IJHMT 81 (2015), 638",
            identifier: "doi:10.1016/j.ijheatmasstransfer.2013.04.015; doi:10.1016/j.ijheatmasstransfer.2014.10.063",
        },
        model: ModelCard::new(
            CARD.name(),
            "1.0.0",
            Ambition::Solid,
            [
                "mean forced heat transfer from the shared physical heated entrance of a smooth straight circular tube",
                "single-phase Newtonian turbulent flow; no roughness or buoyancy correction",
                "bulk transport properties are held constant and the paper's property correction K is exactly one; no wall-viscosity or gas-temperature correction is implied",
                "the length multiplier includes developing hydraulic and thermal boundary layers; L is heated distance from the inlet, never mesh-cell width",
                "Darcy friction factor is [1.8 log10(Re)-1.5]^-2 from corrected Eq. (6); it is not the older fully-developed card's friction fit",
                "Re in [1e4,1e6], Pr in [0.5,200], L/D in [10,1e6] are conservative implementation admissions, not a source-published jointly validated range",
                "the admission excludes the paper's low-Re transition and very-short-tube caveats; inclusion in this numeric box does not establish physical applicability",
            ].into_iter().map(str::to_string).collect(),
            ValidityDomain::unconstrained()
                .with("Re", 1.0e4, 1.0e6)
                .with("Pr", 0.5, 200.0)
                .with("L_over_Dh", 10.0, 1.0e6),
            [
                "transition, rough tubes, non-circular passages, swirl, or strong secondary flow",
                "significant property changes requiring the paper's liquid or gas correction K; re-querying local bulk properties is only a held-property approximation",
                "L/D below ten, a pre-existing unheated hydrodynamic entrance, or a different inlet geometry not represented by the source model",
                "phase change, compressibility, strong buoyancy, non-Newtonian fluids, or physical validation inferred from formula tests",
                "using the mean coefficient independently on every downstream cell, which restarts the thermal entrance and overcounts transfer",
            ].into_iter().map(str::to_string).collect(),
            0.20,
        ),
        discrepancy_basis: DiscrepancyBasis::EngineeringAllowance,
    }
}

/// Mean transfer over one section of a developing smooth turbulent circular tube.
/// Distances use the same physical inlet and diameter as the retained mean card.
#[derive(Debug)]
pub struct TurbulentCircularDuctSection {
    start_over_diameter: f64,
    end_over_diameter: f64,
    upstream: Option<NusseltEvaluation>,
    downstream: NusseltEvaluation,
    nusselt: f64,
}

impl TurbulentCircularDuctSection {
    pub const fn start_over_diameter(&self) -> f64 {
        self.start_over_diameter
    }
    pub const fn end_over_diameter(&self) -> f64 {
        self.end_over_diameter
    }
    /// Original upstream mean, absent only at the exact physical entrance.
    pub const fn upstream(&self) -> Option<&NusseltEvaluation> {
        self.upstream.as_ref()
    }
    pub const fn downstream(&self) -> &NusseltEvaluation {
        &self.downstream
    }
    /// Section-average Nusselt number based on the tube diameter.
    pub const fn nusselt(&self) -> f64 {
        self.nusselt
    }

    /// Move out both physical endpoint ratios, original cumulative evaluations
    /// and the derived section value, in that order. The downstream evaluation
    /// remains an original mean; it is never relabeled as the section value.
    /// Compute the coefficient first to retain the section's complete evidence.
    #[must_use]
    pub fn into_parts(self) -> (f64, f64, Option<NusseltEvaluation>, NusseltEvaluation, f64) {
        (
            self.start_over_diameter,
            self.end_over_diameter,
            self.upstream,
            self.downstream,
            self.nusselt,
        )
    }

    /// Convert using the same bulk-fluid conductivity and tube diameter.
    ///
    /// A section beginning at zero retains the original mean evidence exactly.
    /// A positive-start section retains both mean evaluations, but its model
    /// discrepancy is unbounded: the original mean allowance does not validate
    /// a derived local heat-transfer profile. This primal-only conversion
    /// does not attach input sensitivities.
    pub fn heat_transfer_coefficient(
        &self,
        conductivity: ThermalConductivity,
        diameter: Length,
    ) -> Result<Evidence<HeatTransferCoefficient>, CorrelationError> {
        let mut evidence = self
            .downstream
            .heat_transfer_coefficient(conductivity, diameter)?;
        let Some(upstream) = &self.upstream else {
            return Ok(evidence);
        };
        let value = self.nusselt * conductivity.value() / diameter.value();
        if !value.is_finite() || value <= 0.0 {
            return Err(CorrelationError::NonFiniteResult {
                stage: "turbulent section Nu-to-h conversion",
                value_bits: value.to_bits(),
            });
        }
        evidence.value = HeatTransferCoefficient::new(value);
        evidence.qoi = value;
        evidence.numerical = NumericalCertificate::estimate(value, value);
        evidence.model.discrepancy_rel = f64::INFINITY;
        evidence.model.assumptions.push(
            "section integral of the held-property Gnielinski (2013) developing-tube mean; original mean allowance does not bound section discrepancy".into(),
        );
        evidence.model.assumptions.sort_unstable();
        evidence.model.assumptions.dedup();
        evidence.sensitivity = SensitivitySummary::default();
        let mut values = [0_u8; 32];
        for (bytes, input) in values.chunks_exact_mut(8).zip([
            self.start_over_diameter,
            self.end_over_diameter,
            conductivity.value(),
            diameter.value(),
        ]) {
            bytes.copy_from_slice(&input.to_bits().to_le_bytes());
        }
        evidence.provenance = ProvenanceHash::chain(
            "convection-turbulent-section-nu-to-h-v1",
            &[
                upstream.evidence().provenance,
                self.downstream.evidence().provenance,
                ProvenanceHash::of_bytes(&values),
            ],
        );
        Ok(evidence)
    }

    /// Convert this section's actual coefficient with analytic `Re`, `Pr`, `k`
    /// and `L` partials in `sensitivity.d_qoi`, using the same absolute
    /// convention as [`NusseltEvaluation::heat_transfer_coefficient_with_sensitivities`].
    /// `L` is the diameter. The two inlet-distance/diameter ratios and other
    /// independent inputs stay fixed; no endpoint or total geometry derivative
    /// is implied. Re/Pr partials scale the actual section coefficient, whose
    /// geometry factor differs from the downstream cumulative mean.
    ///
    /// A partial is absent at its source-domain endpoint; an unrelated fixed
    /// length boundary does not remove it. Derivatives describe the admitted
    /// held-property formula and do not validate a local heat-transfer profile
    /// or promote its unbounded positive-start model discrepancy. Invalid
    /// dimensions and nonfinite derivative arithmetic are refused.
    pub fn heat_transfer_coefficient_with_sensitivities(
        &self,
        conductivity: ThermalConductivity,
        diameter: Length,
    ) -> Result<Evidence<HeatTransferCoefficient>, CorrelationError> {
        let mut coefficient = self.heat_transfer_coefficient(conductivity, diameter)?;
        let elasticities = transport_elasticities(
            self.downstream.groups()["Re"],
            self.downstream.groups()["Pr"],
        );
        attach_section_sensitivities(
            &mut coefficient,
            self.nusselt,
            conductivity,
            diameter,
            (self.upstream.as_ref(), &self.downstream),
            elasticities,
        )?;
        Ok(coefficient)
    }

    /// Lower this section coefficient into the existing conduction Robin seam.
    pub fn robin_boundary(
        &self,
        conductivity: ThermalConductivity,
        diameter: Length,
        reference_temperature: Temperature,
    ) -> Result<CorrelationRobinBoundary, CorrelationError> {
        let temperature = reference_temperature.value();
        if !temperature.is_finite() {
            return Err(CorrelationError::InvalidDimensionalInput {
                field: "Robin reference temperature",
                value_bits: temperature.to_bits(),
                requirement: "supply a finite coherent-SI absolute temperature",
            });
        }
        let coefficient = self.heat_transfer_coefficient(conductivity, diameter)?;
        let boundary = fs_conduction::ThermalBc::robin(coefficient.value.value(), temperature)
            .map_err(|_| CorrelationError::NonFiniteResult {
                stage: "conduction Robin boundary admission",
                value_bits: coefficient.value.value().to_bits(),
            })?;
        Ok(CorrelationRobinBoundary {
            coefficient,
            reference_temperature,
            boundary,
        })
    }
}

/// Integrate the developing turbulent mean between distances from one inlet.
///
/// For `s=x/D`, cumulative dimensionless transfer is
/// `s Nu_mean(s) = Nu_infinity (s + s^(1/3))`. Its exact zero-length limit is
/// zero, so a first section returns the original downstream mean without
/// evaluating an inadmissible zero-length card. Positive starts and ends must
/// satisfy the unchanged original card, including `s >= 10`.
///
/// For positive starts `a,b`, rationalizing the cube-root difference gives
/// `Nu_section = Nu_infinity [1 + 1/(a^(2/3)+(ab)^(1/3)+b^(2/3))]`.
/// This is stable even for adjacent floating-point endpoints and integrates
/// to the original cumulative transfer without restarting the inlet per cell.
/// Additivity requires the same Re, Pr, diameter, bulk properties and wall
/// boundary convention. It does not establish a resolved variable-property
/// boundary layer, local experimental validity, or a section adjoint.
pub fn evaluate_turbulent_circular_duct_section(
    reynolds: f64,
    prandtl: f64,
    start_over_diameter: f64,
    end_over_diameter: f64,
) -> Result<TurbulentCircularDuctSection, CorrelationError> {
    if !start_over_diameter.is_finite() || start_over_diameter < 0.0 {
        return Err(CorrelationError::InvalidGroup {
            axis: "section_start_over_Dh",
            value_bits: start_over_diameter.to_bits(),
        });
    }
    if !end_over_diameter.is_finite() || end_over_diameter <= start_over_diameter {
        return Err(CorrelationError::InvalidGroup {
            axis: "section_end_over_Dh",
            value_bits: end_over_diameter.to_bits(),
        });
    }
    let original = |length| {
        evaluate(
            CARD,
            CorrelationInputs::forced(reynolds, prandtl).with_length_ratio(length),
        )
    };
    let downstream = original(end_over_diameter)?;
    let upstream = if start_over_diameter == 0.0 {
        None
    } else {
        Some(original(start_over_diameter)?)
    };
    let nusselt = if upstream.is_none() {
        downstream.evidence().value
    } else {
        section_nusselt(reynolds, prandtl, start_over_diameter, end_over_diameter)
    };
    if !nusselt.is_finite() || nusselt <= 0.0 {
        return Err(CorrelationError::NonFiniteResult {
            stage: "turbulent section cumulative-transfer difference",
            value_bits: nusselt.to_bits(),
        });
    }
    Ok(TurbulentCircularDuctSection {
        start_over_diameter,
        end_over_diameter,
        upstream,
        downstream,
        nusselt,
    })
}
