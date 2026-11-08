//! Ideal-drive feedback along an original caloric T/f path.

use super::*;
use crate::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};
use fs_matdb::{MaterialAnswer, QueryPoint};

impl MaterialOhmicNode {
    /// Bound positive `dP_Joule/dh` over this complete original caloric chart
    /// [W/(J/kg)]. Divide by specimen mass to obtain the thermal `dP/dE` bound.
    ///
    /// Each original caloric interval follows the pinned material source path,
    /// retaining every crossed electrical knot/cell edge and checking support
    /// and admissibility of the COMPLETE source bundle. Native conductors own
    /// every drive and power evaluation, including zero/nonpositive feedback.
    ///
    /// Authored T(h)/f(h) intervals are affine. Sourced sensible intervals keep
    /// their original affine Cp(T); its native minimum bounds the path speed,
    /// rather than replacing nonlinear storage by the h/T endpoint chord.
    /// Latent intervals may have positive electrical feedback at constant T.
    /// Without an owned fraction axis, only temperature changes the source.
    ///
    /// The caller may cache this bound for the SAME curve identity and drive.
    /// It bounds the declared source model with outward rounding, not material
    /// uncertainty. Geometry is fixed; free expansion requires another bound.
    pub fn maximum_positive_joule_specific_enthalpy_slope(
        &self,
        drive: OhmicDrive,
        curve: &EquilibriumEnthalpyPhaseCurve,
    ) -> Result<f64, ConductorError> {
        let original = self.response.original();
        if self.free_expansion.is_some()
            || curve.material_card_identity() != original.material().card_identity()
        {
            return Err(invalid(
                "Joule enthalpy bound requires fixed geometry and the original material card",
            ));
        }
        let inverse = matches!(
            (original.source_basis(), drive),
            (ElectricalPropertyBasis::Resistivity, OhmicDrive::Voltage(_))
                | (
                    ElectricalPropertyBasis::Conductivity,
                    OhmicDrive::Current(_)
                )
        );
        let query_at = |knot: &EnthalpyPhaseKnot| -> Result<QueryPoint, ConductorError> {
            if !(self.temperature_range_k[0]..=self.temperature_range_k[1])
                .contains(&knot.temperature_k)
            {
                return Err(invalid(
                    "caloric chart leaves the Ohmic temperature interval",
                ));
            }
            let mut query = original
                .material()
                .query()
                .with(&self.temperature_axis, knot.temperature_k)
                .map_err(|_| invalid("invalid caloric Ohmic temperature coordinate"))?;
            if let Some(phase) = &self.phase {
                phase.admit_fraction(knot.liquid_mass_fraction)?;
                query = query
                    .with(&phase.axis, knot.liquid_mass_fraction)
                    .map_err(|_| invalid("invalid caloric Ohmic fraction coordinate"))?;
            }
            Ok(query)
        };
        let power_at = |answer: &MaterialAnswer| -> Result<f64, ConductorError> {
            // Retain the source owner's actual coordinates at each crossed knot,
            // rather than reconstructing them from a rounded path parameter.
            let coordinate = |axis: &str| {
                answer
                    .receipt
                    .query_point
                    .iter()
                    .find(|(name, _)| name == axis)
                    .map(|(_, value)| *value)
                    .ok_or_else(|| invalid("Ohmic source path lost an original coordinate"))
            };
            let temperature = coordinate(&self.temperature_axis)?;
            let conductor = if let Some(phase) = &self.phase {
                self.resolve_at_temperature_and_liquid_mass_fraction(
                    temperature,
                    coordinate(&phase.axis)?,
                )?
            } else {
                self.resolve_at_temperature_k(temperature)?
            };
            conductor.joule_power_w(conductor.current_a_for_drive(drive)?)
        };
        let mut bound: f64 = 0.0;
        for (index, pair) in curve.knots().windows(2).enumerate() {
            let start = query_at(&pair[0])?;
            let end = query_at(&pair[1])?;
            let mut electrical = None;
            for property in original.material().properties() {
                let selected = property.answer().receipt.selected;
                let claim = self
                    .response
                    .card()
                    .claims()
                    .claim(selected)
                    .ok_or_else(|| invalid("original Ohmic source claim is absent"))?;
                let path = self
                    .response
                    .card()
                    .claims()
                    .query_path_pinned_typed(&claim.key, &start, &end, selected)
                    .map_err(|source| MaterialStatePointError::Support {
                        property: claim.key.name().to_owned(),
                        source: Box::new(source),
                    })?;
                let requirement = property.requirement().admissibility();
                for (_, answer) in &path.knots {
                    let value = answer.evidence.value.value;
                    let admitted = value.is_finite()
                        && match requirement {
                            ScalarAdmissibility::Finite => true,
                            ScalarAdmissibility::StrictlyPositive => value > 0.0,
                            ScalarAdmissibility::NonNegative => value >= 0.0,
                            ScalarAdmissibility::OpenInterval { lower, upper } => {
                                lower < value && value < upper
                            }
                        };
                    if !admitted {
                        return Err(MaterialStatePointError::OutsideConsumerDomain {
                            property: claim.key.name().to_owned(),
                            value,
                            admissibility: requirement,
                        }
                        .into());
                    }
                }
                if &claim.key == original.source_key() {
                    electrical = Some(path);
                }
            }
            let path =
                electrical.ok_or_else(|| invalid("Ohmic electrical source path is absent"))?;
            let sensible = curve
                .sensible_heat_capacities_on_interval(index)
                .map(|_| {
                    curve
                        .sensible_heat_capacity_bounds_j_kg_k(index)
                        .map(|cp| (pair[1].temperature_k - pair[0].temperature_k, cp[0]))
                        .map_err(|_| invalid("original sensible heat capacity bound is invalid"))
                })
                .transpose()?;
            let mut left_power = power_at(&path.knots[0].1)?;
            for interval in path.knots.windows(2) {
                let right_power = power_at(&interval[1].1)?;
                let width = feedback::minimum_enthalpy_path_width(
                    interval[1].0 - interval[0].0,
                    pair[1].specific_enthalpy_j_kg - pair[0].specific_enthalpy_j_kg,
                    sensible,
                )
                .ok_or(ConductorError::InvalidDerived {
                    quantity: "Ohmic source path enthalpy width",
                })?;
                let slope = feedback::positive_segment_slope_w_k(
                    width,
                    [
                        interval[0].1.evidence.value.value,
                        interval[1].1.evidence.value.value,
                    ],
                    [left_power, right_power],
                    inverse,
                )
                .ok_or(ConductorError::InvalidDerived {
                    quantity: "positive Joule specific enthalpy slope",
                })?;
                bound = bound.max(slope);
                left_power = right_power;
            }
        }
        Ok(bound)
    }
}
