//! Two-film resistance on the same finite-inventory Henry trajectory.
//!
//! Eliminate the interface from J = (D_L/delta_L)(C_i-C_L),
//! J = (D_G/delta_G)(p_G-p_i)/(RT), and C_i = H_s^cp p_i.
//! Thus 1/K_L = delta_L/D_L + H_s^cp RT delta_G/D_G. The dimensionless
//! conversion is H_s^cc = RT H_s^cp, not its reciprocal (IUPAC 2021,
//! https://doi.org/10.1515/pac-2020-0302, section 4.2).
//!
//! This is dilute Fick diffusion, with negligible film storage and no Stefan
//! flow. The gas diffusivity source must describe the selected solute in the
//! declared carrier and remain constant over the accepted pressure/composition
//! domain. A measured liquid diffusivity is never reused as a gas diffusivity.

use super::*;
use fs_qty::semantic::CompositionBasis;

/// Selected solute diffusivity in its declared carrier gas [m²/s].
pub const GAS_FILM_DIFFUSIVITY_PROPERTY: &str = "gas_phase_solute_diffusivity";
/// Canonical semantic mole-fraction axis of the selected gas-film solute.
pub const GAS_FILM_MOLE_FRACTION_AXIS: &str = "gas_mole_fraction";

fn mole_fraction_type() -> QuantitySpec {
    typed(QuantityKind::Composition(CompositionBasis::MoleFraction))
}

/// Explicit gas-phase solute diffusivity in the authored carrier mixture.
/// The original claim must qualify the canonical, semantic mole-fraction axis;
/// runtime queries update it at both gas bulk and gas-side interface.
#[derive(Clone, Debug, PartialEq)]
pub struct GasFilmSource {
    /// Solute identity, which must match the original Henry law's species.
    pub species: SpeciesId,
    /// Immutable gas-phase diffusivity evidence for this solute and carrier.
    pub material: MaterialCard,
    /// Typed temperature, pressure and mole fraction plus fixed source context.
    pub point: QueryPoint,
    /// Explicit original-claim selection; subsequent queries retain that claim.
    pub selection: MaterialPropertySelection,
    /// Authored dilute-Fick limit in (0,1), not a validation of that closure.
    /// For simultaneous exchange it bounds the SUM of exchanging solutes
    /// across all gas-film locations and times, including unequal films.
    pub maximum_mole_fraction: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub(super) struct ResolvedGasFilm {
    pub(super) source: GasFilmSource,
    resolved: ResolvedMaterialStatePoint,
    pub(super) thickness_m: f64,
    pub(super) diffusivity_m2_s: f64,
}

impl ResolvedGasFilm {
    fn point(
        &self,
        temperature_k: f64,
        total_pressure_pa: f64,
        partial_pressure_pa: f64,
    ) -> Result<QueryPoint> {
        let fraction = finite(partial_pressure_pa / total_pressure_pa)?;
        if partial_pressure_pa > 0.0 && fraction == 0.0 {
            return Err(invalid("positive gas-film mole fraction underflows"));
        }
        if !positive(total_pressure_pa)
            || !nonnegative(fraction)
            || fraction > self.source.maximum_mole_fraction
        {
            return Err(DissolvedGasError::GasFilmConcentration {
                actual_mole_fraction: fraction,
                maximum_mole_fraction: self.source.maximum_mole_fraction,
            });
        }
        self.source
            .point
            .clone()
            .with_quantity(
                "temperature",
                typed(QuantityKind::AbsoluteTemperature),
                temperature_k,
            )
            .and_then(|point| {
                point.with_quantity("pressure", typed(QuantityKind::Pressure), total_pressure_pa)
            })
            .and_then(|point| {
                point.with_quantity(GAS_FILM_MOLE_FRACTION_AXIS, mole_fraction_type(), fraction)
            })
            .map_err(|_| invalid("unrepresentable gas-film source conditions"))
    }

    /// Continue the same selected diffusivity claim before preparing a new
    /// temperature-dependent interval. Bulk/interface sampling below still
    /// requires the refreshed value to be constant at that fixed temperature.
    pub(super) fn at_conditions(
        &self,
        temperature_k: f64,
        total_pressure_pa: f64,
        partial_pressure_pa: f64,
    ) -> Result<Self> {
        let point = self.point(temperature_k, total_pressure_pa, partial_pressure_pa)?;
        let resolved = self
            .resolved
            .at_query_point(&self.source.material, &point)?;
        let mut next = self.clone();
        next.diffusivity_m2_s = resolved
            .property(GAS_FILM_DIFFUSIVITY_PROPERTY)
            .unwrap()
            .value_si();
        next.resolved = resolved;
        next.source.point = point;
        Ok(next)
    }

    pub(super) fn sample(
        &self,
        total_pressure_pa: f64,
        partial_pressure_pa: f64,
    ) -> Result<ResolvedMaterialStatePoint> {
        let point = self.point(
            self.source.point.axes()["temperature"],
            total_pressure_pa,
            partial_pressure_pa,
        )?;
        let at = self
            .resolved
            .at_query_point(&self.source.material, &point)?;
        if at
            .property(GAS_FILM_DIFFUSIVITY_PROPERTY)
            .unwrap()
            .value_si()
            != self.diffusivity_m2_s
        {
            return Err(invalid(
                "gas-film diffusivity changed during a fixed-coefficient step",
            ));
        }
        Ok(at)
    }
}

/// Prepared once per interval; both endpoints use these same film weights.
#[derive(Clone, Copy)]
pub(super) struct FilmResponse {
    pub(super) coefficient_m_s: f64,
    pub(super) conductance_m3_s: f64,
    pub(super) liquid_weight: f64,
    pub(super) gas_weight: f64,
    pub(super) liquid_relaxation_ratio: f64,
    pub(super) gas_relaxation_ratio: Option<f64>,
}

impl HenryDissolution {
    /// Add finite gas-side resistance to the original finite-pair law.
    /// The gas and liquid films share the area supplied to `advance`; their
    /// physical thicknesses and original diffusivity sources stay distinct.
    /// Both phases remain at the original exact temperature. Each film's
    /// diffusion/bulk timescale ratio must satisfy the existing step controls.
    /// `DissolutionFilm::gas_side_transfer` must remain None: it is an
    /// alternative description of this same gas resistance.
    ///
    /// Zero diffusivity seals a film. Two sealed films refuse because neither
    /// fixes the common interface trace needed for source admission.
    pub fn with_gas_film(mut self, source: GasFilmSource, thickness_m: f64) -> Result<Self> {
        if source.species != self.species
            || !positive(thickness_m)
            || !positive(source.maximum_mole_fraction)
            || source.maximum_mole_fraction >= 1.0
            || source.point.axes().get("temperature") != Some(&self.temperature_k)
            || source.point.axis_quantities().get("temperature")
                != Some(&typed(QuantityKind::AbsoluteTemperature))
            || source.point.axis_quantities().get("pressure")
                != Some(&typed(QuantityKind::Pressure))
            || !source
                .point
                .axes()
                .get("pressure")
                .is_some_and(|p| positive(*p))
            || source
                .point
                .axis_quantities()
                .get(GAS_FILM_MOLE_FRACTION_AXIS)
                != Some(&mole_fraction_type())
        {
            return Err(invalid(
                "gas film needs the same species/temperature, positive thickness, dilute limit and canonical typed pressure/mole fraction",
            ));
        }
        let resolved = resolve_material_state_point(
            &source.material,
            &source.point,
            &[ScalarPropertyRequirement::try_new(
                GAS_FILM_DIFFUSIVITY_PROPERTY,
                DIFFUSIVITY_DIMS,
                ScalarAdmissibility::NonNegative,
            )?],
            source.selection.clone(),
        )?;
        admit_frozen_sources(&resolved)?;
        let property = resolved.property(GAS_FILM_DIFFUSIVITY_PROPERTY).unwrap();
        let validity = &property.answer().evidence.model.validity;
        if validity.axis_quantities().get(GAS_FILM_MOLE_FRACTION_AXIS)
            != Some(&mole_fraction_type())
            || !validity
                .bound(GAS_FILM_MOLE_FRACTION_AXIS)
                .is_some_and(|(lo, hi)| lo >= 0.0 && hi <= 1.0)
            || validity.axis_quantities().iter().any(|(axis, quantity)| {
                axis != GAS_FILM_MOLE_FRACTION_AXIS
                    && quantity.semantic_type().is_some_and(|semantic| {
                        matches!(semantic.kind(), QuantityKind::Composition(_))
                    })
            })
        {
            return Err(invalid(
                "gas diffusivity claim must qualify its own semantic mole-fraction domain",
            ));
        }
        let diffusivity = property.value_si();
        if self.diffusivity_m2_s == 0.0 && diffusivity == 0.0 {
            return Err(invalid(
                "two sealed films leave the common interface trace undetermined",
            ));
        }
        let film = ResolvedGasFilm {
            source,
            resolved,
            thickness_m,
            diffusivity_m2_s: diffusivity,
        };
        let fraction = film.source.point.axes()[GAS_FILM_MOLE_FRACTION_AXIS];
        if !nonnegative(fraction) || fraction > film.source.maximum_mole_fraction {
            return Err(DissolvedGasError::GasFilmConcentration {
                actual_mole_fraction: fraction,
                maximum_mole_fraction: film.source.maximum_mole_fraction,
            });
        }
        self.gas_film = Some(film);
        Ok(self)
    }

    /// Retained optional gas-film evidence and authored dilution limit.
    pub fn gas_film_source(&self) -> Option<&GasFilmSource> {
        self.gas_film.as_ref().map(|film| &film.source)
    }

    /// Physical gas-film thickness [m], absent for negligible gas resistance.
    pub fn gas_film_thickness_m(&self) -> Option<f64> {
        self.gas_film.as_ref().map(|film| film.thickness_m)
    }

    /// Original resolved gas-phase solute diffusivity [m²/s] at fixed temperature.
    pub fn gas_film_diffusivity_m2_s(&self) -> Option<f64> {
        self.gas_film.as_ref().map(|film| film.diffusivity_m2_s)
    }

    pub(super) fn film_response(
        &self,
        liquid: DissolutionFilm,
        capacity: f64,
    ) -> Result<FilmResponse> {
        let liquid_coefficient = finite(self.diffusivity_m2_s / liquid.thickness_m)?;
        // Both public descriptions supply the SAME gas coefficient on the
        // liquid-concentration basis. The original source alone also supplies
        // a physical gas diffusion time; a bare pressure coefficient does not.
        let (gas_coefficient, gas_diffusion_geometry) =
            match (&self.gas_film, liquid.gas_side_transfer) {
                (Some(_), Some(_)) => {
                    return Err(invalid(
                        "choose either a gas-film source or a gas-side coefficient, not both",
                    ));
                }
                (Some(gas), None) => {
                    let partition = finite(self.henry_mol_m3_pa * R * self.temperature_k)?;
                    if !positive(partition) {
                        return Err(invalid(
                            "positive dimensionless gas/liquid partition underflows",
                        ));
                    }
                    let coefficient = finite((gas.diffusivity_m2_s / gas.thickness_m) / partition)?;
                    if gas.diffusivity_m2_s > 0.0 && !positive(coefficient) {
                        return Err(invalid("positive film coefficient underflows"));
                    }
                    (Some(coefficient), Some((gas.thickness_m, partition)))
                }
                (None, Some(gas)) => {
                    let coefficient = finite(gas.value_si() / self.henry_mol_m3_pa)?;
                    if gas.value_si() > 0.0 && !positive(coefficient) {
                        return Err(invalid(
                            "gas-film liquid-equivalent coefficient is unrepresentable",
                        ));
                    }
                    (Some(coefficient), None)
                }
                (None, None) => (None, None),
            };
        let (coefficient, conductance, liquid_weight, gas_weight) =
            if let Some(gas_coefficient) = gas_coefficient {
                if self.diffusivity_m2_s > 0.0 && !positive(liquid_coefficient) {
                    return Err(invalid("positive film coefficient underflows"));
                }
                if liquid_coefficient == 0.0 && gas_coefficient == 0.0 {
                    return Err(invalid(
                        "two sealed films have no determinate interfacial state",
                    ));
                }
                // Harmonic conductance and complementary weights without
                // overflowing a resistance or losing the smaller weight in 1-w.
                let (coefficient, liquid_weight, gas_weight) =
                    if liquid_coefficient <= gas_coefficient {
                        let ratio = liquid_coefficient / gas_coefficient;
                        let liquid_weight = 1.0 / (1.0 + ratio);
                        (
                            liquid_coefficient * liquid_weight,
                            liquid_weight,
                            ratio * liquid_weight,
                        )
                    } else {
                        let ratio = gas_coefficient / liquid_coefficient;
                        let gas_weight = 1.0 / (1.0 + ratio);
                        (gas_coefficient * gas_weight, ratio * gas_weight, gas_weight)
                    };
                if liquid_coefficient > 0.0
                    && gas_coefficient > 0.0
                    && (liquid_weight == 0.0 || gas_weight == 0.0)
                {
                    return Err(invalid("positive two-film interface weight underflows"));
                }
                if liquid_coefficient > 0.0 && gas_coefficient > 0.0 && !positive(coefficient) {
                    return Err(invalid("positive two-film transfer is unrepresentable"));
                }
                (
                    coefficient,
                    finite(coefficient * liquid.area_m2)?,
                    liquid_weight,
                    gas_weight,
                )
            } else {
                // Preserve the original liquid-only arithmetic and response.
                (
                    liquid_coefficient,
                    finite(self.diffusivity_m2_s * liquid.area_m2 / liquid.thickness_m)?,
                    1.0,
                    0.0,
                )
            };
        let conducting = self.diffusivity_m2_s > 0.0
            && gas_coefficient.is_none_or(|coefficient| coefficient > 0.0);
        if conducting && liquid.area_m2 > 0.0 && !positive(conductance) {
            return Err(invalid("positive film conductance underflows"));
        }
        let gas_ratio = gas_diffusion_geometry
            .map(|(thickness, partition)| {
                finite(liquid.area_m2 * thickness * capacity / partition * gas_weight)
            })
            .transpose()?;
        Ok(FilmResponse {
            coefficient_m_s: coefficient,
            conductance_m3_s: conductance,
            liquid_weight,
            gas_weight,
            liquid_relaxation_ratio: finite(
                liquid.area_m2 * liquid.thickness_m * capacity * liquid_weight,
            )?,
            gas_relaxation_ratio: gas_ratio,
        })
    }
}
