//! One homogeneous freely expanding conductor on the original heat timeline.
//!
//! The chosen closure is isotropic free expansion at an explicitly maintained
//! absolute source pressure. Current volume is m/rho(T), all lengths scale from
//! their ORIGINAL dimensions by cbrt(V/V0), and the same current section supplies
//! electrical resistance and mechanical coefficients. No independent expansion
//! coefficient is also applied. Anisotropy, constrained expansion, multilayer
//! compatibility, expansion inertia and strain-dependent material laws are not
//! inferred. The caller must justify a uniform temperature and quasistatic motion.
//!
//! Thermal cells retain reference geometry and invariant mass. The native ambient
//! conductance hook uses the current exposed area for BOTH heat flow and the CFL;
//! no parallel heat law or explicit integrator is introduced.

use super::{
    Boundary, CellId, ConcentricElasticSection, ConcentricSectionTemperatureLaw, DepositionBudget,
    DepositionSimulation, HeatedConcentricSection, ResolvedHeatCapacityCurve, SectionHeatingError,
    positive,
};
use crate::phase::isobaric::IsobaricHeatStep;
use crate::section::ConcentricMaterialLayer;
use crate::state_point::SPECIFIC_HEAT_CAPACITY_PROPERTY;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Pressure, QuantitySpec, Temperature};

/// Actual dimensions of the accepted, fixed-mass homogeneous specimen.
/// The complete original/current material receipts remain in section().
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SectionDensityGeometry {
    pub mass_kg: f64,
    pub density_kg_m3: f64,
    /// Relative to the original length and both original radii, never accumulated.
    pub isotropic_scale: f64,
    pub length_m: f64,
    pub volume_m3: f64,
    /// Entire cylindrical/tubular surface including the two end annuli.
    pub surface_area_m2: f64,
}

#[derive(Clone)]
pub(super) struct DensityGeometry<'a> {
    caloric: &'a ResolvedHeatCapacityCurve,
    reference_length_m: f64,
    reference_volume_m3: f64,
    reference_surface_m2: f64,
    pressure_pa: f64,
    pub(super) current: SectionDensityGeometry,
}

impl<'a> HeatedConcentricSection<'a> {
    /// Couple temperature, original density, free dimensions, resistance and heat.
    ///
    /// Exactly one homogeneous solid layer is admitted, optionally with a bore.
    /// The original Cp and density must explicitly support the named, typed
    /// pressure coordinate. That value is declared ABSOLUTE pressure, maintained
    /// throughout the trajectory. The original caloric chart owns h(T) and rho(T);
    /// its selected claims and all fixed conditions must match the section.
    /// Density must describe the relaxed material on that isobar. Every retained
    /// property must use the owned temperature/pressure and static external
    /// conditions; unowned stress, strain, rate and field axes are refused.
    ///
    /// The single cell keeps its reference inventory. Current density changes
    /// volume, length, radii, resistance and ambient convection area at every
    /// actual thermal substep. Ideal terminal drive and authored total heat are
    /// retained; ambient htc is per CURRENT area. A support reservoir is refused
    /// because its changing footprint/constraint is not supplied by this model.
    ///
    /// Returned storage is enthalpy. Each advance additionally reports pressure
    /// work and internal-energy change from its actual accepted endpoints.
    /// No independent alpha(T), constrained axial response, skin effect,
    /// radiation, melting or mechanical expansion dynamics is implied.
    pub fn try_new_with_free_density_expansion(
        law: ConcentricSectionTemperatureLaw<'a>,
        caloric: &'a ResolvedHeatCapacityCurve,
        reference_length_m: f64,
        boundary: Boundary,
        pressure_axis: &str,
        budget: DepositionBudget,
    ) -> Result<Self, SectionHeatingError> {
        if boundary.reservoir.is_some() {
            return Err(invalid(
                "free density expansion cannot infer a constrained reservoir footprint",
            ));
        }
        let geometry = DensityGeometry::admit(&law, caloric, reference_length_m, pressure_axis)?;
        Self::build(
            law,
            &[caloric],
            reference_length_m,
            &[boundary],
            Vec::new(),
            budget,
            Some(geometry),
        )
    }

    /// Current physical dimensions for the explicitly selected free-expansion
    /// mode. None means the original fixed-geometry section.
    pub fn density_geometry(&self) -> Option<&SectionDensityGeometry> {
        self.observer
            .density_geometry
            .as_ref()
            .map(|geometry| &geometry.current)
    }
}

impl<'a> DensityGeometry<'a> {
    fn admit(
        law: &ConcentricSectionTemperatureLaw<'_>,
        caloric: &'a ResolvedHeatCapacityCurve,
        length_m: f64,
        pressure_axis: &str,
    ) -> Result<Self, SectionHeatingError> {
        positive(length_m)?;
        let reference = law.reference();
        if reference.layers().len() != 1 || pressure_axis == law.bindings[0].axis {
            return Err(invalid(
                "free density expansion requires one homogeneous layer and a separate pressure axis",
            ));
        }
        let binding = &law.bindings[0];
        let pressure_quantity = binding
            .point
            .axis_quantities()
            .get(pressure_axis)
            .copied()
            .ok_or_else(|| invalid("free expansion needs an explicit pressure quantity"))?;
        if pressure_quantity != QuantitySpec::dimensional(Pressure::DIMS)
            && pressure_quantity
                != QuantitySpec::semantic(SemanticType::new(
                    QuantityKind::Pressure,
                    ValueForm::Static,
                ))
        {
            return Err(invalid(
                "free expansion needs explicitly typed static pressure in pascals",
            ));
        }
        let pressure_pa =
            positive(
                *binding.point.axes().get(pressure_axis).ok_or_else(|| {
                    invalid("free expansion source lacks the maintained pressure")
                })?,
            )?;
        let layer = &reference.layers()[0];
        if binding
            .point
            .axis_quantities()
            .iter()
            .any(|(axis, quantity)| axis != &binding.axis && quantity.dims() == Temperature::DIMS)
            || (binding.axis != "T"
                && binding.point.axes().contains_key("T")
                && !binding.point.axis_quantities().contains_key("T"))
        {
            return Err(invalid(
                "free expansion must own the single material temperature",
            ));
        }
        for property in layer.material.properties() {
            // The named, admitted pressure is explicitly held by this owner.
            // Reuse the shared policy for every other source/applicability axis;
            // dimensional p is not silently treated as a stress coordinate.
            let admitted = match binding.point.axis_quantities().get(&binding.axis) {
                Some(&temperature) => property
                    .has_only_owned_coordinates_and_fixed_external_conditions(&[
                        (&binding.axis, temperature),
                        (pressure_axis, pressure_quantity),
                    ]),
                None => property
                    .has_only_owned_and_fixed_external_conditions(pressure_axis, pressure_quantity),
            };
            if !admitted {
                return Err(invalid(
                    "free expansion cannot freeze an unowned stress, strain, rate or field source",
                ));
            }
        }
        let density = reference.density_source_property(0)?;
        let cp = layer
            .material
            .property(SPECIFIC_HEAT_CAPACITY_PROPERTY)
            .ok_or_else(|| invalid("free expansion section lacks original Cp"))?;
        for property in [cp, density] {
            let claim = binding
                .card
                .claims()
                .claim(property.answer().receipt.selected)
                .ok_or_else(|| invalid("original free-expansion source claim is missing"))?;
            if !claim.validity.bounds().contains_key(pressure_axis) {
                return Err(invalid(
                    "both original Cp and density sources must explicitly support the isobar",
                ));
            }
        }
        let volume = positive(reference.area_m2() * length_m)?;
        let mass = positive(reference.line_density_kg_m() * length_m)?;
        let surface = positive(
            2.0 * core::f64::consts::PI
                * (reference.inner_radius_m() + layer.outer_radius_m)
                * length_m
                + 2.0 * reference.area_m2(),
        )?;
        Ok(Self {
            caloric,
            reference_length_m: length_m,
            reference_volume_m3: volume,
            reference_surface_m2: surface,
            pressure_pa,
            current: SectionDensityGeometry {
                mass_kg: mass,
                density_kg_m3: reference.layer_integrals()[0].density_kg_m3(),
                isotropic_scale: 1.0,
                length_m,
                volume_m3: volume,
                surface_area_m2: surface,
            },
        })
    }

    pub(super) fn refresh(
        &mut self,
        law: &ConcentricSectionTemperatureLaw<'_>,
        current_section: ConcentricElasticSection,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<ConcentricElasticSection, SectionHeatingError> {
        let cell = thermal
            .cells()
            .get(&CellId(0))
            .ok_or_else(|| invalid("free expansion lost its original thermal cell"))?;
        if cell.mass_kg != self.current.mass_kg {
            return Err(invalid(
                "free expansion cannot change the original mass inventory",
            ));
        }
        let state = self
            .caloric
            .curve()
            .state_at_specific_enthalpy(cell.specific_h_j_kg)?;
        let rho = state.bulk_density_kg_m3();
        let material_rho = current_section.layer_integrals()[0].density_kg_m3();
        // These are the SAME selected source evaluated through different
        // floating-point interpolation paths, not independently fitted laws.
        if (rho - material_rho).abs() > 128.0 * f64::EPSILON * rho.max(material_rho) {
            return Err(invalid(
                "caloric and mechanical density disagree at the accepted temperature",
            ));
        }
        let volume_m3 = positive(cell.mass_kg / rho)?;
        let scale = positive((volume_m3 / self.reference_volume_m3).cbrt())?;
        let length_m = positive(self.reference_length_m * scale)?;
        let surface_area_m2 = positive(self.reference_surface_m2 * scale * scale)?;
        let reference = law.reference();
        let inner = reference.inner_radius_m() * scale;
        let layer = ConcentricMaterialLayer {
            outer_radius_m: positive(reference.layers()[0].outer_radius_m * scale)?,
            material: current_section.layers()[0].material.clone(),
        };
        let section = match current_section.isotropic_parameters() {
            Some(parameters) => ConcentricElasticSection::try_new_with_isotropic_parameters(
                inner,
                vec![layer],
                parameters,
            )?,
            None => ConcentricElasticSection::try_new(inner, vec![layer])?,
        };
        self.current = SectionDensityGeometry {
            mass_kg: cell.mass_kg,
            density_kg_m3: rho,
            isotropic_scale: scale,
            length_m,
            volume_m3,
            surface_area_m2,
        };
        Ok(section)
    }

    pub(super) fn account(
        &self,
        initial: &DepositionSimulation<'_>,
        final_state: &DepositionSimulation<'_>,
        supplied_heat_j: f64,
    ) -> Result<IsobaricHeatStep, SectionHeatingError> {
        let chart = self.caloric.curve();
        let initial =
            chart.state_at_specific_enthalpy(initial.cells()[&CellId(0)].specific_h_j_kg)?;
        let state =
            chart.state_at_specific_enthalpy(final_state.cells()[&CellId(0)].specific_h_j_kg)?;
        Ok(chart.account_isobaric_transition(
            initial,
            state,
            self.current.mass_kg,
            self.pressure_pa,
            supplied_heat_j,
        )?)
    }
}

fn invalid(what: &'static str) -> SectionHeatingError {
    SectionHeatingError::Invalid(what)
}

#[cfg(test)]
mod tests;
