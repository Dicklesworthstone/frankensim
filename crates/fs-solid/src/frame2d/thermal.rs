//! Prescribed expansion through the original force-based frame equilibrium.
//!
//! Material-coordinate mass and geometry stay at the initial stress-free
//! reference. Each existing Lobatto station may have its own current elastic
//! coefficients and integrated layer expansion. Axial force and elastic storage
//! come from `ConcentricThermoelasticSection`; the same native flexibility,
//! compatibility, corotational tangent, supports and time integrator consume
//! the response. Opposing layer stresses are never replaced by net force.
//!
//! This is one-way prescribed-temperature small-strain elasticity of concentric
//! circumferentially uniform layers. It supplies no heat equation, reciprocal
//! thermoelastic heat, radial expansion, thermal curvature, or remeshing. The
//! separate linear harmonic/viscosity and resultant-only strength adapters do
//! not admit this spatial thermal field.

use std::sync::Arc;

use fs_exec::Cx;
use fs_material::section::ConcentricElasticSection;
use fs_material::section::thermal::{ConcentricAxialResponse, ConcentricThermoelasticSection};
use fs_material::state_point::MaterialStatePointError;

use super::{ElementCommit, Frame2d, SectionModel};
use crate::{SolidError, fiber::SectionState};

fn invalid(what: &str) -> SolidError {
    SolidError::MaterialRefused { what: what.into() }
}
fn material(error: MaterialStatePointError) -> SolidError {
    invalid(&error.to_string())
}
fn checkpoint(cx: &Cx<'_>) -> Result<(), SolidError> {
    cx.checkpoint().map_err(|_| SolidError::Cancelled)
}

/// Immutable current station response with its original integrated source paths.
#[derive(Clone, Debug)]
pub struct FrameThermalSection {
    source: Arc<ConcentricThermoelasticSection>,
    linear_strain_limit: f64,
}

impl FrameThermalSection {
    /// Current layer material states and exact reference-to-current integrals.
    #[must_use]
    pub fn source(&self) -> &ConcentricThermoelasticSection {
        &self.source
    }

    /// Bound applied to total, thermal and elastic outer-fiber axial strains.
    #[must_use]
    pub const fn linear_strain_limit(&self) -> f64 {
        self.linear_strain_limit
    }

    // No allocation in repeated native residual/tangent evaluations.
    fn evaluate(&self, eps0: f64, kappa: f64) -> Result<(SectionState, f64), SolidError> {
        let (n, axial_energy) = self.source
            .axial_force_and_energy_at_strain(eps0, self.linear_strain_limit)
            .map_err(material)?;
        if !kappa.is_finite() {
            return Err(invalid("thermal frame curvature must be finite"));
        }
        for (layer, path) in self.source.section().layers().iter()
            .zip(self.source.expansion_paths())
        {
            let bending_strain = kappa.abs() * layer.outer_radius_m;
            if !bending_strain.is_finite()
                || eps0.abs() + bending_strain > self.linear_strain_limit
                || (eps0 - path.free_linear_strain()).abs() + bending_strain
                    > self.linear_strain_limit
            {
                return Err(invalid("thermal frame outer-fiber strain requires finite-strain or inelastic escalation"));
            }
        }
        let ea = self.source.section().axial_stiffness_n();
        let ei = self.source.section().bending_stiffness_n_m2();
        let m = ei * kappa;
        let bending_energy = 0.5 * m * kappa;
        let energy = axial_energy + bending_energy;
        if !m.is_finite() || !energy.is_finite()
            || (kappa != 0.0 && (m == 0.0 || bending_energy == 0.0))
        {
            return Err(invalid("thermal frame bending response is unrepresentable"));
        }
        Ok((SectionState { n, m, tangent: [[ea, 0.0], [0.0, ei]] }, energy))
    }

    pub(super) fn respond(&self, eps0: f64, kappa: f64) -> Result<SectionState, SolidError> {
        self.evaluate(eps0, kappa).map(|(response, _)| response)
    }
}

/// Recovery at one accepted native station, in the original member coordinates.
#[derive(Clone, Debug, PartialEq)]
pub struct FrameThermalSectionResponse {
    /// Layer centerline stresses, including self-equilibrated thermal mismatch.
    pub axial: ConcentricAxialResponse,
    /// Native signed bending curvature [1/m].
    pub curvature_per_m: f64,
    /// Native signed section bending moment [N m].
    pub bending_moment_n_m: f64,
    /// Axial plus bending elastic storage per ORIGINAL member length [J/m].
    /// This is stored mechanical energy, not dissipated heat.
    pub strain_energy_j_m: f64,
    /// Minimum/maximum signed axial normal stress in each annulus [Pa].
    /// Includes bending and the complete original layer thermal mismatch.
    pub layer_stress_range_pa: Vec<[f64; 2]>,
}

impl Frame2d {
    /// Bind one original-source thermal section per existing Lobatto station.
    ///
    /// The member must have been built by `add_material_element` from its
    /// stress-free reference material states. All current sections must retain
    /// exactly those cards, property selections, reference temperatures and
    /// layer geometry. Current E(T) reaches native stiffness, but rho(T) never
    /// replaces the original mass or rotary inertia of this closed specimen.
    /// Repeated refreshes also pin every original expansion source/reference.
    ///
    /// The current nodal pose and velocities remain intact. Native compatibility
    /// is solved at that pose to refresh section strains and internal forces;
    /// `static_load` can then find its new supported equilibrium, or the original
    /// time owner can advance the frozen-temperature mechanics. The caller owns
    /// the heat/changed-storage balance. Shear and viscous members refuse until
    /// their additional source operators can be refreshed coherently.
    ///
    /// All inputs are immutable and failures/cancellation publish no candidate.
    pub fn with_thermal_material_sections(
        &self,
        element: usize,
        sections: Vec<ConcentricThermoelasticSection>,
        linear_strain_limit: f64,
        cx: &Cx<'_>,
    ) -> Result<Self, SolidError> {
        checkpoint(cx)?;
        let member = self.elements.get(element)
            .ok_or_else(|| invalid("thermal frame member index is out of range"))?;
        let reference = member.material_section.as_ref()
            .ok_or_else(|| invalid("thermal frame requires an original source-bound reference section"))?;
        if member.transverse_shear.is_some() || member.viscous_section.is_some() {
            return Err(invalid("thermal frame refresh requires an elastic Euler-Bernoulli member without a separate viscous operator"));
        }
        if sections.len() != member.sections.len() {
            return Err(invalid("thermal frame needs exactly one source section per existing Lobatto station"));
        }
        let mut prepared = Vec::with_capacity(sections.len());
        for (point, source) in sections.into_iter().enumerate() {
            checkpoint(cx)?;
            admit_reference(reference, &source)?;
            match &member.sections[point] {
                SectionModel::Elastic { .. } => {},
                SectionModel::ThermalElastic(old) => {
                    if old.source.expansion_paths().iter().zip(source.expansion_paths())
                        .any(|(old, new)| old.reference() != new.reference())
                    {
                        return Err(invalid("thermal frame refresh must retain the original expansion source and stress-free reference"));
                    }
                },
                SectionModel::Fiber(_) => {
                    return Err(invalid("thermal elastic refresh cannot replace an inelastic fiber history"));
                },
                SectionModel::LaminateThermal(_) => {
                    return Err(invalid("concentric thermal refresh cannot replace a laminate source section"));
                },
            }
            let section = FrameThermalSection {
                source: Arc::new(source),
                linear_strain_limit,
            };
            section.respond(member.commit.d[point].0, member.commit.d[point].1)?;
            prepared.push(SectionModel::ThermalElastic(section));
        }
        let mut candidate = self.clone();
        candidate.elements[element].sections = prepared;
        checkpoint(cx)?;
        let trial = candidate.element_trial(&candidate.elements[element], &candidate.u)?;
        candidate.elements[element].commit = ElementCommit {
            q: trial.q,
            d: trial.d,
            alpha: trial.alpha,
        };
        checkpoint(cx)?;
        Ok(candidate)
    }

    /// Current integrated thermal source at an existing station, if one is bound.
    #[must_use]
    pub fn thermal_material_section(
        &self,
        element: usize,
        station: usize,
    ) -> Option<&FrameThermalSection> {
        match self.elements.get(element)?.sections.get(station)? {
            SectionModel::ThermalElastic(section) => Some(section),
            _ => None,
        }
    }

    /// Read complete layer stresses and storage at the accepted station strains.
    /// No source query, equilibrium solve, material update or thermal time elapses.
    pub fn thermal_section_response(
        &self,
        element: usize,
        station: usize,
    ) -> Result<FrameThermalSectionResponse, SolidError> {
        let section = self.thermal_material_section(element, station)
            .ok_or_else(|| invalid("thermal frame recovery requires a bound material station"))?;
        let (eps0, kappa) = self.elements[element].commit.d[station];
        let (response, energy) = section.evaluate(eps0, kappa)?;
        let axial = section.source.response_at_axial_strain(eps0, section.linear_strain_limit)
            .map_err(material)?;
        let layer_stress_range_pa = axial.layers().iter()
            .zip(section.source.section().layers())
            .zip(section.source.section().layer_integrals())
            .map(|((response, layer), integral)| {
                let spread = integral.young_modulus_pa() * (kappa.abs() * layer.outer_radius_m);
                let range = [response.stress_pa() - spread, response.stress_pa() + spread];
                if range.iter().any(|value| !value.is_finite()) {
                    Err(invalid("thermal frame layer bending stress is unrepresentable"))
                } else {
                    Ok(range)
                }
            }).collect::<Result<Vec<_>, _>>()?;
        Ok(FrameThermalSectionResponse {
            axial,
            curvature_per_m: kappa,
            bending_moment_n_m: response.m,
            strain_energy_j_m: energy,
            layer_stress_range_pa,
        })
    }
}

fn admit_reference(
    reference: &ConcentricElasticSection,
    thermal: &ConcentricThermoelasticSection,
) -> Result<(), SolidError> {
    let current = thermal.section();
    if current.inner_radius_m().to_bits() != reference.inner_radius_m().to_bits()
        || current.layers().len() != reference.layers().len()
        || current.isotropic_parameters() != reference.isotropic_parameters()
    {
        return Err(invalid("thermal frame must retain its original section geometry and elastic source roles"));
    }
    for (index, ((original, now), path)) in reference.layers().iter().zip(current.layers())
        .zip(thermal.expansion_paths()).enumerate()
    {
        let stress_free = path.reference().resolved();
        if original.outer_radius_m.to_bits() != now.outer_radius_m.to_bits()
            || original.material.card_identity() != now.material.card_identity()
            || stress_free.card_identity() != original.material.card_identity()
            || stress_free.axis_quantities() != original.material.axis_quantities()
            || stress_free.query_point() != original.material.query_point()
        {
            return Err(invalid("thermal frame expansion must start at the original member's stress-free source state"));
        }
        let old_properties = reference.elastic_source_properties(index).map_err(material)?;
        let new_properties = current.elastic_source_properties(index).map_err(material)?;
        if old_properties.len() != new_properties.len() {
            return Err(invalid("thermal frame elastic property roles changed"));
        }
        for (old, new) in old_properties.into_iter().zip(new_properties) {
            if old.requirement() != new.requirement()
                || old.answer().receipt.selected != new.answer().receipt.selected
                || !new.is_constant_at_fixed_environment()
                || new.angular_frequency_band_rad_s().map_err(material)?.is_some()
            {
                return Err(invalid("thermal frame requires the original selected elastic claims at a fixed temperature without frequency dependence"));
            }
        }
    }
    Ok(())
}
