//! Reversible storage and heat rates at actual native laminate stations.
//!
//! The shared material potential owns all temperature derivatives, including
//! the relaxation of free-sided transverse strains. Native accepted axial
//! strain and curvature supply its mechanical state. Integration uses the
//! member's original reference length and existing Lobatto weights, exactly as
//! its force-based compatibility does. No alternate displacement interpolation,
//! current-density remassing or continuous inter-station accuracy is inferred.
//!
//! Elastic capacity supplements the original caloric capacity and may be
//! negative. The signed heat contributions below are reversible exchange, not
//! dissipation. These queries do not advance a mechanical or thermal timeline
//! and do not infer finite-interval work from endpoint motion.

use super::{Frame2d, FrameLaminateSection, SectionModel, checkpoint, invalid, source_error};
use crate::SolidError;
use fs_exec::Cx;
use fs_material::laminate::thermal::beam::thermoelastic::{
    LaminateBeamThermoelasticState, ResolvedLaminateBeamThermoelastic,
};
use fs_material::laminate::thermal::thermoelastic::AffineLaminaThermoelastic;
use std::sync::Arc;

/// One original layer's reversible contribution per reference member length.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FrameLaminateThermoelasticLayer {
    pub layer: usize,
    pub temperature_k: f64,
    pub source: LaminateBeamThermoelasticState,
    /// Conservative elastic capacity bound over the complete original source
    /// interval and strain envelope. Add the actual caloric capacity before
    /// admitting a thermal solve; this quantity alone is not total capacity.
    pub minimum_elastic_heat_capacity_j_m_k: f64,
}

/// Source potentials at one accepted force-based member station.
#[derive(Clone, Debug, PartialEq)]
pub struct FrameLaminateThermoelasticResponse {
    pub element: usize,
    pub station: usize,
    /// Original x/L from the member's start node.
    pub station_coordinate: f64,
    /// Original length times this station's Lobatto weight [m].
    pub reference_length_weight_m: f64,
    pub axial_strain: f64,
    /// Native sign: positive curvature shortens positive-y fibers.
    pub curvature_per_m: f64,
    pub layers: Vec<FrameLaminateThermoelasticLayer>,
}

/// Reference-Lobatto integral of elastic storage for one original member.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct FrameLaminateThermoelasticEnergy {
    pub free_energy_j: f64,
    pub entropy_j_k: f64,
    pub internal_energy_j: f64,
    /// Elastic contribution for a common temperature increment in all layers
    /// and stations, holding each native total axial strain/curvature fixed.
    /// This is not the effective capacity after free mechanical equilibration.
    pub elastic_heat_capacity_j_k: f64,
    pub minimum_elastic_heat_capacity_j_k: f64,
}

/// Signed reversible layer heat at one original station, per reference length.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FrameLaminateThermoelasticLayerHeatRate {
    pub layer: usize,
    pub temperature_k: f64,
    pub elastic_heat_capacity_j_m_k: f64,
    /// T * (N_T * axial_strain_rate + M_T * curvature_rate) [W/m].
    pub deformation_heat_w_m: f64,
    /// -C_elastic * temperature_rate [W/m], with either sign permitted.
    pub temperature_heat_w_m: f64,
    /// Sum of the preceding two terms; positive supplies the caloric owner.
    pub reversible_heat_w_m: f64,
}

/// Actual station heat contributions and their original integration weight.
#[derive(Clone, Debug, PartialEq)]
pub struct FrameLaminateThermoelasticHeatRate {
    pub element: usize,
    pub station: usize,
    pub reference_length_weight_m: f64,
    pub layers: Vec<FrameLaminateThermoelasticLayerHeatRate>,
}

fn finite(value: f64) -> Result<f64, SolidError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid(
            "laminate thermoelastic aggregate is unrepresentable",
        ))
    }
}

fn accepted(
    frame: &Frame2d,
    element: usize,
    station: usize,
) -> Result<
    (
        &FrameLaminateSection,
        &ResolvedLaminateBeamThermoelastic,
        [f64; 2],
    ),
    SolidError,
> {
    let member = frame
        .elements
        .get(element)
        .ok_or_else(|| invalid("thermoelastic laminate member index is out of range"))?;
    let original = member
        .laminate_thermal
        .as_ref()
        .ok_or_else(|| invalid("thermoelastic query requires its original source-bound member"))?;
    let section = frame
        .laminate_thermal_section(element, station)
        .ok_or_else(|| invalid("thermoelastic query requires an original laminate station"))?;
    if original.thermoelastic.is_none() || !Arc::ptr_eq(&section.binding, original) {
        return Err(invalid(
            "thermoelastic query must retain its original member and station binding",
        ));
    }
    let law = section.thermoelastic().ok_or_else(|| {
        invalid("laminate reversible potentials must be explicitly bound before querying")
    })?;
    let (eps, kappa) = member.commit.d[station];
    section.admit_strain(eps, kappa)?;
    Ok((section, law, [eps, kappa]))
}

impl Frame2d {
    /// Attach smooth reversible laws to the SAME original laminate sources.
    /// Each original layer needs an explicit complete temperature interval;
    /// every selected elasticity/CTE coefficient must be affine there. Gaps,
    /// slope changes and unsupported source directions refuse through the
    /// shared material owner instead of inventing entropy at a transition.
    ///
    /// The native tensor-norm strain limit remains unchanged. The derived
    /// principal engineering-component envelope is sqrt(2) times that limit,
    /// which also covers rotated engineering shear. The current accepted pose
    /// is admitted before publication. Geometry, inertia, strength selections,
    /// forces and motion are preserved. A second binding refuses; subsequent
    /// temperature refreshes retain these exact source laws and intervals.
    pub fn with_laminate_thermoelasticity(
        &self,
        element: usize,
        temperature_ranges_k: &[[f64; 2]],
        cx: &Cx<'_>,
    ) -> Result<Self, SolidError> {
        checkpoint(cx)?;
        let member = self
            .elements
            .get(element)
            .ok_or_else(|| invalid("thermoelastic laminate member index is out of range"))?;
        let original = member
            .laminate_thermal
            .as_ref()
            .ok_or_else(|| invalid("thermoelasticity requires an original laminate member"))?;
        if original.thermoelastic.is_some()
            || temperature_ranges_k.len() != original.source.layers().len()
        {
            return Err(invalid(
                "one original interval per layer and no thermoelastic rebinding required",
            ));
        }
        let mut laws = Vec::with_capacity(temperature_ranges_k.len());
        let limit = finite(std::f64::consts::SQRT_2 * original.linear_strain_limit)?;
        for (source, &range) in original.source.layers().iter().zip(temperature_ranges_k) {
            checkpoint(cx)?;
            laws.push(
                AffineLaminaThermoelastic::new(source, range, limit)
                    .map_err(|error| source_error(error, cx))?,
            );
        }
        let laws = Arc::new(laws);
        let mut prepared = Vec::with_capacity(member.sections.len());
        for (station, section) in member.sections.iter().enumerate() {
            checkpoint(cx)?;
            let SectionModel::LaminateThermal(section) = section else {
                return Err(invalid(
                    "thermoelastic binding must retain every original station",
                ));
            };
            if !Arc::ptr_eq(&section.binding, original) || section.thermoelastic.is_some() {
                return Err(invalid(
                    "thermoelastic binding cannot replace a station source",
                ));
            }
            let current = ResolvedLaminateBeamThermoelastic::try_new(
                &laws,
                &section.resolved,
                &section.prepared,
                || cx.checkpoint().is_err(),
            )
            .map_err(|error| source_error(error, cx))?;
            let (eps, kappa) = member.commit.d[station];
            section.admit_strain(eps, kappa)?;
            current
                .state(eps, kappa, || cx.checkpoint().is_err())
                .map_err(|error| source_error(error, cx))?;
            prepared.push(Arc::new(current));
        }
        let mut binding = (**original).clone();
        binding.thermoelastic = Some(laws);
        let binding = Arc::new(binding);
        let mut candidate = self.clone();
        for (section, prepared) in candidate.elements[element]
            .sections
            .iter_mut()
            .zip(prepared)
        {
            let SectionModel::LaminateThermal(section) = section else {
                unreachable!()
            };
            section.binding = Arc::clone(&binding);
            section.thermoelastic = Some(prepared);
        }
        candidate.elements[element].laminate_thermal = Some(binding);
        checkpoint(cx)?;
        Ok(candidate)
    }

    /// Actual accepted station strain/curvature and original layer potentials.
    /// Includes residual bimetal storage even when net N and M vanish.
    pub fn laminate_thermoelastic_response(
        &self,
        element: usize,
        station: usize,
        cx: &Cx<'_>,
    ) -> Result<FrameLaminateThermoelasticResponse, SolidError> {
        checkpoint(cx)?;
        let (_, law, [eps, kappa]) = accepted(self, element, station)?;
        let mut layers = Vec::with_capacity(law.layer_count());
        for layer in 0..law.layer_count() {
            checkpoint(cx)?;
            layers.push(FrameLaminateThermoelasticLayer {
                layer,
                temperature_k: law
                    .layer_temperature_k(layer)
                    .ok_or_else(|| invalid("original layer temperature is absent"))?,
                source: law
                    .layer_state(layer, eps, kappa)
                    .map_err(|error| source_error(error, cx))?,
                minimum_elastic_heat_capacity_j_m_k: law
                    .layer_minimum_heat_capacity_j_m_k(layer)
                    .map_err(|error| source_error(error, cx))?,
            });
        }
        let member = &self.elements[element];
        checkpoint(cx)?;
        Ok(FrameLaminateThermoelasticResponse {
            element,
            station,
            station_coordinate: member.xi[station],
            reference_length_weight_m: member.l0 * member.w[station],
            axial_strain: eps,
            curvature_per_m: kappa,
            layers,
        })
    }

    /// Integrate accepted elastic potentials with the member's ORIGINAL
    /// Lobatto quadrature. No allocation or material query occurs. This is the
    /// native discretization's storage, not a continuous interpolation claim.
    /// Add caloric capacity to the reported elastic capacity before thermal use.
    pub fn laminate_thermoelastic_energy(
        &self,
        element: usize,
        cx: &Cx<'_>,
    ) -> Result<FrameLaminateThermoelasticEnergy, SolidError> {
        checkpoint(cx)?;
        let member = self
            .elements
            .get(element)
            .ok_or_else(|| invalid("thermoelastic laminate member index is out of range"))?;
        let mut result = FrameLaminateThermoelasticEnergy::default();
        for station in 0..member.sections.len() {
            checkpoint(cx)?;
            let (_, law, [eps, kappa]) = accepted(self, element, station)?;
            let weight = member.l0 * member.w[station];
            for layer in 0..law.layer_count() {
                checkpoint(cx)?;
                let source = law
                    .layer_state(layer, eps, kappa)
                    .map_err(|error| source_error(error, cx))?;
                let minimum = law
                    .layer_minimum_heat_capacity_j_m_k(layer)
                    .map_err(|error| source_error(error, cx))?;
                result.free_energy_j =
                    finite(result.free_energy_j + weight * source.free_energy_j_m)?;
                result.entropy_j_k = finite(result.entropy_j_k + weight * source.entropy_j_m_k)?;
                result.internal_energy_j =
                    finite(result.internal_energy_j + weight * source.internal_energy_j_m)?;
                result.elastic_heat_capacity_j_k =
                    finite(result.elastic_heat_capacity_j_k + weight * source.heat_capacity_j_m_k)?;
                result.minimum_elastic_heat_capacity_j_k =
                    finite(result.minimum_elastic_heat_capacity_j_k + weight * minimum)?;
            }
        }
        checkpoint(cx)?;
        Ok(result)
    }

    /// Instantaneous signed elastic feedback to the original caloric owner.
    /// `generalized_strain_rate` is [TOTAL centerline axial strain rate (1/s),
    /// TOTAL native curvature rate (1/(m s))] at this accepted station, including
    /// thermal deformation. It is not elastic-strain rate or a force/load rate.
    /// `temperature_rates_k_s` retains original layer order. A coupled solver
    /// must supply compatible physical rates from its existing mechanical and
    /// thermal owners; this query does not infer them from endpoint motion.
    ///
    /// Each layer returns Q_rev = T*(N_T*epsdot + M_T*kappadot) - C_el*Tdot
    /// per original length. Multiply by `reference_length_weight_m` for its
    /// quadrature contribution. The same residual may instead put C_el on its
    /// capacity side and use only deformation heat on the supply side; never
    /// apply both representations. Total caloric+elastic capacity must remain
    /// positive. No dissipation, external supply or finite-interval work is
    /// asserted by these signed reversible terms.
    pub fn laminate_thermoelastic_heat_rate(
        &self,
        element: usize,
        station: usize,
        generalized_strain_rate: [f64; 2],
        temperature_rates_k_s: &[f64],
        cx: &Cx<'_>,
    ) -> Result<FrameLaminateThermoelasticHeatRate, SolidError> {
        checkpoint(cx)?;
        let (_, law, [eps, kappa]) = accepted(self, element, station)?;
        if temperature_rates_k_s.len() != law.layer_count()
            || generalized_strain_rate
                .iter()
                .chain(temperature_rates_k_s)
                .any(|x| !x.is_finite())
        {
            return Err(invalid(
                "finite native total-strain rates and one temperature rate per original layer required",
            ));
        }
        let mut layers = Vec::with_capacity(law.layer_count());
        for (layer, &temperature_rate) in temperature_rates_k_s.iter().enumerate() {
            checkpoint(cx)?;
            let state = law
                .layer_state(layer, eps, kappa)
                .map_err(|error| source_error(error, cx))?;
            let temperature_k = law
                .layer_temperature_k(layer)
                .ok_or_else(|| invalid("original layer temperature is absent"))?;
            let [n_t, m_t] = state.generalized_force_temperature_derivative;
            let deformation_heat_w_m = finite(
                temperature_k
                    * n_t.mul_add(generalized_strain_rate[0], m_t * generalized_strain_rate[1]),
            )?;
            let temperature_heat_w_m = finite(-state.heat_capacity_j_m_k * temperature_rate)?;
            layers.push(FrameLaminateThermoelasticLayerHeatRate {
                layer,
                temperature_k,
                elastic_heat_capacity_j_m_k: state.heat_capacity_j_m_k,
                deformation_heat_w_m,
                temperature_heat_w_m,
                reversible_heat_w_m: finite(deformation_heat_w_m + temperature_heat_w_m)?,
            });
        }
        checkpoint(cx)?;
        Ok(FrameLaminateThermoelasticHeatRate {
            element,
            station,
            reference_length_weight_m: self.elements[element].l0
                * self.elements[element].w[station],
            layers,
        })
    }
}
