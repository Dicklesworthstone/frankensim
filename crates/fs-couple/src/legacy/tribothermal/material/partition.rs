//! Explicit short-time heat allocation from the same material and caloric state.
use super::*;
use fs_material::state_point::{
    ResolvedMaterialStatePoint, THERMAL_CONDUCTIVITY_DIMS, THERMAL_CONDUCTIVITY_PROPERTY,
    resolve_material_state_point,
};
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_tribo::SurfaceThermalProperties;

/// Heat-allocation model, chosen explicitly rather than by a material name.
#[derive(Clone, Debug)]
pub enum SlidingHeatPartition<'a> {
    /// Caller-sourced constant fractions, including any named outside channel.
    Prescribed(HeatPartition),
    /// Equal-residence-time, one-dimensional semi-infinite conduction estimate.
    SemiInfiniteEffusivity(SemiInfiniteHeatPartition<'a>),
}

/// Bounded local effusivity approximation for an interface-generated heat pulse.
///
/// Each conductivity is evaluated at its own material's actual temperature.
/// Density and differential heat capacity come from its ORIGINAL enthalpy chart,
/// so a second, independently selected heat-capacity law cannot disagree with
/// thermal storage. The conductivity claims are pinned at initialization.
///
/// This assumes homogeneous isotropic half spaces, equal exposure times and
/// negligible thermal resistance at the heat-generation plane. The two bodies'
/// finite-cell temperatures are local background states; their ordinary signed
/// heat exchange remains the existing thermal network's responsibility. This
/// model does not solve flash temperatures, spatial gradients, moving-source
/// Peclet effects, thermal contact resistance or a third-body layer.
///
/// The caller must bound the distance to EVERY competing boundary (including
/// the footprint edge) by `thermal_depth_m`. The admitted penetration
/// `2 sqrt(alpha_max * age)` must stay below that distance times the declared
/// fraction. `alpha_max` retains the largest diffusivity at accepted endpoints
/// and the current trial; substep refinement and the source's between-endpoint
/// behavior remain caller obligations. Age never resets on a step or restart.
#[derive(Clone, Debug)]
pub struct SemiInfiniteHeatPartition<'a> {
    /// Ordered original bulk materials and their independently typed conditions.
    pub surfaces: [ThermalMaterialPoint<'a>; 2],
    /// Independent source selection for each body's conductivity.
    pub selection: [MaterialPropertySelection; 2],
    /// Distances [m] to the nearest competing thermal boundary on each body.
    pub thermal_depth_m: [f64; 2],
    /// Explicit approximation envelope, strictly between zero and one.
    pub maximum_penetration_fraction: f64,
}

/// Exact physical inputs to one endpoint's heat allocation.
#[derive(Clone, Debug)]
pub struct EffusivityPartitionSample {
    /// Original conductivity claims, query points, units and authority ceilings.
    pub conductivity: [ResolvedMaterialStatePoint; 2],
    /// Current density from the original caloric states, not remassed geometry.
    pub density_kg_m3: [f64; 2],
    /// Reciprocal of the original chart's one-sided dT/dh, in J/(kg K).
    pub specific_heat_capacity_j_kg_k: [f64; 2],
    /// Actual current k/(rho cp), in m²/s.
    pub diffusivity_m2_s: [f64; 2],
    /// Maximum accepted/current endpoint diffusivity retained for the envelope.
    pub maximum_diffusivity_m2_s: [f64; 2],
    /// Current conservative endpoint-sampled diffusion depths [m].
    pub penetration_m: [f64; 2],
    /// Elapsed time since the original contact initialization [s].
    pub contact_age_s: f64,
}

impl SlidingHeatPartition<'_> {
    pub(super) fn other(&self) -> f64 {
        match self {
            Self::Prescribed(partition) => partition.other(),
            Self::SemiInfiniteEffusivity(_) => 0.0,
        }
    }

    pub(super) fn pin(&mut self, sample: Option<&EffusivityPartitionSample>) {
        if let (Self::SemiInfiniteEffusivity(model), Some(sample)) = (self, sample) {
            for side in 0..2 {
                model.selection[side] = MaterialPropertySelection::PinnedByProperty(
                    sample.conductivity[side]
                        .properties()
                        .iter()
                        .map(|property| {
                            (
                                property.requirement().name().to_owned(),
                                property.answer().receipt.selected,
                            )
                        })
                        .collect(),
                );
            }
        }
    }

    pub(super) fn resolve<E>(
        &self,
        data: &SlidingInterface<'_>,
        phases: [EquilibriumPhaseState; 2],
        age_s: f64,
        previous_maximum_diffusivity: [f64; 2],
    ) -> Result<(HeatPartition, Option<EffusivityPartitionSample>), SlidingError<E>> {
        let model = match self {
            Self::Prescribed(partition) => return Ok((*partition, None)),
            Self::SemiInfiniteEffusivity(model) => model,
        };
        if !age_s.is_finite()
            || age_s < 0.0
            || !model.maximum_penetration_fraction.is_finite()
            || model.maximum_penetration_fraction <= 0.0
            || model.maximum_penetration_fraction >= 1.0
            || model
                .thermal_depth_m
                .iter()
                .any(|d| !d.is_finite() || *d <= 0.0)
        {
            return Err(SlidingError::Invalid(
                "semi-infinite heat split needs positive depths, a fraction in (0,1) and nonnegative age",
            ));
        }
        let requirement = ScalarPropertyRequirement::try_new(
            THERMAL_CONDUCTIVITY_PROPERTY,
            THERMAL_CONDUCTIVITY_DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )?;
        let mut conductivity = Vec::with_capacity(2);
        let mut cp = [0.0; 2];
        let mut density = [0.0; 2];
        let mut diffusivity = [0.0; 2];
        let mut maximum_diffusivity = [0.0; 2];
        let mut penetration = [0.0; 2];
        let mut surfaces = Vec::with_capacity(2);
        for side in 0..2 {
            let source = model.surfaces[side];
            // at() also checks the exact original enthalpy-card identity.
            if source.material.content_hash() != data.point.materials[side].content_hash() {
                return Err(SlidingError::Invalid(
                    "effusivity surfaces must preserve ordered original material cards",
                ));
            }
            let point = source.at(phases[side])?;
            let resolved = resolve_material_state_point(
                source.material,
                &point,
                core::slice::from_ref(&requirement),
                model.selection[side].clone(),
            )?;
            let property = resolved
                .property(THERMAL_CONDUCTIVITY_PROPERTY)
                .expect("required conductivity was resolved");
            source.admit_static(property)?;
            let k = property.value_si();
            let slope = data.caloric[side].temperature_derivative_at_specific_enthalpy(
                phases[side].specific_enthalpy_j_kg(),
            )?;
            if slope <= 0.0 {
                return Err(SlidingError::Invalid(
                    "semi-infinite dry heat split needs positive finite differential heat capacity",
                ));
            }
            cp[side] = number(slope.recip())?;
            density[side] = phases[side].bulk_density_kg_m3();
            diffusivity[side] = number(k / density[side] / cp[side])?;
            let surface = SurfaceThermalProperties::new(k, diffusivity[side])?;
            maximum_diffusivity[side] = diffusivity[side].max(previous_maximum_diffusivity[side]);
            penetration[side] = number(2.0 * maximum_diffusivity[side].sqrt() * age_s.sqrt())?;
            if age_s > 0.0 && penetration[side] == 0.0 {
                return Err(SlidingError::Invalid(
                    "thermal penetration is unrepresentable",
                ));
            }
            let allowed = number(model.thermal_depth_m[side] * model.maximum_penetration_fraction)?;
            if allowed <= 0.0 || penetration[side] > allowed {
                return Err(SlidingError::Invalid(
                    "semi-infinite heat split exceeded its diffusion-depth envelope; resolve spatial thermal transport",
                ));
            }
            conductivity.push(resolved);
            surfaces.push(surface);
        }
        let partition = HeatPartition::semi_infinite_surfaces(surfaces[0], surfaces[1])?;
        Ok((
            partition,
            Some(EffusivityPartitionSample {
                conductivity: conductivity
                    .try_into()
                    .expect("exactly two original surfaces"),
                density_kg_m3: density,
                specific_heat_capacity_j_kg_k: cp,
                diffusivity_m2_s: diffusivity,
                maximum_diffusivity_m2_s: maximum_diffusivity,
                penetration_m: penetration,
                contact_age_s: age_s,
            }),
        ))
    }
}
