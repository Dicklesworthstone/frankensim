//! Finite-cell deposition of the original intrinsic rod shear heat.
//!
//! Viscosity is resolved at the accepted typed layer state, while every
//! pre-existing mechanical receipt is retained verbatim. Native interval heat
//! is partitioned by the same frozen layer eta_s J contributions. Thermal trials
//! neither rerun mechanics nor change that partition. Final source support is
//! checked again at the accepted endpoint before either history is published.

use fs_material::section::{
    ConcentricElasticSection, ConcentricMaterialLayer, torsion::ConcentricTorsionalSection,
};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
};
use fs_qty::DynViscosity;

use super::{
    DepositionSimulation, Result, RodLayerCell, RodMotionControls, RodMotionInertia,
    RodNodalWrench, RodThermoelasticReport, RodThermoelasticSolve, RodThermoelasticSystem, Segment,
    StageLoads, check, finite, invalid, material, motion, sum,
};
use crate::rod_motion::{PreparedRodLoads, PreparedRodTorsionalLoss, RodTorsionalHeat};
use fs_exec::CancelGate;
use fs_solid::rod::material_field::motion::PreparedRodMotion;
use fs_solid::rod::material_field::thermal::source::RodMaterialSource;

/// Explicit viscosity property and selection for every original segment/layer.
#[derive(Clone, Copy, Debug)]
pub struct RodTorsionalViscosity<'a> {
    /// Nonnegative dynamic shear viscosity [Pa s], never a damping ratio.
    pub property: &'a str,
    /// Segment-major, then original concentric-layer order; one selection each.
    pub selections: &'a [MaterialPropertySelection],
}

/// One original whole-layer cell receiving its intrinsic shear heat.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RodLayerTorsionalHeat {
    pub destination: RodLayerCell,
    /// Accepted-start source viscosity [Pa s].
    pub viscosity_pa_s: f64,
    /// That layer's original polar area moment J = 2I [m4].
    pub polar_moment_m4: f64,
    /// Nonnegative heat deposited into this layer's bound whole cell [J].
    pub heat_j: f64,
}

/// Source-retaining physical heat partition for one accepted interval.
#[derive(Clone, Debug, PartialEq)]
pub struct RodTorsionalHeating {
    /// Frozen accepted-start sources used by the actual mechanical loss flow.
    pub source_sections: Vec<ConcentricTorsionalSection>,
    /// Original layer order and their already-bound finite caloric destinations.
    pub layers: Vec<RodLayerTorsionalHeat>,
    /// Native total intrinsic mechanical heat [J].
    pub mechanical_heat_j: f64,
    /// Sum of actual per-layer thermal deposits [J].
    pub deposited_heat_j: f64,
    /// Deposited minus native heat [J], independently gated before publication.
    pub partition_residual_j: f64,
}

pub(super) struct PreparedHeating {
    pub(super) loss: PreparedRodTorsionalLoss,
    layers: Vec<RodLayerTorsionalHeat>,
}

pub(super) fn prepare<I>(
    model: &PreparedRodMotion<I>,
    source: &RodMaterialSource<'_>,
    segments: &[Segment<'_>],
    viscosity: RodTorsionalViscosity<'_>,
    gate: &CancelGate,
) -> Result<PreparedHeating> {
    check(gate)?;
    let actual = model.source().sections();
    let count: usize = actual.iter().map(|s| s.elastic().layers().len()).sum();
    if actual.len() != segments.len()
        || actual.len() != source.segments().len()
        || viscosity.selections.len() != count
    {
        return Err(invalid(
            "torsional viscosity selections must match every original rod layer",
        ));
    }
    let requirement = ScalarPropertyRequirement::try_new(
        viscosity.property,
        DynViscosity::DIMS,
        ScalarAdmissibility::NonNegative,
    )
    .map_err(material)?;
    let mut resolved = Vec::with_capacity(actual.len());
    let mut layers = Vec::with_capacity(count);
    for (segment, current) in actual.iter().enumerate() {
        check(gate)?;
        let original = &source.segments()[segment];
        let receivers = &segments[segment].receivers;
        if original.layers.len() != current.elastic().layers().len()
            || receivers.len() != original.layers.len()
        {
            return Err(invalid(
                "torsional source, geometry and caloric layers must agree",
            ));
        }
        let mut enriched = Vec::with_capacity(original.layers.len());
        for (layer, (bound, geometry)) in current
            .elastic()
            .layers()
            .iter()
            .zip(current.elastic().layer_integrals())
            .enumerate()
        {
            check(gate)?;
            let material = bound
                .material
                .with_additional_requirements(
                    original.layers[layer].card,
                    core::slice::from_ref(&requirement),
                    viscosity.selections[layers.len()].clone(),
                )
                .map_err(material)?;
            let eta = material
                .property(viscosity.property)
                .ok_or_else(|| invalid("resolved torsional viscosity is missing"))?
                .value_si();
            let polar = finite(2.0 * geometry.second_moment_m4())?;
            layers.push(RodLayerTorsionalHeat {
                destination: receivers[layer].binding,
                viscosity_pa_s: eta,
                polar_moment_m4: polar,
                heat_j: 0.0,
            });
            enriched.push(ConcentricMaterialLayer {
                outer_radius_m: bound.outer_radius_m,
                material,
            });
        }
        resolved.push(
            ConcentricElasticSection::try_new(current.elastic().inner_radius_m(), enriched)
                .map_err(material)?
                .with_shear_modulus(current.shear_property(), current.shear_correction())
                .map_err(material)?
                .with_torsional_viscosity(viscosity.property)
                .map_err(material)?,
        );
    }
    let loss = PreparedRodTorsionalLoss::new(model, resolved, gate).map_err(motion)?;
    check(gate)?;
    Ok(PreparedHeating { loss, layers })
}

impl PreparedHeating {
    pub(super) fn partition(
        &self,
        native: &RodTorsionalHeat,
        tolerance_j: f64,
        gate: &CancelGate,
    ) -> Result<RodTorsionalHeating> {
        if !tolerance_j.is_finite()
            || tolerance_j < 0.0
            || native.segment_heat_j.len() != self.loss.sections().len()
            || !native.total_heat_j.is_finite()
            || native.total_heat_j < 0.0
        {
            return Err(invalid(
                "native torsional heat and partition tolerance must be admitted",
            ));
        }
        let mut layers = self.layers.clone();
        for layer in &mut layers {
            check(gate)?;
            let segment = layer.destination.segment;
            let heat = native.segment_heat_j[segment];
            let denominator = self.loss.sections()[segment].viscous_stiffness_n_m2_s();
            let numerator = finite(layer.viscosity_pa_s * layer.polar_moment_m4)?;
            if !heat.is_finite() || heat < 0.0 || (heat > 0.0 && denominator <= 0.0) {
                return Err(invalid(
                    "positive torsional heat needs positive sourced shear loss",
                ));
            }
            if heat > 0.0 && numerator > 0.0 {
                layer.heat_j = finite(heat * (numerator / denominator))?;
                if layer.heat_j <= 0.0 {
                    return Err(invalid("torsional layer heat is unrepresentable"));
                }
            }
        }
        let deposited_heat_j = sum(layers.iter().map(|l| l.heat_j))?;
        let partition_residual_j = finite(deposited_heat_j - native.total_heat_j)?;
        if partition_residual_j.abs() > tolerance_j {
            return Err(invalid(
                "layer torsional heat does not match native mechanical loss",
            ));
        }
        Ok(RodTorsionalHeating {
            source_sections: self.loss.sections().to_vec(),
            layers,
            mechanical_heat_j: native.total_heat_j,
            deposited_heat_j,
            partition_residual_j,
        })
    }
}

/// Check the SAME interval viscosity claims at the final material state. This
/// cannot silently select another claim to conceal the original domain failure.
pub(super) fn admit_endpoint<I>(
    model: &PreparedRodMotion<I>,
    source: &RodMaterialSource<'_>,
    segments: &[Segment<'_>],
    viscosity: RodTorsionalViscosity<'_>,
    heating: &RodTorsionalHeating,
    gate: &CancelGate,
) -> Result<()> {
    let mut pinned = Vec::with_capacity(heating.layers.len());
    for section in &heating.source_sections {
        for layer in section.shear().elastic().layers() {
            check(gate)?;
            let property = layer
                .material
                .property(viscosity.property)
                .ok_or_else(|| invalid("interval torsional source lost its selected viscosity"))?;
            pinned.push(MaterialPropertySelection::PinnedByProperty(vec![(
                viscosity.property.to_owned(),
                property.answer().receipt.selected,
            )]));
        }
    }
    prepare(
        model,
        source,
        segments,
        RodTorsionalViscosity {
            property: viscosity.property,
            selections: &pinned,
        },
        gate,
    )?;
    Ok(())
}

impl<'r, 'n, 'a, I: RodMotionInertia> RodThermoelasticSystem<'r, 'n, 'a, I> {
    /// Advance native torsional dissipation and deposit its physical heat into
    /// the original layer cells. Sources/weights are frozen at accepted T for
    /// this interval; the next interval resolves viscosity at its new actual T.
    /// The full mechanical/thermal split remains first order.
    #[allow(clippy::too_many_arguments)]
    pub fn step_torsional(
        &mut self,
        thermal: &mut DepositionSimulation<'_>,
        target_s: f64,
        loads: &[RodNodalWrench],
        viscosity: RodTorsionalViscosity<'_>,
        mechanical_controls: RodMotionControls,
        controls: RodThermoelasticSolve,
        energy_guess_j: &[f64],
        gate: &CancelGate,
    ) -> Result<RodThermoelasticReport> {
        self.step_with_torsion(
            thermal,
            target_s,
            StageLoads::World(loads),
            Some(viscosity),
            mechanical_controls,
            controls,
            energy_guess_j,
            gate,
        )
    }

    /// The same finite-cell shear-heating owner under rotating material-point
    /// forces and body/world couples. No source work or numerical defect becomes
    /// material heat. Failed endpoint source/thermal checks publish neither owner.
    #[allow(clippy::too_many_arguments)]
    pub fn step_spatial_torsional(
        &mut self,
        thermal: &mut DepositionSimulation<'_>,
        target_s: f64,
        loads: &PreparedRodLoads,
        viscosity: RodTorsionalViscosity<'_>,
        mechanical_controls: RodMotionControls,
        controls: RodThermoelasticSolve,
        energy_guess_j: &[f64],
        gate: &CancelGate,
    ) -> Result<RodThermoelasticReport> {
        self.step_with_torsion(
            thermal,
            target_s,
            StageLoads::Spatial(loads),
            Some(viscosity),
            mechanical_controls,
            controls,
            energy_guess_j,
            gate,
        )
    }
}
