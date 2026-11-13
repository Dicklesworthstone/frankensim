//! Finite-cell deposition of the original native rod bending loss.
//!
//! Each accepted-start layer supplies its own sourced viscosity and area moment
//! eta I. Both native bending directions have these same layer weights. Their
//! actual integrated loss is deposited once into the original caloric cells;
//! thermal trials neither repeat mechanics nor repartition heat by cell mass.
//! The selected source claims must also support the accepted endpoint state.

use fs_exec::CancelGate;
use fs_material::section::{ConcentricKelvinVoigtSection, ConcentricMaterialLayer};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
};
use fs_qty::DynViscosity;
use fs_solid::rod::material_field::motion::PreparedRodMotion;
use fs_solid::rod::material_field::thermal::source::RodMaterialSource;

use super::{Result, RodLayerCell, Segment, check, finite, invalid, material, motion, sum};
use crate::rod_motion::{PreparedRodBendingLoss, RodBendingHeat};

/// Explicit uniaxial bending viscosity for every original segment/layer.
#[derive(Clone, Copy, Debug)]
pub struct RodBendingViscosity<'a> {
    /// Nonnegative solid bending viscosity [Pa s], independent of shear loss.
    pub property: &'a str,
    /// Segment-major, then original concentric-layer order; one selection each.
    pub selections: &'a [MaterialPropertySelection],
}

/// One original whole-layer cell receiving its actual intrinsic bending heat.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RodLayerBendingHeat {
    pub destination: RodLayerCell,
    /// Accepted-start source viscosity [Pa s].
    pub viscosity_pa_s: f64,
    /// Original layer area moment I about either circular bending axis [m4].
    pub second_moment_m4: f64,
    /// Nonnegative heat deposited into this layer's bound whole cell [J].
    pub heat_j: f64,
}

/// Source-retaining bending heat partition for one accepted native interval.
#[derive(Clone, Debug, PartialEq)]
pub struct RodBendingHeating {
    /// Frozen accepted-start sources used by the actual mechanical loss flow.
    pub source_sections: Vec<ConcentricKelvinVoigtSection>,
    /// Original layer order and their already-bound finite caloric destinations.
    pub layers: Vec<RodLayerBendingHeat>,
    /// Actual native heat summed over both bending directions [J].
    pub mechanical_heat_j: f64,
    /// Sum of actual per-layer thermal deposits [J].
    pub deposited_heat_j: f64,
    /// Deposited minus native heat [J], gated before either owner is published.
    pub partition_residual_j: f64,
}

pub(super) struct PreparedHeating {
    pub(super) loss: PreparedRodBendingLoss,
    layers: Vec<RodLayerBendingHeat>,
}

pub(super) fn prepare<I>(
    model: &PreparedRodMotion<I>,
    source: &RodMaterialSource<'_>,
    segments: &[Segment<'_>],
    viscosity: RodBendingViscosity<'_>,
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
            "bending viscosity selections must match every original rod layer",
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
                "bending source, geometry and caloric layers must agree",
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
                .ok_or_else(|| invalid("resolved bending viscosity is missing"))?
                .value_si();
            layers.push(RodLayerBendingHeat {
                destination: receivers[layer].binding,
                viscosity_pa_s: eta,
                second_moment_m4: geometry.second_moment_m4(),
                heat_j: 0.0,
            });
            enriched.push(ConcentricMaterialLayer {
                outer_radius_m: bound.outer_radius_m,
                material,
            });
        }
        resolved.push(
            current.elastic().with_material_layers(enriched)
                .map_err(material)?
                .with_kelvin_voigt_bending(viscosity.property)
                .map_err(material)?,
        );
    }
    let loss = PreparedRodBendingLoss::new(model, resolved, gate).map_err(motion)?;
    check(gate)?;
    Ok(PreparedHeating { loss, layers })
}

impl PreparedHeating {
    pub(super) fn partition(
        &self,
        native: &RodBendingHeat,
        tolerance_j: f64,
        gate: &CancelGate,
    ) -> Result<RodBendingHeating> {
        check(gate)?;
        if !tolerance_j.is_finite()
            || tolerance_j < 0.0
            || native.segment_heat_j.len() != self.loss.sections().len()
            || !native.total_heat_j.is_finite()
            || native.total_heat_j < 0.0
        {
            return Err(invalid(
                "native bending heat and partition tolerance must be admitted",
            ));
        }
        let mut layers = self.layers.clone();
        let mut offset = 0;
        for (&heat, section) in native.segment_heat_j.iter().zip(self.loss.sections()) {
            check(gate)?;
            let denominator = section.viscous_stiffness_n_m2_s();
            if !heat.is_finite() || heat < 0.0 || (heat > 0.0 && denominator <= 0.0) {
                return Err(invalid(
                    "positive bending heat needs positive sourced bending loss",
                ));
            }
            let end = offset + section.elastic().layers().len();
            let mut segment_deposit_j = 0.0;
            for (layer, &numerator) in layers[offset..end]
                .iter_mut()
                .zip(section.layer_viscous_stiffness_n_m2_s())
            {
                check(gate)?;
                if heat > 0.0 && numerator > 0.0 {
                    layer.heat_j = finite(heat * (numerator / denominator))?;
                    if layer.heat_j <= 0.0 {
                        return Err(invalid("bending layer heat is unrepresentable"));
                    }
                }
                segment_deposit_j = finite(segment_deposit_j + layer.heat_j)?;
            }
            if finite(segment_deposit_j - heat)?.abs() > tolerance_j {
                return Err(invalid(
                    "layer bending heat does not match its native segment loss",
                ));
            }
            offset = end;
        }
        let deposited_heat_j = sum(layers.iter().map(|l| l.heat_j))?;
        let partition_residual_j = finite(deposited_heat_j - native.total_heat_j)?;
        if partition_residual_j.abs() > tolerance_j {
            return Err(invalid(
                "layer bending heat does not match native mechanical loss",
            ));
        }
        Ok(RodBendingHeating {
            source_sections: self.loss.sections().to_vec(),
            layers,
            mechanical_heat_j: native.total_heat_j,
            deposited_heat_j,
            partition_residual_j,
        })
    }
}

/// Re-admit the interval's selected viscosity claims at the actual endpoint.
/// A new claim cannot replace an out-of-domain original selection.
pub(super) fn admit_endpoint<I>(
    model: &PreparedRodMotion<I>,
    source: &RodMaterialSource<'_>,
    segments: &[Segment<'_>],
    viscosity: RodBendingViscosity<'_>,
    heating: &RodBendingHeating,
    gate: &CancelGate,
) -> Result<()> {
    let mut pinned = Vec::with_capacity(heating.layers.len());
    for section in &heating.source_sections {
        for layer in section.elastic().layers() {
            check(gate)?;
            let property = layer
                .material
                .property(viscosity.property)
                .ok_or_else(|| invalid("interval bending source lost its selected viscosity"))?;
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
        RodBendingViscosity {
            property: viscosity.property,
            selections: &pinned,
        },
        gate,
    )?;
    Ok(())
}
