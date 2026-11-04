//! Actual native axial loss deposited into the original finite layer cells.
//!
//! Perfectly bonded concentric layers share the native axial strain rate.
//! Their heat shares are therefore the SAME eta_i A_i contributions used by
//! the material section, not caloric masses or bending moments. Sources and
//! heat shares are frozen at the accepted start of the interval; mechanics
//! runs once and the original selected claims must support the joint endpoint.

use fs_exec::CancelGate;
use fs_material::section::{ConcentricMaterialLayer, axial::ConcentricAxialViscosity};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
};
use fs_qty::DynViscosity;
use fs_solid::rod::material_field::motion::PreparedRodMotion;
use fs_solid::rod::material_field::thermal::source::RodMaterialSource;

use super::{Result, RodLayerCell, Segment, check, finite, invalid, material, motion, sum};
use crate::rod_motion::{PreparedRodAxialLoss, RodAxialHeat};

/// Explicit solid axial viscosity for every original segment/layer.
#[derive(Clone, Copy, Debug)]
pub struct RodAxialViscosity<'a> {
    /// Nonnegative axial viscosity [Pa s], independently selected from angular loss.
    pub property: &'a str,
    /// Segment-major, then original concentric-layer order; one selection each.
    pub selections: &'a [MaterialPropertySelection],
}

/// One original layer cell receiving its actual intrinsic axial heat.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RodLayerAxialHeat {
    pub destination: RodLayerCell,
    /// Accepted-start source viscosity [Pa s].
    pub viscosity_pa_s: f64,
    /// Original reference annular layer area [m2].
    pub area_m2: f64,
    /// Nonnegative heat deposited into this layer's bound whole cell [J].
    pub heat_j: f64,
}

/// Original physical loss sources and their accepted finite-cell deposits.
#[derive(Clone, Debug, PartialEq)]
pub struct RodAxialHeating {
    /// Frozen accepted-start sources used by the actual mechanical loss flow.
    pub source_sections: Vec<ConcentricAxialViscosity>,
    /// Original layer order and already-bound whole caloric destinations.
    pub layers: Vec<RodLayerAxialHeat>,
    pub mechanical_heat_j: f64,
    pub deposited_heat_j: f64,
    /// Deposited minus native heat [J], gated before either history is published.
    pub partition_residual_j: f64,
}

pub(super) struct PreparedHeating {
    pub(super) loss: PreparedRodAxialLoss,
    layers: Vec<RodLayerAxialHeat>,
}

pub(super) fn prepare<I>(
    model: &PreparedRodMotion<I>,
    source: &RodMaterialSource<'_>,
    segments: &[Segment<'_>],
    viscosity: RodAxialViscosity<'_>,
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
            "axial viscosity selections must match every original rod layer",
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
                "axial source, geometry and caloric layers must agree",
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
                .ok_or_else(|| invalid("resolved axial viscosity is missing"))?
                .value_si();
            layers.push(RodLayerAxialHeat {
                destination: receivers[layer].binding,
                viscosity_pa_s: eta,
                area_m2: geometry.area_m2(),
                heat_j: 0.0,
            });
            enriched.push(ConcentricMaterialLayer {
                outer_radius_m: bound.outer_radius_m,
                material,
            });
        }
        resolved.push(
            current
                .elastic()
                .with_material_layers(enriched)
                .map_err(material)?
                .with_axial_viscosity(viscosity.property)
                .map_err(material)?,
        );
    }
    let loss = PreparedRodAxialLoss::new(model, resolved, gate).map_err(motion)?;
    check(gate)?;
    Ok(PreparedHeating { loss, layers })
}

impl PreparedHeating {
    pub(super) fn partition(
        &self,
        native: &RodAxialHeat,
        tolerance_j: f64,
        gate: &CancelGate,
    ) -> Result<RodAxialHeating> {
        check(gate)?;
        if !tolerance_j.is_finite()
            || tolerance_j < 0.0
            || native.segment_heat_j.len() != self.loss.sections().len()
            || !native.total_heat_j.is_finite()
            || native.total_heat_j < 0.0
        {
            return Err(invalid(
                "native axial heat and partition tolerance must be admitted",
            ));
        }
        let mut layers = self.layers.clone();
        let mut offset = 0;
        for (&heat, section) in native.segment_heat_j.iter().zip(self.loss.sections()) {
            check(gate)?;
            let denominator = section.viscous_rigidity_n_s();
            if !heat.is_finite() || heat < 0.0 || (heat > 0.0 && denominator <= 0.0) {
                return Err(invalid(
                    "positive axial heat needs positive sourced axial loss",
                ));
            }
            let end = offset + section.elastic().layers().len();
            let mut segment_deposit_j = 0.0;
            for (layer, &numerator) in layers[offset..end]
                .iter_mut()
                .zip(section.layer_viscous_rigidity_n_s())
            {
                check(gate)?;
                if heat > 0.0 && numerator > 0.0 {
                    layer.heat_j = finite(heat * (numerator / denominator))?;
                    if layer.heat_j <= 0.0 {
                        return Err(invalid("axial layer heat is unrepresentable"));
                    }
                }
                segment_deposit_j = finite(segment_deposit_j + layer.heat_j)?;
            }
            if finite(segment_deposit_j - heat)?.abs() > tolerance_j {
                return Err(invalid(
                    "layer axial heat does not match its native segment loss",
                ));
            }
            offset = end;
        }
        let deposited_heat_j = sum(layers.iter().map(|l| l.heat_j))?;
        let partition_residual_j = finite(deposited_heat_j - native.total_heat_j)?;
        if partition_residual_j.abs() > tolerance_j {
            return Err(invalid(
                "layer axial heat does not match native mechanical loss",
            ));
        }
        Ok(RodAxialHeating {
            source_sections: self.loss.sections().to_vec(),
            layers,
            mechanical_heat_j: native.total_heat_j,
            deposited_heat_j,
            partition_residual_j,
        })
    }
}

/// Re-admit the SAME selected viscosity claims at the actual endpoint.
/// New loss coefficients apply only to the next accepted interval.
pub(super) fn admit_endpoint<I>(
    model: &PreparedRodMotion<I>,
    source: &RodMaterialSource<'_>,
    segments: &[Segment<'_>],
    viscosity: RodAxialViscosity<'_>,
    heating: &RodAxialHeating,
    gate: &CancelGate,
) -> Result<()> {
    let mut pinned = Vec::with_capacity(heating.layers.len());
    for section in &heating.source_sections {
        for layer in section.elastic().layers() {
            check(gate)?;
            let property = layer
                .material
                .property(viscosity.property)
                .ok_or_else(|| invalid("interval axial source lost its selected viscosity"))?;
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
        RodAxialViscosity {
            property: viscosity.property,
            selections: &pinned,
        },
        gate,
    )?;
    Ok(())
}
