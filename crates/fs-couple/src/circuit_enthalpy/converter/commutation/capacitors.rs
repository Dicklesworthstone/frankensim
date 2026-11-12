//! Source-bound linear device output capacitance.
//!
//! The charge owner uses four constant positive capacitances throughout one
//! physical history. Temperature and phase still determine whether each
//! original source is applicable. A temperature-, field- or voltage-dependent
//! capacitance needs its own charge law and reversible energy/entropy coupling;
//! this adapter never freezes such a law into a sampled scalar.

use super::*;
use fs_qty::Capacitance;

/// An explicit output-capacitance property on the original device package.
/// The capacitor and its semiconductor conduction law must name the same cell.
#[derive(Clone, Debug)]
pub struct DeviceCapacitanceBinding<'a> {
    pub cell: CellId,
    pub key: &'a PropertyKey,
    pub selection: MaterialPropertySelection,
}

/// The original constant capacitance claim and its current domain admission.
#[derive(Clone, Debug)]
pub struct DeviceCapacitanceMaterial {
    cell: CellId,
    key: PropertyKey,
    source: ResolvedMaterialStatePoint,
    capacitance_f: f64,
}

/// Equal = the same material state (card, material, axes, query point, requirement, source
/// axes, selected claim, value bits) behind the same cell, key and capacitance. The receipt's
/// selection-policy tag and hash are left out on purpose, because a refresh at the same state
/// re-resolves under `PinnedByProperty`.
impl PartialEq for DeviceCapacitanceMaterial {
    fn eq(&self, other: &Self) -> bool {
        self.cell == other.cell
            && self.key == other.key
            && self.capacitance_f == other.capacitance_f
            && same_material(&self.source, &other.source)
    }
}

impl DeviceCapacitanceMaterial {
    pub const fn cell(&self) -> CellId {
        self.cell
    }

    pub const fn key(&self) -> &PropertyKey {
        &self.key
    }

    pub const fn source(&self) -> &ResolvedMaterialStatePoint {
        &self.source
    }

    pub const fn capacitance_f(&self) -> f64 {
        self.capacitance_f
    }
}

/// Refresh source validity at the actual caloric state without changing the
/// accepted charge law, source identity or held external coordinates.
pub(super) fn resolve_capacitances(
    cells: &[MaterialCell<'_>],
    bindings: &[DeviceCapacitanceBinding<'_>; 4],
    thermal: &DepositionSimulation<'_>,
    sample: &MaterialNetworkSample,
    retained: Option<&[DeviceCapacitanceMaterial; 4]>,
    cancel: &RefCell<impl FnMut() -> bool>,
) -> Result<[DeviceCapacitanceMaterial; 4], Error> {
    resolve_constant_capacitances(cells, bindings, thermal, sample, retained, false, cancel)
}

/// The intrinsic gate-drain branch may explicitly author zero capacitance.
/// Existing device output capacitance and gate-source storage stay positive.
pub(super) fn resolve_constant_capacitances(
    cells: &[MaterialCell<'_>],
    bindings: &[DeviceCapacitanceBinding<'_>; 4],
    thermal: &DepositionSimulation<'_>,
    sample: &MaterialNetworkSample,
    retained: Option<&[DeviceCapacitanceMaterial; 4]>,
    allow_zero: bool,
    cancel: &RefCell<impl FnMut() -> bool>,
) -> Result<[DeviceCapacitanceMaterial; 4], Error> {
    let mut capacitors = Vec::with_capacity(bindings.len());
    for (index, binding) in bindings.iter().enumerate() {
        check(cancel)?;
        if binding.key.quantity() != QuantitySpec::dimensional(Capacitance::DIMS) {
            return Err(invalid(
                "device output capacitance requires an exact static farad property",
            ));
        }
        let body = cells
            .iter()
            .find(|body| body.cell == binding.cell)
            .ok_or_else(|| invalid("device capacitance needs its original material package"))?;
        let cell = thermal
            .cells()
            .get(&binding.cell)
            .ok_or_else(|| invalid("device capacitor package is inactive"))?;
        if cell.material_id != body.material_id {
            return Err(invalid(
                "device capacitor package changed its original caloric material",
            ));
        }
        let phase = sample
            .cells
            .iter()
            .find(|entry| entry.cell == binding.cell)
            .ok_or_else(|| invalid("device capacitor package lacks its admitted caloric phase"))?
            .phase;
        let query = body.point.at(phase).map_err(point)?;
        let material = if let Some(originals) = retained {
            let original = &originals[index];
            if original.cell != binding.cell || &original.key != binding.key {
                return Err(invalid(
                    "device capacitance cannot change its package or original source key",
                ));
            }
            for (axis, value) in original.source.query_point() {
                if axis != body.point.temperature_axis
                    && Some(axis.as_str()) != body.point.liquid_fraction_axis
                    && query.axes().get(axis).map(|v| v.to_bits()) != Some(value.to_bits())
                {
                    return Err(invalid(
                        "device capacitance cannot change its held original source context",
                    ));
                }
            }
            original
                .source
                .at_query_point(body.point.material, &query)
                .map_err(source)?
        } else {
            let requirement = ScalarPropertyRequirement::try_with_key(
                binding.key,
                if allow_zero {
                    ScalarAdmissibility::NonNegative
                } else {
                    ScalarAdmissibility::StrictlyPositive
                },
            )
            .map_err(source)?;
            resolve_material_state_point(
                body.point.material,
                &query,
                &[requirement],
                binding.selection.clone(),
            )
            .map_err(source)?
        };
        let property = material
            .property_by_key(binding.key)
            .ok_or_else(|| invalid("device capacitance lost its original source property"))?;
        body.point.admit_static(property).map_err(point)?;
        let static_source = match body.point.liquid_fraction_axis {
            Some(axis) => property.has_only_owned_and_fixed_external_conditions(
                axis,
                QuantitySpec::dimensional(Dims::NONE),
            ),
            None => property.is_static_at_fixed_external_conditions(),
        };
        if !static_source
            || !property.source_axes().is_empty()
            || property.answer().receipt.decision
                != fs_matdb::EvaluationDecision::ConstantWithinValidity
        {
            return Err(invalid(
                "device capacitance requires an original constant law with owned applicability; variable charge laws cannot be frozen",
            ));
        }
        let capacitance_f = property.value_si();
        if !capacitance_f.is_finite()
            || if allow_zero {
                capacitance_f < 0.0
            } else {
                capacitance_f <= 0.0
            }
        {
            return Err(invalid(
                "device capacitance must obey its positive-storage or explicit zero-Miller role",
            ));
        }
        if retained.is_some_and(|originals| {
            originals[index].capacitance_f.to_bits() != capacitance_f.to_bits()
        }) {
            return Err(invalid(
                "device capacitance cannot change along its original charge history",
            ));
        }
        capacitors.push(DeviceCapacitanceMaterial {
            cell: binding.cell,
            key: binding.key.clone(),
            source: material,
            capacitance_f,
        });
    }
    check(cancel)?;
    capacitors
        .try_into()
        .map_err(|_| invalid("four original device capacitance bindings required"))
}
