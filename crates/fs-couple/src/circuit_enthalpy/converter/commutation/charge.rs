//! Complete source-bound output charge Q(V) for the original semiconductor.
//!
//! The native charge owner retains the full monotone law, including support and
//! its integral energy. Only voltage is a response coordinate. Temperature and
//! phase may qualify applicability, but cannot change the charge law without
//! reversible parameter work and an explicitly coupled thermal storage model.

use super::*;
use fs_matdb::{ClaimSelection, PropertyLineSlice};
use fs_phs::circuit::capacitance::{CapacitorChargeCurve, MAX_CAPACITOR_KNOTS};
use fs_qty::ElectricCharge;

/// A complete output-charge law in coulombs on an explicitly owned voltage axis.
/// The original material context must already declare that axis in volts.
/// The range includes both negative and positive voltage and the uncharged zero.
#[derive(Clone, Debug)]
pub struct DeviceChargeBinding<'a> {
    pub cell: CellId,
    pub key: &'a PropertyKey,
    pub voltage_axis: &'a str,
    pub voltage_range_v: [f64; 2],
    pub selection: MaterialPropertySelection,
}

/// Original Q(V) source, complete voltage support and exact native charge law.
/// The scalar source is the zero-voltage reference at the actual caloric state;
/// the retained source path qualifies every voltage used by the native owner.
#[derive(Clone, Debug)]
pub struct DeviceChargeMaterial {
    cell: CellId,
    key: PropertyKey,
    voltage_axis: String,
    voltage_range_v: [f64; 2],
    source: ResolvedMaterialStatePoint,
    path: PropertyLineSlice,
    curve: CapacitorChargeCurve,
}

/// Equal = the same material state (card, material, axes, query point, requirement, source
/// axes, selected claim, value bits) behind the same cell, key, voltage axis, support, path and
/// charge curve. The receipt's selection-policy tag and hash are left out of the source on
/// purpose, because a refresh at the same state re-resolves under `PinnedByProperty`; the path is
/// always queried under a pin of the selected claim, so it compares as it is.
impl PartialEq for DeviceChargeMaterial {
    fn eq(&self, other: &Self) -> bool {
        self.cell == other.cell
            && self.key == other.key
            && self.voltage_axis == other.voltage_axis
            && self.voltage_range_v == other.voltage_range_v
            && same_material(&self.source, &other.source)
            && self.path == other.path
            && self.curve == other.curve
    }
}
impl DeviceChargeMaterial {
    pub const fn cell(&self) -> CellId {
        self.cell
    }
    pub const fn key(&self) -> &PropertyKey {
        &self.key
    }
    pub fn voltage_axis(&self) -> &str {
        &self.voltage_axis
    }
    pub const fn voltage_range_v(&self) -> [f64; 2] {
        self.voltage_range_v
    }
    pub const fn source(&self) -> &ResolvedMaterialStatePoint {
        &self.source
    }
    pub const fn source_path(&self) -> &PropertyLineSlice {
        &self.path
    }
    pub const fn curve(&self) -> &CapacitorChargeCurve {
        &self.curve
    }
}

/// Resolve complete original charge laws at the actual thermal proposal.
/// A source point at V=0 alone never authorizes a frozen local capacitance.
pub(super) fn resolve_charge_curves(
    cells: &[MaterialCell<'_>],
    bindings: &[DeviceChargeBinding<'_>; 4],
    thermal: &DepositionSimulation<'_>,
    sample: &MaterialNetworkSample,
    retained: Option<&[DeviceChargeMaterial; 4]>,
    cancel: &RefCell<impl FnMut() -> bool>,
) -> Result<[DeviceChargeMaterial; 4], Error> {
    let voltage_quantity = QuantitySpec::dimensional(Voltage::DIMS);
    let mut devices = Vec::with_capacity(bindings.len());
    for (index, binding) in bindings.iter().enumerate() {
        check(cancel)?;
        let [lower_v, upper_v] = binding.voltage_range_v;
        if binding.key.quantity() != QuantitySpec::dimensional(ElectricCharge::DIMS)
            || !lower_v.is_finite()
            || !upper_v.is_finite()
            || lower_v >= 0.0
            || upper_v <= 0.0
        {
            return Err(invalid(
                "device charge needs an exact coulomb property and finite support on both sides of zero voltage",
            ));
        }
        let body = cells
            .iter()
            .find(|body| body.cell == binding.cell)
            .ok_or_else(|| {
                invalid("device charge needs its original semiconductor material package")
            })?;
        let cell = thermal
            .cells()
            .get(&binding.cell)
            .ok_or_else(|| invalid("device charge package is inactive"))?;
        if cell.material_id != body.material_id {
            return Err(invalid(
                "device charge package changed its original caloric material",
            ));
        }
        if binding.voltage_axis == body.point.temperature_axis
            || Some(binding.voltage_axis) == body.point.liquid_fraction_axis
            || body
                .point
                .context
                .axis_quantities()
                .get(binding.voltage_axis)
                != Some(&voltage_quantity)
            || !body.point.context.axes().contains_key(binding.voltage_axis)
        {
            return Err(invalid(
                "device charge requires an existing distinct voltage coordinate on its original material context",
            ));
        }
        let phase = sample
            .cells
            .iter()
            .find(|entry| entry.cell == binding.cell)
            .ok_or_else(|| invalid("device charge package lacks its admitted caloric phase"))?
            .phase;
        let query = body
            .point
            .at(phase)
            .map_err(point)?
            .with(binding.voltage_axis, 0.0)
            .map_err(|_| invalid("device zero-voltage reference is unrepresentable"))?;
        let material = if let Some(originals) = retained {
            let original = &originals[index];
            if original.cell != binding.cell
                || &original.key != binding.key
                || original.voltage_axis != binding.voltage_axis
                || original.voltage_range_v.map(f64::to_bits)
                    != binding.voltage_range_v.map(f64::to_bits)
            {
                return Err(invalid(
                    "device charge must retain its original package, source key, voltage axis and full support",
                ));
            }
            for (axis, value) in original.source.query_point() {
                if axis != body.point.temperature_axis
                    && Some(axis.as_str()) != body.point.liquid_fraction_axis
                    && axis != binding.voltage_axis
                    && query.axes().get(axis).map(|v| v.to_bits()) != Some(value.to_bits())
                {
                    return Err(invalid(
                        "device charge changed its held original source context",
                    ));
                }
            }
            original
                .source
                .at_query_point(body.point.material, &query)
                .map_err(source)?
        } else {
            let requirement =
                ScalarPropertyRequirement::try_with_key(binding.key, ScalarAdmissibility::Finite)
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
            .ok_or_else(|| invalid("device charge lost its original property"))?;
        body.point
            .admit_static_with_axes(property, &[(binding.voltage_axis, voltage_quantity)])
            .map_err(point)?;
        let owned_source = match body.point.liquid_fraction_axis {
            Some(axis) => property.has_only_owned_coordinates_and_fixed_external_conditions(&[
                (binding.voltage_axis, voltage_quantity),
                (axis, QuantitySpec::dimensional(Dims::NONE)),
            ]),
            None => property.has_only_owned_and_fixed_external_conditions(
                binding.voltage_axis,
                voltage_quantity,
            ),
        };
        let axes = property.source_axes();
        if !owned_source
            || axes.len() != 1
            || axes[0].0 != binding.voltage_axis
            || axes[0].1 != Voltage::DIMS
        {
            return Err(invalid(
                "device charge must depend only on its owned voltage; temperature-dependent or other unowned charge laws require accounted coupling",
            ));
        }
        if property.value_si() != 0.0 {
            return Err(invalid(
                "device charge law must preserve the explicit uncharged datum Q(0)=0",
            ));
        }
        let lower = query
            .clone()
            .with(binding.voltage_axis, lower_v)
            .map_err(|_| invalid("device lower voltage support is unrepresentable"))?;
        let upper = query
            .clone()
            .with(binding.voltage_axis, upper_v)
            .map_err(|_| invalid("device upper voltage support is unrepresentable"))?;
        check(cancel)?;
        let mut path = body
            .point
            .material
            .claims()
            .query_line_typed(
                binding.key,
                binding.voltage_axis,
                &lower,
                &upper,
                ClaimSelection::Pinned(property.answer().receipt.selected),
            )
            .map_err(|error| {
                source(fs_material::state_point::MaterialStatePointError::Support {
                    property: binding.key.name().to_owned(),
                    source: Box::new(error),
                })
            })?;
        // Inserting the exact original zero answer splits one affine source
        // segment without changing its law, integral, slope or source identity.
        // It is queried under the line's pin, as every other knot is: the state
        // point's own answer carries the bind's single-claim receipt, so a
        // refresh at the same state would not equal its bound path.
        if !path.knots.iter().any(|(voltage, _)| *voltage == 0.0) {
            if path.knots.len() >= MAX_CAPACITOR_KNOTS {
                return Err(invalid(
                    "device charge reference exceeds the native knot budget",
                ));
            }
            let zero = body
                .point
                .material
                .claims()
                .query_pinned_typed(binding.key, &query, property.answer().receipt.selected)
                .map_err(|error| {
                    source(fs_material::state_point::MaterialStatePointError::Support {
                        property: binding.key.name().to_owned(),
                        source: Box::new(fs_matdb::PropertySupportError::Evaluation {
                            point: query.clone(),
                            error: Box::new(error),
                        }),
                    })
                })?;
            if zero.evidence.value.value.to_bits() != property.value_si().to_bits() {
                return Err(invalid(
                    "device charge zero knot differs from its original zero answer",
                ));
            }
            path.knots.push((0.0, zero));
            path.knots.sort_by(|a, b| a.0.total_cmp(&b.0));
        }
        if !(3..=MAX_CAPACITOR_KNOTS).contains(&path.knots.len()) {
            return Err(invalid(
                "bounded complete original device charge knots required",
            ));
        }
        let mut knots = Vec::with_capacity(path.knots.len());
        for (voltage, answer) in &path.knots {
            check(cancel)?;
            if answer.receipt.selected != property.answer().receipt.selected {
                return Err(invalid(
                    "device charge path changed its original source claim",
                ));
            }
            knots.push((answer.evidence.value.value, *voltage));
        }
        let curve = if let Some(originals) = retained {
            let original = &originals[index].curve;
            if knots.as_slice() != original.knots() {
                return Err(invalid(
                    "device charge history cannot change its full voltage law without parameter work",
                ));
            }
            // Fresh source/path admission above remains mandatory. An identical
            // original law can share its admitted slope and energy tables.
            original.clone()
        } else {
            CapacitorChargeCurve::try_new(knots).map_err(|_| invalid(
                "original device charge must be strictly monotone with positive differential capacitance and representable integral energy",
            ))?
        };
        devices.push(DeviceChargeMaterial {
            cell: binding.cell,
            key: binding.key.clone(),
            voltage_axis: binding.voltage_axis.to_owned(),
            voltage_range_v: binding.voltage_range_v,
            source: material,
            path,
            curve,
        });
    }
    check(cancel)?;
    devices
        .try_into()
        .map_err(|_| invalid("four original device charge bindings required"))
}
