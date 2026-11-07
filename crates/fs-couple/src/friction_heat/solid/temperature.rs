//! Complete source temperature ownership on the original ordered heat routes.

use super::SolidContactHeatContact;
use fs_contact::interface_binding::thermal_friction::FrictionQueryAxisKind;
use fs_qty::semantic::QuantityKind;
use fs_qty::{QuantitySpec, Temperature};
use fs_thermal::enthalpy::CellId;
use std::collections::BTreeMap;

/// Original ordered side of the contact's finite thermal route.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SolidContactThermalSide {
    A,
    B,
}

/// An additional source temperature supplied by one original finite body cell.
/// The primary friction temperature always remains the solved common interface
/// temperature. Names must already be declared additional Temperature axes of
/// the original friction law; this mapping cannot add or reinterpret a claim.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SolidContactBodyTemperature<'a> {
    pub node: usize,
    pub axis: &'a str,
    pub side: SolidContactThermalSide,
}

pub(super) struct OwnedTemperature<'a> {
    pub node: usize,
    pub axis: &'a str,
    pub cell: CellId,
}

pub(super) fn admit<'a>(
    contacts: &[&SolidContactHeatContact<'_>],
    bindings: &[SolidContactBodyTemperature<'a>],
) -> Result<Vec<OwnedTemperature<'a>>, &'static str> {
    let by_node: BTreeMap<_, _> = contacts
        .iter()
        .map(|&contact| (contact.node, contact))
        .collect();
    let mut owned = BTreeMap::new();
    for binding in bindings {
        let contact = by_node
            .get(&binding.node)
            .ok_or("body temperature must belong to a bound native contact node")?;
        if !contact.friction.law.query_axes().iter().any(|axis| {
            axis.name == binding.axis && axis.kind == FrictionQueryAxisKind::Temperature
        }) || !contact
            .friction
            .conditions
            .axes()
            .contains_key(binding.axis)
        {
            return Err("body temperature requires an existing additional source temperature axis");
        }
        let side = match binding.side {
            SolidContactThermalSide::A => 0,
            SolidContactThermalSide::B => 1,
        };
        if owned
            .insert((binding.node, binding.axis), contact.route.cells()[side])
            .is_some()
        {
            return Err("one thermal body owner per source temperature axis required");
        }
    }
    for contact in contacts {
        let primary = &contact.friction.law.axes().temperature;
        let conditions = &contact.friction.conditions;
        for axis in contact.friction.law.query_axes() {
            if axis.kind == FrictionQueryAxisKind::Temperature
                && !owned.contains_key(&(contact.node, axis.name.as_str()))
            {
                return Err("additional source temperatures require explicit thermal ownership");
            }
        }
        // Even a temperature admitted as a fixed external condition must have
        // a declared evolving axis and body owner here. Dimensional kelvin is
        // recognized at arbitrary names, alongside semantic absolute kelvin.
        for axis in conditions.axes().keys() {
            if axis == primary {
                continue;
            }
            let quantity = conditions.axis_quantities().get(axis).copied();
            let temperature = quantity == Some(QuantitySpec::dimensional(Temperature::DIMS))
                || quantity
                    .and_then(QuantitySpec::semantic_type)
                    .is_some_and(|q| q.kind() == QuantityKind::AbsoluteTemperature)
                || (axis == "T" && quantity.is_none());
            if temperature && !owned.contains_key(&(contact.node, axis.as_str())) {
                return Err("additional source temperatures require explicit thermal ownership");
            }
        }
    }
    // Canonical ordering makes mapping order immaterial to query/cache behavior.
    Ok(owned
        .into_iter()
        .map(|((node, axis), cell)| OwnedTemperature { node, axis, cell })
        .collect())
}
