//! One solved matching-interface pressure shared with original thermal links.
use super::*;
use crate::circuit_enthalpy::conduction::contacts::SolvedContactPressure;
use crate::circuit_enthalpy::contacts::heat::{PreparedHeatNetwork, ThermalContactPressure};
use std::collections::BTreeSet;

/// Explicit correspondence between one native electrical interface and the
/// SAME original material thermal contact. Heat retains its reference face
/// area and authored bulk half-paths; only its source pressure evolves here.
#[derive(Clone, Copy, Debug)]
pub struct SpatialElectricalThermalContact {
    pub contact: usize,
    pub thermal_link: u64,
}

pub(super) struct ThermalPressureRoute<'a> {
    contact: usize,
    link: u64,
    axis: &'a str,
}

pub(super) fn routes<'a, E>(
    solid: &ThermoplasticSolid<'_>,
    thermal: &DepositionSimulation<'_>,
    problem: &SpatialConduction,
    contacts: &[ConductionContactThermal<'_>],
    pressures: &[SpatialElectricalContactPressure<'a>],
    mappings: &[SpatialElectricalThermalContact],
    prepared: &PreparedHeatNetwork<'_, '_>,
    cx: &Cx<'_>,
) -> Result<Vec<ThermalPressureRoute<'a>>, ElectricalSolidSourceError<E>> {
    let invalid = ElectricalSolidSourceError::Invalid;
    if mappings.len() > pressures.len()
        || mappings.len() > prepared.thermal_contact_bindings().len()
    {
        return Err(invalid(
            "joint contact mappings require original solved-pressure and thermal owners",
        ));
    }
    let mut used_contacts = BTreeSet::new();
    let mut used_links = BTreeSet::new();
    let mut routes = Vec::with_capacity(mappings.len());
    for mapping in mappings {
        cx.checkpoint().map_err(|_| DepositionError::Cancelled)?;
        if !used_contacts.insert(mapping.contact) || !used_links.insert(mapping.thermal_link) {
            return Err(invalid(
                "joint electrical/thermal interfaces must map one-to-one",
            ));
        }
        let pressure = pressures
            .iter()
            .find(|p| p.contact == mapping.contact)
            .ok_or(invalid(
                "thermal contact requires the same selected native pressure owner",
            ))?;
        let electrical = contacts
            .iter()
            .find(|c| c.contact == mapping.contact)
            .ok_or(invalid(
                "joint contact lacks original electrical source binding",
            ))?;
        let binding = prepared
            .thermal_contact_bindings()
            .iter()
            .find(|c| c.link == mapping.thermal_link)
            .ok_or(invalid(
                "joint contact lacks original pinned thermal source binding",
            ))?;
        let link = thermal
            .links()
            .get(&mapping.thermal_link)
            .ok_or(invalid("joint thermal contact link is inactive"))?;
        let field_contact = problem.contacts().get(mapping.contact).ok_or(invalid(
            "joint contact is absent from original native field",
        ))?;
        let elements = *problem.contact_cells().get(mapping.contact).ok_or(invalid(
            "joint contact lost its original adjacent tetrahedra",
        ))?;
        let cells = elements.map(|element| solid.cells[element].accepted.id);
        let a = electrical.point;
        let b = binding.point;
        let flexible = |axis: &str| {
            axis == pressure.pressure_axis
                || a.temperature_axes.contains(&axis)
                || a.liquid_fraction_axes.contains(&Some(axis))
                || a.mean_temperature_axis == Some(axis)
        };
        if cells != [link.a, link.b]
            || cells != binding.sides.map(|side| side.cell)
            || a.interface.content_hash() != b.interface.content_hash()
            || a.materials
                .iter()
                .zip(b.materials)
                .any(|(x, y)| x.content_hash() != y.content_hash())
            || a.temperature_axes != b.temperature_axes
            || a.liquid_fraction_axes != b.liquid_fraction_axes
            || a.mean_temperature_axis != b.mean_temperature_axis
            || !b.context.axes().contains_key(pressure.pressure_axis)
            || a.context
                .axes()
                .keys()
                .chain(b.context.axes().keys())
                .any(|axis| {
                    !flexible(axis)
                        && (a.context.axes().get(axis) != b.context.axes().get(axis)
                            || a.context.axis_quantities().get(axis)
                                != b.context.axis_quantities().get(axis))
                })
        {
            return Err(invalid(
                "joint electrical/thermal contact must retain the same ordered interface, cells and physical conditions",
            ));
        }
        for (side, element) in binding.sides.iter().zip(elements) {
            if side.material_id != solid.cells[element].accepted.material_id
                || !thermal.uses_material_law(side.material_id, side.curve as &dyn ThermalLaw)
            {
                return Err(invalid(
                    "joint contact must retain original ordered caloric objects",
                ));
            }
        }
        let mut nodes = field_contact.first_vertices;
        nodes.sort_unstable();
        let face = solid
            .geometry
            .faces()
            .iter()
            .find(|face| face.nodes == nodes)
            .ok_or(invalid(
                "joint thermal contact requires the original exterior reference face",
            ))?;
        if !crate::solid_enthalpy::same_measure(link.area_m2, face.area_m2) {
            return Err(invalid(
                "joint thermal contact area must equal the original reference face area",
            ));
        }
        routes.push(ThermalPressureRoute {
            contact: mapping.contact,
            link: mapping.thermal_link,
            axis: pressure.pressure_axis,
        });
    }
    Ok(routes)
}

pub(super) fn values<'a, E>(
    routes: &[ThermalPressureRoute<'a>],
    pressures: &[SolvedContactPressure<'_>],
    cx: &Cx<'_>,
) -> Result<Vec<ThermalContactPressure<'a>>, ElectricalSolidSourceError<E>> {
    routes
        .iter()
        .map(|route| {
            cx.checkpoint().map_err(|_| DepositionError::Cancelled)?;
            let pressure = pressures
                .iter()
                .find(|p| p.contact == route.contact && p.axis == route.axis)
                .ok_or(ElectricalSolidSourceError::Invalid(
                    "joint thermal link lost its original native pressure sample",
                ))?;
            Ok(ThermalContactPressure {
                link: route.link,
                pressure_axis: route.axis,
                pressure_pa: pressure.pressure_pa,
            })
        })
        .collect()
}
