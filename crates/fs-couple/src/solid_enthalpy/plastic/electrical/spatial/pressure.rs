//! Native compressive traction on original prescribed matching reference faces.
use super::*;
use crate::circuit_enthalpy::conduction::contacts::SolvedContactPressure;
use crate::solid_enthalpy::plastic::contact::compressive_pressure;
use std::collections::BTreeSet;

/// One original electrical interface whose source pressure comes from native
/// normal stress on BOTH adjacent solid tetrahedra. Pressure is positive in
/// compression; the original reference normals define this small-strain scope.
#[derive(Clone, Copy, Debug)]
pub struct SpatialElectricalContactPressure<'a> {
    /// Original canonical electrical-contact index.
    pub contact: usize,
    /// Existing typed Pa source coordinate, distinct from every thermal role.
    pub pressure_axis: &'a str,
    /// Positive finite absolute allowance for the two normal tractions [Pa].
    /// After agreement the first side's actual value is retained, not averaged.
    pub agreement_tolerance_pa: f64,
}

pub(super) struct PressureRoute<'a> {
    binding: SpatialElectricalContactPressure<'a>,
    elements: [usize; 2],
    normals: [[f64; 3]; 2],
}

pub(super) fn routes<'a, E>(
    solid: &ThermoplasticSolid<'_>,
    problem: &SpatialConduction,
    bindings: &[SpatialElectricalContactPressure<'a>],
    prescribed_m: &[(usize, f64)],
    cx: &Cx<'_>,
) -> Result<Vec<PressureRoute<'a>>, ElectricalSolidSourceError<E>> {
    let mut used = BTreeSet::new();
    let mut result = Vec::with_capacity(bindings.len());
    if bindings.len() > problem.contacts().len() {
        return Err(ElectricalSolidSourceError::Invalid(
            "pressure bindings must name original electrical contacts",
        ));
    }
    for &binding in bindings {
        cx.checkpoint().map_err(|_| DepositionError::Cancelled)?;
        if !positive(binding.agreement_tolerance_pa)
            || binding.pressure_axis.trim().is_empty()
            || !used.insert(binding.contact)
        {
            return Err(ElectricalSolidSourceError::Invalid(
                "unique contact pressure role and positive traction agreement tolerance required",
            ));
        }
        let contact =
            problem
                .contacts()
                .get(binding.contact)
                .ok_or(ElectricalSolidSourceError::Invalid(
                    "foreign electrical contact pressure binding",
                ))?;
        let elements = *problem.contact_cells().get(binding.contact).ok_or(
            ElectricalSolidSourceError::Invalid(
                "electrical contact lacks original adjacent tetrahedra",
            ),
        )?;
        let mut normals = [[0.0; 3]; 2];
        for (side, mut nodes) in [contact.first_vertices, contact.second_vertices]
            .into_iter()
            .enumerate()
        {
            nodes.sort_unstable();
            let face = solid
                .geometry
                .faces()
                .iter()
                .find(|face| face.nodes == nodes)
                .ok_or(ElectricalSolidSourceError::Invalid(
                    "pressure requires complete exterior reference faces",
                ))?;
            if !solid
                .tetrahedra
                .get(elements[side])
                .is_some_and(|tet| nodes.iter().all(|node| tet.contains(node)))
                || nodes.iter().any(|node| {
                    (0..3).any(|axis| {
                        !prescribed_m
                            .iter()
                            .any(|(dof, value)| *dof == 3 * node + axis && value.is_finite())
                    })
                })
            {
                return Err(ElectricalSolidSourceError::Invalid(
                    "pressure contact requires its original tetrahedron and every paired-face displacement prescribed",
                ));
            }
            normals[side] = face.outward_normal;
        }
        if contact
            .first_vertices
            .iter()
            .zip(contact.second_vertices)
            .any(|(&a, b)| solid.nodes_m[a] != solid.nodes_m[b])
            || normals[0]
                .iter()
                .zip(normals[1])
                .map(|(a, b)| a * b)
                .sum::<f64>()
                > -1.0 + 1e-12
        {
            return Err(ElectricalSolidSourceError::Invalid(
                "pressure contact requires coincident reference traces with opposite native normals",
            ));
        }
        result.push(PressureRoute {
            binding,
            elements,
            normals,
        });
    }
    Ok(result)
}

pub(super) fn sample<'a, E>(
    routes: &[PressureRoute<'a>],
    stresses: &[[f64; 6]],
    cx: &Cx<'_>,
) -> Result<Vec<SolvedContactPressure<'a>>, ElectricalSolidSourceError<E>> {
    let mut result = Vec::with_capacity(routes.len());
    for route in routes {
        cx.checkpoint().map_err(|_| DepositionError::Cancelled)?;
        let mut pressures = [0.0; 2];
        for (side, pressure) in pressures.iter_mut().enumerate() {
            let stress =
                stresses
                    .get(route.elements[side])
                    .ok_or(ElectricalSolidSourceError::Invalid(
                        "pressure contact lost original native stress",
                    ))?;
            *pressure = compressive_pressure::<Infallible>(stress, route.normals[side])
                .map_err(ElectricalSolidSourceError::Support)?;
        }
        if (pressures[0] - pressures[1]).abs() > route.binding.agreement_tolerance_pa {
            return Err(ElectricalSolidSourceError::Invalid(
                "two-sided native normal tractions do not agree at the electrical contact",
            ));
        }
        result.push(SolvedContactPressure {
            contact: route.binding.contact,
            axis: route.binding.pressure_axis,
            pressure_pa: pressures[0],
        });
    }
    Ok(result)
}
