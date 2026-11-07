//! Bind the actual accepted thermal inventory and measure native solid storage.

use std::collections::BTreeSet;

use super::*;
use fs_matdb::QueryPoint;
use fs_material::state_point::{DENSITY_PROPERTY, ResolvedMaterialStatePoint};
use fs_solid::hyper3::HyperelasticLawField;
use fs_thermal::enthalpy::ThermalLaw;

pub(super) struct Source {
    pub tetrahedra: Vec<[u32; 4]>,
    pub density: Vec<[f64; 4]>,
    pub enthalpy: Vec<f64>,
    pub enthalpy_j: Dd,
    pub elastic_j: Dd,
    pub absolute_elastic_j: Dd,
    pub mass_defect_kg: f64,
    pub volume_defect_m3: f64,
    pub maximum_temperature_defect_k: f64,
    pub maximum_cauchy_stress_pa: f64,
}

pub(super) fn preflight(
    body: &FiniteSolidBody<'_>,
    thermal: &DepositionSimulation<'_>,
    cells: &[CellId],
    domain: MaterialFreeSurfaceDomain3<'_>,
    controls: FiniteSolidLiquidReleaseControls,
) -> Result<(), FiniteSolidLiquidReleaseError> {
    let maximum_nodes = controls
        .transfer
        .mechanical
        .max_tetrahedra
        .checked_mul(4)
        .ok_or_else(|| invalid("source node budget overflows"))?;
    if body.problem.tetrahedra.is_empty()
        || body.problem.tetrahedra.len() > controls.transfer.mechanical.max_tetrahedra
        || body.problem.nodes_m.len() > maximum_nodes
        || cells.len() != body.problem.tetrahedra.len()
        || thermal.cells().len() != cells.len()
        || !body.state.time_s.is_finite()
        || body.state.time_s < 0.0
        || body.state.time_s.to_bits() != thermal.time_s().to_bits()
        || ![
            controls.maximum_cauchy_stress_pa,
            controls.maximum_elastic_storage_j,
            controls.maximum_inertial_force_n,
            controls.source_kinetic_tolerance_j,
            controls.source_centroid_tolerance_m,
            controls.source_temperature_tolerance_k,
            controls.source_reference_volume_tolerance_m3,
        ]
        .into_iter()
        .all(|x| x.is_finite() && x >= 0.0)
    {
        return Err(invalid(
            "release requires bounded original cells, one accepted clock and finite physical accuracy controls",
        ));
    }
    if !body.problem.prescribed_m.is_empty()
        || body.pressure_surface.is_some()
        || body.gravity.is_some()
        || body
            .problem
            .nodal_forces_n
            .iter()
            .any(|(_, force)| *force != 0.0)
        || body.problem.body_force_n_m3 != [0.0; 3]
        || domain.gravity_m_s2 != [0.0; 3]
        || !thermal.pending_events().is_empty()
        || !thermal.links().is_empty()
    {
        return Err(invalid(
            "retained supports, loads, gravity, pressure, thermal links or future events require their own discharge before release",
        ));
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
pub(super) fn source(
    body: &FiniteSolidBody<'_>,
    thermal: &DepositionSimulation<'_>,
    cells: &[CellId],
    displacement: &[[f64; 3]],
    scaling: MaterialFreeSurfaceScaling<'_>,
    law: &LiquidTemperatureLaw<'_>,
    caloric: &EquilibriumEnthalpyPhaseCurve,
    controls: FiniteSolidLiquidReleaseControls,
    cx: &Cx<'_>,
    gate: &CancelGate,
) -> Result<Source, FiniteSolidLiquidReleaseError> {
    let receiver = scaling.material().material();
    if law.card().content_hash() != receiver.card_identity()
        || caloric.material_card_identity() != receiver.card_identity()
    {
        return Err(invalid(
            "solid, thermal chart and liquid transport require one original immutable material card",
        ));
    }
    if let Some(field) = body.reference_density_field() {
        for vertex in field.vertices() {
            check(cx, gate)?;
            admit_density_source(vertex.material(), receiver, law.temperature_axis())?;
        }
    }
    let mut nodal_h = vec![None; body.problem.nodes_m.len()];
    let mut used = BTreeSet::new();
    let mut out = Source {
        tetrahedra: Vec::with_capacity(cells.len()),
        density: Vec::with_capacity(cells.len()),
        enthalpy: Vec::new(),
        enthalpy_j: Dd::ZERO,
        elastic_j: Dd::ZERO,
        absolute_elastic_j: Dd::ZERO,
        mass_defect_kg: 0.0,
        volume_defect_m3: 0.0,
        maximum_temperature_defect_k: 0.0,
        maximum_cauchy_stress_pa: 0.0,
    };
    for (element, (&indices, &id)) in body.problem.tetrahedra.iter().zip(cells).enumerate() {
        check(cx, gate)?;
        let cell = thermal
            .cells()
            .get(&id)
            .ok_or_else(|| invalid("original thermal cell is inactive"))?;
        if !used.insert(id)
            || !thermal.uses_material_law(cell.material_id, caloric as &dyn ThermalLaw)
            || cell.boundary.source_w != 0.0
            || cell.boundary.htc_w_m2_k != 0.0
            || cell.boundary.reservoir.is_some()
        {
            return Err(invalid(
                "unique original cells and exact caloric chart without retained heat ports required",
            ));
        }
        let phase = caloric.state_at_specific_enthalpy(cell.specific_h_j_kg)?;
        if !phase.is_fully_liquid() {
            return Err(invalid(
                "native thermal cells must already be fully liquid; this handoff does not melt a solid",
            ));
        }
        let material = match body.problem.material {
            HyperelasticLawField::Uniform(material) => material,
            HyperelasticLawField::PerElement(materials) => &materials[element],
        };
        admit_density_source(material.density(), receiver, law.temperature_axis())?;
        if !same_current_context(material.law().context(), receiver, law.temperature_axis()) {
            return Err(invalid(
                "current mechanical and liquid transport external source contexts differ",
            ));
        }
        let source_temperature = material
            .law()
            .context()
            .axes()
            .get(law.temperature_axis())
            .copied()
            .ok_or_else(|| invalid("mechanical source has no original held caloric temperature"))?;
        let temperature_defect = (source_temperature - phase.temperature_k()).abs();
        out.maximum_temperature_defect_k = out.maximum_temperature_defect_k.max(temperature_defect);
        if !temperature_defect.is_finite()
            || temperature_defect > controls.source_temperature_tolerance_k
        {
            return Err(invalid(
                "accepted mechanical and thermal temperatures disagree",
            ));
        }
        let response = body
            .prepared
            .element_response_with(element, |node| displacement[node], cx)
            .map_err(FiniteSolidError::from)?;
        let stress_norm = response
            .cauchy_pa
            .iter()
            .fold(0.0_f64, |n, value| n.hypot(*value));
        out.maximum_cauchy_stress_pa = out.maximum_cauchy_stress_pa.max(stress_norm);
        out.elastic_j = out.elastic_j + Dd::from_f64(response.potential_j);
        out.absolute_elastic_j = out.absolute_elastic_j + Dd::from_f64(response.potential_j).abs();
        if !stress_norm.is_finite()
            || stress_norm > controls.maximum_cauchy_stress_pa
            || !finite_dd(out.elastic_j)
            || !finite_dd(out.absolute_elastic_j)
            || Dd::from_f64(controls.maximum_elastic_storage_j).lt(out.absolute_elastic_j)
        {
            return Err(invalid(
                "native solid retains elastic stress or storage; a physical relaxation/disposition owner is required",
            ));
        }
        let density = match body.reference_density_field() {
            Some(field) => indices.map(|node| field.vertices()[node].density_kg_m3()),
            None => [material.reference_density_kg_m3(); 4],
        };
        let mass =
            response.reference_volume_m3 * density.into_iter().map(|rho| rho * 0.25).sum::<f64>();
        out.mass_defect_kg += (cell.mass_kg - mass).abs();
        out.volume_defect_m3 += (cell.geometry.volume_m3 - response.reference_volume_m3).abs();
        out.enthalpy_j =
            out.enthalpy_j + Dd::from_f64(cell.mass_kg) * Dd::from_f64(cell.specific_h_j_kg);
        if !mass.is_finite()
            || mass <= 0.0
            || !out.mass_defect_kg.is_finite()
            || !out.volume_defect_m3.is_finite()
            || !finite_dd(out.enthalpy_j)
            || out.mass_defect_kg > controls.transfer.mechanical.mass_tolerance_kg
            || out.volume_defect_m3 > controls.source_reference_volume_tolerance_m3
        {
            return Err(invalid(
                "thermal cells do not retain the original solid reference mass and volume",
            ));
        }
        let mut native_indices = [0u32; 4];
        for (corner, node) in indices.into_iter().enumerate() {
            native_indices[corner] = u32::try_from(node)
                .map_err(|_| invalid("original topology exceeds native liquid vertex indices"))?;
            if nodal_h[node].is_some_and(|old| old != cell.specific_h_j_kg) {
                return Err(invalid(
                    "incident thermal cells disagree at a shared vertex; an explicit conservative P0-to-P1 thermal projection is required",
                ));
            }
            nodal_h[node] = Some(cell.specific_h_j_kg);
        }
        out.tetrahedra.push(native_indices);
        out.density.push(density);
    }
    out.enthalpy = nodal_h
        .into_iter()
        .map(|h| {
            h.ok_or_else(|| invalid("original mechanical vertex has no accepted thermal inventory"))
        })
        .collect::<Result<_, _>>()?;
    check(cx, gate)?;
    Ok(out)
}

fn same_current_context(
    original: &QueryPoint,
    receiver: &ResolvedMaterialStatePoint,
    temperature_axis: &str,
) -> bool {
    original.axis_quantities() == receiver.axis_quantities()
        && original
            .axes()
            .iter()
            .filter(|(axis, _)| axis.as_str() != temperature_axis)
            .map(|(axis, value)| (axis.as_str(), *value))
            .eq(receiver
                .query_point()
                .iter()
                .filter(|(axis, _)| axis != temperature_axis)
                .map(|(axis, value)| (axis.as_str(), *value)))
}

fn admit_density_source(
    original: &ResolvedMaterialStatePoint,
    receiver: &ResolvedMaterialStatePoint,
    temperature_axis: &str,
) -> Result<(), FiniteSolidLiquidReleaseError> {
    let source = original
        .property(DENSITY_PROPERTY)
        .ok_or_else(|| invalid("original solid matter lost its density source"))?;
    let target = receiver
        .property(DENSITY_PROPERTY)
        .ok_or_else(|| invalid("native liquid lost its density source"))?;
    if original.card_identity() != receiver.card_identity()
        || source.answer().receipt.selected != target.answer().receipt.selected
        || original.axis_quantities() != receiver.axis_quantities()
        || !original
            .query_point()
            .iter()
            .filter(|(axis, _)| axis != temperature_axis)
            .eq(receiver
                .query_point()
                .iter()
                .filter(|(axis, _)| axis != temperature_axis))
    {
        return Err(invalid(
            "original reference matter must retain the liquid's pinned density claim and nonthermal source context",
        ));
    }
    Ok(())
}
