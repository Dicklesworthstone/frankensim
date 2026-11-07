//! Native source moments and independent component charge/work admission.
use super::*;

pub(super) fn solve<'m>(
    electrodes: &SpatialNonlinearElectrodes<'m>,
    sources: SpatialElectrostaticSources<'_>,
    displacement_faces: &[ElectrostaticDisplacementFace],
    boundary: &[ElectrodeControl],
    initial: &[f64],
    settings: NonlinearElectrostaticControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<SpatialNonlinearElectricField<'m>, Error> {
    poll(cancelled)?;
    super::super::solve::controls(settings)?;
    let problem = electrodes.problem;
    let (dofs, potential, references) = electrodes
        .layout
        .coordinates(boundary, initial, cancelled)?;
    let mut prepared = sources::prepare(
        electrodes,
        &references,
        sources,
        displacement_faces,
        settings,
        cancelled,
    )?;
    let volume_charge_c = prepared.volume_charge_c;
    let (potential, evaluation, report) = gauss::iterate(
        problem,
        &dofs,
        &prepared.rhs,
        potential,
        settings,
        cancelled,
    )?;
    let source_response =
        prepared.response(problem, &evaluation, &potential, &references, cancelled)?;

    let mut terminal_voltage_v = Vec::with_capacity(boundary.len());
    let mut terminal_charge_c = Vec::with_capacity(boundary.len());
    let mut nodal_boundary_charge_c = vec![None; potential.len()];
    let mut component_charge = std::mem::take(&mut prepared.component_charge_c);
    let mut component_correction = std::mem::take(&mut prepared.component_correction_c);
    let mut component_energy = source_response.component_polarization_pairing_j.clone();
    let mut component_energy_correction = vec![0.0; boundary.len()];
    for cell in 0..problem.materials.len() {
        poll(cancelled)?;
        let element = problem.element(cell);
        let reference = references.terminal_by_vertex[element.vertices[0] as usize];
        add(
            &mut component_energy[reference],
            &mut component_energy_correction[reference],
            evaluation.cell_energy[cell],
        )?;
        add(
            &mut component_energy[reference],
            &mut component_energy_correction[reference],
            evaluation.cell_coenergy[cell],
        )?;
    }
    let (mut total, mut ct) = (volume_charge_c, 0.0);
    add(&mut total, &mut ct, prepared.sheet_charge_c)?;
    add(&mut total, &mut ct, -prepared.outward_displacement_c)?;
    let (mut pairing, mut cp) = (0.0, 0.0);
    let mut maximum_floating_charge_residual_c = 0.0_f64;
    for (v, &voltage) in potential.iter().enumerate() {
        poll(cancelled)?;
        let work = finite(voltage * prepared.physical_load_c(v)?)?;
        add(&mut pairing, &mut cp, work)?;
        let reference = references.terminal_by_vertex[v];
        add(
            &mut component_energy[reference],
            &mut component_energy_correction[reference],
            -work,
        )?;
    }
    for (terminal, trace) in electrodes.layout.terminals.iter().enumerate() {
        poll(cancelled)?;
        let (mut charge, mut cq) = (0.0, 0.0);
        for &v in trace {
            poll(cancelled)?;
            let reaction = finite(evaluation.nodal[v] - prepared.rhs[v])?;
            nodal_boundary_charge_c[v] = Some(reaction);
            add(&mut charge, &mut cq, reaction)?;
        }
        let charge = finite(charge + cq)?;
        let voltage = potential[trace[0]];
        terminal_voltage_v.push(voltage);
        terminal_charge_c.push(charge);
        add(&mut total, &mut ct, charge)?;
        let reference = references.terminal_by_vertex[trace[0]];
        add(
            &mut component_charge[reference],
            &mut component_correction[reference],
            charge,
        )?;
        add(&mut pairing, &mut cp, voltage * charge)?;
        add(
            &mut component_energy[reference],
            &mut component_energy_correction[reference],
            -voltage * charge,
        )?;
        if let ElectrodeControl::Charge(requested) = boundary[terminal] {
            maximum_floating_charge_residual_c =
                maximum_floating_charge_residual_c.max((charge - requested).abs());
        }
    }
    let mut maximum_free_charge_residual_c = 0.0_f64;
    for v in 0..potential.len() {
        poll(cancelled)?;
        if electrodes.layout.node_terminal[v].is_none() {
            maximum_free_charge_residual_c =
                maximum_free_charge_residual_c.max((evaluation.nodal[v] - prepared.rhs[v]).abs());
        }
    }
    let charge_balance_residual_c = finite(total + ct)?;
    let mut worst_component_charge_residual_c = 0.0_f64;
    for (&charge, &correction) in component_charge.iter().zip(&component_correction) {
        poll(cancelled)?;
        let residual = finite(charge + correction)?;
        if residual.abs() > worst_component_charge_residual_c.abs() {
            worst_component_charge_residual_c = residual;
        }
    }
    let maximum_component_charge_residual_c = worst_component_charge_residual_c.abs();
    let (mut energy_residual, mut energy_correction) = (0.0, 0.0);
    for work in [
        evaluation.energy,
        evaluation.coenergy,
        source_response.polarization_pairing_j,
        -finite(pairing + cp)?,
    ] {
        add(&mut energy_residual, &mut energy_correction, work)?;
    }
    let energy_identity_residual_j = finite(energy_residual + energy_correction)?;
    let mut worst_component_energy_residual_j = 0.0_f64;
    for (&work, &correction) in component_energy.iter().zip(&component_energy_correction) {
        poll(cancelled)?;
        let residual = finite(work + correction)?;
        if residual.abs() > worst_component_energy_residual_j.abs() {
            worst_component_energy_residual_j = residual;
        }
    }
    let maximum_component_energy_identity_residual_j = worst_component_energy_residual_j.abs();
    if charge_balance_residual_c.abs() > settings.charge_tolerance_c
        || maximum_component_charge_residual_c > settings.charge_tolerance_c
        || maximum_free_charge_residual_c > settings.charge_tolerance_c
        || maximum_floating_charge_residual_c > settings.charge_tolerance_c
        || energy_identity_residual_j.abs() > settings.energy_tolerance_j
        || maximum_component_energy_identity_residual_j > settings.energy_tolerance_j
        || source_response.polarization_work_residual_j.abs() > settings.energy_tolerance_j
        || source_response
            .worst_component_polarization_work_residual_j
            .abs()
            > settings.energy_tolerance_j
    {
        let mut worst_energy = energy_identity_residual_j;
        for residual in [
            worst_component_energy_residual_j,
            source_response.polarization_work_residual_j,
            source_response.worst_component_polarization_work_residual_j,
        ] {
            if residual.abs() > worst_energy.abs() {
                worst_energy = residual;
            }
        }
        return Err(ElectrostaticError::Balance {
            charge_residual_c: if charge_balance_residual_c.abs()
                >= maximum_component_charge_residual_c
            {
                charge_balance_residual_c
            } else {
                worst_component_charge_residual_c
            },
            maximum_free_charge_residual_c: maximum_free_charge_residual_c
                .max(maximum_floating_charge_residual_c),
            energy_residual_j: worst_energy,
        }
        .into());
    }
    poll(cancelled)?;
    let (free_charge_density_c_m3, additional_sources) = prepared.retain(source_response);
    Ok(SpatialNonlinearElectricField {
        problem,
        faces: electrodes.faces.clone(),
        electrodes: electrodes.layout.clone(),
        controls: boundary.to_vec(),
        dofs,
        references,
        relative_voltage_v: potential,
        terminal_voltage_v,
        terminal_charge_c,
        nodal_boundary_charge_c,
        free_charge_density_c_m3,
        additional_sources,
        evaluation,
        volume_charge_c,
        charge_balance_residual_c,
        maximum_component_charge_residual_c,
        maximum_free_charge_residual_c,
        maximum_floating_charge_residual_c,
        energy_identity_residual_j,
        maximum_component_energy_identity_residual_j,
        report,
    })
}
