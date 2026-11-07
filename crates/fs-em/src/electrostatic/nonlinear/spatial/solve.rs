//! Exact tetrahedral source moments and independent physical field admission.
use super::*;

fn loads(
    problem: &SpatialNonlinearElectrostatics,
    rho: &[f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Vec<f64>, f64), Error> {
    if rho.len() != problem.materials.len() || rho.iter().any(|v| !v.is_finite()) {
        return Err(invalid(
            "one finite free-charge density per original tetrahedron required",
        ));
    }
    let mut rhs = vec![0.0; problem.mesh().vertex_count];
    let mut correction = vec![0.0; rhs.len()];
    let (mut total, mut ct) = (0.0, 0.0);
    for (cell, &density) in rho.iter().enumerate() {
        poll(cancelled)?;
        let element = problem.element(cell);
        let q = finite(element.volume_m3 * density)?;
        let nodal = finite(0.25 * q)?;
        if density != 0.0 && nodal == 0.0 {
            return Err(invalid("tetrahedral free-charge load underflowed"));
        }
        add(&mut total, &mut ct, q)?;
        for &vertex in element.vertices {
            let v = vertex as usize;
            add(&mut rhs[v], &mut correction[v], nodal)?;
        }
    }
    for (value, c) in rhs.iter_mut().zip(correction) {
        *value = finite(*value + c)?;
    }
    Ok((rhs, finite(total + ct)?))
}

pub(super) fn solve<'m>(
    electrodes: &SpatialNonlinearElectrodes<'m>,
    rho: &[f64],
    boundary: &[ElectrodeControl],
    initial: &[f64],
    settings: NonlinearElectrostaticControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<SpatialNonlinearElectricField<'m>, Error> {
    poll(cancelled)?;
    let problem = electrodes.problem;
    let (dofs, potential, references) = electrodes
        .layout
        .coordinates(boundary, initial, cancelled)?;
    let (rhs, volume_charge_c) = loads(problem, rho, cancelled)?;
    let (potential, evaluation, report) =
        gauss::iterate(problem, &dofs, &rhs, potential, settings, cancelled)?;

    let mut terminal_voltage_v = Vec::with_capacity(boundary.len());
    let mut terminal_charge_c = Vec::with_capacity(boundary.len());
    let mut nodal_boundary_charge_c = vec![None; potential.len()];
    let mut component_charge = vec![0.0; boundary.len()];
    let mut component_correction = vec![0.0; boundary.len()];
    for (cell, &density) in rho.iter().enumerate() {
        poll(cancelled)?;
        let element = problem.element(cell);
        let reference = references.terminal_by_vertex[element.vertices[0] as usize];
        add(
            &mut component_charge[reference],
            &mut component_correction[reference],
            element.volume_m3 * density,
        )?;
    }
    let (mut total, mut ct) = (volume_charge_c, 0.0);
    let (mut pairing, mut cp) = (0.0, 0.0);
    let mut maximum_floating_charge_residual_c = 0.0_f64;
    for (v, &q) in rhs.iter().enumerate() {
        poll(cancelled)?;
        add(&mut pairing, &mut cp, potential[v] * q)?;
    }
    for (terminal, trace) in electrodes.layout.terminals.iter().enumerate() {
        poll(cancelled)?;
        let (mut charge, mut cq) = (0.0, 0.0);
        for &v in trace {
            poll(cancelled)?;
            let reaction = finite(evaluation.nodal[v] - rhs[v])?;
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
                maximum_free_charge_residual_c.max((evaluation.nodal[v] - rhs[v]).abs());
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
    let energy_identity_residual_j =
        finite(evaluation.energy + evaluation.coenergy - (pairing + cp))?;
    if charge_balance_residual_c.abs() > settings.charge_tolerance_c
        || maximum_component_charge_residual_c > settings.charge_tolerance_c
        || maximum_free_charge_residual_c > settings.charge_tolerance_c
        || maximum_floating_charge_residual_c > settings.charge_tolerance_c
        || energy_identity_residual_j.abs() > settings.energy_tolerance_j
    {
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
            energy_residual_j: energy_identity_residual_j,
        }
        .into());
    }
    poll(cancelled)?;
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
        free_charge_density_c_m3: rho.to_vec(),
        evaluation,
        volume_charge_c,
        charge_balance_residual_c,
        maximum_component_charge_residual_c,
        maximum_free_charge_residual_c,
        maximum_floating_charge_residual_c,
        energy_identity_residual_j,
        report,
    })
}
