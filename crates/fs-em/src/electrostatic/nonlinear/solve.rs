//! Shared nonlinear Gauss iteration and independent original field/source audit.
use super::electrodes::NonlinearElectrodes;
use super::*;

pub(super) fn controls(controls: NonlinearElectrostaticControls) -> Result<(), Error> {
    if !positive(controls.absolute_residual_tolerance_c)
        || !positive(controls.relative_residual_tolerance)
        || controls.relative_residual_tolerance >= 1.0
        || !positive(controls.charge_tolerance_c)
        || !positive(controls.energy_tolerance_j)
        || !positive(controls.linear.tolerance)
        || controls.linear.tolerance >= 1.0
    {
        return Err(invalid(
            "positive finite electrical residual, physical and linear tolerances required",
        ));
    }
    Ok(())
}

pub(super) fn solve<'m>(
    electrodes: &NonlinearElectrodes<'m>,
    sources: super::super::PlanarElectrostaticSources<'_>,
    boundary: &[ElectrodeControl],
    initial: &[f64],
    settings: NonlinearElectrostaticControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<NonlinearElectricField<'m>, Error> {
    poll(cancelled)?;
    controls(settings)?;
    let problem = electrodes.problem;
    let (dofs, potential, references) = electrodes.coordinates(boundary, initial, cancelled)?;
    let sources::PreparedSources {
        loads,
        component_charge_c: mut component_charge,
        component_correction_c: mut component_correction,
    } = sources::prepare(electrodes, &references, sources, cancelled)?;
    let rhs = &loads.rhs;
    let volume_charge_c = loads.volume_charge_c;
    let (potential, evaluation, report) =
        gauss::iterate(problem, &dofs, rhs, potential, settings, cancelled)?;
    let mut terminal_voltage_v = Vec::with_capacity(boundary.len());
    let mut terminal_charge_c = Vec::with_capacity(boundary.len());
    let mut nodal_boundary_charge_c = vec![None; potential.len()];
    let (mut total, mut ct) = (volume_charge_c, 0.0);
    add(&mut total, &mut ct, loads.sheet_charge_c)?;
    add(&mut total, &mut ct, -loads.outward_displacement_c)?;
    let (mut pairing, mut cp) = (0.0, 0.0);
    let mut maximum_floating_charge_residual_c: f64 = 0.0;
    for (v, &q) in rhs.iter().enumerate() {
        add(&mut pairing, &mut cp, potential[v] * q)?;
    }
    for (t, trace) in electrodes.layout.terminals.iter().enumerate() {
        poll(cancelled)?;
        let (mut charge, mut cq) = (0.0, 0.0);
        for &v in trace {
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
        if let ElectrodeControl::Charge(requested) = boundary[t] {
            maximum_floating_charge_residual_c =
                maximum_floating_charge_residual_c.max((charge - requested).abs());
        }
    }
    let mut maximum_free_charge_residual_c: f64 = 0.0;
    for v in 0..potential.len() {
        poll(cancelled)?;
        if electrodes.layout.node_terminal[v].is_none() {
            maximum_free_charge_residual_c =
                maximum_free_charge_residual_c.max((evaluation.nodal[v] - rhs[v]).abs());
        }
    }
    let charge_balance_residual_c = finite(total + ct)?;
    let mut worst_component_charge_residual_c: f64 = 0.0;
    for (&charge, &correction) in component_charge.iter().zip(&component_correction) {
        poll(cancelled)?;
        let residual = finite(charge + correction)?;
        if residual.abs() > worst_component_charge_residual_c.abs() {
            worst_component_charge_residual_c = residual;
        }
    }
    let maximum_component_charge_residual_c = worst_component_charge_residual_c.abs();
    let energy_identity_residual_j = finite(
        evaluation.energy + evaluation.coenergy + evaluation.switching_heat_j() - (pairing + cp),
    )?;
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
    let additional_sources = sources::retain(problem, &evaluation, sources, loads, cancelled)?;
    Ok(NonlinearElectricField {
        problem,
        electrodes: electrodes.layout.clone(),
        controls: boundary.to_vec(),
        dofs,
        references,
        relative_voltage_v: potential,
        terminal_voltage_v,
        terminal_charge_c,
        nodal_boundary_charge_c,
        free_charge_density_c_m3: sources.free_charge_density_c_m3.to_vec(),
        additional_sources,
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
