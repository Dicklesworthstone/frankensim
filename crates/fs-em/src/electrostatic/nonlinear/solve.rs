//! Physical nonlinear residual, bounded line search and independent field audit.
use super::electrodes::{Dofs, NonlinearElectrodes};
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

fn loads(
    problem: &PlanarNonlinearElectrostatics,
    rho: &[f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Vec<f64>, f64), Error> {
    if rho.len() != problem.elements.len() || rho.iter().any(|v| !v.is_finite()) {
        return Err(invalid(
            "one finite free-charge density per original cell required",
        ));
    }
    let mut rhs = vec![0.0; problem.mesh.vertices().len()];
    let mut correction = vec![0.0; rhs.len()];
    let (mut total, mut ct) = (0.0, 0.0);
    for (element, &density) in problem.elements.iter().zip(rho) {
        poll(cancelled)?;
        let q = finite(element.volume_m3 * density)?;
        if density != 0.0 && q / 3.0 == 0.0 {
            return Err(invalid("free charge load underflowed"));
        }
        add(&mut total, &mut ct, q)?;
        for v in element.vertices {
            add(&mut rhs[v], &mut correction[v], q / 3.0)?;
        }
    }
    for (v, c) in rhs.iter_mut().zip(correction) {
        *v = finite(*v + c)?;
    }
    Ok((rhs, finite(total + ct)?))
}

fn objective_change(
    problem: &PlanarNonlinearElectrostatics,
    dofs: &Dofs,
    rhs: &[f64],
    before_v: &[f64],
    before: &Evaluation,
    after_v: &[f64],
    after: &Evaluation,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, Error> {
    let (mut sum, mut correction) = (problem.coenergy_change(before, after, cancelled)?, 0.0);
    for (v, &q) in rhs.iter().enumerate() {
        poll(cancelled)?;
        add(&mut sum, &mut correction, -q * (after_v[v] - before_v[v]))?;
    }
    for (i, &v) in dofs.representatives.iter().enumerate() {
        add(
            &mut sum,
            &mut correction,
            -dofs.charge[i] * (after_v[v] - before_v[v]),
        )?;
    }
    Ok(finite(sum + correction)?)
}

pub(super) fn solve<'m>(
    electrodes: &NonlinearElectrodes<'m>,
    rho: &[f64],
    boundary: &[ElectrodeControl],
    initial: &[f64],
    settings: NonlinearElectrostaticControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<NonlinearElectricField<'m>, Error> {
    poll(cancelled)?;
    controls(settings)?;
    let problem = electrodes.problem;
    let (dofs, mut potential, reference_terminal, reference_voltage_v) =
        electrodes.coordinates(boundary, initial, cancelled)?;
    let (rhs, volume_charge_c) = loads(problem, rho, cancelled)?;
    let mut evaluation = Evaluation::new(problem);
    let mut trial = Evaluation::new(problem);
    problem.evaluate(&potential, &mut evaluation, cancelled)?;
    let mut residual = vec![0.0; dofs.len()];
    let mut correction = vec![0.0; dofs.len()];
    dofs.residual(&evaluation, &rhs, &mut residual, &mut correction)?;
    let initial_norm = finite(norm(&residual))?;
    let target = finite(
        settings
            .absolute_residual_tolerance_c
            .max(settings.relative_residual_tolerance * initial_norm)
            .min(settings.charge_tolerance_c / (2.0 * (dofs.len().max(1) as f64).sqrt())),
    )?;
    if !positive(target) {
        return Err(invalid("physical charge target is unrepresentable"));
    }
    let mut report = NonlinearElectrostaticReport {
        newton_iterations: 0,
        linear_iterations: 0,
        backtracks: 0,
        operator_applications: 0,
        residual_norm_c: initial_norm,
        residual_target_c: target,
    };
    let mut proposed = potential.clone();
    while report.residual_norm_c > target {
        poll(cancelled)?;
        if report.newton_iterations == settings.max_newton_iterations {
            return Err(Error::Nonconvergence(report));
        }
        let rhs_newton: Vec<_> = residual.iter().map(|v| -v).collect();
        let (direction, linear) = tangent::direction(
            problem,
            &evaluation.tangent,
            &dofs,
            &rhs_newton,
            settings.linear,
            settings.max_operator_applications,
            &mut report.operator_applications,
            cancelled,
        )?;
        report.linear_iterations = report
            .linear_iterations
            .checked_add(linear.iters)
            .ok_or_else(|| invalid("nonlinear work count overflow"))?;
        let (mut derivative, mut cd) = (0.0, 0.0);
        for (&r, &d) in residual.iter().zip(&direction) {
            add(&mut derivative, &mut cd, r * d)?;
        }
        let derivative = finite(derivative + cd)?;
        if derivative >= 0.0 {
            return Err(Error::Nonconvergence(report));
        }
        let mut alpha = 1.0_f64;
        let mut accepted = false;
        for backtrack in 0..=settings.max_backtracks {
            poll(cancelled)?;
            for v in 0..proposed.len() {
                proposed[v] = if let Some(i) = dofs.index[v] {
                    finite(alpha.mul_add(direction[i], potential[v]))?
                } else {
                    potential[v]
                };
            }
            match problem.evaluate(&proposed, &mut trial, cancelled) {
                Ok(()) => {
                    let delta = objective_change(
                        problem,
                        &dofs,
                        &rhs,
                        &potential,
                        &evaluation,
                        &proposed,
                        &trial,
                        cancelled,
                    )?;
                    if delta <= 1e-4 * alpha * derivative {
                        accepted = true;
                        break;
                    }
                }
                Err(Error::FieldRange { .. }) => {}
                Err(error) => return Err(error),
            }
            if backtrack == settings.max_backtracks {
                break;
            }
            alpha *= 0.5;
            report.backtracks = report
                .backtracks
                .checked_add(1)
                .ok_or_else(|| invalid("line search work count overflow"))?;
        }
        if !accepted {
            return Err(Error::Nonconvergence(report));
        }
        std::mem::swap(&mut potential, &mut proposed);
        std::mem::swap(&mut evaluation, &mut trial);
        dofs.residual(&evaluation, &rhs, &mut residual, &mut correction)?;
        report.residual_norm_c = finite(norm(&residual))?;
        report.newton_iterations += 1;
    }
    // Recompute the accepted original material response before publishing.
    problem.evaluate(&potential, &mut evaluation, cancelled)?;
    dofs.residual(&evaluation, &rhs, &mut residual, &mut correction)?;
    report.residual_norm_c = finite(norm(&residual))?;
    if report.residual_norm_c > target {
        return Err(Error::Nonconvergence(report));
    }
    let mut terminal_voltage_v = Vec::with_capacity(boundary.len());
    let mut terminal_charge_c = Vec::with_capacity(boundary.len());
    let mut nodal_boundary_charge_c = vec![None; potential.len()];
    let (mut total, mut ct) = (volume_charge_c, 0.0);
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
    let energy_identity_residual_j =
        finite(evaluation.energy + evaluation.coenergy - (pairing + cp))?;
    if charge_balance_residual_c.abs() > settings.charge_tolerance_c
        || maximum_free_charge_residual_c > settings.charge_tolerance_c
        || maximum_floating_charge_residual_c > settings.charge_tolerance_c
        || energy_identity_residual_j.abs() > settings.energy_tolerance_j
    {
        return Err(ElectrostaticError::Balance {
            charge_residual_c: charge_balance_residual_c,
            maximum_free_charge_residual_c: maximum_free_charge_residual_c
                .max(maximum_floating_charge_residual_c),
            energy_residual_j: energy_identity_residual_j,
        }
        .into());
    }
    poll(cancelled)?;
    Ok(NonlinearElectricField {
        problem,
        electrodes: electrodes.layout.clone(),
        controls: boundary.to_vec(),
        dofs,
        reference_terminal,
        reference_voltage_v,
        relative_voltage_v: potential,
        terminal_voltage_v,
        terminal_charge_c,
        nodal_boundary_charge_c,
        free_charge_density_c_m3: rho.to_vec(),
        evaluation,
        volume_charge_c,
        charge_balance_residual_c,
        maximum_free_charge_residual_c,
        maximum_floating_charge_residual_c,
        energy_identity_residual_j,
        report,
    })
}
