//! Shared unpublished K/ALE/K proposal for prescribed motion and finite pistons.

use super::*;

/// The enclosing owner admits source controls, prepared scratch and topology.
/// This repeats no source preparation and publishes no physical field.
#[allow(clippy::too_many_arguments)]
pub(in crate::compressible::spatial::mixture) fn attempt_source(
    gas: &MixtureGas,
    grid: &CartesianGrid,
    next: &CartesianGrid,
    boundaries: Boundaries3d,
    reservoirs: &reservoir::Reservoirs,
    motion: MeshMotion,
    states: &[f64],
    work: &mut Workspace,
    dt: f64,
    cfl: f64,
    reconstruction: Reconstruction,
    forcing: forcing::Context<'_, '_>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(stepping::StageStats, ForceChange3d), Error> {
    let (force, control, workspace) = forcing;
    copy_values(states, &mut workspace.kicked, cancelled)?;
    let mut source = forcing::Accumulator::default();
    forcing::kick_field(
        gas,
        grid,
        &mut workspace.kicked,
        force,
        control,
        0.5 * dt,
        &mut source,
        false,
        cancelled,
    )?;
    let mut stats = attempt_with_forcing(
        gas,
        grid,
        next,
        boundaries,
        reservoirs,
        motion,
        &workspace.kicked,
        work,
        dt,
        cfl,
        reconstruction,
        Some(control),
        cancelled,
    )?;
    forcing::kick_field(
        gas,
        next,
        &mut work.candidate,
        force,
        control,
        0.5 * dt,
        &mut source,
        false,
        cancelled,
    )?;
    let final_rate =
        forcing::admit_moving_field(gas, next, &work.candidate, control, motion, cancelled)?;
    stepping::check_cfl::<Error>(final_rate * dt, cfl)?;
    stats.rate = stats.rate.max(final_rate);
    Ok((stats, source.finish()?))
}

pub(in crate::compressible::spatial::mixture) fn account_source(
    report: &mut StepReport,
    gas: &MixtureGas,
    source: ForceChange3d,
) -> Result<(), Error> {
    report.source_change = source;
    for axis in 0..3 {
        report.conservation_residual.momentum[axis] -= source.impulse[axis];
    }
    report.conservation_residual.energy -= source.work;
    report.absolute_energy_residual_j = report.conservation_residual.absolute_energy(gas)?;
    if !report
        .conservation_residual
        .momentum
        .iter()
        .all(|v| v.is_finite())
        || !report.conservation_residual.energy.is_finite()
    {
        return Err(Error::Invalid("moving mixture forced balance overflow"));
    }
    Ok(())
}
