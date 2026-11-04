//! Prescribed body forcing on the original conservative Cartesian owners.
//!
//! A source half-step changes momentum and adds its midpoint kinetic work.
//! Density and composition are untouched. Mechanical work never subtracts
//! formation/internal energies or appears as a second heating source.
//! The surrounding owner supplies K-A-K splitting, caloric admission, CFL,
//! cancellation and publication. Acceleration and force density are distinct;
//! this source splitting is not a hydrostatic well-balanced discretization.

use super::*;
pub(super) use crate::compressible::body_force::kick;
/// Three-dimensional acceleration or force density, held over one interval.
pub type BodyForce3d<'a> = crate::compressible::body_force::BodyForce<'a, 3>;
/// Caller-authored physical Mach admission for every forced stage.
pub type ForceControl3d = crate::compressible::body_force::ForceControl;
/// Accepted three-dimensional body impulse and mechanical work.
pub type ForceChange3d = crate::compressible::body_force::ForceChange<3>;

/// Reused original-state buffer for symmetric forcing; no per-attempt allocation.
/// Numerical capacity is 40 bytes/cell, excluding fixed metadata/allocator overhead.
pub struct ForcingWorkspace3d {
    original: Vec<State>,
}

impl ForcingWorkspace3d {
    pub fn workspace_bytes(&self) -> usize {
        self.original.capacity() * core::mem::size_of::<State>()
    }

    pub(super) fn admit(&self, cells: usize) -> Result<(), Error> {
        if self.original.len() != cells {
            return Err(Error::Invalid(
                "forcing workspace and original field differ",
            ));
        }
        Ok(())
    }
}

impl Euler3d {
    pub fn required_forcing_workspace_bytes(&self) -> Result<usize, Error> {
        self.states
            .len()
            .checked_mul(core::mem::size_of::<State>())
            .ok_or(Error::Invalid("forcing workspace size overflow"))
    }

    /// Allocate only the added source-splitting buffer under an explicit byte cap.
    pub fn forcing_workspace(&self, maximum_bytes: usize) -> Result<ForcingWorkspace3d, Error> {
        if self.required_forcing_workspace_bytes()? > maximum_bytes {
            return Err(Error::Invalid("forcing workspace exceeds its byte budget"));
        }
        let work = ForcingWorkspace3d {
            original: allocated(self.states.len(), State::default())?,
        };
        if work.workspace_bytes() > maximum_bytes {
            return Err(Error::Invalid(
                "forcing workspace capacity exceeds its byte budget",
            ));
        }
        Ok(work)
    }
}

#[allow(clippy::too_many_arguments)]
pub(super) fn attempt(
    gas: &FrozenGas,
    grid: &CartesianGrid,
    boundaries: Boundaries3d,
    reservoirs: &reservoir::Reservoirs,
    states: &[State],
    work: &mut Workspace,
    forcing_work: &mut ForcingWorkspace3d,
    force: BodyForce3d<'_>,
    controls: ForceControl3d,
    dt: f64,
    cfl: f64,
    reconstruction: Reconstruction,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(StageStats, ForceChange3d), Error> {
    for (cell, (source, target)) in states.iter().zip(&mut forcing_work.original).enumerate() {
        if cell % 64 == 0 {
            poll(cancelled)?;
        }
        *target = *source;
    }
    let (first, _) = kick_field(
        gas,
        grid,
        &mut forcing_work.original,
        force,
        controls,
        0.5 * dt,
        false,
        cancelled,
    )?;
    // Original primitive recovery admits every kicked cell and every RK state.
    let mut stats = evolve::attempt(
        gas,
        grid,
        boundaries,
        reservoirs,
        &forcing_work.original,
        work,
        dt,
        cfl,
        reconstruction,
        Some(controls),
        cancelled,
    )?;
    let (second, endpoint_rate) = kick_field(
        gas,
        grid,
        &mut work.candidate,
        force,
        controls,
        0.5 * dt,
        true,
        cancelled,
    )?;
    check_cfl::<Error>(endpoint_rate * dt, cfl)?;
    stats.rate = stats.rate.max(endpoint_rate);
    let source = ForceChange3d {
        impulse: std::array::from_fn(|axis| first.impulse[axis] + second.impulse[axis]),
        work: first.work + second.work,
    };
    if !source
        .impulse
        .into_iter()
        .chain([source.work])
        .all(f64::is_finite)
    {
        return Err(Error::Invalid(
            "unrepresentable extensive body-force source",
        ));
    }
    Ok((stats, source))
}

#[allow(clippy::too_many_arguments)]
fn kick_field(
    gas: &FrozenGas,
    grid: &CartesianGrid,
    states: &mut [State],
    force: BodyForce3d<'_>,
    controls: ForceControl3d,
    duration: f64,
    admit_endpoint: bool,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(ForceChange3d, f64), Error> {
    let mut sum = [0.0; 4];
    let mut correction = [0.0; 4];
    let mut rate = 0.0_f64;
    let volume = grid.cell_volume();
    for (cell, state) in states.iter_mut().enumerate() {
        if cell % 64 == 0 {
            poll(cancelled)?;
        }
        let applied = kick::<Error, 3>(
            force,
            cell,
            state.density,
            state.momentum,
            state.energy,
            duration,
        )?;
        state.momentum = applied.momentum;
        state.energy = applied.energy;
        if admit_endpoint {
            let primitive = state.primitive(gas, &mut *cancelled)?;
            controls.admit_velocity::<Error>(primitive.velocity, primitive.sound_speed)?;
            let local: f64 = (0..3)
                .map(|axis| {
                    (primitive.velocity[axis].abs() + primitive.sound_speed) / grid.spacing[axis]
                })
                .sum();
            rate = rate.max(local);
        }
        let value = [
            applied.impulse_density[0],
            applied.impulse_density[1],
            applied.impulse_density[2],
            applied.work_density,
        ];
        // Same fixed-cell-order compensated arithmetic as the native totals.
        for d in 0..4 {
            let y = value[d] * volume - correction[d];
            let next = sum[d] + y;
            correction[d] = (next - sum[d]) - y;
            sum[d] = next;
        }
    }
    poll(cancelled)?;
    if !sum.into_iter().all(f64::is_finite) || !rate.is_finite() {
        return Err(Error::Invalid("unrepresentable body-force field balance"));
    }
    Ok((
        ForceChange3d {
            impulse: [sum[0], sum[1], sum[2]],
            work: sum[3],
        },
        rate,
    ))
}

#[cfg(test)]
mod tests;
