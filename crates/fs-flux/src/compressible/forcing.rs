//! Conservative prescribed acceleration on the original planar Euler owner.
//!
//! Source semantics and kinetic work share the same dimension-generic kernel
//! as the volumetric frozen and multicomponent gas owners. Transport, retry,
//! source calorics and boundary exchange retain their existing implementations.

use super::body_force::kick;
use super::*;

/// Two-dimensional acceleration or force density, held over one interval.
pub type BodyForce2d<'a> = super::body_force::BodyForce<'a, 2>;
/// Caller-authored physical Mach admission for every forced stage.
pub type ForceControl2d = super::body_force::ForceControl;
/// Accepted two-dimensional body impulse and mechanical work, in SI units.
pub type ForceChange2d = super::body_force::ForceChange<2>;

/// Reusable source-splitting buffer. Stores one original conservative state
/// (32 nominal bytes) per cell, excluding fixed metadata and allocator overhead.
pub struct ForcingWorkspace2d {
    original: Vec<State>,
}

impl ForcingWorkspace2d {
    /// Actual owned numerical capacity in bytes.
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

impl Euler2d {
    /// Minimum numerical storage for the added source-splitting buffer.
    pub fn required_forcing_workspace_bytes(&self) -> Result<usize, Error> {
        self.states
            .len()
            .checked_mul(core::mem::size_of::<State>())
            .ok_or(Error::Invalid("forcing workspace size overflow"))
    }

    /// Allocate the source buffer once, admitting its actual capacity as well
    /// as the requested size. Steps and retries then allocate no source arrays.
    pub fn forcing_workspace(&self, maximum_bytes: usize) -> Result<ForcingWorkspace2d, Error> {
        if self.required_forcing_workspace_bytes()? > maximum_bytes {
            return Err(Error::Invalid("forcing workspace exceeds its byte budget"));
        }
        let workspace = ForcingWorkspace2d {
            original: allocated(self.states.len(), State::default())?,
        };
        if workspace.workspace_bytes() > maximum_bytes {
            return Err(Error::Invalid(
                "forcing workspace capacity exceeds its byte budget",
            ));
        }
        Ok(workspace)
    }
}

#[allow(clippy::too_many_arguments)]
pub(super) fn attempt(
    gas: &FrozenGas,
    grid: &CartesianGrid,
    boundaries: Boundaries2d,
    states: &[State],
    work: &mut Workspace,
    forcing_work: &mut ForcingWorkspace2d,
    force: BodyForce2d<'_>,
    controls: ForceControl2d,
    dt: f64,
    cfl: f64,
    reconstruction: Reconstruction,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(StageStats, ForceChange2d), Error> {
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
    let mut stats = super::attempt(
        gas,
        grid,
        boundaries,
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
    cartesian::check_cfl::<Error>(endpoint_rate * dt, cfl)?;
    stats.rate = stats.rate.max(endpoint_rate);
    let source = ForceChange2d {
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
    force: BodyForce2d<'_>,
    controls: ForceControl2d,
    duration: f64,
    admit_endpoint: bool,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(ForceChange2d, f64), Error> {
    let mut sum = [0.0; 3];
    let mut correction = [0.0; 3];
    let mut rate = 0.0_f64;
    let volume = grid.cell_volume();
    for (cell, state) in states.iter_mut().enumerate() {
        if cell % 64 == 0 {
            poll(cancelled)?;
        }
        let applied = kick::<Error, 2>(
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
            let local: f64 = (0..2)
                .map(|axis| {
                    (primitive.velocity[axis].abs() + primitive.sound_speed) / grid.spacing[axis]
                })
                .sum();
            rate = rate.max(local);
        }
        let value = [
            applied.impulse_density[0],
            applied.impulse_density[1],
            applied.work_density,
        ];
        // Fixed-cell-order compensated arithmetic matches the original totals.
        for d in 0..3 {
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
        ForceChange2d {
            impulse: [sum[0], sum[1]],
            work: sum[2],
        },
        rate,
    ))
}

#[cfg(test)]
mod tests;
