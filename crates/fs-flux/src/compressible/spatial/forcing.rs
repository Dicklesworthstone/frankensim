//! Prescribed body forcing on the original conservative Cartesian owners.
//!
//! A source half-step changes momentum and adds its midpoint kinetic work.
//! Density and composition are untouched. Mechanical work never subtracts
//! formation/internal energies or appears as a second heating source.
//! The surrounding owner supplies K-A-K splitting, caloric admission, CFL,
//! cancellation and publication. Acceleration and force density are distinct;
//! this source splitting is not a hydrostatic well-balanced discretization.

use super::*;
use crate::compressible::cartesian::Failure;

/// Physical source held fixed over one complete attempted interval.
/// Arrays use the original x-fast cell order, including every cell.
#[derive(Clone, Copy, Debug)]
pub enum BodyForce3d<'a> {
    /// Uniform acceleration [m/s²], converted using the current stage density.
    UniformAcceleration([f64; 3]),
    /// Per-cell acceleration [m/s²], converted using current stage densities.
    Acceleration(&'a [[f64; 3]]),
    /// Uniform force per volume [N/m³], e.g. an authored periodic pressure drive.
    UniformForceDensity([f64; 3]),
    /// Per-cell force per volume [N/m³].
    ForceDensity(&'a [[f64; 3]]),
}

impl BodyForce3d<'_> {
    fn value(self, cell: usize) -> [f64; 3] {
        match self {
            Self::UniformAcceleration(value) | Self::UniformForceDensity(value) => value,
            Self::Acceleration(field) | Self::ForceDensity(field) => field[cell],
        }
    }

    pub(super) fn admit<E: Failure>(
        self,
        cells: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<bool, E> {
        crate::compressible::cartesian::poll::<E>(cancelled)?;
        let field = match self {
            Self::UniformAcceleration(value) | Self::UniformForceDensity(value) => {
                if !value.into_iter().all(f64::is_finite) {
                    return Err(E::invalid("body force must be finite"));
                }
                return Ok(value.into_iter().any(|x| x != 0.0));
            }
            Self::Acceleration(field) | Self::ForceDensity(field) => field,
        };
        if field.len() != cells {
            return Err(E::invalid(
                "body force requires one vector per original cell",
            ));
        }
        let mut active = false;
        for (cell, value) in field.iter().enumerate() {
            if cell % 64 == 0 {
                crate::compressible::cartesian::poll::<E>(cancelled)?;
            }
            if !value.iter().all(|x| x.is_finite()) {
                return Err(E::invalid("body force must be finite"));
            }
            active |= value.iter().any(|&x| x != 0.0);
        }
        Ok(active)
    }
}

/// Caller-authored physical admission; no low-Mach limit is imposed on Euler.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ForceControl3d {
    /// Positive finite total-speed/sound-speed bound at admitted stage states.
    pub maximum_mach: f64,
}

impl ForceControl3d {
    pub(super) fn admit<E: Failure>(self) -> Result<(), E> {
        if !self.maximum_mach.is_finite() || self.maximum_mach <= 0.0 {
            return Err(E::invalid(
                "body-force Mach bound must be positive and finite",
            ));
        }
        Ok(())
    }

    pub(super) fn admit_velocity<E: Failure>(
        self,
        velocity: [f64; 3],
        sound_speed: f64,
    ) -> Result<(), E> {
        if !sound_speed.is_finite() || sound_speed <= 0.0 {
            return Err(E::invalid(
                "body-force sound speed must be positive and finite",
            ));
        }
        let mach = velocity
            .into_iter()
            .fold(0.0_f64, |norm, u| norm.hypot(u / sound_speed));
        if !mach.is_finite() || mach > self.maximum_mach {
            return Err(E::invalid(
                "body-force stage exceeds its authored Mach bound",
            ));
        }
        Ok(())
    }
}

/// Extensive body source over an accepted physical interval, separate from walls.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ForceChange3d {
    /// Three-axis impulse [kg m/s].
    pub impulse: [f64; 3],
    /// Mechanical work [J]; internal and chemical energy are not source terms.
    pub work: f64,
}

pub(super) struct Kick3d {
    pub momentum: [f64; 3],
    pub energy: f64,
    pub impulse_density: [f64; 3],
    pub work_density: f64,
}

/// Constant-source flow map. Work uses the represented momentum increment:
/// an increment lost to momentum rounding cannot silently become heat.
pub(super) fn kick<E: Failure>(
    force: BodyForce3d<'_>,
    cell: usize,
    density: f64,
    momentum: [f64; 3],
    energy: f64,
    duration: f64,
) -> Result<Kick3d, E> {
    if !density.is_finite()
        || density <= 0.0
        || !energy.is_finite()
        || !momentum.into_iter().all(f64::is_finite)
        || !duration.is_finite()
        || duration <= 0.0
    {
        return Err(E::invalid("invalid conservative body-force source state"));
    }
    let value = force.value(cell);
    let acceleration = matches!(
        force,
        BodyForce3d::UniformAcceleration(_) | BodyForce3d::Acceleration(_)
    );
    let mut result = Kick3d {
        momentum,
        energy,
        impulse_density: [0.0; 3],
        work_density: 0.0,
    };
    for axis in 0..3 {
        if value[axis] == 0.0 {
            continue;
        }
        let impulse = if acceleration {
            (duration * value[axis]) * density
        } else {
            duration * value[axis]
        };
        let next = momentum[axis] + impulse;
        let represented = next - momentum[axis];
        let midpoint_velocity = 0.5 * (momentum[axis] / density) + 0.5 * (next / density);
        let work = represented * midpoint_velocity;
        if !impulse.is_finite()
            || !next.is_finite()
            || !represented.is_finite()
            || !midpoint_velocity.is_finite()
            || !work.is_finite()
        {
            return Err(E::invalid(
                "unrepresentable conservative body-force impulse/work",
            ));
        }
        result.momentum[axis] = next;
        result.impulse_density[axis] = represented;
        result.work_density += work;
    }
    if result.work_density != 0.0 {
        result.energy += result.work_density;
    }
    if !result.work_density.is_finite() || !result.energy.is_finite() {
        return Err(E::invalid("unrepresentable conservative body-force energy"));
    }
    Ok(result)
}

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
        let applied = kick::<Error>(
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
