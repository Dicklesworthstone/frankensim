//! Shared Cartesian source flow map. Spatial transport stays with each owner.

use super::cartesian::{Failure, poll};

/// Physical source held fixed over one complete attempted interval.
/// Arrays use the original x-fast cell order, including every cell.
#[derive(Clone, Copy, Debug)]
pub enum BodyForce<'a, const D: usize> {
    /// Uniform acceleration [m/s²], converted using the current stage density.
    UniformAcceleration([f64; D]),
    /// Per-cell acceleration [m/s²], converted using current stage densities.
    Acceleration(&'a [[f64; D]]),
    /// Uniform force per volume [N/m³], e.g. an authored periodic pressure drive.
    UniformForceDensity([f64; D]),
    /// Per-cell force per volume [N/m³].
    ForceDensity(&'a [[f64; D]]),
}

impl<const D: usize> BodyForce<'_, D> {
    fn value(self, cell: usize) -> [f64; D] {
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
        poll::<E>(cancelled)?;
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
                poll::<E>(cancelled)?;
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
pub struct ForceControl {
    /// Positive finite total-speed/sound-speed bound at admitted stage states.
    pub maximum_mach: f64,
}

impl ForceControl {
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
        velocity: impl IntoIterator<Item = f64>,
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
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ForceChange<const D: usize> {
    /// Cartesian impulse [kg m/s].
    pub impulse: [f64; D],
    /// Mechanical work [J]; internal and chemical energy are not source terms.
    pub work: f64,
}

impl<const D: usize> Default for ForceChange<D> {
    fn default() -> Self {
        Self {
            impulse: [0.0; D],
            work: 0.0,
        }
    }
}

pub(super) struct Kick<const D: usize> {
    pub momentum: [f64; D],
    pub energy: f64,
    pub impulse_density: [f64; D],
    pub work_density: f64,
}

/// Constant-source flow map. Work uses the represented momentum increment:
/// an increment lost to momentum rounding cannot silently become heat.
pub(super) fn kick<E: Failure, const D: usize>(
    force: BodyForce<'_, D>,
    cell: usize,
    density: f64,
    momentum: [f64; D],
    energy: f64,
    duration: f64,
) -> Result<Kick<D>, E> {
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
        BodyForce::UniformAcceleration(_) | BodyForce::Acceleration(_)
    );
    let mut result = Kick {
        momentum,
        energy,
        impulse_density: [0.0; D],
        work_density: 0.0,
    };
    for axis in 0..D {
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
