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
    /// Held local Coriolis acceleration `-2 angular_velocity cross u`, together
    /// with uniform non-Coriolis acceleration [m/s²]. Angular velocity is a
    /// right-handed vector [rad/s], resolved in the same axes as velocity.
    ///
    /// Two-dimensional states admit only rotation normal to the xy plane;
    /// three-dimensional states admit any axis. This is a local rotating-flow
    /// source, not a coordinate transformation. Supply effective gravity or
    /// other physical acceleration explicitly; centrifugal, Euler and origin
    /// acceleration are not inserted implicitly. The native pressure/transport
    /// operator remains the sole owner of its pressure gradient.
    ///
    /// The source uses an analytic implicit midpoint solve, which preserves
    /// pure-Coriolis kinetic energy to roundoff without a rotation CFL limit.
    /// Its phase is second-order accurate, not exact: the physical inertial
    /// period must still be resolved. Zero spin follows the original
    /// acceleration path with identical arithmetic.
    Coriolis {
        /// Frame angular velocity [rad/s].
        angular_velocity: [f64; 3],
        /// Acceleration excluding Coriolis [m/s²].
        acceleration: [f64; D],
    },
    /// The same Coriolis source with a prescribed non-Coriolis acceleration
    /// for every original cell. Values are held over the attempted interval;
    /// each source half-step uses its current density. An authored centrifugal
    /// or Euler field must respect the physical domain and its boundaries.
    CoriolisField {
        /// Uniform frame angular velocity [rad/s].
        angular_velocity: [f64; 3],
        /// Non-Coriolis acceleration [m/s²], one vector per original cell.
        acceleration: &'a [[f64; D]],
    },
}

impl<const D: usize> BodyForce<'_, D> {
    fn value(self, cell: usize) -> [f64; D] {
        match self {
            Self::UniformAcceleration(value)
            | Self::UniformForceDensity(value)
            | Self::Coriolis {
                acceleration: value,
                ..
            } => value,
            Self::Acceleration(field)
            | Self::ForceDensity(field)
            | Self::CoriolisField {
                acceleration: field,
                ..
            } => field[cell],
        }
    }

    fn angular_velocity(self) -> Option<[f64; 3]> {
        match self {
            Self::Coriolis {
                angular_velocity, ..
            }
            | Self::CoriolisField {
                angular_velocity, ..
            } => Some(angular_velocity),
            _ => None,
        }
    }

    pub(super) fn admit<E: Failure>(
        self,
        cells: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<bool, E> {
        poll::<E>(cancelled)?;
        let rotating = if let Some(omega) = self.angular_velocity() {
            admit_rotation::<E, D>(omega)?;
            omega.into_iter().any(|x| x != 0.0)
        } else {
            false
        };
        let field = match self {
            Self::UniformAcceleration(value)
            | Self::UniformForceDensity(value)
            | Self::Coriolis {
                acceleration: value,
                ..
            } => {
                if !value.into_iter().all(f64::is_finite) {
                    return Err(E::invalid("body force must be finite"));
                }
                return Ok(rotating || value.into_iter().any(|x| x != 0.0));
            }
            Self::Acceleration(field)
            | Self::ForceDensity(field)
            | Self::CoriolisField {
                acceleration: field,
                ..
            } => field,
        };
        if field.len() != cells {
            return Err(E::invalid(
                "body force requires one vector per original cell",
            ));
        }
        let mut active = rotating;
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

/// Held-source flow map. Prescribed forces use their exact constant impulse;
/// Coriolis uses an analytic skew midpoint solve. Work uses the represented
/// momentum increment, except that pure Coriolis supplies exactly zero work:
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
    if let Some(omega) = force.angular_velocity() {
        admit_rotation::<E, D>(omega)?;
        if omega.into_iter().any(|x| x != 0.0) {
            return coriolis_kick::<E, D>(
                omega, value, density, momentum, energy, duration,
            );
        }
    }
    let acceleration = matches!(
        force,
        BodyForce::UniformAcceleration(_)
            | BodyForce::Acceleration(_)
            | BodyForce::Coriolis { .. }
            | BodyForce::CoriolisField { .. }
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

fn admit_rotation<E: Failure, const D: usize>(omega: [f64; 3]) -> Result<(), E> {
    if D != 2 && D != 3 {
        return Err(E::invalid(
            "Coriolis forcing requires two or three Cartesian components",
        ));
    }
    if !omega.into_iter().all(f64::is_finite) {
        return Err(E::invalid("Coriolis angular velocity must be finite"));
    }
    if D == 2 && (omega[0] != 0.0 || omega[1] != 0.0) {
        return Err(E::invalid(
            "planar Coriolis forcing requires spin normal to the xy plane",
        ));
    }
    Ok(())
}

/// Integrate dm/dt = rho*a - 2*omega cross m over one source interval.
/// With w = h*omega, the midpoint impulse solves
/// (I + [w]x) dm = h*rho*a - 2*w cross m.
/// Its inverse is (I - [w]x + w*w^T)/(1 + dot(w,w)).
///
/// Solving for the increment preserves a represented geostrophic equilibrium
/// and avoids subtracting two nearly equal midpoint momenta at small h.
/// The Coriolis part does no work. A pure rotation leaves total energy exactly
/// unchanged; with a drive, the represented kinetic change is the mechanical
/// work (as in the ordinary kick), so rounded-away impulses cannot become heat.
/// This is a local source map, not a well-balanced pressure/rotation scheme.
/// Signs follow MIT 16.07 Lecture L8, "Relative Motion using Rotating Axes".
fn coriolis_kick<E: Failure, const D: usize>(
    omega: [f64; 3],
    acceleration: [f64; D],
    density: f64,
    momentum: [f64; D],
    energy: f64,
    duration: f64,
) -> Result<Kick<D>, E> {
    let w = omega.map(|v| duration * v);
    let denominator = 1.0 + dot3(w, w);
    if !denominator.is_finite()
        || !w.into_iter().all(f64::is_finite)
        || !acceleration.into_iter().all(f64::is_finite)
    {
        return Err(E::invalid("unrepresentable Coriolis source coefficients"));
    }
    let m: [f64; 3] = core::array::from_fn(|axis| {
        if axis < D { momentum[axis] } else { 0.0 }
    });
    let drive: [f64; 3] = core::array::from_fn(|axis| {
        if axis < D {
            (duration * acceleration[axis]) * density
        } else {
            0.0
        }
    });
    let rotation = cross3(w, m);
    let rhs: [f64; 3] =
        core::array::from_fn(|axis| (-2.0_f64).mul_add(rotation[axis], drive[axis]));
    if !drive.into_iter().chain(rhs).all(f64::is_finite) {
        return Err(E::invalid("unrepresentable Coriolis source impulse"));
    }
    if rhs.into_iter().all(|value| value == 0.0) {
        return Ok(Kick {
            momentum,
            energy,
            impulse_density: [0.0; D],
            work_density: 0.0,
        });
    }
    let skew = cross3(w, rhs);
    // w dot (w cross m) is exactly zero analytically; do not introduce its
    // floating-point cancellation residue into the parallel impulse.
    let projection = dot3(w, drive);
    let increment: [f64; 3] = core::array::from_fn(|axis| {
        w[axis].mul_add(projection, rhs[axis] - skew[axis]) / denominator
    });
    if !increment.into_iter().all(f64::is_finite) {
        return Err(E::invalid("unrepresentable Coriolis source impulse"));
    }
    let driven = drive.into_iter().any(|value| value != 0.0);
    let mut result = Kick {
        momentum,
        energy,
        impulse_density: [0.0; D],
        work_density: 0.0,
    };
    for axis in 0..D {
        let next = momentum[axis] + increment[axis];
        let represented = next - momentum[axis];
        let midpoint_velocity = 0.5 * (momentum[axis] / density) + 0.5 * (next / density);
        if !next.is_finite() || !represented.is_finite() || !midpoint_velocity.is_finite() {
            return Err(E::invalid("unrepresentable Coriolis source state"));
        }
        result.momentum[axis] = next;
        result.impulse_density[axis] = represented;
        if driven {
            result.work_density += represented * midpoint_velocity;
        }
    }
    if result.work_density != 0.0 {
        result.energy += result.work_density;
    }
    if !result.work_density.is_finite() || !result.energy.is_finite() {
        return Err(E::invalid("unrepresentable Coriolis source work"));
    }
    Ok(result)
}

fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0].mul_add(b[0], a[1].mul_add(b[1], a[2] * b[2]))
}

fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1].mul_add(b[2], -a[2] * b[1]),
        a[2].mul_add(b[0], -a[0] * b[2]),
        a[0].mul_add(b[1], -a[1] * b[0]),
    ]
}
