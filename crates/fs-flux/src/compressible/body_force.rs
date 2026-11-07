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
    /// Spatially varying local Coriolis coefficients and effective acceleration.
    ///
    /// Each borrowed slice contains either ONE uniform vector or one vector
    /// per original cell in x-fast order. The two shapes are independent, so
    /// uniform gravity does not require another cell-sized allocation. Fields
    /// remain fixed throughout the attempted interval, including retries.
    ///
    /// A traditional beta plane uses `angular_velocity = [0, 0, f(y)/2]`
    /// with `f(y) = f0 + beta*y`; the factor two belongs to the Coriolis
    /// equation. Both signs and a zero equatorial value are admitted. Spatial
    /// states also admit a varying local axis. These are coefficients of a
    /// Cartesian source law, not a deforming reference frame or a spherical
    /// geometry model. Authored fields must respect the domain's boundaries.
    ///
    /// Each cell retains the same work-free midpoint Coriolis update and its
    /// current stage density. No coefficient arrays are copied or allocated
    /// during source admission, stepping or retry.
    CoriolisSpatial {
        /// Local angular velocity coefficients [rad/s], uniform or cellwise.
        angular_velocity: &'a [[f64; 3]],
        /// Non-Coriolis acceleration [m/s²], uniform or cellwise.
        acceleration: &'a [[f64; D]],
    },
    /// Exact held-coefficient Coriolis flow, including the integrated physical
    /// acceleration. Each slice independently broadcasts one vector or supplies
    /// one per original cell, with the same units and geometry as
    /// [`Self::CoriolisSpatial`].
    ///
    /// The analytic rotation removes the midpoint source's inertial-phase
    /// error. Pure rotation supplies zero work; driven work is the represented
    /// kinetic-energy change. No source subcycling, coefficient copies or
    /// per-cell allocations are required. Zero spin retains the original
    /// acceleration arithmetic. The complete transport/source split remains
    /// second order, so accuracy of spatial transport and changing coefficients
    /// still requires time refinement.
    CoriolisExact {
        /// Local angular velocity [rad/s], uniform or cellwise.
        angular_velocity: &'a [[f64; 3]],
        /// Non-Coriolis acceleration [m/s²], uniform or cellwise.
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
            Self::CoriolisSpatial { acceleration, .. }
            | Self::CoriolisExact { acceleration, .. } => {
                acceleration[if acceleration.len() == 1 { 0 } else { cell }]
            }
        }
    }

    fn angular_velocity(self, cell: usize) -> Option<[f64; 3]> {
        match self {
            Self::Coriolis {
                angular_velocity, ..
            }
            | Self::CoriolisField {
                angular_velocity, ..
            } => Some(angular_velocity),
            Self::CoriolisSpatial {
                angular_velocity, ..
            }
            | Self::CoriolisExact {
                angular_velocity, ..
            } => Some(angular_velocity[if angular_velocity.len() == 1 { 0 } else { cell }]),
            _ => None,
        }
    }

    pub(super) fn admit<E: Failure>(
        self,
        cells: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<bool, E> {
        poll::<E>(cancelled)?;
        if let Self::CoriolisSpatial {
            angular_velocity,
            acceleration,
        }
        | Self::CoriolisExact {
            angular_velocity,
            acceleration,
        } = self
        {
            // Validate both shapes before indexing either field. Inspect each
            // authored vector once; broadcast values do not cost a cell loop.
            if (angular_velocity.len() != 1 && angular_velocity.len() != cells)
                || (acceleration.len() != 1 && acceleration.len() != cells)
            {
                return Err(E::invalid(
                    "spatial Coriolis fields require one vector or one per original cell",
                ));
            }
            let mut active = false;
            for (cell, &omega) in angular_velocity.iter().enumerate() {
                if cell % 64 == 0 {
                    poll::<E>(cancelled)?;
                }
                admit_rotation::<E, D>(omega)?;
                active |= omega.into_iter().any(|x| x != 0.0);
            }
            for (cell, value) in acceleration.iter().enumerate() {
                if cell % 64 == 0 {
                    poll::<E>(cancelled)?;
                }
                if !value.iter().all(|x| x.is_finite()) {
                    return Err(E::invalid("body force must be finite"));
                }
                active |= value.iter().any(|&x| x != 0.0);
            }
            return Ok(active);
        }
        let rotating = if let Some(omega) = self.angular_velocity(0) {
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
            Self::CoriolisSpatial { .. } | Self::CoriolisExact { .. } => {
                unreachable!("spatial fields admitted above")
            }
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
/// Coriolis uses the selected midpoint or exact analytic flow. Work uses the
/// represented momentum increment, except that pure Coriolis supplies zero work:
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
    if let Some(omega) = force.angular_velocity(cell) {
        admit_rotation::<E, D>(omega)?;
        if omega.into_iter().any(|x| x != 0.0) {
            if matches!(force, BodyForce::CoriolisExact { .. }) {
                return exact_coriolis_kick::<E, D>(
                    omega, value, density, momentum, energy, duration,
                );
            }
            return coriolis_kick::<E, D>(omega, value, density, momentum, energy, duration);
        }
    }
    let acceleration = matches!(
        force,
        BodyForce::UniformAcceleration(_)
            | BodyForce::Acceleration(_)
            | BodyForce::Coriolis { .. }
            | BodyForce::CoriolisField { .. }
            | BodyForce::CoriolisSpatial { .. }
            | BodyForce::CoriolisExact { .. }
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
    let m: [f64; 3] = core::array::from_fn(|axis| if axis < D { momentum[axis] } else { 0.0 });
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
    let driven = drive.into_iter().any(|value| value != 0.0);
    finish_coriolis_kick::<E, D>(density, momentum, energy, increment, driven)
}

/// Integrate the held affine rotation with phi_1(-2 h [omega]x).
/// Acting on the initial NET source impulse retains a represented geostrophic
/// equilibrium without subtracting two large independently rotated states.
fn exact_coriolis_kick<E: Failure, const D: usize>(
    omega: [f64; 3],
    acceleration: [f64; D],
    density: f64,
    momentum: [f64; D],
    energy: f64,
    duration: f64,
) -> Result<Kick<D>, E> {
    let w = omega.map(|v| duration * v);
    let magnitude = w[0].hypot(w[1]).hypot(w[2]);
    let angle = 2.0 * magnitude;
    if !angle.is_finite() || !acceleration.into_iter().all(f64::is_finite) {
        return Err(E::invalid(
            "unrepresentable exact Coriolis source coefficients",
        ));
    }
    if magnitude == 0.0 {
        // A representable spin may round away only after multiplication by h.
        // Its limiting flow is the ordinary acceleration kick.
        return kick::<E, D>(
            BodyForce::UniformAcceleration(acceleration),
            0,
            density,
            momentum,
            energy,
            duration,
        );
    }
    let axis = w.map(|v| v / magnitude);
    let m = core::array::from_fn(|a| if a < D { momentum[a] } else { 0.0 });
    let drive: [f64; 3] = core::array::from_fn(|a| {
        if a < D {
            (duration * acceleration[a]) * density
        } else {
            0.0
        }
    });
    let rotation = cross3(w, m);
    let rhs: [f64; 3] = core::array::from_fn(|a| (-2.0_f64).mul_add(rotation[a], drive[a]));
    if !drive.into_iter().chain(rhs).all(f64::is_finite) {
        return Err(E::invalid("unrepresentable exact Coriolis source impulse"));
    }
    if rhs.into_iter().all(|value| value == 0.0) {
        return Ok(Kick {
            momentum,
            energy,
            impulse_density: [0.0; D],
            work_density: 0.0,
        });
    }
    let (sinc, cosc) = if angle < 1.0e-4 {
        let square = angle * angle;
        (
            1.0 + square * (-1.0 / 6.0 + square / 120.0),
            angle * (0.5 + square * (-1.0 / 24.0 + square / 720.0)),
        )
    } else {
        // Half-angle identities give both coefficients with one sin_cos and
        // no cancellation in 1-cos(angle), including near complete turns.
        let (sine, cosine) = magnitude.sin_cos();
        (sine * cosine / magnitude, sine * sine / magnitude)
    };
    // Rotation cannot supply axial impulse: use only the physical drive to
    // avoid a spurious axial component from dot(axis, cross(w, m)) roundoff.
    let projection = dot3(axis, drive);
    let skew = cross3(axis, rhs);
    let increment = core::array::from_fn(|a| {
        let parallel = axis[a] * projection;
        sinc.mul_add(rhs[a] - parallel, parallel - cosc * skew[a])
    });
    finish_coriolis_kick::<E, D>(
        density,
        momentum,
        energy,
        increment,
        drive.into_iter().any(|v| v != 0.0),
    )
}

fn finish_coriolis_kick<E: Failure, const D: usize>(
    density: f64,
    momentum: [f64; D],
    energy: f64,
    increment: [f64; 3],
    driven: bool,
) -> Result<Kick<D>, E> {
    if !increment.into_iter().all(f64::is_finite) {
        return Err(E::invalid("unrepresentable Coriolis source impulse"));
    }
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

#[cfg(test)]
mod exact_tests;
