//! Mechanical body sources shared by the original one-dimensional gas owners.
//!
//! A held force changes momentum and total energy together. Its work is the
//! represented impulse times the midpoint velocity, so the exact source flow
//! changes kinetic energy without introducing a heat or chemical-energy source.

/// An axial source held over the complete attempted interval, in cell order.
#[derive(Clone, Copy, Debug)]
pub enum BodyForce<'a> {
    /// Uniform acceleration [m/s²], multiplied by each source stage's density.
    UniformAcceleration(f64),
    /// One acceleration [m/s²] per original finite-volume cell.
    Acceleration(&'a [f64]),
    /// Uniform force per volume [N/m³], independent of the evolving density.
    UniformForceDensity(f64),
    /// One force per volume [N/m³] per original finite-volume cell.
    ForceDensity(&'a [f64]),
}

/// Physical admission for the states entering and leaving each force substep.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ForceControl {
    /// Positive finite bound on axial speed divided by the original sound speed.
    /// This is caller-authored, not an imposed low-Mach approximation.
    pub maximum_mach: f64,
}

/// Extensive source-only exchange, separate from pressure, gravity and wall work.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ForceChange {
    /// Actual represented axial impulse applied to the gas [N s].
    pub impulse: f64,
    /// Conjugate mechanical work, positive into the gas [J].
    pub work: f64,
    /// Sum of actual total-energy increments minus conjugate source work [J].
    /// This exposes representation error; it is never deposited as heat or
    /// subtracted from the full step's conservation residual.
    pub energy_residual: f64,
}

pub(crate) trait Failure {
    fn invalid(reason: &'static str) -> Self;
    fn cancelled() -> Self;
}
impl Failure for crate::dynamics::Error {
    fn invalid(reason: &'static str) -> Self {
        Self::Invalid(reason)
    }
    fn cancelled() -> Self {
        Self::Cancelled
    }
}
impl Failure for crate::multicomponent::Error {
    fn invalid(reason: &'static str) -> Self {
        Self::Invalid(reason)
    }
    fn cancelled() -> Self {
        Self::Cancelled
    }
}

impl BodyForce<'_> {
    fn value(self, cell: usize) -> f64 {
        match self {
            Self::UniformAcceleration(value) | Self::UniformForceDensity(value) => value,
            Self::Acceleration(values) | Self::ForceDensity(values) => values[cell],
        }
    }

    pub(crate) fn admit<E: Failure>(
        self,
        cells: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<bool, E> {
        if cancelled() {
            return Err(E::cancelled());
        }
        let values = match self {
            Self::UniformAcceleration(value) | Self::UniformForceDensity(value) => {
                if !value.is_finite() {
                    return Err(E::invalid("body force must be finite"));
                }
                return Ok(value != 0.0);
            }
            Self::Acceleration(values) | Self::ForceDensity(values) => values,
        };
        if values.len() != cells {
            return Err(E::invalid(
                "body force requires one value per original cell",
            ));
        }
        let mut active = false;
        for (cell, &value) in values.iter().enumerate() {
            if cell % 64 == 0 && cancelled() {
                return Err(E::cancelled());
            }
            if !value.is_finite() {
                return Err(E::invalid("body force must be finite"));
            }
            active |= value != 0.0;
        }
        Ok(active)
    }
}

impl ForceControl {
    pub(crate) fn admit<E: Failure>(self) -> Result<(), E> {
        if !self.maximum_mach.is_finite() || self.maximum_mach <= 0.0 {
            return Err(E::invalid(
                "body-force Mach bound must be positive and finite",
            ));
        }
        Ok(())
    }

    pub(crate) fn admit_velocity<E: Failure>(
        self,
        velocity: f64,
        sound_speed: f64,
    ) -> Result<(), E> {
        let mach = velocity.abs() / sound_speed;
        if !sound_speed.is_finite()
            || sound_speed <= 0.0
            || !mach.is_finite()
            || mach > self.maximum_mach
        {
            return Err(E::invalid(
                "body-force state exceeds its authored Mach bound",
            ));
        }
        Ok(())
    }
}

pub(crate) struct Kick {
    pub(crate) momentum: f64,
    pub(crate) energy: f64,
    impulse: f64,
    work: f64,
    energy_residual: f64,
}

/// Exact held-source map before floating-point representation. Use the actual
/// momentum increment: a nominal impulse that rounds away cannot become heat.
pub(crate) fn kick<E: Failure>(
    source: BodyForce<'_>,
    cell: usize,
    density: f64,
    momentum: f64,
    energy: f64,
    duration: f64,
) -> Result<Kick, E> {
    let value = source.value(cell);
    let authored_impulse = match source {
        BodyForce::UniformAcceleration(_) | BodyForce::Acceleration(_) => {
            (duration * value) * density
        }
        BodyForce::UniformForceDensity(_) | BodyForce::ForceDensity(_) => duration * value,
    };
    let next = if authored_impulse == 0.0 {
        momentum
    } else {
        momentum + authored_impulse
    };
    let impulse = next - momentum;
    // Halving before summing avoids overflowing two individually finite speeds.
    let midpoint = 0.5 * (momentum / density) + 0.5 * (next / density);
    let work = impulse * midpoint;
    let final_energy = if work == 0.0 { energy } else { energy + work };
    let energy_residual = (final_energy - energy) - work;
    if !density.is_finite()
        || density <= 0.0
        || !duration.is_finite()
        || duration <= 0.0
        || ![
            authored_impulse,
            next,
            impulse,
            midpoint,
            work,
            final_energy,
            energy_residual,
        ]
        .iter()
        .all(|x| x.is_finite())
    {
        return Err(E::invalid("unrepresentable body-force impulse or work"));
    }
    Ok(Kick {
        momentum: next,
        energy: final_energy,
        impulse,
        work,
        energy_residual,
    })
}

/// Fixed-order compensated extensive sums span both source substeps.
#[derive(Default)]
pub(crate) struct Accumulator {
    sum: [f64; 3],
    correction: [f64; 3],
}
impl Accumulator {
    pub(crate) fn add<E: Failure>(&mut self, kick: &Kick, volume: f64) -> Result<(), E> {
        for (i, density) in [kick.impulse, kick.work, kick.energy_residual]
            .into_iter()
            .enumerate()
        {
            let value = density * volume;
            let y = value - self.correction[i];
            let next = self.sum[i] + y;
            self.correction[i] = (next - self.sum[i]) - y;
            self.sum[i] = next;
            if !value.is_finite() || !next.is_finite() || !self.correction[i].is_finite() {
                return Err(E::invalid("unrepresentable extensive body-force exchange"));
            }
        }
        Ok(())
    }

    pub(crate) fn finish(self) -> ForceChange {
        ForceChange {
            impulse: self.sum[0],
            work: self.sum[1],
            energy_residual: self.sum[2],
        }
    }
}
