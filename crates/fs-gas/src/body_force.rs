//! Mechanical body sources shared by the original one-dimensional gas owners.
//!
//! A mechanical source changes momentum and total energy together. Its work is
//! the represented impulse times the midpoint velocity, so the exact source flow
//! changes kinetic energy without introducing a heat or chemical-energy source.
//! Linear histories integrate over the original owner's actual physical clock,
//! including shortened retries and the two distinct source half-intervals.

/// An axial source in original cell order. The first four variants are held;
/// ramps specify a linear history at an absolute reference time in seconds.
/// A ramp is evaluated on each actual source interval, never on the originally
/// requested duration. It remains prescribed, independent of the gas state.
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
    /// Uniform acceleration `acceleration + jerk * (t - reference_time)`.
    UniformAccelerationRamp {
        /// Absolute flow-clock reference [s].
        reference_time: f64,
        /// Acceleration at the reference time [m/s²].
        acceleration: f64,
        /// Constant rate of change of acceleration [m/s³].
        jerk: f64,
    },
    /// One linear acceleration history per original finite-volume cell.
    AccelerationRamp {
        /// Absolute flow-clock reference [s].
        reference_time: f64,
        /// Acceleration at the reference time [m/s²], in cell order.
        acceleration: &'a [f64],
        /// Constant rate of change [m/s³], in the same cell order.
        jerk: &'a [f64],
    },
    /// Uniform force density `force_density + rate * (t - reference_time)`.
    UniformForceDensityRamp {
        /// Absolute flow-clock reference [s].
        reference_time: f64,
        /// Force per volume at the reference time [N/m³].
        force_density: f64,
        /// Constant rate of change of force per volume [N/(m³ s)].
        rate: f64,
    },
    /// One linear force-density history per original finite-volume cell.
    ForceDensityRamp {
        /// Absolute flow-clock reference [s].
        reference_time: f64,
        /// Force per volume at the reference time [N/m³], in cell order.
        force_density: &'a [f64],
        /// Constant rate of change [N/(m³ s)], in the same cell order.
        rate: &'a [f64],
    },
}

/// Physical admission at source endpoints and any interior velocity extremum.
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
    // (reference value, time derivative, reference time, is acceleration).
    fn coefficients(self, cell: usize) -> (f64, f64, f64, bool) {
        match self {
            Self::UniformAcceleration(value) => (value, 0.0, 0.0, true),
            Self::Acceleration(values) => (values[cell], 0.0, 0.0, true),
            Self::UniformForceDensity(value) => (value, 0.0, 0.0, false),
            Self::ForceDensity(values) => (values[cell], 0.0, 0.0, false),
            Self::UniformAccelerationRamp {
                reference_time,
                acceleration,
                jerk,
            } => (acceleration, jerk, reference_time, true),
            Self::AccelerationRamp {
                reference_time,
                acceleration,
                jerk,
            } => (acceleration[cell], jerk[cell], reference_time, true),
            Self::UniformForceDensityRamp {
                reference_time,
                force_density,
                rate,
            } => (force_density, rate, reference_time, false),
            Self::ForceDensityRamp {
                reference_time,
                force_density,
                rate,
            } => (force_density[cell], rate[cell], reference_time, false),
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
        let (values, rates) = match self {
            Self::UniformAcceleration(value) | Self::UniformForceDensity(value) => {
                if !value.is_finite() {
                    return Err(E::invalid("body force must be finite"));
                }
                return Ok(value != 0.0);
            }
            Self::Acceleration(values) | Self::ForceDensity(values) => (values, None),
            Self::UniformAccelerationRamp {
                reference_time,
                acceleration: value,
                jerk: rate,
            }
            | Self::UniformForceDensityRamp {
                reference_time,
                force_density: value,
                rate,
            } => {
                if ![reference_time, value, rate].iter().all(|x| x.is_finite()) {
                    return Err(E::invalid("body-force history must be finite"));
                }
                return Ok(value != 0.0 || rate != 0.0);
            }
            Self::AccelerationRamp {
                reference_time,
                acceleration: values,
                jerk: rates,
            }
            | Self::ForceDensityRamp {
                reference_time,
                force_density: values,
                rate: rates,
            } => {
                if !reference_time.is_finite() || rates.len() != cells {
                    return Err(E::invalid(
                        "body-force history requires a finite clock and one rate per original cell",
                    ));
                }
                (values, Some(rates))
            }
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
            let rate = rates.map_or(0.0, |rates| rates[cell]);
            if !value.is_finite() || !rate.is_finite() {
                return Err(E::invalid("body force must be finite"));
            }
            active |= value != 0.0 || rate != 0.0;
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
    /// Represented momentum and energy at a strict interior velocity extremum.
    /// Owners decode this through their original EOS and apply Mach/CFL gates.
    pub(crate) interior: Option<[f64; 2]>,
    impulse: f64,
    work: f64,
    energy_residual: f64,
}

/// Exact prescribed-source map before floating-point representation. Use the
/// actual momentum increment: an impulse that rounds away cannot become heat.
/// Time is [accepted step start, offset within that step, source duration].
/// Keeping the offset separate avoids losing a half-step at a large clock.
pub(crate) fn kick<E: Failure>(
    source: BodyForce<'_>,
    cell: usize,
    density: f64,
    momentum: f64,
    energy: f64,
    [time, offset, duration]: [f64; 3],
) -> Result<Kick, E> {
    if !density.is_finite()
        || density <= 0.0
        || !time.is_finite()
        || !offset.is_finite()
        || offset < 0.0
        || !duration.is_finite()
        || duration <= 0.0
    {
        return Err(E::invalid("unrepresentable body-force interval or density"));
    }
    let (reference, rate, reference_time, acceleration) = source.coefficients(cell);
    // Retain the held-source arithmetic exactly, including zero-rate histories.
    let value = if rate == 0.0 {
        reference
    } else {
        rate.mul_add(offset, rate.mul_add(time - reference_time, reference))
    };
    let mean = if rate == 0.0 {
        value
    } else {
        rate.mul_add(0.5 * duration, value)
    };
    let end = if rate == 0.0 {
        value
    } else {
        rate.mul_add(duration, value)
    };
    if ![value, mean, end].iter().all(|x| x.is_finite()) {
        return Err(E::invalid("unrepresentable body-force history"));
    }
    let integral = duration * mean;
    let authored_impulse = if acceleration {
        integral * density
    } else {
        integral
    };
    let mut result = represented::<E>(density, momentum, energy, authored_impulse)?;
    if rate != 0.0 {
        let extremum_time = -value / rate;
        if extremum_time > 0.0 && extremum_time < duration {
            let integral = (0.5 * value) * extremum_time;
            let impulse = if acceleration {
                integral * density
            } else {
                integral
            };
            let interior = represented::<E>(density, momentum, energy, impulse)?;
            result.interior = Some([interior.momentum, interior.energy]);
        }
    }
    Ok(result)
}

fn represented<E: Failure>(
    density: f64,
    momentum: f64,
    energy: f64,
    authored_impulse: f64,
) -> Result<Kick, E> {
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
    if ![
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
        interior: None,
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
