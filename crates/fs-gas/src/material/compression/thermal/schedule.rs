//! Prescribed thermal forcing during source-owned reacting volume strokes.
//!
//! Bath temperature and total conductance are linear between authored physical
//! time knots. Every original interval holds its own midpoint values; full and
//! half-step probes therefore see different forcing when the program changes.
//! No probe crosses a forcing knot. The same native chemistry/heat and volume
//! operators own all thermodynamics, work, entropy and source admission.
//!
//! A supplied finite body changes only its conductance program. Its temperature,
//! latent inventory and density still follow its exact original enthalpy chart.
//! These programs do not infer moving area, convection, wall gradients or the
//! external energy needed to impose a bath-temperature trajectory.

use super::super::{CompressionInterval, FixedCompressionPath, MAX_COMPRESSION_SUBSTEPS, poll};
use super::adaptive::{
    self, AdaptiveMaterialThermalCompressionBodyStep, AdaptiveMaterialThermalCompressionError,
    AdaptiveMaterialThermalCompressionSolve, AdaptiveMaterialThermalCompressionStep, IntervalTime,
    Partner, Totals,
};
use super::{
    ClosedGasParcel, MaterialCompressionError, MaterialThermalCompressionBodyStep,
    MaterialThermalCompressionReactor, MaterialThermalCompressionSolve,
    MaterialThermalCompressionStep, ThermalInterval,
};
use crate::chemistry::heat::{HeatBath, ThermalConductance, material::FiniteHeatBody};
use crate::material::MaterialKineticsError;
use crate::material::isochoric::MaterialIsochoricReactor;
use fs_qty::{Temperature, Time, Volume};

type Error = MaterialCompressionError;

/// One admitted held-bath value at an exact elapsed physical time.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HeatBathKnot {
    pub elapsed_time: Time,
    pub bath: HeatBath,
}

/// One prescribed total conductance at an exact elapsed physical time.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ConductanceKnot {
    pub elapsed_time: Time,
    pub conductance: ThermalConductance,
}

/// Immutable, bounded linear bath-temperature and conductance history.
///
/// The first knot is at zero; subsequent times strictly increase. A single
/// zero-time knot represents a no-op, otherwise the last knot sets duration.
/// There is no extrapolation. At most 4097 knots are admitted, and original
/// bath admission supplies positive temperatures and nonnegative conductances.
#[derive(Clone, Copy, Debug)]
pub struct HeatBathSchedule<'s> {
    knots: &'s [HeatBathKnot],
    constant: bool,
}

impl<'s> HeatBathSchedule<'s> {
    pub fn new(
        knots: &'s [HeatBathKnot],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        let constant = admit(
            knots,
            |k| k.elapsed_time.value(),
            |k| {
                k.bath.temperature().value().is_finite()
                    && k.bath.temperature().value() > 0.0
                    && k.bath.conductance().value().is_finite()
                    && k.bath.conductance().value() >= 0.0
            },
            |a, b| a.bath == b.bath,
            &mut cancelled,
        )?;
        Ok(Self { knots, constant })
    }

    pub fn duration(self) -> Time {
        self.knots[self.knots.len() - 1].elapsed_time
    }

    fn at(self, time_s: f64) -> Result<HeatBath, Error> {
        let (lower, upper) = bracket(self.knots, |k| k.elapsed_time.value(), time_s)?;
        if lower == upper {
            return Ok(self.knots[lower].bath);
        }
        let a = self.knots[lower];
        let b = self.knots[upper];
        let interpolate = |initial, terminal| {
            linear(
                a.elapsed_time.value(),
                b.elapsed_time.value(),
                time_s,
                initial,
                terminal,
            )
        };
        Ok(HeatBath::new(
            ThermalConductance::new(interpolate(
                a.bath.conductance().value(),
                b.bath.conductance().value(),
            )),
            Temperature::new(interpolate(
                a.bath.temperature().value(),
                b.bath.temperature().value(),
            )),
        )
        .map_err(MaterialKineticsError::from)?)
    }
}

/// Immutable linear conductance history for ONE original finite body.
///
/// Times have the same admission as [`HeatBathSchedule`]. The first conductance
/// must match the supplied body's represented value exactly. Returned bodies
/// carry the exact final conductance, not the last midpoint's approximation.
#[derive(Clone, Copy, Debug)]
pub struct ConductanceSchedule<'s> {
    knots: &'s [ConductanceKnot],
    constant: bool,
}

impl<'s> ConductanceSchedule<'s> {
    pub fn new(
        knots: &'s [ConductanceKnot],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        let constant = admit(
            knots,
            |k| k.elapsed_time.value(),
            |k| k.conductance.value().is_finite() && k.conductance.value() >= 0.0,
            // Preserve the exact endpoint representation used by the next
            // body's initial-coefficient admission, including signed zero.
            |a, b| a.conductance.value().to_bits() == b.conductance.value().to_bits(),
            &mut cancelled,
        )?;
        Ok(Self { knots, constant })
    }

    pub fn duration(self) -> Time {
        self.knots[self.knots.len() - 1].elapsed_time
    }

    fn bind(self, body: FiniteHeatBody<'_>) -> Result<(), Error> {
        if self.knots[0].conductance.value().to_bits() != body.conductance().value().to_bits() {
            return Err(Error::Invalid(
                "conductance program must start at the exact original body conductance",
            ));
        }
        Ok(())
    }

    fn at(self, time_s: f64) -> Result<ThermalConductance, Error> {
        let (lower, upper) = bracket(self.knots, |k| k.elapsed_time.value(), time_s)?;
        let a = self.knots[lower];
        if lower == upper {
            return Ok(a.conductance);
        }
        let b = self.knots[upper];
        Ok(ThermalConductance::new(linear(
            a.elapsed_time.value(),
            b.elapsed_time.value(),
            time_s,
            a.conductance.value(),
            b.conductance.value(),
        )))
    }
}

impl MaterialThermalCompressionReactor {
    /// Resolve a prescribed bath program on the original fixed interval kernel.
    ///
    /// The uniform grid requests `controls.substeps` over the complete stroke;
    /// mandatory forcing knots split that grid further. The complete union is
    /// bounded to 4096 original intervals, including fixed-volume thermal ramps.
    /// Constant programs retain the existing fixed endpoint path exactly.
    pub fn advance_with_heat_bath_schedule(
        &mut self,
        parcel: &mut ClosedGasParcel,
        final_volume: Volume,
        schedule: HeatBathSchedule<'_>,
        controls: MaterialThermalCompressionSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MaterialThermalCompressionStep, Error> {
        if schedule.constant {
            return self.advance_with_heat_bath(
                parcel,
                final_volume,
                schedule.duration(),
                schedule.knots[0].bath,
                controls,
                cancelled,
            );
        }
        let (candidate, physical) = fixed(
            self,
            parcel,
            final_volume,
            schedule.duration(),
            schedule,
            controls,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        *parcel = candidate;
        Ok(physical)
    }

    /// Resolve one original finite body's programmed conductance and volume work.
    ///
    /// Source, mass, pressure, original h/phase and both work accounts survive
    /// every interval. A late refusal publishes neither endpoint, including a
    /// body whose internal candidate already crossed a phase boundary.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_finite_body_schedule<'m>(
        &mut self,
        parcel: &mut ClosedGasParcel,
        final_volume: Volume,
        body: FiniteHeatBody<'m>,
        schedule: ConductanceSchedule<'_>,
        controls: MaterialThermalCompressionSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MaterialThermalCompressionBodyStep<'m>, Error> {
        poll(&mut cancelled)?;
        schedule.bind(body)?;
        if schedule.constant {
            return self.advance_with_finite_body(
                parcel,
                final_volume,
                schedule.duration(),
                body,
                controls,
                cancelled,
            );
        }
        let (candidate, physical) = fixed(
            self,
            parcel,
            final_volume,
            schedule.duration(),
            BodyProgram { body, schedule },
            controls,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        *parcel = candidate;
        Ok(physical)
    }

    /// Adapt midpoint thermal forcing and volume/reaction splitting together.
    ///
    /// Each full/two-half probe uses its own physical midpoint bath. All five
    /// native comparison channels and global accepted/evaluation caps remain
    /// active. A varying fixed-volume bath cannot use the one-call shortcut.
    /// Local comparisons are not a global time-accuracy or entropy certificate.
    pub fn advance_with_heat_bath_schedule_adaptive(
        &mut self,
        parcel: &mut ClosedGasParcel,
        final_volume: Volume,
        schedule: HeatBathSchedule<'_>,
        controls: AdaptiveMaterialThermalCompressionSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AdaptiveMaterialThermalCompressionStep, AdaptiveMaterialThermalCompressionError>
    {
        if schedule.constant {
            return self.advance_with_heat_bath_adaptive(
                parcel,
                final_volume,
                schedule.duration(),
                schedule.knots[0].bath,
                controls,
                cancelled,
            );
        }
        let (candidate, report) = adaptive::advance(
            self,
            parcel,
            final_volume,
            schedule.duration(),
            schedule,
            controls,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        *parcel = candidate;
        Ok(report)
    }

    /// Adapt programmed conductance using actual original-body enthalpy error.
    ///
    /// A latent plateau's constant temperature never hides its energy/phase
    /// difference. Coarse and rejected gas/body states remain private, and their
    /// numerical work still consumes the complete stroke's explicit budget.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_finite_body_schedule_adaptive<'m>(
        &mut self,
        parcel: &mut ClosedGasParcel,
        final_volume: Volume,
        body: FiniteHeatBody<'m>,
        schedule: ConductanceSchedule<'_>,
        controls: AdaptiveMaterialThermalCompressionSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<
        AdaptiveMaterialThermalCompressionBodyStep<'m>,
        AdaptiveMaterialThermalCompressionError,
    > {
        poll(&mut cancelled)?;
        schedule.bind(body)?;
        if schedule.constant {
            return self.advance_with_finite_body_adaptive(
                parcel,
                final_volume,
                schedule.duration(),
                body,
                controls,
                cancelled,
            );
        }
        let (candidate, report) = adaptive::advance(
            self,
            parcel,
            final_volume,
            schedule.duration(),
            BodyProgram { body, schedule },
            controls,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        *parcel = candidate;
        Ok(report)
    }
}

impl Partner for HeatBathSchedule<'_> {
    type Physical = MaterialThermalCompressionStep;

    fn react<F: FnMut() -> bool>(
        self,
        reactor: &mut MaterialIsochoricReactor,
        parcel: &mut ClosedGasParcel,
        duration: Time,
        clock: IntervalTime,
        cancelled: &mut F,
    ) -> Result<(Self, ThermalInterval, [f64; 5]), Error> {
        let (_, physical, body) = self
            .at(clock.midpoint_s())?
            .react(reactor, parcel, duration, clock, cancelled)?;
        Ok((self, physical, body))
    }

    fn is_time_constant(self) -> bool {
        self.constant
    }

    fn next_boundary(self, after_s: f64, duration_s: f64) -> f64 {
        next_boundary(self.knots, |k| k.elapsed_time.value(), after_s, duration_s)
    }

    fn finish(self, gas: MaterialThermalCompressionStep, _: [f64; 5]) -> Self::Physical {
        gas
    }
}

#[derive(Clone, Copy)]
struct BodyProgram<'m, 's> {
    body: FiniteHeatBody<'m>,
    schedule: ConductanceSchedule<'s>,
}

impl<'m> Partner for BodyProgram<'m, '_> {
    type Physical = MaterialThermalCompressionBodyStep<'m>;

    fn react<F: FnMut() -> bool>(
        self,
        reactor: &mut MaterialIsochoricReactor,
        parcel: &mut ClosedGasParcel,
        duration: Time,
        clock: IntervalTime,
        cancelled: &mut F,
    ) -> Result<(Self, ThermalInterval, [f64; 5]), Error> {
        let body = self
            .body
            .with_conductance(self.schedule.at(clock.midpoint_s())?)
            .map_err(MaterialKineticsError::from)?;
        let (body, gas, changes) = body.react(reactor, parcel, duration, clock, cancelled)?;
        // Retain the exact prescribed endpoint coefficient without re-solving
        // the caloric state or adding a fictitious heat/work increment.
        let body = body
            .with_conductance(self.schedule.at(clock.end_s)?)
            .map_err(MaterialKineticsError::from)?;
        Ok((Self { body, ..self }, gas, changes))
    }

    fn is_time_constant(self) -> bool {
        self.schedule.constant
    }

    fn next_boundary(self, after_s: f64, duration_s: f64) -> f64 {
        next_boundary(
            self.schedule.knots,
            |k| k.elapsed_time.value(),
            after_s,
            duration_s,
        )
    }

    fn body_inventory(self) -> Option<[f64; 2]> {
        self.body.body_inventory()
    }

    fn finish(self, gas: MaterialThermalCompressionStep, body: [f64; 5]) -> Self::Physical {
        self.body.finish(gas, body)
    }
}

// Shared admission bounds schedule traversal; the immutable slices also prevent
// a stateful callback from changing the forcing between coarse/fine retries.
fn admit<T>(
    knots: &[T],
    time: impl Fn(&T) -> f64,
    valid: impl Fn(&T) -> bool,
    same: impl Fn(&T, &T) -> bool,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<bool, Error> {
    poll(cancelled)?;
    if knots.is_empty() || knots.len() > MAX_COMPRESSION_SUBSTEPS + 1 || time(&knots[0]) != 0.0 {
        return Err(Error::Invalid(
            "thermal program needs 1..=4097 knots starting at zero time",
        ));
    }
    let mut constant = true;
    let mut previous = 0.0;
    for (index, knot) in knots.iter().enumerate() {
        poll(cancelled)?;
        let current = time(knot);
        if !current.is_finite() || current < 0.0 || index > 0 && current <= previous || !valid(knot)
        {
            return Err(Error::Invalid(
                "thermal knots require advancing finite times and admitted physical values",
            ));
        }
        previous = current;
        constant &= same(knot, &knots[0]);
    }
    Ok(constant)
}

fn bracket<T>(knots: &[T], time: impl Fn(&T) -> f64, at: f64) -> Result<(usize, usize), Error> {
    let last = knots.len() - 1;
    if !at.is_finite() || at < 0.0 || at > time(&knots[last]) {
        return Err(Error::Invalid(
            "thermal program does not extrapolate physical time",
        ));
    }
    let upper = knots.partition_point(|k| time(k) <= at);
    let lower = upper - 1;
    Ok((
        lower,
        if time(&knots[lower]) == at {
            lower
        } else {
            upper
        },
    ))
}

// Evaluate from the nearer endpoint, preserving exact knot values and avoiding
// cancellation when an admitted coefficient falls through many orders of magnitude.
fn linear(start: f64, end: f64, at: f64, initial: f64, terminal: f64) -> f64 {
    let left = at - start;
    let right = end - at;
    if left <= right {
        (terminal - initial).mul_add(left / (end - start), initial)
    } else {
        (initial - terminal).mul_add(right / (end - start), terminal)
    }
}

fn next_boundary<T>(knots: &[T], time: impl Fn(&T) -> f64, after: f64, duration: f64) -> f64 {
    knots
        .get(knots.partition_point(|k| time(k) <= after))
        .map_or(duration, time)
}

// Uniform requested times and mandatory forcing knots form one bounded grid.
// The original path admission and interval/thermal totals own all physical work.
#[allow(clippy::too_many_arguments)]
fn fixed<B: Partner, F: FnMut() -> bool>(
    reactor: &mut MaterialThermalCompressionReactor,
    parcel: &ClosedGasParcel,
    final_volume: Volume,
    duration: Time,
    mut partner: B,
    controls: MaterialThermalCompressionSolve,
    cancelled: &mut F,
) -> Result<(ClosedGasParcel, B::Physical), Error> {
    poll(cancelled)?;
    let path = FixedCompressionPath::new(
        parcel.volume_m3(),
        final_volume,
        duration,
        controls.substeps,
    )?;
    // Admit the known union budget before any chemistry. A forcing knot already
    // on the represented uniform grid costs no extra interval. The bounded grid
    // index is rounded only to identify a candidate; exact time equality decides.
    let mut required = controls.substeps;
    let mut boundary = partner.next_boundary(0.0, path.duration_s);
    while boundary < path.duration_s {
        poll(cancelled)?;
        let index = ((boundary / path.duration_s) * controls.substeps as f64).round();
        let grid_time = path.duration_s * (index / controls.substeps as f64);
        if grid_time != boundary {
            required += 1;
            if required > MAX_COMPRESSION_SUBSTEPS {
                return Err(Error::Invalid(
                    "thermal program and uniform grid exceed 4096 intervals",
                ));
            }
        }
        boundary = partner.next_boundary(boundary, path.duration_s);
    }
    let mut candidate = parcel.clone();
    candidate.set_volume_isentropic(path.initial_volume_m3, controls.volume, &mut *cancelled)?;
    let volume_at = |time: f64| {
        if time == path.duration_s {
            path.final_volume_m3
        } else {
            path.initial_volume_m3 * fs_math::det::exp(path.log_ratio * (time / path.duration_s))
        }
    };
    let mut totals = Totals::default();
    let (mut elapsed, mut grid, mut intervals) = (0.0, 1, 0);
    while elapsed < path.duration_s {
        poll(cancelled)?;
        if intervals == MAX_COMPRESSION_SUBSTEPS {
            return Err(Error::Invalid(
                "thermal program and uniform grid exceed 4096 intervals",
            ));
        }
        let grid_end = if grid == controls.substeps {
            path.duration_s
        } else {
            path.duration_s * (grid as f64 / controls.substeps as f64)
        };
        let end = grid_end.min(partner.next_boundary(elapsed, path.duration_s));
        let clock = IntervalTime {
            start_s: elapsed,
            end_s: end,
        };
        let midpoint = clock.midpoint_s();
        if !(elapsed < midpoint && midpoint < end) {
            return Err(Error::Invalid(
                "thermal interval needs distinct represented start/midpoint/end times",
            ));
        }
        totals.interval(
            &mut reactor.reactor,
            &mut candidate,
            &mut partner,
            CompressionInterval {
                duration_s: end - elapsed,
                midpoint_volume_m3: volume_at(midpoint),
                final_volume_m3: volume_at(end),
            },
            clock,
            controls.volume,
            cancelled,
        )?;
        elapsed = end;
        intervals += 1;
        if end == grid_end {
            grid += 1;
        }
    }
    let physical = totals.finish(
        partner,
        path.initial_volume_m3,
        path.final_volume_m3,
        duration,
    )?;
    Ok((candidate, physical))
}
