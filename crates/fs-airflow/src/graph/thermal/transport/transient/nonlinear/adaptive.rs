//! Resolve prescribed-flow convection in time using the original endpoint solve.
//!
//! A complete native backward-Euler step is compared with two native half
//! steps. Every fluid and wall cell participates, including enthalpy differences
//! on an isothermal latent plateau. Only accepted half steps carry heat credit.
//! This is a local time-discretization estimate, not a rigorous error enclosure.

use std::cell::Cell;

use fs_exec::Cx;
use fs_math::eft::two_sum;
use fs_thermal::enthalpy::backward_euler::BoundaryStepCoefficients;
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation, EnergyAudit};

use super::{
    GraphEnthalpyTransport, ImplicitTransportControls, ImplicitTransportError,
    ImplicitTransportReport, nonnegative, positive, thermal_checkpoint,
};

/// Physical resolution and total work for one immutable candidate window.
#[derive(Clone, Copy, Debug)]
pub struct ImplicitTransportWindowControls {
    /// Full-versus-two-half comparison interval; accepted steps are its halves.
    pub maximum_interval_s: f64,
    /// Refinement floor. A final or schedule-clipped interval may be shorter,
    /// but an unresolved short interval refuses instead of bypassing tolerance.
    pub minimum_interval_s: f64,
    pub maximum_native_steps: usize,
    /// Also bounded by the original thermal owner's max_substeps.
    pub maximum_accepted_steps: usize,
    /// All source endpoint calls, including derivative and discarded trials.
    /// Each native solve retains its own endpoint-evaluation ceiling as well.
    pub maximum_source_evaluations: usize,
    /// Per-cell specific storage-coordinate difference [J/kg], including all
    /// finite wall cells. In a latent plateau this still measures phase energy.
    pub absolute_specific_energy_error_j_kg: f64,
    /// Per-cell original T(h) difference [K]. No heat capacity is substituted.
    pub absolute_temperature_error_k: f64,
    /// Scales actual interval changes in h or T, never absolute energy datums.
    pub relative_state_error: f64,
    /// Bounds the independently measured whole-window balance and the sum of
    /// absolute original accepted-step balance defects separately [J].
    pub energy_tolerance_j: f64,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct ImplicitTransportWindowWork {
    pub native_steps: usize,
    pub source_evaluations: usize,
    pub rejected_intervals: usize,
}

#[derive(Clone, Copy, Debug)]
pub struct ImplicitTransportWindowEstimate {
    pub start_s: f64,
    pub end_s: f64,
    pub maximum_specific_energy_difference_j_kg: f64,
    pub maximum_temperature_difference_k: f64,
    pub error_ratio: f64,
}

/// One accepted original step and its actual physical interval/source payload.
#[derive(Debug)]
pub struct ImplicitTransportWindowStep<T> {
    pub start_s: f64,
    pub end_s: f64,
    pub report: ImplicitTransportReport<T>,
}

#[derive(Debug)]
pub struct ImplicitTransportWindowReport<T> {
    pub steps: Vec<ImplicitTransportWindowStep<T>>,
    /// Accepted comparison intervals only; no rejected trajectory is retained.
    pub estimates: Vec<ImplicitTransportWindowEstimate>,
    pub work: ImplicitTransportWindowWork,
    /// Endpoint storage is measured from original per-cell energy changes.
    /// External heat sums ONLY accepted original native audits.
    pub thermal: EnergyAudit,
    pub sum_absolute_native_energy_defect_j: f64,
}

#[derive(Debug)]
pub enum ImplicitTransportWindowError<E> {
    Invalid(&'static str),
    Thermal(DepositionError),
    Budget {
        resource: &'static str,
        work: ImplicitTransportWindowWork,
    },
    Step {
        start_s: f64,
        end_s: f64,
        work: ImplicitTransportWindowWork,
        source: ImplicitTransportError<E>,
    },
    Underresolved {
        estimate: ImplicitTransportWindowEstimate,
        work: ImplicitTransportWindowWork,
    },
    Energy {
        residual_j: f64,
        sum_absolute_native_defect_j: f64,
        tolerance_j: f64,
    },
}
impl<E> From<DepositionError> for ImplicitTransportWindowError<E> {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error)
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for ImplicitTransportWindowError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adaptive graph convection: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for ImplicitTransportWindowError<E> {}

impl GraphEnthalpyTransport {
    /// Resolve one prescribed-flow nonlinear convection window atomically.
    ///
    /// This repeats `advance_implicit_with_heat`, with its original upwind or
    /// limited spatial transport, caloric laws, source admission and cell-energy
    /// acceptance. No Richardson-extrapolated state is ever published. Native
    /// scheduled pump/inlet slope changes are exact comparison boundaries;
    /// every coarse and half step samples its own original endpoint schedule.
    ///
    /// `resolve` receives the ORIGINAL START of the particular native interval
    /// and its ACTUAL trial endpoint. Duration-dependent conducting films must
    /// use that pair, including rejected comparisons. The callback is pure;
    /// its payload/heat is retained only after native and time-error admission.
    /// All source, numerical and physical refusals are terminal, including in
    /// discarded coarse solves; only a completed time-error comparison retries.
    ///
    /// The generic caloric owner exposes h and T(h), not a phase fraction. Both
    /// are checked independently for every cell, with relative scales based on
    /// changes since the comparison start. A large formation-energy datum, a
    /// stationary latent temperature or cancellation of opposing body changes
    /// cannot hide the transported energy difference.
    ///
    /// Structural/heat-impulse timeline events remain explicit native handoffs:
    /// split this call at their left limit and process them through advance_to.
    /// Events exactly at target remain pending. An interior event refuses
    /// before source work. Fixed-mass geometry, original laws and source limits
    /// stay authoritative. Refine spatial cells and maximum interval as well;
    /// endpoint sampling cannot resolve arbitrary unsampled forcing frequencies.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_implicit_adaptive_with_heat<'m, T, E>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        step: ImplicitTransportControls,
        controls: ImplicitTransportWindowControls,
        resolve: impl Fn(
            &Cx<'_>,
            &DepositionSimulation<'m>,
            &DepositionSimulation<'m>,
        ) -> Result<BoundaryStepCoefficients<T>, E>,
    ) -> Result<
        (DepositionSimulation<'m>, ImplicitTransportWindowReport<T>),
        ImplicitTransportWindowError<E>,
    > {
        thermal_checkpoint(cx)?;
        if !target_s.is_finite()
            || target_s <= thermal.time_s()
            || !(target_s - thermal.time_s()).is_finite()
            || !positive(controls.maximum_interval_s)
            || !positive(controls.minimum_interval_s)
            || controls.minimum_interval_s > controls.maximum_interval_s
            || !(1..=1_000_000).contains(&controls.maximum_native_steps)
            || !(2..=4096).contains(&controls.maximum_accepted_steps)
            || !(1..=1_000_000).contains(&controls.maximum_source_evaluations)
            || !positive(controls.absolute_specific_energy_error_j_kg)
            || !positive(controls.absolute_temperature_error_k)
            || !nonnegative(controls.relative_state_error)
            || controls.relative_state_error >= 1.0
            || !nonnegative(controls.energy_tolerance_j)
            || thermal.cells().is_empty()
            || thermal.cells().len() > step.maximum_cells
            || step.maximum_cells > super::MAX_DISTRIBUTED_FLUID_CELLS
            || thermal
                .next_event_time_s()
                .is_some_and(|event| event < target_s)
        {
            return Err(ImplicitTransportWindowError::Invalid(
                "bounded forward window, explicit state tolerances and no interior native thermal events required",
            ));
        }
        // The original endpoint driver admits all solver and native topology
        // controls before its first physical proposal; these are not replaced.
        let maximum_accepted = controls
            .maximum_accepted_steps
            .min(thermal.budget().max_substeps);
        let mut current = thermal.clone();
        let mut work = ImplicitTransportWindowWork::default();
        let mut steps = Vec::new();
        let mut estimates = Vec::new();
        let mut external = Sum::default();
        let mut absolute_defect = Sum::default();
        let mut proposed = controls.maximum_interval_s;
        while current.time_s() < target_s {
            thermal_checkpoint(cx)?;
            if steps.len() + 2 > maximum_accepted {
                return Err(ImplicitTransportWindowError::Budget {
                    resource: "accepted native thermal steps",
                    work,
                });
            }
            let start = current.time_s();
            let stop = self
                .schedule()
                .and_then(|s| s.next_transition_time_s(start))
                .map_or(target_s, |knot| knot.min(target_s));
            let end = if proposed >= stop - start {
                stop
            } else {
                start + proposed
            };
            let width = end - start;
            let middle = start + 0.5 * width;
            if !end.is_finite() || !(start < middle && middle < end) {
                return Err(ImplicitTransportWindowError::Invalid(
                    "native comparison interval is not representable",
                ));
            }
            let (coarse, _) =
                native_step(self, cx, &current, end, step, controls, &mut work, &resolve)?;
            let (half, first) = native_step(
                self, cx, &current, middle, step, controls, &mut work, &resolve,
            )?;
            let (fine, second) =
                native_step(self, cx, &half, end, step, controls, &mut work, &resolve)?;
            let estimate = compare(cx, &current, &coarse, &fine, controls)?;
            if estimate.error_ratio > 1.0 {
                work.rejected_intervals += 1;
                proposed = 0.5 * width;
                if proposed < controls.minimum_interval_s || start + proposed == start {
                    return Err(ImplicitTransportWindowError::Underresolved { estimate, work });
                }
                continue;
            }
            thermal_checkpoint(cx)?;
            for (a, b, report) in [(start, middle, first), (middle, end, second)] {
                if report.thermal.substeps != 1
                    || report.thermal.events != 0
                    || report.thermal.incoming_mass_kg != 0.0
                    || report.thermal.incoming_energy_j != 0.0
                {
                    return Err(ImplicitTransportWindowError::Invalid(
                        "native implicit interval changed fixed physical inventory",
                    ));
                }
                external.add(report.thermal.external_heat_j)?;
                absolute_defect.add(report.thermal.residual_j.abs())?;
                steps.push(ImplicitTransportWindowStep {
                    start_s: a,
                    end_s: b,
                    report,
                });
            }
            proposed = if estimate.error_ratio <= 0.25 {
                (2.0 * width).min(controls.maximum_interval_s)
            } else {
                width
            };
            current = fine;
            estimates.push(estimate);
        }
        let mut stored = Sum::default();
        for (id, old) in thermal.cells() {
            thermal_checkpoint(cx)?;
            let next = &current.cells()[id];
            stored.add(finite(
                old.mass_kg * (next.specific_h_j_kg - old.specific_h_j_kg),
            )?)?;
        }
        let residual = stored.difference(external)?;
        let defect = absolute_defect.value()?;
        if residual.abs() > controls.energy_tolerance_j || defect > controls.energy_tolerance_j {
            return Err(ImplicitTransportWindowError::Energy {
                residual_j: residual,
                sum_absolute_native_defect_j: defect,
                tolerance_j: controls.energy_tolerance_j,
            });
        }
        thermal_checkpoint(cx)?;
        Ok((
            current,
            ImplicitTransportWindowReport {
                thermal: EnergyAudit {
                    incoming_mass_kg: 0.0,
                    incoming_energy_j: 0.0,
                    external_heat_j: external.value()?,
                    stored_energy_change_j: stored.value()?,
                    residual_j: residual,
                    substeps: steps.len(),
                    events: 0,
                },
                steps,
                estimates,
                work,
                sum_absolute_native_energy_defect_j: defect,
            },
        ))
    }
}

#[allow(clippy::too_many_arguments)]
fn native_step<'m, T, E>(
    transport: &GraphEnthalpyTransport,
    cx: &Cx<'_>,
    initial: &DepositionSimulation<'m>,
    target_s: f64,
    mut step: ImplicitTransportControls,
    controls: ImplicitTransportWindowControls,
    work: &mut ImplicitTransportWindowWork,
    resolve: &impl Fn(
        &Cx<'_>,
        &DepositionSimulation<'m>,
        &DepositionSimulation<'m>,
    ) -> Result<BoundaryStepCoefficients<T>, E>,
) -> Result<(DepositionSimulation<'m>, ImplicitTransportReport<T>), ImplicitTransportWindowError<E>>
{
    thermal_checkpoint(cx)?;
    if work.native_steps == controls.maximum_native_steps {
        return Err(ImplicitTransportWindowError::Budget {
            resource: "native step attempts",
            work: *work,
        });
    }
    let remaining = controls.maximum_source_evaluations - work.source_evaluations;
    if remaining == 0 {
        return Err(ImplicitTransportWindowError::Budget {
            resource: "source evaluations",
            work: *work,
        });
    }
    step.maximum_endpoint_evaluations = step.maximum_endpoint_evaluations.min(remaining);
    work.native_steps += 1;
    let calls = Cell::new(0usize);
    let result = transport.advance_implicit_with_heat(
        cx,
        initial,
        target_s,
        step,
        &vec![0.0; initial.cells().len()],
        |cx, candidate| {
            calls.set(calls.get() + 1);
            resolve(cx, initial, candidate)
        },
    );
    // A failed caloric or transport admission may stop before the callback;
    // this count is source work, while each native solve separately caps all
    // endpoint evaluations. Every refusal is terminal, including budget refusal.
    work.source_evaluations += calls.get();
    result.map_err(|source| ImplicitTransportWindowError::Step {
        start_s: initial.time_s(),
        end_s: target_s,
        work: *work,
        source,
    })
}

fn compare<E>(
    cx: &Cx<'_>,
    initial: &DepositionSimulation<'_>,
    coarse: &DepositionSimulation<'_>,
    fine: &DepositionSimulation<'_>,
    controls: ImplicitTransportWindowControls,
) -> Result<ImplicitTransportWindowEstimate, ImplicitTransportWindowError<E>> {
    if coarse.time_s() != fine.time_s()
        || coarse.cells().len() != initial.cells().len()
        || fine.cells().len() != initial.cells().len()
    {
        return Err(ImplicitTransportWindowError::Invalid(
            "native comparison inventories or times differ",
        ));
    }
    let mut estimate = ImplicitTransportWindowEstimate {
        start_s: initial.time_s(),
        end_s: fine.time_s(),
        maximum_specific_energy_difference_j_kg: 0.0,
        maximum_temperature_difference_k: 0.0,
        error_ratio: 0.0,
    };
    for ((id, old), ((a_id, a), (b_id, b))) in initial
        .cells()
        .iter()
        .zip(coarse.cells().iter().zip(fine.cells()))
    {
        thermal_checkpoint(cx)?;
        if id != a_id
            || id != b_id
            || old.mass_kg != a.mass_kg
            || old.mass_kg != b.mass_kg
            || old.material_id != a.material_id
            || old.material_id != b.material_id
            || old.geometry != a.geometry
            || old.geometry != b.geometry
        {
            return Err(ImplicitTransportWindowError::Invalid(
                "native comparison changed physical storage",
            ));
        }
        let h = finite(b.specific_h_j_kg - a.specific_h_j_kg)?.abs();
        let t0 = initial.temperature_k(*id)?;
        let ta = coarse.temperature_k(*id)?;
        let tb = fine.temperature_k(*id)?;
        let t = finite(tb - ta)?.abs();
        let hs = finite(
            controls.absolute_specific_energy_error_j_kg
                + controls.relative_state_error
                    * (a.specific_h_j_kg - old.specific_h_j_kg)
                        .abs()
                        .max((b.specific_h_j_kg - old.specific_h_j_kg).abs()),
        )?;
        let ts = finite(
            controls.absolute_temperature_error_k
                + controls.relative_state_error * (ta - t0).abs().max((tb - t0).abs()),
        )?;
        estimate.maximum_specific_energy_difference_j_kg =
            estimate.maximum_specific_energy_difference_j_kg.max(h);
        estimate.maximum_temperature_difference_k =
            estimate.maximum_temperature_difference_k.max(t);
        estimate.error_ratio = estimate
            .error_ratio
            .max(finite(h / hs)?)
            .max(finite(t / ts)?);
    }
    Ok(estimate)
}

fn finite<E>(value: f64) -> Result<f64, ImplicitTransportWindowError<E>> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(ImplicitTransportWindowError::Invalid(
            "unrepresentable physical window arithmetic",
        ))
    }
}

#[derive(Clone, Copy, Default)]
struct Sum {
    sum: f64,
    correction: f64,
}
impl Sum {
    fn add<E>(&mut self, value: f64) -> Result<(), ImplicitTransportWindowError<E>> {
        let (sum, error) = two_sum(self.sum, finite(value)?);
        self.sum = finite(sum)?;
        self.correction = finite(self.correction + error)?;
        self.value().map(|_| ())
    }
    fn value<E>(self) -> Result<f64, ImplicitTransportWindowError<E>> {
        finite(self.sum + self.correction)
    }
    fn difference<E>(self, other: Self) -> Result<f64, ImplicitTransportWindowError<E>> {
        let (difference, roundoff) = two_sum(self.sum, -other.sum);
        finite(difference + (roundoff + self.correction - other.correction))
    }
}
