//! Adaptive physical time for the original irreversible brittle motion.
//!
//! Each candidate compares one native backward-Euler step with two half steps
//! from the SAME accepted history. Only the two native half steps may publish:
//! Richardson extrapolation would manufacture displacement, velocity and damage
//! histories that have never passed the original material/energy admission.
//!
//! Local displacement, velocity, damage and fracture-work differences drive
//! refinement. The displacement/velocity relative scales use increments, not
//! accumulated position or uniform velocity. Physical damage jumps, absolute
//! constitutive quadrature defects and numerical kinetic loss have independent
//! bounds. Numerical losses remain separate from fracture and generated heat.
//!
//! This is a first-order, fixed-source/fixed-plane time path. Step doubling is
//! an error indicator, not a global error certificate or proof that an unseen
//! reversal/event was resolved. Split non-affine loads at their physical knots.

use fs_exec::Cx;
use fs_material::brittle::BrittleEnergy;
use fs_math::det;

use super::super::path::{Sum, interpolate};
use super::super::{BrittleTetError, checkpoint, energy_gate};
use super::{
    BrittleTetDynamicControls, BrittleTetDynamicEnergy, BrittleTetDynamicError, BrittleTetMotion,
};
use crate::plastic3::dynamics::InertialEnergy;
use crate::plastic3::{IncrementSettings, SmallStrainTetError};

/// Explicit physical intervals, native solve/storage budgets and local errors.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BrittleTetTimeControls<W = BrittleTetDynamicControls> {
    /// Duration of the first full-vs-two-halves comparison.
    pub initial_interval_s: f64,
    /// Smallest comparison interval, except a shorter final remainder.
    /// Accepted native timesteps are its two representable halves.
    pub minimum_interval_s: f64,
    pub maximum_interval_s: f64,
    /// Compact accepted comparison records; each uses TWO native motion steps.
    pub maximum_intervals: usize,
    /// All native solves, including discarded coarse steps and rejected trials.
    pub maximum_trials: usize,
    /// Positive absolute local comparison tolerances in the stated units.
    pub displacement_tolerance_m: f64,
    pub velocity_tolerance_mps: f64,
    pub damage_tolerance: f64,
    pub fracture_dissipation_tolerance_j: f64,
    /// Nonnegative relative tolerance on local motion/fracture increments.
    pub relative_tolerance: f64,
    /// Maximum damage increase per element in EITHER accepted half step.
    pub maximum_damage_increment: f64,
    /// Physical/kinetic identities apply locally and to the complete duration.
    /// With the default backward-Euler policy, constitutive-defect and
    /// numerical-kinetic-loss allowances are WHOLE duration budgets,
    /// partitioned by each actual native timestep.
    pub work: W,
}

/// Fine-vs-coarse differences before accepting the unmodified native fine state.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct BrittleTetTimeErrorEstimate {
    pub displacement_error_m: f64,
    pub velocity_error_mps: f64,
    pub maximum_damage_error: f64,
    /// Sum of absolute element differences, preventing spatial cancellation.
    pub fracture_dissipation_error_j: f64,
    pub normalized_error: f64,
}

/// Compact accepted pair; its two energy receipts belong to actual fine steps.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BrittleTetTimeRecord<E = BrittleTetDynamicEnergy> {
    pub start_time_s: f64,
    pub midpoint_time_s: f64,
    pub end_time_s: f64,
    pub estimate: BrittleTetTimeErrorEstimate,
    pub maximum_damage_increment: f64,
    pub energy: [E; 2],
}

/// Only a complete physical duration publishes its final original motion.
#[derive(Debug)]
pub struct BrittleTetTimeStep<'a> {
    pub motion: BrittleTetMotion<'a>,
    pub records: Vec<BrittleTetTimeRecord>,
    pub trials: usize,
    pub rejected_intervals: usize,
    pub energy: BrittleTetDynamicEnergy,
    pub absolute_constitutive_defect_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum BrittleTetTimeRefinementCause {
    Endpoint(BrittleTetDynamicError),
    ErrorEstimate(BrittleTetTimeErrorEstimate),
    DamageIncrement { observed: f64, maximum: f64 },
}

/// No refused duration publishes its internal partial progress.
#[derive(Clone, Debug, PartialEq)]
pub enum BrittleTetTimeError {
    Invalid(&'static str),
    Allocation,
    Endpoint(BrittleTetDynamicError),
    WorkLimit {
        accepted_intervals: usize,
        trials: usize,
    },
    MinimumInterval {
        time_s: f64,
        interval_s: f64,
        cause: BrittleTetTimeRefinementCause,
    },
}
impl core::fmt::Display for BrittleTetTimeError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adaptive brittle time path refused: {self:?}")
    }
}
impl std::error::Error for BrittleTetTimeError {}
impl From<BrittleTetDynamicError> for BrittleTetTimeError {
    fn from(error: BrittleTetDynamicError) -> Self {
        Self::Endpoint(error)
    }
}
impl From<BrittleTetError> for BrittleTetTimeError {
    fn from(error: BrittleTetError) -> Self {
        Self::Endpoint(error.into())
    }
}

impl<'a> BrittleTetMotion<'a> {
    /// Advance a complete affine physical load interval with native error control.
    ///
    /// Support coordinates interpolate from the original accepted displacement
    /// to their exact endpoint. Applied-force entries give explicit initial and
    /// final values; duplicates retain native summation semantics. Acceleration
    /// need not start in static equilibrium. Every candidate borrows the same
    /// original material owners and immutable cached reference mass.
    ///
    /// Replays consume the native-solve budget, but discarded coarse/rejected
    /// states consume no accepted motion steps. Each accepted pair consumes two.
    /// Only numerical Newton/linear refusals, the two finite numerical-work
    /// gates, excessive local error and damage jumps refine. Material/domain,
    /// geometry, physical-energy and cancellation refusals remain terminal.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_time(
        &self,
        duration_s: f64,
        prescribed_end_m: &[(usize, f64)],
        force_intervals_n: &[(usize, [f64; 2])],
        settings: IncrementSettings,
        controls: BrittleTetTimeControls,
        cx: &Cx<'_>,
    ) -> Result<BrittleTetTimeStep<'a>, BrittleTetTimeError> {
        let result = drive(
            self,
            duration_s,
            prescribed_end_m,
            force_intervals_n,
            settings,
            controls,
            BackwardEuler,
            cx,
        )?;
        let motion = result.motion;
        let [
            material,
            fracture,
            internal,
            dead,
            support,
            constitutive,
            absolute_constitutive,
            kinetic,
            numerical,
            inertial,
            viscous,
        ] = result.values;
        let storage_change = motion.solid.stored_energy_j - self.solid.stored_energy_j;
        let inertia = InertialEnergy {
            reference_mass_kg: self.reference_mass_kg(),
            kinetic_before_j: self.kinetic_energy_j,
            kinetic_after_j: motion.kinetic_energy_j,
            // Sum native average(v) dot M dv, avoiding large-K subtraction.
            kinetic_change_j: kinetic,
            numerical_loss_j: numerical,
            inertial_work_j: inertial,
            balance_error_j: (inertial - kinetic) - numerical,
        };
        let energy = BrittleTetDynamicEnergy {
            inertia,
            stored_before_j: self.solid.stored_energy_j,
            stored_after_j: motion.solid.stored_energy_j,
            material_work_j: material,
            fracture_dissipation_j: fracture,
            viscous_dissipation_j: viscous,
            material_balance_error_j: ((material - storage_change) - fracture) - viscous,
            internal_endpoint_work_j: internal,
            dead_load_endpoint_work_j: dead,
            support_endpoint_work_j: support,
            equilibrium_work_residual_j: ((dead + support) - internal) - inertial,
            constitutive_endpoint_defect_j: constitutive,
            total_balance_error_j: ((((((dead + support) - storage_change) - fracture) - viscous)
                - kinetic)
                - numerical)
                - constitutive,
        };
        for (account, defect, tolerance) in [
            (
                "time path material work/storage/fracture",
                energy.material_balance_error_j,
                controls.work.work.material_energy_tolerance_j,
            ),
            (
                "time path cumulative fracture loss",
                (motion.solid.fracture_dissipation_j - self.solid.fracture_dissipation_j)
                    - fracture,
                controls.work.work.material_energy_tolerance_j,
            ),
            (
                "time path applied/support/inertia equilibrium",
                energy.equilibrium_work_residual_j,
                controls.work.work.equilibrium_work_tolerance_j,
            ),
            (
                "time path constitutive quadrature",
                absolute_constitutive,
                controls.work.work.maximum_endpoint_work_defect_j,
            ),
            (
                "time path numerical kinetic loss",
                numerical,
                controls.work.maximum_numerical_kinetic_loss_j,
            ),
            (
                "time path native kinetic work identity",
                inertia.balance_error_j,
                controls.work.kinetic_work_tolerance_j,
            ),
            (
                "time path complete kinetic work balance",
                energy.total_balance_error_j,
                controls.work.work.material_energy_tolerance_j,
            ),
        ] {
            energy_gate(account, defect, tolerance)?;
        }
        checkpoint(cx)?;
        Ok(BrittleTetTimeStep {
            motion,
            records: result.records,
            trials: result.trials,
            rejected_intervals: result.rejected_intervals,
            energy,
            absolute_constitutive_defect_j: absolute_constitutive,
        })
    }
}

/// Shared native trial representation; it never fabricates a material state.
pub(super) struct TimeTrial<'a, E> {
    pub(super) motion: BrittleTetMotion<'a>,
    pub(super) element_energy: Vec<BrittleEnergy>,
    pub(super) energy: E,
    pub(super) numerical_tangent_evaluations: usize,
}

/// Private integration policy: physical sampling and energy interpretation stay
/// with the native scheme, while subdivision/history/error logic is shared.
pub(super) trait Scheme {
    type Work: Copy;
    type Energy: Copy;
    const MEAN_FORCES: bool;
    fn validate(&self, work: Self::Work) -> Result<(), BrittleTetTimeError>;
    #[allow(clippy::too_many_arguments)]
    fn advance<'a>(
        &self,
        old: &BrittleTetMotion<'a>,
        dt_s: f64,
        prescribed_m: &[(usize, f64)],
        forces_n: &[(usize, f64)],
        settings: IncrementSettings,
        work: Self::Work,
        part: f64,
        cx: &Cx<'_>,
    ) -> Result<TimeTrial<'a, Self::Energy>, BrittleTetDynamicError>;
    /// Material, fracture, internal, applied, support, signed work defect,
    /// absolute work defect, kinetic change, numerical loss, inertial work,
    /// and surviving-bond physical viscous loss.
    /// The concrete public result retains endpoint-vs-mean semantics.
    fn values(&self, energy: Self::Energy) -> [f64; 11];
    fn work_growth(&self, work: Self::Work, part: f64, energy: [Self::Energy; 2]) -> f64;
    fn retryable(&self, error: &BrittleTetDynamicError) -> bool;
}

pub(super) struct TimeTrajectory<'a, E> {
    pub(super) motion: BrittleTetMotion<'a>,
    pub(super) records: Vec<BrittleTetTimeRecord<E>>,
    pub(super) trials: usize,
    pub(super) rejected_intervals: usize,
    pub(super) values: [f64; 11],
    /// Tangent evaluations from completed native solves, including discarded
    /// coarse/rejected candidates.
    /// A failed native solve does not expose its partial tangent telemetry.
    pub(super) numerical_tangent_evaluations: usize,
}

struct BackwardEuler;
impl Scheme for BackwardEuler {
    type Work = BrittleTetDynamicControls;
    type Energy = BrittleTetDynamicEnergy;
    const MEAN_FORCES: bool = false;

    fn validate(&self, c: Self::Work) -> Result<(), BrittleTetTimeError> {
        if ![
            c.work.material_energy_tolerance_j,
            c.work.equilibrium_work_tolerance_j,
            c.work.maximum_endpoint_work_defect_j,
            c.kinetic_work_tolerance_j,
            c.maximum_numerical_kinetic_loss_j,
        ]
        .iter()
        .all(|v| v.is_finite() && *v >= 0.0)
        {
            return Err(BrittleTetTimeError::Invalid(
                "ordered positive intervals, positive error/solve budgets and finite nonnegative work gates required",
            ));
        }
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    fn advance<'a>(
        &self,
        old: &BrittleTetMotion<'a>,
        dt_s: f64,
        prescribed_m: &[(usize, f64)],
        forces_n: &[(usize, f64)],
        settings: IncrementSettings,
        controls: Self::Work,
        part: f64,
        cx: &Cx<'_>,
    ) -> Result<TimeTrial<'a, Self::Energy>, BrittleTetDynamicError> {
        let local = BrittleTetDynamicControls {
            work: super::super::BrittleTetControls {
                maximum_endpoint_work_defect_j: controls.work.maximum_endpoint_work_defect_j * part,
                ..controls.work
            },
            maximum_numerical_kinetic_loss_j: controls.maximum_numerical_kinetic_loss_j * part,
            ..controls
        };
        let step = old.advance(dt_s, prescribed_m, forces_n, settings, local, cx)?;
        Ok(TimeTrial {
            motion: step.motion,
            element_energy: step.element_energy,
            energy: step.energy,
            numerical_tangent_evaluations: 0,
        })
    }

    fn values(&self, e: Self::Energy) -> [f64; 11] {
        [
            e.material_work_j,
            e.fracture_dissipation_j,
            e.internal_endpoint_work_j,
            e.dead_load_endpoint_work_j,
            e.support_endpoint_work_j,
            e.constitutive_endpoint_defect_j,
            e.constitutive_endpoint_defect_j.abs(),
            e.inertia.kinetic_change_j,
            e.inertia.numerical_loss_j,
            e.inertia.inertial_work_j,
            e.viscous_dissipation_j,
        ]
    }

    fn work_growth(&self, c: Self::Work, part: f64, e: [Self::Energy; 2]) -> f64 {
        growth(
            c.work.maximum_endpoint_work_defect_j * part,
            e.iter()
                .map(|e| e.constitutive_endpoint_defect_j.abs())
                .sum(),
        )
        .min(growth(
            c.maximum_numerical_kinetic_loss_j * part,
            e.iter().map(|e| e.inertia.numerical_loss_j).sum(),
        ))
    }
    fn retryable(&self, error: &BrittleTetDynamicError) -> bool {
        newton_retryable(error)
            || matches!(error, BrittleTetDynamicError::Solid(BrittleTetError::EnergyGate {
                account: "dynamic constitutive endpoint integration defect"
                    | "backward-Euler numerical kinetic loss",
                defect_j, ..
            }) if defect_j.is_finite())
    }
}

#[allow(clippy::too_many_arguments)]
pub(super) fn drive<'a, S: Scheme>(
    original: &BrittleTetMotion<'a>,
    duration_s: f64,
    prescribed_end_m: &[(usize, f64)],
    force_intervals_n: &[(usize, [f64; 2])],
    settings: IncrementSettings,
    controls: BrittleTetTimeControls<S::Work>,
    scheme: S,
    cx: &Cx<'_>,
) -> Result<TimeTrajectory<'a, S::Energy>, BrittleTetTimeError> {
    checkpoint(cx)?;
    validate(&controls)?;
    scheme.validate(controls.work)?;
    let target_time = original.time_s + duration_s;
    if !duration_s.is_finite()
        || duration_s <= 0.0
        || !target_time.is_finite()
        || target_time <= original.time_s
    {
        return Err(BrittleTetTimeError::Invalid(
            "positive representable physical duration required",
        ));
    }
    let mut loads = Loads::new(
        original,
        target_time,
        prescribed_end_m,
        force_intervals_n,
        cx,
    )?;
    let duration = target_time - original.time_s;
    let mut current: Option<BrittleTetMotion<'a>> = None;
    let mut records = buffer(controls.maximum_intervals.min(64))?;
    let mut sums = [Sum::default(); 11];
    let mut interval = controls.initial_interval_s;
    let mut trials = 0;
    let mut rejected = 0;
    let mut tangent_evaluations = 0;
    loop {
        checkpoint(cx)?;
        let accepted = current.as_ref().unwrap_or(original);
        if accepted.time_s == target_time {
            break;
        }
        if records.len() >= controls.maximum_intervals || controls.maximum_trials - trials < 3 {
            return Err(BrittleTetTimeError::WorkLimit {
                accepted_intervals: records.len(),
                trials,
            });
        }
        if accepted.remaining_steps() < 2 {
            return Err(BrittleTetDynamicError::StepBudget {
                accepted_steps: accepted.accepted_steps,
                maximum_steps: accepted.maximum_steps,
            }
            .into());
        }
        let end = (accepted.time_s + interval).min(target_time);
        let span = end - accepted.time_s;
        let midpoint = accepted.time_s + 0.5 * span;
        if !span.is_finite() || span <= 0.0 || midpoint <= accepted.time_s || midpoint >= end {
            return Err(BrittleTetTimeError::Invalid(
                "two increasing native physical timesteps must be representable",
            ));
        }
        let attempted = compare_interval(
            accepted,
            midpoint,
            end,
            &mut loads,
            settings,
            controls,
            &scheme,
            &mut trials,
            &mut tangent_evaluations,
            cx,
        );
        let pair = match attempted {
            Ok(pair) => pair,
            Err(BrittleTetTimeError::Endpoint(error)) if scheme.retryable(&error) => {
                interval = refine(
                    accepted.time_s,
                    span,
                    controls.minimum_interval_s,
                    BrittleTetTimeRefinementCause::Endpoint(error),
                )?;
                rejected += 1;
                continue;
            }
            Err(error) => return Err(error),
        };
        let cause = if pair.estimate.normalized_error > 1.0 {
            Some(BrittleTetTimeRefinementCause::ErrorEstimate(pair.estimate))
        } else if pair.maximum_damage_increment > controls.maximum_damage_increment {
            Some(BrittleTetTimeRefinementCause::DamageIncrement {
                observed: pair.maximum_damage_increment,
                maximum: controls.maximum_damage_increment,
            })
        } else {
            None
        };
        if let Some(cause) = cause {
            interval = refine(accepted.time_s, span, controls.minimum_interval_s, cause)?;
            rejected += 1;
            continue;
        }
        let energies = [pair.first_energy, pair.fine.energy];
        for energy in energies {
            for (sum, value) in sums.iter_mut().zip(scheme.values(energy)) {
                sum.add(value);
                if !sum.value().is_finite() {
                    return Err(BrittleTetTimeError::Invalid(
                        "complete physical-time work must be representable",
                    ));
                }
            }
        }
        if records.len() == records.capacity() {
            records
                .try_reserve(1)
                .map_err(|_| BrittleTetTimeError::Allocation)?;
        }
        records.push(BrittleTetTimeRecord {
            start_time_s: accepted.time_s,
            midpoint_time_s: midpoint,
            end_time_s: end,
            estimate: pair.estimate,
            maximum_damage_increment: pair.maximum_damage_increment,
            energy: energies,
        });
        // The same bounded conservative controller is retained at nonsmooth
        // fracture events; it does not presume a smooth higher-order solution.
        let error_growth = if pair.estimate.normalized_error > 0.0 {
            (0.9 / det::sqrt(pair.estimate.normalized_error)).clamp(0.5, 2.0)
        } else {
            2.0
        };
        let damage_growth = growth(
            controls.maximum_damage_increment,
            pair.maximum_damage_increment,
        );
        let work_growth = scheme.work_growth(controls.work, span / duration, energies);
        interval = (span * error_growth.min(damage_growth).min(work_growth))
            .clamp(controls.minimum_interval_s, controls.maximum_interval_s);
        current = Some(pair.fine.motion);
    }
    let motion = current.ok_or(BrittleTetTimeError::Invalid(
        "no physical-time interval accepted",
    ))?;
    checkpoint(cx)?;
    Ok(TimeTrajectory {
        motion,
        records,
        trials,
        rejected_intervals: rejected,
        values: sums.map(Sum::value),
        numerical_tangent_evaluations: tangent_evaluations,
    })
}

struct Compared<'a, E> {
    fine: TimeTrial<'a, E>,
    first_energy: E,
    estimate: BrittleTetTimeErrorEstimate,
    maximum_damage_increment: f64,
}

#[allow(clippy::too_many_arguments)]
fn compare_interval<'a, S: Scheme>(
    old: &BrittleTetMotion<'a>,
    midpoint: f64,
    end: f64,
    loads: &mut Loads<'_>,
    settings: IncrementSettings,
    controls: BrittleTetTimeControls<S::Work>,
    scheme: &S,
    trials: &mut usize,
    tangent_evaluations: &mut usize,
    cx: &Cx<'_>,
) -> Result<Compared<'a, S::Energy>, BrittleTetTimeError> {
    *trials += 1;
    let coarse = loads.advance(
        old,
        end,
        settings,
        controls.work,
        scheme,
        tangent_evaluations,
        cx,
    )?;
    *trials += 1;
    let first = loads.advance(
        old,
        midpoint,
        settings,
        controls.work,
        scheme,
        tangent_evaluations,
        cx,
    )?;
    *trials += 1;
    let fine = loads.advance(
        &first.motion,
        end,
        settings,
        controls.work,
        scheme,
        tangent_evaluations,
        cx,
    )?;
    let (estimate, maximum_damage_increment) =
        estimate(old, &coarse, &first, &fine, controls, scheme, cx)?;
    Ok(Compared {
        fine,
        first_energy: first.energy,
        estimate,
        maximum_damage_increment,
    })
}

fn estimate<S: Scheme>(
    old: &BrittleTetMotion<'_>,
    coarse: &TimeTrial<'_, S::Energy>,
    first: &TimeTrial<'_, S::Energy>,
    fine: &TimeTrial<'_, S::Energy>,
    c: BrittleTetTimeControls<S::Work>,
    scheme: &S,
    cx: &Cx<'_>,
) -> Result<(BrittleTetTimeErrorEstimate, f64), BrittleTetTimeError> {
    let mut out = BrittleTetTimeErrorEstimate::default();
    for node in 0..old.velocity_mps.len() {
        if node % 256 == 0 {
            checkpoint(cx)?;
        }
        for (kind, a, b, start, absolute) in [
            (
                0,
                coarse.motion.solid.solution.displacement_m[node],
                fine.motion.solid.solution.displacement_m[node],
                old.solid.solution.displacement_m[node],
                c.displacement_tolerance_m,
            ),
            (
                1,
                coarse.motion.velocity_mps[node],
                fine.motion.velocity_mps[node],
                old.velocity_mps[node],
                c.velocity_tolerance_mps,
            ),
        ] {
            let difference = norm(core::array::from_fn(|i| b[i] - a[i]));
            let scale = norm(core::array::from_fn(|i| b[i] - start[i]))
                .max(norm(core::array::from_fn(|i| a[i] - start[i])));
            let ratio = normalized(difference, absolute, c.relative_tolerance, scale)?;
            out.normalized_error = out.normalized_error.max(ratio);
            if kind == 0 {
                out.displacement_error_m = out.displacement_error_m.max(difference);
            } else {
                out.velocity_error_mps = out.velocity_error_mps.max(difference);
            }
        }
    }
    let mut fracture_error = Sum::default();
    let mut maximum_increment = 0.0_f64;
    for element in 0..old.solid.materials.len() {
        if element % 256 == 0 {
            checkpoint(cx)?;
        }
        let before = old.solid.solution.element_states[element].damage();
        let middle = first.motion.solid.solution.element_states[element].damage();
        let after = fine.motion.solid.solution.element_states[element].damage();
        let broad = coarse.motion.solid.solution.element_states[element].damage();
        let difference = (after - broad).abs();
        out.maximum_damage_error = out.maximum_damage_error.max(difference);
        maximum_increment = maximum_increment.max(middle - before).max(after - middle);
        out.normalized_error = out.normalized_error.max(normalized(
            difference,
            c.damage_tolerance,
            c.relative_tolerance,
            (after - before).abs().max((broad - before).abs()),
        )?);
        // Compare original incremental work, not differences of large histories.
        let narrow =
            first.element_energy[element].dissipated + fine.element_energy[element].dissipated;
        let wide = coarse.element_energy[element].dissipated;
        fracture_error.add(old.solid.bands[element].reference_volume_m3() * (narrow - wide).abs());
    }
    out.fracture_dissipation_error_j = fracture_error.value();
    out.normalized_error = out.normalized_error.max(normalized(
        out.fracture_dissipation_error_j,
        c.fracture_dissipation_tolerance_j,
        c.relative_tolerance,
        (scheme.values(first.energy)[1] + scheme.values(fine.energy)[1])
            .max(scheme.values(coarse.energy)[1]),
    )?);
    Ok((out, maximum_increment))
}

fn normalized(
    error: f64,
    absolute: f64,
    relative: f64,
    scale: f64,
) -> Result<f64, BrittleTetTimeError> {
    let denominator = relative.mul_add(scale, absolute);
    if !error.is_finite()
        || !scale.is_finite()
        || !denominator.is_finite()
        || error < 0.0
        || denominator <= 0.0
    {
        return Err(BrittleTetTimeError::Invalid(
            "finite local error normalization required",
        ));
    }
    Ok(error / denominator)
}
fn norm(v: [f64; 3]) -> f64 {
    det::hypot(det::hypot(v[0], v[1]), v[2])
}

struct Loads<'a> {
    start_time_s: f64,
    end_time_s: f64,
    initial_u: &'a [f64],
    prescribed_end: &'a [(usize, f64)],
    force_intervals: &'a [(usize, [f64; 2])],
    prescribed: Vec<(usize, f64)>,
    forces: Vec<(usize, f64)>,
}
impl<'a> Loads<'a> {
    fn new(
        old: &'a BrittleTetMotion<'_>,
        end_time_s: f64,
        prescribed_end: &'a [(usize, f64)],
        force_intervals: &'a [(usize, [f64; 2])],
        cx: &Cx<'_>,
    ) -> Result<Self, BrittleTetTimeError> {
        let initial_u = old.solid.solution.displacement_m.as_flattened();
        if prescribed_end.len() > initial_u.len() || force_intervals.len() > initial_u.len() {
            return Err(BrittleTetTimeError::Invalid(
                "load entries exceed original DOF count",
            ));
        }
        let mut prescribed = buffer(prescribed_end.len())?;
        for &(dof, target) in prescribed_end {
            checkpoint(cx)?;
            let Some(&start) = initial_u.get(dof) else {
                return Err(BrittleTetTimeError::Invalid("prescribed DOF outside mesh"));
            };
            if !target.is_finite() {
                return Err(BrittleTetTimeError::Invalid(
                    "finite prescribed target required",
                ));
            }
            prescribed.push((dof, start));
        }
        let mut forces = buffer(force_intervals.len())?;
        for &(dof, values) in force_intervals {
            checkpoint(cx)?;
            if dof >= initial_u.len() || values.iter().any(|x| !x.is_finite()) {
                return Err(BrittleTetTimeError::Invalid(
                    "finite applied forces on existing DOFs required",
                ));
            }
            forces.push((dof, values[0]));
        }
        Ok(Self {
            start_time_s: old.time_s,
            end_time_s,
            initial_u,
            prescribed_end,
            force_intervals,
            prescribed,
            forces,
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn advance<'m, S: Scheme>(
        &mut self,
        old: &BrittleTetMotion<'m>,
        target: f64,
        settings: IncrementSettings,
        controls: S::Work,
        scheme: &S,
        tangent_evaluations: &mut usize,
        cx: &Cx<'_>,
    ) -> Result<TimeTrial<'m, S::Energy>, BrittleTetTimeError> {
        let duration = self.end_time_s - self.start_time_s;
        let fraction = if target == self.end_time_s {
            1.0
        } else {
            (target - self.start_time_s) / duration
        };
        // Average interval endpoint fractions, avoiding an absolute midpoint
        // clock that can round onto one endpoint at a large accumulated time.
        let force_fraction = if S::MEAN_FORCES {
            0.5 * ((old.time_s - self.start_time_s) / duration) + 0.5 * fraction
        } else {
            fraction
        };
        for (i, (entry, &(dof, end))) in self
            .prescribed
            .iter_mut()
            .zip(self.prescribed_end)
            .enumerate()
        {
            if i % 256 == 0 {
                checkpoint(cx)?;
            }
            entry.1 = interpolate(self.initial_u[dof], end, fraction);
        }
        for (i, (entry, &(_, [start, end]))) in
            self.forces.iter_mut().zip(self.force_intervals).enumerate()
        {
            if i % 256 == 0 {
                checkpoint(cx)?;
            }
            entry.1 = interpolate(start, end, force_fraction);
        }
        let dt = target - old.time_s;
        let step = scheme.advance(
            old,
            dt,
            &self.prescribed,
            &self.forces,
            settings,
            controls,
            dt / duration,
            cx,
        )?;
        *tangent_evaluations = (*tangent_evaluations)
            .checked_add(step.numerical_tangent_evaluations)
            .ok_or(BrittleTetTimeError::Invalid(
                "completed tangent count overflow",
            ))?;
        if step.motion.time_s != target {
            return Err(BrittleTetTimeError::Invalid(
                "native clock did not reach the requested physical time",
            ));
        }
        Ok(step)
    }
}

fn validate<W>(c: &BrittleTetTimeControls<W>) -> Result<(), BrittleTetTimeError> {
    if ![
        c.initial_interval_s,
        c.minimum_interval_s,
        c.maximum_interval_s,
        c.displacement_tolerance_m,
        c.velocity_tolerance_mps,
        c.damage_tolerance,
        c.fracture_dissipation_tolerance_j,
        c.maximum_damage_increment,
    ]
    .iter()
    .all(|x| x.is_finite() && *x > 0.0)
        || c.minimum_interval_s > c.initial_interval_s
        || c.initial_interval_s > c.maximum_interval_s
        || c.maximum_damage_increment > 1.0
        || c.maximum_intervals == 0
        || c.maximum_trials < 3
        || !c.relative_tolerance.is_finite()
        || c.relative_tolerance < 0.0
    {
        return Err(BrittleTetTimeError::Invalid(
            "ordered positive intervals, positive error/solve budgets and finite nonnegative work gates required",
        ));
    }
    Ok(())
}
fn buffer<T>(capacity: usize) -> Result<Vec<T>, BrittleTetTimeError> {
    let mut result = Vec::new();
    result
        .try_reserve_exact(capacity)
        .map_err(|_| BrittleTetTimeError::Allocation)?;
    Ok(result)
}
fn growth(allowance: f64, used: f64) -> f64 {
    if used > 0.0 {
        (0.8 * allowance / used).clamp(0.5, 2.0)
    } else {
        2.0
    }
}
pub(super) fn newton_retryable(error: &BrittleTetDynamicError) -> bool {
    matches!(
        error,
        BrittleTetDynamicError::Solid(BrittleTetError::Solve(
            SmallStrainTetError::NewtonStalled { .. }
                | SmallStrainTetError::LinearSolveFailed { .. }
        ))
    )
}
fn refine(
    time_s: f64,
    interval_s: f64,
    minimum: f64,
    cause: BrittleTetTimeRefinementCause,
) -> Result<f64, BrittleTetTimeError> {
    let smaller = 0.5 * interval_s;
    if smaller < minimum || time_s + smaller <= time_s {
        Err(BrittleTetTimeError::MinimumInterval {
            time_s,
            interval_s,
            cause,
        })
    } else {
        Ok(smaller)
    }
}
