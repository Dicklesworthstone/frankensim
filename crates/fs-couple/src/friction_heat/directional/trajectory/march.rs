//! Time refinement of the original sourced rigid-friction/enthalpy phase.
//!
//! A full phase and two half phases start from the SAME accepted checkpoint.
//! Their differences control body momenta, individual cell energies, curved
//! travel and physical work. Only the two-half-phase result can advance. This
//! estimates the local error of the endpoint coefficient/thermal discretization;
//! it is not an interval enclosure or a global error bound. No extrapolated
//! body, heat or source state is manufactured from the comparison.
//!
//! Poses, geometry and normal load remain fixed over this operation. Constitutive
//! speed is phase-average unsigned speed, approaching instantaneous speed under
//! time refinement. Temperature, path and original source receipts refresh at
//! every accepted half phase. There is no imposed tangential drive or event
//! scheduler; pending thermal events and unsupported source domains still refuse.

use super::*;
use fs_contact::interface_binding::thermal_friction::FrictionQueryAxisKind;

/// Absolute physical error scales, with one relative increment scale.
/// Relative momentum scales use CHANGES from the step origin, not a body's
/// background momentum; thermal scales likewise use cell energy changes.
#[derive(Clone, Copy, Debug)]
pub struct DirectionalMarchTolerance {
    pub linear_momentum_ns: f64,
    pub angular_momentum_nms: f64,
    pub cell_energy_j: f64,
    pub distance_m: f64,
    pub dissipated_j: f64,
    pub relative: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct DirectionalMarchControls {
    pub initial_step_s: f64,
    /// Smallest attempted full comparison interval. The final remainder may
    /// be smaller; each accepted interval consists of two half phases.
    pub minimum_step_s: f64,
    pub maximum_step_s: f64,
    /// Explicit source/physical domain for phase-average unsigned speed. Each
    /// phase intersects `dt * this bracket` with the phase distance bracket.
    /// Refining time therefore cannot inflate a distance endpoint into an
    /// out-of-domain speed. This is caller admission, not a free-speed bound.
    pub mean_speed_bracket_m_s: [f64; 2],
    /// Includes rejected comparisons. Each comparison executes three bounded
    /// phase solves; their own source/Newton/quadrature budgets remain active.
    pub maximum_attempts: usize,
    /// Bounds retained accepted phase receipts as well as accepted solve count.
    pub maximum_accepted_phases: usize,
    pub tolerance: DirectionalMarchTolerance,
    /// Independent whole-span energy gate, not multiplied by phase count.
    pub energy_tolerance_j: f64,
    /// Actual total curved path versus every original owned source path change,
    /// including floating-point accumulation and per-phase closure defects.
    pub source_path_tolerance_m: f64,
}

/// Coarse/fine differences in physical units; these are local error estimates.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct DirectionalMarchDefect {
    pub linear_momentum_ns: f64,
    pub angular_momentum_nms: f64,
    pub cell_energy_j: f64,
    pub distance_m: f64,
    pub net_slip_m: f64,
    pub dissipated_j: f64,
    pub normalized: f64,
}

pub struct AcceptedDirectionalHeatPhase<P> {
    pub start_s: f64,
    pub end_s: f64,
    pub report: HeatedDirectionalTrajectoryReport<P>,
}

pub struct DirectionalHeatMarchReport<P> {
    pub phases: Vec<AcceptedDirectionalHeatPhase<P>>,
    pub attempts: usize,
    pub rejected_attempts: usize,
    /// Includes discarded full-step and rejected-comparison solve work.
    pub thermal_trials: usize,
    pub source_evaluations: usize,
    pub maximum_accepted_normalized_defect: f64,
    pub distance_m: f64,
    pub net_slip_m: [f64; 2],
    pub dissipated_j: f64,
    pub body_heat_j: [f64; 2],
    pub external_heat_j: f64,
    pub kinetic_energy_change_j: f64,
    pub thermal_energy_change_j: f64,
    pub mechanical_energy_residual_j: f64,
    pub thermal_energy_residual_j: f64,
    pub combined_energy_residual_j: f64,
    pub source_path_error_m: f64,
}

#[derive(Debug)]
pub enum DirectionalHeatMarchError<E> {
    Invalid(&'static str),
    Phase(DirectionalTrajectoryHeatError<E>),
    Dynamics(DynamicsError),
    Cancelled,
    AttemptBudget,
    AcceptedPhaseBudget,
    Unresolved {
        step_s: f64,
        defect: DirectionalMarchDefect,
    },
    Energy {
        mechanical_j: f64,
        thermal_j: f64,
        combined_j: f64,
    },
    SourcePath {
        error_m: f64,
        tolerance_m: f64,
    },
}
impl<E: core::fmt::Debug> core::fmt::Display for DirectionalHeatMarchError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "directional friction/heat march: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for DirectionalHeatMarchError<E> {}

/// Advance a complete held-geometry span, publishing all original owners once.
///
/// Source callbacks must be pure candidate evaluations, as for the phase API.
/// Their payloads are retained only for accepted half phases. Source, material,
/// numerical, budget or cancellation failure leaves the caller's entire span
/// untouched, even after successful private substeps. No pending event is
/// consumed. All source coordinates keep their declared domain and authority.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn advance_source_directional_heat<'m, P, E>(
    bodies: &mut [RigidBodyState],
    properties: &[MassProperties],
    history: &mut SourceThermalEllipticState,
    law: &SourceThermalEllipticCoulomb,
    contact: RigidDirectionalContact,
    thermal: &mut DepositionSimulation<'m>,
    route: &FrictionHeatRoute,
    material: DirectionalHeatMaterial<'_>,
    duration_s: f64,
    phase_controls: DirectionalHeatSolve,
    relaxation: EllipticRelaxationControls,
    controls: DirectionalMarchControls,
    interface_temperature_guess_k: f64,
    resolve: impl Fn(
        &DepositionSimulation<'m>,
        [EquilibriumPhaseState; 2],
    ) -> Result<DirectionalHeatSource<P>, E>,
    gate: &CancelGate,
) -> Result<DirectionalHeatMarchReport<P>, DirectionalHeatMarchError<E>> {
    let poll = || {
        if gate.is_requested() {
            Err(DirectionalHeatMarchError::Cancelled)
        } else {
            Ok(())
        }
    };
    poll()?;
    let origin = thermal.time_s();
    let target = origin + duration_s;
    let binding = contact.binding;
    let t = controls.tolerance;
    if bodies.len() != properties.len()
        || bodies.len() > phase_controls.mechanical.maximum_bodies
        || binding.body_a >= bodies.len()
        || binding.body_b >= bodies.len()
        || binding.body_a == binding.body_b
        || thermal.cells().len() > phase_controls.maximum_cells
        || thermal.cells().is_empty()
        || !duration_s.is_finite()
        || duration_s <= 0.0
        || !target.is_finite()
        || target <= origin
        || controls.maximum_attempts == 0
        || controls.maximum_accepted_phases < 2
        || [
            controls.initial_step_s,
            controls.minimum_step_s,
            controls.maximum_step_s,
            t.linear_momentum_ns,
            t.angular_momentum_nms,
            t.cell_energy_j,
            t.distance_m,
            t.dissipated_j,
            interface_temperature_guess_k,
        ]
        .into_iter()
        .any(|x| !x.is_finite() || x <= 0.0)
        || [
            t.relative,
            controls.energy_tolerance_j,
            controls.source_path_tolerance_m,
        ]
        .into_iter()
        .any(|x| !x.is_finite() || x < 0.0)
        || controls.minimum_step_s > controls.initial_step_s
        || controls.initial_step_s > controls.maximum_step_s
        || controls
            .mean_speed_bracket_m_s
            .iter()
            .any(|x| !x.is_finite() || *x < 0.0)
        || controls.mean_speed_bracket_m_s[0] > controls.mean_speed_bracket_m_s[1]
        || phase_controls
            .motion
            .distance_bracket_m
            .iter()
            .any(|x| !x.is_finite() || *x < 0.0)
        || phase_controls.motion.distance_bracket_m[0] > phase_controls.motion.distance_bracket_m[1]
    {
        return Err(DirectionalHeatMarchError::Invalid(
            "bounded original owners, time steps and physical error scales required",
        ));
    }
    let initial = [bodies[binding.body_a], bodies[binding.body_b]];
    let props = [properties[binding.body_a], properties[binding.body_b]];
    let mut pair = contact;
    pair.binding.body_a = 0;
    pair.binding.body_b = 1;
    let mut accepted_bodies = initial;
    let mut accepted_history = history.clone();
    let mut accepted_thermal = thermal.clone();
    let zero_energy = vec![0.0; thermal.cells().len()];
    let mut guess_temperature = interface_temperature_guess_k;
    let mut step_s = controls.initial_step_s;
    let mut attempts = 0;
    let mut rejected = 0;
    let mut trials = 0usize;
    let mut evaluations = 0usize;
    let mut maximum_defect = 0.0_f64;
    let mut phases = Vec::new();
    // Compensated totals use actual small transfers, never differences of large
    // aggregate formation energies or reconstructed average-force work.
    let mut distance = Sum::default();
    let mut slip = [Sum::default(); 2];
    let mut heat = Sum::default();
    let mut external = Sum::default();
    let mut body_heat = [Sum::default(); 2];
    while accepted_thermal.time_s() < target {
        poll()?;
        if attempts >= controls.maximum_attempts {
            return Err(DirectionalHeatMarchError::AttemptBudget);
        }
        if phases.len() > controls.maximum_accepted_phases - 2 {
            return Err(DirectionalHeatMarchError::AcceptedPhaseBudget);
        }
        attempts += 1;
        let start = accepted_thermal.time_s();
        let remaining = target - start;
        let planned = step_s.min(remaining);
        // Avoid leaving a sub-ulp or below-minimum tail after repeated rounded
        // steps. Two balanced intervals preserve the exact requested endpoint
        // without silently lengthening a step past its admitted maximum.
        let planned = if remaining > planned && remaining - planned < controls.minimum_step_s {
            0.5 * remaining
        } else {
            planned
        };
        let end = if planned == remaining {
            target
        } else {
            start + planned
        };
        let half = start + 0.5 * (end - start);
        if !end.is_finite() || half <= start || half >= end {
            return Err(DirectionalHeatMarchError::Invalid(
                "comparison times are not representable",
            ));
        }
        let mut coarse_bodies = accepted_bodies;
        let mut coarse_history = accepted_history.clone();
        let mut coarse_thermal = accepted_thermal.clone();
        let mut fine_bodies = accepted_bodies;
        let mut fine_history = accepted_history.clone();
        let mut fine_thermal = accepted_thermal.clone();
        let mut phase = |b: &mut [RigidBodyState; 2],
                         s: &mut SourceThermalEllipticState,
                         h: &mut DepositionSimulation<'m>,
                         stop: f64,
                         ti: f64| {
            let dt = stop - h.time_s();
            let mut local_controls = phase_controls;
            let speed_distance = controls.mean_speed_bracket_m_s.map(|v| v * dt);
            if speed_distance.iter().any(|x| !x.is_finite())
                || controls
                    .mean_speed_bracket_m_s
                    .iter()
                    .zip(speed_distance)
                    .any(|(v, s)| *v > 0.0 && s == 0.0)
            {
                return Err(DirectionalHeatMarchError::Invalid(
                    "source speed bracket has unrepresentable phase travel",
                ));
            }
            local_controls.motion.distance_bracket_m = [
                phase_controls.motion.distance_bracket_m[0].max(speed_distance[0]),
                phase_controls.motion.distance_bracket_m[1].min(speed_distance[1]),
            ];
            let report = apply_source_directional_trajectory_heat(
                b,
                &props,
                s,
                law,
                pair,
                h,
                route,
                material,
                dt,
                local_controls,
                relaxation,
                &zero_energy,
                ti,
                &resolve,
                gate,
            )
            .map_err(DirectionalHeatMarchError::Phase)?;
            if h.time_s() != stop {
                return Err(DirectionalHeatMarchError::Invalid(
                    "phase endpoints disagree with comparison clock",
                ));
            }
            trials = trials.checked_add(report.thermal_trials).ok_or(
                DirectionalHeatMarchError::Invalid("thermal trial count overflow"),
            )?;
            evaluations = evaluations.checked_add(report.source_evaluations).ok_or(
                DirectionalHeatMarchError::Invalid("source evaluation count overflow"),
            )?;
            Ok::<_, DirectionalHeatMarchError<E>>(report)
        };
        let coarse = phase(
            &mut coarse_bodies,
            &mut coarse_history,
            &mut coarse_thermal,
            end,
            guess_temperature,
        )?;
        let first = phase(
            &mut fine_bodies,
            &mut fine_history,
            &mut fine_thermal,
            half,
            guess_temperature,
        )?;
        let second = phase(
            &mut fine_bodies,
            &mut fine_history,
            &mut fine_thermal,
            end,
            first.heat.interface_temperature_k,
        )?;
        let defect = compare(
            accepted_bodies,
            coarse_bodies,
            fine_bodies,
            &accepted_thermal,
            &coarse_thermal,
            &fine_thermal,
            &coarse,
            &first,
            &second,
            t,
            gate,
        )?;
        let dt = end - start;
        if defect.normalized > 1.0 {
            rejected += 1;
            if dt <= controls.minimum_step_s {
                return Err(DirectionalHeatMarchError::Unresolved { step_s: dt, defect });
            }
            let factor = (0.8 / det::sqrt(defect.normalized)).clamp(0.1, 0.5);
            step_s = (dt * factor).max(controls.minimum_step_s);
            continue;
        }
        poll()?;
        maximum_defect = maximum_defect.max(defect.normalized);
        for report in [&first, &second] {
            distance.add(report.actual_distance_m);
            heat.add(report.heat.dissipated_j);
            external.add(report.heat.external_heat_j);
            for i in 0..2 {
                slip[i].add(report.mechanical.relaxation.slip_increment_m[i]);
                body_heat[i].add(report.heat.body_heat_j[i]);
            }
        }
        checked(distance.value)?;
        checked(heat.value)?;
        checked(external.value)?;
        let path_error = source_path_error(law, history, &fine_history, distance.value)?;
        if path_error > controls.source_path_tolerance_m {
            return Err(DirectionalHeatMarchError::SourcePath {
                error_m: path_error,
                tolerance_m: controls.source_path_tolerance_m,
            });
        }
        phases
            .try_reserve(2)
            .map_err(|_| DirectionalHeatMarchError::Invalid("accepted phase allocation refused"))?;
        guess_temperature = second.heat.interface_temperature_k;
        phases.push(AcceptedDirectionalHeatPhase {
            start_s: start,
            end_s: half,
            report: first,
        });
        phases.push(AcceptedDirectionalHeatPhase {
            start_s: half,
            end_s: end,
            report: second,
        });
        accepted_bodies = fine_bodies;
        accepted_history = fine_history;
        accepted_thermal = fine_thermal;
        let factor = if defect.normalized == 0.0 {
            2.0
        } else {
            (0.8 / det::sqrt(defect.normalized)).clamp(0.5, 2.0)
        };
        step_s = (dt * factor)
            .max(controls.minimum_step_s)
            .min(controls.maximum_step_s);
    }
    let mut kinetic = Sum::default();
    for i in 0..2 {
        let old = initial[i]
            .point_kinematics(props[i], Vec3::ZERO)
            .map_err(DirectionalHeatMarchError::Dynamics)?;
        let new = accepted_bodies[i]
            .point_kinematics(props[i], Vec3::ZERO)
            .map_err(DirectionalHeatMarchError::Dynamics)?;
        let dp = accepted_bodies[i]
            .linear_momentum_world()
            .sub(initial[i].linear_momentum_world());
        let dl = accepted_bodies[i]
            .angular_momentum_body()
            .sub(initial[i].angular_momentum_body());
        kinetic.add(
            dp.dot(
                old.center_of_mass_velocity_world
                    .scale(0.5)
                    .add(new.center_of_mass_velocity_world.scale(0.5)),
            ),
        );
        kinetic.add(
            dl.dot(
                old.angular_velocity_body
                    .scale(0.5)
                    .add(new.angular_velocity_body.scale(0.5)),
            ),
        );
    }
    let mut stored = Sum::default();
    for (index, (id, old)) in thermal.cells().iter().enumerate() {
        if index % 256 == 0 {
            poll()?;
        }
        let new = &accepted_thermal.cells()[id];
        stored.add(old.mass_kg * (new.specific_h_j_kg - old.specific_h_j_kg));
    }
    let mechanical_j = checked(kinetic.value + heat.value)?;
    let thermal_j = checked(stored.value - external.value - heat.value)?;
    let combined_j = checked(kinetic.value + stored.value - external.value)?;
    if [mechanical_j, thermal_j, combined_j]
        .into_iter()
        .any(|x| x.abs() > controls.energy_tolerance_j)
    {
        return Err(DirectionalHeatMarchError::Energy {
            mechanical_j,
            thermal_j,
            combined_j,
        });
    }
    let source_path_error_m = source_path_error(law, history, &accepted_history, distance.value)?;
    let report = DirectionalHeatMarchReport {
        phases,
        attempts,
        rejected_attempts: rejected,
        thermal_trials: trials,
        source_evaluations: evaluations,
        maximum_accepted_normalized_defect: maximum_defect,
        distance_m: distance.value,
        net_slip_m: [checked(slip[0].value)?, checked(slip[1].value)?],
        dissipated_j: heat.value,
        body_heat_j: [checked(body_heat[0].value)?, checked(body_heat[1].value)?],
        external_heat_j: external.value,
        kinetic_energy_change_j: checked(kinetic.value)?,
        thermal_energy_change_j: checked(stored.value)?,
        mechanical_energy_residual_j: mechanical_j,
        thermal_energy_residual_j: thermal_j,
        combined_energy_residual_j: combined_j,
        source_path_error_m,
    };
    poll()?;
    bodies[binding.body_a] = accepted_bodies[0];
    bodies[binding.body_b] = accepted_bodies[1];
    *history = accepted_history;
    *thermal = accepted_thermal;
    Ok(report)
}

#[derive(Clone, Copy, Default)]
struct Sum {
    value: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) {
        let y = value - self.correction;
        let next = self.value + y;
        self.correction = (next - self.value) - y;
        self.value = next;
    }
}
fn checked<E>(x: f64) -> Result<f64, DirectionalHeatMarchError<E>> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(DirectionalHeatMarchError::Invalid(
            "unrepresentable march coordinate or account",
        ))
    }
}
fn norm(v: Vec3) -> f64 {
    det::hypot(det::hypot(v.x, v.y), v.z)
}
fn norm2(v: [f64; 2]) -> f64 {
    det::hypot(v[0], v[1])
}

#[allow(clippy::too_many_arguments)]
fn compare<P, E>(
    old: [RigidBodyState; 2],
    coarse: [RigidBodyState; 2],
    fine: [RigidBodyState; 2],
    old_thermal: &DepositionSimulation<'_>,
    coarse_thermal: &DepositionSimulation<'_>,
    fine_thermal: &DepositionSimulation<'_>,
    c: &HeatedDirectionalTrajectoryReport<P>,
    a: &HeatedDirectionalTrajectoryReport<P>,
    b: &HeatedDirectionalTrajectoryReport<P>,
    tolerance: DirectionalMarchTolerance,
    gate: &CancelGate,
) -> Result<DirectionalMarchDefect, DirectionalHeatMarchError<E>> {
    let mut out = DirectionalMarchDefect::default();
    let mut scale = |error: f64,
                     absolute: f64,
                     relative_scale: f64|
     -> Result<(), DirectionalHeatMarchError<E>> {
        let denominator = checked(absolute + tolerance.relative * relative_scale)?;
        out.normalized = out.normalized.max(checked(error / denominator)?);
        Ok(())
    };
    for i in 0..2 {
        let cp = coarse[i]
            .linear_momentum_world()
            .sub(old[i].linear_momentum_world());
        let fp = fine[i]
            .linear_momentum_world()
            .sub(old[i].linear_momentum_world());
        let ep = checked(norm(
            coarse[i]
                .linear_momentum_world()
                .sub(fine[i].linear_momentum_world()),
        ))?;
        out.linear_momentum_ns = out.linear_momentum_ns.max(ep);
        scale(ep, tolerance.linear_momentum_ns, norm(cp).max(norm(fp)))?;
        let cl = coarse[i]
            .angular_momentum_body()
            .sub(old[i].angular_momentum_body());
        let fl = fine[i]
            .angular_momentum_body()
            .sub(old[i].angular_momentum_body());
        let el = checked(norm(
            coarse[i]
                .angular_momentum_body()
                .sub(fine[i].angular_momentum_body()),
        ))?;
        out.angular_momentum_nms = out.angular_momentum_nms.max(el);
        scale(el, tolerance.angular_momentum_nms, norm(cl).max(norm(fl)))?;
    }
    for (index, (id, old)) in old_thermal.cells().iter().enumerate() {
        if index % 256 == 0 && gate.is_requested() {
            return Err(DirectionalHeatMarchError::Cancelled);
        }
        let c = &coarse_thermal.cells()[id];
        let f = &fine_thermal.cells()[id];
        let error = checked((old.mass_kg * (c.specific_h_j_kg - f.specific_h_j_kg)).abs())?;
        out.cell_energy_j = out.cell_energy_j.max(error);
        scale(
            error,
            tolerance.cell_energy_j,
            (old.mass_kg * (c.specific_h_j_kg - old.specific_h_j_kg))
                .abs()
                .max((old.mass_kg * (f.specific_h_j_kg - old.specific_h_j_kg)).abs()),
        )?;
    }
    let fine_path = checked(a.actual_distance_m + b.actual_distance_m)?;
    out.distance_m = checked((c.actual_distance_m - fine_path).abs())?;
    scale(
        out.distance_m,
        tolerance.distance_m,
        c.actual_distance_m.max(fine_path),
    )?;
    let coarse_slip = c.mechanical.relaxation.slip_increment_m;
    let fine_slip = core::array::from_fn(|i| {
        a.mechanical.relaxation.slip_increment_m[i] + b.mechanical.relaxation.slip_increment_m[i]
    });
    out.net_slip_m = checked(norm2(core::array::from_fn(|i| {
        coarse_slip[i] - fine_slip[i]
    })))?;
    scale(
        out.net_slip_m,
        tolerance.distance_m,
        norm2(coarse_slip).max(norm2(fine_slip)),
    )?;
    let fine_heat = checked(a.heat.dissipated_j + b.heat.dissipated_j)?;
    out.dissipated_j = checked((c.heat.dissipated_j - fine_heat).abs())?;
    scale(
        out.dissipated_j,
        tolerance.dissipated_j,
        c.heat.dissipated_j.max(fine_heat),
    )?;
    Ok(out)
}

fn source_path_error<E>(
    law: &SourceThermalEllipticCoulomb,
    initial: &SourceThermalEllipticState,
    next: &SourceThermalEllipticState,
    distance: f64,
) -> Result<f64, DirectionalHeatMarchError<E>> {
    let mut error = 0.0_f64;
    for i in 0..2 {
        for axis in law.principal_sources()[i].query_axes() {
            if axis.kind == FrictionQueryAxisKind::AccumulatedDistance {
                let old = initial.friction_sources()[i].query().axes()[&axis.name];
                let new = next.friction_sources()[i].query().axes()[&axis.name];
                error = error.max(checked((new - old - distance).abs())?);
            }
        }
    }
    Ok(error)
}
