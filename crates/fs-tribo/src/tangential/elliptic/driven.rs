//! Held-load directional Coulomb flow with an explicit external drive.
//!
//! The effective equation is `dv/dt = a - W F(v)`, with the original full SPD
//! body mobility W. A separately retained resting mode uses the STATIC ellipse
//! to hold `F=W^-1 a`; a moving mode uses kinetic friction until a resolved
//! arrest. Constant-direction segments have exact breakaway/stop/restart times.
//! Other sliding segments use the original native matrix projection to impose
//! Coulomb at the endpoint, then hold that traction over a linear-velocity path.
//! This changing-direction rung is a bounded first-order time discretization.
//! Its endpoint arrests are timestep-resolved, not exact continuous events.
//!
//! Every held traction opposes both endpoints, hence the whole linear path.
//! Its physical heat is traction times actual midpoint displacement, not the
//! endpoint-impulse loss of a fictitious force kick. No numerical loss is added
//! to heat. Source/load/geometry changes require an outer refined transaction.

use super::relaxation::{
    EllipticRelaxationControls, EllipticRelaxationError, EllipticRelaxationInput,
};
use super::*;
use crate::TangentialSolveConfig;

/// Retained physical contact mode, never inferred from a small sliding speed.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EllipticDrivenMode {
    Resting,
    Sliding,
}

#[derive(Clone, Copy, Debug)]
pub struct EllipticDrivenInput {
    pub velocity_m_s: [f64; 2],
    pub mobility_per_kg: [[f64; 2]; 2],
    /// Tangent acceleration from the original external body wrenches alone.
    pub driving_acceleration_m_s2: [f64; 2],
    pub mode: EllipticDrivenMode,
    pub normal_n: f64,
    pub dt_s: f64,
    pub maximum_force_n: f64,
    pub contact_temperature_k: Option<f64>,
    pub nominal_pressure_pa: Option<f64>,
}

#[derive(Clone, Copy, Debug)]
pub struct EllipticDrivenControls {
    /// Accuracy axis for changing-direction endpoint-Coulomb substeps.
    pub maximum_step_s: f64,
    /// Refuse rather than hiding an unresolved direction/event below this step.
    pub minimum_step_s: f64,
    /// Includes rejected/refined projection attempts, not just accepted steps.
    pub maximum_substeps: usize,
    /// Only validates an explicitly retained Resting constraint. Sliding never
    /// becomes Resting through this tolerance. Actual body momenta are not reset.
    pub rest_velocity_tolerance_m_s: f64,
    pub solve: TangentialSolveConfig,
    /// Original zero-drive solver and per-segment unsigned-path quadrature.
    pub relaxation: EllipticRelaxationControls,
}

#[derive(Clone, Debug, PartialEq)]
pub struct EllipticDrivenReport {
    pub state: EllipticCoulombState,
    pub mode: EllipticDrivenMode,
    pub velocity_after_m_s: [f64; 2],
    pub driving_force_n: [f64; 2],
    /// Positive traction impulse; the force on body A has the opposite sign.
    pub traction_impulse_ns: [f64; 2],
    pub slip_increment_m: [f64; 2],
    pub distance_m: f64,
    pub dissipated_j: f64,
    pub drive_work_j: f64,
    pub kinetic_energy_change_j: f64,
    pub energy_residual_j: f64,
    /// First transition from nonzero sliding velocity to rest in this phase.
    pub first_arrest_time_s: Option<f64>,
    /// Meaningful only when first_arrest_time_s is Some. False means an
    /// endpoint-resolved event of the changing-direction time discretization.
    pub arrest_time_is_exact: bool,
    pub substeps: usize,
    pub projection_iterations: usize,
    pub path_quadrature_error_m: f64,
    pub quadrature_evaluations: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub enum EllipticDrivenError {
    Constitutive(TriboError),
    Relaxation(EllipticRelaxationError),
    Invalid(&'static str),
    Cancelled,
    SubstepBudget,
    UnresolvedDirection,
}
impl From<TriboError> for EllipticDrivenError {
    fn from(e: TriboError) -> Self {
        Self::Constitutive(e)
    }
}
impl From<EllipticRelaxationError> for EllipticDrivenError {
    fn from(e: EllipticRelaxationError) -> Self {
        Self::Relaxation(e)
    }
}
impl core::fmt::Display for EllipticDrivenError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "driven directional friction: {self:?}")
    }
}
impl std::error::Error for EllipticDrivenError {}
type Error = EllipticDrivenError;

impl EllipticCoulomb {
    /// Integrate one held-source/geometry/load phase under a held external drive.
    ///
    /// This borrows the original state. Resting is an explicit physical mode,
    /// with its actual body-velocity residual bounded independently by controls.
    /// It is not a small-speed capture rule. The body owner must carry the
    /// returned mode and independently audit actual endpoint momentum/velocity.
    /// General changing-direction motion requires refinement of maximum_step_s;
    /// exact zero-drive relaxation and constant-direction limits are retained.
    #[allow(clippy::too_many_lines)]
    pub fn advance_driven(
        &self,
        state: &EllipticCoulombState,
        input: &EllipticDrivenInput,
        controls: EllipticDrivenControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<EllipticDrivenReport, Error> {
        poll(cancelled)?;
        self.check(state)?;
        controls.solve.validate()?;
        let r = controls.relaxation;
        if [
            input.dt_s,
            controls.maximum_step_s,
            controls.minimum_step_s,
            r.time_tolerance_s,
            r.travel_tolerance_m,
        ]
        .into_iter()
        .any(|x| !x.is_finite() || x <= 0.0)
            || controls.minimum_step_s > controls.maximum_step_s
            || controls.maximum_substeps == 0
            || r.maximum_evaluations == 0
            || r.maximum_root_iterations == 0
            || [
                input.normal_n,
                input.maximum_force_n,
                controls.rest_velocity_tolerance_m_s,
            ]
            .into_iter()
            .any(|x| !x.is_finite() || x < 0.0)
            || input
                .velocity_m_s
                .into_iter()
                .chain(input.driving_acceleration_m_s2)
                .any(|x| !x.is_finite())
            || input
                .contact_temperature_k
                .is_some_and(|x| !x.is_finite() || x <= 0.0)
            || input
                .nominal_pressure_pa
                .is_some_and(|x| !x.is_finite() || x < 0.0)
        {
            return Err(Error::Invalid(
                "finite driven phase, tolerances and work bounds required",
            ));
        }
        let w = input.mobility_per_kg;
        let matrix = SymmetricTangentialCompliance::new(w)?;
        let (normalized, determinant) = matrix.material([1.0, 0.0])?;
        let scale = w[0][0].max(w[1][1]);
        if determinant <= 0.0 || scale <= 0.0 {
            return Err(Error::Invalid("driven motion requires SPD body mobility"));
        }
        let inverse = |v: [f64; 2]| -> Result<[f64; 2], Error> {
            let [a, b, d] = normalized;
            let g = v.map(|x| x / scale);
            vector([
                (d * g[0] - b * g[1]) / determinant,
                (a * g[1] - b * g[0]) / determinant,
            ])
        };
        let drive = vector(input.driving_acceleration_m_s2)?;
        let q = inverse(drive)?;
        let static_caps = capacities(self.static_mu, input.normal_n)?;
        let caps = capacities(self.kinetic_mu, input.normal_n)?;
        let q_material = to_principal(state.axis, q)?;
        let can_hold = inside(q_material, static_caps);
        let kinetic_hold = strictly_inside(q_material, caps);
        let force_limit = input.maximum_force_n.min(controls.solve.maximum_force_n);
        let mut mode = input.mode;
        let initial = match mode {
            EllipticDrivenMode::Resting => {
                if length(input.velocity_m_s) > controls.rest_velocity_tolerance_m_s {
                    return Err(Error::Invalid(
                        "retained rest exceeds actual velocity residual gate",
                    ));
                }
                [0.0; 2]
            }
            EllipticDrivenMode::Sliding => {
                if input.velocity_m_s == [0.0; 2] {
                    return Err(Error::Invalid(
                        "zero velocity requires explicit resting mode",
                    ));
                }
                input.velocity_m_s
            }
        };
        let mut report = EllipticDrivenReport {
            state: state.clone(),
            mode,
            velocity_after_m_s: initial,
            driving_force_n: q,
            traction_impulse_ns: [0.0; 2],
            slip_increment_m: [0.0; 2],
            distance_m: 0.0,
            dissipated_j: 0.0,
            drive_work_j: 0.0,
            kinetic_energy_change_j: 0.0,
            energy_residual_j: 0.0,
            first_arrest_time_s: None,
            arrest_time_is_exact: true,
            substeps: 0,
            projection_iterations: 0,
            path_quadrature_error_m: 0.0,
            quadrature_evaluations: 0,
        };
        if drive == [0.0; 2] {
            let native = self.relax_velocity(
                state,
                &EllipticRelaxationInput {
                    velocity_m_s: initial,
                    mobility_per_kg: w,
                    normal_n: input.normal_n,
                    dt_s: input.dt_s,
                    maximum_force_n: force_limit,
                    contact_temperature_k: input.contact_temperature_k,
                    nominal_pressure_pa: input.nominal_pressure_pa,
                },
                r,
                cancelled,
            )?;
            report.mode = if native.velocity_after_m_s == [0.0; 2] {
                EllipticDrivenMode::Resting
            } else {
                EllipticDrivenMode::Sliding
            };
            report.first_arrest_time_s =
                if initial != [0.0; 2] && report.mode == EllipticDrivenMode::Resting {
                    native.kinetic_phase_end_s
                } else {
                    None
                };
            report.state = native.state;
            report.velocity_after_m_s = native.velocity_after_m_s;
            report.traction_impulse_ns = native.traction_impulse_ns;
            report.slip_increment_m = native.slip_increment_m;
            report.distance_m = native.distance_m;
            report.dissipated_j = native.dissipated_j;
            report.path_quadrature_error_m = native.quadrature_error[3];
            report.quadrature_evaluations = native.evaluations;
            report.substeps = 1;
        } else {
            let mut velocity = initial;
            let mut elapsed = 0.0;
            while elapsed < input.dt_s {
                poll(cancelled)?;
                if report.substeps >= controls.maximum_substeps {
                    return Err(Error::SubstepBudget);
                }
                let remaining = input.dt_s - elapsed;
                if mode == EllipticDrivenMode::Resting && can_hold {
                    if length(q) > force_limit {
                        return Err(Error::Invalid("static reaction exceeds force budget"));
                    }
                    for i in 0..2 {
                        report.traction_impulse_ns[i] =
                            checked(report.traction_impulse_ns[i] + remaining * q[i])?;
                    }
                    report.substeps += 1;
                    break;
                }
                if caps[0].max(caps[1]) > force_limit {
                    return Err(Error::Invalid("kinetic capacity exceeds force budget"));
                }
                let mut h = remaining;
                let mut exact = true;
                let (force, after, path, path_error, evaluations) = if mode
                    == EllipticDrivenMode::Resting
                {
                    // Homogeneity makes the outgoing acceleration constant.
                    // The native endpoint equation finds its associated support.
                    report.substeps += 1;
                    let (force, end, iterations) =
                        project(self, state, input, controls, [0.0; 2], h, force_limit)?;
                    report.projection_iterations += iterations;
                    if end == [0.0; 2] {
                        return Err(Error::Invalid("superstatic drive did not break away"));
                    }
                    (force, end, 0.5 * h * length(end), 0.0, 0)
                } else {
                    let force = from_principal(
                        state.axis,
                        support(caps, to_principal(state.axis, velocity)?)?,
                    )?;
                    let acceleration = vector([
                        drive[0] - multiply(w, force)[0],
                        drive[1] - multiply(w, force)[1],
                    ])?;
                    let speed = length(velocity);
                    let along = dot(velocity.map(|x| x / speed), acceleration);
                    let cross = if length(acceleration) == 0.0 {
                        0.0
                    } else {
                        let a = acceleration.map(|x| x / length(acceleration));
                        (velocity[0] / speed) * a[1] - (velocity[1] / speed) * a[0]
                    };
                    if cross.abs() <= 256.0 * f64::EPSILON {
                        report.substeps += 1;
                        let stop = if along < 0.0 {
                            speed / -along
                        } else {
                            f64::INFINITY
                        };
                        h = h.min(stop);
                        let end = if h == stop {
                            [0.0; 2]
                        } else {
                            vector([
                                velocity[0] + h * acceleration[0],
                                velocity[1] + h * acceleration[1],
                            ])?
                        };
                        (force, end, 0.5 * h * (speed + length(end)), 0.0, 0)
                    } else {
                        exact = false;
                        h = h.min(controls.maximum_step_s);
                        let (force, end) = loop {
                            poll(cancelled)?;
                            if report.substeps >= controls.maximum_substeps {
                                return Err(Error::SubstepBudget);
                            }
                            report.substeps += 1;
                            if elapsed + h <= elapsed || h <= 0.0 {
                                return Err(Error::UnresolvedDirection);
                            }
                            let (force, end, iterations) =
                                project(self, state, input, controls, velocity, h, force_limit)?;
                            report.projection_iterations += iterations;
                            // Outside the kinetic ellipse, a coarse endpoint
                            // projection cannot by itself authorize capture in
                            // the static/kinetic hysteresis band. Refine until
                            // the moving branch or constant-direction event is
                            // resolved; exhaustion remains a refusal. This is
                            // not a theorem excluding all anisotropic arrests.
                            if dot(force, velocity) < 0.0 || (end == [0.0; 2] && !kinetic_hold) {
                                h *= 0.5;
                                if h < controls.minimum_step_s {
                                    return Err(Error::UnresolvedDirection);
                                }
                                continue;
                            }
                            break (force, end);
                        };
                        let (path, error, evaluations) = if end == [0.0; 2] {
                            (0.5 * h * speed, 0.0, 0)
                        } else {
                            super::relaxation::linear_velocity_path(velocity, end, h, r, cancelled)?
                        };
                        (force, end, path, error, evaluations)
                    }
                };
                if elapsed + h <= elapsed || !h.is_finite() {
                    return Err(Error::UnresolvedDirection);
                }
                let displacement = vector([
                    0.5 * h * (velocity[0] + after[0]),
                    0.5 * h * (velocity[1] + after[1]),
                ])?;
                let loss = checked(dot(force, displacement))?;
                if loss < 0.0 || !path.is_finite() || path < 0.0 {
                    return Err(Error::Invalid(
                        "negative or unrepresentable driven work/path",
                    ));
                }
                for i in 0..2 {
                    report.traction_impulse_ns[i] =
                        checked(report.traction_impulse_ns[i] + h * force[i])?;
                    report.slip_increment_m[i] =
                        checked(report.slip_increment_m[i] + displacement[i])?;
                }
                report.distance_m = checked(report.distance_m + path)?;
                report.dissipated_j = checked(report.dissipated_j + loss)?;
                report.path_quadrature_error_m =
                    checked(report.path_quadrature_error_m + path_error)?;
                report.quadrature_evaluations = report
                    .quadrature_evaluations
                    .checked_add(evaluations)
                    .ok_or(Error::Invalid("quadrature count overflow"))?;
                elapsed = if h == remaining {
                    input.dt_s
                } else {
                    checked(elapsed + h)?
                };
                if after == [0.0; 2] {
                    mode = EllipticDrivenMode::Resting;
                    if report.first_arrest_time_s.is_none() {
                        report.first_arrest_time_s = Some(elapsed);
                        report.arrest_time_is_exact = exact;
                    }
                } else {
                    mode = EllipticDrivenMode::Sliding;
                }
                velocity = after;
            }
            report.mode = mode;
            report.velocity_after_m_s = velocity;
            report.state.contact_temperature_k = input.contact_temperature_k;
            for i in 0..2 {
                report.state.slip_m[i] = checked(state.slip_m[i] + report.slip_increment_m[i])?;
            }
        }
        let delta = inverse([
            report.velocity_after_m_s[0] - initial[0],
            report.velocity_after_m_s[1] - initial[1],
        ])?;
        report.kinetic_energy_change_j = checked(dot(
            delta,
            [
                0.5 * (initial[0] + report.velocity_after_m_s[0]),
                0.5 * (initial[1] + report.velocity_after_m_s[1]),
            ],
        ))?;
        report.drive_work_j = checked(dot(q, report.slip_increment_m))?;
        report.energy_residual_j =
            checked(report.kinetic_energy_change_j + report.dissipated_j - report.drive_work_j)?;
        poll(cancelled)?;
        Ok(report)
    }
}

#[allow(clippy::too_many_arguments)]
fn project(
    law: &EllipticCoulomb,
    state: &EllipticCoulombState,
    input: &EllipticDrivenInput,
    controls: EllipticDrivenControls,
    velocity: [f64; 2],
    h: f64,
    force_limit: f64,
) -> Result<([f64; 2], [f64; 2], usize), Error> {
    let h2 = checked(h * h)?;
    if h2 == 0.0 {
        return Err(Error::Invalid("unrepresentable driven compliance"));
    }
    let free = vector([
        h * (velocity[0] + h * input.driving_acceleration_m_s2[0]),
        h * (velocity[1] + h * input.driving_acceleration_m_s2[1]),
    ])?;
    let matrix =
        SymmetricTangentialCompliance::new(input.mobility_per_kg.map(|row| row.map(|x| h2 * x)))?;
    let mut solve = controls.solve;
    solve.maximum_force_n = force_limit;
    let step = law.advance_kinetic_with_compliance(
        state,
        &TangentialStepInput {
            normal_n: input.normal_n,
            slip_increment_m: free,
            dt_s: h,
            compliance_m_per_n: 0.0,
            contact_temperature_k: input.contact_temperature_k,
            nominal_pressure_pa: input.nominal_pressure_pa,
            solve,
        },
        matrix,
    )?;
    Ok((
        step.traction_n,
        vector(step.slip_increment_m.map(|x| x / h))?,
        step.iterations,
    ))
}
fn checked(x: f64) -> Result<f64, Error> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(Error::Invalid("unrepresentable driven quantity"))
    }
}
fn vector(x: [f64; 2]) -> Result<[f64; 2], Error> {
    for value in x {
        checked(value)?;
    }
    Ok(x)
}
fn length(x: [f64; 2]) -> f64 {
    det::hypot(x[0], x[1])
}
fn dot(x: [f64; 2], y: [f64; 2]) -> f64 {
    x[0] * y[0] + x[1] * y[1]
}
fn multiply(w: [[f64; 2]; 2], x: [f64; 2]) -> [f64; 2] {
    [dot(w[0], x), dot(w[1], x)]
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
fn strictly_inside(force: [f64; 2], caps: [f64; 2]) -> bool {
    let mut scaled = [0.0; 2];
    for i in 0..2 {
        if caps[i] == 0.0 {
            if force[i] != 0.0 {
                return false;
            }
        } else {
            scaled[i] = force[i] / caps[i];
        }
    }
    length(scaled) < 1.0
}
