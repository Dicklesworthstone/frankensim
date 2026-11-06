//! Full-tensor viscous torque on a rigid body, with prepared implicit stages.
//!
//! The body-frame law is tau = drive - D (omega - Omega_bath). Two implicit
//! midpoint torque stages surround the existing compensated free rigid-body
//! flow. The stationary entry sets Omega_bath = 0; the rotating-bath entry
//! holds its prescribed spin in space and reports its signed source work.
//! Gyroscopic transport stays in its original owner. Relative-slip heat,
//! applied work and the free flow's numerical energy defect remain separate.
//!
//! The Euler equation I omega' + omega cross (I omega) = tau follows
//! Lynch and Park, Modern Robotics, section 8.2:
//! <https://modernrobotics.northwestern.edu/nu-gm-book-resource/8-2-dynamics-of-a-single-rigid-body-part-1-of-2/>.

use super::{LieStepError, compensated, gyrostat::HeldTorque, gyrostat::InertiaTensor3, held};
use fs_ga::{Mat3, So3, Vec3};
use fs_la::factor::{Cholesky, cholesky};
use fs_math::dd::Dd;

/// Finite-inertia reciprocal coupling between two freely rotating bodies.
#[path = "viscous/pair.rs"]
pub mod pair;

type Vector = [Dd; 3];

/// A fixed full inertia, body-frame viscous tensor and positive time interval.
///
/// Preparation admits symmetric positive-semidefinite D and caches ONE native
/// 3-by-3 Cholesky factor of I + dt D / 4. Subsequent steps are allocation-free
/// and use two Cholesky solves; no per-step eigensolve or numerical Jacobian
/// is introduced. Zero drag needs no factorization.
///
/// This is linear viscous drag relative to a stationary or prescribed rotating
/// environment, not dry friction, quadratic aerodynamic drag or a general
/// feedback controller. D is fixed in the body's material axes; a changing
/// external drag tensor requires a different model/preparation.
#[derive(Debug, Clone)]
pub struct ViscousRigidBodyStepper {
    inertia: InertiaTensor3,
    damping_body: Mat3,
    duration_s: f64,
    half_s: f64,
    factor: Option<Cholesky>,
    factor_scale: f64,
    relaxation_rate_bound_per_s: f64,
}

/// Accepted original-body motion with separate drive, drag and numerical work.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ViscousRigidBodyStep {
    /// Accepted body-to-space attitude.
    pub rotation_after: So3,
    /// Accepted angular momentum about the COM, in material axes [kg m^2/s].
    /// This is the authoritative momentum to pass to the next step.
    pub momentum_after_body: Vec3,
    /// Diagnostic I^-1 momentum_after_body in material axes [rad/s].
    pub angular_velocity_after_body: Vec3,
    /// Sum of the two actual endpoint drive impulses in space [N m s].
    pub drive_angular_impulse_space: Vec3,
    /// Sum of the two viscous reaction impulses in space [N m s].
    pub drag_angular_impulse_space: Vec3,
    /// Actual spatial momentum change minus drive and drag impulses [N m s].
    /// Includes native group arithmetic and final public-state rounding.
    pub angular_momentum_residual_space: Vec3,
    /// Signed drive work at the two implicit-stage mean angular rates [J].
    pub drive_work_j: f64,
    /// Nonnegative sum(dt_stage * relative_spin_mid^T D relative_spin_mid) [J].
    /// Relative spin is omega for the stationary entry and omega - Omega_bath
    /// for the rotating-bath entry. This is mechanical loss, with no automatic
    /// thermal destination. Numerical and free-flow errors are never heat.
    pub dissipated_energy_j: f64,
    /// Actual kinetic change between the supplied and returned momentum [J].
    pub kinetic_energy_change_j: f64,
    /// Kinetic change during the native torque-free drift alone [J].
    pub drift_energy_defect_j: f64,
    /// Kinetic change - drive work - bath work + heat - drift defect [J].
    /// Bath work is zero for `step`; `step_in_rotating_bath` reports it in
    /// the enclosing receipt. Retains final-state rounding; no projection.
    pub work_energy_residual_j: f64,
}

/// A body coupled to an externally prescribed, spatially held bath spin.
///
/// The bath receives the opposite of `motion.drag_angular_impulse_space`.
/// Its maintaining actuator supplies `bath_work_j`; a negative value means
/// it extracts work. The body's drag work is bath work minus relative-slip
/// heat. This is an infinite prescribed source, not a finite-inertia body.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ViscousRotatingBathStep {
    /// Original body motion and balances, including bath work in the residual.
    pub motion: ViscousRigidBodyStep,
    /// Sum(tau_drag_body dot Omega_bath_body * dt_stage) [J].
    pub bath_work_j: f64,
}

impl InertiaTensor3 {
    /// Prepare full-tensor body drag for one positive held interval.
    ///
    /// The largest viscous relaxation rate is bounded numerically by
    /// tr(I^-1 D), the trace of the congruent positive-semidefinite form.
    /// Preparation refuses dt/2 * bound > 2, avoiding unresolved midpoint
    /// sign reversals. Midpoint is passive but not L-stable; this admission
    /// does not replace timestep refinement for decay, phase or tumbling.
    ///
    /// Roundoff-sized asymmetry (64 eps of the largest entry) is averaged.
    /// Negative compensated principal minors refuse; D is never projected onto a
    /// different positive tensor. This is floating-point admission, not an
    /// interval certificate. Reprepare when I, D or dt changes.
    pub fn viscous_stepper(
        self,
        damping_body: Mat3,
        duration_s: f64,
    ) -> Result<ViscousRigidBodyStepper, LieStepError> {
        scalar(duration_s, "viscous rigid-body duration")?;
        let half_s = 0.5 * duration_s;
        if duration_s <= 0.0 || half_s == 0.0 {
            return Err(invalid("positive resolved viscous rigid-body duration"));
        }
        let damping_body = admit_damping(damping_body)?;
        let any_drag = damping_body.max_abs() != 0.0;
        let mut trace = Dd::ZERO;
        if any_drag {
            for i in 0..3 {
                let inverse_column = self.velocity(Vec3::new(
                    damping_body.m[i],
                    damping_body.m[3 + i],
                    damping_body.m[6 + i],
                ))?;
                trace = trace + dd(components(inverse_column)[i]);
            }
        }
        let rate = scalar(
            trace.to_f64() * (1.0 + 256.0 * f64::EPSILON),
            "viscous rigid-body relaxation bound",
        )?;
        if rate < 0.0 || (any_drag && rate == 0.0) {
            return Err(invalid("positive viscous rigid-body relaxation mobility"));
        }
        if scalar(half_s * rate, "viscous rigid-body relaxation number")? > 2.0 {
            return Err(invalid(
                "viscous rigid-body interval does not resolve relaxation",
            ));
        }
        let mut matrix = self.matrix().m;
        for (entry, damping) in matrix.iter_mut().zip(damping_body.m) {
            *entry = scalar(
                (dd(*entry) + dd(half_s) * dd(0.5) * dd(damping)).to_f64(),
                "viscous rigid-body factor",
            )?;
        }
        let factor_scale = matrix.iter().map(|x| x.abs()).fold(0.0_f64, f64::max);
        let factor = if any_drag {
            let normalized = matrix.map(|x| x / factor_scale);
            Some(
                cholesky(&normalized, 3)
                    .map_err(|_| invalid("viscous rigid-body factor not SPD"))?,
            )
        } else {
            None
        };
        Ok(ViscousRigidBodyStepper {
            inertia: self,
            damping_body,
            duration_s,
            half_s,
            factor,
            factor_scale,
            relaxation_rate_bound_per_s: rate,
        })
    }
}

impl ViscousRigidBodyStepper {
    /// Positive prepared interval [s].
    pub const fn duration_s(&self) -> f64 {
        self.duration_s
    }

    /// Numerical upper relaxation-rate estimate used at preparation [1/s].
    pub const fn relaxation_rate_bound_per_s(&self) -> f64 {
        self.relaxation_rate_bound_per_s
    }

    /// Admitted symmetric tensor in the original material axes [N m s/rad].
    pub const fn damping_body(&self) -> Mat3 {
        self.damping_body
    }

    /// Advance angular momentum under held drive and implicit viscous drag.
    ///
    /// Drive torque can be held in body or spatial coordinates. A spatial
    /// drive is transformed at BOTH actual endpoint attitudes. The viscous
    /// torque uses each implicit stage's actual mean angular velocity, never
    /// the stale initial velocity. Work pairs those same efforts and rates.
    ///
    /// The symmetric splitting is second order for this smooth held law.
    /// Torque impulse tails remain in DD through native free rotation.
    /// The fixed number of 3-vector operations has no unbounded iteration.
    /// All state is supplied by value; every refusal leaves it unchanged.
    /// This entry uses a stationary bath. No contact, changing inertia,
    /// internal rotor or adjoint is inferred. Residuals are not error bounds.
    pub fn step(
        &self,
        rotation: So3,
        momentum_body: Vec3,
        drive: HeldTorque,
    ) -> Result<ViscousRigidBodyStep, LieStepError> {
        Ok(self.advance(rotation, momentum_body, drive, None)?.motion)
    }

    /// Couple the body to a bath with constant spatial angular velocity [rad/s].
    ///
    /// Each implicit stage uses tau_drag = -D (omega_mid - R^-1 Omega_space).
    /// The prescribed bath spin is transformed at BOTH accepted endpoint
    /// attitudes. Its work is tau_drag dot Omega_bath times stage duration;
    /// the nonnegative relative-slip quadratic form is the dissipated energy.
    /// A faster bath can therefore accelerate the body while producing heat.
    ///
    /// Uses the same cached factor and symmetric second-order composition.
    /// Exactly zero bath spin or zero damping takes the stationary arithmetic.
    /// Body and bath spin must be resolved by timestep refinement, in addition
    /// to the preparation-time relaxation limit. No finite bath inertia,
    /// changing bath law, translational drag, heat partition or adjoint is
    /// inferred. All supplied state remains unchanged on refusal.
    pub fn step_in_rotating_bath(
        &self,
        rotation: So3,
        momentum_body: Vec3,
        drive: HeldTorque,
        bath_angular_velocity_space: Vec3,
    ) -> Result<ViscousRotatingBathStep, LieStepError> {
        let bath = held::finite_vector(
            bath_angular_velocity_space,
            "viscous spatial bath angular velocity",
        )?;
        let bath = if self.factor.is_none()
            || components(bath_angular_velocity_space)
                .iter()
                .all(|x| *x == 0.0)
        {
            None
        } else {
            Some(bath)
        };
        self.advance(rotation, momentum_body, drive, bath)
    }

    fn advance(
        &self,
        rotation: So3,
        momentum_body: Vec3,
        drive: HeldTorque,
        bath_space: Option<Vector>,
    ) -> Result<ViscousRotatingBathStep, LieStepError> {
        let initial = held::finite_vector(momentum_body, "viscous rigid-body momentum")?;
        let first = self.stage(
            initial,
            drive_at(rotation, drive)?,
            bath_at(rotation, bath_space)?,
        )?;
        let (rotation_after, drifted) =
            self.inertia
                .compensated_momentum_step(rotation, first.momentum, self.duration_s)?;
        let second = self.stage(
            drifted,
            drive_at(rotation_after, drive)?,
            bath_at(rotation_after, bath_space)?,
        )?;
        let momentum_after_body = rounded(second.momentum, "viscous rigid-body output momentum")?;
        let accepted = held::finite_vector(momentum_after_body, "viscous accepted momentum")?;
        let angular_velocity_after_body = rounded(
            velocity(self.inertia, accepted)?,
            "viscous output angular velocity",
        )?;
        let kinetic_change = energy(self.inertia, accepted)? - energy(self.inertia, initial)?;
        let drift_defect = energy(self.inertia, drifted)? - energy(self.inertia, first.momentum)?;
        let drive_work = first.drive_work + second.drive_work;
        let bath_work = first.bath_work + second.bath_work;
        let source_work = if bath_space.is_some() {
            drive_work + bath_work
        } else {
            drive_work
        };
        let heat = first.heat + second.heat;
        let drive_impulse = add(
            to_space(rotation, first.drive_impulse)?,
            to_space(rotation_after, second.drive_impulse)?,
        );
        let drag_impulse = add(
            to_space(rotation, first.drag_impulse)?,
            to_space(rotation_after, second.drag_impulse)?,
        );
        let balance = subtract(
            subtract(
                to_space(rotation_after, accepted)?,
                to_space(rotation, initial)?,
            ),
            add(drive_impulse, drag_impulse),
        );
        let motion = ViscousRigidBodyStep {
            rotation_after,
            momentum_after_body,
            angular_velocity_after_body,
            drive_angular_impulse_space: rounded(drive_impulse, "viscous drive impulse")?,
            drag_angular_impulse_space: rounded(drag_impulse, "viscous drag impulse")?,
            angular_momentum_residual_space: rounded(balance, "viscous angular balance")?,
            drive_work_j: scalar(drive_work.to_f64(), "viscous drive work")?,
            dissipated_energy_j: scalar(heat.to_f64(), "viscous heat")?,
            kinetic_energy_change_j: scalar(kinetic_change.to_f64(), "viscous kinetic change")?,
            drift_energy_defect_j: scalar(drift_defect.to_f64(), "viscous free-flow defect")?,
            work_energy_residual_j: scalar(
                (kinetic_change - source_work + heat - drift_defect).to_f64(),
                "viscous work balance",
            )?,
        };
        Ok(ViscousRotatingBathStep {
            motion,
            bath_work_j: scalar(bath_work.to_f64(), "viscous bath work")?,
        })
    }

    fn stage(
        &self,
        momentum: Vector,
        drive: Vector,
        bath: Option<Vector>,
    ) -> Result<Stage, LieStepError> {
        let before = velocity(self.inertia, momentum)?;
        let mut midpoint = before;
        if let Some(factor) = &self.factor {
            let initial_drag = apply(self.damping_body, relative(before, bath))?;
            let rhs = subtract(drive, initial_drag).map(|x| dd(self.half_s) * x);
            let mut change = [0.0; 3];
            for i in 0..3 {
                change[i] = compensated::divide(rhs[i], dd(self.factor_scale))?.to_f64();
            }
            factor.solve(&mut change);
            for i in 0..3 {
                midpoint[i] = compensated::checked(
                    before[i] + dd(0.5) * dd(change[i]),
                    "viscous midpoint velocity",
                )?;
            }
        }
        let drive_impulse = drive.map(|x| dd(self.half_s) * x);
        let trial_drag = apply(self.damping_body, relative(midpoint, bath))?;
        let impulse = subtract(drive, trial_drag).map(|x| dd(self.half_s) * x);
        let mut after = [Dd::ZERO; 3];
        for i in 0..3 {
            after[i] = compensated::checked(momentum[i] + impulse[i], "viscous kicked momentum")?;
        }
        let after_velocity = velocity(self.inertia, after)?;
        let mean = core::array::from_fn(|i| dd(0.5) * before[i] + dd(0.5) * after_velocity[i]);
        let slip = relative(mean, bath);
        let drag = apply(self.damping_body, slip)?;
        let drag_impulse = drag.map(|x| dd(-self.half_s) * x);
        let actual = subtract(after, momentum);
        for i in 0..3 {
            let residual = scalar(
                (actual[i] - drive_impulse[i] - drag_impulse[i]).to_f64(),
                "viscous impulse residual",
            )?;
            let scale = actual[i]
                .to_f64()
                .abs()
                .max(drive_impulse[i].to_f64().abs())
                .max(drag_impulse[i].to_f64().abs());
            if residual.abs() > 4096.0 * f64::EPSILON * scale {
                return Err(invalid("unresolved implicit viscous torque"));
            }
        }
        let drive_work = dd(self.half_s) * dot(mean, drive);
        let bath_work = bath.map_or(Dd::ZERO, |spin| dot(drag_impulse, spin));
        let source_work = if bath.is_some() {
            drive_work + bath_work
        } else {
            drive_work
        };
        let heat = dd(self.half_s) * dot(slip, drag);
        if scalar(heat.to_f64(), "viscous dissipative work")? < 0.0 {
            return Err(invalid("nonnegative viscous dissipative work"));
        }
        let kinetic_change = dot(actual, mean);
        let work_residual = scalar(
            (kinetic_change - source_work + heat).to_f64(),
            "viscous stage work residual",
        )?;
        // Opposing cross-axis work can cancel in the net dot product.
        // Admit rounding against local component transfers, never background
        // kinetic energy, so such cancellation does not create a false refusal.
        let mut scale = 0.0_f64;
        for i in 0..3 {
            for component_impulse in [actual[i], drive_impulse[i], drag_impulse[i]] {
                scale = scale.max(
                    scalar(
                        (component_impulse * mean[i]).to_f64(),
                        "viscous component work",
                    )?
                    .abs(),
                );
            }
        }
        if let Some(spin) = bath {
            for i in 0..3 {
                for rate in [spin[i], slip[i]] {
                    scale = scale.max(
                        scalar(
                            (drag_impulse[i] * rate).to_f64(),
                            "viscous bath/slip component work",
                        )?
                        .abs(),
                    );
                }
            }
        }
        if work_residual.abs() > 12288.0 * f64::EPSILON * scale {
            return Err(invalid("unresolved implicit viscous torque work"));
        }
        Ok(Stage {
            momentum: after,
            drive_impulse,
            drag_impulse,
            drive_work,
            bath_work,
            heat,
        })
    }
}

struct Stage {
    momentum: Vector,
    drive_impulse: Vector,
    drag_impulse: Vector,
    drive_work: Dd,
    bath_work: Dd,
    heat: Dd,
}

fn admit_damping(input: Mat3) -> Result<Mat3, LieStepError> {
    if !input.m.iter().all(|x| x.is_finite()) {
        return Err(invalid("finite viscous damping tensor"));
    }
    let scale = input.max_abs();
    if scale == 0.0 {
        return Ok(input);
    }
    let mut matrix = input.m.map(|x| x / scale);
    for (i, j) in [(1, 3), (2, 6), (5, 7)] {
        if (matrix[i] - matrix[j]).abs() > 64.0 * f64::EPSILON {
            return Err(invalid("symmetric viscous damping tensor"));
        }
        let average = 0.5 * matrix[i] + 0.5 * matrix[j];
        matrix[i] = average;
        matrix[j] = average;
    }
    // A symmetric 3-by-3 form is PSD iff all principal minors are nonnegative.
    // Keep product tails: an eigensolver's off-diagonal stopping tolerance can
    // otherwise overlook an indefinite zero-diagonal form with tiny coupling.
    for i in 0..3 {
        if matrix[3 * i + i] < 0.0 {
            return Err(invalid("positive-semidefinite viscous damping tensor"));
        }
        for j in i + 1..3 {
            if (matrix[3 * i + i] == 0.0 || matrix[3 * j + j] == 0.0) && matrix[3 * i + j] != 0.0 {
                return Err(invalid("positive-semidefinite viscous damping tensor"));
            }
            let minor = dd(matrix[3 * i + i]) * dd(matrix[3 * j + j])
                - dd(matrix[3 * i + j]) * dd(matrix[3 * i + j]);
            if negative(minor) {
                return Err(invalid("positive-semidefinite viscous damping tensor"));
            }
        }
    }
    let [a, b, c, _, d, e, _, _, f] = matrix.map(dd);
    let determinant = a * (d * f - e * e) - b * (b * f - c * e) + c * (b * e - c * d);
    if negative(determinant) {
        return Err(invalid("positive-semidefinite viscous damping tensor"));
    }
    Ok(Mat3 {
        m: matrix.map(|x| x * scale),
    })
}

fn negative(value: Dd) -> bool {
    value.hi < 0.0 || (value.hi == 0.0 && value.lo < 0.0)
}

fn velocity(inertia: InertiaTensor3, momentum: Vector) -> Result<Vector, LieStepError> {
    let axes = inertia.principal_to_body();
    let principal = compensated::spatial_vector_to_body(axes, momentum)?;
    let moments = components(inertia.principal_moments());
    let mut rates = [Dd::ZERO; 3];
    for i in 0..3 {
        rates[i] = compensated::divide(principal[i], dd(moments[i]))?;
    }
    compensated::spatial_vector_to_body(axes.inverse(), rates)
}

fn drive_at(rotation: So3, drive: HeldTorque) -> Result<Vector, LieStepError> {
    match drive {
        HeldTorque::Body(value) => held::finite_vector(value, "viscous body drive torque"),
        HeldTorque::Space(value) => compensated::spatial_vector_to_body(
            rotation,
            held::finite_vector(value, "viscous spatial drive torque")?,
        ),
    }
}

fn bath_at(rotation: So3, bath_space: Option<Vector>) -> Result<Option<Vector>, LieStepError> {
    bath_space
        .map(|spin| compensated::spatial_vector_to_body(rotation, spin))
        .transpose()
}

fn relative(spin: Vector, bath: Option<Vector>) -> Vector {
    bath.map_or(spin, |value| subtract(spin, value))
}

fn to_space(rotation: So3, value: Vector) -> Result<Vector, LieStepError> {
    compensated::spatial_vector_to_body(rotation.inverse(), value)
}

fn apply(matrix: Mat3, value: Vector) -> Result<Vector, LieStepError> {
    let mut result = [Dd::ZERO; 3];
    for i in 0..3 {
        for (j, component) in value.iter().enumerate() {
            result[i] = result[i] + dd(matrix.m[3 * i + j]) * *component;
        }
        result[i] = compensated::checked(result[i], "viscous tensor product")?;
    }
    Ok(result)
}

fn energy(inertia: InertiaTensor3, momentum: Vector) -> Result<Dd, LieStepError> {
    compensated::checked(
        dd(0.5) * dot(momentum, velocity(inertia, momentum)?),
        "viscous rotational kinetic energy",
    )
}

fn dot(a: Vector, b: Vector) -> Dd {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn add(a: Vector, b: Vector) -> Vector {
    core::array::from_fn(|i| a[i] + b[i])
}

fn subtract(a: Vector, b: Vector) -> Vector {
    core::array::from_fn(|i| a[i] - b[i])
}

fn rounded(value: Vector, context: &'static str) -> Result<Vec3, LieStepError> {
    let [x, y, z] = value.map(Dd::to_f64);
    Ok(Vec3::new(
        scalar(x, context)?,
        scalar(y, context)?,
        scalar(z, context)?,
    ))
}

fn components(value: Vec3) -> [f64; 3] {
    [value.x, value.y, value.z]
}

fn dd(value: f64) -> Dd {
    Dd::from_f64(value)
}

fn scalar(value: f64, context: &'static str) -> Result<f64, LieStepError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid(context))
    }
}

fn invalid(context: &'static str) -> LieStepError {
    LieStepError::InvalidParameter { context }
}
