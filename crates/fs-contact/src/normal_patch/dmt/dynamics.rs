//! Finite-mass motion on the original DMT potential, including adhesive snap.
//!
//! For signed approach x and prescribed spring stage D, the physical equation
//! is m x'' = Q + K(D-x) - F_DMT(x) - c x'. The optional dashpot is referenced
//! to the fixed opposing surface, NOT the moving stage. Q is a held applied
//! force. The source's exact reversible work supplies the scalar discrete
//! gradient; the spring, inertia and dashpot use midpoint effort/velocity.
//! Thus stage work and applied-load work close kinetic + contact + spring
//! storage plus c * mean_velocity^2 * dt. Snap energy becomes motion and only
//! the explicitly supplied dashpot dissipates it. No equilibrium jump occurs.
//!
//! This is the scalar discrete-gradient construction of Gonzalez (1996),
//! J. Nonlinear Science 6, 449-467, applied to the existing source potential.
//! Energy closure is not a trajectory-error estimate: explicit step/range
//! limits and temporal refinement remain necessary, especially during snap.

use super::{DmtError, DmtResponse, SourceDmtSphereContact};
use fs_blake3::{ContentHash, DomainHasher};

/// Explicit physical range, work and numerical admission for every proposal.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DmtDynamicsBudget {
    pub maximum_step_s: f64,
    /// At most this many scalar root evaluations, in addition to two bounds.
    pub maximum_iterations: usize,
    /// Finite lower approach bound; the original source owns the upper bound.
    pub minimum_indentation_m: f64,
    pub maximum_indentation_step_m: f64,
    pub maximum_abs_stage_displacement_m: f64,
    pub maximum_abs_velocity_m_s: f64,
    /// Applied, contact, spring and damping efforts must each satisfy this.
    pub maximum_abs_force_n: f64,
    pub force_absolute_tolerance_n: f64,
    pub force_relative_tolerance: f64,
    pub energy_absolute_tolerance_j: f64,
    pub energy_relative_tolerance: f64,
}

/// A refused proposal leaves its immutable input state available for replay.
#[derive(Clone, Debug, PartialEq)]
pub enum DmtDynamicsError {
    Contact(DmtError),
    Invalid {
        field: &'static str,
    },
    StateMismatch,
    /// Inertia/damping/spring do not dominate the negative discrete curvature.
    /// Refine dt; no arbitrary root or quasistatic branch is substituted.
    NonuniqueStep,
    /// The unique root is outside the admitted spatial interval for this step.
    OutsideStepRange,
    Budget {
        field: &'static str,
        value: f64,
        limit: f64,
    },
    NoConvergence {
        residual_n: f64,
        iterations: usize,
    },
    EnergyBalance {
        residual_j: f64,
        tolerance_j: f64,
    },
}
impl From<DmtError> for DmtDynamicsError {
    fn from(value: DmtError) -> Self {
        Self::Contact(value)
    }
}
impl core::fmt::Display for DmtDynamicsError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "DMT inertial contact refusal: {self:?}")
    }
}
impl std::error::Error for DmtDynamicsError {}

/// Accepted normal history. Its fields cannot be changed through diagnostics.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DmtDynamicsState {
    identity: ContentHash,
    indentation_m: f64,
    velocity_m_s: f64,
    stage_displacement_m: f64,
    time_s: f64,
    steps: u64,
}
impl DmtDynamicsState {
    #[must_use]
    pub const fn indentation_m(self) -> f64 {
        self.indentation_m
    }
    #[must_use]
    pub const fn velocity_m_s(self) -> f64 {
        self.velocity_m_s
    }
    #[must_use]
    pub const fn stage_displacement_m(self) -> f64 {
        self.stage_displacement_m
    }
    #[must_use]
    pub const fn time_s(self) -> f64 {
        self.time_s
    }
    #[must_use]
    pub const fn steps(self) -> u64 {
        self.steps
    }
}

/// One fully admitted successor and its independently inspectable work ports.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DmtDynamicsStep {
    state: DmtDynamicsState,
    /// Actual ENDPOINT response from the original source. This supplies real
    /// evolving area/repulsive load to shear; gap attraction supplies no area.
    pub contact: DmtResponse,
    pub kinetic_energy_j: f64,
    pub spring_energy_j: f64,
    pub total_energy_j: f64,
    /// Exact-work mean of the original signed contact force. Its negative is
    /// the generalized force on the mass; it is not an endpoint reaction.
    pub discrete_contact_force_n: f64,
    /// Midpoint spring force ON the mass, K mean(D-x).
    pub spring_force_n: f64,
    /// Fixed-reference dashpot force ON the mass, -c mean(v).
    pub damping_force_n: f64,
    /// Original contact potential change, using its stable exact-piece work.
    pub contact_work_j: f64,
    /// Work supplied by the prescribed spring stage, mean(K(D-x)) * delta D.
    pub prescribed_stage_work_j: f64,
    /// Work supplied by the held applied force Q, Q * delta x.
    pub applied_load_work_j: f64,
    /// Physical dashpot loss only; never a quasistatic snap-release estimate.
    pub damping_dissipation_j: f64,
    /// Stable storage differences + physical loss - both supplied works.
    pub energy_residual_j: f64,
    /// m delta(v)/dt - Q - spring + contact - damping.
    pub force_residual_n: f64,
    pub iterations: usize,
}
impl DmtDynamicsStep {
    /// Accept by retaining this token. Discarding a proposal is a full rollback.
    #[must_use]
    pub const fn state(&self) -> DmtDynamicsState {
        self.state
    }
}

/// One normal relative coordinate with explicit finite (or reduced) mass.
/// All source conditions stay frozen; the contact geometry evolves with x.
/// Construction owns its source once; stepping is allocation-free.
#[derive(Clone, Debug, PartialEq)]
pub struct DmtInertialContact {
    contact: SourceDmtSphereContact,
    mass_kg: f64,
    stiffness_n_per_m: f64,
    damping_n_s_per_m: f64,
    budget: DmtDynamicsBudget,
    identity: ContentHash,
}

#[derive(Clone, Copy)]
struct Trial {
    contact: DmtResponse,
    contact_force: f64,
    contact_work: f64,
    spring_force: f64,
    mean_velocity: f64,
    velocity: f64,
    residual: f64,
    tolerance: f64,
}

impl DmtInertialContact {
    /// K=0 admits direct held-load control; c=0 retains conservative dynamics.
    /// The fixed-reference damper must be explicitly measured or authored.
    pub fn new(
        contact: SourceDmtSphereContact,
        mass_kg: f64,
        stiffness_n_per_m: f64,
        damping_n_s_per_m: f64,
        budget: DmtDynamicsBudget,
    ) -> Result<Self, DmtDynamicsError> {
        positive(mass_kg, "mass")?;
        nonnegative(stiffness_n_per_m, "spring stiffness")?;
        nonnegative(damping_n_s_per_m, "fixed-reference damping")?;
        for (value, field) in [
            (budget.maximum_step_s, "maximum step"),
            (
                budget.maximum_indentation_step_m,
                "maximum indentation step",
            ),
            (
                budget.maximum_abs_stage_displacement_m,
                "maximum stage displacement",
            ),
            (budget.maximum_abs_velocity_m_s, "maximum velocity"),
            (budget.maximum_abs_force_n, "maximum force"),
            (
                budget.force_absolute_tolerance_n,
                "absolute force tolerance",
            ),
        ] {
            positive(value, field)?;
        }
        for (value, field) in [
            (budget.force_relative_tolerance, "relative force tolerance"),
            (
                budget.energy_relative_tolerance,
                "relative energy tolerance",
            ),
        ] {
            nonnegative(value, field)?;
            if value >= 1.0 {
                return Err(DmtDynamicsError::Invalid { field });
            }
        }
        nonnegative(
            budget.energy_absolute_tolerance_j,
            "absolute energy tolerance",
        )?;
        if !(1..=128).contains(&budget.maximum_iterations)
            || !budget.minimum_indentation_m.is_finite()
            || budget.minimum_indentation_m >= contact.envelope().maximum_indentation_m
        {
            return Err(DmtDynamicsError::Invalid {
                field: "iteration/range budget",
            });
        }
        let mut hash = DomainHasher::new("org.frankensim.fs-contact.dmt-inertial.v1");
        hash.update(contact.identity().as_bytes());
        for x in [mass_kg, stiffness_n_per_m, damping_n_s_per_m] {
            hash.update(&x.to_bits().to_le_bytes());
        }
        Ok(Self {
            contact,
            mass_kg,
            stiffness_n_per_m,
            damping_n_s_per_m,
            budget,
            identity: hash.finalize(),
        })
    }
    #[must_use]
    pub const fn contact(&self) -> &SourceDmtSphereContact {
        &self.contact
    }
    #[must_use]
    pub const fn mass_kg(&self) -> f64 {
        self.mass_kg
    }
    #[must_use]
    pub const fn stiffness_n_per_m(&self) -> f64 {
        self.stiffness_n_per_m
    }
    #[must_use]
    pub const fn damping_n_s_per_m(&self) -> f64 {
        self.damping_n_s_per_m
    }
    #[must_use]
    pub const fn budget(&self) -> DmtDynamicsBudget {
        self.budget
    }

    /// Initialize at explicit finite geometry, velocity and stage, with time 0.
    /// No equilibrium, attachment event or initial impulse is inferred.
    pub fn initialize(
        &self,
        indentation_m: f64,
        velocity_m_s: f64,
        stage_displacement_m: f64,
    ) -> Result<DmtDynamicsState, DmtDynamicsError> {
        let state = DmtDynamicsState {
            identity: self.identity,
            indentation_m,
            velocity_m_s,
            stage_displacement_m,
            time_s: 0.0,
            steps: 0,
        };
        self.check_state(state)?;
        self.total_energy_j(state)?;
        Ok(state)
    }

    /// Complete actual endpoint energy; the contact reference may be negative.
    pub fn total_energy_j(&self, state: DmtDynamicsState) -> Result<f64, DmtDynamicsError> {
        self.check_state(state)?;
        let (kinetic, spring) = self.other_storage(state)?;
        finite(
            kinetic
                + spring
                + self
                    .contact
                    .evaluate(state.indentation_m)?
                    .reversible_energy_j,
            "total stored energy",
        )
    }

    /// Propose one step with a prescribed stage endpoint and a held force Q.
    /// D is interpolated linearly for its work port. The scalar solve refuses
    /// unless 4m/dt^2 + 2c/dt + K > S, where S is the original tail stiffness:
    /// F'_DMT >= -S implies d(DG F)/dx >= -S/2. This makes its root unique
    /// across contact, attraction and separation, without selecting equilibria.
    pub fn advance(
        &self,
        before: DmtDynamicsState,
        dt_s: f64,
        stage_displacement_m: f64,
        applied_force_n: f64,
    ) -> Result<DmtDynamicsStep, DmtDynamicsError> {
        self.check_state(before)?;
        positive(dt_s, "time step")?;
        limit(dt_s, self.budget.maximum_step_s, "time step")?;
        limit(
            stage_displacement_m.abs(),
            self.budget.maximum_abs_stage_displacement_m,
            "stage displacement",
        )?;
        limit(
            applied_force_n.abs(),
            self.budget.maximum_abs_force_n,
            "applied force",
        )?;
        let linear_slope = finite(
            2.0 * (self.mass_kg / dt_s) / dt_s
                + self.damping_n_s_per_m / dt_s
                + 0.5 * self.stiffness_n_per_m,
            "inertial slope",
        )?;
        if linear_slope <= 0.5 * self.contact.tail_stiffness_n_per_m {
            return Err(DmtDynamicsError::NonuniqueStep);
        }
        let x0 = before.indentation_m;
        let lower_x =
            (x0 - self.budget.maximum_indentation_step_m).max(self.budget.minimum_indentation_m);
        let upper_x = (x0 + self.budget.maximum_indentation_step_m)
            .min(self.contact.envelope().maximum_indentation_m);
        if !lower_x.is_finite() || !upper_x.is_finite() || lower_x >= upper_x {
            return Err(DmtDynamicsError::Invalid {
                field: "representable step interval",
            });
        }
        // Solve the velocity INCREMENT, rather than dividing a difference of
        // nearly equal absolute coordinates by dt^2 in the momentum equation.
        // The two exact spatial endpoints define its transformed bracket.
        let mut lo = finite(
            2.0 * ((lower_x - x0) / dt_s - before.velocity_m_s),
            "lower velocity increment",
        )?;
        let mut hi = finite(
            2.0 * ((upper_x - x0) / dt_s - before.velocity_m_s),
            "upper velocity increment",
        )?;
        let trial = |dv, x| self.trial(before, dt_s, stage_displacement_m, applied_force_n, dv, x);
        let lower = trial(lo, lower_x)?;
        let upper = trial(hi, upper_x)?;
        for endpoint in [lower, upper] {
            if endpoint.residual.abs() <= endpoint.tolerance {
                return self.finish(
                    before,
                    dt_s,
                    stage_displacement_m,
                    applied_force_n,
                    endpoint,
                    0,
                );
            }
        }
        if lower.residual > 0.0 || upper.residual < 0.0 {
            return Err(DmtDynamicsError::OutsideStepRange);
        }
        let mut dv = if lo < 0.0 && hi > 0.0 {
            0.0
        } else {
            f64::midpoint(lo, hi)
        };
        let mut last_residual = 0.0;
        for iteration in 1..=self.budget.maximum_iterations {
            let x = finite(
                x0 + dt_s * (before.velocity_m_s + 0.5 * dv),
                "endpoint approach",
            )?;
            if x < lower_x || x > upper_x {
                return Err(DmtDynamicsError::OutsideStepRange);
            }
            let t = trial(dv, x)?;
            last_residual = t.residual;
            if t.residual.abs() <= t.tolerance {
                return self.finish(
                    before,
                    dt_s,
                    stage_displacement_m,
                    applied_force_n,
                    t,
                    iteration,
                );
            }
            if t.residual < 0.0 {
                lo = dv;
            } else {
                hi = dv;
            }
            let dx = x - x0;
            // Exact derivative of the secant force, with its local limit when
            // subtraction would lose its slope. Newton is only an accelerator;
            // every proposal stays bracketed and the original residual admits it.
            let contact_slope = if dx.abs()
                > f64::EPSILON.sqrt()
                    * x.abs()
                        .max(x0.abs())
                        .max(self.contact.interaction_range_m())
            {
                (t.contact.normal_force_n - t.contact_force) / dx
            } else {
                0.25 * (t.contact.one_sided_tangents_n_per_m[0]
                    + t.contact.one_sided_tangents_n_per_m[1])
            };
            let slope = (linear_slope + contact_slope) * (0.5 * dt_s);
            let newton = dv - t.residual / slope;
            let guard = 0.05 * (hi - lo);
            let next = if slope > 0.0 && newton > lo + guard && newton < hi - guard {
                newton
            } else {
                f64::midpoint(lo, hi)
            };
            if next == dv || next < lo || next > hi {
                break;
            }
            dv = next;
        }
        Err(DmtDynamicsError::NoConvergence {
            residual_n: last_residual,
            iterations: self.budget.maximum_iterations,
        })
    }

    fn trial(
        &self,
        before: DmtDynamicsState,
        dt: f64,
        stage: f64,
        load: f64,
        dv: f64,
        x: f64,
    ) -> Result<Trial, DmtDynamicsError> {
        let contact = self.contact.evaluate(x)?;
        let dx = x - before.indentation_m;
        let contact_work = self
            .contact
            .reversible_work_between(before.indentation_m, x)?;
        let contact_force = if dx == 0.0 {
            contact.normal_force_n
        } else {
            contact_work / dx
        };
        let mean_velocity = finite(before.velocity_m_s + 0.5 * dv, "mean velocity")?;
        let velocity = finite(before.velocity_m_s + dv, "endpoint velocity")?;
        let inertia = finite((self.mass_kg / dt) * dv, "inertial force")?;
        let spring_force = finite(
            self.stiffness_n_per_m
                * (0.5 * (before.stage_displacement_m - before.indentation_m) + 0.5 * (stage - x)),
            "spring force",
        )?;
        let damping = finite(self.damping_n_s_per_m * mean_velocity, "damping force")?;
        let residual = finite(
            inertia + contact_force - spring_force + damping - load,
            "force residual",
        )?;
        let scale = inertia
            .abs()
            .max(contact_force.abs())
            .max(spring_force.abs())
            .max(damping.abs())
            .max(load.abs());
        let tolerance = finite(
            self.budget.force_absolute_tolerance_n + self.budget.force_relative_tolerance * scale,
            "force tolerance",
        )?;
        Ok(Trial {
            contact,
            contact_force,
            contact_work,
            spring_force,
            mean_velocity,
            velocity,
            residual,
            tolerance,
        })
    }

    fn finish(
        &self,
        before: DmtDynamicsState,
        dt: f64,
        stage: f64,
        load: f64,
        t: Trial,
        iterations: usize,
    ) -> Result<DmtDynamicsStep, DmtDynamicsError> {
        let time_s = finite(before.time_s + dt, "clock")?;
        if time_s <= before.time_s {
            return Err(DmtDynamicsError::Invalid {
                field: "clock resolution",
            });
        }
        let state = DmtDynamicsState {
            identity: self.identity,
            indentation_m: t.contact.indentation_m,
            velocity_m_s: t.velocity,
            stage_displacement_m: stage,
            time_s,
            steps: before
                .steps
                .checked_add(1)
                .ok_or(DmtDynamicsError::Invalid {
                    field: "step count",
                })?,
        };
        self.check_state(state)?;
        let damping_force_n = -self.damping_n_s_per_m * t.mean_velocity;
        // Audit the actually representable endpoint momentum as well as the
        // increment solved above; never hide a lost velocity update.
        let force_residual_n = finite(
            (self.mass_kg / dt) * (t.velocity - before.velocity_m_s) + t.contact_force
                - t.spring_force
                - damping_force_n
                - load,
            "endpoint force residual",
        )?;
        if force_residual_n.abs() > t.tolerance {
            return Err(DmtDynamicsError::NoConvergence {
                residual_n: force_residual_n,
                iterations,
            });
        }
        for (force, name) in [
            (t.contact_force, "discrete contact force"),
            (t.spring_force, "mean spring force"),
            (damping_force_n, "mean damping force"),
        ] {
            limit(force.abs(), self.budget.maximum_abs_force_n, name)?;
        }
        let dx = state.indentation_m - before.indentation_m;
        let stage_step = stage - before.stage_displacement_m;
        let (kinetic_energy_j, spring_energy_j) = self.other_storage(state)?;
        let (old_kinetic, old_spring) = self.other_storage(before)?;
        let old_contact = self
            .contact
            .evaluate(before.indentation_m)?
            .reversible_energy_j;
        let total_energy_j = finite(
            kinetic_energy_j + spring_energy_j + t.contact.reversible_energy_j,
            "total energy",
        )?;
        let kinetic_change = finite(
            0.5 * self.mass_kg
                * (t.velocity - before.velocity_m_s)
                * (t.velocity + before.velocity_m_s),
            "kinetic work",
        )?;
        let spring_change = finite(t.spring_force * (stage_step - dx), "spring work")?;
        let prescribed_stage_work_j = finite(t.spring_force * stage_step, "stage work")?;
        let applied_load_work_j = finite(load * dx, "applied-load work")?;
        let damping_dissipation_j =
            finite((-damping_force_n * t.mean_velocity) * dt, "damping loss")?;
        let energy_residual_j = finite(
            kinetic_change + t.contact_work + spring_change + damping_dissipation_j
                - prescribed_stage_work_j
                - applied_load_work_j,
            "energy residual",
        )?;
        let scale = [
            kinetic_energy_j,
            spring_energy_j,
            t.contact.reversible_energy_j.abs(),
            old_kinetic,
            old_spring,
            old_contact.abs(),
            prescribed_stage_work_j.abs(),
            applied_load_work_j.abs(),
            damping_dissipation_j,
        ]
        .into_iter()
        .fold(0.0_f64, f64::max);
        let tolerance_j = finite(
            self.budget.energy_absolute_tolerance_j + self.budget.energy_relative_tolerance * scale,
            "energy tolerance",
        )?;
        if energy_residual_j.abs() > tolerance_j {
            return Err(DmtDynamicsError::EnergyBalance {
                residual_j: energy_residual_j,
                tolerance_j,
            });
        }
        Ok(DmtDynamicsStep {
            state,
            contact: t.contact,
            kinetic_energy_j,
            spring_energy_j,
            total_energy_j,
            discrete_contact_force_n: t.contact_force,
            spring_force_n: t.spring_force,
            damping_force_n,
            contact_work_j: t.contact_work,
            prescribed_stage_work_j,
            applied_load_work_j,
            damping_dissipation_j,
            energy_residual_j,
            force_residual_n,
            iterations,
        })
    }
    fn check_state(&self, state: DmtDynamicsState) -> Result<(), DmtDynamicsError> {
        if state.identity != self.identity {
            return Err(DmtDynamicsError::StateMismatch);
        }
        if state.indentation_m < self.budget.minimum_indentation_m {
            return Err(DmtDynamicsError::OutsideStepRange);
        }
        let contact = self.contact.evaluate(state.indentation_m)?;
        limit(
            state.velocity_m_s.abs(),
            self.budget.maximum_abs_velocity_m_s,
            "velocity",
        )?;
        limit(
            state.stage_displacement_m.abs(),
            self.budget.maximum_abs_stage_displacement_m,
            "stage displacement",
        )?;
        for (force, name) in [
            (contact.normal_force_n, "endpoint contact force"),
            (
                self.stiffness_n_per_m * (state.stage_displacement_m - state.indentation_m),
                "endpoint spring force",
            ),
            (
                self.damping_n_s_per_m * state.velocity_m_s,
                "endpoint damping force",
            ),
        ] {
            limit(force.abs(), self.budget.maximum_abs_force_n, name)?;
        }
        Ok(())
    }
    fn other_storage(&self, state: DmtDynamicsState) -> Result<(f64, f64), DmtDynamicsError> {
        let extension = state.stage_displacement_m - state.indentation_m;
        let kinetic = finite(
            (0.5 * self.mass_kg * state.velocity_m_s) * state.velocity_m_s,
            "kinetic energy",
        )?;
        let spring = finite(
            (0.5 * self.stiffness_n_per_m * extension) * extension,
            "spring energy",
        )?;
        Ok((kinetic, spring))
    }
}
fn finite(value: f64, field: &'static str) -> Result<f64, DmtDynamicsError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(DmtDynamicsError::Invalid { field })
    }
}
fn positive(value: f64, field: &'static str) -> Result<(), DmtDynamicsError> {
    if value.is_finite() && value > 0.0 {
        Ok(())
    } else {
        Err(DmtDynamicsError::Invalid { field })
    }
}
fn nonnegative(value: f64, field: &'static str) -> Result<(), DmtDynamicsError> {
    if value.is_finite() && value >= 0.0 {
        Ok(())
    } else {
        Err(DmtDynamicsError::Invalid { field })
    }
}
fn limit(value: f64, bound: f64, field: &'static str) -> Result<(), DmtDynamicsError> {
    finite(value, field)?;
    if value <= bound {
        Ok(())
    } else {
        Err(DmtDynamicsError::Budget {
            field,
            value,
            limit: bound,
        })
    }
}
