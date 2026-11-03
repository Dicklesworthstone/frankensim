//! Simultaneous normal impacts through the original frame inertia.
//!
//! Contact geometry is declared at one fixed configuration. The shared
//! fs-contact owner resolves complementarity; the frame owns the actual
//! momentum update and finite-force acceleration before publishing either.

use super::impulse::{ImpulseContinuation, checked_action};
use super::{Frame2d, FrameImpulse, Rayleigh, SolidError};
use fs_contact::multicontact::{
    MultiImpactConfig, MultiImpactError, MultiImpactResult, NormalContactImpact,
    solve_normal_impacts,
};
use fs_math::dd::Dd;

const ROUNDOFF: f64 = 512.0 * f64::EPSILON;

/// One declared contact at the frame's accepted configuration.
///
/// The original node-major covector is the derivative of physical signed gap:
/// gap_rate = normal_covector dot frame_velocity - obstacle_normal_speed.
/// Translational weights are dimensionless and rotational weights have metres.
/// Multiplying the SAME covector by a positive normal impulse [N s] supplies
/// generalized nodal impulses and moments, preserving virtual work.
///
/// Node contacts use the unit world normal at that node. Linear-reference
/// member contacts may use [Frame2d::reference_transverse_trace]. Contacts
/// between two parts of this frame include both signed traces in one row.
/// The caller owns geometry/event discovery and must supply the gap derivative
/// of that geometry, not independently fitted force and velocity weights.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FrameNormalContact<'a> {
    pub normal_covector: &'a [f64],
    /// Nonnegative gap [m]; penetration is refused, never repaired.
    pub gap: f64,
    /// Nonnegative distance [m] admitting the declared contact set.
    pub activation_tolerance: f64,
    /// Prescribed obstacle speed [m/s], with the same normal orientation.
    /// This is held fixed across the instantaneous event.
    pub obstacle_normal_speed: f64,
    /// Newton restitution in [0,1], subject to GLOBAL relative passivity.
    pub restitution: f64,
}

/// One accepted contact event. All work is measured on actual stored frame
/// velocities, in addition to the shared contact owner's algebraic receipt.
#[derive(Clone, Debug, PartialEq)]
pub struct FrameNormalImpact {
    /// Retains the shared law's determinacy qualification. Redundant contacts
    /// may determine velocity and TOTAL relative work without determining each
    /// impulse or its local obstacle-work allocation.
    pub contact: MultiImpactResult,
    pub impulse: FrameImpulse,
    /// Actual represented endpoint gap rates [m/s], in original contact order,
    /// including geometrically inactive contacts.
    pub relative_normal_speeds_after: Vec<f64>,
    /// Work supplied by prescribed obstacles, sum(obstacle_speed * impulse) [J].
    /// Stationary frame supports do no work; their impulses are in impulse.
    pub obstacle_work_j: f64,
    /// Sum of actual midpoint gap rates paired with normal impulses [J].
    /// Native kinetic change = obstacle work + this relative work.
    pub relative_work_j: f64,
    /// Negative part of relative work [J], without any thermal allocation.
    pub dissipated_energy_j: f64,
    /// Positive part of relative work [J], bounded by the numerical allowance.
    pub passivity_residual_j: f64,
}

/// A refusal publishes neither frame rates nor a new Rayleigh reference.
#[derive(Clone, Debug, PartialEq)]
pub enum FrameNormalImpactError {
    Frame(SolidError),
    Contact(MultiImpactError),
    /// Actual stored velocities did not reproduce the shared contact solution.
    UnresolvedEndpoint {
        index: usize,
        residual: f64,
        allowed: f64,
    },
    /// Original momentum work and measured contact work did not close.
    UnresolvedWork {
        residual: f64,
        allowed: f64,
    },
}
impl core::fmt::Display for FrameNormalImpactError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "frame normal impact: {self:?}")
    }
}
impl core::error::Error for FrameNormalImpactError {}
impl From<SolidError> for FrameNormalImpactError {
    fn from(error: SolidError) -> Self {
        Self::Frame(error)
    }
}
impl From<MultiImpactError> for FrameNormalImpactError {
    fn from(error: MultiImpactError) -> Self {
        Self::Contact(error)
    }
}

fn invalid(what: &str) -> FrameNormalImpactError {
    SolidError::InvalidInput { what: what.into() }.into()
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), FrameNormalImpactError> {
    if cancelled() {
        Err(MultiImpactError::Cancelled.into())
    } else {
        Ok(())
    }
}
fn finite(value: Dd) -> Result<f64, FrameNormalImpactError> {
    let rounded = value.to_f64();
    if value.hi.is_finite() && value.lo.is_finite() && rounded.is_finite() {
        Ok(rounded)
    } else {
        Err(invalid("frame contact response or work is unrepresentable"))
    }
}
fn dot(a: &[f64], b: &[f64]) -> Dd {
    a.iter().zip(b).fold(Dd::ZERO, |sum, (&a, &b)| {
        sum + Dd::from_f64(a) * Dd::from_f64(b)
    })
}

impl Frame2d {
    /// Resolve one simultaneous frictionless impact without advancing time.
    ///
    /// The original free mass supplies the FULL coupled mobility
    /// G M_ff^-1 G^T; contacts are never applied as sequential pairwise kicks.
    /// One factor is reused for streamed contact columns, the native velocity
    /// jump, and finite right-limit acceleration. Extra adapter storage is
    /// O(ndof + contacts²), not a retained inverse or ndof-by-contacts matrix.
    /// Contact count/iterations obey config; cancellation is polled between
    /// columns, inside the contact solve and immediately before publication.
    ///
    /// post_event_load, damping, and continuation_dt have exactly the
    /// [Self::apply_impulse] meaning. In particular dt is not an impact
    /// duration. Fixed supports stay stationary. The existing positive-free-
    /// inertia, untimed-material and fixed-reference-inertia boundaries apply.
    ///
    /// All geometrically active contacts participate, including initially
    /// separating contacts another impulse can drive closed. Declared Newton
    /// restitution may be incompatible with global passivity and then refuses;
    /// restitution is never silently reduced. Redundant contact loads retain
    /// the shared solver's checked-representative qualification.
    ///
    /// Actual stored endpoint speeds, complementarity, relative passivity and
    /// moving-obstacle work must pass BEFORE rates or a damping reference
    /// commit. Displacement, member forces, material history and clocks stay
    /// fixed. This is an event response, not collision detection, persistent
    /// contact, friction, a moving-support integrator or a heat partition.
    pub fn apply_normal_impacts(
        &mut self,
        contacts: &[FrameNormalContact<'_>],
        post_event_load: &[f64],
        damping: Rayleigh,
        continuation_dt: f64,
        config: MultiImpactConfig,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FrameNormalImpact, FrameNormalImpactError> {
        poll(&mut cancelled)?;
        let n = self.ndof();
        let m = contacts.len();
        if m > 64
            || m > config.solver.max_variables
            || !config.solver.gradient_tolerance.is_finite()
            || config.solver.gradient_tolerance < 0.0
            || !config.energy_tolerance.is_finite()
            || config.energy_tolerance < 0.0
        {
            return Err(invalid("frame contact count or numerical tolerances"));
        }
        let mut before = Vec::with_capacity(m);
        let mut law = Vec::with_capacity(m);
        for contact in contacts {
            poll(&mut cancelled)?;
            if contact.normal_covector.len() != n
                || contact.normal_covector.iter().any(|v| !v.is_finite())
                || !contact.obstacle_normal_speed.is_finite()
            {
                return Err(invalid(
                    "frame contact needs one finite gap derivative per DOF and finite obstacle speed",
                ));
            }
            let speed = dot(contact.normal_covector, &self.vel)
                - Dd::from_f64(contact.obstacle_normal_speed);
            before.push(speed);
            law.push(NormalContactImpact {
                gap: contact.gap,
                activation_tolerance: contact.activation_tolerance,
                relative_normal_speed: finite(speed)?,
                restitution: contact.restitution,
            });
        }
        let mass = self.impulse_mass(post_event_load)?;
        let inverse = mass.prepare_solve(&self.fixed)?;
        let mut mobility = vec![0.0; m * m];
        for (j, contact) in contacts.iter().enumerate() {
            poll(&mut cancelled)?;
            let response = inverse.apply(contact.normal_covector);
            checked_action(
                &mass,
                &response,
                contact.normal_covector,
                &self.fixed,
                "frame contact mobility does not resolve original free inertia",
            )?;
            for (i, row) in contacts.iter().enumerate() {
                mobility[i * m + j] = finite(dot(row.normal_covector, &response))?;
            }
        }
        let contact = solve_normal_impacts(&law, &mobility, config, &mut cancelled)?;
        let mut applied = vec![0.0; n];
        for (j, value) in applied.iter_mut().enumerate() {
            *value = finite(contacts.iter().zip(&contact.normal_impulses).fold(
                Dd::ZERO,
                |sum, (row, &impulse)| {
                    sum + Dd::from_f64(row.normal_covector[j]) * Dd::from_f64(impulse)
                },
            ))?;
        }
        poll(&mut cancelled)?;
        let trial = self.trial_impulse(
            &applied,
            ImpulseContinuation {
                load: post_event_load,
                damping,
                dt: continuation_dt,
            },
            &mass,
            &inverse,
        )?;
        let mut speeds_after = Vec::with_capacity(m);
        let mut relative_work = Dd::ZERO;
        let mut relative_allowance = Dd::from_f64(config.energy_tolerance);
        let mut obstacle_work = Dd::ZERO;
        let mut work_allowance = Dd::ZERO;
        for (i, row) in contacts.iter().enumerate() {
            poll(&mut cancelled)?;
            let after =
                dot(row.normal_covector, &trial.velocity) - Dd::from_f64(row.obstacle_normal_speed);
            let speed = finite(after)?;
            speeds_after.push(speed);
            // Scale endpoint rounding by relative motion and the actual jump,
            // never by an unrelated large common world velocity.
            let mut allowed = config.solver.gradient_tolerance
                + ROUNDOFF * law[i].relative_normal_speed.abs()
                + ROUNDOFF * contact.relative_normal_speeds_after[i].abs();
            for (&weight, &change) in row
                .normal_covector
                .iter()
                .zip(&trial.result.velocity_change)
            {
                allowed += ROUNDOFF * (weight * change).abs();
            }
            if !allowed.is_finite() {
                return Err(invalid(
                    "frame contact endpoint allowance is unrepresentable",
                ));
            }
            let mut residual = (speed - contact.relative_normal_speeds_after[i]).abs();
            if row.gap <= row.activation_tolerance {
                let gradient = finite(after - Dd::from_f64(contact.restitution_targets[i]))?;
                residual = residual
                    .max(if contact.normal_impulses[i] > 0.0 {
                        gradient.abs()
                    } else {
                        (-gradient).max(0.0)
                    })
                    .max((-speed).max(0.0));
            }
            if residual > allowed {
                return Err(FrameNormalImpactError::UnresolvedEndpoint {
                    index: i,
                    residual,
                    allowed,
                });
            }
            let impulse = Dd::from_f64(contact.normal_impulses[i]);
            let midpoint = before[i] * Dd::from_f64(0.5) + after * Dd::from_f64(0.5);
            let relative = midpoint * impulse;
            let supplied = Dd::from_f64(row.obstacle_normal_speed) * impulse;
            relative_work = relative_work + relative;
            obstacle_work = obstacle_work + supplied;
            relative_allowance = relative_allowance + relative.abs() * Dd::from_f64(ROUNDOFF);
            work_allowance = work_allowance + (relative + supplied).abs() * Dd::from_f64(ROUNDOFF);
        }
        let relative_work_j = finite(relative_work)?;
        let obstacle_work_j = finite(obstacle_work)?;
        let relative_allowance = finite(
            relative_allowance + Dd::from_f64(contact.relative_work.abs()) * Dd::from_f64(ROUNDOFF),
        )?;
        let relative_residual =
            finite((relative_work - Dd::from_f64(contact.relative_work)).abs())?;
        if relative_residual > relative_allowance {
            return Err(FrameNormalImpactError::UnresolvedWork {
                residual: relative_residual,
                allowed: relative_allowance,
            });
        }
        if relative_work_j > relative_allowance {
            return Err(MultiImpactError::EnergyIncrease {
                work: relative_work_j,
                allowed: relative_allowance,
            }
            .into());
        }
        for (j, &impulse) in applied.iter().enumerate() {
            let midpoint = Dd::from_f64(self.vel[j]) * Dd::from_f64(0.5)
                + Dd::from_f64(trial.velocity[j]) * Dd::from_f64(0.5);
            work_allowance =
                work_allowance + (midpoint * Dd::from_f64(impulse)).abs() * Dd::from_f64(ROUNDOFF);
        }
        let work_allowance = finite(
            work_allowance
                + Dd::from_f64(trial.result.kinetic_energy_change_j.abs()) * Dd::from_f64(ROUNDOFF),
        )?;
        let work_residual = finite(
            (Dd::from_f64(trial.result.kinetic_energy_change_j) - relative_work - obstacle_work)
                .abs(),
        )?;
        if work_residual > work_allowance {
            return Err(FrameNormalImpactError::UnresolvedWork {
                residual: work_residual,
                allowed: work_allowance,
            });
        }
        poll(&mut cancelled)?;
        let impulse = self.accept_impulse(trial);
        Ok(FrameNormalImpact {
            contact,
            impulse,
            relative_normal_speeds_after: speeds_after,
            obstacle_work_j,
            relative_work_j,
            dissipated_energy_j: (-relative_work_j).max(0.0),
            passivity_residual_j: relative_work_j.max(0.0),
        })
    }
}
