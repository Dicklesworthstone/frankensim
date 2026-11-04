//! Coupled planar Coulomb impacts on the original frame inertia.

use super::impulse::{ImpulseContinuation, checked_action};
use super::{Frame2d, FrameImpulse, Rayleigh, SolidError};
use fs_contact::frictional_multicontact::{
    FrictionalContact, FrictionalContactConfig, FrictionalContactError, FrictionalContactResiduals,
    FrictionalContactResult, NormalVelocityLaw, solve_frictional_contacts,
};
use fs_math::dd::Dd;

const ROUNDOFF: f64 = 512.0 * f64::EPSILON;

/// One declared planar contact at the accepted frame configuration.
///
/// Both covectors use original node-major coordinates. They measure velocity
/// at the SAME physical contact, in an orthonormal normal/tangent frame:
/// relative velocity = covector dot frame velocity - obstacle velocity.
/// Translational weights are dimensionless and rotational weights have metres.
/// Multiplying these SAME rows by impulses [N s] supplies conjugate nodal
/// impulses and moments. Internal contacts include both signed traces.
///
/// The caller owns geometry, physical unit directions and event discovery.
/// This interface does not infer an independent friction coefficient, rebuild
/// the contact point, or fit force weights separately from velocity weights.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FrameFrictionalContact<'a> {
    /// Stable physical contact identity; duplicates refuse. Reductions use ID
    /// order, while returned vectors retain the caller's original order.
    pub id: u64,
    pub normal_covector: &'a [f64],
    pub tangent_covector: &'a [f64],
    /// Nonnegative signed gap [m]; penetration is never repaired.
    pub gap: f64,
    /// Nonnegative contact activation distance [m].
    pub activation_tolerance: f64,
    /// Prescribed [normal, tangent] obstacle speeds [m/s], held fixed across
    /// the event. The same orientation is used by both contact covectors.
    pub obstacle_velocity: [f64; 2],
    pub friction_coefficient: f64,
    pub normal_law: NormalVelocityLaw,
}

/// Accepted actual frame response and the original coupled Coulomb receipt.
#[derive(Clone, Debug, PartialEq)]
pub struct FrameFrictionalImpact {
    /// The shared owner uses [normal, tangent, auxiliary-zero] coordinates.
    /// The unused second tangent has EXACTLY zero mobility, pre/post velocity
    /// and accepted impulse. It supplies neither inertia nor compliance.
    ///
    /// This planar embedding conservatively uses the PSD lane; retain its
    /// determinacy and total-work witnesses. A checked representative does
    /// not grant unique local loads or individual source-work allocations.
    pub contact: FrictionalContactResult,
    pub impulse: FrameImpulse,
    /// Actual stored [normal, tangent] relative velocities [m/s], including
    /// geometrically inactive contacts, in original input order.
    pub relative_velocities_after: Vec<[f64; 2]>,
    /// Coulomb residuals measured using the actual stored endpoint velocities.
    pub residuals: FrictionalContactResiduals,
    /// Work supplied by prescribed obstacle motion [J]. Stationary frame
    /// supports have zero work and their reaction is in the native impulse.
    pub obstacle_work_j: f64,
    /// Actual midpoint relative velocity paired with the complete impulse [J].
    /// Native kinetic change = obstacle work + relative work.
    pub relative_work_j: f64,
    /// Negative part of TOTAL relative impact work [J], not frictional heat.
    pub dissipated_energy_j: f64,
    /// Positive part of relative work [J], bounded by the numerical allowance.
    pub passivity_residual_j: f64,
}

/// Every refusal leaves all accepted frame state unchanged.
#[derive(Debug)]
pub enum FrameFrictionalImpactError {
    Frame(SolidError),
    Contact(FrictionalContactError),
    UnresolvedEndpoint {
        index: usize,
        /// 0 is normal, 1 is tangent.
        component: usize,
        residual: f64,
        allowed: f64,
    },
    UnresolvedWork {
        residual: f64,
        allowed: f64,
    },
}
impl core::fmt::Display for FrameFrictionalImpactError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "frame frictional impact: {self:?}")
    }
}
impl core::error::Error for FrameFrictionalImpactError {}
impl From<SolidError> for FrameFrictionalImpactError {
    fn from(error: SolidError) -> Self {
        Self::Frame(error)
    }
}
impl From<FrictionalContactError> for FrameFrictionalImpactError {
    fn from(error: FrictionalContactError) -> Self {
        Self::Contact(error)
    }
}
fn invalid(what: &str) -> FrameFrictionalImpactError {
    SolidError::InvalidInput { what: what.into() }.into()
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), FrameFrictionalImpactError> {
    if cancelled() {
        Err(FrictionalContactError::Cancelled.into())
    } else {
        Ok(())
    }
}
fn finite(value: Dd) -> Result<f64, FrameFrictionalImpactError> {
    let rounded = value.to_f64();
    if value.hi.is_finite() && value.lo.is_finite() && rounded.is_finite() {
        Ok(rounded)
    } else {
        Err(invalid(
            "frame frictional response or work is unrepresentable",
        ))
    }
}
fn dot(a: &[f64], b: &[f64]) -> Dd {
    a.iter().zip(b).fold(Dd::ZERO, |sum, (&a, &b)| {
        sum + Dd::from_f64(a) * Dd::from_f64(b)
    })
}

/// Recheck the shared normal/complementarity and maximum-dissipation law on
/// measured native velocities. This is an endpoint audit, not another solve.
fn law_residuals(
    impulse: [f64; 3],
    after: [Dd; 2],
    target: f64,
    mu: f64,
    config: FrictionalContactConfig,
) -> Result<FrictionalContactResiduals, FrameFrictionalImpactError> {
    let pn = Dd::from_f64(impulse[0]);
    let pt = Dd::from_f64(impulse[1]);
    let gradient = after[0] - Dd::from_f64(target);
    let radius = Dd::from_f64(mu) * pn;
    finite(radius)?;
    let tangential_speed = finite(after[1])?;
    let residuals = FrictionalContactResiduals {
        normal_impulse_violation: (-impulse[0]).max(0.0),
        normal_velocity_violation: (-finite(gradient)?).max(0.0),
        normal_complementarity: finite((pn * gradient).abs())?,
        cone_violation: finite(pt.abs() - radius)?.max(0.0),
        sliding_impulse_residual: if tangential_speed.abs() > config.velocity_tolerance {
            finite((pt + radius * Dd::from_f64(tangential_speed.signum())).abs())?
        } else {
            0.0
        },
        tangential_duality_gap: finite((radius * after[1].abs() + pt * after[1]).abs())?,
    };
    if residuals.normal_impulse_violation > config.impulse_tolerance
        || residuals.normal_velocity_violation > config.velocity_tolerance
        || residuals.normal_complementarity > config.complementarity_tolerance
        || residuals.cone_violation > config.impulse_tolerance
        || residuals.sliding_impulse_residual > config.impulse_tolerance
        || residuals.tangential_duality_gap > config.complementarity_tolerance
    {
        return Err(FrictionalContactError::ResidualNotSatisfied(residuals).into());
    }
    Ok(residuals)
}

impl Frame2d {
    /// Resolve simultaneous planar Signorini-Coulomb impacts without advancing
    /// time. The non-associated normal and tangential laws share ONE solve:
    /// friction is not clipped after an independently selected normal kick.
    ///
    /// The original mass supplies the full coupled contact mobility. Only the
    /// two physical response columns per contact are solved; the shared 3-row
    /// owner's unused second tangent is exactly zero. The accepted auxiliary
    /// impulse and velocity must also be exactly zero. No artificial inertia
    /// regularizes this PSD embedding, and its native determinacy qualification
    /// is retained without inferring unique planar contact loads.
    ///
    /// One free-mass factor serves streamed columns, the combined native
    /// impulse, and finite right-limit acceleration. Additional adapter scratch
    /// is O(ndof + contacts²). Contact count, nonlinear evaluations, iterations,
    /// PSD/rank, work and determinacy admission remain owned by fs-contact.
    /// Cancellation is polled between columns, throughout the shared solve and
    /// immediately before the single publication of accepted frame rates.
    ///
    /// post_event_load, damping and continuation_dt have the
    /// [Self::apply_impulse] meaning. In particular dt is not an impact duration.
    /// Positions, section forces, constitutive history and material clocks
    /// remain fixed. The original positive-free-inertia, stationary-support,
    /// untimed-material and reference-inertia boundaries apply.
    ///
    /// Actual stored endpoint speeds, normal complementarity, Coulomb
    /// sticking/sliding, relative passivity and obstacle work must all pass
    /// before publication. This event does not supply collision detection,
    /// sustained contact, finite-rotation world inertia, or a heat allocation.
    pub fn apply_frictional_impacts(
        &mut self,
        contacts: &[FrameFrictionalContact<'_>],
        post_event_load: &[f64],
        damping: Rayleigh,
        continuation_dt: f64,
        config: FrictionalContactConfig,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FrameFrictionalImpact, FrameFrictionalImpactError> {
        poll(&mut cancelled)?;
        let n = self.ndof();
        let m = contacts.len();
        if m > 21 || m > config.maximum_contacts {
            return Err(invalid(
                "frame Coulomb contact count exceeds the native limit",
            ));
        }
        let mut order: Vec<_> = (0..m).collect();
        order.sort_unstable_by_key(|&i| contacts[i].id);
        let mut before = Vec::with_capacity(m);
        let mut law = Vec::with_capacity(m);
        for contact in contacts {
            poll(&mut cancelled)?;
            let rows = [contact.normal_covector, contact.tangent_covector];
            if rows
                .iter()
                .any(|row| row.len() != n || row.iter().any(|v| !v.is_finite()))
                || contact.obstacle_velocity.iter().any(|v| !v.is_finite())
            {
                return Err(invalid(
                    "frame Coulomb contact needs finite normal/tangent gap derivatives and obstacle velocities",
                ));
            }
            let speeds = [
                dot(rows[0], &self.vel) - Dd::from_f64(contact.obstacle_velocity[0]),
                dot(rows[1], &self.vel) - Dd::from_f64(contact.obstacle_velocity[1]),
            ];
            before.push(speeds);
            law.push(FrictionalContact {
                id: contact.id,
                gap: contact.gap,
                activation_tolerance: contact.activation_tolerance,
                relative_velocity: [finite(speeds[0])?, finite(speeds[1])?, 0.0],
                friction_coefficient: contact.friction_coefficient,
                normal_law: contact.normal_law,
            });
        }
        let mass = self.impulse_mass(post_event_load)?;
        let inverse = mass.prepare_solve(&self.fixed)?;
        let dimension = 3 * m;
        let mut mobility = vec![0.0; dimension * dimension];
        for (j, contact) in contacts.iter().enumerate() {
            for (b, column) in [contact.normal_covector, contact.tangent_covector]
                .into_iter()
                .enumerate()
            {
                poll(&mut cancelled)?;
                let response = inverse.apply(column);
                checked_action(
                    &mass,
                    &response,
                    column,
                    &self.fixed,
                    "frame Coulomb mobility does not resolve original free inertia",
                )?;
                for (i, contact) in contacts.iter().enumerate() {
                    for (a, row) in [contact.normal_covector, contact.tangent_covector]
                        .into_iter()
                        .enumerate()
                    {
                        mobility[(3 * i + a) * dimension + 3 * j + b] =
                            finite(dot(row, &response))?;
                    }
                }
            }
        }
        let contact = solve_frictional_contacts(&law, &mobility, config, &mut cancelled)?;
        if contact.impulses.iter().any(|p| p[2] != 0.0)
            || contact
                .relative_velocities_after
                .iter()
                .any(|v| v[2] != 0.0)
        {
            return Err(invalid(
                "planar Coulomb solve produced an unrepresented auxiliary tangent response",
            ));
        }
        let mut applied = vec![0.0; n];
        for (j, value) in applied.iter_mut().enumerate() {
            let mut sum = Dd::ZERO;
            for &i in &order {
                let rows = [contacts[i].normal_covector, contacts[i].tangent_covector];
                for a in 0..2 {
                    sum = sum + Dd::from_f64(rows[a][j]) * Dd::from_f64(contact.impulses[i][a]);
                }
            }
            *value = finite(sum)?;
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
        let mut velocities_after = vec![[0.0; 2]; m];
        let mut residuals = FrictionalContactResiduals::default();
        let mut relative_work = Dd::ZERO;
        let mut relative_allowance = Dd::from_f64(config.energy_tolerance);
        let mut obstacle_work = Dd::ZERO;
        let mut work_allowance = Dd::ZERO;
        for &i in &order {
            poll(&mut cancelled)?;
            let row = &contacts[i];
            let rows = [row.normal_covector, row.tangent_covector];
            let after = [
                dot(rows[0], &trial.velocity) - Dd::from_f64(row.obstacle_velocity[0]),
                dot(rows[1], &trial.velocity) - Dd::from_f64(row.obstacle_velocity[1]),
            ];
            for a in 0..2 {
                let speed = finite(after[a])?;
                velocities_after[i][a] = speed;
                // Relative/jump scales keep a large common world velocity
                // from granting an unrelated endpoint allowance.
                let mut allowed = config.velocity_tolerance
                    + ROUNDOFF * law[i].relative_velocity[a].abs()
                    + ROUNDOFF * contact.relative_velocities_after[i][a].abs();
                for (&weight, &change) in rows[a].iter().zip(&trial.result.velocity_change) {
                    allowed += ROUNDOFF * (weight * change).abs();
                }
                if !allowed.is_finite() {
                    return Err(invalid(
                        "frame Coulomb endpoint allowance is unrepresentable",
                    ));
                }
                let residual = (speed - contact.relative_velocities_after[i][a]).abs();
                if residual > allowed {
                    return Err(FrameFrictionalImpactError::UnresolvedEndpoint {
                        index: i,
                        component: a,
                        residual,
                        allowed,
                    });
                }
                let impulse = Dd::from_f64(contact.impulses[i][a]);
                let midpoint = before[i][a] * Dd::from_f64(0.5) + after[a] * Dd::from_f64(0.5);
                let relative = midpoint * impulse;
                let supplied = Dd::from_f64(row.obstacle_velocity[a]) * impulse;
                relative_work = relative_work + relative;
                obstacle_work = obstacle_work + supplied;
                relative_allowance = relative_allowance + relative.abs() * Dd::from_f64(ROUNDOFF);
                work_allowance =
                    work_allowance + (relative + supplied).abs() * Dd::from_f64(ROUNDOFF);
            }
            if row.gap <= row.activation_tolerance {
                let r = law_residuals(
                    contact.impulses[i],
                    after,
                    contact.normal_targets[i],
                    row.friction_coefficient,
                    config,
                )?;
                residuals.normal_impulse_violation = residuals
                    .normal_impulse_violation
                    .max(r.normal_impulse_violation);
                residuals.normal_velocity_violation = residuals
                    .normal_velocity_violation
                    .max(r.normal_velocity_violation);
                residuals.normal_complementarity = residuals
                    .normal_complementarity
                    .max(r.normal_complementarity);
                residuals.cone_violation = residuals.cone_violation.max(r.cone_violation);
                residuals.sliding_impulse_residual = residuals
                    .sliding_impulse_residual
                    .max(r.sliding_impulse_residual);
                residuals.tangential_duality_gap = residuals
                    .tangential_duality_gap
                    .max(r.tangential_duality_gap);
            }
        }
        let relative_work_j = finite(relative_work)?;
        let obstacle_work_j = finite(obstacle_work)?;
        let relative_allowance = finite(
            relative_allowance + Dd::from_f64(contact.relative_work.abs()) * Dd::from_f64(ROUNDOFF),
        )?;
        let relative_residual =
            finite((relative_work - Dd::from_f64(contact.relative_work)).abs())?;
        if relative_residual > relative_allowance {
            return Err(FrameFrictionalImpactError::UnresolvedWork {
                residual: relative_residual,
                allowed: relative_allowance,
            });
        }
        if relative_work_j > relative_allowance {
            return Err(FrictionalContactError::EnergyIncrease {
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
            return Err(FrameFrictionalImpactError::UnresolvedWork {
                residual: work_residual,
                allowed: work_allowance,
            });
        }
        poll(&mut cancelled)?;
        let impulse = self.accept_impulse(trial);
        Ok(FrameFrictionalImpact {
            contact,
            impulse,
            relative_velocities_after: velocities_after,
            residuals,
            obstacle_work_j,
            relative_work_j,
            dissipated_energy_j: (-relative_work_j).max(0.0),
            passivity_residual_j: relative_work_j.max(0.0),
        })
    }
}
