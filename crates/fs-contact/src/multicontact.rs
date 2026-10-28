//! Simultaneous frictionless normal impacts with assembled passivity admission.
//!
//! One nonnegative quadratic solve couples all geometrically active contacts.
//! Rebound targets are declared per contact; a candidate that creates net
//! relative-contact work is refused, never repaired by silently changing
//! restitution. Geometry, inertia assembly and the physical state update belong
//! to the caller. This is an instantaneous event, not a contact time stepper.

use fs_solver::orthant::active_set::{
    SpdOrthantConfig, SpdOrthantError, SpdOrthantSolution, solve_spd_orthant,
};

/// Geometry and relative motion for one normal contact, in consistent units.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NormalContactImpact {
    /// Nonnegative signed gap; penetration always refuses.
    pub gap: f64,
    /// Nonnegative proximity threshold admitting this contact to the solve.
    pub activation_tolerance: f64,
    /// Relative normal speed; negative means approach.
    pub relative_normal_speed: f64,
    /// Declared Newton restitution in [0, 1].
    pub restitution: f64,
}

/// Numerical admission for the coupled impact event.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MultiImpactConfig {
    /// Bounded positive-definite orthant solve. Gradient tolerance has speed units.
    pub solver: SpdOrthantConfig,
    /// Nonnegative absolute allowance for numerical positive relative work.
    /// In SI this is joules; any allowed positive work remains in the receipt.
    pub energy_tolerance: f64,
}

impl Default for MultiImpactConfig {
    fn default() -> Self {
        Self {
            solver: SpdOrthantConfig {
                max_variables: 64,
                max_iterations: 1024,
                gradient_tolerance: 1.0e-10,
                symmetry_tolerance: 1.0e-12,
                relative_pivot_tolerance: 1.0e-12,
            },
            energy_tolerance: 1.0e-10,
        }
    }
}

/// One globally admitted normal impulse vector, in original contact order.
#[derive(Clone, Debug, PartialEq)]
pub struct MultiImpactResult {
    /// Nonnegative normal impulses; geometrically inactive entries are zero.
    pub normal_impulses: Vec<f64>,
    /// Independently reconstructed speeds, including inactive contacts.
    pub relative_normal_speeds_after: Vec<f64>,
    /// Nonnegative rebound targets, used only on the active set.
    pub restitution_targets: Vec<f64>,
    /// Active contacts in caller order; indexes the solver's reduced vectors.
    pub active_indices: Vec<usize>,
    /// Sum of midpoint relative speeds paired with their impulses.
    /// For a physical assembly this equals kinetic-energy change minus support work.
    pub relative_work: f64,
    /// Nonnegative lost mechanical energy, excluding numerical positive work.
    /// This is not a thermal heat allocation.
    pub dissipated_energy: f64,
    /// Positive part of relative work, explicitly bounded by the numerical allowance.
    pub passivity_residual: f64,
    /// Numerical solve receipt for the active principal Delassus matrix.
    pub solver: SpdOrthantSolution,
}

/// A refusal returns no partially admitted impact state.
#[derive(Clone, Debug, PartialEq)]
pub enum MultiImpactError {
    /// Invalid shape, cap, or numerical energy allowance.
    Invalid(&'static str),
    /// Contact geometry or restitution is inadmissible.
    InvalidContact {
        /// Original contact index.
        index: usize,
        /// Failed physical input.
        reason: &'static str,
    },
    /// Generic positive-definite solve refused, including rank/curvature ambiguity.
    Solver(SpdOrthantError),
    /// Finite inputs generated unrepresentable arithmetic.
    NonFinite,
    /// Reconstructed endpoint separation or complementarity did not resolve.
    UnresolvedContact {
        /// Original contact index.
        index: usize,
        /// Speed-unit residual.
        residual: f64,
    },
    /// Declared rebound targets would inject energy into the assembled system.
    EnergyIncrease {
        /// Candidate positive relative work.
        work: f64,
        /// Explicit numerical allowance, including arithmetic roundoff.
        allowed: f64,
    },
    /// Caller cancellation before publication.
    Cancelled,
}

impl core::fmt::Display for MultiImpactError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "simultaneous normal impact: {self:?}")
    }
}
impl core::error::Error for MultiImpactError {}
impl From<SpdOrthantError> for MultiImpactError {
    fn from(error: SpdOrthantError) -> Self {
        if matches!(error, SpdOrthantError::Cancelled) {
            Self::Cancelled
        } else {
            Self::Solver(error)
        }
    }
}

/// Resolve simultaneous frictionless normal impacts from a supplied Delassus matrix.
///
/// The row-major n*n matrix is J M^-1 J^T, with n equal to the contact count.
/// Geometry-active contacts, including initially separating ones, share one
/// solve. Their target is -e*min(v_before,0): approaching contacts may rebound;
/// other active contacts must finish nonclosing. The convex solve imposes
/// nonnegative impulse, nonnegative endpoint speed minus target, and normal
/// complementarity.
///
/// Restitution targets can be mutually incompatible with global passivity.
/// The actual relative impulse work is therefore checked after the solve.
/// No sequential pairwise impulses or hidden restitution reduction is used.
/// Zero restitution gives the admitted inelastic projection. Singular or
/// insufficiently resolved active matrices refuse through the generic solver;
/// this API does not choose an arbitrary impulse from a nonunique set.
///
/// Positions are never advanced or repaired. A caller performing a physical
/// update must check its measured endpoint speeds, world momentum and energy
/// against this result. No friction, CCD, persistent contact or adjoint is claimed.
///
/// # Errors
/// Returns an error for inadmissible input, unresolved solve, positive
/// relative work beyond the allowance, or cancellation.
pub fn solve_normal_impacts(
    contacts: &[NormalContactImpact],
    delassus: &[f64],
    config: MultiImpactConfig,
    mut cancelled: impl FnMut() -> bool,
) -> Result<MultiImpactResult, MultiImpactError> {
    poll(&mut cancelled)?;
    let n = contacts.len();
    if n > 64 || n > config.solver.max_variables || n.checked_mul(n) != Some(delassus.len()) {
        return Err(MultiImpactError::Invalid("contact count or Delassus shape"));
    }
    if !config.energy_tolerance.is_finite() || config.energy_tolerance < 0.0 {
        return Err(MultiImpactError::Invalid("energy tolerance"));
    }
    if delassus.iter().any(|value| !value.is_finite()) {
        return Err(MultiImpactError::NonFinite);
    }
    let mut active_indices = Vec::with_capacity(n);
    let mut restitution_targets = vec![0.0; n];
    for (index, contact) in contacts.iter().enumerate() {
        poll(&mut cancelled)?;
        let reason = if !contact.gap.is_finite() || contact.gap < 0.0 {
            Some("finite nonnegative gap required; penetration is not repaired")
        } else if !contact.activation_tolerance.is_finite() || contact.activation_tolerance < 0.0 {
            Some("finite nonnegative activation tolerance required")
        } else if !contact.relative_normal_speed.is_finite() {
            Some("relative normal speed must be finite")
        } else if !contact.restitution.is_finite() || !(0.0..=1.0).contains(&contact.restitution) {
            Some("restitution must lie in [0,1]")
        } else {
            None
        };
        if let Some(reason) = reason {
            return Err(MultiImpactError::InvalidContact { index, reason });
        }
        restitution_targets[index] = -contact.restitution * contact.relative_normal_speed.min(0.0);
        if contact.gap <= contact.activation_tolerance {
            active_indices.push(index);
        }
    }
    let m = active_indices.len();
    let mut hessian = Vec::with_capacity(m * m);
    let mut linear = Vec::with_capacity(m);
    for &i in &active_indices {
        poll(&mut cancelled)?;
        linear.push(finite(
            contacts[i].relative_normal_speed - restitution_targets[i],
        )?);
        for &j in &active_indices {
            hessian.push(delassus[i * n + j]);
        }
    }
    let solver = solve_spd_orthant(&hessian, &linear, config.solver, &mut cancelled)?;
    let mut normal_impulses = vec![0.0; n];
    for (&index, &impulse) in active_indices.iter().zip(&solver.point) {
        normal_impulses[index] = impulse;
    }
    let mut relative_normal_speeds_after = Vec::with_capacity(n);
    for i in 0..n {
        poll(&mut cancelled)?;
        let mut speed = contacts[i].relative_normal_speed;
        for (j, &impulse) in normal_impulses.iter().enumerate() {
            speed = finite(delassus[i * n + j].mul_add(impulse, speed))?;
        }
        relative_normal_speeds_after.push(speed);
    }
    for &index in &active_indices {
        let speed = relative_normal_speeds_after[index];
        let gradient = finite(speed - restitution_targets[index])?;
        let residual = if normal_impulses[index] > 0.0 {
            gradient.abs()
        } else {
            (-gradient).max(0.0)
        }
        .max((-speed).max(0.0));
        if residual > config.solver.gradient_tolerance {
            return Err(MultiImpactError::UnresolvedContact { index, residual });
        }
    }
    let mut relative_work = 0.0;
    let mut work_scale = 0.0;
    for i in 0..n {
        poll(&mut cancelled)?;
        let midpoint = f64::midpoint(
            contacts[i].relative_normal_speed,
            relative_normal_speeds_after[i],
        );
        let work = finite(midpoint * normal_impulses[i])?;
        relative_work = finite(relative_work + work)?;
        work_scale = finite(work_scale + work.abs())?;
    }
    let allowed = finite(config.energy_tolerance + 64.0 * f64::EPSILON * work_scale)?;
    if relative_work > allowed {
        return Err(MultiImpactError::EnergyIncrease {
            work: relative_work,
            allowed,
        });
    }
    poll(&mut cancelled)?;
    Ok(MultiImpactResult {
        normal_impulses,
        relative_normal_speeds_after,
        restitution_targets,
        active_indices,
        relative_work,
        dissipated_energy: (-relative_work).max(0.0),
        passivity_residual: relative_work.max(0.0),
        solver,
    })
}

fn finite(value: f64) -> Result<f64, MultiImpactError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(MultiImpactError::NonFinite)
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), MultiImpactError> {
    if cancelled() {
        Err(MultiImpactError::Cancelled)
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn contact(speed: f64, restitution: f64) -> NormalContactImpact {
        NormalContactImpact {
            gap: 0.0,
            activation_tolerance: 0.0,
            relative_normal_speed: speed,
            restitution,
        }
    }
    fn close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() < 2.0e-10,
            "{actual} != {expected}"
        );
    }

    #[test]
    fn simultaneous_restitution_matches_independent_coupled_solution() {
        let w = [2.5, -1.5, -1.5, 2.5];
        let out = solve_normal_impacts(
            &[contact(-1.2, 0.5), contact(-0.8, 0.5)],
            &w,
            MultiImpactConfig::default(),
            || false,
        )
        .unwrap();
        for (got, want) in out.normal_impulses.iter().zip([1.575, 1.425]) {
            close(*got, want);
        }
        close(out.relative_normal_speeds_after[0], 0.6);
        close(out.relative_normal_speeds_after[1], 0.4);
        close(out.dissipated_energy, 0.7575);
        close(out.passivity_residual, 0.0);
    }

    #[test]
    fn active_separating_contact_can_couple_but_passive_rebound_must_not_gain_energy() {
        let w = [2.5, -1.5, -1.5, 2.5];
        let contacts = [contact(-1.0, 1.0), contact(0.5, 1.0)];
        let err = solve_normal_impacts(&contacts, &w, MultiImpactConfig::default(), || false)
            .unwrap_err();
        match err {
            MultiImpactError::EnergyIncrease { work, .. } => close(work, 0.109375),
            other => panic!("wrong refusal: {other:?}"),
        }
        let out = solve_normal_impacts(
            &[contact(-1.0, 0.0), contact(0.5, 0.0)],
            &w,
            MultiImpactConfig::default(),
            || false,
        )
        .unwrap();
        close(out.normal_impulses[0], 0.4375);
        close(out.normal_impulses[1], 0.0625);
        close(out.relative_normal_speeds_after[0], 0.0);
        close(out.relative_normal_speeds_after[1], 0.0);
        close(out.dissipated_energy, 0.203125);
    }

    #[test]
    fn inactive_contacts_keep_zero_impulse_and_full_cross_response() {
        let mut inactive = contact(-100.0, 1.0);
        inactive.gap = 0.1;
        let out = solve_normal_impacts(
            &[inactive, contact(-1.0, 0.0)],
            &[2.5, -1.5, -1.5, 2.5],
            MultiImpactConfig::default(),
            || false,
        )
        .unwrap();
        assert_eq!(out.active_indices, vec![1]);
        assert_eq!(out.normal_impulses[0], 0.0);
        close(out.normal_impulses[1], 0.4);
        close(out.relative_normal_speeds_after[0], -100.6);
        close(out.dissipated_energy, 0.2);
    }

    #[test]
    fn negative_gaps_and_nonunique_curvature_refuse_without_repair() {
        let mut penetrating = contact(-1.0, 0.0);
        penetrating.gap = -1.0e-15;
        penetrating.activation_tolerance = 1.0e-3;
        assert!(matches!(
            solve_normal_impacts(&[penetrating], &[1.0], MultiImpactConfig::default(), || {
                false
            },),
            Err(MultiImpactError::InvalidContact { index: 0, .. })
        ));
        assert!(matches!(
            solve_normal_impacts(
                &[contact(-1.0, 0.0); 2],
                &[1.0; 4],
                MultiImpactConfig::default(),
                || false,
            ),
            Err(MultiImpactError::Solver(
                SpdOrthantError::NotPositiveDefinite { .. } | SpdOrthantError::Degenerate { .. }
            ))
        ));
    }

    #[test]
    fn empty_noop_and_cancelled_entry_are_explicit() {
        let out = solve_normal_impacts(&[], &[], MultiImpactConfig::default(), || false).unwrap();
        assert!(out.normal_impulses.is_empty());
        assert_eq!(out.relative_work, 0.0);
        assert_eq!(
            solve_normal_impacts(
                &[contact(-1.0, 0.5)],
                &[1.0],
                MultiImpactConfig::default(),
                || true,
            ),
            Err(MultiImpactError::Cancelled)
        );
    }
}
