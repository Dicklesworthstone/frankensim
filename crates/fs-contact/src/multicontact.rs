//! Simultaneous frictionless normal impacts with assembled passivity admission.
//!
//! One nonnegative quadratic solve couples all geometrically active contacts.
//! Rebound targets are declared per contact; a candidate that creates net
//! relative-contact work is refused, never repaired by silently changing
//! restitution. Geometry, inertia assembly and the physical state update belong
//! to the caller. This is an instantaneous event, not a contact time stepper.

use fs_math::dd::Dd;
use fs_solver::orthant::active_set::{
    SpdOrthantConfig, SpdOrthantError, SpdOrthantSolution, solve_spd_orthant,
};
use fs_solver::orthant::{
    GramOrthantConfig, GramOrthantError, GramOrthantSolution, solve_gram_orthant,
};

use crate::frictional_multicontact::{FrictionalContactError, range};

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
    /// Bounded orthant solve. Gradient tolerance has speed units. For redundant
    /// rows, max_iterations bounds Gram-coordinate sweeps instead of SPD solves.
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

/// What the admitted contact mobility determines. A velocity witness never
/// grants physical authority to an individual redundant contact load.
#[derive(Clone, Debug, PartialEq)]
pub enum MultiImpactDeterminacy {
    /// Independent contact rows determine every impulse coordinate.
    UniqueImpulses,
    /// The normal complementarity law determines W_active*p even though p may
    /// have a nullspace gauge. For a complete mechanical Gram matrix this also
    /// determines the body's velocity increment and inactive kinematics;
    /// arbitrary inactive cross rows have no such guarantee. Impulses are one
    /// checked representative, not a local pressure, heat, or moving-support-
    /// work allocation.
    UniqueVelocity(NormalVelocityWitness),
}

/// Numerical admission of a positive-semidefinite normal contact mobility.
/// This is a separated floating-point rank witness, not an interval proof.
#[derive(Clone, Debug, PartialEq)]
pub struct NormalVelocityWitness {
    pub dimension: usize,
    pub rank: usize,
    /// Original W/scale eigen-equation residual.
    pub spectral_residual: f64,
    pub orthogonality_residual: f64,
    /// Maximum entry of (B*B^T-W)/scale for the Gram candidate's factor B.
    /// Every endpoint and work gate still uses the ORIGINAL supplied W.
    pub reconstruction_residual: f64,
    /// Maximum component of pre-impact relative velocity in the admitted
    /// mobility nullspace. This is bounded in speed units before accepting a
    /// work/loss value for an impulse representative.
    pub relative_speed_range_residual: f64,
}

/// The existing numerical owner used by this event, with no invented SPD
/// receipt for a singular contact set. The PSD lane is a bounded coordinate
/// baseline; it does not claim a general conic or Painleve solver.
#[derive(Clone, Debug, PartialEq)]
pub enum MultiImpactSolver {
    PositiveDefinite(SpdOrthantSolution),
    Semidefinite(GramOrthantSolution),
}

impl MultiImpactSolver {
    /// Checked impulse representative in active contact order.
    pub fn point(&self) -> &[f64] {
        match self {
            Self::PositiveDefinite(solution) => &solution.point,
            Self::Semidefinite(solution) => &solution.point,
        }
    }
}

/// One globally admitted normal impulse vector, in original contact order.
#[derive(Clone, Debug, PartialEq)]
pub struct MultiImpactResult {
    pub determinacy: MultiImpactDeterminacy,
    /// Nonnegative normal impulses; geometrically inactive entries are zero.
    pub normal_impulses: Vec<f64>,
    /// Independently reconstructed speeds, including inactive contacts. Each
    /// complete original-matrix response is reduced before rounding to f64.
    pub relative_normal_speeds_after: Vec<f64>,
    /// Nonnegative rebound targets, used only on the active set.
    pub restitution_targets: Vec<f64>,
    /// Active contacts in caller order; indexes the solver's reduced vectors.
    pub active_indices: Vec<usize>,
    /// Sum of midpoint relative speeds paired with their impulses. Products,
    /// cancellation between contact rows and unrounded velocity responses are
    /// retained in double-double arithmetic before this final f64 reduction.
    /// For a physical assembly this equals kinetic-energy change minus support work.
    /// The PSD lane requires pre-impact relative motion to lie in the mobility
    /// range, making this total independent of a nullspace impulse gauge.
    pub relative_work: f64,
    /// Nonnegative lost mechanical energy, excluding numerical positive work.
    /// This is not a thermal heat allocation.
    pub dissipated_energy: f64,
    /// Positive part of relative work, explicitly bounded by the numerical allowance.
    pub passivity_residual: f64,
    /// Numerical candidate receipt; final acceptance uses the original matrix.
    pub solver: MultiImpactSolver,
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
    /// The bounded PSD candidate did not resolve the complementarity problem.
    GramSolver(GramOrthantError),
    /// The admitted velocity response does not suffice to determine source
    /// work/loss. This conservative lane requires pre-impact relative motion
    /// in the mobility range; it does not assign work to an arbitrary gauge.
    UnresolvedWork {
        range_velocity_residual: f64,
        allowed: f64,
    },
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
/// Zero restitution gives the admitted inelastic projection. Independent rows
/// use the existing SPD active-set solver. Redundant PSD rows use the existing
/// Gram-coordinate baseline with checked spectral/rank admission. For two
/// normal-complementarity solutions, monotonicity gives delta_p^T W delta_p
/// <= 0; PSD mobility then fixes W*p, without asserting that p is unique.
/// A redundant result explicitly labels its impulse representative. Indefinite
/// or rank-ambiguous mobility and unresolved original-system residuals refuse.
/// Pre-impact relative motion must lie in the admitted mobility range for the
/// PSD lane, so total relative work cannot depend on an impulse nullspace gauge.
/// Incompatible prescribed support motion conservatively refuses; individual
/// redundant support loads/work remain unassigned even when the total resolves.
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
    let (solver, determinacy) = solve_mobility(
        &hessian,
        &linear,
        |i| contacts[active_indices[i]].relative_normal_speed,
        config.solver,
        &mut cancelled,
    )?;
    let mut normal_impulses = vec![0.0; n];
    for (&index, &impulse) in active_indices.iter().zip(solver.point()) {
        normal_impulses[index] = impulse;
    }
    let mut relative_normal_speeds_after = Vec::with_capacity(n);
    let mut work = Dd::ZERO;
    let mut work_scale = Dd::ZERO;
    for i in 0..n {
        poll(&mut cancelled)?;
        let before = Dd::from_f64(contacts[i].relative_normal_speed);
        let mut speed = before;
        for (j, &impulse) in normal_impulses.iter().enumerate() {
            speed = add_product(speed, delassus[i * n + j], impulse);
        }
        let rounded_speed = finite(speed.to_f64())?;
        finite(speed.lo)?;
        relative_normal_speeds_after.push(rounded_speed);
        if contacts[i].gap <= contacts[i].activation_tolerance {
            let gradient = finite((speed - Dd::from_f64(restitution_targets[i])).to_f64())?;
            let residual = if normal_impulses[i] > 0.0 {
                gradient.abs()
            } else {
                (-gradient).max(0.0)
            }
            .max((-rounded_speed).max(0.0));
            if residual > config.solver.gradient_tolerance {
                return Err(MultiImpactError::UnresolvedContact { index: i, residual });
            }
        }
        // Pair the impulse with the complete original response, before a
        // rounded endpoint or midpoint can erase a small net transfer. Scaling
        // the two endpoints first avoids overflow in their same-sign sum.
        let midpoint = before * Dd::from_f64(0.5) + speed * Dd::from_f64(0.5);
        let term = midpoint * Dd::from_f64(normal_impulses[i]);
        work = work + term;
        work_scale = work_scale + term.abs();
    }
    let relative_work = finite(work.to_f64())?;
    finite(work.lo)?;
    finite(work_scale.lo)?;
    let work_scale = finite(work_scale.to_f64())?;
    let allowed = finite(config.energy_tolerance + 64.0 * f64::EPSILON * work_scale)?;
    if relative_work > allowed {
        return Err(MultiImpactError::EnergyIncrease {
            work: relative_work,
            allowed,
        });
    }
    poll(&mut cancelled)?;
    Ok(MultiImpactResult {
        determinacy,
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

fn solve_mobility(
    matrix: &[f64],
    linear: &[f64],
    relative_speed: impl Fn(usize) -> f64,
    config: SpdOrthantConfig,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(MultiImpactSolver, MultiImpactDeterminacy), MultiImpactError> {
    let original_error = match solve_spd_orthant(matrix, linear, config, &mut *cancelled) {
        Ok(solution) => {
            return Ok((
                MultiImpactSolver::PositiveDefinite(solution),
                MultiImpactDeterminacy::UniqueImpulses,
            ));
        }
        Err(
            error @ (SpdOrthantError::NotPositiveDefinite { .. }
            | SpdOrthantError::Degenerate { .. }),
        ) => error,
        Err(error) => return Err(error.into()),
    };
    // Reuse exactly the same original-matrix PSD, separated-rank, symmetry,
    // eigen-equation and orthogonality admission as simultaneous Coulomb.
    let admission = range::admit(matrix, linear.len(), config, cancelled).map_err(|error| {
        if matches!(error, FrictionalContactError::Cancelled) {
            MultiImpactError::Cancelled
        } else {
            MultiImpactError::Invalid("normal contact mobility admission")
        }
    })?;
    let Some(admission) = admission else {
        return Err(original_error.into());
    };
    let n = linear.len();
    // All normal-complementarity roots share W*p. Their midpoint-work
    // difference is therefore v_before dot delta_p, with delta_p in ker W.
    // Admit only range-compatible relative motion before assigning TOTAL loss.
    // Independently moving redundant supports can violate this even when all
    // candidate impulses are passive and native velocity is unique.
    let velocity_scale = (0..n).map(|i| relative_speed(i).abs()).fold(0.0, f64::max);
    let mut normalized_range_residual = 0.0_f64;
    if velocity_scale > 0.0 {
        for &k in &admission.null {
            poll(cancelled)?;
            let component: f64 = (0..n)
                .map(|i| admission.vectors[i * n + k] * (relative_speed(i) / velocity_scale))
                .sum();
            normalized_range_residual = normalized_range_residual.max(component.abs());
        }
    }
    let relative_speed_range_residual = finite(normalized_range_residual * velocity_scale)?;
    let allowed = finite(
        config.gradient_tolerance + admission.zero_band * (n as f64).sqrt() * velocity_scale,
    )?;
    if relative_speed_range_residual > allowed {
        return Err(MultiImpactError::UnresolvedWork {
            range_velocity_residual: relative_speed_range_residual,
            allowed,
        });
    }
    let root_scale = admission.matrix_scale.sqrt();
    let mut factor = Vec::with_capacity(n);
    for i in 0..n {
        poll(cancelled)?;
        let mut row = Vec::with_capacity(admission.retained.len());
        for &k in &admission.retained {
            row.push(finite(
                admission.vectors[i * n + k] * admission.eigenvalues[k].sqrt() * root_scale,
            )?);
        }
        factor.push(row);
    }
    // The factor supplies a PSD numerical candidate, never a replacement
    // physical mobility. Quantify its reconstruction before solving; the
    // caller below independently enforces the law and work on original W.
    let mut reconstruction_residual = 0.0_f64;
    for i in 0..n {
        poll(cancelled)?;
        for j in 0..n {
            let normalized_entry: f64 = factor[i]
                .iter()
                .zip(&factor[j])
                .map(|(a, b)| (a / root_scale) * (b / root_scale))
                .sum();
            reconstruction_residual = reconstruction_residual
                .max((normalized_entry - matrix[i * n + j] / admission.matrix_scale).abs());
        }
    }
    if !reconstruction_residual.is_finite() || reconstruction_residual > 2.0 * admission.zero_band {
        return Err(original_error.into());
    }
    let solution = solve_gram_orthant(
        &factor,
        linear,
        GramOrthantConfig {
            max_variables: config.max_variables,
            max_columns: 64,
            max_setup_products: 64 * 64 * 64,
            max_sweeps: config.max_iterations,
            gradient_tolerance: config.gradient_tolerance,
        },
        &mut *cancelled,
    )
    .map_err(|error| {
        if matches!(error, GramOrthantError::Cancelled) {
            MultiImpactError::Cancelled
        } else {
            MultiImpactError::GramSolver(error)
        }
    })?;
    Ok((
        MultiImpactSolver::Semidefinite(solution),
        MultiImpactDeterminacy::UniqueVelocity(NormalVelocityWitness {
            dimension: n,
            rank: admission.retained.len(),
            spectral_residual: admission.spectral_residual,
            orthogonality_residual: admission.orthogonality_residual,
            reconstruction_residual,
            relative_speed_range_residual,
        }),
    ))
}

// Preserve the fused owner's representable cancellation when the individual
// product overflows. The common path remains an unscaled error-free product.
pub(crate) fn add_product(sum: Dd, a: f64, b: f64) -> Dd {
    let value = sum + Dd::from_f64(a) * Dd::from_f64(b);
    if value.hi.is_finite() && value.lo.is_finite() {
        value
    } else {
        let (a, b) = if a.abs() >= b.abs() {
            (0.5 * a, b)
        } else {
            (a, 0.5 * b)
        };
        (sum * Dd::from_f64(0.5) + Dd::from_f64(a) * Dd::from_f64(b)) * Dd::from_f64(2.0)
    }
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
    fn full_cross_response_retains_small_motion_between_opposing_large_impulses() {
        // Three unit-mobility contacts stop three independent unit speeds.
        // The fourth, inactive observation has Jacobian [a, 1, -a, a].
        // Its fourth mechanical rate is 1/a, so its initial speed is zero
        // and its actual final speed is one. Neither an impulse iteration
        // order nor the large opposing first/third transfers may erase it.
        let a = 2.0_f64.powi(54);
        let raw = [
            1.0,
            0.0,
            0.0,
            a,
            0.0,
            1.0,
            0.0,
            1.0,
            0.0,
            0.0,
            1.0,
            -a,
            a,
            1.0,
            -a,
            3.0 * a * a,
        ];
        let mut separated = contact(0.0, 0.0);
        separated.gap = 1.0;
        let original = [contact(-1.0, 0.0); 3];
        for order in [[0, 1, 2, 3], [2, 0, 3, 1]] {
            let contacts: Vec<_> = order
                .iter()
                .map(|&i| if i == 3 { separated } else { original[i] })
                .collect();
            let mut matrix = Vec::with_capacity(16);
            for &i in &order {
                for &j in &order {
                    matrix.push(raw[4 * i + j]);
                }
            }
            let result =
                solve_normal_impacts(&contacts, &matrix, MultiImpactConfig::default(), || false)
                    .unwrap();
            for (row, &original_row) in order.iter().enumerate() {
                if original_row == 3 {
                    assert_eq!(result.normal_impulses[row], 0.0);
                    assert_eq!(result.relative_normal_speeds_after[row], 1.0);
                } else {
                    assert_eq!(result.normal_impulses[row], 1.0);
                    assert_eq!(result.relative_normal_speeds_after[row], 0.0);
                }
            }
            assert_eq!(result.relative_work, -1.5);
            assert_eq!(result.dissipated_energy, 1.5);
        }
    }

    #[test]
    fn negative_gaps_and_nonpassive_curvature_refuse_without_repair() {
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
                &[1.0, 2.0, 2.0, 1.0],
                MultiImpactConfig::default(),
                || false,
            ),
            Err(MultiImpactError::Solver(
                SpdOrthantError::NotPositiveDefinite { .. } | SpdOrthantError::Degenerate { .. }
            ))
        ));
    }

    #[test]
    fn simultaneous_counterwork_matches_independent_three_body_kinetic_changes() {
        // Three unit masses with velocities [0, v, v/2] and contact rows
        // [1,-1,0], [0,1,-1]. Near e=sqrt(3)/2, large contact works cancel;
        // a rounded relative midpoint can lose the entire net transfer.
        let v = 2.0_f64.powi(30);
        let e = 0.5 * 3.0_f64.sqrt();
        let mut config = MultiImpactConfig::default();
        config.solver.gradient_tolerance = 1.0e-5;
        for restitution in [e.next_down(), e, e.next_up()] {
            let out = solve_normal_impacts(
                &[contact(-v, restitution), contact(0.5 * v, 0.0)],
                &[2.0, -1.0, -1.0, 2.0],
                config,
                || false,
            )
            .unwrap();
            let p = Dd::from_f64(out.normal_impulses[0]);
            let q = Dd::from_f64(out.normal_impulses[1]);
            // Reconstruct each body's velocity change, then its kinetic
            // change. No contact-space W or returned endpoint enters this
            // independent physical account, and no bulk K is subtracted.
            let impulses = [p, q - p, -q];
            let velocities = [0.0, v, 0.5 * v];
            let mut kinetic_change = Dd::ZERO;
            for (impulse, velocity) in impulses.into_iter().zip(velocities) {
                kinetic_change = kinetic_change
                    + impulse * (Dd::from_f64(velocity) + impulse * Dd::from_f64(0.5));
            }
            let expected = kinetic_change.to_f64();
            assert!(expected.abs() < 1024.0);
            assert!(
                (out.relative_work - expected).abs() < 1.0e-10,
                "contact work {}, native body change {expected}",
                out.relative_work,
            );
            assert_eq!(out.dissipated_energy, (-expected).max(0.0));
        }
    }

    #[test]
    fn redundant_normal_rows_resolve_motion_without_inventing_local_loads() {
        for restitution in [0.0, 0.5, 1.0] {
            let out = solve_normal_impacts(
                &[contact(-1.0, restitution); 2],
                &[1.0; 4],
                MultiImpactConfig::default(),
                || false,
            )
            .unwrap();
            let MultiImpactDeterminacy::UniqueVelocity(witness) = &out.determinacy else {
                panic!("duplicate contact rows cannot claim unique load allocation");
            };
            assert_eq!((witness.dimension, witness.rank), (2, 1));
            assert!(witness.reconstruction_residual < 1e-14);
            assert!(matches!(out.solver, MultiImpactSolver::Semidefinite(_)));
            close(out.normal_impulses.iter().sum(), 1.0 + restitution);
            for speed in out.relative_normal_speeds_after {
                close(speed, restitution);
            }
            close(
                out.dissipated_energy,
                0.5 * (1.0 - restitution * restitution),
            );
            // Both [1+e, 0] and [(1+e)/2, (1+e)/2] satisfy the same law.
            // The exposed representative must not erase that physical gauge.
        }
    }

    #[test]
    fn redundant_four_foot_impact_preserves_native_response_under_permutation_and_units() {
        // Unit mass and roll/pitch inertias. Normal Jacobian at (x,y,0) is
        // [1, y, -x]; four contact rows span only three mechanical coordinates.
        let points = [[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], [-1.0, 1.0]];
        for order in [[0, 1, 2, 3], [2, 0, 3, 1]] {
            for (mass_scale, speed_scale) in [(1.0, 1.0), (4.0, 8.0)] {
                let mut matrix = Vec::new();
                let mut contacts = Vec::new();
                for &i in &order {
                    let [x, y] = points[i];
                    contacts.push(contact((-1.0 + 0.2 * y + 0.15 * x) * speed_scale, 0.5));
                    for &j in &order {
                        matrix.push((1.0 + x * points[j][0] + y * points[j][1]) / mass_scale);
                    }
                }
                let out =
                    solve_normal_impacts(&contacts, &matrix, MultiImpactConfig::default(), || {
                        false
                    })
                    .unwrap();
                let MultiImpactDeterminacy::UniqueVelocity(witness) = &out.determinacy else {
                    panic!("four-foot support has a one-dimensional load gauge");
                };
                assert_eq!(witness.rank, 3);
                let mut increment = [0.0; 3];
                for (k, &i) in order.iter().enumerate() {
                    let [x, y] = points[i];
                    let impulse = out.normal_impulses[k] / (mass_scale * speed_scale);
                    increment[0] += impulse;
                    increment[1] += y * impulse;
                    increment[2] -= x * impulse;
                    assert!(
                        (out.relative_normal_speeds_after[k] / speed_scale
                            + 0.5 * contacts[k].relative_normal_speed / speed_scale)
                            .abs()
                            < 2e-9
                    );
                }
                for (value, expected) in increment.into_iter().zip([1.5, -0.3, 0.225]) {
                    assert!((value - expected).abs() < 2e-9);
                }
                // Independent native kinetic energy, no assumed load split.
                let loss = 0.5 * (1.0 - 0.25) * (1.0 + 0.2 * 0.2 + 0.15 * 0.15);
                assert!(
                    (out.dissipated_energy / (mass_scale * speed_scale * speed_scale) - loss).abs()
                        < 2e-9
                );
            }
        }
    }

    #[test]
    fn redundant_rows_preserve_distinct_restitution_targets_and_global_passivity() {
        let out = solve_normal_impacts(
            &[contact(-1.0, 0.0), contact(-1.0, 1.0)],
            &[1.0; 4],
            MultiImpactConfig::default(),
            || false,
        )
        .unwrap();
        assert_eq!(out.restitution_targets, vec![0.0, 1.0]);
        close(out.normal_impulses[0], 0.0);
        close(out.normal_impulses[1], 2.0);
        close(out.dissipated_energy, 0.0);

        // Duplicate the first row of the known energy-injecting chain. A
        // redundant solve must retain the original global passivity refusal.
        let error = solve_normal_impacts(
            &[contact(-1.0, 1.0), contact(0.5, 1.0), contact(-1.0, 1.0)],
            &[2.5, -1.5, 2.5, -1.5, 2.5, -1.5, 2.5, -1.5, 2.5],
            MultiImpactConfig::default(),
            || false,
        )
        .unwrap_err();
        let MultiImpactError::EnergyIncrease { work, .. } = error else {
            panic!("wrong refusal: {error:?}");
        };
        close(work, 0.109375);
    }

    #[test]
    fn redundant_prescribed_motion_cannot_assign_arbitrary_dissipation() {
        // W*p fixes p1+p2=2 and endpoint speeds [1,0]. Every nonnegative
        // split meets the declared Newton laws, but work=-p2 ranges from -2
        // to 0 J. A unique velocity must not publish one arbitrary loss.
        let error = solve_normal_impacts(
            &[contact(-1.0, 1.0), contact(-2.0, 0.0)],
            &[1.0; 4],
            MultiImpactConfig::default(),
            || false,
        )
        .unwrap_err();
        let MultiImpactError::UnresolvedWork {
            range_velocity_residual,
            allowed,
        } = error
        else {
            panic!("wrong refusal: {error:?}");
        };
        assert!(range_velocity_residual > 0.7);
        assert!(range_velocity_residual > allowed);
    }

    #[test]
    fn weak_positive_modes_are_not_erased_to_obtain_a_redundant_solution() {
        let off_diagonal = 1.0 - 1e-13;
        assert!(matches!(
            solve_normal_impacts(
                &[contact(-1.0, 0.0); 2],
                &[1.0, off_diagonal, off_diagonal, 1.0],
                MultiImpactConfig::default(),
                || false,
            ),
            Err(MultiImpactError::Solver(SpdOrthantError::Degenerate { .. }))
        ));
    }

    #[test]
    fn redundant_solver_cancellation_never_publishes_a_partial_event() {
        let mut polls = 0;
        solve_normal_impacts(
            &[contact(-1.0, 0.5); 2],
            &[1.0; 4],
            MultiImpactConfig::default(),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
        for stop in 1..=polls {
            let mut seen = 0;
            assert_eq!(
                solve_normal_impacts(
                    &[contact(-1.0, 0.5); 2],
                    &[1.0; 4],
                    MultiImpactConfig::default(),
                    || {
                        seen += 1;
                        seen == stop
                    },
                ),
                Err(MultiImpactError::Cancelled),
                "checkpoint {stop}",
            );
        }
    }

    #[test]
    fn maximum_redundant_contact_set_still_uses_one_physical_velocity() {
        let out = solve_normal_impacts(
            &[contact(-1.0, 0.5); 64],
            &[1.0; 64 * 64],
            MultiImpactConfig::default(),
            || false,
        )
        .unwrap();
        let MultiImpactDeterminacy::UniqueVelocity(witness) = out.determinacy else {
            panic!("duplicated rows cannot determine individual loads");
        };
        assert_eq!((witness.dimension, witness.rank), (64, 1));
        close(out.normal_impulses.iter().sum(), 1.5);
        close(out.dissipated_energy, 0.375);
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
