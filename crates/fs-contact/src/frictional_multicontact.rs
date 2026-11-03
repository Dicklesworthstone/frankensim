//! Simultaneous Signorini--Coulomb impulses with the complete Delassus matrix.
//!
//! The normal law is `0 <= p_n _|_ (v_n - target) >= 0`; tangential impulses
//! maximize dissipation in the disk `|p_t| <= mu*p_n` at the FINAL velocity.
//! This is the non-associated Coulomb law, not an associated cone quadratic
//! (which would introduce normal dilation), nor successive normal/friction clips.
//! All rows participate in one semismooth Newton--Krylov natural-map solve,
//! using exact branch Jacobian actions and cached local linearizations. Existing
//! numerical owners provide matrix admission and nonlinear globalization.
//!
//! Independent contact rows use a weighted contraction witness or a positive
//! candidate energy form for impulse determinacy. A separate PSD range witness
//! admits redundant rows only when
//! the coupled law determines velocity despite a possible contact-impulse gauge.
//! A separate range check makes total impulse work independent of that gauge.
//! It does not manufacture per-contact load authority from that gauge. Without
//! a determinacy witness, distinct checked roots report nonuniqueness and an
//! otherwise checked candidate remains explicitly unresolved. General singular
//! manifolds, Painleve existence decisions, poses, CCD, gap repair, persistent
//! time stepping, and classical adjoints are not supplied by this event solve.
//! Total midpoint loss is an impact dissipation receipt, NOT a friction-heat
//! allocation: prescribed motion and material energy partition belong downstream.

use fs_solver::orthant::active_set::{SpdOrthantConfig, SpdOrthantError, solve_spd_orthant};
use fs_solver::{
    NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport, SolverRunProgress,
};
use std::cell::{Cell, RefCell};

mod locked_tangents;
mod natural;
mod normal_gauge;
pub(crate) mod range;

/// Directional impulse ellipses through the original simultaneous Coulomb owner.
pub mod elliptic;

/// What the numerical witness establishes about this accepted response.
#[derive(Clone, Debug, PartialEq)]
pub enum FrictionalContactDeterminacy {
    /// All contact impulse coordinates are determined by the admitted law,
    /// through block contraction or the result's full-rank energy witness.
    UniqueImpulses,
    /// The active contact mobility image `W_active*p` is determined, while
    /// impulses can differ by a mobility-null vector. For an actual complete
    /// mechanical Gram matrix this determines the body's velocity increment
    /// and inactive kinematics; arbitrary inactive cross rows have no such
    /// guarantee. The returned
    /// impulses are one law-satisfying representative, not authoritative local
    /// loads or individual heat/support-work allocations. The result's separate
    /// `impulse_work_witness` establishes TOTAL relative work numerically.
    UniqueVelocity(PsdVelocityWitness),
}

/// Which global sufficient argument establishes the mobility image. Neither
/// criterion assigns physical authority to an individual redundant load.
#[derive(Clone, Debug, PartialEq)]
pub enum PsdVelocityCriterion {
    /// Friction-weighted normal impulses are observable from the mobility
    /// range; their worst-case feedback through tangential velocity is < 1.
    NormalRangeObservable,
    /// A checked solution supplies norm subgradients. Their bilinear normal /
    /// tangential feedback annihilates the mobility nullspace, and the total
    /// energy form is positive on its range. This compares the candidate to
    /// EVERY law-satisfying root, not just its local Newton linearization.
    CandidateFrictionDirection {
        /// `mu_i * s_i`, where `s_i` is a subgradient of tangential speed at the
        /// accepted candidate. In canonical active contact order.
        weighted_directions: Vec<[f64; 2]>,
        /// Maximum entry of `Q_null^T N^T C T Q_range`, normalized by the
        /// largest weighted-direction norm (zero when all directions vanish).
        nullspace_coupling_residual: f64,
        /// Positive lower numerical margin of `I + sym(A)`, with
        /// `A = Lambda^-1/2 Q^T N^T C T Q Lambda^1/2` on the positive range.
        energy_coercivity_lower_bound: f64,
        energy_spectral_residual: f64,
        energy_orthogonality_residual: f64,
        /// Original impulse mismatch `max |p_t + p_n*mu*s|` [N s].
        traction_residual: f64,
        /// `max |mu*|v_t| - (mu*s) dot v_t|` [m/s]. Near-zero velocity uses
        /// the same declared numerical sticking tolerance as the contact law;
        /// at each contact it is bounded by `2*mu*velocity_tolerance` plus
        /// roundoff. Sliding direction is checked relative to mu, even if tiny.
        subgradient_velocity_residual: f64,
    },
}

/// Sufficient numerical determinacy witness on a positive mobility range.
/// When rank equals dimension, a positive candidate energy form determines
/// impulses as well as velocities. A rank-deficient witness determines only
/// the mobility image. Admission is floating point, not an interval proof.
#[derive(Clone, Debug, PartialEq)]
pub struct PsdVelocityWitness {
    pub dimension: usize,
    pub rank: usize,
    pub criterion: PsdVelocityCriterion,
    /// Upper bound kappa in `E <= kappa*E` for two law-satisfying roots,
    /// strictly below one. NormalRangeObservable uses a product of operator
    /// norm bounds; CandidateFrictionDirection uses the negative part of the
    /// minimum symmetric feedback eigenvalue, including numerical margins.
    pub friction_feedback_bound: f64,
    /// Maximum row sum of `diag(mu) * N * Q_null` in the admitted basis. It
    /// must vanish only for NormalRangeObservable; the direction criterion
    /// can prove unique motion while individual normal loads remain free.
    pub weighted_normal_nullspace_residual: f64,
    /// Maximum entry of `W/scale * Q - Q*Lambda`.
    pub spectral_residual: f64,
    /// Maximum entry of `Q^T Q - I`.
    pub orthogonality_residual: f64,
    /// Smallest retained eigenvalue of `W/scale` (positive and separated from
    /// the numerical zero band). `scale` is the maximum absolute row sum of W.
    pub minimum_retained_eigenvalue: f64,
    pub maximum_null_eigenvalue_magnitude: f64,
}

/// Numerical total-work determinacy for a velocity-only contact response.
///
/// Two admissible impulses differ by `z` with `W*z=0`. Their midpoint-work
/// difference is `v_before dot z`; it vanishes when the preimpact velocity is
/// in the mobility range. Exactly frictionless tangent coordinates can be
/// omitted because every law-satisfying impulse is zero there. Structurally
/// locked sliding tangent pairs may instead use the explicit Coulomb-work
/// substitutions below before the same original-mobility range check. This
/// remains a floating admission, not an interval theorem or local heat partition.
#[derive(Clone, Debug, PartialEq)]
pub struct PsdImpulseWorkWitness {
    /// Coordinates of the checked original principal mobility, indexed into
    /// normal-first triples in `active_indices` order. Any omitted coordinate
    /// is a tangent at a contact whose friction coefficient is exactly zero.
    pub impulse_coordinates: Vec<usize>,
    /// Numerically separated rank of that principal mobility.
    pub rank: usize,
    /// Exact-law reductions applied to the incoming work covector before its
    /// range check. Empty means the original preimpact velocity is tested.
    /// These do not modify physical W or identify individual redundant loads.
    pub locked_tangent_substitutions: Vec<LockedTangentWorkSubstitution>,
    /// Work-specific admission for exact locked-slip substitutions. Its bound
    /// applies to EVERY same-motion cone representative, not just the returned
    /// impulse. Empty uses the original numerical velocity-range witness.
    pub locked_tangent_work_admission: Option<LockedTangentWorkAdmission>,
    /// Euclidean residual of the checked work-covector range admission [m/s].
    pub range_velocity_residual: f64,
    /// Allowance of the selected range proof [m/s]. Locked-slip work uses its
    /// energy allowance divided by the admitted impulse-difference bound;
    /// the original velocity-range lane uses its declared velocity tolerance.
    pub allowed_range_velocity_residual: f64,
    /// Original-matrix eigen-equation residual after mobility normalization.
    pub spectral_residual: f64,
    /// Orthogonality residual of the checked range/nullspace basis.
    pub orthogonality_residual: f64,
}

/// An inspectable total-work identity for one exactly locked sliding pair.
///
/// Both original active tangent rows AND columns are exactly zero, so every
/// law root has the same nonzero tangent velocity and p_t = -p_n*mu*v_t/|v_t|.
/// Its work covector [v_n,v_t] can therefore be replaced by
/// [v_n-mu*|v_t|,0,0] when comparing impulse representatives. A zero-slip disk
/// is never given a unique traction direction by this reduction.
#[derive(Clone, Debug, PartialEq)]
pub struct LockedTangentWorkSubstitution {
    /// Contact index in the result's canonical `active_indices` order.
    pub active_contact: usize,
    /// mu*v_t/|v_t| in the original pair of tangent coordinates.
    pub weighted_slip_direction: [f64; 2],
    /// -mu*|v_t|, added to its normal work-covector coordinate [m/s].
    pub normal_velocity_shift: f64,
    /// Retained DD low component of that shift, before adding it to v_n.
    pub normal_velocity_shift_low: f64,
}

/// Numerical work uncertainty for all admitted same-motion Coulomb roots.
///
/// The original-matrix row-space residual r and a positive cone-dual witness
/// give |delta work| <= 2*B*|r|, where B bounds EVERY feasible impulse norm.
/// This is a bounded floating sufficient admission, not exact or interval
/// uniqueness. No bounded cone witness means refusal unless r is exactly zero.
#[derive(Clone, Debug, PartialEq)]
pub struct LockedTangentWorkAdmission {
    /// Scalar active tangent coordinates with mu=0, whose work coefficients
    /// are exactly irrelevant because every admitted impulse is zero there.
    pub zero_friction_coordinates: Vec<usize>,
    /// Upper numerical bound on every admissible impulse's Euclidean norm
    /// [N s]; None means an exactly zero row-space residual needed no bound.
    pub impulse_norm_bound: Option<f64>,
    /// Upper numerical bound on work differences between representatives [J].
    pub maximum_work_difference: f64,
    /// The original selected event's energy allowance, including roundoff [J].
    pub allowed_work_difference: f64,
}

/// Explicit choice of normal endpoint law; these alternatives are not aliases.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum NormalVelocityLaw {
    /// Target `-coefficient*min(v_n_before, 0)`, with coefficient in [0, 1].
    NewtonRestitution { coefficient: f64 },
    /// Caller-owned nonnegative normal endpoint target [m/s]. This does not
    /// authorize penetration repair or an undeclared configuration step.
    Target { speed: f64 },
}

/// One candidate point contact in a caller-owned orthonormal, normal-first frame.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FrictionalContact {
    pub id: u64,
    pub gap: f64,
    pub activation_tolerance: f64,
    /// Relative velocity [normal, tangent 1, tangent 2], m/s. Positive normal
    /// speed separates; positive impulse increases the corresponding speed.
    pub relative_velocity: [f64; 3],
    pub friction_coefficient: f64,
    pub normal_law: NormalVelocityLaw,
}

/// Bounded work and separate numerical/physical acceptance controls.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FrictionalContactConfig {
    /// At most 21 contacts (63 scalar rows), bounding both numerical owners.
    pub maximum_contacts: usize,
    /// Maximum residual, analytic Jacobian and right-preconditioner calls,
    /// shared across all branch probes. Actions use quadratic work; the PSD
    /// lane caches one at-most-63-row dense SVD per Newton point for its
    /// numerical right inverse, with cancellation before/after factorization.
    /// Energy determinacy additionally uses at most two bounded range-energy
    /// eigendecompositions and one 63-by-21 representative SVD per solve;
    /// these fixed-size admission calls are outside the Newton callback count.
    /// Total-work admission reuses the cached range; only a failed full-range
    /// check with exactly frictionless tangents may inspect one smaller
    /// principal mobility using the same bounded eigendecomposition owner.
    pub maximum_evaluations: usize,
    pub maximum_newton_iterations: usize,
    pub velocity_tolerance: f64,
    pub impulse_tolerance: f64,
    /// Normal complementarity and tangential duality-gap tolerance [J].
    pub complementarity_tolerance: f64,
    pub energy_tolerance: f64,
    pub matrix: SpdOrthantConfig,
    pub newton: NewtonKrylovConfig,
}
impl Default for FrictionalContactConfig {
    fn default() -> Self {
        Self {
            maximum_contacts: 21,
            maximum_evaluations: 32768,
            maximum_newton_iterations: 96,
            velocity_tolerance: 1e-9,
            impulse_tolerance: 1e-9,
            complementarity_tolerance: 1e-9,
            energy_tolerance: 1e-9,
            matrix: SpdOrthantConfig::default(),
            newton: NewtonKrylovConfig {
                relative_tolerance: 1e-12,
                // Avoid oversolving early/degenerate linearizations. The
                // analytic action has no difference-scale noise; physical
                // gates and the outer absolute target remain independent.
                forcing_minimum: 1e-6,
                ..Default::default()
            },
        }
    }
}

/// Reconstructed law residuals; no natural-map residual substitutes for these.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct FrictionalContactResiduals {
    pub normal_impulse_violation: f64,
    pub normal_velocity_violation: f64,
    pub normal_complementarity: f64,
    pub cone_violation: f64,
    pub sliding_impulse_residual: f64,
    pub tangential_duality_gap: f64,
}

/// Complete accepted candidate, in the ORIGINAL caller contact order.
#[derive(Clone, Debug, PartialEq)]
pub struct FrictionalContactResult {
    pub determinacy: FrictionalContactDeterminacy,
    pub impulses: Vec<[f64; 3]>,
    pub relative_velocities_after: Vec<[f64; 3]>,
    pub normal_targets: Vec<f64>,
    /// Original indices, ordered by stable contact ID, including separating
    /// geometrically active contacts that another contact may drive closed.
    pub active_indices: Vec<usize>,
    /// TOTAL relative midpoint work, independent of any admitted impulse gauge.
    pub relative_work: f64,
    /// Nonnegative TOTAL impact loss; no automatic heat or friction partition.
    pub dissipated_energy: f64,
    pub passivity_residual: f64,
    pub residuals: FrictionalContactResiduals,
    /// Weighted block maximum-norm impulse contraction bound, strictly below
    /// one when that witness applies. Otherwise one: full-rank energy admission
    /// uses `impulse_energy_witness`, and UniqueVelocity carries its own witness.
    /// This is numerical admission, not an interval proof.
    pub contraction_bound: f64,
    /// Positive weights [normal magnitude, tangential Euclidean norm] for each
    /// contact in active_indices order, allowing independent witness checking.
    pub contraction_weights: Vec<f64>,
    /// Optional global candidate energy proof for independent rows whose block
    /// contraction test failed. Its rank is the entire active dimension, so
    /// equal mobility images imply equal impulses. The original contact law,
    /// passivity checks and Newton solve are unchanged; this is not a branch-
    /// local Jacobian nonsingularity test or a compliance regularization.
    pub impulse_energy_witness: Option<PsdVelocityWitness>,
    /// Present for every `UniqueVelocity` result. Unique impulses already
    /// determine work and need no additional range witness.
    pub impulse_work_witness: Option<PsdImpulseWorkWitness>,
    /// Minimum SPD pivot, or zero for the admitted semidefinite range lane.
    pub minimum_normalized_pivot: f64,
    pub evaluations: usize,
    pub nonlinear: Option<NewtonReport>,
}

/// No variant publishes a partially accepted contact event.
#[derive(Debug)]
pub enum FrictionalContactError {
    Invalid(&'static str),
    PenetratingContact {
        index: usize,
        gap: f64,
    },
    /// Indefinite, rank-ambiguous, or unsupported semidefinite mobility retains
    /// the original SPD refusal. A numerical PSD range without a Coulomb
    /// determinacy witness instead reaches UnresolvedDeterminacy or Nonunique.
    /// None of these outcomes asserts physical nonexistence.
    MobilityNotAdmitted(SpdOrthantError),
    Newton(NewtonError),
    NotConverged(NewtonReport),
    SolverInterrupted(SolverRunProgress),
    Cancelled,
    EvaluationBudget,
    NonFinite,
    ResidualNotSatisfied(FrictionalContactResiduals),
    EnergyIncrease {
        work: f64,
        allowed: f64,
    },
    /// Determined velocity alone does not determine total work. The checked
    /// work covector, including any exact locked-slip reductions, retains an
    /// unresolved mobility-null component. No arbitrary loss is published.
    UnresolvedWork {
        /// Range projection or original-matrix work-covector residual [m/s].
        range_velocity_residual: f64,
        /// Numerical admission allowance [m/s]. Locked sliding covectors
        /// derive this from the all-root impulse bound and energy allowance.
        allowed: f64,
    },
    /// Both witnesses satisfy the physical law and passivity gates and have
    /// distinct mobility images, hence different contact velocities. A change
    /// only in redundant impulse allocation does not trigger this refusal.
    /// Values are in original input order; no arbitrary branch is selected.
    Nonunique {
        first: Vec<[f64; 3]>,
        second: Vec<[f64; 3]>,
    },
    /// A residual-checked passive root is evidence of a candidate, not a
    /// determination. The bounded search did not certify a unique response.
    UnresolvedDeterminacy {
        candidate_impulses: Vec<[f64; 3]>,
    },
}
impl core::fmt::Display for FrictionalContactError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "coupled frictional contact: {self:?}")
    }
}
impl std::error::Error for FrictionalContactError {}

/// Solve one frozen-geometry event using full row-major `3N x 3N` mobility.
/// Entries have units kg^-1; vectors are normal-first, contact-major. Every
/// mutation is private scratch. Cancellation, numerical failure, energy growth,
/// or unresolved determinacy returns no authority to update physical state.
pub fn solve_frictional_contacts(
    contacts: &[FrictionalContact],
    delassus: &[f64],
    config: FrictionalContactConfig,
    mut cancelled: impl FnMut() -> bool,
) -> Result<FrictionalContactResult, FrictionalContactError> {
    use FrictionalContactError as E;
    poll(&mut cancelled)?;
    validate_config(config)?;
    let n = contacts.len();
    if n > config.maximum_contacts
        || delassus.len() != 9 * n * n
        || delassus.iter().any(|v| !v.is_finite())
    {
        return Err(E::Invalid(
            "contact count, mobility shape or finite entries",
        ));
    }
    let mut order: Vec<_> = (0..n).collect();
    order.sort_unstable_by_key(|&i| contacts[i].id);
    if order
        .windows(2)
        .any(|w| contacts[w[0]].id == contacts[w[1]].id)
    {
        return Err(E::Invalid("duplicate contact ID"));
    }
    let mut targets = vec![0.0; n];
    for (i, c) in contacts.iter().enumerate() {
        poll(&mut cancelled)?;
        if !c.gap.is_finite()
            || !c.activation_tolerance.is_finite()
            || c.activation_tolerance < 0.0
            || c.relative_velocity.iter().any(|v| !v.is_finite())
            || !c.friction_coefficient.is_finite()
            || c.friction_coefficient < 0.0
        {
            return Err(E::Invalid("contact geometry, velocity or friction"));
        }
        if c.gap < 0.0 {
            return Err(E::PenetratingContact {
                index: i,
                gap: c.gap,
            });
        }
        targets[i] = match c.normal_law {
            NormalVelocityLaw::NewtonRestitution { coefficient }
                if coefficient.is_finite() && (0.0..=1.0).contains(&coefficient) =>
            {
                -coefficient * c.relative_velocity[0].min(0.0)
            }
            NormalVelocityLaw::Target { speed } if speed.is_finite() && speed >= 0.0 => speed,
            _ => return Err(E::Invalid("normal endpoint law")),
        };
    }
    let active: Vec<_> = order
        .iter()
        .copied()
        .filter(|&i| contacts[i].gap <= contacts[i].activation_tolerance)
        .collect();
    let d = 3 * active.len();
    let mut matrix = Vec::with_capacity(d * d);
    let mut free = Vec::with_capacity(d);
    let mut mu = Vec::with_capacity(active.len());
    for &i in &active {
        free.extend_from_slice(&contacts[i].relative_velocity);
        mu.push(contacts[i].friction_coefficient);
        for a in 0..3 {
            for &j in &active {
                for b in 0..3 {
                    matrix.push(delassus[(3 * i + a) * (3 * n) + 3 * j + b]);
                }
            }
        }
    }
    let has_locked_tangent_pairs = locked_tangents::has_locked_pair(&matrix, d);
    let (minimum_pivot, range_admission) =
        match solve_spd_orthant(&matrix, &vec![0.0; d], config.matrix, &mut cancelled) {
            Ok(admission) => (admission.min_normalized_pivot, None),
            Err(SpdOrthantError::Cancelled) => return Err(E::Cancelled),
            Err(
                error @ (SpdOrthantError::NotPositiveDefinite { .. }
                | SpdOrthantError::Degenerate { .. }),
            ) => {
                if (0..active.len()).any(|k| matrix[(3 * k) * d + 3 * k] <= 0.0) {
                    return Err(E::MobilityNotAdmitted(error));
                }
                let admission = range::admit(&matrix, d, config.matrix, &mut cancelled)?;
                match admission {
                    Some(admission) => (0.0, Some(admission)),
                    None => return Err(E::MobilityNotAdmitted(error)),
                }
            }
            Err(error) => return Err(E::MobilityNotAdmitted(error)),
        };
    let mut steps = Vec::with_capacity(d);
    for k in 0..active.len() {
        let a = 3 * k;
        let rn = 1.0 / matrix[a * d + a];
        // One common step preserves the Euclidean Coulomb disk in both axes.
        let tangent_diagonal = matrix[(a + 1) * d + a + 1] + matrix[(a + 2) * d + a + 2];
        // Any finite positive projection step defines the SAME Coulomb law.
        // A truly locked pair can use the existing normal mobility scale;
        // this is numerical solver scaling, never an entry added to W.
        let rt = if tangent_diagonal == 0.0 {
            rn
        } else {
            2.0 / tangent_diagonal
        };
        if !rn.is_finite() || !rt.is_finite() || rn <= 0.0 || rt <= 0.0 {
            return Err(E::NonFinite);
        }
        steps.extend_from_slice(&[rn, rt, rt]);
    }
    let mut range_witness = match &range_admission {
        Some(admission) => admission.normal_observable_witness(&mu, &mut cancelled)?,
        None => None,
    };
    let witness = if range_admission.is_some() {
        None
    } else {
        contraction_witness(&matrix, &steps, &mu, &mut cancelled)?
    };
    let evaluations = Cell::new(0);
    let cancel = RefCell::new(&mut cancelled);
    let active_targets: Vec<_> = active.iter().map(|&i| targets[i]).collect();
    let range_preconditioner = range_admission.is_some();
    let solve = |seed: Vec<f64>| -> Result<(Vec<f64>, NewtonReport), E> {
        let problem = natural::NaturalMap::new(
            &matrix,
            &free,
            &active_targets,
            &steps,
            &mu,
            &evaluations,
            config.maximum_evaluations,
            range_preconditioner,
            || (**cancel.borrow_mut())(),
        );
        let construction = NewtonKrylovState::new(&problem, seed, config.newton);
        if let Some(e) = problem.take_failure() {
            return Err(e);
        }
        let mut state = construction.map_err(E::Newton)?;
        let report = state.run_controlled(&problem, config.maximum_newton_iterations, &mut || {
            problem.has_failed() || (**cancel.borrow_mut())()
        });
        if let Some(e) = problem.take_failure() {
            return Err(e);
        }
        if report.progress == SolverRunProgress::Paused {
            return Err(E::Cancelled);
        }
        if report.progress != SolverRunProgress::Complete {
            return Err(E::SolverInterrupted(report.progress));
        }
        if !report.report.converged {
            return Err(E::NotConverged(report.report));
        }
        if has_locked_tangent_pairs {
            locked_tangents::complete_tractions(
                &matrix,
                &free,
                &mu,
                &mut state.x,
                &mut **cancel.borrow_mut(),
            )?;
        }
        Ok((state.x, report.report))
    };
    let (point, report) = if d == 0 {
        (vec![], None)
    } else {
        let (point, report) = solve(vec![0.0; d])?;
        (point, Some(report))
    };
    let mut result = audit(
        contacts,
        delassus,
        &order,
        &active,
        &targets,
        &point,
        config,
        &mut **cancel.borrow_mut(),
    )?;
    if range_witness.is_none() {
        if let Some(admission) = &range_admission {
            let velocity: Vec<_> = active
                .iter()
                .flat_map(|&i| result.relative_velocities_after[i])
                .collect();
            range_witness = normal_gauge::direction_witness(
                admission,
                &mu,
                &point,
                &velocity,
                None,
                config,
                &mut **cancel.borrow_mut(),
            )?;
            if range_witness.is_none() {
                let proposed = normal_gauge::coherent_representative(
                    admission,
                    &matrix,
                    &point,
                    &free,
                    config,
                    &mut **cancel.borrow_mut(),
                )?;
                if let Some((mut candidate, directions)) = proposed {
                    if has_locked_tangent_pairs {
                        locked_tangents::complete_tractions(
                            &matrix,
                            &free,
                            &mu,
                            &mut candidate,
                            &mut **cancel.borrow_mut(),
                        )?;
                    }
                    // A numerical gauge proposal never bypasses the original
                    // full matrix, unilateral law, circular cone or work gate.
                    let checked = audit(
                        contacts,
                        delassus,
                        &order,
                        &active,
                        &targets,
                        &candidate,
                        config,
                        &mut **cancel.borrow_mut(),
                    );
                    match checked {
                        Ok(representative) => {
                            let velocity: Vec<_> = active
                                .iter()
                                .flat_map(|&i| representative.relative_velocities_after[i])
                                .collect();
                            range_witness = normal_gauge::direction_witness(
                                admission,
                                &mu,
                                &candidate,
                                &velocity,
                                Some(&directions),
                                config,
                                &mut **cancel.borrow_mut(),
                            )?;
                            if range_witness.is_some() {
                                result = representative;
                            }
                        }
                        Err(E::Cancelled) => return Err(E::Cancelled),
                        _ => {}
                    }
                }
            }
        }
    }
    // A strong normal/tangent coupling can defeat the block norm bound even
    // when the law has exactly one impulse. Reuse the global candidate energy
    // argument on the FULL positive range before refusing that physical solve.
    // E + delta_p_n dot C delta_v_t <= 0 holds against every other Coulomb
    // root; positive energy coercivity forces delta_p=0 for an SPD mobility.
    // The spectral work is skipped on the existing contraction fast path.
    if range_admission.is_none() && witness.is_none() && d > 0 {
        let admission =
            range::admit_positive_definite(&matrix, d, config.matrix, &mut **cancel.borrow_mut())?;
        if let Some(admission) = admission {
            let velocity: Vec<_> = active
                .iter()
                .flat_map(|&i| result.relative_velocities_after[i])
                .collect();
            result.impulse_energy_witness = normal_gauge::direction_witness(
                &admission,
                &mu,
                &point,
                &velocity,
                None,
                config,
                &mut **cancel.borrow_mut(),
            )?;
        }
    }
    if let Some(witness) = range_witness {
        let accepted_point: Vec<_> = active.iter().flat_map(|&i| result.impulses[i]).collect();
        result.impulse_work_witness = Some(range::impulse_work_witness(
            // Only a checked semidefinite range creates range_witness above.
            range_admission
                .as_ref()
                .ok_or(E::Invalid("missing mobility range"))?,
            &matrix,
            &free,
            &mu,
            &accepted_point,
            config,
            &mut **cancel.borrow_mut(),
        )?);
        result.determinacy = FrictionalContactDeterminacy::UniqueVelocity(witness);
    } else if let Some((bound, weights)) = witness {
        result.contraction_bound = bound;
        result.contraction_weights = weights;
    } else if result.impulse_energy_witness.is_none() {
        // A second, genuinely loaded branch probe can reveal the separating /
        // compressive ambiguity that a zero/frictionless pair of seeds misses.
        let mut seed = vec![0.0; d];
        for k in 0..active.len() {
            let a = 3 * k;
            let tangential_speed = free[a + 1].hypot(free[a + 2]);
            seed[a] = steps[a]
                * (free[a].abs()
                    + targets[active[k]]
                    + tangential_speed
                    + config.velocity_tolerance);
            if tangential_speed > 0.0 {
                seed[a + 1] = -mu[k] * seed[a] * (free[a + 1] / tangential_speed);
                seed[a + 2] = -mu[k] * seed[a] * (free[a + 2] / tangential_speed);
            }
        }
        match solve(seed) {
            Ok((other, _)) => {
                match audit(
                    contacts,
                    delassus,
                    &order,
                    &active,
                    &targets,
                    &other,
                    config,
                    &mut **cancel.borrow_mut(),
                ) {
                    Ok(second)
                        if active.iter().any(|&i| {
                            (0..3).any(|a| {
                                (result.relative_velocities_after[i][a]
                                    - second.relative_velocities_after[i][a])
                                    .abs()
                                    > 16.0 * config.velocity_tolerance
                            })
                        }) =>
                    {
                        return Err(E::Nonunique {
                            first: result.impulses,
                            second: second.impulses,
                        });
                    }
                    Err(E::Cancelled) => return Err(E::Cancelled),
                    _ => {}
                }
            }
            Err(E::Cancelled) => return Err(E::Cancelled),
            Err(E::EvaluationBudget) => return Err(E::EvaluationBudget),
            Err(_) => {}
        }
        poll(&mut **cancel.borrow_mut())?;
        return Err(E::UnresolvedDeterminacy {
            candidate_impulses: result.impulses,
        });
    }
    result.minimum_normalized_pivot = minimum_pivot;
    result.evaluations = evaluations.get();
    result.nonlinear = report;
    poll(&mut **cancel.borrow_mut())?;
    Ok(result)
}

/// Bound the natural map in contact normal/tangential block norms. Disk
/// projection is nonexpansive in its argument and 1-Lipschitz in its radius.
/// Hence `delta T <= L delta p`; a positive w with Lw < w proves contraction
/// in max_i |delta p_i|/w_i. This is a sufficient condition, not a failed
/// uniqueness test when no witness is found. It performs no impulse iteration.
fn contraction_witness(
    w: &[f64],
    steps: &[f64],
    mu: &[f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Option<(f64, Vec<f64>)>, FrictionalContactError> {
    let n = mu.len();
    let d = 3 * n;
    let blocks = 2 * n;
    if n == 0 {
        return Ok(Some((0.0, vec![])));
    }
    let mut l = vec![0.0; blocks * blocks];
    for i in 0..n {
        poll(cancelled)?;
        let a = 3 * i;
        let rn = steps[a];
        let rt = steps[a + 1];
        for j in 0..n {
            let b = 3 * j;
            let diagonal = if i == j { 1.0 } else { 0.0 };
            l[(2 * i) * blocks + 2 * j] = (diagonal - rn * w[a * d + b]).abs();
            l[(2 * i) * blocks + 2 * j + 1] = rn * w[a * d + b + 1].hypot(w[a * d + b + 2]);
            l[(2 * i + 1) * blocks + 2 * j] = rt * w[(a + 1) * d + b].hypot(w[(a + 2) * d + b])
                + if i == j { mu[i] } else { 0.0 };
            // Exact 2x2 spectral norm, evaluated by two stable hypotenuses.
            // The Frobenius bound adds a spurious sqrt(2) for the +/-rho
            // diagonal residual of a directional mobility, and can refuse
            // even uncoupled, uniquely determined tangential responses.
            // The common roundoff inflation below keeps this a sufficient
            // Lipschitz certificate without another eigendecomposition.
            let x = diagonal - rt * w[(a + 1) * d + b + 1];
            let y = -rt * w[(a + 1) * d + b + 2];
            let z = -rt * w[(a + 2) * d + b + 1];
            let t = diagonal - rt * w[(a + 2) * d + b + 2];
            l[(2 * i + 1) * blocks + 2 * j + 1] = (0.5 * x + 0.5 * t).hypot(0.5 * z - 0.5 * y)
                + (0.5 * x - 0.5 * t).hypot(0.5 * z + 0.5 * y);
        }
    }
    let mut weights = vec![1.0; blocks];
    for _ in 0..128 {
        poll(cancelled)?;
        let image: Vec<_> = l
            .chunks_exact(blocks)
            .map(|row| fs_solver::dot(row, &weights))
            .collect();
        if image.iter().any(|v| !v.is_finite()) {
            return Ok(None);
        }
        let bound = image
            .iter()
            .zip(&weights)
            .map(|(a, b)| a / b)
            .fold(0.0, f64::max)
            * (1.0 + 256.0 * f64::EPSILON * (blocks as f64));
        if bound < 1.0 - 1e-8 {
            return Ok(Some((bound, weights)));
        }
        weights = image.into_iter().map(|v| 1.0 + v).collect();
    }
    Ok(None)
}

fn audit(
    contacts: &[FrictionalContact],
    w: &[f64],
    order: &[usize],
    active: &[usize],
    targets: &[f64],
    point: &[f64],
    config: FrictionalContactConfig,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<FrictionalContactResult, FrictionalContactError> {
    use FrictionalContactError as E;
    let n = contacts.len();
    let d = 3 * n;
    let mut impulses = vec![[0.0; 3]; n];
    for (k, &i) in active.iter().enumerate() {
        impulses[i].copy_from_slice(&point[3 * k..3 * k + 3]);
    }
    let mut velocities: Vec<_> = contacts.iter().map(|c| c.relative_velocity).collect();
    for &i in order {
        poll(cancelled)?;
        for a in 0..3 {
            for &j in active {
                for b in 0..3 {
                    velocities[i][a] += w[(3 * i + a) * d + 3 * j + b] * impulses[j][b];
                }
            }
        }
    }
    let mut r = FrictionalContactResiduals::default();
    if impulses
        .iter()
        .chain(&velocities)
        .flatten()
        .any(|v| !v.is_finite())
    {
        return Err(E::NonFinite);
    }
    for &i in active {
        let p = impulses[i];
        let v = velocities[i];
        let g = v[0] - targets[i];
        let radius = contacts[i].friction_coefficient * p[0];
        let speed = v[1].hypot(v[2]);
        if !g.is_finite() || !radius.is_finite() || !speed.is_finite() {
            return Err(E::NonFinite);
        }
        r.normal_impulse_violation = r.normal_impulse_violation.max((-p[0]).max(0.0));
        r.normal_velocity_violation = r.normal_velocity_violation.max((-g).max(0.0));
        r.normal_complementarity = r.normal_complementarity.max((p[0] * g).abs());
        r.cone_violation = r.cone_violation.max((p[1].hypot(p[2]) - radius).max(0.0));
        let duality_gap = (radius * speed + p[1] * v[1] + p[2] * v[2]).abs();
        if !duality_gap.is_finite() {
            return Err(E::NonFinite);
        }
        r.tangential_duality_gap = r.tangential_duality_gap.max(duality_gap);
        if speed > config.velocity_tolerance {
            r.sliding_impulse_residual = r
                .sliding_impulse_residual
                .max((p[1] + radius * (v[1] / speed)).hypot(p[2] + radius * (v[2] / speed)));
        }
    }
    let mut work = 0.0;
    let mut work_scale = 0.0;
    for &i in order {
        for a in 0..3 {
            let term =
                impulses[i][a] * (0.5 * contacts[i].relative_velocity[a] + 0.5 * velocities[i][a]);
            work += term;
            work_scale += term.abs();
        }
    }
    if impulses
        .iter()
        .chain(&velocities)
        .flatten()
        .any(|v| !v.is_finite())
        || !work.is_finite()
        || !work_scale.is_finite()
        || [
            r.normal_impulse_violation,
            r.normal_velocity_violation,
            r.normal_complementarity,
            r.cone_violation,
            r.sliding_impulse_residual,
            r.tangential_duality_gap,
        ]
        .iter()
        .any(|v| !v.is_finite())
    {
        return Err(E::NonFinite);
    }
    if r.normal_impulse_violation > config.impulse_tolerance
        || r.normal_velocity_violation > config.velocity_tolerance
        || r.normal_complementarity > config.complementarity_tolerance
        || r.cone_violation > config.impulse_tolerance
        || r.sliding_impulse_residual > config.impulse_tolerance
        || r.tangential_duality_gap > config.complementarity_tolerance
    {
        return Err(E::ResidualNotSatisfied(r));
    }
    let allowed = config.energy_tolerance + 64.0 * f64::EPSILON * work_scale;
    if work > allowed {
        return Err(E::EnergyIncrease { work, allowed });
    }
    poll(cancelled)?;
    Ok(FrictionalContactResult {
        determinacy: FrictionalContactDeterminacy::UniqueImpulses,
        impulses,
        relative_velocities_after: velocities,
        normal_targets: targets.to_vec(),
        active_indices: active.to_vec(),
        relative_work: work,
        dissipated_energy: (-work).max(0.0),
        passivity_residual: work.max(0.0),
        residuals: r,
        contraction_bound: 1.0,
        contraction_weights: vec![],
        impulse_energy_witness: None,
        impulse_work_witness: None,
        minimum_normalized_pivot: 0.0,
        evaluations: 0,
        nonlinear: None,
    })
}

fn validate_config(c: FrictionalContactConfig) -> Result<(), FrictionalContactError> {
    if !(1..=21).contains(&c.maximum_contacts)
        || !(1..=1_000_000).contains(&c.maximum_evaluations)
        || !(1..=256).contains(&c.maximum_newton_iterations)
        || [
            c.velocity_tolerance,
            c.impulse_tolerance,
            c.complementarity_tolerance,
        ]
        .iter()
        .any(|v| !v.is_finite() || *v <= 0.0)
        || !c.energy_tolerance.is_finite()
        || c.energy_tolerance < 0.0
        || c.newton.linear_restart > 64
        || c.newton.max_linear_cycles > 128
    {
        return Err(FrictionalContactError::Invalid("configuration"));
    }
    c.newton.validate().map_err(FrictionalContactError::Newton)
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), FrictionalContactError> {
    if cancelled() {
        Err(FrictionalContactError::Cancelled)
    } else {
        Ok(())
    }
}
