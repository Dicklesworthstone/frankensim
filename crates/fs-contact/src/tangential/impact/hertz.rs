//! Coupled normal restitution and finite Hertz sliding/spin impulses.
//!
//! The normal impulse is an unknown of the same event as the finite patch
//! friction. Each normal trial delegates its complete tangent/spin problem to
//! `fs-tribo`; the resulting normal cross-response closes the original Newton
//! restitution equation. This is the non-associated law, with no normal
//! dilation from an associated cone minimization and no independent spin clamp.
//!
//! The effective Hertz radius and normalized pressure-impulse shape are supplied
//! physical assumptions. This owner does not determine a transient footprint,
//! impact duration, normal compliance, thermal partition or multiple contacts.

use fs_math::dd::Dd;
use fs_solver::orthant::active_set::{
    PreparedSpdMatrix, SpdMatrixAdmission, SpdOrthantConfig, SpdOrthantError, prepare_spd_matrix,
};
use fs_tribo::sliding_patch::contact::impulse::{
    HertzImpulseControls, HertzImpulseError, HertzImpulseInput, HertzImpulsePatch,
    HertzImpulseReport,
};

/// Work-conjugate physical coordinates with positive normal speed separating.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HertzImpactInput {
    /// `[normal, tangent_x, tangent_y, radius * relative_spin]` [m/s].
    pub velocity_before_m_s: [f64; 4],
    /// Full symmetric positive-definite response to
    /// `[Jn, Jtx, Jty, normal_axis_angular_impulse / radius]` [1/kg].
    pub scaled_mobility_per_kg: [[f64; 4]; 4],
    /// Newton coefficient in `[0,1]`; normal target is `-e*min(vn_before,0)`.
    pub restitution_coefficient: f64,
}

/// Bounded nested work and independent physical publication gates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HertzImpactControls {
    /// Original conditional Hertz controls. Their force/couple impulse bounds
    /// also govern intermediate normal trials; no budget is silently enlarged.
    pub patch: HertzImpulseControls,
    /// Original numerical SPD admission, applied to four rows and a three-row
    /// Schur complement. No diagonal shift or compliance repairs a refusal.
    pub matrix: SpdOrthantConfig,
    /// At most this many safeguarded normal iterations, in `1..=512`, after
    /// at most two initial patch trials. An opening event uses one zero trial.
    pub maximum_iterations: usize,
    /// Total native Hertz support evaluations across all normal trials.
    pub maximum_evaluations: usize,
    /// Original inverse-mobility normal impulse reconstruction allowance [N s].
    /// A stalled normal bracket narrower than this still needs its velocity gate.
    pub normal_impulse_tolerance_ns: f64,
    /// Normal endpoint equation allowance [m/s]. No initial closing speed is
    /// rounded to zero by this tolerance.
    pub normal_velocity_tolerance_m_s: f64,
    /// Normal complementarity and independent assembled work closure [J].
    pub energy_tolerance_j: f64,
    /// Numerical normal-trial and returned normal-impulse bound [N s].
    pub maximum_normal_impulse_ns: f64,
    /// Required positive fractional margin for both global existence growth
    /// and the sufficient uniqueness witness, in `[4096*EPSILON,1)`.
    pub minimum_determinacy_margin: f64,
}

impl HertzImpactControls {
    /// Validate the native patch budget and the new outer physical gates.
    pub fn validate(self) -> Result<(), HertzImpactError> {
        self.patch
            .validate()
            .map_err(|error| HertzImpactError::Patch {
                normal_impulse_ns: 0.0,
                error,
            })?;
        if !(1..=512).contains(&self.maximum_iterations)
            || self.maximum_evaluations == 0
            || [
                self.normal_impulse_tolerance_ns,
                self.normal_velocity_tolerance_m_s,
                self.energy_tolerance_j,
            ]
            .iter()
            .any(|value| !value.is_finite() || *value <= 0.0)
            || !self.maximum_normal_impulse_ns.is_finite()
            || self.maximum_normal_impulse_ns < 0.0
            || !self.minimum_determinacy_margin.is_finite()
            || !(4096.0 * f64::EPSILON..1.0).contains(&self.minimum_determinacy_margin)
        {
            return Err(HertzImpactError::Invalid("coupled Hertz controls"));
        }
        Ok(())
    }
}

/// Sufficient floating-point admission of one globally determined response.
/// This is a numerical witness, not an interval certificate or a classification
/// of the general Painleve/nonunique region that this owner refuses.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HertzImpactDeterminacy {
    /// Full original four-coordinate SPD admission.
    pub matrix: SpdMatrixAdmission,
    /// Schur admission when nonzero friction and normal cross-coupling need it.
    pub schur: Option<SpdMatrixAdmission>,
    /// Lower growth estimate `Wnn - mu*|Wnt|`, including numerical margin [1/kg].
    pub normal_growth_lower_bound_per_kg: f64,
    /// `kappa < 1` in `E <= kappa*E` for any two law-satisfying roots.
    pub friction_feedback_upper_bound: f64,
    /// Maximum absolute Schur row sum, with roundoff margin [1/kg].
    pub schur_spectral_upper_bound_per_kg: f64,
    /// Infinity norm of the independently reconstructed Schur inverse residual.
    pub inverse_schur_residual: f64,
}

/// One complete event, still awaiting native body/application audits.
#[derive(Clone, Copy, Debug, PartialEq)]
#[must_use]
pub struct HertzImpactReport {
    /// `[Jn, Jtx, Jty, normal_axis_angular_impulse / radius]` [N s].
    pub scaled_impulse_ns: [f64; 4],
    /// Actual original `v_before + W*p`, without snapping arrest residuals [m/s].
    pub velocity_after_m_s: [f64; 4],
    /// Original conditional tangent/spin solve at the accepted normal impulse.
    /// Its `dissipated_j` is conditional work after a normal trial, NOT the
    /// complete event loss or a separately additive friction heat allocation.
    pub patch: HertzImpulseReport,
    pub normal_target_m_s: f64,
    /// Signed target defect for positive Jn; only a negative gap for Jn=0.
    pub normal_velocity_residual_m_s: f64,
    /// `Jn * (vn_after - normal_target)` [J].
    pub normal_complementarity_j: f64,
    /// Signed total loss from simultaneous physical midpoint work [J]. It is
    /// nonnegative within the declared energy allowance; a tiny negative
    /// roundoff value is retained, never clamped or converted into heat.
    pub dissipated_j: f64,
    /// Independently reconstructed in the complete original inverse mobility [J].
    pub kinetic_energy_change_j: f64,
    pub energy_residual_j: f64,
    /// Signed simultaneous normal midpoint work [J].
    pub normal_work_j: f64,
    /// Signed simultaneous shear-plus-spin midpoint work [J], not a local heat.
    pub friction_work_j: f64,
    pub determinacy: HertzImpactDeterminacy,
    pub iterations: usize,
    pub evaluations: usize,
    pub patch_solves: usize,
    /// Last numerical bracket [N s]. Inexact inner solves and floating
    /// arithmetic make this a diagnostic, not a validated interval enclosure.
    pub normal_impulse_bracket_ns: [f64; 2],
}

/// No refusal publishes an impulse candidate or changes a borrowed input.
#[derive(Clone, Debug, PartialEq)]
pub enum HertzImpactError {
    Invalid(&'static str),
    Mobility(SpdOrthantError),
    /// The unchanged original patch owner refused this named normal trial.
    /// A trial-budget failure is not proof that the final root exceeds it.
    Patch {
        normal_impulse_ns: f64,
        error: HertzImpulseError,
    },
    /// No sufficient global existence/uniqueness admission in this model rung.
    Determinacy {
        normal_growth_lower_bound_per_kg: f64,
        friction_feedback_upper_bound: f64,
    },
    /// The unique crossing lies beyond the allowed normal-trial interval,
    /// within the declared numerical endpoint resolution.
    NormalImpulseBudget {
        maximum_normal_impulse_ns: f64,
        normal_velocity_residual_m_s: f64,
    },
    EvaluationBudget {
        evaluations: usize,
        patch_solves: usize,
    },
    Unresolved {
        iterations: usize,
        normal_impulse_bracket_ns: [f64; 2],
        normal_velocity_residual_m_s: f64,
    },
    Endpoint(&'static str),
    /// Newton restitution and friction can together create energy even with
    /// `e <= 1`. Such a candidate is refused, never clipped or relabeled heat.
    EnergyBalance {
        dissipated_j: f64,
        kinetic_energy_change_j: f64,
        residual_j: f64,
    },
    Cancelled,
}

impl From<SpdOrthantError> for HertzImpactError {
    fn from(error: SpdOrthantError) -> Self {
        if error == SpdOrthantError::Cancelled {
            Self::Cancelled
        } else {
            Self::Mobility(error)
        }
    }
}
impl core::fmt::Display for HertzImpactError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "coupled normal/Hertz impact: {self:?}")
    }
}
impl std::error::Error for HertzImpactError {}

/// Solve normal restitution, sliding and spin with the complete response.
///
/// With `W = [a,b';b,C]`, a trial `Jn` sends `v_t + b*Jn` and `C` to the
/// original Hertz impulse owner. Its returned q must close
/// `vn + a*Jn + b'*q = -e*min(vn,0)`. A sufficient global witness is admitted
/// before searching: a conservative positive normal growth supplies existence,
/// and an energy argument compares the candidate with EVERY possible root.
/// Outside this sufficient region the function refuses rather than selecting
/// an arbitrary root. Opening contact returns a zero impulse.
///
/// Original SPD numerical owners retain both matrix factors; the outer search
/// needs only fixed four-coordinate scratch. Nested support evaluations share
/// one budget. Original endpoint, normal complementarity, inverse-mobility
/// reconstruction, and total passivity gates govern the final return.
pub fn solve_hertz_impact(
    patch: &HertzImpulsePatch,
    input: &HertzImpactInput,
    controls: HertzImpactControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<HertzImpactReport, HertzImpactError> {
    poll(&mut cancelled)?;
    controls.validate()?;
    if input.velocity_before_m_s.iter().any(|v| !v.is_finite())
        || !input.restitution_coefficient.is_finite()
        || !(0.0..=1.0).contains(&input.restitution_coefficient)
    {
        return Err(HertzImpactError::Invalid("velocity or Newton coefficient"));
    }
    let flat: [f64; 16] = core::array::from_fn(|i| input.scaled_mobility_per_kg[i / 4][i % 4]);
    let factor = prepare_spd_matrix(&flat, 4, controls.matrix, &mut cancelled)?;
    let determinacy = determinacy(patch, input, controls, factor.admission(), &mut cancelled)?;
    let target = -input.restitution_coefficient * input.velocity_before_m_s[0].min(0.0);
    let mut work = Work {
        evaluations: 0,
        patch_solves: 0,
    };
    if input.velocity_before_m_s[0] >= 0.0 {
        let trial = evaluate(
            patch,
            input,
            controls,
            target,
            0.0,
            &mut work,
            &mut cancelled,
        )?;
        return finish(
            input,
            controls,
            &factor,
            determinacy,
            trial,
            target,
            work,
            0,
            [0.0; 2],
            &mut cancelled,
        );
    }
    // Closing motion always reaches a nonzero normal trial; no velocity
    // tolerance invents a dead band in the impact law.
    let demand = (Dd::from_f64(target) - Dd::from_f64(input.velocity_before_m_s[0])).to_f64();
    let upper = demand / determinacy.normal_growth_lower_bound_per_kg;
    if !demand.is_finite() || !upper.is_finite() || upper <= 0.0 {
        return Err(HertzImpactError::Invalid("representable normal bracket"));
    }
    let bound = upper.min(controls.maximum_normal_impulse_ns);
    if bound == 0.0 {
        return Err(HertzImpactError::NormalImpulseBudget {
            maximum_normal_impulse_ns: bound,
            normal_velocity_residual_m_s: -demand,
        });
    }
    let first = (demand / input.scaled_mobility_per_kg[0][0]).min(bound);
    if first == 0.0 || !first.is_finite() {
        return Err(HertzImpactError::Invalid("representable normal impulse"));
    }
    let mut trial = evaluate(
        patch,
        input,
        controls,
        target,
        first,
        &mut work,
        &mut cancelled,
    )?;
    let mut bracket = [0.0, bound];
    if admitted_normal(trial, controls) {
        return finish(
            input,
            controls,
            &factor,
            determinacy,
            trial,
            target,
            work,
            0,
            bracket,
            &mut cancelled,
        );
    }
    let mut gaps = [-demand, trial.gap];
    if trial.gap > 0.0 {
        bracket[1] = first;
    } else {
        bracket[0] = first;
        gaps[0] = trial.gap;
        if first != bound {
            trial = evaluate(
                patch,
                input,
                controls,
                target,
                bound,
                &mut work,
                &mut cancelled,
            )?;
            if admitted_normal(trial, controls) {
                return finish(
                    input,
                    controls,
                    &factor,
                    determinacy,
                    trial,
                    target,
                    work,
                    0,
                    bracket,
                    &mut cancelled,
                );
            }
        }
        gaps[1] = trial.gap;
        if trial.gap < 0.0 {
            return Err(if bound == controls.maximum_normal_impulse_ns {
                HertzImpactError::NormalImpulseBudget {
                    maximum_normal_impulse_ns: bound,
                    normal_velocity_residual_m_s: trial.gap,
                }
            } else {
                HertzImpactError::Unresolved {
                    iterations: 0,
                    normal_impulse_bracket_ns: bracket,
                    normal_velocity_residual_m_s: trial.gap,
                }
            });
        }
    }
    for iteration in 1..=controls.maximum_iterations {
        poll(&mut cancelled)?;
        let width = bracket[1] - bracket[0];
        // Bracketed secant when useful; a central safeguard retains progress
        // through changes between sliding and complete patch arrest.
        let fraction = (-gaps[0] / (gaps[1] - gaps[0])).clamp(0.1, 0.9);
        let next = bracket[0] + fraction * width;
        if !next.is_finite() || next == bracket[0] || next == bracket[1] {
            return Err(HertzImpactError::Unresolved {
                iterations: iteration - 1,
                normal_impulse_bracket_ns: bracket,
                normal_velocity_residual_m_s: trial.gap,
            });
        }
        trial = evaluate(
            patch,
            input,
            controls,
            target,
            next,
            &mut work,
            &mut cancelled,
        )?;
        if admitted_normal(trial, controls) {
            return finish(
                input,
                controls,
                &factor,
                determinacy,
                trial,
                target,
                work,
                iteration,
                bracket,
                &mut cancelled,
            );
        }
        let side = usize::from(trial.gap > 0.0);
        bracket[side] = next;
        gaps[side] = trial.gap;
    }
    Err(HertzImpactError::Unresolved {
        iterations: controls.maximum_iterations,
        normal_impulse_bracket_ns: bracket,
        normal_velocity_residual_m_s: trial.gap,
    })
}

// Let K be the unit-normal, unit-mu Hertz wrench set in radius-scaled
// coordinates. Its support is h(t,z)=E|t+z*Jr/a|. Cauchy--Schwarz and the Hertz
// moments E[r]=0, E[|r/a|²]=2/5 give h(t,z)² <= |t|²+(2/5)z² <= |(t,z)|².
// Thus K is inside the unit Euclidean ball, including mixed sliding/spin.
//
// For two roots, normal equality gives dn=-b'*dq/a and dv_t=S*dq with
// S=C-b*b'/a. Monotonicity of the support subgradient gives
// E=dq'*S*dq <= mu*|dn|*|S*dq|
//                <= (mu/a)*sqrt(b'*S^-1*b)*sqrt(lambda_max(S))*E.
// A kappa strictly below one proves dq=dn=0. The row-sum spectral bound and
// independently reconstructed inverse residual make this a sufficient
// floating admission, not an interval proof. No local Jacobian is substituted.
fn determinacy(
    patch: &HertzImpulsePatch,
    input: &HertzImpactInput,
    controls: HertzImpactControls,
    matrix: SpdMatrixAdmission,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<HertzImpactDeterminacy, HertzImpactError> {
    let w = input.scaled_mobility_per_kg;
    let a = w[0][0];
    let mu = patch.friction_coefficient();
    let b: [f64; 3] = core::array::from_fn(|i| 0.5 * w[0][i + 1] + 0.5 * w[i + 1][0]);
    let normal_row = [w[0][1], w[0][2], w[0][3]];
    let growth = (Dd::from_f64(a)
        - Dd::from_f64(mu) * Dd::from_f64(norm(normal_row) * (1.0 + 256.0 * f64::EPSILON))
        - Dd::from_f64(a) * Dd::from_f64(256.0 * f64::EPSILON))
    .to_f64();
    let schur: [[f64; 3]; 3] = core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            (Dd::from_f64(0.5 * w[i + 1][j + 1] + 0.5 * w[j + 1][i + 1])
                - Dd::from_f64(b[i]) * (Dd::from_f64(b[j]) / Dd::from_f64(a)))
            .to_f64()
        })
    });
    let spectral = schur
        .map(|row| row.into_iter().map(f64::abs).sum::<f64>())
        .into_iter()
        .fold(0.0, f64::max)
        * (1.0 + 256.0 * f64::EPSILON);
    let mut out = HertzImpactDeterminacy {
        matrix,
        schur: None,
        normal_growth_lower_bound_per_kg: growth,
        friction_feedback_upper_bound: 0.0,
        schur_spectral_upper_bound_per_kg: spectral,
        inverse_schur_residual: 0.0,
    };
    if !growth.is_finite() || growth / a < controls.minimum_determinacy_margin {
        return Err(HertzImpactError::Determinacy {
            normal_growth_lower_bound_per_kg: growth,
            friction_feedback_upper_bound: f64::INFINITY,
        });
    }
    if !spectral.is_finite() {
        return Err(HertzImpactError::Endpoint("finite Schur spectral witness"));
    }
    if mu == 0.0 || b == [0.0; 3] {
        return Ok(out);
    }
    let flat: [f64; 9] = core::array::from_fn(|i| schur[i / 3][i % 3]);
    let factor = prepare_spd_matrix(&flat, 3, controls.matrix, &mut *cancelled)?;
    out.schur = Some(factor.admission());
    let mut inverse = [[0.0; 3]; 3];
    for column in 0..3 {
        let mut unit = [0.0; 3];
        unit[column] = 1.0;
        let solution = factor.solve(&unit, &mut *cancelled)?;
        for row in 0..3 {
            inverse[row][column] = solution[row];
        }
    }
    let inverse_norm = inverse
        .map(|row| row.into_iter().map(f64::abs).sum::<f64>())
        .into_iter()
        .fold(0.0, f64::max);
    let mut rho: f64 = 0.0;
    for row in 0..3 {
        poll(cancelled)?;
        let mut row_sum = 0.0;
        for column in 0..3 {
            let product = (0..3).fold(Dd::ZERO, |sum, k| {
                sum + Dd::from_f64(schur[row][k]) * Dd::from_f64(inverse[k][column])
            });
            row_sum += (Dd::from_f64(if row == column { 1.0 } else { 0.0 }) - product)
                .to_f64()
                .abs();
        }
        rho = rho.max(row_sum);
    }
    out.inverse_schur_residual = rho;
    if !rho.is_finite() || rho > 1e-8 || !inverse_norm.is_finite() || !spectral.is_finite() {
        return Err(HertzImpactError::Endpoint("Schur inverse witness"));
    }
    let solution = inverse.map(|row| dot(row, b).to_f64());
    let residual =
        core::array::from_fn(|i| (Dd::from_f64(b[i]) - dot(schur[i], solution)).to_f64());
    let quadratic =
        dot(b, solution).to_f64() + norm(b) * (inverse_norm / (1.0 - rho)) * norm(residual);
    let quadratic = quadratic * (1.0 + 256.0 * f64::EPSILON);
    let feedback = mu * (quadratic / a).sqrt() * (spectral / a).sqrt();
    out.friction_feedback_upper_bound = feedback;
    if !feedback.is_finite() || feedback > 1.0 - controls.minimum_determinacy_margin {
        return Err(HertzImpactError::Determinacy {
            normal_growth_lower_bound_per_kg: growth,
            friction_feedback_upper_bound: feedback,
        });
    }
    Ok(out)
}

#[derive(Clone, Copy)]
struct Work {
    evaluations: usize,
    patch_solves: usize,
}
#[derive(Clone, Copy)]
struct Trial {
    impulse: [f64; 4],
    after: [f64; 4],
    patch: HertzImpulseReport,
    radius_m: f64,
    gap: f64,
}

#[allow(clippy::too_many_arguments)]
fn evaluate(
    patch: &HertzImpulsePatch,
    input: &HertzImpactInput,
    controls: HertzImpactControls,
    target: f64,
    normal: f64,
    work: &mut Work,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Trial, HertzImpactError> {
    poll(cancelled)?;
    let remaining = controls.maximum_evaluations - work.evaluations;
    if remaining == 0 {
        return Err(HertzImpactError::EvaluationBudget {
            evaluations: work.evaluations,
            patch_solves: work.patch_solves,
        });
    }
    let w = input.scaled_mobility_per_kg;
    let free: [f64; 3] = core::array::from_fn(|i| {
        (Dd::from_f64(input.velocity_before_m_s[i + 1])
            + Dd::from_f64(w[i + 1][0]) * Dd::from_f64(normal))
        .to_f64()
    });
    let spin = free[2] / patch.radius_m();
    if free.iter().any(|v| !v.is_finite()) || !spin.is_finite() || (free[2] != 0.0 && spin == 0.0) {
        return Err(HertzImpactError::Invalid(
            "representable conditional velocity",
        ));
    }
    let resolved = patch
        .solve(
            &HertzImpulseInput {
                normal_impulse_ns: normal,
                velocity_before_m_s: [free[0], free[1]],
                spin_before_rad_s: spin,
                scaled_mobility_per_kg: core::array::from_fn(|i| {
                    core::array::from_fn(|j| w[i + 1][j + 1])
                }),
            },
            HertzImpulseControls {
                maximum_evaluations: controls.patch.maximum_evaluations.min(remaining),
                ..controls.patch
            },
            &mut *cancelled,
        )
        .map_err(|error| {
            if error == HertzImpulseError::Cancelled {
                HertzImpactError::Cancelled
            } else {
                HertzImpactError::Patch {
                    normal_impulse_ns: normal,
                    error,
                }
            }
        })?;
    work.evaluations += resolved.evaluations;
    work.patch_solves += 1;
    let impulse = [
        normal,
        resolved.impulse_ns[0],
        resolved.impulse_ns[1],
        resolved.spin_impulse_nms / patch.radius_m(),
    ];
    let after = core::array::from_fn(|i| {
        (Dd::from_f64(input.velocity_before_m_s[i]) + dot(w[i], impulse)).to_f64()
    });
    if impulse.iter().chain(&after).any(|x| !x.is_finite()) {
        return Err(HertzImpactError::Endpoint("finite coupled endpoint"));
    }
    let patch_after = [
        resolved.velocity_after_m_s[0],
        resolved.velocity_after_m_s[1],
        patch.radius_m() * resolved.spin_after_rad_s,
    ];
    if (0..3).any(|i| (after[i + 1] - patch_after[i]).abs() > controls.patch.velocity_tolerance_m_s)
    {
        return Err(HertzImpactError::Endpoint(
            "original coupled tangent/spin endpoint",
        ));
    }
    Ok(Trial {
        impulse,
        after,
        patch: resolved,
        radius_m: patch.radius_m(),
        gap: (Dd::from_f64(after[0]) - Dd::from_f64(target)).to_f64(),
    })
}

fn admitted_normal(trial: Trial, controls: HertzImpactControls) -> bool {
    trial.gap.abs() <= controls.normal_velocity_tolerance_m_s
        && (trial.impulse[0] * trial.gap).abs() <= controls.energy_tolerance_j
}

#[allow(clippy::too_many_arguments)]
fn finish(
    input: &HertzImpactInput,
    controls: HertzImpactControls,
    factor: &PreparedSpdMatrix,
    determinacy: HertzImpactDeterminacy,
    trial: Trial,
    target: f64,
    work: Work,
    iterations: usize,
    bracket: [f64; 2],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<HertzImpactReport, HertzImpactError> {
    poll(cancelled)?;
    let normal_residual = if trial.impulse[0] == 0.0 {
        trial.gap.min(0.0)
    } else {
        trial.gap
    };
    let complementarity = trial.impulse[0] * trial.gap;
    if normal_residual.abs() > controls.normal_velocity_tolerance_m_s
        || !complementarity.is_finite()
        || complementarity.abs() > controls.energy_tolerance_j
        || trial.impulse[0] < 0.0
        || trial.impulse[0] > controls.maximum_normal_impulse_ns
    {
        return Err(HertzImpactError::Endpoint("normal complementarity"));
    }
    let midpoint = core::array::from_fn::<_, 4, _>(|i| {
        (Dd::from_f64(input.velocity_before_m_s[i]) + Dd::from_f64(trial.after[i]))
            * Dd::from_f64(0.5)
    });
    let normal_work = Dd::from_f64(trial.impulse[0]) * midpoint[0];
    let friction_work = (1..4).fold(Dd::ZERO, |sum, i| {
        sum + Dd::from_f64(trial.impulse[i]) * midpoint[i]
    });
    let total_work = normal_work + friction_work;
    let change = core::array::from_fn::<_, 4, _>(|i| {
        (Dd::from_f64(trial.after[i]) - Dd::from_f64(input.velocity_before_m_s[i])).to_f64()
    });
    let reconstructed = factor.solve(&change, &mut *cancelled)?;
    let kinetic = (0..4).fold(Dd::ZERO, |sum, i| {
        sum + Dd::from_f64(reconstructed[i]) * midpoint[i]
    });
    let residual = (kinetic - total_work).to_f64();
    let dissipated = -total_work.to_f64();
    if !dissipated.is_finite()
        || !kinetic.to_f64().is_finite()
        || !residual.is_finite()
        || dissipated < -controls.energy_tolerance_j
        || kinetic.to_f64() > controls.energy_tolerance_j
        || residual.abs() > controls.energy_tolerance_j
    {
        return Err(HertzImpactError::EnergyBalance {
            dissipated_j: dissipated,
            kinetic_energy_change_j: kinetic.to_f64(),
            residual_j: residual,
        });
    }
    // The original four-dimensional inverse audit is distinct from the native
    // three-dimensional patch solve and cannot borrow its conditional loss.
    let tangent_error =
        (reconstructed[1] - trial.impulse[1]).hypot(reconstructed[2] - trial.impulse[2]);
    if (reconstructed[0] - trial.impulse[0]).abs() > controls.normal_impulse_tolerance_ns
        || tangent_error > controls.patch.impulse_tolerance_ns
        || (trial.radius_m * (reconstructed[3] - trial.impulse[3])).abs()
            > controls.patch.spin_impulse_tolerance_nms
    {
        return Err(HertzImpactError::Endpoint(
            "original inverse-mobility impulse",
        ));
    }
    poll(cancelled)?;
    Ok(HertzImpactReport {
        scaled_impulse_ns: trial.impulse,
        velocity_after_m_s: trial.after,
        patch: trial.patch,
        normal_target_m_s: target,
        normal_velocity_residual_m_s: normal_residual,
        normal_complementarity_j: complementarity,
        dissipated_j: dissipated,
        kinetic_energy_change_j: kinetic.to_f64(),
        energy_residual_j: residual,
        normal_work_j: normal_work.to_f64(),
        friction_work_j: friction_work.to_f64(),
        determinacy,
        iterations,
        evaluations: work.evaluations,
        patch_solves: work.patch_solves,
        normal_impulse_bracket_ns: bracket,
    })
}

fn dot<const N: usize>(a: [f64; N], b: [f64; N]) -> Dd {
    a.into_iter().zip(b).fold(Dd::ZERO, |sum, (a, b)| {
        sum + Dd::from_f64(a) * Dd::from_f64(b)
    })
}
fn norm(v: [f64; 3]) -> f64 {
    v[0].hypot(v[1]).hypot(v[2])
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), HertzImpactError> {
    if cancelled() {
        Err(HertzImpactError::Cancelled)
    } else {
        Ok(())
    }
}
