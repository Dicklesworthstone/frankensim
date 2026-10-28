//! Held-effort tangential ports, independent of contact geometry and pHS wrappers.
//!
//! `slip_increment_m` on input is the released (zero-traction) relative motion.
//! The isotropic body response is `du = du_free - compliance * F`. The effort
//! returned here enters the interface; the reactions on the bodies have opposite
//! signs. Work is **F dot du**, using the solved slip, never the released motion.
//! Coulomb and Stribeck have zero reversible storage; heat is exported rather
//! than accumulated in an invented thermal reservoir. Future compliant/layer
//! laws supply their own state and Hamiltonian through the same trait.

use crate::{
    DissipationStep, DryFrictionState, DryInterfaceSystemCard, FrictionLaw, FrictionRegime,
    HeatPartition, TriboError, checked_add, checked_mul, nonnegative_finite, positive_finite,
};
use fs_math::det;

/// Full two-axis body compliance around the original constitutive owner.
pub mod matrix;

fn finite(value: f64, field: &'static str) -> Result<f64, TriboError> {
    crate::finite(value, field)?;
    Ok(value)
}

/// Explicit budgets for the deterministic stick test and sliding bisection.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TangentialSolveConfig {
    /// At most 128 scalar root iterations; zero is refused.
    pub max_iterations: usize,
    /// Inclusive traction magnitude ceiling [N].
    pub maximum_force_n: f64,
    /// Absolute traction residual tolerance [N].
    pub force_absolute_tolerance_n: f64,
    /// Relative traction residual tolerance, strictly between zero and one.
    pub force_relative_tolerance: f64,
}

impl TangentialSolveConfig {
    pub(crate) fn validate(self) -> Result<(), TriboError> {
        if !(1..=128).contains(&self.max_iterations) {
            return Err(TriboError::InvalidInput {
                field: "max_iterations",
            });
        }
        positive_finite(self.maximum_force_n, "maximum_force_n")?;
        positive_finite(
            self.force_absolute_tolerance_n,
            "force_absolute_tolerance_n",
        )?;
        positive_finite(self.force_relative_tolerance, "force_relative_tolerance")?;
        if self.force_relative_tolerance >= 1.0 {
            return Err(TriboError::InvalidInput {
                field: "force_relative_tolerance",
            });
        }
        Ok(())
    }

    pub(crate) fn tolerance(self, applied: f64, expected: f64) -> Result<f64, TriboError> {
        finite(
            self.force_absolute_tolerance_n
                + self.force_relative_tolerance * applied.abs().max(expected.abs()),
            "force_tolerance_n",
        )
    }
}

/// One work-conjugate held-effort interval. All quantities use SI units.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TangentialStepInput {
    /// Compressive normal load [N], non-negative.
    pub normal_n: f64,
    /// Released driver-minus-body relative increment [m], in two tangent axes.
    pub slip_increment_m: [f64; 2],
    /// Positive held-effort interval [s].
    pub dt_s: f64,
    /// Caller thermal boundary [K]; mandatory for a domain-bounded card.
    pub contact_temperature_k: Option<f64>,
    /// Caller nominal pressure [Pa]; mandatory for a domain-bounded card.
    pub nominal_pressure_pa: Option<f64>,
    /// Isotropic held-force body response [m/N]; zero means prescribed motion.
    pub compliance_m_per_n: f64,
    /// Caller-owned numerical and effort budgets, without hidden defaults.
    pub solve: TangentialSolveConfig,
}

/// A stateful interface exposing reversible storage and a dissipative work port.
/// Implementations must return a candidate; a refusal must not mutate `s`.
pub trait TangentialLaw {
    /// Checkpoint state, including any slip, thermal and third-body layer state.
    type State: Clone + PartialEq + core::fmt::Debug;
    /// Stable model identity; this is not a material admission receipt.
    fn model_id(&self) -> &str;
    /// Unloaded constitutive state; thermal inputs are never defaulted.
    fn initial_state(&self) -> Self::State;
    /// Reversible interface energy [J], excluding exported heat.
    fn stored_energy_j(&self, s: &Self::State) -> f64;
    /// Resolves a held-effort interval and its storage/work/heat closure.
    fn advance(
        &self,
        s: &Self::State,
        x: &TangentialStepInput,
    ) -> Result<TangentialStep<Self::State>, TriboError>;
}

/// Dry-law checkpoint. No thermal evolution or third-body layer is inferred.
#[derive(Debug, Clone, PartialEq)]
pub struct DryTangentialState {
    slip_m: [f64; 2],
    contact_temperature_k: Option<f64>,
}

impl DryTangentialState {
    /// Accumulated solved slip in the caller's unchanged tangent frame [m].
    #[must_use]
    pub const fn slip_m(&self) -> [f64; 2] {
        self.slip_m
    }
    /// Last explicit thermal boundary [K], with absence preserved.
    #[must_use]
    pub const fn contact_temperature_k(&self) -> Option<f64> {
        self.contact_temperature_k
    }
}

/// Candidate constitutive step and held-effort dissipative port.
#[derive(Debug, Clone, PartialEq)]
pub struct TangentialStep<S> {
    /// Effort entering the interface [N]; body reactions use the opposite sign.
    pub traction_n: [f64; 2],
    /// Actual resolved relative increment [m].
    pub slip_increment_m: [f64; 2],
    /// Candidate checkpoint, published only by the caller's transaction.
    pub state: S,
    /// Static reaction or kinetic sliding branch.
    pub regime: FrictionRegime,
    /// Reversible storage increment [J].
    pub stored_change_j: f64,
    /// Exported dissipative work [J], non-negative.
    pub dissipated_j: f64,
    /// `stored_change_j + dissipated_j - F dot du` [J].
    pub closure_residual_j: f64,
    /// Kinematic equation residual [m], including static projection roundoff.
    pub compatibility_residual_m: [f64; 2],
    /// Scalar root iterations; zero for a static reaction, affine solution or endpoint root.
    pub iterations: usize,
}

impl<S> TangentialStep<S> {
    /// Routes heat into the existing checked work ledger without a new store.
    pub fn dissipation_step(
        &self,
        partition: HeatPartition,
    ) -> Result<DissipationStep, TriboError> {
        partition.validate()?;
        let step = DissipationStep {
            total_work_j: self.dissipated_j,
            surface_a_heat_j: checked_mul(
                self.dissipated_j,
                partition.surface_a,
                "surface_a_heat_j",
            )?,
            surface_b_heat_j: checked_mul(
                self.dissipated_j,
                partition.surface_b,
                "surface_b_heat_j",
            )?,
            other_work_j: checked_mul(self.dissipated_j, partition.other, "other_work_j")?,
        };
        step.validate()?;
        Ok(step)
    }
}

impl TangentialLaw for FrictionLaw {
    type State = DryTangentialState;

    fn model_id(&self) -> &str {
        match self {
            Self::Coulomb { .. } => "dry-coulomb-set-valued",
            Self::Stribeck { .. } => "dry-stribeck-set-valued-benchmark",
            Self::VelocityDependent { .. } => "dry-velocity-dependent-set-valued-benchmark",
            Self::Hyperbolic { .. } => "dry-hyperbolic-set-valued-benchmark",
        }
    }

    fn initial_state(&self) -> Self::State {
        DryTangentialState {
            slip_m: [0.0; 2],
            contact_temperature_k: None,
        }
    }

    fn stored_energy_j(&self, _s: &Self::State) -> f64 {
        0.0
    }

    fn advance(
        &self,
        s: &Self::State,
        x: &TangentialStepInput,
    ) -> Result<TangentialStep<Self::State>, TriboError> {
        self.validate()?;
        advance_with_capacity(
            s,
            x,
            x.normal_n == 0.0,
            checked_mul(self.static_mu(), x.normal_n, "static_limit")?,
            |speed| {
                let mu = if speed == 0.0 {
                    self.zero_speed_kinetic_coefficient()?
                } else {
                    self.kinetic_coefficient(speed)?
                };
                checked_mul(mu, x.normal_n, "traction_n")
            },
            |distance, compliance, dt| {
                self.affine_sliding_traction_n(x.normal_n, distance, compliance, dt)
            },
        )
    }
}

// Share the actual force-capacity solve. Adhesion does not become a fictitious
// normal load or coefficient; dry callers retain their original arithmetic.
pub(crate) fn advance_with_capacity(
    s: &DryTangentialState,
    x: &TangentialStepInput,
    released: bool,
    static_capacity_n: f64,
    kinetic_capacity_n: impl Fn(f64) -> Result<f64, TriboError>,
    affine_capacity_n: impl Fn(f64, f64, f64) -> Result<Option<f64>, TriboError>,
) -> Result<TangentialStep<DryTangentialState>, TriboError> {
    x.solve.validate()?;
    nonnegative_finite(x.normal_n, "normal_force")?;
    nonnegative_finite(x.compliance_m_per_n, "compliance_m_per_n")?;
    positive_finite(x.dt_s, "dt_s")?;
    if let Some(t) = x.contact_temperature_k {
        positive_finite(t, "temperature_kelvin")?;
    }
    if let Some(p) = x.nominal_pressure_pa {
        nonnegative_finite(p, "contact_pressure_pa")?;
    }
    for value in x.slip_increment_m {
        finite(value, "slip_increment_m")?;
    }

    // The scalar x-axis case preserves the modal owner's operation order.
    let distance = if x.slip_increment_m[1] == 0.0 {
        x.slip_increment_m[0].abs()
    } else {
        det::hypot(x.slip_increment_m[0], x.slip_increment_m[1])
    };
    finite(distance, "free_slip_distance_m")?;
    let (magnitude, regime, iterations) = solve_magnitude(
        x,
        distance,
        released,
        static_capacity_n,
        &kinetic_capacity_n,
        &affine_capacity_n,
    )?;
    let direction = if distance == 0.0 {
        [0.0; 2]
    } else {
        [
            x.slip_increment_m[0] / distance,
            x.slip_increment_m[1] / distance,
        ]
    };
    let mut traction_n = [0.0; 2];
    let mut slip_increment_m = [0.0; 2];
    let mut compatibility_residual_m = [0.0; 2];
    let mut slip_m = s.slip_m;
    for i in 0..2 {
        traction_n[i] = if magnitude == 0.0 {
            0.0
        } else {
            checked_mul(magnitude, direction[i], "traction_n")?
        };
        let response = checked_mul(x.compliance_m_per_n, traction_n[i], "body_response_m")?;
        let actual = finite(x.slip_increment_m[i] - response, "resolved_slip_m")?;
        slip_increment_m[i] = if regime == FrictionRegime::Sticking {
            0.0
        } else {
            actual
        };
        compatibility_residual_m[i] =
            finite(slip_increment_m[i] - actual, "compatibility_residual_m")?;
        slip_m[i] = checked_add(s.slip_m[i], slip_increment_m[i], "accumulated_slip_m")?;
    }
    if regime == FrictionRegime::Sliding
        && !released
        && distance > 0.0
        && slip_increment_m == [0.0; 2]
    {
        // A near-static endpoint may satisfy a loose force tolerance but
        // cannot be published as a kinetic branch at exactly zero slip.
        // The modal owner applies the same actual-state refusal after solve.
        return Err(TriboError::InvalidInput {
            field: "zero_sliding_increment",
        });
    }
    let work = checked_add(
        checked_mul(traction_n[0], slip_increment_m[0], "tangent_work_j")?,
        checked_mul(traction_n[1], slip_increment_m[1], "tangent_work_j")?,
        "tangent_work_j",
    )?;
    nonnegative_finite(work, "tangent_work_j")?;
    let state = DryTangentialState {
        slip_m,
        contact_temperature_k: x.contact_temperature_k,
    };
    let stored_change_j = 0.0;
    let dissipated_j = work;
    Ok(TangentialStep {
        traction_n,
        slip_increment_m,
        state,
        regime,
        stored_change_j,
        dissipated_j,
        closure_residual_j: stored_change_j + dissipated_j - work,
        compatibility_residual_m,
        iterations,
    })
}

impl TangentialLaw for DryInterfaceSystemCard {
    type State = DryTangentialState;
    fn model_id(&self) -> &str {
        self.friction_law().model_id()
    }
    fn initial_state(&self) -> Self::State {
        self.friction_law().initial_state()
    }
    fn stored_energy_j(&self, s: &Self::State) -> f64 {
        self.friction_law().stored_energy_j(s)
    }
    fn advance(
        &self,
        s: &Self::State,
        x: &TangentialStepInput,
    ) -> Result<TangentialStep<Self::State>, TriboError> {
        let temperature = x.contact_temperature_k.ok_or(TriboError::InvalidInput {
            field: "temperature_kelvin",
        })?;
        let pressure = x.nominal_pressure_pa.ok_or(TriboError::InvalidInput {
            field: "contact_pressure_pa",
        })?;
        // Refuse unavailable material state before any root work. Slip is checked
        // on the solved motion, not on the released predictor.
        self.applicability()
            .temperature_kelvin()
            .require_contains(temperature, "temperature_kelvin")?;
        self.applicability()
            .contact_pressure_pa()
            .require_contains(pressure, "contact_pressure_pa")?;
        let step = self.friction_law().advance(s, x)?;
        let speed = finite(
            det::hypot(step.slip_increment_m[0], step.slip_increment_m[1]) / x.dt_s,
            "slip_speed_mps",
        )?;
        self.validate_state(DryFrictionState::new(temperature, pressure, speed)?)?;
        Ok(step)
    }
}

impl FrictionLaw {
    /// Closed-form kinetic traction magnitude for an affine body load line.
    ///
    /// The caller must perform the set-valued static-friction test first and
    /// verify its force budget and constitutive residual afterwards. For a
    /// Coulomb law this is `mu_k * normal_n`. For a linear velocity law, solving
    /// `F = N * max(0, mu_0 + slope * (du_free - compliance * F) / dt)`
    /// gives the positive affine root, or zero on the clamped branch.
    ///
    /// Returns `None` for nonlinear laws or an ill-conditioned affine formula;
    /// callers retain their safeguarded root solver in those cases. A candidate
    /// is not a static reaction and does not change nonlinear branch selection.
    pub fn affine_sliding_traction_n(
        &self,
        normal_n: f64,
        free_slip_m: f64,
        compliance_m_per_n: f64,
        dt_s: f64,
    ) -> Result<Option<f64>, TriboError> {
        self.validate()?;
        nonnegative_finite(normal_n, "normal_force")?;
        positive_finite(free_slip_m, "free_slip_distance_m")?;
        nonnegative_finite(compliance_m_per_n, "compliance_m_per_n")?;
        positive_finite(dt_s, "dt_s")?;
        match *self {
            Self::Coulomb { kinetic_mu, .. } => {
                Ok(Some(checked_mul(kinetic_mu, normal_n, "traction_n")?))
            }
            Self::VelocityDependent {
                slope_per_speed, ..
            } => {
                let speed = finite(free_slip_m / dt_s, "slip_speed_mps")?;
                let released_force = checked_mul(self.kinetic_mu(speed)?, normal_n, "traction_n")?;
                if released_force == 0.0 {
                    return Ok(Some(0.0));
                }
                let gain = normal_n * slope_per_speed * (compliance_m_per_n / dt_s);
                let denominator = 1.0 + gain;
                // Near a singular weakening load line, bisection preserves the
                // declared residual tolerance without amplifying roundoff.
                if !gain.is_finite() || denominator <= 64.0 * f64::EPSILON * (1.0 + gain.abs()) {
                    return Ok(None);
                }
                let force = released_force / denominator;
                Ok(force.is_finite().then_some(force))
            }
            Self::Stribeck { .. } | Self::Hyperbolic { .. } => Ok(None),
        }
    }
}

// Operation order, endpoint limits and stopping criterion intentionally match
// fs-couple's modal friction scalar owner. Its in-crate comparison exercises both
// implementations bitwise; no modal/network type enters this constitutive leaf.
fn solve_magnitude(
    x: &TangentialStepInput,
    distance: f64,
    released: bool,
    cap: f64,
    kinetic_capacity_n: &impl Fn(f64) -> Result<f64, TriboError>,
    affine_capacity_n: &impl Fn(f64, f64, f64) -> Result<Option<f64>, TriboError>,
) -> Result<(f64, FrictionRegime, usize), TriboError> {
    if released {
        return Ok((0.0, FrictionRegime::Sliding, 0));
    }
    let compliance = x.compliance_m_per_n;
    if distance == 0.0 {
        return Ok((0.0, FrictionRegime::Sticking, 0));
    }
    if compliance == 0.0 {
        let force = kinetic_capacity_n(finite(distance / x.dt_s, "slip_speed_mps")?)?;
        require_force_budget(force, x.solve)?;
        return Ok((force, FrictionRegime::Sliding, 0));
    }
    let required = finite(distance / compliance, "sticking_force_n")?;
    if required <= cap {
        require_force_budget(required, x.solve)?;
        return Ok((required, FrictionRegime::Sticking, 0));
    }
    let evaluate = |magnitude: f64| -> Result<(f64, f64), TriboError> {
        let slip = finite(distance - compliance * magnitude, "resolved_slip_m")?.max(0.0);
        let expected = kinetic_capacity_n(finite(slip / x.dt_s, "slip_speed_mps")?)?;
        Ok((
            finite(magnitude - expected, "constitutive_residual_n")?,
            x.solve.tolerance(magnitude, expected)?,
        ))
    };
    if let Some(force) = affine_capacity_n(distance, compliance, x.dt_s)? {
        require_force_budget(force, x.solve)?;
        if force < required {
            let (residual, tolerance) = evaluate(force)?;
            if residual.abs() <= tolerance {
                return Ok((force, FrictionRegime::Sliding, 0));
            }
        }
    }
    if evaluate(0.0)?.0 == 0.0 {
        return Ok((0.0, FrictionRegime::Sliding, 0));
    }
    let mut lo = 0.0;
    let mut hi = required.min(x.solve.maximum_force_n);
    let (mut residual, mut tolerance) = evaluate(hi)?;
    if residual.abs() <= tolerance {
        return Ok((hi, FrictionRegime::Sliding, 0));
    }
    if residual < 0.0 {
        return Err(TriboError::InvalidInput {
            field: "friction root outside maximum_force_n",
        });
    }
    let mut used = 0;
    for iteration in 1..=x.solve.max_iterations {
        used = iteration;
        let mid = f64::midpoint(lo, hi);
        (residual, tolerance) = evaluate(mid)?;
        if residual.abs() <= tolerance {
            return Ok((mid, FrictionRegime::Sliding, iteration));
        }
        if mid == lo || mid == hi {
            break;
        }
        if residual < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    Err(TriboError::TangentialSolve {
        residual_n: residual,
        tolerance_n: tolerance,
        iterations: used,
    })
}

fn require_force_budget(force: f64, config: TangentialSolveConfig) -> Result<(), TriboError> {
    if force > config.maximum_force_n {
        return Err(TriboError::OutsideApplicability {
            field: "maximum_force_n",
            value: force,
            minimum: 0.0,
            maximum: config.maximum_force_n,
        });
    }
    Ok(())
}
