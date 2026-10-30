//! Orthotropic Coulomb friction with an associated, maximum-dissipation slip law.
//!
//! In the retained material axes the admissible static traction satisfies
//! `sum((F_i / (N mu_static_i))^2) <= 1`. During sliding the kinetic ellipse
//! supplies the support point in the **actual** slip direction:
//! `F_i = a_i^2 du_i / hypot(a_0 du_0, a_1 du_1)`, `a_i = N mu_kinetic_i`.
//! Force and slip need not be parallel. Independent component clamps or a
//! direction-dependent scalar coefficient do not implement this associated law.
//!
//! With isotropic body compliance C, `du = du_free - C F` reduces to a bounded
//! scalar projection onto the kinetic ellipse. Static feasibility is tested
//! first, including when the kinetic ellipse is smaller. The law allocates no
//! workspace, has zero elastic storage, and exports exactly `F dot du` as heat.
//! See the Abaqus Theory Guide, "Coulomb friction", anisotropic friction,
//! equations (10)--(11), for the ellipse and associated slip convention:
//! <https://docs.software.vt.edu/abaqusv2025/English/SIMACAETHERefMap/simathe-c-coulombfric.htm>.
//!
//! Coefficients are explicit caller data, not calibrated material properties.
//! Their axes are a supplied effective interface texture, not an inferred blend
//! of two surfaces. This rung excludes adhesion, elastic microslip, wear and
//! a thermal evolution equation. The mechanics owner supplies normal load and
//! transports the material axes with the existing objective-history protocol.

use fs_math::det;

use super::transport::{ObjectiveTangentialLaw, TangentialHistoryRotation};
use super::{TangentialLaw, TangentialStep, TangentialStepInput};
use crate::{FrictionRegime, TriboError, finite, nonnegative_finite, positive_finite};

mod compliance;
pub use compliance::SymmetricTangentialCompliance;

/// Continuous held-load relaxation with finite slip-arrest events.
pub mod relaxation;

/// Driven static hold, breakaway and bounded directional sliding evolution.
pub mod driven;

/// Instantaneous directional impulses without advancing physical slip history.
pub mod impulse;

/// Constant static/kinetic principal coefficients and an initial texture axis.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EllipticCoulomb {
    static_mu: [f64; 2],
    kinetic_mu: [f64; 2],
    initial_axis: [f64; 2],
}

/// Retained material direction and accumulated solved slip in the current frame.
#[derive(Clone, Debug, PartialEq)]
pub struct EllipticCoulombState {
    owner: EllipticCoulomb,
    axis: [f64; 2],
    slip_m: [f64; 2],
    contact_temperature_k: Option<f64>,
}

impl EllipticCoulombState {
    /// First material principal axis in the current tangent coordinates.
    #[must_use]
    pub const fn principal_axis(&self) -> [f64; 2] {
        self.axis
    }

    /// Accumulated actual slip, never the released motion predictor [m].
    #[must_use]
    pub const fn slip_m(&self) -> [f64; 2] {
        self.slip_m
    }

    /// Last explicit thermal boundary; this law does not evolve temperature.
    #[must_use]
    pub const fn contact_temperature_k(&self) -> Option<f64> {
        self.contact_temperature_k
    }
}

impl EllipticCoulomb {
    /// Admit nonnegative principal coefficients, with kinetic <= static on each
    /// axis. Zero coefficients admit a frictionless direction or whole plane.
    /// The nonzero finite axis is normalized once, at construction.
    pub fn new(
        static_mu: [f64; 2],
        kinetic_mu: [f64; 2],
        principal_axis: [f64; 2],
    ) -> Result<Self, TriboError> {
        for i in 0..2 {
            nonnegative_finite(static_mu[i], "elliptic_static_mu")?;
            nonnegative_finite(kinetic_mu[i], "elliptic_kinetic_mu")?;
            if kinetic_mu[i] > static_mu[i] {
                return Err(invalid("elliptic_kinetic_above_static"));
            }
            finite(principal_axis[i], "elliptic_principal_axis")?;
        }
        let length = det::hypot(principal_axis[0], principal_axis[1]);
        positive_finite(length, "elliptic_principal_axis_length")?;
        Ok(Self {
            static_mu,
            kinetic_mu,
            initial_axis: principal_axis.map(|x| x / length),
        })
    }

    /// Original ordered principal static coefficients.
    #[must_use]
    pub const fn static_coefficients(&self) -> [f64; 2] {
        self.static_mu
    }

    /// Original ordered principal kinetic coefficients.
    #[must_use]
    pub const fn kinetic_coefficients(&self) -> [f64; 2] {
        self.kinetic_mu
    }

    /// Continue explicit coefficient changes without advancing slip or texture.
    ///
    /// The source owner authorizes these zero-storage coefficient changes.
    /// Native ranges and checkpoint identity are checked atomically.
    pub fn continue_coefficients(
        &self,
        state: &EllipticCoulombState,
        static_mu: [f64; 2],
        kinetic_mu: [f64; 2],
    ) -> Result<(Self, EllipticCoulombState), TriboError> {
        self.check(state)?;
        let mut law = Self::new(static_mu, kinetic_mu, self.initial_axis)?;
        // Constructor normalization must not accumulate into source identity.
        law.initial_axis = self.initial_axis;
        let mut next = state.clone();
        next.owner = law;
        Ok((law, next))
    }

    /// Resolve an already sliding interval on the kinetic ellipse.
    ///
    /// The mechanical owner selects this branch from its accepted motion. The
    /// larger static ellipse cannot capture a moving midpoint merely because
    /// its released reaction lies below the breakaway threshold. Zero solved
    /// slip is still possible inside the kinetic ellipse; a dynamics owner must
    /// independently admit arrest with its impulse and event-time checks.
    /// The returned checkpoint retains this law's ORIGINAL static coefficients.
    pub fn advance_kinetic(
        &self,
        state: &EllipticCoulombState,
        input: &TangentialStepInput,
    ) -> Result<TangentialStep<EllipticCoulombState>, TriboError> {
        self.check(state)?;
        // Both coefficient pairs have already been validated. This temporary
        // capacity choice needs no normalization, allocation or source query.
        let kinetic = Self {
            static_mu: self.kinetic_mu,
            ..*self
        };
        let mut history = state.clone();
        history.owner = kinetic;
        let mut step = kinetic.advance(&history, input)?;
        step.state.owner = *self;
        Ok(step)
    }

    fn check(&self, state: &EllipticCoulombState) -> Result<(), TriboError> {
        if state.owner != *self {
            return Err(invalid("elliptic_foreign_checkpoint"));
        }
        Ok(())
    }
}

impl TangentialLaw for EllipticCoulomb {
    type State = EllipticCoulombState;

    fn model_id(&self) -> &str {
        "dry-elliptic-coulomb-associated-v1"
    }

    fn initial_state(&self) -> Self::State {
        Self::State {
            owner: *self,
            axis: self.initial_axis,
            slip_m: [0.0; 2],
            contact_temperature_k: None,
        }
    }

    fn stored_energy_j(&self, _state: &Self::State) -> f64 {
        0.0
    }

    fn advance(
        &self,
        state: &Self::State,
        input: &TangentialStepInput,
    ) -> Result<TangentialStep<Self::State>, TriboError> {
        self.check(state)?;
        input.solve.validate()?;
        nonnegative_finite(input.normal_n, "normal_n")?;
        nonnegative_finite(input.compliance_m_per_n, "compliance_m_per_n")?;
        positive_finite(input.dt_s, "dt_s")?;
        if let Some(t) = input.contact_temperature_k {
            positive_finite(t, "contact_temperature_k")?;
        }
        if let Some(p) = input.nominal_pressure_pa {
            nonnegative_finite(p, "nominal_pressure_pa")?;
        }
        let free = to_principal(state.axis, input.slip_increment_m)?;
        let static_caps = capacities(self.static_mu, input.normal_n)?;
        let kinetic_caps = capacities(self.kinetic_mu, input.normal_n)?;
        let (force, slip, regime, iterations) = resolve(
            free,
            static_caps,
            kinetic_caps,
            input.compliance_m_per_n,
            input.normal_n > 0.0,
            ProjectionControls::from(input.solve),
        )?;
        let traction_n = from_principal(state.axis, force)?;
        let slip_increment_m = if input.compliance_m_per_n == 0.0 {
            // Preserve the exact prescribed kinematic boundary and source path.
            input.slip_increment_m
        } else {
            from_principal(state.axis, slip)?
        };
        let force_norm = det::hypot(traction_n[0], traction_n[1]);
        if !force_norm.is_finite() || force_norm > input.solve.maximum_force_n {
            return Err(invalid("elliptic_maximum_force_n"));
        }
        let mut next = state.clone();
        next.contact_temperature_k = input.contact_temperature_k;
        let mut compatibility_residual_m = [0.0; 2];
        for i in 0..2 {
            next.slip_m[i] += slip_increment_m[i];
            finite(next.slip_m[i], "elliptic_accumulated_slip_m")?;
            compatibility_residual_m[i] = slip_increment_m[i]
                - (input.slip_increment_m[i] - input.compliance_m_per_n * traction_n[i]);
            finite(
                compatibility_residual_m[i],
                "elliptic_compatibility_residual_m",
            )?;
        }
        // Constitutive work is nonnegative term by term in the material axes.
        // Rotating a degenerate ellipse can put two cancelling contributions
        // in the world-frame dot product, even when physical heat is exactly
        // zero. Retain that isometry roundoff as a work residual, not heat or
        // a spurious negative-dissipation refusal.
        let dissipated_j = force[0] * slip[0] + force[1] * slip[1];
        nonnegative_finite(dissipated_j, "elliptic_dissipated_j")?;
        let actual_work_j =
            traction_n[0] * slip_increment_m[0] + traction_n[1] * slip_increment_m[1];
        finite(actual_work_j, "elliptic_actual_work_j")?;
        let closure_residual_j = dissipated_j - actual_work_j;
        finite(closure_residual_j, "elliptic_closure_residual_j")?;
        Ok(TangentialStep {
            traction_n,
            slip_increment_m,
            state: next,
            regime,
            stored_change_j: 0.0,
            dissipated_j,
            closure_residual_j,
            compatibility_residual_m,
            iterations,
        })
    }
}

impl ObjectiveTangentialLaw for EllipticCoulomb {
    fn transport_state(
        &self,
        state: &Self::State,
        rotation: TangentialHistoryRotation,
    ) -> Result<Self::State, TriboError> {
        self.check(state)?;
        let mut next = state.clone();
        next.axis = rotation.apply(state.axis)?;
        next.slip_m = rotation.apply(state.slip_m)?;
        Ok(next)
    }
}

fn invalid(field: &'static str) -> TriboError {
    TriboError::InvalidInput { field }
}

fn capacities(mu: [f64; 2], normal: f64) -> Result<[f64; 2], TriboError> {
    let result = mu.map(|x| x * normal);
    for x in result {
        finite(x, "elliptic_force_capacity_n")?;
    }
    Ok(result)
}

fn to_principal(axis: [f64; 2], vector: [f64; 2]) -> Result<[f64; 2], TriboError> {
    from_principal([axis[0], -axis[1]], vector)
}

fn from_principal(axis: [f64; 2], vector: [f64; 2]) -> Result<[f64; 2], TriboError> {
    let result = [
        axis[0] * vector[0] - axis[1] * vector[1],
        axis[1] * vector[0] + axis[0] * vector[1],
    ];
    for x in result {
        finite(x, "elliptic_tangent_vector")?;
    }
    Ok(result)
}

fn inside(force: [f64; 2], caps: [f64; 2]) -> bool {
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
    det::hypot(scaled[0], scaled[1]) <= 1.0
}

// Support point of the kinetic ellipse; scaling avoids squaring dimensional
// forces or multiplying dimensional force by displacement before normalization.
fn support(caps: [f64; 2], slip: [f64; 2]) -> Result<[f64; 2], TriboError> {
    let scale = caps[0].max(caps[1]);
    let motion = slip[0].abs().max(slip[1].abs());
    if scale == 0.0 || motion == 0.0 {
        return Ok([0.0; 2]);
    }
    let weighted = std::array::from_fn::<_, 2, _>(|i| (caps[i] / scale) * (slip[i] / motion));
    let norm = det::hypot(weighted[0], weighted[1]);
    if norm == 0.0 {
        return Ok([0.0; 2]);
    }
    let result = std::array::from_fn(|i| caps[i] * (weighted[i] / norm));
    for x in result {
        finite(x, "elliptic_sliding_traction_n")?;
    }
    Ok(result)
}

type Resolution = ([f64; 2], [f64; 2], FrictionRegime, usize);

// Stateless numerical budgets shared by force/displacement and impulse/velocity
// projections. Physical wrappers validate their own quantities and units.
#[derive(Clone, Copy)]
struct ProjectionControls {
    maximum_iterations: usize,
    absolute_tolerance: f64,
    relative_tolerance: f64,
}
impl From<super::TangentialSolveConfig> for ProjectionControls {
    fn from(solve: super::TangentialSolveConfig) -> Self {
        Self {
            maximum_iterations: solve.max_iterations,
            absolute_tolerance: solve.force_absolute_tolerance_n,
            relative_tolerance: solve.force_relative_tolerance,
        }
    }
}
impl ProjectionControls {
    fn tolerance(self, applied: f64, expected: f64) -> Result<f64, TriboError> {
        let result =
            self.absolute_tolerance + self.relative_tolerance * applied.abs().max(expected.abs());
        finite(result, "elliptic_projection_tolerance")?;
        Ok(result)
    }
}

fn resolve(
    free: [f64; 2],
    static_caps: [f64; 2],
    kinetic_caps: [f64; 2],
    c: f64,
    normal_active: bool,
    solve: ProjectionControls,
) -> Result<Resolution, TriboError> {
    if !normal_active {
        return Ok(([0.0; 2], free, FrictionRegime::Sliding, 0));
    }
    if free == [0.0; 2] {
        return Ok(([0.0; 2], [0.0; 2], FrictionRegime::Sticking, 0));
    }
    if c == 0.0 {
        return Ok((
            support(kinetic_caps, free)?,
            free,
            FrictionRegime::Sliding,
            0,
        ));
    }
    let required = free.map(|x| x / c);
    for x in required {
        finite(x, "elliptic_static_reaction_n")?;
    }
    if inside(required, static_caps) {
        return Ok((required, [0.0; 2], FrictionRegime::Sticking, 0));
    }
    // Circular ellipses and principal-axis motion have an exact radial answer.
    // Preserve that common fast path rather than iterating an affine load line.
    if kinetic_caps[0] == kinetic_caps[1] || free[0] == 0.0 || free[1] == 0.0 {
        let force = support(kinetic_caps, free)?;
        let slip = std::array::from_fn(|i| free[i] - c * force[i]);
        for value in slip {
            finite(value, "elliptic_resolved_slip_m")?;
        }
        return Ok((force, slip, FrictionRegime::Sliding, 0));
    }
    // Degenerate ellipses can stick along one direction while slipping freely
    // along the other. The admissible interior reaction then does zero work.
    let projected = std::array::from_fn(|i| {
        if kinetic_caps[i] == 0.0 {
            0.0
        } else {
            required[i]
        }
    });
    if inside(projected, kinetic_caps) {
        let slip = std::array::from_fn(|i| if kinetic_caps[i] == 0.0 { free[i] } else { 0.0 });
        return Ok((projected, slip, FrictionRegime::Sliding, 0));
    }
    let scale = kinetic_caps[0].max(kinetic_caps[1]);
    let distance = free[0].abs().max(free[1].abs());
    let r = kinetic_caps.map(|a| a / scale);
    let d = free.map(|x| x / distance);
    let h = c * (scale / distance);
    finite(h, "elliptic_scaled_compliance")?;
    for i in 0..2 {
        if kinetic_caps[i] > 0.0 && r[i] * r[i] == 0.0 {
            return Err(invalid("elliptic_unresolved_axis_scale"));
        }
    }
    // q = support_work / (scale * distance). Its exact root lies in [0, hi].
    // Each force component decreases monotonically with q. This is a scalar
    // projection, not a nested vector Newton solve or a sampled direction grid.
    let mut lo = 0.0;
    let mut hi = det::hypot(r[0] * d[0], r[1] * d[1]);
    for iteration in 1..=solve.maximum_iterations {
        let q = f64::midpoint(lo, hi);
        let mut force = [0.0; 2];
        let mut slip = [0.0; 2];
        let mut ellipse = [0.0; 2];
        for i in 0..2 {
            if kinetic_caps[i] > 0.0 {
                let denominator = q + h * r[i] * r[i];
                ellipse[i] = r[i] * d[i] / denominator;
                force[i] = kinetic_caps[i] * ellipse[i];
            }
            slip[i] = free[i] - c * force[i];
            finite(force[i], "elliptic_sliding_traction_n")?;
            finite(slip[i], "elliptic_resolved_slip_m")?;
        }
        let expected = support(kinetic_caps, slip)?;
        let residual = det::hypot(force[0] - expected[0], force[1] - expected[1]);
        let tolerance = solve.tolerance(
            det::hypot(force[0], force[1]),
            det::hypot(expected[0], expected[1]),
        )?;
        if residual <= tolerance && (0..2).all(|i| force[i] * slip[i] >= 0.0) && slip != [0.0; 2] {
            return Ok((force, slip, FrictionRegime::Sliding, iteration));
        }
        if q == lo || q == hi {
            break;
        }
        if det::hypot(ellipse[0], ellipse[1]) > 1.0 {
            lo = q;
        } else {
            hi = q;
        }
    }
    Err(invalid("elliptic_sliding_root_unresolved"))
}
