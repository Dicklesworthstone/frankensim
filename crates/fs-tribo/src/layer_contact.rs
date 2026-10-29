//! Isothermal elastic-viscoplastic shear of a retained third-body layer.
//!
//! This is a homogeneous, pressure-independent shear layer, not a measured
//! rosin law. In two tangent axes, total slip is e + p and the continuous law is
//! F = k e, p_dot = (|F|-Y)_+ F/(c |F|). The held effort uses the midpoint
//! elastic strain: F = k(e_old+e_new)/2. With du = du_free-C F, the radial
//! return has a closed form; no velocity regularization or fitted spectrum is
//! introduced. Work splits into delta(k |e|^2/2) + Y |dp| + c |dp|^2/dt.
//!
//! The contact must remain compressive. Opening, geometry, wear and layer
//! advection need their own state transfer; none may silently erase this stored
//! shear energy. Frame-only changes use `tangential::transport`, retaining both
//! elastic and plastic vectors and the original constitutive owner.
//! Thermal owners may use `continue_state` with an unchanged elastic potential,
//! or `continue_state_with_work` when they account explicitly for coefficient
//! work and the mechanical entropy. Neither operation resets accepted strain.

use crate::{FrictionRegime, TangentialLaw, TangentialStep, TangentialStepInput, TriboError,
    checked_add, checked_mul, finite, nonnegative_finite, positive_finite};
use fs_math::det;

/// Frozen coefficients integrated from layer area and thickness by the source
/// owner: k=A G/h [N/m], Y=A tau_y [N], c=A eta_plastic/h [N s/m].
/// These numerical inputs alone carry no measured-material authority.
#[derive(Clone, Debug, PartialEq)]
pub struct ThirdBodyLayerContact {
    stiffness_n_m: f64,
    yield_force_n: f64,
    plastic_resistance_n_s_m: f64,
    temperature_kelvin: f64,
    maximum_elastic_slip_m: f64,
}

/// Accepted elastic and plastic slip in the checkpoint's current tangent frame.
/// The private owner binding prevents reusing stored strain with another law.
#[derive(Clone, Debug, PartialEq)]
pub struct ThirdBodyLayerState {
    elastic_slip_m: [f64; 2],
    plastic_slip_m: [f64; 2],
    owner: ThirdBodyLayerContact,
}
impl ThirdBodyLayerState {
    #[must_use]
    pub const fn elastic_slip_m(&self) -> [f64; 2] { self.elastic_slip_m }
    #[must_use]
    pub const fn plastic_slip_m(&self) -> [f64; 2] { self.plastic_slip_m }
}

/// Thermodynamic mechanical contribution for an admitted affine k(T).
/// The thermal owner adds its independent stress-free caloric internal energy.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThirdBodyLayerThermodynamics {
    pub free_energy_j: f64,
    pub entropy_j_k: f64,
    pub internal_energy_j: f64,
}

/// Retained strain at new coefficients, with the work of changing its potential.
#[derive(Clone, Debug, PartialEq)]
pub struct ThirdBodyLayerContinuation {
    pub state: ThirdBodyLayerState,
    /// psi(T_new,e_old) - psi(T_old,e_old), not dissipated heat [J].
    pub coefficient_free_energy_change_j: f64,
}

impl ThirdBodyLayerContact {
    /// Bind an explicit finite shear domain; source owners derive its slip
    /// limit from thickness times the admitted small elastic shear strain.
    pub fn new(stiffness_n_m: f64, yield_force_n: f64, plastic_resistance_n_s_m: f64,
        temperature_kelvin: f64, maximum_elastic_slip_m: f64) -> Result<Self, TriboError> {
        positive_finite(stiffness_n_m, "layer_stiffness_n_m")?;
        nonnegative_finite(yield_force_n, "layer_yield_force_n")?;
        positive_finite(plastic_resistance_n_s_m, "layer_plastic_resistance_n_s_m")?;
        positive_finite(temperature_kelvin, "layer_temperature_kelvin")?;
        positive_finite(maximum_elastic_slip_m, "layer_maximum_elastic_slip_m")?;
        Ok(Self { stiffness_n_m, yield_force_n, plastic_resistance_n_s_m,
            temperature_kelvin, maximum_elastic_slip_m })
    }
    #[must_use]
    pub const fn stiffness_n_m(&self) -> f64 { self.stiffness_n_m }
    #[must_use]
    pub const fn yield_force_n(&self) -> f64 { self.yield_force_n }
    #[must_use]
    pub const fn plastic_resistance_n_s_m(&self) -> f64 { self.plastic_resistance_n_s_m }
    #[must_use]
    pub const fn temperature_kelvin(&self) -> f64 { self.temperature_kelvin }
    #[must_use]
    pub const fn maximum_elastic_slip_m(&self) -> f64 { self.maximum_elastic_slip_m }

    /// Continue the SAME retained strain while changing only dissipative
    /// coefficients and their temperature. The elastic potential and its domain
    /// must be unchanged: this operation supplies no thermoelastic work.
    ///
    /// Material identity, source-claim continuity, geometry and temperature
    /// evolution belong to the coupled owner. This numerical leaf never selects
    /// a replacement material. Neither signed plastic slip nor reversible
    /// elastic slip is reset, including on a failed continuation.
    pub fn continue_state(&self, previous: &ThirdBodyLayerState)
        -> Result<ThirdBodyLayerState, TriboError> {
        if self.stiffness_n_m.to_bits() != previous.owner.stiffness_n_m.to_bits()
            || self.maximum_elastic_slip_m.to_bits()
                != previous.owner.maximum_elastic_slip_m.to_bits() {
            return Err(TriboError::InvalidInput {
                field: "layer_continuation_requires_unchanged_elastic_potential",
            });
        }
        Ok(self.continue_state_with_work(previous)?.state)
    }

    /// Continue the same retained strain with explicit coefficient work.
    /// A thermal owner must supply an admitted free-energy temperature law and
    /// close `delta(psi + T s)` with its caloric update. This numerical transfer
    /// alone supplies neither heat nor source authority. Geometry/tangent-frame
    /// changes remain outside this operation; the strain limit is unchanged.
    pub fn continue_state_with_work(&self, previous: &ThirdBodyLayerState)
        -> Result<ThirdBodyLayerContinuation, TriboError> {
        if self.maximum_elastic_slip_m.to_bits()
            != previous.owner.maximum_elastic_slip_m.to_bits() {
            return Err(TriboError::InvalidInput {
                field: "layer_continuation_requires_unchanged_elastic_domain",
            });
        }
        self.admit_elastic_slip(previous.elastic_slip_m)?;
        for value in previous.plastic_slip_m {
            finite(value, "layer_continued_plastic_slip_m")?;
        }
        let e = det::hypot(previous.elastic_slip_m[0], previous.elastic_slip_m[1]);
        let coefficient_free_energy_change_j = if self.stiffness_n_m
            == previous.owner.stiffness_n_m { 0.0 } else {
            checked_mul(checked_mul(0.5 * (self.stiffness_n_m
                - previous.owner.stiffness_n_m), e, "layer_coefficient_work_j")?,
                e, "layer_coefficient_work_j")?
        };
        let state = ThirdBodyLayerState {
            elastic_slip_m: previous.elastic_slip_m,
            plastic_slip_m: previous.plastic_slip_m,
            owner: self.clone(),
        };
        Ok(ThirdBodyLayerContinuation { state, coefficient_free_energy_change_j })
    }

    /// `s = -d(psi)/dT` and `U = psi + T s` for the retained affine spring.
    /// The caller owns admission of the complete affine temperature interval;
    /// an arbitrary supplied slope never becomes a measured material claim.
    pub fn thermodynamics(&self, state: &ThirdBodyLayerState,
        stiffness_slope_n_m_k: f64) -> Result<ThirdBodyLayerThermodynamics, TriboError> {
        if state.owner != *self {
            return Err(TriboError::InvalidInput { field: "layer_checkpoint_owner" });
        }
        finite(stiffness_slope_n_m_k, "layer_stiffness_temperature_derivative")?;
        self.admit_elastic_slip(state.elastic_slip_m)?;
        let e = det::hypot(state.elastic_slip_m[0], state.elastic_slip_m[1]);
        let free_energy_j = self.stored_energy_j(state);
        nonnegative_finite(free_energy_j, "layer_free_energy_j")?;
        let entropy_j_k = if stiffness_slope_n_m_k == 0.0 { 0.0 } else {
            checked_mul(checked_mul(-0.5 * stiffness_slope_n_m_k, e,
                "layer_mechanical_entropy_j_k")?, e, "layer_mechanical_entropy_j_k")?
        };
        let internal_energy_j = checked_add(free_energy_j,
            checked_mul(self.temperature_kelvin, entropy_j_k,
                "layer_entropy_energy_j")?, "layer_internal_energy_j")?;
        Ok(ThirdBodyLayerThermodynamics { free_energy_j, entropy_j_k, internal_energy_j })
    }

    fn admit_elastic_slip(&self, slip: [f64; 2]) -> Result<(), TriboError> {
        for value in slip { finite(value, "layer_elastic_slip_m")?; }
        let magnitude = det::hypot(slip[0], slip[1]);
        finite(magnitude, "layer_elastic_slip_norm_m")?;
        if magnitude > self.maximum_elastic_slip_m {
            return Err(TriboError::OutsideApplicability { field: "layer_elastic_slip_m",
                value: magnitude, minimum: 0.0, maximum: self.maximum_elastic_slip_m });
        }
        Ok(())
    }
}

impl crate::tangential::transport::ObjectiveTangentialLaw for ThirdBodyLayerContact {
    fn transport_state(&self, state: &Self::State,
        rotation: crate::tangential::transport::TangentialHistoryRotation)
        -> Result<Self::State, TriboError> {
        if state.owner != *self {
            return Err(TriboError::InvalidInput { field: "layer_checkpoint_owner" });
        }
        self.admit_elastic_slip(state.elastic_slip_m)?;
        let elastic_slip_m = rotation.apply(state.elastic_slip_m)?;
        let plastic_slip_m = rotation.apply(state.plastic_slip_m)?;
        self.admit_elastic_slip(elastic_slip_m)?;
        Ok(ThirdBodyLayerState { elastic_slip_m, plastic_slip_m, owner: state.owner.clone() })
    }
}

impl TangentialLaw for ThirdBodyLayerContact {
    type State = ThirdBodyLayerState;
    fn model_id(&self) -> &str { "retained-isothermal-elastic-viscoplastic-shear-layer-v1" }
    fn initial_state(&self) -> Self::State {
        ThirdBodyLayerState { elastic_slip_m: [0.0; 2], plastic_slip_m: [0.0; 2],
            owner: self.clone() }
    }
    fn stored_energy_j(&self, state: &Self::State) -> f64 {
        // States are constructed only by this law. A foreign checkpoint is not
        // a physical material transfer and must not produce a finite ledger.
        if state.owner != *self { return f64::NAN; }
        let e = det::hypot(state.elastic_slip_m[0], state.elastic_slip_m[1]);
        (0.5 * self.stiffness_n_m * e) * e
    }
    fn advance(&self, state: &Self::State, input: &TangentialStepInput)
        -> Result<TangentialStep<Self::State>, TriboError> {
        input.solve.validate()?;
        if state.owner != *self {
            return Err(TriboError::InvalidInput { field: "layer_checkpoint_owner" });
        }
        self.admit_elastic_slip(state.elastic_slip_m)?;
        for value in state.plastic_slip_m { finite(value, "layer_plastic_slip_m")?; }
        // No separation or reattachment rule is implied by this law.
        positive_finite(input.normal_n, "layer_requires_retained_compressive_contact")?;
        positive_finite(input.dt_s, "dt_s")?;
        nonnegative_finite(input.compliance_m_per_n, "compliance_m_per_n")?;
        let temperature = input.contact_temperature_k
            .ok_or(TriboError::InvalidInput { field: "layer_temperature_kelvin" })?;
        if temperature.to_bits() != self.temperature_kelvin.to_bits() {
            return Err(TriboError::OutsideApplicability { field: "layer_temperature_kelvin",
                value: temperature, minimum: self.temperature_kelvin, maximum: self.temperature_kelvin });
        }
        if let Some(pressure) = input.nominal_pressure_pa {
            nonnegative_finite(pressure, "contact_pressure_pa")?;
        }
        for value in input.slip_increment_m { finite(value, "slip_increment_m")?; }

        let half_k = 0.5 * self.stiffness_n_m;
        positive_finite(half_k, "layer_half_stiffness")?;
        let flow = input.dt_s / self.plastic_resistance_n_s_m;
        positive_finite(flow, "layer_plastic_step_compliance")?;
        let a = checked_add(1.0, checked_mul(half_k, input.compliance_m_per_n,
            "layer_load_line")?, "layer_load_line")?;
        let b = checked_mul(half_k, flow, "layer_return_coefficient")?;
        positive_finite(b, "layer_return_coefficient")?;
        let mut trial = [0.0; 2];
        for (i, value) in trial.iter_mut().enumerate() {
            *value = checked_add(checked_mul(self.stiffness_n_m, state.elastic_slip_m[i],
                "layer_trial_force_n")?, checked_mul(half_k, input.slip_increment_m[i],
                "layer_trial_force_n")?, "layer_trial_force_n")?;
        }
        let norm = det::hypot(trial[0], trial[1]);
        finite(norm, "layer_trial_force_norm_n")?;
        let trial_magnitude = norm / a;
        finite(trial_magnitude, "layer_elastic_trial_force_n")?;
        let yielding = trial_magnitude > self.yield_force_n;
        let magnitude = if yielding {
            // Equivalent to (norm+b*Y)/(a+b), without a possibly enormous b*Y.
            let excess = norm - a * self.yield_force_n;
            checked_add(self.yield_force_n, excess / checked_add(a, b, "layer_return_denominator")?,
                "layer_return_force_n")?
        } else { trial_magnitude };
        finite(magnitude, "layer_return_force_n")?;
        if magnitude > input.solve.maximum_force_n {
            return Err(TriboError::OutsideApplicability { field: "maximum_force_n",
                value: magnitude, minimum: 0.0, maximum: input.solve.maximum_force_n });
        }
        let direction = if norm == 0.0 { [0.0; 2] } else { [trial[0]/norm, trial[1]/norm] };
        let plastic_distance = if yielding {
            checked_mul(flow, magnitude - self.yield_force_n, "layer_plastic_increment_m")?
        } else { 0.0 };
        nonnegative_finite(plastic_distance, "layer_plastic_increment_m")?;
        let mut traction_n = [0.0; 2];
        let mut slip_increment_m = [0.0; 2];
        let mut candidate = state.clone();
        let mut force_residual = 0.0_f64;
        for i in 0..2 {
            traction_n[i] = checked_mul(magnitude, direction[i], "layer_traction_n")?;
            slip_increment_m[i] = input.slip_increment_m[i]
                - checked_mul(input.compliance_m_per_n, traction_n[i], "layer_body_response_m")?;
            finite(slip_increment_m[i], "layer_resolved_slip_m")?;
            let plastic = checked_mul(plastic_distance, direction[i], "layer_plastic_increment_m")?;
            candidate.elastic_slip_m[i] = state.elastic_slip_m[i] + slip_increment_m[i] - plastic;
            finite(candidate.elastic_slip_m[i], "layer_candidate_elastic_slip_m")?;
            candidate.plastic_slip_m[i] = checked_add(state.plastic_slip_m[i], plastic,
                "layer_accumulated_plastic_slip_m")?;
            let midpoint_force = half_k * (state.elastic_slip_m[i] + candidate.elastic_slip_m[i]);
            finite(midpoint_force, "layer_midpoint_force_n")?;
            force_residual = force_residual.max((traction_n[i] - midpoint_force).abs());
        }
        self.admit_elastic_slip(candidate.elastic_slip_m)?;
        let tolerance = input.solve.tolerance(magnitude, magnitude)?;
        if !force_residual.is_finite() || force_residual > tolerance {
            return Err(TriboError::TangentialSolve {
                residual_n: force_residual, tolerance_n: tolerance, iterations: 0 });
        }
        let old_energy = self.stored_energy_j(state);
        let new_energy = self.stored_energy_j(&candidate);
        nonnegative_finite(old_energy, "layer_old_energy_j")?;
        nonnegative_finite(new_energy, "layer_new_energy_j")?;
        let stored_change_j = new_energy - old_energy;
        let yield_work = checked_mul(self.yield_force_n, plastic_distance, "layer_yield_work_j")?;
        // This is actual plastic work, not a numerical energy deficit. Avoid
        // c/dt overflow in the zero-flow elastic branch.
        let viscous_work = if plastic_distance == 0.0 { 0.0 } else {
            checked_mul(plastic_distance / flow, plastic_distance, "layer_viscous_work_j")?
        };
        let dissipated_j = checked_add(yield_work, viscous_work, "layer_dissipated_j")?;
        nonnegative_finite(dissipated_j, "layer_dissipated_j")?;
        let work = checked_add(checked_mul(traction_n[0], slip_increment_m[0], "layer_work_j")?,
            checked_mul(traction_n[1], slip_increment_m[1], "layer_work_j")?, "layer_work_j")?;
        let closure_residual_j = stored_change_j + dissipated_j - work;
        finite(closure_residual_j, "layer_closure_residual_j")?;
        Ok(TangentialStep { traction_n, slip_increment_m, state: candidate,
            regime: if yielding { FrictionRegime::Sliding } else { FrictionRegime::Sticking },
            stored_change_j, dissipated_j, closure_residual_j,
            compatibility_residual_m: [0.0; 2], iterations: 0 })
    }
}
