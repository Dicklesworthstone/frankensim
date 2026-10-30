//! Isothermal elastic-viscoplastic shear of a retained third-body layer.
//!
//! This is a homogeneous, pressure-independent shear layer, not a measured
//! rosin law. In two tangent axes, total slip is e + p and the continuous law is
//! F = k e + sum(k_j r_j), p_dot = (|F|-Y)_+ F/(c |F|), with optional
//! Maxwell branches `r_j_dot = e_dot-r_j/tau_j`. The held effort uses midpoint
//! strains: F = k(e_old+e_new)/2 + sum(k_j*r_j_mid). With du = du_free-C F, the radial
//! return has a closed form; no velocity regularization or fitted spectrum is
//! introduced. Work splits into the change of equilibrium/branch spring storage,
//! plastic work and each branch's `dt k_j |r_j_mid|^2/tau_j` heat. Relaxation
//! therefore occurs below yield and under held total slip, in the SAME law.
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

/// Explicit storage bound on material shear histories per physical contact.
pub const MAX_LAYER_RELAXATION_ARMS: usize = 64;

/// One Maxwell branch in parallel with the equilibrium shear spring. The
/// source owner derives stiffness as A*DeltaG/h; tau is a constitutive time [s].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThirdBodyLayerRelaxationArm {
    pub stiffness_n_m: f64,
    pub relaxation_time_s: f64,
}

/// Physical loss channels of one accepted constitutive interval [J].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThirdBodyLayerDissipation {
    pub plastic_j: f64,
    pub relaxation_j: f64,
}

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
    relaxation: Vec<ThirdBodyLayerRelaxationArm>,
}

/// Accepted elastic and plastic slip in the checkpoint's current tangent frame.
/// The private owner binding prevents reusing stored strain with another law.
#[derive(Clone, Debug, PartialEq)]
pub struct ThirdBodyLayerState {
    elastic_slip_m: [f64; 2],
    plastic_slip_m: [f64; 2],
    relaxation_slip_m: Vec<[f64; 2]>,
    owner: ThirdBodyLayerContact,
}
impl ThirdBodyLayerState {
    #[must_use]
    pub const fn elastic_slip_m(&self) -> [f64; 2] { self.elastic_slip_m }
    #[must_use]
    pub const fn plastic_slip_m(&self) -> [f64; 2] { self.plastic_slip_m }
    /// Spring displacement in each original Maxwell arm, in source order.
    pub fn relaxation_slip_m(&self) -> &[[f64; 2]] { &self.relaxation_slip_m }
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
            temperature_kelvin, maximum_elastic_slip_m, relaxation: Vec::new() })
    }
    /// Bind source-derived linear shear memory before creating its state. The
    /// original `stiffness_n_m` is the EQUILIBRIUM stiffness, not the sum of arms.
    /// Midpoint integration requires dt <= 2*tau for every arm, avoiding an
    /// unresolved alternating relaxation response. Refine dt independently.
    pub fn with_relaxation(mut self, arms: Vec<ThirdBodyLayerRelaxationArm>) -> Result<Self, TriboError> {
        if arms.len() > MAX_LAYER_RELAXATION_ARMS || !self.relaxation.is_empty() {
            return Err(TriboError::InvalidInput { field: "layer_relaxation_count_or_rebinding" });
        }
        let mut total = self.stiffness_n_m;
        for arm in &arms {
            positive_finite(arm.stiffness_n_m, "layer_relaxation_stiffness_n_m")?;
            positive_finite(arm.relaxation_time_s, "layer_relaxation_time_s")?;
            total = checked_add(total, arm.stiffness_n_m, "layer_instantaneous_stiffness_n_m")?;
        }
        self.relaxation = arms;
        Ok(self)
    }
    pub fn relaxation(&self) -> &[ThirdBodyLayerRelaxationArm] { &self.relaxation }
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
            || self.relaxation.len() != previous.owner.relaxation.len()
            || self.relaxation.iter().zip(&previous.owner.relaxation)
                .any(|(a, b)| a.stiffness_n_m.to_bits() != b.stiffness_n_m.to_bits())
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
        if self.relaxation.len() != previous.owner.relaxation.len()
            || previous.relaxation_slip_m.len() != self.relaxation.len()
            || self.maximum_elastic_slip_m.to_bits()
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
        let mut coefficient_free_energy_change_j = if self.stiffness_n_m
            == previous.owner.stiffness_n_m { 0.0 } else {
            checked_mul(checked_mul(0.5 * (self.stiffness_n_m
                - previous.owner.stiffness_n_m), e, "layer_coefficient_work_j")?,
                e, "layer_coefficient_work_j")?
        };
        for ((old, new), slip) in previous.owner.relaxation.iter()
            .zip(&self.relaxation).zip(&previous.relaxation_slip_m) {
            self.admit_elastic_slip(*slip)?;
            let r = det::hypot(slip[0], slip[1]);
            let change = checked_mul(checked_mul(0.5 * (new.stiffness_n_m-old.stiffness_n_m),
                r, "layer_memory_coefficient_work_j")?, r, "layer_memory_coefficient_work_j")?;
            coefficient_free_energy_change_j = checked_add(coefficient_free_energy_change_j,
                change, "layer_coefficient_work_j")?;
        }
        let state = ThirdBodyLayerState {
            elastic_slip_m: previous.elastic_slip_m,
            plastic_slip_m: previous.plastic_slip_m,
            relaxation_slip_m: previous.relaxation_slip_m.clone(),
            owner: self.clone(),
        };
        Ok(ThirdBodyLayerContinuation { state, coefficient_free_energy_change_j })
    }

    /// `s = -d(psi)/dT` and `U = psi + T s` for the retained affine spring.
    /// The caller owns admission of the complete affine temperature interval;
    /// an arbitrary supplied slope never becomes a measured material claim.
    /// Bound Maxwell stiffnesses must be temperature-independent on that
    /// interval; tau may vary because it stores no energy. Source consumers
    /// must enforce this requirement, rather than omit branch entropy.
    pub fn thermodynamics(&self, state: &ThirdBodyLayerState,
        stiffness_slope_n_m_k: f64) -> Result<ThirdBodyLayerThermodynamics, TriboError> {
        self.thermodynamics_with_relaxation(state, stiffness_slope_n_m_k,
            &[0.0; MAX_LAYER_RELAXATION_ARMS][..self.relaxation.len()])
    }

    /// Complete mechanical entropy/internal energy for affine equilibrium AND
    /// Maxwell stiffnesses: s=-k'|e|^2/2-sum(k_j'|r_j|^2/2), U=psi+T*s.
    /// Derivatives must cover every retained arm in source order. The caller
    /// admits their original material laws over the full thermal interval and
    /// includes coefficient work plus delta(T*s) in the finite caloric balance.
    /// No strain is reset and no relaxation loss is relabelled as reversible heat.
    pub fn thermodynamics_with_relaxation(&self, state: &ThirdBodyLayerState,
        stiffness_slope_n_m_k: f64, relaxation_slopes_n_m_k: &[f64])
        -> Result<ThirdBodyLayerThermodynamics, TriboError> {
        if state.owner != *self {
            return Err(TriboError::InvalidInput { field: "layer_checkpoint_owner" });
        }
        if relaxation_slopes_n_m_k.len() != self.relaxation.len()
            || state.relaxation_slip_m.len() != self.relaxation.len() {
            return Err(TriboError::InvalidInput { field: "layer_thermal_memory_derivative_count" });
        }
        finite(stiffness_slope_n_m_k, "layer_stiffness_temperature_derivative")?;
        self.admit_elastic_slip(state.elastic_slip_m)?;
        let e = det::hypot(state.elastic_slip_m[0], state.elastic_slip_m[1]);
        let free_energy_j = self.stored_energy_j(state);
        nonnegative_finite(free_energy_j, "layer_free_energy_j")?;
        let mut entropy_j_k = if stiffness_slope_n_m_k == 0.0 { 0.0 } else {
            checked_mul(checked_mul(-0.5 * stiffness_slope_n_m_k, e,
                "layer_mechanical_entropy_j_k")?, e, "layer_mechanical_entropy_j_k")?
        };
        for (&slope, slip) in relaxation_slopes_n_m_k.iter().zip(&state.relaxation_slip_m) {
            finite(slope, "layer_memory_temperature_derivative")?;
            self.admit_elastic_slip(*slip)?;
            if slope != 0.0 {
                let r = det::hypot(slip[0], slip[1]);
                let entropy = checked_mul(checked_mul(-0.5*slope, r,
                    "layer_memory_entropy_j_k")?, r, "layer_memory_entropy_j_k")?;
                entropy_j_k = checked_add(entropy_j_k, entropy, "layer_mechanical_entropy_j_k")?;
            }
        }
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

    /// Evaluate physical loss independently from the work/storage residual.
    /// Both states belong to the SAME held-coefficient law; any thermal
    /// continuation happens before this mechanical interval.
    pub fn step_dissipation(&self, old: &ThirdBodyLayerState, new: &ThirdBodyLayerState,
        dt_s: f64) -> Result<ThirdBodyLayerDissipation, TriboError> {
        positive_finite(dt_s, "layer_dissipation_dt_s")?;
        if old.owner != *self || new.owner != *self
            || old.relaxation_slip_m.len() != self.relaxation.len()
            || new.relaxation_slip_m.len() != self.relaxation.len() {
            return Err(TriboError::InvalidInput { field: "layer_checkpoint_owner" });
        }
        let distance = det::hypot(new.plastic_slip_m[0]-old.plastic_slip_m[0],
            new.plastic_slip_m[1]-old.plastic_slip_m[1]);
        let plastic_j = if distance == 0.0 { 0.0 } else {
            checked_add(checked_mul(self.yield_force_n, distance, "layer_yield_work_j")?,
                checked_mul(checked_mul(self.plastic_resistance_n_s_m, distance/dt_s,
                    "layer_plastic_viscous_force_n")?, distance, "layer_viscous_work_j")?,
                "layer_plastic_dissipation_j")?
        };
        let mut relaxation_j = 0.0;
        for ((arm, before), after) in self.relaxation.iter().zip(&old.relaxation_slip_m)
            .zip(&new.relaxation_slip_m) {
            let r = det::hypot(0.5*before[0]+0.5*after[0], 0.5*before[1]+0.5*after[1]);
            let heat = checked_mul(checked_mul(checked_mul(dt_s/arm.relaxation_time_s,
                arm.stiffness_n_m, "layer_relaxation_heat_j")?, r,
                "layer_relaxation_heat_j")?, r, "layer_relaxation_heat_j")?;
            relaxation_j = checked_add(relaxation_j, heat, "layer_relaxation_heat_j")?;
        }
        nonnegative_finite(plastic_j, "layer_plastic_dissipation_j")?;
        nonnegative_finite(relaxation_j, "layer_relaxation_heat_j")?;
        Ok(ThirdBodyLayerDissipation { plastic_j, relaxation_j })
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
        let relaxation_slip_m = state.relaxation_slip_m.iter()
            .map(|&slip| rotation.apply(slip)).collect::<Result<Vec<_>, _>>()?;
        self.admit_elastic_slip(elastic_slip_m)?;
        Ok(ThirdBodyLayerState { elastic_slip_m, plastic_slip_m, relaxation_slip_m,
            owner: state.owner.clone() })
    }
}

impl TangentialLaw for ThirdBodyLayerContact {
    type State = ThirdBodyLayerState;
    fn model_id(&self) -> &str {
        if self.relaxation.is_empty() { "retained-isothermal-elastic-viscoplastic-shear-layer-v1" }
        else { "retained-isothermal-viscoelastic-viscoplastic-shear-layer-v1" }
    }
    fn initial_state(&self) -> Self::State {
        ThirdBodyLayerState { elastic_slip_m: [0.0; 2], plastic_slip_m: [0.0; 2],
            relaxation_slip_m: vec![[0.0; 2]; self.relaxation.len()],
            owner: self.clone() }
    }
    fn stored_energy_j(&self, state: &Self::State) -> f64 {
        // States are constructed only by this law. A foreign checkpoint is not
        // a physical material transfer and must not produce a finite ledger.
        if state.owner != *self { return f64::NAN; }
        let e = det::hypot(state.elastic_slip_m[0], state.elastic_slip_m[1]);
        let mut energy = (0.5 * self.stiffness_n_m * e) * e;
        for (arm, slip) in self.relaxation.iter().zip(&state.relaxation_slip_m) {
            let r = det::hypot(slip[0], slip[1]);
            energy += (0.5 * arm.stiffness_n_m * r) * r;
        }
        energy
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

        let mut stiffness = self.stiffness_n_m;
        let mut history = state.elastic_slip_m.map(|e| self.stiffness_n_m*e);
        let mut weights = Vec::with_capacity(self.relaxation.len());
        for (arm, slip) in self.relaxation.iter().zip(&state.relaxation_slip_m) {
            self.admit_elastic_slip(*slip)?;
            let ratio = input.dt_s/arm.relaxation_time_s;
            if !ratio.is_finite() || ratio > 2.0 {
                return Err(TriboError::OutsideApplicability { field: "layer_dt_over_relaxation_time",
                    value: ratio, minimum: 0.0, maximum: 2.0 });
            }
            let weight = 1.0/(1.0+0.5*ratio);
            let effective = checked_mul(arm.stiffness_n_m, weight, "layer_effective_stiffness_n_m")?;
            stiffness = checked_add(stiffness, effective, "layer_effective_stiffness_n_m")?;
            for i in 0..2 { history[i] = checked_add(history[i],
                checked_mul(effective, slip[i], "layer_memory_force_n")?, "layer_history_force_n")?; }
            weights.push(weight);
        }
        let half_k = 0.5 * stiffness;
        positive_finite(half_k, "layer_half_stiffness")?;
        let flow = input.dt_s / self.plastic_resistance_n_s_m;
        positive_finite(flow, "layer_plastic_step_compliance")?;
        let a = checked_add(1.0, checked_mul(half_k, input.compliance_m_per_n,
            "layer_load_line")?, "layer_load_line")?;
        let b = checked_mul(half_k, flow, "layer_return_coefficient")?;
        positive_finite(b, "layer_return_coefficient")?;
        let mut trial = [0.0; 2];
        for (i, value) in trial.iter_mut().enumerate() {
            *value = checked_add(history[i], checked_mul(half_k, input.slip_increment_m[i],
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
            let delta = candidate.elastic_slip_m[i]-state.elastic_slip_m[i];
            let mut midpoint_force = 0.5*self.stiffness_n_m
                * (state.elastic_slip_m[i] + candidate.elastic_slip_m[i]);
            for (j, arm) in self.relaxation.iter().enumerate() {
                let mid = weights[j]*(state.relaxation_slip_m[j][i]+0.5*delta);
                candidate.relaxation_slip_m[j][i] = 2.0*mid-state.relaxation_slip_m[j][i];
                midpoint_force += arm.stiffness_n_m*mid;
            }
            finite(midpoint_force, "layer_midpoint_force_n")?;
            force_residual = force_residual.max((traction_n[i] - midpoint_force).abs());
        }
        self.admit_elastic_slip(candidate.elastic_slip_m)?;
        for slip in &candidate.relaxation_slip_m { self.admit_elastic_slip(*slip)?; }
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
        let relaxation_j = self.step_dissipation(state, &candidate, input.dt_s)?.relaxation_j;
        let dissipated_j = checked_add(checked_add(yield_work, viscous_work, "layer_plastic_dissipation_j")?,
            relaxation_j, "layer_dissipated_j")?;
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
