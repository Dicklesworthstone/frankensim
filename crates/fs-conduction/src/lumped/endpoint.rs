//! Original lumped enthalpy and boundary powers for coupled thermal sources.
//!
//! Preparation retains the original geometry, transport sources, caloric chart
//! and Biot admission. A coupled owner may vary the endpoint ambient temperature
//! and integrated internal heat while solving its other stores. A driven step
//! delegates its implicit source feedback to the original bounded enthalpy
//! solver; read-only samples do not themselves accept a trajectory.

use super::*;
use fs_material::phase::UniformEnthalpyStepInput;

/// One prepared physical interval on the original uniform-body owner.
pub struct LumpedEnthalpyEndpoint<'a> {
    body: LumpedEnthalpyBody<'a>,
    config: LumpedEnthalpyMarchConfig,
    initial: EquilibriumPhaseState,
    maximum_thermal_residual_j: f64,
}

/// Original endpoint calorics and integrated boundary powers, in SI units.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LumpedEnthalpyEndpointSample {
    pub specific_enthalpy_j_kg: f64,
    pub temperature_k: f64,
    pub stored_energy_change_j: f64,
    pub convection_heat_j: f64,
    pub radiation_heat_j: f64,
    /// Actual storage minus the original endpoint net power times duration.
    pub residual_j: f64,
}

/// One implicit proposal with an additional caloric-state-dependent source.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DrivenLumpedEnthalpyStep {
    /// Original caloric state; latent heat remains in its enthalpy chart.
    pub state: EquilibriumPhaseState,
    /// Native storage, convection, radiation and complete source residual.
    pub endpoint: LumpedEnthalpyEndpointSample,
    /// Additional endpoint source power times the physical duration [J].
    pub additional_internal_heat_j: f64,
    /// Source calls, bounded by native bisection iterations plus two endpoints.
    pub source_evaluations: usize,
}

impl LumpedThermalEnvironment {
    /// Prepare repeated coupled residual evaluations on the SAME native wall.
    ///
    /// Geometry, source identity, radiation enclosure, constant convection
    /// coefficient and the original whole-chart Biot ceiling are admitted once.
    /// Endpoint ambient temperature and net integrated internal heat are supplied
    /// to each evaluation. Use the original `advance` for its accepted trajectory;
    /// coupled owners must recheck other stores against that returned wall.
    pub fn prepare_endpoint<'a>(
        &self,
        input: UniformEnthalpyStepInput<'a>,
    ) -> Result<LumpedEnthalpyEndpoint<'a>, ConductionError> {
        if input.initial.phase_curve_identity() != input.curve.identity()
            || !(input.duration_s.is_finite() && input.duration_s > 0.0)
            || !(input.volume_m3.is_finite() && input.volume_m3 > 0.0)
            || !(self.maximum_thermal_residual_j.is_finite()
                && self.maximum_thermal_residual_j >= 0.0)
        {
            return Err(lumped_error(
                "coupled uniform body",
                "thermal proposal needs its initial chart, positive duration/volume and a finite nonnegative energy-residual budget".to_owned(),
            ));
        }
        let body = LumpedEnthalpyBody::try_new_with_transport(
            "coupled uniform body",
            input.mass_kg,
            input.surface_area_m2,
            self.convection_w_per_m2_k,
            input.volume_m3 / input.surface_area_m2,
            self.transport.clone(),
            input.curve,
        )?;
        let config = LumpedEnthalpyMarchConfig {
            initial_specific_enthalpy_j_kg: input.initial.specific_enthalpy_j_kg(),
            ambient_temperature_k: self.temperature_k,
            radiation_temperature_k: self.radiation_temperature_k,
            internal_power_w: input.internal_heat_j / input.duration_s,
            duration_s: input.duration_s,
            maximum_step_s: input.duration_s,
            maximum_steps: 1,
            enthalpy_tolerance_j_kg: self.enthalpy_tolerance_j_kg,
        };
        if input.internal_heat_j != 0.0 && config.internal_power_w == 0.0 {
            return Err(lumped_error(
                &body.name,
                "nonzero prepared source underflows to zero power".to_owned(),
            ));
        }
        validate_enthalpy_config(&body, config)?;
        let maximum_biot = body.maximum_biot(config.radiation_temperature_k);
        let gate = BiotGate::corpus_default();
        if !maximum_biot.is_finite() || maximum_biot > gate.ceiling() {
            return Err(lumped_error(
                &body.name,
                format!(
                    "convection-plus-radiation Biot {maximum_biot} exceeds ceiling {}; escalate to spatial transient enthalpy transport",
                    gate.ceiling()
                ),
            ));
        }
        phase_state(&body, config.initial_specific_enthalpy_j_kg)?;
        Ok(LumpedEnthalpyEndpoint {
            body,
            config,
            initial: input.initial,
            maximum_thermal_residual_j: self.maximum_thermal_residual_j,
        })
    }
}

impl LumpedEnthalpyEndpoint<'_> {
    /// Original whole-surface convection/radiation bound used at preparation.
    /// Coupled owners must add independently occupied heat paths to this bound.
    #[must_use]
    pub fn maximum_biot(&self) -> f64 {
        self.body.maximum_biot(self.config.radiation_temperature_k)
    }

    /// Sufficient implicit step bound for native radiation and an extra source.
    ///
    /// The supplied bound is the maximum positive source-power slope with
    /// respect to TOTAL body energy [W/J], throughout the original chart.
    /// Combining it with the native radiation bound retains a residual slope
    /// of at least half the stored-energy slope. Infinity admits every finite
    /// duration; this bound alone does not certify temporal accuracy.
    pub fn maximum_step_s_for_internal_power(
        &self,
        maximum_positive_feedback_w_per_j: f64,
        mut checkpoint: impl FnMut() -> Result<(), ConductionError>,
    ) -> Result<f64, ConductionError> {
        checkpoint()?;
        if !maximum_positive_feedback_w_per_j.is_finite() || maximum_positive_feedback_w_per_j < 0.0
        {
            return Err(lumped_error(
                &self.body.name,
                "source feedback bound must be finite and nonnegative [W/J]".to_owned(),
            ));
        }
        let native = monotone_enthalpy_step_limit_s(
            &self.body,
            self.config.radiation_temperature_k,
            &mut checkpoint,
        )?;
        let limit = if maximum_positive_feedback_w_per_j == 0.0 {
            native
        } else {
            let source = (0.5 / maximum_positive_feedback_w_per_j).next_down();
            // Harmonic combination avoids overflowing twice the source slope.
            let smaller = native.min(source);
            let larger = native.max(source);
            if larger.is_infinite() {
                smaller
            } else {
                let denominator = (1.0 + (smaller / larger).next_up()).next_up();
                (smaller / denominator).next_down()
            }
        };
        if limit <= 0.0 || limit.is_nan() {
            return Err(lumped_error(
                &self.body.name,
                "combined thermal/source step bound is not representable".to_owned(),
            ));
        }
        checkpoint()?;
        Ok(limit)
    }

    /// Solve the native backward-Euler balance with an additional source [W].
    ///
    /// `source(h, temperature)` must return the same signed power at the same
    /// caloric state and must not mutate an external reservoir. Its positive
    /// slope against total body energy must obey the supplied whole-chart
    /// bound. Passing h retains the phase coordinate on latent plateaus.
    /// Prepared internal power remains included; source and boundary powers
    /// use the SAME endpoint. No state is published or mutated by this proposal.
    pub fn solve_with_internal_power(
        &self,
        maximum_positive_feedback_w_per_j: f64,
        mut source: impl FnMut(f64, f64) -> Result<f64, ConductionError>,
        mut checkpoint: impl FnMut() -> Result<(), ConductionError>,
    ) -> Result<DrivenLumpedEnthalpyStep, ConductionError> {
        let limit = self.maximum_step_s_for_internal_power(
            maximum_positive_feedback_w_per_j,
            &mut checkpoint,
        )?;
        let dt = self.config.duration_s;
        if dt > limit {
            return Err(lumped_error(&self.body.name,
                "combined source/radiative feedback exceeds the monotone implicit step bound; reduce the physical interval".to_owned()));
        }
        let energy_tolerance_j = self
            .maximum_thermal_residual_j
            .min(self.body.mass_kg * self.config.enthalpy_tolerance_j_kg);
        let mut source_evaluations = 0;
        // Cache the sign endpoints, including adjacent-float termination,
        // without allocating or evaluating the accepted source a second time.
        let mut negative = None;
        let mut positive = None;
        let h = bisect_enthalpy_step(
            &self.body,
            self.config.initial_specific_enthalpy_j_kg,
            self.config.enthalpy_tolerance_j_kg,
            energy_tolerance_j,
            |h| {
                checkpoint()?;
                let values = self
                    .body
                    .phase_curve
                    .values_at_specific_enthalpy(h)
                    .map_err(|error| lumped_error(&self.body.name, error.to_string()))?;
                source_evaluations += 1;
                let power = source(h, values.temperature_k)?;
                checkpoint()?;
                let heat = dt * power;
                let total_power = self.config.internal_power_w + power;
                if !power.is_finite()
                    || !heat.is_finite()
                    || !total_power.is_finite()
                    || (power != 0.0 && heat == 0.0)
                {
                    return Err(lumped_error(
                        &self.body.name,
                        "driven source power or heat is not representable".to_owned(),
                    ));
                }
                let sample = self.evaluate_native_power(
                    h,
                    values.temperature_k,
                    self.config.ambient_temperature_k,
                    total_power,
                )?;
                if sample.residual_j <= 0.0 {
                    negative = Some((sample, heat));
                } else {
                    positive = Some((sample, heat));
                }
                Ok(sample.residual_j)
            },
        )?;
        let (endpoint, additional_internal_heat_j) = [negative, positive]
            .into_iter()
            .flatten()
            .find(|(sample, _)| sample.specific_enthalpy_j_kg == h)
            .ok_or_else(|| {
                lumped_error(
                    &self.body.name,
                    "native enthalpy solution has no retained source sample".to_owned(),
                )
            })?;
        if endpoint.residual_j.abs() > energy_tolerance_j {
            return Err(lumped_error(
                &self.body.name,
                "driven endpoint does not close its original energy budget".to_owned(),
            ));
        }
        let state = if h == self.initial.specific_enthalpy_j_kg() {
            self.initial
        } else {
            phase_state(&self.body, h)?
        };
        checkpoint()?;
        Ok(DrivenLumpedEnthalpyStep {
            state,
            endpoint,
            additional_internal_heat_j,
            source_evaluations,
        })
    }

    /// Evaluate native backward-Euler storage and source/boundary powers.
    ///
    /// `storage_increment_j` is a numerical coordinate; the residual uses its
    /// actually represented h increment. `internal_heat_j` is the COMPLETE
    /// integrated internal source for this trial, replacing the preparation
    /// value. It may include a coupled heat debit, whose separate balance is
    /// owned by the caller. No latent heat or finite-reservoir heat is inferred.
    pub fn evaluate(
        &self,
        storage_increment_j: f64,
        ambient_temperature_k: f64,
        internal_heat_j: f64,
    ) -> Result<LumpedEnthalpyEndpointSample, ConductionError> {
        if !storage_increment_j.is_finite() {
            return Err(lumped_error(
                &self.body.name,
                "nonfinite endpoint heat or invalid ambient temperature".to_owned(),
            ));
        }
        let h0 = self.config.initial_specific_enthalpy_j_kg;
        let h = h0 + storage_increment_j / self.body.mass_kg;
        if !h.is_finite() || (storage_increment_j != 0.0 && h == h0) {
            return Err(lumped_error(
                &self.body.name,
                "endpoint energy cannot change its original enthalpy coordinate".to_owned(),
            ));
        }
        self.evaluate_specific_enthalpy(h, ambient_temperature_k, internal_heat_j)
    }

    /// Evaluate the same native powers at an already represented caloric state.
    /// This avoids converting an accepted h to joules and back when composing
    /// another thermal owner. It does not accept or advance that state.
    pub fn evaluate_specific_enthalpy(
        &self,
        specific_enthalpy_j_kg: f64,
        ambient_temperature_k: f64,
        internal_heat_j: f64,
    ) -> Result<LumpedEnthalpyEndpointSample, ConductionError> {
        if !specific_enthalpy_j_kg.is_finite()
            || !internal_heat_j.is_finite()
            || !ambient_temperature_k.is_finite()
            || ambient_temperature_k <= 0.0
        {
            return Err(lumped_error(
                &self.body.name,
                "nonfinite endpoint heat or invalid ambient temperature".to_owned(),
            ));
        }
        let h = specific_enthalpy_j_kg;
        let values = self
            .body
            .phase_curve
            .values_at_specific_enthalpy(h)
            .map_err(|error| lumped_error(&self.body.name, error.to_string()))?;
        let internal_power_w = internal_heat_j / self.config.duration_s;
        if internal_heat_j != 0.0 && internal_power_w == 0.0 {
            return Err(lumped_error(
                &self.body.name,
                "nonzero endpoint source underflows to zero power".to_owned(),
            ));
        }
        self.evaluate_native_power(
            h,
            values.temperature_k,
            ambient_temperature_k,
            internal_power_w,
        )
    }

    fn evaluate_native_power(
        &self,
        h: f64,
        temperature_k: f64,
        ambient_temperature_k: f64,
        internal_power_w: f64,
    ) -> Result<LumpedEnthalpyEndpointSample, ConductionError> {
        let h0 = self.config.initial_specific_enthalpy_j_kg;
        let mut config = self.config;
        config.ambient_temperature_k = ambient_temperature_k;
        config.internal_power_w = internal_power_w;
        let (convection, radiation, net) = enthalpy_power(&self.body, config, temperature_k)?;
        let stored_energy_change_j = self.body.mass_kg * (h - h0);
        let convection_heat_j = config.duration_s * convection;
        let radiation_heat_j = config.duration_s * radiation;
        let residual_j = stored_energy_change_j - config.duration_s * net;
        if ![
            stored_energy_change_j,
            convection_heat_j,
            radiation_heat_j,
            residual_j,
        ]
        .into_iter()
        .all(f64::is_finite)
        {
            return Err(lumped_error(
                &self.body.name,
                "endpoint storage or boundary heat is unrepresentable".to_owned(),
            ));
        }
        Ok(LumpedEnthalpyEndpointSample {
            specific_enthalpy_j_kg: h,
            temperature_k,
            stored_energy_change_j,
            convection_heat_j,
            radiation_heat_j,
            residual_j,
        })
    }
}
