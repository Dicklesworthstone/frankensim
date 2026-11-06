//! Original lumped enthalpy and boundary powers without an inner thermal solve.
//!
//! Preparation retains the original geometry, transport sources, caloric chart
//! and Biot admission. A coupled owner may vary the endpoint ambient temperature
//! and integrated internal heat while solving its other stores. These samples
//! are residual evaluations, not accepted trajectories or substitute wall laws.

use super::*;
use fs_material::phase::UniformEnthalpyStepInput;

/// One prepared physical interval on the original uniform-body owner.
pub struct LumpedEnthalpyEndpoint<'a> {
    body: LumpedEnthalpyBody<'a>,
    config: LumpedEnthalpyMarchConfig,
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
        Ok(LumpedEnthalpyEndpoint { body, config })
    }
}

impl LumpedEnthalpyEndpoint<'_> {
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
        if !storage_increment_j.is_finite()
            || !internal_heat_j.is_finite()
            || !ambient_temperature_k.is_finite()
            || ambient_temperature_k <= 0.0
        {
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
        let values = self
            .body
            .phase_curve
            .values_at_specific_enthalpy(h)
            .map_err(|error| lumped_error(&self.body.name, error.to_string()))?;
        let mut config = self.config;
        config.ambient_temperature_k = ambient_temperature_k;
        config.internal_power_w = internal_heat_j / config.duration_s;
        if internal_heat_j != 0.0 && config.internal_power_w == 0.0 {
            return Err(lumped_error(
                &self.body.name,
                "nonzero endpoint source underflows to zero power".to_owned(),
            ));
        }
        let (convection, radiation, net) =
            enthalpy_power(&self.body, config, values.temperature_k)?;
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
            temperature_k: values.temperature_k,
            stored_energy_change_j,
            convection_heat_j,
            radiation_heat_j,
            residual_j,
        })
    }
}
