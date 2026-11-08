//! Native equilibrium response of the same mixed-flow contact closure.

use super::{Error, MixtureWallStreamStep, PreparedMixtureWallStream, phase_error, positive};
use crate::chemistry::ThermalConductance;
use fs_material::thermochemical_gas::nonideal::flash::equilibrium::acoustics::{
    PHASE_ACOUSTIC_PT_EVALUATIONS, PhaseEquilibriumAcousticControls, PhaseEquilibriumAcoustics,
};
use fs_qty::Temperature;

/// Original outlet and its response to wall temperature at fixed feed and pressure.
///
/// The native nine-PT owner admits a differentiable phase neighborhood, including
/// equilibrium phase redistribution. Its refinement evidence is numerical, not
/// an interval-certified derivative bound. A bubble/dew crossing, critical
/// singularity, unresolved neighborhood or different phase policy refuses.
#[derive(Clone, Debug, PartialEq)]
pub struct MixtureWallStreamResponse {
    /// Unchanged original PH heat, entropy and value-query work receipts.
    pub stream: MixtureWallStreamStep,
    /// Original response at the identical PH endpoint; absent only for zero G.
    pub outlet_response: Option<PhaseEquilibriumAcoustics>,
    /// d(actual stream heat power)/dT_wall at the converged contact [W/K].
    pub heat_power_temperature_derivative_w_k: f64,
    /// dT_out/dT_wall at fixed feed, flow, pressure and conductance.
    pub outlet_temperature_derivative: f64,
}
impl MixtureWallStreamResponse {
    /// Value queries plus the additional native response probes, excluding preparation.
    #[must_use]
    pub fn pt_evaluations(&self) -> usize {
        self.stream.pt_evaluations
            + self
                .outlet_response
                .as_ref()
                .map_or(0, PhaseEquilibriumAcoustics::pt_evaluations)
    }

    /// Native EOS work ceilings with the same scope as pt_evaluations().
    #[must_use]
    pub fn eos_evaluation_upper_bound(&self) -> usize {
        self.stream.eos_evaluation_upper_bound
            + self
                .outlet_response
                .as_ref()
                .map_or(0, PhaseEquilibriumAcoustics::eos_evaluation_upper_bound)
    }
}

impl PreparedMixtureWallStream<'_> {
    /// Solve the original contact and differentiate its equilibrium closure.
    ///
    /// For q = G (T_wall - T_out(q / mass_flow)),
    /// dq/dT_wall = G C / (G + C), with C = mass_flow Cp_eq at the
    /// original outlet. No frozen-phase or constant-Cp caloric path is created.
    /// The returned tangent uses the native response's controlled numerical Cp;
    /// the original enthalpy residual still decides thermal-step acceptance.
    ///
    /// The response policy must equal the prepared PH phase policy exactly.
    /// Nine PT calls are reserved inside the supplied shared allowance before
    /// solving the value. Zero conductance skips them; zero thermal drive at
    /// nonzero conductance still needs the nonzero equilibrium response.
    pub fn evaluate_response(
        &self,
        wall_temperature: Temperature,
        conductance: ThermalConductance,
        response_controls: PhaseEquilibriumAcousticControls,
        remaining_duty_evaluations: usize,
        remaining_pt_evaluations: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MixtureWallStreamResponse, Error> {
        super::poll(&mut cancelled)?;
        response_controls
            .validate()
            .map_err(|e| phase_error(0, 0, e))?;
        if response_controls.phase != self.controls.equilibrium.phase_options() {
            return Err(Error::Invalid(
                "wall-stream response must retain the original PH phase policy",
            ));
        }
        let g = conductance.value();
        if !g.is_finite() || g < 0.0 {
            return Err(Error::Invalid(
                "nonnegative wall-stream conductance required",
            ));
        }
        let response_pt = if g == 0.0 {
            0
        } else {
            PHASE_ACOUSTIC_PT_EVALUATIONS
        };
        let maximum_pt = self
            .controls
            .maximum_pt_evaluations
            .min(remaining_pt_evaluations);
        if maximum_pt < response_pt {
            return Err(Error::WorkLimit {
                duty_evaluations: 0,
                pt_evaluations: 0,
            });
        }
        let stream = self.evaluate(
            wall_temperature,
            conductance,
            remaining_duty_evaluations,
            maximum_pt - response_pt,
            &mut cancelled,
        )?;
        let (outlet_response, derivative, outlet_derivative) = if g == 0.0 {
            (None, 0.0, 0.0)
        } else {
            let response = self
                .inlet
                .model
                .phase_equilibrium_acoustics(
                    Temperature::new(stream.outlet_equilibrium.state().temperature_k()),
                    self.inlet.pressure,
                    response_controls,
                    &mut cancelled,
                )
                .map_err(|e| phase_error(0, stream.pt_evaluations, e))?;
            if response.source() != stream.outlet_equilibrium.state() {
                return Err(Error::Invalid(
                    "wall-stream response changed the original PH endpoint",
                ));
            }
            let capacity = self.inlet.mass_flow.value() * response.isobaric_heat_capacity_j_kg_k();
            if !positive(capacity) {
                return Err(Error::Invalid(
                    "unrepresentable equilibrium heat-capacity flow rate",
                ));
            }
            // Avoid the overflowing product G*C and retain the smaller rate.
            let derivative = if g <= capacity {
                g / (1.0 + g / capacity)
            } else {
                capacity / (1.0 + capacity / g)
            };
            let outlet_derivative = derivative / capacity;
            if !positive(derivative) || !positive(outlet_derivative) {
                return Err(Error::Invalid(
                    "unrepresentable wall-stream equilibrium tangent",
                ));
            }
            (Some(response), derivative, outlet_derivative)
        };
        super::poll(&mut cancelled)?;
        Ok(MixtureWallStreamResponse {
            stream,
            outlet_response,
            heat_power_temperature_derivative_w_k: derivative,
            outlet_temperature_derivative: outlet_derivative,
        })
    }
}
