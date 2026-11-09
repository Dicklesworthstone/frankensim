//! Complete frequency curvature with three original-state solves on one LU.

use super::{
    ScatteringFactor, ScatteringFrequencyResponse, ScatteringModel,
    mixed::{add_i_state, observe, residual},
};
use crate::harmonic::{
    C64, HarmonicError, HarmonicFrequencyCurvatureResponse, HarmonicFrequencyResponse,
    HarmonicParameterDirection, Sum, finite, invalid, mixed as physical, poll,
};
use physical::{SourceJet, TangentState};

/// One original amplitude and its first two partial derivatives per rad/s.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ScatteringComplexFrequencyJet {
    pub value: C64,
    pub frequency: C64,
    pub second_frequency: C64,
}

/// Original signed watts, watts per rad/s and watts per (rad/s)^2.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ScatteringPowerFrequencyJet {
    pub value: f64,
    pub frequency: f64,
    pub second_frequency: f64,
}

/// Independent physical and wave-flux accounts at one selected original port.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ScatteringWaveFrequencyCurvature {
    pub original_port: usize,
    /// Positive output/input reference, held fixed in frequency.
    pub reference: f64,
    pub incident: ScatteringComplexFrequencyJet,
    pub outgoing: ScatteringComplexFrequencyJet,
    pub incident_power: ScatteringPowerFrequencyJet,
    pub outgoing_power: ScatteringPowerFrequencyJet,
    pub physical_power: ScatteringPowerFrequencyJet,
    /// Incident minus outgoing minus independently evaluated conjugate work.
    pub power_residual: ScatteringPowerFrequencyJet,
}

/// Original loaded response and complete first/second frequency derivatives.
///
/// Reference loading is boundary work, never internal material heat. Original
/// descriptor constraints, all material memory and physical source jets remain.
/// A nonlinear equilibrium image retains its signed incremental-work semantics.
#[derive(Clone, Debug, PartialEq)]
pub struct ScatteringFrequencyCurvatureResponse {
    pub physical: HarmonicFrequencyCurvatureResponse,
    /// Actual u=(2a-y)/Z at selected ports and ordinary u elsewhere.
    pub physical_inputs: Vec<ScatteringComplexFrequencyJet>,
    /// Same selected-port order as ScatteringModel::wave_ports().
    pub waves: Vec<ScatteringWaveFrequencyCurvature>,
    /// Original loss minus incident/outgoing flux and ordinary source work.
    pub wave_power_residual: ScatteringPowerFrequencyJet,
}

impl ScatteringFactor<'_> {
    /// Resolve a full loaded frequency curvature using exactly three solves.
    ///
    /// One already closed LU serves the primal, omega and omega-omega systems.
    /// Input slots are incident output waves at selected ports and ordinary
    /// conjugate inputs elsewhere. Optional first/second source derivatives
    /// preserve coherent amplitude and phase dispersion; None means zero.
    /// Fixed positive references supply feedback before every solve. Physical
    /// dispersive loading belongs in retained model states, not in Z(omega).
    ///
    /// The higher original unclosed equations are checked using the actual
    /// second input derivative. No dummy parameter solve, parameter direction,
    /// extra factor or finite frequency perturbation is introduced. The native
    /// curvature budget plus three scattering payload/work bounds and bounded
    /// O(n^2+n*ports) observation workspace are admitted before allocation.
    pub fn response_with_frequency_curvature(
        &self,
        input: &[C64],
        input_angular_derivative: Option<&[C64]>,
        input_angular_second_derivative: Option<&[C64]>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ScatteringFrequencyCurvatureResponse, HarmonicError> {
        poll(&mut cancelled)?;
        let model = &self.model.closed;
        for values in [
            Some(input),
            input_angular_derivative,
            input_angular_second_derivative,
        ]
        .into_iter()
        .flatten()
        {
            if values.len() != model.m || values.iter().any(|&value| !finite(value)) {
                return Err(invalid(
                    "complete finite original frequency source jets required",
                ));
            }
        }
        physical::admit_curvature(model)?;
        self.model.admit_curvature_observation()?;
        let frequency = self.response_with_frequency_derivative(
            input,
            input_angular_derivative,
            &mut cancelled,
        )?;
        let mut rhs = vec![C64::ZERO; model.n];
        let (second, second_input) = self.second_frequency_state(
            &frequency,
            input_angular_second_derivative,
            &mut rhs,
            &mut cancelled,
        )?;
        let actual_input = frequency.response.physical_inputs;
        let actual_frequency = frequency.physical_input_derivative;
        let base_waves = frequency.response.waves;
        let first_waves = frequency.wave_derivatives;
        let wave_residual = frequency.response.wave_power_residual_w;
        let wave_residual_first = frequency.wave_power_residual_derivative;
        let actual = SourceJet {
            input: &actual_input,
            frequency: Some(&actual_frequency),
            second_frequency: Some(&second_input),
            parameter: None,
            mixed: None,
        };
        let ordinary_frequency = HarmonicFrequencyResponse {
            response: frequency.response.physical,
            state_derivative: frequency.state_derivative,
            effort_derivative: frequency.effort_derivative,
            output_derivative: frequency.output_derivative,
            mean_stored_energy_derivative: frequency.mean_stored_energy_derivative,
            port_power_derivative: frequency.port_power_derivative,
            supplied_power_derivative: frequency.supplied_power_derivative,
            dissipated_power_derivative: frequency.dissipated_power_derivative,
            power_residual_derivative: frequency.power_residual_derivative,
            relative_derivative_residual: frequency.relative_derivative_residual,
        };
        let physical =
            physical::observe_curvature(model, actual, ordinary_frequency, second, &mut cancelled)?;
        let first = &physical.frequency_response;
        let mut supply_second = physical::physical_power_curvature_sum(
            model,
            actual,
            &first.response.effort,
            &first.effort_derivative,
            &physical.effort_second_derivative,
            |p| !self.model.selected[p],
            &mut cancelled,
        )?;
        let mut waves = Vec::with_capacity(base_waves.len());
        for (base, first_wave) in base_waves.into_iter().zip(first_waves) {
            poll(&mut cancelled)?;
            let p = base.original_port;
            let incident_second = input_angular_second_derivative.map_or(C64::ZERO, |a| a[p]);
            let incident = ScatteringComplexFrequencyJet {
                value: base.incident,
                frequency: first_wave.incident,
                second_frequency: incident_second,
            };
            let outgoing = ScatteringComplexFrequencyJet {
                value: base.outgoing,
                frequency: first_wave.outgoing,
                second_frequency: physical.output_second_derivative[p] - incident_second,
            };
            if !finite(outgoing.second_frequency) {
                return Err(invalid("outgoing wave curvature is not representable"));
            }
            let incident_second_power = wave_second_power(incident, base.reference)?;
            let outgoing_second_power = wave_second_power(outgoing, base.reference)?;
            let work_second = physical.port_power_second_derivative[p];
            let mut defect = Sum::default();
            defect
                .product(incident_second_power, 1.0)
                .and_then(|()| defect.product(outgoing_second_power, -1.0))
                .and_then(|()| defect.product(work_second, -1.0))
                .and_then(|()| supply_second.product(incident_second_power, 1.0))
                .and_then(|()| supply_second.product(outgoing_second_power, -1.0))
                .ok_or_else(|| invalid("wave power curvature balance is not representable"))?;
            waves.push(ScatteringWaveFrequencyCurvature {
                original_port: p,
                reference: base.reference,
                incident,
                outgoing,
                incident_power: power_jet(
                    base.incident_power_w,
                    first_wave.incident_power_w,
                    incident_second_power,
                ),
                outgoing_power: power_jet(
                    base.outgoing_power_w,
                    first_wave.outgoing_power_w,
                    outgoing_second_power,
                ),
                physical_power: power_jet(
                    base.physical_power_w,
                    first_wave.physical_power_w,
                    work_second,
                ),
                power_residual: power_jet(
                    base.power_residual_w,
                    first_wave.power_residual_w,
                    power_value(defect)?,
                ),
            });
        }
        let mut defect = Sum::default();
        defect
            .product(physical.dissipated_power_second_derivative, 1.0)
            .and_then(|()| defect.add_scaled(supply_second, -1.0))
            .ok_or_else(|| invalid("original loss minus wave curvature is not representable"))?;
        let physical_inputs = (0..model.m)
            .map(|p| ScatteringComplexFrequencyJet {
                value: actual_input[p],
                frequency: actual_frequency[p],
                second_frequency: second_input[p],
            })
            .collect();
        poll(&mut cancelled)?;
        Ok(ScatteringFrequencyCurvatureResponse {
            physical,
            physical_inputs,
            waves,
            wave_power_residual: power_jet(
                wave_residual,
                wave_residual_first,
                power_value(defect)?,
            ),
        })
    }

    /// Shared by the three-solve curvature and five-solve mixed owners.
    pub(super) fn second_frequency_state(
        &self,
        frequency: &ScatteringFrequencyResponse,
        source_second: Option<&[C64]>,
        rhs: &mut [C64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(TangentState, Vec<C64>), HarmonicError> {
        let model = &self.model.closed;
        let source = source_second.map_or_else(|| vec![C64::ZERO; model.m], <[C64]>::to_vec);
        let effective = self.effective_input(&source, cancelled)?;
        for (i, value) in rhs.iter_mut().enumerate() {
            poll(cancelled)?;
            let row = &model.g[i * model.m..(i + 1) * model.m];
            let mut terms = crate::harmonic::input_row_sums(row, |p| effective[p])?;
            if i < model.n_diff {
                add_i_state(&mut terms, frequency.state_derivative[i], 2.0)?;
            }
            *value = crate::harmonic::input_value(terms, C64::ZERO)?;
        }
        let (state, relative) = self.factor.solve_rhs(rhs, cancelled)?;
        let mut second = physical::recover_tangent(
            model,
            state,
            HarmonicParameterDirection::default(),
            &frequency.response.physical.state,
            &frequency.response.physical.effort,
            relative,
            cancelled,
        )?;
        let actual = self.restore_input_tangent(
            source_second,
            &second.output,
            &frequency.physical_input_derivative,
            None,
            cancelled,
        )?;
        second.relative_residual =
            second
                .relative_residual
                .max(residual::original_higher_residual(
                    self.model,
                    frequency.response.physical.omega_rad_s,
                    &second,
                    &actual,
                    &frequency.state_derivative,
                    &frequency.effort_derivative,
                    &frequency.physical_input_derivative,
                    &frequency.state_derivative,
                    2.0,
                    HarmonicParameterDirection::default(),
                    cancelled,
                )?);
        Ok((second, actual))
    }
}

impl ScatteringModel {
    fn admit_curvature_observation(&self) -> Result<(), HarmonicError> {
        let (n, m, k) = (
            self.closed.n as u128,
            self.closed.m as u128,
            self.ports.len() as u128,
        );
        let bytes = 3 * (8 * n * k + 256 * k + 64 * m) + 256 * (n + m + k);
        let work =
            3 * (2 * n * n * (k + 1) + n * m + k * k + 32 * m) + 24 * n * (n + m) + 96 * (m + k);
        if bytes > self.budget.maximum_payload_bytes as u128
            || work > self.budget.maximum_work as u128
        {
            return Err(invalid(
                "insufficient loaded curvature payload or work budget",
            ));
        }
        Ok(())
    }
}

fn wave_second_power(a: ScatteringComplexFrequencyJet, z: f64) -> Result<f64, HarmonicError> {
    Ok(observe::wave_power(
        super::ScatteringComplexJet {
            value: a.value,
            frequency: a.frequency,
            second_frequency: a.second_frequency,
            parameter: C64::ZERO,
            mixed: C64::ZERO,
        },
        z,
        0.0,
    )?
    .second_frequency)
}
fn power_jet(value: f64, frequency: f64, second_frequency: f64) -> ScatteringPowerFrequencyJet {
    ScatteringPowerFrequencyJet {
        value,
        frequency,
        second_frequency,
    }
}
fn power_value(value: Sum) -> Result<f64, HarmonicError> {
    value
        .value()
        .ok_or_else(|| invalid("wave power curvature sum is not representable"))
}
