//! Complete original wave-loaded frequency/parameter jets on one closed LU.
//!
//! The boundary reference is positive and constant in frequency. Its authored
//! parameter derivative changes both feedback and flux normalization. Original
//! material heat excludes this boundary exchange at every derivative order.

use super::{ScatteringFactor, ScatteringModel};
use crate::harmonic::{
    C64, HarmonicError, HarmonicFrequencyParameterResponse, HarmonicFrequencyResponse,
    HarmonicParameterDirection, Sum, finite, invalid, mixed as physical, parameter, poll,
};
use physical::{SourceJet, TangentState};

pub(super) mod observe;
pub(super) mod residual;
pub use observe::{ScatteringComplexJet, ScatteringPowerJet, ScatteringWaveFrequencyParameterJet};

/// Every original physical field and actual wave-boundary quantity through
/// value, omega, omega-omega, theta and omega-theta order.
///
/// The real parameter is authored in the original fixed state/descriptor
/// chart. Angular derivatives are per rad/s. Power derivatives are signed
/// slopes of watts; they are not additional losses or guarantees at a pole.
#[derive(Clone, Debug, PartialEq)]
pub struct ScatteringFrequencyParameterResponse {
    /// Full original states, efforts, outputs, storage and R/T heat jets.
    /// Source work uses the restored physical inputs, never the 2a/Z drive.
    pub physical: HarmonicFrequencyParameterResponse,
    /// Actual u=(2a-y)/Z on selected ports, ordinary u elsewhere.
    /// Original source-port order is retained at all five orders.
    pub physical_inputs: Vec<ScatteringComplexJet>,
    /// Selected ports in the authored wave_ports() order.
    pub waves: Vec<ScatteringWaveFrequencyParameterJet>,
    /// Independently observed original loss minus incident-outgoing flux
    /// and ordinary-port source work, at all five orders.
    pub wave_power_residual: ScatteringPowerJet,
}

impl ScatteringFactor<'_> {
    /// Resolve the full loaded frequency/parameter jet without differentiating
    /// an isolated impedance or refactorizing a perturbed model.
    ///
    /// Source slots are incident output waves on selected ports and ordinary
    /// physical inputs elsewhere. Optional slices independently author their
    /// omega, omega-omega, theta and omega-theta derivatives; None means zero.
    /// All slices contain every original port. The reference derivative is
    /// dZ/dtheta in wave_ports() order. Z has no frequency dependence; physical
    /// dispersive boundary laws belong in the original material/memory states.
    ///
    /// Five solves share this already closed factor: primal, omega, theta,
    /// omega-omega and omega-theta. Both occurrences of a changing G, the
    /// changing Q observation, the full source law and reference feedback are
    /// retained. Original descriptor constraints, memory and incremental T
    /// work remain. Every higher solve also checks the original unclosed rows.
    ///
    /// Native mixed-jet admission applies. Added scattering resources are six
    /// constructor payload bounds plus 512(n+m+k) bytes, and six constructor
    /// work scores plus 48n(n+m)+192(m+k), where k is the selected count.
    /// All inputs and budgets are checked before a response allocation. A
    /// cancellation or failed residual/arithmetic check returns no partial jet.
    #[allow(clippy::too_many_arguments)]
    pub fn response_with_frequency_parameter_derivative(
        &self,
        input: &[C64],
        direction: HarmonicParameterDirection<'_>,
        input_angular_derivative: Option<&[C64]>,
        input_angular_second_derivative: Option<&[C64]>,
        input_parameter_derivative: Option<&[C64]>,
        input_angular_parameter_derivative: Option<&[C64]>,
        reference_derivative: Option<&[f64]>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ScatteringFrequencyParameterResponse, HarmonicError> {
        poll(&mut cancelled)?;
        let model = &self.model.closed;
        let source = SourceJet {
            input,
            frequency: input_angular_derivative,
            second_frequency: input_angular_second_derivative,
            parameter: input_parameter_derivative,
            mixed: input_angular_parameter_derivative,
        };
        for values in [
            Some(source.input),
            source.frequency,
            source.second_frequency,
            source.parameter,
            source.mixed,
        ]
        .into_iter()
        .flatten()
        {
            if values.len() != model.m || values.iter().any(|&value| !finite(value)) {
                return Err(invalid("complete finite original source jets required"));
            }
        }
        if reference_derivative.is_some_and(|values| {
            values.len() != self.model.ports.len() || values.iter().any(|value| !value.is_finite())
        }) {
            return Err(invalid(
                "complete finite selected-reference derivative required",
            ));
        }
        direction.validate(model, &mut cancelled)?;
        physical::admit_mixed(model)?;
        self.model.admit_mixed_observation()?;
        let frequency =
            self.response_with_frequency_derivative(input, source.frequency, &mut cancelled)?;

        // The parameter solve shares the existing primal, including its actual
        // loaded input. The common observation owner restores original Q/G/R/T.
        let action = parameter::generator_action(
            model,
            &direction,
            &frequency.response.physical.state,
            &frequency.response.physical.effort,
            &mut cancelled,
        )?;
        let mut rhs = action.generator;
        super::parameter::add_parameter_source(
            self.model,
            direction,
            &frequency.response.physical.effort,
            &frequency.response.physical_inputs,
            &action.storage,
            source.parameter,
            reference_derivative,
            &mut rhs,
            None,
            &mut cancelled,
        )?;
        let (state, relative) = self.factor.solve_rhs(&rhs, &mut cancelled)?;
        let mut tangent = parameter::observe_state(
            model,
            &frequency.response.physical,
            direction,
            state,
            &action.storage,
            relative,
            &mut cancelled,
        )?;
        let parameter_input = self.restore_input_tangent(
            source.parameter,
            &tangent.output,
            &frequency.response.physical_inputs,
            reference_derivative,
            &mut cancelled,
        )?;
        tangent.relative_residual =
            tangent
                .relative_residual
                .max(self.model.original_parameter_residual(
                    super::parameter::ParameterResidualState {
                        laplace: C64::new(0.0, -frequency.response.physical.omega_rad_s),
                        laplace_derivative: C64::ZERO,
                        state: &frequency.response.physical.state,
                        effort: &frequency.response.physical.effort,
                        input: &frequency.response.physical_inputs,
                        state_derivative: &tangent.state,
                        effort_derivative: &tangent.effort,
                        input_derivative: &parameter_input,
                    },
                    direction,
                    &mut cancelled,
                )?);
        let parameter = parameter::finish_response(
            model,
            direction.ports,
            frequency.response.physical.clone(),
            tangent,
            &frequency.response.physical_inputs,
            Some(&parameter_input),
            &mut cancelled,
        )?;

        let (second, second_input) = self.second_frequency_state(
            &frequency,
            source.second_frequency,
            &mut rhs,
            &mut cancelled,
        )?;

        // Differentiate the complete loaded omega equation. Its coefficient
        // action and explicit source retain G_theta, Q_theta and Z_theta.
        let action = parameter::generator_action(
            model,
            &direction,
            &frequency.state_derivative,
            &frequency.effort_derivative,
            &mut cancelled,
        )?;
        rhs = action.generator;
        super::parameter::add_parameter_source(
            self.model,
            direction,
            &frequency.effort_derivative,
            &frequency.physical_input_derivative,
            &action.storage,
            source.mixed,
            reference_derivative,
            &mut rhs,
            None,
            &mut cancelled,
        )?;
        for (i, value) in rhs[..model.n_diff].iter_mut().enumerate() {
            poll(&mut cancelled)?;
            let mut terms = [Sum::default(); 2];
            add_i_state(&mut terms, parameter.state_derivative[i], 1.0)?;
            *value = crate::harmonic::input_value(terms, *value)?;
        }
        let (state, relative) = self.factor.solve_rhs(&rhs, &mut cancelled)?;
        let mut mixed = physical::recover_tangent(
            model,
            state,
            direction,
            &frequency.state_derivative,
            &frequency.effort_derivative,
            relative,
            &mut cancelled,
        )?;
        let mixed_input = self.restore_input_tangent(
            source.mixed,
            &mixed.output,
            &frequency.physical_input_derivative,
            reference_derivative,
            &mut cancelled,
        )?;
        mixed.relative_residual = mixed
            .relative_residual
            .max(residual::original_higher_residual(
                self.model,
                frequency.response.physical.omega_rad_s,
                &mixed,
                &mixed_input,
                &frequency.state_derivative,
                &frequency.effort_derivative,
                &frequency.physical_input_derivative,
                &parameter.state_derivative,
                1.0,
                direction,
                &mut cancelled,
            )?);

        let actual_input = frequency.response.physical_inputs;
        let actual_frequency = frequency.physical_input_derivative;
        let actual_source = SourceJet {
            input: &actual_input,
            frequency: Some(&actual_frequency),
            second_frequency: Some(&second_input),
            parameter: Some(&parameter_input),
            mixed: Some(&mixed_input),
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
        let curvature = physical::observe_curvature(
            model,
            actual_source,
            ordinary_frequency,
            second,
            &mut cancelled,
        )?;
        let physical = physical::observe_mixed(
            model,
            actual_source,
            direction,
            curvature,
            parameter,
            mixed,
            &mut cancelled,
        )?;
        let (waves, wave_power_residual) = observe::observe_waves(
            self.model,
            source,
            actual_source,
            direction,
            reference_derivative,
            &physical,
            &mut cancelled,
        )?;
        let physical_inputs = (0..model.m)
            .map(|p| ScatteringComplexJet::from_source(actual_source, p))
            .collect();
        poll(&mut cancelled)?;
        Ok(ScatteringFrequencyParameterResponse {
            physical,
            physical_inputs,
            waves,
            wave_power_residual,
        })
    }

    pub(super) fn restore_input_tangent(
        &self,
        source: Option<&[C64]>,
        output: &[C64],
        base_input: &[C64],
        reference_derivative: Option<&[f64]>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<C64>, HarmonicError> {
        let mut input =
            source.map_or_else(|| vec![C64::ZERO; self.model.closed.m], <[C64]>::to_vec);
        for (index, &(p, z)) in self.model.ports.iter().enumerate() {
            poll(cancelled)?;
            let dz = reference_derivative.map_or(0.0, |v| v[index]);
            let mut terms = [Sum::default(); 2];
            for (value, coefficient) in [
                (input[p], 2.0 / z),
                (output[p], -1.0 / z),
                (base_input[p], -dz / z),
            ] {
                terms[0]
                    .product(value.re, coefficient)
                    .and_then(|()| terms[1].product(value.im, coefficient))
                    .ok_or_else(|| invalid("actual loaded input jet is not representable"))?;
            }
            input[p] = crate::harmonic::input_value(terms, C64::ZERO)?;
        }
        Ok(input)
    }
}

impl ScatteringModel {
    pub(in crate::harmonic::scattering) fn admit_mixed_observation(
        &self,
    ) -> Result<(), HarmonicError> {
        let (n, m, k) = (
            self.closed.n as u128,
            self.closed.m as u128,
            self.ports.len() as u128,
        );
        let bytes = 6 * (8 * n * k + 256 * k + 64 * m) + 512 * (n + m + k);
        let work =
            6 * (2 * n * n * (k + 1) + n * m + k * k + 32 * m) + 48 * n * (n + m) + 192 * (m + k);
        if bytes > self.budget.maximum_payload_bytes as u128
            || work > self.budget.maximum_work as u128
        {
            return Err(invalid(
                "insufficient complete loaded jet payload or work budget",
            ));
        }
        Ok(())
    }
}

pub(super) fn add_i_state(terms: &mut [Sum; 2], state: C64, factor: f64) -> Result<(), HarmonicError> {
    terms[0]
        .product(-state.im, factor)
        .and_then(|()| terms[1].product(state.re, factor))
        .ok_or_else(|| invalid("loaded jet right-hand side is not representable"))
}
