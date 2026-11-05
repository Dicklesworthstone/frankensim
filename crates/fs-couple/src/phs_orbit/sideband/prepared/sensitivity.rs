//! Analytic carrier derivatives in the original physical port/state chart.

use super::super::{
    ObservationBudget, PhsPortSidebandCheck, admit_storage, input_work, minimum_observation_work,
    observe_with_budget,
};
use super::{
    C64, FactoredPhsPortSidebands, PhsPortSidebandError, PhsPortSidebandResponse, Result,
    SidebandResponse, admit_forcing, input_forcing, poll,
};
use fs_orbit::sideband::sensitivity::SidebandCarrierDerivativeControls;

/// Complete physical response and its analytic derivative per angular frequency.
///
/// The base orbit, sustaining source, physical maps and sideband indices stay
/// fixed. These are derivatives of the COMPLEX coefficients under exp(+i nu t),
/// not a derivative of the real time waveform or of a changing base branch.
/// No i*t term or derivative of the base clock is added to the coefficients.
#[derive(Clone, Debug)]
pub struct PhsPortSidebandCarrierDerivativeResponse {
    pub response: PhsPortSidebandResponse,
    /// Same state ordering and retained bands as `response.state_response`,
    /// in original state units per angular frequency. Its residuals check the
    /// differentiated original equations, with the authored derivative limits.
    pub state_angular_derivative: SidebandResponse,
    /// Derivative of the original `G^T Hxx(x_base) dx` output coefficients,
    /// in original physical port units per angular frequency, with the same
    /// independent output order and indexing as `response.port_coefficients`.
    pub port_angular_derivative: Vec<C64>,
    /// Independent shifted-grid checks in output units per angular frequency.
    pub derivative_observation_trace: Vec<PhsPortSidebandCheck>,
    /// TOTAL physical observation Hessian actions across response and derivative.
    /// `response.hessian_actions` remains the primal observation's own count.
    pub total_hessian_actions: usize,
}

impl FactoredPhsPortSidebands<'_, '_> {
    /// Solve the physical response and its carrier derivative with this SAME LU.
    ///
    /// `forcing_angular_derivative` is the complete authored original-port
    /// derivative per angular frequency at the same signed sideband indices.
    /// None explicitly holds the probe phasors fixed. G maps both arrays in
    /// their original units. For these complete pHS ODEs P(iw)=iwI, hence its
    /// exact carrier derivative is iI for EVERY original coordinate, including
    /// memory. The shared generic owner solves and audits
    /// `D X_nu = G u_nu - D_nu X` without another factorization or frequency probe.
    ///
    /// The fixed-base physical output derivative is `G^T Hxx(x_base) X_nu`.
    /// Full mechanical momentum states retain the velocity chain rule and all
    /// constitutive memory derivatives automatically. A derivative along a
    /// changing nonlinear base would require additional base/tangent terms and
    /// is not this operation.
    ///
    /// `derivative_envelope_tolerances` independently bound the sampled output
    /// derivative in each original port's units per angular frequency. Source
    /// observation work/actions are CUMULATIVE across the two observations; the
    /// first output and trace remain charged while the second is computed.
    /// The generic derivative controls add their separately bounded workspace
    /// and work to the original core allowances. Live payload is bounded by the
    /// sum of original core, added derivative core, and physical adapter caps.
    ///
    /// # Errors
    /// Invalid original forcing/slope or derivative tolerances; exhausted pair
    /// budgets; generic primal/derivative residual or omitted-band refusal;
    /// unresolved original output or its derivative; source refusal or
    /// cancellation. No partial pair escapes; the same factor remains usable.
    pub fn solve_with_carrier_derivative(
        &mut self,
        forcing: &[C64],
        forcing_angular_derivative: Option<&[C64]>,
        controls: SidebandCarrierDerivativeControls,
        derivative_envelope_tolerances: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PhsPortSidebandCarrierDerivativeResponse> {
        poll(&mut cancelled)?;
        let original = self.source.problem.original();
        let system = original.system;
        let c = &self.source.controls;
        let (n, m, bands) = (
            system.state_dim(),
            system.port_dim(),
            2 * c.sideband.harmonics + 1,
        );
        admit_forcing(system, forcing, c)?;
        if let Some(slope) = forcing_angular_derivative {
            admit_forcing(system, slope, c)?;
        }
        if derivative_envelope_tolerances.len() != m
            || derivative_envelope_tolerances
                .iter()
                .any(|x| !x.is_finite() || *x <= 0.0)
        {
            return Err(PhsPortSidebandError::Invalid(
                "positive original-port output derivative tolerances required",
            ));
        }
        // Both envelopes require their own initial, doubled and quadrupled
        // grids. Admit their unavoidable combined cost before any new source
        // map or derivative preparation. Include the known complete iI fills.
        let mapped_inputs = if forcing_angular_derivative.is_some() {
            2
        } else {
            1
        };
        let mut budget = ObservationBudget {
            actions: 0,
            work: mapped_inputs * input_work(n, m, c.sideband.harmonics)
                + 16 * n as u128 * n as u128 * bands as u128,
            retained_bytes: 0,
        };
        if 28 * c.initial_samples > c.max_hessian_actions {
            return Err(PhsPortSidebandError::Budget(
                "response and carrier derivative exceed the combined Hessian-action cap",
            ));
        }
        if budget.work + 2 * minimum_observation_work(n, m, c.initial_samples)
            > c.max_observation_work as u128
        {
            return Err(PhsPortSidebandError::Budget(
                "response and carrier derivative exceed the combined physical work cap",
            ));
        }
        admit_storage(n, m, c, c.initial_samples)?;
        let pair = {
            let rhs = input_forcing(system, forcing, c.sideband.harmonics, &mut cancelled)?;
            let derivative_rhs = match forcing_angular_derivative {
                Some(slope) => input_forcing(system, slope, c.sideband.harmonics, &mut cancelled)?,
                None => vec![C64::ZERO; n * bands],
            };
            let workspace = self.factor.prepare_carrier_derivative(
                controls,
                |_, out| {
                    if out.len() != n * n {
                        return false;
                    }
                    out.fill(C64::ZERO);
                    for i in 0..n {
                        out[i * n + i] = C64::new(0.0, 1.0);
                    }
                    true
                },
                &mut cancelled,
            )?;
            workspace.solve(&rhs, &derivative_rhs, &mut cancelled)?
        }; // Release mapped inputs and derivative matrices before FFT observation.
        let (port_coefficients, observation_trace, hessian_actions) = observe_with_budget(
            original,
            self.source.orbit,
            &pair.response.coefficients,
            c,
            &c.envelope_tolerances,
            &mut budget,
            &mut cancelled,
        )?;
        // The existing 32*n*bands allowance covers BOTH retained state arrays.
        // Add the first physical output/trace while the second observer's grids
        // and spectra are live; no caller-owned result has been published yet.
        budget.retained_bytes = core::mem::size_of::<C64>() as u128
            * port_coefficients.capacity() as u128
            + core::mem::size_of::<PhsPortSidebandCheck>() as u128
                * observation_trace.capacity() as u128
            + 8 * observation_trace
                .iter()
                .map(|check| {
                    check.omitted_envelope_norm.capacity() as u128
                        + check
                            .sampling_change_norm
                            .as_ref()
                            .map_or(0, |x| x.capacity() as u128)
                })
                .sum::<u128>();
        let (port_angular_derivative, derivative_observation_trace, _) = observe_with_budget(
            original,
            self.source.orbit,
            &pair.angular_derivative.coefficients,
            c,
            derivative_envelope_tolerances,
            &mut budget,
            &mut cancelled,
        )
        .map_err(|error| match error {
            PhsPortSidebandError::ObservationUnresolved { last } => {
                PhsPortSidebandError::DerivativeObservationUnresolved { last }
            }
            other => other,
        })?;
        poll(&mut cancelled)?;
        Ok(PhsPortSidebandCarrierDerivativeResponse {
            response: PhsPortSidebandResponse {
                omega_rad_s: self.source.orbit.omega,
                carrier_angular_frequency: self.carrier_angular_frequency(),
                state_harmonics: c.sideband.harmonics,
                output_harmonics: c.output_harmonics,
                state_dimension: n,
                port_count: m,
                state_response: pair.response,
                port_coefficients,
                base_residual: self.base_residual,
                tangent_sampling_defect: self.tangent_sampling_defect,
                observation_trace,
                hessian_actions,
            },
            state_angular_derivative: pair.angular_derivative,
            port_angular_derivative,
            derivative_observation_trace,
            total_hessian_actions: budget.actions,
        })
    }
}
