//! Resonances of the original coupled spatial structure and pressure fields.
//!
//! The physical pencil, forces, inward flow, displacement and pressure all use
//! `exp(-i omega t)`. The pressure potential is not a pressure observation:
//! `p=-i omega psi` and `p_omega=-i(psi+omega psi_omega)`. The canonical whole
//! pencil owner supplies that kinematic term with its same-factor tangent.

use super::{
    PreparedStructuralInterior, StructuralInterior, StructuralInteriorError,
    StructuralInteriorResponse, invalid, poll, response,
};
mod parameter;
pub use parameter::{
    StructuralInteriorDerivativeFields, StructuralInteriorFrequencyParameterResponse,
    StructuralInteriorParameterDirection, StructuralInteriorPowerParameterJet,
    StructuralInteriorResponseField, StructuralInteriorSourceParameterJet,
};

use fs_exec::Cx;
use fs_math::c64::C64;
use fs_modal::harmonic::resonance::HarmonicTransferSample;
use fs_modal::harmonic::{HarmonicFrequencyControls, HarmonicResponseError};
use fs_phs::harmonic::resonance::{ResonanceSearchError, survey_resonances};
pub use fs_phs::harmonic::resonance::{
    ResonanceSurvey, ResonanceSurveyControls, ResonanceSurveyError,
};

/// Original source arrays, either peak phasors or their derivatives per rad/s
/// when supplied as `source_derivative`. No pressure-to-flow normalization occurs.
#[derive(Clone, Copy, Debug)]
pub struct StructuralInteriorSources<'s> {
    pub structural_force: &'s [C64],
    pub inward_flow_load: &'s [C64],
}

/// Independent derivatives of the original physical work channels [W per rad/s].
/// These slopes are signed even when their underlying dissipated power is positive.
#[derive(Clone, Copy, Debug)]
pub struct StructuralInteriorPowerDerivative {
    pub structural_input_w_per_rad_s: f64,
    pub acoustic_input_w_per_rad_s: f64,
    pub fluid_to_structure_w_per_rad_s: f64,
    pub structural_dissipation_w_per_rad_s: f64,
    pub wall_dissipation_w_per_rad_s: f64,
}

impl StructuralInteriorPowerDerivative {
    #[must_use]
    pub fn structural_balance_residual_w_per_rad_s(self) -> f64 {
        self.structural_input_w_per_rad_s + self.fluid_to_structure_w_per_rad_s
            - self.structural_dissipation_w_per_rad_s
    }

    #[must_use]
    pub fn acoustic_balance_residual_w_per_rad_s(self) -> f64 {
        self.acoustic_input_w_per_rad_s
            - self.fluid_to_structure_w_per_rad_s
            - self.wall_dissipation_w_per_rad_s
    }
}

/// Complete coupled physical response and continuous frequency tangent.
#[derive(Clone, Debug)]
pub struct StructuralInteriorFrequencyResponse {
    pub response: StructuralInteriorResponse,
    pub structural_displacement_derivative_per_rad_s: Vec<C64>,
    pub structural_velocity_derivative_per_rad_s: Vec<C64>,
    pub pressure_potential_derivative_pa_s_per_rad_s: Vec<C64>,
    pub pressure_derivative_pa_per_rad_s: Vec<C64>,
    pub fluid_force_derivative_per_rad_s: Vec<C64>,
    pub displaced_flow_derivative_m3_s_per_rad_s: Vec<C64>,
    pub power_derivative: StructuralInteriorPowerDerivative,
    /// Canonical whole-pencil input minus dissipative-work derivative.
    pub power_balance_residual_derivative_w_per_rad_s: f64,
    /// Native tangent solve and original-coordinate differentiated-row residual.
    /// This is a backward check, not a forward error bound near resonance.
    pub relative_derivative_residual: f64,
}

/// A fixed physical receiver. Structural weights retain original coordinate
/// units and signs; they are linear observation rows, without normalization.
/// Pressure is evaluated in the original hierarchical tetrahedral basis.
#[derive(Clone, Copy, Debug)]
pub enum StructuralInteriorReceiver<'r> {
    StructuralDisplacement(&'r [f64]),
    StructuralVelocity(&'r [f64]),
    PressurePoint {
        tetrahedron: usize,
        barycentric: [f64; 4],
    },
}

impl PreparedStructuralInterior<'_, '_, '_, '_> {
    /// Reuse this complete system's factor for the analytic frequency tangent.
    /// M,C,K, the original interface and real wall admittances stay fixed.
    /// `source_derivative=None` holds BOTH force and flow phasors fixed; Some
    /// supplies both complete original derivatives, including zero arrays for
    /// a source that stays fixed. Source and receiver conventions are unchanged.
    ///
    /// Pressure, traction and displaced-flow derivatives remain physical; the
    /// structural and acoustic work derivatives include every source chain-rule
    /// term. Symmetric damping alone supplies heat, with no gyroscopic debit.
    ///
    /// In addition to the native harmonic tangent budget this adapter admits
    /// `256*n+256` bytes and work score `128*n*n+256*n` for original field,
    /// interface, storage and work recovery, where n counts all coordinates.
    /// These are added numerical payload/work bounds, excluding the prepared
    /// model, caller inputs and allocator metadata. Refusal leaves it reusable.
    pub fn solve_with_frequency_derivative(
        &self,
        structural_force: &[C64],
        inward_flow_load: &[C64],
        source_derivative: Option<StructuralInteriorSources<'_>>,
        controls: HarmonicFrequencyControls,
        cx: &Cx<'_>,
    ) -> Result<StructuralInteriorFrequencyResponse, StructuralInteriorError> {
        poll(cx)?;
        let ns = self.model.structural_dofs();
        let n = ns + self.model.pressure_dofs();
        let controls = native_controls(n, controls)?;
        let force = self.load(structural_force, inward_flow_load)?;
        let slope = source_derivative
            .map(|d| self.load(d.structural_force, d.inward_flow_load))
            .transpose()?;
        let mut harmonic = self.harmonic.solve_with_frequency_derivative(
            &force,
            slope.as_deref(),
            None,
            controls,
            cx,
        )?;
        let potential_derivative = harmonic.displacement_derivative.split_off(ns);
        let pressure_derivative = harmonic.velocity_derivative.split_off(ns);
        let response = self.response(harmonic.response, structural_force, inward_flow_load, cx)?;
        let mut cancelled = || cx.checkpoint().is_err();
        let (force_derivative, flow_derivative) = response::interface(
            &self.model.coupling,
            &harmonic.velocity_derivative,
            &pressure_derivative,
            &mut cancelled,
        )?;
        let power = StructuralInteriorPowerDerivative {
            structural_input_w_per_rad_s: response::work_derivative(
                structural_force,
                source_derivative.map(|d| d.structural_force),
                &response.structural_velocity,
                &harmonic.velocity_derivative,
            ),
            acoustic_input_w_per_rad_s: response::work_derivative(
                inward_flow_load,
                source_derivative.map(|d| d.inward_flow_load),
                &response.pressure_pa,
                &pressure_derivative,
            ),
            fluid_to_structure_w_per_rad_s: response::work_derivative(
                &response.fluid_force,
                Some(&force_derivative),
                &response.structural_velocity,
                &harmonic.velocity_derivative,
            ),
            structural_dissipation_w_per_rad_s: response::dissipation_derivative(
                &self.model.damping,
                n,
                0,
                &response.structural_velocity,
                &harmonic.velocity_derivative,
                &mut cancelled,
            )?,
            wall_dissipation_w_per_rad_s: response::dissipation_derivative(
                &self.model.damping,
                n,
                ns,
                &response.pressure_pa,
                &pressure_derivative,
                &mut cancelled,
            )?,
        };
        if [
            power.structural_input_w_per_rad_s,
            power.acoustic_input_w_per_rad_s,
            power.fluid_to_structure_w_per_rad_s,
            power.structural_dissipation_w_per_rad_s,
            power.wall_dissipation_w_per_rad_s,
        ]
        .iter()
        .any(|v| !v.is_finite())
        {
            return Err(invalid(
                "original physical work derivative is unrepresentable",
            ));
        }
        poll(cx)?;
        Ok(StructuralInteriorFrequencyResponse {
            response,
            structural_displacement_derivative_per_rad_s: harmonic.displacement_derivative,
            structural_velocity_derivative_per_rad_s: harmonic.velocity_derivative,
            pressure_potential_derivative_pa_s_per_rad_s: potential_derivative,
            pressure_derivative_pa_per_rad_s: pressure_derivative,
            fluid_force_derivative_per_rad_s: force_derivative,
            displaced_flow_derivative_m3_s_per_rad_s: flow_derivative,
            power_derivative: power,
            power_balance_residual_derivative_w_per_rad_s: harmonic
                .power_balance_residual_derivative,
            relative_derivative_residual: harmonic.relative_derivative_residual,
        })
    }
}

impl StructuralInterior<'_, '_, '_> {
    /// Find driven receiver peaks resolved by the supplied rad/s grid.
    ///
    /// Every query prepares the same complete original spatial pencil once and
    /// reuses that factor for its analytic tangent. Fixed coherent structural
    /// forces and inward-flow phasors interfere before observation. The physical
    /// receiver uses `exp(-i omega t)` directly; no conjugation changes its phase.
    /// Work retains the separate conjugate effort/flow convention.
    ///
    /// The shared search owns refinement and total query/peak budgets. The model
    /// and tangent budgets apply to each query. Cancellation or a genuine coupled
    /// singularity publishes no partial survey. The result is grid-resolved, not
    /// a free-mode census; it supplies neither decay Q nor a completeness claim.
    #[allow(clippy::too_many_arguments)]
    pub fn survey_resonances(
        &self,
        frequencies_rad_s: &[f64],
        structural_force: &[C64],
        inward_flow_load: &[C64],
        receiver: StructuralInteriorReceiver<'_>,
        frequency_controls: HarmonicFrequencyControls,
        survey_controls: ResonanceSurveyControls,
        cx: &Cx<'_>,
    ) -> Result<ResonanceSurvey, ResonanceSurveyError<StructuralInteriorError>> {
        poll(cx).map_err(ResonanceSearchError::Source)?;
        let mut cancelled = || cx.checkpoint().is_err();
        survey_resonances(
            frequencies_rad_s,
            survey_controls,
            &mut cancelled,
            |omega, _, _| {
                let result = self
                    .prepare_angular(omega, cx)?
                    .solve_with_frequency_derivative(
                        structural_force,
                        inward_flow_load,
                        None,
                        frequency_controls,
                        cx,
                    )?;
                let sample = match receiver {
                    StructuralInteriorReceiver::StructuralDisplacement(weights) => {
                        HarmonicTransferSample::project_real(
                            weights,
                            &result.response.structural_displacement,
                            &result.structural_displacement_derivative_per_rad_s,
                        )?
                    }
                    StructuralInteriorReceiver::StructuralVelocity(weights) => {
                        HarmonicTransferSample::project_real(
                            weights,
                            &result.response.structural_velocity,
                            &result.structural_velocity_derivative_per_rad_s,
                        )?
                    }
                    StructuralInteriorReceiver::PressurePoint {
                        tetrahedron,
                        barycentric,
                    } => HarmonicTransferSample {
                        amplitude: self.interior.pressure_at(
                            &result.response.pressure_pa,
                            tetrahedron,
                            barycentric,
                        )?,
                        angular_derivative: self.interior.pressure_at(
                            &result.pressure_derivative_pa_per_rad_s,
                            tetrahedron,
                            barycentric,
                        )?,
                    },
                };
                Ok((
                    sample.amplitude,
                    sample.angular_derivative,
                    result
                        .response
                        .relative_residual
                        .max(result.relative_derivative_residual),
                ))
            },
        )
    }
}

fn native_controls(
    n: usize,
    controls: HarmonicFrequencyControls,
) -> Result<HarmonicFrequencyControls, StructuralInteriorError> {
    let bytes = n
        .checked_mul(256)
        .and_then(|v| v.checked_add(256))
        .ok_or_else(|| invalid("physical frequency payload extent overflow"))?;
    let maximum_additional_payload_bytes = controls
        .maximum_additional_payload_bytes
        .checked_sub(bytes)
        .ok_or(HarmonicResponseError::Budget {
            quantity: "structural-interior frequency recovery payload bytes",
            requested: bytes,
            cap: controls.maximum_additional_payload_bytes,
        })?;
    let n = u64::try_from(n).map_err(|_| invalid("frequency work extent overflow"))?;
    let work = n
        .checked_mul(n)
        .and_then(|v| v.checked_mul(128))
        .and_then(|v| n.checked_mul(256).and_then(|w| v.checked_add(w)))
        .ok_or_else(|| invalid("physical frequency work extent overflow"))?;
    let maximum_work =
        controls
            .maximum_work
            .checked_sub(work)
            .ok_or(HarmonicResponseError::WorkBudget {
                requested: work,
                cap: controls.maximum_work,
            })?;
    Ok(HarmonicFrequencyControls {
        maximum_additional_payload_bytes,
        maximum_work,
    })
}
