//! Original storage-parameter sensitivities of a changing forced periodic orbit.

use super::super::{
    EffortObservation, ObservationBudget, PhsPortSidebandCheck, PhsStorageParameterAction,
    minimum_observation_work, observe_effort_with_budget,
};
use super::sensitivity::trace_payload;
use super::{
    C64, FactoredPhsPortSidebands, OrbitProblem, PhsPortSidebandError, PreparedPhsPortSidebands,
    Result, SidebandResponse, SourceProblem, poll,
};
use core::cell::Cell;
use fs_math::det;
use fs_orbit::sideband::SidebandError;

/// Additional physical source-resolution and work limits for ONE derivative.
#[derive(Clone, Debug)]
pub struct PhsStorageParameterControls {
    /// Positive absolute RMS tolerances on the projection of
    /// `(J-R) partial_p grad H`, per original equation, in equation units per
    /// parameter unit. Retained-grid change and observed omitted forcing share
    /// each allowance. The complete normalized defect also shares the native
    /// omitted-equation budget with the solved response; it is never hidden by
    /// treating unlisted source harmonics as exactly zero.
    pub forcing_rms_tolerances: Vec<f64>,
    /// Total constitutive parameter-action calls across forcing projection AND
    /// complete output observation. Each needs two successful grid comparisons.
    pub maximum_parameter_actions: usize,
}

/// Local derivatives at fixed drive frequency, phase and authored input phasors.
/// All original states, including material memory, keep their order and units.
#[derive(Clone, Debug)]
pub struct PhsForcedStorageParameterTangent {
    pub omega_rad_s: f64,
    pub state_harmonics: usize,
    pub output_harmonics: usize,
    pub state_dimension: usize,
    pub port_count: usize,
    /// Sideband-major `(k+state_harmonics)*state_dimension+j`; physical state
    /// units per parameter. Both signed Fourier coefficients participate.
    /// Its native residual refers to the resolved finite projection of the
    /// original constitutive parameter forcing; its projection errors follow.
    pub state_parameter_derivative: SidebandResponse,
    /// `(k+output_harmonics)*port_count+p`, in physical output units per
    /// parameter. This is `G^T (Hxx x_p + partial_p grad H)`, including the
    /// DIRECT material term. These are Fourier coefficients, not doubled peaks.
    pub port_parameter_derivative: Vec<C64>,
    /// Grid and omitted-forcing norms in original equation units per parameter.
    pub forcing_projection_trace: Vec<PhsPortSidebandCheck>,
    /// Native omitted-equation defect plus the normalized measured source
    /// sampling/tail defect, conservatively combined by the triangle inequality.
    /// This stays within the native omitted-equation allowance. It remains a
    /// finite-grid check, not an infinite-tail or forward-error bound.
    pub combined_omitted_equation_defect: f64,
    pub base_residual: f64,
    pub tangent_sampling_defect: f64,
    pub observation_trace: Vec<PhsPortSidebandCheck>,
    pub hessian_actions: usize,
    pub parameter_actions: usize,
}

fn require_forced(problem: &SourceProblem<'_>) -> Result<()> {
    if !matches!(problem, SourceProblem::Forced(_)) || problem.autonomous() {
        Err(PhsPortSidebandError::Invalid(
            "storage parameter tangents require an original externally forced source without feedback",
        ))
    } else {
        Ok(())
    }
}

impl PreparedPhsPortSidebands<'_, '_> {
    /// Factor at zero carrier and differentiate the original forced response
    /// with respect to a storage parameter. See the factored entry for the
    /// fixed-model-chart, complete direct-output and source-resolution contract.
    pub fn forced_storage_parameter_tangent(
        &mut self,
        action: &PhsStorageParameterAction<'_>,
        controls: &PhsStorageParameterControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PhsForcedStorageParameterTangent> {
        poll(&mut cancelled)?;
        require_forced(&self.source.problem)?;
        self.factor(0.0, &mut cancelled)?
            .forced_storage_parameter_tangent(action, controls, cancelled)
    }
}

impl FactoredPhsPortSidebands<'_, '_> {
    /// Differentiate a nonlinear forced orbit using this SAME zero-carrier LU.
    ///
    /// Holding J/R/G, the physical state chart, drive clock and input phasors
    /// fixed gives `D x_p = (J-R) partial_p grad H`, where D is the existing
    /// full periodic tangent. The original compensated vector-field owner maps
    /// the complete constitutive action. Its Fourier forcing is independently
    /// resolved by the shared shifted-grid FFT observer before the native solve.
    /// Source projection defects remain explicit and share the native omitted
    /// equation allowance with the response; no unresolved material harmonic is
    /// silently discarded. All original nonlinear and memory states remain.
    ///
    /// Physical output is observed as `G^T (Hxx x_p + partial_p grad H)` on
    /// independently refined grids, including direct material effects invisible
    /// in a state-only sensitivity. The source's existing output tolerances now
    /// have output units per declared parameter unit. Native normalized solve
    /// and omitted-equation limits likewise apply to this derivative. The
    /// constitutive action must cover the SAME complete original storage law;
    /// no parameter-to-model correspondence is inferred from a scalar label.
    ///
    /// Feedback sources refuse: changing effort can also change their input
    /// law, whose direct material derivative is not supplied by this action.
    /// Changing dissipation/interconnection/ports, basis, mass normalization or
    /// drive coefficients is also outside this fixed-chart storage direction.
    /// In particular changing frequency while holding a damping ratio fixed
    /// also changes R and requires more than this storage-only derivative.
    ///
    /// Both projections share the original adapter work/payload/Hessian caps
    /// and one parameter-action cap. The same shifted-grid refinement, native
    /// factor/residual correction and original observation owners are reused.
    /// No perturbed model, second Jacobian or finite-difference solve is formed.
    /// Singular or unresolved derivatives refuse without regularization or a
    /// branch switch. This local finite-harmonic sensitivity does not establish
    /// branch uniqueness, Floquet stability, infinite tails or continuum error.
    /// Failure/cancellation publishes no partial result and leaves the factor
    /// reusable. Source callbacks remain deterministic/read-only.
    pub fn forced_storage_parameter_tangent(
        &mut self,
        action: &PhsStorageParameterAction<'_>,
        controls: &PhsStorageParameterControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PhsForcedStorageParameterTangent> {
        poll(&mut cancelled)?;
        require_forced(&self.source.problem)?;
        if self.carrier_angular_frequency() != 0.0 {
            return Err(PhsPortSidebandError::Invalid(
                "a changing forced base requires a zero-carrier factor",
            ));
        }
        let original = self.source.problem.original();
        let system = original.system;
        let base = self.source.orbit;
        let c = &self.source.controls;
        let (n, m, h) = (system.state_dim(), system.port_dim(), c.sideband.harmonics);
        if controls.forcing_rms_tolerances.len() != n
            || controls
                .forcing_rms_tolerances
                .iter()
                .any(|v| !v.is_finite() || *v <= 0.0)
        {
            return Err(PhsPortSidebandError::Invalid(
                "positive physical parameter-forcing tolerance required for every original equation",
            ));
        }
        if controls.maximum_parameter_actions < 14 * c.initial_samples {
            return Err(PhsPortSidebandError::Budget(
                "parameter forcing and output require two independent grid comparisons each",
            ));
        }
        let source_map_work = 128 * n as u128 * (2 * n + m) as u128;
        let minimum_work = minimum_observation_work(n, n, c.initial_samples)
            + minimum_observation_work(n, m, c.initial_samples)
            + source_map_work * (7 * c.initial_samples) as u128;
        if minimum_work > c.max_observation_work as u128 {
            return Err(PhsPortSidebandError::Budget(
                "parameter forcing and output exceed the combined physical work cap",
            ));
        }
        let actions = Cell::new(0_usize);
        let exhausted = Cell::new(false);
        let counted = |state: &[f64], effort: &mut [f64]| {
            if actions.get() == controls.maximum_parameter_actions {
                exhausted.set(true);
                return false;
            }
            actions.set(actions.get() + 1);
            action(state, effort)
        };
        let mut budget = ObservationBudget {
            actions: 0,
            work: 0,
            retained_bytes: 16 * m as u128,
        };
        let (state_parameter_derivative, forcing_projection_trace, combined_defect) = {
            // No second projection implementation: the SAME observer samples
            // the constitutive derivative alone and maps it through original J-R.
            let mut projection = c.clone();
            projection.output_harmonics = h;
            let zero_input = vec![0.0; m];
            let result = observe_effort_with_budget(
                original,
                base,
                &[],
                &projection,
                &controls.forcing_rms_tolerances,
                n,
                source_map_work,
                &|effort, flow| {
                    system
                        .vector_field_from_effort_into(effort, &zero_input, flow)
                        .is_ok()
                },
                EffortObservation::Parameter(&counted),
                &mut budget,
                &mut cancelled,
            );
            if exhausted.get() {
                return Err(PhsPortSidebandError::Budget(
                    "constitutive parameter-action cap exhausted",
                ));
            }
            let (forcing, trace, _) = result?;
            let last = trace.last().ok_or(PhsPortSidebandError::Invalid(
                "missing parameter forcing assessment",
            ))?;
            let change =
                last.sampling_change_norm
                    .as_ref()
                    .ok_or(PhsPortSidebandError::Invalid(
                        "missing parameter forcing grid comparison",
                    ))?;
            let mut source_defect = 0.0;
            for ((tail, change), scale) in last
                .omitted_envelope_norm
                .iter()
                .zip(change)
                .zip(self.factor.equation_scales())
            {
                poll(&mut cancelled)?;
                source_defect = det::hypot(source_defect, (tail + change) / scale);
            }
            let limit = c.sideband.maximum_omitted_residual;
            if !source_defect.is_finite() || source_defect > limit {
                return Err(SidebandError::OmittedSidebands {
                    measured: source_defect,
                    limit,
                }
                .into());
            }
            let state = self.factor.solve(&forcing, &mut cancelled)?;
            let combined = state.omitted_residual + source_defect;
            if !combined.is_finite() || combined > limit {
                return Err(SidebandError::OmittedSidebands {
                    measured: combined,
                    limit,
                }
                .into());
            }
            (state, trace, combined)
        };
        // Force coefficients and cloned projection controls are now released.
        // Keep their diagnostic history charged while observing the total effort.
        budget.retained_bytes = trace_payload(
            &forcing_projection_trace,
            forcing_projection_trace.capacity(),
        );
        let result = observe_effort_with_budget(
            original,
            base,
            &state_parameter_derivative.coefficients,
            c,
            &c.envelope_tolerances,
            m,
            0,
            &|effort, output| system.output_from_effort_into(effort, output).is_ok(),
            EffortObservation::StateAndParameter(&counted),
            &mut budget,
            &mut cancelled,
        );
        if exhausted.get() {
            return Err(PhsPortSidebandError::Budget(
                "constitutive parameter-action cap exhausted",
            ));
        }
        let (port_parameter_derivative, observation_trace, hessian_actions) = result?;
        poll(&mut cancelled)?;
        Ok(PhsForcedStorageParameterTangent {
            omega_rad_s: base.omega,
            state_harmonics: h,
            output_harmonics: c.output_harmonics,
            state_dimension: n,
            port_count: m,
            state_parameter_derivative,
            port_parameter_derivative,
            forcing_projection_trace,
            combined_omitted_equation_defect: combined_defect,
            base_residual: self.base_residual,
            tangent_sampling_defect: self.tangent_sampling_defect,
            observation_trace,
            hessian_actions,
            parameter_actions: actions.get(),
        })
    }
}
