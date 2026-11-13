//! Original physical-parameter sensitivities of a changing forced periodic orbit.

use super::super::super::{InputKind, PhsFeedbackParameterAction};
use super::super::{
    EffortObservation, ObservationBudget, PhsPortSidebandCheck, PhsStorageParameterAction,
    minimum_observation_work, observe_sampled_effort_with_budget,
};
use super::sensitivity::trace_payload;
use super::{
    C64, FactoredPhsPortSidebands, OrbitProblem, PeriodicPortInput, PhsPortSidebandError,
    PreparedPhsPortSidebands, Result, SidebandResponse, SourceProblem, poll,
};
use core::cell::{Cell, RefCell};
use fs_math::det;
use fs_orbit::sideband::SidebandError;
use fs_phs::PhsStructureParameterDirection;

/// ONE real parameter's original physical derivatives at fixed phase/clock.
/// None holds the corresponding source fixed. The storage action covers every
/// nonlinear and memory coordinate; the SAME port direction changes excitation
/// and observation. The input holds derivatives of authored Fourier coefficients
/// in the original port chart, with the SAME clock and coefficient convention.
/// A feedback source additionally requires its COMPLETE fixed-state parameter
/// action; the source's existing total state Jacobian supplies the state chain.
/// No model/basis/drive derivative is inferred from a scalar parameter label.
#[derive(Clone, Copy, Default)]
pub struct PhsForcedParameterDirection<'a> {
    pub structure: PhsStructureParameterDirection<'a>,
    pub storage: Option<&'a PhsStorageParameterAction<'a>>,
    pub input: Option<&'a PeriodicPortInput>,
    pub feedback: Option<&'a PhsFeedbackParameterAction<'a>>,
}

/// Additional physical source-resolution and work limits for ONE derivative.
#[derive(Clone, Debug)]
pub struct PhsParameterControls {
    /// Positive absolute RMS tolerances on the projection of the complete
    /// direct parameter forcing, in original equation units per parameter unit.
    /// Retained-grid change and observed omitted forcing share
    /// each allowance. The complete normalized defect also shares the native
    /// omitted-equation budget with the solved response; it is never hidden by
    /// treating unlisted source harmonics as exactly zero.
    pub forcing_rms_tolerances: Vec<f64>,
    /// Total real parameter samples across forcing projection AND complete
    /// output observation. Each evaluates the optional storage action once and
    /// at most one original gradient. Forcing also evaluates the optional
    /// feedback parameter action once; both projections need two successful
    /// grid comparisons. Callback-internal work remains source-defined.
    pub maximum_parameter_actions: usize,
}

/// Compatibility name for the storage-only subset of the same physical owner.
pub type PhsStorageParameterControls = PhsParameterControls;

/// Local derivatives at fixed drive frequency, phase and physical state chart.
/// All original states, including material memory, keep their order and units.
#[derive(Clone, Debug)]
pub struct PhsForcedParameterTangent {
    pub omega_rad_s: f64,
    pub state_harmonics: usize,
    pub output_harmonics: usize,
    pub state_dimension: usize,
    pub port_count: usize,
    /// Sideband-major `(k+state_harmonics)*state_dimension+j`; physical state
    /// units per parameter. Both signed Fourier coefficients participate.
    /// Its native residual refers to the resolved finite projection of the
    /// complete original parameter forcing; its projection errors follow.
    pub state_parameter_derivative: SidebandResponse,
    /// `(k+output_harmonics)*port_count+p`, in physical output units per
    /// parameter. This is `G^T (Hxx x_p + partial_p grad H) + G_p^T grad H`,
    /// including DIRECT material/map terms. These are Fourier coefficients,
    /// not doubled peaks.
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

/// Compatibility name for a storage-only direction of the same response owner.
pub type PhsForcedStorageParameterTangent = PhsForcedParameterTangent;

fn require_forced(problem: &SourceProblem<'_>, feedback_action: bool) -> Result<()> {
    if problem.autonomous() {
        return Err(PhsPortSidebandError::Invalid(
            "a fixed-clock parameter tangent requires an original externally forced source",
        ));
    }
    if matches!(problem, SourceProblem::Feedback(_)) != feedback_action {
        return Err(PhsPortSidebandError::Invalid(
            "feedback sources require an explicit complete fixed-state parameter action, and ordinary forced sources exclude one",
        ));
    }
    Ok(())
}

impl PreparedPhsPortSidebands<'_, '_> {
    /// Factor at zero carrier and differentiate the original forced response
    /// with respect to a storage parameter. See the factored entry for the
    /// fixed-model-chart, complete direct-output and source-resolution contract.
    pub fn forced_storage_parameter_tangent(
        &mut self,
        action: &PhsStorageParameterAction<'_>,
        controls: &PhsStorageParameterControls,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PhsForcedStorageParameterTangent> {
        self.forced_parameter_tangent(
            PhsForcedParameterDirection {
                storage: Some(action),
                ..Default::default()
            },
            controls,
            cancelled,
        )
    }

    /// Factor at zero carrier and differentiate complete physical coefficients,
    /// storage and drive using the factored entry's original-source contract.
    pub fn forced_parameter_tangent(
        &mut self,
        direction: PhsForcedParameterDirection<'_>,
        controls: &PhsParameterControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PhsForcedParameterTangent> {
        poll(&mut cancelled)?;
        require_forced(&self.source.problem, direction.feedback.is_some())?;
        self.factor(0.0, &mut cancelled)?
            .forced_parameter_tangent(direction, controls, cancelled)
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
        cancelled: impl FnMut() -> bool,
    ) -> Result<PhsForcedStorageParameterTangent> {
        self.forced_parameter_tangent(
            PhsForcedParameterDirection {
                storage: Some(action),
                ..Default::default()
            },
            controls,
            cancelled,
        )
    }

    /// Differentiate the complete original forced response with this SAME LU.
    ///
    /// At fixed physical chart and clock, `e=grad H` and `h_p=partial_p grad H`
    /// give `D x_p=(J-R)h_p+(J_p-R_p)e+G_p u+G u_p`. The original source-bound
    /// coefficient owner contracts ALL terms before rounding each equation.
    /// Independently refined source projection precedes the existing native
    /// solve; measured source errors share the original equation normalization
    /// and omitted-equation allowance with the solved response.
    ///
    /// Output is `G^T(Hxx x_p+h_p)+G_p^T e`, including direct storage and
    /// reciprocal map changes. The same observer resolves its total waveform,
    /// not the terms separately. Only the real envelope gets direct parameter
    /// terms. Every original nonlinear/material-memory state remains.
    ///
    /// Signed J/R/G directions retain their source symmetry without requiring
    /// positivity. Input derivatives include DC, amplitude and phase changes,
    /// with matching original ports/clock and explicitly resolved bandwidth.
    /// For driven feedback, `u` includes BOTH original feedback and imposed
    /// drive. The required feedback action receives original conjugate flows
    /// and their direct derivative, then supplies the COMPLETE fixed-state
    /// input derivative. The cached total feedback Jacobian already includes
    /// the state chain; it is never added a second time. External drive and
    /// feedback parameter derivatives are combined before the source map.
    /// No finite parameter interval, changing basis/clock or feedback-chain
    /// rule is inferred. Singular/unresolved/cancelled
    /// operations publish no partial result and leave the factor reusable.
    ///
    /// Both projections share work/payload/Hessian and real-parameter-sample
    /// caps. Storage callback/gradient internal work remains source-defined.
    /// This finite-harmonic local slope is not a stability, branch uniqueness,
    /// infinite-tail or near-pole forward-error certificate.
    pub fn forced_parameter_tangent(
        &mut self,
        direction: PhsForcedParameterDirection<'_>,
        controls: &PhsParameterControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PhsForcedParameterTangent> {
        poll(&mut cancelled)?;
        require_forced(&self.source.problem, direction.feedback.is_some())?;
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
        let (periodic, feedback) = match (&self.source.problem, &original.input) {
            (SourceProblem::Forced(_), InputKind::Forced { periodic, .. }) => (*periodic, None),
            (SourceProblem::Feedback(problem), InputKind::Feedback { law, .. }) => (
                problem.external.ok_or(PhsPortSidebandError::Invalid(
                    "original periodic input required",
                ))?,
                Some(*law),
            ),
            _ => {
                return Err(PhsPortSidebandError::Invalid(
                    "original forced input required",
                ));
            }
        };
        let physical = system
            .prepare_parameter_direction(direction.structure)
            .map_err(|_| {
                PhsPortSidebandError::Invalid(
                    "complete finite original skew/symmetric parameter direction required",
                )
            })?;
        if direction.input.is_some_and(|input| {
            input.port_count() != m
                || input.omega_rad_s() != base.omega
                || input.highest_nonzero_harmonic() >= c.initial_samples / 2
        }) || (direction.structure.ports.is_some()
            && periodic.highest_nonzero_harmonic() >= c.initial_samples / 2)
        {
            return Err(PhsPortSidebandError::Invalid(
                "parameter forcing needs original ports/clock and a resolved initial source bandwidth",
            ));
        }
        let storage_only = direction.structure.interconnection.is_none()
            && direction.structure.dissipation.is_none()
            && direction.structure.ports.is_none()
            && direction.input.is_none()
            && feedback.is_none();
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
        let source_map_work =
            128 * n as u128 * (2 * n + m) as u128 * if storage_only { 1 } else { 2 }
                + if feedback.is_some() {
                    192 * n as u128 * m as u128 + 8 * m as u128
                } else {
                    0
                }
                + 64 * m as u128
                    * (if direction.structure.ports.is_some() {
                        periodic.harmonic_count()
                    } else {
                        0
                    } + direction.input.map_or(0, PeriodicPortInput::harmonic_count))
                        as u128;
        let output_map_work = if direction.structure.ports.is_some() {
            64 * n as u128 * m as u128
        } else {
            0
        };
        let minimum_work = minimum_observation_work(n, n, c.initial_samples)
            + minimum_observation_work(n, m, c.initial_samples)
            + (source_map_work + output_map_work) * (7 * c.initial_samples) as u128;
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
            if let Some(action) = direction.storage {
                action(state, effort)
            } else {
                effort.fill(0.0);
                true
            }
        };
        let mut budget = ObservationBudget {
            actions: 0,
            work: 0,
            retained_bytes: if storage_only {
                16 * m as u128
            } else {
                8 * (n as u128 + (3 + if feedback.is_some() { 3 } else { 0 }) * m as u128)
            },
        };
        let (state_parameter_derivative, forcing_projection_trace, combined_defect) = {
            // No second projection implementation: the SAME observer samples
            // the original state/phase and maps the complete direct source.
            let mut projection = c.clone();
            projection.output_harmonics = h;
            let scratch = RefCell::new((
                vec![0.0; if storage_only { 0 } else { n }],
                vec![0.0; m],
                vec![0.0; if storage_only { 0 } else { m }],
                vec![0.0; if feedback.is_some() { m } else { 0 }],
                vec![0.0; if feedback.is_some() { m } else { 0 }],
                vec![0.0; if feedback.is_some() { m } else { 0 }],
            ));
            let result = observe_sampled_effort_with_budget(
                original,
                base,
                &[],
                &projection,
                &controls.forcing_rms_tolerances,
                n,
                source_map_work,
                &|phase, state, effort, _, flow| {
                    let mut scratch = scratch.borrow_mut();
                    let (
                        original_effort,
                        input,
                        input_derivative,
                        flows,
                        flow_derivative,
                        feedback_input,
                    ) = &mut *scratch;
                    if storage_only {
                        return system
                            .vector_field_from_effort_into(effort, input, flow)
                            .is_ok();
                    }
                    if direction.structure.interconnection.is_some()
                        || direction.structure.dissipation.is_some()
                        || feedback.is_some()
                    {
                        original_effort.fill(f64::NAN);
                        system.storage().gradient(state, original_effort);
                    }
                    if direction.structure.ports.is_some()
                        && !periodic.sample_at_frequency_into(phase, 1.0, input)
                    {
                        return false;
                    }
                    if direction.input.is_some_and(|source| {
                        !source.sample_at_frequency_into(phase, 1.0, input_derivative)
                    }) {
                        return false;
                    }
                    if direction.input.is_none() {
                        input_derivative.fill(0.0);
                    }
                    if let Some(law) = feedback {
                        if system
                            .output_from_effort_into(original_effort, flows)
                            .is_err()
                            || physical
                                .output_from_effort_into(original_effort, effort, flow_derivative)
                                .is_err()
                        {
                            return false;
                        }
                        if direction.structure.ports.is_some() {
                            feedback_input.fill(f64::NAN);
                            if !law.input(state, flows, feedback_input)
                                || feedback_input.iter().any(|v| !v.is_finite())
                            {
                                return false;
                            }
                            for (total, value) in input.iter_mut().zip(feedback_input.iter()) {
                                *total += value;
                            }
                        }
                        feedback_input.fill(f64::NAN);
                        if !direction.feedback.is_some_and(|action| {
                            action(state, flows, flow_derivative, feedback_input)
                        }) || feedback_input.iter().any(|v| !v.is_finite())
                        {
                            return false;
                        }
                        for (total, value) in input_derivative.iter_mut().zip(feedback_input.iter())
                        {
                            *total += value;
                        }
                    }
                    physical
                        .vector_field_from_effort_into(
                            original_effort,
                            effort,
                            input,
                            input_derivative,
                            flow,
                        )
                        .is_ok()
                },
                EffortObservation::Parameter(&counted),
                &mut budget,
                &mut cancelled,
            );
            if exhausted.get() {
                return Err(PhsPortSidebandError::Budget(
                    "physical parameter-sample cap exhausted",
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
        ) + if direction.structure.ports.is_some() {
            8 * n as u128
        } else {
            0
        };
        let output_effort = RefCell::new(vec![
            0.0;
            if direction.structure.ports.is_some() {
                n
            } else {
                0
            }
        ]);
        let result = observe_sampled_effort_with_budget(
            original,
            base,
            &state_parameter_derivative.coefficients,
            c,
            &c.envelope_tolerances,
            m,
            output_map_work,
            &|_, state, effort, imaginary, output| {
                if imaginary || direction.structure.ports.is_none() {
                    return system.output_from_effort_into(effort, output).is_ok();
                }
                let mut original_effort = output_effort.borrow_mut();
                original_effort.fill(f64::NAN);
                system.storage().gradient(state, &mut original_effort);
                physical
                    .output_from_effort_into(&original_effort, effort, output)
                    .is_ok()
            },
            EffortObservation::StateAndParameter(&counted),
            &mut budget,
            &mut cancelled,
        );
        if exhausted.get() {
            return Err(PhsPortSidebandError::Budget(
                "physical parameter-sample cap exhausted",
            ));
        }
        let (port_parameter_derivative, observation_trace, hessian_actions) = result?;
        poll(&mut cancelled)?;
        Ok(PhsForcedParameterTangent {
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
