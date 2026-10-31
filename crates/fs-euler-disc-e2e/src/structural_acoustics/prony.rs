//! Causal, source-bound scalar Prony memory in the existing Euler modal audio path.
//!
//! The original homogeneous isotropic basis is an EQUILIBRIUM elastic basis.
//! At fixed Poisson ratio, every Maxwell arm multiplies its complete stiffness,
//! so each retained mass-normalized mode keeps its original shape and gains
//! genuine constitutive memory. Storage dispersion and loss therefore enter
//! the same native held-force propagator; no loss-factor-only oscillator is
//! substituted and no frequency-dependent stiffness is frozen at one sample.
//!
//! The retained mechanical subspace, source/fit band and prepared-solver limits
//! remain explicit. A time-domain transient is the response of the supplied
//! Prony law; a finite-band lowering error does not certify out-of-band forcing.
//! Pressure retains the existing fixed modal radiation transfers, hence the
//! original narrow-band observer approximation. This does not turn them into
//! a broadband radiation model or add acoustic backreaction.

use super::{
    ContentHash, DomainHasher, EulerControlStream, GeneralizedForceReconstructionInput,
    IsotropicElasticStatePoint, LoweredModel, ModalAcousticRadiation, ModalAcousticTimeBudget,
    ModalAcousticTimeModel, ModalLossSpectrum, PhysicalContactForceSampling,
    PhysicalModalAudioModel, PhysicalModalInitialState, PhysicalPressureObserver,
    PhysicalPressureSignal, PointForceProjection, ResolvedElasticDiscProfile, StructuralModalBasis,
    StructuralModalBasisError, fixed_rate_frame_count_with_roundoff_bound, hash_pressure_observer,
    modal_force_for_control_interval, reconstruct_projected_modal_forces,
};
use fs_couple::modal_acoustic_time::{ModalAcousticMode, ModalAcousticState};
use fs_couple::render::schedule::force::coupled::{
    CoupledModalFrame, CoupledModalSystem, ModalCouplingConfig, ModalCouplingError,
    PreparedModalBudget, harmonic::ModalHarmonicResponse,
};
use fs_exec::{CancelGate, Cx};
use fs_material::visco::GeneralizedMaxwell;
use fs_math::{c64::C64, det};
use fs_phs::RelaxationBranch;

/// Explicit physical, allocation and preparation limits of the native runtime.
#[derive(Clone, Copy, Debug)]
pub struct StructuralPronyAudioBudget {
    /// Original per-mode state and pressure limits.
    pub modal: ModalAcousticTimeBudget,
    /// Existing complete-system work, storage, pressure and Nyquist admission.
    pub coupling: ModalCouplingConfig,
    /// Existing exact held-force preparation budget.
    pub prepared: PreparedModalBudget,
    /// Maximum total positive Maxwell arms over all retained modes.
    pub maximum_material_branches: usize,
}

/// Complete immutable source/basis/observer binding for causal material audio.
#[derive(Clone, Copy, Debug)]
pub struct StructuralPronyAudioRequest<'a> {
    /// Original equilibrium structural shapes and force projection.
    pub basis: &'a StructuralModalBasis,
    /// Exact specimen from which that basis was assembled.
    pub specimen: &'a ResolvedElasticDiscProfile,
    /// Original isotropic elastic source, declaring the EQUILIBRIUM Young modulus.
    ///
    /// Its source identity, density and full stiffness must match the specimen;
    /// its Young modulus must equal the Prony equilibrium modulus exactly.
    /// A measured storage modulus at a nonzero reference frequency is refused.
    pub equilibrium: &'a IsotropicElasticStatePoint,
    /// Passive source-authored Prony coefficients and their declared fit band.
    pub material: &'a LoweredModel,
    /// Nonzero identity of the constitutive model and its source evidence.
    pub material_model_identity: ContentHash,
    /// Existing identity-bound modal radiation observation.
    pub radiation: &'a ModalAcousticRadiation,
    /// Fixed output sample rate [Hz].
    pub sample_rate_hz: u32,
    /// Caller-owned limits; no modes or positive branches are dropped to fit.
    pub budget: StructuralPronyAudioBudget,
}

/// A physical point-force projection and the accepted native energy/pressure row.
#[derive(Clone, Debug)]
pub struct StructuralPronyPressureFrame {
    /// Original P1 boundary force projection on every retained mode.
    pub force_projection: PointForceProjection,
    /// Pressure, storage, external work and independent relaxation heat [SI].
    pub acoustic: CoupledModalFrame,
}

/// Refusal before publication of a material binding or complete audio step.
#[derive(Debug)]
pub enum StructuralPronyAudioError {
    /// Source, basis, force projection or radiation binding refused.
    Structural(StructuralModalBasisError),
    /// The existing complete native propagator refused.
    Coupling(ModalCouplingError),
    /// Incompatible constitutive reference, source band or projection extent.
    Invalid {
        /// Exact rejected invariant.
        what: &'static str,
    },
}

impl core::fmt::Display for StructuralPronyAudioError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Structural(error) => error.fmt(formatter),
            Self::Coupling(error) => error.fmt(formatter),
            Self::Invalid { what } => write!(formatter, "structural Prony audio: {what}"),
        }
    }
}
impl std::error::Error for StructuralPronyAudioError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Structural(error) => Some(error),
            Self::Coupling(error) => Some(error),
            Self::Invalid { .. } => None,
        }
    }
}
impl From<StructuralModalBasisError> for StructuralPronyAudioError {
    fn from(error: StructuralModalBasisError) -> Self {
        Self::Structural(error)
    }
}
impl From<ModalCouplingError> for StructuralPronyAudioError {
    fn from(error: ModalCouplingError) -> Self {
        Self::Coupling(error)
    }
}

fn invalid(what: &'static str) -> StructuralPronyAudioError {
    StructuralPronyAudioError::Invalid { what }
}

/// Point-force-driven material audio using the existing Maxwell memory owner.
///
/// The stored request enables a transactional reset to zero mechanical state
/// and fully relaxed material memory. Accepted states and material histories
/// are separately inspectable; neither is silently reset during stepping.
pub struct StructuralPronyAudioModel<'a> {
    request: StructuralPronyAudioRequest<'a>,
    runtime: CoupledModalSystem,
    identity: ContentHash,
}

impl<'a> StructuralPronyAudioModel<'a> {
    /// Bind the original equilibrium shapes and observer to causal Prony memory.
    ///
    /// Every positive arm is retained. The entire interval from each equilibrium
    /// frequency to its instantaneous-stiffness upper bound must lie in the
    /// supplied material band. The native owner separately enforces its
    /// conservative Nyquist, work, state and energy bounds.
    pub fn try_new(
        request: StructuralPronyAudioRequest<'a>,
        gate: &CancelGate,
    ) -> Result<Self, StructuralPronyAudioError> {
        if gate.is_requested() {
            return Err(ModalCouplingError::Cancelled.into());
        }
        validate_source(&request)?;
        let runtime = prepare_runtime(&request, None, gate)?;
        let identity = binding_identity(&request);
        Ok(Self {
            request,
            runtime,
            identity,
        })
    }

    /// Identity of all constitutive coefficients, source bindings and observer.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Original mass-normalized mechanical states in structural-mode order.
    #[must_use]
    pub fn states(&self) -> &[ModalAcousticState] {
        self.runtime.components()[0].states()
    }

    /// Energy-normalized Maxwell histories, ordered by mode then source arm.
    #[must_use]
    pub fn material_memory(&self) -> &[f64] {
        self.runtime.material_memory()
    }

    /// Completed samples only; refused trials consume no time or history.
    #[must_use]
    pub const fn samples_rendered(&self) -> u64 {
        self.runtime.samples_rendered()
    }

    /// Last accepted complete native pressure/energy row.
    #[must_use]
    pub fn last_frame(&self) -> Option<&CoupledModalFrame> {
        self.runtime.last_frame()
    }

    /// Return to zero mechanical state and fully relaxed material memory.
    ///
    /// Rebuilds from the retained immutable sources. Cancellation or a failed
    /// preparation leaves the currently accepted runtime intact.
    pub fn reset(&mut self, gate: &CancelGate) -> Result<(), StructuralPronyAudioError> {
        let replacement = Self::try_new(self.request, gate)?;
        *self = replacement;
        Ok(())
    }

    /// Initialize a cold model at fully RELAXED static equilibrium.
    ///
    /// Mechanical displacement uses the declared equilibrium stiffness. Each
    /// Maxwell history is then initialized at that displacement, so a constant
    /// preload creates neither unrelaxed branch stress nor an artificial ring.
    /// Preparation is transactional; use `reset` before reinitializing a model
    /// that has already accepted a sample.
    pub fn initialize_static_equilibrium(
        &mut self,
        force_n_per_sqrt_kg: &[f64],
        gate: &CancelGate,
    ) -> Result<(), StructuralPronyAudioError> {
        if self.samples_rendered() != 0 {
            return Err(invalid(
                "static initialization requires a cold runtime; reset first",
            ));
        }
        let replacement = prepare_runtime(&self.request, Some(force_n_per_sqrt_kg), gate)?;
        self.runtime = replacement;
        Ok(())
    }

    /// Synthesize an admitted Euler contact-control horizon through Prony memory.
    ///
    /// Reuses the existing closing-else-opening contact projection and the same
    /// band-limited reconstruction of interval force measures as physical modal
    /// audio. The existing reconstruction owns the 48 kHz master clock, filter,
    /// source-bandwidth admission and horizon padding; a different runtime
    /// sample rate refuses. No interval is treated as a new point impulse.
    ///
    /// This method starts a cold source horizon. `Zero` requires zero accepted
    /// state/history; `StaticEquilibriumAtFirstHeldForce` initializes both the
    /// equilibrium coordinates and relaxed memory under the first reconstructed
    /// force. Static initialization retains the established narrow-band observer
    /// subtraction about the CURRENT held force; it is not broadband radiation.
    ///
    /// Timeline, source, force reconstruction and initialization refuse before
    /// advancing time. Once stepping begins, each complete accepted native frame
    /// remains committed if a later frame or cancellation refuses; `reset`
    /// permits deterministic replay of the complete horizon.
    pub fn synthesize_control_stream(
        &mut self,
        controls: &EulerControlStream<'_>,
        initial_state: PhysicalModalInitialState,
        force_reconstruction: GeneralizedForceReconstructionInput,
        maximum_contact_distance_m: f64,
        gate: &CancelGate,
        cx: &Cx<'_>,
    ) -> Result<PhysicalPressureSignal, StructuralPronyAudioError> {
        if gate.is_requested() {
            return Err(ModalCouplingError::Cancelled.into());
        }
        cx.checkpoint()
            .map_err(|_| StructuralModalBasisError::Cancelled)?;
        if self.samples_rendered() != 0 {
            return Err(invalid(
                "control-stream synthesis requires a cold runtime; reset first",
            ));
        }
        if initial_state == PhysicalModalInitialState::Zero
            && (self.states().iter().any(|state| {
                state.displacement_m_sqrt_kg != 0.0 || state.velocity_m_sqrt_kg_per_s != 0.0
            }) || self.material_memory().iter().any(|memory| *memory != 0.0))
        {
            return Err(invalid(
                "zero initialization requires zero mechanical state and memory; reset first",
            ));
        }
        if controls.source().metadata().specimen_profile_identity
            != self.request.basis.profile_identity
        {
            return Err(StructuralModalBasisError::IdentityMismatch {
                what: "control trajectory specimen profile does not match structural basis profile",
            }
            .into());
        }
        let intervals = controls.audio();
        let first = intervals
            .first()
            .ok_or(StructuralModalBasisError::ControlTimeline {
                what: "control stream has no positive-duration audio intervals",
            })?;
        let last = intervals
            .last()
            .expect("nonempty interval slice has a last item");
        for pair in intervals.windows(2) {
            if pair[0].end_time_s.to_bits() != pair[1].start_time_s.to_bits() {
                return Err(StructuralModalBasisError::ControlTimeline {
                    what: "mechanics audio intervals are not exactly contiguous",
                }
                .into());
            }
        }
        let frame_count = fixed_rate_frame_count_with_roundoff_bound(
            first.start_time_s, last.end_time_s, self.request.sample_rate_hz,
            force_reconstruction.clock_roundoff_operation_count,
        ).ok_or(StructuralModalBasisError::ControlTimeline {
            what: "control horizon is not a positive integral number of audio samples within declared clock roundoff",
        })?;
        let mut modal_forces = Vec::with_capacity(intervals.len());
        for interval in intervals {
            if gate.is_requested() {
                return Err(ModalCouplingError::Cancelled.into());
            }
            cx.checkpoint()
                .map_err(|_| StructuralModalBasisError::Cancelled)?;
            modal_forces.push(modal_force_for_control_interval(
                &self.request.basis.mesh,
                &self.request.basis.modes,
                controls,
                interval,
                maximum_contact_distance_m,
            )?);
        }
        let reconstructed =
            reconstruct_projected_modal_forces(intervals, &modal_forces, force_reconstruction, cx)?;
        if reconstructed.sample_rate_hz != self.request.sample_rate_hz
            || reconstructed.frame_count() != frame_count
        {
            return Err(StructuralModalBasisError::ControlTimeline {
                what: "reconstructed specimen force clock does not match modal audio clock",
            }
            .into());
        }
        let mut pressure_pa = Vec::new();
        pressure_pa.try_reserve_exact(frame_count).map_err(|_| {
            StructuralModalBasisError::PressureCapacity {
                requested: frame_count,
            }
        })?;
        let static_transfers: Vec<_> = self
            .request
            .radiation
            .modes
            .iter()
            .map(|mode| mode.observer_pressure_per_modal_velocity)
            .collect();
        if initial_state == PhysicalModalInitialState::StaticEquilibriumAtFirstHeldForce {
            self.initialize_static_equilibrium(
                reconstructed
                    .frame(0)
                    .expect("positive frame count has a first force row"),
                gate,
            )?;
        }
        let mut peak_abs_pressure_pa = 0.0_f64;
        for index in 0..frame_count {
            if index % 64 == 0 {
                cx.checkpoint()
                    .map_err(|_| StructuralModalBasisError::Cancelled)?;
            }
            let force = reconstructed
                .frame(index)
                .expect("reconstructed clock matches the complete pressure horizon");
            let native_pressure = self
                .runtime
                .step_under_gate(force, gate)?
                .observer_pressure_pa;
            let pressure =
                if initial_state == PhysicalModalInitialState::StaticEquilibriumAtFirstHeldForce {
                    self.runtime.components()[0]
                        .observer_pressure_with_transfers_about_static_equilibrium(
                            &static_transfers,
                            force,
                        )
                        .map_err(StructuralModalBasisError::from)?
                } else {
                    native_pressure
                };
            peak_abs_pressure_pa = peak_abs_pressure_pa.max(pressure.abs());
            pressure_pa.push(pressure);
        }
        if gate.is_requested() {
            return Err(ModalCouplingError::Cancelled.into());
        }
        cx.checkpoint()
            .map_err(|_| StructuralModalBasisError::Cancelled)?;
        let observer = PhysicalPressureObserver::BodyFixed(self.request.radiation.observer);
        let mut identity =
            DomainHasher::new("org.frankensim.euler-disc.prony-physical-pressure-signal.v1");
        identity.update(self.identity.as_bytes());
        identity.update(
            controls
                .source()
                .metadata()
                .configuration_identity
                .as_bytes(),
        );
        identity.update(reconstructed.identity.as_bytes());
        identity.update(&self.request.sample_rate_hz.to_le_bytes());
        identity.update(&first.start_time_s.to_bits().to_le_bytes());
        identity.update(&[match initial_state {
            PhysicalModalInitialState::Zero => 0,
            PhysicalModalInitialState::StaticEquilibriumAtFirstHeldForce => 1,
        }]);
        hash_pressure_observer(&mut identity, observer);
        identity.update(&(pressure_pa.len() as u64).to_le_bytes());
        for pressure in &pressure_pa {
            identity.update(&pressure.to_bits().to_le_bytes());
        }
        Ok(PhysicalPressureSignal {
            start_time_s: first.start_time_s,
            sample_rate_hz: self.request.sample_rate_hz,
            pressure_pa,
            peak_abs_pressure_pa,
            contact_force_sampling:
                PhysicalContactForceSampling::IntervalMeasureAtClosingElseOpeningEndpointBandLimitedV1,
            observer,
            structural_basis_identity: self.request.basis.identity,
            radiation_identity: self.request.radiation.identity,
            damping_model_identity: self.request.material_model_identity,
            identity: identity.finalize(),
        })
    }

    /// Advance one actual mass-normalized held-force sample [N/sqrt(kg)].
    ///
    /// The native owner publishes modal states, material histories, pressure
    /// and independent work/heat diagnostics together after all gates pass.
    pub fn step_modal_force(
        &mut self,
        force_n_per_sqrt_kg: &[f64],
        gate: &CancelGate,
    ) -> Result<&CoupledModalFrame, StructuralPronyAudioError> {
        Ok(self.runtime.step_under_gate(force_n_per_sqrt_kg, gate)?)
    }

    /// Project an actual body-frame contact force and advance one audio sample.
    pub fn step_point_force(
        &mut self,
        point_body_m: [f64; 3],
        force_body_n: [f64; 3],
        maximum_distance_m: f64,
        gate: &CancelGate,
    ) -> Result<StructuralPronyPressureFrame, StructuralPronyAudioError> {
        if gate.is_requested() {
            return Err(ModalCouplingError::Cancelled.into());
        }
        let force_projection = self.request.basis.project_point_force(
            point_body_m,
            force_body_n,
            maximum_distance_m,
        )?;
        let acoustic = self
            .runtime
            .step_under_gate(&force_projection.modal_force_n_per_sqrt_kg, gate)?
            .clone();
        Ok(StructuralPronyPressureFrame {
            force_projection,
            acoustic,
        })
    }

    /// Continuous response of the SAME full memory system, in exp(-i omega t).
    ///
    /// This read-only query keeps all branches and checks the material band at
    /// the actual drive frequency. It does not replace the time propagator with
    /// a loss-factor fit or alter accepted states.
    pub fn harmonic_response(
        &self,
        omega_rad_s: f64,
        peak_force_n_per_sqrt_kg: &[C64],
        maximum_work: usize,
        gate: &CancelGate,
    ) -> Result<ModalHarmonicResponse, StructuralPronyAudioError> {
        self.request
            .material
            .modulus_checked(omega_rad_s)
            .map_err(StructuralModalBasisError::from)?;
        Ok(self.runtime.harmonic_response(
            omega_rad_s,
            peak_force_n_per_sqrt_kg,
            maximum_work,
            gate,
        )?)
    }
}

fn validate_source(
    request: &StructuralPronyAudioRequest<'_>,
) -> Result<(), StructuralPronyAudioError> {
    let reference = fs_solid::TetElasticMaterial::from_resolved_elastic_state(request.equilibrium);
    if request.basis.specimen_identity != request.specimen.identity
        || request.basis.material_state_identity != request.specimen.material_state_identity
        || request.equilibrium.identity() != request.specimen.material_state_identity
        || request.equilibrium.resolved().card_identity() != request.specimen.material_card_identity
        || reference.material_state_identity
            != request.specimen.elastic_material.material_state_identity
        || reference.density_kg_m3.to_bits()
            != request.specimen.elastic_material.density_kg_m3.to_bits()
        || reference.stiffness_mandel_pa()
            != request.specimen.elastic_material.stiffness_mandel_pa()
    {
        return Err(invalid(
            "Prony equilibrium source, specimen and structural basis must be identical",
        ));
    }
    if request.material.model.e_inf.to_bits() != request.equilibrium.young_modulus_pa().to_bits() {
        return Err(invalid(
            "Prony equilibrium modulus must equal the original equilibrium Young modulus",
        ));
    }
    if request.material_model_identity == ContentHash([0; 32]) {
        return Err(invalid("material_model_identity must not be zero"));
    }
    Ok(())
}

fn prepare_runtime(
    request: &StructuralPronyAudioRequest<'_>,
    static_force: Option<&[f64]>,
    gate: &CancelGate,
) -> Result<CoupledModalSystem, StructuralPronyAudioError> {
    if gate.is_requested() {
        return Err(ModalCouplingError::Cancelled.into());
    }
    let loss = ModalLossSpectrum {
        structural_basis_identity: request.basis.identity,
        material_state_identity: request.basis.material_state_identity,
        damping_model_identity: request.material_model_identity,
        // Intrinsic loss belongs to native material memory exactly once.
        loss_factors: vec![0.0; request.basis.modes.len()],
    };
    let mut elastic = PhysicalModalAudioModel::try_new(
        request.basis,
        &loss,
        request.radiation,
        request.sample_rate_hz,
        request.budget.modal,
    )?;
    if let Some(force) = static_force {
        elastic
            .runtime
            .initialize_static_equilibrium(force)
            .map_err(StructuralModalBasisError::from)?;
    }
    // Attach AFTER setting equilibrium coordinates. None initializes relaxed
    // memory z_j = sqrt(k_j) q, not zero histories with spurious branch stress.
    bind_prony_memory(elastic.runtime, request.material, request.budget, gate)
}

fn bind_prony_memory(
    elastic: ModalAcousticTimeModel,
    material: &LoweredModel,
    budget: StructuralPronyAudioBudget,
    gate: &CancelGate,
) -> Result<CoupledModalSystem, StructuralPronyAudioError> {
    let arms = project_prony_arms(elastic.modes(), material, budget, gate)?;
    let system = CoupledModalSystem::new_prepared(
        vec![elastic],
        Vec::new(),
        budget.coupling,
        budget.prepared,
        gate,
    )?;
    Ok(system.with_material_relaxation(arms, None, gate)?)
}

fn project_prony_arms(
    modes: &[ModalAcousticMode],
    material: &LoweredModel,
    budget: StructuralPronyAudioBudget,
    gate: &CancelGate,
) -> Result<Vec<RelaxationBranch>, StructuralPronyAudioError> {
    if gate.is_requested() {
        return Err(ModalCouplingError::Cancelled.into());
    }
    let n = modes.len();
    let positive = material
        .model
        .terms
        .iter()
        .filter(|(modulus, _)| *modulus > 0.0)
        .count();
    let branches = n
        .checked_mul(positive)
        .ok_or_else(|| invalid("material branch extent overflow"))?;
    let dimension = n
        .checked_mul(2)
        .ok_or_else(|| invalid("material projection dimension overflow"))?;
    let entries = branches
        .checked_mul(dimension)
        .ok_or_else(|| invalid("material projection storage overflow"))?;
    if n == 0
        || n > budget.coupling.max_modes
        || n > budget.prepared.max_modes
        || branches > budget.maximum_material_branches
        || entries > budget.coupling.max_setup_terms
        || material.model.terms.len() > budget.coupling.max_setup_terms
    {
        return Err(invalid(
            "complete Prony projection exceeds the authored mode, branch or setup budget",
        ));
    }
    let law = GeneralizedMaxwell::new(material.model.e_inf, material.model.terms.clone())
        .map_err(StructuralModalBasisError::from)?;
    if !material.band.0.is_finite()
        || !material.band.1.is_finite()
        || material.band.0 < 0.0
        || material.band.1 <= material.band.0
        || [
            material.sup_rel_err,
            material.sup_rel_storage_err,
            material.sup_rel_loss_err,
        ]
        .iter()
        .any(|error| !error.is_finite() || *error < 0.0)
        || material.verification_points < 2
    {
        return Err(invalid(
            "Prony fit band and measured error metadata must be finite and admissible",
        ));
    }
    let mut instantaneous_ratio = 1.0;
    for &(modulus, _) in &law.terms {
        instantaneous_ratio += modulus / law.e_inf;
    }
    if !instantaneous_ratio.is_finite() || instantaneous_ratio < 1.0 {
        return Err(invalid("instantaneous modulus is not representable"));
    }
    let root_ratio = det::sqrt(instantaneous_ratio);
    let mut arms = Vec::with_capacity(branches);
    for (index, mode) in modes.iter().enumerate() {
        if gate.is_requested() {
            return Err(ModalCouplingError::Cancelled.into());
        }
        if mode.damping_ratio != 0.0 {
            return Err(invalid(
                "intrinsic Prony loss cannot be added to a separately damped base mode",
            ));
        }
        let omega = mode.angular_frequency_rad_s;
        let upper = omega * root_ratio;
        if !upper.is_finite() || omega < material.band.0 || upper > material.band.1 {
            return Err(invalid(
                "equilibrium through instantaneous modal frequencies must lie in the material band",
            ));
        }
        for &(modulus, relaxation_time_s) in &law.terms {
            if modulus == 0.0 {
                continue;
            }
            // q has mass-normalized length units; k_j*q² has units of energy.
            // Scale by the dimensionless modulus ratio before the second omega
            // to avoid an unnecessary square of a dimensional modulus.
            let stiffness = (omega * (modulus / law.e_inf)) * omega;
            if !stiffness.is_finite() || stiffness <= 0.0 {
                return Err(invalid(
                    "positive Prony modal stiffness is not representable",
                ));
            }
            let mut projection = vec![0.0; dimension];
            projection[2 * index] = 1.0;
            arms.push(RelaxationBranch {
                projection,
                stiffness,
                relaxation_time_s,
            });
        }
    }
    Ok(arms)
}

fn binding_identity(request: &StructuralPronyAudioRequest<'_>) -> ContentHash {
    let mut identity = DomainHasher::new("org.frankensim.euler-disc.structural-prony-audio.v1");
    for hash in [
        request.basis.identity,
        request.specimen.identity,
        request.equilibrium.identity(),
        request.material_model_identity,
        request.radiation.identity,
    ] {
        identity.update(hash.as_bytes());
    }
    identity.update(&request.sample_rate_hz.to_le_bytes());
    for value in [
        request.material.model.e_inf,
        request.material.band.0,
        request.material.band.1,
        request.material.sup_rel_err,
        request.material.sup_rel_storage_err,
        request.material.sup_rel_loss_err,
    ] {
        identity.update(&value.to_bits().to_le_bytes());
    }
    identity.update(&(request.material.verification_points as u64).to_le_bytes());
    identity.update(&(request.material.model.terms.len() as u64).to_le_bytes());
    for &(modulus, time) in &request.material.model.terms {
        identity.update(&modulus.to_bits().to_le_bytes());
        identity.update(&time.to_bits().to_le_bytes());
    }
    identity.finalize()
}

#[cfg(test)]
mod tests;
