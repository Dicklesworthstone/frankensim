//! Nonlinear storage coupled to the SAME exact linear and friction owners.
//!
//! Each trial supplies a held reaction -DG(V,q0,q1), plus any explicitly bound
//! passive constitutive strain-rate force, to the original complete
//! linear network. The original scalar or simultaneous bow solve includes that
//! reaction. The fixed point is accepted only when the reaction law AND the combined energy
//! balance close. fs-phs owns DG; this adapter changes neither the exponential
//! linear propagator nor the friction law, its stick/slip branch rule or work
//! port. Only a step whose reaction iteration exhausts its budget is solved
//! again, with each port's branch memory and any joint rescue's stick/slip
//! partition held across trials (`stage_trial`).
//! This is a bounded iterated discrete-gradient composition, not SAV, a unique
//! nonlinear-solution claim, a trajectory certificate or a real-time guarantee.
use super::*;
use crate::stribeck_friction::ContactBranch;
use fs_phs::Storage;
use std::borrow::Cow;

mod equilibrium;
pub use equilibrium::NonlinearPreloadDiagnostics;
pub mod history;
use history::BoundHistory;

/// Displacement-only positive constitutive energy. Contiguous bindings supply
/// additional energy; mapped bindings may explicitly replace selected original
/// diagonal stiffness. Inertia remains with the original linear owner. The
/// constitutive owner admits its own physical validity envelope.
pub trait ModalPotential: Storage + Send + Sync {
    /// Number of original mass-normalized displacement coordinates.
    fn coordinate_count(&self) -> usize;
    /// Reject nonfinite coordinates or motion outside the constitutive domain.
    fn validate_coordinates(&self, q: &[f64]) -> Result<(), ModalCouplingError>;
    /// Whether endpoint admission additionally needs actual modal velocities.
    /// This immutable selection prepares reusable mapped state buffers once;
    /// displacement-only laws incur no velocity gathers. Wrappers must forward
    /// this selection together with `validate_state`.
    fn requires_velocity_validation(&self) -> bool { false }
    /// Admit a complete native endpoint without changing storage or history.
    /// Override together with `requires_velocity_validation` for rate-dependent
    /// validity. The host supplies the original accepted or staged velocities,
    /// never displacement increment / dt from its held-force loss quadrature.
    /// Admission runs at cold binding and before dynamic/static publication.
    fn validate_state(&self, q: &[f64], velocity: &[f64]) -> Result<(), ModalCouplingError> {
        if velocity.len() != self.coordinate_count() || velocity.iter().any(|v| !v.is_finite()) {
            return Err(invalid("constitutive admission needs complete finite modal velocities"));
        }
        self.validate_coordinates(q)
    }
    /// Whether the constitutive reaction below can be nonzero. Elastic laws
    /// keep the original solve path without additional per-trial loss work.
    fn has_constitutive_loss(&self) -> bool { false }
    /// Optional memoryless constitutive reaction in this SAME coordinate basis.
    /// It must be passive: `reaction dot velocity <= 0`. The default is elastic.
    ///
    /// The coupled solve evaluates this at midpoint coordinates and actual
    /// displacement increment / dt. Its held-force work is included once in
    /// network dissipation. Endpoint acceleration observation instead supplies
    /// endpoint coordinates/velocities. A law must not repeat modal or joint
    /// damping already owned by the linear network. No state may be advanced.
    fn dissipative_reaction(
        &self, _q: &[f64], _velocity: &[f64], reaction: &mut [f64],
    ) -> Result<(), ModalCouplingError> {
        reaction.fill(0.0);
        Ok(())
    }
}

/// A contiguous, nonoverlapping part of one original component.
pub struct ModalPotentialBinding {
    /// Original component index, without reordering its states.
    pub component: usize,
    /// First coordinate within that component.
    pub coordinate_start: usize,
    /// Additional conservative storage in the unchanged normalized basis.
    pub potential: Box<dyn ModalPotential>,
}

/// A coordinate in an existing component's mass-normalized basis.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ModalCoordinate {
    pub component: usize,
    pub mode: usize,
}

/// Conservative storage coupling coordinates of several original components.
/// The coordinate list is the constitutive owner's order, not a new basis.
pub struct MappedModalPotentialBinding {
    pub coordinates: Vec<ModalCoordinate>,
    /// Local indices whose original diagonal elastic energy is INCLUDED in
    /// `potential`. That energy and its midpoint force are subtracted from the
    /// nonlinear correction, since the exact linear propagator retains them.
    /// This allows a positive full strain energy with a signed coupling
    /// remainder, without counting linear stiffness twice.
    pub replaces_linear: Vec<usize>,
    pub potential: Box<dyn ModalPotential>,
}

/// Bounded constitutive iteration; nonconvergence requires physical clock/modal
/// refinement, not clipping a reaction or publishing a partial solve.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NonlinearModalConfig {
    /// Maximum complete linear/friction trials, in 1..=128. Static preload
    /// uses the same cap for constitutive evaluations, including backtracking.
    pub max_iterations: usize,
    /// Maximum generalized reaction magnitude [N/sqrt(kg)].
    pub maximum_generalized_force: f64,
    /// Absolute reaction-law residual allowance [N/sqrt(kg)].
    pub force_absolute_tolerance: f64,
    /// Relative reaction-law residual allowance, in (0,1).
    pub force_relative_tolerance: f64,
}
impl NonlinearModalConfig {
    /// Admit finite positive force budgets and a bounded trial count.
    pub fn validate(self) -> Result<(), ModalCouplingError> {
        if !(1..=128).contains(&self.max_iterations)
            || [
                self.maximum_generalized_force,
                self.force_absolute_tolerance,
                self.force_relative_tolerance,
            ]
            .iter()
            .any(|v| !v.is_finite() || *v <= 0.0)
            || self.force_relative_tolerance >= 1.0
        {
            return Err(invalid(
                "nonlinear storage requires finite force/iteration budgets",
            ));
        }
        Ok(())
    }
}

/// One jointly accepted constitutive, friction and complete-network sample.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NonlinearFrictionFrame {
    /// Network energy includes nonlinear storage; external work excludes the
    /// internal geometric reaction.
    pub friction: FrictionModalFrame,
    /// Conservative correction [J], already counted in friction energy. May be
    /// signed when full strain storage replaces selected linear stiffness.
    pub nonlinear_energy_j: f64,
    /// Constitutive strain-rate loss [J], already in network dissipation.
    pub nonlinear_dissipation_j: f64,
    /// Complete linear/friction trials required by this accepted sample.
    pub nonlinear_iterations: usize,
    /// Maximum constitutive force mismatch [N/sqrt(kg)].
    pub reaction_residual: f64,
    /// Absolute-plus-relative allowed mismatch [N/sqrt(kg)].
    pub reaction_tolerance: f64,
}

/// Constitutive report shared by scalar and simultaneous friction.
/// The additional storage is already included in the accepted network energy.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NonlinearFrictionDiagnostics {
    /// Conservative correction [J]; may be signed for linear replacements.
    pub nonlinear_energy_j: f64,
    /// Constitutive strain-rate loss [J], already in network dissipation.
    pub nonlinear_dissipation_j: f64,
    /// Complete network/friction trials required by this accepted sample.
    pub nonlinear_iterations: usize,
    /// Maximum constitutive force mismatch [N/sqrt(kg)].
    pub reaction_residual: f64,
    /// Absolute-plus-relative allowed mismatch [N/sqrt(kg)].
    pub reaction_tolerance: f64,
}
impl NonlinearFrictionFrame {
    pub(super) fn diagnostics(self) -> NonlinearFrictionDiagnostics {
        NonlinearFrictionDiagnostics {
            nonlinear_energy_j: self.nonlinear_energy_j,
            nonlinear_dissipation_j: self.nonlinear_dissipation_j,
            nonlinear_iterations: self.nonlinear_iterations,
            reaction_residual: self.reaction_residual,
            reaction_tolerance: self.reaction_tolerance,
        }
    }
}

// Both adapters stage the ORIGINAL complete network, with no publication. The
// constitutive iteration below owns neither a friction law nor a propagator.
pub(super) trait NonlinearTrial {
    fn core(&self) -> &FrictionModalSystem;
    fn stage(&mut self, forces: &[f64], gate: Option<&CancelGate>)
        -> Result<NonlinearTrialWork, ModalFrictionError>;
    /// Each dry port's McIntyre-Schumacher-Woodhouse memory, which its root
    /// rule starts from: the accepted step's, except inside the held-branch
    /// pass of `stage_trial`, which sets the previous trial's and restores it.
    fn branch_memory(&self) -> Vec<ContactBranch>;
    fn set_branch_memory(&mut self, memory: &[ContactBranch]);
    /// Each port's branch as the last `stage` resolved it, not yet committed.
    fn staged_branches(&self) -> Vec<ContactBranch>;
    /// Held pass only: keeps the joint stick/slip partition the last `stage`
    /// resolved, for the next stage's active-set rescue to try first. Returns,
    /// per port, whether it differs from the partition held before (none differ
    /// on the first call). A trial without a joint rescue holds nothing.
    fn hold_partition(&mut self) -> Vec<bool> { Vec::new() }
    /// Releases the held partition, so later stages rescue as before.
    fn release_partition(&mut self) {}
}
pub(super) struct NonlinearTrialWork {
    pub observer_pressure_pa: f64,
    pub network_energy_j: f64,
    pub drive_work_j: f64,
    pub absolute_drive_work_j: f64,
    pub network_dissipation_j: f64,
    pub friction_dissipation_j: f64,
    /// Signed F*d of the endpoint stick rule not certified as kinetic loss.
    pub declared_slip_work_j: f64,
    pub internal_storage: Option<NonlinearInternalStorage>,
    pub contacts: Option<NonlinearContactWork>,
}
/// Original contact mechanisms may store energy outside the modal coordinates.
/// Their dissipation and supplied work already belong to the trial's friction
/// accounts; only accepted/staged storage must be added to the joint DG balance.
#[derive(Clone, Copy)]
pub(super) struct NonlinearInternalStorage {
    pub before_energy_j: f64,
    pub after_energy_j: f64,
}
/// Already evaluated by the original normal-contact owner at this same trial.
/// For x=N^T q-s, -R N^T delta_q = -delta_U-loss-R delta_s. Thus normal
/// actuator work belongs in the complete balance without changing the DG law.
#[derive(Clone, Copy)]
pub(super) struct NonlinearContactWork {
    pub before_energy_j: f64,
    pub after_energy_j: f64,
    pub dissipation_j: f64,
    pub work_j: f64,
    pub absolute_work_j: f64,
}
pub(super) struct NonlinearCorrection {
    pub observer_pressure_pa: f64,
    pub network_energy_j: f64,
    pub network_dissipation_j: f64,
    pub external_work_j: f64,
    pub energy_residual_j: f64,
    pub energy_tolerance_j: f64,
    pub diagnostics: NonlinearFrictionDiagnostics,
}
struct SingleTrial<'a> {
    inner: &'a mut FrictionModalSystem,
    drive: FrictionDrive,
    frame: Option<FrictionModalFrame>,
}
impl NonlinearTrial for SingleTrial<'_> {
    fn core(&self) -> &FrictionModalSystem { self.inner }
    fn stage(&mut self, forces: &[f64], gate: Option<&CancelGate>)
        -> Result<NonlinearTrialWork, ModalFrictionError>
    {
        let frame = self.inner.stage_inner(forces, self.drive, gate)?;
        self.frame = Some(frame);
        Ok(NonlinearTrialWork {
            observer_pressure_pa: frame.observer_pressure_pa,
            network_energy_j: frame.network_energy_j,
            drive_work_j: frame.drive_work_j,
            absolute_drive_work_j: frame.drive_work_j.abs(),
            network_dissipation_j: frame.network_dissipation_j,
            friction_dissipation_j: frame.friction_dissipation_j,
            declared_slip_work_j: frame.declared_slip_work_j,
            internal_storage: None,
            contacts: None,
        })
    }
    fn branch_memory(&self) -> Vec<ContactBranch> { vec![self.inner.branch] }
    fn set_branch_memory(&mut self, memory: &[ContactBranch]) { self.inner.branch = memory[0]; }
    fn staged_branches(&self) -> Vec<ContactBranch> { vec![self.inner.staged_branch] }
}
/// A reaction iteration that used its whole budget without being accepted.
struct Exhausted {
    ratio: f64,
    energy: (f64, f64),
}
struct BoundPotential {
    coordinates: Vec<usize>,
    contiguous: Option<std::ops::Range<usize>>,
    replaced_stiffness: Vec<f64>,
    state_validation: Option<PotentialStateValidation>,
    potential: Box<dyn ModalPotential>,
}
// Original component/mode addresses survive cold extensions of the network.
// Both mapped buffers are allocated once, including noncontiguous bindings.
struct PotentialStateValidation {
    coordinates: Vec<ModalCoordinate>,
    q: Vec<f64>,
    velocity: Vec<f64>,
}
impl PotentialStateValidation {
    fn gather(&mut self, models: &[ModalAcousticTimeModel]) -> Result<(), ModalCouplingError> {
        for (index, address) in self.coordinates.iter().enumerate() {
            let state = models.get(address.component)
                .and_then(|model| model.states().get(address.mode))
                .ok_or_else(|| invalid("constitutive admission lost its original modal state"))?;
            self.q[index] = state.displacement_m_sqrt_kg;
            self.velocity[index] = state.velocity_m_sqrt_kg_per_s;
        }
        Ok(())
    }
}
impl BoundPotential {
    fn select<'a>(&self, q: &'a [f64]) -> Cow<'a, [f64]> {
        match &self.contiguous {
            Some(range) => Cow::Borrowed(&q[range.clone()]),
            None => Cow::Owned(self.coordinates.iter().map(|&index| q[index]).collect()),
        }
    }
}
pub(super) struct NonlinearFrictionState {
    potentials: Vec<BoundPotential>,
    histories: Vec<BoundHistory>,
    staged_history_sample: Option<u64>,
    config: NonlinearModalConfig,
    last: Option<NonlinearFrictionFrame>,
    pub(super) preload: Option<NonlinearPreloadDiagnostics>,
}
impl FrictionModalSystem {
    /// Bind disjoint additional storage before playback. Existing nonzero initial
    /// motion is preserved. Rendering, gestures, interfaces and observation use
    /// this same owner; refused trials publish no accepted state or diagnostic.
    pub fn with_nonlinear_potentials(
        self,
        bindings: Vec<ModalPotentialBinding>,
        config: NonlinearModalConfig,
        gate: &CancelGate,
    ) -> Result<Self, ModalFrictionError> {
        let mut mapped = Vec::with_capacity(bindings.len());
        for binding in bindings {
            poll(Some(gate))?;
            let count = binding.potential.coordinate_count();
            let component = self.components().get(binding.component)
                .ok_or_else(|| invalid("nonlinear potential names an unknown original component"))?;
            let end = binding.coordinate_start.checked_add(count)
                .ok_or_else(|| invalid("nonlinear component coordinate overflow"))?;
            if count == 0 || end > component.modes().len() {
                return Err(invalid("nonlinear potential must fit its unchanged component basis").into());
            }
            mapped.push(MappedModalPotentialBinding {
                coordinates: (binding.coordinate_start..end).map(|mode| ModalCoordinate {
                    component: binding.component, mode,
                }).collect(),
                replaces_linear: Vec::new(),
                potential: binding.potential,
            });
        }
        self.with_mapped_nonlinear_potentials(mapped, config, gate)
    }

    /// Bind a complete positive strain energy across existing components.
    /// Referenced linear oscillators remain in the exact propagator. Only the
    /// constitutive correction enters the joint friction/contact iteration.
    /// Disjoint cold constituents may join an existing binding under the SAME
    /// solve configuration. Every original storage law remains active in one
    /// combined constitutive iteration. Duplicate coordinates, overlapping
    /// bindings, changed numerical budgets, invalid replacement indices or
    /// nonlinear energy outside the complete network budget refuse.
    pub fn with_mapped_nonlinear_potentials(
        mut self,
        bindings: Vec<MappedModalPotentialBinding>,
        config: NonlinearModalConfig,
        gate: &CancelGate,
    ) -> Result<Self, ModalFrictionError> {
        poll(Some(gate))?;
        config.validate()?;
        if self.samples_rendered() != 0
            || self.nonlinear.as_ref().is_some_and(|state| state.config != config)
            || bindings.is_empty()
            || bindings.len() > self.mode_count()
        {
            return Err(invalid("nonlinear storage requires cold nonempty binding and finite force/iteration budgets").into());
        }
        let mut occupied = vec![false; self.mode_count()];
        let (mut potentials,histories) = self.nonlinear.take()
            .map_or_else(||(Vec::new(),Vec::new()), |state| (state.potentials,state.histories));
        for potential in &potentials {
            for &coordinate in &potential.coordinates { occupied[coordinate] = true; }
        }
        for binding in bindings {
            poll(Some(gate))?;
            let count = binding.potential.coordinate_count();
            if count == 0 || count != binding.coordinates.len() || count > self.mode_count() {
                return Err(invalid("nonlinear coordinate map must match its full constitutive basis").into());
            }
            let mut coordinates = Vec::with_capacity(count);
            let mut replaced_stiffness = vec![0.0; count];
            let mut replaced = vec![false; count];
            for &index in &binding.replaces_linear {
                if index >= count || replaced[index] {
                    return Err(invalid("linear replacement indices must be unique constitutive coordinates").into());
                }
                replaced[index] = true;
            }
            for (local, coordinate) in binding.coordinates.iter().enumerate() {
                let mode = self.components().get(coordinate.component)
                    .and_then(|component| component.modes().get(coordinate.mode))
                    .ok_or_else(|| invalid("nonlinear map names an unknown original modal coordinate"))?;
                let global = self.network.offsets[coordinate.component] + coordinate.mode;
                if occupied[global] {
                    return Err(invalid("overlapping additive nonlinear bindings are not admitted").into());
                }
                occupied[global] = true;
                coordinates.push(global);
                if replaced[local] {
                    replaced_stiffness[local] = finite(mode.angular_frequency_rad_s * mode.angular_frequency_rad_s)?;
                }
            }
            let contiguous = coordinates.windows(2).all(|pair| pair[1] == pair[0] + 1)
                .then(|| coordinates[0]..coordinates[0] + count);
            let state_validation = binding.potential.requires_velocity_validation().then(|| {
                PotentialStateValidation {
                    coordinates: binding.coordinates,
                    q: vec![0.0; count],
                    velocity: vec![0.0; count],
                }
            });
            potentials.push(BoundPotential {
                coordinates,
                contiguous,
                replaced_stiffness,
                state_validation,
                potential: binding.potential,
            });
        }
        let mut state = NonlinearFrictionState {
            potentials,
            histories,
            staged_history_sample: None,
            config,
            last: None,
            preload: None,
        };
        state.validate_states(self.components())?;
        let q = coordinates(self.components());
        let total = finite(self.network.total_energy_j()? + state.potential_energy(&q)? + state.history_energy(&q)?)?;
        limit(
            "linear plus nonlinear energy",
            total,
            self.network.config.maximum_total_energy_j,
        )?;
        self.nonlinear = Some(state);
        poll(Some(gate))?;
        Ok(self)
    }

    /// Cold extension of the original full operator. Existing signed columns,
    /// initial states, material interfaces, nonlinear storage and pressure
    /// weights are preserved. Appending components leaves every original
    /// constitutive coordinate index unchanged; new coordinates participate in
    /// that same nonlinear solve through the complete mechanical network.
    /// No already-playing history can be reset.
    pub(crate) fn with_internal_components(
        self, additional: Vec<ModalAcousticTimeModel>, gate: &CancelGate,
    ) -> Result<Self, ModalFrictionError> {
        self.with_internal_components_and_connections(additional, Vec::new(), &[], gate)
    }

    /// Append cold internal coordinates and their reciprocal laws to the SAME
    /// full operator. Columns use the complete extended coordinate order;
    /// constant-loss indices address only the additional connection list.
    /// Existing physical columns and all original source/history indices stay
    /// unchanged. This is preparation, never a reset of a playing network.
    pub(crate) fn with_internal_components_and_connections(
        self,
        additional: Vec<ModalAcousticTimeModel>,
        connections: Vec<super::super::super::ModalColumnConnection>,
        anchors: &[super::super::super::ConstantLossConnection],
        gate: &CancelGate,
    ) -> Result<Self, ModalFrictionError> {
        poll(Some(gate))?;
        if self.samples_rendered() != 0 || additional.is_empty() {
            return Err(invalid("internal coordinates require cold nonempty preparation").into());
        }
        let count = additional.iter().try_fold(self.mode_count(), |count, model| {
            count.checked_add(model.modes().len()).ok_or_else(|| invalid("internal modal dimension overflow"))
        })?;
        if count > self.network.config.max_modes {
            return Err(invalid("internal coordinates exceed the original complete-network mode budget").into());
        }
        let original_links = self.network.columns.len();
        let link_count = original_links.checked_add(connections.len())
            .ok_or_else(|| invalid("internal connection count overflow"))?;
        if anchors.iter().any(|anchor| anchor.connection >= connections.len()) {
            return Err(invalid("internal constant-loss anchor names an absent new connection").into());
        }
        let mut models = self.network.models.clone();
        models.extend(additional);
        // Admit complete connection work BEFORE allocating the enlarged columns.
        super::super::super::connection_layout(&models, link_count, self.network.config, gate)?;
        let mut columns = Vec::with_capacity(link_count);
        for (column, law) in self.network.columns.iter().zip(&self.network.connections) {
            poll(Some(gate))?;
            let mut column = column.clone();
            column.resize(count, 0.0);
            columns.push(super::super::super::ModalColumnConnection {
                column, stiffness_n_m: law.stiffness_n_m,
                damping_n_s_m: law.damping_n_s_m, rest_extension_m: law.rest_extension_m,
            });
        }
        columns.extend(connections);
        // Original anchors keep their indices; new laws join the same coupled
        // constant-loss basis, without reinterpreting Kelvin--Voigt dashpots.
        let mut constant_loss = self.network.constant_loss_connections();
        constant_loss.extend(anchors.iter().map(|anchor| super::super::super::ConstantLossConnection {
            connection: original_links + anchor.connection,
            reference_omega_rad_s: anchor.reference_omega_rad_s,
        }));
        let network = CoupledModalSystem::new_with_columns_preparation(models, columns,
            &constant_loss, self.network.config, self.network.preparation(), gate)?;
        let network = if self.network.relaxation_branches.is_empty() { network } else {
            let mut branches = self.network.relaxation_branches.clone();
            for branch in &mut branches { branch.projection.resize(2 * count, 0.0); }
            network.with_material_relaxation(branches, Some(self.network.material_memory.clone()), gate)?
        };
        let network=if self.network.plane_stress_material_points.is_empty() {network} else {
            network.with_plane_stress_material_points(self.network.plane_stress_material_points.clone(),gate)?
        };
        let network = if let Some(air) = &self.network.radiation {
            // Cold appended mechanical coordinates leave every original body
            // velocity address intact; the same acoustic history remains at rest.
            network.with_radiation_load(&air.model, &air.source_modes, air.maximum_poles, gate)?
        } else { network };
        let mut next = Self::new(network, self.port, self.config, gate)?;
        if let Some(mut weights) = self.acceleration_pressure {
            weights.resize(count, 0.0);
            next = next.with_acceleration_pressure(weights)?;
        }
        if let Some(interface) = self.interface {
            next = next.with_interface(interface)?;
        }
        // Bound potentials retain their exact original coordinate prefix and
        // any replaced linear stiffness. A larger network must not turn their
        // storage into a diagnostic-only source, or silently lose it entirely.
        next.adhesive = self.adhesive;
        next.nonlinear = self.nonlinear;
        next.validate_material_history_band()?;
        limit("complete initial material/nonlinear energy", next.total_energy_j()?,
            next.network.config.maximum_total_energy_j)?;
        poll(Some(gate))?;
        Ok(next)
    }
    /// Whether an explicit constitutive participant is retained by this owner.
    pub fn has_nonlinear_potentials(&self) -> bool {
        self.nonlinear.is_some()
    }
    /// Last jointly accepted nonlinear sample; refusal leaves it unchanged.
    pub fn nonlinear_frame(&self) -> Option<&NonlinearFrictionFrame> {
        self.nonlinear
            .as_ref()
            .and_then(|state| state.last.as_ref())
    }
    /// Conservative correction [J] at accepted coordinates, without advancing
    /// state. May be signed when a complete strain law replaces linear storage.
    pub fn nonlinear_energy_j(&self) -> Result<f64, ModalCouplingError> {
        match &self.nonlinear {
            Some(state) => {
                let q=coordinates(self.components());
                finite(state.potential_energy(&q)?+state.history_energy(&q)?)
            },
            None => Ok(0.0),
        }
    }
}
impl NonlinearFrictionState {
    pub(super) fn initialize_equilibrium(&mut self, network: &mut CoupledModalSystem,
        external: &[f64], obstacles: Option<&super::multiple::obstacles::ObstacleSet>,
        gate: &CancelGate) -> Result<NonlinearPreloadDiagnostics, ModalFrictionError>
    {
        if !self.histories.is_empty() {
            return Err(invalid("initialize static equilibrium before binding hereditary material history").into());
        }
        equilibrium::initialize(self, network, external, obstacles, gate)
    }
    fn validate_states(&mut self, models: &[ModalAcousticTimeModel]) -> Result<(), ModalCouplingError> {
        for binding in &mut self.potentials {
            let Some(validation) = &mut binding.state_validation else { continue; };
            validation.gather(models)?;
            binding.potential.validate_state(&validation.q, &validation.velocity)?;
        }
        Ok(())
    }
    /// Cold restoration must admit displacement-only laws as well as actual
    /// velocity-dependent limits, while retaining every accepted history.
    pub(super) fn validate_restored_states(&mut self, models: &[ModalAcousticTimeModel],
        linear_memory: &[f64])
        -> Result<(), ModalCouplingError>
    {
        let q = coordinates(models);
        finite(self.potential_energy(&q)? + self.history_energy(&q)?)?;
        self.validate_states(models)?;
        self.validate_history_states(models, linear_memory, false)
    }
    fn potential_energy(&self, q: &[f64]) -> Result<f64, ModalCouplingError> {
        let mut energy = 0.0;
        for binding in &self.potentials {
            let selected = binding.select(q);
            binding.potential.validate_coordinates(&selected)?;
            let value = finite(binding.potential.hamiltonian(&selected))?;
            if value < 0.0 {
                return Err(invalid(
                    "additive nonlinear storage must have a nonnegative energy reference",
                ));
            }
            let reference: f64 = binding.replaced_stiffness.iter().zip(selected.iter())
                .map(|(stiffness, q)| 0.5 * stiffness * q * q).sum();
            energy = finite(energy + finite(value - reference)?)?;
        }
        Ok(energy)
    }
    pub(super) fn stage(
        &mut self,
        inner: &mut FrictionModalSystem,
        external: &[f64],
        drive: FrictionDrive,
        gate: Option<&CancelGate>,
    ) -> Result<NonlinearFrictionFrame, ModalFrictionError> {
        let mut trial = SingleTrial { inner, drive, frame: None };
        let correction = self.stage_trial(&mut trial, external, gate)?;
        let mut frame = trial.frame.ok_or_else(|| invalid("nonlinear solve accepted no mechanical trial"))?;
        frame.observer_pressure_pa = correction.observer_pressure_pa;
        frame.network_energy_j = correction.network_energy_j;
        frame.network_dissipation_j = correction.network_dissipation_j;
        frame.external_work_j = correction.external_work_j;
        frame.energy_residual_j = correction.energy_residual_j;
        frame.energy_tolerance_j = correction.energy_tolerance_j;
        Ok(NonlinearFrictionFrame {
            friction: frame,
            nonlinear_energy_j: correction.diagnostics.nonlinear_energy_j,
            nonlinear_dissipation_j: correction.diagnostics.nonlinear_dissipation_j,
            nonlinear_iterations: correction.diagnostics.nonlinear_iterations,
            reaction_residual: correction.diagnostics.reaction_residual,
            reaction_tolerance: correction.diagnostics.reaction_tolerance,
        })
    }
    // Every trial first chooses each port's root from the ACCEPTED step's
    // branch memory, so an accepted fixed point is the branch rule's own
    // answer at its own forces. A port at the static cone's edge can have no
    // such point: sticking is inadmissible under the reaction of its own stuck
    // response and admissible under that of the sliding one, so the trials
    // alternate until the budget ends. Only then is the step solved again from
    // the same accepted state with each trial's memory set to the previous
    // trial's resolved branches: a port that had to leave its branch keeps the
    // new one while that branch's root exists. A step the joint sweep cannot
    // converge goes to the active-set rescue, which keeps no branch memory and
    // restarts its walk at all-stuck on every trial. So in this pass the
    // rescue first tries the partition the previous trial resolved, kept while
    // it holds a consistent state. A port may change branch once, and partition
    // once, in this pass; a second change, or any failure, returns the first
    // pass's refusal. A step the first pass accepts is never solved again.
    pub(super) fn stage_trial<T: NonlinearTrial>(
        &mut self,
        trial: &mut T,
        external: &[f64],
        gate: Option<&CancelGate>,
    ) -> Result<NonlinearCorrection, ModalFrictionError> {
        let ratio = match self.iterate(trial, external, gate, false)? {
            Ok(accepted) => return Ok(accepted),
            Err(Exhausted { ratio, .. }) if ratio > 1.0 => ratio,
            Err(Exhausted { energy: (residual_j, tolerance_j), .. }) => {
                return Err(ModalCouplingError::EnergyBalance { residual_j, tolerance_j }.into());
            }
        };
        let accepted_memory = trial.branch_memory();
        let held = self.iterate(trial, external, gate, true);
        trial.set_branch_memory(&accepted_memory);
        trial.release_partition();
        if let Ok(Ok(mut accepted)) = held {
            // Disclosed: every trial staged for this step, both passes.
            accepted.diagnostics.nonlinear_iterations = self.config.max_iterations
                .saturating_add(accepted.diagnostics.nonlinear_iterations);
            return Ok(accepted);
        }
        poll(gate)?;
        Err(ModalCouplingError::SolveResidual { relative: ratio, tolerance: 1.0 }.into())
    }
    fn iterate<T: NonlinearTrial>(
        &mut self,
        trial: &mut T,
        external: &[f64],
        gate: Option<&CancelGate>,
        hold: bool,
    ) -> Result<Result<NonlinearCorrection, Exhausted>, ModalFrictionError> {
        self.staged_history_sample = None;
        poll(gate)?;
        let inner = trial.core();
        if external.len() != inner.mode_count() || external.iter().any(|v| !v.is_finite()) {
            return Err(
                invalid("external force must retain every nonlinear bank coordinate").into(),
            );
        }
        self.validate_states(inner.components())?;
        self.validate_history_states(inner.components(), inner.material_memory(), false)?;
        let before_q = coordinates(inner.components());
        let dt = inner.sample_period_s();
        let has_loss = self.potentials.iter().any(|p| p.potential.has_constitutive_loss());
        let before_potential = finite(self.potential_energy(&before_q)?+self.history_energy(&before_q)?)?;
        let before = finite(inner.network.total_energy_j()? + before_potential)?;
        let mut reaction = vec![0.0; inner.mode_count()];
        self.reactions(&before_q, &before_q, &mut reaction, gate)?;
        let mut loss_reaction = if has_loss { vec![0.0; inner.mode_count()] } else { Vec::new() };
        if has_loss {
            self.loss_reactions(&before_q, &velocities(inner.components()), &mut loss_reaction, gate)?;
            self.add_loss_reaction(&mut reaction, &loss_reaction)?;
        }
        let mut history_reaction=if self.histories.is_empty() {Vec::new()} else {vec![0.0;inner.mode_count()]};
        if !self.histories.is_empty() {
            self.history_reactions(&before_q,false,&mut history_reaction,gate)?;
            self.add_loss_reaction(&mut reaction,&history_reaction)?;
        }
        let mut next_reaction = reaction.clone();
        let mut forces = external.to_vec();
        let mut last_ratio = f64::INFINITY;
        let mut last_energy = (0.0, 0.0);
        let mut held: Option<Vec<ContactBranch>> = None;
        let mut changes: Vec<u8> = Vec::new();
        let mut partition_changes: Vec<u8> = Vec::new();
        for iteration in 1..=self.config.max_iterations {
            poll(gate)?;
            for i in 0..forces.len() {
                forces[i] = if reaction[i] == 0.0 {
                    external[i]
                } else {
                    finite(external[i] + reaction[i])?
                };
            }
            if let Some(previous) = &held { trial.set_branch_memory(previous); }
            let frame = trial.stage(&forces, gate)?;
            if hold {
                let staged = trial.staged_branches();
                if let Some(previous) = &held {
                    for ((count, now), before) in changes.iter_mut().zip(&staged).zip(previous) {
                        if now != before { *count = count.saturating_add(1); }
                    }
                    if changes.iter().any(|&count| count > 1) { break; }
                } else {
                    changes = vec![0; staged.len()];
                }
                held = Some(staged);
                let moved = trial.hold_partition();
                partition_changes.resize(moved.len(), 0);
                for (count, moved) in partition_changes.iter_mut().zip(moved) {
                    if moved { *count = count.saturating_add(1); }
                }
                if partition_changes.iter().any(|&count| count > 1) { break; }
            }
            let inner = trial.core();
            let after_q = coordinates(&inner.network.candidates);
            let mut after_potential = self.potential_energy(&after_q)?;
            self.reactions(&before_q, &after_q, &mut next_reaction, gate)?;
            let mut nonlinear_dissipation = if has_loss {
                let midpoint: Vec<_> = before_q.iter().zip(&after_q)
                    .map(|(&a, &b)| f64::midpoint(a, b)).collect();
                let increment: Vec<_> = before_q.iter().zip(&after_q).map(|(&a, &b)| b - a).collect();
                let velocity: Vec<_> = increment.iter().map(|&delta| delta / dt).collect();
                self.loss_reactions(&midpoint, &velocity, &mut loss_reaction, gate)?;
                self.add_loss_reaction(&mut next_reaction, &loss_reaction)?;
                // Work conjugacy uses actual held-force displacement, never a
                // separately integrated endpoint velocity or guessed decay ratio.
                finite(-dot(&loss_reaction, &increment)?)?
            } else { 0.0 };
            if !self.histories.is_empty() {
                let (stored,loss)=self.stage_histories(&before_q,&after_q,dt,&mut history_reaction,
                    inner.network.config.energy_absolute_tolerance_j,
                    inner.network.config.energy_relative_tolerance,gate)?;
                self.add_loss_reaction(&mut next_reaction,&history_reaction)?;
                after_potential=finite(after_potential+stored)?;
                nonlinear_dissipation=finite(nonlinear_dissipation+loss)?;
            }
            if nonlinear_dissipation < 0.0 {
                return Err(invalid("nonlinear constitutive reaction is not passive").into());
            }
            let network_dissipation = finite(frame.network_dissipation_j + nonlinear_dissipation)?;
            let residual = reaction
                .iter()
                .zip(&next_reaction)
                .fold(0.0_f64, |m, (a, b)| m.max((a - b).abs()));
            let scale = reaction
                .iter()
                .chain(&next_reaction)
                .fold(0.0_f64, |m, v| m.max(v.abs()));
            let tolerance = finite(
                self.config.force_absolute_tolerance + self.config.force_relative_tolerance * scale,
            )?;
            last_ratio = finite(residual / tolerance)?;
            let after = finite(frame.network_energy_j + after_potential)?;
            let external_work = dot(external, &inner.network.free_delta)?;
            // Layer and series-drive springs stay outside the modal energy
            // report. Their original accepted and staged storage still enter
            // this complete balance, without duplicating plastic/dashpot loss.
            let (stored_before, stored_after) = match frame.internal_storage {
                Some(storage) => (finite(before + storage.before_energy_j)?,
                    finite(after + storage.after_energy_j)?),
                None => (before, after),
            };
            // Keep obstacle storage separate in the public frame, but include
            // BOTH accepted endpoints and its own loss/work in this transaction.
            let (total_before, total_after) = match frame.contacts {
                Some(contact) => (finite(stored_before + contact.before_energy_j)?,
                    finite(stored_after + contact.after_energy_j)?),
                None => (stored_before, stored_after),
            };
            let balance = if let Some(contact) = frame.contacts {
                finite(total_after - total_before + network_dissipation
                    + frame.friction_dissipation_j + frame.declared_slip_work_j + contact.dissipation_j
                    - external_work - frame.drive_work_j - contact.work_j)?
            } else { finite(
                total_after - total_before + network_dissipation + frame.friction_dissipation_j
                    + frame.declared_slip_work_j
                    - external_work
                    - frame.drive_work_j,
            )? };
            let mut energy_scale = total_before
                .max(total_after)
                .max(network_dissipation)
                .max(frame.friction_dissipation_j)
                .max(frame.declared_slip_work_j.abs())
                .max(external_work.abs())
                .max(frame.absolute_drive_work_j);
            if let Some(contact) = frame.contacts {
                energy_scale = energy_scale.max(contact.dissipation_j).max(contact.absolute_work_j);
            }
            let energy_tolerance = finite(
                inner.network.config.energy_absolute_tolerance_j
                    + inner.network.config.energy_relative_tolerance * energy_scale,
            )?;
            last_energy = (balance, energy_tolerance);
            if residual <= tolerance && balance.abs() <= energy_tolerance {
                // The constitutive limit sees the actual converged endpoint,
                // independently of acceleration observation and loss quadrature.
                self.validate_states(&inner.network.candidates)?;
                self.validate_history_states(&inner.network.candidates,
                    &inner.network.candidate_material_memory, true)?;
                limit(
                    "total linear, nonlinear and contact energy",
                    total_after,
                    inner.network.config.maximum_total_energy_j,
                )?;
                // Actual endpoint pressure uses constitutive force at (q1,v1),
                // not the held discrete-gradient/midpoint reaction.
                let mut observer_pressure_pa = frame.observer_pressure_pa;
                if let Some(weights) = inner.acceleration_pressure_weights() {
                    self.reactions(&after_q, &after_q, &mut next_reaction, gate)?;
                    if has_loss {
                        self.loss_reactions(&after_q, &velocities(&inner.network.candidates),
                            &mut loss_reaction, gate)?;
                        self.add_loss_reaction(&mut next_reaction, &loss_reaction)?;
                    }
                    if !self.histories.is_empty() {
                        self.history_reactions(&after_q,true,&mut history_reaction,gate)?;
                        self.add_loss_reaction(&mut next_reaction,&history_reaction)?;
                    }
                    for ((weight, endpoint), held) in
                        weights.iter().zip(&next_reaction).zip(&reaction)
                    {
                        observer_pressure_pa =
                            finite(observer_pressure_pa + weight * (endpoint - held))?;
                    }
                    limit(
                        "nonlinear acceleration observer pressure",
                        observer_pressure_pa.abs(),
                        inner.network.config.maximum_abs_pressure_pa,
                    )?;
                }
                let accepted = NonlinearCorrection {
                    observer_pressure_pa,
                    network_energy_j: after,
                    network_dissipation_j: network_dissipation,
                    external_work_j: external_work,
                    energy_residual_j: balance,
                    energy_tolerance_j: energy_tolerance,
                    diagnostics: NonlinearFrictionDiagnostics {
                        nonlinear_energy_j: after_potential,
                        nonlinear_dissipation_j: nonlinear_dissipation,
                        nonlinear_iterations: iteration,
                        reaction_residual: residual,
                        reaction_tolerance: tolerance,
                    },
                };
                poll(gate)?;
                self.staged_history_sample = Some(inner.samples_rendered());
                return Ok(Ok(accepted));
            }
            std::mem::swap(&mut reaction, &mut next_reaction);
        }
        Ok(Err(Exhausted { ratio: last_ratio, energy: last_energy }))
    }
    pub(super) fn publish(&mut self, frame: NonlinearFrictionFrame) {
        self.publish_histories();
        self.last = Some(frame);
    }
    fn add_loss_reaction(&self, reaction: &mut [f64], loss: &[f64]) -> Result<(), ModalFrictionError> {
        for (reaction, &loss) in reaction.iter_mut().zip(loss) {
            // Preserve the original elastic arithmetic when no loss is selected.
            if loss != 0.0 { *reaction = finite(*reaction + loss)?; }
            limit("complete nonlinear generalized reaction", reaction.abs(),
                self.config.maximum_generalized_force)?;
        }
        Ok(())
    }
    fn loss_reactions(
        &self, q: &[f64], velocity: &[f64], out: &mut [f64], gate: Option<&CancelGate>,
    ) -> Result<(), ModalFrictionError> {
        out.fill(0.0);
        if q.len() != velocity.len() || q.len() != out.len()
            || q.iter().chain(velocity).any(|x| !x.is_finite()) {
            return Err(invalid("nonlinear loss requires complete finite coordinates and velocities").into());
        }
        for binding in &self.potentials {
            poll(gate)?;
            if !binding.potential.has_constitutive_loss() { continue; }
            let selected_q = binding.select(q);
            let selected_v = binding.select(velocity);
            let mut reaction = vec![0.0; binding.coordinates.len()];
            binding.potential.dissipative_reaction(&selected_q, &selected_v, &mut reaction)?;
            if finite(dot(&reaction, &selected_v)?)? > 0.0 {
                return Err(invalid("nonlinear constitutive reaction supplies energy").into());
            }
            for (&index, value) in binding.coordinates.iter().zip(reaction) {
                out[index] = finite(value)?;
            }
        }
        Ok(())
    }
    fn reactions(
        &self,
        before: &[f64],
        after: &[f64],
        out: &mut [f64],
        gate: Option<&CancelGate>,
    ) -> Result<(), ModalFrictionError> {
        out.fill(0.0);
        for binding in &self.potentials {
            poll(gate)?;
            let before = binding.select(before);
            let after = binding.select(after);
            let gradient = fs_phs::discrete_gradient(
                binding.potential.as_ref(),
                &before,
                &after,
            );
            for (local, (&index, value)) in binding.coordinates.iter().zip(gradient).enumerate() {
                let reference = binding.replaced_stiffness[local] * f64::midpoint(before[local], after[local]);
                out[index] = finite(if binding.replaced_stiffness[local] == 0.0 { -value } else { reference - value })?;
                limit(
                    "nonlinear generalized reaction",
                    out[index].abs(),
                    self.config.maximum_generalized_force,
                )?;
            }
        }
        Ok(())
    }
}
fn coordinates(models: &[ModalAcousticTimeModel]) -> Vec<f64> {
    models
        .iter()
        .flat_map(|model| {
            model
                .states()
                .iter()
                .map(|state| state.displacement_m_sqrt_kg)
        })
        .collect()
}
fn velocities(models: &[ModalAcousticTimeModel]) -> Vec<f64> {
    models.iter().flat_map(|model| model.states().iter()
        .map(|state| state.velocity_m_sqrt_kg_per_s)).collect()
}
