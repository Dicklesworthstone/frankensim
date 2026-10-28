//! Nonlinear storage coupled to the SAME exact linear and friction owners.
//!
//! Each trial supplies a held reaction -DG(V,q0,q1) to the original complete
//! linear network. The original scalar or simultaneous bow solve includes that
//! reaction. The fixed point is accepted only when the reaction law AND the combined energy
//! balance close. fs-phs owns DG; this adapter changes neither the exponential
//! linear propagator nor the friction law, its stick/slip branch or work port.
//! This is a bounded iterated discrete-gradient composition, not SAV, a unique
//! nonlinear-solution claim, a trajectory certificate or a real-time guarantee.
use super::*;
use fs_phs::Storage;
use std::borrow::Cow;

/// Displacement-only positive constitutive energy. Contiguous bindings supply
/// additional energy; mapped bindings may explicitly replace selected original
/// diagonal stiffness. Inertia remains with the original linear owner. The
/// constitutive owner admits its own physical validity envelope.
pub trait ModalPotential: Storage + Send + Sync {
    /// Number of original mass-normalized displacement coordinates.
    fn coordinate_count(&self) -> usize;
    /// Reject nonfinite coordinates or motion outside the constitutive domain.
    fn validate_coordinates(&self, q: &[f64]) -> Result<(), ModalCouplingError>;
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
    /// Maximum complete linear/friction trials, in 1..=128.
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
}
pub(super) struct NonlinearTrialWork {
    pub observer_pressure_pa: f64,
    pub network_energy_j: f64,
    pub drive_work_j: f64,
    pub absolute_drive_work_j: f64,
    pub network_dissipation_j: f64,
    pub friction_dissipation_j: f64,
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
            internal_storage: None,
            contacts: None,
        })
    }
}
struct BoundPotential {
    coordinates: Vec<usize>,
    contiguous: Option<std::ops::Range<usize>>,
    replaced_stiffness: Vec<f64>,
    potential: Box<dyn ModalPotential>,
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
    config: NonlinearModalConfig,
    last: Option<NonlinearFrictionFrame>,
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
        let mut potentials = self.nonlinear.take().map_or_else(Vec::new, |state| state.potentials);
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
            potentials.push(BoundPotential {
                coordinates,
                contiguous,
                replaced_stiffness,
                potential: binding.potential,
            });
        }
        let state = NonlinearFrictionState {
            potentials,
            config,
            last: None,
        };
        let q = coordinates(self.components());
        let total = finite(self.network.total_energy_j()? + state.potential_energy(&q)?)?;
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
        let network = match self.network.prepared_budget() {
            Some(budget) => CoupledModalSystem::new_with_columns_inner(models, columns, self.network.config, Some(budget),
                &constant_loss, gate),
            None if constant_loss.is_empty() => CoupledModalSystem::new_with_columns(models, columns, self.network.config, gate),
            None => return Err(invalid("constant-loss internal connections require prepared propagation").into()),
        }?;
        let network = if self.network.relaxation_branches.is_empty() { network } else {
            let mut branches = self.network.relaxation_branches.clone();
            for branch in &mut branches { branch.projection.resize(2 * count, 0.0); }
            network.with_material_relaxation(branches, Some(self.network.material_memory.clone()), gate)?
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
        next.nonlinear = self.nonlinear;
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
            Some(state) => state.potential_energy(&coordinates(self.components())),
            None => Ok(0.0),
        }
    }
}
impl NonlinearFrictionState {
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
        frame.external_work_j = correction.external_work_j;
        frame.energy_residual_j = correction.energy_residual_j;
        frame.energy_tolerance_j = correction.energy_tolerance_j;
        Ok(NonlinearFrictionFrame {
            friction: frame,
            nonlinear_energy_j: correction.diagnostics.nonlinear_energy_j,
            nonlinear_iterations: correction.diagnostics.nonlinear_iterations,
            reaction_residual: correction.diagnostics.reaction_residual,
            reaction_tolerance: correction.diagnostics.reaction_tolerance,
        })
    }
    pub(super) fn stage_trial<T: NonlinearTrial>(
        &self,
        trial: &mut T,
        external: &[f64],
        gate: Option<&CancelGate>,
    ) -> Result<NonlinearCorrection, ModalFrictionError> {
        poll(gate)?;
        let inner = trial.core();
        if external.len() != inner.mode_count() || external.iter().any(|v| !v.is_finite()) {
            return Err(
                invalid("external force must retain every nonlinear bank coordinate").into(),
            );
        }
        let before_q = coordinates(inner.components());
        let before_potential = self.potential_energy(&before_q)?;
        let before = finite(inner.network.total_energy_j()? + before_potential)?;
        let mut reaction = vec![0.0; inner.mode_count()];
        self.reactions(&before_q, &before_q, &mut reaction, gate)?;
        let mut next_reaction = reaction.clone();
        let mut forces = external.to_vec();
        let mut last_ratio = f64::INFINITY;
        let mut last_energy = (0.0, 0.0);
        for iteration in 1..=self.config.max_iterations {
            poll(gate)?;
            for i in 0..forces.len() {
                forces[i] = if reaction[i] == 0.0 {
                    external[i]
                } else {
                    finite(external[i] + reaction[i])?
                };
            }
            let frame = trial.stage(&forces, gate)?;
            let inner = trial.core();
            let after_q = coordinates(&inner.network.candidates);
            let after_potential = self.potential_energy(&after_q)?;
            self.reactions(&before_q, &after_q, &mut next_reaction, gate)?;
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
                finite(total_after - total_before + frame.network_dissipation_j
                    + frame.friction_dissipation_j + contact.dissipation_j
                    - external_work - frame.drive_work_j - contact.work_j)?
            } else { finite(
                total_after - total_before + frame.network_dissipation_j + frame.friction_dissipation_j
                    - external_work
                    - frame.drive_work_j,
            )? };
            let mut energy_scale = total_before
                .max(total_after)
                .max(frame.network_dissipation_j)
                .max(frame.friction_dissipation_j)
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
                limit(
                    "total linear, nonlinear and contact energy",
                    total_after,
                    inner.network.config.maximum_total_energy_j,
                )?;
                // Actual endpoint pressure uses -grad V(q1), not held DG.
                let mut observer_pressure_pa = frame.observer_pressure_pa;
                if let Some(weights) = inner.acceleration_pressure_weights() {
                    self.reactions(&after_q, &after_q, &mut next_reaction, gate)?;
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
                    external_work_j: external_work,
                    energy_residual_j: balance,
                    energy_tolerance_j: energy_tolerance,
                    diagnostics: NonlinearFrictionDiagnostics {
                        nonlinear_energy_j: after_potential,
                        nonlinear_iterations: iteration,
                        reaction_residual: residual,
                        reaction_tolerance: tolerance,
                    },
                };
                poll(gate)?;
                return Ok(accepted);
            }
            std::mem::swap(&mut reaction, &mut next_reaction);
        }
        if last_ratio <= 1.0 {
            return Err(ModalCouplingError::EnergyBalance {
                residual_j: last_energy.0,
                tolerance_j: last_energy.1,
            }
            .into());
        }
        Err(ModalCouplingError::SolveResidual {
            relative: last_ratio,
            tolerance: 1.0,
        }
        .into())
    }
    pub(super) fn publish(&mut self, frame: NonlinearFrictionFrame) {
        self.last = Some(frame);
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
