use super::super::{
    BodyForce3d, Error, ForceChange3d, ForceControl3d, ForcingWorkspace3d, ReservoirFaceState,
    StepReport, forcing, reservoir_stages::Context,
};
use super::*;
use crate::compressible::mixture::molecular::stabilized::{IntervalFailure, StepWork};
pub use crate::compressible::mixture::molecular::{
    MolecularStepControl, StabilizedMolecularStepControl,
};

mod moving;
pub use moving::MovingMolecularStepReport3d;

#[cfg(feature = "experimental-sgs")]
mod lagrangian_channel;
#[cfg(feature = "experimental-sgs")]
mod lagrangian_species;

/// Complete accepted volumetric multicomponent molecular/advection interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MolecularStepReport3d {
    /// Common physical duration; each molecular half represents half.
    pub elapsed: Time,
    /// Complete rejected trials before the accepted interval.
    pub coupling_retries: usize,
    /// First molecular half interval.
    pub first: MolecularReport,
    /// Central original spatial advection, including supplied body forcing.
    pub advection: StepReport,
    /// Second molecular half interval.
    pub second: MolecularReport,
    /// Original extensive field.
    pub before: Totals,
    /// Final extensive field.
    pub after: Totals,
    /// Complete inward molecular and advective boundary exchange.
    pub boundary_change: Totals,
    /// Body impulse and signed work, separate from boundary exchange.
    pub source_change: ForceChange3d,
    /// `after - before - boundary_change - source_change`.
    pub conservation_residual: Totals,
    /// Physical chemical-plus-kinetic energy residual from conservative differences.
    pub absolute_energy_residual_j: f64,
}

/// Complete spatial mixture split and actual RKL2 integration work.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StabilizedMolecularStepReport3d {
    /// Original physical balances and reports for the accepted complete split.
    pub transport: MolecularStepReport3d,
    /// Both molecular halves and all discarded complete-interval probes.
    pub operator_evaluations: usize,
    /// Largest stage count in the final accepted pair of molecular halves.
    pub maximum_stages_used: usize,
    /// Largest scaled parabolic number in the final accepted molecular halves.
    pub maximum_scaled_parabolic_number: f64,
}

/// Atomic nonreacting multicomponent viscous, thermal and advective flow.
///
/// The prepared molecular operator retains its original source, physical walls
/// and optional SGS configuration. Its requested-duration halves surround the
/// original three-dimensional mixture Euler transaction. No chemistry, frozen
/// composition surrogate or additional viscous-heating source is introduced.
///
/// Shortened or refused phases restart the complete interval from accepted
/// cells. Final cancellation preserves both owners' published state and loads.
/// Symmetric composition permits second-order smooth evolution when the
/// constituent paths retain their order; limiting and first-order retries can
/// reduce it. No entropy, turbulence or experimental-validation claim is added.
pub struct MolecularEuler3d<L: TransportLaw> {
    accepted: MixtureEuler3d,
    trial: MixtureEuler3d,
    molecular: MolecularTransport3d<L>,
    boundary_impulses: Vec<f64>,
    boundary_trial: Vec<f64>,
    maximum_workspace_bytes: usize,
    last_step: Option<MolecularStepReport3d>,
}
impl<L: TransportLaw> MolecularEuler3d<L> {
    /// Both full fields, prepared molecular/SGS buffers and complete boundary
    /// buffers, including incoming excess capacity. Immutable sources, bounded
    /// stack scratch and provider storage are excluded.
    pub fn required_workspace_bytes(
        flow: &MixtureEuler3d,
        molecular: &MolecularTransport3d<L>,
    ) -> Result<usize, DissipationError> {
        flow.workspace_bytes()
            .checked_add(MixtureEuler3d::required_workspace_bytes(
                flow.gas(),
                flow.grid(),
            )?)
            .and_then(|bytes| bytes.checked_add(flow.reservoirs.workspace_bytes()))
            .and_then(|bytes| bytes.checked_add(molecular.workspace_bytes()))
            .and_then(|bytes| {
                bytes.checked_add(
                    flow.boundary_impulses
                        .len()
                        .checked_mul(2 * size_of::<f64>())?,
                )
            })
            .ok_or(DissipationError::Invalid(
                "coupled spatial mixture workspace size overflow",
            ))
    }

    /// Complete coupled storage plus unallocated RKL2 face/boundary scratch.
    /// Add supplied forcing workspace to this allowance before calling `new`.
    pub fn required_stabilized_workspace_bytes(
        flow: &MixtureEuler3d,
        molecular: &MolecularTransport3d<L>,
    ) -> Result<usize, DissipationError> {
        Self::required_workspace_bytes(flow, molecular)?
            .checked_add(molecular.additional_stabilized_workspace_bytes()?)
            .ok_or(DissipationError::Invalid(
                "coupled stabilized spatial mixture workspace size overflow",
            ))
    }

    /// Additional numerical payload for Lagrangian SGS: 216 bytes per cell.
    /// Add this to `required_workspace_bytes` or
    /// `required_stabilized_workspace_bytes` for the molecular pair, plus any
    /// supplied forcing workspace. Actual capacities are admitted after allocation.
    #[cfg(feature = "experimental-sgs")]
    pub fn required_lagrangian_dynamic_sgs_workspace_bytes(
        flow: &MixtureEuler3d,
    ) -> Result<usize, DissipationError> {
        sgs::SubgridTransport::required_lagrangian_bytes(flow)
    }

    /// Bind periodic Lagrangian stress to the original complete nonreacting split.
    ///
    /// The supplied molecular operator must not already select an SGS closure.
    /// Its Maxwell--Stefan/Soret/Dufour, original heat capacities and referenced
    /// partial enthalpies remain the source of molecular and modeled transport.
    /// Fixed turbulent Prandtl/Schmidt numbers do not imply scalar-memory fits.
    ///
    /// One fit samples the accepted physical START state. Both molecular halves,
    /// every SSP/RKL stage and all whole-split retries use its held coefficient;
    /// face density, strain and composition still follow each actual stage.
    /// History/report publication occurs only with complete flow acceptance.
    /// The pathline/coefficient update is first order in time. This constructor
    /// starts a fresh history; no chemistry or reacting-history coupling is added.
    ///
    /// The byte cap covers the complete original pair and all additional SGS
    /// capacities. Detaching the history during a step does not remove its
    /// storage from that cap. No history/filter allocation occurs during stepping.
    #[cfg(feature = "experimental-sgs")]
    pub fn with_lagrangian_dynamic_sgs(
        flow: MixtureEuler3d,
        molecular: MolecularTransport3d<L>,
        config: sgs::MixtureLagrangianSubgridConfig,
        maximum_workspace_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, DissipationError> {
        poll(&mut cancelled)?;
        config.validate()?;
        if molecular.sgs.is_some() {
            return Err(DissipationError::Invalid(
                "Lagrangian mixture attachment requires a molecular operator without SGS",
            ));
        }
        let required = Self::required_workspace_bytes(&flow, &molecular)?
            .checked_add(Self::required_lagrangian_dynamic_sgs_workspace_bytes(
                &flow,
            )?)
            .ok_or(DissipationError::Invalid(
                "coupled Lagrangian mixture workspace size overflow",
            ))?;
        admit_workspace(required, maximum_workspace_bytes)?;
        let sgs = sgs::SubgridTransport::new_lagrangian(&flow, config)?;
        poll(&mut cancelled)?;
        let mut result = Self::new(flow, molecular, maximum_workspace_bytes, &mut cancelled)?;
        result.molecular.sgs = Some(sgs);
        admit_workspace(result.workspace_bytes(), maximum_workspace_bytes)?;
        poll(&mut cancelled)?;
        Ok(result)
    }

    /// Explicit controls for the attached pathline closure; None for other models.
    #[cfg(feature = "experimental-sgs")]
    pub fn lagrangian_dynamic_subgrid_config(&self) -> Option<sgs::MixtureLagrangianSubgridConfig> {
        self.molecular
            .sgs
            .as_ref()
            .and_then(sgs::SubgridTransport::lagrangian_config)
    }

    /// Last committed START-state sample time, distinct from the flow endpoint.
    /// None until a complete interval accepts or when no history model is selected.
    #[cfg(feature = "experimental-sgs")]
    pub fn lagrangian_subgrid_sample_time(&self) -> Option<Time> {
        self.molecular
            .sgs
            .as_ref()
            .and_then(sgs::SubgridTransport::lagrangian_sample_time)
    }

    /// Last committed start-state fit. Its cell diagnostics do not replace
    /// the conservative face's stage-local strain and material admission.
    #[cfg(feature = "experimental-sgs")]
    pub fn lagrangian_subgrid_report(&self) -> Option<sgs::DynamicSgsReport> {
        self.molecular
            .sgs
            .as_ref()
            .and_then(sgs::SubgridTransport::lagrangian_report)
    }

    /// Bind the original field and its already-prepared molecular operator.
    /// The stored complete budget also admits externally supplied forcing scratch,
    /// even for an exactly zero source. Both owners remain immutable to callers.
    pub fn new(
        flow: MixtureEuler3d,
        molecular: MolecularTransport3d<L>,
        maximum_workspace_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, DissipationError> {
        poll(&mut cancelled)?;
        if molecular.gas() != flow.gas()
            || molecular.grid != flow.grid
            || molecular.boundaries != flow.boundaries
        {
            return Err(DissipationError::Invalid(
                "coupled spatial mixture molecular binding differs",
            ));
        }
        reservoir::admit(
            flow.boundaries,
            &flow.reservoirs,
            molecular.pressure_outlets,
        )?;
        let required = Self::required_workspace_bytes(&flow, &molecular)?;
        admit_workspace(required, maximum_workspace_bytes)?;
        let mut states = allocated(flow.states.len(), 0.0)?;
        copy_values(&flow.states, &mut states, &mut cancelled)?;
        let mut trial = MixtureEuler3d::new(
            flow.gas.clone(),
            flow.grid,
            states,
            flow.boundaries,
            MixtureEuler3d::required_workspace_bytes(flow.gas(), flow.grid())?,
            &mut cancelled,
        )?;
        trial.reservoirs = flow.reservoirs.try_clone(&mut cancelled)?;
        trial.workspace_bytes = trial
            .workspace_bytes
            .checked_add(trial.reservoirs.workspace_bytes())
            .ok_or(DissipationError::Invalid(
                "coupled reservoir workspace size overflow",
            ))?;
        copy_physical(&flow, &mut trial, &mut cancelled)?;
        let mut boundary_impulses = allocated(flow.boundary_impulses.len(), 0.0)?;
        copy_values(
            &flow.boundary_impulses,
            &mut boundary_impulses,
            &mut cancelled,
        )?;
        let boundary_trial = allocated(flow.boundary_impulses.len(), 0.0)?;
        let result = Self {
            accepted: flow,
            trial,
            molecular,
            boundary_impulses,
            boundary_trial,
            maximum_workspace_bytes,
            last_step: None,
        };
        admit_workspace(result.workspace_bytes(), maximum_workspace_bytes)?;
        poll(&mut cancelled)?;
        Ok(result)
    }

    /// Complete accepted field and clock. Its Euler metadata describes the
    /// middle phase; this owner's report and face loads describe the whole split.
    pub const fn flow(&self) -> &MixtureEuler3d {
        &self.accepted
    }

    /// Borrow actual accepted-endpoint primitives and viscous face tractions.
    /// The source's Maxwell--Stefan/Fourier/full-stress operator and explicitly
    /// configured SGS/walls are evaluated once under the supplied numerical
    /// controls. Accepted state, clock, closure history and step loads stay put.
    /// Failure returns no view; reusable constitutive caches may change.
    pub fn stress_field(
        &mut self,
        controls: DissipationControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<crate::compressible::spatial::molecular::StressField3d<'_>, DissipationError> {
        poll(&mut cancelled)?;
        #[cfg(feature = "experimental-sgs")]
        if let Some(mut filter) = self
            .molecular
            .sgs
            .as_mut()
            .and_then(sgs::SubgridTransport::take_lagrangian_filter)
        {
            let result = (|| {
                let trial = self
                    .molecular
                    .sgs
                    .as_mut()
                    .expect("selected Lagrangian mixture SGS owner")
                    .prepare_lagrangian_with_walls(
                        &mut filter,
                        &self.accepted,
                        &self.molecular.walls,
                        &mut cancelled,
                    )?;
                let result =
                    self.molecular
                        .prepare_stress_field(&self.accepted, controls, &mut cancelled);
                drop(trial);
                result
            })();
            self.molecular
                .sgs
                .as_mut()
                .expect("retained Lagrangian mixture SGS owner")
                .restore_lagrangian_filter(filter, None);
            result?;
        } else {
            self.molecular
                .prepare_stress_field(&self.accepted, controls, &mut cancelled)?;
        }
        #[cfg(not(feature = "experimental-sgs"))]
        self.molecular
            .prepare_stress_field(&self.accepted, controls, &mut cancelled)?;
        Ok(
            crate::compressible::spatial::molecular::StressField3d::mixture(
                &self.molecular.grid,
                self.molecular.boundaries,
                self.accepted.time.value(),
                &self.molecular.work.primitive,
                &self.molecular.work.flux,
                &self.accepted.states,
                self.accepted.gas.species_count(),
            ),
        )
    }
    /// Atomically refresh all prepared exterior profiles on both split fields.
    /// Prepare their topology/storage on the Euler field before composition.
    /// The callback and original caloric admission run once per face; the
    /// molecular halves and middle HLLC step then use the same held traces.
    /// Cancellation preserves both profiles, cells, clocks and accepted reports.
    pub fn update_reservoir_profiles(
        &mut self,
        sample: impl FnMut(
            BoundarySide,
            usize,
            [f64; 3],
            &mut [f64],
        ) -> Result<super::super::ReservoirFaceState, Error>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<(), Error> {
        super::super::reservoir_profiles::update_pair(
            &mut self.accepted,
            &mut self.trial,
            sample,
            cancelled,
        )
    }
    /// Bound original operator; its published report describes the second half.
    pub const fn molecular(&self) -> &MolecularTransport3d<L> {
        &self.molecular
    }
    /// Actual owned numerical capacities, following state-buffer exchanges.
    pub fn workspace_bytes(&self) -> usize {
        self.accepted.workspace_bytes()
            + self.trial.workspace_bytes()
            + self.molecular.workspace_bytes()
            + (self.boundary_impulses.capacity() + self.boundary_trial.capacity())
                * size_of::<f64>()
    }
    /// Last completely accepted physical interval.
    pub const fn last_step(&self) -> Option<&MolecularStepReport3d> {
        self.last_step.as_ref()
    }
    /// Complete outward exchange through each physical face.
    pub fn boundary_impulses(
        &self,
        side: BoundarySide,
    ) -> impl ExactSizeIterator<Item = BoundaryImpulse<'_>> {
        boundary_view(
            &self.boundary_impulses,
            self.accepted.grid(),
            self.accepted.gas().species_count(),
            side,
        )
    }
    /// Recover both accepted owners and their original wall/source configuration.
    pub fn into_parts(self) -> (MixtureEuler3d, MolecularTransport3d<L>) {
        (self.accepted, self.molecular)
    }

    /// Advance at most the requested duration, restarting both molecular halves
    /// when the central interval shortens. Exhausted whole-step retries return
    /// a `DissipationError::StepFailure` with zero committed elapsed time.
    pub fn step(
        &mut self,
        maximum_duration: Time,
        controls: MolecularStepControl,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MolecularStepReport3d, DissipationError> {
        self.step_impl(maximum_duration, controls, None, None, &[], None, cancelled)
    }

    /// Advance the complete split with prescribed reservoir histories.
    ///
    /// Molecular halves sample physical intervals `[t,t+dt/2]` and
    /// `[t+dt/2,t+dt]`; the central Euler stages sample `[t,t+dt]`. Every
    /// numerical retry replays the same immutable provider at its actual stage
    /// times. One face-evaluation budget covers initial admission, both halves,
    /// advection and all discarded outer attempts. Provider/source refusals and
    /// budget exhaustion are terminal. The callback must fill every canonical
    /// species and be pure and replayable; mutable access is observation scratch.
    ///
    /// Existing prepared profile banks and original molecular physics are
    /// reused. Cells, clock, endpoint profiles, reports and all boundary loads
    /// publish only with a completely accepted split. Smooth prescribed data
    /// retain the constituent second-order timing; discontinuities require
    /// caller-chosen physical steps. Finite wall inertia is explicitly refused.
    pub fn step_with_reservoir_profiles(
        &mut self,
        maximum_duration: Time,
        controls: MolecularStepControl,
        maximum_face_evaluations: usize,
        mut sample: impl FnMut(
            f64,
            BoundarySide,
            usize,
            [f64; 3],
            &mut [f64],
        ) -> Result<ReservoirFaceState, Error>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MolecularStepReport3d, DissipationError> {
        let mut context = Context::new(&mut sample, maximum_face_evaluations);
        self.step_impl(
            maximum_duration,
            controls,
            None,
            None,
            &[],
            Some(&mut context),
            cancelled,
        )
    }

    /// The same physical split with native RKL2 molecular stage abscissae.
    ///
    /// The provider budget spans both molecular halves and the full Euler
    /// stage sequence; the RKL operator budget independently counts every
    /// molecular residual, including endpoint admission and discarded probes.
    pub fn step_stabilized_with_reservoir_profiles(
        &mut self,
        maximum_duration: Time,
        controls: StabilizedMolecularStepControl,
        maximum_face_evaluations: usize,
        mut sample: impl FnMut(
            f64,
            BoundarySide,
            usize,
            [f64; 3],
            &mut [f64],
        ) -> Result<ReservoirFaceState, Error>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<StabilizedMolecularStepReport3d, DissipationError> {
        let mut context = Context::new(&mut sample, maximum_face_evaluations);
        let mut work = StepWork::new(controls.molecular);
        let transport = self.step_impl(
            maximum_duration,
            controls.split(),
            None,
            Some(&mut work),
            &[],
            Some(&mut context),
            cancelled,
        )?;
        Ok(StabilizedMolecularStepReport3d {
            transport,
            operator_evaluations: work.evaluations,
            maximum_stages_used: work.maximum_stages_used,
            maximum_scaled_parabolic_number: work.maximum_scaled_parabolic_number,
        })
    }

    /// Drive the same original three-dimensional split with held acceleration or
    /// force density. Source work stays separate from molecular and wall energy.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_forcing(
        &mut self,
        maximum_duration: Time,
        controls: MolecularStepControl,
        force: BodyForce3d<'_>,
        force_controls: ForceControl3d,
        workspace: &mut ForcingWorkspace3d,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MolecularStepReport3d, DissipationError> {
        self.step_impl(
            maximum_duration,
            controls,
            Some((force, force_controls, workspace)),
            None,
            &[],
            None,
            cancelled,
        )
    }

    /// Compose held body forcing with stage-time reservoirs across the complete
    /// molecular--Euler--molecular split.
    ///
    /// Exact kinetic-work half-kicks surround the central Euler RK2 interval;
    /// acceleration uses each kick's actual mixture density. Original stress,
    /// species, enthalpy, heat and optional SGS operators sample the prescribed
    /// exterior gas throughout the two molecular half intervals. Force work and
    /// impulse retain their separate ledger and final Mach/acoustic admission.
    ///
    /// One provider budget covers every stage and rejected split; it is never
    /// reset on shortening. Numerical retries replay from the accepted state;
    /// source/budget refusal is terminal. Cells, clock, profiles, all loads and
    /// reports publish together. Zero forcing uses the original timed split.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_forcing_and_reservoir_profiles(
        &mut self,
        maximum_duration: Time,
        controls: MolecularStepControl,
        force: BodyForce3d<'_>,
        force_controls: ForceControl3d,
        workspace: &mut ForcingWorkspace3d,
        maximum_face_evaluations: usize,
        mut sample: impl FnMut(
            f64,
            BoundarySide,
            usize,
            [f64; 3],
            &mut [f64],
        ) -> Result<ReservoirFaceState, Error>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MolecularStepReport3d, DissipationError> {
        let mut context = Context::new(&mut sample, maximum_face_evaluations);
        self.step_impl(
            maximum_duration,
            controls,
            Some((force, force_controls, workspace)),
            None,
            &[],
            Some(&mut context),
            cancelled,
        )
    }

    /// The same forced physical split with native RKL2 molecular stage times.
    ///
    /// The face-provider budget spans both halves, central forced RK stages
    /// and all retries. The independent RKL operator budget also counts every
    /// discarded residual. Supplied forcing scratch counts even for zero force.
    /// Prescribed wall velocities are retained; finite inertia and accelerated
    /// walls remain outside the timed-profile composition.
    #[allow(clippy::too_many_arguments)]
    pub fn step_stabilized_with_forcing_and_reservoir_profiles(
        &mut self,
        maximum_duration: Time,
        controls: StabilizedMolecularStepControl,
        force: BodyForce3d<'_>,
        force_controls: ForceControl3d,
        workspace: &mut ForcingWorkspace3d,
        maximum_face_evaluations: usize,
        mut sample: impl FnMut(
            f64,
            BoundarySide,
            usize,
            [f64; 3],
            &mut [f64],
        ) -> Result<ReservoirFaceState, Error>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<StabilizedMolecularStepReport3d, DissipationError> {
        let mut context = Context::new(&mut sample, maximum_face_evaluations);
        let mut work = StepWork::new(controls.molecular);
        let transport = self.step_impl(
            maximum_duration,
            controls.split(),
            Some((force, force_controls, workspace)),
            Some(&mut work),
            &[],
            Some(&mut context),
            cancelled,
        )?;
        Ok(StabilizedMolecularStepReport3d {
            transport,
            operator_evaluations: work.evaluations,
            maximum_stages_used: work.maximum_stages_used,
            maximum_scaled_parabolic_number: work.maximum_scaled_parabolic_number,
        })
    }

    /// Drive existing no-slip walls through their stage-resolved tangential
    /// acceleration history. The endpoint velocities publish only with the
    /// complete accepted fluid interval, including any shortening retries.
    pub fn step_with_wall_acceleration(
        &mut self,
        maximum_duration: Time,
        controls: MolecularStepControl,
        accelerations: &[WallAcceleration3d],
        cancelled: impl FnMut() -> bool,
    ) -> Result<MolecularStepReport3d, DissipationError> {
        self.step_impl(
            maximum_duration,
            controls,
            None,
            None,
            accelerations,
            None,
            cancelled,
        )
    }

    /// Compose held body forcing with accelerating walls. Boundary impulse and
    /// work retain their own ledger; body impulse and work remain source terms.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_forcing_and_wall_acceleration(
        &mut self,
        maximum_duration: Time,
        controls: MolecularStepControl,
        force: BodyForce3d<'_>,
        force_controls: ForceControl3d,
        workspace: &mut ForcingWorkspace3d,
        accelerations: &[WallAcceleration3d],
        cancelled: impl FnMut() -> bool,
    ) -> Result<MolecularStepReport3d, DissipationError> {
        self.step_impl(
            maximum_duration,
            controls,
            Some((force, force_controls, workspace)),
            None,
            accelerations,
            None,
            cancelled,
        )
    }

    /// Apply original RKL2 molecular halves around one conservative spatial
    /// Euler interval. Maxwell--Stefan, Soret/Dufour, enthalpy, heat, full stress
    /// and optional SGS keep their original face operators and signed loads.
    /// One actual-evaluation budget covers both halves and discarded retries.
    ///
    /// This path requires prescribed constant wall velocities. Finite inertia
    /// is refused before source work; accelerating walls retain the original
    /// stage-resolved SSP entry points. The physical macrostep cap remains an
    /// accuracy choice, and the linear stability extension does not establish
    /// nonlinear positivity or entropy decay. All owners publish atomically.
    pub fn step_stabilized(
        &mut self,
        maximum_duration: Time,
        controls: StabilizedMolecularStepControl,
        cancelled: impl FnMut() -> bool,
    ) -> Result<StabilizedMolecularStepReport3d, DissipationError> {
        let mut work = StepWork::new(controls.molecular);
        let transport = self.step_impl(
            maximum_duration,
            controls.split(),
            None,
            Some(&mut work),
            &[],
            None,
            cancelled,
        )?;
        Ok(StabilizedMolecularStepReport3d {
            transport,
            operator_evaluations: work.evaluations,
            maximum_stages_used: work.maximum_stages_used,
            maximum_scaled_parabolic_number: work.maximum_scaled_parabolic_number,
        })
    }

    /// Retain original body impulse/work and final Mach/acoustic admission
    /// between stabilized halves. Supplied forcing scratch counts in full even
    /// when the authored force is zero.
    #[allow(clippy::too_many_arguments)]
    pub fn step_stabilized_with_forcing(
        &mut self,
        maximum_duration: Time,
        controls: StabilizedMolecularStepControl,
        force: BodyForce3d<'_>,
        force_controls: ForceControl3d,
        workspace: &mut ForcingWorkspace3d,
        cancelled: impl FnMut() -> bool,
    ) -> Result<StabilizedMolecularStepReport3d, DissipationError> {
        let mut work = StepWork::new(controls.molecular);
        let transport = self.step_impl(
            maximum_duration,
            controls.split(),
            Some((force, force_controls, workspace)),
            Some(&mut work),
            &[],
            None,
            cancelled,
        )?;
        Ok(StabilizedMolecularStepReport3d {
            transport,
            operator_evaluations: work.evaluations,
            maximum_stages_used: work.maximum_stages_used,
            maximum_scaled_parabolic_number: work.maximum_scaled_parabolic_number,
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn molecular_half(
        &mut self,
        duration: Time,
        controls: DissipationControl,
        stabilization: Option<&mut StepWork>,
        drive: Option<(&WallDrive3d, f64)>,
        previous: Option<&WallInertiaReport3d>,
        wall_work: &WallWork,
        profiles: Option<(&mut Context<'_>, timed::Clock)>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MolecularReport, DissipationError> {
        match stabilization {
            None => self.molecular.advance_driven_with_profiles_unpublished(
                &mut self.trial,
                duration,
                controls,
                self.maximum_workspace_bytes,
                drive,
                previous,
                wall_work,
                profiles,
                cancelled,
            ),
            Some(work) => {
                let report = self
                    .molecular
                    .advance_stabilized_with_profiles_unpublished(
                        &mut self.trial,
                        duration,
                        work.controls,
                        self.maximum_workspace_bytes,
                        &mut work.evaluations,
                        profiles,
                        cancelled,
                    )?;
                work.record(
                    report.maximum_stages_used,
                    report.maximum_scaled_parabolic_number,
                );
                Ok(report.transport)
            }
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn step_impl(
        &mut self,
        maximum_duration: Time,
        controls: MolecularStepControl,
        forcing: Option<forcing::Context<'_, '_>>,
        stabilization: Option<&mut StepWork>,
        accelerations: &[WallAcceleration3d],
        profiles: Option<&mut Context<'_>>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MolecularStepReport3d, DissipationError> {
        poll(&mut cancelled)?;
        #[cfg(feature = "experimental-sgs")]
        if let Some(mut filter) = self
            .molecular
            .sgs
            .as_mut()
            .and_then(sgs::SubgridTransport::take_lagrangian_filter)
        {
            // The trial borrows only this detached filter, leaving the complete
            // conservative owner available. Its capacity is still budgeted.
            // No fallible action follows complete flow/history acceptance.
            let result = (|| {
                let trial = self
                    .molecular
                    .sgs
                    .as_mut()
                    .expect("the Lagrangian mixture owner was just selected")
                    .prepare_lagrangian_with_walls(
                        &mut filter,
                        &self.accepted,
                        &self.molecular.walls,
                        &mut cancelled,
                    )?;
                let fit = trial.reports();
                let report = self.step_prepared(
                    maximum_duration,
                    controls,
                    forcing,
                    stabilization,
                    accelerations,
                    profiles,
                    &mut cancelled,
                )?;
                trial.commit();
                Ok((report, fit))
            })();
            let committed_fit = result.as_ref().ok().map(|(_, fit)| *fit);
            self.molecular
                .sgs
                .as_mut()
                .expect("a coupled step retains its original mixture SGS owner")
                .restore_lagrangian_filter(filter, committed_fit);
            return result.map(|(report, _)| report);
        }
        self.step_prepared(
            maximum_duration,
            controls,
            forcing,
            stabilization,
            accelerations,
            profiles,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn step_prepared(
        &mut self,
        maximum_duration: Time,
        controls: MolecularStepControl,
        mut forcing: Option<forcing::Context<'_, '_>>,
        mut stabilization: Option<&mut StepWork>,
        accelerations: &[WallAcceleration3d],
        mut profiles: Option<&mut Context<'_>>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MolecularStepReport3d, DissipationError> {
        poll(&mut cancelled)?;
        if profiles.is_some() {
            timed::admit(
                &self.accepted,
                timed::Clock::direct(self.accepted.time, maximum_duration.value()),
                maximum_duration.value(),
                self.molecular.inertia.is_some(),
            )?;
            if !accelerations.is_empty() {
                return Err(DissipationError::Invalid(
                    "time-resolved molecular reservoirs do not compose wall acceleration",
                ));
            }
        }
        if stabilization.is_some()
            && (self.molecular.inertia.is_some() || !accelerations.is_empty())
        {
            return Err(DissipationError::Invalid(
                "stabilized molecular transport requires prescribed walls without finite inertia or acceleration",
            ));
        }
        let drive = self
            .molecular
            .prepare_wall_drive(accelerations, &mut cancelled)?;
        let forcing_bytes = forcing
            .as_ref()
            .map_or(0, |(_, _, work)| work.workspace_bytes());
        let active = if let Some((force, control, workspace)) = forcing.as_mut() {
            control.admit::<Error>()?;
            workspace.admit(&self.accepted)?;
            force.admit::<Error>(self.accepted.cell_count(), &mut cancelled)?
        } else {
            false
        };
        if !active {
            forcing = None;
        }
        let bytes =
            self.workspace_bytes()
                .checked_add(forcing_bytes)
                .ok_or(DissipationError::Invalid(
                    "coupled spatial mixture workspace size overflow",
                ))?;
        admit_workspace(bytes, self.maximum_workspace_bytes)?;
        controls.molecular.validate()?;
        if !positive(maximum_duration.value())
            || controls.maximum_coupling_retries > 32
            || !positive(controls.advection.cfl)
            || controls.advection.cfl > 0.45
            || controls.advection.maximum_halvings > 32
        {
            return Err(DissipationError::Invalid(
                "coupled spatial mixture duration or controls",
            ));
        }
        if let Some(work) = stabilization.as_ref() {
            work.controls.validate()?;
            // The original operator admits its own field and scratch. Preserve
            // the accepted field, complete load buffers and all supplied force
            // scratch in the same complete-owner budget before lazy allocation.
            let outside = bytes
                .checked_sub(self.trial.workspace_bytes())
                .and_then(|value| value.checked_sub(self.molecular.workspace_bytes()))
                .ok_or(DissipationError::Invalid(
                    "coupled stabilized spatial mixture workspace size overflow",
                ))?;
            let allowance = self
                .maximum_workspace_bytes
                .min(super::super::MAX_MIXTURE_WORKSPACE_BYTES)
                .checked_sub(outside)
                .ok_or(DissipationError::Invalid(
                    "coupled stabilized spatial mixture total workspace budget",
                ))?;
            self.molecular
                .prepare_stabilized(&self.trial, true, allowance)?;
        }
        let n = self.accepted.gas.species_count();
        let width = n + 4;
        let before = checked_totals(
            &self.accepted.states,
            &self.accepted.grid,
            n,
            &mut cancelled,
        )?;
        let mut rate: f64 = 0.0;
        for row in self.accepted.states.chunks_exact(width) {
            let p = cell(row, n).primitive(&self.accepted.gas, &mut cancelled)?;
            if let Some((_, control, _)) = forcing.as_ref() {
                control.admit_velocity::<Error>(p.velocity, p.sound_speed)?;
            }
            rate = rate.max(forcing::acoustic_rate(p, &self.accepted.grid)?);
        }
        if let Some(context) = profiles.as_deref_mut() {
            context.sample_at(
                self.accepted.time,
                &self.accepted.gas,
                &self.accepted.grid,
                &mut self.trial.reservoirs,
                &mut cancelled,
            )?;
            rate = rate.max(
                self.trial
                    .reservoirs
                    .maximum_rate_at(&self.accepted.grid, true),
            );
        }
        // Evolved source stages and faces remain authoritative; this initial
        // headroom avoids a predictably excessive first molecular interval.
        let mut dt = maximum_duration
            .value()
            .min(0.8 * controls.advection.cfl / rate);
        let mut last = DissipationError::Invalid("no coupled spatial mixture attempt");
        let mut last_dt = dt;
        let wall_work = WallWork::default();
        for retry in 0..=controls.maximum_coupling_retries {
            poll(&mut cancelled)?;
            if let Some(work) = stabilization.as_deref_mut() {
                work.reset_extrema();
            }
            if !positive(0.5 * dt)
                || !(self.accepted.time + dt).is_finite()
                || self.accepted.time + dt <= self.accepted.time
            {
                return Err(DissipationError::Invalid(
                    "unrepresentable coupled spatial mixture increment",
                ));
            }
            copy_physical(&self.accepted, &mut self.trial, &mut cancelled)?;
            let mut shortened = None;
            let mut endpoint_walls = None;
            let attempt = (|| -> Result<MolecularStepReport3d, DissipationError> {
                endpoint_walls = drive.as_ref().map(|d| d.at(dt)).transpose()?;
                let first = self.molecular_half(
                    Time::new(0.5 * dt),
                    controls.molecular,
                    stabilization.as_deref_mut(),
                    drive.as_ref().map(|d| (d, 0.0)),
                    None,
                    &wall_work,
                    profiles.as_deref_mut().map(|context| {
                        (
                            context,
                            timed::Clock {
                                origin_s: self.accepted.time,
                                phase_start_s: 0.0,
                                phase_end_s: 0.5 * dt,
                            },
                        )
                    }),
                    &mut cancelled,
                )?;
                copy_values(
                    self.molecular.trial_impulses(),
                    &mut self.boundary_trial,
                    &mut cancelled,
                )?;
                // Admit the immutable source/context once at the outer owner,
                // preserving native stage checks inside the original Euler step.
                let advection = if let Some(context) = profiles.as_deref_mut() {
                    self.trial.step_impl_with_contact(
                        Time::new(dt),
                        controls.advection,
                        forcing.as_mut().map(|(force, control, workspace)| {
                            (*force, *control, &mut **workspace)
                        }),
                        None,
                        Some(context),
                        &mut cancelled,
                    )?
                } else if let Some((force, control, workspace)) = forcing.as_mut() {
                    self.trial.step_impl(
                        Time::new(dt),
                        controls.advection,
                        Some((*force, *control, &mut **workspace)),
                        &mut cancelled,
                    )?
                } else {
                    self.trial
                        .step(Time::new(dt), controls.advection, &mut cancelled)?
                };
                if advection.elapsed.value() != dt {
                    shortened = Some(advection.elapsed.value());
                    return Err(DissipationError::Invalid(
                        "advection shortened coupled spatial mixture interval",
                    ));
                }
                add_impulses(
                    &self.trial.boundary_impulses,
                    &mut self.boundary_trial,
                    &mut cancelled,
                )?;
                let second = self.molecular_half(
                    Time::new(0.5 * dt),
                    controls.molecular,
                    stabilization.as_deref_mut(),
                    // The molecular clock covers [0,h/2] then [h/2,h].
                    // Advection uses the same fixed geometry and no wall work.
                    drive.as_ref().map(|d| (d, 0.5 * dt)),
                    first.wall_inertia.as_ref(),
                    &wall_work,
                    profiles.as_deref_mut().map(|context| {
                        (
                            context,
                            timed::Clock {
                                origin_s: self.accepted.time,
                                phase_start_s: 0.5 * dt,
                                phase_end_s: dt,
                            },
                        )
                    }),
                    &mut cancelled,
                )?;
                add_impulses(
                    self.molecular.trial_impulses(),
                    &mut self.boundary_trial,
                    &mut cancelled,
                )?;
                if let Some((_, control, _)) = forcing.as_ref() {
                    let final_rate = forcing::admit_field(
                        &self.trial.gas,
                        &self.trial.grid,
                        &self.trial.states,
                        *control,
                        &mut cancelled,
                    )?;
                    stepping::check_cfl::<Error>(final_rate * dt, controls.advection.cfl)?;
                }
                let after =
                    checked_totals(&self.trial.states, &self.trial.grid, n, &mut cancelled)?;
                let outward = sum_rows(&self.boundary_trial, width, 1.0, &mut cancelled)?;
                let boundary_change = Totals::from_row(&outward.map(|value| -value), n);
                let source_change = advection.source_change;
                let mut residual = [0.0; WIDTH_LIMIT];
                for (d, value) in residual[..width].iter_mut().enumerate() {
                    *value =
                        after.component(d) - before.component(d) - boundary_change.component(d);
                }
                for axis in 0..3 {
                    residual[n + axis] -= source_change.impulse[axis];
                }
                residual[n + 3] -= source_change.work;
                if !residual[..width].iter().all(|value| value.is_finite()) {
                    return Err(DissipationError::Invalid(
                        "coupled spatial mixture balance overflow",
                    ));
                }
                let conservation_residual = Totals::from_row(&residual[..width], n);
                let absolute_energy_residual_j =
                    conservation_residual.absolute_energy(&self.trial.gas)?;
                poll(&mut cancelled)?;
                Ok(MolecularStepReport3d {
                    elapsed: Time::new(dt),
                    coupling_retries: retry,
                    first,
                    advection,
                    second,
                    before,
                    after,
                    boundary_change,
                    source_change,
                    conservation_residual,
                    absolute_energy_residual_j,
                })
            })();
            match attempt {
                Ok(report) => {
                    std::mem::swap(&mut self.accepted, &mut self.trial);
                    std::mem::swap(&mut self.boundary_impulses, &mut self.boundary_trial);
                    self.molecular.publish(report.second, endpoint_walls);
                    self.last_step = Some(report);
                    return Ok(report);
                }
                Err(DissipationError::Cancelled) => return Err(DissipationError::Cancelled),
                Err(error) if profiles.as_ref().is_some_and(|context| context.failed()) => {
                    return Err(error);
                }
                Err(error) if stabilization.is_some() && error.is_terminal() => {
                    return Err(error.step_failure(0.0, dt));
                }
                Err(error) => {
                    last = error;
                    last_dt = dt;
                    dt = shortened.unwrap_or(0.5 * dt);
                }
            }
        }
        Err(DissipationError::StepFailure {
            elapsed_time_s: 0.0,
            step_s: last_dt,
            last: Box::new(last),
        })
    }
}

fn copy_physical(
    source: &MixtureEuler3d,
    target: &mut MixtureEuler3d,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), DissipationError> {
    copy_values(&source.states, &mut target.states, cancelled)?;
    copy_values(
        &source.boundary_impulses,
        &mut target.boundary_impulses,
        cancelled,
    )?;
    match (&source.reservoirs.profiles, &mut target.reservoirs.profiles) {
        (Some(source), Some(target)) => target.copy_active_from(source, cancelled)?,
        (None, None) => {}
        _ => {
            return Err(DissipationError::Invalid(
                "coupled mixture reservoir profile layout differs",
            ));
        }
    }
    // A previous accepted or discarded moving split may leave the reusable
    // frame on another geometry. Every transaction starts from the owner.
    target.grid = source.grid;
    target.time = source.time;
    target.last_step = source.last_step;
    Ok(())
}

fn add_impulses(
    source: &[f64],
    target: &mut [f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), DissipationError> {
    for (source, target) in source.chunks(4096).zip(target.chunks_mut(4096)) {
        poll(cancelled)?;
        for (source, target) in source.iter().zip(target) {
            *target += source;
        }
    }
    Ok(())
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
