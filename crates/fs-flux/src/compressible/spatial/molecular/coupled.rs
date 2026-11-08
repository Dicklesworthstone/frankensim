use super::*;
use crate::compressible::mixture::molecular::stabilized::IntervalFailure;
use crate::compressible::molecular::stabilized::StepWork;

#[cfg(feature = "experimental-sgs")]
mod lagrangian_channel;

/// One completely accepted three-dimensional molecular/Euler/molecular step.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MolecularStepReport3d {
    /// Physical duration shared by all stages [s].
    pub elapsed: Time,
    /// Rejected complete trials before acceptance.
    pub coupling_retries: usize,
    /// First molecular half step.
    pub first: MolecularReport3d,
    /// Original Euler transport report between molecular half steps.
    pub advection: StepReport,
    /// Second molecular half step.
    pub second: MolecularReport3d,
    /// Original complete field totals.
    pub before: Totals,
    /// Accepted complete field totals.
    pub after: Totals,
    /// Combined inward molecular and advective boundary transfer.
    pub boundary_change: Totals,
    /// Accepted body impulse/work from the central forced Euler interval.
    pub source_change: ForceChange3d,
    /// after - before - boundary_change - source_change for all five components.
    pub conservation_residual: Totals,
}

/// Complete three-dimensional coupled interval and actual RKL2 work.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StabilizedMolecularStepReport3d {
    /// Original accepted molecular/advection/molecular physical balances.
    pub transport: MolecularStepReport3d,
    /// All molecular evaluations across both halves and discarded whole trials.
    pub operator_evaluations: usize,
    /// Largest stage count in the final accepted pair of molecular halves.
    pub maximum_stages_used: usize,
    /// Largest scaled parabolic number in the final accepted molecular halves.
    pub maximum_scaled_parabolic_number: f64,
}

/// Atomic source-driven three-dimensional compressible viscous/thermal flow.
///
/// A symmetric molecular/advection/molecular composition retains one actual
/// accepted duration. A shorter Euler step or a refused later stage restarts
/// the complete trial from the accepted field. Molecular halves use SSP RK2
/// by default or opt-in RKL2; smooth MUSCL Euler permits second-order splitting,
/// while limiting or
/// first-order fallback reduces order. No entropy, turbulence, physical shock
/// thickness or experimental-accuracy certificate is implied.
pub struct MolecularEuler3d<L> {
    accepted: Euler3d,
    trial: Euler3d,
    molecular: MolecularTransport3d<L>,
    boundary_impulses: Vec<BoundaryImpulse>,
    boundary_trial: Vec<BoundaryImpulse>,
    last_step: Option<MolecularStepReport3d>,
}

impl<L: TransportLaw> MolecularEuler3d<L> {
    /// Enable the original molecular bath law in both split half steps.
    ///
    /// Complete source-admitted uniform or face-resolved Euler reservoirs
    /// supply the held bath traces. Unbound open sides retain zero molecular
    /// flux. This selection survives accepted/trial swaps and forced or
    /// stabilized steps; timed step APIs supply the physical stage resampling.
    /// Changing the selection invalidates the prior complete-step report.
    pub fn with_reservoir_exchange(mut self) -> Self {
        if !self.molecular.reservoir_exchange_enabled() {
            self.last_step = None;
        }
        self.molecular = self.molecular.with_reservoir_exchange();
        self
    }

    /// Whether complete reservoirs supply heat and molecular traction.
    pub const fn reservoir_exchange_enabled(&self) -> bool {
        self.molecular.reservoir_exchange_enabled()
    }

    /// Enable local SGS stress and heat at complete reservoirs in both halves.
    /// Uses the same admitted bath traces as molecular exchange, including
    /// timed-stage sampling. Dynamic and Lagrangian closures refuse.
    /// Changing this selection invalidates the previous complete-step report.
    #[cfg(feature = "experimental-sgs")]
    pub fn with_local_sgs_reservoir_exchange(mut self) -> Result<Self, MolecularError> {
        let changed = !self.molecular.local_sgs_reservoir_exchange_enabled();
        self.molecular = self.molecular.with_local_sgs_reservoir_exchange()?;
        if changed {
            self.last_step = None;
        }
        Ok(self)
    }

    #[cfg(feature = "experimental-sgs")]
    pub fn local_sgs_reservoir_exchange_enabled(&self) -> bool {
        self.molecular.local_sgs_reservoir_exchange_enabled()
    }

    /// Apply the same physical wall to both molecular half steps. Euler keeps
    /// its impermeable pressure boundary; wall work/heat enter the combined
    /// conservative exchange report only through the molecular face flux.
    /// Changing a wall invalidates the previous complete-step report so its
    /// accepted impulse cannot be relabeled with a different wall condition.
    pub fn with_wall(
        mut self,
        side: BoundarySide,
        wall: MolecularWall3d,
    ) -> Result<Self, MolecularError> {
        let changed = self.molecular.wall(side) != wall;
        self.molecular = self.molecular.with_wall(side, wall)?;
        if changed {
            self.last_step = None;
        }
        Ok(self)
    }

    pub fn wall(&self, side: BoundarySide) -> MolecularWall3d {
        self.molecular.wall(side)
    }

    /// Apply the same equilibrium SGS wall stress in both molecular half
    /// steps, retaining the original complete-step wall-work ledger.
    #[cfg(feature = "experimental-sgs")]
    pub fn with_equilibrium_sgs_wall(
        mut self,
        side: BoundarySide,
        config: sgs::EquilibriumWallConfig,
    ) -> Result<Self, MolecularError> {
        let changed = self.molecular.equilibrium_sgs_wall(side) != Some(config);
        self.molecular = self.molecular.with_equilibrium_sgs_wall(side, config)?;
        if changed {
            self.last_step = None;
        }
        Ok(self)
    }

    #[cfg(feature = "experimental-sgs")]
    pub fn equilibrium_sgs_wall(&self, side: BoundarySide) -> Option<sgs::EquilibriumWallConfig> {
        self.molecular.equilibrium_sgs_wall(side)
    }

    /// Apply the explicit equilibrium stress/heat law to both molecular half
    /// steps and retain its heat in the complete-step boundary energy ledger.
    #[cfg(feature = "experimental-sgs")]
    pub fn with_equilibrium_sgs_thermal_wall(
        mut self,
        side: BoundarySide,
        config: sgs::EquilibriumThermalWallConfig,
    ) -> Result<Self, MolecularError> {
        let changed = self.molecular.equilibrium_sgs_thermal_wall(side) != Some(config);
        self.molecular = self.molecular.with_equilibrium_sgs_thermal_wall(side, config)?;
        if changed {
            self.last_step = None;
        }
        Ok(self)
    }

    #[cfg(feature = "experimental-sgs")]
    pub fn equilibrium_sgs_thermal_wall(&self, side: BoundarySide) -> Option<sgs::EquilibriumThermalWallConfig> {
        self.molecular.equilibrium_sgs_thermal_wall(side)
    }

    /// Compose explicit subgrid stress and heat with the original Euler
    /// advection and atomic retry mechanism. Every stage CFL includes SGS.
    #[cfg(feature = "experimental-sgs")]
    pub fn with_sgs(
        flow: Euler3d,
        law: L,
        config: sgs::SubgridConfig,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, MolecularError> {
        config.validate()?;
        let mut result = Self::new(flow, law, &mut cancelled)?;
        result.molecular.sgs = Some(sgs::SubgridTransport::new(&result.accepted, config)?);
        poll(&mut cancelled)?;
        Ok(result)
    }

    /// Fit periodic dynamic stress through both molecular halves and the same
    /// atomic Euler composition. The byte budget covers additional SGS arrays,
    /// excluding original molecular and Euler storage.
    #[cfg(feature = "experimental-sgs")]
    pub fn with_dynamic_sgs(
        flow: Euler3d,
        law: L,
        config: sgs::DynamicSubgridConfig,
        maximum_subgrid_workspace_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, MolecularError> {
        poll(&mut cancelled)?;
        let sgs =
            sgs::SubgridTransport::new_dynamic(&flow, config, maximum_subgrid_workspace_bytes)?;
        poll(&mut cancelled)?;
        let mut result = Self::new(flow, law, &mut cancelled)?;
        result.molecular.sgs = Some(sgs);
        poll(&mut cancelled)?;
        Ok(result)
    }

    /// Exact dynamic controls; `None` for molecular-only or local SGS flow.
    #[cfg(feature = "experimental-sgs")]
    pub fn dynamic_subgrid_config(&self) -> Option<sgs::DynamicSubgridConfig> {
        self.molecular.dynamic_subgrid_config()
    }

    /// Minimum additional numerical payload for Lagrangian SGS: 208 bytes/cell.
    /// Includes both history banks, physical-field staging and held coefficients;
    /// excludes the original flow/molecular workspaces and allocator metadata.
    #[cfg(feature = "experimental-sgs")]
    pub fn required_lagrangian_dynamic_sgs_workspace_bytes(
        flow: &Euler3d,
    ) -> Result<usize, MolecularError> {
        sgs::SubgridTransport::required_lagrangian_bytes(flow)
    }

    /// Add periodic Lagrangian Favre--Germano stress and fixed-Prandtl heat.
    ///
    /// One fit samples the accepted physical START state and is held through
    /// both molecular halves, all SSP/RKL stages and every coupling retry. The
    /// fit uses physical seconds, and its history commits only when the whole
    /// step accepts. Rejected steps cannot age or replace accepted memory.
    /// Face density/strain and all transport stability/ratio checks remain
    /// stage-local. The held coefficient and pathline update are first order
    /// in time, including when the underlying flow integrator is second order.
    ///
    /// The budget covers all additional SGS Vec capacities. No filter/history
    /// allocation occurs during stepping. This constructor starts a fresh
    /// history on the supplied trajectory; it does not import earlier memory.
    #[cfg(feature = "experimental-sgs")]
    pub fn with_lagrangian_dynamic_sgs(
        flow: Euler3d,
        law: L,
        config: sgs::LagrangianSubgridConfig,
        maximum_subgrid_workspace_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, MolecularError> {
        poll(&mut cancelled)?;
        let sgs =
            sgs::SubgridTransport::new_lagrangian(&flow, config, maximum_subgrid_workspace_bytes)?;
        poll(&mut cancelled)?;
        let mut result = Self::new(flow, law, &mut cancelled)?;
        result.molecular.sgs = Some(sgs);
        poll(&mut cancelled)?;
        Ok(result)
    }

    /// Exact Lagrangian controls; other closure choices return None.
    #[cfg(feature = "experimental-sgs")]
    pub fn lagrangian_dynamic_subgrid_config(&self) -> Option<sgs::LagrangianSubgridConfig> {
        self.molecular
            .sgs
            .as_ref()
            .and_then(sgs::SubgridTransport::lagrangian_config)
    }

    /// Last committed history SAMPLE time, at the start of the accepted step.
    /// This is not the flow endpoint. None before the first successful step or
    /// when the selected closure does not use Lagrangian history.
    #[cfg(feature = "experimental-sgs")]
    pub fn lagrangian_subgrid_sample_time(&self) -> Option<Time> {
        self.molecular
            .sgs
            .as_ref()
            .and_then(sgs::SubgridTransport::lagrangian_sample_time)
    }

    /// Last committed start-state fit. Cell-viscosity diagnostics describe that
    /// sample; actual face viscosities and bounds use each stage's face strain.
    #[cfg(feature = "experimental-sgs")]
    pub fn lagrangian_subgrid_report(&self) -> Option<sgs::DynamicSgsReport> {
        self.molecular
            .sgs
            .as_ref()
            .and_then(sgs::SubgridTransport::lagrangian_report)
    }

    /// Additional numerical payload for paired Lagrangian momentum/heat:
    /// 320 bytes/cell, excluding original flow and molecular workspaces.
    #[cfg(feature = "experimental-sgs")]
    pub fn required_lagrangian_dynamic_thermal_sgs_workspace_bytes(
        flow: &Euler3d,
    ) -> Result<usize, MolecularError> {
        sgs::SubgridTransport::required_lagrangian_thermal_bytes(flow)
    }

    /// Fit independent stress and temperature-flux memories at physical step
    /// starts. Both histories publish together only after complete acceptance.
    ///
    /// This shares the fixed-Prandtl path's split/retry/stabilization transaction.
    /// Conservative faces evaluate current density, strain and source Cp; no
    /// turbulent Prandtl number, extra heating source or stage-time memory is
    /// introduced. Step and vorticity queries reuse prepared numerical storage.
    /// The held coefficient/pathline update remains first order in time, and
    /// requires a fully periodic Cartesian domain with resolved smooth fields.
    #[cfg(feature = "experimental-sgs")]
    pub fn with_lagrangian_dynamic_thermal_sgs(
        flow: Euler3d,
        law: L,
        config: sgs::LagrangianThermalSubgridConfig,
        maximum_subgrid_workspace_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, MolecularError> {
        poll(&mut cancelled)?;
        let sgs = sgs::SubgridTransport::new_lagrangian_thermal(
            &flow,
            config,
            maximum_subgrid_workspace_bytes,
        )?;
        poll(&mut cancelled)?;
        let mut result = Self::new(flow, law, &mut cancelled)?;
        result.molecular.sgs = Some(sgs);
        poll(&mut cancelled)?;
        Ok(result)
    }

    #[cfg(feature = "experimental-sgs")]
    pub fn lagrangian_thermal_subgrid_config(&self) -> Option<sgs::LagrangianThermalSubgridConfig> {
        self.molecular
            .sgs
            .as_ref()
            .and_then(sgs::SubgridTransport::lagrangian_thermal_config)
    }

    /// Last committed scalar fit at lagrangian_subgrid_sample_time().
    /// None before acceptance or when independent heat memory is not selected.
    #[cfg(feature = "experimental-sgs")]
    pub fn lagrangian_subgrid_heat_report(&self) -> Option<sgs::DynamicScalarReport> {
        self.molecular
            .sgs
            .as_ref()
            .and_then(sgs::SubgridTransport::lagrangian_heat_report)
    }

    /// Fit momentum and temperature transport independently through both
    /// molecular halves of the existing atomic conservative composition.
    /// The byte budget covers additional SGS arrays, excluding original flow
    /// and molecular workspaces. No fixed turbulent Prandtl number is used.
    #[cfg(feature = "experimental-sgs")]
    pub fn with_dynamic_thermal_sgs(
        flow: Euler3d,
        law: L,
        config: sgs::DynamicThermalSubgridConfig,
        maximum_subgrid_workspace_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, MolecularError> {
        poll(&mut cancelled)?;
        let sgs = sgs::SubgridTransport::new_dynamic_thermal(
            &flow,
            config,
            maximum_subgrid_workspace_bytes,
        )?;
        poll(&mut cancelled)?;
        let mut result = Self::new(flow, law, &mut cancelled)?;
        result.molecular.sgs = Some(sgs);
        poll(&mut cancelled)?;
        Ok(result)
    }

    /// Exact dynamic momentum/heat controls; other closure choices return None.
    #[cfg(feature = "experimental-sgs")]
    pub fn dynamic_thermal_subgrid_config(&self) -> Option<sgs::DynamicThermalSubgridConfig> {
        self.molecular.dynamic_thermal_subgrid_config()
    }

    /// Fit independent channel stress/heat through both original molecular
    /// halves. Explicit lower/upper no-slip walls retain molecular wall heat,
    /// prescribed wall work and the complete conservative interval transaction.
    #[cfg(feature = "experimental-sgs")]
    pub fn with_channel_dynamic_thermal_sgs(
        flow: Euler3d,
        law: L,
        config: sgs::ChannelDynamicThermalSubgridConfig,
        walls: [MolecularWall3d; 2],
        maximum_subgrid_workspace_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, MolecularError> {
        poll(&mut cancelled)?;
        config.validate()?;
        let normal = sgs::channel_axis(&flow.grid, flow.boundaries)?;
        let required =
            MolecularTransport3d::<L>::required_channel_dynamic_thermal_sgs_workspace_bytes(&flow)?;
        if required > maximum_subgrid_workspace_bytes {
            return Err(fs_lbm::sgs3::dynamic::DynamicSgsError::WorkspaceBytes {
                required,
                maximum: maximum_subgrid_workspace_bytes,
            }
            .into());
        }
        let mut result = Self::new(flow, law, &mut cancelled)?;
        for (side, wall) in walls.into_iter().enumerate() {
            result = result.with_wall(BoundarySide::ALL[2 * normal + side], wall)?;
        }
        result.molecular.sgs = Some(sgs::SubgridTransport::new_channel_dynamic_thermal(
            &result.accepted,
            config,
            &result.molecular.walls,
            maximum_subgrid_workspace_bytes,
        )?);
        poll(&mut cancelled)?;
        Ok(result)
    }

    /// Exact channel controls; other closure geometries return None.
    #[cfg(feature = "experimental-sgs")]
    pub fn channel_dynamic_thermal_subgrid_config(
        &self,
    ) -> Option<sgs::ChannelDynamicThermalSubgridConfig> {
        self.molecular.channel_dynamic_thermal_subgrid_config()
    }

    /// Exact local controls; `None` for molecular-only or dynamic SGS flow.
    #[cfg(feature = "experimental-sgs")]
    pub fn subgrid_config(&self) -> Option<sgs::SubgridConfig> {
        self.molecular.subgrid_config()
    }

    /// Bind the identical original source while preserving the complete field.
    /// This owns a second complete Euler field, including any prepared face
    /// profiles. Each profile copy retains its original per-owner byte limit;
    /// [`Self::reservoir_profile_workspace_bytes`] reports their combined cost.
    /// To admit the combined profile storage before allocating it, construct
    /// from an unprofiled field and use [`Self::prepare_reservoir_profiles`].
    pub fn new(
        flow: Euler3d,
        law: L,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, MolecularError> {
        poll(&mut cancelled)?;
        let molecular = MolecularTransport3d::new(&flow, law, &mut cancelled)?;
        let mut cells = allocated(flow.states.len(), State::default())?;
        cells.copy_from_slice(&flow.states);
        let mut trial = Euler3d::new(
            flow.gas.clone(),
            flow.grid,
            cells,
            flow.boundaries,
            &mut cancelled,
        )?;
        trial.reservoirs = flow.reservoirs.try_clone(&mut cancelled)?;
        if let Some(balance) = &flow.hydrostatic {
            trial.bind_hydrostatic_reference(
                balance.axis,
                balance.reference,
                balance.quadrature,
                &mut cancelled,
            )?;
        }
        let faces = flow.boundary_impulses.len();
        let result = Self {
            accepted: flow,
            trial,
            molecular,
            boundary_impulses: allocated(faces, Totals::default())?,
            boundary_trial: allocated(faces, Totals::default())?,
            last_step: None,
        };
        poll(&mut cancelled)?;
        Ok(result)
    }

    /// Accepted field and physical clock. Its Euler report describes the middle
    /// advective stage; this owner's report and loads cover the complete split.
    pub const fn flow(&self) -> &Euler3d {
        &self.accepted
    }

    /// Minimum additional profile bytes across both complete Euler owners.
    pub fn required_reservoir_profile_workspace_bytes(
        &self,
        sides: &[BoundarySide],
    ) -> Result<usize, MolecularError> {
        self.accepted
            .required_reservoir_profile_workspace_bytes(sides)?
            .checked_mul(2)
            .ok_or(MolecularError::Invalid(
                "coupled reservoir profile size overflow",
            ))
    }

    /// Prepare identical held exterior faces in both accepted and trial fields.
    /// The byte allowance covers BOTH owners' active and candidate arrays.
    /// Callbacks run once, on the original grid/source; no resampling occurs
    /// during molecular halves, Euler stages or retries. Molecular exchange is
    /// zero by default and can use these same held traces through the explicit
    /// [`Self::with_reservoir_exchange`] selection.
    pub fn prepare_reservoir_profiles(
        &mut self,
        sides: &[BoundarySide],
        maximum_profile_bytes: usize,
        sample: impl FnMut(BoundarySide, usize, [f64; 3]) -> Result<ReservoirFaceState, Error>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), MolecularError> {
        let profiles = self.accepted.reservoirs.prepare_profiles(
            &self.accepted.gas,
            &self.accepted.grid,
            self.accepted.boundaries,
            sides,
            maximum_profile_bytes / 2,
            sample,
            &mut cancelled,
        )?;
        let trial_profiles = profiles.try_clone(&mut cancelled)?;
        poll(&mut cancelled)?;
        self.accepted.reservoirs.profiles = profiles;
        self.trial.reservoirs.profiles = trial_profiles;
        self.accepted.last_step = None;
        self.trial.last_step = None;
        self.last_step = None;
        Ok(())
    }

    /// Replace both complete fields' held traces in one control transaction.
    /// Every face is sampled once into reusable scratch. A failed sample,
    /// original-source refusal or cancellation preserves both active copies
    /// and all accepted flow state. Successful updates invalidate old reports
    /// but leave the physical clock, cells and historical impulses unchanged.
    pub fn update_reservoir_profiles(
        &mut self,
        sample: impl FnMut(BoundarySide, usize, [f64; 3]) -> Result<ReservoirFaceState, Error>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), MolecularError> {
        self.accepted.reservoirs.profiles.stage(
            &self.accepted.gas,
            &self.accepted.grid,
            sample,
            &mut cancelled,
        )?;
        self.trial
            .reservoirs
            .profiles
            .copy_staged_from(&self.accepted.reservoirs.profiles, &mut cancelled)?;
        poll(&mut cancelled)?;
        self.accepted.reservoirs.profiles.publish();
        self.trial.reservoirs.profiles.publish();
        self.accepted.last_step = None;
        self.trial.last_step = None;
        self.last_step = None;
        Ok(())
    }

    /// Actual profile numerical capacity across accepted and trial owners.
    pub fn reservoir_profile_workspace_bytes(&self) -> usize {
        self.accepted.reservoir_profile_workspace_bytes()
            + self.trial.reservoir_profile_workspace_bytes()
    }

    /// Exact immutable source law, including retained material/collision data.
    pub fn law(&self) -> &L {
        self.molecular.law()
    }

    /// Required additional storage for the accepted field's complete stress curl.
    pub fn required_vorticity_workspace_bytes(&self) -> Result<usize, MolecularError> {
        self.molecular.required_vorticity_workspace_bytes()
    }

    /// Allocate reusable inviscid/stress observations under an explicit byte cap.
    pub fn vorticity_workspace(
        &self,
        maximum_bytes: usize,
    ) -> Result<StressVorticityWorkspace3d, MolecularError> {
        self.molecular.vorticity_workspace(maximum_bytes)
    }

    /// Observe the accepted state with its configured molecular/SGS face stress.
    ///
    /// Includes held wall and enabled reservoir traction; imposed body-force
    /// curl remains excluded. No accepted state, clock or step load is changed.
    /// A Lagrangian closure evaluates an uncommitted sample at this endpoint:
    /// observations neither age accepted history nor replace its public report.
    pub fn vorticity<'a>(
        &mut self,
        work: &'a mut StressVorticityWorkspace3d,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<&'a [StressVorticityBudget3d], MolecularError> {
        poll(&mut cancelled)?;
        #[cfg(feature = "experimental-sgs")]
        if let Some(mut filter) = self
            .molecular
            .sgs
            .as_mut()
            .and_then(sgs::SubgridTransport::take_lagrangian_filter)
        {
            let result = match self
                .molecular
                .sgs
                .as_mut()
                .expect("the Lagrangian SGS owner was just selected")
                .prepare_lagrangian_with_walls(
                    &mut filter,
                    &self.accepted,
                    &self.molecular.walls,
                    &mut cancelled,
                ) {
                Ok(trial) => {
                    let observation =
                        self.molecular
                            .vorticity(&self.accepted, work, &mut cancelled);
                    drop(trial);
                    observation
                }
                Err(error) => Err(error),
            };
            self.molecular
                .sgs
                .as_mut()
                .expect("an observation retains its original SGS owner")
                .restore_lagrangian_filter(filter, None);
            return result;
        }
        self.molecular.vorticity(&self.accepted, work, cancelled)
    }

    /// Last completely accepted physical interval.
    pub const fn last_step(&self) -> Option<&MolecularStepReport3d> {
        self.last_step.as_ref()
    }

    /// Combined outward transfer, including both molecular support impulses.
    pub fn boundary_impulses(&self, side: BoundarySide) -> &[BoundaryImpulse] {
        &self.boundary_impulses[self.accepted.grid.boundary_range(side)]
    }

    /// Advance at most the requested duration, publishing only a complete split.
    pub fn step(
        &mut self,
        maximum_duration: Time,
        controls: MolecularStepControl,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MolecularStepReport3d, MolecularError> {
        self.step_impl(maximum_duration, controls, None, None, None, cancelled)
    }

    /// Compose molecular SSP and Euler with a replayable prescribed exterior.
    ///
    /// The first molecular half observes `[t, t + dt/2]`, native Euler observes
    /// `[t, t + dt]`, and the second half observes `[t + dt/2, t + dt]`.
    /// Each operator samples its own stage times; a shortened Euler step
    /// restarts the complete split from the accepted state. The provider must
    /// therefore permit repeated and nonmonotone time observations, following
    /// [`Euler3d::step_with_reservoir_profiles`].
    ///
    /// One cumulative face-evaluation budget includes the initial exterior CFL
    /// admission, both molecular halves, central Euler and every discarded
    /// trial. The original bath, wall, SGS and caloric operators remain active;
    /// molecular bath exchange requires [`Self::with_reservoir_exchange`].
    /// A successful split publishes the endpoint profile with all accepted
    /// cells, time, loads and reports. Any refusal leaves those data unchanged.
    /// This method does not advance an evolving random-source owner for you.
    pub fn step_with_reservoir_profiles(
        &mut self,
        maximum_duration: Time,
        controls: MolecularStepControl,
        maximum_face_evaluations: usize,
        mut sample: impl FnMut(f64, BoundarySide, usize, [f64; 3]) -> Result<ReservoirFaceState, Error>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MolecularStepReport3d, MolecularError> {
        let mut source = timed::TimedSource::new(&mut sample, maximum_face_evaluations);
        self.step_impl(
            maximum_duration,
            controls,
            None,
            None,
            Some(&mut source),
            cancelled,
        )
    }

    /// Sustain a driven viscous/channel flow with conservative body work.
    /// Molecular halves retain their current source, SGS and wall flux owners.
    /// A shortened or refused central force/advection/force interval retries
    /// the WHOLE molecular split from the accepted state with one duration.
    /// Allocate workspace from flow().forcing_workspace(maximum_bytes).
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_forcing(
        &mut self,
        maximum_duration: Time,
        controls: MolecularStepControl,
        force: BodyForce3d<'_>,
        force_controls: ForceControl3d,
        workspace: &mut ForcingWorkspace3d,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MolecularStepReport3d, MolecularError> {
        force_controls.admit::<Error>()?;
        let active = force.admit::<Error>(self.accepted.states.len(), &mut cancelled)?;
        workspace.admit(self.accepted.states.len())?;
        if !active {
            return self.step(maximum_duration, controls, cancelled);
        }
        self.step_impl(
            maximum_duration,
            controls,
            Some((force, force_controls, workspace)),
            None,
            None,
            cancelled,
        )
    }

    /// Drive the original viscous split with held forcing and a timed exterior.
    ///
    /// The physical stage clocks and one cumulative provider budget follow
    /// [`Self::step_with_reservoir_profiles`]. Central Euler retains its original
    /// conservative force half-kicks, kinetic work and authored Mach admission.
    /// Any shortened central interval replays both molecular halves and source
    /// kicks from the accepted field. Publication includes the endpoint profile
    /// and all source/boundary accounts only after the complete split succeeds.
    /// Zero forcing follows the exact unforced timed path.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_forcing_and_reservoir_profiles(
        &mut self,
        maximum_duration: Time,
        controls: MolecularStepControl,
        force: BodyForce3d<'_>,
        force_controls: ForceControl3d,
        workspace: &mut ForcingWorkspace3d,
        maximum_face_evaluations: usize,
        mut sample: impl FnMut(f64, BoundarySide, usize, [f64; 3]) -> Result<ReservoirFaceState, Error>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MolecularStepReport3d, MolecularError> {
        force_controls.admit::<Error>()?;
        let active = force.admit::<Error>(self.accepted.states.len(), &mut cancelled)?;
        workspace.admit(self.accepted.states.len())?;
        if !active {
            return self.step_with_reservoir_profiles(
                maximum_duration,
                controls,
                maximum_face_evaluations,
                sample,
                cancelled,
            );
        }
        let mut source = timed::TimedSource::new(&mut sample, maximum_face_evaluations);
        self.step_impl(
            maximum_duration,
            controls,
            Some((force, force_controls, workspace)),
            None,
            Some(&mut source),
            cancelled,
        )
    }

    /// Compose RKL2 molecular halves around the original volumetric Euler step.
    ///
    /// One cumulative evaluation budget covers both halves and discarded whole
    /// trials. Native wall, SGS, reservoir, source and hyperbolic CFL admission
    /// remain active. Refusal publishes no field, clock, load or report.
    pub fn step_stabilized(
        &mut self,
        maximum_duration: Time,
        controls: MolecularStepControl,
        stabilization: MolecularStabilization,
        cancelled: impl FnMut() -> bool,
    ) -> Result<StabilizedMolecularStepReport3d, MolecularError> {
        let mut work = StepWork::new(stabilization);
        let transport = self.step_impl(
            maximum_duration,
            controls,
            None,
            Some(&mut work),
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

    /// Use native RKL2 molecular stage times inside the complete prescribed split.
    ///
    /// Physical intervals and atomic publication follow
    /// [`Self::step_with_reservoir_profiles`]. The operator-evaluation budget
    /// and the provider-face budget each span both halves and every rejected
    /// whole split; neither is reset by a smaller trial. Original Legendre
    /// field and boundary weights are retained without per-stage allocation.
    #[allow(clippy::too_many_arguments)]
    pub fn step_stabilized_with_reservoir_profiles(
        &mut self,
        maximum_duration: Time,
        controls: MolecularStepControl,
        stabilization: MolecularStabilization,
        maximum_face_evaluations: usize,
        mut sample: impl FnMut(f64, BoundarySide, usize, [f64; 3]) -> Result<ReservoirFaceState, Error>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<StabilizedMolecularStepReport3d, MolecularError> {
        let mut source = timed::TimedSource::new(&mut sample, maximum_face_evaluations);
        let mut work = StepWork::new(stabilization);
        let transport = self.step_impl(
            maximum_duration,
            controls,
            None,
            Some(&mut work),
            Some(&mut source),
            cancelled,
        )?;
        Ok(StabilizedMolecularStepReport3d {
            transport,
            operator_evaluations: work.evaluations,
            maximum_stages_used: work.maximum_stages_used,
            maximum_scaled_parabolic_number: work.maximum_scaled_parabolic_number,
        })
    }

    /// Stabilize both molecular halves around native force/advection/force.
    /// Body work and wall transfer retain separate original physical accounts.
    #[allow(clippy::too_many_arguments)]
    pub fn step_stabilized_with_forcing(
        &mut self,
        maximum_duration: Time,
        controls: MolecularStepControl,
        stabilization: MolecularStabilization,
        force: BodyForce3d<'_>,
        force_controls: ForceControl3d,
        workspace: &mut ForcingWorkspace3d,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<StabilizedMolecularStepReport3d, MolecularError> {
        force_controls.admit::<Error>()?;
        let active = force.admit::<Error>(self.accepted.states.len(), &mut cancelled)?;
        workspace.admit(self.accepted.states.len())?;
        if !active {
            return self.step_stabilized(maximum_duration, controls, stabilization, cancelled);
        }
        let mut work = StepWork::new(stabilization);
        let transport = self.step_impl(
            maximum_duration,
            controls,
            Some((force, force_controls, workspace)),
            Some(&mut work),
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

    /// RKL2 molecular halves with conservative forcing and stage-time profiles.
    ///
    /// Both cumulative budgets span all phases and discarded whole splits.
    /// Native Legendre weights, bath/wall/SGS fluxes, force work, Mach and CFL
    /// checks are retained. Provider failure or exhausted work publishes no
    /// accepted state, clock, endpoint profile, load or report.
    #[allow(clippy::too_many_arguments)]
    pub fn step_stabilized_with_forcing_and_reservoir_profiles(
        &mut self,
        maximum_duration: Time,
        controls: MolecularStepControl,
        stabilization: MolecularStabilization,
        force: BodyForce3d<'_>,
        force_controls: ForceControl3d,
        workspace: &mut ForcingWorkspace3d,
        maximum_face_evaluations: usize,
        mut sample: impl FnMut(f64, BoundarySide, usize, [f64; 3]) -> Result<ReservoirFaceState, Error>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<StabilizedMolecularStepReport3d, MolecularError> {
        force_controls.admit::<Error>()?;
        let active = force.admit::<Error>(self.accepted.states.len(), &mut cancelled)?;
        workspace.admit(self.accepted.states.len())?;
        if !active {
            return self.step_stabilized_with_reservoir_profiles(
                maximum_duration,
                controls,
                stabilization,
                maximum_face_evaluations,
                sample,
                cancelled,
            );
        }
        let mut source = timed::TimedSource::new(&mut sample, maximum_face_evaluations);
        let mut work = StepWork::new(stabilization);
        let transport = self.step_impl(
            maximum_duration,
            controls,
            Some((force, force_controls, workspace)),
            Some(&mut work),
            Some(&mut source),
            cancelled,
        )?;
        Ok(StabilizedMolecularStepReport3d {
            transport,
            operator_evaluations: work.evaluations,
            maximum_stages_used: work.maximum_stages_used,
            maximum_scaled_parabolic_number: work.maximum_scaled_parabolic_number,
        })
    }

    fn molecular_half(
        &mut self,
        duration: Time,
        controls: MolecularControl,
        stabilization: Option<&mut StepWork>,
        clock: timed::Clock,
        profiles: Option<&mut timed::TimedSource<'_>>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MolecularReport3d, MolecularError> {
        match stabilization {
            None => self.molecular.advance_impl(
                &mut self.trial,
                duration,
                controls,
                clock,
                profiles,
                cancelled,
            ),
            Some(work) => {
                let report = self.molecular.advance_stabilized_with_budget(
                    &mut self.trial,
                    duration,
                    controls,
                    work.controls,
                    &mut work.evaluations,
                    clock,
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

    fn step_impl(
        &mut self,
        maximum_duration: Time,
        controls: MolecularStepControl,
        forcing: Option<(BodyForce3d<'_>, ForceControl3d, &mut ForcingWorkspace3d)>,
        stabilization: Option<&mut StepWork>,
        profiles: Option<&mut timed::TimedSource<'_>>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MolecularStepReport3d, MolecularError> {
        poll(&mut cancelled)?;
        controls.molecular.validate()?;
        if !positive(maximum_duration.value())
            || controls.maximum_coupling_retries > 32
            || !positive(controls.advection.cfl)
            || controls.advection.cfl > 0.45
            || controls.advection.maximum_halvings > 32
        {
            return Err(MolecularError::Invalid(
                "coupled molecular duration or controls",
            ));
        }
        #[cfg(feature = "experimental-sgs")]
        if let Some(mut filter) = self
            .molecular
            .sgs
            .as_mut()
            .and_then(sgs::SubgridTransport::take_lagrangian_filter)
        {
            // The borrowed history trial outlives every internal split retry.
            // Detaching only its owner leaves the held face cache usable and
            // requires no history copy or allocation. No fallible operation
            // follows acceptance of the complete flow and memory transaction.
            let result = (|| {
                let trial = self
                    .molecular
                    .sgs
                    .as_mut()
                    .expect("the Lagrangian SGS owner was just selected")
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
                .expect("a coupled step retains its original SGS owner")
                .restore_lagrangian_filter(filter, committed_fit);
            return result.map(|(report, _)| report);
        }
        self.step_prepared(
            maximum_duration,
            controls,
            forcing,
            stabilization,
            profiles,
            cancelled,
        )
    }

    fn step_prepared(
        &mut self,
        maximum_duration: Time,
        controls: MolecularStepControl,
        mut forcing: Option<(BodyForce3d<'_>, ForceControl3d, &mut ForcingWorkspace3d)>,
        mut stabilization: Option<&mut StepWork>,
        mut profiles: Option<&mut timed::TimedSource<'_>>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MolecularStepReport3d, MolecularError> {
        if let Some(work) = stabilization.as_deref() {
            work.controls.controller(controls.molecular)?;
            self.molecular
                .prepare_stabilized(work.controls.maximum_extra_workspace_bytes)?;
        }
        let before = checked_totals(
            &self.accepted.states,
            self.accepted.grid.cell_volume(),
            &mut cancelled,
        )?;
        // A failed timed split may leave endpoint data in the private owner.
        // Restore the accepted bath before every call, including a later held
        // call. Only reusable private scratch is touched until this copy passes.
        self.trial
            .reservoirs
            .profiles
            .copy_accepted_from(&self.accepted.reservoirs.profiles, &mut cancelled)?;
        self.trial.reservoirs.profiles.publish();
        let reservoirs = timed::stage(
            &mut self.trial.reservoirs,
            profiles.as_deref_mut(),
            &self.accepted.gas,
            &self.accepted.grid,
            self.accepted.time,
            true,
            &mut cancelled,
        )?
        .ok_or(MolecularError::Invalid("missing coupled reservoir view"))?;
        let mut rate = reservoirs.maximum_rate(
            &self.accepted.gas,
            &self.accepted.grid,
            &self.accepted.states,
            self.accepted.hydrostatic.as_ref(),
            &mut cancelled,
        )?;
        for q in &self.accepted.states {
            let p = q.primitive(&self.accepted.gas, &mut cancelled)?;
            if let Some((_, controls, _)) = &forcing {
                controls.admit_velocity::<Error>(p.velocity, p.sound_speed)?;
            }
            let local: f64 = (0..3)
                .map(|a| (p.velocity[a].abs() + p.sound_speed) / self.accepted.grid.spacing[a])
                .sum();
            rate = rate.max(local);
        }
        let mut dt = maximum_duration.value().min(controls.advection.cfl / rate);
        for retry in 0..=controls.maximum_coupling_retries {
            poll(&mut cancelled)?;
            if let Some(work) = stabilization.as_deref_mut() {
                work.reset_extrema();
            }
            if !positive(dt * 0.5) {
                return Err(MolecularError::Invalid("unrepresentable coupled duration"));
            }
            if profiles.is_some() {
                timed::admit_interval(self.accepted.time, dt)?;
            }
            self.trial.states.copy_from_slice(&self.accepted.states);
            self.trial.time = self.accepted.time;
            self.trial.last_step = self.accepted.last_step;
            self.trial
                .boundary_impulses
                .copy_from_slice(&self.accepted.boundary_impulses);
            let mut shortened = None;
            let attempt = (|| {
                let first = self.molecular_half(
                    Time::new(0.5 * dt),
                    controls.molecular,
                    stabilization.as_deref_mut(),
                    timed::Clock::split(self.accepted.time, dt, false),
                    profiles.as_deref_mut(),
                    &mut cancelled,
                )?;
                self.boundary_trial
                    .copy_from_slice(&self.molecular.boundary_impulses);
                let advection = match forcing.as_mut() {
                    Some((force, force_controls, workspace)) => self.trial.step_impl(
                        Time::new(dt),
                        controls.advection,
                        Some((*force, *force_controls, &mut **workspace)),
                        profiles.as_deref_mut(),
                        &mut cancelled,
                    )?,
                    None => self.trial.step_impl(
                        Time::new(dt),
                        controls.advection,
                        None,
                        profiles.as_deref_mut(),
                        &mut cancelled,
                    )?,
                };
                if advection.elapsed.value() != dt {
                    shortened = Some(advection.elapsed.value());
                    return Err(MolecularError::Invalid(
                        "Euler CFL shortened the split interval",
                    ));
                }
                for (i, (total, impulse)) in self
                    .boundary_trial
                    .iter_mut()
                    .zip(&self.trial.boundary_impulses)
                    .enumerate()
                {
                    if i % 64 == 0 {
                        poll(&mut cancelled)?;
                    }
                    *total = add_totals(*total, *impulse);
                }
                let second = self.molecular_half(
                    Time::new(0.5 * dt),
                    controls.molecular,
                    stabilization.as_deref_mut(),
                    timed::Clock::split(self.accepted.time, dt, true),
                    profiles.as_deref_mut(),
                    &mut cancelled,
                )?;
                for (i, (total, impulse)) in self
                    .boundary_trial
                    .iter_mut()
                    .zip(&self.molecular.boundary_impulses)
                    .enumerate()
                {
                    if i % 64 == 0 {
                        poll(&mut cancelled)?;
                    }
                    *total = add_totals(*total, *impulse);
                }
                if let Some((_, force_controls, _)) = &forcing {
                    let mut endpoint_rate = 0.0_f64;
                    for state in &self.trial.states {
                        let primitive = state.primitive(&self.trial.gas, &mut cancelled)?;
                        force_controls
                            .admit_velocity::<Error>(primitive.velocity, primitive.sound_speed)?;
                        let local: f64 = (0..3)
                            .map(|axis| {
                                (primitive.velocity[axis].abs() + primitive.sound_speed)
                                    / self.trial.grid.spacing[axis]
                            })
                            .sum();
                        endpoint_rate = endpoint_rate.max(local);
                    }
                    check_cfl::<Error>(endpoint_rate * dt, controls.advection.cfl)?;
                }
                let after = checked_totals(
                    &self.trial.states,
                    self.accepted.grid.cell_volume(),
                    &mut cancelled,
                )?;
                let outward = sum_arrays(
                    self.boundary_trial.iter().map(|q| q.array()),
                    &mut cancelled,
                )?;
                let boundary_change = Totals::from_array(outward.map(|v| -v));
                let source_change = advection.source_change;
                let source_totals = Totals {
                    mass: 0.0,
                    momentum: source_change.impulse,
                    energy: source_change.work,
                };
                let conservation_residual =
                    balance(before, after, add_totals(boundary_change, source_totals))?;
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
                })
            })();
            match attempt {
                Ok(report) => {
                    std::mem::swap(&mut self.accepted, &mut self.trial);
                    self.boundary_impulses.copy_from_slice(&self.boundary_trial);
                    self.last_step = Some(report);
                    return Ok(report);
                }
                Err(error)
                    if error.is_terminal()
                        || profiles.as_deref().is_some_and(timed::TimedSource::failed) =>
                {
                    return Err(error);
                }
                Err(error) => {
                    if retry == controls.maximum_coupling_retries {
                        return Err(MolecularError::CouplingRetries {
                            attempts: retry + 1,
                            last: Box::new(error),
                        });
                    }
                    dt = shortened.unwrap_or(dt * 0.5);
                }
            }
        }
        unreachable!("bounded coupled retries return on the final attempt")
    }
}
