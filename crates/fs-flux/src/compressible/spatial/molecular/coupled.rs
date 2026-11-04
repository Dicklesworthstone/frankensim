use super::*;

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

/// Atomic source-driven three-dimensional compressible viscous/thermal flow.
///
/// A symmetric molecular/advection/molecular composition retains one actual
/// accepted duration. A shorter Euler step or a refused later stage restarts
/// the complete trial from the accepted field. Both molecular stages use SSP
/// RK2; smooth MUSCL Euler permits second-order splitting, while limiting or
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
        trial.reservoirs = flow.reservoirs.clone();
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

    /// Exact immutable source law, including retained material/collision data.
    pub fn law(&self) -> &L {
        self.molecular.law()
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
        self.step_impl(maximum_duration, controls, None, cancelled)
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
            maximum_duration, controls, Some((force, force_controls, workspace)), cancelled,
        )
    }

    fn step_impl(
        &mut self,
        maximum_duration: Time,
        controls: MolecularStepControl,
        mut forcing: Option<(BodyForce3d<'_>, ForceControl3d, &mut ForcingWorkspace3d)>,
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
        let before = checked_totals(
            &self.accepted.states,
            self.accepted.grid.cell_volume(),
            &mut cancelled,
        )?;
        let mut rate: f64 = 0.0;
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
            if !positive(dt * 0.5) {
                return Err(MolecularError::Invalid("unrepresentable coupled duration"));
            }
            self.trial.states.copy_from_slice(&self.accepted.states);
            self.trial.time = self.accepted.time;
            self.trial.last_step = self.accepted.last_step;
            self.trial
                .boundary_impulses
                .copy_from_slice(&self.accepted.boundary_impulses);
            let mut shortened = None;
            let attempt = (|| {
                let first = self.molecular.advance(
                    &mut self.trial,
                    Time::new(0.5 * dt),
                    controls.molecular,
                    &mut cancelled,
                )?;
                self.boundary_trial
                    .copy_from_slice(&self.molecular.boundary_impulses);
                let advection = match forcing.as_mut() {
                    Some((force, force_controls, workspace)) => self.trial.step_impl(
                        Time::new(dt), controls.advection,
                        Some((*force, *force_controls, &mut **workspace)), &mut cancelled,
                    )?,
                    None => self.trial.step(
                        Time::new(dt), controls.advection, &mut cancelled,
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
                let second = self.molecular.advance(
                    &mut self.trial,
                    Time::new(0.5 * dt),
                    controls.molecular,
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
                        force_controls.admit_velocity::<Error>(
                            primitive.velocity, primitive.sound_speed,
                        )?;
                        let local: f64 = (0..3).map(|axis| {
                            (primitive.velocity[axis].abs() + primitive.sound_speed)
                                / self.trial.grid.spacing[axis]
                        }).sum();
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
                    mass: 0.0, momentum: source_change.impulse, energy: source_change.work,
                };
                let conservation_residual = balance(
                    before, after, add_totals(boundary_change, source_totals),
                )?;
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
                Err(MolecularError::Cancelled) => return Err(MolecularError::Cancelled),
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
