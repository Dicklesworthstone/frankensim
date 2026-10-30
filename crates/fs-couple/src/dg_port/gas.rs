//! Temporal source-gas changes on the original, continuing DG port system.

use super::*;
use fs_duct::dg::DgGasMaterialUpdate;
use fs_material::gas_field::AxialGasField;

/// Native gas-history transfer selected by the existing acoustic owner.
/// None of these choices replaces or discards an unsupported history.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum DgGasRefresh {
    /// Uniform gas without viscothermal gas memory. Wall memory is retained.
    Lossless,
    /// Uniform circular diffusion with its original physical arm histories.
    UniformDiffusion {
        /// Maximum relative error for the new gas on the retained angular band.
        maximum_relative_error: f64,
    },
    /// Spatial circular diffusion on the original weak basis and quadrature.
    SpatialDiffusion {
        /// Maximum relative error for every local gas on the retained band.
        maximum_relative_error: f64,
    },
    /// Uniform molecular/thermal gas with original radial/spatial histories.
    /// Retains the original species fractions and source molecular law;
    /// preserves entropy perturbations, compression and viscous momentum.
    UniformMolecularDiffusion {
        /// New-gas complex and loss admission on the original angular band.
        maximum_relative_error: f64,
    },
    /// Nonuniform molecular/thermal gas on the original weak basis and wall
    /// quadrature. Preserves local reservoir entropy, weak compression and
    /// momentum at unchanged primitive stations and species fractions.
    SpatialMolecularDiffusion {
        /// Every local gas must pass the original complex and loss band.
        maximum_relative_error: f64,
    },
}

/// One accepted source-gas intervention in the complete represented load.
/// Rational exterior storage is not inferred by this account.
#[derive(Clone, Debug, PartialEq)]
pub struct DgPortGasUpdate {
    /// Original native main-bore transfer, including wall/diffusion memory.
    pub duct: DgGasMaterialUpdate,
    /// Native channel transfers in original attachment order. Independently
    /// sourced branch channels retain their gas and have zero parameter work.
    pub terminal_channels: Vec<DgTerminalChannelGasUpdate>,
    /// Bore, finite-branch and finite-channel storage before the change [J].
    pub energy_before_j: f64,
    /// The same represented storage after the change [J].
    pub energy_after_j: f64,
    /// Signed sum of the native bore and channel storage changes [J].
    pub parameter_work_j: f64,
    /// Sum of absolute native interventions charged to the work allowance [J].
    /// Opposite signed transfers cannot conceal an excessive intervention.
    pub absolute_parameter_work_j: f64,
}

/// A source-field change on an existing finite branch. The external port ID
/// identifies the original branch; slice order never reassigns its gas.
#[derive(Clone, Copy, Debug)]
pub struct DgBranchGasInput<'a> {
    pub port: PortId,
    pub field: &'a AxialGasField,
    pub refresh: DgGasRefresh,
}

/// Native transfer on one explicitly selected branch, in attachment order.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DgBranchGasUpdate {
    pub port: PortId,
    pub duct: DgGasMaterialUpdate,
}

/// One atomic material intervention on selected regions of the original load.
/// Unselected ducts keep their original propagators and physical histories.
#[derive(Clone, Debug, PartialEq)]
pub struct DgRegionalGasUpdate {
    pub main: Option<DgGasMaterialUpdate>,
    pub branches: Vec<DgBranchGasUpdate>,
    /// Every original terminal channel, in attachment order. Its local outlet
    /// gas selects the native momentum transfer; unchanged sources report zero.
    pub terminal_channels: Vec<DgTerminalChannelGasUpdate>,
    pub energy_before_j: f64,
    pub energy_after_j: f64,
    /// Signed native main, branch and channel work INTO represented motion [J].
    pub parameter_work_j: f64,
    /// Sum of the absolute native contributions, bounded before publication [J].
    pub absolute_parameter_work_j: f64,
}

fn duct_error(error: fs_duct::DuctError) -> AcousticRealizeError {
    match error {
        fs_duct::DuctError::Cancelled { .. } => AcousticRealizeError::Cancelled,
        error => AcousticRealizeError::Duct(error),
    }
}

fn field_error(error: fs_duct::profile::gas_field::GasFieldDuctError) -> AcousticRealizeError {
    match error {
        fs_duct::profile::gas_field::GasFieldDuctError::Duct(error) => duct_error(error),
        error => AcousticRealizeError::Nonlinear(error.to_string()),
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), AcousticRealizeError> {
    if cancelled() {
        Err(AcousticRealizeError::Cancelled)
    } else {
        Ok(())
    }
}

fn charge_work(
    signed: &mut f64,
    absolute: &mut f64,
    work: f64,
    allowance: f64,
) -> Result<(), AcousticRealizeError> {
    *signed += work;
    *absolute += work.abs();
    if !signed.is_finite() || !absolute.is_finite() || *absolute > allowance {
        return Err(invalid(
            "live regional gas changes exceed the total absolute parameter-work allowance",
        ));
    }
    Ok(())
}

/// Shared preparation for the source wall and nonlinear port adapters. The
/// original native transfer alone owns physical memory and parameter work.
pub(crate) fn prepare_dg_gas(
    duct: &DgDuct,
    state: &DgState,
    propagator: &DgPropagator,
    field: &AxialGasField,
    refresh: DgGasRefresh,
    maximum_absolute_parameter_work_j: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(DgDuct, DgState, DgPropagator, DgGasMaterialUpdate), AcousticRealizeError> {
    poll(cancelled)?;
    let mut state = state.clone();
    propagator.sync_nodal(&mut state).map_err(duct_error)?;
    let mut candidate = duct.clone();
    let report = match refresh {
        DgGasRefresh::Lossless => candidate
            .rebind_uniform_gas_field(
                &mut state,
                field,
                maximum_absolute_parameter_work_j,
                cancelled,
            )
            .map_err(duct_error),
        DgGasRefresh::UniformDiffusion {
            maximum_relative_error,
        } => candidate
            .rebind_uniform_viscothermal_gas_field(
                &mut state,
                field,
                maximum_relative_error,
                maximum_absolute_parameter_work_j,
                cancelled,
            )
            .map_err(duct_error),
        DgGasRefresh::SpatialDiffusion {
            maximum_relative_error,
        } => candidate
            .rebind_spatial_viscothermal_gas_field(
                &mut state,
                field,
                maximum_relative_error,
                maximum_absolute_parameter_work_j,
                cancelled,
            )
            .map_err(field_error),
        DgGasRefresh::UniformMolecularDiffusion {
            maximum_relative_error,
        } => candidate
            .rebind_uniform_molecular_gas_field(
                &mut state,
                field,
                maximum_relative_error,
                maximum_absolute_parameter_work_j,
                cancelled,
            )
            .map_err(field_error),
        DgGasRefresh::SpatialMolecularDiffusion {
            maximum_relative_error,
        } => candidate
            .rebind_spatial_molecular_gas_field(
                &mut state,
                field,
                maximum_relative_error,
                maximum_absolute_parameter_work_j,
                cancelled,
            )
            .map_err(field_error),
    }?;
    if candidate.total_dofs() != duct.total_dofs() {
        return Err(invalid(
            "live gas transfer cannot change the admitted native state count",
        ));
    }
    poll(cancelled)?;
    // Native identity covers the operator, physical mass transforms, geometry
    // and typed boundaries. An unchanged operator keeps its already-certified
    // exponential at the same clock/tolerance after full source/transfer admission.
    // This avoids another matrix exponential for unchanged finite-gas trials.
    let stepping = if candidate.operator_identity() == duct.operator_identity() {
        propagator.clone()
    } else {
        candidate
            .propagator(propagator.dt(), propagator.tolerance())
            .map_err(duct_error)?
    };
    poll(cancelled)?;
    Ok((candidate, state, stepping, report))
}

impl DgPortLoad {
    /// Refresh the MAIN bore gas on this continuing acoustic trajectory.
    ///
    /// The original native transfer preserves weak compression and momentum,
    /// all admitted viscothermal arms, wall motion and hereditary histories.
    /// Pressure, flow, matched boundaries and the exact mean source impedance
    /// follow the new complete source field. The original clock, dissipative
    /// accounts and last accepted frames survive; this consumes no sample.
    ///
    /// Original finite side branches retain their independently supplied gas,
    /// geometry and physical histories. Their junction elimination and the
    /// terminal Schur complement are rebuilt for the new main bore. The causal
    /// exterior load, including a driven mouth, keeps its source and EVERY
    /// pole state: changing the interior does not redefine the exterior gas.
    /// Branch junctions and computed exterior openings must still share the
    /// original equilibrium pressure; an unsupported pressure jump refuses.
    /// Attached slit/circular channels follow their actual local gas. Their
    /// original native owners retain every transverse/radial cell momentum,
    /// refresh viscous response and report kinetic parameter work separately.
    /// Independently supplied branch gaps retain their gas and physical state.
    /// Rigid closures retain their source without inventing fluid storage.
    ///
    /// Actual main-bore and channel storage changes are separately reported and
    /// charged once as parameter work. The sum of absolute interventions must
    /// fit the supplied work allowance. The background owner owes signed work
    /// and any heat/species balance. This does not solve transport or add mean
    /// pressure traction. Every source, state, propagator and algebraic coupling
    /// publishes together after all checks; failure preserves the entire load.
    pub fn rebind_gas_field_controlled(
        &mut self,
        field: &AxialGasField,
        refresh: DgGasRefresh,
        maximum_absolute_parameter_work_j: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DgPortGasUpdate, AcousticRealizeError> {
        let report = self.rebind_regional_gas_fields_controlled(
            Some((field, refresh)),
            &[],
            maximum_absolute_parameter_work_j,
            cancelled,
        )?;
        Ok(DgPortGasUpdate {
            duct: report
                .main
                .expect("the main source was explicitly selected"),
            terminal_channels: report.terminal_channels,
            energy_before_j: report.energy_before_j,
            energy_after_j: report.energy_after_j,
            parameter_work_j: report.parameter_work_j,
            absolute_parameter_work_j: report.absolute_parameter_work_j,
        })
    }

    /// Change the main source and/or named finite branches in ONE continuing load.
    ///
    /// `main = None` retains the original main operator and propagator. Branch
    /// inputs identify existing external port IDs, not positions; unknown or
    /// duplicate IDs refuse before native work. Unlisted branches retain their
    /// gas and their original propagation. Each changed duct uses its own
    /// admitted native uniform/spatial diffusion and wall-memory transfer.
    ///
    /// Junction and computed-exterior equilibrium pressures are checked on the
    /// complete proposed endpoint. This admits coordinated pressure changes
    /// that cannot pass through separately published intermediate regions.
    /// Every finite terminal passage follows its actual NEW local outlet gas,
    /// preserving cell momentum and retaining its original Mach admission.
    /// Branch and terminal elimination are rebuilt once from the final sources.
    /// Exterior pole histories, clocks and accepted interval frames survive.
    ///
    /// The absolute-work allowance covers all changed ducts and channels before
    /// any owner is published. Signed work remains an obligation of the gas or
    /// thermal owner; this prescribed fixed-geometry intervention does not
    /// itself supply heat, transport species or add static pressure traction.
    pub fn rebind_regional_gas_fields_controlled(
        &mut self,
        main: Option<(&AxialGasField, DgGasRefresh)>,
        branch_inputs: &[DgBranchGasInput<'_>],
        maximum_absolute_parameter_work_j: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DgRegionalGasUpdate, AcousticRealizeError> {
        poll(cancelled)?;
        if !maximum_absolute_parameter_work_j.is_finite()
            || maximum_absolute_parameter_work_j < 0.0
            || (main.is_none() && branch_inputs.is_empty())
            || branch_inputs.len() > self.branches.len()
        {
            return Err(invalid(
                "regional gas refresh requires selected original sources and a finite nonnegative work allowance",
            ));
        }
        for (index, input) in branch_inputs.iter().enumerate() {
            poll(cancelled)?;
            if !self
                .branches
                .iter()
                .any(|branch| branch.port() == input.port)
                || branch_inputs[..index]
                    .iter()
                    .any(|old| old.port == input.port)
            {
                return Err(invalid(
                    "regional gas refresh requires unique existing branch port IDs",
                ));
            }
        }
        let energy_before_j = self.stored_energy_j();
        let main_candidate = main
            .map(|(field, refresh)| {
                prepare_dg_gas(
                    &self.duct,
                    &self.state,
                    &self.propagator,
                    field,
                    refresh,
                    maximum_absolute_parameter_work_j,
                    cancelled,
                )
            })
            .transpose()?;
        let (duct, state, propagator) = main_candidate.as_ref().map_or(
            (&self.duct, &self.state, &self.propagator),
            |(duct, state, propagator, _)| (duct, state, propagator),
        );
        let main_report = main_candidate.as_ref().map(|(_, _, _, report)| *report);
        let mut parameter_work_j = 0.0;
        let mut absolute_parameter_work_j = 0.0;
        if let Some(report) = main_report {
            charge_work(
                &mut parameter_work_j,
                &mut absolute_parameter_work_j,
                report.parameter_work_j,
                maximum_absolute_parameter_work_j,
            )?;
        }
        let mut branch_reports = Vec::with_capacity(branch_inputs.len());
        // The common main-only path borrows every original branch. A regional
        // change allocates one candidate set and never rebuilds unselected
        // native propagators or refits their material histories.
        let branch_candidates = if branch_inputs.is_empty() {
            None
        } else {
            let mut candidates = Vec::with_capacity(self.branches.len());
            for branch in &self.branches {
                poll(cancelled)?;
                if let Some(input) = branch_inputs
                    .iter()
                    .find(|input| input.port == branch.port())
                {
                    let (candidate, report) = branch.prepare_gas(
                        input.field,
                        input.refresh,
                        maximum_absolute_parameter_work_j,
                        cancelled,
                    )?;
                    charge_work(
                        &mut parameter_work_j,
                        &mut absolute_parameter_work_j,
                        report.parameter_work_j,
                        maximum_absolute_parameter_work_j,
                    )?;
                    branch_reports.push(DgBranchGasUpdate {
                        port: branch.port(),
                        duct: report,
                    });
                    candidates.push(candidate);
                } else {
                    candidates.push(branch.clone());
                }
            }
            Some(candidates)
        };
        let branches = branch_candidates.as_deref().unwrap_or(&self.branches);
        for branch in branches {
            poll(cancelled)?;
            let source = duct.side_port_source(branch.port()).map_err(field_error)?;
            if source.gas.pressure.to_bits() != branch.duct().inlet_gas().pressure.to_bits() {
                return Err(invalid(
                    "live DG gas update requires equal equilibrium pressure at each retained branch junction",
                ));
            }
        }
        let terminal_gas = |port| -> Result<_, AcousticRealizeError> {
            if let Some(branch) = branches.iter().find(|branch| branch.port() == port) {
                Ok(*branch.duct().outlet_gas())
            } else if duct.outlet() == DgBoundary::Port(port) {
                Ok(*duct.outlet_gas())
            } else {
                Ok(duct.side_port_source(port).map_err(field_error)?.gas)
            }
        };
        if let Some(record) = &self.radiation_record {
            let pressure = record.gas_record().state().pressure.to_bits();
            for &port in self.terminal_ports() {
                poll(cancelled)?;
                if terminal_gas(port)?.pressure.to_bits() != pressure {
                    return Err(invalid(
                        "live DG openings and their retained computed exterior require equal equilibrium pressures",
                    ));
                }
            }
            if self.has_driven_mouth() && duct.inlet_gas().pressure.to_bits() != pressure {
                return Err(invalid(
                    "live DG driven mouth and its retained exterior require equal equilibrium pressures",
                ));
            }
        }
        let mut terminal_channels = Vec::with_capacity(self.terminal_channels.len());
        let mut channel_reports = Vec::with_capacity(self.terminal_channels.len());
        for channel in &self.terminal_channels {
            poll(cancelled)?;
            let gas = terminal_gas(channel.input().port)?;
            let (candidate, report) =
                channel.prepare_gas(&gas, maximum_absolute_parameter_work_j, cancelled)?;
            let work = report
                .transfer
                .map_or(0.0, |transfer| transfer.parameter_work_j);
            charge_work(
                &mut parameter_work_j,
                &mut absolute_parameter_work_j,
                work,
                maximum_absolute_parameter_work_j,
            )?;
            terminal_channels.push(candidate);
            channel_reports.push(report);
        }
        let response = propagator.multiport_response(state).map_err(duct_error)?;
        let branch_coupling = if branches.is_empty() {
            None
        } else {
            Some(BranchCoupling::new(&response, branches)?)
        };
        let response = match &branch_coupling {
            Some(coupling) => coupling.response(response, branches)?,
            None => response,
        };
        let terminal_coupling = match (&self.terminal, &self.terminal_coupling) {
            (Some(terminal), Some(previous)) => Some(TerminalCoupling::with_channels(
                &response,
                self.inlet,
                previous.ports.clone(),
                terminal,
                previous.driven_mouth,
                &terminal_channels,
            )?),
            (None, None) => None,
            _ => {
                return Err(invalid(
                    "live DG gas update lost the original terminal owner",
                ));
            }
        };
        let impedance = match (&self.terminal, &terminal_coupling) {
            (Some(terminal), Some(coupling)) => {
                coupling
                    .relation(response, terminal, &terminal_channels)?
                    .inlet
                    .impedance_pa_s_m3
            }
            _ if branch_coupling.is_none() => propagator.port_response(state, self.inlet)
                .map_err(duct_error)?.impedance_pa_s_m3,
            _ => {
                if response.ports != [self.inlet] {
                    return Err(invalid("live DG inlet-only response lost sealed branch ownership"));
                }
                response.impedance_pa_s_m3[0]
            }
        };
        if !impedance.is_finite() || impedance <= 0.0 {
            return Err(invalid(
                "live DG gas update requires strictly positive exact mean-port impedance",
            ));
        }
        let duct_energy = branches
            .iter()
            .fold(state.stored_energy_j(), |sum, branch| {
                sum + branch.state().stored_energy_j()
            });
        let energy_after_j = terminal_channels
            .iter()
            .fold(duct_energy, |sum, channel| sum + channel.stored_energy_j());
        if !energy_before_j.is_finite() || !energy_after_j.is_finite() {
            return Err(invalid("live DG/channel gas storage cannot be represented"));
        }
        let report = DgRegionalGasUpdate {
            main: main_report,
            branches: branch_reports,
            terminal_channels: channel_reports,
            energy_before_j,
            energy_after_j,
            parameter_work_j,
            absolute_parameter_work_j,
        };
        poll(cancelled)?;
        if let Some((duct, state, propagator, _)) = main_candidate {
            self.duct = duct;
            self.state = state;
            self.propagator = propagator;
        }
        if let Some(branches) = branch_candidates {
            self.branches = branches;
        }
        self.branch_coupling = branch_coupling;
        self.terminal_coupling = terminal_coupling;
        self.terminal_channels = terminal_channels;
        self.impedance_pa_s_m3 = impedance;
        Ok(report)
    }

    /// Signed parameter work accumulated in all represented duct/channel
    /// owners [J]. Pressure work, viscous heat and numerical loss stay separate.
    #[must_use]
    pub fn cumulative_parameter_work_j(&self) -> f64 {
        let ducts = self
            .branches
            .iter()
            .fold(self.state.cumulative_parameter_work_j(), |sum, branch| {
                sum + branch.state().cumulative_parameter_work_j()
            });
        self.terminal_channels.iter().fold(ducts, |sum, channel| {
            sum + channel.cumulative_parameter_work_j()
        })
    }
}

#[cfg(test)]
mod tests;
