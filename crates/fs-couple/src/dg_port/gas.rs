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
    }?;
    if candidate.total_dofs() != duct.total_dofs() {
        return Err(invalid(
            "live gas transfer cannot change the admitted native state count",
        ));
    }
    poll(cancelled)?;
    let stepping = candidate
        .propagator(propagator.dt(), propagator.tolerance())
        .map_err(duct_error)?;
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
    /// An attached finite gap whose local gas would change refuses, because
    /// its transverse momentum transfer is not represented by this operation.
    /// Independently supplied branch gaps whose gas stays equal remain valid.
    ///
    /// Actual main-bore storage change is separately reported and charged once
    /// as parameter work. The owner prescribing the background owes that work
    /// and any heat/species balance. This does not solve transport or add mean
    /// pressure traction. Every source, state, propagator and algebraic coupling
    /// publishes together after all checks; failure preserves the entire load.
    pub fn rebind_gas_field_controlled(
        &mut self,
        field: &AxialGasField,
        refresh: DgGasRefresh,
        maximum_absolute_parameter_work_j: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DgGasMaterialUpdate, AcousticRealizeError> {
        let (duct, state, propagator, report) = prepare_dg_gas(
            &self.duct,
            &self.state,
            &self.propagator,
            field,
            refresh,
            maximum_absolute_parameter_work_j,
            cancelled,
        )?;
        for branch in &self.branches {
            poll(cancelled)?;
            let source = duct.side_port_source(branch.port()).map_err(field_error)?;
            if source.gas.pressure.to_bits() != branch.duct().inlet_gas().pressure.to_bits() {
                return Err(invalid(
                    "live DG gas update requires equal equilibrium pressure at each retained branch junction",
                ));
            }
        }
        let terminal_gas = |port| -> Result<_, AcousticRealizeError> {
            if let Some(branch) = self.branches.iter().find(|branch| branch.port() == port) {
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
        for channel in &self.terminal_channels {
            poll(cancelled)?;
            let gas = terminal_gas(channel.input().port)?;
            if gas != *channel.gas() {
                return Err(invalid(
                    "live bore gas change requires physical transfer of the attached terminal channel gas and momentum",
                ));
            }
        }
        let response = propagator.multiport_response(&state).map_err(duct_error)?;
        let branch_coupling = if self.branches.is_empty() {
            None
        } else {
            Some(BranchCoupling::new(&response, &self.branches)?)
        };
        let response = match &branch_coupling {
            Some(coupling) => coupling.response(response, &self.branches)?,
            None => response,
        };
        let terminal_coupling = match (&self.terminal, &self.terminal_coupling) {
            (Some(terminal), Some(previous)) => Some(TerminalCoupling::with_channels(
                &response,
                self.inlet,
                previous.ports.clone(),
                terminal,
                previous.driven_mouth,
                &self.terminal_channels,
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
                    .relation(response, terminal, &self.terminal_channels)?
                    .inlet
                    .impedance_pa_s_m3
            }
            _ => {
                propagator
                    .port_response(&state, self.inlet)
                    .map_err(duct_error)?
                    .impedance_pa_s_m3
            }
        };
        if !impedance.is_finite() || impedance <= 0.0 {
            return Err(invalid(
                "live DG gas update requires strictly positive exact mean-port impedance",
            ));
        }
        poll(cancelled)?;
        self.duct = duct;
        self.state = state;
        self.propagator = propagator;
        self.branch_coupling = branch_coupling;
        self.terminal_coupling = terminal_coupling;
        self.impedance_pa_s_m3 = impedance;
        Ok(report)
    }
}

#[cfg(test)]
mod tests;
