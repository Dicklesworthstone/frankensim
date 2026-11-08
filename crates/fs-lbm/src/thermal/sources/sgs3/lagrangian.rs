//! Pathline-averaged momentum and heat in the original D3Q19/D3Q7 slab.
//!
//! The shared Favre/Germano owner carries independent signed stress and heat
//! histories along one accepted velocity trajectory. The native source owner
//! still supplies forcing, thermal walls, sensible storage and publication.
//! A prepared token keeps both histories pending through every outer material
//! query; only a successful physical transaction may commit it.

use super::*;
use crate::sgs3::dynamic::lagrangian::scalar::{
    LagrangianDynamicScalar, LagrangianScalarConfig, LagrangianScalarTrial,
};
use crate::sgs3::dynamic::lagrangian::{AdaptiveLagrangianMemory, LagrangianSgsConfig};
use crate::sgs3::dynamic::{DynamicSgsError, DynamicSgsReport, scalar::DynamicScalarReport};

/// Independent stress/heat memories; every time is in lattice timesteps.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermalLagrangianSgsConfig {
    pub momentum: LagrangianSgsConfig,
    pub heat: LagrangianScalarConfig,
    /// Admitted before preparation and after the actual flow/source step.
    pub maximum_mach: f64,
}

impl ThermalLagrangianSgsConfig {
    pub fn validate(self) -> Result<(), ThermalSgsError> {
        self.momentum.validate()?;
        self.heat.validate()?;
        if !self.maximum_mach.is_finite() || self.maximum_mach <= 0.0 || self.maximum_mach > 0.3 {
            return Err(ThermalSgsError::Invalid(
                "SGS Mach bound must lie in (0, 0.3]",
            ));
        }
        Ok(())
    }
}

/// Reused fluid-only sample arrays, full-slab rate arrays and one joint history.
/// This owns a single physical trajectory. Do not reuse it with a different
/// flow state merely because its grid dimensions happen to match.
pub struct ThermalLagrangianSgsWorkspace3 {
    velocity: Vec<[f64; 3]>,
    density: Vec<f64>,
    temperature: Vec<f64>,
    momentum: Vec<f64>,
    thermal: Vec<f64>,
    filter: LagrangianDynamicScalar,
}

impl ThermalLagrangianSgsWorkspace3 {
    /// Select bounded invariant momentum memory before the first accepted fit.
    /// Heat keeps its separately authored scalar memory time.
    pub fn with_adaptive_momentum_memory(
        mut self,
        config: AdaptiveLagrangianMemory,
    ) -> Result<Self, ThermalSgsError> {
        self.filter = self.filter.with_adaptive_momentum_memory(config)?;
        Ok(self)
    }

    pub const fn accepted_time(&self) -> Option<f64> {
        self.filter.accepted_time()
    }

    pub const fn accepted_momentum_memory_time_range(&self) -> Option<[f64; 2]> {
        self.filter.accepted_momentum_memory_time_range()
    }

    /// Most recently prepared rate scratch, valid only after successful
    /// preparation. It is never a molecular input to a later step.
    pub fn relaxation(&self) -> RelaxationField<'_> {
        RelaxationField {
            momentum: &self.momentum,
            thermal: &self.thermal,
        }
    }

    /// Actual numerical capacities; metadata and the native source snapshot
    /// are separate. The original source max_cells still admits its snapshot.
    pub fn workspace_bytes(&self) -> usize {
        (3 * self.velocity.capacity()
            + self.density.capacity()
            + self.temperature.capacity()
            + self.momentum.capacity()
            + self.thermal.capacity())
            * core::mem::size_of::<f64>()
            + self.filter.workspace_bytes()
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermalLagrangianSgsReport {
    pub momentum_fit: DynamicSgsReport,
    pub heat_fit: DynamicScalarReport,
    pub transport: ThermalSgsReport,
    /// Actual min/max time used by this candidate's momentum update.
    pub momentum_memory_time_range: [f64; 2],
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermalLagrangianSgsStep {
    pub heat: ThermalSourceBalance,
    pub sgs: ThermalLagrangianSgsReport,
}

/// Pending histories and their admitted total rates. The lifetime binds only
/// the workspace, not the sampled state, molecular input or source arrays.
/// Dropping a token rejects both histories. Private scratch can still change.
#[must_use = "commit only after the complete consuming physical step succeeds"]
pub struct ThermalLagrangianSgsTrial<'a> {
    history: LagrangianScalarTrial<'a>,
    momentum: &'a [f64],
    thermal: &'a [f64],
    report: ThermalLagrangianSgsReport,
}

impl ThermalLagrangianSgsTrial<'_> {
    pub const fn report(&self) -> ThermalLagrangianSgsReport {
        self.report
    }

    pub fn relaxation(&self) -> RelaxationField<'_> {
        RelaxationField {
            momentum: self.momentum,
            thermal: self.thermal,
        }
    }

    /// Infallible publication after all native and outer material checks.
    pub fn commit(self) {
        self.history.commit();
    }
}

impl ThermalLbm3 {
    pub fn required_thermal_lagrangian_sgs_workspace_bytes(
        &self,
    ) -> Result<usize, ThermalSgsError> {
        let (cells, plane) = self.layout()?;
        let fluid = cells - 2 * plane;
        let history = LagrangianDynamicScalar::required_bytes(fluid)?;
        fluid
            .checked_mul(5 * core::mem::size_of::<f64>())
            .and_then(|a| {
                cells
                    .checked_mul(2 * core::mem::size_of::<f64>())
                    .and_then(|b| a.checked_add(b))
            })
            .and_then(|a| a.checked_add(history))
            .ok_or(ThermalSgsError::Invalid(
                "thermal Lagrangian workspace size overflow",
            ))
    }

    /// Prepare a resolved channel: periodic x/y and the original rigid z walls.
    /// Fixed-temperature walls require two fluid layers; insulating or interior
    /// reconstructed wall derivatives require three when fitting heat.
    pub fn thermal_lagrangian_sgs_workspace(
        &self,
        maximum_bytes: usize,
    ) -> Result<ThermalLagrangianSgsWorkspace3, ThermalSgsError> {
        let (cells, plane) = self.layout()?;
        let fluid = cells - 2 * plane;
        let required = self.required_thermal_lagrangian_sgs_workspace_bytes()?;
        if required > maximum_bytes {
            return Err(DynamicSgsError::WorkspaceBytes {
                required,
                maximum: maximum_bytes,
            }
            .into());
        }
        let result = ThermalLagrangianSgsWorkspace3 {
            filter: LagrangianDynamicScalar::new_channel(
                [self.grid.nx, self.grid.ny, self.grid.nz - 2],
                [1.0; 3],
                2,
                maximum_bytes,
            )?,
            velocity: filled(fluid, [0.0; 3])?,
            density: filled(fluid, 0.0)?,
            temperature: filled(fluid, 0.0)?,
            momentum: filled(cells, 0.0)?,
            thermal: filled(cells, 0.0)?,
        };
        let required = result.workspace_bytes();
        if required > maximum_bytes {
            return Err(DynamicSgsError::WorkspaceBytes {
                required,
                maximum: maximum_bytes,
            }
            .into());
        }
        Ok(result)
    }

    /// Fit independent stress and heat from one accepted Guo-centered rho/u/T
    /// sample, retaining both histories until the consuming transaction accepts.
    /// sample_time_lu is the actual current time in lattice timesteps. This slab
    /// has no authoritative clock; the caller owns that time and the trajectory.
    /// It must increase between accepted fits, including any skipped intervals.
    /// Memory times and departure bounds use that same elapsed lattice time.
    ///
    /// The first sample uses the shared local tangential average; subsequent
    /// samples advect signed contractions without crossing either wall. Current
    /// scalar wall values/derivative policies are identical to the instantaneous
    /// channel fit. Momentum and heat each enforce their own molecular ratio.
    /// Frozen rates apply for one native step. Modeled mechanical dissipation
    /// is diagnostic; this remains low-Mach sensible heat, not total energy.
    pub fn prepare_thermal_lagrangian_sgs<'a>(
        &self,
        config: ThermalLagrangianSgsConfig,
        sample_time_lu: f64,
        molecular: Option<RelaxationField<'_>>,
        workspace: &'a mut ThermalLagrangianSgsWorkspace3,
        gate: &CancelGate,
    ) -> Result<ThermalLagrangianSgsTrial<'a>, ThermalSgsError> {
        poll(gate)?;
        config.validate()?;
        let (cells, plane) = self.layout()?;
        let fluid = cells - 2 * plane;
        let walls = channel::scalar_walls(self, plane)?;
        if workspace.filter.shape() != [self.grid.nx, self.grid.ny, self.grid.nz - 2]
            || workspace.filter.normal_axis() != Some(2)
            || workspace.velocity.len() != fluid
            || workspace.density.len() != fluid
            || workspace.temperature.len() != fluid
            || workspace.momentum.len() != cells
            || workspace.thermal.len() != cells
        {
            return Err(ThermalSgsError::Invalid(
                "Lagrangian SGS workspace does not match the slab",
            ));
        }
        if !self.flow_valid() {
            return Err(ThermalSgsError::Invalid("invalid thermal/flow state"));
        }
        if let Some(field) = molecular {
            field.admit(cells)?;
        }
        for cell in 0..cells {
            if cell % 64 == 0 {
                poll(gate)?;
            }
            workspace.momentum[cell] = molecular.map_or(self.grid.tau[cell], |f| f.momentum[cell]);
            workspace.thermal[cell] = molecular.map_or(self.tau_g, |f| f.thermal[cell]);
            if cell >= plane && cell < cells - plane {
                let active = cell - plane;
                let (rho, velocity) = self.grid.moments3(cell);
                admit_velocity(rho, velocity, cell, config.maximum_mach)?;
                workspace.density[active] = rho;
                workspace.velocity[active] = velocity;
                workspace.temperature[active] = self.gpop[cell].iter().sum();
            }
        }
        let history = workspace.filter.prepare_with_walls(
            config.momentum,
            config.heat,
            sample_time_lu,
            &workspace.density,
            &workspace.velocity,
            [[0.0; 3]; 2],
            &workspace.temperature,
            walls,
            || gate.is_requested(),
        )?;
        let transport = channel::augment_transport(
            self,
            &workspace.density,
            &workspace.velocity,
            &mut workspace.momentum,
            &mut workspace.thermal,
            history.viscosities(),
            history.diffusivities(),
            config.momentum.maximum_eddy_to_molecular,
            config.heat.maximum_diffusivity_to_molecular,
            gate,
        )?;
        let report = ThermalLagrangianSgsReport {
            momentum_fit: history.momentum_report(),
            heat_fit: history.scalar_report(),
            transport,
            momentum_memory_time_range: history.momentum_memory_time_range(),
        };
        Ok(ThermalLagrangianSgsTrial {
            history,
            momentum: &workspace.momentum,
            thermal: &workspace.thermal,
            report,
        })
    }

    /// Advance the existing source/flow/heat owner and both pathline histories.
    /// Native balance, thermal interval, endpoint Mach and cancellation checks
    /// all precede physical publication. A later refusal leaves history at the
    /// last accepted sample, so retry uses the same elapsed physical interval.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_sources_and_lagrangian_sgs(
        &mut self,
        heat_per_step: &[f64],
        forces: &[[f64; 3]],
        sample_time_lu: f64,
        config: ThermalLagrangianSgsConfig,
        molecular: Option<RelaxationField<'_>>,
        workspace: &mut ThermalLagrangianSgsWorkspace3,
        limits: ThermalSourceLimits,
        gate: &CancelGate,
    ) -> Result<ThermalLagrangianSgsStep, ThermalSgsError> {
        let (heat, pending) = step_lagrangian_sgs_mapped(
            self,
            heat_per_step,
            forces,
            sample_time_lu,
            config,
            molecular,
            workspace,
            limits,
            gate,
            Ok,
        )?;
        let sgs = pending.report();
        pending.commit();
        Ok(ThermalLagrangianSgsStep { heat, sgs })
    }
}

// The returned history token deliberately survives the native source step.
// Material compositions keep it until all original endpoint queries and their
// outer publication have succeeded. There is no fallible work after commit.
#[allow(clippy::too_many_arguments)]
pub(crate) fn step_lagrangian_sgs_mapped<'a, R>(
    state: &mut ThermalLbm3,
    heat: &[f64],
    forces: &[[f64; 3]],
    sample_time_lu: f64,
    config: ThermalLagrangianSgsConfig,
    molecular: Option<RelaxationField<'_>>,
    workspace: &'a mut ThermalLagrangianSgsWorkspace3,
    limits: ThermalSourceLimits,
    gate: &CancelGate,
    finish: impl FnOnce(ThermalSourceBalance) -> Result<R, ThermalSourceError>,
) -> Result<(R, ThermalLagrangianSgsTrial<'a>), ThermalSgsError> {
    let pending =
        state.prepare_thermal_lagrangian_sgs(config, sample_time_lu, molecular, workspace, gate)?;
    let mut admission = Ok(());
    let result = source_step_transport_mapped::<ThermalLbm3, 3, 7, _>(
        state,
        heat,
        forces,
        limits,
        Some(pending.relaxation()),
        |candidate, balance| {
            admission = admit_endpoint(candidate, config.maximum_mach, gate);
            if admission.is_err() {
                return Err(ThermalSourceError::Invalid(
                    "Lagrangian SGS endpoint admission failed",
                ));
            }
            finish(balance)
        },
    );
    admission?;
    Ok((result?, pending))
}

#[cfg(test)]
mod tests;
