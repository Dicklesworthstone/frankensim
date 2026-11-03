//! Same-state dynamic momentum and heat on the original D3Q19/D3Q7 slab.
//!
//! Only the two periodic tangent directions are test-filtered. Each z layer
//! fits its own coefficient, with the actual stationary halfway velocity and
//! fixed-temperature wall conditions in normal derivatives. The original
//! source owner measures wall heat and publishes the complete candidate.
//! This remains resolved-wall, low-Mach Boussinesq sensible-temperature
//! transport; modeled mechanical dissipation is diagnostic, not a heat source.

use super::*;
use crate::sgs3::dynamic::channel::{
    ChannelDynamicScalar, ChannelDynamicScalarConfig, ChannelDynamicSgs, ChannelDynamicSgsConfig,
    ScalarWall,
};
use crate::sgs3::dynamic::{DynamicSgsError, DynamicSgsReport, scalar::DynamicScalarReport};

/// Independent plane-averaged fits. No turbulent Prandtl number is prescribed.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermalChannelSgsConfig {
    pub momentum: ChannelDynamicSgsConfig,
    pub heat: ChannelDynamicScalarConfig,
    /// Flow-lattice Mach bound in `(0, 0.3]`, checked before publication.
    pub maximum_mach: f64,
}

impl ThermalChannelSgsConfig {
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

impl From<DynamicSgsError> for ThermalSgsError {
    fn from(error: DynamicSgsError) -> Self {
        match error {
            DynamicSgsError::Cancelled => Self::Cancelled,
            error => Self::Dynamic(error),
        }
    }
}

/// Reusable fluid-only macroscopic/filter buffers and complete relaxation
/// fields. The byte budget excludes the original source owner's temporary
/// complete solver snapshot; its cell budget remains independently enforced.
pub struct ThermalChannelSgsWorkspace3 {
    velocity: Vec<[f64; 3]>,
    density: Vec<f64>,
    temperature: Vec<f64>,
    momentum: Vec<f64>,
    thermal: Vec<f64>,
    filter: ChannelDynamicSgs,
    scalar: ChannelDynamicScalar,
}

impl ThermalChannelSgsWorkspace3 {
    /// Step-local total rates, consumable only after successful preparation.
    /// These must never be supplied as molecular rates on the following step.
    pub fn relaxation(&self) -> RelaxationField<'_> {
        RelaxationField {
            momentum: &self.momentum,
            thermal: &self.thermal,
        }
    }

    /// Actual numerical capacities, excluding fixed metadata and allocator overhead.
    pub fn workspace_bytes(&self) -> usize {
        (3 * self.velocity.capacity()
            + self.density.capacity()
            + self.temperature.capacity()
            + self.momentum.capacity()
            + self.thermal.capacity())
            * core::mem::size_of::<f64>()
            + self.filter.workspace_bytes()
            + self.scalar.workspace_bytes()
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ThermalChannelSgsReport {
    pub momentum_fit: DynamicSgsReport,
    pub heat_fit: DynamicScalarReport,
    pub transport: ThermalSgsReport,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermalChannelSgsStep {
    pub heat: ThermalSourceBalance,
    pub sgs: ThermalChannelSgsReport,
}

fn admit_bytes(required: usize, maximum: usize) -> Result<(), ThermalSgsError> {
    if required > maximum {
        Err(DynamicSgsError::WorkspaceBytes { required, maximum }.into())
    } else {
        Ok(())
    }
}

impl ThermalLbm3 {
    /// Numerical payload for the two nonlocal filters, current rho/u/T and
    /// full step-local rates. Fixed metadata and the source snapshot are separate.
    pub fn required_thermal_channel_sgs_workspace_bytes(&self) -> Result<usize, ThermalSgsError> {
        let (cells, plane) = self.layout()?;
        let fluid = cells - 2 * plane;
        let required = fluid
            .checked_mul(5 * core::mem::size_of::<f64>())
            .and_then(|n| {
                cells
                    .checked_mul(2 * core::mem::size_of::<f64>())
                    .and_then(|m| n.checked_add(m))
            })
            .and_then(|n| {
                ChannelDynamicSgs::required_bytes(fluid)
                    .ok()
                    .and_then(|m| n.checked_add(m))
            })
            .and_then(|n| {
                ChannelDynamicScalar::required_bytes(fluid)
                    .ok()
                    .and_then(|m| n.checked_add(m))
            })
            .ok_or(ThermalSgsError::Invalid(
                "channel SGS workspace size overflow",
            ))?;
        Ok(required)
    }

    /// Allocate once for the existing periodic-x/y, rigid-z slab. At least
    /// three cells in each tangent direction and two fluid layers are required.
    pub fn thermal_channel_sgs_workspace(
        &self,
        maximum_bytes: usize,
    ) -> Result<ThermalChannelSgsWorkspace3, ThermalSgsError> {
        let (cells, plane) = self.layout()?;
        let fluid = cells - 2 * plane;
        admit_bytes(
            self.required_thermal_channel_sgs_workspace_bytes()?,
            maximum_bytes,
        )?;
        let shape = [self.grid.nx, self.grid.ny, self.grid.nz - 2];
        let result = ThermalChannelSgsWorkspace3 {
            filter: ChannelDynamicSgs::new(shape, 1.0, 2, maximum_bytes)?,
            scalar: ChannelDynamicScalar::new(shape, 2, maximum_bytes)?,
            velocity: filled(fluid, [0.0; 3])?,
            density: filled(fluid, 0.0)?,
            temperature: filled(fluid, 0.0)?,
            momentum: filled(cells, 0.0)?,
            thermal: filled(cells, 0.0)?,
        };
        admit_bytes(result.workspace_bytes(), maximum_bytes)?;
        Ok(result)
    }

    /// Prepare independent momentum and heat diffusion from the SAME accepted
    /// Guo-centered velocity, density and temperature. Neither fit wraps across
    /// the plates. Molecular rates come from the caller or original slab;
    /// coefficients are frozen for one source/flow/temperature step.
    ///
    /// Zero viscosity and zero diffusivity each preserve their original rate
    /// bits independently. Preparation allocates nothing and polls cancellation
    /// every 64 cells, including shared filtering. Failed scratch is overwritten
    /// on retry; physical state is never modified here.
    pub fn prepare_thermal_channel_sgs(
        &self,
        config: ThermalChannelSgsConfig,
        molecular: Option<RelaxationField<'_>>,
        workspace: &mut ThermalChannelSgsWorkspace3,
        gate: &CancelGate,
    ) -> Result<ThermalChannelSgsReport, ThermalSgsError> {
        poll(gate)?;
        config.validate()?;
        let (cells, plane) = self.layout()?;
        let fluid = cells - 2 * plane;
        if workspace.filter.shape() != [self.grid.nx, self.grid.ny, self.grid.nz - 2]
            || workspace.velocity.len() != fluid
            || workspace.density.len() != fluid
            || workspace.temperature.len() != fluid
            || workspace.momentum.len() != cells
            || workspace.thermal.len() != cells
        {
            return Err(ThermalSgsError::Invalid(
                "channel SGS workspace does not match the slab",
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
            workspace.momentum[cell] =
                molecular.map_or(self.grid.tau[cell], |field| field.momentum[cell]);
            workspace.thermal[cell] = molecular.map_or(self.tau_g, |field| field.thermal[cell]);
            if cell >= plane && cell < cells - plane {
                let active = cell - plane;
                let (rho, velocity) = self.grid.moments3(cell);
                admit_velocity(rho, velocity, cell, config.maximum_mach)?;
                workspace.density[active] = rho;
                workspace.velocity[active] = velocity;
                workspace.temperature[active] = self.gpop[cell].iter().sum();
            }
        }
        let (momentum_fit, heat_fit) = workspace.filter.prepare_with_scalar(
            config.momentum,
            config.heat,
            &workspace.density,
            &workspace.velocity,
            [[0.0; 3]; 2],
            &workspace.temperature,
            [
                ScalarWall::Value(self.t_bottom),
                ScalarWall::Value(self.t_top),
            ],
            &mut workspace.scalar,
            || gate.is_requested(),
        )?;
        let mut transport = ThermalSgsReport::default();
        for active in 0..fluid {
            if active % 64 == 0 {
                poll(gate)?;
            }
            let cell = active + plane;
            let nu = (workspace.momentum[cell] - 0.5) / 3.0;
            let alpha = (workspace.thermal[cell] - 0.5) / 4.0;
            let eddy = workspace.filter.viscosities()[active];
            let diffusivity = workspace.scalar.diffusivities()[active];
            let ratio = eddy / nu;
            if !ratio.is_finite() || ratio > config.momentum.maximum_eddy_to_molecular {
                return Err(SgsError::ViscosityRatioExceeded {
                    cell,
                    ratio,
                    maximum: config.momentum.maximum_eddy_to_molecular,
                }
                .into());
            }
            let heat_ratio = diffusivity / alpha;
            if !heat_ratio.is_finite() || heat_ratio > config.heat.maximum_diffusivity_to_molecular
            {
                return Err(ThermalSgsError::DiffusivityRatio {
                    cell,
                    ratio: heat_ratio,
                    maximum: config.heat.maximum_diffusivity_to_molecular,
                });
            }
            if eddy != 0.0 {
                workspace.momentum[cell] = 0.5 + 3.0 * (nu + eddy);
            }
            if diffusivity != 0.0 {
                workspace.thermal[cell] = 0.5 + 4.0 * (alpha + diffusivity);
            }
            if !workspace.momentum[cell].is_finite() || !workspace.thermal[cell].is_finite() {
                return Err(ThermalSgsError::Invalid(
                    "total channel SGS transport is not representable",
                ));
            }
            let gradient = self.thermal_sgs_gradient(cell, |i| workspace.velocity[i - plane]);
            let dilation_third = (gradient[0][0] + gradient[1][1] + gradient[2][2]) / 3.0;
            let mut strain2 = 0.0;
            for i in 0..3 {
                for j in 0..3 {
                    let s = 0.5 * (gradient[i][j] + gradient[j][i])
                        - if i == j { dilation_third } else { 0.0 };
                    strain2 += s * s;
                }
            }
            let factor = 2.0 * workspace.density[active] * strain2;
            transport.momentum.fluid_cells += 1;
            transport.momentum.resolved_dissipation_lu += factor * nu;
            transport.momentum.sgs_dissipation_lu += factor * eddy;
            transport.momentum.maximum_eddy_viscosity_lu =
                transport.momentum.maximum_eddy_viscosity_lu.max(eddy);
            transport.momentum.maximum_eddy_to_molecular =
                transport.momentum.maximum_eddy_to_molecular.max(ratio);
            transport.maximum_eddy_diffusivity_lu =
                transport.maximum_eddy_diffusivity_lu.max(diffusivity);
        }
        poll(gate)?;
        if !transport.momentum.resolved_dissipation_lu.is_finite()
            || !transport.momentum.sgs_dissipation_lu.is_finite()
        {
            return Err(SgsError::NonFiniteViscosity.into());
        }
        Ok(ThermalChannelSgsReport {
            momentum_fit,
            heat_fit,
            transport,
        })
    }

    /// Advance through the original transactional source and wall-heat owner.
    /// Arrays include both wall planes; their heating must remain zero. Both
    /// transport fields are borrowed for this step only. The native source and
    /// wall balance, temperature interval, Mach and cancellation checks all
    /// precede publication. The serial collision kernel retains its original
    /// cancellation granularity (checks around, not inside, the kernel).
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_sources_and_channel_sgs(
        &mut self,
        heat_per_step: &[f64],
        forces: &[[f64; 3]],
        config: ThermalChannelSgsConfig,
        molecular: Option<RelaxationField<'_>>,
        workspace: &mut ThermalChannelSgsWorkspace3,
        limits: ThermalSourceLimits,
        gate: &CancelGate,
    ) -> Result<ThermalChannelSgsStep, ThermalSgsError> {
        let (heat, sgs) = step_channel_sgs_mapped(
            self,
            heat_per_step,
            forces,
            config,
            molecular,
            workspace,
            limits,
            gate,
            Ok,
        )?;
        Ok(ThermalChannelSgsStep { heat, sgs })
    }
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn step_channel_sgs_mapped<R>(
    state: &mut ThermalLbm3,
    heat: &[f64],
    forces: &[[f64; 3]],
    config: ThermalChannelSgsConfig,
    molecular: Option<RelaxationField<'_>>,
    workspace: &mut ThermalChannelSgsWorkspace3,
    limits: ThermalSourceLimits,
    gate: &CancelGate,
    finish: impl FnOnce(ThermalSourceBalance) -> Result<R, ThermalSourceError>,
) -> Result<(R, ThermalChannelSgsReport), ThermalSgsError> {
    let report = state.prepare_thermal_channel_sgs(config, molecular, workspace, gate)?;
    let mut admission = Ok(());
    let result = source_step_transport_mapped::<ThermalLbm3, 3, 7, _>(
        state,
        heat,
        forces,
        limits,
        Some(workspace.relaxation()),
        |candidate, balance| {
            admission = (|| {
                let (cells, plane) = candidate.layout()?;
                for cell in plane..cells - plane {
                    if cell % 64 == 0 {
                        poll(gate)?;
                    }
                    let (rho, u) = candidate.grid.moments3(cell);
                    admit_velocity(rho, u, cell, config.maximum_mach)?;
                }
                poll(gate)
            })();
            if admission.is_err() {
                return Err(ThermalSourceError::Invalid(
                    "channel SGS endpoint admission failed",
                ));
            }
            finish(balance)
        },
    );
    admission?;
    Ok((result?, report))
}

#[cfg(test)]
mod tests;
