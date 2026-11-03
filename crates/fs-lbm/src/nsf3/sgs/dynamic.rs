//! Same-state Favre--Germano momentum and temperature transport for NSF3.
//!
//! The original Guo-corrected velocity and current temperature drive both
//! fits. Only the original deviatoric collision rate receives `nu_t`; the
//! original conservative energy faces receive `rho (cp/R) alpha_t`. Thus heat
//! need not vanish when the momentum fit does. This is the shared exploratory
//! low-Mach gradient-diffusion closure, not a complete compressible SGS model.

use super::*;
use crate::sgs3::dynamic::channel::{
    ChannelDynamicScalar, ChannelDynamicScalarConfig, ChannelDynamicSgs, ChannelDynamicSgsConfig,
    ScalarWall,
};
use crate::sgs3::dynamic::scalar::{DynamicScalarConfig, PeriodicDynamicScalar};
use crate::sgs3::dynamic::{DynamicSgsConfig, DynamicSgsError, PeriodicDynamicSgs};
use std::sync::Mutex;

/// Independently fitted momentum and temperature diffusivity. Coefficient and
/// molecular-ratio bounds are admissions, not coefficient clamps.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NsfDynamicThermalSgsConfig {
    pub momentum: DynamicSgsConfig,
    pub heat: DynamicScalarConfig,
}

impl NsfDynamicThermalSgsConfig {
    pub fn validate(self) -> Result<(), NsfError> {
        self.momentum.validate().map_err(NsfError::DynamicSgs)?;
        self.heat.validate().map_err(NsfError::DynamicSgs)
    }
}

/// Independent plane-averaged momentum and heat fits for a resolved channel.
/// The two homogeneous directions are periodic; no test filter crosses a wall.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NsfChannelDynamicThermalSgsConfig {
    pub momentum: ChannelDynamicSgsConfig,
    pub heat: ChannelDynamicScalarConfig,
}

impl NsfChannelDynamicThermalSgsConfig {
    pub fn validate(self) -> Result<(), NsfError> {
        self.momentum.validate().map_err(NsfError::DynamicSgs)?;
        self.heat.validate().map_err(NsfError::DynamicSgs)
    }
}

struct Workspace {
    momentum: PeriodicDynamicSgs,
    heat: PeriodicDynamicScalar,
}

struct ChannelWorkspace {
    momentum: ChannelDynamicSgs,
    heat: ChannelDynamicScalar,
}

enum DynamicState {
    Periodic {
        config: NsfDynamicThermalSgsConfig,
        workspace: Mutex<Workspace>,
    },
    Channel {
        config: NsfChannelDynamicThermalSgsConfig,
        normal_axis: usize,
        workspace: Mutex<ChannelWorkspace>,
    },
}

pub(in crate::nsf3) struct DynamicTransport {
    // Each variant pairs its own controls and filter geometry. One lock around
    // the serial fit preserves Sync for the original parallel cell/face loops.
    state: DynamicState,
    bytes: usize,
}

impl<L: Lattice3<27>> NsfGrid3<L> {
    /// Added reusable filter payload: 200 nominal bytes per cell, excluding
    /// fixed metadata and the existing NSF state/predictor/transport vectors.
    pub fn required_dynamic_sgs_workspace_bytes(&self) -> Result<usize, NsfError> {
        PeriodicDynamicSgs::required_bytes(self.n)
            .and_then(|momentum| {
                momentum
                    .checked_add(PeriodicDynamicScalar::required_bytes(self.n)?)
                    .ok_or(DynamicSgsError::Invalid("NSF dynamic workspace overflow"))
            })
            .map_err(NsfError::DynamicSgs)
    }

    /// Fit periodic subgrid momentum and heat in the existing NSF predictor
    /// and accepted endpoint. Fully periodic, all-fluid domains with at least
    /// three cells per axis are required. The source-owned shear/bulk physics,
    /// force chart, total-energy flux and publication transaction are retained.
    ///
    /// The byte limit covers the two nonlocal filter workspaces. Their arrays
    /// are reused; this does not make the surrounding NSF step allocation-free
    /// or extend its Mach, thermal or stability qualification. Local and dynamic
    /// SGS controls are mutually exclusive. No wall or shock model is supplied.
    pub fn with_dynamic_thermal_sgs(
        mut self,
        config: NsfDynamicThermalSgsConfig,
        maximum_workspace_bytes: usize,
    ) -> Result<Self, NsfError> {
        config.validate()?;
        if self.sgs.is_some() || self.dynamic_sgs.is_some() {
            return Err(NsfError::Input(
                "attach dynamic SGS to a molecular NSF grid",
            ));
        }
        if self.spec.periodic.iter().any(|periodic| !periodic)
            || self.spec.dims.iter().any(|&n| n < 3)
            || self.solid.iter().any(|solid| *solid)
            || self.open.is_some()
            || self.reservoir_faces.iter().flatten().any(Option::is_some)
            || !self.reservoir_apertures.is_empty()
        {
            return Err(NsfError::Input(
                "dynamic NSF SGS requires periodic fluid axes with at least three cells each",
            ));
        }
        let admit = |required| {
            if required > maximum_workspace_bytes {
                Err(NsfError::DynamicSgs(DynamicSgsError::WorkspaceBytes {
                    required,
                    maximum: maximum_workspace_bytes,
                }))
            } else {
                Ok(())
            }
        };
        admit(self.required_dynamic_sgs_workspace_bytes()?)?;
        let momentum = PeriodicDynamicSgs::new(self.spec.dims, 1.0, maximum_workspace_bytes)
            .map_err(NsfError::DynamicSgs)?;
        let heat = PeriodicDynamicScalar::new(self.spec.dims, maximum_workspace_bytes)
            .map_err(NsfError::DynamicSgs)?;
        let bytes = momentum
            .workspace_bytes()
            .checked_add(heat.workspace_bytes())
            .ok_or(NsfError::Input("NSF dynamic workspace capacity overflow"))?;
        admit(bytes)?;
        self.dynamic_sgs = Some(DynamicTransport {
            state: DynamicState::Periodic {
                config,
                workspace: Mutex::new(Workspace { momentum, heat }),
            },
            bytes,
        });
        self.refresh_dynamic_transport()?;
        Ok(self)
    }

    /// Fit momentum and heat separately on each homogeneous channel plane.
    /// Exactly one axis must have halfway walls, and every axis needs at least
    /// three cells. All cells must be fluid. The three wall-normal centers also
    /// support one-sided differentiation of the Favre-filtered adiabatic scalar.
    /// The original stationary or translating tangential wall velocities and
    /// adiabatic/isothermal thermal condition supply the boundary derivatives
    /// on every preparation, including an accepted-time wall-velocity update.
    ///
    /// This resolved-wall exploratory closure adds no wall function. Modeled
    /// heat acts only on interior shared faces; wall conduction remains the
    /// original molecular half-cell law. Moving-wall work uses the original
    /// total stress and energy ledger. Gas sources, bulk relaxation, Guo force
    /// correction and the original temperature/Mach/viscosity bounds remain
    /// authoritative. Open boundaries and interior solids are unsupported.
    /// The reusable numerical workspace costs 200 nominal bytes per fluid cell.
    pub fn with_channel_dynamic_thermal_sgs(
        mut self,
        config: NsfChannelDynamicThermalSgsConfig,
        maximum_workspace_bytes: usize,
    ) -> Result<Self, NsfError> {
        config.validate()?;
        if self.sgs.is_some() || self.dynamic_sgs.is_some() {
            return Err(NsfError::Input(
                "attach dynamic SGS to a molecular NSF grid",
            ));
        }
        let mut walls = self
            .spec
            .periodic
            .iter()
            .enumerate()
            .filter_map(|(axis, periodic)| (!periodic).then_some(axis));
        let normal_axis = walls.next().ok_or(NsfError::Input(
            "channel dynamic NSF SGS requires one wall axis",
        ))?;
        if walls.next().is_some()
            || self.spec.dims.iter().any(|&n| n < 3)
            || self.solid.iter().any(|solid| *solid)
            || self.open.is_some()
            || self.reservoir_faces.iter().flatten().any(Option::is_some)
            || !self.reservoir_apertures.is_empty()
        {
            return Err(NsfError::Input(
                "channel dynamic NSF SGS requires fluid walls and two periodic axes",
            ));
        }
        let admit = |required| {
            if required > maximum_workspace_bytes {
                Err(NsfError::DynamicSgs(DynamicSgsError::WorkspaceBytes {
                    required,
                    maximum: maximum_workspace_bytes,
                }))
            } else {
                Ok(())
            }
        };
        let required = ChannelDynamicSgs::required_bytes(self.n)
            .and_then(|momentum| {
                momentum
                    .checked_add(ChannelDynamicScalar::required_bytes(self.n)?)
                    .ok_or(DynamicSgsError::Invalid("NSF channel workspace overflow"))
            })
            .map_err(NsfError::DynamicSgs)?;
        admit(required)?;
        let momentum =
            ChannelDynamicSgs::new(self.spec.dims, 1.0, normal_axis, maximum_workspace_bytes)
                .map_err(NsfError::DynamicSgs)?;
        let heat = ChannelDynamicScalar::new(self.spec.dims, normal_axis, maximum_workspace_bytes)
            .map_err(NsfError::DynamicSgs)?;
        let bytes = momentum
            .workspace_bytes()
            .checked_add(heat.workspace_bytes())
            .ok_or(NsfError::Input("NSF channel workspace capacity overflow"))?;
        admit(bytes)?;
        self.dynamic_sgs = Some(DynamicTransport {
            state: DynamicState::Channel {
                config,
                normal_axis,
                workspace: Mutex::new(ChannelWorkspace { momentum, heat }),
            },
            bytes,
        });
        self.refresh_dynamic_transport()?;
        Ok(self)
    }

    fn refresh_dynamic_transport(&mut self) -> Result<(), NsfError> {
        let (rho, mom) = self.macro_fields(&self.f);
        let velocity = self.velocity_field(&rho, &mom, &self.force);
        let transport = self.sgs_field(&rho, &velocity, &self.theta, &self.material, None)?;
        let stress = self.stress_field(
            &self.f,
            &self.force,
            &rho,
            &self.theta,
            &self.material,
            &transport,
        );
        let (div, wall, open) = self.flux_divergence(
            &self.e,
            &rho,
            &velocity,
            &self.theta,
            &stress,
            &transport,
            None,
        )?;
        self.sgs_transport = transport;
        self.div = div;
        self.wall_out = wall;
        self.open_out = open;
        Ok(())
    }

    /// Periodic controls; other closure geometries return None.
    pub fn dynamic_thermal_sgs_config(&self) -> Option<NsfDynamicThermalSgsConfig> {
        match &self.dynamic_sgs.as_ref()?.state {
            DynamicState::Periodic { config, .. } => Some(*config),
            DynamicState::Channel { .. } => None,
        }
    }

    /// Channel controls; other closure geometries return None.
    pub fn channel_dynamic_thermal_sgs_config(&self) -> Option<NsfChannelDynamicThermalSgsConfig> {
        match &self.dynamic_sgs.as_ref()?.state {
            DynamicState::Channel { config, .. } => Some(*config),
            DynamicState::Periodic { .. } => None,
        }
    }

    /// Actual retained numerical capacities of the two dynamic filter buffers.
    pub fn dynamic_sgs_workspace_bytes(&self) -> Option<usize> {
        self.dynamic_sgs.as_ref().map(|dynamic| dynamic.bytes)
    }

    pub(super) fn dynamic_sgs_field(
        &self,
        rho: &[f64],
        velocity: &[[f64; 3]],
        theta: &[f64],
        material: &[Option<NsfCellMaterial>],
    ) -> Result<Vec<NsfSgsTransport>, NsfError> {
        let dynamic = self.dynamic_sgs.as_ref().expect("admitted dynamic NSF fit");
        // NSF3 already owns lattice-unit density, physical velocity (including
        // half-force) and temperature. No staging or separate velocity sweep.
        match &dynamic.state {
            DynamicState::Periodic { config, workspace } => {
                let mut workspace = workspace
                    .lock()
                    .map_err(|_| NsfError::Input("dynamic NSF workspace lock poisoned"))?;
                let Workspace { momentum, heat } = &mut *workspace;
                momentum
                    .prepare_with_scalar(
                        config.momentum,
                        config.heat,
                        rho,
                        velocity,
                        theta,
                        heat,
                        || false,
                    )
                    .map_err(NsfError::DynamicSgs)?;
                self.dynamic_transport_field(
                    rho,
                    theta,
                    material,
                    momentum.viscosities(),
                    heat.diffusivities(),
                    config.momentum.maximum_eddy_to_molecular,
                )
            }
            DynamicState::Channel {
                config,
                normal_axis,
                workspace,
            } => {
                let mut workspace = workspace
                    .lock()
                    .map_err(|_| NsfError::Input("channel NSF workspace lock poisoned"))?;
                let ChannelWorkspace { momentum, heat } = &mut *workspace;
                let scalar_wall = match self.spec.wall {
                    WallThermal::Adiabatic => ScalarWall::ZeroGradient,
                    WallThermal::Isothermal(kelvin) => {
                        ScalarWall::Value(self.law.theta_of_kelvin(kelvin))
                    }
                };
                momentum
                    .prepare_with_scalar(
                        config.momentum,
                        config.heat,
                        rho,
                        velocity,
                        self.edge_velocity[*normal_axis],
                        theta,
                        [scalar_wall; 2],
                        heat,
                        || false,
                    )
                    .map_err(NsfError::DynamicSgs)?;
                self.dynamic_transport_field(
                    rho,
                    theta,
                    material,
                    momentum.viscosities(),
                    heat.diffusivities(),
                    config.momentum.maximum_eddy_to_molecular,
                )
            }
        }
    }

    fn dynamic_transport_field(
        &self,
        rho: &[f64],
        theta: &[f64],
        material: &[Option<NsfCellMaterial>],
        eddy_viscosity: &[f64],
        heat_diffusivity: &[f64],
        maximum_eddy_to_molecular: f64,
    ) -> Result<Vec<NsfSgsTransport>, NsfError> {
        (0..self.n)
            .map(|i| {
                let local = material.get(i).and_then(Option::as_ref);
                let molecular =
                    local.map_or_else(|| self.law.mu_lattice(theta[i]), |m| m.mu_lattice) / rho[i];
                if !(molecular.is_finite() && molecular > 0.0) {
                    return Err(NsfError::Sgs(SgsError::InvalidMolecularViscosity {
                        cell: i,
                    }));
                }
                let eddy = eddy_viscosity[i];
                let ratio = eddy / molecular;
                if !ratio.is_finite() || ratio > maximum_eddy_to_molecular {
                    return Err(NsfError::Sgs(SgsError::ViscosityRatioExceeded {
                        cell: i,
                        ratio,
                        maximum: maximum_eddy_to_molecular,
                    }));
                }
                if molecular + eddy > MAX_VERIFIED_LATTICE_VISCOSITY {
                    return Err(NsfError::Viscosity {
                        observed: molecular + eddy,
                        maximum: MAX_VERIFIED_LATTICE_VISCOSITY,
                    });
                }
                let cp_over_r = local.map_or(self.law.gamma / (self.law.gamma - 1.0), |m| {
                    m.state.specific_heat_cp / self.law.gas_constant
                });
                let conductivity = rho[i] * cp_over_r * heat_diffusivity[i];
                if !(conductivity.is_finite() && conductivity >= 0.0) {
                    return Err(NsfError::Input("nonfinite dynamic NSF heat conductivity"));
                }
                // Admit accepted endpoint coefficients from the SAME cached
                // cell source, without another material query. The actual
                // shared-face source is also checked where its flux is formed.
                let molecular_heat = local.map_or_else(
                    || self.law.kappa_lattice(theta[i]),
                    |m| {
                        m.state.thermal_conductivity * self.law.transport_scale()
                            / self.law.gas_constant
                    },
                );
                self.admit_dynamic_heat_ratio(i, conductivity, molecular_heat)?;
                Ok(NsfSgsTransport {
                    eddy_viscosity_lu: eddy,
                    turbulent_conductivity_lu: conductivity,
                })
            })
            .collect()
    }

    // Called after the original shared-face molecular query, before adding
    // modeled conduction. Heat-only fits must reach this gate even when nu_t=0.
    pub(in crate::nsf3) fn admit_dynamic_heat_ratio(
        &self,
        cell: usize,
        modeled: f64,
        molecular: f64,
    ) -> Result<(), NsfError> {
        let Some(dynamic) = &self.dynamic_sgs else {
            return Ok(());
        };
        if modeled == 0.0 {
            return Ok(());
        }
        let ratio = if molecular.is_finite() && molecular > 0.0 {
            modeled / molecular
        } else {
            f64::INFINITY
        };
        let maximum = match &dynamic.state {
            DynamicState::Periodic { config, .. } => config.heat.maximum_diffusivity_to_molecular,
            DynamicState::Channel { config, .. } => config.heat.maximum_diffusivity_to_molecular,
        };
        if !ratio.is_finite() || ratio > maximum {
            Err(NsfError::SgsHeatRatio {
                cell,
                ratio,
                maximum,
            })
        } else {
            Ok(())
        }
    }
}

#[cfg(test)]
mod tests;

#[cfg(test)]
mod channel_tests;
