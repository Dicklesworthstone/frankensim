//! Accepted-time pathline memory in the original conservative NSF step.
//!
//! Every predictor uses the same accepted signed contractions and endpoint
//! epoch. Only the final admitted force/temperature/material state commits
//! history, together with the physical state. Heat uses its independent fit
//! and actual source cp; stress work stays in the original total-energy flux.

use super::*;
use crate::sgs3::dynamic::lagrangian::scalar::{LagrangianScalarConfig, LagrangianScalarTrial};
use crate::sgs3::dynamic::lagrangian::{
    AdaptiveLagrangianMemory, LagrangianAdvection, LagrangianSgsConfig,
};
use std::sync::MutexGuard;

/// Independent momentum and temperature histories in lattice units.
/// `memory_time` and the optional adaptive minimum are lattice timesteps,
/// not seconds. The original grid clock supplies all sample times.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NsfLagrangianThermalSgsConfig {
    pub momentum: LagrangianSgsConfig,
    pub heat: LagrangianScalarConfig,
    pub adaptive_momentum_memory: Option<AdaptiveLagrangianMemory>,
}

impl NsfLagrangianThermalSgsConfig {
    pub fn validate(self) -> Result<(), NsfError> {
        self.momentum.validate().map_err(NsfError::DynamicSgs)?;
        self.heat.validate().map_err(NsfError::DynamicSgs)?;
        if let Some(memory) = self.adaptive_momentum_memory {
            memory.validate().map_err(NsfError::DynamicSgs)?;
            if memory.minimum_time > self.momentum.memory_time {
                return Err(NsfError::Input(
                    "adaptive NSF minimum memory time exceeds momentum memory time",
                ));
            }
        }
        Ok(())
    }
}

impl DynamicTransport {
    pub(in crate::nsf3) fn begin_native_step(&mut self, steps: usize) -> Result<(), NsfError> {
        if matches!(self.state, DynamicState::Lagrangian { .. }) {
            let next = steps
                .checked_add(1)
                .ok_or(NsfError::Input("Lagrangian NSF clock overflow"))?;
            let time = next as f64;
            if time <= steps as f64 || !time.is_finite() {
                return Err(NsfError::Input(
                    "Lagrangian NSF clock lost one-step resolution",
                ));
            }
            self.sample_time = time;
        }
        Ok(())
    }

    pub(in crate::nsf3) fn finish_native_step(&mut self, steps: usize) {
        self.sample_time = steps as f64;
    }

    pub(in crate::nsf3) fn lagrangian_workspace(
        &self,
    ) -> Result<Option<MutexGuard<'_, LagrangianDynamicScalar>>, NsfError> {
        match &self.state {
            DynamicState::Lagrangian { workspace, .. } => workspace
                .lock()
                .map(Some)
                .map_err(|_| NsfError::Input("Lagrangian NSF workspace lock poisoned")),
            _ => Ok(None),
        }
    }
}

impl<L: Lattice3<27>> NsfGrid3<L> {
    /// Numerical payload for both original filters and all accepted/candidate
    /// history banks: 264 nominal bytes per cell, excluding fixed metadata and
    /// the existing fluid/energy state. Regular native steps reuse this storage.
    pub fn required_lagrangian_sgs_workspace_bytes(&self) -> Result<usize, NsfError> {
        LagrangianDynamicScalar::required_bytes(self.n).map_err(NsfError::DynamicSgs)
    }

    /// Attach pathline-averaged Favre--Germano momentum and heat to the native
    /// conservative NSF owner. All-fluid periodic boxes and resolved channels
    /// with two periodic axes are supported; each dimension needs three cells.
    /// Moving tangential walls and original local thermal laws are retained.
    ///
    /// The initial accepted state seeds local spatial averages. Every later
    /// corrector starts from that SAME accepted history; rejected correctors,
    /// endpoint probes and failed steps never age it. The final admitted state
    /// publishes one joint momentum/heat sample at the original integer clock.
    /// Same-time wall/force refreshes re-emit the accepted contractions without
    /// advancing them. Physical wall heat remains molecular; modeled heat acts
    /// on shared interior faces. The original total-energy stress work is the
    /// only mechanical energy transfer: no extra viscous-heat source is added.
    ///
    /// The pathline/source update is first order and bounded by the authored
    /// per-axis departure limit; it does not inherit the energy trapezoid's
    /// temporal order. This is exploratory low-Mach gradient diffusion, not
    /// a complete compressible LES closure, shock model or DNS qualification.
    /// Openings, interior solids and concurrent local/instantaneous SGS refuse.
    /// Reversible thermal intervals fork their own bounded history workspace.
    pub fn with_lagrangian_thermal_sgs(
        self,
        config: NsfLagrangianThermalSgsConfig,
        maximum_workspace_bytes: usize,
    ) -> Result<Self, NsfError> {
        self.with_lagrangian_thermal_sgs_advection(
            config,
            LagrangianAdvection::Euler,
            maximum_workspace_bytes,
        )
    }

    /// Attach the same accepted-time momentum/heat history with an explicit
    /// replayable characteristic policy. Midpoint substeps sample the current
    /// physical force-corrected velocity along each path and honor the current
    /// moving-wall values. They add no full-grid arrays and do not change the
    /// 264-byte nominal payload, independent relaxation times, source sampling,
    /// molecular wall heat, or joint endpoint publication. Only frozen-field
    /// trajectory accuracy improves; the full unsteady closure remains first
    /// order. Reversible native intervals preserve the selected policy.
    pub fn with_lagrangian_thermal_sgs_advection(
        mut self,
        config: NsfLagrangianThermalSgsConfig,
        advection: LagrangianAdvection,
        maximum_workspace_bytes: usize,
    ) -> Result<Self, NsfError> {
        config.validate()?;
        advection.validate().map_err(NsfError::DynamicSgs)?;
        if self.sgs.is_some() || self.dynamic_sgs.is_some() {
            return Err(NsfError::Input(
                "attach Lagrangian SGS to a molecular NSF grid",
            ));
        }
        let mut walls = self
            .spec
            .periodic
            .iter()
            .enumerate()
            .filter_map(|(axis, &periodic)| (!periodic).then_some(axis));
        let normal_axis = walls.next();
        if walls.next().is_some()
            || self.spec.dims.iter().any(|&n| n < 3)
            || self.solid.iter().any(|&solid| solid)
            || self.open.is_some()
            || self.reservoir_faces.iter().flatten().any(Option::is_some)
            || !self.reservoir_apertures.is_empty()
        {
            return Err(NsfError::Input(
                "Lagrangian NSF SGS requires periodic fluid or one resolved wall axis",
            ));
        }
        let mut workspace = match normal_axis {
            Some(axis) => LagrangianDynamicScalar::new_channel(
                self.spec.dims,
                [1.0; 3],
                axis,
                maximum_workspace_bytes,
            ),
            None => LagrangianDynamicScalar::new(self.spec.dims, [1.0; 3], maximum_workspace_bytes),
        }
        .map_err(NsfError::DynamicSgs)?
        .with_advection(advection)
        .map_err(NsfError::DynamicSgs)?;
        if let Some(memory) = config.adaptive_momentum_memory {
            workspace = workspace
                .with_adaptive_momentum_memory(memory)
                .map_err(NsfError::DynamicSgs)?;
        }
        let bytes = workspace.workspace_bytes();
        self.dynamic_sgs = Some(DynamicTransport {
            state: DynamicState::Lagrangian {
                config,
                normal_axis,
                workspace: Arc::new(Mutex::new(workspace)),
            },
            bytes,
            sample_time: self.steps as f64,
        });
        let (rho, mom) = self.macro_fields(&self.f);
        let velocity = self.velocity_field(&rho, &mom, &self.force);
        let retained = self
            .dynamic_sgs
            .clone()
            .expect("attached Lagrangian NSF SGS");
        let mut workspace = retained.lagrangian_workspace()?;
        let (transport, history) = self.prepare_endpoint_sgs(
            workspace.as_deref_mut(),
            &rho,
            &velocity,
            &self.theta,
            &self.material,
            None,
            || false,
        )?;
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
        if let Some(history) = history {
            history.commit();
        }
        self.sgs_transport = transport;
        self.div = div;
        self.wall_out = wall;
        self.open_out = open;
        Ok(self)
    }

    pub fn lagrangian_thermal_sgs_config(&self) -> Option<NsfLagrangianThermalSgsConfig> {
        match &self.dynamic_sgs.as_ref()?.state {
            DynamicState::Lagrangian { config, .. } => Some(*config),
            _ => None,
        }
    }

    /// Characteristic policy retained by this trajectory; other closure
    /// families and molecular-only grids return None.
    pub fn lagrangian_sgs_advection(&self) -> Result<Option<LagrangianAdvection>, NsfError> {
        self.dynamic_sgs
            .as_ref()
            .map(|dynamic| {
                dynamic
                    .lagrangian_workspace()
                    .map(|workspace| workspace.as_ref().map(|workspace| workspace.advection()))
            })
            .unwrap_or(Ok(None))
    }

    /// Last committed model sample in lattice timesteps. For this native
    /// owner it equals `steps()`; rejected candidates never alter it.
    pub fn lagrangian_sgs_accepted_time(&self) -> Result<Option<f64>, NsfError> {
        self.dynamic_sgs
            .as_ref()
            .map(|dynamic| {
                dynamic.lagrangian_workspace().map(|workspace| {
                    workspace
                        .as_ref()
                        .and_then(|workspace| workspace.accepted_time())
                })
            })
            .unwrap_or(Ok(None))
    }

    /// Actual local momentum memory range at the committed sample, in lattice
    /// timesteps. Independent scalar memory remains the authored heat value.
    pub fn lagrangian_sgs_memory_time_range(&self) -> Result<Option<[f64; 2]>, NsfError> {
        self.dynamic_sgs
            .as_ref()
            .map(|dynamic| {
                dynamic.lagrangian_workspace().map(|workspace| {
                    workspace
                        .as_ref()
                        .and_then(|workspace| workspace.accepted_momentum_memory_time_range())
                })
            })
            .unwrap_or(Ok(None))
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn prepare_lagrangian_trial<'a>(
        &self,
        workspace: &'a mut LagrangianDynamicScalar,
        config: NsfLagrangianThermalSgsConfig,
        normal_axis: Option<usize>,
        time: f64,
        rho: &[f64],
        velocity: &[[f64; 3]],
        theta: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<LagrangianScalarTrial<'a>, NsfError> {
        workspace
            .prepare_or_refresh(
                config.momentum,
                config.heat,
                time,
                rho,
                velocity,
                normal_axis.map(|axis| self.edge_velocity[axis]),
                theta,
                normal_axis.map(|axis| self.channel_scalar_walls(axis)),
                cancelled,
            )
            .map_err(NsfError::DynamicSgs)
    }

    #[allow(clippy::too_many_arguments)]
    pub(in crate::nsf3) fn prepare_endpoint_sgs<'a>(
        &self,
        workspace: Option<&'a mut LagrangianDynamicScalar>,
        rho: &[f64],
        velocity: &[[f64; 3]],
        theta: &[f64],
        material: &[Option<NsfCellMaterial>],
        faces: Option<&OpenFaces>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<(Vec<NsfSgsTransport>, Option<LagrangianScalarTrial<'a>>), NsfError> {
        let Some(workspace) = workspace else {
            return self
                .sgs_field(rho, velocity, theta, material, faces)
                .map(|transport| (transport, None));
        };
        let dynamic = self
            .dynamic_sgs
            .as_ref()
            .ok_or(NsfError::Input("missing Lagrangian NSF controls"))?;
        let DynamicState::Lagrangian {
            config,
            normal_axis,
            ..
        } = &dynamic.state
        else {
            return Err(NsfError::Input("Lagrangian NSF workspace/control mismatch"));
        };
        let trial = self.prepare_lagrangian_trial(
            workspace,
            *config,
            *normal_axis,
            dynamic.sample_time,
            rho,
            velocity,
            theta,
            cancelled,
        )?;
        let transport = self.dynamic_transport_field(
            rho,
            theta,
            material,
            trial.viscosities(),
            trial.diffusivities(),
            config.momentum.maximum_eddy_to_molecular,
        )?;
        Ok((transport, Some(trial)))
    }
}

#[cfg(test)]
mod tests;

#[cfg(test)]
mod advection_tests;
