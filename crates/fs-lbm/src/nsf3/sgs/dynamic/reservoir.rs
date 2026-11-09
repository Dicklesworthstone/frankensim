//! Instantaneous transverse Favre--Germano transport at prescribed reservoirs.
//!
//! The native full-face ghost remains the only pressure, mass-flow, material
//! and total-energy boundary authority. This attachment supplies current
//! deviatoric eddy viscosity and independently fitted heat diffusivity; it
//! does not replace a reservoir by an impermeable wall or invent inflow noise.

use super::*;
use crate::sgs3::dynamic::reservoir::{ReservoirDynamicScalar, ReservoirDynamicSgs, ReservoirPrimitive};

pub(super) struct ReservoirWorkspace {
    momentum: ReservoirDynamicSgs,
    heat: ReservoirDynamicScalar,
}

impl<L: Lattice3<27>> NsfGrid3<L> {
    /// Fit momentum and heat on the current native field between two complete
    /// prescribed reservoir faces. Exactly one axis is nonperiodic; both of
    /// its faces must already be installed, every cell must be fluid, and the
    /// transverse axes must be periodic. At least three cells per axis are
    /// required. Characteristic faces and localized apertures are excluded.
    ///
    /// Nine-cell test filtering acts only in the transverse directions, with
    /// squared volume-width ratio 9^(2/3). Domain controls average contractions
    /// on each transverse plane; LocalBox controls average their transverse
    /// nine-cell neighborhood. Normal derivatives use the actual full-spacing
    /// ghost primitive values, including throughflow and reservoir temperature.
    /// Spatial face profiles retain their individual physical ghosts and use
    /// the same transverse Favre weights for their test-level ghost values.
    /// This is an explicit transverse LES model, not a general 3-D open filter.
    ///
    /// Ghost changes on the admitted axis re-fit both closures within the
    /// original force-rebase/energy-flux transaction. Removing a reservoir or
    /// changing boundary topology refuses. Shared-face and reservoir molecular
    /// heat-ratio gates retain their native material state. Adjacent modeled
    /// conductivity is extrapolated to each reservoir face exactly as for the
    /// existing local SGS owner; no exterior eddy viscosity is prescribed.
    /// The numerical filter payload is 200 nominal bytes per fluid cell, with
    /// no halo arrays. This adds no characteristic/radiation condition, inflow
    /// turbulence spectrum, shock model or complete compressible SGS energy.
    pub fn with_reservoir_dynamic_thermal_sgs(
        mut self,
        config: NsfDynamicThermalSgsConfig,
        maximum_workspace_bytes: usize,
    ) -> Result<Self, NsfError> {
        config.validate()?;
        if self.sgs.is_some() || self.dynamic_sgs.is_some() {
            return Err(NsfError::Input(
                "attach reservoir dynamic SGS to a molecular NSF grid",
            ));
        }
        let mut nonperiodic = self
            .spec
            .periodic
            .iter()
            .enumerate()
            .filter_map(|(axis, periodic)| (!periodic).then_some(axis));
        let normal_axis = nonperiodic.next().ok_or(NsfError::Input(
            "reservoir dynamic SGS requires one nonperiodic axis",
        ))?;
        if nonperiodic.next().is_some()
            || self.spec.dims.iter().any(|&n| n < 3)
            || self.solid.iter().any(|solid| *solid)
            || self.open.is_some()
            || !self.reservoir_apertures.is_empty()
            || self.reservoir_faces[normal_axis]
                .iter()
                .any(Option::is_none)
            || self
                .reservoir_faces
                .iter()
                .enumerate()
                .any(|(axis, faces)| axis != normal_axis && faces.iter().any(Option::is_some))
        {
            return Err(NsfError::Input(
                "reservoir dynamic SGS requires two complete reservoirs and periodic transverse fluid axes",
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
        let momentum = ReservoirDynamicSgs::new(
            self.spec.dims,
            [1.0; 3],
            normal_axis,
            maximum_workspace_bytes,
        )
        .map_err(NsfError::DynamicSgs)?;
        let heat =
            ReservoirDynamicScalar::new(self.spec.dims, normal_axis, maximum_workspace_bytes)
                .map_err(NsfError::DynamicSgs)?;
        let bytes = momentum
            .workspace_bytes()
            .checked_add(heat.workspace_bytes())
            .ok_or(NsfError::Input(
                "reservoir dynamic workspace capacity overflow",
            ))?;
        admit(bytes)?;
        self.dynamic_sgs = Some(DynamicTransport {
            state: DynamicState::Reservoir {
                config,
                normal_axis,
                workspace: Arc::new(Mutex::new(ReservoirWorkspace { momentum, heat })),
            },
            bytes,
            sample_time: self.steps as f64,
        });
        self.refresh_dynamic_transport()?;
        Ok(self)
    }

    pub fn reservoir_dynamic_thermal_sgs_config(&self) -> Option<NsfDynamicThermalSgsConfig> {
        match &self.dynamic_sgs.as_ref()?.state {
            DynamicState::Reservoir { config, .. } => Some(*config),
            _ => None,
        }
    }

    // Replacement preserves both complete ghost planes. Original closed
    // dynamic owners and topology removal remain refused by set_open_face.
    pub(in crate::nsf3) fn dynamic_reservoir_face_update_allowed(
        &self,
        axis: usize,
        present: bool,
    ) -> bool {
        present && self.dynamic_sgs.as_ref().is_some_and(|dynamic| {
            matches!(dynamic.state, DynamicState::Reservoir { normal_axis, .. } if normal_axis == axis)
        })
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn reservoir_dynamic_sgs_field(
        &self,
        config: NsfDynamicThermalSgsConfig,
        normal_axis: usize,
        workspace: &Mutex<ReservoirWorkspace>,
        rho: &[f64],
        velocity: &[[f64; 3]],
        theta: &[f64],
        material: &[Option<NsfCellMaterial>],
    ) -> Result<Vec<NsfSgsTransport>, NsfError> {
        let low = self.reservoir_faces[normal_axis][0]
            .ok_or(NsfError::Input("missing admitted low reservoir ghost"))?;
        let high = self.reservoir_faces[normal_axis][1]
            .ok_or(NsfError::Input("missing admitted high reservoir ghost"))?;
        let mut workspace = workspace
            .lock()
            .map_err(|_| NsfError::Input("reservoir NSF workspace lock poisoned"))?;
        let ReservoirWorkspace { momentum, heat } = &mut *workspace;
        if self.reservoir_face_profiles[normal_axis].iter().any(Option::is_some) {
            momentum.prepare_with_profile(
                config.momentum, config.heat, rho, velocity, theta,
                |side, cell| {
                    let ghost = self.reservoir_ghost(cell, normal_axis, side)
                        .expect("both complete reservoir faces were admitted");
                    ReservoirPrimitive { density: ghost.rho, velocity: ghost.u, scalar: ghost.theta }
                }, heat, || false,
            )
        } else {
            momentum.prepare_with_scalar(
                config.momentum,
                config.heat,
                rho,
                velocity,
                [low.u, high.u],
                theta,
                [low.theta, high.theta],
                heat,
                || false,
            )
        }.map_err(NsfError::DynamicSgs)?;
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

#[cfg(test)]
mod tests;
