//! Physical-wall Lagrangian turbulence on the original conservative gas split.

use super::*;

impl<L: TransportLaw> MolecularEuler3d<L> {
    /// Attach pathline momentum history with fixed turbulent Prandtl heat to a
    /// resolved no-slip channel. Exactly one Euler axis reflects and the other
    /// two are periodic; every dimension needs at least three fluid cells.
    /// The explicit lower/upper walls retain their original velocity, work and
    /// isothermal/adiabatic heat laws. Later wall replacements feed the next
    /// accepted-start fit.
    ///
    /// Tangential histories wrap periodically. Normal departures inside the
    /// outer half cells use nearest-fluid-center extension; crossing a physical
    /// wall refuses without publishing history or a partial physical step.
    /// The same complete split, retry and diagnostic transaction owns memory.
    /// This first-order exploratory model adds no wall law or inflow history.
    ///
    /// Additional numerical payload remains 208 bytes/cell. The cap excludes
    /// the original molecular and Euler buffers, as for the periodic entry.
    pub fn with_lagrangian_dynamic_channel_sgs(
        flow: Euler3d,
        law: L,
        config: sgs::LagrangianSubgridConfig,
        walls: [MolecularWall3d; 2],
        maximum_subgrid_workspace_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, MolecularError> {
        poll(&mut cancelled)?;
        let subgrid = sgs::SubgridTransport::new_lagrangian_channel(
            &flow,
            config,
            maximum_subgrid_workspace_bytes,
        )?;
        Self::with_channel_lagrangian_owner(flow, law, subgrid, walls, &mut cancelled)
    }

    /// Attach independent momentum and temperature histories to the same
    /// physical channel. Both fits use current accepted-start wall data and
    /// publish only after complete molecular/Euler/molecular acceptance.
    ///
    /// SSP/RKL stages and whole-split retries hold both coefficients while the
    /// original conservative faces reevaluate current density, strain and
    /// source Cp. Wall work/heat remain molecular unless an equilibrium stress
    /// model is explicitly attached. Diagnostic vorticity queries use the same
    /// walls-aware history trial and leave accepted memory and reports intact.
    ///
    /// The first-order boundary treatment matches
    /// `with_lagrangian_dynamic_channel_sgs`. Additional numerical payload is
    /// 320 bytes/cell, with no new history or staging allocation during stepping.
    pub fn with_lagrangian_dynamic_channel_thermal_sgs(
        flow: Euler3d,
        law: L,
        config: sgs::LagrangianThermalSubgridConfig,
        walls: [MolecularWall3d; 2],
        maximum_subgrid_workspace_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, MolecularError> {
        poll(&mut cancelled)?;
        let subgrid = sgs::SubgridTransport::new_lagrangian_channel_thermal(
            &flow,
            config,
            maximum_subgrid_workspace_bytes,
        )?;
        Self::with_channel_lagrangian_owner(flow, law, subgrid, walls, &mut cancelled)
    }

    /// Physical wall normal for a channel history, or None for other geometries.
    /// Existing Lagrangian config and report getters apply to either geometry.
    pub fn lagrangian_channel_normal_axis(&self) -> Option<usize> {
        self.molecular
            .sgs
            .as_ref()
            .and_then(sgs::SubgridTransport::lagrangian_channel_normal_axis)
    }

    fn with_channel_lagrangian_owner(
        flow: Euler3d,
        law: L,
        subgrid: sgs::SubgridTransport,
        walls: [MolecularWall3d; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, MolecularError> {
        let normal = subgrid
            .lagrangian_channel_normal_axis()
            .ok_or(MolecularError::Invalid(
                "Lagrangian channel requires a physical wall pair",
            ))?;
        if walls
            .iter()
            .any(|wall| !matches!(wall, MolecularWall3d::NoSlip { .. }))
        {
            return Err(MolecularError::Invalid(
                "Lagrangian channel requires two explicit no-slip walls",
            ));
        }
        poll(cancelled)?;
        let mut result = Self::new(flow, law, &mut *cancelled)?;
        for (end, wall) in walls.into_iter().enumerate() {
            result = result.with_wall(BoundarySide::ALL[2 * normal + end], wall)?;
        }
        result.molecular.sgs = Some(subgrid);
        poll(cancelled)?;
        Ok(result)
    }
}
