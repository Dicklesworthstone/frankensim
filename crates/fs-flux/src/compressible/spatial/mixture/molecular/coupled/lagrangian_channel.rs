//! Lagrangian turbulence histories beside the original physical mixture walls.

use super::*;
use crate::compressible::spatial::molecular::sgs::channel_axis;

impl<L: TransportLaw> MolecularEuler3d<L> {
    /// Attach Lagrangian momentum with explicit Prandtl/Schmidt transport to a
    /// channel. Exactly one Euler axis reflects; its supplied lower/upper walls
    /// are no-slip, with finite tangential velocity and original thermal laws.
    /// Every axis requires at least three fluid cells. The additional numerical
    /// payload is the same 216 bytes/cell as the periodic history owner.
    ///
    /// Tangential filtering and pathline interpolation never cross a wall.
    /// Normal departures inside the outer half cell use constant extension of
    /// the nearest fluid history; crossing a physical wall refuses the trial.
    /// This is a first-order, resolved-wall exploration model.
    pub fn with_lagrangian_dynamic_channel_sgs(
        flow: MixtureEuler3d,
        molecular: MolecularTransport3d<L>,
        config: sgs::MixtureLagrangianSubgridConfig,
        walls: [MolecularWall3d; 2],
        maximum_workspace_bytes: usize,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, DissipationError> {
        config.validate()?;
        let extra = Self::required_lagrangian_dynamic_sgs_workspace_bytes(&flow)?;
        Self::with_channel_lagrangian(
            flow,
            molecular,
            walls,
            extra,
            |flow| sgs::SubgridTransport::new_lagrangian_channel(flow, config),
            maximum_workspace_bytes,
            cancelled,
        )
    }

    /// Attach independent Lagrangian stress, temperature and common-species
    /// histories to a physical channel. Uses the same wall admission and
    /// first-order pathline boundary treatment as
    /// `with_lagrangian_dynamic_channel_sgs`.
    ///
    /// One accepted START-state sample sees current wall velocities and
    /// isothermal/adiabatic heat data. Species use interior reconstruction;
    /// impermeability does not impose an invented mass-fraction gradient.
    /// All three coefficients stay fixed through the complete SSP/RKL split
    /// and its retries. The original face owner reevaluates current strain,
    /// material heat capacities, mass-frame closure and partial enthalpies.
    ///
    /// Modeled wall fluxes are zero by default; original molecular wall work,
    /// heat and impermeability remain intact. This starts fresh histories and
    /// adds no inflow, turbulent reaction-rate or SGS-energy closure. Additional
    /// numerical payload is 448 bytes/cell, independent of species count.
    pub fn with_lagrangian_dynamic_channel_species_sgs(
        flow: MixtureEuler3d,
        molecular: MolecularTransport3d<L>,
        config: sgs::MixtureLagrangianSpeciesSubgridConfig,
        walls: [MolecularWall3d; 2],
        maximum_workspace_bytes: usize,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, DissipationError> {
        config.validate()?;
        let extra = Self::required_lagrangian_dynamic_species_sgs_workspace_bytes(&flow)?;
        Self::with_channel_lagrangian(
            flow,
            molecular,
            walls,
            extra,
            |flow| sgs::SubgridTransport::new_lagrangian_channel_species(flow, config),
            maximum_workspace_bytes,
            cancelled,
        )
    }

    /// Physical normal axis of an attached channel history; None for periodic
    /// histories and other SGS models. Coefficient/report getters are shared.
    pub fn lagrangian_channel_normal_axis(&self) -> Option<usize> {
        self.molecular
            .sgs
            .as_ref()
            .and_then(sgs::SubgridTransport::lagrangian_channel_normal_axis)
    }

    #[allow(clippy::too_many_arguments)]
    fn with_channel_lagrangian(
        flow: MixtureEuler3d,
        mut molecular: MolecularTransport3d<L>,
        walls: [MolecularWall3d; 2],
        extra: usize,
        build: impl FnOnce(&MixtureEuler3d) -> Result<sgs::SubgridTransport, DissipationError>,
        maximum_workspace_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, DissipationError> {
        poll(&mut cancelled)?;
        if molecular.sgs.is_some() {
            return Err(DissipationError::Invalid(
                "Lagrangian channel attachment requires a molecular operator without SGS",
            ));
        }
        let normal = channel_axis(&flow.grid, flow.boundaries).map_err(|_| {
            DissipationError::Invalid(
                "Lagrangian channel requires three cells per axis, two periodic axes and one reflecting wall pair",
            )
        })?;
        if walls
            .iter()
            .any(|wall| !matches!(wall, MolecularWall3d::NoSlip { .. }))
        {
            return Err(DissipationError::Invalid(
                "Lagrangian channel requires two explicit no-slip walls",
            ));
        }
        let required = Self::required_workspace_bytes(&flow, &molecular)?
            .checked_add(extra)
            .ok_or(DissipationError::Invalid(
                "coupled Lagrangian channel workspace size overflow",
            ))?;
        admit_workspace(required, maximum_workspace_bytes)?;
        for (end, wall) in walls.into_iter().enumerate() {
            molecular = molecular.with_wall(BoundarySide::ALL[2 * normal + end], wall)?;
        }
        let sgs = build(&flow)?;
        poll(&mut cancelled)?;
        let mut result = Self::new(flow, molecular, maximum_workspace_bytes, &mut cancelled)?;
        result.molecular.sgs = Some(sgs);
        admit_workspace(result.workspace_bytes(), maximum_workspace_bytes)?;
        poll(&mut cancelled)?;
        Ok(result)
    }
}
