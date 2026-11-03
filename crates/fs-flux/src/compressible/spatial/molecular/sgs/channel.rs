//! Plane-channel dynamic stress and independent heat on the original FV owner.
//!
//! Only the two homogeneous periodic directions are test-filtered. The current
//! physical no-slip wall data determine the shared half-cell derivatives.
//! Interior modeled transport uses the original stress/work/heat faces and
//! parabolic estimate; impermeable wall fluxes remain purely molecular.

use super::*;
use fs_lbm::sgs3::dynamic::DynamicSgsError;
use fs_lbm::sgs3::dynamic::channel::{ChannelDynamicScalar, ChannelDynamicSgs, ScalarWall};

/// Independent plane-averaged stress and temperature diffusion controls.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ChannelDynamicThermalSubgridConfig {
    pub momentum: ChannelDynamicSgsConfig,
    pub heat: ChannelDynamicScalarConfig,
}

impl ChannelDynamicThermalSubgridConfig {
    pub fn validate(self) -> Result<(), MolecularError> {
        self.momentum.validate()?;
        self.heat.validate()?;
        Ok(())
    }
}

/// Common fixed Cartesian channel admission for frozen and evolving mixtures.
pub(in crate::compressible::spatial) fn channel_axis(
    grid: &CartesianGrid,
    boundaries: Boundaries3d,
) -> Result<usize, DynamicSgsError> {
    if grid.shape.iter().any(|&n| n < 3) {
        return Err(DynamicSgsError::Invalid(
            "channel transport requires three cells per axis",
        ));
    }
    let mut normal = None;
    for axis in 0..3 {
        match boundaries.axis(axis) {
            AxisBoundaries::Periodic => {}
            AxisBoundaries::Ends {
                lower: BoundaryCondition::Reflecting,
                upper: BoundaryCondition::Reflecting,
            } if normal.is_none() => normal = Some(axis),
            _ => {
                return Err(DynamicSgsError::Invalid(
                    "channel transport requires two periodic axes and one reflecting wall pair",
                ));
            }
        }
    }
    normal.ok_or(DynamicSgsError::Invalid(
        "channel transport requires a physical wall pair",
    ))
}

/// Read current owner wall data, never a stale copy from filter construction.
pub(in crate::compressible::spatial) fn channel_wall_values(
    normal_axis: usize,
    walls: &[MolecularWall3d; 6],
) -> Result<([[f64; 3]; 2], [ScalarWall; 2]), DynamicSgsError> {
    if normal_axis >= 3 {
        return Err(DynamicSgsError::Invalid("invalid channel wall axis"));
    }
    let mut velocity = [[0.0; 3]; 2];
    let mut scalar = [ScalarWall::ZeroGradient; 2];
    for side in 0..2 {
        let MolecularWall3d::NoSlip {
            velocity: u,
            thermal,
        } = walls[2 * normal_axis + side]
        else {
            return Err(DynamicSgsError::Invalid(
                "channel dynamic transport requires two explicit no-slip walls",
            ));
        };
        if u.iter().any(|v| !v.is_finite()) || u[normal_axis] != 0.0 {
            return Err(DynamicSgsError::Invalid(
                "channel wall velocity must be finite and tangential",
            ));
        }
        velocity[side] = u;
        scalar[side] = match thermal {
            WallThermalCondition::Adiabatic => ScalarWall::ZeroGradient,
            WallThermalCondition::Isothermal(t) if positive(t.value()) => {
                ScalarWall::Value(t.value())
            }
            _ => {
                return Err(DynamicSgsError::Invalid(
                    "positive finite channel wall temperature required",
                ));
            }
        };
    }
    Ok((velocity, scalar))
}

impl<L: TransportLaw> MolecularTransport3d<L> {
    /// Added numerical capacity: 256 bytes/cell, excluding original molecular
    /// and Euler buffers, fixed metadata and allocator overhead.
    pub fn required_channel_dynamic_thermal_sgs_workspace_bytes(
        flow: &Euler3d,
    ) -> Result<usize, MolecularError> {
        thermal::required_bytes(flow.states.len())
    }

    /// Fit momentum and heat independently at every original SSP stage and
    /// retry, using the explicit lower/upper no-slip wall pair. The wall axis
    /// is inferred from the unique reflecting pair; the other axes are periodic.
    /// Both stationary/moving tangential walls and adiabatic/isothermal heat
    /// conditions retain their original material and half-cell flux owner.
    ///
    /// This is a resolved-wall, low-Mach exploration closure. It supplies no
    /// coarse-grid wall function, shock model or full compressible SGS energy
    /// equation. All added arrays obey the byte budget and are reused.
    pub fn with_channel_dynamic_thermal_sgs(
        flow: &Euler3d,
        law: L,
        config: ChannelDynamicThermalSubgridConfig,
        walls: [MolecularWall3d; 2],
        maximum_subgrid_workspace_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, MolecularError> {
        poll(&mut cancelled)?;
        config.validate()?;
        let normal = channel_axis(&flow.grid, flow.boundaries)?;
        dynamic::admit_bytes(
            thermal::required_bytes(flow.states.len())?,
            maximum_subgrid_workspace_bytes,
        )?;
        let mut result = Self::new(flow, law, &mut cancelled)?;
        for (side, wall) in walls.into_iter().enumerate() {
            result = result.with_wall(BoundarySide::ALL[2 * normal + side], wall)?;
        }
        result.sgs = Some(SubgridTransport::new_channel_dynamic_thermal(
            flow,
            config,
            &result.walls,
            maximum_subgrid_workspace_bytes,
        )?);
        poll(&mut cancelled)?;
        Ok(result)
    }

    /// Exact channel controls. Other closure geometries return None.
    pub fn channel_dynamic_thermal_subgrid_config(
        &self,
    ) -> Option<ChannelDynamicThermalSubgridConfig> {
        match &self.sgs.as_ref()?.closure {
            SubgridClosure::ChannelDynamicThermal(dynamic) => Some(dynamic.config),
            _ => None,
        }
    }
}

pub(super) struct ChannelDynamicThermalTransport {
    pub(super) config: ChannelDynamicThermalSubgridConfig,
    momentum: ChannelDynamicSgs,
    heat: ChannelDynamicScalar,
    work: thermal::ThermalWorkspace,
}

impl SubgridTransport {
    pub(in crate::compressible::spatial::molecular) fn new_channel_dynamic_thermal(
        flow: &Euler3d,
        config: ChannelDynamicThermalSubgridConfig,
        walls: &[MolecularWall3d; 6],
        maximum_bytes: usize,
    ) -> Result<Self, MolecularError> {
        config.validate()?;
        let normal = channel_axis(&flow.grid, flow.boundaries)?;
        channel_wall_values(normal, walls)?;
        let n = flow.states.len();
        dynamic::admit_bytes(thermal::required_bytes(n)?, maximum_bytes)?;
        let closure = ChannelDynamicThermalTransport {
            config,
            momentum: ChannelDynamicSgs::new_cartesian(
                flow.grid.shape,
                flow.grid.spacing,
                normal,
                maximum_bytes,
            )?,
            heat: ChannelDynamicScalar::new(flow.grid.shape, normal, maximum_bytes)?,
            work: thermal::ThermalWorkspace::new(n)?,
        };
        let dynamic_viscosity = allocated(n, 0.0)?;
        let bytes = closure.momentum.workspace_bytes()
            + closure.heat.workspace_bytes()
            + closure.work.workspace_bytes()
            + dynamic_viscosity.capacity() * core::mem::size_of::<f64>();
        dynamic::admit_bytes(bytes, maximum_bytes)?;
        Ok(Self {
            closure: SubgridClosure::ChannelDynamicThermal(closure),
            dynamic_viscosity,
        })
    }
}

impl ChannelDynamicThermalTransport {
    pub(super) fn prepare(
        &mut self,
        primitive: &[Primitive],
        walls: &[MolecularWall3d; 6],
        dynamic_viscosity: &mut [f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), MolecularError> {
        let (wall_velocity, wall_scalar) = channel_wall_values(self.momentum.normal_axis(), walls)?;
        self.work.stage(primitive, cancelled)?;
        self.momentum.prepare_with_scalar(
            self.config.momentum,
            self.config.heat,
            &self.work.density,
            &self.work.velocity,
            wall_velocity,
            &self.work.temperature,
            wall_scalar,
            &mut self.heat,
            &mut *cancelled,
        )?;
        self.work.store_transport(
            primitive,
            self.momentum.viscosities(),
            self.heat.diffusivities(),
            dynamic_viscosity,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn face_coefficients(
        &self,
        gas: &FrozenGas,
        left: usize,
        right: usize,
        axis: usize,
        face: usize,
        temperature: f64,
        pressure: f64,
        mu_t: f64,
        molecular_mu: f64,
        molecular_conductivity: f64,
    ) -> Result<(f64, f64), MolecularError> {
        self.work.face_coefficients(
            gas,
            left,
            right,
            axis,
            face,
            temperature,
            pressure,
            mu_t,
            molecular_mu,
            molecular_conductivity,
            self.config.momentum.maximum_eddy_to_molecular,
            self.config.heat.maximum_diffusivity_to_molecular,
        )
    }
}
