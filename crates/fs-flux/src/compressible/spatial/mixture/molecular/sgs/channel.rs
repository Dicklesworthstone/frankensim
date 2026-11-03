//! Wall-plane dynamic stress and heat on the evolving-mixture transport owner.
//!
//! The shared channel filter fits momentum and temperature independently in
//! each homogeneous plane, using the actual Cartesian spacings and current
//! no-slip/thermal wall data. Only tangential directions are test filtered;
//! no sample crosses a physical wall. Its discretization and Cabot reference
//! remain with `fs_lbm::sgs3::dynamic::channel`.
//!
//! The original mixture face still owns Cp(T,Y), referenced partial enthalpies,
//! Maxwell--Stefan diffusion, conservative stress work and wall heat exchange.
//! Species SGS mixing uses an explicitly supplied common Schmidt number. This
//! low-Mach, wall-resolved exploration model supplies no wall law, reacting
//! flame closure or missing compressible SGS energy terms.

use super::dynamic::dynamic_error;
use super::*;
use fs_lbm::sgs3::dynamic::channel::{ChannelDynamicScalar, ChannelDynamicSgs, ScalarWall};
pub use fs_lbm::sgs3::dynamic::channel::{ChannelDynamicScalarConfig, ChannelDynamicSgsConfig};

/// Independently fitted channel stress/heat and authored species mixing.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MixtureDynamicChannelThermalSubgridConfig {
    pub momentum: ChannelDynamicSgsConfig,
    pub heat: ChannelDynamicScalarConfig,
    /// Positive finite common species Schmidt number; no species fit is made.
    pub turbulent_schmidt: f64,
}

impl MixtureDynamicChannelThermalSubgridConfig {
    pub fn validate(self) -> Result<(), DissipationError> {
        self.momentum.validate().map_err(dynamic_error)?;
        self.heat.validate().map_err(dynamic_error)?;
        if !positive(self.turbulent_schmidt) {
            return Err(DissipationError::Invalid(
                "positive finite turbulent Schmidt number required",
            ));
        }
        Ok(())
    }
}

impl<L: TransportLaw> MolecularTransport3d<L> {
    /// Complete flow/operator payload: 264 additional nominal numerical bytes
    /// per active fluid cell, as for periodic independent heat fitting.
    pub fn required_workspace_bytes_with_dynamic_channel_thermal_sgs(
        flow: &MixtureEuler3d,
        source: &SpatialTransport<L>,
    ) -> Result<usize, DissipationError> {
        let momentum =
            ChannelDynamicSgs::required_bytes(flow.cell_count()).map_err(dynamic_error)?;
        let heat =
            ChannelDynamicScalar::required_bytes(flow.cell_count()).map_err(dynamic_error)?;
        let staging = flow.cell_count().checked_mul(6 * size_of::<f64>()).ok_or(
            DissipationError::Invalid("dynamic mixture channel workspace overflow"),
        )?;
        Self::required_workspace_bytes_with_sgs(flow, source)?
            .checked_add(momentum)
            .and_then(|bytes| bytes.checked_add(heat))
            .and_then(|bytes| bytes.checked_add(staging))
            .ok_or(DissipationError::Invalid(
                "dynamic mixture channel workspace overflow",
            ))
    }

    /// Fit momentum and heat separately on the current mixture at every SSP
    /// stage/retry. `normal_axis` has two reflecting, explicitly no-slip walls;
    /// the two tangent axes are periodic. `walls` are lower then upper.
    /// At least three fluid cells are required in each tangent direction and
    /// two normal to fixed-temperature walls (three if either is adiabatic).
    /// Wall velocity/temperature may subsequently change through `with_wall`;
    /// preparation always consumes the current original wall data.
    #[allow(clippy::too_many_arguments)]
    pub fn with_dynamic_channel_thermal_sgs(
        flow: &MixtureEuler3d,
        source: SpatialTransport<L>,
        normal_axis: usize,
        walls: [MolecularWall3d; 2],
        config: MixtureDynamicChannelThermalSubgridConfig,
        maximum_workspace_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, DissipationError> {
        poll(&mut cancelled)?;
        config.validate()?;
        admit_domain(flow, normal_axis)?;
        for wall in walls {
            admit_wall(wall, flow.grid.shape[normal_axis])?;
        }
        admit_workspace(
            Self::required_workspace_bytes_with_dynamic_channel_thermal_sgs(flow, &source)?,
            maximum_workspace_bytes,
        )?;
        let mut result = Self::new(flow, source, maximum_workspace_bytes, &mut cancelled)?;
        for (end, wall) in walls.into_iter().enumerate() {
            result = result.with_wall(BoundarySide::ALL[2 * normal_axis + end], wall)?;
        }
        result.sgs =
            Some(SubgridTransport {
                closure: ViscosityClosure::DynamicChannelThermal(
                    DynamicChannelThermalTransport::new(&flow.grid, normal_axis, config)?,
                ),
                dynamic_viscosity: allocated(flow.cell_count(), 0.0)?,
                heat_capacity_v_density: allocated(flow.cell_count(), 0.0)?,
            });
        let bytes = result
            .workspace_bytes()
            .checked_add(flow.workspace_bytes())
            .ok_or(DissipationError::Invalid(
                "dynamic mixture channel workspace overflow",
            ))?;
        admit_workspace(bytes, maximum_workspace_bytes)?;
        poll(&mut cancelled)?;
        Ok(result)
    }

    pub fn dynamic_channel_thermal_subgrid_config(
        &self,
    ) -> Option<MixtureDynamicChannelThermalSubgridConfig> {
        match &self.sgs.as_ref()?.closure {
            ViscosityClosure::DynamicChannelThermal(dynamic) => Some(dynamic.config),
            _ => None,
        }
    }
}

fn admit_domain(flow: &MixtureEuler3d, normal_axis: usize) -> Result<(), DissipationError> {
    if normal_axis >= 3 {
        return Err(DissipationError::Invalid(
            "mixture channel normal axis outside 0..3",
        ));
    }
    for axis in 0..3 {
        let expected = if axis == normal_axis {
            AxisBoundaries::Ends {
                lower: BoundaryCondition::Reflecting,
                upper: BoundaryCondition::Reflecting,
            }
        } else {
            AxisBoundaries::Periodic
        };
        let minimum = if axis == normal_axis { 2 } else { 3 };
        if flow.boundaries.axis(axis) != expected || flow.grid.shape[axis] < minimum {
            return Err(DissipationError::Invalid(
                "mixture channel requires one reflecting normal axis and two resolved periodic axes",
            ));
        }
    }
    Ok(())
}

fn admit_wall(wall: MolecularWall3d, normal_cells: usize) -> Result<(), DissipationError> {
    match wall {
        MolecularWall3d::SlipAdiabatic => Err(DissipationError::Invalid(
            "dynamic mixture channel requires two explicit no-slip walls",
        )),
        MolecularWall3d::NoSlip {
            thermal: WallThermalCondition::Adiabatic,
            ..
        } if normal_cells < 3 => Err(DissipationError::Invalid(
            "adiabatic dynamic mixture channel requires three normal fluid cells",
        )),
        MolecularWall3d::NoSlip { .. } => Ok(()),
    }
}

pub(super) struct DynamicChannelThermalTransport {
    pub(super) config: MixtureDynamicChannelThermalSubgridConfig,
    momentum: ChannelDynamicSgs,
    heat: ChannelDynamicScalar,
    density: Vec<f64>,
    velocity: Vec<[f64; 3]>,
    temperature: Vec<f64>,
    rho_diffusivity: Vec<f64>,
}

impl DynamicChannelThermalTransport {
    fn new(
        grid: &CartesianGrid,
        normal_axis: usize,
        config: MixtureDynamicChannelThermalSubgridConfig,
    ) -> Result<Self, DissipationError> {
        let momentum_bytes =
            ChannelDynamicSgs::required_bytes(grid.len()).map_err(dynamic_error)?;
        let heat_bytes = ChannelDynamicScalar::required_bytes(grid.len()).map_err(dynamic_error)?;
        Ok(Self {
            config,
            momentum: ChannelDynamicSgs::new_cartesian(
                grid.shape,
                grid.spacing,
                normal_axis,
                momentum_bytes,
            )
            .map_err(dynamic_error)?,
            heat: ChannelDynamicScalar::new(grid.shape, normal_axis, heat_bytes)
                .map_err(dynamic_error)?,
            density: allocated(grid.len(), 0.0)?,
            velocity: allocated(grid.len(), [0.0; 3])?,
            temperature: allocated(grid.len(), 0.0)?,
            rho_diffusivity: allocated(grid.len(), 0.0)?,
        })
    }

    pub(super) fn workspace_bytes(&self) -> usize {
        self.momentum.workspace_bytes()
            + self.heat.workspace_bytes()
            + (self.density.capacity()
                + self.temperature.capacity()
                + self.rho_diffusivity.capacity())
                * size_of::<f64>()
            + self.velocity.capacity() * size_of::<[f64; 3]>()
    }

    pub(super) fn admit_wall(&self, wall: MolecularWall3d) -> Result<(), DissipationError> {
        admit_wall(wall, self.momentum.shape()[self.momentum.normal_axis()])
    }

    /// True even for a pure heat fit, so the parent's original source Cv is
    /// refreshed when eddy stress and common-Schmidt species mixing are zero.
    pub(super) fn prepare(
        &mut self,
        primitive: &[Primitive],
        walls: &[MolecularWall3d; 6],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<bool, DissipationError> {
        let mut wall_velocity = [[0.0; 3]; 2];
        let mut scalar_walls = [ScalarWall::ZeroGradient; 2];
        for end in 0..2 {
            let wall = walls[2 * self.momentum.normal_axis() + end];
            self.admit_wall(wall)?;
            let MolecularWall3d::NoSlip { velocity, thermal } = wall else {
                unreachable!()
            };
            wall_velocity[end] = velocity;
            scalar_walls[end] = match thermal {
                WallThermalCondition::Adiabatic => ScalarWall::ZeroGradient,
                WallThermalCondition::Isothermal(value) => ScalarWall::Value(value.value()),
            };
        }
        for (i, p) in primitive.iter().enumerate() {
            if i % 64 == 0 {
                poll(cancelled)?;
            }
            self.density[i] = p.density;
            self.velocity[i] = p.velocity;
            self.temperature[i] = p.temperature;
        }
        let (momentum, heat) = self
            .momentum
            .prepare_with_scalar(
                self.config.momentum,
                self.config.heat,
                &self.density,
                &self.velocity,
                wall_velocity,
                &self.temperature,
                scalar_walls,
                &mut self.heat,
                &mut *cancelled,
            )
            .map_err(dynamic_error)?;
        for (i, (weight, &alpha)) in self
            .rho_diffusivity
            .iter_mut()
            .zip(self.heat.diffusivities())
            .enumerate()
        {
            if i % 64 == 0 {
                poll(cancelled)?;
            }
            *weight = self.density[i] * alpha;
            if !weight.is_finite() {
                return Err(DissipationError::Invalid(
                    "unrepresentable dynamic mixture channel heat coefficient",
                ));
            }
        }
        poll(cancelled)?;
        Ok(momentum.maximum_eddy_viscosity != 0.0 || heat.maximum_diffusivity != 0.0)
    }

    pub(super) fn viscosities(&self) -> &[f64] {
        self.momentum.viscosities()
    }

    pub(super) fn face_heat(&self, left: usize, right: usize) -> (f64, f64) {
        (
            0.5 * self.rho_diffusivity[left] + 0.5 * self.rho_diffusivity[right],
            self.config.heat.maximum_diffusivity_to_molecular,
        )
    }
}
