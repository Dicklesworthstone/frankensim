//! Explicit subgrid stress, heat and species transport on the original mixture.
//!
//! The Maxwell--Stefan flux and its Soret/Dufour terms remain source-owned.
//! This exploration closure adds `mu_t = rho nu_t`, deviatoric stress/work,
//! `k_t = mu_t cp/Pr_t`, and common-Schmidt mass-frame species fluxes
//! `j_i^t = -(mu_t/Sc_t) grad(Y_i)`. The same original referenced partial
//! enthalpies accompany those species fluxes in conservative total energy.
//! There is no separate viscous or chemical heating source.
//! Dynamic stress, heat and species apply their fitted coefficients to the
//! conservative face strain, with separately averaged density and coefficient.
//! A zero cell-centered strain cannot erase transport across resolved face shear.
//! Lagrangian stress, optionally with independent heat/common-species memories,
//! holds the fitted coefficients over one complete physical coupled step and
//! commits all signed pathline histories only with accepted flow.
//! The existing stage-local composition, calorics and face strain remain active.
//! An optional joint Germano fit supplies a common species diffusivity directly,
//! independently of stress and heat; it retains this same mass-frame and
//! referenced-enthalpy flux rather than prescribing a turbulent Schmidt number.
//!
//! Physical wall faces retain the original molecular law unless explicitly
//! given [`EquilibriumWallConfig`] for smooth adiabatic wall stress, or
//! [`EquilibriumThermalWallConfig`] for equilibrium isothermal stress/heat.
//! Other walls require resolved molecular layers; no wall species or turbulent
//! reaction-rate closure is inferred. The existing
//! reactive split integrator can compose this operator without changing its
//! chemistry owner. Isotropic SGS energy, backscatter, pressure dilatation,
//! flame wrinkling and shock/turbulence interaction are not represented.
//!
//! Complete prescribed reservoirs retain molecular exchange by default.
//! [`MolecularTransport3d::with_local_sgs_reservoir_exchange`] explicitly adds
//! local modeled stress, heat and common-species mixing using the same physical
//! bath traces and source energy effort. No exterior dynamic fit is inferred.
//!
//! The heat/species energy split follows CFX Theory Guide 2025 R1, §1.5.7,
//! equations (1-165)--(1-173):
//! https://ansyshelp.ansys.com/public/Views/Secured/corp/v251/en/cfx_thry/BGBHJAGD.html
//! The shared velocity closures retain their original Nicoud--Ducros and
//! Vreman references in `fs-lbm::sgs3`.

use super::super::sum;
use super::*;
use crate::compressible::spatial::molecular::sgs::face_gradient;
use crate::compressible::spatial::molecular::sgs::wall::EquilibriumWallModel;
pub use crate::compressible::spatial::molecular::sgs::{
    ChannelDynamicThermalSubgridConfig, DynamicAveraging, DynamicScalarConfig, DynamicScalarReport,
    DynamicSgsConfig, DynamicSgsReport, DynamicSubgridConfig, DynamicThermalSubgridConfig,
    EquilibriumThermalWallConfig, EquilibriumWallConfig, LagrangianSgsConfig,
    LagrangianSubgridConfig, SgsModel, SgsUse, SubgridConfig, WallHeatConfig,
};
use fs_gas::multicomponent::dissipation::{MaxwellStefanSolveV1, TransportFaceReport};
use fs_gas::multicomponent::{Error as FlowError, MAX_SPECIES};
use fs_lbm::sgs3::{eddy_viscosity, smagorinsky_viscosity, vreman_eddy_viscosity};
use fs_qty::Temperature;
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

mod dynamic;
mod flux;
pub(in crate::compressible::spatial::mixture) mod lagrangian;
mod reservoir;
pub use lagrangian::{
    LagrangianComponentsConfig, LagrangianScalarConfig, MixtureLagrangianSpeciesSubgridConfig,
    MixtureLagrangianSubgridConfig,
};
mod wall;
pub use dynamic::{
    ChannelDynamicComponentsConfig, DynamicComponentsConfig,
    MixtureChannelDynamicSpeciesSubgridConfig, MixtureDynamicSpeciesSubgridConfig,
    MixtureDynamicSubgridConfig, MixtureDynamicThermalSubgridConfig,
};
mod channel;
pub use channel::{
    ChannelDynamicScalarConfig, ChannelDynamicSgsConfig, MixtureDynamicChannelThermalSubgridConfig,
};

/// Explicit stress/heat controls shared with frozen gas, plus species mixing.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MixtureSubgridConfig {
    pub transport: SubgridConfig,
    /// Positive finite Schmidt number shared by all turbulent species fluxes.
    /// Molecular diffusion retains its distinct original binary coefficients.
    pub turbulent_schmidt: f64,
}

impl MixtureSubgridConfig {
    pub fn validate(self) -> Result<(), DissipationError> {
        self.transport.validate().map_err(|_| {
            DissipationError::Invalid("invalid mixture subgrid stress/heat controls")
        })?;
        if !positive(self.turbulent_schmidt) {
            return Err(DissipationError::Invalid(
                "positive finite turbulent Schmidt number required",
            ));
        }
        Ok(())
    }
}

impl<L: TransportLaw> MolecularTransport3d<L> {
    /// Combined numerical bytes including the prepared SGS viscosity/Cv cache.
    pub fn required_workspace_bytes_with_sgs(
        flow: &MixtureEuler3d,
        source: &SpatialTransport<L>,
    ) -> Result<usize, DissipationError> {
        Self::required_workspace_bytes(flow, source)?
            .checked_add(flow.cell_count().checked_mul(2 * size_of::<f64>()).ok_or(
                DissipationError::Invalid("mixture subgrid workspace overflow"),
            )?)
            .ok_or(DissipationError::Invalid(
                "mixture subgrid workspace overflow",
            ))
    }

    /// Prepare explicitly authored SGS on the existing molecular/mixture owner.
    /// Zero model coefficient preserves the molecular field and reports exactly.
    pub fn with_sgs(
        flow: &MixtureEuler3d,
        source: SpatialTransport<L>,
        config: MixtureSubgridConfig,
        maximum_workspace_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, DissipationError> {
        poll(&mut cancelled)?;
        config.validate()?;
        admit_workspace(
            Self::required_workspace_bytes_with_sgs(flow, &source)?,
            maximum_workspace_bytes,
        )?;
        let mut result = Self::new(flow, source, maximum_workspace_bytes, &mut cancelled)?;
        result.sgs = Some(SubgridTransport::new(&flow.grid, config)?);
        admit_workspace(
            result.workspace_bytes() + flow.workspace_bytes(),
            maximum_workspace_bytes,
        )?;
        poll(&mut cancelled)?;
        Ok(result)
    }

    pub fn subgrid_config(&self) -> Option<MixtureSubgridConfig> {
        match &self.sgs.as_ref()?.closure {
            ViscosityClosure::Local { config, .. } => Some(*config),
            ViscosityClosure::Dynamic(_)
            | ViscosityClosure::Lagrangian(_)
            | ViscosityClosure::DynamicThermal(_)
            | ViscosityClosure::DynamicChannelThermal(_)
            | ViscosityClosure::DynamicSpecies(_)
            | ViscosityClosure::ChannelDynamicSpecies(_) => None,
        }
    }
}

enum ViscosityClosure {
    Local {
        config: MixtureSubgridConfig,
        widths: [f64; 3],
        filter_width: f64,
        reservoir_exchange: bool,
    },
    Dynamic(dynamic::DynamicTransport),
    Lagrangian(lagrangian::LagrangianTransport),
    DynamicThermal(dynamic::DynamicThermalTransport),
    DynamicChannelThermal(channel::DynamicChannelThermalTransport),
    DynamicSpecies(dynamic::DynamicSpeciesTransport),
    ChannelDynamicSpecies(dynamic::ChannelDynamicSpeciesTransport),
}

pub(super) struct SubgridTransport {
    closure: ViscosityClosure,
    wall_models: [Option<EquilibriumWallModel>; 6],
    dynamic_coefficient_squared: Vec<f64>,
    heat_capacity_v_density: Vec<f64>,
}

impl SubgridTransport {
    fn new(grid: &CartesianGrid, config: MixtureSubgridConfig) -> Result<Self, DissipationError> {
        let widths = grid
            .spacing
            .map(|h| h * config.transport.filter_width_factor);
        let filter_width = grid.cell_volume().cbrt() * config.transport.filter_width_factor;
        if !widths.into_iter().chain([filter_width]).all(positive) {
            return Err(DissipationError::Invalid(
                "unrepresentable mixture subgrid widths",
            ));
        }
        Ok(Self {
            wall_models: [None; 6],
            closure: ViscosityClosure::Local {
                config,
                widths,
                filter_width,
                reservoir_exchange: false,
            },
            dynamic_coefficient_squared: allocated(grid.len(), 0.0)?,
            heat_capacity_v_density: allocated(grid.len(), 0.0)?,
        })
    }

    pub(super) fn workspace_bytes(&self) -> usize {
        (self.dynamic_coefficient_squared.capacity() + self.heat_capacity_v_density.capacity())
            * size_of::<f64>()
            + match &self.closure {
                ViscosityClosure::Local { .. } => 0,
                ViscosityClosure::Dynamic(dynamic) => dynamic.workspace_bytes(),
                ViscosityClosure::Lagrangian(dynamic) => dynamic.workspace_bytes(),
                ViscosityClosure::DynamicThermal(dynamic) => dynamic.workspace_bytes(),
                ViscosityClosure::DynamicChannelThermal(dynamic) => dynamic.workspace_bytes(),
                ViscosityClosure::DynamicSpecies(dynamic) => dynamic.workspace_bytes(),
                ViscosityClosure::ChannelDynamicSpecies(dynamic) => dynamic.workspace_bytes(),
            }
    }

    pub(super) fn admit_wall(
        &self,
        side: BoundarySide,
        wall: MolecularWall3d,
    ) -> Result<(), DissipationError> {
        if let Some(config) = self.wall_models[side_index(side)] {
            config.admit_wall(wall).map_err(DissipationError::Invalid)?;
        }
        match &self.closure {
            ViscosityClosure::Lagrangian(dynamic) => dynamic.admit_wall(side, wall),
            ViscosityClosure::DynamicChannelThermal(dynamic) => dynamic.admit_wall(wall),
            ViscosityClosure::ChannelDynamicSpecies(_)
                if !matches!(wall, MolecularWall3d::NoSlip { .. }) =>
            {
                Err(DissipationError::Invalid(
                    "channel SGS requires no-slip wall values",
                ))
            }
            _ => Ok(()),
        }
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn prepare(
        &mut self,
        gas: &MixtureGas,
        states: &[f64],
        primitive: &[Primitive],
        _gradient: &[[[f64; 3]; 3]],
        walls: &[MolecularWall3d; 6],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), DissipationError> {
        match &mut self.closure {
            ViscosityClosure::Local { config, .. } => {
                if config.transport.model.coefficient() == 0.0 {
                    return Ok(());
                }
            }
            ViscosityClosure::Dynamic(dynamic) => {
                if dynamic.prepare(primitive, cancelled)? == 0.0 {
                    // A fit can become zero after a nonzero earlier SSP stage.
                    // Preserve molecular face arithmetic without stale stress.
                    for values in self.dynamic_coefficient_squared.chunks_mut(64) {
                        poll(cancelled)?;
                        values.fill(0.0);
                    }
                    return Ok(());
                }
            }
            ViscosityClosure::Lagrangian(dynamic) => {
                if !dynamic.admit_prepared_step()? {
                    return Ok(());
                }
            }
            ViscosityClosure::DynamicThermal(dynamic) => {
                if !dynamic.prepare(primitive, cancelled)? {
                    for values in self.dynamic_coefficient_squared.chunks_mut(64) {
                        poll(cancelled)?;
                        values.fill(0.0);
                    }
                    return Ok(());
                }
            }
            ViscosityClosure::DynamicSpecies(dynamic) => {
                if !dynamic.prepare(states, primitive, cancelled)? {
                    for values in self.dynamic_coefficient_squared.chunks_mut(64) {
                        poll(cancelled)?;
                        values.fill(0.0);
                    }
                    return Ok(());
                }
            }
            ViscosityClosure::DynamicChannelThermal(dynamic) => {
                if !dynamic.prepare(primitive, walls, cancelled)? {
                    for values in self.dynamic_coefficient_squared.chunks_mut(64) {
                        poll(cancelled)?;
                        values.fill(0.0);
                    }
                    return Ok(());
                }
            }
            ViscosityClosure::ChannelDynamicSpecies(dynamic) => {
                if !dynamic.prepare(states, primitive, walls, cancelled)? {
                    for values in self.dynamic_coefficient_squared.chunks_mut(64) {
                        poll(cancelled)?;
                        values.fill(0.0);
                    }
                    return Ok(());
                }
            }
        }
        let n = gas.species_count();
        for (i, (p, state)) in primitive.iter().zip(states.chunks_exact(n + 4)).enumerate() {
            poll(cancelled)?;
            let coefficient = match &self.closure {
                // Local closures supply their authored coefficient at the face.
                ViscosityClosure::Local { .. } => Some(0.0),
                ViscosityClosure::Dynamic(dynamic) => dynamic.coefficient_squared(i),
                ViscosityClosure::Lagrangian(_) => Some(self.dynamic_coefficient_squared[i]),
                ViscosityClosure::DynamicThermal(dynamic) => dynamic.coefficient_squared(i),
                ViscosityClosure::DynamicChannelThermal(dynamic) => dynamic.coefficient_squared(i),
                ViscosityClosure::DynamicSpecies(dynamic) => dynamic.coefficient_squared(i),
                ViscosityClosure::ChannelDynamicSpecies(dynamic) => dynamic.coefficient_squared(i),
            }
            .ok_or(DissipationError::Invalid(
                "dynamic mixture coefficient cache length",
            ))?;
            let mut cv = [0.0; MAX_SPECIES];
            for (k, model) in gas.source().models().iter().enumerate() {
                if k % 16 == 0 {
                    poll(cancelled)?;
                }
                if state[k] != 0.0 {
                    cv[k] = (state[k] / model.molar_mass().value())
                        * (model
                            .evaluate_properties(Temperature::new(p.temperature))
                            .map_err(FlowError::from)?
                            .cp()
                            .value()
                            - R);
                }
            }
            let cv = sum(cv[..n].iter().copied());
            if !positive(cv) {
                return Err(DissipationError::Invalid(
                    "unrepresentable mixture subgrid coefficients",
                ));
            }
            self.dynamic_coefficient_squared[i] = coefficient;
            self.heat_capacity_v_density[i] = cv;
        }
        Ok(())
    }

    /// Add to the source's normal/tangent flux before global-axis restoration.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn add_face(
        &self,
        gas: &MixtureGas,
        states: &[f64],
        primitive: &[Primitive],
        gradient: &[[[f64; 3]; 3]],
        left: usize,
        right: usize,
        axis: usize,
        distance: f64,
        receipt: TransportFaceReport,
        controls: MaxwellStefanSolveV1,
        out: &mut [f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<TransportFaceReport, DissipationError> {
        let dynamic_heat = match &self.closure {
            ViscosityClosure::Lagrangian(dynamic) => dynamic.face_heat_coefficient(left, right),
            ViscosityClosure::DynamicThermal(dynamic) => {
                Some(dynamic.face_heat_coefficient(left, right))
            }
            ViscosityClosure::DynamicChannelThermal(dynamic) => {
                Some(dynamic.face_heat_coefficient(left, right))
            }
            ViscosityClosure::DynamicSpecies(dynamic) => {
                Some(dynamic.face_heat_coefficient(left, right))
            }
            ViscosityClosure::ChannelDynamicSpecies(dynamic) => {
                Some(dynamic.face_heat_coefficient(left, right))
            }
            _ => None,
        };
        let dynamic_species = match &self.closure {
            ViscosityClosure::Lagrangian(dynamic) => dynamic.face_species_coefficient(left, right),
            ViscosityClosure::DynamicSpecies(dynamic) => Some((
                dynamic.face_species_coefficient(left, right),
                dynamic.config.maximum_species_diffusivity_m2_s,
            )),
            ViscosityClosure::ChannelDynamicSpecies(dynamic) => Some((
                dynamic.face_species_coefficient(left, right),
                dynamic.config.maximum_species_diffusivity_m2_s,
            )),
            _ => None,
        };
        let momentum_coefficient = 0.5 * self.dynamic_coefficient_squared[left]
            + 0.5 * self.dynamic_coefficient_squared[right];
        let momentum_disabled = match &self.closure {
            ViscosityClosure::Local { config, .. } => config.transport.model.coefficient() == 0.0,
            _ => momentum_coefficient == 0.0,
        };
        // An exactly disabled model makes no extra caloric calls/arithmetic.
        if momentum_disabled
            && dynamic_heat.is_none_or(|(coefficient, _)| coefficient == 0.0)
            && dynamic_species.is_none_or(|(coefficient, _)| coefficient == 0.0)
        {
            return Ok(receipt);
        }
        // Reconstruct once at the conservative face for all three transports.
        // This owner's cell tensor is derivative-first; transpose to du_i/dx_j.
        let transpose =
            |g: [[f64; 3]; 3]| std::array::from_fn(|i| std::array::from_fn(|j| g[j][i]));
        let g = face_gradient(
            primitive[left].velocity,
            primitive[right].velocity,
            transpose(gradient[left]),
            transpose(gradient[right]),
            axis,
            distance,
        );
        let filter_width = match &self.closure {
            ViscosityClosure::Local { filter_width, .. } => *filter_width,
            ViscosityClosure::Dynamic(dynamic) => dynamic.filter_width(),
            ViscosityClosure::Lagrangian(dynamic) => dynamic.filter_width(),
            ViscosityClosure::DynamicThermal(dynamic) => dynamic.filter_width(),
            ViscosityClosure::DynamicChannelThermal(dynamic) => dynamic.filter_width(),
            ViscosityClosure::DynamicSpecies(dynamic) => dynamic.filter_width(),
            ViscosityClosure::ChannelDynamicSpecies(dynamic) => dynamic.filter_width(),
        };
        let density = 0.5 * primitive[left].density + 0.5 * primitive[right].density;
        // Dynamic momentum, heat and species share |S|. Reuse the normalized
        // invariant once while retaining each coefficient's original length
        // scaling; forming Delta² first would narrow the representable range.
        let shear_rate = if matches!(&self.closure, ViscosityClosure::Local { .. }) {
            0.0
        } else {
            eddy_viscosity(SgsModel::Smagorinsky { cs: 1.0 }, 1.0, g)
                .map_err(|_| DissipationError::Invalid("dynamic mixture face strain or diffusivity"))?
        };
        let face_diffusivity = |coefficient: f64| -> Result<f64, DissipationError> {
            if coefficient == 0.0 {
                return Ok(0.0);
            }
            smagorinsky_viscosity(
                fs_math::det::sqrt(coefficient),
                filter_width,
                shear_rate,
            )
            .map_err(|_| DissipationError::Invalid("dynamic mixture face strain or diffusivity"))
        };
        let mu = density
            * if let ViscosityClosure::Local { config, widths, .. } = &self.closure {
                match config.transport.model {
                    SgsModel::Vreman { cv } => vreman_eddy_viscosity(cv, *widths, g),
                    model => eddy_viscosity(model, filter_width, g),
                }
                .map_err(|_| {
                    DissipationError::Invalid("mixture subgrid face gradient or viscosity")
                })?
            } else {
                face_diffusivity(momentum_coefficient)?
            };
        let dynamic_heat = dynamic_heat
            .map(|(coefficient, maximum)| {
                face_diffusivity(coefficient).map(|alpha| (density * alpha, maximum))
            })
            .transpose()?;
        let dynamic_species = dynamic_species
            .map(|(coefficient, maximum)| {
                let diffusivity = face_diffusivity(coefficient)?;
                if diffusivity > maximum {
                    return Err(DissipationError::Invalid(
                        "dynamic species face diffusivity exceeds its physical admission bound",
                    ));
                }
                Ok(density * diffusivity)
            })
            .transpose()?;
        if mu == 0.0
            && dynamic_heat.is_none_or(|(weight, _)| weight == 0.0)
            && dynamic_species.is_none_or(|weight| weight == 0.0)
        {
            return Ok(receipt);
        }
        poll(cancelled)?;
        let (prandtl, schmidt, maximum_ratio) = match &self.closure {
            ViscosityClosure::Local { config, .. } => (
                Some(config.transport.turbulent_prandtl),
                Some(config.turbulent_schmidt),
                config.transport.maximum_eddy_to_molecular,
            ),
            ViscosityClosure::Dynamic(dynamic) => (
                Some(dynamic.config.transport.turbulent_prandtl),
                Some(dynamic.config.turbulent_schmidt),
                dynamic.config.transport.model.maximum_eddy_to_molecular,
            ),
            ViscosityClosure::Lagrangian(dynamic) => dynamic.transport_controls(),
            ViscosityClosure::DynamicThermal(dynamic) => (
                None,
                Some(dynamic.config.turbulent_schmidt),
                dynamic.config.transport.momentum.maximum_eddy_to_molecular,
            ),
            ViscosityClosure::DynamicSpecies(dynamic) => (
                None,
                None,
                dynamic.config.transport.momentum.maximum_eddy_to_molecular,
            ),
            ViscosityClosure::DynamicChannelThermal(dynamic) => (
                None,
                Some(dynamic.config.turbulent_schmidt),
                dynamic.config.momentum.maximum_eddy_to_molecular,
            ),
            ViscosityClosure::ChannelDynamicSpecies(dynamic) => (
                None,
                None,
                dynamic.config.transport.momentum.maximum_eddy_to_molecular,
            ),
        };
        if mu != 0.0 {
            let ratio = mu / receipt.dynamic_viscosity_pa_s;
            if !ratio.is_finite() || ratio > maximum_ratio {
                return Err(DissipationError::Invalid(
                    "mixture subgrid eddy/molecular viscosity ratio exceeded",
                ));
            }
        }
        let n = gas.species_count();
        let width = n + 4;
        flux::add(
            gas,
            flux::Face {
                partial_densities: [
                    &states[left * width..left * width + n],
                    &states[right * width..right * width + n],
                ],
                primitive: [primitive[left], primitive[right]],
                cv_density: [
                    self.heat_capacity_v_density[left],
                    self.heat_capacity_v_density[right],
                ],
                gradient: g,
                axis,
                distance,
                work_velocity: std::array::from_fn(|d| {
                    0.5 * primitive[left].velocity[d] + 0.5 * primitive[right].velocity[d]
                }),
            },
            mu,
            dynamic_heat,
            dynamic_species,
            prandtl,
            schmidt,
            receipt,
            controls,
            out,
            cancelled,
        )
    }
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}

#[cfg(test)]
mod tests;
