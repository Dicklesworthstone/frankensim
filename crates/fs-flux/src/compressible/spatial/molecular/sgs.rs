//! Local equilibrium subgrid stress and heat transport on the conservative gas.
//!
//! The resolved velocity is momentum/density (the Favre velocity). Each SSP
//! stage recomputes the shared LBM Smagorinsky/WALE/Vreman closure from the full
//! physical velocity gradient at each conservative face. The normal derivative
//! uses the adjacent velocities, exactly as the stress does; transverse
//! derivatives use the two cell traces. Cell-centered central differences alone
//! miss alternating-cell shear and must not select the local face viscosity.
//! Dynamic fits also retain C_s squared independently of cell strain and
//! evaluate `mu_t = rho_face C_s_face^2 Delta^2 |S_face|` on this same face
//! tensor. Density and the nonnegative coefficient are interpolated separately;
//! a zero centered strain therefore cannot erase a resolved face shear.
//! The opt-in Lagrangian fit samples each accepted physical step's start state,
//! holds its coefficient through the complete split, and commits its history
//! only with that step. Its pathline/coefficient time integration is first order;
//! internal SSP/RKL stages and coupling retries never age the memory.
//! Face stress adds `mu_t = rho nu_t` only to the
//! deviatoric viscosity; bulk viscosity remains the original molecular law.
//! Gradient-diffusion heat adds `k_t = mu_t cp(T)/Pr_t`, using the original
//! caloric source at the face. Stress work and heat use the existing shared
//! total-energy flux, so resolved kinetic loss is not counted as heat twice.
//!
//! This is an explicit exploration model for low-Mach, frozen-composition
//! flow, not a closure for shock/turbulence interaction. Isotropic SGS kinetic
//! energy, pressure dilatation, backscatter and species mixing are absent.
//! Reflecting boundaries retain their adiabatic slip semantics unless an
//! explicit wall is authored. [`EquilibriumWallConfig`] optionally supplies
//! smooth adiabatic wall stress when the viscous sublayer is not resolved.
//! [`EquilibriumThermalWallConfig`] explicitly adds the integrated thermal wall
//! resistance for isothermal no-slip walls, using the original face calorics.
//! Complete reservoirs remain molecular-only unless local SGS bath exchange is
//! explicitly enabled. That policy uses the admitted bath trace, half-cell
//! normal differences and the same tangential tensor as the bath stress. It
//! supplies no exterior Germano fit, Lagrangian history or inflow fluctuations.
//! No coefficient or turbulent Prandtl number is chosen implicitly.
//!
//! References: Nicoud & Ducros (1999), doi:10.1023/A:1009995426001;
//! Vreman (2004), doi:10.1063/1.1785131. The stress and heat closure follow
//! Nicoud, "Unsteady flows modeling and computation" (2007), equations
//! (24)-(26): https://imag.umontpellier.fr/~nicoud/PDF/notes_Nicoud_1and2.pdf.

use super::*;
use fs_lbm::sgs3::{SgsConfig, eddy_viscosity, smagorinsky_viscosity, vreman_eddy_viscosity};
pub use fs_lbm::sgs3::{SgsModel, SgsUse};
pub(in crate::compressible::spatial) mod wall;
pub use wall::{EquilibriumThermalWallConfig, EquilibriumWallConfig, WallHeatConfig};
use wall::EquilibriumWallModel;

mod dynamic;
pub use dynamic::DynamicSubgridConfig;
mod lagrangian;
pub use fs_lbm::sgs3::dynamic::lagrangian::LagrangianSgsConfig;
pub use fs_lbm::sgs3::dynamic::lagrangian::scalar::LagrangianScalarConfig;
pub use fs_lbm::sgs3::dynamic::scalar::{DynamicScalarConfig, DynamicScalarReport};
pub use fs_lbm::sgs3::dynamic::{DynamicAveraging, DynamicSgsConfig, DynamicSgsReport};
pub use lagrangian::{LagrangianSubgridConfig, LagrangianThermalSubgridConfig};
mod thermal;
pub use thermal::DynamicThermalSubgridConfig;
mod channel;
pub use channel::ChannelDynamicThermalSubgridConfig;
pub(in crate::compressible::spatial) use channel::{channel_axis, channel_wall_values};
pub use fs_lbm::sgs3::dynamic::channel::{ChannelDynamicScalarConfig, ChannelDynamicSgsConfig};

/// Match the local constitutive tensor to the conservative stress stencil.
/// The input and output tensors are component-first: `g[i][j] = du_i/dx_j`.
pub(in crate::compressible::spatial) fn face_gradient(
    left_velocity: [f64; 3],
    right_velocity: [f64; 3],
    left_gradient: [[f64; 3]; 3],
    right_gradient: [[f64; 3]; 3],
    normal_axis: usize,
    distance: f64,
) -> [[f64; 3]; 3] {
    std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            if j == normal_axis {
                (right_velocity[i] - left_velocity[i]) / distance
            } else {
                0.5 * left_gradient[i][j] + 0.5 * right_gradient[i][j]
            }
        })
    })
}

/// Convert the bath stress stencil to a full component-first velocity tensor.
/// Tangential input rows follow the two derivative axes; their columns are
/// velocity components in normal-frame order `[normal, tangent1, tangent2]`.
pub(in crate::compressible::spatial) fn reservoir_face_gradient(
    left_velocity: [f64; 3],
    right_velocity: [f64; 3],
    tangential_gradient: [[f64; 3]; 2],
    normal_axis: usize,
    distance: f64,
) -> [[f64; 3]; 3] {
    let [a, b] = tangents(normal_axis);
    let mut result = [[0.0; 3]; 3];
    for (column, component) in [normal_axis, a, b].into_iter().enumerate() {
        result[component][normal_axis] =
            (right_velocity[component] - left_velocity[component]) / distance;
        result[component][a] = tangential_gradient[0][column];
        result[component][b] = tangential_gradient[1][column];
    }
    result
}

fn local_eddy_viscosity(
    model: SgsModel,
    widths: [f64; 3],
    filter_width: f64,
    gradient: [[f64; 3]; 3],
) -> Result<f64, MolecularError> {
    match model {
        SgsModel::Vreman { cv } => vreman_eddy_viscosity(cv, widths, gradient),
        model => eddy_viscosity(model, filter_width, gradient),
    }
    .map_err(MolecularError::Subgrid)
}

/// Explicit subgrid controls; enabling the Cargo feature alone changes no flow.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SubgridConfig {
    pub model: SgsModel,
    /// Positive multiplier on the grid widths. Smagorinsky/WALE use the
    /// geometric mean; Vreman retains all three directional widths.
    pub filter_width_factor: f64,
    /// Positive, finite gradient-diffusion turbulent Prandtl number.
    pub turbulent_prandtl: f64,
    /// Maximum face `mu_t/mu`, refused rather than clipped.
    pub maximum_eddy_to_molecular: f64,
    /// Acceptance renders refuse, as in the originating closure owner.
    pub usage: SgsUse,
}

impl SubgridConfig {
    pub fn validate(self) -> Result<(), MolecularError> {
        SgsConfig {
            model: self.model,
            filter_width_lu: self.filter_width_factor,
            maximum_eddy_to_molecular: self.maximum_eddy_to_molecular,
            usage: self.usage,
        }
        .validate()
        .map_err(MolecularError::Subgrid)?;
        if !positive(self.turbulent_prandtl) {
            return Err(MolecularError::Invalid(
                "positive finite turbulent Prandtl number required",
            ));
        }
        Ok(())
    }
}

impl<L: TransportLaw> MolecularTransport3d<L> {
    /// Prepare the same molecular operator with explicitly authored SGS.
    /// The zero-coefficient control is bit-identical to molecular-only flow.
    pub fn with_sgs(
        flow: &Euler3d,
        law: L,
        config: SubgridConfig,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, MolecularError> {
        config.validate()?;
        let mut result = Self::new(flow, law, &mut cancelled)?;
        result.sgs = Some(SubgridTransport::new(flow, config)?);
        poll(&mut cancelled)?;
        Ok(result)
    }

    /// Extend explicitly authored local SGS stress and heat to complete baths.
    ///
    /// This also enables molecular reservoir exchange. Source-admitted bath
    /// traces supply the half-cell normal and profile tangential gradients;
    /// density and caloric properties use the same cell/bath midpoint as the
    /// molecular coefficients. Conservative work uses the prescribed bath
    /// velocity. The authored viscosity-ratio limit and stage diffusion bound
    /// include the modeled bath transport.
    ///
    /// Requires local Smagorinsky, WALE or Vreman SGS. Dynamic and Lagrangian
    /// fits refuse: no exterior filter samples or histories are fabricated.
    /// Unbound sides and pressure-only outlets retain zero diffusive flux.
    pub fn with_local_sgs_reservoir_exchange(mut self) -> Result<Self, MolecularError> {
        let Some(SubgridTransport {
            closure: SubgridClosure::Local {
                reservoir_exchange, ..
            },
            ..
        }) = self.sgs.as_mut()
        else {
            return Err(MolecularError::Invalid(
                "reservoir SGS exchange requires an authored local closure",
            ));
        };
        *reservoir_exchange = true;
        self.reservoir_exchange = true;
        Ok(self)
    }

    /// Whether complete reservoirs use the explicitly selected local SGS law.
    pub fn local_sgs_reservoir_exchange_enabled(&self) -> bool {
        self.sgs
            .as_ref()
            .is_some_and(SubgridTransport::local_reservoir_exchange_enabled)
    }

    /// Exact local controls; `None` for molecular-only or dynamic SGS flow.
    pub fn subgrid_config(&self) -> Option<SubgridConfig> {
        self.sgs.as_ref().and_then(|s| match &s.closure {
            SubgridClosure::Local { config, .. } => Some(*config),
            SubgridClosure::Dynamic(_)
            | SubgridClosure::Lagrangian(_)
            | SubgridClosure::DynamicThermal(_)
            | SubgridClosure::ChannelDynamicThermal(_) => None,
        })
    }
}

pub(super) struct SubgridTransport {
    closure: SubgridClosure,
    wall_models: [Option<EquilibriumWallModel>; 6],
    // Reused at every stage. Store C_s squared independently of cell strain,
    // so the conservative face selects its own viscosity. No hot allocation.
    dynamic_coefficient_squared: Vec<f64>,
}

enum SubgridClosure {
    Local {
        config: SubgridConfig,
        widths: [f64; 3],
        filter_width: f64,
        reservoir_exchange: bool,
    },
    Dynamic(dynamic::DynamicTransport),
    Lagrangian(lagrangian::LagrangianTransport),
    DynamicThermal(thermal::DynamicThermalTransport),
    ChannelDynamicThermal(channel::ChannelDynamicThermalTransport),
}

impl SubgridTransport {
    /// Owned capacities of the closure and of the fitted-coefficient cache; the
    /// six wall models are fixed-size fields.
    pub(super) fn workspace_bytes(&self) -> usize {
        self.dynamic_coefficient_squared.capacity() * core::mem::size_of::<f64>()
            + match &self.closure {
                SubgridClosure::Local { .. } => 0,
                SubgridClosure::Dynamic(dynamic) => dynamic.workspace_bytes(),
                SubgridClosure::Lagrangian(dynamic) => dynamic.workspace_bytes(),
                SubgridClosure::DynamicThermal(dynamic) => dynamic.workspace_bytes(),
                SubgridClosure::ChannelDynamicThermal(dynamic) => dynamic.workspace_bytes(),
            }
    }

    fn local_reservoir_exchange_enabled(&self) -> bool {
        matches!(
            &self.closure,
            SubgridClosure::Local {
                reservoir_exchange: true,
                ..
            }
        )
    }

    pub(super) fn requires_no_slip_wall(&self) -> bool {
        match &self.closure {
            SubgridClosure::ChannelDynamicThermal(_) => true,
            SubgridClosure::Lagrangian(dynamic) => dynamic.normal_axis().is_some(),
            _ => false,
        }
    }

    pub(super) fn new(flow: &Euler3d, config: SubgridConfig) -> Result<Self, MolecularError> {
        let widths = flow.grid.spacing.map(|h| h * config.filter_width_factor);
        let filter_width = flow.grid.cell_volume().cbrt() * config.filter_width_factor;
        if !widths.into_iter().chain([filter_width]).all(positive) {
            return Err(MolecularError::Invalid(
                "unrepresentable subgrid filter widths",
            ));
        }
        Ok(Self {
            wall_models: [None; 6],
            closure: SubgridClosure::Local {
                config,
                widths,
                filter_width,
                reservoir_exchange: false,
            },
            dynamic_coefficient_squared: allocated(flow.states.len(), 0.0)?,
        })
    }

    pub(super) fn prepare(
        &mut self,
        primitive: &[Primitive],
        walls: &[MolecularWall3d; 6],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), MolecularError> {
        match &mut self.closure {
            SubgridClosure::Local { .. } => Ok(()),
            SubgridClosure::Dynamic(dynamic) => {
                dynamic.prepare(primitive, &mut self.dynamic_coefficient_squared, cancelled)
            }
            SubgridClosure::Lagrangian(dynamic) => dynamic.admit_prepared_step(),
            SubgridClosure::DynamicThermal(dynamic) => {
                dynamic.prepare(primitive, &mut self.dynamic_coefficient_squared, cancelled)
            }
            SubgridClosure::ChannelDynamicThermal(dynamic) => dynamic.prepare(
                primitive,
                walls,
                &mut self.dynamic_coefficient_squared,
                cancelled,
            ),
        }
    }

    /// Local bath constitutive law on the exact conservative stress stencil.
    /// Disabled policies return before touching the caloric source or closure.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn reservoir_face_coefficients(
        &self,
        gas: &FrozenGas,
        axis: usize,
        face: usize,
        left: Primitive,
        right: Primitive,
        distance: f64,
        tangential_gradient: [[f64; 3]; 2],
        temperature: f64,
        pressure: f64,
        molecular_mu: f64,
        molecular_conductivity: f64,
    ) -> Result<(f64, f64, f64), MolecularError> {
        let SubgridClosure::Local {
            config,
            widths,
            filter_width,
            reservoir_exchange: true,
        } = &self.closure
        else {
            return Ok((molecular_mu, molecular_conductivity, 0.0));
        };
        if config.model.coefficient() == 0.0 {
            return Ok((molecular_mu, molecular_conductivity, 0.0));
        }
        let gradient = reservoir_face_gradient(
            left.velocity,
            right.velocity,
            tangential_gradient,
            axis,
            distance,
        );
        let nu = local_eddy_viscosity(config.model, *widths, *filter_width, gradient)?;
        let mu_t = (0.5 * left.density + 0.5 * right.density) * nu;
        add_eddy_transport(
            gas,
            axis,
            face,
            temperature,
            pressure,
            mu_t,
            molecular_mu,
            molecular_conductivity,
            config.turbulent_prandtl,
            config.maximum_eddy_to_molecular,
        )
    }

    /// Total viscosity/conductivity plus the separately retained modeled
    /// viscosity. Molecular face quadrature may change without subtracting a
    /// small modeled stress back out of a much larger molecular coefficient.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn face_coefficients(
        &self,
        gas: &FrozenGas,
        left: usize,
        right: usize,
        axis: usize,
        face: usize,
        distance: f64,
        primitive: &[Primitive],
        gradient: &[[[f64; 3]; 3]],
        temperature: f64,
        pressure: f64,
        molecular_mu: f64,
        molecular_conductivity: f64,
    ) -> Result<(f64, f64, f64), MolecularError> {
        if let SubgridClosure::Lagrangian(dynamic) = &self.closure {
            if let Some(config) = dynamic.thermal_config() {
                return dynamic.thermal_face_coefficients(
                    config,
                    gas,
                    left,
                    right,
                    axis,
                    face,
                    distance,
                    &self.dynamic_coefficient_squared,
                    primitive,
                    gradient,
                    temperature,
                    pressure,
                    molecular_mu,
                    molecular_conductivity,
                );
            }
        }
        let (mu_t, rho_alpha) = if let SubgridClosure::Local {
            config,
            widths,
            filter_width,
            ..
        } = &self.closure
        {
            if config.model.coefficient() == 0.0 {
                return Ok((molecular_mu, molecular_conductivity, 0.0));
            }
            let g = face_gradient(
                primitive[left].velocity,
                primitive[right].velocity,
                gradient[left],
                gradient[right],
                axis,
                distance,
            );
            let nu = local_eddy_viscosity(config.model, *widths, *filter_width, g)?;
            (
                (0.5 * primitive[left].density + 0.5 * primitive[right].density) * nu,
                0.0,
            )
        } else {
            let coefficient = 0.5 * self.dynamic_coefficient_squared[left]
                + 0.5 * self.dynamic_coefficient_squared[right];
            let heat_coefficient = match &self.closure {
                SubgridClosure::DynamicThermal(dynamic) => {
                    dynamic.face_heat_coefficient(left, right)
                }
                SubgridClosure::ChannelDynamicThermal(dynamic) => {
                    dynamic.face_heat_coefficient(left, right)
                }
                _ => 0.0,
            };
            if coefficient == 0.0 && heat_coefficient == 0.0 {
                (0.0, 0.0)
            } else {
                let width = match &self.closure {
                    SubgridClosure::Local { filter_width, .. } => *filter_width,
                    SubgridClosure::Dynamic(dynamic) => dynamic.filter_width(),
                    SubgridClosure::Lagrangian(dynamic) => dynamic.filter_width(),
                    SubgridClosure::DynamicThermal(dynamic) => dynamic.filter_width(),
                    SubgridClosure::ChannelDynamicThermal(dynamic) => dynamic.filter_width(),
                };
                let g = face_gradient(
                    primitive[left].velocity,
                    primitive[right].velocity,
                    gradient[left],
                    gradient[right],
                    axis,
                    distance,
                );
                // Resolve the physical strain invariant once, then preserve
                // the scale-safe (sqrt(C) * Delta)^2 grouping for each fit.
                let shear_rate = eddy_viscosity(SgsModel::Smagorinsky { cs: 1.0 }, 1.0, g)
                    .map_err(MolecularError::Subgrid)?;
                // Both independent fits act on the conservative face tensor.
                // Neither a centered strain null nor a zero momentum fit may
                // erase turbulent heat transfer resolved by this face.
                let density = 0.5 * primitive[left].density + 0.5 * primitive[right].density;
                let transport = |coefficient: f64| -> Result<f64, MolecularError> {
                    if coefficient == 0.0 {
                        return Ok(0.0);
                    }
                    let nu = smagorinsky_viscosity(
                        fs_math::det::sqrt(coefficient),
                        width,
                        shear_rate,
                    )
                    .map_err(MolecularError::Subgrid)?;
                    Ok(density * nu)
                };
                (transport(coefficient)?, transport(heat_coefficient)?)
            }
        };
        let (prandtl, maximum_ratio) = match &self.closure {
            SubgridClosure::Local { config, .. } => {
                (config.turbulent_prandtl, config.maximum_eddy_to_molecular)
            }
            SubgridClosure::Dynamic(dynamic) => (
                dynamic.config.turbulent_prandtl,
                dynamic.config.model.maximum_eddy_to_molecular,
            ),
            SubgridClosure::Lagrangian(dynamic) => {
                let config = dynamic
                    .fixed_prandtl_config()
                    .expect("independent Lagrangian heat was evaluated above");
                (
                    config.turbulent_prandtl,
                    config.model.maximum_eddy_to_molecular,
                )
            }
            SubgridClosure::DynamicThermal(dynamic) => {
                return dynamic.face_coefficients(
                    gas,
                    axis,
                    face,
                    temperature,
                    pressure,
                    mu_t,
                    rho_alpha,
                    molecular_mu,
                    molecular_conductivity,
                );
            }
            SubgridClosure::ChannelDynamicThermal(dynamic) => {
                return dynamic.face_coefficients(
                    gas,
                    axis,
                    face,
                    temperature,
                    pressure,
                    mu_t,
                    rho_alpha,
                    molecular_mu,
                    molecular_conductivity,
                );
            }
        };
        add_eddy_transport(
            gas,
            axis,
            face,
            temperature,
            pressure,
            mu_t,
            molecular_mu,
            molecular_conductivity,
            prandtl,
            maximum_ratio,
        )
    }
}

#[allow(clippy::too_many_arguments)]
fn add_eddy_transport(
    gas: &FrozenGas,
    axis: usize,
    face: usize,
    temperature: f64,
    pressure: f64,
    mu_t: f64,
    molecular_mu: f64,
    molecular_conductivity: f64,
    prandtl: f64,
    maximum_ratio: f64,
) -> Result<(f64, f64, f64), MolecularError> {
    // Preserve the original molecular arithmetic and source calls exactly.
    if mu_t == 0.0 {
        return Ok((molecular_mu, molecular_conductivity, 0.0));
    }
    let ratio = mu_t / molecular_mu;
    if !ratio.is_finite() || ratio > maximum_ratio {
        return Err(MolecularError::SubgridViscosityRatio {
            axis,
            face,
            observed: ratio,
            maximum: maximum_ratio,
        });
    }
    let cp = gas
        .source()
        .evaluate_properties(Temperature::new(temperature), Pressure::new(pressure))
        .map_err(Error::from)?
        .mass_specific_cp()
        .value();
    let mu = molecular_mu + mu_t;
    let conductivity = molecular_conductivity + mu_t * (cp / prandtl);
    if !positive(cp) || !mu.is_finite() || !conductivity.is_finite() {
        return Err(MolecularError::Invalid(
            "unrepresentable subgrid stress or heat coefficient",
        ));
    }
    Ok((mu, conductivity, mu_t))
}
