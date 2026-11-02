//! Local equilibrium subgrid stress and heat transport on the conservative gas.
//!
//! The resolved velocity is momentum/density (the Favre velocity). Each SSP
//! stage recomputes the shared LBM Smagorinsky/WALE/Vreman closure from its full
//! physical velocity gradient. Face stress adds `mu_t = rho nu_t` only to the
//! deviatoric viscosity; bulk viscosity remains the original molecular law.
//! Gradient-diffusion heat adds `k_t = mu_t cp(T)/Pr_t`, using the original
//! caloric source at the face. Stress work and heat use the existing shared
//! total-energy flux, so resolved kinetic loss is not counted as heat twice.
//!
//! This is an explicit exploration model for low-Mach, frozen-composition
//! flow, not a closure for shock/turbulence interaction. Isotropic SGS kinetic
//! energy, pressure dilatation, backscatter, species mixing and wall models
//! are absent. Reflecting boundaries retain their adiabatic slip semantics.
//! No coefficient or turbulent Prandtl number is chosen implicitly.
//!
//! References: Nicoud & Ducros (1999), doi:10.1023/A:1009995426001;
//! Vreman (2004), doi:10.1063/1.1785131. The stress and heat closure follow
//! Nicoud, "Unsteady flows modeling and computation" (2007), equations
//! (24)-(26): https://imag.umontpellier.fr/~nicoud/PDF/notes_Nicoud_1and2.pdf.

use super::*;
use fs_lbm::sgs3::{SgsConfig, eddy_viscosity, vreman_eddy_viscosity};
pub use fs_lbm::sgs3::{SgsModel, SgsUse};

mod dynamic;
pub use dynamic::DynamicSubgridConfig;
pub use fs_lbm::sgs3::dynamic::{DynamicAveraging, DynamicSgsConfig, DynamicSgsReport};

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

    /// Exact local controls; `None` for molecular-only or dynamic SGS flow.
    pub fn subgrid_config(&self) -> Option<SubgridConfig> {
        self.sgs.as_ref().and_then(|s| match &s.closure {
            SubgridClosure::Local { config, .. } => Some(*config),
            SubgridClosure::Dynamic(_) => None,
        })
    }
}

pub(super) struct SubgridTransport {
    closure: SubgridClosure,
    // Reused at every stage. No allocation in the per-cell/face hot path.
    dynamic_viscosity: Vec<f64>,
}

enum SubgridClosure {
    Local {
        config: SubgridConfig,
        widths: [f64; 3],
        filter_width: f64,
    },
    Dynamic(dynamic::DynamicTransport),
}

impl SubgridClosure {
    fn heat_and_ratio(&self) -> (f64, f64) {
        match self {
            Self::Local { config, .. } => {
                (config.turbulent_prandtl, config.maximum_eddy_to_molecular)
            }
            Self::Dynamic(dynamic) => (
                dynamic.config.turbulent_prandtl,
                dynamic.config.model.maximum_eddy_to_molecular,
            ),
        }
    }
}

impl SubgridTransport {
    pub(super) fn new(flow: &Euler3d, config: SubgridConfig) -> Result<Self, MolecularError> {
        let widths = flow.grid.spacing.map(|h| h * config.filter_width_factor);
        let filter_width = flow.grid.cell_volume().cbrt() * config.filter_width_factor;
        if !widths.into_iter().chain([filter_width]).all(positive) {
            return Err(MolecularError::Invalid(
                "unrepresentable subgrid filter widths",
            ));
        }
        Ok(Self {
            closure: SubgridClosure::Local {
                config,
                widths,
                filter_width,
            },
            dynamic_viscosity: allocated(flow.states.len(), 0.0)?,
        })
    }

    pub(super) fn prepare(
        &mut self,
        primitive: &[Primitive],
        gradient: &[[[f64; 3]; 3]],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), MolecularError> {
        let (config, widths, filter_width) = match &mut self.closure {
            SubgridClosure::Local {
                config,
                widths,
                filter_width,
            } => (*config, *widths, *filter_width),
            SubgridClosure::Dynamic(dynamic) => {
                return dynamic.prepare(primitive, &mut self.dynamic_viscosity, cancelled);
            }
        };
        if config.model.coefficient() == 0.0 {
            return Ok(());
        }
        for (cell, ((p, g), mu)) in primitive
            .iter()
            .zip(gradient)
            .zip(&mut self.dynamic_viscosity)
            .enumerate()
        {
            if cell % 64 == 0 {
                poll(cancelled)?;
            }
            let nu = match config.model {
                SgsModel::Vreman { cv } => vreman_eddy_viscosity(cv, widths, *g),
                model => eddy_viscosity(model, filter_width, *g),
            }
            .map_err(MolecularError::Subgrid)?;
            *mu = p.density * nu;
            if !mu.is_finite() {
                return Err(MolecularError::Invalid(
                    "unrepresentable subgrid dynamic viscosity",
                ));
            }
        }
        Ok(())
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
        molecular_mu: f64,
        molecular_conductivity: f64,
    ) -> Result<(f64, f64), MolecularError> {
        let mu_t = 0.5 * self.dynamic_viscosity[left] + 0.5 * self.dynamic_viscosity[right];
        // Preserve the original molecular arithmetic and source calls exactly.
        if mu_t == 0.0 {
            return Ok((molecular_mu, molecular_conductivity));
        }
        let ratio = mu_t / molecular_mu;
        let (prandtl, maximum_ratio) = self.closure.heat_and_ratio();
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
        Ok((mu, conductivity))
    }
}
