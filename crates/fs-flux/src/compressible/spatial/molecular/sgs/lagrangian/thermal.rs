//! Independent pathline memory for conservative temperature transport.
//!
//! The scalar coefficient is fitted from the scalar Germano identity and held
//! with the momentum coefficient over the physical interval. Each face uses
//! its current density and full velocity-gradient strain:
//! alpha_t = C_phi,face Delta^2 |S_face|, k_t = rho_face Cp(T_face) alpha_t.
//! No division by momentum viscosity or inferred turbulent Prandtl is needed.
//! The heat fit therefore remains meaningful when the momentum fit is zero.

use super::*;

/// Independent momentum and temperature memories, in physical seconds.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LagrangianThermalSubgridConfig {
    pub momentum: LagrangianSgsConfig,
    pub heat: LagrangianScalarConfig,
}

impl LagrangianThermalSubgridConfig {
    pub fn validate(self) -> Result<(), MolecularError> {
        self.momentum.validate()?;
        self.heat.validate()?;
        Ok(())
    }
}

impl SubgridTransport {
    /// 264 bytes/cell for the paired history/filter, plus 56 for current
    /// density/velocity/temperature and both held face coefficients.
    pub(in crate::compressible::spatial::molecular) fn required_lagrangian_thermal_bytes(
        flow: &Euler3d,
    ) -> Result<usize, MolecularError> {
        lagrangian_bytes(flow.states.len(), true)
    }

    pub(in crate::compressible::spatial::molecular) fn new_lagrangian_thermal(
        flow: &Euler3d,
        config: LagrangianThermalSubgridConfig,
        maximum_bytes: usize,
    ) -> Result<Self, MolecularError> {
        config.validate()?;
        Self::new_lagrangian_owner(flow, LagrangianConfig::Thermal(config), maximum_bytes)
    }

    pub(in crate::compressible::spatial::molecular) fn lagrangian_thermal_config(
        &self,
    ) -> Option<LagrangianThermalSubgridConfig> {
        match &self.closure {
            SubgridClosure::Lagrangian(dynamic) => dynamic.thermal_config(),
            _ => None,
        }
    }

    pub(in crate::compressible::spatial::molecular) fn lagrangian_heat_report(
        &self,
    ) -> Option<DynamicScalarReport> {
        match &self.closure {
            SubgridClosure::Lagrangian(dynamic) => dynamic.last_heat_report,
            _ => None,
        }
    }
}

impl LagrangianTransport {
    pub(in crate::compressible::spatial::molecular::sgs) fn thermal_config(
        &self,
    ) -> Option<LagrangianThermalSubgridConfig> {
        match self.config {
            LagrangianConfig::Thermal(config) => Some(config),
            LagrangianConfig::FixedPrandtl(_) => None,
        }
    }

    #[allow(clippy::too_many_arguments)]
    pub(in crate::compressible::spatial::molecular::sgs) fn thermal_face_coefficients(
        &self,
        config: LagrangianThermalSubgridConfig,
        gas: &FrozenGas,
        left: usize,
        right: usize,
        axis: usize,
        face: usize,
        distance: f64,
        momentum_coefficients: &[f64],
        primitive: &[Primitive],
        gradient: &[[[f64; 3]; 3]],
        temperature: f64,
        pressure: f64,
        molecular_mu: f64,
        molecular_conductivity: f64,
    ) -> Result<(f64, f64), MolecularError> {
        let momentum = 0.5 * momentum_coefficients[left] + 0.5 * momentum_coefficients[right];
        let heat = 0.5 * self.heat_coefficient[left] + 0.5 * self.heat_coefficient[right];
        if momentum == 0.0 && heat == 0.0 {
            return Ok((molecular_mu, molecular_conductivity));
        }
        // Reconstruct once, at the exact conservative face. A zero centered
        // cell strain must not erase a shear resolved by this face stencil.
        let g = face_gradient(
            primitive[left].velocity,
            primitive[right].velocity,
            gradient[left],
            gradient[right],
            axis,
            distance,
        );
        let density = 0.5 * primitive[left].density + 0.5 * primitive[right].density;
        let transport = |coefficient: f64| -> Result<f64, MolecularError> {
            if coefficient == 0.0 {
                return Ok(0.0);
            }
            let value = density
                * eddy_viscosity(
                    SgsModel::Smagorinsky {
                        cs: fs_math::det::sqrt(coefficient),
                    },
                    self.filter_width,
                    g,
                )
                .map_err(MolecularError::Subgrid)?;
            if !value.is_finite() {
                return Err(MolecularError::Invalid(
                    "unrepresentable Lagrangian face transport",
                ));
            }
            Ok(value)
        };
        let mu_t = transport(momentum)?;
        let rho_alpha = transport(heat)?;
        if mu_t != 0.0 {
            let ratio = mu_t / molecular_mu;
            if !ratio.is_finite() || ratio > config.momentum.maximum_eddy_to_molecular {
                return Err(MolecularError::SubgridViscosityRatio {
                    axis,
                    face,
                    observed: ratio,
                    maximum: config.momentum.maximum_eddy_to_molecular,
                });
            }
        }
        let modeled_conductivity = if rho_alpha == 0.0 {
            0.0
        } else {
            let cp = gas
                .source()
                .evaluate_properties(Temperature::new(temperature), Pressure::new(pressure))
                .map_err(Error::from)?
                .mass_specific_cp()
                .value();
            if !positive(cp) {
                return Err(MolecularError::Invalid("invalid Lagrangian heat source Cp"));
            }
            let value = rho_alpha * cp;
            let ratio = value / molecular_conductivity;
            if !ratio.is_finite() || ratio > config.heat.maximum_diffusivity_to_molecular {
                return Err(MolecularError::SubgridThermalDiffusivityRatio {
                    axis,
                    face,
                    observed: ratio,
                    maximum: config.heat.maximum_diffusivity_to_molecular,
                });
            }
            value
        };
        let mu = molecular_mu + mu_t;
        let conductivity = molecular_conductivity + modeled_conductivity;
        if !mu.is_finite() || !conductivity.is_finite() {
            return Err(MolecularError::Invalid(
                "unrepresentable Lagrangian stress or heat coefficient",
            ));
        }
        Ok((mu, conductivity))
    }
}
