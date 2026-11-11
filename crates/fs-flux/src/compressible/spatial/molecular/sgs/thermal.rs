//! Independent dynamic momentum and temperature diffusion on one gas owner.
//!
//! The shared Favre Germano scalar fit retains C_phi independently of the
//! centered cell strain. Each conservative face uses its full velocity tensor:
//! `alpha_t = C_phi,face Delta^2 |S_face|` and
//! `k_t = rho_face alpha_t Cp(T_face)` from the original caloric source.
//! A heat fit can remain active when the momentum fit is zero. No division
//! by the momentum coefficient, inferred Prandtl floor, or second heating
//! source is introduced. The existing total-energy flux and stability estimate
//! consume the resulting conductivity. This remains a low-Mach temperature
//! gradient-diffusion closure; compressible pressure-dilatation and SGS energy
//! transport are not supplied. Fully periodic, uniform Cartesian grids with
//! at least three cells per axis are required; each physical cell width enters
//! both fits, the conservative flux and the directional stability estimate.

use super::*;
use fs_lbm::sgs3::dynamic::PeriodicDynamicSgs;
use fs_lbm::sgs3::dynamic::scalar::PeriodicDynamicScalar;

/// Separately fitted stress and heat, with independent averaging and bounds.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DynamicThermalSubgridConfig {
    pub momentum: DynamicSgsConfig,
    pub heat: DynamicScalarConfig,
}

impl DynamicThermalSubgridConfig {
    pub fn validate(self) -> Result<(), MolecularError> {
        self.momentum.validate()?;
        self.heat.validate()?;
        Ok(())
    }
}

impl<L: TransportLaw> MolecularTransport3d<L> {
    /// Additional numerical payload for independently fitted momentum/heat:
    /// 256 bytes/cell, excluding existing molecular and Euler workspaces.
    pub fn required_dynamic_thermal_sgs_workspace_bytes(
        flow: &Euler3d,
    ) -> Result<usize, MolecularError> {
        required_bytes(flow.states.len())
    }

    /// Fit momentum and temperature from the same physical SSP-stage fields.
    /// Reuses the original conservative stress/work/heat and parabolic bound.
    /// The byte limit covers every added SGS array's allocation capacity.
    pub fn with_dynamic_thermal_sgs(
        flow: &Euler3d,
        law: L,
        config: DynamicThermalSubgridConfig,
        maximum_subgrid_workspace_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, MolecularError> {
        poll(&mut cancelled)?;
        let sgs =
            SubgridTransport::new_dynamic_thermal(flow, config, maximum_subgrid_workspace_bytes)?;
        poll(&mut cancelled)?;
        let mut result = Self::new(flow, law, &mut cancelled)?;
        result.sgs = Some(sgs);
        poll(&mut cancelled)?;
        Ok(result)
    }

    /// Exact independently fitted momentum/heat controls. Fixed-Prandtl and
    /// molecular closures return None; their existing getters stay distinct.
    pub fn dynamic_thermal_subgrid_config(&self) -> Option<DynamicThermalSubgridConfig> {
        match &self.sgs.as_ref()?.closure {
            SubgridClosure::DynamicThermal(dynamic) => Some(dynamic.config),
            _ => None,
        }
    }
}

pub(super) fn required_bytes(cells: usize) -> Result<usize, MolecularError> {
    let fields =
        cells
            .checked_mul(7 * core::mem::size_of::<f64>())
            .ok_or(MolecularError::Invalid(
                "dynamic thermal SGS workspace overflow",
            ))?;
    PeriodicDynamicSgs::required_bytes(cells)?
        .checked_add(PeriodicDynamicScalar::required_bytes(cells)?)
        .and_then(|bytes| bytes.checked_add(fields))
        .ok_or(MolecularError::Invalid(
            "dynamic thermal SGS workspace overflow",
        ))
}

pub(super) struct DynamicThermalTransport {
    pub(super) config: DynamicThermalSubgridConfig,
    momentum: PeriodicDynamicSgs,
    heat: PeriodicDynamicScalar,
    work: ThermalWorkspace,
}

/// Same-stage primitive staging and conservative coefficients shared by the
/// periodic and channel fits. Geometry remains with each paired filter owner.
pub(super) struct ThermalWorkspace {
    pub(super) density: Vec<f64>,
    pub(super) velocity: Vec<[f64; 3]>,
    pub(super) temperature: Vec<f64>,
    heat_coefficient: Vec<f64>,
}

impl SubgridTransport {
    pub(in crate::compressible::spatial::molecular) fn new_dynamic_thermal(
        flow: &Euler3d,
        config: DynamicThermalSubgridConfig,
        maximum_bytes: usize,
    ) -> Result<Self, MolecularError> {
        config.validate()?;
        dynamic::admit_domain(flow)?;
        let n = flow.states.len();
        dynamic::admit_bytes(required_bytes(n)?, maximum_bytes)?;
        let closure = DynamicThermalTransport {
            config,
            momentum: PeriodicDynamicSgs::new_cartesian(
                flow.grid.shape,
                flow.grid.spacing,
                maximum_bytes,
            )?,
            heat: PeriodicDynamicScalar::new(flow.grid.shape, maximum_bytes)?,
            work: ThermalWorkspace::new(n)?,
        };
        let dynamic_coefficient_squared = allocated(n, 0.0)?;
        let bytes = closure.workspace_bytes()
            + dynamic_coefficient_squared.capacity() * core::mem::size_of::<f64>();
        dynamic::admit_bytes(bytes, maximum_bytes)?;
        Ok(Self {
            closure: SubgridClosure::DynamicThermal(closure),
            wall_models: [None; 6],
            dynamic_coefficient_squared,
        })
    }
}

impl DynamicThermalTransport {
    /// Owned capacities of the momentum and heat filters and the staged fields.
    pub(super) fn workspace_bytes(&self) -> usize {
        self.momentum.workspace_bytes() + self.heat.workspace_bytes() + self.work.workspace_bytes()
    }

    pub(super) fn filter_width(&self) -> f64 {
        self.momentum.filter_width()
    }

    pub(super) fn prepare(
        &mut self,
        primitive: &[Primitive],
        dynamic_coefficient_squared: &mut [f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), MolecularError> {
        self.work.stage(primitive, cancelled)?;
        self.momentum.prepare_with_scalar(
            self.config.momentum,
            self.config.heat,
            &self.work.density,
            &self.work.velocity,
            &self.work.temperature,
            &mut self.heat,
            &mut *cancelled,
        )?;
        self.work.store_transport(
            self.momentum.coefficients_squared(),
            self.heat.coefficients(),
            dynamic_coefficient_squared,
            cancelled,
        )
    }

    pub(super) fn face_heat_coefficient(&self, left: usize, right: usize) -> f64 {
        self.work.face_heat_coefficient(left, right)
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn face_coefficients(
        &self,
        gas: &FrozenGas,
        axis: usize,
        face: usize,
        temperature: f64,
        pressure: f64,
        mu_t: f64,
        rho_alpha: f64,
        molecular_mu: f64,
        molecular_conductivity: f64,
    ) -> Result<(f64, f64, f64), MolecularError> {
        ThermalWorkspace::face_coefficients(
            gas,
            axis,
            face,
            temperature,
            pressure,
            mu_t,
            rho_alpha,
            molecular_mu,
            molecular_conductivity,
            self.config.momentum.maximum_eddy_to_molecular,
            self.config.heat.maximum_diffusivity_to_molecular,
        )
    }
}

impl ThermalWorkspace {
    pub(super) fn new(cells: usize) -> Result<Self, MolecularError> {
        Ok(Self {
            density: allocated(cells, 0.0)?,
            velocity: allocated(cells, [0.0; 3])?,
            temperature: allocated(cells, 0.0)?,
            heat_coefficient: allocated(cells, 0.0)?,
        })
    }

    pub(super) fn workspace_bytes(&self) -> usize {
        (self.density.capacity() + self.temperature.capacity() + self.heat_coefficient.capacity())
            * core::mem::size_of::<f64>()
            + self.velocity.capacity() * core::mem::size_of::<[f64; 3]>()
    }

    pub(super) fn stage(
        &mut self,
        primitive: &[Primitive],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), MolecularError> {
        for (cell, p) in primitive.iter().enumerate() {
            if cell % 64 == 0 {
                poll(cancelled)?;
            }
            self.density[cell] = p.density;
            self.velocity[cell] = p.velocity;
            self.temperature[cell] = p.temperature;
        }
        Ok(())
    }

    pub(super) fn store_transport(
        &mut self,
        coefficients: impl ExactSizeIterator<Item = f64>,
        heat_coefficients: impl ExactSizeIterator<Item = f64>,
        dynamic_coefficient_squared: &mut [f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), MolecularError> {
        for (cell, ((heat, fitted), coefficient)) in heat_coefficients
            .zip(coefficients)
            .zip(dynamic_coefficient_squared)
            .enumerate()
        {
            if cell % 64 == 0 {
                poll(cancelled)?;
            }
            *coefficient = fitted;
            self.heat_coefficient[cell] = heat;
        }
        poll(cancelled)?;
        Ok(())
    }

    pub(super) fn face_heat_coefficient(&self, left: usize, right: usize) -> f64 {
        0.5 * self.heat_coefficient[left] + 0.5 * self.heat_coefficient[right]
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn face_coefficients(
        gas: &FrozenGas,
        axis: usize,
        face: usize,
        temperature: f64,
        pressure: f64,
        mu_t: f64,
        rho_alpha: f64,
        molecular_mu: f64,
        molecular_conductivity: f64,
        maximum_viscosity_ratio: f64,
        maximum_diffusivity_ratio: f64,
    ) -> Result<(f64, f64, f64), MolecularError> {
        if mu_t == 0.0 && rho_alpha == 0.0 {
            return Ok((molecular_mu, molecular_conductivity, 0.0));
        }
        if mu_t != 0.0 {
            let ratio = mu_t / molecular_mu;
            if !ratio.is_finite() || ratio > maximum_viscosity_ratio {
                return Err(MolecularError::SubgridViscosityRatio {
                    axis,
                    face,
                    observed: ratio,
                    maximum: maximum_viscosity_ratio,
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
                return Err(MolecularError::Invalid("invalid dynamic heat source Cp"));
            }
            let conductivity = rho_alpha * cp;
            let ratio = conductivity / molecular_conductivity;
            if !ratio.is_finite() || ratio > maximum_diffusivity_ratio {
                return Err(MolecularError::SubgridThermalDiffusivityRatio {
                    axis,
                    face,
                    observed: ratio,
                    maximum: maximum_diffusivity_ratio,
                });
            }
            conductivity
        };
        let mu = molecular_mu + mu_t;
        let conductivity = molecular_conductivity + modeled_conductivity;
        if !mu.is_finite() || !conductivity.is_finite() {
            return Err(MolecularError::Invalid(
                "unrepresentable dynamic stress or heat coefficient",
            ));
        }
        Ok((mu, conductivity, mu_t))
    }
}
