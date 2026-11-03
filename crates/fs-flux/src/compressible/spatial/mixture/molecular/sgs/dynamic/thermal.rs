//! Independently fitted temperature mixing on the original mixture owner.
//!
//! The scalar Germano fit supplies alpha_t directly, without dividing by the
//! momentum coefficient or inventing a Prandtl floor. Temperature is the
//! fitted scalar; this remains a low-Mach gradient-diffusion heat model for
//! evolving mixtures, not a closure for all compressible SGS energy terms.
//! Species mixing retains the explicitly authored common Schmidt number.
//! Momentum and temperature gradients share the actual directional Cartesian
//! cell widths; scalar and momentum test filters use the same physical box.

use super::*;
use fs_lbm::sgs3::dynamic::scalar::PeriodicDynamicScalar;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MixtureDynamicThermalSubgridConfig {
    pub transport: DynamicThermalSubgridConfig,
    /// Positive finite common species Schmidt number; no species fit is made.
    pub turbulent_schmidt: f64,
}

impl MixtureDynamicThermalSubgridConfig {
    pub fn validate(self) -> Result<(), DissipationError> {
        self.transport.validate().map_err(|_| {
            DissipationError::Invalid("invalid dynamic mixture momentum/heat controls")
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
    /// Combined flow and operator numerical payload. Thermal fitting adds a
    /// scalar workspace, temperature staging and rho*alpha_t cache to dynamic
    /// momentum: nominally 264 additional bytes/cell over molecular-only flow.
    pub fn required_workspace_bytes_with_dynamic_thermal_sgs(
        flow: &MixtureEuler3d,
        source: &SpatialTransport<L>,
    ) -> Result<usize, DissipationError> {
        let scalar =
            PeriodicDynamicScalar::required_bytes(flow.cell_count()).map_err(dynamic_error)?;
        let fields = flow.cell_count().checked_mul(2 * size_of::<f64>()).ok_or(
            DissipationError::Invalid("dynamic mixture heat workspace overflow"),
        )?;
        Self::required_workspace_bytes_with_dynamic_sgs(flow, source)?
            .checked_add(scalar)
            .and_then(|bytes| bytes.checked_add(fields))
            .ok_or(DissipationError::Invalid(
                "dynamic mixture heat workspace overflow",
            ))
    }

    /// Fit momentum and heat independently at each original SSP stage. The
    /// original Maxwell--Stefan/species-enthalpy and stress-work owner remains
    /// authoritative. The complete flow/operator capacity obeys the byte bound.
    pub fn with_dynamic_thermal_sgs(
        flow: &MixtureEuler3d,
        source: SpatialTransport<L>,
        config: MixtureDynamicThermalSubgridConfig,
        maximum_workspace_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, DissipationError> {
        poll(&mut cancelled)?;
        config.validate()?;
        admit_domain(flow)?;
        admit_workspace(
            Self::required_workspace_bytes_with_dynamic_thermal_sgs(flow, &source)?,
            maximum_workspace_bytes,
        )?;
        let mut result = Self::new(flow, source, maximum_workspace_bytes, &mut cancelled)?;
        result.sgs = Some(SubgridTransport {
            closure: ViscosityClosure::DynamicThermal(DynamicThermalTransport::new(
                &flow.grid, config,
            )?),
            dynamic_viscosity: allocated(flow.cell_count(), 0.0)?,
            heat_capacity_v_density: allocated(flow.cell_count(), 0.0)?,
        });
        let bytes = result
            .workspace_bytes()
            .checked_add(flow.workspace_bytes())
            .ok_or(DissipationError::Invalid(
                "dynamic mixture heat workspace overflow",
            ))?;
        admit_workspace(bytes, maximum_workspace_bytes)?;
        poll(&mut cancelled)?;
        Ok(result)
    }

    /// Active independently fitted heat/momentum controls, without a fictitious
    /// fixed Prandtl number. Other molecular/SGS producers return None.
    pub fn dynamic_thermal_subgrid_config(&self) -> Option<MixtureDynamicThermalSubgridConfig> {
        match &self.sgs.as_ref()?.closure {
            ViscosityClosure::DynamicThermal(dynamic) => Some(dynamic.config),
            _ => None,
        }
    }
}

pub(in crate::compressible::spatial::mixture::molecular::sgs) struct DynamicThermalTransport {
    pub(in crate::compressible::spatial::mixture::molecular::sgs) config:
        MixtureDynamicThermalSubgridConfig,
    work: ThermalWorkspace,
}

/// Same-stage momentum/heat arrays shared by fixed-Schmidt and fitted-species
/// consumers. Neither consumer owns a second thermodynamic or flux operator.
pub(super) struct ThermalWorkspace {
    pub(super) momentum: PeriodicDynamicSgs,
    pub(super) heat: PeriodicDynamicScalar,
    pub(super) fields: ThermalFields,
}

/// Geometry-independent same-stage fields shared by periodic and channel fits.
pub(super) struct ThermalFields {
    pub(super) density: Vec<f64>,
    pub(super) velocity: Vec<[f64; 3]>,
    pub(super) temperature: Vec<f64>,
    rho_diffusivity: Vec<f64>,
}

impl DynamicThermalTransport {
    fn new(
        grid: &CartesianGrid,
        config: MixtureDynamicThermalSubgridConfig,
    ) -> Result<Self, DissipationError> {
        Ok(Self {
            config,
            work: ThermalWorkspace::new(grid)?,
        })
    }

    pub(in crate::compressible::spatial::mixture::molecular::sgs) fn workspace_bytes(
        &self,
    ) -> usize {
        self.work.workspace_bytes()
    }

    /// True if either independently fitted coefficient is active. A pure heat
    /// fit must still refresh the parent's source Cv and energy coefficients.
    pub(in crate::compressible::spatial::mixture::molecular::sgs) fn prepare(
        &mut self,
        primitive: &[Primitive],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<bool, DissipationError> {
        self.work.fields.load(primitive, cancelled)?;
        let (momentum, heat) = self
            .work
            .momentum
            .prepare_with_scalar(
                self.config.transport.momentum,
                self.config.transport.heat,
                &self.work.fields.density,
                &self.work.fields.velocity,
                &self.work.fields.temperature,
                &mut self.work.heat,
                &mut *cancelled,
            )
            .map_err(dynamic_error)?;
        self.work
            .fields
            .collect_heat(self.work.heat.diffusivities(), cancelled)?;
        Ok(momentum.maximum_eddy_viscosity != 0.0 || heat.maximum_diffusivity != 0.0)
    }

    pub(in crate::compressible::spatial::mixture::molecular::sgs) fn viscosities(&self) -> &[f64] {
        self.work.momentum.viscosities()
    }

    /// Density-weighted thermal diffusivity and its molecular ratio bound.
    pub(in crate::compressible::spatial::mixture::molecular::sgs) fn face_heat(
        &self,
        left: usize,
        right: usize,
    ) -> (f64, f64) {
        (
            self.work.fields.face_heat(left, right),
            self.config.transport.heat.maximum_diffusivity_to_molecular,
        )
    }
}

impl ThermalWorkspace {
    pub(super) fn new(grid: &CartesianGrid) -> Result<Self, DissipationError> {
        let momentum_bytes =
            PeriodicDynamicSgs::required_bytes(grid.len()).map_err(dynamic_error)?;
        let heat_bytes =
            PeriodicDynamicScalar::required_bytes(grid.len()).map_err(dynamic_error)?;
        Ok(Self {
            momentum: PeriodicDynamicSgs::new_cartesian(grid.shape, grid.spacing, momentum_bytes)
                .map_err(dynamic_error)?,
            heat: PeriodicDynamicScalar::new(grid.shape, heat_bytes).map_err(dynamic_error)?,
            fields: ThermalFields::new(grid.len())?,
        })
    }

    pub(super) fn workspace_bytes(&self) -> usize {
        self.momentum.workspace_bytes()
            + self.heat.workspace_bytes()
            + self.fields.workspace_bytes()
    }
}

impl ThermalFields {
    pub(super) fn new(cells: usize) -> Result<Self, DissipationError> {
        Ok(Self {
            density: allocated(cells, 0.0)?,
            velocity: allocated(cells, [0.0; 3])?,
            temperature: allocated(cells, 0.0)?,
            rho_diffusivity: allocated(cells, 0.0)?,
        })
    }

    pub(super) fn workspace_bytes(&self) -> usize {
        (self.density.capacity() + self.temperature.capacity() + self.rho_diffusivity.capacity())
            * size_of::<f64>()
            + self.velocity.capacity() * size_of::<[f64; 3]>()
    }

    pub(super) fn load(
        &mut self,
        primitive: &[Primitive],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), DissipationError> {
        for (i, p) in primitive.iter().enumerate() {
            if i % 64 == 0 {
                poll(cancelled)?;
            }
            self.density[i] = p.density;
            self.velocity[i] = p.velocity;
            self.temperature[i] = p.temperature;
        }
        poll(cancelled)
    }

    pub(super) fn collect_heat(
        &mut self,
        diffusivities: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), DissipationError> {
        for (i, (weight, &alpha)) in self
            .rho_diffusivity
            .iter_mut()
            .zip(diffusivities)
            .enumerate()
        {
            if i % 64 == 0 {
                poll(cancelled)?;
            }
            *weight = self.density[i] * alpha;
            if !weight.is_finite() {
                return Err(DissipationError::Invalid(
                    "unrepresentable dynamic mixture heat coefficient",
                ));
            }
        }
        poll(cancelled)
    }

    pub(super) fn face_heat(&self, left: usize, right: usize) -> f64 {
        0.5 * self.rho_diffusivity[left] + 0.5 * self.rho_diffusivity[right]
    }
}
