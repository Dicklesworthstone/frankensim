//! Explicit subgrid stress, heat and species transport on the original mixture.
//!
//! The Maxwell--Stefan flux and its Soret/Dufour terms remain source-owned.
//! This exploration closure adds `mu_t = rho nu_t`, deviatoric stress/work,
//! `k_t = mu_t cp/Pr_t`, and common-Schmidt mass-frame species fluxes
//! `j_i^t = -(mu_t/Sc_t) grad(Y_i)`. The same original referenced partial
//! enthalpies accompany those species fluxes in conservative total energy.
//! There is no separate viscous or chemical heating source.
//!
//! Physical wall faces retain the original molecular law: modeled wall stress,
//! heat and species flux are zero. Wall-resolved grids are required; no wall
//! function or turbulent reaction-rate closure is inferred. The existing
//! reactive split integrator can compose this operator without changing its
//! chemistry owner. Isotropic SGS energy, backscatter, pressure dilatation,
//! flame wrinkling and shock/turbulence interaction are not represented.
//!
//! The heat/species energy split follows CFX Theory Guide 2025 R1, §1.5.7,
//! equations (1-165)--(1-173):
//! https://ansyshelp.ansys.com/public/Views/Secured/corp/v251/en/cfx_thry/BGBHJAGD.html
//! The shared velocity closures retain their original Nicoud--Ducros and
//! Vreman references in `fs-lbm::sgs3`.

use super::super::sum;
use super::*;
pub use crate::compressible::spatial::molecular::sgs::{SgsModel, SgsUse, SubgridConfig};
use fs_gas::multicomponent::dissipation::{MaxwellStefanSolveV1, TransportFaceReport};
use fs_gas::multicomponent::{Error as FlowError, MAX_SPECIES};
use fs_lbm::sgs3::{eddy_viscosity, vreman_eddy_viscosity};
use fs_qty::Temperature;
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

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
        self.sgs.as_ref().map(|s| s.config)
    }
}

pub(super) struct SubgridTransport {
    config: MixtureSubgridConfig,
    widths: [f64; 3],
    filter_width: f64,
    dynamic_viscosity: Vec<f64>,
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
            config,
            widths,
            filter_width,
            dynamic_viscosity: allocated(grid.len(), 0.0)?,
            heat_capacity_v_density: allocated(grid.len(), 0.0)?,
        })
    }

    pub(super) fn workspace_bytes(&self) -> usize {
        (self.dynamic_viscosity.capacity() + self.heat_capacity_v_density.capacity())
            * size_of::<f64>()
    }

    pub(super) fn prepare(
        &mut self,
        gas: &MixtureGas,
        states: &[f64],
        primitive: &[Primitive],
        gradient: &[[[f64; 3]; 3]],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), DissipationError> {
        if self.config.transport.model.coefficient() == 0.0 {
            return Ok(());
        }
        let n = gas.species_count();
        for (i, ((p, g), state)) in primitive
            .iter()
            .zip(gradient)
            .zip(states.chunks_exact(n + 4))
            .enumerate()
        {
            poll(cancelled)?;
            // The field stores derivative-first; the shared SGS tensor is
            // velocity-component-first. Vreman's directional metric needs this.
            let velocity_gradient =
                std::array::from_fn(|component| std::array::from_fn(|axis| g[axis][component]));
            let nu = match self.config.transport.model {
                SgsModel::Vreman { cv } => {
                    vreman_eddy_viscosity(cv, self.widths, velocity_gradient)
                }
                model => eddy_viscosity(model, self.filter_width, velocity_gradient),
            }
            .map_err(|_| {
                DissipationError::Invalid("mixture subgrid velocity gradient or viscosity")
            })?;
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
            let mu = p.density * nu;
            if !mu.is_finite() || !positive(cv) {
                return Err(DissipationError::Invalid(
                    "unrepresentable mixture subgrid coefficients",
                ));
            }
            self.dynamic_viscosity[i] = mu;
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
        mut receipt: TransportFaceReport,
        controls: MaxwellStefanSolveV1,
        out: &mut [f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<TransportFaceReport, DissipationError> {
        let mu = 0.5 * self.dynamic_viscosity[left] + 0.5 * self.dynamic_viscosity[right];
        // An exactly disabled model makes no extra caloric calls/arithmetic.
        if mu == 0.0 {
            return Ok(receipt);
        }
        poll(cancelled)?;
        let ratio = mu / receipt.dynamic_viscosity_pa_s;
        if !ratio.is_finite() || ratio > self.config.transport.maximum_eddy_to_molecular {
            return Err(DissipationError::Invalid(
                "mixture subgrid eddy/molecular viscosity ratio exceeded",
            ));
        }
        let n = gas.species_count();
        let width = n + 4;
        let pl = primitive[left];
        let pr = primitive[right];
        let temperature = receipt.enthalpy_temperature_k;
        let mut y = [0.0; MAX_SPECIES];
        let mut species_flux = [0.0; MAX_SPECIES];
        let mass_diffusion = mu / self.config.turbulent_schmidt;
        for k in 0..n {
            let yl = states[left * width + k] / pl.density;
            let yr = states[right * width + k] / pr.density;
            y[k] = 0.5 * yl + 0.5 * yr;
            species_flux[k] = -mass_diffusion * ((yr - yl) / distance);
        }
        let pivot = (0..n)
            .max_by(|&a, &b| y[a].total_cmp(&y[b]))
            .expect("admitted species bank");
        // Close a present species, so exactly absent species stay absent. With
        // common Schmidt this changes only summation roundoff, not the model.
        y[pivot] = 1.0 - sum((0..n).filter(|&k| k != pivot).map(|k| y[k]));
        species_flux[pivot] = -sum((0..n).filter(|&k| k != pivot).map(|k| species_flux[k]));
        let mut cp = [0.0; MAX_SPECIES];
        let mut enthalpy = [0.0; MAX_SPECIES];
        for (k, model) in gas.source().models().iter().enumerate() {
            if k % 16 == 0 {
                poll(cancelled)?;
            }
            if y[k] == 0.0 {
                continue;
            }
            let inverse_mass = model.molar_mass().value().recip();
            cp[k] = y[k]
                * inverse_mass
                * model
                    .evaluate_properties(Temperature::new(temperature))
                    .map_err(FlowError::from)?
                    .cp()
                    .value();
            // Evaluate changes within each original caloric source. Subtracting
            // two large formation enthalpies would destroy gauge invariance.
            enthalpy[k] = inverse_mass
                * (model
                    .molar_internal_energy_change(
                        Temperature::new(gas.domain().minimum_temperature),
                        Temperature::new(temperature),
                    )
                    .map_err(FlowError::from)?
                    .value()
                    + R * temperature);
        }
        let cp = sum(cp[..n].iter().copied());
        let conductivity = mu * (cp / self.config.transport.turbulent_prandtl);
        let heat = -conductivity * ((pr.temperature - pl.temperature) / distance);
        let mut energy = heat + sum((0..n).map(|k| enthalpy[k] * species_flux[k]));
        for k in 0..n {
            out[k] += species_flux[k];
        }
        let [a, b] = tangents(axis);
        let divergence = (pr.velocity[axis] - pl.velocity[axis]) / distance
            + 0.5 * gradient[left][a][a]
            + 0.5 * gradient[right][a][a]
            + 0.5 * gradient[left][b][b]
            + 0.5 * gradient[right][b][b];
        for (frame, component) in [axis, a, b].into_iter().enumerate() {
            let normal = (pr.velocity[component] - pl.velocity[component]) / distance;
            let transpose = if component == axis {
                normal - (2.0 / 3.0) * divergence
            } else {
                0.5 * gradient[left][component][axis] + 0.5 * gradient[right][component][axis]
            };
            let stress = mu * (normal + transpose);
            out[n + frame] -= stress;
            energy -= (0.5 * pl.velocity[component] + 0.5 * pr.velocity[component]) * stress;
        }
        out[n + 3] += energy;
        receipt.reduced_heat_flux_w_m2 += heat;
        receipt.mass_frame_residual_kg_m2_s = sum(out[..n].iter().copied());
        let frame_scale = sum(out[..n].iter().map(|j| j.abs()));
        if receipt.mass_frame_residual_kg_m2_s.abs()
            > controls.absolute_mass_flux_tolerance_kg_m2_s
                + controls.relative_tolerance * frame_scale
        {
            return Err(DissipationError::Invalid(
                "mixture subgrid total mass-frame residual",
            ));
        }
        let cv_min = self.heat_capacity_v_density[left].min(self.heat_capacity_v_density[right]);
        let rho_min = pl.density.min(pr.density);
        let gamma = 1.0
            + (pl.pressure / (pl.temperature * self.heat_capacity_v_density[left]))
                .max(pr.pressure / (pr.temperature * self.heat_capacity_v_density[right]));
        let enthalpy_span = (0..n)
            .filter(|&k| y[k] > 0.0)
            .map(|k| (enthalpy[k] - enthalpy[pivot]).abs())
            .fold(0.0, f64::max);
        let species_scale = (mass_diffusion / rho_min)
            * gamma
            * (1.0 + (enthalpy_span / temperature) * (pl.density.max(pr.density) / cv_min));
        // Sum molecular and modeled bounds: max(D_m,D_t) alone underestimates
        // their combined diffusion. Include enthalpy/species coupling and the
        // full 3-D stress row sum, then let the owner apply all inverse lengths.
        receipt.diffusivity_m2_s +=
            ((20.0 / 3.0) * mu / rho_min).max(conductivity / cv_min + species_scale);
        if !positive(cp)
            || !conductivity.is_finite()
            || !mass_diffusion.is_finite()
            || !out.iter().all(|v| v.is_finite())
            || !receipt.reduced_heat_flux_w_m2.is_finite()
            || !receipt.diffusivity_m2_s.is_finite()
            || !receipt.mass_frame_residual_kg_m2_s.is_finite()
            || !frame_scale.is_finite()
            || !gamma.is_finite()
            || !enthalpy_span.is_finite()
            || !species_scale.is_finite()
        {
            return Err(DissipationError::Invalid(
                "unrepresentable mixture subgrid face transport",
            ));
        }
        poll(cancelled)?;
        Ok(receipt)
    }
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}

#[cfg(test)]
mod tests;
