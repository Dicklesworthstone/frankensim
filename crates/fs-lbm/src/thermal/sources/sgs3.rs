//! Explicit SGS momentum and temperature diffusion on the existing D3Q19/D3Q7
//! slab. The shared closure supplies `nu_t`; the declared turbulent Prandtl
//! number supplies `alpha_t=nu_t/Pr_t`. Neither is a new material law.
//!
//! Coefficients use the accepted, Guo-force-centered velocity and are frozen
//! over one source/flow/heat step. Stationary halfway plates use their actual
//! zero velocity in the gradient reconstruction. This remains Boussinesq
//! sensible-temperature transport: modeled kinetic dissipation is diagnostic,
//! not an additional heat source or a total-energy conservation claim.
//! The turbulent Prandtl number is constant; no wall function or thermal wall
//! law substitutes for resolving the boundary layers.

use super::{
    Slab, ThermalLbm3, ThermalSourceBalance, ThermalSourceError, ThermalSourceLimits,
    source_step_transport_mapped, transport::RelaxationField,
};
use crate::sgs3::{SgsConfig, SgsError, SgsReport, eddy_viscosity};
use fs_exec::CancelGate;

/// Plane-averaged dynamic momentum and heat on the original halfway slab.
pub mod channel;

/// Explicit exploration controls; no turbulent Prandtl number is defaulted.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermalSgsConfig {
    /// Shared model, lattice filter width, viscosity-ratio bound and usage.
    pub closure: SgsConfig,
    /// Positive finite `Pr_t=nu_t/alpha_t`.
    pub turbulent_prandtl: f64,
    /// Flow-lattice Mach bound in `(0, 0.3]`, checked before and after the step.
    pub maximum_mach: f64,
}

impl ThermalSgsConfig {
    fn validate(self) -> Result<(), ThermalSgsError> {
        self.closure.validate()?;
        if !self.turbulent_prandtl.is_finite() || self.turbulent_prandtl <= 0.0 {
            return Err(ThermalSgsError::Invalid(
                "turbulent Prandtl number must be positive and finite",
            ));
        }
        if !self.maximum_mach.is_finite() || self.maximum_mach <= 0.0 || self.maximum_mach > 0.3 {
            return Err(ThermalSgsError::Invalid(
                "SGS Mach bound must lie in (0, 0.3]",
            ));
        }
        Ok(())
    }
}

/// Reusable SGS preparation storage. The original thermal owner still holds
/// a tentative complete state during its checked source/transport step.
pub struct ThermalSgsWorkspace3 {
    velocity: Vec<[f64; 3]>,
    density: Vec<f64>,
    momentum: Vec<f64>,
    thermal: Vec<f64>,
}

impl ThermalSgsWorkspace3 {
    /// Total relaxation fields, valid only after successful preparation. Do
    /// not pass these back as molecular input on the following step.
    #[must_use]
    pub fn relaxation(&self) -> RelaxationField<'_> {
        RelaxationField {
            momentum: &self.momentum,
            thermal: &self.thermal,
        }
    }
}

/// Pre-step modeled transport and mechanical dissipation diagnostics.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ThermalSgsReport {
    /// Shared shear-only dissipation report, not integrated heat.
    pub momentum: SgsReport,
    /// Largest modeled temperature diffusivity in lattice units.
    pub maximum_eddy_diffusivity_lu: f64,
}

/// Accepted original source/wall heat balance and pre-step SGS diagnostics.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermalSgsStep {
    pub heat: ThermalSourceBalance,
    pub sgs: ThermalSgsReport,
}

/// Refusals leave the solver unchanged; scratch contents may be overwritten.
#[derive(Clone, Debug, PartialEq)]
pub enum ThermalSgsError {
    Invalid(&'static str),
    Allocation,
    Closure(SgsError),
    Dynamic(crate::sgs3::dynamic::DynamicSgsError),
    DiffusivityRatio {
        cell: usize,
        ratio: f64,
        maximum: f64,
    },
    Source(ThermalSourceError),
    Mach {
        cell: usize,
        mach_squared: f64,
        maximum: f64,
    },
    Cancelled,
}
impl core::fmt::Display for ThermalSgsError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "thermal SGS: {self:?}")
    }
}
impl std::error::Error for ThermalSgsError {}
impl From<SgsError> for ThermalSgsError {
    fn from(error: SgsError) -> Self {
        Self::Closure(error)
    }
}
impl From<ThermalSourceError> for ThermalSgsError {
    fn from(error: ThermalSourceError) -> Self {
        Self::Source(error)
    }
}

fn poll(gate: &CancelGate) -> Result<(), ThermalSgsError> {
    if gate.is_requested() {
        Err(ThermalSgsError::Cancelled)
    } else {
        Ok(())
    }
}

fn filled<T: Clone>(cells: usize, value: T) -> Result<Vec<T>, ThermalSgsError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(cells)
        .map_err(|_| ThermalSgsError::Allocation)?;
    values.resize(cells, value);
    Ok(values)
}

fn admit_velocity(rho: f64, u: [f64; 3], cell: usize, maximum: f64) -> Result<(), ThermalSgsError> {
    let mach_squared = 3.0 * u.iter().map(|v| v * v).sum::<f64>();
    if !rho.is_finite() || rho <= 0.0 || !mach_squared.is_finite() {
        return Err(ThermalSgsError::Invalid(
            "SGS requires finite physical velocity and positive density",
        ));
    }
    if mach_squared > maximum * maximum {
        return Err(ThermalSgsError::Mach {
            cell,
            mach_squared,
            maximum,
        });
    }
    Ok(())
}

impl ThermalLbm3 {
    /// Allocate SGS scratch once, subject to an explicit complete-cell budget.
    pub fn thermal_sgs_workspace(
        &self,
        max_cells: usize,
    ) -> Result<ThermalSgsWorkspace3, ThermalSgsError> {
        let (cells, _) = self.layout()?;
        if cells > max_cells {
            return Err(ThermalSgsError::Invalid(
                "SGS workspace cell budget exceeded",
            ));
        }
        Ok(ThermalSgsWorkspace3 {
            velocity: filled(cells, [0.0; 3])?,
            density: filled(cells, 0.0)?,
            momentum: filled(cells, 0.0)?,
            thermal: filled(cells, 0.0)?,
        })
    }

    /// Prepare total momentum AND heat transport using the shared SGS model.
    /// `molecular` optionally supplies independently authored local relaxation
    /// fields; otherwise the slab's original per-cell tau and scalar tau_g
    /// supply molecular transport. Never supply a previously augmented field.
    ///
    /// Flow uses `nu=(tau-1/2)/3`; D3Q7 uses `alpha=(tau_g-1/2)/4`.
    /// Zero modeled viscosity preserves both original relaxation bits exactly.
    /// Periodic derivatives are centered. At the first/last fluid layer the
    /// quadratic derivative uses the stationary wall at half a cell spacing
    /// and the adjacent fluid center, and is exact for a quadratic profile.
    /// Preparation is allocation-free with fixed cell order and cancellation
    /// polls every 64 cells. It never changes physical state.
    pub fn prepare_thermal_sgs(
        &self,
        config: ThermalSgsConfig,
        molecular: Option<RelaxationField<'_>>,
        workspace: &mut ThermalSgsWorkspace3,
        gate: &CancelGate,
    ) -> Result<ThermalSgsReport, ThermalSgsError> {
        poll(gate)?;
        config.validate()?;
        let (cells, plane) = self.layout()?;
        if [
            workspace.velocity.len(),
            workspace.density.len(),
            workspace.momentum.len(),
            workspace.thermal.len(),
        ]
        .iter()
        .any(|&n| n != cells)
        {
            return Err(ThermalSgsError::Invalid(
                "SGS workspace does not match the slab",
            ));
        }
        if !self.flow_valid() {
            return Err(ThermalSgsError::Invalid("invalid thermal/flow state"));
        }
        if let Some(field) = molecular {
            field.admit(cells)?;
        }
        for cell in 0..cells {
            if cell % 64 == 0 {
                poll(gate)?;
            }
            workspace.momentum[cell] =
                molecular.map_or(self.grid.tau[cell], |field| field.momentum[cell]);
            workspace.thermal[cell] = molecular.map_or(self.tau_g, |field| field.thermal[cell]);
            if cell < plane || cell >= cells - plane {
                workspace.velocity[cell] = [0.0; 3];
                workspace.density[cell] = 0.0;
            } else {
                let (rho, u) = self.grid.moments3(cell);
                admit_velocity(rho, u, cell, config.maximum_mach)?;
                workspace.velocity[cell] = u;
                workspace.density[cell] = rho;
            }
        }
        let mut report = ThermalSgsReport::default();
        for cell in plane..cells - plane {
            if cell % 64 == 0 {
                poll(gate)?;
            }
            let nu = (workspace.momentum[cell] - 0.5) / 3.0;
            let alpha = (workspace.thermal[cell] - 0.5) / 4.0;
            let gradient = self.thermal_sgs_gradient(cell, |i| workspace.velocity[i]);
            let eddy = eddy_viscosity(
                config.closure.model,
                config.closure.filter_width_lu,
                gradient,
            )?;
            let ratio = eddy / nu;
            if !ratio.is_finite() || ratio > config.closure.maximum_eddy_to_molecular {
                return Err(SgsError::ViscosityRatioExceeded {
                    cell,
                    ratio,
                    maximum: config.closure.maximum_eddy_to_molecular,
                }
                .into());
            }
            let diffusivity = eddy / config.turbulent_prandtl;
            if eddy != 0.0 {
                workspace.momentum[cell] = 0.5 + 3.0 * (nu + eddy);
                workspace.thermal[cell] = 0.5 + 4.0 * (alpha + diffusivity);
                if !diffusivity.is_finite()
                    || diffusivity <= 0.0
                    || !workspace.momentum[cell].is_finite()
                    || !workspace.thermal[cell].is_finite()
                {
                    return Err(ThermalSgsError::Invalid(
                        "total SGS transport is not representable",
                    ));
                }
            }
            let dilation_third = (gradient[0][0] + gradient[1][1] + gradient[2][2]) / 3.0;
            let mut strain2 = 0.0;
            for i in 0..3 {
                for j in 0..3 {
                    let s = 0.5 * (gradient[i][j] + gradient[j][i])
                        - if i == j { dilation_third } else { 0.0 };
                    strain2 += s * s;
                }
            }
            let factor = 2.0 * workspace.density[cell] * strain2;
            report.momentum.fluid_cells += 1;
            report.momentum.resolved_dissipation_lu += factor * nu;
            report.momentum.sgs_dissipation_lu += factor * eddy;
            report.momentum.maximum_eddy_viscosity_lu =
                report.momentum.maximum_eddy_viscosity_lu.max(eddy);
            report.momentum.maximum_eddy_to_molecular =
                report.momentum.maximum_eddy_to_molecular.max(ratio);
            report.maximum_eddy_diffusivity_lu =
                report.maximum_eddy_diffusivity_lu.max(diffusivity);
        }
        poll(gate)?;
        if !report.momentum.resolved_dissipation_lu.is_finite()
            || !report.momentum.sgs_dissipation_lu.is_finite()
        {
            return Err(SgsError::NonFiniteViscosity.into());
        }
        Ok(report)
    }

    fn thermal_sgs_gradient(
        &self,
        cell: usize,
        velocity: impl Fn(usize) -> [f64; 3],
    ) -> [[f64; 3]; 3] {
        let g = &self.grid;
        let plane = g.nx * g.ny;
        let (x, y, z) = (cell % g.nx, (cell / g.nx) % g.ny, cell / plane);
        let xm = g.idx(if x == 0 { g.nx - 1 } else { x - 1 }, y, z);
        let xp = g.idx(if x + 1 == g.nx { 0 } else { x + 1 }, y, z);
        let ym = g.idx(x, if y == 0 { g.ny - 1 } else { y - 1 }, z);
        let yp = g.idx(x, if y + 1 == g.ny { 0 } else { y + 1 }, z);
        core::array::from_fn(|component| {
            [
                0.5 * (velocity(xp)[component] - velocity(xm)[component]),
                0.5 * (velocity(yp)[component] - velocity(ym)[component]),
                if z == 1 {
                    velocity(cell)[component] + velocity(cell + plane)[component] / 3.0
                } else if z + 2 == g.nz {
                    -velocity(cell)[component] - velocity(cell - plane)[component] / 3.0
                } else {
                    0.5 * (velocity(cell + plane)[component] - velocity(cell - plane)[component])
                },
            ]
        })
    }

    /// Advance the original source/flow/heat owner with matching turbulent
    /// momentum and temperature diffusion. Heat inputs are prescribed lattice
    /// temperature increments; forces are external force densities, excluding
    /// buoyancy. All inputs include the two wall planes, whose heating is zero.
    ///
    /// The original source, temperature and independent wall-balance checks,
    /// plus endpoint Mach/cancellation admission, precede publication. Neither
    /// total rate persists into later ordinary steps. Scalar `rayleigh`,
    /// `nusselt` and `diffusivity` still describe the molecular setup; the
    /// returned wall exchanges describe the actual turbulent heat step.
    /// Cancellation is checked during preparation and around the unchanged
    /// serial thermal kernel, not within each native collision/stream update.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_sources_and_sgs(
        &mut self,
        heat_per_step: &[f64],
        forces: &[[f64; 3]],
        config: ThermalSgsConfig,
        molecular: Option<RelaxationField<'_>>,
        workspace: &mut ThermalSgsWorkspace3,
        limits: ThermalSourceLimits,
        gate: &CancelGate,
    ) -> Result<ThermalSgsStep, ThermalSgsError> {
        let (heat, sgs) = step_sgs_mapped(
            self,
            heat_per_step,
            forces,
            config,
            molecular,
            workspace,
            limits,
            gate,
            Ok,
        )?;
        Ok(ThermalSgsStep { heat, sgs })
    }
}

// Physical-unit adapters keep their original balance checks inside the SAME
// transaction. The side channel retains the precise SGS endpoint error while
// the original source owner refuses its candidate before publication.
#[allow(clippy::too_many_arguments)]
pub(crate) fn step_sgs_mapped<R>(
    state: &mut ThermalLbm3,
    heat: &[f64],
    forces: &[[f64; 3]],
    config: ThermalSgsConfig,
    molecular: Option<RelaxationField<'_>>,
    workspace: &mut ThermalSgsWorkspace3,
    limits: ThermalSourceLimits,
    gate: &CancelGate,
    finish: impl FnOnce(ThermalSourceBalance) -> Result<R, ThermalSourceError>,
) -> Result<(R, ThermalSgsReport), ThermalSgsError> {
    let report = state.prepare_thermal_sgs(config, molecular, workspace, gate)?;
    let mut admission = Ok(());
    let result = source_step_transport_mapped::<ThermalLbm3, 3, 7, _>(
        state,
        heat,
        forces,
        limits,
        Some(workspace.relaxation()),
        |candidate, balance| {
            admission = (|| {
                let (cells, plane) = candidate.layout()?;
                for cell in plane..cells - plane {
                    if cell % 64 == 0 {
                        poll(gate)?;
                    }
                    let (rho, u) = candidate.grid.moments3(cell);
                    admit_velocity(rho, u, cell, config.maximum_mach)?;
                }
                poll(gate)
            })();
            if admission.is_err() {
                return Err(ThermalSourceError::Invalid("SGS endpoint admission failed"));
            }
            finish(balance)
        },
    );
    admission?;
    Ok((result?, report))
}

#[cfg(test)]
mod tests;
