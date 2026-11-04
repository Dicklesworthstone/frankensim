//! Explicit full-gradient closures; the original stress Smagorinsky is intact.

use super::*;
use crate::d3q19::freesurface3::{
    FreeSurfaceKinematicsError3, FreeSurfaceKinematicsWorkspace3, WetGradientStencil3,
};

impl From<FreeSurfaceKinematicsError3> for FreeSurfaceSgsError3 {
    fn from(error: FreeSurfaceKinematicsError3) -> Self {
        Self::Kinematics(error)
    }
}

impl FreeSurface3 {
    /// Prepare WALE or Vreman from the FULL force-centered liquid gradient,
    /// including rotation. The separately budgeted velocity cache is reused.
    /// Occupied liquid cells need resolved second-order stencils; under-resolved
    /// liquid sheets refuse the whole preparation. Empty interface ghosts keep
    /// their native molecular relaxation, without an invented gradient or liquid
    /// modeled dissipation. No gas velocity, wall velocity or curl is invented.
    /// BGK/TRT and its authored TRT even/odd product remain the native family.
    pub fn prepare_gradient_sgs(
        &self,
        gate: &CancelGate,
        config: SgsConfig,
        kinematics: &mut FreeSurfaceKinematicsWorkspace3,
        workspace: &mut FreeSurfaceSgsWorkspace3,
    ) -> Result<FreeSurfaceSgsReport3, FreeSurfaceSgsError3> {
        check_cancel(gate)?;
        config.validate()?;
        if !matches!(
            config.model,
            SgsModel::Wale { .. } | SgsModel::Vreman { .. }
        ) {
            return Err(FreeSurfaceSgsError3::Invalid(
                "gradient SGS requires WALE or Vreman; use prepare_sgs for stress Smagorinsky",
            ));
        }
        let model = self.collision_model();
        let molecular_tau = match model {
            CollisionModel3::Bgk { tau } | CollisionModel3::Trt { tau, .. } => tau,
            _ => {
                return Err(FreeSurfaceSgsError3::Invalid(
                    "gradient SGS requires native BGK or TRT collision",
                ));
            }
        };
        let molecular_nu = (molecular_tau - 0.5) / 3.0;
        if !molecular_nu.is_finite() || molecular_nu <= 0.0 {
            return Err(SgsError::InvalidMolecularViscosity { cell: 0 }.into());
        }
        let field = self.kinematics(gate, kinematics)?;
        workspace.tau.retain(|&(key, lane), _| {
            self.tiles
                .get(&key)
                .is_some_and(|t| t.cells[lane] != Cell3::Gas)
        });
        let mut report = FreeSurfaceSgsReport3::default();
        for (&key, tile) in &self.tiles {
            check_cancel(gate)?;
            let slot = self.grid.slot_of(key).expect("active");
            for lane in 0..TILE_CELLS {
                if tile.cells[lane] == Cell3::Gas {
                    continue;
                }
                let (x, y, z) = Self::global_of(key, lane);
                let coordinates = [x as usize, y as usize, z as usize];
                let invalid = |reason| FreeSurfaceSgsError3::Cell {
                    coordinates,
                    reason,
                };
                let fill = self.fill(key, lane);
                if tile.cells[lane] == Cell3::Interface && fill == 0.0 {
                    // Empty interfaces still carry populations for native gas
                    // reconstruction and mass exchange. Preserve their original
                    // collision law; no liquid volume supports a modeled eddy.
                    workspace.tau.insert((key, lane), molecular_tau);
                    report.wet_cells += 1;
                    continue;
                }
                let sample = field
                    .sample_at(coordinates[0], coordinates[1], coordinates[2])?
                    .ok_or_else(|| {
                        invalid("gradient SGS requires a liquid sample at every occupied cell")
                    })?;
                if sample.stencils.iter().any(|s| {
                    matches!(
                        s,
                        WetGradientStencil3::ForwardFirstOrder
                            | WetGradientStencil3::BackwardFirstOrder
                    )
                }) {
                    return Err(invalid(
                        "gradient SGS requires second-order liquid stencils",
                    ));
                }
                let gradient = sample.kinematics.gradient();
                let eddy = eddy_viscosity(config.model, config.filter_width_lu, gradient)?;
                let ratio = eddy / molecular_nu;
                if !ratio.is_finite() || ratio > config.maximum_eddy_to_molecular {
                    return Err(SgsError::ViscosityRatioExceeded {
                        cell: report.wet_cells,
                        ratio,
                        maximum: config.maximum_eddy_to_molecular,
                    }
                    .into());
                }
                let tau = if eddy == 0.0 {
                    molecular_tau
                } else {
                    molecular_tau + 3.0 * eddy
                };
                D3q19::with_shear_relaxation(model, tau).map_err(|source| {
                    FreeSurfaceError3::Collision {
                        tile_key: key,
                        lane,
                        source,
                    }
                })?;
                workspace.tau.insert((key, lane), tau);
                // The closures see rotation and dilation; reported power is
                // explicitly deviatoric shear, matching the original SGS path.
                let third_trace =
                    gradient[0][0] / 3.0 + gradient[1][1] / 3.0 + gradient[2][2] / 3.0;
                let mut rate_squared = 0.0;
                for i in 0..3 {
                    for j in 0..3 {
                        let strain = 0.5 * gradient[i][j] + 0.5 * gradient[j][i]
                            - if i == j { third_trace } else { 0.0 };
                        rate_squared += 2.0 * strain * strain;
                    }
                }
                let (rho, _) = self.grid.cell_macros(slot, lane);
                let power = fill * rho * rate_squared;
                if !power.is_finite() {
                    return Err(invalid("unrepresentable liquid shear dissipation"));
                }
                report.wet_cells += 1;
                report.liquid_volume_lu += fill;
                report.resolved_dissipation_lu += power * molecular_nu;
                report.sgs_dissipation_lu += power * eddy;
                report.maximum_eddy_viscosity_lu = report.maximum_eddy_viscosity_lu.max(eddy);
                report.maximum_eddy_to_molecular = report.maximum_eddy_to_molecular.max(ratio);
            }
        }
        check_cancel(gate)?;
        if !report.resolved_dissipation_lu.is_finite() || !report.sgs_dissipation_lu.is_finite() {
            return Err(FreeSurfaceSgsError3::Invalid(
                "SGS dissipation sum overflow",
            ));
        }
        Ok(report)
    }

    /// Advance native collision, gas reconstruction, VOF mass and conversion
    /// using explicit full-gradient WALE/Vreman viscosity. One native candidate
    /// commits only after every closure, boundary, conversion and cancel check.
    /// Zero coefficient preserves original step bits on admitted stencils.
    pub fn step_with_gradient_sgs(
        &mut self,
        gate: &CancelGate,
        config: SgsConfig,
        kinematics: &mut FreeSurfaceKinematicsWorkspace3,
        workspace: &mut FreeSurfaceSgsWorkspace3,
    ) -> Result<FreeSurfaceSgsReport3, FreeSurfaceSgsError3> {
        let report = self.prepare_gradient_sgs(gate, config, kinematics, workspace)?;
        let mut candidate = self.clone();
        if config.model.coefficient() == 0.0 {
            candidate.step_inner(gate, None, None, None)?;
        } else {
            candidate.step_with_material_relaxation_and_pressure(gate, &workspace.tau, None)?;
        }
        check_cancel(gate)?;
        *self = candidate;
        Ok(report)
    }
}
