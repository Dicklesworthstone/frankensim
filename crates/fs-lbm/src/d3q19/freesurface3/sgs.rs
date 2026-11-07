//! Local Smagorinsky transport on the native sparse free surface.
//!
//! The force-corrected, deviatoric non-equilibrium stress gives `tau S_dev`.
//! Its Smagorinsky viscosity is `K/tau`, so the physical root of
//! `tau (tau - tau_molecular) = 3 K` supplies the current collision, without
//! a lagged viscosity or gradients through empty gas cells. The original VOF
//! owner still performs gas reconstruction, mass exchange and phase conversion.
//!
//! The quadratic follows the Smagorinsky BGK construction in Chiappini,
//! Water 12 (2020), 1212, equations (6)-(8), doi:10.3390/w12041212. Here the
//! trace is removed before fitting shear, and the original Guo half-source
//! correction is retained. BGK/TRT retain their tied shear/bulk relaxation;
//! TRT preserves the authored even/odd product as the even viscosity changes.
//! Independent central-moment collision instead retains its molecular bulk
//! rate at every cell while adding modeled viscosity only to deviatoric shear.
//! This is weakly compressible, isothermal exploration:
//! no subgrid free-surface, entrainment, wall damping, resolved heating, or
//! turbulent-flow validation is supplied. The separate `prepare_gradient_sgs`
//! and `step_with_gradient_sgs` path admits WALE/Vreman from the native full
//! liquid velocity gradient; this symmetric stress path remains Smagorinsky.

use std::collections::BTreeMap;

use super::transport::MassLocation3;
use super::{Cell3, FreeSurface3, FreeSurfaceError3, TILE_CELLS, check_cancel};
use crate::d3q19::{CollisionModel3, force_corrected_stress3};
use crate::lattice3::{D3q19, ViscousLattice3};
use crate::sgs3::{SgsConfig, SgsError, SgsModel, eddy_viscosity};
use fs_exec::CancelGate;
use fs_math::det;

mod gradient;
mod motion;

#[derive(Clone, Debug, PartialEq)]
pub enum FreeSurfaceSgsError3 {
    Flow(FreeSurfaceError3),
    Closure(SgsError),
    Kinematics(super::FreeSurfaceKinematicsError3),
    Invalid(&'static str),
    Cell {
        coordinates: [usize; 3],
        reason: &'static str,
    },
}

impl core::fmt::Display for FreeSurfaceSgsError3 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-LBM-FREE-SURFACE-SGS3: {self:?}")
    }
}
impl std::error::Error for FreeSurfaceSgsError3 {}
impl From<FreeSurfaceError3> for FreeSurfaceSgsError3 {
    fn from(error: FreeSurfaceError3) -> Self {
        Self::Flow(error)
    }
}
impl From<SgsError> for FreeSurfaceSgsError3 {
    fn from(error: SgsError) -> Self {
        Self::Closure(error)
    }
}

/// Instantaneous pre-step diagnostics in lattice units. Dissipation is weighted
/// by the native bounded liquid fill, so empty interface ghosts and unlocated
/// conversion carry contribute no liquid dissipation. This is not heat added
/// to the isothermal solver or a discrete energy-conservation certificate.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct FreeSurfaceSgsReport3 {
    /// All collision-bearing fluid/interface cells, including empty ghosts.
    pub wet_cells: usize,
    pub liquid_volume_lu: f64,
    /// Molecular deviatoric shear power `2 rho nu D_dev:D_dev`.
    pub resolved_dissipation_lu: f64,
    /// Molecular volumetric power `zeta (div u)^2`, separate from shear.
    /// BGK/TRT's additional tied eddy bulk damping is not molecular power.
    pub resolved_bulk_dissipation_lu: f64,
    /// Modeled deviatoric power only; no eddy bulk viscosity is authored.
    pub sgs_dissipation_lu: f64,
    pub maximum_eddy_viscosity_lu: f64,
    pub maximum_eddy_to_molecular: f64,
}

/// Reused native relaxation cache, keyed by Morton tile and lane, never sparse
/// slot. Existing wet cells reuse their entries; newly wet cells allocate an
/// entry on the next preparation. Outputs are consumable only after success.
#[derive(Default)]
pub struct FreeSurfaceSgsWorkspace3 {
    tau: BTreeMap<(u64, usize), f64>,
}

impl FreeSurfaceSgsWorkspace3 {
    /// Relaxation in the current scratch cache, valid only after preparation
    /// succeeds. A refused preparation may leave partial scratch. After success,
    /// gas, inactive and out-of-domain cells return None. A step's report/cache
    /// describe its input state; the next preparation admits newly wet cells.
    pub fn relaxation_at(&self, coordinates: [usize; 3]) -> Option<f64> {
        let MassLocation3::Cell(key, lane) = FreeSurface3::transport_location(coordinates)? else {
            return None;
        };
        self.tau.get(&(key, lane)).copied()
    }
}

/// Cancellation and closure refusal leave the native accepted state untouched.
impl FreeSurface3 {
    /// Prepare an explicitly authored Smagorinsky closure at every wet cell.
    /// The native molecular law must be BGK, TRT or independent central moment;
    /// its family, TRT even/odd product and independent molecular bulk/higher
    /// rates are retained. No collision family is replaced. A bound bulk field
    /// must supply every wet cell, including empty interface carriers.
    /// Coefficient zero retains the exact original relaxation. The viscosity
    /// ratio is refused, never clipped, even in a zero-mass interface ghost.
    pub fn prepare_sgs(
        &self,
        gate: &CancelGate,
        config: SgsConfig,
        workspace: &mut FreeSurfaceSgsWorkspace3,
    ) -> Result<FreeSurfaceSgsReport3, FreeSurfaceSgsError3> {
        check_cancel(gate)?;
        config.validate()?;
        if !matches!(config.model, SgsModel::Smagorinsky { .. }) {
            return Err(FreeSurfaceSgsError3::Invalid(
                "free-surface stress SGS supports Smagorinsky only",
            ));
        }
        let model = self.collision_model();
        let molecular_tau = molecular_relaxation(model)?;
        let molecular_nu = (molecular_tau - 0.5) / 3.0;
        if !molecular_nu.is_finite() || molecular_nu <= 0.0 {
            return Err(SgsError::InvalidMolecularViscosity { cell: 0 }.into());
        }
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
                let local_model = self.collision_with_bulk(key, lane, model)?;
                let (x, y, z) = Self::global_of(key, lane);
                let invalid = |reason| FreeSurfaceSgsError3::Cell {
                    coordinates: [x as usize, y as usize, z as usize],
                    reason,
                };
                let f = self.grid.populations(slot, lane);
                let (rho, raw) = self.grid.cell_macros(slot, lane);
                let force = self.force_at_lane(key, lane);
                let u = core::array::from_fn(|a| raw[a] + 0.5 * force[a] / rho);
                if !rho.is_finite()
                    || rho <= 0.0
                    || f.iter().chain(&u).chain(&force).any(|v| !v.is_finite())
                {
                    return Err(invalid(
                        "positive finite density and source-frame state required",
                    ));
                }
                let [xx, yy, zz, xy, xz, yz] = force_corrected_stress3(&f, rho, u, force);
                let trace = xx / 3.0 + yy / 3.0 + zz / 3.0;
                let strain_tau = [
                    [xx - trace, xy, xz],
                    [xy, yy - trace, yz],
                    [xz, yz, zz - trace],
                ]
                .map(|row| row.map(|v| -1.5 * (v / rho)));
                // The shared closure returns K = tau * nu_t from tau*S_dev.
                let k = eddy_viscosity(config.model, config.filter_width_lu, strain_tau)?;
                let (tau, eddy) = relaxation(molecular_tau, k)
                    .ok_or_else(|| invalid("unrepresentable self-consistent SGS relaxation"))?;
                let ratio = eddy / molecular_nu;
                if !ratio.is_finite() || ratio > config.maximum_eddy_to_molecular {
                    return Err(SgsError::ViscosityRatioExceeded {
                        cell: report.wet_cells,
                        ratio,
                        maximum: config.maximum_eddy_to_molecular,
                    }
                    .into());
                }
                // A finite even tau can still round TRT's derived odd tau to
                // one half. Admit the exact same model the native owner uses.
                D3q19::with_shear_relaxation(local_model, tau).map_err(|source| {
                    FreeSurfaceError3::Collision {
                        tile_key: key,
                        lane,
                        source,
                    }
                })?;
                workspace.tau.insert((key, lane), tau);
                let rate_squared: f64 = strain_tau
                    .iter()
                    .flatten()
                    .map(|s| 2.0 * (s / tau).powi(2))
                    .sum();
                let fill = self.fill(key, lane);
                let shear_power = fill * rho * rate_squared;
                // Guo-corrected spherical stress carries tau_bulk*div(u).
                // An independent local bulk rate is never retuned by SGS.
                let bulk_rate = match local_model {
                    CollisionModel3::CentralMomentBulk { bulk_rate, .. } => bulk_rate,
                    _ => 1.0 / tau,
                };
                let divergence = (-4.5 * (trace / rho)) * bulk_rate;
                let bulk_power = (fill * rho * local_model.kinematic_bulk_viscosity())
                    * divergence * divergence;
                if !shear_power.is_finite() || !bulk_power.is_finite() {
                    return Err(invalid("unrepresentable liquid viscous dissipation"));
                }
                report.wet_cells += 1;
                report.liquid_volume_lu += fill;
                report.resolved_dissipation_lu += shear_power * molecular_nu;
                report.resolved_bulk_dissipation_lu += bulk_power;
                report.sgs_dissipation_lu += shear_power * eddy;
                report.maximum_eddy_viscosity_lu = report.maximum_eddy_viscosity_lu.max(eddy);
                report.maximum_eddy_to_molecular = report.maximum_eddy_to_molecular.max(ratio);
            }
        }
        check_cancel(gate)?;
        if !report.resolved_dissipation_lu.is_finite()
            || !report.resolved_bulk_dissipation_lu.is_finite()
            || !report.sgs_dissipation_lu.is_finite()
        {
            return Err(FreeSurfaceSgsError3::Invalid(
                "SGS dissipation sum overflow",
            ));
        }
        Ok(report)
    }

    /// Advance actual free-surface populations, mass and conversion together.
    /// A single private native candidate preserves all accepted state on a
    /// collision, conversion, closure or cancellation refusal. The zero closure
    /// uses the unmodified native collision path and preserves its trajectory.
    pub fn step_with_sgs(
        &mut self,
        gate: &CancelGate,
        config: SgsConfig,
        workspace: &mut FreeSurfaceSgsWorkspace3,
    ) -> Result<FreeSurfaceSgsReport3, FreeSurfaceSgsError3> {
        let (candidate, report) = self.sgs_candidate(gate, config, workspace)?;
        check_cancel(gate)?;
        *self = candidate;
        Ok(report)
    }

    /// Propose the same native SGS step without publishing it. Physical-unit
    /// owners can admit the returned state and diagnostics without a second
    /// population clone. The report and workspace describe the input state.
    /// If a local bulk field is bound, a wet-topology change requires source
    /// rebinding and refuses here; this isothermal native path cannot supply it.
    pub fn sgs_candidate(
        &self,
        gate: &CancelGate,
        config: SgsConfig,
        workspace: &mut FreeSurfaceSgsWorkspace3,
    ) -> Result<(Self, FreeSurfaceSgsReport3), FreeSurfaceSgsError3> {
        let report = self.prepare_sgs(gate, config, workspace)?;
        let mut candidate = self.clone();
        if config.model.coefficient() == 0.0 {
            candidate.step_inner(gate, None, None, None)?;
        } else {
            candidate.step_with_material_relaxation_and_pressure(gate, &workspace.tau, None)?;
        }
        candidate.admit_sgs_bulk_field(gate)?;
        check_cancel(gate)?;
        Ok((candidate, report))
    }

    fn admit_sgs_bulk_field(&self, gate: &CancelGate) -> Result<(), FreeSurfaceSgsError3> {
        if let Some(field) = &self.local_bulk_relaxation {
            self.admit_bulk_relaxations(field, gate)?;
        }
        Ok(())
    }

    /// Allocation-free endpoint gate for the original physical scenario speed.
    pub(crate) fn admit_sgs_speed(
        &self,
        maximum: f64,
        gate: &CancelGate,
    ) -> Result<(), FreeSurfaceSgsError3> {
        check_cancel(gate)?;
        for (&key, tile) in &self.tiles {
            check_cancel(gate)?;
            let slot = self.grid.slot_of(key).expect("active");
            for lane in 0..TILE_CELLS {
                if tile.cells[lane] == Cell3::Gas {
                    continue;
                }
                let (rho, raw) = self.grid.cell_macros(slot, lane);
                let force = self.force_at_lane(key, lane);
                let velocity: [f64; 3] = core::array::from_fn(|a| raw[a] + 0.5 * force[a] / rho);
                let speed = velocity[0].hypot(velocity[1]).hypot(velocity[2]);
                if !speed.is_finite() || speed > maximum {
                    return Err(FreeSurfaceSgsError3::Invalid(
                        "SGS liquid exceeds the original material scenario speed",
                    ));
                }
            }
        }
        check_cancel(gate)?;
        Ok(())
    }
}

fn molecular_relaxation(model: CollisionModel3) -> Result<f64, FreeSurfaceSgsError3> {
    match model {
        CollisionModel3::Bgk { tau } | CollisionModel3::Trt { tau, .. } => Ok(tau),
        CollisionModel3::CentralMomentBulk { shear_rate, .. } => Ok(1.0 / shear_rate),
        _ => Err(FreeSurfaceSgsError3::Invalid(
            "free-surface SGS requires native BGK, TRT or independent central-moment collision",
        )),
    }
}

// Rationalized positive root: no cancellation at small modeled viscosity and
// no tau^2 or 12*K overflow before the square root. K has viscosity*tau units.
fn relaxation(tau: f64, k: f64) -> Option<(f64, f64)> {
    if k == 0.0 {
        return Some((tau, 0.0));
    }
    let stress_scale = 2.0 * det::sqrt(3.0) * det::sqrt(k);
    let scale = tau.max(stress_scale);
    let root = scale * det::sqrt((tau / scale).powi(2) + (stress_scale / scale).powi(2));
    let increment = 3.0 * (k / (0.5 * root + 0.5 * tau));
    let effective = tau + increment;
    let eddy = k / effective;
    (effective.is_finite() && effective > 0.5 && eddy.is_finite()).then_some((effective, eddy))
}

#[cfg(test)]
#[path = "sgs/tests.rs"]
mod tests;
