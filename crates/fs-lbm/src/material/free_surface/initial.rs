//! Supplied physical mass and velocity enter the existing material-bound VOF.

use super::*;
use crate::d3q19::freesurface3::FreeSurfaceInitialCell3;
use fs_exec::Cx;

/// One occupied cell of the same frozen pure liquid as the receiving solver.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialFreeSurfaceInitialCell3 {
    /// Global Cartesian cell coordinates.
    pub coordinates: [usize; 3],
    /// Liquid volume fraction in (0,1]. Full boundary cells retain fill one.
    pub fill_fraction: f64,
    /// Physical initial velocity [m/s], not a population velocity before forcing.
    pub velocity_m_s: [f64; 3],
}

impl<'a> MaterialFreeSurfaceScaling<'a> {
    /// Construct native flow from transferred fill and physical velocity.
    ///
    /// Mass is `fill*rho_source*dx³`; numerical density is initialized at the
    /// source's reference density (zero incompressible pressure perturbation).
    /// This represents one frozen, pure liquid, not a density-changing phase
    /// mixture. The caller owns remap geometry, subcell kinetic-energy loss,
    /// initial pressure compatibility and any discarded constitutive stress.
    /// No accepted trajectory is mutated and no thermal state is inferred.
    ///
    /// # Errors
    /// Refuses the existing domain/material gates, initial speeds above the
    /// admitted scenario bound, native fill/state/population errors, and
    /// cancellation before returning the complete candidate.
    pub fn free_surface3_from_cells(
        self,
        domain: MaterialFreeSurfaceDomain3<'_>,
        initial: &[MaterialFreeSurfaceInitialCell3],
        cx: &Cx<'_>,
    ) -> Result<MaterialFreeSurface3<'a>, MaterialLbmError> {
        self.require_frozen_source()?;
        self.free_surface3_from_cells_inner(domain, initial, cx)
    }

    pub(super) fn free_surface3_from_cells_inner(
        self,
        domain: MaterialFreeSurfaceDomain3<'_>,
        initial: &[MaterialFreeSurfaceInitialCell3],
        cx: &Cx<'_>,
    ) -> Result<MaterialFreeSurface3<'a>, MaterialLbmError> {
        cx.checkpoint()
            .map_err(|_| crate::d3q19::freesurface3::FreeSurfaceError3::Cancelled)?;
        if self.planar_wetting.is_some() && domain.contact != ContactModel3::Neutral {
            return Err(MaterialLbmError::Invalid(
                "a sourced planar angle requires Neutral on all other walls",
            ));
        }
        // Geometry admission is independent of the occupied-state list. The
        // native constructor checks every supplied cell against these walls.
        admit_domain(&domain, &|_, _, _| false)?;
        if initial.len() > domain.max_cells {
            return Err(MaterialLbmError::Invalid(
                "initial liquid exceeds the domain cell budget",
            ));
        }
        let mut gravity = [0.0; 3];
        for (component, physical) in gravity.iter_mut().zip(domain.gravity_m_s2) {
            *component = self.gravity_lattice(physical)?;
        }
        let mut native = Vec::with_capacity(initial.len());
        for cell in initial {
            if cx.checkpoint().is_err() {
                return Err(crate::d3q19::freesurface3::FreeSurfaceError3::Cancelled.into());
            }
            let velocity = cell.velocity_m_s.map(|v| v / self.velocity_scale);
            let speed = velocity[0].hypot(velocity[1]).hypot(velocity[2]);
            if !speed.is_finite() || speed > self.flow.u_lattice {
                return Err(MaterialLbmError::Invalid(
                    "initial free-surface velocity exceeds the admitted scenario speed",
                ));
            }
            native.push(FreeSurfaceInitialCell3 {
                coordinates: cell.coordinates,
                fill_fraction: cell.fill_fraction,
                density: 1.0,
                velocity,
            });
        }
        let mut lattice = FreeSurface3::from_cells(
            domain.cells,
            self.flow.tau,
            gravity,
            self.surface_tension_lattice,
            domain.contact,
            domain.solid_tiles,
            &native,
            cx,
        )?;
        if let Some((wall, wetting)) = self.planar_wetting {
            lattice.set_planar_contact_angle(wall, wetting.contact_angle_radians())?;
        }
        cx.checkpoint()
            .map_err(|_| crate::d3q19::freesurface3::FreeSurfaceError3::Cancelled)?;
        Ok(MaterialFreeSurface3 {
            lattice,
            scaling: self,
        })
    }
}
