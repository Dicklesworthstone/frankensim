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

    /// Initialize liquid fill and velocity with a physical pressure field [Pa].
    ///
    /// Requires [`Self::with_pressure_reference`]. `pressure_pa` supplies the
    /// positive absolute LIQUID pressure at each occupied cell center [m],
    /// including any hydrostatic or capillary jump. The field is sampled once
    /// per supplied cell; it is not retained or imposed during later steps.
    /// Subsequent gas loading uses [`MaterialFreeSurface3::step_with_gas_pressure`]
    /// and the same fixed reference. A common shift of both reference and
    /// pressure leaves the dimensionless initial state unchanged.
    ///
    /// Density is `1 + (p - p_reference)/(rho_source*(dx/dt)^2*cs^2)`.
    /// Its excursion from one must not exceed the caller's finite bound in
    /// `(0,1)`. Native mass is therefore `fill*density*rho_source*dx^3`,
    /// including the weakly compressible numerical density perturbation;
    /// this is not a physical liquid equation of state or a mass-preserving
    /// pressure projection. The source properties stay frozen. The caller
    /// owns initial pressure/boundary compatibility and grid/time convergence;
    /// equilibrium population initialization is not exactly well-balanced.
    ///
    /// # Errors
    /// The ordinary initial-state and source gates apply. Missing reference,
    /// nonpositive/nonfinite pressure, unrepresentable coordinates/conversion,
    /// excessive density excursion, or cancellation refuse the whole candidate.
    pub fn free_surface3_from_cells_with_pressure(
        self,
        domain: MaterialFreeSurfaceDomain3<'_>,
        initial: &[MaterialFreeSurfaceInitialCell3],
        maximum_density_excursion: f64,
        pressure_pa: impl Fn([f64; 3]) -> f64,
        cx: &Cx<'_>,
    ) -> Result<MaterialFreeSurface3<'a>, MaterialLbmError> {
        cx.checkpoint()
            .map_err(|_| crate::d3q19::freesurface3::FreeSurfaceError3::Cancelled)?;
        self.require_frozen_source()?;
        let reference = self.reference_pressure_pa.ok_or(MaterialLbmError::Invalid(
            "initial absolute liquid pressure requires an explicit fixed pressure reference",
        ))?;
        if !maximum_density_excursion.is_finite()
            || maximum_density_excursion <= 0.0
            || maximum_density_excursion >= 1.0
        {
            return Err(MaterialLbmError::Invalid(
                "initial density-excursion bound must be finite and in (0,1)",
            ));
        }
        let dx = self.cell_size_m;
        let pressure_unit = self.pressure_unit_pa;
        self.free_surface3_from_cells_with_density(domain, initial, cx, |coordinates| {
            let position = coordinates.map(|c| (c as f64 + 0.5) * dx);
            let absolute = pressure_pa(position);
            let perturbation = (absolute - reference) / pressure_unit;
            let density = 1.0 + perturbation / CS2;
            if position.iter().any(|x| !x.is_finite())
                || !absolute.is_finite()
                || absolute <= 0.0
                || !density.is_finite()
                || density <= 0.0
                || (density - 1.0).abs() > maximum_density_excursion
            {
                return Err(
                    crate::d3q19::freesurface3::FreeSurfaceError3::InvalidInitialState {
                        coordinates: Some(coordinates),
                        reason: "initial liquid pressure or its density excursion is inadmissible",
                    }
                    .into(),
                );
            }
            Ok(density)
        })
    }

    pub(super) fn free_surface3_from_cells_inner(
        self,
        domain: MaterialFreeSurfaceDomain3<'_>,
        initial: &[MaterialFreeSurfaceInitialCell3],
        cx: &Cx<'_>,
    ) -> Result<MaterialFreeSurface3<'a>, MaterialLbmError> {
        self.free_surface3_from_cells_with_density(domain, initial, cx, |_| Ok(1.0))
    }

    fn free_surface3_from_cells_with_density(
        self,
        domain: MaterialFreeSurfaceDomain3<'_>,
        initial: &[MaterialFreeSurfaceInitialCell3],
        cx: &Cx<'_>,
        density: impl Fn([usize; 3]) -> Result<f64, MaterialLbmError>,
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
                density: density(cell.coordinates)?,
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
