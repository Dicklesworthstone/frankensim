//! Original finite true-absorption cells on the native scattering flights.
//!
//! The grid is scalar, so every Stokes component keeps the same attenuation.
//! Its original constituent source supplies unpolarized LTE only to reciprocal
//! camera paths. Scattering survival is already sampled by the flight owner.

use super::*;
use crate::charts::Ray;
use crate::polarization::Stokes;
use crate::polarization::dielectric::transport::BeamError;
use crate::volumes::transport::GridSegmentTransport;

impl<'a, 'scene, 'optics> CurvedRayleighScene<'a, 'scene, 'optics> {
    /// Bind finite density, true absorption and optional LTE from the original
    /// native cells. A nonempty table has one entry per original scene medium.
    /// Borrowed grids retain their exact palette, source and snapshot identity.
    ///
    /// Each bound host must be homogeneous; at every transported wavelength it
    /// must also be index-one and clear, as required by `AbsorbingGrid`. The grid
    /// owns its absorption, independently of the existing scattering law. A
    /// camera cannot also bind an isothermal/spatial source to that medium.
    ///
    /// The cumulative cell cap applies to one camera/power path or one finite
    /// source replay. Each segment also retains the grid's own native cap. Zero
    /// permits misses/vacuum only. Proved escape tails clip the finite grid and
    /// are counted; no field is extrapolated beyond its native support.
    pub fn with_absorbing_grids(
        mut self,
        cx: &Cx<'_>,
        grids: &'a [Option<&'a AbsorbingGrid<'optics>>],
        max_cell_visits_per_path: u64,
    ) -> Result<Self, CurvedSceneError> {
        cx.checkpoint()?;
        if !grids.is_empty() && grids.len() != self.scene.media.len() {
            return Err(SceneError::InvalidBindings.into());
        }
        for (index, grid) in grids.iter().enumerate() {
            if index % 1024 == 0 {
                cx.checkpoint()?;
            }
            if grid.is_some() && !matches!(self.scene.media[index], SceneMedium::Homogeneous(_)) {
                return Err(scattering_error(PolarizedScatteringError::CarrierMedium));
            }
        }
        self.absorption_grids = grids;
        self.max_absorption_cells = max_cell_visits_per_path;
        Ok(self)
    }

    /// Exact original finite absorption/source owner, if explicitly bound.
    pub fn absorbing_grid(&self, medium: usize) -> Option<&AbsorbingGrid<'optics>> {
        self.absorption_grids.get(medium).copied().flatten()
    }

    pub const fn max_absorption_cell_visits(&self) -> u64 {
        self.max_absorption_cells
    }

    pub(super) fn validate_absorption_host(
        &self,
        medium: usize,
        wavelength_nm: f64,
    ) -> Result<(), CurvedSceneError> {
        let Some(SceneMedium::Homogeneous(material)) = self.scene.media.get(medium) else {
            return Err(scattering_error(PolarizedScatteringError::CarrierMedium));
        };
        // The native grid carrier condition, evaluated on its actual host and
        // wavelength. Geometry/phase-index ownership is never replaced by a
        // grid's density or by equality at a single spatial sample.
        if material
            .ior()
            .eval(wavelength_nm)
            .map_err(BeamError::from)?
            != 1.0
            || material
                .absorption_ref()
                .extinction_per_m(wavelength_nm)
                .map_err(BeamError::from)?
                != 0.0
        {
            return Err(scattering_error(PolarizedScatteringError::CarrierMedium));
        }
        Ok(())
    }

    pub(super) fn absorption_segment<const EMIT: bool>(
        &self,
        cx: &Cx<'_>,
        medium: usize,
        wavelength_nm: f64,
        ray: Ray,
        distance_m: f64,
        remaining_cells: u64,
    ) -> Result<Option<GridSegmentTransport<1>>, CurvedSceneError> {
        let Some(grid) = self.absorbing_grid(medium) else {
            return Ok(None);
        };
        self.validate_absorption_host(medium, wavelength_nm)?;
        grid.transport_observed::<1, EMIT>(
            cx,
            ray.origin,
            ray.dir,
            distance_m,
            &[wavelength_nm],
            remaining_cells.min(u64::from(u32::MAX)) as u32,
        )
        .map(Some)
        .map_err(volume_error)
    }

    pub(super) fn transport_absorbing_segment(
        &self,
        cx: &Cx<'_>,
        path: &mut CurvedPathResult,
        ray: Ray,
        distance_m: f64,
    ) -> Result<(), CurvedSceneError> {
        let remaining = self
            .max_absorption_cells
            .checked_sub(path.absorption_cell_visits)
            .ok_or_else(|| volume_error(VolumeTransportError::Budget))?;
        let segment = if self.radiance_brushing {
            self.absorption_segment::<true>(
                cx,
                path.medium,
                path.beam.wavelength_nm(),
                ray,
                distance_m,
                remaining,
            )?
        } else {
            self.absorption_segment::<false>(
                cx,
                path.medium,
                path.beam.wavelength_nm(),
                ray,
                distance_m,
                remaining,
            )?
        };
        let Some(segment) = segment else {
            return Ok(());
        };
        path.absorption_cell_visits += u64::from(segment.cell_visits);
        if self.radiance_brushing {
            super::super::accumulate_reduced_emission(
                path,
                segment.transport.emission_w_m2_sr_nm[0],
            )?;
        }
        let depth = segment.optical_depth[0];
        path.optical_depth = super::super::checked_add(path.optical_depth, depth)?;
        if depth != 0.0 {
            path.beam = attenuate(path.beam, depth)?;
        }
        Ok(())
    }
}

/// Keep native multiplication when representable, recovering a bright prefix
/// whose separately evaluated survival underflows. A nonzero component that
/// cannot be represented after absorption refuses instead of becoming a false
/// extinction or a changed polarization state.
fn attenuate(beam: PolarizedBeam, depth: f64) -> Result<PolarizedBeam, CurvedSceneError> {
    let survival = det::exp(-depth);
    let scale = |component: f64| -> Result<f64, CurvedSceneError> {
        if component == 0.0 {
            return Ok(component);
        }
        let direct = component * survival;
        let value = if direct != 0.0 {
            direct
        } else {
            det::exp(det::ln(component.abs()) - depth).copysign(component)
        };
        if !value.is_finite() || value == 0.0 {
            return Err(volume_error(VolumeTransportError::Unrepresentable));
        }
        Ok(value)
    };
    let incoming = beam.power();
    PolarizedBeam::try_new(
        beam.ray().origin,
        beam.frame(),
        beam.wavelength_nm(),
        beam.medium(),
        Stokes {
            i: scale(incoming.i)?,
            q: scale(incoming.q)?,
            u: scale(incoming.u)?,
            v: scale(incoming.v)?,
        },
    )
    .map_err(CurvedSceneError::from)
}
