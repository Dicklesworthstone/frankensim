//! Prepared, source-bound optical properties in every original physical cell.
//!
//! Mass absorption/scattering are borrowed m²/kg packets, evaluated by the
//! constitutive owner at the accepted local state. They replace the geometry
//! owner's uniform coefficients. Density remains its original kg/m³ snapshot.
//! No wavelength or temperature interpolation is performed here: the fixed
//! packet and explicit source bounds are immutable for one transport operation.

use super::super::{AbsorbingGrid, SegmentTransport, cells::validate_ray, checkpoint};
use super::deposition::{ScatteringDepositionConfig, ScatteringPowerDeposition};
use super::prepared::{self, PreparedOptics};
use super::{FreeFlightSample, PhaseFunction, VolumeTransportError, unit_variate};
use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::Cx;
use fs_geom::{Point3, Vec3};

/// Borrow one reciprocal phase law, or an indexed palette shared by cells.
/// Large spectral/angular tables are never copied per cell or collision.
#[derive(Clone, Copy)]
pub enum CellPhaseFunction<'a> {
    Uniform(&'a PhaseFunction),
    Indexed {
        palette: &'a [PhaseFunction],
        /// One palette index per original X-major physical cell.
        cell_indices: &'a [usize],
    },
}

impl CellPhaseFunction<'_> {
    fn at(&self, cell: usize) -> &PhaseFunction {
        match self {
            Self::Uniform(phase) => phase,
            Self::Indexed {
                palette,
                cell_indices,
            } => &palette[cell_indices[cell]],
        }
    }

    fn admit(
        &self,
        check: &mut (impl FnMut() -> Result<(), VolumeTransportError> + ?Sized),
        cells: usize,
        bounds: (f64, f64),
        hash: &mut DomainHasher,
    ) -> Result<(), VolumeTransportError> {
        let validate = |phase: &PhaseFunction| {
            phase.validate()?;
            phase.validate_wavelength(bounds.0)?;
            phase.validate_wavelength(bounds.1)
        };
        match self {
            Self::Uniform(phase) => {
                check()?;
                validate(phase)?;
                hash.update(b"uniform-phase");
                phase.hash_into(hash);
            }
            Self::Indexed {
                palette,
                cell_indices,
            } => {
                if cell_indices.len() != cells || palette.is_empty() {
                    return Err(VolumeTransportError::Phase);
                }
                // Admission work stays bounded by the geometry's accepted cell
                // count; unused, arbitrarily large palettes confer no support.
                if palette.len() > cells {
                    return Err(VolumeTransportError::Budget);
                }
                hash.update(b"indexed-phase");
                hash.update(&(palette.len() as u64).to_le_bytes());
                for phase in *palette {
                    check()?;
                    validate(phase)?;
                    phase.hash_into(hash);
                }
                for (cell, &index) in cell_indices.iter().enumerate() {
                    if cell % 1024 == 0 {
                        check()?;
                    }
                    if index >= palette.len() {
                        return Err(VolumeTransportError::Phase);
                    }
                    hash.update(&(index as u64).to_le_bytes());
                }
            }
        }
        Ok(())
    }
}

/// Immutable packet field borrowing the accepted geometry and coefficient
/// buffers. Storage beyond those caller-owned arrays is O(packet width).
///
/// Every original cell has separate absorption/scattering, including cells
/// outside the incident beam: a scattered path can enter any cell. Source
/// identity binds the constitutive data and conditions, not measured-material
/// authority. The declared support is not inferred from packet extrema.
#[derive(Clone, Copy)]
pub struct CellScatteringGrid<'g, 'a, const N: usize> {
    geometry: &'g AbsorbingGrid<'a>,
    wavelengths_nm: [f64; N],
    absorption_m2_kg: &'g [[f64; N]],
    scattering_m2_kg: &'g [[f64; N]],
    phase: CellPhaseFunction<'g>,
    bounds_nm: (f64, f64),
    source_identity: ContentHash,
    identity: ContentHash,
    scattering_lanes: [bool; N],
    extinction_lanes: [bool; N],
    has_emission: bool,
}

impl<'g, 'a, const N: usize> CellScatteringGrid<'g, 'a, N> {
    /// Validate borrowed SI coefficients once, with no allocation or spectra
    /// copied into cells. `geometry` supplies density, metre geometry, optional
    /// kelvin temperatures and DDA work limits; its uniform spectrum is replaced.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        cx: &Cx<'_>,
        geometry: &'g AbsorbingGrid<'a>,
        wavelengths_nm: [f64; N],
        absorption_m2_kg: &'g [[f64; N]],
        scattering_m2_kg: &'g [[f64; N]],
        phase: CellPhaseFunction<'g>,
        wavelength_bounds_nm: (f64, f64),
        source_identity: ContentHash,
    ) -> Result<Self, VolumeTransportError> {
        Self::try_new_with_checkpoint(
            &mut || checkpoint(cx),
            geometry,
            wavelengths_nm,
            absorption_m2_kg,
            scattering_m2_kg,
            phase,
            wavelength_bounds_nm,
            source_identity,
        )
    }

    /// Same admission with the hosting solver's cancellation/budget checkpoint.
    /// The callback's exact refusal is propagated; no partial field is returned.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new_with_checkpoint(
        check: &mut (impl FnMut() -> Result<(), VolumeTransportError> + ?Sized),
        geometry: &'g AbsorbingGrid<'a>,
        wavelengths_nm: [f64; N],
        absorption_m2_kg: &'g [[f64; N]],
        scattering_m2_kg: &'g [[f64; N]],
        phase: CellPhaseFunction<'g>,
        wavelength_bounds_nm: (f64, f64),
        source_identity: ContentHash,
    ) -> Result<Self, VolumeTransportError> {
        check()?;
        if N == 0 {
            return Err(VolumeTransportError::Sample);
        }
        if source_identity == ContentHash([0; 32]) {
            return Err(VolumeTransportError::SourceIdentity);
        }
        let (lo, hi) = wavelength_bounds_nm;
        if !lo.is_finite()
            || !hi.is_finite()
            || lo <= 0.0
            || lo >= hi
            || wavelengths_nm
                .iter()
                .any(|l| !l.is_finite() || !(lo..=hi).contains(l))
        {
            return Err(VolumeTransportError::Wavelength);
        }
        let cells = geometry.cell_count();
        if absorption_m2_kg.len() != cells || scattering_m2_kg.len() != cells {
            return Err(VolumeTransportError::Grid);
        }
        let mut hash = DomainHasher::new("org.frankensim.render.cell-optical-packet.v1");
        hash.update(geometry.content_identity().as_bytes());
        hash.update(source_identity.as_bytes());
        hash.update(&(N as u64).to_le_bytes());
        for value in [lo, hi].iter().chain(wavelengths_nm.iter()) {
            hash.update(&value.to_bits().to_le_bytes());
        }
        phase.admit(check, cells, wavelength_bounds_nm, &mut hash)?;
        let mut scattering_lanes = [false; N];
        let mut extinction_lanes = [false; N];
        let mut has_absorption = false;
        for cell in 0..cells {
            check()?;
            let dense = geometry.grid.data[cell] > 0.0;
            for lane in 0..N {
                let a = absorption_m2_kg[cell][lane];
                let s = scattering_m2_kg[cell][lane];
                if !a.is_finite() || a < 0.0 {
                    return Err(VolumeTransportError::AbsorptionSpectrum);
                }
                if !s.is_finite() || s < 0.0 {
                    return Err(VolumeTransportError::ScatteringSpectrum);
                }
                if !(a + s).is_finite() {
                    return Err(VolumeTransportError::Unrepresentable);
                }
                scattering_lanes[lane] |= dense && s > 0.0;
                extinction_lanes[lane] |= dense && (a > 0.0 || s > 0.0);
                has_absorption |= dense && a > 0.0;
                hash.update(&a.to_bits().to_le_bytes());
                hash.update(&s.to_bits().to_le_bytes());
            }
        }
        check()?;
        Ok(Self {
            geometry,
            wavelengths_nm,
            absorption_m2_kg,
            scattering_m2_kg,
            phase,
            bounds_nm: wavelength_bounds_nm,
            source_identity,
            identity: hash.finalize(),
            scattering_lanes,
            extinction_lanes,
            has_emission: has_absorption && geometry.temperature_k.is_some(),
        })
    }

    #[must_use]
    pub const fn content_identity(&self) -> ContentHash {
        self.identity
    }
    #[must_use]
    pub const fn source_identity(&self) -> ContentHash {
        self.source_identity
    }
    #[must_use]
    pub const fn wavelength_bounds_nm(&self) -> (f64, f64) {
        self.bounds_nm
    }
    #[must_use]
    pub const fn wavelengths_nm(&self) -> &[f64; N] {
        &self.wavelengths_nm
    }
    #[must_use]
    pub const fn absorption_grid(&self) -> &AbsorbingGrid<'a> {
        self.geometry
    }
    #[must_use]
    pub fn has_scattering(&self) -> bool {
        self.scattering_lanes.iter().any(|s| *s)
    }
    #[must_use]
    pub fn phase_ref_at(&self, cell: usize) -> Option<&PhaseFunction> {
        (cell < self.geometry.cell_count()).then(|| self.phase.at(cell))
    }

    /// Exact shadow attenuation through local extinction at the fixed packet.
    pub fn transmittance(
        &self,
        cx: &Cx<'_>,
        origin: Point3,
        direction: Vec3,
        max_distance_m: f64,
    ) -> Result<[f64; N], VolumeTransportError> {
        validate_ray(origin, direction, max_distance_m)?;
        Ok(prepared::transport(
            self,
            &mut || checkpoint(cx),
            origin,
            direction,
            max_distance_m,
            false,
        )?
        .transmittance)
    }

    /// Exact attenuation and local-equilibrium source. The local absorption
    /// creates emission, while the local total extinction attenuates it.
    pub fn transport(
        &self,
        cx: &Cx<'_>,
        origin: Point3,
        direction: Vec3,
        max_distance_m: f64,
    ) -> Result<SegmentTransport<N>, VolumeTransportError> {
        validate_ray(origin, direction, max_distance_m)?;
        prepared::transport(
            self,
            &mut || checkpoint(cx),
            origin,
            direction,
            max_distance_m,
            true,
        )
    }

    /// Native spectral-mixture free flight using each crossed cell's coefficients.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_segment(
        &self,
        cx: &Cx<'_>,
        origin: Point3,
        direction: Vec3,
        max_distance_m: f64,
        u_lane: f64,
        u_distance: f64,
    ) -> Result<FreeFlightSample<N>, VolumeTransportError> {
        validate_ray(origin, direction, max_distance_m)?;
        if !unit_variate(u_lane) || !unit_variate(u_distance) {
            return Err(VolumeTransportError::Sample);
        }
        prepared::sample(
            self,
            &mut || checkpoint(cx),
            origin,
            direction,
            max_distance_m,
            u_lane,
            u_distance,
            true,
        )
    }

    /// Exact spatial removal from the directed beam; only absorption is heat.
    pub fn deposit_beam_power(
        &self,
        cx: &Cx<'_>,
        origin: Point3,
        unit_direction: Vec3,
        max_distance_m: f64,
        incident_power_w: &[f64; N],
        max_deposition_cells: usize,
    ) -> Result<super::super::deposition::BeamPowerDeposition<N>, VolumeTransportError> {
        super::super::deposition::deposit_prepared(
            &mut || checkpoint(cx),
            self.geometry,
            (origin, unit_direction, max_distance_m),
            &self.wavelengths_nm,
            incident_power_w,
            max_deposition_cells,
            (self.identity, self.bounds_nm),
            |cell| self.coefficients(cell),
        )
    }

    /// Reabsorption after scattering, with local collision albedo and phase.
    /// Absorption-only lanes retain exact spatial Beer transport and zero sample
    /// error. Other lanes use the uniform owner's same histories and reductions.
    pub fn deposit_multiple_scattering_power(
        &self,
        cx: &Cx<'_>,
        origin: Point3,
        unit_direction: Vec3,
        incident_power_w: &[f64; N],
        config: ScatteringDepositionConfig,
    ) -> Result<ScatteringPowerDeposition<N>, VolumeTransportError> {
        self.deposit_multiple_scattering_power_with_checkpoint(
            &mut || checkpoint(cx),
            origin,
            unit_direction,
            incident_power_w,
            config,
        )
    }

    /// Hosting-solver checkpoint variant. Per-cell crossings, histories,
    /// allocation budgets and reduction checkpoints propagate the exact error.
    pub fn deposit_multiple_scattering_power_with_checkpoint(
        &self,
        check: &mut (impl FnMut() -> Result<(), VolumeTransportError> + ?Sized),
        origin: Point3,
        unit_direction: Vec3,
        incident_power_w: &[f64; N],
        config: ScatteringDepositionConfig,
    ) -> Result<ScatteringPowerDeposition<N>, VolumeTransportError> {
        super::deposition::deposit_prepared(
            self,
            check,
            origin,
            unit_direction,
            incident_power_w,
            config,
        )
    }
}

impl<const N: usize> PreparedOptics<N> for CellScatteringGrid<'_, '_, N> {
    fn geometry(&self) -> &AbsorbingGrid<'_> {
        self.geometry
    }
    fn identity(&self) -> ContentHash {
        self.identity
    }
    fn bounds(&self) -> (f64, f64) {
        self.bounds_nm
    }
    fn wavelengths(&self) -> &[f64; N] {
        &self.wavelengths_nm
    }
    fn coefficient(
        &self,
        cell: usize,
        lane: usize,
    ) -> Result<(f64, f64, f64), VolumeTransportError> {
        let a = self.absorption_m2_kg[cell][lane];
        let s = self.scattering_m2_kg[cell][lane];
        Ok((a, s, a + s))
    }
    fn coefficients(
        &self,
        cell: usize,
    ) -> Result<([f64; N], [f64; N], [f64; N]), VolumeTransportError> {
        let a = self.absorption_m2_kg[cell];
        let s = self.scattering_m2_kg[cell];
        Ok((a, s, std::array::from_fn(|lane| a[lane] + s[lane])))
    }
    fn phase(&self, cell: usize) -> &PhaseFunction {
        self.phase.at(cell)
    }
    fn has_scattering(&self) -> bool {
        self.has_scattering()
    }
    fn lane_has_scattering(&self, lane: usize) -> bool {
        self.scattering_lanes[lane]
    }
    fn lane_has_extinction(&self, lane: usize) -> bool {
        self.extinction_lanes[lane]
    }
    fn has_emission(&self) -> bool {
        self.has_emission
    }
}

#[cfg(test)]
mod tests;
