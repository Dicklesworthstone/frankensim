//! Exact incoherent slab absorption deposited at its physical normal depths.
//!
//! The original glass remains one homogeneous optical layer. Thermal cells
//! partition its positive forward/backward Beer losses without adding optical
//! interfaces. All internal reflections are retained by the prepared slab
//! owner. This adapter reuses the existing surface observer, source schedules,
//! cell validation, energy accounts and transactional native enthalpy advance.
//!
//! Optical laws and geometry are frozen over their declared receiver support.
//! Absorbed radiation thermalizes promptly; carrier storage and fluorescence
//! require their own physical owners. Emission, cooling, conduction and phase
//! change are not inferred or duplicated by this incident-heating adapter.

use super::coating::PreparedCoatingSource;
use super::response::{
    RESPONSE_ROUNDOFF, hash_geometry, hash_values, mix, response_hash, validate_packet,
};
use super::{
    CheckedSum, DepositionError, DepositionSimulation, SpectralPowerSchedule, SurfaceHeatingBudget,
    SurfaceHeatingError, SurfaceHeatingReceiver, SurfaceHeatingTolerance, SurfaceOpticalHeating,
    SurfacePowerResponse, invalid,
};
use crate::dielectric::{DielectricGlass, DielectricSurface};
use crate::tracer::Material;
use fs_blake3::ContentHash;

/// Frozen spectral slab powers with a complete physical depth partition.
#[derive(Clone, Debug, PartialEq)]
pub struct SlabDepthPowerResponse<const N: usize> {
    total: SurfacePowerResponse<N>,
    depth_edges_m: Box<[f64]>,
    absorbed: Box<[[f64; N]]>,
}

impl<const N: usize> SlabDepthPowerResponse<N> {
    /// Resolve the original glass once per wavelength and integrate each depth
    /// interval analytically. `depth_edges_m` starts at zero and increases
    /// strictly to the slab's positive normal thickness. No interval can be
    /// omitted, overlap, or create an artificial optical interface.
    ///
    /// Receiver order always follows authored depth order. `reverse` illuminates
    /// the last depth first; incident/exit indices and cosine must describe that
    /// actual incidence. Powers are intercepted exterior normal watts, not
    /// radiance and not power already admitted through the first interface.
    /// The caller supplies original source identity and a depth-storage ceiling.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        glass: &DielectricGlass,
        wavelengths_nm: [f64; N],
        depth_edges_m: &[f64],
        incident_eta: [f64; N],
        exit_eta: [f64; N],
        incidence_cosine: f64,
        reverse: bool,
        s_power_fraction: [f64; N],
        source_identity: ContentHash,
        max_depth_cells: usize,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, SurfaceHeatingError> {
        checkpoint()?;
        validate_packet(&wavelengths_nm, &s_power_fraction)?;
        if depth_edges_m.len() < 2 {
            return Err(invalid(
                "slab heating needs a complete positive depth partition",
            ));
        }
        let count = depth_edges_m.len() - 1;
        if count > max_depth_cells {
            return Err(SurfaceHeatingError::Budget);
        }
        if depth_edges_m[0] != 0.0 {
            return Err(invalid(
                "slab depths must increase from zero to finite thickness",
            ));
        }
        for pair in depth_edges_m.windows(2) {
            checkpoint()?;
            if !pair[1].is_finite() || pair[1] <= pair[0] {
                return Err(invalid(
                    "slab depths must increase from zero to finite thickness",
                ));
            }
        }
        let thickness = depth_edges_m[count];
        let mut hash = response_hash(b"incoherent-slab-depth-heat", source_identity)?;
        hash.update(
            Material::Dielectric {
                glass: *glass,
                surface: DielectricSurface::SMOOTH,
            }
            .content_identity()
            .as_bytes(),
        );
        hash_geometry(
            &mut hash,
            incidence_cosine,
            &incident_eta,
            &s_power_fraction,
        );
        hash_values(&mut hash, exit_eta);
        hash.update(&[u8::from(reverse)]);
        hash.update(&(count as u64).to_le_bytes());
        let mut depths = Vec::new();
        depths
            .try_reserve_exact(count + 1)
            .map_err(|_| SurfaceHeatingError::Budget)?;
        let mut absorbed = Vec::new();
        absorbed
            .try_reserve_exact(count)
            .map_err(|_| SurfaceHeatingError::Budget)?;
        for &depth in depth_edges_m {
            checkpoint()?;
            depths.push(depth);
            hash.update(&depth.to_bits().to_le_bytes());
        }
        absorbed.resize(count, [0.0; N]);
        let mut reflected = [0.0; N];
        let mut transmitted = [0.0; N];
        let mut total_absorbed = [0.0; N];
        for lane in 0..N {
            checkpoint()?;
            let profile = glass.prepare_incoherent_slab_absorption_between(
                wavelengths_nm[lane],
                thickness,
                incidence_cosine,
                incident_eta[lane],
                exit_eta[lane],
            )?;
            let total = profile.optical();
            let modes = [total.perpendicular, total.parallel];
            reflected[lane] = mix(modes.map(|p| p.reflectance), s_power_fraction[lane]);
            transmitted[lane] = mix(modes.map(|p| p.transmittance), s_power_fraction[lane]);
            total_absorbed[lane] = mix(modes.map(|p| p.absorptance), s_power_fraction[lane]);
            let mut sum = CheckedSum::default();
            for (cell, edges) in absorbed.iter_mut().zip(depth_edges_m.windows(2)) {
                checkpoint()?;
                let fraction = if reverse {
                    profile.absorption_between_from_back(edges[0], edges[1])?
                } else {
                    profile.absorption_between(edges[0], edges[1])?
                };
                cell[lane] = mix(fraction, s_power_fraction[lane]);
                if !cell[lane].is_finite() || !(0.0..=1.0).contains(&cell[lane]) {
                    return Err(invalid("slab depth absorption must be finite and passive"));
                }
                sum.add(cell[lane])?;
                hash.update(&cell[lane].to_bits().to_le_bytes());
            }
            if (sum.value()? - total_absorbed[lane]).abs() > RESPONSE_ROUNDOFF {
                return Err(invalid(
                    "slab depth heat must close original bulk absorption",
                ));
            }
        }
        let total = SurfacePowerResponse::finish(
            hash,
            wavelengths_nm,
            reflected,
            transmitted,
            total_absorbed,
        )?;
        checkpoint()?;
        Ok(Self {
            total,
            depth_edges_m: depths.into_boxed_slice(),
            absorbed: absorbed.into_boxed_slice(),
        })
    }

    /// Original reflection, transmission and total absorption, counted once.
    pub fn total(&self) -> &SurfacePowerResponse<N> {
        &self.total
    }

    /// Authored normal-depth boundaries in metres, independent of illumination side.
    pub fn depth_edges_m(&self) -> &[f64] {
        &self.depth_edges_m
    }

    /// Absorbed fractions of exterior incident watts, in authored depth order.
    pub fn depth_absorbed_fractions(&self) -> &[[f64; N]] {
        &self.absorbed
    }

    /// Original optical laws, source, geometry and complete depth partition.
    pub fn identity(&self) -> ContentHash {
        self.total.identity()
    }
}

/// One physical slab illumination and its existing depth-cell receivers.
#[derive(Clone, Debug)]
pub struct SlabDepthIlluminationSource<const N: usize> {
    pub source_id: u64,
    /// Exactly one receiver per response interval, in authored depth order.
    /// Several intervals may share a receiver; their heat adds only once.
    /// This mapping and the cells' physical volumes are caller-authored:
    /// a generic thermal receiver does not establish slab shape or coverage.
    pub depth_receivers: Vec<SurfaceHeatingReceiver>,
    pub response: SlabDepthPowerResponse<N>,
    pub power: SpectralPowerSchedule<N>,
}

impl<const N: usize> SurfaceOpticalHeating<N> {
    /// Bind complete depth partitions to the original native thermal observer.
    /// All source and receiver accounts use the existing surface time integrator.
    /// `max_total_depth_cells` bounds retained bindings across the entire source
    /// set before allocation. Refusal or cancellation publishes no partial state.
    pub fn try_new_slab_depths(
        sources: &[SlabDepthIlluminationSource<N>],
        thermal: &DepositionSimulation<'_>,
        budget: SurfaceHeatingBudget,
        max_total_depth_cells: usize,
        tolerance: SurfaceHeatingTolerance,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, SurfaceHeatingError> {
        checkpoint()?;
        if sources.len() > budget.max_sources {
            return Err(SurfaceHeatingError::Budget);
        }
        let mut count = 0_usize;
        for source in sources {
            checkpoint()?;
            if source.depth_receivers.len() != source.response.absorbed.len() {
                return Err(invalid(
                    "one thermal receiver per slab depth interval required",
                ));
            }
            count = count
                .checked_add(source.depth_receivers.len())
                .ok_or(SurfaceHeatingError::Budget)?;
            if count > max_total_depth_cells {
                return Err(SurfaceHeatingError::Budget);
            }
        }
        let mut prepared = Vec::new();
        prepared
            .try_reserve_exact(sources.len())
            .map_err(|_| SurfaceHeatingError::Budget)?;
        for source in sources {
            checkpoint()?;
            let mut bindings = Vec::new();
            bindings
                .try_reserve_exact(source.depth_receivers.len())
                .map_err(|_| SurfaceHeatingError::Budget)?;
            for (&receiver, &fraction) in source
                .depth_receivers
                .iter()
                .zip(source.response.absorbed.iter())
            {
                checkpoint()?;
                bindings.push((receiver, fraction));
            }
            prepared.push(PreparedCoatingSource {
                source_id: source.source_id,
                receiver: None,
                bindings,
                response: source.response.total,
                power: source.power.clone(),
            });
        }
        Self::try_new_partitioned(prepared, thermal, budget, tolerance, checkpoint)
    }
}
