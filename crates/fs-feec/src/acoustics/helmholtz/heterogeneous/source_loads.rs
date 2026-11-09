//! Cached original source columns for real full-pressure multiphysics evolution.

use super::{
    BoundaryVelocitySource, HeterogeneousError, HeterogeneousHelmholtzFe, VolumeVelocitySource,
    boundary_source, invalid, poll, source,
};
use fs_math::dd::Dd;

/// Additional source-column storage and work, independent of evolution limits.
#[derive(Clone, Copy, Debug)]
pub struct PressureSourceLoadLimits {
    /// Cached real columns plus one full pressure-load workspace [bytes].
    /// Original sources, evolution caches and allocator metadata are separate.
    pub maximum_payload_bytes: usize,
    /// Original pressure coordinates times the complete source count.
    pub maximum_products_per_projection: usize,
}

impl Default for PressureSourceLoadLimits {
    fn default() -> Self {
        Self {
            maximum_payload_bytes: 16 * 1024 * 1024,
            maximum_products_per_projection: 1_000_000,
        }
    }
}

/// Complete original real source columns, bound to their fluid material/geometry.
/// This is a load map, not a reduced acoustic state or a transient material law.
pub struct PreparedPressureSourceLoads<'m, 'f, 'c> {
    model: &'m HeterogeneousHelmholtzFe<'f, 'c>,
    boundary_count: usize,
    volume_count: usize,
    /// Row-major original pressure coordinate x source; no physical rescaling.
    loads: Vec<f64>,
}

impl<'f, 'c> HeterogeneousHelmholtzFe<'f, 'c> {
    /// Retain already integrated boundary and volume columns for real drives.
    ///
    /// This map is for complete pressure-potential pencils without eliminated
    /// pressure coordinates. Boundary sources retain their original P2 support
    /// and volume sources their complete spatial shape, including signed and
    /// zero-integral sources. Exact source-owner and count admission precedes
    /// storage allocation, even when later amplitudes will be zero. No callback,
    /// material sample or geometry quadrature is repeated. Linear velocity
    /// shapes can be supplied through the original P2 source by using their
    /// exact corner and midpoint values.
    pub fn prepare_source_loads<'m>(
        &'m self,
        boundary_sources: &[&BoundaryVelocitySource<'_, '_, '_>],
        volume_sources: &[&VolumeVelocitySource<'_, '_, '_>],
        limits: PressureSourceLoadLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedPressureSourceLoads<'m, 'f, 'c>, HeterogeneousError> {
        poll(cancelled)?;
        boundary_source::admit_count(boundary_sources.len())?;
        source::admit_count(volume_sources.len())?;
        let count = boundary_sources.len() + volume_sources.len();
        let entries = self
            .ndof()
            .checked_mul(count)
            .ok_or_else(|| invalid("source-load extent overflow"))?;
        let payload = entries
            .checked_add(self.ndof())
            .and_then(|n| n.checked_mul(8))
            .ok_or_else(|| invalid("source-load payload overflow"))?;
        for (what, requested, maximum) in [
            (
                "source-load products",
                entries,
                limits.maximum_products_per_projection,
            ),
            (
                "source-load payload bytes",
                payload,
                limits.maximum_payload_bytes,
            ),
        ] {
            if requested > maximum {
                return Err(HeterogeneousError::Limit {
                    what,
                    requested,
                    maximum,
                });
            }
        }
        for source in boundary_sources {
            poll(cancelled)?;
            source.validate(self)?;
        }
        for source in volume_sources {
            poll(cancelled)?;
            source.validate(self)?;
        }
        let mut loads = Vec::new();
        loads
            .try_reserve_exact(entries)
            .map_err(|_| invalid("source-load allocation refused"))?;
        for row in 0..self.ndof() {
            poll(cancelled)?;
            loads.extend(boundary_sources.iter().map(|source| source.load[row]));
            loads.extend(volume_sources.iter().map(|source| source.load[row]));
        }
        poll(cancelled)?;
        Ok(PreparedPressureSourceLoads {
            model: self,
            boundary_count: boundary_sources.len(),
            volume_count: volume_sources.len(),
            loads,
        })
    }
}

impl PreparedPressureSourceLoads<'_, '_, '_> {
    pub fn pressure_dofs(&self) -> usize {
        self.model.ndof()
    }
    pub fn boundary_sources(&self) -> usize {
        self.boundary_count
    }
    pub fn volume_sources(&self) -> usize {
        self.volume_count
    }

    /// Real inward boundary velocities [m/s] and volume-flow amplitudes [m³/s]
    /// give the original generalized volume flow [m³/s]. The pressure-potential
    /// equation receives this load directly, without a phasor factor, frequency
    /// or density. Its actual trajectory work is the pressure integral dotted
    /// with the held load; there is no harmonic factor of one half.
    ///
    /// No allocation, material query or geometric integration occurs. Owner,
    /// shape and finite-amplitude admission precede writing. Later cancellation
    /// or overflow may alter this caller-owned scratch; it is never accepted
    /// physical state and must not be published after a refusal.
    pub fn project_into(
        &self,
        model: &HeterogeneousHelmholtzFe<'_, '_>,
        boundary_amplitudes_m_s: &[f64],
        volume_amplitudes_m3_s: &[f64],
        load: &mut [f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), HeterogeneousError> {
        poll(cancelled)?;
        if !core::ptr::eq(self.model, model)
            || boundary_amplitudes_m_s.len() != self.boundary_count
            || volume_amplitudes_m3_s.len() != self.volume_count
            || load.len() != self.pressure_dofs()
            || boundary_amplitudes_m_s
                .iter()
                .chain(volume_amplitudes_m3_s)
                .any(|value| !value.is_finite())
        {
            return Err(invalid("source-load owner, dimensions or finite real amplitudes").into());
        }
        let count = self.boundary_count + self.volume_count;
        for (row, value) in load.iter_mut().enumerate() {
            poll(cancelled)?;
            let mut sum = Dd::from_f64(0.0);
            for (&weight, &amplitude) in self.loads[row * count..(row + 1) * count]
                .iter()
                .zip(boundary_amplitudes_m_s.iter().chain(volume_amplitudes_m3_s))
            {
                sum = sum + Dd::from_f64(weight) * Dd::from_f64(amplitude);
            }
            *value = sum.to_f64();
            if !value.is_finite() {
                return Err(invalid("real source-load accumulation overflow").into());
            }
        }
        poll(cancelled)
    }
}
