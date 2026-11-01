//! Native-grid attenuation over whole emitting patches, before photopic integration.
//!
//! For nonnegative source weights w_lambda, the weighted transmission lies in
//! sum(w_lambda*T_lambda)/sum(w_lambda). Multiplying the original photopic
//! source luminance by this ratio retains its exact represented spectrum and
//! normalization. It is not a gray approximation or centre-ray evaluation.

use super::*;
use crate::photometry::source::{SPECTRAL_BINS, wavelength_nm};
use crate::volumes::transport::enclosure::{
    GRID_BUNDLE_ENCLOSURE_SEMANTICS_VERSION, GRID_EXTINCTION_ENCLOSURE_SEMANTICS_VERSION,
    MAX_PREPARED_EXTINCTION_VALUES, PreparedGridExtinction,
};

struct SourceWeights {
    primitive: usize,
    weights: [Interval; SPECTRAL_BINS],
    sum: Interval,
}

/// One bounded spectral packet per emitter, not per triangle or refinement leaf.
/// Original density, optical palettes and medium bounds remain borrowed.
pub(super) struct PreparedAreaGrid<'field> {
    packet: PreparedGridExtinction<'field, 'field, SPECTRAL_BINS>,
    sources: Vec<SourceWeights>,
    budget: GridEnclosureBudget,
    identity: ContentHash,
}
impl<'field> PreparedAreaGrid<'field> {
    pub(super) fn new(
        medium: PhotometricMedium<'field>, sources: &[PhotometricAreaEmitter],
        budget: GridEnclosureBudget, cx: &Cx<'_>,
    ) -> Result<Self, AreaOcclusionError> {
        cx.checkpoint()?;
        if budget.max_nodes_per_segment == 0 { return Err(AreaOcclusionError::InvalidOptions); }
        let wavelengths: [f64; SPECTRAL_BINS] = std::array::from_fn(wavelength_nm);
        let (packet, identity) = match medium {
            PhotometricMedium::Absorbing(grid) => (
                grid.prepare_extinction_enclosure(cx, &wavelengths, MAX_PREPARED_EXTINCTION_VALUES)?,
                grid.content_identity(),
            ),
            PhotometricMedium::Scattering(grid) => (
                grid.prepare_extinction_enclosure(cx, &wavelengths, MAX_PREPARED_EXTINCTION_VALUES)?,
                grid.content_identity(),
            ),
            PhotometricMedium::Homogeneous(_) => return Err(AreaOcclusionError::InvalidOptions),
        };
        let mut prepared = Vec::new();
        // Source count is already admitted by AreaReferenceBudget. The extra
        // storage is 80 intervals plus one sum and one index per source, never
        // a Cartesian product with the grid's independent optical palettes.
        prepared.try_reserve_exact(sources.len()).map_err(|_| AreaOcclusionError::Allocation)?;
        for source in sources {
            cx.checkpoint()?;
            let weights: [Interval; SPECTRAL_BINS] = std::array::from_fn(|bin| mul(
                point(source.spectrum().bins_w_per_lumen_nm()[bin]),
                point(observer::xyz_at_bin(bin)[1]),
            ));
            let sum = weights.iter().copied().fold(zero(), add);
            if !finite(sum) || sum.lo() <= 0.0 { return Err(AreaReferenceError::NumericalRange.into()); }
            prepared.push(SourceWeights { primitive: source.primitive_index(), weights, sum });
        }
        prepared.sort_unstable_by_key(|source| source.primitive);
        cx.checkpoint()?;
        Ok(Self { packet, sources: prepared, budget, identity })
    }

    pub(super) const fn identity(&self) -> ContentHash { self.identity }

    pub(super) fn hash_into(&self, hash: &mut DomainHasher) {
        hash.update(b"physical-spectral-grid");
        hash.update(&self.identity.0);
        hash.update(&GRID_EXTINCTION_ENCLOSURE_SEMANTICS_VERSION.to_le_bytes());
        hash.update(&GRID_BUNDLE_ENCLOSURE_SEMANTICS_VERSION.to_le_bytes());
        for bin in 0..SPECTRAL_BINS { hash_float(hash, wavelength_nm(bin)); }
        // Spectrum and observer weights already belong to the original source
        // preparation identity. Work budgets do not change the enclosed model.
    }

    pub(super) fn photopic_transfer(
        &self, cx: &Cx<'_>, primitive: usize, origin: Vector, receiver: Vector,
    ) -> Result<Interval, AreaOcclusionError> {
        cx.checkpoint()?;
        let index = self.sources.binary_search_by_key(&primitive, |source| source.primitive)
            .map_err(|_| AreaReferenceError::InvalidGeometry)?;
        let source = &self.sources[index];
        let transmission = self.packet.transmittance_box(cx, origin, receiver, self.budget)?;
        // A common packet range bounds any positive weighted average exactly;
        // in particular an empty/vacuum field retains the original one bits.
        if transmission.iter().all(|t| *t == transmission[0]) { return Ok(transmission[0]); }
        let mut sum = zero();
        for (bin, (&weight, value)) in source.weights.iter().zip(transmission).enumerate() {
            if bin % 16 == 0 { cx.checkpoint()?; }
            sum = add(sum, mul(weight, value));
        }
        let fraction = sum / source.sum;
        let lo = fraction.lo().max(0.0);
        let hi = fraction.hi().min(1.0);
        if !lo.is_finite() || !hi.is_finite() || lo > hi { return Err(AreaReferenceError::NumericalRange.into()); }
        cx.checkpoint()?;
        Ok(Interval::new(lo, hi))
    }
}
