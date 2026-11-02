//! Spectral source terms through the ORIGINAL native grid snapshot.
//! One shared extinction preparation, one observer-weight packet per fixture,
//! and fixed stack packets per query. No density clone, per-cell spectrum copy,
//! per-panel allocation or representative gray optical coefficient is used.
use super::*;
use crate::photometry::source::{SPECTRAL_BINS, LUMINOUS_EFFICACY_LM_W, photopic_weight_nm};
use crate::photometry::transport::PhotometricMedium;
use crate::volumes::transport::enclosure::{
    GridEnclosureBudget, MAX_PREPARED_EXTINCTION_VALUES, PreparedGridExtinction, PreparedGridScattering,
};

pub(super) struct GridScattering<'field> {
    pub(super) packet: PreparedGridScattering<'field, 'field, SPECTRAL_BINS>,
    pub(super) weights: Vec<[Interval; SPECTRAL_BINS]>,
    budget: GridEnclosureBudget,
    identity: ContentHash,
}
impl<'field> GridScattering<'field> {
    pub(super) fn new(scene: &PhotometricScene<'field>,
        packet: PreparedGridExtinction<'field, 'field, SPECTRAL_BINS>, budget: GridEnclosureBudget, cx: &Cx<'_>)
        -> Result<Self, ScatterEnclosureError> {
        cx.checkpoint()?;
        let identity = match &scene.medium {
            PhotometricMedium::Absorbing(grid) => grid.content_identity(),
            PhotometricMedium::Scattering(grid) => grid.content_identity(),
            PhotometricMedium::Homogeneous(_) => return Err(ScatterEnclosureError::HomogeneousOnly),
        };
        // Independently bounded from the native optical coefficient table.
        // At 80 intervals per source this cap admits at most 4096 fixtures.
        if scene.sources.len() > MAX_PREPARED_EXTINCTION_VALUES / SPECTRAL_BINS {
            return Err(ScatterEnclosureError::Budget);
        }
        let packet = packet.into_single_scattering(cx, MAX_PREPARED_EXTINCTION_VALUES)?;
        let mut weights = Vec::new();
        weights.try_reserve_exact(scene.sources.len()).map_err(|_| ScatterEnclosureError::Allocation)?;
        for source in scene.sources {
            cx.checkpoint()?;
            let mut values = [zero(); SPECTRAL_BINS];
            for (bin, (&power, value)) in source.spectrum().bins_w_per_lumen_nm().iter().zip(&mut values).enumerate() {
                *value = mul(mul(point(power), point(photopic_weight_nm(bin))), point(LUMINOUS_EFFICACY_LM_W));
                if !finite(*value) { return Err(ScatterEnclosureError::NumericalRange); }
            }
            weights.push(values);
        }
        cx.checkpoint()?;
        Ok(Self { packet, weights, budget, identity })
    }
    pub(super) fn identity(&self) -> ContentHash { self.identity }
    pub(super) fn active(&self) -> bool { self.packet.has_scattering() }

    // Integrate rho*kappa_s*phase exactly as a retained cell model, enclosing
    // phase variation over the whole camera panel. Attenuation varies within
    // the panel but is bounded for every endpoint by both ray-family queries.
    // Products of nonnegative enclosures remain valid despite their shared
    // spatial dependencies; refinement narrows those dependencies.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn response(&self, cx: &Cx<'_>, source: usize, light: Vector, camera: Vector,
        start: Vector, end: Vector, position: Vector, cosine: Interval)
        -> Result<Interval, ScatterEnclosureError> {
        let scattering = self.packet.scattering_integral_box(cx, start, end, cosine, self.budget)?;
        let weights = &self.weights[source];
        if scattering.iter().zip(weights).all(|(s,w)| s.hi() == 0.0 || w.hi() == 0.0) { return Ok(zero()); }
        let light_t = self.packet.extinction().transmittance_box(cx, light, position, self.budget)?;
        let camera_t = self.packet.extinction().transmittance_box(cx, camera, position, self.budget)?;
        let mut value = zero();
        for bin in 0..SPECTRAL_BINS {
            if bin % 16 == 0 { cx.checkpoint()?; }
            value = add(value, mul(mul(mul(weights[bin], scattering[bin]), light_t[bin]), camera_t[bin]));
            if !finite(value) { return Err(ScatterEnclosureError::NumericalRange); }
        }
        cx.checkpoint()?;
        Ok(value)
    }
}
