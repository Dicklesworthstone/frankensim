//! A positive-operator envelope for volume-only orders in a finite native grid.
//!
//! Let D enclose the box diagonal, s bound rho*kappa_s, P bound the phase
//! density, and a be an available original albedo/phase contraction bound.
//! Along ANY straight chord, integral T*s ds <= min(1, D*s). Consequently
//! q <= min(a, 4*pi*P*min(1,D*s)) bounds one further volume-scattering event.
//! The first-order response to bounded incident point-source irradiance J is
//! at most B1 = J*P*min(1,D*s,a). Higher orders are <= B1*q/(1-q), when q<1.
//! No surface reflection, thermal source, spatial/spectral approximation or
//! physical calibration is certified. The caller must establish finite J.

use super::*;

/// Version of the finite-grid chord, phase and original contraction envelope.
pub const VOLUME_CONTRACTION_ENCLOSURE_VERSION: u32 = 1;

/// Native-grid transport envelope at the originally prepared wavelengths.
/// This is not a scene or source bound; source irradiance must be bounded
/// separately, and lanes with contraction >= 1 cannot supply a geometric tail.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct VolumeContractionEnvelope<const N: usize> {
    bounds_m: [Interval; 3],
    diameter_m: f64,
    scattering_integral_upper: [f64; N],
    phase_peak_per_sr: [f64; N],
    contraction_upper: [f64; N],
}
impl<const N: usize> VolumeContractionEnvelope<N> {
    /// Outward world-space support box of the original physical grid.
    #[must_use]
    pub const fn bounds_m(&self) -> [Interval; 3] { self.bounds_m }
    /// Upper bound on the length of any straight chord of the grid, in metres.
    #[must_use]
    pub const fn diameter_m(&self) -> f64 { self.diameter_m }
    /// Bound on integral T*rho*kappa_s ds along any chord, per wavelength.
    #[must_use]
    pub const fn scattering_integral_upper(&self) -> &[f64; N] { &self.scattering_integral_upper }
    /// Upper phase density over every occupied scattering cell and direction.
    #[must_use]
    pub const fn phase_peak_per_sr(&self) -> &[f64; N] { &self.phase_peak_per_sr }
    /// Supremum-norm bound for one further volume event. Values >=1 are not
    /// contractions; they are retained as such, never rounded below one.
    #[must_use]
    pub const fn contraction_upper(&self) -> &[f64; N] { &self.contraction_upper }

    /// Bound orders two and above for one lane, given a finite nonnegative
    /// upper bound on FIRST-order radiance in that lane. The caller must prove
    /// that first-order bound globally, not substitute a sampled camera value.
    /// Exact zero stays zero; positive underflow never claims exact zero.
    #[must_use]
    pub fn tail_upper(&self, lane: usize, first_order_upper: f64) -> Option<f64> {
        let &q = self.contraction_upper.get(lane)?;
        if !first_order_upper.is_finite() || first_order_upper < 0.0 { return None; }
        if first_order_upper == 0.0 || q == 0.0 { return Some(0.0); }
        if !q.is_finite() || q >= 1.0 { return None; }
        let denominator = Interval::point(1.0) - Interval::point(q);
        if denominator.lo() <= 0.0 { return None; }
        let upper = ((Interval::point(first_order_upper) * Interval::point(q)) / denominator).hi();
        (upper.is_finite() && upper > 0.0).then_some(upper)
    }
}

impl<const N: usize> PreparedGridScattering<'_, '_, N> {
    /// Scan the original field once, for reuse by any number of camera rays.
    /// `max_cells` caps complete-field work before scanning. Native coefficient
    /// and phase owners supply all spectra; no field or palette is copied.
    /// One maximum density per declared scattering material avoids repeated phase
    /// evaluation even when cell indices alternate. This small preparation-only
    /// table is capped by `max_cells` as well; queries do not allocate.
    ///
    /// Returns `None` for thermal media or an unrepresentable envelope. This
    /// optional refusal does not invalidate single scattering. Cancellation and
    /// work exhaustion are errors and never return a partial maximum. Existing
    /// absorption-based contraction is combined with a finite-support bound,
    /// so sufficiently thin conservative media can also have a finite tail.
    pub fn volume_contraction_envelope(&self, cx: &Cx<'_>, max_cells: usize)
        -> Result<Option<VolumeContractionEnvelope<N>>, VolumeTransportError> {
        checkpoint(cx)?;
        let grid = &self.extinction.grid.grid;
        if max_cells == 0 || grid.data.len() > max_cells { return Err(VolumeTransportError::Budget); }
        if self.extinction.grid.has_thermal_emission() { return Ok(None); }
        if grid.dims.iter().any(|&n| n as u128 > (1_u128 << 53)) { return Ok(None); }
        let p = Interval::point;
        let extent: [Interval; 3] = std::array::from_fn(|i| p(grid.cell[i]) * p(grid.dims[i] as f64));
        if extent.iter().any(|v| !finite(*v) || v.lo() <= 0.0) { return Ok(None); }
        let bounds_m = std::array::from_fn(|i| Interval::new(grid.origin[i], (p(grid.origin[i])+extent[i]).hi()));
        if bounds_m.iter().any(|v| !finite(*v)) { return Ok(None); }
        // Scale before squaring, including grids with very large SI extents.
        let scale = extent.iter().map(|v| v.hi()).fold(0.0_f64, f64::max);
        let scaled = extent.map(|v| v / p(scale));
        let squares = scaled.map(|v| nonnegative(v*v));
        let diameter_m = (nonnegative(squares[0]+squares[1]+squares[2]).sqrt()*p(scale)).hi();
        if !diameter_m.is_finite() || diameter_m <= 0.0 { return Ok(None); }
        let mut result = VolumeContractionEnvelope {
            bounds_m, diameter_m, scattering_integral_upper: [0.0; N],
            phase_peak_per_sr: [0.0; N], contraction_upper: [0.0; N],
        };
        let Some(medium) = self.extinction.scattering else { return Ok(Some(result)); };
        if !self.active { return Ok(Some(result)); }
        // This existing envelope inspects occupied original absorption/scattering
        // pairs and retained phase integrals. It must not be rebuilt from an
        // invented gray coefficient or assumed exactly normalized table.
        let original_gain = medium.photometric_scattering_gain_upper(cx)?;
        let models = medium.scattering_palette().map_or(1, |m| m.values().len());
        if models > max_cells { return Err(VolumeTransportError::Budget); }
        let mut densities = Vec::new();
        densities.try_reserve_exact(models).map_err(|_| VolumeTransportError::Budget)?;
        densities.resize(models, 0.0_f64);
        for (cell, &density) in grid.data.iter().enumerate() {
            if cell % 256 == 0 { checkpoint(cx)?; }
            let index = medium.scattering_palette().map_or(0, |m| m.cell_indices()[cell] as usize);
            densities[index] = densities[index].max(density);
        }
        let mut max_scattering = [0.0_f64; N];
        for (index, density) in densities.into_iter().enumerate() {
            checkpoint(cx)?;
            if density == 0.0 { continue; }
            let coefficients = match &self.extinction.coefficients {
                Coefficients::Uniform(_) => self.uniform_scattering.as_ref()
                    .ok_or(VolumeTransportError::ScatteringSpectrum)?,
                Coefficients::Indexed { scattering, .. } => &scattering[index],
            };
            if coefficients.iter().all(|v| v.hi() == 0.0) { continue; }
            let law = match medium.scattering_palette() {
                Some(palette) => &palette.values()[index].phase,
                None => medium.phase_at_cell(0)?,
            };
            let spectral = matches!(law, PhaseFunction::SpectralTable(_) | PhaseFunction::SpectralMixture(_));
            let gray = if spectral { None } else { Some(law.enclosure_cosine(cx, Interval::new(-1.0,1.0))?.hi()) };
            for lane in 0..N {
                checkpoint(cx)?;
                if coefficients[lane].hi() == 0.0 { continue; }
                let peak = match gray {
                    Some(upper) => upper,
                    None => law.enclosure_at(cx, self.extinction.wavelengths_nm[lane], Interval::new(-1.0,1.0))?.hi(),
                };
                let sigma = (p(density)*coefficients[lane]).hi();
                if !sigma.is_finite() { return Ok(None); }
                max_scattering[lane] = max_scattering[lane].max(sigma);
                result.phase_peak_per_sr[lane] = result.phase_peak_per_sr[lane].max(peak);
            }
        }
        let solid_angle = p(4.0) * Interval::new(fs_math::next_down(core::f64::consts::PI), fs_math::next_up(core::f64::consts::PI));
        for lane in 0..N {
            checkpoint(cx)?;
            if max_scattering[lane] == 0.0 { continue; }
            let mut chord = (p(max_scattering[lane])*p(diameter_m)).hi().min(1.0);
            // The existing gain uses phase mass >= 1, so it also bounds the
            // unweighted scattering integral. This is a safe extra restriction.
            if let Some(gain) = original_gain { chord = chord.min(gain); }
            result.scattering_integral_upper[lane] = chord;
            let finite_support_gain = (solid_angle*p(result.phase_peak_per_sr[lane])*p(chord)).hi();
            let gain = original_gain.map_or(finite_support_gain, |q| q.min(finite_support_gain));
            if !gain.is_finite() || gain <= 0.0 { return Ok(None); }
            result.contraction_upper[lane] = gain;
        }
        checkpoint(cx)?;
        Ok(Some(result))
    }
}

#[cfg(test)]
mod tests;
