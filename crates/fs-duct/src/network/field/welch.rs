//! Measured source cross spectra through the original interior p/U transfers.

use super::*;
use fs_fft::welch::WelchCrossSpectrum;

/// One positive-bin contribution at an original spatial probe. These are
/// bin-integrated second moments, not PSD/Hz or finite-shutter first moments.
/// Flow and flux point in the original edge's from-to direction.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NetworkFieldBinStatistics {
    pub mean_square_pressure_pa2: f64,
    pub mean_square_volume_flow_m6_s2: f64,
    /// Cross-section mean axial particle velocity, not a resolved radial field.
    pub mean_square_particle_velocity_m2_s2: f64,
    /// Complex 1/2 E[p conj(U)] [W]. Its real part is transported power;
    /// the imaginary part retains pressure/flow quadrature at standing waves.
    pub pressure_flow_moment_w: C64,
    pub mean_power_w: f64,
    pub mean_intensity_w_m2: f64,
}

impl PreparedNetworkFields {
    /// Observe a recorded interior Welch bin at every existing physical probe.
    /// Exact source labels, pressure/flow SI quantities, order and frequency
    /// are checked against the original loaded network. Each borrowed segment
    /// factor receives sqrt(2*df) once, is propagated coherently in DD, and only
    /// then contributes its half-weighted field moments. Pressure antinodes
    /// therefore need not carry traveling-wave energy or nonzero volume flow.
    ///
    /// No network solve, new probe, material interpolation, source covariance
    /// matrix, independent-realization assumption or cross-frequency coherence
    /// is introduced. The record's leakage/resolution limits remain. Local
    /// signed power is transport plus storage exchange, not dissipated heat.
    ///
    /// Work counts R*(2*P*C + 3*P + 2*C) complex-product equivalents for R
    /// segments, P original probes and C source coordinates, with a 262144
    /// ceiling and at most 4096 segments. Only four DD sums per probe and the
    /// returned readings are allocated, bounded by the original 1024-probe
    /// admission. Factors stream through a fixed 32-source stack buffer.
    /// Cancellation returns no partially published readings.
    pub fn observe_welch_bin(
        &self,
        spectrum: &WelchCrossSpectrum,
        bin: usize,
        recorded_order: &[NetworkSourceCoordinate],
        budget: NetworkCoherenceBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<NetworkFieldBinStatistics>, DuctError> {
        poll(&mut cancelled, 0)?;
        let view = self.ports.recorded_bin(spectrum, bin, recorded_order)?;
        let (r, p, c) = (
            view.segment_count(),
            self.locations.len(),
            view.channel_count(),
        );
        let products = p
            .checked_mul(c)
            .and_then(|v| v.checked_mul(2))
            .and_then(|v| p.checked_mul(3).and_then(|q| v.checked_add(q)))
            .and_then(|v| c.checked_mul(2).and_then(|q| v.checked_add(q)))
            .and_then(|v| v.checked_mul(r))
            .ok_or_else(|| bad("recorded spatial work extent overflow"))?;
        if r > 4096
            || r > budget.maximum_rank
            || products > 262_144
            || products > budget.maximum_complex_products
        {
            return Err(bad(
                "recorded spatial observation exceeds its factor/work budget",
            ));
        }
        let mut sums = Vec::new();
        sums.try_reserve_exact(p)
            .map_err(|_| bad("recorded spatial moment allocation failed"))?;
        sums.resize(p, [Dd::ZERO; 4]);
        let mut excitation = [C64::ZERO; MAX_NETWORK_PORTS];
        for (segment, factor) in view.factors().enumerate() {
            poll(&mut cancelled, segment)?;
            for (out, &value) in excitation[..c].iter_mut().zip(factor) {
                let value = view
                    .peak_factor(value)
                    .map_err(|_| bad("recorded network source peak is unrepresentable"))?;
                *out = C64::new(value.re, value.im);
            }
            for (probe, sum) in sums.iter_mut().enumerate() {
                poll(&mut cancelled, probe)?;
                let [[pr, pi], [qr, qi]] = self.parts(probe, &excitation[..c])?;
                sum[0] = sum[0] + half_product(pr, pr) + half_product(pi, pi);
                sum[1] = sum[1] + half_product(qr, qr) + half_product(qi, qi);
                sum[2] = sum[2] + half_product(pr, qr) + half_product(pi, qi);
                sum[3] = sum[3] + half_product(pi, qr) - half_product(pr, qi);
            }
        }
        let mut result = Vec::new();
        result
            .try_reserve_exact(p)
            .map_err(|_| bad("recorded spatial result allocation failed"))?;
        for (probe, (sum, location)) in sums.into_iter().zip(&self.locations).enumerate() {
            poll(&mut cancelled, probe)?;
            let [pp, qq, pq_re, pq_im] = sum.map(Dd::to_f64);
            let velocity = (qq / location.area_m2) / location.area_m2;
            let intensity = pq_re / location.area_m2;
            if pp < 0.0
                || qq < 0.0
                || [pp, qq, pq_re, pq_im, velocity, intensity]
                    .iter()
                    .any(|v| !v.is_finite())
            {
                return Err(DuctError::Singular);
            }
            result.push(NetworkFieldBinStatistics {
                mean_square_pressure_pa2: pp,
                mean_square_volume_flow_m6_s2: qq,
                mean_square_particle_velocity_m2_s2: velocity,
                pressure_flow_moment_w: C64::new(pq_re, pq_im),
                mean_power_w: pq_re,
                mean_intensity_w_m2: intensity,
            });
        }
        poll(&mut cancelled, r)?;
        Ok(result)
    }
}

// The physical half is applied once, to the larger operand first. This avoids
// both an avoidable overflow before halving and loss of a subnormal small
// operand whose product with the other field component remains representable.
fn half_product(a: Dd, b: Dd) -> Dd {
    if a.hi.abs() >= b.hi.abs() {
        (a * Dd::from_f64(0.5)) * b
    } else {
        a * (b * Dd::from_f64(0.5))
    }
}

#[cfg(test)]
mod tests;
