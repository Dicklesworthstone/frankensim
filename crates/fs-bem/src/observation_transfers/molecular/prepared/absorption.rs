//! Molecular heat from the original retained source bank and positive volume rule.
//!
//! For peak exp(-i omega t) phasors, Q = omega Im(1/K) |p|²/2.
//! The original molecular source supplies K at this bank's frequency. Source
//! amplitudes combine before squaring; independent coherence factors combine
//! afterward. Spatial quadrature is explicit and never normalized to boundary
//! work, which also includes outgoing radiation and may include other walls.

use super::{MolecularObservationError, PreparedMolecularResponse, invalid, poll};
use crate::observation_transfers::{
    HarmonicVolumeQuadrature, MAX_SURFACE_QUADRATURE_POINTS, ReceiverSpec, SourceId,
    TransferError, cap, finite,
    harmonic::{MAX_MODES, MAX_PRODUCTS},
};
use fs_blake3::ContentHash;
use fs_material::relaxation::RelaxationState;
use fs_math::{c64::C64, dd::Dd, det};

/// Explicit preparation/storage and per-query arithmetic admission.
#[derive(Clone, Copy, Debug)]
pub struct MolecularAbsorptionLimits {
    pub maximum_quadrature_points: usize,
    pub maximum_coherence_modes: usize,
    /// Conservative product equivalents: R * (modes * (S + 4) + 4).
    /// The diagonal-source path costs R * (2*S + 4), without an S-by-S factor.
    pub maximum_products: usize,
    /// Retained volume weights, density accumulation/output and this owner.
    /// The original response and caller-owned quadrature/factors are borrowed.
    pub maximum_workspace_bytes: usize,
}

/// Complete-cycle molecular heating of the original receiver volume rule.
/// This is a spectral-line mean, not a power spectral density or a time window.
pub struct MolecularAbsorptionReading<'a> {
    /// One value [W/m³] per original receiver, including zero-volume readouts.
    pub absorbed_power_density_w_m3: &'a [f64],
    /// Physical volume integral [W] of the same densities. No wall loss,
    /// outgoing flux, source work or amplitude attenuation is added to it.
    pub absorbed_power_w: f64,
}

/// Repeated absorption queries on one original molecular solution bank.
///
/// Preparation performs no boundary solve. It retains the actual material and
/// frequency through the immutable response and evaluates its caloric law once.
/// Coherent and rank-r queries cost O(S*R) and O(r*S*R); independent-source
/// queries cost O(S*R). Queries allocate nothing and form no dense covariance.
///
/// Repeated nonnegative volume weights coalesce once. Pressure is scaled by
/// sqrt(omega Im(1/K)/2) before squaring, so weak loss does not require an
/// unrepresentable unweighted pressure norm. Fixed-order DD sums retain small
/// positive contributions. No spatial/error certificate or thermal evolution
/// is inferred; original receiver quadrature errors remain separately available
/// from the response. Failure publishes no reading and the owner is reusable.
pub struct PreparedMolecularAbsorption<'a> {
    response: &'a PreparedMolecularResponse,
    quadrature: &'a [HarmonicVolumeQuadrature],
    weights: Vec<Dd>,
    density_sums: Vec<Dd>,
    densities: Vec<f64>,
    root_loss: f64,
    volume_m3: f64,
    limits: MolecularAbsorptionLimits,
    workspace_bytes: usize,
}

impl PreparedMolecularResponse {
    /// Bind physical volume quadrature to THIS source/material/receiver bank.
    /// No independently supplied density, sound speed or loss fit is accepted.
    pub fn prepare_absorption<'a>(
        &'a self,
        quadrature: &'a [HarmonicVolumeQuadrature],
        limits: MolecularAbsorptionLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedMolecularAbsorption<'a>, MolecularObservationError> {
        poll(&mut cancelled)?;
        cap(quadrature.len(), MAX_SURFACE_QUADRATURE_POINTS)?;
        cap(quadrature.len(), limits.maximum_quadrature_points)?;
        cap(limits.maximum_coherence_modes, MAX_MODES)?;
        cap(limits.maximum_products, MAX_PRODUCTS)?;
        if limits.maximum_coherence_modes == 0 || limits.maximum_products == 0 {
            return Err(invalid("molecular absorption requires a nonzero query allowance"));
        }
        let receivers = self.harmonic().receivers().len();
        let sources = self.sources().len();
        admit_products(receivers, sources, 1, limits.maximum_products)?;
        let workspace_bytes = receivers
            .checked_mul(2 * core::mem::size_of::<Dd>() + core::mem::size_of::<f64>())
            .and_then(|v| v.checked_add(core::mem::size_of::<PreparedMolecularAbsorption<'_>>()))
            .ok_or_else(|| invalid("molecular absorption workspace overflows"))?;
        cap(workspace_bytes, limits.maximum_workspace_bytes)?;
        for sample in quadrature {
            poll(&mut cancelled)?;
            if sample.receiver >= receivers
                || !sample.volume_m3.is_finite()
                || sample.volume_m3 < 0.0
            {
                return Err(invalid("molecular absorption volume quadrature is invalid"));
            }
        }
        let omega = self.harmonic().angular_frequency_rad_s();
        let compliance = self
            .source()
            .caloric_response(omega)
            .map_err(TransferError::from)?
            .bulk_modulus_pa
            .recip();
        if !compliance.im.is_finite() || compliance.im < 0.0 {
            return Err(invalid("molecular absorption compliance is nonpassive"));
        }
        let root_loss = det::sqrt(0.5) * det::sqrt(omega) * det::sqrt(compliance.im);
        if !root_loss.is_finite()
            || root_loss < 0.0
            || (compliance.im > 0.0 && root_loss == 0.0)
        {
            return Err(invalid("molecular absorption coefficient is unrepresentable"));
        }
        let mut weights = allocated(receivers, Dd::ZERO)?;
        for sample in quadrature {
            poll(&mut cancelled)?;
            weights[sample.receiver] =
                weights[sample.receiver] + Dd::from_f64(sample.volume_m3);
        }
        let mut volume = Dd::ZERO;
        for &weight in &weights {
            if !dd_finite(weight) {
                return Err(invalid("molecular absorption volume is unrepresentable"));
            }
            volume = volume + weight;
        }
        let volume_m3 = volume.to_f64();
        if !dd_finite(volume) || !volume_m3.is_finite() || volume_m3 < 0.0 {
            return Err(invalid("molecular absorption total volume is unrepresentable"));
        }
        let density_sums = allocated(receivers, Dd::ZERO)?;
        let densities = allocated(receivers, 0.0)?;
        poll(&mut cancelled)?;
        Ok(PreparedMolecularAbsorption {
            response: self,
            quadrature,
            weights,
            density_sums,
            densities,
            root_loss,
            volume_m3,
            limits,
            workspace_bytes,
        })
    }
}

impl<'a> PreparedMolecularAbsorption<'a> {
    pub fn response(&self) -> &'a PreparedMolecularResponse {
        self.response
    }
    pub fn source(&self) -> &RelaxationState {
        self.response.source()
    }
    pub fn sources(&self) -> &[SourceId] {
        self.response.sources()
    }
    pub fn receivers(&self) -> &[ReceiverSpec] {
        self.response.harmonic().receivers()
    }
    pub fn quadrature(&self) -> &'a [HarmonicVolumeQuadrature] {
        self.quadrature
    }
    pub fn volume_m3(&self) -> f64 {
        self.volume_m3
    }
    pub fn angular_frequency_rad_s(&self) -> f64 {
        self.response.harmonic().angular_frequency_rad_s()
    }
    pub fn producer_hash(&self) -> ContentHash {
        self.response.harmonic().producer_hash()
    }
    pub fn workspace_bytes(&self) -> usize {
        self.workspace_bytes
    }

    /// One coherent source vector in the original source coordinates and phase.
    pub fn evaluate_coherent(
        &mut self,
        rates: &[C64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<MolecularAbsorptionReading<'_>, MolecularObservationError> {
        self.evaluate_coherence_modes(rates, 1, cancelled)
    }

    /// Mode-major, source-minor factors B with peak second moment C = B B*.
    /// Zero modes explicitly mean silence. Admitted factors are never turned
    /// into an incoherent sum of individual source intensities.
    pub fn evaluate_coherence_modes(
        &mut self,
        factors: &[C64],
        modes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MolecularAbsorptionReading<'_>, MolecularObservationError> {
        poll(&mut cancelled)?;
        cap(modes, self.limits.maximum_coherence_modes)?;
        let sources = self.sources().len();
        let entries = modes
            .checked_mul(sources)
            .ok_or_else(|| invalid("molecular absorption factor extent overflows"))?;
        if factors.len() != entries {
            return Err(invalid("molecular absorption factor source order or shape differs"));
        }
        admit_products(
            self.receivers().len(),
            sources,
            modes,
            self.limits.maximum_products,
        )?;
        for &factor in factors {
            poll(&mut cancelled)?;
            if !finite(factor) {
                return Err(invalid("molecular absorption factors must be finite"));
            }
        }
        self.density_sums.fill(Dd::ZERO);
        if self.root_loss != 0.0 {
            let field = self.response.harmonic();
            for (receiver, sum) in self.density_sums.iter_mut().enumerate() {
                poll(&mut cancelled)?;
                for factor in factors.chunks_exact(sources) {
                    poll(&mut cancelled)?;
                    let pressure = field.pressure_field(factor, receiver)?;
                    *sum = *sum + Dd::from_f64(heat_density(pressure, self.root_loss)?);
                }
            }
        }
        self.finish(&mut cancelled)
    }

    /// Diagonal source peak second moments in their original squared SI units.
    /// Uses one retained pressure column per source, with no dense diagonal
    /// factor allocation and no quadratic growth in the source count.
    pub fn evaluate_independent_sources(
        &mut self,
        peak_second_moments: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MolecularAbsorptionReading<'_>, MolecularObservationError> {
        poll(&mut cancelled)?;
        let sources = self.sources().len();
        if peak_second_moments.len() != sources {
            return Err(invalid("independent molecular source powers have the wrong shape"));
        }
        let products = self
            .receivers()
            .len()
            .checked_mul(
                sources
                    .checked_mul(2)
                    .and_then(|v| v.checked_add(4))
                    .ok_or_else(|| invalid("molecular absorption work overflows"))?,
            )
            .ok_or_else(|| invalid("molecular absorption work overflows"))?;
        cap(products, self.limits.maximum_products)?;
        for &power in peak_second_moments {
            poll(&mut cancelled)?;
            if !power.is_finite() || power < 0.0 {
                return Err(invalid("independent molecular source powers must be nonnegative"));
            }
        }
        self.density_sums.fill(Dd::ZERO);
        if self.root_loss != 0.0 {
            let field = self.response.harmonic();
            for (source, &power) in peak_second_moments.iter().enumerate() {
                poll(&mut cancelled)?;
                let scale = det::sqrt(power);
                for (receiver, sum) in self.density_sums.iter_mut().enumerate() {
                    poll(&mut cancelled)?;
                    let pressure = field.pressure_transfer(source, receiver).scale(scale);
                    if !finite(pressure) {
                        return Err(invalid("independent molecular pressure is unrepresentable"));
                    }
                    *sum = *sum + Dd::from_f64(heat_density(pressure, self.root_loss)?);
                }
            }
        }
        self.finish(&mut cancelled)
    }

    fn finish(
        &mut self,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MolecularAbsorptionReading<'_>, MolecularObservationError> {
        let mut power = Dd::ZERO;
        for ((density, &sum), &weight) in self
            .densities
            .iter_mut()
            .zip(&self.density_sums)
            .zip(&self.weights)
        {
            poll(cancelled)?;
            *density = sum.to_f64();
            if !dd_finite(sum) || !density.is_finite() || *density < 0.0 {
                return Err(invalid("molecular absorbed heat density is unrepresentable"));
            }
            let contribution = sum * weight;
            if !dd_finite(contribution)
                || (*density > 0.0 && weight.to_f64() > 0.0 && contribution == Dd::ZERO)
            {
                return Err(invalid("molecular absorbed volume power is unrepresentable"));
            }
            power = power + contribution;
        }
        let absorbed_power_w = power.to_f64();
        if !dd_finite(power) || !absorbed_power_w.is_finite() || absorbed_power_w < 0.0 {
            return Err(invalid("molecular absorbed power is unrepresentable"));
        }
        poll(cancelled)?;
        Ok(MolecularAbsorptionReading {
            absorbed_power_density_w_m3: &self.densities,
            absorbed_power_w,
        })
    }
}

fn heat_density(pressure: C64, root_loss: f64) -> Result<f64, MolecularObservationError> {
    let scaled = pressure.scale(root_loss);
    let heat = scaled.norm_sq();
    if !finite(scaled)
        || !heat.is_finite()
        || (root_loss > 0.0 && pressure != C64::ZERO && heat == 0.0)
    {
        return Err(invalid("molecular heat density is unrepresentable"));
    }
    Ok(heat)
}

fn admit_products(
    receivers: usize,
    sources: usize,
    modes: usize,
    maximum: usize,
) -> Result<(), MolecularObservationError> {
    let products = sources
        .checked_add(4)
        .and_then(|v| v.checked_mul(modes))
        .and_then(|v| v.checked_add(4))
        .and_then(|v| v.checked_mul(receivers))
        .ok_or_else(|| invalid("molecular absorption work overflows"))?;
    cap(products, maximum)?;
    Ok(())
}

fn dd_finite(value: Dd) -> bool {
    value.hi.is_finite() && value.lo.is_finite()
}

fn allocated<T: Clone>(count: usize, value: T) -> Result<Vec<T>, MolecularObservationError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(count)
        .map_err(|_| invalid("molecular absorption allocation failed"))?;
    values.resize(count, value);
    Ok(values)
}

#[cfg(test)]
mod tests;
