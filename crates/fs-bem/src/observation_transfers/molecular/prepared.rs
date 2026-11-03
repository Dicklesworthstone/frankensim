//! Matched source loading and receiver fields from one molecular solve bank.
//!
//! The original boundary owner supplies the pressures AND actual velocities.
//! Their work projections retain source interference and wall/exterior energy
//! separately; the observation bank reads those same solutions. This is a
//! sampled physical load, not a fitted or certified passive realization.

use super::{
    MolecularBemObservation, MolecularObservationError, MolecularReceiverError,
    MolecularReceiverQuadrature, invalid,
};
use crate::helmholtz::molecular::MolecularRadiation;
use crate::observation_transfers::{
    MAX_OBSERVATION_ENTRIES, PreparedHarmonicObservation, SourceId, TransferError, cap, finite,
    harmonic::{self, Sum},
};
use fs_material::relaxation::RelaxationState;
use fs_math::c64::C64;

/// Additional live-storage and projection-work limits for the complete source
/// response. The original `ComplexOptions` still bounds each BEM assembly.
#[derive(Debug, Clone, Copy)]
pub struct MolecularResponseLimits {
    /// Maximum live complex entries outside the original BEM assembly. The
    /// conservative reservation is `3*s*s + 2*s*p + 8*s*r`: three source-work
    /// matrices, the pressure/velocity solution bank, and simultaneous raw and
    /// converted receiver banks. Adaptive receivers add `2*s*r` for four real
    /// discrepancy components per source/receiver; their row storage has a
    /// separate MolecularNearFieldOptions cap. Here s/p/r are sources/panels/receivers.
    pub maximum_retained_entries: usize,
    /// Three weighted products per source pair and panel, including the wall
    /// work projection. Receiver work retains its original owner's admission.
    pub maximum_panel_products: u64,
}

/// Mean work for coherent peak source-rate phasors, in watts. Exterior work
/// includes eventual molecular heat and outgoing flux; it is not the power
/// crossing an arbitrarily distant receiver sphere.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct MolecularResponsePower {
    pub driving_power_w: f64,
    pub exterior_power_w: f64,
    pub absorbed_wall_power_w: f64,
    /// Raw input minus exterior minus wall work. No channel is clipped or
    /// reconstructed to force this residual to zero.
    pub balance_residual_w: f64,
}

/// Immutable complete source-load matrix and matched harmonic receiver bank.
/// Mechanical coordinates retain force/velocity work pairs; aperture-flow
/// coordinates retain pressure/volume-flow pairs. Cross terms are never
/// discarded, symmetrized or repaired to imply continuum passivity.
pub struct PreparedMolecularResponse {
    source: RelaxationState,
    harmonic: PreparedHarmonicObservation,
    // The already-admitted 2*s*p source bank also supplies exact coherent
    // boundary traces for incident/motion interference and interface work.
    boundary_solutions: Vec<MolecularRadiation>,
    receiver_quadrature: Option<MolecularReceiverQuadrature>,
    impedance: Vec<C64>,
    exterior_work: Vec<C64>,
    wall_work: Vec<C64>,
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), MolecularObservationError> {
    if cancelled() {
        Err(TransferError::Cancelled.into())
    } else {
        Ok(())
    }
}

impl MolecularBemObservation<'_> {
    /// Prepare all source reaction rows and receivers from the SAME complete
    /// solution bank. Each existing bounded source batch shares one assembly
    /// and LU; neither additional reaction rows nor receivers re-solve it.
    ///
    /// For source traces b_i, solved pressures p_j and actual velocities v_i,
    /// Z_ij = integral conj(b_i) p_j dA. The exterior and wall work matrices
    /// instead test with v_i and b_i-v_i. Direct wall projection preserves
    /// weak absorption that subtraction of two large total powers can lose.
    /// The three raw matrices remain separately available for diagnostics.
    ///
    /// Source identity, molecular reservoirs, wall laws, positive frequency,
    /// geometric limits and cancellations retain their original admission.
    /// This does not fit a time-domain load or establish spatial convergence.
    pub fn prepare_response(
        &self,
        omega: f64,
        limits: MolecularResponseLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedMolecularResponse, MolecularObservationError> {
        poll(&mut cancelled)?;
        harmonic::admit_frequency(omega, self.boundary.medium)?;
        let n = self.boundary.sources.len();
        let panels = self.boundary.surface.areas().len();
        let receivers = self.boundary.receivers.len();
        let entries = n
            .checked_mul(n)
            .ok_or_else(|| invalid("molecular source matrix extent overflow"))?;
        let samples = n
            .checked_mul(receivers)
            .and_then(|v| v.checked_mul(4))
            .ok_or_else(|| invalid("molecular receiver bank extent overflow"))?;
        let live_entries = entries
            .checked_mul(3)
            .and_then(|v| {
                n.checked_mul(panels)
                    .and_then(|p| p.checked_mul(2))
                    .and_then(|p| v.checked_add(p))
            })
            .and_then(|v| samples.checked_mul(2).and_then(|s| v.checked_add(s)))
            .and_then(|v| {
                v.checked_add(if self.near_field_options.is_some() {
                    samples / 2
                } else {
                    0
                })
            })
            .ok_or_else(|| invalid("molecular source response storage overflow"))?;
        cap(
            live_entries,
            limits.maximum_retained_entries.min(MAX_OBSERVATION_ENTRIES),
        )?;
        let products = (entries as u64)
            .checked_mul(panels as u64)
            .and_then(|v| v.checked_mul(3))
            .ok_or_else(|| invalid("molecular source projection work overflow"))?;
        if products > limits.maximum_panel_products {
            return Err(invalid(
                "molecular source projection exceeds its panel-product budget",
            ));
        }
        let solved = self.solutions(omega, &mut cancelled)?;
        let mut impedance = vec![C64::ZERO; entries];
        let mut exterior_work = vec![C64::ZERO; entries];
        let mut wall_work = vec![C64::ZERO; entries];
        for row in 0..n {
            let drive = &self.boundary.fields[row];
            let velocity = solved[row].velocity_m_s();
            for column in 0..n {
                poll(&mut cancelled)?;
                let pressure = solved[column].pressure_pa();
                let mut sums = [Sum::default(); 3];
                for panel in 0..panels {
                    if panel % 128 == 0 {
                        poll(&mut cancelled)?;
                    }
                    let weighted = pressure[panel].scale(self.boundary.surface.areas()[panel]);
                    sums[0].add(drive[panel].conj() * weighted);
                    sums[1].add(velocity[panel].conj() * weighted);
                    sums[2].add((drive[panel] - velocity[panel]).conj() * weighted);
                }
                let values = sums.map(Sum::finish);
                if values.iter().any(|&z| !finite(z)) {
                    return Err(invalid("nonfinite molecular source work projection"));
                }
                let index = row * n + column;
                impedance[index] = values[0];
                exterior_work[index] = values[1];
                wall_work[index] = values[2];
            }
        }
        let (fields, receiver_quadrature) = self.observation_samples(&solved, &mut cancelled)?;
        let harmonic = PreparedHarmonicObservation::from_samples(
            omega,
            self.boundary.medium,
            &self.boundary.sources,
            &self.boundary.receivers,
            self.operator_hash(),
            fields,
        )?
        .without_radiation_momentum();
        poll(&mut cancelled)?;
        Ok(PreparedMolecularResponse {
            source: self.source().clone(),
            harmonic,
            boundary_solutions: solved,
            receiver_quadrature,
            impedance,
            exterior_work,
            wall_work,
        })
    }
}

impl PreparedMolecularResponse {
    /// Source-column quadrature diagnostics, present only for adaptive rows.
    pub fn receiver_quadrature(&self) -> Option<&MolecularReceiverQuadrature> {
        self.receiver_quadrature.as_ref()
    }

    /// Sum absolute source-rate-weighted quadrature estimates in receiver order.
    /// Coherent field amplitudes still come from the original harmonic bank;
    /// these positive discrepancy estimates do not become field/power terms.
    pub fn receiver_quadrature_for_rates(
        &self,
        rates: &[C64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Option<Vec<MolecularReceiverError>>, MolecularObservationError> {
        poll(&mut cancelled)?;
        self.admit_rates(rates)?;
        let Some(quadrature) = &self.receiver_quadrature else {
            return Ok(None);
        };
        let receivers = self.harmonic.receivers().len();
        let mut output = Vec::with_capacity(receivers);
        for receiver in 0..receivers {
            poll(&mut cancelled)?;
            let mut sums = [Sum::default(); 4];
            for (source, &rate) in rates.iter().enumerate() {
                if source % 128 == 0 {
                    poll(&mut cancelled)?;
                }
                let error = quadrature.source_fields[source * receivers + receiver];
                let magnitude = rate.abs();
                sums[0].add(C64::from_re(magnitude * error.pressure_pa));
                for axis in 0..3 {
                    sums[axis + 1].add(C64::from_re(magnitude * error.particle_velocity_m_s[axis]));
                }
            }
            let errors = sums.map(|sum| sum.finish().re);
            if errors.iter().any(|value| !value.is_finite()) {
                return Err(invalid("molecular coherent receiver discrepancy overflow"));
            }
            output.push(MolecularReceiverError {
                pressure_pa: errors[0],
                particle_velocity_m_s: [errors[1], errors[2], errors[3]],
            });
        }
        poll(&mut cancelled)?;
        Ok(Some(output))
    }

    pub fn source(&self) -> &RelaxationState {
        &self.source
    }

    pub fn sources(&self) -> &[SourceId] {
        self.harmonic.sources()
    }

    /// Matched pressure/particle-velocity transfers, including coherent and
    /// partially coherent observations, intensity and microphone cross spectra.
    /// Evaluation never repeats a boundary solve or attenuation correction.
    pub fn harmonic(&self) -> &PreparedHarmonicObservation {
        &self.harmonic
    }

    /// Number of original boundary panels, in the original surface order.
    pub fn panel_count(&self) -> usize {
        self.boundary_solutions
            .first()
            .map_or(0, |solution| solution.pressure_pa().len())
    }

    /// Combine the retained original pressure and ACTUAL normal-velocity
    /// traces for coherent source rates. The latter includes any passive wall
    /// response; it is not reconstructed from the prescribed drive alone.
    ///
    /// This allocation-free query performs 2*s*p complex products, bounded by
    /// the already-admitted source bank and its cold projection-work ceiling.
    /// No additional assembly, solve, attenuation or pressure-to-velocity
    /// approximation occurs. It preserves the boundary fields needed to form
    /// incident/motion interference directly instead of subtracting powers.
    /// Discard both output buffers on refusal: a late cancellation or overflow
    /// may leave a partial traversal, but never mutates the prepared bank.
    pub fn evaluate_boundary_coherent(
        &self,
        rates: &[C64],
        pressure_pa: &mut [C64],
        normal_velocity_m_s: &mut [C64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), MolecularObservationError> {
        poll(&mut cancelled)?;
        self.admit_rates(rates)?;
        let panels = self.panel_count();
        if pressure_pa.len() != panels || normal_velocity_m_s.len() != panels {
            return Err(invalid(
                "molecular boundary output must retain original panel order and extent",
            ));
        }
        for panel in 0..panels {
            poll(&mut cancelled)?;
            let mut pressure = Sum::default();
            let mut velocity = Sum::default();
            for (source, (&rate, solution)) in
                rates.iter().zip(&self.boundary_solutions).enumerate()
            {
                if source % 128 == 0 {
                    poll(&mut cancelled)?;
                }
                pressure.add(solution.pressure_pa()[panel] * rate);
                velocity.add(solution.velocity_m_s()[panel] * rate);
            }
            let (p, v) = (pressure.finish(), velocity.finish());
            if !finite(p) || !finite(v) {
                return Err(invalid("nonfinite coherent molecular boundary trace"));
            }
            pressure_pa[panel] = p;
            normal_velocity_m_s[panel] = v;
        }
        poll(&mut cancelled)
    }

    /// Complete row-major resisting-effort/source-rate matrix. The force ON a
    /// driven structure is the negative of this load's resisting effort.
    pub fn impedance(&self) -> &[C64] {
        &self.impedance
    }

    /// Pressure/volume-flow impedance [Pa s/m³] is available only when every
    /// original source is a unit-integral aperture-flow coordinate. Its order
    /// is exactly `sources()`, suitable for a labeled duct multiport closure.
    pub fn acoustic_impedance(&self) -> Result<&[C64], MolecularObservationError> {
        if self
            .sources()
            .iter()
            .any(|source| !matches!(source, SourceId::ApertureVolumeFlow(_)))
        {
            return Err(invalid(
                "acoustic impedance requires only aperture-volume-flow coordinates",
            ));
        }
        Ok(&self.impedance)
    }

    /// E_ij = integral conj(v_i) p_j dA, in original source coordinates.
    pub fn exterior_work_matrix(&self) -> &[C64] {
        &self.exterior_work
    }

    /// W_ij = integral conj(b_i-v_i) p_j dA. Rigid walls give exact zeros;
    /// passive walls retain their raw collocation work without eigenvalue repair.
    pub fn wall_work_matrix(&self) -> &[C64] {
        &self.wall_work
    }

    fn admit_rates(&self, rates: &[C64]) -> Result<(), MolecularObservationError> {
        if rates.len() != self.sources().len() || rates.iter().any(|&z| !finite(z)) {
            return Err(invalid(
                "molecular source rates do not match the original finite coordinates",
            ));
        }
        Ok(())
    }

    /// Allocation-free O(sources²) resisting effort. Discard the caller's
    /// output buffer on refusal; it may contain an incomplete traversal.
    pub fn resisting_efforts(
        &self,
        rates: &[C64],
        output: &mut [C64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), MolecularObservationError> {
        poll(&mut cancelled)?;
        self.admit_rates(rates)?;
        if output.len() != rates.len() {
            return Err(invalid("molecular resisting-effort output shape"));
        }
        for (row, out) in self.impedance.chunks_exact(rates.len()).zip(output) {
            poll(&mut cancelled)?;
            let mut sum = Sum::default();
            for (&z, &rate) in row.iter().zip(rates) {
                sum.add(z * rate);
            }
            *out = sum.finish();
            if !finite(*out) {
                return Err(invalid("nonfinite molecular resisting effort"));
            }
        }
        poll(&mut cancelled)
    }

    /// Allocation-free O(sources²) coherent mean work with all self/mutual
    /// terms. Negative channel values remain signed: this query grants no
    /// discrete passivity certificate and performs no positive clipping.
    pub fn power_w(
        &self,
        rates: &[C64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MolecularResponsePower, MolecularObservationError> {
        poll(&mut cancelled)?;
        self.admit_rates(rates)?;
        let n = rates.len();
        let mut powers = [Sum::default(); 3];
        for i in 0..n {
            poll(&mut cancelled)?;
            let mut rows = [Sum::default(); 3];
            for (j, &rate) in rates.iter().enumerate() {
                let index = i * n + j;
                rows[0].add(self.impedance[index] * rate);
                rows[1].add(self.exterior_work[index] * rate);
                rows[2].add(self.wall_work[index] * rate);
            }
            for (power, row) in powers.iter_mut().zip(rows) {
                power.add(rates[i].conj() * row.finish());
            }
        }
        let totals = powers.map(Sum::finish);
        if totals.iter().any(|&value| !finite(value)) {
            return Err(invalid("nonfinite molecular coherent work"));
        }
        let values = totals.map(|sum| 0.5 * sum.re);
        let balance_residual_w = values[0] - values[1] - values[2];
        if values.iter().any(|value| !value.is_finite()) || !balance_residual_w.is_finite() {
            return Err(invalid("nonfinite molecular coherent work"));
        }
        poll(&mut cancelled)?;
        Ok(MolecularResponsePower {
            driving_power_w: values[0],
            exterior_power_w: values[1],
            absorbed_wall_power_w: values[2],
            balance_residual_w,
        })
    }
}

#[cfg(test)]
mod tests;
