//! Finite-band pressure AND particle velocity from the original boundary solve.
//!
//! Native exp(-i omega t) pressure/gradient samples are converted by inviscid
//! momentum, u = grad(p)/(i rho omega), BEFORE the existing Fourier realization.
//! One native solve, frequency cache, FFT ladder and compensated replay serve
//! all four physical channels. No time integration of a replayed gradient,
//! plane-wave substitution, absorption overlay or extra delay is introduced.
//!
//! Every producer must explicitly supply a positive lower observation window.
//! In particular, a steady volume flow can have zero acoustic pressure and
//! nonzero velocity; assigning zero DC to that velocity would invent physics.
//! This owner instead realizes the declared finite observation band, with the
//! original upper taper and reconstruction latency. DC and the omitted band
//! are not reconstructed. Zero input prehistory is an explicit replay condition.
//!
//! Receivers remain at their original fixed positions. Their ideal FirstOrder
//! orientation may change at each accepted frame without resetting the source
//! history. Instantaneous p*u is signed and has no harmonic factor of one half.
//! FIR convergence does not establish source, material or spatial accuracy;
//! observation remains separate from work-conjugate source loading.

use super::{
    ContentHash, DomainHasher, FirKernel, LadderRung, LowerBandWindow, Medium, ObservationFir,
    ObservationTransferBake, ProducerUncertainty, ReceiverSpec, SourceId, TransferDeclaration,
    TransferError, bake_observation_transfers_with_errors, finite, harmonic,
};
use crate::near_field::FirstOrder;
use fs_math::c64::C64;
use std::cell::RefCell;

#[cfg(test)]
mod tests;

/// One real, co-located acoustic field in the original Cartesian frame.
/// Pressure is the acoustic perturbation, not the absolute ambient pressure.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct PressureVelocityFrame {
    pub pressure_pa: f64,
    pub particle_velocity_m_s: [f64; 3],
}

impl PressureVelocityFrame {
    /// Signed instantaneous acoustic intensity [W/m^2]. Near-field storage
    /// and interference may reverse it; no clipping or peak/RMS conversion.
    /// This readout changes no source or observer state if it overflows.
    pub fn instantaneous_intensity_w_m2(self) -> Result<[f64; 3], TransferError> {
        let value = self.particle_velocity_m_s.map(|u| self.pressure_pa * u);
        if !self.pressure_pa.is_finite()
            || self.particle_velocity_m_s.iter().any(|v| !v.is_finite())
            || value.iter().any(|v| !v.is_finite())
        {
            return Err(TransferError::Invalid(
                "nonfinite instantaneous pressure/velocity intensity",
            ));
        }
        Ok(value)
    }
}

/// An immutable typed p/u realization of one original finite-boundary owner.
///
/// The numeric payload stays private: it must never be passed to a consumer
/// interpreting channels 1..3 as pressure gradients. Its field identity is
/// distinct from the original p/grad declaration and includes the explicit
/// lower window and the source-bound microphone calibration medium.
#[derive(Clone, Debug)]
pub struct PressureVelocityBake {
    numeric: ObservationTransferBake,
    medium: Medium,
    lower_window: LowerBandWindow,
    producer_hash: ContentHash,
    identity: ContentHash,
}

impl PressureVelocityBake {
    /// Complete p/u realization identity, including its physical producer.
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Original material/geometry/operator identity, before field realization.
    pub const fn producer_hash(&self) -> ContentHash {
        self.producer_hash
    }

    /// Original real rho/c calibration. For molecular propagation this is the
    /// equilibrium gas, not a replacement for its dispersive wave impedance.
    pub const fn medium(&self) -> Medium {
        self.medium
    }

    pub fn sources(&self) -> &[SourceId] {
        &self.numeric.sources
    }

    pub fn receivers(&self) -> &[ReceiverSpec] {
        &self.numeric.receivers
    }

    pub const fn sample_rate_hz(&self) -> f64 {
        self.numeric.fs_hz
    }

    /// The original conservative geometric prefixes, source-major then receiver.
    /// Dispersion and fractional delays are already present in the same kernels.
    pub fn geometric_delays_s(&self) -> &[f64] {
        &self.numeric.geometric_delays_s
    }

    pub const fn reconstruction_latency_samples(&self) -> usize {
        self.numeric.filter_latency_samples
    }

    pub const fn lower_band_window(&self) -> LowerBandWindow {
        self.lower_window
    }

    pub fn pass_band_hz(&self) -> [f64; 2] {
        [self.lower_window.pass_hz(), self.numeric.pass_band_upper_hz]
    }

    pub const fn upper_taper_stop_hz(&self) -> f64 {
        self.numeric.taper_stop_hz
    }

    pub const fn producer_uncertainty(&self) -> ProducerUncertainty {
        self.numeric.producer_uncertainty
    }

    /// Original measured ladder, now with pressure and the three-component
    /// VELOCITY treated as separate physical quantities. Pa and m/s never
    /// share an error denominator; all three velocity directions share one.
    pub fn ladder(&self) -> &[LadderRung] {
        &self.numeric.ladder
    }

    /// Source-major, receiver-major p [Pa], ux/uy/uz [m/s] per original SI
    /// source coordinate. Real FIR transfer_at uses the native DFT convention.
    /// No p/gradient artifact or mutable kernel set is exposed by this type.
    pub fn kernels(&self) -> &[FirKernel] {
        &self.numeric.kernels
    }

    /// Create one shared source history with explicit zero prehistory. Every
    /// source uses its declared coordinate and sample rate, without resampling.
    /// The wrapper adds at most two complete output buffers (4096 bytes under
    /// the existing 256-scalar output cap); original FIR/history caps remain.
    pub fn observer(&self) -> Result<PressureVelocityObserver<'_>, TransferError> {
        Ok(PressureVelocityObserver {
            bake: self,
            fir: self.numeric.observer()?,
            output: vec![0.0; 4 * self.receivers().len()],
            frames: vec![PressureVelocityFrame::default(); self.receivers().len()],
            accepted_samples: 0,
        })
    }
}

/// Allocation-free p/u replay and changing ideal microphone direction.
/// This borrows its immutable source/receiver/material declaration throughout.
#[derive(Clone)]
pub struct PressureVelocityObserver<'a> {
    bake: &'a PressureVelocityBake,
    fir: ObservationFir,
    output: Vec<f64>,
    frames: Vec<PressureVelocityFrame>,
    accepted_samples: usize,
}

impl<'a> PressureVelocityObserver<'a> {
    pub const fn bake(&self) -> &'a PressureVelocityBake {
        self.bake
    }

    pub const fn accepted_sample_count(&self) -> usize {
        self.accepted_samples
    }

    /// No observation is reported before the first accepted input sample.
    pub fn current_fields(&self) -> Option<&[PressureVelocityFrame]> {
        (self.accepted_samples != 0).then_some(self.frames.as_slice())
    }

    /// Advance exactly one sample in the immutable source order. Cancellation
    /// is admitted at entry; the bounded original FIR then completes atomically.
    /// Invalid input, overflow or entry cancellation changes neither history,
    /// accepted count nor the previously visible complete field. No allocation
    /// or fallible operation follows the native FIR's successful commit.
    pub fn step(
        &mut self,
        input: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<&[PressureVelocityFrame], TransferError> {
        if cancelled() {
            return Err(TransferError::Cancelled);
        }
        let next = self
            .accepted_samples
            .checked_add(1)
            .ok_or(TransferError::Invalid(
                "pressure/velocity sample count overflow",
            ))?;
        self.fir.step_into(input, &mut self.output)?;
        for (frame, row) in self.frames.iter_mut().zip(self.output.chunks_exact(4)) {
            *frame = PressureVelocityFrame {
                pressure_pa: row[0],
                particle_velocity_m_s: [row[1], row[2], row[3]],
            };
        }
        self.accepted_samples = next;
        Ok(&self.frames)
    }

    /// Pa-equivalent ideal microphone output from the last accepted field.
    /// The front axis points from microphone toward its front source. Repeated
    /// or changed patterns do not advance/reset the shared field history.
    /// This changes orientation only: receiver translation and Doppler are not
    /// inferred. In a relaxing gas, real equilibrium rho*c is a calibration,
    /// not an assertion of unity sensitivity to its dispersive plane waves.
    pub fn microphone_pressure_pa(
        &self,
        receiver: usize,
        pattern: FirstOrder,
    ) -> Result<f64, TransferError> {
        let frame = self
            .current_fields()
            .and_then(|fields| fields.get(receiver))
            .ok_or(TransferError::Invalid(
                "microphone requires an accepted field at an original receiver",
            ))?;
        Ok(pattern
            .observe(
                C64::from_re(frame.pressure_pa),
                frame.particle_velocity_m_s.map(C64::from_re),
                self.bake.medium,
            )?
            .re)
    }
}

/// The existing scalar-plus-Cartesian-vector Fourier owner, with typed p/u
/// conversion at the native frequency seam. Only physical producers call this.
/// It is deliberately not a public arbitrary-kernel relabeling constructor.
pub(super) fn bake_with_errors<E: From<TransferError>>(
    mut solve_at: impl FnMut(f64) -> Result<Vec<C64>, E>,
    declaration: &TransferDeclaration,
    medium: Medium,
    lower_window: LowerBandWindow,
    cancelled: impl FnMut() -> bool,
) -> Result<PressureVelocityBake, E> {
    let cancel = RefCell::new(cancelled);
    if (*cancel.borrow_mut())() {
        return Err(TransferError::Cancelled.into());
    }
    let channels = declaration.admit_work()?;
    lower_window.admit_ladder(declaration)?;
    harmonic::admit_frequency(core::f64::consts::TAU * lower_window.stop_hz(), medium)?;
    let numeric = bake_observation_transfers_with_errors::<E>(
        |omega| {
            if (*cancel.borrow_mut())() {
                return Err(TransferError::Cancelled.into());
            }
            let weight = lower_window.weight(omega / core::f64::consts::TAU);
            if weight == 0.0 {
                // The declared observation window removes this bin. Never
                // query or manufacture a DC/material response outside support.
                return Ok(vec![C64::ZERO; channels]);
            }
            harmonic::admit_frequency(omega, medium)?;
            let mut row = solve_at(omega)?;
            if row.len() != channels || row.iter().any(|&value| !finite(value)) {
                return Err(TransferError::Invalid(
                    "incomplete or nonfinite original pressure/gradient field",
                )
                .into());
            }
            // Exactly the original harmonic momentum convention. Density is
            // real even when the material supplies a dispersive complex k.
            let euler = C64::new(0.0, -1.0 / (omega * medium.density));
            for (pair, field) in row.chunks_exact_mut(4).enumerate() {
                if pair.is_multiple_of(64) && (*cancel.borrow_mut())() {
                    return Err(TransferError::Cancelled.into());
                }
                field[0] = field[0].scale(weight);
                for component in &mut field[1..] {
                    *component = (*component * euler).scale(weight);
                }
            }
            Ok(row)
        },
        declaration,
        lower_window.pass_hz(),
        || (*cancel.borrow_mut())(),
    )?;
    let mut hash = DomainHasher::new("fs-bem:pressure-velocity-observation:v1");
    hash.update(numeric.declaration_hash.as_bytes());
    for value in [
        medium.density,
        medium.sound_speed,
        lower_window.stop_hz(),
        lower_window.pass_hz(),
    ] {
        hash.update(&value.to_bits().to_le_bytes());
    }
    Ok(PressureVelocityBake {
        numeric,
        medium,
        lower_window,
        producer_hash: declaration.producer_hash,
        identity: hash.finalize(),
    })
}
