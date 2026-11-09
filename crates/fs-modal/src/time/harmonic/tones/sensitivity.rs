//! Coherent tone-bank derivatives using the original scalar numerical owners.

mod affine;
pub use affine::{ModalAffineTonePullback, ModalAffineToneSourcePullback};

use super::super::sensitivity::AffineForceSensitivity;
use super::{MAX_HARMONIC_TONES, PreparedModalTones, Sum, invalid};
use crate::time::{
    HarmonicForcePhasePullback, LinearForceParameterSensitivity, ModalAcousticState,
    ModalAcousticTimeError, PreparedHarmonicForceSensitivity,
};

/// Additional numerical resources for optional tone-bank sensitivities.
/// The borrowed original bank and caller-owned histories are excluded.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ModalToneSensitivityBudget {
    /// Cache payload plus bounded numerical preparation/query scratch [bytes].
    /// Allocator bookkeeping and compiler stack layout are not measured.
    pub maximum_storage_bytes: usize,
    /// Sum of `8192+512*d` per mode/distinct tone, where d is the original
    /// scalar preparation's number of interval doublings. Affine envelopes
    /// instead use `24576+2048*d`, including quadratic clock moments.
    pub maximum_preparation_work: usize,
    /// One mode's forward or reverse query score: `512*(K+T)+512` for K
    /// distinct and T original tones. Affine envelopes use `1024*(K+T)+1024`.
    /// These scores are not timing guarantees.
    pub maximum_query_work: usize,
}

/// One ORIGINAL tone's independent source-frequency and initial-phase partials.
/// Identical nominal frequencies do not identify their physical parameters.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ModalToneSourceSensitivity {
    pub displacement_by_forcing_angular_frequency: f64,
    pub velocity_by_forcing_angular_frequency: f64,
    pub displacement_by_initial_phase: f64,
    pub velocity_by_initial_phase: f64,
}

/// Original coherent endpoint and local physical derivatives for one mode.
#[derive(Clone, Debug, PartialEq)]
pub struct ModalToneSensitivityFrame {
    /// Exactly the original bank endpoint under its prepared envelope law.
    pub state: ModalAcousticState,
    /// Independent k=omega^2 and physical c partials, with one initial state.
    pub parameters: LinearForceParameterSensitivity,
    sources: [ModalToneSourceSensitivity; MAX_HARMONIC_TONES],
    count: usize,
}
impl ModalToneSensitivityFrame {
    /// Complete original caller order, including cancelled duplicate sources.
    pub fn source_partials(&self) -> &[ModalToneSourceSensitivity] {
        &self.sources[..self.count]
    }
}

/// Scalar-objective derivatives for one original input tone.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ModalToneSourcePullback {
    pub cosine_force: f64,
    pub sine_force: f64,
    pub forcing_angular_frequency: f64,
    /// Common additional phase of this source at the interval start [rad].
    pub initial_phase: f64,
}

/// One mode's original state, mechanics and caller-ordered source cotangents.
#[derive(Clone, Debug, PartialEq)]
pub struct ModalTonePullback {
    pub initial_displacement: f64,
    pub initial_velocity: f64,
    pub stiffness: f64,
    pub damping_rate: f64,
    sources: [ModalToneSourcePullback; MAX_HARMONIC_TONES],
    count: usize,
}
impl ModalTonePullback {
    /// Force/clock/phase cotangents in the original caller order.
    pub fn source_cotangents(&self) -> &[ModalToneSourcePullback] {
        &self.sources[..self.count]
    }
}

/// Optional scalar sensitivity caches paired with their unchanged tone bank.
///
/// Shared k/c derivatives contract the aggregate held-force imbalance and
/// source derivative gains ONCE. Exact duplicate frequencies share caches and
/// coherent primal forces; their independent clock/phase parameters retain
/// the original individual quadratures, even at complete destructive cancellation.
/// Warm queries use fixed-capacity stack results, allocate nothing, and evaluate
/// no transcendental functions. Cancellation publishes no partial result.
///
/// Sum source-frequency cotangents for a shared clock across tones/modes;
/// natural-frequency changes instead use the k/c chain rule. A trajectory also
/// pulls back each source rotation through PreparedModalTones::pullback_force_phase.
/// These are local fixed-duration, fixed-basis motion derivatives. They do not
/// differentiate work/heat, model admission, or nonlinear observation operators.
#[derive(Debug)]
pub struct PreparedModalToneSensitivity<'a> {
    bank: &'a PreparedModalTones,
    coefficients: Box<[PreparedHarmonicForceSensitivity]>,
    affine: Option<Box<[AffineForceSensitivity]>>,
}

impl PreparedModalTones {
    /// Prepare all mode/distinct-tone derivatives under additional explicit
    /// resources. No bank, model or physical state is changed on refusal.
    ///
    /// Affine banks additionally prepare linear-envelope coefficient tangents
    /// and exact quadratic-envelope clock moments. Their queries require both
    /// endpoints through advance_mode_affine/pullback_mode_affine.
    /// Cancellation is polled before allocation and between scalar preparations.
    pub fn prepare_sensitivities(
        &self,
        budget: ModalToneSensitivityBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedModalToneSensitivity<'_>, ModalAcousticTimeError> {
        poll(&mut cancelled)?;
        let count = self
            .mode_count()
            .checked_mul(self.unique.len())
            .ok_or_else(|| invalid("tone sensitivity cache extent overflow"))?;
        let scalar_size = core::mem::size_of::<PreparedHarmonicForceSensitivity>();
        let affine_size = if self.affine {
            core::mem::size_of::<AffineForceSensitivity>()
        } else {
            0
        };
        let mut scratch = core::mem::size_of::<ModalToneSensitivityFrame>()
            + core::mem::size_of::<ModalTonePullback>()
            + core::mem::size_of::<[[f64; 4]; MAX_HARMONIC_TONES]>()
            + core::mem::size_of::<[[f64; MAX_HARMONIC_TONES]; 2]>()
            + core::mem::size_of::<[Sum; 2 * MAX_HARMONIC_TONES]>()
            + core::mem::size_of::<[Sum; 5]>()
            + scalar_size;
        if self.affine {
            scratch += core::mem::size_of::<ModalAffineTonePullback>()
                + core::mem::size_of::<[[f64; MAX_HARMONIC_TONES]; 2]>()
                + affine_size;
        }
        let storage = (count as u128) * (scalar_size + affine_size) as u128
            + scratch as u128
            + core::mem::size_of::<PreparedModalToneSensitivity<'_>>() as u128
            + core::mem::size_of::<Vec<PreparedHarmonicForceSensitivity>>() as u128
            + core::mem::size_of::<Vec<AffineForceSensitivity>>() as u128;
        cap(
            "tone sensitivity storage bytes",
            storage,
            budget.maximum_storage_bytes,
        )?;
        cap(
            "tone sensitivity query work score",
            ((if self.affine { 1024 } else { 512 }) * (self.unique.len() + self.tone_count() + 1))
                as u128,
            budget.maximum_query_work,
        )?;
        let mut work = 0u128;
        for mode in &self.modes {
            for motion in &mode.motion {
                poll(&mut cancelled)?;
                let mut scale = motion
                    .transition
                    .angle
                    .max(motion.transition.drag)
                    .max(motion.forcing_phase_rad);
                let mut doublings = 0;
                while scale > 0.25 {
                    scale *= 0.5;
                    doublings += 1;
                }
                work += if self.affine {
                    24576 + 2048 * doublings
                } else {
                    8192 + 512 * doublings
                };
                cap(
                    "tone sensitivity preparation work score",
                    work,
                    budget.maximum_preparation_work,
                )?;
            }
        }
        poll(&mut cancelled)?;
        let mut coefficients = Vec::with_capacity(count);
        let mut affine = self.affine.then(|| Vec::with_capacity(count));
        for mode in &self.modes {
            for motion in &mode.motion {
                poll(&mut cancelled)?;
                coefficients.push(motion.prepare_sensitivities()?);
                if let Some(extra) = &mut affine {
                    poll(&mut cancelled)?;
                    extra.push(AffineForceSensitivity::prepare(motion)?);
                }
            }
        }
        poll(&mut cancelled)?;
        Ok(PreparedModalToneSensitivity {
            bank: self,
            coefficients: coefficients.into_boxed_slice(),
            affine: affine.map(Vec::into_boxed_slice),
        })
    }

    /// Source-state reverse chain for one ORIGINAL tone. Add these force
    /// cotangents to the mechanical endpoint's force cotangents, and add both
    /// frequency contributions. The original phase rotation is unchanged.
    pub fn pullback_force_phase(
        &self,
        tone_index: usize,
        cosine_force: f64,
        sine_force: f64,
        endpoint_cosine_cotangent: f64,
        endpoint_sine_cotangent: f64,
    ) -> Result<HarmonicForcePhasePullback, ModalAcousticTimeError> {
        let &unique = self
            .source_to_unique
            .get(tone_index)
            .ok_or_else(|| invalid("tone index must belong to the preparation"))?;
        self.modes
            .first()
            .ok_or_else(|| invalid("source phase requires a physical mode"))?
            .motion[unique]
            .pullback_force_phase(
                cosine_force,
                sine_force,
                endpoint_cosine_cotangent,
                endpoint_sine_cotangent,
            )
    }
}

impl PreparedModalToneSensitivity<'_> {
    /// Original endpoint and all local k/c and individual source-clock partials.
    /// The caller owns physical state/observer admission, as for advance_mode.
    pub fn advance_mode(
        &self,
        mode_index: usize,
        state: ModalAcousticState,
        cosine: &[f64],
        sine: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ModalToneSensitivityFrame, ModalAcousticTimeError> {
        poll(&mut cancelled)?;
        if self.bank.affine {
            return Err(invalid(
                "affine sensitivities require both envelope endpoints",
            ));
        }
        let mode = self
            .bank
            .modes
            .get(mode_index)
            .ok_or_else(|| invalid("tone sensitivity mode index must be admitted"))?;
        let (c, s) = self.bank.collect_forces(cosine, sine)?;
        let k = self.bank.unique.len();
        let coefficients = &self.coefficients[mode_index * k..(mode_index + 1) * k];
        let next = mode.advance(state, &c[..k], &s[..k], None)?;
        let mut forcing = [Sum::default(); 5];
        for (i, scalar) in coefficients.iter().enumerate() {
            poll(&mut cancelled)?;
            for (sum, value) in forcing.iter_mut().zip(scalar.parameter_forcing(c[i], s[i])) {
                sum.add(value);
            }
        }
        let parameters = coefficients[0].mechanical_partials(state, forcing.map(Sum::value))?;
        let mut sources = [ModalToneSourceSensitivity::default(); MAX_HARMONIC_TONES];
        for (i, &unique) in self.bank.source_to_unique.iter().enumerate() {
            poll(&mut cancelled)?;
            let [qw, vw, qp, vp] = coefficients[unique].clock_partials(cosine[i], sine[i])?;
            sources[i] = ModalToneSourceSensitivity {
                displacement_by_forcing_angular_frequency: qw,
                velocity_by_forcing_angular_frequency: vw,
                displacement_by_initial_phase: qp,
                velocity_by_initial_phase: vp,
            };
        }
        poll(&mut cancelled)?;
        Ok(ModalToneSensitivityFrame {
            state: next,
            parameters,
            sources,
            count: self.bank.tone_count(),
        })
    }

    /// Pull endpoint cotangents through one complete coherent mode response.
    /// Independent original-tone columns survive nominal duplicate cancellation;
    /// explicitly sum them when the authored source shares frequency or phase.
    #[allow(clippy::too_many_arguments)]
    pub fn pullback_mode(
        &self,
        mode_index: usize,
        state: ModalAcousticState,
        cosine: &[f64],
        sine: &[f64],
        endpoint_displacement_cotangent: f64,
        endpoint_velocity_cotangent: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ModalTonePullback, ModalAcousticTimeError> {
        poll(&mut cancelled)?;
        let (qb, vb) = (endpoint_displacement_cotangent, endpoint_velocity_cotangent);
        if ![qb, vb].into_iter().all(f64::is_finite) {
            return Err(invalid("tone endpoint cotangents must be finite"));
        }
        let frame = self.advance_mode(mode_index, state, cosine, sine, &mut cancelled)?;
        let k = self.bank.unique.len();
        let coefficients = &self.coefficients[mode_index * k..(mode_index + 1) * k];
        let mut inputs = [[0.0; 4]; MAX_HARMONIC_TONES];
        for (i, scalar) in coefficients.iter().enumerate() {
            poll(&mut cancelled)?;
            inputs[i] = scalar.input_cotangents(qb, vb)?;
        }
        let p = frame.parameters;
        let stiffness = qb.mul_add(p.displacement_by_stiffness, vb * p.velocity_by_stiffness);
        let damping_rate = qb.mul_add(
            p.displacement_by_damping_rate,
            vb * p.velocity_by_damping_rate,
        );
        if ![stiffness, damping_rate].into_iter().all(f64::is_finite) {
            return Err(invalid(
                "tone physical parameter cotangents must be representable",
            ));
        }
        let mut sources = [ModalToneSourcePullback::default(); MAX_HARMONIC_TONES];
        for (i, &unique) in self.bank.source_to_unique.iter().enumerate() {
            poll(&mut cancelled)?;
            let p = frame.sources[i];
            let forcing_angular_frequency = qb.mul_add(
                p.displacement_by_forcing_angular_frequency,
                vb * p.velocity_by_forcing_angular_frequency,
            );
            let initial_phase = qb.mul_add(
                p.displacement_by_initial_phase,
                vb * p.velocity_by_initial_phase,
            );
            if ![forcing_angular_frequency, initial_phase]
                .into_iter()
                .all(f64::is_finite)
            {
                return Err(invalid("tone source cotangents must be representable"));
            }
            sources[i] = ModalToneSourcePullback {
                cosine_force: inputs[unique][2],
                sine_force: inputs[unique][3],
                forcing_angular_frequency,
                initial_phase,
            };
        }
        poll(&mut cancelled)?;
        Ok(ModalTonePullback {
            initial_displacement: inputs[0][0],
            initial_velocity: inputs[0][1],
            stiffness,
            damping_rate,
            sources,
            count: self.bank.tone_count(),
        })
    }
}

fn cap(what: &'static str, value: u128, limit: usize) -> Result<(), ModalAcousticTimeError> {
    if value > limit as u128 {
        Err(ModalAcousticTimeError::BudgetExceeded {
            what,
            value: value as f64,
            limit: limit as f64,
        })
    } else {
        Ok(())
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), ModalAcousticTimeError> {
    if cancelled() {
        Err(invalid("tone sensitivity operation cancelled"))
    } else {
        Ok(())
    }
}
