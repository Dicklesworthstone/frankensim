//! Complete endpoint and clock derivatives of coherent affine tone banks.

use super::{
    MAX_HARMONIC_TONES, ModalAcousticState, ModalAcousticTimeError, ModalToneSensitivityFrame,
    ModalToneSourceSensitivity, PreparedModalToneSensitivity, Sum, invalid, poll,
};

/// One original source's interval-start-basis envelope and clock cotangents.
/// Both envelope endpoints share the SAME carrier basis, including the end.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ModalAffineToneSourcePullback {
    pub cosine_start: f64,
    pub sine_start: f64,
    pub cosine_end: f64,
    pub sine_end: f64,
    /// Source Omega at fixed envelope endpoints and initial phase.
    pub forcing_angular_frequency: f64,
    /// Common phase applied to BOTH endpoint pairs at the interval start.
    pub initial_phase: f64,
}

/// Complete one-mode motion pullback under rotating linear envelopes.
#[derive(Clone, Debug, PartialEq)]
pub struct ModalAffineTonePullback {
    pub initial_displacement: f64,
    pub initial_velocity: f64,
    pub stiffness: f64,
    pub damping_rate: f64,
    sources: [ModalAffineToneSourcePullback; MAX_HARMONIC_TONES],
    count: usize,
}

impl ModalAffineTonePullback {
    /// Original caller order, retaining independently controlled duplicates.
    pub fn source_cotangents(&self) -> &[ModalAffineToneSourcePullback] {
        &self.sources[..self.count]
    }
}

impl PreparedModalToneSensitivity<'_> {
    /// Original affine endpoint and complete local k/c, source-frequency and
    /// common-phase partials. The original coalesced endpoint rows supply the
    /// shared mechanical response; each original source keeps its own clock.
    ///
    /// `u=t/h`: F=(1-u)*(C0*cos+S0*sin)+u*(C1*cos+S1*sin). Therefore the
    /// frequency derivative contains u*(1-u) and u² moments, not just a linear
    /// envelope. All queries are fixed-duration/fixed-basis, allocate nothing,
    /// evaluate no transcendentals and publish no partial result on refusal.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_mode_affine(
        &self,
        mode_index: usize,
        state: ModalAcousticState,
        cosine_start: &[f64],
        sine_start: &[f64],
        cosine_end: &[f64],
        sine_end: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ModalToneSensitivityFrame, ModalAcousticTimeError> {
        poll(&mut cancelled)?;
        let extra = self
            .affine
            .as_ref()
            .ok_or_else(|| invalid("affine sensitivities require affine tone preparation"))?;
        let mode = self
            .bank
            .modes
            .get(mode_index)
            .ok_or_else(|| invalid("tone sensitivity mode index must be admitted"))?;
        let [c, s, dc, ds] =
            self.bank
                .collect_affine(cosine_start, sine_start, cosine_end, sine_end)?;
        let k = self.bank.unique.len();
        let span = mode_index * k..(mode_index + 1) * k;
        let coefficients = &self.coefficients[span.clone()];
        let extra = &extra[span];
        let next = mode.advance(state, &c[..k], &s[..k], Some((&dc[..k], &ds[..k])))?;
        let mut forcing = [Sum::default(); 5];
        for (i, (scalar, affine)) in coefficients.iter().zip(extra).enumerate() {
            poll(&mut cancelled)?;
            for (sum, value) in forcing.iter_mut().zip(scalar.parameter_forcing(c[i], s[i])) {
                sum.add(value);
            }
            for (sum, value) in
                forcing
                    .iter_mut()
                    .zip(affine.parameter_forcing(self.bank.duration_s, dc[i], ds[i]))
            {
                sum.add(value);
            }
        }
        let parameters = coefficients[0].mechanical_partials(state, forcing.map(Sum::value))?;
        let mut sources = [ModalToneSourceSensitivity::default(); MAX_HARMONIC_TONES];
        for (i, &unique) in self.bank.source_to_unique.iter().enumerate() {
            poll(&mut cancelled)?;
            let [qw, vw, qp, vp] = extra[unique].clock_partials(
                &coefficients[unique],
                cosine_start[i],
                sine_start[i],
                cosine_end[i],
                sine_end[i],
            )?;
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

    /// Pull back one coherent affine interval, including both envelope knots.
    /// A trajectory sums adjacent endpoint contributions at each authored knot
    /// and pulls the carrier basis through PreparedModalTones::pullback_force_phase.
    /// When every interval's envelopes are instead authored in an absolute-clock
    /// basis, its initial-phase contribution also needs the elapsed-time chain.
    /// Frequency, phase and mechanical parameters are independent here; callers
    /// apply their authored parameter sharing explicitly across sources/modes.
    #[allow(clippy::too_many_arguments)]
    pub fn pullback_mode_affine(
        &self,
        mode_index: usize,
        state: ModalAcousticState,
        cosine_start: &[f64],
        sine_start: &[f64],
        cosine_end: &[f64],
        sine_end: &[f64],
        endpoint_displacement_cotangent: f64,
        endpoint_velocity_cotangent: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ModalAffineTonePullback, ModalAcousticTimeError> {
        poll(&mut cancelled)?;
        let (qb, vb) = (endpoint_displacement_cotangent, endpoint_velocity_cotangent);
        if ![qb, vb].into_iter().all(f64::is_finite) {
            return Err(invalid("tone endpoint cotangents must be finite"));
        }
        let frame = self.advance_mode_affine(
            mode_index,
            state,
            cosine_start,
            sine_start,
            cosine_end,
            sine_end,
            &mut cancelled,
        )?;
        let k = self.bank.unique.len();
        let span = mode_index * k..(mode_index + 1) * k;
        let coefficients = &self.coefficients[span.clone()];
        let extra = &self
            .affine
            .as_ref()
            .expect("affine advance admitted the preparation")[span];
        let state_input = coefficients[0].input_cotangents(qb, vb)?;
        let mut inputs = [[0.0; 4]; MAX_HARMONIC_TONES];
        for (i, (scalar, affine)) in coefficients.iter().zip(extra).enumerate() {
            poll(&mut cancelled)?;
            inputs[i] = affine.input_cotangents(scalar, qb, vb)?;
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
        let mut sources = [ModalAffineToneSourcePullback::default(); MAX_HARMONIC_TONES];
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
                return Err(invalid("affine source cotangents must be representable"));
            }
            let [cosine_start, sine_start, cosine_end, sine_end] = inputs[unique];
            sources[i] = ModalAffineToneSourcePullback {
                cosine_start,
                sine_start,
                cosine_end,
                sine_end,
                forcing_angular_frequency,
                initial_phase,
            };
        }
        poll(&mut cancelled)?;
        Ok(ModalAffineTonePullback {
            initial_displacement: state_input[0],
            initial_velocity: state_input[1],
            stiffness,
            damping_rate,
            sources,
            count: self.bank.tone_count(),
        })
    }
}
