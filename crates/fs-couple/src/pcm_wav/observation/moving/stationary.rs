//! Stationary pressure and particle motion from the original accepted histories.
//!
//! The force impulse is the analytic integral of the SAME cubic Hermite force
//! used by the pressure renderer, with an explicit initial impulse per source.
//! It is prepared once, never reset at callback boundaries or re-integrated on
//! the microphone clock. Volume flow and force fields add coherently before
//! intensity and the existing ideal first-order microphone are evaluated.
//!
//! This is the stationary, homogeneous lossless compact-source model, not a
//! moving-source velocity approximation. Initial field history, compactness,
//! source bandwidth and linear acoustic validity remain producer obligations.
//! Each source must fit the underlying field owner's representable range,
//! including its individual instantaneous flux, before coherent composition.

use super::*;
use fs_aeroac::curle3d::CompactForceObserver;
use fs_bem::{helmholtz::Medium, near_field::FirstOrder};
use fs_math::c64::C64;

/// Co-located instantaneous fields at one physical receiver time.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct StationaryCompactSample {
    /// Original arrival time [s], without filter delay compensation.
    pub receiver_time_s: f64,
    /// Coherent acoustic pressure perturbation [Pa].
    pub pressure_pa: f64,
    /// Coherent Cartesian particle velocity [m/s].
    pub particle_velocity_m_s: [f64; 3],
    /// Signed instantaneous intensity `p u` [W/m²], not a sum of source powers.
    pub instantaneous_intensity_w_m2: [f64; 3],
    /// Pa-equivalent output of the declared ideal first-order microphone.
    pub microphone_pressure_pa: f64,
}

struct Transfer {
    force: CompactForceObserver,
    direction: [f64; 3],
    monopole_pressure: f64,
    monopole_velocity: [f64; 2],
    first_impulse: usize,
}

/// Finite accepted stationary sources on the existing physical receiver clock.
///
/// Preparation reuses `MovingCompactRenderer` admission and its original cubic
/// coefficients. Callbacks allocate nothing; cancellation or refusal publishes
/// neither a partial field nor an advanced clock. The same cursor serves full
/// field blocks and the `PressureRenderer`/decimator/PCM streaming interfaces.
/// Full fields remain unfiltered physical observations; only microphone pressure
/// passes through a downstream decimator when the caller explicitly adds one.
pub struct StationaryCompactRenderer<'a> {
    original: MovingCompactRenderer<'a>,
    transfers: Vec<Transfer>,
    impulses: Vec<[Dd; 3]>,
    scratch: Vec<StationaryCompactSample>,
    microphone: FirstOrder,
}

impl<'a> StationaryCompactRenderer<'a> {
    /// Bind one fixed receiver to zero-velocity source histories.
    ///
    /// `initial_force_impulses_n_s[j]` is the source's force integral at its
    /// FIRST retained endpoint, not at the first requested output. Zero means
    /// zero initial source-induced impulse field; it is never inferred. A
    /// different integration constant changes near-field particle motion while
    /// leaving pressure unchanged. Missing prehistory is not reconstructed.
    ///
    /// The original budget additionally covers cached transfers, all impulse
    /// endpoints and field scratch. Nonstationary sources refuse before any
    /// model is constructed; use the moving pressure owner for their pressure.
    pub fn new(
        observer: MovingCompactObserver,
        sources: &[CompactSourceHistory<'a>],
        initial_force_impulses_n_s: &[[f64; 3]],
        microphone: FirstOrder,
        config: MovingPressureConfig,
        mut budget: MovingPressureBudget,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, RenderError> {
        check_cancel(cancelled)?;
        if sources.is_empty()
            || sources.len() > budget.maximum_sources
            || sources.len() != initial_force_impulses_n_s.len()
        {
            return Err(sizing(
                "stationary fields require one initial force impulse per admitted source",
            ));
        }
        let mut endpoints = 0_usize;
        for (source, impulse) in sources.iter().zip(initial_force_impulses_n_s) {
            check_cancel(cancelled)?;
            if source.motion.velocity_m_s() != [0.0; 3] || impulse.iter().any(|x| !x.is_finite()) {
                return Err(sizing(
                    "stationary fields require fixed sources and finite initial force impulses",
                ));
            }
            endpoints = endpoints
                .checked_add(source.samples.len())
                .ok_or_else(|| sizing("stationary impulse endpoint count overflow"))?;
        }
        if endpoints > budget.maximum_source_samples {
            return Err(sizing(
                "stationary force histories exceed the endpoint budget",
            ));
        }
        let extra = sources
            .len()
            .checked_mul(size_of::<Transfer>())
            .and_then(|n| {
                endpoints
                    .checked_mul(size_of::<[Dd; 3]>())
                    .and_then(|m| n.checked_add(m))
            })
            .and_then(|n| {
                config
                    .max_block
                    .checked_mul(size_of::<StationaryCompactSample>())
                    .and_then(|m| n.checked_add(m))
            })
            .ok_or_else(|| sizing("stationary field payload overflow"))?;
        budget.maximum_payload_bytes = budget
            .maximum_payload_bytes
            .checked_sub(extra)
            .ok_or_else(|| sizing("stationary fields exceed the retained payload budget"))?;
        let original = MovingCompactRenderer::new(observer, sources, config, budget, cancelled)?;
        let mut transfers = Vec::new();
        let mut impulses = Vec::new();
        let mut scratch = Vec::new();
        transfers
            .try_reserve_exact(sources.len())
            .map_err(|_| sizing("cannot reserve stationary field transfers"))?;
        impulses
            .try_reserve_exact(endpoints)
            .map_err(|_| sizing("cannot reserve stationary force impulses"))?;
        scratch
            .try_reserve_exact(config.max_block)
            .map_err(|_| sizing("cannot reserve stationary field scratch"))?;
        scratch.resize(config.max_block, StationaryCompactSample::default());
        let medium = observer.medium();
        for (track, initial) in original.tracks.iter().zip(initial_force_impulses_n_s) {
            check_cancel(cancelled)?;
            let source = track.source.motion.reference_position_m();
            let delta = core::array::from_fn(|i| observer.observer_position_m()[i] - source[i]);
            let radius = norm(delta);
            let inverse = (0.25 / core::f64::consts::PI) / radius;
            let monopole_pressure = inverse * medium.density_kg_m3();
            let monopole_velocity = [inverse / radius, inverse / medium.sound_speed_m_s()];
            if [
                monopole_pressure,
                monopole_velocity[0],
                monopole_velocity[1],
            ]
            .iter()
            .any(|x| !x.is_finite() || *x <= 0.0)
            {
                return Err(sizing(
                    "stationary monopole coefficients are unrepresentable",
                ));
            }
            let force = CompactForceObserver::new(medium, observer.observer_position_m(), source)
                .map_err(acoustic)?;
            // Fail during preparation even when the first block's source is zero.
            force
                .field_at_retarded([0.0; 3], [0.0; 3], [0.0; 3])
                .map_err(acoustic)?;
            transfers.push(Transfer {
                force,
                direction: delta.map(|x| x / radius),
                monopole_pressure,
                monopole_velocity,
                first_impulse: impulses.len(),
            });
            let mut impulse = initial.map(Dd::from_f64);
            impulses.push(impulse);
            for left in 0..track.source.samples.len() - 1 {
                if left % 64 == 0 {
                    check_cancel(cancelled)?;
                }
                let increment = integrated_force(track, &original.intervals, left, 1.0);
                for axis in 0..3 {
                    impulse[axis] = impulse[axis] + increment[axis];
                }
                if impulse.iter().any(|x| !x.to_f64().is_finite()) {
                    return Err(sizing("accepted force impulse is unrepresentable"));
                }
                impulses.push(impulse);
            }
        }
        check_cancel(cancelled)?;
        Ok(Self {
            original,
            transfers,
            impulses,
            scratch,
            microphone,
        })
    }

    /// Original receiver clock and finite observation window.
    pub const fn config(&self) -> MovingPressureConfig {
        self.original.config()
    }
    /// Original ideal microphone, including its physical front-axis convention.
    pub const fn microphone(&self) -> FirstOrder {
        self.microphone
    }
    /// Number of complete published receiver samples.
    pub const fn samples_rendered(&self) -> u64 {
        self.original.samples_rendered()
    }
    /// Number of unobserved samples; no missing tail is synthesized.
    pub fn remaining_samples(&self) -> u64 {
        self.original.remaining_samples()
    }
    /// Largest complete transactional callback.
    pub const fn max_block_len(&self) -> usize {
        self.original.max_block_len()
    }
    /// Same integer-index clock expression as the original pressure owner.
    pub fn receiver_time_at(&self, sample: u64) -> Result<f64, RenderError> {
        self.original.receiver_time_at(sample)
    }

    /// Publish pressure, velocity, intensity and microphone output together.
    pub fn field_block(
        &mut self,
        output: &mut [StationaryCompactSample],
    ) -> Result<(), RenderError> {
        self.field_block_controlled(output, &mut || false)
    }
    /// Atomic full-field callback with bounded cancellation polling.
    pub fn field_block_controlled(
        &mut self,
        output: &mut [StationaryCompactSample],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), RenderError> {
        self.prepare_block(output.len(), cancelled)?;
        output.copy_from_slice(&self.scratch[..output.len()]);
        self.original.completed += output.len() as u64;
        Ok(())
    }
    /// Emit Pa-equivalent microphone pressure through the existing audio path.
    pub fn block(&mut self, output: &mut [f64]) -> Result<(), RenderError> {
        self.block_controlled(output, &mut || false)
    }
    /// Atomic microphone callback; no normalization, limiter or invented gain.
    pub fn block_controlled(
        &mut self,
        output: &mut [f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), RenderError> {
        self.prepare_block(output.len(), cancelled)?;
        for (out, field) in output.iter_mut().zip(&self.scratch) {
            *out = field.microphone_pressure_pa;
        }
        self.original.completed += output.len() as u64;
        Ok(())
    }

    fn prepare_block(
        &mut self,
        count: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), RenderError> {
        if count == 0 {
            return Err(RenderError::EmptyBlock);
        }
        if count > self.max_block_len() {
            return Err(sizing("stationary field block exceeds its capacity"));
        }
        self.validate_sample_count(
            u64::try_from(count).map_err(|_| sizing("stationary field count exceeds u64"))?,
        )?;
        for i in 0..count {
            check_cancel(cancelled)?;
            let time = self.receiver_time_at(self.original.completed + i as u64)?;
            self.scratch[i] = self.observe(time, cancelled)?;
        }
        check_cancel(cancelled)
    }

    fn observe(
        &self,
        time: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<StationaryCompactSample, RenderError> {
        let mut pressure = Dd::ZERO;
        let mut velocity = [Dd::ZERO; 3];
        for (index, (track, transfer)) in
            self.original.tracks.iter().zip(&self.transfers).enumerate()
        {
            if index % 64 == 0 {
                check_cancel(cancelled)?;
            }
            let prepared = prepare_history(
                &self.original.observer,
                &track.source,
                time,
                self.config().maximum_time_residual_s,
                cancelled,
            )?;
            let source_time = prepared.source_time_s();
            let strength = interpolate(track, &self.original.intervals, source_time)?;
            let samples = track.source.samples;
            let left = samples.partition_point(|sample| sample.source_time_s <= source_time) - 1;
            let mut impulse = self.impulses[transfer.first_impulse + left];
            if source_time != samples[left].source_time_s {
                let fraction = (source_time - samples[left].source_time_s)
                    / (samples[left + 1].source_time_s - samples[left].source_time_s);
                let increment = integrated_force(track, &self.original.intervals, left, fraction);
                for axis in 0..3 {
                    impulse[axis] = impulse[axis] + increment[axis];
                }
            }
            let field = transfer
                .force
                .field_at_retarded(
                    impulse.map(Dd::to_f64),
                    strength.force_on_fluid_n,
                    strength.force_rate_on_fluid_n_s,
                )
                .map_err(acoustic)?;
            let radial_velocity = Dd::from_f64(transfer.monopole_velocity[0])
                * Dd::from_f64(strength.volume_flow_m3_s)
                + Dd::from_f64(transfer.monopole_velocity[1])
                    * Dd::from_f64(strength.volume_flow_rate_m3_s2);
            pressure = pressure
                + Dd::from_f64(field.pressure_pa)
                + Dd::from_f64(transfer.monopole_pressure)
                    * Dd::from_f64(strength.volume_flow_rate_m3_s2);
            for axis in 0..3 {
                velocity[axis] = velocity[axis]
                    + Dd::from_f64(field.particle_velocity_m_s[axis])
                    + Dd::from_f64(transfer.direction[axis]) * radial_velocity;
            }
        }
        let pressure_pa = pressure.to_f64();
        let particle_velocity_m_s = velocity.map(Dd::to_f64);
        let instantaneous_intensity_w_m2 = particle_velocity_m_s.map(|v| pressure_pa * v);
        if !pressure_pa.is_finite()
            || particle_velocity_m_s
                .iter()
                .chain(&instantaneous_intensity_w_m2)
                .any(|x| !x.is_finite())
        {
            return Err(sizing(
                "coherent stationary pressure, velocity or intensity is unrepresentable",
            ));
        }
        let medium = self.original.observer.medium();
        let microphone_pressure_pa = self
            .microphone
            .observe(
                C64::new(pressure_pa, 0.0),
                particle_velocity_m_s.map(|v| C64::new(v, 0.0)),
                Medium {
                    density: medium.density_kg_m3(),
                    sound_speed: medium.sound_speed_m_s(),
                },
            )
            .map_err(|_| sizing("stationary first-order microphone output is unrepresentable"))?
            .re;
        Ok(StationaryCompactSample {
            receiver_time_s: time,
            pressure_pa,
            particle_velocity_m_s,
            instantaneous_intensity_w_m2,
            microphone_pressure_pa,
        })
    }
}

impl PressureRenderer for StationaryCompactRenderer<'_> {
    fn samples_rendered(&self) -> u64 {
        Self::samples_rendered(self)
    }
    fn max_block_len(&self) -> usize {
        Self::max_block_len(self)
    }
    fn validate_sample_rate(&self, rate: u32) -> Result<(), RenderError> {
        self.original.validate_sample_rate(rate)
    }
    fn validate_sample_count(&self, samples: u64) -> Result<(), RenderError> {
        self.original.validate_sample_count(samples)
    }
    fn block(&mut self, output: &mut [f64]) -> Result<(), RenderError> {
        Self::block(self, output)
    }
}

// Integral from the left endpoint to the normalized fraction of one original
// cubic. The DD prefix and identical coefficients preserve history across both
// nonuniform source knots and arbitrary receiver callback partitions.
fn integrated_force(
    track: &PreparedTrack<'_, UniformSourceMotion>,
    intervals: &[CubicInterval],
    left: usize,
    fraction: f64,
) -> [Dd; 3] {
    let samples = track.source.samples;
    let h = Dd::from_f64(samples[left + 1].source_time_s - samples[left].source_time_s);
    let u = Dd::from_f64(fraction);
    core::array::from_fn(|axis| {
        let [b, c, d] = intervals[track.first_interval + left].terms[axis + 1].map(Dd::from_f64);
        h * u
            * (Dd::from_f64(samples[left].strength.force_on_fluid_n[axis])
                + u * (b / Dd::from_f64(2.0)
                    + u * (c / Dd::from_f64(3.0) + u * d / Dd::from_f64(4.0))))
    })
}
