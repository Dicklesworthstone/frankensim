//! Moving volume-source fields on the original accepted source and receiver clocks.
//!
//! The shared history owner compiles Q and Qdot once. Every direct or reflected
//! path reads that SAME interpolant at its own retarded source time, then uses
//! the original moving-source potential gradient for all velocity components.
//! Coherent pressure and velocity are summed before intensity and microphone
//! projection. Neither geometric spreading nor Doppler is applied a second time.

use super::*;
use fs_bem::{helmholtz::Medium, near_field::FirstOrder};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::c64::C64;

/// Co-located instantaneous acoustic fields on one physical receiver clock.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct CompactFieldSample {
    /// Original arrival time [s], without filter delay compensation.
    pub receiver_time_s: f64,
    /// Coherent acoustic pressure perturbation [Pa].
    pub pressure_pa: f64,
    /// Coherent Cartesian perturbation velocity [m/s] in the medium's axes.
    /// This is fluid particle motion, not source or microphone translation.
    pub particle_velocity_m_s: [f64; 3],
    /// Signed instantaneous pressure-work flux `p u` [W/m²].
    /// In uniform flow this excludes mean-flow advection of acoustic energy.
    pub instantaneous_intensity_w_m2: [f64; 3],
    /// Pa-equivalent output of the declared ideal first-order microphone.
    pub microphone_pressure_pa: f64,
}

impl CompactFieldSample {
    // Stationary force/volume and moving volume observations use the same
    // physical reduction. Preserve the stationary owner's rounding order.
    pub(super) fn observe(
        time: f64,
        pressure: Dd,
        velocity: [Dd; 3],
        microphone: FirstOrder,
        medium: FluidAcousticState,
    ) -> Result<Self, RenderError> {
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
                "coherent compact pressure, velocity or intensity is unrepresentable",
            ));
        }
        let microphone_pressure_pa = microphone
            .observe(
                C64::new(pressure_pa, 0.0),
                particle_velocity_m_s.map(|v| C64::new(v, 0.0)),
                Medium {
                    density: medium.density_kg_m3(),
                    sound_speed: medium.sound_speed_m_s(),
                },
            )
            .map_err(|_| sizing("compact first-order microphone output is unrepresentable"))?
            .re;
        Ok(Self {
            receiver_time_s: time,
            pressure_pa,
            particle_velocity_m_s,
            instantaneous_intensity_w_m2,
            microphone_pressure_pa,
        })
    }
}

/// Coherent full fields and directional audio from accepted moving Q histories.
///
/// Uniform and curved sources, moving receivers, ideal-wall paths and uniform
/// flow use the SAME motion policies and history admission as
/// [`MovingCompactRenderer`]. All source forces AND force derivatives must be
/// zero at every accepted endpoint: moving force-source velocity requires past
/// worldline information beyond the local retarded strength. Such histories
/// refuse at construction, rather than losing their dipole contribution.
///
/// Near-field, nonradial and source-acceleration terms come from
/// [`PreparedMovingCompactSample::evaluate_monopole_field`]. The ideal pickup
/// uses a fixed front axis in these Cartesian coordinates, including for a
/// moving microphone. It describes an ideal field projection, not the capsule's
/// mechanical transfer function or motion-induced aerodynamic self-noise.
///
/// Callbacks allocate nothing. Fields and microphone pressure share one sample
/// cursor; a refused or cancelled block publishes neither output nor clock.
/// Decimation/PCM use microphone pressure through [`PressureRenderer`]. Full
/// field blocks remain physical observations at the original receiver rate.
/// Source compactness, bandwidth and linearity remain producer obligations;
/// no missing prehistory or tail is synthesized. Each individual path must fit
/// its underlying field owner's representable velocity/flux range.
pub struct MovingMonopoleRenderer<'a, M: CompactMotion = UniformSourceMotion> {
    original: MovingCompactRenderer<'a, M>,
    scratch: Vec<CompactFieldSample>,
    microphone: FirstOrder,
    medium: FluidAcousticState,
}

impl<'a, M: CompactMotion> MovingMonopoleRenderer<'a, M> {
    /// Admit original accepted histories and preallocate full-field callbacks.
    /// The original budget additionally charges all retained field scratch.
    pub fn new(
        observer: M::Observer,
        sources: &[CompactSourceHistory<'a, M>],
        microphone: FirstOrder,
        config: MovingPressureConfig,
        mut budget: MovingPressureBudget,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, RenderError> {
        check_cancel(cancelled)?;
        let extra = config
            .max_block
            .checked_mul(size_of::<CompactFieldSample>())
            .ok_or_else(|| sizing("moving field scratch payload overflows"))?;
        budget.maximum_payload_bytes = budget
            .maximum_payload_bytes
            .checked_sub(extra)
            .ok_or_else(|| sizing("moving fields exceed the retained payload budget"))?;
        let original = MovingCompactRenderer::new(observer, sources, config, budget, cancelled)?;
        // Admission above bounds the endpoint count and validates all values.
        // Zero Hermite endpoint values AND slopes imply zero force throughout.
        for track in &original.tracks {
            for (index, sample) in track.source.samples.iter().enumerate() {
                if index % 64 == 0 {
                    check_cancel(cancelled)?;
                }
                if sample.strength.force_on_fluid_n != [0.0; 3]
                    || sample.strength.force_rate_on_fluid_n_s != [0.0; 3]
                {
                    return Err(sizing(
                        "moving monopole fields require zero force and force-rate histories",
                    ));
                }
            }
        }
        let medium = original.observer_at(config.first_receiver_time_s)?.medium();
        let mut scratch = Vec::new();
        scratch
            .try_reserve_exact(config.max_block)
            .map_err(|_| sizing("cannot reserve moving field scratch"))?;
        scratch.resize(config.max_block, CompactFieldSample::default());
        check_cancel(cancelled)?;
        Ok(Self {
            original,
            scratch,
            microphone,
            medium,
        })
    }

    /// Original receiver owner, including its identity and accepted trajectory.
    pub const fn observer(&self) -> &M::Observer {
        self.original.observer()
    }
    /// Physical receiver position in the medium frame at an admitted time.
    pub fn observer_at(&self, time: f64) -> Result<MovingCompactObserver, RenderError> {
        self.original.observer_at(time)
    }
    /// Original physical receiver clock, bandwidth and finite horizon.
    pub const fn config(&self) -> MovingPressureConfig {
        self.original.config()
    }
    /// Ideal microphone with its fixed physical front-axis convention.
    pub const fn microphone(&self) -> FirstOrder {
        self.microphone
    }
    /// Physical source count; reflected images do not become new sources.
    pub fn source_count(&self) -> usize {
        self.original.source_count()
    }
    /// Complete direct/reflected path count, charged to the original work bound.
    pub const fn source_path_count(&self) -> usize {
        self.original.source_path_count()
    }
    /// Complete published samples on the shared receiver clock.
    pub const fn samples_rendered(&self) -> u64 {
        self.original.samples_rendered()
    }
    /// Unobserved samples within retained physical source coverage.
    pub fn remaining_samples(&self) -> u64 {
        self.original.remaining_samples()
    }
    /// Largest callback fitting the retained scratch buffers.
    pub const fn max_block_len(&self) -> usize {
        self.original.max_block_len()
    }
    /// Original bounded root work over all paths and receiver samples.
    pub const fn maximum_retarded_iterations(&self) -> u64 {
        self.original.maximum_retarded_iterations()
    }
    /// Original integer-index timestamp; block partitioning cannot change it.
    pub fn receiver_time_at(&self, sample: u64) -> Result<f64, RenderError> {
        self.original.receiver_time_at(sample)
    }

    /// Publish coherent pressure, vector velocity, intensity and pickup together.
    pub fn field_block(&mut self, output: &mut [CompactFieldSample]) -> Result<(), RenderError> {
        self.field_block_controlled(output, &mut || false)
    }
    /// Atomic full-field callback with bounded cancellation polling.
    pub fn field_block_controlled(
        &mut self,
        output: &mut [CompactFieldSample],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), RenderError> {
        self.prepare_block(output.len(), cancelled)?;
        output.copy_from_slice(&self.scratch[..output.len()]);
        self.original.completed += output.len() as u64;
        Ok(())
    }
    /// Emit Pa-equivalent ideal microphone pressure through the audio pipeline.
    pub fn block(&mut self, output: &mut [f64]) -> Result<(), RenderError> {
        self.block_controlled(output, &mut || false)
    }
    /// Atomic pickup callback, without material normalization or hidden gain.
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
            return Err(sizing("moving field block exceeds its capacity"));
        }
        self.validate_sample_count(
            u64::try_from(count).map_err(|_| sizing("moving field count exceeds u64"))?,
        )?;
        for i in 0..count {
            check_cancel(cancelled)?;
            let time = self.receiver_time_at(self.original.completed + i as u64)?;
            let mut pressure = Dd::ZERO;
            let mut velocity = [Dd::ZERO; 3];
            for track in &self.original.tracks {
                for path in 0..track.source.motion.path_count() {
                    check_cancel(cancelled)?;
                    let prepared = prepare_path_history(
                        &self.original.observer,
                        &track.source,
                        path,
                        time,
                        self.config().maximum_time_residual_s,
                        cancelled,
                    )?;
                    let strength =
                        interpolate(track, &self.original.intervals, prepared.source_time_s())?;
                    let field = prepared
                        .evaluate_monopole_field(
                            strength.volume_flow_m3_s,
                            strength.volume_flow_rate_m3_s2,
                        )
                        .map_err(acoustic)?;
                    pressure = pressure + Dd::from_f64(field.sample.pressure_pa);
                    for axis in 0..3 {
                        velocity[axis] =
                            velocity[axis] + Dd::from_f64(field.particle_velocity_m_s[axis]);
                    }
                }
            }
            self.scratch[i] = CompactFieldSample::observe(
                time,
                pressure,
                velocity,
                self.microphone,
                self.medium,
            )?;
        }
        check_cancel(cancelled)
    }
}

impl<M: CompactMotion> PressureRenderer for MovingMonopoleRenderer<'_, M> {
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
