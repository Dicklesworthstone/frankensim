//! Passive, integer-delay characteristic propagation with observable storage.
//!
//! This is the pure-delay/resistive-load member of the characteristic runtime,
//! not a fitted-filter energy certificate. Two opposing N-sample delays store
//! pressure waves. With volume-flow impedance Z, E = dt/Z sum(p+^2 + p-^2).
//! At the far end p- = r p+, |r| <= 1; removed wave power is nonnegative.
//! No fractional-delay interpolation, distributed loss, dispersion or radiation
//! model is silently added. See Smith, Physical Audio Signal Processing,
//! Digital Waveguide Theory, acoustic tubes and energy-density waves.

/// Parameters of an initially quiescent bidirectional characteristic line.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WaveguideSpec {
    /// One-way transit time in complete samples (at least one).
    pub one_way_samples: usize,
    /// Pressure/volume-flow impedance [Pa s/m^3].
    pub impedance_pa_s_m3: f64,
    /// Sample interval [s].
    pub time_step_s: f64,
    /// Memoryless terminal pressure reflectance, in [-1, 1].
    pub reflection: f64,
    /// Admission ceiling for requested f64 buffer payloads, not allocator/RSS.
    pub max_memory_bytes: usize,
}

/// Invalid admission or a refused, uncommitted numerical step.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct WaveguideError(pub &'static str);

impl core::fmt::Display for WaveguideError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.write_str(self.0)
    }
}
impl std::error::Error for WaveguideError {}

/// Simultaneous port observations and end-of-step wave storage.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct WaveguideFrame {
    /// Wave arriving at the inlet before this step [Pa].
    pub incoming_pressure_pa: f64,
    /// Inlet pressure and flow use that same incoming wave.
    pub inlet_pressure_pa: f64,
    /// Net flow into the line at the inlet [m^3/s].
    pub inlet_flow_m3_s: f64,
    /// Pressure at the terminal, after one-way propagation [Pa].
    pub terminal_pressure_pa: f64,
    /// Flow into the terminal load [m^3/s].
    pub terminal_flow_m3_s: f64,
    /// Energy in BOTH directions after the shift [J].
    pub stored_energy_j: f64,
    /// Change in actual wave storage [J], not inferred from supplied work.
    pub storage_change_j: f64,
    /// Work entering the inlet over this step [J].
    pub inlet_work_j: f64,
    /// Energy absorbed at the passive terminal over this step [J].
    pub terminal_loss_j: f64,
}
impl WaveguideFrame {
    /// Discrete balance diagnostic [J]; subject to floating-point rounding.
    #[must_use]
    pub fn balance_residual_j(&self) -> f64 {
        self.storage_change_j + self.terminal_loss_j - self.inlet_work_j
    }
}

/// Exact sample shifts with a fixed-order quadratic-storage reduction.
///
/// Stepping allocates nothing. The two energy trees update in O(log N) without
/// accumulating roundoff in a running energy integral. Allocation is admitted
/// before construction; bounded work is not a hard-real-time certification.
#[derive(Debug, Clone)]
pub struct PassiveWaveguide {
    spec: WaveguideSpec,
    waves: Vec<f64>,
    energy: Vec<f64>,
    leaves: usize,
    head: usize,
    sqrt_energy_scale: f64,
}

impl PassiveWaveguide {
    /// Create a zero-state line under an explicit memory-payload budget.
    ///
    /// # Errors
    /// Nonphysical/nonfinite parameters, size overflow, budget or allocation.
    pub fn new(spec: WaveguideSpec) -> Result<Self, WaveguideError> {
        if spec.one_way_samples == 0
            || !spec.impedance_pa_s_m3.is_finite() || spec.impedance_pa_s_m3 <= 0.0
            || !spec.time_step_s.is_finite() || spec.time_step_s <= 0.0
            || !spec.reflection.is_finite() || spec.reflection.abs() > 1.0
        {
            return Err(WaveguideError("positive delay, time and impedance, and passive finite reflectance required"));
        }
        let scale = spec.time_step_s / spec.impedance_pa_s_m3;
        if !scale.is_finite() || scale <= 0.0 {
            return Err(WaveguideError("waveguide energy scale is not representable"));
        }
        let size_error = WaveguideError("waveguide buffer size overflow");
        let leaves = spec.one_way_samples.checked_next_power_of_two().ok_or(size_error)?;
        let wave_len = spec.one_way_samples.checked_mul(2).ok_or(size_error)?;
        let energy_len = leaves.checked_mul(4).ok_or(size_error)?;
        let bytes = wave_len.checked_add(energy_len)
            .and_then(|n| n.checked_mul(core::mem::size_of::<f64>())).ok_or(size_error)?;
        if bytes > spec.max_memory_bytes {
            return Err(WaveguideError("waveguide buffer payload exceeds memory budget"));
        }
        let allocate = |len: usize| -> Result<Vec<f64>, WaveguideError> {
            let mut buffer = Vec::new();
            buffer.try_reserve_exact(len)
                .map_err(|_| WaveguideError("waveguide buffer allocation failed"))?;
            buffer.resize(len, 0.0);
            Ok(buffer)
        };
        Ok(Self {
            spec, waves: allocate(wave_len)?, energy: allocate(energy_len)?,
            leaves, head: 0, sqrt_energy_scale: fs_math::det::sqrt(scale),
        })
    }

    /// Immutable propagation/load declaration.
    #[must_use]
    pub const fn spec(&self) -> &WaveguideSpec { &self.spec }

    /// Incoming pressure to use in the inlet junction's NEXT solve.
    #[must_use]
    pub fn incoming_pressure_pa(&self) -> f64 {
        self.waves[self.spec.one_way_samples + self.head]
    }

    /// Fixed-order sum of the wave-state energies [J].
    #[must_use]
    pub fn stored_energy_j(&self) -> f64 {
        self.energy[1] + self.energy[2 * self.leaves + 1]
    }

    fn wave_energy(&self, pressure: f64) -> f64 {
        let scaled = pressure * self.sqrt_energy_scale;
        scaled * scaled
    }

    // Each direction has a separate tree; preview and commit perform the same
    // additions, in the same order, along the single replacement leaf's path.
    fn replaced_root(&self, offset: usize, value: f64) -> f64 {
        let mut node = self.leaves + self.head;
        let mut sum = value;
        while node > 1 {
            sum = if node % 2 == 0 {
                sum + self.energy[offset + node + 1]
            } else {
                self.energy[offset + node - 1] + sum
            };
            node /= 2;
        }
        sum
    }

    fn replace_leaf(&mut self, offset: usize, value: f64) {
        let mut node = self.leaves + self.head;
        self.energy[offset + node] = value;
        while node > 1 {
            node /= 2;
            self.energy[offset + node] = self.energy[offset + 2 * node]
                + self.energy[offset + 2 * node + 1];
        }
    }

    /// Calculate an entire step without changing either delay or its storage.
    /// Used by coupled solvers to admit all participants before any publication.
    ///
    /// # Errors
    /// Any nonfinite input, wave energy, port observation or work.
    pub fn preview_step(&self, outgoing: f64) -> Result<WaveguideFrame, WaveguideError> {
        let incoming = self.incoming_pressure_pa();
        let incident = self.waves[self.head];
        let reflected = self.spec.reflection * incident;
        let stored = self.replaced_root(0, self.wave_energy(outgoing))
            + self.replaced_root(2 * self.leaves, self.wave_energy(reflected));
        let inlet_pressure = outgoing + incoming;
        let inlet_flow = (outgoing - incoming) / self.spec.impedance_pa_s_m3;
        let frame = WaveguideFrame {
            incoming_pressure_pa: incoming,
            inlet_pressure_pa: inlet_pressure,
            inlet_flow_m3_s: inlet_flow,
            terminal_pressure_pa: incident + reflected,
            terminal_flow_m3_s: (incident - reflected) / self.spec.impedance_pa_s_m3,
            stored_energy_j: stored,
            storage_change_j: stored - self.stored_energy_j(),
            inlet_work_j: inlet_pressure * (inlet_flow * self.spec.time_step_s),
            terminal_loss_j: self.wave_energy(incident)
                * (1.0 - self.spec.reflection) * (1.0 + self.spec.reflection),
        };
        if ![outgoing, reflected, frame.inlet_pressure_pa, frame.inlet_flow_m3_s,
            frame.terminal_pressure_pa, frame.terminal_flow_m3_s, stored,
            frame.storage_change_j, frame.inlet_work_j, frame.terminal_loss_j,
            frame.balance_residual_j()].iter().all(|v| v.is_finite())
        {
            return Err(WaveguideError("waveguide step left the finite set"));
        }
        Ok(frame)
    }

    /// Accept one outgoing wave, shift both directions, and reflect at the end.
    /// The returned incoming wave is the one used BEFORE the shift, not the
    /// next sample's feedback. No output or state is published on a refusal.
    ///
    /// # Errors
    /// Same numerical admission as `preview_step`; no allocation in this call.
    pub fn step(&mut self, outgoing: f64) -> Result<WaveguideFrame, WaveguideError> {
        let frame = self.preview_step(outgoing)?;
        let reflected = self.spec.reflection * self.waves[self.head];
        self.replace_leaf(0, self.wave_energy(outgoing));
        self.replace_leaf(2 * self.leaves, self.wave_energy(reflected));
        self.waves[self.head] = outgoing;
        self.waves[self.spec.one_way_samples + self.head] = reflected;
        self.head = (self.head + 1) % self.spec.one_way_samples;
        Ok(frame)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn spec(samples: usize, reflection: f64) -> WaveguideSpec {
        WaveguideSpec {
            one_way_samples: samples, impedance_pa_s_m3: 1e6,
            time_step_s: 1.0 / 48000.0, reflection, max_memory_bytes: 1 << 20,
        }
    }

    #[test]
    fn pulse_reaches_terminal_then_returns_after_the_full_round_trip() {
        for delay in [1, 3, 32] {
            for r in [-1.0, -0.7, 0.0, 0.4, 1.0] {
                let mut line = PassiveWaveguide::new(spec(delay, r)).unwrap();
                for n in 0..=2 * delay + 1 {
                    let outgoing = if n == 0 { 20.0 } else { 0.0 };
                    let f = line.step(outgoing).unwrap();
                    assert_eq!(f.incoming_pressure_pa, if n == 2 * delay { r * 20.0 } else { 0.0 });
                    assert_eq!(f.terminal_pressure_pa, if n == delay { 20.0 + r * 20.0 } else { 0.0 });
                    assert_eq!(f.inlet_pressure_pa, outgoing + f.incoming_pressure_pa);
                }
                assert_eq!(line.stored_energy_j(), 0.0);
            }
        }
    }

    #[test]
    fn actual_wave_storage_closes_balance_without_a_running_integral() {
        for r in [-1.0, -0.7, 0.0, 0.4, 1.0] {
            let mut line = PassiveWaveguide::new(spec(7, r)).unwrap();
            for n in 0..2000 {
                let input = if n < 1970 { f64::from(n % 31) - 15.0 } else { 0.0 };
                let before = line.stored_energy_j();
                let f = line.step(input).unwrap();
                // Independent direct sum, not the tree or port-work equation.
                let direct = line.waves.iter().map(|p| p * p).sum::<f64>()
                    * line.spec.time_step_s / line.spec.impedance_pa_s_m3;
                let scale = before + f.stored_energy_j + f.inlet_work_j.abs()
                    + f.terminal_loss_j + f64::MIN_POSITIVE;
                assert!((direct - f.stored_energy_j).abs() <= 1e-13 * scale);
                assert!(f.balance_residual_j().abs() <= 1e-13 * scale);
                assert!(f.terminal_loss_j >= 0.0);
                let terminal_work = f.terminal_pressure_pa
                    * f.terminal_flow_m3_s * line.spec.time_step_s;
                assert!((terminal_work - f.terminal_loss_j).abs() <= 1e-13 * scale);
            }
            // Exact zero after complete drain, not accumulated energy drift.
            assert_eq!(line.stored_energy_j(), 0.0);
        }
    }

    #[test]
    fn preview_and_refused_steps_leave_every_state_byte_unchanged() {
        let mut line = PassiveWaveguide::new(spec(3, -0.7)).unwrap();
        line.step(30.0).unwrap();
        let before = line.clone();
        let preview = line.preview_step(7.0).unwrap();
        for bad in [f64::NAN, f64::INFINITY, f64::MAX] {
            assert!(line.step(bad).is_err());
        }
        assert_eq!(line.waves, before.waves);
        assert_eq!(line.energy, before.energy);
        assert_eq!(line.head, before.head);
        assert_eq!(line.step(7.0).unwrap(), preview);
        assert_eq!(line.stored_energy_j(), preview.stored_energy_j);
    }

    #[test]
    fn storage_budget_and_passivity_are_admitted_before_allocation() {
        let mut s = spec(3, 0.5);
        // 6 wave values plus two eight-value energy trees.
        s.max_memory_bytes = 22 * 8;
        assert!(PassiveWaveguide::new(s).is_ok());
        s.max_memory_bytes -= 1;
        assert!(PassiveWaveguide::new(s).is_err());
        for r in [-1.01, 1.01, f64::NAN, f64::INFINITY] {
            assert!(PassiveWaveguide::new(spec(3, r)).is_err());
        }
        for n in [0, usize::MAX] {
            assert!(PassiveWaveguide::new(spec(n, 0.0)).is_err());
        }
    }
}
