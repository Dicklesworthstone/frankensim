//! Moving aperture coupled to an explicit uniform, lossless acoustic tube.
//!
//! The existing midpoint aperture island and fs-vfit characteristic shifts
//! exchange pressure and flow at the same sample. Wave energy, contact energy,
//! terminal absorption and BOTH external supply ports close one accounting
//! window. This is the cylindrical/integer-transit/resistive-termination slice,
//! not a certificate for arbitrary fitted bores or a microphone/radiation model.
//! The requested and represented lengths remain distinct; any transit-time
//! rounding must fit the caller's explicit geometry tolerance.

use super::dynamic::{
    ApertureDrive, ApertureFrame, ApertureProgress, ApertureTerminal, DynamicAperture,
};
use crate::acoustic_realize::AcousticRealizeError;
use fs_exec::CancelGate;
use fs_vfit::waveguide::{PassiveWaveguide, WaveguideFrame, WaveguideSpec};

/// Physical tube and admitted transit-time approximation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct UniformTubeSpec {
    /// Requested centerline length [m].
    pub length_m: f64,
    /// Constant internal radius [m].
    pub radius_m: f64,
    /// Sound speed in the same fluid as the aperture [m/s].
    pub sound_speed_m_s: f64,
    /// Pressure reflectance at the far end, [-1, 1]. This is an explicit
    /// memoryless passive load, not an inferred radiation termination.
    pub terminal_reflection: f64,
    /// Maximum allowed absolute length change from integer sample transit [m].
    pub max_length_error_m: f64,
    /// Payload budget for propagation buffers only (not aperture scratch/RSS).
    pub max_wave_memory_bytes: usize,
}

fn invalid(what: &'static str) -> AcousticRealizeError {
    AcousticRealizeError::InvalidDescription { what }
}

impl UniformTubeSpec {
    /// Compute Z = rho c / (pi r^2) for binding the aperture's load.
    /// Set `DynamicApertureSpec::impedance_pa_s_m3` to this value; construction
    /// refuses mismatches rather than silently changing either participant.
    ///
    /// # Errors
    /// Nonfinite/nonpositive fluid or geometry, or unrepresentable impedance.
    pub fn characteristic_impedance(&self, density_kg_m3: f64) -> Result<f64, AcousticRealizeError> {
        if ![self.length_m, self.radius_m, self.sound_speed_m_s, density_kg_m3]
            .iter().all(|x| x.is_finite() && *x > 0.0)
        {
            return Err(invalid("uniform tube requires positive finite geometry and fluid properties"));
        }
        let area = core::f64::consts::PI * self.radius_m * self.radius_m;
        let z = density_kg_m3 * self.sound_speed_m_s / area;
        if !area.is_finite() || area <= 0.0 || !z.is_finite() || z <= 0.0 {
            return Err(invalid("uniform tube area or impedance is not representable"));
        }
        Ok(z)
    }
}

/// Only externally imposed inputs; the incoming wave always comes from the tube.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct TubeDrive {
    /// Upstream pressure [Pa], held across the aperture midpoint step.
    pub upstream_pressure_pa: f64,
    /// Additional flow supplied at the inlet by a separate body [m^3/s].
    pub body_flow_m3_s: f64,
}

/// One completely accepted aperture-and-propagation step.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct TubeFrame {
    /// Accepted midpoint junction observables and final mechanical state.
    pub aperture: ApertureFrame,
    /// Same-step inlet/terminal observables and final wave storage.
    pub waveguide: WaveguideFrame,
    /// Sum of mechanical, contact and two-direction wave storage [J].
    pub stored_energy_j: f64,
    /// Change of that total storage [J].
    pub storage_change_j: f64,
    /// Jet, mechanical/contact and terminal dissipation [J].
    pub dissipated_energy_j: f64,
    /// Work supplied from upstream, P_upstream (U_bore - U_body) dt [J].
    pub upstream_work_j: f64,
    /// Work supplied by the body-flow source, P_bore U_body dt [J].
    pub body_work_j: f64,
}
impl TubeFrame {
    /// Full discrete balance for this admitted tube/junction [J].
    /// This reports actual storage and separate supplies, not an energy repair.
    #[must_use]
    pub fn balance_residual_j(&self) -> f64 {
        self.storage_change_j + self.dissipated_energy_j
            - self.upstream_work_j - self.body_work_j
    }
}

/// Coupled runtime with one accepted clock and no external access to mutate
/// either participant independently. Cancellation is observed between samples;
/// a refused sample changes neither the mechanics nor the propagation state.
pub struct ApertureTube {
    aperture: DynamicAperture,
    line: PassiveWaveguide,
    spec: UniformTubeSpec,
    represented_length_m: f64,
}

impl ApertureTube {
    /// Couple a new aperture to an initially quiescent physical tube.
    /// Initial aperture motion/contact is allowed and retains its own energy;
    /// a previously advanced aperture is refused to prevent a clock reset.
    ///
    /// # Errors
    /// Geometry/load mismatch, unapproved transit rounding, memory admission,
    /// nonzero aperture clock or nonfinite derived quantities.
    pub fn new(aperture: DynamicAperture, spec: UniformTubeSpec) -> Result<Self, AcousticRealizeError> {
        if aperture.accepted_steps() != 0 {
            return Err(invalid("coupling a quiescent tube requires an unadvanced aperture"));
        }
        let a = aperture.spec();
        let z = spec.characteristic_impedance(a.density_kg_m3)?;
        if a.impedance_pa_s_m3.to_bits() != z.to_bits() {
            return Err(invalid("aperture impedance must equal tube.characteristic_impedance(density)"));
        }
        if !spec.max_length_error_m.is_finite() || spec.max_length_error_m < 0.0 {
            return Err(invalid("uniform tube requires a finite nonnegative length-error allowance"));
        }
        let cell = spec.sound_speed_m_s * a.time_step_s;
        let samples = (spec.length_m / cell).round();
        if !cell.is_finite() || cell <= 0.0 || !samples.is_finite()
            || samples < 1.0 || samples >= usize::MAX as f64
        {
            return Err(invalid("uniform tube transit does not fit a positive integer delay"));
        }
        let represented_length_m = samples * cell;
        if !represented_length_m.is_finite()
            || (represented_length_m - spec.length_m).abs() > spec.max_length_error_m
        {
            return Err(invalid("uniform tube integer transit exceeds the declared length-error allowance"));
        }
        let line = PassiveWaveguide::new(WaveguideSpec {
            one_way_samples: samples as usize,
            impedance_pa_s_m3: z,
            time_step_s: a.time_step_s,
            reflection: spec.terminal_reflection,
            max_memory_bytes: spec.max_wave_memory_bytes,
        }).map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
        Ok(Self { aperture, line, spec, represented_length_m })
    }

    /// Read the actual accepted mechanical state and explicit contact law.
    #[must_use]
    pub const fn aperture(&self) -> &DynamicAperture { &self.aperture }

    /// Requested physical geometry and tolerance, preserved unchanged.
    #[must_use]
    pub const fn spec(&self) -> &UniformTubeSpec { &self.spec }

    /// Geometry actually represented by the admitted integer transit [m].
    #[must_use]
    pub const fn represented_length_m(&self) -> f64 { self.represented_length_m }

    /// One-way propagation delay; a reflection returns after twice this count.
    #[must_use]
    pub fn one_way_samples(&self) -> usize { self.line.spec().one_way_samples }

    /// Total current mechanical/contact/wave energy [J].
    #[must_use]
    pub fn stored_energy_j(&self) -> f64 {
        self.aperture.stored_energy_j() + self.line.stored_energy_j()
    }

    /// Increase the total sample budget without disturbing waves or mechanics.
    ///
    /// # Errors
    /// Same explicit, finite-horizon admission as the aperture.
    pub fn extend_step_budget(&mut self, total: u64) -> Result<(), AcousticRealizeError> {
        self.aperture.extend_step_budget(total)
    }

    /// Solve the aperture against the retained incoming wave, then propagate.
    /// The inlet power is internal: it cancels when the two stores are combined.
    /// A body-flow source has a SEPARATE supply term, never uncredited energy.
    ///
    /// # Errors
    /// Input, budget, solver, propagation or combined observation refusal. No
    /// participant advances until all candidate observations are finite.
    pub fn step(&mut self, drive: TubeDrive) -> Result<TubeFrame, AcousticRealizeError> {
        let aperture = self.aperture.preview_step(ApertureDrive {
            upstream_pressure_pa: drive.upstream_pressure_pa,
            incoming_pressure_pa: self.line.incoming_pressure_pa(),
            body_flow_m3_s: drive.body_flow_m3_s,
        })?;
        let waveguide = self.line.preview_step(aperture.outgoing_pressure_pa)
            .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
        let dt = self.aperture.spec().time_step_s;
        let frame = TubeFrame {
            aperture, waveguide,
            stored_energy_j: aperture.stored_energy_j + waveguide.stored_energy_j,
            storage_change_j: aperture.storage_change_j + waveguide.storage_change_j,
            dissipated_energy_j: aperture.dissipated_energy_j + waveguide.terminal_loss_j,
            upstream_work_j: drive.upstream_pressure_pa
                * ((aperture.bore_flow_m3_s - drive.body_flow_m3_s) * dt),
            body_work_j: aperture.bore_pressure_pa * (drive.body_flow_m3_s * dt),
        };
        if ![frame.stored_energy_j, frame.storage_change_j, frame.dissipated_energy_j,
            frame.upstream_work_j, frame.body_work_j, frame.balance_residual_j()]
            .iter().all(|v| v.is_finite())
        {
            return Err(invalid("coupled tube observation left the finite set"));
        }
        // Repeats the same bounded immutable preview before the line's
        // infallible commit. On any error neither participant has changed.
        self.line.step(aperture.outgoing_pressure_pa)
            .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
        self.aperture.accept_frame(aperture);
        Ok(frame)
    }

    /// Advance a caller-sized block. Cancellation/exhaustion returns a complete
    /// prefix and leaves the remaining output slots untouched. Resume with that
    /// suffix; failed samples may be corrected and retried without a reset.
    ///
    /// # Errors
    /// Shape mismatch or per-step refusal. On refusal, the aperture's accepted
    /// count identifies exactly the completed prefix; both clocks agree.
    pub fn advance_block(
        &mut self, inputs: &[TubeDrive], out: &mut [TubeFrame], gate: &CancelGate,
    ) -> Result<ApertureProgress, AcousticRealizeError> {
        if inputs.len() != out.len() {
            return Err(invalid("coupled tube input/output block lengths must match"));
        }
        for (completed, (drive, slot)) in inputs.iter().zip(out.iter_mut()).enumerate() {
            if gate.is_requested() {
                return Ok(ApertureProgress { completed, terminal: ApertureTerminal::Cancelled });
            }
            if self.aperture.accepted_steps() >= self.aperture.spec().max_steps {
                return Ok(ApertureProgress { completed, terminal: ApertureTerminal::BudgetExhausted });
            }
            *slot = self.step(*drive)?;
        }
        Ok(ApertureProgress { completed: inputs.len(), terminal: ApertureTerminal::Complete })
    }
}
