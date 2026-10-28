//! The existing jet/mouth equations on the original causal DG bore.
//!
//! The exact step-mean source law is E = H + Z Q. E is inlet pressure,
//! or inlet minus exterior mouth pressure for an explicit joint mouth load.
//! At the physical inlet characteristic impedance Zc = rho c / A, it gives
//! wave relation p- = (Z-Zc)/(Z+Zc) p+ + H Zc/(Z+Zc).
//! This adapter changes only the bore owner. It adds no jet, propagator,
//! delayed feedback sample or radiation model.
//!
//! The retained jet card and authored receptivity remain exploration-only.
//! This is not the resolved compressible zone/duct interface of M3.3-W1 and
//! cannot qualify a flute acceptance render. Optional shared rational outlet
//! reaction retains its own physical source and reports; inlet pressure remains
//! an internal observation. Solid wall motion and live bore geometry changes
//! require their own coupled owners.

use super::network::FlueResonator;
use super::{FlueError, FlueGeometry, FlueVoice, MouthStep};
use crate::dg_port::{DgPortLoad, DgTerminalFrame, DgDrivenMouthFrame};
use crate::bernoulli_aperture::radiation::environment::multiport::bake::mixed::DrivenMouthSurface;
use crate::wind_line::WindBoundary;
use fs_aeroac::jetcard::JetCard;
use fs_duct::dg::DgStepOutput;
use fs_exec::CancelGate;

fn invalid(what: &'static str) -> FlueError {
    FlueError::Invalid { what }
}

fn poll(gate: &CancelGate) -> Result<(), FlueError> {
    if gate.is_requested() {
        Err(FlueError::Physics(
            crate::acoustic_realize::AcousticRealizeError::Cancelled,
        ))
    } else {
        Ok(())
    }
}

/// Original geometric/material bore and explicit numerical/source limits.
/// The gas at the actual inlet supplies jet density, viscosity and sound speed.
#[derive(Clone, Debug)]
pub struct FlueDgSpec {
    /// The original bounded cold bore and exact propagator, prepared once.
    pub load: DgPortLoad,
    /// Declared minimum jet speed U_min [m/s]; it sizes the retained time history. Below it, a
    /// jet whose lookup would read motion older than that history carries no deflection, and the
    /// windway's steady flow splits at zero deflection (ruling 8132). A validity edge, not a
    /// stability limit. Zero blowing and the unsteady windway's release keep their own paths.
    pub minimum_jet_speed_m_s: f64,
    /// Incompressible source Mach ceiling, in (0, 0.3].
    pub maximum_jet_mach: f64,
    /// Maximum retained convection-history samples, including interpolation.
    pub max_jet_history_samples: usize,
}

/// One accepted bore step in the same wave coordinates used by the source.
/// Pressure and flow are means/held values for THIS step, not terminal sound.
#[derive(Clone, Debug, PartialEq)]
pub struct FlueDgFrame {
    pub inlet_pressure_pa: f64,
    pub inlet_flow_m3_s: f64,
    /// Independently integrated original DG inlet work for this step [J].
    pub inlet_work_j: f64,
    pub outgoing_pressure_pa: f64,
    pub incoming_pressure_pa: f64,
    /// Actual axial outlet exchange, when loaded. Also included in terminals;
    /// it is a convenience observation, not another energy contribution.
    pub terminal: Option<DgTerminalFrame>,
    /// Every DG terminal exchange in its original relative MIMO channel order. These
    /// are neither inferred heat nor fabricated terminal storage certificates.
    pub terminals: Vec<DgTerminalFrame>,
    /// Original joint exterior mouth channel, when selected. Its pressure and
    /// signed work are separate from the actual bore observations above.
    pub driven_mouth: Option<DgDrivenMouthFrame>,
    /// Original DG inlet plus exterior mouth work minus solved source work [J].
    pub port_work_residual_j: f64,
    /// External ports and summed original bore/branch ledgers. Dense loss is still inferred
    /// where its balance_checked field is false.
    pub duct: DgStepOutput,
}

/// The original DG owner and optional causal terminal under the common flue loop.
pub struct FlueDg {
    load: DgPortLoad,
    characteristic_impedance_pa_s_m3: f64,
    maximum_jet_mach: f64,
    last: Option<FlueDgFrame>,
}

impl FlueDg {
    pub const fn load(&self) -> &DgPortLoad {
        &self.load
    }

    pub const fn last_frame(&self) -> Option<&FlueDgFrame> {
        self.last.as_ref()
    }

    pub const fn characteristic_impedance_pa_s_m3(&self) -> f64 {
        self.characteristic_impedance_pa_s_m3
    }

    pub fn stored_energy_j(&self) -> f64 {
        self.load.stored_energy_j()
    }

    pub(super) fn boundary(&self) -> Result<WindBoundary, FlueError> {
        let relation = self.load.port_response().map_err(FlueError::Physics)?;
        let zc = self.characteristic_impedance_pa_s_m3;
        let denominator = relation.impedance_pa_s_m3 + zc;
        if !denominator.is_finite() || denominator <= 0.0 {
            return Err(invalid("DG characteristic boundary has nonfinite impedance"));
        }
        let boundary = WindBoundary {
            direct_reflection: (relation.impedance_pa_s_m3 - zc) / denominator,
            history_pressure_pa: relation.history_pressure_pa * (zc / denominator),
        };
        if !boundary.direct_reflection.is_finite()
            || !boundary.history_pressure_pa.is_finite()
        {
            return Err(invalid("DG characteristic boundary left the finite domain"));
        }
        Ok(boundary)
    }

    pub(super) fn incoming(&self) -> f64 {
        self.last.as_ref().map_or(0.0, |frame| frame.incoming_pressure_pa)
    }

    pub(super) fn push(&mut self, mouth: MouthStep) -> Result<f64, FlueError> {
        let relation = self.load.port_response().map_err(FlueError::Physics)?;
        let dt = self.load.propagator().dt();
        let mut outgoing_pressure_pa = mouth.outgoing_pa;
        let inlet_flow_m3_s = mouth.flow_m3_s;
        let mut inlet_pressure_pa = mouth.pressure_pa;
        let mut incoming_pressure_pa = mouth.incoming_pa;
        let expected_pressure = relation.history_pressure_pa
            + relation.impedance_pa_s_m3 * inlet_flow_m3_s;
        if ![
            outgoing_pressure_pa,
            inlet_flow_m3_s,
            inlet_pressure_pa,
            incoming_pressure_pa,
        ]
        .iter()
        .all(|value| value.is_finite())
        {
            return Err(invalid("DG source wave implies nonfinite pressure or flow"));
        }
        let inlet = self.load.inlet();
        let trial = self.load.preview_step(inlet_flow_m3_s).map_err(FlueError::Physics)?;
        let output = trial.output();
        let driven_mouth = trial.driven_mouth_frame().copied();
        let actual_inlet = output.port_outputs.iter().find(|port| port.id == inlet)
            .ok_or_else(|| invalid("DG source lost its original inlet observation"))?.conjugate_rate;
        let mouth_pressure = driven_mouth.as_ref().map_or(0.0, |frame| frame.mean_pressure_pa);
        let pressure_scale = relation.history_pressure_pa.abs()
            + (relation.impedance_pa_s_m3 * inlet_flow_m3_s).abs()
            + inlet_pressure_pa.abs() + mouth_pressure.abs();
        let work = inlet_pressure_pa * (inlet_flow_m3_s * dt);
        let port_work_residual_j = trial.source_work_j() - work;
        let work_scale = trial.source_work_j().abs() + work.abs()
            + driven_mouth.as_ref().map_or(0.0, |frame| frame.load_work_j.abs());
        if (expected_pressure - inlet_pressure_pa).abs()
            > 1e-10 * pressure_scale.max(1e-300)
            || (actual_inlet - mouth_pressure - inlet_pressure_pa).abs()
                > 1e-10 * pressure_scale.max(1e-300)
            || !work.is_finite()
            || !port_work_residual_j.is_finite()
            || port_work_residual_j.abs() > 1e-10 * work_scale.max(1e-300)
        {
            return Err(invalid("DG source wave and original mean port work disagree"));
        }
        // The nonlinear solve uses the source pressure jump. Exterior records
        // and launched-wave observations retain actual DG inlet pressure.
        if driven_mouth.is_some() {
            inlet_pressure_pa = actual_inlet;
            outgoing_pressure_pa = 0.5 * (actual_inlet
                + self.characteristic_impedance_pa_s_m3 * inlet_flow_m3_s);
            incoming_pressure_pa = 0.5 * (actual_inlet
                - self.characteristic_impedance_pa_s_m3 * inlet_flow_m3_s);
            if !outgoing_pressure_pa.is_finite() || !incoming_pressure_pa.is_finite() {
                return Err(invalid("joint mouth produces nonfinite physical inlet waves"));
            }
        }
        // All fallible checks precede the common owner's atomic publication.
        let inlet_work_j = trial.inlet_work_j();
        let terminal = trial.terminal_frame().copied();
        let terminals = trial.terminal_frames().to_vec();
        let duct = trial.accept();
        self.last = Some(FlueDgFrame {
            inlet_pressure_pa,
            inlet_flow_m3_s,
            inlet_work_j,
            outgoing_pressure_pa,
            incoming_pressure_pa,
            port_work_residual_j,
            terminal,
            terminals,
            driven_mouth,
            duct,
        });
        Ok(incoming_pressure_pa)
    }

    // The voice owns the declared minimum speed and its state below it (ruling 8132).
    pub(super) fn validate_jet_speed(&self, speed: f64) -> Result<(), FlueError> {
        if !speed.is_finite()
            || speed > self.maximum_jet_mach * self.load.duct().inlet_gas().sound_speed
        {
            return Err(invalid("DG flue exceeds its declared low-Mach jet domain"));
        }
        Ok(())
    }
}

impl FlueVoice {
    /// Compose the original jet, optional windway and mouth laws with the
    /// declared material-gas DG bore. The complete causal wall-memory state
    /// participates in the same-step mouth solve through its exact mean port.
    ///
    /// Gas is read from the actual inlet; a conflicting mouth bore area refuses.
    /// Every original acoustic and diffusion mode remains retained. The shared
    /// one-port owner checks the total DOF allowance and mechanical Nyquist band.
    /// step_block keeps its existing launched-wave observation. Inspect
    /// dg().last_frame() for physical mean inlet pressure, flow and work.
    ///
    /// The source is the original estimate-only jet card, not resolved flow.
    /// No flute acceptance, exterior microphone or solid-wall claim follows.
    #[allow(clippy::too_many_arguments)]
    pub fn from_dg(
        card: Option<&JetCard>,
        displacement_gain: f64,
        seed_amplitude_rel: f64,
        geometry: FlueGeometry,
        source: FlueDgSpec,
        sample_rate_hz: u32,
        voice_seed: u64,
        voice_index: u32,
        gate: &CancelGate,
    ) -> Result<Self, FlueError> {
        poll(gate)?;
        let card = card.ok_or(FlueError::NoCard)?;
        let gas = *source.load.duct().inlet_gas();
        let area = source.load.duct().inlet_area_m2();
        let bound_mouth = source.load.radiation_record().and_then(|record| record.driven_mouth());
        if sample_rate_hz == 0
            || source.load.has_driven_mouth() != bound_mouth.is_some()
            || source.load.state().step_count() != 0
            || source.load.propagator().dt().to_bits()
                != (1.0 / f64::from(sample_rate_hz)).to_bits()
            || !source.minimum_jet_speed_m_s.is_finite()
            || source.minimum_jet_speed_m_s <= 0.0
            || !source.maximum_jet_mach.is_finite()
            || !(0.0 < source.maximum_jet_mach && source.maximum_jet_mach <= 0.3)
            || source.minimum_jet_speed_m_s > source.maximum_jet_mach * gas.sound_speed
            || !(4..=1_048_576).contains(&source.max_jet_history_samples)
            || bound_mouth.is_some_and(|mouth| mouth.size_m
                != [geometry.flue_width_m, geometry.cut_up_m] || geometry.mouth_loss_rel != 0.0)
            || !geometry.bore_area_m2.is_finite()
            || (geometry.bore_area_m2 - area).abs() > 32.0 * f64::EPSILON * area
        {
            return Err(invalid("DG flue needs matching inlet area and bounded physical jet controls"));
        }
        let minimum = source.minimum_jet_speed_m_s;
        let maximum_jet_mach = source.maximum_jet_mach;
        let capacity = source.max_jet_history_samples;
        let load = source.load;
        poll(gate)?;
        let mut voice = Self::with_resonator(
            Some(card),
            displacement_gain,
            seed_amplitude_rel,
            geometry,
            &gas,
            sample_rate_hz,
            voice_seed,
            voice_index,
            minimum,
            capacity,
            |zc| {
                Ok(FlueResonator::Dg(Box::new(FlueDg {
                    load,
                    characteristic_impedance_pa_s_m3: zc,
                    maximum_jet_mach,
                    last: None,
                })))
            },
        )?;
        if bound_mouth.is_some() {
            voice.inertance = 0.0;
            voice.resistance = 0.0;
        }
        poll(gate)?;
        Ok(voice)
    }

    /// Attach one complete mouth/terminal exterior image to a cold DG voice.
    /// The original source, internal separation and material transport remain;
    /// the joint exterior replaces both former independent passive mouth terms.
    pub(crate) fn attach_joint_mouth_radiation(
        &mut self, surface: DrivenMouthSurface, rate: u32, gate: &CancelGate,
    ) -> Result<(), FlueError> {
        if self.poisoned || self.samples != 0 || self.mouth_radiation.is_some()
            || self.geometry.mouth_loss_rel != 0.0
            || surface.size_m != [self.geometry.flue_width_m, self.geometry.cut_up_m]
        {
            return Err(invalid("joint mouth requires a cold actual rectangle, zero legacy loss and no independent self load"));
        }
        let FlueResonator::Dg(duct) = &mut self.line else {
            return Err(invalid("joint mouth source requires the original DG pressure-jump loop"));
        };
        let record = duct.load.radiation_record()
            .ok_or_else(|| invalid("joint mouth requires a computed terminal exterior"))?;
        if record.gas_record().state().pressure.to_bits() != self.gas.pressure.to_bits() {
            return Err(invalid("joint mouth and bore require equal static pressure"));
        }
        let candidate = duct.load.clone().with_driven_mouth(surface, rate, gate)
            .map_err(|error| FlueError::Physics(
                crate::acoustic_realize::AcousticRealizeError::Nonlinear(error.to_string())))?;
        duct.load = candidate;
        self.inertance = 0.0;
        self.resistance = 0.0;
        Ok(())
    }

    /// The accepted DG state and source-bound physical bore; absent on other
    /// propagation paths. The source and returned load have no mutable escape.
    pub fn dg(&self) -> Option<&FlueDg> {
        match &self.line {
            FlueResonator::Dg(duct) => Some(duct),
            _ => None,
        }
    }
}

/// Existing scheduled pressure/PCM interface for explicit mean-inlet observation.
pub mod performance;

#[cfg(test)]
mod tests;
