//! Sampled speed control as a pure Component, with command holds between samples.
//!
//! This adapter delegates the PI, directional torque/power limits and tracking
//! anti-windup to [`SpeedController`]. Its actuator is an ideal torque source:
//! the admitted torque request is the achieved command passed to anti-windup.
//! A real machine's current allocation, voltage limits and thermal policy still
//! belong to their original owners; this adapter does not simulate those.

use core::num::NonZeroUsize;

use fs_exec::Cx;
use fs_exec::solver::snapshot_v2::StateIdentityCharterV2;
use fs_math::dd::Dd;
use fs_ports::StableId;
use fs_ports::component::{Class, Component, PortSpec, PortTangent, Proposal, Receipt};

use crate::ControlError;
use crate::speed::{
    ShaftPowerLimits, SpeedCommand, SpeedController, SpeedControllerParameters,
    SpeedControllerState, SpeedDemand,
};

/// Accepted digital state. The PI state and held command are not physical energy.
///
/// A sample is committed with the first accepted plant substep of its period.
/// Further substeps retain that command and do not run the PI again. All fields
/// are part of the Component codec, including the compensated sample phase.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct SampledSpeedState {
    pi: SpeedControllerState,
    held_torque_nm: f64,
    phase_s: Dd,
    samples: u64,
}

impl SampledSpeedState {
    /// Start at a sample boundary with an authored initial integral contribution.
    pub fn new(integral_torque_nm: f64) -> Result<Self, ControlError> {
        if !integral_torque_nm.is_finite() {
            return Err(invalid("initial integral torque must be finite"));
        }
        Ok(Self {
            pi: SpeedControllerState { integral_torque_nm },
            ..Self::default()
        })
    }

    /// The accepted integral torque contribution [N m].
    #[must_use]
    pub fn integral_torque_nm(&self) -> f64 {
        self.pi.integral_torque_nm
    }

    /// Command held over the most recently accepted interval [N m].
    #[must_use]
    pub fn held_torque_nm(&self) -> f64 {
        self.held_torque_nm
    }

    /// Number of accepted samples, independent of the plant's substep count.
    #[must_use]
    pub fn samples(&self) -> u64 {
        self.samples
    }

    /// Time since the last sample [s], zero at a new sample boundary.
    #[must_use]
    pub fn phase_s(&self) -> f64 {
        self.phase_s.to_f64()
    }
}

/// A fixed-period speed controller driving an ideal torque actuator.
///
/// Scalar signal inputs, in order: `reference` [rad/s], `feedforward` [N m],
/// `rotor_speed` [rad/s], `stator_speed` [rad/s]. Output `torque` [N m] can
/// drive a rotational `fs_system::elements::ControlledSource`.
///
/// Measurements must be the plant's accepted start-of-step values, for example
/// the `velocity` signal of a rotational `Mass::publishing`. The work-conjugate
/// mean speed of a power port is a different quantity. No sampling delay or
/// sensor noise is inferred. The control law runs only at sample boundaries;
/// the output is held across smaller plant steps and event retries.
/// The next sample is a step limit, so no accepted step straddles two commands.
/// Proposals remain pure under retries and checkpoint restoration. This adapter
/// supplies no adaptive error chart or time-order claim; generic library callers
/// must establish an appropriate metric and order for adaptive integration.
///
/// The port map is piecewise smooth. Its declared tangent is the active branch's
/// Newton derivative, taking the limited branch at a clipping boundary; it is
/// not a differentiability claim at the switching surface.
#[derive(Debug, Clone)]
pub struct SampledSpeedController {
    controller: SpeedController,
    limits: ShaftPowerLimits,
    ports: [PortSpec; 5],
}

fn invalid(what: &'static str) -> ControlError {
    ControlError::InvalidSample { what }
}

fn poll(cx: &Cx<'_>) -> Result<(), ControlError> {
    cx.checkpoint().map_err(|_| ControlError::Cancelled)
}

impl SampledSpeedController {
    /// Admit the original sampled law and fixed directional command capabilities.
    pub fn new(
        parameters: SpeedControllerParameters,
        limits: ShaftPowerLimits,
    ) -> Result<Self, ControlError> {
        let controller = SpeedController::new(parameters)?;
        // Use the source owner's admission for all four capabilities.
        controller.demand_controlled(
            &SpeedControllerState::default(),
            [0.0; 2],
            SpeedCommand {
                relative_speed_rad_s: 0.0,
                feedforward_torque_nm: 0.0,
            },
            limits,
            || false,
        )?;
        let input = |name| {
            PortSpec::signal_in(
                StableId::new(name).expect("static control port"),
                NonZeroUsize::MIN,
            )
        };
        Ok(Self {
            controller,
            limits,
            ports: [
                input("reference"),
                input("feedforward"),
                input("rotor_speed"),
                input("stator_speed"),
                PortSpec::signal_out(
                    StableId::new("torque").expect("static control port"),
                    NonZeroUsize::MIN,
                ),
            ],
        })
    }

    fn validate_state(&self, state: &SampledSpeedState) -> Result<(), ControlError> {
        let phase = state.phase_s;
        let period = self.controller.parameters().sample_period_s;
        if !state.pi.integral_torque_nm.is_finite()
            || !state.held_torque_nm.is_finite()
            || state.held_torque_nm > self.limits.positive_torque_limit_nm
            || state.held_torque_nm < -self.limits.negative_torque_limit_nm
            || !phase.hi.is_finite()
            || !phase.lo.is_finite()
            || Dd::from_pair(phase.hi, phase.lo) != phase
            || phase.lt(Dd::ZERO)
            || !phase.lt(Dd::from_f64(period))
            || (state.samples == 0 && (phase != Dd::ZERO || state.held_torque_nm != 0.0))
        {
            return Err(invalid("invalid sampled controller state or phase"));
        }
        Ok(())
    }

    fn remaining(&self, state: &SampledSpeedState) -> f64 {
        (Dd::from_f64(self.controller.parameters().sample_period_s) - state.phase_s).to_f64()
    }

    fn interval(&self, state: &SampledSpeedState, proposed_s: f64) -> f64 {
        let remaining = self.remaining(state);
        let leftover = remaining - proposed_s;
        let resolution = 8.0 * f64::EPSILON * self.controller.parameters().sample_period_s;
        // Honor the caller's strict ceiling without stranding a rounding-sized
        // interval before the sample. Two resolvable shorter steps reach it
        // instead; no extra time is integrated and no sample is taken early.
        if leftover > 0.0 && leftover <= resolution {
            proposed_s.min(0.5 * remaining)
        } else {
            proposed_s.min(remaining)
        }
    }

    fn validate(
        &self,
        state: &SampledSpeedState,
        dt_s: f64,
        inputs: &[f64],
    ) -> Result<(), ControlError> {
        self.validate_state(state)?;
        if !(dt_s.is_finite() && dt_s > 0.0)
            || inputs.len() != 4
            || inputs.iter().any(|x| !x.is_finite())
        {
            return Err(invalid(
                "positive finite step and four finite control signals required",
            ));
        }
        Ok(())
    }

    fn demand(
        &self,
        state: &SampledSpeedState,
        inputs: &[f64],
        cx: &Cx<'_>,
    ) -> Result<SpeedDemand, ControlError> {
        self.controller.demand_controlled(
            &state.pi,
            [inputs[2], inputs[3]],
            SpeedCommand {
                relative_speed_rad_s: inputs[0],
                feedforward_torque_nm: inputs[1],
            },
            self.limits,
            || cx.is_cancel_requested(),
        )
    }

    fn gradient(
        &self,
        state: &SampledSpeedState,
        dt_s: f64,
        inputs: &[f64],
        cx: &Cx<'_>,
    ) -> Result<[f64; 4], ControlError> {
        poll(cx)?;
        self.validate(state, dt_s, inputs)?;
        if dt_s > self.interval(state, dt_s) {
            return Err(invalid("a tangent cannot cross a controller sample"));
        }
        if state.phase_s != Dd::ZERO {
            return Ok([0.0; 4]);
        }
        let demand = self.demand(state, inputs, cx)?;
        let [lower, upper] = demand.torque_bounds_nm;
        let (bound, cap) = if demand.unlimited_torque_nm <= lower {
            (lower, self.limits.negative_torque_limit_nm)
        } else if demand.unlimited_torque_nm >= upper {
            (upper, self.limits.positive_torque_limit_nm)
        } else {
            let parameters = self.controller.parameters();
            let gain = parameters.proportional_gain_nm_s_per_rad
                + parameters.sample_period_s * parameters.integral_gain_nm_per_rad;
            if !gain.is_finite() {
                return Err(ControlError::NumericalRange {
                    what: "sampled control tangent",
                });
            }
            return Ok([gain, 1.0, -gain, gain]);
        };
        // A power-limited endpoint is +/- P/|omega|. Its derivative is
        // -torque/omega in either rotation direction. A torque cap is constant.
        let speed = demand.measured_relative_speed_rad_s;
        let derivative = if speed != 0.0 && bound.abs() < cap {
            -bound / speed
        } else {
            0.0
        };
        if !derivative.is_finite() {
            return Err(ControlError::NumericalRange {
                what: "power-limit control tangent",
            });
        }
        Ok([0.0, 0.0, derivative, -derivative])
    }
}

impl Component for SampledSpeedController {
    type State = SampledSpeedState;
    type Error = ControlError;

    const STATE_CHARTER: Option<StateIdentityCharterV2> = Some(StateIdentityCharterV2 {
        owner: "fs-control::component::SampledSpeedController",
        state_family: "sampled-speed-control",
        schema_grammar: "integral torque Nm; held torque Nm; normalized double-double phase seconds in [0,Ts); accepted sample count u64",
        codec_grammar: "five u64 words: integral bits, held bits, phase hi bits, phase lo bits, sample count",
        codec_version: 1,
    });

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::NonSmooth
    }

    fn stored_energy_j(&self, state: &SampledSpeedState) -> Result<f64, ControlError> {
        self.validate_state(state)?;
        Ok(0.0)
    }

    fn limit_dt_s(&self, state: &SampledSpeedState, proposed_s: f64) -> Result<f64, ControlError> {
        self.validate_state(state)?;
        if !(proposed_s.is_finite() && proposed_s > 0.0) {
            return Err(invalid(
                "positive finite proposed control interval required",
            ));
        }
        Ok(self.interval(state, proposed_s))
    }

    fn next_event_s(&self, state: &SampledSpeedState) -> Option<f64> {
        Some(self.remaining(state))
    }

    fn propose(
        &self,
        old: &SampledSpeedState,
        dt_s: f64,
        inputs: &[f64],
        new: &mut SampledSpeedState,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, ControlError> {
        poll(cx)?;
        self.validate(old, dt_s, inputs)?;
        if outputs.len() != 1 {
            return Err(invalid("one torque output required"));
        }
        let remaining = self.remaining(old);
        let interval = self.interval(old, dt_s);
        if dt_s > interval {
            return Ok(Proposal::Shorten { dt_s: interval });
        }
        let mut next = *old;
        if old.phase_s == Dd::ZERO {
            let demand = self.demand(old, inputs, cx)?;
            let proposal = demand
                .finish_controlled(demand.requested_torque_nm, || cx.is_cancel_requested())?;
            next.pi = proposal.next_state;
            next.held_torque_nm = proposal.achievable_torque_nm;
            next.samples = old
                .samples
                .checked_add(1)
                .ok_or(ControlError::NumericalRange {
                    what: "accepted control sample count",
                })?;
        }
        // Subtraction above rounds the distance once. Taking that exact
        // returned limit lands on the boundary, without a residual microstep.
        next.phase_s = if dt_s == remaining {
            Dd::ZERO
        } else {
            old.phase_s + Dd::from_f64(dt_s)
        };
        self.validate_state(&next)?;
        poll(cx)?;
        outputs[0] = next.held_torque_nm;
        *new = next;
        Ok(Proposal::Step(Receipt::lossless(0.0, 0.0)))
    }

    fn encode_state(&self, state: &SampledSpeedState, words: &mut Vec<u64>) -> bool {
        if self.validate_state(state).is_err() {
            return false;
        }
        words.extend([
            state.pi.integral_torque_nm.to_bits(),
            state.held_torque_nm.to_bits(),
            state.phase_s.hi.to_bits(),
            state.phase_s.lo.to_bits(),
            state.samples,
        ]);
        true
    }

    fn decode_state(&self, words: &[u64]) -> Option<SampledSpeedState> {
        let &[integral, held, hi, lo, samples] = words else {
            return None;
        };
        let state = SampledSpeedState {
            pi: SpeedControllerState {
                integral_torque_nm: f64::from_bits(integral),
            },
            held_torque_nm: f64::from_bits(held),
            phase_s: Dd {
                hi: f64::from_bits(hi),
                lo: f64::from_bits(lo),
            },
            samples,
        };
        self.validate_state(&state).ok().map(|()| state)
    }

    fn identity_words(&self, words: &mut Vec<u64>) {
        let p = self.controller.parameters();
        words.extend(
            [
                p.sample_period_s,
                p.proportional_gain_nm_s_per_rad,
                p.integral_gain_nm_per_rad,
                p.anti_windup_time_s,
                self.limits.positive_torque_limit_nm,
                self.limits.negative_torque_limit_nm,
                self.limits.motoring_power_w,
                self.limits.regeneration_power_w,
            ]
            .map(f64::to_bits),
        );
    }
}

impl PortTangent for SampledSpeedController {
    fn jvp(
        &self,
        old: &SampledSpeedState,
        dt_s: f64,
        inputs: &[f64],
        direction: &[f64],
        output: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), ControlError> {
        if direction.len() != 4 || output.len() != 1 || direction.iter().any(|x| !x.is_finite()) {
            return Err(invalid(
                "four finite control directions and one output required",
            ));
        }
        let gradient = self.gradient(old, dt_s, inputs, cx)?;
        let value = gradient
            .iter()
            .zip(direction)
            .fold(0.0, |sum, (g, d)| sum + g * d);
        if !value.is_finite() {
            return Err(ControlError::NumericalRange {
                what: "sampled control JVP",
            });
        }
        poll(cx)?;
        output[0] = value;
        Ok(())
    }

    fn vjp(
        &self,
        old: &SampledSpeedState,
        dt_s: f64,
        inputs: &[f64],
        output: &[f64],
        inputs_bar: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), ControlError> {
        if output.len() != 1 || inputs_bar.len() != 4 || !output[0].is_finite() {
            return Err(invalid(
                "one finite control cotangent and four inputs required",
            ));
        }
        let values = self.gradient(old, dt_s, inputs, cx)?.map(|g| g * output[0]);
        if values.iter().any(|x| !x.is_finite()) {
            return Err(ControlError::NumericalRange {
                what: "sampled control VJP",
            });
        }
        poll(cx)?;
        inputs_bar.copy_from_slice(&values);
        Ok(())
    }
}
