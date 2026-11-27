//! The gap between a point carried by a planar body and a fixed straight face, as a function of
//! the body's pose: a geometry owner for a contact on a body that moves and turns.
//!
//! The body translates by `s` along a fixed unit axis `u` and turns by `theta` about a pivot that
//! the translation carries. A point fixed in the body, at offset `d` from the pivot at
//! `theta = 0`, is at `P = p0 + s u + R(theta) d`, and its gap to the face through `A` with unit
//! outward normal `n` is
//!
//! ```text
//! g = n . (P - A) = g0 + a s + c cos(theta) + k sin(theta),
//! ```
//!
//! with `g0 = n . (p0 - A)`, `a = n . u`, `c = n . d` and `k = n_y d_x - n_x d_y`. A positive gap
//! is open.
//!
//! Port `s` receives the rate of `s` [m/s] and returns a force [N], port `theta` receives the
//! rate of `theta` [rad/s] and returns a torque [N m], and port `gap` receives the effort `R` [N]
//! of what it is joined to and returns the rate of the gap [m/s]. The pose moves by `dt` times
//! the rates received. Over that move the gap's rate is `J . (ds, dtheta) / dt`, with `J` the
//! exact discrete gradient of `g`: `J_s = a`, and in the product-to-sum form
//!
//! ```text
//! J_theta = (h(theta1) - h(theta0)) / (theta1 - theta0)
//!         = sinc(dtheta / 2) (k cos(theta_m) - c sin(theta_m)),
//! ```
//!
//! with `h = c cos + k sin`, `theta_m` the midpoint and `sinc(x) = sin(x) / x`, so a pose that
//! turns by nothing loses no digit. The efforts are `-R J`: the ports' powers cancel, the element
//! stores nothing, and its receipt is lossless. A contact bonded at `gap` (the end that closes as
//! the gap shrinks) and started at the closure `-g` integrates `-g` at the pose, to rounding, at
//! every step: its closure is the geometry's, not an integral of a frozen slope.
//!
//! No claim: one point and one straight face of unbounded extent. A face with ends, a curved face
//! or a face on a second moving body is not this element. The gap's rate is the step's exact
//! mean, and its end does not follow the midpoint relation (`EndOfStep::Untold`): it is joined
//! to a load line only where the line reads mean flows only, as the contact does.

use fs_exec::Cx;
use fs_exec::solver::snapshot_v2::StateIdentityCharterV2;
use fs_ports::PortKind;
use fs_ports::component::{
    Class, Component, EndOfStep, PortRole, PortSpec, PortTangent, Proposal, Receipt,
};

use super::{ElementError, finite, inputs_finite, port, positive, roundoff};

/// The pose of a [`PoseGap`]'s body: its translation along the axis [m] and its turn [rad].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Pose {
    /// Translation along the axis [m].
    pub s_m: f64,
    /// Turn about the pivot [rad].
    pub theta_rad: f64,
}

/// The gap from a point on a translating, turning planar body to a fixed straight face: see the
/// [module documentation](self).
#[derive(Debug, Clone)]
pub struct PoseGap {
    /// `g0 = n . (p0 - A)` [m].
    offset_m: f64,
    /// `a = n . u`.
    along: f64,
    /// `c = n . d` [m].
    cos_m: f64,
    /// `k = n_y d_x - n_x d_y` [m].
    sin_m: f64,
    ports: [PortSpec; 3],
}

impl PoseGap {
    /// The gap from the point at offset `point_m` from the pivot (at `theta = 0`) of a body whose
    /// pivot is at `pivot_m` when `s = 0` and moves along `axis`, to the face through `face_m`
    /// with outward normal `normal`. `axis` and `normal` are scaled to unit length. Its input
    /// order is `[rate of s, rate of theta, effort at gap]`, its output order `[force at s,
    /// torque at theta, rate of the gap]`.
    ///
    /// # Errors
    /// [`ElementError::Parameter`] when a coordinate is not finite, when `axis` or `normal` has no
    /// finite nonzero length, or when a derived coefficient of the gap is not finite.
    pub fn new(
        pivot_m: [f64; 2],
        axis: [f64; 2],
        point_m: [f64; 2],
        face_m: [f64; 2],
        normal: [f64; 2],
    ) -> Result<Self, ElementError> {
        let coordinates = [pivot_m, point_m, face_m].into_iter().flatten();
        for value in coordinates.chain(axis).chain(normal) {
            finite("pose gap coordinate", value)?;
        }
        let [ux, uy] = unit("pose gap axis length", axis)?;
        let [nx, ny] = unit("pose gap face normal length", normal)?;
        let [dx, dy] = point_m;
        let offset = nx * (pivot_m[0] - face_m[0]) + ny * (pivot_m[1] - face_m[1]);
        Ok(Self {
            offset_m: finite("pose gap offset [m]", offset)?,
            along: nx * ux + ny * uy,
            cos_m: finite("pose gap cosine coefficient [m]", nx * dx + ny * dy)?,
            sin_m: finite("pose gap sine coefficient [m]", ny * dx - nx * dy)?,
            ports: [
                port(
                    "s",
                    Some(PortKind::MechanicalForceVelocity),
                    PortRole::FlowIn,
                ),
                port(
                    "theta",
                    Some(PortKind::RotationalTorqueAngularVelocity),
                    PortRole::FlowIn,
                ),
                port(
                    "gap",
                    Some(PortKind::MechanicalForceVelocity),
                    PortRole::EffortIn,
                ),
            ],
        })
    }

    /// The gap [m] at `pose`: `g0 + a s + c cos(theta) + k sin(theta)`.
    #[must_use]
    pub fn gap_m(&self, pose: &Pose) -> f64 {
        let (sin, cos) = pose.theta_rad.sin_cos();
        self.offset_m + self.along * pose.s_m + self.cos_m * cos + self.sin_m * sin
    }

    /// `J_theta` over a turn `turn` from `theta0`, and its derivative with respect to the turn.
    fn turn_slope(&self, theta0: f64, turn: f64) -> (f64, f64) {
        let half = 0.5 * turn;
        let (sin, cos) = (theta0 + half).sin_cos();
        // h'(theta_m) and h''(theta_m) = -h(theta_m).
        let slope = self.sin_m * cos - self.cos_m * sin;
        let curve = -(self.cos_m * cos + self.sin_m * sin);
        let (sinc, sinc_slope) = sinc(half);
        (sinc * slope, 0.5 * (sinc_slope * slope + sinc * curve))
    }

    fn validate(dt_s: f64, old: &Pose, inputs: &[f64]) -> Result<(), ElementError> {
        positive("pose gap timestep [s]", dt_s)?;
        finite("pose gap translation [m]", old.s_m)?;
        finite("pose gap turn [rad]", old.theta_rad)?;
        inputs_finite(inputs)
    }
}

/// `v` over its length, refused when the length is not finite and positive.
fn unit(name: &'static str, v: [f64; 2]) -> Result<[f64; 2], ElementError> {
    let length = positive(name, v[0].hypot(v[1]))?;
    Ok([v[0] / length, v[1] / length])
}

/// `sin(x) / x` and its derivative, each from its series near zero, where the quotients lose
/// their digits: the first below `1e-4`, the second below `1e-2`, with every term the error
/// needs.
fn sinc(x: f64) -> (f64, f64) {
    let x2 = x * x;
    let value = if x.abs() < 1.0e-4 {
        1.0 - x2 / 6.0
    } else {
        x.sin() / x
    };
    let slope = if x.abs() < 1.0e-2 {
        x * (-1.0 / 3.0 + x2 * (1.0 / 30.0 - x2 / 840.0))
    } else {
        (x * x.cos() - x.sin()) / x2
    };
    (value, slope)
}

impl Component for PoseGap {
    type State = Pose;
    type Error = ElementError;

    const STATE_CHARTER: Option<StateIdentityCharterV2> = Some(StateIdentityCharterV2 {
        owner: "fs-system::elements::PoseGap",
        state_family: "planar-pose",
        schema_grammar: "s_m f64, theta_rad f64; both finite",
        codec_grammar: "two u64 words: the bits of s_m then theta_rad",
        codec_version: 1,
    });

    fn encode_state(&self, state: &Pose, words: &mut Vec<u64>) -> bool {
        if !(state.s_m.is_finite() && state.theta_rad.is_finite()) {
            return false;
        }
        words.extend([state.s_m.to_bits(), state.theta_rad.to_bits()]);
        true
    }

    fn decode_state(&self, words: &[u64]) -> Option<Pose> {
        let &[s, theta] = words else {
            return None;
        };
        let pose = Pose {
            s_m: f64::from_bits(s),
            theta_rad: f64::from_bits(theta),
        };
        (pose.s_m.is_finite() && pose.theta_rad.is_finite()).then_some(pose)
    }

    fn identity_words(&self, words: &mut Vec<u64>) {
        // The four coefficients decide every method; the ports are fixed.
        words.extend([self.offset_m, self.along, self.cos_m, self.sin_m].map(f64::to_bits));
    }

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::Smooth
    }

    fn stored_energy_j(&self, _state: &Pose) -> Result<f64, ElementError> {
        Ok(0.0)
    }

    fn end_of_step(&self, _port: usize) -> EndOfStep {
        EndOfStep::Untold
    }

    fn propose(
        &self,
        old: &Pose,
        dt_s: f64,
        inputs: &[f64],
        new: &mut Pose,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, ElementError> {
        cx.checkpoint()?;
        Self::validate(dt_s, old, inputs)?;
        let [rate_s, rate_theta, effort] = [inputs[0], inputs[1], inputs[2]];
        let turn = dt_s * rate_theta;
        let (slope, _) = self.turn_slope(old.theta_rad, turn);
        let end = Pose {
            s_m: finite("pose gap translation [m]", old.s_m + dt_s * rate_s)?,
            theta_rad: finite("pose gap turn [rad]", old.theta_rad + turn)?,
        };
        let response = [
            finite("pose gap force [N]", -effort * self.along)?,
            finite("pose gap torque [N m]", -effort * slope)?,
            finite(
                "pose gap rate [m/s]",
                self.along * rate_s + slope * rate_theta,
            )?,
        ];
        // Each port's work by its own size: they cancel, so their sum says nothing of rounding.
        let work = inputs
            .iter()
            .zip(&response)
            .fold(0.0, |sum, (input, output)| {
                sum + (dt_s * input * output).abs()
            });
        let allowance = finite("pose gap energy allowance [J]", roundoff(work))?;
        *new = end;
        outputs.copy_from_slice(&response);
        Ok(Proposal::Step(Receipt::lossless(0.0, allowance)))
    }
}

impl PortTangent for PoseGap {
    fn jvp(
        &self,
        old: &Pose,
        dt_s: f64,
        inputs: &[f64],
        d_inputs: &[f64],
        d_outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), ElementError> {
        cx.checkpoint()?;
        Self::validate(dt_s, old, inputs)?;
        inputs_finite(d_inputs)?;
        let [rate_theta, effort] = [inputs[1], inputs[2]];
        let (slope, by_turn) = self.turn_slope(old.theta_rad, dt_s * rate_theta);
        // d J_theta / d rate_theta = dt by_turn.
        let bend = dt_s * by_turn;
        d_outputs[0] = -self.along * d_inputs[2];
        d_outputs[1] = -effort * bend * d_inputs[1] - slope * d_inputs[2];
        d_outputs[2] = self.along * d_inputs[0] + (slope + rate_theta * bend) * d_inputs[1];
        Ok(())
    }

    fn vjp(
        &self,
        old: &Pose,
        dt_s: f64,
        inputs: &[f64],
        bar_outputs: &[f64],
        bar_inputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), ElementError> {
        cx.checkpoint()?;
        Self::validate(dt_s, old, inputs)?;
        inputs_finite(bar_outputs)?;
        let [rate_theta, effort] = [inputs[1], inputs[2]];
        let (slope, by_turn) = self.turn_slope(old.theta_rad, dt_s * rate_theta);
        let bend = dt_s * by_turn;
        bar_inputs[0] = self.along * bar_outputs[2];
        bar_inputs[1] =
            -effort * bend * bar_outputs[1] + (slope + rate_theta * bend) * bar_outputs[2];
        bar_inputs[2] = -self.along * bar_outputs[0] - slope * bar_outputs[1];
        Ok(())
    }
}
