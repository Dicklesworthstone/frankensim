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
//! A face may be given ends ([`PoseGap::bounded`]): a span `[lo, hi]` of the point's coordinate
//! `z = t . (P - A)` along the face's tangent `t` (the normal turned a quarter turn
//! anticlockwise), and the depth of the band behind the face in which a point counts as pressed
//! into it. The face is then a thin plate, and with `x = n . (P - A)` and `y` the distance of `z`
//! past the nearer end,
//!
//! ```text
//! front:   lo <= z <= hi, x >= -depth   g = x
//! behind:  lo <= z <= hi, x <  -depth   g = -x             (under the plate, open)
//! outside: z < lo or z > hi             g = hypot(x, y)    (to the nearer end, open)
//! ```
//!
//! `g` and its gradient are continuous across front and outside where `x >= 0`, and across
//! behind and outside. A step from the front band to behind (pushed through) or across an end
//! inside the band (the end struck sideways) is refused: a flank is another face. A step within
//! the front or within the behind region takes the exact secant above, its sign flipped behind.
//! A step with an end outside takes the discrete gradient
//!
//! ```text
//! J = grad g(q_m) + (dg - grad g(q_m) . dq) W dq / (dq . W dq),    W = diag(1, r^2),
//! ```
//!
//! with `dq` the step's move of `(s, theta)`, `q_m` its midpoint, `dg = g(q1) - g(q0)` and
//! `r = |d|` the point's arm, so that `J . dq = dg` and the closure stays `-g`. Below
//! `|dq|_W = 256 eps l`, with `l` the size of the terms `g` sums at both ends of the step, the
//! correction is rounding and `J = grad g(q_m)`: the point then moves by at most
//! `sqrt(2) |dq|_W`, and where `g` is continuous along the step, `|dg - J . dq|` is at most
//! `2 sqrt(2) |dq|_W`. Such a step has `g >= 0` at both ends, so a contact started at `-g`
//! carries no load over it (`R = 0`), and its `J` matters only to the closure and the tangents.
//!
//! No claim: one point and one straight face, unbounded or a thin plate. A flank, a curved face
//! or a face on a second moving body is not this element. A step's regions are read at its two
//! ends: a point that enters and leaves a region within one step, through the band or out of
//! the span and back, is not seen. The gap's rate is the step's exact mean, and its
//! end does not follow the midpoint relation (`EndOfStep::Untold`): it is joined to a load line
//! only where the line reads mean flows only, as the contact does.

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
    /// `t . (p0 - A)` [m] and `t . u`: the terms of `z` that do not turn. Its turning terms are
    /// `-k cos(theta) + c sin(theta)`.
    tangent: [f64; 2],
    /// The face's ends and the depth of the band behind it, for a face with ends.
    span: Option<Span>,
    ports: [PortSpec; 3],
}

/// A face's ends along its tangent [m] and the depth of the band behind it [m].
#[derive(Debug, Clone, Copy)]
struct Span {
    ends_m: [f64; 2],
    depth_m: f64,
}

/// Where a pose puts the point against a face with ends.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Region {
    Front,
    Behind,
    Outside,
}

/// The gap at a pose: its region, `x`, `g`, the gradient of `g` by `(s, theta)` and its Hessian
/// `[ss, s theta, theta theta]`.
#[derive(Debug, Clone, Copy)]
struct Jet {
    region: Region,
    x_m: f64,
    value: f64,
    grad: [f64; 2],
    hess: [f64; 3],
}

/// A step's secant `J` and its derivative by the rates, `m[i][k] = d J_i / d rate_k`.
#[derive(Debug, Clone, Copy)]
struct Secant {
    j: [f64; 2],
    m: [[f64; 2]; 2],
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
        let along_face = nx * (pivot_m[1] - face_m[1]) - ny * (pivot_m[0] - face_m[0]);
        Ok(Self {
            offset_m: finite("pose gap offset [m]", offset)?,
            along: nx * ux + ny * uy,
            cos_m: finite("pose gap cosine coefficient [m]", nx * dx + ny * dy)?,
            sin_m: finite("pose gap sine coefficient [m]", ny * dx - nx * dy)?,
            tangent: [
                finite("pose gap tangent offset [m]", along_face)?,
                nx * uy - ny * ux,
            ],
            span: None,
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

    /// This gap with its face cut to the span `ends_m = [lo, hi]` along the face's tangent from
    /// the face point, with a band `depth_m` deep behind it: see the
    /// [module documentation](self).
    ///
    /// # Errors
    /// [`ElementError::Parameter`] when an end is not finite, when `lo` is not below `hi`, or
    /// when `depth_m` is not finite and positive.
    pub fn bounded(mut self, ends_m: [f64; 2], depth_m: f64) -> Result<Self, ElementError> {
        let [lo, hi] = [
            finite("pose gap face end [m]", ends_m[0])?,
            finite("pose gap face end [m]", ends_m[1])?,
        ];
        positive("pose gap face span [m]", hi - lo)?;
        let depth_m = positive("pose gap band depth [m]", depth_m)?;
        self.span = Some(Span { ends_m, depth_m });
        Ok(self)
    }

    /// The gap [m] at `pose`: `g0 + a s + c cos(theta) + k sin(theta)`, or for a face with ends
    /// the distance by region in the [module documentation](self).
    #[must_use]
    pub fn gap_m(&self, pose: &Pose) -> f64 {
        match &self.span {
            None => self.normal_m(pose),
            Some(span) => self.jet(span, pose).value,
        }
    }

    /// `x = n . (P - A)` [m] at `pose`.
    fn normal_m(&self, pose: &Pose) -> f64 {
        let (sin, cos) = pose.theta_rad.sin_cos();
        self.offset_m + self.along * pose.s_m + self.cos_m * cos + self.sin_m * sin
    }

    /// The gap's [`Jet`] at `pose` against a face with ends.
    fn jet(&self, span: &Span, pose: &Pose) -> Jet {
        let (sin, cos) = pose.theta_rad.sin_cos();
        // h = c cos + k sin turns x, and h' = k cos - c sin; z turns by -h', and dz/dtheta = h.
        let turned = self.cos_m * cos + self.sin_m * sin;
        let slope = self.sin_m * cos - self.cos_m * sin;
        let x = self.normal_m(pose);
        let z = self.tangent[0] + self.tangent[1] * pose.s_m - slope;
        let [lo, hi] = span.ends_m;
        let front = Jet {
            region: Region::Front,
            x_m: x,
            value: x,
            grad: [self.along, slope],
            hess: [0.0, 0.0, -turned],
        };
        if (lo..=hi).contains(&z) {
            if x >= -span.depth_m {
                return front;
            }
            return Jet {
                region: Region::Behind,
                value: -x,
                grad: [-self.along, -slope],
                hess: [0.0, 0.0, turned],
                ..front
            };
        }
        // Past `hi` the distance is `z - hi`, before `lo` it is `lo - z`.
        let (sign, end) = if z > hi { (1.0, hi) } else { (-1.0, lo) };
        let y = sign * (z - end);
        let [dx, dy] = [[self.along, slope], [sign * self.tangent[1], sign * turned]];
        let g = x.hypot(y);
        let grad = [(x * dx[0] + y * dy[0]) / g, (x * dx[1] + y * dy[1]) / g];
        let outer = |a: usize, b: usize| (dx[a] * dx[b] + dy[a] * dy[b] - grad[a] * grad[b]) / g;
        Jet {
            region: Region::Outside,
            x_m: x,
            value: g,
            grad,
            hess: [
                outer(0, 0),
                outer(0, 1),
                outer(1, 1) + (y * sign * slope - x * turned) / g,
            ],
        }
    }

    /// The sum of the sizes of the terms `g` sums at `pose` [m], which bounds its rounding.
    fn terms_m(&self, span: &Span, pose: &Pose) -> f64 {
        let [lo, hi] = span.ends_m;
        self.offset_m.abs()
            + self.tangent[0].abs()
            + (self.along.abs() + self.tangent[1].abs()) * pose.s_m.abs()
            + 2.0 * (self.cos_m.abs() + self.sin_m.abs())
            + lo.abs().max(hi.abs())
    }

    /// The secant of the step from `old` to `end` at `rates`, or the step's refusal.
    fn secant(
        &self,
        old: &Pose,
        end: &Pose,
        dt_s: f64,
        rates: [f64; 2],
    ) -> Result<Secant, ElementError> {
        let linear = |sign: f64| {
            let (slope, by_turn) = self.turn_slope(old.theta_rad, dt_s * rates[1]);
            Secant {
                j: [sign * self.along, sign * slope],
                m: [[0.0, 0.0], [0.0, sign * dt_s * by_turn]],
            }
        };
        let Some(span) = &self.span else {
            return Ok(linear(1.0));
        };
        let (from, to) = (self.jet(span, old), self.jet(span, end));
        let through = |reason: &str| ElementError::Law {
            reason: format!("pose gap step {reason}: a flank is another face"),
        };
        match (from.region, to.region) {
            (Region::Front, Region::Front) => Ok(linear(1.0)),
            (Region::Behind, Region::Behind) => Ok(linear(-1.0)),
            (Region::Front, Region::Behind) | (Region::Behind, Region::Front) => {
                Err(through("goes between the front band and behind the face"))
            }
            (Region::Front, Region::Outside) if from.x_m < 0.0 => {
                Err(through("crosses an end of the face inside the band"))
            }
            (Region::Outside, Region::Front) if to.x_m < 0.0 => {
                Err(through("crosses an end of the face inside the band"))
            }
            _ => Ok(self.discrete(span, old, [from, to], dt_s, rates)),
        }
    }

    /// The discrete gradient of a step with an end outside: see the
    /// [module documentation](self).
    fn discrete(
        &self,
        span: &Span,
        old: &Pose,
        ends: [Jet; 2],
        dt_s: f64,
        rates: [f64; 2],
    ) -> Secant {
        let step = [dt_s * rates[0], dt_s * rates[1]];
        let mid = self.jet(
            span,
            &Pose {
                s_m: old.s_m + 0.5 * step[0],
                theta_rad: old.theta_rad + 0.5 * step[1],
            },
        );
        let (g, h) = (mid.grad, mid.hess);
        let hessian = [[h[0], h[1]], [h[1], h[2]]];
        let arm2 = self.cos_m * self.cos_m + self.sin_m * self.sin_m;
        let w = [step[0], arm2 * step[1]];
        let norm = step[0] * w[0] + step[1] * w[1];
        let end = Pose {
            s_m: old.s_m + step[0],
            theta_rad: old.theta_rad + step[1],
        };
        let floor = 256.0 * f64::EPSILON * (self.terms_m(span, old) + self.terms_m(span, &end));
        if !(norm > floor * floor) {
            return Secant {
                j: g,
                m: hessian.map(|row| row.map(|entry| 0.5 * dt_s * entry)),
            };
        }
        let [start, at_end] = ends;
        let residual = (at_end.value - start.value) - (g[0] * step[0] + g[1] * step[1]);
        let ratio = residual / norm;
        // d J / d dq = H/2 + w rho^T / N - 2 D w w^T / N^2 + (D / N) W, with
        // rho = grad g(q1) - grad g(q_m) - H dq / 2.
        let rho: [f64; 2] = std::array::from_fn(|k| {
            at_end.grad[k] - g[k] - 0.5 * (hessian[k][0] * step[0] + hessian[k][1] * step[1])
        });
        let weight = [1.0, arm2];
        let m = std::array::from_fn(|i| {
            std::array::from_fn(|k| {
                let diagonal = if i == k { ratio * weight[i] } else { 0.0 };
                let bend = w[i] * (rho[k] - 2.0 * ratio * w[k]) / norm;
                dt_s * (0.5 * hessian[i][k] + bend + diagonal)
            })
        });
        Secant {
            j: [g[0] + ratio * w[0], g[1] + ratio * w[1]],
            m,
        }
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

    /// The step's end pose and its secant at `inputs`, after checking them.
    fn advance(
        &self,
        old: &Pose,
        dt_s: f64,
        inputs: &[f64],
    ) -> Result<(Pose, Secant), ElementError> {
        positive("pose gap timestep [s]", dt_s)?;
        finite("pose gap translation [m]", old.s_m)?;
        finite("pose gap turn [rad]", old.theta_rad)?;
        inputs_finite(inputs)?;
        let rates = [inputs[0], inputs[1]];
        let end = Pose {
            s_m: finite("pose gap translation [m]", old.s_m + dt_s * rates[0])?,
            theta_rad: finite("pose gap turn [rad]", old.theta_rad + dt_s * rates[1])?,
        };
        let secant = self.secant(old, &end, dt_s, rates)?;
        Ok((end, secant))
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
        // The four coefficients decide every method of a face without ends; the ports are fixed.
        words.extend([self.offset_m, self.along, self.cos_m, self.sin_m].map(f64::to_bits));
        if let Some(span) = &self.span {
            let [lo, hi] = span.ends_m;
            let terms = [self.tangent[0], self.tangent[1], lo, hi, span.depth_m];
            words.extend(terms.map(f64::to_bits));
        }
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
        let (end, Secant { j, .. }) = self.advance(old, dt_s, inputs)?;
        let [rate_s, rate_theta, effort] = [inputs[0], inputs[1], inputs[2]];
        let response = [
            finite("pose gap force [N]", -effort * j[0])?,
            finite("pose gap torque [N m]", -effort * j[1])?,
            finite("pose gap rate [m/s]", j[0] * rate_s + j[1] * rate_theta)?,
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
        let (_, Secant { j, m }) = self.advance(old, dt_s, inputs)?;
        inputs_finite(d_inputs)?;
        let [rate_s, rate_theta, effort] = [inputs[0], inputs[1], inputs[2]];
        // The efforts are -R J and the rate J . v, with J a function of the rates v.
        let loaded = |i: usize| -effort * m[i][0] * d_inputs[0] - effort * m[i][1] * d_inputs[1];
        let carried = |k: usize| j[k] + (rate_s * m[0][k] + rate_theta * m[1][k]);
        d_outputs[0] = loaded(0) - j[0] * d_inputs[2];
        d_outputs[1] = loaded(1) - j[1] * d_inputs[2];
        d_outputs[2] = carried(0) * d_inputs[0] + carried(1) * d_inputs[1];
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
        let (_, Secant { j, m }) = self.advance(old, dt_s, inputs)?;
        inputs_finite(bar_outputs)?;
        let [rate_s, rate_theta, effort] = [inputs[0], inputs[1], inputs[2]];
        let [bar_force, bar_torque, bar_rate] = [bar_outputs[0], bar_outputs[1], bar_outputs[2]];
        let loaded = |k: usize| -effort * m[0][k] * bar_force - effort * m[1][k] * bar_torque;
        let carried = |k: usize| j[k] + (rate_s * m[0][k] + rate_theta * m[1][k]);
        bar_inputs[0] = loaded(0) + carried(0) * bar_rate;
        bar_inputs[1] = loaded(1) + carried(1) * bar_rate;
        bar_inputs[2] = -j[0] * bar_force - j[1] * bar_torque;
        Ok(())
    }
}
