//! An ideal dry-friction slider between two ends, with constant static and
//! kinetic coefficients under a constant normal force.
//!
//! It is the element the non-smooth solve tier is exercised with, and
//! nothing more: friction that follows a material, a surface card, a
//! temperature or a slip history belongs to `fs-tribo`, whose law a
//! component of its own puts behind the same tier.
//!
//! # What a port tells the slider, and what it does not
//!
//! The elements of this crate use the midpoint rule: a mass returns its
//! MEAN velocity over the step on its port, and its velocity at the end of
//! the step is twice that mean less the velocity it started with. A
//! prescribed flow is another thing: a boundary's value is held over the
//! step, so its start, its mean and its end are that one value, and it may
//! jump from one step to the next. Friction has to know the slip a step
//! starts with, and a mean alone does not say which of the two an end is.
//!
//! The load line does. An end whose row of the line's matrix is zero does
//! not answer the slider's force at all: its flow is given, and held over
//! the step. The slider keeps the slip its last step ended with. Within a
//! step that slip goes to twice its mean less its start, whether an end
//! answers (the midpoint relation) or is held (it does not move). What a
//! given end can add is a jump between two steps: the slider remembers the
//! flow each given end was given, and when it comes back with another
//! value, the step starts with that much more slip. (The first version of
//! this element did not look for the jump. When a prescribed side jumped,
//! its stored slip was off by the jump for ever, alternating in sign, and
//! friction acted every other step.)
//!
//! # The law of a step
//!
//! `s0` is the slip the step starts with, `mean = r - c e` the load line for
//! the mean slip of the step, with `r` the released mean slip, `c > 0` the
//! compliance of what the slider is joined to, and `e` the friction force
//! it returns (positive along the slip it resists). `N` is the normal
//! force. The slider also keeps whether its last step ended stuck.
//!
//! Two obvious laws fail on step means. "Stuck means the mean slip is zero"
//! makes the end slip `-s0`: a thrown mass never comes to rest, its velocity
//! flips sign every step and nothing is dissipated. "Stuck means the end
//! slip is zero" alone gives a negative loss `e * mean` in a step in which
//! the drive reverses the body. The law takes the end-of-step arrest where
//! its loss is not negative, and the mean elsewhere.
//!
//! With `ea = (r - s0/2) / c` the force that ends the step at zero slip and
//! `sg = sign(s0)`:
//!
//! | state | condition | force `e` | the step ends |
//! |---|---|---|---|
//! | stuck, `s0 = 0` | `abs(r) <= c mu_s N` | `r / c` | stuck (held) |
//! | stuck, `s0 = 0` | otherwise | `mu_k N sign(r)` | sliding (broke away) |
//! | stuck, the side it is stuck to moved (`s0 != 0`) | `0 <= sg ea <= mu_s N` | `ea` | stuck, at zero slip (arrested) |
//! | stuck, `s0 != 0` | else `abs(r) <= c mu_s N` | `r / c` | stuck, zero mean slip (held) |
//! | stuck, `s0 != 0` | else | `mu_k N sign(r)` | sliding (broke away) |
//! | sliding | `0 <= sg ea <= mu_k N` | `ea` | stuck, at zero slip exactly (arrested) |
//! | sliding | `sg ea > mu_k N` | `mu_k N sg` | sliding the same way (slid) |
//! | sliding | `0 <= sg r < abs(s0)/2` | `0` | reversed, slowly (released) |
//! | sliding | `-c mu_k N <= sg r < 0` | `r / c` | reversed (turned) |
//! | sliding | `sg r < -c mu_k N` | `-mu_k N sg` | reversed |
//!
//! - The loss `e * mean` is not negative in any row, and the first law is
//!   exact on step means: what the two ends give up leaves through `heat`.
//! - `abs(e) <= mu_s N` stuck and `<= mu_k N` sliding.
//! - A mass thrown along a surface is arrested with the slider's slip at
//!   exactly zero and all of its kinetic energy given up as heat. (The
//!   mass's own velocity is then zero to the rounding its loop was closed
//!   to: the mass and the slider each keep their own record.)
//! - Stuck, the force that carries the body with a side that moved is taken
//!   up to the STATIC limit: a belt that speeds up step by step carries a
//!   block until the force it needs exceeds `mu_s N`, and a jump is caught
//!   within one step if `mu_s N` allows. Where that force would have a
//!   negative loss, the zero-mean hold is taken instead.
//! - In the sliding state `e` is continuous and does not decrease with `r`.
//!   It jumps at break-away when `mu_s > mu_k`, as the physics does, and in
//!   the stuck state with a moved side it jumps between its two holds.
//! - A step in which the slip reverses is first-order accurate. The
//!   "released" row gives the friction of such a step up altogether: the
//!   drive alone reverses the body there, a friction force that ended the
//!   step at zero slip would have a negative loss, and the row keeps the
//!   force continuous with its neighbours.
//!
//! With no compliance (`c <= 0`: both flows given, or two ends on one rigid
//! thing) the slip is what it is and the force is `mu_k N sign(r)`, none at
//! exactly zero slip.
//!
//! Not covered: a prescribed side whose velocity varies WITHIN a step is
//! seen as its step mean held over the step, so a stuck block follows a
//! ramp as a staircase; and the slider's record of an end that answers is
//! only as good as the midpoint relation is for what is attached there.

use fs_exec::Cx;
use fs_ports::PortKind;
use fs_ports::component::{Class, Component, LoadLine, PortRole, PortSpec, Proposal, Receipt};

use super::{ElementError, heat, inputs_finite, port, positive, roundoff, temperature};

/// A dry-friction slider between ports `a` and `b`, with the loss leaving
/// through port `heat`.
///
/// `a` and `b` receive the velocities [m/s] of the two sides and return
/// the friction force [N] and its opposite; the slip is `a - b`. `heat`
/// receives a temperature [K] and returns the heat rate into the slider
/// [W], which is minus the loss. Its state is [`Slip`]. Added with
/// `SystemBuilder::add_on_line`.
#[derive(Debug, Clone)]
pub struct Slider {
    /// Largest force held at zero slip [N].
    hold_n: f64,
    /// Force while sliding [N].
    slide_n: f64,
    ports: [PortSpec; 3],
}

/// What a [`Slider`] keeps between steps: the slip its last step ended
/// with [m/s], whether that step ended stuck, and, for each end whose flow
/// was given in that step, the flow it was given.
///
/// The slip is the flow port `a` receives less the flow port `b` receives,
/// with the signs the ports receive them in. A member of a common-flow
/// junction receives MINUS the flow of the port that sets it, so a slider
/// whose `a` is joined to a mass moving at `v`, with `b` held at rest,
/// starts with a slip of `-v`. Nothing can check a starting slip against
/// what is attached; a pair that starts stuck needs none, whatever the
/// flow it moves with.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Slip {
    /// Flow at `a` less flow at `b` at the end of the last step.
    ended_m_s: f64,
    stuck: bool,
    /// Per end, the flow it was given in the last step, if it was given
    /// one: a given end that comes back with another value has jumped by
    /// the difference, and the next step starts with that much more slip.
    given_m_s: [Option<f64>; 2],
}

impl Slip {
    /// Stuck: the two ends move together.
    #[must_use]
    pub fn stuck() -> Self {
        Self {
            ended_m_s: 0.0,
            stuck: true,
            given_m_s: [None; 2],
        }
    }

    /// Sliding with this slip (flow at `a` less flow at `b`).
    #[must_use]
    pub fn sliding(slip_m_s: f64) -> Self {
        Self {
            ended_m_s: slip_m_s,
            stuck: false,
            given_m_s: [None; 2],
        }
    }

    /// The slip the last step ended with: flow at `a` less flow at `b`.
    #[must_use]
    pub fn slip_m_s(&self) -> f64 {
        self.ended_m_s
    }

    /// Whether the last step ended stuck.
    #[must_use]
    pub fn is_stuck(&self) -> bool {
        self.stuck
    }
}

impl Default for Slip {
    fn default() -> Self {
        Self::stuck()
    }
}

/// Which row of the law a step took (see the module's table).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SliderRow {
    /// Stuck, and held with zero mean slip.
    Held,
    /// Stuck, and broke away.
    BrokeAway,
    /// Ended the step at zero slip: a sliding body arrested, or a stuck one
    /// carried with a side that moved.
    Arrested,
    /// Sliding on the same way.
    Slid,
    /// Sliding, reversed by the drive alone; no friction in the step.
    Released,
    /// Sliding, reversed with zero mean slip.
    Turned,
    /// Sliding, reversed against friction of the new direction.
    Reversed,
}

impl Slider {
    /// A slider pressed with `normal_n` [N], holding up to
    /// `static_coefficient * normal_n` at rest and resisting with
    /// `kinetic_coefficient * normal_n` while it slides.
    ///
    /// # Errors
    /// A normal force or a coefficient that is not positive and finite, or
    /// a kinetic coefficient above the static one.
    pub fn new(
        normal_n: f64,
        static_coefficient: f64,
        kinetic_coefficient: f64,
    ) -> Result<Self, ElementError> {
        let normal_n = positive("normal_n", normal_n)?;
        let holds = positive("static_coefficient", static_coefficient)?;
        let slides = positive("kinetic_coefficient", kinetic_coefficient)?;
        if slides > holds {
            return Err(ElementError::Parameter {
                name: "kinetic_coefficient",
                value: kinetic_coefficient,
            });
        }
        let kind = Some(PortKind::MechanicalForceVelocity);
        Ok(Self {
            hold_n: positive("static friction force", holds * normal_n)?,
            slide_n: positive("kinetic friction force", slides * normal_n)?,
            ports: [
                port("a", kind, PortRole::FlowIn),
                port("b", kind, PortRole::FlowIn),
                heat("heat", PortRole::EffortIn),
            ],
        })
    }

    /// The force, the mean slip and the row of one step: `start` the slip
    /// the step begins with, `released` the released mean slip, `c` the
    /// compliance of the slip to the force, `stuck` whether the last step
    /// ended stuck (a zero `start` is stuck whatever it says).
    ///
    /// A `c` that is not positive is no compliance: both sides are given,
    /// or what the slider is joined to is rigid to within rounding, which
    /// can leave a hair below zero. The slip is then what it is, and the
    /// force is the sliding one. The mean slip returned is still the point
    /// of the line as given, `released - c * force`, so that the loss is
    /// accounted at the flows the ports carry.
    ///
    /// A force that is a quotient `released / c` is kept inside the limit
    /// its row was chosen by: the quotient can exceed it by a unit in the
    /// last place where the test on `released` passed.
    #[must_use]
    pub fn on_line(&self, start: f64, released: f64, c: f64, stuck: bool) -> (f64, f64, SliderRow) {
        let (hold, slide) = (self.hold_n, self.slide_n);
        if !(c > 0.0) {
            let (force, row) = if released == 0.0 {
                (0.0, SliderRow::Held)
            } else {
                (slide * released.signum(), SliderRow::Slid)
            };
            return (force, released - c * force, row);
        }
        // The zero-mean hold and the break-away, shared by both stuck cases.
        let hold_or_break = || {
            if released.abs() <= c * hold {
                ((released / c).clamp(-hold, hold), 0.0, SliderRow::Held)
            } else {
                let force = slide * released.signum();
                (force, released - c * force, SliderRow::BrokeAway)
            }
        };
        if start == 0.0 {
            return hold_or_break();
        }
        let way = start.signum();
        let arrest = (released - 0.5 * start) / c;
        let along = way * arrest;
        if stuck {
            // The side it is stuck to moved. Carried with it up to the
            // static limit if that costs no negative loss; else held on
            // the mean; else it breaks away.
            return if along >= 0.0 && along <= hold {
                (arrest, 0.5 * start, SliderRow::Arrested)
            } else {
                hold_or_break()
            };
        }
        if along >= 0.0 && along <= slide {
            (arrest, 0.5 * start, SliderRow::Arrested)
        } else if along > slide {
            let force = slide * way;
            (force, released - c * force, SliderRow::Slid)
        } else if way * released >= 0.0 {
            (0.0, released, SliderRow::Released)
        } else if way * released >= -c * slide {
            ((released / c).clamp(-slide, slide), 0.0, SliderRow::Turned)
        } else {
            let force = -slide * way;
            (force, released - c * force, SliderRow::Reversed)
        }
    }
}

impl Component for Slider {
    type State = Slip;
    type Error = ElementError;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::NonSmooth
    }

    fn stored_energy_j(&self, _slip: &Slip) -> Result<f64, ElementError> {
        Ok(0.0)
    }

    fn propose(
        &self,
        old: &Slip,
        dt_s: f64,
        inputs: &[f64],
        new: &mut Slip,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, ElementError> {
        self.propose_on_line(old, dt_s, inputs, &[0.0; 4], new, outputs, cx)
    }
}

impl LoadLine for Slider {
    fn line_ports(&self) -> &[usize] {
        &[0, 1]
    }

    fn propose_on_line(
        &self,
        old: &Slip,
        dt_s: f64,
        inputs: &[f64],
        compliance: &[f64],
        new: &mut Slip,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, ElementError> {
        cx.checkpoint()?;
        inputs_finite(inputs)?;
        inputs_finite(compliance)?;
        let kelvin = temperature(inputs[2])?;
        // An end whose row is zero does not answer the force: its flow is
        // given, and held over the step. If it was given another value in
        // the last step, it has jumped by the difference, and the step
        // starts with that much more slip than the last one ended with.
        let given = [
            compliance[0] == 0.0 && compliance[1] == 0.0,
            compliance[2] == 0.0 && compliance[3] == 0.0,
        ];
        let jumped = |end: usize| match (given[end], old.given_m_s[end]) {
            (true, Some(before)) => inputs[end] - before,
            _ => 0.0,
        };
        // The force acts as (e, -e) on the two ends, so a unit of it takes
        // (1, -1) C (1, -1)^T away from the slip.
        let c = compliance[0] - compliance[1] - compliance[2] + compliance[3];
        let released = inputs[0] - inputs[1];
        // With both ends given the slip is the given one from start to end.
        let start = if given[0] && given[1] {
            released
        } else {
            old.ended_m_s + jumped(0) - jumped(1)
        };
        let (force, mean, row) = self.on_line(start, released, c, old.stuck);
        // An end that answers follows the midpoint relation and a given
        // end keeps its flow, so the slip ends at twice its mean less its
        // start either way.
        *new = Slip {
            ended_m_s: 2.0 * mean - start,
            stuck: matches!(row, SliderRow::Held | SliderRow::Arrested),
            given_m_s: [
                given[0].then_some(inputs[0]),
                given[1].then_some(inputs[1]),
            ],
        };
        let loss_w = force * mean;
        outputs[0] = force;
        outputs[1] = -force;
        outputs[2] = -loss_w;
        let dissipated_j = dt_s * loss_w;
        // The two ends' powers are each as large as their own velocities
        // make them; their sum is rounded at that size, not at the loss's.
        let at_ends = inputs[0].abs()
            + inputs[1].abs()
            + force.abs() * compliance.iter().map(|entry| entry.abs()).sum::<f64>();
        Ok(Proposal::Step(Receipt {
            stored_change_j: 0.0,
            dissipated_j,
            numerical_loss_j: 0.0,
            entropy_generated_j_per_k: dissipated_j / kelvin,
            energy_tolerance_j: roundoff(dt_s * force.abs() * at_ends + dissipated_j),
            residual: 0.0,
            tolerance: 0.0,
        }))
    }
}
