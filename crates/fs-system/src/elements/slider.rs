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
//! prescribed end receives a mean too, the mean of what is prescribed
//! (`FlowSource` delivers the exact step mean of its waveform, and a
//! boundary is handed the step to do the same). Friction has to know the
//! slip a step starts and ends with, and for a prescribed end a mean does
//! not say where the side ends its step.
//!
//! Left to itself, the slider takes the mean for the end: a prescribed
//! side is HELD over its step, so its start, its mean and its end are
//! that one value, and it jumps from one step to the next. That is the
//! truth for a side at a constant speed and for one that does jump. It is
//! not for a side that accelerates: see "Not covered" below.
//!
//! A side that accelerates is TOLD its span. [`Slider::told`] gives an end
//! a signal input of width two, `span_a` or `span_b`: the velocities of
//! what moves that end at the start and at the end of the step, as the
//! mover returns them (`FlowSource::publishing` puts them out; a boundary
//! can prescribe them). Nothing is inferred. The sign is the mover's: the
//! flow the port itself receives is the mean with the OTHER sign, because
//! every connection hands a port minus the flow of what it is joined to.
//! Then the end's step is no longer the midpoint relation's by
//! `2 mean - from - end`, with `from` where its last step ended (the
//! span's start, in the first step), and the step's start slip is the
//! last end slip and that, for `a` less that for `b`. An end that is held
//! is the case `end = mean`, where this is the jump. The rows of the law
//! below are the same either way: only the start slip they are handed
//! differs.
//!
//! An end that ANSWERS the slider's force is taken to end its step by
//! the midpoint relation. Beside a mass the end of a step and twice its
//! mean less its start are one number, and every row below rests on that.
//! Beside a neighbour whose step is not the midpoint rule's they are not:
//! a modal structure, exact for a held force, returns the mean velocity of
//! its step and ends it elsewhere. Held at no mean slip such a neighbour is
//! never brought to rest; its velocity at the port changes sign from step
//! to step and does not decay (`tests/load_line.rs` shows it on two exact
//! oscillators: 0.5 m/s, 200 steps on). So the kernel does not seal one
//! beside the slider. A component says how the flow of a port ends its
//! step (`fs_ports::component::EndOfStep`; the default is the midpoint
//! relation), and one on a load line that needs the end names an input for
//! it (`LoadLine::line_ends`), which the kernel wires and gives a row of
//! the line. The slider names none: it is for neighbours on the midpoint
//! rule.
//!
//! A span on an end that answers the slider's force is refused. So is a
//! span of the port's sign, at its third step: a span on one side of zero
//! belongs to a flow on the other, but one step cannot say that a span on
//! the flow's side is wrong (a side that reverses twice within a step has
//! both its ends on one side and its mean on the other), and the slider
//! may not decide by a step's inputs at all, because the kernel also
//! proposes at inputs it has perturbed. It counts the steps taken in a
//! row with the span on the flow's side and refuses the step after two:
//! the port's sign given for the mover's, or a motion the step does not
//! resolve.
//!
//! The load line says which end is which. An end whose row of the line's
//! matrix is zero does not answer the slider's force at all: its flow is
//! given. The slider keeps the slip its last step ended with. Within a
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
//! | sliding | `-c mu_k N <= sg r < 0` | `r / c` | stuck, zero mean slip (turned) |
//! | sliding | `sg r < -c mu_k N` | `-mu_k N sg` | reversed |
//!
//! - The loss `e * mean` is not negative in any row, and the first law is
//!   exact on step means: what the two ends give up leaves through `heat`.
//! - `abs(e) <= mu_s N` stuck and `<= mu_k N` sliding.
//! - A mass thrown along a surface is arrested with the slider's slip at
//!   exactly zero and all of its kinetic energy given up as heat. (The
//!   mass's own velocity is then zero to the rounding its loop was closed
//!   to: the mass and the slider each keep their own record.)
//! - A pair held at zero mean slip is stuck from whichever state it came.
//!   A turned step holds the mean with a force inside `mu_k N`, so inside
//!   `mu_s N`, as a held one does, and the step after it reads the stuck
//!   rows: carried, held, or broken away at `mu_s N`. The two records
//!   differ by their rounding, alternating in sign (2.2e-14 m/s on the
//!   friction belt of fs-cli), so where the slip is within that rounding
//!   of zero the slider's start slip and the released slip can have
//!   opposite signs, step after step. Each such step holds the mean at
//!   zero: turned from the sliding state, held from the stuck one. While
//!   turned was recorded as sliding, every such step stayed turned, and a
//!   belt caught on the boundary of a step was never stuck again (none of
//!   49 steps there).
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
//! # Not covered
//!
//! A prescribed side that ACCELERATES and is not told its span. Taken as
//! held, it is off by half a step of its acceleration at every step, and
//! no smaller step mends what follows (the numbers are in the crate's
//! contract and asserted in `tests/load_line.rs`):
//!
//! - A side that leaves rest in the first step has a first mean that is
//!   not the mass's velocity, a slip the slider is not told. A pair that
//!   needs 8 N of the 9.8 N it holds breaks away within a few steps, at
//!   every step length (which step hangs on a rounding residue).
//! - A side that changes speed after a step at rest is followed as a
//!   staircase: with the right force, one mean behind, and with the catch
//!   of every step booked as heat, `m a^2 dt^2 / 2` a step, where a stuck
//!   pair loses nothing.
//!
//! The slider cannot refuse such a side: it cannot tell a side that
//! accelerates from one that jumps. A prescribed side that does not stay
//! at one speed is to be told its span.
//!
//! A change of the side's acceleration INSIDE a step, with the span told.
//! The slider is told the step's two ends and its mean, no more. In the
//! step of the change, the force that would end the step on the side has
//! a negative loss where the side speeds up more, and the slider holds
//! on the mean instead; the mass ends that step behind and the next step
//! makes it up, with a force above the needed one by `m` times a quarter
//! of the change of acceleration at most, for one step. A pair that is
//! within that of its static limit breaks away there though the side
//! never asks for more than it holds. Changes in consecutive steps add
//! up, each by its own quarter at most as far as measured. A change of
//! acceleration at the end of a step is followed exactly.
//! Tabulated `FlowSource`s advertise their knots through `next_event_s`,
//! so the system resolves these changes at step boundaries. A raw
//! prescribed boundary must arrange that step cap itself.
//!
//! A starting slip that contradicts the ports. Nothing can check it. A
//! pair declared stuck on a side that already moves keeps a record that
//! is off by that slip and alternates in sign, and its run can end not
//! converged.
//!
//! The slider's record of an end that answers is only as good as the
//! midpoint relation is for what is attached there.

use core::num::NonZeroUsize;

use fs_exec::Cx;
use fs_ports::component::{Class, Component, LoadLine, PortRole, PortSpec, Proposal, Receipt};
use fs_ports::{PortKind, StableId};

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
    ports: Vec<PortSpec>,
    /// Per end that is told its span, where the span starts among the
    /// inputs.
    spans: [Option<usize>; 2],
}

/// An end of a [`Slider`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SliderEnd {
    /// Port `a`.
    A,
    /// Port `b`.
    B,
}

/// What a [`Slider`] keeps between steps: the slip its last step ended
/// with [m/s], whether that step ended stuck, and, for each end whose flow
/// was given in that step, the flow that end ended the step with.
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
    /// Per end whose flow was given in the last step, the flow it ended
    /// that step with, as the port receives it: the flow it was given if
    /// it is held, minus the end of its span if it is told one. The next
    /// step takes it up there.
    given_m_s: [Option<f64>; 2],
    /// Per end that is told its span: how many steps in a row, up to the
    /// last, had the span on one side of zero and the flow the port
    /// receives on that same side.
    mirrored: [u8; 2],
}

impl Slip {
    /// Stuck: the two ends move together.
    #[must_use]
    pub fn stuck() -> Self {
        Self {
            ended_m_s: 0.0,
            stuck: true,
            given_m_s: [None; 2],
            mirrored: [0; 2],
        }
    }

    /// Sliding with this slip (flow at `a` less flow at `b`).
    #[must_use]
    pub fn sliding(slip_m_s: f64) -> Self {
        Self {
            ended_m_s: slip_m_s,
            stuck: false,
            given_m_s: [None; 2],
            mirrored: [0; 2],
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
    /// Sliding, reversed with zero mean slip: held from sliding, with a
    /// force inside the sliding limit, and recorded stuck as a held step
    /// is.
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
            ports: vec![
                port("a", kind, PortRole::FlowIn),
                port("b", kind, PortRole::FlowIn),
                heat("heat", PortRole::EffortIn),
            ],
            spans: [None; 2],
        })
    }

    /// The same slider with a signal input `span_a` or `span_b` of width
    /// two for `end`: the velocities [m/s] of what moves that end at the
    /// start and at the end of the step, as the mover returns them (a
    /// `FlowSource::publishing` puts them out). The flow the port itself
    /// receives is the mean, and with the other sign: every connection
    /// hands a port minus the flow of what it is joined to. For an end
    /// that is prescribed and does not stay at one speed; see the module's
    /// "What a port tells the slider". Asking twice adds it once.
    #[must_use]
    pub fn told(mut self, end: SliderEnd) -> Self {
        let (index, name) = match end {
            SliderEnd::A => (0, "span_a"),
            SliderEnd::B => (1, "span_b"),
        };
        if self.spans[index].is_none() {
            // After `a`, `b`, `heat` and the span that is there already.
            let at = 3 + 2 * self.spans.iter().flatten().count();
            let id = StableId::new(name).expect("element port names are valid stable ids");
            self.ports
                .push(PortSpec::signal_in(id, NonZeroUsize::MIN.saturating_add(1)));
            self.spans[index] = Some(at);
        }
        self
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

/// How many steps in a row a span may lie on the side of zero of its port's
/// flow before the next step is refused.
const MIRRORED_STEPS: u8 = 2;

impl Slider {
    /// The span of end `end` (0 for `a`) in the sign the port receives its
    /// flow in: where it starts and ends the step. `from` and `to` are the
    /// two velocities of the signal, `given` whether the end's row of the
    /// line is zero, `mirrored` the end's count of the state.
    ///
    /// # Errors
    /// A span on an end that answers the slider's force, or one that had
    /// the sign of the port's flow in the last [`MIRRORED_STEPS`] steps.
    fn told_span(
        end: usize,
        given: bool,
        mirrored: u8,
        from: f64,
        to: f64,
    ) -> Result<(f64, f64), ElementError> {
        let name = if end == 0 { "a" } else { "b" };
        if !given {
            return Err(ElementError::Span {
                end: name,
                what: "the end answers the slider's force, so where it ends its step follows \
                       from its mean; a span is for an end that is prescribed",
            });
        }
        // Decided by the steps that were taken, never by this step's
        // inputs: the kernel also proposes at inputs it has perturbed to
        // difference the slider, and a perturbed flow of a slow side
        // changes sign.
        if mirrored >= MIRRORED_STEPS {
            return Err(ElementError::Span {
                end: name,
                what: "in its last two steps it lay on the side of zero of the flow the port \
                       received; a span is minus that flow at the step's two ends, so it is \
                       given in the port's sign, or the step does not resolve the motion",
            });
        }
        Ok((-from, -to))
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
        // Per given end: the flow it ends the step with, and by how much
        // its step is not the midpoint relation's, `2 mean - from - end`
        // with `from` where its last step ended. Held, it ends with the
        // flow it is given, and that is the jump from the last step's.
        // Told its span, it ends where the span ends.
        let mut ends = [None; 2];
        let mut moved = [0.0_f64; 2];
        let mut mirrored = [0_u8; 2];
        for end in 0..2 {
            if let Some(at) = self.spans[end] {
                let (from, to) = (inputs[at], inputs[at + 1]);
                let (starts, finishes) =
                    Self::told_span(end, given[end], old.mirrored[end], from, to)?;
                // The port receives minus what its mover returns: where
                // the span lies on one side of zero, the flow lies on the
                // other. One step in which it does not can be a side that
                // reverses twice within it; steps in a row are counted.
                if from * to > 0.0 && from * inputs[end] > 0.0 {
                    mirrored[end] = old.mirrored[end] + 1;
                }
                let from = old.given_m_s[end].unwrap_or(starts);
                moved[end] = 2.0 * inputs[end] - from - finishes;
                ends[end] = Some(finishes);
            } else if given[end] {
                if let Some(before) = old.given_m_s[end] {
                    moved[end] = inputs[end] - before;
                }
                ends[end] = Some(inputs[end]);
            }
        }
        // The force acts as (e, -e) on the two ends, so a unit of it takes
        // (1, -1) C (1, -1)^T away from the slip.
        let c = compliance[0] - compliance[1] - compliance[2] + compliance[3];
        let released = inputs[0] - inputs[1];
        // With both ends given the slip ends where the two ends end (held:
        // at the given slip, from start to end). Else it goes on from
        // where the last step left it.
        let start = match ends {
            [Some(a), Some(b)] => 2.0 * released - (a - b),
            _ => old.ended_m_s + moved[0] - moved[1],
        };
        let (force, mean, row) = self.on_line(start, released, c, old.stuck);
        // An end that answers follows the midpoint relation, and what a
        // given end does besides is in the start slip, so the slip ends
        // at twice its mean less that start.
        *new = Slip {
            ended_m_s: 2.0 * mean - start,
            // Held at zero mean slip, from either state, is stuck.
            stuck: matches!(
                row,
                SliderRow::Held | SliderRow::Arrested | SliderRow::Turned
            ),
            given_m_s: ends,
            mirrored,
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
