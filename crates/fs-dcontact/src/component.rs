//! A normal contact between two bodies as a Component of the composition
//! kernel (`fs-ports`): [`ContactInterface`].
//!
//! The interface holds one admitted one-point [`Obstacle`] on the unit
//! opening map `[-1]` (as fs-couple's modal impact owner admits it), and its
//! state is the closure `x` [m], the displacement of end `a` towards end
//! `b`: `x' = v_a - v_b`. The penetration is `x - gap`. The contact pushes
//! the ends apart with the reaction `R >= 0`: its effort is `R` at `a` and
//! `-R` at `b`, so the power it receives is `R (v_a - v_b)`.
//!
//! # A step
//!
//! The closure moves by `dt` times the mean closure rate the ends give, and
//! the reaction is the law's discrete gradient over that move: the
//! [`OpeningContactStep`] secant of the power-law potential
//! `w K / (alpha + 1) [x - gap]_+^(alpha + 1)` between the start and end
//! closures (its derivative where they coincide). Its work over the move is
//! the potential's change, so a lossless contact stores what the ends give
//! it.
//!
//! On the load line ([`LoadLine`]) the neighbours' mean flows fall
//! linearly with the efforts, so the end closure is
//! `x(R) = x0 + dt (released_a - released_b) - S R` with
//! `S = dt (c_aa - c_ab - c_ba + c_bb)`, and the reaction solves
//! `R = law(x0, x(R))`. The secant of a convex potential is monotone in the
//! end closure, so `R - law(x(R))` is increasing and has one root in
//! `[0, law(x(0))]`. That bracket is bisected exactly as fs-couple's modal
//! impact owner bisects it (`bracket_reaction`): the root is accepted at the
//! first evaluated reaction whose residual is within
//! `force_absolute_tolerance_n + force_relative_tolerance max(|R|, |law|)`.
//! The search refuses ([`DContactError::Reaction`]) when the iterations run
//! out, when the bisection runs out of representable midpoints first, or
//! when the bracket passes `maximum_force_n`. The penetration the step ends
//! with is refused above `maximum_penetration_m`.
//!
//! Beside a neighbour that tells the end of its flow
//! ([`fs_ports::component::EndOfStep::Told`], such as a modal structure that
//! publishes its ends), an end reads that flow on a signal input
//! ([`ContactInterface::reading_ends`]). The contact's law needs only the
//! mean closure rate, so the end flows are checked to be finite and are not
//! otherwise used.
//!
//! # Not covered
//!
//! - Contact loss: an obstacle with a Hunt–Crossley loss is refused. The
//!   loss is heat, which the contract books only as the component's own
//!   stored thermal energy or through a heat port; this interface has
//!   neither yet.
//! - One contact point on one closure. Distributed obstacles
//!   ([`crate::AffineContactStep`]) are not reached by this interface.
//! - No events: a step that enters or leaves contact is not shortened to the
//!   instant of touch. The secant spans the whole move, as the modal impact
//!   owner's does.
//! - No checkpoint codec.

use core::num::NonZeroUsize;

use fs_exec::Cx;
use fs_math::det;
use fs_ports::component::{Class, Component, LoadLine, PortRole, PortSpec, Proposal, Receipt};
use fs_ports::{PortKind, StableId};

use crate::{DContactError, Obstacle, OpeningContactStep};

/// The budgets of one step's reaction search: the same fields and admission
/// as fs-couple's `ModalContactConfig`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ContactSolve {
    /// Bisection refinements, in `1..=128`. At most two bracket endpoint
    /// trials precede them.
    pub max_iterations: usize,
    /// Ceiling of the normal reaction [N]; it also bounds the bracket.
    pub maximum_force_n: f64,
    /// Ceiling of the penetration a step ends with [m].
    pub maximum_penetration_m: f64,
    /// Absolute tolerance of the reaction's residual [N].
    pub force_absolute_tolerance_n: f64,
    /// Relative tolerance of the reaction's residual, in `(0, 1)`, scaled by
    /// the larger of the trial reaction and the law's.
    pub force_relative_tolerance: f64,
}

impl ContactSolve {
    fn validate(&self) -> Result<(), DContactError> {
        if !(1..=128).contains(&self.max_iterations)
            || [
                self.maximum_force_n,
                self.maximum_penetration_m,
                self.force_absolute_tolerance_n,
                self.force_relative_tolerance,
            ]
            .iter()
            .any(|value| !value.is_finite() || *value <= 0.0)
            || self.force_relative_tolerance >= 1.0
        {
            return Err(DContactError::Parameter {
                what: "contact solve needs 1 to 128 iterations and positive finite force, penetration and residual budgets, the relative one below 1",
            });
        }
        Ok(())
    }

    /// Unfused, as fs-couple's `force_tolerance` rounds it.
    fn tolerance(&self, reaction: f64, expected: f64) -> f64 {
        self.force_absolute_tolerance_n
            + self.force_relative_tolerance * reaction.abs().max(expected.abs())
    }
}

/// The accepted state of a [`ContactInterface`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ContactState {
    /// Closure [m]: the displacement of end `a` towards end `b`; the
    /// penetration is the closure less the obstacle's gap.
    pub closure_m: f64,
    /// The reaction of the step that led here [N]; zero at the start.
    pub reaction_n: f64,
}

impl ContactState {
    /// A contact at `closure_m` that has not yet taken a step.
    #[must_use]
    pub const fn at(closure_m: f64) -> Self {
        Self {
            closure_m,
            reaction_n: 0.0,
        }
    }
}

/// A normal power-law contact between two mechanical ends. See the module
/// documentation.
#[derive(Debug, Clone)]
pub struct ContactInterface {
    law: Obstacle,
    solve: ContactSolve,
    ports: Vec<PortSpec>,
    /// Per end (`a`, `b`): the port of its end input, if it reads one.
    ends: [Option<usize>; 2],
}

fn port_error(_: impl core::fmt::Debug) -> DContactError {
    DContactError::Parameter {
        what: "contact interface port specification",
    }
}

fn finite(value: f64, what: &'static str) -> Result<f64, DContactError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(DContactError::Parameter { what })
    }
}

/// One trial of the reaction: the end closure, the law's reaction there,
/// the residual and its tolerance.
#[derive(Clone, Copy)]
struct Trial {
    reaction: f64,
    closure: f64,
    expected: f64,
    residual: f64,
    tolerance: f64,
}

impl ContactInterface {
    /// A contact with the one-point `law` (opening map `[-1]`, a finite gap,
    /// a nonnegative weight and stiffness, `alpha >= 1`, no internal loss
    /// and a provenance) and the search budgets `solve`.
    ///
    /// # Errors
    /// [`DContactError::Shape`] for an obstacle that is not one point on the
    /// unit opening map; [`DContactError::Parameter`] for nonfinite or
    /// nonphysical coefficients, an empty provenance, a nonzero internal
    /// loss (no heat port; see the module's "Not covered") or budgets
    /// [`ContactSolve`] does not admit.
    pub fn new(law: Obstacle, solve: ContactSolve) -> Result<Self, DContactError> {
        solve.validate()?;
        if law.n_points() != 1
            || law.collocation() != [-1.0]
            || law.gaps().len() != 1
            || law.weights().len() != 1
        {
            return Err(DContactError::Shape {
                what: "a contact interface holds one point on the unit opening map [-1]",
            });
        }
        // The secant's own admission: finite gap, weight, stiffness, alpha.
        OpeningContactStep::new(&law, 0.0)?;
        if law.provenance().trim().is_empty() {
            return Err(DContactError::Parameter {
                what: "a contact interface's law needs a provenance",
            });
        }
        // Admitted nonnegative above.
        if law.internal_loss() > 0.0 {
            return Err(DContactError::Parameter {
                what: "contact loss is heat, and this interface has no heat port: internal_loss must be zero",
            });
        }
        let name = |id: &str| StableId::new(id).map_err(port_error);
        let mechanical = |id: &str| {
            PortSpec::power(
                name(id)?,
                PortKind::MechanicalForceVelocity,
                PortRole::FlowIn,
                NonZeroUsize::MIN,
            )
            .map_err(port_error)
        };
        Ok(Self {
            law,
            solve,
            ports: vec![mechanical("a")?, mechanical("b")?],
            ends: [None; 2],
        })
    }

    /// The interface reading, at each end marked in `ends` (`a`, `b`), the
    /// flow its neighbours tell for the END of the step: a signal input
    /// `a.end` or `b.end` of width one, which [`LoadLine::line_ends`] names,
    /// so that the end joins a neighbour that answers
    /// [`fs_ports::component::EndOfStep::Told`]. The law does not use the
    /// end flows.
    ///
    /// # Errors
    /// No end marked, or end inputs already read.
    pub fn reading_ends(mut self, ends: [bool; 2]) -> Result<Self, DContactError> {
        if ends == [false; 2] || self.ends != [None; 2] {
            return Err(DContactError::Parameter {
                what: "contact interface end inputs: mark at least one end, once",
            });
        }
        for (end, id) in [(0, "a.end"), (1, "b.end")] {
            if ends[end] {
                let id = StableId::new(id).map_err(port_error)?;
                self.ends[end] = Some(self.ports.len());
                self.ports.push(PortSpec::signal_in(id, NonZeroUsize::MIN));
            }
        }
        Ok(self)
    }

    /// The admitted law.
    #[must_use]
    pub const fn law(&self) -> &Obstacle {
        &self.law
    }

    /// The search budgets.
    #[must_use]
    pub const fn solve(&self) -> ContactSolve {
        self.solve
    }

    /// How many ends read an end input.
    fn told(&self) -> usize {
        self.ends.iter().flatten().count()
    }

    /// The penetration at `closure_m` [m]: the closure less the gap, or zero
    /// when the ends are apart.
    #[must_use]
    pub fn penetration_m(&self, closure_m: f64) -> f64 {
        (closure_m - self.law.gaps()[0]).max(0.0)
    }

    /// The potential stored at `closure_m` [J],
    /// `w K / (alpha + 1) [x - gap]_+^(alpha + 1)`.
    fn potential_j(&self, closure_m: f64) -> f64 {
        let penetration = self.penetration_m(closure_m);
        if penetration <= 0.0 {
            return 0.0;
        }
        let beta = self.law.alpha() + 1.0;
        self.law.weights()[0] * self.law.stiffness() * det::pow(penetration, beta) / beta
    }

    /// The law's reaction for a move of the closure from `x0` to `x1` [N]:
    /// the secant, which is the conservative reaction (no loss is admitted).
    fn law_force(step: &OpeningContactStep, x1: f64) -> Result<f64, DContactError> {
        let (elastic, _) = step.coefficients(-x1)?;
        Ok(elastic.max(0.0))
    }

    /// The bracket search of fs-couple's `bracket_reaction`, on
    /// `x(R) = free - slope R`.
    fn settle(
        &self,
        step: &OpeningContactStep,
        free: f64,
        slope: f64,
        cx: &Cx<'_>,
    ) -> Result<Trial, DContactError> {
        let evaluate = |reaction: f64| -> Result<Trial, DContactError> {
            cx.checkpoint().map_err(|_| DContactError::Reaction {
                what: "cancelled",
                value: reaction,
                limit: 0.0,
            })?;
            // Unfused, as fs-couple's `solve_contact` rounds it.
            let closure = finite(free - slope * reaction, "contact trial closure")?;
            let expected = Self::law_force(step, closure)?;
            Ok(Trial {
                reaction,
                closure,
                expected,
                residual: finite(reaction - expected, "contact trial residual")?,
                tolerance: finite(
                    self.solve.tolerance(reaction, expected),
                    "contact trial tolerance",
                )?,
            })
        };
        let at_zero = evaluate(0.0)?;
        if at_zero.residual == 0.0 {
            return Ok(at_zero);
        }
        let mut lo = 0.0;
        // The secant of the convex potential is monotone in the closure, so
        // the reaction at no reaction bounds the root.
        let mut hi = (-at_zero.residual).min(self.solve.maximum_force_n);
        let at_hi = evaluate(hi)?;
        if at_hi.residual.abs() <= at_hi.tolerance {
            return Ok(at_hi);
        }
        if at_hi.residual < 0.0 {
            return Err(DContactError::Reaction {
                what: "the reaction passes maximum_force_n",
                value: hi - at_hi.residual,
                limit: self.solve.maximum_force_n,
            });
        }
        let mut last = at_hi;
        for _ in 1..=self.solve.max_iterations {
            let mid = f64::midpoint(lo, hi);
            let trial = evaluate(mid)?;
            if trial.residual.abs() <= trial.tolerance {
                return Ok(trial);
            }
            // Exact on purpose, as in `bracket_reaction`: the midpoint has
            // rounded onto an end of the bracket.
            #[allow(clippy::float_cmp)]
            let exhausted = mid == lo || mid == hi;
            if exhausted {
                return Err(DContactError::Reaction {
                    what: "the bisection ran out of representable midpoints",
                    value: trial.residual,
                    limit: trial.tolerance,
                });
            }
            if trial.residual < 0.0 {
                lo = mid;
            } else {
                hi = mid;
            }
            last = trial;
        }
        Err(DContactError::Reaction {
            what: "the reaction did not settle within max_iterations",
            value: last.residual,
            limit: last.tolerance,
        })
    }
}

impl Component for ContactInterface {
    type State = ContactState;
    type Error = DContactError;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::Monotone1d
    }

    fn stored_energy_j(&self, state: &ContactState) -> Result<f64, DContactError> {
        finite(state.closure_m, "contact closure")?;
        finite(self.potential_j(state.closure_m), "contact potential")
    }

    fn propose(
        &self,
        old: &ContactState,
        dt_s: f64,
        inputs: &[f64],
        new: &mut ContactState,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, DContactError> {
        // The line's square matrix and one row per end input, all zero.
        let zero = [0.0; 8];
        let compliance = &zero[..2 * (2 + self.told())];
        self.propose_on_line(old, dt_s, inputs, compliance, new, outputs, cx)
    }
}

impl LoadLine for ContactInterface {
    fn line_ports(&self) -> &[usize] {
        &[0, 1]
    }

    fn line_ends(&self) -> &[Option<usize>] {
        &self.ends
    }

    fn propose_on_line(
        &self,
        old: &ContactState,
        dt_s: f64,
        inputs: &[f64],
        compliance: &[f64],
        new: &mut ContactState,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, DContactError> {
        let (&[flow_a, flow_b, ref end_flows @ ..], &[c_aa, c_ab, c_ba, c_bb, ref end_rows @ ..]) =
            (inputs, compliance)
        else {
            return Err(DContactError::Shape {
                what: "contact interface inputs or line",
            });
        };
        if end_flows.len() != self.told() || end_rows.len() != 2 * self.told() {
            return Err(DContactError::Shape {
                what: "contact interface end inputs or their rows",
            });
        }
        let [force_a, force_b] = outputs else {
            return Err(DContactError::Shape {
                what: "contact interface outputs",
            });
        };
        if !(dt_s.is_finite() && dt_s > 0.0) {
            return Err(DContactError::Parameter {
                what: "contact step must be positive and finite",
            });
        }
        for &value in [flow_a, flow_b, c_aa, c_ab, c_ba, c_bb]
            .iter()
            .chain(end_flows)
            .chain(end_rows)
        {
            finite(value, "contact interface inputs or line")?;
        }
        let x0 = finite(old.closure_m, "contact closure")?;
        let step = OpeningContactStep::new(&self.law, -x0)?;
        // The end closure at no reaction, and what a unit of reaction takes
        // from it: the efforts are R at a and -R at b.
        let free = finite(x0 + dt_s * (flow_a - flow_b), "contact free closure")?;
        let slope = finite(dt_s * ((c_aa - c_ab) - (c_ba - c_bb)), "contact line slope")?;
        let trial = self.settle(&step, free, slope, cx)?;
        let penetration = self.penetration_m(trial.closure);
        if penetration > self.solve.maximum_penetration_m {
            return Err(DContactError::Reaction {
                what: "the penetration passes maximum_penetration_m",
                value: penetration,
                limit: self.solve.maximum_penetration_m,
            });
        }
        let moved = trial.closure - x0;
        // The secant's work over the move is the potential's change, with
        // no subtraction of nearly equal energies.
        let stored_change_j = finite(trial.expected * moved, "contact stored change")?;
        let work_j = finite(trial.reaction * moved, "contact work")?;
        let before_j = self.potential_j(x0);
        let after_j = self.potential_j(trial.closure);
        *new = ContactState {
            closure_m: trial.closure,
            reaction_n: trial.reaction,
        };
        *force_a = trial.reaction;
        *force_b = -trial.reaction;
        // The ends' powers are each as large as their own velocities make
        // them; the kernel's sum of them is rounded at that size.
        let at_ends = flow_a.abs()
            + flow_b.abs()
            + trial.reaction.abs() * (c_aa.abs() + c_ab.abs() + c_ba.abs() + c_bb.abs());
        let scale = dt_s * trial.reaction.abs() * at_ends + before_j + after_j;
        Ok(Proposal::Step(Receipt {
            stored_change_j,
            dissipated_j: 0.0,
            // What the settled reaction does beyond the law's secant.
            numerical_loss_j: work_j - stored_change_j,
            entropy_generated_j_per_k: 0.0,
            energy_tolerance_j: 64.0 * f64::EPSILON * scale + f64::MIN_POSITIVE,
            residual: trial.residual.abs(),
            tolerance: trial.tolerance,
        }))
    }
}
