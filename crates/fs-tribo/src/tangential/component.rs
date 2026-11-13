//! Friction laws on the composition contract: a sliding interface under a
//! prescribed slip rate, and a friction interface between two bodies.
//!
//! # Under a prescribed slip
//!
//! A [`TangentialLaw`] owns the traction, the solved slip, the regime, any
//! reversible storage and the dissipation it exports. This adapter drives
//! one with a released slip increment of `rate * dt` at zero compliance
//! (prescribed motion) and with the interface temperature its heat port
//! receives, and hands the dissipation to that port. It adds no heat
//! fraction, no heat capacity and no flash temperature. Where the heat goes
//! is decided by what the port is connected to: a thermal node between two
//! bodies divides it by their conductances and their temperatures.
//!
//! Port `heat` receives the interface temperature [K] and returns the heat
//! rate into the interface [W], which is minus the dissipation rate. Port
//! `traction` is a signal output of two values, the traction entering the
//! interface [N] along the two tangent axes.
//!
//! The law at a temperature. Friction that depends on temperature enters as
//! the caller's source `law_at(T)`, which returns the law to use at the
//! interface temperature of a step; a law that is itself bounded in
//! temperature also receives that temperature as its thermal boundary. The
//! source is called at every evaluation of a step and must be a pure
//! function.
//!
//! Work. No mechanical port is offered here. A dry law is set-valued in the
//! slip rate, and a kernel that differences smooth port maps cannot solve
//! for the slip through a port ([`FrictionInterface`], below, has the ports,
//! on the kernel's non-smooth tier). The motion is prescribed, so the work it does
//! against the traction is the component's own account, as a prescribed
//! force's is: the stored energy is the law's reversible energy less the
//! work received so far, and a step's receipt books the law's closure
//! residual `storage + dissipation - work` as numerical loss.
//!
//! Entropy generated: the dissipation over the interface temperature.
//!
//! Class: `Smooth` in its one input, as far as the source and the law are
//! smooth in temperature at a prescribed slip. No tangent is offered; a
//! kernel differences the map.
//!
//! Refused: an interface temperature or a step that is not positive and
//! finite; whatever the source or the law refuses, unchanged; a step whose
//! work balance is off by more than the declared tolerance, or whose
//! reported stored change disagrees with the law's stored energies by more
//! than it (the admission of the hand-written friction-heat coupling's
//! `tangential_heat_candidate`); a source whose law gives the held state
//! another stored energy than the state was left with, by any amount. A
//! reversible law whose stiffness depends on temperature needs a
//! thermoelastic account that this adapter does not have.
//!
//! One call does a bounded amount of work (the law's own root budget), and
//! [`TriboError`] has no cancellation value, so `cx` is not polled here.
//!
//! # Between two bodies
//!
//! [`FrictionInterface`] puts the same kind of law between two mechanical
//! ports. `a` and `b` receive the velocities [m/s] of the two sides and
//! return the friction force [N] and its opposite; the slip rate is the flow
//! at `a` less the flow at `b`. `heat` receives the interface temperature
//! [K] and returns the heat rate into the interface [W], minus the loss. The
//! law comes from `law_at(T)`, as above.
//!
//! A dry law holds any force up to its limit at zero slip, so its force is
//! not a function of the slip. The component is therefore added on the
//! kernel's non-smooth tier (`fs_ports::component::LoadLine`, `add_on_line`
//! in `fs-system`), which hands it the line its neighbours present:
//! `mean slip = released - c force`.
//!
//! The bodies of a kernel step on the midpoint rule: a port carries the MEAN
//! velocity of a step, and the slip a step ends with is twice its mean less
//! the slip it started with. With `s0` the slip a step starts with, `K` the
//! law's kinetic traction, and `along` the force that would end the step at
//! zero slip, `(released - s0 / 2) / c` taken along `s0`:
//!
//! | the step starts | condition | what the law is asked | the step |
//! |---|---|---|---|
//! | at rest | the load is held at both ends of the step, or at neither | its own return map on the line, for the slip that ends the step at none: released increment `dt (released - s0 / 2)`, compliance `dt c` | holds and ends with no slip, or slides at the traction the law solves |
//! | at rest | the load is held at the step's start and not at its end | its static verdict at loads between the two | is not taken: the component asks for the step that ends where the load reaches what the law holds |
//! | sliding | `along > K` | its traction at a prescribed mean slip; the mean solves `m = released - c K(m)` | slides on the same way |
//! | sliding | `along = K`, to the stop's band | its traction at the mean `s0 / 2` | ends at the stop under that traction, at rest; what is left of the slip, inside the tolerance, is kept for the next step |
//! | sliding | `along < K` | its traction at the mean `s0 / 2` | is not taken: the slip would pass zero inside it, and the component asks for the step that ends at the stop (`Proposal::Shorten`) |
//!
//! - The law's static map is asked only of a pair at rest. Asked of a
//!   sliding pair it would hold the MEAN slip of the step at zero: the slip
//!   would end at `-s0`, and the heat of the stopping distance would be
//!   lost. (The hand-written `tribothermal::dynamics::collinear` of
//!   `fs-couple` says the same of its static probe.)
//! - A step never carries the slip through zero. Arrest and reversal go
//!   through the kernel's event path: the step is cut at the stop, and the
//!   next one starts at rest, where the static map decides whether the pair
//!   holds or slides the other way.
//! - Rest is a mark in the state. What is near enough to rest is decided
//!   by momentum: `dt / (2 c)` is the inertia the slip moves against, and
//!   [`InterfaceLoad::impulse_tolerance_n_s`] is the relative momentum
//!   below which there is none. A starting slip below it is rest, and the
//!   stop's band is the law's own force tolerance plus that tolerance over
//!   `dt`. (It is the impulse tolerance of the hand-written module. With
//!   the law's force tolerance alone, 1e-12, a run of 47 steps ended its
//!   cut step with a slip of 7e-15 m/s that the band did not take, and the
//!   step after it could not be located.)
//! - A stop is taken under the law's own kinetic traction. It leaves a
//!   slip inside the tolerance, which the state keeps; the first step at
//!   rest holds the pair with the force that ends that step at no slip,
//!   and from then on the slip is exactly zero.
//! - A break-away is located as a stop is. The static map decides for a
//!   whole step from the step's mean, so under a load that rises through
//!   a step it would let go at the start of the first step whose middle is
//!   beyond the limit, up to half a step early or late. The line gives the
//!   mean load on the slip in every step, `(released - s0) / c`, and the
//!   state keeps the last step's with that step's length. The load at this
//!   step's two ends is read off the straight line through the two means,
//!   which is exact for a load that is linear in time whatever the two
//!   lengths. Where the law holds the load at the start and not at the
//!   end, the step is cut where its verdict changes; the cut step is held,
//!   and the next one starts at the limit and lets go from its start.
//!   What the law holds is asked of the law, by bisection on its own
//!   static verdict inside its root budget: no capacity is restated here,
//!   so a law whose holding depends on its state or its load is asked the
//!   same way. How near its limit it is asked depends on what its verdicts
//!   cost it, by its own count of iterations. Where they cost nothing (a
//!   closed form, as Coulomb's) the instant is narrowed to rounding. A
//!   law that searches for its verdict is asked no nearer than the impulse
//!   tolerance needs: just beyond what it holds its sliding root is next
//!   to zero, its search for that root takes more with every halving, and
//!   narrowed to rounding it was refused by the law's own budget (found
//!   on a bowed string). An instant within the impulse tolerance of either end of
//!   the step is that end (half the load's rate times the time squared is
//!   the momentum of letting go that much early or late). Under a load
//!   that curves the first cut lands near the instant and the cut step is
//!   located again from its own line, until that holds.
//!
//!   The line through two means was chosen over carrying the load's end
//!   from step to step as `end = 2 mean - start`: where the load curves,
//!   the error of that record changes sign with every step and never
//!   decays.
//! - Every row is wide, or has the slope of its neighbours. A kernel that
//!   differences this component moves a flow by `sqrt(eps)` of its scale,
//!   far more than the stop's band is wide. The force that would end a
//!   step at exactly no slip is stiff in the released flow (slope `1 / c`),
//!   and returned only inside the band it would be differenced as the
//!   sliding rows' slope of zero: a step that needs one more update at
//!   such an iterate does not converge. The component as first landed did
//!   that, and a block on a belt at 0.137 m/s was refused at its third
//!   catch. The stop therefore has the sliding rows' traction, and the
//!   stiff force is returned only at rest, where its row is as wide as the
//!   static limit.
//! - The heat is the work of the two ends: the force times the mean slip
//!   that crossed, less what the law stores. It is the law's own
//!   dissipation to the tolerance of the row (the law is asked at the slip
//!   a row aims at, and the line settles a tolerance away). The law's own
//!   accounts for the step it was asked for are admitted as above, within
//!   `work_tolerance_j`.
//!
//! Where the stop is found. A step that ends at zero slip has the mean
//! `s0 / 2` whatever its length, so the kinetic traction of the stop is
//! known before the stop is. The stop is then located from the line of the
//! step that was asked for, at `dt (abs(s0) / 2) / (c K - (released - s0))`
//! with the last term taken along `s0`, as if what the interface is joined
//! to were bodies under held forces: for those `released - s0` and `c` grow
//! in proportion to the step. For those it is exact: one cut, and the cut
//! step ends at rest. Beside anything else (a spring) the located instant is
//! only as good as the step it is located from. The cut step ends near the
//! stop. Beyond it, it is cut again, from a line that is now close. Short of
//! it, it is taken, sliding, and the stop is left to the next step, which is
//! asked at full length again. Measured, in
//! `crates/fs-couple/tests/tribothermal/dynamics/collinear/graph.rs`: beside
//! a linear spring two or three cuts in the step that holds the stop; beside
//! springs with a cubic term in steps of 0.07 s, up to five cuts in that
//! step, or, with two bodies, three further steps of 4e-7 s and less. The
//! stop is reached in every case.
//!
//! What it needs of the kernel. A perturbed proposal that asks to shorten
//! must not shorten the step, and the component must first be asked at the
//! released flows of a step (`fs-system`, "Requests to shorten inside a
//! loop"). Without the first a stop was cut a second time, short of itself;
//! without the second a pair that starts sliding with no starting flows
//! given to its ports was cut without end.
//!
//! A prescribed end (its row of the line is zero) is held over the step,
//! and a change of its flow from one step to the next is a jump: the step
//! starts with that much more slip, and the pair slides until the other
//! side is caught. With both ends prescribed there is no line, and the
//! component is the prescribed interface above with its force on two ports.
//!
//! Class `NonSmooth`. Entropy generated: the heat over the interface
//! temperature. `cx` is not polled, as above.
//!
//! Refused, besides what the prescribed interface refuses: a law that holds
//! a slip it is prescribed, or that keeps slip where it holds a pair at rest
//! (`friction_interface_law_holds_a_prescribed_slip`,
//! `friction_interface_static_map_kept_slip`). The rows above are for laws
//! that are set-valued at rest and slide at any prescribed slip, the dry
//! laws; a law with a compliance of its own is a function of the slip and
//! needs no line. Also refused: a traction against the slip it is asked at,
//! and a sliding root not found within the law's budget
//! (`TriboError::TangentialSolve`).
//!
//! Not covered:
//!
//! - One tangent axis. The law's second axis is given no slip.
//! - A prescribed side that ACCELERATES. Held over each step it is a
//!   staircase of jumps: a pair that static friction would carry is
//!   dragged under the kinetic traction after every jump and the step is
//!   cut where it is caught, a fraction `m a / K` of the step before. So
//!   each step asked for is followed by a train of steps that shrink by
//!   that fraction, down to the impulse tolerance: more than ten steps
//!   for one, and heat where a carried pair makes none (a test of the
//!   graph file above has the numbers). The kernel's ideal slider can be
//!   told the span of such a side; this component cannot.
//! - A break-away in a run's first step, which has no step before it to
//!   draw the line through, and one under a load that jumps (it is beyond
//!   the limit from the step's start). Both are decided by the step's
//!   mean, and the pair lets go from the step's start.
//! - A traction that falls with speed can meet the line at more than one
//!   speed. While sliding, the root is bracketed between half the starting
//!   slip and the released speed, and a sign bracket is not a proof that it
//!   is the only one. From rest the law's own map chooses.
//! - The normal load is a constant of the component.
//!
//! # Beside a neighbour that tells the end of its flow
//!
//! The table above rests on the midpoint rule: a step's slip ends at twice
//! its mean less its start. A neighbour whose step is exact for a held force
//! (a modal string, `ModalStructure` of `fs-nlmodal`) ends its step
//! elsewhere and tells where, on a signal output
//! (`fs_ports::component::EndOfStep::Told`). Beside it the midpoint record
//! is not the neighbour's own: on a bowed string (48 kHz, 64 modes) a pair
//! held at no mean slip ended its steps with a slip that changed sign every
//! step and did not decay, 5e-3 m/s from the start and 0.045 m/s through a
//! stick.
//!
//! [`FrictionInterface::reading_ends`] gives an end a signal input, `a.end`
//! or `b.end`, that [`LoadLine::line_ends`] names. The kernel wires the
//! neighbours' told ends to it and hands the component one more row of the
//! line per such end: `end = released_end - row . (F, -F)`. An end that
//! reads none must be given (its row of the line zero), and it ends the step
//! at the flow it is given. With `r_e` the released end slip and `c_e` what
//! a unit of the force takes from it, the slip ends the step at
//! `r_e - c_e F`; the mean slip is `r - c F` (`r` the released slip, `c` the
//! line's compliance, as above), and the two ends do the work
//! `dt F (r - c F)`.
//!
//! | the step starts | condition | what the law is asked | the step |
//! |---|---|---|---|
//! | at rest | the load is held at both ends of the step, or at neither | its own return map on the END line: released increment `dt r_e`, compliance `dt c_e` | holds and ends with no slip, or lets go at the traction the law solves and ends at `r_e - c_e F` |
//! | at rest | the load is held at the step's start and not at its end | its static verdict at loads between the two | is not taken: the component asks for the step that ends where the load reaches what the law holds |
//! | sliding | `along > K` | its traction at a prescribed mean slip; the mean solves `m = r - c K(m)` | slides on, and ends at `r_e - c_e K` |
//! | sliding | the stop step's mean is not along the slip, and the law does not hold `along` | its own return map on the END line, for its verdict at `along`; then its traction at a prescribed mean slip, from a mean at rest | slides on, and ends at `r_e - c_e K` |
//! | sliding | the stop step's mean is not along the slip, the law holds `along`, and `along < K(s0)` | its own return map on the END line, for its verdict at `along`; then its traction `K(s0)` at the starting slip | is not taken: the component asks for the step that ends where the end line stops under `K(s0)`, which only places the cut |
//! | sliding | `along = K`, to the stop's band | its traction at the stop step's mean | ends at the stop under that traction, at rest; the leftover `r_e - c_e K`, inside the tolerance, is kept |
//! | sliding | `along < K` | its traction at the stop step's mean | is not taken: the component asks for the step that ends at the stop |
//!
//! Here `way` is the sign of the starting slip `s0`, `along = way r_e / c_e`
//! the force that ends the step at no slip, `way r - c along` the stop
//! step's mean speed, and `K` the law's kinetic traction at that speed.
//!
//! - A pair at rest is held at no END slip, as `fs-couple`'s hand-written
//!   bowed string holds it. Beside a string the step's mean slip is then not
//!   zero, and its work `dt F (r - c F)` has either sign. It is booked as
//!   `numerical_loss_j`, signed: a defect of the step, with no heat and no
//!   entropy, which the ledger shows on the component.
//! - A sliding step asks the law at its MEAN slip, as above, so its heat is
//!   the law's own dissipation. The hand-written module asks its law at the
//!   end slip, so the two differ in sliding steps; a graph is compared with
//!   the module by its verdicts and figures (Helmholtz or not, the slip
//!   fraction, the period), not to rounding.
//! - A step that lets go from rest has the traction the law solves on the
//!   end line, and its heat is the work of the two ends less what the law
//!   stores, not the law's dissipation at its end slip.
//! - Rest is a starting slip below the impulse tolerance in
//!   `dt abs(s0) / c_e`; for a mass `dt / c_e` is its inertia, as
//!   `dt / (2 c)` is above. A break-away is located as above, from the load
//!   `(r_e - s0) / c_e` of each step (for a mass the push on the slip, as
//!   `(released - s0) / c` is above), with the law's verdicts asked on the
//!   end line. A stop is located on the end line, at
//!   `dt abs(s0) / (c_e K - way (r_e - s0))`, which is exact for bodies under
//!   held forces. A stop's leftover needs nothing more: the released end is
//!   the neighbour's own, which holds it.
//! - Refused, by field: an end that reads no end input and is not given,
//!   beside one that does (`friction_interface_untold_end`: a body that
//!   ends its step by the midpoint relation would need a start of its own,
//!   which the state does not keep); an end line that takes nothing from the
//!   slip where the mean line does (`friction_interface_end_compliance`); a
//!   sliding pair whose stop step's mean is not along its slip, while its law
//!   holds `along` and `along` is not below `K(s0)` by the stop's band, or
//!   whose law's traction on the end line would stop the mean
//!   (`friction_interface_stop_against_the_slip`); a step that lets go from
//!   rest with a heat below minus the work tolerance
//!   (`friction_interface_let_go_heat`).

use core::fmt;
use core::marker::PhantomData;
use core::num::NonZeroUsize;

use fs_exec::Cx;
use fs_ports::component::{Class, Component, LoadLine, PortRole, PortSpec, Proposal, Receipt};
use fs_ports::{PortKind, StableId};

use super::{TangentialLaw, TangentialSolveConfig, TangentialStep, TangentialStepInput};
use crate::{FrictionRegime, TriboError};

/// The load and motion a [`SlidingInterface`] is held to.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PrescribedSlip {
    /// Compressive normal load [N], non-negative.
    pub normal_n: f64,
    /// Relative slip rate along the two tangent axes [m/s].
    pub slip_rate_m_s: [f64; 2],
    /// Nominal pressure [Pa], for a law that is bounded in it.
    pub nominal_pressure_pa: Option<f64>,
    /// The law's numerical and effort budgets.
    pub solve: TangentialSolveConfig,
    /// Largest admitted defect of a step's work balance [J].
    pub work_tolerance_j: f64,
}

/// State of a [`SlidingInterface`]: the law's own checkpoint and the two
/// terms of the stored energy.
#[derive(Debug, Clone, PartialEq)]
pub struct SlidingState<S> {
    tangential: S,
    /// The law's reversible energy at `tangential` [J].
    reversible_j: f64,
    /// Work the prescribed motion has done on the interface so far [J].
    work_j: f64,
}

impl<S> SlidingState<S> {
    /// The law's checkpoint: its slip, its thermal boundary, its layers.
    #[must_use]
    pub const fn tangential(&self) -> &S {
        &self.tangential
    }

    /// Work the prescribed motion has done on the interface so far [J].
    #[must_use]
    pub const fn work_j(&self) -> f64 {
        self.work_j
    }
}

/// A sliding interface whose slip rate and normal load are prescribed, with
/// the law `law_at` returns for the interface temperature of each step.
pub struct SlidingInterface<L, F> {
    law_at: F,
    drive: PrescribedSlip,
    ports: [PortSpec; 2],
    law: PhantomData<fn() -> L>,
}

impl<L, F> fmt::Debug for SlidingInterface<L, F> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("SlidingInterface")
            .field("drive", &self.drive)
            .finish_non_exhaustive()
    }
}

fn refused(field: &'static str) -> TriboError {
    TriboError::InvalidInput { field }
}

fn finite(value: f64, field: &'static str) -> Result<f64, TriboError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(refused(field))
    }
}

/// One evaluated step: the candidate state, the heat rate into the
/// interface [W], the traction [N] and the receipt.
type Stepped<S> = (SlidingState<S>, f64, [f64; 2], Receipt);

impl<L, F> SlidingInterface<L, F>
where
    L: TangentialLaw,
    F: Fn(f64) -> Result<L, TriboError>,
{
    /// An interface held to `drive`, with the law `law_at` gives at an
    /// interface temperature [K].
    ///
    /// # Errors
    /// A normal load or nominal pressure that is negative or not finite, a
    /// slip rate that is not finite, a work tolerance that is negative or
    /// not finite, or solve budgets the law's owner refuses.
    pub fn new(law_at: F, drive: PrescribedSlip) -> Result<Self, TriboError> {
        drive.solve.validate()?;
        crate::nonnegative_finite(drive.normal_n, "normal_force")?;
        if let Some(pressure_pa) = drive.nominal_pressure_pa {
            crate::nonnegative_finite(pressure_pa, "contact_pressure_pa")?;
        }
        for rate in drive.slip_rate_m_s {
            finite(rate, "slip_rate_m_s")?;
        }
        crate::nonnegative_finite(drive.work_tolerance_j, "work_tolerance_j")?;
        let name = |id: &str| StableId::new(id).map_err(|_| refused("port_identifier"));
        Ok(Self {
            law_at,
            drive,
            ports: [
                PortSpec::heat(name("heat")?, PortRole::EffortIn)
                    .map_err(|_| refused("heat_port"))?,
                PortSpec::signal_out(
                    name("traction")?,
                    NonZeroUsize::new(2).ok_or(refused("traction_port"))?,
                ),
            ],
            law: PhantomData,
        })
    }

    /// The load and motion this interface is held to.
    #[must_use]
    pub const fn drive(&self) -> &PrescribedSlip {
        &self.drive
    }

    /// The unloaded state of the law at `interface_temperature_k`, with no
    /// work received: the state to hand to a kernel with this component.
    ///
    /// # Errors
    /// A temperature that is not positive and finite, the source's refusal,
    /// or a reversible energy that is not finite.
    pub fn start(&self, interface_temperature_k: f64) -> Result<SlidingState<L::State>, TriboError> {
        let law = self.law_at_k(interface_temperature_k)?;
        let tangential = law.initial_state();
        let reversible_j = finite(law.stored_energy_j(&tangential), "stored_energy_j")?;
        Ok(SlidingState {
            tangential,
            reversible_j,
            work_j: 0.0,
        })
    }

    fn law_at_k(&self, interface_temperature_k: f64) -> Result<L, TriboError> {
        if !(interface_temperature_k.is_finite() && interface_temperature_k > 0.0) {
            return Err(refused("interface_temperature_k"));
        }
        (self.law_at)(interface_temperature_k)
    }

    /// The law's step from `old` over `dt_s` at the interface temperature
    /// `inputs[0]`, admitted as the friction-heat coupling admits one.
    fn step(
        &self,
        old: &SlidingState<L::State>,
        dt_s: f64,
        inputs: &[f64],
    ) -> Result<Stepped<L::State>, TriboError> {
        let &[kelvin] = inputs else {
            return Err(refused("sliding_interface_inputs"));
        };
        if !(dt_s.is_finite() && dt_s > 0.0) {
            return Err(refused("dt_s"));
        }
        let law = self.law_at_k(kelvin)?;
        let tolerance_j = self.drive.work_tolerance_j;
        // The held state must be worth to this law exactly what it was left
        // worth: otherwise changing the law has created or destroyed
        // energy. No tolerance is given here. The receipt's stored change
        // is formed from this value and the kernel re-measures it from the
        // one the state holds, so any difference the check let through
        // would be a broken audit instead of this refusal.
        let before_j = law.stored_energy_j(&old.tangential);
        #[allow(clippy::float_cmp)]
        if before_j != old.reversible_j {
            return Err(refused("sliding_interface_law_changed_stored_energy"));
        }
        let step = law.advance(
            &old.tangential,
            &TangentialStepInput {
                normal_n: self.drive.normal_n,
                slip_increment_m: self.drive.slip_rate_m_s.map(|rate| rate * dt_s),
                dt_s,
                contact_temperature_k: Some(kelvin),
                nominal_pressure_pa: self.drive.nominal_pressure_pa,
                // Prescribed motion: the bodies do not give under traction.
                compliance_m_per_n: 0.0,
                solve: self.drive.solve,
            },
        )?;
        // Work is traction times the solved slip, never the released motion.
        let work_j = step.traction_n[0] * step.slip_increment_m[0]
            + step.traction_n[1] * step.slip_increment_m[1];
        let after_j = law.stored_energy_j(&step.state);
        let storage_j = after_j - before_j;
        let closure_j = storage_j + step.dissipated_j - work_j;
        let work_total_j = old.work_j + work_j;
        let heat_w = -step.dissipated_j / dt_s;
        if !(work_j.is_finite()
            && storage_j.is_finite()
            && closure_j.is_finite()
            && work_total_j.is_finite()
            && heat_w.is_finite()
            && step.dissipated_j >= 0.0
            && step.stored_change_j.is_finite()
            && step.closure_residual_j.is_finite()
            && closure_j.abs() <= tolerance_j
            && (storage_j - step.stored_change_j).abs() <= tolerance_j)
        {
            return Err(refused("sliding_interface_work_balance"));
        }
        let scale = step.dissipated_j
            + work_j.abs()
            + before_j.abs()
            + after_j.abs()
            + old.work_j.abs()
            + work_total_j.abs();
        let receipt = Receipt {
            stored_change_j: storage_j - work_j,
            dissipated_j: step.dissipated_j,
            // What the law's own closure leaves of `storage + heat - work`.
            numerical_loss_j: -closure_j,
            entropy_generated_j_per_k: step.dissipated_j / kelvin,
            energy_tolerance_j: 64.0 * f64::EPSILON * scale + f64::MIN_POSITIVE,
            residual: closure_j.abs(),
            tolerance: tolerance_j,
        };
        let state = SlidingState {
            tangential: step.state,
            reversible_j: after_j,
            work_j: work_total_j,
        };
        Ok((state, heat_w, step.traction_n, receipt))
    }
}

impl<L, F> Component for SlidingInterface<L, F>
where
    L: TangentialLaw,
    F: Fn(f64) -> Result<L, TriboError>,
{
    type State = SlidingState<L::State>;
    type Error = TriboError;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::Smooth
    }

    fn stored_energy_j(&self, state: &Self::State) -> Result<f64, TriboError> {
        finite(state.reversible_j - state.work_j, "stored_energy_j")
    }

    fn propose(
        &self,
        old: &Self::State,
        dt_s: f64,
        inputs: &[f64],
        new: &mut Self::State,
        outputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<Proposal, TriboError> {
        let [heat_w, traction_n @ ..] = outputs else {
            return Err(refused("sliding_interface_outputs"));
        };
        let (state, heat, traction, receipt) = self.step(old, dt_s, inputs)?;
        if traction_n.len() != traction.len() {
            return Err(refused("sliding_interface_outputs"));
        }
        *new = state;
        *heat_w = heat;
        traction_n.copy_from_slice(&traction);
        Ok(Proposal::Step(receipt))
    }
}

/// The load a [`FrictionInterface`] is pressed with, and the law's budgets.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct InterfaceLoad {
    /// Compressive normal load [N], non-negative.
    pub normal_n: f64,
    /// Nominal pressure [Pa], for a law that is bounded in it.
    pub nominal_pressure_pa: Option<f64>,
    /// The law's numerical and effort budgets.
    pub solve: TangentialSolveConfig,
    /// Relative momentum [N s] below which the pair is at rest: a slip
    /// whose removal over the step needs less than this is none, and a
    /// step is taken to end at rest when the kinetic traction would leave
    /// less. Positive.
    pub impulse_tolerance_n_s: f64,
    /// Largest admitted defect of the law's own work balance for a step
    /// [J], as for [`PrescribedSlip`].
    pub work_tolerance_j: f64,
}

/// State of a [`FrictionInterface`]: the law's own checkpoint, its
/// reversible energy, the slip rate the last step ended with, whether it
/// ended at rest, and, for each end whose flow was given in that step, the
/// flow it was given.
///
/// The slip rate is the flow port `a` receives less the flow port `b`
/// receives, with the signs the ports receive them in: a member of a
/// common-flow junction receives MINUS the flow of the port that sets it.
#[derive(Debug, Clone, PartialEq)]
pub struct FrictionState<S> {
    tangential: S,
    /// The law's reversible energy at `tangential` [J].
    reversible_j: f64,
    /// Flow at `a` less flow at `b` at the end of the last step [m/s]. At
    /// rest it is none, or what a stop has left of it.
    slip_m_s: f64,
    /// Whether the last step ended at rest: held, or stopped with a slip
    /// below the impulse tolerance.
    resting: bool,
    /// The mean load on the slip over the last step [N] and that step's
    /// length [s], where the step was on a line: what the bodies would
    /// have pushed the slip with had the interface returned no force.
    load: Option<(f64, f64)>,
    given_m_s: [Option<f64>; 2],
}

impl<S> FrictionState<S> {
    /// The law's checkpoint: its slip, its thermal boundary, its layers.
    #[must_use]
    pub const fn tangential(&self) -> &S {
        &self.tangential
    }

    /// The slip rate the last step ended with: flow at `a` less flow at
    /// `b` [m/s]. After a step that ended at a stop it is what the stop
    /// left, below the impulse tolerance in momentum; the next step at rest
    /// takes it out, and it is then exactly zero.
    #[must_use]
    pub const fn slip_m_s(&self) -> f64 {
        self.slip_m_s
    }

    /// Whether the last step ended at rest: held, or stopped.
    #[must_use]
    pub const fn at_rest(&self) -> bool {
        self.resting
    }

    /// The flow each given end is given in this step, and by how much it
    /// has jumped since the last one.
    ///
    /// An end whose row of the line is zero does not answer the force: its
    /// flow is given, and held over the step. Given another value than in
    /// the last step, it has jumped, and the step starts with that much
    /// more slip than the last one ended with.
    fn given_ends(&self, flows: [f64; 2], given: [bool; 2]) -> ([Option<f64>; 2], [f64; 2]) {
        let mut ends = [None; 2];
        let mut jumped = [0.0_f64; 2];
        for end in 0..2 {
            if given[end] {
                if let Some(before) = self.given_m_s[end] {
                    jumped[end] = flows[end] - before;
                }
                ends[end] = Some(flows[end]);
            }
        }
        (ends, jumped)
    }
}

/// A friction interface between two mechanical ports, with the law `law_at`
/// returns for the interface temperature of each step. See the module's
/// "Between two bodies".
pub struct FrictionInterface<L, F> {
    law_at: F,
    load: InterfaceLoad,
    ports: Vec<PortSpec>,
    /// Per end (`a`, `b`): the port of its end input, if it reads one.
    ends: [Option<usize>; 2],
    law: PhantomData<fn() -> L>,
}

impl<L, F> fmt::Debug for FrictionInterface<L, F> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("FrictionInterface")
            .field("load", &self.load)
            .finish_non_exhaustive()
    }
}

/// How a taken step's work is booked.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Booked {
    /// What the ends give up less what the law stores is heat, and rounding
    /// below zero is numerical loss.
    Heat,
    /// As `Heat`, for a step that lets go from rest on the end line: a heat
    /// below minus the work tolerance is refused.
    LetGo,
    /// Held at no end slip: the work less the storage is numerical loss,
    /// signed, with no heat.
    Held,
}

/// What a step on the line comes to.
enum OnLine<S> {
    /// The step is taken: the force returned at `a` [N], the mean slip
    /// rate [m/s], the slip rate the step ends with [m/s], and the law's
    /// own step.
    Taken {
        force_n: f64,
        mean_m_s: f64,
        ends_m_s: f64,
        /// Whether the step ends at rest.
        resting: bool,
        law: TangentialStep<S>,
        booked: Booked,
    },
    /// An event lies inside the step, this long after its start [s]: the
    /// slip reaches zero, or the load on a pair at rest reaches what the
    /// law holds.
    EventAfter(f64),
}

impl<L, F> FrictionInterface<L, F>
where
    L: TangentialLaw,
    F: Fn(f64) -> Result<L, TriboError>,
{
    /// An interface pressed with `load`, with the law `law_at` gives at an
    /// interface temperature [K].
    ///
    /// # Errors
    /// A normal load or nominal pressure that is negative or not finite, a
    /// work tolerance that is negative or not finite, or solve budgets the
    /// law's owner refuses.
    pub fn new(law_at: F, load: InterfaceLoad) -> Result<Self, TriboError> {
        load.solve.validate()?;
        crate::nonnegative_finite(load.normal_n, "normal_force")?;
        if let Some(pressure_pa) = load.nominal_pressure_pa {
            crate::nonnegative_finite(pressure_pa, "contact_pressure_pa")?;
        }
        crate::nonnegative_finite(load.work_tolerance_j, "work_tolerance_j")?;
        crate::positive_finite(load.impulse_tolerance_n_s, "impulse_tolerance_n_s")?;
        let name = |id: &str| StableId::new(id).map_err(|_| refused("port_identifier"));
        let mechanical = |id: &str| {
            PortSpec::power(
                name(id)?,
                PortKind::MechanicalForceVelocity,
                PortRole::FlowIn,
                NonZeroUsize::MIN,
            )
            .map_err(|_| refused("mechanical_port"))
        };
        Ok(Self {
            law_at,
            load,
            ports: vec![
                mechanical("a")?,
                mechanical("b")?,
                PortSpec::heat(name("heat")?, PortRole::EffortIn)
                    .map_err(|_| refused("heat_port"))?,
            ],
            ends: [None; 2],
            law: PhantomData,
        })
    }

    /// The interface reading, at each end marked in `ends` (`a`, `b`), the
    /// flow its neighbours tell for the END of the step: a signal input
    /// `a.end` or `b.end` of width one after `heat`, which
    /// [`LoadLine::line_ends`] names. See the module's "Beside a neighbour
    /// that tells the end of its flow". An end whose flow a kernel
    /// prescribes takes none.
    ///
    /// # Errors
    /// No end marked, or end inputs already read.
    pub fn reading_ends(mut self, ends: [bool; 2]) -> Result<Self, TriboError> {
        if ends == [false; 2] || self.ends != [None; 2] {
            return Err(refused("friction_interface_end_inputs"));
        }
        for (end, id) in [(0, "a.end"), (1, "b.end")] {
            if ends[end] {
                let id = StableId::new(id).map_err(|_| refused("port_identifier"))?;
                self.ends[end] = Some(self.ports.len());
                self.ports.push(PortSpec::signal_in(id, NonZeroUsize::MIN));
            }
        }
        Ok(self)
    }

    /// How many ends read an end input.
    fn told(&self) -> usize {
        self.ends.iter().flatten().count()
    }

    /// Where an end reads an end input: the released slip at the step's end
    /// and what a unit of the force takes from it, from the end inputs'
    /// flows and their rows of the line (of the lengths the ends need). An
    /// end that reads none must be given, and ends the step at its flow.
    fn end_line(
        &self,
        flows: [f64; 2],
        given: [bool; 2],
        end_flows: &[f64],
        end_rows: &[f64],
    ) -> Result<Option<(f64, f64)>, TriboError> {
        if self.told() == 0 {
            return Ok(None);
        }
        let mut told = end_flows.iter().zip(end_rows.as_chunks::<2>().0);
        let mut line = [(0.0_f64, 0.0_f64); 2];
        for (((end, reads), given), flow) in line.iter_mut().zip(self.ends).zip(given).zip(flows) {
            *end = if reads.is_some() {
                let (&released, row) = told
                    .next()
                    .ok_or_else(|| refused("friction_interface_inputs"))?;
                // Its row of the line taken along the efforts (F, -F).
                (released, row[0] - row[1])
            } else if given {
                (flow, 0.0)
            } else {
                return Err(refused("friction_interface_untold_end"));
            };
        }
        let [(released_a, taken_a), (released_b, taken_b)] = line;
        Ok(Some((released_a - released_b, taken_a - taken_b)))
    }

    /// The load this interface is pressed with.
    #[must_use]
    pub const fn load(&self) -> &InterfaceLoad {
        &self.load
    }

    /// The unloaded state of the law at `interface_temperature_k`, sliding
    /// at `slip_m_s` (the flow port `a` receives less the flow port `b`
    /// receives; zero is rest): the state to hand to a kernel with this
    /// component.
    ///
    /// # Errors
    /// A temperature that is not positive and finite, a slip rate that is
    /// not finite, the source's refusal, or a reversible energy that is not
    /// finite.
    pub fn start(
        &self,
        interface_temperature_k: f64,
        slip_m_s: f64,
    ) -> Result<FrictionState<L::State>, TriboError> {
        let law = self.law_at_k(interface_temperature_k)?;
        let tangential = law.initial_state();
        let reversible_j = finite(law.stored_energy_j(&tangential), "stored_energy_j")?;
        Ok(FrictionState {
            tangential,
            reversible_j,
            slip_m_s: finite(slip_m_s, "slip_m_s")?,
            resting: slip_m_s == 0.0,
            load: None,
            given_m_s: [None; 2],
        })
    }

    fn law_at_k(&self, interface_temperature_k: f64) -> Result<L, TriboError> {
        if !(interface_temperature_k.is_finite() && interface_temperature_k > 0.0) {
            return Err(refused("interface_temperature_k"));
        }
        (self.law_at)(interface_temperature_k)
    }

    /// The law's step from `old` over `dt_s` at a mean slip rate, as the
    /// law resolves it against a compliance [m/N] (zero: the slip is
    /// prescribed).
    fn advance(
        &self,
        law: &L,
        old: &L::State,
        dt_s: f64,
        kelvin: f64,
        released_m_s: f64,
        compliance_m_per_n: f64,
    ) -> Result<TangentialStep<L::State>, TriboError> {
        law.advance(
            old,
            &TangentialStepInput {
                normal_n: self.load.normal_n,
                slip_increment_m: [dt_s * released_m_s, 0.0],
                dt_s,
                contact_temperature_k: Some(kelvin),
                nominal_pressure_pa: self.load.nominal_pressure_pa,
                compliance_m_per_n,
                solve: self.load.solve,
            },
        )
    }

    /// The law's kinetic step at a prescribed mean slip rate that is not
    /// zero.
    fn sliding(
        &self,
        law: &L,
        old: &L::State,
        dt_s: f64,
        kelvin: f64,
        mean_m_s: f64,
    ) -> Result<TangentialStep<L::State>, TriboError> {
        let step = self.advance(law, old, dt_s, kelvin, mean_m_s, 0.0)?;
        if step.regime != FrictionRegime::Sliding {
            return Err(refused("friction_interface_law_holds_a_prescribed_slip"));
        }
        Ok(step)
    }

    /// Where, in a step at rest, the load on the pair reaches what the law
    /// holds: the fraction of the step, if the load is held at the step's
    /// start, is not at its end, and the instant between lies inside the
    /// step by more than the impulse tolerance at either end.
    ///
    /// `load_n` is the step's mean load and `last` the last step's with
    /// its length. The load at the step's two ends is read off the
    /// straight line through the two means, which is exact for a load
    /// that is linear in time whatever the two lengths. (Not
    /// `end = 2 mean - start` carried from step to step: where the load
    /// curves, the error of that record changes sign every step and never
    /// decays.) What the law holds is asked of the law, by its own static
    /// map at loads between the two ends; no capacity is restated here.
    #[allow(clippy::too_many_arguments)]
    fn lets_go_inside(
        &self,
        law: &L,
        old: &L::State,
        dt_s: f64,
        kelvin: f64,
        c: f64,
        load_n: f64,
        last: (f64, f64),
    ) -> Result<Option<f64>, TriboError> {
        let (last_n, last_s) = last;
        let rate_n_s = (load_n - last_n) / (0.5 * (dt_s + last_s));
        let (from_n, to_n) = (load_n - 0.5 * dt_s * rate_n_s, load_n + 0.5 * dt_s * rate_n_s);
        // The law's verdict at a load, and what the law says it cost: its
        // own count of iterations, none where the verdict is a closed form.
        let holds = |force_n: f64| -> Result<(bool, usize), TriboError> {
            let step = self.advance(law, old, dt_s, kelvin, c * force_n, dt_s * c)?;
            Ok((step.regime == FrictionRegime::Sticking, step.iterations))
        };
        let (held_at_start, mut cost) = holds(from_n)?;
        if !held_at_start {
            return Ok(None);
        }
        let (held_at_end, at_end) = holds(to_n)?;
        if held_at_end {
            return Ok(None);
        }
        cost = cost.max(at_end);
        // The load passes the limit at a constant rate, so letting go a
        // time `t` early or late is a momentum of half the rate times
        // `t` squared.
        let swing_n_s = dt_s * (to_n - from_n).abs();
        let tolerance_n_s = self.load.impulse_tolerance_n_s;
        // The largest fraction of the step at which the law still holds,
        // narrowed by halving. How far: the instant matters to the impulse
        // tolerance and no further, and a law that searches for its
        // verdict is asked no nearer its limit than that. Just beyond what
        // such a law holds its sliding root is next to zero, and its own
        // search for it takes more with every halving: narrowed to
        // rounding, the stable-branch law on a bowed string's line was
        // asked 51 times in one step, answered the last of them in 123 and
        // 125 of the 128 iterations it may spend, and refused the next. A
        // law whose verdicts cost it nothing is narrowed to rounding,
        // which then costs nothing either.
        let (mut held, mut gone) = (0.0_f64, 1.0_f64);
        let mut asked = 2_usize;
        while gone - held > 4.0 * f64::EPSILON && asked < self.load.solve.max_iterations {
            let width = gone - held;
            if cost > 0 && 0.5 * width * width * swing_n_s <= tolerance_n_s {
                break;
            }
            asked += 1;
            let middle = 0.5 * (held + gone);
            let (held_there, there) = holds(from_n + middle * (to_n - from_n))?;
            cost = cost.max(there);
            if held_there {
                held = middle;
            } else {
                gone = middle;
            }
        }
        // An instant within the tolerance of a step's end is that end.
        let (early, late) = (0.5 * held * held * swing_n_s, 0.5 * (1.0 - held) * (1.0 - held) * swing_n_s);
        Ok((early > tolerance_n_s && late > tolerance_n_s).then_some(held))
    }

    /// One step on the load line `mean = released - c force`, from the slip
    /// rate `start`; `resting`: the last step ended at rest and no given
    /// end has jumped since; `last`: the last step's mean load and length.
    /// See the module's table.
    #[allow(clippy::too_many_arguments)]
    fn on_line(
        &self,
        law: &L,
        old: &L::State,
        dt_s: f64,
        kelvin: f64,
        start: f64,
        resting: bool,
        last: Option<(f64, f64)>,
        released: f64,
        c: f64,
    ) -> Result<OnLine<L::State>, TriboError> {
        if !(c > 0.0) {
            // No compliance: both flows are given, or what the interface is
            // joined to is rigid to within rounding. The slip is what it
            // is. With none at all nothing says what the pair holds, and
            // the law returns no traction.
            let step = if released == 0.0 {
                self.advance(law, old, dt_s, kelvin, 0.0, 0.0)?
            } else {
                self.sliding(law, old, dt_s, kelvin, released)?
            };
            let force_n = step.traction_n[0];
            let mean_m_s = released - c * force_n;
            return Ok(OnLine::Taken {
                force_n,
                mean_m_s,
                ends_m_s: 2.0 * mean_m_s - start,
                resting: mean_m_s == 0.0,
                law: step,
                booked: Booked::Heat,
            });
        }
        // `dt / (2 c)` is the inertia the slip moves against, so a force
        // times `dt` is the relative momentum it gives or takes. A slip
        // below the tolerance in that measure is rest.
        let impulse_n_s = self.load.impulse_tolerance_n_s;
        if resting || dt_s * (0.5 * start.abs() / c) <= impulse_n_s {
            // At rest: the law's own static return map on the line, asked
            // for the slip that ends the step at none. With no slip to
            // start with that is the released flow itself; after a stop it
            // takes out what the stop left. The law holds with no slip
            // beyond that, or slides at the traction it solves.
            //
            // The map decides for the whole step from the step's mean.
            // Where the load reaches what the law holds inside the step,
            // the step is cut there first: the cut step is held, and the
            // next one lets go from its start.
            if let Some(last) = last {
                let load_n = (released - start) / c;
                if let Some(fraction) =
                    self.lets_go_inside(law, old, dt_s, kelvin, c, load_n, last)?
                {
                    let shorter_s = fraction * dt_s;
                    if shorter_s > 0.0 && shorter_s < dt_s {
                        return Ok(OnLine::EventAfter(shorter_s));
                    }
                }
            }
            let aimed = released - 0.5 * start;
            let step = self.advance(law, old, dt_s, kelvin, aimed, dt_s * c)?;
            let force_n = step.traction_n[0];
            return match step.regime {
                FrictionRegime::Sticking => {
                    // Exact on purpose: held is no slip at all.
                    #[allow(clippy::float_cmp)]
                    if step.slip_increment_m != [0.0; 2] {
                        return Err(refused("friction_interface_static_map_kept_slip"));
                    }
                    Ok(OnLine::Taken {
                        force_n,
                        mean_m_s: 0.5 * start,
                        ends_m_s: 0.0,
                        resting: true,
                        law: step,
                        booked: Booked::Heat,
                    })
                }
                FrictionRegime::Sliding => {
                    let mean_m_s = released - c * force_n;
                    Ok(OnLine::Taken {
                        force_n,
                        mean_m_s,
                        ends_m_s: 2.0 * mean_m_s - start,
                        resting: false,
                        law: step,
                        booked: Booked::Heat,
                    })
                }
            };
        }
        // Sliding. A step that ends at zero slip has the mean `start / 2`
        // whatever its length, so the kinetic traction of the stop is the
        // law's at that mean, and `along` is the force that would end this
        // step there.
        let way = start.signum();
        let stop = self.sliding(law, old, dt_s, kelvin, 0.5 * start)?;
        let kinetic_n = way * stop.traction_n[0];
        if !(kinetic_n >= 0.0) {
            return Err(refused("friction_interface_traction_drives_the_slip"));
        }
        let along_n = way * (released - 0.5 * start) / c;
        let band_n = self.load.solve.tolerance(along_n, kinetic_n)? + impulse_n_s / dt_s;
        if along_n > kinetic_n + band_n {
            let lower = 0.5 * start.abs();
            return self.slides_on(
                law, old, dt_s, kelvin, start, released, c, kinetic_n, lower, None,
            );
        }
        if along_n < kinetic_n - band_n {
            // The kinetic traction stops the slip inside the step. Where,
            // by the line of this step: for bodies under held forces what
            // the released flow adds to the starting slip, and the
            // compliance, both grow in proportion to the step.
            let shorter_s = dt_s * (0.5 * start.abs()) / (c * kinetic_n - way * (released - start));
            if shorter_s > 0.0 && shorter_s < dt_s {
                return Ok(OnLine::EventAfter(shorter_s));
            }
        }
        // The step ends at the stop, to the band, under the law's kinetic
        // traction. What that leaves of the slip is below the impulse
        // tolerance; it is kept, and the first step at rest takes it out.
        //
        // The force that would end this step at exactly no slip,
        // `along`, is not returned here. It is stiff in the released flow
        // (its slope is `1 / c`) and it would be returned only inside the
        // band, which is far narrower than the step by which a kernel
        // differences this component: the kernel would be handed the
        // sliding rows' slope of zero for a row of slope `1 / c`, and a
        // step that needs one more update at such an iterate does not
        // converge. (As first landed, the belt of the graph tests at
        // 0.137 m/s and 0.25 ms was refused at its third catch.)
        let force_n = stop.traction_n[0];
        let mean_m_s = released - c * force_n;
        Ok(OnLine::Taken {
            force_n,
            mean_m_s,
            ends_m_s: 2.0 * mean_m_s - start,
            resting: true,
            law: stop,
            booked: Booked::Heat,
        })
    }

    /// One step on the line beside a neighbour that tells the end of its
    /// flow: the mean line `mean = released - c force` and the end line
    /// `end = released_end - c_end force`. See the module's "Beside a
    /// neighbour that tells the end of its flow".
    #[allow(clippy::too_many_arguments)]
    fn on_end_line(
        &self,
        law: &L,
        old: &L::State,
        dt_s: f64,
        kelvin: f64,
        start: f64,
        resting: bool,
        last: Option<(f64, f64)>,
        released: f64,
        c: f64,
        (released_end, c_end): (f64, f64),
    ) -> Result<OnLine<L::State>, TriboError> {
        // `dt / c_end` is the inertia the slip moves against.
        let impulse_n_s = self.load.impulse_tolerance_n_s;
        if resting || dt_s * (start.abs() / c_end) <= impulse_n_s {
            // At rest: the break-away located from the loads of this step
            // and the last, then the law's own return map on the END line,
            // for the slip that ends the step at none. What a stop left is
            // in the released end already: it is the neighbour's own.
            if let Some(last) = last {
                let load_n = (released_end - start) / c_end;
                if let Some(fraction) =
                    self.lets_go_inside(law, old, dt_s, kelvin, c_end, load_n, last)?
                {
                    let shorter_s = fraction * dt_s;
                    if shorter_s > 0.0 && shorter_s < dt_s {
                        return Ok(OnLine::EventAfter(shorter_s));
                    }
                }
            }
            let step = self.advance(law, old, dt_s, kelvin, released_end, dt_s * c_end)?;
            let force_n = step.traction_n[0];
            let mean_m_s = released - c * force_n;
            return match step.regime {
                FrictionRegime::Sticking => {
                    // Exact on purpose: held is no slip at all.
                    #[allow(clippy::float_cmp)]
                    if step.slip_increment_m != [0.0; 2] {
                        return Err(refused("friction_interface_static_map_kept_slip"));
                    }
                    Ok(OnLine::Taken {
                        force_n,
                        mean_m_s,
                        ends_m_s: 0.0,
                        resting: true,
                        law: step,
                        booked: Booked::Held,
                    })
                }
                FrictionRegime::Sliding => Ok(OnLine::Taken {
                    force_n,
                    mean_m_s,
                    ends_m_s: released_end - c_end * force_n,
                    resting: false,
                    law: step,
                    booked: Booked::LetGo,
                }),
            };
        }
        // Sliding. `along` ends the step at no slip; the step that does has
        // the mean speed `stop`, and the law's kinetic traction there is
        // the stop's, asked at a prescribed mean as above.
        let way = start.signum();
        let along_n = way * released_end / c_end;
        let stop_m_s = way * released - c * along_n;
        if !(stop_m_s > 0.0) {
            return self.against_the_slip(
                law,
                old,
                dt_s,
                kelvin,
                start,
                released,
                c,
                (released_end, c_end),
            );
        }
        let stop = self.sliding(law, old, dt_s, kelvin, way * stop_m_s)?;
        let kinetic_n = way * stop.traction_n[0];
        if !(kinetic_n >= 0.0) {
            return Err(refused("friction_interface_traction_drives_the_slip"));
        }
        let band_n = self.load.solve.tolerance(along_n, kinetic_n)? + impulse_n_s / dt_s;
        if along_n > kinetic_n + band_n {
            return self.slides_on(
                law,
                old,
                dt_s,
                kelvin,
                start,
                released,
                c,
                kinetic_n,
                stop_m_s,
                Some((released_end, c_end)),
            );
        }
        if along_n < kinetic_n - band_n {
            // The end slip passes zero inside the step. Where, by the end
            // line of this step: exact for bodies under held forces.
            let shorter_s = dt_s * start.abs() / (c_end * kinetic_n - way * (released_end - start));
            if shorter_s > 0.0 && shorter_s < dt_s {
                return Ok(OnLine::EventAfter(shorter_s));
            }
        }
        // The step ends at the stop, to the band, under the law's kinetic
        // traction (the wide row, as above); the leftover is kept.
        let force_n = stop.traction_n[0];
        Ok(OnLine::Taken {
            force_n,
            mean_m_s: released - c * force_n,
            ends_m_s: released_end - c_end * force_n,
            resting: true,
            law: stop,
            booked: Booked::Heat,
        })
    }

    /// A sliding step beside a told end whose stop step's mean would run
    /// against the slip: `way r - c along <= 0`.
    ///
    /// The law's own return map on the end line says whether it holds `along`
    /// at all. Where it does not, the end cannot stop in this step, and the
    /// step slides on, its mean bracketed from rest. Seen on a bowed string,
    /// `along` 0.726 N against the 0.5465 N its law holds, and on a mode on a
    /// belt, 5.0094 N against 5 N, with stop speeds of -0.0208 and -3.1e-6 m/s.
    ///
    /// Where the law holds `along`, its kinetic traction at the STARTING slip,
    /// which is along the slip, says where the end line stops. Below that
    /// traction, `along` cannot hold the step whole, and the step is cut there.
    /// The starting traction only places the cut: the shorter step the
    /// kernel proposes next is a step of its own, taken by the rows of
    /// [`Self::on_end_line`], which ask the law at that step's mean. Seen on
    /// the bow at its 7th capture: `along` 0.1411 N, which the law holds; a
    /// stop speed of -2.29e-3 m/s; a traction of 0.4985 N at the starting
    /// 0.0182 m/s; and an end line that stops 0.0719 of the way through the
    /// step. Anything else is refused.
    #[allow(clippy::too_many_arguments)]
    fn against_the_slip(
        &self,
        law: &L,
        old: &L::State,
        dt_s: f64,
        kelvin: f64,
        start: f64,
        released: f64,
        c: f64,
        (released_end, c_end): (f64, f64),
    ) -> Result<OnLine<L::State>, TriboError> {
        let way = start.signum();
        let along_n = way * released_end / c_end;
        let reach = self.advance(law, old, dt_s, kelvin, released_end, dt_s * c_end)?;
        if reach.regime == FrictionRegime::Sliding {
            let kinetic_n = way * reach.traction_n[0];
            if !(kinetic_n > 0.0 && way * released - c * kinetic_n > 0.0) {
                return Err(refused("friction_interface_stop_against_the_slip"));
            }
            return self.slides_on(
                law,
                old,
                dt_s,
                kelvin,
                start,
                released,
                c,
                kinetic_n,
                0.0,
                Some((released_end, c_end)),
            );
        }
        let starting_n = way * self.sliding(law, old, dt_s, kelvin, start)?.traction_n[0];
        let band_n = self.load.solve.tolerance(along_n, starting_n)?
            + self.load.impulse_tolerance_n_s / dt_s;
        if along_n < starting_n - band_n {
            let shorter_s =
                dt_s * start.abs() / (c_end * starting_n - way * (released_end - start));
            if shorter_s > 0.0 && shorter_s < dt_s {
                return Ok(OnLine::EventAfter(shorter_s));
            }
        }
        Err(refused("friction_interface_stop_against_the_slip"))
    }

    /// A step that slides on the way it started: the mean speed along the
    /// slip solves `speed = way released - c K(speed)`, with `K` the law's
    /// kinetic traction. It is bracketed between `lower`, the stop step's
    /// mean speed, where the traction falls short of what the line needs,
    /// and the released speed, where the line needs none. The step ends by
    /// the midpoint relation, or on the end line `end_line` where an end is
    /// told.
    #[allow(clippy::too_many_arguments)]
    fn slides_on(
        &self,
        law: &L,
        old: &L::State,
        dt_s: f64,
        kelvin: f64,
        start: f64,
        released: f64,
        c: f64,
        kinetic_n: f64,
        lower: f64,
        end_line: Option<(f64, f64)>,
    ) -> Result<OnLine<L::State>, TriboError> {
        let way = start.signum();
        let taken = |step: TangentialStep<L::State>| {
            let force_n = step.traction_n[0];
            let mean_m_s = released - c * force_n;
            OnLine::Taken {
                force_n,
                mean_m_s,
                ends_m_s: match end_line {
                    None => 2.0 * mean_m_s - start,
                    Some((released_end, c_end)) => released_end - c_end * force_n,
                },
                resting: false,
                law: step,
                booked: Booked::Heat,
            }
        };
        let solve = self.load.solve;
        // The law's traction at a speed less what the line needs there.
        let at = |speed: f64, needed_n: f64| -> Result<_, TriboError> {
            let step = self.sliding(law, old, dt_s, kelvin, way * speed)?;
            let traction_n = way * step.traction_n[0];
            let residual_n = finite(traction_n - needed_n, "constitutive_residual_n")?;
            Ok((step, residual_n, solve.tolerance(traction_n, needed_n)?))
        };
        // First the speed the stop's traction leaves: the answer wherever
        // the traction does not depend on the speed.
        let (mut lower, mut upper) = (lower, way * released);
        let first = way * released - c * kinetic_n;
        let (step, mut residual_n, mut tolerance_n) = at(first, kinetic_n)?;
        if residual_n.abs() <= tolerance_n {
            return Ok(taken(step));
        }
        let along_n = (way * released - lower) / c;
        let (mut lower_n, mut upper_n) = (kinetic_n - along_n, f64::NAN);
        if residual_n < 0.0 {
            (lower, lower_n) = (first, residual_n);
        } else {
            (upper, upper_n) = (first, residual_n);
        }
        let mut used = 1_usize;
        if upper_n.is_nan() {
            // At the released speed the line needs no force.
            let (step, residual, tolerance) = at(upper, 0.0)?;
            used += 1;
            if residual.abs() <= tolerance {
                return Ok(taken(step));
            }
            (upper_n, residual_n, tolerance_n) = (residual, residual, tolerance);
        }
        while used < solve.max_iterations {
            used += 1;
            let width = upper - lower;
            // A secant step where it keeps a sixteenth of the bracket on
            // either side; else the middle.
            let secant = lower - lower_n * width / (upper_n - lower_n);
            let speed = if secant.is_finite()
                && secant > lower + width / 16.0
                && secant < upper - width / 16.0
            {
                secant
            } else {
                lower + 0.5 * width
            };
            if !(speed > lower && speed < upper) {
                break;
            }
            let (step, residual, tolerance) = at(speed, (way * released - speed) / c)?;
            if residual.abs() <= tolerance {
                return Ok(taken(step));
            }
            (residual_n, tolerance_n) = (residual, tolerance);
            if residual < 0.0 {
                (lower, lower_n) = (speed, residual);
            } else {
                (upper, upper_n) = (speed, residual);
            }
        }
        Err(TriboError::TangentialSolve {
            residual_n,
            tolerance_n,
            iterations: used,
        })
    }
}

impl<L, F> Component for FrictionInterface<L, F>
where
    L: TangentialLaw,
    F: Fn(f64) -> Result<L, TriboError>,
{
    type State = FrictionState<L::State>;
    type Error = TriboError;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::NonSmooth
    }

    fn stored_energy_j(&self, state: &Self::State) -> Result<f64, TriboError> {
        finite(state.reversible_j, "stored_energy_j")
    }

    fn propose(
        &self,
        old: &Self::State,
        dt_s: f64,
        inputs: &[f64],
        new: &mut Self::State,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, TriboError> {
        // The line's square matrix and one row per end input, all zero.
        let zero = [0.0; 8];
        let compliance = &zero[..2 * (2 + self.told())];
        self.propose_on_line(old, dt_s, inputs, compliance, new, outputs, cx)
    }
}

impl<L, F> LoadLine for FrictionInterface<L, F>
where
    L: TangentialLaw,
    F: Fn(f64) -> Result<L, TriboError>,
{
    fn line_ports(&self) -> &[usize] {
        &[0, 1]
    }

    fn line_ends(&self) -> &[Option<usize>] {
        &self.ends
    }

    // One pass from the inputs to the receipt, as before the end line joined
    // it: its reading and its dispatch took the pass past 100 lines.
    #[allow(clippy::too_many_lines)]
    fn propose_on_line(
        &self,
        old: &Self::State,
        dt_s: f64,
        inputs: &[f64],
        compliance: &[f64],
        new: &mut Self::State,
        outputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<Proposal, TriboError> {
        // After the temperature, the end inputs; after the line's square
        // matrix, their rows.
        let (
            &[flow_a, flow_b, kelvin, ref end_flows @ ..],
            &[c_aa, c_ab, c_ba, c_bb, ref end_rows @ ..],
        ) = (inputs, compliance)
        else {
            return Err(refused("friction_interface_inputs"));
        };
        if end_flows.len() != self.told() || end_rows.len() != 2 * self.told() {
            return Err(refused("friction_interface_inputs"));
        }
        let [force_a, force_b, heat_w] = outputs else {
            return Err(refused("friction_interface_outputs"));
        };
        if !(dt_s.is_finite() && dt_s > 0.0) {
            return Err(refused("dt_s"));
        }
        for value in [flow_a, flow_b, c_aa, c_ab, c_ba, c_bb] {
            finite(value, "friction_interface_inputs")?;
        }
        for &value in end_flows.iter().chain(end_rows) {
            finite(value, "friction_interface_inputs")?;
        }
        let law = self.law_at_k(kelvin)?;
        // As for the prescribed interface: the held state must be worth to
        // this law exactly what it was left worth.
        let before_j = law.stored_energy_j(&old.tangential);
        #[allow(clippy::float_cmp)]
        if before_j != old.reversible_j {
            return Err(refused("friction_interface_law_changed_stored_energy"));
        }
        let given = [c_aa == 0.0 && c_ab == 0.0, c_ba == 0.0 && c_bb == 0.0];
        let (ends, jumped) = old.given_ends([flow_a, flow_b], given);
        // The force acts as (e, -e) on the two ends, so a unit of it takes
        // (1, -1) C (1, -1)^T away from the slip.
        let c = c_aa - c_ab - c_ba + c_bb;
        let released = flow_a - flow_b;
        let start = if given[0] && given[1] {
            released
        } else {
            old.slip_m_s + jumped[0] - jumped[1]
        };
        // At rest for this step too, unless a given end has jumped.
        #[allow(clippy::float_cmp)]
        let resting = old.resting && jumped == [0.0; 2];
        let end_line = self.end_line([flow_a, flow_b], given, end_flows, end_rows)?;
        let on_line = match end_line {
            Some(end_line) if c > 0.0 => {
                if !(end_line.1 > 0.0) {
                    return Err(refused("friction_interface_end_compliance"));
                }
                self.on_end_line(
                    &law,
                    &old.tangential,
                    dt_s,
                    kelvin,
                    start,
                    resting,
                    old.load,
                    released,
                    c,
                    end_line,
                )?
            }
            _ => self.on_line(
                &law,
                &old.tangential,
                dt_s,
                kelvin,
                start,
                resting,
                old.load,
                released,
                c,
            )?,
        };
        let (force_n, mean_m_s, ends_m_s, resting, step, booked) =
            match on_line {
                OnLine::EventAfter(shorter_s) => return Ok(Proposal::Shorten { dt_s: shorter_s }),
                OnLine::Taken {
                    force_n,
                    mean_m_s,
                    ends_m_s,
                    resting,
                    law,
                    booked,
                } => (force_n, mean_m_s, ends_m_s, resting, law, booked),
            };
        // The law's own accounts for the step it was asked for, admitted
        // as the prescribed interface admits them: its storage, its
        // dissipation and the work of its traction over its slip close, and
        // its reported stored change is that of its stored energies.
        let after_j = law.stored_energy_j(&step.state);
        let storage_j = after_j - before_j;
        let law_work_j = step.traction_n[0] * step.slip_increment_m[0]
            + step.traction_n[1] * step.slip_increment_m[1];
        let closure_j = storage_j + step.dissipated_j - law_work_j;
        let tolerance_j = self.load.work_tolerance_j;
        // What the two ends give up is the force times the mean slip that
        // crossed; what the law does not store of it leaves through `heat`.
        // (It is what the ports carried, which is the law's own dissipation
        // only to the tolerance of the row: the law is asked at the slip a
        // row aims at, and the line settles a tolerance away.)
        let work_j = dt_s * force_n * mean_m_s;
        let heat_j = work_j - storage_j;
        // A step at rest that takes out what a stop left moves the pair by
        // that much against the force that holds it, either way: work of
        // at most the force times `c` times the impulse tolerance, which
        // is no heat where it is negative.
        let taken_out_j = force_n.abs() * c.max(0.0) * self.load.impulse_tolerance_n_s;
        if booked == Booked::LetGo && !(heat_j >= -tolerance_j) {
            return Err(refused("friction_interface_let_go_heat"));
        }
        // Held at no end slip, the work is a defect of the step, either way.
        let held = booked == Booked::Held;
        if !(after_j.is_finite()
            && work_j.is_finite()
            && ends_m_s.is_finite()
            && step.dissipated_j >= 0.0
            && closure_j.abs() <= tolerance_j
            && (storage_j - step.stored_change_j).abs() <= tolerance_j
            && (held || heat_j >= -(tolerance_j + taken_out_j)))
        {
            return Err(refused("friction_interface_work_balance"));
        }
        let dissipated_j = if held { 0.0 } else { heat_j.max(0.0) };
        *new = FrictionState {
            tangential: step.state,
            reversible_j: after_j,
            slip_m_s: ends_m_s,
            resting,
            load: (c > 0.0).then(|| match end_line {
                Some((released_end, c_end)) => ((released_end - start) / c_end, dt_s),
                None => ((released - start) / c, dt_s),
            }),
            given_m_s: ends,
        };
        *force_a = force_n;
        *force_b = -force_n;
        *heat_w = -dissipated_j / dt_s;
        // The two ends' powers are each as large as their own velocities
        // make them; their sum is rounded at that size, not at the loss's.
        let at_ends = flow_a.abs()
            + flow_b.abs()
            + force_n.abs() * (c_aa.abs() + c_ab.abs() + c_ba.abs() + c_bb.abs());
        let scale = dt_s * force_n.abs() * at_ends + dissipated_j + before_j.abs() + after_j.abs();
        Ok(Proposal::Step(Receipt {
            stored_change_j: storage_j,
            dissipated_j,
            // A loss that rounding left below zero is no heat.
            numerical_loss_j: heat_j - dissipated_j,
            entropy_generated_j_per_k: dissipated_j / kelvin,
            energy_tolerance_j: 64.0 * f64::EPSILON * scale + f64::MIN_POSITIVE,
            residual: closure_j.abs(),
            tolerance: tolerance_j,
        }))
    }
}
