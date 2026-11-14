//! The Hertz patch on the composition contract: one circular patch between
//! two bodies, its joint slip and spin friction from the owner's law
//! ([`HertzWrenchContact`]), and its heat.
//!
//! Port `wrench` (power, width 3, receives flow) carries the patch's relative
//! motion in the owner's radius-scaled coordinates: the slip velocity of side
//! a against side b along the two tangent axes, and the patch radius `a`
//! times their relative spin about the normal, `[v_x, v_y, a w]` [m/s]. It
//! returns the effort that takes power out of that motion, minus the wrench
//! acting on side a, `-[F_x, F_y, M / a]` [N], so the power it receives is
//! the work the patch absorbs. All three coefficients are a force and a
//! velocity, so the port is one mechanical port. Which body ports they come
//! from, through which arms and axes, is the system's map, not this
//! component's.
//!
//! Port `heat` receives the interface temperature [K] and returns the heat
//! rate into the interface [W], minus the owner's dissipation, as
//! [`FrictionInterface`](crate::tangential::component::FrictionInterface)'s
//! does. Where the heat goes is decided by what the port is connected to.
//!
//! The law at a temperature is the caller's `law_at(T)`, a pure function,
//! asked at every evaluation of a step. It may change the static and kinetic
//! coefficients with temperature. The normal load and the radius are the
//! component's own ([`HertzPatchLoad`]), and a law with another load or
//! radius is refused: the radius scales the port's third coefficient.
//!
//! # A step
//!
//! Dry friction holds any wrench inside its static set at no motion, so the
//! effort is not a function of the flow. The component is added on the
//! kernel's non-smooth tier ([`LoadLine`], `add_on_line` in `fs-system`),
//! which hands it the line its neighbours present, `mean = released - C
//! effort`. The owner is asked once, as the module's mechanical owner asks it
//! (`fs_mbd::tangential_friction::hertz::PreparedRigidHertzPair`):
//!
//! - its free displacement is `dt released` and its compliance `dt C_s`, with
//!   `C_s` the symmetric part of the line (below);
//! - its phase is `Kinetic` when any coefficient of the motion the step
//!   starts with is not zero, and `StaticAdmissible` otherwise (the module's
//!   `accepted_phase`): a moving pair is never captured by the larger static
//!   set;
//! - the effort returned is the owner's opposing wrench, and the step's mean
//!   motion, which the port carries, is the line's at that effort;
//! - the motion the step ends with is twice the owner's mean motion less the
//!   start: the midpoint rule of the kernel's bodies, and the module's own
//!   endpoint relation. A pair the owner holds at rest therefore ends the
//!   step exactly at rest. A coefficient whose row of the line is zero is
//!   given, and is held over the step.
//!
//! The heat is the owner's own dissipation, `dissipated_energy_j`: the native
//! Hertz support at the owner's actual motion. What the port carried, `dt
//! effort . mean`, less that heat is booked as numerical loss, and a step
//! whose loss exceeds the owner's energy tolerance (with the rounding of the
//! port's work) is refused (`hertz_patch_work_balance`). Entropy generated:
//! the heat over the interface temperature. Stored energy: none, since the
//! owner's rung keeps no elastic memory.
//!
//! The symmetric part. The owner's step is a convex projection and takes a
//! symmetric response; it refuses an asymmetry above its
//! `relative_matrix_tolerance`, at most `1e-8` of the largest entry. A kernel
//! line need not be symmetric. One differenced from a neighbour without a
//! declared tangent carries rounding of about the square root of the rounding
//! unit (`fs-system`, `system/line.rs`). A spinning body's line also carries
//! its gyroscopic coupling: a torque across the spin axis precesses the body
//! within the step, which is physics, not rounding. The antisymmetric part
//! `K` does no work (`e . K e = 0`), but it moves: the line's motion over the
//! step is `dt (released - C_s e - K e)`. The owner is handed `C_s` and never
//! sees a non-symmetric response. The step is taken on the full line: the
//! owner steps with the free motion `dt released`; where `dt K e` at its
//! effort is not exactly zero, it steps again with `dt (released - K e)` at
//! the effort `e` of its previous step, until its force and moment change by
//! no more than its own `force_tolerance_n` and `moment_tolerance_nm`. The
//! step then ends on the full line to the owner's tolerances. Not settled
//! after the owner's `maximum_iterations` further steps, or a `dt K e` that
//! is not finite, is refused (`hertz_patch_line_skew`) and writes nothing.
//! A symmetric line takes the owner's one step, bit for bit, as does a line
//! whose `K` the effort does not reach: a pure spin beside a spinning
//! neighbour, whose gyroscopic coupling lies across the spin axis. A slide
//! beside a strongly precessing neighbour costs more owner steps, each
//! change about `dt |K| |dF/d free|` times the last; where that factor is
//! not below one the steps do not settle and the step is refused.
//! [`HertzPatchState::owner_steps`] counts them.
//!
//! Refused, besides what the source and the owner refuse:
//!
//! - an interface temperature or a step that is not positive and finite, and
//!   released flows or a line that are not finite;
//! - a law whose normal load or radius is not the component's, to the bit
//!   (`hertz_patch_law_geometry`);
//! - a step whose wrench would assist the motion at either of its ends by
//!   more than the owner's energy tolerance: the slip or the spin reverses
//!   inside the step (`hertz_patch_reversal`). These are the module's own
//!   four endpoint powers (slip and spin, at the start and at the end), and
//!   its refusal. The component does not locate the reversal or the arrest.
//!
//! The owner's refusals carry their numbers in [`HertzWrenchError`], which
//! [`TriboError`] has no place for. `Law(e)` is passed on as `e`, and
//! `Invalid` keeps its field. The others are named by field
//! (`hertz_patch_effort_budget`, `hertz_patch_static_boundary`,
//! `hertz_patch_unresolved`, `hertz_patch_cancelled`), and
//! [`HertzWrenchContact::advance_phase`] at the same input gives their
//! numbers. The owner asks `cx` whether cancellation is requested
//! ([`Cx::is_cancel_requested`]) before its solve; a step it cancels is
//! refused as `hertz_patch_cancelled`, since `TriboError` has no
//! cancellation value.
//!
//! Not covered:
//!
//! - Arrest and reversal: a step in which the slip or the spin would pass
//!   through zero is refused, not shortened. `FrictionInterface` locates
//!   stops, on one axis.
//! - The normal load and the radius are constants of the component. No normal
//!   contact law is stepped here.
//! - Rolling, partial slip and elastic shear memory: the owner's rung has
//!   none.
//!
//! Class `NonSmooth`.

use core::array::from_fn;
use core::fmt;
use core::num::NonZeroUsize;

use fs_exec::Cx;
use fs_ports::component::{Class, Component, LoadLine, PortRole, PortSpec, Proposal, Receipt};
use fs_ports::{PortKind, StableId};

use super::contact::{
    HertzWrenchConfig, HertzWrenchContact, HertzWrenchError, HertzWrenchInput, HertzWrenchPhase,
    HertzWrenchStep,
};
use crate::TriboError;

/// The patch a [`HertzPatchInterface`] is held to, and the owner's budgets.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HertzPatchLoad {
    /// Compressive normal load [N], non-negative, held over the run.
    pub normal_n: f64,
    /// Circular Hertz radius [m], positive: the scale of the port's third
    /// coefficient, `a` times the relative spin.
    pub radius_m: f64,
    /// The owner's numerical tolerances and effort budgets.
    pub solve: HertzWrenchConfig,
}

/// State of a [`HertzPatchInterface`]: the relative motion the last step
/// ended with, and how many owner steps it took.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HertzPatchState {
    motion_m_s: [f64; 3],
    owner_steps: usize,
}

impl HertzPatchState {
    /// The relative motion the last step ended with, `[v_x, v_y, a w]`
    /// [m/s], side a against side b, as port `wrench` receives it.
    #[must_use]
    pub const fn motion_m_s(&self) -> [f64; 3] {
        self.motion_m_s
    }

    /// How many steps of the owner the last step took: one on a symmetric
    /// line, more where the line's antisymmetric part moves the step (the
    /// module's "The symmetric part"); none before the first step.
    #[must_use]
    pub const fn owner_steps(&self) -> usize {
        self.owner_steps
    }
}

/// One circular Hertz patch between two bodies, with the law `law_at` returns
/// for the interface temperature of each step. See the module.
pub struct HertzPatchInterface<F> {
    law_at: F,
    load: HertzPatchLoad,
    ports: Vec<PortSpec>,
}

impl<F> fmt::Debug for HertzPatchInterface<F> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("HertzPatchInterface")
            .field("load", &self.load)
            .finish_non_exhaustive()
    }
}

fn refused(field: &'static str) -> TriboError {
    TriboError::InvalidInput { field }
}

/// The owner's response over the step, `dt C_s`, from the line `C`: its
/// symmetric part. See the module's "The symmetric part".
fn symmetric_part(line: &[[f64; 3]; 3], dt_s: f64) -> [[f64; 3]; 3] {
    from_fn(|i| from_fn(|j| dt_s * (0.5 * line[i][j] + 0.5 * line[j][i])))
}

/// The motion the line's antisymmetric part `K` adds over the step at
/// `effort`, `dt K effort` [m]. See the module's "The symmetric part".
fn skew_motion_m(line: &[[f64; 3]; 3], dt_s: f64, effort: [f64; 3]) -> [f64; 3] {
    from_fn(|i| {
        dt_s * (0..3)
            .map(|j| (0.5 * line[i][j] - 0.5 * line[j][i]) * effort[j])
            .sum::<f64>()
    })
}

/// The rounding of a step's energy balance: the port's work `dt effort .
/// mean` is a sum of terms each as large as its own flow makes it, and it is
/// rounded at that size, not at the loss's.
fn balance_rounding_j(
    dt_s: f64,
    released: [f64; 3],
    line: &[[f64; 3]; 3],
    effort: [f64; 3],
    heat_j: f64,
) -> f64 {
    let at_port: f64 = (0..3)
        .map(|i| {
            effort[i].abs()
                * (released[i].abs()
                    + (0..3)
                        .map(|j| line[i][j].abs() * effort[j].abs())
                        .sum::<f64>())
        })
        .sum();
    64.0 * f64::EPSILON * (dt_s * at_port + heat_j) + f64::MIN_POSITIVE
}

/// The step's mean motion on the line at `effort`, which the port carries,
/// and the motion it ends with: twice the owner's mean motion less `start`,
/// or the mean for a coefficient whose row of the line is zero (given).
fn line_motion(
    released: [f64; 3],
    line: &[[f64; 3]; 3],
    effort: [f64; 3],
    step: &HertzWrenchStep,
    start: [f64; 3],
    dt_s: f64,
    radius_m: f64,
) -> Result<([f64; 3], [f64; 3]), TriboError> {
    let mut mean = released;
    for (row, flow) in line.iter().zip(&mut mean) {
        for (c, e) in row.iter().zip(effort) {
            *flow -= c * e;
        }
    }
    let moved = [
        step.slip_displacement_m[0] / dt_s,
        step.slip_displacement_m[1] / dt_s,
        radius_m * step.spin_angle_rad / dt_s,
    ];
    let end: [f64; 3] = from_fn(|i| {
        if line[i].iter().all(|&c| c == 0.0) {
            mean[i]
        } else {
            2.0 * moved[i] - start[i]
        }
    });
    if mean.iter().chain(&end).all(|v| v.is_finite()) {
        Ok((mean, end))
    } else {
        Err(refused("hertz_patch_motion"))
    }
}

/// The owner's refusal, in the crate's error. See the module.
fn owner_refused(error: HertzWrenchError) -> TriboError {
    match error {
        HertzWrenchError::Law(error) => error,
        HertzWrenchError::Invalid(field) => refused(field),
        HertzWrenchError::Cancelled => refused("hertz_patch_cancelled"),
        HertzWrenchError::EffortBudget { .. } => refused("hertz_patch_effort_budget"),
        HertzWrenchError::StaticBoundaryUnresolved { .. } => refused("hertz_patch_static_boundary"),
        HertzWrenchError::Unresolved { .. } => refused("hertz_patch_unresolved"),
    }
}

impl<F> HertzPatchInterface<F>
where
    F: Fn(f64) -> Result<HertzWrenchContact, TriboError>,
{
    /// A patch held to `load`, with the law `law_at` gives at an interface
    /// temperature [K].
    ///
    /// # Errors
    /// A normal load that is negative or not finite, a radius that is not
    /// positive and finite, or budgets the owner refuses.
    pub fn new(law_at: F, load: HertzPatchLoad) -> Result<Self, TriboError> {
        load.solve.validate().map_err(owner_refused)?;
        crate::nonnegative_finite(load.normal_n, "normal_force")?;
        crate::positive_finite(load.radius_m, "contact_radius_m")?;
        let name = |id: &str| StableId::new(id).map_err(|_| refused("port_identifier"));
        let ports = vec![
            PortSpec::power(
                name("wrench")?,
                PortKind::MechanicalForceVelocity,
                PortRole::FlowIn,
                NonZeroUsize::MIN.saturating_add(2),
            )
            .map_err(|_| refused("wrench_port"))?,
            PortSpec::heat(name("heat")?, PortRole::EffortIn).map_err(|_| refused("heat_port"))?,
        ];
        Ok(Self {
            law_at,
            load,
            ports,
        })
    }

    /// The patch this interface is held to.
    #[must_use]
    pub const fn load(&self) -> &HertzPatchLoad {
        &self.load
    }

    /// The state of a pair that starts with the relative motion `motion_m_s`,
    /// `[v_x, v_y, a w]` [m/s], side a against side b (zero is rest): the
    /// state to hand to a kernel with this component.
    ///
    /// # Errors
    /// A motion that is not finite.
    pub fn start(&self, motion_m_s: [f64; 3]) -> Result<HertzPatchState, TriboError> {
        if !motion_m_s.iter().all(|v| v.is_finite()) {
            return Err(refused("motion_m_s"));
        }
        Ok(HertzPatchState {
            motion_m_s,
            owner_steps: 0,
        })
    }

    /// The owner's step on the full line, its effort and how many owner steps
    /// it took. The owner steps on `C_s` with the free motion `dt released`;
    /// where `dt K e` at its effort is not zero, it steps again with
    /// `dt (released - K e)` at the effort `e` of its previous step, until its
    /// force and moment change by no more than its own `force_tolerance_n`
    /// and `moment_tolerance_nm`. Refused when they have not settled after
    /// the owner's `maximum_iterations` further steps. See the module's "The
    /// symmetric part".
    fn step_on_line(
        &self,
        law: &HertzWrenchContact,
        released: [f64; 3],
        line: &[[f64; 3]; 3],
        dt_s: f64,
        phase: HertzWrenchPhase,
        cx: &Cx<'_>,
    ) -> Result<(HertzWrenchStep, [f64; 3], usize), TriboError> {
        let HertzPatchLoad {
            radius_m, solve, ..
        } = self.load;
        let symmetric = symmetric_part(line, dt_s);
        // The first step carries no effort: its free motion is the released
        // flow's, as on a symmetric line.
        let owner = |carried: Option<[f64; 3]>| {
            let mut free = released.map(|flow| dt_s * flow);
            if let Some(carried) = carried {
                for (free, shift) in free.iter_mut().zip(skew_motion_m(line, dt_s, carried)) {
                    *free -= shift;
                }
            }
            let input = HertzWrenchInput {
                free_displacement_m: [free[0], free[1]],
                free_spin_angle_rad: free[2] / radius_m,
                scaled_compliance_m_per_n: symmetric,
            };
            let step = law
                .advance_phase(&input, phase, solve, || cx.is_cancel_requested())
                .map_err(owner_refused)?;
            let effort = [
                -step.force_n[0],
                -step.force_n[1],
                -step.spin_moment_nm / radius_m,
            ];
            Ok::<_, TriboError>((step, effort))
        };
        let (mut step, mut effort) = owner(None)?;
        let shift = skew_motion_m(line, dt_s, effort);
        if shift.iter().all(|&m| m == 0.0) {
            return Ok((step, effort, 1));
        }
        if !shift.iter().all(|m| m.is_finite()) {
            return Err(refused("hertz_patch_line_skew"));
        }
        for repeat in 1..=solve.maximum_iterations {
            let (next, next_effort) = owner(Some(effort))?;
            let settled = (0..2)
                .all(|i| (next.force_n[i] - step.force_n[i]).abs() <= solve.force_tolerance_n)
                && (next.spin_moment_nm - step.spin_moment_nm).abs() <= solve.moment_tolerance_nm;
            (step, effort) = (next, next_effort);
            if settled {
                return Ok((step, effort, 1 + repeat));
            }
        }
        Err(refused("hertz_patch_line_skew"))
    }

    fn law_at_k(&self, interface_temperature_k: f64) -> Result<HertzWrenchContact, TriboError> {
        if !(interface_temperature_k.is_finite() && interface_temperature_k > 0.0) {
            return Err(refused("interface_temperature_k"));
        }
        let law = (self.law_at)(interface_temperature_k)?;
        let parameters = law.parameters();
        if parameters.normal_force_n.to_bits() != self.load.normal_n.to_bits()
            || parameters.contact_radius_m.to_bits() != self.load.radius_m.to_bits()
        {
            return Err(refused("hertz_patch_law_geometry"));
        }
        Ok(law)
    }
}

impl<F> Component for HertzPatchInterface<F>
where
    F: Fn(f64) -> Result<HertzWrenchContact, TriboError>,
{
    type State = HertzPatchState;
    type Error = TriboError;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::NonSmooth
    }

    fn stored_energy_j(&self, _state: &Self::State) -> Result<f64, TriboError> {
        Ok(0.0)
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
        self.propose_on_line(old, dt_s, inputs, &[0.0; 9], new, outputs, cx)
    }
}

impl<F> LoadLine for HertzPatchInterface<F>
where
    F: Fn(f64) -> Result<HertzWrenchContact, TriboError>,
{
    fn line_ports(&self) -> &[usize] {
        &[0]
    }

    fn propose_on_line(
        &self,
        old: &Self::State,
        dt_s: f64,
        inputs: &[f64],
        compliance: &[f64],
        new: &mut Self::State,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, TriboError> {
        let &[flow_x, flow_y, flow_spin, kelvin] = inputs else {
            return Err(refused("hertz_patch_inputs"));
        };
        let &[c00, c01, c02, c10, c11, c12, c20, c21, c22] = compliance else {
            return Err(refused("hertz_patch_inputs"));
        };
        let [effort_x, effort_y, effort_spin, heat_w] = outputs else {
            return Err(refused("hertz_patch_outputs"));
        };
        if !(dt_s.is_finite() && dt_s > 0.0) {
            return Err(refused("dt_s"));
        }
        let released = [flow_x, flow_y, flow_spin];
        let line = [[c00, c01, c02], [c10, c11, c12], [c20, c21, c22]];
        if !released.iter().chain(compliance).all(|v| v.is_finite()) {
            return Err(refused("hertz_patch_inputs"));
        }
        let law = self.law_at_k(kelvin)?;
        let HertzPatchLoad {
            radius_m, solve, ..
        } = self.load;
        // The module's `accepted_phase`: an actual nonzero rate is moving.
        let phase = if old.motion_m_s.iter().any(|&v| v != 0.0) {
            HertzWrenchPhase::Kinetic
        } else {
            HertzWrenchPhase::StaticAdmissible
        };
        let (step, effort, owner_steps) =
            self.step_on_line(&law, released, &line, dt_s, phase, cx)?;
        let start = old.motion_m_s;
        let (mean, end) = line_motion(released, &line, effort, &step, start, dt_s, radius_m)?;
        // The wrench acting on side a, and its power along the slip and the
        // spin at the step's two ends: none may assist the motion.
        let acting = effort.map(|e| -e);
        let assisting_j = dt_s
            * [start, end]
                .into_iter()
                .flat_map(|m| [acting[0] * m[0] + acting[1] * m[1], acting[2] * m[2]])
                .fold(0.0_f64, f64::max);
        if !(assisting_j <= solve.energy_tolerance_j) {
            return Err(refused("hertz_patch_reversal"));
        }
        let work_j = dt_s * (effort[0] * mean[0] + effort[1] * mean[1] + effort[2] * mean[2]);
        let heat_j = step.dissipated_energy_j;
        let loss_j = work_j - heat_j;
        let rounding_j = balance_rounding_j(dt_s, released, &line, effort, heat_j);
        let tolerance_j = solve.energy_tolerance_j + rounding_j;
        if !(work_j.is_finite() && loss_j.abs() <= tolerance_j) {
            return Err(refused("hertz_patch_work_balance"));
        }
        *new = HertzPatchState {
            motion_m_s: end,
            owner_steps,
        };
        *effort_x = effort[0];
        *effort_y = effort[1];
        *effort_spin = effort[2];
        *heat_w = -heat_j / dt_s;
        Ok(Proposal::Step(Receipt {
            stored_change_j: 0.0,
            dissipated_j: heat_j,
            numerical_loss_j: loss_j,
            entropy_generated_j_per_k: heat_j / kelvin,
            energy_tolerance_j: rounding_j,
            residual: loss_j.abs(),
            tolerance: tolerance_j,
        }))
    }
}
