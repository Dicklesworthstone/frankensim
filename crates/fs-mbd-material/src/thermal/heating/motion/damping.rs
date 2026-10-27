//! Axial viscous material interfaces between two continuing thermal rigid bodies.
//!
//! Native point impulses remove mechanical work, the original thermal networks
//! receive it once, and the original interface source resolves c(T) at the
//! proposed endpoint. Existing Newton and enthalpy residual owners close that
//! feedback before either physical history is returned.
//!
//! This is a first-order impulse/free-flight/configuration split. The connector
//! has no spring storage, inertia, contact active set or independent heat store.
//! Its fixed area and explicit heat partition describe an admitted lumped viscous
//! interface. Expansion parameter work remains separate; dissipative feedback
//! does not supply the missing thermoelastic free-energy model.

use super::*;
use fs_matdb::InterfaceSystemCard;
use fs_material::state_point::{MaterialStatePointError, ResolvedLinearInterface};
use fs_mbd::{ActionReactionImpulseReceipt, apply_equal_and_opposite_impulse_at_body_points};
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::{
    NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport, SolverRunProgress,
};
use fs_thermal::enthalpy::backward_euler::{BackwardEulerStep, StepCoefficients};
use std::cell::{Cell, RefCell};

/// A material attachment and heat receiver in one body's original coordinates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DamperAttachment {
    /// The same region supplies the surface identity and receiving mass.
    pub cell: CellId,
    /// Reference body coordinates [m]; source expansion moves the attachment.
    pub reference_point_m: fs_ga::Vec3,
}

/// Explicit source, geometry, temperature convention and heat routing.
#[derive(Clone, Debug)]
pub struct MaterialDamper<'a> {
    pub card: &'a InterfaceSystemCard,
    /// A sourced pure Kelvin-Voigt dashpot: k=0, c>0, no constant-loss anchor.
    pub source: ResolvedLinearInterface,
    /// Ordered exactly like the interface's two material surface roles.
    pub attachments: [DamperAttachment; 2],
    /// 0 or 1: whose receiving cell supplies the interface-law temperature.
    /// This is an explicit isothermal-layer approximation, never a mean of T's.
    pub temperature_side: usize,
    /// Fraction of physical dissipation deposited into side 0; side 1 gets the
    /// exact remainder. No effusivity or flash-temperature law is inferred.
    pub first_heat_fraction: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DampingSolve {
    pub maximum_cells: usize,
    /// All source/mechanical proposals, including finite-difference stencils.
    pub maximum_trials: usize,
    pub maximum_newton_iterations: usize,
    pub difference_step_j: f64,
    pub cell_energy_tolerance_j: f64,
    pub thermal_energy_tolerance_j: f64,
    pub total_energy_tolerance_j: f64,
    /// Upper bound in (0,2] on c*dt*(1/m_eff,a+1/m_eff,b), refusing scalar
    /// slip reversal from an excessively coarse midpoint damping substep.
    pub maximum_relaxation_step: f64,
    pub newton: NewtonKrylovConfig,
}

#[derive(Debug)]
pub enum DampingError {
    Motion(HeatedMotionError),
    Material(MaterialStatePointError),
    Thermal(DepositionError),
    Newton(NewtonError),
    NotConverged(NewtonReport),
    Interrupted(SolverRunProgress),
    Invalid(&'static str),
    Budget(&'static str),
}
impl core::fmt::Display for DampingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "material rigid-body damping: {self:?}")
    }
}
impl std::error::Error for DampingError {}
impl From<HeatedMotionError> for DampingError {
    fn from(e: HeatedMotionError) -> Self {
        Self::Motion(e)
    }
}
impl From<PrincipalFrameError> for DampingError {
    fn from(e: PrincipalFrameError) -> Self {
        Self::Motion(e.into())
    }
}
impl From<DynamicsError> for DampingError {
    fn from(e: DynamicsError) -> Self {
        Self::Motion(e.into())
    }
}
impl From<MaterialStatePointError> for DampingError {
    fn from(e: MaterialStatePointError) -> Self {
        Self::Material(e)
    }
}
impl From<DepositionError> for DampingError {
    fn from(e: DepositionError) -> Self {
        Self::Thermal(e)
    }
}
impl From<NewtonError> for DampingError {
    fn from(e: NewtonError) -> Self {
        Self::Newton(e)
    }
}

/// Native impulse, independent heat/storage accounts and fresh source evaluation.
#[derive(Debug)]
pub struct DampingAudit {
    pub impulse: ActionReactionImpulseReceipt,
    pub source: ResolvedLinearInterface,
    pub thermal: [EnergyAudit; 2],
    pub heat_by_body_j: [f64; 2],
    /// c*dt*v_mid^2, independently checked against native impulse work.
    pub dissipated_heat_j: f64,
    pub impulse_energy_residual_j: f64,
    pub external_heat_j: f64,
    /// Quasi-static inertia work, never deposited as dissipation/caloric heat.
    pub parameter_work_j: [f64; 2],
    /// Signed native free-flight drift, never supplied to the thermal solver.
    pub numerical_energy_change_j: [f64; 2],
    pub combined_residual_j: f64,
    pub maximum_cell_residual_j: f64,
    pub nonlinear: NewtonReport,
    pub trials: usize,
}

/// Two accepted thermal/mechanical histories and one source-bound viscous link.
///
/// Each original thermal network retains its material regions, conduction and
/// reservoirs. This exchanges mechanical work and dissipative heat; it adds no
/// direct heat-conduction link between the networks. Clone retains all accepted
/// enthalpy, configuration, motion, source and clock state.
#[derive(Clone)]
pub struct DampedRigidPair<'a> {
    bodies: [HeatedRigidBody<'a>; 2],
    damper: MaterialDamper<'a>,
}

impl<'a> DampedRigidPair<'a> {
    pub fn new(
        bodies: [HeatedRigidBody<'a>; 2],
        mut damper: MaterialDamper<'a>,
    ) -> Result<Self, DampingError> {
        if bodies[0].thermal.time_s() != bodies[1].thermal.time_s()
            || damper.temperature_side > 1
            || !damper.first_heat_fraction.is_finite()
            || !(0.0..=1.0).contains(&damper.first_heat_fraction)
        {
            return Err(DampingError::Invalid(
                "equal clocks and explicit valid temperature/heat routing required",
            ));
        }
        for (body, attachment) in bodies.iter().zip(damper.attachments) {
            attachment_region(body, attachment)?;
        }
        let a = attachment_region(&bodies[0], damper.attachments[0])?;
        let b = attachment_region(&bodies[1], damper.attachments[1])?;
        damper
            .source
            .require_surface_materials(a.card.id(), b.card.id())?;
        damper.source = refresh(&damper, &bodies[damper.temperature_side].thermal)?;
        Ok(Self { bodies, damper })
    }

    pub const fn bodies(&self) -> &[HeatedRigidBody<'a>; 2] {
        &self.bodies
    }

    pub const fn damper(&self) -> &MaterialDamper<'a> {
        &self.damper
    }

    pub fn into_parts(self) -> ([HeatedRigidBody<'a>; 2], MaterialDamper<'a>) {
        (self.bodies, self.damper)
    }

    /// Solve one bounded interval; continue from the returned pair.
    ///
    /// Unknowns are cell-energy increments [J]: body 0's sorted CellId order,
    /// then body 1's. Every trial restarts from the same accepted bodies.
    /// Pending start/interior events refuse; process each body's current-time
    /// events before reconstructing the pair. Target events remain pending.
    ///
    /// The axial central impulse uses actual starting attachments and proposed
    /// endpoint c(T). Native free flight advances both bodies; source expansion
    /// then updates full inertia at conserved world momenta. This split is first
    /// order and requires time refinement independently of Newton convergence.
    pub fn step_to(
        &self,
        target_s: f64,
        controls: DampingSolve,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<(Self, DampingAudit), DampingError> {
        let cancel = RefCell::new(cancelled);
        let checkpoint = || {
            if (*cancel.borrow_mut())() {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        };
        checkpoint()?;
        controls.newton.validate()?;
        let count = self
            .bodies
            .iter()
            .map(|b| b.thermal.cells().len())
            .sum::<usize>();
        let dt = target_s - self.bodies[0].thermal.time_s();
        if controls.maximum_cells == 0
            || count > controls.maximum_cells
            || initial_energy_increments_j.len() != count
            || controls.maximum_trials == 0
            || controls.maximum_newton_iterations == 0
            || !controls.maximum_relaxation_step.is_finite()
            || controls.maximum_relaxation_step <= 0.0
            || controls.maximum_relaxation_step > 2.0
            || [
                controls.cell_energy_tolerance_j,
                controls.thermal_energy_tolerance_j,
                controls.total_energy_tolerance_j,
            ]
            .iter()
            .any(|v| !v.is_finite() || *v < 0.0)
            || self
                .bodies
                .iter()
                .any(|b| dt > b.motion.controls.maximum_step_s)
        {
            return Err(DampingError::Invalid(
                "invalid paired thermal unknowns, step, budget or energy controls",
            ));
        }
        let first = BackwardEulerStep::new(&self.bodies[0].thermal, target_s, checkpoint)?;
        let second = BackwardEulerStep::new(&self.bodies[1].thermal, target_s, checkpoint)?;
        let split = first.cell_ids().len();
        let trials = Cell::new(0usize);
        let evaluate = |x: &[f64]| {
            if x.len() != count {
                return Err(DampingError::Invalid("wrong energy vector length"));
            }
            first.evaluate(
                &x[..split],
                |a| {
                    let b = second.evaluate(
                        &x[split..],
                        |b| {
                            checkpoint()?;
                            if trials.get() >= controls.maximum_trials {
                                return Err(DampingError::Budget("interface proposals"));
                            }
                            trials.set(trials.get() + 1);
                            let source = refresh(
                                &self.damper,
                                if self.damper.temperature_side == 0 { a } else { b },
                            )?;
                            let impulse = self.impulse(&source, dt, controls)?;
                            let qa = number(impulse.heat_j * self.damper.first_heat_fraction)?;
                            let qb = number(impulse.heat_j - qa)?;
                            Ok::<_, DampingError>(StepCoefficients {
                                sources_w: BTreeMap::from([(
                                    self.damper.attachments[1].cell,
                                    number(qb / dt)?,
                                )]),
                                conductances_w_k: BTreeMap::new(),
                                payload: (source, impulse, [qa, qb]),
                            })
                        },
                        checkpoint,
                    )?;
                    let qa = b.payload().2[0];
                    Ok::<_, DampingError>(StepCoefficients {
                        sources_w: BTreeMap::from([(
                            self.damper.attachments[0].cell,
                            number(qa / dt)?,
                        )]),
                        conductances_w_k: BTreeMap::new(),
                        payload: b,
                    })
                },
                checkpoint,
            )
        };
        let problem =
            FallibleFiniteDifference::new(count, controls.difference_step_j, |x: &[f64]| {
                let trial = evaluate(x)?;
                let mut residual = trial.residual_j().to_vec();
                residual.extend_from_slice(trial.payload().residual_j());
                Ok::<_, DampingError>(residual)
            })
            .map_err(numerical)?;
        let construction = NewtonKrylovState::new(
            &problem,
            initial_energy_increments_j.to_vec(),
            controls.newton,
        );
        if let Some(e) = problem.take_failure() {
            return Err(numerical(e));
        }
        let mut solver = construction?;
        let result =
            solver.run_controlled(&problem, controls.maximum_newton_iterations, &mut || {
                (*cancel.borrow_mut())() || problem.has_failed()
            });
        if let Some(e) = problem.take_failure() {
            return Err(numerical(e));
        }
        checkpoint()?;
        if result.progress != SolverRunProgress::Complete {
            return Err(DampingError::Interrupted(result.progress));
        }
        if !result.report.converged {
            return Err(DampingError::NotConverged(result.report));
        }
        let trial = evaluate(&solver.x)?;
        let maximum_cell_residual_j = trial
            .maximum_cell_residual_j()
            .max(trial.payload().maximum_cell_residual_j());
        let (ta, trial, aa) = trial.accept(
            controls.cell_energy_tolerance_j,
            controls.thermal_energy_tolerance_j,
            checkpoint,
        )?;
        let (tb, (source, impulse, heat), ab) = trial.accept(
            controls.cell_energy_tolerance_j,
            controls.thermal_energy_tolerance_j,
            checkpoint,
        )?;
        // Ask the original rate owner for authored/reservoir supply without
        // internal damping. Do not derive supply from a closed combined total.
        let mut external = 0.0;
        for thermal in [&ta, &tb] {
            external = number(
                external
                    + number(
                        dt * thermal.external_step_power_w(
                            &BTreeMap::new(),
                            &BTreeMap::new(),
                            checkpoint,
                        )?,
                    )?,
            )?;
        }
        if number(aa.external_heat_j + ab.external_heat_j - external - impulse.heat_j)?.abs()
            > controls.total_energy_tolerance_j
        {
            return Err(DampingError::Invalid(
                "deposited damping heat differs from the independent thermal supply account",
            ));
        }
        let mut next = self.clone();
        let states = [
            impulse.receipt.first.state_after,
            impulse.receipt.second.state_after,
        ];
        for (i, thermal) in [ta, tb].into_iter().enumerate() {
            checkpoint()?;
            let motion = &mut next.bodies[i].motion;
            motion.begin()?;
            motion.state = states[i];
            motion.step(
                self.bodies[i].thermal.time_s(),
                dt,
                &BTreeMap::new(),
                &thermal,
                &mut || checkpoint(),
            )?;
            next.bodies[i].thermal = thermal;
        }
        let parameter = next.bodies.each_ref().map(|b| b.motion.work.parameter_j);
        let numerical_change = next.bodies.each_ref().map(|b| b.motion.work.numerical_j);
        let mut combined =
            number(aa.stored_energy_change_j + ab.stored_energy_change_j - external)?;
        for i in 0..2 {
            combined = number(
                combined
                    + mechanical_change(
                        self.bodies[i].diagnostics()?,
                        next.bodies[i].diagnostics()?,
                    )?
                    - parameter[i]
                    - numerical_change[i],
            )?;
        }
        if combined.abs() > controls.total_energy_tolerance_j {
            return Err(DampingError::Invalid(
                "mechanical and caloric storage do not close the physical heat/work accounts",
            ));
        }
        next.damper.source = source.clone();
        checkpoint()?;
        Ok((
            next,
            DampingAudit {
                impulse: impulse.receipt,
                source,
                thermal: [aa, ab],
                heat_by_body_j: heat,
                dissipated_heat_j: impulse.heat_j,
                impulse_energy_residual_j: impulse.residual_j,
                external_heat_j: external,
                parameter_work_j: parameter,
                numerical_energy_change_j: numerical_change,
                combined_residual_j: combined,
                maximum_cell_residual_j,
                nonlinear: result.report,
                trials: trials.get(),
            },
        ))
    }

    fn impulse(
        &self,
        source: &ResolvedLinearInterface,
        dt: f64,
        controls: DampingSolve,
    ) -> Result<Dissipation, DampingError> {
        let arms = [
            attachment_arm(&self.bodies[0], self.damper.attachments[0])?,
            attachment_arm(&self.bodies[1], self.damper.attachments[1])?,
        ];
        let a = &self.bodies[0];
        let b = &self.bodies[1];
        let ka = a
            .state()
            .point_kinematics(a.principal_frame().properties(), arms[0])?;
        let kb = b
            .state()
            .point_kinematics(b.principal_frame().properties(), arms[1])?;
        let separation = kb.point_world.sub(ka.point_world);
        let ea = a.state().directional_effective_mass(
            a.principal_frame().properties(),
            arms[0],
            separation,
        )?;
        let eb = b.state().directional_effective_mass(
            b.principal_frame().properties(),
            arms[1],
            separation,
        )?;
        let inverse_mass = number(ea.inverse_mass_kg_inverse + eb.inverse_mass_kg_inverse)?;
        let c = source.damping_n_s_m();
        let rate = number(c * inverse_mass)?;
        source.require_angular_frequency_rad_s(0.0)?;
        source.require_angular_frequency_rad_s(rate)?;
        let relaxation = number(rate * dt)?;
        if relaxation > controls.maximum_relaxation_step {
            return Err(DampingError::Budget(
                "refine physical timestep for the resolved interface relaxation rate",
            ));
        }
        let velocity = number(
            kb.point_velocity_world
                .sub(ka.point_velocity_world)
                .dot(ea.unit_direction_world),
        )?;
        let midpoint_velocity = number(velocity / (1.0 + 0.5 * relaxation))?;
        let magnitude = number(-c * dt * midpoint_velocity)?;
        let receipt = apply_equal_and_opposite_impulse_at_body_points(
            a.state(),
            a.principal_frame().properties(),
            arms[0],
            b.state(),
            b.principal_frame().properties(),
            arms[1],
            ea.unit_direction_world.scale(magnitude),
        )?;
        let heat_j = number(c * dt * midpoint_velocity * midpoint_velocity)?;
        let residual_j = number(receipt.kinetic_energy_change_j + heat_j)?;
        if heat_j < 0.0
            || residual_j.abs() > controls.total_energy_tolerance_j
            || receipt.work_energy_residual_j.abs() > controls.total_energy_tolerance_j
            || number(receipt.impulse_work_j + heat_j)?.abs() > controls.total_energy_tolerance_j
        {
            return Err(DampingError::Invalid(
                "native impulse loss and physical dashpot heating disagree",
            ));
        }
        Ok(Dissipation { receipt, heat_j, residual_j })
    }
}

struct Dissipation {
    receipt: ActionReactionImpulseReceipt,
    heat_j: f64,
    residual_j: f64,
}

fn attachment_region<'b, 'a>(
    body: &'b HeatedRigidBody<'a>,
    attachment: DamperAttachment,
) -> Result<&'b ThermalRegion<'a>, DampingError> {
    if ![
        attachment.reference_point_m.x,
        attachment.reference_point_m.y,
        attachment.reference_point_m.z,
    ]
    .into_iter()
    .all(f64::is_finite)
    {
        return Err(DampingError::Invalid("nonfinite material attachment"));
    }
    body.motion
        .expansion
        .regions
        .iter()
        .find(|r| r.cell_id == attachment.cell)
        .ok_or(DampingError::Invalid(
            "attachment must name an original material/thermal region",
        ))
}

fn attachment_arm(
    body: &HeatedRigidBody<'_>,
    attachment: DamperAttachment,
) -> Result<fs_mbd::Vec3, DampingError> {
    let binding = attachment_region(body, attachment)?;
    let index = body
        .reference()
        .regions()
        .iter()
        .position(|r| r.key() == binding.key)
        .ok_or(DampingError::Invalid("missing original attachment region"))?;
    let region = &body.configuration().regions()[index];
    let point = region.fixed_point_m()
        + (attachment.reference_point_m - region.fixed_point_m())
            .scale(1.0 + region.expansion().free_linear_strain());
    Ok(body.principal_frame().arm_for_body_point(point)?)
}

fn refresh(
    damper: &MaterialDamper<'_>,
    thermal: &DepositionSimulation<'_>,
) -> Result<ResolvedLinearInterface, DampingError> {
    let mut point = damper.source.resolved().query();
    if !point.axes().contains_key("T") {
        return Err(DampingError::Invalid(
            "interface temperature requires the original absolute T coordinate",
        ));
    }
    point = point
        .with(
            "T",
            thermal.temperature_k(damper.attachments[damper.temperature_side].cell)?,
        )
        .map_err(|_| DampingError::Invalid("invalid interface temperature coordinate"))?;
    let source = damper.source.at_query_point(damper.card, &point)?;
    if source.stiffness_n_m() != 0.0
        || source.damping_n_s_m() <= 0.0
        || source.constant_loss_reference_omega_rad_s().is_some()
    {
        return Err(DampingError::Invalid(
            "a pure sourced viscous interface is required; spring/constant-loss storage is unsupported",
        ));
    }
    Ok(source)
}

fn number(value: f64) -> Result<f64, DampingError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(DampingError::Invalid(
            "unrepresentable physical damping arithmetic",
        ))
    }
}

fn numerical(error: FiniteDifferenceError<DampingError>) -> DampingError {
    match error {
        FiniteDifferenceError::Source(e) => e,
        FiniteDifferenceError::Invalid(s) => DampingError::Invalid(s),
    }
}
