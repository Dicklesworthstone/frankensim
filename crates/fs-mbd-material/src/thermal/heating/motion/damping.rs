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

pub use super::loads::MaterialAttachment as DamperAttachment;
use super::loads::{attachment_arm, attachment_region};
use super::*;
use fs_matdb::InterfaceSystemCard;
use fs_material::phase::EquilibriumPhaseState;
use fs_material::state_point::{
    MaterialPropertySelection, MaterialStatePointError, ResolvedLinearInterface,
};
use fs_material::thermal::interface::{
    AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY, ThermalInterfacePoint,
};
use fs_material::thermal::source::ThermalMaterialError;
use fs_mbd::{ActionReactionImpulseReceipt, apply_equal_and_opposite_impulse_at_body_points};
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::{
    NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport, SolverRunProgress,
};
use fs_thermal::enthalpy::backward_euler::{BackwardEulerStep, StepCoefficients};
use std::cell::{Cell, RefCell};

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
    /// Ordered thermal source, physical context or bulk-state identity refusal.
    Interface(ThermalMaterialError),
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
impl From<ThermalMaterialError> for DampingError {
    fn from(e: ThermalMaterialError) -> Self {
        Self::Interface(e)
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
    /// Signed internal conductive heat from side 0 into side 1 [J]. This
    /// transfer is separate from both generated damping heat and external heat.
    pub contact_heat_to_second_j: f64,
    /// Endpoint A/R'' [W/K], or None for the explicit insulated approximation.
    pub contact_conductance_w_k: Option<f64>,
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
/// reservoirs. An optional source-bound contact exchanges heat between those
/// same networks. Clone retains all accepted enthalpy, configuration, motion,
/// source, contact-coordinate binding and clock state.
#[derive(Clone)]
pub struct DampedRigidPair<'a> {
    bodies: [HeatedRigidBody<'a>; 2],
    damper: MaterialDamper<'a>,
    thermal_contact: Option<ThermalInterfacePoint<'a>>,
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
        Ok(Self {
            bodies,
            damper,
            thermal_contact: None,
        })
    }

    pub const fn bodies(&self) -> &[HeatedRigidBody<'a>; 2] {
        &self.bodies
    }

    pub const fn damper(&self) -> &MaterialDamper<'a> {
        &self.damper
    }

    /// Bind negligible-storage conduction through the SAME original interface.
    ///
    /// The complete damper source must already contain the area-specific
    /// resistance claim. It is admitted by the shared thermal-interface owner
    /// with its exact existing claim pin, never selected a second time. The
    /// original mechanical area supplies A/R''. Both regions are lumped thermal
    /// nodes at their own temperatures; no bulk half-cell resistance is inferred.
    /// The source's existing `T` axis must explicitly be one of those side
    /// temperatures, as selected by `MaterialDamper::temperature_side`.
    ///
    /// Pressure, gap, surface state and every other external coordinate remain
    /// fixed. The caller supplies an interface law valid for this continuing
    /// connection; no collision/contact active set or changing contact area is
    /// inferred from the moving attachment points.
    pub fn with_thermal_contact(
        mut self,
        point: ThermalInterfacePoint<'a>,
    ) -> Result<Self, DampingError> {
        if point.interface.content_hash() != self.damper.card.content_hash()
            || point.temperature_axes[self.damper.temperature_side] != "T"
        {
            return Err(DampingError::Invalid(
                "heat and damping require the same original interface and explicit side T convention",
            ));
        }
        for (i, material) in point.materials.into_iter().enumerate() {
            let region = attachment_region(&self.bodies[i], self.damper.attachments[i])?;
            if material.content_hash() != region.card.content_hash() {
                return Err(DampingError::Invalid(
                    "thermal interface must retain both actual receiving material cards",
                ));
            }
        }
        let states = self.contact_states(&self.bodies[0].thermal, &self.bodies[1].thermal)?;
        let query = point.at(states)?;
        let original = self.damper.source.resolved().query();
        let dynamic: Vec<_> = point
            .temperature_axes
            .into_iter()
            .chain(point.liquid_fraction_axes.into_iter().flatten())
            .chain(point.mean_temperature_axis)
            .collect();
        if original.axes().iter().any(|(axis, value)| {
            !dynamic.contains(&axis.as_str()) && query.axes().get(axis) != Some(value)
        }) {
            return Err(DampingError::Invalid(
                "thermal contact cannot change the interface's fixed source context",
            ));
        }
        let source = require_viscous(
            self.damper
                .source
                .at_query_point(self.damper.card, &query)?,
        )?;
        let resistance = source
            .resolved()
            .property(AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY)
            .ok_or(DampingError::Invalid(
                "the original interface bundle needs a unique thermal resistance",
            ))?;
        let thermal = point.resolve_resistance(
            states,
            MaterialPropertySelection::PinnedByProperty(vec![(
                AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY.into(),
                resistance.answer().receipt.selected,
            )]),
        )?;
        let admitted = &thermal.properties()[0];
        if resistance.requirement().quantity() != admitted.requirement().quantity()
            || resistance.value_si().to_bits() != admitted.value_si().to_bits()
        {
            return Err(DampingError::Invalid(
                "thermal resistance must use the original complete interface bundle",
            ));
        }
        self.damper.source = source;
        self.thermal_contact = Some(point);
        contact_conductance(&self.damper.source)?;
        Ok(self)
    }

    /// Original ordered two-temperature contact binding, if conduction is active.
    pub const fn thermal_contact(&self) -> Option<ThermalInterfacePoint<'a>> {
        self.thermal_contact
    }

    /// Decompose without losing the optional physical interface binding.
    pub fn into_parts(
        self,
    ) -> (
        [HeatedRigidBody<'a>; 2],
        MaterialDamper<'a>,
        Option<ThermalInterfacePoint<'a>>,
    ) {
        (self.bodies, self.damper, self.thermal_contact)
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
                            let (source, contact) = self.resolve_interface(a, b, dt)?;
                            let impulse = self.impulse(&source, dt, controls)?;
                            let qa = number(impulse.heat_j * self.damper.first_heat_fraction)?;
                            let qb = number(impulse.heat_j - qa)?;
                            Ok::<_, DampingError>(StepCoefficients {
                                sources_w: BTreeMap::from([(
                                    self.damper.attachments[1].cell,
                                    number(number(qb + contact.heat_to_second_j)? / dt)?,
                                )]),
                                conductances_w_k: BTreeMap::new(),
                                payload: (source, impulse, [qa, qb], contact),
                            })
                        },
                        checkpoint,
                    )?;
                    let qa = number(b.payload().2[0] - b.payload().3.heat_to_second_j)?;
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
        let (tb, (source, impulse, heat, contact), ab) = trial.accept(
            controls.cell_energy_tolerance_j,
            controls.thermal_energy_tolerance_j,
            checkpoint,
        )?;
        // Ask the original rate owner for authored/reservoir supply without
        // internal damping or contact transfer. Do not derive supply from a
        // closed combined total.
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
                contact_heat_to_second_j: contact.heat_to_second_j,
                contact_conductance_w_k: contact.conductance_w_k,
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

    fn contact_states(
        &self,
        first: &DepositionSimulation<'_>,
        second: &DepositionSimulation<'_>,
    ) -> Result<[EquilibriumPhaseState; 2], DampingError> {
        let state = |i: usize, thermal: &DepositionSimulation<'_>| {
            let attachment = self.damper.attachments[i];
            let region = attachment_region(&self.bodies[i], attachment)?;
            let cell = thermal
                .cells()
                .get(&attachment.cell)
                .ok_or(DampingError::Invalid("thermal contact cell is absent"))?;
            region
                .caloric
                .curve()
                .state_at_specific_enthalpy(cell.specific_h_j_kg)
                .map_err(|e| DampingError::Motion(HeatingError::Caloric(e).into()))
        };
        Ok([state(0, first)?, state(1, second)?])
    }

    fn resolve_interface(
        &self,
        first: &DepositionSimulation<'_>,
        second: &DepositionSimulation<'_>,
        dt: f64,
    ) -> Result<(ResolvedLinearInterface, ContactExchange), DampingError> {
        if let Some(point) = self.thermal_contact {
            let states = self.contact_states(first, second)?;
            let query = point.at(states)?;
            let source = require_viscous(
                self.damper
                    .source
                    .at_query_point(self.damper.card, &query)?,
            )?;
            let exchange = Self::contact_exchange(&source, states, dt)?;
            Ok((source, exchange))
        } else {
            Ok((
                refresh(
                    &self.damper,
                    if self.damper.temperature_side == 0 {
                        first
                    } else {
                        second
                    },
                )?,
                ContactExchange {
                    conductance_w_k: None,
                    heat_to_second_j: 0.0,
                },
            ))
        }
    }

    fn contact_exchange(
        source: &ResolvedLinearInterface,
        states: [EquilibriumPhaseState; 2],
        dt: f64,
    ) -> Result<ContactExchange, DampingError> {
        let conductance = contact_conductance(source)?;
        let difference = number(states[0].temperature_k() - states[1].temperature_k())?;
        let heat = number(dt * number(conductance * difference)?)?;
        Ok(ContactExchange {
            conductance_w_k: Some(conductance),
            heat_to_second_j: heat,
        })
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
        Ok(Dissipation {
            receipt,
            heat_j,
            residual_j,
        })
    }
}

struct Dissipation {
    receipt: ActionReactionImpulseReceipt,
    heat_j: f64,
    residual_j: f64,
}

struct ContactExchange {
    conductance_w_k: Option<f64>,
    heat_to_second_j: f64,
}

fn contact_conductance(source: &ResolvedLinearInterface) -> Result<f64, DampingError> {
    let resistance = source
        .resolved()
        .property(AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY)
        .ok_or(DampingError::Invalid("missing original thermal resistance"))?
        .value_si();
    let conductance = number(source.area_m2() / resistance)?;
    if resistance <= 0.0 || conductance <= 0.0 {
        return Err(DampingError::Invalid(
            "positive representable contact resistance and conductance required",
        ));
    }
    Ok(conductance)
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
    require_viscous(damper.source.at_query_point(damper.card, &point)?)
}

fn require_viscous(
    source: ResolvedLinearInterface,
) -> Result<ResolvedLinearInterface, DampingError> {
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
