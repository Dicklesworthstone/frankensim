//! One native momentum/contact/caloric solve for a whole retained contact island.
//!
//! The unknowns are the original cell energy increments, a massless temperature
//! for every original interface, and all shear/spin impulses. Every trial first
//! applies ALL impulses to the original bodies, then evaluates every retained
//! return map at that same midpoint motion and caloric endpoint. Shared inertia,
//! shared cell heat, and affine contact free energy therefore enter once.
//!
//! Only physical microslip heats the common nodes. Each original surface keeps
//! its own reversible -T ds exchange; numerical return-map loss stays separate.
//! This is a first-order constitutive step at fixed original poses and supplied
//! positive normal patches, not normal dynamics or moving thermal geometry.

use super::*;
use crate::friction_heat::FrictionHeatTolerances;
use crate::friction_heat::network::{
    FrictionHeatNetworkInput, FrictionHeatNetworkStep, InterfaceHeatInput, InterfaceHeatReport,
};

/// One accepted original contact and its physical native body binding.
pub struct RigidCaloricNetworkContact<'s, 'a> {
    pub contact: &'s CaloricPartialSlip<'a>,
    /// Original A/B surface order, including the existing thermal link order.
    pub bodies: [usize; 2],
    pub arms_body_m: [Vec3; 2],
    pub frame: TangentFrame,
    pub creepage_reference_speed_mps: f64,
    pub ownership: &'s GeneralizedWorkOwnership,
    /// Shear impulses and torsional impulse divided by the declared length.
    pub impulse_guess_n_s: [f64; 3],
    pub interface_temperature_guess_k: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RigidCaloricNetworkSolve {
    /// These are joint-trial and joint-Newton budgets for the entire island.
    pub contact: ResolvedCaloricContactSolve,
    pub maximum_bodies: usize,
    pub maximum_contacts: usize,
    /// Converts impulse coordinates and residuals [N s] to solver energies [J].
    /// This is numerical conditioning, never an additional physical velocity.
    pub impulse_coordinate_speed_mps: f64,
    pub torsional_coordinate_length_m: f64,
    pub geometry_tolerance_m: f64,
    pub linear_momentum_tolerance_n_s: f64,
    pub angular_momentum_tolerance_nm_s: f64,
    pub mechanical_energy_tolerance_j: f64,
}

pub struct CaloricNetworkInterfaceSource {
    pub link_id: u64,
    pub conductivity_w_m_k: [f64; 2],
}

pub struct RigidCaloricNetworkSource<S> {
    /// Exactly the callback's ascending original link ID order.
    pub interfaces: Vec<CaloricNetworkInterfaceSource>,
    /// OTHER sources only. Contact heat, reversible exchange and the bound
    /// series links have exclusive owners inside this step.
    pub transport: StepCoefficients<S>,
}

pub struct RigidCaloricInterfaceReport {
    pub heat: InterfaceHeatReport,
    pub contact: PartialSlipStep,
    pub material: Option<SlidingMaterialSample>,
    pub thermodynamics: ContactThermodynamics,
    pub contact_internal_energy_change_j: f64,
    pub reversible_heat_j: [f64; 2],
    pub physical_heat_j: [f64; 2],
    pub deposited_heat_j: [f64; 2],
    pub conduction_heat_j: [f64; 2],
    pub applied_impulse_on_a_n_s: Vec3,
    pub applied_spin_impulse_on_a_nm_s: Vec3,
    pub contact_impulse_residual_n_s: Vec3,
    pub contact_spin_impulse_residual_nm_s: Vec3,
    pub midpoint_relative_velocity_mps: [f64; 2],
    pub midpoint_relative_spin_rad_s: f64,
}

pub struct RigidCaloricNetworkReport<S> {
    /// Ascending original thermal link ID, matching the returned contacts.
    pub interfaces: Vec<RigidCaloricInterfaceReport>,
    /// The original finite network's storage is counted exactly once.
    pub thermal: EnergyAudit,
    pub external_heat_j: f64,
    pub source: S,
    pub nonlinear: NewtonReport,
    pub trials: usize,
    pub maximum_cell_residual_j: f64,
    pub body_momentum_residual_n_s: Vec<Vec3>,
    pub body_angular_momentum_residual_nm_s: Vec<Vec3>,
    pub linear_momentum_change_n_s: Vec3,
    pub world_angular_momentum_change_nm_s: Vec3,
    pub kinetic_energy_change_j: f64,
    pub mechanical_energy_residual_j: f64,
    /// Delta(K + sum U_contact + H_caloric) + sum D_numerical - Q_external.
    pub total_energy_residual_j: f64,
}

/// Immutable original owners in; all accepted owners out together.
pub struct RigidCaloricNetworkCandidate<'a, 't, S> {
    bodies: Vec<RigidBodyState>,
    contacts: Vec<CaloricPartialSlip<'a>>,
    thermal: DepositionSimulation<'t>,
    report: RigidCaloricNetworkReport<S>,
}
impl<'a, 't, S> RigidCaloricNetworkCandidate<'a, 't, S> {
    pub fn bodies(&self) -> &[RigidBodyState] {
        &self.bodies
    }
    pub fn contacts(&self) -> &[CaloricPartialSlip<'a>] {
        &self.contacts
    }
    pub const fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub const fn report(&self) -> &RigidCaloricNetworkReport<S> {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        Vec<RigidBodyState>,
        Vec<CaloricPartialSlip<'a>>,
        DepositionSimulation<'t>,
        RigidCaloricNetworkReport<S>,
    ) {
        (self.bodies, self.contacts, self.thermal, self.report)
    }
}

impl<E> From<DepositionError> for RigidCaloricContactError<E> {
    fn from(error: DepositionError) -> Self {
        Self::Contact(error.into())
    }
}

struct Staged<S> {
    bodies: Vec<RigidBodyState>,
    mechanical: Vec<MechanicalEndpoint>,
    impulses: Vec<Vec3>,
    couples: Vec<Vec3>,
    slip: Vec<[f64; 2]>,
    spin: Vec<f64>,
    impulse_residual_j: Vec<f64>,
    external_heat_j: f64,
    source: S,
}

/// Resolve every native contact, source and thermal cell at one common endpoint.
///
/// Caller enumeration does not choose a mechanical or thermal update order:
/// contacts are canonicalized by original link identity. No intermediate owner
/// is published. The callback is evaluated once per joint trial, with the whole
/// proposed original network and all common-interface temperatures. Affine
/// stiffness and any pinned friction source retain their original identities,
/// temperature axes, domains, accumulated travel and regional history.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn step_rigid_caloric_network<'a, 't, S, E>(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    contacts: &[RigidCaloricNetworkContact<'_, 'a>],
    thermal: &DepositionSimulation<'t>,
    target_s: f64,
    controls: RigidCaloricNetworkSolve,
    energy_guess_j: &[f64],
    resolve: impl Fn(
        &DepositionSimulation<'t>,
        &[u64],
        &[f64],
    ) -> Result<RigidCaloricNetworkSource<S>, E>,
    cancelled: impl FnMut() -> bool,
) -> Result<RigidCaloricNetworkCandidate<'a, 't, S>, RigidCaloricContactError<E>> {
    let cancel = RefCell::new(cancelled);
    let poll = || {
        if (*cancel.borrow_mut())() {
            Err(DepositionError::Cancelled)
        } else {
            Ok(())
        }
    };
    poll()?;
    let cc = controls.contact.contact;
    cc.newton.validate()?;
    let dt = target_s - thermal.time_s();
    if bodies.is_empty()
        || bodies.len() != properties.len()
        || bodies.len() > controls.maximum_bodies
        || contacts.is_empty()
        || contacts.len() > controls.maximum_contacts
        || thermal.cells().len() > cc.maximum_cells
        || energy_guess_j.len() != thermal.cells().len()
        || energy_guess_j.iter().any(|x| !x.is_finite())
        || cc.maximum_contact_trials == 0
        || cc.maximum_newton_iterations == 0
        || !dt.is_finite()
        || dt <= 0.0
        || [
            cc.difference_step_j,
            controls.contact.temperature_scale_j_per_k,
            controls.impulse_coordinate_speed_mps,
            controls.torsional_coordinate_length_m,
            controls.geometry_tolerance_m,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x <= 0.0)
        || [
            cc.cell_energy_tolerance_j,
            cc.thermal_energy_tolerance_j,
            cc.total_energy_tolerance_j,
            controls.contact.interface_energy_tolerance_j,
            controls.linear_momentum_tolerance_n_s,
            controls.angular_momentum_tolerance_nm_s,
            controls.mechanical_energy_tolerance_j,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x < 0.0)
    {
        return Err(RigidCaloricContactError::Invalid(
            "bounded synchronized native caloric island and finite controls required",
        ));
    }
    for input in contacts {
        poll()?;
        let contact = input.contact;
        if input.bodies[0] == input.bodies[1]
            || input.bodies.iter().any(|i| *i >= bodies.len())
            || contact.time_s != thermal.time_s()
            || temperatures(thermal, contact.curves, contact.cells)? != contact.temperature_k
            || !input.interface_temperature_guess_k.is_finite()
            || input.interface_temperature_guess_k <= 0.0
            || !input.creepage_reference_speed_mps.is_finite()
            || input.creepage_reference_speed_mps <= 0.0
            || input.impulse_guess_n_s.iter().any(|x| !x.is_finite())
            || (contact.material_friction.is_some() && input.impulse_guess_n_s[2] != 0.0)
        {
            return Err(RigidCaloricContactError::Invalid(
                "each original contact needs valid ordered bodies, clock, caloric objects and guesses",
            ));
        }
        contact.admit_rigid_geometry(
            &input.bodies.map(|i| bodies[i]),
            input.bodies.map(|i| properties[i]),
            input.arms_body_m,
            thermal,
            input.frame,
            controls.geometry_tolerance_m,
        )?;
    }
    let mut ordered: Vec<_> = contacts.iter().collect();
    ordered.sort_by_key(|input| input.contact.resolved_heat_route().unwrap().link_id());
    let mut patches = std::collections::BTreeSet::new();
    for input in &ordered {
        if !patches.insert(input.contact.checkpoint.normal_patch().patch_id()) {
            return Err(RigidCaloricContactError::Invalid(
                "one retained history and heat owner per physical contact required",
            ));
        }
    }
    let routes: Vec<_> = ordered
        .iter()
        .map(|input| input.contact.resolved_heat_route().unwrap())
        .collect();
    let heat_step =
        FrictionHeatNetworkStep::new(thermal, target_s, &routes, controls.maximum_contacts, poll)?;
    let cells = heat_step.cell_ids().len();
    let heat_dimension =
        cells
            .checked_add(contacts.len())
            .ok_or(RigidCaloricContactError::Invalid(
                "native caloric coordinate count overflow",
            ))?;
    let mut dimension = heat_dimension;
    let mut coordinates = Vec::with_capacity(contacts.len());
    let mut before = Vec::with_capacity(contacts.len());
    for input in &ordered {
        let count = if input.contact.material_friction.is_some() {
            2
        } else {
            3
        };
        coordinates.push((dimension, count));
        dimension = dimension
            .checked_add(count)
            .ok_or(RigidCaloricContactError::Invalid(
                "native caloric coordinate count overflow",
            ))?;
        before.push([
            bodies[input.bodies[0]]
                .point_kinematics(properties[input.bodies[0]], input.arms_body_m[0])?,
            bodies[input.bodies[1]]
                .point_kinematics(properties[input.bodies[1]], input.arms_body_m[1])?,
        ]);
    }
    let speed_scale = controls.impulse_coordinate_speed_mps;
    let length = controls.torsional_coordinate_length_m;
    let trials = Cell::new(0usize);
    let evaluate = |x: &[f64]| {
        poll()?;
        if trials.get() >= cc.maximum_contact_trials {
            return Err(RigidCaloricContactError::MotionBudget);
        }
        trials.set(trials.get() + 1);
        if x.len() != dimension || x.iter().any(|x| !x.is_finite()) {
            return Err(RigidCaloricContactError::Invalid(
                "complete finite joint iterate required",
            ));
        }
        let mut next = bodies.to_vec();
        let mut impulses = Vec::with_capacity(ordered.len());
        let mut couples = Vec::with_capacity(ordered.len());
        for (i, input) in ordered.iter().enumerate() {
            poll()?;
            let (start, count) = coordinates[i];
            let impulse = vector(input.frame.longitudinal())
                .scale(x[start] / speed_scale)
                .add(vector(input.frame.lateral()).scale(x[start + 1] / speed_scale));
            let couple = vector(input.frame.normal()).scale(if count == 3 {
                length * x[start + 2] / speed_scale
            } else {
                0.0
            });
            for side in 0..2 {
                let body = input.bodies[side];
                let sign = if side == 0 { 1.0 } else { -1.0 };
                let kicked = next[body]
                    .apply_impulse_at_body_point(
                        properties[body],
                        input.arms_body_m[side],
                        impulse.scale(sign),
                    )?
                    .state_after;
                next[body] = RigidBodyState::new(
                    kicked.pose(),
                    kicked.linear_momentum_world(),
                    kicked.angular_momentum_body().add(
                        kicked
                            .pose()
                            .orientation()
                            .rotate_world_to_body(couple.scale(sign)),
                    ),
                )?;
            }
            impulses.push(impulse);
            couples.push(couple);
        }
        let temperatures: Vec<_> = ordered
            .iter()
            .enumerate()
            .map(|(i, input)| {
                input.interface_temperature_guess_k
                    + x[cells + i] / controls.contact.temperature_scale_j_per_k
            })
            .collect();
        heat_step.evaluate(
            &x[..cells],
            &temperatures,
            |endpoint, temperatures| {
                let RigidCaloricNetworkSource {
                    interfaces,
                    transport,
                } = resolve(endpoint, heat_step.link_ids(), temperatures).map_err(|e| {
                    RigidCaloricContactError::Contact(CaloricContactError::Source(e))
                })?;
                if interfaces.len() != ordered.len()
                    || interfaces
                        .iter()
                        .zip(heat_step.link_ids())
                        .any(|(source, id)| source.link_id != *id)
                {
                    return Err(RigidCaloricContactError::Invalid(
                        "resolved interfaces must retain canonical original link identities",
                    ));
                }
                let StepCoefficients {
                    mut sources_w,
                    conductances_w_k,
                    payload,
                } = transport;
                // Measure actual outside supply before adding any contact's
                // reversible exchange; shared cells must not repeat this sum.
                let external_heat_j = finite(
                    dt * endpoint.external_step_power_w(&sources_w, &conductances_w_k, poll)?,
                )?;
                let mut mechanical = Vec::with_capacity(ordered.len());
                let mut slip = Vec::with_capacity(ordered.len());
                let mut spin = Vec::with_capacity(ordered.len());
                let mut heat = Vec::with_capacity(ordered.len());
                let mut impulse_residual_j = Vec::with_capacity(dimension - heat_dimension);
                for (i, input) in ordered.iter().enumerate() {
                    poll()?;
                    let after = [
                        next[input.bodies[0]]
                            .point_kinematics(properties[input.bodies[0]], input.arms_body_m[0])?,
                        next[input.bodies[1]]
                            .point_kinematics(properties[input.bodies[1]], input.arms_body_m[1])?,
                    ];
                    let relative = before[i][0]
                        .point_velocity_world
                        .add(after[0].point_velocity_world)
                        .sub(before[i][1].point_velocity_world)
                        .sub(after[1].point_velocity_world)
                        .scale(0.5);
                    let relative_spin = before[i][0]
                        .angular_velocity_world
                        .add(after[0].angular_velocity_world)
                        .sub(before[i][1].angular_velocity_world)
                        .sub(after[1].angular_velocity_world)
                        .scale(0.5)
                        .dot(vector(input.frame.normal()));
                    let common_spin = before[i][0]
                        .angular_velocity_world
                        .add(after[0].angular_velocity_world)
                        .add(before[i][1].angular_velocity_world)
                        .add(after[1].angular_velocity_world)
                        .scale(0.25)
                        .dot(vector(input.frame.normal()));
                    if input.contact.material_friction.is_some() && relative_spin != 0.0 {
                        return Err(RigidCaloricContactError::Invalid(
                            "material friction has no torsional-speed mapping",
                        ));
                    }
                    let local_slip = [
                        relative.dot(vector(input.frame.longitudinal())),
                        relative.dot(vector(input.frame.lateral())),
                    ];
                    let m = input.contact.mechanical_endpoint(
                        endpoint,
                        input.frame,
                        PartialSlipKinematics {
                            creepage: local_slip.map(|v| v / input.creepage_reference_speed_mps),
                            rolling_speed_mps: input.creepage_reference_speed_mps,
                            torsional_spin_rad_per_s: relative_spin,
                            contact_spin_rad_per_s: common_spin,
                            dt_s: dt,
                        },
                        input.ownership,
                        cc,
                    )?;
                    let (start, count) = coordinates[i];
                    let expected = [
                        dt * m.contact.tangent_force_n[0],
                        dt * m.contact.tangent_force_n[1],
                        dt * m.contact.torsional_moment_nm / length,
                    ];
                    for axis in 0..count {
                        impulse_residual_j
                            .push(finite(x[start + axis] - speed_scale * expected[axis])?);
                    }
                    for (id, heat_j) in input.contact.cells.into_iter().zip(m.reversible_heat_j) {
                        let power = finite(heat_j / dt)?;
                        if heat_j != 0.0 && power == 0.0 {
                            return Err(RigidCaloricContactError::Invalid(
                                "reversible contact heat vanished in power conversion",
                            ));
                        }
                        let source = sources_w.entry(id).or_insert(0.0);
                        *source = finite(*source + power)?;
                    }
                    heat.push(InterfaceHeatInput {
                        dissipated_j: m.contact.dissipation.heat_j,
                        conductivity_w_m_k: interfaces[i].conductivity_w_m_k,
                    });
                    mechanical.push(m);
                    slip.push(local_slip);
                    spin.push(relative_spin);
                }
                Ok(FrictionHeatNetworkInput {
                    interfaces: heat,
                    transport: StepCoefficients {
                        sources_w,
                        conductances_w_k,
                        payload: Staged {
                            bodies: next,
                            mechanical,
                            impulses,
                            couples,
                            slip,
                            spin,
                            impulse_residual_j,
                            external_heat_j,
                            source: payload,
                        },
                    },
                })
            },
            poll,
        )
    };
    let problem = FallibleFiniteDifference::new(dimension, cc.difference_step_j, |x: &[f64]| {
        evaluate(x).map(|trial| {
            let mut residual = trial.residual_j().to_vec();
            residual.extend_from_slice(&trial.payload().impulse_residual_j);
            residual
        })
    })
    .map_err(rigid_difference_error)?;
    let mut initial = energy_guess_j.to_vec();
    initial.resize(heat_dimension, 0.0);
    for (i, input) in ordered.iter().enumerate() {
        initial.extend(
            input.impulse_guess_n_s[..coordinates[i].1]
                .iter()
                .map(|x| x * speed_scale),
        );
    }
    let construction = NewtonKrylovState::new(&problem, initial, cc.newton);
    if let Some(error) = problem.take_failure() {
        return Err(rigid_difference_error(error));
    }
    let mut solver = construction?;
    let result = solver.run_controlled(&problem, cc.maximum_newton_iterations, &mut || {
        (*cancel.borrow_mut())() || problem.has_failed()
    });
    if let Some(error) = problem.take_failure() {
        return Err(rigid_difference_error(error));
    }
    poll()?;
    if result.progress != SolverRunProgress::Complete {
        return Err(RigidCaloricContactError::Interrupted(result.progress));
    }
    if !result.report.converged {
        return Err(RigidCaloricContactError::NotConverged(result.report));
    }
    let trial = evaluate(&solver.x)?;
    let maximum_cell_residual_j = trial.residual_j()[..cells]
        .iter()
        .copied()
        .map(f64::abs)
        .fold(0.0, f64::max);
    let (next_thermal, staged, heat) = trial.accept(
        FrictionHeatTolerances {
            cell_j: cc.cell_energy_tolerance_j,
            thermal_j: cc.thermal_energy_tolerance_j,
            interface_j: controls.contact.interface_energy_tolerance_j,
            total_j: cc.total_energy_tolerance_j,
        },
        poll,
    )?;
    let mut expected_linear = vec![Vec3::ZERO; bodies.len()];
    let mut expected_angular = vec![Vec3::ZERO; bodies.len()];
    let mut work = 0.0;
    let mut internal_change = 0.0;
    let mut numerical = 0.0;
    for (i, input) in ordered.iter().enumerate() {
        poll()?;
        let m = &staged.mechanical[i];
        let impulse = vector(m.contact.wrench.force_n).scale(dt);
        let couple = vector(m.contact.wrench.torque_nm).scale(dt);
        // Aggregate body residuals alone can hide cancelling contact-law
        // errors in an island. Audit every original wrench before summing it.
        if norm(staged.impulses[i].sub(impulse)) > controls.linear_momentum_tolerance_n_s
            || norm(staged.couples[i].sub(couple)) > controls.angular_momentum_tolerance_nm_s
        {
            return Err(RigidCaloricContactError::Invalid(
                "an original contact impulse does not match its retained constitutive wrench",
            ));
        }
        for side in 0..2 {
            let body = input.bodies[side];
            let sign = if side == 0 { 1.0 } else { -1.0 };
            expected_linear[body] = expected_linear[body].add(impulse.scale(sign));
            expected_angular[body] = expected_angular[body].add(
                bodies[body]
                    .pose()
                    .orientation()
                    .rotate_body_to_world(input.arms_body_m[side])
                    .cross(impulse)
                    .add(couple)
                    .scale(sign),
            );
        }
        work = finite(work + m.contact.generalized_work.work_into_interface_j)?;
        internal_change = finite(internal_change + m.internal_change)?;
        numerical = finite(numerical + m.contact.generalized_work.numerical_dissipation_j)?;
    }
    let mut linear = Vec::with_capacity(bodies.len());
    let mut angular = Vec::with_capacity(bodies.len());
    let mut dp = Vec3::ZERO;
    let mut dl = Vec3::ZERO;
    let mut kinetic = 0.0;
    for (i, before) in bodies.iter().enumerate() {
        poll()?;
        let after = &staged.bodies[i];
        let actual = after
            .linear_momentum_world()
            .sub(before.linear_momentum_world());
        let actual_angular = before.pose().orientation().rotate_body_to_world(
            after
                .angular_momentum_body()
                .sub(before.angular_momentum_body()),
        );
        linear.push(actual.sub(expected_linear[i]));
        angular.push(actual_angular.sub(expected_angular[i]));
        dp = dp.add(actual);
        dl = dl
            .add(actual_angular)
            .add(before.pose().position_world().cross(actual));
        kinetic = finite(
            kinetic + after.kinetic_energy(properties[i])?
                - before.kinetic_energy(properties[i])?,
        )?;
    }
    let mechanical = finite(kinetic + work)?;
    let total = finite(
        kinetic + internal_change + heat.thermal.stored_energy_change_j + numerical
            - staged.external_heat_j,
    )?;
    if linear
        .iter()
        .chain([&dp])
        .any(|v| !v.is_finite() || norm(*v) > controls.linear_momentum_tolerance_n_s)
        || angular
            .iter()
            .chain([&dl])
            .any(|v| !v.is_finite() || norm(*v) > controls.angular_momentum_tolerance_nm_s)
        || mechanical.abs() > controls.mechanical_energy_tolerance_j
        || total.abs() > cc.total_energy_tolerance_j
    {
        return Err(RigidCaloricContactError::Invalid(
            "joint native momentum, work or original caloric energy failed balance",
        ));
    }
    let mut accepted_contacts = Vec::with_capacity(ordered.len());
    let mut reports = Vec::with_capacity(ordered.len());
    for (i, (mut m, h)) in staged
        .mechanical
        .into_iter()
        .zip(heat.interfaces)
        .enumerate()
    {
        poll()?;
        let fraction_a =
            finite(h.conductances_w_k[0] / finite(h.conductances_w_k[0] + h.conductances_w_k[1])?)?;
        let partition = HeatPartition::new(fraction_a, 1.0 - fraction_a, 0.0).map_err(|_| {
            RigidCaloricContactError::Invalid("unrepresentable conduction-derived heat share")
        })?;
        if let Some(sample) = &mut m.material {
            sample.partition = partition;
        }
        let physical_a = finite(m.contact.dissipation.heat_j * fraction_a)?;
        let physical_heat_j = [
            physical_a,
            finite(m.contact.dissipation.heat_j - physical_a)?,
        ];
        let deposited_heat_j = [
            finite(physical_heat_j[0] + m.reversible_heat_j[0])?,
            finite(physical_heat_j[1] + m.reversible_heat_j[1])?,
        ];
        let conduction_heat_j = [
            finite(h.body_heat_j[0] - physical_heat_j[0])?,
            finite(h.body_heat_j[1] - physical_heat_j[1])?,
        ];
        let mut accepted = ordered[i].contact.clone();
        accepted.checkpoint = m.contact.checkpoint.clone();
        accepted.temperature_k = m.temperature_k;
        accepted.thermodynamics = m.thermodynamics;
        accepted.time_s = target_s;
        if let (Some(binding), Some(sample)) = (&mut accepted.material_friction, &m.material) {
            binding.distance_m = sample.accumulated_distance_m;
            binding.time_s = target_s;
        }
        accepted_contacts.push(accepted);
        reports.push(RigidCaloricInterfaceReport {
            contact_impulse_residual_n_s: staged.impulses[i]
                .sub(vector(m.contact.wrench.force_n).scale(dt)),
            contact_spin_impulse_residual_nm_s: staged.couples[i]
                .sub(vector(m.contact.wrench.torque_nm).scale(dt)),
            heat: h,
            contact: m.contact,
            material: m.material,
            thermodynamics: m.thermodynamics,
            contact_internal_energy_change_j: m.internal_change,
            reversible_heat_j: m.reversible_heat_j,
            physical_heat_j,
            deposited_heat_j,
            conduction_heat_j,
            applied_impulse_on_a_n_s: staged.impulses[i],
            applied_spin_impulse_on_a_nm_s: staged.couples[i],
            midpoint_relative_velocity_mps: staged.slip[i],
            midpoint_relative_spin_rad_s: staged.spin[i],
        });
    }
    poll()?;
    Ok(RigidCaloricNetworkCandidate {
        bodies: staged.bodies,
        contacts: accepted_contacts,
        thermal: next_thermal,
        report: RigidCaloricNetworkReport {
            interfaces: reports,
            thermal: heat.thermal,
            external_heat_j: staged.external_heat_j,
            source: staged.source,
            nonlinear: result.report,
            trials: trials.get(),
            maximum_cell_residual_j,
            body_momentum_residual_n_s: linear,
            body_angular_momentum_residual_nm_s: angular,
            linear_momentum_change_n_s: dp,
            world_angular_momentum_change_nm_s: dl,
            kinetic_energy_change_j: kinetic,
            mechanical_energy_residual_j: mechanical,
            total_energy_residual_j: total,
        },
    })
}
