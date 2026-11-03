//! Native Maxwell motion, interval-admitted flight and simultaneous sphere impacts.
//!
//! Collision geometry follows the explicitly declared straight interpolation of
//! native endpoint coordinates. Original interval sweeps cover that entire path;
//! they do not enclose the unknown continuum flow. Every shortened retry starts
//! from the same accepted history and recomputes the original material interval.
//! Positions, velocities and Maxwell states are never repaired or interpolated.

/// Loaded sustained contact through native force, contact and material drift phases.
pub mod sustained;

use super::multiple::{
    InelasticMixedContactControls, InelasticMixedContactError, InelasticMixedContactReport,
    visco_mixed_contact_candidate,
};
use super::surface::{
    FiniteSolidContactSurface, FiniteSolidSphere, FiniteSolidSurfaceControls,
    FiniteSolidSurfaceError,
};
use crate::finite_solid::contact::step::{
    FiniteSolidSpherePathControls, FiniteSolidSpherePathReport, FiniteSolidSphereTrajectoryError,
    admit_path_at_positions,
    free::{FreeFlightError, RigidFreeFlight, rigid_free_flight},
};
use fs_exec::Cx;
use fs_mbd::{Gravity, MassProperties, RigidBodyState, StepReceipt, Vec3, Wrench};
use fs_solid::hyper3::{
    HyperTetSettings,
    visco::{
        FiniteViscoTetEnergy, FiniteViscoTetError,
        dynamics::{
            FiniteViscoMotionControls, FiniteViscoTetMotion, InertialEnergy, SolidMotionMomentum,
        },
    },
};

/// Native solver controls and whole-call allowances; none expands on a retry.
#[derive(Clone, Copy, Debug)]
pub struct ViscoSphereTrajectoryControls {
    pub maximum_step_s: f64,
    pub maximum_steps: usize,
    pub maximum_rejections: usize,
    pub path: FiniteSolidSpherePathControls,
    pub surface: FiniteSolidSurfaceControls,
    pub impact: InelasticMixedContactControls,
    pub settings: HyperTetSettings,
    pub solid: FiniteViscoMotionControls,
    /// Gates on accumulated residual magnitudes as well as final residuals.
    pub linear_tolerance_kg_m_s: f64,
    pub angular_tolerance_kg_m2_s: f64,
    pub energy_tolerance_j: f64,
    /// Sum of BE kinetic loss and absolute material/rigid numerical defects.
    pub integration_tolerance_j: f64,
    /// Sum of native BE angular-defect magnitudes. This is never physical torque.
    pub angular_integration_tolerance_kg_m2_s: f64,
}

/// Accepted native interval; rejected trials contribute neither history nor loss.
#[derive(Clone, Debug)]
pub struct ViscoSphereTrajectoryStep {
    /// Exact duration used by both native owners, retained for deterministic replay.
    pub duration_s: f64,
    pub elapsed_s: f64,
    pub maximum_contact_gap_m: f64,
    pub impact: Option<InelasticMixedContactReport>,
    pub material: FiniteViscoTetEnergy,
    pub inertia: InertialEnergy,
    pub momentum: SolidMotionMomentum,
    pub rigid: Vec<StepReceipt>,
    pub path: FiniteSolidSpherePathReport,
}

/// Physical losses and numerical errors remain separate throughout the call.
#[derive(Clone, Debug, Default)]
pub struct ViscoSphereTrajectoryReport {
    pub duration_s: f64,
    pub rejected_steps: usize,
    pub steps: Vec<ViscoSphereTrajectoryStep>,
    pub kinetic_change_j: f64,
    pub stored_change_j: f64,
    pub source_work_j: f64,
    pub support_work_j: f64,
    pub relaxation_loss_j: f64,
    pub unallocated_contact_loss_j: f64,
    pub numerical_kinetic_loss_j: f64,
    pub material_endpoint_defect_j: f64,
    pub rigid_drift_defect_j: f64,
    pub absolute_integration_defect_j: f64,
    pub work_energy_residual_j: f64,
    pub linear_momentum_change_kg_m_s: Vec3,
    pub angular_momentum_change_kg_m2_s: Vec3,
    pub external_impulse_kg_m_s: Vec3,
    pub external_angular_impulse_kg_m2_s: Vec3,
    pub support_impulse_kg_m_s: Vec3,
    pub support_angular_impulse_kg_m2_s: Vec3,
    pub angular_integration_defect_kg_m2_s: Vec3,
    pub absolute_angular_integration_defect_kg_m2_s: f64,
    pub linear_momentum_residual_kg_m_s: Vec3,
    pub angular_momentum_residual_kg_m2_s: Vec3,
    pub accumulated_linear_residual_kg_m_s: f64,
    pub accumulated_angular_residual_kg_m2_s: f64,
    pub accumulated_energy_residual_j: f64,
}

/// Adopt all native owners together after the complete requested interval passes.
#[derive(Debug)]
pub struct ViscoSphereTrajectoryCandidate<'a> {
    pub solid: FiniteViscoTetMotion<'a>,
    pub rigid_bodies: Vec<RigidBodyState>,
    pub report: ViscoSphereTrajectoryReport,
}

/// No error publishes a partially advanced solid, rigid body or material state.
#[derive(Debug)]
pub enum ViscoSphereTrajectoryError {
    Invalid(&'static str),
    Solid(FiniteViscoTetError),
    Surface(FiniteSolidSurfaceError),
    Impact(InelasticMixedContactError<FiniteViscoTetError>),
    Rigid(FreeFlightError),
    /// Original path, work-budget, unsupported sphere-pair or clock refusal.
    Path(FiniteSolidSphereTrajectoryError),
    Cancelled,
    Audit(Box<ViscoSphereTrajectoryReport>),
}
impl core::fmt::Display for ViscoSphereTrajectoryError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native viscoelastic sphere trajectory: {self:?}")
    }
}
impl std::error::Error for ViscoSphereTrajectoryError {}

/// Advance a real viscoelastic body and COM-centred spheres through transient impacts.
///
/// Held velocity DOFs must initially be stationary; their accepted displacement
/// is held for this call. Dead nodal/reference-volume loads remain fixed. Rigid
/// gravity belongs to the rigid owner; solid gravity must be supplied through
/// its original reference-acceleration load projection without double counting.
///
/// The positive activation band localizes contact without moving witnesses.
/// All active surface rows enter one original simultaneous impact solve. An
/// unresolved swept path shortens and recomputes the original native interval.
/// Each successful interval continues the same Maxwell branches and time.
/// A contact first entering the band at the final endpoint is resolved at the
/// beginning of the next continuation; that endpoint remains nonpenetrating.
///
/// This is an isothermal, first-order numerical trajectory, not true-flow CCD
/// or sustained-contact dynamics. Offset spheres, sphere-pair encounters without
/// a declared pair law, unresolved grazing/Zeno behavior and exhausted budgets
/// refuse. Physical relaxation/contact loss and numerical kinetic, angular and
/// endpoint-work defects never become inferred temperature changes.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn visco_sphere_trajectory_candidate<'a>(
    solid: &FiniteViscoTetMotion<'a>,
    surface: &FiniteSolidContactSurface,
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    spheres: &[FiniteSolidSphere],
    held_velocity_dofs: &[usize],
    nodal_forces_n: &[(usize, f64)],
    body_force_n_m3: [f64; 3],
    wrenches: &[Wrench],
    gravity: Gravity,
    duration_s: f64,
    origin_m: [f64; 3],
    controls: ViscoSphereTrajectoryControls,
    cx: &Cx<'_>,
) -> Result<ViscoSphereTrajectoryCandidate<'a>, ViscoSphereTrajectoryError> {
    use ViscoSphereTrajectoryError as E;
    poll(cx)?;
    validate(
        bodies, properties, spheres, wrenches, duration_s, origin_m, controls,
    )?;
    let ndofs = solid.velocity_mps().len() * 3;
    if held_velocity_dofs.len() > ndofs {
        return Err(E::Invalid("held velocities exceed the original solid DOFs"));
    }
    let mut held = held_velocity_dofs.to_vec();
    held.sort_unstable();
    if held.windows(2).any(|w| w[0] == w[1])
        || held
            .iter()
            .any(|&d| d >= ndofs || solid.velocity_mps()[d / 3][d % 3] != 0.0)
    {
        return Err(E::Invalid(
            "distinct valid stationary support DOFs required",
        ));
    }
    let prescribed: Vec<_> = held
        .iter()
        .map(|&d| (d, solid.solid().displacement_m()[d / 3][d % 3]))
        .collect();
    let mut spheres = spheres.to_vec();
    spheres.sort_unstable_by_key(|s| s.id);
    let start_time = solid.solid().time_s();
    let end_time = start_time + duration_s;
    if !end_time.is_finite() || end_time <= start_time {
        return Err(E::Path(FiniteSolidSphereTrajectoryError::NoProgress {
            elapsed_s: 0.0,
        }));
    }
    let mut current: Option<FiniteViscoTetMotion<'a>> = None;
    let mut rigid = bodies.to_vec();
    let mut report = ViscoSphereTrajectoryReport {
        duration_s,
        ..Default::default()
    };
    let mut elapsed = 0.0;
    // Subtracting an existing material age loses precision at that age's ULP,
    // even when the requested short interval itself is represented accurately.
    // Both native owners still receive the exact same retained per-step dt.
    let clock_tolerance = 8.0 * f64::EPSILON * start_time.abs().max(end_time.abs()).max(duration_s);
    let mut geometry = surface
        .visco_sphere_separations(solid, &rigid, &spheres, controls.surface, cx)
        .map_err(E::Surface)?;
    let mut contacts = geometry.admitted_inelastic_contacts().map_err(E::Surface)?;
    while elapsed < duration_s {
        poll(cx)?;
        if report.steps.len() == controls.maximum_steps {
            return Err(E::Path(FiniteSolidSphereTrajectoryError::StepBudget {
                elapsed_s: elapsed,
            }));
        }
        let accepted = current.as_ref().unwrap_or(solid);
        let maximum_gap = geometry
            .contacts
            .iter()
            .map(|w| w.gap_m)
            .fold(0.0, f64::max);
        let event = if contacts.is_empty() {
            None
        } else {
            Some(
                visco_mixed_contact_candidate(
                    accepted,
                    &held,
                    &rigid,
                    properties,
                    &contacts,
                    origin_m,
                    controls.impact,
                    cx,
                )
                .map_err(E::Impact)?,
            )
        };
        let post = event.as_ref().map_or(accepted, |e| &e.solid);
        let post_rigid = event
            .as_ref()
            .map_or(rigid.as_slice(), |e| e.rigid_bodies.as_slice());
        let start = positions(surface, post, cx)?;
        let mut dt = controls.maximum_step_s.min(duration_s - elapsed);
        let (step, momentum, flights, rigid_next, path, geometry_next, contacts_next) = loop {
            poll(cx)?;
            if !dt.is_finite()
                || dt <= 0.0
                || elapsed + dt <= elapsed
                || post.solid().time_s() + dt <= post.solid().time_s()
            {
                return Err(E::Path(FiniteSolidSphereTrajectoryError::NoProgress {
                    elapsed_s: elapsed,
                }));
            }
            let (step, momentum) = post
                .advance_with_momentum(
                    dt,
                    &prescribed,
                    nodal_forces_n,
                    body_force_n_m3,
                    controls.settings,
                    controls.solid,
                    origin_m,
                    cx,
                )
                .map_err(E::Solid)?;
            let represented = step.motion.solid().time_s() - post.solid().time_s();
            check_clock(dt, represented, clock_tolerance)?;
            let flights: Vec<_> = post_rigid
                .iter()
                .zip(properties)
                .zip(wrenches)
                .map(|((&old, &property), &wrench)| {
                    rigid_free_flight(old, property, wrench, gravity, dt, vector(origin_m), cx)
                })
                .collect::<Result<_, _>>()
                .map_err(E::Rigid)?;
            let next: Vec<_> = flights.iter().map(|f| f.receipt.state_after).collect();
            let end = positions(surface, &step.motion, cx)?;
            match admit_path_at_positions(
                surface,
                post.solid().reference_tetrahedra(),
                &start,
                &end,
                post_rigid,
                &next,
                &spheres,
                false,
                controls.path,
                cx,
            )
            .map_err(E::Path)?
            {
                Ok(path) => {
                    let geometry = surface
                        .visco_sphere_separations(
                            &step.motion,
                            &next,
                            &spheres,
                            controls.surface,
                            cx,
                        )
                        .map_err(E::Surface)?;
                    let contacts = geometry.admitted_inelastic_contacts().map_err(E::Surface)?;
                    break (step, momentum, flights, next, path, geometry, contacts);
                }
                Err(refusal) => {
                    if report.rejected_steps == controls.maximum_rejections {
                        return Err(E::Path(FiniteSolidSphereTrajectoryError::RejectionBudget {
                            elapsed_s: elapsed,
                            last: refusal,
                        }));
                    }
                    report.rejected_steps += 1;
                    dt *= if refusal.start() > 0.0 {
                        0.9 * refusal.start()
                    } else {
                        0.5
                    };
                }
            }
        };
        let impact = event.map(|e| e.report);
        accumulate(
            &mut report,
            &step.material_energy,
            step.inertia,
            momentum,
            &flights,
            impact.as_ref(),
        );
        elapsed += dt;
        if elapsed > duration_s {
            return Err(E::Path(FiniteSolidSphereTrajectoryError::NoProgress {
                elapsed_s: elapsed,
            }));
        }
        check_clock(
            elapsed,
            step.motion.solid().time_s() - start_time,
            clock_tolerance,
        )?;
        report.steps.push(ViscoSphereTrajectoryStep {
            duration_s: dt,
            elapsed_s: elapsed,
            maximum_contact_gap_m: maximum_gap,
            impact,
            material: step.material_energy,
            inertia: step.inertia,
            momentum,
            rigid: flights.into_iter().map(|f| f.receipt).collect(),
            path,
        });
        current = Some(step.motion);
        rigid = rigid_next;
        // The endpoint query already describes the next accepted configuration;
        // its zero-duration impulse changes no positions. Reuse the BVH result.
        geometry = geometry_next;
        contacts = contacts_next;
        finish(&mut report);
        if !admitted(&report, controls) {
            return Err(E::Audit(Box::new(report)));
        }
    }
    let solid = current.expect("positive interval accepted a native material step");
    check_clock(
        duration_s,
        solid.solid().time_s() - start_time,
        clock_tolerance,
    )?;
    poll(cx)?;
    Ok(ViscoSphereTrajectoryCandidate {
        solid,
        rigid_bodies: rigid,
        report,
    })
}

fn accumulate(
    r: &mut ViscoSphereTrajectoryReport,
    material: &FiniteViscoTetEnergy,
    inertia: InertialEnergy,
    momentum: SolidMotionMomentum,
    flights: &[RigidFreeFlight],
    impact: Option<&InelasticMixedContactReport>,
) {
    r.kinetic_change_j += inertia.kinetic_change_j;
    r.stored_change_j += material.stored_after_j - material.stored_before_j;
    r.source_work_j += material.dead_load_endpoint_work_j;
    r.support_work_j += material.support_endpoint_work_j;
    r.relaxation_loss_j += material.dissipated_j;
    r.numerical_kinetic_loss_j += inertia.numerical_loss_j;
    r.material_endpoint_defect_j += material.endpoint_work_defect_j;
    r.absolute_integration_defect_j +=
        inertia.numerical_loss_j + material.endpoint_work_defect_j.abs();
    r.linear_momentum_change_kg_m_s = r
        .linear_momentum_change_kg_m_s
        .add(vector(momentum.change.linear_kg_m_s));
    r.angular_momentum_change_kg_m2_s = r
        .angular_momentum_change_kg_m2_s
        .add(vector(momentum.change.angular_kg_m2_s));
    r.external_impulse_kg_m_s = r
        .external_impulse_kg_m_s
        .add(vector(momentum.applied.linear_kg_m_s));
    r.external_angular_impulse_kg_m2_s = r
        .external_angular_impulse_kg_m2_s
        .add(vector(momentum.applied.angular_kg_m2_s));
    r.support_impulse_kg_m_s = r
        .support_impulse_kg_m_s
        .add(vector(momentum.support.linear_kg_m_s));
    r.support_angular_impulse_kg_m2_s = r
        .support_angular_impulse_kg_m2_s
        .add(vector(momentum.support.angular_kg_m2_s));
    r.angular_integration_defect_kg_m2_s = r
        .angular_integration_defect_kg_m2_s
        .add(vector(momentum.angular_integration_defect_kg_m2_s));
    r.absolute_angular_integration_defect_kg_m2_s +=
        norm(vector(momentum.angular_integration_defect_kg_m2_s));
    r.accumulated_linear_residual_kg_m_s += norm(vector(momentum.linear_residual_kg_m_s));
    r.accumulated_angular_residual_kg_m2_s += norm(vector(momentum.angular_residual_kg_m2_s));
    r.accumulated_energy_residual_j +=
        material.equilibrium_work_residual_j.abs() + inertia.balance_error_j.abs();
    for f in flights {
        r.kinetic_change_j += f.kinetic_change_j;
        r.source_work_j += f.source_work_j;
        r.rigid_drift_defect_j += f.drift_defect_j;
        r.absolute_integration_defect_j += f.drift_defect_j.abs();
        r.linear_momentum_change_kg_m_s = r.linear_momentum_change_kg_m_s.add(f.linear_change);
        r.angular_momentum_change_kg_m2_s = r.angular_momentum_change_kg_m2_s.add(f.angular_change);
        r.external_impulse_kg_m_s = r.external_impulse_kg_m_s.add(f.applied);
        r.external_angular_impulse_kg_m2_s =
            r.external_angular_impulse_kg_m2_s.add(f.applied_angular);
        r.accumulated_linear_residual_kg_m_s += norm(f.linear_residual);
        r.accumulated_angular_residual_kg_m2_s += norm(f.angular_residual);
        r.accumulated_energy_residual_j += f.work_residual_j.abs();
    }
    if let Some(i) = impact {
        let support = vector(i.solid_impulse.support.linear_kg_m_s);
        let moment = vector(i.solid_impulse.support.angular_kg_m2_s);
        r.kinetic_change_j += i.kinetic_change_j;
        r.support_work_j += i.support_work_j;
        r.unallocated_contact_loss_j += i.contact.dissipated_energy;
        r.linear_momentum_change_kg_m_s = r
            .linear_momentum_change_kg_m_s
            .add(support)
            .add(i.linear_momentum_residual_kg_m_s);
        r.angular_momentum_change_kg_m2_s = r
            .angular_momentum_change_kg_m2_s
            .add(moment)
            .add(i.angular_momentum_residual_kg_m2_s);
        r.support_impulse_kg_m_s = r.support_impulse_kg_m_s.add(support);
        r.support_angular_impulse_kg_m2_s = r.support_angular_impulse_kg_m2_s.add(moment);
        r.accumulated_linear_residual_kg_m_s +=
            i.absolute_linear_residual_kg_m_s + norm(i.linear_momentum_residual_kg_m_s);
        r.accumulated_angular_residual_kg_m2_s +=
            i.absolute_angular_residual_kg_m2_s + norm(i.angular_momentum_residual_kg_m2_s);
        r.accumulated_energy_residual_j += i.absolute_work_residual_j
            + i.work_energy_residual_j.abs()
            + i.dissipation_residual_j.abs()
            + (i.relative_impulse_work_j - i.contact.relative_work).abs();
    }
}

fn finish(r: &mut ViscoSphereTrajectoryReport) {
    r.work_energy_residual_j = r.kinetic_change_j
        + r.stored_change_j
        + r.relaxation_loss_j
        + r.unallocated_contact_loss_j
        + r.numerical_kinetic_loss_j
        + r.material_endpoint_defect_j
        - r.source_work_j
        - r.support_work_j
        - r.rigid_drift_defect_j;
    r.linear_momentum_residual_kg_m_s = r
        .linear_momentum_change_kg_m_s
        .sub(r.external_impulse_kg_m_s)
        .sub(r.support_impulse_kg_m_s);
    r.angular_momentum_residual_kg_m2_s = r
        .angular_momentum_change_kg_m2_s
        .sub(r.external_angular_impulse_kg_m2_s)
        .sub(r.support_angular_impulse_kg_m2_s)
        .sub(r.angular_integration_defect_kg_m2_s);
}

fn admitted(r: &ViscoSphereTrajectoryReport, c: ViscoSphereTrajectoryControls) -> bool {
    [
        r.kinetic_change_j,
        r.stored_change_j,
        r.source_work_j,
        r.support_work_j,
        r.relaxation_loss_j,
        r.unallocated_contact_loss_j,
        r.numerical_kinetic_loss_j,
        r.material_endpoint_defect_j,
        r.rigid_drift_defect_j,
        r.work_energy_residual_j,
        r.absolute_integration_defect_j,
        r.absolute_angular_integration_defect_kg_m2_s,
        r.accumulated_linear_residual_kg_m_s,
        r.accumulated_angular_residual_kg_m2_s,
        r.accumulated_energy_residual_j,
    ]
    .iter()
    .all(|v| v.is_finite())
        && [
            r.linear_momentum_change_kg_m_s,
            r.angular_momentum_change_kg_m2_s,
            r.external_impulse_kg_m_s,
            r.external_angular_impulse_kg_m2_s,
            r.support_impulse_kg_m_s,
            r.support_angular_impulse_kg_m2_s,
            r.angular_integration_defect_kg_m2_s,
            r.linear_momentum_residual_kg_m_s,
            r.angular_momentum_residual_kg_m2_s,
        ]
        .iter()
        .all(|v| v.is_finite())
        && r.accumulated_linear_residual_kg_m_s <= c.linear_tolerance_kg_m_s
        && r.accumulated_angular_residual_kg_m2_s <= c.angular_tolerance_kg_m2_s
        && r.accumulated_energy_residual_j <= c.energy_tolerance_j
        && norm(r.linear_momentum_residual_kg_m_s) <= c.linear_tolerance_kg_m_s
        && norm(r.angular_momentum_residual_kg_m2_s) <= c.angular_tolerance_kg_m2_s
        && r.work_energy_residual_j.abs() <= c.energy_tolerance_j
        && r.absolute_integration_defect_j <= c.integration_tolerance_j
        && r.absolute_angular_integration_defect_kg_m2_s <= c.angular_integration_tolerance_kg_m2_s
}

#[allow(clippy::too_many_arguments)]
fn validate(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    spheres: &[FiniteSolidSphere],
    wrenches: &[Wrench],
    duration: f64,
    origin: [f64; 3],
    c: ViscoSphereTrajectoryControls,
) -> Result<(), ViscoSphereTrajectoryError> {
    if bodies.len() != properties.len()
        || bodies.len() != wrenches.len()
        || bodies.len() != spheres.len()
        || bodies.len() > c.surface.maximum_spheres
        || c.impact.maximum_rigid_bodies == 0
        || bodies.len() > c.impact.maximum_rigid_bodies
        || !(1..=21).contains(&c.impact.contact.maximum_contacts)
        || !(1..=65_536).contains(&c.maximum_steps)
        || c.maximum_rejections > 65_536
        || !c.path.valid()
        || origin.iter().any(|v| !v.is_finite())
        || !c.surface.activation_tolerance_m.is_finite()
        || c.surface.activation_tolerance_m <= 0.0
        || !c.impact.geometry_tolerance_m.is_finite()
        || c.surface.activation_tolerance_m > c.impact.geometry_tolerance_m
        || [
            duration,
            c.maximum_step_s,
            c.linear_tolerance_kg_m_s,
            c.angular_tolerance_kg_m2_s,
            c.energy_tolerance_j,
            c.integration_tolerance_j,
            c.angular_integration_tolerance_kg_m2_s,
        ]
        .iter()
        .any(|v| !v.is_finite() || *v <= 0.0)
    {
        return Err(ViscoSphereTrajectoryError::Invalid(
            "bounded matching sphere scene and positive physical/time gates required",
        ));
    }
    let mut occupied = vec![false; bodies.len()];
    for sphere in spheres {
        if sphere.rigid_body >= bodies.len()
            || occupied[sphere.rigid_body]
            || sphere.center_body_m != Vec3::ZERO
        {
            return Err(ViscoSphereTrajectoryError::Invalid(
                "exactly one COM-centred sphere per rigid body required",
            ));
        }
        occupied[sphere.rigid_body] = true;
    }
    Ok(())
}

fn positions(
    surface: &FiniteSolidContactSurface,
    solid: &FiniteViscoTetMotion<'_>,
    cx: &Cx<'_>,
) -> Result<Vec<[f64; 3]>, ViscoSphereTrajectoryError> {
    let state = solid.solid();
    surface
        .current_positions_from_reference(
            state.reference_nodes_m(),
            state.reference_tetrahedra(),
            state.displacement_m(),
            cx,
        )
        .map_err(ViscoSphereTrajectoryError::Surface)
}
fn check_clock(
    requested: f64,
    represented: f64,
    tolerance: f64,
) -> Result<(), ViscoSphereTrajectoryError> {
    if (represented - requested).abs() > tolerance {
        Err(ViscoSphereTrajectoryError::Path(
            FiniteSolidSphereTrajectoryError::TimeResolution {
                requested_step_s: requested,
                represented_step_s: represented,
            },
        ))
    } else {
        Ok(())
    }
}
fn vector(v: [f64; 3]) -> Vec3 {
    Vec3::new(v[0], v[1], v[2])
}
fn norm(v: Vec3) -> f64 {
    v.x.hypot(v.y).hypot(v.z)
}
fn poll(cx: &Cx<'_>) -> Result<(), ViscoSphereTrajectoryError> {
    cx.checkpoint()
        .map_err(|_| ViscoSphereTrajectoryError::Cancelled)
}
