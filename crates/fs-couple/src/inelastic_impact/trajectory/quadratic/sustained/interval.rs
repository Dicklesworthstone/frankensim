//! One private original force/contact/drift trial; accepted owners are borrowed.

use super::{
    AcceptedInterval, Contact, Error, MixedViscoSustainedSpherePathReport,
    MixedViscoSustainedSphereTrajectoryControls, MixedViscoSustainedSphereTrajectoryStep,
    PreparedQuadraticSolidVolume, QuadraticSolidContactSurface, QuadraticSolidSurfaceError,
    QuadraticSolidVolumePathVerdict, QuadraticSphereSceneError, QuadraticSphereScenePathVerdict,
    Refusal, check_clock, poll, positions, vector,
};
use crate::finite_solid::contact::step::sustained::rigid::{rigid_contact_drift, rigid_force_kick};
use crate::inelastic_impact::{
    multiple::mixed_visco_mixed_contact_candidate, surface::scene::FiniteSolidSphereScene,
};
use fs_exec::Cx;
use fs_material::MaterialError;
use fs_mbd::{DynamicsError, Gravity, MassProperties, RigidBodyState, Wrench};
use fs_solid::hyper3::mixed::{
    MixedUpError,
    visco::{MixedViscoTetError, dynamics::MixedViscoTetMotion},
};

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub(super) fn advance<'a>(
    solid: &MixedViscoTetMotion<'a>,
    surface: &QuadraticSolidContactSurface,
    prepared_volume: &PreparedQuadraticSolidVolume<'_>,
    scene: &FiniteSolidSphereScene,
    start: &[[f64; 3]],
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    held: &[usize],
    contacts: &[Contact],
    nodal_forces: &[(usize, f64)],
    body_force: [f64; 3],
    wrenches: &[Wrench],
    gravity: Gravity,
    dt: f64,
    origin: [f64; 3],
    c: MixedViscoSustainedSphereTrajectoryControls,
    allowance: f64,
    clock_tolerance: f64,
    cx: &Cx<'_>,
) -> Result<Result<AcceptedInterval<'a>, Refusal>, Error> {
    let (kicked, kick) = solid
        .force_kick_candidate(
            dt,
            held,
            nodal_forces,
            body_force,
            origin,
            c.impact.solid_impulse,
            cx,
        )
        .map_err(Error::Solid)?;
    let mut kicks = Vec::with_capacity(bodies.len());
    let mut kicked_bodies = Vec::with_capacity(bodies.len());
    for ((&body, &property), &wrench) in bodies.iter().zip(properties).zip(wrenches) {
        poll(cx)?;
        let kick = rigid_force_kick(body, property, wrench, gravity, dt).map_err(Error::Rigid)?;
        kicked_bodies.push(kick.state_after);
        kicks.push(kick);
    }
    let (contacted, contacted_bodies, impact) = if contacts.is_empty() {
        (kicked, kicked_bodies, None)
    } else {
        let event = mixed_visco_mixed_contact_candidate(
            &kicked,
            held,
            &kicked_bodies,
            properties,
            contacts,
            origin,
            c.impact,
            cx,
        )
        .map_err(Error::Impact)?;
        (event.solid, event.rigid_bodies, Some(event.report))
    };
    let mut drift_controls = c.solid;
    drift_controls.maximum_integration_defect_j =
        drift_controls.maximum_integration_defect_j.min(allowance);
    let (next, solid_drift) =
        match contacted.contact_drift_candidate(&kick, c.settings, drift_controls, cx) {
            Ok(value) => value,
            Err(error) if refinable(&error) => return Ok(Err(Refusal::SolidTrial(error))),
            Err(error) => return Err(Error::Solid(error)),
        };
    check_clock(
        dt,
        next.solid().time_s() - solid.solid().time_s(),
        clock_tolerance,
    )
    .map_err(Error::Path)?;
    let end = positions(surface, &next, cx)?;
    let volume = prepared_volume
        .affine_path_at_positions(start, &end, c.volume, cx)
        .map_err(Error::Volume)?;
    if matches!(
        volume.verdict,
        QuadraticSolidVolumePathVerdict::Unresolved { .. }
    ) {
        return Ok(Err(Refusal::Volume(volume)));
    }
    let mut next_bodies = Vec::with_capacity(bodies.len());
    let mut rigid = Vec::with_capacity(bodies.len());
    let mut absolute_defect = solid_drift.absolute_integration_defect_j;
    for (index, (((&body, &property), &wrench), kick)) in contacted_bodies
        .iter()
        .zip(properties)
        .zip(wrenches)
        .zip(kicks)
        .enumerate()
    {
        poll(cx)?;
        let step =
            match rigid_contact_drift(kick, body, property, wrench, gravity, dt, vector(origin)) {
                Ok(value) => value,
                Err(DynamicsError::RotationalSolveFailed) => {
                    return Ok(Err(Refusal::RotationalSolve { rigid_body: index }));
                }
                Err(error) => return Err(Error::Rigid(error)),
            };
        let (linear, angular, energy) = step.phase_residuals();
        if [
            linear,
            angular,
            energy,
            step.kinetic_drift_defect_j,
            step.force_displacement_work_j,
            step.force_splitting_work_j,
            step.torque_impulse_work_j,
        ]
        .iter()
        .any(|v| !v.is_finite())
            || linear > c.impact.linear_tolerance_kg_m_s
            || angular > c.impact.angular_tolerance_kg_m2_s
            || energy > c.impact.energy_tolerance_j
        {
            return Err(Error::PhaseAudit {
                rigid_body: index,
                linear_residual_kg_m_s: linear,
                angular_residual_kg_m2_s: angular,
                energy_residual_j: energy,
            });
        }
        absolute_defect += step.force_splitting_work_j.abs() + step.kinetic_drift_defect_j.abs();
        next_bodies.push(step.drift.state_after);
        rigid.push(step);
    }
    if !absolute_defect.is_finite() {
        return Err(Error::Invalid(
            "unrepresentable absolute quadratic split defect",
        ));
    }
    if absolute_defect > allowance {
        return Ok(Err(Refusal::Integration {
            defect_j: absolute_defect,
            allowance_j: allowance,
        }));
    }
    let contact = scene
        .affine_quadratic_sphere_paths_at_positions(
            surface,
            start,
            &end,
            &contacted_bodies,
            &next_bodies,
            c.path,
            cx,
        )
        .map_err(Error::Scene)?;
    if matches!(
        contact.verdict,
        QuadraticSphereScenePathVerdict::PossibleContact { .. }
    ) {
        return Ok(Err(Refusal::Path(contact)));
    }
    let geometry = match scene.quadratic_contacts_at_positions(
        surface,
        &end,
        &next_bodies,
        c.path.initial,
        cx,
    ) {
        Ok(query) => query,
        Err(
            error @ QuadraticSphereSceneError::Surface(
                QuadraticSolidSurfaceError::UnresolvedGap { .. }
                | QuadraticSolidSurfaceError::UnresolvedActivation { .. },
            ),
        ) => return Ok(Err(Refusal::Endpoint(error))),
        Err(error) => return Err(Error::Scene(error)),
    };
    poll(cx)?;
    Ok(Ok(AcceptedInterval {
        solid: next,
        rigid: next_bodies,
        positions: end,
        geometry,
        step: MixedViscoSustainedSphereTrajectoryStep {
            elapsed_s: 0.0,
            duration_s: dt,
            maximum_contact_gap_m: 0.0,
            contact_keys: Vec::new(),
            incoming_normal_speeds_m_s: Vec::new(),
            normal_targets_m_s: Vec::new(),
            solid_kick: kick.into(),
            solid_drift,
            rigid,
            impact,
            path: MixedViscoSustainedSpherePathReport { volume, contact },
        },
    }))
}

fn refinable(error: &MixedViscoTetError) -> bool {
    matches!(
        error,
        MixedViscoTetError::EnergyGate { .. }
            | MixedViscoTetError::Native(MixedUpError::Material {
                error: MaterialError::State { .. },
                ..
            })
            | MixedViscoTetError::Native(MixedUpError::NewtonStalled { .. })
    )
}
