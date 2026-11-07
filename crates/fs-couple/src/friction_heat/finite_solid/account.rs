//! Original feature/patch admission and the native tangential work account.

use super::*;
use crate::finite_solid::contact::mixed::FiniteSolidContactCounterpart;
use fs_contact::interface_binding::thermal_friction::FrictionQueryAxisKind;
use fs_mbd::Vec3;
use fs_qty::semantic::QuantityKind;
use fs_qty::{Pressure, QuantitySpec, Temperature};
use fs_solid::hyper3::HyperelasticLawField;

fn invalid(s: &'static str) -> DepositionError {
    DepositionError::InvalidInput(s)
}

pub(super) fn source(c: &FiniteSolidHeatContact<'_>) -> Result<(), DepositionError> {
    let law = &c.binding.source.law.principal_sources()[0];
    let point = law.initial_source().query();
    if law.axes().nominal_pressure.is_some()
        || law
            .query_axes()
            .iter()
            .any(|a| a.kind == FrictionQueryAxisKind::NominalPressure)
        || point
            .axis_quantities()
            .values()
            .any(|q| q.dims() == Pressure::DIMS)
    {
        return Err(invalid(
            "finite-solid heat has no sustained source pressure closure",
        ));
    }
    for name in point.axes().keys() {
        let quantity = point.axis_quantities().get(name).copied();
        let temperature = quantity == Some(QuantitySpec::dimensional(Temperature::DIMS))
            || quantity
                .and_then(QuantitySpec::semantic_type)
                .is_some_and(|q| q.kind() == QuantityKind::AbsoluteTemperature)
            || (name == "T" && quantity.is_none());
        if temperature
            && !c.material.point.is_thermal_axis(name, Temperature::DIMS)
            && c.material.interface_temperature_axis != Some(name.as_str())
        {
            return Err(invalid(
                "every source temperature requires its original thermal owner",
            ));
        }
    }
    Ok(())
}

fn vector(v: [f64; 3]) -> Vec3 {
    Vec3::new(v[0], v[1], v[2])
}
fn norm(v: Vec3) -> f64 {
    fs_math::det::hypot(fs_math::det::hypot(v.x, v.y), v.z)
}

pub(super) fn patch(
    c: &FiniteSolidHeatContact<'_>,
    normal: Vec3,
    witnesses: [Vec3; 2],
    controls: FiniteSolidHeatControls,
) -> Result<(), DepositionError> {
    let geometry = c.route.geometry;
    let centre = [geometry[0].start_m, geometry[0].end_m]
        .into_iter()
        .find(|p| *p == geometry[1].start_m || *p == geometry[1].end_m)
        .ok_or(invalid("original thermal face centre missing"))?;
    let other = if geometry[0].start_m == centre {
        geometry[0].end_m
    } else {
        geometry[0].start_m
    };
    let centre = vector(centre);
    let a = vector(other).sub(centre);
    let face_normal = a.scale(1.0 / norm(a));
    if !face_normal.is_finite() || norm(normal.sub(face_normal)) > controls.normal_tolerance {
        return Err(invalid(
            "native contact changed its ordered thermal face normal",
        ));
    }
    for point in witnesses {
        let delta = point.sub(centre);
        let depth = delta.dot(face_normal);
        let radius = norm(delta.sub(face_normal.scale(depth)));
        if !depth.is_finite()
            || !radius.is_finite()
            || depth.abs() > controls.geometry_tolerance_m
            || radius > c.patch_radius_m + controls.geometry_tolerance_m
        {
            return Err(invalid(
                "native witness left its original finite thermal patch",
            ));
        }
    }
    Ok(())
}

pub(super) fn scene<E>(
    motion: &FiniteSolidHeatMotion<'_, '_>,
    solid: &FiniteSolidBody<'_>,
    bodies: &[RigidBodyState],
    contacts: &[&FiniteSolidHeatContact<'_>],
    controls: FiniteSolidHeatControls,
    cx: &Cx<'_>,
) -> Result<(), FiniteSolidHeatError<E>> {
    let query = motion
        .scene
        .contacts(
            motion.surface,
            solid,
            bodies,
            controls.mechanical.surface,
            cx,
        )
        .map_err(FiniteSolidHeatError::Geometry)?;
    if query.contacts.len() != contacts.len() {
        return Err(FiniteSolidHeatError::Invalid(
            "complete fixed native thermal contact set required",
        ));
    }
    let rigid_point = |p: fs_mbd::frictional_impact::RigidImpactPoint| {
        let body = bodies.get(p.body).ok_or(FiniteSolidHeatError::Invalid(
            "thermal witness body missing",
        ))?;
        let property = motion
            .properties
            .get(p.body)
            .ok_or(FiniteSolidHeatError::Invalid(
                "thermal witness inertia missing",
            ))?;
        body.point_kinematics(*property, p.arm_body)
            .map(|k| k.point_world)
            .map_err(FiniteSolidHeatError::Rigid)
    };
    for (native, key) in query.contacts.iter().zip(&query.keys) {
        cx.checkpoint().map_err(|_| DepositionError::Cancelled)?;
        let c = contacts
            .iter()
            .find(|c| c.contact_id == native.id && c.key == *key)
            .ok_or(FiniteSolidHeatError::Invalid(
                "native scene feature has no original thermal route",
            ))?;
        let first = rigid_point(native.point)?;
        let second = match native.counterpart {
            FiniteSolidContactCounterpart::Rigid(p) => rigid_point(p)?,
            FiniteSolidContactCounterpart::Solid(p) => {
                let material = match solid.problem().material {
                    HyperelasticLawField::Uniform(m) => m,
                    HyperelasticLawField::PerElement(m) => m.get(p.tetrahedron()).ok_or(
                        FiniteSolidHeatError::Invalid("thermal solid witness material missing"),
                    )?,
                };
                if material.density().card_identity()
                    != c.material.point.materials[1].content_hash()
                {
                    return Err(FiniteSolidHeatError::Invalid(
                        "thermal counterpart is not the original solid material",
                    ));
                }
                vector(
                    solid
                        .point_position_m(p, cx)
                        .map_err(FiniteSolidHeatError::Solid)?,
                )
            }
        };
        patch(c, native.normal_world, [first, second], controls)?;
    }
    Ok(())
}

pub(super) fn work(
    report: &FiniteSolidSustainedSourceTrajectoryReport,
    contacts: &[&FiniteSolidHeatContact<'_>],
    cx: &Cx<'_>,
) -> Result<(Vec<f64>, f64, f64), DepositionError> {
    let mut heat = vec![0.0; contacts.len()];
    let (mut eligible, mut exported) = (0.0, 0.0);
    for step in &report.steps {
        cx.checkpoint().map_err(|_| DepositionError::Cancelled)?;
        let event = step
            .event_data
            .as_ref()
            .ok_or(invalid("fixed thermal footprint lost native contact"))?;
        if event.work.len() != contacts.len() {
            return Err(invalid("incomplete native contact work"));
        }
        let mut seen = vec![false; contacts.len()];
        for work in &event.work {
            let i = contacts
                .iter()
                .position(|c| c.contact_id == work.contact_id && c.key == work.key)
                .ok_or(invalid(
                    "unmapped native work cannot be assigned to a heat cell",
                ))?;
            if seen[i]
                || work.authority != FiniteSolidFrictionHeatAuthority::UniqueSustainedImpulse
                || !event.motion.iter().any(|m| {
                    m.contact_id == work.contact_id && m.key == work.key && m.motion.is_some()
                })
            {
                return Err(invalid(
                    "only unique sustained passive tangential work may heat the fixed patch",
                ));
            }
            seen[i] = true;
            let loss = work
                .friction_heat_candidate_j
                .ok_or(invalid("native tangential heat authority missing"))?;
            let deposited = contacts[i].heat_fraction * loss;
            if !loss.is_finite()
                || loss < 0.0
                || (loss > 0.0 && contacts[i].heat_fraction > 0.0 && deposited == 0.0)
            {
                return Err(DepositionError::Arithmetic);
            }
            heat[i] += deposited;
            exported += loss - deposited;
            eligible += loss;
        }
    }
    let unconverted = report.unallocated_contact_loss_j - eligible;
    if report.steps.is_empty()
        || heat
            .iter()
            .chain([&exported, &eligible, &unconverted])
            .any(|x| !x.is_finite())
    {
        return Err(DepositionError::Arithmetic);
    }
    Ok((heat, exported, unconverted))
}
