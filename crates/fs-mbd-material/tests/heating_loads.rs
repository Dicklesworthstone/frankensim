//! G1/G3: native loads retain material geometry and the continuing heat history.
#[path = "heating/common.rs"]
mod common;
use common::*;
use fs_ga::{Se3, Twist, Vec3};
use fs_matdb::{ClaimSet, MaterialCard, PropertyValue};
use fs_material::phase::{ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_mbd::{Gravity, RigidBodyIntegrator};
use fs_mbd_material::ResolvedBody;
use fs_mbd_material::principal::transfer::ConfigurationTransferControls;
use fs_mbd_material::thermal::heating::motion::loads::{
    MaterialAttachment, MaterialImpulse, MaterialImpulseControls,
};
use fs_mbd_material::thermal::heating::motion::{HeatedMotionControls, HeatedRigidBody};
use fs_thermal::enthalpy::CellId;

fn motion<'a>(
    reference: &'a ResolvedBody,
    material: &'a MaterialCard,
    chart: &'a ResolvedHeatCapacityCurve,
    watts: f64,
) -> HeatedRigidBody<'a> {
    HeatedRigidBody::new(
        heated(reference, material, chart, watts, vec![]),
        Se3::identity(),
        Twist::new(Vec3::default(), Vec3::default()),
        RigidBodyIntegrator::new(Gravity::ZERO),
        HeatedMotionControls {
            maximum_step_s: 0.0625,
            relative_tensor_tolerance: 1e-12,
            transfer: ConfigurationTransferControls {
                relative_momentum_tolerance: 1e-12,
                maximum_parameter_work_j: 10.,
            },
            maximum_numerical_energy_change_j: 1e-7,
        },
        || Ok(()),
    )
    .unwrap()
}
fn controls() -> MaterialImpulseControls {
    MaterialImpulseControls {
        maximum_impulses: 4,
        work_residual_tolerance_j: 1e-10,
    }
}
fn attachment(y: f64) -> MaterialAttachment {
    MaterialAttachment {
        cell: CellId(0),
        reference_point_m: Vec3::new(0., y, 0.),
    }
}
fn load(y: f64, impulse: f64) -> MaterialImpulse {
    MaterialImpulse {
        attachment: attachment(y),
        impulse_world: fs_mbd::Vec3::new(impulse, 0., 0.),
    }
}
fn same(a: &HeatedRigidBody<'_>, b: &HeatedRigidBody<'_>) {
    assert_eq!(a.state(), b.state());
    assert_eq!(a.configuration(), b.configuration());
    assert_eq!(a.thermal().cells(), b.thermal().cells());
    assert_eq!(a.thermal().time_s(), b.thermal().time_s());
}

#[test]
fn off_centre_external_impulses_keep_work_heat_and_continued_motion_distinct() {
    let material = card(10., 0., 500.);
    let reference = reference(&material, 1);
    let chart = caloric(&material, 1e5, SolidLiquidPhase::Solid);
    let original = motion(&reference, &material, &chart, 20.);
    let saved = original.clone();
    let (loaded, audit) = original
        .apply_material_impulses(&[load(0.25, 2.), load(0.25, -1.)], controls(), || false)
        .unwrap();
    // Unit cube, rho=2: mass=2 kg, Izz=1/3 kg m². The second event extracts work.
    close(loaded.state().linear_momentum_world().x, 1.);
    close(
        loaded.diagnostics().unwrap().angular_momentum_world.z,
        -0.25,
    );
    close(audit.external_work_j, 1. / 4. + 3. / 32.);
    close(audit.kinetic_energy_change_j, audit.external_work_j);
    close(audit.work_energy_residual_j, 0.);
    close(audit.receipts()[0].work.impulse_work_j, 1.375);
    close(audit.receipts()[1].work.impulse_work_j, -1.03125);
    assert_eq!(
        audit.receipts()[0].state_after,
        audit.receipts()[1].state_before
    );
    assert_eq!(loaded.state().pose(), original.state().pose());
    assert_eq!(loaded.thermal().cells(), original.thermal().cells());
    close(audit.time_s, 0.);
    let point = loaded.material_point_kinematics(attachment(0.25)).unwrap();
    close(point.point_world.y, 0.25);
    close(point.point_velocity_world.x, 0.6875);
    let (later, heating) = loaded.advance_to(0.5, config(), || Ok(())).unwrap();
    close(later.state().pose().position_world().x, 0.25);
    close(later.thermal().temperature_k(CellId(0)).unwrap(), 300.5);
    close(later.thermal().energy_j(), 10.);
    close(heating.thermal.external_heat_j, 10.);
    close(heating.combined_residual_j, 0.);
    same(&original, &saved);
}

#[test]
fn source_expansion_moves_the_load_arm_and_source_density_changes_the_response() {
    let material = card(10., 1e-2, 500.);
    let reference = reference(&material, 1);
    let chart = caloric(&material, 1e5, SolidLiquidPhase::Solid);
    let cold = motion(&reference, &material, &chart, 200.);
    let (hot, _) = cold.advance_to(1., config(), || Ok(())).unwrap();
    let (loaded, audit) = hot
        .apply_material_impulses(&[load(0.25, 1.)], controls(), || false)
        .unwrap();
    let stretch: f64 = 1.1;
    let angular_impulse = -0.25 * stretch;
    close(
        audit.receipts()[0].point_kinematics_before.point_world.y,
        0.25 * stretch,
    );
    close(
        loaded.diagnostics().unwrap().angular_momentum_world.z,
        angular_impulse,
    );
    close(
        loaded
            .material_point_kinematics(attachment(0.25))
            .unwrap()
            .angular_velocity_world
            .z,
        angular_impulse / (stretch * stretch / 3.),
    );
    close(
        audit.external_work_j,
        0.25 + angular_impulse.powi(2) / (2. * stretch * stretch / 3.),
    );
    assert_eq!(loaded.thermal().cells(), hot.thermal().cells());
    close(loaded.thermal().time_s(), 1.);
    let mut velocities = Vec::new();
    for density in [2., 8.] {
        let source = card(10., 0., 500.);
        let mut claims = ClaimSet::new();
        for (_, original) in source.claims().claims_ordered() {
            let mut claim = original.clone();
            if claim.key.name() == "density" {
                let PropertyValue::Curve { knots, .. } = &mut claim.value else {
                    panic!("fixture curve")
                };
                for (_, value) in knots {
                    *value = density;
                }
            }
            claims.insert_claim(claim).unwrap();
        }
        let source = MaterialCard::assemble(source.id().clone(), claims, vec![]).unwrap();
        let reference = self::reference(&source, 1);
        let chart = caloric(&source, 1e5, SolidLiquidPhase::Solid);
        let body = motion(&reference, &source, &chart, 0.);
        let (loaded, audit) = body
            .apply_material_impulses(&[load(0., 1.)], controls(), || false)
            .unwrap();
        let velocity = loaded
            .material_point_kinematics(attachment(0.))
            .unwrap()
            .point_velocity_world
            .x;
        close(velocity, 1. / density);
        close(audit.external_work_j, 0.5 / density);
        velocities.push(velocity);
    }
    close(velocities[0] / velocities[1], 4.);
}

#[test]
fn a_late_bad_load_or_cancellation_cannot_publish_a_partial_event_batch() {
    let material = card(10., 0., 500.);
    let reference = reference(&material, 1);
    let chart = caloric(&material, 1e5, SolidLiquidPhase::Solid);
    let original = motion(&reference, &material, &chart, 20.);
    let saved = original.clone();
    let valid = load(0.25, 1.);
    let mut foreign = valid;
    foreign.attachment.cell = CellId(99);
    let mut invalid_point = valid;
    invalid_point.attachment.reference_point_m.x = f64::NAN;
    let mut invalid_impulse = valid;
    invalid_impulse.impulse_world.x = f64::INFINITY;
    for invalid in [foreign, invalid_point, invalid_impulse] {
        assert!(
            original
                .apply_material_impulses(&[valid, invalid], controls(), || false)
                .is_err()
        );
        same(&original, &saved);
    }
    let mut calls = 0;
    original
        .apply_material_impulses(&[valid, valid], controls(), || {
            calls += 1;
            false
        })
        .unwrap();
    let mut now = 0;
    assert!(
        original
            .apply_material_impulses(&[valid, valid], controls(), || {
                now += 1;
                now == calls
            })
            .is_err()
    );
    assert_eq!(now, calls);
    same(&original, &saved);
    let mut limited = controls();
    limited.maximum_impulses = 1;
    assert!(
        original
            .apply_material_impulses(&[valid, valid], limited, || false)
            .is_err()
    );
    limited = controls();
    limited.work_residual_tolerance_j = f64::NAN;
    assert!(
        original
            .apply_material_impulses(&[valid], limited, || false)
            .is_err()
    );
    let (empty, audit) = original
        .apply_material_impulses(&[], controls(), || false)
        .unwrap();
    same(&empty, &original);
    assert!(audit.receipts().is_empty());
    close(audit.external_work_j, 0.);
}
