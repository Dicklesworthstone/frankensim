//! G1/G3: source heat changes continuing momentum-form motion, not a new preset.
#[path = "heating/common.rs"]
mod common;
use common::*;
use fs_ga::{Se3, So3, Twist, Vec3};
use fs_material::phase::SolidLiquidPhase;
use fs_mbd::{Gravity, RigidBodyIntegrator};
use fs_mbd_material::principal::transfer::ConfigurationTransferControls;
use fs_mbd_material::thermal::heating::HeatedBody;
use fs_mbd_material::thermal::heating::motion::{
    HeatedMotionControls, HeatedMotionError, HeatedRigidBody,
};
use fs_thermal::enthalpy::*;

fn controls() -> HeatedMotionControls {
    HeatedMotionControls {
        maximum_step_s: 0.0625,
        relative_tensor_tolerance: 1e-12,
        transfer: ConfigurationTransferControls {
            relative_momentum_tolerance: 1e-12,
            maximum_parameter_work_j: 10.,
        },
        maximum_numerical_energy_change_j: 1e-7,
    }
}
fn same(a: fs_mbd::Vec3, b: fs_mbd::Vec3) {
    close(a.x, b.x);
    close(a.y, b.y);
    close(a.z, b.z);
}
fn angular_speed(body: &HeatedRigidBody<'_>) -> f64 {
    body.principal_frame()
        .properties()
        .angular_velocity_body_checked(body.state().angular_momentum_body())
        .unwrap()
        .norm_squared()
        .sqrt()
}

#[test]
fn heat_capacity_swap_changes_the_live_spin_with_unchanged_mass_and_world_momentum() {
    let mut speeds = Vec::new();
    for cp in [10., 20.] {
        let material = card(cp, 1e-3, 500.);
        let reference = reference(&material, 1);
        let chart = caloric(&material, 1e5, SolidLiquidPhase::Solid);
        let body = heated(&reference, &material, &chart, 40., vec![]);
        let moving = HeatedRigidBody::new(
            body,
            Se3::identity(),
            Twist::new(Vec3::new(0., 0., 5.), Vec3::new(1., -2., 0.5)),
            RigidBodyIntegrator::new(Gravity::ZERO),
            controls(),
            || Ok(()),
        )
        .unwrap();
        let initial = moving.diagnostics().unwrap();
        let (hot, audit) = moving.advance_to(2., config(), || Ok(())).unwrap();
        let stretch = 1. + 1e-3 * 80. / (reference.inertia().mass() * cp);
        let expected_spin = 5. / (stretch * stretch);
        close(angular_speed(&hot), expected_spin);
        close(hot.thermal().time_s(), 2.);
        close(
            hot.thermal().temperature_k(CellId(0)).unwrap(),
            300. + (stretch - 1.) / 1e-3,
        );
        assert_eq!(
            hot.configuration().inertia().mass(),
            reference.inertia().mass()
        );
        assert_eq!(
            hot.state().linear_momentum_world(),
            moving.state().linear_momentum_world()
        );
        same(
            hot.diagnostics().unwrap().angular_momentum_world,
            initial.angular_momentum_world,
        );
        same(
            hot.state().pose().position_world(),
            fs_mbd::Vec3::new(2., -4., 1.),
        );
        close(audit.thermal.external_heat_j, 80.);
        close(
            audit.parameter_work_j,
            initial.rotational_kinetic_energy * (1. / (stretch * stretch) - 1.),
        );
        assert!(audit.parameter_work_j < 0.);
        close(audit.absolute_parameter_work_j, -audit.parameter_work_j);
        close(audit.combined_residual_j, 0.);
        close(audit.numerical_energy_change_j, 0.);
        assert_eq!(audit.thermal.substeps, 32);
        close(moving.thermal().time_s(), 0.);
        close(angular_speed(&moving), 5.);
        speeds.push(angular_speed(&hot));
    }
    assert!(speeds[0] < speeds[1]);
}

#[test]
fn nonuniform_expansion_moves_the_body_origin_without_teleporting_its_mass_centre() {
    let material = card(10., 1e-3, 500.);
    let reference = reference(&material, 2);
    let chart = caloric(&material, 1e5, SolidLiquidPhase::Solid);
    let mut warm = region(&material, &chart, "a", 0, 40.);
    warm.fixed_point_m = Vec3::new(-1., -2., 1.);
    let body = HeatedBody::new(
        &reference,
        vec![warm, region(&material, &chart, "b", 1, 0.)],
        vec![],
        budget(),
        config(),
        || Ok(()),
    )
    .unwrap();
    let moving = HeatedRigidBody::new(
        body,
        Se3::identity(),
        Twist::new(Vec3::default(), Vec3::new(1., 0., 0.)),
        RigidBodyIntegrator::new(Gravity::ZERO),
        controls(),
        || Ok(()),
    )
    .unwrap();
    let (hot, audit) = moving.advance_to(1., config(), || Ok(())).unwrap();
    same(
        hot.state().pose().position_world(),
        fs_mbd::Vec3::new(1.5, 0., 0.),
    );
    let pose = hot.world_from_body().unwrap();
    for (a, b) in pose
        .rotation()
        .matrix()
        .m
        .into_iter()
        .zip(So3::identity().matrix().m)
    {
        close(a, b);
    }
    let centre = pose.translation()
        + pose
            .rotation()
            .rotate(hot.configuration().inertia().center_of_mass())
            .unwrap();
    close(centre.x, 1.5);
    close(centre.y, 0.);
    close(centre.z, 0.);
    assert!(pose.translation().y.abs() > 1e-4);
    assert!(hot.configuration().inertia().inertia_com().m[1].abs() > 1e-4);
    close(angular_speed(&hot), 0.);
    close(audit.parameter_work_j, 0.);
    close(audit.combined_residual_j, 0.);
}

#[test]
fn cooling_restart_and_refusals_keep_thermal_and_motion_histories_together() {
    let material = card(10., 1e-3, 500.);
    let reference = reference(&material, 1);
    let chart = caloric(&material, 1e5, SolidLiquidPhase::Solid);
    let events = vec![TimedEvent {
        time_s: 1.,
        event: Event::Boundary(CellId(0), boundary(-40.)),
    }];
    let body = heated(&reference, &material, &chart, 40., events);
    let moving = HeatedRigidBody::new(
        body.clone(),
        Se3::identity(),
        Twist::new(Vec3::new(0., 0., 5.), Vec3::default()),
        RigidBodyIntegrator::new(Gravity::ZERO),
        controls(),
        || Ok(()),
    )
    .unwrap();
    let (whole, full_audit) = moving.advance_to(2., config(), || Ok(())).unwrap();
    let (hot, heat_audit) = moving.advance_to(1., config(), || Ok(())).unwrap();
    let (cool, cool_audit) = hot.advance_to(2., config(), || Ok(())).unwrap();
    assert_eq!(whole.state(), cool.state());
    assert_eq!(whole.thermal().cells(), cool.thermal().cells());
    close(angular_speed(&whole), 5.);
    assert!(angular_speed(&hot) < 5.);
    assert!(heat_audit.parameter_work_j < 0. && cool_audit.parameter_work_j > 0.);
    close(
        full_audit.parameter_work_j,
        heat_audit.parameter_work_j + cool_audit.parameter_work_j,
    );
    assert!(full_audit.absolute_parameter_work_j > 1e-3);
    close(full_audit.combined_residual_j, 0.);
    let saved = moving.state();
    assert!(
        moving
            .advance_to(1., config(), || Err(DepositionError::Cancelled))
            .is_err()
    );
    assert_eq!(moving.state(), saved);
    let mut constrained = controls();
    constrained.transfer.maximum_parameter_work_j = 1e-5;
    let limited = HeatedRigidBody::new(
        body,
        Se3::identity(),
        Twist::new(Vec3::new(0., 0., 5.), Vec3::default()),
        RigidBodyIntegrator::new(Gravity::ZERO),
        constrained,
        || Ok(()),
    )
    .unwrap();
    assert!(limited.advance_to(2., config(), || Ok(())).is_err());
    assert_eq!(limited.state(), saved);
    close(limited.thermal().energy_j(), 0.);
}

#[test]
fn unsupported_hot_intermediate_state_refuses_the_whole_motion_candidate() {
    let material = card(10., 1e-3, 301.);
    let reference = reference(&material, 1);
    let chart = caloric(&material, 1e5, SolidLiquidPhase::Solid);
    let body = heated(
        &reference,
        &material,
        &chart,
        40.,
        vec![TimedEvent {
            time_s: 1.,
            event: Event::Boundary(CellId(0), boundary(-40.)),
        }],
    );
    let moving = HeatedRigidBody::new(
        body,
        Se3::identity(),
        Twist::new(Vec3::new(0., 0., 5.), Vec3::default()),
        RigidBodyIntegrator::new(Gravity::ZERO),
        controls(),
        || Ok(()),
    )
    .unwrap();
    assert!(matches!(
        moving.advance_to(2., config(), || Ok(())),
        Err(HeatedMotionError::Heating(_))
    ));
    close(moving.thermal().time_s(), 0.);
    close(angular_speed(&moving), 5.);
}

#[test]
fn unchanged_material_configuration_retains_the_original_integrator_bit_for_bit() {
    let material = card(10., 0., 500.);
    let reference = reference(&material, 1);
    let chart = caloric(&material, 1e5, SolidLiquidPhase::Solid);
    let body = heated(&reference, &material, &chart, 0., vec![]);
    let integrator =
        RigidBodyIntegrator::new(Gravity::new(fs_mbd::Vec3::new(0., -1., 0.)).unwrap());
    let moving = HeatedRigidBody::new(
        body,
        Se3::identity(),
        Twist::new(Vec3::new(0.2, 0.3, -0.5), Vec3::new(1., 0., 0.)),
        integrator,
        controls(),
        || Ok(()),
    )
    .unwrap();
    let mut expected = moving.state();
    for _ in 0..16 {
        expected = integrator
            .step(
                expected,
                moving.principal_frame().properties(),
                fs_mbd::Wrench::ZERO,
                0.0625,
            )
            .unwrap()
            .state_after;
    }
    let (actual, audit) = moving.advance_to(1., config(), || Ok(())).unwrap();
    assert_eq!(actual.state(), expected);
    assert_eq!(audit.parameter_work_j, 0.);
    close(audit.combined_residual_j, 0.);
}

#[derive(Clone)]
struct GeneratedHeat {
    work_j: f64,
}
impl DepositionObserver for GeneratedHeat {
    type Error = DepositionError;
    fn sources(
        &mut self,
        _: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<std::collections::BTreeMap<CellId, f64>, Self::Error> {
        Ok(std::collections::BTreeMap::from([(CellId(0), 40.)]))
    }
    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        Ok(())
    }
    fn step(
        &mut self,
        _: f64,
        dt: f64,
        _: &std::collections::BTreeMap<CellId, f64>,
        _: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        self.work_j += 40. * dt;
        Ok(())
    }
}

#[test]
fn generated_heat_uses_the_same_motion_steps_and_is_counted_once() {
    let material = card(10., 1e-3, 500.);
    let reference = reference(&material, 1);
    let chart = caloric(&material, 1e5, SolidLiquidPhase::Solid);
    let body = heated(&reference, &material, &chart, 0., vec![]);
    let moving = HeatedRigidBody::new(
        body,
        Se3::identity(),
        Twist::new(Vec3::new(0., 0., 5.), Vec3::default()),
        RigidBodyIntegrator::new(Gravity::ZERO),
        controls(),
        || Ok(()),
    )
    .unwrap();
    let source = GeneratedHeat { work_j: 0. };
    let (heated, source_after, audit) = moving
        .advance_with(&source, 1., config(), || Ok(()))
        .unwrap();
    close(heated.thermal().temperature_k(CellId(0)).unwrap(), 302.);
    close(angular_speed(&heated), 5. / 1.002_f64.powi(2));
    close(source_after.work_j, 40.);
    close(audit.thermal.external_heat_j, 40.);
    close(audit.thermal.stored_energy_change_j, 40.);
    close(source.work_j, 0.);
    close(audit.combined_residual_j, 0.);
}
