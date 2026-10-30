//! G1/G3/G4: original material heat changes buoyancy and continuing motion.
#[path = "heating/common.rs"]
mod common;
use common::*;
use fs_ga::{Se3, So3, Twist, Vec3};
use fs_material::{fluid_transport::FluidTransportState, phase::SolidLiquidPhase};
use fs_mbd::{Gravity, RigidBodyIntegrator};
use fs_mbd_material::hydrostatic::{DisplacementMesh, HydrostaticMedium};
use fs_mbd_material::principal::transfer::ConfigurationTransferControls;
use fs_mbd_material::thermal::heating::HeatedBody;
use fs_mbd_material::thermal::heating::motion::hydrostatic::ThermalDisplacementRegion;
use fs_mbd_material::thermal::heating::motion::{HeatedMotionControls, HeatedRigidBody};
use fs_thermal::enthalpy::*;

const DT: f64 = 0.0625;
fn gravity() -> Gravity {
    Gravity::new(fs_mbd::Vec3::new(0., 0., -10.)).unwrap()
}
fn controls() -> HeatedMotionControls {
    HeatedMotionControls {
        maximum_step_s: DT,
        relative_tensor_tolerance: 1e-12,
        transfer: ConfigurationTransferControls {
            relative_momentum_tolerance: 1e-12,
            maximum_parameter_work_j: 10.,
        },
        maximum_numerical_energy_change_j: 1e-3,
    }
}
fn fluid(rho: f64) -> FluidTransportState {
    FluidTransportState::try_new(rho, 4000., 0.5, 0.001).unwrap()
}
fn envelope(key: &str, x: f64) -> ThermalDisplacementRegion {
    let cube = DisplacementMesh::cuboid(Vec3::new(1., 1., 1.), || false).unwrap();
    ThermalDisplacementRegion {
        region: key.into(),
        envelope: DisplacementMesh::new(
            cube.cells()
                .iter()
                .map(|cell| cell.map(|p| p + Vec3::new(x, 0., 0.)))
                .collect(),
            6,
            || false,
        )
        .unwrap(),
    }
}
fn moving(body: HeatedBody<'_>, z: f64, controls: HeatedMotionControls) -> HeatedRigidBody<'_> {
    HeatedRigidBody::new(
        body,
        Se3::from_parts(So3::identity(), Vec3::new(0., 0., z)).unwrap(),
        Twist::zero(),
        RigidBodyIntegrator::new(gravity()),
        controls,
        || Ok(()),
    )
    .unwrap()
}

#[test]
fn expansion_data_alone_changes_displacement_and_actual_flotation_velocity() {
    let reservoir = fluid(4.);
    let medium = HydrostaticMedium::new(&reservoir, gravity(), Vec3::default()).unwrap();
    for alpha in [0., 0.01] {
        let material = card(10., alpha, 500.);
        let source = reference(&material, 1);
        let chart = caloric(&material, 1e5, SolidLiquidPhase::Solid);
        let body = moving(
            heated(&source, &material, &chart, 40., vec![]),
            0.,
            controls(),
        );
        let parts = [envelope("a", 0.)];
        let (hot, audit) = body
            .advance_hydrostatic(&medium, &parts, 2. * DT, config(), || Ok(()))
            .unwrap();
        // Initial half-submerged cube is exactly neutral. First heat update
        // changes its volume. The SECOND native kick sees that real buoyancy.
        let stretch_first = 1. + alpha * 40. * DT / (2. * 10.);
        let velocity = DT * 10. * (stretch_first.powi(3) - 1.);
        close(
            hot.state()
                .center_of_mass_velocity_world(hot.principal_frame().properties())
                .unwrap()
                .z,
            velocity,
        );
        let stretch = 1. + alpha * 40. * (2. * DT) / (2. * 10.);
        let z = hot.world_from_body().unwrap().translation().z;
        let depth = stretch / 2. - z;
        close(
            audit.load_after.displaced().volume_m3(),
            stretch * stretch * depth,
        );
        close(
            audit.load_after.pressure_potential_j(),
            20. * stretch * stretch * depth * depth,
        );
        close(
            audit.load_after.wrench().force_world.z,
            40. * stretch * stretch * depth,
        );
        close(
            hot.configuration().inertia().mass(),
            source.inertia().mass(),
        );
        close(audit.thermal.external_heat_j, 40. * 2. * DT);
        close(audit.combined_residual_j, 0.);
        close(audit.parameter_work_j, audit.displacement_parameter_work_j);
        assert_eq!(audit.displacement_parameter_work_j > 0., alpha > 0.);
        close(body.thermal().time_s(), 0.);
    }
}

#[test]
fn independently_heated_regions_shift_the_actual_center_of_buoyancy_and_torque() {
    let material = card(10., 0.01, 500.);
    // A gap keeps different free strains geometrically compatible. The
    // declared rigid connector has negligible material/displacement volume.
    let point = point(300., 1e5);
    let source = fs_mbd_material::resolve_body(
        &[("a", 0.), ("b", 2.)].map(|(key, x)| fs_mbd_material::RegionBinding {
            key,
            material: &material,
            point: &point,
            selection: fs_material::state_point::MaterialPropertySelection::SingleClaimOnly,
            geometry: fs_mbd_material::VolumeMoments::cuboid(Vec3::new(1., 1., 1.))
                .unwrap()
                .transformed(Se3::from_parts(So3::identity(), Vec3::new(x, 0., 0.)).unwrap())
                .unwrap(),
        }),
        2,
        || false,
    )
    .unwrap();
    let chart = caloric(&material, 1e5, SolidLiquidPhase::Solid);
    let body = HeatedBody::new(
        &source,
        vec![
            region(&material, &chart, "a", 0, 40.),
            region(&material, &chart, "b", 1, 0.),
        ],
        vec![],
        budget(),
        config(),
        || Ok(()),
    )
    .unwrap();
    let body = moving(body, -2., controls());
    let reservoir = fluid(2.);
    let medium = HydrostaticMedium::new(&reservoir, gravity(), Vec3::default()).unwrap();
    // Reverse authoring order proves each part follows its named source region.
    let parts = [envelope("b", 2.), envelope("a", 0.)];
    let (hot, audit) = body
        .advance_hydrostatic(&medium, &parts, DT, config(), || Ok(()))
        .unwrap();
    let stretch = 1. + 0.01 * 40. * DT / 20.;
    close(
        audit.load_after.displaced().volume_m3(),
        1. + stretch.powi(3),
    );
    close(
        audit.load_after.center_of_buoyancy_world().unwrap().x,
        2. / (1. + stretch.powi(3)),
    );
    let torque_world = hot
        .state()
        .pose()
        .orientation()
        .rotate_body_to_world(audit.load_after.wrench().torque_body);
    close(torque_world.y, 20. * (stretch.powi(3) - 1.));
    close(audit.combined_residual_j, 0.);
    close(hot.thermal().temperature_k(CellId(1)).unwrap(), 300.);
    assert_eq!(hot.configuration().inertia().mass(), 4.);
}

#[test]
fn source_fluid_density_reaches_the_native_heated_body_step() {
    let material = card(10., 0., 500.);
    let source = reference(&material, 1);
    let chart = caloric(&material, 1e5, SolidLiquidPhase::Solid);
    let mut limits = controls();
    limits.maximum_numerical_energy_change_j = 2.;
    let body = moving(heated(&source, &material, &chart, 0., vec![]), 0., limits);
    let parts = [envelope("a", 0.)];
    for density in [4., 8.] {
        let reservoir = fluid(density);
        let medium = HydrostaticMedium::new(&reservoir, gravity(), Vec3::default()).unwrap();
        let (next, audit) = body
            .advance_hydrostatic(&medium, &parts, DT, config(), || Ok(()))
            .unwrap();
        let expected_v = DT * (density * 10. * 0.5 / 2. - 10.);
        close(
            next.state()
                .center_of_mass_velocity_world(next.principal_frame().properties())
                .unwrap()
                .z,
            expected_v,
        );
        close(audit.parameter_work_j, 0.);
        close(audit.combined_residual_j, 0.);
    }
}

#[test]
fn restart_and_late_refusals_preserve_one_accepted_thermal_mechanical_history() {
    let material = card(10., 0.01, 500.);
    let source = reference(&material, 1);
    let chart = caloric(&material, 1e5, SolidLiquidPhase::Solid);
    let thermal = heated(&source, &material, &chart, 40., vec![]);
    let body = moving(thermal.clone(), 0., controls());
    let reservoir = fluid(4.);
    let medium = HydrostaticMedium::new(&reservoir, gravity(), Vec3::default()).unwrap();
    let parts = [envelope("a", 0.)];
    let (whole, _) = body
        .advance_hydrostatic(&medium, &parts, 4. * DT, config(), || Ok(()))
        .unwrap();
    let (half, _) = body
        .advance_hydrostatic(&medium, &parts, 2. * DT, config(), || Ok(()))
        .unwrap();
    let (split, _) = half
        .advance_hydrostatic(&medium, &parts, 4. * DT, config(), || Ok(()))
        .unwrap();
    assert_eq!(whole.state(), split.state());
    assert_eq!(whole.thermal().cells(), split.thermal().cells());
    let mut calls = 0;
    body.advance_hydrostatic(&medium, &parts, DT, config(), || {
        calls += 1;
        Ok(())
    })
    .unwrap();
    let mut visited = 0;
    assert!(
        body.advance_hydrostatic(&medium, &parts, DT, config(), || {
            visited += 1;
            if visited == calls {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        })
        .is_err()
    );
    let mut tight = controls();
    tight.transfer.maximum_parameter_work_j = 1e-8;
    let limited = moving(thermal, 0., tight);
    assert!(
        limited
            .advance_hydrostatic(&medium, &parts, DT, config(), || Ok(()))
            .is_err()
    );
    assert!(
        body.advance_hydrostatic(&medium, &[envelope("missing", 0.)], DT, config(), || Ok(()))
            .is_err()
    );
    assert!(
        body.advance_hydrostatic(
            &medium,
            &[envelope("a", 0.), envelope("a", 0.)],
            DT,
            config(),
            || Ok(())
        )
        .is_err()
    );
    close(body.thermal().time_s(), 0.);
    close(limited.thermal().energy_j(), 0.);
    assert_eq!(limited.state(), body.state());
}
