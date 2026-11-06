//! Native articulated trajectories with geometry-dependent hydrostatic forcing.
use super::*;
use fs_mbd::articulated::motion::{FreeFloatingMotionError, FreeFloatingMotionState};
use fs_mbd_material::hydrostatic::{
    ArticulatedHydrostaticMotionError, FiniteHydrostaticLayer, FiniteHydrostaticTank,
    FiniteLayeredHydrostaticTank, FiniteTankControls, HydrostaticStepControls,
    HydrostaticTankGeometry,
};

fn controls(dt: f64, energy: f64) -> HydrostaticStepControls {
    HydrostaticStepControls {
        duration_s: dt,
        maximum_energy_defect_j: energy,
    }
}

fn still(z: f64, q: Vec<f64>) -> FreeFloatingMotionState {
    FreeFloatingMotionState {
        base: FreeFloatingBaseState::stationary(pose(z)),
        joint_velocity: vec![0.0; q.len()],
        joint_position: q,
    }
}

fn unit_box_source(density: f64) -> ResolvedArticulated {
    let card = density_card("hydrostatic-heave", density, density);
    let p = point(300.0);
    let regions = [RegionBinding {
        key: "solid",
        geometry: VolumeMoments::cuboid(Vec3::new(1.0, 1.0, 1.0)).unwrap(),
        material: &card,
        point: &p,
        selection: MaterialPropertySelection::SingleClaimOnly,
    }];
    resolve_articulated(
        &[MaterialLink {
            name: "floating-box",
            parent: None,
            parent_from_child_zero: Se3::identity(),
            joint: JointModel::FIXED,
            regions: &regions,
        }],
        1,
        || false,
    )
    .unwrap()
}

fn geometry() -> HydrostaticTankGeometry {
    HydrostaticTankGeometry {
        world_from_tank: Se3::identity(),
        inside_dimensions_m: Vec3::new(6.0, 2.0, 3.0),
    }
}
fn volume_controls() -> FiniteTankControls {
    FiniteTankControls {
        maximum_volume_defect_m3: 1e-13,
        maximum_volume_evaluations: 96,
    }
}

#[test]
fn native_trajectory_rises_floats_or_sinks_with_actual_density_ratio() {
    let water = fluid(1000.0);
    let medium = HydrostaticMedium::new(&water, gravity(), origin()).unwrap();
    for (density, acceleration) in [(500.0, 10.0), (1000.0, 0.0), (2000.0, -5.0)] {
        let source = mechanism(density, true);
        let meshes = [envelope(), envelope()];
        let prepared =
            PreparedArticulatedHydrostatics::new(&source, &meshes, 12, 1e-12, || false).unwrap();
        let state = still(-2.0, vec![0.31]);
        let dt = 0.1;
        let step = prepared
            .step_free_floating(&medium, &state, &[0.0], controls(dt, 1e-8), || false)
            .unwrap();
        let next = &step.dynamics.state_after;
        let dz = 0.5 * acceleration * dt * dt;
        assert!((next.base.world_from_base.translation().z - (-2.0 + dz)).abs() < 1e-11);
        assert!((next.base.twist_body.linear.z - acceleration * dt).abs() < 1e-10);
        assert!(
            next.base
                .twist_body
                .angular
                .dot(next.base.twist_body.angular)
                < 1e-20
        );
        assert!((next.joint_position[0] - 0.31).abs() < 1e-11);
        assert!(next.joint_velocity[0].abs() < 1e-10);
        close(step.dynamics.gravity_work, -0.32 * density * 10.0 * dz);
        close(
            step.body_gravitational_potential_change_j,
            -step.dynamics.gravity_work,
        );
        assert!(step.energy_defect_j.abs() < 1e-8);
        assert!(step.dynamics.linear_momentum_residual_world.norm() < 1e-10);
    }
}

#[test]
fn moving_immersion_is_sampled_at_native_midpoint_and_endpoint() {
    let source = unit_box_source(500.0);
    let meshes = [DisplacementMesh::cuboid(Vec3::new(1.0, 1.0, 1.0), || false).unwrap()];
    let prepared =
        PreparedArticulatedHydrostatics::new(&source, &meshes, 6, 1e-12, || false).unwrap();
    let water = fluid(1000.0);
    let medium = HydrostaticMedium::new(&water, gravity(), origin()).unwrap();
    let mut state = still(0.1, vec![]);
    state.base.twist_body = Twist::new(origin(), Vec3::new(0.0, 0.0, 0.3));
    let dt = 0.04;
    let step = prepared
        .step_free_floating(&medium, &state, &[], controls(dt, 1.0), || false)
        .unwrap();
    // m=500, waterplane=1, so z''=-20 z. Native midpoint geometry uses
    // z_h=z0+v0*dt/2; the accepted COM displacement uses its final momentum.
    let middle = 0.1 + 0.3 * dt / 2.0;
    let acceleration = -20.0 * middle;
    let velocity = 0.3 + acceleration * dt;
    let endpoint = 0.1 + 0.5 * (0.3 + velocity) * dt;
    close(step.load_before.link_loads[0].displaced_mass_kg(), 400.0);
    close(
        step.load_midpoint.link_loads[0].displaced_mass_kg(),
        1000.0 * (0.5 - middle),
    );
    close(
        step.load_after.link_loads[0].displaced_mass_kg(),
        1000.0 * (0.5 - endpoint),
    );
    close(
        step.dynamics
            .state_after
            .base
            .world_from_base
            .translation()
            .z,
        endpoint,
    );
    close(step.dynamics.state_after.base.twist_body.linear.z, velocity);
    let physical_energy_change = 250.0 * (velocity * velocity - 0.3_f64.powi(2))
        + 5000.0 * (endpoint * endpoint - 0.1_f64.powi(2));
    assert!((step.energy_defect_j - physical_energy_change).abs() < 1e-9);
}

#[test]
fn partial_immersion_heave_converges_to_the_analytic_period() {
    let source = unit_box_source(500.0);
    let meshes = [DisplacementMesh::cuboid(Vec3::new(1.0, 1.0, 1.0), || false).unwrap()];
    let prepared =
        PreparedArticulatedHydrostatics::new(&source, &meshes, 6, 1e-12, || false).unwrap();
    let water = fluid(1000.0);
    let medium = HydrostaticMedium::new(&water, gravity(), origin()).unwrap();
    let omega = 20.0_f64.sqrt();
    let mut errors = Vec::new();
    let mut energy_errors = Vec::new();
    for count in [25, 50, 100] {
        let mut state = still(0.1, vec![]);
        for _ in 0..count {
            let step = prepared
                .step_free_floating(
                    &medium,
                    &state,
                    &[],
                    controls(1.0 / count as f64, 1.0),
                    || false,
                )
                .unwrap();
            state = step.dynamics.state_after;
        }
        let z = state.base.world_from_base.translation().z;
        let v = state.base.twist_body.linear.z;
        errors.push((z - 0.1 * omega.cos()).abs() + (v + 0.1 * omega * omega.sin()).abs() / omega);
        energy_errors.push((250.0 * v * v + 5000.0 * z * z - 50.0).abs());
    }
    assert!(errors[2] < 8e-5, "{errors:?}");
    for pair in errors.windows(2) {
        assert!((3.8..4.2).contains(&(pair[0] / pair[1])), "{errors:?}");
    }
    for pair in energy_errors.windows(2) {
        assert!(
            (3.8..4.2).contains(&(pair[0] / pair[1])),
            "{energy_errors:?}"
        );
    }
}

#[test]
fn shared_finite_waterline_changes_the_linked_heave_frequency() {
    let source = mechanism(500.0, true);
    let meshes = [envelope(), envelope()];
    let prepared =
        PreparedArticulatedHydrostatics::new(&source, &meshes, 12, 1e-12, || false).unwrap();
    let water = fluid(1000.0);
    let tank =
        FiniteHydrostaticTank::new(&water, gravity(), geometry(), 9.6, volume_controls()).unwrap();
    // Two 80 kg links, each waterplane=0.4 m2. A=12 m2 and a_total=0.8:
    // omega^2=(rho*g*a/m) A/(A-a_total), stronger than prescribed-level 50.
    let equilibrium = (9.6 + 0.16) / 12.0;
    let omega = (50.0_f64 * 12.0 / 11.2).sqrt();
    let mut errors = Vec::new();
    for count in [25, 50] {
        let mut state = still(equilibrium + 0.025, vec![0.0]);
        for _ in 0..count {
            let step = prepared
                .step_free_floating_in_finite_tank(
                    &tank,
                    &state,
                    &[0.0],
                    controls(0.5 / count as f64, 0.1),
                    || false,
                )
                .unwrap();
            let after = &step.load_after.tank_load;
            close(
                12.0 * after.surface.surface_height_m - after.surface.displaced_volume_m3,
                tank.liquid_volume_m3(),
            );
            state = step.dynamics.state_after;
        }
        let z = state.base.world_from_base.translation().z - equilibrium;
        let v = state.base.twist_body.linear.z;
        errors.push(
            (z - 0.025 * (omega * 0.5).cos()).abs()
                + (v + 0.025 * omega * (omega * 0.5).sin()).abs() / omega,
        );
        assert!(state.joint_position[0].abs() < 1e-9);
        assert!(state.joint_velocity[0].abs() < 1e-8);
    }
    assert!((3.8..4.2).contains(&(errors[0] / errors[1])), "{errors:?}");
    assert!(errors[1] < 4e-5, "{errors:?}");
    close(tank.liquid_mass_kg(), 9600.0);
}

#[test]
fn finite_layered_motion_keeps_one_inventory_and_one_energy_account() {
    let source = mechanism(500.0, true);
    let meshes = [envelope(), envelope()];
    let prepared =
        PreparedArticulatedHydrostatics::new(&source, &meshes, 12, 1e-12, || false).unwrap();
    let dense = fluid(1200.0);
    let light = fluid(800.0);
    let tank = FiniteLayeredHydrostaticTank::new(
        &[
            FiniteHydrostaticLayer {
                source: &dense,
                liquid_volume_m3: 9.36,
            },
            FiniteHydrostaticLayer {
                source: &light,
                liquid_volume_m3: 1.44,
            },
        ],
        gravity(),
        geometry(),
        volume_controls(),
        || false,
    )
    .unwrap();
    let mut state = still(0.9, vec![0.0]);
    state.base.twist_body = Twist::new(origin(), Vec3::new(0.0, 0.0, 0.1));
    let dt = 0.01;
    let step = prepared
        .step_free_floating_in_finite_layered_tank(&tank, &state, &[0.0], controls(dt, 0.1), || {
            false
        })
        .unwrap();
    let center_mid = 0.9 + 0.1 * dt / 2.0;
    let bottom = center_mid - 0.2;
    let lower = (9.36 - 0.8 * bottom) / 11.2;
    let upper = (10.8 - 0.8 * bottom) / 11.2;
    let displaced_mass_per_link = 400.0 * 0.4 * (lower - bottom) + 800.0 * 0.4 * (upper - bottom);
    let acceleration = 10.0 * (displaced_mass_per_link / 80.0 - 1.0);
    close(step.load_midpoint.tank_load.surface.surface_height_m, upper);
    close(
        step.load_midpoint.tank_load.bodies[0].displaced_mass_kg(),
        displaced_mass_per_link,
    );
    close(
        step.dynamics.state_after.base.twist_body.linear.z,
        0.1 + acceleration * dt,
    );
    close(
        step.fluid_potential_change_j,
        step.load_after
            .tank_load
            .surface
            .liquid_gravitational_energy_change_j
            - step
                .load_before
                .tank_load
                .surface
                .liquid_gravitational_energy_change_j,
    );
    let masses: Vec<_> = tank.layer_masses_kg().collect();
    close(masses[0], 1200.0 * 9.36);
    close(masses[1], 800.0 * 1.44);
    assert!(step.energy_defect_j.abs() < 0.1);
}

#[test]
fn internal_actuation_supplies_work_and_retains_free_base_recoil() {
    let source = mechanism(1000.0, true);
    let meshes = [envelope(), envelope()];
    let prepared =
        PreparedArticulatedHydrostatics::new(&source, &meshes, 12, 1e-12, || false).unwrap();
    let water = fluid(1000.0);
    let medium = HydrostaticMedium::new(&water, gravity(), origin()).unwrap();
    let state = still(-2.0, vec![0.31]);
    let torque = 5.0;
    let step = prepared
        .step_free_floating(&medium, &state, &[torque], controls(0.02, 1e-7), || false)
        .unwrap();
    assert!(step.dynamics.state_after.joint_velocity[0] > 0.0);
    assert!(step.dynamics.state_after.base.twist_body.angular.y < 0.0);
    assert!(step.dynamics.actuator_work > 1e-6);
    assert!(
        (step.dynamics.actuator_work
            - torque * (step.dynamics.state_after.joint_position[0] - 0.31))
            .abs()
            < 1e-10
    );
    assert!(step.dynamics.external_linear_impulse_world.norm() < 1e-10);
    assert!(step.dynamics.external_angular_impulse_world.norm() < 1e-9);
    assert!(step.energy_defect_j.abs() < 1e-7);
}

#[test]
fn energy_and_endpoint_admission_and_cancellation_are_transactional() {
    let source = unit_box_source(500.0);
    let meshes = [DisplacementMesh::cuboid(Vec3::new(1.0, 1.0, 1.0), || false).unwrap()];
    let prepared =
        PreparedArticulatedHydrostatics::new(&source, &meshes, 6, 1e-12, || false).unwrap();
    let water = fluid(1000.0);
    let medium = HydrostaticMedium::new(&water, gravity(), origin()).unwrap();
    let state = still(0.1, vec![]);
    let original = state.clone();
    assert!(matches!(
        prepared.step_free_floating(&medium, &state, &[], controls(0.3, 1e-6), || false),
        Err(ArticulatedHydrostaticMotionError::EnergyDefect { .. })
    ));
    let mut polls = 0;
    let expected = prepared
        .step_free_floating(&medium, &state, &[], controls(0.001, 1e-3), || {
            polls += 1;
            false
        })
        .unwrap();
    for stop in [1, polls / 2, polls] {
        let mut seen = 0;
        let result =
            prepared.step_free_floating(&medium, &state, &[], controls(0.001, 1e-3), || {
                seen += 1;
                seen == stop
            });
        assert!(matches!(
            result,
            Err(ArticulatedHydrostaticMotionError::Load(
                ArticulatedHydrostaticError::Physics(HydrostaticError::Cancelled)
            )) | Err(ArticulatedHydrostaticMotionError::Motion(
                FreeFloatingMotionError::Cancelled
            ))
        ));
        assert_eq!(
            prepared
                .step_free_floating(&medium, &state, &[], controls(0.001, 1e-3), || false,)
                .unwrap(),
            expected
        );
    }
    assert_eq!(state, original);

    let source = mechanism(2000.0, true);
    let meshes = [envelope(), envelope()];
    let prepared =
        PreparedArticulatedHydrostatics::new(&source, &meshes, 12, 1e-12, || false).unwrap();
    let tank =
        FiniteHydrostaticTank::new(&water, gravity(), geometry(), 9.6, volume_controls()).unwrap();
    let mut sinking = still(0.35, vec![0.0]);
    sinking.base.twist_body = Twist::new(origin(), Vec3::new(0.0, 0.0, -1.0));
    let before = sinking.clone();
    // Start z=.35 and midpoint z=.25 both clear the .2 m lower hull arm.
    // The candidate endpoint z=.05 crosses the floor and must not publish.
    assert!(matches!(
        prepared.step_free_floating_in_finite_tank(
            &tank,
            &sinking,
            &[0.0],
            controls(0.2, 1.0),
            || false,
        ),
        Err(ArticulatedHydrostaticMotionError::Load(
            ArticulatedHydrostaticError::Physics(HydrostaticError::Invalid(_))
        ))
    ));
    assert_eq!(sinking, before);
    assert!(
        prepared
            .step_free_floating_in_finite_tank(
                &tank,
                &sinking,
                &[0.0],
                controls(0.05, 1e-8),
                || false,
            )
            .is_ok()
    );
}
