//! Shared finite-liquid inventories coupled through native articulated dynamics.
use super::*;
use fs_mbd_material::hydrostatic::{
    ArticulatedFiniteTankLoad, FiniteHydrostaticLayer, FiniteHydrostaticTank,
    FiniteLayeredHydrostaticTank, FiniteTankControls, HydrostaticTankGeometry,
};

fn controls() -> FiniteTankControls {
    FiniteTankControls {
        maximum_volume_defect_m3: 1e-13,
        maximum_volume_evaluations: 96,
    }
}

fn geometry(height: f64) -> HydrostaticTankGeometry {
    HydrostaticTankGeometry {
        world_from_tank: Se3::identity(),
        inside_dimensions_m: Vec3::new(6.0, 2.0, height),
    }
}

fn sliders() -> ResolvedArticulated {
    let card = density_card("finite-tank-sliders", 1500.0, 1500.0);
    let p = point(300.0);
    let regions = [RegionBinding {
        key: "solid",
        geometry: VolumeMoments::cuboid(Vec3::new(1.0, 0.4, 0.4))
            .unwrap()
            .transformed(shifted(0.75))
            .unwrap(),
        material: &card,
        point: &p,
        selection: MaterialPropertySelection::SingleClaimOnly,
    }];
    let joint = JointModel::prismatic(Vec3::new(0.0, 0.0, 1.0), None).unwrap();
    resolve_articulated(
        &[
            MaterialLink {
                name: "lower-slider",
                parent: None,
                parent_from_child_zero: Se3::identity(),
                joint,
                regions: &regions,
            },
            MaterialLink {
                name: "upper-slider",
                parent: Some(0),
                parent_from_child_zero: shifted(1.5),
                joint,
                regions: &regions,
            },
        ],
        2,
        || false,
    )
    .unwrap()
}

// Both boxes intersect the waterline: each horizontal section has area 0.4.
// Their lower faces are at b[0], b[1], in the same 12 m2 tank.
fn level(volume: f64, b: [f64; 2]) -> f64 {
    (volume - 0.4 * (b[0] + b[1])) / (12.0 - 0.8)
}

// Integrate z over the actual liquid = tank minus the two submerged boxes.
// This is independent of the pressure-depth/cross-displacement energy formula.
fn liquid_moment(h: f64, b: [f64; 2]) -> f64 {
    6.0 * h * h - b.iter().map(|z| 0.2 * (h * h - z * z)).sum::<f64>()
}

#[test]
fn one_slider_changes_neighbor_buoyancy_and_native_holding_force() {
    let source = sliders();
    let meshes = [envelope(), envelope()];
    let prepared =
        PreparedArticulatedHydrostatics::new(&source, &meshes, 12, 1e-12, || false).unwrap();
    let water = fluid(1000.0);
    let tank =
        FiniteHydrostaticTank::new(&water, gravity(), geometry(2.0), 9.6, controls()).unwrap();
    let base = pose(0.9);
    let initial = prepared
        .load_finite_tank(&tank, base, &[0.0, 0.0], || false)
        .unwrap();
    let q = [0.0, 0.05];
    let load = prepared
        .load_finite_tank(&tank, base, &q, || false)
        .unwrap();
    let b = [0.7, 0.75];
    let h = level(9.6, b);
    close(load.tank_load.surface.surface_height_m, h);
    close(
        initial.tank_load.surface.surface_height_m - h,
        0.4 * 0.05 / 11.2,
    );
    assert!(load.external_body_wrenches[0].force.z < initial.external_body_wrenches[0].force.z);
    let forces = b.map(|z| 10000.0 * 0.4 * (h - z));
    for (wrench, force) in load.external_body_wrenches.iter().zip(forces) {
        close(wrench.force.z, force);
        close(wrench.torque.y, -0.75 * force);
    }
    close(
        load.tank_load.surface.displaced_volume_m3,
        0.4 * (2.0 * h - b[0] - b[1]),
    );
    close(
        load.tank_load.surface.liquid_gravitational_energy_change_j,
        10000.0 * (liquid_moment(h, b) - 6.0 * 0.8_f64.powi(2)),
    );
    assert!(load.tank_load.surface.liquid_volume_defect_m3.abs() <= 1e-13);
    close(tank.liquid_mass_kg(), 9600.0);
    close(tank.liquid_volume_m3(), 9.6);

    // q0 lifts both 240 kg links; q1 lifts only the child. The analytically
    // balanced efforts must hold both original native prismatic coordinates.
    let holding = [4800.0 - forces[0] - forces[1], 2400.0 - forces[1]];
    let result = forward_dynamics(
        source.model(),
        BaseState::stationary(base),
        &q,
        &[0.0; 2],
        &holding,
        Vec3::new(0.0, 0.0, -10.0),
        &load.external_body_wrenches,
    )
    .unwrap();
    assert!(
        result
            .generalized_acceleration
            .iter()
            .all(|a| a.abs() < 1e-10)
    );
}

#[test]
fn layered_articulation_retains_one_stack_energy_and_distinct_displaced_masses() {
    let source = sliders();
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
        geometry(2.0),
        controls(),
        || false,
    )
    .unwrap();
    let load = prepared
        .load_finite_layered_tank(&tank, pose(0.9), &[0.0, 0.05], || false)
        .unwrap();
    let b = [0.7, 0.75];
    let lower = level(9.36, b);
    let upper = level(10.8, b);
    close(load.tank_load.surface.surface_height_m, upper);
    for (index, z) in b.into_iter().enumerate() {
        let lower_volume = 0.4 * (lower - z);
        let upper_volume = 0.4 * (upper - z);
        let mass = 400.0 * lower_volume + 800.0 * upper_volume;
        let first_moment =
            400.0 * lower_volume * (lower + z) / 2.0 + 800.0 * upper_volume * (upper + z) / 2.0;
        close(load.tank_load.bodies[index].displaced_mass_kg(), mass);
        close(
            load.tank_load.bodies[index].displaced().volume_m3(),
            upper_volume,
        );
        close(
            load.tank_load.bodies[index]
                .center_of_buoyancy_world()
                .unwrap()
                .z,
            first_moment / mass,
        );
        close(load.external_body_wrenches[index].force.z, 10.0 * mass);
        close(load.external_body_wrenches[index].torque.y, -7.5 * mass);
    }
    close(
        load.tank_load.surface.liquid_gravitational_energy_change_j,
        4000.0 * (liquid_moment(lower, b) - 6.0 * 0.78_f64.powi(2))
            + 8000.0 * (liquid_moment(upper, b) - 6.0 * 0.9_f64.powi(2)),
    );
    let masses: Vec<_> = tank.layer_masses_kg().collect();
    close(masses[0], 1200.0 * 9.36);
    close(masses[1], 800.0 * 1.44);
    close(tank.liquid_volume_m3(), 10.8);
}

fn check_virtual_work(
    prepared: &PreparedArticulatedHydrostatics<'_>,
    base: Se3,
    evaluate: impl Fn(&[f64]) -> ArticulatedFiniteTankLoad,
) {
    let q = [0.18, -0.29];
    let load = evaluate(&q);
    let holding = inverse_dynamics(
        prepared.source().model(),
        BaseState::stationary(base),
        &q,
        &[0.0; 2],
        &[0.0; 2],
        origin(),
        &load.external_body_wrenches,
    )
    .unwrap();
    let step = 2e-6;
    for axis in 0..2 {
        let mut plus = q;
        let mut minus = q;
        plus[axis] += step;
        minus[axis] -= step;
        let energy = |q: &[f64]| {
            evaluate(q)
                .tank_load
                .surface
                .liquid_gravitational_energy_change_j
        };
        let derivative = (energy(&plus) - energy(&minus)) / (2.0 * step);
        assert!(
            (derivative - holding.generalized_force[axis]).abs() < 5e-5,
            "joint {axis}: liquid gradient {derivative}, holding {}",
            holding.generalized_force[axis],
        );
    }
}

#[test]
fn moving_common_surfaces_preserve_hinge_virtual_work() {
    let source = mechanism(1500.0, false);
    let meshes = [envelope(), envelope()];
    let prepared =
        PreparedArticulatedHydrostatics::new(&source, &meshes, 12, 1e-12, || false).unwrap();
    let water = fluid(1000.0);
    let uniform =
        FiniteHydrostaticTank::new(&water, gravity(), geometry(3.0), 10.8, controls()).unwrap();
    let base = pose(1.05);
    check_virtual_work(&prepared, base, |q| {
        prepared
            .load_finite_tank(&uniform, base, q, || false)
            .unwrap()
    });
    let dense = fluid(1200.0);
    let light = fluid(800.0);
    let layered = FiniteLayeredHydrostaticTank::new(
        &[
            FiniteHydrostaticLayer {
                source: &dense,
                liquid_volume_m3: 9.0,
            },
            FiniteHydrostaticLayer {
                source: &light,
                liquid_volume_m3: 1.8,
            },
        ],
        gravity(),
        geometry(3.0),
        controls(),
        || false,
    )
    .unwrap();
    check_virtual_work(&prepared, base, |q| {
        prepared
            .load_finite_layered_tank(&layered, base, q, || false)
            .unwrap()
    });
}

#[test]
fn free_floating_finite_tank_keeps_weight_and_buoyancy_separate() {
    let water = fluid(1000.0);
    let tank =
        FiniteHydrostaticTank::new(&water, gravity(), geometry(4.0), 36.0, controls()).unwrap();
    for (density, acceleration) in [(500.0, 10.0), (1000.0, 0.0), (2000.0, -5.0)] {
        let source = mechanism(density, true);
        let meshes = [envelope(), envelope()];
        let prepared =
            PreparedArticulatedHydrostatics::new(&source, &meshes, 12, 1e-12, || false).unwrap();
        let base_pose = pose(1.5);
        let load = prepared
            .load_finite_tank(&tank, base_pose, &[0.31], || false)
            .unwrap();
        close(
            load.tank_load.surface.surface_height_m,
            (36.0 + 0.32) / 12.0,
        );
        let result = free_floating_forward_dynamics(
            source.model(),
            FreeFloatingBaseState::stationary(base_pose),
            &[0.31],
            &[0.0],
            &[0.0],
            Vec3::new(0.0, 0.0, -10.0),
            &load.external_body_wrenches,
        )
        .unwrap();
        let base = result.base_spatial_acceleration_body;
        assert!((base.linear.z - acceleration).abs() < 1e-10);
        assert!(base.angular.dot(base.angular) < 1e-20);
        assert!(result.generalized_acceleration[0].abs() < 1e-10);
    }
}

#[test]
fn overflow_floor_and_late_cancellation_return_no_shared_load() {
    let source = sliders();
    let meshes = [envelope(), envelope()];
    let prepared =
        PreparedArticulatedHydrostatics::new(&source, &meshes, 12, 1e-12, || false).unwrap();
    let water = fluid(1000.0);
    let tank =
        FiniteHydrostaticTank::new(&water, gravity(), geometry(2.0), 9.6, controls()).unwrap();
    let short =
        FiniteHydrostaticTank::new(&water, gravity(), geometry(0.805), 9.6, controls()).unwrap();
    assert!(matches!(
        prepared.load_finite_tank(&short, pose(0.9), &[0.0; 2], || false),
        Err(ArticulatedHydrostaticError::Physics(
            HydrostaticError::TankOverflow { .. }
        ))
    ));
    assert!(matches!(
        prepared.load_finite_tank(&tank, pose(0.1), &[0.0; 2], || false),
        Err(ArticulatedHydrostaticError::Physics(
            HydrostaticError::Invalid(_)
        ))
    ));
    assert!(matches!(
        prepared.load_finite_tank(&tank, pose(0.9), &[f64::NAN, 0.0], || false),
        Err(ArticulatedHydrostaticError::Kinematics(_))
    ));
    let mut polls = 0;
    let expected = prepared
        .load_finite_tank(&tank, pose(0.9), &[0.0; 2], || {
            polls += 1;
            false
        })
        .unwrap();
    for stop in [1, polls / 2, polls] {
        let mut seen = 0;
        assert!(matches!(
            prepared.load_finite_tank(&tank, pose(0.9), &[0.0; 2], || {
                seen += 1;
                seen == stop
            }),
            Err(ArticulatedHydrostaticError::Physics(
                HydrostaticError::Cancelled
            ))
        ));
        assert_eq!(
            prepared
                .load_finite_tank(&tank, pose(0.9), &[0.0; 2], || false)
                .unwrap(),
            expected,
        );
    }
}
