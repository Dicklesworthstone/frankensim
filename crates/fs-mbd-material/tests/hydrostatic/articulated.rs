//! G0/G3 physical articulated buoyancy through the original dynamics owners.
use super::*;
use fs_mbd::articulated::{
    BaseState, FreeFloatingBaseState, JointModel, forward_dynamics,
    free_floating_forward_dynamics, inverse_dynamics,
};
use fs_mbd_material::articulated::{MaterialLink, ResolvedArticulated, resolve_articulated};
use fs_mbd_material::hydrostatic::{
    ArticulatedHydrostaticError, DensityProfileHydrostaticMedium, HydrostaticDensitySample,
    HydrostaticLayer, LayeredHydrostaticMedium, LinearDensityHydrostaticMedium,
    MAX_DISPLACEMENT_CELLS, PreparedArticulatedHydrostatics, PrescribedHydrostaticReservoir,
};

fn shifted(x: f64) -> Se3 {
    Se3::from_parts(So3::identity(), Vec3::new(x, 0.0, 0.0)).unwrap()
}

fn hinge() -> JointModel {
    JointModel::revolute(Vec3::new(0.0, 1.0, 0.0), None).unwrap()
}

fn envelope() -> DisplacementMesh {
    let mesh = DisplacementMesh::cuboid(Vec3::new(1.0, 0.4, 0.4), || false).unwrap();
    let cells = mesh
        .cells()
        .iter()
        .map(|cell| (*cell).map(|v| v + Vec3::new(0.75, 0.0, 0.0)))
        .collect();
    DisplacementMesh::new(cells, 6, || false).unwrap()
}

fn mechanism(density: f64, floating: bool) -> ResolvedArticulated {
    let card = density_card("articulated-solid", density, density);
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
    resolve_articulated(
        &[
            MaterialLink {
                name: "root",
                parent: None,
                parent_from_child_zero: Se3::identity(),
                joint: if floating { JointModel::FIXED } else { hinge() },
                regions: &regions,
            },
            MaterialLink {
                name: "child",
                parent: Some(0),
                parent_from_child_zero: shifted(1.5),
                joint: hinge(),
                regions: &regions,
            },
        ],
        2,
        || false,
    )
    .unwrap()
}

#[test]
fn submerged_hollow_hinge_uses_displacement_and_link_origin_torque() {
    // Material COM x=0.6, sealed hull centroid x=0.75, pivot x=0.
    // This distinguishes all three origins and two independently owned volumes.
    let card = density_card("compact-ballast", 3000.0, 3000.0);
    let p = point(300.0);
    let regions = [RegionBinding {
        key: "ballast",
        geometry: VolumeMoments::cuboid(Vec3::new(0.5, 0.3, 0.3))
            .unwrap()
            .transformed(shifted(0.6))
            .unwrap(),
        material: &card,
        point: &p,
        selection: MaterialPropertySelection::SingleClaimOnly,
    }];
    let source = resolve_articulated(
        &[MaterialLink {
            name: "sealed-hinge",
            parent: None,
            parent_from_child_zero: Se3::identity(),
            joint: hinge(),
            regions: &regions,
        }],
        1,
        || false,
    )
    .unwrap();
    let meshes = [envelope()];
    let prepared =
        PreparedArticulatedHydrostatics::new(&source, &meshes, 6, 1e-12, || false).unwrap();
    let water = fluid(1000.0);
    let medium = HydrostaticMedium::new(&water, gravity(), origin()).unwrap();
    let angle: f64 = 0.37;
    let base_pose = pose(-2.0);
    let load = prepared.load(&medium, base_pose, &[angle], || false).unwrap();
    close(load.link_loads[0].displaced_mass_kg(), 160.0);
    close(
        load.external_body_wrenches[0].torque.y,
        -160.0 * 10.0 * 0.75 * angle.cos(),
    );
    let hold = inverse_dynamics(
        source.model(),
        BaseState::stationary(base_pose),
        &[angle],
        &[0.0],
        &[0.0],
        Vec3::new(0.0, 0.0, -10.0),
        &load.external_body_wrenches,
    )
    .unwrap();
    // Upward buoyancy and downward weight have different arms about the hinge.
    let expected = (160.0 * 0.75 - 135.0 * 0.6) * 10.0 * angle.cos();
    close(hold.generalized_force[0], expected);
    let acceleration = forward_dynamics(
        source.model(),
        BaseState::stationary(base_pose),
        &[angle],
        &[0.0],
        &[expected],
        Vec3::new(0.0, 0.0, -10.0),
        &load.external_body_wrenches,
    )
    .unwrap();
    assert!(acceleration.generalized_acceleration[0].abs() < 1e-11);
}

fn check_virtual_work<R: PrescribedHydrostaticReservoir>(
    prepared: &PreparedArticulatedHydrostatics<'_>,
    medium: &R,
) {
    let q = [0.18, -0.29];
    let base_pose = pose(0.08);
    let load = prepared.load(medium, base_pose, &q, || false).unwrap();
    for link in &load.link_loads {
        assert!(link.displaced().volume_m3() > 0.0);
        assert!(link.displaced().volume_m3() < 0.16);
    }
    // With gravity disabled only for this virtual-work projection, RNEA's
    // holding effort is the positive pressure-potential coordinate derivative.
    let holding = inverse_dynamics(
        prepared.source().model(),
        BaseState::stationary(base_pose),
        &q,
        &[0.0; 2],
        &[0.0; 2],
        origin(),
        &load.external_body_wrenches,
    )
    .unwrap();
    let h = 2e-6;
    for axis in 0..2 {
        let mut plus = q;
        let mut minus = q;
        plus[axis] += h;
        minus[axis] -= h;
        let energy_plus = prepared
            .load(medium, base_pose, &plus, || false)
            .unwrap()
            .pressure_potential_j;
        let energy_minus = prepared
            .load(medium, base_pose, &minus, || false)
            .unwrap()
            .pressure_potential_j;
        let derivative = (energy_plus - energy_minus) / (2.0 * h);
        assert!(
            (derivative - holding.generalized_force[axis]).abs() < 2e-5,
            "joint {axis}: dU/dq={derivative}, holding={}",
            holding.generalized_force[axis],
        );
    }
}

#[test]
fn all_prescribed_media_preserve_partial_immersion_joint_virtual_work() {
    let source = mechanism(1500.0, false);
    let meshes = [envelope(), envelope()];
    let prepared =
        PreparedArticulatedHydrostatics::new(&source, &meshes, 12, 1e-12, || false).unwrap();
    let water = fluid(1000.0);
    let brine = fluid(1200.0);
    let homogeneous = HydrostaticMedium::new(&water, gravity(), origin()).unwrap();
    check_virtual_work(&prepared, &homogeneous);
    let linear =
        LinearDensityHydrostaticMedium::new(&water, 500.0, gravity(), origin()).unwrap();
    check_virtual_work(&prepared, &linear);
    let layered = LayeredHydrostaticMedium::new(
        None,
        &[
            HydrostaticLayer { depth_m: 0.0, source: &water },
            HydrostaticLayer { depth_m: 0.15, source: &brine },
        ],
        gravity(),
        origin(),
        || false,
    )
    .unwrap();
    check_virtual_work(&prepared, &layered);
    let profile = DensityProfileHydrostaticMedium::new(
        &[
            HydrostaticDensitySample { depth_m: 0.0, source: &water },
            HydrostaticDensitySample { depth_m: 0.3, source: &brine },
        ],
        gravity(),
        origin(),
        || false,
    )
    .unwrap();
    check_virtual_work(&prepared, &profile);
}

#[test]
fn free_floating_linked_bodies_rise_float_or_sink_without_double_gravity() {
    let water = fluid(1000.0);
    let medium = HydrostaticMedium::new(&water, gravity(), origin()).unwrap();
    let meshes = [envelope(), envelope()];
    for (density, acceleration) in [(500.0, 10.0), (1000.0, 0.0), (2000.0, -5.0)] {
        let source = mechanism(density, true);
        let prepared =
            PreparedArticulatedHydrostatics::new(&source, &meshes, 12, 1e-12, || false).unwrap();
        let base_pose = pose(-2.0);
        let load = prepared.load(&medium, base_pose, &[0.31], || false).unwrap();
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
        assert!(base.linear.x.abs() < 1e-10 && base.linear.y.abs() < 1e-10);
        assert!(base.angular.dot(base.angular) < 1e-20);
        assert!(result.generalized_acceleration[0].abs() < 1e-10);
    }
}

#[test]
fn moving_the_world_frame_preserves_link_wrenches_and_pressure_potential() {
    let source = mechanism(1500.0, false);
    let meshes = [envelope(), envelope()];
    let prepared =
        PreparedArticulatedHydrostatics::new(&source, &meshes, 12, 1e-12, || false).unwrap();
    let water = fluid(1000.0);
    let medium = HydrostaticMedium::new(&water, gravity(), origin()).unwrap();
    let q = [0.18, -0.29];
    let base_pose = pose(0.08);
    let original = prepared.load(&medium, base_pose, &q, || false).unwrap();
    let rotation = So3::exp(So3Tangent::new(Vec3::new(0.3, -0.2, 0.4))).unwrap();
    let transform = Se3::from_parts(rotation, Vec3::new(7.0, -4.0, 3.0)).unwrap();
    let g = rotation.rotate(Vec3::new(0.0, 0.0, -10.0)).unwrap();
    let moved_medium = HydrostaticMedium::new(
        &water,
        Gravity::new(fs_mbd::Vec3::new(g.x, g.y, g.z)).unwrap(),
        transform.translation(),
    )
    .unwrap();
    let moved = prepared
        .load(&moved_medium, transform.compose(base_pose).unwrap(), &q, || false)
        .unwrap();
    close(original.pressure_potential_j, moved.pressure_potential_j);
    for (a, b) in original.external_body_wrenches.iter().zip(&moved.external_body_wrenches) {
        for (a, b) in a.to_array().into_iter().zip(b.to_array()) {
            assert!((a - b).abs() <= 1e-9 + 2e-12 * a.abs().max(b.abs()));
        }
    }
}

#[test]
fn admission_and_late_cancellation_do_not_publish_partial_link_loads() {
    let source = mechanism(1500.0, false);
    let meshes = [envelope(), envelope()];
    for budget in [0, 11, MAX_DISPLACEMENT_CELLS + 1] {
        assert!(
            PreparedArticulatedHydrostatics::new(&source, &meshes, budget, 1e-12, || false)
                .is_err()
        );
    }
    assert!(PreparedArticulatedHydrostatics::new(&source, &meshes[..1], 12, 1e-12, || false).is_err());
    assert!(PreparedArticulatedHydrostatics::new(&source, &meshes, 12, f64::NAN, || false).is_err());
    let prepared =
        PreparedArticulatedHydrostatics::new(&source, &meshes, 12, 1e-12, || false).unwrap();
    assert!(std::ptr::eq(prepared.source(), &source));
    assert_eq!(prepared.cell_count(), 12);
    let water = fluid(1000.0);
    let medium = HydrostaticMedium::new(&water, gravity(), origin()).unwrap();
    assert!(matches!(
        prepared.load(&medium, pose(0.08), &[0.18], || false),
        Err(ArticulatedHydrostaticError::Kinematics(_)),
    ));
    let mut polls = 0;
    let expected = prepared
        .load(&medium, pose(0.08), &[0.18, -0.29], || {
            polls += 1;
            false
        })
        .unwrap();
    for stop in [1, polls / 2, polls] {
        let mut count = 0;
        assert!(matches!(
            prepared.load(&medium, pose(0.08), &[0.18, -0.29], || {
                count += 1;
                count == stop
            }),
            Err(ArticulatedHydrostaticError::Physics(HydrostaticError::Cancelled)),
        ));
    }
    assert_eq!(
        prepared.load(&medium, pose(0.08), &[0.18, -0.29], || false).unwrap(),
        expected,
    );
}
