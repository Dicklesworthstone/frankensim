//! G1/G3: original Cosserat mechanics on unequal physical material intervals.
use super::*;
use fs_ga::{So3, So3Tangent, Vec3};
use fs_solid::rod::material_field::motion::RodRootBoundary;
use fs_solid::rod::material_field::reference::RodReferenceGeometry;

fn arc(lengths: Vec<f64>, curvature: f64) -> RodReferenceGeometry {
    let mut station = 0.0;
    let mut positions = vec![[0.0; 3]];
    let mut frames = vec![So3::identity()];
    for &length in &lengths {
        station += length;
        let angle = curvature * station;
        positions.push([
            angle.sin() / curvature,
            (1.0 - angle.cos()) / curvature,
            0.0,
        ]);
        frames.push(So3::exp(So3Tangent::new(Vec3::new(0.0, 0.0, angle))).unwrap());
    }
    RodReferenceGeometry::try_nonuniform(lengths, positions, frames).unwrap()
}

#[test]
fn g1_each_original_interval_has_its_actual_axial_and_torsional_compliance() {
    let lengths = vec![0.125, 0.375, 0.5];
    let sections = vec![section(1.0), section(2.0), section(0.5)];
    let mut rod =
        SpatialMaterialRod::try_nonuniform_straight(lengths.clone(), sections.clone()).unwrap();
    let (force, moment) = (0.5, 0.01);
    let gate = CancelGate::new();
    rod.solve_static(
        &TipLoad {
            force: [force, 0.0, 0.0],
            moment: [moment, 0.0, 0.0],
        },
        2,
        1e-10,
        &gate,
    )
    .unwrap();
    let mut extension = 0.0;
    let mut twist = 0.0;
    let mut energy = 0.0;
    for (i, (&length, source)) in lengths.iter().zip(&sections).enumerate() {
        let ea = source.elastic().axial_stiffness_n();
        let gj = source.torsional_stiffness_n_m2();
        let (gamma, curvature) = rod.segment_strains(i).unwrap();
        near(gamma[0], force / ea, 2e-7);
        near(curvature[0], moment / gj, 2e-7);
        let (n, m) = rod.segment_resultants(i).unwrap();
        near(n[0], force, 2e-7);
        near(m[0], moment, 2e-7);
        extension += length * force / ea;
        twist += length * moment / gj;
        energy += 0.5 * length * (force * force / ea + moment * moment / gj);
    }
    near(rod.positions_m()[3][0] - 1.0, extension, 2e-7);
    near(
        rod.frames()[3]
            .body_minus(rod.frames()[0])
            .unwrap()
            .angular
            .x,
        twist,
        2e-7,
    );
    near(rod.energy_j().unwrap(), energy, 2e-7);
    assert_eq!(rod.reference_segment_lengths_m(), lengths);
    assert_eq!(
        rod.reference_material_coordinates_m(),
        &[0.0, 0.125, 0.5, 1.0]
    );
    let model = rod
        .clone()
        .into_motion_model(RodRootBoundary::Clamped, &gate)
        .unwrap();
    let response = model.evaluate(&model.initial_pose(), &gate).unwrap();
    near(response.force_n[0][0], force, 2e-7);
    near(response.force_n[3][0], -force, 2e-7);
    near(response.torque_body_nm[3][0], -moment, 2e-7);
    assert_eq!(response.potential_energy_j, rod.energy_j().unwrap());
}

#[test]
fn g1_curved_nonuniform_refinement_converges_to_the_same_loaded_circular_arc() {
    let source = section(1.0);
    let moment = 0.1;
    let ei = source.elastic().bending_stiffness_n_m2();
    let curvature = 0.7 + moment / ei;
    let exact = [
        curvature.sin() / curvature,
        (1.0 - curvature.cos()) / curvature,
    ];
    let mut errors = Vec::new();
    for refinement in [1, 2, 4] {
        let lengths: Vec<_> = [0.25, 0.75]
            .into_iter()
            .flat_map(|length| std::iter::repeat_n(length / refinement as f64, refinement))
            .collect();
        let geometry = arc(lengths.clone(), 0.7);
        let mut rod =
            SpatialMaterialRod::try_from_reference(geometry, vec![source.clone(); lengths.len()])
                .unwrap();
        assert_eq!(rod.energy_j().unwrap(), 0.0);
        for segment in 0..lengths.len() {
            assert_eq!(rod.segment_strains(segment).unwrap(), ([0.0; 3], [0.0; 3]));
        }
        rod.solve_static(
            &TipLoad {
                force: [0.0; 3],
                moment: [0.0, 0.0, moment],
            },
            2,
            1e-10,
            &CancelGate::new(),
        )
        .unwrap();
        near(rod.energy_j().unwrap(), 0.5 * moment * moment / ei, 2e-6);
        let tip = rod.positions_m().last().unwrap();
        errors.push((tip[0] - exact[0]).hypot(tip[1] - exact[1]));
    }
    for pair in errors.windows(2) {
        assert!((3.7..4.3).contains(&(pair[0] / pair[1])), "{errors:?}");
    }
    assert!(errors[2] < 5e-5, "{errors:?}");
}

#[test]
fn g1_reference_inventory_first_moments_and_all_rotary_axes_follow_each_interval() {
    let lengths = vec![0.125, 0.375, 0.5];
    let radii = [0.02_f64, 0.04, 0.03];
    let sections: Vec<_> = radii
        .iter()
        .map(|&radius| {
            let mut layers = section(1.0).elastic().layers().to_vec();
            layers[0].outer_radius_m = radius;
            ConcentricElasticSection::try_new(0.0, layers)
                .unwrap()
                .with_shear_modulus("shear_modulus", 0.8)
                .unwrap()
        })
        .collect();
    let gate = CancelGate::new();
    let rod = SpatialMaterialRod::try_nonuniform_straight(lengths.clone(), sections)
        .unwrap()
        .with_gravity([0.0, -9.0, 0.0], &gate)
        .unwrap();
    assert_eq!(rod.gravity_potential_j().unwrap(), 0.0);
    let stations = rod.reference_material_coordinates_m().to_vec();
    let model = rod.into_motion_model(RodRootBoundary::Free, &gate).unwrap();
    let mut nodal_mass = [0.0; 4];
    let mut nodal_diameter_inertia = [0.0; 4];
    let mut first_moment = 0.0;
    for i in 0..3 {
        let mass = 1000.0 * core::f64::consts::PI * radii[i].powi(2) * lengths[i];
        let diameter = 1000.0 * core::f64::consts::PI * radii[i].powi(4) / 4.0 * lengths[i];
        near(
            model.source().reference_mass().segment_mass_kg()[i],
            mass,
            2e-15,
        );
        for node in [i, i + 1] {
            nodal_mass[node] += 0.5 * mass;
            nodal_diameter_inertia[node] += 0.5 * diameter;
        }
        first_moment += 0.5 * mass * (stations[i] + stations[i + 1]);
    }
    let mut actual_moment = 0.0;
    for i in 0..4 {
        let inertia = model.inertia()[i];
        near(inertia.mass_kg(), nodal_mass[i], 2e-15);
        actual_moment += inertia.mass_kg() * stations[i];
        for (actual, expected) in inertia
            .principal_kg_m2()
            .into_iter()
            .zip([2.0, 1.0, 1.0].map(|factor| factor * nodal_diameter_inertia[i]))
        {
            near(actual, expected, 3e-15);
        }
    }
    near(actual_moment, first_moment, 2e-15);
    let mut pose = model.initial_pose();
    for point in &mut pose.positions_m {
        point[1] += 0.125;
    }
    let loads = model.body_loads(&pose, &gate).unwrap();
    near(
        loads.gravity_potential_j,
        9.0 * 0.125 * nodal_mass.iter().sum::<f64>(),
        2e-15,
    );
    for (force, mass) in loads.force_n.iter().zip(nodal_mass) {
        near(force[1], -9.0 * mass, 2e-15);
    }
}

#[test]
fn g3_nonuniform_finite_rotation_preserves_objective_energy_and_virtual_work() {
    let gate = CancelGate::new();
    let model = SpatialMaterialRod::try_from_reference(
        arc(vec![0.125, 0.375, 0.5], 0.7),
        vec![section(1.0), section(2.0), section(0.5)],
    )
    .unwrap()
    .into_motion_model(RodRootBoundary::Free, &gate)
    .unwrap();
    let mut pose = model.initial_pose();
    pose.positions_m[1][2] = 0.002;
    pose.frames[2] = pose.frames[2]
        .body_plus(So3Tangent::new(Vec3::new(0.02, 0.03, -0.01)))
        .unwrap();
    let before = model.evaluate(&pose, &gate).unwrap();
    let rotation = So3::exp(So3Tangent::new(Vec3::new(0.3, -0.2, 0.4))).unwrap();
    for (point, frame) in pose.positions_m.iter_mut().zip(&mut pose.frames) {
        let moved = rotation
            .rotate(Vec3::new(point[0], point[1], point[2]))
            .unwrap();
        *point = [moved.x + 0.2, moved.y - 0.3, moved.z + 0.1];
        *frame = rotation.compose(*frame).unwrap();
    }
    let after = model.evaluate(&pose, &gate).unwrap();
    near(after.potential_energy_j, before.potential_energy_j, 3e-12);
    for i in 0..4 {
        let force = before.force_n[i];
        let expected = rotation
            .rotate(Vec3::new(force[0], force[1], force[2]))
            .unwrap();
        for (actual, expected) in after.force_n[i]
            .into_iter()
            .zip([expected.x, expected.y, expected.z])
        {
            assert!((actual - expected).abs() < 1e-10);
        }
        for c in 0..3 {
            assert!((after.torque_body_nm[i][c] - before.torque_body_nm[i][c]).abs() < 1e-10);
        }
    }
}

#[test]
fn g3_uniform_intervals_keep_the_original_equilibrium_and_inventory_bits() {
    let sections = vec![section(1.0); 4];
    let mut original = SpatialMaterialRod::try_straight(1.0, sections.clone()).unwrap();
    let mut explicit =
        SpatialMaterialRod::try_nonuniform_straight(vec![0.25; 4], sections).unwrap();
    let load = TipLoad {
        force: [0.5, 0.001, 0.0],
        moment: [0.002, 0.0, 0.0005],
    };
    let gate = CancelGate::new();
    assert_eq!(
        original.solve_static(&load, 2, 1e-8, &gate).unwrap(),
        explicit.solve_static(&load, 2, 1e-8, &gate).unwrap()
    );
    assert_eq!(original.positions_m(), explicit.positions_m());
    assert_eq!(original.frames(), explicit.frames());
    assert_eq!(
        original.energy_j().unwrap().to_bits(),
        explicit.energy_j().unwrap().to_bits()
    );
    assert_eq!(original.reference_mass(), explicit.reference_mass());
}

#[test]
fn g0_material_station_roundoff_cannot_create_false_straight_prestress() {
    let material = section(1.0);
    for lengths in [vec![1.0, 3e-16], vec![1e16, 3.0]] {
        assert!(
            SpatialMaterialRod::try_nonuniform_straight(lengths, vec![material.clone(); 2],)
                .is_err()
        );
    }
    // A large grading ratio is admissible when every authored interval is
    // actually representable; neither its length nor inventory is rounded away.
    for lengths in [vec![1.0, 2.0_f64.powi(-50)], vec![1e16, 4.0]] {
        let rod =
            SpatialMaterialRod::try_nonuniform_straight(lengths.clone(), vec![material.clone(); 2])
                .unwrap();
        assert_eq!(rod.reference_segment_lengths_m(), lengths);
        assert_eq!(rod.energy_j().unwrap(), 0.0);
        for segment in 0..2 {
            assert_eq!(rod.segment_strains(segment).unwrap(), ([0.0; 3], [0.0; 3]));
        }
    }
}

#[test]
fn g0_g4_invalid_intervals_and_cancellation_publish_no_configuration() {
    let material = section(1.0);
    for lengths in [
        vec![],
        vec![0.0],
        vec![-1.0],
        vec![f64::NAN],
        vec![f64::INFINITY],
        vec![f64::MAX, f64::MAX],
        vec![1.0, f64::MIN_POSITIVE],
        vec![1.0; 65],
    ] {
        assert!(
            SpatialMaterialRod::try_nonuniform_straight(
                lengths.clone(),
                vec![material.clone(); lengths.len()]
            )
            .is_err()
        );
    }
    assert!(
        SpatialMaterialRod::try_nonuniform_straight(vec![0.25, 0.75], vec![material.clone()])
            .is_err()
    );
    assert!(
        RodReferenceGeometry::try_nonuniform(
            vec![0.125, 0.875],
            vec![[0.0; 3], [0.25, 0.0, 0.0], [1.0, 0.0, 0.0]],
            vec![So3::identity(); 3]
        )
        .is_err()
    );
    let mut rod =
        SpatialMaterialRod::try_nonuniform_straight(vec![0.25, 0.75], vec![material; 2]).unwrap();
    let before = rod.positions_m().to_vec();
    let gate = CancelGate::new();
    gate.request();
    assert_eq!(
        rod.solve_static(&TipLoad::default(), 1, 1e-9, &gate),
        Err(SpatialRodError::Cancelled)
    );
    assert_eq!(rod.positions_m(), before);
    assert!(rod.into_motion_model(RodRootBoundary::Free, &gate).is_err());
}
