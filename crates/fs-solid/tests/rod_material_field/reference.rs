//! Intrinsic shape tests use synthetic geometry/materials, not bow calibration.
use super::*;
use fs_ga::{So3, So3Tangent, Vec3};
use fs_math::det;
use fs_solid::rod::material_field::reference::RodReferenceGeometry;
use fs_solid::rod::material_field::supports::RodSupport;

fn rotation(angle: f64) -> So3 {
    So3::exp(So3Tangent::new(Vec3::new(0.0, 0.0, angle))).unwrap()
}

fn arc(count: usize, length: f64, curvature: f64) -> RodReferenceGeometry {
    let mut positions = Vec::new();
    let mut frames = Vec::new();
    for node in 0..=count {
        let angle = curvature * length * node as f64 / count as f64;
        positions.push([
            det::sin(angle) / curvature,
            (1.0 - det::cos(angle)) / curvature,
            0.0,
        ]);
        frames.push(rotation(angle));
    }
    RodReferenceGeometry::try_new(length, positions, frames).unwrap()
}

#[test]
fn g0_intrinsic_camber_and_twist_do_not_create_phantom_stress() {
    let curved = arc(4, 1.0, 0.7);
    assert!(
        curved
            .intrinsic_strains()
            .iter()
            .all(|(_, k)| (k[2] - 0.7).abs() < 1e-13)
    );
    let positions = curved.positions_m().to_vec();
    let rod = SpatialMaterialRod::try_from_reference(curved, vec![section(1.0); 4]).unwrap();
    assert_eq!(rod.positions_m(), positions);
    assert_eq!(rod.energy_j().unwrap(), 0.0);
    for segment in 0..4 {
        assert_eq!(rod.segment_strains(segment).unwrap(), ([0.0; 3], [0.0; 3]));
        assert_eq!(
            rod.segment_resultants(segment).unwrap(),
            ([0.0; 3], [0.0; 3])
        );
    }
    let frames = (0..=4)
        .map(|i| So3::exp(So3Tangent::new(Vec3::new(0.2 * i as f64, 0.0, 0.0))).unwrap())
        .collect();
    let twisted = RodReferenceGeometry::try_new(
        1.0,
        (0..=4).map(|i| [0.25 * i as f64, 0.0, 0.0]).collect(),
        frames,
    )
    .unwrap();
    let rod = SpatialMaterialRod::try_from_reference(twisted, vec![section(1.0); 4]).unwrap();
    assert_eq!(rod.energy_j().unwrap(), 0.0);
    assert!(
        rod.reference_geometry()
            .unwrap()
            .intrinsic_strains()
            .iter()
            .all(|(_, k)| (k[0] - 0.8).abs() < 1e-13)
    );
}

#[test]
fn g1_curved_material_rod_converges_to_circular_arc_under_end_moment() {
    let source = section(1.0);
    let moment = 0.1;
    let ei = source.elastic().bending_stiffness_n_m2();
    let curvature = 0.7 + moment / ei;
    let exact = [
        det::sin(curvature) / curvature,
        (1.0 - det::cos(curvature)) / curvature,
    ];
    let mut errors = Vec::new();
    for count in [2, 4, 8] {
        let mut rod = SpatialMaterialRod::try_from_reference(
            arc(count, 1.0, 0.7),
            vec![source.clone(); count],
        )
        .unwrap();
        rod.solve_static(
            &TipLoad {
                force: [0.0; 3],
                moment: [0.0, 0.0, moment],
            },
            2,
            2e-7,
            &CancelGate::new(),
        )
        .unwrap();
        near(rod.energy_j().unwrap(), 0.5 * moment * moment / ei, 2e-4);
        for segment in 0..count {
            near(rod.segment_resultants(segment).unwrap().1[2], moment, 2e-4);
        }
        let tip = rod.positions_m()[count];
        errors.push(det::hypot(tip[0] - exact[0], tip[1] - exact[1]));
    }
    for pair in errors.windows(2) {
        assert!(
            pair[0] / pair[1] > 3.5 && pair[0] / pair[1] < 4.5,
            "{errors:?}"
        );
    }
    assert!(errors[2] < 3e-5, "{errors:?}");
}

#[test]
fn g3_reference_geometry_preserves_material_sensitivity_and_gravity_datum() {
    let geometry = arc(2, 1.0, 0.7);
    let gate = CancelGate::new();
    let mut energies = Vec::new();
    for scale in [1.0, 2.0] {
        let rod = SpatialMaterialRod::try_from_reference(geometry.clone(), vec![section(scale); 2])
            .unwrap()
            .with_gravity([0.0, -9.81, 0.0], &gate)
            .unwrap();
        assert_eq!(rod.gravity_potential_j().unwrap(), 0.0);
        let mut tip = rod.positions_m()[2];
        tip[1] += 0.001;
        let changed = rod
            .with_supports(
                &[RodSupport {
                    node: 2,
                    position_m: tip.map(Some),
                    frame: Some(rod.frames()[2]),
                }],
                &gate,
            )
            .unwrap();
        energies.push(changed.energy_j().unwrap());
        near(
            changed.gravity_potential_j().unwrap(),
            rod.reference_mass().nodal_mass_kg()[2] * 9.81 * 0.001,
            1e-12,
        );
        assert_eq!(rod.energy_j().unwrap(), 0.0);
    }
    assert!(energies[0] > 0.0);
    near(energies[1], 2.0 * energies[0], 1e-13);
}

#[test]
fn g3_superposed_rigid_motion_preserves_relative_strain_and_energy() {
    let original = arc(2, 1.0, 0.7);
    let rotation = rotation(0.9);
    let transform = |p: [f64; 3]| {
        let p = rotation.rotate(Vec3::new(p[0], p[1], p[2])).unwrap();
        [p.x + 0.3, p.y - 0.2, p.z + 0.1]
    };
    let moved = RodReferenceGeometry::try_new(
        1.0,
        original
            .positions_m()
            .iter()
            .copied()
            .map(transform)
            .collect(),
        original
            .frames()
            .iter()
            .map(|&q| rotation.compose(q).unwrap())
            .collect(),
    )
    .unwrap();
    let a = SpatialMaterialRod::try_from_reference(original, vec![section(1.0); 2]).unwrap();
    let b = SpatialMaterialRod::try_from_reference(moved, vec![section(1.0); 2]).unwrap();
    let mut tip = a.positions_m()[2];
    tip[1] += 0.001;
    let gate = CancelGate::new();
    let a = a
        .with_supports(
            &[RodSupport {
                node: 2,
                position_m: tip.map(Some),
                frame: None,
            }],
            &gate,
        )
        .unwrap();
    let b = b
        .with_supports(
            &[RodSupport {
                node: 2,
                position_m: transform(tip).map(Some),
                frame: None,
            }],
            &gate,
        )
        .unwrap();
    near(a.energy_j().unwrap(), b.energy_j().unwrap(), 2e-11);
    for segment in 0..2 {
        let (ag, ak) = a.segment_strains(segment).unwrap();
        let (bg, bk) = b.segment_strains(segment).unwrap();
        for (x, y) in ag.iter().chain(&ak).zip(bg.iter().chain(&bk)) {
            assert!((x - y).abs() < 2e-14);
        }
    }
}

#[test]
fn g0_bad_reference_data_refuses_without_reinterpreting_a_specimen() {
    assert!(
        RodReferenceGeometry::try_new(
            1.0,
            vec![[0.0; 3], [2.0, 0.0, 0.0]],
            vec![So3::identity(); 2],
        )
        .is_err()
    );
    assert!(
        RodReferenceGeometry::try_new(1.0, vec![[0.0; 3]; 2], vec![So3::identity(); 2]).is_err()
    );
    assert!(
        RodReferenceGeometry::try_new(1.0, vec![[0.0; 3], [1.0, 0.0, 0.0]], vec![So3::identity()])
            .is_err()
    );
    assert!(
        RodReferenceGeometry::try_new(
            1.0,
            vec![[0.0; 3], [f64::NAN, 0.0, 0.0]],
            vec![So3::identity(); 2]
        )
        .is_err()
    );
    assert!(
        RodReferenceGeometry::try_new(
            1.0,
            vec![[0.0; 3], [1.0, 0.0, 0.0]],
            vec![So3::identity(), rotation(core::f64::consts::PI)]
        )
        .is_err()
    );
    assert!(SpatialMaterialRod::try_from_reference(arc(2, 1.0, 0.7), vec![section(1.0)]).is_err());
}
