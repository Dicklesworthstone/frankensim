//! G1/G3/G4 native principal tension under axial, bending and torsional loads.
use super::*;
use fs_ga::{So3Tangent, Vec3};
use fs_solid::rod::material_field::reference::RodReferenceGeometry;

fn area() -> f64 {
    core::f64::consts::PI * RADIUS.powi(2)
}
fn inertia() -> f64 {
    core::f64::consts::PI * RADIUS.powi(4) / 4.0
}
fn polar() -> f64 {
    2.0 * inertia()
}

// Independent closed-form eigenvalue of [[sigma,tau],[tau,0]]. The test values
// are moderate; the material owner's scaled implementation is not reused here.
fn principal(normal: f64, shear: f64) -> f64 {
    0.5 * (normal + (normal * normal + 4.0 * shear * shear).sqrt())
}

fn rotated_vector(rotation: So3, vector: [f64; 3]) -> [f64; 3] {
    let v = rotation
        .rotate(Vec3::new(vector[0], vector[1], vector[2]))
        .unwrap();
    [v.x, v.y, v.z]
}

#[test]
fn g1_native_torque_has_principal_tension_when_the_axial_check_is_zero() {
    let gate = CancelGate::new();
    let source = section(31, 300.0, &[(1e6, 100.0, 1e-5)]).0;
    let capacities = [capacity(&source)];
    for torque in [-0.02, 0.02] {
        let mut rod = SpatialMaterialRod::try_straight(1.0, vec![source.clone()]).unwrap();
        rod.solve_static_world(
            &WorldTipLoad {
                force: [0.0; 3],
                moment: [torque, 0.0, 0.0],
            },
            1,
            1e-9,
            &gate,
        )
        .unwrap();
        let axial = rod
            .assess_tensile_strength(&capacities, &[], 1, &gate)
            .unwrap();
        near(axial.utilization(), 0.0, 1e-9);
        let report = rod
            .assess_axial_bending_torsion_strength(&capacities, &[], 1, &gate)
            .unwrap();
        let sample = &report.segment_samples()[0];
        let expected = torque.abs() * RADIUS / polar();
        assert_eq!(sample.segment(), 0);
        assert_eq!(sample.capacity(), &capacities[0]);
        assert_eq!(sample.shear_section(), &rod.sections()[0]);
        near(sample.moment_body_n_m()[0], torque, 1e-8);
        near(
            sample.assessment().layers()[0].maximum_torsional_shear_pa,
            expected,
            2e-7,
        );
        near(
            sample.assessment().layers()[0].principal_tensile_stress_pa,
            expected,
            2e-7,
        );
        near(report.utilization(), expected / 100.0, 2e-7);
        assert!(report.reaches_limit());
        assert_eq!(report.governing_member(), RodTensileMember::Segment(0));
        assert_eq!(report.governing_layer(), 0);
    }
}

#[test]
fn g1_layered_native_torsion_uses_each_original_shear_modulus() {
    let gate = CancelGate::new();
    let source = section(39, 300.0, &[(1e6, 100.0, 1e-5), (2e6, 200.0, 1e-5)]).0;
    let capacities = [capacity(&source)];
    let mut rod = SpatialMaterialRod::try_straight(1.0, vec![source]).unwrap();
    let torque = 0.04;
    rod.solve_static_world(
        &WorldTipLoad {
            force: [0.0; 3],
            moment: [torque, 0.0, 0.0],
        },
        1,
        1e-9,
        &gate,
    )
    .unwrap();
    let report = rod
        .assess_axial_bending_torsion_strength(&capacities, &[], 2, &gate)
        .unwrap();
    // Equal areas give inner/outer polar moments J/4 and 3J/4. Common twist
    // shares torque through G_i J_i, not through one averaged shear stress.
    let twist = torque / ((0.4e6 / 4.0 + 0.8e6 * 3.0 / 4.0) * polar());
    for (index, (shear, limit)) in [(0.4e6, 100.0), (0.8e6, 200.0)].into_iter().enumerate() {
        let expected = shear * twist * RADIUS * ((index + 1) as f64 / 2.0).sqrt();
        let layer = &report.segment_samples()[0].assessment().layers()[index];
        near(layer.maximum_torsional_shear_pa, expected, 2e-7);
        near(layer.principal_tensile_stress_pa, expected, 2e-7);
        near(layer.tensile_utilization, expected / limit, 2e-7);
    }
    assert_eq!(report.governing_layer(), 1);
    assert!(report.reaches_limit());
}

#[test]
fn g1_g3_axial_torsion_retains_compression_and_is_objective_under_world_rotation() {
    let gate = CancelGate::new();
    let source = section(32, 300.0, &[(1e6, 100.0, 1e-5)]).0;
    let capacities = [capacity(&source)];
    let rotations = [
        So3::identity(),
        So3::exp(So3Tangent::new(Vec3::new(0.3, -0.4, 0.2))).unwrap(),
    ];
    for axial_force in [-0.5, 0.5] {
        let expected = principal(axial_force / area(), 0.02 * RADIUS / polar());
        for rotation in rotations {
            let reference = RodReferenceGeometry::try_new(
                1.0,
                vec![[0.0; 3], rotated_vector(rotation, [1.0, 0.0, 0.0])],
                vec![rotation; 2],
            )
            .unwrap();
            let mut rod =
                SpatialMaterialRod::try_from_reference(reference, vec![source.clone()]).unwrap();
            rod.solve_static_world(
                &WorldTipLoad {
                    force: rotated_vector(rotation, [axial_force, 0.0, 0.0]),
                    moment: rotated_vector(rotation, [0.02, 0.0, 0.0]),
                },
                1,
                1e-9,
                &gate,
            )
            .unwrap();
            let report = rod
                .assess_axial_bending_torsion_strength(&capacities, &[], 1, &gate)
                .unwrap();
            let layer = &report.segment_samples()[0].assessment().layers()[0];
            near(layer.minimum_axial_stress_pa, axial_force / area(), 2e-6);
            near(layer.maximum_axial_stress_pa, axial_force / area(), 2e-6);
            near(layer.principal_tensile_stress_pa, expected, 2e-6);
            near(report.utilization(), expected / 100.0, 2e-6);
            assert!(layer.principal_tensile_stress_pa > 0.0);
            assert_eq!(report.reaches_limit(), axial_force > 0.0);
        }
    }
}

#[test]
fn g3_native_bending_and_torsion_combine_at_the_same_section_surface() {
    let gate = CancelGate::new();
    let source = section(33, 300.0, &[(1e6, 1000.0, 1e-5)]).0;
    let capacities = [capacity(&source)];
    let mut utilizations = Vec::new();
    for bending in [[0.0, 0.1], [0.06, 0.08], [-0.08, 0.06]] {
        let mut rod = SpatialMaterialRod::try_straight(1.0, vec![source.clone()]).unwrap();
        rod.solve_static(
            &TipLoad {
                force: [0.0; 3],
                moment: [0.02, bending[0], bending[1]],
            },
            1,
            1e-9,
            &gate,
        )
        .unwrap();
        let report = rod
            .assess_axial_bending_torsion_strength(&capacities, &[], 1, &gate)
            .unwrap();
        let sample = &report.segment_samples()[0];
        let moment = sample.moment_body_n_m();
        let normal =
            sample.force_body_n()[0] / area() + moment[1].hypot(moment[2]) * RADIUS / inertia();
        let shear = moment[0].abs() * RADIUS / polar();
        let layer = &sample.assessment().layers()[0];
        near(layer.maximum_axial_stress_pa, normal, 1e-12);
        near(layer.maximum_torsional_shear_pa, shear, 1e-12);
        near(
            layer.principal_tensile_stress_pa,
            principal(normal, shear),
            1e-12,
        );
        assert!(layer.principal_tensile_stress_pa > layer.maximum_axial_stress_pa);
        utilizations.push(report.utilization());
    }
    for &value in &utilizations[1..] {
        near(value, utilizations[0], 3e-6);
    }
}

#[test]
fn g1_thermal_mismatch_and_native_torque_share_original_stresses_and_strain_domain() {
    let gate = CancelGate::new();
    let (source, thermal) = section(34, 400.0, &[(1e6, 1000.0, 1e-5), (1e6, 1000.0, 3e-5)]);
    let capacities = [capacity(&source)];
    let mut rod = SpatialMaterialRod::try_straight(1.0, vec![source])
        .unwrap()
        .with_thermal_sections(vec![thermal], 0.004, &gate)
        .unwrap();
    rod.solve_static_world(
        &WorldTipLoad {
            force: [0.0; 3],
            moment: [0.1, 0.0, 0.0],
        },
        1,
        1e-9,
        &gate,
    )
    .unwrap();
    {
        let report = rod
            .assess_axial_bending_torsion_strength(&capacities, &[], 2, &gate)
            .unwrap();
        let sample = &report.segment_samples()[0];
        near(sample.force_body_n()[0], 0.0, 1e-8);
        near(sample.total_axial_strain(), 0.002, 1e-8);
        for (index, normal) in [1000.0, -1000.0].into_iter().enumerate() {
            let radius = RADIUS * ((index + 1) as f64 / 2.0).sqrt();
            let shear = 0.1 * radius / polar();
            let layer = &sample.assessment().layers()[index];
            near(layer.maximum_axial_stress_pa, normal, 2e-7);
            near(layer.maximum_torsional_shear_pa, shear, 2e-7);
            near(
                layer.principal_tensile_stress_pa,
                principal(normal, shear),
                2e-7,
            );
            near(
                layer.tensile_utilization,
                principal(normal, shear) / 500.0,
                2e-7,
            );
        }
        assert_eq!(report.governing_layer(), 0);
        assert!(report.utilization() > 2.0);
    }
    // Axial/thermal strains still satisfy the source's bound, but native outer
    // torsional gamma now exceeds it. The combined assessment must refuse.
    rod.solve_static_world(
        &WorldTipLoad {
            force: [0.0; 3],
            moment: [0.4, 0.0, 0.0],
        },
        1,
        1e-9,
        &gate,
    )
    .unwrap();
    assert!(rod.thermal_axial_response(0).unwrap().is_some());
    let before = rod.positions_m().to_vec();
    let frames = rod.frames().to_vec();
    let energy = rod.energy_j().unwrap();
    assert!(matches!(
        rod.assess_axial_bending_torsion_strength(&capacities, &[], 2, &gate),
        Err(SpatialRodError::Material(_))
    ));
    assert_eq!(rod.positions_m(), before);
    assert_eq!(rod.frames(), frames);
    assert_eq!(rod.energy_j().unwrap(), energy);
}

#[test]
fn g3_combined_report_marks_transverse_shear_without_inventing_its_stress() {
    let gate = CancelGate::new();
    let source = section(35, 300.0, &[(1e6, 100.0, 1e-5)]).0;
    let capacities = [capacity(&source)];
    let rod = SpatialMaterialRod::try_straight(1.0, vec![source])
        .unwrap()
        .with_supports(
            &[RodSupport {
                node: 1,
                position_m: [Some(1.0), Some(0.001), Some(0.0)],
                frame: Some(So3::identity()),
            }],
            &gate,
        )
        .unwrap();
    let report = rod
        .assess_axial_bending_torsion_strength(&capacities, &[], 1, &gate)
        .unwrap();
    let sample = &report.segment_samples()[0];
    let force = rod.segment_resultants(0).unwrap().0;
    assert_eq!(sample.unassessed_transverse_force_n(), [force[1], force[2]]);
    assert!(sample.unassessed_transverse_force_n()[0] > 1.0);
    assert!(sample.has_unassessed_transverse_shear());
    assert!(report.has_unassessed_transverse_shear());
    // This is explicitly a partial criterion, not a "safe material" verdict.
    assert_eq!(
        sample.assessment().layers()[0].maximum_torsional_shear_pa,
        0.0
    );
    assert_eq!(report.utilization(), 0.0);
}

#[test]
fn g3_g4_strands_keep_their_axial_law_and_combined_refusals_are_transactional() {
    let gate = CancelGate::new();
    let source = section(36, 300.0, &[(1e6, 1e9, 1e-5)]).0;
    let strand = section(37, 300.0, &[(1e6, 1e4, 1e-5)]).0;
    let capacities = [capacity(&source)];
    let foreign = [capacity(&section(38, 300.0, &[(1e6, 1e9, 1e-5)]).0)];
    let link_capacities = [capacity(&strand), capacity(&strand)];
    let bare = SpatialMaterialRod::try_straight(1.0, vec![source])
        .unwrap()
        .with_supports(
            &[RodSupport {
                node: 1,
                position_m: [Some(1.0), Some(0.0), Some(0.0)],
                frame: Some(So3::identity()),
            }],
            &gate,
        )
        .unwrap();
    for count in [1, 7] {
        let links = [0.99, 1.01].map(|natural_length| {
            MaterialTensionLink::try_new(
                strand.elastic().clone(),
                [
                    TensionLinkAttachment {
                        node: 0,
                        offset_body_m: [0.0; 3],
                    },
                    TensionLinkAttachment {
                        node: 1,
                        offset_body_m: [0.0; 3],
                    },
                ],
                natural_length,
                count,
                0.05,
            )
            .unwrap()
        });
        let rod = bare.with_tension_links(links.into(), &gate).unwrap();
        let before = rod.positions_m().to_vec();
        let energy = rod.energy_j().unwrap();
        let report = rod
            .assess_axial_bending_torsion_strength(&capacities, &link_capacities, 3, &gate)
            .unwrap();
        let axial = rod
            .assess_tensile_strength(&capacities, &link_capacities, 3, &gate)
            .unwrap();
        for index in 0..2 {
            assert_eq!(&report.link_samples()[index], &axial.samples()[index + 1]);
        }
        let expected = 1e6 * ((1.0 - 0.99) / 0.99) / 1e4;
        near(report.utilization(), expected, 1e-14);
        assert_eq!(report.governing_member(), RodTensileMember::TensionLink(0));
        assert_eq!(report.link_samples()[1].assessment().utilization(), 0.0);
        assert!(!report.has_unassessed_transverse_shear());
        assert!(
            rod.assess_axial_bending_torsion_strength(&capacities, &link_capacities, 2, &gate)
                .is_err()
        );
        assert!(
            rod.assess_axial_bending_torsion_strength(&capacities, &[], 3, &gate)
                .is_err()
        );
        assert!(
            rod.assess_axial_bending_torsion_strength(&foreign, &link_capacities, 3, &gate)
                .is_err()
        );
        let cancelled = CancelGate::new();
        cancelled.request();
        assert_eq!(
            rod.assess_axial_bending_torsion_strength(&capacities, &link_capacities, 3, &cancelled),
            Err(SpatialRodError::Cancelled)
        );
        assert_eq!(rod.positions_m(), before);
        assert_eq!(rod.energy_j().unwrap(), energy);
    }
}
