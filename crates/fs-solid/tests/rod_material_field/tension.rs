//! G1/G3 manufactured tensile-link checks; no instrument calibration.
use super::*;
use fs_ga::{So3, So3Tangent, Vec3};
use fs_math::det;
use fs_solid::rod::material_field::reference::RodReferenceGeometry;
use fs_solid::rod::material_field::supports::RodSupport;
use fs_solid::rod::material_field::tension::{MaterialTensionLink, TensionLinkAttachment};

fn link(scale: f64, natural: f64, offset: [f64; 3], end: usize) -> MaterialTensionLink {
    let material = section(scale).elastic().layers()[0].material.clone();
    let section = ConcentricElasticSection::try_new(
        0.0,
        vec![ConcentricMaterialLayer {
            outer_radius_m: 0.005,
            material,
        }],
    )
    .unwrap();
    MaterialTensionLink::try_new(
        section,
        [
            TensionLinkAttachment {
                node: 0,
                offset_body_m: offset,
            },
            TensionLinkAttachment {
                node: end,
                offset_body_m: offset,
            },
        ],
        natural,
        3,
        0.1,
    )
    .unwrap()
}
fn rotation(angle: f64) -> So3 {
    So3::exp(So3Tangent::new(Vec3::new(0.0, 0.0, angle))).unwrap()
}

#[test]
fn g1_shared_axial_equilibrium_obeys_actual_rod_and_link_compliances() {
    let gate = CancelGate::new();
    let source = section(1.0);
    for scale in [0.5, 2.0] {
        let link = link(scale, 0.99, [0.0; 3], 2);
        let kr = source.elastic().axial_stiffness_n();
        let kl = link.section().axial_stiffness_n() * f64::from(link.parallel_strands())
            / link.natural_length_m();
        let expected_length = (kr + kl * 0.99) / (kr + kl);
        let mut rod = SpatialMaterialRod::try_straight(1.0, vec![source.clone(); 2])
            .unwrap()
            .with_tension_links(vec![link], &gate)
            .unwrap()
            .with_supports(
                &(1..=2)
                    .map(|node| RodSupport {
                        node,
                        position_m: [None, Some(0.0), Some(0.0)],
                        frame: Some(So3::identity()),
                    })
                    .collect::<Vec<_>>(),
                &gate,
            )
            .unwrap();
        rod.solve_static(&TipLoad::default(), 1, 1e-7, &gate)
            .unwrap();
        let state = &rod.tension_link_states().unwrap()[0];
        near(rod.positions_m()[2][0], expected_length, 1e-8);
        near(state.tension_n, kl * (expected_length - 0.99), 3e-5);
        near(
            rod.segment_resultants(0).unwrap().0[0],
            -state.tension_n,
            3e-5,
        );
        near(
            rod.energy_j().unwrap(),
            0.5 * kr * (expected_length - 1.0).powi(2)
                + 0.5 * kl * (expected_length - 0.99).powi(2),
            3e-5,
        );
    }
}

#[test]
fn g1_offset_link_energy_supplies_the_same_support_force_and_body_moment() {
    let gate = CancelGate::new();
    let bare = SpatialMaterialRod::try_straight(1.0, vec![section(1.0)])
        .unwrap()
        .with_supports(
            &[RodSupport {
                node: 1,
                position_m: [Some(1.0), Some(0.002), Some(0.0)],
                frame: Some(rotation(0.1)),
            }],
            &gate,
        )
        .unwrap();
    let linked = bare
        .with_tension_links(vec![link(1.0, 0.98, [0.0, 0.02, 0.0], 1)], &gate)
        .unwrap();
    let state = &linked.tension_link_states().unwrap()[0];
    near(
        linked.energy_j().unwrap() - bare.energy_j().unwrap(),
        state.stored_energy_j,
        1e-12,
    );
    for c in 0..3 {
        assert_eq!(state.force_n[0][c], -state.force_n[1][c]);
    }
    let a = bare.support_reactions(&TipLoad::default(), &gate).unwrap();
    let b = linked
        .support_reactions(&TipLoad::default(), &gate)
        .unwrap();
    for c in 0..3 {
        let reaction = b.nodes[0].force_n[c].unwrap() - a.nodes[0].force_n[c].unwrap();
        assert!((reaction + state.force_n[1][c]).abs() < 2e-6);
        let reaction =
            b.nodes[0].moment_body_nm.unwrap()[c] - a.nodes[0].moment_body_nm.unwrap()[c];
        assert!((reaction + state.moment_body_nm[1][c]).abs() < 2e-6);
    }
    let root_a = bare.root_support_reaction(&gate).unwrap();
    let root_b = linked.root_support_reaction(&gate).unwrap();
    assert!(
        (root_b.moment_body_nm.unwrap()[2] - root_a.moment_body_nm.unwrap()[2]
            + state.moment_body_nm[0][2])
            .abs()
            < 2e-6
    );
}

#[test]
fn g3_intrinsic_camber_and_stick_material_determine_tension_together() {
    let gate = CancelGate::new();
    let curvature = 0.2;
    let positions: Vec<_> = (0..=2)
        .map(|i| {
            let a = curvature * i as f64 / 2.0;
            [
                det::sin(a) / curvature,
                (1.0 - det::cos(a)) / curvature,
                0.0,
            ]
        })
        .collect();
    let geometry = RodReferenceGeometry::try_new(
        1.0,
        positions.clone(),
        (0..=2)
            .map(|i| rotation(curvature * i as f64 / 2.0))
            .collect(),
    )
    .unwrap();
    let mut results = Vec::new();
    for stiffness in [0.5, 2.0] {
        let mut rod =
            SpatialMaterialRod::try_from_reference(geometry.clone(), vec![section(stiffness); 2])
                .unwrap()
                .with_tension_links(vec![link(1.0, 0.998, [0.0, -0.02, 0.0], 2)], &gate)
                .unwrap();
        let initial = rod.tension_link_states().unwrap()[0].tension_n;
        rod.solve_static(&TipLoad::default(), 1, 2e-7, &gate)
            .unwrap();
        let tension = rod.tension_link_states().unwrap()[0].tension_n;
        assert!(tension > 0.0 && tension < initial);
        assert_ne!(rod.positions_m(), positions);
        assert_eq!(
            rod.reference_geometry().unwrap().positions_m(),
            geometry.positions_m()
        );
        assert_eq!(
            rod.reference_geometry().unwrap().frames(),
            geometry.frames()
        );
        results.push(tension);
    }
    assert!(
        results[0] < results[1],
        "more compliant stick releases more preload: {results:?}"
    );
}

#[test]
fn g1_link_reference_mass_and_rotating_gravity_offset_are_counted_once() {
    let gate = CancelGate::new();
    let gravity = [0.0, -9.81, 0.0];
    let offset = [0.03, 0.02, 0.0];
    let link = link(1.0, 1.02, offset, 1);
    let mass = link.reference_mass_kg();
    let rod = SpatialMaterialRod::try_straight(1.0, vec![section(1.0)]).unwrap();
    let linked = rod
        .with_gravity(gravity, &gate)
        .unwrap()
        .with_tension_links(vec![link.clone()], &gate)
        .unwrap();
    near(
        linked.total_reference_mass_kg(),
        rod.reference_mass().total_mass_kg() + mass,
        1e-14,
    );
    near(
        linked
            .gravity_nodal_forces_n()
            .unwrap()
            .iter()
            .map(|f| f[1])
            .sum(),
        -9.81 * linked.total_reference_mass_kg(),
        1e-14,
    );
    assert_eq!(linked.gravity_potential_j().unwrap(), 0.0);
    let frame = rotation(0.3);
    let rotated = linked
        .with_supports(
            &[RodSupport {
                node: 1,
                position_m: [Some(1.0), Some(0.0), Some(0.0)],
                frame: Some(frame),
            }],
            &gate,
        )
        .unwrap();
    assert_eq!(rotated.tension_link_states().unwrap()[0].tension_n, 0.0);
    let rotated_offset = frame
        .rotate(Vec3::new(offset[0], offset[1], offset[2]))
        .unwrap();
    near(
        rotated.gravity_potential_j().unwrap(),
        0.5 * mass * 9.81 * (rotated_offset.y - offset[1]),
        1e-13,
    );
    let expected = -0.5 * mass * 9.81 * (offset[0] * det::cos(0.3) - offset[1] * det::sin(0.3));
    near(
        rotated.gravity_nodal_moments_body_nm().unwrap().unwrap()[1][2],
        expected,
        1e-13,
    );
    let replaced = linked.with_tension_links(vec![link], &gate).unwrap();
    assert_eq!(
        replaced.gravity_nodal_forces_n(),
        linked.gravity_nodal_forces_n()
    );
    let removed = linked.with_tension_links(Vec::new(), &gate).unwrap();
    assert_eq!(
        removed.gravity_nodal_forces_n(),
        rod.with_gravity(gravity, &gate)
            .unwrap()
            .gravity_nodal_forces_n()
    );
}

#[test]
fn g0_slack_source_domain_and_cancelled_replacement_keep_the_old_specimen() {
    let gate = CancelGate::new();
    let bare = SpatialMaterialRod::try_straight(1.0, vec![section(1.0)]).unwrap();
    let slack = bare
        .with_tension_links(vec![link(1.0, 1.01, [0.0; 3], 1)], &gate)
        .unwrap();
    let state = &slack.tension_link_states().unwrap()[0];
    assert_eq!(state.tension_n, 0.0);
    assert_eq!(state.stored_energy_j, 0.0);
    assert_eq!(slack.energy_j().unwrap(), bare.energy_j().unwrap());
    assert!(
        slack
            .with_tension_links(vec![link(1.0, 0.5, [0.0; 3], 1)], &gate)
            .is_err()
    );
    assert!(
        slack
            .with_tension_links(vec![link(1.0, 1.0, [0.0; 3], 2)], &gate)
            .is_err()
    );
    assert!(
        slack
            .with_tension_links(vec![link(1.0, 1.0, [0.0; 3], 1); 65], &gate)
            .is_err()
    );
    gate.request();
    assert!(matches!(
        slack.with_tension_links(Vec::new(), &gate),
        Err(SpatialRodError::Cancelled)
    ));
    assert_eq!(slack.tension_links().len(), 1);
    assert_eq!(slack.positions_m(), bare.positions_m());
}
