//! Native support motion closes an open crack before mixed contact carries load.

use super::*;
use fs_solid::brittle::solid::BrittleTetError;
use fs_solid::brittle::solid::fracture::contact::approach::{
    BrittleFacePathApproachControls, BrittleFacePathApproachError,
};
use fs_solid::brittle::solid::path::{BrittleTetPathControls, BrittleTetPathError};

const GAP: f64 = 1.0 / 1024.0;

fn approach_controls() -> BrittleFacePathApproachControls {
    BrittleFacePathApproachControls {
        path: BrittleTetPathControls {
            initial_fraction: 0.5,
            minimum_fraction: 1.0 / 64.0,
            maximum_fraction: 0.5,
            maximum_accepted_steps: 16,
            maximum_trials: 32,
            maximum_damage_increment: 0.1,
            initial_equilibrium_tolerance_n: 1e-10,
            work: work_controls(),
        },
        activation_tolerance_m: 1e-12,
        fraction_tolerance: 1e-8,
        maximum_prefix_trials: 64,
    }
}

fn with_open<R>(gap: f64, run: impl FnOnce(&BrittleTetSolid<'_>, &Cx<'_>) -> R) -> R {
    with_cx(false, |cx| {
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, budget(), 8, cx).unwrap();
        let materials = sources(&geometry);
        let failed = failed(&geometry, &materials, cx);
        let cut = topology(
            &geometry,
            &[InteriorSeam {
                id: 951,
                nodes: [0, 1, 2],
                side_a_element: 0,
            }],
            cx,
        );
        let target =
            TetLoadGeometry3::new(cut.nodes_m(), cut.tetrahedra(), budget(), 8, cx).unwrap();
        let split = failed
            .split_fracture(&cut, &target, fracture_controls(), cx)
            .unwrap();
        let open = split
            .solid
            .advance(
                &translated_upper(&split.solid, gap),
                &[],
                settings(),
                work_controls(),
                cx,
            )
            .unwrap();
        run(&open.solid, cx)
    })
}

fn translated_upper(solid: &BrittleTetSolid<'_>, translation: f64) -> Vec<(usize, f64)> {
    let mut displacement = solid.solution().displacement_m.clone();
    for &node in &solid.geometry().tetrahedra()[1] {
        displacement[node][2] += translation;
    }
    prescribed(&displacement)
}

fn minimum_gap(solid: &BrittleTetSolid<'_>) -> f64 {
    let face = solid.released_faces()[0];
    (0..3)
        .map(|vertex| {
            let a = solid.solution().displacement_m[face.paired_nodes[0][vertex]];
            let b = solid.solution().displacement_m[face.paired_nodes[1][vertex]];
            (0..3).fold(0.0, |sum, axis| {
                face.normal[axis].mul_add(b[axis] - a[axis], sum)
            })
        })
        .fold(f64::INFINITY, f64::min)
}

#[test]
fn open_support_path_reaches_native_loaded_recontact_without_projection() {
    with_open(GAP, |open, cx| {
        let original = open.solution().clone();
        let loss = open.fracture_dissipation_j();
        let endpoint = translated_upper(open, -3.0 * GAP);
        let numerical = approach_controls();
        let approached = open
            .advance_toward_released_face_contact(&endpoint, &[], settings(), numerical, cx)
            .unwrap();
        let contact = approached.contact.unwrap();
        assert_eq!(contact.face_id, 951);
        assert_eq!(contact.vertex, 0);
        assert!(contact.gap_m >= 0.0 && contact.gap_m <= numerical.activation_tolerance_m);
        assert!(approached.fraction <= 1.0 / 3.0);
        assert!(contact.bracket_end_fraction >= 1.0 / 3.0);
        assert!(contact.bracket_end_fraction - approached.fraction <= numerical.fraction_tolerance);
        assert!(approached.prefix_trials > 1);
        assert!(approached.prefix_trials <= numerical.maximum_prefix_trials);
        assert_eq!(approached.advance.records.len(), 2);
        assert_eq!(approached.advance.records[1].fraction, 1.0);
        for (dof, &(index, end)) in endpoint.iter().enumerate() {
            assert_eq!(dof, index);
            let start = original.displacement_m.as_flattened()[dof];
            let actual = approached
                .advance
                .solid
                .solution()
                .displacement_m
                .as_flattened()[dof];
            close(
                actual,
                (1.0 - approached.fraction) * start + approached.fraction * end,
                1e-17,
            );
        }
        close(
            approached.advance.solid.fracture_dissipation_j(),
            loss,
            1e-12,
        );
        close(approached.advance.energy.fracture_dissipation_j, 0.0, 1e-12);
        close(approached.advance.energy.material_work_j, 0.0, 1e-12);
        close(
            approached.advance.absolute_endpoint_work_defect_j,
            0.0,
            1e-12,
        );

        let (held, force) = supported_load(&approached.advance.solid);
        let loaded = approached
            .advance
            .solid
            .advance_released_face_equilibrium(&held, &force, &[0.5; 3], controls(), cx)
            .unwrap();
        for (i, row) in loaded.rows.iter().enumerate() {
            close(row.normal_force_n, (i + 1) as f64 / 8.0, 1e-10);
            close(row.gap_m, 0.0, controls().activation_tolerance_m);
        }
        close(loaded.native.solid.fracture_dissipation_j(), loss, 1e-12);
        close(
            loaded.absolute_contact_work_j,
            0.0,
            controls().contact_work_tolerance_j,
        );
        assert_eq!(open.solution(), &original);
        assert_eq!(open.fracture_dissipation_j(), loss);
    });
}

#[test]
fn a_fully_open_interval_keeps_the_native_path_and_its_receipts() {
    with_open(GAP, |open, cx| {
        let endpoint = translated_upper(open, -GAP / 4.0);
        let numerical = approach_controls();
        let direct = open
            .advance_path(&endpoint, &[], settings(), numerical.path, cx)
            .unwrap();
        let approached = open
            .advance_toward_released_face_contact(&endpoint, &[], settings(), numerical, cx)
            .unwrap();
        assert_eq!(approached.fraction, 1.0);
        assert_eq!(approached.prefix_trials, 1);
        assert_eq!(approached.contact, None);
        assert_eq!(approached.advance.solid.solution(), direct.solid.solution());
        assert_eq!(approached.advance.records, direct.records);
        assert_eq!(approached.advance.energy, direct.energy);
        assert_eq!(
            approached.advance.absolute_endpoint_work_defect_j,
            direct.absolute_endpoint_work_defect_j
        );
        close(
            minimum_gap(&approached.advance.solid),
            3.0 * GAP / 4.0,
            1e-16,
        );
    });
}

#[test]
fn approach_never_publishes_the_original_negative_roundoff_allowance() {
    with_open(GAP, |open, cx| {
        let mut numerical = approach_controls();
        numerical.activation_tolerance_m = 1e-14;
        numerical.fraction_tolerance = 1e-10;
        let endpoint = translated_upper(open, -GAP - fracture_controls().opening_tolerance_m / 4.0);
        // The original released-face policy deliberately admits this endpoint.
        // Event localization must retain an actual nonpenetrating native side.
        let direct = open
            .advance_path(&endpoint, &[], settings(), numerical.path, cx)
            .unwrap();
        assert!(minimum_gap(&direct.solid) < 0.0);
        let approached = open
            .advance_toward_released_face_contact(&endpoint, &[], settings(), numerical, cx)
            .unwrap();
        let contact = approached.contact.unwrap();
        assert!(approached.fraction < 1.0);
        assert!(minimum_gap(&approached.advance.solid) >= 0.0);
        assert!(contact.gap_m <= numerical.activation_tolerance_m);
        assert!(contact.bracket_end_fraction - approached.fraction <= numerical.fraction_tolerance);
        assert_eq!(
            open.released_faces()[0].opening_tolerance_m,
            fracture_controls().opening_tolerance_m
        );
    });
}

#[test]
fn native_refusals_and_contact_work_limits_keep_the_open_history() {
    with_open(GAP, |open, cx| {
        let original = open.solution().clone();
        let endpoint = translated_upper(open, -3.0 * GAP);
        let mut numerical = approach_controls();
        numerical.maximum_prefix_trials = 1;
        assert!(matches!(
            open.advance_toward_released_face_contact(&endpoint, &[], settings(), numerical, cx),
            Err(BrittleFacePathApproachError::WorkLimit {
                prefix_trials: 1,
                ..
            })
        ));

        let free_dof = 3 * open.released_faces()[0].paired_nodes[1][0] + 2;
        let held: Vec<_> = endpoint
            .iter()
            .copied()
            .filter(|(dof, _)| *dof != free_dof)
            .collect();
        assert!(matches!(
            open.advance_toward_released_face_contact(
                &held,
                &[(free_dof, [1.0, 1.0])],
                settings(),
                approach_controls(),
                cx
            ),
            Err(BrittleFacePathApproachError::Path(
                BrittleTetPathError::InitialEquilibrium { .. }
            ))
        ));
        let mut outside_domain = endpoint.clone();
        let apex = open.geometry().tetrahedra()[1]
            .iter()
            .copied()
            .find(|&node| open.geometry().nodes_m()[node][2] == 1.0)
            .unwrap();
        outside_domain
            .iter_mut()
            .find(|(dof, _)| *dof == 3 * apex + 2)
            .unwrap()
            .1 += 0.125;
        // This trial both closes the face and exceeds the actual material
        // strain domain. Preserve the native constitutive/solver refusal;
        // bisection must never relabel it as a contact event.
        let native_refusal = open
            .advance_path(
                &outside_domain,
                &[],
                settings(),
                approach_controls().path,
                cx,
            )
            .unwrap_err();
        assert!(matches!(
            &native_refusal,
            BrittleTetPathError::Endpoint(BrittleTetError::Material { .. })
                | BrittleTetPathError::Endpoint(BrittleTetError::Solve(
                    fs_solid::plastic3::SmallStrainTetError::InvalidInput { .. }
                ))
        ));
        assert_eq!(
            open.advance_toward_released_face_contact(
                &outside_domain,
                &[],
                settings(),
                approach_controls(),
                cx
            )
            .unwrap_err(),
            BrittleFacePathApproachError::Path(native_refusal)
        );
        with_cx(true, |cancelled| {
            assert!(matches!(
                open.advance_toward_released_face_contact(
                    &endpoint,
                    &[],
                    settings(),
                    approach_controls(),
                    cancelled
                ),
                Err(BrittleFacePathApproachError::Path(
                    BrittleTetPathError::Endpoint(BrittleTetError::Cancelled)
                ))
            ));
        });
        assert_eq!(open.solution(), &original);
        let retry = open
            .advance_toward_released_face_contact(
                &endpoint,
                &[],
                settings(),
                approach_controls(),
                cx,
            )
            .unwrap();
        assert!(retry.contact.is_some());
        assert_eq!(open.solution(), &original);
    });
    with_open(0.0, |touching, cx| {
        assert!(matches!(
            touching.advance_toward_released_face_contact(
                &translated_upper(touching, 0.0),
                &[],
                settings(),
                approach_controls(),
                cx
            ),
            Err(BrittleFacePathApproachError::AlreadyActive { .. })
        ));
    });
}
