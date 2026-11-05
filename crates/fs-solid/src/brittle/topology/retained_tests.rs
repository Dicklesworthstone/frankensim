//! G0/G1/G3: additional cuts retain physical cohesive points and their owners.
use super::*;
use crate::cohesive::p1::{MaterialFacePair, P1Interfaces, TraceBudget};
use fs_material::cohesive::source::ResolvedBilinearCohesive;

const NODES: [[f64; 3]; 5] = [
    [0.0, 0.0, -1.0],
    [0.0, 0.0, 1.0],
    [1.0, 0.0, 0.0],
    [-0.5, 0.9, 0.0],
    [-0.5, -0.9, 0.0],
];
const TETS: [[usize; 4]; 3] = [[0, 1, 2, 3], [0, 1, 3, 4], [0, 1, 4, 2]];
const VOLUMES: [f64; 3] = [0.3; 3];
const FIRST: InteriorSeam = InteriorSeam {
    id: 19,
    nodes: [3, 0, 1],
    side_a_element: 1,
};
const MORE: [InteriorSeam; 2] = [
    // Original vertex 4 keeps its ID while 0/1 are duplicated on A. Sorting
    // the new IDs would therefore rotate this material frame and point order.
    InteriorSeam {
        id: 27,
        nodes: [0, 1, 4],
        side_a_element: 1,
    },
    InteriorSeam {
        id: 38,
        nodes: [2, 1, 0],
        side_a_element: 0,
    },
];
fn budget() -> FractureTopologyBudget {
    FractureTopologyBudget {
        maximum_output_nodes: 32,
        maximum_face_entries: 12,
        maximum_seams: 3,
    }
}
fn trace_budget() -> TraceBudget {
    TraceBudget {
        maximum_nodes: 32,
        maximum_tetrahedra: 3,
        maximum_interfaces: 3,
        minimum_face_quality: 1e-10,
    }
}
fn first() -> TetFractureTopology<'static> {
    TetFractureTopology::prepare(&NODES, &TETS, &VOLUMES, &[FIRST], budget(), &mut || false)
        .unwrap()
}
fn law() -> BilinearCohesive {
    BilinearCohesive::new(1e9, 2e5, 100., 2e9, 0.01).unwrap()
}
fn trace_at<'a>(
    nodes: &[[f64; 3]],
    tetrahedra: &[[usize; 4]],
    pairs: &[FacePair<'a>],
    ordered: bool,
) -> P1Interfaces<'a> {
    let free: Vec<_> = (0..3 * nodes.len()).collect();
    if ordered {
        P1Interfaces::new_ordered_controlled(
            nodes,
            tetrahedra,
            &free,
            pairs,
            trace_budget(),
            &mut || false,
        )
    } else {
        P1Interfaces::new_controlled(nodes, tetrahedra, &free, pairs, trace_budget(), &mut || {
            false
        })
    }
    .unwrap()
}
fn trace<'a>(topology: &TetFractureTopology<'_>, law: &'a BilinearCohesive) -> P1Interfaces<'a> {
    let pairs: Vec<_> = topology
        .seams()
        .iter()
        .map(|seam| seam.face_pair(law))
        .collect();
    trace_at(topology.nodes_m(), topology.tetrahedra(), &pairs, true)
}
fn opening(topology: &TetFractureTopology<'_>, map: &P1Interfaces<'_>, amount: f64) -> Vec<f64> {
    let mut u = vec![0.0; map.full_dimension()];
    let face = topology.seams()[0];
    for (a, b) in face.side_a.into_iter().zip(face.side_b) {
        if a != b {
            for c in 0..3 {
                u[3 * b + c] = amount * map.points()[0].axes[0][c];
            }
        }
    }
    u
}

#[test]
fn successive_cuts_retain_evolved_point_histories_and_their_physical_jumps() {
    let law = law();
    let old = first();
    let old_map = trace(&old, &law);
    let mut history = vec![law.bonded_point(); old_map.points().len()];
    let mut u = Vec::new();
    for amount in [0.0003, 0.0033, 0.0017] {
        u = opening(&old, &old_map, amount);
        for (point, jump) in history
            .iter_mut()
            .zip(old_map.jumps(&u, &mut || false).unwrap())
        {
            *point = point.trial(jump, 1e-7).unwrap().point;
        }
    }
    assert!(history.iter().all(|point| point.damage() > 0.0));
    assert!(history.iter().any(|point| point.failed()));
    let saved = history.clone();
    let first_cut = old
        .refine_controlled(&MORE[..1], budget(), &mut || false)
        .unwrap();
    let middle_map = trace(first_cut.topology(), &law);
    first_cut
        .require_retained_interfaces_controlled(&old_map, &middle_map, &mut || false)
        .unwrap();
    let middle = first_cut
        .transfer_motion_controlled(&u, &vec![0.0; u.len()], &[2.0, 3.0, 5.0], &mut || false)
        .unwrap();
    let second_cut = first_cut
        .topology()
        .refine_controlled(&MORE[1..], budget(), &mut || false)
        .unwrap();
    let next_map = trace(second_cut.topology(), &law);
    second_cut
        .require_retained_interfaces_controlled(&middle_map, &next_map, &mut || false)
        .unwrap();
    let next = second_cut
        .transfer_motion_controlled(
            &middle.displacement_m,
            &middle.velocity_m_s,
            &[2.0, 3.0, 5.0],
            &mut || false,
        )
        .unwrap();
    let old_jumps = old_map.jumps(&u, &mut || false).unwrap();
    let next_jumps = next_map.jumps(&next.displacement_m, &mut || false).unwrap();
    assert_eq!(&next_jumps[..old_jumps.len()], old_jumps);
    for (point, jump) in history.iter().zip(next_jumps) {
        let retained = point.trial(jump, 1e-7).unwrap();
        assert_eq!(retained.point, *point);
        assert_eq!(retained.energy.work, 0.0);
        assert_eq!(retained.energy.dissipated, 0.0);
        assert!(core::ptr::eq(point.law(), &law));
    }
    assert_eq!(history, saved);
}

#[test]
fn reordered_roles_ids_incomplete_maps_and_unrelated_geometry_refuse() {
    let law = law();
    let old = first();
    let old_map = trace(&old, &law);
    let change = old
        .refine_controlled(&MORE, budget(), &mut || false)
        .unwrap();
    let next = change.topology();
    let pairs: Vec<_> = next
        .seams()
        .iter()
        .map(|seam| seam.face_pair(&law))
        .collect();
    let correct = trace(next, &law);
    let canonical = trace_at(next.nodes_m(), next.tetrahedra(), &pairs, false);
    assert_ne!(canonical.points(), correct.points());
    assert!(
        change
            .require_retained_interfaces_controlled(&old_map, &canonical, &mut || false)
            .is_err()
    );
    let mut reversed = pairs.clone();
    let first_reversed = &mut reversed[0];
    std::mem::swap(&mut first_reversed.side_a, &mut first_reversed.side_b);
    let reversed = trace_at(next.nodes_m(), next.tetrahedra(), &reversed, true);
    assert!(
        change
            .require_retained_interfaces_controlled(&old_map, &reversed, &mut || false)
            .is_err()
    );
    let mut renamed = pairs.clone();
    renamed[0].id += 100;
    let renamed = trace_at(next.nodes_m(), next.tetrahedra(), &renamed, true);
    assert!(
        change
            .require_retained_interfaces_controlled(&old_map, &renamed, &mut || false)
            .is_err()
    );
    let incomplete = trace_at(next.nodes_m(), next.tetrahedra(), &pairs[..2], true);
    assert!(
        change
            .require_retained_interfaces_controlled(&old_map, &incomplete, &mut || false)
            .is_err()
    );
    // Both supplied maps agree with EACH OTHER after translation, but neither
    // belongs to the actual topology. Pairwise agreement alone is insufficient.
    let shift = |nodes: &[[f64; 3]]| {
        nodes
            .iter()
            .map(|x| [x[0] + 0.125, x[1] - 0.25, x[2] + 0.5])
            .collect::<Vec<_>>()
    };
    let old_pairs: Vec<_> = old
        .seams()
        .iter()
        .map(|seam| seam.face_pair(&law))
        .collect();
    let shifted_old = trace_at(&shift(old.nodes_m()), old.tetrahedra(), &old_pairs, true);
    let shifted_next = trace_at(&shift(next.nodes_m()), next.tetrahedra(), &pairs, true);
    assert!(
        change
            .require_retained_interfaces_controlled(&shifted_old, &shifted_next, &mut || false)
            .is_err()
    );
    change
        .require_retained_interfaces_controlled(&old_map, &correct, &mut || false)
        .unwrap();
}

fn source() -> ResolvedBilinearCohesive {
    use fs_evidence::ValidityDomain;
    use fs_matdb::{
        ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialStateId, PropertyClaim,
        PropertyKey, PropertyValue, Provenance, QueryPoint, SurfaceSpec, SystemContext,
        UncertaintyModel,
    };
    use fs_material::cohesive::source::*;
    use fs_material::state_point::{MaterialPropertySelection, interface_stiffness_per_area_dims};
    use fs_qty::Pressure;
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        (
            COHESIVE_STIFFNESS_PROPERTY,
            interface_stiffness_per_area_dims(),
            1e9,
        ),
        (COHESIVE_STRENGTH_PROPERTY, Pressure::DIMS, 2e5),
        (
            COHESIVE_FRACTURE_ENERGY_PROPERTY,
            FRACTURE_ENERGY_PER_AREA_DIMS,
            100.,
        ),
        (
            COHESIVE_COMPRESSION_STIFFNESS_PROPERTY,
            interface_stiffness_per_area_dims(),
            2e9,
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained(),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured retained-interface test; no measured claim".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let surface = |chemistry: &str| SurfaceSpec {
        material: MaterialStateId {
            chemistry: chemistry.into(),
            phase: "solid".into(),
            process: "test".into(),
            revision: 0,
        },
        texture_frame: "fixed reference frame".into(),
    };
    let card = InterfaceSystemCard::assemble(
        surface("A"),
        surface("B"),
        SystemContext {
            medium: "explicit cohesive seam".into(),
            third_body: None,
            environment: "fixed".into(),
            history: "virgin".into(),
        },
        claims,
        vec![],
    )
    .unwrap();
    resolve_bilinear_cohesive(
        &card,
        &QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly,
        0.01,
    )
    .unwrap()
}
fn sourced_trace<'a>(
    topology: &TetFractureTopology<'_>,
    source: &'a ResolvedBilinearCohesive,
) -> P1Interfaces<'a> {
    let pairs: Vec<_> = topology
        .seams()
        .iter()
        .map(|seam| MaterialFacePair {
            id: seam.id,
            side_a: seam.side_a,
            side_b: seam.side_b,
            source,
            surface_materials: [
                &source.surfaces()[0].material,
                &source.surfaces()[1].material,
            ],
        })
        .collect();
    let free: Vec<_> = (0..3 * topology.nodes_m().len()).collect();
    P1Interfaces::from_interfaces_ordered_controlled(
        topology.nodes_m(),
        topology.tetrahedra(),
        &free,
        &pairs,
        trace_budget(),
        &mut || false,
    )
    .unwrap()
}

#[test]
fn equal_coefficients_do_not_replace_original_law_or_interface_source_owners() {
    let law = law();
    let cloned = law.clone();
    assert_eq!(law, cloned);
    let old = first();
    let change = old
        .refine_controlled(&MORE, budget(), &mut || false)
        .unwrap();
    assert!(
        change
            .require_retained_interfaces_controlled(
                &trace(&old, &law),
                &trace(change.topology(), &cloned),
                &mut || false
            )
            .is_err()
    );
    let source = source();
    let clone = source.clone();
    assert_eq!(source, clone);
    let previous = sourced_trace(&old, &source);
    let next = sourced_trace(change.topology(), &source);
    change
        .require_retained_interfaces_controlled(&previous, &next, &mut || false)
        .unwrap();
    assert!(
        change
            .require_retained_interfaces_controlled(
                &previous,
                &sourced_trace(change.topology(), &clone),
                &mut || false
            )
            .is_err()
    );
    let raw = trace(change.topology(), source.law());
    assert!(core::ptr::eq(
        previous.patches()[0].law(),
        raw.patches()[0].law()
    ));
    assert!(
        change
            .require_retained_interfaces_controlled(&previous, &raw, &mut || false)
            .is_err()
    );
}

#[test]
fn cancellation_at_admission_middle_and_publication_leaves_histories_reusable() {
    let law = law();
    let old = first();
    let previous = trace(&old, &law);
    let change = old
        .refine_controlled(&MORE, budget(), &mut || false)
        .unwrap();
    let next = trace(change.topology(), &law);
    let u = opening(&old, &previous, 0.0011);
    let history: Vec<_> = previous
        .jumps(&u, &mut || false)
        .unwrap()
        .into_iter()
        .map(|jump| law.bonded_point().trial(jump, 1e-7).unwrap().point)
        .collect();
    let saved = history.clone();
    let old_points = previous.points().to_vec();
    let new_points = next.points().to_vec();
    let mut polls = 0;
    change
        .require_retained_interfaces_controlled(&previous, &next, &mut || {
            polls += 1;
            false
        })
        .unwrap();
    assert!(polls >= 3);
    for stop in [1, polls / 2, polls] {
        let mut count = 0;
        assert_eq!(
            change.require_retained_interfaces_controlled(&previous, &next, &mut || {
                count += 1;
                count == stop
            }),
            Err(FractureTopologyError::Cancelled)
        );
        assert_eq!(history, saved);
        assert_eq!(previous.points(), old_points);
        assert_eq!(next.points(), new_points);
        change
            .require_retained_interfaces_controlled(&previous, &next, &mut || false)
            .unwrap();
    }
}
