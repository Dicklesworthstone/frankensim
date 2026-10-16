//! Claim-integrity regressions for the mesh certificates
//! (`docs/CLAIM_INTEGRITY.md`).
//!
//! Every test here encodes the minimal repro from one filed defect and
//! fails against the pre-fix behavior:
//!
//! - `frankensim-extreal-program-f85xj.2.14` — vertex-sharing face pairs
//!   were skipped wholesale, so a face piercing a neighbour it touched
//!   was reported `proven_free`;
//! - `.2.15` — non-finite coordinates passed BOTH certificates: NaN is
//!   dropped by the AABB `min`/`max` broad phase and falsifies every
//!   red-flag comparison in the manifold path;
//! - `.2.16` — a GLOBAL inward-orientation failure was published as a
//!   localized `MisorientedEdge { edge: [0, 0] }`, a self-loop that is
//!   not an edge of any mesh;
//! - `.2.17` — `oriented` documented outwardness but omitted it whenever
//!   no probe was supplied, and the report kept no record of which
//!   meaning the flag carried.

use fs_geom::Point3;
use fs_rep_mesh::{Soup, shapes};
use fs_topo::{
    IntersectKind, ManifoldDefect, SelfIntersectRefusal, manifold_certificate,
    self_intersection_certificate,
};

fn soup(positions: &[[f64; 3]], triangles: &[[u32; 3]]) -> Soup {
    Soup {
        positions: positions
            .iter()
            .map(|p| Point3::new(p[0], p[1], p[2]))
            .collect(),
        triangles: triangles.to_vec(),
    }
}

fn reversed(mut s: Soup) -> Soup {
    for t in &mut s.triangles {
        t.swap(1, 2);
    }
    s
}

/// The closed unit corner tetrahedron, outward-oriented. Every one of its
/// six face pairs shares an edge, so it is the sharpest available check
/// that legitimate adjacency still passes.
fn tetrahedron() -> Soup {
    soup(
        &[
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ],
        &[[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]],
    )
}

// ---------------------------------------------------------------- .2.14

/// A face that shares one vertex with its neighbour and pierces it is a
/// genuine self-intersection. Before the fix the pair was excluded by
/// `if ti.iter().any(|v| tj.contains(v)) { continue; }`, so
/// `pairs_tested == 0` and `proven_free() == true`.
#[test]
fn vertex_sharing_pair_that_pierces_is_not_proven_free() {
    // T0 spans the z = 0 corner triangle; T1 is a blade through it,
    // hinged on the shared vertex 0. The blade's edge 3→4 crosses
    // (0.3, 0.2, 0), which is strictly inside T0 because 0.3 + 0.2 < 1.
    let pierced = soup(
        &[
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.3, 0.2, -1.0],
            [0.3, 0.2, 1.0],
        ],
        &[[0, 1, 2], [0, 3, 4]],
    );
    let report = self_intersection_certificate(&pierced);
    assert!(
        report.admitted(),
        "the repro is finite, in range and non-degenerate: {:?}",
        report.refusals
    );
    assert_eq!(
        report.shared_feature_pairs_tested, 1,
        "the vertex-sharing pair must be DECIDED, not skipped"
    );
    assert!(
        !report.proven_free(),
        "a blade piercing the face it shares a vertex with is a self-intersection"
    );
    assert!(
        report
            .intersections
            .iter()
            .any(|&(a, b, k)| (a, b) == (0, 1) && k == IntersectKind::Crossing),
        "the piercing is strict, so it localizes as Crossing on pair (0, 1): {:?}",
        report.intersections
    );
}

/// The other half of the claim: deciding shared-feature pairs must not
/// manufacture false FAILs on the adjacency that legitimately touches.
#[test]
fn legitimate_adjacency_still_passes_and_is_actually_tested() {
    let report = self_intersection_certificate(&tetrahedron());
    assert_eq!(
        report.shared_feature_pairs_tested, 6,
        "all six tetrahedron face pairs share an edge and must be decided"
    );
    assert!(
        report.proven_free(),
        "a tetrahedron meets itself only along its shared edges: {:?}",
        report.intersections
    );

    // Vertex-sharing (not edge-sharing) adjacency, at scale.
    let sphere = shapes::icosphere(Point3::new(0.0, 0.0, 0.0), 1.0, 3);
    let report = self_intersection_certificate(&sphere);
    assert!(
        report.shared_feature_pairs_tested > 0,
        "the icosphere's adjacent pairs must reach the narrow phase"
    );
    assert!(
        report.proven_free(),
        "a clean icosphere stays PROVEN free: {:?}",
        report.intersections
    );
}

/// Coplanar shared-feature pairs go through the exact 2D corner-cone and
/// hinge-side tests rather than the non-coplanar chord argument.
#[test]
fn coplanar_shared_feature_pairs_are_decided_exactly() {
    // Shared vertex, cones disjoint: they meet only at the origin.
    let disjoint = soup(
        &[
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, -1.0, 0.0],
            [-1.0, -1.0, 0.0],
        ],
        &[[0, 1, 2], [0, 3, 4]],
    );
    assert!(
        self_intersection_certificate(&disjoint).proven_free(),
        "coplanar wedges meeting only at the shared corner are clean"
    );

    // Shared vertex, cones overlapping: a genuine coplanar overlap.
    let overlapping = soup(
        &[
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
            [-1.0, 1.0, 0.0],
        ],
        &[[0, 1, 2], [0, 3, 4]],
    );
    assert!(
        !self_intersection_certificate(&overlapping).proven_free(),
        "overlapping coplanar corner cones are a self-intersection"
    );

    // Shared edge, opposite sides: a flat hinge, clean.
    let hinge = soup(
        &[
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
        ],
        &[[0, 1, 2], [0, 1, 3]],
    );
    assert!(
        self_intersection_certificate(&hinge).proven_free(),
        "two coplanar faces hinged on a shared edge with opposite apexes are clean"
    );

    // Shared edge, same side: positive-area overlap.
    let folded = soup(
        &[
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.2, 0.5, 0.0],
        ],
        &[[0, 1, 2], [0, 1, 3]],
    );
    assert!(
        !self_intersection_certificate(&folded).proven_free(),
        "two coplanar faces on the SAME side of their shared edge overlap"
    );
}

// ---------------------------------------------------------------- .2.15

/// NaN coordinates are dropped by `f64::min`/`f64::max`, so the offending
/// triangle got the inverted box `lo = +inf, hi = -inf` and every pair
/// containing it was culled before the exact narrow phase. The result was
/// `pairs_tested == 0` and `proven_free() == true` on all-NaN geometry.
#[test]
fn non_finite_coordinates_refuse_the_self_intersection_certificate() {
    let nan = f64::NAN;
    let bad = soup(
        &[
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.2, nan, -1.0],
            [0.4, nan, 1.0],
            [0.3, nan, 0.0],
        ],
        &[[0, 1, 2], [3, 4, 5]],
    );
    let report = self_intersection_certificate(&bad);
    assert!(!report.admitted(), "a NaN soup must not be admitted");
    assert!(
        !report.proven_free(),
        "a refused soup proves nothing about self-intersection"
    );
    assert_eq!(
        report.refusals,
        vec![
            SelfIntersectRefusal::NonFiniteVertex { vertex: 3 },
            SelfIntersectRefusal::NonFiniteVertex { vertex: 4 },
            SelfIntersectRefusal::NonFiniteVertex { vertex: 5 },
        ],
        "the refusal names every offending vertex, in index order"
    );
    assert_eq!(report.pairs_tested, 0, "no predicate may run on NaN input");
}

/// Out-of-range indices used to panic inside `positions[..]`; they are an
/// admission refusal now, so the certificate stays total.
#[test]
fn out_of_range_indices_refuse_the_self_intersection_certificate() {
    let bad = soup(
        &[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        &[[0, 1, 2], [0, 1, 7]],
    );
    let report = self_intersection_certificate(&bad);
    assert!(!report.proven_free());
    assert_eq!(
        report.refusals,
        vec![SelfIntersectRefusal::VertexIndexOutOfRange { face: 1, vertex: 7 }]
    );
}

/// An exactly-degenerate face has no supporting plane, so every
/// plane-separation argument in the narrow phase is vacuous on it.
#[test]
fn exactly_degenerate_faces_refuse_rather_than_pass() {
    let sliver = soup(
        &[
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ],
        &[[0, 1, 3], [0, 1, 2]],
    );
    let report = self_intersection_certificate(&sliver);
    assert!(
        !report.proven_free(),
        "a zero-area face cannot be argued about"
    );
    assert_eq!(
        report.refusals,
        vec![SelfIntersectRefusal::DegenerateFace { face: 1 }]
    );
}

/// `manifold_certificate` never read positions in a way NaN could fail:
/// `n.norm() < 1e-30`, `den > 1e-30` and `(w - 1.0).abs() > 0.5` are all
/// FALSE for NaN, so an all-NaN combinatorial sphere certified.
#[test]
fn non_finite_coordinates_refuse_the_manifold_certificate() {
    let mut nan_sphere = shapes::icosphere(Point3::new(0.0, 0.0, 0.0), 1.0, 1);
    for p in &mut nan_sphere.positions {
        *p = Point3::new(f64::NAN, f64::NAN, f64::NAN);
    }
    let report = manifold_certificate(&nan_sphere, Some(Point3::new(0.0, 0.0, 0.0)));
    assert!(
        !report.certified(),
        "NaN geometry cannot certify as a closed oriented manifold surface"
    );
    assert!(
        report
            .defects
            .iter()
            .any(|d| matches!(d, ManifoldDefect::NonFiniteVertex { .. })),
        "the refusal must be named: {:?}",
        report.defects
    );
    assert_eq!(
        report.outward, None,
        "no outwardness verdict may be minted from a NaN winding sum"
    );
    // The combinatorial verdicts are honest and survive: they never read
    // a coordinate.
    assert!(report.manifold && report.closed && report.consistently_oriented);
}

/// A non-finite probe must not silently pass the winding comparison
/// either.
#[test]
fn non_finite_probe_yields_no_outwardness_verdict() {
    let sphere = shapes::icosphere(Point3::new(0.0, 0.0, 0.0), 1.0, 1);
    let report = manifold_certificate(&sphere, Some(Point3::new(f64::NAN, 0.0, 0.0)));
    assert_eq!(report.outward, None);
    assert!(!report.certified());
    assert!(
        report
            .defects
            .iter()
            .any(|d| matches!(d, ManifoldDefect::IndeterminateWinding { .. })),
        "{:?}",
        report.defects
    );
}

// ---------------------------------------------------------------- .2.16

/// Reversing every triangle keeps the surface consistently oriented and
/// leaves no misoriented edge — the edge-use census proves `dir_sum == 0`
/// everywhere before the winding probe runs. The old code nonetheless
/// fabricated `MisorientedEdge { edge: [0, 0] }`, a self-loop on vertex 0
/// that no mesh contains, inviting a local repair at a non-existent edge.
#[test]
fn inward_orientation_is_reported_globally_not_as_a_fabricated_edge() {
    let inward = reversed(shapes::icosphere(Point3::new(0.0, 0.0, 0.0), 1.0, 2));
    let probe = Point3::new(0.0, 0.0, 0.0);
    let report = manifold_certificate(&inward, Some(probe));

    assert!(
        !report
            .defects
            .iter()
            .any(|d| matches!(d, ManifoldDefect::MisorientedEdge { .. })),
        "no edge is misoriented on a uniformly reversed surface: {:?}",
        report.defects
    );
    assert!(
        report.consistently_oriented,
        "uniform reversal preserves combinatorial orientation consistency"
    );
    let inward_defect = report
        .defects
        .iter()
        .find_map(|d| match d {
            ManifoldDefect::InwardOrientation { probe, winding } => Some((*probe, *winding)),
            _ => None,
        })
        .unwrap_or_else(|| panic!("expected a global InwardOrientation: {:?}", report.defects));
    assert_eq!(
        inward_defect.0, probe,
        "the defect carries the actual probe"
    );
    assert!(
        (inward_defect.1 + 1.0).abs() < 1e-6,
        "a reversed closed surface winds -1 at an interior point, got {}",
        inward_defect.1
    );
    assert_eq!(report.outward, Some(false));
    assert!(!report.certified());
}

// ---------------------------------------------------------------- .2.17

/// Without a probe there is NO outwardness evidence, and the report must
/// say so instead of publishing a flag whose documentation promises it.
#[test]
fn an_unprobed_report_makes_no_outwardness_claim() {
    let inward = reversed(shapes::icosphere(Point3::new(0.0, 0.0, 0.0), 1.0, 2));
    let report = manifold_certificate(&inward, None);

    assert_eq!(
        report.outward, None,
        "no probe ran, so outwardness is unknown — not true"
    );
    assert!(
        report.combinatorially_certified(),
        "the combinatorial claim is genuine and stays available"
    );
    assert!(
        !report.certified(),
        "certified() must not pass an entirely inward-facing surface"
    );
}

/// The positive control: a probed clean sphere certifies, and the same
/// sphere without a probe reports the weaker claim rather than the strong
/// one.
#[test]
fn a_probed_clean_surface_still_certifies() {
    let sphere = shapes::icosphere(Point3::new(0.0, 0.0, 0.0), 1.0, 2);
    let probed = manifold_certificate(&sphere, Some(Point3::new(0.0, 0.0, 0.0)));
    assert!(probed.certified());
    assert_eq!(probed.outward, Some(true));
    assert!(probed.defects.is_empty(), "{:?}", probed.defects);

    let unprobed = manifold_certificate(&sphere, None);
    assert!(unprobed.combinatorially_certified());
    assert!(!unprobed.certified());
    assert_eq!(unprobed.outward, None);

    let tetra = manifold_certificate(&tetrahedron(), Some(Point3::new(0.2, 0.2, 0.2)));
    assert!(tetra.certified(), "{:?}", tetra.defects);
}
