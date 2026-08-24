//! Feature-complex / CCD-candidate battery (bead rjnd, part 3).
//!
//! - gf-001 G0: complex construction — feature counts, canonical edge
//!   dedup, deterministic order.
//! - gf-002 G0: BVH enumeration equals brute-force inflated-box
//!   overlap exactly (same set, same order after sorting), and a pair
//!   that must collide within the motion window is present.
//! - gf-003 G3: translating both complexes identically preserves the
//!   candidate set (indices), and widening a motion bound never
//!   removes candidates (monotone superset).
//! - gf-004 G5: identical inputs replay identical candidate vectors.
//! - gf-005 G0/G4: index/position/inflation/cap/cancellation refusals
//!   fail closed with the named typed error.
//! Aggregate outcomes use canonical fs-obs events; evaluated cases carry
//! the shared execution seed and constructor-only gf-001 uses zero.

use asupersync::types::Budget;
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_query::{
    Feature, FeatureComplex, MAX_COMPLEX_FEATURES, QueryError, ccd_candidates,
};

const EXECUTION_SEED: u64 = 0xFEA7;

fn verdict(case: &str, pass: bool, detail: &str, seed: u64) {
    let mut emitter = fs_obs::Emitter::new("fs-query/features", case);
    let event = emitter.emit(
        if pass {
            fs_obs::Severity::Info
        } else {
            fs_obs::Severity::Error
        },
        fs_obs::EventKind::ConformanceCase {
            suite: "fs-query/features".to_string(),
            case: case.to_string(),
            pass,
            detail: detail.to_string(),
            seed,
        },
        None,
    );
    fs_obs::lint_failure_record(&event).expect("feature-candidate verdict must be replayable");
    let line = event.to_jsonl();
    fs_obs::validate_line(&line).expect("feature verdict must use the fs-obs wire schema");
    println!("{line}");
    assert!(pass, "case {case}: {detail}");
}

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: EXECUTION_SEED,
                kernel_id: 14,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

/// A tetrahedron with one vertex at `origin`, edge scale `s`.
fn tetra(origin: [f64; 3], s: f64) -> (Vec<[f64; 3]>, Vec<[u32; 3]>) {
    let [x, y, z] = origin;
    let positions = vec![[x, y, z], [x + s, y, z], [x, y + s, z], [x, y, z + s]];
    let triangles = vec![[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]];
    (positions, triangles)
}

fn shifted(positions: &[[f64; 3]], d: [f64; 3]) -> Vec<[f64; 3]> {
    positions
        .iter()
        .map(|p| [p[0] + d[0], p[1] + d[1], p[2] + d[2]])
        .collect()
}

/// Brute-force oracle: inflate each side's feature boxes and test all
/// pairs. Mirrors the documented semantics exactly.
fn brute_pairs(
    a_pos: &[[f64; 3]],
    a_tri: &[[u32; 3]],
    b_pos: &[[f64; 3]],
    b_tri: &[[u32; 3]],
    motion_a: f64,
    motion_b: f64,
) -> Vec<(usize, usize)> {
    let feature_points = |pos: &[[f64; 3]], tri: &[[u32; 3]]| -> Vec<Vec<[f64; 3]>> {
        let mut edges: Vec<(u32, u32)> = Vec::new();
        for t in tri {
            for (a, b) in [(t[0], t[1]), (t[1], t[2]), (t[0], t[2])] {
                edges.push((a.min(b), a.max(b)));
            }
        }
        edges.sort_unstable();
        edges.dedup();
        let mut out: Vec<Vec<[f64; 3]>> = Vec::new();
        for p in pos {
            out.push(vec![*p]);
        }
        for (a, b) in edges {
            out.push(vec![pos[a as usize], pos[b as usize]]);
        }
        for t in tri {
            out.push(vec![
                pos[t[0] as usize],
                pos[t[1] as usize],
                pos[t[2] as usize],
            ]);
        }
        out
    };
    let boxes = |sets: &[Vec<[f64; 3]>], pad: f64| -> Vec<([f64; 3], [f64; 3])> {
        sets.iter()
            .map(|points| {
                let mut min = [f64::INFINITY; 3];
                let mut max = [f64::NEG_INFINITY; 3];
                for p in points {
                    for a in 0..3 {
                        min[a] = min[a].min(p[a]);
                        max[a] = max[a].max(p[a]);
                    }
                }
                for a in 0..3 {
                    min[a] = (min[a].next_down() - pad).next_down();
                    max[a] = (max[a].next_up() + pad).next_up();
                }
                (min, max)
            })
            .collect()
    };
    let ba = boxes(&feature_points(a_pos, a_tri), motion_a);
    let bb = boxes(&feature_points(b_pos, b_tri), motion_b);
    let mut pairs = Vec::new();
    for (i, (amin, amax)) in ba.iter().enumerate() {
        for (j, (bmin, bmax)) in bb.iter().enumerate() {
            let overlap = (0..3).all(|k| amin[k] <= bmax[k] && bmin[k] <= amax[k]);
            if overlap {
                pairs.push((i, j));
            }
        }
    }
    pairs
}

#[test]
fn gf_001_complex_construction_is_canonical() {
    let (pos, tri) = tetra([0.0, 0.0, 0.0], 1.0);
    let complex = FeatureComplex::from_triangles(&pos, &tri).expect("tetra complex");
    // 4 vertices + 6 unique edges + 4 faces.
    assert_eq!(complex.len(), 14);
    assert_eq!(complex.feature(0), Some(Feature::Vertex(0)));
    assert_eq!(complex.feature(4), Some(Feature::Edge(0, 1)));
    assert_eq!(complex.feature(10), Some(Feature::Face(0)));
    let again = FeatureComplex::from_triangles(&pos, &tri).expect("rebuild");
    let same = (0..complex.len()).all(|i| complex.feature(i) == again.feature(i));
    verdict(
        "gf-001",
        same,
        "tetra complex: 4+6+4 features in deterministic canonical order",
        0,
    );
}

#[test]
fn gf_002_bvh_matches_brute_force_and_finds_the_hit() {
    let (a_pos, a_tri) = tetra([0.0, 0.0, 0.0], 1.0);
    // Second tetra separated by a 0.3 gap along +x: within a 0.2+0.2
    // motion window the near features MUST be candidates.
    let (b_pos, b_tri) = tetra([1.3, 0.0, 0.0], 1.0);
    let a = FeatureComplex::from_triangles(&a_pos, &a_tri).expect("a");
    let b = FeatureComplex::from_triangles(&b_pos, &b_tri).expect("b");
    let pairs = with_cx(|cx| ccd_candidates(&a, &b, 0.2, 0.2, 100_000, cx)).expect("candidates");
    let mut oracle = brute_pairs(&a_pos, &a_tri, &b_pos, &b_tri, 0.2, 0.2);
    oracle.sort_unstable();
    assert_eq!(
        pairs, oracle,
        "BVH enumeration must equal the brute-force oracle exactly"
    );
    // Vertex 1 of A ([1,0,0]) and vertex 0 of B ([1.3,0,0]) close to a
    // 0.3 gap: with 0.4 total inflation they must appear.
    assert!(
        pairs.contains(&(1, 0)),
        "the closest vertex pair must be a candidate"
    );
    // And with a tiny window the far corners must NOT be candidates.
    let tight = with_cx(|cx| ccd_candidates(&a, &b, 0.01, 0.01, 100_000, cx)).expect("tight");
    assert!(
        tight.is_empty(),
        "a 0.02 total window cannot bridge a 0.3 gap, got {} pairs",
        tight.len()
    );
    verdict(
        "gf-002",
        true,
        &format!(
            "BVH == brute force ({} pairs at 0.4 window; 0 at 0.02 window)",
            pairs.len()
        ),
        EXECUTION_SEED,
    );
}

#[test]
fn gf_003_translation_invariance_and_monotone_windows() {
    let (a_pos, a_tri) = tetra([0.0, 0.0, 0.0], 1.0);
    let (b_pos, b_tri) = tetra([1.1, 0.2, -0.1], 1.0);
    let d = [3.5, -2.25, 0.625];
    let a = FeatureComplex::from_triangles(&a_pos, &a_tri).expect("a");
    let b = FeatureComplex::from_triangles(&b_pos, &b_tri).expect("b");
    let a_shift = FeatureComplex::from_triangles(&shifted(&a_pos, d), &a_tri).expect("a shifted");
    let b_shift = FeatureComplex::from_triangles(&shifted(&b_pos, d), &b_tri).expect("b shifted");
    let (base, moved, narrow, wide) = with_cx(|cx| {
        (
            ccd_candidates(&a, &b, 0.15, 0.15, 100_000, cx).expect("base"),
            ccd_candidates(&a_shift, &b_shift, 0.15, 0.15, 100_000, cx).expect("moved"),
            ccd_candidates(&a, &b, 0.05, 0.05, 100_000, cx).expect("narrow"),
            ccd_candidates(&a, &b, 0.45, 0.45, 100_000, cx).expect("wide"),
        )
    });
    // Dyadic translation: box arithmetic is exact, so the candidate
    // sets agree exactly. (A general translation could flip borderline
    // overlaps within one ulp; the semantics stay conservative.)
    assert_eq!(base, moved, "dyadic rigid translation preserves candidates");
    for pair in &narrow {
        assert!(
            base.contains(pair),
            "widening the window must never remove candidate {pair:?}"
        );
    }
    for pair in &base {
        assert!(
            wide.contains(pair),
            "widening the window must never remove candidate {pair:?}"
        );
    }
    verdict(
        "gf-003",
        true,
        &format!(
            "translation-invariant ({} pairs); windows monotone {} ⊆ {} ⊆ {}",
            base.len(),
            narrow.len(),
            base.len(),
            wide.len()
        ),
        EXECUTION_SEED,
    );
}

#[test]
fn gf_004_replay_is_identical() {
    let (a_pos, a_tri) = tetra([0.0, 0.0, 0.0], 1.0);
    let (b_pos, b_tri) = tetra([0.9, 0.1, 0.2], 1.0);
    let a = FeatureComplex::from_triangles(&a_pos, &a_tri).expect("a");
    let b = FeatureComplex::from_triangles(&b_pos, &b_tri).expect("b");
    let first = with_cx(|cx| ccd_candidates(&a, &b, 0.1, 0.1, 100_000, cx)).expect("first");
    let second = with_cx(|cx| ccd_candidates(&a, &b, 0.1, 0.1, 100_000, cx)).expect("second");
    assert_eq!(first, second);
    verdict(
        "gf-004",
        true,
        &format!("{} candidate pairs replay identically", first.len()),
        EXECUTION_SEED,
    );
}

#[test]
fn gf_005_refusals_fail_closed() {
    let (pos, tri) = tetra([0.0, 0.0, 0.0], 1.0);
    let complex = FeatureComplex::from_triangles(&pos, &tri).expect("valid");

    let bad_index = FeatureComplex::from_triangles(&pos, &[[0, 1, 9]]);
    assert!(matches!(
        bad_index,
        Err(QueryError::InvalidBoundaryIndex { index: 9, .. })
    ));
    let degenerate = FeatureComplex::from_triangles(&pos, &[[0, 0, 1]]);
    assert!(matches!(
        degenerate,
        Err(QueryError::InvalidBoundaryIndex { .. })
    ));
    let bad_pos = FeatureComplex::from_triangles(&[[0.0, f64::NAN, 0.0]], &[]);
    assert!(matches!(
        bad_pos,
        Err(QueryError::InvalidPointSample { .. })
    ));
    let oversized_triangles = vec![[0, 0, 0]; MAX_COMPLEX_FEATURES + 1];
    assert!(
        matches!(
            FeatureComplex::from_triangles(&[], &oversized_triangles),
            Err(QueryError::FeatureComplexTooLarge { features, max })
                if features == MAX_COMPLEX_FEATURES + 1 && max == MAX_COMPLEX_FEATURES
        ),
        "the feature ceiling must refuse before malformed triangle validation or edge allocation"
    );

    let (nan_inflation, negative_inflation, capped, cancelled_ok) = with_cx(|cx| {
        (
            ccd_candidates(&complex, &complex, f64::NAN, 0.0, 100, cx),
            ccd_candidates(&complex, &complex, 0.1, -0.1, 100, cx),
            ccd_candidates(&complex, &complex, 0.1, 0.1, 3, cx),
            ccd_candidates(&complex, &complex, 0.1, 0.1, 100_000, cx),
        )
    });
    assert!(matches!(
        nan_inflation,
        Err(QueryError::FeatureInvalidInflation { .. })
    ));
    assert!(matches!(
        negative_inflation,
        Err(QueryError::FeatureInvalidInflation { .. })
    ));
    assert!(matches!(
        capped,
        Err(QueryError::FeatureTooManyPairs { max: 3 })
    ));
    assert!(cancelled_ok.is_ok(), "self-overlap enumerates fine");

    let gate = CancelGate::new();
    gate.request();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let cancelled = pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: EXECUTION_SEED,
                kernel_id: 15,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        ccd_candidates(&complex, &complex, 0.1, 0.1, 100_000, &cx)
    });
    assert!(matches!(cancelled, Err(QueryError::Cancelled)));
    verdict(
        "gf-005",
        true,
        "index/position/inflation/cap/cancellation all refuse typed",
        EXECUTION_SEED,
    );
}
