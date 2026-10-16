//! fs-rep-voxel conformance (the wqd.6 bead). Acceptance: morphology
//! matches brute force; the Euclidean DT is EXACT vs the O(n²) reference;
//! point-cloud queries match brute force and normals recover analytic
//! surfaces; lattice graphs round-trip through FrankenNetworkx with
//! attributes preserved, degenerates refuse structurally, and realization
//! behaves like a watertight solid's level set.

use fs_geom::Chart;
use fs_rep_voxel::{
    LatticeGraph, LatticeNode, OccupancyChart, OccupancyField, PointCloud, Strut, VoxelError,
};

const TEST_DT_BUDGET: usize = 1_000_000;

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-rep-voxel/conformance\",\"case\":\"{case}\",\"verdict\":\"pass\",\
         \"detail\":\"{detail}\"}}"
    );
}

fn lcg(seed: &mut u64) -> f64 {
    *seed = seed
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    ((*seed >> 11) as f64) / (1u64 << 53) as f64
}

/// A solid ball of voxels.
fn ball(radius: i32) -> OccupancyField {
    let mut f = OccupancyField::new(0.1, [0.0; 3]).expect("field");
    for x in -radius..=radius {
        for y in -radius..=radius {
            for z in -radius..=radius {
                if x * x + y * y + z * z <= radius * radius {
                    f.set([x, y, z]);
                }
            }
        }
    }
    f
}

/// Brute-force 6-connected dilation of an active set.
fn brute_dilate(
    cells: &std::collections::BTreeSet<[i32; 3]>,
) -> std::collections::BTreeSet<[i32; 3]> {
    let mut out = cells.clone();
    for c in cells {
        for d in [
            [-1, 0, 0],
            [1, 0, 0],
            [0, -1, 0],
            [0, 1, 0],
            [0, 0, -1],
            [0, 0, 1],
        ] {
            out.insert([c[0] + d[0], c[1] + d[1], c[2] + d[2]]);
        }
    }
    out
}

fn active_set(f: &OccupancyField) -> std::collections::BTreeSet<[i32; 3]> {
    f.grid.iter_active().map(|(c, _)| c).collect()
}

#[test]
fn rv_001_morphology_matches_brute_force_and_algebra() {
    let mut f = ball(4);
    // A thin spur that opening must remove.
    for k in 5..9 {
        f.set([k, 0, 0]);
    }
    let before = active_set(&f);
    // Dilation vs brute force.
    let mut dilated = f.clone();
    dilated.dilate(1);
    assert_eq!(
        active_set(&dilated),
        brute_dilate(&before),
        "dilate == brute force"
    );
    // Erosion is dilation's dual on the fixtures: erode(dilate(A)) ⊇ A.
    let mut closed = f.clone();
    closed.close(1);
    let after_close = active_set(&closed);
    assert!(
        before.is_subset(&after_close),
        "closing must never remove original voxels"
    );
    // Opening ⊆ A, and it removes the 1-voxel spur but keeps the ball.
    let mut opened = f.clone();
    opened.open(1);
    let after_open = active_set(&opened);
    assert!(after_open.is_subset(&before), "opening must not add voxels");
    assert!(
        !after_open.contains(&[7, 0, 0]),
        "the spur must be opened away"
    );
    assert!(
        after_open.contains(&[0, 0, 0]),
        "the ball core survives opening"
    );
    // Boolean algebra: (A ∪ B) \ B ⊆ A, A ∩ B ⊆ A.
    let b = ball(3);
    let mut u = f.clone();
    u.union(&b).expect("matching frames");
    u.subtract(&b).expect("matching frames");
    assert!(active_set(&u).is_subset(&before));
    let mut i = f.clone();
    i.intersect(&b).expect("matching frames");
    assert!(active_set(&i).is_subset(&active_set(&b)));
    println!("{}", f.stats_json());
    verdict(
        "rv-001",
        "dilate == brute force; open/close/boolean algebra laws hold; spur removed",
    );
}

#[test]
fn rv_002_euclidean_dt_is_exact_vs_reference() {
    // Sparse scattered seeds + a slab: exactness on non-trivial topology.
    let mut f = OccupancyField::new(1.0, [0.0; 3]).expect("field");
    let mut seed = 0xD7_0001u64;
    let mut seeds = Vec::new();
    for _ in 0..24 {
        let c = [
            (lcg(&mut seed) * 20.0) as i32,
            (lcg(&mut seed) * 14.0) as i32,
            (lcg(&mut seed) * 10.0) as i32,
        ];
        f.set(c);
        seeds.push(c);
    }
    for x in 3..9 {
        for y in 2..5 {
            f.set([x, y, 7]);
            seeds.push([x, y, 7]);
        }
    }
    seeds.sort_unstable();
    seeds.dedup();
    let dt = fs_rep_voxel::euclidean_dt(&f, TEST_DT_BUDGET)
        .expect("admissible box")
        .expect("nonempty");
    let dt_min = dt.min();
    let dt_dims = dt.dims();
    // O(n²) reference over the whole box, gated at the "EXACT" the case name
    // promises — BIT-FOR-BIT, not a tolerance. `brute` is an exact-integer
    // squared distance (i32 products); the DT stores the SAME exact-integer
    // squared distance (guarded <= 2^53) and returns `det::sqrt(sq) *
    // voxel_size`, where `det::sqrt` IS the correctly-rounded IEEE sqrt. So an
    // exact DT reproduces `sqrt(brute) * voxel_size` bit-for-bit; a wrong
    // nearest seed, a lower-envelope off-by-one, or a voxel-scaling slip flips
    // at least one bit. (The old `< 1e-9` bound was ~7 orders looser than the
    // integer sqrt-gap of ~0.02 here, so it "verified EXACT" while gating
    // nothing near a real, discrete DT error.)
    let voxel = dt.voxel_size();
    for x in dt_min[0]..dt_min[0] + i32::try_from(dt_dims[0]).expect("fits") {
        for y in dt_min[1]..dt_min[1] + i32::try_from(dt_dims[1]).expect("fits") {
            for z in dt_min[2]..dt_min[2] + i32::try_from(dt_dims[2]).expect("fits") {
                let brute = seeds
                    .iter()
                    .map(|s| {
                        let d = [s[0] - x, s[1] - y, s[2] - z];
                        f64::from(d[0] * d[0] + d[1] * d[1] + d[2] * d[2])
                    })
                    .fold(f64::INFINITY, f64::min);
                let got = dt.distance([x, y, z]).expect("in box");
                let brute_dist = brute.sqrt() * voxel;
                assert_eq!(
                    got.to_bits(),
                    brute_dist.to_bits(),
                    "DT not bit-exact at ({x},{y},{z}): {got} vs {brute_dist} (sq={brute})"
                );
            }
        }
    }
    // Triangle inequality (1-Lipschitz in the voxel metric).
    let mut s2 = 0xD7_0002u64;
    for _ in 0..200 {
        let p = [
            dt_min[0] + (lcg(&mut s2) * dt_dims[0] as f64) as i32,
            dt_min[1] + (lcg(&mut s2) * dt_dims[1] as f64) as i32,
            dt_min[2] + (lcg(&mut s2) * dt_dims[2] as f64) as i32,
        ];
        let q = [p[0] + 1, p[1], p[2]];
        if let (Some(dp), Some(dq)) = (dt.distance(p), dt.distance(q)) {
            assert!(
                (dp - dq).abs() <= 1.0 + 1e-9,
                "1-Lipschitz violated: |{dp} - {dq}| > 1"
            );
        }
    }
    verdict(
        "rv-002",
        "DT exact vs O(n^2) on scattered+slab fixture; 1-Lipschitz",
    );
}

#[test]
fn rv_003_point_cloud_queries_and_normals() {
    let mut seed = 0xC1_0003u64;
    // Random cloud for query correctness.
    let pts: Vec<[f64; 3]> = (0..400)
        .map(|_| {
            [
                lcg(&mut seed) * 4.0,
                lcg(&mut seed) * 4.0,
                lcg(&mut seed) * 4.0,
            ]
        })
        .collect();
    let cloud = PointCloud::new(pts.clone(), 0.5).expect("cloud");
    let dist2 = |a: [f64; 3], b: [f64; 3]| {
        (a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)
    };
    for _ in 0..24 {
        let q = [
            lcg(&mut seed) * 4.0,
            lcg(&mut seed) * 4.0,
            lcg(&mut seed) * 4.0,
        ];
        // Radius query vs brute force.
        let r = 0.4 + lcg(&mut seed) * 0.8;
        let mut brute: Vec<usize> = (0..pts.len())
            .filter(|&i| dist2(pts[i], q) <= r * r)
            .collect();
        let mut got = cloud.radius_query(q, r);
        brute.sort_unstable();
        got.sort_unstable();
        assert_eq!(got, brute, "radius query must match brute force");
        // kNN vs brute force (as a set of distances — ties may reorder).
        let k = 7;
        let knn = cloud.knn(q, k);
        assert_eq!(knn.len(), k);
        let mut brute_d: Vec<f64> = (0..pts.len()).map(|i| dist2(pts[i], q)).collect();
        brute_d.sort_by(f64::total_cmp);
        let worst_knn = knn.iter().map(|&i| dist2(pts[i], q)).fold(0.0, f64::max);
        assert!(
            worst_knn <= brute_d[k - 1] + 1e-12,
            "kNN must return the k nearest"
        );
    }
    // Normals on a sampled sphere: PCA + propagation gives outward normals.
    let mut sphere_pts = Vec::new();
    let n_lat = 24;
    for i in 1..n_lat {
        let theta = core::f64::consts::PI * f64::from(i) / f64::from(n_lat);
        let n_lon = 2 * n_lat;
        for j in 0..n_lon {
            let phi = 2.0 * core::f64::consts::PI * f64::from(j) / f64::from(n_lon);
            sphere_pts.push([
                theta.sin() * phi.cos(),
                theta.sin() * phi.sin(),
                theta.cos(),
            ]);
        }
    }
    let mut sphere = PointCloud::new(sphere_pts.clone(), 0.3).expect("sphere cloud");
    sphere.estimate_normals(8).expect("normals");
    let normals = sphere.normals.as_ref().expect("present");
    let mut aligned = 0usize;
    for (p, n) in sphere_pts.iter().zip(normals) {
        let radial = p; // unit sphere: outward normal == position
        let cos = p[0] * n[0] + p[1] * n[1] + p[2] * n[2];
        let _ = radial;
        if cos > 0.9 {
            aligned += 1;
        }
    }
    let frac = aligned as f64 / sphere_pts.len() as f64;
    if frac <= 0.97 {
        let mut hist = [0usize; 5]; // cos buckets: <-0.9, -0.9..0, 0..0.9, >0.9, other
        for (p, n) in sphere_pts.iter().zip(normals) {
            let cos = p[0] * n[0] + p[1] * n[1] + p[2] * n[2];
            let b: usize = match cos {
                c if c > 0.9 => 3,
                c if c > 0.0 => 2,
                c if c > -0.9 => 1,
                _ => 0,
            };
            hist[b] += 1;
        }
        panic!("normals must align outward: frac {frac:.3}, hist {hist:?}");
    }
    // Degenerate queries refuse structurally.
    assert!(PointCloud::new(Vec::new(), 0.5).is_err());
    assert!(
        PointCloud::new(vec![[0.0; 3]; 4], 0.5)
            .expect("tiny")
            .clone()
            .estimate_normals(8)
            .is_err()
    );
    verdict(
        "rv-003",
        "radius/kNN match brute force; sphere normals >97% outward-aligned",
    );
}

/// A small ground-structure fixture (tetrahedral cell).
fn tetra_fixture() -> (Vec<LatticeNode>, Vec<Strut>) {
    let nodes = vec![
        LatticeNode {
            pos: [0.0, 0.0, 0.0],
            radius: 0.02,
        },
        LatticeNode {
            pos: [1.0, 0.0, 0.0],
            radius: 0.02,
        },
        LatticeNode {
            pos: [0.5, 0.9, 0.0],
            radius: 0.03,
        },
        LatticeNode {
            pos: [0.5, 0.3, 0.8],
            radius: 0.02,
        },
    ];
    let struts = vec![
        Strut {
            a: 0,
            b: 1,
            radius: 0.05,
        },
        Strut {
            a: 0,
            b: 2,
            radius: 0.04,
        },
        Strut {
            a: 1,
            b: 2,
            radius: 0.04,
        },
        Strut {
            a: 0,
            b: 3,
            radius: 0.03,
        },
        Strut {
            a: 1,
            b: 3,
            radius: 0.03,
        },
        Strut {
            a: 2,
            b: 3,
            radius: 0.03,
        },
    ];
    (nodes, struts)
}

#[test]
fn rv_004_lattice_round_trip_degenerates_and_realization() {
    let (nodes, struts) = tetra_fixture();
    let lattice = LatticeGraph::new(nodes.clone(), struts.clone()).expect("valid");
    // FrankenNetworkx round-trip preserves attributes exactly.
    let g = lattice.to_fnx();
    let back = LatticeGraph::from_fnx(&g).expect("round trip");
    assert_eq!(back.nodes, lattice.nodes, "node attributes preserved");
    let mut orig_struts = lattice.struts.clone();
    let mut back_struts = back.struts.clone();
    let key = |s: &Strut| (s.a.min(s.b), s.a.max(s.b));
    orig_struts.sort_by_key(&key);
    back_struts.sort_by_key(&key);
    assert_eq!(back_struts.len(), orig_struts.len());
    for (a, b) in orig_struts.iter().zip(&back_struts) {
        assert_eq!(key(a), key(b));
        assert!(
            (a.radius - b.radius).abs() < 1e-15,
            "strut radius preserved"
        );
    }
    // Degenerates refuse with the offending element named.
    let coincident = LatticeGraph::new(
        vec![
            LatticeNode {
                pos: [0.0; 3],
                radius: 0.01,
            },
            LatticeNode {
                pos: [0.0; 3],
                radius: 0.01,
            },
        ],
        Vec::new(),
    );
    assert!(
        matches!(coincident, Err(fs_rep_voxel::VoxelError::Lattice { ref what }) if what.contains("coincident"))
    );
    let zero_len = LatticeGraph::new(
        vec![LatticeNode {
            pos: [0.0; 3],
            radius: 0.01,
        }],
        vec![Strut {
            a: 0,
            b: 0,
            radius: 0.01,
        }],
    );
    assert!(
        matches!(zero_len, Err(fs_rep_voxel::VoxelError::Lattice { ref what }) if what.contains("zero-length"))
    );
    // Realization: level-set behavior. Strut midpoints are inside…
    for s in &struts {
        let (a, b) = (nodes[s.a].pos, nodes[s.b].pos);
        let mid = [
            f64::midpoint(a[0], b[0]),
            f64::midpoint(a[1], b[1]),
            f64::midpoint(a[2], b[2]),
        ];
        assert!(lattice.sdf(mid) < 0.0, "strut midpoint must be inside");
    }
    // …far points are outside, and a probe crossing one strut sees
    // exactly one sign change in and one out (a closed surface).
    assert!(lattice.sdf([10.0, 10.0, 10.0]) > 1.0);
    let mut signs = Vec::new();
    let steps = 400;
    for k in 0..=steps {
        let t = f64::from(k) / f64::from(steps);
        let p = [0.5, -0.5 + t, 0.0]; // crosses the 0-1 and 0-2/1-2 plane region
        signs.push(lattice.sdf(p) < 0.0);
    }
    let transitions = signs.windows(2).filter(|w| w[0] != w[1]).count();
    assert!(
        transitions >= 2 && transitions % 2 == 0,
        "a straight probe must enter and exit a watertight solid an even number of \
         times (got {transitions})"
    );
    for line in lattice.realization_receipts() {
        println!("{line}");
    }
    verdict(
        "rv-004",
        "fnx round-trip exact; coincident/zero-length refuse; realization level set \
         closed along probes",
    );
}

#[test]
fn rv_005_occupancy_chart_contract() {
    use asupersync::types::Budget;
    use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
    use fs_geom::{Differentiability, Point3};
    let field = ball(6); // voxel_size 0.1, centered at origin voxel
    let chart = OccupancyChart::try_new(field, TEST_DT_BUDGET).expect("admissible chart");
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 1,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        // Center is inside with negative distance; far point positive.
        let center = chart.eval(Point3::new(0.05, 0.05, 0.05), &cx);
        assert!(center.signed_distance < 0.0, "ball center is inside");
        let far = chart.eval(Point3::new(5.0, 0.0, 0.0), &cx);
        assert!(far.signed_distance > 0.0, "far point is outside");
        // The analytic ball radius is ~0.65 (6.5 voxels x 0.1): the chart
        // distance at a far probe matches |p| - R within a voxel diagonal.
        let probe = Point3::new(2.0, 0.0, 0.0);
        let s = chart.eval(probe, &cx);
        let analytic = 2.05 - 0.65; // voxel-center probe to ball surface
        assert!(
            (s.signed_distance - analytic).abs() < 0.2,
            "chart distance {} vs analytic {analytic}",
            s.signed_distance
        );
        // Center-metric error is estimate-grade: the piecewise-constant
        // representative has no continuity or Lipschitz authority.
        assert_eq!(s.error.kind, fs_evidence::NumericalKind::Estimate);
        assert!(s.lipschitz.is_none());
        assert_eq!(chart.differentiability(), Differentiability::Unknown);
        // Support box contains the ball.
        let support = chart.support();
        assert!(chart.inside(Point3::new(0.05, 0.05, 0.05), &cx));
        let _ = support;
        for invalid in [
            Point3::new(f64::NAN, 0.0, 0.0),
            Point3::new(0.0, f64::INFINITY, 0.0),
            Point3::new(0.0, 0.0, f64::MAX),
        ] {
            let sample = chart.eval(invalid, &cx);
            assert!(sample.signed_distance.is_nan());
            assert!(sample.gradient.is_none());
            assert!(sample.lipschitz.is_none());
            assert_eq!(sample.error.kind, fs_evidence::NumericalKind::NoClaim);
        }

        // G0 regression for the authority bug: the abstract solid is the
        // closed unit cube. Floor semantics choose the empty cell at x = 1,
        // whose center is one unit from the occupied center. The returned
        // representative is therefore +1 even though the true signed
        // distance on the cube face is exactly zero. A half-diagonal band
        // falsely excluded zero; the chart must not label any such band a
        // rigorous enclosure, continuous field, or Lipschitz field.
        let mut one_voxel = OccupancyField::new(1.0, [0.0; 3]).expect("one-voxel frame");
        one_voxel.set([0, 0, 0]);
        let one_voxel = OccupancyChart::try_new(one_voxel, 27).expect("one-voxel chart");
        let face = one_voxel.eval(Point3::new(1.0, 0.5, 0.5), &cx);
        assert_eq!(face.signed_distance.to_bits(), 1.0f64.to_bits());
        assert_eq!(face.error.kind, fs_evidence::NumericalKind::Estimate);
        assert!(
            face.error.lo <= 0.0 && face.error.hi >= 0.0,
            "resolution estimate must include the analytic face distance: {:?}",
            face.error
        );
        assert!(face.lipschitz.is_none());
        assert_eq!(one_voxel.differentiability(), Differentiability::Unknown);

        let just_inside_x = f64::from_bits(1.0f64.to_bits() - 1);
        let just_inside = one_voxel.eval(Point3::new(just_inside_x, 0.5, 0.5), &cx);
        assert_eq!(just_inside.signed_distance.to_bits(), (-1.0f64).to_bits());
        assert_eq!(
            (face.signed_distance - just_inside.signed_distance).to_bits(),
            2.0f64.to_bits(),
            "center metric must remain explicitly discontinuous at the voxel face"
        );
    });
    verdict(
        "rv-005",
        "occupancy chart: center metric is estimate-only, with no false continuity, \
         Lipschitz, or enclosure authority",
    );
}

fn assert_frame_ops_refuse_without_mutation(
    base: &OccupancyField,
    other: &OccupancyField,
    mismatch: &str,
) {
    let before = active_set(base);
    let mut candidate = base.clone();
    assert!(matches!(
        candidate.union(other),
        Err(VoxelError::FrameMismatch { .. })
    ));
    assert_eq!(
        active_set(&candidate),
        before,
        "union mutated on {mismatch} mismatch"
    );

    let mut candidate = base.clone();
    assert!(matches!(
        candidate.intersect(other),
        Err(VoxelError::FrameMismatch { .. })
    ));
    assert_eq!(
        active_set(&candidate),
        before,
        "intersection mutated on {mismatch} mismatch"
    );

    let mut candidate = base.clone();
    assert!(matches!(
        candidate.subtract(other),
        Err(VoxelError::FrameMismatch { .. })
    ));
    assert_eq!(
        active_set(&candidate),
        before,
        "subtraction mutated on {mismatch} mismatch"
    );
}

#[test]
fn rv_006_frames_and_world_coordinates_fail_closed() {
    let mut base = OccupancyField::new(1.0, [0.0; 3]).expect("base frame");
    base.set([1, 2, 3]);
    base.set([4, 5, 6]);
    assert_eq!(base.voxel_size().to_bits(), 1.0f64.to_bits());
    assert_eq!(base.origin().map(f64::to_bits), [0.0f64.to_bits(); 3]);
    let conversion = OccupancyField::new(2.0, [10.0, 20.0, 30.0]).expect("conversion frame");
    assert_eq!(
        conversion.voxel_of([9.999, 20.0, 34.0]),
        Ok([-1, 0, 2]),
        "finite coordinate conversion keeps floor boundary semantics"
    );
    for point in [
        [f64::NAN, 0.0, 0.0],
        [0.0, f64::NEG_INFINITY, 0.0],
        [0.0, 0.0, f64::MAX],
    ] {
        assert!(matches!(
            conversion.voxel_of(point),
            Err(VoxelError::WorldCoordinateOutOfRange { .. })
        ));
    }

    for (mismatch, other) in [
        (
            "voxel size",
            OccupancyField::new(0.5, [0.0; 3]).expect("different size"),
        ),
        (
            "origin",
            OccupancyField::new(1.0, [1.0, 0.0, 0.0]).expect("different origin"),
        ),
    ] {
        assert_frame_ops_refuse_without_mutation(&base, &other, mismatch);
    }

    assert!(matches!(
        OccupancyField::new(1.0, [f64::NAN, 0.0, 0.0]),
        Err(VoxelError::Parameters { .. })
    ));
    assert!(matches!(
        OccupancyField::new(1.0, [0.0, f64::INFINITY, 0.0]),
        Err(VoxelError::Parameters { .. })
    ));
    let empty = OccupancyField::new(1.0, [0.0; 3]).expect("empty frame");
    assert!(matches!(
        OccupancyChart::try_new(empty, TEST_DT_BUDGET),
        Err(VoxelError::EmptyOccupancy {
            operation: "occupancy chart construction"
        })
    ));

    verdict(
        "rv-006",
        "frames are immutable; booleans preserve receivers on mismatch; world coordinate \
         conversion and empty charts fail closed",
    );
}

#[test]
fn rv_007_dense_work_fails_closed() {
    // The full i32 span must be computed in i64/u128, then rejected by
    // the explicit budget rather than overflowing signed subtraction.
    let mut extrema = OccupancyField::new(1.0, [0.0; 3]).expect("extrema frame");
    extrema.set([i32::MIN, 0, 0]);
    extrema.set([i32::MAX, 0, 0]);
    assert!(matches!(
        fs_rep_voxel::euclidean_dt(&extrema, 1_024),
        Err(VoxelError::VoxelBudgetExceeded {
            required: 4_294_967_296,
            maximum: 1_024,
            ..
        })
    ));

    let mut volume = OccupancyField::new(1.0, [0.0; 3]).expect("volume frame");
    volume.set([0, 0, 0]);
    volume.set([10, 10, 10]);
    assert!(matches!(
        fs_rep_voxel::euclidean_dt(&volume, 1_000),
        Err(VoxelError::VoxelBudgetExceeded {
            required: 1_331,
            maximum: 1_000,
            ..
        })
    ));

    // The explicit allocation budget is not permission to leave the
    // numerical range where integer squared distances are exact in f64.
    let mut too_wide = OccupancyField::new(1.0, [0.0; 3]).expect("wide frame");
    too_wide.set([0, 0, 0]);
    too_wide.set([1 << 26, 0, 0]);
    assert!(matches!(
        fs_rep_voxel::euclidean_dt(&too_wide, (1 << 26) + 1),
        Err(VoxelError::ExactnessRangeExceeded {
            max_squared_distance: 4_503_599_627_370_496,
            maximum: 4_503_599_627_370_495,
            ..
        })
    ));

    // A single active cell needs a 3^3 complement scan. The constructor
    // must enforce that budget before populating the complement field.
    let mut one = OccupancyField::new(1.0, [0.0; 3]).expect("one-cell frame");
    one.set([0, 0, 0]);
    assert!(matches!(
        OccupancyChart::try_new(one, 26),
        Err(VoxelError::VoxelBudgetExceeded {
            required: 27,
            maximum: 26,
            ..
        })
    ));

    for boundary in [i32::MIN, i32::MAX] {
        let mut edge = OccupancyField::new(1.0, [0.0; 3]).expect("edge frame");
        edge.set([boundary, 0, 0]);
        assert!(matches!(
            OccupancyChart::try_new(edge, 1_000),
            Err(VoxelError::CoordinateRange {
                operation: "occupancy complement halo",
                axis: 0,
                halo: 1,
                ..
            })
        ));
    }

    verdict(
        "rv-007",
        "DT and complement halo refuse coordinate overflow, excess volume, and numeric \
         inexactness before dense work",
    );
}

#[test]
fn rv_008_chart_support_is_exact_voxel_union() {
    // Chart support is the union of voxel cubes, not centers padded by a
    // full edge length. This fixture has exactly representable bounds.
    let mut support_field = OccupancyField::new(2.0, [10.0, 20.0, 30.0]).expect("support frame");
    support_field.set([-1, 2, 0]);
    support_field.set([1, -2, 3]);
    let support = OccupancyChart::try_new(support_field, TEST_DT_BUDGET)
        .expect("support chart")
        .support();
    assert_eq!(
        [support.min.x, support.min.y, support.min.z].map(f64::to_bits),
        [8.0, 16.0, 30.0].map(f64::to_bits)
    );
    assert_eq!(
        [support.max.x, support.max.y, support.max.z].map(f64::to_bits),
        [14.0, 26.0, 38.0].map(f64::to_bits)
    );

    verdict(
        "rv-008",
        "chart support is the exact union of occupied voxel cubes",
    );
}
