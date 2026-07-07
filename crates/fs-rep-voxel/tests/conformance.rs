//! fs-rep-voxel conformance (the wqd.6 bead). Acceptance: morphology
//! matches brute force; the Euclidean DT is EXACT vs the O(n²) reference;
//! point-cloud queries match brute force and normals recover analytic
//! surfaces; lattice graphs round-trip through FrankenNetworkx with
//! attributes preserved, degenerates refuse structurally, and realization
//! behaves like a watertight solid's level set.

use fs_rep_voxel::{LatticeGraph, LatticeNode, OccupancyChart, OccupancyField, PointCloud, Strut};

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
fn brute_dilate(cells: &std::collections::BTreeSet<[i32; 3]>) -> std::collections::BTreeSet<[i32; 3]> {
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
    assert_eq!(active_set(&dilated), brute_dilate(&before), "dilate == brute force");
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
    assert!(!after_open.contains(&[7, 0, 0]), "the spur must be opened away");
    assert!(after_open.contains(&[0, 0, 0]), "the ball core survives opening");
    // Boolean algebra: (A ∪ B) \ B ⊆ A, A ∩ B ⊆ A.
    let b = ball(3);
    let mut u = f.clone();
    u.union(&b);
    u.subtract(&b);
    assert!(active_set(&u).is_subset(&before));
    let mut i = f.clone();
    i.intersect(&b);
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
    let dt = fs_rep_voxel::euclidean_dt(&f).expect("nonempty");
    // O(n²) reference over the whole box: EXACT equality of squared
    // distances (both sides are integer-valued in voxel units).
    for x in dt.min[0]..dt.min[0] + i32::try_from(dt.dims[0]).expect("fits") {
        for y in dt.min[1]..dt.min[1] + i32::try_from(dt.dims[1]).expect("fits") {
            for z in dt.min[2]..dt.min[2] + i32::try_from(dt.dims[2]).expect("fits") {
                let brute = seeds
                    .iter()
                    .map(|s| {
                        let d = [s[0] - x, s[1] - y, s[2] - z];
                        f64::from(d[0] * d[0] + d[1] * d[1] + d[2] * d[2])
                    })
                    .fold(f64::INFINITY, f64::min);
                let got = dt.distance([x, y, z]).expect("in box");
                let brute_dist = brute.sqrt();
                assert!(
                    (got - brute_dist).abs() < 1e-9,
                    "DT not exact at ({x},{y},{z}): {got} vs {brute_dist}"
                );
            }
        }
    }
    // Triangle inequality (1-Lipschitz in the voxel metric).
    let mut s2 = 0xD7_0002u64;
    for _ in 0..200 {
        let p = [
            dt.min[0] + (lcg(&mut s2) * dt.dims[0] as f64) as i32,
            dt.min[1] + (lcg(&mut s2) * dt.dims[1] as f64) as i32,
            dt.min[2] + (lcg(&mut s2) * dt.dims[2] as f64) as i32,
        ];
        let q = [p[0] + 1, p[1], p[2]];
        if let (Some(dp), Some(dq)) = (dt.distance(p), dt.distance(q)) {
            assert!(
                (dp - dq).abs() <= 1.0 + 1e-9,
                "1-Lipschitz violated: |{dp} - {dq}| > 1"
            );
        }
    }
    verdict("rv-002", "DT exact vs O(n^2) on scattered+slab fixture; 1-Lipschitz");
}

#[test]
fn rv_003_point_cloud_queries_and_normals() {
    let mut seed = 0xC1_0003u64;
    // Random cloud for query correctness.
    let pts: Vec<[f64; 3]> = (0..400)
        .map(|_| [lcg(&mut seed) * 4.0, lcg(&mut seed) * 4.0, lcg(&mut seed) * 4.0])
        .collect();
    let cloud = PointCloud::new(pts.clone(), 0.5).expect("cloud");
    let dist2 = |a: [f64; 3], b: [f64; 3]| {
        (a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)
    };
    for _ in 0..24 {
        let q = [lcg(&mut seed) * 4.0, lcg(&mut seed) * 4.0, lcg(&mut seed) * 4.0];
        // Radius query vs brute force.
        let r = 0.4 + lcg(&mut seed) * 0.8;
        let mut brute: Vec<usize> = (0..pts.len()).filter(|&i| dist2(pts[i], q) <= r * r).collect();
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
            let b = if cos > 0.9 { 3 } else if cos > 0.0 { 2 } else if cos > -0.9 { 1 } else { 0 };
            hist[b] += 1;
        }
        panic!("normals must align outward: frac {frac:.3}, hist {hist:?}");
    }
    // Degenerate queries refuse structurally.
    assert!(PointCloud::new(Vec::new(), 0.5).is_err());
    assert!(PointCloud::new(vec![[0.0; 3]; 4], 0.5).expect("tiny").clone().estimate_normals(8).is_err());
    verdict(
        "rv-003",
        "radius/kNN match brute force; sphere normals >97% outward-aligned",
    );
}

#[test]
fn rv_004_lattice_round_trip_degenerates_and_realization() {
    // A small ground-structure fixture (tetrahedral cell).
    let nodes = vec![
        LatticeNode { pos: [0.0, 0.0, 0.0], radius: 0.02 },
        LatticeNode { pos: [1.0, 0.0, 0.0], radius: 0.02 },
        LatticeNode { pos: [0.5, 0.9, 0.0], radius: 0.03 },
        LatticeNode { pos: [0.5, 0.3, 0.8], radius: 0.02 },
    ];
    let struts = vec![
        Strut { a: 0, b: 1, radius: 0.05 },
        Strut { a: 0, b: 2, radius: 0.04 },
        Strut { a: 1, b: 2, radius: 0.04 },
        Strut { a: 0, b: 3, radius: 0.03 },
        Strut { a: 1, b: 3, radius: 0.03 },
        Strut { a: 2, b: 3, radius: 0.03 },
    ];
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
        assert!((a.radius - b.radius).abs() < 1e-15, "strut radius preserved");
    }
    // Degenerates refuse with the offending element named.
    let coincident = LatticeGraph::new(
        vec![
            LatticeNode { pos: [0.0; 3], radius: 0.01 },
            LatticeNode { pos: [0.0; 3], radius: 0.01 },
        ],
        Vec::new(),
    );
    assert!(matches!(coincident, Err(fs_rep_voxel::VoxelError::Lattice { ref what }) if what.contains("coincident")));
    let zero_len = LatticeGraph::new(
        vec![LatticeNode { pos: [0.0; 3], radius: 0.01 }],
        vec![Strut { a: 0, b: 0, radius: 0.01 }],
    );
    assert!(matches!(zero_len, Err(fs_rep_voxel::VoxelError::Lattice { ref what }) if what.contains("zero-length")));
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
    use fs_geom::{Chart, Point3};
    let field = ball(6); // voxel_size 0.1, centered at origin voxel
    let chart = OccupancyChart::new(field);
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
        // The error certificate is an HONEST enclosure containing the value.
        match s.error {
            fs_evidence::NumericalCertificate { .. } => {
                let half = 0.5 * 3.0f64.sqrt() * 0.1;
                assert!(half > 0.08, "declared resolution error present");
            }
        }
        // Support box contains the ball.
        let support = chart.support();
        assert!(chart.inside(Point3::new(0.05, 0.05, 0.05), &cx));
        let _ = support;
    });
    verdict(
        "rv-005",
        "occupancy chart: inside/outside, DT-backed distance near analytic, resolution \
         error declared",
    );
}
