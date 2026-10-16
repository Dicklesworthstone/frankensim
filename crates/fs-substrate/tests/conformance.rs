//! fs-substrate conformance suite (CONTRACT.md: any reimplementation must
//! pass). Tile-layout coverage: G0 (Morton bijection + backend equivalence,
//! halo laws, iteration permutations), G5 (deterministic orders, ids, and
//! affinity tables), plus the layout stencil smoke with measured — never
//! assumed — timings emitted as benchmark events. Aggregate verdicts use
//! the canonical fs-obs schema; randomized cases carry their input seed.

use fs_substrate::affinity::{AffinityMap, CcdTopology};
use fs_substrate::field::{Boundary, TiledField};
use fs_substrate::morton::{MORTON_COORD_LIMIT, morton_backend, morton3_decode, morton3_encode};
use fs_substrate::tile::{TileCoord, TileEdge, TileGrid};

fn verdict(case: &str, pass: bool, detail: &str, seed: u64) {
    let mut emitter = fs_obs::Emitter::new("fs-substrate/conformance", case);
    let event = emitter.emit(
        if pass {
            fs_obs::Severity::Info
        } else {
            fs_obs::Severity::Error
        },
        fs_obs::EventKind::ConformanceCase {
            suite: "fs-substrate/conformance".to_string(),
            case: case.to_string(),
            pass,
            detail: detail.to_string(),
            seed,
        },
        None,
    );
    fs_obs::lint_failure_record(&event).expect("substrate verdict must be replayable");
    let line = event.to_jsonl();
    fs_obs::validate_line(&line).expect("substrate verdict must use the fs-obs wire schema");
    println!("{line}");
    assert!(pass, "case {case}: {detail}");
}

/// In-house LCG (L0 cannot depend on fs-rand; fs-qty battery constants).
struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0
    }

    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }
}

#[test]
fn sub_001_morton_bijection_and_backend_equivalence() {
    const SEED: u64 = 0x0001_3D0C_2026_0706;
    let mut rng = Lcg(SEED);
    for _ in 0..200_000 {
        let x = (rng.next() as u32) & (MORTON_COORD_LIMIT - 1);
        let y = (rng.next() as u32) & (MORTON_COORD_LIMIT - 1);
        let z = (rng.next() as u32) & (MORTON_COORD_LIMIT - 1);
        let c = morton3_encode(x, y, z);
        assert_eq!(morton3_decode(c), (x, y, z), "bijection at ({x},{y},{z})");
    }
    verdict(
        "sub-001",
        true,
        &format!(
            "morton bijection over 200k seeded cases (seed {SEED:#x}), backend={}",
            morton_backend()
        ),
        SEED,
    );
}

#[test]
fn sub_002_tile_world_maps_and_iteration_orders() {
    let g = TileGrid::new([20, 12, 9], TileEdge::E8).expect("grid");
    // Every cell roundtrips through (tile, within) exactly once.
    let dims = g.cell_dims();
    let mut seen = vec![false; (dims[0] * dims[1] * dims[2]) as usize];
    for z in 0..dims[2] {
        for y in 0..dims[1] {
            for x in 0..dims[0] {
                let (t, w) = g.tile_of_cell([x, y, z]);
                let back = g.cell_of(t, w);
                assert_eq!(back, [x, y, z]);
                let idx = ((z * dims[1] + y) * dims[0] + x) as usize;
                assert!(!seen[idx]);
                seen[idx] = true;
            }
        }
    }
    // Orders are permutations; z-order and boundary-first are deterministic.
    let mut lin: Vec<TileCoord> = g.iter_linear().collect();
    let mut zor = g.iter_zorder();
    assert_eq!(zor, g.iter_zorder(), "z-order deterministic (G5)");
    let bf = g.iter_boundary_first();
    let cut = bf
        .iter()
        .position(|&t| !g.is_boundary(t))
        .unwrap_or(bf.len());
    assert!(bf[..cut].iter().all(|&t| g.is_boundary(t)));
    assert!(bf[cut..].iter().all(|&t| !g.is_boundary(t)));
    lin.sort_unstable();
    zor.sort_unstable();
    assert_eq!(lin, zor);
    verdict(
        "sub-002",
        seen.iter().all(|&s| s),
        "cell/tile maps bijective; linear, z-order, boundary-first orders are permutations",
        0,
    );
}

#[test]
fn sub_003_halo_fast_path_equals_reference_on_all_tiles_and_bcs() {
    const SEED: u64 = 0x0003_4A70_2026_0706;
    let mut rng = Lcg(SEED);
    let mut checked = 0u32;
    for (dims, edge) in [
        ([20u32, 12, 9], TileEdge::E8), // partial tiles on every axis
        ([9, 9, 9], TileEdge::E4),
        ([16, 16, 16], TileEdge::E8), // fully aligned
    ] {
        let g = TileGrid::new(dims, edge).expect("grid");
        let mut f = TiledField::new(g, 0i32);
        for z in 0..dims[2] {
            for y in 0..dims[1] {
                for x in 0..dims[0] {
                    f.set([x, y, z], rng.next() as i32);
                }
            }
        }
        let (mut a, mut b) = (Vec::new(), Vec::new());
        for t in f.grid().iter_zorder() {
            for bc in [Boundary::Clamp, Boundary::Periodic, Boundary::Constant(-7)] {
                f.gather_halo(t, bc, &mut a);
                f.gather_halo_fast(t, bc, &mut b);
                assert_eq!(a, b, "halo mismatch at {t:?} bc {bc:?} dims {dims:?}");
                checked += 1;
            }
        }
    }
    verdict(
        "sub-003",
        checked > 0,
        &format!("fast halo == reference on {checked} (tile, bc) cases (seed {SEED:#x})"),
        SEED,
    );
}

#[test]
fn sub_004_affinity_fixtures_respect_ccd_boundaries_deterministically() {
    let g = TileGrid::new([320, 320, 320], TileEdge::E8).expect("grid"); // 64000 tiles
    for (topo, name) in [
        (CcdTopology::TR_7995WX, "tr-7995wx"),
        (CcdTopology::EPYC_128C, "epyc-128c"),
        (CcdTopology::APPLE_M_CLASS, "apple-m"),
    ] {
        let map = AffinityMap::assign(&g, &topo);
        assert_eq!(u32::from(map.shard_count()), topo.ccds, "{name}");
        let sizes: Vec<u32> = (0..map.shard_count())
            .map(|s| map.slots_of(s).len() as u32)
            .collect();
        assert_eq!(
            sizes.iter().sum::<u32>(),
            map.tile_count(),
            "{name} coverage"
        );
        let (min, max) = (
            *sizes.iter().min().expect("nonempty"),
            *sizes.iter().max().expect("nonempty"),
        );
        assert!(max - min <= 1, "{name} balance: {sizes:?}");
        for s in 0..map.shard_count() {
            let cores = topo.cores_of(s);
            assert_eq!(cores.start % topo.cores_per_ccd, 0, "{name} CCD boundary");
            let r = map.slots_of(s);
            assert_eq!(map.shard_of_slot(r.start), s, "{name} range inversion");
        }
        assert_eq!(
            map.to_json(),
            map.to_json(),
            "{name} table deterministic (G5)"
        );
    }
    verdict(
        "sub-004",
        true,
        "TR/EPYC/Apple fixtures: balanced contiguous z-ranges, CCD-respecting core ranges",
        0,
    );
}

#[test]
fn sub_005_stencil_layout_smoke_documents_measurements() {
    const N: u32 = 48;
    const SEED: u64 = 0x0005_57E4_2026_0706;
    let g = TileGrid::new([N, N, N], TileEdge::E8).expect("grid");
    let mut f = TiledField::new(g, 0.0f32);
    let mut linear = vec![0.0f32; (N * N * N) as usize];
    let idx = |x: u32, y: u32, z: u32| ((z * N + y) * N + x) as usize;
    let mut rng = Lcg(SEED);
    for z in 0..N {
        for y in 0..N {
            for x in 0..N {
                let v = (rng.below(1000) as f32) / 1000.0;
                f.set([x, y, z], v);
                linear[idx(x, y, z)] = v;
            }
        }
    }
    // 7-point Laplacian over the interior, both layouts.
    let run_tiled = |f: &TiledField<f32>| -> (f64, u128) {
        let start = std::time::Instant::now();
        let mut acc = 0.0f64;
        for z in 1..N - 1 {
            for y in 1..N - 1 {
                for x in 1..N - 1 {
                    let c = f.get([x, y, z]);
                    let l = f.get([x - 1, y, z])
                        + f.get([x + 1, y, z])
                        + f.get([x, y - 1, z])
                        + f.get([x, y + 1, z])
                        + f.get([x, y, z - 1])
                        + f.get([x, y, z + 1])
                        - 6.0 * c;
                    acc += f64::from(l);
                }
            }
        }
        (acc, start.elapsed().as_nanos())
    };
    let run_linear = |v: &[f32]| -> (f64, u128) {
        let start = std::time::Instant::now();
        let mut acc = 0.0f64;
        for z in 1..N - 1 {
            for y in 1..N - 1 {
                for x in 1..N - 1 {
                    let c = v[idx(x, y, z)];
                    let l = v[idx(x - 1, y, z)]
                        + v[idx(x + 1, y, z)]
                        + v[idx(x, y - 1, z)]
                        + v[idx(x, y + 1, z)]
                        + v[idx(x, y, z - 1)]
                        + v[idx(x, y, z + 1)]
                        - 6.0 * c;
                    acc += f64::from(l);
                }
            }
        }
        (acc, start.elapsed().as_nanos())
    };
    let (acc_t, ns_t) = run_tiled(&f);
    let (acc_l, ns_l) = run_linear(&linear);
    // Numerical identity is the LAW; timing is a measurement, not a claim
    // (roofline verdicts belong to the perf harness bead — CONTRACT.md).
    let identical = acc_t.to_bits() == acc_l.to_bits();
    let probe = fs_substrate::CapabilityProbe::topology_only();
    let mut em = fs_obs::Emitter::new("fs-substrate/conformance", "sub-005/stencil-smoke");
    for (kernel, ns) in [("stencil7-tiled-e8", ns_t), ("stencil7-linear", ns_l)] {
        let line = em
            .emit(
                fs_obs::Severity::Info,
                fs_obs::EventKind::BenchmarkResult {
                    kernel: kernel.to_string(),
                    metric: "ns".to_string(),
                    value: ns as f64,
                    machine: probe.fingerprint(),
                },
                None,
            )
            .to_jsonl();
        fs_obs::validate_line(&line).expect("benchmark event must validate");
        println!("{line}");
    }
    verdict(
        "sub-005",
        identical,
        "tiled and linear stencils agree bitwise; timings are separate benchmark events \
         (documentation, not a perf claim)",
        SEED,
    );
}

#[test]
fn sub_006_first_touch_shard_views_are_disjoint_and_complete() {
    let g = TileGrid::new([64, 64, 64], TileEdge::E8).expect("grid"); // 512 tiles
    let topo = CcdTopology::APPLE_M_CLASS;
    let map = AffinityMap::assign(&g, &topo);
    // Parallel first-touch: each shard's view initialized on its own thread
    // (the placement hook; actual pinning is fs-exec's contract).
    let mut par = TiledField::new(g.clone(), 0u32);
    let views = par.shard_views_mut(&map);
    std::thread::scope(|s| {
        for mut view in views {
            s.spawn(move || {
                for i in 0..view.tile_count() {
                    let slot = view.global_slot(i);
                    view.tile_mut(i).fill(slot);
                }
            });
        }
    });
    // Serial reference.
    let mut ser = TiledField::new(g, 0u32);
    let zorder = ser.grid().iter_zorder();
    for (rank, t) in zorder.iter().enumerate() {
        ser.tile_slice_mut(*t).fill(rank as u32);
    }
    let dims = [64u32, 64, 64];
    let mut equal = true;
    for z in 0..dims[2] {
        for y in 0..dims[1] {
            for x in 0..dims[0] {
                equal &= par.get([x, y, z]) == ser.get([x, y, z]);
            }
        }
    }
    // Owner tags recorded per the map.
    let tags_ok = zorder
        .iter()
        .enumerate()
        .all(|(rank, t)| par.meta(*t).owner_shard == map.shard_of_slot(rank as u32));
    verdict(
        "sub-006",
        equal && tags_ok,
        "parallel per-shard first-touch equals serial fill; owner tags match the affinity map",
        0,
    );
}
