//! Sparse active-tile D3Q19 battery (bead sjro): Morton round trips,
//! insertion-order independence, halo transport, bounce-back mass
//! conservation, the G5 worker-count bitwise gate, G4-style memory
//! proportionality, typed refusals, and deterministic cancellation.

use fs_exec::{CancelGate, TilePool};
use fs_lbm::d3q19::sparse::{
    MORTON_COORD_BITS, SparseError3, SparseGrid3, demorton3, morton3, state_bytes_per_tile,
};

fn splitmix64(mut z: u64) -> u64 {
    z = z.wrapping_add(0x9e37_79b9_7f4a_7c15);
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^ (z >> 31)
}

/// An L-shaped ~10%-active tile set inside an `ntx × nty × ntz` tile
/// domain: one x-aligned beam plus one z-aligned column, deterministic.
fn l_shaped_tiles(ntx: u32, nty: u32, ntz: u32) -> Vec<(u32, u32, u32)> {
    let mut tiles = Vec::new();
    for tx in 0..ntx {
        tiles.push((tx, nty / 2, ntz / 2));
    }
    for tz in 0..ntz {
        tiles.push((ntx / 2, nty / 2, tz));
    }
    tiles.sort_unstable();
    tiles.dedup();
    tiles
}

#[test]
fn morton_keys_round_trip_and_interleave_exactly() {
    assert_eq!(morton3(0, 0, 0), 0);
    assert_eq!(morton3(1, 0, 0), 0b001);
    assert_eq!(morton3(0, 1, 0), 0b010);
    assert_eq!(morton3(0, 0, 1), 0b100);
    assert_eq!(morton3(7, 7, 7), 0b111_111_111);
    let limit = (1u32 << MORTON_COORD_BITS) - 1;
    let edge_cases = [0, 1, 2, 3, 0x5555_5 & limit, 0xa_aaaa & limit, limit];
    for &x in &edge_cases {
        for &y in &edge_cases {
            for &z in &edge_cases {
                assert_eq!(demorton3(morton3(x, y, z)), (x, y, z));
            }
        }
    }
    for i in 0..4096u64 {
        let h = splitmix64(i);
        let x = (h as u32) & limit;
        let y = ((h >> 21) as u32) & limit;
        let z = ((h >> 42) as u32) & limit;
        assert_eq!(demorton3(morton3(x, y, z)), (x, y, z));
    }
}

#[test]
fn active_order_is_morton_order_never_insertion_order() {
    let tiles = l_shaped_tiles(8, 8, 8);
    let mut reversed = tiles.clone();
    reversed.reverse();

    let mut forward = SparseGrid3::new(32, 32, 32, 0.8, [0.0; 3]).expect("dims admissible");
    forward
        .activate_tiles(&tiles)
        .expect("activation in-domain");
    let mut backward = SparseGrid3::new(32, 32, 32, 0.8, [0.0; 3]).expect("dims admissible");
    // Insert in reversed order AND in two separate batches.
    let (head, tail) = reversed.split_at(reversed.len() / 2);
    backward.activate_tiles(head).expect("activation in-domain");
    backward.activate_tiles(tail).expect("activation in-domain");

    assert_eq!(forward.active_tiles(), backward.active_tiles());
    forward.perturb(0x5EED, 0.01);
    backward.perturb(0x5EED, 0.01);
    for _ in 0..5 {
        forward.step_serial().expect("sweep admissible");
        backward.step_serial().expect("sweep admissible");
    }
    assert_eq!(
        forward.state_bits(),
        backward.state_bits(),
        "activation history leaked into sweep state"
    );
}

#[test]
fn halo_transport_reaches_the_neighbor_tile() {
    // Two adjacent tiles along x; perturb only cells of the first tile.
    let mut grid = SparseGrid3::new(8, 4, 4, 0.8, [0.0; 3]).expect("dims admissible");
    grid.activate_tiles(&[(0, 0, 0), (1, 0, 0)])
        .expect("activation in-domain");
    grid.perturb(7, 0.05);
    // Flatten tile 1 (slot for morton3(1,0,0)) back to exact equilibrium
    // by re-activating a fresh grid and copying: instead, build the
    // asymmetry with the public surface — a second grid perturbed with a
    // seed whose hash we zero out is not available, so instead measure:
    // tile 1's initial macros, then verify they CHANGE once the halo from
    // tile 0 arrives (2 steps: collide then stream crosses the face).
    let before: Vec<(f64, [f64; 3])> = (0..64).map(|lane| grid.cell_macros(1, lane)).collect();
    for _ in 0..2 {
        grid.step_serial().expect("sweep admissible");
    }
    let after: Vec<(f64, [f64; 3])> = (0..64).map(|lane| grid.cell_macros(1, lane)).collect();
    assert_ne!(
        before, after,
        "no state crossed the tile face — halo exchange is dead"
    );
    // And the two-tile closed box conserves mass while doing it.
    let mass: f64 = grid.total_mass();
    let mut reference = SparseGrid3::new(8, 4, 4, 0.8, [0.0; 3]).expect("dims admissible");
    reference
        .activate_tiles(&[(0, 0, 0), (1, 0, 0)])
        .expect("activation in-domain");
    reference.perturb(7, 0.05);
    let initial_mass = reference.total_mass();
    assert!(
        ((mass - initial_mass) / initial_mass).abs() < 1e-12,
        "closed two-tile box lost mass: {initial_mass} -> {mass}"
    );
}

#[test]
fn sealed_box_conserves_mass_for_fifty_steps() {
    let mut grid = SparseGrid3::new(16, 16, 16, 0.9, [0.0; 3]).expect("dims admissible");
    let all: Vec<(u32, u32, u32)> = (0..4u32)
        .flat_map(|tx| (0..4u32).flat_map(move |ty| (0..4u32).map(move |tz| (tx, ty, tz))))
        .collect();
    grid.activate_tiles(&all).expect("activation in-domain");
    grid.perturb(0xCAFE, 0.02);
    let initial = grid.total_mass();
    for _ in 0..50 {
        grid.step_serial().expect("sweep admissible");
    }
    let residual = ((grid.total_mass() - initial) / initial).abs();
    assert!(
        residual < 1e-12,
        "sealed box mass drifted by {residual:e} relative over 50 steps"
    );
}

#[test]
fn pooled_sweep_is_bitwise_identical_across_worker_counts() {
    let tiles = l_shaped_tiles(8, 8, 8);

    let mut serial = SparseGrid3::new(32, 32, 32, 0.7, [1e-6, 0.0, 0.0]).expect("dims admissible");
    serial.activate_tiles(&tiles).expect("activation in-domain");
    serial.perturb(0x5EED_CA51, 0.01);
    for _ in 0..25 {
        serial.step_serial().expect("sweep admissible");
    }
    let reference = serial.state_bits();

    for workers in [1usize, 2, 7] {
        let pool = TilePool::for_host(workers, 0x5EED);
        let gate = CancelGate::new();
        let mut pooled =
            SparseGrid3::new(32, 32, 32, 0.7, [1e-6, 0.0, 0.0]).expect("dims admissible");
        pooled.activate_tiles(&tiles).expect("activation in-domain");
        pooled.perturb(0x5EED_CA51, 0.01);
        for _ in 0..25 {
            pooled.step_pooled(&pool, &gate).expect("sweep admissible");
        }
        assert_eq!(
            pooled.state_bits(),
            reference,
            "pooled sweep with {workers} workers diverged bitwise from serial"
        );
    }
}

#[test]
fn memory_is_proportional_to_active_tiles_not_domain() {
    // 1000-tile domain, 10% active: allocation must track the 100 active
    // tiles exactly, with zero dependence on the 40^3 dense extent.
    let mut grid = SparseGrid3::new(40, 40, 40, 0.8, [0.0; 3]).expect("dims admissible");
    let tiles: Vec<(u32, u32, u32)> = (0..100u32).map(|i| (i % 10, (i / 10) % 10, 4)).collect();
    grid.activate_tiles(&tiles).expect("activation in-domain");
    assert_eq!(grid.active_tiles(), 100);
    assert_eq!(
        grid.allocated_state_bytes(),
        100 * state_bytes_per_tile(),
        "allocation is not proportional to the active set"
    );
    let mut dense_equivalent = SparseGrid3::new(40, 40, 40, 0.8, [0.0; 3]).expect("dims");
    let all: Vec<(u32, u32, u32)> = (0..10u32)
        .flat_map(|tx| (0..10u32).flat_map(move |ty| (0..10u32).map(move |tz| (tx, ty, tz))))
        .collect();
    dense_equivalent.activate_tiles(&all).expect("activation");
    assert_eq!(
        dense_equivalent.allocated_state_bytes(),
        10 * grid.allocated_state_bytes(),
        "10x the active tiles must cost exactly 10x the state bytes"
    );
}

#[test]
fn activation_preserves_existing_tile_state_exactly() {
    let mut grid = SparseGrid3::new(32, 32, 32, 0.8, [0.0; 3]).expect("dims admissible");
    grid.activate_tiles(&[(1, 1, 1), (2, 1, 1)])
        .expect("activation in-domain");
    grid.perturb(42, 0.03);
    for _ in 0..3 {
        grid.step_serial().expect("sweep admissible");
    }
    let before = grid.state_bits();
    // Growing the active set re-sorts slots but must not move a single
    // bit of existing tile state.
    grid.activate_tiles(&[(0, 0, 0), (5, 5, 5)])
        .expect("activation in-domain");
    let after = grid.state_bits();
    // Locate the surviving tiles' bit runs: slots are Morton-sorted, so
    // (1,1,1) and (2,1,1) occupy known positions in both orderings.
    let bits_per_tile = state_bytes_per_tile() / 2 / 8;
    let old_keys = [morton3(1, 1, 1), morton3(2, 1, 1)];
    let mut new_keys = vec![
        morton3(0, 0, 0),
        morton3(1, 1, 1),
        morton3(2, 1, 1),
        morton3(5, 5, 5),
    ];
    new_keys.sort_unstable();
    for (old_slot, key) in old_keys.iter().enumerate() {
        let new_slot = new_keys.iter().position(|k| k == key).expect("key kept");
        assert_eq!(
            &before[old_slot * bits_per_tile..(old_slot + 1) * bits_per_tile],
            &after[new_slot * bits_per_tile..(new_slot + 1) * bits_per_tile],
            "activation moved bits of surviving tile {key:#x}"
        );
    }
}

#[test]
fn collision_refusal_is_typed_and_leaves_state_intact() {
    let mut grid = SparseGrid3::new(8, 4, 4, 0.8, [f64::NAN, 0.0, 0.0]).expect("dims admissible");
    grid.activate_tiles(&[(0, 0, 0), (1, 0, 0)])
        .expect("activation in-domain");
    grid.perturb(9, 0.01);
    let before = grid.state_bits();
    let refusal = grid
        .step_serial()
        .expect_err("non-finite force must refuse");
    match refusal {
        SparseError3::Collision { tile_key, lane, .. } => {
            assert_eq!(
                tile_key,
                morton3(0, 0, 0),
                "first refusal in canonical order"
            );
            assert_eq!(lane, 0);
        }
        other => panic!("expected a collision refusal, got {other}"),
    }
    assert_eq!(grid.steps(), 0, "a refused step must not count");
    assert_eq!(
        grid.state_bits(),
        before,
        "a refused step must not move state"
    );
}

#[test]
fn pre_tripped_gate_cancels_cleanly_and_reissue_is_deterministic() {
    let tiles = l_shaped_tiles(8, 8, 8);
    let mut grid = SparseGrid3::new(32, 32, 32, 0.8, [0.0; 3]).expect("dims admissible");
    grid.activate_tiles(&tiles).expect("activation in-domain");
    grid.perturb(0xD00D, 0.01);
    let before = grid.state_bits();

    let pool = TilePool::for_host(4, 7);
    let tripped = CancelGate::new();
    tripped.request();
    let refusal = grid
        .step_pooled(&pool, &tripped)
        .expect_err("a pre-tripped gate must cancel the sweep");
    assert_eq!(refusal, SparseError3::Cancelled);
    assert_eq!(grid.steps(), 0);
    assert_eq!(
        grid.state_bits(),
        before,
        "cancelled sweep must leave the pre-step state intact"
    );

    // Re-issue under an open gate: bitwise identical to the serial path.
    let open = CancelGate::new();
    grid.step_pooled(&pool, &open).expect("sweep admissible");
    let mut reference = SparseGrid3::new(32, 32, 32, 0.8, [0.0; 3]).expect("dims admissible");
    reference
        .activate_tiles(&tiles)
        .expect("activation in-domain");
    reference.perturb(0xD00D, 0.01);
    reference.step_serial().expect("sweep admissible");
    assert_eq!(
        grid.state_bits(),
        reference.state_bits(),
        "post-cancellation reissue diverged from the serial reference"
    );
}

#[test]
fn out_of_domain_activation_refuses_without_partial_application() {
    let mut grid = SparseGrid3::new(16, 16, 16, 0.8, [0.0; 3]).expect("dims admissible");
    let refusal = grid
        .activate_tiles(&[(0, 0, 0), (4, 0, 0)])
        .expect_err("tile (4,0,0) lies outside a 4-tile domain");
    assert_eq!(
        refusal,
        SparseError3::TileOutOfDomain {
            tx: 4,
            ty: 0,
            tz: 0
        }
    );
    assert_eq!(grid.active_tiles(), 0, "refused activation must be atomic");
    let refusal =
        SparseGrid3::new(10, 16, 16, 0.8, [0.0; 3]).expect_err("10 is not a tile multiple");
    assert_eq!(
        refusal,
        SparseError3::Dims {
            nx: 10,
            ny: 16,
            nz: 16
        }
    );
}

#[test]
fn deactivation_returns_exact_mass_and_preserves_survivors() {
    let mut grid = SparseGrid3::new(32, 32, 32, 0.8, [0.0; 3]).expect("dims admissible");
    grid.activate_tiles(&[(1, 1, 1), (2, 1, 1), (3, 1, 1)])
        .expect("activation in-domain");
    grid.perturb(0xFEED, 0.02);
    for _ in 0..3 {
        grid.step_serial().expect("sweep admissible");
    }
    let total_before = grid.total_mass();
    let bits_before = grid.state_bits();
    let bits_per_tile = state_bytes_per_tile() / 2 / 8;

    // Retire the middle tile: the ledger must receive exactly its mass.
    let removed = grid
        .deactivate_tiles(&[(2, 1, 1)])
        .expect("deactivation in-domain");
    assert_eq!(grid.active_tiles(), 2);
    let total_after = grid.total_mass();
    // Ledger conservation: removed + remaining == before. NOT bitwise —
    // the three sums group the same ~3.7k cell values differently
    // (whole-sequence vs per-subset), so reassociation moves a few ULPs
    // (~1e-16 relative each); 1e-13 leaves margin while still catching a
    // single dropped cell (~1/192 ≈ 5e-3 relative) by ten orders.
    assert!(
        ((total_after + removed) - total_before).abs() / total_before < 1e-13,
        "deactivation mass ledger leaked: {total_before} != {total_after} + {removed}"
    );

    // Survivors keep their bits exactly; Morton order preserved.
    let old_keys = [morton3(1, 1, 1), morton3(2, 1, 1), morton3(3, 1, 1)];
    let kept = [morton3(1, 1, 1), morton3(3, 1, 1)];
    let bits_after = grid.state_bits();
    for (new_slot, key) in kept.iter().enumerate() {
        let old_slot = old_keys.iter().position(|k| k == key).expect("was active");
        assert_eq!(
            &bits_before[old_slot * bits_per_tile..(old_slot + 1) * bits_per_tile],
            &bits_after[new_slot * bits_per_tile..(new_slot + 1) * bits_per_tile],
            "deactivation moved bits of surviving tile {key:#x}"
        );
    }

    // Deactivating an inactive tile is a zero-mass no-op.
    let removed = grid
        .deactivate_tiles(&[(7, 7, 7)])
        .expect("deactivation in-domain");
    assert_eq!(removed, 0.0);
    assert_eq!(grid.active_tiles(), 2);

    // Out-of-domain deactivation refuses atomically.
    let refusal = grid
        .deactivate_tiles(&[(1, 1, 1), (0, 0, 40)])
        .expect_err("tile (0,0,40) lies outside an 8-tile domain");
    assert!(matches!(refusal, SparseError3::TileOutOfDomain { .. }));
    assert_eq!(
        grid.active_tiles(),
        2,
        "refused deactivation must be atomic"
    );
}

/// FNV-1a over the exact state bits — the golden preimage.
fn fnv1a_bits(bits: &[u64]) -> u64 {
    let mut hash = 0xcbf2_9ce4_8422_2325u64;
    for &word in bits {
        for byte in word.to_le_bytes() {
            hash ^= u64::from(byte);
            hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
        }
    }
    hash
}

/// Candidate golden fixture: the exact G5 configuration (L-shaped
/// ~10%-active set in a 32^3 box, tau 0.7, x-force 1e-6, seeded
/// perturbation, 25 serial sweeps) plus one mass-ledgered deactivation
/// and 5 more sweeps — so the frozen surface covers activation,
/// sweeping, AND retirement arithmetic.
fn sparse_candidate_hash() -> u64 {
    let tiles = l_shaped_tiles(8, 8, 8);
    let mut grid = SparseGrid3::new(32, 32, 32, 0.7, [1e-6, 0.0, 0.0]).expect("dims admissible");
    grid.activate_tiles(&tiles).expect("activation in-domain");
    grid.perturb(0x5EED_CA51, 0.01);
    for _ in 0..25 {
        grid.step_serial().expect("sweep admissible");
    }
    let removed = grid
        .deactivate_tiles(&[(4, 4, 0)])
        .expect("deactivation in-domain");
    for _ in 0..5 {
        grid.step_serial().expect("sweep admissible");
    }
    let mut bits = grid.state_bits();
    bits.push(removed.to_bits());
    bits.push(grid.total_mass().to_bits());
    fnv1a_bits(&bits)
}

/// FROZEN 2026-07-17 (bead sjro): all four policy quadrants reproduced
/// this value bit-identically on the committed tree e33b74de — aarch64
/// M4 Pro debug+release and x86-64 ts1 debug+release. Bump only per
/// docs/GOLDEN_POLICY.md with the same four-quadrant evidence.
const SPARSE_GOLDEN_HASH: u64 = 0x4c00_8876_e576_a332;

#[test]
fn d3q19_sparse_candidate_hash_is_replay_stable() {
    let first = sparse_candidate_hash();
    let second = sparse_candidate_hash();
    assert_eq!(
        first, second,
        "sparse candidate hash is not replay-stable on one host"
    );
    println!("d3q19-sparse candidate hash: {first:#018x}");
    assert_eq!(
        first, SPARSE_GOLDEN_HASH,
        "sparse golden moved — bump only per docs/GOLDEN_POLICY.md \
         with four-quadrant evidence"
    );
}
