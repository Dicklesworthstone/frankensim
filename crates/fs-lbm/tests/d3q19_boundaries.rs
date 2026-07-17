//! D3Q19 boundary battery (bead 40p2): deterministic per-tile link masks,
//! voxelized solids, regularized velocity/pressure faces, stationary and
//! moving halfway bounce-back, closed-domain leak, lid-cavity circulation,
//! pressure-driven duct accuracy, and a boundary-surface golden.

use fs_lbm::d3q19::{CollisionModel3, D3Q19_BOUNDARY_BIT_SEMANTICS_VERSION, E3, Q3, TILE};
use fs_lbm::{
    BoundaryGrid3, BoundarySpec3, CS2, Face3, FaceBoundary3, duct_analytic, equilibrium3,
};

fn verdict(name: &str, pass: bool, details: &str) {
    println!("{{\"test\":\"{name}\",\"pass\":{pass},\"details\":\"{details}\"}}");
    assert!(pass, "{name}: {details}");
}

fn mask(directions: &[usize]) -> u32 {
    directions
        .iter()
        .fold(0u32, |bits, direction| bits | (1u32 << direction))
}

fn nonequilibrium_stress(populations: [f64; Q3], rho: f64, velocity: [f64; 3]) -> [[f64; 3]; 3] {
    let equilibrium = equilibrium3(rho, velocity);
    let mut stress = [[0.0; 3]; 3];
    for q in 0..Q3 {
        let e = [f64::from(E3[q].0), f64::from(E3[q].1), f64::from(E3[q].2)];
        let nonequilibrium = populations[q] - equilibrium[q];
        for row in 0..3 {
            for column in 0..3 {
                stress[row][column] += e[row] * e[column] * nonequilibrium;
            }
        }
    }
    stress
}

/// lbm3bc-001: wall and SDF-solid links are represented by exact masks in
/// tile-major order, and enumeration is stable by `(tile,lane,direction)`.
#[test]
fn link_masks_cover_planar_and_voxelized_boundaries() {
    assert_eq!(D3Q19_BOUNDARY_BIT_SEMANTICS_VERSION, 1);
    let grid = BoundaryGrid3::new(8, 8, 8, 0.8, [0.0; 3], BoundarySpec3::periodic_duct_z());
    // Hand-enumerated from the pinned table; this oracle does not reuse the
    // production predicate that constructs the masks.
    let positive_x = mask(&[1, 7, 9, 11, 13]);
    let negative_x = mask(&[2, 8, 10, 12, 14]);
    let positive_y = mask(&[3, 7, 10, 15, 17]);
    let negative_y = mask(&[4, 8, 9, 16, 18]);
    let positive_z = mask(&[5, 11, 14, 15, 18]);
    let negative_z = mask(&[6, 12, 13, 16, 17]);
    assert_eq!(grid.link_mask(0, 3, 3), positive_x);
    assert_eq!(grid.link_mask(7, 3, 3), negative_x);
    assert_eq!(grid.link_mask(3, 0, 3), positive_y);
    assert_eq!(grid.link_mask(3, 7, 3), negative_y);
    assert_eq!(grid.link_mask(0, 0, 3), positive_x | positive_y);
    assert_eq!(grid.link_mask(3, 3, 3), 0);
    let cavity = BoundaryGrid3::new(
        8,
        8,
        8,
        0.8,
        [0.0; 3],
        BoundarySpec3::lid_cavity([0.02, 0.0, 0.0]),
    );
    assert_eq!(cavity.link_mask(3, 3, 0), positive_z);
    assert_eq!(cavity.link_mask(3, 3, 7), negative_z);
    assert_eq!(grid.tile_link_masks().len(), 8);
    assert!(
        grid.tile_link_masks()
            .iter()
            .all(|tile| std::ptr::from_ref(tile).addr().is_multiple_of(128))
    );

    let links = grid.boundary_links();
    assert!(links.windows(2).all(|pair| {
        (pair[0].tile, pair[0].lane, pair[0].direction)
            < (pair[1].tile, pair[1].lane, pair[1].direction)
    }));

    let mut voxel = BoundaryGrid3::new(8, 8, 8, 0.8, [0.0; 3], BoundarySpec3::periodic());
    // One complete x=3 slab, sampled at cell centers.
    voxel.voxelize_sdf(|point| (point[0] - 3.5).abs() - 0.25);
    assert!(voxel.is_solid(3, 4, 4));
    assert!(!voxel.is_solid(4, 4, 4));
    assert_eq!(voxel.link_mask(4, 4, 4), positive_x);

    let mut isolated = BoundaryGrid3::new(8, 8, 8, 0.8, [0.0; 3], BoundarySpec3::periodic());
    isolated.voxelize_sdf(|point| {
        let dx = point[0] - 3.5;
        let dy = point[1] - 3.5;
        let dz = point[2] - 3.5;
        dx.mul_add(dx, dy.mul_add(dy, dz * dz)) - 0.01
    });
    let mut actual = isolated
        .boundary_links()
        .into_iter()
        .map(|link| (link.cell, link.direction))
        .collect::<Vec<_>>();
    let mut expected = vec![
        ([4, 3, 3], 1),
        ([2, 3, 3], 2),
        ([3, 4, 3], 3),
        ([3, 2, 3], 4),
        ([3, 3, 4], 5),
        ([3, 3, 2], 6),
        ([4, 4, 3], 7),
        ([2, 2, 3], 8),
        ([4, 2, 3], 9),
        ([2, 4, 3], 10),
        ([4, 3, 4], 11),
        ([2, 3, 2], 12),
        ([4, 3, 2], 13),
        ([2, 3, 4], 14),
        ([3, 4, 4], 15),
        ([3, 2, 2], 16),
        ([3, 4, 2], 17),
        ([3, 2, 4], 18),
    ];
    actual.sort_unstable();
    expected.sort_unstable();
    assert_eq!(actual, expected);
    verdict(
        "lbm3bc-001-link-masks",
        true,
        &format!(
            "{} planar links; voxel slab mask {positive_x:#07x}",
            links.len()
        ),
    );
}

/// lbm3bc-001b (G0/G3): the general grid exposes collision selection without
/// changing the legacy BGK constructor. Equal-rate central moments follow the
/// BGK macroscopic trajectory to solve roundoff, while a split-rate reduced
/// cumulant grid remains conservative and refuses unverified forcing.
#[test]
fn collision_model_selection_is_backward_compatible_and_conservative() {
    let tau = 0.81;
    let rate = 1.0 / tau;
    let central_model = CollisionModel3::CentralMoment {
        second_order_rate: rate,
        higher_order_rate: rate,
    };
    let mut bgk = BoundaryGrid3::new(8, 8, 8, tau, [0.0; 3], BoundarySpec3::periodic());
    let mut central = BoundaryGrid3::with_collision_model(
        8,
        8,
        8,
        central_model,
        [0.0; 3],
        BoundarySpec3::periodic(),
    );
    assert_eq!(bgk.collision_model(), CollisionModel3::Bgk { tau });
    assert_eq!(central.collision_model(), central_model);
    assert_eq!(bgk.viscosity().to_bits(), ((tau - 0.5) / 3.0).to_bits());
    assert!((central.viscosity() - bgk.viscosity()).abs() < f64::EPSILON);

    bgk.perturb(0xC011_1DE3, 1e-3);
    central.perturb(0xC011_1DE3, 1e-3);
    bgk.run(2);
    central.run(2);
    let mut max_macro_delta = 0.0_f64;
    for z in 0..8 {
        for y in 0..8 {
            for x in 0..8 {
                max_macro_delta =
                    max_macro_delta.max((bgk.density(x, y, z) - central.density(x, y, z)).abs());
                let bgk_velocity = bgk.velocity(x, y, z);
                let central_velocity = central.velocity(x, y, z);
                for axis in 0..3 {
                    max_macro_delta =
                        max_macro_delta.max((bgk_velocity[axis] - central_velocity[axis]).abs());
                }
            }
        }
    }
    assert!(
        max_macro_delta < 5e-12,
        "equal-rate central/BGK grid trajectories differ by {max_macro_delta:.3e}"
    );

    let reduced_model = CollisionModel3::ReducedCumulant {
        second_order_rate: rate,
        third_order_rate: 1.35,
        fourth_order_rate: 1.72,
    };
    let mut reduced = BoundaryGrid3::with_collision_model(
        8,
        8,
        8,
        reduced_model,
        [0.0; 3],
        BoundarySpec3::periodic(),
    );
    reduced.perturb(0xC011_1DE3, 1e-3);
    let mass_before = reduced.total_mass();
    reduced.run(4);
    assert!((reduced.total_mass() - mass_before).abs() < 1e-11);

    let forced_moment_grid = std::panic::catch_unwind(|| {
        BoundaryGrid3::with_collision_model(
            8,
            8,
            8,
            central_model,
            [1e-7, 0.0, 0.0],
            BoundarySpec3::periodic(),
        )
    });
    assert!(forced_moment_grid.is_err());
}

/// lbm3bc-002: regularized faces impose their declared macroscopic moments;
/// pressure and density remain finite during a velocity-inlet/pressure-outlet
/// transient.
#[test]
fn regularized_velocity_and_pressure_faces_impose_moments() {
    let inlet = [0.0, 0.0, 0.02];
    let mut grid = BoundaryGrid3::new(
        8,
        8,
        16,
        0.8,
        [0.0; 3],
        BoundarySpec3::velocity_pressure_duct_z(inlet, 1.0),
    );
    grid.run(200);

    let positive_z = mask(&[5, 11, 14, 15, 18]);
    assert_eq!(grid.open_link_mask(2, 3, 0), positive_z);
    // At the inlet/x-min/y-min rim, links whose pull crosses only the open
    // face stay open; diagonals that also cross a side wall are wall-owned.
    let open_only_positive_z = mask(&[5, 14, 18]);
    assert_eq!(grid.open_link_mask(0, 0, 0), open_only_positive_z);
    assert_eq!(
        grid.link_mask(0, 0, 0) & positive_z,
        positive_z & !open_only_positive_z
    );

    let mut max_inlet_error = 0.0_f64;
    let mut max_outlet_density_error = 0.0_f64;
    for y in 1..7 {
        for x in 1..7 {
            let u = grid.velocity(x, y, 0);
            max_inlet_error = max_inlet_error.max(
                u.iter()
                    .zip(inlet)
                    .map(|(actual, target)| (actual - target).abs())
                    .fold(0.0, f64::max),
            );
            max_outlet_density_error =
                max_outlet_density_error.max((grid.density(x, y, 15) - 1.0).abs());
        }
    }
    let finite = (0..16).all(|z| {
        (0..8).all(|y| {
            (0..8).all(|x| {
                let rho = grid.density(x, y, z);
                let u = grid.velocity(x, y, z);
                rho.is_finite() && rho > 0.0 && u.iter().all(|value| value.is_finite())
            })
        })
    });
    let boundary_rho = grid.density(2, 3, 0);
    let boundary_velocity = grid.velocity(2, 3, 0);
    let interior_rho = grid.density(2, 3, 1);
    let interior_velocity = grid.velocity(2, 3, 1);
    let inlet_density_error = (boundary_rho - interior_rho).abs();
    let boundary_stress =
        nonequilibrium_stress(grid.populations(2, 3, 0), boundary_rho, boundary_velocity);
    let interior_stress =
        nonequilibrium_stress(grid.populations(2, 3, 1), interior_rho, interior_velocity);
    let stress_norm = interior_stress
        .iter()
        .flatten()
        .map(|value| value.abs())
        .fold(0.0, f64::max);
    let stress_error = boundary_stress
        .iter()
        .flatten()
        .zip(interior_stress.iter().flatten())
        .map(|(boundary, interior)| (boundary - interior).abs())
        .fold(0.0, f64::max);

    let outlet_rho = grid.density(2, 3, 15);
    let outlet_velocity = grid.velocity(2, 3, 15);
    let outlet_interior_rho = grid.density(2, 3, 14);
    let outlet_interior_velocity = grid.velocity(2, 3, 14);
    let outlet_velocity_error = outlet_velocity
        .iter()
        .zip(outlet_interior_velocity)
        .map(|(boundary, interior)| (boundary - interior).abs())
        .fold(0.0, f64::max);
    let outlet_stress =
        nonequilibrium_stress(grid.populations(2, 3, 15), outlet_rho, outlet_velocity);
    let outlet_interior_stress = nonequilibrium_stress(
        grid.populations(2, 3, 14),
        outlet_interior_rho,
        outlet_interior_velocity,
    );
    let outlet_stress_error = outlet_stress
        .iter()
        .flatten()
        .zip(outlet_interior_stress.iter().flatten())
        .map(|(boundary, interior)| (boundary - interior).abs())
        .fold(0.0, f64::max);
    verdict(
        "lbm3bc-002-regularized-faces",
        finite
            && max_inlet_error < 2e-12
            && max_outlet_density_error < 2e-12
            && inlet_density_error < 2e-12
            && outlet_velocity_error < 2e-12
            && stress_norm > 1e-8
            && stress_error < 2e-12
            && outlet_stress_error < 2e-12,
        &format!(
            "inlet u/rho err {max_inlet_error:.3e}/{inlet_density_error:.3e}; outlet rho/u err {max_outlet_density_error:.3e}/{outlet_velocity_error:.3e}; stress norm {stress_norm:.3e}; inlet/outlet stress err {stress_error:.3e}/{outlet_stress_error:.3e}"
        ),
    );
}

#[test]
fn regularized_faces_are_axis_generic() {
    let wall = FaceBoundary3::stationary_wall();
    let target_x = [0.01, 0.002, -0.003];
    let target_y = [-0.002, 0.01, 0.003];
    let cases = [
        (
            BoundarySpec3::new([
                FaceBoundary3::Velocity { velocity: target_x },
                FaceBoundary3::Pressure { density: 1.0 },
                wall,
                wall,
                wall,
                wall,
            ]),
            [0, 3, 3],
            [7, 3, 3],
            target_x,
        ),
        (
            BoundarySpec3::new([
                wall,
                wall,
                FaceBoundary3::Velocity { velocity: target_y },
                FaceBoundary3::Pressure { density: 1.0 },
                wall,
                wall,
            ]),
            [3, 0, 3],
            [3, 7, 3],
            target_y,
        ),
    ];
    for (boundaries, inlet, outlet, target) in cases {
        let mut grid = BoundaryGrid3::new(8, 8, 8, 0.8, [0.0; 3], boundaries);
        grid.step();
        let actual = grid.velocity(inlet[0], inlet[1], inlet[2]);
        assert!(
            actual
                .iter()
                .zip(target)
                .all(|(actual, target)| (*actual - target).abs() < 2e-15)
        );
        assert!((grid.density(outlet[0], outlet[1], outlet[2]) - 1.0).abs() < 2e-15);
    }
}

/// lbm3bc-003: a closed no-slip duct has no systematic leak. The default
/// fixture is short; the release acceptance below carries the full 10k-step
/// requirement.
#[test]
fn stationary_halfway_bounce_back_has_no_systematic_leak() {
    let mut grid = BoundaryGrid3::new(
        8,
        8,
        8,
        0.8,
        [0.0, 0.0, 1e-6],
        BoundarySpec3::periodic_duct_z(),
    );
    grid.perturb(0x40_02, 1e-3);
    let initial = grid.total_mass();
    grid.run(1_000);
    let relative_leak = (grid.total_mass() - initial).abs() / initial;

    let mut voxel = BoundaryGrid3::new(8, 8, 8, 0.8, [0.0; 3], BoundarySpec3::periodic());
    voxel.voxelize_sdf(|point| {
        let dx = point[0] - 4.0;
        let dy = point[1] - 4.0;
        let dz = point[2] - 4.0;
        dx.mul_add(dx, dy.mul_add(dy, dz * dz)) - 2.0
    });
    voxel.perturb(0x5D_F0, 1e-3);
    let voxel_initial = voxel.total_mass();
    voxel.run(1_000);
    let voxel_leak = (voxel.total_mass() - voxel_initial).abs() / voxel_initial;
    verdict(
        "lbm3bc-003-closed-leak",
        relative_leak < 1e-11 && voxel_leak < 1e-11,
        &format!(
            "planar leak {relative_leak:.3e}; voxel-obstacle leak {voxel_leak:.3e} over 1000 steps"
        ),
    );
}

/// lbm3bc-004: the moving-lid sign is correct and the primary mid-plane
/// circulation is clockwise for a +x lid.
#[test]
fn moving_lid_one_step_matches_the_link_momentum_oracle() {
    let lid_velocity = 0.04;
    let mut grid = BoundaryGrid3::new(
        8,
        8,
        8,
        0.8,
        [0.0; 3],
        BoundarySpec3::lid_cavity([lid_velocity, 0.0, 0.0]),
    );
    grid.step();
    // The two x-bearing incoming lid diagonals each contribute u_lid/6
    // with opposite population signs, hence total x momentum u_lid/3.
    let observed = grid.velocity(4, 7, 4)[0];
    let expected = lid_velocity / 3.0;
    verdict(
        "lbm3bc-004a-lid-one-step",
        (observed - expected).abs() < 2e-15,
        &format!("observed {observed:.16e}; expected {expected:.16e}"),
    );
}

#[test]
fn moving_wall_open_rim_keeps_wall_links_stationary() {
    let wall = FaceBoundary3::stationary_wall();
    let lid_velocity = 0.04;
    let boundaries = BoundarySpec3::new([
        wall,
        wall,
        wall,
        FaceBoundary3::Wall {
            velocity: [lid_velocity, 0.0, 0.0],
        },
        FaceBoundary3::Velocity { velocity: [0.0; 3] },
        FaceBoundary3::Pressure { density: 1.0 },
    ]);
    let mut grid = BoundaryGrid3::new(8, 8, 8, 0.8, [0.0; 3], boundaries);
    grid.step();

    let rim_ux = grid.velocity(4, 7, 0)[0];
    let lid_ux = grid.velocity(4, 7, 4)[0];
    verdict(
        "lbm3bc-004b-moving-open-rim",
        rim_ux.abs() < 2e-15 && (lid_ux - lid_velocity / 3.0).abs() < 2e-15,
        &format!("rim ux {rim_ux:.3e}; unobstructed lid ux {lid_ux:.16e}"),
    );
}

#[test]
fn moving_lid_drives_the_expected_primary_vortex() {
    let mut grid = BoundaryGrid3::new(
        8,
        8,
        8,
        0.8,
        [0.0; 3],
        BoundarySpec3::lid_cavity([0.04, 0.0, 0.0]),
    );
    let initial = grid.total_mass();
    grid.run(1_500);

    let top_ux = grid.velocity(4, 7, 4)[0];
    let lower_ux = grid.velocity(4, 1, 4)[0];
    let dvdx = 0.5 * (grid.velocity(5, 4, 4)[1] - grid.velocity(3, 4, 4)[1]);
    let dudy = 0.5 * (grid.velocity(4, 5, 4)[0] - grid.velocity(4, 3, 4)[0]);
    let omega_z = dvdx - dudy;
    let relative_leak = (grid.total_mass() - initial).abs() / initial;
    verdict(
        "lbm3bc-004-lid-vortex",
        top_ux > 0.0 && lower_ux < 0.0 && omega_z < 0.0 && relative_leak < 1e-11,
        &format!(
            "top ux {top_ux:.4e}; lower ux {lower_ux:.4e}; omega_z {omega_z:.4e}; leak {relative_leak:.2e}"
        ),
    );
}

/// A small pressure-driven duct catches sign or reconstruction regressions in
/// the normal test lane. The full 3% gate is the ignored release fixture.
#[test]
fn pressure_driven_duct_has_the_poiseuille_shape() {
    let (rho_in, rho_out) = (1.000_05, 0.999_95);
    let mut grid = BoundaryGrid3::new(
        8,
        8,
        16,
        0.933_012_701_892_219_3,
        [0.0; 3],
        BoundarySpec3::pressure_duct_z(rho_in, rho_out),
    );
    grid.run(4_000);
    let section = grid.velocity_section_z(8, 2);
    let center = section[4 * 8 + 4];
    let near_wall = section[4 * 8];
    let symmetry = (section[3 * 8 + 2] - section[4 * 8 + 5]).abs() / center.abs();
    verdict(
        "lbm3bc-005-poiseuille-shape",
        center > 0.0 && near_wall > 0.0 && center > 2.0 * near_wall && symmetry < 2e-4,
        &format!("center {center:.4e}; near wall {near_wall:.4e}; symmetry {symmetry:.3e}"),
    );
}

/// Frozen boundary golden (fs-lbm:d3q19-boundary-bits, registered in
/// golden-couplings.json). Frozen 2026-07-17 after all four policy quadrants
/// reproduced the value bit-identically on a committed tree: aarch64 M4 Pro
/// debug + release and x86-64 ts1 debug + release, at the d99f5a60/2919144
/// boundary (empty fs-lbm diff between them). Bump only with a semantic
/// change to the D3Q19 boundary bits, per docs/GOLDEN_POLICY.md.
const BOUNDARY_GOLDEN_HASH: u64 = 0x6683_0bf8_5cbb_3750;

/// Preimage for the boundary golden: link masks, densities, velocities, and
/// total mass across a voxelized-sphere lid cavity and a velocity/pressure
/// duct, hashed FNV-1a over exact bit patterns.
fn boundary_candidate_hash() -> u64 {
    let mut accumulator = 0xcbf2_9ce4_8422_2325u64;
    let mut feed = |bytes: &[u8]| {
        for byte in bytes {
            accumulator ^= u64::from(*byte);
            accumulator = accumulator.wrapping_mul(0x0000_0100_0000_01b3);
        }
    };
    let mut grid = BoundaryGrid3::new(
        8,
        8,
        8,
        0.75,
        [0.0; 3],
        BoundarySpec3::lid_cavity([0.025, 0.0, 0.0]),
    );
    let center = [3.5, 3.5, 3.5];
    grid.voxelize_sdf(|point| {
        let dx = point[0] - center[0];
        let dy = point[1] - center[1];
        let dz = point[2] - center[2];
        dx.mul_add(dx, dy.mul_add(dy, dz * dz)) - 0.75
    });
    grid.perturb(0x40_02_BC, 1e-4);
    grid.run(64);
    for (x, y, z) in [(0, 0, 0), (2, 5, 3), (6, 7, 4), (4, 1, 6)] {
        feed(&grid.link_mask(x, y, z).to_le_bytes());
        feed(&grid.density(x, y, z).to_bits().to_le_bytes());
        for component in grid.velocity(x, y, z) {
            feed(&component.to_bits().to_le_bytes());
        }
    }
    feed(&grid.total_mass().to_bits().to_le_bytes());

    let mut open = BoundaryGrid3::new(
        8,
        8,
        8,
        0.8,
        [0.0; 3],
        BoundarySpec3::velocity_pressure_duct_z([0.0, 0.0, 0.015], 1.0),
    );
    open.perturb(0x40_02_0B, 1e-4);
    open.run(32);
    for (x, y, z) in [(2, 3, 0), (5, 4, 7), (3, 3, 4)] {
        feed(&open.density(x, y, z).to_bits().to_le_bytes());
        for component in open.velocity(x, y, z) {
            feed(&component.to_bits().to_le_bytes());
        }
    }
    feed(&open.total_mass().to_bits().to_le_bytes());
    accumulator
}

#[test]
fn d3q19_boundary_candidate_hash_is_replay_stable() {
    let first = boundary_candidate_hash();
    let second = boundary_candidate_hash();
    println!("{{\"test\":\"lbm3bc-006-candidate\",\"hash\":\"{first:#018x}\"}}");
    assert_ne!(first, 0);
    assert_eq!(
        first, second,
        "D3Q19 boundary candidate hash is not replay-stable"
    );
    assert_eq!(
        first, BOUNDARY_GOLDEN_HASH,
        "d3q19 boundary bits changed: {first:#018x} vs {BOUNDARY_GOLDEN_HASH:#018x} — bump \
         only with a semantic boundary change and a same-commit registry re-pin \
         (docs/GOLDEN_POLICY.md)"
    );
}

/// Full no-slip leak gate from the bead: relative mass leak below 1e-11 over
/// 10,000 steps. Run explicitly in release.
#[test]
#[ignore = "release acceptance: 10k-step closed-domain leak (bead 40p2)"]
fn acceptance_no_slip_leak_10k() {
    let mut grid = BoundaryGrid3::new(
        8,
        8,
        16,
        0.8,
        [0.0, 0.0, 1e-6],
        BoundarySpec3::periodic_duct_z(),
    );
    grid.perturb(0x1E_A5, 1e-3);
    let initial = grid.total_mass();
    grid.run(10_000);
    let relative_leak = (grid.total_mass() - initial).abs() / initial;

    let mut voxel = BoundaryGrid3::new(8, 8, 8, 0.8, [0.0; 3], BoundarySpec3::periodic());
    voxel.voxelize_sdf(|point| {
        let dx = point[0] - 4.0;
        let dy = point[1] - 4.0;
        let dz = point[2] - 4.0;
        dx.mul_add(dx, dy.mul_add(dy, dz * dz)) - 2.0
    });
    voxel.perturb(0x1E_A6, 1e-3);
    let voxel_initial = voxel.total_mass();
    voxel.run(10_000);
    let voxel_leak = (voxel.total_mass() - voxel_initial).abs() / voxel_initial;
    verdict(
        "lbm3bc-acc-leak-10k",
        relative_leak < 1e-11 && voxel_leak < 1e-11,
        &format!("planar leak {relative_leak:.3e}; voxel leak {voxel_leak:.3e}"),
    );
}

/// Full pressure-boundary Poiseuille gate. The pressure gradient is derived
/// from `p = c_s^2 rho`; the mid-duct section must match the rectangular-duct
/// analytic series to 3% away from the inlet/outlet planes.
#[test]
#[ignore = "release acceptance: pressure-driven Poiseuille 3% gate (bead 40p2)"]
fn acceptance_pressure_poiseuille_within_three_percent() {
    // The frozen D3Q19 core uses the same 32x32 cross-section for its honest
    // full-rim 3% gate: the magic relaxation time cancels the leading
    // halfway-wall slip, while the remaining corner defect converges below
    // the bar only at the release resolution.
    let (nx, ny, nz) = (32, 32, 64);
    let (rho_in, rho_out) = (1.000_02, 0.999_98);
    let tau = 0.933_012_701_892_219_3;
    let mut grid = BoundaryGrid3::new(
        nx,
        ny,
        nz,
        tau,
        [0.0; 3],
        BoundarySpec3::pressure_duct_z(rho_in, rho_out),
    );
    grid.run(20_000);
    let section = grid.velocity_section_z(nz / 2, 2);
    let viscosity = grid.viscosity();
    let pressure_acceleration = CS2 * (rho_in - rho_out) / ((nz - 1) as f64);
    let mut max_relative = 0.0_f64;
    let mut max_interior = 0.0_f64;
    let mut max_rim = 0.0_f64;
    let mut argmax = (0usize, 0usize);
    for y in 0..ny {
        for x in 0..nx {
            let expected = duct_analytic(pressure_acceleration, viscosity, nx, ny, x, y);
            let relative = (section[y * nx + x] - expected).abs() / expected.abs();
            if x == 0 || y == 0 || x + 1 == nx || y + 1 == ny {
                max_rim = max_rim.max(relative);
            } else {
                max_interior = max_interior.max(relative);
            }
            if relative > max_relative {
                max_relative = relative;
                argmax = (x, y);
            }
        }
    }
    verdict(
        "lbm3bc-acc-pressure-poiseuille",
        max_relative < 0.03,
        &format!(
            "max rel {max_relative:.4} at {argmax:?}; interior {max_interior:.4}; rim {max_rim:.4}"
        ),
    );
}

#[test]
#[should_panic(expected = "periodic boundaries must be paired")]
fn periodic_faces_must_be_paired() {
    let wall = FaceBoundary3::stationary_wall();
    let invalid = BoundarySpec3::new([FaceBoundary3::Periodic, wall, wall, wall, wall, wall]);
    let _ = BoundaryGrid3::new(8, 8, 8, 0.8, [0.0; 3], invalid);
}

#[test]
#[should_panic(expected = "must be tangential")]
fn moving_wall_rejects_normal_velocity() {
    let wall = FaceBoundary3::stationary_wall();
    let invalid = BoundarySpec3::new([
        wall,
        wall,
        wall,
        FaceBoundary3::Wall {
            velocity: [0.0, 0.01, 0.0],
        },
        wall,
        wall,
    ]);
    let _ = BoundaryGrid3::new(8, 8, 8, 0.8, [0.0; 3], invalid);
}

#[test]
#[should_panic(expected = "currently require zero body force")]
fn open_boundaries_reject_body_force_until_guo_moments_are_reconstructed() {
    let _ = BoundaryGrid3::new(
        8,
        8,
        8,
        0.8,
        [0.0, 0.0, 1e-6],
        BoundarySpec3::pressure_duct_z(1.001, 0.999),
    );
}

#[test]
fn voxel_topology_commits_atomically_and_then_becomes_immutable() {
    use std::panic::{AssertUnwindSafe, catch_unwind};

    let mut grid = BoundaryGrid3::new(8, 8, 8, 0.8, [0.0; 3], BoundarySpec3::periodic());
    let masks_before = grid.tile_link_masks().to_vec();
    let failed = catch_unwind(AssertUnwindSafe(|| {
        grid.voxelize_sdf(|point| if point[2] > 4.0 { f64::NAN } else { 1.0 });
    }));
    assert!(failed.is_err());
    assert_eq!(grid.tile_link_masks(), masks_before);
    assert!((0..8).all(|z| (0..8).all(|y| (0..8).all(|x| !grid.is_solid(x, y, z)))));

    let all_solid = catch_unwind(AssertUnwindSafe(|| grid.voxelize_sdf(|_| -1.0)));
    assert!(all_solid.is_err());
    assert_eq!(grid.tile_link_masks(), masks_before);

    grid.voxelize_sdf(|point| {
        let dx = point[0] - 4.0;
        let dy = point[1] - 4.0;
        let dz = point[2] - 4.0;
        dx.mul_add(dx, dy.mul_add(dy, dz * dz)) - 1.0
    });
    let masks_after = grid.tile_link_masks().to_vec();
    let second = catch_unwind(AssertUnwindSafe(|| grid.voxelize_sdf(|_| 1.0)));
    assert!(second.is_err());
    assert_eq!(grid.tile_link_masks(), masks_after);
    verdict(
        "lbm3bc-007-topology-atomic",
        true,
        "failed samples left state unchanged; committed topology refused mutation",
    );
}

#[test]
fn obstructed_open_neighbor_rejects_topology_atomically() {
    use std::panic::{AssertUnwindSafe, catch_unwind};

    let mut grid = BoundaryGrid3::new(
        8,
        8,
        8,
        0.8,
        [0.0; 3],
        BoundarySpec3::pressure_duct_z(1.001, 0.999),
    );
    let masks_before = grid.tile_link_masks().to_vec();
    let failed = catch_unwind(AssertUnwindSafe(|| {
        grid.voxelize_sdf(|point| {
            if (point[2] - 1.5).abs() < 0.25 {
                -1.0
            } else {
                1.0
            }
        });
    }));
    assert!(failed.is_err());
    assert_eq!(grid.tile_link_masks(), masks_before);
    assert!((0..8).all(|z| (0..8).all(|y| (0..8).all(|x| !grid.is_solid(x, y, z)))));
    grid.voxelize_sdf(|_| 1.0);
}

#[test]
fn rejected_perturbation_leaves_state_and_topology_unlocked() {
    use std::panic::{AssertUnwindSafe, catch_unwind};

    let mut grid = BoundaryGrid3::new(8, 8, 8, 0.8, [0.0; 3], BoundarySpec3::periodic());
    let before = grid.populations(0, 0, 0);
    let failed = catch_unwind(AssertUnwindSafe(|| grid.perturb(1, 2.0)));
    assert!(failed.is_err());
    assert!(
        grid.populations(0, 0, 0)
            .into_iter()
            .zip(before)
            .all(|(after, before)| after.to_bits() == before.to_bits())
    );
    grid.voxelize_sdf(|_| 1.0);
}

#[test]
fn face_normals_and_tile_size_are_pinned() {
    assert_eq!(TILE, 4);
    assert_eq!(Face3::XMin.normal(), [-1, 0, 0]);
    assert_eq!(Face3::YMax.normal(), [0, 1, 0]);
    assert_eq!(Face3::ZMax.normal(), [0, 0, 1]);
}
