//! Battery for the D2Q9 lattice Boltzmann core (fs-lbm). Covers the
//! equilibrium moments, mass conservation, the load-bearing physical check
//! (steady Poiseuille channel flow matches the analytic parabola), and the
//! lattice-scaling assistant's stability bookkeeping.

use fs_lbm::core2::{
    CurvedMovingWallLink2, MovingWallMomentumExchange2, PartialSaturationCell2, VelocityPressureX2,
};
use fs_lbm::{
    CS2, Cell, Color, Grid, Lbm, MACH_LIMIT, Q, equilibrium, plan_scaling, poiseuille_analytic,
};

fn d2q9_nonequilibrium_stress(
    populations: &[f64; Q],
    rho: f64,
    velocity: [f64; 2],
) -> [[f64; 2]; 2] {
    const EX: [f64; Q] = [0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0];
    const EY: [f64; Q] = [0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0];
    let equilibrium = equilibrium(rho, velocity[0], velocity[1]);
    let mut stress = [[0.0; 2]; 2];
    for q in 0..Q {
        let nonequilibrium = populations[q] - equilibrium[q];
        let e = [EX[q], EY[q]];
        for row in 0..2 {
            for column in 0..2 {
                stress[row][column] += e[row] * e[column] * nonequilibrium;
            }
        }
    }
    stress
}

fn d2q9_active_raw_momentum(grid: &Grid) -> [f64; 2] {
    const EX: [f64; Q] = [0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0];
    const EY: [f64; Q] = [0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0];
    let mut momentum = [0.0; 2];
    for (populations, flag) in grid.f.iter().zip(&grid.flags) {
        if !matches!(*flag, Cell::Fluid | Cell::Interface) {
            continue;
        }
        for q in 0..Q {
            momentum[0] += EX[q] * populations[q];
            momentum[1] += EY[q] * populations[q];
        }
    }
    momentum
}

fn d2q9_curved_wall_links(
    grid: &Grid,
    fraction: f64,
    wall_velocities: &[[f64; 2]],
) -> Vec<[Option<CurvedMovingWallLink2>; Q]> {
    let mut links = vec![[None; Q]; grid.nx * grid.ny];
    for y in 0..grid.ny {
        for x in 0..grid.nx {
            let destination = grid.idx(x, y);
            if !matches!(grid.flags[destination], Cell::Fluid | Cell::Interface) {
                continue;
            }
            for q in 0..Q {
                let Some(source) = grid.source(x, y, q) else {
                    continue;
                };
                if grid.flags[source] == Cell::Wall {
                    links[destination][q] = Some(CurvedMovingWallLink2::new(
                        fraction,
                        wall_velocities[source],
                    ));
                }
            }
        }
    }
    links
}

const MOVING_COUETTE_NX: usize = 8;
const MOVING_COUETTE_NY: usize = 10;
const MOVING_COUETTE_SHEAR_SPEED: f64 = 0.04;

#[derive(Debug, Clone, PartialEq)]
struct MovingCouetteObservation {
    relative_profile: Vec<f64>,
    normalized_profile_error: f64,
    normalized_transverse_velocity: f64,
    top: MovingWallMomentumExchange2,
    bottom: MovingWallMomentumExchange2,
}

fn moving_couette_eta(y: usize) -> f64 {
    (y as f64 - 0.5) / (MOVING_COUETTE_NY - 2) as f64
}

#[allow(clippy::too_many_lines)]
fn observe_moving_couette(frame_speed: f64) -> MovingCouetteObservation {
    let mut grid = Grid::uniform(MOVING_COUETTE_NX, MOVING_COUETTE_NY, 0.8);
    grid.periodic_x = true;
    grid.periodic_y = false;

    let mut wall_velocities = vec![[0.0; 2]; MOVING_COUETTE_NX * MOVING_COUETTE_NY];
    let mut both_walls = vec![false; wall_velocities.len()];
    let mut top_wall = vec![false; wall_velocities.len()];
    let mut bottom_wall = vec![false; wall_velocities.len()];
    for x in 0..MOVING_COUETTE_NX {
        let bottom = grid.idx(x, 0);
        let top = grid.idx(x, MOVING_COUETTE_NY - 1);
        grid.flags[bottom] = Cell::Wall;
        grid.flags[top] = Cell::Wall;
        wall_velocities[bottom] = [frame_speed, 0.0];
        wall_velocities[top] = [frame_speed + MOVING_COUETTE_SHEAR_SPEED, 0.0];
        both_walls[bottom] = true;
        both_walls[top] = true;
        bottom_wall[bottom] = true;
        top_wall[top] = true;
    }

    // Start from the analytic linear velocity field. Its missing viscous
    // nonequilibrium relaxes during the warm-up; the observation is taken only
    // after more than six momentum-diffusion times across the fluid height.
    for y in 1..MOVING_COUETTE_NY - 1 {
        let velocity = frame_speed + MOVING_COUETTE_SHEAR_SPEED * moving_couette_eta(y);
        for x in 0..MOVING_COUETTE_NX {
            let cell = grid.idx(x, y);
            grid.f[cell] = equilibrium(1.0, velocity, 0.0);
        }
    }

    let mut scratch = Vec::new();
    let origin = [
        0.5 * (MOVING_COUETTE_NX - 1) as f64,
        0.5 * (MOVING_COUETTE_NY - 1) as f64,
    ];
    for _ in 0..4_000 {
        let _receipt = grid.step_with_moving_wall_momentum(
            &mut scratch,
            &both_walls,
            &wall_velocities,
            origin,
        );
    }

    let mut relative_profile = Vec::with_capacity(MOVING_COUETTE_NY - 2);
    let mut normalized_profile_error = 0.0_f64;
    let mut normalized_transverse_velocity = 0.0_f64;
    for y in 1..MOVING_COUETTE_NY - 1 {
        let mut mean_velocity = [0.0; 2];
        for x in 0..MOVING_COUETTE_NX {
            let velocity = grid.moments(grid.idx(x, y)).u;
            mean_velocity[0] += velocity[0];
            mean_velocity[1] += velocity[1];
        }
        mean_velocity[0] /= MOVING_COUETTE_NX as f64;
        mean_velocity[1] /= MOVING_COUETTE_NX as f64;
        let relative_velocity = mean_velocity[0] - frame_speed;
        let analytic = MOVING_COUETTE_SHEAR_SPEED * moving_couette_eta(y);
        normalized_profile_error = normalized_profile_error
            .max((relative_velocity - analytic).abs() / MOVING_COUETTE_SHEAR_SPEED);
        normalized_transverse_velocity =
            normalized_transverse_velocity.max(mean_velocity[1].abs() / MOVING_COUETTE_SHEAR_SPEED);
        relative_profile.push(relative_velocity);
    }

    // Probe both walls from one identical post-collision state so differences
    // are selection effects only, not consecutive time-step drift.
    let mut post = Vec::new();
    grid.collide_into(&mut post);
    let mut top_probe = grid.clone();
    let top =
        top_probe.stream_from_with_moving_wall_momentum(&post, &top_wall, &wall_velocities, origin);
    let mut bottom_probe = grid;
    let bottom = bottom_probe.stream_from_with_moving_wall_momentum(
        &post,
        &bottom_wall,
        &wall_velocities,
        origin,
    );

    MovingCouetteObservation {
        relative_profile,
        normalized_profile_error,
        normalized_transverse_velocity,
        top,
        bottom,
    }
}

fn moving_wall_receipt_balance_error(receipt: &MovingWallMomentumExchange2) -> f64 {
    (0..2).fold(0.0_f64, |error, axis| {
        error.max(
            (receipt.wall_impulse[axis] + receipt.fluid_population_impulse[axis]
                - receipt.wall_velocity_mass_impulse[axis])
                .abs(),
        )
    })
}

fn normalized_pair_residual(left: f64, right: f64) -> f64 {
    (left - right).abs() / left.abs().max(right.abs())
}

#[test]
fn the_equilibrium_recovers_its_moments() {
    let (rho, ux, uy) = (1.0, 0.05, -0.02);
    let f = equilibrium(rho, ux, uy);
    let sum: f64 = f.iter().sum();
    assert!((sum - rho).abs() < 1e-12); // density
    // momentum: Σ eₓ fᵢ = ρ uₓ (D2Q9 velocities).
    let ex = [0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0];
    let ey = [0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0];
    let mx: f64 = f.iter().zip(ex).map(|(fi, e)| fi * e).sum();
    let my: f64 = f.iter().zip(ey).map(|(fi, e)| fi * e).sum();
    assert!((mx - rho * ux).abs() < 1e-12 && (my - rho * uy).abs() < 1e-12);
}

#[test]
fn mass_is_conserved() {
    let mut lbm = Lbm::channel(6, 12, 0.8, 1e-4);
    let m0 = lbm.total_mass();
    lbm.run(200);
    // Mass is conserved BY CONSTRUCTION (collision, forcing, streaming, and
    // bounce-back all preserve Σf), so the only drift is summation roundoff:
    // measured 9.38e-13 over 200 steps on mass 72, BIT-IDENTICAL on aarch64
    // (M4 Pro) and x86-64 (Threadripper 5975WX). Gate at 1e-11 (~10x that
    // roundoff floor) so the CONTRACT's "mass is conserved" claim is verified
    // to roundoff and a future systematic per-step leak is actually caught —
    // the old 1e-9 bound was ~1000x loose and would have passed a real
    // ~5e-12/step leak.
    assert!((lbm.total_mass() - m0).abs() < 1e-11, "mass drifted");
    assert!((m0 - f64::from(6 * 12)).abs() < 1e-9); // unit density
}

#[test]
fn poiseuille_flow_matches_the_analytic_parabola() {
    let (nx, ny, tau, gx) = (4, 25, 0.8, 1e-5);
    let mut lbm = Lbm::channel(nx, ny, tau, gx);
    lbm.run(20_000); // reach steady state
    let profile = lbm.x_velocity_profile();
    let nu = lbm.viscosity();
    // the profile matches the analytic parabola at every row (halfway
    // bounce-back resolves the quadratic exactly).
    let mut max_rel = 0.0_f64;
    for (y, &u) in profile.iter().enumerate() {
        let a = poiseuille_analytic(gx, nu, ny, y);
        max_rel = max_rel.max((u - a).abs() / a.abs());
    }
    assert!(max_rel < 0.03, "profile off by {max_rel:.4}");
    // and it is a parabola: symmetric with its peak at the centre.
    let mid = ny / 2;
    assert!(profile[mid] > profile[0] && profile[mid] > profile[ny - 1]);
    assert!((profile[1] - profile[ny - 2]).abs() / profile[mid] < 0.02);
}

#[test]
fn the_scaling_assistant_derives_tau_and_flags_stability() {
    // Re 100, L 40 lu, u 0.05 -> nu 0.02, tau 0.56, low Mach -> stable.
    let plan = plan_scaling(100.0, 40.0, 0.05);
    assert!((plan.viscosity - 0.02).abs() < 1e-12);
    assert!((plan.tau - 0.56).abs() < 1e-12);
    assert!(plan.stable && plan.tau_margin > 0.0);
    assert!(plan.mach < MACH_LIMIT);
    // a comfortably-stable plan is verified-color.
    assert!(matches!(plan.color(), Color::Verified { .. }));
}

#[test]
fn the_scaling_assistant_rejects_a_high_mach_plan() {
    // too large a lattice velocity breaks the low-Mach (incompressible) regime.
    let plan = plan_scaling(100.0, 20.0, 0.25);
    assert!(plan.mach > MACH_LIMIT);
    assert!(!plan.stable);
    // an unstable plan is not verified-color.
    assert!(matches!(plan.color(), Color::Estimated { .. }));
}

#[test]
#[should_panic(expected = "must be positive")]
fn the_scaling_assistant_rejects_nonsense() {
    let _ = plan_scaling(-1.0, 20.0, 0.1);
}

#[test]
fn the_solver_is_deterministic() {
    let mut a = Lbm::channel(4, 10, 0.7, 1e-4);
    let mut b = Lbm::channel(4, 10, 0.7, 1e-4);
    a.run(100);
    b.run(100);
    assert_eq!(a.x_velocity_profile(), b.x_velocity_profile());
}

#[test]
fn d2q9_wall_momentum_has_exact_link_sign_and_obstacle_selection() {
    let mut grid = Grid::uniform(3, 3, 0.8);
    grid.periodic_x = false;
    grid.periodic_y = false;
    grid.flags.fill(Cell::Gas);
    let fluid = grid.idx(1, 1);
    let left_wall = grid.idx(0, 1);
    let right_wall = grid.idx(2, 1);
    let top_wall = grid.idx(1, 2);
    grid.flags[fluid] = Cell::Fluid;
    grid.flags[left_wall] = Cell::Wall;
    grid.flags[right_wall] = Cell::Wall;
    grid.flags[top_wall] = Cell::Wall;

    let mut post = vec![[0.0; Q]; grid.nx * grid.ny];
    // D2Q9 directions 1 and 3 point east and west. An east-going population
    // of 0.25 transfers +2f = +0.5 x-momentum to the right wall.
    post[fluid][1] = 0.25;
    // The west-going population would transfer -0.75 to the left wall.
    post[fluid][3] = 0.375;
    // A north-going population of 0.125 transfers +0.25 y-momentum to the
    // top wall, independently pinning the lift-axis sign.
    post[fluid][2] = 0.125;

    let mut right_only = vec![false; grid.nx * grid.ny];
    right_only[right_wall] = true;
    let right_receipt = grid.stream_from_with_wall_momentum(&post, &right_only);
    assert_eq!(right_receipt.measured_links, 1);
    assert_eq!(right_receipt.wall_impulse[0].to_bits(), 0.5f64.to_bits());
    assert_eq!(right_receipt.wall_impulse[1].to_bits(), 0.0f64.to_bits());
    assert_eq!(grid.f[fluid][3].to_bits(), 0.25f64.to_bits());

    let mut both = right_only;
    both[left_wall] = true;
    let both_receipt = grid.stream_from_with_wall_momentum(&post, &both);
    assert_eq!(both_receipt.measured_links, 2);
    assert_eq!(both_receipt.wall_impulse[0].to_bits(), (-0.25f64).to_bits());
    assert_eq!(both_receipt.wall_impulse[1].to_bits(), 0.0f64.to_bits());
    assert_eq!(grid.f[fluid][1].to_bits(), 0.375f64.to_bits());

    let mut top_only = vec![false; grid.nx * grid.ny];
    top_only[top_wall] = true;
    let top_receipt = grid.stream_from_with_wall_momentum(&post, &top_only);
    assert_eq!(top_receipt.measured_links, 1);
    assert_eq!(top_receipt.wall_impulse[0].to_bits(), 0.0f64.to_bits());
    assert_eq!(top_receipt.wall_impulse[1].to_bits(), 0.25f64.to_bits());
    assert_eq!(grid.f[fluid][4].to_bits(), 0.125f64.to_bits());
}

#[test]
fn d2q9_wall_momentum_is_zero_at_rest_and_replays_bitwise() {
    let mut resting = Grid::uniform(5, 5, 0.8);
    let resting_wall = resting.idx(2, 2);
    resting.flags[resting_wall] = Cell::Wall;
    let mut resting_mask = vec![false; resting.nx * resting.ny];
    resting_mask[resting_wall] = true;
    let resting_post = resting.f.clone();
    let resting_receipt = resting.stream_from_with_wall_momentum(&resting_post, &resting_mask);
    assert_eq!(resting_receipt.measured_links, 8);
    assert!(resting_receipt.wall_impulse[0].abs() < 1e-15);
    assert!(resting_receipt.wall_impulse[1].abs() < 1e-15);

    let mut first = Grid::uniform(7, 7, 0.8);
    let wall = first.idx(3, 3);
    let upstream = first.idx(2, 3);
    first.flags[wall] = Cell::Wall;
    first.f[upstream] = equilibrium(1.0, 0.04, 0.01);
    let mut second = first.clone();
    let mut legacy = first.clone();
    let mut mask = vec![false; first.nx * first.ny];
    mask[wall] = true;
    let (mut first_scratch, mut second_scratch) = (Vec::new(), Vec::new());
    let mut legacy_scratch = Vec::new();
    let mut observed_nonzero_impulse = false;

    for _ in 0..12 {
        let first_receipt = first.step_with_wall_momentum(&mut first_scratch, &mask);
        let second_receipt = second.step_with_wall_momentum(&mut second_scratch, &mask);
        legacy.step(&mut legacy_scratch);
        observed_nonzero_impulse |= first_receipt.wall_impulse != [0.0, 0.0];
        assert_eq!(first_receipt.measured_links, second_receipt.measured_links);
        assert_eq!(
            first_receipt.wall_impulse.map(f64::to_bits),
            second_receipt.wall_impulse.map(f64::to_bits)
        );
    }
    assert!(observed_nonzero_impulse);
    for (first_cell, (second_cell, legacy_cell)) in
        first.f.iter().zip(second.f.iter().zip(&legacy.f))
    {
        assert_eq!(first_cell.map(f64::to_bits), second_cell.map(f64::to_bits));
        assert_eq!(first_cell.map(f64::to_bits), legacy_cell.map(f64::to_bits));
    }
}

#[test]
fn d2q9_moving_wall_exchange_pins_relative_force_torque_work_and_balance() {
    // G0/G3: one independently enumerable link pins direction conventions,
    // the moving bounce correction, relative-velocity force, and the complete
    // moving-mass momentum balance.
    let mut grid = Grid::uniform(3, 3, 0.8);
    grid.periodic_x = false;
    grid.periodic_y = false;
    grid.flags.fill(Cell::Gas);
    let fluid = grid.idx(1, 1);
    let wall = grid.idx(2, 1);
    grid.flags[fluid] = Cell::Fluid;
    grid.flags[wall] = Cell::Wall;

    let mut post = vec![[0.0; Q]; grid.nx * grid.ny];
    let outgoing = 0.25;
    post[fluid][1] = outgoing;
    let mut measured = vec![false; grid.nx * grid.ny];
    measured[wall] = true;
    let mut wall_velocities = vec![[0.0; 2]; grid.nx * grid.ny];
    let wall_velocity = [0.03, 0.02];
    wall_velocities[wall] = wall_velocity;
    let origin = [0.5, 0.25];

    let receipt =
        grid.stream_from_with_moving_wall_momentum(&post, &measured, &wall_velocities, origin);

    // Pull direction q=3 points west from this right-hand wall into fluid.
    // Its D2Q9 weight is 1/9 and the post-collision density is 0.25.
    let incoming = outgoing + 2.0 * (1.0 / 9.0) * outgoing * -wall_velocity[0] / CS2;
    let expected_wall = [
        (1.0 - wall_velocity[0]) * outgoing - (-1.0 - wall_velocity[0]) * incoming,
        (0.0 - wall_velocity[1]) * outgoing - (0.0 - wall_velocity[1]) * incoming,
    ];
    let expected_fluid = [-(incoming + outgoing), 0.0];
    let expected_mass_change = incoming - outgoing;
    let expected_mass_impulse = [
        wall_velocity[0] * expected_mass_change,
        wall_velocity[1] * expected_mass_change,
    ];
    let link_offset = [1.5 - origin[0], 1.0 - origin[1]];
    let expected_wall_torque =
        link_offset[0] * expected_wall[1] - link_offset[1] * expected_wall[0];
    let expected_work = expected_wall[0] * wall_velocity[0] + expected_wall[1] * wall_velocity[1];

    assert_eq!(receipt.measured_links, 1);
    assert!((grid.f[fluid][3] - incoming).abs() < 1e-15);
    for axis in 0..2 {
        assert!((receipt.wall_impulse[axis] - expected_wall[axis]).abs() < 1e-15);
        assert!((receipt.fluid_population_impulse[axis] - expected_fluid[axis]).abs() < 1e-15);
        assert!(
            (receipt.wall_velocity_mass_impulse[axis] - expected_mass_impulse[axis]).abs() < 1e-15
        );
        assert!(
            (receipt.wall_impulse[axis] + receipt.fluid_population_impulse[axis]
                - receipt.wall_velocity_mass_impulse[axis])
                .abs()
                < 1e-15
        );
    }
    assert!((receipt.fluid_mass_change - expected_mass_change).abs() < 1e-15);
    assert!((receipt.wall_angular_impulse - expected_wall_torque).abs() < 1e-15);
    assert!(
        (receipt.wall_angular_impulse + receipt.fluid_population_angular_impulse
            - receipt.wall_velocity_mass_angular_impulse)
            .abs()
            < 1e-15
    );
    assert!((receipt.wall_work - expected_work).abs() < 1e-15);

    let shifted_origin = [1.0, -0.5];
    let shifted = grid.stream_from_with_moving_wall_momentum(
        &post,
        &measured,
        &wall_velocities,
        shifted_origin,
    );
    let origin_shift = [shifted_origin[0] - origin[0], shifted_origin[1] - origin[1]];
    let translated_torque = receipt.wall_angular_impulse
        - origin_shift[0] * receipt.wall_impulse[1]
        + origin_shift[1] * receipt.wall_impulse[0];
    assert!((shifted.wall_angular_impulse - translated_torque).abs() < 1e-15);
}

#[test]
#[allow(clippy::too_many_lines)]
fn d2q9_moving_couette_measures_common_boost_residual() {
    // G1/G3/G5: this is one deliberately narrow, low-Mach envelope. It checks
    // a manufactured linear Couette profile, compares two inertial frames,
    // and retains separate top/bottom exchange receipts. It is not a general
    // Galilean-invariance claim for arbitrary geometry or transient motion.
    let reference = observe_moving_couette(-0.01);
    let boosted = observe_moving_couette(0.01);
    let replay = observe_moving_couette(-0.01);
    assert_eq!(reference, replay, "fixed-frame deck must replay exactly");

    for observation in [&reference, &boosted] {
        assert!(
            observation.normalized_profile_error < 0.03,
            "linear Couette profile residual exceeded 3% of wall-speed difference: {}",
            observation.normalized_profile_error
        );
        assert!(
            observation.normalized_transverse_velocity < 1.0e-10,
            "transverse velocity exceeded the symmetric-channel roundoff band: {}",
            observation.normalized_transverse_velocity
        );
        assert_eq!(
            observation.top.measured_links,
            3 * MOVING_COUETTE_NX,
            "top wall must contribute one normal and two diagonal links per column"
        );
        assert_eq!(
            observation.bottom.measured_links,
            3 * MOVING_COUETTE_NX,
            "bottom wall must contribute one normal and two diagonal links per column"
        );
        assert!(
            moving_wall_receipt_balance_error(&observation.top) < 1.0e-12,
            "top moving-wall receipt violated its discrete momentum balance"
        );
        assert!(
            moving_wall_receipt_balance_error(&observation.bottom) < 1.0e-12,
            "bottom moving-wall receipt violated its discrete momentum balance"
        );

        let shear_force_scale = observation.top.wall_impulse[0]
            .abs()
            .max(observation.bottom.wall_impulse[0].abs());
        assert!(
            shear_force_scale > 1.0e-8,
            "Couette deck must resolve nonzero wall shear"
        );
        let shear_balance = (observation.top.wall_impulse[0] + observation.bottom.wall_impulse[0])
            .abs()
            / shear_force_scale;
        assert!(
            shear_balance < 0.03,
            "top and bottom shear impulses failed to balance within 3%: {shear_balance}"
        );

        let normal_force_scale = observation.top.wall_impulse[1]
            .abs()
            .max(observation.bottom.wall_impulse[1].abs());
        assert!(
            normal_force_scale > 0.0,
            "wall pressure impulse must be nonzero"
        );
        let normal_balance = (observation.top.wall_impulse[1] + observation.bottom.wall_impulse[1])
            .abs()
            / normal_force_scale;
        assert!(
            normal_balance < 1.0e-12,
            "symmetric wall-pressure impulses failed to cancel: {normal_balance}"
        );
    }

    let profile_frame_residual = reference
        .relative_profile
        .iter()
        .zip(&boosted.relative_profile)
        .fold(0.0_f64, |residual, (&reference_u, &boosted_u)| {
            residual.max((reference_u - boosted_u).abs() / MOVING_COUETTE_SHEAR_SPEED)
        });
    let top_impulse_frame_residual =
        normalized_pair_residual(reference.top.wall_impulse[0], boosted.top.wall_impulse[0]);
    let bottom_impulse_frame_residual = normalized_pair_residual(
        reference.bottom.wall_impulse[0],
        boosted.bottom.wall_impulse[0],
    );
    let impulse_frame_residual = top_impulse_frame_residual.max(bottom_impulse_frame_residual);

    assert!(
        profile_frame_residual < 0.03,
        "common boost changed the relative Couette profile by at least 3%: {profile_frame_residual}"
    );
    assert!(
        impulse_frame_residual < 0.03,
        "common boost changed the boundary-relative shear impulse by at least 3%: \
         {impulse_frame_residual}"
    );

    println!(
        concat!(
            "{{\"deck\":\"d2q9-moving-couette-common-boost\",",
            "\"common_boost\":0.02,\"profile_residual\":{:.17e},",
            "\"wall_impulse_residual\":{:.17e},",
            "\"reference_profile_error\":{:.17e},\"boosted_profile_error\":{:.17e},",
            "\"links_per_wall\":{}}}"
        ),
        profile_frame_residual,
        impulse_frame_residual,
        reference.normalized_profile_error,
        boosted.normalized_profile_error,
        reference.top.measured_links
    );
}

#[test]
#[allow(clippy::too_many_lines)]
fn d2q9_partial_saturation_pins_noble_torczynski_collision_and_exchange() {
    const OPPOSITE: [usize; Q] = [0, 3, 4, 1, 2, 7, 8, 5, 6];
    const EX: [f64; Q] = [0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0];
    const EY: [f64; Q] = [0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0];

    // G0/G3: independently enumerate every population of one nonsymmetric
    // cell so neither the fluid nor solid collision can pass by cancellation.
    let mut grid = Grid::uniform(3, 3, 0.8);
    let coupled = grid.idx(1, 1);
    let mut populations = equilibrium(1.1, 0.04, -0.02);
    populations[1] += 0.006;
    populations[3] -= 0.002;
    populations[5] += 0.003;
    populations[7] -= 0.001;
    grid.f[coupled] = populations;

    let solid_fraction = 0.4;
    let solid_velocity = [-0.01, 0.015];
    let descriptor = PartialSaturationCell2::new(solid_fraction, solid_velocity);
    assert_eq!(
        descriptor.solid_fraction().to_bits(),
        solid_fraction.to_bits()
    );
    assert_eq!(
        descriptor.solid_velocity().map(f64::to_bits),
        solid_velocity.map(f64::to_bits)
    );
    let mut cells = vec![PartialSaturationCell2::default(); grid.nx * grid.ny];
    cells[coupled] = descriptor;

    let moments = grid.moments(coupled);
    let fluid_equilibrium = equilibrium(moments.rho, moments.u[0], moments.u[1]);
    let solid_equilibrium = equilibrium(moments.rho, solid_velocity[0], solid_velocity[1]);
    let relaxation_margin = grid.tau[coupled] - 0.5;
    let coupling_weight =
        solid_fraction * relaxation_margin / ((1.0 - solid_fraction) + relaxation_margin);
    let origin = [0.25, -0.5];
    let mut post = vec![[123.0; Q]; 2];
    let receipt = grid.collide_into_with_partial_saturation(&mut post, &cells, origin);

    let mut ordinary_post = Vec::new();
    grid.collide_into(&mut ordinary_post);
    let mut observed_mass_change = 0.0;
    let mut observed_fluid_impulse = [0.0; 2];
    for q in 0..Q {
        let fluid_collision = (fluid_equilibrium[q] - populations[q]) / grid.tau[coupled];
        let solid_collision = populations[OPPOSITE[q]] - populations[q] + solid_equilibrium[q]
            - fluid_equilibrium[OPPOSITE[q]];
        let independently_blended = populations[q]
            + (1.0 - coupling_weight) * fluid_collision
            + coupling_weight * solid_collision;
        assert!(
            (post[coupled][q] - independently_blended).abs() < 5.0e-15,
            "population {q} did not match the independently blended PSC formula"
        );

        let observed_correction = post[coupled][q] - ordinary_post[coupled][q];
        observed_mass_change += observed_correction;
        observed_fluid_impulse[0] += EX[q] * observed_correction;
        observed_fluid_impulse[1] += EY[q] * observed_correction;
    }

    let analytic_fluid_impulse = [
        coupling_weight * moments.rho * (solid_velocity[0] - moments.u[0]),
        coupling_weight * moments.rho * (solid_velocity[1] - moments.u[1]),
    ];
    assert!(receipt.fluid_population_mass_change.abs() < 2.0e-14);
    assert!((observed_mass_change - receipt.fluid_population_mass_change).abs() < 2.0e-15);
    for axis in 0..2 {
        assert!(
            (receipt.fluid_population_impulse[axis] - observed_fluid_impulse[axis]).abs() < 2.0e-15
        );
        assert!(
            (receipt.fluid_population_impulse[axis] - analytic_fluid_impulse[axis]).abs() < 2.0e-14
        );
        assert!(
            (receipt.solid_impulse[axis] + receipt.fluid_population_impulse[axis]).abs() < 2.0e-15
        );
    }
    assert_eq!(receipt.coupled_cells, 1);
    assert!((receipt.coupling_weight_sum - coupling_weight).abs() < 1.0e-15);

    let offset = [1.0 - origin[0], 1.0 - origin[1]];
    let expected_solid_torque =
        offset[0] * receipt.solid_impulse[1] - offset[1] * receipt.solid_impulse[0];
    let expected_work =
        receipt.solid_impulse[0] * solid_velocity[0] + receipt.solid_impulse[1] * solid_velocity[1];
    assert!((receipt.solid_angular_impulse - expected_solid_torque).abs() < 2.0e-15);
    assert!(
        (receipt.solid_angular_impulse + receipt.fluid_population_angular_impulse).abs() < 2.0e-15
    );
    assert!((receipt.solid_work - expected_work).abs() < 2.0e-15);

    let shifted_origin = [1.25, 0.0];
    let mut shifted_post = Vec::new();
    let shifted =
        grid.collide_into_with_partial_saturation(&mut shifted_post, &cells, shifted_origin);
    for (unshifted, shifted) in post.iter().zip(&shifted_post) {
        assert_eq!(unshifted.map(f64::to_bits), shifted.map(f64::to_bits));
    }
    let origin_shift = [shifted_origin[0] - origin[0], shifted_origin[1] - origin[1]];
    let translated_torque = receipt.solid_angular_impulse
        - origin_shift[0] * receipt.solid_impulse[1]
        + origin_shift[1] * receipt.solid_impulse[0];
    assert!((shifted.solid_angular_impulse - translated_torque).abs() < 2.0e-15);
}

#[test]
fn d2q9_partial_saturation_preserves_endpoints_and_replays() {
    // The zero-fraction seam is operation-tree compatible with ordinary BGK.
    let mut base = Grid::uniform(4, 4, 0.8);
    let perturbed = base.idx(2, 1);
    base.f[perturbed] = equilibrium(1.05, 0.03, -0.01);
    base.f[perturbed][1] += 0.002;
    base.f[perturbed][3] -= 0.001;
    let zero_cells = vec![PartialSaturationCell2::default(); base.nx * base.ny];
    let mut ordinary_post = Vec::new();
    base.collide_into(&mut ordinary_post);
    let mut partial_post = vec![[77.0; Q]; 1];
    let zero_receipt =
        base.collide_into_with_partial_saturation(&mut partial_post, &zero_cells, [0.0; 2]);
    assert_eq!(zero_receipt.coupled_cells, 0);
    assert_eq!(zero_receipt.coupling_weight_sum.to_bits(), 0.0f64.to_bits());
    for (ordinary, partial) in ordinary_post.iter().zip(&partial_post) {
        assert_eq!(ordinary.map(f64::to_bits), partial.map(f64::to_bits));
    }

    let mut ordinary_step = base.clone();
    let mut partial_step = base.clone();
    let (mut ordinary_scratch, mut partial_scratch) = (Vec::new(), vec![[91.0; Q]; 2]);
    ordinary_step.step(&mut ordinary_scratch);
    let step_receipt =
        partial_step.step_with_partial_saturation(&mut partial_scratch, &zero_cells, [0.0; 2]);
    assert_eq!(step_receipt.coupled_cells, 0);
    for (ordinary, partial) in ordinary_step.f.iter().zip(&partial_step.f) {
        assert_eq!(ordinary.map(f64::to_bits), partial.map(f64::to_bits));
    }
    for (ordinary, partial) in ordinary_scratch.iter().zip(&partial_scratch) {
        assert_eq!(ordinary.map(f64::to_bits), partial.map(f64::to_bits));
    }

    // At full saturation, a stationary equilibrium is a fixed point and the
    // nonlinear Noble-Torczynski weight reaches exactly one.
    let full_grid = Grid::uniform(3, 2, 1.0);
    let full_cells = vec![PartialSaturationCell2::new(1.0, [0.0; 2]); 6];
    let mut full_post = Vec::new();
    let full_receipt =
        full_grid.collide_into_with_partial_saturation(&mut full_post, &full_cells, [1.0, 0.5]);
    assert_eq!(full_receipt.coupled_cells, 6);
    assert_eq!(full_receipt.coupling_weight_sum.to_bits(), 6.0f64.to_bits());
    assert!(full_receipt.fluid_population_mass_change.abs() < 1.0e-15);
    assert!(full_receipt.fluid_population_impulse == [0.0; 2]);
    assert!(full_receipt.solid_impulse == [0.0; 2]);
    for (before, after) in full_grid.f.iter().zip(&full_post) {
        assert_eq!(before.map(f64::to_bits), after.map(f64::to_bits));
    }

    // G5: nonzero saturation and moving-solid forcing replay bit-for-bit in
    // fixed row-major cell/direction order.
    let mut first = Grid::uniform(5, 4, 0.8);
    let first_cell = first.idx(1, 1);
    let second_cell = first.idx(3, 2);
    first.f[first_cell] = equilibrium(1.0, 0.025, -0.01);
    first.f[second_cell] = equilibrium(0.98, -0.015, 0.005);
    let mut replay_cells = vec![PartialSaturationCell2::default(); first.nx * first.ny];
    replay_cells[first_cell] = PartialSaturationCell2::new(0.25, [-0.01, 0.02]);
    replay_cells[second_cell] = PartialSaturationCell2::new(0.55, [0.015, -0.01]);
    let mut second = first.clone();
    let (mut first_scratch, mut second_scratch) = (Vec::new(), Vec::new());
    let mut observed_nonzero_impulse = false;
    for _ in 0..12 {
        let first_receipt =
            first.step_with_partial_saturation(&mut first_scratch, &replay_cells, [2.0, 1.5]);
        let second_receipt =
            second.step_with_partial_saturation(&mut second_scratch, &replay_cells, [2.0, 1.5]);
        observed_nonzero_impulse |= first_receipt.solid_impulse != [0.0; 2];
        assert_eq!(first_receipt, second_receipt);
        for (first_population, second_population) in first.f.iter().zip(&second.f) {
            assert_eq!(
                first_population.map(f64::to_bits),
                second_population.map(f64::to_bits)
            );
        }
    }
    assert!(observed_nonzero_impulse);
}

#[test]
fn d2q9_partial_saturation_refuses_bad_requests_atomically() {
    let invalid_descriptors = [
        (f64::NAN, [0.0, 0.0]),
        (-0.1, [0.0, 0.0]),
        (1.1, [0.0, 0.0]),
        (0.0, [0.01, 0.0]),
        (0.5, [f64::NAN, 0.0]),
        (0.5, [0.2, 0.0]),
    ];
    for (solid_fraction, solid_velocity) in invalid_descriptors {
        let refusal = std::panic::catch_unwind(|| {
            let _ = PartialSaturationCell2::new(solid_fraction, solid_velocity);
        });
        assert!(refusal.is_err());
    }

    let mut grid = Grid::uniform(3, 3, 0.8);
    let original_populations = grid.f.clone();
    let sentinel = vec![[71.0; Q]; 2];
    let mut scratch = sentinel.clone();
    let short = vec![PartialSaturationCell2::default(); grid.nx * grid.ny - 1];
    let short_refusal = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        let _ = grid.step_with_partial_saturation(&mut scratch, &short, [0.0; 2]);
    }));
    assert!(short_refusal.is_err());
    assert_eq!(scratch, sentinel);
    assert_eq!(grid.f, original_populations);

    let wall = grid.idx(1, 1);
    grid.flags[wall] = Cell::Wall;
    let mut wall_cells = vec![PartialSaturationCell2::default(); grid.nx * grid.ny];
    wall_cells[wall] = PartialSaturationCell2::new(0.25, [0.0; 2]);
    let mut wall_post = sentinel.clone();
    let wall_refusal = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        let _ = grid.collide_into_with_partial_saturation(&mut wall_post, &wall_cells, [0.0; 2]);
    }));
    assert!(wall_refusal.is_err());
    assert_eq!(wall_post, sentinel);

    grid.flags[wall] = Cell::Fluid;
    grid.g = [1.0e-5, 0.0];
    let zero_cells = vec![PartialSaturationCell2::default(); grid.nx * grid.ny];
    let mut forced_post = sentinel.clone();
    let force_refusal = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        let _ = grid.collide_into_with_partial_saturation(&mut forced_post, &zero_cells, [0.0; 2]);
    }));
    assert!(force_refusal.is_err());
    assert_eq!(forced_post, sentinel);

    grid.g = [0.0; 2];
    let mut bad_origin_post = sentinel.clone();
    let origin_refusal = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        let _ = grid.collide_into_with_partial_saturation(
            &mut bad_origin_post,
            &zero_cells,
            [f64::NAN, 0.0],
        );
    }));
    assert!(origin_refusal.is_err());
    assert_eq!(bad_origin_post, sentinel);
    assert_eq!(grid.f, original_populations);

    // Zero saturation still traverses ordinary BGK collision. A finite state
    // whose finite momentum produces an overflowing equilibrium must refuse
    // before publishing the otherwise-unrelated caller scratch buffer.
    let mut pathological = Grid::uniform(2, 2, 0.8);
    let pathological_cell = pathological.idx(0, 0);
    let mut extreme = [0.0; Q];
    extreme[1] = 1.0e200;
    extreme[3] = -1.0e200;
    extreme[8] = 1.0;
    pathological.f[pathological_cell] = extreme;
    let pathological_before = pathological.f.clone();
    let pathological_cells =
        vec![PartialSaturationCell2::default(); pathological.nx * pathological.ny];
    let mut pathological_post = sentinel.clone();
    let proposal_refusal = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        let _ = pathological.collide_into_with_partial_saturation(
            &mut pathological_post,
            &pathological_cells,
            [0.0; 2],
        );
    }));
    assert!(proposal_refusal.is_err());
    assert_eq!(pathological_post, sentinel);
    assert_eq!(pathological.f, pathological_before);
}

#[test]
fn d2q9_zero_velocity_moving_api_matches_stationary_api_and_replays() {
    // G5: zero velocity is an exact compatibility boundary, and a nonzero
    // moving-wall step replays bit-for-bit in fixed cell/link order.
    let mut legacy = Grid::uniform(5, 5, 0.8);
    let wall = legacy.idx(2, 2);
    let upstream = legacy.idx(1, 2);
    legacy.flags[wall] = Cell::Wall;
    legacy.f[upstream] = equilibrium(1.0, 0.03, -0.01);
    let post = legacy.f.clone();
    let mut moving = legacy.clone();
    let mut measured = vec![false; legacy.nx * legacy.ny];
    measured[wall] = true;
    let stationary_velocities = vec![[0.0; 2]; legacy.nx * legacy.ny];

    let stationary = legacy.stream_from_with_wall_momentum(&post, &measured);
    let through_moving = moving.stream_from_with_moving_wall_momentum(
        &post,
        &measured,
        &stationary_velocities,
        [0.0; 2],
    );
    assert_eq!(
        stationary.wall_impulse.map(f64::to_bits),
        through_moving.wall_impulse.map(f64::to_bits)
    );
    for (legacy_cell, moving_cell) in legacy.f.iter().zip(&moving.f) {
        assert_eq!(legacy_cell.map(f64::to_bits), moving_cell.map(f64::to_bits));
    }
    assert_eq!(through_moving.fluid_mass_change.to_bits(), 0.0f64.to_bits());
    assert_eq!(
        through_moving.wall_velocity_mass_impulse.map(f64::to_bits),
        [0, 0]
    );

    let mut first = Grid::uniform(7, 7, 0.8);
    let moving_wall = first.idx(3, 3);
    first.flags[moving_wall] = Cell::Wall;
    let mut second = first.clone();
    let mut moving_velocities = vec![[0.0; 2]; first.nx * first.ny];
    moving_velocities[moving_wall] = [0.01, -0.02];
    let mut moving_mask = vec![false; first.nx * first.ny];
    moving_mask[moving_wall] = true;
    let (mut first_scratch, mut second_scratch) = (Vec::new(), Vec::new());
    for _ in 0..8 {
        let first_receipt = first.step_with_moving_wall_momentum(
            &mut first_scratch,
            &moving_mask,
            &moving_velocities,
            [3.0, 3.0],
        );
        let second_receipt = second.step_with_moving_wall_momentum(
            &mut second_scratch,
            &moving_mask,
            &moving_velocities,
            [3.0, 3.0],
        );
        assert_eq!(first_receipt, second_receipt);
        for (first_cell, second_cell) in first.f.iter().zip(&second.f) {
            assert_eq!(first_cell.map(f64::to_bits), second_cell.map(f64::to_bits));
        }
    }
}

#[test]
fn d2q9_moving_wall_refuses_bad_fields_before_advancing() {
    let mut grid = Grid::uniform(3, 3, 0.8);
    let wall = grid.idx(1, 1);
    grid.flags[wall] = Cell::Wall;
    let original = grid.f.clone();
    let measured = vec![false; grid.nx * grid.ny];
    let mut scratch = Vec::new();

    let mut non_wall_motion = vec![[0.0; 2]; grid.nx * grid.ny];
    non_wall_motion[grid.idx(0, 0)] = [0.01, 0.0];
    let non_wall_refusal = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        let _ = grid.step_with_moving_wall_momentum(
            &mut scratch,
            &measured,
            &non_wall_motion,
            [0.0; 2],
        );
    }));
    assert!(non_wall_refusal.is_err());
    assert_eq!(grid.f, original);
    assert!(scratch.is_empty());

    let mut too_fast = vec![[0.0; 2]; grid.nx * grid.ny];
    too_fast[wall] = [0.2, 0.0];
    let speed_refusal = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        let _ = grid.step_with_moving_wall_momentum(&mut scratch, &measured, &too_fast, [0.0; 2]);
    }));
    assert!(speed_refusal.is_err());
    assert_eq!(grid.f, original);
    assert!(scratch.is_empty());

    let valid_motion = {
        let mut velocities = vec![[0.0; 2]; grid.nx * grid.ny];
        velocities[wall] = [0.01, 0.0];
        velocities
    };
    let bad_origin = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        let _ = grid.step_with_moving_wall_momentum(
            &mut scratch,
            &measured,
            &valid_motion,
            [f64::NAN, 0.0],
        );
    }));
    assert!(bad_origin.is_err());
    assert_eq!(grid.f, original);
    assert!(scratch.is_empty());

    let zero_post = vec![[0.0; Q]; grid.nx * grid.ny];
    let density_refusal = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        let _ = grid.stream_from_with_moving_wall_momentum(
            &zero_post,
            &measured,
            &valid_motion,
            [0.0; 2],
        );
    }));
    assert!(density_refusal.is_err());
    assert_eq!(grid.f, original);
}

#[test]
fn d2q9_curved_moving_wall_pins_both_linear_bfl_branches_and_receipts() {
    // G0/G3: independently enumerate the short-link and long-link linear BFL
    // rules on one east-wall link. The receipt uses the actual wall
    // intersection, not the halfway midpoint.
    let mut base = Grid::uniform(5, 3, 0.8);
    base.periodic_x = false;
    base.periodic_y = false;
    base.flags.fill(Cell::Gas);
    let far = base.idx(1, 1);
    let fluid = base.idx(2, 1);
    let wall = base.idx(3, 1);
    base.flags[far] = Cell::Fluid;
    base.flags[fluid] = Cell::Fluid;
    base.flags[wall] = Cell::Wall;

    let mut post = vec![[0.0; Q]; base.nx * base.ny];
    let outgoing = 0.25;
    let prior_opposite = 0.08;
    let far_outgoing = 0.125;
    post[fluid][0] = 0.20;
    post[fluid][1] = outgoing;
    post[fluid][3] = prior_opposite;
    post[far][1] = far_outgoing;
    let rho_post = post[fluid].iter().sum::<f64>();
    let wall_velocity = [0.03, 0.02];
    let mut wall_velocities = vec![[0.0; 2]; base.nx * base.ny];
    wall_velocities[wall] = wall_velocity;
    let mut measured = vec![false; base.nx * base.ny];
    measured[wall] = true;
    let origin = [0.4, 0.25];
    let moving_correction = 2.0 * (1.0 / 9.0) * rho_post * -wall_velocity[0] / CS2;

    let short_fraction = 0.25;
    let short_links = d2q9_curved_wall_links(&base, short_fraction, &wall_velocities);
    let mut short = base.clone();
    let short_receipt =
        short.stream_from_with_curved_moving_wall_momentum(&post, &measured, &short_links, origin);
    let short_incoming = 2.0 * short_fraction * outgoing
        + (1.0 - 2.0 * short_fraction) * far_outgoing
        + moving_correction;
    assert!((short.f[fluid][3] - short_incoming).abs() < 1e-15);

    let long_fraction = 0.75;
    let long_links = d2q9_curved_wall_links(&base, long_fraction, &wall_velocities);
    let mut long = base.clone();
    let long_receipt =
        long.stream_from_with_curved_moving_wall_momentum(&post, &measured, &long_links, origin);
    let denominator = 2.0 * long_fraction;
    let long_incoming = (outgoing + moving_correction) / denominator
        + ((denominator - 1.0) / denominator) * prior_opposite;
    assert!((long.f[fluid][3] - long_incoming).abs() < 1e-15);

    let assert_receipt = |receipt: MovingWallMomentumExchange2, incoming: f64, fraction: f64| {
        let e = [-1.0, 0.0];
        let expected_wall = [
            (-e[0] - wall_velocity[0]) * outgoing - (e[0] - wall_velocity[0]) * incoming,
            (-e[1] - wall_velocity[1]) * outgoing - (e[1] - wall_velocity[1]) * incoming,
        ];
        let expected_fluid = [e[0] * (incoming + outgoing), 0.0];
        let mass_change = incoming - outgoing;
        let expected_mass_impulse = [
            wall_velocity[0] * mass_change,
            wall_velocity[1] * mass_change,
        ];
        let link_offset = [2.0 + fraction - origin[0], 1.0 - origin[1]];
        let expected_wall_torque =
            link_offset[0] * expected_wall[1] - link_offset[1] * expected_wall[0];
        let expected_work =
            expected_wall[0] * wall_velocity[0] + expected_wall[1] * wall_velocity[1];

        assert_eq!(receipt.measured_links, 1);
        for axis in 0..2 {
            assert!((receipt.wall_impulse[axis] - expected_wall[axis]).abs() < 1e-15);
            assert!((receipt.fluid_population_impulse[axis] - expected_fluid[axis]).abs() < 1e-15);
            assert!(
                (receipt.wall_velocity_mass_impulse[axis] - expected_mass_impulse[axis]).abs()
                    < 1e-15
            );
            assert!(
                (receipt.wall_impulse[axis] + receipt.fluid_population_impulse[axis]
                    - receipt.wall_velocity_mass_impulse[axis])
                    .abs()
                    < 1e-15
            );
        }
        assert!((receipt.fluid_mass_change - mass_change).abs() < 1e-15);
        assert!((receipt.wall_angular_impulse - expected_wall_torque).abs() < 1e-15);
        assert!(
            (receipt.wall_angular_impulse + receipt.fluid_population_angular_impulse
                - receipt.wall_velocity_mass_angular_impulse)
                .abs()
                < 1e-15
        );
        assert!((receipt.wall_work - expected_work).abs() < 1e-15);
    };
    assert_receipt(short_receipt, short_incoming, short_fraction);
    assert_receipt(long_receipt, long_incoming, long_fraction);
}

#[test]
fn d2q9_curved_halfway_is_exactly_legacy_compatible_and_replays() {
    // G3/G5: q=1/2 is the exact compatibility seam, including the receipt;
    // a non-halfway collide-stream trajectory also replays bit-for-bit.
    let mut legacy = Grid::uniform(5, 5, 0.8);
    let wall = legacy.idx(2, 2);
    legacy.flags[wall] = Cell::Wall;
    let perturbed = legacy.idx(1, 2);
    legacy.f[perturbed] = equilibrium(1.1, 0.03, -0.01);
    let post = legacy.f.clone();
    let mut curved = legacy.clone();
    let mut measured = vec![false; legacy.nx * legacy.ny];
    measured[wall] = true;
    let mut wall_velocities = vec![[0.0; 2]; legacy.nx * legacy.ny];
    wall_velocities[wall] = [0.01, -0.02];
    let halfway_links = d2q9_curved_wall_links(&legacy, 0.5, &wall_velocities);

    let legacy_receipt = legacy.stream_from_with_moving_wall_momentum(
        &post,
        &measured,
        &wall_velocities,
        [1.25, -0.5],
    );
    let curved_receipt = curved.stream_from_with_curved_moving_wall_momentum(
        &post,
        &measured,
        &halfway_links,
        [1.25, -0.5],
    );
    assert_eq!(legacy_receipt, curved_receipt);
    for (legacy_cell, curved_cell) in legacy.f.iter().zip(&curved.f) {
        assert_eq!(legacy_cell.map(f64::to_bits), curved_cell.map(f64::to_bits));
    }

    let mut stationary_legacy = Grid::uniform(5, 5, 0.8);
    let stationary_wall = stationary_legacy.idx(2, 2);
    stationary_legacy.flags[stationary_wall] = Cell::Wall;
    let stationary_post = stationary_legacy.f.clone();
    let mut stationary_curved = stationary_legacy.clone();
    let stationary_velocities = vec![[0.0; 2]; stationary_legacy.nx * stationary_legacy.ny];
    let stationary_links = d2q9_curved_wall_links(&stationary_legacy, 0.5, &stationary_velocities);
    let mut stationary_measured = vec![false; stationary_legacy.nx * stationary_legacy.ny];
    stationary_measured[stationary_wall] = true;
    let stationary_legacy_receipt = stationary_legacy.stream_from_with_moving_wall_momentum(
        &stationary_post,
        &stationary_measured,
        &stationary_velocities,
        [0.0; 2],
    );
    let stationary_curved_receipt = stationary_curved.stream_from_with_curved_moving_wall_momentum(
        &stationary_post,
        &stationary_measured,
        &stationary_links,
        [0.0; 2],
    );
    assert_eq!(stationary_legacy_receipt, stationary_curved_receipt);
    for (legacy_cell, curved_cell) in stationary_legacy.f.iter().zip(&stationary_curved.f) {
        assert_eq!(legacy_cell.map(f64::to_bits), curved_cell.map(f64::to_bits));
    }

    let mut first = Grid::uniform(7, 7, 0.8);
    let replay_wall = first.idx(3, 3);
    first.flags[replay_wall] = Cell::Wall;
    let mut second = first.clone();
    let mut replay_velocities = vec![[0.0; 2]; first.nx * first.ny];
    replay_velocities[replay_wall] = [0.01, -0.015];
    let replay_links = d2q9_curved_wall_links(&first, 0.25, &replay_velocities);
    let mut replay_measured = vec![false; first.nx * first.ny];
    replay_measured[replay_wall] = true;
    let (mut first_scratch, mut second_scratch) = (Vec::new(), Vec::new());
    for _ in 0..6 {
        let first_receipt = first.step_with_curved_moving_wall_momentum(
            &mut first_scratch,
            &replay_measured,
            &replay_links,
            [3.0, 3.0],
        );
        let second_receipt = second.step_with_curved_moving_wall_momentum(
            &mut second_scratch,
            &replay_measured,
            &replay_links,
            [3.0, 3.0],
        );
        assert_eq!(first_receipt, second_receipt);
        for (first_cell, second_cell) in first.f.iter().zip(&second.f) {
            assert_eq!(first_cell.map(f64::to_bits), second_cell.map(f64::to_bits));
        }
    }
}

#[test]
fn d2q9_curved_wall_admission_is_strict_and_atomic() {
    let valid = CurvedMovingWallLink2::new(0.25, [0.01, -0.02]);
    assert_eq!(valid.fraction().to_bits(), 0.25f64.to_bits());
    assert_eq!(valid.wall_velocity(), [0.01, -0.02]);
    for bad_fraction in [0.0, -0.1, 1.1, f64::NAN] {
        assert!(
            std::panic::catch_unwind(|| {
                let _ = CurvedMovingWallLink2::new(bad_fraction, [0.0; 2]);
            })
            .is_err()
        );
    }
    assert!(
        std::panic::catch_unwind(|| {
            let _ = CurvedMovingWallLink2::new(0.5, [0.2, 0.0]);
        })
        .is_err()
    );

    let mut grid = Grid::uniform(3, 1, 0.8);
    grid.periodic_x = false;
    grid.periodic_y = false;
    grid.flags.fill(Cell::Gas);
    let far = grid.idx(0, 0);
    let fluid = grid.idx(1, 0);
    let wall = grid.idx(2, 0);
    grid.flags[fluid] = Cell::Fluid;
    grid.flags[wall] = Cell::Wall;
    let original = grid.f.clone();
    let mut measured = vec![false; grid.nx * grid.ny];
    measured[wall] = true;
    let mut scratch = Vec::new();
    let missing = vec![[None; Q]; grid.nx * grid.ny];
    let missing_refusal = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        let _ =
            grid.step_with_curved_moving_wall_momentum(&mut scratch, &measured, &missing, [0.0; 2]);
    }));
    assert!(missing_refusal.is_err());
    assert_eq!(grid.f, original);
    assert!(scratch.is_empty());

    let mut spurious = missing.clone();
    spurious[fluid][3] = Some(CurvedMovingWallLink2::new(0.5, [0.0; 2]));
    spurious[fluid][0] = Some(CurvedMovingWallLink2::new(0.5, [0.0; 2]));
    let spurious_refusal = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        let _ = grid.step_with_curved_moving_wall_momentum(
            &mut scratch,
            &measured,
            &spurious,
            [0.0; 2],
        );
    }));
    assert!(spurious_refusal.is_err());
    assert_eq!(grid.f, original);
    assert!(scratch.is_empty());

    let mut short = missing.clone();
    short[fluid][3] = Some(CurvedMovingWallLink2::new(0.25, [0.0; 2]));
    let mut post = vec![[0.0; Q]; grid.nx * grid.ny];
    post[fluid][1] = 0.25;
    let no_far_refusal = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        let _ =
            grid.stream_from_with_curved_moving_wall_momentum(&post, &measured, &short, [0.0; 2]);
    }));
    assert!(no_far_refusal.is_err());
    assert_eq!(grid.f, original);

    grid.flags[far] = Cell::Fluid;
    let original_with_far = grid.f.clone();
    post[far][1] = f64::NAN;
    let bad_far_refusal = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        let _ =
            grid.stream_from_with_curved_moving_wall_momentum(&post, &measured, &short, [0.0; 2]);
    }));
    assert!(bad_far_refusal.is_err());
    assert_eq!(grid.f, original_with_far);

    let mut long = short;
    long[fluid][3] = Some(CurvedMovingWallLink2::new(0.75, [0.0; 2]));
    post[far][1] = 0.125;
    post[fluid][3] = f64::NAN;
    let bad_local_refusal = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        let _ =
            grid.stream_from_with_curved_moving_wall_momentum(&post, &measured, &long, [0.0; 2]);
    }));
    assert!(bad_local_refusal.is_err());
    assert_eq!(grid.f, original_with_far);
}

#[test]
fn d2q9_wall_topology_transition_initializes_fresh_cells_and_closes_receipt() {
    // G0/G3: move one lattice wall cell into an adjacent fluid cell. The
    // vacated cell has seven unique surviving one-ring donors, all at the same
    // equilibrium, while the covered cell carries independently known mass
    // and momentum.
    let mut first = Grid::uniform(5, 3, 0.8);
    first.periodic_x = false;
    first.periodic_y = false;
    let vacated = first.idx(1, 1);
    let covered = first.idx(2, 1);
    first.tau.fill(0.9);
    first.fext.fill([0.02, -0.01]);
    first.flags[vacated] = Cell::Wall;
    first.f[covered] = equilibrium(1.2, 0.04, -0.01);
    let mut second = first.clone();
    let before_mass = first.total_mass();
    let before_momentum = d2q9_active_raw_momentum(&first);
    let mut next_walls = vec![false; first.nx * first.ny];
    next_walls[covered] = true;

    let first_receipt = first.transition_wall_topology(&next_walls);
    let second_receipt = second.transition_wall_topology(&next_walls);

    assert_eq!(first_receipt, second_receipt);
    assert_eq!(first.flags, second.flags);
    assert_eq!(first.f, second.f);
    assert_eq!(first.tau, second.tau);
    assert_eq!(first.fext, second.fext);
    assert_eq!(first_receipt.covered_fluid_cells, 1);
    assert_eq!(first_receipt.fresh_fluid_cells, 1);
    assert_eq!(first_receipt.fresh_donor_samples, 7);
    assert_eq!(first.flags[vacated], Cell::Fluid);
    assert_eq!(first.flags[covered], Cell::Wall);
    assert!(first.f[covered].iter().all(|population| *population == 0.0));
    assert!((first_receipt.removed_mass - 1.2).abs() < 1e-12);
    assert!((first_receipt.fresh_mass - 1.0).abs() < 1e-12);
    assert!((first_receipt.removed_momentum[0] - 1.2 * 0.04).abs() < 1e-12);
    assert!((first_receipt.removed_momentum[1] - 1.2 * -0.01).abs() < 1e-12);
    assert!(first_receipt.fresh_momentum[0].abs() < 1e-15);
    assert!(first_receipt.fresh_momentum[1].abs() < 1e-15);
    assert!((first.tau[vacated] - 0.9).abs() < 1e-15);
    assert!((first.fext[vacated][0] - 0.02).abs() < 1e-15);
    assert!((first.fext[vacated][1] + 0.01).abs() < 1e-15);

    let after_mass = first.total_mass();
    let after_momentum = d2q9_active_raw_momentum(&first);
    assert!((after_mass - before_mass - first_receipt.net_mass_change).abs() < 1e-12);
    for axis in 0..2 {
        assert!(
            (after_momentum[axis]
                - before_momentum[axis]
                - first_receipt.net_momentum_change[axis])
                .abs()
                < 1e-12
        );
        assert!(
            (first_receipt.fresh_momentum[axis]
                - first_receipt.removed_momentum[axis]
                - first_receipt.net_momentum_change[axis])
                .abs()
                < 1e-15
        );
    }

    // Reapplying the committed mask is a deterministic no-op receipt.
    let committed = first.clone();
    let no_op = first.transition_wall_topology(&next_walls);
    assert_eq!(no_op.covered_fluid_cells, 0);
    assert_eq!(no_op.fresh_fluid_cells, 0);
    assert_eq!(no_op.fresh_donor_samples, 0);
    assert_eq!(no_op.removed_mass.to_bits(), 0.0f64.to_bits());
    assert_eq!(no_op.fresh_mass.to_bits(), 0.0f64.to_bits());
    assert_eq!(first.flags, committed.flags);
    assert_eq!(first.f, committed.f);
    assert_eq!(first.tau, committed.tau);
    assert_eq!(first.fext, committed.fext);
}

#[test]
fn d2q9_wall_topology_transition_refuses_without_donors_atomically() {
    let mut grid = Grid::uniform(3, 3, 0.8);
    grid.flags.fill(Cell::Wall);
    let original_flags = grid.flags.clone();
    let original_populations = grid.f.clone();
    let original_tau = grid.tau.clone();
    let original_external_force = grid.fext.clone();
    let mut next_walls = vec![true; grid.nx * grid.ny];
    let fresh = grid.idx(1, 1);
    next_walls[fresh] = false;

    let refusal = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        let _ = grid.transition_wall_topology(&next_walls);
    }));
    assert!(refusal.is_err());
    assert_eq!(grid.flags, original_flags);
    assert_eq!(grid.f, original_populations);
    assert_eq!(grid.tau, original_tau);
    assert_eq!(grid.fext, original_external_force);

    let mut mixed = Grid::uniform(3, 3, 0.8);
    mixed.flags[0] = Cell::Gas;
    let mixed_original = mixed.clone();
    let no_walls = vec![false; mixed.nx * mixed.ny];
    let mixed_refusal = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        let _ = mixed.transition_wall_topology(&no_walls);
    }));
    assert!(mixed_refusal.is_err());
    assert_eq!(mixed.flags, mixed_original.flags);
    assert_eq!(mixed.f, mixed_original.f);
    assert_eq!(mixed.tau, mixed_original.tau);
    assert_eq!(mixed.fext, mixed_original.fext);
}

#[test]
fn d2q9_wall_momentum_refuses_invalid_masks_before_advancing() {
    let mut grid = Grid::uniform(3, 3, 0.8);
    let original = grid.f.clone();
    let mut scratch = Vec::new();
    let non_wall_mask = vec![true; grid.nx * grid.ny];
    let refusal = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        let _ = grid.step_with_wall_momentum(&mut scratch, &non_wall_mask);
    }));
    assert!(refusal.is_err());
    assert_eq!(grid.f, original);
    assert!(scratch.is_empty());

    let short_mask = vec![false; grid.nx * grid.ny - 1];
    let short_refusal = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        let _ = grid.step_with_wall_momentum(&mut scratch, &short_mask);
    }));
    assert!(short_refusal.is_err());
    assert_eq!(grid.f, original);
    assert!(scratch.is_empty());
}

#[test]
fn d2q9_regularized_x_faces_impose_moments_and_preserve_stress() {
    let boundary = VelocityPressureX2::new([0.02, 0.003], 1.0);
    assert_eq!(boundary.inlet_velocity(), [0.02, 0.003]);
    assert_eq!(boundary.outlet_density().to_bits(), 1.0f64.to_bits());

    let mut plain = Grid::uniform(9, 8, 0.8);
    plain.periodic_x = false;
    // Seed a density/momentum-conserving xx nonequilibrium on the initial
    // inlet column so the copied-stress oracle cannot pass on equilibrium
    // zeros alone.
    for y in 0..plain.ny {
        let inlet = plain.idx(0, y);
        plain.f[inlet][0] -= 0.002;
        plain.f[inlet][1] += 0.001;
        plain.f[inlet][3] += 0.001;
    }
    let mut measured = plain.clone();
    let measured_mask = vec![false; measured.nx * measured.ny];
    let (mut plain_scratch, mut measured_scratch) = (Vec::new(), Vec::new());
    for _ in 0..4 {
        plain.step_velocity_pressure_x(&mut plain_scratch, boundary);
        let receipt = measured.step_velocity_pressure_x_with_wall_momentum(
            &mut measured_scratch,
            boundary,
            &measured_mask,
        );
        assert_eq!(receipt.measured_links, 0);
        assert_eq!(receipt.wall_impulse.map(f64::to_bits), [0, 0]);
    }
    for (plain_cell, measured_cell) in plain.f.iter().zip(&measured.f) {
        assert_eq!(
            plain_cell.map(f64::to_bits),
            measured_cell.map(f64::to_bits)
        );
    }

    let mut max_stress = 0.0f64;
    for y in 0..plain.ny {
        let inlet = plain.idx(0, y);
        let inlet_source = plain.idx(1, y);
        let outlet_source = plain.idx(plain.nx - 2, y);
        let outlet = plain.idx(plain.nx - 1, y);
        let inlet_moments = plain.moments(inlet);
        let inlet_source_moments = plain.moments(inlet_source);
        let outlet_source_moments = plain.moments(outlet_source);
        let outlet_moments = plain.moments(outlet);
        assert!((inlet_moments.rho - inlet_source_moments.rho).abs() < 2e-12);
        assert!(
            inlet_moments
                .u
                .into_iter()
                .zip(boundary.inlet_velocity())
                .all(|(actual, target)| (actual - target).abs() < 2e-12)
        );
        assert!((outlet_moments.rho - boundary.outlet_density()).abs() < 2e-12);
        assert!(
            outlet_moments
                .u
                .into_iter()
                .zip(outlet_source_moments.u)
                .all(|(actual, source)| (actual - source).abs() < 2e-12)
        );

        let inlet_stress =
            d2q9_nonequilibrium_stress(&plain.f[inlet], inlet_moments.rho, inlet_moments.u);
        let inlet_source_stress = d2q9_nonequilibrium_stress(
            &plain.f[inlet_source],
            inlet_source_moments.rho,
            inlet_source_moments.u,
        );
        let outlet_stress =
            d2q9_nonequilibrium_stress(&plain.f[outlet], outlet_moments.rho, outlet_moments.u);
        let outlet_source_stress = d2q9_nonequilibrium_stress(
            &plain.f[outlet_source],
            outlet_source_moments.rho,
            outlet_source_moments.u,
        );
        for row in 0..2 {
            for column in 0..2 {
                max_stress = max_stress.max(inlet_source_stress[row][column].abs());
                assert!(
                    (inlet_stress[row][column] - inlet_source_stress[row][column]).abs() < 2e-12
                );
                assert!(
                    (outlet_stress[row][column] - outlet_source_stress[row][column]).abs() < 2e-12
                );
            }
        }
    }
    assert!(max_stress > 1e-10, "stress-copy oracle must be non-vacuous");
}

#[test]
fn d2q9_regularized_x_faces_refuse_invalid_setup_before_advancing() {
    assert!(std::panic::catch_unwind(|| VelocityPressureX2::new([0.2, 0.0], 1.0)).is_err());
    assert!(std::panic::catch_unwind(|| VelocityPressureX2::new([0.02, 0.0], 0.0)).is_err());

    let boundary = VelocityPressureX2::new([0.02, 0.0], 1.0);
    let mut grid = Grid::uniform(5, 4, 0.8);
    let original = grid.f.clone();
    let mut scratch = Vec::new();
    let periodic_refusal = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        grid.step_velocity_pressure_x(&mut scratch, boundary);
    }));
    assert!(periodic_refusal.is_err());
    assert_eq!(grid.f, original);
    assert!(scratch.is_empty());

    grid.periodic_x = false;
    grid.g = [1e-6, 0.0];
    let force_refusal = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        grid.step_velocity_pressure_x(&mut scratch, boundary);
    }));
    assert!(force_refusal.is_err());
    assert_eq!(grid.f, original);
    assert!(scratch.is_empty());

    grid.g = [0.0; 2];
    let blocked_interior = grid.idx(1, 2);
    grid.flags[blocked_interior] = Cell::Wall;
    let topology_refusal = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        grid.step_velocity_pressure_x(&mut scratch, boundary);
    }));
    assert!(topology_refusal.is_err());
    assert_eq!(grid.f, original);
    assert!(scratch.is_empty());
}
