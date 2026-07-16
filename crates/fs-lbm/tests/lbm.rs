//! Battery for the D2Q9 lattice Boltzmann core (fs-lbm). Covers the
//! equilibrium moments, mass conservation, the load-bearing physical check
//! (steady Poiseuille channel flow matches the analytic parabola), and the
//! lattice-scaling assistant's stability bookkeeping.

use fs_lbm::core2::VelocityPressureX2;
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
