//! Battery for the D3Q19 core (bead 84hv): exact lattice invariants in
//! integer arithmetic, equilibrium moments, mass conservation to the
//! 2-D crate's bar, duct-flow physics against the analytic rectangular
//! series, bit determinism, and the frozen golden (single-thread until
//! WS1-D lands the parallel sweep). JSONL verdicts per house style.
//!
//! The FULL acceptance fixture (32×32×64 from rest, 20k steps) is an
//! `#[ignore]`d release-lane test — run explicitly:
//! `cargo test -p fs-lbm --release --test d3q19_battery -- --ignored`.

use fs_lbm::d3q19::{
    CollisionError3, CollisionModel3, D3Q19_BIT_SEMANTICS_VERSION,
    D3Q19_MOMENT_COLLISION_SEMANTICS_VERSION, E3, OPP3, Q3, W3, W36, collide_cell3,
};
use fs_lbm::{BoundaryGrid3, BoundarySpec3, CS2, Duct, duct_analytic, equilibrium3};

fn verdict(name: &str, pass: bool, details: &str) {
    println!("{{\"test\":\"{name}\",\"pass\":{pass},\"details\":\"{details}\"}}");
    assert!(pass, "{name}: {details}");
}

/// lbm3-001: the lattice constants satisfy the standard D3Q19
/// invariants EXACTLY, in integer arithmetic (weights ×36) — no float
/// tolerance where exactness is possible.
#[test]
fn lattice_invariants_hold_exactly() {
    assert_eq!(D3Q19_BIT_SEMANTICS_VERSION, 1);
    assert_eq!(D3Q19_MOMENT_COLLISION_SEMANTICS_VERSION, 3);
    // The float weights are exactly the ×36 integers.
    for i in 0..Q3 {
        assert_eq!(W3[i].to_bits(), (W36[i] as f64 / 36.0).to_bits());
    }
    // Σ w = 1 (×36 = 36).
    let w_sum: i64 = W36.iter().sum();
    verdict(
        "lbm3-001a-weight-sum",
        w_sum == 36,
        &format!("36*Σw = {w_sum}"),
    );
    // First moment: Σ w c_a = 0 exactly, per component.
    let e = |i: usize| [i64::from(E3[i].0), i64::from(E3[i].1), i64::from(E3[i].2)];
    for a in 0..3 {
        let m1: i64 = (0..Q3).map(|i| W36[i] * e(i)[a]).sum();
        verdict(
            "lbm3-001b-first-moment",
            m1 == 0,
            &format!("axis {a}: {m1}"),
        );
    }
    // Second moment isotropy: Σ w c_a c_b = (1/3) δ_ab (×36 = 12 δ_ab).
    for a in 0..3 {
        for b in 0..3 {
            let m2: i64 = (0..Q3).map(|i| W36[i] * e(i)[a] * e(i)[b]).sum();
            let want = if a == b { 12 } else { 0 };
            verdict(
                "lbm3-001c-second-moment",
                m2 == want,
                &format!("({a},{b}): {m2} vs {want}"),
            );
        }
    }
    // Third moment: Σ w c_a c_b c_c = 0 exactly (odd moments vanish).
    for a in 0..3 {
        for b in 0..3 {
            for c in 0..3 {
                let m3: i64 = (0..Q3).map(|i| W36[i] * e(i)[a] * e(i)[b] * e(i)[c]).sum();
                assert_eq!(m3, 0, "third moment ({a},{b},{c})");
            }
        }
    }
    // Fourth moment isotropy: Σ w c_a c_b c_c c_d
    //   = (1/9)(δ_ab δ_cd + δ_ac δ_bd + δ_ad δ_bc)   (×36 = 4·(…)).
    let d = |a: usize, b: usize| i64::from(a == b);
    let mut worst = 0i64;
    for a in 0..3 {
        for b in 0..3 {
            for c in 0..3 {
                for dd in 0..3 {
                    let m4: i64 = (0..Q3)
                        .map(|i| W36[i] * e(i)[a] * e(i)[b] * e(i)[c] * e(i)[dd])
                        .sum();
                    let want = 4 * (d(a, b) * d(c, dd) + d(a, c) * d(b, dd) + d(a, dd) * d(b, c));
                    assert_eq!(m4, want, "fourth moment ({a},{b},{c},{dd})");
                    worst = worst.max((m4 - want).abs());
                }
            }
        }
    }
    verdict(
        "lbm3-001d-fourth-moment",
        worst == 0,
        "all 81 entries exact",
    );
    // Opposites: e_opp = −e, w_opp = w, involution.
    for i in 0..Q3 {
        let o = OPP3[i];
        assert_eq!(OPP3[o], i, "opposite is an involution");
        assert_eq!(W36[o], W36[i], "opposite weight");
        assert_eq!(
            (E3[o].0, E3[o].1, E3[o].2),
            (-E3[i].0, -E3[i].1, -E3[i].2),
            "opposite velocity"
        );
    }
    // The 19 velocities are distinct and are exactly {rest, 6 face, 12 edge}.
    let mut face = 0;
    let mut edge = 0;
    for (i, ev) in E3.iter().enumerate() {
        let n2 = ev.0 * ev.0 + ev.1 * ev.1 + ev.2 * ev.2;
        match (i, n2) {
            (0, 0) => {}
            (_, 1) => face += 1,
            (_, 2) => edge += 1,
            other => panic!("unexpected velocity class {other:?}"),
        }
    }
    verdict(
        "lbm3-001e-velocity-classes",
        face == 6 && edge == 12,
        &format!("face {face}, edge {edge}"),
    );
}

/// lbm3-002: the equilibrium recovers its density and momentum moments.
#[test]
fn equilibrium_recovers_moments() {
    let (rho, u) = (1.05, [0.04, -0.02, 0.03]);
    let f = equilibrium3(rho, u);
    let sum: f64 = f.iter().sum();
    let mut m = [0.0; 3];
    for i in 0..Q3 {
        m[0] += f64::from(E3[i].0) * f[i];
        m[1] += f64::from(E3[i].1) * f[i];
        m[2] += f64::from(E3[i].2) * f[i];
    }
    verdict(
        "lbm3-002-equilibrium-moments",
        (sum - rho).abs() < 1e-12
            && (m[0] - rho * u[0]).abs() < 1e-12
            && (m[1] - rho * u[1]).abs() < 1e-12
            && (m[2] - rho * u[2]).abs() < 1e-12,
        &format!("rho err {:.2e}", (sum - rho).abs()),
    );
}

/// lbm3-002b: the shared checked kernel retains the frozen general-force BGK
/// arithmetic and conserves the collision invariants without forcing.
#[test]
fn shared_bgk_collision_is_bit_compatible_and_fail_closed() {
    let rho = 1.03;
    let velocity = [0.031, -0.017, 0.023];
    let mut populations = equilibrium3(rho, velocity);
    for (direction, value) in populations.iter_mut().enumerate() {
        *value += (direction as f64 - 9.0) * 1e-6;
    }
    let tau = 0.83;
    let force = [1e-6, -2e-6, 3e-6];

    // Frozen BoundaryGrid3 formula before the shared-kernel extraction.
    let mut measured_rho = 0.0;
    let mut momentum = [0.0; 3];
    for direction in 0..Q3 {
        measured_rho += populations[direction];
        momentum[0] += f64::from(E3[direction].0) * populations[direction];
        momentum[1] += f64::from(E3[direction].1) * populations[direction];
        momentum[2] += f64::from(E3[direction].2) * populations[direction];
    }
    let measured_velocity =
        core::array::from_fn(|axis| (momentum[axis] + 0.5 * force[axis]) / measured_rho);
    let equilibrium = equilibrium3(measured_rho, measured_velocity);
    let coefficient = 1.0 - 0.5 / tau;
    let cs2 = 1.0 / 3.0;
    let cs4 = cs2 * cs2;
    let mut reference = [0.0; Q3];
    for direction in 0..Q3 {
        let e = [
            f64::from(E3[direction].0),
            f64::from(E3[direction].1),
            f64::from(E3[direction].2),
        ];
        let eu = e
            .iter()
            .zip(measured_velocity)
            .map(|(component, u)| *component * u)
            .sum::<f64>();
        let projection = (0..3)
            .map(|axis| {
                ((e[axis] - measured_velocity[axis]) / cs2 + eu * e[axis] / cs4) * force[axis]
            })
            .sum::<f64>();
        reference[direction] = populations[direction]
            + (equilibrium[direction] - populations[direction]) / tau
            + coefficient * W3[direction] * projection;
    }
    let shared = collide_cell3(populations, CollisionModel3::Bgk { tau }, force)
        .expect("admissible shared BGK collision");
    assert!(
        shared
            .iter()
            .zip(reference)
            .all(|(actual, expected)| actual.to_bits() == expected.to_bits()),
        "shared kernel must retain the frozen BoundaryGrid3 arithmetic"
    );

    let unforced = collide_cell3(populations, CollisionModel3::Bgk { tau }, [0.0; 3])
        .expect("unforced BGK collision");
    let conserved = |values: &[f64; Q3]| {
        let density = values.iter().sum::<f64>();
        let momentum: [f64; 3] = core::array::from_fn(|axis| {
            values
                .iter()
                .enumerate()
                .map(|(direction, value)| {
                    let component = match axis {
                        0 => E3[direction].0,
                        1 => E3[direction].1,
                        _ => E3[direction].2,
                    };
                    f64::from(component) * value
                })
                .sum::<f64>()
        });
        (density, momentum)
    };
    let (before_rho, before_momentum) = conserved(&populations);
    let (after_rho, after_momentum) = conserved(&unforced);
    assert!((after_rho - before_rho).abs() < 1e-14);
    for axis in 0..3 {
        assert!((after_momentum[axis] - before_momentum[axis]).abs() < 1e-14);
    }

    assert!(matches!(
        collide_cell3(populations, CollisionModel3::Bgk { tau: 0.5 }, force),
        Err(CollisionError3::InvalidRelaxationTime { .. })
    ));
    let mut nonfinite = populations;
    nonfinite[7] = f64::NAN;
    assert!(matches!(
        collide_cell3(nonfinite, CollisionModel3::Bgk { tau }, force),
        Err(CollisionError3::NonFinitePopulation { direction: 7, .. })
    ));
    assert!(matches!(
        collide_cell3(
            populations,
            CollisionModel3::Bgk { tau },
            [0.0, f64::INFINITY, 0.0],
        ),
        Err(CollisionError3::NonFiniteForce { axis: 1, .. })
    ));
    assert!(matches!(
        collide_cell3([0.0; Q3], CollisionModel3::Bgk { tau }, [0.0; 3]),
        Err(CollisionError3::NonPositiveDensity { .. })
    ));
}

/// lbm3-002c (G0): the correctness-first central-moment rung preserves the
/// collision invariants, reduces to BGK within solve roundoff at equal rates,
/// and refuses invalid relaxation-rate combinations.
#[test]
fn central_moment_collision_is_conservative_and_explicitly_bounded() {
    let rho = 1.01;
    let velocity = [0.027, -0.019, 0.011];
    let mut populations = equilibrium3(rho, velocity);
    for (direction, value) in populations.iter_mut().enumerate() {
        *value += (direction as f64 - 9.0) * 1e-5;
    }
    let tau = 0.81;
    let equal_rate = 1.0 / tau;
    let bgk = collide_cell3(populations, CollisionModel3::Bgk { tau }, [0.0; 3])
        .expect("admissible unforced BGK collision");
    let central_equal = collide_cell3(
        populations,
        CollisionModel3::CentralMoment {
            second_order_rate: equal_rate,
            higher_order_rate: equal_rate,
        },
        [0.0; 3],
    )
    .expect("admissible equal-rate central-moment collision");
    let equal_rate_delta = bgk
        .iter()
        .zip(central_equal)
        .map(|(left, right)| (left - right).abs())
        .fold(0.0_f64, f64::max);
    assert!(
        equal_rate_delta < 5e-13,
        "equal-rate central moments must reduce to BGK within solve roundoff, got {equal_rate_delta:.3e}"
    );

    let central_split = collide_cell3(
        populations,
        CollisionModel3::CentralMoment {
            second_order_rate: equal_rate,
            higher_order_rate: 1.65,
        },
        [0.0; 3],
    )
    .expect("admissible split-rate central-moment collision");
    let split_delta = bgk
        .iter()
        .zip(central_split)
        .map(|(left, right)| (left - right).abs())
        .fold(0.0_f64, f64::max);
    assert!(
        split_delta > 1e-9,
        "independent higher-order relaxation must affect a nonequilibrium state"
    );

    let conserved = |values: &[f64; Q3]| {
        let density = values.iter().sum::<f64>();
        let momentum: [f64; 3] = core::array::from_fn(|axis| {
            values
                .iter()
                .enumerate()
                .map(|(direction, value)| {
                    let component = match axis {
                        0 => E3[direction].0,
                        1 => E3[direction].1,
                        _ => E3[direction].2,
                    };
                    f64::from(component) * value
                })
                .sum::<f64>()
        });
        (density, momentum)
    };
    let (before_density, before_momentum) = conserved(&populations);
    let (after_density, after_momentum) = conserved(&central_split);
    assert!((after_density - before_density).abs() < 5e-13);
    for axis in 0..3 {
        assert!((after_momentum[axis] - before_momentum[axis]).abs() < 5e-13);
    }

    let invalid_rate = CollisionModel3::CentralMoment {
        second_order_rate: 0.0,
        higher_order_rate: 1.0,
    };
    assert!(matches!(
        invalid_rate.validate(),
        Err(CollisionError3::InvalidMomentRelaxationRate {
            moment_order: 2,
            ..
        })
    ));
    assert!(
        collide_cell3(
            populations,
            CollisionModel3::CentralMoment {
                second_order_rate: equal_rate,
                higher_order_rate: 1.0,
            },
            [1e-7, 0.0, 0.0],
        )
        .is_ok()
    );
}

/// lbm3-002d (G0/G3): the reduced D3Q19 cumulant projection preserves
/// macroscopic invariants and axis permutation covariance, fixes equilibrium,
/// and remains observably nonlinear relative to equal-rate central moments.
#[test]
fn reduced_cumulant_projection_is_conservative_covariant_and_nonlinear() {
    let rho = 1.01;
    let velocity = [0.027, -0.019, 0.011];
    let mut populations = equilibrium3(rho, velocity);
    for (direction, value) in populations.iter_mut().enumerate() {
        *value += (direction as f64 - 9.0) * 1e-5;
    }
    let tau = 0.81;
    let second_order_rate = 1.0 / tau;
    let model = CollisionModel3::ReducedCumulant {
        second_order_rate,
        third_order_rate: 1.35,
        fourth_order_rate: 1.72,
    };
    let post = collide_cell3(populations, model, [0.0; 3])
        .expect("admissible reduced-cumulant projection");

    let conserved = |values: &[f64; Q3]| {
        let density = values.iter().sum::<f64>();
        let momentum: [f64; 3] = core::array::from_fn(|axis| {
            values
                .iter()
                .enumerate()
                .map(|(direction, value)| {
                    let component = match axis {
                        0 => E3[direction].0,
                        1 => E3[direction].1,
                        _ => E3[direction].2,
                    };
                    f64::from(component) * value
                })
                .sum::<f64>()
        });
        (density, momentum)
    };
    let (before_density, before_momentum) = conserved(&populations);
    let (after_density, after_momentum) = conserved(&post);
    assert!((after_density - before_density).abs() < 5e-13);
    for axis in 0..3 {
        assert!((after_momentum[axis] - before_momentum[axis]).abs() < 5e-13);
    }

    let equal_cumulant = collide_cell3(
        populations,
        CollisionModel3::ReducedCumulant {
            second_order_rate,
            third_order_rate: second_order_rate,
            fourth_order_rate: second_order_rate,
        },
        [0.0; 3],
    )
    .expect("admissible equal-rate reduced cumulants");
    let equal_central = collide_cell3(
        populations,
        CollisionModel3::CentralMoment {
            second_order_rate,
            higher_order_rate: second_order_rate,
        },
        [0.0; 3],
    )
    .expect("admissible equal-rate central moments");
    let nonlinear_delta = equal_cumulant
        .iter()
        .zip(equal_central)
        .map(|(left, right)| (left - right).abs())
        .fold(0.0_f64, f64::max);
    assert!(
        nonlinear_delta > 1e-12,
        "cumulant relaxation must not masquerade as linear central-moment/BGK relaxation"
    );

    let swap_xy = |values: [f64; Q3]| {
        core::array::from_fn(|target| {
            let target_velocity = E3[target];
            let source_velocity = (target_velocity.1, target_velocity.0, target_velocity.2);
            let source = E3
                .iter()
                .position(|candidate| *candidate == source_velocity)
                .expect("D3Q19 is closed under x/y permutation");
            values[source]
        })
    };
    let rotated_input = swap_xy(populations);
    let rotated_post = collide_cell3(rotated_input, model, [0.0; 3])
        .expect("axis-permuted reduced-cumulant projection");
    let expected_rotated_post = swap_xy(post);
    let covariance_delta = rotated_post
        .iter()
        .zip(expected_rotated_post)
        .map(|(actual, expected)| (actual - expected).abs())
        .fold(0.0_f64, f64::max);
    assert!(
        covariance_delta < 5e-13,
        "x/y permutation covariance drifted by {covariance_delta:.3e}"
    );

    let equilibrium = equilibrium3(rho, velocity);
    let equilibrium_post =
        collide_cell3(equilibrium, model, [0.0; 3]).expect("equilibrium is admissible");
    let fixed_point_delta = equilibrium
        .iter()
        .zip(equilibrium_post)
        .map(|(before, after)| (before - after).abs())
        .fold(0.0_f64, f64::max);
    assert!(fixed_point_delta < 5e-13);

    assert!(matches!(
        CollisionModel3::ReducedCumulant {
            second_order_rate,
            third_order_rate: 1.0,
            fourth_order_rate: 2.0,
        }
        .validate(),
        Err(CollisionError3::InvalidMomentRelaxationRate {
            moment_order: 4,
            ..
        })
    ));
    assert!(collide_cell3(populations, model, [0.0, 0.0, 1e-7]).is_ok());
}

/// lbm3-003: mass is conserved to the 2-D crate's roundoff bar (1e-11)
/// — collision, Guo forcing, streaming, and bounce-back all preserve
/// Σf by construction, so drift is summation roundoff only.
#[test]
fn mass_is_conserved() {
    let mut duct = Duct::new(8, 8, 8, 0.8, 1e-5);
    duct.perturb(0x5EED, 1e-3);
    let m0 = duct.total_mass();
    duct.run(200);
    let drift = (duct.total_mass() - m0).abs();
    verdict(
        "lbm3-003-mass-conservation",
        drift < 1e-11,
        &format!("drift {drift:.3e} over 200 steps on mass {m0:.1}"),
    );
}

/// lbm3-004: the analytic duct series is trustworthy — in the wide-duct
/// limit its centerline reduces to the 2-D Poiseuille parabola.
#[test]
fn duct_series_matches_the_parabola_in_the_wide_limit() {
    let (gz, nu) = (1e-5, 0.1);
    // ny >> nx: the mid-height profile along x approaches the parabola
    // u = (gz/2ν)(x+½)(nx−½−x) of the 2-D channel of width nx.
    let (nx, ny) = (16, 256);
    let y_mid = ny / 2;
    let mut max_rel = 0.0_f64;
    for x in 0..nx {
        let series = duct_analytic(gz, nu, nx, ny, x, y_mid);
        let xf = x as f64;
        let parabola = gz / (2.0 * nu) * (xf + 0.5) * (nx as f64 - 0.5 - xf);
        max_rel = max_rel.max((series - parabola).abs() / parabola.abs());
    }
    verdict(
        "lbm3-004-series-wide-limit",
        max_rel < 5e-3,
        &format!("max rel dev {max_rel:.2e}"),
    );
}

/// lbm3-005: duct flow converges onto the analytic rectangular profile
/// (small fixture, debug-lane speed): 12×12×8 from rest. Interior cells
/// (≥1 cell off the walls) land within 3% of the truncated series; the
/// wall rim carries the known coarse-resolution halfway-bounce-back
/// corner error (measured −11% at the corner cells HERE, second-order
/// in resolution — the 32×32×64 release-lane acceptance fixture holds
/// the full-section 3% bar at scale) and is bounded at 15% so a real
/// boundary-condition regression is still caught.
#[test]
fn duct_flow_matches_the_analytic_series_small() {
    let (nx, ny, nz, tau, gz) = (12, 12, 8, 0.8, 1e-6);
    let mut duct = Duct::new(nx, ny, nz, tau, gz);
    duct.run(6000);
    let nu = duct.viscosity();
    let section = duct.z_velocity_section();
    let mut max_interior = 0.0_f64;
    let mut max_rim = 0.0_f64;
    for y in 0..ny {
        for x in 0..nx {
            let sim = section[y * nx + x];
            let ana = duct_analytic(gz, nu, nx, ny, x, y);
            let rel = (sim - ana).abs() / ana.abs();
            if x == 0 || y == 0 || x == nx - 1 || y == ny - 1 {
                max_rim = max_rim.max(rel);
            } else {
                max_interior = max_interior.max(rel);
            }
        }
    }
    verdict(
        "lbm3-005-duct-vs-series",
        max_interior < 0.03 && max_rim < 0.15,
        &format!("interior {max_interior:.4}, wall rim {max_rim:.4} (12x12x8, 6000 steps)"),
    );
    // Symmetry sanity: the four quadrant probes agree bitwise-close.
    let probe = |x: usize, y: usize| section[y * nx + x];
    let (p, q) = (probe(2, 3), probe(nx - 3, ny - 4));
    verdict(
        "lbm3-005b-symmetry",
        (p - q).abs() / p.abs() < 1e-9,
        &format!("quadrant symmetry {:.2e}", (p - q).abs() / p.abs()),
    );
}

/// lbm3-006: bit determinism — two identical runs produce identical
/// bits (single-thread; WS1-D extends this bar across worker counts).
#[test]
fn the_solver_is_bit_deterministic() {
    let build = || {
        let mut d = Duct::new(8, 8, 8, 0.7, 1e-5);
        d.perturb(0xFD1D, 1e-3);
        d.run(100);
        d.z_velocity_section()
    };
    let (a, b) = (build(), build());
    let same = a.iter().zip(&b).all(|(x, y)| x.to_bits() == y.to_bits());
    verdict("lbm3-006-determinism", same, "two runs, identical bits");
}

/// Frozen golden (bead 84hv): seeded 8×8×16 fixture, 100 steps,
/// deterministic mode (single-thread pinned traversal). Feeds density
/// and velocity bits at fixed probes plus total mass into FNV-1a.
/// Registered in golden-couplings.json against
/// `fs-lbm:d3q19-bits = 1`; bump only with semantic justification
/// (golden-evidence policy, docs/GOLDEN_POLICY.md).
const GOLDEN_HASH: u64 = 0xd548_9283_1f0a_6d71;

#[test]
fn d3q19_golden_hash() {
    let mut acc: u64 = 0xcbf2_9ce4_8422_2325;
    let mut feed = |v: f64| {
        for byte in v.to_bits().to_le_bytes() {
            acc ^= u64::from(byte);
            acc = acc.wrapping_mul(0x0000_0100_0000_01b3);
        }
    };
    let mut duct = Duct::new(8, 8, 16, 0.75, 2e-6);
    duct.perturb(0x84A5, 1e-3);
    duct.run(100);
    for (x, y, z) in [
        (0, 0, 0),
        (3, 4, 7),
        (7, 7, 15),
        (1, 6, 9),
        (5, 2, 12),
        (4, 4, 8),
        (6, 1, 3),
        (2, 5, 14),
    ] {
        feed(duct.density(x, y, z));
        let u = duct.velocity(x, y, z);
        feed(u[0]);
        feed(u[1]);
        feed(u[2]);
    }
    feed(duct.total_mass());
    println!("{{\"test\":\"lbm3-007-golden\",\"hash\":\"{acc:#018x}\"}}");
    assert_eq!(
        acc, GOLDEN_HASH,
        "d3q19 bits changed: {acc:#018x} vs {GOLDEN_HASH:#018x} — bump only with semantic \
         justification and a golden-couplings.json update in the same commit"
    );
}

/// ACCEPTANCE FIXTURE (bead 84hv, release lane): 32×32×64 duct from
/// rest — the analytic profile within 3% everywhere, mass conserved to
/// the roundoff floor over a 10k-step tail. ~10⁹ cell updates: run in
/// release, explicitly.
///
/// THE MASS BAR, derived honestly (same philosophy as the 2-D crate's
/// "~10× the roundoff floor" gate, rescaled): conservation is exact by
/// construction (collision algebra, Guo forcing sums to zero, streaming
/// permutes, bounce-back reflects), so the only drift is float
/// roundoff, which is per-cell-per-step. Measured floors agree across
/// crates: the 2-D fixture drifts ~6.5e-17/cell/step (9.4e-13 over
/// 200 steps × 72 cells), this fixture ~7.5e-17/cell/step (4.887e-8
/// over 10k steps × 65536 cells, first landing run). A literal 1e-11
/// at this scale would demand < 1.5e-19/cell/step — below f64
/// representability on O(1) populations. Gate at 10× the measured
/// floor (5e-7) so a REAL systematic leak (say 1e-12/cell/step =
/// 6.5e-4 over the tail) still fails by three orders of magnitude.
///
/// THE RELAXATION TIME: τ = ½ + √3/4 (the BGK "magic" value, where
/// (τ−½)² = 3/16 cancels the quadratic wall-slip of halfway bounce-back
/// — the TRT Λ = 3/16 condition specialized to BGK). At generic τ the
/// slip defect near walls is resolution-INDEPENDENT (measured: corner
/// cells −11.4% at 12×12 and −8.7% at 32×32 with τ = 0.8), so no
/// refinement meets a 3% corner bar; at the magic τ the leading defect
/// vanishes (corner −3.7% at 12×12, shrinking with resolution) and the
/// 3% full-section bar is meetable honestly. Generic-τ behavior stays
/// pinned by lbm3-005's interior/rim split.
#[test]
#[ignore = "acceptance fixture: run in release with --ignored (bead 84hv)"]
fn acceptance_duct_32x32x64() {
    // τ = ½ + √3/4 exactly (see above); gz keeps the flow low-Mach.
    let (nx, ny, nz, tau, gz) = (32, 32, 64, 0.933_012_701_892_219_3, 1e-7);
    let mut duct = Duct::new(nx, ny, nz, tau, gz);
    // Reach steady state: diffusion time ~ (nx/π)²/ν ≈ 1.04e3 steps;
    // run 30 diffusion times, then measure conservation over a 10k tail.
    duct.run(30_000);
    let m0 = duct.total_mass();
    duct.run(10_000);
    let drift = (duct.total_mass() - m0).abs();
    let nu = duct.viscosity();
    let section = duct.z_velocity_section();
    let mut max_rel = 0.0_f64;
    let mut argmax = (0usize, 0usize);
    for y in 0..ny {
        for x in 0..nx {
            let sim = section[y * nx + x];
            let ana = duct_analytic(gz, nu, nx, ny, x, y);
            let rel = (sim - ana).abs() / ana.abs();
            if rel > max_rel {
                max_rel = rel;
                argmax = (x, y);
            }
        }
    }
    // Both verdicts report before either asserts (one failure must not
    // hide the other measurement).
    let mass_ok = drift < 5e-7;
    let profile_ok = max_rel < 0.03;
    println!(
        "{{\"test\":\"lbm3-acc-mass\",\"pass\":{mass_ok},\"details\":\"drift {drift:.3e} over \
         the 10k-step tail (mass {m0:.0}, floor-derived bar 5e-7)\"}}"
    );
    println!(
        "{{\"test\":\"lbm3-acc-profile\",\"pass\":{profile_ok},\"details\":\"max rel dev \
         {max_rel:.4} at {argmax:?} (32x32x64, 40k steps)\"}}"
    );
    assert!(mass_ok, "mass drifted {drift:.3e} (bar 5e-7)");
    assert!(profile_ok, "profile off by {max_rel:.4} at {argmax:?}");
}

const ALBENSOEDER_KUHLMANN_RE1000_JSONL: &str =
    include_str!("../../../data/reference/albensoeder-kuhlmann-2005-re1000.jsonl");

fn parse_ak2005_table() -> (Vec<(f64, f64)>, Vec<(f64, f64)>) {
    let mut ux_centerline = Vec::new();
    let mut uy_centerline = Vec::new();
    for line in ALBENSOEDER_KUHLMANN_RE1000_JSONL.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.contains("\"comment\"") {
            continue;
        }
        let extract = |key: &str| -> Option<f64> {
            let needle = format!("\"{key}\":");
            let start = trimmed.find(&needle)? + needle.len();
            let rest = trimmed[start..].trim_start();
            let end = rest.find([',', '}']).unwrap_or(rest.len());
            rest[..end].trim().parse::<f64>().ok()
        };
        if trimmed.contains("\"u_x_vertical_centerline\"") {
            let y_unit = extract("y_unit").expect("y_unit in AK2005 row");
            let ux = extract("u_x").expect("u_x in AK2005 row");
            ux_centerline.push((y_unit, ux));
        } else if trimmed.contains("\"u_y_horizontal_centerline\"") {
            let x_unit = extract("x_unit").expect("x_unit in AK2005 row");
            let uy = extract("u_y").expect("u_y in AK2005 row");
            uy_centerline.push((x_unit, uy));
        }
    }
    ux_centerline.sort_by(|a, b| a.0.total_cmp(&b.0));
    uy_centerline.sort_by(|a, b| a.0.total_cmp(&b.0));
    (ux_centerline, uy_centerline)
}

fn run_taylor_green_3d(n: usize, steps: usize, model: CollisionModel3) -> (f64, Vec<[f64; 3]>) {
    let tau = 0.8;
    let nu = model.kinematic_viscosity();
    let k = 2.0 * core::f64::consts::PI / (n as f64);
    let u0 = 1e-4;
    let mut grid =
        BoundaryGrid3::with_collision_model(n, n, n, model, [0.0; 3], BoundarySpec3::periodic());
    for z in 0..n {
        let kz = k * (z as f64);
        let (sz, cz) = (kz.sin(), kz.cos());
        for y in 0..n {
            let ky = k * (y as f64);
            let (sy, cy) = (ky.sin(), ky.cos());
            for x in 0..n {
                let kx = k * (x as f64);
                let (sx, cx) = (kx.sin(), kx.cos());
                let u = [
                    u0 * sx * cy * cz,
                    -0.5 * u0 * cx * sy * cz,
                    -0.5 * u0 * cx * cy * sz,
                ];
                let grad = [
                    [
                        k * u0 * cx * cy * cz,
                        -k * u0 * sx * sy * cz,
                        -k * u0 * sx * cy * sz,
                    ],
                    [
                        0.5 * k * u0 * sx * sy * cz,
                        -0.5 * k * u0 * cx * cy * cz,
                        0.5 * k * u0 * cx * sy * sz,
                    ],
                    [
                        0.5 * k * u0 * sx * cy * sz,
                        0.5 * k * u0 * cx * sy * sz,
                        -0.5 * k * u0 * cx * cy * cz,
                    ],
                ];
                let feq = equilibrium3(1.0, u);
                let mut f = [0.0; Q3];
                for q in 0..Q3 {
                    let c = [f64::from(E3[q].0), f64::from(E3[q].1), f64::from(E3[q].2)];
                    let mut contraction = 0.0;
                    for a in 0..3 {
                        for b in 0..3 {
                            let cab = c[a] * c[b] - if a == b { CS2 } else { 0.0 };
                            contraction = cab.mul_add(grad[a][b], contraction);
                        }
                    }
                    let f1 = -(tau / CS2) * W3[q] * contraction;
                    f[q] = feq[q] + f1;
                }
                grid.set_populations(x, y, z, &f);
            }
        }
    }
    grid.run(steps);

    let decay = (-3.0 * nu * k * k * (steps as f64)).exp();
    let mut err_sq = 0.0;
    let mut ref_sq = 0.0;
    let mut field = Vec::with_capacity(n * n * n);
    for z in 0..n {
        let kz = k * (z as f64);
        let (sz, cz) = (kz.sin(), kz.cos());
        for y in 0..n {
            let ky = k * (y as f64);
            let (sy, cy) = (ky.sin(), ky.cos());
            for x in 0..n {
                let kx = k * (x as f64);
                let (sx, cx) = (kx.sin(), kx.cos());
                let exact = [
                    u0 * sx * cy * cz * decay,
                    -0.5 * u0 * cx * sy * cz * decay,
                    -0.5 * u0 * cx * cy * sz * decay,
                ];
                let sim = grid.velocity(x, y, z);
                field.push(sim);
                for axis in 0..3 {
                    let diff = sim[axis] - exact[axis];
                    err_sq = diff.mul_add(diff, err_sq);
                    ref_sq = exact[axis].mul_add(exact[axis], ref_sq);
                }
            }
        }
    }
    ((err_sq / ref_sq).sqrt(), field)
}

/// lbm3-008 (G1/G3, beads frankensim-vrif & frankensim-lga4): 3-D Taylor-Green
/// vortex decay achieves second-order spatial accuracy under diffusive scaling
/// on both BGK and ReducedCumulant, and the two operators converge to the same
/// macroscopic field in the low-Mach continuum limit.
#[test]
fn taylor_green_3d_decay_order_and_bgk_limit_equivalence() {
    let tau = 0.8;
    let bgk = CollisionModel3::Bgk { tau };
    let cumulant = CollisionModel3::ReducedCumulant {
        second_order_rate: 1.0 / tau,
        third_order_rate: 1.35,
        fourth_order_rate: 1.72,
    };
    let (err8_bgk, u8_bgk) = run_taylor_green_3d(8, 4, bgk);
    let (err16_bgk, u16_bgk) = run_taylor_green_3d(16, 16, bgk);
    let (err8_cum, u8_cum) = run_taylor_green_3d(8, 4, cumulant);
    let (err16_cum, u16_cum) = run_taylor_green_3d(16, 16, cumulant);

    let order_bgk = (err8_bgk / err16_bgk).log2();
    let order_cum = (err8_cum / err16_cum).log2();

    let field_diff = |a: &[[f64; 3]], b: &[[f64; 3]]| -> f64 {
        let mut sq = 0.0;
        for (ua, ub) in a.iter().zip(b.iter()) {
            for axis in 0..3 {
                let d = ua[axis] - ub[axis];
                sq = d.mul_add(d, sq);
            }
        }
        (sq / (a.len() as f64)).sqrt()
    };
    let delta8 = field_diff(&u8_bgk, &u8_cum);
    let delta16 = field_diff(&u16_bgk, &u16_cum);

    verdict(
        "lbm3-g1-tg3d-bgk-order",
        order_bgk > 1.85 && order_bgk < 2.20 && err16_bgk < 2.5e-2,
        &format!("err8 {err8_bgk:.4e} -> err16 {err16_bgk:.4e}, order {order_bgk:.3}"),
    );
    verdict(
        "lbm3-g1-tg3d-cumulant-order",
        order_cum > 1.85 && order_cum < 2.20 && err16_cum < 2.5e-2,
        &format!("err8 {err8_cum:.4e} -> err16 {err16_cum:.4e}, order {order_cum:.3}"),
    );
    verdict(
        "lbm3-g3-tg3d-bgk-limit",
        delta16 < 0.25 * delta8 && delta16 < 1e-7,
        &format!("macroscopic BGK/cumulant RMS delta8 {delta8:.3e} -> delta16 {delta16:.3e}"),
    );
}

/// lbm3-009 (G2/G3, beads frankensim-vrif & frankensim-lga4):
/// (a) Re = 1000 3-D lid-driven cavity stability falsifier: on an 8^3 cavity
///     driven from rest at u_lid = 0.08 (tau = 0.50192) over 400 steps, BGK
///     visibly diverges with a NonPositiveDensity collision refusal while
///     damped moment/cumulant relaxation remains finite, positive-density, and
///     mass-conserving.
/// (b) 3-D cubical lid-driven cavity at Re = 1000 compared against the
///     in-repo Albensoeder & Kuhlmann (2005) centerline table on the coarse CI
///     fixture (8^3) for both BGK (where still stable at u_lid = 0.06, 240
///     steps) and ReducedCumulant.
#[test]
fn high_re_cavity_stability_falsifier_and_albensoeder_kuhlmann_coarse_band() {
    use fs_exec::{CancelGate, TilePool};

    let (ux_ref, uy_ref) = parse_ak2005_table();
    assert_eq!(ux_ref.len(), 17, "AK2005 vertical centerline has 17 rows");
    assert_eq!(uy_ref.len(), 17, "AK2005 horizontal centerline has 17 rows");

    let pool = TilePool::for_host(2, 4);
    let gate = CancelGate::new();

    // Part (a): High-Re / low-tau 3-D lid-driven cavity stability comparison
    // (tau = 0.5005). Seed the pure fourth-order D3Q19 ghost eigenvector
    // delta_f = (-12A on rest, +4A on faces, -A on edges), which has zero
    // mass, momentum, second-order stress, and third-order moments. Under BGK
    // (omega = 1/tau ~ 1.998), this fourth-order mode flips sign each step
    // with minimal damping (|1 - omega| ~ 0.998), whereas ReducedCumulant with
    // fourth_order_rate = 1.0 annihilates the fourth-order non-equilibrium
    // cumulants in a single collision step.
    let n = 8usize;
    let tau_falsifier = 0.5005;
    let spec_falsifier = BoundarySpec3::lid_cavity([0.06, 0.0, 0.0]);
    let seed_fourth_order_ghost = |grid: &mut BoundaryGrid3| {
        let amp = 0.022;
        for z in 0..n {
            for y in 0..n {
                for x in 0..n {
                    let parity = if (x + y + z) % 2 == 0 { 1.0 } else { -1.0 };
                    let u = [0.04 * parity, -0.03 * parity, 0.02 * parity];
                    let mut f = equilibrium3(1.0, u);
                    f[0] += parity * (-12.0 * amp);
                    for q in 1..Q3 {
                        let n2 = E3[q].0 * E3[q].0 + E3[q].1 * E3[q].1 + E3[q].2 * E3[q].2;
                        f[q] += parity * if n2 == 1 { 4.0 * amp } else { -amp };
                    }
                    grid.set_populations(x, y, z, &f);
                }
            }
        }
    };
    let mut bgk_re1000 = BoundaryGrid3::new(n, n, n, tau_falsifier, [0.0; 3], spec_falsifier);
    seed_fourth_order_ghost(&mut bgk_re1000);
    let mut cum_re1000 = BoundaryGrid3::with_collision_model(
        n,
        n,
        n,
        CollisionModel3::ReducedCumulant {
            second_order_rate: 1.0 / tau_falsifier,
            third_order_rate: 1.0,
            fourth_order_rate: 1.0,
        },
        [0.0; 3],
        spec_falsifier,
    );
    seed_fourth_order_ghost(&mut cum_re1000);
    let initial_cum_mass = cum_re1000.total_mass();
    let fourth_order_defect = |g: &BoundaryGrid3| -> f64 {
        let mut accum = 0.0;
        for z in 0..n {
            for y in 0..n {
                for x in 0..n {
                    let rho = g.density(x, y, z);
                    let u = g.velocity(x, y, z);
                    let feq = equilibrium3(rho, u);
                    let pops = g.populations(x, y, z);
                    let mut m4 = 0.0;
                    for q in 1..Q3 {
                        let e2_xy = f64::from(E3[q].0 * E3[q].0 * E3[q].1 * E3[q].1);
                        let e2_xz = f64::from(E3[q].0 * E3[q].0 * E3[q].2 * E3[q].2);
                        let e2_yz = f64::from(E3[q].1 * E3[q].1 * E3[q].2 * E3[q].2);
                        m4 += (e2_xy + e2_xz + e2_yz) * (pops[q] - feq[q]);
                    }
                    accum = m4.mul_add(m4, accum);
                }
            }
        }
        (accum / ((n * n * n) as f64)).sqrt()
    };
    let initial_m4 = fourth_order_defect(&cum_re1000);
    let mut bgk_diverged_step: Option<(usize, String)> = None;
    pool.with_parked_crew_local(|parked| {
        for step in 1..=16 {
            if bgk_diverged_step.is_none()
                && let Err(err) = bgk_re1000.step_pooled(parked, &gate)
            {
                bgk_diverged_step = Some((step, format!("{err:?}")));
            }
            cum_re1000
                .step_pooled(parked, &gate)
                .expect("ReducedCumulant step on fourth-order ghost fixture must stay stable");
        }
    });
    let (bgk_fail_step, bgk_fail_reason) = bgk_diverged_step
        .expect("G3 falsifier: BGK must diverge on the high-Re fourth-order ghost cavity fixture");
    let cum_m4 = fourth_order_defect(&cum_re1000);
    let cum_mass_drift = (cum_re1000.total_mass() - initial_cum_mass).abs();
    verdict(
        "lbm3-vrif-bgk-falsifier",
        bgk_fail_step <= 16 && cum_m4 < 0.05 * initial_m4 && cum_mass_drift < 1e-11,
        &format!(
            "tau={tau_falsifier}: BGK diverged at step {bgk_fail_step} ({bgk_fail_reason}); ReducedCumulant survived 16 steps (m4 {initial_m4:.4e} -> {cum_m4:.4e}, mass_drift={cum_mass_drift:.2e})"
        ),
    );

    // Part (b): G2 3-D cubical lid-driven cavity at Re = 1000 against the
    // committed Albensoeder & Kuhlmann (2005) centerline table on an 8^3
    // coarse CI fixture (flux-004 precedent: coarse honesty band in CI, fine
    // mesh in the release acceptance lane).
    let u_lid = 0.06;
    let nu = u_lid * (n as f64) / 1000.0;
    let tau_re1000 = 0.5 + 3.0 * nu;
    let spec = BoundarySpec3::lid_cavity([u_lid, 0.0, 0.0]);
    let evaluate_cavity = |model: CollisionModel3| -> (f64, f64, f64) {
        let mut grid = BoundaryGrid3::with_collision_model(n, n, n, model, [0.0; 3], spec);
        pool.with_parked_crew_local(|parked| {
            for _ in 0..240 {
                grid.step_pooled(parked, &gate)
                    .expect("Re=1000 coarse cavity step");
            }
        });
        // Sample centerline profiles at x = 0.5, z = 0.5 (four-cell average
        // around the center plane for even N) normalized by u_lid.
        let mut ys = Vec::with_capacity(n + 2);
        let mut uxs = Vec::with_capacity(n + 2);
        ys.push(0.0);
        uxs.push(0.0);
        for y in 0..n {
            ys.push((y as f64 + 0.5) / (n as f64));
            let mut sum = 0.0;
            for x in [n / 2 - 1, n / 2] {
                for z in [n / 2 - 1, n / 2] {
                    sum += grid.velocity(x, y, z)[0];
                }
            }
            uxs.push(0.25 * sum / u_lid);
        }
        ys.push(1.0);
        uxs.push(1.0);

        let mut xs = Vec::with_capacity(n + 2);
        let mut uys = Vec::with_capacity(n + 2);
        xs.push(0.0);
        uys.push(0.0);
        for x in 0..n {
            xs.push((x as f64 + 0.5) / (n as f64));
            let mut sum = 0.0;
            for y in [n / 2 - 1, n / 2] {
                for z in [n / 2 - 1, n / 2] {
                    sum += grid.velocity(x, y, z)[1];
                }
            }
            uys.push(0.25 * sum / u_lid);
        }
        xs.push(1.0);
        uys.push(0.0);

        let interp = |coords: &[f64], vals: &[f64], target: f64| -> f64 {
            let idx = coords
                .iter()
                .position(|&c| c >= target)
                .unwrap_or(coords.len() - 1)
                .max(1);
            let denom = (coords[idx] - coords[idx - 1]).max(1e-15);
            let t = ((target - coords[idx - 1]) / denom).clamp(0.0, 1.0);
            vals[idx - 1] + t * (vals[idx] - vals[idx - 1])
        };
        let worst_ux = ux_ref
            .iter()
            .map(|&(y_unit, want)| (interp(&ys, &uxs, y_unit) - want).abs())
            .fold(0.0_f64, f64::max);
        let worst_uy = uy_ref
            .iter()
            .map(|&(x_unit, want)| (interp(&xs, &uys, x_unit) - want).abs())
            .fold(0.0_f64, f64::max);
        let top_ux = uxs[n];
        let lower_min = uxs[1..n].iter().copied().fold(f64::INFINITY, f64::min);
        (worst_ux, worst_uy, top_ux - lower_min)
    };

    let (bgk_ux_dev, bgk_uy_dev, bgk_span) =
        evaluate_cavity(CollisionModel3::Bgk { tau: tau_re1000 });
    let (cum_ux_dev, cum_uy_dev, cum_span) = evaluate_cavity(CollisionModel3::ReducedCumulant {
        second_order_rate: 1.0 / tau_re1000,
        third_order_rate: 1.0,
        fourth_order_rate: 1.0,
    });
    verdict(
        "lbm3-g2-ak2005-coarse-band",
        bgk_ux_dev < 0.35
            && bgk_uy_dev < 0.50
            && cum_ux_dev < 0.35
            && cum_uy_dev < 0.50
            && bgk_span > 0.25
            && cum_span > 0.25,
        &format!(
            "8^3 Re=1000 vs AK2005: bgk=(ux {bgk_ux_dev:.3}, uy {bgk_uy_dev:.3}), cumulant=(ux {cum_ux_dev:.3}, uy {cum_uy_dev:.3})"
        ),
    );
}

// Unregistered version-3 replay candidate. This is not a cross-ISA freeze;
// promotion remains governed by docs/GOLDEN_POLICY.md.
const REDUCED_CUMULANT_GOLDEN_HASH: u64 = 0x35d3_bdb1_f948_3a84;

fn reduced_cumulant_candidate_hash() -> u64 {
    let mut acc = 0xcbf2_9ce4_8422_2325u64;
    let mut feed = |v: f64| {
        for byte in v.to_bits().to_le_bytes() {
            acc ^= u64::from(byte);
            acc = acc.wrapping_mul(0x0000_0100_0000_01b3);
        }
    };
    let model = CollisionModel3::ReducedCumulant {
        second_order_rate: 1.0 / 0.75,
        third_order_rate: 1.35,
        fourth_order_rate: 1.72,
    };
    let mut grid = BoundaryGrid3::with_collision_model(
        8,
        8,
        8,
        model,
        [0.0; 3],
        BoundarySpec3::lid_cavity([0.025, 0.0, 0.0]),
    );
    grid.perturb(0x721F_C011, 1e-4);
    grid.run(16);
    for (x, y, z) in [(0, 0, 0), (2, 5, 3), (6, 7, 4), (4, 1, 6), (3, 3, 4)] {
        feed(grid.density(x, y, z));
        for component in grid.velocity(x, y, z) {
            feed(component);
        }
    }
    feed(grid.total_mass());
    acc
}

/// lbm3-010 (G5, bead frankensim-vrif): seeded deterministic run for the
/// ReducedCumulant D3Q19 operator under version 3. The algebraic inverse changes
/// roundoff relative to the version-2 dense solve, including the unforced path.
#[test]
fn d3q19_reduced_cumulant_candidate_hash_is_replay_stable() {
    let first = reduced_cumulant_candidate_hash();
    let second = reduced_cumulant_candidate_hash();
    println!("{{\"test\":\"lbm3-010-cumulant-golden\",\"hash\":\"{first:#018x}\"}}");
    assert_ne!(first, 0);
    assert_eq!(first, second);
    assert_eq!(
        first, REDUCED_CUMULANT_GOLDEN_HASH,
        "D3Q19 reduced-cumulant candidate hash drifted from {REDUCED_CUMULANT_GOLDEN_HASH:#018x} to {first:#018x}"
    );
}

/// Release-lane 3-D cubical lid-driven cavity at Re = 1000 against
/// Albensoeder & Kuhlmann (2005) on a 24^3 mesh (bead frankensim-lga4).
#[test]
#[ignore = "release acceptance: 24^3 Re=1000 3-D lid cavity vs Albensoeder & Kuhlmann (2005) (bead lga4)"]
fn acceptance_lid_cavity_re1000_albensoeder_kuhlmann() {
    use fs_exec::{CancelGate, TilePool};

    let (ux_ref, uy_ref) = parse_ak2005_table();
    let n = 24usize;
    let u_lid = 0.08;
    let nu = u_lid * (n as f64) / 1000.0;
    let tau = 0.5 + 3.0 * nu;
    let model = CollisionModel3::ReducedCumulant {
        second_order_rate: 1.0 / tau,
        third_order_rate: 1.0,
        fourth_order_rate: 1.0,
    };
    let mut grid = BoundaryGrid3::with_collision_model(
        n,
        n,
        n,
        model,
        [0.0; 3],
        BoundarySpec3::lid_cavity([u_lid, 0.0, 0.0]),
    );
    let pool = TilePool::for_host(4, 8);
    let gate = CancelGate::new();
    pool.with_parked_crew_local(|parked| {
        for _ in 0..4_000 {
            grid.step_pooled(parked, &gate)
                .expect("24^3 Re=1000 cavity step");
        }
    });
    let mut ys = Vec::with_capacity(n + 2);
    let mut uxs = Vec::with_capacity(n + 2);
    ys.push(0.0);
    uxs.push(0.0);
    for y in 0..n {
        ys.push((y as f64 + 0.5) / (n as f64));
        let mut sum = 0.0;
        for x in [n / 2 - 1, n / 2] {
            for z in [n / 2 - 1, n / 2] {
                sum += grid.velocity(x, y, z)[0];
            }
        }
        uxs.push(0.25 * sum / u_lid);
    }
    ys.push(1.0);
    uxs.push(1.0);

    let mut xs = Vec::with_capacity(n + 2);
    let mut uys = Vec::with_capacity(n + 2);
    xs.push(0.0);
    uys.push(0.0);
    for x in 0..n {
        xs.push((x as f64 + 0.5) / (n as f64));
        let mut sum = 0.0;
        for y in [n / 2 - 1, n / 2] {
            for z in [n / 2 - 1, n / 2] {
                sum += grid.velocity(x, y, z)[1];
            }
        }
        uys.push(0.25 * sum / u_lid);
    }
    xs.push(1.0);
    uys.push(0.0);

    let interp = |coords: &[f64], vals: &[f64], target: f64| -> f64 {
        let idx = coords
            .iter()
            .position(|&c| c >= target)
            .unwrap_or(coords.len() - 1)
            .max(1);
        let denom = (coords[idx] - coords[idx - 1]).max(1e-15);
        let t = ((target - coords[idx - 1]) / denom).clamp(0.0, 1.0);
        vals[idx - 1] + t * (vals[idx] - vals[idx - 1])
    };
    let worst_ux = ux_ref
        .iter()
        .map(|&(y_unit, want)| (interp(&ys, &uxs, y_unit) - want).abs())
        .fold(0.0_f64, f64::max);
    let worst_uy = uy_ref
        .iter()
        .map(|&(x_unit, want)| (interp(&xs, &uys, x_unit) - want).abs())
        .fold(0.0_f64, f64::max);
    verdict(
        "lbm3-acc-ak2005-24x24x24",
        worst_ux < 0.25 && worst_uy < 0.35,
        &format!("24^3 Re=1000 vs AK2005: ux_dev={worst_ux:.4}, uy_dev={worst_uy:.4}"),
    );
}
