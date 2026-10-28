//! IRP .130 G0/G1/G2/G4/G5. Native release execution is central.
//! The axial Gaussian fixture tests rotational error of a translationally
//! invariant 3-D jet cross-section, not nozzle/edge-tone qualification.

use fs_exec::{CancelGate, Cx, TilePool};
use fs_lbm::d3q27::*;
use fs_lbm::lattice3::{
    AcousticScalingPlan, D3q19, Grid3, GridBudget3, GridError3, acoustic_scaling,
};
use fs_material::gas::{GasSpec, GasState};
use fs_math::det;
use std::time::Instant;

const TAU: f64 = core::f64::consts::TAU;
const SEED: u64 = 0xD327_0130;
// Independent exact-product oracle: input theta=1/4, stress rates3/2
// produce theta=3/8. All higher cumulants vanish; tensor inverse is dyadic.
// FNV1a over post populations in E27 order, calculated independently of Rust.
const D3Q27_GOLDEN_HASH: u64 = 0xb256_c77f_0c5c_1bf7;

fn rates(shear: f64, higher: f64) -> Rates27 {
    Rates27 {
        omega: [
            shear, 1.0, higher, higher, higher, higher, higher, higher, higher, higher,
        ],
    }
}
fn config(model: CollisionModel27) -> CollisionConfig27 {
    CollisionConfig27 {
        model,
        maximum_mach: 0.2,
    }
}
fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: SEED,
                kernel_id: 130,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        f(&cx)
    })
}
fn gas() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap()
}
fn gas_identity(g: &GasState) -> String {
    let bytes: Vec<_> = [
        g.temperature,
        g.pressure,
        g.density,
        g.sound_speed,
        g.dynamic_viscosity,
    ]
    .into_iter()
    .flat_map(|v| v.to_bits().to_le_bytes())
    .collect();
    fs_blake3::hash_domain("org.frankensim.lbm.gas-fixture.v1", &bytes).to_hex()
}
fn emit(lane: &str, stage: &str, lattice: &str, verdict: &str, detail: &str, start: Instant) {
    let g = gas();
    println!(
        "{{\"suite\":\"fs-lbm/d3q27\",\"lane\":\"{lane}\",\"stage\":\"{stage}\",\"lattice\":\"{lattice}\",\"verdict\":\"{verdict}\",\"seed\":{SEED},\"inputs\":{{\"gas_state\":\"dry-air-USSA1976\",\"T_K\":{},\"p_Pa\":{},\"nu_m2_s\":{},\"c_m_s\":{},\"blake3\":\"{}\"}},\"detail\":{detail},\"elapsed_s\":{}}}",
        g.temperature,
        g.pressure,
        g.dynamic_viscosity / g.density,
        g.sound_speed,
        gas_identity(&g),
        start.elapsed().as_secs_f64()
    );
}
fn random(seed: &mut u64) -> f64 {
    *seed = seed
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    (*seed >> 11) as f64 / ((1u64 << 53) as f64)
}
fn conserved(f: &[f64; Q27]) -> [f64; 4] {
    core::array::from_fn(|a| {
        f.iter()
            .enumerate()
            .map(|(q, v)| {
                if a == 0 {
                    *v
                } else {
                    *v * f64::from(E27[q][a - 1])
                }
            })
            .sum()
    })
}

#[test]
fn d3q27_lattice_invariants_hold_exactly() {
    assert_eq!(D3Q27_BIT_SEMANTICS_VERSION, 3);
    assert_eq!(W216.iter().sum::<i64>(), 216);
    for q in 0..Q27 {
        assert_eq!(W27[q].to_bits(), (W216[q] as f64 / 216.0).to_bits());
        assert_eq!(OPP27[OPP27[q]], q);
        assert_eq!(E27[OPP27[q]], E27[q].map(|v| -v));
        assert_eq!(W216[OPP27[q]], W216[q]);
        assert!(E27[..q].iter().all(|e| *e != E27[q]));
    }
    for a in 0..3 {
        for b in 0..3 {
            let second: i64 = (0..Q27)
                .map(|q| W216[q] * i64::from(E27[q][a] * E27[q][b]))
                .sum();
            assert_eq!(second, if a == b { 72 } else { 0 });
            for c in 0..3 {
                for d in 0..3 {
                    let fourth: i64 = (0..Q27)
                        .map(|q| W216[q] * i64::from(E27[q][a] * E27[q][b] * E27[q][c] * E27[q][d]))
                        .sum();
                    let delta = |i, j| i64::from(i == j);
                    assert_eq!(
                        fourth,
                        24 * (delta(a, b) * delta(c, d)
                            + delta(a, c) * delta(b, d)
                            + delta(a, d) * delta(b, c))
                    );
                }
            }
        }
    }
    // The mixed sixth moment absent from D3Q19 is present exactly.
    assert_eq!(
        (0..Q27)
            .map(|q| W216[q] * i64::from(E27[q].iter().map(|e| e * e).product::<i32>()))
            .sum::<i64>(),
        8
    );
}

#[test]
fn d3q27_equilibrium_recovers_moments() {
    for form in [Equilibrium27::Quadratic, Equilibrium27::Factorized] {
        for rho in [0.7, 1.0, 1.3] {
            let u = [0.03, -0.04, 0.02];
            let f = equilibrium27(rho, u, form);
            let actual = conserved(&f);
            assert!((actual[0] / rho - 1.0).abs() < 1e-14);
            for a in 0..3 {
                assert!((actual[a + 1] - rho * u[a]).abs() < 1e-14 * rho);
                for b in 0..3 {
                    let stress: f64 = (0..Q27)
                        .map(|q| f[q] * f64::from(E27[q][a] * E27[q][b]))
                        .sum();
                    let exact = rho * (u[a] * u[b] + if a == b { CS2_27 } else { 0.0 });
                    assert!((stress - exact).abs() < 1e-14 * rho);
                }
            }
        }
    }
}

#[test]
fn cumulant_collision_conserves_mass_momentum() {
    let start = Instant::now();
    let mut seed = SEED;
    for _ in 0..10_000 {
        let f = core::array::from_fn(|q| W27[q] * (0.99 + 0.02 * random(&mut seed)));
        let before = conserved(&f);
        for model in [
            CollisionModel27::CentralMoment {
                rates: rates(1.8, 0.7),
            },
            CollisionModel27::Cumulant {
                rates: rates(1.8, 0.7),
            },
        ] {
            for force in [[0.0; 3], [1e-6, -2e-6, 0.5e-6]] {
                let out = collide_cell27(f, config(model), force).unwrap();
                let after = conserved(&out);
                for a in 0..4 {
                    let expected = before[a] + if a > 0 { force[a - 1] } else { 0.0 };
                    assert!(
                        (after[a] - expected).abs() <= 1e-14 * before[0],
                        "a{a}: {:?} -> {:?}",
                        before,
                        after
                    );
                }
            }
        }
    }
    emit(
        "unit",
        "conservation",
        "D3Q27",
        "pass",
        "{\"states\":10000,\"tolerance\":1e-14}",
        start,
    );
}

// Independent product probabilities and their analytic velocity derivative.
// These do not call the production equilibrium, moment transform or recurrence.
fn factorized_force_oracle(rho: f64, u: [f64; 3], force: [f64; 3]) -> ([f64; Q27], [f64; Q27]) {
    let factors: [[(f64, f64); 3]; Q27] = core::array::from_fn(|q| {
        core::array::from_fn(|a| {
            let c = f64::from(E27[q][a]);
            if c == 0.0 {
                (2.0 / 3.0 - u[a] * u[a], -2.0 * u[a])
            } else {
                (0.5 * (1.0 / 3.0 + u[a] * u[a] + c * u[a]), u[a] + 0.5 * c)
            }
        })
    });
    let eq = core::array::from_fn(|q| rho * factors[q].iter().map(|v| v.0).product::<f64>());
    let source = core::array::from_fn(|q| {
        (0..3)
            .map(|a| {
                force[a] * factors[q][a].1 * factors[q][(a + 1) % 3].0 * factors[q][(a + 2) % 3].0
            })
            .sum()
    });
    (eq, source)
}

fn moment_about(f: &[f64; Q27], u: [f64; 3], powers: [usize; 3]) -> f64 {
    (0..Q27)
        .map(|q| {
            f[q] * (0..3)
                .map(|a| {
                    let c = f64::from(E27[q][a]) - u[a];
                    match powers[a] {
                        0 => 1.0,
                        1 => c,
                        2 => c * c,
                        _ => unreachable!(),
                    }
                })
                .product::<f64>()
        })
        .sum()
}

#[test]
fn factorized_bgk_force_matches_maxwellian_moments_in_moving_frames() {
    // G0/G3: f=feq-S/2 represents a physical Maxwellian in the trapezoidal
    // source frame. Its postcollision state is feq+S/2 at every relaxation
    // time. A quadratic source passes mass/impulse/stress but fails the higher
    // moments when the selected equilibrium is factorized and u is nonzero.
    for rho in [0.8, 1.4] {
        for u in [[0.0; 3], [0.065, -0.04, 0.025], [-0.035, 0.06, 0.045]] {
            for tau in [0.5005, 0.8, 1.7] {
                for direction in [
                    [8e-4, 0.0, 0.0],
                    [0.0, -5e-4, 0.0],
                    [0.0, 0.0, 3e-4],
                    [7e-4, -4e-4, 2e-4],
                ] {
                    for sign in [-1.0, 1.0] {
                        let force = direction.map(|v| sign * v);
                        let (eq, source) = factorized_force_oracle(rho, u, force);
                        let before = core::array::from_fn(|q| eq[q] - 0.5 * source[q]);
                        let post = collide_cell27(
                            before,
                            config(CollisionModel27::Bgk {
                                tau,
                                equilibrium: Equilibrium27::Factorized,
                            }),
                            force,
                        )
                        .unwrap();
                        for q in 0..Q27 {
                            let expected = eq[q] + 0.5 * source[q];
                            assert!(
                                (post[q] - expected).abs() < 3e-14 * rho,
                                "population {q}: got {}, expected {expected}; u={u:?}, force={force:?}, tau={tau}",
                                post[q]
                            );
                        }
                        // Independent Gaussian moment identities exercise the
                        // complete tensor, not only the conserved quantities.
                        for z in 0..=2 {
                            for y in 0..=2 {
                                for x in 0..=2 {
                                    let powers = [x, y, z];
                                    let odd = powers.iter().filter(|&&p| p == 1).count();
                                    let gaussian = powers
                                        .iter()
                                        .filter(|&&p| p == 2)
                                        .fold(1.0, |v, _| v * CS2_27);
                                    let expected = match odd {
                                        0 => rho * gaussian,
                                        1 => {
                                            0.5 * force
                                                [powers.iter().position(|&p| p == 1).unwrap()]
                                                * gaussian
                                        }
                                        _ => 0.0,
                                    };
                                    assert!(
                                        (moment_about(&post, u, powers) - expected).abs()
                                            < 5e-14 * rho,
                                        "central moment {powers:?}; u={u:?}, force={force:?}, tau={tau}"
                                    );
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

#[test]
fn central_moment_force_matches_factorized_equilibrium_derivative() {
    // G0/G3: the transformed population is f = feq - S/2. For an equilibrium
    // physical population, collision must give feq + S/2 for ANY grouped rates.
    // S=(F/rho).d_u feq has the Gaussian central moments in Fei et al. (2018),
    // equations5,20--21: cs²F at third order and cs⁴F at fifth, with zero stress.
    for rho in [0.8, 1.4] {
        for u in [[0.0; 3], [0.065, -0.04, 0.025]] {
            for omega in [
                [0.65, 1.7, 0.4, 1.5, 0.7, 0.9, 1.3, 1.1, 0.25, 1.6],
                [1.85, 0.6, 1.8, 0.3, 1.2, 1.4, 0.8, 1.7, 1.75, 0.5],
            ] {
                for direction in [
                    [8e-4, 0.0, 0.0],
                    [0.0, -5e-4, 0.0],
                    [0.0, 0.0, 3e-4],
                    [7e-4, -4e-4, 2e-4],
                ] {
                    for sign in [-1.0, 1.0] {
                        let force = direction.map(|v| sign * v);
                        let (eq, source) = factorized_force_oracle(rho, u, force);
                        let before = core::array::from_fn(|q| eq[q] - 0.5 * source[q]);
                        let post = collide_cell27(
                            before,
                            config(CollisionModel27::CentralMoment {
                                rates: Rates27 { omega },
                            }),
                            force,
                        )
                        .unwrap();
                        for q in 0..Q27 {
                            let expected = eq[q] + 0.5 * source[q];
                            assert!(
                                (post[q] - expected).abs() < 3e-14 * rho,
                                "population {q}: got {}, expected {expected}; force={force:?}, rates={omega:?}",
                                post[q]
                            );
                        }
                        // All independent moments include mass, impulse, each
                        // shear/normal stress, and all third/fifth force terms.
                        for z in 0..=2 {
                            for y in 0..=2 {
                                for x in 0..=2 {
                                    let powers = [x, y, z];
                                    let odd = powers.iter().filter(|&&p| p == 1).count();
                                    let gaussian = powers
                                        .iter()
                                        .filter(|&&p| p == 2)
                                        .fold(1.0, |v, _| v * CS2_27);
                                    let expected = match odd {
                                        0 => rho * gaussian,
                                        1 => {
                                            0.5 * force
                                                [powers.iter().position(|&p| p == 1).unwrap()]
                                                * gaussian
                                        }
                                        _ => 0.0,
                                    };
                                    let observed = moment_about(&post, u, powers);
                                    assert!(
                                        (observed - expected).abs() < 5e-14 * rho,
                                        "central moment {powers:?}: got {observed}, expected {expected}; force={force:?}"
                                    );
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

#[test]
fn central_moment_force_preserves_independent_shear_and_bulk_relaxation() {
    // A Hermite stress perturbation has independently prescribed normal/shear
    // stresses and no mass, momentum or odd moments. The force may not change
    // their viscous rates or create an additional third/fifth nonequilibrium.
    let rho = 1.2;
    let stress = [
        [5e-4, 2e-4, -1e-4],
        [2e-4, -1e-4, 1.5e-4],
        [-1e-4, 1.5e-4, 3e-4],
    ];
    let trace = stress[0][0] + stress[1][1] + stress[2][2];
    for omega in [
        [1.7, 0.9, 0.6, 0.8, 0.7, 1.2, 1.1, 0.65, 0.75, 0.85],
        [0.7, 1.6, 1.4, 0.5, 1.1, 0.8, 1.2, 1.5, 1.8, 0.4],
    ] {
        for sign in [-1.0, 1.0] {
            let force = [1e-3, -5e-4, 2.5e-4].map(|v| sign * v);
            let (eq, source) = factorized_force_oracle(rho, [0.0; 3], force);
            let before = core::array::from_fn(|q| {
                let e = E27[q].map(f64::from);
                let mut perturbation = 0.0;
                for a in 0..3 {
                    for b in 0..3 {
                        perturbation +=
                            (e[a] * e[b] - if a == b { CS2_27 } else { 0.0 }) * stress[a][b];
                    }
                }
                eq[q] + W27[q] * perturbation / (2.0 * CS2_27 * CS2_27) - 0.5 * source[q]
            });
            let post = collide_cell27(
                before,
                config(CollisionModel27::CentralMoment {
                    rates: Rates27 { omega },
                }),
                force,
            )
            .unwrap();
            for a in 0..3 {
                for b in a..3 {
                    let mut powers = [0; 3];
                    powers[a] += 1;
                    powers[b] += 1;
                    let isotropic = if a == b { trace / 3.0 } else { 0.0 };
                    let expected = if a == b { rho * CS2_27 } else { 0.0 }
                        + (1.0 - omega[0]) * (stress[a][b] - isotropic)
                        + (1.0 - omega[1]) * isotropic;
                    assert!((moment_about(&post, [0.0; 3], powers) - expected).abs() < 5e-14 * rho);
                }
            }
            for (axis, powers) in [
                (0, [1, 2, 0]),
                (0, [1, 0, 2]),
                (1, [2, 1, 0]),
                (1, [0, 1, 2]),
                (2, [2, 0, 1]),
                (2, [0, 2, 1]),
                (0, [1, 2, 2]),
                (1, [2, 1, 2]),
                (2, [2, 2, 1]),
            ] {
                let gaussian = powers
                    .iter()
                    .filter(|&&p| p == 2)
                    .fold(1.0, |v, _| v * CS2_27);
                let expected = 0.5 * force[axis] * gaussian;
                assert!((moment_about(&post, [0.0; 3], powers) - expected).abs() < 5e-14 * rho);
            }
        }
    }
}

fn central(f: &[f64; Q27], powers: [usize; 3]) -> f64 {
    let c = conserved(f);
    let u = [c[1] / c[0], c[2] / c[0], c[3] / c[0]];
    (0..Q27)
        .map(|q| {
            f[q] / c[0]
                * (0..3)
                    .map(|a| det::powi(f64::from(E27[q][a]) - u[a], powers[a] as i32))
                    .product::<f64>()
        })
        .sum()
}
/// Independent set-partition definition (Bell(6)=203), sharing neither the
/// production recurrence nor its separable central transform.
fn cumulant_oracle(f: &[f64; Q27], powers: [usize; 3]) -> f64 {
    let axes: Vec<_> = (0..3)
        .flat_map(|a| std::iter::repeat_n(a, powers[a]))
        .collect();
    fn visit(
        f: &[f64; Q27],
        axes: &[usize],
        i: usize,
        blocks: &mut Vec<[usize; 3]>,
        out: &mut f64,
    ) {
        if i == axes.len() {
            let k = blocks.len();
            let factorial = (1..k).product::<usize>() as f64;
            let sign = if k % 2 == 1 { 1.0 } else { -1.0 };
            *out += sign * factorial * blocks.iter().map(|p| central(f, *p)).product::<f64>();
            return;
        }
        for b in 0..blocks.len() {
            blocks[b][axes[i]] += 1;
            visit(f, axes, i + 1, blocks, out);
            blocks[b][axes[i]] -= 1;
        }
        let mut block = [0; 3];
        block[axes[i]] = 1;
        blocks.push(block);
        visit(f, axes, i + 1, blocks, out);
        blocks.pop();
    }
    let mut result = 0.0;
    visit(f, &axes, 0, &mut Vec::new(), &mut result);
    result
}

#[test]
fn full_cumulants_through_sixth_order_match_independent_partitions() {
    let mut seed = SEED;
    let r = Rates27 {
        omega: [1.7, 0.9, 0.6, 0.8, 0.7, 1.2, 1.1, 0.65, 0.75, 0.85],
    };
    for _ in 0..32 {
        let f = core::array::from_fn(|q| W27[q] * (0.8 + 0.4 * random(&mut seed)));
        let post =
            collide_cell27(f, config(CollisionModel27::Cumulant { rates: r }), [0.0; 3]).unwrap();
        for (p, rate) in [
            ([1, 1, 1], r.omega[4]),
            ([2, 1, 1], r.omega[7]),
            ([1, 2, 1], r.omega[7]),
            ([1, 1, 2], r.omega[7]),
            ([2, 2, 1], r.omega[8]),
            ([2, 1, 2], r.omega[8]),
            ([1, 2, 2], r.omega[8]),
            ([2, 2, 2], r.omega[9]),
        ] {
            let before = cumulant_oracle(&f, p);
            let actual = cumulant_oracle(&post, p);
            assert!(
                (actual - (1.0 - rate) * before).abs() < 2e-14,
                "{p:?}: {actual} vs {}",
                (1.0 - rate) * before
            );
        }
        let cm = collide_cell27(
            f,
            config(CollisionModel27::CentralMoment { rates: r }),
            [0.0; 3],
        )
        .unwrap();
        assert!(
            cm.iter().zip(post).any(|(a, b)| (*a - b).abs() > 1e-9),
            "central moments must not alias full cumulants"
        );
    }
}

#[test]
fn geier_diagonal_cubic_correction_matches_publication() {
    let r = Rates27 {
        omega: [1.7, 0.9, 0.6, 0.8, 0.7, 1.2, 1.1, 0.65, 0.75, 0.85],
    };
    let mut seed = SEED;
    let equilibrium = equilibrium27(1.1, [0.07, -0.05, 0.02], Equilibrium27::Factorized);
    for _ in 0..32 {
        let f = core::array::from_fn(|q| equilibrium[q] * (0.98 + 0.04 * random(&mut seed)));
        let state = conserved(&f);
        let rho = state[0];
        let u = [state[1] / rho, state[2] / rho, state[3] / rho];
        let c = [
            central(&f, [2, 0, 0]) * rho,
            central(&f, [0, 2, 0]) * rho,
            central(&f, [0, 0, 2]) * rho,
        ];
        let trace = c.iter().sum::<f64>();
        let dx = -r.omega[0] / (2.0 * rho) * (2.0 * c[0] - c[1] - c[2])
            - r.omega[1] / (2.0 * rho) * (trace - rho);
        let dy = dx + 3.0 * r.omega[0] / (2.0 * rho) * (c[0] - c[1]);
        let dz = dx + 3.0 * r.omega[0] / (2.0 * rho) * (c[0] - c[2]);
        let correction = [u[0] * u[0] * dx, u[1] * u[1] * dy, u[2] * u[2] * dz];
        let expected_xy = (1.0 - r.omega[0]) * (c[0] - c[1])
            - 3.0 * rho * (1.0 - r.omega[0] / 2.0) * (correction[0] - correction[1]);
        let expected_xz = (1.0 - r.omega[0]) * (c[0] - c[2])
            - 3.0 * rho * (1.0 - r.omega[0] / 2.0) * (correction[0] - correction[2]);
        let expected_trace = rho * r.omega[1] + (1.0 - r.omega[1]) * trace
            - 3.0 * rho * (1.0 - r.omega[1] / 2.0) * correction.iter().sum::<f64>();
        for model in [
            CollisionModel27::CentralMoment { rates: r },
            CollisionModel27::Cumulant { rates: r },
        ] {
            let post = collide_cell27(f, config(model), [0.0; 3]).unwrap();
            let p = [
                central(&post, [2, 0, 0]) * rho,
                central(&post, [0, 2, 0]) * rho,
                central(&post, [0, 0, 2]) * rho,
            ];
            assert!((p[0] - p[1] - expected_xy).abs() < 1e-14 * rho);
            assert!((p[0] - p[2] - expected_xz).abs() < 1e-14 * rho);
            assert!((p.iter().sum::<f64>() - expected_trace).abs() < 1e-14 * rho);
        }
        assert!(
            (expected_xy - (1.0 - r.omega[0]) * (c[0] - c[1])).abs() > 1e-10,
            "uncorrected collision must be distinguishable"
        );
    }
}

#[test]
fn refusals_cap_and_cap_plus_one() {
    let rest = equilibrium27(1.0, [0.0; 3], Equilibrium27::Factorized);
    for tau in [0.5, 0.0, f64::NAN, f64::INFINITY] {
        assert!(
            collide_cell27(
                rest,
                config(CollisionModel27::Bgk {
                    tau,
                    equilibrium: Equilibrium27::Quadratic
                }),
                [0.0; 3]
            )
            .is_err()
        );
    }
    collide_cell27(
        rest,
        config(CollisionModel27::Bgk {
            tau: 0.5 + 1e-12,
            equilibrium: Equilibrium27::Quadratic,
        }),
        [0.0; 3],
    )
    .unwrap();
    for i in 0..10 {
        for bad in [0.0, 2.0, f64::NAN, f64::INFINITY] {
            let mut r = rates(1.0, 1.0);
            r.omega[i] = bad;
            assert!(
                collide_cell27(
                    rest,
                    config(CollisionModel27::Cumulant { rates: r }),
                    [0.0; 3]
                )
                .is_err()
            );
        }
    }
    let mut r = rates(1.0, 1.0);
    r.omega = [2.0 - 1e-12; 10];
    r.validate().unwrap();
    collide_cell27(
        rest,
        config(CollisionModel27::Cumulant { rates: r }),
        [0.0; 3],
    )
    .unwrap();
    let mut bad = rest;
    bad[26] = f64::NAN;
    assert!(
        collide_cell27(
            bad,
            config(CollisionModel27::Cumulant { rates: r }),
            [0.0; 3]
        )
        .is_err()
    );
    assert!(
        collide_cell27(
            rest,
            config(CollisionModel27::Cumulant { rates: r }),
            [f64::NAN, 0.0, 0.0]
        )
        .is_err()
    );
    let fast = equilibrium27(1.0, [0.2, 0.0, 0.0], Equilibrium27::Factorized);
    assert!(matches!(
        collide_cell27(
            fast,
            config(CollisionModel27::Cumulant { rates: r }),
            [0.0; 3]
        ),
        Err(CollisionError27::Mach { .. })
    ));
    let make = |dims| {
        Grid27::try_new(
            dims,
            config(CollisionModel27::Cumulant {
                rates: rates(1.0, 1.0),
            }),
            [0.0; 3],
            [true; 3],
            GridBudget3 {
                maximum_cells: 8,
                maximum_mach: 0.2,
            },
            |_| false,
            |_| rest,
        )
    };
    make([2, 2, 2]).unwrap();
    assert!(matches!(
        make([3, 3, 1]),
        Err(GridError3::Cap {
            requested: 9,
            maximum: 8
        })
    ));
    assert!(make([0, 2, 2]).is_err());
    let g = gas();
    for dx in [200e-6, 100e-6, 50e-6] {
        let p = AcousticScalingPlan::from_gas(&g, dx, 0.02 * g.sound_speed, 0.2).unwrap();
        assert!((p.dt_s * g.sound_speed / dx - CS2_27.sqrt()).abs() < 1e-14);
        assert!(
            (p.lattice.viscosity * dx * dx / p.dt_s / (g.dynamic_viscosity / g.density) - 1.0)
                .abs()
                < 1e-14
        );
    }
    assert!(acoustic_scaling(0.0, g.sound_speed, 1e-4, 1.0, 0.2).is_err());
    assert!(AcousticScalingPlan::from_gas(&g, 0.0, 1.0, 0.2).is_err());
    assert!(AcousticScalingPlan::from_gas(&g, 1e-4, 0.21 * g.sound_speed, 0.2).is_err());
}

fn seeded_grid() -> Grid27 {
    let c = config(CollisionModel27::Cumulant {
        rates: rates(1.7, 0.8),
    });
    Grid27::try_new(
        [8, 8, 4],
        c,
        [1e-7, 0.0, 0.0],
        [true, false, true],
        GridBudget3 {
            maximum_cells: 256,
            maximum_mach: 0.2,
        },
        |p| p == [3, 3, 2],
        |p| {
            let u = [
                0.005 * det::sin(TAU * p[1] as f64 / 8.0),
                0.002 * det::cos(TAU * p[0] as f64 / 8.0),
                0.001,
            ];
            equilibrium27(1.0, u, Equilibrium27::Factorized)
        },
    )
    .unwrap()
}
fn bits<L: fs_lbm::lattice3::Lattice3<Q>, const Q: usize>(grid: &Grid3<L, Q>) -> Vec<u64> {
    grid.populations()
        .iter()
        .flatten()
        .map(|v| v.to_bits())
        .collect()
}

#[test]
fn d3q27_bitwise_worker_count_invariance() {
    let mut reference = seeded_grid();
    with_cx(|cx| {
        for _ in 0..12 {
            reference.step_serial(cx).unwrap();
        }
    });
    for workers in [1, 2, 7] {
        let pool = TilePool::for_host(workers, SEED);
        let gate = CancelGate::new();
        let mut grid = seeded_grid();
        for _ in 0..12 {
            grid.step_pooled(&pool, &gate).unwrap();
        }
        assert_eq!(bits(&grid), bits(&reference), "workers={workers}");
    }
    let mut grid = seeded_grid();
    let before = bits(&grid);
    let gate = CancelGate::new();
    gate.request();
    assert!(matches!(
        grid.step_pooled(&TilePool::for_host(1, SEED), &gate),
        Err(GridError3::Cancelled)
    ));
    assert_eq!(bits(&grid), before);
    assert_eq!(grid.completed_steps(), 0);
}

#[test]
fn cancellation_between_passes_keeps_published_state_and_reissues_exactly() {
    use fs_exec::{KernelRunner, RunError, RunReport, TileKernel};
    use std::sync::atomic::{AtomicUsize, Ordering};
    struct CancelAfterCollision {
        pool: TilePool,
        calls: AtomicUsize,
    }
    impl KernelRunner for CancelAfterCollision {
        fn workers(&self) -> usize {
            self.pool.workers()
        }
        fn run_with_gate<K: TileKernel>(
            &self,
            kernel: &K,
            gate: &CancelGate,
        ) -> (Result<K::Out, RunError>, RunReport) {
            let result = self.pool.run_with_gate(kernel, gate);
            if self.calls.fetch_add(1, Ordering::AcqRel) == 0 {
                gate.request();
            }
            result
        }
    }
    let runner = CancelAfterCollision {
        pool: TilePool::for_host(2, SEED),
        calls: AtomicUsize::new(0),
    };
    let gate = CancelGate::new();
    let mut grid = seeded_grid();
    let before = bits(&grid);
    assert!(matches!(
        grid.step_pooled(&runner, &gate),
        Err(GridError3::Cancelled)
    ));
    assert_eq!(bits(&grid), before);
    assert_eq!(grid.completed_steps(), 0);
    let fresh = CancelGate::new();
    grid.step_pooled(&runner.pool, &fresh).unwrap();
    let mut reference = seeded_grid();
    with_cx(|cx| reference.step_serial(cx).unwrap());
    assert_eq!(bits(&grid), bits(&reference));
}

#[test]
fn d3q27_golden_hash() {
    let theta = 0.25;
    let input = core::array::from_fn(|q| {
        (0..3)
            .map(|a| {
                if E27[q][a] == 0 {
                    1.0 - theta
                } else {
                    0.5 * theta
                }
            })
            .product()
    });
    let post = collide_cell27(
        input,
        config(CollisionModel27::Cumulant {
            rates: Rates27 {
                omega: [1.5, 1.5, 0.6, 0.7, 0.8, 0.9, 1.1, 1.2, 1.3, 1.4],
            },
        }),
        [0.0; 3],
    )
    .unwrap();
    let mut hash = 0xcbf29ce484222325u64;
    for (q, v) in post.into_iter().enumerate() {
        let exact: f64 = (0..3)
            .map(|a| {
                if E27[q][a] == 0 {
                    5.0 / 8.0
                } else {
                    3.0 / 16.0
                }
            })
            .product();
        assert_eq!(
            v.to_bits(),
            exact.to_bits(),
            "independent dyadic collision oracle q{q}"
        );
        for byte in v.to_bits().to_le_bytes() {
            hash = (hash ^ u64::from(byte)).wrapping_mul(0x100000001b3);
        }
    }
    println!(
        "{{\"suite\":\"fs-lbm/d3q27\",\"lane\":\"golden\",\"stage\":\"full-cumulant-product-state\",\"verdict\":\"pass\",\"hash\":\"{hash:016x}\"}}"
    );
    assert_eq!(
        hash, D3Q27_GOLDEN_HASH,
        "new semantic surface; no D3Q19 golden bump"
    );
}

fn shear_phasor(grid: &Grid27) -> fs_math::c64::C64 {
    let n = grid.dimensions()[0];
    let mut p = fs_math::c64::C64::ZERO;
    for x in 0..n {
        let v = grid.macroscopic(x).unwrap().1[1];
        let phase = TAU * x as f64 / n as f64;
        p.re += v * det::cos(phase);
        p.im -= v * det::sin(phase);
    }
    p
}

fn shear_decay(n: usize, cumulant: bool, higher: f64) -> (f64, f64, AcousticScalingPlan) {
    let g = gas();
    let dx = 1.6e-3 / n as f64;
    let scaling = AcousticScalingPlan::from_gas(&g, dx, 0.02 * g.sound_speed, 0.2).unwrap();
    let r = rates(1.0 / scaling.lattice.tau, higher);
    let model = if cumulant {
        CollisionModel27::Cumulant { rates: r }
    } else {
        CollisionModel27::CentralMoment { rates: r }
    };
    let mean = scaling.lattice.u_lattice;
    let amplitude = 1e-4 * CS2_27.sqrt();
    let mut grid = Grid27::try_new(
        [n, 1, 1],
        config(model),
        [0.0; 3],
        [true; 3],
        GridBudget3 {
            maximum_cells: n,
            maximum_mach: 0.2,
        },
        |_| false,
        |p| {
            equilibrium27(
                1.0,
                [
                    mean,
                    amplitude * det::sin(TAU * p[0] as f64 / n as f64),
                    0.0,
                ],
                Equilibrium27::Factorized,
            )
        },
    )
    .unwrap();
    with_cx(|cx| {
        // Remove kinetic start-up before estimating the hydrodynamic pole.
        // 16N² warm-up steps are declared, not fitted to the observed result.
        for _ in 0..16 * n * n {
            grid.step_serial(cx).unwrap();
        }
        let before = shear_phasor(&grid);
        for _ in 0..20 * n {
            grid.step_serial(cx).unwrap();
        }
        let after = shear_phasor(&grid);
        let seconds = 20.0 * n as f64 * scaling.dt_s;
        let decay = -det::ln(after.abs() / before.abs()) / seconds;
        let rotation = after / before;
        let phase_speed = -det::atan2(rotation.im, rotation.re) / (TAU / 1.6e-3 * seconds);
        (decay, phase_speed, scaling)
    })
}

#[test]
fn collision_dissipation_ladder() {
    use fs_ladder::convergence::{ConvergenceEvaluator, ConvergencePlan, MeshRung};
    let start = Instant::now();
    let g = gas();
    let nu = g.dynamic_viscosity / g.density;
    let exact = nu * (TAU / 1.6e-3).powi(2);
    for cumulant in [false, true] {
        let model = if cumulant {
            "Cumulant-Geier2015-full"
        } else {
            "CentralMoment-full"
        };
        let mut by_rate = Vec::new();
        for higher in [0.6, 1.0, 1.4] {
            let mut plan = ConvergencePlan::new("shear-wave-decay-rate", 2.0);
            let mut values = Vec::new();
            let mut errors = Vec::new();
            let mut rungs = Vec::new();
            for (ordinal, n) in [8, 16, 32].into_iter().enumerate() {
                let rung_start = Instant::now();
                let (qoi, phase, scaling) = shear_decay(n, cumulant, higher);
                let err = (qoi / exact - 1.0).abs();
                assert!(qoi.is_finite() && qoi > 0.0 && err.is_finite());
                let mut rung = MeshRung::new(
                    ordinal,
                    format!("D3Q27-{n}"),
                    scaling.dx_m,
                    "m",
                    n * Q27,
                    "shear-wave-decay-rate",
                    qoi,
                    "s^-1",
                );
                // Explicit updates have no iterative algebraic solve residual.
                rung.solver_residual = 0.0;
                rung.budget_consumed_s = rung_start.elapsed().as_secs_f64();
                plan = plan.with_rung(rung);
                let record = format!(
                    "{{\"dx_m\":{},\"dt_s\":{},\"tau\":{},\"qoi\":{qoi},\"err\":{err},\"oracle_s_inverse\":{exact},\"phase_speed_m_s\":{phase},\"phase_oracle_m_s\":{},\"warmup_steps\":{},\"fit_steps\":{},\"collision\":{{\"model\":\"{model}\",\"rates\":{:?}}}}}",
                    scaling.dx_m,
                    scaling.dt_s,
                    scaling.lattice.tau,
                    0.02 * g.sound_speed,
                    16 * n * n,
                    20 * n,
                    rates(1.0 / scaling.lattice.tau, higher).omega
                );
                emit(
                    "dissipation-ladder",
                    "rung",
                    "D3Q27",
                    "pass",
                    &record,
                    rung_start,
                );
                values.push(qoi);
                errors.push(err);
                rungs.push(record);
            }
            let result = ConvergenceEvaluator::evaluate(&plan);
            let order = result
                .observed_order
                .expect("three-rung fit must resolve an order");
            let error_orders = [
                det::ln(errors[0] / errors[1]) / det::ln(2.0),
                det::ln(errors[1] / errors[2]) / det::ln(2.0),
            ];
            let pass = order >= 1.8 && error_orders.iter().all(|&p| p >= 1.8);
            emit(
                "dissipation-ladder",
                "order",
                "D3Q27",
                if pass { "pass" } else { "fail" },
                &format!(
                    "{{\"collision\":{{\"model\":\"{model}\",\"free_higher_rate\":{higher}}},\"rungs\":[{}],\"observed_order\":{order},\"error_orders\":{error_orders:?},\"tolerance_minimum_order\":1.8}}",
                    rungs.join(",")
                ),
                start,
            );
            assert!(
                pass,
                "{model} rate{higher}: apparent order{order}, analytic-error orders{error_orders:?}, errors{errors:?}"
            );
            by_rate.push(values);
        }
        let spread: Vec<_> = (0..3)
            .map(|i| {
                let min = by_rate.iter().map(|r| r[i]).fold(f64::INFINITY, f64::min);
                let max = by_rate
                    .iter()
                    .map(|r| r[i])
                    .fold(f64::NEG_INFINITY, f64::max);
                max - min
            })
            .collect();
        let pass = spread.iter().all(|s| s.is_finite() && *s > 1e-10 * exact)
            && spread[0] > spread[1]
            && spread[1] > spread[2];
        emit(
            "dissipation-ladder",
            "terminal",
            "D3Q27",
            if pass { "pass" } else { "fail" },
            &format!(
                "{{\"collision\":{{\"model\":\"{model}\",\"higher_rate_axis\":[0.6,1,1.4]}},\"rate_spread_per_rung\":{spread:?},\"spread_qoi_unit\":\"s^-1\",\"liveness_floor_relative\":1e-10,\"tolerance\":\"strictly shrinking\"}}"
            ),
            start,
        );
        assert!(
            pass,
            "{model} rate axis must be live and vanishing: {spread:?}"
        );
    }
}

/// Fourier interpolation avoids adding a bilinear interpolation's own m=4
/// error to BOTH lattices and thereby hiding the lattice discrimination.
fn ring_m4(mut velocity: Vec<f64>, n: usize, radius_cells: f64) -> f64 {
    let mut spectrum: Vec<_> = velocity
        .drain(..)
        .map(|v| fs_fft::C64::new(v, 0.0))
        .collect();
    fs_fft::FftNd::new(&[n, n]).forward(&mut spectrum);
    let mut mean = 0.0;
    let mut fourth = [0.0; 2];
    for angle in 0..128 {
        let theta = TAU * angle as f64 / 128.0;
        let x = n as f64 / 2.0 + radius_cells * det::cos(theta);
        let y = n as f64 / 2.0 + radius_cells * det::sin(theta);
        let mut value = 0.0;
        for ky in 0..n {
            for kx in 0..n {
                let sx = if kx < n / 2 {
                    kx as isize
                } else {
                    kx as isize - n as isize
                };
                let sy = if ky < n / 2 {
                    ky as isize
                } else {
                    ky as isize - n as isize
                };
                let phase = TAU * (sx as f64 * x + sy as f64 * y) / n as f64;
                let z = spectrum[ky * n + kx];
                value += z.re * det::cos(phase) - z.im * det::sin(phase);
            }
        }
        value /= (n * n) as f64;
        mean += value;
        fourth[0] += value * det::cos(4.0 * theta);
        fourth[1] -= value * det::sin(4.0 * theta);
    }
    assert!(mean.is_finite() && mean > 0.0);
    2.0 * det::hypot(fourth[0], fourth[1]) / mean
}

fn axial_gaussian(n: usize, dt: f64, dx: f64, u_m_s: f64) -> impl Fn([usize; 3]) -> [f64; 3] {
    move |p| {
        let radius = 0.4e-3;
        let r = |axis: usize| {
            let x = (p[axis] as f64 - n as f64 / 2.0) * dx;
            (-1..=1)
                .map(|image| {
                    let d = (x + f64::from(image) * n as f64 * dx) / radius;
                    det::exp(-0.5 * d * d)
                })
                .sum::<f64>()
        };
        [0.0, 0.0, u_m_s * dt / dx * r(0) * r(1)]
    }
}

fn mean_jet<L: fs_lbm::lattice3::Lattice3<Q>, const Q: usize>(
    mut grid: Grid3<L, Q>,
    steps: usize,
) -> Vec<f64> {
    let mut mean = vec![0.0; grid.populations().len()];
    let mut samples = 0;
    with_cx(|cx| {
        for step in 0..steps {
            grid.step_serial(cx).unwrap();
            if step >= steps * 3 / 4 {
                for (cell, mean) in mean.iter_mut().enumerate() {
                    *mean += grid.macroscopic(cell).unwrap().1[2];
                }
                samples += 1;
            }
        }
    });
    for v in &mut mean {
        *v /= samples as f64;
    }
    mean
}

#[test]
fn rotational_invariance_d3q27_vs_d3q19_control() {
    let start = Instant::now();
    let g = gas();
    let nu = g.dynamic_viscosity / g.density;
    let radius = 0.4e-3;
    let reynolds = 200.0;
    let speed = reynolds * nu / (2.0 * radius);
    let nu_lattice = 0.001;
    let tau = 0.5 + 3.0 * nu_lattice;
    let mut control = Vec::new();
    let mut bgk27 = Vec::new();
    let mut cm = Vec::new();
    let mut full = Vec::new();
    // Four rungs: at tau - 1/2 = 0.003 the central-moment operator resets the
    // fourth-order moments to their equilibrium while the stress moments are
    // over-relaxed, so its anisotropy converges more slowly than the cumulant's.
    // Measured at 1b576cbe5, CM m4/m0 was 0.031, 0.0175 and 0.0062 on n = 16, 32, 64.
    for n in [16, 32, 64, 128] {
        let rung_start = Instant::now();
        let dx = 8.0 * radius / n as f64;
        // Diffusive scaling here isolates rotational incompressible error;
        // the SEPARATE dissipation ladder keeps physical acoustic c fixed.
        let dt = nu_lattice * dx * dx / nu;
        let steps = (0.04 * radius * radius / (nu * dt)).ceil() as usize;
        let budget = GridBudget3 {
            maximum_cells: n * n,
            maximum_mach: 0.2,
        };
        let initial = axial_gaussian(n, dt, dx, speed);
        let d19 = Grid3::<D3q19, 19>::try_new(
            [n, n, 1],
            fs_lbm::d3q19::CollisionModel3::Bgk { tau },
            [0.0; 3],
            [true; 3],
            budget,
            |_| false,
            |p| fs_lbm::equilibrium3(1.0, initial(p)),
        )
        .unwrap();
        let a19 = ring_m4(mean_jet(d19, steps), n, 1.5 * radius / dx);
        control.push(a19);
        let mut measured = Vec::new();
        for (model, form) in [
            (
                CollisionModel27::Bgk {
                    tau,
                    equilibrium: Equilibrium27::Quadratic,
                },
                Equilibrium27::Quadratic,
            ),
            (
                CollisionModel27::CentralMoment {
                    rates: rates(1.0 / tau, 1.0),
                },
                Equilibrium27::Factorized,
            ),
            (
                CollisionModel27::Cumulant {
                    rates: rates(1.0 / tau, 1.0),
                },
                Equilibrium27::Factorized,
            ),
        ] {
            let grid = Grid27::try_new(
                [n, n, 1],
                config(model),
                [0.0; 3],
                [true; 3],
                budget,
                |_| false,
                |p| equilibrium27(1.0, initial(p), form),
            )
            .unwrap();
            measured.push(ring_m4(mean_jet(grid, steps), n, 1.5 * radius / dx));
        }
        bgk27.push(measured[0]);
        cm.push(measured[1]);
        full.push(measured[2]);
        emit(
            "d3q27-invariance",
            "rung",
            "D3Q27-vs-D3Q19",
            "pass",
            &format!(
                "{{\"dx_m\":{dx},\"dt_s\":{dt},\"tau\":{tau},\"steps\":{steps},\"Re\":{reynolds},\"physical_speed_m_s\":{speed},\"lattice_Ma\":{},\"radial_sample_r_over_R\":1.5,\"mean_window\":\"last quarter\",\"m4_over_m0\":{{\"D3Q19_BGK\":{a19},\"D3Q27_BGK\":{},\"D3Q27_CM\":{},\"D3Q27_Cumulant\":{}}},\"collision\":{{\"stress_tau\":{tau},\"higher_rates\":[1,1,1,1,1,1,1,1],\"bulk_rate\":1}},\"scaling\":\"diffusive; no acoustic-speed claim\"}}",
                speed * dt / dx / CS2_27.sqrt(),
                measured[0],
                measured[1],
                measured[2]
            ),
            rung_start,
        );
    }
    // Compare like-for-like BGK polynomial equilibria to isolate the lattice;
    // CM and full cumulant must independently meet the same final m4 band.
    let ratio = control[3] / bgk27[3];
    let discriminate = ratio.is_finite() && ratio >= 2.0;
    let converges = |a: &[f64]| {
        a.iter().all(|v| v.is_finite() && *v > 0.0)
            && a[3] <= 0.005
            && a.windows(2).all(|w| w[0] > w[1])
    };
    let pass = discriminate && converges(&bgk27) && converges(&cm) && converges(&full);
    emit(
        "d3q27-invariance",
        "terminal",
        "D3Q27-vs-D3Q19",
        if !discriminate {
            "refused"
        } else if pass {
            "pass"
        } else {
            "fail"
        },
        &format!(
            "{{\"control_ratio\":{ratio},\"D3Q19_BGK\":{control:?},\"D3Q27_BGK\":{bgk27:?},\"D3Q27_CM\":{cm:?},\"D3Q27_Cumulant\":{full:?},\"tolerance\":{{\"finest_m4_over_m0_max\":0.005,\"control_ratio_min\":2}},\"Re\":{reynolds},\"nondiscriminating_action\":\"refuse fixture; raise and record Re within the explicit Mach envelope\"}}"
        ),
        start,
    );
    assert!(
        discriminate,
        "fixture cannot distinguish lattices: ratio={ratio}; raise/record Re, never call this a pass"
    );
    assert!(
        pass,
        "rotational ladders must decrease within0.5%: BGK{bgk27:?}, CM{cm:?}, full{full:?}"
    );
}
