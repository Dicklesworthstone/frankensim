//! Physical initial velocity, source-corrected stress and the first real impulse.

use fs_exec::{CancelGate, Cx};
use fs_lbm::d3q19::{CollisionModel3, equilibrium3};
use fs_lbm::d3q27::{
    CollisionConfig27, CollisionModel27, D3q27, Equilibrium27, Rates27, equilibrium27,
};
use fs_lbm::lattice3::sources::ForcedGrid3;
use fs_lbm::lattice3::{D3q19, EquilibriumLattice3, Grid3, GridBudget3, GridError3, Lattice3};

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let arenas = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    arenas.scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            fs_exec::StreamKey {
                seed: 0x1A17,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        ))
    })
}

fn budget() -> GridBudget3 {
    GridBudget3 {
        maximum_cells: 24,
        maximum_mach: 0.2,
    }
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 4e-15,
        "{actual:.17e} != {expected:.17e}"
    );
}

fn physical_moments<L: Lattice3<Q>, const Q: usize>(
    populations: &[f64; Q],
    force: [f64; 3],
) -> (f64, [f64; 3], [[f64; 3]; 3]) {
    let rho = populations.iter().sum::<f64>();
    let velocity: [f64; 3] = core::array::from_fn(|a| {
        (populations
            .iter()
            .zip(L::E)
            .map(|(f, e)| f * f64::from(e[a]))
            .sum::<f64>()
            + 0.5 * force[a])
            / rho
    });
    let stress = core::array::from_fn(|a| {
        core::array::from_fn(|b| {
            populations
                .iter()
                .zip(L::E)
                .map(|(f, e)| f * f64::from(e[a]) * f64::from(e[b]))
                .sum::<f64>()
                + 0.5 * (velocity[a] * force[b] + force[a] * velocity[b])
                - rho * velocity[a] * velocity[b]
        })
    });
    (rho, velocity, stress)
}

fn exercise<L: EquilibriumLattice3<Q>, const Q: usize>(
    collision: L::Collision,
    equilibrium: impl Fn(f64, [f64; 3]) -> [f64; Q],
) {
    let density = 1.7;
    let velocity = [0.041, -0.023, 0.017];
    let force = [0.0004, -0.0003, 0.0002];
    let make_equilibrium = |force| {
        Grid3::<L, Q>::try_new_equilibrium(
            [4, 3, 2],
            collision,
            force,
            [true; 3],
            budget(),
            |_| false,
            |_| (density, velocity),
        )
        .unwrap()
    };
    let zero = make_equilibrium([0.0; 3]);
    for populations in zero.populations() {
        assert_eq!(
            populations.map(f64::to_bits),
            equilibrium(density, velocity).map(f64::to_bits)
        );
    }
    let mut uniform = make_equilibrium(force);
    for (cell, populations) in uniform.populations().iter().enumerate() {
        let (rho, u, stress) = physical_moments::<L, Q>(populations, force);
        near(rho, density);
        let observed = uniform.macroscopic(cell).unwrap();
        near(observed.0, density);
        for a in 0..3 {
            near(u[a], velocity[a]);
            near(observed.1[a], velocity[a]);
            for b in 0..3 {
                near(stress[a][b], if a == b { density * L::CS2 } else { 0.0 });
            }
        }
    }
    with_cx(&CancelGate::new_clock_free(), |cx| {
        uniform.step_serial(cx).unwrap();
        for cell in 0..24 {
            let (rho, u) = uniform.macroscopic(cell).unwrap();
            near(rho, density);
            for a in 0..3 {
                near(rho * u[a], density * velocity[a] + force[a]);
            }
        }
        assert_eq!(uniform.completed_steps(), 1);

        let original = Grid3::<L, Q>::try_new(
            [4, 3, 2],
            collision,
            [0.0; 3],
            [true; 3],
            budget(),
            |p| p == [1, 1, 0],
            |p| {
                let mut f = equilibrium(0.9 + 0.1 * p[0] as f64, velocity);
                // Independent shear and normal nonequilibrium stresses, both
                // with zero mass and momentum. Re-equilibration must fail.
                f[7] += 0.0003;
                f[8] += 0.0003;
                f[9] -= 0.0003;
                f[10] -= 0.0003;
                f[0] -= 0.002;
                f[1] += 0.001;
                f[2] += 0.001;
                f
            },
        )
        .unwrap();
        let before = original.populations().to_vec();
        let forces: Vec<_> = (0..24)
            .map(|cell| {
                if cell % 5 == 0 {
                    [0.0; 3]
                } else {
                    force.map(|f| f * (1.0 + cell as f64 / 24.0))
                }
            })
            .collect();
        let flow = ForcedGrid3::try_new_preserving_velocity(original, cx, forces.clone()).unwrap();
        assert_eq!(flow.forces(), forces);
        assert_eq!(flow.completed_steps(), 0);
        for cell in 0..24 {
            let old = physical_moments::<L, Q>(&before[cell], [0.0; 3]);
            let new = physical_moments::<L, Q>(&flow.populations()[cell], forces[cell]);
            near(new.0, old.0);
            for a in 0..3 {
                near(new.1[a], old.1[a]);
                for b in 0..3 {
                    near(new.2[a][b], old.2[a][b]);
                }
            }
            if forces[cell] == [0.0; 3] {
                assert_eq!(
                    flow.populations()[cell].map(f64::to_bits),
                    before[cell].map(f64::to_bits)
                );
            }
        }
        let bits = zero
            .populations()
            .iter()
            .map(|f| f.map(f64::to_bits))
            .collect::<Vec<_>>();
        let flow = ForcedGrid3::try_new_preserving_velocity(zero, cx, vec![[0.0; 3]; 24]).unwrap();
        assert_eq!(
            flow.populations()
                .iter()
                .map(|f| f.map(f64::to_bits))
                .collect::<Vec<_>>(),
            bits
        );
    });
}

#[test]
fn native_initial_velocity_stress_and_first_step_are_force_consistent() {
    for model in [
        CollisionModel3::Bgk { tau: 0.9 },
        CollisionModel3::Trt {
            tau: 1.2,
            magic: 3.0 / 16.0,
        },
        CollisionModel3::CentralMoment {
            second_order_rate: 1.1,
            higher_order_rate: 1.6,
        },
        CollisionModel3::CentralMomentBulk {
            shear_rate: 1.1,
            bulk_rate: 0.9,
            higher_order_rate: 1.5,
        },
        CollisionModel3::ReducedCumulant {
            second_order_rate: 1.2,
            third_order_rate: 1.1,
            fourth_order_rate: 0.7,
        },
    ] {
        exercise::<D3q19, 19>(model, equilibrium3);
    }
    for model in [
        CollisionModel27::Bgk {
            tau: 0.9,
            equilibrium: Equilibrium27::Quadratic,
        },
        CollisionModel27::Bgk {
            tau: 0.9,
            equilibrium: Equilibrium27::Factorized,
        },
        CollisionModel27::Trt {
            tau: 1.2,
            magic: 3.0 / 16.0,
            equilibrium: Equilibrium27::Quadratic,
        },
        CollisionModel27::Trt {
            tau: 1.2,
            magic: 3.0 / 16.0,
            equilibrium: Equilibrium27::Factorized,
        },
        CollisionModel27::CentralMoment {
            rates: Rates27 { omega: [1.2; 10] },
        },
        CollisionModel27::Cumulant {
            rates: Rates27 { omega: [1.2; 10] },
        },
    ] {
        let form = match model {
            CollisionModel27::Bgk { equilibrium, .. }
            | CollisionModel27::Trt { equilibrium, .. } => equilibrium,
            _ => Equilibrium27::Factorized,
        };
        exercise::<D3q27, 27>(
            CollisionConfig27 {
                model,
                maximum_mach: 0.2,
            },
            |rho, u| equilibrium27(rho, u, form),
        );
    }
}

#[test]
fn factorized_initialization_retains_all_product_source_moments() {
    let rho = 1.3;
    let u = [0.062, -0.041, 0.029];
    let force = [0.003, -0.002, 0.001];
    for model in [
        CollisionModel27::Bgk {
            tau: 0.9,
            equilibrium: Equilibrium27::Factorized,
        },
        CollisionModel27::Trt {
            tau: 0.9,
            magic: 3.0 / 16.0,
            equilibrium: Equilibrium27::Factorized,
        },
    ] {
        let grid = Grid3::<D3q27, 27>::try_new_equilibrium(
            [1; 3],
            CollisionConfig27 {
                model,
                maximum_mach: 0.2,
            },
            force,
            [true; 3],
            budget(),
            |_| false,
            |_| (rho, u),
        )
        .unwrap();
        for index in 0..27 {
            let powers = [index % 3, index / 3 % 3, index / 9];
            // Independent tensor-moment oracle: M0=1, M1=u, M2=cs²+u².
            let moments: [f64; 3] = core::array::from_fn(|a| match powers[a] {
                0 => 1.0,
                1 => u[a],
                _ => 1.0 / 3.0 + u[a] * u[a],
            });
            let derivative: [f64; 3] = core::array::from_fn(|a| match powers[a] {
                0 => 0.0,
                1 => 1.0,
                _ => 2.0 * u[a],
            });
            let source: f64 = (0..3)
                .map(|a| force[a] * derivative[a] * moments[(a + 1) % 3] * moments[(a + 2) % 3])
                .sum();
            let expected = rho * moments.iter().product::<f64>() - 0.5 * source;
            let actual = grid.populations()[0]
                .iter()
                .zip(D3q27::E)
                .map(|(f, e)| {
                    f * (0..3)
                        .map(|a| f64::from(e[a]).powi(powers[a]))
                        .product::<f64>()
                })
                .sum();
            near(actual, expected);
        }
    }
}

#[test]
fn admission_uses_physical_velocity_and_both_declared_mach_limits() {
    let model = CollisionModel3::Bgk { tau: 0.9 };
    let force = [0.4, -0.3, 0.2];
    let initial = || {
        Grid3::<D3q19, 19>::try_new_equilibrium(
            [1; 3],
            model,
            [0.0; 3],
            [true; 3],
            budget(),
            |_| false,
            |_| (1.2, [0.0; 3]),
        )
        .unwrap()
    };
    let grid = Grid3::<D3q19, 19>::try_new_equilibrium(
        [1; 3],
        model,
        force,
        [true; 3],
        budget(),
        |_| false,
        |_| (1.2, [0.0; 3]),
    )
    .unwrap();
    for v in grid.macroscopic(0).unwrap().1 {
        near(v, 0.0);
    }
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let flow = ForcedGrid3::try_new_preserving_velocity(initial(), cx, vec![force]).unwrap();
        for v in flow.macroscopic(0).unwrap().1 {
            near(v, 0.0);
        }
    });
    // The stored raw speed exceeds the Mach cap above; no elapsed step is
    // claimed. The physical velocity, not that fictitious speed, is admitted.
    let collision = CollisionConfig27 {
        model: CollisionModel27::Bgk {
            tau: 0.9,
            equilibrium: Equilibrium27::Factorized,
        },
        maximum_mach: 0.1,
    };
    assert!(matches!(
        Grid3::<D3q27, 27>::try_new_equilibrium(
            [1; 3],
            collision,
            [0.0; 3],
            [true; 3],
            budget(),
            |_| false,
            |_| (1.0, [0.08, 0.0, 0.0]),
        ),
        Err(GridError3::Mach { maximum: 0.1, .. })
    ));
    let raw = Grid3::<D3q27, 27>::try_new(
        [1; 3],
        collision,
        [0.0; 3],
        [true; 3],
        budget(),
        |_| false,
        |_| equilibrium27(1.0, [0.08, 0.0, 0.0], Equilibrium27::Factorized),
    )
    .unwrap();
    with_cx(&CancelGate::new_clock_free(), |cx| {
        assert!(matches!(
            ForcedGrid3::try_new_preserving_velocity(raw, cx, vec![[0.0; 3]]),
            Err(GridError3::Mach { maximum: 0.1, .. })
        ));
    });
}

#[test]
fn invalid_initial_states_and_source_attachment_return_no_owner() {
    let collision = CollisionModel3::Bgk { tau: 0.9 };
    for (rho, u) in [
        (0.0, [0.0; 3]),
        (-1.0, [0.0; 3]),
        (f64::NAN, [0.0; 3]),
        (f64::INFINITY, [0.0; 3]),
        (1.0, [f64::NAN, 0.0, 0.0]),
        (1.0, [f64::INFINITY, 0.0, 0.0]),
        (1.0, [0.3, 0.0, 0.0]),
    ] {
        assert!(
            Grid3::<D3q19, 19>::try_new_equilibrium(
                [2, 2, 1],
                collision,
                [0.0; 3],
                [true; 3],
                budget(),
                |_| false,
                |p| if p == [1, 1, 0] {
                    (rho, u)
                } else {
                    (1.0, [0.0; 3])
                },
            )
            .is_err()
        );
    }
    let make = |force| {
        Grid3::<D3q19, 19>::try_new_equilibrium(
            [2, 2, 1],
            collision,
            force,
            [true; 3],
            budget(),
            |p| p == [1, 1, 0],
            |_| (1.0, [0.01, -0.02, 0.03]),
        )
        .unwrap()
    };
    with_cx(&CancelGate::new_clock_free(), |cx| {
        for forces in [
            vec![[0.0; 3]; 3],
            vec![[0.0; 3], [f64::NAN, 0.0, 0.0], [0.0; 3], [0.0; 3]],
            vec![[0.0; 3], [0.0; 3], [0.0; 3], [0.001, 0.0, 0.0]],
        ] {
            assert!(ForcedGrid3::try_new_preserving_velocity(make([0.0; 3]), cx, forces).is_err());
        }
        assert!(
            ForcedGrid3::try_new_preserving_velocity(
                make([0.001, 0.0, 0.0]),
                cx,
                vec![[0.0; 3]; 4]
            )
            .is_err()
        );
        let mut advanced = make([0.0; 3]);
        advanced.step_serial(cx).unwrap();
        assert!(ForcedGrid3::try_new_preserving_velocity(advanced, cx, vec![[0.0; 3]; 4]).is_err());
    });
    let gate = CancelGate::new_clock_free();
    gate.request();
    with_cx(&gate, |cx| {
        assert!(matches!(
            ForcedGrid3::try_new_preserving_velocity(make([0.0; 3]), cx, vec![[0.0; 3]; 4]),
            Err(GridError3::Cancelled)
        ));
    });
}
