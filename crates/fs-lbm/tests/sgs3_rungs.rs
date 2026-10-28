//! G0/G1/G4/G5 exploration controls for the explicit SGS owner. The TGV
//! ladder tests vanishing closure on an analytic resolved field and short
//! fixed-Re evolution, not a Re1600 turbulent-DNS dissipation benchmark.

use fs_exec::{CancelGate, Cx, TilePool};
use fs_lbm::d3q19::CollisionModel3;
use fs_lbm::d3q27::{
    CollisionConfig27, CollisionModel27, D3q27, Equilibrium27, Rates27, equilibrium27,
};
use fs_lbm::lattice3::{D3q19, Grid3, GridBudget3, GridError3, Lattice3};
use fs_lbm::sgs3::{
    SgsConfig, SgsError, SgsLattice3, SgsModel, SgsStepError, SgsUse, eddy_viscosity,
};
use fs_math::det;

const SEED: u64 = 0x0535_003;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let arenas = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    arenas.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: SEED,
                kernel_id: 136,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn config(model: SgsModel) -> SgsConfig {
    SgsConfig {
        model,
        filter_width_lu: 1.0,
        maximum_eddy_to_molecular: 1000.0,
        usage: SgsUse::Exploration,
    }
}

fn bgk(tau: f64) -> CollisionConfig27 {
    CollisionConfig27 {
        model: CollisionModel27::Bgk {
            tau,
            equilibrium: Equilibrium27::Factorized,
        },
        maximum_mach: 0.2,
    }
}

fn tgv<L: Lattice3<Q>, const Q: usize>(n: usize, collision: L::Collision) -> Grid3<L, Q> {
    Grid3::try_new(
        [n; 3],
        collision,
        [0.0; 3],
        [true; 3],
        GridBudget3 {
            maximum_cells: n * n * n,
            maximum_mach: 0.2,
        },
        |_| false,
        |p| {
            let [x, y, z] = p.map(|v| core::f64::consts::TAU * (v as f64 + 0.5) / n as f64);
            let u = [
                0.02 * det::sin(x) * det::cos(y) * det::cos(z),
                -0.02 * det::cos(x) * det::sin(y) * det::cos(z),
                0.0,
            ];
            let speed2: f64 = u.iter().map(|v| v * v).sum();
            core::array::from_fn(|q| {
                let eu: f64 = (0..3).map(|a| f64::from(L::E[q][a]) * u[a]).sum();
                L::W[q]
                    * (1.0 + eu / L::CS2 + 0.5 * eu * eu / (L::CS2 * L::CS2)
                        - 0.5 * speed2 / L::CS2)
            })
        },
    )
    .unwrap()
}

fn bits<L: Lattice3<Q>, const Q: usize>(grid: &Grid3<L, Q>) -> Vec<u64> {
    grid.populations()
        .iter()
        .flatten()
        .map(|v| v.to_bits())
        .collect()
}

fn zero_control<L: SgsLattice3<Q>, const Q: usize>(collision: L::Collision) {
    let mut reference = tgv::<L, Q>(8, collision);
    let mut controlled = tgv::<L, Q>(8, collision);
    let mut workspace = controlled.sgs_workspace().unwrap();
    let allocation = workspace.collisions().as_ptr();
    with_cx(|cx| {
        for model in [
            SgsModel::Smagorinsky { cs: 0.0 },
            SgsModel::Wale { cw: 0.0 },
        ] {
            for _ in 0..3 {
                reference.step_serial(cx).unwrap();
                let report = controlled
                    .step_serial_with_sgs(cx, config(model), None, &mut workspace)
                    .unwrap();
                assert_eq!(report.sgs_dissipation_lu, 0.0);
                assert_eq!(bits(&reference), bits(&controlled));
                assert_eq!(allocation, workspace.collisions().as_ptr());
            }
        }
    });
}

#[test]
fn sgs_zero_constant_is_bitwise_identity() {
    zero_control::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.8 });
    zero_control::<D3q19, 19>(CollisionModel3::CentralMoment {
        second_order_rate: 1.3,
        higher_order_rate: 1.1,
    });
    zero_control::<D3q27, 27>(bgk(0.8));
    for model in [
        CollisionModel27::CentralMoment {
            rates: Rates27 { omega: [1.3; 10] },
        },
        CollisionModel27::Cumulant {
            rates: Rates27 { omega: [1.3; 10] },
        },
    ] {
        zero_control::<D3q27, 27>(CollisionConfig27 {
            model,
            maximum_mach: 0.2,
        });
    }
}

#[test]
fn halfway_wall_gradients_match_quadratic_poiseuille_flow() {
    let ny = 8;
    let nu = 0.01;
    let mut rates = Rates27 { omega: [1.1; 10] };
    rates.omega[0] = 1.0 / (0.5 + 3.0 * nu);
    let collision = CollisionConfig27 {
        model: CollisionModel27::Cumulant { rates },
        maximum_mach: 0.2,
    };
    // One solid layer on either side makes the internal-solid and external
    // halfway-wall reconstructions exercise the same physical distance.
    for internal_wall in [false, true] {
        let height = if internal_wall { ny + 2 } else { ny };
        let grid = Grid3::<D3q27, 27>::try_new(
            [3, height, 3],
            collision,
            [0.0; 3],
            [true, false, true],
            GridBudget3 {
                maximum_cells: 3 * height * 3,
                maximum_mach: 0.2,
            },
            |p| internal_wall && (p[1] == 0 || p[1] + 1 == height),
            |p| {
                let y = p[1] as f64 + 0.5 - if internal_wall { 1.0 } else { 0.0 };
                equilibrium27(
                    1.0,
                    [0.0001 * y * (ny as f64 - y), 0.0, 0.0],
                    Equilibrium27::Factorized,
                )
            },
        )
        .unwrap();
        let mut workspace = grid.sgs_workspace().unwrap();
        with_cx(|cx| {
            grid.prepare_sgs(
                cx,
                config(SgsModel::Smagorinsky { cs: 0.2 }),
                None,
                &mut workspace,
            )
            .unwrap();
        });
        for y in 0..ny {
            let at = y + usize::from(internal_wall);
            let cell = (height + at) * 3 + 1;
            let updated = workspace.collisions()[cell];
            let CollisionModel27::Cumulant { rates: total } = updated.model else {
                panic!("model changed");
            };
            assert_eq!(total.omega[1..], rates.omega[1..]);
            let observed = (1.0 / total.omega[0] - 0.5) / 3.0 - nu;
            let derivative = 0.0001 * (ny as f64 - 2.0 * (y as f64 + 0.5));
            let expected = 0.2 * 0.2 * derivative.abs();
            assert!(
                (observed - expected).abs() < 2e-16,
                "y={y} got={observed} expected={expected}"
            );
        }
    }
}

#[test]
fn sgs_contribution_vanishes_under_refinement() {
    let mut ratios = Vec::new();
    for n in [8, 16, 32] {
        // Fixed U/(k nu)=1600, with acoustic dt proportional to dx. The
        // lattice molecular viscosity grows with n, while the filter remains
        // one lattice cell; its physical width therefore decreases as 1/n.
        let nu = 0.02 * n as f64 / (core::f64::consts::TAU * 1600.0);
        let mut grid = tgv::<D3q27, 27>(n, bgk(0.5 + 3.0 * nu));
        let mut workspace = grid.sgs_workspace().unwrap();
        let settings = config(SgsModel::Smagorinsky { cs: 0.12 });
        let initial = with_cx(|cx| {
            grid.prepare_sgs(cx, settings, None, &mut workspace)
                .unwrap()
        });
        ratios.push(initial.dissipation_ratio());
        let mass_before: f64 = grid
            .populations()
            .iter()
            .map(|f| f.iter().sum::<f64>())
            .sum();
        with_cx(|cx| {
            for _ in 0..n / 4 {
                let report = grid
                    .step_serial_with_sgs(cx, settings, None, &mut workspace)
                    .unwrap();
                assert!(report.resolved_dissipation_lu > 0.0);
                assert!(report.sgs_dissipation_lu > 0.0);
            }
        });
        let mass_after: f64 = grid
            .populations()
            .iter()
            .map(|f| f.iter().sum::<f64>())
            .sum();
        assert!((mass_after / mass_before - 1.0).abs() < 2e-12);
    }
    for pair in ratios.windows(2) {
        assert!(pair[1] < pair[0]);
        let order = (pair[0] / pair[1]).ln() / 2.0_f64.ln();
        assert!(
            (order - 2.0).abs() <= 0.3,
            "ratios={ratios:?} order={order}"
        );
    }
}

#[test]
fn refusals_are_typed_and_leave_populations_unchanged() {
    let mut grid = tgv::<D3q27, 27>(8, bgk(0.8));
    let original = bits(&grid);
    let mut workspace = grid.sgs_workspace().unwrap();
    let settings = config(SgsModel::Wale { cw: 0.5 });
    let report = with_cx(|cx| {
        grid.prepare_sgs(cx, settings, None, &mut workspace)
            .unwrap()
    });
    let exact = SgsConfig {
        maximum_eddy_to_molecular: report.maximum_eddy_to_molecular,
        ..settings
    };
    with_cx(|cx| grid.prepare_sgs(cx, exact, None, &mut workspace).unwrap());
    let exceeded = SgsConfig {
        maximum_eddy_to_molecular: exact.maximum_eddy_to_molecular * (1.0 - 1e-12),
        ..settings
    };
    with_cx(|cx| {
        assert!(matches!(
            grid.step_serial_with_sgs(cx, exceeded, None, &mut workspace),
            Err(SgsStepError::Closure(
                SgsError::ViscosityRatioExceeded { .. }
            ))
        ));
        let acceptance = SgsConfig {
            usage: SgsUse::AcceptanceRender,
            ..settings
        };
        assert!(matches!(
            grid.step_serial_with_sgs(cx, acceptance, None, &mut workspace),
            Err(SgsStepError::Closure(SgsError::AcceptanceRenderForbidden))
        ));
        assert!(
            grid.prepare_sgs(cx, settings, Some(&[bgk(0.8)]), &mut workspace)
                .is_err()
        );
    });
    for cs in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(eddy_viscosity(SgsModel::Smagorinsky { cs }, 1.0, [[0.0; 3]; 3]).is_err());
    }
    assert!(eddy_viscosity(settings.model, 0.0, [[0.0; 3]; 3]).is_err());
    assert!(eddy_viscosity(settings.model, 1.0, [[f64::NAN; 3]; 3]).is_err());
    assert_eq!(bits(&grid), original);
    assert_eq!(grid.completed_steps(), 0);
}

#[test]
fn sgs_bitwise_worker_count_invariance() {
    let mut reference = tgv::<D3q27, 27>(6, bgk(0.8));
    let settings = config(SgsModel::Wale { cw: 0.5 });
    let mut workspace = reference.sgs_workspace().unwrap();
    let mut reports = Vec::new();
    with_cx(|cx| {
        for _ in 0..3 {
            reports.push(
                reference
                    .step_serial_with_sgs(cx, settings, None, &mut workspace)
                    .unwrap(),
            );
        }
    });
    for workers in [1, 2, 4] {
        let pool = TilePool::for_host(workers, SEED);
        let gate = CancelGate::new_clock_free();
        let mut grid = tgv::<D3q27, 27>(6, bgk(0.8));
        let mut workspace = grid.sgs_workspace().unwrap();
        for expected in &reports {
            let actual = grid
                .step_pooled_with_sgs(&pool, &gate, settings, None, &mut workspace)
                .unwrap();
            assert_eq!(&actual, expected);
        }
        assert_eq!(bits(&grid), bits(&reference), "workers={workers}");
    }
}

#[test]
fn cancelled_sgs_step_can_be_retried() {
    let mut grid = tgv::<D3q27, 27>(4, bgk(0.8));
    let original = bits(&grid);
    let mut workspace = grid.sgs_workspace().unwrap();
    let pool = TilePool::for_host(2, SEED);
    let gate = CancelGate::new_clock_free();
    gate.request();
    let settings = config(SgsModel::Smagorinsky { cs: 0.17 });
    assert!(matches!(
        grid.step_pooled_with_sgs(&pool, &gate, settings, None, &mut workspace),
        Err(SgsStepError::Grid(GridError3::Cancelled))
    ));
    assert_eq!(bits(&grid), original);
    assert_eq!(grid.completed_steps(), 0);
    let fresh = CancelGate::new_clock_free();
    grid.step_pooled_with_sgs(&pool, &fresh, settings, None, &mut workspace)
        .unwrap();
    assert_eq!(grid.completed_steps(), 1);
}
