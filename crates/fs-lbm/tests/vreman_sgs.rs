//! Vreman uses the original lattice collision, streaming and SGS workspace.
//! These are short numerical consistency tests, not turbulent-DNS validation.
use fs_exec::{CancelGate, Cx};
use fs_lbm::d3q19::CollisionModel3;
use fs_lbm::d3q27::{CollisionConfig27, CollisionModel27, D3q27, Rates27};
use fs_lbm::lattice3::{D3q19, Grid3, GridBudget3, Lattice3};
use fs_lbm::sgs3::{SgsConfig, SgsLattice3, SgsModel, SgsUse};
use fs_math::det;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let arenas = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    arenas.scope(|arena| {
        let cx = Cx::new(&gate, arena,
            fs_exec::StreamKey { seed: 0x0535_004, kernel_id: 136, tile: 0, iteration: 0 },
            fs_exec::Budget::INFINITE, fs_exec::ExecMode::Deterministic);
        f(&cx)
    })
}

fn config(cv: f64) -> SgsConfig {
    SgsConfig { model: SgsModel::Vreman { cv }, filter_width_lu: 1.0,
        maximum_eddy_to_molecular: 1000.0, usage: SgsUse::Exploration }
}

fn vortex<L: Lattice3<Q>, const Q: usize>(collision: L::Collision) -> Grid3<L, Q> {
    Grid3::try_new([8; 3], collision, [0.0; 3], [true; 3],
        GridBudget3 { maximum_cells: 512, maximum_mach: 0.2 }, |_| false, |p| {
            let [x, y, z] = p.map(|v| core::f64::consts::TAU * (v as f64 + 0.5) / 8.0);
            let u = [0.02 * det::sin(x) * det::cos(y) * det::cos(z),
                -0.02 * det::cos(x) * det::sin(y) * det::cos(z), 0.0];
            let speed2: f64 = u.iter().map(|v| v * v).sum();
            core::array::from_fn(|q| {
                let eu: f64 = (0..3).map(|a| f64::from(L::E[q][a]) * u[a]).sum();
                L::W[q] * (1.0 + eu / L::CS2 + 0.5 * eu * eu / (L::CS2 * L::CS2)
                    - 0.5 * speed2 / L::CS2)
            })
        }).unwrap()
}

fn bits<L: Lattice3<Q>, const Q: usize>(grid: &Grid3<L, Q>) -> Vec<u64> {
    grid.populations().iter().flatten().map(|x| x.to_bits()).collect()
}

fn zero<L: SgsLattice3<Q>, const Q: usize>(collision: L::Collision) {
    let mut ordinary = vortex::<L, Q>(collision);
    let mut disabled = vortex::<L, Q>(collision);
    let mut workspace = disabled.sgs_workspace().unwrap();
    let allocation = workspace.collisions().as_ptr();
    with_cx(|cx| {
        for _ in 0..4 {
            ordinary.step_serial(cx).unwrap();
            let report = disabled.step_serial_with_sgs(cx, config(0.0), None, &mut workspace).unwrap();
            assert_eq!(report.sgs_dissipation_lu, 0.0);
            assert_eq!(bits(&ordinary), bits(&disabled));
            assert_eq!(allocation, workspace.collisions().as_ptr());
        }
    });
}

#[test]
fn vreman_zero_is_exact_on_both_stencils_and_multiple_collision_models() {
    zero::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.8 });
    zero::<D3q19, 19>(CollisionModel3::CentralMoment {
        second_order_rate: 1.3, higher_order_rate: 1.1,
    });
    zero::<D3q27, 27>(CollisionConfig27 {
        model: CollisionModel27::Cumulant { rates: Rates27 { omega: [1.3; 10] } },
        maximum_mach: 0.2,
    });
}

#[test]
fn nonzero_vreman_changes_shear_preserves_bulk_and_conserves_mass() {
    let rates = Rates27 { omega: [1.3; 10] };
    let collision = CollisionConfig27 {
        model: CollisionModel27::Cumulant { rates }, maximum_mach: 0.2,
    };
    let mut grid = vortex::<D3q27, 27>(collision);
    let mut ordinary = vortex::<D3q27, 27>(collision);
    let mass = |g: &Grid3<D3q27, 27>| g.populations().iter().flatten().sum::<f64>();
    let before = mass(&grid);
    let mut workspace = grid.sgs_workspace().unwrap();
    with_cx(|cx| {
        let report = grid.prepare_sgs(cx, config(0.07), None, &mut workspace).unwrap();
        assert!(report.sgs_dissipation_lu > 0.0);
        let mut changed = false;
        for c in workspace.collisions() {
            let CollisionModel27::Cumulant { rates: total } = c.model else { panic!("collision changed") };
            assert_eq!(total.omega[1..], rates.omega[1..]);
            changed |= total.omega[0] < rates.omega[0];
        }
        assert!(changed);
        for _ in 0..4 {
            grid.step_serial_with_sgs(cx, config(0.07), None, &mut workspace).unwrap();
            ordinary.step_serial(cx).unwrap();
        }
    });
    assert!((mass(&grid) / before - 1.0).abs() < 3e-12);
    assert_ne!(bits(&grid), bits(&ordinary));
}
