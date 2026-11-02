use super::*;
use crate::d3q19::CollisionModel3;
use crate::d3q27::{CollisionConfig27, CollisionModel27, D3q27, Equilibrium27};
use crate::lattice3::{D3q19, GridBudget3};
use crate::sgs3::SgsUse;
use crate::sgs3::dynamic::DynamicAveraging;
use fs_exec::TilePool;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let arenas = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    arenas.scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 93,
                kernel_id: 136,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        ))
    })
}

fn config() -> DynamicSgsConfig {
    DynamicSgsConfig {
        averaging: DynamicAveraging::LocalBox,
        maximum_coefficient_squared: 10.0,
        maximum_eddy_to_molecular: 100.0,
        usage: SgsUse::Exploration,
    }
}

fn grid<L: SgsLattice3<Q>, const Q: usize>(collision: L::Collision, varying: bool) -> Grid3<L, Q> {
    Grid3::try_new(
        [7; 3],
        collision,
        [0.0; 3],
        [true; 3],
        GridBudget3 {
            maximum_cells: 343,
            maximum_mach: 0.2,
        },
        |_| false,
        |p| {
            let [x, y, z] = p.map(|i| core::f64::consts::TAU * (i as f64 + 0.5) / 7.0);
            let u = if varying {
                [
                    0.003 * x.sin() + 0.002 * (y + z).sin() + 0.001 * (2.0 * x - y).cos(),
                    0.002 * y.cos() + 0.003 * (z + x).sin() + 0.001 * (x - 2.0 * z).sin(),
                    0.002 * z.sin() + 0.001 * (x + y).cos() + 0.003 * (2.0 * y - z).sin(),
                ]
            } else {
                [0.01, -0.005, 0.003]
            };
            let rho = if varying { 1.0 + 0.002 * x.sin() } else { 1.0 };
            let speed2 = u.iter().map(|v| v * v).sum::<f64>();
            core::array::from_fn(|q| {
                let eu = (0..3).map(|a| f64::from(L::E[q][a]) * u[a]).sum::<f64>();
                rho * L::W[q]
                    * (1.0 + eu / L::CS2 + 0.5 * eu * eu / (L::CS2 * L::CS2)
                        - 0.5 * speed2 / L::CS2)
            })
        },
    )
    .unwrap()
}

fn zero<L: SgsLattice3<Q>, const Q: usize>(collision: L::Collision) {
    let mut original = grid::<L, Q>(collision, false);
    let mut dynamic = grid::<L, Q>(collision, false);
    let mut work = dynamic.dynamic_sgs_workspace(usize::MAX).unwrap();
    with_cx(|cx| {
        for _ in 0..3 {
            original.step_serial(cx).unwrap();
            let report = dynamic
                .step_serial_with_dynamic_sgs(cx, config(), None, &mut work)
                .unwrap();
            assert_eq!(report.fit.maximum_eddy_viscosity, 0.0);
            assert_eq!(report.dissipation.sgs_dissipation_lu, 0.0);
            assert_eq!(original.populations(), dynamic.populations());
        }
    });
}

#[test]
fn dynamic_uniform_velocity_preserves_both_original_lattices() {
    zero::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.8 });
    zero::<D3q27, 27>(CollisionConfig27 {
        model: CollisionModel27::Bgk {
            tau: 0.8,
            equilibrium: Equilibrium27::Factorized,
        },
        maximum_mach: 0.2,
    });
}

#[test]
fn dynamic_actual_collision_field_and_pooled_step_agree() {
    let law = CollisionModel3::Bgk { tau: 0.8 };
    let mut serial = grid::<D3q19, 19>(law, true);
    let mut prepared = grid::<D3q19, 19>(law, true);
    let mut pooled = grid::<D3q19, 19>(law, true);
    let mut work = serial.dynamic_sgs_workspace(usize::MAX).unwrap();
    let mut parallel = pooled.dynamic_sgs_workspace(usize::MAX).unwrap();
    with_cx(|cx| {
        let report = serial
            .prepare_dynamic_sgs(cx, config(), None, &mut work)
            .unwrap();
        assert!(report.fit.maximum_eddy_viscosity > 0.0);
        for (c, nu_t) in work.collisions().iter().zip(work.eddy_viscosities()) {
            let nu = D3q19::molecular_viscosity(*c).unwrap();
            assert!((nu - (0.8 - 0.5) / 3.0 - nu_t).abs() < 1e-16);
        }
        prepared
            .step_serial_with_collisions(cx, work.collisions())
            .unwrap();
        let accepted = serial
            .step_serial_with_dynamic_sgs(cx, config(), None, &mut work)
            .unwrap();
        assert_eq!(report, accepted);
        let pool = TilePool::for_host(2, 93);
        let gate = CancelGate::new_clock_free();
        let parallel_report = pooled
            .step_pooled_with_dynamic_sgs(&pool, &gate, config(), None, &mut parallel)
            .unwrap();
        assert_eq!(accepted, parallel_report);
        assert_eq!(serial.populations(), prepared.populations());
        assert_eq!(serial.populations(), pooled.populations());
    });
}

#[test]
fn dynamic_refused_candidate_retains_original_populations_and_clock() {
    let mut flow = grid::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.8 }, true);
    assert!(matches!(
        flow.dynamic_sgs_workspace(1),
        Err(SgsStepError::Dynamic(
            DynamicSgsError::WorkspaceBytes { .. }
        ))
    ));
    let mut work = flow.dynamic_sgs_workspace(usize::MAX).unwrap();
    let before = flow.populations().to_vec();
    let clock = flow.completed_steps;
    with_cx(|cx| {
        let refused = DynamicSgsConfig {
            maximum_coefficient_squared: 0.0,
            ..config()
        };
        assert!(matches!(
            flow.step_serial_with_dynamic_sgs(cx, refused, None, &mut work),
            Err(SgsStepError::Dynamic(DynamicSgsError::Coefficient { .. }))
        ));
        assert_eq!(flow.populations(), before);
        assert_eq!(flow.completed_steps, clock);
        let refused = DynamicSgsConfig {
            maximum_eddy_to_molecular: 0.0,
            ..config()
        };
        assert!(matches!(
            flow.step_serial_with_dynamic_sgs(cx, refused, None, &mut work),
            Err(SgsStepError::Closure(
                SgsError::ViscosityRatioExceeded { .. }
            ))
        ));
        assert_eq!(flow.populations(), before);
        assert_eq!(flow.completed_steps, clock);
    });
    flow.periodic[0] = false;
    assert!(matches!(
        flow.dynamic_sgs_workspace(usize::MAX),
        Err(SgsStepError::Dynamic(DynamicSgsError::Invalid(_)))
    ));
    flow.periodic[0] = true;
    flow.solid[0] = true;
    assert!(matches!(
        flow.dynamic_sgs_workspace(usize::MAX),
        Err(SgsStepError::Dynamic(DynamicSgsError::Invalid(_)))
    ));
}
