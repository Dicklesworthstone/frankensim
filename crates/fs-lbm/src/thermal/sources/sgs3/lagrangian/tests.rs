//! Native source/flow/heat composition and joint-history publication.

use super::*;
use crate::d3q19::equilibrium3;
use crate::sgs3::SgsUse;
use crate::sgs3::dynamic::channel::ScalarWall;
use crate::sgs3::dynamic::lagrangian::scalar::{LagrangianDynamicScalar, LagrangianScalarConfig};
use crate::sgs3::dynamic::lagrangian::{AdaptiveLagrangianMemory, LagrangianSgsConfig};
use crate::thermal::boundary::ThermalWall;

fn config() -> ThermalLagrangianSgsConfig {
    ThermalLagrangianSgsConfig {
        momentum: LagrangianSgsConfig {
            memory_time: 2.0,
            maximum_departure_cells: 1.0,
            maximum_coefficient_squared: 100.0,
            maximum_eddy_to_molecular: 1000.0,
            usage: SgsUse::Exploration,
        },
        heat: LagrangianScalarConfig {
            memory_time: 0.4,
            maximum_departure_cells: 1.0,
            maximum_coefficient: 100.0,
            maximum_diffusivity_to_molecular: 1000.0,
            usage: SgsUse::Exploration,
        },
        maximum_mach: 0.2,
    }
}
fn adaptive() -> AdaptiveLagrangianMemory {
    AdaptiveLagrangianMemory {
        coefficient: 0.05,
        minimum_time: 0.1,
    }
}
fn limits(cells: usize) -> ThermalSourceLimits {
    ThermalSourceLimits {
        max_cells: cells,
        temperature_range: [-0.2, 1.2],
        balance_tolerance: 2e-12,
    }
}
fn bits(s: &ThermalLbm3) -> Vec<u64> {
    s.grid
        .f
        .iter()
        .flatten()
        .chain(s.gpop.iter().flatten())
        .chain(s.grid.fext.iter().flatten())
        .chain(&s.grid.tau)
        .chain([&s.tau_g, &s.t_bottom, &s.t_top, &s.gbeta, &s.t_ref])
        .map(|x| x.to_bits())
        .collect()
}
fn noise(i: usize, salt: u64) -> f64 {
    let mut x = (i as u64 + 1).wrapping_mul(0x9e37_79b9_7f4a_7c15) ^ salt;
    x = (x ^ (x >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    x ^= x >> 31;
    (x >> 11) as f64 / (1u64 << 53) as f64 * 2.0 - 1.0
}
fn seeded() -> ThermalLbm3 {
    let mut state = ThermalLbm3::slab(5, 4, 5, 0.8, 0.9, 0.0);
    for active in 0..100 {
        let cell = active + 20;
        let h = (active / 20) as f64 + 0.5;
        let envelope = 4.0 * h * (5.0 - h) / 25.0;
        let u = core::array::from_fn(|d| 0.02 * envelope * noise(active, 31 * d as u64 + 7));
        let force = [1e-5, -2e-5, 3e-5];
        state.grid.f[cell] = equilibrium3(1.0, u);
        crate::thermal::forcing::reframe3(&mut state.grid.f[cell], u, [0.0; 3], force);
        state.grid.fext[cell] = force;
        let temperature = 1.0 - h / 5.0 + 0.05 * envelope * noise(active, 101);
        let previous: f64 = state.gpop[cell].iter().sum();
        state.gpop[cell][0] += temperature - previous;
    }
    state
}
fn fields(state: &ThermalLbm3) -> (Vec<f64>, Vec<[f64; 3]>, Vec<f64>) {
    let (cells, plane) = state.layout().unwrap();
    let mut rho = Vec::new();
    let mut velocity = Vec::new();
    let mut temperature = Vec::new();
    for cell in plane..cells - plane {
        let (r, u) = state.grid.moments3(cell);
        rho.push(r);
        velocity.push(u);
        temperature.push(state.gpop[cell].iter().sum());
    }
    (rho, velocity, temperature)
}
fn total_rates(
    momentum: &[f64],
    thermal: &[f64],
    eddy: &[f64],
    diffusivity: &[f64],
) -> (Vec<f64>, Vec<f64>) {
    let plane = (momentum.len() - eddy.len()) / 2;
    let mut a = momentum.to_vec();
    let mut b = thermal.to_vec();
    for active in 0..eddy.len() {
        let cell = active + plane;
        if eddy[active] != 0.0 {
            a[cell] = 0.5 + 3.0 * ((a[cell] - 0.5) / 3.0 + eddy[active]);
        }
        if diffusivity[active] != 0.0 {
            b[cell] = 0.5 + 4.0 * ((b[cell] - 0.5) / 4.0 + diffusivity[active]);
        }
    }
    (a, b)
}

#[test]
fn evolving_native_steps_match_joint_history_and_original_transport() {
    let gate = CancelGate::new_clock_free();
    for use_adaptive in [false, true] {
        let mut actual = seeded();
        let mut reference = actual.snapshot().0;
        let n = actual.gpop.len();
        let mut workspace = actual.thermal_lagrangian_sgs_workspace(usize::MAX).unwrap();
        let mut joint =
            LagrangianDynamicScalar::new_channel([5, 4, 5], [1.0; 3], 2, usize::MAX).unwrap();
        if use_adaptive {
            workspace = workspace.with_adaptive_momentum_memory(adaptive()).unwrap();
            joint = joint.with_adaptive_momentum_memory(adaptive()).unwrap();
        }
        let pointer = workspace.relaxation().thermal.as_ptr();
        let original_tau = actual.grid.tau.clone();
        let mut heat_only = false;
        for epoch in 0..3 {
            let time = epoch as f64;
            let momentum: Vec<_> = (0..n)
                .map(|i| 0.75 + 0.002 * (i % 5 + epoch) as f64)
                .collect();
            let thermal: Vec<_> = (0..n)
                .map(|i| 0.85 + 0.002 * (i % 7 + epoch) as f64)
                .collect();
            let (rho, velocity, temperature) = fields(&reference);
            let trial = joint
                .prepare_with_walls(
                    config().momentum,
                    config().heat,
                    time,
                    &rho,
                    &velocity,
                    [[0.0; 3]; 2],
                    &temperature,
                    [ScalarWall::Value(1.0), ScalarWall::Value(0.0)],
                    || false,
                )
                .unwrap();
            heat_only |= trial
                .viscosities()
                .iter()
                .zip(trial.diffusivities())
                .any(|(&nu, &alpha)| nu == 0.0 && alpha > 0.0);
            let fits = (
                trial.momentum_report(),
                trial.scalar_report(),
                trial.momentum_memory_time_range(),
            );
            let (total_momentum, total_thermal) = total_rates(
                &momentum,
                &thermal,
                trial.viscosities(),
                trial.diffusivities(),
            );
            let mut heat = vec![0.0001; n];
            heat[..20].fill(0.0);
            heat[n - 20..].fill(0.0);
            let force = vec![[1e-5 * (time + 1.0), -1e-5, 0.0]; n];
            let expected = reference
                .step_with_sources_and_transport(
                    &heat,
                    &force,
                    RelaxationField {
                        momentum: &total_momentum,
                        thermal: &total_thermal,
                    },
                    limits(n),
                )
                .unwrap();
            trial.commit();
            let result = actual
                .step_with_sources_and_lagrangian_sgs(
                    &heat,
                    &force,
                    time,
                    config(),
                    Some(RelaxationField {
                        momentum: &momentum,
                        thermal: &thermal,
                    }),
                    &mut workspace,
                    limits(n),
                    &gate,
                )
                .unwrap();
            assert_eq!(result.heat, expected);
            assert_eq!(result.sgs.momentum_fit, fits.0);
            assert_eq!(result.sgs.heat_fit, fits.1);
            assert_eq!(result.sgs.momentum_memory_time_range, fits.2);
            assert_eq!(workspace.relaxation().momentum, total_momentum);
            assert_eq!(workspace.relaxation().thermal, total_thermal);
            assert_eq!(bits(&actual), bits(&reference));
            assert_eq!(actual.last_wall_exchange(), reference.last_wall_exchange());
            assert_eq!(workspace.accepted_time(), Some(time));
            assert_eq!(pointer, workspace.relaxation().thermal.as_ptr());
            assert_eq!(actual.grid.tau, original_tau);
            assert_eq!(actual.tau_g, 0.9);
        }
        assert!(
            heat_only,
            "independently positive heat transport must survive a zero momentum fit"
        );
    }
}

#[test]
fn quiescent_history_preserves_original_sources_and_wall_laws_bitwise() {
    let gate = CancelGate::new_clock_free();
    for walls in [
        [ThermalWall::Temperature(1.0), ThermalWall::Temperature(0.0)],
        [ThermalWall::Adiabatic, ThermalWall::Adiabatic],
        [
            ThermalWall::HeatFlux(0.0002),
            ThermalWall::Convective {
                transfer: 0.03,
                ambient: 0.1,
            },
        ],
    ] {
        let mut actual = ThermalLbm3::slab(4, 3, 4, 0.8, 0.9, 0.0);
        actual.set_thermal_walls(walls).unwrap();
        let mut reference = actual.snapshot().0;
        let n = actual.gpop.len();
        let mut work = actual.thermal_lagrangian_sgs_workspace(usize::MAX).unwrap();
        let momentum: Vec<_> = (0..n).map(|i| 0.7 + 0.01 * (i % 5) as f64).collect();
        let thermal: Vec<_> = (0..n).map(|i| 0.8 + 0.01 * (i % 7) as f64).collect();
        let mut heat = vec![0.0001; n];
        heat[..12].fill(0.0);
        heat[n - 12..].fill(0.0);
        let force = vec![[0.0; 3]; n];
        for time in [0.0, 1.0] {
            let molecular = RelaxationField {
                momentum: &momentum,
                thermal: &thermal,
            };
            let expected = reference
                .step_with_sources_and_transport(&heat, &force, molecular, limits(n))
                .unwrap();
            let result = actual
                .step_with_sources_and_lagrangian_sgs(
                    &heat,
                    &force,
                    time,
                    config(),
                    Some(molecular),
                    &mut work,
                    limits(n),
                    &gate,
                )
                .unwrap();
            assert_eq!(result.heat, expected);
            assert_eq!(result.sgs.transport.momentum.maximum_eddy_viscosity_lu, 0.0);
            assert_eq!(result.sgs.transport.maximum_eddy_diffusivity_lu, 0.0);
            assert_eq!(bits(&actual), bits(&reference));
            assert_eq!(actual.last_wall_exchange(), reference.last_wall_exchange());
            assert_eq!(work.accepted_time(), Some(time));
        }
    }
}

#[test]
fn physical_force_frame_and_active_scalar_walls_reach_the_shared_fit() {
    let gate = CancelGate::new_clock_free();
    for (walls, scalar_walls) in [
        (
            [ThermalWall::Temperature(0.4), ThermalWall::Temperature(0.2)],
            [ScalarWall::Value(0.4), ScalarWall::Value(0.2)],
        ),
        (
            [ThermalWall::Adiabatic, ThermalWall::Adiabatic],
            [ScalarWall::ZeroGradient; 2],
        ),
        (
            [
                ThermalWall::HeatFlux(0.0002),
                ThermalWall::Convective {
                    transfer: 0.03,
                    ambient: 0.1,
                },
            ],
            [ScalarWall::Interior; 2],
        ),
    ] {
        let mut forced = seeded();
        forced.set_thermal_walls(walls).unwrap();
        let mut unforced = forced.snapshot().0;
        let n = forced.gpop.len();
        for cell in 20..n - 20 {
            let u = unforced.grid.moments3(cell).1;
            let previous = unforced.grid.fext[cell];
            crate::thermal::forcing::reframe3(&mut unforced.grid.f[cell], u, previous, [0.0; 3]);
            unforced.grid.fext[cell] = [0.0; 3];
        }
        let (rho, velocity, temperature) = fields(&forced);
        let mut joint =
            LagrangianDynamicScalar::new_channel([5, 4, 5], [1.0; 3], 2, usize::MAX).unwrap();
        let expected = joint
            .prepare_with_walls(
                config().momentum,
                config().heat,
                0.0,
                &rho,
                &velocity,
                [[0.0; 3]; 2],
                &temperature,
                scalar_walls,
                || false,
            )
            .unwrap();
        let (momentum, thermal) = total_rates(
            &forced.grid.tau,
            &vec![forced.tau_g; n],
            expected.viscosities(),
            expected.diffusivities(),
        );
        let mut a = forced.thermal_lagrangian_sgs_workspace(usize::MAX).unwrap();
        let mut b = unforced
            .thermal_lagrangian_sgs_workspace(usize::MAX)
            .unwrap();
        let first = forced
            .prepare_thermal_lagrangian_sgs(config(), 0.0, None, &mut a, &gate)
            .unwrap();
        let second = unforced
            .prepare_thermal_lagrangian_sgs(config(), 0.0, None, &mut b, &gate)
            .unwrap();
        assert_eq!(first.relaxation().momentum, momentum);
        assert_eq!(first.relaxation().thermal, thermal);
        for (&x, &y) in first
            .relaxation()
            .momentum
            .iter()
            .zip(second.relaxation().momentum)
            .chain(
                first
                    .relaxation()
                    .thermal
                    .iter()
                    .zip(second.relaxation().thermal),
            )
        {
            assert!(
                (x - y).abs() < 2e-13,
                "force chart changed physical transport: {x} vs {y}"
            );
        }
        drop(first);
        drop(second);
        assert_eq!(a.accepted_time(), None);
        assert_eq!(b.accepted_time(), None);
    }
    let mut thin = ThermalLbm3::slab(3, 3, 2, 0.8, 0.9, 0.0);
    thin.set_thermal_walls([ThermalWall::HeatFlux(0.0002), ThermalWall::Temperature(0.0)])
        .unwrap();
    let mut work = thin.thermal_lagrangian_sgs_workspace(usize::MAX).unwrap();
    assert!(matches!(
        thin.prepare_thermal_lagrangian_sgs(config(), 0.0, None, &mut work, &gate),
        Err(ThermalSgsError::Dynamic(_))
    ));
    assert_eq!(work.accepted_time(), None);
}

#[test]
fn dropped_capped_source_and_cancelled_trials_keep_warmed_memory_and_retry() {
    let gate = CancelGate::new_clock_free();
    let mut state = seeded();
    let mut reference = state.snapshot().0;
    let n = state.gpop.len();
    let bytes = state
        .required_thermal_lagrangian_sgs_workspace_bytes()
        .unwrap();
    assert!(state.thermal_lagrangian_sgs_workspace(bytes - 1).is_err());
    let mut work = state
        .thermal_lagrangian_sgs_workspace(bytes)
        .unwrap()
        .with_adaptive_momentum_memory(adaptive())
        .unwrap();
    let mut control = reference
        .thermal_lagrangian_sgs_workspace(bytes)
        .unwrap()
        .with_adaptive_momentum_memory(adaptive())
        .unwrap();
    assert_eq!(work.workspace_bytes(), bytes);
    let mut heat = vec![0.0; n];
    let force = vec![[0.0; 3]; n];
    for (s, w) in [(&mut state, &mut work), (&mut reference, &mut control)] {
        s.step_with_sources_and_lagrangian_sgs(
            &heat,
            &force,
            0.0,
            config(),
            None,
            w,
            limits(n),
            &gate,
        )
        .unwrap();
    }
    let before = bits(&state);
    let exchange = state.last_wall_exchange();
    let range = work.accepted_momentum_memory_time_range();
    drop(
        state
            .prepare_thermal_lagrangian_sgs(config(), 1.0, None, &mut work, &gate)
            .unwrap(),
    );
    let mut capped = config();
    capped.heat.maximum_diffusivity_to_molecular = 0.0;
    assert!(matches!(
        state.step_with_sources_and_lagrangian_sgs(
            &heat,
            &force,
            1.0,
            capped,
            None,
            &mut work,
            limits(n),
            &gate
        ),
        Err(ThermalSgsError::DiffusivityRatio { .. })
    ));
    heat[30] = 100.0;
    assert!(matches!(
        state.step_with_sources_and_lagrangian_sgs(
            &heat,
            &force,
            1.0,
            config(),
            None,
            &mut work,
            limits(n),
            &gate
        ),
        Err(ThermalSgsError::Source(_))
    ));
    heat[30] = 0.0;
    assert!(
        state
            .prepare_thermal_lagrangian_sgs(config(), 0.0, None, &mut work, &gate)
            .is_err()
    );
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        state.prepare_thermal_lagrangian_sgs(config(), 1.0, None, &mut work, &cancelled),
        Err(ThermalSgsError::Cancelled)
    ));
    assert_eq!(bits(&state), before);
    assert_eq!(state.last_wall_exchange(), exchange);
    assert_eq!(work.accepted_time(), Some(0.0));
    assert_eq!(work.accepted_momentum_memory_time_range(), range);
    let actual = state
        .step_with_sources_and_lagrangian_sgs(
            &heat,
            &force,
            1.0,
            config(),
            None,
            &mut work,
            limits(n),
            &gate,
        )
        .unwrap();
    let expected = reference
        .step_with_sources_and_lagrangian_sgs(
            &heat,
            &force,
            1.0,
            config(),
            None,
            &mut control,
            limits(n),
            &gate,
        )
        .unwrap();
    assert_eq!(actual.heat, expected.heat);
    assert_eq!(actual.sgs, expected.sgs);
    assert_eq!(bits(&state), bits(&reference));
    assert_eq!(state.last_wall_exchange(), reference.last_wall_exchange());
}

#[test]
fn endpoint_mach_failure_preserves_wall_heat_flow_and_history_together() {
    let gate = CancelGate::new_clock_free();
    let mut state = ThermalLbm3::slab(3, 3, 4, 0.8, 0.9, 0.0);
    state
        .set_thermal_walls([
            ThermalWall::HeatFlux(0.0002),
            ThermalWall::Convective {
                transfer: 0.03,
                ambient: 0.1,
            },
        ])
        .unwrap();
    let mut reference = state.snapshot().0;
    let n = state.gpop.len();
    let heat = vec![0.0; n];
    let quiet = vec![[0.0; 3]; n];
    let mut work = state.thermal_lagrangian_sgs_workspace(usize::MAX).unwrap();
    let mut control = reference
        .thermal_lagrangian_sgs_workspace(usize::MAX)
        .unwrap();
    for (s, w) in [(&mut state, &mut work), (&mut reference, &mut control)] {
        s.step_with_sources_and_lagrangian_sgs(
            &heat,
            &quiet,
            0.0,
            config(),
            None,
            w,
            limits(n),
            &gate,
        )
        .unwrap();
    }
    let before = bits(&state);
    let exchange = state.last_wall_exchange();
    let mut strict = config();
    strict.maximum_mach = 0.001;
    assert!(matches!(
        state.step_with_sources_and_lagrangian_sgs(
            &heat,
            &vec![[0.002, 0.0, 0.0]; n],
            1.0,
            strict,
            None,
            &mut work,
            limits(n),
            &gate
        ),
        Err(ThermalSgsError::Mach { .. })
    ));
    assert_eq!(bits(&state), before);
    assert_eq!(state.last_wall_exchange(), exchange);
    assert_eq!(work.accepted_time(), Some(0.0));
    let actual = state
        .step_with_sources_and_lagrangian_sgs(
            &heat,
            &quiet,
            1.0,
            config(),
            None,
            &mut work,
            limits(n),
            &gate,
        )
        .unwrap();
    let expected = reference
        .step_with_sources_and_lagrangian_sgs(
            &heat,
            &quiet,
            1.0,
            config(),
            None,
            &mut control,
            limits(n),
            &gate,
        )
        .unwrap();
    assert_eq!(actual.heat, expected.heat);
    assert_eq!(actual.sgs, expected.sgs);
    assert_eq!(bits(&state), bits(&reference));
    assert_eq!(state.last_wall_exchange(), reference.last_wall_exchange());
}
