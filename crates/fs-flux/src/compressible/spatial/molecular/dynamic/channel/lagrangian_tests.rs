use super::*;
use crate::compressible::spatial::molecular::sgs::{
    LagrangianScalarConfig, LagrangianSgsConfig, LagrangianSubgridConfig,
    LagrangianThermalSubgridConfig, SubgridTransport,
};
use fs_lbm::sgs3::dynamic::lagrangian::scalar::LagrangianDynamicScalar;

fn memory() -> LagrangianThermalSubgridConfig {
    LagrangianThermalSubgridConfig {
        momentum: LagrangianSgsConfig {
            memory_time: 0.01,
            maximum_departure_cells: 1.0,
            maximum_coefficient_squared: 100.0,
            maximum_eddy_to_molecular: 1000.0,
            usage: SgsUse::Exploration,
        },
        heat: LagrangianScalarConfig {
            memory_time: 0.017,
            maximum_departure_cells: 1.0,
            maximum_coefficient: 100.0,
            maximum_diffusivity_to_molecular: 1000.0,
            usage: SgsUse::Exploration,
        },
    }
}

fn physical(flow: &Euler3d) -> (Vec<f64>, Vec<[f64; 3]>, Vec<f64>) {
    let p: Vec<_> = flow
        .states()
        .iter()
        .map(|state| state.primitive(flow.gas(), || false).unwrap())
        .collect();
    (
        p.iter().map(|p| p.density).collect(),
        p.iter().map(|p| p.velocity).collect(),
        p.iter().map(|p| p.temperature).collect(),
    )
}

fn channel_owner(normal: usize, mu: f64, k: f64) -> MolecularEuler3d<ConstantTransport> {
    let flow = channel(normal, true, [0.17, 0.11, 0.13]);
    let source = law(flow.gas(), mu, 0.003, k);
    let bytes = MolecularEuler3d::<ConstantTransport>::required_lagrangian_dynamic_thermal_sgs_workspace_bytes(&flow).unwrap();
    assert_eq!(bytes, 320 * flow.states.len());
    let result = MolecularEuler3d::with_lagrangian_dynamic_channel_thermal_sgs(
        flow,
        source,
        memory(),
        walls(normal, false),
        bytes,
        || false,
    )
    .unwrap();
    assert_eq!(result.lagrangian_channel_normal_axis(), Some(normal));
    result
}

#[test]
fn g1_channel_lagrangian_uses_current_wall_data_and_preserves_molecular_wall_fluxes() {
    for normal in 0..3 {
        let flow = channel(normal, true, [0.17, 0.11, 0.13]);
        let mut active = configured(&flow, normal, walls(normal, false));
        let mut base = configured(&flow, normal, walls(normal, false));
        active.sgs = Some(
            SubgridTransport::new_lagrangian_channel_thermal(
                &flow,
                memory(),
                320 * flow.states.len(),
            )
            .unwrap(),
        );
        let replacement = MolecularWall3d::NoSlip {
            velocity: core::array::from_fn(|axis| if axis == normal { 0.0 } else { 0.07 }),
            thermal: WallThermalCondition::Adiabatic,
        };
        active = active
            .with_wall(BoundarySide::ALL[2 * normal + 1], replacement)
            .unwrap();
        base = base
            .with_wall(BoundarySide::ALL[2 * normal + 1], replacement)
            .unwrap();

        let (rho, velocity, temperature) = physical(&flow);
        let (wall_velocity, scalar_walls) = channel_wall_values(normal, &active.walls).unwrap();
        let mut reference = LagrangianDynamicScalar::new_channel(
            flow.grid.shape,
            flow.grid.spacing,
            normal,
            usize::MAX,
        )
        .unwrap();
        let expected = reference
            .prepare_with_walls(
                memory().momentum,
                memory().heat,
                flow.time().value(),
                &rho,
                &velocity,
                wall_velocity,
                &temperature,
                scalar_walls,
                || false,
            )
            .unwrap();
        let mut filter = active
            .sgs
            .as_mut()
            .unwrap()
            .take_lagrangian_filter()
            .unwrap();
        let trial = active
            .sgs
            .as_mut()
            .unwrap()
            .prepare_lagrangian_with_walls(&mut filter, &flow, &active.walls, &mut || false)
            .unwrap();
        assert_eq!(
            trial.reports(),
            (expected.momentum_report(), Some(expected.scalar_report()))
        );
        let base_rate = prepare(&mut base, &flow);
        assert!(prepare(&mut active, &flow) >= base_rate);
        assert_ne!(
            active.work.flux, base.work.flux,
            "interior modeled transport must be active"
        );
        for cell in 0..flow.states.len() {
            let c = flow.grid.coordinates(cell);
            if c[normal] != 0 {
                continue;
            }
            for end in [0, flow.grid.shape[normal]] {
                let mut boundary = c;
                boundary[normal] = end;
                let face = flow.grid.cell_face(normal, boundary);
                assert_eq!(active.work.flux[normal][face], base.work.flux[normal][face]);
            }
        }
        drop(trial);
        active
            .sgs
            .as_mut()
            .unwrap()
            .restore_lagrangian_filter(filter, None);
        assert_eq!(active.sgs.as_ref().unwrap().lagrangian_sample_time(), None);
    }
}

#[test]
fn g4_channel_physical_start_samples_observation_and_rkl_replay_preserve_both_histories() {
    let make = || channel_owner(2, 2.0, 4.0);
    let mut observed = make();
    let mut control = make();
    let mut history = LagrangianDynamicScalar::new_channel(
        observed.flow().grid.shape,
        observed.flow().grid.spacing,
        2,
        usize::MAX,
    )
    .unwrap();
    for duration in [0.0001, 0.0002] {
        let start = observed.flow().time();
        let (rho, velocity, temperature) = physical(observed.flow());
        let wall = BoundarySide::ALL.map(|side| observed.wall(side));
        let (wall_velocity, scalar_walls) = channel_wall_values(2, &wall).unwrap();
        let trial = history
            .prepare_with_walls(
                memory().momentum,
                memory().heat,
                start.value(),
                &rho,
                &velocity,
                wall_velocity,
                &temperature,
                scalar_walls,
                || false,
            )
            .unwrap();
        let report = observed
            .step(Time::new(duration), MolecularStepControl::default(), || {
                false
            })
            .unwrap();
        control
            .step(Time::new(duration), MolecularStepControl::default(), || {
                false
            })
            .unwrap();
        assert_eq!(observed.lagrangian_subgrid_sample_time(), Some(start));
        assert_eq!(
            observed.lagrangian_subgrid_report(),
            Some(trial.momentum_report())
        );
        assert_eq!(
            observed.lagrangian_subgrid_heat_report(),
            Some(trial.scalar_report())
        );
        for residual in report.conservation_residual.array() {
            close(residual, 0.0, 3e-12);
        }
        trial.commit();
    }
    let before = observed.flow().states().to_vec();
    let before_time = observed.flow().time();
    let before_sample = observed.lagrangian_subgrid_sample_time();
    let before_momentum = observed.lagrangian_subgrid_report();
    let before_heat = observed.lagrangian_subgrid_heat_report();
    let before_report = observed.last_step().copied();
    let before_loads = BoundarySide::ALL.map(|side| observed.boundary_impulses(side).to_vec());
    let mut work = observed.vorticity_workspace(usize::MAX).unwrap();
    let first = observed.vorticity(&mut work, || false).unwrap().to_vec();
    assert_eq!(observed.vorticity(&mut work, || false).unwrap(), first);
    assert_eq!(observed.lagrangian_subgrid_sample_time(), before_sample);
    assert_eq!(observed.lagrangian_subgrid_report(), before_momentum);
    assert_eq!(observed.lagrangian_subgrid_heat_report(), before_heat);
    let mut polls = 0;
    let expected = control
        .step_stabilized(
            Time::new(0.004),
            MolecularStepControl::default(),
            MolecularStabilization::default(),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    assert!(expected.maximum_stages_used > 2);
    let mut calls = 0;
    assert_eq!(
        observed.step_stabilized(
            Time::new(0.004),
            MolecularStepControl::default(),
            MolecularStabilization::default(),
            || {
                calls += 1;
                calls == polls
            },
        ),
        Err(MolecularError::Cancelled)
    );
    assert_eq!(observed.flow().states(), before);
    assert_eq!(observed.flow().time(), before_time);
    assert_eq!(observed.last_step().copied(), before_report);
    assert_eq!(observed.lagrangian_subgrid_sample_time(), before_sample);
    assert_eq!(observed.lagrangian_subgrid_report(), before_momentum);
    assert_eq!(observed.lagrangian_subgrid_heat_report(), before_heat);
    for (i, side) in BoundarySide::ALL.into_iter().enumerate() {
        assert_eq!(observed.boundary_impulses(side), before_loads[i]);
    }
    assert_eq!(
        observed
            .step_stabilized(
                Time::new(0.004),
                MolecularStepControl::default(),
                MolecularStabilization::default(),
                || false,
            )
            .unwrap(),
        expected
    );
    assert_eq!(observed.flow().states(), control.flow().states());
    assert_eq!(observed.lagrangian_subgrid_sample_time(), Some(before_time));
    assert_eq!(
        observed.lagrangian_subgrid_report(),
        control.lagrangian_subgrid_report()
    );
    assert_eq!(
        observed.lagrangian_subgrid_heat_report(),
        control.lagrangian_subgrid_heat_report()
    );
    // A later sample detects an accidentally committed diagnostic or failed trial.
    assert_eq!(
        observed
            .step(Time::new(0.0001), MolecularStepControl::default(), || false)
            .unwrap(),
        control
            .step(Time::new(0.0001), MolecularStepControl::default(), || false)
            .unwrap(),
    );
    assert_eq!(observed.flow().states(), control.flow().states());
}

#[test]
fn g0_channel_lagrangian_budget_geometry_walls_and_zero_fit_are_explicit() {
    let make = || {
        field(
            gas(0.0),
            CartesianGrid::new(3, 3, 3, 0.1, 0.2, 0.3).unwrap(),
            boundaries(1),
            |_| [1.0, 0.0, 0.0, 0.0, 2.0],
        )
    };
    let wall = [MolecularWall3d::NoSlip {
        velocity: [0.0; 3],
        thermal: WallThermalCondition::Adiabatic,
    }; 2];
    let flow = make();
    let source = law(flow.gas(), 0.002, 0.003, 0.004);
    assert!(matches!(
        MolecularEuler3d::with_lagrangian_dynamic_channel_thermal_sgs(
            flow,
            source,
            memory(),
            wall,
            320 * 27 - 1,
            || false,
        ),
        Err(MolecularError::DynamicSubgrid(
            DynamicSgsError::WorkspaceBytes { .. }
        ))
    ));
    let flow = make();
    let source = law(flow.gas(), 0.002, 0.003, 0.004);
    assert!(
        MolecularEuler3d::with_lagrangian_dynamic_channel_thermal_sgs(
            flow,
            source,
            memory(),
            [MolecularWall3d::SlipAdiabatic; 2],
            320 * 27,
            || false,
        )
        .is_err()
    );
    let flow = cube(3, |_| [1.0, 0.0, 0.0, 0.0, 2.0]);
    let source = law(flow.gas(), 0.002, 0.003, 0.004);
    assert!(
        MolecularEuler3d::with_lagrangian_dynamic_channel_thermal_sgs(
            flow,
            source,
            memory(),
            wall,
            usize::MAX,
            || false,
        )
        .is_err()
    );

    let flow = make();
    let source = law(flow.gas(), 0.002, 0.003, 0.004);
    let mut thermal = MolecularEuler3d::with_lagrangian_dynamic_channel_thermal_sgs(
        flow,
        source,
        memory(),
        wall,
        320 * 27,
        || false,
    )
    .unwrap();
    let flow = make();
    let source = law(flow.gas(), 0.002, 0.003, 0.004);
    let fixed_config = LagrangianSubgridConfig {
        model: memory().momentum,
        turbulent_prandtl: 0.8,
    };
    assert_eq!(
        MolecularEuler3d::<ConstantTransport>::required_lagrangian_dynamic_sgs_workspace_bytes(
            &flow
        )
        .unwrap(),
        208 * 27
    );
    let mut fixed = MolecularEuler3d::with_lagrangian_dynamic_channel_sgs(
        flow,
        source,
        fixed_config,
        wall,
        208 * 27,
        || false,
    )
    .unwrap();
    assert_eq!(
        fixed.lagrangian_dynamic_subgrid_config(),
        Some(fixed_config)
    );
    let flow = make();
    let source = law(flow.gas(), 0.002, 0.003, 0.004);
    let mut base = MolecularEuler3d::new(flow, source, || false).unwrap();
    for (end, value) in wall.into_iter().enumerate() {
        base = base.with_wall(BoundarySide::ALL[2 + end], value).unwrap();
    }
    for _ in 0..2 {
        let expected = base
            .step(Time::new(0.001), MolecularStepControl::default(), || false)
            .unwrap();
        assert_eq!(
            thermal
                .step(Time::new(0.001), MolecularStepControl::default(), || false)
                .unwrap(),
            expected
        );
        assert_eq!(
            fixed
                .step(Time::new(0.001), MolecularStepControl::default(), || false)
                .unwrap(),
            expected
        );
        assert_eq!(thermal.flow().states(), base.flow().states());
        assert_eq!(fixed.flow().states(), base.flow().states());
    }
    assert_eq!(
        thermal
            .lagrangian_subgrid_heat_report()
            .unwrap()
            .maximum_coefficient,
        0.0
    );
    assert_eq!(fixed.lagrangian_subgrid_heat_report(), None);
    assert!(
        thermal
            .with_wall(BoundarySide::YMin, MolecularWall3d::SlipAdiabatic)
            .is_err()
    );
}
