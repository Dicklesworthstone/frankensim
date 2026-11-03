use super::*;
use crate::compressible::spatial::molecular::sgs::{
    ChannelDynamicScalarConfig, ChannelDynamicSgsConfig, ChannelDynamicThermalSubgridConfig,
    channel_wall_values,
};
use fs_lbm::sgs3::dynamic::channel::{ChannelDynamicScalar, ChannelDynamicSgs};

fn settings() -> ChannelDynamicThermalSubgridConfig {
    ChannelDynamicThermalSubgridConfig {
        momentum: ChannelDynamicSgsConfig {
            maximum_coefficient_squared: 100.0,
            maximum_eddy_to_molecular: 1000.0,
            usage: SgsUse::Exploration,
        },
        heat: ChannelDynamicScalarConfig {
            maximum_coefficient: 100.0,
            maximum_diffusivity_to_molecular: 1000.0,
            usage: SgsUse::Exploration,
        },
    }
}

fn boundaries(normal: usize) -> Boundaries3d {
    let mut pairs = [AxisBoundaries::Periodic; 3];
    pairs[normal] = WALL;
    Boundaries3d {
        x: pairs[0],
        y: pairs[1],
        z: pairs[2],
    }
}

fn boost(normal: usize) -> [f64; 3] {
    core::array::from_fn(|axis| {
        if axis == normal {
            0.0
        } else {
            0.01 * (axis + 1) as f64
        }
    })
}

fn walls(normal: usize, adiabatic: bool) -> [MolecularWall3d; 2] {
    [3.0, 2.0].map(|t| MolecularWall3d::NoSlip {
        velocity: boost(normal),
        thermal: if adiabatic {
            WallThermalCondition::Adiabatic
        } else {
            WallThermalCondition::Isothermal(Temperature::new(t))
        },
    })
}

fn noise(i: usize, salt: u64) -> f64 {
    let mut h = ((i as u64 + 1).wrapping_mul(0x9e37_79b9_7f4a_7c15)) ^ salt;
    h = (h ^ (h >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    h = (h ^ (h >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    h ^= h >> 31;
    (h >> 11) as f64 / 9_007_199_254_740_992.0 * 2.0 - 1.0
}

fn channel(normal: usize, variable_density: bool, spacing: [f64; 3]) -> Euler3d {
    let shape = [5, 4, 5];
    let grid = CartesianGrid::new(
        shape[0], shape[1], shape[2], spacing[0], spacing[1], spacing[2],
    )
    .unwrap();
    let gas = gas(if variable_density { 0.2 } else { 0.0 });
    let states = (0..shape.into_iter().product::<usize>())
        .map(|cell| {
            let p = grid.coordinates(cell);
            let eta = (p[normal] as f64 + 0.5) / shape[normal] as f64;
            let envelope = 4.0 * eta * (1.0 - eta);
            let rho = if variable_density {
                1.0 + 0.15 * noise(cell, 911)
            } else {
                1.0
            };
            let u = core::array::from_fn(|axis| {
                boost(normal)[axis] + 0.02 * envelope * noise(cell, 31 * axis as u64 + 7)
            });
            let temperature = 3.0 - eta + 0.05 * envelope * noise(cell, 101);
            State::from_primitive(&gas, rho, u, rho * temperature).unwrap()
        })
        .collect();
    Euler3d::new(gas, grid, states, boundaries(normal), || false).unwrap()
}

fn core_fit(
    op: &MolecularTransport3d<ConstantTransport>,
    normal: usize,
) -> (ChannelDynamicSgs, ChannelDynamicScalar) {
    let p = &op.work.primitive;
    let rho: Vec<_> = p.iter().map(|p| p.density).collect();
    let u: Vec<_> = p.iter().map(|p| p.velocity).collect();
    let t: Vec<_> = p.iter().map(|p| p.temperature).collect();
    let mut momentum =
        ChannelDynamicSgs::new_cartesian(op.grid.shape, op.grid.spacing, normal, usize::MAX)
            .unwrap();
    let mut heat = ChannelDynamicScalar::new(op.grid.shape, normal, usize::MAX).unwrap();
    let (uw, tw) = channel_wall_values(normal, &op.walls).unwrap();
    momentum
        .prepare_with_scalar(
            settings().momentum,
            settings().heat,
            &rho,
            &u,
            uw,
            &t,
            tw,
            &mut heat,
            || false,
        )
        .unwrap();
    (momentum, heat)
}

fn configured(
    flow: &Euler3d,
    normal: usize,
    wall: [MolecularWall3d; 2],
) -> MolecularTransport3d<ConstantTransport> {
    let mut result =
        MolecularTransport3d::new(flow, law(flow.gas(), 0.002, 0.003, 0.004), || false).unwrap();
    for (side, wall) in wall.into_iter().enumerate() {
        result = result
            .with_wall(BoundarySide::ALL[2 * normal + side], wall)
            .unwrap();
    }
    result
}

fn assert_faces(
    flow: &Euler3d,
    normal: usize,
    active: &mut MolecularTransport3d<ConstantTransport>,
    base: &mut MolecularTransport3d<ConstantTransport>,
) -> (ChannelDynamicSgs, ChannelDynamicScalar) {
    let base_rate = prepare(base, flow);
    let active_rate = prepare(active, flow);
    assert!(active_rate >= base_rate);
    let (momentum, heat) = core_fit(active, normal);
    let axis = (normal + 1) % 3;
    for cell in 0..flow.states.len() {
        let mut c = flow.grid.coordinates(cell);
        let face = flow.grid.cell_face(axis, c);
        c[axis] = (c[axis] + flow.grid.shape[axis] - 1) % flow.grid.shape[axis];
        let left = flow.grid.index(c);
        let pl = active.work.primitive[left];
        let pr = active.work.primitive[cell];
        let mu = 0.5 * pl.density * momentum.viscosities()[left]
            + 0.5 * pr.density * momentum.viscosities()[cell];
        let rho_alpha = 0.5 * pl.density * heat.diffusivities()[left]
            + 0.5 * pr.density * heat.diffusivities()[cell];
        let ts = tangents(axis);
        let derivative = (pr.velocity[axis] - pl.velocity[axis]) / flow.grid.spacing[axis];
        let transverse: f64 = ts
            .into_iter()
            .map(|t| {
                0.5 * active.work.gradient[left][t][t] + 0.5 * active.work.gradient[cell][t][t]
            })
            .sum();
        let normal_stress = (4.0 / 3.0) * mu * derivative - (2.0 / 3.0) * mu * transverse;
        let shear: [f64; 2] = ts.map(|t| {
            mu * ((pr.velocity[t] - pl.velocity[t]) / flow.grid.spacing[axis]
                + 0.5 * active.work.gradient[left][axis][t]
                + 0.5 * active.work.gradient[cell][axis][t])
        });
        let cp = flow
            .gas
            .source()
            .evaluate_properties(
                Temperature::new(0.5 * pl.temperature + 0.5 * pr.temperature),
                Pressure::new(0.5 * pl.pressure + 0.5 * pr.pressure),
            )
            .unwrap()
            .mass_specific_cp()
            .value();
        let energy = -normal_stress * (0.5 * pl.velocity[axis] + 0.5 * pr.velocity[axis])
            - shear[0] * (0.5 * pl.velocity[ts[0]] + 0.5 * pr.velocity[ts[0]])
            - shear[1] * (0.5 * pl.velocity[ts[1]] + 0.5 * pr.velocity[ts[1]])
            - rho_alpha * cp * (pr.temperature - pl.temperature) / flow.grid.spacing[axis];
        let a = active.work.flux[axis][face];
        let b = base.work.flux[axis][face];
        close(a.normal_momentum - b.normal_momentum, -normal_stress, 2e-13);
        for t in 0..2 {
            close(
                a.transverse_momentum[t] - b.transverse_momentum[t],
                -shear[t],
                2e-13,
            );
        }
        close(a.energy - b.energy, energy, 2e-13);
        if flow.grid.coordinates(cell)[normal] == 0 {
            let mut wall = flow.grid.coordinates(cell);
            for side in [0, flow.grid.shape[normal]] {
                wall[normal] = side;
                let face = flow.grid.cell_face(normal, wall);
                assert_eq!(active.work.flux[normal][face], base.work.flux[normal][face]);
            }
        }
    }
    (momentum, heat)
}

#[test]
fn g1_channel_dynamic_fv_uses_current_walls_geometry_source_heat_and_molecular_wall_faces() {
    for normal in 0..3 {
        for adiabatic in [false, true] {
            let flow = channel(normal, true, [0.17, 0.11, 0.13]);
            let wall = walls(normal, adiabatic);
            let mut base = configured(&flow, normal, wall);
            let mut active = MolecularTransport3d::with_channel_dynamic_thermal_sgs(
                &flow,
                law(flow.gas(), 0.002, 0.003, 0.004),
                settings(),
                wall,
                256 * flow.states.len(),
                || false,
            )
            .unwrap();
            assert_eq!(
                active.channel_dynamic_thermal_subgrid_config(),
                Some(settings())
            );
            assert_eq!(active.dynamic_thermal_subgrid_config(), None);
            assert_eq!(active.dynamic_subgrid_config(), None);
            let _ = assert_faces(&flow, normal, &mut active, &mut base);
            let changed = MolecularWall3d::NoSlip {
                velocity: core::array::from_fn(|axis| if axis == normal { 0.0 } else { 0.07 }),
                thermal: WallThermalCondition::Isothermal(Temperature::new(2.8)),
            };
            active = active
                .with_wall(BoundarySide::ALL[2 * normal + 1], changed)
                .unwrap();
            base = base
                .with_wall(BoundarySide::ALL[2 * normal + 1], changed)
                .unwrap();
            let _ = assert_faces(&flow, normal, &mut active, &mut base);
        }
    }
}

#[test]
fn g1_channel_dynamic_fv_heat_survives_zero_momentum_and_reaches_conservative_update() {
    let mut flow = channel(2, false, [1.0; 3]);
    let mut base_flow = channel(2, false, [1.0; 3]);
    let wall = walls(2, false);
    let mut base = configured(&base_flow, 2, wall);
    let mut active = MolecularTransport3d::with_channel_dynamic_thermal_sgs(
        &flow,
        law(flow.gas(), 0.002, 0.003, 0.004),
        settings(),
        wall,
        256 * flow.states.len(),
        || false,
    )
    .unwrap();
    let (momentum, heat) = assert_faces(&flow, 2, &mut active, &mut base);
    let mut heat_only = 0;
    for cell in 0..flow.states.len() {
        if momentum.viscosities()[cell] == 0.0 && heat.diffusivities()[cell] > 0.0 {
            let face = flow.grid.cell_face(0, flow.grid.coordinates(cell));
            let a = active.work.flux[0][face];
            let b = base.work.flux[0][face];
            assert_eq!(a.normal_momentum, b.normal_momentum);
            assert_eq!(a.transverse_momentum, b.transverse_momentum);
            assert_ne!(a.energy.to_bits(), b.energy.to_bits());
            heat_only += 1;
        }
    }
    assert!(heat_only > 0);
    let original = flow.totals();
    let report = active
        .advance(
            &mut flow,
            Time::new(0.001),
            MolecularControl::default(),
            || false,
        )
        .unwrap();
    base.advance(
        &mut base_flow,
        Time::new(0.001),
        MolecularControl::default(),
        || false,
    )
    .unwrap();
    assert_ne!(flow.states(), base_flow.states());
    close(flow.totals().mass, original.mass, 3e-13);
    for value in report.conservation_residual.array() {
        close(value, 0.0, 3e-12);
    }
}

fn coupled() -> MolecularEuler3d<ConstantTransport> {
    let flow = channel(2, false, [0.2; 3]);
    let source = law(flow.gas(), 0.002, 0.003, 0.004);
    MolecularEuler3d::with_channel_dynamic_thermal_sgs(
        flow,
        source,
        settings(),
        walls(2, false),
        256 * 100,
        || false,
    )
    .unwrap()
}

#[test]
fn g0_g4_g5_channel_dynamic_coupled_zero_fit_energy_and_cancelled_trials_preserve_native_owner() {
    let uniform = || {
        field(
            gas(0.0),
            CartesianGrid::new(4, 4, 4, 0.2, 0.2, 0.2).unwrap(),
            boundaries(2),
            |_| [1.0, 0.01, 0.02, 0.0, 2.0],
        )
    };
    let wall = [MolecularWall3d::NoSlip {
        velocity: boost(2),
        thermal: WallThermalCondition::Adiabatic,
    }; 2];
    let flow = uniform();
    let source = law(flow.gas(), 0.002, 0.003, 0.004);
    let mut base = MolecularEuler3d::new(flow, source, || false).unwrap();
    for side in [BoundarySide::ZMin, BoundarySide::ZMax] {
        base = base.with_wall(side, wall[0]).unwrap();
    }
    let flow = uniform();
    let source = law(flow.gas(), 0.002, 0.003, 0.004);
    let mut active = MolecularEuler3d::with_channel_dynamic_thermal_sgs(
        flow,
        source,
        settings(),
        wall,
        256 * 64,
        || false,
    )
    .unwrap();
    assert_eq!(
        active.channel_dynamic_thermal_subgrid_config(),
        Some(settings())
    );
    for _ in 0..3 {
        let a = active
            .step(Time::new(0.001), MolecularStepControl::default(), || false)
            .unwrap();
        let b = base
            .step(Time::new(0.001), MolecularStepControl::default(), || false)
            .unwrap();
        assert_eq!(a, b);
        assert_eq!(active.flow().states(), base.flow().states());
    }
    let mut reference = coupled();
    reference
        .step(Time::new(0.001), MolecularStepControl::default(), || false)
        .unwrap();
    let mut calls = 0;
    let expected = reference
        .step(Time::new(0.001), MolecularStepControl::default(), || {
            calls += 1;
            false
        })
        .unwrap();
    for value in expected.conservation_residual.array() {
        close(value, 0.0, 3e-12);
    }
    assert!(expected.boundary_change.energy.abs() > 0.0);
    for cutoff in [1, calls / 2, calls] {
        let mut candidate = coupled();
        candidate
            .step(Time::new(0.001), MolecularStepControl::default(), || false)
            .unwrap();
        let states = candidate.flow().states().to_vec();
        let time = candidate.flow().time();
        let previous = candidate.last_step().copied();
        let impulses = BoundarySide::ALL.map(|side| candidate.boundary_impulses(side).to_vec());
        let mut seen = 0;
        assert_eq!(
            candidate.step(Time::new(0.001), MolecularStepControl::default(), || {
                seen += 1;
                seen == cutoff
            }),
            Err(MolecularError::Cancelled)
        );
        assert_eq!(candidate.flow().states(), states);
        assert_eq!(candidate.flow().time(), time);
        assert_eq!(candidate.last_step().copied(), previous);
        for (i, side) in BoundarySide::ALL.into_iter().enumerate() {
            assert_eq!(candidate.boundary_impulses(side), impulses[i]);
        }
        assert_eq!(
            candidate
                .step(Time::new(0.001), MolecularStepControl::default(), || false)
                .unwrap(),
            expected
        );
        assert_eq!(candidate.flow().states(), reference.flow().states());
    }
}

#[test]
fn g0_channel_dynamic_fv_admission_bounds_and_wall_rebinding_are_explicit() {
    let mut flow = channel(2, false, [1.0; 3]);
    let bytes = MolecularTransport3d::<ConstantTransport>::required_channel_dynamic_thermal_sgs_workspace_bytes(&flow).unwrap();
    assert_eq!(bytes, 256 * flow.states.len());
    let build = |flow: &Euler3d, config, wall, budget| {
        MolecularTransport3d::with_channel_dynamic_thermal_sgs(
            flow,
            law(flow.gas(), 0.002, 0.003, 0.004),
            config,
            wall,
            budget,
            || false,
        )
    };
    assert!(matches!(
        build(&flow, settings(), walls(2, false), bytes - 1),
        Err(MolecularError::DynamicSubgrid(
            DynamicSgsError::WorkspaceBytes { .. }
        ))
    ));
    assert!(
        build(
            &flow,
            settings(),
            [MolecularWall3d::SlipAdiabatic; 2],
            bytes
        )
        .is_err()
    );
    assert!(
        build(
            &cube(4, |_| [1.0, 0.0, 0.0, 0.0, 2.0]),
            settings(),
            walls(2, false),
            usize::MAX
        )
        .is_err()
    );
    let short = field(
        gas(0.0),
        CartesianGrid::new(3, 3, 2, 1.0, 1.0, 1.0).unwrap(),
        boundaries(2),
        |_| [1.0, 0.0, 0.0, 0.0, 2.0],
    );
    assert!(build(&short, settings(), walls(2, false), usize::MAX).is_err());
    let op = build(&flow, settings(), walls(2, false), bytes).unwrap();
    assert!(
        op.with_wall(BoundarySide::ZMin, MolecularWall3d::SlipAdiabatic)
            .is_err()
    );
    let mut forbidden = settings();
    forbidden.heat.usage = SgsUse::AcceptanceRender;
    assert!(build(&flow, forbidden, walls(2, false), bytes).is_err());
    forbidden = settings();
    forbidden.heat.maximum_diffusivity_to_molecular = 0.0;
    let mut op = build(&flow, forbidden, walls(2, false), bytes).unwrap();
    let before = flow.states().to_vec();
    assert!(matches!(
        op.advance(
            &mut flow,
            Time::new(0.001),
            MolecularControl::default(),
            || false
        ),
        Err(MolecularError::SubgridThermalDiffusivityRatio { .. })
    ));
    assert_eq!(flow.states(), before);
    assert_eq!(flow.time().value(), 0.0);
    for side in BoundarySide::ALL {
        assert!(
            op.boundary_impulses(side)
                .iter()
                .all(|v| *v == Totals::default())
        );
    }
}
