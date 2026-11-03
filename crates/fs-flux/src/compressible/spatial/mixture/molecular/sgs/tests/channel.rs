use super::*;
use fs_lbm::sgs3::dynamic::channel::{ChannelDynamicScalar, ChannelDynamicSgs, ScalarWall};

fn config() -> MixtureDynamicChannelThermalSubgridConfig {
    MixtureDynamicChannelThermalSubgridConfig {
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
        turbulent_schmidt: 0.8,
    }
}

fn boundaries(normal: usize) -> Boundaries3d {
    let mut axes = [AxisBoundaries::Periodic; 3];
    axes[normal] = WALL;
    Boundaries3d {
        x: axes[0],
        y: axes[1],
        z: axes[2],
    }
}

fn walls(normal: usize, mixed: bool) -> [MolecularWall3d; 2] {
    let mut velocity = [[0.03, -0.01, 0.02], [-0.02, 0.015, -0.01]];
    velocity[0][normal] = 0.0;
    velocity[1][normal] = 0.0;
    [
        MolecularWall3d::NoSlip {
            velocity: velocity[0],
            thermal: WallThermalCondition::Isothermal(Temperature::new(2.0)),
        },
        MolecularWall3d::NoSlip {
            velocity: velocity[1],
            thermal: if mixed {
                WallThermalCondition::Adiabatic
            } else {
                WallThermalCondition::Isothermal(Temperature::new(1.0))
            },
        },
    ]
}

fn noise(i: usize, salt: u64) -> f64 {
    let mut h = ((i as u64 + 1).wrapping_mul(0x9e37_79b9_7f4a_7c15)) ^ salt;
    h = (h ^ (h >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    h = (h ^ (h >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    h ^= h >> 31;
    (h >> 11) as f64 / 9_007_199_254_740_992.0 * 2.0 - 1.0
}

fn channel_field(normal: usize, formation: f64, variable_density: bool) -> MixtureEuler3d {
    let grid = CartesianGrid::new(5, 4, 5, 0.3, 0.7, 1.1).unwrap();
    field(material(formation), grid, boundaries(normal), |point| {
        let coordinates = std::array::from_fn(|axis| (point[axis] / grid.spacing[axis]) as usize);
        let i = grid.index(coordinates);
        let eta = point[normal] / (grid.shape[normal] as f64 * grid.spacing[normal]);
        let envelope = 4.0 * eta * (1.0 - eta);
        let rho = if variable_density {
            1.0 + 0.15 * noise(i, 911)
        } else {
            1.0
        };
        let first = 0.45 + 0.1 * noise(i, 51);
        (
            [rho * first, rho * (1.0 - first), 0.0],
            std::array::from_fn(|axis| 0.02 * envelope * noise(i, 31 * axis as u64 + 7)),
            2.0 - eta + 0.05 * envelope * noise(i, 101),
        )
    })
}

fn prepared(
    flow: &MixtureEuler3d,
    normal: usize,
    walls: [MolecularWall3d; 2],
) -> MolecularTransport3d<ConstantTransport> {
    let source = original(flow);
    let bytes = MolecularTransport3d::required_workspace_bytes_with_dynamic_channel_thermal_sgs(
        flow, &source,
    )
    .unwrap();
    let base = MolecularTransport3d::required_workspace_bytes(flow, &source).unwrap();
    assert_eq!(bytes - base, 264 * flow.cell_count());
    let result = MolecularTransport3d::with_dynamic_channel_thermal_sgs(
        flow,
        source,
        normal,
        walls,
        config(),
        bytes,
        || false,
    )
    .unwrap();
    assert_eq!(result.workspace_bytes() + flow.workspace_bytes(), bytes);
    assert_eq!(
        result.dynamic_channel_thermal_subgrid_config(),
        Some(config())
    );
    assert_eq!(result.dynamic_thermal_subgrid_config(), None);
    assert_eq!(result.dynamic_subgrid_config(), None);
    assert_eq!(result.subgrid_config(), None);
    result
}

fn molecular(
    flow: &MixtureEuler3d,
    normal: usize,
    walls: [MolecularWall3d; 2],
) -> MolecularTransport3d<ConstantTransport> {
    let mut op = operator(flow, None);
    for (end, wall) in walls.into_iter().enumerate() {
        op = op
            .with_wall(BoundarySide::ALL[2 * normal + end], wall)
            .unwrap();
    }
    op
}

#[test]
fn channel_mixture_faces_keep_source_enthalpy_cp_and_molecular_wall_heat_work() {
    let mut active_species = 0;
    let mut active_heat = 0;
    for normal in 0..3 {
        for mixed in [false, true] {
            let flow = channel_field(normal, 0.0, true);
            let walls = walls(normal, mixed);
            let mut op = prepared(&flow, normal, walls);
            let mut base = molecular(&flow, normal, walls);
            // An original molecular wall can still control the global bound.
            assert!(faces(&mut op, &flow).rate >= faces(&mut base, &flow).rate);
            let rho: Vec<_> = op.work.primitive.iter().map(|p| p.density).collect();
            let velocity: Vec<_> = op.work.primitive.iter().map(|p| p.velocity).collect();
            let temperature: Vec<_> = op.work.primitive.iter().map(|p| p.temperature).collect();
            let mut core = ChannelDynamicSgs::new_cartesian(
                flow.grid.shape,
                flow.grid.spacing,
                normal,
                usize::MAX,
            )
            .unwrap();
            let mut heat = ChannelDynamicScalar::new(flow.grid.shape, normal, usize::MAX).unwrap();
            let wall_velocity = walls.map(|w| match w {
                MolecularWall3d::NoSlip { velocity, .. } => velocity,
                _ => unreachable!(),
            });
            let scalar_walls = walls.map(|w| match w {
                MolecularWall3d::NoSlip {
                    thermal: WallThermalCondition::Isothermal(t),
                    ..
                } => ScalarWall::Value(t.value()),
                _ => ScalarWall::ZeroGradient,
            });
            core.prepare_with_scalar(
                config().momentum,
                config().heat,
                &rho,
                &velocity,
                wall_velocity,
                &temperature,
                scalar_walls,
                &mut heat,
                || false,
            )
            .unwrap();
            for axis in 0..3 {
                for left in 0..flow.cell_count() {
                    let mut c = flow.grid.coordinates(left);
                    if axis == normal && c[axis] + 1 == flow.grid.shape[axis] {
                        continue;
                    }
                    c[axis] = (c[axis] + 1) % flow.grid.shape[axis];
                    let right = flow.grid.index(c);
                    let face = flow.grid.cell_face(axis, c);
                    let delta =
                        |k| op.work.flux[axis][7 * face + k] - base.work.flux[axis][7 * face + k];
                    let mu = 0.5 * rho[left] * core.viscosities()[left]
                        + 0.5 * rho[right] * core.viscosities()[right];
                    let t = 0.5 * temperature[left] + 0.5 * temperature[right];
                    let mut cp = 0.0;
                    let mut species_energy = 0.0;
                    for (k, model) in flow.gas().source().models().iter().enumerate() {
                        let yl = flow.states()[7 * left + k] / rho[left];
                        let yr = flow.states()[7 * right + k] / rho[right];
                        let j = -(mu / config().turbulent_schmidt) * (yr - yl)
                            / flow.grid.spacing[axis];
                        close(delta(k), j, 2e-14);
                        let inverse_mass = model.molar_mass().value().recip();
                        cp += 0.5
                            * (yl + yr)
                            * inverse_mass
                            * model
                                .evaluate_properties(Temperature::new(t))
                                .unwrap()
                                .cp()
                                .value();
                        let h = inverse_mass
                            * (model
                                .molar_internal_energy_change(
                                    Temperature::new(flow.gas().domain().minimum_temperature),
                                    Temperature::new(t),
                                )
                                .unwrap()
                                .value()
                                + R * t);
                        species_energy += h * j;
                        active_species += usize::from(j.abs() > 1e-9);
                    }
                    let weight = 0.5 * rho[left] * heat.diffusivities()[left]
                        + 0.5 * rho[right] * heat.diffusivities()[right];
                    let q = -weight * cp * (temperature[right] - temperature[left])
                        / flow.grid.spacing[axis];
                    let work: f64 = (0..3)
                        .map(|a| 0.5 * (velocity[left][a] + velocity[right][a]) * delta(3 + a))
                        .sum();
                    close(delta(6), q + species_energy + work, 5e-13);
                    active_heat += usize::from(q.abs() > 1e-9);
                    assert_eq!(delta(2), 0.0);
                }
            }
            for cell in 0..flow.cell_count() {
                let mut c = flow.grid.coordinates(cell);
                if c[normal] != 0 {
                    continue;
                }
                for end in [0, flow.grid.shape[normal]] {
                    c[normal] = end;
                    let f = flow.grid.cell_face(normal, c);
                    assert_eq!(
                        &op.work.flux[normal][7 * f..7 * f + 7],
                        &base.work.flux[normal][7 * f..7 * f + 7]
                    );
                }
            }
        }
    }
    // Plane averaging may reject backscatter in every plane of one fixture.
    // Across the rotated wall cases both modeled energy mechanisms must act.
    assert!(
        active_species > 0 && active_heat > 0,
        "active species faces={active_species}, active heat faces={active_heat}"
    );
}

#[test]
fn channel_mixture_heat_remains_active_on_planes_with_zero_momentum_fit() {
    let flow = channel_field(2, 0.0, false);
    let walls = [
        MolecularWall3d::NoSlip {
            velocity: [0.0; 3],
            thermal: WallThermalCondition::Isothermal(Temperature::new(2.0)),
        },
        MolecularWall3d::NoSlip {
            velocity: [0.0; 3],
            thermal: WallThermalCondition::Isothermal(Temperature::new(1.0)),
        },
    ];
    let mut op = prepared(&flow, 2, walls);
    let mut base = molecular(&flow, 2, walls);
    faces(&mut op, &flow);
    faces(&mut base, &flow);
    let sgs = op.sgs.as_ref().unwrap();
    let ViscosityClosure::DynamicChannelThermal(thermal) = &sgs.closure else {
        unreachable!()
    };
    let mut active = 0;
    for left in 0..flow.cell_count() {
        let mut c = flow.grid.coordinates(left);
        c[0] = (c[0] + 1) % flow.grid.shape[0];
        let right = flow.grid.index(c);
        if sgs.dynamic_viscosity[left] + sgs.dynamic_viscosity[right] != 0.0
            || thermal.face_heat(left, right).0 == 0.0
        {
            continue;
        }
        let face = flow.grid.cell_face(0, c);
        assert_eq!(
            &op.work.flux[0][7 * face..7 * face + 6],
            &base.work.flux[0][7 * face..7 * face + 6]
        );
        assert_ne!(
            op.work.flux[0][7 * face + 6],
            base.work.flux[0][7 * face + 6]
        );
        assert!(sgs.heat_capacity_v_density[left] > 0.0);
        active += 1;
    }
    assert!(active > 0);
}

#[test]
fn channel_mixture_advances_with_wall_energy_balance_gauge_and_current_wall_data() {
    let mut flow = channel_field(1, 0.0, true);
    let mut shifted = channel_field(1, 1e12, true);
    let old_walls = walls(1, false);
    let current = walls(1, true);
    let mut op = prepared(&flow, 1, old_walls)
        .with_wall(BoundarySide::YMax, current[1])
        .unwrap();
    let mut expected = prepared(&shifted, 1, current);
    let bytes = op.workspace_bytes() + flow.workspace_bytes();
    let controls = DissipationControl::default();
    let report = op
        .advance(&mut flow, Time::new(0.02), controls, || false)
        .unwrap();
    let shifted_report = expected
        .advance(&mut shifted, Time::new(0.02), controls, || false)
        .unwrap();
    balance(&report);
    assert_eq!(
        report.conservation_residual,
        shifted_report.conservation_residual
    );
    assert!(report.boundary_change.energy != 0.0);
    assert!(report.boundary_change.momentum.iter().any(|v| *v != 0.0));
    assert_eq!(flow.states(), shifted.states());
    assert_eq!(op.boundary_impulses, expected.boundary_impulses);
    assert!(flow.states().chunks_exact(7).all(|row| row[2] == 0.0));
    assert_eq!(op.workspace_bytes() + flow.workspace_bytes(), bytes);
}

#[test]
fn channel_mixture_refusal_and_cancelled_refit_preserve_accepted_state_and_retry() {
    let mut flow = channel_field(2, 0.0, true);
    let walls = walls(2, true);
    let source = original(&flow);
    let bytes = MolecularTransport3d::required_workspace_bytes_with_dynamic_channel_thermal_sgs(
        &flow, &source,
    )
    .unwrap();
    assert!(
        MolecularTransport3d::with_dynamic_channel_thermal_sgs(
            &flow,
            source,
            2,
            walls,
            config(),
            bytes - 1,
            || false,
        )
        .is_err()
    );
    assert!(
        prepared(&flow, 2, walls)
            .with_wall(BoundarySide::ZMin, MolecularWall3d::SlipAdiabatic)
            .is_err()
    );
    let small = field(
        material(0.0),
        CartesianGrid::new(3, 3, 2, 0.2, 0.3, 0.4).unwrap(),
        boundaries(2),
        |_| ([0.6, 0.4, 0.0], [0.0; 3], 1.5),
    );
    let fixed = self::walls(2, false);
    assert!(
        prepared(&small, 2, fixed)
            .with_wall(BoundarySide::ZMax, walls[1])
            .is_err()
    );
    for normal in [0, 3] {
        assert!(
            MolecularTransport3d::with_dynamic_channel_thermal_sgs(
                &flow,
                original(&flow),
                normal,
                walls,
                config(),
                bytes,
                || false,
            )
            .is_err()
        );
    }
    let mut op = prepared(&flow, 2, walls);
    let controls = DissipationControl::default();
    let report = op
        .advance(&mut flow, Time::new(0.002), controls, || false)
        .unwrap();
    let state = flow.states().to_vec();
    let loads = op.boundary_impulses.clone();
    let mut polls = 0;
    assert_eq!(
        op.advance(&mut flow, Time::new(0.002), controls, || {
            polls += 1;
            polls == 1000
        }),
        Err(DissipationError::Cancelled)
    );
    assert_eq!(flow.states(), state);
    assert_eq!(op.boundary_impulses, loads);
    assert_eq!(op.last_step(), Some(&report));
    let mut clean = MixtureEuler3d::new(
        flow.gas().clone(),
        flow.grid,
        state.clone(),
        flow.boundaries,
        flow.workspace_bytes(),
        || false,
    )
    .unwrap();
    let mut fresh = prepared(&clean, 2, walls);
    let a = op
        .advance(&mut flow, Time::new(0.002), controls, || false)
        .unwrap();
    let b = fresh
        .advance(&mut clean, Time::new(0.002), controls, || false)
        .unwrap();
    assert_eq!(a, b);
    assert_eq!(flow.states(), clean.states());
    let state = flow.states().to_vec();
    let ViscosityClosure::DynamicChannelThermal(thermal) = &mut op.sgs.as_mut().unwrap().closure
    else {
        unreachable!()
    };
    thermal.config.heat.maximum_diffusivity_to_molecular = f64::MIN_POSITIVE;
    assert_eq!(
        op.advance(&mut flow, Time::new(0.002), controls, || false),
        Err(DissipationError::Invalid(
            "mixture subgrid heat/molecular conductivity ratio exceeded"
        ))
    );
    assert_eq!(flow.states(), state);
    assert_eq!(op.last_step(), Some(&a));
}
