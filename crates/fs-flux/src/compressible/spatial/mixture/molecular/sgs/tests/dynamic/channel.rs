use super::*;
use fs_gas::multicomponent::dissipation::SoretDufourTransport;
use fs_lbm::sgs3::dynamic::channel::components::ChannelDynamicComponents;
use fs_lbm::sgs3::dynamic::channel::{
    ChannelDynamicScalar, ChannelDynamicScalarConfig, ChannelDynamicSgs, ChannelDynamicSgsConfig,
};
use fs_thermochem::MassSpecificEnergyQuantityV1;

fn config() -> MixtureChannelDynamicSpeciesSubgridConfig {
    MixtureChannelDynamicSpeciesSubgridConfig {
        transport: ChannelDynamicThermalSubgridConfig {
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
        },
        species: ChannelDynamicComponentsConfig {
            maximum_coefficient: 100.0,
            usage: SgsUse::Exploration,
        },
        maximum_species_diffusivity_m2_s: 10.0,
    }
}

fn boundaries(normal: usize) -> Boundaries3d {
    let mut boundaries = PERIODIC;
    match normal {
        0 => boundaries.x = WALL,
        1 => boundaries.y = WALL,
        _ => boundaries.z = WALL,
    }
    boundaries
}

fn walls() -> [MolecularWall3d; 2] {
    [1.1, 1.0].map(|temperature| MolecularWall3d::NoSlip {
        velocity: [0.0; 3],
        thermal: WallThermalCondition::Isothermal(Temperature::new(temperature)),
    })
}

fn noise(cell: usize, salt: u64) -> f64 {
    let mut h = ((cell as u64 + 1).wrapping_mul(0x9e37_79b9_7f4a_7c15)) ^ salt;
    h = (h ^ (h >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    h = (h ^ (h >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    h ^= h >> 31;
    (h >> 11) as f64 / 9_007_199_254_740_992.0 * 2.0 - 1.0
}

fn channel(
    normal: usize,
    spacing: [f64; 3],
    variable_density: bool,
    formation: f64,
) -> MixtureEuler3d {
    let grid = CartesianGrid::new(5, 4, 5, spacing[0], spacing[1], spacing[2]).unwrap();
    field(material(formation), grid, boundaries(normal), |point| {
        let c = std::array::from_fn(|axis| (point[axis] / spacing[axis]) as usize);
        let cell = grid.index(c);
        let eta = (c[normal] as f64 + 0.5) / grid.shape[normal] as f64;
        let envelope = 4.0 * eta * (1.0 - eta);
        let phi = 1.0 - eta + 0.05 * envelope * noise(cell, 101);
        let fraction = 0.25 + 0.5 * phi;
        let rho = if variable_density {
            1.0 + 0.15 * noise(cell, 911)
        } else {
            1.0
        };
        (
            [rho * fraction, rho * (1.0 - fraction), 0.0],
            std::array::from_fn(|j| 0.02 * envelope * noise(cell, 31 * j as u64 + 7)),
            1.0 + 0.1 * phi,
        )
    })
}

fn prepared<L: TransportLaw>(
    flow: &MixtureEuler3d,
    source: SpatialTransport<L>,
    cfg: MixtureChannelDynamicSpeciesSubgridConfig,
    wall: [MolecularWall3d; 2],
) -> MolecularTransport3d<L> {
    let base = MolecularTransport3d::required_workspace_bytes(flow, &source).unwrap();
    let bytes = MolecularTransport3d::required_workspace_bytes_with_channel_dynamic_species_sgs(
        flow, &source,
    )
    .unwrap();
    assert_eq!(bytes - base, 352 * flow.cell_count());
    let op = MolecularTransport3d::with_channel_dynamic_species_sgs(
        flow,
        source,
        cfg,
        wall,
        bytes,
        || false,
    )
    .unwrap();
    assert_eq!(op.workspace_bytes() + flow.workspace_bytes(), bytes);
    assert_eq!(op.channel_dynamic_species_subgrid_config(), Some(cfg));
    assert_eq!(op.dynamic_species_subgrid_config(), None);
    assert_eq!(op.dynamic_thermal_subgrid_config(), None);
    assert_eq!(op.dynamic_subgrid_config(), None);
    assert_eq!(op.subgrid_config(), None);
    op
}

fn baseline<L: TransportLaw>(
    flow: &MixtureEuler3d,
    source: SpatialTransport<L>,
    normal: usize,
    wall: [MolecularWall3d; 2],
) -> MolecularTransport3d<L> {
    let bytes = MolecularTransport3d::required_workspace_bytes(flow, &source).unwrap();
    let mut op = MolecularTransport3d::new(flow, source, bytes, || false).unwrap();
    for (side, value) in wall.into_iter().enumerate() {
        op = op
            .with_wall(BoundarySide::ALL[2 * normal + side], value)
            .unwrap();
    }
    op
}

fn stage<L: TransportLaw>(op: &mut MolecularTransport3d<L>, flow: &MixtureEuler3d) -> Stage {
    crate::compressible::spatial::mixture::molecular::spatial(
        &mut op.source,
        &op.grid,
        op.boundaries,
        &op.walls,
        flow.states(),
        &mut op.work.primitive,
        &mut op.work.gradient,
        &mut op.work.flux,
        op.sgs.as_mut(),
        DissipationControl::default(),
        &mut || false,
    )
    .unwrap()
}

fn check_shared(
    op: &MolecularTransport3d<ConstantTransport>,
    flow: &MixtureEuler3d,
    normal: usize,
) {
    let (wall_velocity, heat_walls) =
        crate::compressible::spatial::molecular::sgs::channel_wall_values(normal, &op.walls)
            .unwrap();
    let rho: Vec<_> = op.work.primitive.iter().map(|p| p.density).collect();
    let velocity: Vec<_> = op.work.primitive.iter().map(|p| p.velocity).collect();
    let temperature: Vec<_> = op.work.primitive.iter().map(|p| p.temperature).collect();
    let mut momentum =
        ChannelDynamicSgs::new_cartesian(flow.grid.shape, flow.grid.spacing, normal, usize::MAX)
            .unwrap();
    let mut heat = ChannelDynamicScalar::new(flow.grid.shape, normal, usize::MAX).unwrap();
    let mut species =
        ChannelDynamicComponents::new(flow.grid.shape, normal, 3, usize::MAX).unwrap();
    let cfg = config();
    momentum
        .prepare_with_scalar_and_components(
            cfg.transport.momentum,
            cfg.transport.heat,
            cfg.species,
            &rho,
            &velocity,
            wall_velocity,
            &temperature,
            heat_walls,
            &mut heat,
            |species, cell| flow.states()[7 * cell + species] / rho[cell],
            &mut species,
            || false,
        )
        .unwrap();
    let sgs = op.sgs.as_ref().unwrap();
    let ViscosityClosure::ChannelDynamicSpecies(active) = &sgs.closure else {
        unreachable!()
    };
    for cell in 0..flow.cell_count() {
        assert_eq!(
            sgs.dynamic_viscosity[cell],
            rho[cell] * momentum.viscosities()[cell]
        );
        assert_eq!(
            active.face_heat(cell, cell).0,
            rho[cell] * heat.diffusivities()[cell]
        );
        assert_eq!(
            active.face_species(cell, cell),
            rho[cell] * species.diffusivities()[cell]
        );
    }
}

fn check_wall_flux<A: TransportLaw, B: TransportLaw>(
    active: &MolecularTransport3d<A>,
    base: &MolecularTransport3d<B>,
    flow: &MixtureEuler3d,
    normal: usize,
) {
    for cell in 0..flow.cell_count() {
        let c = flow.grid.coordinates(cell);
        if c[normal] != 0 {
            continue;
        }
        for side in [0, flow.grid.shape[normal]] {
            let mut face_c = c;
            face_c[normal] = side;
            let face = flow.grid.cell_face(normal, face_c) * 7;
            assert_eq!(
                &active.work.flux[normal][face..face + 7],
                &base.work.flux[normal][face..face + 7]
            );
            assert!(
                active.work.flux[normal][face..face + 3]
                    .iter()
                    .all(|&j| j == 0.0)
            );
        }
    }
}

#[test]
fn g1_channel_mixture_uses_current_walls_and_rectangular_favre_species_on_all_axes() {
    for normal in 0..3 {
        let flow = channel(normal, [0.09, 0.14, 0.20], true, 0.0);
        let mut wall = walls();
        let mut velocity = [0.0; 3];
        velocity[(normal + 1) % 3] = 0.004;
        wall[1] = MolecularWall3d::NoSlip {
            velocity,
            thermal: WallThermalCondition::Adiabatic,
        };
        let mut active = prepared(&flow, original(&flow), config(), wall);
        let mut base = baseline(&flow, original(&flow), normal, wall);
        stage(&mut active, &flow);
        stage(&mut base, &flow);
        check_shared(&active, &flow, normal);
        check_wall_flux(&active, &base, &flow, normal);

        velocity[(normal + 1) % 3] = -0.015;
        let replacement = MolecularWall3d::NoSlip {
            velocity,
            thermal: WallThermalCondition::Adiabatic,
        };
        active = active
            .with_wall(BoundarySide::ALL[2 * normal], replacement)
            .unwrap();
        base = base
            .with_wall(BoundarySide::ALL[2 * normal], replacement)
            .unwrap();
        stage(&mut active, &flow);
        stage(&mut base, &flow);
        check_shared(&active, &flow, normal);
        check_wall_flux(&active, &base, &flow, normal);
    }
}

#[test]
fn g1_channel_species_and_referenced_enthalpy_remain_active_at_zero_eddy_viscosity() {
    let flow = channel(2, [1.0; 3], false, 0.0);
    let mut active = prepared(&flow, original(&flow), config(), walls());
    let mut base = baseline(&flow, original(&flow), 2, walls());
    assert!(stage(&mut active, &flow).rate > stage(&mut base, &flow).rate);
    let sgs = active.sgs.as_ref().unwrap();
    let ViscosityClosure::ChannelDynamicSpecies(model) = &sgs.closure else {
        unreachable!()
    };
    let mut nonzero = 0;
    for axis in 0..2 {
        for left in 0..flow.cell_count() {
            let mut c = flow.grid.coordinates(left);
            c[axis] = (c[axis] + 1) % flow.grid.shape[axis];
            let right = flow.grid.index(c);
            if sgs.dynamic_viscosity[left] != 0.0 || sgs.dynamic_viscosity[right] != 0.0 {
                continue;
            }
            let pl = active.work.primitive[left];
            let pr = active.work.primitive[right];
            let j = -model.face_species(left, right)
                * (flow.states()[7 * right] / pr.density - flow.states()[7 * left] / pl.density)
                / flow.grid.spacing[axis];
            let face = flow.grid.cell_face(axis, c) * 7;
            let delta: [f64; 7] = std::array::from_fn(|k| {
                active.work.flux[axis][face + k] - base.work.flux[axis][face + k]
            });
            close(delta[0], j, 4e-14);
            close(delta[1], -j, 4e-14);
            assert_eq!(delta[2], 0.0);
            assert_eq!(&delta[3..6], &[0.0; 3]);
            let t = 0.5 * pl.temperature + 0.5 * pr.temperature;
            let y = 0.5 * flow.states()[7 * left] / pl.density
                + 0.5 * flow.states()[7 * right] / pr.density;
            let cp = y * (3.5 + 0.2 * t) + (1.0 - y) * 0.5 * (4.5 + 0.3 * t);
            let h: Vec<_> = flow
                .gas()
                .source()
                .models()
                .iter()
                .take(2)
                .map(|model| {
                    (model
                        .molar_internal_energy_change(
                            Temperature::new(flow.gas().domain().minimum_temperature),
                            Temperature::new(t),
                        )
                        .unwrap()
                        .value()
                        + R * t)
                        / model.molar_mass().value()
                })
                .collect();
            let heat = -model.face_heat(left, right).0 * cp * (pr.temperature - pl.temperature)
                / flow.grid.spacing[axis];
            close(delta[6], heat + (h[0] - h[1]) * j, 4e-14);
            nonzero += usize::from(j.abs() > 1e-10);
        }
    }
    assert!(
        nonzero > 0,
        "fixture must expose species transport on a plane with clipped momentum fit"
    );
}

fn cross_source(
    flow: &MixtureEuler3d,
) -> SpatialTransport<SoretDufourTransport<ConstantTransport>> {
    let heats = flow
        .gas()
        .source()
        .models()
        .iter()
        .enumerate()
        .map(|(k, model)| {
            (
                model.species().clone(),
                MassSpecificEnergyQuantityV1::new([0.2, -0.3, 0.0][k]),
            )
        })
        .collect();
    let law = SoretDufourTransport::new(
        original(flow).law().clone(),
        "channel-cross-source",
        1,
        TransportDomain {
            minimum_temperature_k: 0.001,
            maximum_temperature_k: 50.0,
            minimum_pressure_pa: 0.001,
            maximum_pressure_pa: 100.0,
        },
        heats,
    )
    .unwrap();
    SpatialTransport::new(flow.gas().clone(), law, || false).unwrap()
}

#[test]
fn g3_channel_mixture_preserves_cross_diffusion_wall_flux_and_conservative_source_gauge() {
    let mut flow = channel(2, [0.09, 0.14, 0.20], true, 0.0);
    let mut shifted = channel(2, [0.09, 0.14, 0.20], true, 1e12);
    let mut active = prepared(&flow, cross_source(&flow), config(), walls());
    let mut shifted_op = prepared(&shifted, cross_source(&shifted), config(), walls());
    let mut molecular = baseline(&flow, cross_source(&flow), 2, walls());
    let mut plain = prepared(&flow, original(&flow), config(), walls());
    let mut plain_molecular = baseline(&flow, original(&flow), 2, walls());
    stage(&mut active, &flow);
    stage(&mut molecular, &flow);
    stage(&mut plain, &flow);
    stage(&mut plain_molecular, &flow);
    check_wall_flux(&active, &molecular, &flow, 2);
    let mut cross_effect = false;
    for axis in 0..3 {
        for face in 0..active.work.flux[axis].len() / 7 {
            for component in [0, 1, 6] {
                let i = 7 * face + component;
                let cross = molecular.work.flux[axis][i] - plain_molecular.work.flux[axis][i];
                close(
                    active.work.flux[axis][i] - plain.work.flux[axis][i],
                    cross,
                    3e-14,
                );
                cross_effect |= cross.abs() > 1e-8;
            }
        }
    }
    assert!(cross_effect, "source thermodiffusion must be active");
    let bytes = active.workspace_bytes() + flow.workspace_bytes();
    let controls = DissipationControl::default();
    let report = active
        .advance(&mut flow, Time::new(0.001), controls, || false)
        .unwrap();
    balance(&report);
    shifted_op
        .advance(&mut shifted, Time::new(0.001), controls, || false)
        .unwrap();
    assert_eq!(flow.states(), shifted.states());
    assert!(flow.states().chunks_exact(7).all(|row| row[2] == 0.0));
    assert_eq!(active.workspace_bytes() + flow.workspace_bytes(), bytes);
}

#[test]
fn g0_channel_zero_fit_clears_all_transport_and_keeps_molecular_trajectory() {
    let flow = channel(2, [1.0; 3], false, 0.0);
    let mut active = prepared(&flow, original(&flow), config(), walls());
    stage(&mut active, &flow);
    let adiabatic = MolecularWall3d::NoSlip {
        velocity: [0.0; 3],
        thermal: WallThermalCondition::Adiabatic,
    };
    for side in [4, 5] {
        active = active
            .with_wall(BoundarySide::ALL[side], adiabatic)
            .unwrap();
    }
    let uniform = || {
        field(material(0.0), flow.grid, boundaries(2), |_| {
            ([0.5, 0.5, 0.0], [0.0; 3], 1.0)
        })
    };
    let mut actual = uniform();
    let mut expected = uniform();
    let mut base = baseline(&expected, original(&expected), 2, [adiabatic; 2]);
    stage(&mut active, &actual);
    stage(&mut base, &expected);
    assert_eq!(active.work.flux, base.work.flux);
    let sgs = active.sgs.as_ref().unwrap();
    assert!(sgs.dynamic_viscosity.iter().all(|&mu| mu == 0.0));
    let ViscosityClosure::ChannelDynamicSpecies(model) = &sgs.closure else {
        unreachable!()
    };
    assert!(
        (0..actual.cell_count())
            .all(|cell| model.face_heat(cell, cell).0 == 0.0
                && model.face_species(cell, cell) == 0.0)
    );
    let controls = DissipationControl::default();
    let report = active
        .advance(&mut actual, Time::new(0.001), controls, || false)
        .unwrap();
    let reference = base
        .advance(&mut expected, Time::new(0.001), controls, || false)
        .unwrap();
    assert_eq!(report, reference);
    assert_eq!(actual.states(), expected.states());
}

#[test]
fn g4_channel_species_admission_late_cancellation_and_retry_keep_accepted_state() {
    let mut flow = channel(2, [1.0; 3], false, 0.0);
    let source = original(&flow);
    let bytes = MolecularTransport3d::required_workspace_bytes_with_channel_dynamic_species_sgs(
        &flow, &source,
    )
    .unwrap();
    assert!(
        MolecularTransport3d::with_channel_dynamic_species_sgs(
            &flow,
            source,
            config(),
            walls(),
            bytes - 1,
            || false
        )
        .is_err()
    );
    for (shape, boundary) in [
        ([3, 3, 3], PERIODIC),
        ([3, 3, 3], WALLS),
        ([3, 3, 2], boundaries(2)),
    ] {
        let grid = CartesianGrid::new(shape[0], shape[1], shape[2], 1.0, 1.0, 1.0).unwrap();
        let bad = field(material(0.0), grid, boundary, |_| {
            ([0.5, 0.5, 0.0], [0.0; 3], 1.0)
        });
        assert!(
            MolecularTransport3d::with_channel_dynamic_species_sgs(
                &bad,
                original(&bad),
                config(),
                walls(),
                usize::MAX,
                || false
            )
            .is_err()
        );
    }
    assert!(
        MolecularTransport3d::with_channel_dynamic_species_sgs(
            &flow,
            original(&flow),
            config(),
            [MolecularWall3d::SlipAdiabatic; 2],
            bytes,
            || false
        )
        .is_err()
    );
    let rebound = prepared(&flow, original(&flow), config(), walls());
    assert!(
        rebound
            .with_wall(BoundarySide::ALL[4], MolecularWall3d::SlipAdiabatic)
            .is_err()
    );
    for bound in [-1.0, f64::INFINITY, f64::NAN] {
        let mut invalid = config();
        invalid.maximum_species_diffusivity_m2_s = bound;
        assert!(invalid.validate().is_err());
    }
    let mut forbidden = config();
    forbidden.species.usage = SgsUse::AcceptanceRender;
    assert!(forbidden.validate().is_err());

    let mut cfg = config();
    cfg.maximum_species_diffusivity_m2_s = 0.0;
    let mut active = prepared(&flow, original(&flow), cfg, walls());
    let before = flow.states().to_vec();
    let controls = DissipationControl::default();
    assert!(
        active
            .advance(&mut flow, Time::new(0.001), controls, || false)
            .is_err()
    );
    assert_eq!(flow.states(), before);
    assert_eq!(active.last_step(), None);
    let ViscosityClosure::ChannelDynamicSpecies(species) =
        &mut active.sgs.as_mut().unwrap().closure
    else {
        unreachable!()
    };
    species.config = config();
    let mut reference_flow = channel(2, [1.0; 3], false, 0.0);
    let mut reference = prepared(
        &reference_flow,
        original(&reference_flow),
        config(),
        walls(),
    );
    let mut polls = 0;
    let report = reference
        .advance(&mut reference_flow, Time::new(0.001), controls, || {
            polls += 1;
            false
        })
        .unwrap();
    let loads = active.boundary_impulses.clone();
    let mut cancelled_polls = 0;
    assert_eq!(
        active.advance(&mut flow, Time::new(0.001), controls, || {
            cancelled_polls += 1;
            cancelled_polls == polls - 1
        }),
        Err(DissipationError::Cancelled)
    );
    assert_eq!(flow.states(), before);
    assert_eq!(active.boundary_impulses, loads);
    assert_eq!(active.last_step(), None);
    let retry = active
        .advance(&mut flow, Time::new(0.001), controls, || false)
        .unwrap();
    assert_eq!(retry, report);
    assert_eq!(flow.states(), reference_flow.states());
}
