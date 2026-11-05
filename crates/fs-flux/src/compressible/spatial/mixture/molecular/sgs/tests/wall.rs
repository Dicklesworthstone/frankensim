use super::*;
use fs_lbm::conjugate::turbulence::law_of_the_wall_ratios;

fn wall_config() -> EquilibriumWallConfig {
    EquilibriumWallConfig {
        maximum_wall_distance_plus: 300.0,
        maximum_mach_number: 0.3,
    }
}

fn wall(velocity: [f64; 3]) -> MolecularWall3d {
    MolecularWall3d::NoSlip {
        velocity,
        thermal: WallThermalCondition::Adiabatic,
    }
}

fn channel(normal: usize, velocity: [f64; 3]) -> MixtureEuler3d {
    let mut axes = [AxisBoundaries::Periodic; 3];
    axes[normal] = WALL;
    field(
        material(0.0),
        CartesianGrid::new(4, 4, 4, 0.2, 0.3, 0.4).unwrap(),
        Boundaries3d {
            x: axes[0],
            y: axes[1],
            z: axes[2],
        },
        |_| ([0.4, 0.6, 0.0], velocity, 2.0),
    )
}

fn source(flow: &MixtureEuler3d) -> SpatialTransport<ConstantTransport> {
    let law = ConstantTransport::new(
        flow.gas().source().clone(),
        "equilibrium-mixture-sgs-wall",
        1,
        TransportDomain {
            minimum_temperature_k: 0.001,
            maximum_temperature_k: 50.0,
            minimum_pressure_pa: 0.001,
            maximum_pressure_pa: 100.0,
        },
        TransportCoefficients {
            dynamic_viscosity: DynViscosity::new(1e-4),
            bulk_viscosity: DynViscosity::new(0.0),
            thermal_conductivity: ThermalConductivityQuantityV1::new(0.0),
        },
        vec![0.0, 0.01, 0.012, 0.01, 0.0, 0.013, 0.012, 0.013, 0.0],
    )
    .unwrap();
    SpatialTransport::new(flow.gas().clone(), law, || false).unwrap()
}

fn momentum() -> ChannelDynamicSgsConfig {
    ChannelDynamicSgsConfig {
        maximum_coefficient_squared: 100.0,
        maximum_eddy_to_molecular: 1000.0,
        usage: SgsUse::Exploration,
    }
}

fn heat() -> ChannelDynamicScalarConfig {
    ChannelDynamicScalarConfig {
        maximum_coefficient: 100.0,
        maximum_diffusivity_to_molecular: 1000.0,
        usage: SgsUse::Exploration,
    }
}

fn prepared(
    flow: &MixtureEuler3d,
    normal: usize,
    velocity: [f64; 3],
    kind: usize,
) -> MolecularTransport3d<ConstantTransport> {
    let source = source(flow);
    match kind {
        0 => {
            let bytes =
                MolecularTransport3d::required_workspace_bytes_with_sgs(flow, &source).unwrap();
            MolecularTransport3d::with_sgs(
                flow,
                source,
                config(SgsModel::Smagorinsky { cs: 0.0 }),
                bytes,
                || false,
            )
            .unwrap()
            .with_wall(BoundarySide::ALL[2 * normal], wall(velocity))
            .unwrap()
            .with_wall(BoundarySide::ALL[2 * normal + 1], wall(velocity))
            .unwrap()
        }
        1 => {
            let bytes =
                MolecularTransport3d::required_workspace_bytes_with_dynamic_channel_thermal_sgs(
                    flow, &source,
                )
                .unwrap();
            MolecularTransport3d::with_dynamic_channel_thermal_sgs(
                flow,
                source,
                normal,
                [wall(velocity); 2],
                MixtureDynamicChannelThermalSubgridConfig {
                    momentum: momentum(),
                    heat: heat(),
                    turbulent_schmidt: 0.8,
                },
                bytes,
                || false,
            )
            .unwrap()
        }
        _ => {
            let bytes =
                MolecularTransport3d::required_workspace_bytes_with_channel_dynamic_species_sgs(
                    flow, &source,
                )
                .unwrap();
            MolecularTransport3d::with_channel_dynamic_species_sgs(
                flow,
                source,
                MixtureChannelDynamicSpeciesSubgridConfig {
                    transport: ChannelDynamicThermalSubgridConfig {
                        momentum: momentum(),
                        heat: heat(),
                    },
                    species: ChannelDynamicComponentsConfig {
                        maximum_coefficient: 100.0,
                        usage: SgsUse::Exploration,
                    },
                    maximum_species_diffusivity_m2_s: 10.0,
                },
                [wall(velocity); 2],
                bytes,
                || false,
            )
            .unwrap()
        }
    }
}

#[test]
fn equilibrium_sgs_wall_mixture_faces_preserve_species_and_native_work_for_local_and_dynamic_models()
 {
    for axis in 0..3 {
        let [a, b] = tangents(axis);
        let mut velocity = [0.0; 3];
        velocity[a] = 0.2;
        velocity[b] = -0.1;
        let wall_velocity = velocity.map(|u| 0.25 * u);
        let flow = channel(axis, velocity);
        let distance = 0.5 * flow.grid.spacing[axis];
        let ratio =
            law_of_the_wall_ratios((0.15_f64.powi(2) + 0.075_f64.powi(2)).sqrt() * distance / 1e-4)
                .1;
        let mut base = prepared(&flow, axis, wall_velocity, 0);
        let base_rate = faces(&mut base, &flow).rate;
        for kind in 0..3 {
            let mut op = prepared(&flow, axis, wall_velocity, kind);
            let bytes = op.workspace_bytes();
            for upper in [false, true] {
                let side = BoundarySide::ALL[2 * axis + usize::from(upper)];
                op = op.with_equilibrium_sgs_wall(side, wall_config()).unwrap();
                assert_eq!(op.equilibrium_sgs_wall(side), Some(wall_config()));
            }
            assert!(faces(&mut op, &flow).rate > base_rate);
            assert_eq!(op.workspace_bytes(), bytes);
            for upper in [false, true] {
                let sign = if upper { 1.0 } else { -1.0 };
                let mut face = [1; 3];
                face[axis] = if upper { 4 } else { 0 };
                let index = op.grid.cell_face(axis, face) * 7;
                let flux = &op.work.flux[axis][index..index + 7];
                assert_eq!(&flux[..3], &[0.0; 3]);
                close(flux[3 + axis], 0.0, 1e-14);
                for t in [a, b] {
                    let momentum =
                        1e-4 * ratio * sign * (velocity[t] - wall_velocity[t]) / distance;
                    close(flux[3 + t], momentum, 2e-14);
                }
                close(
                    flux[6],
                    flux[3 + a] * wall_velocity[a] + flux[3 + b] * wall_velocity[b],
                    2e-14,
                );
            }
        }
    }
}

#[test]
fn equilibrium_sgs_wall_mixture_braking_conserves_species_and_total_energy_with_wall_work() {
    for speed in [0.0, 0.1] {
        let mut base = channel(1, [0.4, 0.0, 0.0]);
        let mut active = channel(1, [0.4, 0.0, 0.0]);
        let mut base_op = prepared(&base, 1, [speed, 0.0, 0.0], 0);
        let mut op = prepared(&active, 1, [speed, 0.0, 0.0], 0);
        for side in [BoundarySide::YMin, BoundarySide::YMax] {
            op = op.with_equilibrium_sgs_wall(side, wall_config()).unwrap();
        }
        let before = active.totals();
        base_op
            .advance(
                &mut base,
                Time::new(0.02),
                DissipationControl::default(),
                || false,
            )
            .unwrap();
        let report = op
            .advance(
                &mut active,
                Time::new(0.02),
                DissipationControl::default(),
                || false,
            )
            .unwrap();
        assert!(active.totals().momentum[0] < base.totals().momentum[0]);
        close(
            report.boundary_change.energy,
            speed * report.boundary_change.momentum[0],
            1e-13,
        );
        close(
            active.totals().energy - before.energy,
            report.boundary_change.energy,
            2e-13,
        );
        for i in 0..3 {
            close(
                active.totals().species_masses[i],
                before.species_masses[i],
                1e-13,
            );
        }
        assert!(active.states().chunks_exact(7).all(|row| row[2] == 0.0));
        balance(&report);
    }
}

#[test]
fn equilibrium_sgs_wall_mixture_refusal_preserves_species_loads_and_report() {
    let mut flow = channel(1, [0.4, 0.0, 0.0]);
    let mut op = prepared(&flow, 1, [0.0; 3], 0);
    for side in [BoundarySide::YMin, BoundarySide::YMax] {
        op = op.with_equilibrium_sgs_wall(side, wall_config()).unwrap();
    }
    op.advance(
        &mut flow,
        Time::new(0.01),
        DissipationControl::default(),
        || false,
    )
    .unwrap();
    let states = flow.states().to_vec();
    let loads = op.boundary_impulses.clone();
    let report = op.last_step().cloned();
    op = op
        .with_equilibrium_sgs_wall(
            BoundarySide::YMax,
            EquilibriumWallConfig {
                maximum_mach_number: 1e-6,
                ..wall_config()
            },
        )
        .unwrap();
    assert!(
        op.advance(
            &mut flow,
            Time::new(0.01),
            DissipationControl::default(),
            || false
        )
        .is_err()
    );
    assert_eq!(flow.states(), states);
    assert_eq!(op.boundary_impulses, loads);
    assert_eq!(op.last_step(), report.as_ref());
    assert!(
        op.with_wall(
            BoundarySide::YMin,
            MolecularWall3d::NoSlip {
                velocity: [0.0; 3],
                thermal: WallThermalCondition::Isothermal(Temperature::new(2.0)),
            }
        )
        .is_err()
    );
}
