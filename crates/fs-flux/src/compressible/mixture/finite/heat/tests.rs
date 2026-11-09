//! Physical and transactional regressions on both original Cartesian owners.
use super::*;
use crate::compressible::mixture::{
    self as m2,
    tests::{close, standard},
};
use crate::compressible::spatial::{self as spatial, mixture as m3};
use crate::compressible::{
    self as planar, AxisBoundaries, BoundaryCondition, Reconstruction, StepControl,
};

enum Flow {
    Planar(m2::MixtureEuler2d, planar::BoundarySide),
    Spatial(m3::MixtureEuler3d, spatial::BoundarySide),
}

#[derive(Debug, PartialEq)]
enum Snapshot {
    Planar(
        Vec<f64>,
        Plenum,
        Option<WallThermal>,
        Option<m2::FiniteReservoirReport>,
        Time,
        Option<m2::StepReport>,
        usize,
    ),
    Spatial(
        Vec<f64>,
        Plenum,
        Option<WallThermal>,
        Option<m3::FiniteReservoirReport>,
        Time,
        Option<m3::StepReport>,
        usize,
    ),
}

impl Flow {
    fn new(dimensions: usize, slot: usize, thermal: Option<WallThermal>) -> Self {
        let gas = standard(0.2);
        let open = AxisBoundaries::Ends {
            lower: if slot % 2 == 0 {
                BoundaryCondition::Transmissive
            } else {
                BoundaryCondition::Reflecting
            },
            upper: if slot % 2 == 1 {
                BoundaryCondition::Transmissive
            } else {
                BoundaryCondition::Reflecting
            },
        };
        let mut axes = [AxisBoundaries::Periodic; 3];
        axes[slot / 2] = open;
        if dimensions == 2 {
            let side = planar::SIDES[slot];
            let grid = planar::CartesianGrid::new(2, 2, 0.2, 0.3, 0.4).unwrap();
            let mut row = [0.0; 6];
            m2::state_into(&gas, &[0.8, 0.2, 0.0], [0.0; 2], 1.0, &mut row).unwrap();
            let states = row.repeat(4);
            let bytes =
                m2::MixtureEuler2d::required_workspace_bytes_with_reservoirs(&gas, &grid, 1)
                    .unwrap();
            let flow = m2::MixtureEuler2d::new_with_reservoirs(
                gas,
                grid,
                states,
                planar::Boundaries2d {
                    x: axes[0],
                    y: axes[1],
                },
                &[m2::ReservoirBoundary {
                    side,
                    partial_densities: &[0.2, 0.8, 0.0],
                    velocity: [0.0; 2],
                    pressure_pa: 1.0,
                }],
                bytes,
                || false,
            )
            .unwrap();
            let bytes = flow
                .required_workspace_bytes_with_finite_reservoirs(1)
                .unwrap();
            let flow = flow
                .with_finite_reservoirs(
                    &[m2::FiniteReservoir {
                        side,
                        volume_m3: 0.06,
                    }],
                    bytes,
                    || false,
                )
                .unwrap();
            let flow = if let Some(thermal) = thermal {
                let bytes = flow
                    .required_workspace_bytes_with_finite_reservoir_heat(1)
                    .unwrap();
                let flow = flow
                    .with_finite_reservoir_heat(
                        &[m2::FiniteReservoirHeat { side, thermal }],
                        WallControl::default(),
                        bytes,
                        || false,
                    )
                    .unwrap();
                assert_eq!(flow.workspace_bytes(), bytes);
                flow
            } else {
                flow
            };
            Self::Planar(flow, side)
        } else {
            let side = spatial::BoundarySide::ALL[slot];
            let grid = spatial::CartesianGrid::new(2, 2, 2, 0.2, 0.3, 0.4).unwrap();
            let mut row = [0.0; 7];
            m3::state_into(&gas, &[0.8, 0.2, 0.0], [0.0; 3], 1.0, &mut row).unwrap();
            let states = row.repeat(8);
            let bytes =
                m3::MixtureEuler3d::required_workspace_bytes_with_reservoirs(&gas, &grid, 1)
                    .unwrap();
            let flow = m3::MixtureEuler3d::new_with_reservoirs(
                gas,
                grid,
                states,
                spatial::Boundaries3d {
                    x: axes[0],
                    y: axes[1],
                    z: axes[2],
                },
                &[m3::ReservoirBoundary {
                    side,
                    partial_densities: &[0.2, 0.8, 0.0],
                    velocity: [0.0; 3],
                    pressure_pa: 1.0,
                }],
                bytes,
                || false,
            )
            .unwrap();
            let bytes = flow
                .required_workspace_bytes_with_finite_reservoirs(1)
                .unwrap();
            let flow = flow
                .with_finite_reservoirs(
                    &[m3::FiniteReservoir {
                        side,
                        volume_m3: 0.06,
                    }],
                    bytes,
                    || false,
                )
                .unwrap();
            let flow = if let Some(thermal) = thermal {
                let bytes = flow
                    .required_workspace_bytes_with_finite_reservoir_heat(1)
                    .unwrap();
                let flow = flow
                    .with_finite_reservoir_heat(
                        &[m3::FiniteReservoirHeat { side, thermal }],
                        WallControl::default(),
                        bytes,
                        || false,
                    )
                    .unwrap();
                assert_eq!(flow.workspace_bytes(), bytes);
                flow
            } else {
                flow
            };
            Self::Spatial(flow, side)
        }
    }

    fn step(
        &mut self,
        dt: f64,
        reconstruction: Reconstruction,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), Error> {
        let control = StepControl {
            reconstruction,
            ..StepControl::default()
        };
        let (elapsed, halvings, energy_residual) = match self {
            Self::Planar(flow, _) => {
                let r = flow.step(Time::new(dt), control, cancelled)?;
                (
                    r.elapsed.value(),
                    r.halvings,
                    r.conservation_residual.energy,
                )
            }
            Self::Spatial(flow, _) => {
                let r = flow.step(Time::new(dt), control, cancelled)?;
                (
                    r.elapsed.value(),
                    r.halvings,
                    r.conservation_residual.energy,
                )
            }
        };
        assert_eq!(halvings, 0);
        close(elapsed, dt, 1e-15);
        close(energy_residual, 0.0, 2e-13);
        Ok(())
    }

    fn tank(&self) -> &Plenum {
        match self {
            Self::Planar(flow, side) => flow.finite_reservoir(*side).unwrap(),
            Self::Spatial(flow, side) => flow.finite_reservoir(*side).unwrap(),
        }
    }

    fn thermal(&self) -> Option<WallThermal> {
        match self {
            Self::Planar(flow, side) => flow.finite_reservoir_thermal(*side),
            Self::Spatial(flow, side) => flow.finite_reservoir_thermal(*side),
        }
    }

    fn heat(&self) -> FiniteReservoirHeatReport {
        match self {
            Self::Planar(flow, side) => flow.finite_reservoir_report(*side).unwrap().heat.unwrap(),
            Self::Spatial(flow, side) => flow.finite_reservoir_report(*side).unwrap().heat.unwrap(),
        }
    }

    fn pressure(&self) -> f64 {
        match self {
            Self::Planar(flow, side) => flow.reservoir(*side).unwrap().pressure_pa,
            Self::Spatial(flow, side) => flow.reservoir(*side).unwrap().pressure_pa,
        }
    }

    fn states(&self) -> &[f64] {
        match self {
            Self::Planar(flow, _) => flow.states(),
            Self::Spatial(flow, _) => flow.states(),
        }
    }

    fn field_mass(&self) -> f64 {
        match self {
            Self::Planar(flow, _) => flow.totals().species_masses.iter().sum(),
            Self::Spatial(flow, _) => flow.totals().species_masses.iter().sum(),
        }
    }

    fn storage(&self) -> ([f64; 3], f64) {
        let (species, energy) = match self {
            Self::Planar(flow, _) => (flow.totals().species_masses, flow.totals().energy),
            Self::Spatial(flow, _) => (flow.totals().species_masses, flow.totals().energy),
        };
        let tank = self.tank();
        let q = tank.state();
        let wall = match self.thermal() {
            Some(WallThermal::Lumped {
                temperature_k,
                heat_capacity_j_per_k,
                ..
            }) => temperature_k * heat_capacity_j_per_k,
            _ => 0.0,
        };
        (
            std::array::from_fn(|s| species[s] + q.partial_densities[s] * tank.volume_m3()),
            energy + q.energy * tank.volume_m3() + wall,
        )
    }

    fn snapshot(&self) -> Snapshot {
        match self {
            Self::Planar(flow, side) => Snapshot::Planar(
                flow.states().to_vec(),
                flow.finite_reservoir(*side).unwrap().clone(),
                flow.finite_reservoir_thermal(*side),
                flow.finite_reservoir_report(*side).copied(),
                flow.time(),
                flow.last_step().copied(),
                flow.workspace_bytes(),
            ),
            Self::Spatial(flow, side) => Snapshot::Spatial(
                flow.states().to_vec(),
                flow.finite_reservoir(*side).unwrap().clone(),
                flow.finite_reservoir_thermal(*side),
                flow.finite_reservoir_report(*side).copied(),
                flow.time(),
                flow.last_step().copied(),
                flow.workspace_bytes(),
            ),
        }
    }
}

fn wall(temperature_k: f64, conductance_w_per_k: f64) -> WallThermal {
    WallThermal::Lumped {
        temperature_k,
        heat_capacity_j_per_k: 0.7,
        conductance_w_per_k,
    }
}

#[test]
fn g0_cartesian_finite_wall_heating_and_cooling_drive_pressure_and_conserve_all_axes() {
    for dimensions in [2, 3] {
        for slot in 0..2 * dimensions {
            for reconstruction in [Reconstruction::FirstOrder, Reconstruction::Muscl] {
                for temperature in [0.6, 5.0] {
                    let mut flow = Flow::new(dimensions, slot, Some(wall(temperature, 4.0)));
                    let initial = flow.storage();
                    let field_mass = flow.field_mass();
                    for _ in 0..3 {
                        flow.step(0.0002, reconstruction, &mut || false).unwrap();
                        let now = flow.storage();
                        for s in 0..3 {
                            close(now.0[s], initial.0[s], 3e-13);
                        }
                        close(now.1, initial.1, 3e-11);
                        let heat = flow.heat();
                        assert_eq!(heat.after, flow.thermal().unwrap());
                        close(heat.source.elapsed.value(), 0.0002, 1e-15);
                        close(
                            heat.source.heat_to_gas_j + heat.source.wall_energy_change_j,
                            0.0,
                            3e-11,
                        );
                        close(heat.gas_wall_energy_residual_j, 0.0, 3e-11);
                        assert_eq!(heat.source.external_heat_to_gas_j, 0.0);
                        assert_eq!(heat.source.mechanical_work_to_gas_j, 0.0);
                        assert_eq!(heat.source.wall_impulse_ns, 0.0);
                        assert!(if temperature > 2.0 {
                            heat.source.heat_to_gas_j > 0.0
                        } else {
                            heat.source.heat_to_gas_j < 0.0
                        });
                        assert_eq!(flow.tank().state().partial_densities[2], 0.0);
                    }
                    assert!(if temperature > 2.0 {
                        flow.pressure() > 1.0
                    } else {
                        flow.pressure() < 1.0
                    });
                    assert!(if temperature > 2.0 {
                        flow.field_mass() > field_mass
                    } else {
                        flow.field_mass() < field_mass
                    });
                }
            }
        }
    }
}

#[test]
fn g0_cartesian_prescribed_bath_heat_is_external_and_adiabatic_paths_preserve_bits() {
    for dimensions in [2, 3] {
        for reconstruction in [Reconstruction::FirstOrder, Reconstruction::Muscl] {
            let mut bath = Flow::new(
                dimensions,
                1,
                Some(WallThermal::Reservoir {
                    temperature_k: 5.0,
                    conductance_w_per_k: 4.0,
                }),
            );
            let initial = bath.storage();
            let mut external = 0.0;
            for _ in 0..3 {
                bath.step(0.0002, reconstruction, &mut || false).unwrap();
                let heat = bath.heat();
                external += heat.source.external_heat_to_gas_j;
                assert!(heat.source.external_heat_to_gas_j > 0.0);
                assert_eq!(heat.source.wall_energy_change_j, 0.0);
                assert_eq!(heat.before, heat.after);
                close(bath.storage().1, initial.1 + external, 3e-11);
            }
            for thermal in [WallThermal::Adiabatic, wall(5.0, 0.0)] {
                let mut plain = Flow::new(dimensions, 0, None);
                let mut bound = Flow::new(dimensions, 0, Some(thermal));
                for _ in 0..3 {
                    plain.step(0.0002, reconstruction, &mut || false).unwrap();
                    bound.step(0.0002, reconstruction, &mut || false).unwrap();
                    assert_eq!(bound.states(), plain.states());
                    assert_eq!(bound.tank(), plain.tank());
                    assert_eq!(bound.heat().source.heat_to_gas_j, 0.0);
                    assert_eq!(bound.thermal(), Some(thermal));
                }
            }
        }
    }
}

#[test]
fn g4_cartesian_heated_tank_final_poll_keeps_wall_gas_receipts_and_deterministic_retry() {
    for dimensions in [2, 3] {
        let make = || {
            let mut flow = Flow::new(dimensions, 2 * dimensions - 1, Some(wall(5.0, 4.0)));
            flow.step(0.0002, Reconstruction::Muscl, &mut || false)
                .unwrap();
            flow
        };
        let mut probe = make();
        let mut polls = 0;
        probe
            .step(0.0002, Reconstruction::Muscl, &mut || {
                polls += 1;
                false
            })
            .unwrap();
        let mut flow = make();
        let before = flow.snapshot();
        let mut calls = 0;
        assert_eq!(
            flow.step(0.0002, Reconstruction::Muscl, &mut || {
                calls += 1;
                calls == polls
            }),
            Err(Error::Cancelled)
        );
        assert_eq!(calls, polls);
        assert_eq!(flow.snapshot(), before);
        flow.step(0.0002, Reconstruction::Muscl, &mut || false)
            .unwrap();
        assert_eq!(flow.snapshot(), probe.snapshot());
    }
}

#[test]
fn g2_finite_wall_candidates_and_split_restores_reuse_prepared_storage() {
    let gas = standard(0.2);
    let initial = gas.conserved(&[0.2, 0.8, 0.0], 0.0, 1.0).unwrap();
    let mut owner = FiniteReservoirs::<3, 6>::empty();
    owner
        .bind(&gas, 0, 0.06, initial.as_cell(), &mut || false)
        .unwrap();
    owner
        .bind_heat(&gas, 0, wall(5.0, 4.0), WallControl::default(), &mut || {
            false
        })
        .unwrap();
    let pointers = |owner: &FiniteReservoirs<3, 6>| {
        let mut values = [
            owner.accepted[0]
                .as_ref()
                .unwrap()
                .state()
                .partial_densities
                .as_ptr(),
            owner.trial[0]
                .as_ref()
                .unwrap()
                .state()
                .partial_densities
                .as_ptr(),
            owner.heat_base[0]
                .as_ref()
                .unwrap()
                .state()
                .partial_densities
                .as_ptr(),
        ];
        values.sort_unstable();
        values
    };
    let buffers = pointers(&owner);
    let bytes = owner.workspace_bytes();
    let before = owner.try_clone(&gas, &mut || false).unwrap();
    let faces = [1, 0, 0, 0, 0, 0];
    let zero = [0.0; 7];
    for _ in 0..3 {
        owner.begin_heat(&gas, 0.0002, &mut || false).unwrap();
        owner.stage(&gas, &zero, faces, 2.0, &mut || false).unwrap();
        owner.stage(&gas, &zero, faces, 1.0, &mut || false).unwrap();
        owner.finish_heat(&gas, 0.0002, &mut || false).unwrap();
        owner.publish();
        assert_eq!(pointers(&owner), buffers);
        assert_eq!(owner.workspace_bytes(), bytes);
    }
    let warmed = owner.inventory(0).unwrap().clone();
    let report = owner.report(0).copied();
    let thermal = owner.thermal(0);
    let mut probe = owner.try_clone(&gas, &mut || false).unwrap();
    let mut polls = 0;
    probe
        .copy_accepted_from(&gas, &before, &mut || {
            polls += 1;
            false
        })
        .unwrap();
    let mut calls = 0;
    assert_eq!(
        owner.copy_accepted_from(&gas, &before, &mut || {
            calls += 1;
            calls == polls
        }),
        Err(Error::Cancelled)
    );
    assert_eq!(owner.inventory(0), Some(&warmed));
    assert_eq!(owner.report(0).copied(), report);
    assert_eq!(owner.thermal(0), thermal);
    owner
        .copy_accepted_from(&gas, &before, &mut || false)
        .unwrap();
    assert_eq!(owner.inventory(0), before.inventory(0));
    assert_eq!(owner.thermal(0), before.thermal(0));
    assert_eq!(owner.report(0), None);
    assert_eq!(pointers(&owner), buffers);
    assert_eq!(owner.workspace_bytes(), bytes);
}

#[test]
fn g4_finite_tank_heat_binding_respects_total_budget_and_original_source_admission() {
    for dimensions in [2, 3] {
        let flow = Flow::new(dimensions, 0, None);
        match flow {
            Flow::Planar(flow, side) => {
                let bytes = flow
                    .required_workspace_bytes_with_finite_reservoir_heat(1)
                    .unwrap();
                assert!(
                    flow.with_finite_reservoir_heat(
                        &[m2::FiniteReservoirHeat {
                            side,
                            thermal: wall(5.0, 4.0)
                        }],
                        WallControl::default(),
                        bytes - 1,
                        || false,
                    )
                    .is_err()
                );
            }
            Flow::Spatial(flow, side) => {
                let bytes = flow
                    .required_workspace_bytes_with_finite_reservoir_heat(1)
                    .unwrap();
                assert!(
                    flow.with_finite_reservoir_heat(
                        &[m3::FiniteReservoirHeat {
                            side,
                            thermal: wall(5.0, 4.0)
                        }],
                        WallControl::default(),
                        bytes - 1,
                        || false,
                    )
                    .is_err()
                );
            }
        }
    }
    let gas = standard(0.2);
    let initial = gas.conserved(&[0.2, 0.8, 0.0], 0.0, 1.0).unwrap();
    let mut owner = FiniteReservoirs::<2, 4>::empty();
    assert!(
        owner
            .bind_heat(&gas, 0, wall(5.0, 4.0), WallControl::default(), &mut || {
                false
            })
            .is_err()
    );
    owner
        .bind(&gas, 0, 0.06, initial.as_cell(), &mut || false)
        .unwrap();
    let bytes = owner.workspace_bytes();
    for thermal in [
        wall(-1.0, 4.0),
        WallThermal::Caloric {
            mass_kg: 1.0,
            specific_energy_j_kg: 1.0,
            conductance_w_per_k: 4.0,
        },
    ] {
        assert!(matches!(
            owner.bind_heat(&gas, 0, thermal, WallControl::default(), &mut || false),
            Err(Error::Thermal(_))
        ));
        assert_eq!(owner.thermal(0), None);
        assert_eq!(owner.workspace_bytes(), bytes);
        assert_eq!(owner.inventory(0).unwrap().state(), initial.as_cell());
    }
    owner
        .bind_heat(&gas, 0, wall(5.0, 4.0), WallControl::default(), &mut || {
            false
        })
        .unwrap();
    assert!(
        owner
            .bind_heat(&gas, 0, wall(5.0, 4.0), WallControl::default(), &mut || {
                false
            })
            .is_err()
    );
}

#[test]
fn g1_heated_pressure_cfl_retries_restart_both_tank_and_wall_from_accepted() {
    for dimensions in [2, 3] {
        let make = || Flow::new(dimensions, 0, Some(wall(40.0, 1000.0)));
        let mut attempted = make();
        let control = StepControl {
            reconstruction: Reconstruction::FirstOrder,
            ..StepControl::default()
        };
        let (dt, halvings) = match &mut attempted {
            Flow::Planar(flow, _) => {
                let r = flow.step(Time::new(0.1), control, || false).unwrap();
                assert!(r.maximum_cfl <= control.cfl * (1.0 + 1e-12));
                (r.elapsed.value(), r.halvings)
            }
            Flow::Spatial(flow, _) => {
                let r = flow.step(Time::new(0.1), control, || false).unwrap();
                assert!(r.maximum_cfl <= control.cfl * (1.0 + 1e-12));
                (r.elapsed.value(), r.halvings)
            }
        };
        assert!(
            halvings > 0,
            "heating must tighten the accepted cold-state wave bound"
        );
        let mut direct = make();
        direct
            .step(dt, Reconstruction::FirstOrder, &mut || false)
            .unwrap();
        assert_eq!(attempted.states(), direct.states());
        assert_eq!(attempted.tank(), direct.tank());
        assert_eq!(attempted.thermal(), direct.thermal());
        assert_eq!(attempted.heat(), direct.heat());
    }
}
