use super::*;
use crate::thermal::boundary::conduction::{WallCell, WallConductionControls, WallConductionLink};
use crate::thermal::open::immersed::ImmersedThermalBody3;

fn retained(field: &ThermalOpenField3) -> (Vec<u64>, String, usize, Option<ThermalOpenBalance3>) {
    (
        bits(field),
        format!("{:?}", field.history),
        field.completed_steps(),
        field.last_balance(),
    )
}

fn finite_state(field: &ThermalOpenField3, face: usize) -> Vec<[u64; 3]> {
    field
        .wall_states(face)
        .unwrap()
        .iter()
        .map(|wall| {
            let ThermalWall::FiniteSolid(wall) = wall else {
                panic!("finite wall")
            };
            [
                wall.temperature.to_bits(),
                wall.heat_capacity.to_bits(),
                wall.resistance.to_bits(),
            ]
        })
        .collect()
}

#[test]
fn transient_temperature_flux_and_ambient_use_original_wall_exchange_without_reset() {
    let dims = [5, 4, 3];
    let cells: usize = dims.into_iter().product();
    let area = (dims[0] * dims[2]) as f64;
    let mut boundaries = ThermalOpenBoundaries3 {
        faces: [ThermalOpenFace3::Periodic; 6],
        internal_wall: ThermalWall::Adiabatic,
    };
    boundaries.faces[2] = ThermalOpenFace3::Wall(ThermalWall::Temperature(0.5));
    boundaries.faces[3] = ThermalOpenFace3::Wall(ThermalWall::Adiabatic);
    let mut field = field(dims, [0.02, 0.0, 0.0], 0.5, 0.83, boundaries);
    let heat = vec![0.0; cells];
    step(&mut field, &heat, None);
    for law in [
        ThermalWall::Temperature(0.8),
        ThermalWall::HeatFlux(0.002),
        ThermalWall::HeatFlux(-0.001),
        ThermalWall::Convective {
            transfer: 0.05,
            ambient: 1.2,
        },
        ThermalWall::Convective {
            transfer: 0.03,
            ambient: -0.5,
        },
        ThermalWall::Adiabatic,
    ] {
        let before = retained(&field);
        let energy: f64 = field.populations().iter().flatten().sum();
        field.set_wall_boundary(2, law).unwrap();
        assert_eq!(retained(&field), before);
        assert_eq!(field.boundaries().faces[2], ThermalOpenFace3::Wall(law));
        assert!(
            field
                .wall_states(2)
                .unwrap()
                .iter()
                .all(|wall| *wall == law)
        );
        let balance = step(&mut field, &heat, None);
        let wall = balance.walls[2].unwrap();
        assert!(wall.solid.is_none());
        let expected = match law {
            ThermalWall::Temperature(temperature) => {
                near(wall.mean_temperature, temperature, 2e-15);
                assert!(wall.heat_in > 0.0);
                wall.heat_in
            }
            ThermalWall::HeatFlux(flux) => area * flux,
            ThermalWall::Convective { transfer, ambient } => {
                area * transfer * (ambient - wall.mean_temperature)
            }
            ThermalWall::Adiabatic => 0.0,
            ThermalWall::FiniteSolid(_) => unreachable!(),
        };
        near(wall.heat_in, expected, 2e-14);
        near(
            field.populations().iter().flatten().sum::<f64>() - energy,
            expected,
            3e-13,
        );
        near(balance.combined_defect, 0.0, 2e-13);
    }
    let before = (retained(&field), field.boundaries(), field.walls.clone());
    for (face, law) in [
        (0, ThermalWall::HeatFlux(0.1)),
        (6, ThermalWall::Adiabatic),
        (2, ThermalWall::Temperature(11.0)),
        (2, ThermalWall::HeatFlux(f64::NAN)),
        (
            2,
            ThermalWall::Convective {
                transfer: -0.1,
                ambient: 0.5,
            },
        ),
        (
            2,
            ThermalWall::FiniteSolid(FiniteThermalWall {
                temperature: 0.5,
                heat_capacity: 1.0,
                resistance: 0.0,
                heat_per_step: 0.0,
            }),
        ),
    ] {
        assert!(field.set_wall_boundary(face, law).is_err());
        assert_eq!(
            (retained(&field), field.boundaries(), field.walls.clone()),
            before
        );
    }
}

#[test]
fn transient_internal_flux_heats_only_unassigned_links_and_preserves_immersed_storage() {
    let dims = [8, 4, 3];
    let cells: usize = dims.into_iter().product();
    let assigned = 2 + 8 * (1 + 4);
    let unassigned = 5 + 8 * (1 + 4);
    let mut solid = vec![false; cells];
    solid[assigned] = true;
    solid[unassigned] = true;
    let mut field = ThermalOpenField3::try_new::<ThermalOpenError3>(
        dims,
        &solid,
        &vec![0.5; cells],
        0.83,
        ThermalOpenBoundaries3 {
            faces: [ThermalOpenFace3::Periodic; 6],
            internal_wall: ThermalWall::Adiabatic,
        },
        limits(),
        0,
        |_| Ok([0.0; 3]),
        || Ok(()),
    )
    .unwrap();
    field
        .set_immersed_thermal_bodies(
            &[ImmersedThermalBody3 {
                solid_cells: &[assigned],
                solid: FiniteThermalWall {
                    temperature: 0.5,
                    heat_capacity: 0.8,
                    resistance: 0.3,
                    heat_per_step: 0.0,
                },
            }],
            64 * 1024,
        )
        .unwrap();
    let heat = vec![0.0; cells];
    for flux in [0.003, -0.001, 0.0] {
        let before = retained(&field);
        let body = field.immersed_body(0).unwrap();
        let energy: f64 = field.populations().iter().flatten().sum::<f64>()
            + body.solid.heat_capacity * body.solid.temperature;
        field
            .set_internal_wall_boundary(ThermalWall::HeatFlux(flux))
            .unwrap();
        assert_eq!(retained(&field), before);
        assert_eq!(field.immersed_body(0).unwrap(), body);
        let balance = step(&mut field, &heat, None);
        let body = field.immersed_body(0).unwrap();
        near(
            field.populations().iter().flatten().sum::<f64>()
                + body.solid.heat_capacity * body.solid.temperature
                - energy,
            6.0 * flux,
            3e-13,
        );
        near(
            balance.internal_wall.unwrap().heat_in - body.last_exchange.unwrap().heat_in,
            6.0 * flux,
            2e-14,
        );
        near(balance.combined_defect, 0.0, 2e-13);
    }
    let before = (retained(&field), field.boundaries(), field.immersed_body(0));
    assert!(
        field
            .set_internal_wall_boundary(ThermalWall::FiniteSolid(
                field.immersed_body(0).unwrap().solid,
            ))
            .is_err()
    );
    assert_eq!(
        (retained(&field), field.boundaries(), field.immersed_body(0)),
        before
    );
}

#[test]
fn spatial_wall_power_preserves_conduction_and_closes_the_original_solid_fluid_balance() {
    let mut field = heated_fixture();
    let cells = field.populations().len();
    let plane = field.wall_states(2).unwrap().len();
    field
        .set_wall_conduction(
            1,
            &[WallConductionLink {
                a: WallCell { side: 0, link: 0 },
                b: WallCell {
                    side: 1,
                    link: plane - 1,
                },
                conductance: 3.0,
            }],
            WallConductionControls {
                max_nodes: 2,
                max_links: 1,
                max_factor_bytes: 1024 * 1024,
                balance_tolerance: 1e-10,
            },
        )
        .unwrap();
    let declaration = field.boundaries();
    let heat = vec![0.0; cells];
    for pulse in [1.0, -0.5, 0.0] {
        let before = retained(&field);
        let physical = [finite_state(&field, 2), finite_state(&field, 3)];
        let mut lower = vec![0.0; plane];
        let mut upper = vec![0.0; plane];
        lower[0] = pulse * 0.02;
        upper[plane - 1] = pulse * -0.004;
        field.set_finite_wall_heat(2, &lower).unwrap();
        field.set_finite_wall_heat(3, &upper).unwrap();
        assert_eq!(retained(&field), before);
        assert_eq!([finite_state(&field, 2), finite_state(&field, 3)], physical);
        assert_eq!(field.boundaries(), declaration);
        assert!(field.wall_conduction(1).is_some());
        let energy_before: f64 = field.populations().iter().flatten().sum();
        let balance = step(&mut field, &heat, None);
        let bottom = balance.walls[2].unwrap().solid.unwrap();
        let top = balance.walls[3].unwrap().solid.unwrap();
        let power = lower.iter().chain(&upper).sum::<f64>();
        near(bottom.source + top.source, power, 2e-15);
        near(bottom.conducted_in + top.conducted_in, 0.0, 2e-14);
        near(
            field.populations().iter().flatten().sum::<f64>() - energy_before
                + bottom.stored_change
                + top.stored_change,
            power + balance.reservoirs.iter().map(|r| r.heat_in).sum::<f64>(),
            3e-13,
        );
        near(balance.combined_defect, 0.0, 2e-13);
        if pulse == 1.0 {
            // Independent two-node implicit conduction solve. Initially the
            // normal outgoing populations see uniform fluid temperature 0.4.
            let k = 1.0 / 4.4;
            let diagonal = 2.7 + k + 3.0;
            let rhs_a = 2.7 * 1.2 + 0.02 + k * 0.4;
            let rhs_b = 2.7 * 0.7 - 0.004 + k * 0.4;
            let determinant = diagonal * diagonal - 9.0;
            let expected_a = (diagonal * rhs_a + 3.0 * rhs_b) / determinant;
            let expected_b = (diagonal * rhs_b + 3.0 * rhs_a) / determinant;
            let ThermalWall::FiniteSolid(a) = field.wall_states(2).unwrap()[0] else {
                panic!("finite cell")
            };
            let ThermalWall::FiniteSolid(b) = field.wall_states(3).unwrap()[plane - 1] else {
                panic!("finite cell")
            };
            near(a.temperature, expected_a, 3e-14);
            near(b.temperature, expected_b, 3e-14);
        }
    }
    let before = (retained(&field), field.walls.clone());
    for source in [vec![0.0; plane - 1], {
        let mut values = vec![0.001; plane];
        values[plane - 1] = f64::NAN;
        values
    }] {
        assert!(field.set_finite_wall_heat(2, &source).is_err());
        assert_eq!((retained(&field), field.walls.clone()), before);
    }
    assert!(
        field
            .set_wall_boundary(2, ThermalWall::Temperature(0.5))
            .is_err()
    );
    assert_eq!((retained(&field), field.walls.clone()), before);
}

#[test]
fn wall_power_rejects_blocked_slots_on_all_noncubic_faces_before_changing_any_source() {
    let dims = [5, 4, 3];
    let cells: usize = dims.into_iter().product();
    let mut solid = vec![false; cells];
    solid[0] = true;
    solid[cells - 1] = true;
    let wall = ThermalWall::FiniteSolid(FiniteThermalWall {
        temperature: 0.5,
        heat_capacity: 2.0,
        resistance: 0.3,
        heat_per_step: 0.0,
    });
    let mut field = ThermalOpenField3::try_new::<ThermalOpenError3>(
        dims,
        &solid,
        &vec![0.5; cells],
        0.83,
        ThermalOpenBoundaries3 {
            faces: [ThermalOpenFace3::Wall(wall); 6],
            internal_wall: ThermalWall::Adiabatic,
        },
        limits(),
        0,
        |_| Ok([0.0; 3]),
        || Ok(()),
    )
    .unwrap();
    let mut total_source = 0.0;
    for face in 0..6 {
        let plane = cells / dims[face / 2];
        let blocked = if face % 2 == 0 { 0 } else { plane - 1 };
        let mut source = vec![0.0001; plane];
        let before = field.walls.clone();
        assert!(field.set_finite_wall_heat(face, &source).is_err());
        assert_eq!(field.walls, before);
        source[blocked] = 0.0;
        field.set_finite_wall_heat(face, &source).unwrap();
        total_source += source.iter().sum::<f64>();
    }
    let balance = step(&mut field, &vec![0.0; cells], None);
    let walls: Vec<_> = balance
        .walls
        .into_iter()
        .flatten()
        .map(|w| w.solid.unwrap())
        .collect();
    near(walls.iter().map(|w| w.source).sum(), total_source, 2e-15);
    near(
        balance.stored_change + walls.iter().map(|w| w.stored_change).sum::<f64>(),
        total_source,
        2e-13,
    );
}
