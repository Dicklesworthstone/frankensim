use super::*;
use crate::thermal::boundary::{FiniteThermalWall, ThermalWall};
use crate::thermal::open::{ThermalOpenBalance3, ThermalOpenBoundaries3, ThermalOpenFace3};

fn controls() -> WallConductionControls {
    WallConductionControls {
        max_nodes: 64,
        max_links: 256,
        max_factor_bytes: 1024 * 1024,
        balance_tolerance: 1e-9,
    }
}
fn field(dims: [usize; 3], solids: &[usize]) -> ThermalOpenField3 {
    let cells: usize = dims.into_iter().product();
    let mut mask = vec![false; cells];
    for &cell in solids {
        mask[cell] = true;
    }
    ThermalOpenField3::try_new::<ThermalOpenError3>(
        dims,
        &mask,
        &vec![0.2; cells],
        0.85,
        ThermalOpenBoundaries3 {
            faces: [ThermalOpenFace3::Periodic; 6],
            internal_wall: ThermalWall::Adiabatic,
        },
        ThermalOpenLimits3 {
            maximum_cells: 4096,
            maximum_storage_bytes: 2 * 1024 * 1024,
            temperature_range: [-2.0, 5.0],
            balance_tolerance: 1e-9,
        },
        0,
        |_| Ok([0.0; 3]),
        || Ok(()),
    )
    .unwrap()
}
fn solid(temperature: f64, capacity: f64, source: f64) -> FiniteThermalWall {
    FiniteThermalWall {
        temperature,
        heat_capacity: capacity,
        resistance: 0.7,
        heat_per_step: source,
    }
}
fn edge(a: usize, b: usize, conductance: f64) -> ImmersedConductionLink3 {
    ImmersedConductionLink3 { a, b, conductance }
}
fn step(field: &mut ThermalOpenField3) -> ThermalOpenBalance3 {
    let heat = vec![0.0; field.populations().len()];
    let trial = field
        .begin_step::<ThermalOpenError3>(&heat, None, || Ok(()))
        .unwrap();
    let balance = field
        .finish_step::<ThermalOpenError3>(&heat, trial, || Ok(()))
        .unwrap();
    field.publish(field.completed_steps() + 1, balance);
    balance
}
fn fluid_heat(field: &ThermalOpenField3) -> f64 {
    field.populations().iter().flatten().sum()
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.17e} != {b:.17e}");
}

fn pair(g: f64) -> ThermalOpenField3 {
    let mut field = field([5, 4, 3], &[27, 28]);
    field
        .set_immersed_thermal_network(
            &[
                ImmersedThermalBody3 {
                    solid_cells: &[27],
                    solid: solid(1.4, 4.0, 0.01),
                },
                ImmersedThermalBody3 {
                    solid_cells: &[28],
                    solid: solid(-0.2, 2.0, -0.005),
                },
            ],
            &[edge(0, 1, g)],
            64 * 1024,
            controls(),
        )
        .unwrap();
    field
}

#[test]
fn simultaneous_graph_and_fluid_endpoint_matches_independent_two_by_two_balance() {
    // The large conductance has an explicit solid timestep limit far below
    // one lattice interval; the implicit endpoint remains finite and passive.
    for g in [0.8, 10_000.0] {
        let mut field = pair(g);
        let k = 5.0 / 4.7;
        let aa = 4.0 + k + g;
        let bb = 2.0 + k + g;
        let ra = 4.0 * 1.4 + 0.01 + k * 0.2;
        let rb = 2.0 * -0.2 - 0.005 + k * 0.2;
        let determinant = aa * bb - g * g;
        let expected = [
            (bb * ra + g * rb) / determinant,
            (aa * rb + g * ra) / determinant,
        ];
        let before = fluid_heat(&field);
        let result = step(&mut field);
        let a = field.immersed_body(0).unwrap();
        let b = field.immersed_body(1).unwrap();
        near(a.solid.temperature, expected[0], 2e-11);
        near(b.solid.temperature, expected[1], 2e-11);
        let qa = a.last_exchange.unwrap();
        let qb = b.last_exchange.unwrap();
        near(qa.heat_in, k * (a.solid.temperature - 0.2), 2e-15);
        near(qb.heat_in, k * (b.solid.temperature - 0.2), 2e-15);
        near(
            qa.conducted_in,
            g * (b.solid.temperature - a.solid.temperature),
            0.0,
        );
        assert_eq!(qa.conducted_in, -qb.conducted_in);
        near(
            fluid_heat(&field) - before + qa.stored_change + qb.stored_change,
            0.005,
            2e-11,
        );
        near(result.combined_defect, 0.0, 2e-11);
        for body in [a, b] {
            assert!((-0.2..=1.4).contains(&body.solid.temperature));
        }
    }
}

fn buried_core(source: f64) -> ThermalOpenField3 {
    let mut outer = Vec::new();
    for z in 1..=3 {
        for y in 1..=3 {
            for x in 1..=3 {
                let cell = x + 5 * (y + 5 * z);
                if cell != 62 {
                    outer.push(cell);
                }
            }
        }
    }
    let mut cells = outer.clone();
    cells.push(62);
    let mut field = field([5, 5, 5], &cells);
    field
        .set_immersed_thermal_network(
            &[
                ImmersedThermalBody3 {
                    solid_cells: &[62],
                    solid: solid(1.5, 4.0, source),
                },
                ImmersedThermalBody3 {
                    solid_cells: &outer,
                    solid: solid(0.5, 26.0, 0.0),
                },
            ],
            &[edge(0, 1, 1.2)],
            128 * 1024,
            controls(),
        )
        .unwrap();
    field
}

#[test]
fn buried_source_conducts_through_original_solid_storage_before_heating_the_fluid() {
    let mut field = buried_core(0.02);
    let before = fluid_heat(&field);
    let k = 54.0 / 4.7;
    let g = 1.2;
    let aa = 4.0 + g;
    let bb = 26.0 + k + g;
    let ra = 4.0 * 1.5 + 0.02;
    let rb = 26.0 * 0.5 + k * 0.2;
    let determinant = aa * bb - g * g;
    let expected = [
        (bb * ra + g * rb) / determinant,
        (aa * rb + g * ra) / determinant,
    ];
    let balance = step(&mut field);
    let core = field.immersed_body(0).unwrap();
    let shell = field.immersed_body(1).unwrap();
    assert_eq!(core.interface_links, 0);
    assert_eq!(shell.interface_links, 54);
    near(core.solid.temperature, expected[0], 2e-15);
    near(shell.solid.temperature, expected[1], 2e-15);
    let inner = core.last_exchange.unwrap();
    let outer = shell.last_exchange.unwrap();
    assert_eq!(inner.heat_in, 0.0);
    assert_eq!(inner.source, 0.02);
    assert!(inner.conducted_in < 0.0);
    assert_eq!(inner.conducted_in, -outer.conducted_in);
    near(inner.stored_change, 0.02 + inner.conducted_in, 2e-15);
    near(
        fluid_heat(&field) - before + inner.stored_change + outer.stored_change,
        0.02,
        1e-13,
    );
    near(
        balance.internal_wall.unwrap().solid.unwrap().source,
        0.02,
        0.0,
    );
    near(balance.combined_defect, 0.0, 1e-13);
}

#[test]
fn buried_core_and_shell_reach_the_whole_capacity_weighted_equilibrium() {
    let mut field = buried_core(0.0);
    let initial = fluid_heat(&field) + 4.0 * 1.5 + 26.0 * 0.5;
    let expected = (98.0 * 0.2 + 4.0 * 1.5 + 26.0 * 0.5) / (98.0 + 4.0 + 26.0);
    for _ in 0..600 {
        step(&mut field);
    }
    let a = field.immersed_body(0).unwrap();
    let b = field.immersed_body(1).unwrap();
    near(a.solid.temperature, expected, 1e-10);
    near(b.solid.temperature, expected, 1e-10);
    for cell in 0..125 {
        if let Some(temperature) = field.temperature(cell) {
            near(temperature, expected, 1e-10);
        }
    }
    near(
        fluid_heat(&field) + 4.0 * a.solid.temperature + 26.0 * b.solid.temperature,
        initial,
        3e-11,
    );
}

#[test]
fn late_failure_preserves_graph_body_states_receipts_and_exact_retry() {
    let mut actual = pair(0.8);
    let mut control = pair(0.8);
    step(&mut actual);
    step(&mut control);
    let bodies = [actual.immersed_body(0), actual.immersed_body(1)];
    let populations = actual.populations().to_vec();
    let last = actual.last_balance();
    let heat = [0.0; 60];
    let candidate = actual
        .begin_step::<ThermalOpenError3>(&heat, None, || Ok(()))
        .unwrap();
    let mut calls = 0;
    assert!(matches!(
        actual.finish_step::<ThermalOpenError3>(&heat, candidate, || {
            calls += 1;
            if calls == 4 {
                Err(ThermalOpenError3::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(ThermalOpenError3::Cancelled)
    ));
    assert_eq!([actual.immersed_body(0), actual.immersed_body(1)], bodies);
    assert_eq!(actual.populations(), populations);
    assert_eq!(actual.last_balance(), last);
    assert_eq!(actual.completed_steps(), 1);
    assert_eq!(step(&mut actual), step(&mut control));
    assert_eq!(actual.populations(), control.populations());
    assert_eq!(
        [actual.immersed_body(0), actual.immersed_body(1)],
        [control.immersed_body(0), control.immersed_body(1)]
    );
}

#[test]
fn invalid_graphs_and_orphan_interior_nodes_preserve_the_previous_plan() {
    let mut field = pair(0.8);
    step(&mut field);
    let before = [field.immersed_body(0), field.immersed_body(1)];
    let nodes = [
        ImmersedThermalBody3 {
            solid_cells: &[27],
            solid: solid(1.0, 4.0, 0.0),
        },
        ImmersedThermalBody3 {
            solid_cells: &[28],
            solid: solid(0.5, 2.0, 0.0),
        },
    ];
    for links in [
        vec![edge(0, 0, 1.0)],
        vec![edge(0, 2, 1.0)],
        vec![edge(0, 1, -1.0)],
        vec![edge(0, 1, f64::INFINITY)],
        vec![edge(0, 1, 0.5), edge(1, 0, 0.5)],
    ] {
        assert!(
            field
                .set_immersed_thermal_network(&nodes, &links, 64 * 1024, controls())
                .is_err()
        );
        assert_eq!([field.immersed_body(0), field.immersed_body(1)], before);
    }
    let too_small = WallConductionControls {
        max_factor_bytes: 1,
        ..controls()
    };
    assert!(
        field
            .set_immersed_thermal_network(&nodes, &[edge(0, 1, 0.8)], 64 * 1024, too_small)
            .is_err()
    );
    assert_eq!([field.immersed_body(0), field.immersed_body(1)], before);
    let mut buried = buried_core(0.0);
    let before = buried.immersed_body(0);
    assert!(
        buried
            .set_immersed_thermal_network(
                &[ImmersedThermalBody3 {
                    solid_cells: &[62],
                    solid: solid(1.5, 4.0, 0.0)
                },],
                &[],
                64 * 1024,
                controls()
            )
            .is_err()
    );
    assert_eq!(buried.immersed_body(0), before);
}
