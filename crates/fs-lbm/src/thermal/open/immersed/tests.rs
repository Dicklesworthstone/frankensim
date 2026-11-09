use super::*;
use crate::thermal::boundary::ThermalWall;
use crate::thermal::open::{ThermalOpenBalance3, ThermalOpenBoundaries3};

fn fixture(solid_cells: &[usize], source: f64, varying: bool) -> ThermalOpenField3 {
    let dims = [5, 4, 3];
    let mut solid = vec![false; 60];
    for &cell in solid_cells {
        solid[cell] = true;
    }
    let initial: Vec<_> = (0..60)
        .map(|cell| 0.2 + if varying { 0.001 * cell as f64 } else { 0.0 })
        .collect();
    let mut field = ThermalOpenField3::try_new::<ThermalOpenError3>(
        dims,
        &solid,
        &initial,
        0.85,
        ThermalOpenBoundaries3 {
            faces: [ThermalOpenFace3::Periodic; 6],
            internal_wall: ThermalWall::Adiabatic,
        },
        ThermalOpenLimits3 {
            maximum_cells: 60,
            maximum_storage_bytes: 1024 * 1024,
            temperature_range: [-2.0, 3.0],
            balance_tolerance: 1e-11,
        },
        0,
        |_| Ok([0.0; 3]),
        || Ok(()),
    )
    .unwrap();
    field
        .set_immersed_thermal_bodies(
            &[ImmersedThermalBody3 {
                solid_cells,
                solid: FiniteThermalWall {
                    temperature: 1.4,
                    heat_capacity: 4.0,
                    resistance: 0.7,
                    heat_per_step: source,
                },
            }],
            64 * 1024,
        )
        .unwrap();
    field
}

fn step(field: &mut ThermalOpenField3) -> ThermalOpenBalance3 {
    let heat = [0.0; 60];
    let candidate = field
        .begin_step::<ThermalOpenError3>(&heat, None, || Ok(()))
        .unwrap();
    let balance = field
        .finish_step::<ThermalOpenError3>(&heat, candidate, || Ok(()))
        .unwrap();
    field.publish(field.completed_steps() + 1, balance);
    balance
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.17e} != {b:.17e}");
}
fn fluid_heat(field: &ThermalOpenField3) -> f64 {
    field.populations().iter().flatten().sum()
}

#[test]
fn two_voxel_body_uses_one_capacity_and_source_across_all_ten_faces() {
    // Adjacent cells (2,1,1) and (3,1,1) have ten exposed normal faces.
    // These independently enumerated neighbors also check D3Q7 directions.
    let neighbors = [26, 22, 32, 7, 47, 29, 23, 33, 8, 48];
    let trace_sum: f64 = neighbors
        .map(|cell| 0.2 + 0.001 * f64::from(cell))
        .iter()
        .sum();
    let k = 1.0 / 4.7;
    let expected = (4.0 * 1.4 + 0.06 + k * trace_sum) / (4.0 + 10.0 * k);
    // An unassigned obstacle exercises the stateless fallback between managed
    // links without advancing the shared body's streaming cursor.
    let mut field = fixture(&[0, 27, 28], 0.06, true);
    let solid = field.immersed_body(0).unwrap().solid;
    field
        .set_immersed_thermal_bodies(
            &[ImmersedThermalBody3 {
                solid_cells: &[27, 28],
                solid,
            }],
            64 * 1024,
        )
        .unwrap();
    let before = fluid_heat(&field);
    let balance = step(&mut field);
    let body = field.immersed_body(0).unwrap();
    assert_eq!(body.interface_links, 10);
    near(body.solid.temperature, expected, 3e-15);
    let exchange = body.last_exchange.unwrap();
    assert_eq!(exchange.source, 0.06);
    near(exchange.heat_in, k * (10.0 * expected - trace_sum), 3e-15);
    near(exchange.stored_change, 4.0 * (expected - 1.4), 3e-15);
    near(
        fluid_heat(&field) - before + exchange.stored_change,
        0.06,
        1e-14,
    );
    assert_eq!(balance.internal_wall.unwrap().solid.unwrap(), exchange);
    near(balance.combined_defect, 0.0, 5e-15);
}

#[test]
fn finite_body_and_closed_fluid_reach_the_capacity_weighted_equilibrium() {
    let mut field = fixture(&[27, 28], 0.0, false);
    let initial = fluid_heat(&field) + 4.0 * 1.4;
    let expected = (58.0 * 0.2 + 4.0 * 1.4) / 62.0;
    for _ in 0..1000 {
        let balance = step(&mut field);
        near(balance.combined_defect, 0.0, 1e-13);
    }
    let body = field.immersed_body(0).unwrap();
    near(body.solid.temperature, expected, 2e-12);
    for cell in 0..60 {
        if let Some(temperature) = field.temperature(cell) {
            near(temperature, expected, 2e-12);
        }
    }
    near(
        fluid_heat(&field) + 4.0 * body.solid.temperature,
        initial,
        3e-12,
    );
}

#[test]
fn neighboring_hot_and_cold_bodies_keep_separate_inventories_and_sources() {
    let mut field = fixture(&[27, 28], 0.0, false);
    let hot = FiniteThermalWall {
        temperature: 1.4,
        heat_capacity: 4.0,
        resistance: 0.7,
        heat_per_step: 0.01,
    };
    let cold = FiniteThermalWall {
        temperature: -0.2,
        heat_capacity: 2.0,
        resistance: 0.7,
        heat_per_step: -0.005,
    };
    field
        .set_immersed_thermal_bodies(
            &[
                ImmersedThermalBody3 {
                    solid_cells: &[27],
                    solid: hot,
                },
                ImmersedThermalBody3 {
                    solid_cells: &[28],
                    solid: cold,
                },
            ],
            64 * 1024,
        )
        .unwrap();
    let before = fluid_heat(&field);
    let balance = step(&mut field);
    let mut stored = 0.0;
    for (i, original) in [hot, cold].into_iter().enumerate() {
        let body = field.immersed_body(i).unwrap();
        assert_eq!(body.interface_links, 5);
        let k = 5.0 / 4.7;
        let expected =
            (original.heat_capacity * original.temperature + original.heat_per_step + k * 0.2)
                / (original.heat_capacity + k);
        near(body.solid.temperature, expected, 2e-15);
        let exchange = body.last_exchange.unwrap();
        assert_eq!(exchange.source, original.heat_per_step);
        assert_eq!(exchange.heat_in > 0.0, i == 0);
        near(exchange.defect(), 0.0, 4e-15);
        stored += exchange.stored_change;
    }
    near(fluid_heat(&field) - before + stored, 0.005, 1e-14);
    near(balance.combined_defect, 0.0, 5e-15);
}

#[test]
fn periodic_seams_and_authored_cell_order_preserve_shared_exchange() {
    // Straddle each seam of the noncubic grid. The shared face is internal,
    // including z where both exterior faces lead to the same fluid voxel.
    for other in [4, 15, 40] {
        let mut actual = fixture(&[0, other], 0.002, true);
        let mut reversed = fixture(&[other, 0], 0.002, true);
        assert_eq!(actual.immersed_body(0).unwrap().interface_links, 10);
        for _ in 0..12 {
            assert_eq!(step(&mut actual), step(&mut reversed));
            assert_eq!(actual.populations(), reversed.populations());
            assert_eq!(actual.immersed_body(0), reversed.immersed_body(0));
        }
    }
}

#[test]
fn extreme_capacities_preserve_equilibrium_and_refuse_unresolved_storage() {
    for capacity in [f64::from_bits(1), 1e-24, 1.0, 1e24, f64::MAX] {
        let mut field = fixture(&[27, 28], 0.0, false);
        for (cell, populations) in field.populations.iter_mut().enumerate() {
            if !field.solid[cell] {
                *populations = D3Q7.equilibrium(0.5, [0.0; 3]);
            }
        }
        field
            .set_immersed_thermal_bodies(
                &[ImmersedThermalBody3 {
                    solid_cells: &[27, 28],
                    solid: FiniteThermalWall {
                        temperature: 0.5,
                        heat_capacity: capacity,
                        resistance: 0.7,
                        heat_per_step: 0.0,
                    },
                }],
                64 * 1024,
            )
            .unwrap();
        let exchange = step(&mut field).internal_wall.unwrap().solid.unwrap();
        assert_eq!(
            field.immersed_body(0).unwrap().solid.temperature.to_bits(),
            0.5_f64.to_bits()
        );
        assert_eq!(exchange.stored_change, 0.0);
        assert_eq!(exchange.heat_in, 0.0);
        assert_eq!(exchange.defect(), 0.0);
    }
    let mut field = fixture(&[27, 28], 0.0, false);
    let solid = FiniteThermalWall {
        heat_capacity: f64::MAX,
        ..field.immersed_body(0).unwrap().solid
    };
    field
        .set_immersed_thermal_bodies(
            &[ImmersedThermalBody3 {
                solid_cells: &[27, 28],
                solid,
            }],
            64 * 1024,
        )
        .unwrap();
    // A finite outgoing heat cannot be silently accepted when the required
    // temperature decrement lies below resolution of the huge capacity.
    let before = field.populations().to_vec();
    assert!(matches!(
        field.begin_step::<ThermalOpenError3>(&[0.0; 60], None, || Ok(())),
        Err(ThermalOpenError3::Balance { .. })
    ));
    assert_eq!(field.populations(), before);
    assert_eq!(field.immersed_body(0).unwrap().solid, solid);
    assert_eq!(field.completed_steps(), 0);
}

#[test]
fn late_refusal_and_cancellation_keep_body_storage_receipts_and_exact_retry() {
    let mut actual = fixture(&[27, 28], 0.001, true);
    let mut control = fixture(&[27, 28], 0.001, true);
    step(&mut actual);
    step(&mut control);
    let accepted = actual.immersed_body(0);
    let populations = actual.populations().to_vec();
    let balance = actual.last_balance();
    let history = format!("{:?}", actual.history);
    let heat = [0.0; 60];
    let candidate = actual
        .begin_step::<ThermalOpenError3>(&heat, None, || Ok(()))
        .unwrap();
    actual.velocities[0][0] = 0.25;
    assert!(matches!(
        actual.finish_step::<ThermalOpenError3>(&heat, candidate, || Ok(())),
        Err(ThermalOpenError3::Velocity { .. })
    ));
    actual.velocities[0][0] = 0.0;
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
    assert_eq!(actual.immersed_body(0), accepted);
    assert_eq!(actual.populations(), populations);
    assert_eq!(actual.last_balance(), balance);
    assert_eq!(actual.completed_steps(), 1);
    assert_eq!(format!("{:?}", actual.history), history);
    assert_eq!(step(&mut actual), step(&mut control));
    assert_eq!(actual.populations(), control.populations());
    assert_eq!(actual.immersed_body(0), control.immersed_body(0));
}

#[test]
fn ownership_budget_and_temperature_refusals_do_not_replace_the_body() {
    let mut field = fixture(&[27, 28], 0.0, false);
    let original = field.immersed_body(0).unwrap();
    for cells in [&[27, 27][..], &[27, 1][..], &[60][..], &[][..]] {
        assert!(
            field
                .set_immersed_thermal_bodies(
                    &[ImmersedThermalBody3 {
                        solid_cells: cells,
                        solid: original.solid,
                    }],
                    64 * 1024
                )
                .is_err()
        );
        assert_eq!(field.immersed_body(0), Some(original));
    }
    let declaration = ImmersedThermalBody3 {
        solid_cells: &[27, 28],
        solid: original.solid,
    };
    assert!(
        field
            .set_immersed_thermal_bodies(&[declaration], 0)
            .is_err()
    );
    assert!(field.set_immersed_body_heat(0, f64::NAN).is_err());
    assert_eq!(field.immersed_body(0), Some(original));
    field.set_immersed_body_heat(0, 100.0).unwrap();
    assert!(matches!(
        field.begin_step::<ThermalOpenError3>(&[0.0; 60], None, || Ok(())),
        Err(ThermalOpenError3::BoundaryTemperature { .. })
    ));
    near(field.immersed_body(0).unwrap().solid.temperature, 1.4, 0.0);
    assert!(field.immersed_body(0).unwrap().last_exchange.is_none());
    field.set_immersed_body_heat(0, 0.001).unwrap();
    near(
        step(&mut field)
            .internal_wall
            .unwrap()
            .solid
            .unwrap()
            .source,
        0.001,
        0.0,
    );
}
