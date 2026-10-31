//! G1/G3/G4: directional material transport, bracket selection and finite inventory.
use super::*;
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::pressure::direction::{
    LiquidFlowDirection, LiquidFlowSelection,
};

fn circulation_budget() -> PressureDrivenBudget {
    PressureDrivenBudget {
        mass_flow_bracket_kg_per_s: (0.0002, 0.02),
        max_flow_evaluations: 128,
        ..pressure_budget()
    }
}

// Independent constant-Cp/HTC solution with source rho(T), in actual traversal.
fn reference(inlet: f64, reverse: bool, authored_drive: f64) -> (f64, f64, f64) {
    let walls = if reverse {
        [300.0, 400.0]
    } else {
        [400.0, 300.0]
    };
    let drive = if reverse {
        -authored_drive
    } else {
        authored_drive
    };
    let trial = |mass: f64| {
        let decay = (-0.549 / (mass * 2000.0)).exp();
        let middle = walls[0] + (inlet - walls[0]) * decay;
        let outlet = walls[1] + (middle - walls[1]) * decay;
        let down = 1000.0 - 0.5 * (0.5 * (inlet + middle) - 300.0);
        let up = 1000.0 - 0.5 * (0.5 * (middle + outlet) - 300.0);
        (
            1e7 * mass * (1.0 / down + 1.0 / up) - 9.81 * (down - up) - drive,
            outlet,
            mass * 2000.0 * (outlet - inlet),
        )
    };
    let (mut lo, mut hi) = circulation_budget().mass_flow_bracket_kg_per_s;
    assert!(trial(lo).0 < 0.0 && trial(hi).0 > 0.0);
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        if trial(mid).0 < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let mass = 0.5 * (lo + hi);
    let (_, outlet, heat) = trial(mass);
    (mass, outlet, heat)
}

#[test]
fn reverse_circulation_recomputes_upstream_material_heat_and_authored_signs() {
    with_card_path(
        &buoyant_card(),
        LiquidPropertyTemperature::FilmMean,
        |cx, path| {
            let solve = |direction| {
                path.solve_gravity_circulation(
                    cx,
                    &[400.0, 300.0],
                    &resistances(),
                    Pressure::new(0.0),
                    &closed_gravity(),
                    LiquidFlowSelection::Direction(direction),
                    local(),
                    circulation_budget(),
                )
                .unwrap()
            };
            let forward = solve(LiquidFlowDirection::Forward);
            let reverse = solve(LiquidFlowDirection::Reverse);
            let point = reverse.operating_point();
            let expected = reference(300.0, true, 0.0);
            close(reverse.signed_mass_flow_kg_per_s(), -expected.0);
            close(point.thermal.regions[1].outlet.temperature_k(), expected.1);
            close(point.thermal.total_heat_rate_w, expected.2);
            assert_eq!(point.thermal.regions[0].region, "last");
            assert_eq!(point.thermal.regions[1].region, "first");
            assert_eq!(point.thermal.regions[0].inlet, path.inlet());
            assert_eq!(
                point.thermal.regions[1].inlet,
                point.thermal.regions[0].outlet
            );
            assert!(core::ptr::eq(
                point.thermal.source,
                forward.operating_point().thermal.source
            ));
            assert!(point.hydraulics.total_dissipated_power_w > 0.0);
            assert!(point.thermal.total_heat_rate_w > 0.0);
            assert!(reverse.signed_hydrostatic_gain_pa() < 0.0);
            close(
                reverse.signed_pressure_drop_pa(),
                reverse.signed_hydrostatic_gain_pa(),
            );
            close(reverse.signed_pressure_residual_pa(), 0.0);
            assert_eq!(
                point.hydrostatic.as_ref().unwrap().gravitational_power_w(),
                0.0
            );
            assert!(
                (point.mass_flow_kg_per_s - forward.operating_point().mass_flow_kg_per_s).abs()
                    > 1e-5
            );
            assert_ne!(
                point.thermal.total_heat_rate_w,
                forward.operating_point().thermal.total_heat_rate_w
            );
            for (hydraulic, thermal) in point.hydraulics.regions.iter().zip(&point.thermal.regions)
            {
                assert_eq!(hydraulic.region, thermal.region);
                let bulk = 0.5 * (thermal.inlet.temperature_k() + thermal.outlet.temperature_k());
                close(
                    hydraulic.material.density_kg_m3(),
                    1000.0 - 0.5 * (bulk - 300.0),
                );
            }
        },
    );
}

#[test]
fn orientation_search_exposes_ambiguity_and_counts_both_pressure_brackets() {
    with_card_path(
        &buoyant_card(),
        LiquidPropertyTemperature::BulkMean,
        |cx, path| {
            let solve = |drive, budget| {
                path.solve_gravity_circulation(
                    cx,
                    &[400.0, 300.0],
                    &resistances(),
                    Pressure::new(drive),
                    &closed_gravity(),
                    LiquidFlowSelection::OnlyBracketed,
                    local(),
                    budget,
                )
            };
            assert!(
                matches!(solve(0.0,circulation_budget()),Err(PressureDrivenError::AmbiguousCirculationBrackets {
            forward_mass_flow_kg_per_s,reverse_mass_flow_kg_per_s,
        }) if forward_mass_flow_kg_per_s>0.0 && reverse_mass_flow_kg_per_s<0.0)
            );
            for (drive, direction) in [
                (200.0, LiquidFlowDirection::Forward),
                (-200.0, LiquidFlowDirection::Reverse),
            ] {
                let answer = solve(drive, circulation_budget()).unwrap();
                assert_eq!(answer.direction(), direction);
                close(answer.signed_driving_pressure_pa(), drive);
                close(
                    answer.signed_pressure_drop_pa(),
                    drive + answer.signed_hydrostatic_gain_pa(),
                );
                assert_eq!(
                    answer.total_flow_evaluations(),
                    answer.operating_point().flow_evaluations + 2
                );
                assert!(
                    answer.total_flow_evaluations() <= circulation_budget().max_flow_evaluations
                );
            }
            let forward = solve(200.0, circulation_budget()).unwrap();
            let exhausted = PressureDrivenBudget {
                max_flow_evaluations: forward.operating_point().flow_evaluations,
                ..circulation_budget()
            };
            assert!(
                matches!(solve(200.0,exhausted),Err(PressureDrivenError::CirculationWorkBudget {
            performed,maximum,
        }) if performed==maximum && performed==exhausted.max_flow_evaluations)
            );
        },
    );
}

#[test]
fn unsupported_orientation_is_not_mistaken_for_an_absent_circulation_root() {
    let card = card_with_transport_ceiling(false, 1.0, 320.0);
    with_card_path(&card, LiquidPropertyTemperature::BulkMean, |cx, path| {
        let solve = |selection| {
            path.solve_gravity_circulation(
                cx,
                &[400.0, 250.0],
                &resistances(),
                Pressure::new(-50.0),
                &closed_gravity(),
                selection,
                local(),
                circulation_budget(),
            )
        };
        // Cooling first keeps reverse-flow bulk states inside the source domain.
        assert!(solve(LiquidFlowSelection::Direction(LiquidFlowDirection::Reverse)).is_ok());
        // Heating first at the low-flow endpoint exceeds the same pinned source.
        assert!(matches!(
            solve(LiquidFlowSelection::OnlyBracketed),
            Err(PressureDrivenError::Thermal(_))
        ));
    });
    with_path(
        false,
        1.0,
        LiquidPropertyTemperature::BulkMean,
        |cx, path| {
            assert!(matches!(
                path.solve_gravity_circulation(
                    cx,
                    &[400.0, 300.0],
                    &resistances(),
                    Pressure::new(0.0),
                    &closed_gravity(),
                    LiquidFlowSelection::OnlyBracketed,
                    local(),
                    circulation_budget()
                ),
                Err(PressureDrivenError::NoCirculationBracket { .. })
            ));
        },
    );
}

#[test]
fn finite_inventory_reverses_without_reauthoring_and_preserves_energy_and_rollback() {
    with_card_path(
        &buoyant_card(),
        LiquidPropertyTemperature::FilmMean,
        |cx, path| {
            let controls = GravityReservoirControls {
                pressure: circulation_budget(),
                ..gravity_reservoir_controls()
            };
            let mut tank = path
                .reservoir(cx, 0.1, Temperature::new(300.0), 0.0)
                .unwrap();
            for (end, drive, reverse) in [(10.0, 200.0, false), (20.0, -200.0, true)] {
                let initial = tank.clone();
                let (mut lo, mut hi) = controls.reservoir.temperature_bracket_k;
                for _ in 0..80 {
                    let mid = 0.5 * (lo + hi);
                    let residual = 200.0 * (mid - initial.state().temperature_k())
                        - 10.0 * reference(mid, reverse, drive).2;
                    if residual < 0.0 {
                        lo = mid;
                    } else {
                        hi = mid;
                    }
                }
                let answer = path
                    .advance_reservoir_circulation(
                        cx,
                        &mut tank,
                        end,
                        0.0,
                        &[400.0, 300.0],
                        &resistances(),
                        Pressure::new(drive),
                        &closed_gravity(),
                        LiquidFlowSelection::OnlyBracketed,
                        controls,
                    )
                    .unwrap();
                let point = answer.circuit.operating_point();
                close(tank.state().temperature_k(), 0.5 * (lo + hi));
                assert_eq!(answer.circuit.signed_mass_flow_kg_per_s() < 0.0, reverse);
                assert_eq!(
                    point.thermal.regions[0].region,
                    if reverse { "last" } else { "first" }
                );
                assert_eq!(point.thermal.regions[0].inlet, tank.state());
                assert_eq!(answer.initial, initial.state());
                assert_eq!(answer.final_state, tank.state());
                assert_eq!(answer.end_time_s, end);
                close(answer.stored_enthalpy_change_j, answer.returned_heat_j);
                assert!(
                    answer.energy_residual_j.abs()
                        <= controls.reservoir.energy_absolute_tolerance_j
                );
                assert!(
                    answer.total_flow_evaluations
                        >= answer.inlet_evaluations * 2 + point.flow_evaluations
                );
                assert!(answer.total_flow_evaluations <= controls.max_total_flow_evaluations);
            }
            let original = tank.clone();
            // Both zero-pump orientation brackets close, so no branch is silently chosen.
            assert!(matches!(
                path.advance_reservoir_circulation(
                    cx,
                    &mut tank,
                    30.0,
                    0.0,
                    &[400.0, 300.0],
                    &resistances(),
                    Pressure::new(0.0),
                    &closed_gravity(),
                    LiquidFlowSelection::OnlyBracketed,
                    controls
                ),
                Err(RecirculationError::Pressure(
                    PressureDrivenError::AmbiguousCirculationBrackets { .. }
                ))
            ));
            assert_eq!(tank, original);
            assert_eq!(path.inlet().temperature_k(), 300.0);
        },
    );
}
