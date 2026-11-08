//! G1/G3/G4: generated stacks close each original local contact flux.
use super::*;
use fs_conduction::layered::{GeneratedTemperatureContactWallOptions, NonlinearGenerationOptions};

fn generated_options() -> GeneratedTemperatureContactWallOptions {
    let mut solver = options();
    solver.temperature_tolerance_k = 1e-11;
    solver.flux_absolute_tolerance_w_per_m2 = 1e-11;
    solver.flux_relative_tolerance = 1e-13;
    GeneratedTemperatureContactWallOptions {
        generation: NonlinearGenerationOptions {
            temperature_interval_k: [280.0, 500.0],
            solver,
        },
        max_contact_iterations: 128,
        relaxation: 0.5,
    }
}

fn heater_wall(reflected: bool) -> LayeredWall {
    let heater = curve(
        "nonlinear heated solid",
        1.0,
        vec![(280.0, 280.0 / 300.0), (500.0, 500.0 / 300.0)],
    );
    let sink = PlanarThermalLayer::new(
        "unheated solid",
        1.0,
        ConductivityTable::declared(1.0).unwrap(),
    )
    .unwrap();
    LayeredWall::new(
        if reflected {
            vec![sink, heater]
        } else {
            vec![heater, sink]
        },
        vec![1.0],
        [0.0; 2],
    )
    .unwrap()
}

fn heater_contact(sign: f64) -> (InterfaceSystemCard, QueryPoint) {
    let a = (90_000.0 + 1_800.0 * sign).sqrt();
    let b = 300.0 + sign;
    let resistance = (a - b) / sign;
    let mean = 0.5 * (a + b);
    // The actual contact is supported, but the equal 300 K reservoirs are
    // deliberately outside this law. Heating/cooling must not seed equilibrium.
    let (low, high, reference) = if sign > 0.0 {
        (301.0, 320.0, 310.0)
    } else {
        (280.0, 299.0, 290.0)
    };
    (
        contact_card(
            "T",
            vec![
                (low, resistance + 0.01 * (low - mean)),
                (high, resistance + 0.01 * (high - mean)),
            ],
            ValidityDomain::unconstrained().with("T", low, high),
        ),
        mean_point(reference),
    )
}

fn local_binding<'a>(
    card: &'a InterfaceSystemCard,
    point: &QueryPoint,
    index: usize,
) -> LayeredTemperatureContact<'a> {
    LayeredTemperatureContact::new(
        index,
        "original generated-wall joint",
        card,
        point,
        TemperatureContactCoordinates::Mean("T".into()),
        ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
    )
    .unwrap()
}

#[test]
fn g1_generated_equal_reservoirs_use_actual_local_flux_with_both_source_signs_and_reflection() {
    for sign in [1.0, -1.0] {
        let (card, point) = heater_contact(sign);
        for reflected in [false, true] {
            let wall = heater_wall(reflected);
            let original = wall.clone();
            let sources = if reflected {
                [0.0, 8.0 * sign]
            } else {
                [8.0 * sign, 0.0]
            };
            let bindings = [local_binding(&card, &point, 0)];
            let answer = with_cx(|cx| {
                wall.solve_temperature_dependent_with_generation_and_contacts(
                    cx,
                    [300.0; 2],
                    3.0,
                    &sources,
                    &bindings,
                    generated_options(),
                )
            })
            .unwrap();
            // Independent Kirchhoff construction for the heater:
            // q(x)=sign*(-7+8x), T(x)^2=90000+600*sign*(7x-4x^2).
            // The unheated solid has |q|=1; total prescribed heat is sign*24 W.
            let native = &answer.conduction;
            let boundary = if reflected {
                [-sign, 7.0 * sign]
            } else {
                [-7.0 * sign, sign]
            };
            close(native.boundary_fluxes_w_per_m2[0], boundary[0]);
            close(native.boundary_fluxes_w_per_m2[1], boundary[1]);
            close(native.reservoir_heat_rates_w[0], -3.0 * boundary[0]);
            close(native.reservoir_heat_rates_w[1], 3.0 * boundary[1]);
            close(native.generated_heat_rate_w, 24.0 * sign);
            close(native.energy_residual_w, 0.0);
            assert_eq!(native.volumetric_generation_w_m3, sources);
            let heated = usize::from(reflected);
            let row = &native.layer_temperatures[heated];
            let extremum = (90_000.0 + 1_837.5 * sign).sqrt();
            close(
                if sign > 0.0 {
                    row.maximum_temperature_k
                } else {
                    row.minimum_temperature_k
                },
                extremum,
            );
            close(
                if sign > 0.0 {
                    row.maximum_position_m
                } else {
                    row.minimum_position_m
                },
                if reflected { 0.125 } else { 0.875 },
            );
            for x in [0.0, 0.25, 0.5, 0.875, 1.0] {
                let coordinate = if reflected { 1.0 - x } else { x };
                close(
                    with_cx(|cx| native.temperature_at(cx, heated, x)).unwrap(),
                    (90_000.0 + 600.0 * sign * (7.0 * coordinate - 4.0 * coordinate * coordinate))
                        .sqrt(),
                );
            }
            let state = &answer.contacts[0];
            let q = if reflected { -sign } else { sign };
            close(state.solved_flux_w_per_m2, q);
            close(state.physical_flux_w_per_m2, q);
            close(state.continuity_residual_w_per_m2, 0.0);
            assert!(state.flux_residual_w_per_m2.abs() <= 1.1e-11);
            let a = (90_000.0 + 1_800.0 * sign).sqrt();
            let b = 300.0 + sign;
            close(state.temperatures_k[0], if reflected { b } else { a });
            close(state.temperatures_k[1], if reflected { a } else { b });
            close(state.physical.value_m2_k_per_w(), (a - b) / sign);
            assert_eq!(
                state.physical.receipt().selected,
                bindings[0].reference_resistance().receipt().selected
            );
            card.claims()
                .verify_receipt(state.applied.receipt())
                .unwrap();
            card.claims()
                .verify_receipt(state.physical.receipt())
                .unwrap();
            let query = &state.physical.receipt().query_point;
            close(
                query.iter().find(|(axis, _)| axis == "T").unwrap().1,
                0.5 * (a + b),
            );
            assert!(query.contains(&("pressure".into(), 100_000.0)));
            assert!(answer.contact_iterations > 1);
            assert!(answer.total_shooting_iterations >= native.iterations);
            assert_eq!(wall, original);
        }
    }
}

#[test]
fn g1_zero_net_generation_retains_opposite_local_contact_fluxes_and_binding_order() {
    let wall = LayeredWall::new(
        (0..3)
            .map(|i| {
                PlanarThermalLayer::new(
                    format!("layer {i}"),
                    1.0,
                    ConductivityTable::declared(1.0).unwrap(),
                )
                .unwrap()
            })
            .collect(),
        vec![0.4, 0.6],
        [0.0; 2],
    )
    .unwrap();
    let card = contact_card(
        "T",
        vec![(290.0, 0.92), (320.0, 1.22)],
        ValidityDomain::unconstrained().with("T", 290.0, 320.0),
    );
    let bindings = [binding(&card, 1), binding(&card, 0)];
    let answer = with_cx(|cx| {
        wall.solve_temperature_dependent_with_generation_and_contacts(
            cx,
            [300.0; 2],
            2.0,
            &[2.0, -4.0, 2.0],
            &bindings,
            generated_options(),
        )
    })
    .unwrap();
    // Independent exact field: faces [300,299], [297,297], [299,300].
    // q at the two contacts is +2 and -2 even though both boundary fluxes and
    // net generation vanish. The middle sink has a 296.5 K interior minimum.
    for (row, expected) in answer.conduction.layer_temperatures.iter().zip([
        [300.0, 299.0],
        [297.0, 297.0],
        [299.0, 300.0],
    ]) {
        close(row.face_temperatures_k[0], expected[0]);
        close(row.face_temperatures_k[1], expected[1]);
    }
    for (state, (index, flux)) in answer.contacts.iter().zip([(1, -2.0), (0, 2.0)]) {
        assert_eq!(state.interface_index, index);
        close(state.solved_flux_w_per_m2, flux);
        close(state.physical_flux_w_per_m2, flux);
        close(state.physical.value_m2_k_per_w(), 1.0);
        close(state.continuity_residual_w_per_m2, 0.0);
    }
    close(
        answer.conduction.layer_temperatures[1].minimum_temperature_k,
        296.5,
    );
    close(
        answer.conduction.layer_temperatures[1].minimum_position_m,
        0.5,
    );
    for value in answer.conduction.boundary_fluxes_w_per_m2 {
        close(value, 0.0);
    }
    close(answer.conduction.generated_heat_rate_w, 0.0);
    close(answer.conduction.energy_residual_w, 0.0);
    assert!(answer.contact_iterations > 1);
}

#[test]
fn g3_generated_contact_path_preserves_zero_source_and_fixed_source_bulk_results() {
    let options = generated_options();
    let wall = wall();
    let card = mean_card(300.0, 450.0);
    let bindings = [binding(&card, 0)];
    for reservoirs in [[355.0; 2], [400.0, 300.0], [300.0, 400.0]] {
        let ordinary = with_cx(|cx| {
            wall.solve_temperature_dependent_with_contacts(
                cx,
                reservoirs,
                2.0,
                &bindings,
                TemperatureContactWallOptions {
                    solver: options.generation.solver,
                    max_contact_iterations: options.max_contact_iterations,
                    relaxation: options.relaxation,
                },
            )
        })
        .unwrap();
        let generated = with_cx(|cx| {
            wall.solve_temperature_dependent_with_generation_and_contacts(
                cx, reservoirs, 2.0, &[0.0; 3], &bindings, options,
            )
        })
        .unwrap();
        for flux in generated.conduction.boundary_fluxes_w_per_m2 {
            close(flux, ordinary.conduction.heat_flux_w_per_m2);
        }
        close(generated.conduction.generated_heat_rate_w, 0.0);
        close(
            generated.contacts[0].physical.value_m2_k_per_w(),
            ordinary.contacts[0].physical.value_m2_k_per_w(),
        );
        if reservoirs[0] == reservoirs[1] {
            assert_eq!(generated.contact_iterations, 1);
            close(generated.contacts[0].applied.value_m2_k_per_w(), 1.3);
            assert!(
                generated.contacts[0]
                    .physical
                    .receipt()
                    .query_point
                    .contains(&("T".into(), 355.0))
            );
        }
    }
    let wall = LayeredWall::new(wall.layers().to_vec(), vec![0.4, 0.5], [0.1, 0.2]).unwrap();
    let constant = contact_card(
        "T",
        vec![(300.0, 0.4), (450.0, 0.4)],
        ValidityDomain::unconstrained().with("T", 300.0, 450.0),
    );
    let expected = with_cx(|cx| {
        wall.solve_temperature_dependent_with_generation(
            cx,
            [380.0, 300.0],
            2.0,
            &[2.0, -1.0, 3.0],
            options.generation,
        )
    })
    .unwrap();
    for bindings in [vec![], vec![binding(&constant, 0)]] {
        let answer = with_cx(|cx| {
            wall.solve_temperature_dependent_with_generation_and_contacts(
                cx,
                [380.0, 300.0],
                2.0,
                &[2.0, -1.0, 3.0],
                &bindings,
                options,
            )
        })
        .unwrap();
        assert_eq!(answer.conduction, expected);
        assert_eq!(answer.contact_iterations, 1);
    }
}

#[test]
fn g0_g4_generated_contacts_preserve_source_extremum_budgets_and_cancellation_refusals() {
    let wall = heater_wall(false);
    let (card, point) = heater_contact(1.0);
    let bindings = [local_binding(&card, &point, 0)];
    for sources in [vec![], vec![8.0], vec![8.0, f64::NAN]] {
        assert!(
            with_cx(
                |cx| wall.solve_temperature_dependent_with_generation_and_contacts(
                    cx,
                    [300.0; 2],
                    1.0,
                    &sources,
                    &bindings,
                    generated_options(),
                )
            )
            .is_err()
        );
    }
    let mut budget = generated_options();
    budget.max_contact_iterations = 1;
    assert!(
        with_cx(
            |cx| wall.solve_temperature_dependent_with_generation_and_contacts(
                cx,
                [300.0; 2],
                1.0,
                &[8.0, 0.0],
                &bindings,
                budget,
            )
        )
        .is_err()
    );
    let mut native_budget = generated_options();
    native_budget.generation.solver.max_iterations = 1;
    assert!(
        with_cx(
            |cx| wall.solve_temperature_dependent_with_generation_and_contacts(
                cx,
                [300.0; 2],
                1.0,
                &[8.0, 0.0],
                &bindings,
                native_budget,
            )
        )
        .is_err()
    );
    let mut narrow = generated_options();
    // Both exact faces fit, but the original heated-layer maximum does not.
    narrow.generation.temperature_interval_k[1] = 302.99;
    assert!(
        with_cx(
            |cx| wall.solve_temperature_dependent_with_generation_and_contacts(
                cx,
                [300.0; 2],
                1.0,
                &[8.0, 0.0],
                &bindings,
                narrow,
            )
        )
        .is_err()
    );
    let narrow_card = mean_card(309.0, 311.0);
    assert!(
        with_cx(
            |cx| wall.solve_temperature_dependent_with_generation_and_contacts(
                cx,
                [300.0; 2],
                1.0,
                &[8.0, 0.0],
                &[binding(&narrow_card, 0)],
                generated_options(),
            )
        )
        .is_err()
    );
    assert!(matches!(
        with_cancelled_cx(
            |cx| wall.solve_temperature_dependent_with_generation_and_contacts(
                cx,
                [300.0; 2],
                1.0,
                &[8.0, 0.0],
                &bindings,
                generated_options(),
            )
        ),
        Err(ConductionError::Cancelled { .. })
    ));
}
