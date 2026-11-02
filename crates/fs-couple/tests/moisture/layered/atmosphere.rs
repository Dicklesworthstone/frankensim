//! Finite native vapor, material partition and complete stack inventories.
use super::super::atmosphere as air;
use super::*;
use fs_couple::moisture::atmosphere::AtmosphereMoistureError;
use fs_couple::moisture::layered::atmosphere::LayeredMoistureParticipant;

const CAPACITY: [f64; 2] = [0.032, 0.012];
const SORPTION: [f64; 2] = [0.2, 0.4];

fn stack(activity: [f64; 2], cx: &Cx<'_>) -> MoistureLayeredSlab {
    MoistureLayeredSlab::new(
        vec![
            layer(
                material(800.0, vec![(0.0, 0.0), (1.0, SORPTION[0])], [0.01; 2]),
                0.2,
                vec![SORPTION[0] * activity[0]],
                cx,
            ),
            layer(
                material(300.0, vec![(0.0, 0.0), (1.0, SORPTION[1])], [0.02; 2]),
                0.1,
                vec![SORPTION[1] * activity[1]],
                cx,
            ),
        ],
        256,
        cx,
    )
    .unwrap()
}

/// Independent two-equation elimination: a_k = u_k + v_k*RH_end.
fn response(initial: [f64; 2], dt: f64) -> ([f64; 2], [f64; 2]) {
    let half = [0.016, 0.048];
    let join = 0.012;
    let a = CAPACITY[0] + dt * (half[0] + join);
    let d = CAPACITY[1] + dt * (half[1] + join);
    let b = dt * join;
    let determinant = a * d - b * b;
    let stored = [CAPACITY[0] * initial[0], CAPACITY[1] * initial[1]];
    (
        [
            (d * stored[0] + b * stored[1]) / determinant,
            (b * stored[0] + a * stored[1]) / determinant,
        ],
        [
            dt * (d * half[0] + b * half[1]) / determinant,
            dt * (b * half[0] + a * half[1]) / determinant,
        ],
    )
}

#[test]
fn g3_finite_vapor_and_two_materials_match_coupled_adsorption_desorption_and_continuation() {
    with_cx(|cx| {
        for (humidity, initial) in [(0.8, [0.2, 0.1]), (0.1, [0.8, 0.5])] {
            let mut atmosphere = air::atmosphere(humidity, 0.1);
            let mut moisture = stack(initial, cx);
            let capacity = atmosphere.vapor_mass_kg().unwrap() / humidity;
            let carrier = air::carrier_mass(&atmosphere);
            for _ in 0..2 {
                let dt = 0.5;
                let initial: [f64; 2] = moisture.activity().try_into().unwrap();
                let rh0 = atmosphere.relative_humidity().unwrap();
                let (u, v) = response(initial, dt);
                let expected_rh = (capacity * rh0
                    + (0..2)
                        .map(|k| CAPACITY[k] * (initial[k] - u[k]))
                        .sum::<f64>())
                    / (capacity + (0..2).map(|k| CAPACITY[k] * v[k]).sum::<f64>());
                let old_water = atmosphere.vapor_mass_kg().unwrap() + moisture.water_mass_kg();
                let old_mass = atmosphere.gas().mass_kg() + moisture.wet_mass_kg();
                let candidate = atmosphere
                    .advance_layered(
                        &moisture,
                        dt,
                        [true; 2],
                        config(),
                        nonlinear(),
                        air::controls(),
                        cx,
                    )
                    .unwrap();
                near(
                    candidate.atmosphere().relative_humidity().unwrap(),
                    expected_rh,
                    3e-10,
                );
                for k in 0..2 {
                    near(
                        candidate.moisture().activity()[k],
                        u[k] + v[k] * expected_rh,
                        3e-10,
                    );
                    near(
                        candidate.step().sorption[k]
                            .as_ref()
                            .unwrap()
                            .property("sorption")
                            .unwrap()
                            .value_si(),
                        SORPTION[k] * expected_rh,
                        3e-10,
                    );
                    assert_eq!(
                        candidate.step().sorption[k]
                            .as_ref()
                            .unwrap()
                            .card_identity(),
                        moisture.layers()[k].slab().card().content_hash()
                    );
                }
                near(
                    candidate.atmosphere().vapor_mass_kg().unwrap()
                        + candidate.moisture().water_mass_kg(),
                    old_water,
                    2e-13,
                );
                near(
                    candidate.atmosphere().gas().mass_kg() + candidate.moisture().wet_mass_kg(),
                    old_mass,
                    2e-13,
                );
                near(air::carrier_mass(candidate.atmosphere()), carrier, 2e-13);
                near(
                    candidate.step().water_to_solid_kg,
                    candidate.moisture().water_mass_kg() - moisture.water_mass_kg(),
                    2e-13,
                );
                near(
                    candidate.step().water_to_solid_kg,
                    candidate
                        .step()
                        .transport
                        .boundary_water_kg
                        .iter()
                        .sum::<f64>(),
                    2e-13,
                );
                assert_eq!(
                    candidate.moisture().layers()[0].slab().time_s(),
                    moisture.time_s() + dt
                );
                let (next_atmosphere, next_moisture, _) = candidate.into_parts();
                atmosphere = next_atmosphere;
                moisture = next_moisture;
            }
        }
    });
}

#[test]
fn g1_common_activity_and_sealed_internal_redistribution_retain_the_exact_native_gas() {
    with_cx(|cx| {
        let atmosphere = air::atmosphere(0.5, 0.1);
        let equilibrium = stack([atmosphere.relative_humidity().unwrap(); 2], cx);
        let accepted = atmosphere
            .advance_layered(
                &equilibrium,
                1.0,
                [true; 2],
                config(),
                nonlinear(),
                air::controls(),
                cx,
            )
            .unwrap();
        assert_eq!(accepted.atmosphere(), &atmosphere);
        assert_eq!(accepted.moisture().activity(), equilibrium.activity());
        assert_eq!(accepted.step().water_to_solid_kg, 0.0);
        assert_eq!(accepted.step().iterations, 0);

        let gradient = stack([0.8, 0.1], cx);
        let accepted = atmosphere
            .advance_layered(
                &gradient,
                1.0,
                [false; 2],
                config(),
                nonlinear(),
                air::controls(),
                cx,
            )
            .unwrap();
        assert_eq!(accepted.atmosphere(), &atmosphere);
        assert!(accepted.moisture().activity()[0] < gradient.activity()[0]);
        assert!(accepted.moisture().activity()[1] > gradient.activity()[1]);
        near(
            accepted.moisture().water_mass_kg(),
            gradient.water_mass_kg(),
            2e-13,
        );
        assert!(accepted.step().sorption.iter().all(Option::is_none));
    });
}

#[test]
fn g3_shared_vapor_balances_opposing_layered_exchanges_independent_of_participant_order() {
    with_cx(|cx| {
        let atmosphere = air::atmosphere(0.5, 0.2);
        let wet = stack([0.8; 2], cx);
        let dry = stack([0.1; 2], cx);
        let participant = |id, moisture| LayeredMoistureParticipant {
            id,
            moisture,
            exposed: [true; 2],
            transport: config(),
            nonlinear: nonlinear(),
        };
        let candidate = atmosphere
            .advance_shared_layered(
                &[participant(9, &wet), participant(2, &dry)],
                0.5,
                4,
                air::controls(),
                cx,
            )
            .unwrap();
        let reversed = atmosphere
            .advance_shared_layered(
                &[participant(2, &dry), participant(9, &wet)],
                0.5,
                4,
                air::controls(),
                cx,
            )
            .unwrap();
        assert_eq!(candidate.atmosphere(), reversed.atmosphere());
        assert_eq!(
            candidate.solids().iter().map(|s| s.id).collect::<Vec<_>>(),
            [2, 9]
        );
        for (left, right) in candidate.solids().iter().zip(reversed.solids()) {
            assert_eq!(left.moisture.activity(), right.moisture.activity());
        }
        assert!(candidate.solids()[0].moisture.water_mass_kg() > dry.water_mass_kg());
        assert!(candidate.solids()[1].moisture.water_mass_kg() < wet.water_mass_kg());
        near(
            candidate.atmosphere().vapor_mass_kg().unwrap()
                + candidate
                    .solids()
                    .iter()
                    .map(|s| s.moisture.water_mass_kg())
                    .sum::<f64>(),
            atmosphere.vapor_mass_kg().unwrap() + wet.water_mass_kg() + dry.water_mass_kg(),
            2e-13,
        );
        near(
            air::carrier_mass(candidate.atmosphere()),
            air::carrier_mass(&atmosphere),
            2e-13,
        );
        near(candidate.step().total_mass_residual_kg, 0.0, 2e-13);
    });
}

#[test]
fn g0_finite_layered_exchange_refuses_insufficient_work_and_changing_unowned_source_context() {
    with_cx(|cx| {
        let atmosphere = air::atmosphere(0.8, 0.1);
        let saved = atmosphere.clone();
        let moisture = stack([0.0; 2], cx);
        let mut controls = air::controls();
        controls.maximum_iterations = 1;
        assert!(matches!(
            atmosphere.advance_layered(
                &moisture,
                0.5,
                [true; 2],
                config(),
                nonlinear(),
                controls,
                cx
            ),
            Err(AtmosphereMoistureError::NotConverged { .. })
        ));
        assert_eq!(atmosphere, saved);
        assert_eq!(moisture.activity(), &[0.0; 2]);

        let card = material(500.0, vec![(0.0, 0.0), (1.0, 0.2)], [0.01; 2]);
        let sorption_point = isotherm(&card)
            .source_point()
            .clone()
            .with_quantity("pressure", qty(QuantityKind::Pressure), 101_325.0)
            .unwrap();
        let sorption = isotherm_at(&card, sorption_point);
        let slab = MoistureSlab::new_kirchhoff(
            card,
            point(0.1)
                .with_quantity("pressure", qty(QuantityKind::Pressure), 101_325.0)
                .unwrap(),
            "moisture_diffusivity",
            MaterialPropertySelection::SingleClaimOnly,
            MoistureSlabGeometry {
                length_m: 0.2,
                cross_section_m2: 0.001,
                maximum_moisture: 0.6,
                maximum_cells: 256,
            },
            vec![0.1],
            nonlinear(),
            cx,
        )
        .unwrap();
        let held_pressure = MoistureLayeredSlab::new(
            vec![MoistureLayer::new(slab, sorption, cx).unwrap()],
            256,
            cx,
        )
        .unwrap();
        assert!(matches!(
            atmosphere.advance_layered(
                &held_pressure,
                0.5,
                [true; 2],
                config(),
                nonlinear(),
                air::controls(),
                cx
            ),
            Err(AtmosphereMoistureError::Moisture(MoistureError::Invalid(_)))
        ));
        assert_eq!(atmosphere, saved);
        assert_eq!(held_pressure.layers()[0].slab().moisture(), &[0.1]);
    });
}
