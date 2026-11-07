//! G1/G3 finite gas with original coated bodies and independent physical balances.
use super::*;
use crate::atmosphere::{atmosphere as gas, carrier_mass, controls as gas_controls};
use fs_couple::moisture::{
    atmosphere::AtmosphereMoistureError,
    radial::layered::atmosphere::LayeredRadialMoistureParticipant,
};

const RHO: [f64; 2] = [300.0, 900.0];
const SLOPE: [f64; 2] = [0.2, 0.4];
const DIFFUSION: [f64; 2] = [0.001, 0.0005];
const INTERFACE_H: f64 = 0.03;
const EXTERIOR_H: f64 = 0.04;

fn coated(
    symmetry: RadialMoistureSymmetry,
    activity: [f64; 2],
    cx: &Cx<'_>,
) -> MoistureLayeredRadial {
    MoistureLayeredRadial::new(
        (0..2)
            .map(|i| {
                radial_layer(
                    linear_card(RHO[i], SLOPE[i], DIFFUSION[i]),
                    symmetry,
                    if i == 0 { vec![0.0, 0.03] } else { vec![0.03, 0.05] },
                    vec![activity[i]],
                    cx,
                )
            })
            .collect(),
        2,
        cx,
    )
    .unwrap()
    .with_transfer(
        [
            MoistureActivityTransfer::IdealContact,
            MoistureActivityTransfer::ConductanceKgM2S(EXTERIOR_H),
        ],
        vec![MoistureActivityTransfer::ConductanceKgM2S(INTERFACE_H)],
        cx,
    )
    .unwrap()
}

fn participant(id: u64, moisture: &MoistureLayeredRadial) -> LayeredRadialMoistureParticipant<'_> {
    LayeredRadialMoistureParticipant {
        id,
        moisture,
        exposed: [false, true],
        transport: config(),
        nonlinear: nonlinear(),
    }
}

fn cell_water(before: &MoistureRadial, after: &MoistureRadial) -> f64 {
    before
        .moisture()
        .iter()
        .zip(after.moisture())
        .zip(before.dry_masses_kg())
        .map(|((&old, &new), &dry)| dry * (new - old))
        .sum()
}

fn actual_change(before: &MoistureLayeredRadial, after: &MoistureLayeredRadial) -> f64 {
    before
        .layers()
        .iter()
        .zip(after.layers())
        .map(|(old, new)| cell_water(old.radial(), new.radial()))
        .sum()
}

/// Independent physical storage and constant-conductivity curved resistances.
/// Returns C_i=B_i*s_i, core/coating half paths, complete interface and exterior.
fn coefficients(symmetry: RadialMoistureSymmetry) -> ([f64; 2], [f64; 2], f64, f64, f64) {
    let capacity = [
        RHO[0] * measures(symmetry, 0.0, 0.03).0 * SLOPE[0],
        RHO[1] * measures(symmetry, 0.03, 0.05).0 * SLOPE[1],
    ];
    let half = [
        RHO[0] * DIFFUSION[0] * SLOPE[0] / measures(symmetry, 0.015, 0.03).1,
        RHO[1] * DIFFUSION[1] * SLOPE[1] / measures(symmetry, 0.03, 0.04).1,
    ];
    let join = 1.0
        / (1.0 / half[0]
            + 1.0 / (area(symmetry, 0.03) * INTERFACE_H)
            + 1.0 / half[1]);
    let outer_bulk =
        RHO[1] * DIFFUSION[1] * SLOPE[1] / measures(symmetry, 0.04, 0.05).1;
    let outer =
        1.0 / (1.0 / outer_bulk + 1.0 / (area(symmetry, 0.05) * EXTERIOR_H));
    (capacity, half, join, outer, outer_bulk)
}

/// Eliminate two implicit cell equations to obtain a_i'=u_i+v_i*RH'.
fn response(
    symmetry: RadialMoistureSymmetry,
    initial: [f64; 2],
    dt: f64,
) -> ([f64; 2], [f64; 2]) {
    let (c, _, join, outer, _) = coefficients(symmetry);
    let a = c[0] + dt * join;
    let b = dt * join;
    let d = c[1] + dt * (join + outer);
    let determinant = a * d - b * b;
    let stored = [c[0] * initial[0], c[1] * initial[1]];
    (
        [
            (d * stored[0] + b * stored[1]) / determinant,
            (b * stored[0] + a * stored[1]) / determinant,
        ],
        [
            b * dt * outer / determinant,
            a * dt * outer / determinant,
        ],
    )
}

#[test]
fn finite_coated_radial_matches_analytic_absorption_drying_and_continuation() {
    with_cx(|cx| {
        for symmetry in [
            RadialMoistureSymmetry::Sphere,
            RadialMoistureSymmetry::Cylinder { length_m: 0.1 },
        ] {
            for (initial_rh, initial) in [(0.8, [0.1, 0.2]), (0.15, [0.8, 0.6])] {
                let mut atmosphere = gas(initial_rh, 0.5);
                let mut body = coated(symmetry, initial, cx);
                let gas_capacity = atmosphere.vapor_mass_kg().unwrap() / initial_rh;
                let (capacity, halves, join, outer, outer_bulk) = coefficients(symmetry);
                let carrier = carrier_mass(&atmosphere);
                for interval in 0..3 {
                    let dt = 2.0;
                    let initial: [f64; 2] = body.activity().try_into().unwrap();
                    let rh = atmosphere.relative_humidity().unwrap();
                    let (u, v) = response(symmetry, initial, dt);
                    let expected_rh = (gas_capacity * rh
                        + (0..2)
                            .map(|i| capacity[i] * (initial[i] - u[i]))
                            .sum::<f64>())
                        / (gas_capacity
                            + (0..2).map(|i| capacity[i] * v[i]).sum::<f64>());
                    let expected: [f64; 2] =
                        std::array::from_fn(|i| u[i] + v[i] * expected_rh);
                    let exchange = dt * outer * (expected_rh - expected[1]);
                    let internal = dt * join * (expected[0] - expected[1]);
                    let surface = expected[1] + exchange / (dt * outer_bulk);
                    let total_water =
                        atmosphere.vapor_mass_kg().unwrap() + body.water_mass_kg();
                    let total_mass = atmosphere.gas().mass_kg() + body.wet_mass_kg();
                    let candidate = atmosphere
                        .advance_layered_radial(
                            &body,
                            dt,
                            [false, true],
                            config(),
                            nonlinear(),
                            gas_controls(),
                            cx,
                        )
                        .unwrap();
                    let next = candidate.moisture();
                    let report = candidate.step();
                    near(
                        candidate.atmosphere().relative_humidity().unwrap(),
                        expected_rh,
                        4e-11,
                    );
                    for i in 0..2 {
                        near(next.activity()[i], expected[i], 4e-12);
                        near(
                            next.layers()[i].radial().moisture()[0],
                            SLOPE[i] * expected[i],
                            3e-12,
                        );
                    }
                    near(report.water_to_solid_kg, exchange, 2e-12);
                    near(report.transport.boundary_water_kg[1], exchange, 2e-12);
                    assert_eq!(report.transport.boundary_water_kg[0], 0.0);
                    near(
                        report.transport.interfaces[0].water_left_to_right_kg,
                        internal,
                        2e-12,
                    );
                    near(
                        report.transport.interfaces[0].relative_activity,
                        expected[0] - internal / (dt * halves[0]),
                        4e-12,
                    );
                    near(
                        report.transport.interfaces[0].right_relative_activity,
                        expected[1] + internal / (dt * halves[1]),
                        4e-12,
                    );
                    near(
                        cell_water(body.layers()[0].radial(), next.layers()[0].radial()),
                        -internal,
                        2e-12,
                    );
                    near(
                        cell_water(body.layers()[1].radial(), next.layers()[1].radial()),
                        internal + exchange,
                        2e-12,
                    );
                    near(actual_change(&body, next), exchange, 2e-12);
                    near(
                        report.transport.surface_activities[1].unwrap(),
                        surface,
                        4e-12,
                    );
                    let source = report.sorption[1].as_ref().unwrap();
                    near(source.query().axes()["relative_humidity"], surface, 4e-12);
                    assert!((surface - expected_rh).abs() > 1e-3);
                    assert_eq!(
                        source.property("sorption").unwrap().answer().receipt.selected,
                        body.layers()[1].sorption().selected_claim(),
                    );
                    assert!(report.sorption[0].is_none());
                    near(
                        candidate.atmosphere().vapor_mass_kg().unwrap()
                            + next.water_mass_kg(),
                        total_water,
                        3e-12,
                    );
                    near(
                        candidate.atmosphere().gas().mass_kg() + next.wet_mass_kg(),
                        total_mass,
                        3e-12,
                    );
                    near(carrier_mass(candidate.atmosphere()), carrier, 3e-13);
                    near(
                        report.water_balance_residual_kg,
                        0.0,
                        gas_controls().mass_tolerance_kg,
                    );
                    assert_eq!(
                        report.transport.interval_s,
                        [f64::from(interval) * dt, f64::from(interval + 1) * dt],
                    );
                    assert_eq!(exchange.is_sign_positive(), initial_rh == 0.8);
                    let (next_air, next_body, _) = candidate.into_parts();
                    atmosphere = next_air;
                    body = next_body;
                }
            }
        }
    });
}

#[test]
fn finite_coated_radial_shared_inventory_matches_coupled_solution_and_canonical_order() {
    with_cx(|cx| {
        let symmetries = [
            RadialMoistureSymmetry::Sphere,
            RadialMoistureSymmetry::Cylinder { length_m: 0.1 },
        ];
        let initial = [[0.1, 0.2], [0.8, 0.7]];
        let dry = coated(symmetries[0], initial[0], cx);
        let wet = coated(symmetries[1], initial[1], cx);
        let atmosphere = gas(0.4, 1.0);
        let capacity = atmosphere.vapor_mass_kg().unwrap() / 0.4;
        let dt = 2.0;
        let response = [
            response(symmetries[0], initial[0], dt),
            response(symmetries[1], initial[1], dt),
        ];
        let storage = [coefficients(symmetries[0]).0, coefficients(symmetries[1]).0];
        let mut numerator = capacity * 0.4;
        let mut denominator = capacity;
        for body in 0..2 {
            for cell in 0..2 {
                numerator += storage[body][cell]
                    * (initial[body][cell] - response[body].0[cell]);
                denominator += storage[body][cell] * response[body].1[cell];
            }
        }
        let expected_rh = numerator / denominator;
        let forward = atmosphere
            .advance_shared_layered_radial(
                &[participant(7, &dry), participant(2, &wet)],
                dt,
                4,
                gas_controls(),
                cx,
            )
            .unwrap();
        let reverse = atmosphere
            .advance_shared_layered_radial(
                &[participant(2, &wet), participant(7, &dry)],
                dt,
                4,
                gas_controls(),
                cx,
            )
            .unwrap();
        assert_eq!(forward.atmosphere(), reverse.atmosphere());
        assert_eq!(
            forward.solids().iter().map(|b| b.id).collect::<Vec<_>>(),
            vec![2, 7],
        );
        near(
            forward.atmosphere().relative_humidity().unwrap(),
            expected_rh,
            4e-11,
        );
        let mut change = 0.0;
        for (a, b) in forward.solids().iter().zip(reverse.solids()) {
            assert_eq!(a.moisture.activity(), b.moisture.activity());
            let index = if a.id == 7 { 0 } else { 1 };
            for cell in 0..2 {
                near(
                    a.moisture.activity()[cell],
                    response[index].0[cell] + response[index].1[cell] * expected_rh,
                    4e-12,
                );
            }
            let original = if a.id == 7 { &dry } else { &wet };
            let local = actual_change(original, &a.moisture);
            assert_eq!(local.is_sign_positive(), a.id == 7);
            near(local, a.transport.boundary_water_kg[1], 2e-12);
            change += local;
            assert_eq!(a.moisture.time_s(), dt);
        }
        near(forward.step().water_to_solids_kg, change, 2e-12);
        near(
            atmosphere.vapor_mass_kg().unwrap()
                - forward.atmosphere().vapor_mass_kg().unwrap(),
            change,
            2e-12,
        );
        near(forward.step().total_mass_residual_kg, 0.0, 1e-13);
        near(
            carrier_mass(forward.atmosphere()),
            carrier_mass(&atmosphere),
            3e-13,
        );
        assert!(matches!(
            atmosphere.advance_shared_layered_radial(
                &[participant(7, &dry), participant(7, &wet)],
                dt, 4, gas_controls(), cx,
            ),
            Err(AtmosphereMoistureError::Invalid(_)),
        ));
        assert!(matches!(
            atmosphere.advance_shared_layered_radial(
                &[participant(7, &dry), participant(2, &wet)],
                dt, 3, gas_controls(), cx,
            ),
            Err(AtmosphereMoistureError::Invalid(_)),
        ));
        assert!(matches!(
            atmosphere.advance_shared_layered_radial(
                &[participant(7, &dry), participant(2, &forward.solids()[0].moisture)],
                dt, 4, gas_controls(), cx,
            ),
            Err(AtmosphereMoistureError::Invalid(_)),
        ));
        assert_eq!(dry.time_s(), 0.0);
        assert_eq!(wet.time_s(), 0.0);
    });
}

#[test]
fn finite_coated_radial_wets_degenerate_materials_with_exact_local_and_gas_accounts() {
    with_cx(|cx| {
        let body = MoistureLayeredRadial::new(
            vec![
                radial_layer(
                    material(
                        500.0,
                        vec![(0.0, 0.0), (1.0, 0.25)],
                        [0.0, 0.006],
                    ),
                    RadialMoistureSymmetry::Sphere,
                    vec![0.0, 0.03],
                    vec![0.0],
                    cx,
                ),
                radial_layer(
                    material(
                        750.0,
                        vec![(0.0, 0.0), (1.0, 0.5)],
                        [0.0, 0.015],
                    ),
                    RadialMoistureSymmetry::Sphere,
                    vec![0.03, 0.05],
                    vec![0.0],
                    cx,
                ),
            ],
            2,
            cx,
        )
        .unwrap();
        let atmosphere = gas(0.8, 0.2);
        let saved = atmosphere.clone();
        let accepted = atmosphere
            .advance_layered_radial(
                &body, 0.01, [false, true], config(), nonlinear(), gas_controls(), cx,
            )
            .unwrap();
        let next = accepted.moisture();
        let report = &accepted.step().transport;
        assert!(next.activity()[0] > 0.0 && next.activity()[1] > next.activity()[0]);
        assert!(accepted.atmosphere().relative_humidity().unwrap() > 0.0);
        assert!(accepted.atmosphere().relative_humidity().unwrap() < 0.8);
        let internal = report.interfaces[0].water_left_to_right_kg;
        assert!(internal < 0.0);
        near(next.layers()[0].radial().water_mass_kg(), -internal, 3e-12);
        near(
            next.layers()[1].radial().water_mass_kg(),
            internal + report.boundary_water_kg[1],
            3e-12,
        );
        near(next.water_mass_kg(), report.boundary_water_kg[1], 3e-12);
        near(
            accepted.atmosphere().vapor_mass_kg().unwrap() + next.water_mass_kg(),
            atmosphere.vapor_mass_kg().unwrap(),
            3e-12,
        );
        assert_eq!(
            report.sorption_sources,
            body.layers().iter().map(|l| l.sorption().identity()).collect::<Vec<_>>(),
        );
        assert_eq!(atmosphere, saved);
        assert_eq!(contents(&body), vec![0.0, 0.0]);
        assert_eq!(body.time_s(), 0.0);
    });
}

fn restricted_body(cx: &Cx<'_>) -> MoistureLayeredRadial {
    let layers = (0..2)
        .map(|i| {
            let card = linear_card(RHO[i], SLOPE[i], DIFFUSION[i]);
            let full = radial_layer(
                card.clone(),
                RadialMoistureSymmetry::Sphere,
                if i == 0 { vec![0.0, 0.03] } else { vec![0.03, 0.05] },
                vec![if i == 0 { 0.3 } else { 0.6 }],
                cx,
            );
            let source = resolve_material_state_point(
                &card,
                full.sorption().source_point(),
                &[ScalarPropertyRequirement::try_with_quantity(
                    "sorption",
                    qty(QuantityKind::Moisture(MoistureBasis::DryMass)),
                    ScalarAdmissibility::NonNegative,
                )
                .unwrap()],
                MaterialPropertySelection::SingleClaimOnly,
            )
            .unwrap();
            let limited =
                MonotoneMoistureSorption::from_binding(&card, &source, "sorption", [0.2, 0.8])
                    .unwrap();
            MoistureRadialLayer::new(full.radial().clone(), limited, cx).unwrap()
        })
        .collect();
    MoistureLayeredRadial::new(layers, 2, cx).unwrap()
}

#[test]
fn finite_coated_radial_respects_source_support_sealed_surfaces_and_atomic_refusals() {
    with_cx(|cx| {
        let atmosphere = gas(0.5, 0.2);
        let saved = atmosphere.clone();
        let restricted = restricted_body(cx);
        assert!(matches!(
            atmosphere.advance_layered_radial(
                &restricted, 0.2, [false, true], config(), nonlinear(), gas_controls(), cx,
            ),
            Err(AtmosphereMoistureError::Invalid(_)),
        ));
        let closed = restricted
            .clone()
            .with_transfer(
                [
                    MoistureActivityTransfer::IdealContact,
                    MoistureActivityTransfer::Sealed,
                ],
                vec![MoistureActivityTransfer::IdealContact],
                cx,
            )
            .unwrap();
        let humid = gas(0.9, 0.2); // Outside the body's chart, but physically sealed.
        let accepted = humid
            .advance_layered_radial(
                &closed, 0.2, [false, true], config(), nonlinear(), gas_controls(), cx,
            )
            .unwrap();
        assert_eq!(accepted.atmosphere(), &humid);
        assert!(accepted.moisture().activity()[0] > closed.activity()[0]);
        assert!(accepted.moisture().activity()[1] < closed.activity()[1]);
        assert!(accepted.step().sorption.iter().all(Option::is_none));
        near(accepted.step().water_to_solid_kg, 0.0, 1e-13);
        let body = coated(RadialMoistureSymmetry::Sphere, [0.0, 0.0], cx);
        assert!(matches!(
            atmosphere.advance_layered_radial(
                &body, 0.2, [true, true], config(), nonlinear(), gas_controls(), cx,
            ),
            Err(AtmosphereMoistureError::Invalid(_)),
        ));
        let mut limits = gas_controls();
        limits.maximum_iterations = 1;
        assert!(matches!(
            atmosphere.advance_layered_radial(
                &body, 2.0, [false, true], config(), nonlinear(), limits, cx,
            ),
            Err(AtmosphereMoistureError::NotConverged { .. }),
        ));
        assert_eq!(atmosphere, saved);
        assert_eq!(contents(&body), vec![0.0, 0.0]);
        assert_eq!(body.time_s(), 0.0);

        let card = linear_card(500.0, 0.25, 0.001);
        let point_sorption = isotherm(&card)
            .source_point()
            .clone()
            .with_quantity("pressure", qty(QuantityKind::Pressure), 101_325.0)
            .unwrap();
        let sorption = isotherm_at(&card, point_sorption);
        let field = MoistureRadial::new(
            card,
            point(0.1)
                .with_quantity("pressure", qty(QuantityKind::Pressure), 101_325.0)
                .unwrap(),
            "moisture_diffusivity",
            MaterialPropertySelection::SingleClaimOnly,
            MoistureRadialGeometry {
                symmetry: RadialMoistureSymmetry::Sphere,
                faces_m: vec![0.0, 0.05],
                maximum_moisture: 0.6,
                maximum_cells: 1,
            },
            vec![0.1],
            cx,
        )
        .unwrap();
        let held_pressure = MoistureLayeredRadial::new(
            vec![MoistureRadialLayer::new(field, sorption, cx).unwrap()],
            1,
            cx,
        )
        .unwrap();
        assert!(matches!(
            atmosphere.advance_layered_radial(
                &held_pressure, 0.2, [false, true], config(), nonlinear(), gas_controls(), cx,
            ),
            Err(AtmosphereMoistureError::Moisture(MoistureError::Invalid(_))),
        ));
        assert_eq!(atmosphere, saved);
        assert_eq!(contents(&held_pressure), vec![0.1]);
    });
}
