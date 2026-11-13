//! G1/G3: independent coupled inventories across unlike native material layers.
#[path = "pair/adaptive.rs"]
mod adaptive;

use super::*;
use fs_couple::moisture::atmosphere::pair::{
    MoistureAtmospherePort, PairedAtmosphereMoistureSolve,
};

fn controls() -> PairedAtmosphereMoistureSolve {
    PairedAtmosphereMoistureSolve {
        atmospheres: [air::controls(); 2],
        maximum_evaluations: 10_000,
        mass_tolerance_kg: 5e-12,
    }
}

/// Eliminate both finite vapor inventories from the four backward-Euler rows.
/// B_i*(r_i-r_i0) = -dt*G_i*(r_i-a_i) leaves K_i = dt*G_i*B_i/(B_i+dt*G_i).
/// The remaining two material-activity equations have independent capacities.
fn endpoint(
    initial: [f64; 2],
    rh: [f64; 2],
    vapor_capacity: [f64; 2],
    dt: f64,
    boundary: [f64; 2],
    join: f64,
) -> ([f64; 2], [f64; 2], [f64; 2]) {
    let weight = std::array::from_fn::<_, 2, _>(|i| {
        dt * boundary[i] * vapor_capacity[i] / (vapor_capacity[i] + dt * boundary[i])
    });
    let diagonal = [
        CAPACITY[0] + dt * join + weight[0],
        CAPACITY[1] + dt * join + weight[1],
    ];
    let off = dt * join;
    let rhs = [
        CAPACITY[0] * initial[0] + weight[0] * rh[0],
        CAPACITY[1] * initial[1] + weight[1] * rh[1],
    ];
    let determinant = diagonal[0] * diagonal[1] - off * off;
    let activity = [
        (diagonal[1] * rhs[0] + off * rhs[1]) / determinant,
        (off * rhs[0] + diagonal[0] * rhs[1]) / determinant,
    ];
    let transfer = std::array::from_fn(|i| weight[i] * (rh[i] - activity[i]));
    let final_rh = std::array::from_fn(|i| rh[i] - transfer[i] / vapor_capacity[i]);
    (activity, final_rh, transfer)
}

#[test]
fn g1_heterogeneous_layers_and_finite_chambers_match_independent_elimination() {
    with_cx(|cx| {
        let ports = [
            Some(MoistureAtmospherePort::First),
            Some(MoistureAtmospherePort::Second),
        ];
        for finite_transfer in [false, true] {
            let mut moisture = stack([0.3, 0.65], cx);
            let (boundary, join) = if finite_transfer {
                moisture = moisture
                    .with_transfer(
                        [
                            MoistureActivityTransfer::ConductanceKgM2S(40.0),
                            MoistureActivityTransfer::ConductanceKgM2S(60.0),
                        ],
                        vec![MoistureActivityTransfer::ConductanceKgM2S(30.0)],
                        cx,
                    )
                    .unwrap();
                // Actual area is 0.001 m². Bulk half-cell and authored surface/
                // interface resistances add; unlike materials are never averaged.
                (
                    [
                        1.0 / (1.0 / 0.016 + 1.0 / 0.040),
                        1.0 / (1.0 / 0.048 + 1.0 / 0.060),
                    ],
                    1.0 / (1.0 / 0.012 + 1.0 / 0.030),
                )
            } else {
                ([0.016, 0.048], 0.012)
            };
            let identity = moisture.definition_identity();
            let mut gas = [air::atmosphere(0.85, 0.3), air::atmosphere(0.1, 0.7)];
            let capacity = gas
                .each_ref()
                .map(|g| g.vapor_mass_kg().unwrap() / g.relative_humidity().unwrap());
            let carrier = gas.each_ref().map(air::carrier_mass);
            for _ in 0..2 {
                let dt = 0.4;
                let initial: [f64; 2] = moisture.activity().try_into().unwrap();
                let rh = gas.each_ref().map(|g| g.relative_humidity().unwrap());
                let expected = endpoint(initial, rh, capacity, dt, boundary, join);
                let before_mass =
                    gas[0].gas().mass_kg() + gas[1].gas().mass_kg() + moisture.wet_mass_kg();
                let candidate = gas[0]
                    .advance_paired_layered(
                        &gas[1],
                        &moisture,
                        dt,
                        &ports,
                        config(),
                        nonlinear(),
                        controls(),
                        cx,
                    )
                    .unwrap();
                let report = candidate.step();
                assert_eq!(report.face_atmospheres, ports);
                assert!(report.field_evaluations > 1);
                assert!(report.field_evaluations <= controls().maximum_evaluations);
                assert_eq!(candidate.moisture().definition_identity(), identity);
                for i in 0..2 {
                    near(candidate.moisture().activity()[i], expected.0[i], 2e-10);
                    near(
                        candidate.moisture().layers()[i].slab().moisture()[0],
                        SORPTION[i] * expected.0[i],
                        1e-10,
                    );
                    near(
                        candidate.atmospheres()[i].relative_humidity().unwrap(),
                        expected.1[i],
                        2e-10,
                    );
                    near(
                        report.atmospheres[i].water_to_solids_kg,
                        expected.2[i],
                        2e-12,
                    );
                    near(report.transport.boundary_water_kg[i], expected.2[i], 2e-12);
                    near(
                        air::carrier_mass(&candidate.atmospheres()[i]),
                        carrier[i],
                        2e-13,
                    );
                    assert_eq!(
                        report.transport.sorption_sources[i],
                        moisture.layers()[i].sorption().identity()
                    );
                    assert_eq!(
                        report.transport.transport_end[i].card_identity(),
                        moisture.layers()[i].slab().card().content_hash()
                    );
                    assert_eq!(
                        candidate.atmospheres()[i].source_catalog(),
                        gas[i].source_catalog()
                    );
                    let surface = if finite_transfer {
                        expected.1[i] - expected.2[i] / (dt * [0.040, 0.060][i])
                    } else {
                        expected.1[i]
                    };
                    near(
                        report.transport.surface_activities[i].unwrap(),
                        surface,
                        2e-10,
                    );
                }
                near(
                    report.transport.interfaces[0].water_left_to_right_kg,
                    dt * join * (expected.0[0] - expected.0[1]),
                    2e-12,
                );
                near(
                    report.water_balance_residual_kg,
                    0.0,
                    controls().mass_tolerance_kg,
                );
                near(
                    report.total_mass_residual_kg,
                    0.0,
                    controls().mass_tolerance_kg,
                );
                near(
                    candidate.atmospheres()[0].gas().mass_kg()
                        + candidate.atmospheres()[1].gas().mass_kg()
                        + candidate.moisture().wet_mass_kg(),
                    before_mass,
                    controls().mass_tolerance_kg,
                );
                let (next_gas, next_field, _) = candidate.into_parts();
                gas = next_gas;
                moisture = next_field;
            }
        }
    });
}

#[test]
fn g3_zero_net_uptake_preserves_throughflow_and_sealed_or_failed_chambers() {
    with_cx(|cx| {
        let gas = [air::atmosphere(0.8, 1.0), air::atmosphere(0.2, 1.0)];
        let saved = gas.clone();
        let rh = gas.each_ref().map(|g| g.relative_humidity().unwrap());
        let capacity = gas
            .each_ref()
            .map(|g| g.vapor_mass_kg().unwrap() / g.relative_humidity().unwrap());
        let dt = 0.5;
        let solve = |a| endpoint([a; 2], rh, capacity, dt, [0.016, 0.048], 0.012);
        let change = |a| {
            let endpoint = solve(a).0;
            CAPACITY[0] * (endpoint[0] - a) + CAPACITY[1] * (endpoint[1] - a)
        };
        // A constant initial activity whose independently computed net uptake is
        // zero still transmits water from the wet chamber to the dry chamber.
        let initial = change(0.0) / (change(0.0) - change(1.0));
        let moisture = stack([initial; 2], cx);
        let ports = [
            Some(MoistureAtmospherePort::First),
            Some(MoistureAtmospherePort::Second),
        ];
        let candidate = gas[0]
            .advance_paired_layered(
                &gas[1],
                &moisture,
                dt,
                &ports,
                config(),
                nonlinear(),
                controls(),
                cx,
            )
            .unwrap();
        near(
            candidate.step().water_change_kg,
            0.0,
            controls().mass_tolerance_kg,
        );
        let transfers = candidate.step().transport.boundary_water_kg;
        assert!(transfers[0] > 1e-5 && transfers[1] < -1e-5);
        near(transfers[0], -transfers[1], controls().mass_tolerance_kg);
        assert!(
            candidate.atmospheres()[0].vapor_mass_kg().unwrap() < gas[0].vapor_mass_kg().unwrap()
        );
        assert!(
            candidate.atmospheres()[1].vapor_mass_kg().unwrap() > gas[1].vapor_mass_kg().unwrap()
        );

        let reversed = gas[1]
            .advance_paired_layered(
                &gas[0],
                &moisture,
                dt,
                &[ports[1], ports[0]],
                config(),
                nonlinear(),
                controls(),
                cx,
            )
            .unwrap();
        for i in 0..2 {
            near(
                reversed.moisture().activity()[i],
                candidate.moisture().activity()[i],
                2e-10,
            );
            near(
                reversed.atmospheres()[1 - i].vapor_mass_kg().unwrap(),
                candidate.atmospheres()[i].vapor_mass_kg().unwrap(),
                2e-12,
            );
        }

        let mut limited = controls();
        limited.maximum_evaluations = 2;
        assert!(matches!(
            gas[0].advance_paired_layered(
                &gas[1],
                &moisture,
                dt,
                &ports,
                config(),
                nonlinear(),
                limited,
                cx,
            ),
            Err(AtmosphereMoistureError::EvaluationBudget {
                evaluations: 2,
                maximum: 2
            })
        ));
        assert_eq!(gas, saved);
        assert_eq!(moisture.activity(), &[initial; 2]);
        assert_eq!(moisture.time_s(), 0.0);

        let sealed = moisture
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
        let candidate = gas[0]
            .advance_paired_layered(
                &gas[1],
                &sealed,
                dt,
                &ports,
                config(),
                nonlinear(),
                controls(),
                cx,
            )
            .unwrap();
        assert_eq!(&candidate.atmospheres()[1], &gas[1]);
        assert_eq!(candidate.step().face_atmospheres, [ports[0], None]);
        assert_eq!(candidate.step().transport.boundary_water_kg[1], 0.0);
        assert_eq!(candidate.step().transport.surface_activities[1], None);
    });
}
