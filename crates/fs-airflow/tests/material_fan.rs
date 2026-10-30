//! G1/G3: one source gas changes fan pressure, losses and transported heat.
//! Manufactured gas/curve data test the coupling, not fan performance accuracy.
#[path = "finite_gas_convection/source.rs"]
#[allow(dead_code)]
mod gas_source;

use fs_airflow::composite::{compose_parallel, compose_series};
use fs_airflow::graph::sensitivity::{GraphFanSensitivityError, GraphSensitivityConfig};
use fs_airflow::graph::thermal::medium::fan::{FanGasDomain, MaterialFanError};
use fs_airflow::graph::thermal::medium::{FlowHeatGeometry, FrozenFlowMedium};
use fs_airflow::graph::{GraphBranch, GraphFanConfig, GraphSolveConfig, LossGraph};
use fs_airflow::{
    EnclosureNetwork, FanArrangement, FanBank, FanCurve, FanPoint, LeakageElement, LossElement,
    LossNetwork, LossResistance, SourceProvenance, ToleranceBasis, solve_operating_point,
};
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_convection::{CorrelationId, ThermalDirection};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::gas::{GasSpec, GasState};
use fs_qty::{Area, Length, Pressure, Temperature, VolumetricFlowRate};
use fs_solver::CheckedCgConfig;

fn bank(speed: f64) -> FanBank {
    FanBank::new(
        FanCurve::new(
            "reference gas fan",
            vec![
                FanPoint::new(VolumetricFlowRate::new(0.0), Pressure::new(800.0)),
                FanPoint::new(VolumetricFlowRate::new(0.04), Pressure::new(0.0)),
            ],
            SourceProvenance::new("manufactured static fan curve", "material-fan-fixture"),
            0.0,
            ToleranceBasis::Analytic,
            VolumetricFlowRate::new(0.0),
            (0.5, 1.5),
        )
        .unwrap(),
        1,
        FanArrangement::Series,
        speed,
    )
    .unwrap()
}
fn domain() -> FanGasDomain {
    FanGasDomain {
        density_ratio: [0.1, 5.0],
        reynolds_ratio: [0.1, 5.0],
        maximum_pressure_fraction: 0.05,
    }
}
fn close(a: f64, b: f64, relative: f64) {
    assert!(
        (a - b).abs() <= relative * a.abs().max(b.abs()).max(1e-20),
        "{a:e} != {b:e}"
    );
}
fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 91,
                kernel_id: 839,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}
fn graph_config() -> GraphSolveConfig {
    GraphSolveConfig {
        max_sweeps: 128,
        max_node_iterations: 80,
        absolute_flow_tolerance: VolumetricFlowRate::new(1e-12),
        relative_flow_tolerance: 1e-12,
    }
}
fn fan_config() -> GraphFanConfig {
    GraphFanConfig {
        max_iterations: 128,
        absolute_flow_tolerance: VolumetricFlowRate::new(1e-12),
        absolute_pressure_tolerance: Pressure::new(1e-8),
    }
}
fn sensitivity() -> GraphSensitivityConfig {
    GraphSensitivityConfig {
        minimum_absolute_flow: VolumetricFlowRate::new(1e-12),
        linear: CheckedCgConfig {
            tolerance: 1e-12,
            max_iterations: 64,
            max_corrections: 2,
        },
    }
}
fn graph(resistance: f64) -> LossGraph {
    LossGraph::new(
        2,
        vec![GraphBranch {
            from: 0,
            to: 1,
            loss: LossElement::new(
                "analytic loss",
                LossResistance::new(resistance),
                0.0,
                SourceProvenance::new("analytic resistance", "material-fan-loss"),
                ToleranceBasis::Analytic,
            )
            .unwrap(),
        }],
    )
    .unwrap()
}

#[test]
fn original_mixture_drives_fan_network_and_heat_with_no_density_reentry() {
    let reference_model = gas_source::model(0.75);
    let reference = reference_model
        .evaluate(&gas_source::point(400.0, 100000.0))
        .unwrap();
    let reference_medium = FrozenFlowMedium::ThermochemicalGas(&reference);
    let source = bank(1.0);
    let saved_source = source.clone();
    let mut results = Vec::new();
    for (fraction, temperature, pressure) in [
        (0.75, 400.0, 100000.0),
        (0.75, 400.0, 200000.0),
        (0.25, 400.0, 100000.0),
        (0.75, 500.0, 100000.0),
    ] {
        let model = gas_source::model(fraction);
        let state = model
            .evaluate(&gas_source::point(temperature, pressure))
            .unwrap();
        let medium = FrozenFlowMedium::ThermochemicalGas(&state);
        let fan = medium
            .bind_fan_bank(&source, reference_medium, domain())
            .unwrap();
        assert_eq!(fan.curve(), source.curve());
        assert_eq!(fan.reference_medium(), reference_medium);
        let FrozenFlowMedium::ThermochemicalGas(retained) = fan.inlet_medium() else {
            panic!("lost original thermochemical source");
        };
        assert!(core::ptr::eq(retained, &state));
        assert_eq!(retained.source_packs(), model.source_packs());
        let ratio = state.state().density / reference.state().density;
        close(
            fan.gas_density_binding().unwrap().density_ratio(),
            ratio,
            1e-14,
        );

        let primary = medium.orifice_loss("primary", Area::new(4e-4)).unwrap();
        let leakage = medium.orifice_loss("leakage", Area::new(4e-4)).unwrap();
        let network = EnclosureNetwork::new(
            LossNetwork::Element(primary.element().clone()),
            LeakageElement::new(leakage.element().clone()),
        );
        // The existing interval-Newton producer receives the bound bank.
        let point = solve_operating_point(&fan, &network).unwrap();
        let q = point.flow.value.value();
        let resistance = network.equivalent_resistance().value();
        let slope = ratio * 800.0 / 0.04;
        let expected = 2.0 * ratio * 800.0
            / (slope + (slope * slope + 4.0 * resistance * ratio * 800.0).sqrt());
        close(q, expected, 2e-8);
        close(
            point.pressure.value.value(),
            ratio * (800.0 - 20000.0 * q),
            2e-8,
        );
        let channel = medium
            .channel(
                "wall",
                VolumetricFlowRate::new(0.5 * q),
                Temperature::new(temperature),
                FlowHeatGeometry {
                    flow_area: Area::new(0.01),
                    characteristic_length: Length::new(0.2),
                    heated_length: Length::new(0.4),
                    wetted_area: Area::new(0.02),
                    aspect_ratio: None,
                },
                CorrelationId::FlatPlateLaminarAverage,
                ThermalDirection::HeatingFluid,
            )
            .unwrap();
        let capacity = 0.5 * q * state.state().density * state.state().specific_heat_cp;
        close(channel.path().capacity_rate_w_per_k(), capacity, 1e-14);
        let thermal = channel.path().march(&[temperature + 30.0]).unwrap();
        close(
            thermal.total_heat_rate_w,
            capacity * (thermal.outlet_temperature_k - temperature),
            2e-12,
        );
        results.push((
            q,
            point.pressure.value.value(),
            capacity,
            thermal.total_heat_rate_w,
        ));
    }
    // A same-speed fan and its quadratic losses both scale with density. The
    // original pressure-only source fixture cannot reproduce this invariant.
    for result in &results[1..] {
        close(result.0, results[0].0, 2e-8);
    }
    close(results[1].1, 2.0 * results[0].1, 2e-8);
    close(results[1].2, 2.0 * results[0].2, 2e-8);
    assert_ne!(results[0].2, results[2].2);
    assert_ne!(results[0].3, results[1].3);
    assert_eq!(source, saved_source);
}

#[test]
fn density_factor_reaches_graph_operating_point_and_original_curve_derivatives() {
    let reference = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 100000.0).unwrap();
    let inlet = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 200000.0).unwrap();
    let medium = FrozenFlowMedium::DeclaredGas(&inlet);
    let fan = medium
        .bind_fan_bank(
            &bank(1.0),
            FrozenFlowMedium::DeclaredGas(&reference),
            domain(),
        )
        .unwrap();
    let resistance = 2e6;
    let network = graph(resistance);
    with_cx(|cx| {
        let solved = network
            .solve_with_fan(0, 1, &fan, graph_config(), fan_config(), cx)
            .unwrap();
        let slope = 40000.0;
        let q = 3200.0 / (slope + (slope * slope + 4.0 * resistance * 1600.0).sqrt());
        close(solved.flow.value(), q, 2e-9);
        assert_eq!(solved.fan.gas_density_binding(), fan.gas_density_binding());
        let linear = network
            .solve_with_fan_linearized(0, 1, &fan, graph_config(), fan_config(), sensitivity(), cx)
            .unwrap();
        let speed = linear
            .tangent(&[0.0], 1.0, &[Pressure::new(0.0); 2], cx)
            .unwrap();
        close(speed.fan_flow.value(), q, 2e-9);
        let sample = linear
            .tangent(&[0.0], 0.0, &[Pressure::new(1.0), Pressure::new(0.0)], cx)
            .unwrap();
        let expected = 2.0 * (1.0 - q / 0.04) / (2.0 * resistance * q + slope);
        close(sample.fan_flow.value(), expected, 2e-9);
    });
}

#[test]
fn material_density_scales_series_parallel_composition_once_and_preserves_identity() {
    let air = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 100000.0).unwrap();
    let dense = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 200000.0).unwrap();
    let source = bank(1.0);
    let a = FrozenFlowMedium::DeclaredGas(&dense)
        .bind_fan_bank(&source, FrozenFlowMedium::DeclaredGas(&dense), domain())
        .unwrap();
    let b = FrozenFlowMedium::DeclaredGas(&dense)
        .bind_fan_bank(&source, FrozenFlowMedium::DeclaredGas(&air), domain())
        .unwrap();
    let members = [a.bank().clone(), b.bank().clone()];
    let reversed = [members[1].clone(), members[0].clone()];
    let series = compose_series(&members).unwrap();
    let parallel = compose_parallel(&members).unwrap();
    close(
        series
            .pressure_at(VolumetricFlowRate::new(0.01))
            .unwrap()
            .value(),
        1800.0,
        1e-13,
    );
    // At 600 Pa: q_a=.01, q_b=.025, so q_parallel=.035.
    close(
        parallel
            .pressure_at(VolumetricFlowRate::new(0.035))
            .unwrap()
            .value(),
        600.0,
        1e-13,
    );
    assert_eq!(series, compose_series(&reversed).unwrap());
    assert_eq!(parallel, compose_parallel(&reversed).unwrap());
    for composed in [&series, &parallel] {
        assert!(composed.has_material_density());
        assert!(
            FrozenFlowMedium::DeclaredGas(&dense)
                .bind_fan_bank(composed, FrozenFlowMedium::DeclaredGas(&air), domain())
                .is_err()
        );
    }

    // Same density with different T/p and viscosity changes the physical input
    // identity even where the nominal static-pressure curve is identical.
    let hot = GasState::try_new(&GasSpec::dry_air_ussa1976(), 600.0, 400000.0).unwrap();
    let c = FrozenFlowMedium::DeclaredGas(&hot)
        .bind_fan_bank(&source, FrozenFlowMedium::DeclaredGas(&hot), domain())
        .unwrap();
    let d = FrozenFlowMedium::DeclaredGas(&hot)
        .bind_fan_bank(&source, FrozenFlowMedium::DeclaredGas(&air), domain())
        .unwrap();
    let same_pressure = compose_series(&[c.bank().clone(), d.bank().clone()]).unwrap();
    assert_eq!(series.curve().points(), same_pressure.curve().points());
    assert_ne!(series.curve().source(), same_pressure.curve().source());
    for unsupported in [
        vec![a.bank().clone(), c.bank().clone()],
        vec![a.bank().clone(), source.clone()],
    ] {
        assert!(matches!(
            compose_series(&unsupported),
            Err(fs_airflow::AirflowError::MaterialFanComposition { .. })
        ));
        assert!(matches!(
            compose_parallel(&unsupported),
            Err(fs_airflow::AirflowError::MaterialFanComposition { .. })
        ));
    }
    let nested = compose_series(&[series, a.bank().clone()]).unwrap();
    assert_eq!(nested.inlet_gas(), Some(&dense));
    close(
        nested
            .pressure_at(VolumetricFlowRate::new(0.01))
            .unwrap()
            .value(),
        2400.0,
        1e-13,
    );
    let narrow = FrozenFlowMedium::DeclaredGas(&dense)
        .bind_fan_bank(
            &source,
            FrozenFlowMedium::DeclaredGas(&dense),
            FanGasDomain {
                maximum_pressure_fraction: 0.005,
                ..domain()
            },
        )
        .unwrap();
    let pair = [narrow.bank().clone(), narrow.bank().clone()];
    assert!(matches!(
        compose_series(&pair),
        Err(fs_airflow::AirflowError::MaterialFanComposition { .. })
    ));
    assert!(compose_parallel(&pair).is_ok());
}

#[test]
fn material_refusals_leave_original_curve_available_and_never_apply_twice() {
    let mut air = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 100000.0).unwrap();
    let dense = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 200000.0).unwrap();
    let source = bank(1.0);
    let saved = source.clone();
    let medium = FrozenFlowMedium::DeclaredGas(&dense);
    let reference = FrozenFlowMedium::DeclaredGas(&air);
    for unsupported in [
        FanGasDomain {
            density_ratio: [0.5, 1.5],
            ..domain()
        },
        FanGasDomain {
            reynolds_ratio: [0.5, 1.5],
            ..domain()
        },
        FanGasDomain {
            maximum_pressure_fraction: 0.001,
            ..domain()
        },
    ] {
        assert!(matches!(
            medium.bind_fan_bank(&source, reference, unsupported),
            Err(MaterialFanError::OutsideDomain { .. })
        ));
    }
    let admitted = medium.bind_fan_bank(&source, reference, domain()).unwrap();
    assert!(matches!(
        medium.bind_fan_bank(&admitted, reference, domain()),
        Err(MaterialFanError::Invalid(_))
    ));
    air.density *= 2.0;
    assert!(matches!(
        FrozenFlowMedium::DeclaredGas(&air).bind_fan_bank(
            &source,
            FrozenFlowMedium::DeclaredGas(&dense),
            domain()
        ),
        Err(MaterialFanError::Material(_))
    ));
    assert_eq!(source, saved);
}

#[test]
fn speed_derivative_requires_interior_material_similarity() {
    let air = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 100000.0).unwrap();
    let medium = FrozenFlowMedium::DeclaredGas(&air);
    let fan = medium
        .bind_fan_bank(
            &bank(1.0),
            medium,
            FanGasDomain {
                reynolds_ratio: [1.0, 2.0],
                ..domain()
            },
        )
        .unwrap();
    let network = graph(2e6);
    with_cx(|cx| {
        assert!(
            network
                .solve_with_fan(0, 1, &fan, graph_config(), fan_config(), cx)
                .is_ok()
        );
        assert!(matches!(
            network.solve_with_fan_linearized(
                0,
                1,
                &fan,
                graph_config(),
                fan_config(),
                sensitivity(),
                cx
            ),
            Err(GraphFanSensitivityError::NonDifferentiableFan(_))
        ));
    });
}
