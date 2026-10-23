//! Model-level regressions: synthetic coolant cards, not material measurements.
#[path = "finite_gas_convection/source.rs"]
#[allow(dead_code)]
mod gas_source;

use fs_airflow::graph::thermal::medium::{FlowHeatGeometry, FrozenFlowMedium, MediumFlowError};
use fs_airflow::LossNetwork;
use fs_convection::{CorrelationId, ThermalDirection};
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
use fs_material::gas::{ConductivityModel, GasSpec, GasState};
use fs_material::liquid::{ResolvedLiquidState, resolve_liquid_state};
use fs_material::state_point::MaterialPropertySelection;
use fs_qty::semantic::{Composition, CompositionBasis, QuantityKind, SemanticType, ValueForm};
use fs_qty::{Area, Dims, Length, QuantitySpec, Temperature, VolumetricFlowRate};

fn liquid(name: &str, rho: f64, cp: f64, mu: f64, k: f64) -> ResolvedLiquidState {
    let mut claims = ClaimSet::new();
    for (name, value, dims, quantity) in [
        ("density", rho, Dims([-3,1,0,0,0,0]), QuantitySpec::dimensional(Dims([-3,1,0,0,0,0]))),
        ("specific_heat_capacity", cp, Dims([2,0,-2,-1,0,0]), QuantitySpec::dimensional(Dims([2,0,-2,-1,0,0]))),
        ("dynamic_viscosity", mu, Dims([-1,1,-1,0,0,0]), QuantitySpec::dimensional(Dims([-1,1,-1,0,0,0]))),
        ("thermal_conductivity", k, Dims([1,1,-3,-1,0,0]), QuantitySpec::semantic(
            SemanticType::new(QuantityKind::ThermalConductivity, ValueForm::Static))),
    ] {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::with_quantity(name, quantity),
            value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained(),
            uncertainty: UncertaintyModel::RelativeHalfWidth { fraction: 0.01, confidence: 0.95 },
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: Provenance { source: "synthetic coolant wiring fixture".into(),
                license: "test-only".into(), artifact: None },
        }).unwrap();
    }
    let card = MaterialCard::assemble(MaterialStateId {
        chemistry: name.into(), phase: "liquid".into(), process: "synthetic constant coefficients".into(), revision: 0,
    }, claims, vec![]).unwrap();
    resolve_liquid_state(&card, &QueryPoint::new(), MaterialPropertySelection::SingleClaimOnly).unwrap()
}

fn geometry() -> FlowHeatGeometry {
    FlowHeatGeometry { flow_area: Area::new(core::f64::consts::PI * 0.01 * 0.01 / 4.0), characteristic_length: Length::new(0.01),
        heated_length: Length::new(2.0), wetted_area: Area::new(core::f64::consts::PI * 0.01 * 2.0), aspect_ratio: None }
}
fn near(a: f64, b: f64) {
    assert!((a-b).abs() <= 2e-12 * a.abs().max(b.abs()), "{a} != {b}");
}

#[test]
fn liquid_replacement_updates_pressure_loading_robin_and_advected_heat_together() {
    let first = liquid("synthetic-a", 1000.0, 4000.0, 0.001, 0.6);
    let second = liquid("synthetic-b", 900.0, 1800.0, 0.0015, 0.3);
    let a = FrozenFlowMedium::Liquid(&first);
    let b = FrozenFlowMedium::Liquid(&second);
    let flow = VolumetricFlowRate::new(1e-6);
    let ca = a.channel("wall", flow, Temperature::new(300.0), geometry(),
        CorrelationId::CircularDuctLaminarCwt, ThermalDirection::HeatingFluid).unwrap();
    let cb = b.channel("wall", flow, Temperature::new(300.0), geometry(),
        CorrelationId::CircularDuctLaminarCwt, ThermalDirection::HeatingFluid).unwrap();
    assert_eq!(ca.medium(), a);
    assert_eq!(cb.medium(), b);
    near(ca.path().capacity_rate_w_per_k(), 4.0);
    near(cb.path().capacity_rate_w_per_k(), 1.62);
    near(cb.coefficient().value.value(), ca.coefficient().value.value() * 0.5);
    assert_ne!(ca.correlation(), cb.correlation());
    let ra = LossNetwork::Element(a.orifice_loss("opening", Area::new(1e-4)).unwrap().element().clone()).equivalent_resistance().value();
    let rb = LossNetwork::Element(b.orifice_loss("opening", Area::new(1e-4)).unwrap().element().clone()).equivalent_resistance().value();
    near(rb, 0.9 * ra);
    let transport = b.transport_air().unwrap();
    near(transport.density.value(), second.density_kg_m3());
    near(transport.specific_heat_j_kg_k, second.specific_heat_j_kg_k());
    let robin = cb.robin_boundary(Temperature::new(310.0)).unwrap();
    near(robin.coefficient().value.value(), cb.coefficient().value.value());
    match robin.boundary_condition() {
        fs_conduction::ThermalBc::Robin { htc, t_ref } => {
            near(htc.at(0), cb.path().segments()[0].htc_w_per_m2_k());
            near(t_ref.at(0), 310.0);
        }
        other => panic!("wrong boundary: {other:?}"),
    }
    for channel in [&ca, &cb] {
        let march = channel.path().march(&[330.0]).unwrap();
        near(march.total_heat_rate_w, channel.path().capacity_rate_w_per_k() * (march.outlet_temperature_k - 300.0));
        assert!(march.outlet_temperature_k > 300.0 && march.outlet_temperature_k < 330.0);
    }
    assert_ne!(ca.path().march(&[330.0]).unwrap().outlet_temperature_k,
        cb.path().march(&[330.0]).unwrap().outlet_temperature_k);
}

#[test]
fn gas_composition_changes_a_real_heat_channel_without_air_defaults() {
    let air_spec = GasSpec::dry_air_ussa1976();
    let other = GasSpec { molar_mass: 0.004, gamma: 5.0/3.0, sutherland_beta: 1e-6,
        sutherland_s: 80.0, conductivity: ConductivityModel::Eucken };
    let air = GasState::try_new(&air_spec, 300.0, 101325.0).unwrap();
    let composition = Composition::new(CompositionBasis::MoleFraction, vec![0.5,0.5]).unwrap();
    let mixture = GasState::try_new_mixture(&[air_spec,other], &composition, 300.0, 101325.0, None).unwrap();
    let ca = FrozenFlowMedium::DeclaredGas(&air).channel("wall", VolumetricFlowRate::new(1e-4),
        Temperature::new(300.0), geometry(), CorrelationId::CircularDuctLaminarCwt, ThermalDirection::HeatingFluid).unwrap();
    let cm = FrozenFlowMedium::DeclaredGas(&mixture).channel("wall", VolumetricFlowRate::new(1e-4),
        Temperature::new(300.0), geometry(), CorrelationId::CircularDuctLaminarCwt, ThermalDirection::HeatingFluid).unwrap();
    assert_ne!(ca.coefficient().value.value(), cm.coefficient().value.value());
    assert_ne!(ca.path().capacity_rate_w_per_k(), cm.path().capacity_rate_w_per_k());
    near(cm.path().capacity_rate_w_per_k(), mixture.density * 1e-4 * mixture.specific_heat_cp);
    near(cm.coefficient().value.value(), 3.66 * mixture.thermal_conductivity / 0.01);
}

#[test]
fn viscosity_and_flow_reach_correlation_admission_not_an_unrelated_knob() {
    let water_like = liquid("synthetic-a", 1000.0, 4000.0, 0.001, 0.6);
    let medium = FrozenFlowMedium::Liquid(&water_like);
    let make = |flow| medium.channel("wall", VolumetricFlowRate::new(flow), Temperature::new(300.0),
        geometry(), CorrelationId::CircularDuctHausen, ThermalDirection::HeatingFluid).unwrap();
    let a = make(1e-6);
    let b = make(2e-6);
    assert!(b.coefficient().value.value() > a.coefficient().value.value());
    near(b.path().capacity_rate_w_per_k(), 2.0 * a.path().capacity_rate_w_per_k());
    let refused = medium.channel("wall", VolumetricFlowRate::new(1e-6), Temperature::new(300.0),
        geometry(), CorrelationId::Gnielinski, ThermalDirection::HeatingFluid);
    assert!(matches!(refused, Err(MediumFlowError::Correlation(_))));
    let viscous = liquid("synthetic-viscous", 1000.0, 4000.0, 1.0, 0.6);
    assert!(matches!(FrozenFlowMedium::Liquid(&viscous).channel("wall", VolumetricFlowRate::new(1e-6),
        Temperature::new(300.0), geometry(), CorrelationId::CircularDuctHausen, ThermalDirection::HeatingFluid),
        Err(MediumFlowError::Correlation(_))));
    assert!(matches!(medium.channel("wall", VolumetricFlowRate::new(1e-6), Temperature::new(300.0),
        geometry(), CorrelationId::CircularDuctLaminarChf, ThermalDirection::HeatingFluid),
        Err(MediumFlowError::UnsupportedCorrelation(_))));
}

#[test]
fn invalid_geometry_flow_and_material_refuse_without_a_partial_channel() {
    let mut gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 101325.0).unwrap();
    for q in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        assert!(FrozenFlowMedium::DeclaredGas(&gas).channel("wall", VolumetricFlowRate::new(q),
            Temperature::new(300.0), geometry(), CorrelationId::CircularDuctLaminarCwt,
            ThermalDirection::HeatingFluid).is_err());
    }
    let bad_geometry = FlowHeatGeometry { flow_area: Area::new(0.0), ..geometry() };
    assert!(FrozenFlowMedium::DeclaredGas(&gas).channel("wall", VolumetricFlowRate::new(1e-4),
        Temperature::new(300.0), bad_geometry, CorrelationId::CircularDuctLaminarCwt,
        ThermalDirection::HeatingFluid).is_err());
    gas.thermal_conductivity = f64::NAN;
    assert!(FrozenFlowMedium::DeclaredGas(&gas).transport_air().is_err());
}

#[test]
fn thermochemical_pack_changes_reach_solved_flow_and_series_heat_with_original_sources() {
    use fs_airflow::graph::thermal::medium::FlowHeatSection;
    use fs_airflow::{EnclosureNetwork, FanArrangement, FanBank, FanCurve, FanPoint,
        LeakageElement, SourceProvenance, ToleranceBasis, SHARP_EDGED_ORIFICE_CD,
        solve_operating_point};
    use fs_material::fluid_transport::FluidTransport;
    use fs_qty::Pressure;

    // An authored regulated pressure source, not a measured fan rescaled to
    // another atmosphere. Its unchanged boundary law supplies 2000 Pa.
    let drive = FanBank::new(FanCurve::new("controlled pressure fixture", vec![
        FanPoint::new(VolumetricFlowRate::new(0.0), Pressure::new(2000.0)),
        FanPoint::new(VolumetricFlowRate::new(0.05), Pressure::new(2000.0)),
    ], SourceProvenance::new("manufactured regulated source", "constant-pressure-v1"),
        0.0, ToleranceBasis::Analytic, VolumetricFlowRate::new(0.0), (1.0, 1.0))
        .unwrap(), 1, FanArrangement::Series, 1.0).unwrap();
    let patch = FlowHeatGeometry { flow_area: Area::new(0.01),
        characteristic_length: Length::new(0.2), heated_length: Length::new(0.4),
        wetted_area: Area::new(0.02), aspect_ratio: None };
    let sections = ["upstream", "downstream"].map(|region| FlowHeatSection {
        region, geometry: patch, correlation: CorrelationId::FlatPlateLaminarAverage,
        direction: ThermalDirection::HeatingFluid,
    });
    let run = |fraction, pressure| {
        let model = gas_source::model(fraction);
        let state = model.evaluate(&gas_source::point(400.0, pressure)).unwrap();
        let medium = FrozenFlowMedium::ThermochemicalGas(&state);
        let transport = state.transport_state().unwrap();
        assert_eq!(medium.transport_state().unwrap(), transport);
        let primary = medium.orifice_loss("primary", Area::new(4e-4)).unwrap();
        let leakage = medium.orifice_loss("leakage", Area::new(4e-4)).unwrap();
        assert_eq!(primary.medium(), medium);
        let network = EnclosureNetwork::new(LossNetwork::Element(primary.element().clone()),
            LeakageElement::new(leakage.element().clone()));
        let solved = solve_operating_point(&drive, &network).unwrap();
        let branch = solved.branches.iter().find(|b| b.path == "primary").unwrap();
        let q = branch.flow.value.value();
        let expected_q = SHARP_EDGED_ORIFICE_CD * 4e-4
            * (2.0 * 2000.0 / state.state().density).sqrt();
        assert!((q / expected_q - 1.0).abs() < 1e-8);
        let series = medium.channel_series(branch.flow.value, Temperature::new(400.0),
            &sections).unwrap();
        assert_eq!(series.medium(), medium);
        let retained = match series.medium() {
            FrozenFlowMedium::ThermochemicalGas(original) => original,
            _ => panic!("lost thermochemical source identity"),
        };
        assert!(core::ptr::eq(retained, &state));
        assert_eq!(retained.source_packs(), model.source_packs());
        assert_eq!(retained.source_packs().len(), 2);
        let capacity = state.state().density * q * state.state().specific_heat_cp;
        near(series.path().capacity_rate_w_per_k(), capacity);
        let mut expected_outlet = 400.0;
        for (exchange, wall) in series.exchanges().iter().zip([420.0, 430.0]) {
            assert_eq!(exchange.medium(), medium);
            assert_eq!(*exchange.transport_state(), transport);
            let h = exchange.coefficient().value.value();
            // Independent constant-wall-temperature segment solution.
            expected_outlet = wall + (expected_outlet - wall)
                * (-h * patch.wetted_area.value() / capacity).exp();
            let robin = exchange.robin_boundary(Temperature::new(expected_outlet)).unwrap();
            near(robin.coefficient().value.value(), h);
        }
        let result = series.path().march(&[420.0, 430.0]).unwrap();
        near(result.outlet_temperature_k, expected_outlet);
        near(result.total_heat_rate_w, capacity * (expected_outlet - 400.0));
        (q, capacity, result.outlet_temperature_k, result.total_heat_rate_w,
            series.exchanges()[0].coefficient().value.value())
    };
    let reference = run(0.75, 100000.0);
    let pressurized = run(0.75, 200000.0);
    let substituted = run(0.25, 100000.0);
    assert!((pressurized.0 / reference.0 - 1.0 / 2.0_f64.sqrt()).abs() < 1e-8);
    assert!((pressurized.1 / reference.1 - 2.0_f64.sqrt()).abs() < 1e-8);
    assert_ne!(reference.2, pressurized.2);
    assert_ne!(reference.0, substituted.0);
    assert_ne!(reference.1, substituted.1);
    assert_ne!(reference.3, substituted.3);
    assert_ne!(reference.4, substituted.4);
}

#[test]
fn thermochemical_channel_retains_correlation_refusal_without_replacing_the_medium() {
    let model = gas_source::model(0.75);
    let state = model.evaluate(&gas_source::point(400.0, 100000.0)).unwrap();
    let saved = state.clone();
    let medium = FrozenFlowMedium::ThermochemicalGas(&state);
    assert!(matches!(medium.channel("wall", VolumetricFlowRate::new(1.0),
        Temperature::new(400.0), geometry(), CorrelationId::CircularDuctLaminarCwt,
        ThermalDirection::HeatingFluid), Err(MediumFlowError::Correlation(_))));
    assert_eq!(state, saved);
    // A supported state remains usable after refusal; no fallback gas or
    // retained failed channel changes its transport or source authority.
    let channel = medium.channel("wall", VolumetricFlowRate::new(1e-4),
        Temperature::new(400.0), geometry(), CorrelationId::CircularDuctLaminarCwt,
        ThermalDirection::HeatingFluid).unwrap();
    assert_eq!(channel.medium(), medium);
    assert_eq!(state, saved);
}
