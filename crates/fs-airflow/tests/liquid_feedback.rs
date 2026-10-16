//! Synthetic temperature laws through actual material/correlation/caloric owners.
use fs_airflow::graph::thermal::medium::{FlowHeatGeometry, FlowHeatSection, FrozenFlowMedium};
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::{
    LiquidFeedbackBudget, LiquidHeatPath, LiquidHeatSection, LiquidPropertyTemperature,
};
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_convection::{CorrelationId, ThermalDirection};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel};
use fs_material::liquid_temperature::{LiquidTemperatureLaw, LiquidViscosityBasis};
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::state_point::MaterialPropertySelection;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Area, Dims, Length, MassFlowRate, QuantitySpec, Temperature, VolumetricFlowRate};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| f(&Cx::new(
        &gate, arena, StreamKey { seed: 32, kernel_id: 735, tile: 0, iteration: 0 },
        Budget::INFINITE, ExecMode::Deterministic)))
}
fn point(t: f64) -> QueryPoint {
    QueryPoint::new().with("T", t).unwrap().with("pressure", 100000.0).unwrap()
}
fn card(variable_transport: bool, variable_cp: bool, viscosity_scale: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, initial, final_value, varying) in [
        ("density", Dims([-3,1,0,0,0,0]), 1000.0, 850.0, variable_transport),
        ("specific_heat_capacity", Dims([2,0,-2,-1,0,0]), 1500.0, 3500.0, variable_cp),
        ("dynamic_viscosity", Dims([-1,1,-1,0,0,0]), 0.002*viscosity_scale, 0.001*viscosity_scale, variable_transport),
        ("thermal_conductivity", Dims([1,1,-3,-1,0,0]), 0.4, 0.8, variable_transport),
    ] {
        let quantity = if name == "thermal_conductivity" {
            QuantitySpec::semantic(SemanticType::new(QuantityKind::ThermalConductivity, ValueForm::Static))
        } else { QuantitySpec::dimensional(dims) };
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::with_quantity(name, quantity),
            value: if varying { PropertyValue::Curve { abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
                knots: vec![(250.0, initial), (450.0, final_value)], dims }
            } else { PropertyValue::Scalar { value: initial, dims } },
            validity: ValidityDomain::unconstrained().with("T", 250.0, 450.0).with("pressure", 100000.0, 100000.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: if varying { InterpolationPolicy::LinearInside } else { InterpolationPolicy::ConstantWithinValidity },
            observations: vec![], provenance: Provenance { source: "synthetic feedback equations, not coolant measurements".into(),
                license: "test-only".into(), artifact: None },
        }).unwrap();
    }
    MaterialCard::assemble(MaterialStateId { chemistry: "synthetic coolant".into(), phase: "liquid".into(),
        process: "test-only".into(), revision: 0 }, claims, vec![]).unwrap()
}
fn law(card: &MaterialCard) -> LiquidTemperatureLaw<'_> {
    LiquidTemperatureLaw::try_new(card, &point(300.0), "T", LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly).unwrap()
}
fn caloric(card: &MaterialCard) -> ResolvedHeatCapacityCurve {
    let cp = PropertyKey::with_quantity("specific_heat_capacity", QuantitySpec::dimensional(Dims([2,0,-2,-1,0,0])));
    let rho = PropertyKey::with_quantity("density", QuantitySpec::dimensional(Dims([-3,1,0,0,0,0])));
    ResolvedHeatCapacityCurve::try_new(card, HeatCapacitySourceQuery {
        heat_capacity: &cp, density: &rho, temperature_axis: "T", lower: &point(250.0), upper: &point(450.0),
        selection: SelectionPolicy::SingleClaimOnly,
    }, SolidLiquidPhase::Liquid, 0.0, 0.5).unwrap()
}
fn section(region: &str) -> LiquidHeatSection<'_> {
    LiquidHeatSection { region, geometry: FlowHeatGeometry {
        flow_area: Area::new(1e-4), characteristic_length: Length::new(0.01), heated_length: Length::new(2.0),
        wetted_area: Area::new(0.008), aspect_ratio: None,
    }, correlation: CorrelationId::CircularDuctLaminarCwt }
}
fn budget() -> LiquidFeedbackBudget {
    LiquidFeedbackBudget { temperature_tolerance_k: 1e-9, coefficient_relative_tolerance: 1e-10, max_region_iterations: 96 }
}
fn close(a: f64, b: f64) {
    assert!((a-b).abs() < 2e-9*(1.0+a.abs().max(b.abs())), "{a} != {b}");
}

#[test]
fn constant_properties_match_existing_enthalpy_and_robin_path() {
    let card = card(false, false, 1.0);
    let law = law(&card);
    let source = caloric(&card);
    with_cx(|cx| {
        let sections = [section("first"), section("last")];
        let path = LiquidHeatPath::try_new(cx, &law, &source, MassFlowRate::new(0.001), Temperature::new(300.0),
            &sections, LiquidPropertyTemperature::BulkMean).unwrap();
        let old_sections: Vec<_> = sections.iter().map(|s| FlowHeatSection { region: s.region, geometry: s.geometry,
            correlation: s.correlation, direction: ThermalDirection::HeatingFluid }).collect();
        let old = FrozenFlowMedium::Liquid(law.reference_state()).channel_series(VolumetricFlowRate::new(1e-6),
            Temperature::new(300.0), &old_sections).unwrap();
        let walls = [370.0, 280.0];
        let expected = old.march_liquid_enthalpy(cx, &source, "T", &walls).unwrap();
        let result = path.march(cx, &walls, budget()).unwrap();
        close(result.total_heat_rate_w, expected.total_heat_rate_w());
        for (actual, expected) in result.regions.iter().zip(expected.regions()) {
            close(actual.outlet.temperature_k(), expected.outlet().temperature_k());
            close(actual.boundary.reference_temperature().value(), expected.reference_temperature_k());
            close(actual.boundary.coefficient().value.value(), expected.exchange().coefficient().value.value());
        }
        assert_eq!(result, path.march(cx, &walls, budget()).unwrap());
    });
}

#[test]
fn data_only_transport_changes_reach_heat_and_mass_conservation() {
    let card = card(true, true, 1.0);
    let law = law(&card);
    let source = caloric(&card);
    with_cx(|cx| {
        let path = LiquidHeatPath::try_new(cx, &law, &source, MassFlowRate::new(0.001), Temperature::new(300.0),
            &[section("first"), section("last")], LiquidPropertyTemperature::BulkMean).unwrap();
        let walls = [410.0, 280.0];
        let result = path.march(cx, &walls, budget()).unwrap();
        let mut previous = path.inlet();
        for (region, wall) in result.regions.iter().zip(walls) {
            assert_eq!(region.inlet, previous);
            close(region.volume_flow_m3_s * region.material.density_kg_m3(), 0.001);
            close(region.velocity_m_s, region.volume_flow_m3_s / 1e-4);
            close(region.property_temperature_k, 0.5*(region.inlet.temperature_k()+region.outlet.temperature_k()));
            close(region.heat_rate_w, 0.001*(region.outlet.specific_enthalpy_j_kg()-region.inlet.specific_enthalpy_j_kg()));
            let h = region.boundary.coefficient().value.value();
            close(h, 3.66*region.material.thermal_conductivity_w_m_k()/0.01);
            close(region.heat_rate_w, h*0.008*(wall-region.boundary.reference_temperature().value()));
            assert_eq!(region.material.material().card_identity(), card.content_hash());
            previous = region.outlet;
        }
        assert_ne!(result.regions[0].material.dynamic_viscosity_pa_s(), law.reference_state().dynamic_viscosity_pa_s());
        assert_ne!(result.regions[0].boundary.coefficient().value.value(), result.regions[1].boundary.coefficient().value.value());
        assert!(result.regions[0].heat_rate_w > 0.0 && result.regions[1].heat_rate_w < 0.0);
        close(result.total_heat_rate_w, 0.001*(previous.specific_enthalpy_j_kg()-path.inlet().specific_enthalpy_j_kg()));
    });
}

#[test]
fn film_choice_is_explicit_and_viscosity_changes_turbulent_heat_transfer() {
    let a = card(true, false, 1.0);
    let b = card(true, false, 0.75);
    with_cx(|cx| {
        let run = |card: &MaterialCard, convention| {
            let law = law(card);
            let source = caloric(card);
            let mut section = section("wall");
            section.correlation = CorrelationId::DittusBoelter;
            section.geometry.wetted_area = Area::new(0.08);
            let path = LiquidHeatPath::try_new(cx, &law, &source, MassFlowRate::new(0.4), Temperature::new(300.0),
                &[section], convention).unwrap();
            let result = path.march(cx, &[420.0], budget()).unwrap();
            result.regions[0].clone()
        };
        let bulk = run(&a, LiquidPropertyTemperature::BulkMean);
        let film = run(&a, LiquidPropertyTemperature::FilmMean);
        let thinner = run(&b, LiquidPropertyTemperature::BulkMean);
        assert!(film.property_temperature_k > bulk.property_temperature_k);
        assert!(thinner.correlation.groups()["Re"] > bulk.correlation.groups()["Re"]);
        assert!(thinner.outlet.temperature_k() > bulk.outlet.temperature_k());
        assert_ne!(film.boundary.coefficient().value.value(), bulk.boundary.coefficient().value.value());
    });
}

#[test]
fn source_mismatches_invalid_budgets_and_domain_exits_refuse() {
    let card = card(true, true, 1.0);
    let other = self::card(false, true, 1.0);
    let law = law(&card);
    let source = caloric(&card);
    with_cx(|cx| {
        assert!(LiquidHeatPath::try_new(cx, &law, &caloric(&other), MassFlowRate::new(0.001), Temperature::new(300.0),
            &[section("wall")], LiquidPropertyTemperature::BulkMean).is_err());
        assert!(LiquidHeatPath::try_new(cx, &law, &source, MassFlowRate::new(0.001), Temperature::new(300.0),
            &[section("same"), section("same")], LiquidPropertyTemperature::BulkMean).is_err());
        let path = LiquidHeatPath::try_new(cx, &law, &source, MassFlowRate::new(0.001), Temperature::new(300.0),
            &[section("wall")], LiquidPropertyTemperature::BulkMean).unwrap();
        for wall in [0.0, f64::NAN, 249.0, 451.0] { assert!(path.march(cx, &[wall], budget()).is_err()); }
        assert!(path.march(cx, &[350.0, 360.0], budget()).is_err());
        assert!(path.march(cx, &[350.0], LiquidFeedbackBudget { max_region_iterations: 0, ..budget() }).is_err());
        assert!(path.march(cx, &[350.0], LiquidFeedbackBudget { max_region_iterations: 1, ..budget() }).is_err());
        let no_heat = path.march(cx, &[300.0], budget()).unwrap();
        assert_eq!(no_heat.total_heat_rate_w, 0.0);
        assert_eq!(no_heat.regions[0].inlet, no_heat.regions[0].outlet);
    });
}

#[path = "liquid_feedback/coupled.rs"]
mod coupled;
