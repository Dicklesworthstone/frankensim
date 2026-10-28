//! Synthetic-source caloric integration, not measured coolant validation.
use fs_airflow::graph::thermal::medium::{FlowHeatGeometry, FlowHeatSection, FrozenFlowMedium};
use fs_airflow::graph::thermal::medium::series::enthalpy::LiquidEnthalpyError;
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_convection::{CorrelationId, ThermalDirection};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::liquid::{ResolvedLiquidState, resolve_liquid_state};
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::state_point::MaterialPropertySelection;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Area, Dims, Length, QuantitySpec, Temperature, VolumetricFlowRate};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(&gate, arena, StreamKey { seed: 31, kernel_id: 735, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic))
    })
}
fn point(t: f64, pressure: f64) -> QueryPoint {
    QueryPoint::new().with("T", t).unwrap().with("pressure", pressure).unwrap()
}
fn card(variable: bool, name: &str) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (key, value, dims) in [
        ("density",1000.0,Dims([-3,1,0,0,0,0])),
        ("specific_heat_capacity",2000.0,Dims([2,0,-2,-1,0,0])),
        ("dynamic_viscosity",0.001,Dims([-1,1,-1,0,0,0])),
        ("thermal_conductivity",0.6,Dims([1,1,-3,-1,0,0])),
    ] {
        let changing = variable && key == "specific_heat_capacity";
        let quantity = if key == "thermal_conductivity" {
            QuantitySpec::semantic(SemanticType::new(QuantityKind::ThermalConductivity, ValueForm::Static))
        } else { QuantitySpec::dimensional(dims) };
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::with_quantity(key, quantity),
            value: if changing {
                PropertyValue::Curve { abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
                    knots: vec![(250.0,1000.0),(350.0,2000.0),(450.0,3000.0)], dims }
            } else { PropertyValue::Scalar { value, dims } },
            validity: ValidityDomain::unconstrained().with("T",250.0,450.0).with("pressure",1.0,2.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: if changing { InterpolationPolicy::LinearInside }
                else { InterpolationPolicy::ConstantWithinValidity },
            observations: vec![],
            provenance: Provenance { source: "synthetic sensible coolant, not a measurement".into(), license: "test-only".into(), artifact: None },
        }).unwrap();
    }
    MaterialCard::assemble(MaterialStateId { chemistry: name.into(), phase: "liquid".into(),
        process: "synthetic fixture".into(), revision: 0 }, claims, vec![]).unwrap()
}
fn liquid(card: &MaterialCard, pressure: f64) -> ResolvedLiquidState {
    resolve_liquid_state(card, &point(300.0,pressure), MaterialPropertySelection::SingleClaimOnly).unwrap()
}
fn caloric(card: &MaterialCard, pressure: f64, phase: SolidLiquidPhase) -> ResolvedHeatCapacityCurve {
    let cp = PropertyKey::with_quantity("specific_heat_capacity", QuantitySpec::dimensional(Dims([2,0,-2,-1,0,0])));
    let rho = PropertyKey::with_quantity("density", QuantitySpec::dimensional(Dims([-3,1,0,0,0,0])));
    ResolvedHeatCapacityCurve::try_new(card, HeatCapacitySourceQuery {
        heat_capacity: &cp, density: &rho, temperature_axis: "T",
        lower: &point(250.0,pressure), upper: &point(450.0,pressure),
        selection: SelectionPolicy::SingleClaimOnly,
    }, phase, 0.0, 0.5).unwrap()
}
fn section(region: &str, area: f64) -> FlowHeatSection<'_> {
    FlowHeatSection { region, geometry: FlowHeatGeometry {
        flow_area: Area::new(1e-4), characteristic_length: Length::new(0.01),
        heated_length: Length::new(16.0), wetted_area: Area::new(area), aspect_ratio: None,
    }, correlation: CorrelationId::CircularDuctLaminarCwt, direction: ThermalDirection::HeatingFluid }
}
fn close(a: f64,b: f64) { assert!((a-b).abs() < 2e-9*(1.0+a.abs().max(b.abs())), "{a} != {b}"); }

#[test]
fn constant_cp_matches_existing_march_and_its_actual_robin_rows() {
    let card = card(false,"constant");
    let liquid = liquid(&card,1.0);
    let source = caloric(&card,1.0,SolidLiquidPhase::Liquid);
    let path = FrozenFlowMedium::Liquid(&liquid).channel_series(VolumetricFlowRate::new(1e-6),
        Temperature::new(300.0), &[section("first",0.004),section("last",0.003)]).unwrap();
    with_cx(|cx| {
        let actual = path.march_liquid_enthalpy(cx,&source,"T",&[370.0,280.0]).unwrap();
        let legacy = path.path().march(&[370.0,280.0]).unwrap();
        close(actual.total_heat_rate_w(),legacy.total_heat_rate_w);
        close(actual.outlet().temperature_k(),legacy.outlet_temperature_k);
        assert_eq!(actual.source(),&source);
        for (actual,expected) in actual.regions().iter().zip(&legacy.segments) {
            close(actual.reference_temperature_k(),expected.reference_temperature_k);
            close(actual.heat_rate_w(),expected.heat_rate_w);
            let boundary = actual.robin_boundary().unwrap();
            match boundary.boundary_condition() {
                fs_conduction::ThermalBc::Robin { htc, t_ref } => {
                    close(htc.at(0),actual.exchange().coefficient().value.value());
                    close(t_ref.at(0),actual.reference_temperature_k());
                }
                other => panic!("unexpected boundary {other:?}"),
            }
        }
    });
}

// Independent integral of cp(T)/(wall-T) for cp(T)=10*T-1500.
fn continuous_cp_reference(inlet: f64, wall: f64, ua_over_mass: f64) -> f64 {
    let mut low = 0.0;
    let mut high = 1.0;
    for _ in 0..80 {
        let fraction = 0.5*(low+high);
        let t = inlet + fraction*(wall-inlet);
        let integral = (10.0*wall-1500.0)*(-(1.0-fraction).ln()) - 10.0*(t-inlet);
        if integral < ua_over_mass { low = fraction; } else { high = fraction; }
    }
    inlet + 0.5*(low+high)*(wall-inlet)
}

#[test]
fn variable_cp_heating_and_cooling_follow_source_enthalpy_not_frozen_cp() {
    let card = card(true,"varying");
    let liquid = liquid(&card,1.0);
    let source = caloric(&card,1.0,SolidLiquidPhase::Liquid);
    assert!(source.curve().knots().len() > 100);
    for (inlet,wall) in [(300.0,420.0),(420.0,270.0)] {
        let path = FrozenFlowMedium::Liquid(&liquid).channel_series(VolumetricFlowRate::new(1e-6),
            Temperature::new(inlet), &[section("wall",0.008)]).unwrap();
        with_cx(|cx| {
            let actual = path.march_liquid_enthalpy(cx,&source,"T",&[wall]).unwrap();
            let ua = path.path().segments()[0].htc_w_per_m2_k()*0.008;
            let expected = continuous_cp_reference(inlet,wall,ua/path.mass_flow_kg_per_s());
            // Includes the declared 0.5 J/kg Cp-to-enthalpy chord approximation.
            assert!((actual.outlet().temperature_k()-expected).abs() < 0.002);
            let legacy = path.path().march(&[wall]).unwrap();
            assert!((actual.outlet().temperature_k()-legacy.outlet_temperature_k).abs() > 0.5);
            close(actual.total_heat_rate_w(),path.mass_flow_kg_per_s()*
                (actual.outlet().specific_enthalpy_j_kg()-actual.regions()[0].inlet().specific_enthalpy_j_kg()));
        });
    }
}

#[test]
fn exact_source_knot_inlet_uses_the_correct_cooling_or_heating_piece() {
    let card = card(true,"knots");
    let liquid = liquid(&card,1.0);
    let source = caloric(&card,1.0,SolidLiquidPhase::Liquid);
    let t = source.curve().knots()[source.curve().knots().len()/2].temperature_k;
    for wall in [280.0,420.0,t] {
        let path = FrozenFlowMedium::Liquid(&liquid).channel_series(VolumetricFlowRate::new(1e-6),
            Temperature::new(t), &[section("wall",0.004)]).unwrap();
        with_cx(|cx| {
            let march = path.march_liquid_enthalpy(cx,&source,"T",&[wall]).unwrap();
            if wall == t { assert_eq!(march.total_heat_rate_w(),0.0); }
            else {
                assert_eq!(march.total_heat_rate_w().is_sign_positive(),wall>t);
                assert!(march.outlet().temperature_k() >= wall.min(t));
                assert!(march.outlet().temperature_k() <= wall.max(t));
            }
        });
    }
}

#[test]
fn source_chart_refinement_of_wall_regions_preserves_energy_and_outlet() {
    let card = card(true,"split");
    let liquid = liquid(&card,1.0);
    let source = caloric(&card,1.0,SolidLiquidPhase::Liquid);
    let medium = FrozenFlowMedium::Liquid(&liquid);
    let whole = medium.channel_series(VolumetricFlowRate::new(1e-6),Temperature::new(300.0),&[section("all",0.008)]).unwrap();
    with_cx(|cx| {
        let expected = whole.march_liquid_enthalpy(cx,&source,"T",&[420.0]).unwrap();
        for count in [2,8] {
            let names:Vec<_> = (0..count).map(|i| format!("wall-{i}")).collect();
            let pieces:Vec<_> = names.iter().map(|n| section(n,0.008/count as f64)).collect();
            let path = medium.channel_series(VolumetricFlowRate::new(1e-6),Temperature::new(300.0),&pieces).unwrap();
            let actual = path.march_liquid_enthalpy(cx,&source,"T",&vec![420.0;count]).unwrap();
            close(actual.outlet().temperature_k(),expected.outlet().temperature_k());
            close(actual.total_heat_rate_w(),expected.total_heat_rate_w());
            for pair in actual.regions().windows(2) { assert_eq!(pair[0].outlet(),pair[1].inlet()); }
        }
    });
}

#[test]
fn every_robin_region_closes_its_own_enthalpy_transfer() {
    let card = card(true,"balance");
    let liquid = liquid(&card,1.0);
    let source = caloric(&card,1.0,SolidLiquidPhase::Liquid);
    let path = FrozenFlowMedium::Liquid(&liquid).channel_series(VolumetricFlowRate::new(1e-6),
        Temperature::new(300.0),&[section("heat",0.008),section("cool",0.012)]).unwrap();
    with_cx(|cx| {
        let march = path.march_liquid_enthalpy(cx,&source,"T",&[420.0,270.0]).unwrap();
        for (region,wall) in march.regions().iter().zip([420.0,270.0]) {
            let boundary = region.robin_boundary().unwrap();
            let q = boundary.coefficient().value.value()*region.exchange().geometry().wetted_area.value()
                *(wall-region.reference_temperature_k());
            close(q,region.heat_rate_w());
        }
        close(march.total_heat_rate_w(),march.mass_flow_kg_per_s()*
            (march.outlet().specific_enthalpy_j_kg()-march.regions()[0].inlet().specific_enthalpy_j_kg()));
        assert!(march.regions()[0].heat_rate_w()>0.0);
        assert!(march.regions()[1].heat_rate_w()<0.0);
        assert_eq!(path.march_liquid_enthalpy(cx,&source,"T",&[420.0,270.0]).unwrap(),march);
    });
}

#[test]
fn mismatched_card_fixed_context_axis_or_phase_refuses() {
    let card = card(true,"one");
    let other_card = self::card(true,"two");
    let liquid = liquid(&card,1.0);
    let path = FrozenFlowMedium::Liquid(&liquid).channel_series(VolumetricFlowRate::new(1e-6),
        Temperature::new(300.0),&[section("wall",0.008)]).unwrap();
    with_cx(|cx| {
        for bad in [caloric(&other_card,1.0,SolidLiquidPhase::Liquid),
            caloric(&card,2.0,SolidLiquidPhase::Liquid),caloric(&card,1.0,SolidLiquidPhase::Solid)] {
            assert!(matches!(path.march_liquid_enthalpy(cx,&bad,"T",&[420.0]),Err(LiquidEnthalpyError::Invalid(_))));
        }
        let good = caloric(&card,1.0,SolidLiquidPhase::Liquid);
        assert!(path.march_liquid_enthalpy(cx,&good,"pressure",&[420.0]).is_err());
        assert!(path.march_liquid_enthalpy(cx,&good,"missing",&[420.0]).is_err());
    });
}

#[test]
fn unsupported_late_fluid_temperature_returns_no_partial_march() {
    let card = card(true,"bounded");
    let liquid = liquid(&card,1.0);
    let source = caloric(&card,1.0,SolidLiquidPhase::Liquid);
    let path = FrozenFlowMedium::Liquid(&liquid).channel_series(VolumetricFlowRate::new(1e-6),
        Temperature::new(300.0),&[section("first",0.004),section("last",1.0)]).unwrap();
    with_cx(|cx| {
        let before = path.march_liquid_enthalpy(cx,&source,"T",&[350.0,320.0]).unwrap();
        assert!(path.march_liquid_enthalpy(cx,&source,"T",&[350.0,500.0]).is_err());
        assert!(path.march_liquid_enthalpy(cx,&source,"T",&[350.0,200.0]).is_err());
        assert!(path.march_liquid_enthalpy(cx,&source,"T",&[350.0,f64::NAN]).is_err());
        assert!(path.march_liquid_enthalpy(cx,&source,"T",&[350.0]).is_err());
        assert_eq!(path.march_liquid_enthalpy(cx,&source,"T",&[350.0,320.0]).unwrap(),before);
        // A wall outside the source does not require extrapolating the fluid.
        let weak = FrozenFlowMedium::Liquid(&liquid).channel_series(VolumetricFlowRate::new(1e-6),
            Temperature::new(300.0),&[section("weak",0.00001)]).unwrap();
        assert!(weak.march_liquid_enthalpy(cx,&source,"T",&[500.0]).is_ok());
    });
}

#[path = "liquid_enthalpy_flow/coupled.rs"]
mod coupled;

fn heating_resistance(factor: f64) -> fs_airflow::graph::thermal::medium::series::hydraulic::MaterialFlowResistance {
    fs_airflow::graph::thermal::medium::series::hydraulic::MaterialFlowResistance::try_new(
        factor, 0.0, Area::new(1e-4), Length::new(0.01), (0.0, 2300.0),
    ).unwrap()
}

#[test]
fn hydraulic_power_heats_the_liquid_while_the_original_wall_receives_rejected_heat() {
    let card = card(false,"dissipative");
    let liquid = liquid(&card,1.0);
    let source = caloric(&card,1.0,SolidLiquidPhase::Liquid);
    let path = FrozenFlowMedium::Liquid(&liquid).channel_series(VolumetricFlowRate::new(1e-6),
        Temperature::new(300.0), &[section("wall",0.008)]).unwrap();
    with_cx(|cx| {
        let result = path.march_liquid_enthalpy_with_dissipation(cx,&source,"T",&[300.0],
            &[heating_resistance(2e16)]).unwrap();
        let ua = path.path().segments()[0].htc_w_per_m2_k()*0.008;
        let expected = 300.0 + 20.0/ua*(1.0-(-ua/2.0).exp());
        close(result.outlet().temperature_k(),expected);
        close(result.hydraulic_heat_rate_w(),20.0);
        close(result.enthalpy_gain_w(),2.0*(expected-300.0));
        close(result.balance_residual_w(),0.0);
        assert!(result.total_heat_rate_w()<0.0);
        let region = &result.regions()[0];
        close(ua*(300.0-region.reference_temperature_k()),result.total_heat_rate_w());
        close(region.hydraulic_loss().unwrap().dissipated_power_w,20.0);
        close(region.hydraulic_heat_rate_w(),20.0);
        let legacy = path.march_liquid_enthalpy(cx,&source,"T",&[300.0]).unwrap();
        assert_eq!(legacy.total_heat_rate_w(),0.0);
        assert_eq!(legacy.hydraulic_heat_rate_w(),0.0);
    });
}

#[test]
fn distributed_hydraulic_heating_preserves_variable_cp_energy_under_subdivision() {
    let card = card(true,"dissipative-refinement");
    let liquid = liquid(&card,1.0);
    let source = caloric(&card,1.0,SolidLiquidPhase::Liquid);
    with_cx(|cx| {
        let mut expected = None;
        for count in [1,2,8] {
            let names:Vec<_> = (0..count).map(|i| format!("heated-{i}")).collect();
            let sections:Vec<_> = names.iter().map(|name| section(name,0.008/count as f64)).collect();
            let path = FrozenFlowMedium::Liquid(&liquid).channel_series(VolumetricFlowRate::new(1e-6),
                Temperature::new(300.0),&sections).unwrap();
            let result = path.march_liquid_enthalpy_with_dissipation(cx,&source,"T",&vec![320.0;count],
                &vec![heating_resistance(2e16/count as f64);count]).unwrap();
            close(result.hydraulic_heat_rate_w(),20.0);
            close(result.enthalpy_gain_w(),result.mass_flow_kg_per_s()*
                (result.outlet().specific_enthalpy_j_kg()-result.regions()[0].inlet().specific_enthalpy_j_kg()));
            close(result.balance_residual_w(),0.0);
            if let Some((temperature,wall_heat)) = expected {
                close(result.outlet().temperature_k(),temperature);
                close(result.total_heat_rate_w(),wall_heat);
            } else { expected = Some((result.outlet().temperature_k(),result.total_heat_rate_w())); }
        }
    });
}

#[test]
fn hydraulic_heating_rejects_geometry_reynolds_and_caloric_domain_mismatches() {
    use fs_airflow::graph::thermal::medium::series::hydraulic::MaterialFlowResistance;
    let card = card(false,"dissipative-refusal");
    let liquid = liquid(&card,1.0);
    let source = caloric(&card,1.0,SolidLiquidPhase::Liquid);
    let path = FrozenFlowMedium::Liquid(&liquid).channel_series(VolumetricFlowRate::new(1e-6),
        Temperature::new(300.0), &[section("wall",0.008)]).unwrap();
    with_cx(|cx| {
        assert!(path.march_liquid_enthalpy_with_dissipation(cx,&source,"T",&[300.0],&[]).is_err());
        for loss in [
            MaterialFlowResistance::try_new(2e16,0.0,Area::new(2e-4),Length::new(0.01),(0.0,2300.0)).unwrap(),
            MaterialFlowResistance::try_new(2e16,0.0,Area::new(1e-4),Length::new(0.01),(0.0,1.0)).unwrap(),
            heating_resistance(1e20),
        ] {
            assert!(path.march_liquid_enthalpy_with_dissipation(cx,&source,"T",&[300.0],&[loss]).is_err());
        }
        assert_eq!(path.path().inlet_temperature_k(),300.0);
    });
}

#[test]
fn advected_energy_remains_independent_of_large_cancelling_heat_reports() {
    let card = card(false,"cancelling-heat");
    let liquid = liquid(&card,1.0);
    let source = caloric(&card,1.0,SolidLiquidPhase::Liquid);
    let loss = heating_resistance(2.0_f64.powi(60)/1e-15);
    let power = loss.evaluate(&liquid,fs_qty::MassFlowRate::new(0.001)).unwrap().dissipated_power_w;
    let htc = 3.66*0.6/0.01;
    let path = FrozenFlowMedium::Liquid(&liquid).channel_series(VolumetricFlowRate::new(1e-6),
        Temperature::new(299.999999), &[section("wall",power/(50.0*htc))]).unwrap();
    with_cx(|cx| {
        let report = path.march_liquid_enthalpy_with_dissipation(cx,&source,"T",&[250.0],&[loss]).unwrap();
        let actual = report.mass_flow_kg_per_s()*(report.outlet().specific_enthalpy_j_kg()
            - report.regions()[0].inlet().specific_enthalpy_j_kg());
        assert!(actual>1e-7);
        assert_eq!(report.total_heat_rate_w()+report.hydraulic_heat_rate_w(),0.0);
        assert_eq!(report.enthalpy_gain_w(),actual);
        assert_eq!(report.balance_residual_w(),actual);
    });
}
