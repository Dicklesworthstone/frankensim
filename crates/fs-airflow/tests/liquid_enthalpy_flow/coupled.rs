//! The solid is a two-node conducting body with independently deposited power.
use super::*;
use fs_airflow::{AirflowError, conjugate::{ConjugateConfig, Relaxation, SolidRegionState, solve_conjugate_iqn}};
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::LiquidSolidResponse;
use fs_convection::CorrelationRobinBoundary;
use fs_couple::iqn_ils::IqnIlsConfig;

const AREAS: [f64;2] = [0.008,0.006];
const NAMES: [&str;2] = ["upstream","downstream"];

fn config() -> ConjugateConfig {
    ConjugateConfig { temperature_tolerance_k: 1e-10, balance_tolerance_w: 1e-8,
        balance_relative_tolerance: 1e-10, max_iterations: 128,
        relaxation: Relaxation::Fixed { omega: 0.5 }, ..ConjugateConfig::default() }
}
fn solid(rows: &[CorrelationRobinBoundary], power: [f64;2]) -> LiquidSolidResponse<[f64;2]> {
    let ua = [rows[0].coefficient().value.value()*AREAS[0],rows[1].coefficient().value.value()*AREAS[1]];
    let refs = [rows[0].reference_temperature().value(),rows[1].reference_temperature().value()];
    let g = 0.7;
    let diag = [ua[0]+g,ua[1]+g];
    let rhs = [power[0]+ua[0]*refs[0],power[1]+ua[1]*refs[1]];
    let determinant = diag[0]*diag[1]-g*g;
    let walls = [(diag[1]*rhs[0]+g*rhs[1])/determinant,(g*rhs[0]+diag[0]*rhs[1])/determinant];
    let regions = (0..2).map(|i| SolidRegionState {
        region: NAMES[i].into(), area_m2: AREAS[i], mean_wall_temperature_k: walls[i],
        heat_rate_w: ua[i]*(walls[i]-refs[i]), mean_reference_temperature_k: Some(refs[i]),
    }).collect();
    LiquidSolidResponse { regions, state: walls }
}
fn sections() -> [FlowHeatSection<'static>;2] {
    [section(NAMES[0],AREAS[0]),section(NAMES[1],AREAS[1])]
}

#[test]
fn coupled_constant_cp_matches_the_existing_iqn_solid_air_path() {
    let card = card(false,"coupled-constant");
    let liquid = liquid(&card,1.0);
    let source = caloric(&card,1.0,SolidLiquidPhase::Liquid);
    let path = FrozenFlowMedium::Liquid(&liquid).channel_series(VolumetricFlowRate::new(1e-6),
        Temperature::new(300.0),&sections()).unwrap();
    with_cx(|cx| {
        let mut calls = 0;
        let actual = path.solve_liquid_enthalpy_iqn(cx,&source,"T",&config(),IqnIlsConfig::default(),|_,rows| {
            calls+=1; Ok(solid(rows,[30.0,10.0]))
        }).unwrap();
        assert_eq!(calls,actual.iterations);
        let legacy = solve_conjugate_iqn(cx,path.path(),&config(),IqnIlsConfig::default(),|_,refs| {
            let rows:Vec<_> = path.exchanges().iter().zip(refs).map(|(e,&t)| e.robin_boundary(Temperature::new(t)).unwrap()).collect();
            Ok(solid(&rows,[30.0,10.0]).regions)
        }).unwrap();
        close(actual.fluid.outlet().temperature_k(),legacy.march.outlet_temperature_k);
        close(actual.fluid.total_heat_rate_w(),40.0);
        for (a,b) in actual.solid.regions.iter().zip(&legacy.solid) { close(a.mean_wall_temperature_k,b.mean_wall_temperature_k); }
        assert_eq!(actual.solid.state,solid(&actual.boundary_rows,[30.0,10.0]).state);
        assert!(actual.maximum_reference_residual_k<=config().temperature_tolerance_k);
        assert!(actual.maximum_heat_imbalance_w<config().balance_tolerance_w);
    });
}

#[test]
fn coupled_temperature_dependent_coolant_closes_the_original_solid_power() {
    let card = card(true,"coupled-varying");
    let liquid = liquid(&card,1.0);
    let source = caloric(&card,1.0,SolidLiquidPhase::Liquid);
    let path = FrozenFlowMedium::Liquid(&liquid).channel_series(VolumetricFlowRate::new(1e-6),
        Temperature::new(300.0),&sections()).unwrap();
    for sign in [1.0,-1.0] {
        with_cx(|cx| {
            let powers = [30.0*sign,10.0*sign];
            let actual = path.solve_liquid_enthalpy_iqn(cx,&source,"T",&config(),IqnIlsConfig::default(),|_,rows| Ok(solid(rows,powers))).unwrap();
            close(actual.fluid.total_heat_rate_w(),powers.iter().sum());
            // Independent integrated source law: h(T)=1000*(T-250)+5*(T-250)^2.
            let h_out = 62500.0+sign*40.0/path.mass_flow_kg_per_s();
            let expected = 250.0+((1e6+20.0*h_out).sqrt()-1000.0)/10.0;
            assert!((actual.fluid.outlet().temperature_k()-expected).abs()<0.002);
            assert_eq!(actual.solid.state,solid(&actual.boundary_rows,powers).state);
            assert_eq!(actual.fluid.source(),&source);
            for (body,fluid) in actual.solid.regions.iter().zip(actual.fluid.regions()) {
                close(body.heat_rate_w,fluid.heat_rate_w());
                assert_eq!(body.mean_wall_temperature_k,actual.solid.state[if body.region==NAMES[0] {0}else{1}]);
            }
        });
    }
}

#[test]
fn opposite_regional_heat_errors_cannot_cancel_into_acceptance() {
    let card = card(true,"opposing-errors");
    let liquid = liquid(&card,1.0);
    let source = caloric(&card,1.0,SolidLiquidPhase::Liquid);
    let path = FrozenFlowMedium::Liquid(&liquid).channel_series(VolumetricFlowRate::new(1e-6),
        Temperature::new(300.0),&sections()).unwrap();
    with_cx(|cx| {
        let result = path.solve_liquid_enthalpy_iqn(cx,&source,"T",&config(),IqnIlsConfig::default(),|_,rows| {
            let mut response = solid(rows,[30.0,10.0]);
            response.regions[0].heat_rate_w+=1.0;
            response.regions[1].heat_rate_w-=1.0;
            Ok(response)
        });
        assert!(matches!(result,Err(LiquidEnthalpyError::Flow(fs_airflow::graph::thermal::medium::MediumFlowError::Airflow(
            AirflowError::ConjugateBalanceUnclosed {..})))));
    });
}

#[test]
fn actual_sent_boundary_references_and_area_are_required_from_the_solid() {
    let card = card(false,"wrong-solid-inputs");
    let liquid = liquid(&card,1.0);
    let source = caloric(&card,1.0,SolidLiquidPhase::Liquid);
    let path = FrozenFlowMedium::Liquid(&liquid).channel_series(VolumetricFlowRate::new(1e-6),
        Temperature::new(300.0),&sections()).unwrap();
    with_cx(|cx| {
        for mutation in 0..4 {
            let mut calls = 0;
            let result = path.solve_liquid_enthalpy_iqn(cx,&source,"T",&config(),IqnIlsConfig::default(),|_,rows| {
                calls+=1;
                let mut response = solid(rows,[30.0,10.0]);
                match mutation {
                    0 => response.regions[1].mean_reference_temperature_k=None,
                    1 => response.regions[1].mean_reference_temperature_k=Some(200.0),
                    2 => response.regions[1].area_m2*=2.0,
                    _ => response.regions.swap(0,1),
                }
                Ok(response)
            });
            assert!(result.is_err()); assert_eq!(calls,1);
        }
    });
}

#[test]
fn malformed_policy_runs_no_solid_and_iteration_exhaustion_publishes_no_solution() {
    let card = card(false,"limits");
    let liquid = liquid(&card,1.0);
    let source = caloric(&card,1.0,SolidLiquidPhase::Liquid);
    let path = FrozenFlowMedium::Liquid(&liquid).channel_series(VolumetricFlowRate::new(1e-6),
        Temperature::new(300.0),&sections()).unwrap();
    with_cx(|cx| {
        let mut calls = 0;
        for (cfg,iqn) in [(ConjugateConfig {max_iterations:0,..config()},IqnIlsConfig::default()),
            (config(),IqnIlsConfig {max_history:0,..IqnIlsConfig::default()})] {
            assert!(path.solve_liquid_enthalpy_iqn(cx,&source,"T",&cfg,iqn,|_,rows| {
                calls+=1; Ok(solid(rows,[30.0,10.0]))
            }).is_err());
        }
        let outside = FrozenFlowMedium::Liquid(&liquid).channel_series(VolumetricFlowRate::new(1e-6),
            Temperature::new(500.0),&sections()).unwrap();
        assert!(outside.solve_liquid_enthalpy_iqn(cx,&source,"T",&config(),IqnIlsConfig::default(),|_,rows| {
            calls+=1; Ok(solid(rows,[30.0,10.0]))
        }).is_err());
        assert_eq!(calls,0);
        let result = path.solve_liquid_enthalpy_iqn(cx,&source,"T",&ConjugateConfig {max_iterations:1,..config()},
            IqnIlsConfig::default(),|_,rows| {calls+=1;Ok(solid(rows,[30.0,10.0]))});
        assert_eq!(calls,1);
        assert!(matches!(result,Err(LiquidEnthalpyError::Flow(fs_airflow::graph::thermal::medium::MediumFlowError::Airflow(
            AirflowError::ConjugateNotConverged {..})))));
    });
}
