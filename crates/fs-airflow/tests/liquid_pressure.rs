//! Synthetic source data through production material, heat and pressure owners.
use fs_airflow::graph::thermal::medium::{FlowHeatGeometry};
use fs_airflow::graph::thermal::medium::series::hydraulic::MaterialFlowResistance;
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::{
    LiquidHeatPath, LiquidHeatSection, LiquidPropertyTemperature, LiquidFeedbackBudget,
};
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::pressure::{
    LiquidResistanceSection, PressureDrivenBudget, PressureDrivenError,
};
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_convection::CorrelationId;
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel};
use fs_material::liquid_temperature::{LiquidTemperatureLaw, LiquidViscosityBasis};
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::state_point::MaterialPropertySelection;
use fs_qty::{Area, Dims, Length, MassFlowRate, Pressure, QuantitySpec, Temperature};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};

fn area() -> Area { Area::new(core::f64::consts::PI*0.01*0.01/4.0) }
fn geometry() -> FlowHeatGeometry {
    FlowHeatGeometry { flow_area: area(), characteristic_length: Length::new(0.01),
        heated_length: Length::new(2.0), wetted_area: Area::new(0.003), aspect_ratio: None }
}
fn point(t: f64) -> QueryPoint {
    QueryPoint::new().with("T",t).unwrap().with("pressure",100000.0).unwrap()
}
fn card(variable: bool, viscosity_scale: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name,value,dims) in [
        ("density",1000.0,Dims([-3,1,0,0,0,0])),
        ("specific_heat_capacity",2000.0,Dims([2,0,-2,-1,0,0])),
        ("thermal_conductivity",0.5,Dims([1,1,-3,-1,0,0])),
        ("dynamic_viscosity",0.002*viscosity_scale,Dims([-1,1,-1,0,0,0])),
    ] {
        let varying = name == "dynamic_viscosity" && variable;
        let quantity = if name == "thermal_conductivity" {
            QuantitySpec::semantic(SemanticType::new(QuantityKind::ThermalConductivity,ValueForm::Static))
        } else { QuantitySpec::dimensional(dims) };
        claims.insert_claim(PropertyClaim { key: PropertyKey::with_quantity(name,quantity),
            value: if varying { PropertyValue::Curve { abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
                knots: vec![(250.0,0.0025*viscosity_scale),(450.0,0.0005*viscosity_scale)],dims }
            } else { PropertyValue::Scalar { value,dims } },
            validity: ValidityDomain::unconstrained().with("T",250.0,450.0).with("pressure",100000.0,100000.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: if varying { InterpolationPolicy::LinearInside } else { InterpolationPolicy::ConstantWithinValidity },
            observations: vec![], provenance: Provenance { source: "synthetic hydraulic thermal fixture".into(),
                license: "test-only".into(),artifact: None },
        }).unwrap();
    }
    MaterialCard::assemble(MaterialStateId { chemistry: "synthetic-liquid".into(), phase: "liquid".into(),
        process: "numerical model test".into(),revision: 0 },claims,vec![]).unwrap()
}
fn with_path<R>(variable: bool, viscosity_scale: f64, convention: LiquidPropertyTemperature,
    run: impl FnOnce(&Cx<'_>, &LiquidHeatPath<'_, '_>) -> R) -> R {
    let card = card(variable,viscosity_scale);
    let law = LiquidTemperatureLaw::try_new(&card,&point(300.0),"T",LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly).unwrap();
    let cp = PropertyKey::with_quantity("specific_heat_capacity",QuantitySpec::dimensional(Dims([2,0,-2,-1,0,0])));
    let rho = PropertyKey::with_quantity("density",QuantitySpec::dimensional(Dims([-3,1,0,0,0,0])));
    let source = ResolvedHeatCapacityCurve::try_new(&card,HeatCapacitySourceQuery {
        heat_capacity: &cp,density: &rho,temperature_axis: "T",lower: &point(250.0),upper: &point(450.0),
        selection: SelectionPolicy::SingleClaimOnly,
    },SolidLiquidPhase::Liquid,0.0,0.5).unwrap();
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(&gate,arena,StreamKey { seed: 91,kernel_id: 740,tile: 0,iteration: 0 },
            Budget::INFINITE,ExecMode::Deterministic);
        let sections = [LiquidHeatSection { region: "first",geometry: geometry(),correlation: CorrelationId::CircularDuctLaminarCwt },
            LiquidHeatSection { region: "last",geometry: geometry(),correlation: CorrelationId::CircularDuctLaminarCwt }];
        let path = LiquidHeatPath::try_new(&cx,&law,&source,MassFlowRate::new(0.001),Temperature::new(300.0),
            &sections,convention).unwrap();
        run(&cx,&path)
    })
}
fn resistances() -> [LiquidResistanceSection<'static>; 2] {
    let resistance = MaterialFlowResistance::try_new(5e9,0.0,area(),Length::new(0.01),(0.0,2300.0)).unwrap();
    [LiquidResistanceSection { region: "first",resistance },LiquidResistanceSection { region: "last",resistance }]
}
fn local() -> LiquidFeedbackBudget {
    LiquidFeedbackBudget { temperature_tolerance_k: 1e-9,coefficient_relative_tolerance: 1e-9,max_region_iterations: 64 }
}
fn pressure_budget() -> PressureDrivenBudget {
    PressureDrivenBudget { mass_flow_bracket_kg_per_s: (0.0002,0.008),pressure_absolute_tolerance_pa: 1e-9,
        pressure_relative_tolerance: 1e-9,max_flow_evaluations: 64,maximum_neglected_hydraulic_power_w: 0.01 }
}
fn close(a: f64,b: f64) { assert!((a-b).abs()<3e-8*(1.0+a.abs().max(b.abs())),"{a} != {b}"); }

#[test]
fn pressure_selects_mass_instead_of_reusing_template_flow() {
    for scale in [1.0,2.0] {
        with_path(false,scale,LiquidPropertyTemperature::BulkMean,|cx,path| {
            let result = path.solve_pressure_driven(cx,&[300.0,300.0],&resistances(),Pressure::new(50.0),local(),pressure_budget()).unwrap();
            close(result.mass_flow_kg_per_s,50.0/(20000.0*scale));
            assert_ne!(result.mass_flow_kg_per_s,path.mass_flow_kg_per_s());
            close(result.hydraulics.total_pressure_drop_pa,50.0);
            close(result.thermal.total_heat_rate_w,0.0);
            assert!(result.flow_evaluations<=pressure_budget().max_flow_evaluations);
        });
    }
}
#[test]
fn heating_changes_viscosity_flow_and_the_downstream_heat_trajectory_together() {
    with_path(true,1.0,LiquidPropertyTemperature::BulkMean,|cx,path| {
        let solve = |walls: &[f64]| path.solve_pressure_driven(cx,walls,&resistances(),Pressure::new(50.0),local(),pressure_budget()).unwrap();
        let cold = solve(&[300.0,300.0]);
        let hot = solve(&[400.0,400.0]);
        assert!(hot.mass_flow_kg_per_s>cold.mass_flow_kg_per_s);
        close(hot.hydraulics.total_pressure_drop_pa,50.0);
        close(hot.thermal.regions[1].inlet.temperature_k(),hot.thermal.regions[0].outlet.temperature_k());
        close(hot.thermal.total_heat_rate_w,hot.mass_flow_kg_per_s*2000.0*(hot.thermal.regions[1].outlet.temperature_k()-300.0));
        for r in &hot.hydraulics.regions {
            close(r.loss.mass_flow_kg_per_s,hot.mass_flow_kg_per_s);
            assert!(r.material.dynamic_viscosity_pa_s()<0.002);
        }
        let again = solve(&[400.0,400.0]);
        assert_eq!(hot.mass_flow_kg_per_s,again.mass_flow_kg_per_s);
        assert_eq!(hot.thermal,again.thermal);
        assert_eq!(hot.hydraulics,again.hydraulics);
    });
}
#[test]
fn hydraulic_material_uses_bulk_not_heat_transfer_film_temperature() {
    with_path(true,1.0,LiquidPropertyTemperature::FilmMean,|cx,path| {
        let result = path.solve_pressure_driven(cx,&[400.0,380.0],&resistances(),Pressure::new(50.0),local(),pressure_budget()).unwrap();
        for (loss,heat) in result.hydraulics.regions.iter().zip(&result.thermal.regions) {
            let bulk = 0.5*(heat.inlet.temperature_k()+heat.outlet.temperature_k());
            close(loss.material.dynamic_viscosity_pa_s(),0.005-1e-5*bulk);
            assert!(loss.material.dynamic_viscosity_pa_s()>heat.material.dynamic_viscosity_pa_s());
            let hydraulic_receipt = &loss.material.material().property("dynamic_viscosity").unwrap().answer().receipt;
            let heat_receipt = &heat.material.material().property("dynamic_viscosity").unwrap().answer().receipt;
            assert_eq!(hydraulic_receipt.selected,heat_receipt.selected);
            assert!(hydraulic_receipt.query_point.contains(&("pressure".into(),100000.0)));
        }
    });
}
#[test]
fn unsupported_bracket_geometry_and_neglected_heat_are_not_silently_accepted() {
    with_path(false,1.0,LiquidPropertyTemperature::BulkMean,|cx,path| {
        let solve = |r: &[LiquidResistanceSection<'_>],p,b| path.solve_pressure_driven(cx,&[350.0,340.0],r,Pressure::new(p),local(),b);
        assert!(matches!(solve(&resistances(),100000.0,pressure_budget()),Err(PressureDrivenError::NotBracketed { .. })));
        let mut reversed = resistances(); reversed.swap(0,1);
        assert!(matches!(solve(&reversed,50.0,pressure_budget()),Err(PressureDrivenError::Invalid(_))));
        assert!(matches!(solve(&resistances(),50.0,PressureDrivenBudget { maximum_neglected_hydraulic_power_w: 0.0,
            ..pressure_budget() }),Err(PressureDrivenError::NeglectedDissipation { .. })));
        assert!(matches!(solve(&resistances(),50.0,PressureDrivenBudget { max_flow_evaluations: 2,
            ..pressure_budget() }),Err(PressureDrivenError::NotConverged { evaluations: 2,.. })));
        assert!(solve(&resistances(),-1.0,pressure_budget()).is_err());
    });
}

fn solid_controls() -> fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::pressure::solid::PressureSolidControls {
    use fs_airflow::conjugate::{ConjugateConfig,Relaxation};
    use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::pressure::solid::PressureSolidControls;
    PressureSolidControls { pressure: pressure_budget(),local: local(),acceleration: Default::default(),max_total_solid_solves: 1024,
        coupling: ConjugateConfig { temperature_tolerance_k: 1e-7,max_iterations: 64,balance_tolerance_w: 1e-5,
            balance_relative_tolerance: 1e-7,relaxation: Relaxation::Fixed { omega: 1.0 } } }
}

#[test]
fn pressure_and_temperature_close_around_the_actual_fem_solid_solver() {
    use fs_conduction::{ConductionMesh,ConductionProblem,ConductivityModel,InitialGuess,LinearConfig,
        Nonlinearity,ScalarField,SolveConfig,StopRule,ThermalBoundaryBuilder,solve};
    use fs_conduction::fixtures::{box_grid,on_box_face};
    use fs_rep_mesh::TetComplex;
    use fs_airflow::conjugate::SolidRegionState;
    use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::LiquidSolidResponse;
    with_path(true,1.0,LiquidPropertyTemperature::BulkMean,|cx,path| {
        let (topology,positions) = box_grid([2,1,1],[0.06,0.1,0.01]);
        let mesh = ConductionMesh::new(TetComplex::from_tets(positions.len(),topology.tets),positions).unwrap();
        let material = ConductivityModel::isotropic_declared(5.0).unwrap();
        let source = ScalarField::nodal("uniform 20W",mesh.vertex_count(),vec![20.0/(0.06*0.1*0.01);mesh.vertex_count()]).unwrap();
        let mut callbacks = 0;
        let mut controls = solid_controls();
        // Pressure closure must not request more accuracy than the inner FEM/thermal solve.
        controls.pressure.pressure_relative_tolerance = 1e-7;
        let answer = path.solve_pressure_coupled_iqn(cx,&[325.0,330.0],&resistances(),Pressure::new(50.0),&controls,|cx,rows| {
            callbacks += 1;
            let mut builder = ThermalBoundaryBuilder::new(&mesh);
            for (index,name) in ["first","last"].iter().enumerate() {
                builder = builder.region(name,|face| face.outward_normal[2]>0.99 && on_box_face(face.centroid[2],0.01)
                    && if index == 0 { face.centroid[0]<0.03 } else { face.centroid[0]>=0.03 },
                    rows[index].boundary_condition().clone()).unwrap();
            }
            let boundary = builder.adiabatic_remainder().finish().unwrap();
            let solution = solve(cx,ConductionProblem { mesh: &mesh,boundary: &boundary,material: &material,
                source: &source,element_materials: None },SolveConfig {
                nonlinearity: Nonlinearity::FixedPoint { relaxation: 1.0,max_backtracks: 8 },
                stop: StopRule { residual_rtol: 1e-10,residual_atol: 1e-18,step_atol: 0.0,max_iterations: 12 },
                linear: LinearConfig { tolerance: 1e-10,max_iterations: 10000,restart: 64 },initial: InitialGuess::Uniform(330.0),
            }).unwrap();
            let regions = solution.report.robin_fluxes.iter().map(SolidRegionState::from_robin_flux).collect();
            Ok(LiquidSolidResponse { regions,state: solution })
        }).unwrap();
        assert_eq!(callbacks,answer.total_solid_solves);
        assert!(callbacks<=controls.max_total_solid_solves);
        close(answer.operating_point.hydraulics.total_pressure_drop_pa,50.0);
        let result = answer.operating_point.thermal;
        assert!((result.fluid.total_heat_rate_w-20.0).abs()<2e-5);
        assert!((result.solid.state.report.energy.robin_out_w-20.0).abs()<2e-5);
        assert!(result.maximum_heat_imbalance_w<1e-5);
        for (row,region) in result.boundary_state.regions.iter().zip(&result.solid.regions) {
            close(row.boundary.reference_temperature().value(),region.mean_reference_temperature_k.unwrap());
        }
    });
}

#[test]
fn shared_solid_budget_is_not_reset_by_outer_pressure_trials() {
    use fs_airflow::conjugate::SolidRegionState;
    use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::LiquidSolidResponse;
    with_path(false,1.0,LiquidPropertyTemperature::BulkMean,|cx,path| {
        let mut controls = solid_controls(); controls.max_total_solid_solves = 1;
        let mut called = 0;
        let result = path.solve_pressure_coupled_iqn(cx,&[300.0,300.0],&resistances(),Pressure::new(50.0),&controls,|_,rows| {
            called += 1;
            Ok(LiquidSolidResponse { state: called,regions: rows.iter().enumerate().map(|(i,r)| SolidRegionState {
                region: ["first","last"][i].into(),area_m2: 0.003,mean_wall_temperature_k: 300.0,
                heat_rate_w: 0.0,mean_reference_temperature_k: Some(r.reference_temperature().value()),
            }).collect() })
        });
        // The first isothermal flow trial converges in one solid call. The
        // SECOND pressure trial must refuse before making another callback.
        assert_eq!(called,1);
        assert!(matches!(result,Err(PressureDrivenError::SolidSolveBudget { performed: 1,maximum: 1 })));
    });
}

#[test]
fn malformed_hydraulic_assignment_and_solid_reply_cannot_publish_a_solution() {
    use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::LiquidSolidResponse;
    with_path(false,1.0,LiquidPropertyTemperature::BulkMean,|cx,path| {
        let mut called = 0;
        let mut bad = resistances(); bad[0].region = "foreign-wall";
        let result = path.solve_pressure_coupled_iqn(cx,&[300.0,300.0],&bad,Pressure::new(50.0),&solid_controls(),|_,_| {
            called += 1; Ok(LiquidSolidResponse { state: (),regions: vec![] })
        });
        assert!(matches!(result,Err(PressureDrivenError::Invalid(_))));
        assert_eq!(called,0);
        let result = path.solve_pressure_coupled_iqn(cx,&[300.0,300.0],&resistances(),Pressure::new(50.0),&solid_controls(),|_,_| {
            called += 1; Ok(LiquidSolidResponse { state: (),regions: vec![] })
        });
        assert!(matches!(result,Err(PressureDrivenError::Thermal(_))));
        assert_eq!(called,1);
    });
}
