//! Synthetic source data through production material, heat and pressure owners.
use fs_airflow::graph::thermal::medium::{FlowHeatGeometry};
use fs_airflow::graph::thermal::medium::series::hydraulic::MaterialFlowResistance;
use fs_airflow::graph::thermal::medium::series::hydraulic::pressure::SeriesGravity;
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::recirculation::{
    GravityReservoirControls, RecirculationError, ReservoirStepControls,
};
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
    card_with_transport_ceiling(variable,viscosity_scale,450.0)
}
fn card_with_transport_ceiling(variable: bool, viscosity_scale: f64, ceiling: f64) -> MaterialCard {
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
            validity: ValidityDomain::unconstrained().with("T",250.0,
                if name == "dynamic_viscosity" || name == "thermal_conductivity" { ceiling } else { 450.0 })
                .with("pressure",100000.0,100000.0),
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
    with_card_path(&card,convention,run)
}
fn with_card_path<R>(card: &MaterialCard, convention: LiquidPropertyTemperature,
    run: impl FnOnce(&Cx<'_>, &LiquidHeatPath<'_, '_>) -> R) -> R {
    with_card_path_gate(card,convention,|cx,path,_| run(cx,path))
}
fn with_card_path_gate<R>(card: &MaterialCard, convention: LiquidPropertyTemperature,
    run: impl FnOnce(&Cx<'_>, &LiquidHeatPath<'_, '_>, &CancelGate) -> R) -> R {
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
        run(&cx,&path,&gate)
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

fn buoyant_card() -> MaterialCard {
    let base = card(false,1.0);
    let mut claims = ClaimSet::new();
    for (_,claim) in base.claims().claims_ordered() {
        let mut claim = claim.clone();
        if claim.key.name() == "density" {
            claim.value = PropertyValue::Curve { abscissa: "T".into(),abscissa_dims: Temperature::DIMS,
                knots: vec![(250.0,1025.0),(450.0,925.0)],dims: Dims([-3,1,0,0,0,0]) };
            claim.interpolation = InterpolationPolicy::LinearInside;
        }
        claims.insert_claim(claim).unwrap();
    }
    MaterialCard::assemble(base.id().clone(),claims,vec![]).unwrap()
}
#[test]
fn source_density_closes_zero_pump_buoyant_flow_and_the_original_heat_path() {
    let card = buoyant_card();
    with_card_path(&card,LiquidPropertyTemperature::FilmMean,|cx,path| {
        // First leg descends while heating; the warmer second leg ascends.
        // Returning to the 300 K inlet represents a prescribed thermal reservoir.
        let gravity = SeriesGravity::try_new([0.0,0.0,-9.81],
            vec![[0.0,0.0,1.0],[0.0;3],[0.0,0.0,1.0]]).unwrap();
        let result = path.solve_pressure_gravity_driven(cx,&[400.0,300.0],&resistances(),
            Pressure::new(0.0),&gravity,local(),pressure_budget()).unwrap();
        // Independent constant-Cp, constant-HTC path and scalar pressure balance.
        // Density and resistance are evaluated at each bulk mean, never the film.
        let reference = |mass:f64| {
            let e = (-3.66*0.5/0.01*0.003/(mass*2000.0)).exp();
            let t1 = 400.0-100.0*e;
            let t2 = 300.0+(t1-300.0)*e;
            let r1 = 1000.0-0.5*(0.5*(300.0+t1)-300.0);
            let r2 = 1000.0-0.5*(0.5*(t1+t2)-300.0);
            (1e7*mass*(1.0/r1+1.0/r2)-9.81*(r1-r2),t2)
        };
        let (mut lo,mut hi) = pressure_budget().mass_flow_bracket_kg_per_s;
        for _ in 0..80 {
            let mid = 0.5*(lo+hi);
            if reference(mid).0<0.0 { lo=mid; } else { hi=mid; }
        }
        let mass = 0.5*(lo+hi);
        close(result.mass_flow_kg_per_s,mass);
        close(result.thermal.regions[1].outlet.temperature_k(),reference(mass).1);
        close(result.thermal.total_heat_rate_w,mass*2000.0*(reference(mass).1-300.0));
        assert!(result.thermal.total_heat_rate_w>0.0);
        assert_eq!(result.driving_pressure_pa,0.0);
        assert_eq!(result.gravity.as_ref(),Some(&gravity));
        let hydrostatic = result.hydrostatic.as_ref().unwrap();
        close(result.hydraulics.total_pressure_drop_pa,hydrostatic.pressure_gain_pa());
        close(result.pressure_residual_pa,0.0);
        assert_eq!(hydrostatic.gravitational_power_w(),0.0);
        for (hydraulic,thermal) in result.hydraulics.regions.iter().zip(&result.thermal.regions) {
            let bulk = 0.5*(thermal.inlet.temperature_k()+thermal.outlet.temperature_k());
            close(hydraulic.material.density_kg_m3(),1000.0-0.5*(bulk-300.0));
            assert_ne!(hydraulic.material.density_kg_m3(),thermal.material.density_kg_m3());
        }
        let no_contrast = path.solve_pressure_gravity_driven(cx,&[300.0,300.0],&resistances(),
            Pressure::new(0.0),&gravity,local(),pressure_budget());
        assert!(matches!(no_contrast,Err(PressureDrivenError::NotBracketed { .. })));
        let reverse = SeriesGravity::try_new([0.0,0.0,9.81],gravity.vertices_m().to_vec()).unwrap();
        assert!(matches!(path.solve_pressure_gravity_driven(cx,&[400.0,300.0],&resistances(),
            Pressure::new(0.0),&reverse,local(),pressure_budget()),Err(PressureDrivenError::NotBracketed { .. })));
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

fn heated_resistances(factor: f64) -> [LiquidResistanceSection<'static>;2] {
    let resistance = MaterialFlowResistance::try_new(factor,0.0,area(),Length::new(0.01),(0.0,2300.0)).unwrap();
    [LiquidResistanceSection { region: "first",resistance }, LiquidResistanceSection { region: "last",resistance }]
}

// Independent closed form for constant Cp/HTC and affine bulk viscosity.
// Eliminating mu(T_bulk) from the exponential regional solution gives one
// linear equation in each outlet. No production thermal or pressure solve.
fn heated_reference(mass: f64, variable: bool) -> (f64,f64,f64) {
    let ua = 3.66*0.5/0.01*0.003;
    let e = (-ua/(mass*2000.0)).exp();
    let f = 1.0-e;
    let q = mass/1000.0;
    let c = 5e15*q*q;
    let (a,b) = if variable { (0.005,1e-5) } else { (0.002,0.0) };
    let (mut inlet,mut pressure,mut power) = (300.0,0.0,0.0);
    for _ in 0..2 {
        let outlet = (e*inlet+f*300.0+f*c/ua*(a-0.5*b*inlet))/(1.0+0.5*f*c*b/ua);
        let viscosity = a-0.5*b*(inlet+outlet);
        pressure += 5e15*viscosity*q;
        power += c*viscosity;
        inlet = outlet;
    }
    (inlet,pressure,power)
}

#[test]
fn hydraulic_heat_and_temperature_dependent_viscosity_close_together() {
    for variable in [false,true] {
        with_path(variable,1.0,LiquidPropertyTemperature::BulkMean,|cx,path| {
            let report = path.march_with_dissipation(cx,&[300.0,300.0],&heated_resistances(5e15),local()).unwrap();
            let (outlet,pressure,power) = heated_reference(0.001,variable);
            close(report.thermal.regions[1].outlet.temperature_k(),outlet);
            close(report.hydraulics.total_pressure_drop_pa,pressure);
            close(report.hydraulics.total_dissipated_power_w,power);
            close(report.enthalpy_gain_w,2.0*(outlet-300.0));
            close(report.enthalpy_gain_w,report.thermal.total_heat_rate_w+power);
            close(report.balance_residual_w,0.0);
            assert!(report.thermal.total_heat_rate_w<0.0);
            for (region,loss) in report.thermal.regions.iter().zip(&report.hydraulics.regions) {
                let wall = region.boundary.coefficient().value.value()*0.003
                    *(300.0-region.boundary.reference_temperature().value());
                close(wall,region.heat_rate_w);
                close(region.heat_rate_w+loss.loss.dissipated_power_w,0.001*
                    (region.outlet.specific_enthalpy_j_kg()-region.inlet.specific_enthalpy_j_kg()));
            }
        });
    }
}

#[test]
fn pressure_driven_heating_matches_an_independent_coupled_operating_point() {
    with_path(true,1.0,LiquidPropertyTemperature::BulkMean,|cx,path| {
        let mut lo = 0.0008;
        let mut hi = 0.0013;
        for _ in 0..80 {
            let mass = 0.5*(lo+hi);
            if heated_reference(mass,true).1 < 2e7 { lo=mass; } else { hi=mass; }
        }
        let expected_mass = 0.5*(lo+hi);
        let expected = heated_reference(expected_mass,true);
        let report = path.solve_pressure_driven_with_dissipation(cx,&[300.0,300.0],
            &heated_resistances(5e15),Pressure::new(2e7),local(),PressureDrivenBudget {
                mass_flow_bracket_kg_per_s: (0.0008,0.0013), maximum_neglected_hydraulic_power_w: 0.0,
                ..pressure_budget()
            }).unwrap();
        assert!((report.mass_flow_kg_per_s/expected_mass-1.0).abs()<2e-8);
        assert!(report.mass_flow_kg_per_s>0.001);
        close(report.thermal.regions[1].outlet.temperature_k(),expected.0);
        close(report.hydraulics.total_pressure_drop_pa,2e7);
        close(report.hydraulics.total_dissipated_power_w,expected.2);
        close(report.thermal.total_heat_rate_w+expected.2,
            report.mass_flow_kg_per_s*2000.0*(expected.0-300.0));
        assert!(path.solve_pressure_driven_with_dissipation(cx,&[300.0,300.0],
            &heated_resistances(5e15),Pressure::new(2e7),local(),pressure_budget()).is_err());
    });
}

#[test]
fn heated_outlet_prediction_stays_inside_narrow_transport_support() {
    let card = card_with_transport_ceiling(false,1.0,340.0);
    for convention in [LiquidPropertyTemperature::BulkMean,LiquidPropertyTemperature::FilmMean] {
        with_card_path(&card,convention,|cx,path| {
            let budget = LiquidFeedbackBudget { max_region_iterations: 1,..local() };
            let result = path.march_with_dissipation(cx,&[300.0,300.0],&heated_resistances(5e15),budget).unwrap();
            close(result.thermal.regions[1].outlet.temperature_k(),heated_reference(0.001,false).0);
            // Cp reaches450 K, but actual property queries above340 K still
            // fail through the original source owner, without extrapolation.
            assert!(matches!(path.march_with_dissipation(cx,&[300.0,300.0],&heated_resistances(5e16),budget),
                Err(PressureDrivenError::Thermal(
                    fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::LiquidFeedbackError::Material(_)))));
        });
    }
}

mod porous;
mod direction;
#[path = "liquid_pressure/transient.rs"]
mod transient;

mod dissipative_solid;

fn gravity_reservoir_controls() -> GravityReservoirControls {
    GravityReservoirControls {
        local: local(), pressure: pressure_budget(), max_total_flow_evaluations: 4096,
        reservoir: ReservoirStepControls {
            temperature_bracket_k: (300.0,315.0), energy_absolute_tolerance_j: 1e-6,
            max_inlet_evaluations: 64, maximum_neglected_hydraulic_energy_j: 0.1,
        },
    }
}
fn closed_gravity() -> SeriesGravity {
    SeriesGravity::try_new([0.0,0.0,-9.81],
        vec![[0.0,0.0,1.0],[0.0;3],[0.0,0.0,1.0]]).unwrap()
}

// Independent constant-Cp/constant-HTC equations with the source's linear rho(T).
// Solves only scalar algebra, never calling the production heat/pressure owners.
fn thermosiphon_reference(inlet: f64) -> (f64,f64) {
    let trial = |mass:f64| {
        let e = (-0.549/(mass*2000.0)).exp();
        let middle = 400.0+(inlet-400.0)*e;
        let outlet = 300.0+(middle-300.0)*e;
        let down = 1000.0-0.5*(0.5*(inlet+middle)-300.0);
        let up = 1000.0-0.5*(0.5*(middle+outlet)-300.0);
        let residual = 1e7*mass*(1.0/down+1.0/up)-9.81*(down-up);
        (residual,mass*2000.0*(outlet-inlet))
    };
    let (mut lo,mut hi) = pressure_budget().mass_flow_bracket_kg_per_s;
    assert!(trial(lo).0<0.0 && trial(hi).0>0.0);
    for _ in 0..80 {
        let mid = 0.5*(lo+hi);
        if trial(mid).0<0.0 { lo=mid; } else { hi=mid; }
    }
    let mass = 0.5*(lo+hi);
    (mass,trial(mass).1)
}

#[test]
fn finite_thermosiphon_inventory_closes_energy_and_recomputes_its_circulation() {
    with_card_path(&buoyant_card(),LiquidPropertyTemperature::FilmMean,|cx,path| {
        let gravity = closed_gravity();
        let controls = gravity_reservoir_controls();
        let mut tank = path.reservoir(cx,0.1,Temperature::new(300.0),0.0).unwrap();
        let original_source = tank.source();
        let mut previous_flow = thermosiphon_reference(300.0).0;
        for (end,external) in [(10.0,0.0),(20.0,-100.0),(30.0,100.0)] {
            let initial = tank.clone();
            let duration = end-initial.time_s();
            let mut lo = controls.reservoir.temperature_bracket_k.0;
            let mut hi = controls.reservoir.temperature_bracket_k.1;
            for _ in 0..80 {
                let mid = 0.5*(lo+hi);
                let residual = 200.0*(mid-initial.state().temperature_k())
                    -duration*thermosiphon_reference(mid).1-external;
                if residual<0.0 { lo=mid; } else { hi=mid; }
            }
            let expected = 0.5*(lo+hi);
            let answer = path.advance_reservoir_with_gravity(cx,&mut tank,end,external,
                &[400.0,300.0],&resistances(),Pressure::new(0.0),&gravity,controls).unwrap();
            close(tank.state().temperature_k(),expected);
            close(answer.circuit.mass_flow_kg_per_s,thermosiphon_reference(expected).0);
            close(answer.returned_heat_j,duration*thermosiphon_reference(expected).1);
            assert!(answer.energy_residual_j.abs()<=controls.reservoir.energy_absolute_tolerance_j);
            close(answer.stored_enthalpy_change_j,answer.returned_heat_j+external);
            assert_eq!(answer.initial,initial.state());
            assert_eq!(answer.final_state,tank.state());
            assert_eq!(answer.start_time_s,initial.time_s());
            assert_eq!(answer.end_time_s,end);
            assert_eq!(tank.time_s(),end);
            assert_eq!(tank.mass_kg(),initial.mass_kg());
            assert!(core::ptr::eq(tank.source(),original_source));
            assert!(core::ptr::eq(answer.circuit.thermal.source,original_source));
            assert_eq!(answer.circuit.thermal.regions[0].inlet,tank.state());
            assert_eq!(answer.circuit.driving_pressure_pa,0.0);
            assert_eq!(answer.circuit.hydrostatic.as_ref().unwrap().gravitational_power_w(),0.0);
            close(answer.circuit.hydraulics.total_pressure_drop_pa,
                answer.circuit.hydrostatic.as_ref().unwrap().pressure_gain_pa());
            assert!(answer.circuit.mass_flow_kg_per_s<previous_flow);
            assert!(answer.total_flow_evaluations>answer.circuit.flow_evaluations);
            assert!(answer.total_flow_evaluations<=controls.max_total_flow_evaluations);
            assert_eq!(answer.total_solid_solves,0);
            assert!(answer.neglected_hydraulic_energy_j<=controls.reservoir.maximum_neglected_hydraulic_energy_j);
            previous_flow = answer.circuit.mass_flow_kg_per_s;
            // Repeating from the same committed history reproduces the entire state.
            let mut replay = initial;
            let replayed = path.advance_reservoir_with_gravity(cx,&mut replay,end,external,
                &[400.0,300.0],&resistances(),Pressure::new(0.0),&gravity,controls).unwrap();
            assert_eq!(tank,replay);
            assert_eq!(answer.circuit.thermal,replayed.circuit.thermal);
            assert_eq!(answer.total_flow_evaluations,replayed.total_flow_evaluations);
        }
        assert_eq!(path.inlet().temperature_k(),300.0);
    });
}

#[test]
fn thermosiphon_step_shares_actual_pressure_work_and_rolls_back_failures() {
    with_card_path_gate(&buoyant_card(),LiquidPropertyTemperature::BulkMean,|cx,path,gate| {
        let gravity = closed_gravity();
        let controls = gravity_reservoir_controls();
        let mut tank = path.reservoir(cx,0.1,Temperature::new(300.0),0.0).unwrap();
        let original = tank.clone();
        let first = path.solve_pressure_gravity_driven(cx,&[400.0,300.0],&resistances(),
            Pressure::new(0.0),&gravity,controls.local,controls.pressure).unwrap();
        let limited = GravityReservoirControls { max_total_flow_evaluations: first.flow_evaluations, ..controls };
        assert!(matches!(path.advance_reservoir_with_gravity(cx,&mut tank,10.0,0.0,
            &[400.0,300.0],&resistances(),Pressure::new(0.0),&gravity,limited),
            Err(RecirculationError::WorkBudget { flow_evaluations, solid_solves: 0 })
                if flow_evaluations == first.flow_evaluations));
        assert_eq!(tank,original);

        // Same elevation alone is insufficient: a reservoir closes the full path.
        let open = SeriesGravity::try_new([0.0,0.0,-9.81],
            vec![[0.0,0.0,1.0],[0.0;3],[1.0,0.0,1.0]]).unwrap();
        assert!(matches!(path.advance_reservoir_with_gravity(cx,&mut tank,10.0,0.0,
            &[400.0,300.0],&resistances(),Pressure::new(0.0),&open,controls),
            Err(RecirculationError::Invalid(_))));
        assert_eq!(tank,original);
        for step in [
            ReservoirStepControls { temperature_bracket_k: (300.0,301.0), ..controls.reservoir },
            ReservoirStepControls { temperature_bracket_k: (100.0,315.0), ..controls.reservoir },
            ReservoirStepControls { maximum_neglected_hydraulic_energy_j: 0.0, ..controls.reservoir },
            ReservoirStepControls { max_inlet_evaluations: 2, ..controls.reservoir },
        ] {
            assert!(path.advance_reservoir_with_gravity(cx,&mut tank,10.0,0.0,
                &[400.0,300.0],&resistances(),Pressure::new(0.0),&gravity,
                GravityReservoirControls { reservoir: step, ..controls }).is_err());
            assert_eq!(tank,original);
        }
        gate.request();
        assert!(path.advance_reservoir_with_gravity(cx,&mut tank,10.0,0.0,
            &[400.0,300.0],&resistances(),Pressure::new(0.0),&gravity,controls).is_err());
        assert_eq!(tank,original);
    });
}

#[test]
fn finite_inventory_does_not_invent_circulation_without_density_contrast() {
    with_path(false,1.0,LiquidPropertyTemperature::BulkMean,|cx,path| {
        let mut tank = path.reservoir(cx,0.1,Temperature::new(300.0),0.0).unwrap();
        let original = tank.clone();
        assert!(matches!(path.advance_reservoir_with_gravity(cx,&mut tank,10.0,0.0,
            &[400.0,300.0],&resistances(),Pressure::new(0.0),&closed_gravity(),
            gravity_reservoir_controls()),Err(RecirculationError::Pressure(
                PressureDrivenError::NotBracketed { .. }))));
        assert_eq!(tank,original);
    });
}
