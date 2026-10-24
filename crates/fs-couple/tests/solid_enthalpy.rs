//! Caloric/material data reaches actual 3-D deformation, stress and vibration.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::ContentHash;
use fs_couple::solid_enthalpy::*;
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel};
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve,
    HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::state_point::{MaterialPropertySelection, INVERSE_TEMPERATURE_DIMS,
    SPECIFIC_HEAT_CAPACITY_DIMS, resolve_isotropic_elastic_state_point};
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_qty::{Density, Dims, Pressure, QuantitySpec, Temperature};
use fs_solid::linear3::*;
use fs_solid::loads3::VolumeLoad3;
use fs_thermal::enthalpy::*;

#[path = "solid_enthalpy/directional.rs"]
mod directional;
#[path = "solid_enthalpy/valve.rs"]
mod valve;

const NODES: [[f64;3];4] = [[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]];
const TETS: [[usize;4];1] = [[0,1,2,3]];
const FREE: [usize;6] = [0,1,2,4,5,8];
fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate=CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| f(&Cx::new(&gate,arena,
        StreamKey {seed:3,kernel_id:17,tile:0,iteration:0},Budget::INFINITE,ExecMode::Deterministic)))
}
fn point(t:f64)->QueryPoint {
    QueryPoint::new().with_quantity("T",QuantitySpec::dimensional(Temperature::DIMS),t).unwrap()
        .with("lot",7.0).unwrap()
}
fn material(young:f64,alpha:f64,cp:f64,nu:f64,split:bool)->(MaterialCard,ResolvedHeatCapacityCurve) {
    let heat=PropertyKey::new("specific_heat_capacity",SPECIFIC_HEAT_CAPACITY_DIMS);
    let density=PropertyKey::new("density",Density::DIMS);
    let mut claims=ClaimSet::new();
    let mut add=|key:PropertyKey,a:f64,b:f64,low:f64,high:f64| {
        let dims=key.quantity().dims();
        claims.insert_claim(PropertyClaim {key,
            value:PropertyValue::Curve {abscissa:"T".into(),abscissa_dims:Temperature::DIMS,
                knots:vec![(300.,a),(500.,b)],dims},
            validity:ValidityDomain::unconstrained().with_quantity("T",
                QuantitySpec::dimensional(Temperature::DIMS),low,high).with("lot",7.,7.),
            uncertainty:UncertaintyModel::Unstated,interpolation:InterpolationPolicy::LinearInside,
            observations:vec![],provenance:Provenance {source:"synthetic thermoelastic fixture".into(),
                license:"CC0-1.0".into(),artifact:None}}).unwrap();
    };
    add(heat.clone(),cp,cp,300.,500.);
    add(density.clone(),6.,2.,300.,500.);
    add(PropertyKey::new("young_modulus",Pressure::DIMS),young,0.6*young,300.,if split {349.} else {500.});
    if split { add(PropertyKey::new("young_modulus",Pressure::DIMS),2.*young,1.2*young,350.,500.); }
    add(PropertyKey::new("poisson_ratio",Dims([0;6])),nu,nu,300.,500.);
    add(PropertyKey::new("linear_thermal_expansion_coefficient",INVERSE_TEMPERATURE_DIMS),
        alpha,3.*alpha,300.,500.);
    let card=MaterialCard::assemble(MaterialStateId {chemistry:"authored specimen".into(),phase:"solid".into(),
        process:"synthetic".into(),revision:0},claims,vec![]).unwrap();
    let caloric=ResolvedHeatCapacityCurve::try_new(&card,HeatCapacitySourceQuery {
        heat_capacity:&heat,density:&density,temperature_axis:"T",lower:&point(300.),upper:&point(500.),
        selection:SelectionPolicy::SingleClaimOnly},SolidLiquidPhase::Solid,0.,1e-8).unwrap();
    (card,caloric)
}
fn config()->SolidEquilibriumConfig {
    SolidEquilibriumConfig {static_solve:TetStaticSolveConfig::default(),mesh_frame:ContentHash([1;32]),
        maximum_nodal_displacement_m:0.1}
}
fn thermal_config()->StepConfig { StepConfig {max_dt_s:0.125,energy_tolerance_j:1e-8} }
fn initial<'a>(curve:&'a EquilibriumEnthalpyPhaseCurve, cells:&[(u64,f64)],mass:f64)
    ->DepositionSimulation<'a> {
    let events=cells.iter().enumerate().map(|(i,&(id,power))| TimedEvent {time_s:0.,
        event:Event::Birth(DepositedCell {id:CellId(id),segment_id:id,body_id:77,material_id:MaterialId(0),
            mass_kg:mass,specific_h_j_kg:0.,geometry:CellGeometry {
                start_m:[0.,0.,i as f64],end_m:[1.,0.,i as f64],volume_m3:1./6.,surface_area_m2:3.},
            boundary:Boundary {ambient_k:300.,htc_w_m2_k:0.,reservoir:None,source_w:power}})}).collect();
    DepositionSimulation::new(events,&[(MaterialId(0),curve as &dyn ThermalLaw)],DepositionBudget {
        max_cells:8,max_links:8,max_events:16,max_substeps:1000}).unwrap()
        .advance_to(0.,thermal_config(),||Ok(())).unwrap().0
}
fn source<'a>(id:u64,card:&'a MaterialCard,curve:&'a EquilibriumEnthalpyPhaseCurve,context:&'a QueryPoint)
    ->TetCellMaterial<'a> {
    TetCellMaterial {cell:CellId(id),curve,point:ThermalMaterialPoint {
        material:card,context,temperature_axis:"T",liquid_fraction_axis:None},
        mechanical_selection:MaterialPropertySelection::SingleClaimOnly,
        stiffness_model:ThermalStiffnessModel::Isotropic,
        expansion_model:ThermalExpansionModel::Isotropic,
        expansion_selection:MaterialPropertySelection::SingleClaimOnly}
}
fn near(a:f64,b:f64,tol:f64) { assert!((a-b).abs()<=tol,"{a} != {b}, tolerance {tol}"); }

#[test]
fn caloric_and_expansion_data_drive_free_motion_and_constrained_stress_without_changing_mass() {
    with_cx(|cx| {
        for (young,alpha,cp) in [(1000.,1e-5,10.),(1000.,2e-5,10.),(2000.,1e-5,20.)] {
            let (card,caloric)=material(young,alpha,cp,0.25,false);
            let thermal=initial(caloric.curve(),&[(30,1000.)],1.);
            let context=point(300.);
            let model=SolidEnthalpy::try_new(&thermal,&NODES,&TETS,
                vec![source(30,&card,caloric.curve(),&context)],TetAssemblyBudget::standard(),0.1,config().mesh_frame,cx).unwrap();
            let loads=model.geometry().assemble(&[],VolumeLoad3::None,&[],cx).unwrap();
            let free=model.advance_to(&thermal,1.,thermal_config(),&loads,&FREE,config(),cx).unwrap();
            let rise=1000./cp;
            let strain=alpha*(rise+0.005*rise*rise); // exact integral of the authored alpha(T)
            let modulus=young*(1.-0.002*rise);
            near(free.thermal.temperature_k(CellId(30)).unwrap(),300.+rise,1e-10);
            near(free.thermal_audit.external_heat_j,1000.,1e-8);
            near(free.state.samples()[0].reference_elastic.property("density").unwrap().value_si(),6.,1e-10);
            assert!(free.state.samples()[0].elastic.property("density").is_none());
            near(free.state.materials()[0].density_kg_m3,6.,1e-12);
            near(free.state.problem(&FREE).assemble(cx).unwrap().total_mass_kg,1.,1e-12);
            for (node,p) in NODES.iter().enumerate() { for c in 0..3 {
                near(free.equilibrium.solution.displacement.displacement_m()[node][c],strain*p[c],1e-10);
                near(free.equilibrium.deformed_mesh.nodes_m()[node][c],(1.+strain)*p[c],1e-10);
            }}
            for s in free.equilibrium.solution.stress.elements()[0].stress_mandel_pa() { near(*s,0.,1e-7); }
            let fixed:Vec<_>=(0..12).collect();
            let clamped=model.advance_to(&thermal,1.,thermal_config(),&loads,&fixed,config(),cx).unwrap();
            let stress=clamped.equilibrium.solution.stress.elements()[0].stress_mandel_pa();
            for &s in &stress[..3] { near(s,-modulus*strain/(1.-2.*0.25),1e-8); }
            assert_eq!(clamped.equilibrium.deformed_mesh.nodes_m(),&NODES);
            assert_eq!(free.state.samples()[0].elastic.card_identity(),card.content_hash());
            assert_eq!(thermal.time_s(),0.);
            assert_eq!(thermal.cells()[&CellId(30)].specific_h_j_kg,0.);
        }
    });
}

#[test]
fn spatial_heating_and_element_order_produce_compatible_nonuniform_deformation() {
    with_cx(|cx| {
        let nodes=[[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.],[0.,0.,-1.]];
        let tets=[[0,1,2,3],[0,2,1,4]];
        let (card,caloric)=material(1000.,1e-5,10.,0.,false);
        let thermal=initial(caloric.curve(),&[(30,1000.),(10,500.)],1.);
        let context=point(300.);
        let model=SolidEnthalpy::try_new(&thermal,&nodes,&tets,vec![
            source(30,&card,caloric.curve(),&context),source(10,&card,caloric.curve(),&context)],
            TetAssemblyBudget::standard(),0.1,config().mesh_frame,cx).unwrap();
        let loads=model.geometry().assemble(&[],VolumeLoad3::None,&[],cx).unwrap();
        let next=model.advance_to(&thermal,1.,thermal_config(),&loads,&FREE,config(),cx).unwrap();
        let (hot,cold)=(1e-5*150.,1e-5*62.5);
        let (eh,ec)=(800.,900.);
        // nu=0: shared-face in-plane strain is stiffness-weighted; each free
        // apex expands by its OWN thermal strain normal to that face.
        let shared=(eh*hot+ec*cold)/(eh+ec);
        let u=next.equilibrium.solution.displacement.displacement_m();
        near(u[1][0],shared,1e-10); near(u[2][1],shared,1e-10);
        near(u[3][2],hot,1e-10); near(u[4][2],-cold,1e-10);
        for (i,(e,eps)) in [(eh,hot),(ec,cold)].into_iter().enumerate() {
            let stress=next.equilibrium.solution.stress.elements()[i].stress_mandel_pa();
            near(stress[0],e*(shared-eps),1e-8);
            near(stress[1],e*(shared-eps),1e-8); near(stress[2],0.,1e-8);
        }
        assert_eq!(next.state.samples().iter().map(|s|s.cell).collect::<Vec<_>>(),vec![CellId(30),CellId(10)]);
        near(next.thermal_audit.stored_energy_change_j,1500.,1e-8);
        near(next.state.problem(&FREE).assemble(cx).unwrap().total_mass_kg,2.,1e-12);
    });
}

#[test]
fn heated_loaded_equilibrium_reaches_the_native_prestress_pencil_with_exact_force_binding() {
    with_cx(|cx| {
        let (card,caloric)=material(1000.,1e-5,10.,0.,false);
        let thermal=initial(caloric.curve(),&[(30,1000.)],1.);
        let context=point(300.);
        let model=SolidEnthalpy::try_new(&thermal,&NODES,&TETS,
            vec![source(30,&card,caloric.curve(),&context)],TetAssemblyBudget::standard(),0.1,config().mesh_frame,cx).unwrap();
        let force=[(11,0.01)];
        let loads=model.geometry().assemble(&[],VolumeLoad3::None,&force,cx).unwrap();
        let fixed=[0,1,2,3,6,8,9];
        let next=model.advance_to(&thermal,1.,thermal_config(),&loads,&fixed,config(),cx).unwrap();
        let problem=next.state.problem(&fixed);
        let tangent=problem.assemble_thermal_tangent_with_nodal_forces(
            &next.equilibrium.solution.displacement,next.state.thermal_strains(),loads.nodal_forces_n(),
            config().mesh_frame,0.1,cx).unwrap();
        let i=tangent.pencil().free_dofs.iter().position(|&d|d==4).unwrap();
        near(tangent.geometric_stiffness().get(i,i),-800.*0.0015/6.,1e-8);
        near(tangent.stress().elements()[0].stress_mandel_pa()[2],6.*0.01,1e-8);
        near(tangent.pencil().total_mass_kg,1.,1e-12);
        assert_eq!(tangent.stress().identity(),next.equilibrium.solution.stress.identity());
        assert!(matches!(problem.assemble_thermal_tangent_with_nodal_forces(
            &next.equilibrium.solution.displacement,next.state.thermal_strains(),&[(11,0.02)],
            config().mesh_frame,0.1,cx),Err(TetElasticError::StressRecoveryInputMismatch)));
    });
}

#[test]
fn selected_elastic_source_cannot_switch_when_heating_crosses_its_support() {
    with_cx(|cx| {
        let (card,caloric)=material(1000.,1e-5,10.,0.25,true);
        let thermal=initial(caloric.curve(),&[(30,1000.)],1.);
        let context=point(300.);
        let model=SolidEnthalpy::try_new(&thermal,&NODES,&TETS,
            vec![source(30,&card,caloric.curve(),&context)],TetAssemblyBudget::standard(),0.1,config().mesh_frame,cx).unwrap();
        let loads=model.geometry().assemble(&[],VolumeLoad3::None,&[],cx).unwrap();
        let initial_claim=model.resolve(&thermal,cx).unwrap().samples()[0].elastic
            .property("young_modulus").unwrap().answer().receipt.selected;
        let independent=resolve_isotropic_elastic_state_point(&card,&point(400.),
            MaterialPropertySelection::SingleClaimOnly).unwrap();
        assert_ne!(initial_claim,independent.resolved().property("young_modulus").unwrap().answer().receipt.selected);
        assert!(matches!(model.advance_to(&thermal,1.,thermal_config(),&loads,&FREE,config(),cx),
            Err(SolidEnthalpyError::Material(_))));
        let supported=model.advance_to(&thermal,0.25,thermal_config(),&loads,&FREE,config(),cx).unwrap();
        assert_eq!(initial_claim,supported.state.samples()[0].elastic
            .property("young_modulus").unwrap().answer().receipt.selected);
        assert_eq!(thermal.time_s(),0.);
        assert_eq!(thermal.cells()[&CellId(30)].specific_h_j_kg,0.);
    });
}

#[test]
fn foreign_law_mass_melt_and_geometric_domain_refusals_leave_the_initial_state_unchanged() {
    with_cx(|cx| {
        let (card,caloric)=material(1000.,1e-5,10.,0.25,false);
        let curve=caloric.curve();
        let thermal=initial(curve,&[(30,1000.)],1.);
        let context=point(300.);
        let clone=curve.clone();
        assert!(SolidEnthalpy::try_new(&thermal,&NODES,&TETS,
            vec![source(30,&card,&clone,&context)],TetAssemblyBudget::standard(),0.1,config().mesh_frame,cx).is_err());
        let wrong_mass=initial(curve,&[(30,1000.)],2.);
        assert!(SolidEnthalpy::try_new(&wrong_mass,&NODES,&TETS,
            vec![source(30,&card,curve,&context)],TetAssemblyBudget::standard(),0.1,config().mesh_frame,cx).is_err());
        let model=SolidEnthalpy::try_new(&thermal,&NODES,&TETS,
            vec![source(30,&card,curve,&context)],TetAssemblyBudget::standard(),0.1,config().mesh_frame,cx).unwrap();
        let loads=model.geometry().assemble(&[],VolumeLoad3::None,&[],cx).unwrap();
        let small=SolidEquilibriumConfig {maximum_nodal_displacement_m:1e-5,..config()};
        assert!(model.advance_to(&thermal,1.,thermal_config(),&loads,&FREE,small,cx).is_err());
        let melt=EquilibriumEnthalpyPhaseCurve::try_new(card.content_hash(),vec![
            EnthalpyPhaseKnot {specific_enthalpy_j_kg:0.,temperature_k:300.,liquid_mass_fraction:0.,bulk_density_kg_m3:6.},
            EnthalpyPhaseKnot {specific_enthalpy_j_kg:1000.,temperature_k:400.,liquid_mass_fraction:0.,bulk_density_kg_m3:4.},
            EnthalpyPhaseKnot {specific_enthalpy_j_kg:1100.,temperature_k:400.,liquid_mass_fraction:1.,bulk_density_kg_m3:3.},
        ]).unwrap();
        let molten=initial(&melt,&[(30,1050.)],1.);
        let melt_model=SolidEnthalpy::try_new(&molten,&NODES,&TETS,
            vec![source(30,&card,&melt,&context)],TetAssemblyBudget::standard(),0.1,config().mesh_frame,cx).unwrap();
        let loads=melt_model.geometry().assemble(&[],VolumeLoad3::None,&[],cx).unwrap();
        assert!(matches!(melt_model.advance_to(&molten,1.,thermal_config(),&loads,&FREE,config(),cx),
            Err(SolidEnthalpyError::Invalid(_))));
        assert_eq!(thermal.time_s(),0.); assert_eq!(molten.time_s(),0.);
        assert_eq!(thermal.cells()[&CellId(30)].specific_h_j_kg,0.);
        assert_eq!(molten.cells()[&CellId(30)].specific_h_j_kg,0.);
    });
}

#[test]
fn intermediate_melting_cannot_be_hidden_by_a_fully_solid_final_state() {
    with_cx(|cx| {
        let (card,_)=material(1000.,1e-5,10.,0.25,false);
        let curve=EquilibriumEnthalpyPhaseCurve::try_new(card.content_hash(),vec![
            EnthalpyPhaseKnot {specific_enthalpy_j_kg:0.,temperature_k:300.,liquid_mass_fraction:0.,bulk_density_kg_m3:6.},
            EnthalpyPhaseKnot {specific_enthalpy_j_kg:1000.,temperature_k:400.,liquid_mass_fraction:0.,bulk_density_kg_m3:4.},
            EnthalpyPhaseKnot {specific_enthalpy_j_kg:1100.,temperature_k:400.,liquid_mass_fraction:1.,bulk_density_kg_m3:3.},
        ]).unwrap();
        let seed=initial(&curve,&[(30,1050.)],1.);
        let cell=seed.cells()[&CellId(30)].clone();
        let cooling=Boundary {source_w:-1050.,..cell.boundary};
        let thermal=DepositionSimulation::new(vec![
            TimedEvent {time_s:0.,event:Event::Birth(cell)},
            TimedEvent {time_s:1.,event:Event::Boundary(CellId(30),cooling)},
        ],&[(MaterialId(0),&curve as &dyn ThermalLaw)],DepositionBudget {
            max_cells:1,max_links:0,max_events:2,max_substeps:100}).unwrap()
            .advance_to(0.,thermal_config(),||Ok(())).unwrap().0;
        let hot=thermal.advance_to(1.,thermal_config(),||Ok(())).unwrap().0;
        assert!(!curve.state_at_specific_enthalpy(hot.cells()[&CellId(30)].specific_h_j_kg)
            .unwrap().is_fully_solid());
        let cold=thermal.advance_to(2.,thermal_config(),||Ok(())).unwrap().0;
        assert!(curve.state_at_specific_enthalpy(cold.cells()[&CellId(30)].specific_h_j_kg)
            .unwrap().is_fully_solid());
        let context=point(300.);
        let model=SolidEnthalpy::try_new(&thermal,&NODES,&TETS,
            vec![source(30,&card,&curve,&context)],TetAssemblyBudget::standard(),0.1,config().mesh_frame,cx).unwrap();
        let loads=model.geometry().assemble(&[],VolumeLoad3::None,&[],cx).unwrap();
        assert!(matches!(model.advance_to(&thermal,2.,thermal_config(),&loads,&FREE,config(),cx),
            Err(SolidEnthalpyError::Invalid(_))));
        assert_eq!(thermal.time_s(),0.);
        assert_eq!(thermal.cells()[&CellId(30)].specific_h_j_kg,0.);
    });
}
