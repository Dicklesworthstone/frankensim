//! Source data and solved heat drive the existing signed vibration pencil.
use super::*;
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_matdb::SelectionPolicy;
use fs_modal::eigh_gen_dense;
use fs_solid::linear3::{TetElasticMaterial, TetLinearElasticProblem, TetMaterialField,
    TetThermalStrainField, TetThermalStrainState};
use fs_solid::prestress::PrestressedTetAssembly;
use fs_mbd_material::deformation::thermoelastic::ThermalBody;
use fs_mbd_material::thermal::heating::ThermalRegion;
use fs_mbd_material::thermal::heating::constrained::{ConstrainedHeatedBody, ConstrainedRegion};
use fs_thermal::enthalpy::{Boundary, CellGeometry, CellId, DepositionBudget,
    Event, MaterialId, StepConfig, TimedEvent};

const AXIAL: [usize; 7] = [0,1,2,3,6,8,9];
fn spectrum(p: &PrestressedTetAssembly) -> Vec<f64> {
    let a=p.pencil();
    eigh_gen_dense(&a.stiffness.to_dense(),&a.mass.to_dense(),a.free_dofs.len()).unwrap()
        .into_iter().map(|mode|mode.lambda).collect()
}

// Independently construct the old solid API from the retained source receipts.
// Production callers never need to repeat this binding or copy coefficients.
fn direct(body:&ThermalBody<'_>, mesh:&TetrahedralMass<'_>, fixed:&[usize], cx:&Cx<'_>) -> PrestressedTetAssembly {
    let mut ordered:Vec<_>=body.sources().iter().zip(mesh.reference().regions())
        .zip(mesh.region_elements()).map(|((source,region),&e)| {
            let property=|key|source.stiffness().property(key).unwrap().value_si();
            let material=TetElasticMaterial::try_new(
                region.material().property(DENSITY_PROPERTY).unwrap().value_si(),
                property(YOUNG_MODULUS_PROPERTY),property(POISSON_RATIO_PROPERTY),
                *source.material_state_identity().as_bytes()).unwrap();
            let a=source.expansion().free_linear_strain();
            let point=|source:&ResolvedMaterialStatePoint|source.query_point().iter()
                .find(|(axis,_)|axis=="T").unwrap().1;
            let thermal=TetThermalStrainState::try_new(point(source.stiffness()),
                point(source.expansion().reference().resolved()),[a,a,a,0.,0.,0.],
                *source.material_state_identity().as_bytes(),source.expansion().identity(),
                config().maximum_absolute_mandel_strain).unwrap();
            (e,material,thermal)
        }).collect();
    ordered.sort_by_key(|row|row.0);
    let materials:Vec<_>=ordered.iter().map(|r|r.1.clone()).collect();
    let strains:Vec<_>=ordered.iter().map(|r|r.2).collect();
    let p=TetLinearElasticProblem {nodes_m:mesh.nodes_m(),tetrahedra:mesh.tetrahedra(),
        materials:TetMaterialField::PerElement(&materials),fixed_dofs:fixed,budget:config().assembly};
    p.assemble_thermal_tangent(body.displacement(),TetThermalStrainField::PerElement(&strains),
        config().mesh_frame,config().maximum_absolute_mandel_strain,cx).unwrap()
}

#[test]
fn data_only_expansion_swap_changes_the_actual_spectrum_not_reference_mass() {
    let mut spectra=Vec::new(); let mut masses=Vec::new();
    for alpha in [0.0,1e-4,2e-4] {
        let cards=[card("specimen",6.,[1000.,1000.],alpha)]; let reference=reference(&cards);
        let mesh=TetrahedralMass::bind(&reference,&NODES[..4],&TETS[..1],&KEYS[..1],4,||false).unwrap();
        let body=with_cx(|cx|mesh.solve_thermal(&inputs(&cards,&[400.]),&AXIAL,config(),cx)).unwrap();
        let pencil=with_cx(|cx|body.initial_stress_pencil(cx)).unwrap();
        let expected=with_cx(|cx|direct(&body,&mesh,&AXIAL,cx));
        assert_eq!(pencil.pencil().stiffness.to_dense(),expected.pencil().stiffness.to_dense());
        assert_eq!(pencil.identity(),expected.identity());
        assert_eq!(pencil.stress(),body.stress());
        assert!(body.sources()[0].stiffness().property(DENSITY_PROPERTY).is_none());
        near(pencil.stress().elements()[0].stress_mandel_pa()[0],-1000.*alpha*100.);
        masses.push(pencil.pencil().mass.to_dense()); spectra.push(spectrum(&pencil));
    }
    assert_eq!(masses[0],masses[1]); assert_eq!(masses[1],masses[2]);
    for i in 0..spectra[0].len() {assert!(spectra[2][i]<spectra[1][i] && spectra[1][i]<spectra[0][i]);}
}

#[test]
fn heterogeneous_source_and_element_order_are_retained_without_new_queries() {
    let cards=[card("a",6.,[1000.,600.],1e-4),card("b",12.,[2000.,1600.],2e-4)];
    let reference=reference(&cards);
    let mesh=TetrahedralMass::bind(&reference,&NODES,&[TETS[1],TETS[0]],&["b","a"],5,||false).unwrap();
    let mut input=inputs(&cards,&[375.,350.]);
    let body=with_cx(|cx|mesh.solve_thermal(&input,&SUPPORT,config(),cx)).unwrap();
    let a=with_cx(|cx|body.initial_stress_pencil(cx)).unwrap();
    let b=with_cx(|cx|direct(&body,&mesh,&SUPPORT,cx));
    assert_eq!(a.identity(),b.identity()); assert_eq!(a.pencil().stiffness.to_dense(),b.pencil().stiffness.to_dense());
    input.reverse();
    let reordered=with_cx(|cx|mesh.solve_thermal(&input,&SUPPORT,config(),cx)).unwrap();
    assert_eq!(a.identity(),with_cx(|cx|reordered.initial_stress_pencil(cx)).unwrap().identity());
    near(body.sources()[0].stiffness().property(YOUNG_MODULUS_PROPERTY).unwrap().value_si(),700.);
    near(body.sources()[1].stiffness().property(YOUNG_MODULUS_PROPERTY).unwrap().value_si(),1800.);
}

#[test]
fn snapshot_supports_and_cancellation_cannot_silently_release_the_solid() {
    let cards=[card("a",6.,[1000.,1000.],1e-4)]; let reference=reference(&cards);
    let mesh=TetrahedralMass::bind(&reference,&NODES[..4],&TETS[..1],&KEYS[..1],4,||false).unwrap();
    let all:Vec<_>=(0..12).collect();
    let fixed=with_cx(|cx|mesh.solve_thermal(&inputs(&cards,&[400.]),&all,config(),cx)).unwrap();
    assert!(matches!(with_cx(|cx|fixed.initial_stress_pencil(cx)),Err(ThermalBodyError::Solid(TetElasticError::NoFreeDofs))));
    let body=with_cx(|cx|mesh.solve_thermal(&inputs(&cards,&[400.]),&AXIAL,config(),cx)).unwrap();
    let before=body.clone(); let gate=CancelGate::new(); gate.request();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx=Cx::new(&gate,arena,StreamKey {seed:3,kernel_id:17,tile:0,iteration:0},Budget::INFINITE,ExecMode::Deterministic);
        assert!(matches!(body.initial_stress_pencil(&cx),Err(ThermalBodyError::Solid(TetElasticError::Cancelled))));
    });
    assert_eq!(body,before);
}

// A separate heating fixture gives density and Cp full caloric-chart support.
// All comparisons keep geometry, drive, stiffness and expansion identical.
fn heat_card(cp:f64) -> MaterialCard {
    let mut claims=ClaimSet::new();
    for (key,dims,value) in [(DENSITY_PROPERTY,Density::DIMS,6.),
        (SPECIFIC_HEAT_CAPACITY_PROPERTY,SPECIFIC_HEAT_CAPACITY_DIMS,cp),
        (YOUNG_MODULUS_PROPERTY,Pressure::DIMS,1000.),(POISSON_RATIO_PROPERTY,Dims([0;6]),0.25),
        (LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,INVERSE_TEMPERATURE_DIMS,0.001)] {
        claims.insert_claim(PropertyClaim {key:PropertyKey::new(key,dims),value:PropertyValue::Curve {
            abscissa:"T".into(),abscissa_dims:Temperature::DIMS,knots:vec![(300.,value),(400.,value)],dims},
            validity:ValidityDomain::unconstrained().with("T",300.,400.).with("pressure",1e5,1e5),
            uncertainty:UncertaintyModel::Unstated,interpolation:InterpolationPolicy::LinearInside,
            observations:vec![],provenance:Provenance {source:"synthetic thermal-vibration data; not measurements".into(),
                license:"CC0-1.0".into(),artifact:None}}).unwrap();
    }
    MaterialCard::assemble(MaterialStateId {chemistry:"specimen".into(),phase:"solid".into(),
        process:"synthetic fixture".into(),revision:0},claims,vec![]).unwrap()
}
fn caloric(card:&MaterialCard) -> ResolvedHeatCapacityCurve {
    let point=|t|QueryPoint::new().with("T",t).unwrap().with("pressure",1e5).unwrap();
    ResolvedHeatCapacityCurve::try_new(card,HeatCapacitySourceQuery {
        heat_capacity:&PropertyKey::new(SPECIFIC_HEAT_CAPACITY_PROPERTY,SPECIFIC_HEAT_CAPACITY_DIMS),
        density:&PropertyKey::new(DENSITY_PROPERTY,Density::DIMS),temperature_axis:"T",
        lower:&point(300.),upper:&point(400.),selection:SelectionPolicy::SingleClaimOnly,
    },SolidLiquidPhase::Solid,0.,1e-8).unwrap()
}
fn boundary(source_w:f64)->Boundary {Boundary {ambient_k:300.,htc_w_m2_k:0.,reservoir:None,source_w}}
fn step()->StepConfig {StepConfig {max_dt_s:0.25,energy_tolerance_j:1e-8}}
fn heated<'a>(reference:&'a ResolvedBody,card:&'a MaterialCard,caloric:&'a ResolvedHeatCapacityCurve,
    cx:&Cx<'_>)->ConstrainedHeatedBody<'a> {
    let mesh=TetrahedralMass::bind(reference,&NODES[..4],&TETS[..1],&KEYS[..1],4,||false).unwrap();
    let storage=ThermalRegion {key:"a",card,caloric,cell_id:CellId(0),material_id:MaterialId(0),
        segment_id:0,body_id:0,geometry:CellGeometry {start_m:[0.,0.,0.],end_m:[0.5,0.5,0.5],
            volume_m3:reference.regions()[0].geometry().volume_m3(),surface_area_m2:2.4},
        initial_specific_enthalpy_j_kg:0.,boundary:boundary(100.),
        expansion_selection:MaterialPropertySelection::SingleClaimOnly,fixed_point_m:Vec3::default()};
    ConstrainedHeatedBody::new(mesh,vec![ConstrainedRegion {storage,stiffness_selection:MaterialPropertySelection::SingleClaimOnly}],
        AXIAL.to_vec(),config(),vec![TimedEvent {time_s:1.,event:Event::Boundary(CellId(0),boundary(-100.))}],
        DepositionBudget {max_cells:1,max_links:0,max_events:2,max_substeps:100},step(),cx).unwrap()
}

#[test]
fn heat_capacity_alone_changes_solved_temperature_stress_and_eigenvalues() {
    let mut spectra=Vec::new(); let mut matrices=Vec::new();
    for cp in [10.,20.] {
        let cards=[heat_card(cp)]; let reference=reference(&cards); let curve=caloric(&cards[0]);
        let old=with_cx(|cx|heated(&reference,&cards[0],&curve,cx));
        let initial=with_cx(|cx|old.configuration().initial_stress_pencil(cx)).unwrap();
        let (hot,report)=with_cx(|cx|old.advance_to(1.,step(),4,cx)).unwrap();
        assert_eq!(report.shape_solves,4); near(report.thermal.external_heat_j,100.);
        near(hot.thermal().temperature_k(CellId(0)).unwrap(),300.+100./cp);
        let p=with_cx(|cx|hot.configuration().initial_stress_pencil(cx)).unwrap();
        near(p.stress().elements()[0].stress_mandel_pa()[0],-100./cp);
        assert_eq!(p.pencil().mass.to_dense(),initial.pencil().mass.to_dense());
        assert_eq!(p.stress().solution_identity(),hot.configuration().displacement().identity());
        assert_eq!(old.thermal().time_s(),0.);
        assert_eq!(initial.identity(),with_cx(|cx|old.configuration().initial_stress_pencil(cx)).unwrap().identity());
        matrices.push(p.pencil().mass.to_dense()); spectra.push(spectrum(&p));
    }
    assert_eq!(matrices[0],matrices[1]);
    for (hot,cooler) in spectra[0].iter().zip(&spectra[1]) {assert!(hot<cooler);}
}

#[test]
fn cooling_and_aligned_restart_refresh_the_pencil_without_compounding_stress() {
    let cards=[heat_card(10.)]; let reference=reference(&cards); let curve=caloric(&cards[0]);
    let old=with_cx(|cx|heated(&reference,&cards[0],&curve,cx));
    let (whole,_)=with_cx(|cx|old.advance_to(1.,step(),4,cx)).unwrap();
    let (half,_)=with_cx(|cx|old.advance_to(0.5,step(),2,cx)).unwrap();
    let (resumed,_)=with_cx(|cx|half.advance_to(1.,step(),2,cx)).unwrap();
    let a=with_cx(|cx|whole.configuration().initial_stress_pencil(cx)).unwrap();
    let b=with_cx(|cx|resumed.configuration().initial_stress_pencil(cx)).unwrap();
    assert_eq!(a.identity(),b.identity()); assert_eq!(a.pencil().stiffness.to_dense(),b.pencil().stiffness.to_dense());
    let (cool,report)=with_cx(|cx|whole.advance_to(2.,step(),4,cx)).unwrap();
    near(report.thermal.external_heat_j,-100.); near(cool.thermal().temperature_k(CellId(0)).unwrap(),300.);
    let reset=with_cx(|cx|cool.configuration().initial_stress_pencil(cx)).unwrap();
    let initial=with_cx(|cx|old.configuration().initial_stress_pencil(cx)).unwrap();
    assert_eq!(reset.pencil().stiffness.to_dense(),initial.pencil().stiffness.to_dense());
    assert_eq!(reset.pencil().mass.to_dense(),initial.pencil().mass.to_dense());
    assert_eq!(a.identity(),with_cx(|cx|whole.configuration().initial_stress_pencil(cx)).unwrap().identity());
}
