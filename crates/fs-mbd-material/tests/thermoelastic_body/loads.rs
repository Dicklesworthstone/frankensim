//! Material density, dead load and heating share one accepted 3-D state.
use super::*;
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_matdb::SelectionPolicy;
use fs_mbd_material::thermal::heating::ThermalRegion;
use fs_mbd_material::thermal::heating::constrained::{ConstrainedHeatedBody, ConstrainedRegion};
use fs_modal::eigh_gen_dense;
use fs_thermal::enthalpy::{Boundary, CellGeometry, CellId, DepositionBudget, Event, MaterialId,
    StepConfig, TimedEvent};

#[test]
fn fixed_force_material_swap_changes_shape_and_retains_the_loaded_tangent() {
    for young in [500.,1000.] {
        let cards=[card("a",6.,[1000.,young],1e-4)];
        let original=reference(&cards);
        let mesh=TetrahedralMass::bind(&original,&NODES[..4],&TETS[..1],&KEYS[..1],4,||false).unwrap();
        // The second load is on a free DOF. The first is applied directly at a
        // support and must change its reaction without changing the strain.
        let forces=[(0,2.),(3,1.)];
        let body=with_cx(|cx|mesh.solve_thermal_with_nodal_forces(
            &inputs(&cards,&[400.]),&SUPPORT,&forces,config(),cx)).unwrap();
        let u=body.displacement().displacement_m();
        // Unit P1 tetrahedron: V=1/6 and dN1/dx=1, so sigma_xx=6 F.
        near(u[1][0],0.01+6./young);
        near(u[2][1],0.01-0.25*6./young);
        near(u[3][2],0.01-0.25*6./young);
        near(body.stress().elements()[0].stress_mandel_pa()[0],6.);
        near(body.stress().elements()[0].stress_mandel_pa()[1],0.);
        near(body.stress().elements()[0].stress_mandel_pa()[2],0.);
        near(body.support_reactions_n()[0][0],-3.);
        assert_eq!(body.nodal_forces_n(),&forces[..]);
        assert_eq!(body.body().inertia().mass(),original.inertia().mass());
        let pencil=with_cx(|cx|body.initial_stress_pencil(cx)).unwrap();
        assert_eq!(pencil.stress(),body.stress());
        assert_eq!(pencil.pencil().free_dofs,vec![3,6,7,9,10,11]);
        // Only N1 has an x gradient among the free nodes. V*sigma_xx=1.
        for (i,value) in pencil.geometric_stiffness().to_dense().iter().enumerate() {
            near(*value,if i==0 {1.}else{0.});
        }
    }
}

#[test]
fn density_drives_weight_and_thermal_deformation_keeps_that_reference_matter() {
    for rho in [6.,12.] {
        let cards=[card("a",rho,[1000.,1000.],1e-4)];
        let original=reference(&cards);
        let mesh=TetrahedralMass::bind(&original,&NODES[..4],&TETS[..1],&KEYS[..1],4,||false).unwrap();
        let forces=with_cx(|cx|mesh.gravity_nodal_forces([0.,0.,-4.],config().assembly,4,cx)).unwrap();
        for &(dof,force) in &forces {
            assert_eq!(dof%3,2);
            near(force,-rho/6.);
        }
        assert_eq!(forces.len(),4);
        let body=with_cx(|cx|mesh.solve_thermal_with_nodal_forces(
            &inputs(&cards,&[400.]),&SUPPORT,&forces,config(),cx)).unwrap();
        near(body.stress().elements()[0].stress_mandel_pa()[2],-rho);
        near(body.displacement().displacement_m()[3][2],0.01-rho/1000.);
        near(body.support_reactions_n().iter().map(|f|f[2]).sum(),4.*original.inertia().mass());
        assert_eq!(body.body().inertia().mass(),original.inertia().mass());
        assert!(body.sources()[0].stiffness().property(DENSITY_PROPERTY).is_none());
        assert_ne!(body.body().regions()[0].density_kg_m3(),rho);
        assert_eq!(forces,with_cx(|cx|mesh.gravity_nodal_forces([0.,0.,-4.],config().assembly,4,cx)).unwrap());
    }
}

#[test]
fn heterogeneous_gravity_uses_the_original_region_to_element_map() {
    let cards=[card("a",6.,[1000.,1000.],0.),card("b",12.,[2000.,2000.],0.)];
    let original=reference(&cards);
    let mesh=TetrahedralMass::bind(&original,&NODES,&[TETS[1],TETS[0]],&["b","a"],5,||false).unwrap();
    let forces=with_cx(|cx|mesh.gravity_nodal_forces([0.,0.,-4.],config().assembly,8,cx)).unwrap();
    assert_eq!(forces.len(),5);
    for ((dof,value),(expected_dof,expected_value)) in forces.iter()
        .zip([(2,-3.),(5,-3.),(8,-3.),(11,-1.),(14,-2.)]) {
        assert_eq!(*dof,expected_dof); near(*value,expected_value);
    }
}

fn heating_card(cp:f64)->MaterialCard {
    let mut claims=ClaimSet::new();
    for (key,dims,lo,hi) in [
        (DENSITY_PROPERTY,Density::DIMS,6.,6.),
        (SPECIFIC_HEAT_CAPACITY_PROPERTY,SPECIFIC_HEAT_CAPACITY_DIMS,cp,cp),
        (YOUNG_MODULUS_PROPERTY,Pressure::DIMS,1000.,500.),
        (POISSON_RATIO_PROPERTY,Dims([0;6]),0.25,0.25),
        (LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,INVERSE_TEMPERATURE_DIMS,1e-4,1e-4),
    ] {
        claims.insert_claim(PropertyClaim {key:PropertyKey::new(key,dims),
            value:PropertyValue::Curve {abscissa:"T".into(),abscissa_dims:Temperature::DIMS,
                knots:vec![(300.,lo),(400.,hi)],dims},
            validity:ValidityDomain::unconstrained().with("T",300.,400.).with("pressure",1e5,1e5),
            uncertainty:UncertaintyModel::Unstated,interpolation:InterpolationPolicy::LinearInside,
            observations:vec![],provenance:Provenance {source:"synthetic loaded thermal specimen; not measurements".into(),
                license:"CC0-1.0".into(),artifact:None}}).unwrap();
    }
    MaterialCard::assemble(MaterialStateId {chemistry:"a".into(),phase:"solid".into(),
        process:"synthetic fixture".into(),revision:0},claims,vec![]).unwrap()
}
fn caloric(card:&MaterialCard)->ResolvedHeatCapacityCurve {
    let point=|t|QueryPoint::new().with("T",t).unwrap().with("pressure",1e5).unwrap();
    ResolvedHeatCapacityCurve::try_new(card,HeatCapacitySourceQuery {
        heat_capacity:&PropertyKey::new(SPECIFIC_HEAT_CAPACITY_PROPERTY,SPECIFIC_HEAT_CAPACITY_DIMS),
        density:&PropertyKey::new(DENSITY_PROPERTY,Density::DIMS),temperature_axis:"T",
        lower:&point(300.),upper:&point(400.),selection:SelectionPolicy::SingleClaimOnly,
    },SolidLiquidPhase::Solid,0.,1e-8).unwrap()
}
fn boundary(w:f64)->Boundary {Boundary {ambient_k:300.,htc_w_m2_k:0.,reservoir:None,source_w:w}}
fn step()->StepConfig {StepConfig {max_dt_s:0.25,energy_tolerance_j:1e-8}}
fn heated<'a>(reference:&'a ResolvedBody,card:&'a MaterialCard,curve:&'a ResolvedHeatCapacityCurve,
    cx:&Cx<'_>)->ConstrainedHeatedBody<'a> {
    let mesh=TetrahedralMass::bind(reference,&NODES[..4],&TETS[..1],&KEYS[..1],4,||false).unwrap();
    let storage=ThermalRegion {key:"a",card,caloric:curve,cell_id:CellId(0),material_id:MaterialId(0),
        segment_id:0,body_id:0,geometry:CellGeometry {start_m:[0.,0.,0.],end_m:[0.5,0.5,0.5],
            volume_m3:reference.regions()[0].geometry().volume_m3(),surface_area_m2:2.4},
        initial_specific_enthalpy_j_kg:0.,boundary:boundary(100.),
        expansion_selection:MaterialPropertySelection::SingleClaimOnly,fixed_point_m:Vec3::default()};
    ConstrainedHeatedBody::new_with_nodal_forces(mesh,
        vec![ConstrainedRegion {storage,stiffness_selection:MaterialPropertySelection::SingleClaimOnly,
            stiffness_model:ThermalStiffnessModel::Isotropic,
            expansion_model:fs_mbd_material::deformation::thermoelastic::ThermalExpansionModel::Isotropic}],
        SUPPORT.to_vec(),vec![(0,2.),(3,1.)],config(),
        vec![TimedEvent {time_s:1.,event:Event::Boundary(CellId(0),boundary(-100.))}],
        DepositionBudget {max_cells:1,max_links:0,max_events:2,max_substeps:100},step(),cx).unwrap()
}
fn spectrum(body:&ConstrainedHeatedBody<'_>)->Vec<f64> {
    let pencil=with_cx(|cx|body.configuration().initial_stress_pencil(cx)).unwrap();
    let a=pencil.pencil();
    eigh_gen_dense(&a.stiffness.to_dense(),&a.mass.to_dense(),a.free_dofs.len()).unwrap()
        .into_iter().map(|mode|mode.lambda).collect()
}

#[test]
fn heat_capacity_changes_loaded_shape_and_modes_on_the_existing_heat_clock() {
    let mut displacement=Vec::new(); let mut spectra=Vec::new();
    for cp in [10.,20.] {
        let cards=[heating_card(cp)]; let original=reference(&cards); let curve=caloric(&cards[0]);
        let initial=with_cx(|cx|heated(&original,&cards[0],&curve,cx));
        let (hot,report)=with_cx(|cx|initial.advance_to(1.,step(),4,cx)).unwrap();
        assert_eq!(report.shape_solves,4);
        near(report.thermal.external_heat_j,100.);
        let delta_t=100./cp; let young=1000.-5.*delta_t;
        near(hot.thermal().temperature_k(CellId(0)).unwrap(),300.+delta_t);
        let accepted=hot.configuration();
        near(accepted.displacement().displacement_m()[1][0],1e-4*delta_t+6./young);
        near(accepted.stress().elements()[0].stress_mandel_pa()[0],6.);
        near(accepted.support_reactions_n()[0][0],-3.);
        assert_eq!(accepted.nodal_forces_n(),initial.nodal_forces_n());
        assert_eq!(accepted.body().inertia().mass(),original.inertia().mass());
        let direct=with_cx(|cx|hot.mesh().solve_thermal_with_nodal_forces(
            &inputs(&cards,&[300.+delta_t]),&SUPPORT,hot.nodal_forces_n(),config(),cx)).unwrap();
        assert_eq!(accepted,&direct);
        displacement.push(accepted.displacement().displacement_m()[1][0]); spectra.push(spectrum(&hot));

        let (half,_)=with_cx(|cx|initial.advance_to(0.5,step(),2,cx)).unwrap();
        let (resumed,_)=with_cx(|cx|half.advance_to(1.,step(),2,cx)).unwrap();
        assert_eq!(hot.configuration(),resumed.configuration());
        let (cooled,audit)=with_cx(|cx|hot.advance_to(2.,step(),4,cx)).unwrap();
        near(audit.thermal.external_heat_j,-100.);
        assert_eq!(cooled.configuration(),initial.configuration());
        assert!(with_cx(|cx|initial.advance_to(1.,step(),3,cx)).is_err());
        near(initial.thermal().time_s(),0.);
        assert_eq!(initial.nodal_forces_n(),&[(0,2.),(3,1.)][..]);
    }
    assert!(displacement[0]>displacement[1]);
    for (hot,cooler) in spectra[0].iter().zip(&spectra[1]) {assert!(hot<cooler);}
}

#[test]
fn unloaded_identity_and_refusal_preserve_the_accepted_snapshot() {
    let cards=[card("a",6.,[1000.,700.],1e-4)]; let original=reference(&cards);
    let mesh=TetrahedralMass::bind(&original,&NODES[..4],&TETS[..1],&KEYS[..1],4,||false).unwrap();
    let states=inputs(&cards,&[400.]);
    let plain=with_cx(|cx|mesh.solve_thermal(&states,&SUPPORT,config(),cx)).unwrap();
    let empty=with_cx(|cx|mesh.solve_thermal_with_nodal_forces(&states,&SUPPORT,&[],config(),cx)).unwrap();
    assert_eq!(plain,empty);
    for forces in [vec![(3,1.),(3,2.)],vec![(12,1.)],vec![(3,f64::NAN)]] {
        assert!(matches!(with_cx(|cx|mesh.solve_thermal_with_nodal_forces(
            &states,&SUPPORT,&forces,config(),cx)),Err(ThermalBodyError::Solid(TetElasticError::InvalidNodalForce {..}))));
    }
    let before=plain.clone(); let gate=CancelGate::new(); gate.request();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx=Cx::new(&gate,arena,StreamKey {seed:3,kernel_id:17,tile:0,iteration:0},Budget::INFINITE,ExecMode::Deterministic);
        assert!(mesh.gravity_nodal_forces([0.,0.,-4.],config().assembly,4,&cx).is_err());
        assert!(mesh.solve_thermal_with_nodal_forces(&states,&SUPPORT,&[(3,1.)],config(),&cx).is_err());
    });
    assert_eq!(plain,before);
}
