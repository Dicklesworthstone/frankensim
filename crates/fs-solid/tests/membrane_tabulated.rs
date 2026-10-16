//! Actual source -> nonlinear equilibrium -> prestressed vibration chain.
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet,InterpolationPolicy,MaterialCard,MaterialStateId,PropertyClaim,
    PropertyKey,PropertyValue,Provenance,QueryPoint};
use fs_qty::Dims;
use fs_material::sheet::{SheetCapabilities,SheetLaw,resolve_sheet,ResolvedSheet,MEMBRANE_DIMS};
use fs_material::sheet::tabulated::{CURVE_PROPERTIES,COUPLING_PROPERTY,STRAIN_AXIS};
use fs_material::state_point::MaterialPropertySelection;
use fs_solid::membrane::{MembraneBudget,MembraneDomain,MembranePatch,MembraneSolveSettings};
use fs_solid::membrane::vibration::VibrationSettings;
use fs_solid::mesh2::Mesh2;
use fs_solver::nonlinear::NewtonKrylovConfig;

fn card(hardening:f64, mass:f64, policy:InterpolationPolicy, restrict:bool)->MaterialCard {
    let mut claims=ClaimSet::new();
    let mut add=|key:&str,value:PropertyValue, interpolation, validity| {
        let dims=value.dims();
        claims.insert_claim(PropertyClaim {key:PropertyKey::new(key,dims),value,interpolation,validity,
            uncertainty:fs_matdb::UncertaintyModel::Unstated,observations:vec![],
            provenance:Provenance {source:"synthetic signed Green-strain law; not a measured fabric".into(),license:"test".into(),artifact:None}}).unwrap();
    };
    let valid=ValidityDomain::unconstrained().with("T",290.0,310.0);
    for (key,value,dims) in [("areal_density",mass,Dims([-2,1,0,0,0,0])),
        ("thickness",0.001,Dims([1,0,0,0,0,0])),(COUPLING_PROPERTY,0.0,MEMBRANE_DIMS)] {
        add(key,PropertyValue::Scalar {value,dims},InterpolationPolicy::ConstantWithinValidity,valid.clone());
    }
    for (i,key) in CURVE_PROPERTIES.iter().enumerate() {
        let knots=if i==2 {vec![(-0.2,-800.0),(0.0,0.0),(0.2,800.0)]}
            else {vec![(-0.2,-1600.0),(0.0,0.0),(0.025,200.0),(0.1,200.0+0.075*hardening)]};
        let validity=if restrict {valid.clone().with(STRAIN_AXIS,-0.02,0.035)} else {valid.clone()};
        add(key,PropertyValue::Curve {abscissa:STRAIN_AXIS.into(),abscissa_dims:Dims::NONE,knots,dims:MEMBRANE_DIMS},policy,validity);
    }
    MaterialCard::assemble(MaterialStateId {chemistry:"synthetic characterized sheet".into(),phase:"solid".into(),
        process:"fixed specimen".into(),revision:0},claims,vec![]).unwrap()
}
fn sheet(card:&MaterialCard)->ResolvedSheet {
    resolve_sheet(card,&QueryPoint::new().with("T",300.0).unwrap(),
        SheetCapabilities {membrane:true,bending:false},MaterialPropertySelection::SingleClaimOnly).unwrap()
}
fn patch(card:&MaterialCard)->MembranePatch {
    let mesh=Mesh2 {nodes:vec![[0.0,0.0],[0.1,0.0],[0.0,0.1],[-0.1,0.0],[0.0,-0.1]],
        elems:vec![vec![0,1,2],vec![0,2,3],vec![0,3,4],vec![0,4,1]],patches:vec![]};
    MembranePatch::try_new(mesh,vec![sheet(card);4],vec![0.0;4],MembraneDomain {
        minimum_area_ratio:0.01,maximum_abs_strain:0.25},MembraneBudget {max_nodes:5,max_triangles:4},||false).unwrap()
}
fn settings()->MembraneSolveSettings {
    MembraneSolveSettings {newton:NewtonKrylovConfig {absolute_tolerance:1e-9,relative_tolerance:1e-11,
        linear_restart:8,..Default::default()},max_iterations:80,max_evaluations:2000,max_krylov_restart:8}
}
fn fixed(p:&MembranePatch)->Vec<(usize,f64)> {
    let mut bc=vec![(0,0.0),(1,0.0)];
    for (n,x) in p.mesh().nodes.iter().enumerate().skip(1) {
        bc.extend([(3*n,0.02*x[0]),(3*n+1,0.02*x[1]),(3*n+2,0.0)]);
    }
    bc
}
#[test]
fn source_curves_and_receipts_survive_without_equivalent_scalar_moduli() {
    let c=card(20000.0,0.2,InterpolationPolicy::LinearInside,false);let s=sheet(&c);
    assert_eq!(s.material().card_identity(),c.content_hash());
    assert!(s.material().property("warp_stiffness").is_none());
    assert!(s.membrane().is_err());assert!(s.bending().is_err());
    let law=s.membrane_law().unwrap();
    assert!((law.response([0.04,0.0,0.0]).unwrap().tangent_n_m[0]-20000.0).abs()<1e-9);
    let other=card(8000.0,0.2,InterpolationPolicy::LinearInside,false);
    assert!(ResolvedSheet::try_from_material_with_card(s.material().clone(),&other).is_err());
    let detached=ResolvedSheet::try_from_material(s.material().clone()).unwrap();
    assert!(detached.membrane_law().is_err()); // no original curve card, no fabricated law
}
#[test]
fn hardening_changes_actual_equilibrium_and_the_exact_prestressed_tangent() {
    let p=patch(&card(20000.0,0.2,InterpolationPolicy::LinearInside,false));
    let mut loads=[[0.0;3];5];loads[0][2]=-53.12;
    let v=p.solve_vibration(&loads,&fixed(&p),&[[0.0;3];5],settings(),
        VibrationSettings {max_free_dofs:1,relative_residual_tolerance:1e-8},||false).unwrap();
    // At h=.02: E=.0402, S=200+20000*(E-.025)=504, gamma=.04.
    // F=4h(S+G*h²/L²)=53.12; dF/dh=4S+4C*h²/L²+12G*h²/L²=7136.
    assert!((v.equilibrium().displacement_m[0][2]+0.02).abs()<1e-9);
    let expected=7136.0/p.nodal_mass_kg()[0];
    assert!((v.modes()[0].lambda/expected-1.0).abs()<1e-8);
    let linear=patch(&card(8000.0,0.2,InterpolationPolicy::LinearInside,false));
    let other=linear.solve(&loads,&fixed(&linear),&[[0.0;3];5],settings(),||false).unwrap();
    assert!(other.displacement_m[0][2]< -0.02);
    assert_eq!(p.reference_mass_kg(),linear.reference_mass_kg());
}
#[test]
fn actual_source_support_and_interpolation_authority_control_the_solver() {
    let limited=card(20000.0,0.2,InterpolationPolicy::LinearInside,true);
    assert!(sheet(&limited).membrane_law().unwrap().response([0.04,0.0,0.0]).is_err());
    let exact=card(20000.0,0.2,InterpolationPolicy::TabulatedOnly,false);
    assert!(resolve_sheet(&exact,&QueryPoint::new().with("T",300.0).unwrap(),
        SheetCapabilities {membrane:true,bending:false},MaterialPropertySelection::SingleClaimOnly).is_err());
    let p=patch(&limited);let mut loads=[[0.0;3];5];loads[0][2]=-53.12;
    assert!(p.solve(&loads,&fixed(&p),&[[0.0;3];5],settings(),||false).is_err());
    assert!(p.solve(&loads,&fixed(&p),&[[0.0;3];5],settings(),||true).is_err());
}

#[test]
fn a_converged_kink_equilibrium_cannot_mint_a_unique_vibration_basis() {
    let p=patch(&card(20000.0,0.2,InterpolationPolicy::LinearInside,false));
    let stretch=1.05_f64.sqrt();let mut bc=vec![(0,0.0),(1,0.0)];
    for (n,x) in p.mesh().nodes.iter().enumerate().skip(1) {
        bc.extend([(3*n,(stretch-1.0)*x[0]),(3*n+1,(stretch-1.0)*x[1]),(3*n+2,0.0)]);
    }
    let zero=[[0.0;3];5];
    assert!(p.solve(&zero,&bc,&zero,settings(),||false).is_ok());
    assert!(p.solve_vibration(&zero,&bc,&zero,settings(),
        VibrationSettings {max_free_dofs:1,relative_residual_tolerance:1e-8},||false).is_err());
}
