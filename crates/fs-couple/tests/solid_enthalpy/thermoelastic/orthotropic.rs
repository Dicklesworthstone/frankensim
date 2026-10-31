//! The original anisotropic source, finite cells, actual FEM and inertia share
//! one thermal potential. Authored constants verify mechanics, not calibration.
use super::*;
use fs_couple::solid_enthalpy::{TetCellMaterial,ThermalExpansionModel,ThermalStiffnessModel};
use fs_material::state_point::{ORTHOTROPIC_YOUNG_MODULUS_PROPERTIES,
    ORTHOTROPIC_POISSON_RATIO_PROPERTIES,ORTHOTROPIC_SHEAR_MODULUS_PROPERTIES,
    ORTHOTROPIC_THERMAL_EXPANSION_COEFFICIENT_PROPERTIES};
use fs_solid::loads3::transient::TransientConfig3;

const IDENTITY:[[f64;3];3]=[[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]];
const TURN:[[f64;3];3]=[[0.8,-0.6,0.0],[0.6,0.8,0.0],[0.0,0.0,1.0]];
const E:[f64;3]=[1000.0,1600.0,2500.0];
const ET:[f64;3]=[2.0,-1.0,3.0];
const ALPHA:[f64;3]=[0.0003,-0.0001,0.0005];
const ALPHAT:[f64;3]=[0.000001,0.0000002,-0.0000005];

fn fixture(scale:f64,varying:bool,axial_only:bool)->(MaterialCard,EquilibriumEnthalpyPhaseCurve) {
    let mut claims=ClaimSet::new();
    let mut add=|name:&str,dims,value:f64,slope:f64| {
        claims.insert_claim(PropertyClaim {key:PropertyKey::new(name,dims),
            value:PropertyValue::Curve {abscissa:"T".into(),abscissa_dims:Temperature::DIMS,
                knots:vec![(280.0,value-20.0*slope),(360.0,value+60.0*slope)],dims},
            validity:ValidityDomain::unconstrained().with_quantity("T",quantity(),280.0,360.0),
            uncertainty:UncertaintyModel::Unstated,interpolation:InterpolationPolicy::LinearInside,
            observations:vec![],provenance:Provenance {source:"authored reciprocal orthotropic fixture".into(),
                license:"test-only".into(),artifact:None}}).unwrap();
    };
    add("density",Density::DIMS,6.0,0.0);
    add("specific_heat_capacity",SPECIFIC_HEAT_CAPACITY_DIMS,10.0,0.0);
    for i in 0..3 {
        add(ORTHOTROPIC_YOUNG_MODULUS_PROPERTIES[i],Pressure::DIMS,E[i],if varying {ET[i]} else {0.0});
        add(ORTHOTROPIC_POISSON_RATIO_PROPERTIES[i],Dims::NONE,0.0,0.0);
        add(ORTHOTROPIC_SHEAR_MODULUS_PROPERTIES[i],Pressure::DIMS,[300.0,600.0,700.0][i],0.0);
        let alpha=if axial_only {if i==0 {0.001} else {0.0}} else {ALPHA[i]};
        add(ORTHOTROPIC_THERMAL_EXPANSION_COEFFICIENT_PROPERTIES[i],INVERSE_TEMPERATURE_DIMS,
            scale*alpha,if varying {scale*ALPHAT[i]} else {0.0});
    }
    let card=MaterialCard::assemble(MaterialStateId {chemistry:"synthetic orthotropic source".into(),
        phase:"solid".into(),process:"numerical reciprocal fixture".into(),revision:0},claims,vec![]).unwrap();
    let curve=EquilibriumEnthalpyPhaseCurve::try_single_phase(card.content_hash(),SolidLiquidPhase::Solid,
        vec![EnthalpyPhaseKnot {temperature_k:280.0,specific_enthalpy_j_kg:0.0,liquid_mass_fraction:0.0,bulk_density_kg_m3:6.0},
            EnthalpyPhaseKnot {temperature_k:360.0,specific_enthalpy_j_kg:800.0,liquid_mass_fraction:0.0,bulk_density_kg_m3:6.0}]).unwrap();
    (card,curve)
}
fn oriented<'a>(id:u64,card:&'a MaterialCard,curve:&'a EquilibriumEnthalpyPhaseCurve,
    point:&'a QueryPoint,stiffness:[[f64;3];3],expansion:[[f64;3];3])->TetCellMaterial<'a> {
    TetCellMaterial {
        stiffness_model:ThermalStiffnessModel::OrientedOrthotropic {
            principal_to_mesh:stiffness,orientation_identity:ContentHash([6;32])},
        expansion_model:ThermalExpansionModel::OrientedOrthotropic {
            principal_to_mesh:expansion,principal_frame:ContentHash([7;32])},
        ..source(id,card,curve,point)
    }
}
fn rotate(x:[f64;3],q:[[f64;3];3])->[f64;3] {
    core::array::from_fn(|i|(0..3).map(|j|q[i][j]*x[j]).sum())
}
fn tensor_rotate(v:[f64;6],q:[[f64;3];3])->[f64;6] {
    let r=core::f64::consts::FRAC_1_SQRT_2;
    let a=[[v[0],r*v[3],r*v[5]],[r*v[3],v[1],r*v[4]],[r*v[5],r*v[4],v[2]]];
    let b:[[f64;3];3]=core::array::from_fn(|i|core::array::from_fn(|j|
        (0..3).map(|p|(0..3).map(|k|q[i][p]*a[p][k]*q[j][k]).sum::<f64>()).sum()));
    [b[0][0],b[1][1],b[2][2],b[0][1]/r,b[1][2]/r,b[2][0]/r]
}

// Independent principal nu=0 Helmholtz potential with integrated affine alpha.
fn element(t:f64,u:f64,scale:f64)->([f64;3],f64) {
    let delta=t-300.0;
    let e=core::array::from_fn::<_,3,_>(|i|E[i]+ET[i]*delta);
    let free=core::array::from_fn::<_,3,_>(|i|scale*(ALPHA[i]*delta+0.5*ALPHAT[i]*delta*delta));
    let elastic=[u-free[0],-free[1],-free[2]];
    let stress=core::array::from_fn(|i|e[i]*elastic[i]);
    let psi=0.5*VOLUME*(0..3).map(|i|e[i]*elastic[i]*elastic[i]).sum::<f64>();
    let entropy=VOLUME*(0..3).map(|i|scale*(ALPHA[i]+ALPHAT[i]*delta)*stress[i]
        -0.5*ET[i]*elastic[i]*elastic[i]).sum::<f64>();
    (stress,psi+t*entropy)
}

#[test]
fn finite_heating_material_swaps_and_rotated_constraints_use_original_anisotropic_stress() {
    with_cx(|cx| {
        let mut constrained_temperatures=Vec::new();
        for scale in [1.0,2.0] {for clamped in [false,true] {
            let (card,curve)=fixture(scale,true,false);
            let heat=thermal(&[&curve],&[10.0],&[200.0],0.0);let context=context(300.0);
            let expected_t=if clamped {root(280.0,359.0,|t|10.0*(t-300.0)+element(t,0.0,scale).1-20.0)} else {302.0};
            let mut original_stress=None;
            for q in [IDENTITY,TURN] {
                let nodes=NODES.map(|x|rotate(x,q));
                let source=SolidEnthalpy::try_new(&heat,&nodes,&TETS,
                    vec![oriented(30,&card,&curve,&context,q,q)],TetAssemblyBudget::standard(),
                    0.05,mechanical().mesh_frame,cx).unwrap();
                let loads=source.geometry().assemble(&[],VolumeLoad3::None,&[],cx).unwrap();
                // Free rigid-mode pins are defined in the global frame. They
                // remove only rigid motion; the compatible free strain is read
                // from the original tensor stress report, not an assumed gauge.
                let fixed=if clamped {(0..12).collect::<Vec<_>>()} else {FREE.to_vec()};
                let body=ThermoelasticSolid::try_new(&source,&heat,&loads,&fixed,&[[280.0,360.0]],
                    mechanical(),solve(),cx).unwrap();
                let next=body.advance(&heat,&(),0.1,&loads,thermal_controls(),solve(),cx).unwrap();
                near(next.thermal.temperature_k(IDS[0]).unwrap(),expected_t,3e-8);
                near(next.report.external_heat_j,20.0,1e-9);
                near(next.report.combined_residual_j,0.0,1e-9);
                let stress=next.solid.equilibrium().solution.stress.elements()[0].stress_mandel_pa();
                let principal=if clamped {let sigma=element(expected_t,0.0,scale).0;
                    [sigma[0],sigma[1],sigma[2],0.0,0.0,0.0]} else {[0.0;6]};
                for (actual,expected) in stress.iter().zip(tensor_rotate(principal,q)) {near(*actual,expected,2e-7);}
                if let Some(first)=original_stress {for (actual,expected) in stress.iter().zip(tensor_rotate(first,q)) {near(*actual,expected,2e-7);}}
                else {original_stress=Some(*stress);}
                near(next.thermal.energy_j()-heat.energy_j()+next.report.internal_energy_change_j,20.0,1e-9);
                assert_eq!(next.solid.state().samples()[0].elastic.card_identity(),card.content_hash());
                near(next.solid.state().samples()[0].reference_mass_kg,1.0,1e-14);
                let acoustic=next.solid.adiabatic_state(cx).unwrap();
                assert_eq!(acoustic.responses()[0].isothermal_stiffness_mandel_pa(),
                    *next.solid.state().materials()[0].stiffness_mandel_pa());
                if clamped && q==IDENTITY {constrained_temperatures.push(expected_t);}
            }
        }}
        assert!((constrained_temperatures[1]-constrained_temperatures[0]).abs()>1e-4);
    });
}

#[test]
fn original_inertial_trajectory_and_finite_heat_match_independent_orthotropic_endpoint_equations() {
    with_cx(|cx| {
        let (card,curve)=fixture(1.0,true,false);let heat=thermal(&[&curve],&[10.0],&[50.0],0.0);
        let context=context(300.0);
        let source=SolidEnthalpy::try_new(&heat,&NODES,&TETS,
            vec![oriented(30,&card,&curve,&context,IDENTITY,IDENTITY)],TetAssemblyBudget::standard(),
            0.05,mechanical().mesh_frame,cx).unwrap();
        let initial=source.geometry().assemble(&[],VolumeLoad3::None,&[(3,0.8)],cx).unwrap();
        let applied=source.geometry().assemble(&[],VolumeLoad3::None,&[(3,1.2)],cx).unwrap();
        let fixed:Vec<_>=(0..12).filter(|d|*d!=3).collect();
        let body=ThermoelasticSolid::try_new(&source,&heat,&initial,&fixed,&[[280.0,360.0]],mechanical(),solve(),cx).unwrap();
        let u0=0.8/(VOLUME*E[0]);let velocity=0.08;let dt:f64=0.02;
        let evaluate=|t:f64| {
            let delta=t-300.0;let k0=VOLUME*E[0];let k1=VOLUME*(E[0]+ET[0]*delta);
            let free=ALPHA[0]*delta+0.5*ALPHAT[0]*delta*delta;
            let u=((0.2/dt.powi(2)-0.5*k0)*u0+0.2/dt*velocity+1.0+0.5*k1*free)
                /(0.2/dt.powi(2)+0.5*k1);
            let before=element(300.0,u0,1.0);let after=element(t,u,1.0);
            let work=0.5*VOLUME*(before.0[0]+after.0[0])*(u-u0);
            let q=work-(after.1-before.1);
            (10.0*(t-300.0)-50.0*dt-q,[u,2.0*(u-u0)/dt-velocity,q])
        };
        let expected_t=root(280.0,359.0,|t|evaluate(t).0);let expected=evaluate(expected_t).1;
        let body=body.begin_transient(&initial,&[[0.0;3],[velocity,0.0,0.0],[0.0;3],[0.0;3]],
            TransientConfig3 {maximum_steps:100,maximum_displacement_m:0.05,
                absolute_force_tolerance_n:1e-11,absolute_energy_tolerance_j:1e-10,
                relative_energy_tolerance:1e-8,..TransientConfig3::default()},cx).unwrap();
        let next=body.advance(&heat,&(),dt,&applied,thermal_controls(),solve(),cx).unwrap();
        let motion=next.solid.nodal_motion(cx).unwrap();
        near(next.thermal.temperature_k(IDS[0]).unwrap(),expected_t,3e-8);
        near(motion[0][1][0],expected[0],2e-10);near(motion[1][1][0],expected[1],2e-9);
        near(next.report.reversible_heat_j,expected[2],2e-9);
        near(next.report.combined_residual_j,0.0,1e-9);
        near(next.report.kinetic_energy_change_j,0.05*(expected[1]*expected[1]-velocity*velocity),2e-10);
        assert_eq!(body.step(),0);assert_eq!(next.solid.step(),1);assert_eq!(heat.time_s(),0.0);
        next.solid.adiabatic_prestressed_state(cx).unwrap();
        let observer=Steps {count:0,cancel_at:1};
        assert!(body.advance(&heat,&observer,dt,&applied,thermal_controls(),solve(),cx).is_err());
        assert_eq!(body.step(),0);assert_eq!(heat.time_s(),0.0);assert_eq!(observer.count,0);
    });
}

#[test]
fn nonuniform_finite_heating_bends_the_actual_anisotropic_supported_mesh() {
    with_cx(|cx| {
        let nodes=[[0.,0.,0.],[1.,0.,0.],[1.,1.,0.],[0.,1.,0.],
            [0.,0.,1.],[1.,0.,1.],[1.,1.,1.],[0.,1.,1.]];
        let tets=[[0,1,2,6],[0,2,3,6],[0,3,7,6],[0,7,4,6],[0,4,5,6],[0,5,1,6]];
        let (card,curve)=fixture(1.0,false,true);let context=context(300.0);
        let fixed:Vec<_>=[0,3,4,7].into_iter().flat_map(|node|(0..3).map(move |k|3*node+k)).collect();
        let mut deflections=Vec::new();
        for gradient in [-1.0,0.0,1.0] {
            let events:Vec<_>=tets.iter().enumerate().map(|(i,tet)| {
                let z=tet.iter().map(|&n|nodes[n][2]).sum::<f64>()/4.0;
                TimedEvent {time_s:0.0,event:Event::Birth(DepositedCell {id:CellId(100+i as u64),
                    segment_id:i as u64,body_id:3,material_id:MaterialId(0),mass_kg:1.0,specific_h_j_kg:200.0,
                    geometry:CellGeometry {start_m:[0.0;3],end_m:[1.0,0.0,0.0],volume_m3:VOLUME,surface_area_m2:3.0},
                    boundary:Boundary {ambient_k:300.0,htc_w_m2_k:0.0,reservoir:None,source_w:10.0*(1.0+gradient*(z-0.5))}})}
            }).collect();
            let heat=DepositionSimulation::new(events,&[(MaterialId(0),&curve as &dyn ThermalLaw)],
                DepositionBudget {max_cells:6,max_links:0,max_events:6,max_substeps:1000}).unwrap()
                .advance_to(0.0,thermal_controls(),||Ok(())).unwrap().0;
            let cells=(0..6).map(|i|oriented(100+i,&card,&curve,&context,IDENTITY,IDENTITY)).collect();
            let source=SolidEnthalpy::try_new(&heat,&nodes,&tets,cells,TetAssemblyBudget::standard(),
                0.05,mechanical().mesh_frame,cx).unwrap();
            let loads=source.geometry().assemble(&[],VolumeLoad3::None,&[],cx).unwrap();
            let body=ThermoelasticSolid::try_new(&source,&heat,&loads,&fixed,&[[280.0,360.0];6],mechanical(),solve(),cx).unwrap();
            let next=body.advance(&heat,&(),1.0,&loads,thermal_controls(),solve(),cx).unwrap();
            let u=next.solid.equilibrium().solution.displacement.displacement_m();
            let deflection=[1,2,5,6].into_iter().map(|i|u[i][2]).sum::<f64>()/4.0;
            deflections.push(deflection);
            near(next.report.external_heat_j,60.0,1e-8);near(next.report.combined_residual_j,0.0,1e-9);
            if gradient==0.0 {near(deflection,0.0,1e-10);}
            else {
                assert!(next.report.cells.iter().any(|cell|cell.reversible_heat_j.abs()>1e-7));
                assert!(next.solid.equilibrium().solution.stress.elements().iter().any(|cell|
                    cell.stress_mandel_pa().iter().any(|s|s.abs()>1e-3)));
            }
        }
        assert!(deflections[0]>1e-6 && deflections[2]<-1e-6);
        near(deflections[0],-deflections[2],1e-9);
        // Independent dense B^T C B assembly and six-temperature Newton
        // balance (not the runtime's heat fixed point) gives these values.
        // Its maximum local energy residual is 2.72e-13 J. This checks the
        // resolved six-tet fixture, not a continuum beam-convergence claim.
        near(deflections[0],1.0920071031292661e-4,3e-11);
        near(deflections[2],-1.0920069677871665e-4,3e-11);
    });
}
