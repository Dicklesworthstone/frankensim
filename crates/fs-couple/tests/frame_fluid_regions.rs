//! G1/G3: public regional-fluid path, analytic integrals and actual motion.
use fs_couple::air_path::cylinder::{CylinderFlowResolution,frame::{FrameFluidBinding,
    FrameFluidBudget,FrameFluidRegion,MaterialFluidFrame,FrameFluidError}};
use fs_couple::frame_vibration::FrameVibrationBudget;
use fs_couple::modal_acoustic_time::ModalAcousticState;
use fs_exec::CancelGate;
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet,InterpolationPolicy,MaterialCard,MaterialStateId,PropertyClaim,
    PropertyKey,PropertyValue,Provenance,QueryPoint,UncertaintyModel};
use fs_material::fluid_transport::{FluidTransport,FluidTransportState};
use fs_material::section::{ConcentricElasticSection,ConcentricMaterialLayer};
use fs_material::state_point::{MaterialPropertySelection,ScalarAdmissibility,
    ScalarPropertyRequirement,resolve_material_state_point};
use fs_qty::{Density,DynViscosity,Pressure};
use fs_solid::frame2d::{Frame2d,FrameMassFormulation,FrameViscousSection,Geometry};

fn specimen(rotations:bool)->Frame2d {
    let mut claims=ClaimSet::new();let mut requirements=Vec::new();
    for (name,value,dims) in [("density",1000.0,Density::DIMS),
        ("young_modulus",2e6,Pressure::DIMS),("eta",500.0,DynViscosity::DIMS)] {
        let mut validity=ValidityDomain::unconstrained().with("T",300.0,400.0);
        if name=="eta" {validity=validity.with("omega",0.0,1e7);}
        claims.insert_claim(PropertyClaim {key:PropertyKey::new(name,dims),
            value:PropertyValue::Scalar {value,dims},validity,
            uncertainty:UncertaintyModel::Unstated,interpolation:InterpolationPolicy::ConstantWithinValidity,
            observations:vec![],provenance:Provenance {source:"synthetic regional-fluid test".into(),
                license:"CC0-1.0".into(),artifact:None}}).unwrap();
        requirements.push(ScalarPropertyRequirement::try_new(name,dims,ScalarAdmissibility::Finite).unwrap());
    }
    let card=MaterialCard::assemble(MaterialStateId {chemistry:"manufactured beam".into(),
        phase:"solid".into(),process:"constant reference law".into(),revision:0},claims,vec![]).unwrap();
    let point=QueryPoint::new().with("T",300.0).unwrap().with("omega",1.0).unwrap();
    let material=resolve_material_state_point(&card,&point,&requirements,
        MaterialPropertySelection::SingleClaimOnly).unwrap();
    let section=ConcentricElasticSection::try_new(0.0,vec![ConcentricMaterialLayer {outer_radius_m:0.02,material}]).unwrap();
    let mut frame=Frame2d::with_mass_formulation(Geometry::Linear,FrameMassFormulation::ConsistentReference);
    frame.add_node(0.0,0.0);frame.add_node(1.0,0.0);
    frame.fix(0,[true;3]);frame.fix(1,[true,false,!rotations]);
    frame.add_material_viscous_element(0,1,3,FrameViscousSection::try_new(section,"eta","eta").unwrap()).unwrap();
    frame
}
fn fluid(rho:f64,mu:f64)->FluidTransportState {FluidTransportState::try_new(rho,4000.0,0.5,mu).unwrap()}
fn resolution()->CylinderFlowResolution {CylinderFlowResolution {
    minimum_wavenumber_radius:0.01,maximum_wavenumber_radius:10.0,intervals:3}}
fn dry_budget()->FrameVibrationBudget {FrameVibrationBudget {max_frame_dofs:6,max_dense_entries:36}}
fn budget()->FrameFluidBudget {FrameFluidBudget {max_states:64,max_dense_entries:4096,
    max_projection_terms:4096,mass_inverse_tolerance:1e-10}}
fn build<'a,F:FluidTransport>(frame:&'a Frame2d,m:usize,regions:Vec<FrameFluidRegion<F>>)
    ->Result<MaterialFluidFrame<'a,F>,FrameFluidError> {
    MaterialFluidFrame::new_with_regions(frame,m,regions,resolution(),dry_budget(),budget(),&CancelGate::new_clock_free())
}
fn region(span:[f64;2],rho:f64,mu:f64)->FrameFluidRegion<FluidTransportState> {
    FrameFluidRegion {element:0,span,fluid:fluid(rho,mu)}
}
fn near(a:f64,b:f64,tol:f64) {assert!((a-b).abs()<=tol*(1.0+a.abs()+b.abs()),"{a:e} != {b:e}");}
fn rest<F>(sim:&MaterialFluidFrame<'_,F>)->Vec<f64> {
    sim.state_from_motion(&vec![ModalAcousticState::default();sim.basis().modes().len()],
        &vec![0.0;sim.memory_masses_kg().len()]).unwrap()
}
// Analytic monomial antiderivative, independent of production Gauss nodes.
fn gram(sim:&MaterialFluidFrame<'_,FluidTransportState>,i:usize,j:usize,span:[f64;2])->f64 {
    let polynomial=|index:usize| {
        let u=&sim.basis().modes()[index].displacement;
        [0.0,0.0,3.0*u[4]-u[5],-2.0*u[4]+u[5]]
    };
    let (p,q)=(polynomial(i),polynomial(j));let mut integral=0.0;
    for a in 0..4 {for b in 0..4 {
        let power=(a+b+1) as i32;
        integral+=p[a]*q[b]*(span[1].powi(power)-span[0].powi(power))/f64::from(power);
    }}
    integral
}

#[test]
fn whole_region_keeps_existing_operators_and_material_geometry() {
    let frame=specimen(false);let source=fluid(400.0,0.1);let gate=CancelGate::new_clock_free();
    let full=MaterialFluidFrame::new(&frame,1,vec![FrameFluidBinding {element:0,fluid:source}],
        resolution(),dry_budget(),budget(),&gate).unwrap();
    let region=build(&frame,1,vec![FrameFluidRegion {element:0,span:[0.0,1.0],fluid:source}]).unwrap();
    assert_eq!(full.system().structure(),region.system().structure());
    assert_eq!(full.modal_mass_matrix(),region.modal_mass_matrix());
    assert_eq!(full.memory_masses_kg(),region.memory_masses_kg());
    assert_eq!(region.members()[0].span(),[0.0,1.0]);
    assert_eq!(region.members()[0].length_m(),1.0);
    assert_eq!(region.members()[0].flow().fluid(),&source);
    assert!(core::ptr::eq(region.frame(),&frame));
}

#[test]
fn partial_coverage_integrates_original_shapes_not_a_length_fraction() {
    let frame=specimen(false);
    let primitive=|x:f64|9.0*x.powi(5)/5.0-2.0*x.powi(6)+4.0*x.powi(7)/7.0;
    for span in [[0.0,0.5],[0.5,1.0],[0.2,0.75]] {
        let sim=build(&frame,1,vec![region(span,400.0,0.1)]).unwrap();
        let added=0.4*(primitive(span[1])-primitive(span[0]))/(13.0/35.0);
        near(sim.modal_mass_matrix()[0],1.0+added,1e-11);
        assert!((added-0.4*(span[1]-span[0])).abs()>1e-3);
        near(sim.members()[0].length_m(),span[1]-span[0],1e-14);
        assert_eq!(sim.basis().ndof(),6);
        assert_eq!(sim.members()[0].flow().radius_m(),frame.material_section(0).unwrap().outer_radius_m());
    }
}

#[test]
fn mixed_regions_keep_cross_mass_and_complex_reciprocal_force() {
    let frame=specimen(true);let spans=[[0.37,1.0],[0.0,0.37]];
    let sim=build(&frame,2,vec![region(spans[0],400.0,0.1),region(spans[1],1000.0,0.01)]).unwrap();
    assert_eq!(sim.members()[0].span(),spans[0]);assert_eq!(sim.members()[1].span(),spans[1]);
    let n=sim.system().state_dim();let (_,r,_)=sim.system().structure();let m=2;
    assert!(sim.modal_mass_matrix()[1].abs()>1e-3);
    for i in 0..m {for j in 0..m {
        let diagonal=if i==j {1.0} else {0.0};
        let expected_mass=diagonal+sim.members().iter().map(|region|
            region.flow().added_mass_per_length_kg_m()*gram(&sim,i,j,region.span())).sum::<f64>();
        near(sim.modal_mass_matrix()[i*m+j],expected_mass,1e-10);
        for omega in [1.0,10.0,100.0] {
            let mut re=r[(m+i)*n+m+j];let mut im=omega*sim.modal_mass_matrix()[i*m+j];
            if i==j {im-=sim.basis().modes()[i].eigenvalue_per_s2/omega;}
            for (k,mass) in sim.memory_masses_kg().iter().enumerate() {
                let idx=2*m+k;let c=r[idx*n+idx];let wm=omega*mass;
                let product=r[(m+i)*n+idx]*r[idx*n+m+j];
                re-=product*c/(c*c+wm*wm);im+=product*wm/(c*c+wm*wm);
            }
            // Uniform eta/E gives the exact original solid matrix in dry modes.
            let mut expected_re=if i==j {500.0/2e6*sim.basis().modes()[i].eigenvalue_per_s2} else {0.0};
            let mut expected_im=if i==j {omega-sim.basis().modes()[i].eigenvalue_per_s2/omega} else {0.0};
            for region in sim.members() {
                let z=region.flow().impedance_per_length(omega).unwrap();let g=gram(&sim,i,j,region.span());
                expected_re+=g*z.re;expected_im+=g*z.im;
            }
            near(re,expected_re,2e-9);near(im,expected_im,2e-9);
        }
    }}
}

#[test]
fn identical_partition_preserves_native_motion_and_total_fluid_loss() {
    let frame=specimen(false);let gate=CancelGate::new_clock_free();
    let full=build(&frame,1,vec![region([0.0,1.0],400.0,0.1)]).unwrap();
    let split=build(&frame,1,vec![region([0.0,0.37],400.0,0.1),region([0.37,1.0],400.0,0.1)]).unwrap();
    let (mut a,mut b)=(rest(&full),rest(&split));
    for _ in 0..40 {
        let force=[0.0,0.0,0.0,0.0,0.01,0.0];
        let ra=full.step_candidate(&a,&force,0.001,1e-9,&gate).unwrap();
        let rb=split.step_candidate(&b,&force,0.001,1e-9,&gate).unwrap();
        near(ra.native.x[0],rb.native.x[0],1e-10);near(ra.native.x[1],rb.native.x[1],1e-10);
        near(ra.fluid_member_dissipation_j[0],rb.fluid_member_dissipation_j.iter().sum(),1e-11);
        a=ra.native.x;b=rb.native.x;
    }
    assert_eq!(split.members().len(),2);assert_eq!(split.memory_masses_kg().len(),24);
}

#[test]
fn regional_sources_reach_prepared_motion_and_physical_losses() {
    let frame=specimen(false);let gate=CancelGate::new_clock_free();let mut positions=Vec::new();
    for span in [[0.0,0.5],[0.5,1.0]] {
        let sim=build(&frame,1,vec![region(span,400.0,0.1)]).unwrap();
        let mut prepared=sim.prepare(&gate).unwrap();let mut x=rest(&sim);
        for _ in 0..20 {
            let force=[0.0,0.0,0.0,0.0,0.01,0.0];
            let dense=sim.step_candidate(&x,&force,0.001,1e-9,&gate).unwrap();
            let proposal=prepared.step_candidate(&x,&force,0.001,1e-9,&gate).unwrap();
            for (a,b) in dense.native.x.iter().zip(&proposal.state) {near(*a,*b,1e-9);}
            assert!(proposal.fluid_member_dissipation_j[0]>0.0);
            near(proposal.energy_residual_j,0.0,1e-9);x=proposal.state;
        }
        positions.push(x[0]*sim.basis().modes()[0].displacement[4]);
    }
    assert!(positions[1]<0.9*positions[0]);
}

#[test]
fn overlaps_bad_spans_and_capacity_refuse_before_source_evaluation() {
    use std::cell::Cell;
    struct Counted<'a>(&'a Cell<usize>);
    impl FluidTransport for Counted<'_> {
        fn transport_state(&self)->Result<FluidTransportState,fs_material::MaterialError> {
            self.0.set(self.0.get()+1);Ok(fluid(400.0,0.1))
        }
    }
    let frame=specimen(false);let calls=Cell::new(0);let gate=CancelGate::new_clock_free();
    let r=|span|FrameFluidRegion {element:0,span,fluid:Counted(&calls)};
    for span in [[-0.1,0.5],[0.2,1.1],[0.5,0.5],[0.7,0.2],[f64::NAN,1.0],
        [0.0,f64::INFINITY],[1.0-f64::EPSILON,1.0]] {
        assert!(build(&frame,1,vec![r(span)]).is_err());
    }
    assert!(build(&frame,1,vec![r([0.5,1.0]),r([0.0,0.6])]).is_err());
    assert!(build(&frame,1,vec![r([0.0,1.0]),r([0.0,1.0])]).is_err());
    assert!(MaterialFluidFrame::new_with_regions(&frame,1,vec![r([0.0,0.5]),r([0.5,1.0])],
        resolution(),dry_budget(),FrameFluidBudget {max_states:25,..budget()},&gate).is_err());
    assert_eq!(calls.get(),0);
    let cancelled=CancelGate::new_clock_free();cancelled.request();
    assert!(matches!(MaterialFluidFrame::new_with_regions(&frame,1,vec![r([0.0,1.0])],
        resolution(),dry_budget(),budget(),&cancelled),Err(FrameFluidError::Cancelled)));
    assert_eq!(calls.get(),0);
    let sim=build(&frame,1,vec![r([0.5,1.0]),r([0.0,0.5])]).unwrap();
    assert_eq!(calls.get(),2);assert_eq!(sim.members()[0].span(),[0.5,1.0]);
    assert_eq!(frame.displacements(),&[0.0;6]);
}
