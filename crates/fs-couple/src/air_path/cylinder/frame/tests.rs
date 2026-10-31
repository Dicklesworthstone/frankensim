//! G1/G3: manufactured source data, native dynamics and independent physical matrices.
use super::*;
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
use fs_material::fluid_transport::FluidTransportState;
use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer};
use fs_material::state_point::{MaterialPropertySelection, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_material_state_point};
use fs_qty::{Density, DynViscosity, Pressure};
use fs_solid::frame2d::{FrameMassFormulation, FrameViscousSection, Geometry, SectionModel};

fn section(rho: f64, young: f64, radius: f64, eta: f64) -> ConcentricElasticSection {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (name, value, dims) in [("density",rho,Density::DIMS),
        ("young_modulus",young,Pressure::DIMS), ("axial_eta",eta,DynViscosity::DIMS),
        ("bending_eta",eta,DynViscosity::DIMS)] {
        let mut validity = ValidityDomain::unconstrained().with("T",300.0,400.0);
        if name.ends_with("_eta") { validity = validity.with("omega",0.0,1e7); }
        claims.insert_claim(PropertyClaim { key: PropertyKey::new(name,dims),
            value: PropertyValue::Scalar { value,dims }, validity,
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity, observations: vec![],
            provenance: Provenance { source: "manufactured fluid/frame test".into(),
                license: "CC0-1.0".into(), artifact: None },
        }).unwrap();
        requirements.push(ScalarPropertyRequirement::try_new(name,dims,ScalarAdmissibility::Finite).unwrap());
    }
    let card = MaterialCard::assemble(MaterialStateId { chemistry: "synthetic beam".into(),
        phase: "solid".into(), process: "test-only constant elastic/viscous law".into(), revision: 0 },
        claims,vec![]).unwrap();
    let point = QueryPoint::new().with("T",300.0).unwrap().with("omega",1.0).unwrap();
    let material = resolve_material_state_point(&card,&point,&requirements,
        MaterialPropertySelection::SingleClaimOnly).unwrap();
    ConcentricElasticSection::try_new(0.0,vec![ConcentricMaterialLayer { outer_radius_m: radius,material }]).unwrap()
}
fn frame(rho: f64, young: f64, eta: f64, two: bool) -> Frame2d {
    let mut frame = Frame2d::with_mass_formulation(Geometry::Linear,FrameMassFormulation::ConsistentReference);
    frame.add_node(0.0,0.0);
    frame.add_node(if two {0.5} else {1.0},0.0);
    if two { frame.add_node(1.0,0.0); }
    frame.fix(0,[true;3]);
    frame.fix(1,[true,false,!two]);
    if two { frame.fix(2,[true,false,false]); }
    for element in 0..(if two {2} else {1}) {
        let section = section(rho,young,0.02,eta);
        let section = FrameViscousSection::try_new(section,"axial_eta","bending_eta").unwrap();
        frame.add_material_viscous_element(element,element+1,3,section).unwrap();
    }
    frame
}
fn fluid(rho: f64, mu: f64) -> FluidTransportState {
    // Deliberately authored numerical transport; no measured-liquid claim.
    FluidTransportState::try_new(rho,4000.0,0.5,mu).unwrap()
}
fn budget() -> FrameFluidBudget {
    FrameFluidBudget { max_states: 64,max_dense_entries: 4096,max_projection_terms: 4096,
        mass_inverse_tolerance: 1e-10 }
}
fn frame_budget() -> FrameVibrationBudget {
    FrameVibrationBudget { max_frame_dofs: 12,max_dense_entries: 144 }
}
fn resolution() -> CylinderFlowResolution {
    CylinderFlowResolution { minimum_wavenumber_radius: 0.01,maximum_wavenumber_radius: 10.0,intervals: 3 }
}
fn model(frame: &Frame2d, modes: usize, member: usize, rho: f64, mu: f64)
    -> MaterialFluidFrame<'_,FluidTransportState> {
    MaterialFluidFrame::new(frame,modes,vec![FrameFluidBinding {element:member,fluid:fluid(rho,mu)}],
        resolution(),frame_budget(),budget(),&CancelGate::new_clock_free()).unwrap()
}
fn near(a: f64,b: f64,tolerance: f64) {
    assert!((a-b).abs() <= tolerance*(1.0+a.abs()+b.abs()),"{a:e} != {b:e}");
}
fn zero(model: &MaterialFluidFrame<'_,FluidTransportState>) -> Vec<f64> {
    model.state_from_motion(&vec![ModalAcousticState::default();model.basis().modes().len()],
        &vec![0.0;model.memory_masses_kg().len()]).unwrap()
}

#[test]
fn original_material_geometry_drives_all_added_inertia_and_fluid_sources_survive() {
    let frame = frame(1000.0,2e6,500.0,false);
    let sim = model(&frame,1,0,400.0,0.1);
    assert!(core::ptr::eq(sim.frame(),&frame));
    assert_eq!(sim.members()[0].element(),0);
    assert_eq!(sim.members()[0].flow().radius_m(),frame.material_section(0).unwrap().outer_radius_m());
    assert_eq!(sim.members()[0].flow().fluid(),&fluid(400.0,0.1));
    // All transverse mass comes from one homogeneous circular member. Fluid
    // and solid have the SAME Hermite spatial inertia, so the ratio is exact.
    near(sim.modal_mass_matrix()[0],1.0+400.0/1000.0,2e-10);
    assert_eq!(sim.memory_masses_kg().len(),12);
    assert_eq!(sim.system().state_dim(),14);
    assert!(sim.mass_inverse_residual() < 1e-12);
}

#[test]
fn partial_immersion_keeps_full_modal_added_mass_not_independent_retuning() {
    let frame = frame(1000.0,2e6,500.0,true);
    let sim = model(&frame,4,1,400.0,0.1);
    let length = 0.5;
    let line_mass = 400.0*core::f64::consts::PI*0.02_f64.powi(2);
    let l = length;
    let local = [[156.0,22.0*l,54.0,-13.0*l], [22.0*l,4.0*l*l,13.0*l,-3.0*l*l],
        [54.0,13.0*l,156.0,-22.0*l], [-13.0*l,-3.0*l*l,-22.0*l,4.0*l*l]];
    let dofs = [4,5,7,8];
    let mut off_diagonal = 0.0_f64;
    for i in 0..4 { for j in 0..4 {
        let mut added = 0.0;
        for a in 0..4 { for b in 0..4 {
            added += sim.basis().modes()[i].displacement[dofs[a]]*local[a][b]
                *sim.basis().modes()[j].displacement[dofs[b]]*line_mass*length/420.0;
        }}
        near(sim.modal_mass_matrix()[i*4+j],added+if i==j {1.0} else {0.0},5e-9);
        if i!=j {off_diagonal=off_diagonal.max(added.abs());}
    }}
    assert!(off_diagonal>1e-3);
    let motion: Vec<_> = (0..4).map(|i| ModalAcousticState {
        displacement_m_sqrt_kg: (i+1) as f64*1e-5,
        velocity_m_sqrt_kg_per_s: (i+1) as f64*-2e-4 }).collect();
    let state=sim.state_from_motion(&motion,&[0.01;12]).unwrap();
    let recovered=sim.modal_states(&state).unwrap();
    for (a,b) in motion.iter().zip(recovered) {
        near(a.displacement_m_sqrt_kg,b.displacement_m_sqrt_kg,1e-12);
        near(a.velocity_m_sqrt_kg_per_s,b.velocity_m_sqrt_kg_per_s,1e-12);
    }
}

#[test]
fn coupled_native_history_has_separate_solid_fluid_loss_and_actual_nodal_work() {
    let frame=frame(1000.0,2e6,500.0,false);
    let before=frame.displacements().to_vec();
    let sim=model(&frame,1,0,400.0,0.1);
    let mut state=sim.state_from_motion(&[ModalAcousticState {
        displacement_m_sqrt_kg:1e-5,velocity_m_sqrt_kg_per_s:2e-4}],&[1e-4;12]).unwrap();
    let initial=sim.system().hamiltonian(&state);
    let (mut work,mut solid,mut liquid)=(0.0,0.0,0.0);
    for step in 0..50 {
        let force=if step<20 {0.01} else {-0.003};
        let proposal=sim.step_candidate(&state,&[0.0,0.0,0.0,0.0,force,0.0],0.001,1e-9,
            &CancelGate::new_clock_free()).unwrap();
        let delta=sim.basis().modes()[0].displacement[4]*(proposal.native.x[0]-state[0]);
        near(proposal.nodal_work_j,force*delta,1e-12);
        work+=proposal.nodal_work_j;
        solid+=proposal.material_port_dissipation_j.iter().sum::<f64>();
        liquid+=proposal.fluid_member_dissipation_j.iter().sum::<f64>();
        state=proposal.native.x;
    }
    assert!(solid>0.0 && liquid>0.0);
    near(sim.system().hamiltonian(&state)-initial+solid+liquid,work,1e-9);
    assert_eq!(frame.displacements(),before.as_slice());
    assert!(state[2..].iter().any(|p| *p!=0.0));
}

#[test]
fn complex_elimination_recovers_original_impedance_including_cross_mode_terms() {
    let frame=frame(1000.0,2e6,500.0,true);
    let sim=model(&frame,4,1,400.0,0.1);
    let m=4;
    let n=sim.system().state_dim();
    let (_,r,_)=sim.system().structure();
    for omega in [1.0,10.0,100.0] {
        for i in 0..m { for j in 0..m {
            let mut re=r[(m+i)*n+m+j];
            let mut im=omega*sim.modal_mass_matrix()[i*m+j];
            if i==j { im-=sim.basis().modes()[i].eigenvalue_per_s2/omega; }
            for k in 0..sim.memory_masses_kg().len() {
                let slot=2*m+k;
                let c=r[slot*n+slot];
                let wm=omega*sim.memory_masses_kg()[k];
                let product=r[(m+i)*n+slot]*r[slot*n+m+j];
                re-=product*c/(c*c+wm*wm);
                im+=product*wm/(c*c+wm*wm);
            }
            let mut expected_re=0.0;
            let mut expected_im=if i==j {omega-sim.basis().modes()[i].eigenvalue_per_s2/omega} else {0.0};
            for port in &sim.material_ports {
                expected_re+=port.resistance*port.projection[i]*port.projection[j];
            }
            for station in &sim.stations {
                let z=sim.members[station.member].flow.impedance_per_length(omega).unwrap();
                let factor=station.weight_m*station.projection[i]*station.projection[j];
                expected_re+=factor*z.re;
                expected_im+=factor*z.im;
            }
            near(re,expected_re,2e-9);
            near(im,expected_im,2e-9);
        }}
    }
}

#[test]
fn changing_only_source_density_modulus_or_fluid_changes_the_actual_motion() {
    let frames=[frame(1000.0,2e6,500.0,false),frame(2000.0,2e6,500.0,false),
        frame(1000.0,8e6,500.0,false)];
    let mut response=Vec::new();
    for (source,rho,mu) in [(&frames[0],400.0,0.1),(&frames[1],400.0,0.1),
        (&frames[2],400.0,0.1),(&frames[0],1600.0,0.1),(&frames[0],400.0,1.0)] {
        let sim=model(source,1,0,rho,mu);
        let mut state=zero(&sim);
        for _ in 0..20 {
            state=sim.step_candidate(&state,&[0.0,0.0,0.0,0.0,0.01,0.0],0.001,1e-9,
                &CancelGate::new_clock_free()).unwrap().native.x;
        }
        response.push(state[0]*sim.basis().modes()[0].displacement[4]);
    }
    assert!(response.iter().all(|x| x.is_finite() && *x>0.0));
    for changed in &response[1..] { assert!((changed-response[0]).abs()>1e-12); }
    assert!(response[1]<response[0] && response[3]<response[0]);
}

#[test]
fn dry_limit_preserves_material_damping_and_refusals_preserve_all_input_state() {
    let frame=frame(1000.0,2e6,500.0,false);
    let gate=CancelGate::new_clock_free();
    let dry=MaterialFluidFrame::<FluidTransportState>::new(&frame,1,vec![],resolution(),frame_budget(),budget(),&gate).unwrap();
    assert_eq!(dry.modal_mass_matrix(),&[1.0]);
    assert_eq!(dry.system().state_dim(),2);
    assert_eq!(dry.memory_masses_kg().len(),0);
    let state=dry.state_from_motion(&[ModalAcousticState {displacement_m_sqrt_kg:1e-4,
        velocity_m_sqrt_kg_per_s:1e-3}],&[]).unwrap();
    let proposal=dry.step_candidate(&state,&[0.0;6],0.001,1e-9,&gate).unwrap();
    assert!(proposal.material_port_dissipation_j.iter().sum::<f64>()>0.0);
    assert!(proposal.fluid_member_dissipation_j.is_empty());
    let saved=state.clone();
    assert!(dry.step_candidate(&state,&[0.0;6],1e-9,1e-9,&gate).is_err());
    assert!(dry.step_candidate(&state,&[f64::NAN;6],0.001,1e-9,&gate).is_err());
    gate.request();
    assert!(matches!(dry.step_candidate(&state,&[0.0;6],0.001,1e-9,&gate),Err(FrameFluidError::Cancelled)));
    assert_eq!(state,saved);
    let bindings=|| vec![FrameFluidBinding {element:0,fluid:fluid(400.0,0.1)}];
    for bad in [FrameFluidBudget {max_states:13,..budget()},
        FrameFluidBudget {max_dense_entries:195,..budget()},
        FrameFluidBudget {max_projection_terms:0,..budget()}] {
        assert!(MaterialFluidFrame::new(&frame,1,bindings(),resolution(),frame_budget(),bad,
            &CancelGate::new_clock_free()).is_err());
    }
    let mut duplicate=bindings(); duplicate.extend(bindings());
    assert!(MaterialFluidFrame::new(&frame,1,duplicate,resolution(),frame_budget(),budget(),
        &CancelGate::new_clock_free()).is_err());
    let sim=model(&frame,1,0,400.0,0.1);
    assert!(sim.state_from_motion(&[ModalAcousticState::default()],&[0.0;11]).is_err());
    let mut raw=Frame2d::new(Geometry::Linear);
    raw.add_node(0.0,0.0); raw.add_node(1.0,0.0);
    raw.fix(0,[true;3]); raw.fix(1,[true,false,true]);
    raw.add_element(0,1,3,1.0,&|| SectionModel::Elastic {ea:200.0,ei:8.0}).unwrap();
    assert!(MaterialFluidFrame::new(&raw,1,bindings(),resolution(),frame_budget(),budget(),
        &CancelGate::new_clock_free()).is_err());
}

mod update;
