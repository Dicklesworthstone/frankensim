use super::*;
use fs_couple::bernoulli_aperture::plate::closure::PlateClosureSpec;
use fs_couple::bernoulli_aperture::dynamic::{DynamicAperture,ApertureState,ApertureDrive};
use fs_couple::bernoulli_aperture::tube::{ApertureTube,UniformTubeSpec,TubeDrive,TubeFrame};

fn profile(r: &PlateApertureReduction) -> PlateClosureSpec {
    PlateClosureSpec {
        nodal_rest_gap_m:r.chart().mesh.nodes.iter().map(|p| r.options().rest_opening_m*(0.05+1.9*p.1/0.01)).collect(),
        lay_triangles:(0..r.chart().mesh.tris.len()).collect(),
        stiffness_pa_per_m_alpha:1e12,alpha:2.0,internal_loss_s_per_m:0.5,
        provenance:"illustrative asymmetric plate lay; not measured cane".into(),max_penetration_m:0.0002,
    }
}
fn valve(r: PlateApertureReduction, p: PlateClosureSpec, state: ApertureState) -> DynamicAperture {
    DynamicAperture::from_plate_with_closure(r,p,1.2,1e6,1e-5,2048,state).unwrap()
}
fn rest(r: &PlateApertureReduction) -> ApertureState {
    ApertureState { opening_m:r.options().rest_opening_m,opening_velocity_m_s:0.0 }
}
fn normal_force(law: &fs_dcontact::Obstacle,q0:f64,q1:f64,vm:f64) -> f64 {
    law.collocation().iter().zip(law.gaps()).zip(law.weights()).map(|((&b,&c),&w)| {
        let (a,z) = ((b*q0-c).max(0.0),(b*q1-c).max(0.0));
        let e = if b==0.0 || q0==q1 { w*law.stiffness()*a*a }
            else if a>0.0 && z>0.0 { w*law.stiffness()*(a*a+a*z+z*z)/3.0 }
            else { w*law.stiffness()*(z.powi(3)-a.powi(3))/(3.0*b*(q1-q0)) };
        -b*(e*(1.0+law.internal_loss()*b*vm)).max(0.0)
    }).sum()
}
// Prescribe a next mechanical state, then independently derive the mouth
// pressure needed by momentum and characteristic continuity for that state.
fn manufactured(model:&DynamicAperture,next:f64) -> ApertureDrive {
    let s=model.spec();let old=model.state();let dt=s.time_step_s;
    let vm=(next-old.opening_m)/dt;let v1=2.0*vm-old.opening_velocity_m_s;
    let r=model.plate_reduction().unwrap();let area=r.pressure_area_m2();
    let k=r.stiffness_n_m();let m=r.mass_kg();let damping=2.0*s.damping_ratio*(k*m).sqrt();
    let force=normal_force(model.contact_law(),old.opening_m,next,vm);
    let dp=(force-m*(v1-old.opening_velocity_m_s)/dt-k*(f64::midpoint(old.opening_m,next)-r.options().rest_opening_m)-damping*vm)/area;
    let gap_area=model.plate_closure().unwrap().open_area_m2(f64::midpoint(old.opening_m,next)).unwrap();
    let jet=gap_area*dp.signum()*(2.0*dp.abs()/s.density_kg_m3).sqrt();
    let (incoming,body)=(75.0,2e-7);
    let outgoing=incoming+s.impedance_pa_s_m3*(jet-area*vm+body);
    ApertureDrive { upstream_pressure_pa:dp+outgoing+incoming,incoming_pressure_pa:incoming,body_flow_m3_s:body }
}
#[test]
fn physical_triangle_area_and_signed_gap_map_enter_the_original_contact_potential() {
    let r=reduction(4e9,900.0);let p=profile(&r);let c=r.compile_closure(p.clone()).unwrap();
    near(c.contact_law().weights().iter().sum(),0.025*0.01,1e-12);
    assert_eq!(c.contact_law().provenance(),p.provenance);
    let q=0.5*r.options().rest_opening_m;
    for (i,&node) in c.contact_nodes().iter().enumerate() {
        let expected=p.nodal_rest_gap_m[node]+r.shape_per_opening()[node][0]*(q-r.options().rest_opening_m);
        near(c.nodal_gap_m(node,q).unwrap(),expected,1e-12);
        let penetration=c.contact_law().collocation()[i]*q-c.contact_law().gaps()[i];
        assert!((penetration+expected).abs()<1e-18);
    }
    let observation=c.probe(q).unwrap();
    assert!(observation.active_lay_points>0 && observation.active_lay_points<c.contact_nodes().len());
    assert!(observation.active_lay_area_m2>0.0 && observation.active_lay_area_m2<0.025*0.01);
}
#[test]
fn partial_slit_closure_uses_spatial_open_area_not_positive_mean_or_clipped_endpoints() {
    let r=reduction(4e9,900.0);let c=r.compile_closure(profile(&r)).unwrap();
    let q=0.0; // zero mean-coordinate opening can leave an asymmetric slit open.
    let exact=c.open_area_m2(q).unwrap();assert!(exact>0.0);
    let mut quadrature=0.0;let mut clipped_trapezoid=0.0;
    for &[a,b] in &r.options().slit_edges {
        let (x,y)=(r.chart().mesh.nodes[a],r.chart().mesh.nodes[b]);let l=(x.0-y.0).hypot(x.1-y.1);
        let (ga,gb)=(c.nodal_gap_m(a,q).unwrap(),c.nodal_gap_m(b,q).unwrap());
        clipped_trapezoid+=l*0.5*(ga.max(0.0)+gb.max(0.0));
        for i in 0..20000 { let t=(i as f64+0.5)/20000.0;quadrature+=l*((1.0-t)*ga+t*gb).max(0.0)/20000.0; }
    }
    near(exact,quadrature,2e-8);
    assert!(clipped_trapezoid>=exact);
}
#[test]
fn spatial_flow_and_distributed_contact_share_the_implicit_pressure_motion_solve() {
    let r=reduction(4e9,900.0);let p=profile(&r);
    let initial=ApertureState { opening_m:0.5*r.options().rest_opening_m,opening_velocity_m_s:-0.05 };
    let mut model=valve(r,p,initial);let target=initial.opening_m-1e-6;
    let input=manufactured(&model,target);
    let frame=model.step(input).unwrap();
    assert!((frame.state.opening_m-target).abs()<1e-12);
    assert!((frame.state.opening_velocity_m_s-(-0.15)).abs()<1e-7);
    let profile=model.plate_closure().unwrap();
    let area=profile.open_area_m2(f64::midpoint(initial.opening_m,frame.state.opening_m)).unwrap();
    near(frame.midpoint_opening_m*model.spec().aperture.width_m,area,1e-12);
    assert!(profile.probe(frame.state.opening_m).unwrap().active_lay_points>0);
    assert!(frame.flow_residual_m3_s.abs()<1e-12);
    let scale=frame.stored_energy_j+frame.storage_change_j.abs()+frame.dissipated_energy_j+frame.pressure_work_j.abs();
    assert!(frame.balance_residual_j().abs()<=3e-10*scale);
}
#[test]
fn penetration_domain_failure_leaves_all_state_available_for_exact_retry() {
    let r=reduction(4e9,900.0);let mut p=profile(&r);p.max_penetration_m=1e-7;
    let initial=rest(&r);let mut a=valve(r.clone(),p.clone(),initial);let mut b=valve(r,p,initial);
    let target=-initial.opening_m;
    let rejected=manufactured(&a,target);
    let error=a.step(rejected).unwrap_err();
    assert!(matches!(error,fs_couple::acoustic_realize::AcousticRealizeError::InvalidDescription {
        what:"plate closure exceeds its declared penetration allowance"
    }),"{error:?}");
    assert_eq!(a.state(),initial);assert_eq!(a.accepted_steps(),0);
    let gentle=ApertureDrive { upstream_pressure_pa:0.01,..ApertureDrive::default() };
    for _ in 0..32 { assert_eq!(a.step(gentle).unwrap(),b.step(gentle).unwrap()); }
}
#[test]
fn incomplete_profile_and_duplicated_physical_area_never_acquire_defaults() {
    let r=reduction(4e9,900.0);let base=profile(&r);
    let mut cases=Vec::new();
    let mut p=base.clone();p.nodal_rest_gap_m.pop();cases.push(p);
    let mut p=base.clone();p.nodal_rest_gap_m[0]=f64::NAN;cases.push(p);
    let mut p=base.clone();p.lay_triangles[0]=usize::MAX;cases.push(p);
    let mut p=base.clone();p.lay_triangles[1]=p.lay_triangles[0];cases.push(p);
    let mut p=base.clone();p.max_penetration_m=0.0;cases.push(p);
    let mut p=base.clone();p.stiffness_pa_per_m_alpha=-1.0;cases.push(p);
    let mut p=base;p.provenance.clear();cases.push(p);
    for p in cases { assert!(r.compile_closure(p).is_err()); }
}
#[test]
fn tube_feedback_keeps_profiled_contact_and_survives_cancelled_blocks() {
    let r=reduction(4e9,900.0);let p=profile(&r);let initial=rest(&r);
    let ts=UniformTubeSpec { length_m:343.0*1e-5*32.0,radius_m:0.007,sound_speed_m_s:343.0,
        terminal_reflection:-0.8,max_length_error_m:1e-12,max_wave_memory_bytes:1<<20 };
    let build=|| {
        let valve=DynamicAperture::from_plate_with_closure(r.clone(),p.clone(),1.2,
            ts.characteristic_impedance(1.2).unwrap(),1e-5,513,initial).unwrap();
        ApertureTube::new(valve,ts).unwrap()
    };
    let (mut a,mut b)=(build(),build());
    let inputs:Vec<_>=(0..513).map(|i| TubeDrive { upstream_pressure_pa:if i<128 {5.0}else{0.0},body_flow_m3_s:0.0 }).collect();
    let expected:Vec<_>=inputs.iter().map(|&d|a.step(d).unwrap()).collect();
    let mut actual=vec![TubeFrame::default();513];
    let cancelled=CancelGate::new_clock_free();cancelled.request();
    let progress=b.advance_block(&inputs,&mut actual,&cancelled).unwrap();
    assert_eq!(progress.completed,0);assert_eq!(b.aperture().state(),initial);
    let active=CancelGate::new_clock_free();
    for (input,out) in inputs.chunks(37).zip(actual.chunks_mut(37)) { b.advance_block(input,out,&active).unwrap(); }
    assert_eq!(actual,expected);
    assert!(b.aperture().plate_closure().is_some());assert!(b.aperture().plate_reduction().is_some());
}
