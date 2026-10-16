//! Physical spatial integrals, endpoint work and the existing coupled time owner.
use super::*;
use super::super::{ImpactConfig, ImpactSystem};
use fs_exec::CancelGate;

fn support()->TranslatingSupport {TranslatingSupport {mass_kg:0.01,stiffness_n_m:12.0,
    damping_n_s_m:0.02,initial_position_m:0.0,initial_velocity_m_s:0.0,
    maximum_travel_m:0.02,maximum_slope:0.2}}
fn wire(n:usize)->WireSpan {WireSpan {endpoints_m:[[0.0,0.0],[0.2,0.0]],linear_density_kg_m:0.01,
    tension_n:0.8,bending_stiffness_n_m2:1e-7,damping_per_s:vec![2.0;n]}}
fn law()->Option<StringStretching> {Some(StringStretching {axial_rigidity_n:100.0,maximum_slope:0.2})}
fn near(a:f64,b:f64,tol:f64) {assert!((a-b).abs()<tol*(a.abs()+b.abs()).max(1e-20),"{a} != {b}");}
fn model()->SupportedStrings {SupportedStrings::with_endpoint_motion(&[wire(3),wire(2)],
    &[law(),law()],&[StringEndMotion::Second,StringEndMotion::First],support()).unwrap()}

#[test]
fn fixed_nut_and_moving_bridge_retain_full_lift_inertia_stretch_and_endpoint_work() {
    for end in [StringEndMotion::First,StringEndMotion::Second] {for n in [1,3,12] {
        let (x,v)=(0.002,0.03);
        let source=SupportedStrings::with_endpoint_motion(&[wire(n)],&[law()],&[end],
            TranslatingSupport {initial_position_m:x,initial_velocity_m_s:v,..support()}).unwrap();
        let body=source.body();let z:Vec<_>=body.initial.iter().map(|s|s.displacement_m_sqrt_kg).collect();
        let kinetic:f64=body.initial.iter().map(|s|0.5*s.velocity_m_sqrt_kg_per_s.powi(2)).sum();
        near(kinetic,0.5*(0.01+0.002/3.0)*v*v,1e-14);
        near(source.potential(&z),0.5*12.0*x*x+0.8*x*x/(2.0*0.2)+100.0*x.powi(4)/(8.0*0.2_f64.powi(3)),1e-13);
        let observed=source.string_observation(0,&z).unwrap();
        near(observed.additional_strain,0.5*(x/0.2).powi(2),1e-12);
        near(observed.tension_n,0.8+100.0*observed.additional_strain,1e-14);
        near(observed.slope_bound,x/0.2,1e-12);
        for station in [0.0,0.23,1.0] {
            let row=source.absolute_shape(0,station).unwrap();
            let displacement:f64=row.iter().zip(&z).map(|(b,z)|b*z).sum();
            near(displacement,end.lift(station)*x,1e-12);
            let speed:f64=row.iter().zip(&body.initial).map(|(b,s)|b*s.velocity_m_sqrt_kg_per_s).sum();
            near(speed,end.lift(station)*v,1e-12);
        }
        assert_eq!(source.total_mass_kg(),0.012); // physical mass, not lift inertia
    }}
    let original=SupportedStrings::new(&[wire(3)],law(),support()).unwrap();
    let explicit=SupportedStrings::with_endpoint_motion(&[wire(3)],&[law()],&[StringEndMotion::Both],support()).unwrap();
    assert_eq!(original.initial_state(),explicit.initial_state());
    assert_eq!(original.clock_bounds(),explicit.clock_bounds());
    let z=[2e-5,-3e-5,1e-5,1e-4];
    assert_eq!(original.potential(&z),explicit.potential(&z));
}

#[test]
fn coupled_gradient_and_hessian_include_bridge_tension_and_each_wire_cross_term() {
    let s=model();let n=s.mode_count();
    let z:Vec<_>=(0..n).map(|i|2e-5*(i as f64-1.2)).collect();
    let d:Vec<_>=(0..n).map(|i|0.003*(i as f64+0.7)).collect();
    let (mut g,mut hv)=(vec![0.0;n],vec![0.0;n]);
    s.gradient(&z,&mut g);s.hessian_vector(&z,&d,&mut hv);
    for i in 0..n {
        let h=1e-9;let mut plus=z.clone();let mut minus=z.clone();plus[i]+=h;minus[i]-=h;
        near(g[i],(s.potential(&plus)-s.potential(&minus))/(2.0*h),2e-6);
    }
    let h=1e-7;let plus:Vec<_>=z.iter().zip(&d).map(|(z,d)|z+h*d).collect();
    let minus:Vec<_>=z.iter().zip(&d).map(|(z,d)|z-h*d).collect();
    let (mut a,mut b)=(vec![0.0;n],vec![0.0;n]);s.gradient(&plus,&mut a);s.gradient(&minus,&mut b);
    for i in 0..n {near(hv[i],(a[i]-b[i])/(2.0*h),2e-7);}
    // Independent midpoint spatial integration of T w'^2/2 + EI w''^2/2
    // plus EA/(8L)*(integral w'^2)^2. The lift's derivative is not zero.
    let x=z[n-1]*s.force_weight();let mut energy=0.5*support().stiffness_n_m*x*x;
    for (idx,w) in s.0.wires.iter().enumerate() {
        let l=w.length_m();let cells=8192;let dx=l/cells as f64;
        let mut slope2=0.0;let mut curvature2=0.0;
        for j in 0..cells {
            let station=(j as f64+0.5)*dx;
            let mut slope=if idx==0 {x/l} else {-x/l};let mut curvature=0.0;
            for (mode,i) in s.0.ranges[idx].clone().enumerate() {
                let k=(mode+1) as f64*core::f64::consts::PI/l;
                let amplitude=(2.0/(w.linear_density_kg_m*l)).sqrt()*(z[i]-s.0.translation[i]*x);
                slope+=amplitude*k*(k*station).cos();curvature-=amplitude*k*k*(k*station).sin();
            }
            slope2+=slope*slope*dx;curvature2+=curvature*curvature*dx;
        }
        energy+=0.5*w.tension_n*slope2+0.5*w.bending_stiffness_n_m2*curvature2+100.0/(8.0*l)*slope2*slope2;
    }
    near(s.potential(&z),energy,1e-10);
}

#[test]
fn stationary_drag_matches_spatial_absolute_velocity_loss_not_relative_bending_loss() {
    let base=model();
    let drag=vec![AbsoluteStringDrag {modal_per_s:vec![3.0;3],remainder_per_s:3.0},
        AbsoluteStringDrag {modal_per_s:vec![7.0;2],remainder_per_s:7.0}];
    let s=base.clone().with_absolute_drag(drag.clone()).unwrap();let n=s.mode_count();
    let v:Vec<_>=(0..n).map(|i|0.007*(i as f64-1.4)).collect();
    let mut a=vec![0.0;4*n*n];let mut b=a.clone();
    s.add_resistance(0,2*n,&mut a).unwrap();base.add_resistance(0,2*n,&mut b).unwrap();
    let mut power=0.0;
    for i in 0..n {for j in 0..n {power+=v[i]*(a[(2*i+1)*2*n+2*j+1]-b[(2*i+1)*2*n+2*j+1])*v[j];}}
    let speed=v[n-1]*s.force_weight();let mut spatial=0.0;
    for (idx,w) in s.0.wires.iter().enumerate() {
        let l=w.length_m();let cells=8192;let dx=l/cells as f64;
        for j in 0..cells {
            let t=(j as f64+0.5)/cells as f64;let mut local=s.0.end_motion[idx].lift(t)*speed;
            for (mode,i) in s.0.ranges[idx].clone().enumerate() {
                local+=(2.0/(w.linear_density_kg_m*l)).sqrt()*((mode+1) as f64*core::f64::consts::PI*t).sin()
                    *(v[i]-s.0.translation[i]*speed);
            }
            spatial+=drag[idx].remainder_per_s*w.linear_density_kg_m*local*local*dx;
        }
    }
    near(power,spatial,2e-8);assert!(power>0.0);
    assert!(s.clone().with_absolute_drag(drag.clone()).is_err());
    let mut invalid=drag.clone();invalid[1].modal_per_s.pop();assert!(base.clone().with_absolute_drag(invalid).is_err());
    let mut invalid=drag;invalid[0].remainder_per_s=-1.0;assert!(base.with_absolute_drag(invalid).is_err());
}

#[test]
fn contact_and_slope_admission_use_the_actual_endpoint_lift() {
    let s=model();let line=LineContact::uniform(0.2,8,0.00002,1e6,1.5,0.05,"fixture".into()).unwrap();
    let shapes=vec![vec![2.0];8];let n=s.mode_count()+1;
    let contact=s.contact(0,&line,&shapes,0..1,1,n).unwrap();
    let x=0.001;let mut q=vec![0.0;n];q[0]=0.0001;
    for i in 0..s.support_coordinate() {q[1+i]=s.0.translation[i]*x;}
    q[n-1]=x/s.force_weight();
    for (row,station) in contact.collocation().chunks_exact(n).zip(&line.stations_m) {
        let gap:f64=row.iter().zip(&q).map(|(b,q)|b*q).sum();
        near(gap,2.0*q[0]-station/0.2*x,1e-12);
    }
    let small=TranslatingSupport {initial_position_m:0.001,maximum_slope:0.001,..support()};
    assert!(SupportedStrings::with_endpoint_motion(&[wire(3)],&[law()],&[StringEndMotion::Second],small).is_err());
    assert!(SupportedStrings::new(&[wire(3)],law(),small).is_ok()); // rigid translation has no slope
    assert!(SupportedStrings::with_endpoint_motion(&[wire(3)],&[law()],&[],support()).is_err());
}

#[test]
fn one_plucked_span_drives_the_shared_bridge_and_unplayed_string_in_the_existing_owner() {
    let gate=CancelGate::new_clock_free();
    let mut end_states=Vec::new();
    for ea in [100.0,1000.0] {
        let s=SupportedStrings::with_endpoint_motion(&[wire(3),wire(2)],
            &[Some(StringStretching {axial_rigidity_n:ea,maximum_slope:0.2}),law()],
            &[StringEndMotion::Second,StringEndMotion::First],support()).unwrap();
        let mut initial=vec![ModalAcousticState::default();5];initial[0].displacement_m_sqrt_kg=1e-4;
        let body=s.body_with_relative_states(&initial).unwrap();
        let cfg=ImpactConfig {dt_s:2e-5,max_steps:512,maximum_energy_j:1.0,
            energy_absolute_tolerance_j:1e-10,energy_relative_tolerance:1e-7,maximum_generalized_force:1e4};
        let mut run=ImpactSystem::new(vec![body],vec![],vec![],vec![],cfg).unwrap().prepare_analytic().unwrap();
        let initial_energy=run.stored_energy_j();let mut dissipated=0.0;
        for _ in 0..256 {
            let f=run.step(&[0.0;6],&gate).unwrap();dissipated+=f.dissipated_energy_j;
            assert!((f.stored_energy_j+dissipated-initial_energy).abs()<1e-8);
        }
        let state=run.state().to_vec();let x=state[10]*s.force_weight();
        assert!(x.abs()>1e-12);
        assert!((state[6]-s.0.translation[3]*x).abs()>1e-12);
        let before=state.clone();assert!(run.step(&[f64::NAN;6],&gate).is_err());assert_eq!(run.state(),before);
        end_states.push(state);
    }
    assert_ne!(end_states[0],end_states[1]);
}