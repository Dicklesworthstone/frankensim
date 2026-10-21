//! Independent spatial mechanics and real axial work on a common termination.
use super::*;
use super::super::{ImpactConfig,ImpactSystem};
use fs_exec::CancelGate;

fn wire(n:usize)->WireSpan {WireSpan {endpoints_m:[[0.0,0.0],[0.2,0.0]],linear_density_kg_m:0.012,
    tension_n:0.8,bending_stiffness_n_m2:1e-7,damping_per_s:vec![2.0;n]}}
fn rail()->TranslatingSupport {TranslatingSupport {mass_kg:0.01,stiffness_n_m:12.0,damping_n_s_m:0.02,
    initial_position_m:0.0,initial_velocity_m_s:0.0,maximum_travel_m:0.01,maximum_slope:0.2}}
fn law()->Option<StringStretching> {Some(StringStretching {axial_rigidity_n:100.0,maximum_slope:0.2})}
fn axial(modes:Vec<usize>)->LongitudinalSupportModes {LongitudinalSupportModes {
    solid_linear_density_kg_m:0.01,modes,maximum_strain:0.02}}
fn geometry(direction:[f64;3])->StringSupportGeometry {
    StringSupportGeometry::try_new([[0.0;3],[0.2,0.0,0.0]],[0.0,1.0,0.0],direction).unwrap()
}
fn near(a:f64,b:f64,tol:f64) {assert!((a-b).abs()<=tol*(a.abs()+b.abs()).max(1e-18),"{a} != {b}");}

#[test]
fn g1_oblique_axial_support_matches_spatial_mass_strain_force_and_hessian() {
    for ends in [StringEndMotion::First,StringEndMotion::Second] {
        let (x,v)=(2e-5,0.03);let w=wire(2);let l=w.length_m();
        let support=TranslatingSupport {initial_position_m:x,initial_velocity_m_s:v,..rail()};
        let s=SupportedStrings::with_longitudinal_motion(&[w.clone()],&[law()],&[ends],
            &[geometry([0.6,0.8,0.0])],&[axial(vec![1,2,3,4])],support).unwrap();
        let relative:Vec<_>=(0..8).map(|i|ModalAcousticState {
            displacement_m_sqrt_kg:if i<4 {2e-7*(i as f64-1.2)} else {1e-8*(i as f64-4.7)},
            velocity_m_sqrt_kg_per_s:0.001*(i as f64-2.3),
        }).collect();
        let body=s.body_with_relative_states(&relative).unwrap();
        let z:Vec<_>=body.initial.iter().map(|q|q.displacement_m_sqrt_kg).collect();
        let velocity:Vec<_>=body.initial.iter().map(|q|q.velocity_m_sqrt_kg_per_s).collect();
        let trans=(2.0/(w.linear_density_kg_m*l)).sqrt();let solid=(2.0/(0.01*l)).sqrt();
        let sign=if ends==StringEndMotion::First {-1.0} else {1.0};
        let cells=8192;let dx=l/cells as f64;
        let (mut kinetic,mut energy,mut gas,mut mean)=(0.5*support.mass_kg*v*v,0.5*support.stiffness_n_m*x*x,0.0,0.0);
        for cell in 0..cells {
            let t=(cell as f64+0.5)/cells as f64;let lift=ends.lift(t);
            let (mut slope2,mut curvature2,mut speed2)=(0.0,0.0,0.0);
            for plane in 0..2 {
                let projection=if plane==0 {0.8} else {0.0};
                let (mut slope,mut curvature,mut speed)=(projection*sign*x/l,0.0,projection*lift*v);
                for mode in 0..2 {
                    let k=(mode+1) as f64*core::f64::consts::PI/l;let q=relative[plane*2+mode];
                    slope+=trans*k*q.displacement_m_sqrt_kg*(k*t*l).cos();
                    curvature-=trans*k*k*q.displacement_m_sqrt_kg*(k*t*l).sin();
                    speed+=trans*q.velocity_m_sqrt_kg_per_s*(k*t*l).sin();
                }
                slope2+=slope*slope;curvature2+=curvature*curvature;speed2+=speed*speed;
            }
            let (mut strain,mut axial_speed)=(0.6*sign*x/l,0.6*lift*v);
            for mode in 0..4 {
                let k=(mode+1) as f64*core::f64::consts::PI/l;let q=relative[4+mode];
                strain+=solid*k*q.displacement_m_sqrt_kg*(k*t*l).cos();
                axial_speed+=solid*q.velocity_m_sqrt_kg_per_s*(k*t*l).sin();
            }
            strain+=0.5*slope2;mean+=strain*dx/l;
            kinetic+=0.5*(w.linear_density_kg_m*speed2+0.01*axial_speed*axial_speed)*dx;
            energy+=(0.5*w.tension_n*slope2+0.5*w.bending_stiffness_n_m2*curvature2+50.0*strain*strain)*dx;
            gas+=3.0*w.linear_density_kg_m*speed2*dx;
        }
        near(0.5*velocity.iter().map(|v|v*v).sum::<f64>(),kinetic,2e-8);
        near(s.potential(&z),energy,2e-10);
        near(s.longitudinal_observation(0,&z).unwrap().mean_dynamic_strain,mean,1e-10);
        near(s.total_mass_kg(),support.mass_kg+l*(0.012*0.64+0.01*0.36),1e-14);
        let n=s.mode_count();let mut original=vec![0.0;4*n*n];let mut loaded=original.clone();
        s.add_resistance(0,2*n,&mut original).unwrap();
        s.clone().with_absolute_drag(vec![AbsoluteStringDrag {modal_per_s:vec![3.0;4],remainder_per_s:3.0}])
            .unwrap().add_resistance(0,2*n,&mut loaded).unwrap();
        let power=(0..n).flat_map(|i|(0..n).map(move|j|(i,j))).map(|(i,j)|
            velocity[i]*(loaded[(2*i+1)*2*n+2*j+1]-original[(2*i+1)*2*n+2*j+1])*velocity[j]).sum();
        near(power,gas,2e-8);
        let mut gradient=vec![0.0;n];s.gradient(&z,&mut gradient);
        for i in 0..n {
            let mut a=z.clone();let mut b=z.clone();a[i]+=1e-10;b[i]-=1e-10;
            near(gradient[i],(s.potential(&a)-s.potential(&b))/2e-10,2e-5);
        }
        let d:Vec<_>=(0..n).map(|i|0.001*(i as f64-3.2)).collect();
        let mut actual=vec![0.0;n];s.hessian_vector(&z,&d,&mut actual);
        let plus:Vec<_>=z.iter().zip(&d).map(|(q,d)|q+1e-8*d).collect();
        let minus:Vec<_>=z.iter().zip(&d).map(|(q,d)|q-1e-8*d).collect();
        let (mut gp,mut gm)=(vec![0.0;n],vec![0.0;n]);s.gradient(&plus,&mut gp);s.gradient(&minus,&mut gm);
        for i in 0..n {near(actual[i],(gp[i]-gm[i])/2e-8,2e-6);}
        for station in [0.0,0.37,1.0] {
            let row=s.absolute_axial_shape(0,station).unwrap();
            let sine=if station==0.0 || station==1.0 {0.0} else {(1..=4).map(|mode|
                solid*relative[3+mode].displacement_m_sqrt_kg*(mode as f64*core::f64::consts::PI*station).sin()).sum()};
            let expected=0.6*ends.lift(station)*x+sine;
            near(row.iter().zip(&z).map(|(b,q)|b*q).sum(),expected,1e-12);
        }
    }
}

#[test]
fn g3_axial_waves_move_the_common_support_and_unplayed_span_with_closed_work() {
    let s=SupportedStrings::with_longitudinal_motion(&[wire(1),wire(1)],&[law(),law()],
        &[StringEndMotion::Second,StringEndMotion::First],&[geometry([1.0,0.0,0.0]);2],
        &[axial(vec![1,2]),axial(vec![1,2])],rail()).unwrap();
    let mut initial=vec![ModalAcousticState::default();8];initial[2].displacement_m_sqrt_kg=1e-7;
    let row=s.absolute_axial_shape(0,0.37).unwrap();let n=s.mode_count();
    let config=ImpactConfig {dt_s:2e-6,max_steps:256,maximum_energy_j:1.0,
        energy_absolute_tolerance_j:1e-12,energy_relative_tolerance:1e-9,maximum_generalized_force:1e4};
    let mut run=ImpactSystem::new(vec![s.body_with_relative_states(&initial).unwrap()],vec![],vec![],vec![],config)
        .unwrap().prepare_analytic().unwrap();
    let before=run.stored_energy_j();let (mut work,mut loss,mut bridge,mut neighbor)=(0.0,0.0,0.0_f64,0.0_f64);
    for step in 0..256 {
        let force=if step<16 {0.001} else {0.0};
        let old:Vec<_>=run.state().chunks_exact(2).map(|q|q[0]).collect();
        let frame=run.step(&row.iter().map(|b|b*force).collect::<Vec<_>>(),&CancelGate::new_clock_free()).unwrap();
        let new:Vec<_>=run.state().chunks_exact(2).map(|q|q[0]).collect();
        let port_work=force*row.iter().zip(old.iter().zip(&new)).map(|(b,(a,z))|b*(z-a)).sum::<f64>();
        near(frame.supplied_work_j,port_work,1e-10);
        work+=frame.supplied_work_j;loss+=frame.dissipated_energy_j;
        assert!((frame.stored_energy_j+loss-before-work).abs()<1e-10);
        let x=new[n-1]*s.force_weight();bridge=bridge.max(x.abs());
        neighbor=neighbor.max((new[6]-s.0.translation[6]*x).abs());
        assert_eq!([new[0],new[1],new[4],new[5]],[0.0;4]);
    }
    assert!(bridge>1e-12 && neighbor>1e-14,"axial support/neighbor path must be live");
    assert!(StringSupportGeometry::try_new([[0.0;3],[0.2,0.0,0.0]],[1.0,0.0,0.0],[1.0,0.0,0.0]).is_err());
    assert!(SupportedStrings::with_transverse_motion(&[wire(1)],&[law()],&[StringEndMotion::Second],
        &[StringTransverseMotion::Spatial {geometry:geometry([1.0,0.0,0.0])}],rail()).is_err());
}
