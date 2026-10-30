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

#[test]
fn g1_axial_viscous_power_is_the_spatial_strain_rate_integral_with_endpoint_work() {
    for ends in [StringEndMotion::First,StringEndMotion::Second,StringEndMotion::Both] {
        let w=wire(2);let l=w.length_m();let viscosity=0.06;let support_speed=0.03;
        let base=SupportedStrings::with_longitudinal_motion(&[w.clone()],&[law()],&[ends],
            &[geometry([0.6,0.8,0.0])],&[axial(vec![1,3])],
            TranslatingSupport {initial_velocity_m_s:support_speed,..rail()}).unwrap();
        let loaded=base.clone().with_axial_viscosity(vec![viscosity]).unwrap();
        assert_eq!(loaded.axial_viscosity_n_s(),Some(&[viscosity][..]));
        assert_eq!(loaded.total_mass_kg(),base.total_mass_kg());
        assert_eq!(loaded.clock_bounds(),base.clock_bounds());
        let n=base.mode_count();let dim=2*n;
        let (mut original,mut damped)=(vec![0.0;dim*dim],vec![0.0;dim*dim]);
        base.add_resistance(0,dim,&mut original).unwrap();
        loaded.add_resistance(0,dim,&mut damped).unwrap();
        for i in 0..dim {for j in 0..dim {assert_eq!(damped[i*dim+j],damped[j*dim+i]);}}
        for moving_internal in [false,true] {
            let mut relative=vec![ModalAcousticState::default();n-1];
            if moving_internal {
                // Arbitrary nonzero transverse velocity must not be reused as
                // small-strain axial viscosity. Solid and gas-loaded masses differ.
                for (i,state) in relative.iter_mut().enumerate() {
                    state.velocity_m_sqrt_kg_per_s=0.001*(i as f64-2.3);
                }
            }
            let body=loaded.body_with_relative_states(&relative).unwrap();
            let velocity:Vec<_>=body.initial.iter().map(|q|q.velocity_m_sqrt_kg_per_s).collect();
            let power: f64=(0..n).map(|i|(0..n).map(|j|velocity[i]
                *(damped[(2*i+1)*dim+2*j+1]-original[(2*i+1)*dim+2*j+1])*velocity[j]).sum::<f64>()).sum();
            let sign=match ends {StringEndMotion::First=>-1.0,StringEndMotion::Second=>1.0,StringEndMotion::Both=>0.0};
            let solid_normalization=(2.0/(0.01*l)).sqrt();
            let cells=4096;let mut spatial=0.0;
            for cell in 0..cells {
                let x=l*(cell as f64+0.5)/cells as f64;
                let mut strain_rate=sign*0.6*support_speed/l;
                for (local,mode) in [1,3].into_iter().enumerate() {
                    let k=mode as f64*core::f64::consts::PI/l;
                    strain_rate+=solid_normalization*k*relative[4+local].velocity_m_sqrt_kg_per_s*(k*x).cos();
                }
                spatial+=viscosity*strain_rate*strain_rate*l/cells as f64;
            }
            let scale=spatial.max(viscosity*support_speed*support_speed/l);
            assert!((power-spatial).abs()<2e-12*scale,"{power} != {spatial}");
            if !moving_internal {
                near(spatial,viscosity/l*(sign*0.6*support_speed).powi(2),2e-13);
            }
        }
        let zero=base.clone().with_axial_viscosity(vec![0.0]).unwrap();
        let mut zero_resistance=vec![0.0;dim*dim];
        zero.add_resistance(0,dim,&mut zero_resistance).unwrap();
        assert_eq!(original,zero_resistance,"explicit zero preserves the elastic operator");
        assert!(loaded.with_axial_viscosity(vec![viscosity]).is_err());
        for bad in [vec![],vec![-1.0],vec![f64::NAN],vec![f64::MAX]] {
            assert!(base.clone().with_axial_viscosity(bad).is_err());
        }
    }
    let transverse=SupportedStrings::new(&[wire(1)],law(),rail()).unwrap();
    assert!(transverse.with_axial_viscosity(vec![1.0]).is_err());
}

#[test]
fn g1_finite_strain_viscosity_matches_spatial_power_and_the_full_support_tangent() {
    for ends in [StringEndMotion::First,StringEndMotion::Second,StringEndMotion::Both] {
        let w=wire(2);let l=w.length_m();let rigidity=0.06;
        let support=TranslatingSupport {initial_position_m:2e-4,initial_velocity_m_s:0.03,..rail()};
        let base=SupportedStrings::with_longitudinal_motion(&[w.clone()],&[law()],&[ends],
            &[geometry([0.6,0.48,0.64])],&[axial(vec![1,2,3,4])],support).unwrap();
        let loaded=base.clone().with_finite_strain_axial_viscosity(vec![rigidity]).unwrap();
        let n=loaded.mode_count();let dim=2*n;
        let (mut before,mut after)=(vec![0.0;dim*dim],vec![0.0;dim*dim]);
        base.add_resistance(0,dim,&mut before).unwrap();
        loaded.add_resistance(0,dim,&mut after).unwrap();
        assert_eq!(before,after,"finite-strain viscosity must not also be a linear resistance");
        let relative:Vec<_>=(0..8).map(|i|ModalAcousticState {
            displacement_m_sqrt_kg:if i<4 {2e-6*(i as f64-1.2)} else {1e-8*(i as f64-4.7)},
            velocity_m_sqrt_kg_per_s:if i<4 {0.001*(i as f64-2.3)} else {1e-6*(i as f64-6.2)},
        }).collect();
        let body=loaded.body_with_relative_states(&relative).unwrap();
        let z:Vec<_>=body.initial.iter().map(|s|s.displacement_m_sqrt_kg).collect();
        let v:Vec<_>=body.initial.iter().map(|s|s.velocity_m_sqrt_kg_per_s).collect();
        let mut force=vec![0.0;n];assert!(loaded.viscous_reaction(&z,&v,&mut force));
        let power=-force.iter().zip(&v).map(|(f,v)|f*v).sum::<f64>();
        let sign=match ends {StringEndMotion::First=>-1.0,StringEndMotion::Second=>1.0,StringEndMotion::Both=>0.0};
        let transverse=(2.0/(w.linear_density_kg_m*l)).sqrt();let solid=(2.0/(0.01*l)).sqrt();
        let cells=512;let mut spatial=0.0;
        for cell in 0..cells {
            let x=l*(cell as f64+0.5)/cells as f64;
            let mut rate=sign*0.6*support.initial_velocity_m_s/l;
            for mode in 0..4 {
                let k=(mode+1) as f64*core::f64::consts::PI/l;
                rate+=solid*k*relative[4+mode].velocity_m_sqrt_kg_per_s*(k*x).cos();
            }
            for (plane,projection) in [0.48,0.64].into_iter().enumerate() {
                let mut slope=sign*projection*support.initial_position_m/l;
                let mut slope_rate=sign*projection*support.initial_velocity_m_s/l;
                for mode in 0..2 {
                    let k=(mode+1) as f64*core::f64::consts::PI/l;let state=relative[2*plane+mode];
                    slope+=transverse*k*state.displacement_m_sqrt_kg*(k*x).cos();
                    slope_rate+=transverse*k*state.velocity_m_sqrt_kg_per_s*(k*x).cos();
                }
                rate+=slope*slope_rate;
            }
            spatial+=rigidity*rate*rate*l/cells as f64;
        }
        assert!(power>0.0);near(power,spatial,2e-12);
        // Both product-rule terms and the kinetic pullback depend on the
        // support component. Probe the complete normalized coordinate direction.
        let dz:Vec<_>=(0..n).map(|i|0.0001*(i as f64-3.1)).collect();
        let dv:Vec<_>=(0..n).map(|i|0.02*(i as f64+1.7)).collect();
        let mut tangent=vec![0.0;n];assert!(loaded.viscous_reaction_tangent(&z,&v,&dz,&dv,&mut tangent));
        let h=1e-7;
        let plus:Vec<_>=z.iter().zip(&dz).map(|(z,d)|z+h*d).collect();
        let minus:Vec<_>=z.iter().zip(&dz).map(|(z,d)|z-h*d).collect();
        let vplus:Vec<_>=v.iter().zip(&dv).map(|(v,d)|v+h*d).collect();
        let vminus:Vec<_>=v.iter().zip(&dv).map(|(v,d)|v-h*d).collect();
        let (mut fp,mut fm)=(vec![0.0;n],vec![0.0;n]);
        assert!(loaded.viscous_reaction(&plus,&vplus,&mut fp));
        assert!(loaded.viscous_reaction(&minus,&vminus,&mut fm));
        for i in 0..n {near(tangent[i],(fp[i]-fm[i])/(2.0*h),2e-6);}
        let zero=base.clone().with_finite_strain_axial_viscosity(vec![0.0]).unwrap();
        assert!(!zero.has_finite_strain_axial_viscosity());
        assert!(zero.viscous_reaction(&z,&v,&mut force));assert_eq!(force,vec![0.0;n]);
        if ends==StringEndMotion::Both {
            let mut translation=vec![0.0;n];translation[n-1]=1.0;
            for (value,c) in translation.iter_mut().zip(&base.0.translation) {*value=c*base.force_weight();}
            assert!(loaded.viscous_reaction(&z,&translation,&mut force));
            assert_eq!(force,vec![0.0;n],"rigid endpoint/string translation has no axial strain rate");
        }
        assert!(loaded.clone().with_axial_viscosity(vec![rigidity]).is_err());
        assert!(loaded.with_finite_strain_axial_viscosity(vec![rigidity]).is_err());
    }
}

#[test]
fn g1_transverse_motion_dissipates_axially_at_finite_slope_with_no_axial_velocity() {
    let base=SupportedStrings::with_longitudinal_motion(&[wire(1)],&[law()],&[StringEndMotion::Second],
        &[geometry([0.6,0.8,0.0])],&[axial(vec![1,2])],rail()).unwrap();
    let loaded=base.clone().with_finite_strain_axial_viscosity(vec![0.06]).unwrap();
    let relative=vec![ModalAcousticState {displacement_m_sqrt_kg:1e-5,
        velocity_m_sqrt_kg_per_s:0.001},ModalAcousticState::default(),
        ModalAcousticState::default(),ModalAcousticState::default()];
    let body=loaded.body_with_relative_states(&relative).unwrap();
    let z:Vec<_>=body.initial.iter().map(|s|s.displacement_m_sqrt_kg).collect();
    let v:Vec<_>=body.initial.iter().map(|s|s.velocity_m_sqrt_kg_per_s).collect();
    let mut force=vec![0.0;loaded.mode_count()];assert!(loaded.viscous_reaction(&z,&v,&mut force));
    assert!(force[0]*v[0]<0.0,"geometric strain rate must reach the transverse motion");
    assert!(force[3].abs()>0.0,"the same viscous stress reacts on the second axial harmonic");
    let w=wire(1);let scale=(2.0/(w.linear_density_kg_m*w.length_m())).sqrt()
        *core::f64::consts::PI/w.length_m();
    let a=scale*relative[0].displacement_m_sqrt_kg;
    let b=scale*relative[0].velocity_m_sqrt_kg_per_s;
    // epsilon_dot=a*b*cos^2(pi*x/L), integral cos^4 = 3L/8.
    near(-force.iter().zip(&v).map(|(f,v)|f*v).sum::<f64>(),
        0.06*wire(1).length_m()*3.0/8.0*(a*b).powi(2),1e-13);
}
