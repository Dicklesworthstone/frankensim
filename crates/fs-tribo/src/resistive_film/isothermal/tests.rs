use super::*;
use super::super::{FilmCell,FilmChannel,FilmLimits,GapPort};
fn gas()->AmbientGas {AmbientGas{pressure_pa:100000.0,temperature_k:293.15,specific_gas_constant_j_kg_k:287.05}}
fn gap(h:f64,b:&[f64])->GapPort {GapPort{reference_m:h,closure:b.to_vec()}}
fn model(drained:bool)->IsothermalFilm {
    let a=gap(0.001,&[1.0,-1.0]);let b=gap(0.0012,&[0.7,-0.7]);
    let geometry=ResistiveFilm::new(vec![FilmCell{area_m2:0.002,gap:a.clone()},
        FilmCell{area_m2:0.003,gap:b}],vec![
        FilmChannel{from:0,to:Some(1),width_m:0.04,length_m:0.02,gap:gap(0.0011,&[0.85,-0.85])},
        FilmChannel{from:0,to:None,width_m:0.03,length_m:0.01,
            gap:if drained{a}else{gap(0.0,&[0.0;2])}}],2,1.8e-5,
        FilmLimits{minimum_cell_gap_m:1e-6,maximum_gap_m:0.01,maximum_pressure_pa:300000.0}).unwrap();
    IsothermalFilm::at_ambient(geometry,&[0.0;2],gas()).unwrap()
}
fn close(a:f64,b:f64,tol:f64) {assert!((a-b).abs()<=tol*(a.abs()+b.abs()).max(1e-8),"{a} != {b}");}

#[test]
fn sealed_pocket_has_positive_mass_and_reciprocal_compression_storage() {
    let m=model(false);let z=m.initial_state();let initial=m.observe(&[0.0;2],z).unwrap();
    assert_eq!(initial.free_energy_j,0.0);assert_eq!(initial.minimum_pressure_pa,gas().pressure_pa);
    let q=[0.0002,0.0];let state=m.state(&q,z).unwrap();
    close(state.pressure[0],125000.0,1e-13);close(state.observation.mass_kg,initial.mass_kg,1e-14);
    assert!(state.observation.free_energy_j>0.0);
    let mut gq=[0.0;2];let mut gz=[0.0;2];m.gradient_into(&q,z,&mut gq,&mut gz).unwrap();
    assert!(gq[0]>0.0);assert_eq!(gq[0],-gq[1]);
    let translated=m.observe(&[0.03,0.03],z).unwrap();assert_eq!(translated.free_energy_j,0.0);
}
#[test]
fn physical_gradient_recovers_compressible_poiseuille_mass_flux() {
    let m=model(false);let mut z=m.initial_state().to_vec();z[0]*=1.1;z[1]*=0.9;
    let mut gq=[0.0;2];let mut gz=[0.0;2];let q=[0.0;2];let mut flow=[0.0;2];
    m.gradient_into(&q,&z,&mut gq,&mut gz).unwrap();let power=m.flow_into(&q,&z,&gz,&mut flow).unwrap();
    let g=m.geometry.factors[0]*0.0011_f64.powi(3);
    let rt=gas().temperature_k*gas().specific_gas_constant_j_kg_k;
    let expected=g*(110000.0_f64.powi(2)-90000.0_f64.powi(2))/(2.0*rt);
    close(flow[0]*m.scales[0]/rt,expected,1e-12);
    close(flow[1]*m.scales[1]/rt,-expected,1e-12);
    close(dot(&gz,&flow),power,1e-12);assert!(power>0.0);
}
#[test]
fn storage_gradient_and_hessian_agree_with_independent_differences() {
    let m=model(true);let q=[0.0001,-0.00003];let mut z=m.initial_state().to_vec();z[1]*=0.95;
    let mut gq=[0.0;2];let mut gz=[0.0;2];m.gradient_into(&q,&z,&mut gq,&mut gz).unwrap();
    for i in 0..4 {
        let mut qa=q;let mut qb=q;let mut za=z.clone();let mut zb=z.clone();
        let h=if i<2{1e-9}else{1e-6};
        if i<2{qa[i]+=h;qb[i]-=h;}else{za[i-2]+=h;zb[i-2]-=h;}
        let fd=(m.observe(&qa,&za).unwrap().free_energy_j-m.observe(&qb,&zb).unwrap().free_energy_j)/(2.0*h);
        close(fd,if i<2{gq[i]}else{gz[i-2]},1e-6);
    }
    let dq=[0.0001,-0.0002];let dz=[0.03,-0.04];let mut hq=[0.0;2];let mut hz=[0.0;2];
    m.hessian_into(&q,&z,&dq,&dz,&mut hq,&mut hz).unwrap();
    let eps=1e-5;let mut qa=q;let mut qb=q;let mut za=z.clone();let mut zb=z.clone();
    for i in 0..2{qa[i]+=eps*dq[i];qb[i]-=eps*dq[i];za[i]+=eps*dz[i];zb[i]-=eps*dz[i];}
    let (mut ga,mut gb,mut ma,mut mb)=([0.0;2],[0.0;2],[0.0;2],[0.0;2]);
    m.gradient_into(&qa,&za,&mut ga,&mut ma).unwrap();m.gradient_into(&qb,&zb,&mut gb,&mut mb).unwrap();
    for i in 0..2{close(hq[i],(ga[i]-gb[i])/(2.0*eps),1e-7);close(hz[i],(ma[i]-mb[i])/(2.0*eps),1e-7);}
    assert!(dot(&dq,&hq)+dot(&dz,&hz)>=0.0);
}
#[test]
fn arbitrary_discrete_effort_is_passive_and_sealed_exchange_conserves_mass() {
    for drained in [false,true] {
        let m=model(drained);let z=m.initial_state();let mut out=[0.0;2];
        for e in [[0.01,-0.02],[0.0,0.0],[-0.03,-0.07]] {
            let power=m.flow_into(&[0.0;2],z,&e,&mut out).unwrap();
            close(dot(&e,&out),power,1e-12);assert!(power>=0.0);
            if !drained{close(m.scales[0]*out[0],-m.scales[1]*out[1],1e-12);}
        }
    }
}
#[test]
fn flow_tangent_includes_density_gap_and_discrete_effort_near_equal_pressure() {
    let m=model(true);
    for q in [[0.0;2],[0.0001,-0.00003]] {
        let z=m.initial_state();let e=[0.01,-0.02];let dq=[0.0001,0.00004];let dz=[0.02,-0.04];let de=[0.03,0.02];
        let mut tangent=[0.0;2];m.flow_tangent_into(&q,z,&e,&dq,&dz,&de,&mut tangent).unwrap();
        let eps=1e-5;let mut qa=q;let mut qb=q;let mut za=z.to_vec();let mut zb=z.to_vec();let mut ea=e;let mut eb=e;
        for i in 0..2{qa[i]+=eps*dq[i];qb[i]-=eps*dq[i];za[i]+=eps*dz[i];zb[i]-=eps*dz[i];ea[i]+=eps*de[i];eb[i]-=eps*de[i];}
        let mut a=[0.0;2];let mut b=[0.0;2];m.flow_into(&qa,&za,&ea,&mut a).unwrap();m.flow_into(&qb,&zb,&eb,&mut b).unwrap();
        for i in 0..2{close(tangent[i],(a[i]-b[i])/(2.0*eps),2e-7);}
    }
}
#[test]
fn channels_can_close_and_reopen_without_discarding_trapped_mass() {
    let mut m=model(true);m.geometry.channels[0].gap=gap(0.0,&[0.0;2]);
    m.geometry.channels[1].gap=gap(0.0001,&[1.0,-1.0]);
    let z=m.initial_state();let mut flow=[7.0;2];
    m.flow_into(&[0.0002,0.0],z,&[0.02;2],&mut flow).unwrap();assert_eq!(flow,[0.0;2]);
    let closed=m.observe(&[0.0002,0.0],z).unwrap();assert!(closed.maximum_pressure_pa>100000.0);
    m.flow_into(&[0.0;2],z,&[0.02;2],&mut flow).unwrap();assert!(flow[0]>0.0);assert_eq!(flow[1],0.0);
    assert_eq!(m.observe(&[0.0;2],z).unwrap().mass_kg,closed.mass_kg);
}
#[test]
fn invalid_mass_pressure_and_collapsed_gap_leave_outputs_untouched() {
    let m=model(true);let z=m.initial_state();let mut qout=[21.0;2];let mut zout=[22.0;2];
    for (q,z) in [([0.002,0.0],z.to_vec()),([0.0;2],vec![0.0,z[1]]),
        ([0.0;2],vec![f64::NAN,z[1]]),([0.0;2],vec![10.0*z[0],z[1]])] {
        assert!(m.gradient_into(&q,&z,&mut qout,&mut zout).is_err());
        assert_eq!(qout,[21.0;2]);assert_eq!(zout,[22.0;2]);
    }
    assert!(m.flow_into(&[0.0;2],z,&[f64::NAN;2],&mut zout).is_err());assert_eq!(zout,[22.0;2]);
    assert!(IsothermalFilm::at_ambient(m.geometry.clone(),&[0.0;2],AmbientGas{pressure_pa:0.0,..gas()}).is_err());
}
