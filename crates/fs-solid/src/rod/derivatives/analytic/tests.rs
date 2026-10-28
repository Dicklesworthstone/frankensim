use super::*;
use crate::rod::{RodSection, TipLoad};
use crate::rod::derivatives::{RodEnergy, UniformEnergy, perturb};
use std::cell::Cell;

fn specimen(angle:f64)->Rod {
    let mut r=Rod::straight(0.5,1,RodSection {ea:40.0,ga:20.0,gj:3.0,ei:2.0});
    r.positions[1]=[0.57,0.04,-0.03];
    r.frames[0]=So3::exp(So3Tangent::new(Vec3::new(0.2,-0.1,0.3))).unwrap();
    r.frames[1]=r.frames[0].body_plus(So3Tangent::new(Vec3::new(0.36*angle,-0.48*angle,0.8*angle))).unwrap();
    r
}
fn response(r:&Rod,tangent:bool)->PairResponse {
    let kin=Kinematics::new(r,0).unwrap();
    let (g,k)=kin.strains();let s=r.section;
    kin.linearize(Vec3::new(s.ea*g[0],s.ga*g[1],s.ga*g[2]),
        Vec3::new(s.gj*k[0],s.ei*k[1],s.ei*k[2]),
        [s.ea,s.ga,s.ga],[s.gj,s.ei,s.ei],tangent).unwrap()
}

#[test]
fn g1_analytical_force_and_body_tangent_match_independent_probes() {
    for angle in [0.0,1e-7,0.02,0.4,1.5,3.0] {
        let rod=specimen(angle);
        let exact=response(&rod,true);
        let matrix=exact.tangent.unwrap();
        for col in 0..12 {
            let step=2e-6;
            let mut plus=rod.clone();let mut minus=rod.clone();
            perturb(&mut plus,col/6,col%6,step).unwrap();
            perturb(&mut minus,col/6,col%6,-step).unwrap();
            let energy_derivative=(plus.energy().unwrap()-minus.energy().unwrap())/(2.0*step);
            let expected=exact.gradient[col/6][col%6];
            assert!((energy_derivative-expected).abs()<3e-6*expected.abs().max(1.0));
            let upper=response(&plus,false);let lower=response(&minus,false);
            for row in 0..12 {
                let probe=(upper.gradient[row/6][row%6]-lower.gradient[row/6][row%6])/(2.0*step);
                assert!((matrix[row][col]-probe).abs()<5e-6*matrix[row][col].abs().max(1.0),
                    "angle {angle}, ({row},{col}): {} != {probe}",matrix[row][col]);
            }
        }
    }
}

#[test]
fn g0_body_torque_basis_terms_are_retained_not_symmetrized() {
    let rod=specimen(1.2);
    let exact=response(&rod,true);let k=exact.tangent.unwrap();
    for row in 0..12 {for col in 0..12 {
        let expected=if row/6==col/6 && row%6>=3 && col%6>=3 {
            let force=exact.gradient[row/6];
            hat(Vec3::new(force[3],force[4],force[5])).m[3*(row%6-3)+col%6-3]
        } else {0.0};
        assert!((k[row][col]-k[col][row]-expected).abs()<2e-11);
    }}
    let f0=vector(exact.gradient[0][..3].try_into().unwrap());
    let f1=vector(exact.gradient[1][..3].try_into().unwrap());
    assert!(vadd(f0,f1).norm()<1e-12);
    let m0=rod.frames[0].rotate(vector(exact.gradient[0][3..].try_into().unwrap())).unwrap();
    let m1=rod.frames[1].rotate(vector(exact.gradient[1][3..].try_into().unwrap())).unwrap();
    let total=vadd(vadd(vector(rod.positions[0]).cross(f0),vector(rod.positions[1]).cross(f1)),vadd(m0,m1));
    assert!(total.norm()<2e-11,"{total:?}");
}

#[test]
fn g0_reference_tangent_has_six_rigid_modes_and_exact_axial_stiffness() {
    let rod=Rod::straight(0.5,1,RodSection {ea:40.0,ga:20.0,gj:3.0,ei:2.0});
    let exact=response(&rod,true);let k=exact.tangent.unwrap();
    assert!(exact.gradient.iter().flatten().all(|v| *v==0.0));
    assert_eq!(k[0][0],80.0);assert_eq!(k[0][6],-80.0);
    for mode in 0..6 {
        let mut d=[0.0;12];
        if mode<3 {d[mode]=1.0;d[6+mode]=1.0;}
        else {
            let mut axis=[0.0;3];axis[mode-3]=1.0;
            for node in 0..2 {
                let translation=vector(axis).cross(vector(rod.positions[node]));
                d[6*node..6*node+3].copy_from_slice(&array(translation));
                d[6*node+3..6*node+6].copy_from_slice(&axis);
            }
        }
        for row in &k {assert!(row.iter().zip(d).map(|(a,b)|a*b).sum::<f64>().abs()<1e-12);}
    }
}

#[test]
fn g3_translated_world_coordinates_do_not_need_representable_fd_steps() {
    let mut rod=Rod::straight(1.0,4,RodSection {ea:40.0,ga:20.0,gj:3.0,ei:2.0});
    let before=UniformEnergy.linearize(&rod,&(0..24).collect::<Vec<_>>(),true,&||Ok(())).unwrap().unwrap();
    for p in &mut rod.positions {for x in p {*x+=2.0_f64.powi(40);}}
    let after=UniformEnergy.linearize(&rod,&(0..24).collect::<Vec<_>>(),true,&||Ok(())).unwrap().unwrap();
    assert_eq!(before.gradient,after.gradient);assert_eq!(before.tangent,after.tangent);
    // This is precisely a configuration where the old translation probe is lost.
    assert_eq!(rod.positions[1][0]+3e-6,rod.positions[1][0]);
}

#[test]
fn g1_original_newton_consumes_analytical_tangent_without_probe_fallback() {
    struct ExactOnly<'a>(&'a Cell<usize>);
    impl RodEnergy<SolidError> for ExactOnly<'_> {
        fn value(&self,r:&Rod)->Result<f64,SolidError>{r.energy()}
        fn linearize(&self,r:&Rod,f:&[usize],t:bool,c:&impl Fn()->Result<(),SolidError>)
            ->Result<Option<Linearization>,SolidError> {
            if t {self.0.set(self.0.get()+1);}
            UniformEnergy.linearize(r,f,t,c)
        }
        fn derivative(&self,_:&Rod,_:&Rod,_:&Rod,_:usize,_:usize,_:f64)->Result<f64,SolidError>{
            panic!("native Newton must not use finite-difference force probes")
        }
    }
    let mut rod=Rod::straight(1.0,4,RodSection {ea:100.0,ga:50.0,gj:3.0,ei:2.0});
    let count=Cell::new(0);
    rod.solve_static_with_energy(&TipLoad {force:[2.0,0.0,0.0],moment:[0.0;3]},
        1,1e-10,&ExactOnly(&count),&||Ok(()),None,1.0).unwrap();
    assert!(count.get()>0);
    for (node,p) in rod.positions.iter().enumerate(){assert!((p[0]-1.02*node as f64/4.0).abs()<1e-12);}
}

#[test]
fn g4_analytical_assembly_cancellation_and_log_cut_refuse_without_mutation() {
    let mut rod=specimen(0.2);let before=rod.positions.clone();
    let gate=||Err(SolidError::InvalidInput {what:"cancelled fixture".into()});
    assert!(UniformEnergy.linearize(&rod,&[],true,&gate).is_err());
    assert_eq!(before,rod.positions);
    rod.frames[1]=rod.frames[0].body_plus(So3Tangent::new(Vec3::new(core::f64::consts::PI,0.0,0.0))).unwrap();
    assert!(rod.energy().is_ok());
    assert!(UniformEnergy.linearize(&rod,&[],false,&||Ok(())).is_err());
    assert_eq!(before,rod.positions);
}
