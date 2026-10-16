use super::*;

fn diagonal(value:f64)->Tangent6 {
    core::array::from_fn(|i|core::array::from_fn(|j|if i==j {value}else{0.0}))
}
fn close(a:f64,b:f64,tol:f64) {assert!((a-b).abs()<=tol*(1.0+a.abs().max(b.abs())),"{a} != {b}");}
fn cclose(a:C64,b:C64,tol:f64) {close(a.re,b.re,tol);close(a.im,b.im,tol);}
fn model(scale:f64)->TensorProny {
    let mut e=diagonal(3.0);e[0][2]=0.7;e[1][3]=-0.2;e[4][1]=0.3;
    let mut a=diagonal(2.0);a[0][2]=1.3;a[2][4]=-0.7;
    let mut b=diagonal(1.0);b[1][2]=-1.5;b[5][0]=0.5;
    let scaled=|mut m:Tangent6| {for row in &mut m {for v in row {*v*=det::sqrt(scale);}}m};
    TensorProny::new_mandel(gram(&scaled(e)).unwrap(),vec![
        TensorPronyBranch::from_gram_rows(scaled(a),0.01).unwrap(),
        TensorPronyBranch::from_gram_rows(scaled(b),0.03).unwrap()],0.1).unwrap()
}

#[test]
fn isotropic_plane_stress_and_section_recover_the_closed_form() {
    let young=12e9;let nu=0.25;let tau=0.02;
    let c=crate::IsotropicElastic::new(young,nu,0.1).unwrap().stiffness();
    let branch=TensorPronyBranch::from_stiffness_mandel(c.map(|r|r.map(|v|1.5*v)),tau).unwrap();
    let law=TensorProny::from_tensor_shear(c,vec![branch],0.1).unwrap();
    for omega in [0.0,1.0,50.0,1e5] {
        let z=omega*tau;let em=C64::new(young*(1.0+1.5*z*z/(1.0+z*z)),young*1.5*z/(1.0+z*z));
        let ps=law.plane_stress(omega).unwrap();let d=ps.engineering();
        let expected=[[1.0,nu,0.0],[nu,1.0,0.0],[0.0,0.0,0.5*(1.0-nu)]];
        for i in 0..3 {for j in 0..3 {cclose(d[i][j],em.scale(expected[i][j]/(1.0-nu*nu)),2e-13);}}
        let input=[C64::new(0.002,0.001),C64::new(-0.001,0.0002),C64::new(0.0005,-0.0003)];
        let strain=ps.strain_tensor(input).unwrap();
        cclose(strain[2],(input[0]+input[1]).scale(-nu/(1.0-nu)),2e-13);
        cclose(strain[3],input[2].scale(0.5),2e-13);assert_eq!(strain[4],C64::ZERO);assert_eq!(strain[5],C64::ZERO);
        let section=ps.section(0.003).unwrap();let twice=ps.section(0.006).unwrap();
        for i in 0..3 {for j in 0..3 {
            cclose(section.membrane[i][j],d[i][j].scale(0.003),2e-14);
            cclose(section.bending[i][j],d[i][j].scale(0.003_f64.powi(3)/12.0),2e-14);
            cclose(twice.membrane[i][j],section.membrane[i][j].scale(2.0),2e-14);
            cclose(twice.bending[i][j],section.bending[i][j].scale(8.0),2e-14);
        }}
        for h in [0.0,-1.0,f64::NAN,f64::INFINITY,1e308] {assert!(ps.section(h).is_err());}
    }
}

#[test]
fn full_complex_elimination_closes_traction_and_dissipative_work() {
    let law=model(1.0);
    for omega in [0.0,1.0,50.0,1e4] {
        let (storage,loss)=law.modulus_mandel(omega).unwrap();let ps=law.plane_stress(omega).unwrap();
        let input=[C64::new(0.002,0.0007),C64::new(-0.001,0.0003),C64::new(0.0005,-0.0001)];
        let physical=ps.strain_tensor(input).unwrap();
        let e:[C64;6]=core::array::from_fn(|i|physical[i].scale(WEIGHTS[i]));
        let stress:[C64;6]=core::array::from_fn(|i|(0..6).fold(C64::ZERO,|sum,j|
            sum+C64::new(storage[i][j],loss[i][j])*e[j]));
        let reduced=complex_apply(&ps.engineering(),&input);
        for i in 0..3 {
            cclose(stress[OUT_OF_PLANE[i]],C64::ZERO,2e-13);
            cclose(stress[IN_PLANE[i]].scale(1.0/WEIGHTS[IN_PLANE[i]]),reduced[i],2e-13);
        }
        let work=(0..6).fold(C64::ZERO,|sum,i|sum+e[i].conj()*stress[i]);
        let planar=(0..3).fold(C64::ZERO,|sum,i|sum+input[i].conj()*reduced[i]);
        cclose(work,planar,2e-13);assert!(work.im>=-1e-15);
        if omega>0.0 {assert!(work.im>0.0);}
        let d=ps.engineering();for i in 0..3 {for j in 0..3 {cclose(d[i][j],d[j][i],2e-13);}}
    }
    // Condensing individual branches gives the WRONG operator even though
    // each tensor is SPD. This witness detects that tempting implementation.
    let omega=50.0;let actual=law.plane_stress(omega).unwrap().engineering();
    let mut naive=TensorProny::new_mandel(*law.equilibrium_mandel(),vec![],0.1)
        .unwrap().plane_stress(0.0).unwrap().engineering();
    for branch in law.branches() {
        let projected=TensorProny::new_mandel(*branch.stiffness_mandel(),vec![],0.1)
            .unwrap().plane_stress(0.0).unwrap().engineering();
        let z=omega*branch.tau_s;let weight=C64::new(z*z/(1.0+z*z),z/(1.0+z*z));
        for i in 0..3 {for j in 0..3 {naive[i][j]=naive[i][j]+projected[i][j]*weight;}}
    }
    assert!((actual[0][0]-naive[0][0]).abs()>1e-3);
}

#[test]
fn rotation_preserves_stress_work_heat_and_full_frequency_response() {
    let law=model(1.0);let q=crate::tensor::rotation([0.0,0.0,1.0],0.7);
    let rotated=law.rotated(q).unwrap();let mut a=law.history();let mut b=rotated.history();
    for k in 0..20 {
        let e=core::array::from_fn(|i|0.002*det::sin((k*6+i) as f64));
        let actual=a.advance(e,0.007,||false).unwrap();
        let changed=b.advance(crate::tensor::rotate(&e,&q),0.007,||false).unwrap();
        let expected=crate::tensor::rotate(&actual.stress,&q);
        for i in 0..6 {close(changed.stress[i],expected[i],2e-12);}
        close(actual.work,changed.work,2e-12);close(actual.stored_energy,changed.stored_energy,2e-12);
        close(actual.dissipated,changed.dissipated,2e-12);
    }
    let e=[0.002,-0.001,0.0007,0.003,-0.002,0.001];let er=crate::tensor::rotate(&e,&q);
    for omega in [0.0,20.0,2000.0] {
        let (a,b)=law.modulus_mandel(omega).unwrap();let (ar,br)=rotated.modulus_mandel(omega).unwrap();
        for (c,cr) in [(a,ar),(b,br)] {
            let expected=crate::tensor::rotate(&tensor(&apply(&c,&mandel(&e))),&q);
            let actual=tensor(&apply(&cr,&mandel(&er)));
            for i in 0..6 {close(actual[i],expected[i],2e-12);}
        }
    }
    for bad in [[[0.0;3];3],[[-1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]],[[f64::NAN;3];3]] {
        assert!(law.rotated(bad).is_err());
    }
}

#[test]
fn plane_stress_is_scale_covariant_and_preserves_strain_admission() {
    let expected=model(1.0).plane_stress(20.0).unwrap().engineering();
    for scale in [1e-180,1e180] {
        let ps=model(scale).plane_stress(20.0).unwrap();let actual=ps.engineering();
        for i in 0..3 {for j in 0..3 {
            cclose(C64::new(actual[i][j].re/scale,actual[i][j].im/scale),expected[i][j],2e-12);
        }}
        assert!(ps.strain_tensor([C64::new(f64::NAN,0.0);3]).is_err());
        assert!(ps.strain_tensor([C64::ONE;3]).is_err());
    }
}
