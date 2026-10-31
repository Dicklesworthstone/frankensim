use super::*;

fn diagonal(value: f64) -> Tangent6 { core::array::from_fn(|i| core::array::from_fn(|j| if i==j {value} else {0.0})) }
fn model() -> TensorProny {
    let mut a = [[0.0;6];6]; a[0] = [2.0,-1.0,0.0,3.0,0.5,0.0];
    let mut b = [[0.0;6];6]; b[2] = [0.0,1.0,2.0,0.0,-1.0,0.5];
    TensorProny::new_mandel(diagonal(12.0), vec![
        TensorPronyBranch::from_gram_rows(a,0.01).unwrap(),
        TensorPronyBranch::from_gram_rows(b,0.3).unwrap()],0.1).unwrap()
}
fn columns() -> Vec<Voigt> { vec![[1.0,0.0,-0.2,0.3,0.0,0.0],[-0.1,0.5,0.0,-0.4,0.2,0.0]] }
fn budget() -> TensorProjectionBudget { TensorProjectionBudget {max_coordinates:4,max_scalar_terms:10_000} }
fn close(a:f64,b:f64) { assert!((a-b).abs() < 2e-12*a.abs().max(b.abs()).max(1e-12),"{a} vs {b}"); }
fn projected_stiffness(p:&ProjectedTensorProny,omega:f64)->(Vec<f64>,Vec<f64>) {
    let n=p.coordinates(); let mut re=p.equilibrium_stiffness().to_vec(); let mut im=vec![0.0;n*n];
    for b in p.branches() {
        let z=omega*b.relaxation_time_s(); let den=1.0+z*z;
        for i in 0..n {for j in 0..n {
            let value=b.rows().iter().fold(0.0,|s,row|s+row[i]*row[j]);
            re[i*n+j]+=value*z*z/den; im[i*n+j]+=value*z/den;
        }}
    }
    (re,im)
}

#[test]
fn projected_tensor_preserves_full_complex_work_and_rank_deficient_branches() {
    let m=model(); let c=columns(); let p=m.project_strain(&c,0.25,budget(),||false).unwrap();
    assert_eq!(p.coordinates(),2); assert_eq!(p.branches().len(),2);
    assert!(p.branches()[0].rows()[1..].iter().flatten().all(|x|*x==0.0));
    for omega in [0.0,0.1,10.0,100.0,10_000.0] {
        let (storage,loss)=m.modulus_mandel(omega).unwrap();
        let (re,im)=projected_stiffness(&p,omega);
        for i in 0..2 {for j in 0..2 {
            close(re[2*i+j],0.25*dot(&mandel(&c[i]),&apply(&storage,&mandel(&c[j]))));
            close(im[2*i+j],0.25*dot(&mandel(&c[i]),&apply(&loss,&mandel(&c[j]))));
        }}
    }
    let x=[0.002,-0.003]; let strain=p.strain(&x).unwrap();
    let mut history=m.history(); let report=history.advance(strain,0.0,||false).unwrap();
    let force=p.conjugate_force(&report.stress).unwrap();
    close(0.5*(force[0]*x[0]+force[1]*x[1]),0.25*report.stored_energy);
    assert_ne!(p.equilibrium_stiffness()[1],0.0,"cross terms must survive");
}

#[test]
fn projected_factor_histories_match_the_original_material_on_reversing_strains() {
    let m=model(); let p=m.project_strain(&columns(),0.25,budget(),||false).unwrap();
    let mut history=m.history(); let mut z=vec![[0.0;6];p.branches().len()];
    let mut old=[0.0;2];
    for step in 0..100 {
        let x=[0.004*det::sin(step as f64*0.37),0.002*det::cos(step as f64*0.29)];
        let dt=if step==0 {0.0} else {0.003};
        let native=history.advance(p.strain(&x).unwrap(),dt,||false).unwrap();
        let mut force:Vec<_>=(0..2).map(|i|(0..2).map(|j|p.equilibrium_stiffness()[2*i+j]*x[j]).sum::<f64>()).collect();
        let mut stored=0.5*(force[0]*x[0]+force[1]*x[1]); let mut heat=0.0;
        for (z,branch) in z.iter_mut().zip(p.branches()) {
            for (k,row) in branch.rows().iter().enumerate() {
                let delta=row.iter().enumerate().fold(0.0,|s,(i,a)|s+a*(x[i]-old[i]));
                let (next,_,loss)=ramp(z[k],delta,dt/branch.relaxation_time_s());
                z[k]=next;heat+=loss;stored+=0.5*next*next;
                for i in 0..2 {force[i]+=row[i]*next;}
            }
        }
        let expected=p.conjugate_force(&native.stress).unwrap();
        for i in 0..2 {close(force[i],expected[i]);}
        close(stored,0.25*native.stored_energy); close(heat,0.25*native.dissipated);
        old=x;
    }
}

#[test]
fn projection_preserves_shear_work_and_volume_scaling() {
    let m=TensorProny::new_mandel(diagonal(2.0),vec![],0.1).unwrap();
    let c=[[0.0,0.0,0.0,1.0,0.0,0.0]];
    let p=m.project_strain(&c,0.5,budget(),||false).unwrap();
    close(p.equilibrium_stiffness()[0],2.0); // 2*V*C_xyxy, not V*C_xyxy.
    let large=m.project_strain(&c,2.0,budget(),||false).unwrap();
    close(large.equilibrium_stiffness()[0],4.0*p.equilibrium_stiffness()[0]);
    let stress=[0.0,0.0,0.0,7.0,0.0,0.0];
    close(p.conjugate_force(&stress).unwrap()[0],7.0);
    assert!(p.strain(&[0.08]).is_err(),"Frobenius shear norm exceeds 0.1");
}

#[test]
fn projection_refuses_domain_work_and_cancellation_without_mutating_the_law() {
    let m=model(); let saved=m.clone(); let c=columns();
    for v in [0.0,-1.0,f64::NAN,f64::INFINITY] {assert!(m.project_strain(&c,v,budget(),||false).is_err());}
    assert!(m.project_strain(&[],1.0,budget(),||false).is_err());
    assert!(m.project_strain(&c,1.0,TensorProjectionBudget {max_scalar_terms:1,..budget()},||false).is_err());
    let mut bad=c.clone();bad[0][0]=f64::NAN;assert!(m.project_strain(&bad,1.0,budget(),||false).is_err());
    let mut polls=0;let baseline=m.project_strain(&c,1.0,budget(),||{polls+=1;false}).unwrap();
    for stop in 1..=polls {
        let mut at=0;assert!(matches!(m.project_strain(&c,1.0,budget(),||{at+=1;at==stop}),Err(TensorPronyError::Cancelled)));
    }
    assert_eq!(m,saved);
    assert_eq!(baseline,m.project_strain(&c,1.0,budget(),||false).unwrap());
    assert!(baseline.strain(&[f64::NAN,0.0]).is_err());
    assert!(baseline.strain(&[0.0]).is_err());
    assert!(baseline.conjugate_force(&[f64::NAN;6]).is_err());
}

#[test]
fn thermal_projection_retains_signed_full_tensor_derivatives_and_original_strain_work() {
    let m = model();
    let columns = columns();
    let volume = 0.25;
    let p = m.project_strain(&columns,volume,budget(),||false).unwrap();
    let mut equilibrium = diagonal(-0.7);
    equilibrium[0][3] = 1.3; equilibrium[3][0] = 1.3;
    let mut first = diagonal(0.2);
    first[0][0] = -1.1; first[1][4] = -0.8; first[4][1] = -0.8;
    let mut second = diagonal(-0.3);
    second[2][5] = 0.9; second[5][2] = 0.9;
    let derivative = source::thermal::TensorTemperatureDerivative {
        equilibrium_pa_k: equilibrium, branches_pa_k: vec![first,second],
        interval_k: [280.0,320.0],
    };
    let projected = p.thermal_stiffness_derivatives(&derivative,budget(),||false).unwrap();
    let pairs = [(equilibrium,projected.equilibrium()),
        (first,projected.branches()[0].as_slice()),
        (second,projected.branches()[1].as_slice())];
    // Independent physical Voigt contraction: both shear factors sqrt(2)
    // remain, including every normal/shear cross term and both tangent axes.
    let weights = [1.0,1.0,1.0,2.0_f64.sqrt(),2.0_f64.sqrt(),2.0_f64.sqrt()];
    for (matrix, reduced) in pairs {
        for x in [[0.004,-0.003],[-0.002,0.005]] {
            let strain: [f64;6] = std::array::from_fn(|i|
                weights[i]*(columns[0][i]*x[0]+columns[1][i]*x[1]));
            let direct = volume*(0..6).map(|i| (0..6)
                .map(|j| strain[i]*matrix[i][j]*strain[j]).sum::<f64>()).sum::<f64>();
            let reduced_work = reduced[0]*x[0]*x[0] + 2.0*reduced[1]*x[0]*x[1]
                + reduced[3]*x[1]*x[1];
            close(direct,reduced_work);
            assert_eq!(reduced[1],reduced[2]);
        }
    }
    let viscous = m.with_kelvin_viscosity_mandel(diagonal(0.2)).unwrap()
        .project_strain(&columns,volume,budget(),||false).unwrap();
    assert_eq!(projected,viscous.thermal_stiffness_derivatives(&derivative,budget(),||false).unwrap());
    assert!(p.thermal_stiffness_derivatives(&derivative,
        TensorProjectionBudget {max_scalar_terms:503,..budget()},||false).is_err());
    assert_eq!(projected,p.thermal_stiffness_derivatives(&derivative,
        TensorProjectionBudget {max_scalar_terms:504,..budget()},||false).unwrap());
    assert!(matches!(p.thermal_stiffness_derivatives(&derivative,budget(),||true),
        Err(TensorPronyError::Cancelled)));
    let mut incomplete = derivative; incomplete.branches_pa_k.pop();
    assert!(p.thermal_stiffness_derivatives(&incomplete,budget(),||false).is_err());
}
