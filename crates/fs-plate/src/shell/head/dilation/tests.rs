use super::*;
use crate::shell::{head::TensionedDiskSpec, profile::ProfileBudget};
fn budget()->MembraneReductionBudget {
    MembraneReductionBudget {max_modes:8,max_nodes:100,max_facet_pairs:2000,
        max_solve_entries:100_000,relative_tolerance:1e-6}
}
fn fixture()->(TensionedDisk,Vec<ModePair>) {
    let f=TensionedDisk::new(TensionedDiskSpec {radius_m:0.17,thickness_m:0.000254,
        young_pa:4e9,poisson:0.38,density_kg_m3:1390.0,tension_n_m:3000.0,
        radial_intervals:2,azimuths:8},ProfileBudget {max_nodes:100,max_triangles:100,max_feature_evaluations:0}).unwrap();
    let n=f.model.free;let mut k=vec![0.0;n*n];let mut m=k.clone();
    for i in 0..n {for j in 0..n {k[i*n+j]=f.model.k.get(i,j);m[i*n+j]=f.model.m.get(i,j);}}
    let modes=fs_modal::eigh_gen_dense(&k,&m,n).unwrap().into_iter().take(3).collect();(f,modes)
}
fn solve(n:usize,k:&[f64],rhs:&[Vec<f64>])->Result<Vec<Vec<f64>>,PlateError> {
    let mut id=vec![0.0;n*n];for i in 0..n {id[i*n+i]=1.0;}
    let pairs=fs_modal::eigh_gen_dense(k,&id,n).unwrap();
    Ok(rhs.iter().map(|b|{let mut u=vec![0.0;n];for p in &pairs {
        assert!(p.lambda>0.0);let a=p.phi.iter().zip(b).map(|(p,b)|p*b).sum::<f64>()/p.lambda;
        for (u,p) in u.iter_mut().zip(&p.phi) {*u+=a*p;}
    }u}).collect())
}
fn reduced(f:&TensionedDisk,m:&[ModePair])->DilationReduction {DilationReduction::new(f,m,0.8,4000.0,budget(),solve).unwrap()}
#[test]
fn fixed_rim_limit_keeps_original_head_but_reports_its_reciprocal_reaction() {
    let (f,m)=fixture();let d=reduced(&f,&m);let q=[1e-4,-2e-4,0.6e-4,0.0];
    assert_eq!(d.potential(&q),d.membrane.potential(&q[..3]));
    let mut original=[0.0;3];let mut joint=[0.0;4];
    d.membrane.gradient(&q[..3],&mut original);d.gradient(&q,&mut joint);
    assert_eq!(&joint[..3],&original);assert!(joint[3]>0.0);
    assert_eq!(d.baseline_tension_n_m(&q),f.spec.tension_n_m);
}
#[test]
fn rim_stiffness_comes_from_actual_polygon_area_and_no_extra_inplane_solve() {
    let (f,m)=fixture();let mut calls=0;
    let d=DilationReduction::new(&f,&m,0.8,4000.0,budget(),|n,k,b|{calls+=1;solve(n,k,b)}).unwrap();
    assert_eq!(calls,1);
    let area=f.mass_kg/(f.spec.density_kg_m3*f.spec.thickness_m);
    let expected=(2.0*f.spec.young_pa*f.spec.thickness_m/(1.0-f.spec.poisson)
        *area/f.spec.radius_m.powi(2)+4000.0)/0.8;
    assert!((d.rim_omega.powi(2)-expected).abs()<1e-12*expected);
    let q=[0.0,0.0,0.0,0.8_f64.sqrt()*2e-5];let mut g=[0.0;4];d.gradient(&q,&mut g);
    assert!((g[3]-expected*q[3]).abs()<1e-12*(expected*q[3]).abs());
}
#[test]
fn imposed_edge_extension_matches_the_actual_prestress_pencil_increment() {
    let (f,m)=fixture();let d=reduced(&f,&m);let mut a=[0.0;16];let mut b=a;let mut scratch=[[0.0;3];4];
    d.hessian_matrix(&[0.0;4],&mut a,&mut scratch);
    let s=2e-5;let q=[0.0,0.0,0.0,s/d.force_weight()];d.hessian_matrix(&q,&mut b,&mut scratch);
    let tension=f.spec.young_pa*f.spec.thickness_m/(1.0-f.spec.poisson)*s/f.spec.radius_m;
    let mut geometric=[[0.0;3];3];
    for t in &f.mesh.tris {
        let p=t.map(|i|f.mesh.nodes[i]);let a2=(p[1].0-p[0].0)*(p[2].1-p[0].1)-(p[2].0-p[0].0)*(p[1].1-p[0].1);
        let grad:[[f64;2];3]=core::array::from_fn(|i|[(p[(i+1)%3].1-p[(i+2)%3].1)/a2,(p[(i+2)%3].0-p[(i+1)%3].0)/a2]);
        let gradients:Vec<[f64;2]>=m.iter().map(|m|core::array::from_fn(|c|(0..3).map(|i|
            grad[i][c]*f.model.dof_map[3*t[i]].map_or(0.0,|k|m.phi[k])).sum())).collect();
        for i in 0..3 {for j in 0..3 {geometric[i][j]+=0.5*a2*tension*(gradients[i][0]*gradients[j][0]+gradients[i][1]*gradients[j][1]);}}
    }
    let scale=geometric.iter().flatten().fold(1.0_f64,|s,v|s.max(v.abs()));
    for i in 0..3 {for j in 0..3 {assert!((b[4*i+j]-a[4*i+j]-geometric[i][j]).abs()<1e-7*scale);}}
}
#[test]
fn complete_energy_gradient_and_mixed_tangent_agree_without_allocated_hot_scratch() {
    let (f,m)=fixture();let d=reduced(&f,&m);let q=[1e-4,-1.3e-4,0.7e-4,3e-5];
    let mut g=[0.0;4];d.gradient(&q,&mut g);let mut a=[0.0;16];let mut scratch=[[0.0;3];4];d.hessian_matrix(&q,&mut a,&mut scratch);
    for col in 0..4 {
        let h=1e-8;let mut plus=q;let mut minus=q;plus[col]+=h;minus[col]-=h;
        let fd=(d.potential(&plus)-d.potential(&minus))/(2.0*h);
        assert!((g[col]-fd).abs()<1e-6*g[col].abs().max(1.0));
        let mut gp=[0.0;4];let mut gm=gp;d.gradient(&plus,&mut gp);d.gradient(&minus,&mut gm);
        let mut direction=[0.0;4];direction[col]=1.0;let mut product=[0.0;4];d.hessian_vector(&q,&direction,&mut product);
        for row in 0..4 {
            let scale=a[4*row+col].abs().max(1.0);
            assert!((product[row]-a[4*row+col]).abs()<1e-10*scale);
            assert!(((gp[row]-gm[row])/(2.0*h)-a[4*row+col]).abs()<1e-6*scale);
            assert_eq!(a[4*row+col],a[4*col+row]);
        }
    }
    assert!(a[3].abs()>0.0);assert!(d.potential(&q)>0.0);
}
#[test]
fn mass_normalization_changes_coordinates_not_physical_energy_or_force_work() {
    let (f,m)=fixture();let a=reduced(&f,&m);let b=DilationReduction::new(&f,&m,3.2,4000.0,budget(),solve).unwrap();
    let q=[1e-4,2e-4,-1e-4,1e-5/a.force_weight()];let mut p=q;p[3]=1e-5/b.force_weight();
    assert!((a.potential(&q)-b.potential(&p)).abs()<1e-12*a.potential(&q));
    let mut ga=[0.0;4];let mut gb=ga;a.gradient(&q,&mut ga);b.gradient(&p,&mut gb);
    assert!((ga[3]/a.force_weight()-gb[3]/b.force_weight()).abs()<1e-10*ga[3].abs());
    assert!((5.0*a.force_weight()*q[3]-5.0e-5).abs()<1e-18);
}
#[test]
fn invalid_actuator_material_and_output_dimensions_refuse() {
    let (mut f,m)=fixture();
    for mass in [0.0,-1.0,f64::NAN,f64::INFINITY] {assert!(DilationReduction::new(&f,&m,mass,0.0,budget(),solve).is_err());}
    let d=reduced(&f,&m);let mut out=[0.0;4];d.gradient(&[0.0;3],&mut out);assert!(out.iter().all(|v|v.is_nan()));
    assert!(d.clock_bounds(0.0).is_err());assert!(d.clock_bounds(f64::NAN).is_err());
    let c=d.clock_bounds(0.001).unwrap();assert!(c[..3].iter().zip(d.membrane.omegas()).all(|(a,b)|a>b));
    f.section.d[0]*=1.1;assert!(DilationReduction::new(&f,&m,1.0,0.0,budget(),solve).is_err());
}
