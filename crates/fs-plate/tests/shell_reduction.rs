use fs_plate::{PlateSection, ModePair};
use fs_plate::shell::{ShellMesh,ShellSupport,ShellModel,assemble_shell_sections};
use fs_plate::shell::reduction::{ShellReduction,ReductionBudget};
fn fixture()->(ShellMesh,Vec<PlateSection>,ShellModel,Vec<ModePair>) {
    let mesh=ShellMesh::new(vec![[0.0,0.0,0.0],[0.1,0.0,0.01],[0.0,0.1,0.02],[0.08,0.09,0.035]],vec![[0,1,3],[0,3,2]]).unwrap();
    let sections=vec![PlateSection::isotropic(112.6e9,0.342,0.001,8607.0).unwrap(),
        PlateSection::isotropic(112.6e9,0.342,0.0007,8607.0).unwrap()];
    let model=assemble_shell_sections(&mesh,&sections,&[0,1,2],ShellSupport::Clamped).unwrap();
    let n=model.free;let mut k=vec![0.0;n*n];let mut m=k.clone();
    for i in 0..n {for j in 0..n {k[i*n+j]=model.k.get(i,j);m[i*n+j]=model.m.get(i,j);}}
    let modes=fs_modal::eigh_gen_dense(&k,&m,n).unwrap();
    (mesh,sections,model,modes)
}
fn budget()->ReductionBudget {ReductionBudget{max_modes:12,max_facet_modes:100,relative_tolerance:1e-6}}
#[test]
fn nonlinear_gradient_matches_independent_metric_energy_and_finite_differences() {
    let (mesh,sections,model,modes)=fixture();
    let reduced=ShellReduction::new(&mesh,&sections,&model,&modes,budget()).unwrap();
    let n=modes.len();let q:Vec<_>=(0..n).map(|i|1e-5*((i*7+3) as f64).sin()).collect();
    let mut force=vec![0.0;n];reduced.gradient(&q,&mut force);
    assert!(reduced.potential(&q)>0.0);
    for i in 0..n {
        let h=1e-9;let mut plus=q.clone();let mut minus=q.clone();plus[i]+=h;minus[i]-=h;
        let fd=(reduced.potential(&plus)-reduced.potential(&minus))/(2.0*h);
        assert!((fd-force[i]).abs()<=2e-6*fd.abs().max(force[i].abs()).max(1.0));
    }
    // Independent full nodal reconstruction and local metric, subtracting only
    // the linear membrane already present in the independently assembled K.
    let mut u=vec![0.0;6*mesh.nodes.len()];
    for (dof,free) in model.dof_map.iter().enumerate() {if let Some(free)=free {
        u[dof]=modes.iter().zip(&q).map(|(m,q)|m.phi[*free]*q).sum();
    }}
    let linear=0.5*modes.iter().zip(&q).map(|(m,q)|m.lambda*q*q).sum::<f64>();
    let mut correction=0.0;
    for (e,tri) in mesh.tris.iter().enumerate() {
        let g=mesh.facet(e).unwrap();let mut dx=[0.0;3];let mut dy=[0.0;3];
        for a in 0..3 {for c in 0..3 {dx[c]+=g.gradient[a][0]*u[6*tri[a]+c];dy[c]+=g.gradient[a][1]*u[6*tri[a]+c];}}
        let dot=|a:[f64;3],b:[f64;3]|a.iter().zip(b).map(|(x,y)|x*y).sum::<f64>();
        let old=[dot(g.frame[0],dx),dot(g.frame[1],dy),dot(g.frame[0],dy)+dot(g.frame[1],dx)];
        let new=[old[0]+0.5*dot(dx,dx),old[1]+0.5*dot(dy,dy),old[2]+dot(dx,dy)];
        let s=&sections[e];let energy=|v:[f64;3]| {
            0.5*g.area_m2*12.0/(s.thickness*s.thickness)*(0..3).map(|i|(0..3).map(|j|v[i]*s.d[3*i+j]*v[j]).sum::<f64>()).sum::<f64>()
        };
        correction+=energy(new)-energy(old);
    }
    assert!((reduced.potential(&q)-(linear+correction)).abs()<1e-6*(linear+correction).abs());
}
#[test]
fn small_amplitude_recovers_linear_pencil_and_curvature_breaks_sign_symmetry() {
    let (mesh,sections,model,modes)=fixture();let r=ShellReduction::new(&mesh,&sections,&model,&modes,budget()).unwrap();
    let mut q=vec![0.0;modes.len()];q[0]=1e-10;let mut g=q.clone();r.gradient(&q,&mut g);
    assert!((g[0]/q[0]-modes[0].lambda).abs()<2e-5*modes[0].lambda);
    q[0]=1e-4;let positive=r.potential(&q);q[0]=-q[0];let negative=r.potential(&q);
    assert!((positive-negative).abs()>1e-6*positive.max(negative),"curved shell must retain its cubic potential");
}
#[test]
fn force_projection_is_reciprocal_and_reduction_does_not_accept_invented_frequencies() {
    let (mesh,sections,model,mut modes)=fixture();let r=ShellReduction::new(&mesh,&sections,&model,&modes,budget()).unwrap();
    let b=r.point_port(0,[0.2,0.3,0.5],[0.0,0.0,1.0]).unwrap();
    let v:Vec<_>=(0..modes.len()).map(|i|0.01*i as f64).collect();let force=12.0;
    let generalized=b.iter().zip(&v).map(|(b,v)|b*force*v).sum::<f64>();
    let velocity=b.iter().zip(&v).map(|(b,v)|b*v).sum::<f64>();
    assert!((generalized-force*velocity).abs()<1e-12*generalized.abs().max(1.0));
    assert!(r.point_port(0,[0.0,0.0,1.0],[0.0,0.0,2.0]).is_err());
    modes[0].lambda*=2.0;
    assert!(ShellReduction::new(&mesh,&sections,&model,&modes,budget()).is_err());
    let mut cap=budget();cap.max_facet_modes=1;
    assert!(ShellReduction::new(&mesh,&sections,&model,&modes,cap).is_err());
    let mut g=vec![0.0;1];r.gradient(&[],&mut g);assert!(g[0].is_nan());
}
