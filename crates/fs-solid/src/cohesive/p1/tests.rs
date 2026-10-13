use super::*;
use crate::cohesive::{Budget, CohesiveSystem, Config};
use fs_solver::{LinearOp, NewtonKrylovConfig};

fn mesh() -> (Vec<[f64; 3]>, Vec<[usize; 4]>) {
    (vec![[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,-1.],
        [0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]], vec![[0,1,2,3],[4,5,6,7]])
}
fn law() -> BilinearCohesive { BilinearCohesive::new(1e9,2e5,100.,2e9,0.01).unwrap() }
fn pair(law: &BilinearCohesive) -> FacePair<'_> {
    FacePair { id: 41, side_a: [2,0,1], side_b: [5,6,4], law }
}
fn budget() -> TraceBudget { TraceBudget { maximum_nodes: 100, maximum_tetrahedra: 100,
    maximum_interfaces: 10, minimum_face_quality: 1e-10 } }
fn close(a: f64, b: f64, tol: f64) { assert!((a-b).abs() <= tol, "{a:e} != {b:e}"); }
fn build<'a>(nodes: &[[f64;3]], tets: &[[usize;4]], free: &[usize], pairs: &[FacePair<'a>]) -> Result<P1Interfaces<'a>,Error> {
    P1Interfaces::new_controlled(nodes,tets,free,pairs,budget(),&mut || false)
}

#[test]
fn triangle_quadrature_matches_consistent_interface_matrix_and_permutations() {
    let (nodes,tets) = mesh(); let law = law(); let free: Vec<_> = (0..24).collect();
    let faces = build(&nodes,&tets,&free,&[pair(&law)]).unwrap();
    let mut permutation = pair(&law); permutation.side_a = [1,2,0]; permutation.side_b = [6,4,5];
    let again = build(&nodes,&tets,&free,&[permutation]).unwrap();
    assert_eq!(faces.points(),again.points());
    close(faces.points().iter().map(|p| p.area_weight_m2).sum(),0.5,1e-15);
    for point in faces.points() { assert_eq!(point.axes,[[0.,0.,1.],[1.,0.,0.],[0.,1.,0.]]); }
    // Independent analytic triangle mass entries: A*(1+delta_ij)/12.
    for j in 0..3 { for c in 0..3 {
        let mut u = vec![0.0;24]; u[3*(4+j)+c] = 1.0;
        let jumps = faces.jumps(&u,&mut || false).unwrap();
        let f = faces.full_forces(&jumps,&mut || false).unwrap();
        for i in 0..3 { for d in 0..3 {
            let expected = if c == d { 0.5*(if i == j { 2.0 } else { 1.0 })/12.0 } else { 0.0 };
            close(f[3*(4+i)+d],expected,2e-16);
            close(f[3*i+d],-expected,2e-16);
        } }
    } }
}

#[test]
fn jump_force_actions_preserve_virtual_work_force_moment_and_rigid_motion() {
    let (nodes,tets) = mesh(); let law = law();
    let faces = build(&nodes,&tets,&(0..24).collect::<Vec<_>>(),&[pair(&law)]).unwrap();
    let u: Vec<_> = (0..24).map(|i| (i as f64*0.7).sin()*0.01).collect();
    let traction = [[17.,-31.,12.],[-20.,11.,30.],[6.,4.,-7.]];
    let jumps = faces.jumps(&u,&mut || false).unwrap();
    let force = faces.full_forces(&traction,&mut || false).unwrap();
    let work: f64 = force.iter().zip(&u).map(|(a,b)| a*b).sum();
    let local: f64 = jumps.iter().zip(traction).zip(faces.points()).map(|((j,t),p)| p.area_weight_m2*dot(*j,t)).sum();
    close(work,local,1e-15);
    let mut resultant = [0.0;3]; let mut moment = [0.0;3]; let mut rigid = Vec::new();
    for (node,f) in nodes.iter().zip(force.chunks_exact(3)) {
        let m = cross(*node,[f[0],f[1],f[2]]);
        for c in 0..3 { resultant[c] += f[c]; moment[c] += m[c]; }
        let motion = cross([0.003,0.001,-0.002],*node);
        rigid.extend(motion.map(|v| v+0.01));
    }
    for v in resultant.into_iter().chain(moment) { close(v,0.,1e-14); }
    for jump in faces.jumps(&rigid,&mut || false).unwrap() { for v in jump { close(v,0.,1e-17); } }
}

#[test]
fn prescribed_and_free_motions_reconstruct_the_identical_full_jump() {
    let (nodes,tets) = mesh(); let law = law(); let free = [0,4,8,12,16,20];
    let faces = build(&nodes,&tets,&free,&[pair(&law)]).unwrap();
    let reduced: Vec<_> = (0..free.len()).map(|i| i as f64*0.001).collect();
    let mut fixed: Vec<_> = (0..24).map(|i| i as f64*0.0001).collect();
    for &dof in &free { fixed[dof] = 0.0; }
    let full = faces.expand(&reduced,&fixed,&mut || false).unwrap();
    let exact = faces.jumps(&full,&mut || false).unwrap();
    let offsets = faces.prescribed_jumps(&fixed,&mut || false).unwrap();
    for ((patch,offset),expected) in faces.patches.iter().zip(offsets).zip(exact) {
        for (a,b) in patch.jump(&reduced,offset).unwrap().into_iter().zip(expected) { close(a,b,1e-17); }
    }
    fixed[free[0]] = 1.0;
    assert!(faces.prescribed_jumps(&fixed,&mut || false).is_err());
}

#[test]
fn fully_prescribed_faces_still_fracture_and_keep_compressive_reclosure() {
    struct Identity;
    impl LinearOp for Identity { fn n(&self)->usize { 6 } fn apply(&self,x:&[f64],y:&mut[f64]) { y.copy_from_slice(x); } }
    let (nodes,tets) = mesh(); let law = law();
    let faces = build(&nodes,&tets,&[9,10,11,21,22,23],&[pair(&law)]).unwrap();
    assert!(faces.patches().iter().all(|p| p.terms.is_empty()));
    let op = Identity;
    let system = CohesiveSystem::new(&op,faces.patches().to_vec(),Budget { max_dofs:6,max_patches:3,
        max_jump_terms:0,max_krylov_basis_entries:4096 }).unwrap();
    let config = Config { newton:NewtonKrylovConfig::default(),max_newton_iterations:10,
        max_operator_evaluations:100,force_tolerance_n:1e-8,local_energy_tolerance_j_m2:1e-7,
        total_cohesive_energy_tolerance_j:1e-7,maximum_displacement_m:0.01 };
    let state = system.bonded_state(); let mut fixed = vec![0.;24];
    for node in [4,5,6] { fixed[3*node+2] = 0.002; }
    let offsets = faces.prescribed_jumps(&fixed,&mut || false).unwrap();
    let result = state.solve_controlled(&[0.;6],&offsets,config,&mut || false).unwrap();
    assert!(result.state.points().iter().all(|p| p.failed()));
    assert!(state.points().iter().all(|p| p.damage() == 0.0));
    close(result.cohesive_energy_j.dissipated,50.,1e-12);
    for node in [4,5,6] { fixed[3*node+2] = -1e-5; }
    let offsets = faces.prescribed_jumps(&fixed,&mut || false).unwrap();
    let closed = result.state.solve_controlled(&[0.;6],&offsets,config,&mut || false).unwrap();
    for p in closed.state.points() { assert!(p.failed()); close(p.traction_pa()[0],-2e4,1e-10); }
}

#[test]
fn invalid_geometry_topology_controls_and_cancellation_refuse_without_panicking() {
    let (nodes,tets) = mesh(); let law = law(); let p = pair(&law); let free:Vec<_>=(0..24).collect();
    assert!(build(&nodes,&tets,&free,&[p,p]).is_err());
    let mut bad = nodes.clone(); bad[4][0] = 1e-12;
    assert!(build(&bad,&tets,&free,&[p]).is_err());
    bad = nodes.clone(); bad[7][2] = -1.0;
    assert!(build(&bad,&tets,&free,&[p]).is_err());
    bad = nodes.clone(); bad[1] = bad[0];
    assert!(build(&bad,&tets,&free,&[p]).is_err());
    bad = nodes.clone(); bad[7][0] = f64::NAN;
    assert!(build(&bad,&tets,&free,&[p]).is_err());
    let mut repeated = tets.clone(); repeated.push(tets[0]);
    assert!(build(&nodes,&repeated,&free,&[p]).is_err());
    assert!(build(&nodes,&tets,&[2,1],&[p]).is_err());
    assert!(build(&nodes,&tets,&[24],&[p]).is_err());
    let mut cap = budget(); cap.maximum_interfaces = 0;
    assert!(P1Interfaces::new_controlled(&nodes,&tets,&free,&[p],cap,&mut || false).is_err());
    for stop in [1,5,8] {
        let mut polls = 0;
        let result = P1Interfaces::new_controlled(&nodes,&tets,&free,&[p],budget(),&mut || { polls+=1; polls==stop });
        assert!(matches!(result,Err(Error::Cancelled)));
    }
}
