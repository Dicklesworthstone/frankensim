use super::*;
use fs_material::{IsotropicElastic, visco_tensor::TensorPronyBranch};

fn specimen(square: bool) -> TensorMembrane {
    let nodes=if square {vec![[0.0,0.0,0.0],[1.0,0.0,0.0],[1.0,1.0,0.0],[0.0,1.0,0.0]]}
        else {vec![[0.0,0.0,0.0],[1.0,0.0,0.0],[0.0,1.0,0.0]]};
    let tris=if square {vec![[0,1,2],[0,2,3]]} else {vec![[0,1,2]]};
    let c=IsotropicElastic::new(1000.0,0.25,0.2).unwrap().stiffness();
    let material=MembraneFacetMaterial { law:TensorProny::from_tensor_shear(c,
        vec![TensorPronyBranch::from_stiffness_mandel(c.map(|r|r.map(|v|0.4*v)),0.02).unwrap()],0.2).unwrap(),
        thickness_m:0.1,density_kg_m3:10.0 };
    let mut fixed=vec![0,1,4];fixed.extend((0..nodes.len()).map(|n|3*n+2));
    let nf=tris.len();
    TensorMembrane::new(ShellMesh::new(nodes,tris).unwrap(),vec![material;nf],&fixed,
        TensorProjectionBudget {max_coordinates:128,max_scalar_terms:100_000_000},
        TensorMaterialBudget {max_state_coordinates:128,max_material_points:128},
        &CancelGate::new_clock_free()).unwrap()
}
fn dot(a:[f64;3],b:[f64;3])->f64 {a.into_iter().zip(b).map(|(a,b)|a*b).sum()}
fn cross(a:[f64;3],b:[f64;3])->[f64;3] {
    [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]
}
fn close(a:f64,b:f64) {assert!((a-b).abs()<1e-11*(1.0+a.abs().max(b.abs())),"{a} != {b}");}
fn interpolate(values:&[[f64;3];3],n:[f64;3])->[f64;3] {
    core::array::from_fn(|axis|(0..3).map(|i|n[i]*values[i][axis]).sum())
}

#[test]
fn linear_surface_and_edge_loads_preserve_resultants_moments_and_virtual_work() {
    let mesh=specimen(false);let gate=CancelGate::new_clock_free();
    let traction=[[[1.0,2.0,-3.0],[-4.0,5.0,6.0],[7.0,-8.0,9.0]]];
    let edge=MembraneEdgeTraction {nodes:[0,1],values_n_m:[[2.0,-3.0,4.0],[-1.0,2.0,5.0]]};
    let force=mesh.distributed_nodal_forces(&traction,&[edge],[0.0;3],&gate).unwrap();
    let velocity=[[0.0;3],[0.2,0.0,0.0],[-0.1,0.3,0.0]];
    let points:[[f64;3];3]=mesh.mesh().nodes.clone().try_into().unwrap();
    let mut expected_force=[0.0;3];let mut expected_moment=[0.0;3];let mut work=0.0;
    // Independent exact degree-two triangle quadrature, not the nodal formula.
    for n in [[0.5,0.5,0.0],[0.0,0.5,0.5],[0.5,0.0,0.5]] {
        let t=interpolate(&traction[0],n);let x=interpolate(&points,n);let v=interpolate(&velocity,n);
        let moment=cross(x,t);let weight=0.5/3.0;
        for a in 0..3 {expected_force[a]+=weight*t[a];expected_moment[a]+=weight*moment[a];}
        work+=weight*dot(t,v);
    }
    // Independent two-point Gauss rule on the unit-length boundary edge.
    let half=0.5/fs_math::det::sqrt(3.0);
    for s in [0.5-half,0.5+half] {
        let t=core::array::from_fn(|a|(1.0-s)*edge.values_n_m[0][a]+s*edge.values_n_m[1][a]);
        let v=core::array::from_fn(|a|(1.0-s)*velocity[0][a]+s*velocity[1][a]);
        let x=[s,0.0,0.0];let moment=cross(x,t);
        for a in 0..3 {expected_force[a]+=0.5*t[a];expected_moment[a]+=0.5*moment[a];}
        work+=0.5*dot(t,v);
    }
    let mut moment=[0.0;3];
    for (x,f) in points.into_iter().zip(&force) {let m=cross(x,*f);for a in 0..3 {moment[a]+=m[a];}}
    for a in 0..3 {close(force.iter().map(|f|f[a]).sum(),expected_force[a]);close(moment[a],expected_moment[a]);}
    close(force.iter().copied().zip(velocity).map(|(f,v)|dot(f,v)).sum(),work);
    let u=mesh.restrict_state(&velocity).unwrap();let free=mesh.project_nodal_forces(&force).unwrap();
    close(free.iter().zip(u).map(|(f,u)|f*u).sum(),work);
    let reverse=MembraneEdgeTraction {nodes:[1,0],values_n_m:[edge.values_n_m[1],edge.values_n_m[0]]};
    assert_eq!(force,mesh.distributed_nodal_forces(&traction,&[reverse],[0.0;3],&gate).unwrap());
}

#[test]
fn pressure_sign_body_acceleration_and_boundary_admission_remain_physical() {
    let mesh=specimen(true);let gate=CancelGate::new_clock_free();
    let zero=vec![[[0.0;3];3];2];let acceleration=[2.0,-3.0,4.0];
    let gravity=mesh.distributed_nodal_forces(&zero,&[],acceleration,&gate).unwrap();
    for (f,m) in gravity.iter().zip(mesh.nodal_masses_kg()) {
        for a in 0..3 {close(f[a],m*acceleration[a]);}
    }
    let pressure=[[1.0,2.0,3.0],[4.0,5.0,6.0]];
    let f=mesh.pressure_nodal_forces(&pressure,&gate).unwrap();
    close(f.iter().map(|f|f[2]).sum(),-3.5);
    for f in f {close(f[0],0.0);close(f[1],0.0);}
    let edge=MembraneEdgeTraction {nodes:[0,1],values_n_m:[[1.0;3];2]};
    assert!(mesh.distributed_nodal_forces(&zero,&[edge,edge],[0.0;3],&gate).is_err());
    for nodes in [[0,2],[0,99],[1,1]] {
        assert!(mesh.distributed_nodal_forces(&zero,&[MembraneEdgeTraction {nodes,..edge}],[0.0;3],&gate).is_err());
    }
    assert!(mesh.distributed_nodal_forces(&zero[..1],&[],[0.0;3],&gate).is_err());
    assert!(mesh.pressure_nodal_forces(&[[f64::NAN;3];2],&gate).is_err());
    let cancelled=CancelGate::new_clock_free();cancelled.request();
    assert!(mesh.distributed_nodal_forces(&zero,&[],[0.0;3],&cancelled).is_err());
}

#[test]
fn distributed_loading_drives_the_same_history_and_rejects_foreign_bodies_atomically() {
    let mesh=specimen(false);let other=specimen(false);let gate=CancelGate::new_clock_free();
    let nodal=mesh.distributed_nodal_forces(&[[[0.01,0.0,0.0];3]],&[],[0.0;3],&gate).unwrap();
    let free=mesh.project_nodal_forces(&nodal).unwrap();let zeros=[[0.0;3];3];
    let mut a=mesh.history(&zeros,&zeros,InitialMemory::Relaxed,&gate).unwrap();
    let mut b=a.clone();
    for _ in 0..8 {
        let ra=mesh.advance_nodal(&mut a,&nodal,0.002,&gate).unwrap();
        let rb=b.advance(&free,0.002,&gate).unwrap();
        assert_eq!(a.state(),b.state());assert_eq!(ra.supplied.to_bits(),rb.supplied.to_bits());
        assert_eq!(ra.dissipated.to_bits(),rb.dissipated.to_bits());
    }
    assert!(a.displacement_m().iter().any(|u|u.abs()>1e-7));
    assert!(a.dissipated_heat_j()>0.0);close(a.balance_residual_j(),0.0);
    let before=a.state().to_vec();let time=a.time_s();let heat=a.dissipated_heat_j();
    assert!(other.advance_nodal(&mut a,&nodal,0.002,&gate).is_err());
    assert!(mesh.advance_nodal(&mut a,&[[f64::NAN;3];3],0.002,&gate).is_err());
    assert_eq!(a.state(),before);assert_eq!(a.time_s(),time);assert_eq!(a.dissipated_heat_j(),heat);
}
