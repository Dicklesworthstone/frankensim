use super::*;
use fs_material::{IsotropicElastic, visco_tensor::TensorPronyBranch};

fn projection() -> TensorProjectionBudget {
    TensorProjectionBudget { max_coordinates: 128, max_scalar_terms: 100_000_000 }
}
fn dynamics() -> TensorMaterialBudget {
    TensorMaterialBudget { max_state_coordinates: 128, max_material_points: 128 }
}
fn material(memory: bool, poisson: f64) -> MembraneFacetMaterial {
    let c = IsotropicElastic::new(1000.0, poisson, 0.2).unwrap().stiffness();
    let branches = if memory {
        vec![TensorPronyBranch::from_stiffness_mandel(c.map(|r| r.map(|v| 0.4*v)), 0.02).unwrap()]
    } else { vec![] };
    MembraneFacetMaterial {
        law: TensorProny::from_tensor_shear(c, branches, 0.2).unwrap(),
        thickness_m: 0.1, density_kg_m3: 10.0,
    }
}
fn triangle() -> ShellMesh {
    ShellMesh::new(vec![[0.0,0.0,0.0],[1.0,0.0,0.0],[0.0,1.0,0.0]], vec![[0,1,2]]).unwrap()
}
const FIXED: [usize;6] = [0,1,2,4,5,8];
fn close(a: f64, b: f64, tol: f64) {
    assert!((a-b).abs() <= tol*(1.0+a.abs().max(b.abs())), "{a:e} != {b:e}");
}

#[test]
fn affine_patch_has_original_plane_stress_energy_mass_and_cartesian_work() {
    let gate = CancelGate::new_clock_free();
    for memory in [false,true] {
        let mesh = TensorMembrane::new(triangle(), vec![material(memory,0.25)], &FIXED,
            projection(), dynamics(), &gate).unwrap();
        assert_eq!(mesh.body().coordinates(),3);
        assert_eq!(mesh.body().state_dimension(), if memory {12} else {6});
        assert_eq!(mesh.dof_map(), &[None,None,None,Some(0),None,None,Some(1),Some(2),None]);
        for mass in mesh.nodal_masses_kg() { close(*mass, 1.0/6.0, 1e-14); }
        let u = [[0.0;3],[0.001,0.0,0.0],[0.0006,-0.0002,0.0]];
        let relaxed = mesh.history(&u, &[[0.0;3];3], InitialMemory::Relaxed, &gate).unwrap();
        let strain = relaxed.plane_stress_strain(0).unwrap();
        for (got, expected) in strain.into_iter().zip([0.001,-0.0002,-0.0008/3.0,0.0003,0.0,0.0]) {
            close(got,expected,1e-12);
        }
        let energy = 0.5*0.05*(1000.0/0.9375*(1e-6+4e-8-1e-7)+400.0*3.6e-7);
        close(relaxed.stored_energy_j(),energy,1e-12);
        let instant = mesh.history(&u, &[[0.0;3];3], InitialMemory::Instantaneous, &gate).unwrap();
        close(instant.stored_energy_j(),energy*if memory {1.4} else {1.0},1e-12);
        for s in instant.plane_stress_normal_stress_pa(0).unwrap() { close(s,0.0,1e-10); }
        assert_eq!(mesh.expand_coordinates(relaxed.displacement_m()).unwrap(),u);
        let force = [[19.0,-8.0,7.0],[2.0,3.0,4.0],[5.0,6.0,7.0]];
        let f = mesh.project_nodal_forces(&force).unwrap();
        let nodal_work: f64 = u.iter().flatten().zip(force.iter().flatten()).map(|(a,b)|a*b).sum();
        let free_work: f64 = relaxed.displacement_m().iter().zip(f).map(|(a,b)|a*b).sum();
        close(nodal_work,free_work,1e-15);
    }
}

#[test]
fn oriented_membrane_motion_and_heat_follow_the_original_geometry() {
    let gate = CancelGate::new_clock_free();
    let rot = |v:[f64;3]| [-v[1],v[0],v[2]];
    let mut a = material(true,0.25);
    // A directional, non-proportional branch makes accidental global/local
    // mixing visible; material coefficients are in the actual facet frame.
    let mut gram = [[0.0;6];6]; gram[0]=[8.0,2.0,1.0,3.0,0.0,0.0];
    a.law = TensorProny::new_mandel(*a.law.equilibrium_mandel(),
        vec![TensorPronyBranch::from_gram_rows(gram,0.02).unwrap()],0.2).unwrap();
    let mesh = triangle();
    let rotated = ShellMesh::new(mesh.nodes.iter().copied().map(rot).collect(),mesh.tris.clone()).unwrap();
    let original = TensorMembrane::new(mesh,vec![a.clone()],&FIXED,projection(),dynamics(),&gate).unwrap();
    let turned = TensorMembrane::new(rotated,vec![a],&[0,1,2,3,5,8],projection(),dynamics(),&gate).unwrap();
    let zeros = [[0.0;3];3];
    let mut x = original.history(&zeros,&zeros,InitialMemory::Relaxed,&gate).unwrap();
    let mut y = turned.history(&zeros,&zeros,InitialMemory::Relaxed,&gate).unwrap();
    let load = [[0.0;3],[0.01,0.0,0.0],[0.003,-0.002,0.0]];
    let f = original.project_nodal_forces(&load).unwrap();
    let g = turned.project_nodal_forces(&load.map(rot)).unwrap();
    let unloaded = [0.0;3];
    for step in 0..12 {
        x.advance(if step<4 {&f} else {&unloaded},0.002,&gate).unwrap();
        y.advance(if step<4 {&g} else {&unloaded},0.002,&gate).unwrap();
        let a = original.expand_coordinates(x.displacement_m()).unwrap();
        let b = turned.expand_coordinates(y.displacement_m()).unwrap();
        for (a,b) in a.into_iter().map(rot).zip(b) { for (a,b) in a.into_iter().zip(b) {close(a,b,1e-10);} }
        close(x.stored_energy_j(),y.stored_energy_j(),1e-10);
        close(x.dissipated_heat_j(),y.dissipated_heat_j(),1e-10);
        close(x.balance_residual_j(),0.0,1e-10);
        for p in y.plane_stress_normal_stress_pa(0).unwrap() { close(p,0.0,1e-10); }
    }
    assert!(x.dissipated_heat_j()>0.0);
    assert!(x.displacement_m().iter().any(|x|x.abs()>1e-7));
    let saved = x.state().to_vec(); let time=x.time_s(); let work=x.supplied_work_j();
    let cancelled=CancelGate::new_clock_free(); cancelled.request();
    assert!(x.advance(&f,0.002,&cancelled).is_err());
    assert_eq!(x.state(),saved); assert_eq!(x.time_s(),time); assert_eq!(x.supplied_work_j(),work);
}

fn strip(cells: usize) -> (ShellMesh, Vec<usize>) {
    let nodes: Vec<_> = (0..=cells).flat_map(|i| [[i as f64/cells as f64,0.0,0.0],
        [i as f64/cells as f64,0.1,0.0]]).collect();
    let tris = (0..cells).flat_map(|i| [[2*i,2*i+2,2*i+3],[2*i,2*i+3,2*i+1]]).collect();
    let mut fixed: Vec<_> = (0..nodes.len()).flat_map(|i| [3*i+1,3*i+2]).collect();
    fixed.extend([0,3]);
    (ShellMesh::new(nodes,tris).unwrap(),fixed)
}

#[test]
fn three_mesh_rungs_resolve_the_axial_strip_frequency_without_retuning() {
    let gate = CancelGate::new_clock_free();
    let exact = core::f64::consts::FRAC_PI_2*fs_math::det::sqrt(1000.0/10.0);
    let mut previous = 1.0;
    for cells in [1,2,4] {
        let (geometry,fixed)=strip(cells);
        let facets=geometry.tris.len();
        let mesh=TensorMembrane::new(geometry,vec![material(false,0.0);facets],&fixed,
            projection(),dynamics(),&gate).unwrap();
        let n=mesh.body().coordinates(); assert_eq!(n,2*cells);
        let mass:Vec<_>=mesh.dof_map().iter().enumerate().filter_map(|(i,d)|d.map(|_|mesh.nodal_masses_kg()[i/3])).collect();
        let mut operator=vec![0.0;n*n];
        for point in mesh.body().plane_stress_points() {
            for i in 0..n {for j in 0..n {
                operator[i*n+j]+=point.energy_matrix()[i*n+j]/fs_math::det::sqrt(mass[i]*mass[j]);
            }}
        }
        let (eigenvalues,_)=fs_la::eigen::jacobi_eigh(&operator,n);
        assert!(eigenvalues.iter().all(|v|*v>0.0));
        let lowest=eigenvalues.into_iter().fold(f64::INFINITY,f64::min);
        let error=(fs_math::det::sqrt(lowest)/exact-1.0).abs();
        assert!(error<0.4*previous,"cells={cells}: error={error}, previous={previous}");
        if cells==4 {assert!(error<0.012);}
        previous=error;
        close(mesh.nodal_masses_kg().iter().sum(),0.1,1e-13);
    }
}

#[test]
fn incomplete_geometry_constraints_materials_and_budgets_refuse_before_publication() {
    let gate=CancelGate::new_clock_free();
    let build=|mesh,materials,fixed:&[usize],p,d|TensorMembrane::new(mesh,materials,fixed,p,d,&gate);
    assert!(build(triangle(),vec![],&FIXED,projection(),dynamics()).is_err());
    assert!(build(triangle(),vec![material(false,0.25)],&[9],projection(),dynamics()).is_err());
    assert!(build(triangle(),vec![material(false,0.25)],&[0,0],projection(),dynamics()).is_err());
    assert!(build(triangle(),vec![material(true,0.25)],&FIXED,projection(),
        TensorMaterialBudget {max_state_coordinates:11,..dynamics()}).is_err());
    assert!(build(triangle(),vec![material(true,0.25)],&FIXED,
        TensorProjectionBudget {max_scalar_terms:1,..projection()},dynamics()).is_err());
    let mut bad_mesh=triangle();bad_mesh.tris[0][2]=99;
    assert!(matches!(build(bad_mesh,vec![material(false,0.25)],&FIXED,projection(),dynamics()),Err(MembraneError::Geometry(_))));
    let mut isolated=triangle();isolated.nodes.push([2.0,2.0,0.0]);
    assert!(build(isolated,vec![material(false,0.25)],&FIXED,projection(),dynamics()).is_err());
    let mut duplicate=triangle();duplicate.tris.push([2,1,0]);
    assert!(build(duplicate,vec![material(false,0.25);2],&FIXED,projection(),dynamics()).is_err());
    for value in [0.0,-1.0,f64::INFINITY,f64::NAN] {
        let mut m=material(false,0.25);m.thickness_m=value;
        assert!(build(triangle(),vec![m],&FIXED,projection(),dynamics()).is_err());
    }
    let mesh=build(triangle(),vec![material(false,0.25)],&FIXED,projection(),dynamics()).unwrap();
    let mut state=[[0.0;3];3];state[0][0]=0.1;
    assert!(mesh.history(&state,&[[0.0;3];3],InitialMemory::Relaxed,&gate).is_err());
    assert!(mesh.project_nodal_forces(&[[f64::NAN;3];3]).is_err());
    let cancelled=CancelGate::new_clock_free();cancelled.request();
    assert!(matches!(TensorMembrane::new(triangle(),vec![material(false,0.25)],&FIXED,
        projection(),dynamics(),&cancelled),Err(MembraneError::Dynamics(TensorMaterialError::Cancelled))));
}


#[test]
fn consistent_membrane_mass_matches_velocity_quadrature_and_retains_material_ringdown() {
    let gate=CancelGate::new_clock_free();
    let mesh=TensorMembrane::new_with_consistent_mass(triangle(),vec![material(true,0.25)],
        &FIXED,projection(),dynamics(),&gate).unwrap();
    let mass=mesh.body().mass_matrix_kg();
    let expected=[1.0/12.0,1.0/24.0,0.0,1.0/24.0,1.0/12.0,0.0,0.0,0.0,1.0/12.0];
    for (got,want) in mass.iter().zip(expected) {close(*got,want,1e-15);}
    let velocity=[[0.0;3],[0.003,0.0,0.0],[-0.002,0.004,0.0]];
    let mut history=mesh.history(&[[0.0;3];3],&velocity,InitialMemory::Relaxed,&gate).unwrap();
    // Three barycentric quadrature points integrate the squared P1 velocity
    // exactly, independently of the closed-form element mass in the constructor.
    let mut kinetic=0.0;
    for shape in [[2.0/3.0,1.0/6.0,1.0/6.0],[1.0/6.0,2.0/3.0,1.0/6.0],
        [1.0/6.0,1.0/6.0,2.0/3.0]] {
        let v:[f64;3]=core::array::from_fn(|axis|(0..3).map(|node|shape[node]*velocity[node][axis]).sum());
        kinetic+=0.5*(0.5/3.0)*v.iter().map(|v|v*v).sum::<f64>();
    }
    close(history.stored_energy_j(),kinetic,1e-15);
    let free=mesh.project_nodal_forces(&velocity).unwrap();
    for i in 0..3 {close(history.state()[3+i],(0..3).map(|j|mass[3*i+j]*free[j]).sum(),1e-15);}
    // Full row sums still describe physical nodal load weights, including the
    // fixed vertex; they are not substituted for the reduced consistent mass.
    for row_sum in mesh.nodal_masses_kg() {close(*row_sum,1.0/6.0,1e-15);}
    for tick in 0..16 {
        let force=if tick<4 {[0.01,-0.003,0.002]} else {[0.0;3]};
        history.advance(&force,0.002,&gate).unwrap();
        close(history.balance_residual_j(),0.0,1e-11);
        for traction in history.plane_stress_normal_stress_pa(0).unwrap() {close(traction,0.0,1e-10);}
    }
    assert!(history.dissipated_heat_j()>0.0);
    assert!(history.displacement_m().iter().any(|v|v.abs()>1e-6));
}

#[test]
fn consistent_inertia_strip_modes_converge_to_the_continuum_frequency() {
    let gate=CancelGate::new_clock_free();
    let exact=core::f64::consts::FRAC_PI_2*fs_math::det::sqrt(1000.0/10.0);
    let mut previous=1.0;
    for cells in [1,2,4] {
        let (geometry,fixed)=strip(cells);let facets=geometry.tris.len();
        let mesh=TensorMembrane::new_with_consistent_mass(geometry,vec![material(false,0.0);facets],
            &fixed,projection(),dynamics(),&gate).unwrap();
        let n=mesh.body().coordinates();let mass=mesh.body().mass_matrix_kg();
        let factor=fs_la::factor::cholesky(&mass,n).unwrap();
        // Generalized eigenvalues from M^-1/2 K M^-1/2 using an independent
        // symmetric eigendecomposition of the original physical mass matrix.
        let (values,vectors)=fs_la::eigen::jacobi_eigh(&mass,n);
        assert!(values.iter().all(|v|*v>0.0));
        let mut inverse_sqrt=vec![0.0;n*n];
        for i in 0..n {for j in 0..n {inverse_sqrt[i*n+j]=(0..n).map(|k|
            vectors[i*n+k]*vectors[j*n+k]/fs_math::det::sqrt(values[k])).sum();}}
        let mut stiffness=vec![0.0;n*n];
        for point in mesh.body().plane_stress_points() {
            for (a,b) in stiffness.iter_mut().zip(point.energy_matrix()) {*a+=b;}
        }
        let mut operator=vec![0.0;n*n];
        for i in 0..n {for j in 0..n {for a in 0..n {for b in 0..n {
            operator[i*n+j]+=inverse_sqrt[i*n+a]*stiffness[a*n+b]*inverse_sqrt[b*n+j];
        }}}}
        let (eigenvalues,modes)=fs_la::eigen::jacobi_eigh(&operator,n);
        let (index,&lowest)=eigenvalues.iter().enumerate().min_by(|a,b|a.1.total_cmp(b.1)).unwrap();
        let omega=fs_math::det::sqrt(lowest);let error=(omega/exact-1.0).abs();
        assert!(error<0.4*previous,"cells={cells}: error={error}, previous={previous}");
        if cells==4 {assert!(error<0.012);}
        previous=error;
        // Check the generalized eigenpair against untransformed M and K too.
        let mode:Vec<_>=(0..n).map(|i|(0..n).map(|k|inverse_sqrt[i*n+k]*modes[k*n+index]).sum::<f64>()).collect();
        let mut acceleration:Vec<_>=(0..n).map(|i|(0..n).map(|j|stiffness[i*n+j]*mode[j]).sum::<f64>()).collect();
        factor.solve(&mut acceleration);
        for (actual,q) in acceleration.into_iter().zip(mode) {close(actual,lowest*q,1e-9);}
    }
}
