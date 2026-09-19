use fs_plate::{ModePair, PlateError};
use fs_plate::shell::head::{TensionedDisk, TensionedDiskSpec};
use fs_plate::shell::head::nonlinear::{MembraneReduction, MembraneReductionBudget};
use fs_plate::shell::profile::ProfileBudget;

fn budget() -> MembraneReductionBudget {
    MembraneReductionBudget { max_modes: 8, max_nodes: 100,
        max_facet_pairs: 2000, max_solve_entries: 100_000, relative_tolerance: 1e-6 }
}
fn fixture() -> (TensionedDisk, Vec<ModePair>, Vec<usize>) {
    let disk = TensionedDisk::new(TensionedDiskSpec { radius_m: 0.17,
        thickness_m: 0.000254, young_pa: 4e9, poisson: 0.38,
        density_kg_m3: 1390.0, tension_n_m: 3000.0, radial_intervals: 2, azimuths: 8 },
        ProfileBudget { max_nodes: 100, max_triangles: 100, max_feature_evaluations: 0 }).unwrap();
    let n = disk.model.free;
    let mut k = vec![0.0;n*n]; let mut m = k.clone();
    for i in 0..n { for j in 0..n {
        k[i*n+j] = disk.model.k.get(i,j); m[i*n+j] = disk.model.m.get(i,j);
    }}
    let mut modes = fs_modal::eigh_gen_dense(&k,&m,n).unwrap(); modes.truncate(3);
    let rim = (0..disk.mesh.nodes.len()).filter(|&i| disk.model.dof_map[3*i].is_none()).collect();
    (disk,modes,rim)
}
// Independent test solver: spectral inverse from the existing modal owner,
// rather than duplicating the production host's Cholesky solve.
fn solve(n: usize, k: &[f64], rhs: &[Vec<f64>]) -> Result<Vec<Vec<f64>>,PlateError> {
    let mut identity = vec![0.0;n*n]; for i in 0..n { identity[i*n+i] = 1.0; }
    let eigen = fs_modal::eigh_gen_dense(k,&identity,n).unwrap();
    assert!(eigen.iter().all(|e| e.lambda > 0.0));
    Ok(rhs.iter().map(|b| {
        let mut u = vec![0.0;n];
        for e in &eigen {
            let amplitude = e.phi.iter().zip(b).map(|(p,b)|p*b).sum::<f64>()/e.lambda;
            for (u,p) in u.iter_mut().zip(&e.phi) { *u += amplitude*p; }
        }
        u
    }).collect())
}
fn reduced(d: &TensionedDisk, modes: &[ModePair], rim: &[usize]) -> MembraneReduction {
    MembraneReduction::new(&d.mesh,&d.section,&d.model,modes,rim,budget(),solve).unwrap()
}
#[test]
fn stretching_is_positive_quartic_and_gradient_is_energy_derivative() {
    let (d,modes,rim) = fixture(); let r = reduced(&d,&modes,&rim);
    let q = [3e-4,-1.1e-4,0.7e-4]; let twice = q.map(|v|2.0*v); let minus = q.map(|v|-v);
    let e = r.stretching_energy(&q); assert!(e > 0.0);
    assert!((r.stretching_energy(&twice)-16.0*e).abs() < 1e-12*e);
    assert_eq!(r.potential(&q),r.potential(&minus));
    let mut g = [0.0;3]; r.gradient(&q,&mut g);
    for i in 0..3 {
        let h = 1e-8; let mut a = q; let mut b = q; a[i] += h; b[i] -= h;
        let fd = (r.potential(&a)-r.potential(&b))/(2.0*h);
        assert!((g[i]-fd).abs() < 2e-6*g[i].abs().max(fd.abs()).max(1e-8));
    }
    let tiny = [1e-10,0.0,0.0]; r.gradient(&tiny,&mut g);
    assert!((g[0]/tiny[0]-modes[0].lambda).abs() < 1e-6*modes[0].lambda);
    r.gradient(&q,&mut g);
    let linear_force = modes.iter().zip(q).map(|(m,q)|m.lambda*q*q).sum::<f64>();
    let radial_force = q.iter().zip(g).map(|(q,g)|q*g).sum::<f64>();
    assert!((radial_force-linear_force-4.0*e).abs() < 1e-6*linear_force);
    assert!(r.maximum_slope(&q) > 0.0 && r.maximum_slope(&q) < 1.0);
    assert!(r.solve_residual() < budget().relative_tolerance);
}
#[test]
fn interior_relaxation_is_not_artificially_locked() {
    let (d,modes,rim) = fixture(); let relaxed = reduced(&d,&modes,&rim);
    let all: Vec<_> = (0..d.mesh.nodes.len()).collect();
    let locked = MembraneReduction::new(&d.mesh,&d.section,&d.model,&modes,&all,budget(),
        |_,_,_| panic!("all-fixed geometry must not call a linear solver")).unwrap();
    let q = [3e-4,-1.1e-4,0.7e-4];
    assert!(relaxed.stretching_energy(&q) < 0.999*locked.stretching_energy(&q));
    assert!(relaxed.stretching_energy(&q) > 0.0);
}
#[test]
fn condensed_energy_matches_independent_nodal_membrane_reconstruction() {
    let (d,modes,rim) = fixture(); let q = [3e-4,-1.1e-4,0.7e-4];
    let mut columns = Vec::new();
    let r = MembraneReduction::new(&d.mesh,&d.section,&d.model,&modes,&rim,budget(),|n,k,rhs| {
        let u = solve(n,k,rhs)?; columns = u.clone(); Ok(u)
    }).unwrap();
    let mut free = 0; let mut u = vec![[0.0;2];d.mesh.nodes.len()];
    for (node,uv) in u.iter_mut().enumerate() { if !rim.contains(&node) {
        let mut pair = 0;
        for i in 0..3 { for j in i..3 {
            for c in 0..2 { uv[c] += columns[pair][free+c]*q[i]*q[j]; } pair += 1;
        }}
        free += 2;
    }}
    let w: Vec<_> = (0..d.mesh.nodes.len()).map(|node| d.model.dof_map[3*node]
        .map_or(0.0, |p| modes.iter().zip(q).map(|(m,q)|m.phi[p]*q).sum::<f64>())).collect();
    let mut energy = 0.0;
    for tri in &d.mesh.tris {
        let (x0,y0)=d.mesh.nodes[tri[0]]; let (x1,y1)=d.mesh.nodes[tri[1]]; let (x2,y2)=d.mesh.nodes[tri[2]];
        let a2=(x1-x0)*(y2-y0)-(x2-x0)*(y1-y0);
        let b=[(y1-y2)/a2,(y2-y0)/a2,(y0-y1)/a2];
        let c=[(x2-x1)/a2,(x0-x2)/a2,(x1-x0)/a2];
        let mut wx=0.0; let mut wy=0.0; let mut strain=[0.0;3];
        for i in 0..3 {
            wx+=b[i]*w[tri[i]]; wy+=c[i]*w[tri[i]];
            strain[0]+=b[i]*u[tri[i]][0]; strain[1]+=c[i]*u[tri[i]][1];
            strain[2]+=c[i]*u[tri[i]][0]+b[i]*u[tri[i]][1];
        }
        strain[0]+=0.5*wx*wx; strain[1]+=0.5*wy*wy; strain[2]+=wx*wy;
        for i in 0..3 { for j in 0..3 {
            energy+=0.25*a2*12.0/d.section.thickness.powi(2)*strain[i]*d.section.d[3*i+j]*strain[j];
        }}
    }
    assert!((energy-r.stretching_energy(&q)).abs() < 2e-12*energy);
}
#[test]
fn refuses_wrong_basis_budget_supports_and_fake_static_solutions() {
    let (d,mut modes,rim) = fixture();
    let mut cap = budget(); cap.max_solve_entries = 1;
    assert!(MembraneReduction::new(&d.mesh,&d.section,&d.model,&modes,&rim,cap,
        |_,_,_| panic!("budget admission must precede the solve")).is_err());
    assert!(MembraneReduction::new(&d.mesh,&d.section,&d.model,&modes,&rim,budget(),
        |n,_,rhs| Ok(vec![vec![0.0;n];rhs.len()])).is_err());
    assert!(MembraneReduction::new(&d.mesh,&d.section,&d.model,&modes,&[0,0],budget(),solve).is_err());
    modes[0].lambda *= 2.0;
    assert!(MembraneReduction::new(&d.mesh,&d.section,&d.model,&modes,&rim,budget(),solve).is_err());
    modes[0].lambda *= 0.5;
    let r = reduced(&d,&modes,&rim); let mut out = [0.0;2]; r.gradient(&[0.0;3],&mut out);
    assert!(out.iter().all(|v|v.is_nan())); assert!(r.potential(&[f64::NAN;3]).is_nan());
    assert!(r.maximum_slope(&[]).is_nan());
}
