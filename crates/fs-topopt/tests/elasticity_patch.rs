//! G0/G3 physical patch tests, independent of the assembled stiffness oracle.
use fs_solver::LinearOp;
use fs_topopt::DensityElasticity;

fn close(a: f64, b: f64, scale: f64) {
    assert!((a-b).abs() < 2e-11*scale.max(b.abs()), "{a:e} != {b:e}");
}

#[test]
fn g0_all_three_simple_shears_use_mu_not_twice_mu() {
    let (mesh, positions) = fs_feec::kuhn_cube(2);
    let young = 30.0; let nu = 0.25; let gamma = 0.03;
    let mu = young/(2.0*(1.0+nu));
    let el = DensityElasticity::new(&mesh, &positions, young, nu, &|_| false);
    for (component, axis, shear) in [(0,1,3), (1,2,4), (2,0,5)] {
        let mut u = vec![0.0; el.n()];
        for (i,p) in positions.iter().enumerate() { u[3*i+component] = gamma*p[axis]; }
        for stress in el.cell_stress(&u) {
            for (i, value) in stress.into_iter().enumerate() {
                close(value, if i==shear { mu*gamma } else { 0.0 }, 1.0);
            }
        }
        // Unit cube; cell_energies stores u_e^T K_e u_e, twice stored energy.
        close(el.cell_energies(&u).iter().sum(), mu*gamma*gamma, 1.0);
        let mut ku = vec![0.0; el.n()]; el.apply(&u, &mut ku);
        close(u.iter().zip(ku).map(|(u,f)|u*f).sum(), mu*gamma*gamma, 1.0);
    }
}

#[test]
fn g0_affine_stress_and_nodal_work_match_continuum_hooke_law() {
    let (mesh, mut positions) = fs_feec::kuhn_cube(2);
    for p in &mut positions { let [x,y,z]=*p; *p=[1.7*x+0.2*y,0.9*y+0.1*z,1.3*z]; }
    let geo = fs_feec::element_geometry(&mesh,&positions);
    let young = 37.0; let nu = 0.27;
    let mu = young/(2.0*(1.0+nu));
    let lambda = young*nu/((1.0+nu)*(1.0-2.0*nu));
    let a = [[0.01,0.02,-0.03],[-0.04,0.05,0.01],[0.02,-0.01,-0.02]];
    let strain: [[f64;3];3] = core::array::from_fn(|i|core::array::from_fn(|j|0.5*(a[i][j]+a[j][i])));
    let trace = strain[0][0]+strain[1][1]+strain[2][2];
    let stress: [[f64;3];3] = core::array::from_fn(|i|core::array::from_fn(|j|
        2.0*mu*strain[i][j]+if i==j {lambda*trace} else {0.0}));
    let sigma = [stress[0][0],stress[1][1],stress[2][2],stress[0][1],stress[1][2],stress[2][0]];
    let el = DensityElasticity::new(&mesh,&positions,young,nu,&|_|false);
    let u: Vec<_> = positions.iter().flat_map(|p|(0..3).map(move|i|
        a[i][0]*p[0]+a[i][1]*p[1]+a[i][2]*p[2])).collect();
    let mut expected_force = vec![0.0;el.n()];
    for (e,tet) in mesh.tets.iter().enumerate() {
        for (a,&node) in tet.iter().enumerate() { for i in 0..3 {
            expected_force[3*node as usize+i] += geo.vol_signed[e].abs()
                *(0..3).map(|j|stress[i][j]*geo.grads[e][a][j]).sum::<f64>();
        }}
    }
    let mut ku = vec![0.0;el.n()]; el.apply(&u,&mut ku);
    for (&got,&expected) in ku.iter().zip(&expected_force) {close(got,expected,1.0);}
    for got in el.cell_stress(&u) { for (&got,&expected) in got.iter().zip(&sigma) {close(got,expected,1.0);} }
    let density: f64 = (0..3).flat_map(|i|(0..3).map(move|j|strain[i][j]*stress[i][j])).sum();
    close(el.cell_energies(&u).iter().sum(),density*geo.vol_signed.iter().map(|v|v.abs()).sum::<f64>(),1.0);
}

#[test]
fn g3_rigid_motion_has_no_stress_and_mass_keeps_its_original_units() {
    let (mesh, positions) = fs_feec::kuhn_cube(2);
    let el = DensityElasticity::new(&mesh,&positions,30.0,0.25,&|_|false);
    let u: Vec<_> = positions.iter().flat_map(|p|[1.0-p[1],2.0+p[0],-0.5]).collect();
    for s in el.cell_stress(&u) {for x in s {close(x,0.0,1.0);}}
    close(el.cell_energies(&u).iter().sum(),0.0,1.0);
    let t: Vec<_> = positions.iter().flat_map(|_|[1.0,2.0,-0.5]).collect();
    let mut mt = vec![0.0;el.n()]; el.apply_mass(&vec![7.0;el.cells()],&t,&mut mt);
    close(t.iter().zip(mt).map(|(u,f)|u*f).sum(),7.0*5.25,1.0);
}
