//! Initial stress reaches the ACTUAL consistent-mass vibration problem.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::ContentHash;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_modal::eigh_gen_dense;
use fs_solid::linear3::*;
use fs_solid::prestress::PrestressedTetAssembly;

const NODES: [[f64; 3]; 4] = [[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]];
const TETS: [[usize; 4]; 1] = [[0,1,2,3]];
// Prevent axial expansion but allow lateral Poisson response and transverse modes.
const FIXED: [usize; 7] = [0,1,2,3,6,8,9];
fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| f(&Cx::new(&gate, arena,
        StreamKey { seed:3, kernel_id:17, tile:0, iteration:0 },
        Budget::INFINITE, ExecMode::Deterministic)))
}
fn material() -> TetElasticMaterial { TetElasticMaterial::try_new(6.,1000.,0.25,[7;32]).unwrap() }
fn strain(a: f64) -> TetThermalStrainState {
    TetThermalStrainState::try_new(300.+a*1000.,300.,[a,a,a,0.,0.,0.],
        [7;32],ContentHash([9;32]),0.1).unwrap()
}
fn problem<'a>(mat: &'a TetElasticMaterial, tets: &'a [[usize;4]], fixed: &'a [usize]) -> TetLinearElasticProblem<'a> {
    TetLinearElasticProblem { nodes_m:&NODES, tetrahedra:tets,
        materials:TetMaterialField::Uniform(mat), fixed_dofs:fixed, budget:TetAssemblyBudget::standard() }
}
fn assembled(a: f64, cx: &Cx<'_>) -> PrestressedTetAssembly {
    let mat=material(); let thermal=strain(a); let p=problem(&mat,&TETS,&FIXED);
    let field=TetThermalStrainField::Uniform(&thermal);
    let u=p.solve_thermal_displacement(field,TetStaticSolveConfig::default(),cx).unwrap();
    p.assemble_thermal_tangent(&u,field,ContentHash([1;32]),0.1,cx).unwrap()
}
fn close(a:f64,b:f64) { assert!((a-b).abs()<2e-10*a.abs().max(b.abs()).max(1.),"{a} != {b}"); }
fn eigen(p: &PrestressedTetAssembly) -> Vec<f64> {
    let a=p.pencil();
    eigh_gen_dense(&a.stiffness.to_dense(),&a.mass.to_dense(),a.free_dofs.len()).unwrap()
        .into_iter().map(|m|m.lambda).collect()
}

#[test]
fn zero_stress_preserves_the_original_pencil_and_repeat_identity() {
    with_cx(|cx| {
        let a=assembled(0.,cx); let mat=material(); let b=problem(&mat,&TETS,&FIXED).assemble(cx).unwrap();
        assert_eq!(a.pencil().stiffness.to_dense(),b.stiffness.to_dense());
        assert_eq!(a.pencil().mass.to_dense(),b.mass.to_dense());
        assert_eq!(a.pencil().free_dofs,b.free_dofs);
        assert!(a.geometric_stiffness().to_dense().iter().all(|&v|v==0.));
        assert_eq!(a.identity(),assembled(0.,cx).identity());
    });
}

#[test]
fn thermal_restraint_stress_changes_modes_with_the_correct_sign_and_unchanged_mass() {
    with_cx(|cx| {
        let cold=assembled(-0.01,cx); let base=assembled(0.,cx); let hot=assembled(0.01,cx);
        for (a,p) in [(-0.01,&cold),(0.01,&hot)] {
            let s=p.stress().elements()[0].stress_mandel_pa();
            close(s[0],-1000.*a);
            for &v in &s[1..] {close(v,0.);}
            let i=p.pencil().free_dofs.iter().position(|&d|d==4).unwrap();
            close(p.geometric_stiffness().get(i,i),-1000.*a/6.);
            assert_eq!(p.pencil().mass.to_dense(),base.pencil().mass.to_dense());
            assert_eq!(p.pencil().total_mass_kg,1.);
        }
        let (a,b,c)=(eigen(&hot),eigen(&base),eigen(&cold));
        assert!(a.iter().zip(&b).any(|(a,b)|*a < *b-1.));
        for ((a,b),c) in a.iter().zip(b).zip(c) {assert!(*a<=b+1e-8 && b<=c+1e-8);}
        assert_ne!(hot.identity(),cold.identity());
    });
}

#[test]
fn geometric_energy_matches_independent_affine_gradient_contraction() {
    with_cx(|cx| {
        let p=assembled(0.01,cx); let a=p.pencil();
        let full=[[0.3,-0.1,0.7],[0.2,0.9,-0.4],[-0.3,0.2,0.6],[0.1,-0.6,0.8]];
        let mut admitted=[[0.;3];4]; let mut q=vec![0.;a.free_dofs.len()];
        for (i,&d) in a.free_dofs.iter().enumerate() {q[i]=full[d/3][d%3]; admitted[d/3][d%3]=q[i];}
        let mut gq=vec![0.;q.len()]; p.geometric_stiffness().spmv(&q,&mut gq);
        let actual:f64=q.iter().zip(gq).map(|(x,y)|x*y).sum();
        // For the unit reference tet, du/dx = u(node 1)-u(node 0).
        // sigma_xx=-10 Pa, volume=1/6 m3, all other stress components zero.
        let expected=-10./6.*(0..3).map(|i|(admitted[1][i]-admitted[0][i]).powi(2)).sum::<f64>();
        close(actual,expected);
        for i in 0..q.len() {for j in 0..q.len() {
            assert_eq!(p.geometric_stiffness().get(i,j).to_bits(),p.geometric_stiffness().get(j,i).to_bits());
        }}
    });
}

#[test]
fn element_orientation_does_not_change_the_physical_operator() {
    with_cx(|cx| {
        let a=assembled(0.01,cx); let mat=material(); let thermal=strain(0.01);
        for tet in [[0,2,1,3],[3,2,0,1]] {
            let tets=[tet]; let p=problem(&mat,&tets,&FIXED); let field=TetThermalStrainField::Uniform(&thermal);
            let u=p.solve_thermal_displacement(field,TetStaticSolveConfig::default(),cx).unwrap();
            let b=p.assemble_thermal_tangent(&u,field,ContentHash([1;32]),0.1,cx).unwrap();
            for (x,y) in a.pencil().stiffness.to_dense().iter().zip(b.pencil().stiffness.to_dense()) {close(*x,y);}
            for (x,y) in a.pencil().mass.to_dense().iter().zip(b.pencil().mass.to_dense()) {close(*x,y);}
        }
    });
}

#[test]
fn stale_material_thermal_or_support_inputs_cannot_reuse_a_solution() {
    with_cx(|cx| {
        let mat=material(); let t=strain(0.01); let p=problem(&mat,&TETS,&FIXED);
        let field=TetThermalStrainField::Uniform(&t);
        let u=p.solve_thermal_displacement(field,TetStaticSolveConfig::default(),cx).unwrap();
        let changed=strain(0.02);
        assert!(matches!(p.assemble_thermal_tangent(&u,TetThermalStrainField::Uniform(&changed),ContentHash([1;32]),0.1,cx),
            Err(TetElasticError::StressRecoveryInputMismatch)));
        let fixed=[0,1,2,3,4,6,8,9];
        assert!(problem(&mat,&TETS,&fixed).assemble_thermal_tangent(&u,field,ContentHash([1;32]),0.1,cx).is_err());
        let other=TetElasticMaterial::try_new(6.,2000.,0.25,[7;32]).unwrap();
        assert!(problem(&other,&TETS,&FIXED).assemble_thermal_tangent(&u,field,ContentHash([1;32]),0.1,cx).is_err());
        assert!(p.assemble_thermal_tangent(&u,field,ContentHash([1;32]),0.001,cx).is_err());
    });
}

#[test]
fn all_fixed_and_cancelled_requests_do_not_publish_a_pencil() {
    let (mat,t)=(material(),strain(0.01)); let fixed:Vec<_>=(0..12).collect();
    with_cx(|cx| {
        let p=problem(&mat,&TETS,&fixed); let field=TetThermalStrainField::Uniform(&t);
        let u=p.solve_thermal_displacement(field,TetStaticSolveConfig::default(),cx).unwrap();
        assert!(matches!(p.assemble_thermal_tangent(&u,field,ContentHash([1;32]),0.1,cx),Err(TetElasticError::NoFreeDofs)));
    });
    let p=problem(&mat,&TETS,&FIXED); let field=TetThermalStrainField::Uniform(&t);
    let u=with_cx(|cx|p.solve_thermal_displacement(field,TetStaticSolveConfig::default(),cx)).unwrap();
    let gate=CancelGate::new(); gate.request();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx=Cx::new(&gate,arena,StreamKey {seed:3,kernel_id:17,tile:0,iteration:0},Budget::INFINITE,ExecMode::Deterministic);
        assert!(matches!(p.assemble_thermal_tangent(&u,field,ContentHash([1;32]),0.1,&cx),Err(TetElasticError::Cancelled)));
    });
}


#[test]
fn thermal_and_dead_nodal_load_share_one_equilibrium_and_stress_identity() {
    with_cx(|cx| {
        let mat=material(); let thermal=strain(0.01); let p=problem(&mat,&TETS,&FIXED);
        let field=TetThermalStrainField::Uniform(&thermal);
        let forces=[(4usize,0.25)];
        let u=p.solve_thermal_displacement_with_nodal_forces(
            field,&forces,TetStaticSolveConfig::default(),cx).unwrap();
        assert!(u.displacement_m()[1][1] > 0.0);
        let stress=p.recover_thermal_stress_with_nodal_forces(
            &u,field,&forces,ContentHash([1;32]),0.1,cx).unwrap();
        close(stress.nodal_internal_force_n()[1][1],0.25);
        assert!(matches!(
            p.recover_thermal_stress_with_nodal_forces(
                &u,field,&[(4,0.5)],ContentHash([1;32]),0.1,cx),
            Err(TetElasticError::StressRecoveryInputMismatch)
        ));
    });
}

#[test]
fn fixed_dof_load_preserves_displacement_and_malformed_loads_refuse() {
    with_cx(|cx| {
        let mat=material(); let thermal=strain(0.01); let p=problem(&mat,&TETS,&FIXED);
        let field=TetThermalStrainField::Uniform(&thermal);
        let base=p.solve_thermal_displacement(field,TetStaticSolveConfig::default(),cx).unwrap();
        let loaded=p.solve_thermal_displacement_with_nodal_forces(
            field,&[(0,0.5)],TetStaticSolveConfig::default(),cx).unwrap();
        assert_eq!(base.displacement_m(),loaded.displacement_m());
        for bad in [&[(4usize,1.0),(4,2.0)][..],&[(5usize,1.0),(4,2.0)][..]] {
            assert!(matches!(
                p.solve_thermal_displacement_with_nodal_forces(
                    field,bad,TetStaticSolveConfig::default(),cx),
                Err(TetElasticError::InvalidNodalForce{..})
            ));
        }
        assert!(matches!(
            p.solve_thermal_displacement_with_nodal_forces(
                field,&[(99,f64::NAN)],TetStaticSolveConfig::default(),cx),
            Err(TetElasticError::InvalidNodalForce{..})
        ));
    });
}

#[path = "thermal_prestress/column.rs"]
mod column;

#[path = "thermal_prestress/constitutive.rs"]
mod constitutive;
