#![cfg(feature = "certified-speculation")]
//! Actual affine forcing must not be certified as its cell-average substitute.
use fs_verify::interval::Iv;
use fs_verify::tet::{AffineSourceTetProblem, BoundaryCondition as Bc, BoundaryFace,
    ConductivityTensor, FluxBudget, TensorTetProblem, TetError, affine_source_energy_bound,
    affine_source_goal_bound, affine_source_mean_bound, tensor_energy_bound, tensor_mean_bound};
use std::collections::BTreeMap;

const VERTICES: [[f64;3];4] = [[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]];
const TETS: [[usize;4];1] = [[0,1,2,3]];
const IDENTITY: ConductivityTensor = [[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]];
fn contains(x: Iv, truth: f64) { assert!(x.lo <= truth && truth <= x.hi, "{x:?} misses {truth}"); }
fn fixed_tet() -> Vec<BoundaryFace> {
    [[0,1,2],[0,1,3],[0,2,3],[1,2,3]].into_iter().map(|vertices|
        BoundaryFace { vertices, condition: Bc::Dirichlet([0.0;3]) }).collect()
}

#[test]
fn zero_mean_source_has_nonzero_certified_flux_defect_including_tensor_cross_terms() {
    let boundary = fixed_tet();
    // Exact correction b=-lambda_0*(x,y,z); no RT0 mean load remains.
    for (k, exact) in [(IDENTITY, 1.0/420.0),
        ([[2.,1.,0.],[1.,2.,0.],[0.,0.,4.]], 1.0/1008.0)] {
        let tensors = [k]; let source = [[-3.,1.,1.,1.]];
        let problem = AffineSourceTetProblem { vertices: &VERTICES, tets: &TETS,
            conductivity: &tensors, source: &source, boundary: &boundary };
        let result = affine_source_energy_bound(&problem, &[0.;4], FluxBudget::default(), || true).unwrap();
        contains(result.majorant_squared, exact);
        assert!(result.majorant_squared.hi-result.majorant_squared.lo < 1e-12);
        assert!(result.energy_error_upper > 0.03);
        let averaged = TensorTetProblem { vertices: &VERTICES, tets: &TETS,
            conductivity: &tensors, source: &[0.], boundary: &boundary };
        assert!(tensor_energy_bound(&averaged, &[0.;4], FluxBudget::default(), || true)
            .unwrap().energy_error_upper < 1e-12);
        for flux in result.outward_flux_integrals[0] { contains(flux, 0.0); }
    }
}

fn cube(n: usize, robin: bool) -> (Vec<[f64;3]>, Vec<[usize;4]>, Vec<BoundaryFace>) {
    let id = |i, j, k| (k*(n+1)+j)*(n+1)+i;
    let mut vertices = Vec::new();
    for k in 0..=n { for j in 0..=n { for i in 0..=n {
        vertices.push([i as f64/n as f64, j as f64/n as f64, k as f64/n as f64]);
    }}}
    let mut tets = Vec::new();
    for k in 0..n { for j in 0..n { for i in 0..n {
        for order in [[0,1,2],[0,2,1],[1,0,2],[1,2,0],[2,0,1],[2,1,0]] {
            let mut c = [i,j,k]; let mut tet = [id(i,j,k);4];
            for m in 0..3 { c[order[m]]+=1; tet[m+1]=id(c[0],c[1],c[2]); }
            tets.push(tet);
        }
    }}}
    let mut counts = BTreeMap::<[usize;3], usize>::new();
    for t in &tets { for opposite in 0..4 {
        let mut f = [0;3]; let mut j=0;
        for (i, v) in t.iter().enumerate() { if i!=opposite { f[j]=*v; j+=1; } }
        f.sort_unstable(); *counts.entry(f).or_default()+=1;
    }}
    let boundary = counts.into_iter().filter(|(_, count)| *count==1).map(|(face, _)| {
        let condition = if face.iter().all(|&v| vertices[v][0]==0.) { Bc::Dirichlet([2.;3]) }
            else if face.iter().all(|&v| vertices[v][0]==1.) {
                if robin { Bc::Robin { h: 2., reference: [1.;3] } } else { Bc::Dirichlet([2.;3]) }
            } else { Bc::Neumann(0.) };
        BoundaryFace { vertices: face, condition }
    }).collect();
    (vertices,tets,boundary)
}

#[test]
fn manufactured_cubic_temperature_mean_and_true_energy_error_are_enclosed() {
    // T=2+x-x^3, f=6*x, with either fixed T=2 at x=1 or h=2,Tref=1.
    // Mean is exactly 9/4. The source is genuinely nonconstant in every tet.
    for robin in [false,true] { for n in [1,2] {
        let (vertices,tets,boundary) = cube(n,robin);
        let source: Vec<_> = tets.iter().map(|t| t.map(|i| 6.*vertices[i][0])).collect();
        let conductivity = vec![IDENTITY;tets.len()];
        let problem = AffineSourceTetProblem { vertices:&vertices, tets:&tets,
            source:&source, conductivity:&conductivity, boundary:&boundary };
        let u: Vec<_> = vertices.iter().map(|p| 2.+p[0]-p[0]*p[0]*p[0]).collect();
        // A conforming, deliberately inexact unit-source dual is sufficient.
        let z: Vec<_> = vertices.iter().map(|p| p[0]*(1.-p[0])/2.).collect();
        let result = affine_source_mean_bound(&problem,&u,&z,FluxBudget::default(),||true).unwrap();
        contains(result.enclosure,2.25);
        let h=1./n as f64;
        let energy_squared: f64 = (0..n).map(|i| {
            let m=(i as f64+0.5)*h;
            // integral (T' - interpolation slope)^2 over this slab.
            3.*m*m*h.powi(3)+h.powi(5)/20.
        }).sum();
        assert!(result.integral.primal.energy_error_upper >= energy_squared.sqrt());
        assert!(result.enclosure.hi-result.enclosure.lo < 8.0);
        for (flux,t) in result.integral.primal.outward_flux_integrals.iter().zip(&tets) {
            let sum=flux.iter().copied().fold(Iv::zero(),Iv::add);
            let exact_source_mean: f64=t.iter().map(|&v| 6.*vertices[v][0]/4.).sum();
            contains(sum,exact_source_mean*h.powi(3)/6.);
        }
    }}
}

#[test]
fn constant_sources_reproduce_original_energy_and_mean_bits() {
    let (vertices,tets,boundary)=cube(2,true);
    let k=vec![IDENTITY;tets.len()]; let constant=vec![2.;tets.len()]; let affine=vec![[2.;4];tets.len()];
    let legacy=TensorTetProblem { vertices:&vertices,tets:&tets,conductivity:&k,source:&constant,boundary:&boundary };
    let lifted=AffineSourceTetProblem { vertices:&vertices,tets:&tets,conductivity:&k,source:&affine,boundary:&boundary };
    let u=vec![2.;vertices.len()];let z=vec![0.;vertices.len()];
    let a=tensor_mean_bound(&legacy,&u,&z,FluxBudget::default(),||true).unwrap();
    let b=affine_source_mean_bound(&lifted,&u,&z,FluxBudget::default(),||true).unwrap();
    assert_eq!(a.enclosure,b.enclosure);assert_eq!(a.candidate_mean,b.candidate_mean);
    assert_eq!(a.integral.residual_correction,b.integral.residual_correction);
    assert_eq!(a.integral.primal.majorant_squared,b.integral.primal.majorant_squared);
    assert_eq!(a.integral.primal.outward_flux_integrals,b.integral.primal.outward_flux_integrals);
}

#[test]
fn source_order_and_actual_goal_load_are_preserved() {
    let boundary:Vec<_>=fixed_tet().into_iter().map(|f|BoundaryFace { condition:Bc::Robin {h:1.,reference:[0.;3]}, ..f }).collect();
    let k=[IDENTITY];let source=[[0.,4.,0.,0.]];
    let p=AffineSourceTetProblem {vertices:&VERTICES,tets:&TETS,conductivity:&k,source:&source,boundary:&boundary};
    let goal=affine_source_goal_bound(&p,&[0.;4],&[0.,1.,0.,0.],&[1.],FluxBudget::default(),||true).unwrap();
    contains(goal.residual_correction,1./15.);
    assert!(goal.residual_correction.lo>1./24.);
    let permuted=[[2,0,3,1]];let f=[[0.,0.,0.,4.]];
    let other=AffineSourceTetProblem {tets:&permuted,source:&f,..p};
    let replay=affine_source_goal_bound(&other,&[0.;4],&[0.,1.,0.,0.],&[1.],FluxBudget::default(),||true).unwrap();
    contains(replay.residual_correction,1./15.);
}

#[test]
fn budgets_bad_data_and_late_cancellation_never_return_a_partial_bound() {
    let boundary=fixed_tet();let k=[IDENTITY];let source=[[-3.,1.,1.,1.]];
    let p=AffineSourceTetProblem {vertices:&VERTICES,tets:&TETS,conductivity:&k,source:&source,boundary:&boundary};
    let budget=FluxBudget::default();
    assert!(matches!(affine_source_energy_bound(&p,&[0.;4],FluxBudget {max_cells:0,..budget},||true),Err(TetError::Budget)));
    let mut polls=0;let before=affine_source_energy_bound(&p,&[0.;4],budget,||{polls+=1;true}).unwrap();
    for stop in [1,polls/2,polls] {
        let mut count=0;
        assert!(matches!(affine_source_energy_bound(&p,&[0.;4],budget,||{count+=1;count!=stop}),Err(TetError::Cancelled)));
    }
    let bad=AffineSourceTetProblem {source:&[],..p};
    assert!(affine_source_energy_bound(&bad,&[0.;4],budget,||true).is_err());
    let bad_source=[[0.,f64::NAN,0.,0.]];
    assert!(affine_source_energy_bound(&AffineSourceTetProblem {source:&bad_source,..p},&[0.;4],budget,||true).is_err());
    assert!(affine_source_energy_bound(&p,&[1.;4],budget,||true).is_err());
    assert_eq!(before.majorant_squared,affine_source_energy_bound(&p,&[0.;4],budget,||true).unwrap().majorant_squared);
}
