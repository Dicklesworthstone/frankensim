use super::*;

const COUPLED: ConductivityTensor = [[2.0, -1.0, -1.0], [-1.0, 2.0, 0.0], [-1.0, 0.0, 2.0]];

struct Fixture {
    vertices: Vec<[f64; 3]>,
    tets: Vec<[usize; 4]>,
    conductivity: Vec<ConductivityTensor>,
    source: Vec<f64>,
    boundary: Vec<BoundaryFace>,
}
impl Fixture {
    fn problem(&self) -> TensorTetProblem<'_> {
        TensorTetProblem { vertices: &self.vertices, tets: &self.tets, conductivity: &self.conductivity,
            source: &self.source, boundary: &self.boundary }
    }
}
fn mixed_tet() -> Fixture {
    // Exact solution u=x; q=-K grad u=(-2,1,1). The sloping face has
    // ZERO flux, so all declared normal fluxes are exact rational values.
    Fixture {
        vertices: vec![[0.0; 3], [1.0,0.0,0.0], [0.0,1.0,0.0], [0.0,0.0,1.0]],
        tets: vec![[0,1,2,3]], conductivity: vec![COUPLED], source: vec![0.0],
        boundary: vec![
            BoundaryFace { vertices: [0,2,3], condition: BoundaryCondition::Dirichlet([0.0;3]) },
            BoundaryFace { vertices: [0,1,3], condition: BoundaryCondition::Neumann(-1.0) },
            BoundaryFace { vertices: [0,1,2], condition: BoundaryCondition::Neumann(-1.0) },
            BoundaryFace { vertices: [1,2,3], condition: BoundaryCondition::Neumann(0.0) },
        ],
    }
}
fn slab(n: usize) -> (Fixture, Vec<f64>) {
    let mut vertices = Vec::new();
    for z in 0..=n { for y in 0..=n { for x in 0..=n {
        vertices.push([x as f64/n as f64, y as f64/n as f64, z as f64/n as f64]);
    }}}
    let index = |p: [usize; 3]| p[0]+(n+1)*(p[1]+(n+1)*p[2]);
    let mut tets = Vec::new();
    for z in 0..n { for y in 0..n { for x in 0..n {
        for axes in [[0,1,2], [0,2,1], [1,0,2], [1,2,0], [2,0,1], [2,1,0]] {
            let mut p = [x,y,z];
            let mut tet = [index(p); 4];
            for (i,d) in axes.into_iter().enumerate() { p[d] += 1; tet[i+1] = index(p); }
            tets.push(tet);
        }
    }}}
    let mut counts = BTreeMap::<[usize;3], usize>::new();
    for tet in &tets { for skip in 0..4 {
        let mut face = [0;3];
        let mut j = 0;
        for (i,&v) in tet.iter().enumerate() { if i != skip { face[j] = v; j += 1; } }
        face.sort_unstable();
        *counts.entry(face).or_default() += 1;
    }}
    let u: Vec<_> = vertices.iter().map(|p| p[0]*(1.0-p[0])).collect();
    let boundary = counts.into_iter().filter(|(_,c)| *c==1).map(|(face,_)| {
        let end = face.iter().all(|&i| vertices[i][0]==0.0) || face.iter().all(|&i| vertices[i][0]==1.0);
        BoundaryFace { vertices: face, condition: if end {
            BoundaryCondition::Dirichlet([0.0;3])
        } else { BoundaryCondition::Neumann(0.0) } }
    }).collect();
    let ne = tets.len();
    (Fixture { vertices, tets, conductivity: vec![[[4.0,0.0,0.0],[0.0,2.0,0.0],[0.0,0.0,1.0]];ne],
        source: vec![8.0;ne], boundary }, u)
}

#[test]
fn full_inverse_metric_contains_the_exact_cross_terms_at_extreme_scales() {
    for s in [1e-150, 1.0, 1e150] {
        let k = [[2.0,-1.0,0.0],[-1.0,2.0,0.0],[0.0,0.0,1.0]].map(|r| r.map(|v| s*v));
        let t = Tensor::new(k).unwrap();
        let q = t.inverse_quadratic([Iv::point(1.0), Iv::point(1.0), Iv::zero()]);
        let exact = 2.0/s;
        assert!(q.lo <= exact && exact <= q.hi, "{q:?} versus {exact}");
        assert!(q.hi < exact*(1.0+1e-12));
    }
}

#[test]
fn off_diagonal_fluxes_recover_an_affine_mixed_boundary_patch() {
    let fixture = mixed_tet();
    let u = [0.0,1.0,0.0,0.0];
    let b = tensor_energy_bound(&fixture.problem(), &u, FluxBudget::default(), || true).unwrap();
    assert!(b.energy_error_upper < 1e-7, "{b:?}");
    // Replacing K by its diagonal changes the PDE, not just the evidence label.
    let mut wrong = mixed_tet();
    wrong.conductivity[0] = [[2.0,0.0,0.0],[0.0,2.0,0.0],[0.0,0.0,2.0]];
    let b = tensor_energy_bound(&wrong.problem(), &u, FluxBudget::default(), || true).unwrap();
    assert!(b.energy_error_upper > 0.1);
}

#[test]
fn anisotropic_energy_covers_an_exact_algebraic_error_with_no_proposal_iterations() {
    let fixture = mixed_tet();
    let b = tensor_energy_bound(&fixture.problem(), &[0.0;4],
        FluxBudget { max_iterations: 0, ..FluxBudget::default() }, || true).unwrap();
    // u=x on a volume-1/6 tetrahedron: ||u-0||_K^2=Kxx/6=1/3.
    let exact = (1.0_f64/3.0).sqrt();
    assert!(b.energy_error_upper >= exact, "{b:?}");
    assert!(b.energy_error_upper < exact*(1.0+1e-10));
    assert_eq!(b.proposal_iterations,0);
}

#[test]
fn anisotropic_energy_and_mean_bounds_refine_on_the_manufactured_slab() {
    let mut energy = f64::INFINITY;
    let mut width = f64::INFINITY;
    for n in [1,2,4] {
        let (fixture,u) = slab(n);
        let z: Vec<_> = u.iter().map(|v| v/8.0).collect();
        let b = tensor_mean_bound(&fixture.problem(), &u, &z, FluxBudget::default(), || true).unwrap();
        let exact_energy = 2.0/(3.0_f64.sqrt()*n as f64);
        assert!(b.integral.primal.energy_error_upper >= exact_energy, "{b:?}");
        assert!(b.integral.primal.energy_error_upper < 0.8*energy);
        assert!(b.enclosure.lo <= 1.0/6.0 && b.enclosure.hi >= 1.0/6.0, "{b:?}");
        let next_width = b.enclosure.hi-b.enclosure.lo;
        assert!(next_width < 0.6*width, "{next_width} vs {width}");
        energy = b.integral.primal.energy_error_upper;
        width = next_width;
    }
}

#[test]
fn tensor_goal_retains_the_residual_and_signed_functional() {
    let fixture = mixed_tet();
    let u = [0.0,1.0,0.0,0.0];
    let z = [0.0,0.25,0.0,0.0];
    let b = tensor_goal_bound(&fixture.problem(), &u, &z, &[-2.0], FluxBudget::default(), || true).unwrap();
    // Integral -2*x on the reference tet equals -1/12.
    assert!(b.enclosure.lo <= -1.0/12.0 && b.enclosure.hi >= -1.0/12.0, "{b:?}");
    assert!(b.enclosure.hi-b.enclosure.lo < 1e-7);
    assert!(b.residual_correction.lo <= 0.0 && b.residual_correction.hi >= 0.0);
    let (fixture, u) = slab(4);
    let zero = vec![0.0;u.len()];
    let z: Vec<_> = u.iter().map(|v| v/8.0).collect();
    let b = tensor_mean_bound(&fixture.problem(), &zero, &z, FluxBudget::default(), || true).unwrap();
    assert!(b.integral.residual_correction.lo > 0.14);
    assert!(b.enclosure.lo <= 1.0/6.0 && b.enclosure.hi >= 1.0/6.0);
}

#[test]
fn rotating_geometry_and_conductivity_preserves_the_physical_bound() {
    let mut fixture = mixed_tet();
    let u = [0.0;4];
    let original = tensor_energy_bound(&fixture.problem(), &u, FluxBudget::default(), || true).unwrap();
    let axes = [2,0,1]; // proper, exact rotation; no irrational coordinates
    for p in &mut fixture.vertices { *p = axes.map(|i| p[i]); }
    for k in &mut fixture.conductivity { *k = axes.map(|i| axes.map(|j| k[i][j])); }
    let rotated = tensor_energy_bound(&fixture.problem(), &u, FluxBudget::default(), || true).unwrap();
    assert!((rotated.energy_error_upper-original.energy_error_upper).abs() < 1e-11);
}

#[test]
fn anisotropic_robin_trace_energy_is_not_lost() {
    let mut fixture = mixed_tet();
    for face in &mut fixture.boundary {
        face.condition = BoundaryCondition::Robin { h: 2.0, reference: [1.0;3] };
    }
    let b = tensor_energy_bound(&fixture.problem(), &[0.0;4], FluxBudget::default(), || true).unwrap();
    // Exact solution is one: only Robin trace energy remains.
    let exact = (3.0+3.0_f64.sqrt()).sqrt();
    assert!(b.energy_error_upper >= exact);
    assert!(b.energy_error_upper < exact*(1.0+1e-10));
}

#[test]
fn invalid_or_unprovable_tensors_and_resource_interruptions_return_no_bound() {
    let mut nonsymmetric = COUPLED;
    nonsymmetric[0][1] = 0.0;
    let mut nonfinite = COUPLED;
    nonfinite[2][2] = f64::NAN;
    for k in [nonsymmetric, nonfinite,
        [[1.0,2.0,0.0],[2.0,1.0,0.0],[0.0,0.0,1.0]],
        [[1.0,1.0,0.0],[1.0,1.0+f64::EPSILON,0.0],[0.0,0.0,1.0]],
        [[0.0;3];3],
    ] { assert!(Tensor::new(k).is_err()); }
    let fixture = mixed_tet();
    assert_eq!(tensor_energy_bound(&fixture.problem(), &[0.0;4], FluxBudget::default(), || false).unwrap_err(), TetError::Cancelled);
    assert_eq!(tensor_energy_bound(&fixture.problem(), &[0.0;4], FluxBudget { max_cells: 0, max_iterations: 0 }, || true).unwrap_err(), TetError::Budget);
    let mut calls = 0;
    assert_eq!(tensor_energy_bound(&fixture.problem(), &[0.0;4], FluxBudget::default(), || {
        calls += 1; calls < 5
    }).unwrap_err(), TetError::Cancelled);
    let wrong = TensorTetProblem { conductivity: &[], ..fixture.problem() };
    assert!(matches!(tensor_energy_bound(&wrong, &[0.0;4], FluxBudget::default(), || true), Err(TetError::Invalid(_))));
}
