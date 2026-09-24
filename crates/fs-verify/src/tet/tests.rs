use super::*;

struct Fixture {
    vertices: Vec<[f64; 3]>,
    tets: Vec<[usize; 4]>,
    conductivity: Vec<f64>,
    source: Vec<f64>,
    boundary: Vec<BoundaryFace>,
}
impl Fixture {
    fn problem(&self) -> TetProblem<'_> {
        TetProblem { vertices: &self.vertices, tets: &self.tets,
            conductivity: &self.conductivity, source: &self.source, boundary: &self.boundary }
    }
}
fn cube(n: usize, source: f64, value: impl Fn(f64) -> f64) -> (Fixture, Vec<f64>) {
    let mut vertices = Vec::new();
    for z in 0..=n { for y in 0..=n { for x in 0..=n {
        vertices.push([x as f64/n as f64, y as f64/n as f64, z as f64/n as f64]);
    }}}
    let index = |p: [usize; 3]| p[0]+(n+1)*(p[1]+(n+1)*p[2]);
    let mut tets = Vec::new();
    for z in 0..n { for y in 0..n { for x in 0..n {
        for permutation in [[0,1,2], [0,2,1], [1,0,2], [1,2,0], [2,0,1], [2,1,0]] {
            let mut p = [x,y,z];
            let mut tet = [index(p); 4];
            for (i, d) in permutation.into_iter().enumerate() { p[d] += 1; tet[i+1] = index(p); }
            tets.push(tet);
        }
    }}}
    let mut incidence = BTreeMap::<[usize; 3], usize>::new();
    for tet in &tets {
        for opposite in 0..4 {
            let mut face = [0; 3];
            let mut j = 0;
            for (i, &v) in tet.iter().enumerate() {
                if i != opposite { face[j] = v; j += 1; }
            }
            face.sort_unstable();
            *incidence.entry(face).or_default() += 1;
        }
    }
    let u: Vec<_> = vertices.iter().map(|p| value(p[0])).collect();
    let boundary = incidence.into_iter().filter(|(_, count)| *count == 1)
        .map(|(face, _)| {
            let prescribed = face.iter().all(|&i| vertices[i][0] == 0.0)
                || face.iter().all(|&i| vertices[i][0] == 1.0);
            BoundaryFace { vertices: face, condition: if prescribed {
                BoundaryCondition::Dirichlet(face.map(|i| u[i]))
            } else { BoundaryCondition::Neumann(0.0) } }
        }).collect();
    let count = tets.len();
    (Fixture { vertices, tets, conductivity: vec![1.0; count], source: vec![source; count], boundary }, u)
}

#[test]
fn affine_patch_has_only_outward_rounding_width() {
    let (fixture, u) = cube(2, 0.0, |x| 300.0+2.0*x);
    let bound = energy_bound(&fixture.problem(), &u, FluxBudget::default(), || true).unwrap();
    assert!(bound.energy_error_upper < 1e-7, "{bound:?}");
    assert!(bound.energy_error_upper >= 0.0);
}

#[test]
fn quadratic_slab_bounds_true_energy_error_and_shrinks_under_refinement() {
    let mut previous = f64::INFINITY;
    for n in [1, 2, 4] {
        let (fixture, u) = cube(n, 2.0, |x| x*(1.0-x));
        let bound = energy_bound(&fixture.problem(), &u, FluxBudget::default(), || true).unwrap();
        // Exact continuum solution x(1-x); exact interpolation error h/sqrt(3).
        let exact = 1.0/(3.0_f64.sqrt()*n as f64);
        assert!(bound.energy_error_upper >= exact, "underbound n={n}: {bound:?}");
        assert!(bound.energy_error_upper < 1.6*exact, "fixture effectivity n={n}: {bound:?}");
        assert!(bound.energy_error_upper < 0.6*previous);
        previous = bound.energy_error_upper;
        for (e, flux) in bound.outward_flux_integrals.iter().enumerate() {
            let sum = flux.iter().copied().fold(Iv::zero(), Iv::add);
            // Each dyadic cube has six equal-volume tets. This check is a
            // diagnostic; the proof uses the correlated forest expressions.
            let exact_load = 2.0/(6.0*(n*n*n) as f64);
            assert!(sum.lo <= exact_load && exact_load <= sum.hi, "cell {e}");
        }
    }
}

#[test]
fn robin_trace_defect_cannot_be_dropped() {
    let (mut fixture, u) = cube(2, 0.0, |_| 0.0);
    for face in &mut fixture.boundary {
        face.condition = BoundaryCondition::Robin { h: 2.0, reference: [1.0; 3] };
    }
    // Exact solution is 1. Candidate is 0. All error is the Robin trace:
    // ||e||_a^2 = 2 * surface_area(unit cube) = 12.
    let bound = energy_bound(&fixture.problem(), &u, FluxBudget::default(), || true).unwrap();
    assert!(bound.energy_error_upper >= 12.0_f64.sqrt());
    assert!(bound.energy_error_upper < 12.0_f64.sqrt()*(1.0+1e-10));
}

#[test]
fn proposal_exhaustion_does_not_turn_a_residual_into_equilibrium() {
    let (fixture, u) = cube(4, 2.0, |x| x*(1.0-x));
    let budget = FluxBudget { max_iterations: 0, ..FluxBudget::default() };
    let bound = energy_bound(&fixture.problem(), &u, budget, || true).unwrap();
    assert_eq!(bound.proposal_iterations, 0);
    assert!(bound.energy_error_upper >= 1.0/(4.0*3.0_f64.sqrt()));
    assert!(bound.energy_error_upper < 0.5);
}

#[test]
fn outward_neumann_sign_is_physical_not_an_absolute_value() {
    let (mut fixture, u) = cube(2, 0.0, |x| x);
    for face in &mut fixture.boundary {
        if face.vertices.iter().all(|&i| fixture.vertices[i][0] == 1.0) {
            face.condition = BoundaryCondition::Neumann(-1.0);
        }
    }
    let exact = energy_bound(&fixture.problem(), &u, FluxBudget::default(), || true).unwrap();
    assert!(exact.energy_error_upper < 1e-7);
    for face in &mut fixture.boundary {
        if matches!(face.condition, BoundaryCondition::Neumann(q) if q == -1.0) {
            face.condition = BoundaryCondition::Neumann(1.0);
        }
    }
    // Opposite load has solution -x; the retained +x candidate is wrong by 2 in energy norm.
    let wrong = energy_bound(&fixture.problem(), &u, FluxBudget::default(), || true).unwrap();
    assert!(wrong.energy_error_upper >= 2.0);
}

#[test]
fn invalid_trace_unanchored_domain_and_resource_stops_return_no_bound() {
    let (mut fixture, mut u) = cube(1, 0.0, |x| x);
    u[0] = 7.0;
    assert!(matches!(energy_bound(&fixture.problem(), &u, FluxBudget::default(), || true), Err(TetError::Invalid(_))));
    u[0] = 0.0;
    assert_eq!(energy_bound(&fixture.problem(), &u, FluxBudget::default(), || false).unwrap_err(), TetError::Cancelled);
    assert_eq!(energy_bound(&fixture.problem(), &u, FluxBudget { max_cells: 1, max_iterations: 0 }, || true).unwrap_err(), TetError::Budget);
    for face in &mut fixture.boundary { face.condition = BoundaryCondition::Neumann(0.0); }
    assert!(matches!(energy_bound(&fixture.problem(), &u, FluxBudget::default(), || true), Err(TetError::Unsupported(_))));
    fixture.conductivity[0] = f64::NAN;
    assert!(matches!(energy_bound(&fixture.problem(), &u, FluxBudget::default(), || true), Err(TetError::Invalid(_))));
}

#[test]
fn duplicate_cells_and_missing_boundary_are_not_silently_admitted() {
    let (mut fixture, u) = cube(1, 0.0, |x| x);
    fixture.tets[1] = fixture.tets[0];
    assert!(energy_bound(&fixture.problem(), &u, FluxBudget::default(), || true).is_err());
    let (mut fixture, u) = cube(1, 0.0, |x| x);
    fixture.boundary.pop();
    assert!(energy_bound(&fixture.problem(), &u, FluxBudget::default(), || true).is_err());
}
