use super::*;

struct Slab {
    vertices: Vec<[f64; 3]>,
    tets: Vec<[usize; 4]>,
    conductivity: Vec<ConductivityTensor>,
    source: Vec<[f64; 4]>,
    boundary: Vec<BoundaryFace>,
}
impl Slab {
    fn new(source: impl Fn(f64) -> f64) -> Self {
        let vertices: Vec<_> = (0..8).map(|i|
            [(i & 1) as f64, ((i >> 1) & 1) as f64, ((i >> 2) & 1) as f64]).collect();
        let tets: Vec<_> = [[0,1,2], [0,2,1], [1,0,2], [1,2,0], [2,0,1], [2,1,0]]
            .into_iter().map(|axes| {
                let mut tet = [0; 4];
                for (i, axis) in axes.into_iter().enumerate() { tet[i+1] = tet[i] | (1 << axis); }
                tet
            }).collect();
        let mut incidence = BTreeMap::<[usize; 3], usize>::new();
        for tet in &tets {
            for opposite in 0..4 {
                let mut face = [0; 3]; let mut slot = 0;
                for (i, &v) in tet.iter().enumerate() {
                    if i != opposite { face[slot] = v; slot += 1; }
                }
                face.sort_unstable(); *incidence.entry(face).or_default() += 1;
            }
        }
        let boundary = incidence.into_iter().filter(|(_, n)| *n == 1).map(|(face, _)| BoundaryFace {
            vertices: face,
            condition: if face.iter().all(|&v| vertices[v][0] == 1.0) {
                BoundaryCondition::Dirichlet([300.0; 3])
            } else { BoundaryCondition::Neumann(0.0) },
        }).collect();
        let loads = tets.iter().map(|tet| tet.map(|v| source(vertices[v][0]))).collect();
        Self { vertices, conductivity: vec![[[2.0,0.0,0.0],[0.0,3.0,0.0],[0.0,0.0,4.0]]; tets.len()],
            tets, source: loads, boundary }
    }
    fn problem(&self) -> AffineSourceTetProblem<'_> {
        AffineSourceTetProblem { vertices: &self.vertices, tets: &self.tets,
            conductivity: &self.conductivity, source: &self.source, boundary: &self.boundary }
    }
    fn range(&self) -> Result<Iv, TetError> {
        self.problem().slab_temperature_range(&vec![300.0; self.vertices.len()], 0, false,
            FluxBudget::default(), || true)
    }
}

#[test]
fn signed_heating_and_cooling_enclose_exact_continuum_extrema() {
    for (source, expected) in [(4.0, [300.0,301.0]), (-4.0, [299.0,300.0]), (0.0,[300.0,300.0])] {
        let slab = Slab::new(|_| source);
        let range = slab.range().unwrap();
        assert!(range.lo <= expected[0] && range.hi >= expected[1]);
        assert!((range.lo-expected[0]).abs() < 1e-10 && (range.hi-expected[1]).abs() < 1e-10);
        // This unconverged field is not the analytic solution for nonzero f.
        assert_eq!(range, slab.problem().slab_temperature_range(&[300.0;8],0,false,
            FluxBudget { max_cells: 6, max_iterations: 0 },||true).unwrap());
    }
}

#[test]
fn robin_end_and_affine_source_keep_their_boundary_and_source_allowances() {
    let mut slab = Slab::new(|_| 4.0);
    for face in &mut slab.boundary {
        if matches!(face.condition, BoundaryCondition::Dirichlet(_)) {
            face.condition = BoundaryCondition::Robin { h:4.0, reference:[300.0;3] };
        }
    }
    // T=302-x^2: q(1)=4=4*(T(1)-300). Omitting L/h misses the maximum.
    let range = slab.range().unwrap();
    assert!(range.lo <= 300.0 && range.hi >= 302.0 && range.hi < 302.0+1e-10);
    let affine = Slab::new(|x| 12.0*x);
    // T=301-x^3; the actual P1-source supremum is 12, not its mean 6.
    let range = affine.range().unwrap();
    assert!(range.hi >= 303.0 && range.hi < 303.0+1e-10);
    let mut sign_changing = Slab::new(|x| 8.0*x-4.0);
    let first = sign_changing.range().unwrap();
    assert!(first.lo <= 299.0 && first.hi >= 301.0);
    sign_changing.source.iter_mut().for_each(|row| row.fill(0.0));
    assert_eq!(sign_changing.range().unwrap(), Iv::point(300.0));
}

#[test]
fn rotated_reversed_and_tangentially_heterogeneous_slabs_preserve_the_bound() {
    let mut slab = Slab::new(|_| 4.0);
    for (i, k) in slab.conductivity.iter_mut().enumerate() {
        k[1][1] += i as f64; k[1][2] = 0.5; k[2][1] = 0.5;
    }
    let expected = slab.range().unwrap();
    for p in &mut slab.vertices { p.swap(0,2); p[2] = 1.0-p[2]; }
    for k in &mut slab.conductivity {
        *k = [2,1,0].map(|i| [2,1,0].map(|j| k[i][j]));
    }
    let actual = slab.problem().slab_temperature_range(&[300.0;8],2,true,
        FluxBudget::default(),||true).unwrap();
    assert_eq!(actual, expected);
    assert!(slab.problem().slab_temperature_range(&[300.0;8],2,false,
        FluxBudget::default(),||true).is_err());
}

#[test]
fn unsupported_physics_and_geometry_never_become_maximum_certificates() {
    let mut slab = Slab::new(|_|4.0);
    slab.conductivity[0][0][0] = 2.5;
    assert!(matches!(slab.range(),Err(TetError::Unsupported(_))));
    let mut slab = Slab::new(|_|4.0);
    slab.conductivity.iter_mut().for_each(|k| { k[0][1]=0.25; k[1][0]=0.25; });
    assert!(matches!(slab.range(),Err(TetError::Unsupported(_))));
    let mut slab = Slab::new(|_|4.0);
    slab.boundary.iter_mut().for_each(|f| f.condition = BoundaryCondition::Neumann(0.0));
    assert!(slab.range().is_err());
    let mut slab = Slab::new(|_|4.0);
    let i = slab.boundary.iter().position(|f| matches!(f.condition,BoundaryCondition::Neumann(_))).unwrap();
    slab.boundary[i].condition = BoundaryCondition::Neumann(1.0);
    assert!(matches!(slab.range(),Err(TetError::Unsupported(_))));
    let mut slab = Slab::new(|_|4.0);
    slab.vertices[7][0] = 0.875;
    assert!(slab.range().is_err());
}

#[test]
fn invalid_fields_work_limits_and_cancellation_return_no_range() {
    let slab = Slab::new(|_|4.0); let problem = slab.problem();
    let b = FluxBudget::default();
    assert!(matches!(problem.slab_temperature_range(&[300.0;8],3,false,b,||true),Err(TetError::Invalid(_))));
    assert!(problem.slab_temperature_range(&[300.0;7],0,false,b,||true).is_err());
    assert!(problem.slab_temperature_range(&[299.0;8],0,false,b,||true).is_err());
    assert!(matches!(problem.slab_temperature_range(&[300.0;8],0,false,
        FluxBudget { max_cells: 5, max_iterations: 0 },||true),Err(TetError::Budget)));
    assert!(matches!(problem.slab_temperature_range(&[300.0;8],0,false,b,||false),Err(TetError::Cancelled)));
    let mut polls = 0;
    assert!(matches!(problem.slab_temperature_range(&[300.0;8],0,false,b,|| {
        polls += 1; polls < 10
    }),Err(TetError::Cancelled)));
}
