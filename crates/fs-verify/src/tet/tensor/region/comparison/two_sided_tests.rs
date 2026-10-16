use super::*;

struct Slab {
    vertices: Vec<[f64; 3]>,
    tets: Vec<[usize; 4]>,
    conductivity: Vec<ConductivityTensor>,
    source: Vec<[f64; 4]>,
    boundary: Vec<BoundaryFace>,
}
impl Slab {
    fn new(source: impl Fn(f64) -> f64, h: [f64; 2]) -> Self {
        let vertices: Vec<_> = (0..8).map(|i|
            [(i & 1) as f64, ((i >> 1) & 1) as f64, ((i >> 2) & 1) as f64]).collect();
        let tets: Vec<_> = [[0,1,2],[0,2,1],[1,0,2],[1,2,0],[2,0,1],[2,1,0]]
            .into_iter().map(|axes| {
                let mut tet = [0;4];
                for (i, axis) in axes.into_iter().enumerate() { tet[i+1] = tet[i] | (1 << axis); }
                tet
            }).collect();
        let mut incidence = BTreeMap::<[usize;3], usize>::new();
        for tet in &tets { for opposite in 0..4 {
            let mut face = [0;3]; let mut slot = 0;
            for (i, &v) in tet.iter().enumerate() {
                if i != opposite { face[slot] = v; slot += 1; }
            }
            face.sort_unstable(); *incidence.entry(face).or_default() += 1;
        }}
        let boundary = incidence.into_iter().filter(|(_, n)| *n == 1).map(|(face, _)| {
            let end = [0,1].into_iter().find(|&e| face.iter().all(|&v| vertices[v][0] == e as f64));
            BoundaryFace { vertices: face, condition: match end {
                Some(i) if h[i].is_infinite() => BoundaryCondition::Dirichlet([300.0;3]),
                Some(i) => BoundaryCondition::Robin { h: h[i], reference: [300.0;3] },
                None => BoundaryCondition::Neumann(0.0),
            }}
        }).collect();
        let source = tets.iter().map(|tet| tet.map(|v| source(vertices[v][0]))).collect();
        let conductivity = vec![[[2.0,0.0,0.0],[0.0,3.0,0.0],[0.0,0.0,4.0]];tets.len()];
        Self { vertices, tets, conductivity, source, boundary }
    }
    fn problem(&self) -> AffineSourceTetProblem<'_> {
        AffineSourceTetProblem { vertices: &self.vertices, tets: &self.tets,
            conductivity: &self.conductivity, source: &self.source, boundary: &self.boundary }
    }
    fn range(&self) -> Result<Iv, TetError> {
        self.problem().slab_temperature_range_two_sided(&[300.0;8],0,FluxBudget::default(),||true)
    }
}

#[test]
fn both_ends_enclose_exact_symmetric_asymmetric_and_mixed_boundary_extrema() {
    // Exact stationary points of -2*T''=4 with the indicated end laws:
    // RR(4,4): x=1/2; RR(2,8): x=1/3; RD(4,D): x=1/3; DD: x=1/2.
    for (h, rise) in [([4.0,4.0],0.75),([2.0,8.0],7.0/9.0),
        ([4.0,f64::INFINITY],4.0/9.0),([f64::INFINITY,f64::INFINITY],0.25)]
    {
        for sign in [-1.0,1.0] {
            let slab = Slab::new(|_| sign*4.0,h);
            let range = slab.range().unwrap();
            let exact = 300.0+sign*rise;
            assert!(range.lo <= exact && exact <= range.hi, "{h:?} {sign}: {range:?} vs {exact}");
            let endpoint = if sign > 0.0 { range.hi } else { range.lo };
            assert!((endpoint-exact).abs() < 2e-11);
            let old = slab.problem().slab_temperature_range(&[300.0;8],0,false,
                FluxBudget::default(),||true).unwrap();
            assert!(range.hi-range.lo < old.hi-old.lo, "both ends must actually improve the bound");
        }
    }
}

#[test]
fn weakest_patch_and_all_reference_temperatures_are_retained() {
    let mut slab = Slab::new(|_|4.0,[2.0,8.0]);
    let expected = slab.range().unwrap();
    let left = slab.boundary.iter().position(|f| f.vertices.iter().all(|&v| slab.vertices[v][0]==0.0)).unwrap();
    // One stronger patch cannot replace the other patch's weaker coefficient.
    slab.boundary[left].condition = BoundaryCondition::Robin { h:100.0,reference:[300.0;3] };
    assert_eq!(slab.range().unwrap(),expected);
    slab.boundary[left].condition = BoundaryCondition::Dirichlet([300.0;3]);
    assert_eq!(slab.range().unwrap(),expected,"one Dirichlet patch is not an all-Dirichlet end");
    let mut slab = Slab::new(|_|0.0,[2.0,8.0]);
    let side = slab.boundary.iter().position(|f| matches!(f.condition,BoundaryCondition::Neumann(_))).unwrap();
    slab.boundary[side].condition = BoundaryCondition::Robin { h:3.0,reference:[280.0,320.0,340.0] };
    assert_eq!(slab.range().unwrap(),Iv { lo:280.0,hi:340.0 });
    let slab = Slab::new(|x|8.0*x-4.0,[4.0,4.0]);
    let range = slab.range().unwrap();
    assert!(range.lo <= 299.25 && range.hi >= 300.75,
        "a zero-mean signed source cannot be replaced by zero heating");
}

#[test]
fn two_end_bound_is_independent_of_candidate_axis_and_end_orientation() {
    let mut slab = Slab::new(|_|4.0,[2.0,8.0]);
    for (i,k) in slab.conductivity.iter_mut().enumerate() {
        k[1][1] += i as f64; k[1][2]=0.5; k[2][1]=0.5;
    }
    let expected = slab.range().unwrap();
    assert_eq!(expected,slab.problem().slab_temperature_range_two_sided(&[123.0;8],0,
        FluxBudget::default(),||true).unwrap(),"nodal extrema do not enter the comparison");
    for p in &mut slab.vertices { p.swap(0,2); p[2]=1.0-p[2]; }
    for k in &mut slab.conductivity { *k = [2,1,0].map(|i| [2,1,0].map(|j| k[i][j])); }
    let rotated = slab.problem().slab_temperature_range_two_sided(&[300.0;8],2,
        FluxBudget::default(),||true).unwrap();
    assert_eq!(expected,rotated);
}

#[test]
fn partial_cooling_and_numerical_or_work_failures_never_produce_a_range() {
    let mut slab = Slab::new(|_|4.0,[4.0,4.0]);
    let left = slab.boundary.iter().position(|f| f.vertices.iter().all(|&v| slab.vertices[v][0]==0.0)).unwrap();
    slab.boundary[left].condition=BoundaryCondition::Neumann(0.0);
    assert!(matches!(slab.range(),Err(TetError::Unsupported(_))));
    assert!(slab.problem().slab_temperature_range(&[300.0;8],0,false,
        FluxBudget::default(),||true).is_ok(),"the original one-ended fallback remains valid");
    let mut slab = Slab::new(|_|4.0,[4.0,4.0]);
    slab.conductivity[0][0][0]=2.5;
    assert!(matches!(slab.range(),Err(TetError::Unsupported(_))));
    let slab = Slab::new(|_|4.0,[4.0,4.0]); let p=slab.problem(); let b=FluxBudget::default();
    assert!(matches!(p.slab_temperature_range_two_sided(&[300.0;8],3,b,||true),Err(TetError::Invalid(_))));
    assert!(p.slab_temperature_range_two_sided(&[300.0;7],0,b,||true).is_err());
    assert!(matches!(p.slab_temperature_range_two_sided(&[300.0;8],0,
        FluxBudget {max_cells:5,max_iterations:0},||true),Err(TetError::Budget)));
    assert!(matches!(p.slab_temperature_range_two_sided(&[300.0;8],0,b,||false),Err(TetError::Cancelled)));
    let mut count=0;
    assert!(matches!(p.slab_temperature_range_two_sided(&[300.0;8],0,b,|| {
        count+=1; count<30
    }),Err(TetError::Cancelled)));
    let slab = Slab::new(|_|4.0,[f64::MIN_POSITIVE/10.0;2]);
    assert!(matches!(slab.range(),Err(TetError::Unbounded)));
}
