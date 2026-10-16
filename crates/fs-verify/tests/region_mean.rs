#![cfg(feature = "certified-speculation")]
//! Exact regional functionals; selection must not carve out a different PDE.
use std::collections::BTreeMap;
use fs_verify::tet::{AffineSourceTetProblem, BoundaryCondition, BoundaryFace,
    ConductivityTensor, FluxBudget, RegionMeanSelection, TetError, affine_source_mean_bound};

struct Fixture {
    points: Vec<[f64; 3]>, tets: Vec<[usize; 4]>, tensors: Vec<ConductivityTensor>,
    source: Vec<[f64; 4]>, boundary: Vec<BoundaryFace>, candidate: Vec<f64>,
}
impl Fixture {
    fn problem(&self) -> AffineSourceTetProblem<'_> {
        AffineSourceTetProblem { vertices: &self.points, tets: &self.tets,
            conductivity: &self.tensors, source: &self.source, boundary: &self.boundary }
    }
    fn left(&self, end: f64) -> RegionMeanSelection {
        let cells: Vec<_> = self.tets.iter().enumerate().filter_map(|(e, tet)|
            tet.iter().all(|&v| self.points[v][0] <= end).then_some(e)).collect();
        RegionMeanSelection::new(self.tets.len(), &cells, FluxBudget::default(), || true).unwrap()
    }
}
fn slab(xs: &[f64], heated: bool) -> Fixture {
    let mut points = Vec::new();
    for z in [0.0, 1.0] { for y in [0.0, 1.0] { for &x in xs { points.push([x,y,z]); } } }
    let index = |p: [usize; 3]| p[0] + xs.len() * (p[1] + 2*p[2]);
    let mut tets = Vec::new();
    for x in 0..xs.len()-1 {
        for permutation in [[0,1,2], [0,2,1], [1,0,2], [1,2,0], [2,0,1], [2,1,0]] {
            let mut p = [x,0,0]; let mut tet = [index(p);4];
            for (i,d) in permutation.into_iter().enumerate() { p[d]+=1; tet[i+1]=index(p); }
            tets.push(tet);
        }
    }
    let candidate: Vec<_> = points.iter().map(|p|
        if heated { 300.0 + p[0] - p[0]*p[0]*p[0] } else { 300.0 + 2.0*p[0] }).collect();
    let source = tets.iter().map(|tet| tet.map(|v|
        if heated { 6.0*points[v][0] } else { 0.0 })).collect();
    let mut incidence = BTreeMap::<[usize;3], usize>::new();
    for tet in &tets { for i in 0..4 {
        let mut face: Vec<_> = tet.iter().enumerate().filter_map(|(j,&v)|(i!=j).then_some(v)).collect();
        face.sort_unstable(); *incidence.entry([face[0],face[1],face[2]]).or_default()+=1;
    } }
    let boundary = incidence.into_iter().filter(|(_,n)|*n==1).map(|(vertices,_)| {
        let end = vertices.iter().all(|&v| points[v][0]==0.0)
            || vertices.iter().all(|&v| points[v][0]==1.0);
        BoundaryFace { vertices, condition: if end {
            BoundaryCondition::Dirichlet(vertices.map(|v|candidate[v]))
        } else { BoundaryCondition::Neumann(0.0) } }
    }).collect();
    Fixture { tensors: vec![[[1.0,0.0,0.0],[0.0,2.0,0.5],[0.0,0.5,3.0]];tets.len()],
        points,tets,source,boundary,candidate }
}
fn contains(interval: fs_verify::interval::Iv, truth: f64) {
    assert!(interval.lo<=truth && truth<=interval.hi, "{interval:?} excludes {truth}");
}

#[test]
fn region_volume_not_whole_volume_or_equal_cell_averaging() {
    let f=slab(&[0.0,0.125,0.5,1.0],false);
    let selected=f.left(0.5);
    let bound=f.problem().region_mean_bound(&f.candidate,&vec![0.0;f.points.len()],
        &selected,FluxBudget::default(),||true).unwrap();
    contains(bound.region_volume,0.5); contains(bound.integral.domain_volume,1.0);
    contains(bound.candidate_mean,300.5); contains(bound.enclosure,300.5);
    assert!(bound.enclosure.hi-bound.enclosure.lo<1e-5);
    let all=RegionMeanSelection::new(f.tets.len(),&(0..f.tets.len()).collect::<Vec<_>>(),
        FluxBudget::default(),||true).unwrap();
    let whole=f.problem().region_mean_bound(&f.candidate,&vec![0.0;f.points.len()],
        &all,FluxBudget::default(),||true).unwrap();
    contains(whole.enclosure,301.0);
    assert!(bound.enclosure.hi<whole.enclosure.lo);
}

#[test]
fn region_supported_dual_and_inexact_fields_keep_the_correct_continuum_goal() {
    let f=slab(&[0.0,0.25,0.5,0.75,1.0],true);
    let selected=f.left(0.5);
    // Exact dual for integral over 0<=x<=1/2, with kx=1 and zero end values.
    let dual:Vec<_>=f.points.iter().map(|p| if p[0]<=0.5 {
        -0.5*p[0]*p[0]+0.375*p[0]
    } else {0.125*(1.0-p[0])}).collect();
    for candidate in [&f.candidate, &vec![300.0;f.points.len()]] {
        let bound=f.problem().region_mean_bound(candidate,&dual,&selected,FluxBudget::default(),||true).unwrap();
        contains(bound.enclosure,300.21875);
        assert_eq!(bound.integral.dual.cell_majorant_squared_upper.len(),f.tets.len());
        assert!(bound.integral.dual.energy_error_upper>0.0);
    }
    // A wrong/unfinished dual is not trusted: its missing region-source response
    // is included in the same full-domain bound rather than changing the goal.
    let bound=f.problem().region_mean_bound(&vec![300.0;f.points.len()],&vec![0.0;f.points.len()],
        &selected,FluxBudget::default(),||true).unwrap();
    contains(bound.enclosure,300.21875);
    assert!(bound.enclosure.hi-bound.enclosure.lo>0.4);
}

#[test]
fn selecting_all_cells_preserves_existing_mean_bits_and_reordering_replays() {
    let f=slab(&[0.0,0.25,0.5,0.75,1.0],true);
    let dual:Vec<_>=f.points.iter().map(|p|0.5*p[0]*(1.0-p[0])).collect();
    let indices:Vec<_>=(0..f.tets.len()).rev().collect();
    let selection=RegionMeanSelection::new(f.tets.len(),&indices,FluxBudget::default(),||true).unwrap();
    let region=f.problem().region_mean_bound(&f.candidate,&dual,&selection,FluxBudget::default(),||true).unwrap();
    let whole=affine_source_mean_bound(&f.problem(),&f.candidate,&dual,FluxBudget::default(),||true).unwrap();
    assert_eq!(region.enclosure,whole.enclosure);
    assert_eq!(region.candidate_mean,whole.candidate_mean);
    assert_eq!(region.region_volume,whole.integral.domain_volume);
    assert_eq!(selection.cells(),&(0..f.tets.len()).collect::<Vec<_>>());
}

#[test]
fn bad_selection_and_cancelled_or_unresolved_domains_return_no_bound() {
    let budget=FluxBudget::default();
    for cells in [vec![],vec![0,0],vec![2],vec![usize::MAX]] {
        assert!(RegionMeanSelection::new(2,&cells,budget,||true).is_err());
    }
    assert_eq!(RegionMeanSelection::new(2,&[0],budget,||false),Err(TetError::Cancelled));
    let f=slab(&[0.0,0.5,1.0],false);
    let wrong=RegionMeanSelection::new(1,&[0],budget,||true).unwrap();
    assert!(f.problem().region_mean_bound(&f.candidate,&vec![0.0;f.points.len()],&wrong,budget,||true).is_err());
    let selected=f.left(0.5);
    assert!(matches!(f.problem().region_mean_bound(&f.candidate,&vec![0.0;f.points.len()],
        &selected,budget,||false),Err(TetError::Cancelled)));
    let small=FluxBudget {max_cells:1,..budget};
    assert!(matches!(f.problem().region_mean_bound(&f.candidate,&vec![0.0;f.points.len()],
        &selected,small,||true),Err(TetError::Budget)));
    let mut invalid=f; invalid.tensors[invalid.tets.len()-1][0][0]=-1.0;
    assert!(invalid.problem().region_mean_bound(&invalid.candidate,&vec![0.0;invalid.points.len()],
        &selected,budget,||true).is_err(),"unselected bad material must not be hidden");
}
