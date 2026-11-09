//! Native constitutive refresh without reconstructing geometry or constraints.
use std::cell::Cell;
use std::ops::ControlFlow;

use fs_blake3::ContentHash;
use fs_cutfem::elastic3::adaptive::AdaptiveElasticity3;
use fs_cutfem::elastic3::{CutElasticity3, ElasticMaterialField3, ElasticityError3,
    ElasticityFormulation3, ElasticityOptions3, MaterialRefreshLimits3};
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use fs_material::elastic::acoustic::ElasticAcousticState;
use fs_solver::LinearOp;

struct Slab<'a>(&'a Cell<usize>);
impl Slab<'_> {
    fn touch(&self) { self.0.set(self.0.get() + 1); }
}
impl CutSdf3 for Slab<'_> {
    fn value(&self, p: [f64; 3]) -> f64 { self.touch(); p[2] - 0.73 }
    fn enclose(&self, lo: [f64; 3], hi: [f64; 3]) -> Interval {
        self.touch(); Interval::new(lo[2] - 0.73, hi[2] - 0.73)
    }
    fn derivative_enclose(&self, _: [f64; 3], _: [f64; 3], axis: HeightAxis) -> Interval {
        self.touch(); let value = if axis == HeightAxis::Z { 1.0 } else { 0.0 };
        Interval::new(value, value)
    }
}
fn proceed() -> ControlFlow<()> { ControlFlow::Continue(()) }
fn original() -> ElasticAcousticState {
    ElasticAcousticState::from_bulk_and_shear(3.0, 20.0, 8.0, ContentHash([41; 32]), 0.1).unwrap()
}
fn updated(bounds: HexCell, formulation: ElasticityFormulation3) -> ElasticAcousticState {
    if formulation == ElasticityFormulation3::BBar {
        return ElasticAcousticState::from_bulk_and_shear(3.0,
            120.0 + 10.0 * bounds.lo()[0], 2.0 + bounds.lo()[2], ContentHash([42; 32]), 0.07).unwrap();
    }
    let mut tensor = std::array::from_fn(|i|
        std::array::from_fn(|j| if i == j { 30.0 + 5.0 * i as f64 + bounds.lo()[0] } else { 0.0 }));
    tensor[0][1] = 4.0; tensor[1][0] = 4.0;
    tensor[0][3] = 1.5; tensor[3][0] = 1.5;
    tensor[2][5] = -1.1; tensor[5][2] = -1.1;
    ElasticAcousticState::try_new_mandel(tensor, 3.0, ContentHash([43; 32]), 0.07).unwrap()
}
fn geometry<T>(action: impl FnOnce(&mut QuadratureControl3<'_>) -> T) -> T {
    let mut callback = |_| proceed();
    let mut control = QuadratureControl3::new(
        QuadratureOptions3 { depth: 0, ..Default::default() }, &mut callback).unwrap();
    action(&mut control)
}
fn raw(field: ElasticMaterialField3<'_>, formulation: ElasticityFormulation3,
    calls: &Cell<usize>) -> CutElasticity3 {
    geometry(|control| CutElasticity3::build_with_material_field_and_embedded_dirichlet(
        HexCell::try_new([0.0; 3], [1.0; 3]).unwrap(), [2, 1, 2], &Slab(calls), field,
        &|p| { calls.set(calls.get() + 1); p[0] == 0.0 },
        &|_, _| { calls.set(calls.get() + 1); true },
        ElasticityOptions3 { formulation, ..Default::default() },
        Default::default(), Default::default(), control).unwrap())
}
fn adaptive(field: ElasticMaterialField3<'_>, formulation: ElasticityFormulation3,
    calls: &Cell<usize>) -> AdaptiveElasticity3 {
    let tree = Octree3::uniform(1, 4, 4096).unwrap();
    let leaf = *tree.leaves().iter().find(|cell| cell.index() == [0, 0, 1]).unwrap();
    let tree = tree.refined(&[leaf], proceed).unwrap();
    geometry(|control| AdaptiveElasticity3::build_with_material_field_and_embedded_dirichlet(
        HexCell::try_new([0.0; 3], [1.0; 3]).unwrap(), &tree, &Slab(calls), field,
        &|p| { calls.set(calls.get() + 1); p[0] == 0.0 },
        &|_, _| { calls.set(calls.get() + 1); true },
        ElasticityOptions3 { formulation, ..Default::default() },
        Default::default(), Default::default(), control).unwrap())
}
fn action(op: &impl LinearOp, input: &[f64]) -> Vec<f64> {
    let mut result = vec![0.0; op.n()]; op.apply(input, &mut result); result
}
fn assert_same_operator(a: &impl LinearOp, b: &impl LinearOp) {
    assert_eq!(a.n(), b.n());
    let mut basis = vec![0.0; a.n()];
    for i in 0..a.n() {
        basis[i] = 1.0;
        assert_eq!(action(a, &basis), action(b, &basis), "native column {i}");
        basis[i] = 0.0;
    }
}
fn strain(_: usize, p: [f64; 3]) -> Result<[f64; 6], ElasticityError3> {
    Ok([0.001 * (1.0 + p[0]), -0.0003, 0.0002, 0.0001 * p[2], 0.0, 0.0])
}

#[test]
fn cartesian_refresh_matches_fresh_bulk_nitsche_ghost_and_thermal_loads() {
    for formulation in [ElasticityFormulation3::Displacement, ElasticityFormulation3::BBar] {
        let calls = Cell::new(0);
        let mut op = raw(ElasticMaterialField3::Uniform(&original()), formulation, &calls);
        let scales: Vec<_> = (0..op.cells()).map(|i| 0.3 + 0.1 * i as f64).collect();
        op.set_scales(&scales).unwrap();
        let nodes = op.nodes().as_ptr();
        let inventory = op.reference_cell_masses_kg().unwrap().to_vec();
        let bounds = op.volume_bounds();
        let geometry_calls = calls.get();
        let material_calls = Cell::new(0);
        let query = |bounds| { material_calls.set(material_calls.get() + 1); Ok(updated(bounds, formulation)) };
        op.refresh_material_field(ElasticMaterialField3::ByCell(&query), Default::default(), proceed).unwrap();
        assert_eq!(material_calls.get(), op.cells());
        assert_eq!(calls.get(), geometry_calls, "no SDF, support or clamp callback may repeat");
        assert_eq!(op.nodes().as_ptr(), nodes, "original geometry storage must survive");
        assert_eq!(op.reference_cell_masses_kg().unwrap(), inventory);
        assert_eq!(op.volume_bounds(), bounds);
        assert_eq!(op.scales(), scales);
        let mut fresh = raw(ElasticMaterialField3::ByCell(&|bounds| Ok(updated(bounds, formulation))),
            formulation, &calls);
        fresh.set_scales(&scales).unwrap();
        assert_same_operator(&op, &fresh);
        for cell in 0..op.cells() { assert_eq!(op.cell_material_state(cell), fresh.cell_material_state(cell)); }
        let motion = |p: [f64; 3], _: [f64; 3]| [0.0001 * p[1], 0.0, 0.0001];
        assert_eq!(op.prescribed_displacement_load(&motion, proceed).unwrap(),
            fresh.prescribed_displacement_load(&motion, proceed).unwrap());
        let field = op.prepare_eigenstrain(&strain, 100_000, proceed).unwrap();
        let other = fresh.prepare_eigenstrain(&strain, 100_000, proceed).unwrap();
        assert_eq!(op.eigenstrain_load(&field, proceed).unwrap(), fresh.eigenstrain_load(&other, proceed).unwrap());
        let u = vec![0.0; op.n()];
        assert_eq!(op.eigenstrain_response(&field, &u, 0.1, 100_000, proceed).unwrap(),
            fresh.eigenstrain_response(&other, &u, 0.1, 100_000, proceed).unwrap());
    }
}

#[test]
fn adaptive_refresh_keeps_exact_hanging_constraints_and_matches_fresh_assembly() {
    for formulation in [ElasticityFormulation3::Displacement, ElasticityFormulation3::BBar] {
        let calls = Cell::new(0);
        let mut op = adaptive(ElasticMaterialField3::Uniform(&original()), formulation, &calls);
        assert!(op.physical_nodes().len() > op.nodes().len());
        let scales: Vec<_> = (0..op.cells()).map(|i| 0.3 + 0.02 * (i % 9) as f64).collect();
        op.set_scales(&scales).unwrap();
        let masters = op.nodes().as_ptr(); let physical = op.physical_nodes().as_ptr();
        let leaves = op.leaves().as_ptr(); let edges = op.filter_edges().as_ptr();
        let geometry_calls = calls.get();
        op.refresh_material_field(ElasticMaterialField3::ByCell(&|bounds| Ok(updated(bounds, formulation))),
            Default::default(), proceed).unwrap();
        assert_eq!(calls.get(), geometry_calls);
        assert_eq!(op.nodes().as_ptr(), masters); assert_eq!(op.physical_nodes().as_ptr(), physical);
        assert_eq!(op.leaves().as_ptr(), leaves); assert_eq!(op.filter_edges().as_ptr(), edges);
        assert_eq!(op.scales(), scales);
        let mut fresh = adaptive(ElasticMaterialField3::ByCell(&|bounds| Ok(updated(bounds, formulation))),
            formulation, &calls);
        fresh.set_scales(&scales).unwrap();
        assert_same_operator(&op, &fresh);
        let field = op.prepare_eigenstrain(&strain, 100_000, proceed).unwrap();
        let other = fresh.prepare_eigenstrain(&strain, 100_000, proceed).unwrap();
        assert_eq!(op.eigenstrain_load(&field, proceed).unwrap(), fresh.eigenstrain_load(&other, proceed).unwrap());
        let coefficients: Vec<_> = op.nodes().iter().flat_map(|p| [0.001 * p[0], 0.0, 0.0]).collect();
        assert_eq!(op.physical_displacements(&coefficients).unwrap(), fresh.physical_displacements(&coefficients).unwrap());
    }
}

#[test]
fn only_success_rotates_native_preparations_and_reference_mass_remains_unchanged() {
    let calls = Cell::new(0);
    let mut op = raw(ElasticMaterialField3::Uniform(&original()), ElasticityFormulation3::Displacement, &calls);
    let field = op.prepare_eigenstrain(&strain, 100_000, proceed).unwrap();
    let mass = op.prepare_material_mass(100_000, proceed).unwrap();
    let body = op.prepare_material_body_load(&|_| [0.0, -9.81, 0.0], 100_000, proceed).unwrap();
    let occupancy = vec![0.61; op.cells()];
    let mass_before = mass.diagonal(&occupancy, proceed).unwrap();
    let inventory = op.reference_cell_masses_kg().unwrap().to_vec();
    let rejected = ElasticAcousticState::from_bulk_and_shear(3.1, 120.0, 2.0,
        ContentHash([55; 32]), 0.1).unwrap();
    assert!(op.refresh_material_field(ElasticMaterialField3::Uniform(&rejected),
        Default::default(), proceed).is_err());
    op.validate_eigenstrain(&field).unwrap();
    op.validate_material_mass(&mass).unwrap();
    op.validate_material_body_load(&body).unwrap();
    op.refresh_material_field(ElasticMaterialField3::Uniform(&updated(
        HexCell::try_new([0.0; 3], [1.0; 3]).unwrap(), ElasticityFormulation3::Displacement)),
        Default::default(), proceed).unwrap();
    assert!(op.validate_eigenstrain(&field).is_err());
    assert!(op.validate_material_mass(&mass).is_err());
    assert!(op.validate_material_body_load(&body).is_err());
    assert_eq!(op.reference_cell_masses_kg().unwrap(), inventory);
    let new_mass = op.prepare_material_mass(100_000, proceed).unwrap();
    assert_eq!(new_mass.diagonal(&occupancy, proceed).unwrap(), mass_before);
    op.validate_material_mass(&new_mass).unwrap();
}

#[test]
fn allowance_and_material_refusals_happen_before_geometry_or_operator_changes() {
    let calls = Cell::new(0);
    let mut op = raw(ElasticMaterialField3::Uniform(&original()), ElasticityFormulation3::BBar, &calls);
    let u: Vec<_> = (0..op.n()).map(|i| 0.001 * (i % 7) as f64).collect();
    let before = action(&op, &u);
    let field = op.prepare_eigenstrain(&strain, 100_000, proceed).unwrap();
    let geometry_calls = calls.get(); let queries = Cell::new(0);
    let query = |bounds| { queries.set(queries.get() + 1); Ok(updated(bounds, ElasticityFormulation3::BBar)) };
    for limits in [
        MaterialRefreshLimits3 { max_cells: op.cells() - 1, ..Default::default() },
        MaterialRefreshLimits3 { max_points: 0, ..Default::default() },
        MaterialRefreshLimits3 { max_scratch_bytes: 1, ..Default::default() },
    ] {
        assert!(op.refresh_material_field(ElasticMaterialField3::ByCell(&query), limits, proceed).is_err());
        assert_eq!(queries.get(), 0, "preflight must precede material callbacks");
        assert_eq!(action(&op, &u), before); op.validate_eigenstrain(&field).unwrap();
    }
    // A generic tensor, even if well-conditioned, must not acquire a B-bar fit.
    assert!(op.refresh_material_field(ElasticMaterialField3::ByCell(&|bounds|
        Ok(updated(bounds, ElasticityFormulation3::Displacement))), Default::default(), proceed).is_err());
    assert_eq!(action(&op, &u), before); op.validate_eigenstrain(&field).unwrap();
    // Refuse late in a spatial source so already-staged cells cannot escape.
    let queries = Cell::new(0);
    let cells = op.cells();
    let bad = |bounds| {
        queries.set(queries.get() + 1);
        if queries.get() == cells { Err(ElasticityError3::Invalid("missing thermal source")) }
        else { Ok(updated(bounds, ElasticityFormulation3::BBar)) }
    };
    assert!(op.refresh_material_field(ElasticMaterialField3::ByCell(&bad), Default::default(), proceed).is_err());
    assert_eq!(queries.get(), cells);
    assert_eq!(action(&op, &u), before); op.validate_eigenstrain(&field).unwrap();
    assert_eq!(calls.get(), geometry_calls);
}

#[test]
fn cancellation_including_the_final_publish_poll_keeps_all_old_values_and_owners() {
    let calls = Cell::new(0);
    let material = original();
    let query = |bounds| Ok(updated(bounds, ElasticityFormulation3::Displacement));
    let mut successful = raw(ElasticMaterialField3::Uniform(&material), ElasticityFormulation3::Displacement, &calls);
    let mut polls = 0usize;
    successful.refresh_material_field(ElasticMaterialField3::ByCell(&query), Default::default(), || {
        polls += 1; proceed()
    }).unwrap();
    assert!(polls > successful.cells() * 3);
    for stop in [1, polls / 2, polls - 1, polls] {
        let mut op = raw(ElasticMaterialField3::Uniform(&material), ElasticityFormulation3::Displacement, &calls);
        let scales: Vec<_> = (0..op.cells()).map(|i| 0.4 + 0.1 * i as f64).collect();
        op.set_scales(&scales).unwrap();
        let field = op.prepare_eigenstrain(&strain, 100_000, proceed).unwrap();
        let mass = op.prepare_material_mass(100_000, proceed).unwrap();
        let body = op.prepare_material_body_load(&|_| [0.0, -1.0, 0.0], 100_000, proceed).unwrap();
        let u: Vec<_> = (0..op.n()).map(|i| 0.001 * (i % 7) as f64).collect();
        let applied = action(&op, &u); let inventory = op.reference_cell_masses_kg().unwrap().to_vec();
        let mut count = 0usize;
        assert_eq!(op.refresh_material_field(ElasticMaterialField3::ByCell(&query), Default::default(), || {
            count += 1; if count == stop { ControlFlow::Break(()) } else { proceed() }
        }), Err(ElasticityError3::Cancelled));
        assert_eq!(action(&op, &u), applied); assert_eq!(op.scales(), scales);
        assert_eq!(op.reference_cell_masses_kg().unwrap(), inventory);
        for index in 0..op.cells() { assert_eq!(op.cell_material_state(index), Some(&material)); }
        op.validate_eigenstrain(&field).unwrap(); op.validate_material_mass(&mass).unwrap();
        op.validate_material_body_load(&body).unwrap();
    }
}

#[test]
fn legacy_without_reference_inventory_cannot_be_thermally_rematerialized() {
    let calls = Cell::new(0);
    let mut op = geometry(|control| CutElasticity3::build(
        HexCell::try_new([0.0; 3], [1.0; 3]).unwrap(), [1; 3], &Slab(&calls),
        &IsotropicElastic::new(20.0, 0.25, 0.1).unwrap(), &|p| p[0] == 0.0,
        Default::default(), control).unwrap());
    assert!(matches!(op.refresh_material_field(ElasticMaterialField3::Uniform(&original()),
        Default::default(), proceed), Err(ElasticityError3::Invalid("material refresh requires original native reference density"))));
}
