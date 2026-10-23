use super::*;
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_rep_mesh::TetComplex;

const DATA: &str = "FS_SCAN_MATERIALS_V1
material 20 2 6 150
phase_knot 20 -200 280 0
phase_knot 20 100 310 0
phase_knot 20 200 310 1
phase_knot 20 900 380 1
material 10 4 12 200
phase_knot 10 -200 300 0
phase_knot 10 100 330 0
phase_knot 10 300 330 1
phase_knot 10 900 390 1
element 0 20
element 1 10
initial 1 160
";
const CAP: ScanMaterialBudget = ScanMaterialBudget {
    max_bytes: 4096, max_materials: 2, max_phase_knots: 8,
};
fn mesh() -> ConductionMesh {
    ConductionMesh::new(TetComplex::from_tets(8, vec![[0,1,2,3], [4,5,6,7]]),
        vec![[0.,0.,0.], [1.,0.,0.], [0.,1.,0.], [0.,0.,1.],
            [3.,0.,0.], [4.,0.,0.], [3.,1.,0.], [3.,0.,1.]]).unwrap()
}
fn with_cx(f: impl FnOnce(&CancelGate, &Cx<'_>)) {
    let gate = CancelGate::new();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        f(&gate, &Cx::new(&gate, arena,
            StreamKey { seed: 76, kernel_id: 51, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic));
    });
}
#[test]
fn material_ids_bind_native_masses_charts_and_nonuniform_initial_state() {
    with_cx(|_, cx| {
        let mesh = mesh();
        let map = ScanMaterialMap::from_text(cx, &mesh, DATA, CAP).unwrap();
        assert!(std::ptr::eq(map.mesh(), &mesh));
        assert_eq!(map.material_count(), 2);
        assert_eq!(map.element_material_ids(), &[1,0]);
        assert_eq!(map.initial_specific_enthalpy_j_kg(), &[150.,160.,150.,150.,200.,200.,200.,200.]);
        assert!((map.reference_mass_kg() - 3.).abs() < 1e-14);
        let references = map.references();
        assert_eq!(references[0].reference_density_kg_m3, 12.);
        assert_eq!(references[1].reference_density_kg_m3, 6.);
        let prepared = HeterogeneousEnthalpyBackwardEuler::new(cx, &mesh, &references,
            map.element_material_ids(), EnthalpyBudget { max_vertices: 8, max_elements: 2 }).unwrap();
        assert_eq!(prepared.vertex_material_ids(), &[1,1,1,1,0,0,0,0]);
        for (&mass, expected) in prepared.reference_nodal_masses_kg().iter()
            .zip([0.25,0.25,0.25,0.25,0.5,0.5,0.5,0.5]) {
            assert!((mass-expected).abs() < 1e-14);
        }
        let changed = ScanMaterialMap::from_text(cx, &mesh,
            &DATA.replace("material 10 4 12 200", "material 10 4 6 200"), CAP).unwrap();
        assert!((changed.reference_mass_kg()-2.).abs() < 1e-14);
        assert_ne!(changed.identity(), map.identity());
    });
}
#[test]
fn malformed_or_underbudget_material_maps_cannot_silently_fall_back() {
    with_cx(|gate, cx| {
        let mesh = mesh();
        for bad in [DATA.replace("element 1 10\n", ""),
            DATA.replace("element 1 10", "element 1 30"),
            DATA.replace("initial 1 160", "initial 1 1000"),
            DATA.replace("initial 1 160", "initial 8 160"),
            DATA.replace("material 20 2 6 150", "material 20 NaN 6 150"),
            DATA.replace("phase_knot 20 200 310 1", "phase_knot 20 50 310 1"),
            DATA.replace("phase_knot 10 300 330 1", "phase_knot 10 300 330 100"),
            format!("{DATA}element 0 20\n"), format!("{DATA}initial 1 160\n"),
            format!("{DATA}unsupported 1\n"), DATA.replace("element 1 10", "element 1 20")]
        { assert!(ScanMaterialMap::from_text(cx, &mesh, &bad, CAP).is_err(), "{bad}"); }
        for cap in [ScanMaterialBudget { max_bytes: DATA.len()-1, ..CAP },
            ScanMaterialBudget { max_materials: 1, ..CAP },
            ScanMaterialBudget { max_phase_knots: 7, ..CAP }]
        { assert!(ScanMaterialMap::from_text(cx, &mesh, DATA, cap).is_err()); }
        gate.request();
        assert!(ScanMaterialMap::from_text(cx, &mesh, DATA, CAP).is_err());
    });
}
#[test]
fn conflicting_caloric_ids_at_a_shared_vertex_remain_a_native_refusal() {
    with_cx(|_, cx| {
        let mesh = ConductionMesh::new(TetComplex::from_tets(7, vec![[0,1,2,3], [0,4,6,5]]),
            vec![[0.,0.,0.], [1.,0.,0.], [0.,1.,0.], [0.,0.,1.],
                [-1.,0.,0.], [0.,-1.,0.], [0.,0.,-1.]]).unwrap();
        let error = ScanMaterialMap::from_text(cx, &mesh, DATA, CAP).err().unwrap().to_string();
        assert!(error.contains("vertex 0 belongs to materials"), "{error}");
    });
}
