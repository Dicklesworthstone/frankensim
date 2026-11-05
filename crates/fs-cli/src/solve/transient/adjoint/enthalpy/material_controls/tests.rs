//! Live contact support must not weaken the original fixed-owner contract.
use super::*;
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_conduction::{ScalarField, ThermalBoundaryBuilder, fixtures::box_grid};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};

#[test]
fn fixed_material_controls_reject_an_equal_but_replaced_interface_owner() {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(&gate, arena,
            StreamKey { seed: 71, kernel_id: 825, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic);
        let (complex, positions) = box_grid([1, 1, 1], [1.0, 1.0, 1.0]);
        let mesh = ConductionMesh::new(complex, positions).unwrap();
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .adiabatic_remainder().finish().unwrap();
        let material = ConductivityModel::isotropic_declared(2.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        let problem = ConductionProblem { mesh: &mesh, boundary: &boundary,
            material: &material, element_materials: None, source: &source };
        // Even the empty contact layout has an exact original operator owner.
        // Equal data in a distinct allocation must not substitute for it.
        let interfaces = ThermalInterfaces::new(&mesh, &boundary, Vec::new()).unwrap();
        let replacement = interfaces.clone();
        assert_eq!(interfaces, replacement);
        let temperature = vec![300.0; mesh.vertex_count()];
        let dual = vec![0.0; mesh.vertex_count()];
        let mut controls = Controls::new(&cx, &ProjectSpec::default(), problem,
            Some(&interfaces), &vec![0; mesh.element_count()], &BTreeMap::new(),
            &temperature, None).unwrap();
        for substituted in [Some(&replacement), None] {
            assert!(controls.add_capacity(&cx, problem, substituted,
                &temperature, &dual, None).is_err());
            assert!(controls.fragment(0).is_ok(), "a refusal must not advance the prefix");
        }
        controls.add_capacity(&cx, problem, Some(&interfaces),
            &temperature, &dual, None).unwrap();
        assert!(controls.fragment(1).is_ok());
    });
}
