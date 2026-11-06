//! G1/G3: sourced failed bands continue on actual split P1 geometry.
//! Prescribed manufactured load paths verify state/energy transfer, not a
//! calibrated material or a quasi-static stability claim through softening.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::brittle::source::*;
use fs_material::cohesive::source::FRACTURE_ENERGY_PER_AREA_DIMS;
use fs_material::state_point::{
    MaterialPropertySelection, POISSON_RATIO_PROPERTY, YOUNG_MODULUS_PROPERTY,
};
use fs_qty::{Dims, Pressure};
use fs_solid::brittle::TetCrackBandGeometry;
use fs_solid::brittle::solid::fracture::{BrittleTetFractureControls, BrittleTetFractureError};
use fs_solid::brittle::solid::{BrittleTetControls, BrittleTetError, BrittleTetSolid};
use fs_solid::brittle::topology::{FractureTopologyBudget, InteriorSeam, TetFractureTopology};
use fs_solid::linear3::TetAssemblyBudget;
use fs_solid::loads3::TetLoadGeometry3;
use fs_solid::plastic3::IncrementSettings;

const YOUNG: f64 = 1024.0;
const GC: f64 = 1.0 / 1024.0;
const FAILED_STRAIN: [f64; 3] = [1.0 / 1024.0, 1.0 / 2048.0, 1.0 / 128.0];
const NODES: [[f64; 3]; 5] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, -1.0],
    [0.0, 0.0, 1.0],
];
const TETS: [[usize; 4]; 2] = [[0, 2, 1, 3], [0, 1, 2, 4]];

fn with_cx<R>(cancelled: bool, run: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    if cancelled {
        gate.request();
    }
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        run(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 931,
                kernel_id: 931,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

fn budget() -> TetAssemblyBudget {
    TetAssemblyBudget {
        maximum_nodes: 12,
        maximum_tetrahedra: 4,
        maximum_free_dofs: 36,
        minimum_scaled_jacobian: 1e-10,
    }
}

fn controls() -> BrittleTetControls {
    BrittleTetControls {
        material_energy_tolerance_j: 1e-10,
        equilibrium_work_tolerance_j: 1e-10,
        // The native owner separately reports this endpoint quadrature error.
        // The prescribed rupture path intentionally spans a softening event.
        maximum_endpoint_work_defect_j: 0.05,
    }
}

fn release_controls() -> BrittleTetFractureControls {
    BrittleTetFractureControls {
        normal_alignment_tolerance: 1e-13,
        traction_tolerance_pa: 1e-10,
        energy_tolerance_j: 1e-12,
        opening_tolerance_m: 1e-13,
    }
}

fn source_card() -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, YOUNG),
        (POISSON_RATIO_PROPERTY, Dims::NONE, 0.0),
        (BRITTLE_TENSILE_STRENGTH_PROPERTY, Pressure::DIMS, 1.0),
        (BRITTLE_SHEAR_STRENGTH_PROPERTY, Pressure::DIMS, 1.25),
        (
            BRITTLE_MODE_I_FRACTURE_ENERGY_PROPERTY,
            FRACTURE_ENERGY_PER_AREA_DIMS,
            GC,
        ),
        (
            BRITTLE_MODE_II_FRACTURE_ENERGY_PROPERTY,
            FRACTURE_ENERGY_PER_AREA_DIMS,
            2.0 * GC,
        ),
        (BRITTLE_BK_EXPONENT_PROPERTY, Dims::NONE, 1.5),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained(),
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured native brittle fragment fixture".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured brittle solid".into(),
            phase: "solid".into(),
            process: "fracture transfer verification".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn materials(
    card: &MaterialCard,
    geometry: &TetLoadGeometry3<'_>,
    normal: [f64; 3],
) -> Vec<ResolvedBrittleCrackBand> {
    (0..geometry.tetrahedra().len())
        .map(|element| {
            let band = TetCrackBandGeometry::from_element(geometry, element, normal).unwrap();
            resolve_brittle_crack_band(
                card,
                &QueryPoint::new(),
                MaterialPropertySelection::SingleClaimOnly,
                band.material_settings(0.03, 1e-9),
            )
            .unwrap()
        })
        .collect()
}

fn prescribed(u: &[[f64; 3]]) -> Vec<(usize, f64)> {
    u.as_flattened().iter().copied().enumerate().collect()
}

fn affine(nodes: &[[f64; 3]], diagonal: [f64; 3]) -> Vec<[f64; 3]> {
    nodes
        .iter()
        .map(|p| core::array::from_fn(|i| diagonal[i] * p[i]))
        .collect()
}

fn loaded<'a>(
    geometry: &'a TetLoadGeometry3<'a>,
    sources: &'a [ResolvedBrittleCrackBand],
    diagonal: [f64; 3],
    cx: &Cx<'_>,
) -> BrittleTetSolid<'a> {
    BrittleTetSolid::try_new(geometry, sources, budget(), cx)
        .unwrap()
        .advance(
            &prescribed(&affine(geometry.nodes_m(), diagonal)),
            &[],
            IncrementSettings::default(),
            controls(),
            cx,
        )
        .unwrap()
        .solid
}

fn topology<'g>(
    geometry: &'g TetLoadGeometry3<'_>,
    face: [usize; 3],
    side_a: usize,
    id: u64,
    cx: &Cx<'_>,
) -> TetFractureTopology<'g> {
    TetFractureTopology::new(
        geometry,
        &[InteriorSeam {
            id,
            nodes: face,
            side_a_element: side_a,
        }],
        FractureTopologyBudget {
            maximum_output_nodes: 12,
            maximum_face_entries: 16,
            maximum_seams: 1,
        },
        cx,
    )
    .unwrap()
}

fn exact_material_transfer(old: &BrittleTetSolid<'_>, next: &BrittleTetSolid<'_>) {
    assert!(core::ptr::eq(
        old.materials().as_ptr(),
        next.materials().as_ptr()
    ));
    assert_eq!(
        old.solution().element_states,
        next.solution().element_states
    );
    for (a, b) in old
        .solution()
        .element_strain
        .iter()
        .zip(&next.solution().element_strain)
    {
        assert_eq!(a.map(f64::to_bits), b.map(f64::to_bits));
    }
    for (a, b) in old
        .solution()
        .element_stress
        .iter()
        .zip(&next.solution().element_stress)
    {
        assert_eq!(a.map(f64::to_bits), b.map(f64::to_bits));
    }
    for (a, b) in old
        .solution()
        .element_states
        .iter()
        .zip(&next.solution().element_states)
    {
        assert!(core::ptr::eq(
            a.cohesive_point().law(),
            b.cohesive_point().law()
        ));
        assert_eq!(
            a.cohesive_point().jump_m().map(f64::to_bits),
            b.cohesive_point().jump_m().map(f64::to_bits)
        );
        assert_eq!(a.dissipated_j_m3().to_bits(), b.dissipated_j_m3().to_bits());
    }
    assert_eq!(
        old.stored_energy_j().to_bits(),
        next.stored_energy_j().to_bits()
    );
    assert_eq!(
        old.fracture_dissipation_j().to_bits(),
        next.fracture_dissipation_j().to_bits()
    );
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:e} != {b:e}");
}

#[test]
fn completed_native_bands_split_with_exact_history_and_continue_as_independent_fragments() {
    with_cx(false, |cx| {
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, budget(), 8, cx).unwrap();
        let sources = materials(&source_card(), &geometry, [0.0, 0.0, 1.0]);
        let old = loaded(&geometry, &sources, FAILED_STRAIN, cx);
        assert!(old.solution().element_states.iter().all(|s| s.failed()));
        let cut = topology(&geometry, [0, 1, 2], 0, 931, cx);
        let target =
            TetLoadGeometry3::new(cut.nodes_m(), cut.tetrahedra(), budget(), 8, cx).unwrap();
        let step = old
            .split_fracture(&cut, &target, release_controls(), cx)
            .unwrap();
        exact_material_transfer(&old, &step.solid);
        assert!(step.solid.solution().increments.is_empty());
        assert_eq!(step.energy.storage_defect_j, 0.0);
        assert_eq!(step.energy.fracture_history_defect_j, 0.0);
        assert_eq!(cut.element_fragments(), &[0, 1]);
        assert_eq!(step.solid.released_faces().len(), 1);
        for (new, &original) in cut.original_nodes().iter().enumerate() {
            assert_eq!(
                step.solid.solution().displacement_m[new],
                old.solution().displacement_m[original]
            );
        }
        let volume: f64 = geometry.element_volumes_m3().iter().sum();
        let expected_storage =
            0.5 * YOUNG * (FAILED_STRAIN[0].powi(2) + FAILED_STRAIN[1].powi(2)) * volume;
        close(step.solid.stored_energy_j(), expected_storage, 1e-15);
        let expected_loss: f64 = old
            .band_geometry()
            .iter()
            .map(|b| GC * b.projected_area_m2())
            .sum();
        close(step.solid.fracture_dissipation_j(), expected_loss, 1e-15);

        // Distinct rigid normal translations were impossible before node
        // splitting. The native volume solver now accepts them independently.
        let mut separated = step.solid.solution().displacement_m.clone();
        for (node, u) in separated.iter_mut().enumerate() {
            u[2] += if cut.tetrahedra()[0].contains(&node) {
                -1.0 / 512.0
            } else {
                1.0 / 512.0
            };
        }
        let continued = step
            .solid
            .advance(
                &prescribed(&separated),
                &[],
                IncrementSettings::default(),
                controls(),
                cx,
            )
            .unwrap();
        exact_material_transfer(&step.solid, &continued.solid);
        assert_eq!(continued.solid.solution().displacement_m, separated);
        assert_eq!(
            continued.solid.released_faces(),
            step.solid.released_faces()
        );
        close(continued.energy.material_work_j, 0.0, 1e-15);
        close(continued.energy.fracture_dissipation_j, 0.0, 1e-15);
        close(continued.energy.support_endpoint_work_j, 0.0, 1e-15);
        for face in continued.solid.released_faces() {
            for i in 0..3 {
                let a = separated[face.paired_nodes[0][i]];
                let b = separated[face.paired_nodes[1][i]];
                close(
                    (0..3).map(|k| face.normal[k] * (b[k] - a[k])).sum(),
                    1.0 / 256.0,
                    1e-15,
                );
            }
        }
    });
}

#[test]
fn partial_compressed_misaligned_or_foreign_geometry_and_cancellation_refuse_atomically() {
    with_cx(false, |cx| {
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, budget(), 8, cx).unwrap();
        let card = source_card();
        let sources = materials(&card, &geometry, [0.0, 0.0, 1.0]);
        let cut = topology(&geometry, [0, 1, 2], 0, 932, cx);
        let target =
            TetLoadGeometry3::new(cut.nodes_m(), cut.tetrahedra(), budget(), 8, cx).unwrap();
        let partial = loaded(&geometry, &sources, [0.0, 0.0, 1.0 / 512.0], cx);
        assert!(
            partial
                .solution()
                .element_states
                .iter()
                .all(|s| s.damage() > 0.0 && !s.failed())
        );
        assert!(
            partial
                .split_fracture(&cut, &target, release_controls(), cx)
                .is_err()
        );
        let failed = loaded(&geometry, &sources, FAILED_STRAIN, cx);
        let compressed = failed
            .advance(
                &prescribed(&affine(&NODES, [0.0, 0.0, -1.0 / 1024.0])),
                &[],
                IncrementSettings::default(),
                controls(),
                cx,
            )
            .unwrap()
            .solid;
        assert!(
            compressed
                .solution()
                .element_states
                .iter()
                .all(|s| s.failed())
        );
        let loose = BrittleTetFractureControls {
            traction_tolerance_pa: 1e100,
            ..release_controls()
        };
        assert!(compressed.split_fracture(&cut, &target, loose, cx).is_err());
        let wrong_sources = materials(&card, &geometry, [1.0, 0.0, 0.0]);
        let wrong_plane = loaded(&geometry, &wrong_sources, [1.0 / 128.0, 0.0, 0.0], cx);
        assert!(
            wrong_plane
                .solution()
                .element_states
                .iter()
                .all(|s| s.failed())
        );
        assert!(
            wrong_plane
                .split_fracture(&cut, &target, loose, cx)
                .is_err()
        );
        assert!(
            failed
                .split_fracture(&cut, &geometry, release_controls(), cx)
                .is_err()
        );
        let shifted_nodes = NODES.map(|p| [p[0] + 3.0, p[1], p[2]]);
        let foreign = TetLoadGeometry3::new(&shifted_nodes, &TETS, budget(), 8, cx).unwrap();
        let foreign_cut = topology(&foreign, [0, 1, 2], 0, 933, cx);
        let foreign_target = TetLoadGeometry3::new(
            foreign_cut.nodes_m(),
            foreign_cut.tetrahedra(),
            budget(),
            8,
            cx,
        )
        .unwrap();
        assert!(
            failed
                .split_fracture(&foreign_cut, &foreign_target, release_controls(), cx)
                .is_err()
        );
        with_cx(true, |cancelled| {
            assert!(matches!(
                failed.split_fracture(&cut, &target, release_controls(), cancelled),
                Err(BrittleTetFractureError::Solid(BrittleTetError::Cancelled))
            ));
        });
        assert!(failed.released_faces().is_empty());
        let first = failed
            .split_fracture(&cut, &target, release_controls(), cx)
            .unwrap();
        let retry = failed
            .split_fracture(&cut, &target, release_controls(), cx)
            .unwrap();
        exact_material_transfer(&failed, &first.solid);
        exact_material_transfer(&first.solid, &retry.solid);
        assert_eq!(
            first.solid.solution().displacement_m,
            retry.solid.solution().displacement_m
        );
        assert_eq!(first.solid.released_faces(), retry.solid.released_faces());
    });
}

#[test]
fn sequential_partial_front_cuts_retain_prior_opening_history_and_closure_guard() {
    with_cx(false, |cx| {
        let nodes = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.5, 0.5, -1.0],
            [0.5, 0.5, 1.0],
        ];
        let tets = [[0, 2, 1, 4], [0, 1, 2, 5], [0, 3, 2, 4], [0, 2, 3, 5]];
        let geometry = TetLoadGeometry3::new(&nodes, &tets, budget(), 16, cx).unwrap();
        let sources = materials(&source_card(), &geometry, [0.0, 0.0, 1.0]);
        let old = loaded(&geometry, &sources, FAILED_STRAIN, cx);
        let first_cut = topology(&geometry, [0, 1, 2], 0, 941, cx);
        let first_geometry = TetLoadGeometry3::new(
            first_cut.nodes_m(),
            first_cut.tetrahedra(),
            budget(),
            16,
            cx,
        )
        .unwrap();
        let first = old
            .split_fracture(&first_cut, &first_geometry, release_controls(), cx)
            .unwrap();
        exact_material_transfer(&old, &first.solid);
        let first_face = first.solid.released_faces()[0];
        let free_vertex = (0..3)
            .find(|&i| first_face.paired_nodes[0][i] != first_face.paired_nodes[1][i])
            .unwrap();
        assert_eq!(
            (0..3)
                .filter(|&i| first_face.paired_nodes[0][i] != first_face.paired_nodes[1][i])
                .count(),
            1
        );
        let mut opened = first.solid.solution().displacement_m.clone();
        opened[first_face.paired_nodes[1][free_vertex]][2] += 1.0 / 4096.0;
        let opened = first
            .solid
            .advance(
                &prescribed(&opened),
                &[],
                IncrementSettings::default(),
                controls(),
                cx,
            )
            .unwrap()
            .solid;
        assert_eq!(
            opened.fracture_dissipation_j(),
            old.fracture_dissipation_j()
        );
        let second_tet = first_geometry.tetrahedra()[2];
        let second_cut = topology(
            &first_geometry,
            [second_tet[0], second_tet[2], second_tet[1]],
            2,
            942,
            cx,
        );
        let second_geometry = TetLoadGeometry3::new(
            second_cut.nodes_m(),
            second_cut.tetrahedra(),
            budget(),
            16,
            cx,
        )
        .unwrap();
        let second = opened
            .split_fracture(&second_cut, &second_geometry, release_controls(), cx)
            .unwrap();
        exact_material_transfer(&opened, &second.solid);
        assert_eq!(second.solid.released_faces().len(), 2);
        let retained = second.solid.released_faces()[0];
        assert_eq!(retained.id, first_face.id);
        assert_eq!(retained.elements, first_face.elements);
        assert_eq!(retained.normal, first_face.normal);
        let current = &second.solid.solution().displacement_m;
        let a = current[retained.paired_nodes[0][free_vertex]];
        let b = current[retained.paired_nodes[1][free_vertex]];
        close(
            (0..3).map(|c| retained.normal[c] * (b[c] - a[c])).sum(),
            1.0 / 4096.0,
            1e-15,
        );
        for vertex in 0..3 {
            for side in 0..2 {
                assert_eq!(
                    second_cut.original_nodes()[retained.paired_nodes[side][vertex]],
                    first_face.paired_nodes[side][vertex]
                );
            }
        }
        let mut closing = current.clone();
        closing[retained.paired_nodes[1][free_vertex]][2] -= 1.0 / 2048.0;
        assert!(matches!(
            second.solid.advance(
                &prescribed(&closing),
                &[],
                IncrementSettings::default(),
                controls(),
                cx
            ),
            Err(BrittleTetError::ReleasedFaceClosure { face: 941, .. })
        ));
        let retry = second
            .solid
            .advance(
                &prescribed(current),
                &[],
                IncrementSettings::default(),
                controls(),
                cx,
            )
            .unwrap();
        exact_material_transfer(&second.solid, &retry.solid);
        assert_eq!(second.solid.released_faces(), retry.solid.released_faces());
    });
}
