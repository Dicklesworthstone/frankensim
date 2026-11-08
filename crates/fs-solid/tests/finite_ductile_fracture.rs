//! G0/G3: explicit finite-support failure, native topology transfer and motion.
//! Synthetic constants exercise mechanics; they do not qualify a real alloy.
#![deny(unsafe_code)]

#[path = "finite_ductile_fracture/contact.rs"]
mod contact;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::ductile::evolution::{DuctileJ2Settings, DuctileWorkEvolution};
use fs_material::ductile::{DuctileDomain, TabulatedDuctileLocus};
use fs_material::plastic::energy::EnergyTolerances;
use fs_material::plastic::finite::ductile::FiniteDuctileJ2;
use fs_material::plastic::finite::{FiniteJ2, FiniteJ2Limits, FiniteJ2Parameters};
use fs_solid::brittle::TetCrackBandGeometry;
use fs_solid::brittle::topology::{FractureTopologyBudget, InteriorSeam, TetFractureTopology};
use fs_solid::ductile::{DuctileStepConditions, DuctileTemperatures};
use fs_solid::hyper3::plastic::FinitePlasticStepControls;
use fs_solid::hyper3::plastic::ductile::dynamics::{
    FiniteDuctileMotionControls, FiniteDuctileTetMotion,
};
use fs_solid::hyper3::plastic::ductile::fracture::{
    FiniteDuctileFractureControls, FiniteDuctileFractureError,
};
use fs_solid::hyper3::plastic::ductile::{FiniteDuctileTetError, FiniteDuctileTetSolid};
use fs_solid::loads3::TetLoadGeometry3;
use fs_solid::{IncrementSettings, TetAssemblyBudget};

const NODES: [[f64; 3]; 5] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, -1.0],
    [0.0, 0.0, 1.0],
];
const TETS: [[usize; 4]; 2] = [[0, 2, 1, 3], [0, 1, 2, 4]];
type Model = FiniteDuctileJ2<TabulatedDuctileLocus>;

fn with_cx<R>(f: impl FnOnce(&CancelGate, &Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 211,
                kernel_id: 908,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&gate, &cx)
    })
}
fn bands(geometry: &TetLoadGeometry3<'_>) -> Vec<TetCrackBandGeometry> {
    (0..geometry.tetrahedra().len())
        .map(|element| {
            TetCrackBandGeometry::from_element(geometry, element, [0.0, 0.0, 1.0]).unwrap()
        })
        .collect()
}
fn model(width: f64, terminal: bool) -> Model {
    let native = FiniteJ2::new(
        FiniteJ2Parameters {
            shear_modulus_pa: 500.0,
            bulk_modulus_pa: 1000.0,
            yield_stress_pa: 3.0,
            hardening_modulus_pa: 0.0,
        },
        FiniteJ2Limits {
            minimum_jacobian: 0.25,
            maximum_elastic_log_strain: 0.5,
            maximum_plastic_log_strain: 2.0,
            maximum_accumulated_plastic_strain: 3.0,
        },
    )
    .unwrap();
    let locus = TabulatedDuctileLocus::new(
        vec![(-20.0, 0.002), (20.0, 0.002)],
        DuctileDomain {
            triaxiality: [-20.0, 20.0],
            strain_rate_s_inv: [1e-10, 10.0],
            temperature_k: [300.0, 300.0],
        },
    )
    .unwrap();
    let model = FiniteDuctileJ2::new(
        native,
        locus,
        DuctileJ2Settings {
            characteristic_length_m: width,
            post_initiation_plastic_work_j_m2: 0.005,
            maximum_damage_increment: 1.0,
            energy_tolerances: EnergyTolerances {
                energy_j_m3: 1e-9,
                yield_pa: 1e-8,
            },
        },
    )
    .unwrap();
    if terminal {
        model.with_evolution(DuctileWorkEvolution::FiniteSupportQuadratic)
    } else {
        model
    }
}
fn conditions() -> DuctileStepConditions<'static> {
    DuctileStepConditions {
        duration_s: 0.1,
        temperatures_k: DuctileTemperatures::Uniform(300.0),
    }
}
fn controls() -> FinitePlasticStepControls {
    FinitePlasticStepControls {
        equilibrium_work_tolerance_j: 1e-8,
        maximum_endpoint_work_defect_j: 1.0,
        yield_tolerance_pa: 1e-8,
    }
}
fn motion_controls() -> FiniteDuctileMotionControls {
    FiniteDuctileMotionControls {
        material: controls(),
        maximum_numerical_kinetic_loss_j: 1.0,
    }
}
fn release_controls() -> FiniteDuctileFractureControls {
    FiniteDuctileFractureControls {
        traction_tolerance_pa: 1e-10,
        energy_tolerance_j: 1e-10,
        opening_tolerance_m: 1e-12,
        normal_alignment_tolerance: 1e-10,
    }
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() <= 3e-9 * (1.0 + b.abs()), "{a:e} != {b:e}");
}
fn failed_f() -> [f64; 9] {
    // Positive dilation avoids even a sub-ulp negative native log(J).
    // Ry(pi/2) makes the original +Z facet normal become current +X.
    [
        0.0,
        0.0,
        (-0.015_f64).exp(),
        0.0,
        (-0.015_f64).exp(),
        0.0,
        -(0.06_f64).exp(),
        0.0,
        0.0,
    ]
}
fn affine(nodes: &[[f64; 3]], f: [f64; 9]) -> Vec<(usize, f64)> {
    nodes
        .iter()
        .enumerate()
        .flat_map(|(node, x)| {
            (0..3).map(move |i| {
                (
                    3 * node + i,
                    (0..3)
                        .map(|j| (f[3 * i + j] - if i == j { 1.0 } else { 0.0 }) * x[j])
                        .sum(),
                )
            })
        })
        .collect()
}
fn failed_solid<'g, 'm>(
    geometry: &'g TetLoadGeometry3<'g>,
    materials: &'m [&'m Model],
    bands: &'g [TetCrackBandGeometry],
    cx: &Cx<'_>,
) -> FiniteDuctileTetSolid<'g, 'm, TabulatedDuctileLocus> {
    let mut solid = FiniteDuctileTetSolid::try_new(
        geometry,
        materials,
        bands,
        TetAssemblyBudget::standard(),
        cx,
    )
    .unwrap();
    // Rotate the stress-free virgin body first. Combining a large rotation
    // with newly created bulk tension would manufacture a large endpoint
    // work-quadrature defect unrelated to the failure transfer under test.
    let rotation = [0.0, 0.0, 1.0, 0.0, 1.0, 0.0, -1.0, 0.0, 0.0];
    solid = solid
        .advance(
            &affine(geometry.nodes_m(), rotation),
            &[],
            [0.0; 3],
            conditions(),
            IncrementSettings::default(),
            controls(),
            cx,
        )
        .unwrap()
        .solid;
    // Resolve real positive-integrity plastic work and stored-energy release
    // before terminal failure. A single virgin-to-zero-integrity jump cannot
    // provide this history: endpoint plastic work and old active storage vanish.
    for step in 1..=32 {
        let strain = 0.05 * step as f64 / 32.0;
        let longitudinal = (0.01 + strain).exp();
        let transverse = (0.01 - 0.5 * strain).exp();
        let f = [
            0.0,
            0.0,
            transverse,
            0.0,
            transverse,
            0.0,
            -longitudinal,
            0.0,
            0.0,
        ];
        solid = solid
            .advance(
                &affine(geometry.nodes_m(), f),
                &[],
                [0.0; 3],
                conditions(),
                IncrementSettings::default(),
                controls(),
                cx,
            )
            .unwrap()
            .solid;
    }
    solid
}
fn topology<'g>(geometry: &'g TetLoadGeometry3<'_>, cx: &Cx<'_>) -> TetFractureTopology<'g> {
    TetFractureTopology::new(
        geometry,
        &[InteriorSeam {
            id: 71,
            nodes: [0, 1, 2],
            side_a_element: 0,
        }],
        FractureTopologyBudget {
            maximum_output_nodes: 8,
            maximum_face_entries: 8,
            maximum_seams: 1,
        },
        cx,
    )
    .unwrap()
}

#[test]
fn terminal_split_retains_every_native_history_and_both_consistent_density_modes() {
    with_cx(|_, cx| {
        let geometry =
            TetLoadGeometry3::new(&NODES, &TETS, TetAssemblyBudget::standard(), 8, cx).unwrap();
        let original_bands = bands(&geometry);
        let material = model(original_bands[0].band_width_m(), true);
        let materials = [&material, &material];
        let topology = topology(&geometry, cx);
        let target = TetLoadGeometry3::new(
            topology.nodes_m(),
            topology.tetrahedra(),
            TetAssemblyBudget::standard(),
            8,
            cx,
        )
        .unwrap();
        let target_bands = bands(&target);
        let density = [2.0, 3.0, 5.0, 7.0, 11.0];
        let velocity: Vec<_> = NODES
            .iter()
            .map(|x| [0.2 + 0.1 * x[0], -0.1 + 0.02 * x[1], 0.05 + 0.03 * x[2]])
            .collect();
        for graded in [false, true] {
            let solid = failed_solid(&geometry, &materials, &original_bands, cx);
            assert!(solid.element_states().iter().all(|s| s.is_fully_failed()));
            assert!(
                solid
                    .element_states()
                    .iter()
                    .all(|s| s.compression_stored_energy_j_m3() == 0.0)
            );
            let histories = solid.element_states().to_vec();
            let deformation = solid.element_deformation().to_vec();
            let stress = solid.element_first_piola_pa().to_vec();
            let native = if graded {
                FiniteDuctileTetMotion::try_new_with_nodal_density(solid, &density, &velocity, cx)
            } else {
                FiniteDuctileTetMotion::try_new(solid, &[6.0, 12.0], &velocity, cx)
            }
            .unwrap();
            let cut = native
                .split_fracture(&topology, &target, &target_bands, release_controls(), cx)
                .unwrap();
            let split = cut.motion.motion();
            assert_eq!(split.solid().element_states(), histories);
            assert_eq!(split.solid().element_deformation(), deformation);
            assert_eq!(split.solid().element_first_piola_pa(), stress);
            assert_eq!(native.solid().element_states(), histories);
            assert!(std::ptr::eq(split.solid().materials()[0], &material));
            assert_eq!(cut.inventory_before, cut.inventory_after);
            assert_eq!(cut.fragments.len(), 2);
            near(
                cut.inventory_after.mass_kg,
                if graded { 19.0 / 12.0 } else { 3.0 },
            );
            assert_eq!(native.reference_mass_kg(), split.reference_mass_kg());
            assert_eq!(cut.energy.stored_before_j, cut.energy.stored_after_j);
            assert_eq!(
                cut.energy.plastic_dissipation_before_j,
                cut.energy.plastic_dissipation_after_j
            );
            assert_eq!(
                cut.energy.damage_dissipation_before_j,
                cut.energy.damage_dissipation_after_j
            );
            assert_eq!(
                cut.energy.endpoint_defect_before_j,
                cut.energy.endpoint_defect_after_j
            );
            near(cut.energy.kinetic_defect_j, 0.0);
            assert!(cut.energy.plastic_dissipation_after_j > 0.0);
            assert!(cut.energy.damage_dissipation_after_j > 0.0);
            let face = cut.motion.released_faces()[0];
            assert_eq!(face.id, 71);
            assert_eq!(face.elements, [0, 1]);
            near(face.normal[0], 1.0);
            near(face.normal[1], 0.0);
            near(face.normal[2], 0.0);
            for (next, &previous) in topology.original_nodes().iter().enumerate() {
                assert_eq!(
                    split.solid().displacement_m()[next],
                    native.solid().displacement_m()[previous]
                );
                assert_eq!(split.velocity_mps()[next], velocity[previous]);
                if graded {
                    assert_eq!(
                        split.reference_nodal_density_kg_m3().unwrap()[next],
                        density[previous]
                    );
                }
                for i in 0..3 {
                    assert_eq!(
                        target.nodes_m()[next][i] + split.solid().displacement_m()[next][i],
                        NODES[previous][i] + native.solid().displacement_m()[previous][i]
                    );
                }
            }
            // A nonuniform acceleration probes the full original consistent
            // operator, including its off-diagonal graded-density entries.
            let acceleration: Vec<_> = (0..5)
                .map(|i| {
                    [
                        0.13 * (i as f64 + 1.0),
                        -0.17 * (i as f64),
                        0.07 - 0.03 * i as f64,
                    ]
                })
                .collect();
            let copied: Vec<_> = topology
                .original_nodes()
                .iter()
                .map(|&i| acceleration[i])
                .collect();
            let original = native
                .reference_acceleration_load_n(&acceleration, cx)
                .unwrap();
            let next = split.reference_acceleration_load_n(&copied, cx).unwrap();
            let mut returned = vec![0.0; 15];
            for (dof, force) in next {
                returned[3 * topology.original_nodes()[dof / 3] + dof % 3] += force;
            }
            for (dof, force) in original {
                near(returned[dof], force);
            }
        }
    });
}

#[test]
fn separated_native_fragments_accelerate_and_refuse_uncontrolled_closure() {
    with_cx(|_, cx| {
        let geometry =
            TetLoadGeometry3::new(&NODES, &TETS, TetAssemblyBudget::standard(), 8, cx).unwrap();
        let original_bands = bands(&geometry);
        let material = model(original_bands[0].band_width_m(), true);
        let materials = [&material, &material];
        let solid = failed_solid(&geometry, &materials, &original_bands, cx);
        let native =
            FiniteDuctileTetMotion::try_new(solid, &[6.0, 12.0], &[[0.0; 3]; 5], cx).unwrap();
        let topology = topology(&geometry, cx);
        let target = TetLoadGeometry3::new(
            topology.nodes_m(),
            topology.tetrahedra(),
            TetAssemblyBudget::standard(),
            8,
            cx,
        )
        .unwrap();
        let target_bands = bands(&target);
        let cut = native
            .split_fracture(&topology, &target, &target_bands, release_controls(), cx)
            .unwrap();
        let released = cut.motion;
        let mut acceleration = vec![[0.0; 3]; target.nodes_m().len()];
        for (element, tet) in target.tetrahedra().iter().enumerate() {
            for &node in tet {
                acceleration[node][0] = if element == 0 { -0.1 } else { 0.1 };
            }
        }
        let loads = released
            .motion()
            .reference_acceleration_load_n(&acceleration, cx)
            .unwrap();
        let moved = released
            .advance(
                &[],
                &loads,
                [0.0; 3],
                conditions(),
                IncrementSettings::default(),
                motion_controls(),
                cx,
            )
            .unwrap();
        for (node, a) in acceleration.iter().enumerate() {
            near(moved.motion.motion().velocity_mps()[node][0], 0.1 * a[0]);
            near(
                moved.motion.motion().solid().displacement_m()[node][0]
                    - released.motion().solid().displacement_m()[node][0],
                0.01 * a[0],
            );
        }
        assert_eq!(moved.motion.released_faces(), released.released_faces());
        near(moved.material_energy.plastic_dissipation_j, 0.0);
        near(moved.material_energy.damage_release_j, 0.0);
        near(moved.balance_error_j, 0.0);
        let closing: Vec<_> = loads.iter().map(|&(dof, force)| (dof, -force)).collect();
        let history = released.motion().solid().element_states().to_vec();
        assert!(matches!(
            released.advance(
                &[],
                &closing,
                [0.0; 3],
                conditions(),
                IncrementSettings::default(),
                motion_controls(),
                cx
            ),
            Err(FiniteDuctileFractureError::Face { .. })
        ));
        assert_eq!(released.motion().solid().element_states(), history);
        assert_eq!(released.motion().velocity_mps(), &[[0.0; 3]; 8]);
    });
}

#[test]
fn incomplete_compressed_wrong_band_and_cancelled_transfers_are_transactional() {
    with_cx(|gate, cx| {
        let geometry =
            TetLoadGeometry3::new(&NODES, &TETS, TetAssemblyBudget::standard(), 8, cx).unwrap();
        let original_bands = bands(&geometry);
        let material = model(original_bands[0].band_width_m(), true);
        let materials = [&material, &material];
        let topology = topology(&geometry, cx);
        let target = TetLoadGeometry3::new(
            topology.nodes_m(),
            topology.tetrahedra(),
            TetAssemblyBudget::standard(),
            8,
            cx,
        )
        .unwrap();
        let target_bands = bands(&target);
        let virgin = FiniteDuctileTetSolid::try_new(
            &geometry,
            &materials,
            &original_bands,
            TetAssemblyBudget::standard(),
            cx,
        )
        .unwrap();
        let virgin =
            FiniteDuctileTetMotion::try_new(virgin, &[6.0, 12.0], &[[0.0; 3]; 5], cx).unwrap();
        assert!(matches!(
            virgin.split_fracture(&topology, &target, &target_bands, release_controls(), cx),
            Err(FiniteDuctileFractureError::Face { .. })
        ));
        assert!(
            virgin
                .solid()
                .element_states()
                .iter()
                .all(|s| s.native_state().alpha == 0.0)
        );
        let exponential = model(original_bands[0].band_width_m(), false);
        let exponential_materials = [&exponential, &exponential];
        let attenuated = failed_solid(&geometry, &exponential_materials, &original_bands, cx);
        assert!(
            attenuated
                .element_states()
                .iter()
                .all(|s| s.damage() > 0.99 && !s.is_fully_failed())
        );
        let attenuated =
            FiniteDuctileTetMotion::try_new(attenuated, &[6.0, 12.0], &[[0.0; 3]; 5], cx).unwrap();
        assert!(matches!(
            attenuated.split_fracture(&topology, &target, &target_bands, release_controls(), cx),
            Err(FiniteDuctileFractureError::Face { .. })
        ));
        let failed = failed_solid(&geometry, &materials, &original_bands, cx);
        let compressed_f = failed_f().map(|v| 0.98 * v);
        let compressed = failed
            .advance(
                &affine(&NODES, compressed_f),
                &[],
                [0.0; 3],
                conditions(),
                IncrementSettings::default(),
                controls(),
                cx,
            )
            .unwrap()
            .solid;
        assert!(
            compressed
                .element_states()
                .iter()
                .all(|s| s.is_fully_failed() && s.compression_stored_energy_j_m3() > 0.0)
        );
        let compressed =
            FiniteDuctileTetMotion::try_new(compressed, &[6.0, 12.0], &[[0.0; 3]; 5], cx).unwrap();
        let preserved = compressed.solid().element_states().to_vec();
        assert!(matches!(
            compressed.split_fracture(&topology, &target, &target_bands, release_controls(), cx),
            Err(FiniteDuctileFractureError::Face { .. })
        ));
        assert_eq!(compressed.solid().element_states(), preserved);
        let failed =
            FiniteDuctileTetMotion::try_new(failed, &[6.0, 12.0], &[[0.0; 3]; 5], cx).unwrap();
        let wrong: Vec<_> = (0..2)
            .map(|element| {
                TetCrackBandGeometry::from_element(&target, element, [1.0, 0.0, 0.0]).unwrap()
            })
            .collect();
        let preserved = failed.solid().element_states().to_vec();
        assert!(matches!(
            failed.split_fracture(&topology, &target, &wrong, release_controls(), cx),
            Err(FiniteDuctileFractureError::Invalid(_))
        ));
        assert_eq!(failed.solid().element_states(), preserved);
        gate.request();
        assert!(matches!(
            failed.split_fracture(&topology, &target, &target_bands, release_controls(), cx),
            Err(FiniteDuctileFractureError::Native(
                FiniteDuctileTetError::Solve(fs_solid::plastic3::SmallStrainTetError::Cancelled)
            ))
        ));
        assert_eq!(failed.solid().element_states(), preserved);
    });
}
