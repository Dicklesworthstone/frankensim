//! G1/G3: thermal tensile failure continues as native separated components.
//! Manufactured source data verify reference inventory and work, not material qualification.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::brittle::source::expansion::ResolvedThermalBrittle;
use fs_material::brittle::source::expansion::dynamics::resolve_thermal_brittle_dynamic_material;
use fs_material::brittle::source::{
    BRITTLE_BK_EXPONENT_PROPERTY, BRITTLE_MODE_I_FRACTURE_ENERGY_PROPERTY,
    BRITTLE_MODE_II_FRACTURE_ENERGY_PROPERTY, BRITTLE_SHEAR_STRENGTH_PROPERTY,
    BRITTLE_TENSILE_STRENGTH_PROPERTY,
};
use fs_material::cohesive::source::FRACTURE_ENERGY_PER_AREA_DIMS;
use fs_material::state_point::{
    DENSITY_PROPERTY, INVERSE_TEMPERATURE_DIMS, LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
    MaterialPropertySelection, POISSON_RATIO_PROPERTY, YOUNG_MODULUS_PROPERTY,
};
use fs_qty::{Density, Dims, Pressure, QuantitySpec, Temperature};
use fs_solid::brittle::TetCrackBandGeometry;
use fs_solid::brittle::solid::expansion::{ThermalBrittleTetError, ThermalBrittleTetSolid};
use fs_solid::brittle::solid::fracture::{BrittleTetFractureControls, BrittleTetFractureError};
use fs_solid::brittle::solid::{BrittleTetControls, BrittleTetError};
use fs_solid::brittle::topology::{FractureTopologyBudget, InteriorSeam, TetFractureTopology};
use fs_solid::linear3::TetAssemblyBudget;
use fs_solid::loads3::TetLoadGeometry3;
use fs_solid::plastic3::IncrementSettings;

const NODES: [[f64; 3]; 5] = [
    [0., 0., 0.],
    [1., 0., 0.],
    [0., 1., 0.],
    [0., 0., -1.],
    [0., 0., 1.],
];
const TETS: [[usize; 4]; 2] = [[0, 2, 1, 3], [0, 1, 2, 4]];
const DENSITY: f64 = 6.0;
const GC: f64 = 0.0012;

fn with_cx<R>(run: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        run(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 598,
                kernel_id: 598,
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
        maximum_nodes: 8,
        maximum_tetrahedra: 2,
        maximum_free_dofs: 24,
        minimum_scaled_jacobian: 1e-10,
    }
}

fn settings() -> IncrementSettings {
    IncrementSettings {
        max_newton_iterations: 30,
        relative_tolerance: 1e-11,
        max_backtracks: 18,
        linear_tolerance: 1e-13,
        max_linear_iterations: 120,
    }
}

fn work_controls() -> BrittleTetControls {
    BrittleTetControls {
        material_energy_tolerance_j: 1e-10,
        equilibrium_work_tolerance_j: 1e-10,
        maximum_endpoint_work_defect_j: 1.0,
    }
}

fn fracture_controls() -> BrittleTetFractureControls {
    BrittleTetFractureControls {
        normal_alignment_tolerance: 1e-12,
        traction_tolerance_pa: 1e-10,
        energy_tolerance_j: 1e-10,
        opening_tolerance_m: 1e-12,
    }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}"
    );
}

fn source_card() -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, 1200.0),
        (POISSON_RATIO_PROPERTY, Dims::NONE, 0.0),
        (DENSITY_PROPERTY, Density::DIMS, DENSITY),
        (BRITTLE_TENSILE_STRENGTH_PROPERTY, Pressure::DIMS, 1.2),
        (BRITTLE_SHEAR_STRENGTH_PROPERTY, Pressure::DIMS, 1.5),
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
        (
            LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
            INVERSE_TEMPERATURE_DIMS,
            1e-4,
        ),
    ] {
        let varies = name == DENSITY_PROPERTY;
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: if varies {
                    PropertyValue::Curve {
                        abscissa: "T".into(),
                        abscissa_dims: Temperature::DIMS,
                        knots: vec![(300.0, 0.5 * value), (400.0, value)],
                        dims,
                    }
                } else {
                    PropertyValue::Scalar { value, dims }
                },
                validity: if varies {
                    ValidityDomain::unconstrained().with_quantity(
                        "T",
                        QuantitySpec::dimensional(Temperature::DIMS),
                        300.0,
                        400.0,
                    )
                } else {
                    ValidityDomain::unconstrained()
                },
                interpolation: if varies {
                    InterpolationPolicy::LinearInside
                } else {
                    InterpolationPolicy::ConstantWithinValidity
                },
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured thermal brittle fragment fixture".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured brittle fragments".into(),
            phase: "solid".into(),
            process: "thermal fracture continuation verification".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn reference_sources(
    card: &MaterialCard,
    geometry: &TetLoadGeometry3<'_>,
) -> Vec<ResolvedThermalBrittle> {
    let query = QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), 400.0)
        .unwrap();
    (0..2)
        .map(|element| {
            let band = TetCrackBandGeometry::from_element(geometry, element, [0., 0., 1.]).unwrap();
            resolve_thermal_brittle_dynamic_material(
                card,
                &query,
                &query,
                MaterialPropertySelection::SingleClaimOnly,
                band.material_settings(0.03, 1e-9),
            )
            .unwrap()
        })
        .collect()
}

fn at_temperature(
    sources: &[ResolvedThermalBrittle],
    card: &MaterialCard,
    kelvin: f64,
) -> Vec<ResolvedThermalBrittle> {
    sources
        .iter()
        .map(|s| s.at_temperature(card, kelvin).unwrap())
        .collect()
}

fn prescribed(displacement: &[[f64; 3]]) -> Vec<(usize, f64)> {
    displacement
        .as_flattened()
        .iter()
        .copied()
        .enumerate()
        .collect()
}

fn topology<'g>(geometry: &'g TetLoadGeometry3<'_>, cx: &Cx<'_>) -> TetFractureTopology<'g> {
    TetFractureTopology::new(
        geometry,
        &[InteriorSeam {
            id: 598,
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
fn restrained_thermal_fracture_preserves_native_mechanical_history_and_continues_after_release() {
    with_cx(|cx| {
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, budget(), 8, cx).unwrap();
        let card = source_card();
        let reference = reference_sources(&card, &geometry);
        let cold = at_temperature(&reference, &card, 300.0);
        let warm = at_temperature(&reference, &card, 350.0);
        let initial = ThermalBrittleTetSolid::try_new(&geometry, &reference, budget(), cx).unwrap();
        let cooled = initial
            .advance(
                &cold,
                &prescribed(&[[0.0; 3]; 5]),
                &[],
                settings(),
                work_controls(),
                cx,
            )
            .unwrap()
            .solid;
        for (total, state) in cooled
            .solution()
            .element_strain
            .iter()
            .zip(&cooled.solution().element_states)
        {
            assert_eq!(*total, [0.0; 6]);
            assert!(state.failed());
            for axis in 0..3 {
                close(state.strain()[axis], 0.01, 1e-15);
            }
        }
        // Both V/h=1/2 bands retain their original fracture inventory.
        // Undamaged transverse strains still store (2V)*E*0.01^2 = .04 J.
        close(cooled.stored_energy_j(), 0.04, 1e-13);
        close(cooled.fracture_dissipation_j(), GC, 1e-14);
        let cut = topology(&geometry, cx);
        let target =
            TetLoadGeometry3::new(cut.nodes_m(), cut.tetrahedra(), budget(), 8, cx).unwrap();
        let released = cooled
            .split_fracture(&cut, &target, fracture_controls(), cx)
            .unwrap();
        assert!(core::ptr::eq(
            released.solid.materials().as_ptr(),
            cold.as_ptr()
        ));
        assert_eq!(
            released.solid.solution().element_strain,
            cooled.solution().element_strain
        );
        assert_eq!(
            released.solid.solution().element_states,
            cooled.solution().element_states
        );
        assert_eq!(
            released.solid.solution().element_stress,
            cooled.solution().element_stress
        );
        assert_eq!(released.solid.band_geometry(), cooled.band_geometry());
        assert!(released.solid.solution().increments.is_empty());
        assert_eq!(released.solid.released_faces().len(), 1);
        assert_eq!(released.solid.released_faces()[0].id, 598);
        assert_eq!(released.energy.storage_defect_j, 0.0);
        assert_eq!(released.energy.fracture_history_defect_j, 0.0);
        for (node, &old) in cut.original_nodes().iter().enumerate() {
            assert_eq!(
                released.solid.solution().displacement_m[node],
                cooled.solution().displacement_m[old]
            );
        }

        // A later source change remains thermal: total strain follows free
        // expansion while the retained mechanical history stays unchanged.
        let delta_free = warm[0].free_linear_strain() - cold[0].free_linear_strain();
        let expansion: Vec<_> = target
            .nodes_m()
            .iter()
            .map(|x| x.map(|q| delta_free * q))
            .collect();
        let continued = released
            .solid
            .advance(
                &warm,
                &prescribed(&expansion),
                &[],
                settings(),
                work_controls(),
                cx,
            )
            .unwrap();
        assert!(core::ptr::eq(
            continued.solid.materials().as_ptr(),
            warm.as_ptr()
        ));
        assert_eq!(
            continued.solid.released_faces(),
            released.solid.released_faces()
        );
        for (total, state) in continued
            .solid
            .solution()
            .element_strain
            .iter()
            .zip(&continued.solid.solution().element_states)
        {
            for axis in 0..3 {
                close(total[axis], delta_free, 1e-15);
                close(state.strain()[axis], 0.01, 1e-15);
            }
            assert!(state.failed());
        }
        close(
            continued.solid.stored_energy_j(),
            cooled.stored_energy_j(),
            1e-13,
        );
        close(continued.solid.fracture_dissipation_j(), GC, 1e-14);
        close(
            continued.energy.mechanical_work_j + continued.energy.eigenstrain_work_j,
            0.0,
            1e-13,
        );
        close(continued.energy.total_balance_error_j, 0.0, 1e-13);
    });
}

#[test]
fn unfinished_loaded_and_closing_thermal_faces_refuse_without_losing_original_history() {
    with_cx(|cx| {
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, budget(), 8, cx).unwrap();
        let card = source_card();
        let reference = reference_sources(&card, &geometry);
        let cold = at_temperature(&reference, &card, 300.0);
        let warm = at_temperature(&reference, &card, 350.0);
        let initial = ThermalBrittleTetSolid::try_new(&geometry, &reference, budget(), cx).unwrap();
        let partial = initial
            .advance(
                &warm,
                &prescribed(&[[0.0; 3]; 5]),
                &[],
                settings(),
                work_controls(),
                cx,
            )
            .unwrap()
            .solid;
        assert!(
            partial
                .solution()
                .element_states
                .iter()
                .all(|state| !state.failed())
        );
        let cut = topology(&geometry, cx);
        let target =
            TetLoadGeometry3::new(cut.nodes_m(), cut.tetrahedra(), budget(), 8, cx).unwrap();
        assert!(matches!(
            partial.split_fracture(&cut, &target, fracture_controls(), cx),
            Err(BrittleTetFractureError::Solid(
                BrittleTetError::Material { .. }
            ))
        ));
        assert!(partial.released_faces().is_empty());

        let cooled = initial
            .advance(
                &cold,
                &prescribed(&[[0.0; 3]; 5]),
                &[],
                settings(),
                work_controls(),
                cx,
            )
            .unwrap()
            .solid;
        let compression = NODES.map(|x| [0.0, 0.0, -0.011 * x[2]]);
        let loaded = cooled
            .advance(
                &cold,
                &prescribed(&compression),
                &[],
                settings(),
                work_controls(),
                cx,
            )
            .unwrap()
            .solid;
        assert!(
            loaded
                .solution()
                .element_states
                .iter()
                .all(|state| state.failed())
        );
        assert!(
            loaded
                .solution()
                .element_stress
                .iter()
                .all(|stress| stress[2] < 0.0)
        );
        assert!(matches!(
            loaded.split_fracture(&cut, &target, fracture_controls(), cx),
            Err(BrittleTetFractureError::Solid(
                BrittleTetError::Material { .. }
            ))
        ));
        close(loaded.fracture_dissipation_j(), GC, 1e-14);

        let split = cooled
            .split_fracture(&cut, &target, fracture_controls(), cx)
            .unwrap()
            .solid;
        let snapshot = split.solution().clone();
        let delta_free = warm[0].free_linear_strain() - cold[0].free_linear_strain();
        let mut closing: Vec<_> = target
            .nodes_m()
            .iter()
            .map(|x| x.map(|q| delta_free * q))
            .collect();
        for &node in &target.tetrahedra()[1] {
            closing[node][2] -= 1e-4;
        }
        assert!(matches!(
            split.advance(
                &warm,
                &prescribed(&closing),
                &[],
                settings(),
                work_controls(),
                cx,
            ),
            Err(ThermalBrittleTetError::Solid(
                BrittleTetError::ReleasedFaceClosure { face: 598, .. }
            ))
        ));
        assert_eq!(split.solution().displacement_m, snapshot.displacement_m);
        assert_eq!(split.solution().element_states, snapshot.element_states);
        assert_eq!(split.solution().element_strain, snapshot.element_strain);
        assert!(core::ptr::eq(split.materials().as_ptr(), cold.as_ptr()));
        close(split.fracture_dissipation_j(), GC, 1e-14);
    });
}
