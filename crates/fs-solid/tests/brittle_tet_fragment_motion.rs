//! G1/G3: actual failed bulk bodies continue through native fragment dynamics.
//! Manufactured material data verify mechanics, not experimental qualification.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::brittle::source::dynamics::resolve_brittle_dynamic_material;
use fs_material::brittle::source::{
    BRITTLE_BK_EXPONENT_PROPERTY, BRITTLE_MODE_I_FRACTURE_ENERGY_PROPERTY,
    BRITTLE_MODE_II_FRACTURE_ENERGY_PROPERTY, BRITTLE_SHEAR_STRENGTH_PROPERTY,
    BRITTLE_TENSILE_STRENGTH_PROPERTY, ResolvedBrittleCrackBand,
};
use fs_material::cohesive::source::FRACTURE_ENERGY_PER_AREA_DIMS;
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, POISSON_RATIO_PROPERTY, YOUNG_MODULUS_PROPERTY,
};
use fs_qty::{Density, Dims, Pressure, QuantitySpec, Temperature};
use fs_solid::brittle::TetCrackBandGeometry;
use fs_solid::brittle::solid::dynamics::fracture::BrittleTetFractureMotionError;
use fs_solid::brittle::solid::dynamics::midpoint::BrittleTetMidpointControls;
use fs_solid::brittle::solid::dynamics::{
    BrittleTetDynamicControls, BrittleTetDynamicError, BrittleTetMotion,
};
use fs_solid::brittle::solid::fracture::BrittleTetFractureControls;
use fs_solid::brittle::solid::{BrittleTetControls, BrittleTetError, BrittleTetSolid};
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
                seed: 86,
                kernel_id: 86,
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

fn midpoint_controls() -> BrittleTetMidpointControls {
    BrittleTetMidpointControls {
        material_energy_tolerance_j: 1e-10,
        equilibrium_work_tolerance_j: 1e-10,
        kinetic_work_tolerance_j: 1e-10,
        total_energy_tolerance_j: 1e-10,
        difference_step_strain: Some(1e-7),
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
    for (name, dims, value, varies) in [
        (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, 1200.0, true),
        (POISSON_RATIO_PROPERTY, Dims::NONE, 0.0, false),
        (DENSITY_PROPERTY, Density::DIMS, DENSITY, true),
        (
            BRITTLE_TENSILE_STRENGTH_PROPERTY,
            Pressure::DIMS,
            1.2,
            false,
        ),
        (BRITTLE_SHEAR_STRENGTH_PROPERTY, Pressure::DIMS, 1.5, false),
        (
            BRITTLE_MODE_I_FRACTURE_ENERGY_PROPERTY,
            FRACTURE_ENERGY_PER_AREA_DIMS,
            GC,
            false,
        ),
        (
            BRITTLE_MODE_II_FRACTURE_ENERGY_PROPERTY,
            FRACTURE_ENERGY_PER_AREA_DIMS,
            2.0 * GC,
            false,
        ),
        (BRITTLE_BK_EXPONENT_PROPERTY, Dims::NONE, 1.5, false),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: if varies {
                    PropertyValue::Curve {
                        abscissa: "T".into(),
                        abscissa_dims: Temperature::DIMS,
                        knots: vec![(300.0, value), (600.0, 0.5 * value)],
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
                        600.0,
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
                    source: "manufactured brittle fragment motion fixture".into(),
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
            process: "native numerical verification".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn sources(card: &MaterialCard, geometry: &TetLoadGeometry3<'_>) -> Vec<ResolvedBrittleCrackBand> {
    (0..2)
        .map(|element| {
            let band = TetCrackBandGeometry::from_element(geometry, element, [0., 0., 1.]).unwrap();
            resolve_brittle_dynamic_material(
                card,
                &QueryPoint::new()
                    .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), 300.0)
                    .unwrap(),
                MaterialPropertySelection::SingleClaimOnly,
                band.material_settings(0.03, 1e-9),
            )
            .unwrap()
        })
        .collect()
}

fn failed<'a>(
    geometry: &'a TetLoadGeometry3<'a>,
    materials: &'a [ResolvedBrittleCrackBand],
    cx: &Cx<'_>,
) -> BrittleTetSolid<'a> {
    let prescribed: Vec<_> = NODES
        .iter()
        .enumerate()
        .flat_map(|(node, x)| {
            (0..3).map(move |c| (3 * node + c, if c == 2 { 0.009 * x[2] } else { 0.0 }))
        })
        .collect();
    let solid = BrittleTetSolid::try_new(geometry, materials, budget(), cx)
        .unwrap()
        .advance(
            &prescribed,
            &[],
            settings(),
            BrittleTetControls {
                material_energy_tolerance_j: 1e-10,
                equilibrium_work_tolerance_j: 1e-10,
                maximum_endpoint_work_defect_j: 1.0,
            },
            cx,
        )
        .unwrap()
        .solid;
    assert!(solid.solution().element_states.iter().all(|s| s.failed()));
    // Each original band has V/h = 1/2 m2. This historical fracture energy
    // remains charged to those two material volumes when the seam is cut.
    close(solid.fracture_dissipation_j(), GC, 1e-15);
    solid
}

fn topology<'g>(geometry: &'g TetLoadGeometry3<'_>, cx: &Cx<'_>) -> TetFractureTopology<'g> {
    TetFractureTopology::new(
        geometry,
        &[InteriorSeam {
            id: 86,
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

fn independent_kinetic(geometry: &TetLoadGeometry3<'_>, velocities: &[[f64; 3]]) -> f64 {
    // Degree-two tetrahedral Gauss rule, independent of the mass assembly.
    let a = 0.5854101966249685;
    let b = 0.1381966011250105;
    let mut energy = 0.0;
    for (tet, &volume) in geometry
        .tetrahedra()
        .iter()
        .zip(geometry.element_volumes_m3())
    {
        for q in 0..4 {
            let v: [f64; 3] = core::array::from_fn(|c| {
                (0..4)
                    .map(|i| (if i == q { a } else { b }) * velocities[tet[i]][c])
                    .sum()
            });
            energy += DENSITY * volume / 8.0 * v.iter().map(|x| x * x).sum::<f64>();
        }
    }
    energy
}

fn independent_momentum(geometry: &TetLoadGeometry3<'_>, velocities: &[[f64; 3]]) -> [f64; 3] {
    let mut momentum = [0.0; 3];
    for (tet, &volume) in geometry
        .tetrahedra()
        .iter()
        .zip(geometry.element_volumes_m3())
    {
        for &node in tet {
            for c in 0..3 {
                momentum[c] += DENSITY * volume / 4.0 * velocities[node][c];
            }
        }
    }
    momentum
}

#[test]
fn failed_original_body_splits_and_continues_with_native_mass_and_midpoint_momentum() {
    with_cx(false, |cx| {
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, budget(), 8, cx).unwrap();
        let card = source_card();
        let materials = sources(&card, &geometry);
        let solid = failed(&geometry, &materials, cx);
        let velocity = NODES.map(|x| {
            [
                0.003 + 0.004 * x[0],
                -0.002 + 0.001 * x[1],
                0.001 + 0.0005 * x[2],
            ]
        });
        let initial = BrittleTetMotion::try_from_solid(&solid, &velocity, 1.25, 5, cx).unwrap();
        let old = initial
            .advance_midpoint(0.01, &[], &[], settings(), midpoint_controls(), cx)
            .unwrap()
            .motion;
        assert!(old.solid().stored_energy_j() > 0.0);
        let cuts = topology(&geometry, cx);
        let target =
            TetLoadGeometry3::new(cuts.nodes_m(), cuts.tetrahedra(), budget(), 8, cx).unwrap();
        let handoff = old
            .split_fracture(&cuts, &target, fracture_controls(), 1e-12, cx)
            .unwrap();
        assert_eq!(handoff.inventory_before, handoff.inventory_after);
        assert_eq!(handoff.fragments.len(), 2);
        assert_eq!(
            handoff.energy.stored_before_j,
            handoff.energy.stored_after_j
        );
        assert_eq!(
            handoff.energy.fracture_dissipation_before_j,
            handoff.energy.fracture_dissipation_after_j
        );
        close(
            handoff.kinetic_after_j,
            independent_kinetic(&target, handoff.motion.velocity_mps()),
            1e-18,
        );
        close(handoff.motion.reference_mass_kg(), 2.0, 1e-15);
        assert_eq!(handoff.motion.accepted_steps(), 1);
        assert_eq!(handoff.motion.remaining_steps(), 4);
        assert_eq!(handoff.motion.time_s(), old.time_s());
        assert!(core::ptr::eq(
            handoff.motion.solid().materials().as_ptr(),
            materials.as_ptr()
        ));
        for (node, &original) in cuts.original_nodes().iter().enumerate() {
            assert_eq!(
                handoff.motion.velocity_mps()[node],
                old.velocity_mps()[original]
            );
            assert_eq!(
                handoff.motion.solid().solution().displacement_m[node],
                old.solid().solution().displacement_m[original]
            );
        }
        let momentum = independent_momentum(&target, handoff.motion.velocity_mps());
        let initial_total = handoff.motion.kinetic_energy_j()
            + handoff.motion.solid().stored_energy_j()
            + handoff.motion.solid().fracture_dissipation_j();
        let mut force = vec![0.0; 3 * target.nodes_m().len()];
        let acceleration = 0.1;
        for (element, tet) in target.tetrahedra().iter().enumerate() {
            for &node in tet {
                force[3 * node + 2] += (if element == 0 { -1.0 } else { 1.0 })
                    * DENSITY
                    * target.element_volumes_m3()[element]
                    * acceleration
                    / 4.0;
            }
        }
        let forces: Vec<_> = force
            .into_iter()
            .enumerate()
            .filter(|(_, f)| *f != 0.0)
            .collect();
        let mut motion = handoff.motion;
        let mut work = 0.0;
        for step_index in 1..=4 {
            let step = motion
                .advance_midpoint(0.02, &[], &forces, settings(), midpoint_controls(), cx)
                .unwrap();
            work += step.energy.external_mean_work_j();
            assert_eq!(step.energy.inertia.numerical_loss_j, 0.0);
            close(
                step.motion.kinetic_energy_j(),
                independent_kinetic(&target, step.motion.velocity_mps()),
                1e-17,
            );
            close(
                step.motion.kinetic_energy_j()
                    + step.motion.solid().stored_energy_j()
                    + step.motion.solid().fracture_dissipation_j(),
                initial_total + work,
                1e-11,
            );
            for (actual, expected) in independent_momentum(&target, step.motion.velocity_mps())
                .into_iter()
                .zip(momentum)
            {
                close(actual, expected, 1e-11);
            }
            let face = &cuts.seams()[0];
            for vertex in 0..3 {
                let u = &step.motion.solid().solution().displacement_m;
                close(
                    u[face.side_b[vertex]][2] - u[face.side_a[vertex]][2],
                    acceleration * (step_index as f64 * 0.02).powi(2),
                    1e-11,
                );
            }
            close(step.motion.solid().fracture_dissipation_j(), GC, 1e-14);
            assert!(
                step.motion
                    .solid()
                    .solution()
                    .element_states
                    .iter()
                    .all(|s| s.failed())
            );
            motion = step.motion;
        }
        assert_eq!(motion.accepted_steps(), 5);
        assert_eq!(motion.remaining_steps(), 0);
        assert!(matches!(
            motion.advance_midpoint(0.02, &[], &forces, settings(), midpoint_controls(), cx),
            Err(BrittleTetDynamicError::StepBudget {
                accepted_steps: 5,
                maximum_steps: 5
            })
        ));
    });
}

#[test]
fn split_motion_retains_thermal_history_frozen_density_and_geometric_opening_guard() {
    with_cx(false, |cx| {
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, budget(), 8, cx).unwrap();
        let card = source_card();
        let cold = sources(&card, &geometry);
        let hot: Vec<_> = cold
            .iter()
            .map(|source| source.at_temperature(&card, "T", 600.0).unwrap())
            .collect();
        let solid = failed(&geometry, &cold, cx);
        let old =
            BrittleTetMotion::try_from_solid(&solid, &[[0.003, -0.002, 0.001]; 5], 2.0, 4, cx)
                .unwrap();
        let heated = old
            .advance_midpoint_at_temperature(
                0.01,
                &hot,
                "T",
                &[],
                &[],
                settings(),
                midpoint_controls(),
                cx,
            )
            .unwrap()
            .mechanical
            .motion;
        assert_eq!(hot[0].reference_density_kg_m3().unwrap(), 0.5 * DENSITY);
        let cuts = topology(&geometry, cx);
        let target =
            TetLoadGeometry3::new(cuts.nodes_m(), cuts.tetrahedra(), budget(), 8, cx).unwrap();
        let split = heated
            .split_fracture(&cuts, &target, fracture_controls(), 1e-12, cx)
            .unwrap()
            .motion;
        assert_eq!(split.reference_density_kg_m3(), &[DENSITY; 2]);
        assert!(core::ptr::eq(
            split.reference_density_kg_m3().as_ptr(),
            old.reference_density_kg_m3().as_ptr()
        ));
        close(split.reference_mass_kg(), 2.0, 1e-15);
        assert_eq!(split.accepted_steps(), 1);
        assert_eq!(split.solid().released_faces().len(), 1);
        let cooled = split
            .advance_midpoint_at_temperature(
                0.01,
                &cold,
                "T",
                &[],
                &[],
                settings(),
                midpoint_controls(),
                cx,
            )
            .unwrap();
        close(cooled.total_balance_error_j, 0.0, 1e-12);
        let motion = cooled.mechanical.motion;
        assert_eq!(motion.solid().released_faces().len(), 1);
        assert_eq!(motion.accepted_steps(), 2);
        close(
            motion.solid().fracture_dissipation_j(),
            solid.fracture_dissipation_j(),
            1e-14,
        );
        let snapshot = motion.solid().solution().clone();
        let velocity = motion.velocity_mps().to_vec();
        let time = motion.time_s();
        let prescribed = |shift: f64| {
            let mut u = snapshot.displacement_m.clone();
            for &node in &target.tetrahedra()[1] {
                u[node][2] += shift;
            }
            u.into_iter().flatten().enumerate().collect::<Vec<_>>()
        };
        let refused = motion.advance_midpoint(
            0.01,
            &prescribed(-1e-4),
            &[],
            settings(),
            midpoint_controls(),
            cx,
        );
        assert!(matches!(
            refused,
            Err(BrittleTetDynamicError::Solid(
                BrittleTetError::ReleasedFaceClosure { face: 86, .. }
            ))
        ));
        // Backward Euler has its own accepted-solid path. Allow its explicit
        // numerical kinetic loss so this must refuse at the geometric guard.
        let refused_backward_euler = motion.advance(
            0.01,
            &prescribed(-1e-4),
            &[],
            settings(),
            BrittleTetDynamicControls {
                work: BrittleTetControls {
                    material_energy_tolerance_j: 1e-9,
                    equilibrium_work_tolerance_j: 1e-9,
                    maximum_endpoint_work_defect_j: 1e-6,
                },
                kinetic_work_tolerance_j: 1e-9,
                maximum_numerical_kinetic_loss_j: 1.0,
            },
            cx,
        );
        assert!(matches!(
            refused_backward_euler,
            Err(BrittleTetDynamicError::Solid(
                BrittleTetError::ReleasedFaceClosure { face: 86, .. }
            ))
        ));
        assert_eq!(
            motion.solid().solution().displacement_m,
            snapshot.displacement_m
        );
        assert_eq!(
            motion.solid().solution().element_states,
            snapshot.element_states
        );
        assert_eq!(motion.velocity_mps(), velocity);
        assert_eq!(motion.time_s(), time);
        assert_eq!(motion.remaining_steps(), 2);
        let opened = motion
            .advance_midpoint(
                0.01,
                &prescribed(1e-4),
                &[],
                settings(),
                midpoint_controls(),
                cx,
            )
            .unwrap();
        assert_eq!(opened.motion.accepted_steps(), 3);
        close(opened.motion.solid().fracture_dissipation_j(), GC, 1e-14);
    });
}

#[test]
fn fracture_transfer_refusals_are_atomic_and_cutting_does_not_reset_exhausted_steps() {
    with_cx(false, |cx| {
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, budget(), 8, cx).unwrap();
        let card = source_card();
        let materials = sources(&card, &geometry);
        let cuts = topology(&geometry, cx);
        let target =
            TetLoadGeometry3::new(cuts.nodes_m(), cuts.tetrahedra(), budget(), 8, cx).unwrap();
        let intact = BrittleTetSolid::try_new(&geometry, &materials, budget(), cx).unwrap();
        let intact_motion =
            BrittleTetMotion::try_from_solid(&intact, &[[0.0; 3]; 5], 0.0, 1, cx).unwrap();
        assert!(matches!(
            intact_motion.split_fracture(&cuts, &target, fracture_controls(), 1e-12, cx),
            Err(BrittleTetFractureMotionError::Fracture(_))
        ));
        let solid = failed(&geometry, &materials, cx);
        let start =
            BrittleTetMotion::try_from_solid(&solid, &[[0.003, -0.002, 0.001]; 5], 3.0, 1, cx)
                .unwrap();
        let old = start
            .advance_midpoint(0.01, &[], &[], settings(), midpoint_controls(), cx)
            .unwrap()
            .motion;
        let snapshot = old.solid().solution().clone();
        for tolerance in [-1.0, f64::NAN] {
            assert!(matches!(
                old.split_fracture(&cuts, &target, fracture_controls(), tolerance, cx),
                Err(BrittleTetFractureMotionError::Motion(
                    BrittleTetDynamicError::Invalid(_)
                ))
            ));
        }
        with_cx(true, |cancelled| {
            assert!(matches!(
                old.split_fracture(&cuts, &target, fracture_controls(), 1e-12, cancelled),
                Err(BrittleTetFractureMotionError::Motion(
                    BrittleTetDynamicError::Solid(BrittleTetError::Cancelled)
                ))
            ));
        });
        let a = old
            .split_fracture(&cuts, &target, fracture_controls(), 1e-12, cx)
            .unwrap();
        let b = old
            .split_fracture(&cuts, &target, fracture_controls(), 1e-12, cx)
            .unwrap();
        assert_eq!(a.inventory_after, b.inventory_after);
        assert_eq!(
            a.motion.solid().solution().displacement_m,
            b.motion.solid().solution().displacement_m
        );
        assert_eq!(
            a.motion.solid().solution().element_states,
            b.motion.solid().solution().element_states
        );
        assert_eq!(a.motion.velocity_mps(), b.motion.velocity_mps());
        assert_eq!(a.motion.accepted_steps(), 1);
        assert_eq!(a.motion.remaining_steps(), 0);
        assert_eq!(a.motion.time_s(), old.time_s());
        assert!(matches!(
            a.motion
                .advance_midpoint(0.01, &[], &[], settings(), midpoint_controls(), cx),
            Err(BrittleTetDynamicError::StepBudget {
                accepted_steps: 1,
                maximum_steps: 1
            })
        ));
        assert_eq!(
            old.solid().solution().displacement_m,
            snapshot.displacement_m
        );
        assert_eq!(
            old.solid().solution().element_states,
            snapshot.element_states
        );
    });
}
