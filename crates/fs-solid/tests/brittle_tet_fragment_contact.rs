//! G1/G3: native simultaneous impact of retained brittle fragments.
//! Manufactured source data; exact simplex mass and Gauss inventories provide
//! independent mechanics oracles without treating a tetrahedron as a point mass.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_contact::multicontact::{MultiImpactConfig, MultiImpactDeterminacy};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::brittle::source::dynamics::resolve_brittle_dynamic_material;
use fs_material::brittle::source::*;
use fs_material::cohesive::source::FRACTURE_ENERGY_PER_AREA_DIMS;
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, POISSON_RATIO_PROPERTY, YOUNG_MODULUS_PROPERTY,
};
use fs_qty::{Density, Dims, Pressure};
use fs_solid::brittle::TetCrackBandGeometry;
use fs_solid::brittle::solid::dynamics::fracture::contact::{
    BrittleFaceImpactControls, BrittleFaceImpactError, SolidImpulseControls,
};
use fs_solid::brittle::solid::dynamics::midpoint::BrittleTetMidpointControls;
use fs_solid::brittle::solid::dynamics::{BrittleTetDynamicError, BrittleTetMotion};
use fs_solid::brittle::solid::fracture::BrittleTetFractureControls;
use fs_solid::brittle::solid::{BrittleTetControls, BrittleTetSolid};
use fs_solid::brittle::topology::{FractureTopologyBudget, InteriorSeam, TetFractureTopology};
use fs_solid::linear3::TetAssemblyBudget;
use fs_solid::loads3::TetLoadGeometry3;
use fs_solid::plastic3::IncrementSettings;
use fs_solver::CheckedCgConfig;

#[path = "brittle_tet_fragment_contact/approach.rs"]
mod approach;
#[path = "brittle_tet_fragment_contact/friction.rs"]
mod friction;
#[path = "brittle_tet_fragment_contact/quasistatic.rs"]
mod quasistatic;
#[path = "brittle_tet_fragment_contact/sustained.rs"]
mod sustained;
#[path = "brittle_tet_fragment_contact/time.rs"]
mod time;

const NODES: [[f64; 3]; 5] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, -1.0],
    [0.0, 0.0, 1.0],
];
const TETS: [[usize; 4]; 2] = [[0, 2, 1, 3], [0, 1, 2, 4]];
const SPEED: f64 = 1.0 / 16.0;
const DRIFT: [f64; 2] = [1.0 / 32.0, -1.0 / 64.0];
const APPROACH_DT: f64 = 1.0 / 64.0;
const ORIGIN: [f64; 3] = [0.125, -0.25, 0.5];

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
                seed: 951,
                kernel_id: 951,
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
        maximum_endpoint_work_defect_j: 0.05,
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
        normal_alignment_tolerance: 1e-13,
        traction_tolerance_pa: 1e-10,
        energy_tolerance_j: 1e-12,
        opening_tolerance_m: 1e-12,
    }
}

fn impact_controls() -> BrittleFaceImpactControls {
    let mut impact = MultiImpactConfig::default();
    impact.solver.gradient_tolerance = 1e-12;
    BrittleFaceImpactControls {
        impact,
        impulse: SolidImpulseControls {
            solve: CheckedCgConfig {
                tolerance: 1e-13,
                max_iterations: 64,
                max_corrections: 2,
            },
            maximum_workspace_components: 8192,
            linear_tolerance_kg_m_s: 1e-10,
            angular_tolerance_kg_m2_s: 1e-10,
            work_tolerance_j: 1e-10,
        },
        activation_tolerance_m: 1e-12,
        tangential_tolerance_m: 1e-12,
        velocity_tolerance_mps: 1e-9,
        energy_tolerance_j: 1e-10,
        maximum_contact_rows: 12,
    }
}

fn sources(geometry: &TetLoadGeometry3<'_>) -> Vec<ResolvedBrittleCrackBand> {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, 1024.0),
        (POISSON_RATIO_PROPERTY, Dims::NONE, 0.0),
        (DENSITY_PROPERTY, Density::DIMS, 6.0),
        (BRITTLE_TENSILE_STRENGTH_PROPERTY, Pressure::DIMS, 1.0),
        (BRITTLE_SHEAR_STRENGTH_PROPERTY, Pressure::DIMS, 1.25),
        (
            BRITTLE_MODE_I_FRACTURE_ENERGY_PROPERTY,
            FRACTURE_ENERGY_PER_AREA_DIMS,
            1.0 / 1024.0,
        ),
        (
            BRITTLE_MODE_II_FRACTURE_ENERGY_PROPERTY,
            FRACTURE_ENERGY_PER_AREA_DIMS,
            1.0 / 512.0,
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
                    source: "manufactured brittle-fragment contact fixture".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured brittle solid".into(),
            phase: "solid".into(),
            process: "native coupled contact verification".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    (0..geometry.tetrahedra().len())
        .map(|element| {
            let band =
                TetCrackBandGeometry::from_element(geometry, element, [0.0, 0.0, 1.0]).unwrap();
            resolve_brittle_dynamic_material(
                &card,
                &QueryPoint::new(),
                MaterialPropertySelection::SingleClaimOnly,
                band.material_settings(0.03, 1e-9),
            )
            .unwrap()
        })
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

fn failed<'a>(
    geometry: &'a TetLoadGeometry3<'a>,
    materials: &'a [ResolvedBrittleCrackBand],
    cx: &Cx<'_>,
) -> BrittleTetSolid<'a> {
    let displacement: Vec<_> = geometry
        .nodes_m()
        .iter()
        .map(|x| [0.0, 0.0, x[2] / 128.0])
        .collect();
    let solid = BrittleTetSolid::try_new(geometry, materials, budget(), cx)
        .unwrap()
        .advance(
            &prescribed(&displacement),
            &[],
            settings(),
            work_controls(),
            cx,
        )
        .unwrap()
        .solid;
    assert!(solid.solution().element_states.iter().all(|s| s.failed()));
    assert_eq!(solid.stored_energy_j(), 0.0);
    solid
}

fn topology<'g>(
    geometry: &'g TetLoadGeometry3<'_>,
    seams: &[InteriorSeam],
    cx: &Cx<'_>,
) -> TetFractureTopology<'g> {
    TetFractureTopology::new(
        geometry,
        seams,
        FractureTopologyBudget {
            maximum_output_nodes: 12,
            maximum_face_entries: 16,
            maximum_seams: 2,
        },
        cx,
    )
    .unwrap()
}

// Start at a finite opening and take one actual constant-velocity midpoint
// interval to exactly touching pairs. Thus impact must retain a USED physical
// step budget, not merely a fresh initializer's zero count.
fn arriving<'a>(solid: &BrittleTetSolid<'a>, cx: &Cx<'_>) -> BrittleTetMotion<'a> {
    let mut opened = solid.solution().displacement_m.clone();
    let mut target = opened.clone();
    let mut velocity = vec![[0.0; 3]; opened.len()];
    for (element, tet) in solid.geometry().tetrahedra().iter().enumerate() {
        let approach = if element % 2 == 0 { SPEED } else { -SPEED };
        for &node in tet {
            velocity[node] = [DRIFT[0], DRIFT[1], approach];
            opened[node][2] = solid.solution().displacement_m[node][2] - approach * APPROACH_DT;
            target[node][0] = solid.solution().displacement_m[node][0] + DRIFT[0] * APPROACH_DT;
            target[node][1] = solid.solution().displacement_m[node][1] + DRIFT[1] * APPROACH_DT;
        }
    }
    let open = solid
        .advance(&prescribed(&opened), &[], settings(), work_controls(), cx)
        .unwrap()
        .solid;
    let motion = BrittleTetMotion::try_from_solid(&open, &velocity, 2.0, 2, cx).unwrap();
    let arrived = motion
        .advance_midpoint(
            APPROACH_DT,
            &prescribed(&target),
            &[],
            settings(),
            midpoint_controls(),
            cx,
        )
        .unwrap()
        .motion;
    assert_eq!(arrived.velocity_mps(), velocity);
    assert_eq!(arrived.accepted_steps(), 1);
    assert_eq!(arrived.remaining_steps(), 1);
    arrived
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:e} != {b:e}");
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn inventory(motion: &BrittleTetMotion<'_>) -> ([f64; 3], [f64; 3], f64) {
    // Positive degree-two tetrahedral quadrature integrates position cross
    // velocity and squared velocity, including all consistent-mass cross terms.
    let a = 0.5854101966249685;
    let b = 0.1381966011250105;
    let geometry = motion.solid().geometry();
    let mut linear = [0.0; 3];
    let mut angular = [0.0; 3];
    let mut kinetic = 0.0;
    for (element, tet) in geometry.tetrahedra().iter().enumerate() {
        let weight = motion.reference_density_kg_m3()[element]
            * geometry.element_volumes_m3()[element]
            / 4.0;
        for q in 0..4 {
            let mut position = [0.0; 3];
            let mut velocity = [0.0; 3];
            for (i, &node) in tet.iter().enumerate() {
                let shape = if i == q { a } else { b };
                for c in 0..3 {
                    position[c] += shape
                        * (geometry.nodes_m()[node][c]
                            + motion.solid().solution().displacement_m[node][c]
                            - ORIGIN[c]);
                    velocity[c] += shape * motion.velocity_mps()[node][c];
                }
            }
            let moment = cross(position, velocity);
            for c in 0..3 {
                linear[c] += weight * velocity[c];
                angular[c] += weight * moment[c];
                kinetic += 0.5 * weight * velocity[c] * velocity[c];
            }
        }
    }
    (linear, angular, kinetic)
}

fn retained(before: &BrittleTetMotion<'_>, after: &BrittleTetMotion<'_>) {
    assert_eq!(
        before.solid().solution().displacement_m,
        after.solid().solution().displacement_m
    );
    assert_eq!(
        before.solid().solution().element_strain,
        after.solid().solution().element_strain
    );
    assert_eq!(
        before.solid().solution().element_stress,
        after.solid().solution().element_stress
    );
    assert_eq!(
        before.solid().solution().element_states,
        after.solid().solution().element_states
    );
    assert_eq!(
        before.solid().released_faces(),
        after.solid().released_faces()
    );
    assert!(core::ptr::eq(
        before.solid().materials().as_ptr(),
        after.solid().materials().as_ptr()
    ));
    assert_eq!(
        before.solid().stored_energy_j(),
        after.solid().stored_energy_j()
    );
    assert_eq!(
        before.solid().fracture_dissipation_j(),
        after.solid().fracture_dissipation_j()
    );
    assert_eq!(
        before.reference_density_kg_m3(),
        after.reference_density_kg_m3()
    );
    assert_eq!(before.reference_mass_kg(), after.reference_mass_kg());
    assert_eq!(before.time_s().to_bits(), after.time_s().to_bits());
    assert_eq!(before.accepted_steps(), after.accepted_steps());
    assert_eq!(before.remaining_steps(), after.remaining_steps());
}

#[test]
fn simultaneous_restitution_matches_consistent_mass_and_continues_without_consuming_an_event_step()
{
    with_cx(false, |cx| {
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, budget(), 8, cx).unwrap();
        let material = sources(&geometry);
        let solid = failed(&geometry, &material, cx);
        let cut = topology(
            &geometry,
            &[InteriorSeam {
                id: 951,
                nodes: [0, 1, 2],
                side_a_element: 0,
            }],
            cx,
        );
        let target =
            TetLoadGeometry3::new(cut.nodes_m(), cut.tetrahedra(), budget(), 8, cx).unwrap();
        let split = solid
            .split_fracture(&cut, &target, fracture_controls(), cx)
            .unwrap();
        let motion = arriving(&split.solid, cx);
        let before = inventory(&motion);
        close(
            before.2,
            SPEED * SPEED + DRIFT[0].powi(2) + DRIFT[1].powi(2),
            1e-15,
        );
        for e in [0.0, 0.5, 1.0] {
            let hit = motion
                .impact_released_faces(&[e], &[], ORIGIN, impact_controls(), cx)
                .unwrap();
            retained(&motion, &hit.motion);
            assert_eq!(hit.rows.len(), 3);
            assert!(matches!(
                hit.contact.determinacy,
                MultiImpactDeterminacy::UniqueImpulses
            ));
            // Each tetrahedron has m=1 and scalar M=(I+ones)/20.
            // Its inverse is 20I-4ones. The three paired face rows therefore
            // have W=40I-8ones, whose uniform-mode eigenvalue is 16.
            let impulse = (1.0 + e) * SPEED / 8.0;
            for (&p, &speed) in hit
                .contact
                .normal_impulses
                .iter()
                .zip(&hit.contact.relative_normal_speeds_after)
            {
                close(p, impulse, 1e-12);
                close(speed, 2.0 * e * SPEED, 1e-11);
            }
            for (element, tet) in target.tetrahedra().iter().enumerate() {
                let sign = if element == 0 { 1.0 } else { -1.0 };
                for &node in tet {
                    let expected = if target.nodes_m()[node][2] == 0.0 {
                        -sign * e * SPEED
                    } else {
                        sign * (2.5 + 1.5 * e) * SPEED
                    };
                    close(hit.motion.velocity_mps()[node][2], expected, 1e-11);
                    assert_eq!(hit.motion.velocity_mps()[node][0], DRIFT[0]);
                    assert_eq!(hit.motion.velocity_mps()[node][1], DRIFT[1]);
                }
            }
            let after = inventory(&hit.motion);
            for c in 0..3 {
                close(after.0[c], before.0[c], 1e-12);
                close(after.1[c], before.1[c], 1e-12);
            }
            let expected_loss = 3.0 * (1.0 - e * e) * SPEED * SPEED / 8.0;
            close(hit.contact.dissipated_energy, expected_loss, 1e-12);
            close(after.2, before.2 - expected_loss, 1e-12);
            close(hit.motion.kinetic_energy_j(), after.2, 1e-12);
            close(hit.impulse.kinetic_change_j, -expected_loss, 1e-12);
            assert_eq!(hit.impulse.support_work_j, 0.0);
            assert_eq!(hit.absolute_contact_couple_kg_m2_s, 0.0);
            close(hit.relative_work_error_j, 0.0, 1e-12);
            close(hit.energy_balance_error_j, 0.0, 1e-12);
            // At e=0 the faces stop approaching but internal velocity remains;
            // a rigid point-mass shortcut would incorrectly remove all of it.
            if e == 0.0 {
                close(
                    after.2 - DRIFT[0].powi(2) - DRIFT[1].powi(2),
                    5.0 * SPEED * SPEED / 8.0,
                    1e-12,
                );
            }
            let next = hit
                .motion
                .advance_midpoint(1.0 / 1024.0, &[], &[], settings(), midpoint_controls(), cx)
                .unwrap();
            assert_eq!(next.motion.accepted_steps(), 2);
            assert_eq!(next.motion.remaining_steps(), 0);
            assert_eq!(
                next.motion.solid().fracture_dissipation_j(),
                motion.solid().fracture_dissipation_j()
            );
            assert_eq!(
                next.motion.solid().released_faces(),
                motion.solid().released_faces()
            );
            assert!(matches!(
                next.motion.advance_midpoint(
                    1.0 / 1024.0,
                    &[],
                    &[],
                    settings(),
                    midpoint_controls(),
                    cx
                ),
                Err(BrittleTetDynamicError::StepBudget { .. })
            ));
        }
    });
}

#[test]
fn redundant_shared_edge_contacts_determine_native_velocity_without_unique_row_impulses() {
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
        let material = sources(&geometry);
        let solid = failed(&geometry, &material, cx);
        let cut = topology(
            &geometry,
            &[
                InteriorSeam {
                    id: 951,
                    nodes: [0, 1, 2],
                    side_a_element: 0,
                },
                InteriorSeam {
                    id: 952,
                    nodes: [0, 2, 3],
                    side_a_element: 2,
                },
            ],
            cx,
        );
        let target =
            TetLoadGeometry3::new(cut.nodes_m(), cut.tetrahedra(), budget(), 16, cx).unwrap();
        let split = solid
            .split_fracture(&cut, &target, fracture_controls(), cx)
            .unwrap();
        let motion = arriving(&split.solid, cx);
        let e = 0.5;
        let hit = motion
            .impact_released_faces(&[e, e], &[], ORIGIN, impact_controls(), cx)
            .unwrap();
        retained(&motion, &hit.motion);
        assert_eq!(hit.rows.len(), 6);
        let MultiImpactDeterminacy::UniqueVelocity(witness) = &hit.contact.determinacy else {
            panic!("two repeated shared-edge rows must not acquire unique impulse authority");
        };
        assert_eq!(witness.dimension, 6);
        assert_eq!(witness.rank, 4);
        let mut physical_impulse = vec![0.0; target.nodes_m().len()];
        for (row, &p) in hit.rows.iter().zip(&hit.contact.normal_impulses) {
            physical_impulse[row.nodes[1]] += p;
        }
        let d = (1.0 + e) * SPEED;
        for node in 0..target.nodes_m().len() {
            if physical_impulse[node] == 0.0 {
                continue;
            }
            let p = target.nodes_m()[node];
            // Sum M*delta_v from the TWO original tetrahedra in this body.
            // Shared diagonal vertices receive twice a corner's impulse.
            let expected = if p[0] == p[1] { d / 4.0 } else { d / 8.0 };
            close(physical_impulse[node], expected, 1e-10);
        }
        for (element, tet) in target.tetrahedra().iter().enumerate() {
            let sign = if element % 2 == 0 { 1.0 } else { -1.0 };
            for &node in tet {
                let expected = if target.nodes_m()[node][2] == 0.0 {
                    -sign * e * SPEED
                } else {
                    sign * (2.5 + 1.5 * e) * SPEED
                };
                close(hit.motion.velocity_mps()[node][2], expected, 1e-9);
            }
        }
        let before = inventory(&motion);
        let after = inventory(&hit.motion);
        for c in 0..3 {
            close(after.0[c], before.0[c], 1e-10);
            close(after.1[c], before.1[c], 1e-10);
        }
        let loss = 3.0 * (1.0 - e * e) * SPEED * SPEED / 4.0;
        close(hit.contact.dissipated_energy, loss, 1e-10);
        close(after.2, before.2 - loss, 1e-10);
        close(hit.energy_balance_error_j, 0.0, 1e-10);
    });
}

#[test]
fn held_moving_fragment_supplies_explicit_support_work_through_the_same_mass_operator() {
    with_cx(false, |cx| {
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, budget(), 8, cx).unwrap();
        let material = sources(&geometry);
        let solid = failed(&geometry, &material, cx);
        let cut = topology(
            &geometry,
            &[InteriorSeam {
                id: 951,
                nodes: [0, 1, 2],
                side_a_element: 0,
            }],
            cx,
        );
        let target =
            TetLoadGeometry3::new(cut.nodes_m(), cut.tetrahedra(), budget(), 8, cx).unwrap();
        let split = solid
            .split_fracture(&cut, &target, fracture_controls(), cx)
            .unwrap();
        let motion = arriving(&split.solid, cx);
        let mut held: Vec<_> = target.tetrahedra()[0]
            .iter()
            .flat_map(|node| (0..3).map(move |c| 3 * node + c))
            .collect();
        held.sort_unstable();
        let e = 0.5;
        let hit = motion
            .impact_released_faces(&[e], &held, ORIGIN, impact_controls(), cx)
            .unwrap();
        retained(&motion, &hit.motion);
        for &dof in &held {
            assert_eq!(
                hit.motion.velocity_mps().as_flattened()[dof].to_bits(),
                motion.velocity_mps().as_flattened()[dof].to_bits()
            );
        }
        // Only B responds: W=20I-4ones, uniform-mode eigenvalue eight.
        let p = (1.0 + e) * SPEED / 4.0;
        for (&actual, &speed) in hit
            .contact
            .normal_impulses
            .iter()
            .zip(&hit.contact.relative_normal_speeds_after)
        {
            close(actual, p, 1e-12);
            close(speed, 2.0 * e * SPEED, 1e-11);
        }
        let support_work = 3.0 * p * SPEED;
        let loss = 3.0 * (1.0 - e * e) * SPEED * SPEED / 4.0;
        close(hit.impulse.support_work_j, support_work, 1e-12);
        close(hit.contact.dissipated_energy, loss, 1e-12);
        close(hit.impulse.kinetic_change_j, support_work - loss, 1e-12);
        assert!(
            hit.impulse.kinetic_change_j > 0.0,
            "prescribed support motion must supply its measured work"
        );
        let before = inventory(&motion);
        let after = inventory(&hit.motion);
        for c in 0..3 {
            close(
                after.0[c] - before.0[c],
                hit.impulse.support.linear_kg_m_s[c],
                1e-12,
            );
            close(
                after.1[c] - before.1[c],
                hit.impulse.support.angular_kg_m2_s[c],
                1e-12,
            );
        }
        close(hit.impulse.support.linear_kg_m_s[2], 3.0 * p, 1e-12);
        close(after.2 - before.2, support_work - loss, 1e-12);
        close(hit.relative_work_error_j, 0.0, 1e-12);
        close(hit.energy_balance_error_j, 0.0, 1e-12);
    });
}

#[test]
fn unsupported_contact_geometry_budgets_and_cancellation_leave_deterministic_retry() {
    with_cx(false, |cx| {
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, budget(), 8, cx).unwrap();
        let material = sources(&geometry);
        let solid = failed(&geometry, &material, cx);
        let cut = topology(
            &geometry,
            &[InteriorSeam {
                id: 951,
                nodes: [0, 1, 2],
                side_a_element: 0,
            }],
            cx,
        );
        let target =
            TetLoadGeometry3::new(cut.nodes_m(), cut.tetrahedra(), budget(), 8, cx).unwrap();
        let split = solid
            .split_fracture(&cut, &target, fracture_controls(), cx)
            .unwrap();
        let motion = arriving(&split.solid, cx);
        let saved_velocity = motion.velocity_mps().to_vec();
        let saved_history = motion.solid().solution().element_states.clone();
        for offset in [
            [1e-4, 0.0, 0.0],
            [0.0, 0.0, -2.0_f64.powi(-45)],
            [0.0, 0.0, 1e-4],
        ] {
            let mut displacement = motion.solid().solution().displacement_m.clone();
            for &node in &target.tetrahedra()[1] {
                for c in 0..3 {
                    displacement[node][c] += offset[c];
                }
            }
            let shifted = motion
                .solid()
                .advance(
                    &prescribed(&displacement),
                    &[],
                    settings(),
                    work_controls(),
                    cx,
                )
                .unwrap()
                .solid;
            let candidate = BrittleTetMotion::try_from_solid(
                &shifted,
                motion.velocity_mps(),
                motion.time_s(),
                1,
                cx,
            )
            .unwrap();
            let result =
                candidate.impact_released_faces(&[0.5], &[], ORIGIN, impact_controls(), cx);
            if offset[2] > 0.0 {
                assert!(matches!(result, Err(BrittleFaceImpactError::Invalid(_))));
            } else {
                assert!(matches!(
                    result,
                    Err(BrittleFaceImpactError::Geometry { face: 951, .. })
                ));
            }
            assert_eq!(candidate.velocity_mps(), saved_velocity);
            assert_eq!(candidate.solid().solution().displacement_m, displacement);
        }
        for restitution in [&[][..], &[1.1][..], &[f64::NAN][..]] {
            assert!(
                motion
                    .impact_released_faces(restitution, &[], ORIGIN, impact_controls(), cx)
                    .is_err()
            );
        }
        for invalid in [
            BrittleFaceImpactControls {
                maximum_contact_rows: 2,
                ..impact_controls()
            },
            BrittleFaceImpactControls {
                activation_tolerance_m: -1.0,
                ..impact_controls()
            },
            BrittleFaceImpactControls {
                energy_tolerance_j: f64::NAN,
                ..impact_controls()
            },
        ] {
            assert!(
                motion
                    .impact_released_faces(&[0.5], &[], ORIGIN, invalid, cx)
                    .is_err()
            );
        }
        assert!(
            motion
                .impact_released_faces(&[0.5], &[0, 0], ORIGIN, impact_controls(), cx)
                .is_err()
        );
        with_cx(true, |cancelled| {
            assert!(
                motion
                    .impact_released_faces(&[0.5], &[], ORIGIN, impact_controls(), cancelled)
                    .is_err()
            );
        });
        assert_eq!(motion.velocity_mps(), saved_velocity);
        assert_eq!(motion.solid().solution().element_states, saved_history);
        assert_eq!(motion.accepted_steps(), 1);
        let first = motion
            .impact_released_faces(&[0.5], &[], ORIGIN, impact_controls(), cx)
            .unwrap();
        let retry = motion
            .impact_released_faces(&[0.5], &[], ORIGIN, impact_controls(), cx)
            .unwrap();
        assert_eq!(first.motion.velocity_mps(), retry.motion.velocity_mps());
        assert_eq!(first.contact, retry.contact);
        retained(&motion, &first.motion);
        retained(&motion, &retry.motion);
    });
}
