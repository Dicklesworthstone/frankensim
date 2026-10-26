//! G1/G3/G4: actual nonlinear/circuit/induction fields feed native mechanics.
//! Source coefficients are synthetic, not measured actuator material data.

#[path = "../nonlinear_eddy/common.rs"]
mod source;

use super::{AMBIENT_MU, DEPTH, common, force_limits, frame, near, specimen};
use fs_em::{
    planar::{
        PlanarMagnetostatics,
        eddy::{
            EddyConductor, EddyControls, PlanarEddy,
            nonlinear::{NonlinearEddyError, PlanarNonlinearEddy},
        },
        nonlinear::{IsotropicMagneticMaterial, PlanarNonlinearMagnetostatics},
    },
    winding::{
        nonlinear::{
            NonlinearWindingCircuit, NonlinearWindingConfig, NonlinearWindingInitial,
            NonlinearWindings,
        },
        ohmic::WindingConductor,
    },
};
use fs_ga::{Mat3, Se3, Vec3};
use fs_matdb::{ClaimSet, MaterialCard, MaterialStateId};
use fs_material::{
    conductor::{ElectricalPropertyBasis, ResolvedConductorState},
    field_storage::{MAGNETIC_PERMEABILITY_DIMS, MAGNETIC_PERMEABILITY_PROPERTY},
    state_point::{
        MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
        resolve_material_state_point,
    },
};
use fs_mbd::articulated::{
    ArticulatedModel, FreeFloatingBaseState, JointModel, Link, SpatialInertia,
    free_floating_forward_dynamics,
};
use fs_phs::circuit::{Branch, CircuitGraph};
use fs_rep_mesh::TriComplex2;

#[path = "../nonlinear_winding_circuit/common.rs"]
mod winding;

fn geometry(n: usize) -> (TriComplex2, Vec<usize>, Vec<usize>, Vec<f64>) {
    let (seed, enclosed, body) = specimen(n, AMBIENT_MU, false, DEPTH);
    let mesh = seed.mesh().clone();
    let mut profile = vec![0.0; mesh.faces().len()];
    for &cell in &body {
        profile[cell] = 4.0;
    }
    // Equal area return away from the enclosed body's contour. Both legs are
    // explicit original mesh regions, so axial current closes without a sink.
    for (cell, face) in mesh.faces().iter().enumerate() {
        let center: [f64; 2] = [0, 1].map(|axis| {
            face.iter()
                .map(|&v| mesh.vertices()[v as usize][axis])
                .sum::<f64>()
                / 3.0
        });
        if center[0] > 0.75 && center[1].abs() < 0.5 {
            profile[cell] = -4.0;
        }
    }
    (mesh, enclosed, body, profile)
}
fn potential([x, y]: [f64; 2], gradient: f64) -> f64 {
    0.3 * x + gradient * (x * x - y * y)
}
fn boundary(mesh: &TriComplex2, gradient: f64) -> Vec<(usize, f64)> {
    mesh.boundary_trace()
        .unwrap()
        .vertices
        .into_iter()
        .map(|v| (v as usize, potential(mesh.vertices()[v as usize], gradient)))
        .collect()
}
fn problem(
    mesh: TriComplex2,
    profile: &[f64],
    card: &MaterialCard,
) -> PlanarNonlinearMagnetostatics {
    let nonlinear = IsotropicMagneticMaterial::Reversible(source::law(card));
    let ambient = IsotropicMagneticMaterial::Linear(common::material(1.0));
    let materials = profile
        .iter()
        .map(|p| {
            if *p == 0.0 {
                ambient.clone()
            } else {
                nonlinear.clone()
            }
        })
        .collect();
    PlanarNonlinearMagnetostatics::new(mesh, materials, common::limits(), &mut || false).unwrap()
}
fn history<'m>(
    problem: &'m PlanarNonlinearMagnetostatics,
    profile: &[f64],
    current: f64,
    conductors: Vec<EddyConductor>,
    gradient: f64,
) -> PlanarNonlinearEddy<'m> {
    let j: Vec<_> = profile.iter().map(|p| p * current).collect();
    let bc = boundary(problem.mesh(), gradient);
    let initial: Vec<_> = problem
        .mesh()
        .vertices()
        .iter()
        .map(|&p| potential(p, gradient))
        .collect();
    let field = problem
        .solve(&j, &bc, &initial, source::controls().nonlinear, &mut || {
            false
        })
        .unwrap();
    let insulators: Vec<_> = (0..profile.len())
        .filter(|&cell| !conductors.iter().any(|s| s.triangle == cell))
        .collect();
    PlanarNonlinearEddy::new_with_impressed_current(
        &field,
        conductors,
        &insulators,
        &bc,
        &j,
        0.0,
        source::controls(),
        &mut || false,
    )
    .unwrap()
}
fn coil_sources(card: &MaterialCard, profile: &[f64]) -> Vec<WindingConductor> {
    let conductor = source::conductor(card, 300.0, ElectricalPropertyBasis::Resistivity);
    profile
        .iter()
        .enumerate()
        .filter(|(_, p)| **p != 0.0)
        .map(|(triangle, _)| WindingConductor {
            triangle,
            conductor: conductor.clone(),
        })
        .collect()
}
fn circuit<'m>(
    problem: &'m PlanarNonlinearMagnetostatics,
    profile: &[f64],
    card: &MaterialCard,
    voltage: f64,
) -> NonlinearWindingCircuit<'m> {
    NonlinearWindingCircuit::new(
        history(problem, profile, 0.0, vec![], 0.1),
        CircuitGraph {
            node_count: 3,
            branches: vec![
                (1, 0, Branch::VoltageSource { port: 0 }),
                (1, 2, Branch::Inductor { henries: 123.0 }),
                (2, 0, Branch::Resistor { ohms: 1.0 }),
            ],
            transformers: vec![],
        },
        NonlinearWindings {
            branch_indices: vec![1],
            profiles_per_m2: vec![profile.to_vec()],
            end_inductance_h: vec![0.0],
        },
        coil_sources(card, profile),
        NonlinearWindingInitial {
            currents_a: vec![0.0],
            remaining_fluxes_wb: vec![],
            charges_c: vec![],
            external_inputs: vec![voltage],
        },
        NonlinearWindingConfig {
            max_windings: 1,
            max_profile_values: profile.len(),
            current_balance_tolerance_turns: 1e-10,
        },
        &mut || false,
    )
    .unwrap()
}

#[test]
fn linear_history_limit_matches_existing_force_and_independent_volume_lorentz_refinement() {
    let mut previous = f64::INFINITY;
    for n in [8, 16, 32] {
        let (mesh, enclosed, body, profile) = geometry(n);
        let material = common::material(1.0);
        let linear = PlanarMagnetostatics::new(
            mesh.clone(),
            vec![material.clone(); profile.len()],
            common::limits(),
            &mut || false,
        )
        .unwrap();
        let nonlinear = PlanarNonlinearMagnetostatics::new(
            mesh,
            vec![IsotropicMagneticMaterial::Linear(material); profile.len()],
            common::limits(),
            &mut || false,
        )
        .unwrap();
        let j: Vec<_> = profile.iter().map(|p| 0.5 * p).collect();
        let bc = boundary(linear.mesh(), 0.0);
        let reference = linear
            .solve(&j, &bc, common::config(), &mut || false)
            .unwrap();
        let induction = PlanarEddy::new_with_impressed_current(
            &reference,
            vec![],
            &(0..profile.len()).collect::<Vec<_>>(),
            &bc,
            &j,
            0.0,
            EddyControls {
                linear: common::config(),
                energy_tolerance_j: 1e-9,
                current_balance_tolerance_a: 1e-8,
            },
            &mut || false,
        )
        .unwrap();
        let evolving = history(&nonlinear, &profile, 0.5, vec![], 0.0);
        let origin = [0.0, 0.0, DEPTH / 2.0];
        let original = reference
            .maxwell_wrench(&enclosed, &frame(), origin, force_limits(), &mut || false)
            .unwrap();
        for load in [
            induction
                .maxwell_wrench(&enclosed, &frame(), origin, force_limits(), &mut || false)
                .unwrap(),
            evolving
                .maxwell_wrench(&enclosed, &frame(), origin, force_limits(), &mut || false)
                .unwrap(),
        ] {
            for (actual, expected) in load.to_array().into_iter().zip(original.to_array()) {
                near(actual, expected);
            }
        }
        // Independent volume J cross B, not a boundary stress or trace call.
        // In this uniform linear medium the body's net magnetic force is its
        // integrated Lorentz load. Both discrete formulations converge to it.
        let mut lorentz = [0.0; 2];
        for cell in body {
            let volume = linear.mesh().face_measure(cell).unwrap();
            let b = reference.magnetic_flux_density_t()[cell];
            lorentz[0] -= j[cell] * b[1] * volume;
            lorentz[1] += j[cell] * b[0] * volume;
        }
        let error = (original.force.x - lorentz[0]).hypot(original.force.y - lorentz[1]);
        assert!(
            error < 0.65 * previous,
            "n={n}, error={error}, previous={previous}"
        );
        previous = error;
        if n == 32 {
            assert!(error < 0.05 * lorentz[0].hypot(lorentz[1]));
        }
        assert_eq!(evolving.time_s(), 0.0);
    }
}

#[test]
fn source_saturation_and_circuit_voltage_change_native_force_and_body_acceleration() {
    let mechanics = ArticulatedModel::new(vec![Link::new(
        "body",
        None,
        Se3::identity(),
        JointModel::FIXED,
        SpatialInertia::new(
            2.0,
            Vec3::new(0.0, 0.0, 0.0),
            Mat3 {
                m: [0.2, 0.0, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0, 0.2],
            },
        )
        .unwrap(),
    )])
    .unwrap();
    let mut loads = Vec::new();
    for saturating in [false, true] {
        let card = source::source(
            if saturating {
                &[(0.0, 0.0), (0.1, 0.025), (4.0, 3.925)]
            } else {
                &[(0.0, 0.0), (4.0, 1.0)]
            },
            8.0,
        );
        let (mesh, enclosed, _, profile) = geometry(16);
        let problem = problem(mesh, &profile, &card);
        for voltage in [0.4, 0.8] {
            let original = circuit(&problem, &profile, &card, voltage);
            let old_a = original.field().potential_wb_m().to_vec();
            let accepted = original
                .step_to_candidate(
                    0.1,
                    vec![],
                    original.conductors().to_vec(),
                    &[voltage],
                    winding::controls(),
                    &mut || false,
                )
                .unwrap();
            let load = accepted
                .circuit()
                .field()
                .maxwell_wrench(
                    &enclosed,
                    &frame(),
                    [0.0, 0.0, DEPTH / 2.0],
                    force_limits(),
                    &mut || false,
                )
                .unwrap();
            let dynamics = free_floating_forward_dynamics(
                &mechanics,
                FreeFloatingBaseState::stationary(Se3::identity()),
                &[],
                &[],
                &[],
                Vec3::new(0.0, 0.0, 0.0),
                &[load],
            )
            .unwrap();
            near(
                dynamics.base_spatial_acceleration_body.linear.x,
                load.force.x / 2.0,
            );
            near(
                dynamics.base_spatial_acceleration_body.linear.y,
                load.force.y / 2.0,
            );
            near(
                dynamics.base_spatial_acceleration_body.angular.z,
                load.torque.z / 0.2,
            );
            assert!(accepted.circuit().currents_a()[0] > 0.0);
            assert_eq!(accepted.circuit().time_s(), 0.1);
            assert_eq!(original.field().potential_wb_m(), old_a);
            assert_eq!(original.time_s(), 0.0);
            loads.push((load.force.x, accepted.circuit().currents_a()[0]));
        }
    }
    assert!((loads[0].0 - loads[1].0).abs() > 1e-5);
    assert!((loads[2].0 - loads[3].0).abs() > 1e-5);
    assert!((loads[0].0 - loads[2].0).abs() > 1e-5);
    assert!((loads[0].1 - loads[2].1).abs() > 1e-5);
}

#[test]
fn decaying_eddy_history_keeps_its_actual_force_instead_of_resolving_static_equilibrium() {
    let card = source::source(&[(0.0, 0.0), (0.1, 0.025), (4.0, 3.925)], 8.0);
    let (mesh, enclosed, _, profile) = geometry(16);
    let problem = problem(mesh, &profile, &card);
    let conductors: Vec<_> = coil_sources(&card, &profile)
        .into_iter()
        .map(|s| EddyConductor {
            triangle: s.triangle,
            region: if profile[s.triangle] > 0.0 { 0 } else { 1 },
            conductor: s.conductor,
        })
        .collect();
    let original = history(&problem, &profile, 0.5, conductors, 0.0);
    let next = original
        .step_to_candidate(
            0.01,
            original.conductors().to_vec(),
            &vec![0.0; profile.len()],
            source::controls(),
            &mut || false,
        )
        .unwrap();
    assert!(next.step().heat_j > 0.0);
    let actual = next
        .state()
        .maxwell_wrench(
            &enclosed,
            &frame(),
            [0.0, 0.0, DEPTH / 2.0],
            force_limits(),
            &mut || false,
        )
        .unwrap();
    let static_replacement = history(&problem, &profile, 0.0, vec![], 0.0);
    let replaced = static_replacement
        .maxwell_wrench(
            &enclosed,
            &frame(),
            [0.0, 0.0, DEPTH / 2.0],
            force_limits(),
            &mut || false,
        )
        .unwrap();
    assert!((actual.force.x - replaced.force.x).abs() > 1e-5);
    let saved = next.state().potential_wb_m().to_vec();
    let mut checkpoints = 0;
    next.state()
        .maxwell_wrench(
            &enclosed,
            &frame(),
            [0.0, 0.0, DEPTH / 2.0],
            force_limits(),
            &mut || {
                checkpoints += 1;
                false
            },
        )
        .unwrap();
    let mut seen = 0;
    assert!(matches!(
        next.state().maxwell_wrench(
            &enclosed,
            &frame(),
            [0.0, 0.0, DEPTH / 2.0],
            force_limits(),
            &mut || {
                seen += 1;
                seen == checkpoints
            }
        ),
        Err(NonlinearEddyError::Cancelled)
    ));
    assert_eq!(next.state().potential_wb_m(), saved);
    assert_eq!(next.state().time_s(), 0.01);
    assert_eq!(original.time_s(), 0.0);
}

#[test]
fn conducting_nonlinear_or_impressed_contours_refuse_even_when_recovered_b_is_finite() {
    let (mesh, enclosed, body, profile) = geometry(8);
    let card = source::source(&[(0.0, 0.0), (4.0, 1.0)], 8.0);
    let nonlinear = problem(mesh.clone(), &profile, &card);
    let original = history(&nonlinear, &profile, 0.5, vec![], 0.0);
    assert!(
        original
            .maxwell_wrench(&body, &frame(), [0.0; 3], force_limits(), &mut || false)
            .is_err()
    );
    assert!(
        original
            .maxwell_wrench(
                &enclosed,
                &frame(),
                [0.0; 3],
                fs_em::force::MaxwellForceLimits {
                    max_cells: 1,
                    max_edges: 1
                },
                &mut || false
            )
            .is_err()
    );

    // The entire exterior is sourced linear material, with an explicit eddy
    // conductor touching the otherwise valid contour. Zero instantaneous J or
    // a finite B sample must not erase that conductor's physical support.
    let mu = common::source(
        fs_matdb::PropertyValue::Scalar {
            value: 1.0,
            dims: MAGNETIC_PERMEABILITY_DIMS,
        },
        fs_evidence::ValidityDomain::unconstrained(),
    );
    let mut claims = ClaimSet::new();
    for (_, claim) in mu
        .claims()
        .claims_ordered()
        .chain(card.claims().claims_ordered())
    {
        claims.insert_claim(claim.clone()).unwrap();
    }
    let linear_card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic ambient".into(),
            phase: "solid".into(),
            process: "test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let snapshot = resolve_material_state_point(
        &linear_card,
        &source::point(300.0),
        &[
            ScalarPropertyRequirement::try_new(
                MAGNETIC_PERMEABILITY_PROPERTY,
                MAGNETIC_PERMEABILITY_DIMS,
                ScalarAdmissibility::StrictlyPositive,
            )
            .unwrap(),
            ElectricalPropertyBasis::Resistivity.requirement().unwrap(),
        ],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let conductor = ResolvedConductorState::try_from_material_with_basis(
        snapshot.clone(),
        ElectricalPropertyBasis::Resistivity,
    )
    .unwrap();
    let linear = PlanarNonlinearMagnetostatics::new(
        mesh,
        vec![IsotropicMagneticMaterial::Linear(snapshot); profile.len()],
        common::limits(),
        &mut || false,
    )
    .unwrap();
    let boundary_cell = linear
        .mesh()
        .trace_for_faces(enclosed.iter().copied())
        .unwrap()
        .edges[0]
        .source_face;
    let conducting = history(
        &linear,
        &profile,
        0.0,
        vec![EddyConductor {
            triangle: boundary_cell,
            region: 0,
            conductor,
        }],
        0.0,
    );
    assert!(
        conducting
            .maxwell_wrench(&enclosed, &frame(), [0.0; 3], force_limits(), &mut || false)
            .is_err()
    );
    let mut imposed = vec![0.0; profile.len()];
    imposed[boundary_cell] = 1.0;
    let another = (0..profile.len()).find(|&i| i != boundary_cell).unwrap();
    imposed[another] = -1.0;
    let impressed = history(&linear, &imposed, 1.0, vec![], 0.0);
    assert!(
        impressed
            .maxwell_wrench(&enclosed, &frame(), [0.0; 3], force_limits(), &mut || false)
            .is_err()
    );
    assert_eq!(conducting.time_s(), 0.0);
    assert_eq!(impressed.time_s(), 0.0);
}
