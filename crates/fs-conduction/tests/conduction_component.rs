//! The transient finite-element body on the composition kernel.
//!
//! `fs_conduction::component::ConductionBody` puts the backward-Euler owner
//! behind heat ports on named boundary patches. Here it passes the
//! conformance kit, a heater through the kernel steps it exactly as the
//! owner steps itself, and a lump joined to it through a conductance resolved
//! from a contact interface card meets the step condensed from the owner's
//! own responses and reaches the two bodies' shared equilibrium. A body whose
//! conductivity depends on temperature passes the same kit and heater on the
//! owner's Newton endpoint. Run with
//! `cargo test -p fs-conduction --features component --test conduction_component`.
#![cfg(feature = "component")]
// Exact float comparisons here assert bitwise agreement with the owner.
#![allow(clippy::float_cmp)]

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_conduction::component::{ConductionBody, PatchPort};
use fs_conduction::fixtures::{box_grid, on_box_face};
use fs_conduction::transient::VolumetricHeatCapacity;
use fs_conduction::transient::backward_euler::{BackwardEuler, NonlinearStepConfig, StepConfig};
use fs_conduction::{
    AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS, AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY,
    ConductionError, ConductionMesh, ConductionProblem, ConductivityModel, ConductivityTable,
    InterfaceResistance, LinearConfig, ResistanceOrigin, ScalarField, ThermalBc, ThermalBoundary,
    ThermalBoundaryBuilder,
};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, SurfaceSpec,
    SystemContext, UncertaintyModel,
};
use fs_ports::component::conformance::{Case, Observed, check_component, check_port_tangent};
use fs_ports::component::{Component, Proposal};
use fs_system::SystemBuilder;
use fs_system::elements::{Conductance, HeatCapacity, Heater, SourceAccount, Waveform};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let key = StreamKey {
            seed: 65,
            kernel_id: 1109,
            tile: 0,
            iteration: 0,
        };
        f(&Cx::new(
            &gate,
            arena,
            key,
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

/// A 30 x 20 x 20 mm block of 72 tetrahedra: steel-like `k` and `rho c`.
const LENGTH_M: f64 = 0.03;
const CONDUCTIVITY_W_M_K: f64 = 15.0;
const RHO_C_J_M3_K: f64 = 3.6e6;

/// The owner's step configuration. The Krylov tolerance is tight enough that
/// the solve's noise stays below the kernel's interface tolerance.
fn config() -> StepConfig {
    StepConfig {
        linear: LinearConfig {
            tolerance: 1.0e-13,
            max_iterations: 2000,
            restart: 20,
        },
        energy_tolerance_j: 1.0e-8,
    }
}

/// The owner's inputs, leaked: a body in a System lives as long as it does.
#[derive(Clone, Copy)]
struct Parts {
    mesh: &'static ConductionMesh,
    boundary: &'static ThermalBoundary,
    material: &'static ConductivityModel,
    source: &'static ScalarField,
}

fn leak<T>(value: T) -> &'static T {
    Box::leak(Box::new(value))
}

fn mesh() -> &'static ConductionMesh {
    let (complex, positions) = box_grid([3, 2, 2], [LENGTH_M, 0.02, 0.02]);
    leak(ConductionMesh::new(complex, positions).unwrap())
}

/// Two insulated end faces, `left` (x = 0) and `right` (x = L), and an
/// insulated remainder; each end has its own law.
fn boundary(mesh: &ConductionMesh, left: ThermalBc, right: ThermalBc) -> &'static ThermalBoundary {
    leak(
        ThermalBoundaryBuilder::new(mesh)
            .region("left", |face| on_box_face(face.centroid[0], 0.0), left)
            .unwrap()
            .region(
                "right",
                |face| on_box_face(face.centroid[0], LENGTH_M),
                right,
            )
            .unwrap()
            .adiabatic_remainder()
            .finish()
            .unwrap(),
    )
}

fn parts() -> Parts {
    let mesh = mesh();
    Parts {
        mesh,
        boundary: boundary(mesh, ThermalBc::adiabatic(), ThermalBc::adiabatic()),
        material: leak(ConductivityModel::isotropic_declared(CONDUCTIVITY_W_M_K).unwrap()),
        source: leak(ScalarField::Uniform(0.0)),
    }
}

fn problem(parts: Parts) -> ConductionProblem<'static> {
    ConductionProblem {
        mesh: parts.mesh,
        boundary: parts.boundary,
        material: parts.material,
        element_materials: None,
        source: parts.source,
    }
}

fn solver(cx: &Cx<'_>, mesh: &'static ConductionMesh) -> BackwardEuler<'static> {
    BackwardEuler::uniform(
        cx,
        mesh,
        VolumetricHeatCapacity::declared(RHO_C_J_M3_K).unwrap(),
    )
    .unwrap()
}

fn try_body(
    cx: &Cx<'_>,
    parts: Parts,
    ports: &[&str],
) -> Result<ConductionBody<'static>, ConductionError> {
    let ports: Vec<PatchPort<'_>> = ports.iter().map(|&name| PatchPort::Mean(name)).collect();
    ConductionBody::new(solver(cx, parts.mesh), problem(parts), &ports, config())
}

fn body(cx: &Cx<'_>, parts: Parts, ports: &[&str]) -> ConductionBody<'static> {
    try_body(cx, parts, ports).unwrap()
}

/// 300 K at the left end rising to 360 K at the right.
fn initial(parts: Parts) -> Vec<f64> {
    parts
        .mesh
        .positions()
        .iter()
        .map(|p| 300.0 + 2000.0 * p[0])
        .collect()
}

/// One step of the body as a Component, returning the end field and the
/// ports' temperatures.
fn step(
    body: &ConductionBody<'_>,
    old: &Vec<f64>,
    dt_s: f64,
    heat_w: &[f64],
    cx: &Cx<'_>,
) -> (Vec<f64>, Vec<f64>) {
    let mut new = Vec::new();
    let mut means = vec![0.0; heat_w.len()];
    let proposal = body
        .propose(old, dt_s, heat_w, &mut new, &mut means, cx)
        .unwrap();
    assert!(matches!(proposal, Proposal::Step(_)));
    (new, means)
}

fn largest_gap(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0, f64::max)
}

#[test]
fn the_body_conforms_and_its_tangent_is_the_derivative_of_its_step() {
    with_cx(|cx| {
        let parts = parts();
        let body = body(cx, parts, &["left", "right"]);
        let old = initial(parts);
        let stale = vec![1000.0; old.len()];
        let capacity_j_k = body.stored_energy_j(&vec![1.0; old.len()]).unwrap();
        // The declared rho c over the block's volume, to rounding.
        let volume_m3 = LENGTH_M * 0.02 * 0.02;
        assert!((capacity_j_k - RHO_C_J_M3_K * volume_m3).abs() <= 1.0e-12 * capacity_j_k);
        for (dt_s, heat_w) in [
            (0.5, [12.0, -3.0]),
            (2.0, [0.0, 0.0]),
            (10.0, [-40.0, 25.0]),
            (1.0e-3, [5.0, 5.0]),
            (4.0, [30.0, -30.0]),
        ] {
            let case = Case {
                state: &old,
                stale: &stale,
                dt_s,
                inputs: &heat_w,
            };
            let observed = check_component(&body, &case, cx)
                .unwrap_or_else(|failure| panic!("{dt_s} s, {heat_w:?} W: {failure}"));
            let Observed::Stepped {
                receipt,
                audit,
                outputs,
            } = observed
            else {
                panic!("{dt_s} s, {heat_w:?} W: the body asked for a shorter step");
            };
            // All the heat that came in is stored: backward Euler with
            // Neumann data conserves 1ᵀ C T to its solve and rounding.
            let heat_j = dt_s * (heat_w[0] + heat_w[1]);
            assert!(
                receipt.numerical_loss_j.abs() <= 1.0e-9 * (1.0 + heat_j.abs()),
                "{dt_s} s, {heat_w:?} W: {receipt:?}"
            );
            assert!((audit.port_energy_in_j - heat_j).abs() <= 1.0e-12 * (1.0 + heat_j.abs()));
            // Insulated, the field relaxes: the cold end warms, the hot end cools.
            if heat_w == [0.0, 0.0] {
                assert!(outputs[0] > 300.0 && outputs[1] < 360.0, "{outputs:?}");
            }
            // The kit differences at a step of cbrt(eps) of the largest input,
            // and the patch temperatures are rounded at eps of 400 K: below
            // that the difference cannot see the slope (at 1 ms it is 1e-4 K/W).
            let difference_step_w =
                heat_w.iter().fold(1.0_f64, |m, q| m.max(q.abs())) * f64::EPSILON.cbrt();
            let floor_k_w = 8.0 * f64::EPSILON * 400.0 / difference_step_w;
            check_port_tangent(&body, &case, 1.0e-6, floor_k_w, cx)
                .unwrap_or_else(|failure| panic!("tangent at {dt_s} s, {heat_w:?} W: {failure}"));
        }
        // Equal heat in and out stores nothing: the field only relaxes.
        let (end, _) = step(&body, &old, 4.0, &[30.0, -30.0], cx);
        let stored = body.stored_energy_j(&end).unwrap() - body.stored_energy_j(&old).unwrap();
        assert!(
            stored.abs() <= 1.0e-9 * body.stored_energy_j(&old).unwrap(),
            "{stored} J"
        );
    });
}

/// A trace port hands the body one load per vertex of its patch. Given the
/// consistent loads of a uniform flux, `A_f / 3` of each face's share on each
/// of its vertices, it is the mean port's step: the claim, declared before
/// the first run, is that the two end fields agree within 1e-9 K.
#[test]
fn a_trace_port_conforms_and_steps_the_body_as_the_mean_port_under_a_uniform_flux() {
    with_cx(|cx| {
        let parts = parts();
        let ports = [PatchPort::Trace("left"), PatchPort::Mean("right")];
        let traced =
            ConductionBody::new(solver(cx, parts.mesh), problem(parts), &ports, config()).unwrap();
        let vertices = traced.trace_vertices(0).unwrap().to_vec();
        assert!(traced.trace_vertices(1).is_none());
        // The left end of the 3 x 2 x 2 grid: nine vertices, in ascending order, all at x = 0.
        assert_eq!(vertices.len(), 9);
        assert!(vertices.windows(2).all(|pair| pair[0] < pair[1]));
        assert!(
            vertices
                .iter()
                .all(|&v| parts.mesh.positions()[v][0] == 0.0)
        );
        let old = initial(parts);
        let stale = vec![1000.0; old.len()];
        // The consistent loads of one watt spread uniformly over the left end.
        let area_m2 = traced.patch_area_m2(0);
        let mut uniform = [0.0; 9];
        for face in parts.mesh.boundary() {
            if on_box_face(face.centroid[0], 0.0) {
                for vertex in face.vertices {
                    let point = vertices.binary_search(&(vertex as usize)).unwrap();
                    uniform[point] += face.area / 3.0 / area_m2;
                }
            }
        }
        assert!((uniform.iter().sum::<f64>() - 1.0).abs() <= 8.0 * f64::EPSILON);
        let meaned = body(cx, parts, &["left", "right"]);
        for (dt_s, left_w, right_w) in [(0.5, 12.0, -3.0), (4.0, -40.0, 25.0)] {
            let mut inputs: Vec<f64> = uniform.iter().map(|share| share * left_w).collect();
            inputs.push(right_w);
            let (end, readings) = step(&traced, &old, dt_s, &inputs, cx);
            let (expected, means) = step(&meaned, &old, dt_s, &[left_w, right_w], cx);
            let gap_k = largest_gap(&end, &expected);
            eprintln!("{dt_s} s, {left_w} W: trace against mean port, largest gap {gap_k:e} K");
            assert!(gap_k <= 1.0e-9, "{dt_s} s: {gap_k} K");
            // A trace port returns the end temperature at each of its points.
            for (reading, &vertex) in readings.iter().zip(&vertices) {
                assert_eq!(reading.to_bits(), end[vertex].to_bits());
            }
            assert!((readings[9] - means[1]).abs() <= 1.0e-9);
            // Loads that are not a uniform flux: the kit's laws and the tangent.
            for (point, load) in inputs[..9].iter_mut().enumerate() {
                *load += 0.3 * left_w * (point as f64 - 4.0) / 9.0;
            }
            let case = Case {
                state: &old,
                stale: &stale,
                dt_s,
                inputs: &inputs,
            };
            let observed = check_component(&traced, &case, cx)
                .unwrap_or_else(|failure| panic!("{dt_s} s: {failure}"));
            let Observed::Stepped { receipt, audit, .. } = observed else {
                panic!("{dt_s} s: the body asked for a shorter step");
            };
            let heat_j = dt_s * inputs.iter().sum::<f64>();
            assert!(receipt.numerical_loss_j.abs() <= 1.0e-9 * (1.0 + heat_j.abs()));
            assert!((audit.port_energy_in_j - heat_j).abs() <= 1.0e-12 * (1.0 + heat_j.abs()));
            let difference_step_w =
                inputs.iter().fold(1.0_f64, |m, q| m.max(q.abs())) * f64::EPSILON.cbrt();
            let floor_k_w = 8.0 * f64::EPSILON * 400.0 / difference_step_w;
            check_port_tangent(&traced, &case, 1.0e-6, floor_k_w, cx)
                .unwrap_or_else(|failure| panic!("tangent at {dt_s} s: {failure}"));
        }
    });
}

#[test]
fn a_heater_through_the_kernel_steps_the_body_as_the_owner_does() {
    const POWER_W: f64 = 20.0;
    const DT_S: f64 = 0.5;
    with_cx(|cx| {
        let parts = parts();
        let mut builder = SystemBuilder::new();
        let node = builder
            .add_with_tangent("body", body(cx, parts, &["left"]), initial(parts))
            .unwrap();
        let heater = builder
            .add_with_tangent(
                "heater",
                Heater::new(Waveform::constant(POWER_W).unwrap()).unwrap(),
                SourceAccount::default(),
            )
            .unwrap();
        let heat = builder.port(heater, "heat").unwrap();
        builder
            .bond(builder.port(node, "left").unwrap(), heat)
            .unwrap();
        let mut system = builder.seal().unwrap();
        // The bond is a loop, so the heater's temperature input needs a
        // starting value.
        assert!(system.guess(heat, &[300.0]));

        // The owner driven directly: its left end prescribes the flux -P / A.
        let area_m2 = body(cx, parts, &["left"]).patch_area_m2(0);
        let heated = Parts {
            boundary: boundary(
                parts.mesh,
                ThermalBc::neumann(-POWER_W / area_m2).unwrap(),
                ThermalBc::adiabatic(),
            ),
            ..parts
        };
        let owner = solver(cx, parts.mesh);
        let mut field = initial(parts);
        let start_j = body(cx, parts, &["left"]).stored_energy_j(&field).unwrap();
        for index in 0..8 {
            system.step(DT_S, cx).unwrap();
            field = owner
                .advance(cx, problem(heated), None, &field, DT_S, config())
                .unwrap()
                .temperature;
            let kernel = system.state::<Vec<f64>>(node).unwrap();
            assert!(
                kernel
                    .iter()
                    .zip(&field)
                    .all(|(a, b)| a.to_bits() == b.to_bits()),
                "step {index}: the kernel's body differs from the owner by {} K",
                largest_gap(kernel, &field)
            );
        }
        // The heater's 80 J are in the body and the window closes.
        let stored_j = body(cx, parts, &["left"]).stored_energy_j(&field).unwrap() - start_j;
        assert!(
            (stored_j - 8.0 * DT_S * POWER_W).abs() <= 1.0e-9 * start_j,
            "{stored_j} J"
        );
        let ledger = system.ledger();
        assert!(
            ledger.energy_residual_j().abs() <= 1.0e-9 * start_j,
            "{}",
            ledger.energy_residual_j()
        );
        // The heater dissipates its power at the body's patch temperature.
        assert!(ledger.entropy_generated_j_per_k() > 0.0);
    });
}

/// The declared contact between the lump and the body's left patch: one
/// area-specific resistance claim on an ordered interface card.
const CONTACT_M2_K_W: f64 = 8.0e-4;

fn contact_card() -> InterfaceSystemCard {
    let surface = |chemistry: &str, texture_frame: &str| SurfaceSpec {
        material: MaterialStateId {
            chemistry: chemistry.to_string(),
            phase: "solid".to_string(),
            process: "as-fixtured".to_string(),
            revision: 0,
        },
        texture_frame: texture_frame.to_string(),
    };
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(
                AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY,
                AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS,
            ),
            value: PropertyValue::Scalar {
                value: CONTACT_M2_K_W,
                dims: AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS,
            },
            validity: ValidityDomain::unconstrained(),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: Vec::new(),
            provenance: Provenance {
                source: "conduction_component lump-to-body declared contact fixture".to_string(),
                license: "internal-test-use".to_string(),
                artifact: None,
            },
        })
        .unwrap();
    InterfaceSystemCard::assemble(
        surface("lump", "interface-normal-plus-x"),
        surface("body", "interface-normal-minus-x"),
        SystemContext {
            medium: "dry".to_string(),
            third_body: Some("declared-contact-layer".to_string()),
            environment: "vacuum".to_string(),
            history: "unaged".to_string(),
        },
        claims,
        Vec::new(),
    )
    .unwrap()
}

/// The wall's conductance from the card: the resolved claim, its receipt
/// replayed against the card, one series term `R'' / A` over the patch whose
/// origin cites the card and that receipt.
fn card_wall_w_k(area_m2: f64) -> f64 {
    let card = contact_card();
    let contact = InterfaceResistance::from_card(
        "lump-body",
        &card,
        &QueryPoint::new(),
        SelectionPolicy::SingleClaimOnly,
    )
    .unwrap();
    card.claims().verify_receipt(contact.receipt()).unwrap();
    assert_eq!(contact.value_m2_k_per_w(), CONTACT_M2_K_W);
    let term = contact.term_for_area("lump-body", area_m2).unwrap();
    let ResistanceOrigin::InterfaceCard {
        card_identity,
        receipt,
    } = term.origin()
    else {
        panic!("the term cites {:?}, not the card", term.origin());
    };
    assert_eq!(*card_identity, card.content_hash());
    assert_eq!(receipt, contact.receipt());
    1.0 / term.value_k_per_w()
}

#[test]
fn a_lump_through_a_card_contact_meets_the_condensed_step_and_the_shared_equilibrium() {
    const LUMP_J_K: f64 = 20.0;
    const LUMP_START_K: f64 = 420.0;
    const DT_S: f64 = 2.0;
    const REFERENCE_W: f64 = 10.0;
    // 3000 s: some fifty times the slowest joint mode (about 60 s).
    const STEPS: usize = 1500;
    with_cx(|cx| {
        let parts = parts();
        // The owner's own responses, read through an unconnected copy.
        let probe = body(cx, parts, &["left"]);
        // 8e-4 m^2 K/W over the 20 x 20 mm patch: 0.5 W/K.
        let wall_w_k = card_wall_w_k(probe.patch_area_m2(0));
        assert!((wall_w_k - 0.5).abs() <= 1.0e-12, "{wall_w_k} W/K");
        let mut builder = SystemBuilder::new();
        let node = builder
            .add_with_tangent("body", body(cx, parts, &["left"]), initial(parts))
            .unwrap();
        let lump = builder
            .add_with_tangent("lump", HeatCapacity::new(LUMP_J_K).unwrap(), LUMP_START_K)
            .unwrap();
        let wall = builder
            .add_with_tangent("wall", Conductance::new(wall_w_k).unwrap(), ())
            .unwrap();
        let (a, b) = (
            builder.port(wall, "a").unwrap(),
            builder.port(wall, "b").unwrap(),
        );
        builder
            .bond(builder.port(lump, "heat").unwrap(), a)
            .unwrap();
        builder
            .bond(b, builder.port(node, "left").unwrap())
            .unwrap();
        let mut system = builder.seal().unwrap();
        assert!(system.guess(a, &[LUMP_START_K]));
        assert!(system.guess(b, &[300.0]));

        let body_j_k = probe
            .stored_energy_j(&vec![1.0; parts.mesh.vertex_count()])
            .unwrap();
        let start = initial(parts);
        let start_j = LUMP_J_K * LUMP_START_K + probe.stored_energy_j(&start).unwrap();
        let equilibrium_k = start_j / (LUMP_J_K + body_j_k);
        let mut worst_k = 0.0_f64;
        for index in 0..STEPS {
            let lump_k = *system.state::<f64>(lump).unwrap();
            let field = system.state::<Vec<f64>>(node).unwrap().clone();
            // The body's patch temperature is affine in its heat: T = a + r Q
            // from two owner steps. The lump's is its mean over the step,
            // T_L - dt Q / (2 C_L), so the wall's heat Q = G (T_L - T) is
            // Q = G (T_L - a) / (1 + G (dt / (2 C_L) + r)).
            let (_, at_zero) = step(&probe, &field, DT_S, &[0.0], cx);
            let (_, at_reference) = step(&probe, &field, DT_S, &[REFERENCE_W], cx);
            let slope_k_w = (at_reference[0] - at_zero[0]) / REFERENCE_W;
            let heat_w = wall_w_k * (lump_k - at_zero[0])
                / (1.0 + wall_w_k * (0.5 * DT_S / LUMP_J_K + slope_k_w));
            let (expected, _) = step(&probe, &field, DT_S, &[heat_w], cx);

            system.step(DT_S, cx).unwrap();
            let lump_gap =
                (*system.state::<f64>(lump).unwrap() - (lump_k - DT_S * heat_w / LUMP_J_K)).abs();
            let body_gap = largest_gap(system.state::<Vec<f64>>(node).unwrap(), &expected);
            worst_k = worst_k.max(lump_gap).max(body_gap);
            assert!(
                lump_gap <= 1.0e-10 && body_gap <= 1.0e-10,
                "step {index}: {lump_gap} K, {body_gap} K"
            );
        }
        // Closed: the two bodies hold the energy they started with, and both
        // reach the temperature that shares it.
        let lump_k = *system.state::<f64>(lump).unwrap();
        let field = system.state::<Vec<f64>>(node).unwrap();
        let end_j = LUMP_J_K * lump_k + probe.stored_energy_j(field).unwrap();
        assert!(
            (end_j - start_j).abs() <= 1.0e-9 * start_j,
            "{} J",
            end_j - start_j
        );
        assert!(
            (lump_k - equilibrium_k).abs() <= 1.0e-10,
            "lump {lump_k} K, shared {equilibrium_k} K"
        );
        let spread = field
            .iter()
            .map(|t| (t - equilibrium_k).abs())
            .fold(0.0, f64::max);
        assert!(spread <= 1.0e-10, "body {spread} K from {equilibrium_k} K");
        assert!(system.ledger().energy_residual_j().abs() <= 1.0e-9 * start_j);
        eprintln!(
            "condensed step: largest gap {worst_k:e} K over {STEPS} steps; lump {lump_k} K, body within {spread:e} K of the shared {equilibrium_k} K"
        );
    });
}

#[test]
fn the_body_refuses_what_it_cannot_carry_through_its_ports() {
    with_cx(|cx| {
        let parts = parts();
        let refused = |parts: Parts, ports: &[&str], what: &str| match try_body(cx, parts, ports) {
            Err(ConductionError::Config {
                parameter: "conduction body",
                what: why,
            }) => {
                assert!(why.contains(what), "{why:?} does not say {what:?}");
            }
            Err(other) => panic!("{what}: {other:?}"),
            Ok(_) => panic!("{what}: accepted"),
        };
        let with_right = |right: ThermalBc| Parts {
            boundary: boundary(parts.mesh, ThermalBc::adiabatic(), right),
            ..parts
        };
        refused(
            with_right(ThermalBc::dirichlet(300.0).unwrap()),
            &["left"],
            "region right is a dirichlet row",
        );
        refused(
            with_right(ThermalBc::robin(40.0, 300.0).unwrap()),
            &["left"],
            "region right is a robin row",
        );
        refused(
            with_right(ThermalBc::neumann(5.0).unwrap()),
            &["left"],
            "region right is a neumann row",
        );
        // A port's own region may not carry a flux of its own either.
        refused(
            with_right(ThermalBc::neumann(5.0).unwrap()),
            &["left", "right"],
            "region right",
        );
        refused(
            Parts {
                source: leak(ScalarField::Uniform(1.0e4)),
                ..parts
            },
            &["left"],
            "volumetric source",
        );
        let curve = ConductivityTable::declared_curve(vec![(250.0, 10.0), (450.0, 20.0)]).unwrap();
        refused(
            Parts {
                material: leak(ConductivityModel::isotropic(curve)),
                ..parts
            },
            &["left"],
            "conductivity",
        );
        refused(parts, &["top"], "port top names no boundary region");
        refused(parts, &["left", "left"], "port left is named twice");
        refused(parts, &[], "at least one port");
        // A solver prepared on another mesh of the same shape.
        let other = solver(cx, mesh());
        match ConductionBody::new(other, problem(parts), &[PatchPort::Mean("left")], config()) {
            Err(ConductionError::Config { what, .. }) => {
                assert!(what.contains("same mesh"), "{what}");
            }
            other => panic!("another mesh: {other:?}"),
        }
        // A step the owner refuses comes back as the owner's error.
        let body = body(cx, parts, &["left"]);
        let mut new = Vec::new();
        let mut out = [0.0];
        assert!(
            body.propose(&vec![300.0; 3], 1.0, &[1.0], &mut new, &mut out, cx)
                .is_err()
        );
        assert!(
            body.propose(&initial(parts), 1.0, &[f64::NAN], &mut new, &mut out, cx)
                .is_err()
        );
    });
}

/// A conductivity linear in temperature, 10 W/m/K at 250 K to 20 W/m/K at
/// 450 K, with no knot in between: the tangent sees no kink.
fn dependent_parts() -> Parts {
    let curve = ConductivityTable::declared_curve(vec![(250.0, 10.0), (450.0, 20.0)]).unwrap();
    Parts {
        material: leak(ConductivityModel::isotropic(curve)),
        ..parts()
    }
}

/// The owner's Newton policy, tighter than its default so that where Newton
/// stops lies far below what the kit's difference quotient can see.
fn newton() -> NonlinearStepConfig {
    NonlinearStepConfig {
        residual_rtol: 1.0e-12,
        residual_atol_j: 1.0e-12,
        ..NonlinearStepConfig::default()
    }
}

fn try_dependent(
    cx: &Cx<'_>,
    parts: Parts,
    ports: &[&str],
    config: StepConfig,
    policy: NonlinearStepConfig,
) -> Result<ConductionBody<'static>, ConductionError> {
    let ports: Vec<PatchPort<'_>> = ports.iter().map(|&name| PatchPort::Mean(name)).collect();
    ConductionBody::temperature_dependent(
        solver(cx, parts.mesh),
        problem(parts),
        &ports,
        config,
        policy,
    )
}

/// The owner's step configuration with its Krylov solves taken to 1e-15.
/// A k(T) body's step Jacobian is not symmetric, so the tangent's forward
/// and adjoint solves agree only as far as each converges, not to rounding
/// as for a constant conductivity, where both apply one symmetric
/// polynomial of the same matrix. The kit's transpose law allows rounding
/// only: at 1e-13 the first run broke it at 4 s by 5e-13 relative.
fn converged() -> StepConfig {
    let mut config = config();
    config.linear.tolerance = 1.0e-15;
    config
}

/// Claims, declared before the first run: the k(T) body passes the
/// conformance kit, and its port tangent is the difference quotient of its
/// step within 1e-6 relative, above a floor of the patch rounding plus
/// 1e-9 K for where Newton stops. The kit runs on `converged`.
#[test]
fn a_temperature_dependent_body_conforms_and_its_tangent_differentiates_its_newton_step() {
    with_cx(|cx| {
        let parts = dependent_parts();
        let body = try_dependent(cx, parts, &["left", "right"], converged(), newton()).unwrap();
        let old = initial(parts);
        let stale = vec![1000.0; old.len()];
        for (dt_s, heat_w) in [
            (0.5, [12.0, -3.0]),
            (2.0, [0.0, 0.0]),
            (10.0, [-40.0, 25.0]),
            (4.0, [30.0, -30.0]),
        ] {
            let case = Case {
                state: &old,
                stale: &stale,
                dt_s,
                inputs: &heat_w,
            };
            let observed = check_component(&body, &case, cx)
                .unwrap_or_else(|failure| panic!("{dt_s} s, {heat_w:?} W: {failure}"));
            let Observed::Stepped { receipt, audit, .. } = observed else {
                panic!("{dt_s} s, {heat_w:?} W: the body asked for a shorter step");
            };
            let heat_j = dt_s * (heat_w[0] + heat_w[1]);
            assert!(
                receipt.numerical_loss_j.abs() <= 1.0e-9 * (1.0 + heat_j.abs()),
                "{dt_s} s, {heat_w:?} W: {receipt:?}"
            );
            assert!((audit.port_energy_in_j - heat_j).abs() <= 1.0e-12 * (1.0 + heat_j.abs()));
            let difference_step_w =
                heat_w.iter().fold(1.0_f64, |m, q| m.max(q.abs())) * f64::EPSILON.cbrt();
            let floor_k_w = (8.0 * f64::EPSILON * 400.0 + 1.0e-9) / difference_step_w;
            check_port_tangent(&body, &case, 1.0e-6, floor_k_w, cx)
                .unwrap_or_else(|failure| panic!("tangent at {dt_s} s, {heat_w:?} W: {failure}"));
        }
    });
}

/// Claim, declared before the first run: a heater through the kernel steps
/// the k(T) body bit for bit as the owner's `advance_nonlinear` steps
/// itself, and the window closes.
#[test]
fn a_heater_through_the_kernel_steps_a_temperature_dependent_body_as_the_owner_by_newton() {
    const POWER_W: f64 = 20.0;
    const DT_S: f64 = 0.5;
    with_cx(|cx| {
        let parts = dependent_parts();
        let body = try_dependent(cx, parts, &["left"], config(), newton()).unwrap();
        let mut builder = SystemBuilder::new();
        let node = builder
            .add_with_tangent(
                "body",
                try_dependent(cx, parts, &["left"], config(), newton()).unwrap(),
                initial(parts),
            )
            .unwrap();
        let heater = builder
            .add_with_tangent(
                "heater",
                Heater::new(Waveform::constant(POWER_W).unwrap()).unwrap(),
                SourceAccount::default(),
            )
            .unwrap();
        let heat = builder.port(heater, "heat").unwrap();
        builder
            .bond(builder.port(node, "left").unwrap(), heat)
            .unwrap();
        let mut system = builder.seal().unwrap();
        assert!(system.guess(heat, &[300.0]));
        let area_m2 = body.patch_area_m2(0);
        let heated = Parts {
            boundary: boundary(
                parts.mesh,
                ThermalBc::neumann(-POWER_W / area_m2).unwrap(),
                ThermalBc::adiabatic(),
            ),
            ..parts
        };
        let owner = solver(cx, parts.mesh);
        let mut field = initial(parts);
        let start_j = body.stored_energy_j(&field).unwrap();
        let mut newton_steps = 0;
        for index in 0..8 {
            system.step(DT_S, cx).unwrap();
            let solved = owner
                .advance_nonlinear(cx, problem(heated), None, &field, DT_S, config(), newton())
                .unwrap();
            newton_steps += solved.nonlinear_iterations;
            field = solved.step.temperature;
            let kernel = system.state::<Vec<f64>>(node).unwrap();
            assert!(
                kernel
                    .iter()
                    .zip(&field)
                    .all(|(a, b)| a.to_bits() == b.to_bits()),
                "step {index}: the kernel's body differs from the owner by {} K",
                largest_gap(kernel, &field)
            );
        }
        let stored_j = body.stored_energy_j(&field).unwrap() - start_j;
        eprintln!(
            "k(T) body, 8 steps of {DT_S} s at {POWER_W} W: {newton_steps} Newton updates, stored {stored_j} J"
        );
        assert!(newton_steps >= 8, "Newton did no work: {newton_steps}");
        assert!(
            (stored_j - 8.0 * DT_S * POWER_W).abs() <= 1.0e-9 * start_j,
            "{stored_j} J"
        );
        let ledger = system.ledger();
        assert!(
            ledger.energy_residual_j().abs() <= 1.0e-9 * start_j,
            "{}",
            ledger.energy_residual_j()
        );
    });
}

#[test]
fn a_body_is_put_behind_ports_on_the_path_its_conductivity_needs() {
    with_cx(|cx| {
        let says =
            |result: Result<ConductionBody<'static>, ConductionError>, what: &str| match result {
                Err(error) => {
                    let why = error.to_string();
                    assert!(why.contains(what), "{why:?} does not say {what:?}");
                }
                Ok(_) => panic!("{what}: accepted"),
            };
        says(
            try_body(cx, dependent_parts(), &["left"]),
            "ConductionBody::temperature_dependent",
        );
        says(
            try_dependent(cx, parts(), &["left"], config(), newton()),
            "ConductionBody::new",
        );
        says(
            try_dependent(
                cx,
                dependent_parts(),
                &["left"],
                config(),
                NonlinearStepConfig {
                    max_iterations: 0,
                    ..newton()
                },
            ),
            "positive Newton budget",
        );
        // An endpoint beyond the curve's 450 K is the owner's refusal.
        let body = try_dependent(cx, dependent_parts(), &["left"], config(), newton()).unwrap();
        let old = vec![440.0; dependent_parts().mesh.vertex_count()];
        let mut new = Vec::new();
        let mut out = [0.0];
        let refused = body.propose(&old, 10.0, &[2000.0], &mut new, &mut out, cx);
        eprintln!("440 K and 2000 W for 10 s: {refused:?}");
        assert!(refused.is_err());
        assert!(new.is_empty());
    });
}
