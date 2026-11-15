//! `conduction_body`: fs-conduction's backward-Euler finite-element body,
//! lowered from a scene.
//!
//! The check is the one-dimensional closed form of fs-conduction
//! tests/conduction_pair.rs, from a scene. Two equal slabs,
//! examples/system/two_slabs.json, start at 400 K and 300 K with insulated
//! outer ends and exchange heat through a 0.5 W/K conductance (a contact of
//! 8e-4 m^2 K/W over the 20 x 20 mm end, Bi = 5). The claim, declared before
//! the first run: the energy that has crossed the wall by t* = 0.1 L^2 /
//! alpha converges to the closed form at order at least
//! [`CLAIMED_ORDER`] between 4 and 8 cells along each slab, with the step
//! t* / (2 n^2); the pair's stored energy and the ledger close; the wall's
//! last rate is its conductance times the two faces' reported temperatures.
//! Then the refusals.
//!
//! `trace_contact`: the same two slabs joined trace to trace,
//! examples/system/two_slabs_trace.json (equal cross grids, the owner's
//! matching contact) and two_slabs_mapped.json (2 x 2 against 3 x 3, the
//! owner's common refinement), through 8e-4 m^2 K/W declared in the scene.
//! The reference is fs-conduction's own contact between the slabs on one
//! mesh ([`owner_contact`]), stepped with the scene's step and budgets. The
//! claims, declared before the first run, on 4 and 8 cells along each slab:
//! the energy that crossed the contact within [`CLAIMED_GAP_J`]; the least
//! and greatest temperature of each slab and of each contact face within
//! [`CLAIMED_GAP_K`]; the pair's stored energy and the ledger closed to 1e-9
//! of what the pair holds. Then the examples through the binary and the
//! refusals.

use std::f64::consts::{FRAC_PI_2, PI};

use super::*;
use fs_cli::SystemScene;
use fs_conduction::fixtures::{box_grid, on_box_face};
use fs_conduction::interface::{NonmatchingOptions, NonmatchingSurface};
use fs_conduction::transient::VolumetricHeatCapacity;
use fs_conduction::transient::backward_euler::{BackwardEuler, StepConfig};
use fs_conduction::{
    AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS, AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY, BoundaryFace,
    ConductionMesh, ConductionProblem, ConductivityModel, InterfaceFacePair, InterfaceResistance,
    InterfaceSurface, LinearConfig, ScalarField, ThermalBoundary, ThermalBoundaryBuilder,
    ThermalInterfaces,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, SurfaceSpec,
    SystemContext, UncertaintyModel,
};
use fs_rep_mesh::TetComplex;

const LENGTH_M: f64 = 0.03;
const CROSS_M: f64 = 0.02;
const RHO_C_J_M3_K: f64 = 3.6e6;
const CONDUCTANCE_W_K: f64 = 0.5;
/// Bi = 2 L / (k R''), R'' = A / G, for the half-problem about the wall.
const BIOT: f64 = 2.0 * LENGTH_M * CONDUCTANCE_W_K / (15.0 * CROSS_M * CROSS_M);
const HOT_K: f64 = 400.0;
const COLD_K: f64 = 300.0;
/// The comparison time in units of L^2 / alpha.
const FOURIER: f64 = 0.1;
/// The least order between the two rungs that meets the claim.
const CLAIMED_ORDER: f64 = 1.8;

/// The excess heat the hot slab holds over the mean at the start [J].
fn excess_j() -> f64 {
    RHO_C_J_M3_K * LENGTH_M * CROSS_M * CROSS_M * 0.5 * (HOT_K - COLD_K)
}

/// The energy that has crossed the wall by t*, from the closed form: each
/// slab is the plane wall insulated at its outer end with a Robin end of
/// `h = 2 G / A`, still holding `sum C_n sin(b_n) / b_n exp(-b_n^2 Fo)` of its
/// excess, `b_n tan b_n = Bi`, `C_n = 4 sin b_n / (2 b_n + sin 2 b_n)`.
fn closed_form_crossed_j() -> f64 {
    let mut held = 0.0;
    for n in 0..64 {
        // b tan b - Bi rises from -Bi to +infinity on (n pi, n pi + pi/2).
        let start = f64::from(n) * PI;
        let (mut lo, mut hi) = (start, start + FRAC_PI_2);
        loop {
            let mid = 0.5 * (lo + hi);
            if mid <= lo || mid >= hi {
                break;
            }
            if mid * mid.tan() < BIOT {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        let b = 0.5 * (lo + hi);
        let term =
            4.0 * b.sin() / (2.0 * b + (2.0 * b).sin()) * b.sin() / b * (-b * b * FOURIER).exp();
        held += term;
        if term.abs() < 1.0e-17 {
            break;
        }
    }
    excess_j() * (1.0 - held)
}

fn probe(scene: &SystemScene, path: &str) -> f64 {
    scene
        .probe(path)
        .unwrap()
        .unwrap_or_else(|| panic!("{path} is not reported"))
}

/// The example with `cells` along each slab and its step t* / (2 n^2).
fn two_slabs(cells: u32) -> String {
    let text = std::fs::read_to_string(example("two_slabs.json")).unwrap();
    assert_eq!(text.matches("\"cells\": [4, 2, 2]").count(), 2);
    assert!(text.contains("\"step_s\": 0.675,") && text.contains("\"duration_s\": 21.6,"));
    let step_s = 21.6 / f64::from(2 * cells * cells);
    text.replace(
        "\"cells\": [4, 2, 2]",
        &format!("\"cells\": [{cells}, 2, 2]"),
    )
    .replace("\"step_s\": 0.675,", &format!("\"step_s\": {step_s},"))
}

#[test]
fn two_conduction_bodies_from_a_scene_converge_to_the_closed_form_at_second_order() {
    let exact_j = closed_form_crossed_j();
    // The pair's stored energy, rho c V (T_a + T_b).
    let held_j = RHO_C_J_M3_K * LENGTH_M * CROSS_M * CROSS_M * (HOT_K + COLD_K);
    eprintln!("closed form: {exact_j:.15e} J crossed (Bi {BIOT})");
    let mut errors = Vec::new();
    for cells in [4, 8] {
        let text = two_slabs(cells);
        scene_context(|cx, _| {
            let mut scene = SystemScene::from_text(&text, "two_slabs.json").unwrap();
            let mut crossed_j = 0.0;
            let mut time_s = scene.time_s();
            while scene.advance(f64::MAX, cx).unwrap() {
                let now_s = scene.time_s();
                crossed_j += probe(&scene, "wall.heat_rate_w") * (now_s - time_s);
                time_s = now_s;
            }
            let (face_a, face_b, rate) = (
                probe(&scene, "a.x_max.temperature_k"),
                probe(&scene, "b.x_min.temperature_k"),
                probe(&scene, "wall.heat_rate_w"),
            );
            let wall = rate / (CONDUCTANCE_W_K * (face_a - face_b)) - 1.0;
            let ledger = scene.ledger();
            let (stored, residual) = (ledger.stored_change_j(), ledger.energy_residual_j());
            eprintln!(
                "n {cells}: {} steps to {time_s} s; crossed {crossed_j:.15e} J; faces {face_a} K, {face_b} K, wall rate / G dT - 1 = {wall:e}; a {} .. {} K, b {} .. {} K; stored change {stored:e} J, ledger residual {residual:e} J",
                scene.steps(),
                probe(&scene, "a.min_temperature_k"),
                probe(&scene, "a.max_temperature_k"),
                probe(&scene, "b.min_temperature_k"),
                probe(&scene, "b.max_temperature_k"),
            );
            assert!(
                (time_s / 21.6 - 1.0).abs() <= 1e-12,
                "the run ends at t*: {time_s}"
            );
            // The wall reads the faces the bodies return on the settled load
            // line.
            assert!(wall.abs() <= 1e-9, "wall rate off its faces by {wall:e}");
            // The pair is closed: nothing enters, so what one slab gives the
            // other takes, up to the owner's linear solves.
            assert!(
                stored.abs() <= 1e-9 * held_j,
                "the pair changed by {stored:e} J"
            );
            assert!(
                residual.abs() <= 1e-9 * held_j,
                "ledger residual {residual:e} J"
            );
            errors.push((crossed_j - exact_j).abs() / excess_j());
        });
    }
    let order = (errors[0] / errors[1]).log2();
    eprintln!("relative errors {errors:?}; order {order:.3}");
    assert!(
        order >= CLAIMED_ORDER,
        "order {order:.3} below {CLAIMED_ORDER}"
    );
}

#[test]
fn the_two_slab_example_runs_through_the_binary_and_reports_its_faces() {
    let path = example("two_slabs.json");
    let (code, json, stderr) = run_binary(&["--json", "system", &path]);
    assert_eq!(code, 0, "{stderr}");
    let (face_a, face_b) = (
        value(&json, &element("a"), "x_max.temperature_k"),
        value(&json, &element("b"), "x_min.temperature_k"),
    );
    let rate = value(&json, &element("wall"), "heat_rate_w");
    assert!(
        COLD_K < face_b && face_b < face_a && face_a < HOT_K && rate > 0.0,
        "faces {face_a} K, {face_b} K, rate {rate} W"
    );
    let (code, text, stderr) = run_binary(&["system", &path]);
    assert_eq!(code, 0, "{stderr}");
    assert!(
        text.contains("conduction bodies step fs-conduction"),
        "{text}"
    );
}

#[test]
fn conduction_body_scenes_refuse_bad_boxes_ports_materials_and_solves_by_name() {
    let base = two_slabs(4);
    for (name, from, to, reason) in [
        (
            "face",
            "\"ports\": [\"x_max\"]",
            "\"ports\": [\"w_max\"]",
            "ports name box faces: x_min, x_max",
        ),
        (
            "twice",
            "\"ports\": [\"x_max\"]",
            "\"ports\": [\"x_max\", \"x_max\"]",
            "port x_max is named twice",
        ),
        (
            "no ports",
            "\"ports\": [\"x_max\"]",
            "\"ports\": []",
            "ports is a nonempty array of box faces",
        ),
        (
            "fraction",
            "\"cells\": [4, 2, 2]",
            "\"cells\": [4.5, 2, 2]",
            "cells are whole numbers from 1 to 64",
        ),
        (
            "too many",
            "\"cells\": [4, 2, 2]",
            "\"cells\": [64, 64, 2]",
            "at most 4096 cells",
        ),
        (
            "box",
            "\"box_m\": [0.03, 0.02, 0.02]",
            "\"box_m\": [0.03, 0.02]",
            "box_m is an array of three numbers",
        ),
        (
            "capacity",
            "\"volumetric_heat_capacity_j_m3_k\": 3600000.0",
            "\"volumetric_heat_capacity_j_m3_k\": -3600000.0",
            "is not positive",
        ),
        (
            "solve",
            ",\n      \"solve\": {\"energy_tolerance_j\": 1e-8, \"tolerance\": 1e-13, \"max_iterations\": 4000, \"restart\": 20}",
            "",
            "missing solve (energy_tolerance_j)",
        ),
        (
            "energy",
            "\"energy_tolerance_j\": 1e-8, ",
            "",
            "missing energy_tolerance_j",
        ),
        (
            "tolerance",
            "\"tolerance\": 1e-13",
            "\"tolerance\": 1.5",
            "tolerance must lie in (0, 1)",
        ),
        (
            "restart",
            "\"restart\": 20",
            "\"restart\": 0",
            "restart is a whole number from 1 to 1000",
        ),
        (
            "field",
            "\"temperature_k\": 400.0,",
            "\"temperature_k\": 400.0, \"emissivity\": 0.5,",
            "unsupported field \"emissivity\"",
        ),
    ] {
        assert!(base.contains(from), "{name}: test fixture changed");
        let (code, stderr) = refused(name, &base.replacen(from, to, 1));
        assert_eq!(code, 3, "{name}: {stderr}");
        assert!(stderr.contains(reason), "{name}: {stderr}");
    }
}

/// The contact of the two trace examples [m^2 K/W]: the 0.5 W/K wall of
/// two_slabs.json over its 20 x 20 mm face.
const CONTACT_M2_K_W: f64 = CROSS_M * CROSS_M / CONDUCTANCE_W_K;
/// How far a trace scene may lie from the owner's contact on one mesh.
const CLAIMED_GAP_J: f64 = 1.0e-9;
const CLAIMED_GAP_K: f64 = 1.0e-9;

/// What a trace scene and the owner's contact are compared in.
struct Crossed {
    /// The energy that crossed the contact from slab a to slab b [J].
    energy_j: f64,
    /// The least and the greatest end temperature of slab a, of slab b, of
    /// a's contact face and of b's [K].
    spans_k: [(f64, f64); 4],
}

fn span(values: impl Iterator<Item = f64>) -> (f64, f64) {
    values.fold((f64::INFINITY, f64::NEG_INFINITY), |(low, high), t| {
        (low.min(t), high.max(t))
    })
}

/// The examples' resistance as the owner's, through a declared card.
fn owner_resistance() -> InterfaceResistance {
    let side = |chemistry: &str| SurfaceSpec {
        material: MaterialStateId {
            chemistry: chemistry.to_string(),
            phase: "solid".to_string(),
            process: "as-fixtured".to_string(),
            revision: 0,
        },
        texture_frame: "interface-normal-x".to_string(),
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
                source: "system_command trace-contact reference fixture".to_string(),
                license: "internal-test-use".to_string(),
                artifact: None,
            },
        })
        .unwrap();
    let card = InterfaceSystemCard::assemble(
        side("slab-a"),
        side("slab-b"),
        SystemContext {
            medium: "dry".to_string(),
            third_body: None,
            environment: "unspecified".to_string(),
            history: "unaged".to_string(),
        },
        claims,
        Vec::new(),
    )
    .unwrap();
    InterfaceResistance::from_card(
        "slab-contact",
        &card,
        &QueryPoint::new(),
        SelectionPolicy::SingleClaimOnly,
    )
    .unwrap()
}

/// The owner's contact over the plane x = L of the two-slab mesh: the
/// matching contact of its coincident faces, or the non-matching one with
/// the refinement budgets of two_slabs_mapped.json.
fn owner_interfaces(
    cx: &fs_exec::Cx<'_>,
    mesh: &ConductionMesh,
    boundary: &ThermalBoundary,
    matching: bool,
) -> ThermalInterfaces {
    if matching {
        let pairs = ThermalInterfaces::coincident_face_pairs(mesh)
            .unwrap()
            .into_iter()
            .map(|pair| {
                if mesh.boundary()[pair.side_a].outward_normal[0] > 0.0 {
                    pair
                } else {
                    InterfaceFacePair {
                        side_a: pair.side_b,
                        side_b: pair.side_a,
                    }
                }
            })
            .collect();
        let surface = InterfaceSurface::new("slab-contact", pairs, owner_resistance()).unwrap();
        return ThermalInterfaces::new(mesh, boundary, vec![surface]).unwrap();
    }
    // The faces on the contact plane: a's look along +x, b's along -x.
    let side = |sign: f64| -> Vec<usize> {
        let on_plane = |face: &BoundaryFace| {
            on_box_face(face.centroid[0], LENGTH_M) && sign * face.outward_normal[0] > 0.0
        };
        let faces = mesh.boundary().iter().enumerate();
        faces
            .filter(|(_, face)| on_plane(face))
            .map(|(slot, _)| slot)
            .collect()
    };
    let surface = NonmatchingSurface::new(
        "slab-contact",
        side(1.0),
        side(-1.0),
        owner_resistance(),
        NonmatchingOptions {
            plane_tolerance_m: 1.0e-12,
            coverage_relative_tolerance: 1.0e-10,
            max_pair_tests: 100_000,
            max_overlap_triangles: 10_000,
        },
    )
    .unwrap();
    ThermalInterfaces::with_nonmatching(cx, mesh, boundary, Vec::new(), vec![surface]).unwrap()
}

/// The reference: fs-conduction's own contact between the two slabs on one
/// mesh, stepped by its backward Euler with the scene's step and solve
/// budgets, as fs-conduction tests/conduction_pair.rs builds it. Equal cross
/// grids take the matching contact, others the non-matching one.
fn owner_contact(cx: &fs_exec::Cx<'_>, cells: u32, cross_a: usize, cross_b: usize) -> Crossed {
    let extent = [LENGTH_M, CROSS_M, CROSS_M];
    let along = usize::try_from(cells).unwrap();
    let (slab_a, mut positions) = box_grid([along, cross_a, cross_a], extent);
    let (slab_b, positions_b) = box_grid([along, cross_b, cross_b], extent);
    let per_slab = positions.len();
    let offset = u32::try_from(per_slab).unwrap();
    let mut tets = slab_a.tets;
    tets.extend(slab_b.tets.iter().map(|tet| tet.map(|v| v + offset)));
    positions.extend(positions_b.iter().map(|&[x, y, z]| [x + LENGTH_M, y, z]));
    let mesh =
        ConductionMesh::new(TetComplex::from_tets(positions.len(), tets), positions).unwrap();
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .adiabatic_remainder()
        .finish()
        .unwrap();
    let interfaces = owner_interfaces(cx, &mesh, &boundary, cross_a == cross_b);
    let solver = BackwardEuler::uniform(
        cx,
        &mesh,
        VolumetricHeatCapacity::declared(RHO_C_J_M3_K).unwrap(),
    )
    .unwrap();
    let material = ConductivityModel::isotropic_declared(15.0).unwrap();
    let source = ScalarField::Uniform(0.0);
    let problem = ConductionProblem {
        mesh: &mesh,
        boundary: &boundary,
        material: &material,
        element_materials: None,
        source: &source,
    };
    let config = StepConfig {
        linear: LinearConfig {
            tolerance: 1.0e-13,
            max_iterations: 4000,
            restart: 20,
        },
        energy_tolerance_j: 1.0e-8,
    };
    let mut field: Vec<f64> = (0..mesh.vertex_count())
        .map(|v| if v < per_slab { HOT_K } else { COLD_K })
        .collect();
    let steps = 2 * cells * cells;
    let dt_s = 21.6 / f64::from(steps);
    let mut energy_j = 0.0;
    for _ in 0..steps {
        let step = solver
            .advance(cx, problem, Some(&interfaces), &field, dt_s, config)
            .unwrap();
        assert_eq!(step.contact_fluxes.len(), 1);
        energy_j += dt_s * step.contact_fluxes[0].heat_rate_a_to_b_w;
        field = step.temperature;
    }
    // Both contact faces lie in the plane x = L; a slab's vertices tell them apart.
    let on_plane = |v: &usize| on_box_face(mesh.positions()[*v][0], LENGTH_M);
    Crossed {
        energy_j,
        spans_k: [
            span(field[..per_slab].iter().copied()),
            span(field[per_slab..].iter().copied()),
            span((0..per_slab).filter(on_plane).map(|v| field[v])),
            span((per_slab..field.len()).filter(on_plane).map(|v| field[v])),
        ],
    }
}

/// The trace example `name` with `cells` along each slab and its step
/// t* / (2 n^2).
fn trace_scene(name: &str, cells: u32) -> String {
    let text = std::fs::read_to_string(example(name)).unwrap();
    assert_eq!(text.matches("\"cells\": [4, ").count(), 2);
    assert!(text.contains("\"step_s\": 0.675,") && text.contains("\"duration_s\": 21.6,"));
    let step_s = 21.6 / f64::from(2 * cells * cells);
    text.replace("\"cells\": [4, ", &format!("\"cells\": [{cells}, "))
        .replace("\"step_s\": 0.675,", &format!("\"step_s\": {step_s},"))
}

/// The trace example `name`, whose slab b has `cross_b` cells across, on 4
/// and 8 cells along each slab against the owner's contact on one mesh.
fn against_the_owner(name: &str, cross_b: usize) {
    // The pair's stored energy, rho c V (T_a + T_b).
    let held_j = RHO_C_J_M3_K * LENGTH_M * CROSS_M * CROSS_M * (HOT_K + COLD_K);
    for cells in [4_u32, 8] {
        let text = trace_scene(name, cells);
        scene_context(|cx, _| {
            let owner = owner_contact(cx, cells, 2, cross_b);
            let mut scene = SystemScene::from_text(&text, name).unwrap();
            let mut crossed_j = 0.0;
            let mut time_s = scene.time_s();
            while scene.advance(f64::MAX, cx).unwrap() {
                let now_s = scene.time_s();
                crossed_j += probe(&scene, "wall.heat_rate_w") * (now_s - time_s);
                time_s = now_s;
            }
            let pair = |low: &str, high: &str| (probe(&scene, low), probe(&scene, high));
            let spans_k = [
                pair("a.min_temperature_k", "a.max_temperature_k"),
                pair("b.min_temperature_k", "b.max_temperature_k"),
                pair("a.x_max.min_temperature_k", "a.x_max.max_temperature_k"),
                pair("b.x_min.min_temperature_k", "b.x_min.max_temperature_k"),
            ];
            let gap_j = crossed_j - owner.energy_j;
            let gap_k = spans_k
                .iter()
                .zip(&owner.spans_k)
                .map(|(scene, owner)| (scene.0 - owner.0).abs().max((scene.1 - owner.1).abs()))
                .fold(0.0, f64::max);
            let ledger = scene.ledger();
            let (stored, residual) = (ledger.stored_change_j(), ledger.energy_residual_j());
            eprintln!(
                "{name} n {cells}: {} steps to {time_s} s; crossed {crossed_j:.15e} owner {:.15e} J; scene - owner {gap_j:e} J; largest span difference {gap_k:e} K; faces a {:?} K, b {:?} K; {} interface points, {} W/K; stored change {stored:e} J, ledger residual {residual:e} J",
                scene.steps(),
                owner.energy_j,
                spans_k[2],
                spans_k[3],
                probe(&scene, "wall.interface_points"),
                probe(&scene, "wall.conductance_w_k"),
            );
            assert_eq!(scene.steps(), u64::from(2 * cells * cells));
            assert!(
                gap_j.abs() <= CLAIMED_GAP_J,
                "n {cells}: the scene is {gap_j:e} J from the owner's contact"
            );
            assert!(
                gap_k <= CLAIMED_GAP_K,
                "n {cells}: the scene is {gap_k:e} K from the owner's contact"
            );
            assert!(
                stored.abs() <= 1e-9 * held_j,
                "the pair changed by {stored:e} J"
            );
            assert!(
                residual.abs() <= 1e-9 * held_j,
                "ledger residual {residual:e} J"
            );
        });
    }
}

#[test]
fn two_bodies_trace_to_trace_from_a_scene_are_the_owner_contact_on_one_mesh() {
    against_the_owner("two_slabs_trace.json", 2);
}

#[test]
fn two_bodies_whose_faces_do_not_match_from_a_scene_are_the_owner_nonmatching_contact() {
    against_the_owner("two_slabs_mapped.json", 3);
}

#[test]
fn the_trace_examples_run_through_the_binary_and_report_the_declared_resistance() {
    for (name, points) in [
        ("two_slabs_trace.json", 24.0),
        ("two_slabs_mapped.json", 144.0),
    ] {
        let path = example(name);
        let (code, json, stderr) = run_binary(&["--json", "system", &path]);
        assert_eq!(code, 0, "{name}: {stderr}");
        let wall = element("wall");
        let (rate, total) = (
            value(&json, &wall, "heat_rate_w"),
            value(&json, &wall, "conductance_w_k"),
        );
        eprintln!("{name}: wall rate {rate} W, conductance {total} W/K");
        assert!(rate > 0.0, "{name}: rate {rate} W");
        // The interface points cover the face: their conductances sum to A / R''.
        assert!(
            (total / CONDUCTANCE_W_K - 1.0).abs() <= 1e-12,
            "{name}: {total} W/K"
        );
        // Both are printed exactly: a count, and the number the scene gives.
        let exactly = |key: &str, expected: f64| {
            assert_eq!(
                value(&json, &wall, key).to_bits(),
                expected.to_bits(),
                "{name}: {key}"
            );
        };
        exactly("interface_points", points);
        exactly("resistance_m2_k_w", 8.0e-4);
        assert!(
            json.contains(
                "\"declared\":{\"between\":[\"a.x_max\",\"b.x_min\"],\"resistance_m2_k_w\":0.0008,\"source\":\"declared for this example, not measured"
            ),
            "{name}: {json}"
        );
        let (code, text, stderr) = run_binary(&["system", &path]);
        assert_eq!(code, 0, "{name}: {stderr}");
        assert!(
            text.contains("Contact resistances (declared in the scene, not measured):")
                && text.contains("wall: 8e-4 m^2 K/W between a.x_max and b.x_min <- declared")
                && text.contains("a trace contact carries heat between two traces"),
            "{name}: {text}"
        );
    }
}

/// A refusal case: its name, the changes to the example in order, and what
/// the refusal says.
type Refused<'a> = (&'a str, Vec<(&'a str, &'a str)>, &'a str);

/// Every case changes the example `name`, is refused with exit code 3 and
/// says why.
fn refuse_each(name: &str, cases: Vec<Refused<'_>>) {
    let base = std::fs::read_to_string(example(name)).unwrap();
    for (case, changes, reason) in cases {
        let mut text = base.clone();
        for (from, to) in changes {
            assert!(text.contains(from), "{case}: test fixture changed");
            text = text.replacen(from, to, 1);
        }
        let (code, stderr) = refused(case, &text);
        assert_eq!(code, 3, "{case}: {stderr}");
        assert!(stderr.contains(reason), "{case}: {stderr}");
    }
}

#[test]
fn trace_contact_scenes_refuse_bad_ends_wiring_and_placement_by_name() {
    let ends = "\"between\": [\"a.x_max\", \"b.x_min\"]";
    let cases = vec![
        (
            "both",
            vec![(
                "\"traces\": [\"x_max\"]",
                "\"ports\": [\"x_max\"], \"traces\": [\"x_max\"]",
            )],
            "face x_max is under ports and under traces",
        ),
        (
            "mean",
            vec![("\"traces\": [\"x_max\"]", "\"ports\": [\"x_max\"]")],
            "between names a.x_max, which is not one of a's traces",
        ),
        (
            "missing",
            vec![(ends, "\"between\": [\"c.x_max\", \"b.x_min\"]")],
            "no element is named \"c\"",
        ),
        (
            "one body",
            vec![(ends, "\"between\": [\"a.x_max\", \"a.x_max\"]")],
            "two traces of one body",
        ),
        (
            "not a body",
            vec![(ends, "\"between\": [\"wall.a\", \"b.x_min\"]")],
            "wall is not a conduction_body",
        ),
        (
            "one end",
            vec![(ends, "\"between\": [\"a.x_max\"]")],
            "between is [\"element.face\", \"element.face\"], two traces",
        ),
        (
            "wiring",
            vec![
                (
                    "\"traces\": [\"x_min\"]",
                    "\"traces\": [\"x_min\", \"x_max\"]",
                ),
                (
                    "{\"join\": [\"wall.b\", \"b.x_min\"]}",
                    "{\"join\": [\"wall.b\", \"b.x_max\"]}",
                ),
            ],
            "wall.b is the contact's side at b.x_min; join it to that trace alone",
        ),
        (
            "apart",
            vec![(
                "\"offset_m\": [0.03, 0.0, 0.0]",
                "\"offset_m\": [0.06, 0.0, 0.0]",
            )],
            "the traces do not match",
        ),
    ];
    refuse_each("two_slabs_trace.json", cases);
}

#[test]
fn mapped_trace_contact_scenes_refuse_bad_placement_and_budgets_by_name() {
    let cases = vec![
        (
            "apart mapped",
            vec![(
                "\"offset_m\": [0.03, 0.0, 0.0]",
                "\"offset_m\": [0.06, 0.0, 0.0]",
            )],
            "contact vertices leave the declared plane tolerance",
        ),
        (
            "not mapped",
            vec![(
                ",\n      \"mapped\": {\"plane_tolerance_m\": 1e-12, \"coverage_relative_tolerance\": 1e-10, \"max_pair_tests\": 100000, \"max_overlap_triangles\": 10000}",
                "",
            )],
            "`mapped_contact` joins traces that do not match (without mapped the two faces must have the same points and triangles; mapped gives the refinement's budgets for faces that do not)",
        ),
        (
            "budget",
            vec![("\"max_pair_tests\": 100000", "\"max_pair_tests\": 0")],
            "max_pair_tests is a whole number from 1 to",
        ),
        (
            "no budget",
            vec![(", \"max_overlap_triangles\": 10000", "")],
            "missing max_overlap_triangles",
        ),
    ];
    refuse_each("two_slabs_mapped.json", cases);
}

#[test]
fn trace_contact_scenes_refuse_bad_resistances_sources_and_fields_by_name() {
    let cases = vec![
        (
            "resistance",
            vec![(
                "\"resistance_m2_k_w\": 0.0008",
                "\"resistance_m2_k_w\": -0.0008",
            )],
            "resistance_m2_k_w must be positive and finite",
        ),
        (
            "source",
            vec![(
                ",\n      \"source\": \"declared for this example, not measured: the 0.5 W/K wall of two_slabs.json over its 20 x 20 mm face\"",
                "",
            )],
            "source says where resistance_m2_k_w comes from",
        ),
        (
            "field",
            vec![(
                "\"resistance_m2_k_w\": 0.0008,",
                "\"resistance_m2_k_w\": 0.0008, \"emissivity\": 0.5,",
            )],
            "unsupported field \"emissivity\"",
        ),
    ];
    refuse_each("two_slabs_trace.json", cases);
}
