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
//!
//! A body whose conductivity is a material's: examples/system/
//! steel316_heated_bar.json, a 40 x 20 x 20 mm bar of the NASA/TP-2010-216435
//! 316 stainless card from 290 K, 20 W into one end for 60 s. The reference
//! is fs-conduction's backward Euler on the same box with the source's
//! printed k(T) table, stepped by its Newton endpoint ([`owner_bar`]). The
//! claims, declared before the first run: the compiled pack holds the
//! printed table; the field's least and greatest temperature equal the
//! owner's bit for bit after each of the 30 steps; the 1200 J put in is
//! stored and the ledger closes, each to 1e-9 of the bar's start energy; the
//! hot end passes the table's 300 K knot and stays below its 400 K end. A
//! step that would carry the bar past 400 K is refused and changes nothing.
//! Then the example through the binary with its report's no-claim sentence,
//! that sentence for a scene that also has a card-backed contact, and the
//! refusals.
//!
//! Two bodies of that material trace to trace: examples/system/
//! steel316_blocks_trace.json, the trace examples' slabs and contact with
//! the card's k(T) and the bar's heat capacity, from 390 K and 290 K, and
//! the same with slab b 3 x 3 across and the contact mapped. The reference
//! is [`owner_contact`] with the printed table, stepped by the owner's Newton
//! endpoint with the contact on one mesh. The claims, declared before the
//! first run, on 4 and 8 cells along each slab: the energy crossed within
//! [`CLAIMED_STEEL_GAP_J`] and the eight least and greatest temperatures
//! within [`CLAIMED_STEEL_GAP_K`] of the owner's, and the pair's stored
//! energy and the ledger closed to 1e-9 of what the pair holds. The first
//! run started slab b at 300 K, a knot of the table, and the scene refused
//! its first step (the owner has no tangent where a whole element sits on a
//! knot); that start is now a refusal test.
//!
//! A body whose heat capacity and density are the card's too:
//! examples/system/steel316_enthalpy_bar.json, the heated bar with no
//! declared heat capacity and its mass fixed at 290 K. The reference is
//! fs-conduction's paired material prepared from the same pack's card with
//! the scene's keys, context and span, and its enthalpy owner on the same
//! box stepped by `advance` ([`owner_enthalpy_bar`]). The claims, declared
//! before the first run: the compiled pack holds the printed Cp and density
//! tables; the field's least and greatest temperature equal the owner's bit
//! for bit after each of the 30 steps; the reported density is the owner's
//! and the printed table's at 290 K to 1e-12; the 1200 J put in is stored
//! and the ledger closes, each within 1e-6 J; the hot end passes the 300 K
//! knot and stays below 400 K. Then the example through the binary and the
//! refusals.
//!
//! A contact that follows its card: the trace examples' slabs, a from 380 K
//! and b, with half the heat capacity, from 300 K, through a card authored
//! for the test (not measured) whose resistance falls with the mean of the
//! two traces' temperatures, read from a pack the test writes beside the
//! scene. The reference is fs-conduction's own temperature contact from the
//! same card on one mesh ([`owner_law_contact`]), stepped by its Newton
//! endpoint. The claims, declared before the first run, are on
//! [`a_contact_that_follows_its_card_from_a_scene_is_the_owner_temperature_contact`].
//! Then examples/system/two_blocks_atkins_following.json, the Atkins example
//! with `follows`: that card's resistance does not vary with temperature
//! inside its span, so the example is the read-once one in every number and
//! differs in what it refuses, a state off the card at any step. Then the
//! refusals.

use std::f64::consts::{FRAC_PI_2, PI};

use super::*;
use fs_cli::SystemScene;
use fs_conduction::fixtures::{box_grid, on_box_face};
use fs_conduction::interface::temperature::{
    ContactIterationConfig, TemperatureContact, TemperatureContactCoordinates, TemperatureContacts,
};
use fs_conduction::interface::{NonmatchingOptions, NonmatchingSurface};
use fs_conduction::transient::VolumetricHeatCapacity;
use fs_conduction::transient::backward_euler::{BackwardEuler, NonlinearStepConfig, StepConfig};
use fs_conduction::transient::enthalpy::material::{
    ExactMaterialEnthalpySource, MaterialEnthalpyTransport,
};
use fs_conduction::transient::enthalpy::{EnthalpyBudget, EnthalpyStepConfig};
use fs_conduction::{
    AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS, AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY, BoundaryFace,
    ConductionMesh, ConductionProblem, ConductivityModel, ConductivityTable, InterfaceFacePair,
    InterfaceResistance, InterfaceSurface, LinearConfig, ScalarField, ThermalBc, ThermalBoundary,
    ThermalBoundaryBuilder, ThermalInterfaces,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialCard, MaterialStateId,
    NormalizedInterfacePack, NormalizedMaterialCardPack, NormalizedPack, ObservationDataset,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy,
    SurfaceSpec, SystemContext, UncertaintyModel,
};
use fs_material::phase::{HeatCapacitySourceQuery, SolidLiquidPhase};
use fs_material::state_point::SPECIFIC_HEAT_CAPACITY_DIMS;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, QuantitySpec};
use fs_rep_mesh::TetComplex;
use fs_solver::NewtonKrylovConfig;

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

/// The two slabs' bulk: one conductivity, the owner's Newton policy when it
/// varies with temperature, one volumetric heat capacity, and the slabs'
/// start temperatures.
struct Bulk {
    conductivity: ConductivityModel,
    newton: Option<NonlinearStepConfig>,
    rho_c_j_m3_k: f64,
    start_k: [f64; 2],
}

/// The slab examples' bulk: 15 W/m K and 3.6e6 J/m^3 K, from 400 K and 300 K.
fn declared_bulk() -> Bulk {
    Bulk {
        conductivity: ConductivityModel::isotropic_declared(15.0).unwrap(),
        newton: None,
        rho_c_j_m3_k: RHO_C_J_M3_K,
        start_k: [HOT_K, COLD_K],
    }
}

/// The reference: fs-conduction's own contact between the two slabs on one
/// mesh, stepped by its backward Euler with the scene's step and solve
/// budgets, as fs-conduction tests/conduction_pair.rs builds it, and by its
/// Newton endpoint (`advance_nonlinear`) when the bulk's conductivity varies.
/// Equal cross grids take the matching contact, others the non-matching one.
fn owner_contact(
    cx: &fs_exec::Cx<'_>,
    bulk: &Bulk,
    cells: u32,
    cross_a: usize,
    cross_b: usize,
) -> Crossed {
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
        VolumetricHeatCapacity::declared(bulk.rho_c_j_m3_k).unwrap(),
    )
    .unwrap();
    let source = ScalarField::Uniform(0.0);
    let problem = ConductionProblem {
        mesh: &mesh,
        boundary: &boundary,
        material: &bulk.conductivity,
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
        .map(|v| bulk.start_k[usize::from(v >= per_slab)])
        .collect();
    let steps = 2 * cells * cells;
    let dt_s = 21.6 / f64::from(steps);
    let mut energy_j = 0.0;
    for _ in 0..steps {
        let step = match bulk.newton {
            None => solver.advance(cx, problem, Some(&interfaces), &field, dt_s, config),
            Some(policy) => solver
                .advance_nonlinear(cx, problem, Some(&interfaces), &field, dt_s, config, policy)
                .map(|solution| solution.step),
        }
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
    let text = |cells| trace_scene(name, cells);
    let claimed = [CLAIMED_GAP_J, CLAIMED_GAP_K];
    against(name, &text, &declared_bulk(), cross_b, claimed);
}

/// The scene `text(cells)`, read as the example `name`, against the owner's
/// contact on one mesh with `bulk` and slab b `cross_b` cells across, on 4
/// and 8 cells along each slab: the energy crossed within `claimed[0]` J,
/// the eight least and greatest temperatures within `claimed[1]` K.
fn against(
    name: &str,
    text: &dyn Fn(u32) -> String,
    bulk: &Bulk,
    cross_b: usize,
    claimed: [f64; 2],
) {
    // The pair's stored energy, rho c V (T_a + T_b).
    let held_j =
        bulk.rho_c_j_m3_k * LENGTH_M * CROSS_M * CROSS_M * (bulk.start_k[0] + bulk.start_k[1]);
    for cells in [4_u32, 8] {
        let text = text(cells);
        scene_context(|cx, _| {
            let owner = owner_contact(cx, bulk, cells, 2, cross_b);
            let mut scene = SystemScene::from_text(&text, example(name)).unwrap();
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
                gap_j.abs() <= claimed[0],
                "n {cells}: the scene is {gap_j:e} J from the owner's contact"
            );
            assert!(
                gap_k <= claimed[1],
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
        ("steel316_blocks_trace.json", 24.0),
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

/// The source the pack of examples/system/two_blocks_atkins.json was
/// compiled from, and the pack beside the example.
const ATKINS_SOURCE: &str =
    "data/matdb/seed-v1/atkins-1964-al6061t6-lathe-0p30um-vacuum-interface/interface.tsv";
const ATKINS_PACK: &str = "atkins1964/al6061t6-lathe-0p30um-vacuum.fsintpk";
/// The card example's source of resistance, as its file writes it.
const CARD_FIELDS: &str = "\"pack\": \"atkins1964/al6061t6-lathe-0p30um-vacuum.fsintpk\",\n      \"state\": [\n        {\"name\": \"normal_pressure\", \"value\": \"313 kPa\", \"kind\": \"legacy\"},\n        {\"name\": \"T\", \"value\": \"297.5 K\", \"kind\": \"legacy\"}\n      ]";
/// How far a card-backed contact may lie from the printed table.
const CLAIMED_CARD_GAP: f64 = 1.0e-12;

/// The table the seed source prints, from its `curve` row: the joint's
/// resistance [m^2 K/W] at each tabulated pressure [Pa].
fn printed_knots() -> Vec<(f64, f64)> {
    let path = format!("{}/../../{ATKINS_SOURCE}", env!("CARGO_MANIFEST_DIR"));
    let text = std::fs::read_to_string(path).unwrap();
    let row = text
        .lines()
        .find(|line| line.starts_with("curve\tcontact-resistance\t"))
        .expect("the source has the resistance curve");
    let fields: Vec<&str> = row.split('\t').collect();
    assert_eq!(
        (fields[4], fields[5], fields[6], fields[8]),
        ("normal_pressure", "Pa", "m^2*K/W", "linear")
    );
    fields[7]
        .split(',')
        .map(|knot| {
            let (pascal, resistance) = knot.split_once(':').unwrap();
            (pascal.parse().unwrap(), resistance.parse().unwrap())
        })
        .collect()
}

/// The card example with its changes made, in order.
fn card_scene(case: &str, changes: &[(&str, &str)]) -> String {
    let mut text = std::fs::read_to_string(example("two_blocks_atkins.json")).unwrap();
    for (from, to) in changes {
        assert!(text.contains(from), "{case}: test fixture changed");
        text = text.replacen(from, to, 1);
    }
    text
}

/// A scene beside the card example, run to its end: the energy that crossed
/// its contact and the temperatures a trace scene is compared in.
fn crossed(cx: &fs_exec::Cx<'_>, text: &str) -> Crossed {
    let mut scene = SystemScene::from_text(text, example("two_blocks_atkins.json")).unwrap();
    let mut energy_j = 0.0;
    let mut time_s = scene.time_s();
    while scene.advance(f64::MAX, cx).unwrap() {
        let now_s = scene.time_s();
        energy_j += probe(&scene, "wall.heat_rate_w") * (now_s - time_s);
        time_s = now_s;
    }
    let pair = |low: &str, high: &str| (probe(&scene, low), probe(&scene, high));
    Crossed {
        energy_j,
        spans_k: [
            pair("a.min_temperature_k", "a.max_temperature_k"),
            pair("b.min_temperature_k", "b.max_temperature_k"),
            pair("a.x_max.min_temperature_k", "a.x_max.max_temperature_k"),
            pair("b.x_min.min_temperature_k", "b.x_min.max_temperature_k"),
        ],
    }
}

#[test]
fn a_card_backed_contact_reports_the_printed_resistance_at_and_between_the_tabulated_pressures() {
    let knots = printed_knots();
    // The compiled card holds the printed table and interpolates it linearly.
    let pack =
        NormalizedInterfacePack::from_bytes(&std::fs::read(example(ATKINS_PACK)).unwrap()).unwrap();
    let claims = pack
        .card()
        .claims_for(AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY);
    assert_eq!(claims.len(), 1);
    assert_eq!(claims[0].1.interpolation, InterpolationPolicy::LinearInside);
    let PropertyValue::Curve {
        knots: compiled, ..
    } = &claims[0].1.value
    else {
        panic!("the resistance is not a curve");
    };
    assert_eq!(compiled, &knots, "the pack is not the source's table");
    // Two tabulated pressures, and the pressure midway between the first two.
    let between = (
        0.5 * (knots[0].0 + knots[1].0),
        0.5 * (knots[0].1 + knots[1].1),
    );
    for (pascal, printed) in [knots[1], knots[2], between] {
        let text = card_scene("pressure", &[("\"313 kPa\"", &format!("\"{pascal} Pa\""))]);
        let scene = SystemScene::from_text(&text, example("two_blocks_atkins.json")).unwrap();
        let (resistance, conductance) = (
            probe(&scene, "wall.resistance_m2_k_w"),
            probe(&scene, "wall.conductance_w_k"),
        );
        eprintln!(
            "{pascal} Pa: card {resistance:e} m^2 K/W, printed {printed:e}; {conductance} W/K, A / R'' {}",
            CROSS_M * CROSS_M / printed
        );
        assert!(
            (resistance / printed - 1.0).abs() <= CLAIMED_CARD_GAP,
            "{pascal} Pa: {resistance:e} against the printed {printed:e}"
        );
        assert!(
            (conductance * printed / (CROSS_M * CROSS_M) - 1.0).abs() <= CLAIMED_CARD_GAP,
            "{pascal} Pa: {conductance} W/K"
        );
    }
}

#[test]
fn a_card_backed_contact_is_the_declared_contact_at_the_printed_resistance() {
    let printed = printed_knots()[1];
    assert_eq!(printed.0.to_bits(), 313_000.0_f64.to_bits());
    let declared = card_scene(
        "declared",
        &[(
            CARD_FIELDS,
            &format!(
                "\"resistance_m2_k_w\": {:e},\n      \"source\": \"the card's printed value at 313 kPa\"",
                printed.1
            ),
        )],
    );
    scene_context(|cx, _| {
        let card = crossed(cx, &card_scene("card", &[]));
        let typed = crossed(cx, &declared);
        let gap_j = card.energy_j - typed.energy_j;
        let gap_k = card
            .spans_k
            .iter()
            .zip(&typed.spans_k)
            .map(|(card, typed)| (card.0 - typed.0).abs().max((card.1 - typed.1).abs()))
            .fold(0.0, f64::max);
        eprintln!(
            "card {:.15e} J, declared {:.15e} J; card - declared {gap_j:e} J; largest span difference {gap_k:e} K; faces a {:?} K, b {:?} K",
            card.energy_j, typed.energy_j, card.spans_k[2], card.spans_k[3]
        );
        assert!(card.energy_j > 0.0, "nothing crossed: {} J", card.energy_j);
        assert!(gap_j.abs() <= CLAIMED_GAP_J, "{gap_j:e} J apart");
        assert!(gap_k <= CLAIMED_GAP_K, "{gap_k:e} K apart");
    });
}

#[test]
fn the_card_example_runs_through_the_binary_and_reports_its_card() {
    let path = example("two_blocks_atkins.json");
    let (code, json, stderr) = run_binary(&["--json", "system", &path]);
    assert_eq!(code, 0, "{stderr}");
    assert!(value(&json, &element("wall"), "heat_rate_w") > 0.0);
    for member in [
        "\"card\":{\"between\":[\"a.x_max\",\"b.x_min\"],\"resistance_m2_k_w\":0.00031606,\"pack_id\":\"atkins-1964-al6061t6-lathe-0p30um-vacuum-interface\",\"pack_hash\":\"",
        "\"resistance_receipt_hash\":\"",
        "\"selected_claim\":\"",
        "\"state\":{\"T\":297.5,\"normal_pressure\":313000.0}}",
    ] {
        assert!(json.contains(member), "{member} not in {json}");
    }
    assert!(!json.contains("\"declared\":"), "{json}");
    let (code, text, stderr) = run_binary(&["system", &path]);
    assert_eq!(code, 0, "{stderr}");
    assert!(
        text.contains("Contact resistances from interface cards")
            && text.contains(
                "wall: 3.1606e-4 m^2 K/W between a.x_max and b.x_min <- atkins-1964-al6061t6-lathe-0p30um-vacuum-interface"
            )
            && text.contains("an interface card's value read once at the declared state")
            && !text.contains("declared in the scene, not measured"),
        "{text}"
    );
}

/// A refusal case of the card example: its name, the pack whose bytes stand
/// beside the scene, the changes in order, and what the refusal says.
type CardRefused<'a> = (&'a str, &'a str, Vec<(&'a str, &'a str)>, &'a str);

/// A changed copy of the card example, written beside a copy of `pack`
/// under the name the example gives its pack, and refused: the exit code
/// and the diagnostic.
fn card_refused(case: &str, pack: &str, changes: &[(&str, &str)]) -> (i32, String) {
    let path = scratch(case, &card_scene(case, changes));
    let beside = path.parent().unwrap().join("atkins1964");
    std::fs::create_dir_all(&beside).unwrap();
    std::fs::copy(
        example(pack),
        beside.join("al6061t6-lathe-0p30um-vacuum.fsintpk"),
    )
    .unwrap();
    let (code, stdout, stderr) = run_binary(&["system", path.to_str().unwrap()]);
    assert!(stdout.is_empty(), "{stdout}");
    (code, stderr)
}

#[test]
fn card_backed_contact_scenes_refuse_two_sources_states_off_the_card_and_other_packs_by_name() {
    let source = "\"pack\": \"atkins1964/al6061t6-lathe-0p30um-vacuum.fsintpk\",";
    let pressure =
        "        {\"name\": \"normal_pressure\", \"value\": \"313 kPa\", \"kind\": \"legacy\"},\n";
    let one_of_two = "from a card (pack with state): one of the two, and not both";
    // matdb's own refusal of a state off the card, through the owner.
    let outside = "extrapolation refused: the query point is outside the validity of every one";
    let no_source = format!(",\n      {CARD_FIELDS}");
    let cases: Vec<CardRefused<'_>> = vec![
        (
            "both",
            ATKINS_PACK,
            vec![(
                source,
                "\"resistance_m2_k_w\": 0.0008, \"source\": \"typed\", \"pack\": \"atkins1964/al6061t6-lathe-0p30um-vacuum.fsintpk\",",
            )],
            one_of_two,
        ),
        ("neither", ATKINS_PACK, vec![(&no_source, "")], one_of_two),
        (
            "low pressure",
            ATKINS_PACK,
            vec![("\"313 kPa\"", "\"50 kPa\"")],
            outside,
        ),
        (
            "hot",
            ATKINS_PACK,
            vec![("\"297.5 K\"", "\"350 K\"")],
            outside,
        ),
        (
            "no pressure",
            ATKINS_PACK,
            vec![(pressure, "")],
            "fs-matdb refused \"area-specific-thermal-contact-resistance\"",
        ),
        (
            "absent",
            ATKINS_PACK,
            vec![("vacuum.fsintpk\",", "absent.fsintpk\",")],
            "elements.wall.pack:",
        ),
        (
            "material pack",
            "steel1018/steel-1018.fsmcdpk",
            vec![],
            "elements.wall.pack:",
        ),
        (
            "selection",
            ATKINS_PACK,
            vec![(
                "\"state\": [",
                "\"selection\": \"best\",\n      \"state\": [",
            )],
            "elements.wall.selection:",
        ),
    ];
    for (case, pack, changes, reason) in cases {
        let (code, stderr) = card_refused(case, pack, &changes);
        eprintln!("{case}: exit {code}: {}", stderr.trim_end());
        assert_eq!(code, 3, "{case}: {stderr}");
        assert!(stderr.contains(reason), "{case}: {stderr}");
    }
}

/// The trace example with the stick, tip, head and case of struck_mode.json
/// beside its slabs: a scene with a thermal and a mechanical contact.
fn slabs_beside_a_strike() -> String {
    let slabs = std::fs::read_to_string(example("two_slabs_trace.json")).unwrap();
    let strike = std::fs::read_to_string(example("struck_mode.json")).unwrap();
    let part = |from: &str, to: &str| {
        let start = strike.find(from).expect("the part begins") + from.len();
        let end = start + strike[start..].find(to).expect("the part ends");
        &strike[start..end]
    };
    let elements = part("\"elements\": {\n", "\n  },\n  \"connections\"");
    let connections = part("\"connections\": [\n", "\n  ]\n}");
    let both = slabs
        .replacen(
            "\"elements\": {\n",
            &format!("\"elements\": {{\n{elements},\n"),
            1,
        )
        .replacen(
            "\"connections\": [\n",
            &format!("\"connections\": [\n{connections},\n"),
            1,
        );
    assert!(both.contains("\"contact_interface\"") && both.contains("\"trace_contact\""));
    both
}

#[test]
fn a_scene_with_a_contact_or_a_card_does_not_say_it_has_none() {
    // What follows "no material-card resolution" in the list of what a scene
    // has none of, up to "phase change".
    let plain = ", contact, friction,";
    let mechanical = ", mechanical contact, friction,";
    let carded = " except a trace contact's interface card, mechanical contact, friction,";
    let struck = ", friction,";
    let read = |name: &str| std::fs::read_to_string(example(name)).unwrap();
    for (name, text, listed) in [
        ("two_slabs.json", read("two_slabs.json"), plain),
        ("heated_body.json", read("heated_body.json"), plain),
        (
            "two_slabs_trace.json",
            read("two_slabs_trace.json"),
            mechanical,
        ),
        (
            "two_slabs_mapped.json",
            read("two_slabs_mapped.json"),
            mechanical,
        ),
        (
            "two_blocks_atkins.json",
            read("two_blocks_atkins.json"),
            carded,
        ),
        (
            "two_blocks_atkins_following.json",
            read("two_blocks_atkins_following.json"),
            carded,
        ),
        ("struck_mode.json", read("struck_mode.json"), struck),
        ("two_slabs_trace.json", slabs_beside_a_strike(), struck),
    ] {
        let scene = SystemScene::from_text(&text, example(name)).unwrap();
        let clause = format!(
            "; no material-card resolution{listed} phase change or step-size error control;"
        );
        for report in [scene.render_text(), scene.render_json()] {
            let said: Vec<String> = report
                .match_indices("no material-card resolution")
                .map(|(at, _)| report[at..].chars().take(120).collect())
                .collect();
            assert!(report.contains(&clause), "{name}: {said:?}");
            assert_eq!(said.len(), 1, "{name}: {said:?}");
        }
        eprintln!("{name}: {clause}");
    }
}

/// The stainless example, its pack and the seed source the pack is compiled
/// from.
const STEEL_EXAMPLE: &str = "steel316_heated_bar.json";
const STEEL_PACK: &str = "steel316/stainless-316-nasa.fsmcdpk";
const STEEL_SOURCE: &str = "data/matdb/seed-v1/stainless-316-nasa-tp216435/properties.tsv";
/// The example's heat rate into x_min [W], its step [s] and its steps, and
/// its declared volumetric heat capacity [J/(m^3 K)]: the table's density
/// at 300 K, 7960.74 kg/m^3, times its specific heat there, 467.69 J/(kg K),
/// as the scene writes it.
const STEEL_POWER_W: f64 = 20.0;
const STEEL_STEP_S: f64 = 2.0;
const STEEL_STEPS: u64 = 30;
const STEEL_RHO_C_J_M3_K: f64 = 3_723_158.490_6;
const STEEL_START_K: f64 = 290.0;

/// The conductivity table the seed source prints, from its `curve` row:
/// [K, W/(m K)].
fn printed_conductivity() -> Vec<(f64, f64)> {
    let path = format!("{}/../../{STEEL_SOURCE}", env!("CARGO_MANIFEST_DIR"));
    let text = std::fs::read_to_string(path).unwrap();
    let row = text
        .lines()
        .find(|line| line.starts_with("curve\tthermal-conductivity-273-400k\t"))
        .expect("the source has the conductivity curve");
    let fields: Vec<&str> = row.split('\t').collect();
    assert_eq!(fields[3], "thermal_conductivity");
    assert_eq!(
        (fields[4], fields[5], fields[6], fields[8]),
        ("temperature", "K", "W/m/K", "linear")
    );
    fields[7]
        .split(',')
        .map(|knot| {
            let (kelvin, conductivity) = knot.split_once(':').unwrap();
            (kelvin.parse().unwrap(), conductivity.parse().unwrap())
        })
        .collect()
}

/// The reference: fs-conduction's backward Euler on the example's box with
/// `printed` as its conductivity, stepped by its Newton endpoint under the
/// owner's default policy, with the scene's step and budgets and the heat
/// rate as a uniform flux over x_min, whose area is the faces' summed in
/// boundary order, as a port sums them. The field's least and greatest
/// temperature after each step.
fn owner_bar(cx: &fs_exec::Cx<'_>, printed: Vec<(f64, f64)>) -> Vec<(f64, f64)> {
    let (complex, positions) = box_grid([8, 2, 2], [0.04, 0.02, 0.02]);
    let mesh = ConductionMesh::new(complex, positions).unwrap();
    let heated = |face: &BoundaryFace| on_box_face(face.centroid[0], 0.0);
    let area_m2: f64 = mesh
        .boundary()
        .iter()
        .filter(|face| heated(face))
        .map(|face| face.area)
        .sum();
    let flux = ThermalBc::neumann(-STEEL_POWER_W / area_m2).unwrap();
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region("x_min", heated, flux)
        .unwrap()
        .adiabatic_remainder()
        .finish()
        .unwrap();
    let capacity = VolumetricHeatCapacity::declared(STEEL_RHO_C_J_M3_K).unwrap();
    let solver = BackwardEuler::uniform(cx, &mesh, capacity).unwrap();
    let material =
        ConductivityModel::isotropic(ConductivityTable::declared_curve(printed).unwrap());
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
    let mut field = vec![STEEL_START_K; mesh.vertex_count()];
    let mut spans = Vec::new();
    for _ in 0..STEEL_STEPS {
        let policy = NonlinearStepConfig::default();
        let solved = solver
            .advance_nonlinear(cx, problem, None, &field, STEEL_STEP_S, config, policy)
            .unwrap();
        field = solved.step.temperature;
        spans.push(span(field.iter().copied()));
    }
    spans
}

#[test]
fn a_body_of_a_card_material_is_the_owner_stepping_the_printed_conductivity_by_newton() {
    let printed = printed_conductivity();
    let pack = NormalizedMaterialCardPack::from_bytes(&std::fs::read(example(STEEL_PACK)).unwrap())
        .unwrap();
    let claims = pack.card().claims().claims_for("thermal-conductivity");
    assert_eq!(claims.len(), 1);
    let PropertyValue::Curve {
        knots: compiled, ..
    } = &claims[0].1.value
    else {
        panic!("the conductivity is not a curve");
    };
    assert_eq!(compiled, &printed, "the pack is not the source's table");
    let text = std::fs::read_to_string(example(STEEL_EXAMPLE)).unwrap();
    // The bar's energy at the start, rho c V T, and the heat put in.
    let start_j = STEEL_RHO_C_J_M3_K * 0.04 * 0.02 * 0.02 * STEEL_START_K;
    let put_in_j = STEEL_POWER_W * 60.0;
    scene_context(|cx, _| {
        let owner = owner_bar(cx, printed);
        let mut scene = SystemScene::from_text(&text, example(STEEL_EXAMPLE)).unwrap();
        let mut index = 0;
        while scene.advance(f64::MAX, cx).unwrap() {
            let low = probe(&scene, "bar.min_temperature_k");
            let high = probe(&scene, "bar.max_temperature_k");
            assert!(
                low.to_bits() == owner[index].0.to_bits()
                    && high.to_bits() == owner[index].1.to_bits(),
                "step {index}: the scene's {low} .. {high} K against the owner's {:?} K",
                owner[index]
            );
            index += 1;
        }
        let ledger = scene.ledger();
        let (stored, residual) = (ledger.stored_change_j(), ledger.energy_residual_j());
        let (low, high) = owner[owner.len() - 1];
        eprintln!(
            "{} steps; bar {low} .. {high} K, x_min {} K; stored change {stored} J of {put_in_j} J in, ledger residual {residual:e} J",
            scene.steps(),
            probe(&scene, "bar.x_min.temperature_k"),
        );
        assert_eq!(scene.steps(), STEEL_STEPS);
        assert_eq!(index, owner.len());
        assert!(
            (stored - put_in_j).abs() <= 1e-9 * start_j,
            "stored {stored} J"
        );
        assert!(
            residual.abs() <= 1e-9 * start_j,
            "ledger residual {residual:e} J"
        );
        // The run crosses the table's interior knot and stays on the table.
        assert!(300.0 < high && high < 400.0, "{high} K");
    });
}

#[test]
fn the_material_body_example_runs_through_the_binary_and_reports_its_source() {
    let path = example(STEEL_EXAMPLE);
    let (code, json, stderr) = run_binary(&["--json", "system", &path]);
    assert_eq!(code, 0, "{stderr}");
    assert!(value(&json, &element("bar"), "x_min.temperature_k") > STEEL_START_K);
    for member in [
        "\"material_bindings\":{\"bar\":[\"steel316\"]}",
        "\"has_conductivity\":true",
        "\"property\":\"thermal-conductivity\"",
    ] {
        assert!(json.contains(member), "{member} not in {json}");
    }
    let (code, text, stderr) = run_binary(&["system", &path]);
    assert_eq!(code, 0, "{stderr}");
    assert!(text.contains("bar <- steel316"), "{text}");
    for said in [
        "source resolution is limited to declared pure-phase thermal laws; no contact, friction,",
        "with one declared conductivity or a scene material's k(T) curve, stepped by \
         fs-conduction's Newton endpoint and refused off the curve's span, and one declared \
         volumetric heat capacity",
        "no internal source, temperature-dependent heat capacity (a body takes only a \
         material's k(T)), imported mesh, radiation or contact trace",
    ] {
        assert!(text.contains(said), "{said} not in {text}");
    }
    assert!(
        !text.contains("no material-card resolution")
            && !text.contains("temperature-dependent conductivity"),
        "{text}"
    );
}

/// The card example with the stainless bar, its material and its heater
/// beside the two blocks: a card-backed contact and a material body.
fn blocks_beside_a_bar() -> String {
    let blocks = std::fs::read_to_string(example("two_blocks_atkins.json")).unwrap();
    let bar = std::fs::read_to_string(example(STEEL_EXAMPLE)).unwrap();
    let part = |from: &str, to: &str| {
        let start = bar.find(from).expect("the part begins");
        let end = start + bar[start..].find(to).expect("the part ends");
        &bar[start..end]
    };
    let materials = part("  \"materials\": {", "  \"elements\"");
    let element = part("    \"bar\": {", "\n  },\n  \"connections\"");
    let heater = part("    {\"prescribe\"", "\n  ]");
    let elements = format!("{materials}  \"elements\": {{\n{element},\n");
    let connections = format!("  \"connections\": [\n{heater},\n");
    blocks
        .replacen("  \"elements\": {\n", &elements, 1)
        .replacen("  \"connections\": [\n", &connections, 1)
}

#[test]
fn a_scene_with_a_material_body_and_a_card_backed_contact_names_both_sources() {
    let text = blocks_beside_a_bar();
    let scene = SystemScene::from_text(&text, example("two_blocks_atkins.json")).unwrap();
    let report = scene.render_text();
    assert!(
        report.contains(
            "source resolution is limited to declared pure-phase thermal laws and trace \
             contacts' interface cards; no mechanical contact, friction,"
        ) && !report.contains("no material-card resolution"),
        "{report}"
    );
}

#[test]
fn material_body_scenes_refuse_mixed_sources_off_span_starts_and_bad_policies_by_name() {
    let base = std::fs::read_to_string(example(STEEL_EXAMPLE)).unwrap();
    let material = "\"material\": \"steel316\",";
    let newton = |policy: &str| format!("{material} \"newton\": {policy},");
    let cases = [
        (
            "both",
            material,
            "\"material\": \"steel316\", \"conductivity_w_m_k\": 15.0,".to_owned(),
            "the conductivity is conductivity_w_m_k or a material's, not both",
        ),
        (
            "unknown",
            material,
            "\"material\": \"steel317\",".to_owned(),
            "unknown material \"steel317\"",
        ),
        (
            "no conductivity",
            ",\n      \"conductivity_property\": \"thermal-conductivity\",\n      \"conductivity_kind\": \"thermal-conductivity\"",
            String::new(),
            "has no conductivity",
        ),
        (
            "cold",
            "\"temperature_k\": 290.0,",
            "\"temperature_k\": 250.0,".to_owned(),
            "temperature_k 250 K lies outside material \"steel316\"'s conductivity, 273 K to 400 K",
        ),
        (
            "newton field",
            material,
            newton("{\"damping\": 0.5}"),
            "unsupported field \"damping\"",
        ),
        (
            "newton budget",
            material,
            newton("{\"max_iterations\": 0}"),
            "max_iterations is a whole number from 1 to 1000",
        ),
        (
            "newton rtol",
            material,
            newton("{\"residual_rtol\": 2.0}"),
            "positive Newton budget, rtol in (0,1)",
        ),
    ];
    for (case, from, to, reason) in cases {
        assert!(base.contains(from), "{case}: test fixture changed");
        let text = base.replacen(from, &to, 1);
        let Err(error) = SystemScene::from_text(&text, example(STEEL_EXAMPLE)) else {
            panic!("{case}: admitted");
        };
        eprintln!("{case}: exit {}: {}", error.exit_code, error.message);
        assert!(error.message.contains(reason), "{case}: {}", error.message);
    }
    // A declared conductivity takes no Newton policy.
    let slabs = two_slabs(4);
    let declared = "\"temperature_k\": 400.0,";
    assert!(slabs.contains(declared));
    let (code, stderr) = refused(
        "newton declared",
        &slabs.replacen(declared, "\"temperature_k\": 400.0, \"newton\": {},", 1),
    );
    assert_eq!(code, 3, "{stderr}");
    assert!(
        stderr.contains("a body with conductivity_w_m_k takes none"),
        "{stderr}"
    );
    // 2000 W carries the heated end past the table's 400 K in the first
    // step: the owner refuses the step and the scene does not move.
    let hot = base.replacen("\"heat_rate_w\": 20.0", "\"heat_rate_w\": 2000.0", 1);
    let mut scene = SystemScene::from_text(&hot, example(STEEL_EXAMPLE)).unwrap();
    let before = scene.render_json();
    scene_context(|cx, _| {
        let error = scene.advance(f64::MAX, cx).unwrap_err();
        eprintln!("2000 W: exit {}: {}", error.exit_code, error.message);
        assert_eq!(error.exit_code, fs_cli::exit::REFUSED);
        assert!(error.message.contains("LineSearchFailed"), "{error}");
        assert_eq!(scene.render_json(), before);
        assert_eq!(scene.steps(), 0);
    });
}

/// Two 316 stainless slabs trace to trace: the trace examples' slabs and
/// contact with the card's k(T), from 390 K and 290 K.
const STEEL_TRACE: &str = "steel316_blocks_trace.json";
/// Both sides stop each Newton solve at 1e-10 J plus 1e-10 of the step's
/// first residual, so they agree to that resolution, not to rounding: the
/// claims on the energy crossed [J] and on the eight temperatures [K].
const CLAIMED_STEEL_GAP_J: f64 = 1.0e-6;
const CLAIMED_STEEL_GAP_K: f64 = 1.0e-6;

/// The stainless slabs' bulk: the source's printed k(T) table under the
/// owner's default Newton policy, the bar's declared heat capacity.
fn steel_bulk() -> Bulk {
    let table = ConductivityTable::declared_curve(printed_conductivity()).unwrap();
    Bulk {
        conductivity: ConductivityModel::isotropic(table),
        newton: Some(NonlinearStepConfig::default()),
        rho_c_j_m3_k: STEEL_RHO_C_J_M3_K,
        start_k: [390.0, 290.0],
    }
}

/// The stainless trace example with slab b 3 x 3 across and the contact
/// mapped with the budgets of two_slabs_mapped.json.
fn steel_mapped(cells: u32) -> String {
    let text = trace_scene(STEEL_TRACE, cells);
    let b = format!("\"cells\": [{cells}, 2, 2]");
    let at = text.rfind(&b).expect("slab b's cells");
    let wall = "\"resistance_m2_k_w\": 0.0008,";
    assert!(text.contains(wall));
    let mapped = "\n      \"mapped\": {\"plane_tolerance_m\": 1e-12, \"coverage_relative_tolerance\": \
                  1e-10, \"max_pair_tests\": 100000, \"max_overlap_triangles\": 10000},";
    let (before, after) = (&text[..at], &text[at + b.len()..]);
    let text = format!("{before}\"cells\": [{cells}, 3, 3]{after}");
    text.replacen(wall, &format!("{wall}{mapped}"), 1)
}

#[test]
fn two_bodies_of_a_card_material_trace_to_trace_are_the_owner_newton_contact_on_one_mesh() {
    let text = |cells| trace_scene(STEEL_TRACE, cells);
    let claimed = [CLAIMED_STEEL_GAP_J, CLAIMED_STEEL_GAP_K];
    against(STEEL_TRACE, &text, &steel_bulk(), 2, claimed);
}

#[test]
fn two_bodies_of_a_card_material_on_grids_that_do_not_match_are_the_owner_newton_contact() {
    let claimed = [CLAIMED_STEEL_GAP_J, CLAIMED_STEEL_GAP_K];
    against(STEEL_TRACE, &steel_mapped, &steel_bulk(), 3, claimed);
}

/// Slab b from 300 K, a knot of the printed table. The step the kernel
/// solves through the contact asks the body for its tangent, and the owner
/// refuses a tangent where all four vertices of an element sit on a knot,
/// where k(T) has no two-sided slope: the first step is refused and the
/// scene does not move. The owner alone steps the same field
/// (`owner_contact`).
#[test]
fn a_card_material_slab_started_on_a_knot_of_its_table_is_refused_in_a_contact() {
    let text = trace_scene(STEEL_TRACE, 4);
    let start = "\"temperature_k\": 290.0,";
    assert_eq!(text.matches(start).count(), 1);
    let on_knot = text.replacen(start, "\"temperature_k\": 300.0,", 1);
    let mut scene = SystemScene::from_text(&on_knot, example(STEEL_TRACE)).unwrap();
    let before = scene.render_json();
    scene_context(|cx, _| {
        let bulk = Bulk {
            start_k: [390.0, 300.0],
            ..steel_bulk()
        };
        let owner = owner_contact(cx, &bulk, 4, 2, 2);
        assert!(owner.energy_j > 0.0);
        let error = scene.advance(f64::MAX, cx).unwrap_err();
        eprintln!("300 K: exit {}: {}", error.exit_code, error.message);
        assert_eq!(error.exit_code, fs_cli::exit::REFUSED);
        assert!(
            error.message.contains("no unique two-sided derivative"),
            "{error}"
        );
        assert_eq!(scene.render_json(), before);
        assert_eq!(scene.steps(), 0);
    });
}

/// The stainless bar with its heat capacity and density from the card as
/// well as its conductivity.
const STEEL_ENTHALPY_EXAMPLE: &str = "steel316_enthalpy_bar.json";
/// The temperature at which the example's box holds the bar's mass [K].
const STEEL_REFERENCE_K: f64 = 290.0;
/// The two energy claims [J]: the owner's 1e-8 J gate on each of the 30
/// steps, with room for rounding.
const CLAIMED_ENTHALPY_GAP_J: f64 = 1.0e-6;

/// The table of the seed source's `curve` row `id`: [K, value].
fn printed_curve(id: &str) -> Vec<(f64, f64)> {
    let path = format!("{}/../../{STEEL_SOURCE}", env!("CARGO_MANIFEST_DIR"));
    let text = std::fs::read_to_string(path).unwrap();
    let head = format!("curve\t{id}\t");
    let row = text
        .lines()
        .find(|line| line.starts_with(&head))
        .unwrap_or_else(|| panic!("the source has no curve {id}"));
    let fields: Vec<&str> = row.split('\t').collect();
    assert_eq!(
        (fields[4], fields[5], fields[8]),
        ("temperature", "K", "linear")
    );
    fields[7]
        .split(',')
        .map(|knot| {
            let (kelvin, value) = knot.split_once(':').unwrap();
            (kelvin.parse().unwrap(), value.parse().unwrap())
        })
        .collect()
}

/// The point the example's materials section reads the card at: its five
/// source flags, each a dimensionless number, and the absolute temperature.
fn steel_point(temperature_k: f64) -> QueryPoint {
    let flag = QuantitySpec::dimensional(Dims([0; 6]));
    let mut point = QueryPoint::new();
    for (name, value) in [
        ("source-generic-316-identity", 1.0),
        ("source-uns-designation-known", 0.0),
        ("source-heat-treatment-known", 0.0),
        ("source-product-form-known", 0.0),
        ("source-pressure-known", 0.0),
    ] {
        point = point.with_quantity(name, flag, value).unwrap();
    }
    let absolute = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    point
        .with_quantity("temperature", absolute, temperature_k)
        .unwrap()
}

/// What the reference reports of the enthalpy bar.
struct EnthalpyBar {
    /// The least and greatest nodal temperature after each step [K].
    spans: Vec<(f64, f64)>,
    /// The density the masses were fixed by [kg/m^3].
    density_kg_m3: f64,
}

/// The reference: fs-conduction's paired material prepared from `card`
/// with the example's keys, context and span (273 K to 400 K, 0 J/kg at
/// 273 K) and the density read at 290 K, its enthalpy owner on the
/// example's box, and the heat rate as a uniform flux over x_min, stepped
/// by `advance` from the chart's enthalpy of 290 K under the Newton policy
/// the scene's defaults give (32 updates, 1e-10 relative, 1e-10 J) and the
/// scene's 1e-8 J gate.
fn owner_enthalpy_bar(cx: &fs_exec::Cx<'_>, card: &MaterialCard) -> EnthalpyBar {
    let heat_capacity = PropertyKey::new("specific-heat-capacity", SPECIFIC_HEAT_CAPACITY_DIMS);
    let density = PropertyKey::new("density", Density::DIMS);
    let conductivity = PropertyKey::with_quantity(
        "thermal-conductivity",
        QuantitySpec::semantic(SemanticType::new(
            QuantityKind::ThermalConductivity,
            ValueForm::Static,
        )),
    );
    let (lower, upper) = (steel_point(273.0), steel_point(400.0));
    let reference = steel_point(STEEL_REFERENCE_K);
    let transport = MaterialEnthalpyTransport::from_isotropic_card_exact(
        card,
        ExactMaterialEnthalpySource {
            caloric: HeatCapacitySourceQuery {
                heat_capacity: &heat_capacity,
                density: &density,
                temperature_axis: "temperature",
                lower: &lower,
                upper: &upper,
                selection: SelectionPolicy::SingleClaimOnly,
            },
            conductivity: &conductivity,
            reference: &reference,
            phase: SolidLiquidPhase::Solid,
            lower_specific_enthalpy_j_kg: 0.0,
        },
    )
    .unwrap();
    let (complex, positions) = box_grid([8, 2, 2], [0.04, 0.02, 0.02]);
    let mesh = ConductionMesh::new(complex, positions).unwrap();
    let heated = |face: &BoundaryFace| on_box_face(face.centroid[0], 0.0);
    let area_m2: f64 = mesh
        .boundary()
        .iter()
        .filter(|face| heated(face))
        .map(|face| face.area)
        .sum();
    let flux = ThermalBc::neumann(-STEEL_POWER_W / area_m2).unwrap();
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region("x_min", heated, flux)
        .unwrap()
        .adiabatic_remainder()
        .finish()
        .unwrap();
    let budget = EnthalpyBudget {
        max_vertices: mesh.vertex_count(),
        max_elements: mesh.element_count(),
    };
    let solver = transport.prepare_uniform(cx, &mesh, budget).unwrap();
    let source = ScalarField::Uniform(0.0);
    let problem = ConductionProblem {
        mesh: &mesh,
        boundary: &boundary,
        material: transport.conductivity(),
        element_materials: None,
        source: &source,
    };
    let config = EnthalpyStepConfig {
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1.0e-10,
            relative_tolerance: 1.0e-10,
            ..NewtonKrylovConfig::default()
        },
        max_newton_iterations: 32,
        energy_tolerance_j: 1.0e-8,
    };
    let start = transport
        .caloric()
        .state_at_temperature(STEEL_START_K)
        .unwrap()
        .specific_enthalpy_j_kg();
    let mut field = vec![start; mesh.vertex_count()];
    let mut spans = Vec::new();
    for _ in 0..STEEL_STEPS {
        let step = solver
            .advance(cx, problem, None, &field, STEEL_STEP_S, config)
            .unwrap();
        spans.push(span(step.temperature.iter().copied()));
        field = step.specific_enthalpy_j_kg;
    }
    EnthalpyBar {
        spans,
        density_kg_m3: transport.reference_density_kg_m3(),
    }
}

/// The printed table's value at `t` by linear interpolation.
fn interpolate(table: &[(f64, f64)], t: f64) -> f64 {
    let i = table.iter().position(|&(k, _)| k >= t).unwrap().max(1);
    let ((t0, v0), (t1, v1)) = (table[i - 1], table[i]);
    v0 + (t - t0) / (t1 - t0) * (v1 - v0)
}

#[test]
fn a_body_of_a_card_heat_capacity_is_the_owner_enthalpy_field_on_the_same_card() {
    let pack = NormalizedMaterialCardPack::from_bytes(&std::fs::read(example(STEEL_PACK)).unwrap())
        .unwrap();
    for (property, id) in [
        ("specific-heat-capacity", "specific-heat-capacity-273-400k"),
        ("density", "density-273-400k"),
    ] {
        let claims = pack.card().claims().claims_for(property);
        assert_eq!(claims.len(), 1, "{property}");
        let PropertyValue::Curve {
            knots: compiled, ..
        } = &claims[0].1.value
        else {
            panic!("{property} is not a curve");
        };
        assert_eq!(compiled, &printed_curve(id), "the pack's {property}");
    }
    let text = std::fs::read_to_string(example(STEEL_ENTHALPY_EXAMPLE)).unwrap();
    let put_in_j = STEEL_POWER_W * 60.0;
    scene_context(|cx, _| {
        let owner = owner_enthalpy_bar(cx, pack.card());
        let mut scene = SystemScene::from_text(&text, example(STEEL_ENTHALPY_EXAMPLE)).unwrap();
        let mut index = 0;
        while scene.advance(f64::MAX, cx).unwrap() {
            let low = probe(&scene, "bar.min_temperature_k");
            let high = probe(&scene, "bar.max_temperature_k");
            assert!(
                low.to_bits() == owner.spans[index].0.to_bits()
                    && high.to_bits() == owner.spans[index].1.to_bits(),
                "step {index}: the scene's {low} .. {high} K against the owner's {:?} K",
                owner.spans[index]
            );
            index += 1;
        }
        let density = probe(&scene, "bar.reference_density_kg_m3");
        let printed = interpolate(&printed_curve("density-273-400k"), STEEL_REFERENCE_K);
        let ledger = scene.ledger();
        let (stored, residual) = (ledger.stored_change_j(), ledger.energy_residual_j());
        let (low, high) = owner.spans[owner.spans.len() - 1];
        eprintln!(
            "{} steps; bar {low} .. {high} K, x_min {} K; density {density} kg/m^3 (printed table at 290 K {printed}); stored change {stored} J of {put_in_j} J in, ledger residual {residual:e} J",
            scene.steps(),
            probe(&scene, "bar.x_min.temperature_k"),
        );
        assert_eq!(scene.steps(), STEEL_STEPS);
        assert_eq!(index, owner.spans.len());
        assert_eq!(density.to_bits(), owner.density_kg_m3.to_bits());
        assert!(
            (density - printed).abs() <= 1e-12 * printed,
            "{density} against {printed}"
        );
        assert!(
            (stored - put_in_j).abs() <= CLAIMED_ENTHALPY_GAP_J,
            "stored {stored} J"
        );
        assert!(
            residual.abs() <= CLAIMED_ENTHALPY_GAP_J,
            "ledger residual {residual:e} J"
        );
        // The run crosses the chart's interior knot and stays on the chart.
        assert!(300.0 < high && high < 400.0, "{high} K");
    });
}

#[test]
fn the_card_heat_capacity_example_runs_through_the_binary_and_says_what_the_card_gives() {
    let path = example(STEEL_ENTHALPY_EXAMPLE);
    let (code, json, stderr) = run_binary(&["--json", "system", &path]);
    assert_eq!(code, 0, "{stderr}");
    assert!(value(&json, &element("bar"), "x_min.temperature_k") > STEEL_START_K);
    assert!(value(&json, &element("bar"), "reference_density_kg_m3") > 7960.0);
    assert!(
        json.contains("\"material_bindings\":{\"bar\":[\"steel316\"]}"),
        "{json}"
    );
    let (code, text, stderr) = run_binary(&["system", &path]);
    assert_eq!(code, 0, "{stderr}");
    assert!(text.contains("bar <- steel316"), "{text}");
    for said in [
        "and one declared volumetric heat capacity or a material card's Cp(T) and density, \
         whose specific enthalpy fs-conduction's enthalpy owner stores on the card's chart, \
         with masses fixed by the card's density at the body's reference temperature and no \
         thermal expansion, each port",
        "no internal source, imported mesh, radiation or contact trace",
    ] {
        assert!(text.contains(said), "{said} not in {text}");
    }
    assert!(
        !text.contains("temperature-dependent heat capacity"),
        "{text}"
    );
}

#[test]
fn card_heat_capacity_scenes_refuse_missing_references_off_chart_states_and_contacts_by_name() {
    let base = std::fs::read_to_string(example(STEEL_ENTHALPY_EXAMPLE)).unwrap();
    let reference = "\"reference_temperature_k\": 290.0,";
    let cases = [
        (
            "no reference",
            reference,
            String::new(),
            "needs reference_temperature_k, the temperature at which box_m holds its mass",
        ),
        (
            "declared too",
            reference,
            format!("{reference} \"volumetric_heat_capacity_j_m3_k\": 3723158.4906,"),
            "a body with volumetric_heat_capacity_j_m3_k takes none",
        ),
        (
            "reference off the span",
            reference,
            "\"reference_temperature_k\": 250.0,".to_owned(),
            "reference temperature must be supported",
        ),
        (
            "start off the chart",
            "\"temperature_k\": 290.0,",
            "\"temperature_k\": 250.0,".to_owned(),
            "temperature_k 250 K:",
        ),
        (
            "no conductivity",
            ",\n      \"conductivity_property\": \"thermal-conductivity\",\n      \"conductivity_kind\": \"thermal-conductivity\"",
            String::new(),
            "the declaration names no conductivity_property",
        ),
        (
            "both",
            "\"material\": \"steel316\",",
            "\"material\": \"steel316\", \"conductivity_w_m_k\": 15.0,".to_owned(),
            "the conductivity is conductivity_w_m_k or a material's, not both",
        ),
    ];
    for (case, from, to, reason) in cases {
        assert!(base.contains(from), "{case}: test fixture changed");
        let text = base.replacen(from, &to, 1);
        let Err(error) = SystemScene::from_text(&text, example(STEEL_ENTHALPY_EXAMPLE)) else {
            panic!("{case}: admitted");
        };
        eprintln!("{case}: exit {}: {}", error.exit_code, error.message);
        assert!(error.message.contains(reason), "{case}: {}", error.message);
    }
    // A trace contact joins bodies of declared heat capacity: slab a of the
    // stainless trace example with its heat capacity from the card.
    let slabs = trace_scene(STEEL_TRACE, 4);
    let declared = "\"volumetric_heat_capacity_j_m3_k\": 3723158.4906,";
    assert_eq!(slabs.matches(declared).count(), 2);
    let text = slabs.replacen(declared, reference, 1);
    let Err(error) = SystemScene::from_text(&text, example(STEEL_TRACE)) else {
        panic!("contact: admitted");
    };
    eprintln!("contact: exit {}: {}", error.exit_code, error.message);
    assert!(
        error
            .message
            .contains("fs-conduction's contacts join bodies with volumetric_heat_capacity_j_m3_k"),
        "{}",
        error.message
    );
    // 2000 W carries the heated end off the chart in the first step: the
    // owner refuses the step and the scene does not move.
    let hot = base.replacen("\"heat_rate_w\": 20.0", "\"heat_rate_w\": 2000.0", 1);
    let mut scene = SystemScene::from_text(&hot, example(STEEL_ENTHALPY_EXAMPLE)).unwrap();
    let before = scene.render_json();
    scene_context(|cx, _| {
        let error = scene.advance(f64::MAX, cx).unwrap_err();
        eprintln!("2000 W: exit {}: {}", error.exit_code, error.message);
        assert_eq!(error.exit_code, fs_cli::exit::REFUSED);
        assert_eq!(scene.render_json(), before);
        assert_eq!(scene.steps(), 0);
    });
}

/// The card of the following-contact tests, authored for them and not
/// measured: its resistance [m^2 K/W] at three means of the two traces'
/// temperatures [K], linear between them and with no answer outside.
const AUTHORED_KNOTS: [(f64, f64); 3] = [(280.0, 1.4e-3), (345.0, 8.0e-4), (430.0, 5.0e-4)];
/// The state the authored scenes declare for the card [K].
const AUTHORED_REFERENCE_K: f64 = 320.0;
/// The authored scenes' slabs: a from 380 K with the trace examples' heat
/// capacity, b from 300 K with half of it, so that the interface's means
/// pass the card's middle knot; 12 steps of 1.8 s.
const LAW_START_K: [f64; 2] = [380.0, 300.0];
const LAW_RHO_C_J_M3_K: [f64; 2] = [3.6e6, 1.8e6];
const LAW_STEP_S: f64 = 1.8;
const LAW_STEPS: u32 = 12;
/// How far a scene whose contact follows its card may lie from the owner's
/// temperature contact on one mesh.
const CLAIMED_LAW_GAP_J: f64 = 1.0e-7;
const CLAIMED_LAW_GAP_K: f64 = 1.0e-7;
/// The least the following must move one of the eight reported temperatures
/// from the same scene with its card read once.
const CLAIMED_LAW_MOVE_K: f64 = 0.1;
/// The example whose contact follows the Atkins card.
const FOLLOWING_EXAMPLE: &str = "two_blocks_atkins_following.json";

/// One resistance claim as a compiled interface pack, authored for these
/// tests: a scene reads its bytes, the owner's contact takes its card.
fn authored_pack(
    id: &str,
    value: PropertyValue,
    validity: ValidityDomain,
) -> NormalizedInterfacePack {
    let source = fs_blake3::hash_bytes(id.as_bytes());
    let provenance = || Provenance {
        source: format!(
            "authored for the system-scene following-contact tests ({id}), not measured"
        ),
        license: "internal-test-use".to_string(),
        artifact: Some(source),
    };
    let mut claims = ClaimSet::new();
    let observation = claims
        .register_observation(ObservationDataset {
            specimen: id.to_string(),
            method: "authored numbers, no instrument".to_string(),
            artifact: source,
            caveats: "not measurements".to_string(),
            provenance: provenance(),
        })
        .unwrap();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(
                AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY,
                AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS,
            ),
            value,
            validity,
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![observation],
            provenance: provenance(),
        })
        .unwrap();
    let claims = NormalizedPack::new(
        id,
        "authored-in-test",
        source,
        "internal-test-use",
        claims,
        Vec::new(),
        Vec::new(),
    )
    .unwrap();
    let side = |chemistry: &str| SurfaceSpec {
        material: MaterialStateId {
            chemistry: chemistry.to_string(),
            phase: "solid".to_string(),
            process: "as-fixtured".to_string(),
            revision: 0,
        },
        texture_frame: "interface-normal-x".to_string(),
    };
    let context = SystemContext {
        medium: "dry".to_string(),
        third_body: None,
        environment: "unspecified".to_string(),
        history: "unaged".to_string(),
    };
    NormalizedInterfacePack::new(side("slab-a"), side("slab-b"), context, claims).unwrap()
}

/// The authored card: [`AUTHORED_KNOTS`] over the plain coordinate `T`.
fn authored_curve() -> NormalizedInterfacePack {
    authored_pack(
        "authored-resistance-over-mean-temperature",
        PropertyValue::Curve {
            abscissa: "T".to_string(),
            abscissa_dims: fs_conduction::TEMPERATURE_DIMS,
            knots: AUTHORED_KNOTS.to_vec(),
            dims: AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS,
        },
        ValidityDomain::unconstrained().with("T", AUTHORED_KNOTS[0].0, AUTHORED_KNOTS[2].0),
    )
}

/// [`AUTHORED_KNOTS`] by hand at `kelvin`.
fn authored_resistance(kelvin: f64) -> f64 {
    let segment = usize::from(kelvin > AUTHORED_KNOTS[1].0);
    let ((t0, r0), (t1, r1)) = (AUTHORED_KNOTS[segment], AUTHORED_KNOTS[segment + 1]);
    r0 + (r1 - r0) * (kelvin - t0) / (t1 - t0)
}

/// The trace example `name` as the authored scene: the slabs of
/// [`LAW_START_K`] and [`LAW_RHO_C_J_M3_K`], [`LAW_STEPS`] steps, and the
/// contact from the authored pack at [`AUTHORED_REFERENCE_K`], following its
/// `T` or (`follows` false) read once. Written with the pack beside it;
/// its path.
fn authored_scene(case: &str, name: &str, follows: bool) -> std::path::PathBuf {
    let mut text = std::fs::read_to_string(example(name)).unwrap();
    let follows = if follows {
        ",\n      \"follows\": \"T\""
    } else {
        ""
    };
    let card = format!(
        "\"pack\": \"authored/curve.fsintpk\",\n      \"state\": [{{\"name\": \"T\", \"value\": \"{AUTHORED_REFERENCE_K} K\", \"kind\": \"legacy\"}}]{follows}"
    );
    for (from, to) in [
        ("\"temperature_k\": 400.0,", "\"temperature_k\": 380.0,"),
        (
            "\"volumetric_heat_capacity_j_m3_k\": 3600000.0,\n      \"temperature_k\": 300.0,",
            "\"volumetric_heat_capacity_j_m3_k\": 1800000.0,\n      \"temperature_k\": 300.0,",
        ),
        ("\"step_s\": 0.675,", "\"step_s\": 1.8,"),
        (
            "\"resistance_m2_k_w\": 0.0008,\n      \"source\": \"declared for this example, not measured: the 0.5 W/K wall of two_slabs.json over its 20 x 20 mm face\"",
            card.as_str(),
        ),
    ] {
        assert_eq!(
            text.matches(from).count(),
            1,
            "{case}: test fixture changed"
        );
        text = text.replacen(from, to, 1);
    }
    let path = scratch(case, &text);
    let beside = path.parent().unwrap().join("authored");
    std::fs::create_dir_all(&beside).unwrap();
    std::fs::write(beside.join("curve.fsintpk"), authored_curve().to_bytes()).unwrap();
    path
}

/// A scene file run to its end: what a trace scene is compared in.
fn crossed_at(cx: &fs_exec::Cx<'_>, path: &std::path::Path) -> (Crossed, SystemScene) {
    let text = std::fs::read_to_string(path).unwrap();
    let mut scene = SystemScene::from_text(&text, path).unwrap();
    let mut energy_j = 0.0;
    let mut time_s = scene.time_s();
    while scene.advance(f64::MAX, cx).unwrap() {
        let now_s = scene.time_s();
        energy_j += probe(&scene, "wall.heat_rate_w") * (now_s - time_s);
        time_s = now_s;
    }
    let pair = |low: &str, high: &str| (probe(&scene, low), probe(&scene, high));
    let spans_k = [
        pair("a.min_temperature_k", "a.max_temperature_k"),
        pair("b.min_temperature_k", "b.max_temperature_k"),
        pair("a.x_max.min_temperature_k", "a.x_max.max_temperature_k"),
        pair("b.x_min.min_temperature_k", "b.x_min.max_temperature_k"),
    ];
    (Crossed { energy_j, spans_k }, scene)
}

/// The largest difference between two runs' eight reported temperatures [K].
fn span_gap_k(one: &Crossed, other: &Crossed) -> f64 {
    one.spans_k
        .iter()
        .zip(&other.spans_k)
        .map(|(one, other)| (one.0 - other.0).abs().max((one.1 - other.1).abs()))
        .fold(0.0, f64::max)
}

/// Both authored slabs in one mesh with unshared faces at x = L, slab a
/// first: the mesh, its elements' vertices, and how many vertices and
/// elements are slab a's.
fn law_mesh(cross_b: usize) -> (ConductionMesh, Vec<[u32; 4]>, [usize; 2]) {
    let extent = [LENGTH_M, CROSS_M, CROSS_M];
    let (slab_a, mut positions) = box_grid([4, 2, 2], extent);
    let (slab_b, positions_b) = box_grid([4, cross_b, cross_b], extent);
    let of_a = [positions.len(), slab_a.tets.len()];
    let offset = u32::try_from(positions.len()).unwrap();
    let mut tets = slab_a.tets;
    tets.extend(slab_b.tets.iter().map(|tet| tet.map(|v| v + offset)));
    positions.extend(positions_b.iter().map(|&[x, y, z]| [x + LENGTH_M, y, z]));
    let complex = TetComplex::from_tets(positions.len(), tets.clone());
    (ConductionMesh::new(complex, positions).unwrap(), tets, of_a)
}

/// The owner's temperature contact from `card` on that mesh, at the scene's
/// reference point and over the mean of the two traces: its pointwise
/// contact on equal cross grids, its non-matching one, with the mapped
/// example's budgets, on others.
fn owner_law<'c>(
    mesh: &ConductionMesh,
    card: &'c InterfaceSystemCard,
    cross_b: usize,
) -> TemperatureContact<'c> {
    let point = QueryPoint::new().with("T", AUTHORED_REFERENCE_K).unwrap();
    let mean = TemperatureContactCoordinates::Mean("T".to_string());
    let policy = SelectionPolicy::SingleClaimOnly;
    if cross_b == 2 {
        let pairs = ThermalInterfaces::coincident_face_pairs(mesh).unwrap();
        let pairs = pairs.into_iter().map(|pair| {
            if mesh.boundary()[pair.side_a].outward_normal[0] > 0.0 {
                pair
            } else {
                InterfaceFacePair {
                    side_a: pair.side_b,
                    side_b: pair.side_a,
                }
            }
        });
        return TemperatureContact::new_pointwise(
            "wall",
            pairs.collect(),
            card,
            &point,
            mean,
            policy,
        )
        .unwrap();
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
    let options = NonmatchingOptions {
        plane_tolerance_m: 1.0e-12,
        coverage_relative_tolerance: 1.0e-10,
        max_pair_tests: 100_000,
        max_overlap_triangles: 10_000,
    };
    let (a, b) = (side(1.0), side(-1.0));
    TemperatureContact::new_nonmatching("wall", a, b, card, &point, mean, policy, options).unwrap()
}

/// The reference of the authored scenes: both slabs in one mesh and the
/// owner's temperature contact from the authored card, stepped by the
/// owner's `advance_capacity_newton` with the scenes' step and budgets. The
/// energy is what slab a's field lost, `rho c` times the integral of its P1
/// temperature.
fn owner_law_contact(cx: &fs_exec::Cx<'_>, cross_b: usize) -> Crossed {
    let (mesh, tets, [per_slab, elements_a]) = law_mesh(cross_b);
    let held_by_a = |field: &[f64]| -> f64 {
        let integral: f64 = (0..elements_a)
            .map(|e| {
                let mean = tets[e].iter().map(|&v| field[v as usize]).sum::<f64>() / 4.0;
                mesh.element_volume(e) * mean
            })
            .sum();
        LAW_RHO_C_J_M3_K[0] * integral
    };
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .adiabatic_remainder()
        .finish()
        .unwrap();
    let pack = authored_curve();
    let faces = 10_000;
    let set = TemperatureContacts::with_nonmatching(
        cx,
        &mesh,
        &boundary,
        Vec::new(),
        Vec::new(),
        vec![owner_law(&mesh, pack.card(), cross_b)],
        faces,
        3 * faces,
    )
    .unwrap();
    let capacities: Vec<VolumetricHeatCapacity> = (0..mesh.element_count())
        .map(|element| {
            VolumetricHeatCapacity::declared(LAW_RHO_C_J_M3_K[usize::from(element >= elements_a)])
                .unwrap()
        })
        .collect();
    let engine = BackwardEuler::per_element(cx, &mesh, &capacities).unwrap();
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
    let newton = NonlinearStepConfig {
        residual_rtol: 1.0e-12,
        residual_atol_j: 1.0e-11,
        ..Default::default()
    };
    let coupling = ContactIterationConfig {
        max_iterations: 64,
        relaxation: 0.5,
        residual_tolerance_j: 1.0e-8,
    };
    let mut field: Vec<f64> = (0..mesh.vertex_count())
        .map(|v| LAW_START_K[usize::from(v >= per_slab)])
        .collect();
    let start_j = held_by_a(&field);
    for _ in 0..LAW_STEPS {
        let step = set
            .advance_capacity_newton(
                cx, &engine, problem, &field, LAW_STEP_S, config, newton, coupling, faces,
            )
            .unwrap();
        field = step.conduction.temperature;
    }
    let on_plane = |v: &usize| on_box_face(mesh.positions()[*v][0], LENGTH_M);
    Crossed {
        energy_j: start_j - held_by_a(&field),
        spans_k: [
            span(field[..per_slab].iter().copied()),
            span(field[per_slab..].iter().copied()),
            span((0..per_slab).filter(on_plane).map(|v| field[v])),
            span((per_slab..field.len()).filter(on_plane).map(|v| field[v])),
        ],
    }
}

/// The authored scene made from the trace example `name`, whose slab b has
/// `cross_b` cells across, against the owner's temperature contact, against
/// itself with the card read once, and its reported resistances by hand.
fn following_against_the_owner(name: &str, cross_b: usize) {
    scene_context(|cx, _| {
        let owner = owner_law_contact(cx, cross_b);
        let (scene, ran) = crossed_at(cx, &authored_scene("following.json", name, true));
        let (once, _) = crossed_at(cx, &authored_scene("once.json", name, false));
        let gap_j = scene.energy_j - owner.energy_j;
        let gap_k = span_gap_k(&scene, &owner);
        let moved_k = span_gap_k(&scene, &once);
        let reported = |key: &str| probe(&ran, &format!("wall.{key}"));
        let (declared, low, high) = (
            reported("resistance_m2_k_w"),
            reported("min_resistance_m2_k_w"),
            reported("max_resistance_m2_k_w"),
        );
        // The card falls with the mean, and every interface point's mean
        // lies between the means of the two faces' least and of their
        // greatest temperatures.
        let (face_a, face_b) = (scene.spans_k[2], scene.spans_k[3]);
        let bracket = (
            authored_resistance(0.5 * (face_a.1 + face_b.1)),
            authored_resistance(0.5 * (face_a.0 + face_b.0)),
        );
        eprintln!(
            "{name} following its card: {} steps; crossed {:.15e} owner {:.15e} J, scene - owner {gap_j:e} J; largest span difference {gap_k:e} K; read once crossed {:.15e} J and lies {moved_k:e} K away; resistance declared {declared:e}, over the points {low:e} to {high:e}, by hand between {:e} and {:e} m^2 K/W; faces a {face_a:?} K, b {face_b:?} K",
            ran.steps(),
            scene.energy_j,
            owner.energy_j,
            once.energy_j,
            bracket.0,
            bracket.1,
        );
        assert_eq!(ran.steps(), u64::from(LAW_STEPS));
        assert!(
            gap_j.abs() <= CLAIMED_LAW_GAP_J,
            "{gap_j:e} J from the owner"
        );
        assert!(gap_k <= CLAIMED_LAW_GAP_K, "{gap_k:e} K from the owner");
        assert!(
            moved_k >= CLAIMED_LAW_MOVE_K,
            "the card moved {moved_k:e} K"
        );
        let by_hand = authored_resistance(AUTHORED_REFERENCE_K);
        assert!((declared / by_hand - 1.0).abs() <= 1e-12, "{declared:e}");
        let slack = 1e-12 * by_hand;
        assert!(
            bracket.0 - slack <= low && low <= high && high <= bracket.1 + slack,
            "{low:e} to {high:e} outside {bracket:?}"
        );
        assert!(
            high < declared,
            "the contact stayed at its declared {declared:e}"
        );
    });
}

/// Claims, declared before the first run: over 12 steps the scene whose
/// contact follows the authored card crosses the energy slab a loses in the
/// owner's pointwise temperature contact on one mesh within
/// [`CLAIMED_LAW_GAP_J`] and reports its eight temperatures within
/// [`CLAIMED_LAW_GAP_K`]; the same scene with the card read once lies at
/// least [`CLAIMED_LAW_MOVE_K`] away; the declared resistance is the card's
/// at 320 K by hand, and the least and greatest over the interface points
/// lie inside the bracket the two faces' reported temperatures give by hand.
#[test]
fn a_contact_that_follows_its_card_from_a_scene_is_the_owner_temperature_contact() {
    following_against_the_owner("two_slabs_trace.json", 2);
}

/// The same claims on cross grids that do not match, against the owner's
/// non-matching temperature contact.
#[test]
fn a_mapped_contact_that_follows_its_card_is_the_owner_nonmatching_temperature_contact() {
    following_against_the_owner("two_slabs_mapped.json", 3);
}

/// Claims, declared before the first run. The Atkins card's resistance does
/// not vary with `T` inside its span, so the example whose contact follows
/// it crosses the read-once example's energy within [`CLAIMED_GAP_J`] and
/// reports its temperatures within [`CLAIMED_GAP_K`], and its least and
/// greatest resistance are the declared one within 1e-12. Its reports say
/// that it follows. With block a from 330 K the two faces' mean is 310 K,
/// above the card's 305.85 K: the first step is refused with matdb's text
/// and the scene does not move.
#[test]
fn the_following_example_is_the_read_once_example_and_refuses_a_state_off_its_card() {
    let path = example(FOLLOWING_EXAMPLE);
    let text = std::fs::read_to_string(&path).unwrap();
    assert_eq!(
        text.replacen(",\n      \"follows\": \"T\"", "", 1),
        std::fs::read_to_string(example("two_blocks_atkins.json")).unwrap(),
        "the example is the read-once example plus follows"
    );
    scene_context(|cx, _| {
        let (following, ran) = crossed_at(cx, std::path::Path::new(&path));
        let once = crossed(cx, &card_scene("read once", &[]));
        let (gap_j, gap_k) = (
            following.energy_j - once.energy_j,
            span_gap_k(&following, &once),
        );
        let reported = |key: &str| probe(&ran, &format!("wall.{key}"));
        let (declared, low, high) = (
            reported("resistance_m2_k_w"),
            reported("min_resistance_m2_k_w"),
            reported("max_resistance_m2_k_w"),
        );
        eprintln!(
            "following {:.15e} J, read once {:.15e} J; following - read once {gap_j:e} J; largest span difference {gap_k:e} K; resistance declared {declared:e}, over the points {low:e} to {high:e} m^2 K/W",
            following.energy_j, once.energy_j
        );
        assert!(following.energy_j > 0.0);
        assert!(gap_j.abs() <= CLAIMED_GAP_J, "{gap_j:e} J apart");
        assert!(gap_k <= CLAIMED_GAP_K, "{gap_k:e} K apart");
        for value in [low, high] {
            assert!((value / declared - 1.0).abs() <= 1e-12, "{value:e}");
        }
        let hot = text.replacen("\"temperature_k\": 305.0,", "\"temperature_k\": 330.0,", 1);
        assert_ne!(hot, text, "test fixture changed");
        let mut scene = SystemScene::from_text(&hot, &path).unwrap();
        let before = scene.render_json();
        let error = scene.advance(f64::MAX, cx).unwrap_err();
        eprintln!("330 K: exit {}: {}", error.exit_code, error.message);
        assert_eq!(error.exit_code, fs_cli::exit::REFUSED);
        assert!(
            error
                .message
                .starts_with("step 1 at t = 0 s: component `wall` refused in Propose")
                && error.message.contains(
                    "extrapolation refused: the query point is outside the validity of every one"
                ),
            "{}",
            error.message
        );
        assert_eq!((scene.steps(), scene.render_json()), (0, before));
    });
    let (code, json, stderr) = run_binary(&["--json", "system", &path]);
    assert_eq!(code, 0, "{stderr}");
    assert!(
        json.contains("\"state\":{\"T\":297.5,\"normal_pressure\":313000.0},\"follows\":[\"T\"]}")
            && json.contains("\"min_resistance_m2_k_w\":")
            && json.contains("asked at every interface point of every step"),
        "{json}"
    );
    let (code, report, stderr) = run_binary(&["system", &path]);
    assert_eq!(code, 0, "{stderr}");
    assert!(
        report.contains("Contact resistances that follow an interface card")
            && report.contains(
                "wall: 3.1606e-4 m^2 K/W between a.x_max and b.x_min <- atkins-1964-al6061t6-lathe-0p30um-vacuum-interface, following T"
            )
            && !report.contains("each read once at its declared state"),
        "{report}"
    );
}

/// A following contact that cannot be bound is refused before any step,
/// each case with exit code 3: by the scene for the form of `follows` and
/// for a declared number beside it, by the owner for coordinates it cannot
/// follow, and by the law for a surface over the two temperatures.
#[test]
fn following_contact_scenes_refuse_bad_forms_coordinates_and_surfaces_by_name() {
    let follows = "\"follows\": \"T\"";
    let coordinate =
        "each contact temperature coordinate needs a positive absolute-temperature value";
    for (case, to, reason) in [
        ("number", "\"follows\": 3", "elements.wall.follows:"),
        (
            "one name in a list",
            "\"follows\": [\"T\"]",
            "elements.wall.follows:",
        ),
        (
            "the same name twice",
            "\"follows\": [\"T\", \"T\"]",
            "contact temperature coordinates must be nonblank and distinct",
        ),
        (
            "no such coordinate",
            "\"follows\": \"interface_T\"",
            coordinate,
        ),
        ("a pressure", "\"follows\": \"normal_pressure\"", coordinate),
    ] {
        let text = std::fs::read_to_string(example(FOLLOWING_EXAMPLE)).unwrap();
        assert!(text.contains(follows), "{case}: test fixture changed");
        let path = scratch(case, &text.replacen(follows, to, 1));
        let beside = path.parent().unwrap().join("atkins1964");
        std::fs::create_dir_all(&beside).unwrap();
        std::fs::copy(
            example(ATKINS_PACK),
            beside.join("al6061t6-lathe-0p30um-vacuum.fsintpk"),
        )
        .unwrap();
        let (code, stdout, stderr) = run_binary(&["system", path.to_str().unwrap()]);
        eprintln!("{case}: exit {code}: {}", stderr.trim_end());
        assert!(stdout.is_empty(), "{case}: {stdout}");
        assert_eq!(code, 3, "{case}: {stderr}");
        assert!(stderr.contains(reason), "{case}: {stderr}");
    }
    // A declared number has no card to follow.
    let declared = std::fs::read_to_string(example("two_slabs_trace.json")).unwrap();
    let both = declared.replacen(
        "\"resistance_m2_k_w\": 0.0008,",
        "\"resistance_m2_k_w\": 0.0008, \"follows\": \"T\",",
        1,
    );
    assert_ne!(both, declared, "test fixture changed");
    let (code, stderr) = refused("a declared number", &both);
    eprintln!("a declared number: exit {code}: {}", stderr.trim_end());
    assert_eq!(code, 3, "{stderr}");
    assert!(stderr.contains("one of the two, and not both"), "{stderr}");
    // A surface over the two traces' temperatures: the law is not built.
    let corners = [
        [280.0, 280.0],
        [430.0, 280.0],
        [280.0, 430.0],
        [430.0, 430.0],
    ];
    let surface = authored_pack(
        "authored-resistance-over-two-temperatures",
        PropertyValue::Surface2D {
            axes: ["T".to_string(), "T_b".to_string()],
            axis_dims: [fs_conduction::TEMPERATURE_DIMS; 2],
            samples: corners
                .into_iter()
                .map(|at| (at, 1.0e-3 - 1.0e-6 * (at[0] - 280.0)))
                .collect(),
            triangles: vec![[0, 1, 2], [1, 3, 2]],
            dims: AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS,
        },
        ValidityDomain::unconstrained()
            .with("T", 280.0, 430.0)
            .with("T_b", 280.0, 430.0),
    );
    let path = authored_scene("surface.json", "two_slabs_trace.json", true);
    let text = std::fs::read_to_string(&path).unwrap().replacen(
        "\"kind\": \"legacy\"}]",
        "\"kind\": \"legacy\"}, {\"name\": \"T_b\", \"value\": \"320 K\", \"kind\": \"legacy\"}]",
        1,
    );
    std::fs::write(&path, text).unwrap();
    std::fs::write(
        path.parent().unwrap().join("authored/curve.fsintpk"),
        surface.to_bytes(),
    )
    .unwrap();
    let (code, stdout, stderr) = run_binary(&["system", path.to_str().unwrap()]);
    eprintln!("a surface: exit {code}: {}", stderr.trim_end());
    assert!(stdout.is_empty(), "{stdout}");
    assert_eq!(code, 3, "{stderr}");
    assert!(
        stderr.contains(
            "a surface over a temperature coordinate, and any other interpolation, is not built"
        ),
        "{stderr}"
    );
}
