//! Two finite-element bodies on the composition kernel, joined through a
//! contact interface card.
//!
//! Two equal slabs start at 400 K and 300 K with insulated outer ends and
//! exchange heat through a Conductance whose value is the card's contact
//! resistance over the patch.
//!
//! The claim, declared before the first run: the energy that has crossed the
//! contact by `t* = 0.1 L^2 / alpha`, compared with the one-dimensional
//! closed form (the series below), converges at order 2, at least
//! [`CLAIMED_ORDER`] between the last two rungs, on three paths: (a) the
//! pair on the kernel with matching cross grids, (b) the pair with
//! non-matching cross grids, and (c) the owner's own matching-P1 contact on
//! one mesh at the same steps. Every path is reported before any is
//! asserted.
//!
//! Two more paths join the two bodies trace to trace through the owner's own
//! contact as a `ConductanceMatrix`: (d) on matching cross grids, which is
//! path (c)'s problem, and (e) on cross grids that do not match, which is the
//! owner's non-matching contact on one mesh. Their claims are stated at
//! their tests. Run with
//! `cargo test -p fs-conduction --features component --test conduction_pair -- --nocapture`.
#![cfg(feature = "component")]

use std::f64::consts::{FRAC_PI_2, PI};

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_conduction::component::{ConductionBody, PatchPort, TraceContact};
use fs_conduction::fixtures::{box_grid, on_box_face};
use fs_conduction::interface::{NonmatchingOptions, NonmatchingSurface};
use fs_conduction::transient::VolumetricHeatCapacity;
use fs_conduction::transient::backward_euler::{BackwardEuler, StepConfig};
use fs_conduction::{
    AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS, AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY, BoundaryFace,
    ConductionError, ConductionMesh, ConductionProblem, ConductivityModel, InterfaceFacePair,
    InterfaceResistance, InterfaceSurface, LinearConfig, ScalarField, ThermalBc, ThermalBoundary,
    ThermalBoundaryBuilder, ThermalInterfaces,
};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, SurfaceSpec,
    SystemContext, UncertaintyModel,
};
use fs_mms::{RefinementLadder, fit_order};
use fs_ports::component::Component;
use fs_ports::component::conformance::{Case, check_component, check_port_tangent};
use fs_rep_mesh::TetComplex;
use fs_system::SystemBuilder;
use fs_system::elements::{Conductance, ConductanceMatrix, InterfacePoint};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let key = StreamKey {
            seed: 66,
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

/// Each slab: 30 mm along x, 20 x 20 mm across, steel-like `k` and `rho c`.
const LENGTH_M: f64 = 0.03;
const CROSS_M: f64 = 0.02;
const CONDUCTIVITY_W_M_K: f64 = 15.0;
const RHO_C_J_M3_K: f64 = 3.6e6;
/// The card's contact: Bi = 2 L / (k R'') = 5 for the half-problem.
const CONTACT_M2_K_W: f64 = 8.0e-4;
const HOT_K: f64 = 400.0;
const COLD_K: f64 = 300.0;
/// The comparison time, in units of L^2 / alpha.
const FOURIER: f64 = 0.1;
/// Cells along x on each rung; the step is t* / (2 n^2), so backward
/// Euler's first-order error is second order in h.
const RUNGS: [usize; 4] = [2, 4, 8, 16];
/// The least order between the last two rungs that meets the claim.
const CLAIMED_ORDER: f64 = 1.8;
/// Path (d) against path (c): the most the crossed energy may differ, and
/// the most any vertex's end temperature may.
const CLAIMED_GAP_J: f64 = 1.0e-9;
const CLAIMED_GAP_K: f64 = 1.0e-9;

fn config() -> StepConfig {
    StepConfig {
        linear: LinearConfig {
            tolerance: 1.0e-13,
            max_iterations: 4000,
            restart: 20,
        },
        energy_tolerance_j: 1.0e-8,
    }
}

fn leak<T>(value: T) -> &'static T {
    Box::leak(Box::new(value))
}

fn end_time_s() -> f64 {
    FOURIER * LENGTH_M * LENGTH_M * RHO_C_J_M3_K / CONDUCTIVITY_W_M_K
}

/// The energy that has crossed the contact by `t*`, from the closed form.
///
/// For equal slabs the field is antisymmetric about the contact around the
/// mean temperature, so each slab is the plane wall insulated at its outer
/// end with a Robin end of `h = 2 / R''`: the excess it still holds is
/// `sum C_n sin(b_n) / b_n exp(-b_n^2 Fo)` of the initial one, with
/// `b_n tan b_n = Bi` and `C_n = 4 sin b_n / (2 b_n + sin 2 b_n)`.
fn closed_form_crossed_j() -> f64 {
    let biot = 2.0 * LENGTH_M / (CONDUCTIVITY_W_M_K * CONTACT_M2_K_W);
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
            if mid * mid.tan() < biot {
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
    let excess_j = RHO_C_J_M3_K * LENGTH_M * CROSS_M * CROSS_M * 0.5 * (HOT_K - COLD_K);
    excess_j * (1.0 - held)
}

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
                source: "conduction_pair slab-to-slab declared contact fixture".to_string(),
                license: "internal-test-use".to_string(),
                artifact: None,
            },
        })
        .unwrap();
    InterfaceSystemCard::assemble(
        surface("slab-a", "interface-normal-plus-x"),
        surface("slab-b", "interface-normal-minus-x"),
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

fn contact() -> InterfaceResistance {
    let card = contact_card();
    let contact = InterfaceResistance::from_card(
        "slab-contact",
        &card,
        &QueryPoint::new(),
        SelectionPolicy::SingleClaimOnly,
    )
    .unwrap();
    card.claims().verify_receipt(contact.receipt()).unwrap();
    contact
}

fn material() -> &'static ConductivityModel {
    leak(ConductivityModel::isotropic_declared(CONDUCTIVITY_W_M_K).unwrap())
}

fn source() -> &'static ScalarField {
    leak(ScalarField::Uniform(0.0))
}

/// One slab body and its mesh: `cells` along x, `cross` by `cross` across,
/// its one port on the end at `port_x`.
fn slab_body(
    cx: &Cx<'_>,
    cells: usize,
    cross: usize,
    port: PatchPort<'_>,
    port_x: f64,
) -> (ConductionBody<'static>, &'static ConductionMesh) {
    let (complex, positions) = box_grid([cells, cross, cross], [LENGTH_M, CROSS_M, CROSS_M]);
    let mesh = leak(ConductionMesh::new(complex, positions).unwrap());
    let boundary = leak(
        ThermalBoundaryBuilder::new(mesh)
            .region(
                port.region(),
                |face| on_box_face(face.centroid[0], port_x),
                ThermalBc::adiabatic(),
            )
            .unwrap()
            .adiabatic_remainder()
            .finish()
            .unwrap(),
    );
    let solver = BackwardEuler::uniform(
        cx,
        mesh,
        VolumetricHeatCapacity::declared(RHO_C_J_M3_K).unwrap(),
    )
    .unwrap();
    let problem = ConductionProblem {
        mesh,
        boundary,
        material: material(),
        element_materials: None,
        source: source(),
    };
    (
        ConductionBody::new(solver, problem, &[port], config()).unwrap(),
        mesh,
    )
}

/// The largest spread of a field within one plane x = const: zero when the
/// field is a function of x alone.
fn cross_spread_k(positions: &[[f64; 3]], field: &[f64]) -> f64 {
    let mut planes: Vec<(f64, f64, f64)> = Vec::new();
    for (p, &t) in positions.iter().zip(field) {
        // Box-grid planes share their coordinate bit for bit.
        match planes
            .iter_mut()
            .find(|(x, _, _)| x.to_bits() == p[0].to_bits())
        {
            Some((_, low, high)) => {
                *low = low.min(t);
                *high = high.max(t);
            }
            None => planes.push((p[0], t, t)),
        }
    }
    planes
        .iter()
        .map(|(_, low, high)| high - low)
        .fold(0.0, f64::max)
}

/// What one path reports: the energy that crossed by `t*`, the end field's
/// largest spread across a plane x = const, and the end fields of A and B.
struct Crossed {
    energy_j: f64,
    spread_k: f64,
    end: [Vec<f64>; 2],
}

/// The pair on the kernel: A's `right` end bonded to B's `left` end through
/// the card's conductance. Checks the pair's energy and the ledger close.
fn kernel_pair_crossed(cx: &Cx<'_>, cells: usize, cross_a: usize, cross_b: usize) -> Crossed {
    let (a, mesh_a) = slab_body(cx, cells, cross_a, PatchPort::Mean("right"), LENGTH_M);
    let (b, mesh_b) = slab_body(cx, cells, cross_b, PatchPort::Mean("left"), 0.0);
    let wall_w_k = 1.0
        / contact()
            .term_for_area("slab-contact", a.patch_area_m2(0))
            .unwrap()
            .value_k_per_w();
    let (start_a, start_b) = (
        vec![HOT_K; mesh_a.vertex_count()],
        vec![COLD_K; mesh_b.vertex_count()],
    );
    let start_a_j = a.stored_energy_j(&start_a).unwrap();
    let start_j = start_a_j + b.stored_energy_j(&start_b).unwrap();
    let mut builder = SystemBuilder::new();
    let node_a = builder.add_with_tangent("a", a, start_a).unwrap();
    let node_b = builder.add_with_tangent("b", b, start_b).unwrap();
    let wall = builder
        .add_with_tangent("wall", Conductance::new(wall_w_k).unwrap(), ())
        .unwrap();
    let (wall_a, wall_b) = (
        builder.port(wall, "a").unwrap(),
        builder.port(wall, "b").unwrap(),
    );
    builder
        .bond(builder.port(node_a, "right").unwrap(), wall_a)
        .unwrap();
    builder
        .bond(wall_b, builder.port(node_b, "left").unwrap())
        .unwrap();
    let mut system = builder.seal().unwrap();
    assert!(system.guess(wall_a, &[HOT_K]));
    assert!(system.guess(wall_b, &[COLD_K]));
    let steps = 2 * cells * cells;
    let dt_s = end_time_s() / steps as f64;
    for _ in 0..steps {
        system.step(dt_s, cx).unwrap();
    }
    // The same bodies again, unconnected, read the end fields' energies.
    let (probe_a, _) = slab_body(cx, cells, cross_a, PatchPort::Mean("right"), LENGTH_M);
    let (probe_b, _) = slab_body(cx, cells, cross_b, PatchPort::Mean("left"), 0.0);
    let (end_a, end_b) = (
        system.state::<Vec<f64>>(node_a).unwrap(),
        system.state::<Vec<f64>>(node_b).unwrap(),
    );
    let end_a_j = probe_a.stored_energy_j(end_a).unwrap();
    let end_j = end_a_j + probe_b.stored_energy_j(end_b).unwrap();
    assert!(
        (end_j - start_j).abs() <= 1.0e-9 * start_j,
        "the pair lost {} J",
        start_j - end_j
    );
    assert!(system.ledger().energy_residual_j().abs() <= 1.0e-9 * start_j);
    Crossed {
        energy_j: start_a_j - end_a_j,
        spread_k: cross_spread_k(mesh_a.positions(), end_a)
            .max(cross_spread_k(mesh_b.positions(), end_b)),
        end: [end_a.clone(), end_b.clone()],
    }
}

/// What the trace pair reports besides its crossing: the energy A's own
/// store lost, the entropy the contact generated (in all, and the least in
/// one step), what the kernel's nodal sum over the contact's two ports makes
/// of the same steps, and the kernel's refusal if it refused a step.
struct TracePair {
    crossed: Crossed,
    lost_by_a_j: f64,
    generated_j_k: f64,
    least_generated_j_k: f64,
    nodal_j_k: f64,
    refused: Option<String>,
}

/// The owner's geometry tolerances and budgets for traces that do not match.
fn refinement() -> NonmatchingOptions {
    NonmatchingOptions {
        plane_tolerance_m: 1.0e-12,
        coverage_relative_tolerance: 1.0e-10,
        max_pair_tests: 100_000,
        max_overlap_triangles: 10_000,
    }
}

/// The element of a contact between two traces.
fn wall_of(joint: &TraceContact) -> ConductanceMatrix {
    let points = joint.points.iter().map(|point| InterfacePoint {
        a: point.a,
        b: point.b,
        conductance_w_k: point.conductance_w_k,
    });
    ConductanceMatrix::new(joint.width_a, joint.width_b, points.collect()).unwrap()
}

/// The pair on the kernel, trace to trace: A's `right` end and B's `left`
/// end are field heat ports, joined through the owner's contact from the
/// same card as a `ConductanceMatrix`: the matching contact where the cross
/// grids are equal, the mapped one where they are not. The crossing is the
/// heat the kernel booked out of A's port, step by step, as the owner books
/// its own.
fn trace_pair_crossed(cx: &Cx<'_>, cells: usize, cross_a: usize, cross_b: usize) -> TracePair {
    let shift = [LENGTH_M, 0.0, 0.0];
    let (a, mesh_a) = slab_body(cx, cells, cross_a, PatchPort::Trace("right"), LENGTH_M);
    let (b, mesh_b) = slab_body(cx, cells, cross_b, PatchPort::Trace("left"), 0.0);
    let joint = if cross_a == cross_b {
        a.matching_contact(0, &b, 0, shift, &contact())
    } else {
        a.mapped_contact(cx, 0, &b, 0, shift, &contact(), refinement())
    }
    .unwrap();
    let (points_a, points_b) = (joint.width_a, joint.width_b);
    let wall = wall_of(&joint);
    let wall_ports = wall.ports().to_vec();
    let (start_a, start_b) = (
        vec![HOT_K; mesh_a.vertex_count()],
        vec![COLD_K; mesh_b.vertex_count()],
    );
    let start_a_j = a.stored_energy_j(&start_a).unwrap();
    let start_j = start_a_j + b.stored_energy_j(&start_b).unwrap();
    let mut builder = SystemBuilder::new();
    let node_a = builder.add_with_tangent("a", a, start_a).unwrap();
    let node_b = builder.add_with_tangent("b", b, start_b).unwrap();
    let wall = builder.add_with_tangent("wall", wall, ()).unwrap();
    let (wall_a, wall_b) = (
        builder.port(wall, "a").unwrap(),
        builder.port(wall, "b").unwrap(),
    );
    builder
        .bond(builder.port(node_a, "right").unwrap(), wall_a)
        .unwrap();
    builder
        .bond(wall_b, builder.port(node_b, "left").unwrap())
        .unwrap();
    let mut system = builder.seal().unwrap();
    assert!(system.guess(wall_a, &vec![HOT_K; points_a]));
    assert!(system.guess(wall_b, &vec![COLD_K; points_b]));
    let steps = 2 * cells * cells;
    let dt_s = end_time_s() / steps as f64;
    let (mut least_generated_j_k, mut refused) = (f64::INFINITY, None);
    let (mut generated_j_k, mut nodal_j_k) = (0.0, 0.0);
    for step in 0..steps {
        if let Err(error) = system.step(dt_s, cx) {
            refused = Some(format!("step {step}: {error:?}"));
            break;
        }
        let account = &system.ledger().nodes[2];
        assert_eq!(account.name, "wall");
        let total = account.entropy_generated_j_per_k.value();
        least_generated_j_k = least_generated_j_k.min(total - generated_j_k);
        generated_j_k = total;
        // Minus the entropy the kernel's nodal rule booked into the contact.
        for (spec, port) in wall_ports.iter().zip([wall_a, wall_b]) {
            let (kelvin, loads) = (system.last_input(port), system.last_output(port));
            nodal_j_k -= spec.entropy_in_j_per_k(dt_s, kelvin, loads);
        }
    }
    let (probe_a, _) = slab_body(cx, cells, cross_a, PatchPort::Trace("right"), LENGTH_M);
    let (probe_b, _) = slab_body(cx, cells, cross_b, PatchPort::Trace("left"), 0.0);
    let (end_a, end_b) = (
        system.state::<Vec<f64>>(node_a).unwrap(),
        system.state::<Vec<f64>>(node_b).unwrap(),
    );
    let end_a_j = probe_a.stored_energy_j(end_a).unwrap();
    let end_j = end_a_j + probe_b.stored_energy_j(end_b).unwrap();
    assert!(
        (end_j - start_j).abs() <= 1.0e-9 * start_j,
        "the pair lost {} J",
        start_j - end_j
    );
    assert!(system.ledger().energy_residual_j().abs() <= 1.0e-9 * start_j);
    let account = &system.ledger().nodes[0];
    assert_eq!(account.name, "a");
    TracePair {
        crossed: Crossed {
            energy_j: -account.port_energy_in_j.value(),
            spread_k: cross_spread_k(mesh_a.positions(), end_a)
                .max(cross_spread_k(mesh_b.positions(), end_b)),
            end: [end_a.clone(), end_b.clone()],
        },
        lost_by_a_j: start_a_j - end_a_j,
        generated_j_k,
        least_generated_j_k,
        nodal_j_k,
        refused,
    }
}

/// The owner's own contact: both slabs in one mesh with unshared traces at
/// x = L, bound by the same card. Equal cross grids are its matching-P1
/// contact on coincident faces; unequal ones its non-matching contact on
/// the common refinement of the two traces.
fn owner_contact_crossed(cx: &Cx<'_>, cells: usize, cross_a: usize, cross_b: usize) -> Crossed {
    let (slab_a, mut positions) = box_grid([cells, cross_a, cross_a], [LENGTH_M, CROSS_M, CROSS_M]);
    let (slab_b, positions_b) = box_grid([cells, cross_b, cross_b], [LENGTH_M, CROSS_M, CROSS_M]);
    let per_slab = positions.len();
    let offset = u32::try_from(per_slab).unwrap();
    let mut tets = slab_a.tets;
    tets.extend(slab_b.tets.iter().map(|tet| tet.map(|v| v + offset)));
    positions.extend(positions_b.iter().map(|&[x, y, z]| [x + LENGTH_M, y, z]));
    let mesh =
        ConductionMesh::new(TetComplex::from_tets(positions.len(), tets), positions).unwrap();
    let boundary: ThermalBoundary = ThermalBoundaryBuilder::new(&mesh)
        .adiabatic_remainder()
        .finish()
        .unwrap();
    let interfaces = if cross_a == cross_b {
        let pairs = ThermalInterfaces::coincident_face_pairs(&mesh)
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
        let surface = InterfaceSurface::new("slab-contact", pairs, contact()).unwrap();
        ThermalInterfaces::new(&mesh, &boundary, vec![surface]).unwrap()
    } else {
        // The faces on the contact plane: A's look along +x, B's along -x.
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
            contact(),
            refinement(),
        )
        .unwrap();
        ThermalInterfaces::with_nonmatching(cx, &mesh, &boundary, Vec::new(), vec![surface])
            .unwrap()
    };
    let solver = BackwardEuler::uniform(
        cx,
        &mesh,
        VolumetricHeatCapacity::declared(RHO_C_J_M3_K).unwrap(),
    )
    .unwrap();
    let problem = ConductionProblem {
        mesh: &mesh,
        boundary: &boundary,
        material: material(),
        element_materials: None,
        source: source(),
    };
    let mut field: Vec<f64> = (0..mesh.vertex_count())
        .map(|v| if v < per_slab { HOT_K } else { COLD_K })
        .collect();
    let steps = 2 * cells * cells;
    let dt_s = end_time_s() / steps as f64;
    let mut crossed_j = 0.0;
    for _ in 0..steps {
        let step = solver
            .advance(cx, problem, Some(&interfaces), &field, dt_s, config())
            .unwrap();
        assert_eq!(step.contact_fluxes.len(), 1);
        crossed_j += dt_s * step.contact_fluxes[0].heat_rate_a_to_b_w;
        field = step.temperature;
    }
    // Each slab on its own: the contact plane holds both traces, a jump apart.
    let (positions_a, positions_b) = mesh.positions().split_at(per_slab);
    let (field_a, field_b) = field.split_at(per_slab);
    Crossed {
        energy_j: crossed_j,
        spread_k: cross_spread_k(positions_a, field_a).max(cross_spread_k(positions_b, field_b)),
        end: [field_a.to_vec(), field_b.to_vec()],
    }
}

/// The relative error of each rung against the closed form, the fitted order,
/// and the slope of the last two rungs.
fn ladder(errors: &[f64]) -> (f64, f64) {
    let hs: Vec<f64> = RUNGS.iter().map(|&n| LENGTH_M / n as f64).collect();
    let fit = fit_order(&RefinementLadder::new(hs, errors.to_vec()).unwrap());
    let last = (errors[2] / errors[3]).log2();
    (fit.observed, last)
}

#[test]
fn two_bodies_through_a_card_contact_converge_at_second_order_with_the_owner_contact() {
    with_cx(|cx| {
        let exact_j = closed_form_crossed_j();
        let excess_j = RHO_C_J_M3_K * LENGTH_M * CROSS_M * CROSS_M * 0.5 * (HOT_K - COLD_K);
        eprintln!("closed form: {exact_j:.15e} J crossed of the {excess_j:.15e} J excess");
        let mut paths = [Vec::new(), Vec::new(), Vec::new()];
        for &cells in &RUNGS {
            let matching = kernel_pair_crossed(cx, cells, 2, 2);
            let nonmatching = kernel_pair_crossed(cx, cells, 2, 3);
            let owner = owner_contact_crossed(cx, cells, 2, 2);
            eprintln!(
                "n {cells}: crossed (a) {:.15e} (b) {:.15e} (c) {:.15e} J; (a) - (c) {:e} J; cross-plane spread (a) {:e} (b) {:e} (c) {:e} K",
                matching.energy_j,
                nonmatching.energy_j,
                owner.energy_j,
                matching.energy_j - owner.energy_j,
                matching.spread_k,
                nonmatching.spread_k,
                owner.spread_k
            );
            for (errors, crossed) in paths.iter_mut().zip([matching, nonmatching, owner]) {
                errors.push((crossed.energy_j - exact_j).abs() / excess_j);
            }
        }
        let mut missed = Vec::new();
        for (name, errors) in [
            "(a) matching pair",
            "(b) non-matching pair",
            "(c) owner contact",
        ]
        .iter()
        .zip(&paths)
        {
            let (fitted, last) = ladder(errors);
            eprintln!(
                "{name}: relative errors {errors:?}; fitted order {fitted:.3}, last two rungs {last:.3}"
            );
            if !(last >= CLAIMED_ORDER) {
                missed.push(format!("{name} {last:.3}"));
            }
        }
        assert!(missed.is_empty(), "below order {CLAIMED_ORDER}: {missed:?}");
    });
}

/// The pair trace to trace against the owner's own contact on one mesh, rung
/// by rung, every rung reported: the rungs that missed a claim.
fn against_the_owner(cx: &Cx<'_>, path: &str, cross_b: usize) -> Vec<String> {
    let mut missed = Vec::new();
    for &cells in &RUNGS {
        let pair = trace_pair_crossed(cx, cells, 2, cross_b);
        let owner = owner_contact_crossed(cx, cells, 2, cross_b);
        let gap_j = pair.crossed.energy_j - owner.energy_j;
        let mut gap_k = 0.0_f64;
        for (traced, owned) in pair.crossed.end.iter().zip(&owner.end) {
            for (t, o) in traced.iter().zip(owned) {
                gap_k = gap_k.max((t - o).abs());
            }
        }
        eprintln!(
            "n {cells}: crossed {path} {:.15e} owner {:.15e} J; {path} - owner {gap_j:e} J; A's store lost {:.15e} J; largest vertex difference {gap_k:e} K; cross-plane spread {path} {:e} owner {:e} K; entropy the contact generated {:e} J/K, least in a step {:e} J/K; the nodal sum over its ports {:e} J/K, {:e} J/K from the generation; refused: {:?}",
            pair.crossed.energy_j,
            owner.energy_j,
            pair.lost_by_a_j,
            pair.crossed.spread_k,
            owner.spread_k,
            pair.generated_j_k,
            pair.least_generated_j_k,
            pair.nodal_j_k,
            pair.nodal_j_k - pair.generated_j_k,
            pair.refused
        );
        if pair.refused.is_some() || !(gap_j.abs() <= CLAIMED_GAP_J) || !(gap_k <= CLAIMED_GAP_K) {
            missed.push(format!(
                "n {cells}: {gap_j:e} J, {gap_k:e} K, {:?}",
                pair.refused
            ));
        }
    }
    missed
}

/// Trace to trace through the owner's contact, the pair is the owner's own
/// contact problem on one mesh: the same unknowns under the same equations,
/// solved in another order.
///
/// The claims, declared before the first run: on every rung the energy that
/// crossed agrees with path (c)'s within [`CLAIMED_GAP_J`], every vertex's
/// end temperature agrees with (c)'s within [`CLAIMED_GAP_K`], and the kernel
/// refuses no step, its entropy check included. No order is claimed: the gap
/// is the two solves' tolerances on every rung. Through the patch means,
/// path (a), the crossing differed from (c)'s by 4.6e-2 J on the coarsest
/// rung and 2.9e-5 J on the finest. The entropy the contact generates and
/// the kernel's nodal sum over its two ports are reported with no claim on
/// their difference. Every rung is reported before any is asserted.
#[test]
fn two_bodies_trace_to_trace_through_the_card_are_the_owner_contact_on_one_mesh() {
    with_cx(|cx| {
        let missed = against_the_owner(cx, "(d)", 2);
        assert!(missed.is_empty(), "path (d) against path (c): {missed:?}");
    });
}

/// The same where the cross grids do not match, 2 x 2 against 3 x 3: path
/// (e), through the mapped contact, against the owner's non-matching contact
/// on one mesh at the same steps. The same claims, declared before the first
/// run.
#[test]
fn two_bodies_whose_traces_do_not_match_are_the_owner_nonmatching_contact_on_one_mesh() {
    with_cx(|cx| {
        let missed = against_the_owner(cx, "(e)", 3);
        assert!(
            missed.is_empty(),
            "path (e) against the owner's non-matching contact: {missed:?}"
        );
    });
}

#[test]
fn a_matching_contact_is_the_owner_face_conductance_and_other_traces_are_refused() {
    with_cx(|cx| {
        let shift = [LENGTH_M, 0.0, 0.0];
        let (a, mesh_a) = slab_body(cx, 2, 2, PatchPort::Trace("right"), LENGTH_M);
        let (b, mesh_b) = slab_body(cx, 2, 2, PatchPort::Trace("left"), 0.0);
        let joint = a.matching_contact(0, &b, 0, shift, &contact()).unwrap();
        // Three quadrature points on each of the patch's eight faces.
        assert_eq!(
            (joint.width_a, joint.width_b, joint.points.len()),
            (9, 9, 24)
        );
        let matrix = wall_of(&joint).matrix();
        // The owner's matching-P1 assembly by hand: on every face of the
        // patch `A / (6 R'')` between a vertex and itself and `A / (12 R'')`
        // between two of its vertices, with the sign of the two sides.
        let (trace_a, trace_b) = (a.trace_vertices(0).unwrap(), b.trace_vertices(0).unwrap());
        let across: Vec<usize> = trace_a
            .iter()
            .map(|&vertex| {
                let [_, y, z] = mesh_a.positions()[vertex];
                let same = |w: &&usize| {
                    let [_, v, w] = mesh_b.positions()[**w];
                    v.to_bits() == y.to_bits() && w.to_bits() == z.to_bits()
                };
                trace_b.iter().position(|w| same(&w)).unwrap()
            })
            .collect();
        let mut owner = vec![0.0; 18 * 18];
        for face in mesh_a.boundary() {
            if !on_box_face(face.centroid[0], LENGTH_M) {
                continue;
            }
            let points = face
                .vertices
                .map(|v| trace_a.binary_search(&(v as usize)).unwrap());
            for (i, &row) in points.iter().enumerate() {
                for (j, &column) in points.iter().enumerate() {
                    let mass = if i == j {
                        face.area / 6.0
                    } else {
                        face.area / 12.0
                    };
                    let value = (1.0 / CONTACT_M2_K_W) * mass;
                    let (far_row, far_column) = (9 + across[row], 9 + across[column]);
                    owner[row * 18 + column] += value;
                    owner[far_row * 18 + far_column] += value;
                    owner[row * 18 + far_column] -= value;
                    owner[far_row * 18 + column] -= value;
                }
            }
        }
        // Either entry is a sum of like-signed terms: at most 24 of
        // `G w_i w_j`, three roundings each, against at most 8 of the
        // owner's, two roundings each. 40 roundoff units of the entry hold
        // both sums, and where the owner has no entry there is none.
        let mut worst = 0.0_f64;
        for (entry, expected) in matrix.iter().zip(&owner) {
            assert!(
                (entry - expected).abs() <= 40.0 * f64::EPSILON * expected.abs(),
                "{entry} against {expected} W/K"
            );
            if expected.abs() > 0.0 {
                worst = worst.max(((entry - expected) / expected).abs());
            }
        }
        eprintln!(
            "matching contact against the owner's face conductances: largest relative difference {worst:e} over 324 entries"
        );
        let refused = |other: &ConductionBody<'_>, offset_m: [f64; 3], what: &str| match a
            .matching_contact(0, other, 0, offset_m, &contact())
        {
            Err(ConductionError::Config { what: why, .. }) => {
                assert!(why.contains(what), "{why:?} does not say {what:?}");
            }
            other => panic!("{what}: {other:?}"),
        };
        let (mean, _) = slab_body(cx, 2, 2, PatchPort::Mean("left"), 0.0);
        refused(&mean, shift, "two trace ports");
        // B left where its own mesh has it: its trace is a slab's length away.
        refused(&b, [0.0; 3], "no point of the second trace");
        // A finer cross grid shares the corners only.
        let (fine, _) = slab_body(cx, 2, 3, PatchPort::Trace("left"), 0.0);
        refused(&fine, shift, "do not match");
        // B's own right end coincides with A's, facing the same way.
        let (overlapping, _) = slab_body(cx, 2, 2, PatchPort::Trace("right"), LENGTH_M);
        refused(&overlapping, [0.0; 3], "do not face each other");
    });
}

#[test]
fn a_mapped_contact_covers_the_patch_is_its_own_transpose_and_keeps_a_uniform_field() {
    with_cx(|cx| {
        let shift = [LENGTH_M, 0.0, 0.0];
        let (a, _) = slab_body(cx, 2, 2, PatchPort::Trace("right"), LENGTH_M);
        let (b, _) = slab_body(cx, 2, 3, PatchPort::Trace("left"), 0.0);
        let joint = a
            .mapped_contact(cx, 0, &b, 0, shift, &contact(), refinement())
            .unwrap();
        assert_eq!((joint.width_a, joint.width_b), (9, 16));
        // The quadrature covers the patch once: its conductances sum to
        // `A / R''`, within the owner's coverage tolerance on each face.
        let total: f64 = joint.points.iter().map(|point| point.conductance_w_k).sum();
        let patch = CROSS_M * CROSS_M / CONTACT_M2_K_W;
        assert!((total - patch).abs() <= 1.0e-9 * patch, "{total} W/K");
        eprintln!(
            "mapped contact: {} quadrature points, {total} W/K of {patch} W/K",
            joint.points.len()
        );
        let wall = wall_of(&joint);
        let matrix = wall.matrix();
        // K is symmetric and a uniform temperature moves no heat, each to
        // the rounding of the sums of an entry's or a row's terms.
        for row in 0..25 {
            let entries = &matrix[row * 25..row * 25 + 25];
            let size: f64 = entries.iter().map(|entry| entry.abs()).sum();
            assert!(entries.iter().sum::<f64>().abs() <= 64.0 * f64::EPSILON * size);
            for column in 0..row {
                let mirror = matrix[column * 25 + row];
                assert!((entries[column] - mirror).abs() <= 64.0 * f64::EPSILON * mirror.abs());
            }
        }
        // The kit, and the tangent with its transpose (the dot-product
        // test), at a field with no pattern of the grids.
        let kelvin: Vec<f64> = (0..25)
            .map(|i: i32| 300.0 + 13.0 * f64::from((i * 7) % 11))
            .collect();
        let case = Case {
            state: &(),
            stale: &(),
            dt_s: 0.5,
            inputs: &kelvin,
        };
        check_component(&wall, &case, cx).unwrap_or_else(|failure| panic!("{failure}"));
        check_port_tangent(&wall, &case, 1.0e-8, 1.0e-9, cx)
            .unwrap_or_else(|failure| panic!("tangent: {failure}"));
        // A mean port has no trace, and traces a slab apart share no plane.
        let (mean, _) = slab_body(cx, 2, 3, PatchPort::Mean("left"), 0.0);
        match a.mapped_contact(cx, 0, &mean, 0, shift, &contact(), refinement()) {
            Err(ConductionError::Config { what, .. }) => assert!(what.contains("two trace ports")),
            other => panic!("a mean port: {other:?}"),
        }
        match a.mapped_contact(cx, 0, &b, 0, [0.0; 3], &contact(), refinement()) {
            Err(ConductionError::Interface { what, .. }) => eprintln!("traces apart: {what}"),
            other => panic!("traces a slab apart: {other:?}"),
        }
    });
}
