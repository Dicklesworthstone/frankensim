//! A card's contact resistance at the trace temperatures, on the composition kernel.
//!
//! `fs_conduction::component::TraceResistanceLaw` is an interface card's `R''` at the
//! temperatures of two traces, as the factor a `ConductanceMatrix` of fs-system follows. Two
//! slabs joined trace to trace through it on the kernel are compared here with the owner's own
//! temperature contact, both slabs in one mesh stepped by `advance_capacity_newton`: its
//! pointwise contact where the traces match and its non-matching contact where they do not.
//!
//! The cards are authored for these tests and are no material's. The main one lets `R''` fall
//! from 1.4e-3 m² K/W at 280 K through 8e-4 at 345 K to 5e-4 at 430 K, over the mean of the two
//! trace temperatures, and has no answer outside 280 K to 430 K.
//!
//! `cargo test -p fs-conduction --features component --test trace_resistance_law`.
#![cfg(feature = "component")]

use std::sync::Arc;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_conduction::component::{ConductionBody, PatchPort, TraceContact, TraceResistanceLaw};
use fs_conduction::fixtures::{box_grid, on_box_face};
use fs_conduction::interface::NonmatchingOptions;
use fs_conduction::interface::temperature::{
    ContactIterationConfig, TemperatureContact, TemperatureContactCoordinates, TemperatureContacts,
};
use fs_conduction::transient::VolumetricHeatCapacity;
use fs_conduction::transient::backward_euler::{BackwardEuler, NonlinearStepConfig, StepConfig};
use fs_conduction::{
    AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS, AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY, BoundaryFace,
    ConductionError, ConductionMesh, ConductionProblem, ConductivityModel, InterfaceFacePair,
    LinearConfig, ScalarField, ThermalBc, ThermalBoundaryBuilder, ThermalInterfaces,
};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, SurfaceSpec,
    SystemContext, UncertaintyModel,
};
use fs_ports::component::Component;
use fs_rep_mesh::TetComplex;
use fs_system::SystemBuilder;
use fs_system::elements::{ConductanceLaw, ConductanceMatrix, GuardRow, InterfacePoint};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let key = StreamKey {
            seed: 67,
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

/// Each slab: 30 mm along x in 4 cells, 20 x 20 mm across. Slab B stores half as much per
/// kelvin as slab A, so the mean of the two traces moves as heat crosses.
const LENGTH_M: f64 = 0.03;
const CROSS_M: f64 = 0.02;
const CELLS: usize = 4;
const CONDUCTIVITY_W_M_K: f64 = 15.0;
const RHO_C_A_J_M3_K: f64 = 3.6e6;
const RHO_C_B_J_M3_K: f64 = 1.8e6;
/// Slab A starts at 380 K where y = 0 and 30 K hotter where y = 20 mm; slab B at 300 K. The
/// means at the contact then run from 340 K to 355 K, on both sides of the card's knot.
const HOT_K: f64 = 380.0;
const TILT_K: f64 = 30.0;
const COLD_K: f64 = 300.0;
const STEPS: usize = 12;
const STEP_S: f64 = 1.8;
/// The card's `R''` [m² K/W] over the mean trace temperature [K].
const KNOTS: [(f64, f64); 3] = [(280.0, 1.4e-3), (345.0, 8.0e-4), (430.0, 5.0e-4)];
const REFERENCE_K: f64 = 320.0;
const NAME: &str = "slab-contact";
/// Kernel against owner: the most any vertex may differ after any step, and the most the
/// energy that crossed may.
const CLAIMED_GAP_K: f64 = 1.0e-7;
const CLAIMED_GAP_J: f64 = 1.0e-7;
/// The least the law must move some vertex from the same contact at its constant reference
/// resistance, so that a law the element ignored could not pass.
const CLAIMED_MOVE_K: f64 = 0.1;

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

fn newton() -> NonlinearStepConfig {
    NonlinearStepConfig {
        residual_rtol: 1.0e-12,
        residual_atol_j: 1.0e-11,
        ..Default::default()
    }
}

fn coupling() -> ContactIterationConfig {
    ContactIterationConfig {
        max_iterations: 64,
        relaxation: 0.5,
        residual_tolerance_j: 1.0e-8,
    }
}

fn refinement() -> NonmatchingOptions {
    NonmatchingOptions {
        plane_tolerance_m: 1.0e-12,
        coverage_relative_tolerance: 1.0e-10,
        max_pair_tests: 100_000,
        max_overlap_triangles: 10_000,
    }
}

fn leak<T>(value: T) -> &'static T {
    Box::leak(Box::new(value))
}

fn kelvin() -> fs_qty::QuantitySpec {
    use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
    fs_qty::QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}

/// A card with the one resistance claim `value`.
fn card_of(
    value: PropertyValue,
    interpolation: InterpolationPolicy,
    validity: ValidityDomain,
) -> InterfaceSystemCard {
    let surface = |chemistry: &str| SurfaceSpec {
        material: MaterialStateId {
            chemistry: chemistry.to_string(),
            phase: "solid".to_string(),
            process: "as-fixtured".to_string(),
            revision: 0,
        },
        texture_frame: "interface normal".to_string(),
    };
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(
                AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY,
                AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS,
            ),
            value,
            validity,
            uncertainty: UncertaintyModel::Unstated,
            interpolation,
            observations: Vec::new(),
            provenance: Provenance {
                source: "authored for trace_resistance_law, not measurements".to_string(),
                license: "internal-test-use".to_string(),
                artifact: None,
            },
        })
        .unwrap();
    InterfaceSystemCard::assemble(
        surface("slab-a"),
        surface("slab-b"),
        SystemContext {
            medium: "dry".to_string(),
            third_body: None,
            environment: "vacuum".to_string(),
            history: "unaged".to_string(),
        },
        claims,
        Vec::new(),
    )
    .unwrap()
}

/// A linear curve of the resistance over `abscissa`.
fn curve(abscissa: &str, abscissa_dims: fs_qty::Dims, knots: &[(f64, f64)]) -> PropertyValue {
    PropertyValue::Curve {
        abscissa: abscissa.to_string(),
        abscissa_dims,
        knots: knots.to_vec(),
        dims: AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS,
    }
}

/// The main card: [`KNOTS`] over the axis `T`.
fn card() -> InterfaceSystemCard {
    card_of(
        curve("T", fs_conduction::TEMPERATURE_DIMS, &KNOTS),
        InterpolationPolicy::LinearInside,
        ValidityDomain::unconstrained().with("T", KNOTS[0].0, KNOTS[2].0),
    )
}

fn reference() -> QueryPoint {
    QueryPoint::new().with("T", REFERENCE_K).unwrap()
}

fn mean() -> TemperatureContactCoordinates {
    TemperatureContactCoordinates::Mean("T".to_string())
}

/// The main card's law over the mean of the two traces.
fn law() -> TraceResistanceLaw {
    TraceResistanceLaw::bind(
        NAME,
        &card(),
        &reference(),
        mean(),
        SelectionPolicy::SingleClaimOnly,
    )
    .unwrap()
}

/// [`KNOTS`] by hand: the resistance at `kelvin` and the slope of segment `segment`.
fn by_hand(kelvin: f64, segment: usize) -> (f64, f64) {
    let ((x0, y0), (x1, y1)) = (KNOTS[segment], KNOTS[segment + 1]);
    let slope = (y1 - y0) / (x1 - x0);
    (y0 + (kelvin - x0) * slope, slope)
}

fn close(found: f64, expected: f64, relative: f64, what: &str) {
    assert!(
        (found - expected).abs() <= relative * expected.abs(),
        "{what}: {found:e} against {expected:e}"
    );
}

/// The law as the element's hook: the owner's refusal becomes the element's reason.
#[derive(Debug)]
struct CardLaw(TraceResistanceLaw);

impl ConductanceLaw for CardLaw {
    fn factor(&self, kelvin_a: f64, kelvin_b: f64) -> Result<(f64, [f64; 2]), String> {
        self.0
            .factor(kelvin_a, kelvin_b)
            .map_err(|error| error.to_string())
    }

    fn identity_words(&self, words: &mut Vec<u64>) {
        self.0.identity_words(words);
    }
}

/// How the two traces are joined on the kernel.
#[derive(Debug, Clone, Copy, PartialEq)]
enum Wall {
    /// The contact at the law's reference resistance, with no law.
    Constant,
    /// The contact following the law, guarded at its facing pairs.
    Law,
    /// The contact following the law with no guard rows.
    Unguarded,
}

fn wall_of(joint: &TraceContact, law: &TraceResistanceLaw, wall: Wall) -> ConductanceMatrix {
    let points = joint.points.iter().map(|point| InterfacePoint {
        a: point.a,
        b: point.b,
        conductance_w_k: point.conductance_w_k,
    });
    let element = ConductanceMatrix::new(joint.width_a, joint.width_b, points.collect()).unwrap();
    if wall == Wall::Constant {
        return element;
    }
    let element = element.following(Arc::new(CardLaw(law.clone())));
    if wall == Wall::Unguarded {
        return element;
    }
    let alone = |point: usize| [(point, 1.0), (point, 0.0), (point, 0.0)];
    let rows = joint.facing.iter().map(|&[a, b]| GuardRow {
        a: alone(a),
        b: alone(b),
    });
    element.guarded_at(rows.collect()).unwrap()
}

/// One slab body and its mesh, `cross` by `cross` cells across, its one trace port on the end
/// at `port_x`.
fn slab_body(
    cx: &Cx<'_>,
    cross: usize,
    rho_c_j_m3_k: f64,
    port: &'static str,
    port_x: f64,
) -> (ConductionBody<'static>, &'static ConductionMesh) {
    let (complex, positions) = box_grid([CELLS, cross, cross], [LENGTH_M, CROSS_M, CROSS_M]);
    let mesh = leak(ConductionMesh::new(complex, positions).unwrap());
    let boundary = leak(
        ThermalBoundaryBuilder::new(mesh)
            .region(
                port,
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
        VolumetricHeatCapacity::declared(rho_c_j_m3_k).unwrap(),
    )
    .unwrap();
    let problem = ConductionProblem {
        mesh,
        boundary,
        material: leak(ConductivityModel::isotropic_declared(CONDUCTIVITY_W_M_K).unwrap()),
        element_materials: None,
        source: leak(ScalarField::Uniform(0.0)),
    };
    let body = ConductionBody::new(solver, problem, &[PatchPort::Trace(port)], config()).unwrap();
    (body, mesh)
}

/// Slab A's start: hotter towards y = 20 mm.
fn tilted(position: [f64; 3]) -> f64 {
    HOT_K + TILT_K * position[1] / CROSS_M
}

/// What a run of the pair leaves: both fields after every step it took, both fields as it
/// ended, the energy the kernel booked out of slab A's port [J], and the first refusal.
struct Run {
    fields: Vec<[Vec<f64>; 2]>,
    kept: [Vec<f64>; 2],
    crossed_j: f64,
    refused: Option<String>,
}

/// The pair on the kernel, trace to trace: the matching contact where slab B has slab A's 2 x 2
/// cross grid, the mapped one where it does not.
fn kernel_run(
    cx: &Cx<'_>,
    cross_b: usize,
    wall: Wall,
    start_a: fn([f64; 3]) -> f64,
    steps: usize,
) -> Run {
    let shift = [LENGTH_M, 0.0, 0.0];
    let (a, mesh_a) = slab_body(cx, 2, RHO_C_A_J_M3_K, "right", LENGTH_M);
    let (b, mesh_b) = slab_body(cx, cross_b, RHO_C_B_J_M3_K, "left", 0.0);
    let law = law();
    let joint = if cross_b == 2 {
        a.matching_contact(0, &b, 0, shift, law.reference_resistance())
    } else {
        a.mapped_contact(
            cx,
            0,
            &b,
            0,
            shift,
            law.reference_resistance(),
            refinement(),
        )
    }
    .unwrap();
    assert_eq!(joint.facing.len(), if cross_b == 2 { 9 } else { 0 });
    let element = wall_of(&joint, &law, wall);
    let field_a: Vec<f64> = mesh_a.positions().iter().map(|&p| start_a(p)).collect();
    let field_b = vec![COLD_K; mesh_b.vertex_count()];
    let trace = |body: &ConductionBody<'_>, field: &[f64]| -> Vec<f64> {
        let vertices = body.trace_vertices(0).unwrap();
        vertices.iter().map(|&vertex| field[vertex]).collect()
    };
    let (guess_a, guess_b) = (trace(&a, &field_a), trace(&b, &field_b));
    let mut builder = SystemBuilder::new();
    let node_a = builder.add_with_tangent("a", a, field_a).unwrap();
    let node_b = builder.add_with_tangent("b", b, field_b).unwrap();
    let element = builder.add_with_tangent("wall", element, ()).unwrap();
    let (wall_a, wall_b) = (
        builder.port(element, "a").unwrap(),
        builder.port(element, "b").unwrap(),
    );
    builder
        .bond(builder.port(node_a, "right").unwrap(), wall_a)
        .unwrap();
    builder
        .bond(wall_b, builder.port(node_b, "left").unwrap())
        .unwrap();
    let mut system = builder.seal().unwrap();
    assert!(system.guess(wall_a, &guess_a));
    assert!(system.guess(wall_b, &guess_b));
    let (mut fields, mut refused) = (Vec::new(), None);
    for step in 0..steps {
        if let Err(error) = system.step(STEP_S, cx) {
            refused = Some(format!("step {step}: {error:?}"));
            break;
        }
        fields.push([node_a, node_b].map(|node| system.state::<Vec<f64>>(node).unwrap().clone()));
    }
    let account = &system.ledger().nodes[0];
    assert_eq!(account.name, "a");
    Run {
        fields,
        kept: [node_a, node_b].map(|node| system.state::<Vec<f64>>(node).unwrap().clone()),
        crossed_j: -account.port_energy_in_j.value(),
        refused,
    }
}

/// Both slabs in one mesh with unshared traces at x = 30 mm: the mesh, its start field, and
/// how many vertices and elements are slab A's (they come first).
fn one_mesh(
    cross_b: usize,
    start_a: fn([f64; 3]) -> f64,
) -> (ConductionMesh, Vec<f64>, [usize; 2]) {
    let (slab_a, mut positions) = box_grid([CELLS, 2, 2], [LENGTH_M, CROSS_M, CROSS_M]);
    let (slab_b, positions_b) = box_grid([CELLS, cross_b, cross_b], [LENGTH_M, CROSS_M, CROSS_M]);
    let of_a = [positions.len(), slab_a.tets.len()];
    let offset = u32::try_from(positions.len()).unwrap();
    let mut field: Vec<f64> = positions.iter().map(|&p| start_a(p)).collect();
    field.extend(positions_b.iter().map(|_| COLD_K));
    let mut tets = slab_a.tets;
    tets.extend(slab_b.tets.iter().map(|tet| tet.map(|v| v + offset)));
    positions.extend(positions_b.iter().map(|&[x, y, z]| [x + LENGTH_M, y, z]));
    let mesh =
        ConductionMesh::new(TetComplex::from_tets(positions.len(), tets), positions).unwrap();
    (mesh, field, of_a)
}

/// The owner's contact law on that mesh, from the same card at the same reference point: its
/// pointwise contact on coincident faces where the cross grids are equal, its non-matching
/// contact on the two sides of the plane where they are not.
fn owner_contact<'c>(
    mesh: &ConductionMesh,
    card: &'c InterfaceSystemCard,
    cross_b: usize,
) -> TemperatureContact<'c> {
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
    let (point, selection) = (reference(), SelectionPolicy::SingleClaimOnly);
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
        TemperatureContact::new_pointwise(NAME, pairs.collect(), card, &point, mean(), selection)
    } else {
        let (a, b) = (side(1.0), side(-1.0));
        TemperatureContact::new_nonmatching(
            NAME,
            a,
            b,
            card,
            &point,
            mean(),
            selection,
            refinement(),
        )
    }
    .unwrap()
}

/// The owner's own temperature contact, stepped by `advance_capacity_newton`: both fields after
/// every step, or the owner's refusal of a step.
fn owner_run(
    cx: &Cx<'_>,
    cross_b: usize,
    start_a: fn([f64; 3]) -> f64,
    steps: usize,
) -> Result<Vec<[Vec<f64>; 2]>, ConductionError> {
    let (mesh, mut field, [per_slab, elements_a]) = one_mesh(cross_b, start_a);
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .adiabatic_remainder()
        .finish()
        .unwrap();
    let card = card();
    let contact = owner_contact(&mesh, &card, cross_b);
    let faces = 10_000;
    let (fixed, fixed_nonmatching) = (Vec::new(), Vec::new());
    let set = TemperatureContacts::with_nonmatching(
        cx,
        &mesh,
        &boundary,
        fixed,
        fixed_nonmatching,
        vec![contact],
        faces,
        3 * faces,
    )
    .unwrap();
    let capacities: Vec<VolumetricHeatCapacity> = (0..mesh.element_count())
        .map(|element| {
            let rho_c = if element < elements_a {
                RHO_C_A_J_M3_K
            } else {
                RHO_C_B_J_M3_K
            };
            VolumetricHeatCapacity::declared(rho_c).unwrap()
        })
        .collect();
    let engine = BackwardEuler::per_element(cx, &mesh, &capacities).unwrap();
    let material = ConductivityModel::isotropic_declared(CONDUCTIVITY_W_M_K).unwrap();
    let source = ScalarField::Uniform(0.0);
    let problem = ConductionProblem {
        mesh: &mesh,
        boundary: &boundary,
        material: &material,
        element_materials: None,
        source: &source,
    };
    let mut fields = Vec::new();
    for _ in 0..steps {
        let step = set.advance_capacity_newton(
            cx,
            &engine,
            problem,
            &field,
            STEP_S,
            config(),
            newton(),
            coupling(),
            faces,
        )?;
        field = step.conduction.temperature;
        let (field_a, field_b) = field.split_at(per_slab);
        fields.push([field_a.to_vec(), field_b.to_vec()]);
    }
    Ok(fields)
}

/// The largest difference between two runs' fields over every step and vertex [K].
fn gap_k(one: &[[Vec<f64>; 2]], other: &[[Vec<f64>; 2]]) -> f64 {
    assert_eq!(one.len(), other.len());
    let mut gap = 0.0_f64;
    for (ours, theirs) in one.iter().zip(other) {
        for (field, owned) in ours.iter().zip(theirs) {
            assert_eq!(field.len(), owned.len());
            for (t, o) in field.iter().zip(owned) {
                gap = gap.max((t - o).abs());
            }
        }
    }
    gap
}

/// The pair through the law on the kernel against the owner's temperature contact on one mesh,
/// after every step, and against the same pair at the constant reference resistance.
fn against_the_owner(cx: &Cx<'_>, cross_b: usize) {
    let pair = kernel_run(cx, cross_b, Wall::Law, tilted, STEPS);
    assert_eq!(pair.refused, None);
    let owner = owner_run(cx, cross_b, tilted, STEPS).unwrap();
    let constant = kernel_run(cx, cross_b, Wall::Constant, tilted, STEPS);
    assert_eq!(constant.refused, None);
    // What slab A lost in the owner's run, by slab A's own stored energy.
    let (probe, mesh_a) = slab_body(cx, 2, RHO_C_A_J_M3_K, "right", LENGTH_M);
    let start: Vec<f64> = mesh_a.positions().iter().map(|&p| tilted(p)).collect();
    let owner_crossed_j = probe.stored_energy_j(&start).unwrap()
        - probe.stored_energy_j(&owner[STEPS - 1][0]).unwrap();
    let gap_j = pair.crossed_j - owner_crossed_j;
    let gap = gap_k(&pair.fields, &owner);
    let moved = gap_k(&pair.fields, &constant.fields);
    eprintln!(
        "cross grid of B {cross_b} x {cross_b}: {STEPS} steps; crossed kernel {:.15e} owner {owner_crossed_j:.15e} J, kernel - owner {gap_j:e} J; largest vertex difference over all steps {gap:e} K; the law moves a vertex by up to {moved:e} K from the contact at its reference resistance, which crossed {:.15e} J",
        pair.crossed_j, constant.crossed_j
    );
    assert!(gap <= CLAIMED_GAP_K, "{gap:e} K");
    assert!(gap_j.abs() <= CLAIMED_GAP_J, "{gap_j:e} J");
    assert!(moved >= CLAIMED_MOVE_K, "{moved:e} K");
}

/// Claims, declared before the first run: after every one of 12 steps every vertex of both
/// slabs is within [`CLAIMED_GAP_K`] of the owner's pointwise contact on one mesh, the energy
/// that crossed is within [`CLAIMED_GAP_J`] of what slab A lost there, the kernel refuses no
/// step, and the law moves some vertex by at least [`CLAIMED_MOVE_K`] from the constant contact.
#[test]
fn two_slabs_through_the_law_are_the_owner_pointwise_temperature_contact_on_one_mesh() {
    with_cx(|cx| against_the_owner(cx, 2));
}

/// The same claims where the cross grids do not match, 2 x 2 against 3 x 3: the mapped contact
/// through the law against the owner's non-matching temperature contact on one mesh.
#[test]
fn two_slabs_whose_traces_do_not_match_are_the_owner_nonmatching_temperature_contact() {
    with_cx(|cx| against_the_owner(cx, 3));
}

/// Slab A's start with the trace corner at the origin of the contact plane at 600 K and every
/// other vertex at 330 K.
fn hot_corner(position: [f64; 3]) -> f64 {
    let at_corner = [LENGTH_M, 0.0, 0.0];
    if (0..3).all(|axis| on_box_face(position[axis], at_corner[axis])) {
        600.0
    } else {
        330.0
    }
}

/// A trace corner the card has no answer for, in the state a step starts from, is refused by
/// both with one text, although every quadrature point is inside the card: the corner's mean
/// with slab B's 300 K is 450 K, above the card's 430 K, and a quadrature point reads at most
/// 2/3 of the corner, a mean of 405 K. Both ask the card at the start state first: the owner as
/// its Newton solve's first iterate, the kernel as its interface solve's starting value.
#[test]
fn a_trace_corner_outside_the_card_is_refused_by_the_kernel_and_by_the_owner_alike() {
    with_cx(|cx| {
        let owner = owner_run(cx, 2, hot_corner, 1).expect_err("the corner is outside the card");
        let owner = owner.to_string();
        let law = law()
            .factor(600.0, COLD_K)
            .expect_err("450 K is outside the card");
        eprintln!("the owner's refusal: {owner}");
        assert_eq!(law.to_string(), owner);
        assert!(owner.contains("fs-matdb refused"), "{owner}");
        // The kernel refuses with that text as the element's reason, and keeps both fields.
        let pair = kernel_run(cx, 2, Wall::Law, hot_corner, 1);
        let refused = pair.refused.expect("the guard row at the corner refuses");
        eprintln!("the kernel's refusal: {refused}");
        assert!(refused.contains(&format!("{owner:?}")), "{refused}");
        assert!(pair.fields.is_empty());
        let (_, mesh_a) = slab_body(cx, 2, RHO_C_A_J_M3_K, "right", LENGTH_M);
        let start: Vec<f64> = mesh_a.positions().iter().map(|&p| hot_corner(p)).collect();
        let bits = |field: &[f64]| field.iter().map(|t| t.to_bits()).collect::<Vec<_>>();
        assert_eq!(bits(&pair.kept[0]), bits(&start));
        assert!(pair.kept[1].iter().all(|t| t.to_bits() == COLD_K.to_bits()));
        // Without the guard rows the law is asked at the quadrature points alone, which are
        // inside the card, and the step is taken: the rows are what refuse.
        let unguarded = kernel_run(cx, 2, Wall::Unguarded, hot_corner, 1);
        assert_eq!(unguarded.refused, None);
        assert_eq!(unguarded.fields.len(), 1);
    });
}

/// The law's numbers against the card's knots by hand: the factor, its two slopes inside a
/// segment and at the knots, the refusals, and the other shapes of claim.
#[test]
fn the_law_is_the_card_over_its_reference_with_the_owner_newton_slope_at_a_knot() {
    let law = law();
    let (reference, _) = by_hand(REFERENCE_K, 0);
    close(
        law.reference_resistance().value_m2_k_per_w(),
        reference,
        1.0e-14,
        "the reference resistance",
    );
    // Inside a segment, at a knot (the right segment's slope), and at either end (the inward
    // segment's): `(T_a, T_b, the segment whose slope the law takes)`.
    for (kelvin_a, kelvin_b, segment) in [
        (320.0, 320.0, 0),
        (330.0, 290.0, 0),
        (400.0, 310.0, 1),
        (425.0, 415.0, 1),
        (390.0, 300.0, 1),
        (290.0, 270.0, 0),
        (460.0, 400.0, 1),
    ] {
        let at = 0.5 * kelvin_a + 0.5 * kelvin_b;
        let (resistance, slope) = by_hand(at, segment);
        let (factor, partials) = law.factor(kelvin_a, kelvin_b).unwrap();
        let what = format!("{kelvin_a} K and {kelvin_b} K");
        close(factor, reference / resistance, 1.0e-13, &what);
        let expected = -0.5 * slope / resistance * (reference / resistance);
        close(partials[0], expected, 1.0e-12, &what);
        close(partials[1], expected, 1.0e-12, &what);
    }
    // Away from the knots the slopes are those of central differences.
    let step_k = 1.0e-3;
    for (kelvin_a, kelvin_b) in [(330.0, 290.0), (400.0, 310.0)] {
        let (_, partials) = law.factor(kelvin_a, kelvin_b).unwrap();
        let by_a = (law.factor(kelvin_a + step_k, kelvin_b).unwrap().0
            - law.factor(kelvin_a - step_k, kelvin_b).unwrap().0)
            / (2.0 * step_k);
        let by_b = (law.factor(kelvin_a, kelvin_b + step_k).unwrap().0
            - law.factor(kelvin_a, kelvin_b - step_k).unwrap().0)
            / (2.0 * step_k);
        close(partials[0], by_a, 1.0e-7, "the slope by T_a");
        close(partials[1], by_b, 1.0e-7, "the slope by T_b");
    }
    // Outside the card there is no answer, and a temperature is positive.
    for (kelvin_a, kelvin_b, said) in [
        (470.0, 400.0, "fs-matdb refused"),
        (280.0, 270.0, "fs-matdb refused"),
        (
            -1.0,
            300.0,
            "absolute temperatures must be finite and positive",
        ),
    ] {
        let refusal = law.factor(kelvin_a, kelvin_b).unwrap_err().to_string();
        assert!(refusal.contains(said), "{refusal}");
    }
    // The identity is the card, the claim, the coordinates and the reference point.
    let words = |law: &TraceResistanceLaw| {
        let mut words = Vec::new();
        law.identity_words(&mut words);
        words
    };
    let elsewhere = TraceResistanceLaw::bind(
        NAME,
        &card(),
        &QueryPoint::new().with("T", 330.0).unwrap(),
        mean(),
        SelectionPolicy::SingleClaimOnly,
    )
    .unwrap();
    assert_eq!(words(&law).len(), 4);
    assert_eq!(words(&law), words(&self::law()));
    assert_ne!(words(&law), words(&elsewhere));
}

/// The other shapes of claim: a curve over side A's temperature alone, a curve over a
/// coordinate that is not a temperature (the shape of the committed pressure cards, which
/// state a temperature span only as validity), and a surface over both temperatures.
#[test]
fn the_law_carries_one_sided_and_level_claims_and_refuses_a_temperature_surface() {
    let sides = || TemperatureContactCoordinates::TwoSided(["TA", "TB"].map(String::from));
    let both = |a: f64, b: f64| {
        QueryPoint::new()
            .with_quantity("TA", kelvin(), a)
            .unwrap()
            .with_quantity("TB", kelvin(), b)
            .unwrap()
    };
    let span = || {
        ValidityDomain::unconstrained()
            .with_quantity("TA", kelvin(), 280.0, 430.0)
            .with_quantity("TB", kelvin(), 250.0, 450.0)
    };
    let bind = |card: &InterfaceSystemCard, point: &QueryPoint| {
        TraceResistanceLaw::bind(NAME, card, point, sides(), SelectionPolicy::SingleClaimOnly)
    };
    // Over side A's temperature alone: side B's has no slope.
    let one_sided = card_of(
        curve("TA", fs_qty::Temperature::DIMS, &KNOTS),
        InterpolationPolicy::LinearInside,
        span(),
    );
    let law = bind(&one_sided, &both(REFERENCE_K, REFERENCE_K)).unwrap();
    let (reference, _) = by_hand(REFERENCE_K, 0);
    let (resistance, slope) = by_hand(400.0, 1);
    let (factor, partials) = law.factor(400.0, 300.0).unwrap();
    close(factor, reference / resistance, 1.0e-13, "one-sided factor");
    close(
        partials[0],
        -slope / resistance * (reference / resistance),
        1.0e-12,
        "one-sided slope",
    );
    assert_eq!(partials[1].to_bits(), 0.0_f64.to_bits());
    // Over the pressure, with the temperatures only in the validity: the factor is one with no
    // slope inside the span, and outside it there is no answer.
    let pressure = card_of(
        curve(
            "normal_pressure",
            fs_qty::Pressure::DIMS,
            &[(1.0e5, 2.0e-3), (1.0e7, 5.0e-4)],
        ),
        InterpolationPolicy::LinearInside,
        span().with("normal_pressure", 1.0e5, 1.0e7),
    );
    let point = both(REFERENCE_K, REFERENCE_K)
        .with("normal_pressure", 2.0e6)
        .unwrap();
    let law = bind(&pressure, &point).unwrap();
    let (factor, partials) = law.factor(400.0, 300.0).unwrap();
    assert_eq!(
        [factor, partials[0], partials[1]].map(f64::to_bits),
        [1.0, 0.0, 0.0].map(f64::to_bits)
    );
    let refusal = law.factor(440.0, 300.0).unwrap_err().to_string();
    assert!(refusal.contains("fs-matdb refused"), "{refusal}");
    // A surface over both temperatures is not built.
    let surface = card_of(
        PropertyValue::Surface2D {
            axes: ["TA".to_string(), "TB".to_string()],
            axis_dims: [fs_qty::Temperature::DIMS; 2],
            samples: [
                [280.0, 250.0],
                [430.0, 250.0],
                [280.0, 450.0],
                [430.0, 450.0],
            ]
            .into_iter()
            .map(|at| (at, 1.0e-3 + 1.0e-6 * (at[0] - 280.0)))
            .collect(),
            triangles: vec![[0, 1, 2], [1, 3, 2]],
            dims: AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS,
        },
        InterpolationPolicy::LinearInside,
        span(),
    );
    let refusal = bind(&surface, &both(REFERENCE_K, REFERENCE_K)).unwrap_err();
    let refusal = refusal.to_string();
    assert!(refusal.contains("is not built"), "{refusal}");
}
