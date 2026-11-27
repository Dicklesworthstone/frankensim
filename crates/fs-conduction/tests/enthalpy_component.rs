//! The nonlinear enthalpy body of one material card on the composition kernel.
//!
//! `fs_conduction::component::EnthalpyBody` puts the owner
//! `EnthalpyBackwardEuler`, prepared from a card by
//! `MaterialEnthalpyTransport::prepare_uniform`, behind heat ports on named
//! boundary patches. The card is synthetic: density 8000 kg/m³, Cp rising
//! linearly from 450 J/kg/K at 250 K to 550 J/kg/K at 450 K, k from 12 to
//! 18 W/m/K over the same span, so the body's specific enthalpy above 250 K
//! is `h(T) = 450 u + u² / 4` with `u = T - 250 K`, and its inverse
//! `T(h) = 250 + 2 (sqrt(202500 + h) - 450)` is known in closed form. Here the
//! body passes the conformance kit, its ports read the closed-form `T(h)` of
//! its own state, a trace port steps it as a mean port does under a uniform
//! flux, a heater through the kernel steps it exactly as the owner steps
//! itself, and it refuses what it cannot carry. Its state codec reads what
//! it wrote, and its identity follows the chart, the ports and the step
//! policies. Its stage residual passes the kit, and a loop of the body, a
//! wall and a lump solved monolithically meets the same loop solved by
//! interface Newton. Run with
//! `cargo test -p fs-conduction --features component --test enthalpy_component`.
#![cfg(feature = "component")]
// Exact float comparisons here assert bitwise agreement with the owner.
#![allow(clippy::float_cmp)]

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_conduction::component::{EnthalpyBody, PatchPort};
use fs_conduction::fixtures::{box_grid, on_box_face};
use fs_conduction::material::LiquidMassFractionConductivity;
use fs_conduction::transient::enthalpy::material::{
    ExactMaterialEnthalpySource, MaterialEnthalpyTransport,
};
use fs_conduction::transient::enthalpy::{
    EnthalpyBackwardEuler, EnthalpyBudget, EnthalpyStepConfig,
};
use fs_conduction::{
    ConductionError, ConductionMesh, ConductionProblem, LinearConfig, ScalarField, ThermalBc,
    ThermalBoundary, ThermalBoundaryBuilder,
};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::phase::{HeatCapacitySourceQuery, SolidLiquidPhase};
use fs_ports::component::conformance::{
    Case, Observed, check_component, check_port_tangent, check_stage_residual,
};
use fs_ports::component::{Component, Proposal};
use fs_qty::{Density, Dims, QuantitySpec, Temperature};
use fs_solver::NewtonKrylovConfig;
use fs_system::elements::{Conductance, HeatCapacity, Heater, SourceAccount, Waveform};
use fs_system::{BlockSolveMode, BlockSolverPolicy, SystemBuilder};

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

const CP_DIMS: Dims = Dims([2, 0, -2, -1, 0, 0]);
const K_DIMS: Dims = Dims([1, 1, -3, -1, 0, 0]);
/// A 30 x 20 x 20 mm block of 72 tetrahedra on 36 vertices.
const LENGTH_M: f64 = 0.03;
/// The owner's joule gate on every step.
const ENERGY_TOLERANCE_J: f64 = 1.0e-8;

fn leak<T>(value: T) -> &'static T {
    Box::leak(Box::new(value))
}

fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "temperature",
            QuantitySpec::dimensional(Temperature::DIMS),
            t,
        )
        .unwrap()
}

/// The synthetic card: each property a straight line in temperature over
/// 250 K to 450 K, with no knot between.
fn card() -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, lo, hi) in [
        ("rho", Density::DIMS, 8000.0, 8000.0),
        ("cp", CP_DIMS, 450.0, 550.0),
        ("k", K_DIMS, 12.0, 18.0),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "temperature".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(250.0, lo), (450.0, hi)],
                    dims,
                },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "temperature",
                    QuantitySpec::dimensional(Temperature::DIMS),
                    250.0,
                    450.0,
                ),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: Vec::new(),
                provenance: Provenance {
                    source: "synthetic thermal fixture for the enthalpy body; not measured material data"
                        .into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic enthalpy body fixture".into(),
            phase: "solid".into(),
            process: "numerical reference".into(),
            revision: 0,
        },
        claims,
        Vec::new(),
    )
    .unwrap()
}

/// The card's storage and transport: Cp integrated from 0 J/kg at 250 K,
/// density frozen at 300 K, the scalar k(T) as isotropic conduction.
fn transport() -> &'static MaterialEnthalpyTransport {
    transport_from(0.0)
}

/// The same with the chart's specific enthalpy at 250 K given.
fn transport_from(lower_specific_enthalpy_j_kg: f64) -> &'static MaterialEnthalpyTransport {
    let card = card();
    leak(
        MaterialEnthalpyTransport::from_isotropic_card_exact(
            &card,
            ExactMaterialEnthalpySource {
                caloric: HeatCapacitySourceQuery {
                    heat_capacity: &PropertyKey::new("cp", CP_DIMS),
                    density: &PropertyKey::new("rho", Density::DIMS),
                    temperature_axis: "temperature",
                    lower: &point(250.0),
                    upper: &point(450.0),
                    selection: SelectionPolicy::SingleClaimOnly,
                },
                conductivity: &PropertyKey::new("k", K_DIMS),
                reference: &point(300.0),
                phase: SolidLiquidPhase::Solid,
                lower_specific_enthalpy_j_kg,
            },
        )
        .unwrap(),
    )
}

/// The card's specific enthalpy [J/kg] at `t` [K], in closed form.
fn h_of(t: f64) -> f64 {
    let u = t - 250.0;
    450.0 * u + 0.25 * u * u
}

/// The card's temperature [K] at specific enthalpy `h` [J/kg], in closed form.
fn t_of(h: f64) -> f64 {
    250.0 + 2.0 * ((202_500.0 + h).sqrt() - 450.0)
}

#[derive(Clone, Copy)]
struct Parts {
    transport: &'static MaterialEnthalpyTransport,
    mesh: &'static ConductionMesh,
    boundary: &'static ThermalBoundary,
    source: &'static ScalarField,
}

fn mesh() -> &'static ConductionMesh {
    let (complex, positions) = box_grid([3, 2, 2], [LENGTH_M, 0.02, 0.02]);
    leak(ConductionMesh::new(complex, positions).unwrap())
}

/// Two end faces, `left` (x = 0) and `right` (x = L), each with its own
/// law, and an insulated remainder.
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
        transport: transport(),
        mesh,
        boundary: boundary(mesh, ThermalBc::adiabatic(), ThermalBc::adiabatic()),
        source: leak(ScalarField::Uniform(0.0)),
    }
}

fn problem(parts: Parts) -> ConductionProblem<'static> {
    ConductionProblem {
        mesh: parts.mesh,
        boundary: parts.boundary,
        material: parts.transport.conductivity(),
        element_materials: None,
        source: parts.source,
    }
}

fn solver(cx: &Cx<'_>, parts: Parts) -> EnthalpyBackwardEuler<'static, 'static> {
    parts
        .transport
        .prepare_uniform(
            cx,
            parts.mesh,
            EnthalpyBudget {
                max_vertices: 36,
                max_elements: 72,
            },
        )
        .unwrap()
}

/// The owner's Newton policy: its solve stops at 1e-10 J, far below what
/// the kit's difference quotient can see, and every step closes its energy
/// to 1e-8 J.
fn config() -> EnthalpyStepConfig {
    EnthalpyStepConfig {
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-10,
            relative_tolerance: 1e-12,
            linear_restart: 20,
            max_linear_cycles: 8,
            forcing_minimum: 1e-12,
            forcing_maximum: 1e-4,
            ..NewtonKrylovConfig::default()
        },
        max_newton_iterations: 32,
        energy_tolerance_j: ENERGY_TOLERANCE_J,
    }
}

/// Each transpose solve of the tangent.
fn adjoint() -> LinearConfig {
    LinearConfig {
        tolerance: 1.0e-13,
        max_iterations: 2000,
        restart: 20,
    }
}

fn try_body(
    cx: &Cx<'_>,
    parts: Parts,
    ports: &[&str],
) -> Result<EnthalpyBody<'static, 'static>, ConductionError> {
    let ports: Vec<PatchPort<'_>> = ports.iter().map(|&name| PatchPort::Mean(name)).collect();
    EnthalpyBody::new(
        solver(cx, parts),
        problem(parts),
        &ports,
        config(),
        adjoint(),
    )
}

fn body(cx: &Cx<'_>, parts: Parts, ports: &[&str]) -> EnthalpyBody<'static, 'static> {
    try_body(cx, parts, ports).unwrap()
}

/// 310 K at the left end rising to 350 K at the right, as specific enthalpy.
fn initial(parts: Parts) -> Vec<f64> {
    parts
        .mesh
        .positions()
        .iter()
        .map(|p| h_of(310.0 + 40.0 * p[0] / LENGTH_M))
        .collect()
}

/// The area mean over the face `x = at` of the closed-form `T(h)`, by the
/// P1 rule on each boundary triangle.
fn patch_mean_k(parts: Parts, h: &[f64], at: f64) -> f64 {
    let (mut area, mut integral) = (0.0, 0.0);
    for face in parts.mesh.boundary() {
        if on_box_face(face.centroid[0], at) {
            let sum: f64 = face.vertices.iter().map(|&v| t_of(h[v as usize])).sum();
            area += face.area;
            integral += face.area * sum / 3.0;
        }
    }
    integral / area
}

/// Claims, declared before the first run: the body passes the conformance
/// kit on each case; the heat it reports in, `dt` times the port powers, is
/// stored to the owner's 1e-8 J gate plus rounding; its port tangent is the
/// difference quotient of its step within 1e-6 relative, above a floor of
/// the patch rounding plus 1e-9 K for where Newton stops; and its ports
/// read the closed-form `T(h)` of the state it returns, area-averaged over
/// each end, within 1e-9 K; and its stage residual passes the kit at that
/// state within 1e-8 J and 1e-6 relative.
#[test]
fn an_enthalpy_body_of_a_card_conforms_and_its_tangent_is_the_derivative_of_its_step() {
    with_cx(|cx| {
        let parts = parts();
        let body = body(cx, parts, &["left", "right"]);
        let old = initial(parts);
        let stale = vec![1.0e5; old.len()];
        // The stored energy is the inventory of the card's frozen density:
        // a field of 1 J/kg holds the block's mass in joules.
        let mass_kg = body.stored_energy_j(&vec![1.0; old.len()]).unwrap();
        let volume_m3 = LENGTH_M * 0.02 * 0.02;
        assert!((mass_kg - 8000.0 * volume_m3).abs() <= 1.0e-12 * mass_kg);
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
            let Observed::Stepped {
                receipt,
                audit,
                outputs,
            } = observed
            else {
                panic!("{dt_s} s, {heat_w:?} W: the body asked for a shorter step");
            };
            let heat_j = dt_s * (heat_w[0] + heat_w[1]);
            assert!(
                receipt.numerical_loss_j.abs()
                    <= ENERGY_TOLERANCE_J + 1.0e-12 * (1.0 + heat_j.abs()),
                "{dt_s} s, {heat_w:?} W: {receipt:?}"
            );
            assert!((audit.port_energy_in_j - heat_j).abs() <= 1.0e-12 * (1.0 + heat_j.abs()));
            // The ports read the card's T(h) of the state the body returns.
            let mut new = Vec::new();
            let mut means = [0.0; 2];
            let proposal = body
                .propose(&old, dt_s, &heat_w, &mut new, &mut means, cx)
                .unwrap();
            assert!(matches!(proposal, Proposal::Step(_)));
            assert_eq!(means.as_slice(), outputs.as_slice());
            check_stage_residual(&body, &case, &new, 1.0e-6, 1.0e-8, cx)
                .unwrap_or_else(|failure| panic!("stage at {dt_s} s, {heat_w:?} W: {failure}"));
            for (port, at) in [(0, 0.0), (1, LENGTH_M)] {
                let expected = patch_mean_k(parts, &new, at);
                eprintln!(
                    "{dt_s} s, {heat_w:?} W: port {port} reads {} K, closed form {expected} K",
                    means[port]
                );
                assert!(
                    (means[port] - expected).abs() <= 1.0e-9,
                    "port {port}: {} K against {expected} K",
                    means[port]
                );
            }
            // Insulated, the field relaxes: the cold end warms, the hot end cools.
            if heat_w == [0.0, 0.0] {
                assert!(means[0] > 310.0 && means[1] < 350.0, "{means:?}");
            }
            let difference_step_w =
                heat_w.iter().fold(1.0_f64, |m, q| m.max(q.abs())) * f64::EPSILON.cbrt();
            let floor_k_w = (8.0 * f64::EPSILON * 400.0 + 1.0e-9) / difference_step_w;
            check_port_tangent(&body, &case, 1.0e-6, floor_k_w, cx)
                .unwrap_or_else(|failure| panic!("tangent at {dt_s} s, {heat_w:?} W: {failure}"));
        }
    });
}

/// Claims, declared before the first run: given the consistent loads of a
/// uniform flux, `A_f / 3` of each face's share on each of its vertices, a
/// trace port steps the body as the mean port does, the two end fields'
/// closed-form temperatures within 1e-9 K; a trace port returns the end
/// `T(h)` at each of its vertices, and its stage residual passes the kit at
/// that end as above; under loads that are not a uniform flux the traced
/// body passes the kit and the tangent check as above.
#[test]
fn a_trace_port_of_an_enthalpy_body_conforms_and_steps_it_as_the_mean_port_under_a_uniform_flux() {
    with_cx(|cx| {
        let parts = parts();
        let ports = [PatchPort::Trace("left"), PatchPort::Mean("right")];
        let traced = EnthalpyBody::new(
            solver(cx, parts),
            problem(parts),
            &ports,
            config(),
            adjoint(),
        )
        .unwrap();
        let vertices = traced.trace_vertices(0).unwrap().to_vec();
        assert_eq!(vertices.len(), 9);
        assert!(traced.trace_vertices(1).is_none());
        let old = initial(parts);
        let stale = vec![1.0e5; old.len()];
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
            let (mut end, mut readings) = (Vec::new(), vec![0.0; 10]);
            traced
                .propose(&old, dt_s, &inputs, &mut end, &mut readings, cx)
                .unwrap();
            let (mut expected, mut means) = (Vec::new(), [0.0; 2]);
            meaned
                .propose(
                    &old,
                    dt_s,
                    &[left_w, right_w],
                    &mut expected,
                    &mut means,
                    cx,
                )
                .unwrap();
            let gap_k = end
                .iter()
                .zip(&expected)
                .map(|(a, b)| (t_of(*a) - t_of(*b)).abs())
                .fold(0.0, f64::max);
            eprintln!("{dt_s} s, {left_w} W: trace against mean port, largest gap {gap_k:e} K");
            assert!(gap_k <= 1.0e-9, "{dt_s} s: {gap_k} K");
            for (reading, &vertex) in readings.iter().zip(&vertices) {
                assert!((reading - t_of(end[vertex])).abs() <= 1.0e-9);
            }
            assert!((readings[9] - means[1]).abs() <= 1.0e-9);
            let uniform_case = Case {
                state: &old,
                stale: &stale,
                dt_s,
                inputs: &inputs,
            };
            check_stage_residual(&traced, &uniform_case, &end, 1.0e-6, 1.0e-8, cx)
                .unwrap_or_else(|failure| panic!("stage at {dt_s} s: {failure}"));
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
            assert!(
                receipt.numerical_loss_j.abs()
                    <= ENERGY_TOLERANCE_J + 1.0e-12 * (1.0 + heat_j.abs())
            );
            assert!((audit.port_energy_in_j - heat_j).abs() <= 1.0e-12 * (1.0 + heat_j.abs()));
            let difference_step_w =
                inputs.iter().fold(1.0_f64, |m, q| m.max(q.abs())) * f64::EPSILON.cbrt();
            let floor_k_w = (8.0 * f64::EPSILON * 400.0 + 1.0e-9) / difference_step_w;
            check_port_tangent(&traced, &case, 1.0e-6, floor_k_w, cx)
                .unwrap_or_else(|failure| panic!("tangent at {dt_s} s: {failure}"));
        }
    });
}

/// Claims, declared before the first run: a heater through the kernel
/// steps the body bit for bit as the owner's `advance` steps itself under
/// the heater's flux; the 80 J it puts in are stored to eight steps of the
/// owner's gate; the ledger closes; the heater generates entropy.
#[test]
fn a_heater_through_the_kernel_steps_an_enthalpy_body_as_the_owner_does() {
    const POWER_W: f64 = 20.0;
    const DT_S: f64 = 0.5;
    with_cx(|cx| {
        let parts = parts();
        let reference = body(cx, parts, &["left"]);
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
        assert!(system.guess(heat, &[310.0]));

        // The owner driven directly: its left end prescribes the flux -P / A.
        let area_m2 = reference.patch_area_m2(0);
        let heated = Parts {
            boundary: boundary(
                parts.mesh,
                ThermalBc::neumann(-POWER_W / area_m2).unwrap(),
                ThermalBc::adiabatic(),
            ),
            ..parts
        };
        let owner = solver(cx, parts);
        let mut field = initial(parts);
        let start_j = reference.stored_energy_j(&field).unwrap();
        for index in 0..8 {
            system.step(DT_S, cx).unwrap();
            field = owner
                .advance(cx, problem(heated), None, &field, DT_S, config())
                .unwrap()
                .specific_enthalpy_j_kg;
            let kernel = system.state::<Vec<f64>>(node).unwrap();
            assert!(
                kernel
                    .iter()
                    .zip(&field)
                    .all(|(a, b)| a.to_bits() == b.to_bits()),
                "step {index}: the kernel's body differs from the owner's"
            );
        }
        let stored_j = reference.stored_energy_j(&field).unwrap() - start_j;
        eprintln!(
            "enthalpy body, 8 steps of {DT_S} s at {POWER_W} W: stored {stored_j} J, left end {} K",
            patch_mean_k(parts, &field, 0.0)
        );
        assert!(
            (stored_j - 8.0 * DT_S * POWER_W).abs() <= 8.0 * ENERGY_TOLERANCE_J + 1.0e-12 * 80.0,
            "{stored_j} J"
        );
        let ledger = system.ledger();
        assert!(
            ledger.energy_residual_j().abs() <= 1.0e-9 * start_j,
            "{}",
            ledger.energy_residual_j()
        );
        assert!(ledger.entropy_generated_j_per_k() > 0.0);
    });
}

/// Claims, declared before the first run: with every node offering its
/// stage residual, the kernel's default policy solves the loop of a 20 J/K
/// lump at 420 K, a 0.5 W/K wall and the body's left end monolithically,
/// and the body's receipt there has a zero residual. Over 50 steps of 2 s
/// the body's closed-form temperatures and the lump stay within 1e-8 K of
/// the same loop solved by interface Newton: the owner's Newton stops at
/// 1e-10 J, at most 7e-10 K a step at Cp >= 450 J/kg/K on a vertex of at
/// least 3.3e-4 kg (asserted), and the body's slowest decay keeps at most 14
/// such steps. On every step the
/// body's monolithic receipt books the entropy its own `propose` books at
/// the converged port heat within 1e-8 relative.
#[test]
fn a_monolithic_solve_of_the_lump_loop_meets_interface_newton_and_the_owner_receipt() {
    const DT_S: f64 = 2.0;
    with_cx(|cx| {
        let parts = parts();
        let probe = body(cx, parts, &["left"]);
        let lightest_kg = solver(cx, parts)
            .reference_nodal_masses_kg()
            .iter()
            .fold(f64::INFINITY, |m, &v| m.min(v));
        assert!(lightest_kg >= 3.3e-4, "{lightest_kg} kg");
        let build = || {
            let mut builder = SystemBuilder::new();
            let body = body(cx, parts, &["left"]);
            let node = builder
                .add_with_tangent_and_stage_residual("body", body, initial(parts))
                .unwrap();
            let capacity = HeatCapacity::new(20.0).unwrap();
            let lump = builder
                .add_with_tangent_and_stage_residual("lump", capacity, 420.0)
                .unwrap();
            let conductance = Conductance::new(0.5).unwrap();
            let wall = builder
                .add_with_tangent_and_stage_residual("wall", conductance, ())
                .unwrap();
            let (a, b) = (
                builder.port(wall, "a").unwrap(),
                builder.port(wall, "b").unwrap(),
            );
            builder
                .bond(builder.port(lump, "heat").unwrap(), a)
                .unwrap();
            let left = builder.port(node, "left").unwrap();
            builder.bond(b, left).unwrap();
            let mut system = builder.seal().unwrap();
            assert!(system.guess(a, &[420.0]) && system.guess(b, &[310.0]));
            (system, node, lump, left)
        };
        let (mut monolithic, node, lump, left) = build();
        let (mut interface, ..) = build();
        interface.set_solver_policy(BlockSolverPolicy::InterfaceNewton);
        assert_eq!(
            monolithic.block_modes(),
            [BlockSolveMode::MonolithicNewtonKrylov]
        );
        assert_eq!(interface.block_modes(), [BlockSolveMode::InterfaceNewton]);
        let (mut worst_k, mut worst_entropy) = (0.0_f64, 0.0_f64);
        for index in 0..50 {
            let old = monolithic.state::<Vec<f64>>(node).unwrap().clone();
            monolithic.step(DT_S, cx).unwrap();
            interface.step(DT_S, cx).unwrap();
            let lumps = [&monolithic, &interface].map(|system| *system.state::<f64>(lump).unwrap());
            let fields =
                [&monolithic, &interface].map(|system| system.state::<Vec<f64>>(node).unwrap());
            let gap_k = fields[0]
                .iter()
                .zip(fields[1])
                .map(|(a, b)| (t_of(*a) - t_of(*b)).abs())
                .fold((lumps[0] - lumps[1]).abs(), f64::max);
            let heat_w = monolithic.last_input(left).to_vec();
            let (mut new, mut reading) = (Vec::new(), [0.0]);
            let proposal = probe.propose(&old, DT_S, &heat_w, &mut new, &mut reading, cx);
            let Ok(Proposal::Step(owner)) = proposal else {
                panic!("step {index}: {proposal:?}");
            };
            let receipt = &monolithic.last_receipts()[node.index()];
            assert_eq!(receipt.residual, 0.0, "step {index}: {receipt:?}");
            let booked = receipt.entropy_generated_j_per_k;
            let owner_j_k = owner.entropy_generated_j_per_k;
            let entropy = (booked - owner_j_k).abs() / owner_j_k;
            worst_k = worst_k.max(gap_k);
            worst_entropy = worst_entropy.max(entropy);
            assert!(
                gap_k <= 1.0e-8 && entropy <= 1.0e-8,
                "step {index}: {gap_k} K; {booked} against {owner_j_k} J/K"
            );
        }
        eprintln!(
            "monolithic against interface Newton over 50 steps: largest gap {worst_k:e} K, entropy {worst_entropy:e} relative"
        );
    });
}

/// Claims, declared before the first run: the conduction body's steady slab
/// (10 W in at the left, out through a 0.5 W/K wall into a 1e11 J/K lump,
/// 500 s steps, steady at a step storing at most 1e-9 of `Q dt`, within 100
/// steps) as an enthalpy body. Its own step from there generates
/// `dt Q (1/T_right - 1/T_left)` of its port temperatures within 1e-6: the
/// card's k(T) makes the discrete field not exactly planar, and a face's
/// spread enters only at second order. The ledger books it within 1e-8.
#[test]
fn a_steady_enthalpy_slab_generates_its_heat_times_the_difference_of_its_end_reciprocals() {
    const POWER_W: f64 = 10.0;
    const DT_S: f64 = 500.0;
    with_cx(|cx| {
        let parts = parts();
        let probe = body(cx, parts, &["left", "right"]);
        let mut builder = SystemBuilder::new();
        let ends = body(cx, parts, &["left", "right"]);
        let node = builder
            .add_with_tangent("body", ends, initial(parts))
            .unwrap();
        let heater = Heater::new(Waveform::constant(POWER_W).unwrap()).unwrap();
        let heater = builder
            .add_with_tangent("heater", heater, SourceAccount::default())
            .unwrap();
        let wall = builder
            .add_with_tangent("wall", Conductance::new(0.5).unwrap(), ())
            .unwrap();
        let lump = builder
            .add_with_tangent("lump", HeatCapacity::new(1.0e11).unwrap(), 300.0)
            .unwrap();
        let heat = builder.port(heater, "heat").unwrap();
        let a = builder.port(wall, "a").unwrap();
        let b = builder.port(wall, "b").unwrap();
        for (one, other) in [
            (builder.port(node, "left").unwrap(), heat),
            (builder.port(node, "right").unwrap(), a),
            (b, builder.port(lump, "heat").unwrap()),
        ] {
            builder.bond(one, other).unwrap();
        }
        let mut system = builder.seal().unwrap();
        assert!(system.guess(heat, &[370.0]) && system.guess(a, &[320.0]));
        assert!(system.guess(b, &[300.0]));
        let booked =
            |system: &fs_system::System| system.ledger().nodes[0].entropy_generated_j_per_k.value();
        let (mut steps, mut stored_j, mut last_j_per_k) = (0, f64::INFINITY, 0.0);
        while stored_j.abs() > 1.0e-9 * POWER_W * DT_S {
            assert!(steps < 100, "not steady after {steps} steps: {stored_j} J");
            let old = system.state::<Vec<f64>>(node).unwrap().clone();
            let before = booked(&system);
            system.step(DT_S, cx).unwrap();
            let field = system.state::<Vec<f64>>(node).unwrap();
            stored_j = probe.stored_energy_j(field).unwrap() - probe.stored_energy_j(&old).unwrap();
            (steps, last_j_per_k) = (steps + 1, booked(&system) - before);
        }
        assert_eq!(system.ledger().nodes[0].name, "body");
        let field = system.state::<Vec<f64>>(node).unwrap();
        let (mut end, mut ends_k) = (Vec::new(), [0.0; 2]);
        let heat_w = [POWER_W, -POWER_W];
        let proposal = probe.propose(field, DT_S, &heat_w, &mut end, &mut ends_k, cx);
        let Ok(Proposal::Step(receipt)) = proposal else {
            panic!("the body did not step");
        };
        let generated = receipt.entropy_generated_j_per_k;
        let expected = DT_S * POWER_W * (1.0 / ends_k[1] - 1.0 / ends_k[0]);
        let gap = (generated - expected).abs() / expected;
        let booked_gap = (last_j_per_k - generated).abs() / expected;
        eprintln!(
            "steady after {steps} steps: {generated} J/K a step, dt Q (1/Tc - 1/Th) {expected} J/K, {gap:e} relative; the ledger's last step {booked_gap:e} from it"
        );
        assert!(gap <= 1.0e-6 && booked_gap <= 1.0e-8);
    });
}

/// A Newton policy whose absolute floor, 1e-30 J, is far below its relative
/// target, 1e-10 of the step's first residual: the owner stops on the
/// relative one.
fn relative_config() -> EnthalpyStepConfig {
    let mut config = config();
    config.newton.absolute_tolerance = 1.0e-30;
    config.newton.relative_tolerance = 1.0e-10;
    config
}

/// The receipt's tolerance is the target the owner stopped at. The first
/// scene of a card's heat capacity (fs-cli, 1e-10 relative and 1e-10 J)
/// found the receipt carrying the absolute tolerance alone: the kernel
/// refused a step the owner had accepted at 3.1e-10 J.
#[test]
fn an_enthalpy_body_whose_newton_stops_on_its_relative_tolerance_steps_through_the_kernel() {
    const POWER_W: f64 = 20.0;
    const DT_S: f64 = 0.5;
    with_cx(|cx| {
        let parts = parts();
        let make = || {
            EnthalpyBody::new(
                solver(cx, parts),
                problem(parts),
                &[PatchPort::Mean("left")],
                relative_config(),
                adjoint(),
            )
            .unwrap()
        };
        let mut new = Vec::new();
        let mut out = [0.0];
        let Proposal::Step(receipt) = make()
            .propose(&initial(parts), DT_S, &[POWER_W], &mut new, &mut out, cx)
            .unwrap()
        else {
            panic!("the body did not step");
        };
        eprintln!(
            "relative Newton target: residual {:e} J at tolerance {:e} J",
            receipt.residual, receipt.tolerance
        );
        assert!(receipt.tolerance > 1.0e-30, "{receipt:?}");
        assert!(receipt.residual <= receipt.tolerance, "{receipt:?}");
        let mut builder = SystemBuilder::new();
        let node = builder
            .add_with_tangent("body", make(), initial(parts))
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
        assert!(system.guess(heat, &[310.0]));
        for _ in 0..8 {
            system.step(DT_S, cx).unwrap();
        }
        assert_eq!(system.steps(), 8);
    });
}

#[test]
fn an_enthalpy_body_refuses_what_it_cannot_carry_and_a_step_off_its_chart() {
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
            Parts {
                source: leak(ScalarField::Uniform(1.0e4)),
                ..parts
            },
            &["left"],
            "volumetric source",
        );
        refused(parts, &["top"], "port top names no boundary region");
        refused(parts, &["left", "left"], "port left is named twice");
        refused(parts, &[], "at least one port");
        // A solver prepared on another mesh of the same shape.
        let other = Parts {
            mesh: mesh(),
            ..parts
        };
        match EnthalpyBody::new(
            solver(cx, other),
            problem(parts),
            &[PatchPort::Mean("left")],
            config(),
            adjoint(),
        ) {
            Err(ConductionError::Config { what, .. }) => {
                assert!(what.contains("same mesh"), "{what}");
            }
            other => panic!("another mesh: {other:?}"),
        }
        // A step the owner refuses comes back as the owner's error and
        // publishes nothing: a field of the wrong length, a heat that is
        // not a number, and 1e5 W for one second, which would carry the
        // block far past the chart's 450 K end.
        let body = body(cx, parts, &["left"]);
        let old = initial(parts);
        let mut new = Vec::new();
        let mut out = [0.0];
        assert!(
            body.propose(&vec![h_of(300.0); 3], 1.0, &[1.0], &mut new, &mut out, cx)
                .is_err()
        );
        assert!(
            body.propose(&old, 1.0, &[f64::NAN], &mut new, &mut out, cx)
                .is_err()
        );
        let off_chart = body.propose(&old, 1.0, &[1.0e5], &mut new, &mut out, cx);
        eprintln!("1e5 W for 1 s: {off_chart:?}");
        assert!(off_chart.is_err());
        assert!(new.is_empty());
    });
}

/// Claims, declared before the first run. An enthalpy body writes its state
/// as the vertex count and every specific enthalpy's bits and reads the
/// same bits back, zero included. It refuses its words with another count,
/// a word short, and a value that is NaN or infinite. Its identity is four
/// words, equal for two bodies built the same way and different, from that
/// and from each other, for a chart that starts at another specific
/// enthalpy, another port set, another step policy and another policy of
/// the tangent's transpose solves. A body on a solver with more than its
/// mesh, masses and one chart, here a phase multiplier on the conductivity,
/// writes no words and reads none.
#[test]
fn an_enthalpy_body_reads_the_state_it_wrote_and_its_identity_follows_its_chart_and_policies() {
    with_cx(|cx| {
        let parts = parts();
        let ends = [PatchPort::Mean("left"), PatchPort::Mean("right")];
        let build = |parts: Parts, ports: &[PatchPort<'_>], config, adjoint| {
            EnthalpyBody::new(solver(cx, parts), problem(parts), ports, config, adjoint).unwrap()
        };
        let identity = |body: &EnthalpyBody<'_, '_>| {
            let mut words = Vec::new();
            body.identity_words(&mut words);
            words
        };
        let first = build(parts, &ends, config(), adjoint());
        assert!(<EnthalpyBody<'static, 'static> as Component>::STATE_CHARTER.is_some());

        let mut field = initial(parts);
        field[0] = 0.0;
        let mut words = Vec::new();
        assert!(first.encode_state(&field, &mut words));
        assert_eq!(
            (words.len(), words[0]),
            (field.len() + 1, field.len() as u64)
        );
        let read = first.decode_state(&words).expect("the body's own words");
        assert_eq!(read.len(), field.len());
        assert!(
            read.iter()
                .zip(&field)
                .all(|(a, b)| a.to_bits() == b.to_bits())
        );
        let refused = |change: &dyn Fn(&mut Vec<u64>)| {
            let mut changed = words.clone();
            change(&mut changed);
            first.decode_state(&changed).is_none()
        };
        assert!(refused(&|words| words[0] -= 1), "another count");
        assert!(
            refused(&|words| words.truncate(words.len() - 1)),
            "a word short"
        );
        for value in [f64::NAN, f64::NEG_INFINITY] {
            assert!(refused(&|words| words[2] = value.to_bits()), "{value} J/kg");
        }

        let shifted = Parts {
            transport: transport_from(1.0),
            ..parts
        };
        let looser = LinearConfig {
            tolerance: 1.0e-12,
            ..adjoint()
        };
        let identities = [
            ("the first body", identity(&first)),
            (
                "another chart",
                identity(&build(shifted, &ends, config(), adjoint())),
            ),
            (
                "another port set",
                identity(&build(parts, &ends[..1], config(), adjoint())),
            ),
            (
                "another step policy",
                identity(&build(parts, &ends, relative_config(), adjoint())),
            ),
            (
                "another transpose policy",
                identity(&build(parts, &ends, config(), looser)),
            ),
        ];
        assert_eq!(
            identity(&build(parts, &ends, config(), adjoint())),
            identities[0].1
        );
        for (index, (name, words)) in identities.iter().enumerate() {
            assert_eq!(words.len(), 4, "{name}");
            for (other, earlier) in &identities[..index] {
                assert_ne!(words, earlier, "{name} has the identity of {other}");
            }
        }

        // A phase multiplier on the conductivity is in no identity, so a
        // body on a solver with one writes nothing and reads nothing.
        let law = LiquidMassFractionConductivity::declared(1.0, 0.5).unwrap();
        let phased = solver(cx, parts)
            .with_phase_conductivity(cx, vec![law; parts.mesh.element_count()])
            .unwrap();
        let phased = EnthalpyBody::new(phased, problem(parts), &ends, config(), adjoint()).unwrap();
        let mut none = Vec::new();
        assert!(!phased.encode_state(&field, &mut none) && none.is_empty());
        assert!(phased.decode_state(&words).is_none());
    });
}
