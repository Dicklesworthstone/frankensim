//! The transient finite-element body on the composition kernel.
//!
//! `fs_conduction::component::ConductionBody` puts the backward-Euler owner
//! behind heat ports on named boundary patches. Here it passes the
//! conformance kit, a heater through the kernel steps it exactly as the
//! owner steps itself, and a lump joined to it through a conductance meets
//! the step condensed from the owner's own responses and reaches the two
//! bodies' shared equilibrium. Run with
//! `cargo test -p fs-conduction --features component --test conduction_component`.
#![cfg(feature = "component")]
// Exact float comparisons here assert bitwise agreement with the owner.
#![allow(clippy::float_cmp)]

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_conduction::component::ConductionBody;
use fs_conduction::fixtures::{box_grid, on_box_face};
use fs_conduction::transient::VolumetricHeatCapacity;
use fs_conduction::transient::backward_euler::{BackwardEuler, StepConfig};
use fs_conduction::{
    ConductionError, ConductionMesh, ConductionProblem, ConductivityModel, ConductivityTable,
    LinearConfig, ScalarField, ThermalBc, ThermalBoundary, ThermalBoundaryBuilder,
};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
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
    ConductionBody::new(solver(cx, parts.mesh), problem(parts), ports, config())
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

#[test]
fn a_lump_through_a_conductance_meets_the_condensed_step_and_the_shared_equilibrium() {
    const LUMP_J_K: f64 = 20.0;
    const LUMP_START_K: f64 = 420.0;
    const WALL_W_K: f64 = 0.5;
    const DT_S: f64 = 2.0;
    const REFERENCE_W: f64 = 10.0;
    // 3000 s: some fifty times the slowest joint mode (about 60 s).
    const STEPS: usize = 1500;
    with_cx(|cx| {
        let parts = parts();
        // The owner's own responses, read through an unconnected copy.
        let probe = body(cx, parts, &["left"]);
        let mut builder = SystemBuilder::new();
        let node = builder
            .add_with_tangent("body", body(cx, parts, &["left"]), initial(parts))
            .unwrap();
        let lump = builder
            .add_with_tangent("lump", HeatCapacity::new(LUMP_J_K).unwrap(), LUMP_START_K)
            .unwrap();
        let wall = builder
            .add_with_tangent("wall", Conductance::new(WALL_W_K).unwrap(), ())
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
            let heat_w = WALL_W_K * (lump_k - at_zero[0])
                / (1.0 + WALL_W_K * (0.5 * DT_S / LUMP_J_K + slope_k_w));
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
        match ConductionBody::new(other, problem(parts), &["left"], config()) {
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
