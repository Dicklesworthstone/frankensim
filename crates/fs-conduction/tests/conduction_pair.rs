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
//! asserted. Run with
//! `cargo test -p fs-conduction --features component --test conduction_pair -- --nocapture`.
#![cfg(feature = "component")]

use std::f64::consts::{FRAC_PI_2, PI};

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_conduction::component::ConductionBody;
use fs_conduction::fixtures::{box_grid, on_box_face};
use fs_conduction::transient::VolumetricHeatCapacity;
use fs_conduction::transient::backward_euler::{BackwardEuler, StepConfig};
use fs_conduction::{
    AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS, AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY,
    ConductionMesh, ConductionProblem, ConductivityModel, InterfaceFacePair, InterfaceResistance,
    InterfaceSurface, LinearConfig, ScalarField, ThermalBc, ThermalBoundary,
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
use fs_rep_mesh::TetComplex;
use fs_system::SystemBuilder;
use fs_system::elements::Conductance;

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
    port: &str,
    port_x: f64,
) -> (ConductionBody<'static>, &'static ConductionMesh) {
    let (complex, positions) = box_grid([cells, cross, cross], [LENGTH_M, CROSS_M, CROSS_M]);
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

/// What one path reports: the energy that crossed by `t*`, and the end
/// field's largest spread across a plane x = const.
struct Crossed {
    energy_j: f64,
    spread_k: f64,
}

/// The pair on the kernel: A's `right` end bonded to B's `left` end through
/// the card's conductance. Checks the pair's energy and the ledger close.
fn kernel_pair_crossed(cx: &Cx<'_>, cells: usize, cross_a: usize, cross_b: usize) -> Crossed {
    let (a, mesh_a) = slab_body(cx, cells, cross_a, "right", LENGTH_M);
    let (b, mesh_b) = slab_body(cx, cells, cross_b, "left", 0.0);
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
    let (probe_a, _) = slab_body(cx, cells, cross_a, "right", LENGTH_M);
    let (probe_b, _) = slab_body(cx, cells, cross_b, "left", 0.0);
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
    }
}

/// The owner's own matching-P1 contact: both slabs in one mesh with
/// coincident, unshared traces at x = L, bound by the same card.
fn owner_contact_crossed(cx: &Cx<'_>, cells: usize, cross: usize) -> Crossed {
    let (slab, slab_positions) = box_grid([cells, cross, cross], [LENGTH_M, CROSS_M, CROSS_M]);
    let per_slab = slab_positions.len();
    let offset = u32::try_from(per_slab).unwrap();
    let mut tets = slab.tets.clone();
    tets.extend(slab.tets.iter().map(|tet| tet.map(|v| v + offset)));
    let mut positions = slab_positions.clone();
    positions.extend(slab_positions.iter().map(|&[x, y, z]| [x + LENGTH_M, y, z]));
    let mesh =
        ConductionMesh::new(TetComplex::from_tets(positions.len(), tets), positions).unwrap();
    let boundary: ThermalBoundary = ThermalBoundaryBuilder::new(&mesh)
        .adiabatic_remainder()
        .finish()
        .unwrap();
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
    let interfaces = ThermalInterfaces::new(
        &mesh,
        &boundary,
        vec![InterfaceSurface::new("slab-contact", pairs, contact()).unwrap()],
    )
    .unwrap();
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
            let owner = owner_contact_crossed(cx, cells, 2);
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
