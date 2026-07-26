//! End-to-end conjugate exchange: a solved fan operating point drives a
//! validity-gated convection card, whose coefficient couples a real
//! `fs-conduction` FEM solve to a stream-wise air path until the interface is
//! consistent.
//!
//! This is the fixture the E05 vertical did not have. `plate_fin_e2e.rs`
//! builds the same chain ONE WAY: it computes a bulk air temperature from the
//! *declared* base power, hands that one number to every Robin row, and never
//! feeds the solved surface back. That is why its own comment says "This is
//! not a conjugate fluid solve." Here the solid's surface temperature sets how
//! much heat the air takes, the air's temperature sets how much heat the solid
//! loses, and the loop closes on a declared criterion.
//!
//! What this fixture is NOT: a validated electronics-cooling model, a
//! commercial fan, a real duct, or any claim about a physical enclosure. The
//! fan curve and loss network are synthetic. Its value is that every link is
//! the production surface — a certified operating-point root, a card that
//! refuses outside its domain, a real assembled FEM solve — rather than a
//! stub, so the coupling is exercised against code that can disagree with it.

use fs_airflow::conjugate::{
    AirPath, AirSegment, ConjugateConfig, Relaxation, SolidRegionState, solve_conjugate,
};
use fs_airflow::{
    EnclosureNetwork, FanArrangement, FanBank, FanCurve, FanPoint, LeakageElement,
    LossElement, LossNetwork, LossResistance, OperatingPoint, SourceProvenance, ToleranceBasis,
    solve_operating_point,
};
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_conduction::fixtures::{box_grid, on_box_face};
use fs_conduction::{
    BoundaryFace, ConductionMesh, ConductionProblem, ConductivityModel, InitialGuess, LinearConfig,
    Nonlinearity, ScalarField, SolveConfig, StopRule, ThermalBc, ThermalBoundary,
    ThermalBoundaryBuilder, solve,
};
use fs_convection::{CorrelationId, ThermalConductivity, evaluate};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_qty::{Area, Density, DynViscosity, Length, Pressure, VolumetricFlowRate};
use fs_rep_mesh::TetComplex;

// --- Solid ------------------------------------------------------------------

/// Stream-wise length, m. `L/D_h = 60` clears the cards' inclusive
/// fully-developed floor of 50 with 20% margin.
const BAR_LENGTH_M: f64 = 0.24;
const BAR_WIDTH_M: f64 = 0.02;
const BAR_HEIGHT_M: f64 = 0.004;
const CELLS_X: usize = 8;
const CELLS_Y: usize = 2;
const CELLS_Z: usize = 2;
/// Stream-wise Robin regions. `CELLS_X` is a multiple of this so no cell
/// straddles a region boundary.
const SEGMENTS: usize = 4;

/// A deliberately POOR conductor (alumina/FR4 class, not aluminium). An
/// aluminium bar at this power is isothermal to within millikelvin, which
/// would hide the stream-wise wall gradient the coupling produces — the
/// fixture would still converge and still balance, but it would no longer
/// demonstrate anything the old one-way model could not.
const SOLID_CONDUCTIVITY_W_M_K: f64 = 5.0;
/// Total dissipated power, W. Sized for a ~25 K air rise at the solved flow.
const TOTAL_POWER_W: f64 = 1.2;

// --- Air --------------------------------------------------------------------

const AIR_DENSITY_KG_M3: f64 = 1.2;
const AIR_DYNAMIC_VISCOSITY_PA_S: f64 = 1.8e-5;
const AIR_THERMAL_CONDUCTIVITY_W_M_K: f64 = 0.026;
const AIR_SPECIFIC_HEAT_J_KG_K: f64 = 1_005.0;
const AIR_PRANDTL: f64 = 0.72;
const INLET_TEMPERATURE_K: f64 = 300.0;
const HYDRAULIC_DIAMETER_M: f64 = 0.004;
const CHANNEL_BRANCH: &str = "duct";

/// The cards' inclusive fully-developed floor. Checked at compile time so a
/// later geometry edit cannot quietly walk the fixture out of the domain and
/// leave the runtime admission as the only guard.
const _: () = assert!(BAR_LENGTH_M / HYDRAULIC_DIAMETER_M >= 50.0);

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0x0000_C0FF_EE00_0002,
                kernel_id: 72,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn synthetic_source(identifier: &str) -> SourceProvenance {
    SourceProvenance::new(
        "synthetic conjugate-coupling fixture; not a manufacturer curve",
        identifier,
    )
}

fn fan_curve() -> FanCurve {
    FanCurve::new(
        "synthetic-conjugate-fan",
        vec![
            FanPoint::new(VolumetricFlowRate::new(0.0), Pressure::new(30.0)),
            FanPoint::new(VolumetricFlowRate::new(30.0e-6), Pressure::new(25.0)),
            FanPoint::new(VolumetricFlowRate::new(60.0e-6), Pressure::new(12.0)),
            FanPoint::new(VolumetricFlowRate::new(90.0e-6), Pressure::new(0.0)),
        ],
        synthetic_source("conjugate-fan-v1"),
        0.0,
        ToleranceBasis::EngineeringAllowance,
        VolumetricFlowRate::new(5.0e-6),
        (0.7, 3.0),
    )
    .expect("synthetic fan curve")
}

fn loss(name: &str, resistance_pa_s2_m6: f64) -> LossElement {
    LossElement::new(
        name,
        LossResistance::new(resistance_pa_s2_m6),
        0.0,
        synthetic_source(&format!("conjugate-loss-{name}-v1")),
        ToleranceBasis::EngineeringAllowance,
    )
    .expect("synthetic quadratic loss")
}

fn airflow_network() -> EnclosureNetwork {
    EnclosureNetwork::new(
        LossNetwork::Element(loss(CHANNEL_BRANCH, 1.0e10)),
        LeakageElement::new(loss("declared-case-leakage", 1.6e11)),
    )
}

fn operating_point(speed_ratio: f64) -> OperatingPoint {
    let fan = FanBank::new(fan_curve(), 1, FanArrangement::Series, speed_ratio)
        .expect("speed inside the declared fan-law domain");
    solve_operating_point(&fan, &airflow_network()).expect("unique synthetic operating point")
}

fn flow_area_m2() -> f64 {
    core::f64::consts::PI * HYDRAULIC_DIAMETER_M * HYDRAULIC_DIAMETER_M / 4.0
}

fn bar_mesh() -> ConductionMesh {
    let (complex, positions) = box_grid(
        [CELLS_X, CELLS_Y, CELLS_Z],
        [BAR_LENGTH_M, BAR_WIDTH_M, BAR_HEIGHT_M],
    );
    let complex = TetComplex::from_tets(positions.len(), complex.tets);
    ConductionMesh::new(complex, positions).expect("bar conduction mesh")
}

fn region_name(index: usize) -> String {
    format!("duct-wall-{index}")
}

/// The cooled trace: the top face, cut into `SEGMENTS` stream-wise bands.
fn in_segment(face: &BoundaryFace, index: usize) -> bool {
    if face.outward_normal[2] < 0.99 || !on_box_face(face.centroid[2], BAR_HEIGHT_M) {
        return false;
    }
    let span = BAR_LENGTH_M / SEGMENTS as f64;
    let low = index as f64 * span;
    let high = low + span;
    face.centroid[0] > low && face.centroid[0] < high
}

fn thermal_boundary(
    mesh: &ConductionMesh,
    htc_w_m2_k: f64,
    references_k: &[f64],
) -> ThermalBoundary {
    let mut builder = ThermalBoundaryBuilder::new(mesh);
    for (index, &reference) in references_k.iter().enumerate() {
        builder = builder
            .region(
                &region_name(index),
                |face| in_segment(face, index),
                ThermalBc::robin(htc_w_m2_k, reference).expect("admissible Robin row"),
            )
            .expect("named stream-wise Robin region");
    }
    builder
        .adiabatic_remainder()
        .finish()
        .expect("complete duct-wall boundary partition")
}

fn uniform_source(mesh: &ConductionMesh) -> ScalarField {
    let volume = BAR_LENGTH_M * BAR_WIDTH_M * BAR_HEIGHT_M;
    let density = TOTAL_POWER_W / volume;
    ScalarField::nodal(
        "conjugate bar source",
        mesh.vertex_count(),
        vec![density; mesh.vertex_count()],
    )
    .expect("uniform volumetric source")
}

fn solve_config() -> SolveConfig {
    SolveConfig {
        nonlinearity: Nonlinearity::FixedPoint {
            relaxation: 1.0,
            max_backtracks: 8,
        },
        stop: StopRule {
            residual_rtol: 1.0e-10,
            residual_atol: 1.0e-20,
            step_atol: 0.0,
            max_iterations: 12,
        },
        linear: LinearConfig {
            // A purely-Robin domain has NO Dirichlet row pinning it, so the
            // achievable true relative residual floors out around 1e-12 and a
            // tighter gate would report "the mesh is unpinned" as "the solve
            // failed". This is the exact trap the solid-side prerequisite
            // (frankensim-xcbm1) measured; keep ~100x headroom over the floor.
            tolerance: 1.0e-10,
            max_iterations: 60_000,
            restart: 80,
        },
        initial: InitialGuess::Uniform(320.0),
    }
}

/// Everything the fixture derives from the solved operating point, before any
/// coupling happens.
struct DuctFixture {
    mesh: ConductionMesh,
    path: AirPath,
    htc_w_m2_k: f64,
    capacity_rate_w_per_k: f64,
    reynolds: f64,
}

fn duct_fixture(speed_ratio: f64) -> DuctFixture {
    let point = operating_point(speed_ratio);
    let handoff = point
        .correlation_handoff(
            CHANNEL_BRANCH,
            Area::new(flow_area_m2()),
            Density::new(AIR_DENSITY_KG_M3),
            DynViscosity::new(AIR_DYNAMIC_VISCOSITY_PA_S),
            Length::new(HYDRAULIC_DIAMETER_M),
            AIR_PRANDTL,
        )
        .expect("dimensioned duct handoff");

    let inputs = handoff
        .correlation_inputs
        .with_length_ratio(BAR_LENGTH_M / HYDRAULIC_DIAMETER_M);
    // The coupled wall is very nearly isothermal per segment, so the
    // constant-wall-temperature limit is the matched card. It refuses outside
    // its declared Re / Pr / (L/D_h) domain, so admission here is a gate the
    // fixture passes rather than an assumption it makes.
    let nusselt = evaluate(CorrelationId::CircularDuctLaminarCwt, inputs)
        .expect("fully-developed laminar CWT card admits this duct");
    assert!(nusselt.evidence().model.in_domain);

    let htc = nusselt
        .heat_transfer_coefficient(
            ThermalConductivity::new(AIR_THERMAL_CONDUCTIVITY_W_M_K),
            Length::new(HYDRAULIC_DIAMETER_M),
        )
        .expect("Nu to h")
        .value
        .value();

    let mass_flow = AIR_DENSITY_KG_M3 * handoff.branch_flow.value.value();
    let segment_area = BAR_LENGTH_M * BAR_WIDTH_M / SEGMENTS as f64;
    let segments = (0..SEGMENTS)
        .map(|index| {
            AirSegment::new(&region_name(index), segment_area, htc).expect("valid segment")
        })
        .collect();
    let path = AirPath::new(
        INLET_TEMPERATURE_K,
        mass_flow,
        AIR_SPECIFIC_HEAT_J_KG_K,
        segments,
    )
    .expect("valid duct air path");

    DuctFixture {
        mesh: bar_mesh(),
        capacity_rate_w_per_k: path.capacity_rate_w_per_k(),
        path,
        htc_w_m2_k: htc,
        reynolds: handoff.reynolds,
    }
}

/// The solid side of one exchange: rebuild the Robin rows against the current
/// air references, solve, and hand back the per-region decomposition.
fn solid_exchange(
    cx: &Cx<'_>,
    fixture: &DuctFixture,
    references_k: &[f64],
) -> (Vec<SolidRegionState>, f64) {
    let boundary = thermal_boundary(&fixture.mesh, fixture.htc_w_m2_k, references_k);
    let material = ConductivityModel::isotropic_declared(SOLID_CONDUCTIVITY_W_M_K)
        .expect("declared solid conductivity");
    let source = uniform_source(&fixture.mesh);
    let solution = solve(
        cx,
        ConductionProblem {
            mesh: &fixture.mesh,
            boundary: &boundary,
            material: &material,
            source: &source,
        },
        solve_config(),
    )
    .expect("conjugate inner conduction solve");
    let states = solution
        .report
        .robin_fluxes
        .iter()
        .map(SolidRegionState::from_robin_flux)
        .collect();
    (states, solution.report.energy.robin_out_w)
}

fn config() -> ConjugateConfig {
    ConjugateConfig {
        // Well below the ~1e-3 K resolution any downstream QoI would use, and
        // far above the inner solve's own noise floor.
        temperature_tolerance_k: 1.0e-9,
        max_iterations: 60,
        relaxation: Relaxation::Fixed { omega: 1.0 },
    }
}

#[test]
fn the_duct_fixture_sits_inside_every_declared_card_domain() {
    let fixture = duct_fixture(1.0);
    assert!(
        fixture.reynolds > 1.0 && fixture.reynolds < 2_300.0,
        "Re {} left the laminar card domain",
        fixture.reynolds
    );
    // L/D_h is a compile-time property of the fixture constants; see the
    // const assertion above.
    // The exchange is strongly coupled: NTU near 1 is where the air rise is a
    // first-order effect and the exponential reference law matters. A fixture
    // at NTU ≈ 0 would converge trivially and prove nothing.
    let total_ntu: f64 = fixture
        .path
        .segments()
        .iter()
        .map(|s| s.ntu(fixture.capacity_rate_w_per_k))
        .sum();
    assert!(
        (0.5..10.0).contains(&total_ntu),
        "fixture NTU {total_ntu} is outside the interesting band"
    );
}

#[test]
fn the_coupled_duct_converges_and_closes_its_interface_balance() {
    let fixture = duct_fixture(1.0);
    let mut robin_out_w = 0.0;
    let solution = with_cx(|cx| {
        solve_conjugate(cx, &fixture.path, &config(), |cx, references| {
            let (states, total) = solid_exchange(cx, &fixture, references);
            robin_out_w = total;
            Ok(states)
        })
        .expect("the coupled duct converges")
    })
    .with_decomposition_cross_check(robin_out_w);

    // Clause: converges with a logged iteration history.
    assert!(
        solution.iterations >= 2,
        "a one-iteration 'convergence' would mean the solid never responded"
    );
    assert_eq!(solution.history.len(), solution.iterations);
    assert!(
        solution
            .history
            .last()
            .expect("nonempty")
            .max_reference_change_k
            <= config().temperature_tolerance_k
    );

    // Clause: the flux-balance invariant holds. This is the WIRING falsifier
    // described in the module docs, not a conservation proof: at the fixed
    // point T_ref,eff makes it an identity for a uniform h.
    assert!(
        solution.balance.max_region_imbalance_w < 1.0e-6 * TOTAL_POWER_W.max(1.0),
        "per-region imbalance {} W",
        solution.balance.max_region_imbalance_w
    );

    // The decomposition cross-check: fs-conduction accumulated robin_out_w in
    // its own loop over every Robin face. If the four declared regions did not
    // partition the cooled trace, this disagrees even though the per-region
    // identity above still closes.
    let decomposition = solution
        .balance
        .decomposition_residual_w
        .expect("cross-check attached");
    assert!(
        decomposition.abs() < 1.0e-9,
        "declared regions do not partition the Robin boundary: {decomposition} W"
    );

    // Global energy: everything dissipated leaves through the air.
    assert!(
        (solution.balance.air_total_w - TOTAL_POWER_W).abs() < 1.0e-3 * TOTAL_POWER_W,
        "air carried {} W of {TOTAL_POWER_W} W dissipated",
        solution.balance.air_total_w
    );
    assert!(!solution.worst_recorded_imbalance_w.is_nan());
}

#[test]
fn stream_wise_air_heating_tilts_the_wall_temperature_downstream() {
    // The result the one-way model structurally cannot produce. With one
    // declared bulk air temperature every segment sees the same coolant, so a
    // uniformly heated symmetric bar comes out symmetric. Here the downstream
    // coolant is hotter, so the downstream wall must be hotter.
    let fixture = duct_fixture(1.0);
    let solution = with_cx(|cx| {
        solve_conjugate(cx, &fixture.path, &config(), |cx, references| {
            Ok(solid_exchange(cx, &fixture, references).0)
        })
        .expect("converges")
    });

    let walls: Vec<f64> = solution
        .solid
        .iter()
        .map(|state| state.mean_wall_temperature_k)
        .collect();
    for window in walls.windows(2) {
        assert!(
            window[1] > window[0],
            "wall temperature did not rise downstream: {walls:?}"
        );
    }
    // And the air rise is a first-order effect, not a rounding artifact.
    let rise = solution.march.outlet_temperature_k - INLET_TEMPERATURE_K;
    assert!(
        rise > 5.0,
        "air rise {rise} K is too small for this fixture to be a conjugate demonstration"
    );
    assert!(
        walls[SEGMENTS - 1] - walls[0] > 1.0,
        "stream-wise wall tilt {} K is too small to be distinguishable",
        walls[SEGMENTS - 1] - walls[0]
    );

    // Every segment's air stays below its own wall: heat flows the right way
    // along the whole path.
    for (segment, state) in solution.march.segments.iter().zip(&solution.solid) {
        assert!(segment.outlet_temperature_k < state.mean_wall_temperature_k);
        assert!(state.heat_rate_w > 0.0);
    }
}

#[test]
fn a_faster_fan_cools_the_duct_and_flattens_its_tilt() {
    // A metamorphic check across the whole chain: more flow means more
    // capacity rate, so both the peak wall temperature and the stream-wise
    // tilt must fall. This crosses fan curve, operating point, correlation
    // card, air path, and FEM solve, so it is sensitive to a sign or wiring
    // error anywhere in the chain.
    let run = |speed_ratio: f64| {
        let fixture = duct_fixture(speed_ratio);
        let solution = with_cx(|cx| {
            solve_conjugate(cx, &fixture.path, &config(), |cx, references| {
                Ok(solid_exchange(cx, &fixture, references).0)
            })
            .expect("converges")
        });
        let walls: Vec<f64> = solution
            .solid
            .iter()
            .map(|state| state.mean_wall_temperature_k)
            .collect();
        let peak = walls.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let tilt = walls[SEGMENTS - 1] - walls[0];
        let rise = solution.march.outlet_temperature_k - INLET_TEMPERATURE_K;
        (peak, tilt, rise)
    };

    let (slow_peak, slow_tilt, slow_rise) = run(0.8);
    let (fast_peak, fast_tilt, fast_rise) = run(1.2);
    assert!(
        fast_peak < slow_peak,
        "a faster fan must cool the duct: {fast_peak} vs {slow_peak}"
    );
    assert!(
        fast_tilt < slow_tilt,
        "a faster fan must flatten the stream-wise tilt: {fast_tilt} vs {slow_tilt}"
    );
    assert!(
        fast_rise < slow_rise,
        "a faster fan must reduce the air rise: {fast_rise} vs {slow_rise}"
    );
}

#[test]
fn the_coupled_exchange_replays_bitwise() {
    // G5: the driver has no RNG and no thread-order dependence, so two runs of
    // the same fixture must agree bit for bit, including the whole history.
    let fixture = duct_fixture(1.0);
    let run = || {
        with_cx(|cx| {
            solve_conjugate(cx, &fixture.path, &config(), |cx, references| {
                Ok(solid_exchange(cx, &fixture, references).0)
            })
            .expect("converges")
        })
    };
    let first = run();
    let second = run();
    assert_eq!(first.iterations, second.iterations);
    for (a, b) in first
        .reference_temperatures_k
        .iter()
        .zip(&second.reference_temperatures_k)
    {
        assert_eq!(a.to_bits(), b.to_bits());
    }
    for (a, b) in first.history.iter().zip(&second.history) {
        assert_eq!(
            a.max_reference_change_k.to_bits(),
            b.max_reference_change_k.to_bits()
        );
        assert_eq!(
            a.air_outlet_temperature_k.to_bits(),
            b.air_outlet_temperature_k.to_bits()
        );
        assert_eq!(a.solid_heat_rate_w.to_bits(), b.solid_heat_rate_w.to_bits());
    }
}

#[test]
fn a_declared_ambient_and_the_coupled_answer_disagree_materially() {
    // The falsifier for the whole bead: if the conjugate answer matched the
    // old one-way declared-ambient answer, none of this machinery would be
    // earning its place. Solve the same bar with ONE bulk reference — the
    // arithmetic mid-rise the previous fixture used — and show the peak wall
    // temperature differs by more than any tolerance in this file.
    let fixture = duct_fixture(1.0);
    let coupled = with_cx(|cx| {
        solve_conjugate(cx, &fixture.path, &config(), |cx, references| {
            Ok(solid_exchange(cx, &fixture, references).0)
        })
        .expect("converges")
    });
    let coupled_peak = coupled
        .solid
        .iter()
        .map(|state| state.mean_wall_temperature_k)
        .fold(f64::NEG_INFINITY, f64::max);

    // The one-way model: bulk air = inlet + half the sensible rise implied by
    // the DECLARED power, applied uniformly to every region.
    let declared_bulk =
        INLET_TEMPERATURE_K + TOTAL_POWER_W / (2.0 * fixture.capacity_rate_w_per_k);
    let one_way = with_cx(|cx| {
        solid_exchange(cx, &fixture, &[declared_bulk; SEGMENTS])
    });
    let one_way_peak = one_way
        .0
        .iter()
        .map(|state| state.mean_wall_temperature_k)
        .fold(f64::NEG_INFINITY, f64::max);

    assert!(
        (coupled_peak - one_way_peak).abs() > 1.0,
        "coupled peak {coupled_peak} K and declared-ambient peak {one_way_peak} K are \
         indistinguishable — this fixture would not demonstrate conjugate coupling"
    );
    // The one-way model under-predicts the hot end, because it gives the
    // downstream wall coolant that is colder than the air actually reaching it.
    assert!(
        coupled_peak > one_way_peak,
        "the declared-ambient model should be optimistic at the hot end"
    );
}

#[test]
#[ignore = "diagnostic: prints the fixture's operating numbers"]
fn diagnostic_fixture_numbers() {
    let fixture = duct_fixture(1.0);
    let total_ntu: f64 = fixture
        .path
        .segments()
        .iter()
        .map(|s| s.ntu(fixture.capacity_rate_w_per_k))
        .sum();
    let solution = with_cx(|cx| {
        solve_conjugate(cx, &fixture.path, &config(), |cx, references| {
            Ok(solid_exchange(cx, &fixture, references).0)
        })
        .expect("converges")
    });
    let walls: Vec<f64> = solution
        .solid
        .iter()
        .map(|s| s.mean_wall_temperature_k)
        .collect();
    let declared_bulk =
        INLET_TEMPERATURE_K + TOTAL_POWER_W / (2.0 * fixture.capacity_rate_w_per_k);
    let one_way = with_cx(|cx| {
        solid_exchange(cx, &fixture, &[declared_bulk; SEGMENTS])
    });
    let one_way_peak = one_way
        .0
        .iter()
        .map(|s| s.mean_wall_temperature_k)
        .fold(f64::NEG_INFINITY, f64::max);
    println!("Re                = {:.1}", fixture.reynolds);
    println!("h                 = {:.3} W/m2K", fixture.htc_w_m2_k);
    println!("mdot*cp           = {:.5} W/K", fixture.capacity_rate_w_per_k);
    println!("total NTU         = {total_ntu:.4}");
    println!("iterations        = {}", solution.iterations);
    println!("air rise          = {:.3} K", solution.march.outlet_temperature_k - INLET_TEMPERATURE_K);
    println!("walls             = {walls:?}");
    println!("wall tilt         = {:.3} K", walls[SEGMENTS - 1] - walls[0]);
    println!("coupled peak      = {:.3} K", walls.iter().copied().fold(f64::NEG_INFINITY, f64::max));
    println!("one-way peak      = {one_way_peak:.3} K");
    println!("air total         = {:.6} W (declared {TOTAL_POWER_W})", solution.balance.air_total_w);
    println!("max region imbal  = {:.3e} W", solution.balance.max_region_imbalance_w);
    println!("residual history  = {:?}", solution.history.iter().map(|h| h.max_reference_change_k).collect::<Vec<_>>());
}
