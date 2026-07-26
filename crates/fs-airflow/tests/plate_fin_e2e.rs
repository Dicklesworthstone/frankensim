//! G0/G3 end-to-end plate-fin airflow-to-conduction wiring fixture.
//!
//! The fixture is intentionally synthetic. It proves that a fan-curve
//! operating point drives a dimensioned rectangular-channel handoff, that
//! validity-gated convection rows reach separately named heatsink surfaces,
//! and that declared base-to-fin bondlines participate in the actual
//! conduction solve. It does not validate a commercial fan, a plate-fin
//! array, commercial hardware, or conjugate fluid physics. Its nominal
//! fin-flank coefficient executes one narrow table slice for a smooth
//! simultaneously developing rectangular duct; the wrong-`D_h` falsifier
//! deliberately enters the separately marked engineering bridge.

use fs_airflow::{
    EnclosureNetwork, FanArrangement, FanBank, FanCurve, FanPoint, LeakageElement, LossElement,
    LossNetwork, LossResistance, SourceProvenance, ToleranceBasis, solve_operating_point,
};
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_conduction::fixtures::{box_grid, on_box_face};
use fs_conduction::{
    AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS, AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY,
    ConductionMesh, ConductionProblem, ConductionSolution, ConductivityModel, InitialGuess,
    InterfaceFacePair, InterfaceResistance, InterfaceSurface, LinearConfig, Nonlinearity,
    ScalarField, SolveConfig, StopRule, ThermalBc, ThermalBoundary, ThermalBoundaryBuilder,
    ThermalInterfaces, solve_with_interfaces,
};
use fs_convection::{
    CorrelationError, CorrelationId, EvaluationSourceRegion, ThermalConductivity, evaluate,
};
use fs_evidence::{ProvenanceHash, ValidityDomain};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance as MaterialProvenance, QueryPoint, SelectionPolicy,
    SurfaceSpec, SystemContext, UncertaintyModel,
};
use fs_qty::{Area, Density, DynViscosity, Length, Pressure, Temperature, VolumetricFlowRate};
use fs_rep_mesh::TetComplex;

const GRID_PITCH_M: f64 = 1.0 / 2_048.0;
// 174 pitches = 87/1024 m, exactly representable. This length is chosen to
// maximise the SMALLER of the fixture's two validity margins; see
// FULLY_DEVELOPED_FLOOR and GRAETZ_TABLE_FLOOR below. The former 5/64 put
// L/D_h at exactly 50.0, sitting on a bound one ULP from refusal.
const CHANNEL_LENGTH_M: f64 = 174.0 * GRID_PITCH_M;
const BASE_WIDTH_M: f64 = 9.0 * GRID_PITCH_M;
const BASE_THICKNESS_M: f64 = GRID_PITCH_M;
const FIN_THICKNESS_M: f64 = GRID_PITCH_M;
const FIN_HEIGHT_M: f64 = 4.0 * GRID_PITCH_M;
const FIN_Y_START_CELLS: [usize; 3] = [1, 4, 7];
const FIN_COUNT: usize = FIN_Y_START_CELLS.len();
const INTERNAL_CHANNEL_COUNT: usize = FIN_COUNT - 1;
const CELLS_X: usize = 4;
const CELLS_FIN_Z: usize = 2;

const SOLID_CONDUCTIVITY_W_M_K: f64 = 205.0;
const BASE_SOURCE_W_M3: f64 = 5.0e6;
const BONDLINE_RESISTANCE_M2_K_W: f64 = 2.0e-5;

const AIR_DENSITY_KG_M3: f64 = 1.2;
const AIR_DYNAMIC_VISCOSITY_PA_S: f64 = 1.8e-5;
const AIR_THERMAL_CONDUCTIVITY_W_M_K: f64 = 0.026;
const AIR_SPECIFIC_HEAT_J_KG_K: f64 = 1_005.0;
const AIR_PRANDTL: f64 = 0.72;
const INLET_TEMPERATURE_K: f64 = 300.0;

const CHANNEL_BRANCH: &str = "fin-channels";
const SPEED_RATIOS: [f64; 3] = [0.75, 1.0, 1.25];

/// The `L_over_Dh` lower bound the fully-developed CWT/CHF cards declare
/// (`fs-convection`'s `fully_developed()` domain). It is INCLUSIVE, so a
/// fixture sitting exactly on it is admitted -- and is one ULP from refusal.
const FULLY_DEVELOPED_FLOOR: f64 = 50.0;
/// The Graetz number below which the developing card leaves its Shah-London
/// table slice for the declared engineering bridge
/// (`RECTANGULAR_DEVELOPING_FIRST_TABLE_GRAETZ`).
const GRAETZ_TABLE_FLOOR: f64 = 10.0;
/// How far above BOTH floors this fixture is required to sit.
///
/// The two bounds pull in opposite directions -- `Gz = Re·Pr/(L/D_h)`, so
/// lengthening the channel buys margin on the fully-developed floor and
/// spends it on the Graetz floor. With this fixture's slowest solved
/// Reynolds number the admissible window is roughly `L/D_h ∈ [50, 59.3]`,
/// and the equal-margin optimum is `sqrt(Re·Pr·50/10) ≈ 54.5`. The chosen
/// length lands at `L/D_h = 54.375`: 8.75% clear of the floor and 9.15%
/// clear of the table cliff.
const VALIDITY_MARGIN: f64 = 0.08;

/// Boundary triangles each named Robin trace must own, derived from the mesh
/// parameters rather than pinned as a magic number. Existence-only coverage
/// (`faces > 0`) is not enough: a predicate that silently matched HALF the fin
/// flanks passes every other assertion in this file, because the energy
/// closure compares the region sum against a total that drops by the same
/// amount and no absolute temperature is pinned.
const CHANNEL_GAP_CELLS: usize = 2;
const EXPECTED_FLANK_FACES: usize = INTERNAL_CHANNEL_COUNT * 2 * CELLS_X * CELLS_FIN_Z * 2;
const EXPECTED_FLOOR_FACES: usize = INTERNAL_CHANNEL_COUNT * CHANNEL_GAP_CELLS * CELLS_X * 2;

#[derive(Debug, Clone, Copy)]
struct ChannelGeometry {
    flow_area_m2: f64,
    wetted_perimeter_m: f64,
    hydraulic_diameter_m: f64,
    aspect_ratio: f64,
}

impl ChannelGeometry {
    fn from_mesh_constants() -> Self {
        let gap = 2.0 * GRID_PITCH_M;
        let one_channel_area = gap * FIN_HEIGHT_M;
        // A plate-fin passage is open at the top in this fixture. Its wetted
        // perimeter is the base floor plus two fin flanks, not the perimeter
        // of an invented shroud.
        //
        // MODEL MISMATCH, stated rather than left for a reader to derive: this
        // D_h describes a THREE-SIDED OPEN channel, while every card it feeds
        // is declared for a FOUR-WALL CLOSED duct with its condition on all
        // four walls. The `aspect_ratio` handed to those cards, 0.5, is the
        // closed 2p x 4p rectangle -- whose own D_h would be 2.667p, not the
        // 3.2p computed here. So the (D_h, aspect_ratio) PAIR describes no
        // single duct. The fixture additionally applies CWT to the flanks and
        // CHF to the floor of the SAME passage, which are two mutually
        // exclusive four-wall idealizations held at once.
        //
        // That is defensible for a WIRING fixture -- it exercises the solved
        // flow -> Re -> card -> h -> Robin -> solve chain, and the falsifier
        // tests below move real solved quantities. It is NOT a claim that
        // these h values describe this passage, and nothing here may be read
        // as validating the correlation against plate-fin hardware.
        let one_channel_wetted_perimeter = gap + 2.0 * FIN_HEIGHT_M;
        let flow_area_m2 = INTERNAL_CHANNEL_COUNT as f64 * one_channel_area;
        let wetted_perimeter_m = INTERNAL_CHANNEL_COUNT as f64 * one_channel_wetted_perimeter;
        Self {
            flow_area_m2,
            wetted_perimeter_m,
            hydraulic_diameter_m: 4.0 * flow_area_m2 / wetted_perimeter_m,
            aspect_ratio: gap / FIN_HEIGHT_M,
        }
    }
}

struct PlateFinMesh {
    mesh: ConductionMesh,
    base_vertex_count: usize,
}

#[derive(Debug)]
struct PlateFinRun {
    speed_ratio: f64,
    branch_flow_m3_s: f64,
    reynolds: f64,
    developing_cwt_provenance: ProvenanceHash,
    developing_cwt_source_region: EvaluationSourceRegion,
    chf_provenance: ProvenanceHash,
    developing_cwt_nusselt: f64,
    developing_cwt_h_w_m2_k: f64,
    cwt_limit_h_w_m2_k: f64,
    chf_h_w_m2_k: f64,
    mean_bulk_temperature_k: f64,
    mean_base_temperature_k: f64,
    hydraulic_diameter_m: f64,
    solution: ConductionSolution,
}

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0x0000_A1F1_0A7E_0000,
                kernel_id: 63,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
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
            // The slender, disconnected blocks joined by finite bondline
            // conductance leave this operator badly conditioned, and the
            // achievable true relative residual is SENSITIVE TO THE GEOMETRY:
            // measured 1.98e-11 / 1.82e-11 / 1.85e-11 at speed ratios
            // 0.75 / 1.00 / 1.25 on the current mesh, where the previous
            // 5/64 channel reached ~1.1e-11.
            //
            // The old 2.0e-11 gate was set just above that older floor, which
            // left the slowest run 1.2% under its own fail-closed threshold.
            // Perturbing the geometry by as little as 0.1% pushed it over
            // (observed 2.007e-11) -- so the gate was reporting "the mesh
            // changed" as "the solve failed". Keep roughly 5x headroom over
            // the measured worst case instead; this is still two orders
            // tighter than the 1e-8 energy-closure assertion it feeds.
            tolerance: 1.0e-10,
            max_iterations: 60_000,
            restart: 80,
        },
        initial: InitialGuess::Uniform(320.0),
    }
}

fn synthetic_source(identifier: &str) -> SourceProvenance {
    SourceProvenance::new(
        "Synthetic plate-fin G0/G3 fixture; not manufacturer performance data",
        identifier,
    )
}

fn fan_curve() -> FanCurve {
    FanCurve::new(
        "synthetic-plate-fin-fan",
        vec![
            FanPoint::new(VolumetricFlowRate::new(0.0), Pressure::new(30.0)),
            FanPoint::new(VolumetricFlowRate::new(30.0e-6), Pressure::new(25.0)),
            FanPoint::new(VolumetricFlowRate::new(60.0e-6), Pressure::new(12.0)),
            FanPoint::new(VolumetricFlowRate::new(90.0e-6), Pressure::new(0.0)),
        ],
        synthetic_source("plate-fin-fan-v1"),
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
        synthetic_source(&format!("plate-fin-loss-{name}-v1")),
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

fn operating_point(speed_ratio: f64) -> fs_airflow::OperatingPoint {
    let fan = FanBank::new(fan_curve(), 1, FanArrangement::Series, speed_ratio)
        .expect("speed is inside the declared fan-law domain");
    solve_operating_point(&fan, &airflow_network()).expect("unique synthetic operating point")
}

fn plate_fin_mesh() -> PlateFinMesh {
    let (base, mut positions) = box_grid(
        [CELLS_X, 9, 1],
        [CHANNEL_LENGTH_M, BASE_WIDTH_M, BASE_THICKNESS_M],
    );
    let base_vertex_count = positions.len();
    let mut tets = base.tets;

    for start_cell in FIN_Y_START_CELLS {
        let (fin, fin_positions) = box_grid(
            [CELLS_X, 1, CELLS_FIN_Z],
            [CHANNEL_LENGTH_M, FIN_THICKNESS_M, FIN_HEIGHT_M],
        );
        let offset = u32::try_from(positions.len()).expect("fixture vertex count fits u32");
        tets.extend(
            fin.tets
                .into_iter()
                .map(|tet| tet.map(|vertex| vertex + offset)),
        );
        let y_offset = start_cell as f64 * GRID_PITCH_M;
        positions.extend(
            fin_positions
                .into_iter()
                .map(|[x, y, z]| [x, y + y_offset, z + BASE_THICKNESS_M]),
        );
    }

    let complex = TetComplex::from_tets(positions.len(), tets);
    PlateFinMesh {
        mesh: ConductionMesh::new(complex, positions).expect("plate-fin conduction mesh"),
        base_vertex_count,
    }
}

fn bond_surface(material: &str, texture: &str) -> SurfaceSpec {
    SurfaceSpec {
        material: MaterialStateId {
            chemistry: material.to_string(),
            phase: "solid".to_string(),
            process: "as-fixtured".to_string(),
            revision: 0,
        },
        texture_frame: texture.to_string(),
    }
}

fn bondline_card() -> InterfaceSystemCard {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(
                AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY,
                AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS,
            ),
            value: PropertyValue::Scalar {
                value: BONDLINE_RESISTANCE_M2_K_W,
                dims: AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS,
            },
            validity: ValidityDomain::unconstrained(),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: Vec::new(),
            provenance: MaterialProvenance {
                source: "synthetic plate-fin bondline fixture".to_string(),
                license: "internal-test-use".to_string(),
                artifact: None,
            },
        })
        .expect("bondline claim inserts");
    InterfaceSystemCard::assemble(
        bond_surface("aluminium-base", "bond-normal-plus-z"),
        bond_surface("aluminium-fin", "bond-normal-minus-z"),
        SystemContext {
            medium: "dry".to_string(),
            third_body: Some("declared-synthetic-bondline".to_string()),
            environment: "enclosure-air".to_string(),
            history: "unaged".to_string(),
        },
        claims,
        Vec::new(),
    )
    .expect("bondline card assembles")
}

fn fin_index_for_y(y: f64) -> usize {
    FIN_Y_START_CELLS
        .iter()
        .position(|&start| {
            let low = start as f64 * GRID_PITCH_M;
            let high = low + FIN_THICKNESS_M;
            y >= low && y <= high
        })
        .expect("every coincident face belongs to one fin footprint")
}

fn thermal_interfaces(mesh: &ConductionMesh, boundary: &ThermalBoundary) -> ThermalInterfaces {
    let card = bondline_card();
    let resistance = InterfaceResistance::from_card(
        "plate-fin-bondlines",
        &card,
        &QueryPoint::new(),
        SelectionPolicy::SingleClaimOnly,
    )
    .expect("card-backed bondline resistance");
    let mut per_fin = vec![Vec::<InterfaceFacePair>::new(); FIN_COUNT];

    for pair in ThermalInterfaces::coincident_face_pairs(mesh).expect("coincident fin footprints") {
        let side_a_is_base = mesh.boundary()[pair.side_a].outward_normal[2] > 0.0;
        let oriented = if side_a_is_base {
            pair
        } else {
            InterfaceFacePair {
                side_a: pair.side_b,
                side_b: pair.side_a,
            }
        };
        assert!(
            mesh.boundary()[oriented.side_a].outward_normal[2] > 0.99
                && mesh.boundary()[oriented.side_b].outward_normal[2] < -0.99,
            "bondline orientation must run from base to fin"
        );
        let y = mesh.boundary()[oriented.side_a].centroid[1];
        per_fin[fin_index_for_y(y)].push(oriented);
    }

    let surfaces = per_fin
        .into_iter()
        .enumerate()
        .map(|(fin, pairs)| {
            assert_eq!(
                pairs.len(),
                2 * CELLS_X,
                "one matching triangle pair for each fin-footprint triangle"
            );
            InterfaceSurface::new(format!("fin-{fin}-bondline"), pairs, resistance.clone())
                .expect("declared fin bondline")
        })
        .collect();
    ThermalInterfaces::new(mesh, boundary, surfaces).expect("complete bondline binding")
}

fn is_internal_fin_flank(face: &fs_conduction::BoundaryFace) -> bool {
    if face.centroid[2] <= BASE_THICKNESS_M || face.outward_normal[1].abs() < 0.99 {
        return false;
    }
    let internal_flank_y = [
        2.0 * GRID_PITCH_M,
        4.0 * GRID_PITCH_M,
        5.0 * GRID_PITCH_M,
        7.0 * GRID_PITCH_M,
    ];
    internal_flank_y
        .into_iter()
        .any(|target| on_box_face(face.centroid[1], target))
}

fn is_internal_channel_floor(face: &fs_conduction::BoundaryFace) -> bool {
    if face.outward_normal[2] < 0.99 || !on_box_face(face.centroid[2], BASE_THICKNESS_M) {
        return false;
    }
    let y = face.centroid[1];
    (y > 2.0 * GRID_PITCH_M && y < 4.0 * GRID_PITCH_M)
        || (y > 5.0 * GRID_PITCH_M && y < 7.0 * GRID_PITCH_M)
}

fn thermal_boundary(
    mesh: &ConductionMesh,
    fin_flank: ThermalBc,
    channel_floor: ThermalBc,
) -> ThermalBoundary {
    ThermalBoundaryBuilder::new(mesh)
        .region("fin-flanks", is_internal_fin_flank, fin_flank)
        .expect("named fin-flank Robin region")
        .region(
            "base-channel-floors",
            is_internal_channel_floor,
            channel_floor,
        )
        .expect("named base-floor Robin region")
        .adiabatic_remainder()
        .finish()
        .expect("complete plate-fin boundary partition")
}

fn base_source(mesh: &PlateFinMesh) -> ScalarField {
    let mut values = vec![0.0; mesh.mesh.vertex_count()];
    values[..mesh.base_vertex_count].fill(BASE_SOURCE_W_M3);
    ScalarField::nodal("plate-fin base source", mesh.mesh.vertex_count(), values)
        .expect("base-only source")
}

fn declared_base_power_w() -> f64 {
    BASE_SOURCE_W_M3 * CHANNEL_LENGTH_M * BASE_WIDTH_M * BASE_THICKNESS_M
}

fn solve_plate_fin(
    mesh: &PlateFinMesh,
    geometry: ChannelGeometry,
    speed_ratio: f64,
    hydraulic_diameter_scale: f64,
) -> PlateFinRun {
    let point = operating_point(speed_ratio);
    let hydraulic_diameter_m = geometry.hydraulic_diameter_m * hydraulic_diameter_scale;
    let handoff = point
        .correlation_handoff(
            CHANNEL_BRANCH,
            Area::new(geometry.flow_area_m2),
            Density::new(AIR_DENSITY_KG_M3),
            DynViscosity::new(AIR_DYNAMIC_VISCOSITY_PA_S),
            Length::new(hydraulic_diameter_m),
            AIR_PRANDTL,
        )
        .expect("dimensioned channel-flow handoff");
    let correlation_inputs = handoff
        .correlation_inputs
        .with_aspect_ratio(geometry.aspect_ratio)
        .with_length_ratio(CHANNEL_LENGTH_M / hydraulic_diameter_m);
    let cwt_limit = evaluate(CorrelationId::RectangularDuctLaminarCwt, correlation_inputs)
        .expect("constant-wall-temperature rectangular limit admits");
    let developing_cwt = evaluate(
        CorrelationId::RectangularDuctLaminarCwtDevelopingPr072,
        correlation_inputs,
    )
    .expect("developing constant-wall-temperature table/bridge card admits");
    let chf_limit = evaluate(CorrelationId::RectangularDuctLaminarChf, correlation_inputs)
        .expect("constant-heat-flux rectangular limit admits");

    assert!(cwt_limit.evidence().model.in_domain);
    assert!(developing_cwt.evidence().model.in_domain);
    assert!(chf_limit.evidence().model.in_domain);
    let branch_flow_m3_s = handoff.branch_flow.value.value();
    // One-way fixture closure: for a declared base power and constant cp, the
    // channel-mean bulk temperature is the midpoint of the inlet/outlet
    // sensible-heat rise. This is not a conjugate fluid solve.
    let mean_bulk_temperature_k = INLET_TEMPERATURE_K
        + declared_base_power_w()
            / (2.0 * AIR_DENSITY_KG_M3 * AIR_SPECIFIC_HEAT_J_KG_K * branch_flow_m3_s);
    let fin_coupling = developing_cwt
        .robin_boundary(
            ThermalConductivity::new(AIR_THERMAL_CONDUCTIVITY_W_M_K),
            Length::new(hydraulic_diameter_m),
            Temperature::new(mean_bulk_temperature_k),
        )
        .expect("fin-flank Robin lowering");
    let base_coupling = chf_limit
        .robin_boundary(
            ThermalConductivity::new(AIR_THERMAL_CONDUCTIVITY_W_M_K),
            Length::new(hydraulic_diameter_m),
            Temperature::new(mean_bulk_temperature_k),
        )
        .expect("base-floor Robin lowering");
    let developing_cwt_h_w_m2_k = fin_coupling.coefficient().value.value();
    let cwt_limit_h_w_m2_k = cwt_limit
        .heat_transfer_coefficient(
            ThermalConductivity::new(AIR_THERMAL_CONDUCTIVITY_W_M_K),
            Length::new(hydraulic_diameter_m),
        )
        .expect("CWT limit coefficient")
        .value
        .value();
    let chf_h_w_m2_k = base_coupling.coefficient().value.value();
    let boundary = thermal_boundary(
        &mesh.mesh,
        fin_coupling.boundary_condition().clone(),
        base_coupling.boundary_condition().clone(),
    );
    let interfaces = thermal_interfaces(&mesh.mesh, &boundary);
    let material = ConductivityModel::isotropic_declared(SOLID_CONDUCTIVITY_W_M_K)
        .expect("declared aluminium conductivity");
    let source = base_source(mesh);
    let solution = with_cx(|cx| {
        solve_with_interfaces(
            cx,
            ConductionProblem {
                mesh: &mesh.mesh,
                boundary: &boundary,
                material: &material,
                source: &source,
            },
            &interfaces,
            solve_config(),
        )
        .expect("plate-fin conduction solve")
    });
    let mean_base_temperature_k = solution.temperature[..mesh.base_vertex_count]
        .iter()
        .sum::<f64>()
        / mesh.base_vertex_count as f64;

    PlateFinRun {
        speed_ratio,
        branch_flow_m3_s,
        reynolds: handoff.reynolds,
        developing_cwt_provenance: developing_cwt.evidence().provenance,
        developing_cwt_source_region: developing_cwt.source_region(),
        chf_provenance: chf_limit.evidence().provenance,
        developing_cwt_nusselt: developing_cwt.evidence().value,
        developing_cwt_h_w_m2_k,
        cwt_limit_h_w_m2_k,
        chf_h_w_m2_k,
        mean_bulk_temperature_k,
        mean_base_temperature_k,
        hydraulic_diameter_m,
        solution,
    }
}

#[test]
fn solved_fan_rates_drive_a_declared_plate_fin_thermal_chain() {
    let mesh = plate_fin_mesh();
    let geometry = ChannelGeometry::from_mesh_constants();
    assert_eq!(geometry.aspect_ratio.to_bits(), 0.5f64.to_bits());

    // The fully-developed CWT/CHF cards bound L/D_h at 50.0 INCLUSIVE, and
    // `evaluate().expect()` already enforces it. Asserting `>= 50.0` here
    // therefore proved nothing the solve did not already prove -- and the
    // fixture used to sit on that bound bit-exactly, so one ULP of drift in
    // D_h turned three passing tests into three panics. Assert the MARGIN
    // instead, which is the property that actually keeps the fixture alive.
    let length_over_hydraulic_diameter = CHANNEL_LENGTH_M / geometry.hydraulic_diameter_m;
    assert!(
        length_over_hydraulic_diameter >= FULLY_DEVELOPED_FLOOR * (1.0 + VALIDITY_MARGIN),
        "L/D_h = {length_over_hydraulic_diameter} leaves no margin above the \
         fully-developed floor {FULLY_DEVELOPED_FLOOR}"
    );

    let runs: Vec<_> = SPEED_RATIOS
        .into_iter()
        .map(|speed| solve_plate_fin(&mesh, geometry, speed, 1.0))
        .collect();

    assert!(runs.iter().all(|run| {
        run.developing_cwt_source_region == EvaluationSourceRegion::TableSliceInterpolation
    }));

    // The other side of the squeeze. The source-region assertion above already
    // FAILS if Gz crosses the table floor, but it fails AT the cliff; this
    // fails while there is still room to move, which is the difference between
    // a diagnosable fixture and a surprising one. The slowest run is the
    // binding case because Gz rises with Reynolds.
    let slowest_graetz = runs
        .iter()
        .map(|run| run.reynolds * AIR_PRANDTL / length_over_hydraulic_diameter)
        .fold(f64::INFINITY, f64::min);
    assert!(
        slowest_graetz >= GRAETZ_TABLE_FLOOR * (1.0 + VALIDITY_MARGIN),
        "slowest solved Gz = {slowest_graetz} leaves no margin above the \
         table-slice floor {GRAETZ_TABLE_FLOOR}"
    );

    for pair in runs.windows(2) {
        assert!(pair[1].speed_ratio > pair[0].speed_ratio);
        assert!(pair[1].branch_flow_m3_s > pair[0].branch_flow_m3_s);
        assert!(pair[1].reynolds > pair[0].reynolds);
        assert_ne!(
            pair[1].developing_cwt_provenance,
            pair[0].developing_cwt_provenance
        );
        assert_ne!(pair[1].chf_provenance, pair[0].chf_provenance);
        assert!(
            pair[1].developing_cwt_nusselt > pair[0].developing_cwt_nusselt,
            "the table-slice interpolation must increase Nu with solved Reynolds"
        );
        assert!(
            pair[1].developing_cwt_h_w_m2_k > pair[0].developing_cwt_h_w_m2_k,
            "fixed geometry and fluid properties make the developing Nu change reach h"
        );
        assert!(
            pair[1].mean_bulk_temperature_k < pair[0].mean_bulk_temperature_k,
            "more solved mass flow lowers the declared channel-mean bulk temperature"
        );
        assert!(
            pair[1].mean_base_temperature_k < pair[0].mean_base_temperature_k,
            "the solved airflow change must reach the solid temperature field"
        );

        // The assertion above is NECESSARY but not sufficient, and on its own
        // it wears physics phrasing over a t_ref-wiring check: both Robin rows
        // take `mean_bulk_temperature_k` as their reference, so a base
        // temperature that merely tracked the declared bulk shift would
        // satisfy it while h did nothing at all.
        //
        // Separate the two paths. The excess of the base over the bulk is the
        // part the CONVECTIVE COEFFICIENT owns: it depends on the flank h and
        // not on the reference temperature. Requiring it to fall is a claim
        // about the developing card reaching the solid, and it is falsifiable
        // in exactly the case the phrasing above is not -- freeze h and this
        // difference becomes constant, because the problem is linear in
        // `t_ref` and the mesh, power, and floor card are all unchanged.
        let previous_excess_k = pair[0].mean_base_temperature_k - pair[0].mean_bulk_temperature_k;
        let next_excess_k = pair[1].mean_base_temperature_k - pair[1].mean_bulk_temperature_k;
        assert!(
            next_excess_k < previous_excess_k,
            "a higher solved flank h must shrink the base-over-bulk excess: \
             {previous_excess_k} K -> {next_excess_k} K"
        );
        assert!(
            previous_excess_k > 0.0 && next_excess_k > 0.0,
            "a heated base sits above its coolant: {previous_excess_k} K, {next_excess_k} K"
        );

        // The companion Shah-London CWT/CHF rows are fully-developed LIMITS.
        // Reynolds is only a validity axis for those rows, so they provide a
        // falsifier against accidentally attributing the developing-card
        // behavior to the shared Nu-to-h lowering.
        assert_eq!(
            pair[1].cwt_limit_h_w_m2_k.to_bits(),
            pair[0].cwt_limit_h_w_m2_k.to_bits()
        );
        assert_eq!(
            pair[1].chf_h_w_m2_k.to_bits(),
            pair[0].chf_h_w_m2_k.to_bits()
        );
    }

    for run in &runs {
        assert_eq!(run.solution.report.interface_fluxes.len(), FIN_COUNT);
        for (fin, flux) in run.solution.report.interface_fluxes.iter().enumerate() {
            assert_eq!(flux.interface, format!("fin-{fin}-bondline"));
            assert!(
                flux.heat_rate_a_to_b_w > 0.0,
                "base-to-fin bondline must conduct positive heat"
            );
        }

        let rows = &run.solution.report.robin_fluxes;
        assert_eq!(rows.len(), 2);
        assert_eq!(rows[0].region, "fin-flanks");
        assert_eq!(rows[1].region, "base-channel-floors");
        assert!(rows.iter().all(|row| row.heat_rate_w > 0.0));

        // Coverage by COUNT, not existence. See EXPECTED_*_FACES above for
        // why `faces > 0` cannot catch a half-matching region predicate.
        assert_eq!(
            rows[0].faces, EXPECTED_FLANK_FACES,
            "fin-flank trace must own every flank triangle"
        );
        assert_eq!(
            rows[1].faces, EXPECTED_FLOOR_FACES,
            "channel-floor trace must own every floor triangle"
        );

        // The binding that D_h's derivation actually needs: the hand-computed
        // wetted perimeter must be the perimeter the MESH exposes to the air.
        // Comparing D_h against `4A/P` instead only restates the initializer
        // in `from_mesh_constants` -- same tokens, same operands -- and cannot
        // fail. This can: it fails if either region predicate misses a face,
        // if the perimeter formula counts a surface the mesh does not have,
        // or if the channel length and the mesh extent disagree.
        let meshed_wetted_area_m2 = rows[0].area_m2 + rows[1].area_m2;
        let declared_wetted_area_m2 = CHANNEL_LENGTH_M * geometry.wetted_perimeter_m;
        let wetted_area_error =
            (meshed_wetted_area_m2 - declared_wetted_area_m2).abs() / declared_wetted_area_m2;
        assert!(
            wetted_area_error < 1.0e-12,
            "meshed wetted area {meshed_wetted_area_m2} m² vs the perimeter D_h was \
             built from {declared_wetted_area_m2} m², relative {wetted_area_error}"
        );
        // The report reconstructs each uniform h as Σ(A_f h)/ΣA_f. Bound
        // the two reductions and final division by a small multiple of the
        // participating face count instead of requiring a bitwise round trip.
        let developing_cwt_htc_error = (rows[0].mean_htc_w_per_m2_k - run.developing_cwt_h_w_m2_k)
            .abs()
            / run.developing_cwt_h_w_m2_k;
        let developing_cwt_reduction_bound = 4.0 * rows[0].faces as f64 * f64::EPSILON;
        assert!(
            developing_cwt_htc_error <= developing_cwt_reduction_bound,
            "fin-flank mean h relative error {developing_cwt_htc_error} exceeds \
             {developing_cwt_reduction_bound}"
        );
        let chf_htc_error =
            (rows[1].mean_htc_w_per_m2_k - run.chf_h_w_m2_k).abs() / run.chf_h_w_m2_k;
        let chf_reduction_bound = 4.0 * rows[1].faces as f64 * f64::EPSILON;
        assert!(
            chf_htc_error <= chf_reduction_bound,
            "channel-floor mean h relative error {chf_htc_error} exceeds {chf_reduction_bound}"
        );

        // This sum is deliberately only a WIRING falsifier: it catches a
        // dropped or misnamed region, but both values come from the same
        // boundary quadrature and therefore do not independently validate
        // the convection physics.
        let attributed_robin_w: f64 = rows.iter().map(|row| row.heat_rate_w).sum();
        let robin_total_w = run.solution.report.energy.robin_out_w;
        let attribution_error = (attributed_robin_w - robin_total_w).abs() / robin_total_w.abs();
        assert!(
            attribution_error < 1.0e-12,
            "named traces {attributed_robin_w} W vs total {robin_total_w} W"
        );
        let source_error = (run.solution.report.energy.source_w - declared_base_power_w()).abs()
            / declared_base_power_w();
        assert!(source_error < 1.0e-12);
        assert!(run.solution.report.energy.relative_closure() < 1.0e-8);
    }
}

#[test]
fn a_hydraulic_diameter_misbinding_changes_h_and_the_solved_temperature() {
    let mesh = plate_fin_mesh();
    let geometry = ChannelGeometry::from_mesh_constants();
    let nominal = solve_plate_fin(&mesh, geometry, 1.0, 1.0);
    let wrong = solve_plate_fin(&mesh, geometry, 1.0, 0.5);

    assert_eq!(
        nominal.branch_flow_m3_s.to_bits(),
        wrong.branch_flow_m3_s.to_bits(),
        "the falsifier changes only the downstream geometry handoff"
    );
    assert_eq!(
        nominal.mean_bulk_temperature_k.to_bits(),
        wrong.mean_bulk_temperature_k.to_bits()
    );
    assert_eq!(
        wrong.hydraulic_diameter_m.to_bits(),
        (0.5 * nominal.hydraulic_diameter_m).to_bits()
    );
    assert_eq!(
        nominal.developing_cwt_source_region,
        EvaluationSourceRegion::TableSliceInterpolation
    );
    assert_eq!(
        wrong.developing_cwt_source_region,
        EvaluationSourceRegion::DeclaredEngineeringBridge
    );
    assert_eq!(
        wrong.cwt_limit_h_w_m2_k.to_bits(),
        (2.0 * nominal.cwt_limit_h_w_m2_k).to_bits()
    );
    assert_eq!(
        wrong.chf_h_w_m2_k.to_bits(),
        (2.0 * nominal.chf_h_w_m2_k).to_bits()
    );
    assert!(
        wrong.developing_cwt_nusselt < nominal.developing_cwt_nusselt,
        "halving D_h must reduce Gz and the developing Nusselt number"
    );
    assert!(
        wrong.developing_cwt_h_w_m2_k > nominal.developing_cwt_h_w_m2_k
            && wrong.developing_cwt_h_w_m2_k < 2.0 * nominal.developing_cwt_h_w_m2_k,
        "the wrong-D_h declared-bridge coefficient must change through both Nu(Gz) and 1/D_h"
    );
    assert_ne!(
        wrong.developing_cwt_provenance,
        nominal.developing_cwt_provenance
    );
    assert!(
        wrong.mean_base_temperature_k < nominal.mean_base_temperature_k - 1.0,
        "wrong D_h must materially move the answer: nominal={} K wrong={} K",
        nominal.mean_base_temperature_k,
        wrong.mean_base_temperature_k
    );
}

#[test]
fn a_fan_solved_out_of_domain_flow_refuses_instead_of_extrapolating() {
    let geometry = ChannelGeometry::from_mesh_constants();
    let point = operating_point(2.75);
    let handoff = point
        .correlation_handoff(
            CHANNEL_BRANCH,
            Area::new(geometry.flow_area_m2),
            Density::new(AIR_DENSITY_KG_M3),
            DynViscosity::new(AIR_DYNAMIC_VISCOSITY_PA_S),
            Length::new(geometry.hydraulic_diameter_m),
            AIR_PRANDTL,
        )
        .expect("high-speed fan point still has a dimensioned handoff");
    assert!(
        handoff.reynolds > 2_300.0,
        "fixture must cross the laminar card ceiling"
    );
    let inputs = handoff
        .correlation_inputs
        .with_aspect_ratio(geometry.aspect_ratio)
        .with_length_ratio(CHANNEL_LENGTH_M / geometry.hydraulic_diameter_m);

    for correlation in [
        CorrelationId::RectangularDuctLaminarCwt,
        CorrelationId::RectangularDuctLaminarChf,
        CorrelationId::RectangularDuctLaminarCwtDevelopingPr072,
    ] {
        let error = evaluate(correlation, inputs).expect_err("laminar extrapolation must refuse");
        assert!(
            matches!(error, CorrelationError::OutOfDomain { .. }),
            "expected structured domain refusal, got {error:?}"
        );
        if let CorrelationError::OutOfDomain {
            correlation: refused,
            violations,
        } = error
        {
            assert_eq!(refused, correlation);
            assert!(
                violations.iter().any(|violation| {
                    violation.axis == "Re" && violation.value.is_some_and(|re| re > violation.high)
                }),
                "refusal must retain the solved Reynolds violation: {violations:?}"
            );
        }
    }
}
