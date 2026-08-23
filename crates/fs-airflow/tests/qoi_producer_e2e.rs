//! No-mock thermal-QoI producer end-to-end lane (bead
//! `frankensim-extreal-program-f85xj.5.10`).
//!
//! Driven by `scripts/ci/thermal_qoi_producer_e2e.sh`. The battery runs the
//! REAL production path end to end -- a solved fan/system operating point and
//! a real `fs-conduction` FEM solve feed `extract_thermal_qois` and the
//! operating-envelope audit -- and emits a canonical, hash-chained JSONL
//! event stream. No evidence object is fabricated anywhere in this file:
//! every QoI value is produced by the public producer API from solver output.
//!
//! Modes (env-selected so plain `cargo test` stays green workspace-wide):
//! - `TQOI_MODE=produce`: run the case battery and emit events to
//!   `$TQOI_ARTIFACT_DIR/events.ndjson` (and, with `TQOI_EMIT_STDOUT=1`, as
//!   one base64 stdout line for remote-run retrieval).
//! - `TQOI_MODE=verify`: re-read the emitted stream and re-verify schema,
//!   strict ordinal order, hash-chain continuity (BLAKE3 via `fs-blake3`),
//!   required case coverage, and per-record term completeness. This is the
//!   solver-free checker half of the lane.
//!
//! Tamper evidence: every event carries `prev` (the previous event's digest)
//! and `digest` = BLAKE3 hex over the event's canonical bytes WITHOUT the
//! digest field. Any byte flip breaks the chain at verification. Integrity is
//! not authenticity: the chain detects edits AFTER production; it does not
//! mint authority.
//!
//! Every f64 crosses the artifact boundary as its IEEE bit pattern, so the
//! stream carries no float-formatting nondeterminism at all.

use std::fmt::Write as _;

use fs_airflow::qoi::{
    FanPowerSpec, JunctionRegion, QoiError, SafetyFactorAuthority, SurfaceRegion,
    ThermalOutputAuditError, ThermalQoi, ThermalQoiCardUse, ThermalQoiDeclarations,
    ThermalRequirement, extract_thermal_qois,
};
use fs_airflow::{
    EnclosureNetwork, FanArrangement, FanBank, FanCurve, FanPoint, LeakageElement, LossElement,
    LossNetwork, LossResistance, SourceProvenance, ToleranceBasis, solve_operating_point,
};
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::Blake3;
use fs_conduction::ThermalBoundary;
use fs_conduction::bc::{ThermalBc, ThermalBoundaryBuilder};
use fs_conduction::fixtures::unit_cube;
use fs_conduction::{
    ConductionError, ConductionMesh, ConductionProblem, ConductionSolution, ConductivityModel,
    InitialGuess, LinearConfig, LinearSolveEvidence, Nonlinearity, ResidualClaim, ScalarField,
    SolveConfig, StopRule, solve,
};
use fs_evidence::ModelCard;
use fs_evidence::uncertainty::{EngineeringUncertaintyKind, TermValue};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_qty::{Pressure, Temperature, VolumetricFlowRate};

const SCHEMA: &str = "org.frankensim.fs-airflow.qoi-producer-e2e.v1";

/// One deterministic-repeat arm's full output tuple.
type RepeatJourney = Result<(Vec<f64>, [f64; 7], [String; 2]), String>;
/// When set, `produce` additionally prints the whole canonical stream as one
/// base64 line prefixed with this marker, so a remote runner (RCH) can
/// retrieve the artifact from captured stdout alone.
const STDOUT_MARKER: &str = "TQOI_EVENTS_B64:";

/// Minimal standard base64 (with padding).
fn base64_encode(input: &[u8]) -> String {
    const ALPHABET: &[u8; 64] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    let mut out = String::with_capacity(input.len().div_ceil(3) * 4);
    for chunk in input.chunks(3) {
        let byte0 = u32::from(chunk[0]);
        let byte1 = u32::from(chunk.get(1).copied().unwrap_or(0));
        let byte2 = u32::from(chunk.get(2).copied().unwrap_or(0));
        let word = (byte0 << 16) | (byte1 << 8) | byte2;
        out.push(ALPHABET[(word >> 18) as usize & 63] as char);
        out.push(ALPHABET[(word >> 12) as usize & 63] as char);
        out.push(if chunk.len() > 1 {
            ALPHABET[(word >> 6) as usize & 63] as char
        } else {
            '='
        });
        out.push(if chunk.len() > 2 {
            ALPHABET[word as usize & 63] as char
        } else {
            '='
        });
    }
    out
}

/// Case coverage demanded per profile; mirrored by the driving script.
const CASES_PR: &[&str] = &[
    "five_families_happy",
    "absent_requirement_refuses",
    "tampered_conduction_refuses",
];
const CASES_FULL: &[&str] = &[
    "five_families_happy",
    "absent_requirement_refuses",
    "narrowed_card_demotes_all",
    "tampered_conduction_refuses",
    "cancellation_refuses_before_production",
    "deterministic_repeat_bitwise",
];
const CASES_RECOVERY: &[&str] = &["cancellation_refuses_before_production"];

// ------------------------------------------------------------------ emitter

/// Canonical JSONL emitter with a BLAKE3 hash chain.
struct Emitter {
    body: String,
    ordinal: u64,
    prev: Option<String>,
}

impl Emitter {
    fn new() -> Self {
        Self {
            body: String::new(),
            ordinal: 0,
            prev: None,
        }
    }

    /// Emit one event. Values must already be JSON-safe.
    fn emit(&mut self, event: &str, fields: &[(&str, String)]) {
        let mut core = String::new();
        let prev_json = match &self.prev {
            Some(hash) => format!("\"{hash}\""),
            None => "null".to_string(),
        };
        let _ = write!(
            core,
            "{{\"schema\":\"{SCHEMA}\",\"event\":\"{event}\",\"ordinal\":{},\"prev\":{prev_json}",
            self.ordinal
        );
        for (key, value) in fields {
            let _ = write!(core, ",\"{key}\":{value}");
        }
        core.push('}');

        let mut hasher = Blake3::new();
        hasher.update(core.as_bytes());
        let digest = hasher.finalize().to_hex();

        // Insert the digest before the closing brace; the on-disk line stays
        // valid JSON with deterministic field order, and verification can
        // strip the suffix to recompute the digest from the core bytes.
        self.body.push_str(&core[..core.len() - 1]);
        let _ = writeln!(self.body, ",\"digest\":\"{digest}\"}}");

        self.prev = Some(digest);
        self.ordinal += 1;
    }

    fn write_to(self, directory: &std::path::Path) -> Result<(), String> {
        std::fs::write(directory.join("events.ndjson"), &self.body)
            .map_err(|error| format!("write events.ndjson: {error}"))
    }
}

fn json_str(value: &str) -> String {
    let mut out = String::with_capacity(value.len() + 2);
    out.push('"');
    for character in value.chars() {
        match character {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            control if (control as u32) < 0x20 => {
                let _ = write!(out, "\\u{:04x}", control as u32);
            }
            other => out.push(other),
        }
    }
    out.push('"');
    out
}

/// f64 values cross the artifact boundary as IEEE bit patterns.
fn bits(value: f64) -> String {
    format!("\"{:016x}\"", value.to_bits())
}

fn bits_list(values: &[f64]) -> String {
    let parts: Vec<String> = values.iter().map(|value| bits(*value)).collect();
    format!("[{}]", parts.join(","))
}
fn term_value_json(value: &TermValue) -> String {
    match value {
        TermValue::IntervalBound { lower, upper } => format!(
            "{{\"state\":\"interval\",\"lower\":{},\"upper\":{}}}",
            bits(*lower),
            bits(*upper)
        ),
        TermValue::Unknown { reason } => {
            format!("{{\"state\":\"unknown\",\"reason\":{}}}", json_str(reason))
        }
        TermValue::Negligible { justification } => {
            format!(
                "{{\"state\":\"negligible\",\"why\":{}}}",
                json_str(justification)
            )
        }
        TermValue::Distribution(term) => format!(
            "{{\"state\":\"distribution\",\"mean\":{},\"std\":{},\"half_width\":{},\"level\":{}}}",
            bits(term.mean),
            bits(term.standard_deviation),
            bits(term.conservative_half_width),
            bits(term.level)
        ),
        TermValue::Ensemble(term) => format!(
            "{{\"state\":\"ensemble\",\"members\":{},\"half_width\":{}}}",
            term.member_count,
            bits(term.conservative_half_width)
        ),
        TermValue::CorrelatedBlock(_) => "{\"state\":\"correlated\"}".to_string(),
        _ => "{\"state\":\"other\"}".to_string(),
    }
}

fn refusal_json(code: &str, error: &dyn std::fmt::Display) -> String {
    format!(
        "{{\"code\":{},\"message\":{}}}",
        json_str(code),
        json_str(&error.to_string())
    )
}

// ----------------------------------------------------------------- fixtures

fn source(id: &str) -> SourceProvenance {
    SourceProvenance::new("retained synthetic producer-e2e source", id)
}

fn fan_curve() -> FanCurve {
    FanCurve::new(
        "producer-e2e-fan",
        vec![
            FanPoint::new(VolumetricFlowRate::new(0.0), Pressure::new(160.0)),
            FanPoint::new(VolumetricFlowRate::new(0.04), Pressure::new(130.0)),
            FanPoint::new(VolumetricFlowRate::new(0.08), Pressure::new(70.0)),
            FanPoint::new(VolumetricFlowRate::new(0.12), Pressure::new(0.0)),
        ],
        source("producer-e2e-fan-v1"),
        0.08,
        ToleranceBasis::EngineeringAllowance,
        VolumetricFlowRate::new(0.01),
        (0.7, 1.3),
    )
    .expect("valid fan fixture")
}

fn loss(name: &str, resistance: f64, uncertainty: f64) -> LossElement {
    LossElement::new(
        name,
        LossResistance::new(resistance),
        uncertainty,
        source(&format!("producer-e2e-loss-{name}")),
        ToleranceBasis::EngineeringAllowance,
    )
    .expect("valid loss fixture")
}

fn network() -> EnclosureNetwork {
    let primary = LossNetwork::series(vec![
        LossNetwork::Element(loss("inlet", 40_000.0, 0.10)),
        LossNetwork::Element(loss("heatsink", 30_000.0, 0.12)),
        LossNetwork::Element(loss("outlet", 12_000.0, 0.08)),
    ])
    .expect("series network");
    EnclosureNetwork::new(
        primary,
        LeakageElement::new(loss("leakage", 180_000.0, 0.25)),
    )
}

fn linear_config() -> LinearConfig {
    LinearConfig {
        tolerance: 1.0e-10,
        max_iterations: 60_000,
        restart: 80,
    }
}

fn solve_config() -> SolveConfig {
    SolveConfig {
        initial: InitialGuess::Uniform(320.0),
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
        linear: linear_config(),
    }
}

fn with_cx_gate<R>(gate: &CancelGate, seed: u64, kernel: u64, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            gate,
            arena,
            StreamKey {
                seed,
                kernel_id: kernel,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn ambient_boundary(mesh: &ConductionMesh) -> ThermalBoundary {
    ThermalBoundaryBuilder::new(mesh)
        .region(
            "ambient",
            |_| true,
            ThermalBc::dirichlet(300.0).expect("ambient"),
        )
        .expect("all-face ambient region")
        .finish()
        .expect("boundary")
}

/// A real FEM conduction solve on the unit-cube tet mesh: uniform declared
/// conductivity, uniform volumetric source, ambient Dirichlet on every
/// boundary face. The returned report comes from the production solver.
fn solve_real_conduction() -> Result<(ConductionMesh, ConductionSolution), ConductionError> {
    let (complex, positions) = unit_cube(1);
    let mesh = ConductionMesh::new(complex, positions).expect("unit cube mesh");
    let material = ConductivityModel::isotropic_declared(205.0).expect("declared aluminium");
    let heat_source = ScalarField::uniform("heater", 5.0e6).expect("bounded uniform source");
    let boundary = ambient_boundary(&mesh);
    with_cx_gate(&CancelGate::new(), 0x0000_C0DE_0510_0001, 71, |cx| {
        solve(
            cx,
            ConductionProblem {
                mesh: &mesh,
                boundary: &boundary,
                material: &material,
                source: &heat_source,
            },
            solve_config(),
        )
        .map(|solution| (mesh, solution))
    })
}

fn declarations(mesh: &ConductionMesh) -> (JunctionRegion, SurfaceRegion, FanPowerSpec) {
    let junction = JunctionRegion::try_new("package", vec![7, 0, 6]).expect("junction region");
    let surface =
        SurfaceRegion::try_new("case", (0..mesh.boundary().len()).rev().collect::<Vec<_>>())
            .expect("surface region");
    let power = FanPowerSpec::try_new(0.72, 0.04, source("producer-e2e-efficiency-v1"))
        .expect("efficiency");
    (junction, surface, power)
}

fn requirement_at(effective_limit_k: f64, factor: f64) -> ThermalRequirement {
    ThermalRequirement::try_new(
        Temperature::new(effective_limit_k),
        SafetyFactorAuthority::try_new(factor, source("derating-policy-v1")).expect("factor"),
        source("component-datasheet-limit-v1"),
    )
    .expect("requirement")
}

fn extract_parts(
    mesh: &ConductionMesh,
    solution: &ConductionSolution,
    operating: &fs_airflow::OperatingPoint,
    requirement: Option<&ThermalRequirement>,
) -> Result<fs_airflow::qoi::ThermalQoiSet, QoiError> {
    let (junction, surface, power) = declarations(mesh);
    extract_thermal_qois(
        mesh,
        solution,
        operating,
        &ThermalQoiDeclarations {
            junction_region: &junction,
            surface_region: &surface,
            fan_power: &power,
            requirement,
            discretization: None,
        },
    )
}

fn regime_point(flow_m3_s: f64) -> fs_regime::OperatingPoint {
    fs_regime::OperatingPoint {
        id: "solved".to_string(),
        groups: [
            ("flow_m3_s".to_string(), flow_m3_s),
            ("speed_ratio".to_string(), 1.0),
        ]
        .into_iter()
        .collect(),
    }
}

fn card_uses(
    qois: &fs_airflow::qoi::ThermalQoiSet,
    cards: &[&ModelCard],
) -> Vec<ThermalQoiCardUse> {
    qois.budgets()
        .into_iter()
        .map(|budget| ThermalQoiCardUse {
            qoi: budget.qoi().to_string(),
            model_cards: cards.iter().map(|card| card.name.clone()).collect(),
            override_acknowledgement: None,
        })
        .collect()
}
/// Local accessor trait so the generic emitter can read any `fs-qty`
/// dimensioned value without depending on a crate trait that does not exist.
trait HasValue {
    fn qty_value(&self) -> f64;
}

impl<const M: i8, const KG: i8, const S: i8, const K: i8, const A: i8, const MOL: i8> HasValue
    for fs_qty::Qty<M, KG, S, K, A, MOL>
{
    fn qty_value(&self) -> f64 {
        fs_qty::Qty::value(*self)
    }
}

/// Emit one qoi_record event carrying its complete eight-term budget state.
fn emit_qoi<T: HasValue>(emitter: &mut Emitter, name: &str, unit: &str, qoi: &ThermalQoi<T>) {
    let budget = &qoi.uncertainty;
    let mut terms = String::new();
    for (ordinal, kind) in EngineeringUncertaintyKind::ALL.iter().enumerate() {
        if ordinal > 0 {
            terms.push(',');
        }
        let _ = write!(
            terms,
            "{{\"kind\":{},\"term\":{}}}",
            json_str(kind.name()),
            term_value_json(budget.term(*kind).value())
        );
    }
    let numerical_kind = match qoi.evidence.numerical.kind {
        fs_evidence::NumericalKind::NoClaim => "no_claim".to_string(),
        other => format!("{other:?}"),
    };
    let identity = budget.content_id();
    emitter.emit(
        "qoi_record",
        &[
            ("qoi", json_str(name)),
            ("identity", json_str(&identity.to_hex())),
            ("numerical_kind", json_str(&numerical_kind)),
            ("unit", json_str(unit)),
            ("value_bits", bits(qoi.evidence.value.qty_value())),
            ("terms", format!("[{terms}]")),
        ],
    );
}

/// Emit the raw inputs the script-side independent oracle recomputes from.
fn emit_inputs(
    emitter: &mut Emitter,
    temperatures: &[f64],
    operating: &fs_airflow::OperatingPoint,
) {
    let fan_points: Vec<String> = fan_curve()
        .points()
        .iter()
        .map(|point| {
            format!(
                "{{\"q\":{},\"p\":{}}}",
                bits(point.flow.value()),
                bits(point.pressure.value())
            )
        })
        .collect();
    let branch_flows: Vec<String> = operating
        .branches
        .iter()
        .map(|branch| {
            format!(
                "{{\"path\":{},\"q\":{}}}",
                json_str(&branch.path),
                bits(branch.flow.value.value())
            )
        })
        .collect();
    emitter.emit(
        "input_artifacts",
        &[
            ("temperature_bits", bits_list(temperatures)),
            ("junction_vertices", "[7,0,6]".to_string()),
            ("limit_k_bits", bits(380.0)),
            ("efficiency_bits", bits(0.72)),
            (
                "resistances",
                bits_list(&[40_000.0, 30_000.0, 12_000.0, 180_000.0]),
            ),
            ("fan_points", format!("[{}]", fan_points.join(","))),
            ("operating_flow_bits", bits(operating.flow.value.value())),
            (
                "operating_pressure_bits",
                bits(operating.pressure.value.value()),
            ),
            ("branch_flows", format!("[{}]", branch_flows.join(","))),
            ("leakage_fraction_bits", bits(operating.leakage_fraction)),
        ],
    );
}

// ------------------------------------------------------------------- modes

fn expected_cases(profile: &str) -> &'static [&'static str] {
    match profile {
        "pr" => CASES_PR,
        "full" => CASES_FULL,
        "recovery" => CASES_RECOVERY,
        other => panic!("unknown TQOI_PROFILE {other}"),
    }
}

#[test]
fn producer_e2e_lane() {
    let Ok(mode) = std::env::var("TQOI_MODE") else {
        // Workspace-suite default: the lane only runs under its driver.
        eprintln!("qoi_producer_e2e: skipped (set TQOI_MODE=produce|verify)");
        return;
    };
    let directory = std::env::var("TQOI_ARTIFACT_DIR")
        .map(std::path::PathBuf::from)
        .expect("TQOI_ARTIFACT_DIR must be set for produce/verify modes");
    match mode.as_str() {
        "produce" => produce(&directory).expect("produce battery failed"),
        "verify" => verify(&directory).expect("verify failed"),
        other => panic!("unknown TQOI_MODE {other}"),
    }
}

fn solve_network_operating_point() -> fs_airflow::OperatingPoint {
    let fan = FanBank::new(fan_curve(), 1, FanArrangement::Series, 1.0).expect("fan bank");
    solve_operating_point(&fan, &network()).expect("operating point")
}

#[allow(clippy::too_many_lines)]
fn produce(directory: &std::path::Path) -> Result<(), String> {
    let profile = std::env::var("TQOI_PROFILE").unwrap_or_else(|_| "pr".to_string());
    let cases = expected_cases(&profile);
    let mut emitter = Emitter::new();

    emitter.emit(
        "run_begin",
        &[
            ("profile", json_str(&profile)),
            (
                "source_revision",
                json_str(
                    std::env::var("TQOI_SOURCE_REVISION")
                        .unwrap_or_else(|_| "unknown".into())
                        .as_str(),
                ),
            ),
        ],
    );

    let requires = |name: &str| cases.contains(&name);

    // ------------------------------------------------------ happy five-set
    if requires("five_families_happy") {
        let operating = solve_network_operating_point();
        let (mesh, solution) = solve_real_conduction().map_err(|error| error.to_string())?;
        emit_inputs(&mut emitter, &solution.temperature, &operating);

        let set = extract_parts(
            &mesh,
            &solution,
            &operating,
            Some(&requirement_at(380.0, 1.25)),
        )
        .map_err(|error| format!("happy extraction refused: {error}"))?;

        emit_qoi(
            &mut emitter,
            "junction_maximum",
            "kelvin",
            &set.junction_maximum.qoi,
        );
        emit_qoi(&mut emitter, "pressure_drop", "pascal", &set.pressure_drop);
        emit_qoi(&mut emitter, "fan_power", "watt", &set.fan_power);
        emit_qoi(
            &mut emitter,
            "uniformity_mean",
            "kelvin",
            &set.uniformity.mean_temperature,
        );
        emit_qoi(
            &mut emitter,
            "uniformity_spread",
            "kelvin",
            &set.uniformity.spread,
        );
        emit_qoi(
            &mut emitter,
            "face_mean_std",
            "kelvin",
            &set.uniformity.face_mean_standard_deviation,
        );
        emit_qoi(
            &mut emitter,
            "thermal_margin",
            "kelvin",
            &set.thermal_margin,
        );

        // Envelope audit against the fan's own published model card at the
        // SOLVED operating point.
        let card = fan_curve().model_card();
        let point = regime_point(operating.flow.value.value());
        let uses = card_uses(&set, &[&card]);
        let audited = set
            .clone()
            .audit_operating_envelope(
                std::slice::from_ref(&card),
                std::slice::from_ref(&point),
                &uses,
            )
            .map_err(|error: ThermalOutputAuditError| format!("envelope audit refused: {error}"))?;
        emitter.emit(
            "envelope_audit",
            &[
                ("receipts", audited.audit.receipts.len().to_string()),
                (
                    "first_coverage",
                    json_str(&format!("{:?}", audited.audit.receipts[0].coverage)),
                ),
                (
                    "any_demoted",
                    audited
                        .audit
                        .receipts
                        .iter()
                        .any(fs_regime::OutputClaimReceipt::demoted)
                        .to_string(),
                ),
            ],
        );
    }

    // --------------------------------------------- absent requirement refuses
    if requires("absent_requirement_refuses") {
        let operating = solve_network_operating_point();
        let (mesh, solution) = solve_real_conduction().map_err(|error| error.to_string())?;
        let error = extract_parts(&mesh, &solution, &operating, None)
            .expect_err("margin without a declared requirement must refuse");
        emitter.emit(
            "refusal",
            &[
                ("case", json_str("absent_requirement")),
                ("refusal", refusal_json("qoi-missing-requirement", &error)),
            ],
        );
    }

    // ------------------------------------------- narrowed card demotes all
    if requires("narrowed_card_demotes_all") {
        let operating = solve_network_operating_point();
        let (mesh, solution) = solve_real_conduction().map_err(|error| error.to_string())?;
        let set = extract_parts(
            &mesh,
            &solution,
            &operating,
            Some(&requirement_at(380.0, 1.25)),
        )
        .map_err(|error| format!("pre-narrow extraction refused: {error}"))?;
        let card = fan_curve().model_card();
        let mut narrowed = card.clone();
        narrowed.validity = narrowed.validity.with("flow_m3_s", 0.0, 0.005);
        let point = regime_point(operating.flow.value.value());

        let wide_uses = card_uses(&set, &[&card]);
        let narrow_uses = card_uses(&set, &[&narrowed]);
        let wide = set
            .clone()
            .audit_operating_envelope(
                std::slice::from_ref(&card),
                std::slice::from_ref(&point),
                &wide_uses,
            )
            .map_err(|error: ThermalOutputAuditError| format!("wide audit refused: {error}"))?;
        let tight = set
            .clone()
            .audit_operating_envelope(
                std::slice::from_ref(&narrowed),
                std::slice::from_ref(&point),
                &narrow_uses,
            )
            .map_err(|error: ThermalOutputAuditError| format!("narrow audit refused: {error}"))?;

        let wide_in = wide.audit.receipts.iter().all(|receipt| !receipt.demoted());
        let all_demoted = tight
            .audit
            .receipts
            .iter()
            .all(fs_regime::OutputClaimReceipt::demoted);
        if !wide_in || !all_demoted {
            return Err(format!(
                "narrowing did not demote: wide_in={wide_in} all_demoted={all_demoted}"
            ));
        }
        emitter.emit(
            "demotion",
            &[
                ("case", json_str("narrowed_card_demotes_all")),
                (
                    "wide_coverage",
                    json_str(&format!("{:?}", wide.audit.receipts[0].coverage)),
                ),
                (
                    "narrow_coverage",
                    json_str(&format!("{:?}", tight.audit.receipts[0].coverage)),
                ),
                ("receipt_count", tight.audit.receipts.len().to_string()),
            ],
        );
    }

    // ------------------------------------- tampered conduction report refuses
    if requires("tampered_conduction_refuses") {
        let operating = solve_network_operating_point();
        let (mesh, solution) = solve_real_conduction().map_err(|error| error.to_string())?;

        // (a) receipts claimed but none retained: self-contradictory.
        let mut contradictory = solution.clone();
        contradictory.report.material_receipts = 0;
        let contradiction = extract_parts(
            &mesh,
            &contradictory,
            &operating,
            Some(&requirement_at(380.0, 1.25)),
        )
        .err()
        .ok_or("zero-receipt contradiction was accepted")?;
        if !matches!(contradiction, QoiError::InvalidInput { .. }) {
            return Err(format!(
                "unexpected contradiction variant {contradiction:?}"
            ));
        }

        // (b) an algebraically unconverged linear record is unsupported.
        let mut unconverged = solution.clone();
        unconverged.report.linear.push(LinearSolveEvidence {
            nonlinear_iteration: 1,
            method: "pcg",
            iterations: 500,
            reported: ResidualClaim::RecursiveEstimate(1.0e-11),
            true_relative_residual: 3.2e-2,
            converged_true: false,
            stall: None,
        });
        let unconv_error = extract_parts(
            &mesh,
            &unconverged,
            &operating,
            Some(&requirement_at(380.0, 1.25)),
        )
        .err()
        .ok_or("unconverged linear record was accepted")?;
        if !matches!(unconv_error, QoiError::InvalidInput { .. }) {
            return Err(format!("unexpected unconverged variant {unconv_error:?}"));
        }

        emitter.emit(
            "refusal",
            &[
                ("case", json_str("tampered_conduction")),
                ("refusal", refusal_json("qoi-invalid-input", &contradiction)),
                (
                    "refusal_second",
                    refusal_json("qoi-invalid-input", &unconv_error),
                ),
            ],
        );
    }

    // ------------------------------- cancellation refuses before production
    if requires("cancellation_refuses_before_production") {
        let (complex, positions) = unit_cube(1);
        let mesh = ConductionMesh::new(complex, positions).expect("unit cube mesh");
        let material = ConductivityModel::isotropic_declared(205.0).expect("material");
        let heat_source = ScalarField::uniform("heater", 5.0e6).expect("uniform source");
        let boundary = ambient_boundary(&mesh);
        let gate = CancelGate::new();
        gate.request();
        let outcome = with_cx_gate(&gate, 0x0000_C0DE_0510_0002, 72, |cx| {
            solve(
                cx,
                ConductionProblem {
                    mesh: &mesh,
                    boundary: &boundary,
                    material: &material,
                    source: &heat_source,
                },
                solve_config(),
            )
        });
        let cancelled = outcome
            .err()
            .ok_or("pre-cancelled solve completed; the cancellation drill is vacuous")?;
        let code = match &cancelled {
            ConductionError::Cancelled { stage, .. } => format!("conduction-cancelled:{stage}"),
            other => return Err(format!("expected cancellation, got {other:?}")),
        };
        emitter.emit(
            "refusal",
            &[
                ("case", json_str("cancellation")),
                (
                    "refusal",
                    format!("{{\"code\":{},\"message\":null}}", json_str(&code)),
                ),
            ],
        );
    }

    // ------------------------------------------ deterministic repeat bitwise
    if requires("deterministic_repeat_bitwise") {
        let journey_once = || -> RepeatJourney {
            let operating = solve_network_operating_point();
            let (mesh, solution) = solve_real_conduction().map_err(|error| error.to_string())?;
            let set = extract_parts(
                &mesh,
                &solution,
                &operating,
                Some(&requirement_at(380.0, 1.25)),
            )
            .map_err(|error| format!("repeat arm refused: {error}"))?;
            Ok((
                solution.temperature.clone(),
                [
                    set.junction_maximum.qoi.evidence.value.value(),
                    set.pressure_drop.evidence.value.value(),
                    set.fan_power.evidence.value.value(),
                    set.uniformity.mean_temperature.evidence.value.value(),
                    set.uniformity.spread.evidence.value.value(),
                    set.uniformity
                        .face_mean_standard_deviation
                        .evidence
                        .value
                        .value(),
                    set.thermal_margin.evidence.value.value(),
                ],
                [
                    set.junction_maximum.qoi.uncertainty.content_id().to_hex(),
                    set.thermal_margin.uncertainty.content_id().to_hex(),
                ],
            ))
        };
        let first = journey_once()?;
        let second = journey_once()?;
        for (ordinal, (left, right)) in first.1.iter().zip(second.1.iter()).enumerate() {
            if left.to_bits() != right.to_bits() {
                return Err(format!("repeat arm {ordinal} diverged bitwise"));
            }
        }
        if first.2 != second.2 {
            return Err("repeat rebound a budget identity".to_string());
        }
        emitter.emit(
            "determinism",
            &[
                ("case", json_str("deterministic_repeat_bitwise")),
                ("values_compared", first.1.len().to_string()),
                ("identities_compared", first.2.len().to_string()),
                ("temperature_root_bits", bits_list(&first.0)),
            ],
        );
    }

    emitter.emit(
        "run_end",
        &[
            ("cases", cases.len().to_string()),
            ("ok", "true".to_string()),
        ],
    );
    let stdout_copy = if std::env::var("TQOI_EMIT_STDOUT").as_deref() == Ok("1") {
        Some(base64_encode(emitter.body.as_bytes()))
    } else {
        None
    };
    emitter.write_to(directory)?;
    if let Some(payload) = stdout_copy {
        println!("{STDOUT_MARKER}{payload}");
    }
    Ok(())
}

/// Solver-free integrity check of an emitted stream: schema constant, strict
/// ordinal order, hash-chain continuity, required case coverage, and
/// completeness of the seven-record happy set.
fn verify(directory: &std::path::Path) -> Result<(), String> {
    let profile = std::env::var("TQOI_PROFILE").unwrap_or_else(|_| "pr".to_string());
    let cases = expected_cases(&profile);
    let body = std::fs::read_to_string(directory.join("events.ndjson"))
        .map_err(|error| format!("read events: {error}"))?;
    let lines: Vec<&str> = body.lines().filter(|line| !line.is_empty()).collect();
    if lines.is_empty() {
        return Err("empty event stream".to_string());
    }

    let mut previous: Option<String> = None;
    let mut seen_cases: Vec<String> = Vec::new();
    let mut ended_ok = false;

    for (index, line) in lines.iter().enumerate() {
        let (core, digest) = verify_event_integrity(index, line, &mut previous)?;
        if core.contains("\"case\":\"") {
            let marker = "\"case\":\"";
            let start = core.find(marker).expect("checked") + marker.len();
            let name: String = core[start..].chars().take_while(|c| *c != '"').collect();
            seen_cases.push(name);
        }
        if core.contains("\"event\":\"run_end\"") && core.contains("\"ok\":true") {
            ended_ok = true;
        }
        previous = Some(digest);
    }
    if !ended_ok {
        return Err("stream has no successful run_end".to_string());
    }
    for required in cases {
        match *required {
            "five_families_happy" => {
                let records = lines
                    .iter()
                    .filter(|line| line.contains("\"event\":\"qoi_record\""))
                    .count();
                if records < 7 {
                    return Err("happy case emitted fewer than seven qoi_record events".into());
                }
                if !lines
                    .iter()
                    .any(|line| line.contains("\"event\":\"envelope_audit\""))
                {
                    return Err("happy case emitted no envelope_audit event".into());
                }
            }
            "narrowed_card_demotes_all" => {
                if !lines
                    .iter()
                    .any(|line| line.contains("\"event\":\"demotion\""))
                {
                    return Err("narrowing case emitted no demotion event".into());
                }
            }
            "deterministic_repeat_bitwise" => {
                if !lines
                    .iter()
                    .any(|line| line.contains("\"event\":\"determinism\""))
                {
                    return Err("repeat case emitted no determinism event".into());
                }
            }
            other => {
                if !seen_cases.iter().any(|seen| seen == other) {
                    return Err(format!("required case {other} absent from stream"));
                }
            }
        }
    }
    Ok(())
}

/// Recompute one event's BLAKE3 digest over its core bytes, verify the
/// prev-chain linkage and schema constant, and check ordinal order.
/// Returns the core (digest-stripped) bytes and the verified digest.
fn verify_event_integrity(
    index: usize,
    line: &str,
    previous: &mut Option<String>,
) -> Result<(String, String), String> {
    let digest_marker = ",\"digest\":\"";
    let split = line
        .rfind(digest_marker)
        .ok_or_else(|| format!("event {index}: no digest field"))?;
    let core = &line[..split];
    let digest = &line[split + digest_marker.len()..line.len() - 2];

    let mut hasher = Blake3::new();
    hasher.update(core.as_bytes());
    let recomputed = hasher.finalize().to_hex();
    if recomputed != digest {
        return Err(format!(
            "event {index}: digest mismatch (chain broken or tampered)"
        ));
    }

    let prev_marker = ",\"prev\":";
    let prev_start = core
        .find(prev_marker)
        .ok_or_else(|| format!("event {index}: no prev field"))?
        + prev_marker.len();
    let prev_tail = &core[prev_start..];
    let prev_value = if let Some(rest) = prev_tail.strip_prefix('"') {
        rest.split('"').next().unwrap_or_default().to_string()
    } else {
        prev_tail
            .split([',', '}'])
            .next()
            .unwrap_or_default()
            .trim()
            .to_string()
    };
    match (&*previous, prev_value.as_str()) {
        (None, "null") => {}
        (Some(expected), actual) if expected == actual => {}
        _ => return Err(format!("event {index}: prev chain mismatch")),
    }

    if !core.contains(&format!("\"schema\":\"{SCHEMA}\"")) {
        return Err(format!("event {index}: wrong schema"));
    }
    let ordinal_marker = "\"ordinal\":";
    let ordinal_start = core.find(ordinal_marker).ok_or("no ordinal")? + ordinal_marker.len();
    let ordinal: u64 = core[ordinal_start..]
        .split([',', '}'])
        .next()
        .and_then(|raw| raw.parse().ok())
        .ok_or_else(|| format!("event {index}: bad ordinal"))?;
    if ordinal != index as u64 {
        return Err(format!("event {index}: ordinal {ordinal} out of order"));
    }
    Ok((core.to_string(), digest.to_string()))
}
