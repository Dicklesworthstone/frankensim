//! The production cooling fidelity-graph instance (bead f85xj.10.4).
//!
//! This module assembles THE cooling vertical's fidelity graph from real,
//! carded model surfaces that already ship in this workspace, runs bounded
//! in-process probe campaigns against those real implementations, and fits
//! the resulting evidence into canonical [`fs_ladder`] edges through
//! [`fit_fidelity_campaign`]. It is the instance the doctrine machinery
//! (f85xj.10.1–10.3) exists to serve: nodes are model implementations bound
//! to retained model cards, edges are contextual evidence fitted from paired
//! executions, gaps are explicit acquisition demands, and cost never creates
//! authority.
//!
//! Honesty boundaries:
//!
//! - Every node's [`ModelCardRef`] is minted from the retained card's
//!   canonical ledger-row bytes under [`COOLING_CARD_DOMAIN`], so a consumer
//!   holding the card can recompute and verify the binding. That is exact
//!   binding, not provenance or promotion authority.
//! - Probe costs are wall-clock measurements of the in-process evaluations
//!   and are machine-dependent by nature; they flow into the fitted cost
//!   artifacts and therefore into edge/graph identities. The
//!   [`CoolingFidelityInstance::structural_digest`] deliberately excludes
//!   them: it pins the instance SHAPE (nodes, edge contexts, gaps) so a
//!   structural change fails a golden while a re-measured cost does not.
//! - The thermal-LBM natural-convection node, the resolved 3-D PCB rung, and
//!   the lumped-vs-marched transient pairing remain explicit
//!   [`CampaignGap`]s, not silently absent edges.

use std::collections::BTreeMap;
use std::time::Instant;

use fs_blake3::{ContentHash, hash_domain};
use fs_conduction::fixtures::box_grid;
use fs_conduction::material::ConductivityModel;
use fs_conduction::solve::element_heat_flux;
use fs_conduction::{
    ConductionError, ConductionMesh, EMISSIVITY_DIMS, GrayDiffuseEnclosure,
    LinearizedSurfaceRadiation, RadiationSurface, STEFAN_BOLTZMANN_W_M2_K4,
    SURFACE_EMISSIVITY_PROPERTY, SurfaceEmissivity, ViewFactorMatrix,
};
use fs_convection::{
    CorrelationError, CorrelationId, CorrelationInputs, ThermalDirection, correlation_catalog,
    evaluate,
};
use fs_evidence::{Ambition, ModelCard, ValidityDomain as EvidenceValidityDomain};
use fs_exec::Cx;
use fs_ladder::{
    ClosedInterval, CostModelRef, CostRelationRef, DiscrepancyModelRef, DiscrepancyReference,
    EdgeId, FidelityGraph, FidelityNode, GraphError, ModelCardRef, ModelId, QoiId, RegimeAxis,
    TransferRef,
};
use fs_matdb::{
    ClaimSet, CopperCoverage, InterpolationPolicy, MaterialCard, MaterialStateId,
    PCB_THERMAL_CONDUCTIVITY_DIMS, PcbConductivityDatum, PcbHomogenizationError, PcbLayer,
    PcbPrincipalFrame, PcbScaleSeparation, PcbStackup, PropertyClaim, PropertyKey, PropertyValue,
    Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};

use crate::fidelity_campaign::{
    CampaignAuthority, CampaignError, CampaignGap, CampaignRun, EdgeProbeCampaign, FittedCampaign,
    RunPartition, fit_fidelity_campaign,
};

/// Campaign/instance name recorded in artifacts.
pub const COOLING_INSTANCE_NAME: &str = "cooling-fidelity-instance-v1";
/// Domain separator for node model identities.
pub const COOLING_MODEL_DOMAIN: &str = "org.frankensim.fs-plan.cooling-instance.model.v1";
/// Domain separator for card-reference minting over canonical card bytes.
pub const COOLING_CARD_DOMAIN: &str = "org.frankensim.fs-plan.cooling-instance.model-card.v1";
/// Domain separator for declared model-build identities.
pub const COOLING_BUILD_DOMAIN: &str = "org.frankensim.fs-plan.cooling-instance.model-build.v1";
/// Domain separator for the cost-independent structural digest.
pub const COOLING_STRUCTURE_DOMAIN: &str = "org.frankensim.fs-plan.cooling-instance.structure.v1";

/// Query tolerance that separates the near-anchor radiation regime (fitted
/// discrepancy ≈ 5–11%) from the far-cold regime (≈ 15–21%) on the shipped
/// probe sweeps. Demos use it; nothing in the fit depends on it.
pub const RADIATION_REGIME_SPLIT_TOLERANCE: f64 = 0.125;

/// One retained node card and its exact graph binding.
#[derive(Debug, Clone)]
pub struct CoolingNodeCard {
    /// Graph model identity.
    pub model: ModelId,
    /// Human-auditable node label (also the model-identity preimage).
    pub label: String,
    /// The retained model card the node binds to.
    pub card: ModelCard,
    /// `hash_domain(COOLING_CARD_DOMAIN, card.to_ledger_row_json())`.
    pub card_ref: ModelCardRef,
}

/// Per-fitted-edge numeric summary retained for query resolvers and demos.
///
/// The canonical evidence is the fitted artifacts inside the campaign; this
/// summary re-exposes the numbers a [`fs_ladder::EdgeEvidenceResolver`] needs
/// without re-parsing artifact bytes.
#[derive(Debug, Clone)]
pub struct CoolingEdgeSummary {
    /// Fitted edge identity in the assembled graph.
    pub edge: EdgeId,
    /// Cheap endpoint.
    pub source: ModelId,
    /// Finer endpoint.
    pub target: ModelId,
    /// Compared quantity.
    pub qoi: QoiId,
    /// Exact cost-model reference stored on the edge.
    pub cost_ref: CostModelRef,
    /// Exact discrepancy-model reference stored on the edge.
    pub discrepancy_ref: DiscrepancyModelRef,
    /// Mean measured source cost over the campaign runs, seconds.
    pub mean_source_cost_s: f64,
    /// Mean measured target cost over the campaign runs, seconds.
    pub mean_target_cost_s: f64,
    /// Maximum observed relative discrepancy over fit runs.
    pub fitted_max_rel: f64,
}

/// The assembled production instance.
#[derive(Debug, Clone)]
pub struct CoolingFidelityInstance {
    /// Fitted campaign: graph, edges, artifacts, gaps, freshness authority.
    pub campaign: FittedCampaign,
    /// Retained node cards, sorted by model identity.
    pub cards: Vec<CoolingNodeCard>,
    /// Per-edge numeric summaries, in fitted-edge order.
    pub summaries: Vec<CoolingEdgeSummary>,
    /// Cost-independent digest of the instance shape (nodes, edge contexts,
    /// gaps). Golden-pinned by the instance battery.
    pub structural_digest: ContentHash,
}

/// Structured assembly refusal.
#[derive(Debug)]
pub enum CoolingInstanceError {
    /// A real model surface refused its probe fixture.
    Conduction(ConductionError),
    /// A correlation card refused its probe inputs.
    Correlation(CorrelationError),
    /// PCB homogenization refused its stackup.
    Homogenization(PcbHomogenizationError),
    /// A material-card query refused.
    MatDb(fs_matdb::MatDbError),
    /// Graph construction refused.
    Graph(GraphError),
    /// Campaign fitting refused.
    Campaign(CampaignError),
    /// A fitted edge lost the reference shape this module just built.
    Inconsistent(&'static str),
}

impl std::fmt::Display for CoolingInstanceError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Conduction(error) => write!(f, "cooling instance conduction probe: {error}"),
            Self::Correlation(error) => write!(f, "cooling instance correlation probe: {error}"),
            Self::Homogenization(error) => write!(f, "cooling instance PCB stackup: {error}"),
            Self::MatDb(error) => write!(f, "cooling instance material card: {error}"),
            Self::Graph(error) => write!(f, "cooling instance graph: {error}"),
            Self::Campaign(error) => write!(f, "cooling instance campaign: {error}"),
            Self::Inconsistent(what) => write!(f, "cooling instance inconsistency: {what}"),
        }
    }
}

impl std::error::Error for CoolingInstanceError {}

impl From<ConductionError> for CoolingInstanceError {
    fn from(error: ConductionError) -> Self {
        Self::Conduction(error)
    }
}
impl From<CorrelationError> for CoolingInstanceError {
    fn from(error: CorrelationError) -> Self {
        Self::Correlation(error)
    }
}
impl From<PcbHomogenizationError> for CoolingInstanceError {
    fn from(error: PcbHomogenizationError) -> Self {
        Self::Homogenization(error)
    }
}
impl From<fs_matdb::MatDbError> for CoolingInstanceError {
    fn from(error: fs_matdb::MatDbError) -> Self {
        Self::MatDb(error)
    }
}
impl From<GraphError> for CoolingInstanceError {
    fn from(error: GraphError) -> Self {
        Self::Graph(error)
    }
}
impl From<CampaignError> for CoolingInstanceError {
    fn from(error: CampaignError) -> Self {
        Self::Campaign(error)
    }
}

/// Node labels, fixed and load-bearing: they are the model-identity
/// preimages and the vocabulary the demos and product wiring use.
pub mod labels {
    /// Fully developed laminar constant-wall-temperature duct correlation.
    pub const DUCT_LAMINAR_CWT: &str = "convection.duct-laminar-cwt";
    /// Hausen developing laminar duct correlation.
    pub const DUCT_LAMINAR_HAUSEN: &str = "convection.duct-laminar-hausen";
    /// Churchill–Chu natural-convection vertical-plate correlation.
    pub const NATURAL_CHURCHILL_CHU: &str = "convection.natural-churchill-chu";
    /// Linearized surface-to-ambient radiation Robin row.
    pub const RADIATION_LINEARIZED: &str = "radiation.linearized-robin";
    /// Gray-diffuse enclosure radiosity exchange.
    pub const RADIATION_GRAY_DIFFUSE: &str = "radiation.gray-diffuse-enclosure";
    /// Bulk isotropic FR4 board treatment.
    pub const PCB_BULK_ISOTROPIC: &str = "pcb.bulk-isotropic-fr4";
    /// Receipt-backed homogenized stackup board treatment.
    pub const PCB_HOMOGENIZED: &str = "pcb.homogenized-stackup";
    /// Thermal-LBM natural convection (gap target; no probe surface here yet).
    pub const THERMAL_LBM_NATURAL: &str = "thermal-lbm.natural-convection";
    /// Resolved 3-D per-region PCB stackup (gap target; not implemented).
    pub const PCB_RESOLVED_3D: &str = "pcb.resolved-stackup-3d";
    /// Lumped-capacitance transient (gap source; probes not retained here).
    pub const TRANSIENT_LUMPED: &str = "transient.lumped";
    /// Marched transient finite-element reference (gap target).
    pub const TRANSIENT_MARCHED: &str = "transient.full-march";
}

fn model_id(label: &str) -> ModelId {
    ModelId::new(hash_domain(COOLING_MODEL_DOMAIN, label.as_bytes()))
}

/// Mint the exact card binding for a retained card.
pub fn card_ref_for(card: &ModelCard) -> ModelCardRef {
    ModelCardRef::new(hash_domain(
        COOLING_CARD_DOMAIN,
        card.to_ledger_row_json().as_bytes(),
    ))
}

fn catalog_card(id: CorrelationId) -> Result<ModelCard, CoolingInstanceError> {
    correlation_catalog()
        .into_iter()
        .find(|card| card.id == id)
        .map(|card| card.model)
        .ok_or(CoolingInstanceError::Inconsistent(
            "correlation catalog lost a card this instance names",
        ))
}

fn radiation_linearized_card() -> ModelCard {
    ModelCard::new(
        labels::RADIATION_LINEARIZED,
        "1",
        Ambition::Solid,
        vec![
            "gray surface, hemispherical-total emissivity from a material card".to_string(),
            "surface-to-ambient exchange; surroundings act as a blackbody sink".to_string(),
            "flux linearized at the declared mean temperature (h = 4 eps sigma T_m^3)".to_string(),
        ],
        EvidenceValidityDomain::unconstrained().with("surface_temperature_k", 295.0, 345.0),
        vec![
            "relative error grows away from the linearization anchor; the far-cold band of the \
             probe sweep reaches ~20%"
                .to_string(),
        ],
        0.21,
    )
}

fn radiation_gray_diffuse_card() -> ModelCard {
    ModelCard::new(
        labels::RADIATION_GRAY_DIFFUSE,
        "1",
        Ambition::Solid,
        vec![
            "gray-diffuse opaque surfaces over an admitted view-factor matrix".to_string(),
            "radiosity balance solved exactly for the enclosure".to_string(),
            "no participating media, spectral, or specular behavior".to_string(),
        ],
        EvidenceValidityDomain::unconstrained().with("surface_temperature_k", 250.0, 1000.0),
        vec![
            "does not generate view factors; the admitted matrix's evidence bounds the claim"
                .to_string(),
        ],
        0.02,
    )
}

fn pcb_bulk_card() -> ModelCard {
    ModelCard::new(
        labels::PCB_BULK_ISOTROPIC,
        "1",
        Ambition::Solid,
        vec!["board treated as bulk isotropic FR4 (k = 0.25 W/m-K)".to_string()],
        EvidenceValidityDomain::unconstrained().with("T", 250.0, 400.0),
        vec![
            "ignores copper planes entirely; in-plane conduction is underestimated by more than \
             an order of magnitude on plane-heavy stackups"
                .to_string(),
        ],
        1.0,
    )
}

fn pcb_homogenized_card() -> ModelCard {
    ModelCard::new(
        labels::PCB_HOMOGENIZED,
        "1",
        Ambition::Solid,
        vec![
            "series/parallel laminate homogenization from receipt-backed material cards"
                .to_string(),
            "scale separation declared and admitted before homogenizing".to_string(),
        ],
        EvidenceValidityDomain::unconstrained().with("T", 250.0, 400.0),
        vec![
            "principal-frame effective tensor; resolved per-feature detail is out of scope"
                .to_string(),
        ],
        0.15,
    )
}

fn emissivity_material_card(name: &str, emissivity: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS),
            value: PropertyValue::Scalar {
                value: emissivity,
                dims: EMISSIVITY_DIMS,
            },
            validity: EvidenceValidityDomain::unconstrained().with("T", 250.0, 500.0),
            uncertainty: UncertaintyModel::HalfWidth {
                half_width: 0.01,
                confidence: 0.95,
            },
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: Vec::new(),
            provenance: Provenance {
                source: format!("cooling-instance declared surface finish: {name}"),
                license: "internal-declared".to_string(),
                artifact: None,
            },
        })
        .expect("one well-formed emissivity claim");
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: name.to_string(),
            phase: "solid".to_string(),
            process: "declared-finish".to_string(),
            revision: 0,
        },
        claims,
        Vec::new(),
    )
    .expect("one-claim emissivity card assembles")
}

fn pcb_material_card(chemistry: &str, conductivity: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new("thermal_conductivity", PCB_THERMAL_CONDUCTIVITY_DIMS),
            value: PropertyValue::Scalar {
                value: conductivity,
                dims: PCB_THERMAL_CONDUCTIVITY_DIMS,
            },
            validity: EvidenceValidityDomain::unconstrained().with("T", 250.0, 400.0),
            uncertainty: UncertaintyModel::HalfWidth {
                half_width: 0.0,
                confidence: 0.95,
            },
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: Vec::new(),
            provenance: Provenance {
                source: format!("cooling-instance declared laminate material: {chemistry}"),
                license: "internal-declared".to_string(),
                artifact: None,
            },
        })
        .expect("one well-formed conductivity claim");
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: chemistry.to_string(),
            phase: "solid".to_string(),
            process: "laminate".to_string(),
            revision: 0,
        },
        claims,
        Vec::new(),
    )
    .expect("one-claim conductivity card assembles")
}

fn measured<T>(f: impl FnOnce() -> T) -> (T, f64) {
    let start = Instant::now();
    let value = f();
    (value, start.elapsed().as_secs_f64().max(1.0e-9))
}

fn probe_run(
    label: &str,
    partition: RunPartition,
    params: BTreeMap<String, f64>,
    problem_size: f64,
    qois: (f64, f64, Option<f64>),
    costs: (f64, f64),
) -> CampaignRun {
    CampaignRun {
        run_id: hash_domain(COOLING_MODEL_DOMAIN, format!("run:{label}").as_bytes()),
        case_id: label.to_string(),
        partition,
        params,
        problem_size,
        source_qoi: qois.0,
        target_qoi: qois.1,
        reference_qoi: qois.2,
        source_cost_s: costs.0,
        target_cost_s: costs.1,
    }
}

fn regime(
    entries: &[(&str, f64, f64)],
) -> Result<BTreeMap<RegimeAxis, ClosedInterval>, GraphError> {
    entries
        .iter()
        .map(|(axis, lower, upper)| {
            Ok((
                RegimeAxis::new(*axis)?,
                ClosedInterval::new(*lower, *upper)?,
            ))
        })
        .collect()
}

fn params(entries: &[(&str, f64)]) -> BTreeMap<String, f64> {
    entries
        .iter()
        .map(|(name, value)| ((*name).to_string(), *value))
        .collect()
}

const HOT_EMISSIVITY: f64 = 0.8;
const SINK_EMISSIVITY: f64 = 0.99;
const AMBIENT_K: f64 = 300.0;
const ANCHOR_MEAN_K: f64 = 320.0;

/// Closed-form two-parallel-plate gray exchange: the independent reference
/// the radiation probes are audited against.
fn parallel_plate_exchange_w_m2(hot_k: f64, sink_k: f64) -> f64 {
    let factor = 1.0 / (1.0 / HOT_EMISSIVITY + 1.0 / SINK_EMISSIVITY - 1.0);
    factor * STEFAN_BOLTZMANN_W_M2_K4 * (hot_k.powi(4) - sink_k.powi(4))
}

struct RadiationProbe {
    linearized: LinearizedSurfaceRadiation,
    enclosure: GrayDiffuseEnclosure,
}

impl RadiationProbe {
    fn build() -> Result<Self, CoolingInstanceError> {
        let hot = SurfaceEmissivity::from_card(
            "hot-plate",
            &emissivity_material_card("hot-plate-finish", HOT_EMISSIVITY),
            ANCHOR_MEAN_K,
            SelectionPolicy::SingleClaimOnly,
        )?;
        let linearized =
            LinearizedSurfaceRadiation::new("hot-plate", hot, ANCHOR_MEAN_K, AMBIENT_K, 25.0)?;

        let (complex, positions) = box_grid([2, 2, 2], [1.0, 1.0, 1.0]);
        let mesh = ConductionMesh::new(complex, positions)?;
        let hot_surface = RadiationSurface::new(
            &mesh,
            "hot-plate",
            |face| face.centroid[0].abs() < 1.0e-9,
            SurfaceEmissivity::from_card(
                "hot-plate",
                &emissivity_material_card("hot-plate-finish", HOT_EMISSIVITY),
                ANCHOR_MEAN_K,
                SelectionPolicy::SingleClaimOnly,
            )?,
        )?;
        let area = hot_surface.area_m2();
        let sink_surface = RadiationSurface::new(
            &mesh,
            "surroundings",
            |face| (face.centroid[0] - 1.0).abs() < 1.0e-9,
            SurfaceEmissivity::from_card(
                "surroundings",
                &emissivity_material_card("near-black-surroundings", SINK_EMISSIVITY),
                AMBIENT_K,
                SelectionPolicy::SingleClaimOnly,
            )?,
        )?;
        let matrix = ViewFactorMatrix::infinite_parallel_plates(area)?;
        let enclosure = GrayDiffuseEnclosure::new(vec![hot_surface, sink_surface], matrix)?;
        Ok(Self {
            linearized,
            enclosure,
        })
    }

    fn campaign(
        &self,
        cx: &Cx<'_>,
        source: ModelId,
        target: ModelId,
        bin_name: &str,
        bin_k: (f64, f64),
        fit_temperatures: [f64; 4],
        held_out_temperatures: [f64; 2],
    ) -> Result<EdgeProbeCampaign, CoolingInstanceError> {
        let mut runs = Vec::new();
        let all = fit_temperatures
            .iter()
            .map(|t| (*t, RunPartition::Fit))
            .chain(
                held_out_temperatures
                    .iter()
                    .map(|t| (*t, RunPartition::HeldOut)),
            );
        for (temperature, partition) in all {
            let (point, source_cost) = measured(|| self.linearized.evaluate(temperature));
            let point = point?;
            let (report, target_cost) =
                measured(|| self.enclosure.solve(cx, &[temperature, AMBIENT_K]));
            let report = report?;
            runs.push(probe_run(
                &format!("radiation-{bin_name}-{temperature:.0}K"),
                partition,
                params(&[
                    ("surface_temperature_k", temperature),
                    ("ambient_temperature_k", AMBIENT_K),
                    ("emissivity", HOT_EMISSIVITY),
                ]),
                1.0,
                (
                    point.linearized_outward_flux_w_m2,
                    report.net_outward_flux_w_m2[0],
                    Some(parallel_plate_exchange_w_m2(temperature, AMBIENT_K)),
                ),
                (source_cost, target_cost),
            ));
        }
        Ok(EdgeProbeCampaign {
            source,
            target,
            qoi: QoiId::new("outward-radiative-flux").map_err(CoolingInstanceError::Graph)?,
            qoi_unit: "W/m2".to_string(),
            transfer: TransferRef::new(hash_domain(
                COOLING_MODEL_DOMAIN,
                format!("transfer:radiation:{bin_name}").as_bytes(),
            )),
            regime_bin: regime(&[
                ("surface_temperature_k", bin_k.0, bin_k.1),
                ("ambient_temperature_k", AMBIENT_K, AMBIENT_K),
                ("emissivity", HOT_EMISSIVITY, HOT_EMISSIVITY),
            ])?,
            runs,
        })
    }
}

fn homogenized_stackup() -> Result<fs_matdb::PcbHomogenizedConductivity, CoolingInstanceError> {
    let copper = pcb_material_card("C11000-copper", 400.0);
    let fr4 = pcb_material_card("FR4-dielectric", 0.25);
    let point = QueryPoint::new().with("T", 300.0)?;
    let copper_datum = PcbConductivityDatum::from_card(
        &copper,
        "thermal_conductivity",
        &point,
        SelectionPolicy::SingleClaimOnly,
    )?;
    let fr4_datum = PcbConductivityDatum::from_card(
        &fr4,
        "thermal_conductivity",
        &point,
        SelectionPolicy::SingleClaimOnly,
    )?;
    let provenance = |what: &str| Provenance {
        source: format!("cooling-instance declared stackup: {what}"),
        license: "internal-declared".to_string(),
        artifact: None,
    };
    let plane = PcbLayer::new(
        "L1-plane",
        0.2e-3,
        copper_datum.clone(),
        fr4_datum.clone(),
        CopperCoverage::new("coverage/L1", 0.95, 0.90, 1.0, provenance("L1 coverage"))?,
    )?;
    let core = PcbLayer::new(
        "core",
        0.8e-3,
        copper_datum,
        fr4_datum,
        CopperCoverage::new(
            "coverage/core",
            0.05,
            0.0,
            0.10,
            provenance("core coverage"),
        )?,
    )?;
    Ok(PcbStackup::new(
        "cooling-instance-pcb",
        vec![plane, core],
        PcbPrincipalFrame::default(),
        PcbScaleSeparation::new(25.0e-6, 0.05)?,
    )?
    .homogenize()?)
}

/// Mean |x-flux| over the board mesh under an imposed linear x-gradient.
///
/// Deliberately not a solve: the probe isolates the constitutive treatment
/// (the axis on which the two PCB rungs differ) under an identical unit
/// thermal loading, so the paired QoI difference is exactly the model-form
/// difference.
fn board_mean_x_flux(
    mesh: &ConductionMesh,
    positions: &[[f64; 3]],
    material: &ConductivityModel,
    gradient_k_per_m: f64,
) -> Result<f64, CoolingInstanceError> {
    let temperature: Vec<f64> = positions
        .iter()
        .map(|p| AMBIENT_K + gradient_k_per_m * p[0])
        .collect();
    let flux = element_heat_flux(mesh, material, &temperature)?;
    if flux.is_empty() {
        return Err(CoolingInstanceError::Inconsistent(
            "board fixture produced no elements",
        ));
    }
    Ok(flux.iter().map(|value| value[0].abs()).sum::<f64>() / flux.len() as f64)
}

fn pcb_campaign(
    source: ModelId,
    target: ModelId,
) -> Result<EdgeProbeCampaign, CoolingInstanceError> {
    let (complex, positions) = fs_conduction::fixtures::unit_cube(2);
    let mesh = ConductionMesh::new(complex, positions.clone())?;
    let bulk = ConductivityModel::isotropic_declared(0.25)?;
    let homogenized = ConductivityModel::from_pcb_homogenization(&homogenized_stackup()?)?;
    let gradients = [5.0, 10.0, 15.0, 20.0, 25.0, 30.0];
    let mut runs = Vec::new();
    for (index, gradient) in gradients.into_iter().enumerate() {
        let (source_qoi, source_cost) =
            measured(|| board_mean_x_flux(&mesh, &positions, &bulk, gradient));
        let (target_qoi, target_cost) =
            measured(|| board_mean_x_flux(&mesh, &positions, &homogenized, gradient));
        runs.push(probe_run(
            &format!("pcb-gradient-{gradient:.0}K"),
            if index < 4 {
                RunPartition::Fit
            } else {
                RunPartition::HeldOut
            },
            params(&[("x_gradient_k_per_m", gradient), ("T", 300.0)]),
            1.0,
            (source_qoi?, target_qoi?, None),
            (source_cost, target_cost),
        ));
    }
    Ok(EdgeProbeCampaign {
        source,
        target,
        qoi: QoiId::new("board-mean-x-flux").map_err(CoolingInstanceError::Graph)?,
        qoi_unit: "W/m2".to_string(),
        transfer: TransferRef::new(hash_domain(
            COOLING_MODEL_DOMAIN,
            b"transfer:pcb:bulk-to-homogenized",
        )),
        regime_bin: regime(&[("x_gradient_k_per_m", 5.0, 30.0), ("T", 300.0, 300.0)])?,
        runs,
    })
}

fn convection_campaign(
    source: ModelId,
    target: ModelId,
) -> Result<EdgeProbeCampaign, CoolingInstanceError> {
    let cases: [(f64, f64); 6] = [
        (500.0, 0.7),
        (800.0, 1.0),
        (1_100.0, 2.0),
        (1_400.0, 4.0),
        (1_700.0, 7.0),
        (2_000.0, 10.0),
    ];
    let length_ratio = 100.0;
    let mut runs = Vec::new();
    for (index, (reynolds, prandtl)) in cases.into_iter().enumerate() {
        let inputs = CorrelationInputs::forced(reynolds, prandtl)
            .with_direction(ThermalDirection::CoolingFluid)
            .with_length_ratio(length_ratio);
        let (source_eval, source_cost) =
            measured(|| evaluate(CorrelationId::CircularDuctLaminarCwt, inputs));
        let source_eval = source_eval?;
        let (target_eval, target_cost) =
            measured(|| evaluate(CorrelationId::CircularDuctHausen, inputs));
        let target_eval = target_eval?;
        runs.push(probe_run(
            &format!("duct-cwt-vs-hausen-{index}"),
            if index < 4 {
                RunPartition::Fit
            } else {
                RunPartition::HeldOut
            },
            params(&[
                ("Re", reynolds),
                ("Pr", prandtl),
                ("Pe", reynolds * prandtl),
                ("L_over_Dh", length_ratio),
            ]),
            reynolds,
            (
                source_eval.evidence().value,
                target_eval.evidence().value,
                None,
            ),
            (source_cost, target_cost),
        ));
    }
    Ok(EdgeProbeCampaign {
        source,
        target,
        qoi: QoiId::new("nusselt-number").map_err(CoolingInstanceError::Graph)?,
        qoi_unit: "1".to_string(),
        transfer: TransferRef::new(hash_domain(
            COOLING_MODEL_DOMAIN,
            b"transfer:convection:cwt-to-hausen",
        )),
        regime_bin: regime(&[
            ("Re", 500.0, 2_000.0),
            ("Pr", 0.6, 10.0),
            ("Pe", 300.0, 20_000.0),
            ("L_over_Dh", length_ratio, length_ratio),
        ])?,
        runs,
    })
}

fn structural_digest(
    cards: &[CoolingNodeCard],
    campaigns: &[EdgeProbeCampaign],
    gaps: &[CampaignGap],
) -> ContentHash {
    let mut text = String::new();
    text.push_str(COOLING_INSTANCE_NAME);
    text.push('\n');
    for card in cards {
        text.push_str(&format!(
            "node|{}|{}|{}\n",
            card.label, card.model, card.card_ref
        ));
    }
    for campaign in campaigns {
        text.push_str(&format!(
            "edge|{}|{}|{}|{}",
            campaign.source, campaign.target, campaign.qoi, campaign.qoi_unit
        ));
        for (axis, interval) in &campaign.regime_bin {
            text.push_str(&format!(
                "|{axis}:{}:{}",
                interval.lower(),
                interval.upper()
            ));
        }
        text.push('\n');
    }
    for gap in gaps {
        text.push_str(&format!(
            "gap|{}|{}|{}|{}\n",
            gap.source, gap.target, gap.qoi, gap.reason
        ));
    }
    hash_domain(COOLING_STRUCTURE_DOMAIN, text.as_bytes())
}

/// Assemble the production cooling fidelity-graph instance.
///
/// Runs every probe in-process at fixture scale under the caller's `Cx`,
/// fits the campaign, and returns the instance with retained cards, per-edge
/// summaries, and the cost-independent structural digest.
///
/// # Errors
///
/// Refuses when any real model surface refuses its probe fixture, when graph
/// construction refuses, or when the campaign fitter refuses.
pub fn assemble_cooling_instance(
    cx: &Cx<'_>,
) -> Result<CoolingFidelityInstance, CoolingInstanceError> {
    let node = |label: &str, card: ModelCard| CoolingNodeCard {
        model: model_id(label),
        label: label.to_string(),
        card_ref: card_ref_for(&card),
        card,
    };
    let mut cards = vec![
        node(
            labels::DUCT_LAMINAR_CWT,
            catalog_card(CorrelationId::CircularDuctLaminarCwt)?,
        ),
        node(
            labels::DUCT_LAMINAR_HAUSEN,
            catalog_card(CorrelationId::CircularDuctHausen)?,
        ),
        node(
            labels::NATURAL_CHURCHILL_CHU,
            catalog_card(CorrelationId::ChurchillChuVerticalPlate)?,
        ),
        node(labels::RADIATION_LINEARIZED, radiation_linearized_card()),
        node(
            labels::RADIATION_GRAY_DIFFUSE,
            radiation_gray_diffuse_card(),
        ),
        node(labels::PCB_BULK_ISOTROPIC, pcb_bulk_card()),
        node(labels::PCB_HOMOGENIZED, pcb_homogenized_card()),
    ];
    cards.sort_by_key(|card| card.model);

    let mut graph =
        FidelityGraph::new(COOLING_INSTANCE_NAME).map_err(CoolingInstanceError::Graph)?;
    for card in &cards {
        graph
            .add_node(
                FidelityNode::new(card.model, card.card_ref, card.label.clone())
                    .map_err(CoolingInstanceError::Graph)?,
            )
            .map_err(CoolingInstanceError::Graph)?;
    }

    let radiation = RadiationProbe::build()?;
    let campaigns = vec![
        radiation.campaign(
            cx,
            model_id(labels::RADIATION_LINEARIZED),
            model_id(labels::RADIATION_GRAY_DIFFUSE),
            "near-anchor",
            (316.0, 332.0),
            [321.0, 324.0, 327.0, 330.0],
            [323.0, 329.0],
        )?,
        radiation.campaign(
            cx,
            model_id(labels::RADIATION_LINEARIZED),
            model_id(labels::RADIATION_GRAY_DIFFUSE),
            "far-cold",
            (302.0, 314.0),
            [303.0, 306.0, 309.0, 312.0],
            [304.0, 310.0],
        )?,
        pcb_campaign(
            model_id(labels::PCB_BULK_ISOTROPIC),
            model_id(labels::PCB_HOMOGENIZED),
        )?,
        convection_campaign(
            model_id(labels::DUCT_LAMINAR_CWT),
            model_id(labels::DUCT_LAMINAR_HAUSEN),
        )?,
    ];
    let gaps = vec![
        CampaignGap {
            source: model_id(labels::NATURAL_CHURCHILL_CHU),
            target: model_id(labels::THERMAL_LBM_NATURAL),
            qoi: QoiId::new("nusselt-number").map_err(CoolingInstanceError::Graph)?,
            reason: "the retained correlation is a vertical external plate while the workspace \
                     thermal-LBM fixture is a horizontal periodic Rayleigh-Benard slab; no \
                     shared-validity geometry has a retained paired execution set yet"
                .to_string(),
        },
        CampaignGap {
            source: model_id(labels::PCB_HOMOGENIZED),
            target: model_id(labels::PCB_RESOLVED_3D),
            qoi: QoiId::new("board-mean-x-flux").map_err(CoolingInstanceError::Graph)?,
            reason: "fs-conduction has no per-region material map, so a resolved 3-D stackup \
                     rung cannot execute; acquiring it requires the per-element material \
                     frontend before any probe campaign"
                .to_string(),
        },
        CampaignGap {
            source: model_id(labels::TRANSIENT_LUMPED),
            target: model_id(labels::TRANSIENT_MARCHED),
            qoi: QoiId::new("node-temperature").map_err(CoolingInstanceError::Graph)?,
            reason: "the lumped-vs-marched pairing exists as an fs-conduction test but no \
                     paired probe set is retained through this instance yet"
                .to_string(),
        },
    ];

    let digest = structural_digest(&cards, &campaigns, &gaps);

    let mut model_builds = BTreeMap::new();
    for card in &cards {
        model_builds.insert(
            card.model,
            hash_domain(
                COOLING_BUILD_DOMAIN,
                format!(
                    "{}|fs-plan {}|declared-in-process-probe",
                    card.label,
                    env!("CARGO_PKG_VERSION")
                )
                .as_bytes(),
            ),
        );
    }
    // Gap endpoints that are not graph nodes still name model identities;
    // they carry no build entries because nothing executed them.
    let authority = CampaignAuthority {
        corpus: digest,
        corpus_version: 1,
        machine_fingerprint: b"in-process-probe/wall-clock-costs-are-machine-relative".to_vec(),
        model_builds,
    };

    // The fitter re-sorts edges by identity, but it preserves each input
    // campaign's distinct TransferRef verbatim on the fitted edge, so the
    // transfer reference is the exact join key back to probe metadata.
    let mut by_transfer: BTreeMap<TransferRef, (QoiId, f64, f64)> = BTreeMap::new();
    for campaign in &campaigns {
        let n = campaign.runs.len() as f64;
        by_transfer.insert(
            campaign.transfer,
            (
                campaign.qoi.clone(),
                campaign.runs.iter().map(|r| r.source_cost_s).sum::<f64>() / n,
                campaign.runs.iter().map(|r| r.target_cost_s).sum::<f64>() / n,
            ),
        );
    }

    let campaign = fit_fidelity_campaign(COOLING_INSTANCE_NAME, graph, authority, campaigns, gaps)
        .map_err(CoolingInstanceError::Campaign)?;

    let mut summaries = Vec::new();
    for fitted in &campaign.edges {
        let CostRelationRef::Model(cost_ref) = fitted.edge.cost() else {
            return Err(CoolingInstanceError::Inconsistent(
                "fitted edge lost its native cost-model reference",
            ));
        };
        let DiscrepancyReference::Model(discrepancy_ref) = fitted.edge.discrepancy() else {
            return Err(CoolingInstanceError::Inconsistent(
                "fitted edge lost its native discrepancy-model reference",
            ));
        };
        let Some((qoi, mean_source, mean_target)) = by_transfer.get(&fitted.edge.transfer()) else {
            return Err(CoolingInstanceError::Inconsistent(
                "fitted edge carries a transfer reference this instance never issued",
            ));
        };
        summaries.push(CoolingEdgeSummary {
            edge: fitted.edge.id(),
            source: fitted.edge.source(),
            target: fitted.edge.target(),
            qoi: qoi.clone(),
            cost_ref,
            discrepancy_ref,
            mean_source_cost_s: *mean_source,
            mean_target_cost_s: *mean_target,
            fitted_max_rel: fitted.discrepancy_band.max_observed_rel,
        });
    }

    Ok(CoolingFidelityInstance {
        campaign,
        cards,
        summaries,
        structural_digest: digest,
    })
}
