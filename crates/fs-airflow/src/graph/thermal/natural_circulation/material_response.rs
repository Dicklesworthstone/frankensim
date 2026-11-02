//! Material-updated drag and forced convection on an unchanged native graph.
//!
//! Original material owners remain attached to every numerical response.
//! Reference density owns mass/volume conversion and inertial drag; the actual
//! source viscosity updates linear drag and actual transport supplies Re, Pr,
//! conductivity and h through the existing correlation owner. This is a
//! representative-state closure, not axial integration of viscosity or a
//! second property database. The nonlinear caller selects each material's
//! bulk/film state and retains authority over phase and source applicability.

use fs_convection::{CorrelationId, HeatTransferCoefficient, NusseltEvaluation, ThermalDirection};
use fs_evidence::Evidence;
use fs_exec::Cx;
use fs_material::MaterialError;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_qty::MassFlowRate;

use crate::conjugate::AirSegment;
use crate::graph::thermal::medium::FlowHeatGeometry;
use crate::graph::thermal::medium::convection::{
    ForcedConvectionSettings, MaterialConvectionError, evaluate_forced_from_material,
};
use crate::graph::thermal::medium::series::hydraulic::pressure::BoussinesqDensity;
use crate::graph::thermal::medium::series::hydraulic::{
    MaterialFlowResistance, ResistanceError, ResistanceEvaluation,
};
use crate::graph::thermal::transport::BranchThermalModel;
use crate::graph::{
    GraphBranch, GraphError, LinearLossElement, LinearLossResistance, LossGraph, PassiveBranch,
};
use crate::{AirflowError, LossElement, LossResistance, SourceProvenance, ToleranceBasis};

/// One exchange segment in the branch's original declared order. Geometry
/// and correlation selection are physical inputs, not fitted output targets.
#[derive(Debug, Clone)]
pub struct MaterialExchangeLaw {
    pub region: String,
    pub geometry: FlowHeatGeometry,
    pub correlation: CorrelationId,
    /// The caller updates this convention when wall/bulk heating reverses.
    pub direction: ThermalDirection,
}

/// Original geometric resistance calibration and its heat-exchange segments.
/// The source describes the geometric law, not a detached coefficient at some
/// previous fluid state. Uncertainty is retained, not propagated by this adapter.
#[derive(Debug, Clone)]
pub struct MaterialBranchLaw {
    pub resistance: MaterialFlowResistance,
    pub source: SourceProvenance,
    pub uncertainty_rel: f64,
    pub tolerance_basis: ToleranceBasis,
    /// Empty means an explicitly adiabatic branch.
    pub exchange: Vec<MaterialExchangeLaw>,
}

/// The actual source material and the native material-resistance response.
/// Detaching the graph's scalar coefficients does not transfer this authority.
#[derive(Debug)]
pub struct MaterialBranchEvaluation<M> {
    pub medium: M,
    pub evaluation: ResistanceEvaluation,
    pub linear_volume_resistance_pa_s_m3: f64,
    pub quadratic_volume_resistance_pa_s2_m6: f64,
}

/// Original property-state owner and correlation evidence for one h value.
#[derive(Debug)]
pub struct MaterialSegmentEvaluation<M> {
    pub medium: M,
    pub region: String,
    pub geometry: FlowHeatGeometry,
    pub transport: FluidTransportState,
    pub nusselt: NusseltEvaluation,
    pub coefficient: Evidence<HeatTransferCoefficient>,
}

/// A new coefficient state on exactly the template graph's topology and names.
/// Gravity is deliberately unbound: the coupled owner binds its integrated
/// density profile, never these representative property-state densities.
#[derive(Debug)]
pub struct MaterialGraphResponse<M> {
    pub graph: LossGraph,
    pub models: Vec<BranchThermalModel>,
    pub momentum: Vec<MaterialBranchEvaluation<M>>,
    /// Branch-major, then original declared segment order, even at reverse flow.
    pub heat: Vec<Vec<MaterialSegmentEvaluation<M>>>,
}

#[derive(Debug)]
pub enum MaterialResponseError {
    Invalid(&'static str),
    Branch {
        branch: usize,
        reason: &'static str,
    },
    Material {
        branch: usize,
        segment: Option<usize>,
        source: MaterialError,
    },
    Resistance {
        branch: usize,
        source: ResistanceError,
    },
    Convection {
        branch: usize,
        segment: usize,
        source: MaterialConvectionError,
    },
    Graph(GraphError),
    Airflow(AirflowError),
    Interrupted,
}
impl core::fmt::Display for MaterialResponseError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "material graph response: {self:?}")
    }
}
impl std::error::Error for MaterialResponseError {}
impl From<GraphError> for MaterialResponseError {
    fn from(error: GraphError) -> Self {
        Self::Graph(error)
    }
}
impl From<AirflowError> for MaterialResponseError {
    fn from(error: AirflowError) -> Self {
        Self::Airflow(error)
    }
}

/// Refresh one native coefficient state from explicit current flows/materials.
///
/// Samples each original material exactly once. Branch bulk states determine
/// viscosity; each independently supplied segment state supplies the existing
/// forced-correlation inputs and conductivity. The same snapshot is retained
/// with its original owner, and all original Reynolds/correlation limits apply.
/// Original segment order is retained for negative flows; the transport owner
/// alone reverses traversal. No source, regime, velocity or h fallback exists.
///
/// This does not solve for flow, infer a film convention, update heat capacity
/// in the conservation equation, add dissipative heat or establish convergence.
/// The coupled caller must rebuild and check the final coefficient response.
#[allow(clippy::too_many_arguments)]
pub fn evaluate_material_graph<M: FluidTransport>(
    cx: &Cx<'_>,
    template: &LossGraph,
    density: BoussinesqDensity,
    laws: &[MaterialBranchLaw],
    signed_volume_flows_m3_s: &[f64],
    bulk_materials: Vec<M>,
    segment_materials: Vec<Vec<M>>,
) -> Result<MaterialGraphResponse<M>, MaterialResponseError> {
    poll(cx)?;
    let count = template.branches().len();
    if count == 0
        || laws.len() != count
        || signed_volume_flows_m3_s.len() != count
        || bulk_materials.len() != count
        || segment_materials.len() != count
    {
        return Err(MaterialResponseError::Invalid(
            "one law, flow and material set per original branch required",
        ));
    }
    let mut graph_branches = Vec::with_capacity(count);
    let mut models = Vec::with_capacity(count);
    let mut momentum = Vec::with_capacity(count);
    let mut heat = Vec::with_capacity(count);
    for (branch, (medium, segments)) in bulk_materials
        .into_iter()
        .zip(segment_materials)
        .enumerate()
    {
        poll(cx)?;
        let law = &laws[branch];
        let edge = &template.branches()[branch];
        let q = signed_volume_flows_m3_s[branch];
        if !q.is_finite()
            || segments.len() != law.exchange.len()
            || !law.uncertainty_rel.is_finite()
            || !(0.0..1.0).contains(&law.uncertainty_rel)
            || law.source.citation.trim().is_empty()
            || law.source.identifier.trim().is_empty()
        {
            return Err(MaterialResponseError::Branch {
                branch,
                reason: "invalid flow, material arity or geometric-law authority",
            });
        }
        let state = medium
            .transport_state()
            .map_err(|source| MaterialResponseError::Material {
                branch,
                segment: None,
                source,
            })?;
        admit_density(branch, state, density)?;
        let mass = density.reference_kg_m3() * q;
        if !mass.is_finite() || (q != 0.0 && mass == 0.0) {
            return Err(MaterialResponseError::Branch {
                branch,
                reason: "reference mass flow is not representable",
            });
        }
        let evaluation = law
            .resistance
            .evaluate_snapshot(state, MassFlowRate::new(mass), density.reference_kg_m3())
            .map_err(|source| MaterialResponseError::Resistance { branch, source })?;
        let (linear, quadratic) = law
            .resistance
            .volume_coefficients(state, density.reference_kg_m3())
            .map_err(|source| MaterialResponseError::Resistance { branch, source })?;
        let linear_loss = if linear > 0.0 {
            Some(LinearLossElement::new(
                edge.loss.name.clone(),
                LinearLossResistance::new(linear),
                law.uncertainty_rel,
                law.source.clone(),
                law.tolerance_basis,
            )?)
        } else {
            None
        };
        let replacement = if quadratic > 0.0 {
            let loss = LossElement::new(
                edge.loss.name.clone(),
                LossResistance::new(quadratic),
                law.uncertainty_rel,
                law.source.clone(),
                law.tolerance_basis,
            )?;
            let branch = GraphBranch {
                from: edge.from,
                to: edge.to,
                loss,
            };
            match linear_loss {
                Some(linear_loss) => PassiveBranch::LinearQuadratic {
                    branch,
                    linear_loss,
                },
                None => PassiveBranch::Quadratic(branch),
            }
        } else {
            PassiveBranch::Linear {
                from: edge.from,
                to: edge.to,
                loss: linear_loss.ok_or(MaterialResponseError::Branch {
                    branch,
                    reason: "material resistance must retain an active mechanism",
                })?,
            }
        };
        graph_branches.push(replacement);
        momentum.push(MaterialBranchEvaluation {
            medium,
            evaluation,
            linear_volume_resistance_pa_s_m3: linear,
            quadratic_volume_resistance_pa_s2_m6: quadratic,
        });
        let mut thermal_segments = Vec::with_capacity(segments.len());
        let mut responses = Vec::with_capacity(segments.len());
        for (segment, medium) in segments.into_iter().enumerate() {
            poll(cx)?;
            let exchange = &law.exchange[segment];
            let geometry = exchange.geometry;
            if exchange.region.trim().is_empty()
                || !positive(geometry.flow_area.value())
                || !positive(geometry.characteristic_length.value())
                || !positive(geometry.heated_length.value())
                || !positive(geometry.wetted_area.value())
                || matches!(
                    exchange.correlation,
                    CorrelationId::CircularDuctLaminarChf
                        | CorrelationId::RectangularDuctLaminarChf
                )
            {
                return Err(MaterialResponseError::Branch {
                    branch,
                    reason: "constant-wall exchange needs explicit positive geometry and a compatible correlation",
                });
            }
            let speed = q.abs() / geometry.flow_area.value();
            let length_ratio =
                geometry.heated_length.value() / geometry.characteristic_length.value();
            if !positive(speed) || !positive(length_ratio) {
                return Err(MaterialResponseError::Branch {
                    branch,
                    reason: "forced-convection speed and length ratio must be positive and representable",
                });
            }
            let state =
                medium
                    .transport_state()
                    .map_err(|source| MaterialResponseError::Material {
                        branch,
                        segment: Some(segment),
                        source,
                    })?;
            admit_density(branch, state, density)?;
            let convection = evaluate_forced_from_material(
                exchange.correlation,
                &state,
                speed,
                geometry.characteristic_length,
                ForcedConvectionSettings {
                    length_over_hydraulic_diameter: Some(length_ratio),
                    aspect_ratio: geometry.aspect_ratio,
                    direction: exchange.direction,
                },
            )
            .map_err(|source| MaterialResponseError::Convection {
                branch,
                segment,
                source,
            })?;
            let coefficient = convection.coefficient().clone();
            let nusselt = convection.nusselt().clone();
            thermal_segments.push(AirSegment::new(
                &exchange.region,
                geometry.wetted_area.value(),
                coefficient.value.value(),
            )?);
            responses.push(MaterialSegmentEvaluation {
                medium,
                region: exchange.region.clone(),
                geometry,
                transport: state,
                nusselt,
                coefficient,
            });
        }
        models.push(if thermal_segments.is_empty() {
            BranchThermalModel::Adiabatic
        } else {
            BranchThermalModel::Exchange(thermal_segments)
        });
        heat.push(responses);
    }
    poll(cx)?;
    Ok(MaterialGraphResponse {
        graph: LossGraph::new_mixed(template.node_count(), graph_branches)?,
        models,
        momentum,
        heat,
    })
}

fn admit_density(
    branch: usize,
    state: FluidTransportState,
    density: BoussinesqDensity,
) -> Result<(), MaterialResponseError> {
    let departure = (state.density_kg_m3() - density.reference_kg_m3()) / density.reference_kg_m3();
    if !departure.is_finite() || departure.abs() > density.maximum_relative_departure() {
        Err(MaterialResponseError::Branch {
            branch,
            reason: "material density exceeds the declared Boussinesq allowance",
        })
    } else {
        Ok(())
    }
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn poll(cx: &Cx<'_>) -> Result<(), MaterialResponseError> {
    cx.checkpoint()
        .map_err(|_| MaterialResponseError::Interrupted)
}

#[cfg(test)]
mod tests;
