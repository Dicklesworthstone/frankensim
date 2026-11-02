//! Spatial material drag and heat exchange on the original graph topology.
//!
//! Each authored hydraulic section retains its own bulk material and heat-film
//! states. Series pressure drops add at the SAME reference volume flow:
//! A = sum(a_i * mu_i), B = sum(b_i * rho0). Heat segments remain in declared
//! spatial order; the native transport owner alone reverses traversal.
//!
//! Section factors already include their physical length/calibration. Splitting
//! a channel must partition those factors and wetted areas; it does not reset a
//! developing-flow entrance length. This is a piecewise frozen-property closure,
//! not an accuracy certificate for the caller's spatial samples.

use core::ops::Range;

use fs_exec::Cx;
use fs_material::fluid_transport::FluidTransport;

use super::material_response::{
    MaterialBranchEvaluation, MaterialBranchLaw, MaterialGraphResponse, MaterialResponseError,
    MaterialSegmentEvaluation, evaluate_material_graph,
};
use crate::graph::thermal::medium::series::MAX_MATERIAL_FLOW_SECTIONS;
use crate::graph::thermal::medium::series::hydraulic::pressure::BoussinesqDensity;
use crate::graph::thermal::transport::BranchThermalModel;
use crate::graph::{
    GraphBranch, GraphError, LinearLossElement, LinearLossResistance, LossGraph, PassiveBranch,
};
use crate::{AirflowError, LossElement, LossResistance, SourceProvenance, ToleranceBasis};

/// The original section law and its original source owners, before series
/// aggregation. An adiabatic section has an empty heat vector and range.
#[derive(Debug)]
pub struct DistributedMaterialSection<M> {
    pub law: MaterialBranchLaw,
    pub momentum: MaterialBranchEvaluation<M>,
    pub heat: Vec<MaterialSegmentEvaluation<M>>,
    /// Indices in this branch's concatenated native exchange model.
    pub exchange_range: Range<usize>,
}

/// Applied coefficients on exactly the original nodes, endpoints and branches.
/// Gravity remains unbound for the coupled owner's integrated-density profile.
#[derive(Debug)]
pub struct DistributedMaterialGraphResponse<M> {
    pub graph: LossGraph,
    pub models: Vec<BranchThermalModel>,
    /// Branch-major, then original declared hydraulic-section order.
    pub sections: Vec<Vec<DistributedMaterialSection<M>>>,
}

#[derive(Debug)]
pub enum DistributedResponseError {
    Invalid(&'static str),
    Branch {
        branch: usize,
        reason: &'static str,
    },
    /// Location in the original graph; the nested compiler error retains its
    /// own flattened section index and local heat-segment diagnostics.
    Section {
        branch: usize,
        section: usize,
        source: MaterialResponseError,
    },
    Compiler(MaterialResponseError),
}
impl core::fmt::Display for DistributedResponseError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "distributed material response: {self:?}")
    }
}
impl std::error::Error for DistributedResponseError {}
impl From<GraphError> for DistributedResponseError {
    fn from(error: GraphError) -> Self {
        Self::Compiler(MaterialResponseError::Graph(error))
    }
}
impl From<AirflowError> for DistributedResponseError {
    fn from(error: AirflowError) -> Self {
        Self::Compiler(MaterialResponseError::Airflow(error))
    }
}

/// Evaluate source-backed spatial sections and collapse their series drag onto
/// the same native graph. The material compiler is called once on a temporary
/// batch graph; that graph is NEVER pressure-solved or returned as a flow network.
///
/// Each original branch needs at least one hydraulic section. Every section
/// receives the branch's signed reference Q and an independently supplied bulk
/// material; every authored exchange receives its own supplied property state.
/// All shape/count checks precede material queries. Total hydraulic sections
/// AND total exchange segments are each bounded by MAX_MATERIAL_FLOW_SECTIONS.
///
/// The returned exchange models can be used directly in CirculationResponse.
/// Its callback can sample section materials from the actual segment means.
/// Adiabatic sections add no thermal degree of freedom: their representative
/// state, axial placement and any spatial refinement remain caller decisions.
///
/// Native summed coefficients carry an explicit EngineeringAllowance equal to
/// the largest supplied relative allowance among sections contributing to that
/// mechanism. This conservative positive-sum allowance assumes no statistical
/// independence and does not include material uncertainty or spatial error.
/// Original classifications and citations remain on every retained section.
///
/// Reference density owns continuity and inertia. Original material/correlation
/// support, Boussinesq departure bounds and forced-speed requirements all apply.
/// No averaging of material data, omitted-power deposition, new flow solve or
/// entrance-length inference is performed. The coupled owner must re-evaluate
/// the final solved state and enforce its own residual/work budgets.
#[allow(clippy::too_many_arguments)]
pub fn evaluate_distributed_material_graph<M: FluidTransport>(
    cx: &Cx<'_>,
    template: &LossGraph,
    density: BoussinesqDensity,
    laws: &[Vec<MaterialBranchLaw>],
    signed_volume_flows_m3_s: &[f64],
    bulk_materials: Vec<Vec<M>>,
    segment_materials: Vec<Vec<Vec<M>>>,
) -> Result<DistributedMaterialGraphResponse<M>, DistributedResponseError> {
    poll(cx)?;
    let count = template.branches().len();
    if count == 0
        || laws.len() != count
        || signed_volume_flows_m3_s.len() != count
        || bulk_materials.len() != count
        || segment_materials.len() != count
    {
        return Err(DistributedResponseError::Invalid(
            "one ordered section set and signed flow per original graph branch required",
        ));
    }
    let mut section_count = 0_usize;
    let mut heat_count = 0_usize;
    for (branch, sections) in laws.iter().enumerate() {
        poll(cx)?;
        if sections.is_empty()
            || !signed_volume_flows_m3_s[branch].is_finite()
            || bulk_materials[branch].len() != sections.len()
            || segment_materials[branch].len() != sections.len()
        {
            return Err(DistributedResponseError::Branch {
                branch,
                reason: "nonempty sections, finite flow and matching material counts required",
            });
        }
        section_count = bounded_count(section_count, sections.len())?;
        for (section, law) in sections.iter().enumerate() {
            poll(cx)?;
            heat_count = bounded_count(heat_count, law.exchange.len())?;
            if segment_materials[branch][section].len() != law.exchange.len() {
                return Err(DistributedResponseError::Branch {
                    branch,
                    reason: "each section needs exactly one material per authored exchange",
                });
            }
        }
    }

    let mut batch = Vec::with_capacity(section_count);
    let mut flat_laws = Vec::with_capacity(section_count);
    let mut flows = Vec::with_capacity(section_count);
    let mut bulk = Vec::with_capacity(section_count);
    let mut films = Vec::with_capacity(section_count);
    let mut addresses = Vec::with_capacity(section_count);
    for (branch, (branch_bulk, branch_films)) in bulk_materials
        .into_iter()
        .zip(segment_materials)
        .enumerate()
    {
        for (section, (medium, heat)) in branch_bulk.into_iter().zip(branch_films).enumerate() {
            poll(cx)?;
            // Copy the already valid template law solely for the compiler's
            // topology/name input. No placeholder resistance is fabricated.
            let name = format!("material-section-{branch}-{section}");
            let mut edge = template.branches()[branch].clone();
            edge.loss.name = name.clone();
            let compiled = if let Some(original) = &template.linear_losses()[branch] {
                let mut linear_loss = original.as_ref().clone();
                linear_loss.name = name;
                if edge.loss.resistance.value() > 0.0 {
                    PassiveBranch::LinearQuadratic {
                        branch: edge,
                        linear_loss,
                    }
                } else {
                    PassiveBranch::Linear {
                        from: edge.from,
                        to: edge.to,
                        loss: linear_loss,
                    }
                }
            } else {
                PassiveBranch::Quadratic(edge)
            };
            batch.push(compiled);
            flat_laws.push(laws[branch][section].clone());
            flows.push(signed_volume_flows_m3_s[branch]);
            bulk.push(medium);
            films.push(heat);
            addresses.push((branch, section));
        }
    }
    let batch = LossGraph::new_mixed(template.node_count(), batch)?;
    let MaterialGraphResponse {
        graph: compiled_graph,
        models,
        momentum,
        heat,
    } = evaluate_material_graph(cx, &batch, density, &flat_laws, &flows, bulk, films)
        .map_err(|source| locate(source, &addresses))?;
    // Neither batch topology represents the physical series connection. Only
    // independently evaluated section coefficients and thermal rows are used.
    drop(compiled_graph);
    drop(batch);

    let mut evaluated = flat_laws.into_iter().zip(momentum).zip(heat).zip(models);
    let mut graph_branches = Vec::with_capacity(count);
    let mut models = Vec::with_capacity(count);
    let mut retained = Vec::with_capacity(count);
    for (branch, original) in template.branches().iter().enumerate() {
        let mut linear = PositiveSum::default();
        let mut quadratic = PositiveSum::default();
        let mut linear_allowance = 0.0_f64;
        let mut quadratic_allowance = 0.0_f64;
        let mut segments = Vec::new();
        let mut sections = Vec::with_capacity(laws[branch].len());
        let mut citations = String::from("Series sum of original material-section laws: ");
        let mut identifiers = String::from("material-section-series:");
        for section in 0..laws[branch].len() {
            poll(cx)?;
            let (((law, momentum), heat), model) = evaluated.next().ok_or(
                DistributedResponseError::Invalid("material compiler returned incomplete sections"),
            )?;
            let a = momentum.linear_volume_resistance_pa_s_m3;
            let b = momentum.quadratic_volume_resistance_pa_s2_m6;
            linear.add(a, branch)?;
            quadratic.add(b, branch)?;
            if a > 0.0 {
                linear_allowance = linear_allowance.max(law.uncertainty_rel);
            }
            if b > 0.0 {
                quadratic_allowance = quadratic_allowance.max(law.uncertainty_rel);
            }
            if section > 0 {
                citations.push_str("; ");
                identifiers.push(';');
            }
            citations.push_str(&law.source.citation);
            // Length-delimited identifiers preserve the original boundaries.
            identifiers.push_str(&format!(
                "{}:{}",
                law.source.identifier.len(),
                law.source.identifier
            ));
            let start = segments.len();
            if let BranchThermalModel::Exchange(mut local) = model {
                segments.append(&mut local);
            }
            sections.push(DistributedMaterialSection {
                law,
                momentum,
                heat,
                exchange_range: start..segments.len(),
            });
        }
        let source = SourceProvenance::new(citations, identifiers);
        let a = linear.finish(branch)?;
        let b = quadratic.finish(branch)?;
        let linear_loss = if a > 0.0 {
            Some(LinearLossElement::new(
                original.loss.name.clone(),
                LinearLossResistance::new(a),
                linear_allowance,
                source.clone(),
                ToleranceBasis::EngineeringAllowance,
            )?)
        } else {
            None
        };
        let edge = if b > 0.0 {
            let edge = GraphBranch {
                from: original.from,
                to: original.to,
                loss: LossElement::new(
                    original.loss.name.clone(),
                    LossResistance::new(b),
                    quadratic_allowance,
                    source,
                    ToleranceBasis::EngineeringAllowance,
                )?,
            };
            match linear_loss {
                Some(linear_loss) => PassiveBranch::LinearQuadratic {
                    branch: edge,
                    linear_loss,
                },
                None => PassiveBranch::Quadratic(edge),
            }
        } else {
            PassiveBranch::Linear {
                from: original.from,
                to: original.to,
                loss: linear_loss.ok_or(DistributedResponseError::Branch {
                    branch,
                    reason: "series branch must retain an active resistance mechanism",
                })?,
            }
        };
        graph_branches.push(edge);
        models.push(if segments.is_empty() {
            BranchThermalModel::Adiabatic
        } else {
            BranchThermalModel::Exchange(segments)
        });
        retained.push(sections);
    }
    poll(cx)?;
    Ok(DistributedMaterialGraphResponse {
        graph: LossGraph::new_mixed(template.node_count(), graph_branches)?,
        models,
        sections: retained,
    })
}

fn bounded_count(current: usize, additional: usize) -> Result<usize, DistributedResponseError> {
    current.checked_add(additional).filter(|&count| count <= MAX_MATERIAL_FLOW_SECTIONS)
        .ok_or(DistributedResponseError::Invalid(
            "total hydraulic sections and total heat segments must each fit MAX_MATERIAL_FLOW_SECTIONS",
        ))
}

fn locate(source: MaterialResponseError, addresses: &[(usize, usize)]) -> DistributedResponseError {
    let flat = match &source {
        MaterialResponseError::Branch { branch, .. }
        | MaterialResponseError::Material { branch, .. }
        | MaterialResponseError::Resistance { branch, .. }
        | MaterialResponseError::Convection { branch, .. } => Some(*branch),
        _ => None,
    };
    match flat.and_then(|index| addresses.get(index)).copied() {
        Some((branch, section)) => DistributedResponseError::Section {
            branch,
            section,
            source,
        },
        None => DistributedResponseError::Compiler(source),
    }
}

// Compensated positive sums preserve small spatial contributions without
// claiming an outward-rounded enclosure or allowing an overflowed branch.
#[derive(Default)]
struct PositiveSum {
    sum: f64,
    correction: f64,
}
impl PositiveSum {
    fn add(&mut self, value: f64, branch: usize) -> Result<(), DistributedResponseError> {
        let next = self.sum + value;
        let error = if self.sum >= value {
            (self.sum - next) + value
        } else {
            (value - next) + self.sum
        };
        self.correction += error;
        self.sum = next;
        if !self.sum.is_finite() || !self.correction.is_finite() {
            return Err(DistributedResponseError::Branch {
                branch,
                reason: "summed material resistance is not representable",
            });
        }
        Ok(())
    }
    fn finish(self, branch: usize) -> Result<f64, DistributedResponseError> {
        let value = self.sum + self.correction;
        if value.is_finite() && value >= 0.0 {
            Ok(value)
        } else {
            Err(DistributedResponseError::Branch {
                branch,
                reason: "summed material resistance is not representable",
            })
        }
    }
}
fn poll(cx: &Cx<'_>) -> Result<(), DistributedResponseError> {
    cx.checkpoint()
        .map_err(|_| DistributedResponseError::Compiler(MaterialResponseError::Interrupted))
}

#[cfg(test)]
mod tests;
