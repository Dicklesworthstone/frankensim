//! Original section hydraulic losses deposited into finite liquid storage.
//!
//! This opt-in owns the conversion of actual native `dp * Q` into fluid heat.
//! The original source law must have constant density over its caloric span.
//! Every original section is assumed to deposit its entire irreversible loss
//! into its own represented fluid cell. A porous skeleton or separate body
//! receiving a share needs an explicit partition model outside this adapter.
//!
//! Inertia, gravity, pressure-port work, pump inefficiency and numerical losses
//! are not thermal sources here. Original source and native mechanical-work
//! applicability caps still apply. The enclosing native thermal/momentum driver
//! owns integration and acceptance; evaluate this at each actual candidate.

use super::*;
use crate::graph::thermal::medium::series::enthalpy::admit_dissipative_density;
use fs_material::liquid::ResolvedLiquidState;
use fs_material::phase::ResolvedHeatCapacityCurve;

/// Original source refusal, conflicting ownership or a mechanical/thermal
/// transfer mismatch. A nonlinear solve tolerance cannot relax this transfer.
#[derive(Debug)]
pub enum GraphHeatingError {
    Source(SourceGraphError),
    CellSourceAlreadyOwned(CellId),
    TransferBalance {
        mechanical_loss_j: f64,
        thermal_heat_j: f64,
        residual_j: f64,
        allowed_j: f64,
    },
}
impl core::fmt::Display for GraphHeatingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite liquid graph hydraulic heating: {self:?}")
    }
}
impl std::error::Error for GraphHeatingError {}

/// Local native drag heat in original branch/section order, including zeros.
#[derive(Clone, Copy, Debug)]
pub struct GraphHeatingSample {
    pub cell: CellId,
    pub branch: usize,
    pub section: usize,
    pub power_w: f64,
    pub heat_j: f64,
}

/// Prepared borrowing of the original liquid graph source and caloric binding.
/// No new material, flow law, time history or thermal integrator is created.
pub struct LiquidGraphHeating<'s, 'm> {
    core: &'s BoundGraphCore<'m>,
}

/// Additive native cell powers plus their original source and transfer receipt.
/// The source sections share the exact immutable epoch used by GraphSourceState.
#[derive(Debug)]
pub struct ResolvedLiquidGraphHeating {
    sources_w: BTreeMap<CellId, f64>,
    samples: Vec<GraphHeatingSample>,
    sections: Arc<Vec<Vec<DistributedMaterialSection<ResolvedLiquidState>>>>,
    admission: GraphSourceAdmission,
    duration_s: f64,
    mechanical_loss_j: f64,
    thermal_heat_j: f64,
    transfer_residual_j: f64,
}
impl ResolvedLiquidGraphHeating {
    /// Additional powers for the native thermal residual. Zero entries are
    /// omitted; pre-existing authored boundary sources remain additive.
    pub fn sources_w(&self) -> &BTreeMap<CellId, f64> {
        &self.sources_w
    }
    pub fn samples(&self) -> &[GraphHeatingSample] {
        &self.samples
    }
    pub fn source_sections(&self) -> &[Vec<DistributedMaterialSection<ResolvedLiquidState>>] {
        &self.sections
    }
    /// Original source support and conservative old/new/average speed bounds.
    /// Its envelope powers are applicability limits, never deposited heat.
    pub fn admission(&self) -> &GraphSourceAdmission {
        &self.admission
    }
    pub fn duration_s(&self) -> f64 {
        self.duration_s
    }
    pub fn mechanical_loss_j(&self) -> f64 {
        self.mechanical_loss_j
    }
    pub fn thermal_heat_j(&self) -> f64 {
        self.thermal_heat_j
    }
    pub fn transfer_residual_j(&self) -> f64 {
        self.transfer_residual_j
    }

    /// Apply once to StepCoefficients' additional cell-power map.
    ///
    /// A pre-existing key refuses rather than silently replacing another
    /// source. An application composing multiple physical sources can sum the
    /// read-only source maps explicitly. Cancellation or conflicting ownership
    /// leaves the destination unchanged.
    pub fn apply_to(
        &self,
        cx: &Cx<'_>,
        destination: &mut BTreeMap<CellId, f64>,
    ) -> Result<(), GraphHeatingError> {
        for &cell in self.sources_w.keys() {
            poll(cx).map_err(GraphHeatingError::Source)?;
            if destination.contains_key(&cell) {
                return Err(GraphHeatingError::CellSourceAlreadyOwned(cell));
            }
        }
        poll(cx).map_err(GraphHeatingError::Source)?;
        destination.extend(self.sources_w.iter().map(|(&cell, &power)| (cell, power)));
        Ok(())
    }
}

impl<'s, 'm> LiquidGraphHeating<'s, 'm> {
    pub(super) fn new(
        cx: &Cx<'_>,
        core: &'s BoundGraphCore<'m>,
        caloric: &ResolvedHeatCapacityCurve,
    ) -> Result<Self, GraphHeatingError> {
        poll(cx).map_err(GraphHeatingError::Source)?;
        admit_dissipative_density(caloric, core.reference.density_kg_m3())
            .map_err(|error| GraphHeatingError::Source(SourceGraphError::LiquidCaloric(error)))?;
        poll(cx).map_err(GraphHeatingError::Source)?;
        Ok(Self { core })
    }

    /// Evaluate inside the native inertial driver's exchange callback, using
    /// its actual average flows and same-temperature original source payload.
    /// The hydraulic loss debit must equal the local heat transfer to roundoff.
    pub fn evaluate_inertial(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'_>,
        sources: &GraphSourceState<ResolvedLiquidState>,
        report: &InertialGraphReport,
    ) -> Result<ResolvedLiquidGraphHeating, GraphHeatingError> {
        self.evaluate(cx, thermal, sources, GraphFlow::Inertial(report))
    }

    /// Deposit the solved endpoint's actual drag over the declared thermal
    /// interval. This is the native quasi-steady model, without inertia history.
    pub fn evaluate_quasi_steady(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'_>,
        sources: &GraphSourceState<ResolvedLiquidState>,
        hydraulics: &GraphSolution,
        duration_s: f64,
    ) -> Result<ResolvedLiquidGraphHeating, GraphHeatingError> {
        self.evaluate(
            cx,
            thermal,
            sources,
            GraphFlow::QuasiSteady { hydraulics, duration_s },
        )
    }

    fn evaluate(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'_>,
        sources: &GraphSourceState<ResolvedLiquidState>,
        flow: GraphFlow<'_>,
    ) -> Result<ResolvedLiquidGraphHeating, GraphHeatingError> {
        let source_error = GraphHeatingError::Source;
        // This retains exact source identity and thermal state, original
        // resistance support at all native rates, topology and all work caps.
        let admission = self.core.admit_flow(cx, thermal, sources, flow).map_err(source_error)?;
        let density = self.core.reference.density_kg_m3();
        for material in &sources.node_materials {
            poll(cx).map_err(source_error)?;
            let transport = material.transport_state()
                .map_err(|error| source_error(SourceGraphError::Material(error)))?;
            if transport.density_kg_m3() != density {
                return Err(source_error(SourceGraphError::Invalid(
                    "hydraulic heat requires the original constant liquid density",
                )));
            }
        }

        let duration_s = flow.duration_s();
        let section_count: usize = self.core.inventory.branch_cells().iter().map(Vec::len).sum();
        let mut samples = Vec::with_capacity(section_count);
        let mut sources_w = BTreeMap::new();
        let mut thermal_heat_j = 0.0;
        let mut native_branch_loss_j = 0.0;
        for (branch, cells) in self.core.inventory.branch_cells().iter().enumerate() {
            poll(cx).map_err(source_error)?;
            let branch_power_w = flow.pressure_drop(branch) * flow.average(branch);
            let branch_heat_j = duration_s * branch_power_w;
            if !nonnegative(branch_power_w) || !nonnegative(branch_heat_j)
                || (branch_power_w > 0.0 && branch_heat_j == 0.0)
            {
                return Err(source_error(SourceGraphError::Invalid(
                    "native hydraulic heat is unrepresentable",
                )));
            }
            native_branch_loss_j += branch_heat_j;
            for (section, &cell) in cells.iter().enumerate() {
                poll(cx).map_err(source_error)?;
                let source = &sources.sections[branch][section];
                if source.momentum.evaluation.transport.density_kg_m3() != density {
                    return Err(source_error(SourceGraphError::Cell {
                        cell,
                        reason: "hydraulic heat requires the original constant liquid density",
                    }));
                }
                let power_w = admission.actual_losses[branch][section].dissipated_power_w;
                let heat_j = duration_s * power_w;
                if !nonnegative(power_w) || !nonnegative(heat_j)
                    || (power_w > 0.0 && heat_j == 0.0)
                {
                    return Err(source_error(SourceGraphError::Cell {
                        cell,
                        reason: "local hydraulic heat is unrepresentable",
                    }));
                }
                thermal_heat_j += heat_j;
                if power_w > 0.0 {
                    // The original distributed inventory admits unique pipe
                    // cells; no junction receives fictitious branch loss heat.
                    sources_w.insert(cell, power_w);
                }
                samples.push(GraphHeatingSample { cell, branch, section, power_w, heat_j });
            }
        }
        let mechanical_loss_j = match flow {
            GraphFlow::Inertial(report) => report.hydraulic_loss_j,
            GraphFlow::QuasiSteady { .. } => native_branch_loss_j,
        };
        if [mechanical_loss_j, native_branch_loss_j, thermal_heat_j]
            .into_iter().any(|energy| !nonnegative(energy))
        {
            return Err(source_error(SourceGraphError::Invalid(
                "finite nonnegative native hydraulic debit and thermal credit required",
            )));
        }
        let allowed_j = 256.0 * f64::EPSILON
            * (1.0 + section_count as f64 + admission.actual_losses.len() as f64)
            * mechanical_loss_j.max(native_branch_loss_j).max(thermal_heat_j)
                .max(f64::MIN_POSITIVE);
        for credit in [native_branch_loss_j, thermal_heat_j] {
            let residual_j = credit - mechanical_loss_j;
            if !residual_j.is_finite() || residual_j.abs() > allowed_j {
                return Err(GraphHeatingError::TransferBalance {
                    mechanical_loss_j,
                    thermal_heat_j: credit,
                    residual_j,
                    allowed_j,
                });
            }
        }
        poll(cx).map_err(source_error)?;
        Ok(ResolvedLiquidGraphHeating {
            sources_w,
            samples,
            sections: Arc::clone(&sources.sections),
            admission,
            duration_s,
            mechanical_loss_j,
            thermal_heat_j,
            transfer_residual_j: thermal_heat_j - mechanical_loss_j,
        })
    }
}
