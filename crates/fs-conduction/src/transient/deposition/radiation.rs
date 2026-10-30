//! Closed-enclosure radiation on the shared finite-cell enthalpy history.
//! Reuse the original material resolver, radiosity equations, thermal residual
//! and nonlinear solver. Photons transfer heat between cells; they create no
//! second thermal store and no net external source in a closed enclosure.

use std::cell::Cell;
use std::collections::BTreeMap;
use fs_exec::Cx;
use fs_material::phase::{EquilibriumEnthalpyPhaseCurve, EquilibriumPhaseState, PhaseStateError};
use fs_material::state_point::{MaterialPropertySelection, ResolvedMaterialStatePoint};
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_solver::{NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport, SolverRunProgress};
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_thermal::enthalpy::{CellId, MaterialId, ThermalLaw, DepositionSimulation, DepositionError, EnergyAudit};
use fs_thermal::enthalpy::backward_euler::{BackwardEulerStep, StepCoefficients};
use crate::{ConductionError, GrayDiffuseEnclosure, RadiosityReport, SurfaceEmissivity};

/// One original enclosure surface and its existing finite thermal cell.
/// Several disjoint surface traces may belong to one uniform cell.
pub struct RadiatingCell<'a> {
    /// Original enclosure surface/matrix index, not a thermal cell index.
    pub surface: usize,
    /// Active receiving cell. Passive nonradiating cells need no binding.
    pub cell: CellId,
    /// Exact material-table identity retained by the cell.
    pub material_id: MaterialId,
    /// SAME borrowed caloric object used to initialize the thermal timeline.
    pub curve: &'a EquilibriumEnthalpyPhaseCurve,
    /// Same-card thermal/phase coordinates, preserving other physical conditions.
    pub material: ThermalMaterialPoint<'a>,
    /// Explicit total-emissivity property selection, including exact pins.
    pub selection: MaterialPropertySelection,
}

/// Complete finite-body binding of a closed, opaque gray-diffuse enclosure.
/// Geometry/view factors stay fixed; emissivities are re-resolved at each trial.
pub struct MaterialRadiation<'a> {
    /// Original admitted trace geometry and reciprocal view factors. Its initial
    /// emissivities are replaced completely, never used as missing-data defaults.
    pub enclosure: &'a GrayDiffuseEnclosure,
    /// Exactly one active-cell binding per surface, in any order.
    pub cells: &'a [RadiatingCell<'a>],
    /// Absolute net-enclosure and deposited-power discrepancy allowance [W].
    pub power_tolerance_w: f64,
    /// Absolute largest original radiosity equation residual [W/m2].
    pub radiosity_tolerance_w_m2: f64,
}

/// Original source and solver failures; none returns a partial thermal state.
#[derive(Debug)]
pub enum RadiationCellError {
    /// Incomplete binding, geometry, controls or unrepresentable arithmetic.
    Invalid(&'static str),
    /// Original caloric chart refused a phase state.
    Phase(PhaseStateError),
    /// Original material/radiosity/geometry refusal.
    Radiation(ConductionError),
    /// Original thermal state, topology or energy refusal.
    Thermal(DepositionError),
    /// Original Newton construction refusal.
    Newton(NewtonError),
    /// Original nonlinear driver did not converge.
    NotConverged(NewtonReport),
    /// Numerical interruption other than explicit cancellation.
    Interrupted(SolverRunProgress),
    /// All-trial radiation work budget, including the final fresh evaluation.
    EvaluationBudget,
    /// Cancellation before publication.
    Cancelled,
}
impl From<ConductionError> for RadiationCellError {
    fn from(e: ConductionError) -> Self {
        match e { ConductionError::Cancelled { .. } => Self::Cancelled, other => Self::Radiation(other) }
    }
}
impl From<DepositionError> for RadiationCellError {
    fn from(e: DepositionError) -> Self {
        if e == DepositionError::Cancelled { Self::Cancelled } else { Self::Thermal(e) }
    }
}
impl core::fmt::Display for RadiationCellError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result { write!(f,"finite-cell radiation: {self:?}") }
}
impl std::error::Error for RadiationCellError {}
fn checkpoint(cx: &Cx<'_>) -> Result<(), DepositionError> {
    cx.checkpoint().map_err(|_| DepositionError::Cancelled)
}
fn numerical(error: FiniteDifferenceError<RadiationCellError>) -> RadiationCellError {
    match error { FiniteDifferenceError::Source(e) => e, FiniteDifferenceError::Invalid(s) => RadiationCellError::Invalid(s) }
}
fn nonnegative(x:f64)->bool { x.is_finite() && x>=0.0 }

/// Original phase and full material query for one surface in matrix order.
#[derive(Debug, Clone, PartialEq)]
pub struct RadiatingCellSample {
    /// Original enclosure index.
    pub surface: usize,
    /// Receiving finite thermal cell.
    pub cell: CellId,
    /// Actual caloric endpoint, not averaged temperatures.
    pub phase: EquilibriumPhaseState,
    /// Complete original total-emissivity source bundle and receipts.
    pub material: ResolvedMaterialStatePoint,
}

/// A tentative material/heat evaluation, with no thermal history advanced.
#[derive(Debug, Clone, PartialEq)]
pub struct RadiationCellSources {
    sources_w: BTreeMap<CellId,f64>,
    surfaces: Vec<RadiatingCellSample>,
    enclosure: GrayDiffuseEnclosure,
    radiosity: RadiosityReport,
}
impl RadiationCellSources {
    /// Signed watts INTO each receiving cell, accumulated over its surfaces.
    pub fn sources_w(&self)->&BTreeMap<CellId,f64> { &self.sources_w }
    /// Actual accepted material/phase evaluations in original matrix order.
    pub fn surfaces(&self)->&[RadiatingCellSample] { &self.surfaces }
    /// Exact original radiosity solve, including its independent residuals.
    pub const fn radiosity(&self)->&RadiosityReport { &self.radiosity }
    /// Geometry plus the fresh emissivity snapshot used by that solve.
    pub const fn enclosure(&self)->&GrayDiffuseEnclosure { &self.enclosure }
    /// Sum of radiation watts handed to cells. Ideally zero in a closed enclosure;
    /// retain its roundoff rather than secretly redistributing or clipping heat.
    pub fn net_heat_w(&self)->f64 { self.sources_w.values().copied().sum() }

    /// Conservative outgoing-rate scale for this frozen emissivity snapshot.
    ///
    /// The original gray-diffuse equations give
    /// P = F [I-(I-E)F]^-1 E and q_i = A_i epsilon_i (b_i-(P b)_i),
    /// where b_i = sigma T_i^4. For exactly closed rows, P is nonnegative
    /// with unit row sums, so each outgoing secant sum is at most
    /// 4 sigma A_i epsilon_i Tmax^3. Multiple surfaces of one cell add.
    /// The same net-radiation equations are described by the MOOSE heat-transfer
    /// theory: <https://mooseframework.inl.gov/modules/heat_transfer/>.
    ///
    /// Actual admitted row defects are retained, not normalized. For row sum s,
    /// m = epsilon + (1-epsilon)(1-s) must be positive. Then the nonnegative
    /// inverse obeys [I-(I-E)F]^-1 E 1 <= c 1, c=max(epsilon/m).
    /// Multiplying each surface bound by max(1,s*c) also bounds its positive
    /// diagonal loss and off-diagonal row scale. A nonpositive margin refuses
    /// when this sufficient reflection bound cannot be established.
    ///
    /// Pass these rates to the native thermal additional-conductance limiter;
    /// they are NOT physical heat links or additional source powers. The
    /// original radiation source and its closed-energy gate still own heat.
    /// This coefficient-frozen estimate is not a temporal-error bound, a bound
    /// on d(epsilon)/dT, or a maximum-principle/passivity claim for defective
    /// view factors or independent forcing. All temperatures and emissivities
    /// must be refreshed on each accepted explicit substep.
    pub fn conductance_bounds_w_k(
        &self,
        mut check: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, RadiationCellError> {
        check()?;
        let invalid = RadiationCellError::Invalid;
        let mut maximum_temperature = 0.0_f64;
        let mut amplification = 0.0_f64;
        for (sample, surface) in self.surfaces.iter().zip(self.enclosure.surfaces()) {
            check()?;
            maximum_temperature = maximum_temperature.max(sample.phase.temperature_k());
            let epsilon = surface.emissivity().value();
            let row_sum = self.radiosity.view_factor_row_sums[sample.surface];
            // This form preserves epsilon for an exactly closed row even when
            // 1-epsilon rounds to one for a strongly reflecting surface.
            let margin = epsilon + (1.0 - epsilon) * (1.0 - row_sum);
            let ratio = epsilon / margin;
            if !margin.is_finite() || margin <= 0.0 || !ratio.is_finite() || ratio <= 0.0 {
                return Err(invalid("radiation timestep needs a positive bounded reflection margin"));
            }
            amplification = amplification.max(ratio);
        }
        let emission_slope = ((4.0 * crate::STEFAN_BOLTZMANN_W_M2_K4
            * maximum_temperature) * maximum_temperature) * maximum_temperature;
        if !emission_slope.is_finite() || emission_slope <= 0.0 {
            return Err(invalid("radiation temperature-rate scale is unrepresentable"));
        }
        let mut bounds = BTreeMap::new();
        for (sample, surface) in self.surfaces.iter().zip(self.enclosure.surfaces()) {
            check()?;
            let row_scale = self.radiosity.view_factor_row_sums[sample.surface] * amplification;
            let conductance = ((surface.area_m2() * surface.emissivity().value())
                * emission_slope) * row_scale.max(1.0);
            if !row_scale.is_finite() || !conductance.is_finite() || conductance <= 0.0 {
                return Err(invalid("radiation conductance bound is unrepresentable"));
            }
            let sum = bounds.entry(sample.cell).or_insert(0.0);
            *sum += conductance;
            if !sum.is_finite() {
                return Err(invalid("combined radiative surface bound is unrepresentable"));
            }
        }
        check()?;
        Ok(bounds)
    }

    /// Add radiation to another physical owner's source map, preserving its
    /// independent sources. All additions are checked before modifying the map.
    /// Use once per trial; never also deposit the returned surface heats.
    pub fn add_to(&self, destination:&mut BTreeMap<CellId,f64>)->Result<(),RadiationCellError> {
        let mut staged=Vec::with_capacity(self.sources_w.len());
        for (&id,&heat) in &self.sources_w {
            let prior=destination.get(&id).copied().unwrap_or(0.0);
            let total=prior+heat;
            if !prior.is_finite() || !total.is_finite() || (heat!=0.0 && total==prior) {
                return Err(RadiationCellError::Invalid("radiation disappears or overflows in source addition"));
            }
            staged.push((id,total));
        }
        for (id,value) in staged { destination.insert(id,value); }
        Ok(())
    }
}

impl MaterialRadiation<'_> {
    /// Resolve each surface at its receiving cell's actual T/phase and run the
    /// original radiosity solver. The caller must justify mapping mesh traces to
    /// uniform cells and physical non-overlap. Summed radiative area must fit each
    /// cell's exposed area after existing thermal contacts/reservoir coverage.
    /// Radiation and ambient convection may coexist on that exposed surface:
    /// they are distinct mechanisms, not duplicate convective boundaries.
    ///
    /// Every surface must name an active cell and its actual immutable caloric
    /// object. Future cells, missing emissivity, and foreign charts refuse even
    /// at equilibrium. The original 256-surface bound applies before allocation.
    pub fn resolve(&self,cx:&Cx<'_>,thermal:&DepositionSimulation<'_>)
        ->Result<RadiationCellSources,RadiationCellError> {
        checkpoint(cx)?;
        let invalid=RadiationCellError::Invalid;
        let n=self.enclosure.surfaces().len();
        if self.cells.len()!=n || !nonnegative(self.power_tolerance_w)
            || !nonnegative(self.radiosity_tolerance_w_m2) {
            return Err(invalid("complete radiative bindings and finite nonnegative watt tolerances required"));
        }
        let mut ordered=vec![None;n]; let mut areas=BTreeMap::new();
        for binding in self.cells {
            checkpoint(cx)?;
            if binding.surface>=n || ordered[binding.surface].is_some() {
                return Err(invalid("duplicate or foreign enclosure surface index"));
            }
            let cell=thermal.cells().get(&binding.cell).ok_or(invalid("inactive radiating cell"))?;
            if cell.material_id!=binding.material_id
                || !thermal.uses_material_law(binding.material_id,binding.curve as &dyn ThermalLaw) {
                return Err(invalid("radiating cell must retain its original caloric object and material ID"));
            }
            let area=areas.entry(binding.cell).or_insert(0.0);
            *area+=self.enclosure.surfaces()[binding.surface].area_m2();
            if !area.is_finite() || *area>thermal.exposed_area_m2(binding.cell)? {
                return Err(invalid("radiating traces exceed the cell's exposed area"));
            }
            ordered[binding.surface]=Some(binding);
        }
        let mut emissivities=Vec::with_capacity(n); let mut temperatures=Vec::with_capacity(n);
        let mut samples=Vec::with_capacity(n);
        for (surface,binding) in ordered.into_iter().enumerate() {
            checkpoint(cx)?;
            let binding=binding.ok_or(invalid("unbound enclosure surface"))?;
            let h=thermal.cells()[&binding.cell].specific_h_j_kg;
            let phase=binding.curve.state_at_specific_enthalpy(h).map_err(RadiationCellError::Phase)?;
            let (emissivity,material)=SurfaceEmissivity::from_thermal_state(
                self.enclosure.surfaces()[surface].name(),binding.material,phase,binding.selection.clone())?;
            temperatures.push(phase.temperature_k()); emissivities.push(emissivity);
            samples.push(RadiatingCellSample {surface,cell:binding.cell,phase,material});
        }
        let enclosure=self.enclosure.with_emissivities(emissivities)?;
        let report=enclosure.solve(cx,&temperatures)?;
        if report.net_outward_heat_w.iter().chain(&report.radiosity_w_m2).chain(&report.irradiation_w_m2)
            .any(|v|!v.is_finite()) || !report.enclosure_energy_closure_w.is_finite()
            || !nonnegative(report.linear_residual_max_w_m2)
            || report.linear_residual_max_w_m2>self.radiosity_tolerance_w_m2
            || report.enclosure_energy_closure_w.abs()>self.power_tolerance_w {
            return Err(invalid("original radiosity equations or closed-enclosure power balance did not close"));
        }
        let mut sources=BTreeMap::new();
        for (sample,&outward) in samples.iter().zip(&report.net_outward_heat_w) {
            checkpoint(cx)?;
            let value=sources.entry(sample.cell).or_insert(0.0);
            *value-=outward;
            if !value.is_finite() { return Err(invalid("accumulated radiation cell power is nonfinite")); }
        }
        let total:f64=sources.values().copied().sum();
        if !total.is_finite() || total.abs()>self.power_tolerance_w
            || (total+report.enclosure_energy_closure_w).abs()>self.power_tolerance_w {
            return Err(invalid("surface-to-cell power partition did not close"));
        }
        checkpoint(cx)?;
        Ok(RadiationCellSources {sources_w:sources,surfaces:samples,enclosure,radiosity:report})
    }
}

/// Numerical work and independent physical gates for ONE implicit interval.
#[derive(Debug,Clone,Copy,PartialEq)]
pub struct RadiationStepSolve {
    /// Maximum active thermal cells, including nonradiating cells.
    pub maximum_cells:usize,
    /// All radiosity/material evaluations, including stencils and final admission.
    pub maximum_evaluations:usize,
    /// Original Newton outer work limit.
    pub maximum_iterations:usize,
    /// Maximum component finite-difference perturbation in cell joules.
    pub difference_step_j:f64,
    /// Absolute largest local thermal energy residual [J].
    pub cell_tolerance_j:f64,
    /// Absolute total thermal and nonradiative-external energy residual [J].
    pub total_tolerance_j:f64,
    /// Original numerical solver policy; residual coordinates are joules.
    pub newton:NewtonKrylovConfig,
}
/// Original source, nonlinear and physical accounts for the returned endpoint.
pub struct RadiationStepReport {
    /// Actual endpoint radiative source; its heat is already in the cell state.
    pub radiation:RadiationCellSources,
    /// Original thermal account (external_heat_j includes roundoff-net radiation).
    pub thermal:EnergyAudit,
    /// Authored sources, convection and reservoirs EXCLUDING internal radiation [J].
    pub external_heat_j:f64,
    /// Actual thermal storage change minus nonradiative external heat [J].
    pub energy_residual_j:f64,
    /// Fresh largest local residual admitted after Newton [J].
    pub maximum_cell_residual_j:f64,
    /// Original numerical convergence report; not a time-error certificate.
    pub nonlinear:NewtonReport,
    /// Actual material/radiosity work, including probes and final evaluation.
    pub evaluations:usize,
}
impl MaterialRadiation<'_> {
    /// Advance finite-cell enthalpy with endpoint-evaluated radiation, conduction,
    /// authored heat and reservoirs. Reuses BackwardEulerStep and Newton--Krylov;
    /// no new heat equation, phase model or radiation law is supplied here.
    /// Events inside the interval refuse; process start events first and target
    /// events after joint acceptance. Geometry, view factors and mass stay fixed.
    ///
    /// Unknowns/guesses are cell energy increments [J] in sorted CellId order.
    /// Trials never advance accepted history. A fresh final source, local thermal
    /// gates, and a closed-enclosure total-energy gate precede publication.
    /// First-order time accuracy in general; convergence does not prove temporal
    /// accuracy, uniqueness or unconditional stability. Refine independently.
    pub fn step_candidate<'m>(&self,cx:&Cx<'_>,thermal:&DepositionSimulation<'m>,target_s:f64,
        controls:RadiationStepSolve,guess_j:&[f64])
        ->Result<(DepositionSimulation<'m>,RadiationStepReport),RadiationCellError> {
        checkpoint(cx)?;
        if controls.maximum_cells==0 || thermal.cells().len()>controls.maximum_cells
            || controls.maximum_evaluations==0 || controls.maximum_iterations==0
            || guess_j.len()!=thermal.cells().len() || !nonnegative(controls.cell_tolerance_j)
            || !nonnegative(controls.total_tolerance_j) {
            return Err(RadiationCellError::Invalid("bounded thermal unknowns/work and finite joule tolerances required"));
        }
        controls.newton.validate().map_err(RadiationCellError::Newton)?;
        let step=BackwardEulerStep::new(thermal,target_s,||checkpoint(cx))?;
        let evaluations=Cell::new(0usize);
        let evaluate=|x:&[f64]| step.evaluate(x,|state| {
            if evaluations.get()>=controls.maximum_evaluations { return Err(RadiationCellError::EvaluationBudget); }
            evaluations.set(evaluations.get()+1);
            let source=self.resolve(cx,state)?;
            Ok::<_,RadiationCellError>(StepCoefficients { sources_w:source.sources_w.clone(),
                conductances_w_k:BTreeMap::new(), payload:source })
        },||checkpoint(cx));
        let problem=FallibleFiniteDifference::new(step.cell_ids().len(),controls.difference_step_j,
            |x:&[f64]|evaluate(x).map(|t|t.residual_j().to_vec())).map_err(numerical)?;
        let construction=NewtonKrylovState::new(&problem,guess_j.to_vec(),controls.newton);
        if let Some(e)=problem.take_failure() { return Err(numerical(e)); }
        let mut solver=construction.map_err(RadiationCellError::Newton)?;
        let run=solver.run_controlled(&problem,controls.maximum_iterations,
            &mut ||cx.checkpoint().is_err() || problem.has_failed());
        if let Some(e)=problem.take_failure() { return Err(numerical(e)); }
        checkpoint(cx)?;
        if run.progress!=SolverRunProgress::Complete { return Err(RadiationCellError::Interrupted(run.progress)); }
        if !run.report.converged { return Err(RadiationCellError::NotConverged(run.report)); }
        let trial=evaluate(&solver.x)?;
        let maximum_cell_residual_j=trial.maximum_cell_residual_j();
        let (next,radiation,audit)=trial.accept(controls.cell_tolerance_j,controls.total_tolerance_j,||checkpoint(cx))?;
        let external_heat_j=audit.external_heat_j-step.duration_s()*radiation.net_heat_w();
        let residual=audit.stored_energy_change_j-external_heat_j;
        if !external_heat_j.is_finite() || !residual.is_finite() || residual.abs()>controls.total_tolerance_j {
            return Err(RadiationCellError::Invalid("finite thermal stores and nonradiative supply do not close"));
        }
        checkpoint(cx)?;
        Ok((next,RadiationStepReport {radiation,thermal:audit,external_heat_j,energy_residual_j:residual,
            maximum_cell_residual_j,nonlinear:run.report,evaluations:evaluations.get()}))
    }
}
