//! Spatial Maxwell loss and reversible heat on the original gravity/dead-load beam.
//!
//! Every axial-cell/layer has its own thermal destination and original material
//! history. The existing bridge owner integrates compatible displacement and
//! independent load work. The shared enthalpy residual and Newton--Krylov owner
//! solve all cell energies together; no section is advanced sequentially in time.
//! Fixed thermal geometry and the CaloricPoint thermodynamic restrictions apply.
//! This path does not add drawing tension, moving contact, inertia or birth events.

use std::cell::{Cell, RefCell};
use std::collections::{BTreeMap, BTreeSet};
use fs_material::thermovisco::{Energy, thermal::CaloricPoint};
use fs_solid::bridge as solid;
use fs_solver::{NewtonError, NewtonKrylovState, NewtonReport, SolverRunProgress};
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation, EnergyAudit};
use fs_thermal::enthalpy::backward_euler::{BackwardEulerStep, StepCoefficients};
use super::{LayerThermalBinding, StripEnthalpyConfig, StripEnthalpyError, ThermalCoefficients};
use crate::circuit_enthalpy::material::{MaterialNetworkSample, PreparedMaterialNetwork};
use super::transport::TransportError;
use super::exchange::{LayerThermalExchange, layer_exchange, deposit_heat_j};

/// Atomic spatial tip births after a coupled interval.
pub mod growth;

/// Shared finite electrical storage, Joule heat and spatial material mechanics.
pub mod electrical;

/// Mechanical tuple and its actual source-bound thermal destination.
#[derive(Clone, Copy, Debug)]
pub struct BridgeThermalBinding<'a> {
    pub axial_cell: u64,
    pub section: LayerThermalBinding<'a>,
}

/// Original section/nonlinear budgets plus the integrated beam's deflection gate.
#[derive(Clone, Copy, Debug)]
pub struct BridgeEnthalpyConfig {
    pub feedback: StripEnthalpyConfig,
    pub maximum_deflection_m: f64,
}

#[derive(Debug)]
pub enum BridgeEnthalpyError<E> {
    Coupled(StripEnthalpyError<E>),
    Mechanical(solid::Error),
}
impl<E> From<DepositionError> for BridgeEnthalpyError<E> {
    fn from(e: DepositionError) -> Self { Self::Coupled(e.into()) }
}
impl<E> From<NewtonError> for BridgeEnthalpyError<E> {
    fn from(e: NewtonError) -> Self { Self::Coupled(e.into()) }
}
impl<E> From<solid::Error> for BridgeEnthalpyError<E> {
    fn from(e: solid::Error) -> Self {
        if e == solid::Error::Cancelled { Self::Coupled(StripEnthalpyError::Cancelled) }
        else { Self::Mechanical(e) }
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for BridgeEnthalpyError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "bridge/enthalpy coupling: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for BridgeEnthalpyError<E> {}
fn invalid<E>(s: &'static str) -> BridgeEnthalpyError<E> {
    BridgeEnthalpyError::Coupled(StripEnthalpyError::Invalid(s))
}
fn finite(value: f64) -> Result<f64, DepositionError> {
    if value.is_finite() { Ok(value) } else { Err(DepositionError::Arithmetic) }
}
fn same(a: f64, b: f64) -> bool {
    a.is_finite() && b.is_finite() && a > 0.0 && b > 0.0
        && (a-b).abs() <= 64.0*f64::EPSILON*a.abs().max(b.abs())
}

// The native bridge already integrated every point's free-energy/work/loss
// receipt. Reuse the original quadrature thermodynamics for both this path and
// the circuit-coupled path; geometry, histories and thermal destinations coincide.
struct BridgeThermalExchange {
    layers: BTreeMap<(u64,u64),LayerThermalExchange>,
    mechanical_internal_energy_change_j: f64,
    mechanical_entropy_energy_change_j: f64,
    reversible_heat_j: f64,
    dissipated_j: f64,
}
fn bridge_exchange<E>(
    old_heat:&DepositionSimulation<'_>, next_heat:&DepositionSimulation<'_>,
    old:&solid::Bridge<'_>, next:&solid::Bridge<'_>,
    bound:&BTreeMap<(u64,u64),&BridgeThermalBinding<'_>>,
    energy:&BTreeMap<(u64,u64),Energy>,
    mut checkpoint:impl FnMut()->Result<(),DepositionError>,
) -> Result<BridgeThermalExchange,StripEnthalpyError<E>> {
    if old.cells().len()!=next.cells().len() || energy.len()!=bound.len() {
        return Err(StripEnthalpyError::Invalid("thermal interval changed bridge correspondence"));
    }
    let mut result=BridgeThermalExchange {layers:BTreeMap::new(),
        mechanical_internal_energy_change_j:0.0,mechanical_entropy_energy_change_j:0.0,
        reversible_heat_j:0.0,dissipated_j:0.0};
    for (before,after) in old.cells().iter().zip(next.cells()) {
        checkpoint()?;
        if before.id()!=after.id() || before.length_m()!=after.length_m()
            || before.section().layers().len()!=after.section().layers().len() {
            return Err(StripEnthalpyError::Invalid("thermal interval changed reference section geometry"));
        }
        for (&layer_id,layer) in before.section().layers() {
            checkpoint()?;
            let id=(before.id(),layer_id);
            let binding=bound.get(&id).ok_or(StripEnthalpyError::Invalid("unbound thermal layer"))?;
            let following=after.section().layers().get(&layer_id)
                .ok_or(StripEnthalpyError::Invalid("thermal interval removed material history"))?;
            let cell=binding.section.cell;
            let initial=old_heat.cells().get(&cell)
                .ok_or(StripEnthalpyError::Invalid("missing original caloric cell"))?;
            let endpoint=next_heat.cells().get(&cell)
                .ok_or(StripEnthalpyError::Invalid("missing proposed caloric cell"))?;
            let account=*energy.get(&id)
                .ok_or(StripEnthalpyError::Invalid("missing original point energy account"))?;
            let exchange=layer_exchange(layer,following,initial,endpoint,&binding.section,
                before.length_m(),account)?;
            result.mechanical_internal_energy_change_j=finite(result.mechanical_internal_energy_change_j
                +exchange.mechanical_internal_energy_change_j)?;
            result.mechanical_entropy_energy_change_j=finite(result.mechanical_entropy_energy_change_j
                +exchange.mechanical_entropy_energy_change_j)?;
            result.reversible_heat_j=finite(result.reversible_heat_j+exchange.reversible_heat_j)?;
            result.dissipated_j=finite(result.dissipated_j+exchange.dissipated_j)?;
            result.layers.insert(id,exchange);
        }
    }
    if result.layers.len()!=bound.len() {
        return Err(StripEnthalpyError::Invalid("incomplete spatial thermal exchange"));
    }
    Ok(result)
}

pub struct BridgeEnthalpyReport<S> {
    pub mechanical: solid::Report,
    /// Original extensive free-energy/work/loss accounts [J].
    pub layer_energy_j: BTreeMap<(u64,u64), Energy>,
    /// Includes dashpot heat and signed reversible exchange; do not add them again.
    pub thermal: EnergyAudit,
    pub source: S,
    pub nonlinear: NewtonReport,
    pub structural_trials: usize,
    pub maximum_cell_residual_j: f64,
    /// Signed total material heat, physical dashpot loss plus reversible exchange.
    pub deposited_heat_j: f64,
    pub reversible_heat_j: f64,
    pub mechanical_internal_energy_change_j: f64,
    pub mechanical_entropy_energy_change_j: f64,
    pub non_mechanical_heat_j: f64,
    /// Mechanical internal + caloric storage + numerical loss - load work - other heat.
    /// Gravity is included in the original load work, not counted again as heat.
    pub total_energy_residual_j: f64,
}
pub struct BridgeEnthalpyCandidate<'t,'m,S> {
    thermal: DepositionSimulation<'t>,
    bridge: solid::Bridge<'m>,
    report: BridgeEnthalpyReport<S>,
}
impl<'t,'m,S> BridgeEnthalpyCandidate<'t,'m,S> {
    pub fn thermal(&self) -> &DepositionSimulation<'t> { &self.thermal }
    pub fn bridge(&self) -> &solid::Bridge<'m> { &self.bridge }
    pub fn report(&self) -> &BridgeEnthalpyReport<S> { &self.report }
    pub fn into_parts(self) -> (DepositionSimulation<'t>,solid::Bridge<'m>,BridgeEnthalpyReport<S>) {
        (self.thermal,self.bridge,self.report)
    }
}

/// One immutable, spatially coupled interval. Every born section layer needs a
/// distinct active thermal cell. Extra passive thermal cells may exchange heat.
/// Energy guesses use sorted thermal CellId order, not axial or layer order.
/// The resolver supplies independent heat and conductances, never Maxwell heat.
/// Split at events/load jumps; endpoint events remain pending on success.
#[allow(clippy::too_many_arguments)]
pub fn step_candidate<'t,'m,S,E>(
    thermal: &DepositionSimulation<'t>, bridge: &solid::Bridge<'m>,
    bindings: &[BridgeThermalBinding<'_>], target_s: f64, gravity_m_s2: f64,
    additional_downward_n_m: &[f64], config: BridgeEnthalpyConfig,
    initial_energy_increments_j: &[f64],
    resolve: impl Fn(&DepositionSimulation<'t>) -> Result<ThermalCoefficients<S>,E>,
    cancelled: impl FnMut() -> bool,
) -> Result<BridgeEnthalpyCandidate<'t,'m,S>,BridgeEnthalpyError<E>> {
    let cancel=RefCell::new(cancelled);
    let checkpoint=|| if (*cancel.borrow_mut())() {Err(DepositionError::Cancelled)} else {Ok(())};
    checkpoint()?;
    let c=config.feedback;
    c.newton.validate()?;
    if thermal.time_s()!=bridge.time_s() || bindings.is_empty() || c.maximum_cells==0
        || thermal.cells().len()>c.maximum_cells || bridge.cells().len()>c.maximum_cells
        || bindings.len()>thermal.cells().len() || initial_energy_increments_j.len()!=thermal.cells().len()
        || c.maximum_structural_trials==0 || c.maximum_newton_iterations==0
        || !config.maximum_deflection_m.is_finite() || config.maximum_deflection_m<=0.0
        || [c.cell_energy_tolerance_j,c.thermal_energy_tolerance_j,c.total_energy_tolerance_j]
            .iter().any(|v|!v.is_finite() || *v<0.0) {
        return Err(invalid("complete bounded common-time states and finite energy controls required"));
    }
    let bound=admit_bindings::<E>(thermal,bridge,bindings,checkpoint)?;
    let step=BackwardEulerStep::new(thermal,target_s,checkpoint)?;
    let dt=step.duration_s(); let trials=Cell::new(0usize);
    let mechanics=solid::Config {section:c.mechanics,maximum_deflection_m:config.maximum_deflection_m};
    let evaluate=|x:&[f64]| step.evaluate(x,|endpoint| {
        checkpoint()?;
        if trials.get()>=c.maximum_structural_trials {
            return Err(BridgeEnthalpyError::Coupled(StripEnthalpyError::StructuralBudget));
        }
        trials.set(trials.get()+1);
        let mut temperatures=BTreeMap::new();
        for (&id,b) in &bound { checkpoint()?; temperatures.insert(id,endpoint.temperature_k(b.section.cell)?); }
        let (next,mechanical,layer_energy)=bridge.advance_with_cell_energy(target_s,&temperatures,
            gravity_m_s2,additional_downward_n_m,mechanics,&mut ||(*cancel.borrow_mut())())?;
        let exchange=bridge_exchange(thermal,endpoint,bridge,&next,&bound,&layer_energy,checkpoint)
            .map_err(BridgeEnthalpyError::Coupled)?;
        let coefficients=resolve(endpoint)
            .map_err(|e|BridgeEnthalpyError::Coupled(StripEnthalpyError::Source(e)))?;
        let other=finite(dt*endpoint.external_step_power_w(&coefficients.sources_w,
            &coefficients.conductances_w_k,checkpoint)?)?;
        let mut sources=coefficients.sources_w; let mut deposited=0.0;
        for (id,account) in &exchange.layers {
            checkpoint()?;
            deposited=finite(deposited+deposit_heat_j(&mut sources,bound[id].section.cell,
                account.total_material_heat_j,dt,c.total_energy_tolerance_j)
                .map_err(BridgeEnthalpyError::Coupled)?)?;
        }
        if (exchange.dissipated_j-mechanical.energy_j.dissipated).abs()>c.total_energy_tolerance_j
            || (deposited-exchange.dissipated_j-exchange.reversible_heat_j).abs()>c.total_energy_tolerance_j {
            return Err(invalid("spatial dashpot/reversible heat differs from the whole bridge account"));
        }
        Ok::<_,BridgeEnthalpyError<E>>(StepCoefficients {sources_w:sources,
            conductances_w_k:coefficients.conductances_w_k,
            payload:(next,mechanical,layer_energy,coefficients.source,deposited,other,exchange)})
    },checkpoint);
    let numerical=|e|match e {FiniteDifferenceError::Source(e)=>e,
        FiniteDifferenceError::Invalid(s)=>invalid(s)};
    let problem=FallibleFiniteDifference::new(step.cell_ids().len(),c.difference_step_j,
        |x:&[f64]|evaluate(x).map(|trial|trial.residual_j().to_vec())).map_err(numerical)?;
    let construction=NewtonKrylovState::new(&problem,initial_energy_increments_j.to_vec(),c.newton);
    if let Some(e)=problem.take_failure(){return Err(numerical(e));}
    let mut solver=construction?;
    let result=solver.run_controlled(&problem,c.maximum_newton_iterations,
        &mut ||(*cancel.borrow_mut())() || problem.has_failed());
    if let Some(e)=problem.take_failure(){return Err(numerical(e));}
    checkpoint()?;
    if result.progress!=SolverRunProgress::Complete {
        return Err(BridgeEnthalpyError::Coupled(StripEnthalpyError::Interrupted(result.progress)));
    }
    if !result.report.converged {
        return Err(BridgeEnthalpyError::Coupled(StripEnthalpyError::NotConverged(result.report)));
    }
    let trial=evaluate(&solver.x)?;
    let maximum_cell_residual_j=trial.maximum_cell_residual_j();
    let (thermal,(next,mechanical,layer_energy_j,source,deposited_heat_j,non_mechanical_heat_j,exchange),audit)=
        trial.accept(c.cell_energy_tolerance_j,c.thermal_energy_tolerance_j,checkpoint)?;
    let total=finite(exchange.mechanical_internal_energy_change_j
        +audit.stored_energy_change_j+mechanical.energy_j.numerical_dissipation
        -mechanical.load_work_j-non_mechanical_heat_j)?;
    if total.abs()>c.total_energy_tolerance_j
        || finite(audit.external_heat_j-deposited_heat_j-non_mechanical_heat_j)?.abs()>c.total_energy_tolerance_j {
        return Err(invalid("independent bridge work and local/global thermal energy do not close"));
    }
    checkpoint()?;
    Ok(BridgeEnthalpyCandidate {thermal,bridge:next,report:BridgeEnthalpyReport {
        mechanical,layer_energy_j,thermal:audit,source,nonlinear:result.report,
        structural_trials:trials.get(),maximum_cell_residual_j,deposited_heat_j,
        reversible_heat_j:exchange.reversible_heat_j,
        mechanical_internal_energy_change_j:exchange.mechanical_internal_energy_change_j,
        mechanical_entropy_energy_change_j:exchange.mechanical_entropy_energy_change_j,
        non_mechanical_heat_j,total_energy_residual_j:total}})
}

/// The same source, mass and history checks serve intervals and spatial births.
fn admit_bindings<'b,'m,E>(
    thermal:&DepositionSimulation<'_>, bridge:&solid::Bridge<'_>,
    bindings:&'b [BridgeThermalBinding<'m>],
    mut checkpoint:impl FnMut()->Result<(),DepositionError>,
) -> Result<BTreeMap<(u64,u64),&'b BridgeThermalBinding<'m>>,BridgeEnthalpyError<E>> {
    checkpoint()?;
    if thermal.time_s()!=bridge.time_s() || bindings.is_empty()
        || bindings.len()>thermal.cells().len() {
        return Err(invalid("complete born bindings at one physical time required"));
    }
    let mut axial=BTreeMap::new(); let mut expected=0usize;
    for cell in bridge.cells() {
        checkpoint()?;
        expected=expected.checked_add(cell.section().layers().len())
            .ok_or_else(||invalid("layer count overflow"))?;
        axial.insert(cell.id(),cell);
    }
    if expected!=bindings.len() {return Err(invalid("one binding per born axial-cell/layer required"));}
    let mut bound=BTreeMap::new(); let mut destinations=BTreeSet::new(); let mut masses=BTreeMap::new();
    for binding in bindings {
        checkpoint()?;
        let b=binding.section;
        let cell=axial.get(&binding.axial_cell).ok_or_else(||invalid("unknown axial cell"))?;
        let layer=cell.section().layers().get(&b.layer).ok_or_else(||invalid("unknown section layer"))?;
        let heat=thermal.cells().get(&b.cell).ok_or_else(||invalid("inactive thermal destination"))?;
        if bound.insert((binding.axial_cell,b.layer),binding).is_some() || !destinations.insert(b.cell)
            || !thermal.uses_material_law(heat.material_id,b.caloric.curve()) {
            return Err(invalid("duplicate destination or incompatible caloric object"));
        }
        let density=b.material.reference_density_kg_m3(b.caloric)
            .map_err(|e|BridgeEnthalpyError::Coupled(StripEnthalpyError::Binding(e)))?;
        let volume=finite(layer.geometry.area_m2*cell.length_m())?;
        if !same(volume,heat.geometry.volume_m3) || !same(finite(volume*density)?,heat.mass_kg) {
            return Err(invalid("layer thermal mass/volume do not match the material reference"));
        }
        let mass=masses.entry(binding.axial_cell).or_insert(0.0);
        *mass=finite(*mass+heat.mass_kg)?;
        for point in layer.material_points() {
            if !std::ptr::eq(point.law(),b.material.law()) {
                return Err(invalid("layer does not retain its declared mechanical source"));
            }
            CaloricPoint::try_new(point.clone(),b.caloric.curve(),heat.specific_h_j_kg,
                heat.mass_kg,volume,heat.geometry.surface_area_m2)
                .map_err(|e|BridgeEnthalpyError::Coupled(StripEnthalpyError::CaloricAdmission(e)))?;
        }
    }
    for cell in bridge.cells() {
        checkpoint()?;
        if !same(cell.mass_kg(),*masses.get(&cell.id()).ok_or_else(||invalid("unbound bridge mass"))?) {
            return Err(invalid("gravity and heat storage must use the same reference mass"));
        }
    }
    checkpoint()?;
    Ok(bound)
}

/// Same spatial solve with existing source-pinned material conductivities.
/// No electrical branches are ignored and no new material selection occurs in
/// Newton trials. The original network admits every thermal cell/link, including
/// passive bodies; its full source snapshots remain in the accepted report.
#[allow(clippy::too_many_arguments)]
pub fn step_transport_candidate<'t,'m>(
    thermal: &DepositionSimulation<'t>, bridge: &solid::Bridge<'m>,
    bindings: &[BridgeThermalBinding<'_>], transport: &PreparedMaterialNetwork<'_, '_>,
    target_s: f64, gravity_m_s2: f64, additional_downward_n_m: &[f64],
    config: BridgeEnthalpyConfig, guess_j: &[f64], cancelled: impl FnMut() -> bool,
) -> Result<BridgeEnthalpyCandidate<'t,'m,MaterialNetworkSample>,BridgeEnthalpyError<TransportError>> {
    let cancel=RefCell::new(cancelled);
    step_candidate(thermal,bridge,bindings,target_s,gravity_m_s2,additional_downward_n_m,config,guess_j,
        |state| {
            let resolved=transport.resolve(state,||(*cancel.borrow_mut())()).map_err(TransportError::Network)?;
            if !resolved.resistors.is_empty() {
                return Err(TransportError::Invalid("electrical branches require their circuit energy owner"));
            }
            for binding in bindings {
                if (*cancel.borrow_mut())() {
                    return Err(TransportError::Network(DepositionError::Cancelled.into()));
                }
                let b=binding.section;
                let index=resolved.source.cells.binary_search_by_key(&b.cell,|s|s.cell)
                    .map_err(|_|TransportError::Invalid("bridge layer lacks a transport cell"))?;
                if let Some(material)=&resolved.source.cells[index].conductivity {
                    b.material.admit_state_point(material)
                        .map_err(|error|TransportError::Context {layer:b.layer,error})?;
                }
            }
            Ok(ThermalCoefficients {sources_w:resolved.sources_w,
                conductances_w_k:resolved.conductances_w_k,source:resolved.source})
        },||(*cancel.borrow_mut())())
}
