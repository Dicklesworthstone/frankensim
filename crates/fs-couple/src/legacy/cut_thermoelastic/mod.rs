//! Reciprocal finite heating of the original Cartesian or adaptive CutFEM solid.
//!
//! The existing caloric timeline, material heat observer and native elasticity
//! solve remain the owners. Every Picard trial starts from accepted heat and
//! observer history. Local source-potential exchange supplies reversible heat;
//! ghost/weak-form work is a separate signed numerical account, never heating.
//!
//! This first path holds reference geometry, inventory, stiffness scales and
//! homogeneous strong supports fixed, with independent nodal loads. Embedded
//! supports, topology changes, finite deformation and phase loss of solid
//! support refuse. Endpoint work is trapezoidal and needs timestep refinement.
//! It is not an exact reversible-entropy or continuum-accuracy claim.

mod binding;
pub mod dynamic;
mod heat;
mod native;
mod observer;

pub use binding::CutThermoelasticCell3;
pub use heat::CutThermoelasticHeatError3;
pub use native::CutElasticityMut3;

use std::convert::Infallible;
use fs_cutfem::elastic3::thermal::{
    ThermoelasticCellExchange3, ThermoelasticEquilibrium3,
    ThermoelasticEquilibriumSettings3, ThermoelasticExchangeLimits3,
};
use fs_cutfem::elastic3::{ElasticMaterialField3, ElasticityError3, MaterialRefreshDecision3};
use fs_exec::Cx;
use fs_material::elastic::thermal::ThermoelasticError;
use fs_material::phase::PhaseStateError;
use fs_material::thermal::source::ThermalMaterialError;
use fs_math::dd::Dd;
use fs_solid::cutfront::dynamics3::CutElasticityRef3;
use fs_thermal::enthalpy::{
    CellId, DepositionError, DepositionObserver, DepositionSimulation, EnergyAudit, StepConfig,
};

use binding::Bindings;
use observer::ReversibleHeat;

pub const MAX_CUT_THERMOELASTIC_ITERATIONS: usize = 256;

/// Independent work, allocation and physical-account admission controls.
/// No unbounded default is supplied.
#[derive(Clone, Copy, Debug)]
pub struct CutThermoelasticConfig3 {
    pub maximum_cells: usize,
    pub maximum_iterations: usize,
    /// Relaxation of local reversible joules, in (0,1].
    pub heat_relaxation: f64,
    /// Sum of absolute cell applied/source heat mismatches [J].
    pub heat_tolerance_j: f64,
    /// Separate mechanical and combined heat/work residual allowances [J].
    pub energy_tolerance_j: f64,
    /// Absolute signed work of numerical stabilization over this interval [J].
    /// This is an approximation allowance, not permission to heat the body.
    pub maximum_numerical_work_j: f64,
    pub exchange: ThermoelasticExchangeLimits3,
}

#[derive(Debug)]
pub enum CutThermoelasticError3<E = Infallible> {
    Invalid(&'static str),
    Native(ElasticityError3),
    Source(ThermoelasticError),
    Initialization(fs_time::galpha::initialization::InitialSolveError),
    Dynamics(fs_time::galpha::TimeSolveError),
    Phase(PhaseStateError),
    State(ThermalMaterialError),
    Thermal(DepositionError),
    Observer(E),
    NotConverged {
        iterations: usize,
        heat_defect_j: f64,
        combined_residual_j: f64,
    },
    NumericalWork { work_j: f64, allowance_j: f64 },
    Cancelled,
}
impl<E> From<ElasticityError3> for CutThermoelasticError3<E> {
    fn from(error: ElasticityError3) -> Self {
        match error { ElasticityError3::Cancelled => Self::Cancelled, other => Self::Native(other) }
    }
}
impl<E> From<fs_time::galpha::initialization::InitialSolveError> for CutThermoelasticError3<E> {
    fn from(error: fs_time::galpha::initialization::InitialSolveError) -> Self {
        match error {
            fs_time::galpha::initialization::InitialSolveError::Cancelled => Self::Cancelled,
            other => Self::Initialization(other),
        }
    }
}
impl<E> From<fs_time::galpha::TimeSolveError> for CutThermoelasticError3<E> {
    fn from(error: fs_time::galpha::TimeSolveError) -> Self {
        match error {
            fs_time::galpha::TimeSolveError::Cancelled => Self::Cancelled,
            other => Self::Dynamics(other),
        }
    }
}
impl<E> From<DepositionError> for CutThermoelasticError3<E> {
    fn from(error: DepositionError) -> Self {
        if error == DepositionError::Cancelled { Self::Cancelled } else { Self::Thermal(error) }
    }
}
impl<E> From<ThermoelasticError> for CutThermoelasticError3<E> {
    fn from(error: ThermoelasticError) -> Self { Self::Source(error) }
}
impl<E> From<PhaseStateError> for CutThermoelasticError3<E> {
    fn from(error: PhaseStateError) -> Self { Self::Phase(error) }
}
impl<E> From<ThermalMaterialError> for CutThermoelasticError3<E> {
    fn from(error: ThermalMaterialError) -> Self { Self::State(error) }
}
impl<E: core::fmt::Debug> core::fmt::Display for CutThermoelasticError3<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native cut-solid/finite heat: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for CutThermoelasticError3<E> {}

fn poll<E>(cx: &Cx<'_>) -> Result<(), CutThermoelasticError3<E>> {
    cx.checkpoint().map_err(|_| CutThermoelasticError3::Cancelled)
}
fn checkpoint(cx: &Cx<'_>) -> std::ops::ControlFlow<()> {
    if cx.checkpoint().is_err() { std::ops::ControlFlow::Break(()) }
    else { std::ops::ControlFlow::Continue(()) }
}
fn finite<E>(value: f64) -> Result<f64, CutThermoelasticError3<E>> {
    if value.is_finite() { Ok(value) }
    else { Err(CutThermoelasticError3::Invalid("unrepresentable cut-solid heat or work")) }
}
fn controls<E>(config: CutThermoelasticConfig3) -> Result<(), CutThermoelasticError3<E>> {
    if config.maximum_cells == 0 || config.maximum_iterations == 0
        || config.maximum_iterations > MAX_CUT_THERMOELASTIC_ITERATIONS
        || !config.heat_relaxation.is_finite() || config.heat_relaxation <= 0.0
        || config.heat_relaxation > 1.0
        || [config.heat_tolerance_j, config.energy_tolerance_j]
            .iter().any(|x| !x.is_finite() || *x <= 0.0)
        || !config.maximum_numerical_work_j.is_finite()
        || config.maximum_numerical_work_j < 0.0
    {
        return Err(CutThermoelasticError3::Invalid("bounded work and finite positive coupling tolerances required"));
    }
    Ok(())
}
fn forces<E>(native: CutElasticityRef3<'_>, f: &[f64], maximum: usize)
    -> Result<(), CutThermoelasticError3<E>>
{
    if native::n(native) > maximum || f.len() != native::n(native)
        || f.iter().any(|x| !x.is_finite())
    {
        return Err(CutThermoelasticError3::Invalid("bounded complete independent nodal force vector required"));
    }
    Ok(())
}

struct History {
    time_s: f64,
    temperatures_k: Vec<f64>,
    specific_h_j_kg: Vec<f64>,
    nodal_forces_n: Vec<f64>,
    equilibrium: ThermoelasticEquilibrium3,
}

/// Accepted compatible solid history paired with the original caloric owner.
/// The native operator is borrowed for each operation and is never copied.
/// Retained native preparations reject a foreign or subsequently changed owner.
pub struct CutThermoelasticSolid3<'m> {
    bindings: Bindings<'m>,
    history: History,
    scales: Vec<f64>,
    fixed: Vec<bool>,
    settings: ThermoelasticEquilibriumSettings3,
}

/// A local physical source account in original active-cell order.
#[derive(Clone, Debug)]
pub struct CutThermoelasticCellStep3 {
    pub cell: CellId,
    pub temperature_k: [f64; 2],
    /// Original source exchange integrated at retained native bulk points.
    pub exchange: ThermoelasticCellExchange3,
}

/// Closed accounts from the interval just published to all four owners.
/// Numerical work is signed work spent by the discrete stabilization; it is
/// not a material energy, dissipation law or extra caloric source.
#[derive(Debug)]
pub struct CutThermoelasticReport3 {
    pub iterations: usize,
    /// Original native audit includes applied reversible heat exactly once.
    pub thermal_audit: EnergyAudit,
    pub cells: Vec<CutThermoelasticCellStep3>,
    pub external_work_j: f64,
    pub native_work_j: f64,
    pub physical_bulk_work_j: f64,
    pub numerical_work_j: f64,
    pub internal_energy_change_j: f64,
    pub reversible_heat_j: f64,
    pub applied_heat_j: f64,
    /// Authored external heat, excluding reversible heat and enclosure net.
    pub external_heat_j: f64,
    pub closed_radiation_net_j: f64,
    /// Native resistance work plus represented kinetic change (when inertial)
    /// minus actual external work. Numerical timestep defects are not heat.
    pub mechanical_residual_j: f64,
    pub thermal_residual_j: f64,
    pub combined_residual_j: f64,
    pub heat_defect_j: f64,
}

struct Accepted<'t, O> {
    thermal: DepositionSimulation<'t>,
    observer: O,
    history: History,
    report: CutThermoelasticReport3,
}
enum Trial<'t, O> {
    Accepted(Accepted<'t, O>),
    Retry { heat_j: Vec<f64>, heat_defect_j: f64, combined_residual_j: f64 },
}

impl<'m> CutThermoelasticSolid3<'m> {
    /// Bind an existing native thermal material owner and solve its initial
    /// loaded equilibrium without changing it or the caloric timeline.
    /// Cells follow the operator's documented active-cell order. Native cut
    /// volume, original density/mass, complete material source and caloric
    /// object must agree; there is no conservative-looking rescaling.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        native: CutElasticityRef3<'_>, thermal: &DepositionSimulation<'_>,
        cells: Vec<CutThermoelasticCell3<'m>>, nodal_forces_n: &[f64],
        settings: ThermoelasticEquilibriumSettings3, config: CutThermoelasticConfig3,
        cx: &Cx<'_>,
    ) -> Result<Self, CutThermoelasticError3> {
        poll(cx)?; controls(config)?;
        forces(native, nodal_forces_n, settings.maximum_native_dofs)?;
        let bindings = Bindings::new(native, thermal, cells, config.maximum_cells, cx)?;
        let (temperatures_k, specific_h_j_kg) = bindings.admit(thermal, cx)?;
        let equilibrium = native::solve(native, nodal_forces_n, settings, cx)?;
        let history = History {
            time_s: thermal.time_s(), temperatures_k, specific_h_j_kg,
            nodal_forces_n: nodal_forces_n.to_vec(), equilibrium,
        };
        let scales = native::scales(native).to_vec();
        let fixed = native::fixed(native).to_vec();
        poll(cx)?;
        Ok(Self { bindings, history, scales, fixed, settings })
    }

    pub fn time_s(&self) -> f64 { self.history.time_s }
    pub fn equilibrium(&self) -> &ThermoelasticEquilibrium3 { &self.history.equilibrium }
    pub fn temperatures_k(&self) -> &[f64] { &self.history.temperatures_k }
    pub fn nodal_forces_n(&self) -> &[f64] { &self.history.nodal_forces_n }

    /// Jointly advance the existing native material, compatible equilibrium,
    /// finite thermal timeline and caller's observer history.
    ///
    /// Every rejection leaves all four accepted owners unchanged. A successful
    /// return has already replaced all four; dropping the report cannot split
    /// their epochs. The observer must obey the original DepositionObserver
    /// Clone contract: independent candidate history, without shared mutable
    /// external side effects. Arbitrary interior-sharing callbacks cannot be
    /// rolled back by either this adapter or the native caloric integrator.
    #[allow(clippy::too_many_arguments)]
    pub fn advance<'t, O: DepositionObserver>(
        &mut self, mut native: CutElasticityMut3<'_>,
        thermal: &mut DepositionSimulation<'t>, observer: &mut O,
        target_s: f64, nodal_forces_n: &[f64], thermal_config: StepConfig,
        config: CutThermoelasticConfig3, cx: &Cx<'_>,
    ) -> Result<CutThermoelasticReport3, CutThermoelasticError3<O::Error>> {
        let accepted = self.attempt(&mut native, thermal, observer, target_s,
            nodal_forces_n, thermal_config, config, |_| 0.0, cx)?;
        // No callbacks, validation or allocation after the native commit.
        // Retain displaced values until EVERY owner is replaced: an observer's
        // custom Drop may panic, but must not leave half the accepted epoch old.
        let old_history = core::mem::replace(&mut self.history, accepted.history);
        let old_thermal = core::mem::replace(thermal, accepted.thermal);
        let old_observer = core::mem::replace(observer, accepted.observer);
        drop((old_history, old_thermal, old_observer));
        Ok(accepted.report)
    }

    #[allow(clippy::too_many_arguments)]
    fn attempt<'t, O: DepositionObserver>(
        &self, native: &mut CutElasticityMut3<'_>, thermal: &DepositionSimulation<'t>,
        observer: &O, target_s: f64, nodal_forces_n: &[f64],
        thermal_config: StepConfig, config: CutThermoelasticConfig3,
        radiation_net: impl Fn(&O) -> f64, cx: &Cx<'_>,
    ) -> Result<Accepted<'t, O>, CutThermoelasticError3<O::Error>> {
        poll(cx)?; controls(config)?;
        let view = native.as_ref();
        forces(view, nodal_forces_n, self.settings.maximum_native_dofs)?;
        if self.bindings.cells.len() > config.maximum_cells
            || native::scales(view) != self.scales || native::fixed(view) != self.fixed
        {
            return Err(CutThermoelasticError3::Invalid("native scales, supports and bounded original cells must remain fixed"));
        }
        native::validate(view, &self.history.equilibrium)?;
        self.bindings.admit_native(view, cx)?;
        let (old_t, old_h) = self.bindings.admit(thermal, cx)?;
        if thermal.time_s().to_bits() != self.history.time_s.to_bits()
            || old_t.iter().zip(&self.history.temperatures_k).any(|(a,b)| a.to_bits()!=b.to_bits())
            || old_h.iter().zip(&self.history.specific_h_j_kg).any(|(a,b)| a.to_bits()!=b.to_bits())
        {
            return Err(CutThermoelasticError3::Invalid("accepted native and caloric histories are no longer paired"));
        }
        let dt = target_s - thermal.time_s();
        if !target_s.is_finite() || !dt.is_finite() || dt <= 0.0
            || !thermal_config.energy_tolerance_j.is_finite()
            || thermal_config.energy_tolerance_j < 0.0
            || thermal_config.energy_tolerance_j > config.energy_tolerance_j
        {
            return Err(CutThermoelasticError3::Invalid("positive finite interval and compatible caloric energy tolerance required"));
        }
        let mut applied = vec![0.0; self.bindings.cells.len()];
        let mut last = [f64::INFINITY; 2];
        for iteration in 1..=config.maximum_iterations {
            poll(cx)?;
            let coupled = ReversibleHeat::new(&self.bindings, observer.clone(), &applied, dt, cx)?;
            let (thermal_next, coupled_next, audit) = thermal.advance_observed(
                &coupled, target_s, thermal_config,
                || cx.checkpoint().map_err(|_| DepositionError::Cancelled),
            )?;
            // Drop a disposable observer clone BEFORE native candidate staging.
            // No user-owned destructor may run between material commit and publication.
            drop(coupled);
            let observer_next = coupled_next.inner;
            let applied_actual = coupled_next.applied_j;
            let radiation = finite(radiation_net(&observer_next))?;
            let (temperatures, enthalpies) = self.bindings.admit(&thermal_next, cx)?;
            let materials = self.bindings.materials(&temperatures, cx)?;
            let field = |bounds| self.bindings.material_at(bounds, &materials);
            let field = ElasticMaterialField3::ThermoelasticByCell(&field);
            let decision = native.refresh(field, self.settings.refresh, |candidate| -> Result<_, CutThermoelasticError3<O::Error>> {
                let equilibrium = native::solve(candidate, nodal_forces_n, self.settings, cx)?;
                let exchanges = native::exchange(candidate, &self.bindings.laws,
                    &self.history.temperatures_k, self.history.equilibrium.solution.coefficients(),
                    equilibrium.solution.coefficients(), self.settings.maximum_kinematic_strain,
                    config.exchange, cx)?;
                let (report, proposed) = self.account(&equilibrium, exchanges, &temperatures,
                    nodal_forces_n, &applied_actual, audit, radiation, iteration, config, cx)?;
                let converged = report.heat_defect_j <= config.heat_tolerance_j
                    && report.thermal_residual_j.abs() <= config.energy_tolerance_j
                    && report.combined_residual_j.abs() <= config.energy_tolerance_j;
                if converged {
                    // An unconverged trial can have larger numerical work than
                    // the final coupled state; this is an acceptance gate.
                    if report.numerical_work_j.abs() > config.maximum_numerical_work_j {
                        return Err(CutThermoelasticError3::NumericalWork {
                            work_j: report.numerical_work_j,
                            allowance_j: config.maximum_numerical_work_j,
                        });
                    }
                    let history = History {
                        time_s: target_s, temperatures_k: temperatures,
                        specific_h_j_kg: enthalpies, nodal_forces_n: nodal_forces_n.to_vec(),
                        equilibrium,
                    };
                    Ok(MaterialRefreshDecision3::Accept(Trial::Accepted(Accepted {
                        thermal: thermal_next, observer: observer_next, history, report,
                    })))
                } else {
                    // Candidate preparations are stale after rollback and are not
                    // returned as accepted state. Only local heat scratch survives.
                    drop((equilibrium, observer_next, thermal_next));
                    Ok(MaterialRefreshDecision3::Reject(Trial::Retry {
                        heat_j: proposed, heat_defect_j: report.heat_defect_j,
                        combined_residual_j: report.combined_residual_j,
                    }))
                }
            }, || checkpoint(cx))?;
            match decision {
                MaterialRefreshDecision3::Accept(Trial::Accepted(accepted)) => return Ok(accepted),
                MaterialRefreshDecision3::Reject(Trial::Retry { heat_j, heat_defect_j, combined_residual_j }) => {
                    last = [heat_defect_j, combined_residual_j];
                    for (i, (used, target)) in applied.iter_mut().zip(heat_j).enumerate() {
                        if i % 64 == 0 { poll(cx)?; }
                        *used = finite(config.heat_relaxation.mul_add(target - *used, *used))?;
                    }
                }
                _ => unreachable!("native material decision preserves its supplied payload"),
            }
        }
        Err(CutThermoelasticError3::NotConverged {
            iterations: config.maximum_iterations,
            heat_defect_j: last[0], combined_residual_j: last[1],
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn account<E>(
        &self, next: &ThermoelasticEquilibrium3, exchanges: Vec<ThermoelasticCellExchange3>,
        temperatures: &[f64], nodal_forces: &[f64], applied: &[Dd],
        audit: EnergyAudit, radiation: f64, iteration: usize,
        config: CutThermoelasticConfig3, cx: &Cx<'_>,
    ) -> Result<(CutThermoelasticReport3, Vec<f64>), CutThermoelasticError3<E>> {
        let before = &self.history.equilibrium;
        let mut external = Dd::ZERO;
        let mut mechanical_residual = Dd::ZERO;
        for (i, (&u0, &u1)) in before.solution.coefficients().iter()
            .zip(next.solution.coefficients()).enumerate()
        {
            if i % 192 == 0 { poll(cx)?; }
            let du = Dd::from_f64(u1) - Dd::from_f64(u0);
            external = external + Dd::from_f64(0.5)
                * (Dd::from_f64(self.history.nodal_forces_n[i]) + Dd::from_f64(nodal_forces[i])) * du;
            // Native full residual includes original ghost terms and supported
            // rows. Fixed displacement increments vanish exactly.
            mechanical_residual = mechanical_residual + Dd::from_f64(0.5)
                * (Dd::from_f64(before.equilibrium.residual[i])
                    + Dd::from_f64(next.equilibrium.residual[i])) * du;
        }
        account_exchange(&self.bindings, &self.history.temperatures_k,
            temperatures, exchanges, applied, external, external + mechanical_residual,
            Dd::ZERO, audit, radiation, iteration, config, cx)
            .map(|(report, heat, _)| (report, heat))
    }
}

#[allow(clippy::too_many_arguments)]
fn account_exchange<E>(
    bindings: &Bindings<'_>, before_temperatures: &[f64], temperatures: &[f64],
    exchanges: Vec<ThermoelasticCellExchange3>, applied: &[Dd], external: Dd,
    native_work: Dd, kinetic_change: Dd, audit: EnergyAudit, radiation: f64,
    iteration: usize, config: CutThermoelasticConfig3, cx: &Cx<'_>,
) -> Result<(CutThermoelasticReport3, Vec<f64>, Dd), CutThermoelasticError3<E>> {
    if exchanges.len() != bindings.cells.len()
        || audit.incoming_mass_kg != 0.0 || audit.incoming_energy_j != 0.0
    {
        return Err(CutThermoelasticError3::Invalid("fixed original cell inventory required throughout the interval"));
    }
    let mut bulk = Dd::ZERO;
    let mut internal = Dd::ZERO;
    let mut heat = Dd::ZERO;
    let mut used = Dd::ZERO;
    let mut defect = Dd::ZERO;
    let mut cells = Vec::with_capacity(exchanges.len());
    let mut proposed = Vec::with_capacity(exchanges.len());
    for (i, exchange) in exchanges.into_iter().enumerate() {
        poll(cx)?;
        let q = finite(exchange.retained_reversible_heat_j().to_f64())?;
        bulk = bulk + exchange.retained_mechanical_work_j();
        internal = internal + exchange.retained_internal_energy_change_j();
        heat = heat + exchange.retained_reversible_heat_j();
        used = used + applied[i];
        defect = defect + Dd::from_f64((applied[i] - exchange.retained_reversible_heat_j()).to_f64().abs());
        proposed.push(q);
        cells.push(CutThermoelasticCellStep3 {
            cell: bindings.cells[i].cell,
            temperature_k: [before_temperatures[i], temperatures[i]],
            exchange,
        });
    }
    let numerical = native_work - bulk;
    let numerical_work_j = finite(numerical.to_f64())?;
    let mechanical_residual_j = finite((kinetic_change + native_work - external).to_f64())?;
    if mechanical_residual_j.abs() > config.energy_tolerance_j {
        return Err(CutThermoelasticError3::Invalid("native work and represented kinetic change fail the mechanical energy allowance"));
    }
    // Use the native compensated balance, not subtraction of rounded
    // absolute caloric totals. Internal radiation roundoff is not supply.
    let base = Dd::from_f64(audit.residual_j) + Dd::from_f64(radiation);
    let thermal_residual = base + used - heat;
    let combined = base + used + internal + kinetic_change - external + numerical;
    let report = CutThermoelasticReport3 {
        iterations: iteration, thermal_audit: audit, cells,
        external_work_j: finite(external.to_f64())?,
        native_work_j: finite(native_work.to_f64())?,
        physical_bulk_work_j: finite(bulk.to_f64())?, numerical_work_j,
        internal_energy_change_j: finite(internal.to_f64())?,
        reversible_heat_j: finite(heat.to_f64())?,
        applied_heat_j: finite(used.to_f64())?,
        external_heat_j: finite((Dd::from_f64(audit.external_heat_j) - used
            - Dd::from_f64(radiation)).to_f64())?,
        closed_radiation_net_j: radiation, mechanical_residual_j,
        thermal_residual_j: finite(thermal_residual.to_f64())?,
        combined_residual_j: finite(combined.to_f64())?,
        heat_defect_j: finite(defect.to_f64())?,
    };
    poll(cx)?;
    Ok((report, proposed, internal))
}
