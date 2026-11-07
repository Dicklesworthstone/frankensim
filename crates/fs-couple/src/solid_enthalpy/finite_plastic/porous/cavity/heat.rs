//! Heat exchange between source gas and retained finite thermal porous walls.
//!
//! Each interval first balances gas heat at the ACTUAL accepted cavity volume,
//! then uses that physical gas as a new isentropic reference for native finite
//! wall motion. The heated gas temperature is held during the original explicit
//! wall heat substeps. This is a first-order heat-before-motion split; refine
//! the complete physical interval for time accuracy.
//!
//! No gas is mapped back to undeformed volume. Original material calorics,
//! plastic compaction, optional inertia and current-area pressure remain owned
//! by their native solvers. Every nonlinear probe starts from the same accepted
//! gas, wall, thermal clock and independently cloned observer history.
//!
//! Supplied total conductances [W/K] are explicit fixed links. Optional source
//! contacts reuse original material conduction/convection at held patch geometry.
//! Neither path infers moving area or pore-fluid flow.
//! Native solid transport and all original observer hooks remain active.
//! Original ideal or homogeneous Peng--Robinson parcels supply rigid heating,
//! direct caloric work, phase-domain admission and true acoustic stiffness.

mod gas;
mod observer;
pub mod source;
pub mod march;

pub use gas::ThermalCavityGas;
pub use source::{CavityWallSourceContact, CavityWallSourceReceipt};

use super::{
    FiniteThermalPorousCavitiesEnergy, FiniteThermalPorousCavityEquilibrium,
    FiniteThermalPorousSolid,
};
use crate::iqn_ils::driver::CouplingMethod;
use crate::solid_enthalpy::finite_plastic::{
    FiniteThermoplasticConfig, FiniteThermoplasticError, FiniteThermoplasticReport,
};
use fs_exec::Cx;
use fs_material::thermochemical_gas::energy::{GasEnergyError, GasEnergySolve};
use fs_material::thermochemical_gas::energy::entropy::GasEntropySolve;
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_math::dd::Dd;
use fs_solid::hyper3::{HyperTetError, HyperTetSettings};
use fs_solid::hyper3::plastic::porous::thermal::FiniteThermoplasticControls;
use fs_solid::hyper3::pressure::{ClosedTetPressureSurface, TetPressure};
use fs_solid::hyper3::pressure::cavity::{
    AdiabaticCavity, CavityGasEpoch, ClosedTetCavity,
};
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::{
    NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport, SolverRunProgress,
};
use fs_thermal::enthalpy::{
    CellId, DepositionError, DepositionObserver, DepositionSimulation, StepConfig,
};
use std::cell::Cell;
use std::collections::BTreeSet;
use std::convert::Infallible;
use std::sync::Arc;

use observer::WallObserver;

/// One original ideal or Peng--Robinson inventory and its exact prepared chamber.
#[derive(Clone, Debug)]
pub struct ThermalPorousGasCavity<G = ClosedGasParcel> {
    pub cavity: Arc<ClosedTetCavity>,
    pub gas: G,
    pub exterior_pressure_pa: f64,
    pub entropy: GasEntropySolve,
}

/// Explicit fixed total conductance between one cavity and an original wall cell.
/// Distinct links may share a cell; all their conductances enter its timestep bound.
#[derive(Clone, Copy, Debug)]
pub struct CavityWallHeatLink {
    pub id: u64,
    pub cavity: usize,
    pub cell: CellId,
    pub conductance_w_k: f64,
}

/// Bounded nonlinear heat closure; native mechanical/thermal budgets are separate.
#[derive(Clone, Copy, Debug)]
pub struct CavityHeatSolve {
    pub gas: GasEnergySolve,
    pub maximum_trials: usize,
    pub maximum_newton_iterations: usize,
    /// Positive finite difference increment in the gas heat coordinates [J].
    pub difference_step_j: f64,
    pub newton: NewtonKrylovConfig,
    /// Independently bounds each gas, the interface and complete first law [J].
    pub energy_tolerance_j: f64,
}

/// Original observer or native transport refusal while evaluating a heat link.
#[derive(Debug)]
pub enum CavityHeatObserverError<E> {
    Invalid(&'static str),
    Thermal(DepositionError),
    Source(crate::solid_enthalpy::gas::SolidGasError),
    Observer(E),
}
impl<E> From<DepositionError> for CavityHeatObserverError<E> {
    fn from(error: DepositionError) -> Self { Self::Thermal(error) }
}

/// A refused proposal leaves every accepted owner unchanged.
#[derive(Debug)]
pub enum CavityHeatError<E = Infallible> {
    Invalid(&'static str),
    Binding(FiniteThermoplasticError),
    Geometry(HyperTetError),
    Gas(GasEnergyError),
    Source(crate::solid_enthalpy::gas::SolidGasError),
    Feedback(FiniteThermoplasticError<CavityHeatObserverError<E>>),
    Newton(NewtonError),
    Interrupted(SolverRunProgress),
    NotConverged(NewtonReport),
    TrialBudget,
    Cancelled,
}
impl<E: core::fmt::Debug> core::fmt::Display for CavityHeatError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "thermal porous cavity heat refused: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for CavityHeatError<E> {}

fn poll<E>(cx: &Cx<'_>) -> Result<(), CavityHeatError<E>> {
    cx.checkpoint().map_err(|_| CavityHeatError::Cancelled)
}
fn number<E>(value: f64) -> Result<f64, CavityHeatError<E>> {
    if value.is_finite() { Ok(value) }
    else { Err(CavityHeatError::Invalid("unrepresentable cavity heat or energy account")) }
}
fn difference_error<E>(
    error: FiniteDifferenceError<CavityHeatError<E>>,
) -> CavityHeatError<E> {
    match error {
        FiniteDifferenceError::Source(error) => error,
        FiniteDifferenceError::Invalid(message) => CavityHeatError::Invalid(message),
    }
}

/// Per-link physical heat integral on the accepted native thermal substeps.
#[derive(Clone, Debug)]
pub struct CavityWallHeatReceipt {
    pub id: u64,
    pub cavity: usize,
    pub cell: CellId,
    /// Positive means heat left this gas and entered this wall cell.
    pub heat_into_solid_j: f64,
    pub last_cell_temperature_k: f64,
    pub last_power_into_solid_w: f64,
}

/// One gas's split heat and work, with no absolute formation-energy subtraction.
#[derive(Clone, Debug)]
pub struct CavityGasHeatBalance {
    pub rigid_heat_stage_energy_change_j: f64,
    pub adiabatic_energy_change_j: f64,
    pub total_energy_change_j: f64,
    pub external_heat_j: f64,
    pub heat_into_solid_j: f64,
    pub heat_residual_j: f64,
    /// Temperature held during native wall heat integration, before wall motion.
    pub heat_stage_temperature_k: f64,
}

/// Original native work/heat receipts and independently checked internal exchange.
pub struct CavityHeatBalance<G = ClosedGasParcel> {
    pub native: FiniteThermoplasticReport<
        FiniteThermalPorousCavityEquilibrium<CavityGasEpoch<G>>,
        FiniteThermalPorousCavitiesEnergy,
    >,
    pub gases: Vec<CavityGasHeatBalance>,
    pub links: Vec<CavityWallHeatReceipt>,
    /// Source-owned patches on the same accepted native thermal substeps.
    pub source_contacts: Vec<CavityWallSourceReceipt>,
    pub total_gas_energy_change_j: f64,
    /// ACTUAL integrated source added to the native wall cells [J].
    pub gas_to_solid_heat_j: f64,
    /// Native external heat after removing this internal gas-wall transfer [J].
    pub external_solid_heat_j: f64,
    /// Integrated sum of absolute cell source-representation discrepancies [J].
    /// Opposing cells or substeps cannot conceal lost interface heat.
    pub interface_heat_defect_j: f64,
    /// Native combined residual + rigid gas heat + wall film heat - external gas heat.
    /// Plastic heat, gas pressure work and numerical work defects are not duplicated.
    pub combined_residual_j: f64,
}

/// Retained physical gas inventories, compaction history and explicit heat links.
pub struct ThermalPorousGasSystem<'a, G = ClosedGasParcel> {
    solid: FiniteThermalPorousSolid<'a>,
    cavities: Vec<ThermalPorousGasCavity<G>>,
    links: Vec<CavityWallHeatLink>,
    source_contacts: Option<Arc<source::PreparedContacts<'a, G>>>,
}

/// Commit thermal, solid, gas and original observer histories together.
pub struct CavityHeatCandidate<'a, 't, O, G = ClosedGasParcel> {
    pub thermal: DepositionSimulation<'t>,
    pub system: ThermalPorousGasSystem<'a, G>,
    pub observer: O,
    pub balance: CavityHeatBalance<G>,
    pub nonlinear: NewtonReport,
    pub trials: usize,
}

struct Trial<'a, 't, O, G> {
    thermal: DepositionSimulation<'t>,
    solid: FiniteThermalPorousSolid<'a>,
    observer: O,
    gases: Vec<G>,
    balance: CavityHeatBalance<G>,
}

impl<'a, G: ThermalCavityGas> ThermalPorousGasSystem<'a, G> {
    /// Bind 1..=64 independent cavities and at most 16,384 explicit heat links.
    ///
    /// Gas volumes must equal actual accepted geometric volumes. Cavity faces
    /// cannot overlap, and every link names a bound original structural cell.
    /// Empty links explicitly mean adiabatic walls. Subsequent states retain the
    /// exact cavity preparations, sources, inventories and link identities.
    pub fn try_new(
        solid: FiniteThermalPorousSolid<'a>,
        thermal: &DepositionSimulation<'_>,
        cavities: Vec<ThermalPorousGasCavity<G>>,
        mut links: Vec<CavityWallHeatLink>,
        cx: &Cx<'_>,
    ) -> Result<Self, CavityHeatError> {
        poll(cx)?;
        if cavities.is_empty() || cavities.len() > 64 || links.len() > 16_384
            || thermal.time_s().to_bits() != solid.time_s().to_bits()
        {
            return Err(CavityHeatError::Invalid(
                "bounded nonempty cavities, bounded links and identical accepted clocks required",
            ));
        }
        solid.binding().enthalpies::<Infallible>(thermal, true, cx)
            .map_err(CavityHeatError::Binding)?;
        let mechanics = solid.mechanics();
        let displacement = mechanics.displacement_m().as_flattened();
        let mut faces = BTreeSet::new();
        for entry in &cavities {
            poll(cx)?;
            if !entry.cavity.surface().matches_reference(
                mechanics.reference_nodes_m(), mechanics.reference_tetrahedra(),
            ) {
                return Err(CavityHeatError::Invalid(
                    "cavity must use the exact original solid reference mesh",
                ));
            }
            for &face in entry.cavity.solid_face_indices() {
                if !faces.insert(face) {
                    return Err(CavityHeatError::Invalid("independent gas cavity faces overlap"));
                }
            }
            let epoch = CavityGasEpoch::try_new(
                Arc::clone(&entry.cavity), displacement, entry.gas.clone(), cx,
            ).map_err(CavityHeatError::Geometry)?;
            AdiabaticCavity {
                cavity: entry.cavity.as_ref(),
                reference: &epoch,
                exterior_pressure_pa: entry.exterior_pressure_pa,
                controls: entry.entropy,
            }.evaluate(displacement, cx).map_err(CavityHeatError::Geometry)?;
        }
        let cells: BTreeSet<_> = solid.cells().iter().map(|cell| cell.id).collect();
        links.sort_by_key(|link| link.id);
        for (i, link) in links.iter().enumerate() {
            poll(cx)?;
            if link.cavity >= cavities.len() || !cells.contains(&link.cell)
                || !link.conductance_w_k.is_finite() || link.conductance_w_k < 0.0
                || i > 0 && links[i - 1].id == link.id
            {
                return Err(CavityHeatError::Invalid(
                    "unique heat link IDs, original cells, valid gas indices and nonnegative conductances required",
                ));
            }
        }
        poll(cx)?;
        Ok(Self { solid, cavities, links, source_contacts: None })
    }

    pub fn solid(&self) -> &FiniteThermalPorousSolid<'a> { &self.solid }
    pub fn cavities(&self) -> &[ThermalPorousGasCavity<G>] { &self.cavities }
    pub fn links(&self) -> &[CavityWallHeatLink] { &self.links }
    pub fn time_s(&self) -> f64 { self.solid.time_s() }

    /// Advance physical heat, actual finite cavity pressure and retained compaction.
    ///
    /// Supply one external heat integral and one initial rigid gas heat guess per
    /// cavity, in caller cavity order. Every nonlinear/Jacobian probe is bounded
    /// and disposable. The original observer must clone independent mutable
    /// history; it must not perform irreversible writes during a proposal.
    ///
    /// All native feedback/source/pressure/inertia gates still apply. Passing
    /// these balances establishes neither time accuracy nor a calibrated film
    /// model. The caller supplies the conductances and refines the time interval.
    #[allow(clippy::too_many_arguments)]
    pub fn advance<'t, O: DepositionObserver>(
        &self,
        thermal: &DepositionSimulation<'t>,
        observer: &O,
        target_s: f64,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        pressure: Option<(&ClosedTetPressureSurface, TetPressure<'_>)>,
        settings: HyperTetSettings,
        mechanical_controls: FiniteThermoplasticControls,
        thermal_config: StepConfig,
        feedback: FiniteThermoplasticConfig,
        method: CouplingMethod,
        external_gas_heat_j: &[f64],
        initial_heat_j: &[f64],
        controls: CavityHeatSolve,
        cx: &Cx<'_>,
    ) -> Result<CavityHeatCandidate<'a, 't, O, G>, CavityHeatError<O::Error>> {
        poll(cx)?;
        let count = self.cavities.len();
        let tolerance = controls.energy_tolerance_j;
        if !target_s.is_finite() || target_s <= self.time_s()
            || thermal.time_s().to_bits() != self.time_s().to_bits()
            || external_gas_heat_j.len() != count || initial_heat_j.len() != count
            || external_gas_heat_j.iter().chain(initial_heat_j).any(|q| !q.is_finite())
            || !tolerance.is_finite() || tolerance <= 0.0
            || !feedback.energy_tolerance_j.is_finite()
            || feedback.energy_tolerance_j <= 0.0 || feedback.energy_tolerance_j > tolerance
            || controls.maximum_trials == 0 || controls.maximum_newton_iterations == 0
        {
            return Err(CavityHeatError::Invalid(
                "forward interval, one finite heat per gas, bounded solve and compatible energy gates required",
            ));
        }
        self.solid.binding().enthalpies::<Infallible>(thermal, true, cx)
            .map_err(CavityHeatError::Binding)?;
        for entry in &self.cavities {
            let gas_tolerance = number(entry.gas.mass_kg() * controls.gas.energy_tolerance_j_kg)?;
            if gas_tolerance <= 0.0 || gas_tolerance > tolerance {
                return Err(CavityHeatError::Invalid("gas caloric inversion exceeds enclosing energy gate"));
            }
        }
        // Admit immutable accepted geometry once. All caloric probes share its
        // compact anchor; the native mechanical solver still admits each NEW
        // trial displacement through the original complete-solid owner.
        let mut anchors = Vec::with_capacity(count);
        for entry in &self.cavities {
            anchors.push(CavityGasEpoch::try_new(
                Arc::clone(&entry.cavity),
                self.solid.mechanics().displacement_m().as_flattened(),
                entry.gas.clone(), cx,
            ).map_err(CavityHeatError::Geometry)?);
        }
        let trials = Cell::new(0usize);
        let evaluate = |heat: &[f64]| {
            poll(cx)?;
            if heat.len() != count || heat.iter().any(|q| !q.is_finite()) {
                return Err(CavityHeatError::Invalid("one finite rigid heat coordinate per cavity required"));
            }
            if trials.get() >= controls.maximum_trials { return Err(CavityHeatError::TrialBudget); }
            trials.set(trials.get() + 1);
            let mut epochs = Vec::with_capacity(count);
            let mut rigid_energy = Vec::with_capacity(count);
            let mut temperatures = Vec::with_capacity(count);
            for ((entry, anchor), &q) in self.cavities.iter().zip(&anchors).zip(heat) {
                poll(cx)?;
                let mut gas = entry.gas.clone();
                let heating = gas.add_heat(
                    q, controls.gas, || cx.checkpoint().is_err(),
                ).map_err(CavityHeatError::Gas)?;
                temperatures.push(gas.temperature_k());
                rigid_energy.push(heating.stored_energy_change_j);
                epochs.push(anchor.with_gas_at_anchor(gas, cx)
                    .map_err(CavityHeatError::Geometry)?);
            }
            let sources: Vec<_> = self.cavities.iter().zip(&epochs).map(|(entry, epoch)| {
                AdiabaticCavity {
                    cavity: entry.cavity.as_ref(), reference: epoch,
                    exterior_pressure_pa: entry.exterior_pressure_pa,
                    controls: entry.entropy,
                }
            }).collect();
            let source = self.source_contacts.as_ref()
                .map(|prepared| prepared.flux(&epochs,
                    &mut || cx.checkpoint().map_err(|_| DepositionError::Cancelled),
                )).transpose()
                .map_err(CavityHeatError::Source)?;
            let wall = WallObserver::new(&self.links, &temperatures, source, observer.clone());
            let next = self.solid.advance_with_cavities_and_coupling(
                thermal, &wall, target_s, prescribed_m, nodal_forces_n, body_force_n_m3,
                &sources, pressure, settings, mechanical_controls, thermal_config,
                feedback, method, cx,
            ).map_err(CavityHeatError::Feedback)?;
            let (thermal, solid, wall, native) = next.into_parts();
            let mut gases = Vec::with_capacity(count);
            let mut balances = Vec::with_capacity(count);
            let mut total_rigid = Dd::ZERO;
            let mut total_external = Dd::ZERO;
            let mut total_gas = Dd::ZERO;
            for i in 0..count {
                poll(cx)?;
                let gas = native.equilibrium.cavities[i].gas.gas();
                let adiabatic = native.material_energy.cavity_energy[i].gas_stored_energy_change_j;
                let total = number((Dd::from_f64(rigid_energy[i]) + Dd::from_f64(adiabatic)).to_f64())?;
                let film = number(wall.gas_heat_j[i].to_f64())?;
                let residual = number((Dd::from_f64(rigid_energy[i])
                    - Dd::from_f64(external_gas_heat_j[i]) + wall.gas_heat_j[i]).to_f64())?;
                balances.push(CavityGasHeatBalance {
                    rigid_heat_stage_energy_change_j: rigid_energy[i],
                    adiabatic_energy_change_j: adiabatic, total_energy_change_j: total,
                    external_heat_j: external_gas_heat_j[i], heat_into_solid_j: film,
                    heat_residual_j: residual, heat_stage_temperature_k: temperatures[i],
                });
                gases.push(gas.clone());
                total_rigid = total_rigid + Dd::from_f64(rigid_energy[i]);
                total_external = total_external + Dd::from_f64(external_gas_heat_j[i]);
                total_gas = total_gas + Dd::from_f64(total);
            }
            let film = number(wall.applied_heat_j.to_f64())?;
            let external_solid = number((Dd::from_f64(native.external_heat_j)
                - wall.applied_heat_j).to_f64())?;
            let combined = number((Dd::from_f64(native.combined_residual_j)
                + total_rigid + wall.applied_heat_j - total_external).to_f64())?;
            let balance = CavityHeatBalance {
                native, gases: balances, links: wall.receipts()?,
                source_contacts: wall.source_receipts()?,
                total_gas_energy_change_j: number(total_gas.to_f64())?,
                gas_to_solid_heat_j: film, external_solid_heat_j: external_solid,
                interface_heat_defect_j: number(wall.interface_defect_j.to_f64())?,
                combined_residual_j: combined,
            };
            poll(cx)?;
            Ok(Trial { thermal, solid, observer: wall.inner, gases, balance })
        };
        let problem = FallibleFiniteDifference::new(count, controls.difference_step_j, |heat: &[f64]| {
            let trial = evaluate(heat)?;
            Ok::<_, CavityHeatError<O::Error>>(trial.balance.gases.iter()
                .map(|gas| gas.heat_residual_j).collect())
        }).map_err(difference_error)?;
        let construction = NewtonKrylovState::new(&problem, initial_heat_j.to_vec(), controls.newton);
        if let Some(error) = problem.take_failure() { return Err(difference_error(error)); }
        let mut solver = construction.map_err(CavityHeatError::Newton)?;
        let result = solver.run_controlled(
            &problem, controls.maximum_newton_iterations,
            &mut || cx.checkpoint().is_err() || problem.has_failed(),
        );
        if let Some(error) = problem.take_failure() { return Err(difference_error(error)); }
        poll(cx)?;
        if result.progress != SolverRunProgress::Complete {
            return Err(CavityHeatError::Interrupted(result.progress));
        }
        if !result.report.converged { return Err(CavityHeatError::NotConverged(result.report)); }
        let trial = evaluate(&solver.x)?;
        if trial.balance.gases.iter().any(|gas| gas.heat_residual_j.abs() > tolerance)
            || trial.balance.interface_heat_defect_j > tolerance
            || trial.balance.combined_residual_j.abs() > tolerance
        {
            return Err(CavityHeatError::Invalid("gas, interface or complete energy balance did not close"));
        }
        let cavities = self.cavities.iter().zip(trial.gases).map(|(entry, gas)| {
            ThermalPorousGasCavity {
                cavity: Arc::clone(&entry.cavity), gas,
                exterior_pressure_pa: entry.exterior_pressure_pa, entropy: entry.entropy,
            }
        }).collect();
        poll(cx)?;
        Ok(CavityHeatCandidate {
            thermal: trial.thermal,
            system: Self {
                solid: trial.solid, cavities, links: self.links.clone(),
                source_contacts: self.source_contacts.clone(),
            },
            observer: trial.observer, balance: trial.balance,
            nonlinear: result.report, trials: trials.get(),
        })
    }
}
