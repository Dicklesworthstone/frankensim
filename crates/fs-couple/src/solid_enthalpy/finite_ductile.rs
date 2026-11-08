//! Reciprocal caloric feedback through the original finite ductile trajectory.
//!
//! The original nonpotential finite solid/momentum solve supplies one complete
//! endpoint and compact element energy receipts. Its actual physical duration
//! and caloric endpoint temperatures drive initiation; rejected coupling trials
//! never advance plastic/damage history, velocity, thermal cells or observers.
//! The shared finite thermoplastic feedback owner supplies the original heat
//! timeline, transport, events, bounded Picard/IQN-ILS iteration and work gates.
//!
//! Heat fractions are explicit fractions of ALREADY-ACCOUNTED plastic
//! dissipation and damage release, not Taylor--Quinney fractions of total plastic
//! work. Native energetic hardening is already retained in mechanical storage.
//! Unheated fractions remain in separate irreversible nonthermal inventories;
//! they are not silently discarded, identified as measured fracture-surface
//! energy, or later released without another physical model. Signed material
//! work defects and numerical kinetic losses NEVER enter either heat fraction.
//!
//! The mechanical backbone remains the supplied immutable fixed-coefficient
//! finite J2 law. Temperature changes initiation, without inventing thermal
//! expansion, coefficient transport or reversible entropy. Caloric storage uses
//! its explicit stress-free zero-pressure reference. Thermal cells/links retain
//! reference geometry and independently supplied conductances. This is not a
//! deformed Fourier operator, phase-change solid, crack/contact or fragment model.

use std::collections::BTreeSet;
use std::fmt::Debug;
use std::sync::Arc;

use crate::iqn_ils::driver::CouplingMethod;
use fs_exec::Cx;
use fs_material::ductile::evolution::DifferentiableDuctileInitiationLaw;
use fs_material::phase::ResolvedHeatCapacityCurve;
use fs_material::plastic::finite::ductile::FiniteDuctileEnergy;
use fs_solid::ductile::{DuctileMaterialBinding, DuctileStepConditions, DuctileTemperatures};
use fs_solid::hyper3::plastic::FinitePlasticLawField;
use fs_solid::hyper3::plastic::ductile::dynamics::{
    FiniteDuctileMotionControls, FiniteDuctileTetMotion,
};
use fs_solid::hyper3::plastic::ductile::{FiniteDuctileTetEnergy, FiniteDuctileTetSolid};
use fs_solid::plastic3::{IncrementRecord, IncrementSettings};
use fs_thermal::enthalpy::{
    CellId, DepositedCell, DepositionObserver, DepositionSimulation, StepConfig, ThermalLaw,
};

use super::finite_plastic::FiniteThermoplasticReport;
use super::finite_plastic::feedback::{self, CaloricMaterial, MechanicalEnergy};
pub use super::finite_plastic::{
    CaloricReference, FiniteThermoplasticConfig, FiniteThermoplasticError,
    ReferenceCaloricMaterial, ReferenceCaloricQuery,
};

/// Original finite ductile sheet histories coupled to their caloric cells.
pub mod sheet;

/// Explicit physical destinations for the native irreversible losses.
/// Each fraction must lie in [0,1]. The complementary fraction is retained in
/// the corresponding nonthermal inventory, with no implied later heat release.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DuctileHeatPartition {
    pub plastic_dissipation_heat_fraction: f64,
    pub damage_release_heat_fraction: f64,
}

impl DuctileHeatPartition {
    fn valid(self) -> bool {
        (0.0..=1.0).contains(&self.plastic_dissipation_heat_fraction)
            && (0.0..=1.0).contains(&self.damage_release_heat_fraction)
    }
}

/// One exact reference caloric source and explicitly authored heat partition.
#[derive(Clone, Debug)]
pub struct FiniteDuctileThermalCell<'a> {
    pub cell: CellId,
    pub caloric: ReferenceCaloricMaterial<'a>,
    pub heat_partition: DuctileHeatPartition,
}

impl CaloricMaterial for FiniteDuctileThermalCell<'_> {
    fn caloric_source(&self) -> &ResolvedHeatCapacityCurve {
        self.caloric.caloric_source()
    }
}

/// Irreversibly retained noncaloric energy in one original cell, in joules.
/// This has no recoverable constitutive force, entropy or release law here.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct DuctileNonthermalEnergy {
    pub plastic_j: f64,
    pub damage_j: f64,
}

/// Physical per-interval partition, in original cell order and joules.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DuctileCellHeat {
    pub plastic_heat_j: f64,
    pub damage_heat_j: f64,
    pub nonthermal_plastic_j: f64,
    pub nonthermal_damage_j: f64,
}

/// The unchanged native mechanical account and explicit physical destinations.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteDuctileThermalEnergy<M = FiniteDuctileTetEnergy> {
    pub mechanical: M,
    pub cells: Vec<DuctileCellHeat>,
    pub plastic_heat_j: f64,
    pub damage_heat_j: f64,
    pub nonthermal_plastic_j: f64,
    pub nonthermal_damage_j: f64,
}

impl MechanicalEnergy for FiniteDuctileThermalEnergy {
    fn internal_energy_change_j(&self) -> f64 {
        self.mechanical.stored_after_j - self.mechanical.stored_before_j
    }
    fn endpoint_work_defect_j(&self) -> f64 {
        self.mechanical.endpoint_work_defect_j
    }
    fn dead_load_endpoint_work_j(&self) -> f64 {
        self.mechanical.dead_load_endpoint_work_j
    }
    fn support_endpoint_work_j(&self) -> f64 {
        self.mechanical.support_endpoint_work_j
    }
    fn nonthermal_loss_j(&self) -> f64 {
        self.nonthermal_plastic_j + self.nonthermal_damage_j
    }
}

pub type FiniteDuctileThermalReport =
    FiniteThermoplasticReport<IncrementRecord, FiniteDuctileThermalEnergy>;

#[derive(Debug)]
enum NativeMechanics<'g, 'm, L>
where
    L: DifferentiableDuctileInitiationLaw + Clone + Debug + PartialEq,
{
    Static(FiniteDuctileTetSolid<'g, 'm, L>),
    Moving(FiniteDuctileTetMotion<'g, 'm, L>),
}

impl<'g, 'm, L> NativeMechanics<'g, 'm, L>
where
    L: DifferentiableDuctileInitiationLaw + Clone + Debug + PartialEq,
{
    fn solid(&self) -> &FiniteDuctileTetSolid<'g, 'm, L> {
        match self {
            Self::Static(solid) => solid,
            Self::Moving(motion) => motion.solid(),
        }
    }
}

struct TrialState<'g, 'm, L>
where
    L: DifferentiableDuctileInitiationLaw + Clone + Debug + PartialEq,
{
    mechanics: NativeMechanics<'g, 'm, L>,
    nonthermal: Vec<DuctileNonthermalEnergy>,
}

/// One accepted original specimen, caloric mapping and nonthermal inventory.
#[derive(Debug)]
pub struct FiniteDuctileThermalSolid<'g, 'm, 's, L>
where
    L: DifferentiableDuctileInitiationLaw + Clone + Debug + PartialEq,
{
    mechanics: NativeMechanics<'g, 'm, L>,
    sources: Arc<[FiniteDuctileThermalCell<'s>]>,
    cells: Vec<DepositedCell>,
    nonthermal: Vec<DuctileNonthermalEnergy>,
    time_s: f64,
}

impl<'g, 'm, 's, L> FiniteDuctileThermalSolid<'g, 'm, 's, L>
where
    L: DifferentiableDuctileInitiationLaw + DuctileMaterialBinding + Clone + Debug + PartialEq,
{
    /// Bind original cells before irreversible mechanical evolution. Elastic
    /// preload is admissible, but previous plastic/damage losses cannot be
    /// adopted without their original caloric/nonthermal history. Continue from
    /// the complete returned coupled owner after the first accepted interval.
    pub fn try_new(
        solid: FiniteDuctileTetSolid<'g, 'm, L>,
        thermal: &DepositionSimulation<'_>,
        sources: Vec<FiniteDuctileThermalCell<'s>>,
        cx: &Cx<'_>,
    ) -> Result<Self, FiniteThermoplasticError> {
        poll(cx)?;
        if sources.len() != solid.element_states().len()
            || solid.element_states().iter().any(|state| {
                state.native_state().alpha != 0.0
                    || state.damage() != 0.0
                    || state.dissipated_plastic_j_m3() != 0.0
                    || state.dissipated_damage_j_m3() != 0.0
            })
        {
            return Err(FiniteThermoplasticError::Invalid(
                "one caloric source per original point and no unaccounted prior plastic/damage losses required",
            ));
        }
        let mut used = BTreeSet::new();
        let mut cells = Vec::with_capacity(sources.len());
        for (i, source) in sources.iter().enumerate() {
            poll(cx)?;
            let cell =
                thermal
                    .cells()
                    .get(&source.cell)
                    .ok_or(FiniteThermoplasticError::Invalid(
                        "finite ductile caloric cell is inactive",
                    ))?;
            let temperature = source
                .caloric
                .state_at_specific_enthalpy(cell.specific_h_j_kg)
                .map_err(FiniteThermoplasticError::Caloric)?
                .temperature_k();
            let law = solid.materials()[i].initiation_law();
            let span = law.domain().temperature_k;
            let volume = solid.bands()[i].reference_volume_m3();
            if !used.insert(source.cell)
                || !source.heat_partition.valid()
                || !thermal.uses_material_law(
                    cell.material_id,
                    source.caloric.caloric_source().curve() as &dyn ThermalLaw,
                )
                || temperature != source.caloric.reference_temperature_k()
                || temperature < span[0]
                || temperature > span[1]
                || !law.matches_material(source.caloric.reference_state())
                || !super::same_measure(cell.geometry.volume_m3, volume)
                || !super::same_measure(
                    cell.mass_kg,
                    volume * source.caloric.reference_density_kg_m3(),
                )
            {
                return Err(FiniteThermoplasticError::Invalid(
                    "unique original cell/chart, reference inventory, explicit heat partition and compatible initiation source/temperature required",
                ));
            }
            cells.push(cell.clone());
        }
        poll(cx)?;
        Ok(Self {
            mechanics: NativeMechanics::Static(solid),
            nonthermal: vec![DuctileNonthermalEnergy::default(); cells.len()],
            sources: sources.into(),
            cells,
            time_s: thermal.time_s(),
        })
    }

    /// Bind the original sourced caloric reference density to native consistent
    /// inertia. Temperature, deformation and damage never replace that inventory.
    pub fn with_inertia(
        mut self,
        velocity_mps: &[[f64; 3]],
        cx: &Cx<'_>,
    ) -> Result<Self, FiniteThermoplasticError> {
        poll(cx)?;
        let densities: Vec<_> = self
            .sources
            .iter()
            .map(|source| source.caloric.reference_density_kg_m3())
            .collect();
        self.mechanics = match self.mechanics {
            NativeMechanics::Static(solid) => NativeMechanics::Moving(
                FiniteDuctileTetMotion::try_new(solid, &densities, velocity_mps, cx)
                    .map_err(FiniteThermoplasticError::DuctileSolid)?,
            ),
            NativeMechanics::Moving(_) => {
                return Err(FiniteThermoplasticError::Invalid(
                    "finite ductile inertia is already bound to its original velocity and mass",
                ));
            }
        };
        Ok(self)
    }

    pub fn mechanics(&self) -> &FiniteDuctileTetSolid<'g, 'm, L> {
        self.mechanics.solid()
    }
    pub fn velocity_mps(&self) -> Option<&[[f64; 3]]> {
        match &self.mechanics {
            NativeMechanics::Static(_) => None,
            NativeMechanics::Moving(motion) => Some(motion.velocity_mps()),
        }
    }
    pub fn reference_mass_kg(&self) -> Option<f64> {
        match &self.mechanics {
            NativeMechanics::Static(_) => None,
            NativeMechanics::Moving(motion) => Some(motion.reference_mass_kg()),
        }
    }
    pub fn sources(&self) -> &[FiniteDuctileThermalCell<'s>] {
        &self.sources
    }
    pub fn cells(&self) -> &[DepositedCell] {
        &self.cells
    }
    pub fn nonthermal_energy_j(&self) -> &[DuctileNonthermalEnergy] {
        &self.nonthermal
    }
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }

    fn binding(&self) -> feedback::ThermalBinding<'_, '_, FiniteDuctileThermalCell<'s>> {
        feedback::ThermalBinding {
            cells: &self.cells,
            materials: FinitePlasticLawField::PerElement(&self.sources),
        }
    }

    /// Propose one actual physical interval. Accept thermal, mechanical,
    /// nonthermal and observer histories together through the returned candidate.
    #[allow(clippy::too_many_arguments)]
    pub fn advance<'t, O: DepositionObserver>(
        &self,
        thermal: &DepositionSimulation<'t>,
        observer: &O,
        target_s: f64,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        settings: IncrementSettings,
        mechanical_controls: FiniteDuctileMotionControls,
        thermal_config: StepConfig,
        config: FiniteThermoplasticConfig,
        cx: &Cx<'_>,
    ) -> Result<
        FiniteDuctileThermalCandidate<'g, 'm, 's, 't, L, O>,
        FiniteThermoplasticError<O::Error>,
    > {
        self.advance_with_coupling(
            thermal,
            observer,
            target_s,
            prescribed_m,
            nodal_forces_n,
            body_force_n_m3,
            settings,
            mechanical_controls,
            thermal_config,
            config,
            CouplingMethod::RelaxedPicard,
            cx,
        )
    }

    /// Use the original bounded Picard/IQN-ILS enthalpy iteration. Its secants
    /// belong to this interval and never advance constitutive or observer time.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_coupling<'t, O: DepositionObserver>(
        &self,
        thermal: &DepositionSimulation<'t>,
        observer: &O,
        target_s: f64,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        settings: IncrementSettings,
        mechanical_controls: FiniteDuctileMotionControls,
        thermal_config: StepConfig,
        config: FiniteThermoplasticConfig,
        method: CouplingMethod,
        cx: &Cx<'_>,
    ) -> Result<
        FiniteDuctileThermalCandidate<'g, 'm, 's, 't, L, O>,
        FiniteThermoplasticError<O::Error>,
    > {
        let next = feedback::advance_with_coupling(
            self.binding(),
            self.time_s,
            thermal,
            observer,
            target_s,
            thermal_config,
            config,
            method,
            cx,
            |temperatures, duration_s| {
                self.propose(
                    temperatures,
                    duration_s,
                    prescribed_m,
                    nodal_forces_n,
                    body_force_n_m3,
                    settings,
                    mechanical_controls,
                    config.energy_tolerance_j,
                    cx,
                )
            },
        )?;
        Ok(FiniteDuctileThermalCandidate {
            thermal: next.thermal,
            solid: Self {
                mechanics: next.state.mechanics,
                sources: Arc::clone(&self.sources),
                cells: next.cells,
                nonthermal: next.state.nonthermal,
                time_s: target_s,
            },
            observer: next.observer,
            report: next.report,
        })
    }

    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    fn propose<E>(
        &self,
        temperatures: &[f64],
        duration_s: f64,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        settings: IncrementSettings,
        controls: FiniteDuctileMotionControls,
        tolerance_j: f64,
        cx: &Cx<'_>,
    ) -> Result<
        feedback::MechanicalProposal<
            TrialState<'g, 'm, L>,
            IncrementRecord,
            FiniteDuctileThermalEnergy,
        >,
        FiniteThermoplasticError<E>,
    > {
        poll(cx)?;
        let conditions = DuctileStepConditions {
            duration_s,
            temperatures_k: DuctileTemperatures::PerPoint(temperatures),
        };
        let (mechanics, equilibrium, mechanical, element_energy, inertia) = match &self.mechanics {
            NativeMechanics::Static(solid) => {
                let step = solid
                    .advance(
                        prescribed_m,
                        nodal_forces_n,
                        body_force_n_m3,
                        conditions,
                        settings,
                        controls.material,
                        cx,
                    )
                    .map_err(FiniteThermoplasticError::DuctileSolid)?;
                (
                    NativeMechanics::Static(step.solid),
                    step.equilibrium,
                    step.energy,
                    step.element_energy,
                    None,
                )
            }
            NativeMechanics::Moving(motion) => {
                let step = motion
                    .advance(
                        prescribed_m,
                        nodal_forces_n,
                        body_force_n_m3,
                        conditions,
                        settings,
                        controls,
                        cx,
                    )
                    .map_err(FiniteThermoplasticError::DuctileSolid)?;
                (
                    NativeMechanics::Moving(step.motion),
                    step.equilibrium,
                    step.material_energy,
                    step.element_energy,
                    Some(step.inertia),
                )
            }
        };
        if element_energy.len() != self.cells.len() {
            return Err(FiniteThermoplasticError::Invalid(
                "complete native finite ductile element work receipts required",
            ));
        }
        let mut cells = Vec::with_capacity(self.cells.len());
        let mut heat_j = Vec::with_capacity(self.cells.len());
        let mut nonthermal = Vec::with_capacity(self.cells.len());
        let mut totals = [Sum::default(); 4];
        for (i, native) in element_energy.iter().enumerate() {
            poll(cx)?;
            let partition = self.sources[i].heat_partition;
            let volume = self.mechanics().bands()[i].reference_volume_m3();
            let cell = split(native, partition, volume)?;
            for (sum, value) in totals.iter_mut().zip([
                cell.plastic_heat_j,
                cell.damage_heat_j,
                cell.nonthermal_plastic_j,
                cell.nonthermal_damage_j,
            ]) {
                sum.add(value);
            }
            heat_j.push(finite(cell.plastic_heat_j + cell.damage_heat_j)?);
            nonthermal.push(DuctileNonthermalEnergy {
                plastic_j: accumulate(self.nonthermal[i].plastic_j, cell.nonthermal_plastic_j)?,
                damage_j: accumulate(self.nonthermal[i].damage_j, cell.nonthermal_damage_j)?,
            });
            cells.push(cell);
        }
        let [
            plastic_heat_j,
            damage_heat_j,
            nonthermal_plastic_j,
            nonthermal_damage_j,
        ] = totals.map(|sum| sum.value());
        for residual in [
            plastic_heat_j + nonthermal_plastic_j - mechanical.plastic_dissipation_j,
            damage_heat_j + nonthermal_damage_j - mechanical.damage_release_j,
        ] {
            if finite(residual)?.abs() > tolerance_j {
                return Err(FiniteThermoplasticError::Invalid(
                    "element heat/nonthermal partition does not reproduce original native dissipation",
                ));
            }
        }
        poll(cx)?;
        Ok(feedback::MechanicalProposal {
            state: TrialState {
                mechanics,
                nonthermal,
            },
            equilibrium,
            energy: FiniteDuctileThermalEnergy {
                mechanical,
                cells,
                plastic_heat_j,
                damage_heat_j,
                nonthermal_plastic_j,
                nonthermal_damage_j,
            },
            inertia,
            heat_j,
        })
    }
}

/// The original thermal/observer endpoint and complete coupled structural state.
pub struct FiniteDuctileThermalCandidate<'g, 'm, 's, 't, L, O>
where
    L: DifferentiableDuctileInitiationLaw + Clone + Debug + PartialEq,
{
    thermal: DepositionSimulation<'t>,
    solid: FiniteDuctileThermalSolid<'g, 'm, 's, L>,
    observer: O,
    report: FiniteDuctileThermalReport,
}

impl<'g, 'm, 's, 't, L, O> FiniteDuctileThermalCandidate<'g, 'm, 's, 't, L, O>
where
    L: DifferentiableDuctileInitiationLaw + Clone + Debug + PartialEq,
{
    pub const fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub const fn solid(&self) -> &FiniteDuctileThermalSolid<'g, 'm, 's, L> {
        &self.solid
    }
    pub const fn observer(&self) -> &O {
        &self.observer
    }
    pub const fn report(&self) -> &FiniteDuctileThermalReport {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        FiniteDuctileThermalSolid<'g, 'm, 's, L>,
        O,
        FiniteDuctileThermalReport,
    ) {
        (self.thermal, self.solid, self.observer, self.report)
    }
}

fn split<E>(
    native: &FiniteDuctileEnergy,
    fraction: DuctileHeatPartition,
    volume_m3: f64,
) -> Result<DuctileCellHeat, FiniteThermoplasticError<E>> {
    let plastic = nonnegative(volume_m3 * native.plastic_dissipation_j_m3)?;
    let damage = nonnegative(volume_m3 * native.damage_release_j_m3)?;
    Ok(DuctileCellHeat {
        plastic_heat_j: nonnegative(fraction.plastic_dissipation_heat_fraction * plastic)?,
        damage_heat_j: nonnegative(fraction.damage_release_heat_fraction * damage)?,
        nonthermal_plastic_j: nonnegative(
            (1.0 - fraction.plastic_dissipation_heat_fraction) * plastic,
        )?,
        nonthermal_damage_j: nonnegative((1.0 - fraction.damage_release_heat_fraction) * damage)?,
    })
}

fn accumulate<E>(old: f64, increment: f64) -> Result<f64, FiniteThermoplasticError<E>> {
    let next = nonnegative(old + increment)?;
    if increment > 0.0 && next == old {
        Err(FiniteThermoplasticError::Invalid(
            "nonthermal increment is below retained-history resolution",
        ))
    } else {
        Ok(next)
    }
}
fn poll<E>(cx: &Cx<'_>) -> Result<(), FiniteThermoplasticError<E>> {
    cx.checkpoint()
        .map_err(|_| FiniteThermoplasticError::Cancelled)
}
fn finite<E>(value: f64) -> Result<f64, FiniteThermoplasticError<E>> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(FiniteThermoplasticError::Invalid(
            "finite ductile heat/work arithmetic is unrepresentable",
        ))
    }
}
fn nonnegative<E>(value: f64) -> Result<f64, FiniteThermoplasticError<E>> {
    if value.is_finite() && value >= 0.0 {
        Ok(value)
    } else {
        Err(FiniteThermoplasticError::Invalid(
            "finite nonnegative physical ductile energy required",
        ))
    }
}

#[derive(Clone, Copy, Default)]
struct Sum {
    total: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) {
        let next = self.total + value;
        self.correction += if self.total.abs() >= value.abs() {
            (self.total - next) + value
        } else {
            (value - next) + self.total
        };
        self.total = next;
    }
    fn value(self) -> f64 {
        self.total + self.correction
    }
}
