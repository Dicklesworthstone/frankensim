//! Reciprocal caloric feedback for progressive finite ductile P1/Q1 sheets.
//!
//! Each original element retains one caloric cell and all native quadrature
//! histories. Actual cell temperatures feed the supported initiation law;
//! original point volumes integrate the native plastic/damage heat fractions.
//! Every Picard/IQN proposal restarts from the SAME accepted mechanics,
//! thermal state and observer. Only the complete converged candidate escapes.
//!
//! The mechanical coefficients remain immutable. This adds no temperature law
//! for stiffness, yield, viscosity, expansion or post-initiation work. Thermal
//! cells and conductances retain their original geometry; through-thickness
//! gradients, deformed Fourier transport, melting, contact and fragments are
//! separate models. Optional inertia remains the original in-plane membrane
//! inertia, with numerical kinetic loss excluded from material heating.

use super::*;
use fs_solid::hyper2d::plastic::ductile::dynamics::{
    FiniteDuctileSheetMotion, FiniteDuctileSheetMotionControls,
};
use fs_solid::hyper2d::plastic::ductile::{FiniteDuctileSheet, FiniteDuctileSheetEnergy};

/// Original sheet work plus the shared explicit plastic/damage destinations.
pub type FiniteDuctileSheetThermalEnergy = FiniteDuctileThermalEnergy<FiniteDuctileSheetEnergy>;
pub type FiniteDuctileSheetThermalReport =
    FiniteThermoplasticReport<IncrementRecord, FiniteDuctileSheetThermalEnergy>;

impl MechanicalEnergy for FiniteDuctileSheetThermalEnergy {
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

#[derive(Debug)]
enum NativeSheet<'g, 'm, L>
where
    L: DifferentiableDuctileInitiationLaw + Clone + Debug + PartialEq,
{
    Static(FiniteDuctileSheet<'g, 'm, L>),
    Moving(FiniteDuctileSheetMotion<'g, 'm, L>),
}
impl<'g, 'm, L> NativeSheet<'g, 'm, L>
where
    L: DifferentiableDuctileInitiationLaw + Clone + Debug + PartialEq,
{
    fn sheet(&self) -> &FiniteDuctileSheet<'g, 'm, L> {
        match self {
            Self::Static(sheet) => sheet,
            Self::Moving(motion) => motion.sheet(),
        }
    }
}

struct TrialState<'g, 'm, L>
where
    L: DifferentiableDuctileInitiationLaw + Clone + Debug + PartialEq,
{
    mechanics: NativeSheet<'g, 'm, L>,
    nonthermal: Vec<DuctileNonthermalEnergy>,
}

/// Accepted original sheet, caloric mapping and irreversible nonthermal inventory.
#[derive(Debug)]
pub struct FiniteDuctileThermalSheet<'g, 'm, 's, L>
where
    L: DifferentiableDuctileInitiationLaw + Clone + Debug + PartialEq,
{
    mechanics: NativeSheet<'g, 'm, L>,
    sources: Arc<[FiniteDuctileThermalCell<'s>]>,
    point_counts: Arc<[usize]>,
    cells: Vec<DepositedCell>,
    nonthermal: Vec<DuctileNonthermalEnergy>,
    time_s: f64,
}

impl<'g, 'm, 's, L> FiniteDuctileThermalSheet<'g, 'm, 's, L>
where
    L: DifferentiableDuctileInitiationLaw + DuctileMaterialBinding + Clone + Debug + PartialEq,
{
    /// Bind one exact original caloric cell per element at the accepted time.
    /// Elastic preload is allowed; prior plastic/damage loss requires its
    /// original caloric/nonthermal history and cannot be adopted here.
    pub fn try_new(
        sheet: FiniteDuctileSheet<'g, 'm, L>,
        thermal: &DepositionSimulation<'_>,
        sources: Vec<FiniteDuctileThermalCell<'s>>,
        cx: &Cx<'_>,
    ) -> Result<Self, FiniteThermoplasticError> {
        poll(cx)?;
        if sources.len() != sheet.materials().len()
            || thermal.time_s() != sheet.physical_time_s()
            || sheet.quadrature_states().iter().any(|state| {
                state.native_state().alpha != 0.0
                    || state.damage() != 0.0
                    || state.dissipated_plastic_j_m3() != 0.0
                    || state.dissipated_damage_j_m3() != 0.0
            })
        {
            return Err(FiniteThermoplasticError::Invalid(
                "one caloric cell per sheet element, matching clocks and no unaccounted prior losses required",
            ));
        }
        let mut used = BTreeSet::new();
        let mut cells = Vec::with_capacity(sources.len());
        let mut point_counts = Vec::with_capacity(sources.len());
        let mut covered = 0;
        for (element, source) in sources.iter().enumerate() {
            poll(cx)?;
            let cell =
                thermal
                    .cells()
                    .get(&source.cell)
                    .ok_or(FiniteThermoplasticError::Invalid(
                        "finite ductile sheet caloric cell is inactive",
                    ))?;
            let temperature = source
                .caloric
                .state_at_specific_enthalpy(cell.specific_h_j_kg)
                .map_err(FiniteThermoplasticError::Caloric)?
                .temperature_k();
            let law = sheet.materials()[element].initiation_law();
            let span = law.domain().temperature_k;
            let volume =
                sheet
                    .reference_volume_m3(element)
                    .ok_or(FiniteThermoplasticError::Invalid(
                        "original ductile sheet element is absent",
                    ))?;
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
                || !super::super::same_measure(cell.geometry.volume_m3, volume)
                || !super::super::same_measure(
                    cell.mass_kg,
                    volume * source.caloric.reference_density_kg_m3(),
                )
            {
                return Err(FiniteThermoplasticError::Invalid(
                    "unique sheet cells/charts, original inventory, heat partition and compatible initiation source/temperature required",
                ));
            }
            let mut count = 0;
            while sheet.quadrature_volume_m3(element, count).is_some() {
                count += 1;
            }
            if count == 0 {
                return Err(FiniteThermoplasticError::Invalid(
                    "empty ductile sheet quadrature",
                ));
            }
            covered += count;
            point_counts.push(count);
            cells.push(cell.clone());
        }
        if covered != sheet.quadrature_states().len() {
            return Err(FiniteThermoplasticError::Invalid(
                "incomplete original sheet quadrature mapping",
            ));
        }
        poll(cx)?;
        Ok(Self {
            mechanics: NativeSheet::Static(sheet),
            sources: sources.into(),
            point_counts: point_counts.into(),
            nonthermal: vec![DuctileNonthermalEnergy::default(); cells.len()],
            cells,
            time_s: thermal.time_s(),
        })
    }

    /// Bind the SAME source reference densities to the original consistent mass.
    /// Physical thickness and every material history stay with the native sheet.
    pub fn with_inertia(
        mut self,
        velocity_mps: &[[f64; 2]],
        maximum_assembly_bytes: usize,
        cx: &Cx<'_>,
    ) -> Result<Self, FiniteThermoplasticError> {
        poll(cx)?;
        let densities: Vec<_> = self
            .sources
            .iter()
            .map(|source| source.caloric.reference_density_kg_m3())
            .collect();
        self.mechanics = match self.mechanics {
            NativeSheet::Static(sheet) => NativeSheet::Moving(
                FiniteDuctileSheetMotion::try_new(
                    sheet,
                    &densities,
                    velocity_mps,
                    maximum_assembly_bytes,
                    cx,
                )
                .map_err(FiniteThermoplasticError::DuctileSheet)?,
            ),
            NativeSheet::Moving(_) => {
                return Err(FiniteThermoplasticError::Invalid(
                    "ductile sheet inertia is already bound to its original velocity and mass",
                ));
            }
        };
        Ok(self)
    }

    pub fn mechanics(&self) -> &FiniteDuctileSheet<'g, 'm, L> {
        self.mechanics.sheet()
    }
    pub fn velocity_mps(&self) -> Option<&[[f64; 2]]> {
        match &self.mechanics {
            NativeSheet::Static(_) => None,
            NativeSheet::Moving(motion) => Some(motion.velocity_mps()),
        }
    }
    pub fn reference_mass_kg(&self) -> Option<f64> {
        match &self.mechanics {
            NativeSheet::Static(_) => None,
            NativeSheet::Moving(motion) => Some(motion.reference_mass_kg()),
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

    /// One physical interval with absolute supports and held nodal loads.
    /// Each feedback proposal uses the accepted sheet, never a preceding
    /// Picard/IQN candidate. Only the returned complete endpoint is accepted.
    #[allow(clippy::too_many_arguments)]
    pub fn advance<'t, O: DepositionObserver>(
        &self,
        thermal: &DepositionSimulation<'t>,
        observer: &O,
        target_s: f64,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        settings: IncrementSettings,
        controls: FiniteDuctileSheetMotionControls,
        thermal_config: StepConfig,
        config: FiniteThermoplasticConfig,
        cx: &Cx<'_>,
    ) -> Result<
        FiniteDuctileThermalSheetCandidate<'g, 'm, 's, 't, L, O>,
        FiniteThermoplasticError<O::Error>,
    > {
        self.advance_with_coupling(
            thermal,
            observer,
            target_s,
            prescribed_m,
            nodal_forces_n,
            settings,
            controls,
            thermal_config,
            config,
            CouplingMethod::RelaxedPicard,
            cx,
        )
    }

    /// Reuse the original bounded feedback driver, thermal timeline and observers.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_coupling<'t, O: DepositionObserver>(
        &self,
        thermal: &DepositionSimulation<'t>,
        observer: &O,
        target_s: f64,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        settings: IncrementSettings,
        controls: FiniteDuctileSheetMotionControls,
        thermal_config: StepConfig,
        config: FiniteThermoplasticConfig,
        method: CouplingMethod,
        cx: &Cx<'_>,
    ) -> Result<
        FiniteDuctileThermalSheetCandidate<'g, 'm, 's, 't, L, O>,
        FiniteThermoplasticError<O::Error>,
    > {
        let next = feedback::advance_with_coupling(
            feedback::ThermalBinding {
                cells: &self.cells,
                materials: FinitePlasticLawField::PerElement(&self.sources),
            },
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
                    settings,
                    controls,
                    config.energy_tolerance_j,
                    cx,
                )
            },
        )?;
        if next.state.mechanics.sheet().physical_time_s() != target_s {
            return Err(FiniteThermoplasticError::Invalid(
                "mechanical and caloric sheet clocks differ",
            ));
        }
        Ok(FiniteDuctileThermalSheetCandidate {
            thermal: next.thermal,
            sheet: Self {
                mechanics: next.state.mechanics,
                sources: Arc::clone(&self.sources),
                point_counts: Arc::clone(&self.point_counts),
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
        settings: IncrementSettings,
        controls: FiniteDuctileSheetMotionControls,
        tolerance_j: f64,
        cx: &Cx<'_>,
    ) -> Result<
        feedback::MechanicalProposal<
            TrialState<'g, 'm, L>,
            IncrementRecord,
            FiniteDuctileSheetThermalEnergy,
        >,
        FiniteThermoplasticError<E>,
    > {
        poll(cx)?;
        if temperatures.len() != self.sources.len() {
            return Err(FiniteThermoplasticError::Invalid(
                "one endpoint temperature per sheet cell required",
            ));
        }
        let mut point_temperatures = Vec::with_capacity(self.mechanics().quadrature_states().len());
        for (&temperature, &count) in temperatures.iter().zip(self.point_counts.iter()) {
            poll(cx)?;
            point_temperatures.extend(std::iter::repeat_n(temperature, count));
        }
        let conditions = DuctileStepConditions {
            duration_s,
            temperatures_k: DuctileTemperatures::PerPoint(&point_temperatures),
        };
        let (mechanics, equilibrium, mechanical, accounts, inertia) = match &self.mechanics {
            NativeSheet::Static(sheet) => {
                let step = sheet
                    .advance(
                        prescribed_m,
                        nodal_forces_n,
                        conditions,
                        settings,
                        controls.material,
                        cx,
                    )
                    .map_err(FiniteThermoplasticError::DuctileSheet)?;
                (
                    NativeSheet::Static(step.sheet),
                    step.equilibrium,
                    step.energy,
                    step.point_energy,
                    None,
                )
            }
            NativeSheet::Moving(motion) => {
                let step = motion
                    .advance(
                        prescribed_m,
                        nodal_forces_n,
                        [0.0; 2],
                        conditions,
                        settings,
                        controls,
                        cx,
                    )
                    .map_err(FiniteThermoplasticError::DuctileSheet)?;
                (
                    NativeSheet::Moving(step.motion),
                    step.equilibrium,
                    step.material_energy,
                    step.point_energy,
                    Some(step.inertia),
                )
            }
        };
        if accounts.len() != point_temperatures.len() {
            return Err(FiniteThermoplasticError::Invalid(
                "complete original sheet point work receipts required",
            ));
        }
        let mut cells = Vec::with_capacity(self.cells.len());
        let mut heat_j = Vec::with_capacity(self.cells.len());
        let mut nonthermal = Vec::with_capacity(self.cells.len());
        let mut totals = [Sum::default(); 4];
        let mut offset = 0;
        for (element, &count) in self.point_counts.iter().enumerate() {
            poll(cx)?;
            let mut cell_sums = [Sum::default(); 4];
            for point in 0..count {
                let volume = self
                    .mechanics()
                    .quadrature_volume_m3(element, point)
                    .ok_or(FiniteThermoplasticError::Invalid(
                        "original sheet point mapping changed",
                    ))?;
                let cell = split(
                    &accounts[offset],
                    self.sources[element].heat_partition,
                    volume,
                )?;
                for (sum, value) in cell_sums.iter_mut().zip([
                    cell.plastic_heat_j,
                    cell.damage_heat_j,
                    cell.nonthermal_plastic_j,
                    cell.nonthermal_damage_j,
                ]) {
                    sum.add(value);
                }
                offset += 1;
            }
            let values = cell_sums.map(Sum::value);
            for (sum, value) in totals.iter_mut().zip(values) {
                sum.add(value);
            }
            let cell = DuctileCellHeat {
                plastic_heat_j: values[0],
                damage_heat_j: values[1],
                nonthermal_plastic_j: values[2],
                nonthermal_damage_j: values[3],
            };
            heat_j.push(finite(cell.plastic_heat_j + cell.damage_heat_j)?);
            nonthermal.push(DuctileNonthermalEnergy {
                plastic_j: accumulate(
                    self.nonthermal[element].plastic_j,
                    cell.nonthermal_plastic_j,
                )?,
                damage_j: accumulate(self.nonthermal[element].damage_j, cell.nonthermal_damage_j)?,
            });
            cells.push(cell);
        }
        let [
            plastic_heat_j,
            damage_heat_j,
            nonthermal_plastic_j,
            nonthermal_damage_j,
        ] = totals.map(Sum::value);
        for residual in [
            plastic_heat_j + nonthermal_plastic_j - mechanical.plastic_dissipation_j,
            damage_heat_j + nonthermal_damage_j - mechanical.damage_release_j,
        ] {
            if finite(residual)?.abs() > tolerance_j {
                return Err(FiniteThermoplasticError::Invalid(
                    "sheet heat/nonthermal partition misses original native dissipation",
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

/// Original thermal/observer endpoint and the complete coupled sheet checkpoint.
pub struct FiniteDuctileThermalSheetCandidate<'g, 'm, 's, 't, L, O>
where
    L: DifferentiableDuctileInitiationLaw + Clone + Debug + PartialEq,
{
    thermal: DepositionSimulation<'t>,
    sheet: FiniteDuctileThermalSheet<'g, 'm, 's, L>,
    observer: O,
    report: FiniteDuctileSheetThermalReport,
}
impl<'g, 'm, 's, 't, L, O> FiniteDuctileThermalSheetCandidate<'g, 'm, 's, 't, L, O>
where
    L: DifferentiableDuctileInitiationLaw + Clone + Debug + PartialEq,
{
    pub const fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub const fn sheet(&self) -> &FiniteDuctileThermalSheet<'g, 'm, 's, L> {
        &self.sheet
    }
    pub const fn observer(&self) -> &O {
        &self.observer
    }
    pub const fn report(&self) -> &FiniteDuctileSheetThermalReport {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        FiniteDuctileThermalSheet<'g, 'm, 's, L>,
        O,
        FiniteDuctileSheetThermalReport,
    ) {
        (self.thermal, self.sheet, self.observer, self.report)
    }
}
