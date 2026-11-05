//! Reciprocal caloric feedback for the original finite-plastic sheet.
//!
//! Each P1/Q1 element is bound to one original caloric cell. Native quadrature
//! volumes integrate plastic loss, creep loss and signed reversible heat into
//! that cell; hardening storage and endpoint work defects are never heat.
//! The existing thermal timeline supplies transport, boundary sources, events
//! and substeps. Its actual endpoint temperatures drive the SAME accepted
//! finite-strain plane-stress histories on every feedback proposal.
//!
//! Materials are the original source-admitted families. The caller owns them
//! and constructs the native sheet first; this owner borrows those families
//! through the sheet instead of copying a caloric law or inventing a mixture.
//! Reference area times original thickness determines each fixed cell mass.
//!
//! This is an element-uniform thermal network and quasistatic membrane model.
//! Conductances retain their independently supplied reference interpretation;
//! they are not recomputed from the evolving thickness or deformed geometry.
//! Bending, through-thickness gradients, contact and inertial motion require
//! their own mechanical/thermal owners. Failed proposals change no accepted
//! sheet, thermal simulation or observer history.

use std::collections::BTreeSet;

use fs_exec::Cx;
use fs_solid::hyper2d::plastic::thermal::{
    FiniteThermoplasticElementEnergy, FiniteThermoplasticSheet,
};
use fs_solid::hyper3::plastic::{FinitePlasticLawField, FinitePlasticStepControls};
use fs_solid::plastic3::{IncrementRecord, IncrementSettings};
use fs_thermal::enthalpy::{
    CellId, DepositedCell, DepositionObserver, DepositionSimulation, StepConfig, ThermalLaw,
};

use super::feedback::{self, MechanicalProposal, ThermalBinding};
use super::{
    FiniteThermoplasticConfig, FiniteThermoplasticError, FiniteThermoplasticMaterial,
    FiniteThermoplasticReport, finite, poll,
};

/// Complete accepted sheet with one original caloric cell per element.
#[derive(Debug)]
pub struct FiniteThermoplasticSheetCoupling<'m, 's> {
    sheet: FiniteThermoplasticSheet<'m, FiniteThermoplasticMaterial<'s>>,
    cells: Vec<DepositedCell>,
    time_s: f64,
}

impl<'m, 's> FiniteThermoplasticSheetCoupling<'m, 's> {
    /// Bind an existing native sheet to the SAME accepted caloric state.
    /// Cell IDs follow original element order and must be unique and active.
    /// All material-point temperatures in an element must equal its cell's
    /// temperature; an element mean cannot replace a nonuniform accepted field.
    pub fn try_new(
        thermal: &DepositionSimulation<'_>,
        sheet: FiniteThermoplasticSheet<'m, FiniteThermoplasticMaterial<'s>>,
        cell_ids: &[CellId],
        cx: &Cx<'_>,
    ) -> Result<Self, FiniteThermoplasticError> {
        poll(cx)?;
        if cell_ids.len() != sheet.temperatures_k().len()
            || thermal.time_s() != sheet.physical_time_s()
        {
            return Err(FiniteThermoplasticError::Invalid(
                "one caloric cell per sheet element and matching accepted clocks required",
            ));
        }
        let mut used = BTreeSet::new();
        let mut cells = Vec::with_capacity(cell_ids.len());
        let families = sheet.families();
        let mut point_offset = 0;
        for (element, &id) in cell_ids.iter().enumerate() {
            poll(cx)?;
            let cell = thermal
                .cells()
                .get(&id)
                .ok_or(FiniteThermoplasticError::Invalid(
                    "finite sheet caloric cell is inactive",
                ))?;
            let material = match families {
                FinitePlasticLawField::Uniform(material) => material,
                FinitePlasticLawField::PerElement(materials) => &materials[element],
            };
            if !used.insert(id)
                || !thermal.uses_material_law(
                    cell.material_id,
                    material.caloric_source().curve() as &dyn ThermalLaw,
                )
            {
                return Err(FiniteThermoplasticError::Invalid(
                    "unique sheet cells and exact original caloric objects required",
                ));
            }
            let volume = sheet
                .reference_volume_m3(element)
                .expect("retained sheet element index");
            if !super::super::same_measure(cell.geometry.volume_m3, volume)
                || !super::super::same_measure(
                    cell.mass_kg,
                    volume * material.reference_density_kg_m3(),
                )
            {
                return Err(FiniteThermoplasticError::Invalid(
                    "sheet cell inventory must match original area, thickness and source density",
                ));
            }
            let temperature = material
                .state_at_specific_enthalpy(cell.specific_h_j_kg)
                .map_err(FiniteThermoplasticError::Caloric)?
                .temperature_k();
            let mut point = 0;
            while sheet.quadrature_volume_m3(element, point).is_some() {
                if sheet.quadrature_temperatures_k().get(point_offset) != Some(&temperature) {
                    return Err(FiniteThermoplasticError::Invalid(
                        "each accepted sheet point must match its original caloric cell temperature",
                    ));
                }
                point += 1;
                point_offset += 1;
            }
            cells.push(cell.clone());
        }
        if point_offset != sheet.quadrature_temperatures_k().len() {
            return Err(FiniteThermoplasticError::Invalid(
                "complete native sheet quadrature mapping required",
            ));
        }
        poll(cx)?;
        Ok(Self {
            sheet,
            cells,
            time_s: thermal.time_s(),
        })
    }

    pub const fn mechanics(
        &self,
    ) -> &FiniteThermoplasticSheet<'m, FiniteThermoplasticMaterial<'s>> {
        &self.sheet
    }

    pub const fn time_s(&self) -> f64 {
        self.time_s
    }

    /// Original cell identities and accepted caloric inventory, in element order.
    pub fn cells(&self) -> &[DepositedCell] {
        &self.cells
    }

    /// Advance mechanical and caloric state over the same actual interval.
    /// Every feedback proposal starts from this accepted sheet and the supplied
    /// accepted thermal simulation. Only a complete converged candidate returns.
    /// The observer must preserve its native independent-clone history contract.
    #[allow(clippy::too_many_arguments)]
    pub fn advance<'t, O: DepositionObserver>(
        &self,
        thermal: &DepositionSimulation<'t>,
        observer: &O,
        target_s: f64,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        settings: IncrementSettings,
        mechanical_controls: FinitePlasticStepControls,
        thermal_config: StepConfig,
        config: FiniteThermoplasticConfig,
        cx: &Cx<'_>,
    ) -> Result<FiniteThermoplasticSheetCandidate<'m, 's, 't, O>, FiniteThermoplasticError<O::Error>>
    {
        let binding = ThermalBinding {
            cells: &self.cells,
            materials: self.sheet.families(),
        };
        let candidate = feedback::advance(
            binding,
            self.time_s,
            thermal,
            observer,
            target_s,
            thermal_config,
            config,
            cx,
            |temperatures, duration_s| {
                let step = self
                    .sheet
                    .advance(
                        temperatures,
                        prescribed_m,
                        nodal_forces_n,
                        duration_s,
                        settings,
                        mechanical_controls,
                        cx,
                    )
                    .map_err(FiniteThermoplasticError::Sheet)?;
                let heat_j = self.caloric_heat_j(&step.quadrature_energy, cx)?;
                Ok(MechanicalProposal {
                    state: step.sheet,
                    equilibrium: step.equilibrium,
                    energy: step.energy,
                    inertia: None,
                    heat_j,
                })
            },
        )?;
        Ok(FiniteThermoplasticSheetCandidate {
            thermal: candidate.thermal,
            sheet: Self {
                sheet: candidate.state,
                cells: candidate.cells,
                time_s: target_s,
            },
            observer: candidate.observer,
            report: candidate.report,
        })
    }

    fn caloric_heat_j<E>(
        &self,
        accounts: &[FiniteThermoplasticElementEnergy],
        cx: &Cx<'_>,
    ) -> Result<Vec<f64>, FiniteThermoplasticError<E>> {
        let mut result = Vec::with_capacity(self.cells.len());
        let mut offset = 0;
        for element in 0..self.cells.len() {
            poll(cx)?;
            let mut heat = 0.0;
            let mut point = 0;
            while let Some(volume) = self.sheet.quadrature_volume_m3(element, point) {
                let account = accounts
                    .get(offset)
                    .ok_or(FiniteThermoplasticError::Invalid(
                        "one caloric account per native sheet quadrature point required",
                    ))?;
                heat = finite(heat + volume * account.caloric_heat_j_m3())?;
                offset += 1;
                point += 1;
            }
            result.push(heat);
        }
        if offset != accounts.len() {
            return Err(FiniteThermoplasticError::Invalid(
                "extra native sheet caloric quadrature account",
            ));
        }
        Ok(result)
    }
}

/// Shared physical balance checks with the original sheet equilibrium report.
pub type FiniteThermoplasticSheetReport = FiniteThermoplasticReport<IncrementRecord>;

/// A complete accepted thermal, sheet and observer endpoint.
pub struct FiniteThermoplasticSheetCandidate<'m, 's, 't, O> {
    thermal: DepositionSimulation<'t>,
    sheet: FiniteThermoplasticSheetCoupling<'m, 's>,
    observer: O,
    report: FiniteThermoplasticSheetReport,
}

impl<'m, 's, 't, O> FiniteThermoplasticSheetCandidate<'m, 's, 't, O> {
    pub const fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }

    pub const fn sheet(&self) -> &FiniteThermoplasticSheetCoupling<'m, 's> {
        &self.sheet
    }

    pub const fn observer(&self) -> &O {
        &self.observer
    }

    pub const fn report(&self) -> &FiniteThermoplasticSheetReport {
        &self.report
    }

    /// Adopt all accepted owners together after inspecting their physical report.
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        FiniteThermoplasticSheetCoupling<'m, 's>,
        O,
        FiniteThermoplasticSheetReport,
    ) {
        (self.thermal, self.sheet, self.observer, self.report)
    }
}
