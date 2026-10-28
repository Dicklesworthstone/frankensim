//! Source-resolved mu(T), induction and caloric enthalpy on one accepted clock.
//! All trials restart from original A and H. The native magnetic endpoint-effort
//! step exposes parameter work P; this owner deposits Q_rev=−P−delta(T*S_field)
//! beside (and separately from) its actual spatial Ohmic heat. Consequently the
//! independent audit uses U_field=psi+T*S and original caloric energy exactly once.

use super::super::{
    CircuitEnthalpyReport, CircuitEnthalpySolve, DrivenCoefficients, cell_driver_candidate, finite,
    magnetic::{
        Error, MagneticThermalNetwork,
        caloric::{MagneticEnergyCell, PreparedMagneticCaloric},
        check, invalid, material, same_material,
    },
    material::MaterialNetwork,
};
use super::{EddyThermalSample, deposit_heat, induction};
use fs_em::planar::eddy::{EddyConductor, EddyControls, EddyStep, PlanarEddy, PlanarEddyCandidate};
use fs_feec::weighted::WeightedAssemblyLimits;
use fs_thermal::enthalpy::DepositionSimulation;
use std::{cell::RefCell, collections::BTreeMap};

pub use super::super::magnetic::caloric::{MagneticCaloricSample, MagneticThermalLaw};
#[derive(Clone, Debug, PartialEq)]
pub struct EddyCaloricSource {
    pub material: EddyThermalSample,
    pub fields: Vec<MagneticCaloricSample>,
    /// Includes temperature-independent fixed magnetic regions as psi=U.
    pub field_internal_energy_j: [f64; 2],
    pub reversible_heat_j: f64,
}
pub struct EddyCaloricReport {
    /// Native shared thermal residual plus native spatial induction report.
    /// The driver's work account removes its explicit parameter term; the
    /// independent physical-U audit below closes its thermal destination.
    pub coupled: CircuitEnthalpyReport<EddyCaloricSource, EddyStep>,
    pub external_heat_j: f64,
    pub total_energy_residual_j: f64,
}
pub struct EddyCaloricCandidate<'t, 'm> {
    thermal: DepositionSimulation<'t>,
    field: PlanarEddy<'m>,
    report: EddyCaloricReport,
}
impl<'t, 'm> EddyCaloricCandidate<'t, 'm> {
    pub fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub fn field(&self) -> &PlanarEddy<'m> {
        &self.field
    }
    pub fn report(&self) -> &EddyCaloricReport {
        &self.report
    }
    /// Both accepted histories must publish together; heat is already deposited.
    pub fn into_parts(self) -> (DepositionSimulation<'t>, PlanarEddy<'m>, EddyCaloricReport) {
        (self.thermal, self.field, self.report)
    }
}
impl MagneticThermalNetwork<'_> {
    /// Evolve local mu(T), sigma(T), A and caloric enthalpy together using the
    /// existing spatial field CG and existing enthalpy Newton residual. All
    /// changing magnetic regions (including insulators) need explicit derivative
    /// laws and whole-volume caloric cells. Fixed external regions need mu with
    /// no temperature dependence. Density, geometry and source choices stay fixed.
    ///
    /// The B envelope is enforced at BOTH accepted endpoints and the complete
    /// temperature interval has positive combined caloric+field capacity. This
    /// is first-order thermal coupling; time refinement remains necessary and
    /// no unconditional entropy/passivity or physical-validation claim is made.
    #[allow(clippy::too_many_arguments)]
    pub fn step_eddy_thermal_candidate<'t, 'm>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        field: &PlanarEddy<'m>,
        impressed_jz_a_m2: &[f64],
        bindings: &[MagneticThermalLaw],
        controls: CircuitEnthalpySolve,
        field_controls: EddyControls,
        assembly: WeightedAssemblyLimits,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<EddyCaloricCandidate<'t, 'm>, Error> {
        let cancel = RefCell::new(cancelled);
        check(&cancel)?;
        super::super::admit_driver_unknowns(thermal, controls, initial_energy_increments_j)?;
        if field.time_s().to_bits() != thermal.time_s().to_bits()
            || self.cells.len() != thermal.cells().len()
        {
            return Err(invalid(
                "matching field/thermal clocks and complete original cells required",
            ));
        }
        let network = MaterialNetwork {
            cells: self.cells,
            resistors: &[],
            links: self.links,
        };
        let prepared = network
            .prepare_with_bases(thermal, &[], || (*cancel.borrow_mut())())
            .map_err(material)?;
        let triangles: Vec<_> = field.conductors().iter().map(|s| s.triangle).collect();
        let routes = self.routes(thermal, field.problem(), &triangles, &cancel)?;
        let native = field
            .field(&mut || (*cancel.borrow_mut())())
            .map_err(induction)?;
        let originals: Vec<_> = field.problem().materials().iter().collect();
        let caloric = PreparedMagneticCaloric::new(
            self,
            thermal,
            field.problem().mesh(),
            &originals,
            native.magnetic_flux_density_t(),
            bindings,
            &routes,
            &cancel,
        )?;
        let conductors: BTreeMap<_, _> = field
            .conductors()
            .iter()
            .map(|s| (s.triangle, &s.conductor))
            .collect();
        let regions: BTreeMap<_, _> = field
            .conductors()
            .iter()
            .map(|s| (s.triangle, s.region))
            .collect();
        let initial_sample = prepared
            .resolve(thermal, || (*cancel.borrow_mut())())
            .map_err(material)?;
        let (accepted_sources, _) = self.materials_at_thermal(
            thermal,
            field.problem(),
            &conductors,
            &initial_sample.source,
            &routes,
            &cancel,
        )?;
        if accepted_sources.len() != field.conductors().len()
            || accepted_sources
                .iter()
                .zip(field.conductors())
                .any(|((triangle, c), s)| {
                    *triangle != s.triangle
                        || c.source_basis() != s.conductor.source_basis()
                        || c.source_key() != s.conductor.source_key()
                        || !same_material(c.material(), s.conductor.material())
                })
        {
            return Err(invalid(
                "accepted conductor sources must match actual thermal state",
            ));
        }
        let dt = target_s - thermal.time_s();
        let (thermal, field, report) = cell_driver_candidate(
            thermal,
            target_s,
            controls,
            initial_energy_increments_j,
            |trial, _| {
                check(&cancel)?;
                let mut coefficients = prepared
                    .resolve(trial, || (*cancel.borrow_mut())())
                    .map_err(material)?;
                let (sources, magnetic) = self.materials_at_thermal(
                    trial,
                    field.problem(),
                    &conductors,
                    &coefficients.source,
                    &routes,
                    &cancel,
                )?;
                caloric.validate_sources(trial, &magnetic, &cancel)?;
                let proposal = field
                    .step_with_materials_candidate(
                        target_s,
                        magnetic.iter().map(|s| s.material.clone()).collect(),
                        sources
                            .into_iter()
                            .map(|(triangle, conductor)| EddyConductor {
                                triangle,
                                region: regions[&triangle],
                                conductor,
                            })
                            .collect(),
                        impressed_jz_a_m2,
                        assembly,
                        field_controls,
                        &mut || (*cancel.borrow_mut())(),
                    )
                    .map_err(induction)?;
                let exchange = caloric.deposit(
                    trial,
                    proposal
                        .step()
                        .magnetic_cells
                        .iter()
                        .map(|cell| MagneticEnergyCell {
                            triangle: cell.triangle,
                            flux_density_t: cell.flux_density_t,
                            free_energy_j: cell.free_energy_j,
                            parameter_work_j: cell.parameter_work_j,
                        }),
                    &mut coefficients.sources_w,
                    dt,
                    controls.total_energy_tolerance_j,
                    &cancel,
                )?;
                Ok(DrivenCoefficients {
                    circuit: proposal,
                    sources_w: coefficients.sources_w,
                    conductances_w_k: coefficients.conductances_w_k,
                    circuit_trials: 1,
                    source: EddyCaloricSource {
                        material: EddyThermalSample {
                            network: coefficients.source,
                            magnetic,
                        },
                        fields: exchange.fields,
                        field_internal_energy_j: exchange.field_internal_energy_j,
                        reversible_heat_j: exchange.reversible_heat_j,
                    },
                })
            },
            |proposal: &PlanarEddyCandidate<'_>, dt, checkpoint| {
                deposit_heat(proposal, &routes, dt, checkpoint)
            },
            || (*cancel.borrow_mut())(),
        )?;
        let external_heat_j =
            finite(report.non_electrical_heat_j - report.source.reversible_heat_j)?;
        let total_energy_residual_j = finite(
            report.source.field_internal_energy_j[1] - report.source.field_internal_energy_j[0]
                + report.thermal.stored_energy_change_j
                - report.electrical.supplied_work_j
                - external_heat_j,
        )?;
        if total_energy_residual_j.abs() > controls.total_energy_tolerance_j {
            return Err(invalid(
                "physical magnetic and caloric internal energy do not close external supply",
            ));
        }
        check(&cancel)?;
        Ok(EddyCaloricCandidate {
            thermal,
            field,
            report: EddyCaloricReport {
                coupled: report,
                external_heat_j,
                total_energy_residual_j,
            },
        })
    }
}
