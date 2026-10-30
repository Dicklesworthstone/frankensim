//! Spatial magnetic induction and heat on the original finite enthalpy clock.
//!
//! Every Newton trial restarts from the accepted original A_z. Conductivity is
//! refreshed at the actual thermal endpoint; the same candidate's sigma*E^2
//! heat enters the existing cells exactly once. Magnetic energy, external field
//! source work, conduction, boundary heat and total enthalpy share one audit.
//! Geometry and magnetic permeability remain fixed; no core-loss map is added.

/// Reversible nonlinear induction with physical Joule heating and a separately
/// audited backward-Euler numerical loss on the same finite thermal timeline.
pub mod nonlinear;
/// Reversible permeability/temperature feedback on these same native owners.
pub mod thermal;

use super::{
    CircuitEnthalpyError, CircuitEnthalpyReport, CircuitEnthalpySolve, DrivenCoefficients,
    DriverAccount, ThermalDriver, cell_driver_candidate, finite,
    magnetic::{
        Error, MagneticThermalError, MagneticThermalNetwork, MagneticThermalSample, MagneticRoutes, check, invalid,
        material, same_material,
    },
    material::{MaterialNetwork, MaterialNetworkSample},
};
use fs_em::planar::eddy::{
    EddyConductor, EddyControls, EddyError, EddyStep, PlanarEddy, PlanarEddyCandidate,
};
use fs_thermal::enthalpy::{CellId,DepositionError,DepositionSimulation};
use std::{cell::RefCell, collections::BTreeMap};

#[derive(Clone, Debug, PartialEq)]
pub struct EddyThermalSample<S = MaterialNetworkSample> {
    pub network: S,
    pub magnetic: Vec<MagneticThermalSample>,
}
pub struct EddyEnthalpyCandidate<'t, 'm> {
    thermal: DepositionSimulation<'t>,
    field: PlanarEddy<'m>,
    report: CircuitEnthalpyReport<EddyThermalSample, EddyStep>,
}
impl<'t, 'm> EddyEnthalpyCandidate<'t, 'm> {
    pub fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub fn field(&self) -> &PlanarEddy<'m> {
        &self.field
    }
    pub fn report(&self) -> &CircuitEnthalpyReport<EddyThermalSample, EddyStep> {
        &self.report
    }
    /// Publish BOTH histories together; returned eddy heat is already deposited.
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        PlanarEddy<'m>,
        CircuitEnthalpyReport<EddyThermalSample, EddyStep>,
    ) {
        (self.thermal, self.field, self.report)
    }
}
impl<'m> ThermalDriver for PlanarEddyCandidate<'m> {
    type Model = PlanarEddy<'m>;
    type Step = EddyStep;
    fn account(&self) -> DriverAccount<'_> {
        let s = self.step();
        DriverAccount {
            resistors: &[],
            contacts: &[],
            heat_j: s.heat_j,
            numerical_dissipation_j: 0.0,
            stored_change_j: s.stored_after_j - s.stored_before_j - s.parameter_work_j,
            supplied_j: s.supplied_work_j,
        }
    }
    fn into_parts(self) -> (Self::Model, Self::Step) {
        PlanarEddyCandidate::into_parts(self)
    }
}
fn induction(e: EddyError) -> Error {
    match e {
        EddyError::Cancelled => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(MagneticThermalError::Induction(other)),
    }
}

// The fixed-mu and thermal-mu paths deposit exactly the same native sigma*E²
// integrals. Neither owner may replace them with averaged-current heating.
fn deposit_heat(
    proposal: &PlanarEddyCandidate<'_>,
    routes: &MagneticRoutes,
    dt: f64,
    checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
) -> Result<BTreeMap<CellId, f64>, Error> {
    deposit_cell_heat(&proposal.step().cells, routes, dt, checkpoint)
}

// Identical physical sigma E^2 routing for midpoint and nonlinear induction.
// No magnetic numerical dissipation is present in these admitted cell records.
fn deposit_cell_heat(
    losses: &[fs_em::planar::eddy::EddyCellHeat],
    routes: &MagneticRoutes,
    dt: f64,
    checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
) -> Result<BTreeMap<CellId, f64>, Error> {
    if losses.len() != routes.conductors.len() {
        return Err(invalid("native eddy heat changed triangle correspondence"));
    }
    let mut powers=BTreeMap::new();
    for loss in losses {
        checkpoint()?;
        let cell=routes.conductors.get(&loss.source.triangle)
            .ok_or_else(||invalid("eddy heat lacks an original receiving cell"))?;
        let power=finite(loss.heat_j/dt)?;
        if loss.heat_j<0.0 || loss.heat_j!=0.0 && power==0.0 {
            return Err(invalid("eddy heat disappeared in power conversion"));
        }
        let previous=*powers.get(cell).unwrap_or(&0.0);
        let combined=finite(previous+power)?;
        if power!=0.0 && combined==previous {
            return Err(invalid("local eddy heat disappeared in receiving-cell sum"));
        }
        powers.insert(*cell,combined);
    }
    Ok(powers)
}

impl MagneticThermalNetwork<'_> {
    /// Continue the SAME magnetic and thermal histories on one accepted interval.
    /// All conducting triangles require actual receiving cells. Insulating
    /// triangles may be thermally mapped or explicitly fixed. Original borrowed
    /// caloric laws/cards/source choices, fixed coordinates and volumes stay pinned.
    /// The existing network computes conduction and authored boundary sources;
    /// midpoint induction deposits its actual local integrated heat, not J_avg².
    /// Impressed current is the ENDPOINT profile; the original field's retained
    /// source varies linearly to it and supplies the audited midpoint work.
    ///
    /// First-order coupled accuracy in general: thermal rates are backward Euler
    /// and conductivity uses the endpoint. Refine time and numerical tolerances
    /// independently. No motion, phase transition, variable-mu storage work,
    /// hysteresis, finite winding ends or empirical material validation is inferred.
    #[allow(clippy::too_many_arguments)]
    pub fn step_eddy_candidate<'t, 'm>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        field: &PlanarEddy<'m>,
        impressed_jz_a_m2: &[f64],
        controls: CircuitEnthalpySolve,
        field_controls: EddyControls,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<EddyEnthalpyCandidate<'t, 'm>, Error> {
        let cancel = RefCell::new(cancelled);
        check(&cancel)?;
        super::admit_driver_unknowns(thermal, controls, initial_energy_increments_j)?;
        if field.time_s().to_bits() != thermal.time_s().to_bits()
            || self.cells.len() != thermal.cells().len()
        {
            return Err(invalid(
                "matching magnetic/thermal clocks and complete thermal bindings required",
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
        let initial = prepared
            .resolve(thermal, || (*cancel.borrow_mut())())
            .map_err(material)?;
        let (sources, _) = self.materials_at(
            thermal,
            field.problem(),
            &conductors,
            &initial.source,
            &routes,
            &cancel,
        )?;
        if sources.len() != field.conductors().len()
            || sources
                .iter()
                .zip(field.conductors())
                .any(|((cell, c), s)| {
                    *cell != s.triangle
                        || c.source_basis() != s.conductor.source_basis()
                        || c.source_key() != s.conductor.source_key()
                        || !same_material(c.material(), s.conductor.material())
                })
        {
            return Err(invalid(
                "accepted conductor sources must match the actual thermal state",
            ));
        }
        let (thermal, field, report) = cell_driver_candidate(
            thermal,
            target_s,
            controls,
            initial_energy_increments_j,
            |state, _| {
                check(&cancel)?;
                let coefficients = prepared
                    .resolve(state, || (*cancel.borrow_mut())())
                    .map_err(material)?;
                let (sources, magnetic) = self.materials_at(
                    state,
                    field.problem(),
                    &conductors,
                    &coefficients.source,
                    &routes,
                    &cancel,
                )?;
                let sources = sources
                    .into_iter()
                    .map(|(triangle, conductor)| EddyConductor {
                        triangle,
                        region: regions[&triangle],
                        conductor,
                    })
                    .collect();
                let proposal = field
                    .step_to_candidate(
                        target_s,
                        sources,
                        impressed_jz_a_m2,
                        field_controls,
                        &mut || (*cancel.borrow_mut())(),
                    )
                    .map_err(induction)?;
                Ok(DrivenCoefficients {
                    circuit: proposal,
                    sources_w: coefficients.sources_w,
                    conductances_w_k: coefficients.conductances_w_k,
                    source: EddyThermalSample {
                        network: coefficients.source,
                        magnetic,
                    },
                    circuit_trials: 1,
                })
            },
            |proposal: &PlanarEddyCandidate<'_>, dt, checkpoint| {
                deposit_heat(proposal,&routes,dt,checkpoint)
            },
            || (*cancel.borrow_mut())(),
        )?;
        check(&cancel)?;
        Ok(EddyEnthalpyCandidate {
            thermal,
            field,
            report,
        })
    }
}
