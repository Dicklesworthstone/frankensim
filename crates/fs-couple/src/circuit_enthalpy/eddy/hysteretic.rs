//! Spatial hysteresis, eddy currents and finite material heat on one timeline.
//!
//! The native field owner retains every vector polarization and computes
//! pinning and Joule heat separately. Both physical losses reach their original
//! thermal cells exactly once; backward-Euler numerical dissipation is not heat.
//! Each thermal trial restarts from the SAME accepted field/history. Constant
//! magnetic fits retain their complete source-bounded temperature support,
//! while electrical conductivity responds to the actual endpoint temperature.
//!
//! This is fixed-geometry, fixed-fit, first-order coupling. Changing hysteresis
//! coefficients needs a separate entropy/parameter-work/history law. No fitted
//! core-loss map is added to either resolved physical loss.

use super::{
    CircuitEnthalpyError, CircuitEnthalpyReport, CircuitEnthalpySolve, DrivenCoefficients,
    DriverAccount, EddyThermalSample, ThermalDriver, cell_driver_candidate, deposit_cell_heat,
    finite,
};
use crate::circuit_enthalpy::{
    magnetic::{
        Error, MagneticRoutes, MagneticThermalError, MagneticThermalNetwork, MagneticThermalSample,
        check, invalid, material, same_material,
    },
    material::{MaterialNetwork, MaterialNetworkError, MaterialNetworkSample},
};
use fs_em::planar::{
    eddy::{
        EddyConductor,
        hysteretic::{
            HystereticEddyControls, HystereticEddyError, HystereticEddyStep, PlanarHystereticEddy,
            PlanarHystereticEddyCandidate,
        },
    },
    nonlinear::hysteretic::PlanarHystereticMagnetostatics,
};
use fs_material::{
    conductor::ResolvedConductorState, field_storage::admit_source_temperature_refresh,
};
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation};
use std::{cell::RefCell, collections::BTreeMap};

/// Atomic field/polarization/thermal candidate and the original physical audit.
pub struct HystereticEddyEnthalpyCandidate<'t, 'm> {
    thermal: DepositionSimulation<'t>,
    field: PlanarHystereticEddy<'m>,
    report: CircuitEnthalpyReport<EddyThermalSample, HystereticEddyStep>,
}
impl<'t, 'm> HystereticEddyEnthalpyCandidate<'t, 'm> {
    pub fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub fn field(&self) -> &PlanarHystereticEddy<'m> {
        &self.field
    }
    /// The shared driver's `deposited_joule_j` counter contains total deposited
    /// physical heat here: native Joule PLUS pinning heat. The original field
    /// receipt retains both contributions and numerical dissipation separately.
    pub fn report(&self) -> &CircuitEnthalpyReport<EddyThermalSample, HystereticEddyStep> {
        &self.report
    }
    /// Native Joule plus pinning heat already deposited into finite cells [J].
    /// This excludes the field integrator's separately reported numerical loss.
    pub fn deposited_heat_j(&self) -> f64 {
        self.report.deposited_joule_j
    }
    /// Publish both histories together. All reported physical heat is already
    /// deposited; numerical dissipation has never entered a thermal cell.
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        PlanarHystereticEddy<'m>,
        CircuitEnthalpyReport<EddyThermalSample, HystereticEddyStep>,
    ) {
        (self.thermal, self.field, self.report)
    }
}

impl<'m> ThermalDriver for PlanarHystereticEddyCandidate<'m> {
    type Model = PlanarHystereticEddy<'m>;
    type Step = HystereticEddyStep;
    fn account(&self) -> DriverAccount<'_> {
        let step = self.step();
        DriverAccount {
            resistors: &[],
            contacts: &[],
            heat_j: step.heat_j,
            numerical_dissipation_j: step.numerical_dissipation_j,
            stored_change_j: step.magnetic.stored_energy_change_j,
            supplied_j: step.supplied_work_j,
        }
    }
    fn into_parts(self) -> (Self::Model, Self::Step) {
        PlanarHystereticEddyCandidate::into_parts(self)
    }
}

fn induction(error: HystereticEddyError) -> Error {
    match error {
        HystereticEddyError::Cancelled => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(MagneticThermalError::NonlinearInduction(other)),
    }
}

fn deposit_heat(
    proposal: &PlanarHystereticEddyCandidate<'_>,
    routes: &MagneticRoutes,
    dt: f64,
    checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
) -> Result<BTreeMap<CellId, f64>, Error> {
    let mut powers = deposit_cell_heat(&proposal.step().cells, routes, dt, checkpoint)?;
    if proposal.step().magnetic.cells.len() != routes.magnetic.len() {
        return Err(invalid(
            "hysteresis heat changed magnetic triangle correspondence",
        ));
    }
    for loss in &proposal.step().magnetic.cells {
        checkpoint()?;
        let heat = finite(loss.irreversible_heat_j)?;
        if heat < 0.0 {
            return Err(invalid("hysteresis heat must be nonnegative"));
        }
        if heat == 0.0 {
            continue;
        }
        let cell = routes
            .magnetic
            .get(&loss.triangle)
            .copied()
            .flatten()
            .ok_or_else(|| {
                invalid("every dissipating magnetic triangle needs an actual thermal receiver")
            })?;
        let power = finite(heat / dt)?;
        let previous = *powers.get(&cell).unwrap_or(&0.0);
        let combined = finite(previous + power)?;
        if power == 0.0 || combined == previous {
            return Err(invalid("pinning heat disappeared in receiving-cell power"));
        }
        powers.insert(cell, combined);
    }
    Ok(powers)
}

impl MagneticThermalNetwork<'_> {
    fn hysteretic_materials_at(
        &self,
        thermal: &DepositionSimulation<'_>,
        problem: &PlanarHystereticMagnetostatics,
        conductors: &BTreeMap<usize, &ResolvedConductorState>,
        sample: &MaterialNetworkSample,
        routes: &MagneticRoutes,
        cancelled: &RefCell<impl FnMut() -> bool>,
    ) -> Result<
        (
            Vec<(usize, ResolvedConductorState)>,
            Vec<MagneticThermalSample>,
        ),
        Error,
    > {
        let originals: Vec<_> = problem.materials().iter().map(|s| s.material()).collect();
        let (sources, mut magnetic) = self.materials_at_using(
            thermal,
            &originals,
            conductors,
            sample,
            routes,
            cancelled,
            |_, original, snapshot, binding| {
                if original.card_identity() != snapshot.card_identity()
                    || original.material() != snapshot.material()
                {
                    return Err(invalid(
                        "magnetic fit, conductor and thermal cell must use the same original card",
                    ));
                }
                for property in original.properties() {
                    binding
                        .point
                        .admit_static(property)
                        .map_err(|e| material(MaterialNetworkError::Material(e)))?;
                }
                Ok(())
            },
        )?;
        for sample in &mut magnetic {
            check(cancelled)?;
            let original = &problem.materials()[sample.triangle];
            let refreshed = original
                .at_query_point(&sample.material.query())
                .map_err(|e| CircuitEnthalpyError::Source(MagneticThermalError::Magnetic(e)))?;
            admit_source_temperature_refresh(original.material(), refreshed.material())
                .map_err(|e| CircuitEnthalpyError::Source(MagneticThermalError::Magnetic(e)))?;
            // Electrical bundles may select different property subsets. Keep
            // the complete fitted magnetic claims in the magnetic receipt.
            sample.material = refreshed.material().clone();
        }
        Ok((sources, magnetic))
    }

    /// Continue the same spatial vector magnetic memory and finite heat cells.
    /// Conducting triangles require actual receivers. Insulating magnetic
    /// triangles also require receivers whenever they dissipate pinning heat;
    /// an explicitly fixed region cannot silently absorb an unaccounted loss.
    ///
    /// The full caloric temperature interval must be covered by each original
    /// constant fit. Actual conductivity is refreshed from the same original
    /// source at each thermal trial. All field and thermal acceptance gates
    /// pass before either history can be retained. Refine time for this
    /// first-order fixed-fit magnetoquasistatic/thermal model.
    #[allow(clippy::too_many_arguments)]
    pub fn step_hysteretic_eddy_candidate<'t, 'm>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        field: &PlanarHystereticEddy<'m>,
        impressed_jz_a_m2: &[f64],
        controls: CircuitEnthalpySolve,
        field_controls: HystereticEddyControls,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<HystereticEddyEnthalpyCandidate<'t, 'm>, Error> {
        let cancel = RefCell::new(cancelled);
        check(&cancel)?;
        crate::circuit_enthalpy::admit_driver_unknowns(
            thermal,
            controls,
            initial_energy_increments_j,
        )?;
        if field.time_s().to_bits() != thermal.time_s().to_bits()
            || self.cells.len() != thermal.cells().len()
        {
            return Err(invalid(
                "matching hysteretic field/thermal clocks and complete cell bindings required",
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
        let routes = self.routes_for_mesh(thermal, field.problem().mesh(), &triangles, &cancel)?;
        for (&triangle, &cell) in &routes.magnetic {
            check(&cancel)?;
            if let Some(cell) = cell {
                let binding = self.cells.iter().find(|s| s.cell == cell).ok_or_else(|| {
                    invalid("magnetic receiver lacks its original thermal binding")
                })?;
                let knots = binding.curve.knots();
                let interval = [knots[0].temperature_k, knots[knots.len() - 1].temperature_k];
                field.problem().materials()[triangle]
                    .admit_temperature_interval(binding.point.temperature_axis, interval)
                    .map_err(|e| CircuitEnthalpyError::Source(MagneticThermalError::Magnetic(e)))?;
            }
        }
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
        let (sources, _) = self.hysteretic_materials_at(
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
                .any(|((triangle, conductor), original)| {
                    *triangle != original.triangle
                        || conductor.source_basis() != original.conductor.source_basis()
                        || conductor.source_key() != original.conductor.source_key()
                        || !same_material(conductor.material(), original.conductor.material())
                })
        {
            return Err(invalid(
                "accepted hysteretic conductor sources must match the actual thermal state",
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
                let (sources, magnetic) = self.hysteretic_materials_at(
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
            |proposal, dt, checkpoint| deposit_heat(proposal, &routes, dt, checkpoint),
            || (*cancel.borrow_mut())(),
        )?;
        check(&cancel)?;
        Ok(HystereticEddyEnthalpyCandidate {
            thermal,
            field,
            report,
        })
    }
}
