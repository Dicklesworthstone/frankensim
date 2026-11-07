//! Spatial magnetic memory and physical loss channels on one thermal clock.
//!
//! Every thermal trial restarts the SAME accepted potential and vector play
//! history. The native field solve computes resolved Ohmic, pinning and intrinsic
//! viscous switching heat; each reaches its original finite cells exactly once.
//! The endpoint integration remainder remains numerical work, never heat.
//! Source-backed conductivity changes with
//! temperature, while every fitted magnetic coefficient must remain constant
//! within its explicitly supported temperature interval.

use super::{
    CircuitEnthalpyError, CircuitEnthalpyReport, CircuitEnthalpySolve, CircuitThermalCoefficients,
    DrivenCoefficients, DriverAccount, EddyThermalSample, ThermalDriver, cell_driver_candidate,
    deposit_cell_heat, finite, native_error,
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
            HystereticEddyControls, HystereticEddyStep, PlanarHystereticEddy,
            PlanarHystereticEddyCandidate,
        },
    },
};
use fs_material::{
    conductor::spatial::SpatialConductorMaterial,
    field_storage::{
        admit_source_temperature_refresh, hysteresis::source::ResolvedVectorHysteresis,
    },
};
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation};
use std::{cell::RefCell, collections::BTreeMap};

/// One accepted vector-memory field and finite caloric history, published together.
pub struct HystereticEddyEnthalpyCandidate<'t, 'm, S = MaterialNetworkSample> {
    thermal: DepositionSimulation<'t>,
    field: PlanarHystereticEddy<'m>,
    pub(in crate::circuit_enthalpy) report:
        CircuitEnthalpyReport<EddyThermalSample<S>, HystereticEddyStep>,
}

impl<'t, 'm, S> HystereticEddyEnthalpyCandidate<'t, 'm, S> {
    pub fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub fn field(&self) -> &PlanarHystereticEddy<'m> {
        &self.field
    }
    pub fn report(&self) -> &CircuitEnthalpyReport<EddyThermalSample<S>, HystereticEddyStep> {
        &self.report
    }
    /// Native Joule, pinning and viscous heat deposited into finite cells [J].
    /// Numerical dissipation is excluded from this physical heat total.
    pub fn deposited_heat_j(&self) -> f64 {
        self.report.deposited_joule_j
    }
    /// Physical Joule, pinning and viscous heat is already deposited. Numerical loss
    /// has not been deposited and must not be added to any thermal cell.
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        PlanarHystereticEddy<'m>,
        CircuitEnthalpyReport<EddyThermalSample<S>, HystereticEddyStep>,
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
            // The constitutive owner retains a stable increment beside a large
            // magnetic bias; subtracting endpoint totals would discard it.
            stored_change_j: step.magnetic.stored_energy_change_j,
            supplied_j: step.supplied_work_j,
        }
    }
    fn into_parts(self) -> (Self::Model, Self::Step) {
        PlanarHystereticEddyCandidate::into_parts(self)
    }
}

// Switching also heats electrically insulating magnetic material. Conductor-only
// Joule records cannot stand in for the complete original magnetic cell list.
fn deposit_heat(
    proposal: &PlanarHystereticEddyCandidate<'_>,
    routes: &MagneticRoutes,
    dt: f64,
    checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
) -> Result<BTreeMap<CellId, f64>, Error> {
    let mut powers = deposit_cell_heat(&proposal.step().cells, routes, dt, checkpoint)?;
    deposit_constitutive_heat(
        &proposal.step().magnetic.cells,
        routes,
        dt,
        &mut powers,
        checkpoint,
    )?;
    Ok(powers)
}

// Field-driven and circuit-driven magnetic histories deposit the identical
// original per-cell dry plus viscous receipt, including insulating cells.
// Use the total physical receipt, never the half-viscous variational objective.
pub(in crate::circuit_enthalpy) fn deposit_constitutive_heat(
    losses: &[fs_em::planar::nonlinear::hysteretic::PlanarHystereticCellStep],
    routes: &MagneticRoutes,
    dt: f64,
    powers: &mut BTreeMap<CellId, f64>,
    checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
) -> Result<(), Error> {
    if losses.len() != routes.magnetic.len() {
        return Err(invalid(
            "constitutive heat changed original magnetic cell correspondence",
        ));
    }
    for loss in losses {
        checkpoint()?;
        let receiver = routes
            .magnetic
            .get(&loss.triangle)
            .ok_or_else(|| invalid("constitutive heat lacks an original magnetic region"))?;
        let heat = finite(loss.irreversible_heat_j)?;
        if heat < 0.0 {
            return Err(invalid("negative native magnetic constitutive heat"));
        }
        if heat == 0.0 {
            continue;
        }
        let cell = receiver.ok_or_else(|| {
            invalid("every dissipating magnetic triangle needs an actual thermal heat destination")
        })?;
        let power = finite(heat / dt)?;
        let previous = *powers.get(&cell).unwrap_or(&0.0);
        let combined = finite(previous + power)?;
        if power == 0.0 || combined == previous {
            return Err(invalid(
                "local constitutive heat disappeared in receiving-cell sum",
            ));
        }
        powers.insert(cell, combined);
    }
    Ok(())
}

impl MagneticThermalNetwork<'_> {
    // Source form and the complete supported interval do not depend on a
    // thermal Newton trial. Admit them once for both field and winding drivers.
    pub(in crate::circuit_enthalpy) fn admit_hysteretic_temperature_support(
        &self,
        materials: &[ResolvedVectorHysteresis],
        routes: &MagneticRoutes,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<(), Error> {
        check(cancel)?;
        let bindings: BTreeMap<_, _> = self.cells.iter().map(|cell| (cell.cell, cell)).collect();
        for (&triangle, &cell) in &routes.magnetic {
            check(cancel)?;
            let Some(cell) = cell else { continue };
            let binding = bindings
                .get(&cell)
                .ok_or_else(|| invalid("magnetic receiver lacks its original thermal binding"))?;
            let law = materials
                .get(triangle)
                .ok_or_else(|| invalid("magnetic material changed triangle correspondence"))?;
            for property in law.material().properties() {
                binding
                    .point
                    .admit_static(property)
                    .map_err(|e| material(MaterialNetworkError::Material(e)))?;
            }
            let knots = binding.curve.knots();
            let (Some(first), Some(last)) = (knots.first(), knots.last()) else {
                return Err(invalid(
                    "magnetic receiver lacks its original caloric interval",
                ));
            };
            let axis = binding.point.temperature_axis;
            let mut interval = [first.temperature_k, last.temperature_k];
            if let Some((_, temperature)) = law
                .material()
                .query_point()
                .iter()
                .find(|(name, _)| name == axis)
            {
                // The retained fit can predate this thermal endpoint. Cover
                // its original query as well as every admitted trial endpoint.
                interval[0] = interval[0].min(*temperature);
                interval[1] = interval[1].max(*temperature);
            }
            law.admit_temperature_interval(axis, interval)
                .map_err(|e| CircuitEnthalpyError::Source(MagneticThermalError::Magnetic(e)))?;
        }
        check(cancel)?;
        Ok(())
    }

    // Each field triangle retains its own fitted magnetic source. Electrical
    // bundles may select other properties from that SAME card, so re-resolve
    // the magnetic requirements independently rather than inspecting sigma's
    // receipt as if it established the pinning or saturation coefficients.
    pub(in crate::circuit_enthalpy) fn hysteretic_materials_at<C: SpatialConductorMaterial>(
        &self,
        thermal: &DepositionSimulation<'_>,
        materials: &[ResolvedVectorHysteresis],
        conductors: &BTreeMap<usize, &C>,
        sample: &MaterialNetworkSample,
        routes: &MagneticRoutes,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<
        (
            Vec<(usize, C)>,
            Vec<MagneticThermalSample>,
        ),
        Error,
    > {
        let originals: Vec<_> = materials.iter().map(|m| m.material()).collect();
        let (sources, mut magnetic) = self.materials_at_using(
            thermal, &originals, conductors, sample, routes, cancel,
            |_, original, snapshot, binding| {
                if original.card_identity() != snapshot.card_identity()
                    || original.material() != snapshot.material()
                {
                    return Err(invalid("hysteretic and thermal cells must retain the same original material card"));
                }
                for ((axis, before), (_, after)) in original.query_point().iter()
                    .zip(snapshot.query_point())
                {
                    if before.to_bits() != after.to_bits() && axis != binding.point.temperature_axis {
                        return Err(invalid("only the owned thermal temperature may change a hysteretic source query"));
                    }
                }
                Ok(())
            },
        )?;
        // A thermal/electrical bundle can select different properties from the
        // same card. Keep the complete fitted magnetic claims at the actual
        // endpoint in the magnetic receipt, resolving them only once per cell.
        for sample in &mut magnetic {
            check(cancel)?;
            let original = &materials[sample.triangle];
            let error = |e| CircuitEnthalpyError::Source(MagneticThermalError::Magnetic(e));
            let refreshed = original
                .at_query_point(&sample.material.query())
                .map_err(error)?;
            admit_source_temperature_refresh(original.material(), refreshed.material())
                .map_err(error)?;
            sample.material = refreshed.material().clone();
        }
        Ok((sources, magnetic))
    }

    /// Continue vector hysteresis, eddy diffusion and finite heat storage at
    /// one endpoint. Native source conductivity responds to the solved cell
    /// temperature. All thermal trials restart the same accepted field and
    /// vector polarizations, and a refusal publishes neither state.
    ///
    /// Every heated magnetic triangle needs an explicit original thermal
    /// receiver, including insulating core cells. Fixed regions may dissipate
    /// no pinning or viscous heat. Fits, geometry and solid caloric charts remain fixed;
    /// changing magnetic coefficients, magnetic entropy, phase transitions and
    /// empirical core-loss maps are not supplied by this first-order model.
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
        let network = MaterialNetwork {
            cells: self.cells,
            resistors: &[],
            links: self.links,
        };
        let prepared = network
            .prepare_with_bases(thermal, &[], || (*cancel.borrow_mut())())
            .map_err(material)?;
        self.step_hysteretic_eddy_resolved_candidate(
            thermal,
            target_s,
            field,
            impressed_jz_a_m2,
            controls,
            field_controls,
            initial_energy_increments_j,
            |trial| {
                prepared
                    .resolve(trial, || (*cancel.borrow_mut())())
                    .map_err(material)
            },
            |sample| sample,
            &cancel,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub(in crate::circuit_enthalpy) fn step_hysteretic_eddy_resolved_candidate<
        't,
        'm,
        S,
        E: From<MagneticThermalError>,
    >(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        field: &PlanarHystereticEddy<'m>,
        impressed_jz_a_m2: &[f64],
        controls: CircuitEnthalpySolve,
        field_controls: HystereticEddyControls,
        initial_energy_increments_j: &[f64],
        resolve: impl Fn(
            &DepositionSimulation<'_>,
        ) -> Result<CircuitThermalCoefficients<S>, CircuitEnthalpyError<E>>,
        snapshot: impl Fn(&S) -> &MaterialNetworkSample,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<HystereticEddyEnthalpyCandidate<'t, 'm, S>, CircuitEnthalpyError<E>> {
        check(cancel).map_err(native_error)?;
        crate::circuit_enthalpy::admit_driver_unknowns(
            thermal,
            controls,
            initial_energy_increments_j,
        )?;
        if field.time_s().to_bits() != thermal.time_s().to_bits()
            || self.cells.len() != thermal.cells().len()
        {
            return Err(native_error(invalid(
                "matching hysteretic field/thermal clocks and complete cell bindings required",
            )));
        }
        let triangles: Vec<_> = field.conductors().iter().map(|s| s.triangle).collect();
        let routes = self
            .routes_for_mesh(thermal, field.problem().mesh(), &triangles, cancel)
            .map_err(native_error)?;
        self.admit_hysteretic_temperature_support(field.problem().materials(), &routes, cancel)
            .map_err(native_error)?;
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
        let initial = resolve(thermal)?;
        let (sources, _) = self
            .hysteretic_materials_at(
                thermal,
                field.problem().materials(),
                &conductors,
                snapshot(&initial.source),
                &routes,
                cancel,
            )
            .map_err(native_error)?;
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
            return Err(native_error(invalid(
                "accepted hysteretic conductor sources must match the actual thermal state",
            )));
        }
        let (thermal, field, report) = cell_driver_candidate(
            thermal,
            target_s,
            controls,
            initial_energy_increments_j,
            |state, _| {
                check(cancel).map_err(native_error)?;
                let coefficients = resolve(state)?;
                let (sources, magnetic) = self
                    .hysteretic_materials_at(
                        state,
                        field.problem().materials(),
                        &conductors,
                        snapshot(&coefficients.source),
                        &routes,
                        cancel,
                    )
                    .map_err(native_error)?;
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
                    .map_err(|e| native_error(super::nonlinear::induction(e)))?;
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
            |proposal, dt, checkpoint| {
                deposit_heat(proposal, &routes, dt, checkpoint).map_err(native_error)
            },
            || (*cancel.borrow_mut())(),
        )?;
        check(cancel).map_err(native_error)?;
        Ok(HystereticEddyEnthalpyCandidate {
            thermal,
            field,
            report,
        })
    }
}
