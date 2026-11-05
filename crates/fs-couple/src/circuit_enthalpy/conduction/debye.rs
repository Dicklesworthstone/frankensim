//! Distributed Debye loss and finite caloric storage on one physical interval.
//!
//! The original field history supplies polarization (and admitted Ohmic) heat.
//! Its total local loss enters the existing enthalpy residual exactly once;
//! instantaneous and polarization storage remain electrical energy. Constant
//! source coefficients allow one field solve to serve all thermal trial states.
//! The original claims are nevertheless re-queried at every trial temperature.

use super::*;
use fs_em::electrostatic::{
    ElectrostaticControls, ElectrostaticError,
    dispersion::{
        DebyeFieldStep, DielectricFieldCandidate, DielectricFieldState, PlanarDebyeDielectric,
        PreparedDebyeStep,
    },
};
use fs_material::{
    conductor::ELECTRICAL_CONDUCTIVITY_DIMS,
    field_storage::{FieldStorageError, admit_source_temperature_refresh},
    state_point::{ResolvedMaterialStatePoint, ResolvedScalarProperty},
};

mod heat;

/// The original thermal source results and complete electrical source snapshots
/// revalidated at the accepted thermal endpoint. The immutable field retains its
/// original reference snapshots and coefficients; these receipts establish that
/// the same laws remain valid at the actual receiving temperatures.
#[derive(Debug)]
pub struct DebyeThermalSource<S = MaterialNetworkSample> {
    pub thermal: S,
    /// One source bundle per original electrical cell, in original order.
    /// An explicitly leaky cell retains its conductivity claim in this bundle.
    pub dielectric: Vec<ResolvedMaterialStatePoint>,
}

/// Joint electrical/thermal energy account. The native field step retains the
/// individual loss mechanisms; only their sum is deposited. The legacy shared
/// `coupled.deposited_joule_j` field denotes ALL deposited dielectric losses here.
pub struct DebyeEnthalpyReport<S = MaterialNetworkSample, const D: usize = 2> {
    pub coupled: CircuitEnthalpyReport<DebyeThermalSource<S>, DebyeFieldStep<D>>,
    /// Authored/environmental heat, excluding field loss and enclosure exchange.
    pub external_heat_j: f64,
    /// Computed roundoff net of internal closed-enclosure radiation.
    pub closed_radiation_net_j: f64,
    /// Field storage + caloric storage - electrical work - external heat.
    pub total_energy_residual_j: f64,
}

/// Joint candidate; consume both histories together to continue. A refusal or
/// cancellation leaves input potential, polarization, charge and enthalpy intact.
pub struct DebyeEnthalpyCandidate<
    't,
    S = MaterialNetworkSample,
    F = PlanarDebyeDielectric,
    const D: usize = 2,
> {
    thermal: DepositionSimulation<'t>,
    field: F,
    report: DebyeEnthalpyReport<S, D>,
}
impl<'t, S, F, const D: usize> DebyeEnthalpyCandidate<'t, S, F, D> {
    pub fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub fn field(&self) -> &F {
        &self.field
    }
    pub fn report(&self) -> &DebyeEnthalpyReport<S, D> {
        &self.report
    }
    pub fn into_parts(self) -> (DepositionSimulation<'t>, F, DebyeEnthalpyReport<S, D>) {
        (self.thermal, self.field, self.report)
    }
}

// A single physical field candidate is borrowed through thermal Newton trials.
// Only joint acceptance clones its endpoint; there is no repeated field solve,
// repeated history advance, local-pole replay or per-stencil state allocation.
struct FieldDriver<'a, F, const D: usize>(&'a DielectricFieldCandidate<F, D>);
impl<F: DielectricFieldState, const D: usize> ThermalDriver for FieldDriver<'_, F, D> {
    type Model = F;
    type Step = DebyeFieldStep<D>;
    fn account(&self) -> DriverAccount<'_> {
        let step = self.0.step();
        DriverAccount {
            resistors: &[],
            contacts: &[],
            heat_j: step.heat_j,
            numerical_dissipation_j: 0.0,
            stored_change_j: step.stored_change_j,
            supplied_j: step.supplied_work_j,
        }
    }
    fn into_parts(self) -> (Self::Model, Self::Step) {
        (self.0.field().clone(), self.0.step().clone())
    }
}
fn field(error: ElectrostaticError) -> Error {
    match error {
        ElectrostaticError::Cancelled => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(ConductionThermalError::Debye(other)),
    }
}
fn storage(error: FieldStorageError) -> Error {
    field(ElectrostaticError::Material(error))
}

fn admit_constant_source(
    binding: &MaterialCell<'_>,
    property: &ResolvedScalarProperty,
) -> Result<(), Error> {
    binding
        .point
        .admit_static(property)
        .map_err(|e| material(MaterialNetworkError::Material(e)))?;
    if !property.source_axes().is_empty()
        || property.answer().receipt.decision
            != fs_matdb::EvaluationDecision::ConstantWithinValidity
    {
        return Err(invalid(
            "dielectric thermal bridge requires validity-wide constant storage, relaxation and DC conductivity",
        ));
    }
    Ok(())
}

fn deposit<F: DielectricFieldState, const D: usize>(
    proposal: &FieldDriver<'_, F, D>,
    dt_s: f64,
    routes: &[CellId],
    checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
) -> Result<BTreeMap<CellId, f64>, Error> {
    let heat = &proposal.0.step().cell_heat_j;
    if heat.len() != routes.len() {
        return Err(invalid(
            "Debye heat lost its original electrical-cell correspondence",
        ));
    }
    let mut powers = BTreeMap::new();
    for (&heat_j, &cell) in heat.iter().zip(routes) {
        checkpoint()?;
        if !heat_j.is_finite() || heat_j < 0.0 {
            return Err(invalid(
                "Debye physical heat must be finite and nonnegative",
            ));
        }
        let power = finite(heat_j / dt_s)?;
        let previous = *powers.get(&cell).unwrap_or(&0.0);
        let combined = finite(previous + power)?;
        if heat_j > 0.0 && (power == 0.0 || combined == previous) {
            return Err(invalid(
                "Debye heat disappeared in its receiving thermal cell",
            ));
        }
        powers.insert(cell, combined);
    }
    Ok(powers)
}

impl ConductionThermalNetwork<'_> {
    /// Planar convenience entry point: solve the original field once, then
    /// advance the original thermal bodies with its physical loss. It delegates
    /// to [`Self::step_dielectric_field_candidate`], with the same source,
    /// fixed-volume, material ownership and fully solid chart requirements.
    ///
    /// `prepared` is the field owner's reusable time-step operator. Its interval
    /// must end at `target_s`; field and thermal input clocks must agree. Native
    /// midpoint polarization work and heat are integrated exactly once, and
    /// native backward-Euler thermal transport supplies first-order coupled
    /// time accuracy. `maximum_circuit_trials` bounds thermal material residual
    /// evaluations; the field itself is solved once. `controls.circuit` is unused.
    #[allow(clippy::too_many_arguments)]
    pub fn step_debye_candidate<'t>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        dielectric: &PlanarDebyeDielectric,
        prepared: &PreparedDebyeStep,
        endpoint_voltage_v: &[(usize, f64)],
        field_controls: ElectrostaticControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<DebyeEnthalpyCandidate<'t>, Error> {
        let cancel = RefCell::new(cancelled);
        check(&cancel)?;
        super::super::admit_driver_unknowns(thermal, controls, initial_energy_increments_j)?;
        if dielectric.time_s().to_bits() != thermal.time_s().to_bits()
            || (dielectric.time_s() + prepared.time_step_s()).to_bits() != target_s.to_bits()
        {
            return Err(invalid(
                "matching Debye/thermal clocks and prepared interval required",
            ));
        }
        let proposal = dielectric
            .step_candidate(prepared, endpoint_voltage_v, field_controls, &mut || {
                (*cancel.borrow_mut())()
            })
            .map_err(field)?;
        self.step_dielectric_field_candidate(
            thermal,
            &proposal,
            controls,
            initial_energy_increments_j,
            || (*cancel.borrow_mut())(),
        )
    }

    /// Advance finite thermal bodies using one native dielectric field
    /// candidate, without publishing that field separately. Planar,
    /// axisymmetric and tetrahedral fields share this same residual and
    /// physical-volume admission. `ConductionRegion::triangle` indexes the
    /// original native electrical cell, including a tetrahedron in 3-D.
    /// Several electrical cells may cover one complete receiving thermal cell;
    /// passive thermal cells and source-resolved Fourier links remain supported.
    ///
    /// The sealed native candidate supplies its original physical interval,
    /// independent source work, stable storage change and separate polarization
    /// and Ohmic losses. Its total local heat enters the thermal residual once.
    /// Accepted field and caloric histories publish together, after full energy
    /// checks; cancellation or source refusal publishes neither history.
    ///
    /// Storage, relaxation and every explicitly selected DC conductivity must
    /// be validity-wide constant scalar sources, including a zero-conductivity
    /// claim. All original material/source identities and their validity at the
    /// actual receiving temperature are rechecked for every thermal trial;
    /// nonthermal context stays fixed. Fully solid constant-density caloric
    /// charts are required. No phase polarization, changing geometry or
    /// temperature-dependent electrical constitutive law is inferred.
    ///
    /// The existing field candidate is borrowed throughout thermal Newton;
    /// `maximum_circuit_trials` bounds thermal material residual evaluations,
    /// not repeated field solves. `controls.circuit` is unused. Native midpoint
    /// field work and backward-Euler thermal transport retain their independent
    /// time approximations; coupled transport is first order in general.
    pub fn step_dielectric_field_candidate<'t, F: DielectricFieldState, const D: usize>(
        &self,
        thermal: &DepositionSimulation<'t>,
        proposal: &DielectricFieldCandidate<F, D>,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<DebyeEnthalpyCandidate<'t, MaterialNetworkSample, F, D>, Error> {
        let cancel = RefCell::new(cancelled);
        check(&cancel)?;
        super::super::admit_driver_unknowns(thermal, controls, initial_energy_increments_j)?;
        let network = MaterialNetwork {
            cells: self.cells,
            resistors: &[],
            links: self.links,
        };
        let thermal_sources = network
            .prepare_with_bases(thermal, &[], || (*cancel.borrow_mut())())
            .map_err(material)?;
        self.step_dielectric_field_resolved_candidate(
            thermal,
            proposal,
            controls,
            initial_energy_increments_j,
            |trial| {
                thermal_sources
                    .resolve(trial, || (*cancel.borrow_mut())())
                    .map_err(material)
            },
            |sample| sample,
            |_| 0.0,
            &cancel,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn step_dielectric_field_resolved_candidate<'t, S, F: DielectricFieldState, const D: usize>(
        &self,
        thermal: &DepositionSimulation<'t>,
        proposal: &DielectricFieldCandidate<F, D>,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        resolve: impl Fn(&DepositionSimulation<'t>) -> Result<CircuitThermalCoefficients<S>, Error>,
        sample: impl Fn(&S) -> &MaterialNetworkSample,
        radiation_net_w: impl Fn(&S) -> f64,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<DebyeEnthalpyCandidate<'t, S, F, D>, Error> {
        check(cancel)?;
        super::super::admit_driver_unknowns(thermal, controls, initial_energy_increments_j)?;
        let dielectric = proposal.field();
        let target_s = proposal.step().interval_s[1];
        let cells = dielectric.materials().len();
        if proposal.step().interval_s[0].to_bits() != thermal.time_s().to_bits()
            || dielectric.time_s().to_bits() != target_s.to_bits()
            || dielectric.stored_energy_j().to_bits() != proposal.step().stored_after_j.to_bits()
            || self.cells.len() != thermal.cells().len()
            || dielectric.cell_volumes_m3().len() != cells
            || dielectric.conductivity_s_m().len() != cells
            || dielectric.conductivity_properties().len() != cells
        {
            return Err(invalid(
                "matching native dielectric interval, stored state and complete thermal/material cells required",
            ));
        }
        let routes = self.admit_volume_routes(thermal, dielectric.cell_volumes_m3(), cancel)?;
        let bindings: BTreeMap<_, _> = self
            .cells
            .iter()
            .map(|binding| (binding.cell, binding))
            .collect();
        for (index, (&receiver, original)) in routes.iter().zip(dielectric.materials()).enumerate()
        {
            check(cancel)?;
            let binding = bindings
                .get(&receiver)
                .ok_or_else(|| invalid("Debye receiver lacks original caloric material"))?;
            if binding
                .curve
                .knots()
                .iter()
                .any(|knot| knot.liquid_mass_fraction != 0.0)
            {
                return Err(invalid(
                    "Debye heat requires a fully solid fixed-volume caloric chart",
                ));
            }
            let conductivity = dielectric.conductivity_s_m()[index];
            let conductivity_name = dielectric.conductivity_properties()[index].as_deref();
            if !conductivity.is_finite()
                || conductivity < 0.0
                || conductivity_name.is_none() && conductivity != 0.0
            {
                return Err(invalid(
                    "native leakage must retain its nonnegative conductivity source",
                ));
            }
            let names = std::iter::once(original.high_frequency_property()).chain(
                original.pole_properties().iter().flat_map(|pole| {
                    [
                        pole.permittivity_increment.as_str(),
                        pole.relaxation_time.as_str(),
                    ]
                }),
            );
            for name in names {
                let property = original
                    .material()
                    .property(name)
                    .ok_or_else(|| invalid("Debye law lost an original source role"))?;
                admit_constant_source(binding, property)?;
            }
            if let Some(name) = conductivity_name {
                let key = fs_matdb::PropertyKey::new(name, ELECTRICAL_CONDUCTIVITY_DIMS);
                let property = original
                    .material()
                    .property_by_key(&key)
                    .ok_or_else(|| invalid("native leakage lost its original conductivity key"))?;
                admit_constant_source(binding, property)?;
            }
        }
        let initial = resolve(thermal)?;
        self.debye_materials_at(dielectric, sample(&initial.source), &routes, cancel)?;
        let (thermal, field, coupled) = cell_driver_candidate(
            thermal,
            target_s,
            controls,
            initial_energy_increments_j,
            |state, _| {
                check(cancel)?;
                let coefficients = resolve(state)?;
                if !coefficients.resistors.is_empty() {
                    return Err(invalid(
                        "Debye field cannot silently omit graph electrical work",
                    ));
                }
                let dielectric = self.debye_materials_at(
                    dielectric,
                    sample(&coefficients.source),
                    &routes,
                    cancel,
                )?;
                Ok(DrivenCoefficients {
                    circuit: FieldDriver(proposal),
                    sources_w: coefficients.sources_w,
                    conductances_w_k: coefficients.conductances_w_k,
                    source: DebyeThermalSource {
                        thermal: coefficients.source,
                        dielectric,
                    },
                    circuit_trials: 1,
                })
            },
            |driver, dt, checkpoint| deposit(driver, dt, &routes, checkpoint),
            || (*cancel.borrow_mut())(),
        )?;
        let closed_radiation_net_j = finite(
            (target_s - proposal.step().interval_s[0]) * radiation_net_w(&coupled.source.thermal),
        )?;
        let external_heat_j = finite(coupled.non_electrical_heat_j - closed_radiation_net_j)?;
        // Preserve a resolved displacement-work contribution beside large DC
        // power. Reuse the existing compensated thermal-account sum; each
        // source component comes from the SAME native gauge references.
        let mut total = super::super::entropy::Sum::default();
        for term in [
            coupled.electrical.stored_change_j,
            coupled.thermal.stored_energy_change_j,
            -coupled.electrical.supplied_displacement_work_j,
            -coupled.electrical.supplied_conduction_work_j,
            -coupled.non_electrical_heat_j,
            closed_radiation_net_j,
        ] {
            total.add(term)?;
        }
        let total_energy_residual_j = total.total()?;
        if total_energy_residual_j.abs() > controls.total_energy_tolerance_j {
            return Err(invalid(
                "Debye and caloric storage must close external work without enclosure heat",
            ));
        }
        check(cancel)?;
        Ok(DebyeEnthalpyCandidate {
            thermal,
            field,
            report: DebyeEnthalpyReport {
                coupled,
                external_heat_j,
                closed_radiation_net_j,
                total_energy_residual_j,
            },
        })
    }

    fn debye_materials_at<F: DielectricFieldState>(
        &self,
        field: &F,
        sample: &MaterialNetworkSample,
        routes: &[CellId],
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<Vec<ResolvedMaterialStatePoint>, Error> {
        let bindings: BTreeMap<_, _> = self
            .cells
            .iter()
            .map(|binding| (binding.cell, binding))
            .collect();
        let phases: BTreeMap<_, _> = sample
            .cells
            .iter()
            .map(|cell| (cell.cell, cell.phase))
            .collect();
        field
            .materials()
            .iter()
            .zip(routes)
            .map(|(original, receiver)| {
                check(cancel)?;
                let binding = bindings
                    .get(receiver)
                    .ok_or_else(|| invalid("Debye receiver lacks original material"))?;
                let phase = *phases
                    .get(receiver)
                    .ok_or_else(|| invalid("Debye receiver lacks actual endpoint phase"))?;
                let point = binding
                    .point
                    .at(phase)
                    .map_err(|e| material(MaterialNetworkError::Material(e)))?;
                let snapshot = original
                    .material()
                    .at_query_point(binding.point.material, &point)
                    .map_err(|e| {
                        material(MaterialNetworkError::Material(
                            ThermalMaterialError::Material(e),
                        ))
                    })?;
                admit_source_temperature_refresh(original.material(), &snapshot)
                    .map_err(storage)?;
                Ok(snapshot)
            })
            .collect()
    }
}
