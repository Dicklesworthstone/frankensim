//! Finite circuits and loaded 3D plastic solids on ONE enthalpy residual.
//! Native electrical and tetrahedral proposals always start from accepted history.
//! Only physical plastic/Joule heat and admitted reversible exchange enter heat;
//! hardening storage and backward-Euler loss retain their original accounts.
use super::*;
use super::contact::{ContactError, ContactSample, ThermalSupportContact};
use std::cell::Cell;
use crate::circuit_enthalpy::{CircuitEnthalpyError, CircuitEnthalpyInput,
    CircuitEnthalpyReport, CircuitEnthalpySolve, CircuitThermalCoefficients};
use crate::circuit_enthalpy::contacts::{self, ContactThermalCoefficients, ContactThermalSample};
use crate::circuit_enthalpy::contacts::heat::{HeatNetworkError, HeatNetworkSample, PreparedHeatNetwork};
use crate::circuit_enthalpy::contacts::heat::solid::SolidHeatError;
use crate::material_circuit::{MaterialCircuit, MaterialResistor};
use fs_material::state_point::{MaterialPropertySelection, ResolvedMaterialStatePoint};
use crate::material_circuit::contacts::SeriesElectricalContact;
use fs_matdb::ClaimId;
use fs_material::electrical_interface::{
    ELECTRICAL_CONTACT_RESISTANCE_PROPERTY, resolve_electrical_interface_state,
};

/// Electrical response of an existing whole-face mechanical/thermal support.
/// The link fixes both material surfaces, area, pressure and thermal coordinates;
/// this binding only adds its electrical branch, source claim and heat split.
/// The explicit split is a model input, not inferred from either bulk material.
#[derive(Clone, Copy, Debug)]
pub struct ElectricalSupportContact {
    pub link: u64,
    pub branch: usize,
    pub claim: ClaimId,
    pub heat_fraction_to_first: f64,
}

impl ElectricalSupportContact {
    fn resolve(
        &self,
        support: &ThermalSupportContact<'_>,
        sample: &ContactSample,
    ) -> Result<SeriesElectricalContact, ElectricalSolidSourceError> {
        use ElectricalSolidSourceError as Error;
        if self.link != support.link
            || self.link != sample.link
            || !self.heat_fraction_to_first.is_finite()
            || !(0.0..=1.0).contains(&self.heat_fraction_to_first)
        {
            return Err(Error::Invalid(
                "electrical support needs its original link and an explicit heat partition",
            ));
        }
        let point = support
            .point
            .at(sample.phases)
            .map_err(|e| Error::Support(ContactError::State(e)))?
            .with(support.pressure_axis, sample.compression_pa)
            .map_err(|_| Error::Invalid("electrical support pressure coordinate refused"))?;
        let resolved = resolve_electrical_interface_state(
            support.point.interface,
            &point,
            MaterialPropertySelection::PinnedByProperty(vec![(
                ELECTRICAL_CONTACT_RESISTANCE_PROPERTY.to_owned(),
                self.claim,
            )]),
        )
        .map_err(Error::ElectricalSupport)?;
        let claim = support
            .point
            .interface
            .claims()
            .claim(self.claim)
            .ok_or(Error::Invalid("pinned electrical support claim is absent"))?;
        if !claim.value.axes().all(|(axis, dims)| {
            (axis == support.pressure_axis && dims == fs_qty::Pressure::DIMS)
                || support.point.is_response_axis(axis, dims)
        }) {
            return Err(Error::Invalid(
                "electrical support law requires solved pressure, thermal or typed fixed external coordinates",
            ));
        }
        let property = resolved
            .resolved()
            .property(ELECTRICAL_CONTACT_RESISTANCE_PROPERTY)
            .ok_or(Error::Invalid("electrical support property is absent"))?;
        if property
            .angular_frequency_band_rad_s()
            .map_err(|_| Error::Invalid("electrical support frequency descriptor refused"))?
            .is_some_and(|(lower, _)| lower > 0.0)
        {
            return Err(Error::Invalid(
                "electrical support source excludes static response",
            ));
        }
        SeriesElectricalContact::new(self.link, self.branch, resolved, sample.area_m2)
            .map_err(Error::Circuit)
    }
}

/// The original circuit/enthalpy controls plus the work required by native solids.
#[derive(Clone, Copy, Debug)]
pub struct ElectricalSolidSolve {
    pub coupled: CircuitEnthalpySolve,
    /// Includes every stencil and the final fresh endpoint evaluation.
    pub maximum_mechanical_trials: usize,
    /// Native J2 yield-surface/history tolerance, Pa.
    pub yield_tolerance_pa: f64,
}

#[derive(Debug)]
pub enum ElectricalSolidSourceError {
    Heat(HeatNetworkError),
    ElectricalSupport(fs_material::state_point::MaterialStatePointError),
    Circuit(crate::material_circuit::MaterialCircuitError),
    Mechanical(ThermoplasticError),
    Support(ContactError<Infallible>),
    Thermal(DepositionError),
    Invalid(&'static str),
    MechanicalBudget,
}
impl From<DepositionError> for ElectricalSolidSourceError {
    fn from(error: DepositionError) -> Self { Self::Thermal(error) }
}
impl core::fmt::Display for ElectricalSolidSourceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "electrical/solid source: {self:?}")
    }
}
impl std::error::Error for ElectricalSolidSourceError {}

#[derive(Debug)]
pub enum ElectricalSolidError {
    Binding(SolidHeatError<Infallible>),
    Mechanical(ThermoplasticError),
    Support(ContactError<Infallible>),
    Heat(HeatNetworkError),
    Coupled(CircuitEnthalpyError<ElectricalSolidSourceError>),
    Invalid(&'static str),
    Cancelled,
}
impl core::fmt::Display for ElectricalSolidError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "electrical/solid/enthalpy: {self:?}")
    }
}
impl std::error::Error for ElectricalSolidError {}

/// Original source/electrical/thermal reports and independent closed-system work.
pub struct ElectricalSolidReport {
    /// This subsystem report treats mechanical heat and radiative net roundoff
    /// as callback supply. The CLOSED account below removes those internal terms.
    pub circuit_thermal: CircuitEnthalpyReport<ContactThermalSample<HeatNetworkSample>>,
    /// Solved-pressure thermal interfaces, including their original bulk
    /// conductance and actual series heat. Not added again as source power.
    pub support_contacts: Vec<ContactSample>,
    /// Original native J2 energy increments, in tetrahedron connectivity order.
    pub element_energy: Vec<EnergyIncrement>,
    pub mechanical_trials: usize,
    pub physical_plastic_heat_j: f64,
    pub reversible_heat_j: f64,
    pub deposited_mechanical_heat_j: f64,
    pub stored_mechanical_change_j: f64,
    pub thermoelastic_storage_change_j: f64,
    pub numerical_loss_j: f64,
    /// Native nodal dead-force AND prescribed-reaction work, counted once.
    pub mechanical_work_j: f64,
    pub mechanical_residual_j: f64,
    pub closed_radiation_net_j: f64,
    /// Independent authored heat/reservoir/convection input, not internal loss.
    pub external_heat_j: f64,
    /// Closed caloric residual, including physical Joule/plastic/reversible heat.
    pub thermal_residual_j: f64,
    /// Electrical + elastic/hardening + thermoelastic + caloric + numerical
    /// changes, minus electrical supply, mechanical work and true external heat.
    pub combined_residual_j: f64,
}

pub struct ElectricalSolidCandidate<'m, 't> {
    thermal: DepositionSimulation<'t>,
    solid: ThermoplasticSolid<'m>,
    circuit: MaterialCircuit,
    report: ElectricalSolidReport,
}
impl<'m, 't> ElectricalSolidCandidate<'m, 't> {
    pub const fn thermal(&self) -> &DepositionSimulation<'t> { &self.thermal }
    pub const fn solid(&self) -> &ThermoplasticSolid<'m> { &self.solid }
    pub const fn circuit(&self) -> &MaterialCircuit { &self.circuit }
    pub const fn report(&self) -> &ElectricalSolidReport { &self.report }
    pub fn into_parts(self) -> (DepositionSimulation<'t>, ThermoplasticSolid<'m>,
        MaterialCircuit, ElectricalSolidReport) {
        (self.thermal, self.solid, self.circuit, self.report)
    }
}

struct Trial {
    mechanical: MechanicalProposal,
    material: HeatNetworkSample,
    support_contacts: Vec<ContactSample>,
    applied_heat_j: f64,
    external_heat_j: f64,
    radiation_net_j: f64,
}
fn number(value: f64) -> Result<f64, ElectricalSolidSourceError> {
    if value.is_finite() { Ok(value) } else {
        Err(ElectricalSolidSourceError::Invalid("unrepresentable electrical/solid heat or work"))
    }
}

impl<'m> ThermoplasticSolid<'m> {
    /// Resolve support current sources from the accepted mechanical/thermal state.
    /// Use these in `MaterialCircuit::new_with_contacts` or its explicit
    /// `rebind_contacts_candidate` before the first coupled interval. Branches
    /// are checked by that circuit owner; face/source/pressure and heat-partition
    /// admission are shared with the evolving solve below.
    ///
    /// This does not advance time, modify the circuit or infer an ideal switch.
    /// The prescribed DOF list identifies the same fixed whole-face supports
    /// used during evolution; its target values do not author the pressure.
    pub fn resolve_electrical_supports(
        &self,
        thermal: &DepositionSimulation<'_>,
        supports: &[ThermalSupportContact<'_>],
        electrical: &[ElectricalSupportContact],
        prescribed_m: &[(usize, f64)],
        cx: &Cx<'_>,
    ) -> Result<Vec<SeriesElectricalContact>, ElectricalSolidError> {
        if self.time_s != thermal.time_s() {
            return Err(ElectricalSolidError::Invalid(
                "electrical support requires matching accepted clocks",
            ));
        }
        self.enthalpies::<Infallible>(thermal, true, cx)
            .map_err(ElectricalSolidError::Mechanical)?;
        let routes = self
            .support_contact_routes::<Infallible>(thermal, supports, prescribed_m, cx)
            .map_err(|e| match e {
                ContactError::Cancelled => ElectricalSolidError::Cancelled,
                other => ElectricalSolidError::Support(other),
            })?;
        if electrical.len() > supports.len() {
            return Err(ElectricalSolidError::Invalid(
                "electrical support count exceeds physical supports",
            ));
        }
        let mut seen = BTreeSet::new();
        let mut sources = Vec::with_capacity(electrical.len());
        for binding in electrical {
            cx.checkpoint()
                .map_err(|_| ElectricalSolidError::Cancelled)?;
            let index = supports.iter().position(|s| s.link == binding.link).ok_or(
                ElectricalSolidError::Invalid("electrical contact has no solved support"),
            )?;
            if !seen.insert(binding.link) {
                return Err(ElectricalSolidError::Invalid(
                    "duplicate electrical support",
                ));
            }
            let (element, normal) = routes[index];
            let sample = supports[index]
                .sample::<Infallible>(thermal, &self.solution.element_stress[element], normal, 0.0)
                .map_err(ElectricalSolidError::Support)?;
            sources.push(
                binding
                    .resolve(&supports[index], &sample)
                    .map_err(|e| ElectricalSolidError::Coupled(CircuitEnthalpyError::Source(e)))?,
            );
        }
        cx.checkpoint()
            .map_err(|_| ElectricalSolidError::Cancelled)?;
        Ok(sources)
    }

    /// Jointly propose finite electrical storage, loaded native 3D plasticity,
    /// thermal contact/radiation and the SAME caloric field.
    ///
    /// Uses the existing contact-capable circuit/enthalpy Newton solver. Its
    /// callback evaluates this solid's original mechanics and caloric exchange;
    /// no nested thermal solve or post-hoc mechanical refresh occurs. The final
    /// callback's actual equilibrium and history are the returned solid.
    ///
    /// R(T), k(T), emissivity and yield see one endpoint. Joule heat is deposited
    /// only by the circuit owner. The native J2 owner supplies physical plastic
    /// and reversible heat, not hardening or numerical loss. Electrical and heat
    /// routes come from the prepared network, not a separately authored mapping.
    ///
    /// Fixed geometry, C/L, electrical paths and source conditions during this
    /// interval; no Lorentz force, piezoresistivity, melt or circuit topology
    /// change is inferred. Mechanical-coordinate-dependent transport refuses.
    /// Split at events/load discontinuities. Target events remain pending.
    /// Endpoint thermal/mechanical rates and midpoint currents give generally
    /// first-order coupled accuracy; refine intervals independently.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_electrical_candidate<'t>(
        &self, thermal: &DepositionSimulation<'t>, circuit: &MaterialCircuit,
        prepared: &PreparedHeatNetwork<'_, '_>, ports: &[f64], target_s: f64,
        prescribed_m: &[(usize, f64)], nodal_forces_n: &[(usize, f64)],
        settings: IncrementSettings, controls: ElectricalSolidSolve,
        guess_j: &[f64], cx: &Cx<'_>,
    ) -> Result<ElectricalSolidCandidate<'m, 't>, ElectricalSolidError> {
        self.advance_electrical_with_supports_candidate(thermal, circuit, prepared, &[],
            ports, target_s, prescribed_m, nodal_forces_n, settings, controls, guess_j, cx)
    }

    /// Solve finite electrical heating AND load-dependent support heat transfer.
    ///
    /// Each support's pressure is -n.sigma.n from the SAME native mechanical
    /// endpoint, not an independently supplied context value. Its pinned
    /// interface resistance is placed in series with the prepared network's
    /// fresh bulk conductivity. The existing circuit/enthalpy residual then
    /// couples that flux back to resistivity, yield and deformation.
    ///
    /// A support link must not also have a fixed thermal-contact correction or
    /// electrical contact in the prepared network. The latter would freeze an
    /// interface state while its pressure changes; this method does not infer
    /// a pressure-dependent electrical contact law. Other links may retain
    /// their authored electrical/thermal interfaces and radiation.
    ///
    /// Fixed whole-face prescribed supports and compressive traction only.
    /// Geometry, reference mass and electrical topology remain fixed. Empty
    /// support bindings use the original powered calculation unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_electrical_with_supports_candidate<'t>(
        &self, thermal: &DepositionSimulation<'t>, circuit: &MaterialCircuit,
        prepared: &PreparedHeatNetwork<'_, '_>, support_contacts: &[ThermalSupportContact<'_>],
        ports: &[f64], target_s: f64, prescribed_m: &[(usize, f64)], nodal_forces_n: &[(usize, f64)],
        settings: IncrementSettings, controls: ElectricalSolidSolve,
        guess_j: &[f64], cx: &Cx<'_>,
    ) -> Result<ElectricalSolidCandidate<'m, 't>, ElectricalSolidError> {
        self.advance_electrical_with_pressure_contacts_candidate(thermal, circuit, prepared,
            support_contacts, &[], ports, target_s, prescribed_m, nodal_forces_n,
            settings, controls, guess_j, cx)
    }

    /// Close support pressure, contact current, Joule heat and solid response
    /// through the SAME endpoint residual.
    ///
    /// Each electrical binding names one of `support_contacts`; both responses
    /// use its exact ordered interface card, whole-face area, solved compressive
    /// traction and actual two-sided temperatures. The pinned electrical claim
    /// supplies R''(pressure, thermal state), reduced to R = R'' / area by the
    /// original circuit owner. That owner's contact heat enters the two original
    /// caloric cells exactly once using the explicit partition.
    ///
    /// Initialize the accepted circuit contacts with `resolve_electrical_supports`.
    /// Existing non-support contacts retain their prepared sources and heat routes.
    /// The prepared network must omit these solved-pressure interface corrections;
    /// no thermal or electrical resistance may be added twice.
    ///
    /// Fixed full-face compression and original circuit topology only: tensile
    /// traction, unsupported source states and failed coupled balances refuse
    /// atomically. No changing contact area, opening, arc, electromagnetic force,
    /// piezoresistive bulk law or contact-temperature field is inferred.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_electrical_with_pressure_contacts_candidate<'t>(
        &self,
        thermal: &DepositionSimulation<'t>,
        circuit: &MaterialCircuit,
        prepared: &PreparedHeatNetwork<'_, '_>,
        support_contacts: &[ThermalSupportContact<'_>],
        electrical_supports: &[ElectricalSupportContact],
        ports: &[f64],
        target_s: f64,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        settings: IncrementSettings,
        controls: ElectricalSolidSolve,
        guess_j: &[f64],
        cx: &Cx<'_>,
    ) -> Result<ElectricalSolidCandidate<'m, 't>, ElectricalSolidError> {
        self.advance_electrical_with_motion_candidate(thermal, circuit, prepared,
            support_contacts, electrical_supports, ports, target_s, prescribed_m,
            nodal_forces_n, settings, controls, guess_j, None, cx)
            .map(|(candidate, _)| candidate)
    }

    // All electrical, thermal and support-source admission stays shared with
    // the quasistatic path. Only the native mechanical residual gains inertia.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn advance_electrical_with_motion_candidate<'t>(
        &self, thermal: &DepositionSimulation<'t>, circuit: &MaterialCircuit,
        prepared: &PreparedHeatNetwork<'_, '_>, support_contacts: &[ThermalSupportContact<'_>],
        electrical_supports: &[ElectricalSupportContact], ports: &[f64], target_s: f64,
        prescribed_m: &[(usize, f64)], nodal_forces_n: &[(usize, f64)],
        settings: IncrementSettings, controls: ElectricalSolidSolve, guess_j: &[f64],
        velocity_mps: Option<&[[f64; 3]]>, cx: &Cx<'_>,
    ) -> Result<(ElectricalSolidCandidate<'m, 't>, Option<dynamics::Motion>), ElectricalSolidError> {
        let checkpoint = || cx.checkpoint().map_err(|_| DepositionError::Cancelled);
        checkpoint().map_err(|_| ElectricalSolidError::Cancelled)?;
        let dt = target_s - self.time_s;
        let tolerance = controls.coupled.total_energy_tolerance_j;
        if self.time_s != thermal.time_s() || !positive(dt) || !target_s.is_finite()
            || !positive(tolerance) || !positive(controls.yield_tolerance_pa)
            || controls.maximum_mechanical_trials == 0 {
            return Err(ElectricalSolidError::Invalid("matching clocks and bounded positive mechanical controls required"));
        }
        self.admit_heat_network::<Infallible>(thermal, prepared, cx)
            .map_err(ElectricalSolidError::Binding)?;
        self.enthalpies::<Infallible>(thermal, true, cx).map_err(ElectricalSolidError::Mechanical)?;
        let support_routes = self.support_contact_routes::<Infallible>(
            thermal, support_contacts, prescribed_m, cx).map_err(|e| match e {
                ContactError::Cancelled => ElectricalSolidError::Cancelled,
                other => ElectricalSolidError::Support(other),
            })?;
        let destinations = prepared.resistor_destinations();
        let mut contacts = prepared.heat_destinations(thermal).map_err(ElectricalSolidError::Heat)?;
        let routes: BTreeMap<_, _> = destinations.iter().map(|d| (d.branch, d.cell)).collect();
        // Check immutable source dependencies once. PreparedHeatNetwork retains
        // these exact claims, and still re-evaluates their support and values at
        // every endpoint; do not repeat the structural source queries per stencil.
        let initial = prepared.resolve(cx, thermal).map_err(ElectricalSolidError::Heat)?;
        for support in support_contacts {
            if initial.bulk.source.thermal_contacts.iter().any(|c| c.link == support.link)
                || contacts.iter().any(|c| c.id == support.link) {
                return Err(ElectricalSolidError::Invalid(
                    "a solved-pressure support cannot also have a fixed thermal or electrical interface"));
            }
        }
        if electrical_supports.len() > support_contacts.len() {
            return Err(ElectricalSolidError::Invalid(
                "electrical support count exceeds physical supports",
            ));
        }
        for binding in electrical_supports {
            let support = support_contacts
                .iter()
                .find(|s| s.link == binding.link)
                .ok_or(ElectricalSolidError::Invalid(
                    "electrical contact has no solved support",
                ))?;
            let link = &thermal.links()[&support.link];
            if contacts.iter().any(|c| c.id == binding.link)
                || !routes
                    .get(&binding.branch)
                    .is_some_and(|id| *id == link.a || *id == link.b)
            {
                return Err(ElectricalSolidError::Invalid(
                    "electrical support needs a unique contact and a branch of one receiving material cell",
                ));
            }
            contacts.push(contacts::ContactHeatCells {
                id: binding.link,
                cells: [link.a, link.b],
                heat_fraction_to_first: binding.heat_fraction_to_first,
            });
        }
        // Check the same accepted source pins before allocating nonlinear trials.
        let initial_supports = self.resolve_electrical_supports(
            thermal,
            support_contacts,
            electrical_supports,
            prescribed_m,
            cx,
        )?;
        let mut initial_contacts = initial.contacts;
        initial_contacts.extend(initial_supports);
        contacts::admit_sources::<ElectricalSolidSourceError>(circuit, &initial_contacts, || {
            cx.checkpoint().is_err()
        })
        .map_err(ElectricalSolidError::Coupled)?;
        self.admit_electrical_transport(thermal, &initial.bulk.source,
            &initial.bulk.resistors, &routes, cx).map_err(|error|
                ElectricalSolidError::Coupled(CircuitEnthalpyError::Source(error)))?;
        let mechanics = ThermoplasticConfig {
            // Only the energy/yield controls are consumed by the native endpoint
            // evaluator; the shared Newton owner supplies all coupling iteration.
            maximum_iterations: 1, temperature_tolerance_k: f64::EPSILON,
            energy_tolerance_j: tolerance, yield_tolerance_pa: controls.yield_tolerance_pa,
            enthalpy_relaxation: 1.0,
        };
        let trials = Cell::new(0usize);
        let result = contacts::step_candidate(thermal, target_s,
            CircuitEnthalpyInput { circuit, ports, destinations: &destinations }, &contacts,
            controls.coupled, guess_j, |endpoint| {
                checkpoint()?;
                if trials.get() >= controls.maximum_mechanical_trials {
                    return Err(ElectricalSolidSourceError::MechanicalBudget);
                }
                trials.set(trials.get() + 1);
                let ContactThermalCoefficients { bulk: mut coefficients, mut contacts } =
                    prepared.resolve(cx, endpoint).map_err(ElectricalSolidSourceError::Heat)?;
                let h = self.enthalpies::<Infallible>(endpoint, false, cx)
                    .map_err(ElectricalSolidSourceError::Mechanical)?;
                let mechanical = self.mechanics_with_motion::<Infallible>(&h, prescribed_m,
                    nodal_forces_n, settings, mechanics, dt,
                    velocity_mps, cx)
                    .map_err(ElectricalSolidSourceError::Mechanical)?;
                let mut support_samples = Vec::with_capacity(support_contacts.len());
                for (binding, &(element, normal)) in support_contacts.iter().zip(&support_routes) {
                    checkpoint()?;
                    let bulk = *coefficients.conductances_w_k.get(&binding.link)
                        .ok_or(ElectricalSolidSourceError::Invalid("support lacks original bulk conductance"))?;
                    let sample = binding.sample::<Infallible>(endpoint,
                        &mechanical.solution.element_stress[element], normal, dt)
                        .and_then(|s| s.with_bulk(bulk, dt))
                        .map_err(ElectricalSolidSourceError::Support)?;
                    coefficients.conductances_w_k.insert(binding.link, sample.conductance_w_k);
                    support_samples.push(sample);
                }
                for binding in electrical_supports {
                    checkpoint()?;
                    let support = support_contacts
                        .iter()
                        .find(|s| s.link == binding.link)
                        .ok_or(ElectricalSolidSourceError::Invalid(
                            "electrical support disappeared",
                        ))?;
                    let sample = support_samples
                        .iter()
                        .find(|s| s.link == binding.link)
                        .ok_or(ElectricalSolidSourceError::Invalid(
                            "solved support pressure is absent",
                        ))?;
                    contacts.push(binding.resolve(support, sample)?);
                }
                let radiation_net_j = number(dt * coefficients.source.radiation_net_power_w())?;
                let external_heat_j = number(number(dt * endpoint.external_step_power_w(
                    &coefficients.sources_w, &coefficients.conductances_w_k, checkpoint)?)? - radiation_net_j)?;
                let mut applied_heat_j = 0.0;
                for (bound, heat) in self.cells.iter().zip(&mechanical.heat_j) {
                    checkpoint()?;
                    let watts = number(*heat / dt)?;
                    let previous = coefficients.sources_w.get(&bound.accepted.id).copied().unwrap_or(0.0);
                    let combined = number(previous + watts)?;
                    if *heat != 0.0 && (watts == 0.0 || combined == previous) {
                        return Err(ElectricalSolidSourceError::Invalid("mechanical heat disappeared from its receiving cell"));
                    }
                    if watts != 0.0 { coefficients.sources_w.insert(bound.accepted.id, combined); }
                    applied_heat_j = number(applied_heat_j + number(dt * watts)?)?;
                }
                let heat_error = number(applied_heat_j - mechanical.physical_heat_j - mechanical.reversible_heat_j)?;
                if heat_error.abs() > tolerance {
                    return Err(ElectricalSolidSourceError::Invalid("native caloric heat and deposited mechanical heat disagree"));
                }
                checkpoint()?;
                Ok(ContactThermalCoefficients { contacts, bulk: CircuitThermalCoefficients {
                    resistors: coefficients.resistors, sources_w: coefficients.sources_w,
                    conductances_w_k: coefficients.conductances_w_k,
                    source: Trial { mechanical, material: coefficients.source,
                        support_contacts: support_samples, applied_heat_j,
                        external_heat_j, radiation_net_j },
                } })
            }, || cx.checkpoint().is_err());
        // Cancellation in a nested owner is still cancellation, not a usable
        // numerical failure or a partial circuit/solid result.
        checkpoint().map_err(|_| ElectricalSolidError::Cancelled)?;
        let candidate = result.map_err(ElectricalSolidError::Coupled)?;
        let (thermal, circuit, r) = candidate.into_parts();
        let trial = r.source.material;
        let mechanical = trial.mechanical;
        let thermal_residual_j = r.thermal.stored_energy_change_j - r.thermal.incoming_energy_j
            - trial.external_heat_j - trial.applied_heat_j - r.deposited_joule_j;
        let combined_residual_j = r.electrical.stored_after_j - r.electrical.stored_before_j
            + mechanical.stored_change_j + mechanical.thermoelastic_storage_change_j
            + r.thermal.stored_energy_change_j + mechanical.numerical_loss_j
            + mechanical.motion.as_ref().map_or(0.0, |motion|
                motion.energy.kinetic_change_j + motion.energy.numerical_loss_j)
            - r.electrical.native.supplied - mechanical.external_work_j
            - trial.external_heat_j - r.thermal.incoming_energy_j;
        let source_residual_j = r.thermal.external_heat_j - trial.external_heat_j
            - trial.radiation_net_j - trial.applied_heat_j - r.deposited_joule_j;
        if [thermal_residual_j, combined_residual_j, source_residual_j].iter()
            .any(|r| !r.is_finite() || r.abs() > tolerance) {
            return Err(ElectricalSolidError::Invalid("closed electrical, native solid and caloric balances do not agree"));
        }
        let mut solid = self.clone();
        solid.time_s = target_s;
        solid.materials = mechanical.materials;
        solid.solution = mechanical.solution;
        for cell in &mut solid.cells { cell.accepted = thermal.cells()[&cell.accepted.id].clone(); }
        let report = ElectricalSolidReport {
            support_contacts: trial.support_contacts,
            circuit_thermal: CircuitEnthalpyReport {
                electrical: r.electrical, thermal: r.thermal,
                source: ContactThermalSample { material: trial.material, destinations: r.source.destinations },
                nonlinear: r.nonlinear, maximum_cell_residual_j: r.maximum_cell_residual_j,
                circuit_trials: r.circuit_trials, deposited_joule_j: r.deposited_joule_j,
                non_electrical_heat_j: r.non_electrical_heat_j, total_energy_residual_j: r.total_energy_residual_j,
            },
            element_energy: mechanical.energy, mechanical_trials: trials.get(),
            physical_plastic_heat_j: mechanical.physical_heat_j, reversible_heat_j: mechanical.reversible_heat_j,
            deposited_mechanical_heat_j: trial.applied_heat_j, stored_mechanical_change_j: mechanical.stored_change_j,
            thermoelastic_storage_change_j: mechanical.thermoelastic_storage_change_j,
            numerical_loss_j: mechanical.numerical_loss_j, mechanical_work_j: mechanical.external_work_j,
            mechanical_residual_j: mechanical.mechanical_residual_j,
            closed_radiation_net_j: trial.radiation_net_j, external_heat_j: trial.external_heat_j,
            thermal_residual_j, combined_residual_j,
        };
        checkpoint().map_err(|_| ElectricalSolidError::Cancelled)?;
        Ok((ElectricalSolidCandidate { thermal, solid, circuit, report }, mechanical.motion))
    }

    // The heat network already binds exact fixed conditions. Reuse the original
    // material guard as well: k/rho_e/emissivity cannot freeze a strain, stress or
    // frequency coordinate while the native solid actually changes that state.
    fn admit_electrical_transport(
        &self, thermal: &DepositionSimulation<'_>, sample: &HeatNetworkSample,
        resistors: &[MaterialResistor], routes: &BTreeMap<usize, CellId>, cx: &Cx<'_>,
    ) -> Result<(), ElectricalSolidSourceError> {
        let admit = |id: CellId, state: &ResolvedMaterialStatePoint| -> Result<(), ElectricalSolidSourceError> {
            let Some(bound) = self.cells.iter().find(|b| b.accepted.id == id) else { return Ok(()); };
            let material = bound.source.thermal_material_point().material;
            let h = thermal.cells()[&id].specific_h_j_kg;
            for property in state.properties() {
                poll::<Infallible>(cx).map_err(ElectricalSolidSourceError::Mechanical)?;
                let selected = property.answer().receipt.selected;
                let claim = material.claims().claim(selected)
                    .ok_or(ElectricalSolidSourceError::Invalid("transport claim is not part of this structural material"))?;
                bound.source.pin_thermal_transport(h, &claim.key,
                    MaterialPropertySelection::PinnedByProperty(vec![(claim.key.name().to_owned(), selected)]))
                    .map_err(|error| ElectricalSolidSourceError::Mechanical(ThermoplasticError::Material(error)))?;
            }
            Ok(())
        };
        for cell in &sample.material.cells {
            if let Some(state) = &cell.conductivity { admit(cell.cell, state)?; }
        }
        if let Some(radiation) = &sample.radiation {
            for surface in radiation.surfaces() { admit(surface.cell, &surface.material)?; }
        }
        for resistor in resistors {
            let cell = *routes.get(&resistor.branch)
                .ok_or(ElectricalSolidSourceError::Invalid("material resistor has no original receiving cell"))?;
            admit(cell, resistor.conductor.material())?;
        }
        Ok(())
    }
}
