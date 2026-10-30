//! Finite atmosphere heat exchange on the original material/phase/flow timeline.
//!
//! Condensed cells retain their original total enthalpy (including latent heat).
//! The closed fixed-volume gas retains NASA-9 INTERNAL energy, not cp*T or h.
//! A prescribed stationary film exchanges equal/opposite joules between them.
//! The resulting gas temperature and EOS pressure re-evaluate the same source
//! transport, while condensed h controls phase, conductivity and native momentum.

use super::*;
use fs_conduction::transient::deposition::{EnergyAudit, StepConfig};
use fs_material::thermochemical_gas::ResolvedThermochemicalGasState;
use fs_material::thermochemical_gas::energy::{GasEnergyError, GasEnergySolve, IsochoricGasParcel};

/// Caller-resolved stationary exposed patch. One entry per receiving cell.
/// Area cannot exceed that cell's exposed area after existing contacts/bed.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteGasContact {
    pub cell: CellId,
    pub area_m2: f64,
    /// Condensed centroid-to-interface distance and gas-film thickness [m].
    /// These are physical inputs, never inferred from a material or gas name.
    pub distances_m: [f64; 2],
}

/// Source state and coefficients actually used on one exposed patch.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteGasContactSample {
    pub contact: FiniteGasContact,
    pub specific_enthalpy_j_kg: f64,
    pub condensed_temperature_k: f64,
    pub condensed_conductivity_w_m_k: f64,
    pub gas_temperature_k: f64,
    /// Mean-temperature film at the finite parcel's actual EOS pressure.
    pub film: ResolvedThermochemicalGasState,
    pub conductance_w_k: f64,
    /// Signed watts entering the condensed cell and leaving the gas.
    pub heat_into_condensed_w: f64,
}

/// Accepted exchange ledger. Gas storage belongs to a separate internal-energy
/// owner; its heat is included in `thermal.external_heat_j` only because that
/// ledger describes the condensed subsystem. It is internal to their union.
#[derive(Clone, Debug)]
pub struct FiniteGasExchangeReport {
    pub thermal: EnergyAudit,
    pub heat_into_condensed_j: f64,
    pub gas_internal_energy_change_j: f64,
    /// Actual gas storage change plus its opposite supplied condensed heat.
    pub gas_balance_defect_j: f64,
    /// Actual combined storage minus incoming enthalpy and all other heat.
    /// Optional mechanical work/kinetic accounts still belong to `flow`.
    pub combined_thermal_balance_defect_j: f64,
    /// Thermal external heat with this internal exchange removed. Existing
    /// molecular/Joule source accounts retain their original interpretation.
    pub other_thermal_external_heat_j: f64,
    pub substeps: usize,
    pub last_applied: Vec<FiniteGasContactSample>,
    pub endpoint: Vec<FiniteGasContactSample>,
}

#[derive(Clone, Debug, PartialEq)]
pub enum FiniteGasFlowError<E> {
    Flow(MaterialFlowError<E>),
    Gas(GasEnergyError),
}
impl<E> From<DepositionError> for FiniteGasFlowError<E> {
    fn from(error: DepositionError) -> Self {
        Self::Flow(error.into())
    }
}
impl<E> From<MaterialFlowError<E>> for FiniteGasFlowError<E> {
    fn from(error: MaterialFlowError<E>) -> Self {
        Self::Flow(error)
    }
}
impl<E> From<GasEnergyError> for FiniteGasFlowError<E> {
    fn from(error: GasEnergyError) -> Self {
        Self::Gas(error)
    }
}
impl<E: core::fmt::Display> core::fmt::Display for FiniteGasFlowError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Flow(error) => write!(f, "finite gas material flow: {error}"),
            Self::Gas(error) => write!(f, "finite gas material flow: {error}"),
        }
    }
}
impl<E: std::error::Error + 'static> std::error::Error for FiniteGasFlowError<E> {}

fn finite(value: f64) -> Result<f64, DepositionError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(DepositionError::Arithmetic)
    }
}

fn same_sources(a: &ResolvedThermochemicalGasState, b: &ResolvedThermochemicalGasState) -> bool {
    a.components().len() == b.components().len()
        && a.components().iter().zip(b.components()).all(|(a, b)| {
            a.card_identity() == b.card_identity()
                && a.properties().len() == b.properties().len()
                && a.properties().iter().zip(b.properties()).all(|(a, b)| {
                    a.requirement() == b.requirement()
                        && a.answer().receipt.selected == b.answer().receipt.selected
                })
        })
}

#[derive(Clone)]
struct GasFlow<'a, 'f, F, R, O, T> {
    flow: MaterialFlowFeedback<'f, F, R, O, T>,
    gas: IsochoricGasParcel,
    contacts: &'a [FiniteGasContact],
    controls: GasEnergySolve,
    gas_temperature_slope: f64,
    reference: Vec<FiniteGasContactSample>,
    applied: Vec<FiniteGasContactSample>,
    endpoint: Vec<FiniteGasContactSample>,
    applied_time_s: Option<f64>,
    heat_into_condensed_j: f64,
    substeps: usize,
}

impl<F: Clone, R, O: DepositionObserver, T: FlowTransportSnapshot> GasFlow<'_, '_, F, R, O, T>
where
    R: Clone + FnMut(usize, &DepositedCell, f64) -> Result<F, MaterialError>,
{
    fn sample(
        &self,
        thermal: &DepositionSimulation<'_>,
        samples: &[FlowMaterialSample<F, T>],
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Vec<FiniteGasContactSample>, FiniteGasFlowError<O::Error>> {
        let mut result = Vec::with_capacity(self.contacts.len());
        let by_cell: BTreeMap<_, _> = samples.iter().map(|sample| (sample.cell, sample)).collect();
        let gas_t = self.gas.state().state().temperature;
        if gas_t < self.controls.minimum_temperature_k
            || gas_t > self.controls.maximum_temperature_k
        {
            return Err(invalid("finite gas is outside its admitted caloric bracket").into());
        }
        for (index, contact) in self.contacts.iter().enumerate() {
            check_gate(self.flow.gate, check)?;
            let sample = by_cell.get(&contact.cell).ok_or(invalid(
                "finite gas patch needs a mapped material-flow cell",
            ))?;
            let cell = thermal
                .cells()
                .get(&contact.cell)
                .ok_or(invalid("finite gas patch has no active thermal receiver"))?;
            if sample.specific_enthalpy_j_kg.to_bits() != cell.specific_h_j_kg.to_bits()
                || sample.temperature_k.to_bits() != thermal.temperature_k(contact.cell)?.to_bits()
                || contact.area_m2 > thermal.exposed_area_m2(contact.cell)?
                || cell.boundary.htc_w_m2_k != 0.0
            {
                return Err(invalid("finite gas needs the actual exposed source epoch without duplicate ambient convection").into());
            }
            let k = self.flow.effective_conductivity(sample)?;
            let film_t = 0.5 * gas_t + 0.5 * sample.temperature_k;
            let point = self
                .gas
                .point()
                .clone()
                .with("temperature", film_t)
                .map_err(|_| invalid("unrepresentable finite gas film coordinates"))?;
            let film = self
                .gas
                .model()
                .evaluate(&point)
                .map_err(GasEnergyError::from)?;
            if self
                .reference
                .get(index)
                .is_some_and(|start| !same_sources(&start.film, &film))
            {
                return Err(
                    invalid("finite gas film left its original transport source claims").into(),
                );
            }
            let conductance_w_k = series_conductance_w_k(
                contact.area_m2,
                contact.distances_m,
                [k, film.state().thermal_conductivity],
            )?;
            let heat_into_condensed_w = finite(conductance_w_k * (gas_t - sample.temperature_k))?;
            result.push(FiniteGasContactSample {
                contact: *contact,
                specific_enthalpy_j_kg: sample.specific_enthalpy_j_kg,
                condensed_temperature_k: sample.temperature_k,
                condensed_conductivity_w_m_k: k,
                gas_temperature_k: gas_t,
                film,
                conductance_w_k,
                heat_into_condensed_w,
            });
        }
        Ok(result)
    }
}

impl<F: Clone, R, O: DepositionObserver, T: FlowTransportSnapshot> DepositionObserver
    for GasFlow<'_, '_, F, R, O, T>
where
    R: Clone + FnMut(usize, &DepositedCell, f64) -> Result<F, MaterialError>,
{
    type Error = FiniteGasFlowError<O::Error>;

    fn begin(&mut self) -> Result<(), Self::Error> {
        self.flow.begin()?;
        self.applied_time_s = None;
        Ok(())
    }

    fn requires_step_interval(&self) -> bool {
        self.flow.requires_step_interval()
    }

    fn sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        let mut sources = self.flow.sources(thermal, check)?;
        let epoch = self
            .flow
            .applied()
            .ok_or(invalid("finite gas lacks an applied flow epoch"))?;
        if epoch.time_s.to_bits() != thermal.time_s().to_bits() {
            return Err(invalid("stale finite gas flow epoch").into());
        }
        self.applied = self.sample(thermal, &epoch.samples, check)?;
        for sample in &self.applied {
            if sample.heat_into_condensed_w != 0.0 {
                let value = sources.entry(sample.contact.cell).or_insert(0.0);
                *value = finite(*value + sample.heat_into_condensed_w)?;
            }
        }
        self.applied_time_s = Some(thermal.time_s());
        Ok(sources)
    }

    fn conductances(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<u64, f64>, Self::Error> {
        self.flow.conductances(thermal, check).map_err(Into::into)
    }

    fn reservoir_conductances(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.flow.reservoir_conductances(thermal, check).map_err(Into::into)
    }

    fn heat_transfers(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Vec<HeatTransfer>, Self::Error> {
        self.flow.heat_transfers(thermal, check).map_err(Into::into)
    }

    fn limit_timestep_s(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        proposed_s: f64,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<f64, Self::Error> {
        let dt = self.flow.limit_timestep_s(thermal, proposed_s, check)?;
        if !dt.is_finite()
            || dt <= 0.0
            || dt > proposed_s
            || self.applied_time_s != Some(thermal.time_s())
        {
            return Err(invalid("invalid finite gas timestep/source epoch").into());
        }
        let slope = self.flow.input.caloric.maximum_temperature_slope();
        let mut gas_g = 0.0;
        let mut decay: f64 = 0.0;
        for sample in &self.applied {
            check_gate(self.flow.gate, check)?;
            gas_g = finite(gas_g + sample.conductance_w_k)?;
            decay = decay.max(finite(
                slope * sample.conductance_w_k / thermal.cells()[&sample.contact.cell].mass_kg,
            )?);
        }
        decay = decay.max(finite(
            self.gas_temperature_slope * gas_g / self.gas.mass_kg(),
        )?);
        if decay == 0.0 {
            return Ok(dt);
        }
        // Add the new row decay to the existing conduction/advection limit.
        // Independent minima would allow the summed operator to overshoot.
        let film_dt = 0.8 / decay;
        let limited = if dt <= film_dt {
            dt / (1.0 + dt / film_dt)
        } else {
            film_dt / (1.0 + film_dt / dt)
        };
        if !limited.is_finite() || limited <= 0.0 {
            return Err(DepositionError::Arithmetic.into());
        }
        Ok(limited)
    }

    fn step_sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        start_s: f64,
        dt_s: f64,
        end_s: f64,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.flow
            .step_sources(thermal, start_s, dt_s, end_s, check)
            .map_err(Into::into)
    }

    fn event(
        &mut self,
        event: &Event,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<(), Self::Error> {
        self.flow.event(event, thermal)?;
        self.endpoint = self.sample(thermal, self.flow.endpoint(), &mut || Ok(()))?;
        Ok(())
    }

    fn step(
        &mut self,
        start_s: f64,
        dt: f64,
        old: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        if self.applied_time_s != Some(start_s) {
            return Err(invalid("finite gas has no matching accepted heat interval").into());
        }
        let mut watts = 0.0;
        for sample in &self.applied {
            watts = finite(watts + sample.heat_into_condensed_w)?;
        }
        let heat = finite(dt * watts)?;
        let old_gas = self.gas.state().clone();
        let gate = self.flow.gate;
        let mut checkpoint_error = None;
        let result = self.gas.add_heat(-heat, self.controls, || {
            if let Err(error) = check_gate(gate, check) {
                checkpoint_error = Some(error);
                true
            } else {
                false
            }
        });
        if let Some(error) = checkpoint_error {
            return Err(error.into());
        }
        result?;
        if !same_sources(&old_gas, self.gas.state()) {
            return Err(
                invalid("finite gas storage left its original transport source claims").into(),
            );
        }
        self.flow.step(start_s, dt, old, thermal, check)?;
        self.endpoint = self.sample(thermal, self.flow.endpoint(), check)?;
        self.heat_into_condensed_j = finite(self.heat_into_condensed_j + heat)?;
        self.substeps += 1;
        self.applied_time_s = None;
        Ok(())
    }
}

/// Advance a finite atmosphere and the ORIGINAL material/phase/MAC timeline.
///
/// Source k at each condensed h/phase and gas-film mean temperature forms the
/// two-layer Fourier conductance. Its gas-side state uses the parcel's actual
/// evolving pressure. One frozen source epoch drives each explicit substep;
/// combined gas/condensed/advection decay restricts dt. The accepted opposite
/// gas heat uses the existing nonlinear NASA-9 internal-energy inversion.
/// Latent heat, exact solid occupancy, liquid viscosity, declared Darcy drag,
/// finite momentum and optional molecular heating remain with `flow`.
///
/// The gas is well mixed with frozen composition, fixed mass and volume. This
/// needs a well-mixed bulk / quasisteady film justification and time refinement;
/// there is no empirical convection correlation or inferred film thickness.
/// Interface area and all reference volumes remain fixed. Gas pressure is a
/// thermodynamic observation here, not an imposed MAC traction or piston work.
/// No evaporation, species exchange, moving solid skeleton, volume change,
/// reaction or radiation is inferred. Existing moving-freeze refusals remain.
///
/// Patches must use exposed mapped cells with no existing ambient convection;
/// bed/internal links still retain their own surface coverage. Actual source
/// support is checked at every film/parcel and final endpoint. The returned
/// thermal/flow/gas states are one candidate: any refusal or cancellation leaves
/// ALL inputs untouched, including late gas-domain and phase-momentum failures.
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
pub fn advance_finite_gas_with_material_flow<'t, 'f, F, R, O, T>(
    thermal: &DepositionSimulation<'t>,
    target_s: f64,
    config: StepConfig,
    flow: &MaterialFlowFeedback<'f, F, R, O, T>,
    gas: &IsochoricGasParcel,
    contacts: &[FiniteGasContact],
    controls: GasEnergySolve,
) -> Result<
    (
        DepositionSimulation<'t>,
        MaterialFlowFeedback<'f, F, R, O, T>,
        IsochoricGasParcel,
        FiniteGasExchangeReport,
    ),
    FiniteGasFlowError<O::Error>,
>
where
    F: Clone,
    T: FlowTransportSnapshot,
    R: Clone + FnMut(usize, &DepositedCell, f64) -> Result<F, MaterialError>,
    O: DepositionObserver,
{
    check_gate(flow.gate, &mut || Ok(()))?;
    flow.check_map(thermal, &mut || Ok(()))?;
    if contacts.is_empty() || contacts.len() > flow.input.max_connections {
        return Err(invalid("finite gas needs bounded exposed contacts").into());
    }
    let mut seen = std::collections::BTreeSet::new();
    for contact in contacts {
        check_gate(flow.gate, &mut || Ok(()))?;
        if !seen.insert(contact.cell)
            || !contact.area_m2.is_finite()
            || contact.area_m2 < 0.0
            || contact
                .distances_m
                .iter()
                .any(|d| !d.is_finite() || *d <= 0.0)
        {
            return Err(invalid(
                "finite gas needs unique cells, finite area and positive film distances",
            )
            .into());
        }
    }
    let slope = gas.model().maximum_isochoric_temperature_slope(controls)?;
    let mut observer = GasFlow {
        flow: flow.clone(),
        gas: gas.clone(),
        contacts,
        controls,
        gas_temperature_slope: slope,
        reference: Vec::new(),
        applied: Vec::new(),
        endpoint: Vec::new(),
        applied_time_s: None,
        heat_into_condensed_j: 0.0,
        substeps: 0,
    };
    observer.reference = observer.sample(thermal, flow.endpoint(), &mut || Ok(()))?;
    observer.endpoint = observer.reference.clone();
    let (next, accepted, audit) = thermal.advance_observed(&observer, target_s, config, || {
        check_gate(flow.gate, &mut || Ok(()))
    })?;
    let stored = finite(
        gas.mass_kg()
            * (accepted.gas.state().specific_internal_energy_j_kg()
                - gas.state().specific_internal_energy_j_kg()),
    )?;
    let other_heat = finite(audit.external_heat_j - accepted.heat_into_condensed_j)?;
    let gas_defect = finite(stored + accepted.heat_into_condensed_j)?;
    let gas_tolerance =
        finite(accepted.substeps as f64 * gas.mass_kg() * controls.energy_tolerance_j_kg)?;
    if gas_defect.abs() > gas_tolerance {
        return Err(GasEnergyError::Invalid(
            "accumulated finite gas exchange exceeded its absolute energy budget",
        )
        .into());
    }
    let combined_defect =
        finite(audit.stored_energy_change_j + stored - audit.incoming_energy_j - other_heat)?;
    if combined_defect.abs() > finite(config.energy_tolerance_j + gas_tolerance)? {
        return Err(DepositionError::EnergyBalance(combined_defect).into());
    }
    let report = FiniteGasExchangeReport {
        other_thermal_external_heat_j: other_heat,
        thermal: audit,
        heat_into_condensed_j: accepted.heat_into_condensed_j,
        gas_internal_energy_change_j: stored,
        gas_balance_defect_j: gas_defect,
        combined_thermal_balance_defect_j: combined_defect,
        substeps: accepted.substeps,
        last_applied: accepted.applied,
        endpoint: accepted.endpoint,
    };
    check_gate(flow.gate, &mut || Ok(()))?;
    Ok((next, accepted.flow, accepted.gas, report))
}
