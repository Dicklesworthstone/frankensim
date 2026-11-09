//! Coherent Kelvin/Maxwell heat enters the original finite-cell thermal solve.
//!
//! The original solid observer integrates every material point over the actual
//! thermal interval once, including carrier and beat terms. Only its physical
//! dashpot heat becomes a caloric source. Native nonlinear conduction, contacts,
//! reservoirs and radiation still resolve at every proposed thermal endpoint.
//! The joint account encloses constitutive storage and caloric storage; supplied
//! stress work is not an independently established mechanical force-work balance.

use super::{ElementThermalCell, Receiver, SolidThermalError};
use crate::circuit_enthalpy::{
    CircuitEnthalpyError, CircuitEnthalpySolve, DrivenCoefficients, DriverAccount, ThermalDriver,
    admit_driver_unknowns, cell_driver_candidate,
    contacts::heat::{HeatAssembly, HeatNetworkError, closed_balance},
    material::MaterialCell,
};
use fs_exec::Cx;
use fs_fft::harmonic_window::HarmonicWindow;
use fs_material::{
    thermal::source::ThermalMaterialError,
    visco_tensor::source::{ResolvedTensorMaterial, TensorSourceError},
};
use fs_math::{c64::C64, dd::Dd};
use fs_modal::harmonic::HarmonicResponseError;
use fs_solid::linear3::{
    TetElasticError,
    visco::{
        ViscoElementError, ViscoHarmonicError, ViscoHarmonicMaterialSpectrumField,
        ViscoHarmonicMaterialSpectrumObserver,
    },
};
use fs_solver::NewtonReport;
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation, EnergyAudit};
use std::{
    collections::{BTreeMap, BTreeSet},
    convert::Infallible,
};

/// Original source, quadrature, receiving-cell and heat-provider failures.
#[derive(Debug)]
pub enum HarmonicSolidHeatError<E> {
    Observation(ViscoHarmonicError),
    Source(TensorSourceError),
    Material(ThermalMaterialError),
    Receiver(SolidThermalError<Infallible>),
    Heat(E),
}
impl<E: core::fmt::Debug> core::fmt::Display for HarmonicSolidHeatError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "coherent tensor-solid heat: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for HarmonicSolidHeatError<E> {}
pub type Error<E> = CircuitEnthalpyError<HarmonicSolidHeatError<E>>;

fn poll<E>(cx: &Cx<'_>) -> Result<(), Error<E>> {
    cx.checkpoint().map_err(|_| CircuitEnthalpyError::Cancelled)
}
fn heat<E: From<HeatNetworkError>>(error: HeatNetworkError) -> Error<E> {
    match error {
        HeatNetworkError::Cancelled => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(HarmonicSolidHeatError::Heat(other.into())),
    }
}
fn receiver<E>(error: SolidThermalError<Infallible>) -> Error<E> {
    match error {
        SolidThermalError::Cancelled => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(HarmonicSolidHeatError::Receiver(other)),
    }
}
fn observation<E>(error: ViscoHarmonicError) -> Error<E> {
    match error {
        ViscoHarmonicError::Element(ViscoElementError::Elastic(TetElasticError::Cancelled))
        | ViscoHarmonicError::Harmonic(HarmonicResponseError::Cancelled) => {
            CircuitEnthalpyError::Cancelled
        }
        other => CircuitEnthalpyError::Source(HarmonicSolidHeatError::Observation(other)),
    }
}
fn same_inventory(a: f64, b: f64) -> bool {
    a.is_finite()
        && b.is_finite()
        && a > 0.0
        && b > 0.0
        && (a - b).abs() <= 64.0 * f64::EPSILON * a.max(b)
}
fn finite(value: f64) -> Result<f64, DepositionError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(DepositionError::Arithmetic)
    }
}
fn balance<E>(residual: f64, tolerance: f64) -> Result<(), Error<E>> {
    if finite(residual)?.abs() > tolerance {
        Err(CircuitEnthalpyError::Invalid(
            "harmonic constitutive heat, its caloric handoff or the joint energy account failed",
        ))
    } else {
        Ok(())
    }
}

/// Accepted original quadrature and finite-cell thermal accounting.
pub struct HarmonicSolidHeatReport<S> {
    pub material: ViscoHarmonicMaterialSpectrumField,
    /// In the original mesh-element order, each Kelvin/Maxwell loss appears once.
    pub element_heat_j: Vec<f64>,
    /// Native thermal audit includes material heat in external_heat_j.
    pub thermal: EnergyAudit,
    pub source: S,
    pub nonlinear: NewtonReport,
    pub maximum_cell_residual_j: f64,
    pub thermal_trials: usize,
    pub deposited_heat_j: f64,
    pub heat_handoff_residual_j: f64,
    /// Excludes the internal material heat and numerical closed-radiation net.
    pub other_external_heat_j: f64,
    /// Thermal storage change minus material heat minus other external heat.
    pub caloric_energy_residual_j: f64,
    /// Constitutive plus caloric storage change minus stress work and external
    /// thermal heat. This does not include kinetic energy or applied-force work.
    pub constitutive_thermal_energy_residual_j: f64,
}

/// Unpublished caloric endpoint; the immutable harmonic source never advances.
pub struct HarmonicSolidHeatCandidate<'t, S> {
    thermal: DepositionSimulation<'t>,
    report: HarmonicSolidHeatReport<S>,
}
impl<'t, S> HarmonicSolidHeatCandidate<'t, S> {
    pub fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub fn report(&self) -> &HarmonicSolidHeatReport<S> {
        &self.report
    }
    /// The returned thermal state already includes every reported material joule.
    pub fn into_parts(self) -> (DepositionSimulation<'t>, HarmonicSolidHeatReport<S>) {
        (self.thermal, self.report)
    }
}

// Share the native cell residual and final fresh acceptance without introducing
// a second thermal integrator or cloning the quadrature field for Newton trials.
struct ConstitutiveTransfer<'a>(&'a ViscoHarmonicMaterialSpectrumField);
impl ThermalDriver for ConstitutiveTransfer<'_> {
    type Model = ();
    type Step = ();
    fn account(&self) -> DriverAccount<'_> {
        DriverAccount {
            resistors: &[],
            contacts: &[],
            heat_j: self.0.material_heat_j,
            numerical_dissipation_j: 0.0,
            stored_change_j: self.0.stored_energy_change_j,
            supplied_j: self.0.stress_work_j,
        }
    }
    fn into_parts(self) -> ((), ()) {
        ((), ())
    }
}

// Source axes, rather than equality at two samples, establish independence of
// the evolving temperature. The canonical refresh then admits exact selections,
// all held coordinates and the complete requested interval of applicability.
fn admit_source<E>(
    source: &ResolvedTensorMaterial,
    binding: &MaterialCell<'_>,
    interval: [f64; 2],
    cx: &Cx<'_>,
) -> Result<(), Error<E>> {
    let axis = binding.point.temperature_axis;
    let state = source.state();
    if state.card_identity() != binding.point.material.content_hash()
        || state.axis_quantities() != binding.point.context.axis_quantities()
        || !state
            .query_point()
            .iter()
            .filter(|(a, _)| a != axis)
            .map(|(a, v)| (a, *v))
            .eq(binding
                .point
                .context
                .axes()
                .iter()
                .filter(|(a, _)| a.as_str() != axis)
                .map(|(a, v)| (a, *v)))
    {
        return Err(CircuitEnthalpyError::Invalid(
            "harmonic heat requires the original material card and complete nonthermal query",
        ));
    }
    for state in core::iter::once(state).chain(source.kelvin_source().map(|k| k.state())) {
        for property in state.properties() {
            poll(cx)?;
            if !property.has_only_fixed_external_conditions(Some(axis), true)
                || property.source_axes().iter().any(|(a, _)| a == axis)
            {
                return Err(CircuitEnthalpyError::Invalid(
                    "settled harmonic heating requires temperature-independent tensor, relaxation, Kelvin and density source laws",
                ));
            }
        }
    }
    for temperature in interval {
        poll(cx)?;
        let point = state.query().with(axis, temperature).map_err(|_| {
            CircuitEnthalpyError::Invalid("invalid harmonic source temperature interval")
        })?;
        let refreshed = source
            .at_query_point(binding.point.material, &point)
            .map_err(|e| CircuitEnthalpyError::Source(HarmonicSolidHeatError::Source(e)))?;
        if refreshed.law() != source.law()
            || refreshed.density_kg_m3().to_bits() != source.density_kg_m3().to_bits()
        {
            return Err(CircuitEnthalpyError::Invalid(
                "harmonic material coefficients or reference density changed with temperature",
            ));
        }
    }
    Ok(())
}

/// Deposit exact coherent material heat over `[thermal.time_s(), target_s]`.
///
/// Each original mesh element needs one unique whole-volume, whole-reference-
/// mass receiver in the prepared heat network. Its exact tensor and optional
/// Kelvin sources must be independent of the evolving temperature throughout
/// the explicit fully-solid interval. The original laws still admit every
/// frequency and the complete coherent strain waveform. Additional passive
/// thermal cells, nonlinear transport, contacts and radiation remain active.
///
/// The observer runs once, outside Newton. Exact interval heat is an interval
/// mean source on the existing backward-Euler thermal residual; refine thermal
/// durations to resolve beat-driven transport. Recovered spring energy, stress
/// work and numerical residual are never deposited. No startup history,
/// temperature-dependent mechanics, thermal expansion or dynamic force balance
/// is inferred. Process start events first and split at source/topology changes.
/// Cancellation or a late source refusal changes neither accepted owner.
#[allow(clippy::too_many_arguments)]
pub fn heat_assembly_candidate<'t, A: HeatAssembly + ?Sized>(
    provider: &A,
    cx: &Cx<'_>,
    thermal: &DepositionSimulation<'t>,
    target_s: f64,
    observer: &ViscoHarmonicMaterialSpectrumObserver<'_>,
    displacements: &[&[C64]],
    cells: &[ElementThermalCell<'_>],
    controls: CircuitEnthalpySolve,
    initial_energy_increments_j: &[f64],
) -> Result<HarmonicSolidHeatCandidate<'t, A::Sample>, Error<A::Error>> {
    poll(cx)?;
    admit_driver_unknowns(thermal, controls, initial_energy_increments_j)?;
    let dt = finite(target_s - thermal.time_s())?;
    let window = HarmonicWindow::new(thermal.time_s(), dt).map_err(|_| {
        CircuitEnthalpyError::Invalid(
            "harmonic heat needs a positive representable thermal interval",
        )
    })?;
    let solid = observer.source();
    let sources = solid
        .material_sources()
        .ok_or(CircuitEnthalpyError::Invalid(
            "harmonic self-heating requires original tensor and density source bindings",
        ))?;
    let network = provider
        .prepared_heat_network()
        .thermal_material_network()
        .map_err(heat)?;
    if cells.len() != solid.element_count()
        || sources.len() != cells.len()
        || network.cells.len() != thermal.cells().len()
    {
        return Err(CircuitEnthalpyError::Invalid(
            "complete original element and caloric bindings required",
        ));
    }
    let network_cells: BTreeMap<_, _> = network.cells.iter().map(|b| (b.cell, b)).collect();
    if network_cells.len() != network.cells.len() {
        return Err(CircuitEnthalpyError::Invalid(
            "duplicate heat-network cell binding",
        ));
    }
    let mut used = BTreeSet::new();
    let mut receivers = Vec::with_capacity(cells.len());
    for (element, (&binding, source)) in cells.iter().zip(sources).enumerate() {
        poll(cx)?;
        let [low, high] = binding.temperature_interval_k;
        if !used.insert(binding.cell)
            || !low.is_finite()
            || !high.is_finite()
            || low <= 0.0
            || high < low
        {
            return Err(CircuitEnthalpyError::Invalid(
                "unique original-element receivers and finite positive temperature intervals required",
            ));
        }
        let network_cell =
            *network_cells
                .get(&binding.cell)
                .ok_or(CircuitEnthalpyError::Invalid(
                    "element is missing from the original heat network",
                ))?;
        if !core::ptr::eq(binding.curve, network_cell.curve) {
            return Err(CircuitEnthalpyError::Invalid(
                "element and heat network must borrow the same caloric law",
            ));
        }
        let original = thermal
            .cells()
            .get(&binding.cell)
            .ok_or(CircuitEnthalpyError::Invalid(
                "inactive harmonic heat receiver",
            ))?
            .clone();
        if !solid
            .element_reference_mass_kg(element)
            .is_some_and(|m| same_inventory(m, original.mass_kg))
            || !solid
                .element_reference_volume_m3(element)
                .is_some_and(|v| same_inventory(v, original.geometry.volume_m3))
        {
            return Err(CircuitEnthalpyError::Invalid(
                "harmonic heat needs the complete original element reference mass and volume",
            ));
        }
        let admitted = Receiver { binding, original };
        admitted.admit(thermal).map_err(receiver)?;
        admit_source(&source.material, network_cell, [low, high], cx)?;
        receivers.push(admitted);
    }
    let resolve = |state: &DepositionSimulation<'_>| {
        poll(cx)?;
        for (receiver, source) in receivers.iter().zip(sources) {
            poll(cx)?;
            receiver.admit(state).map_err(self::receiver)?;
            let binding = network_cells[&receiver.binding.cell];
            let cell = &state.cells()[&receiver.binding.cell];
            let phase = binding
                .curve
                .state_at_specific_enthalpy(cell.specific_h_j_kg)
                .map_err(|e| self::receiver(SolidThermalError::Phase(e)))?;
            let point = binding
                .point
                .at(phase)
                .map_err(|e| CircuitEnthalpyError::Source(HarmonicSolidHeatError::Material(e)))?;
            let refreshed = source
                .material
                .at_query_point(binding.point.material, &point)
                .map_err(|e| CircuitEnthalpyError::Source(HarmonicSolidHeatError::Source(e)))?;
            if refreshed.law() != source.material.law()
                || refreshed.density_kg_m3().to_bits() != source.material.density_kg_m3().to_bits()
            {
                return Err(CircuitEnthalpyError::Invalid(
                    "thermal proposal changed the prepared harmonic material",
                ));
            }
        }
        let coefficients = provider
            .resolve_heat(cx, state)
            .map_err(|e| e.map_source(HarmonicSolidHeatError::Heat))?;
        if !coefficients.contacts.is_empty() || !coefficients.bulk.resistors.is_empty() {
            return Err(CircuitEnthalpyError::Invalid(
                "harmonic material heat cannot discard electrical work",
            ));
        }
        Ok(coefficients.bulk)
    };
    resolve(thermal)?;
    let observed = observer
        .observe_window(displacements, window, cx)
        .map_err(observation)?;
    balance(
        observed.energy_balance_residual_j,
        controls.total_energy_tolerance_j,
    )?;
    let mut element_heat = vec![Dd::ZERO; cells.len()];
    for point in &observed.points {
        poll(cx)?;
        let sum = element_heat
            .get_mut(point.element)
            .ok_or(CircuitEnthalpyError::Invalid(
                "harmonic heat lost its original element",
            ))?;
        *sum =
            *sum + Dd::from_f64(point.volume_m3) * Dd::from_f64(point.material.material_heat_j_m3);
    }
    let mut total = Dd::ZERO;
    let mut powers = BTreeMap::new();
    let mut element_heat_j = Vec::with_capacity(cells.len());
    for (binding, energy) in cells.iter().zip(element_heat) {
        poll(cx)?;
        let heat = finite(energy.to_f64())?;
        let power = finite(heat / dt)?;
        if heat < 0.0 || heat > 0.0 && power == 0.0 {
            return Err(CircuitEnthalpyError::Invalid(
                "unrepresentable physical harmonic material heat",
            ));
        }
        total = total + energy;
        powers.insert(binding.cell, power);
        element_heat_j.push(heat);
    }
    balance(
        (total - Dd::from_f64(observed.material_heat_j)).to_f64(),
        controls.total_energy_tolerance_j,
    )?;
    let (next, (), report) = cell_driver_candidate(
        thermal,
        target_s,
        controls,
        initial_energy_increments_j,
        |state, _remaining| {
            let coefficients = resolve(state)?;
            Ok(DrivenCoefficients {
                circuit: ConstitutiveTransfer(&observed),
                sources_w: coefficients.sources_w,
                conductances_w_k: coefficients.conductances_w_k,
                source: coefficients.source,
                circuit_trials: 1,
            })
        },
        |_: &ConstitutiveTransfer<'_>, _, checkpoint| {
            checkpoint()?;
            Ok(powers.clone())
        },
        || cx.checkpoint().is_err(),
    )?;
    let sample = report.source.as_ref();
    closed_balance(
        dt,
        report.total_energy_residual_j,
        sample,
        controls.total_energy_tolerance_j,
    )
    .map_err(heat)?;
    let radiation_net = finite(dt * sample.radiation_net_power_w())?;
    let other_external_heat_j = finite(report.non_electrical_heat_j - radiation_net)?;
    let combined = finite(report.total_energy_residual_j + radiation_net)?;
    let handoff = finite(
        (Dd::from_f64(report.deposited_joule_j) - Dd::from_f64(observed.material_heat_j)).to_f64(),
    )?;
    let caloric = finite(
        (Dd::from_f64(report.thermal.stored_energy_change_j)
            - Dd::from_f64(report.deposited_joule_j)
            - Dd::from_f64(other_external_heat_j))
        .to_f64(),
    )?;
    balance(handoff, controls.total_energy_tolerance_j)?;
    balance(caloric, controls.total_energy_tolerance_j)?;
    poll(cx)?;
    Ok(HarmonicSolidHeatCandidate {
        thermal: next,
        report: HarmonicSolidHeatReport {
            material: observed,
            element_heat_j,
            thermal: report.thermal,
            source: report.source,
            nonlinear: report.nonlinear,
            maximum_cell_residual_j: report.maximum_cell_residual_j,
            thermal_trials: report.circuit_trials,
            deposited_heat_j: report.deposited_joule_j,
            heat_handoff_residual_j: handoff,
            other_external_heat_j,
            caloric_energy_residual_j: caloric,
            constitutive_thermal_energy_residual_j: combined,
        },
    })
}
