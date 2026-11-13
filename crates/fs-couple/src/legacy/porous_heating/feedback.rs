//! Quasisteady harmonic absorption refreshed by native cell temperatures.
//!
//! Each thermally lumped specimen cell is one homogeneous acoustic sublayer.
//! Its actual caloric-law temperature reevaluates the ORIGINAL fluid source
//! and PorousFluidLaw. Density, sound speed, viscosity, conductivity and pore
//! resistivity change together; no independent temperature multiplier is fitted.
//! Changed layers then participate in ONE coherent two-sided stack solve.
//!
//! Pressure and composition are prescribed isobaric background conditions, not
//! a finite fluid inventory. Thermal mass is never replaced by rho(T)*pore volume.
//! The original thermal law must represent the specimen's storage under this
//! background, with a justified common fluid/skeleton temperature. Reservoir mass
//! transport, mean pressure evolution and constituent heat fractions are absent.
//!
//! Every heat rate uses the accepted START temperature of its native explicit
//! step. A changed endpoint prepares the next start and admits the final state;
//! it never retroactively changes deposited heat. The coherent drive and real
//! exterior wave reference stay fixed. Acoustic settling must be fast compared
//! with thermal evolution: startup/storage transients and parametric acoustic
//! work are outside this cycle-averaged quasisteady approximation. Native
//! conduction CFL is unchanged; independently refine max_dt for thermal feedback.

use super::*;
use fs_duct::compose::ComposeOptions;
use fs_duct::porous::{
    PassivePorousTermination, PorousLayer, PorousLayerPower, PorousLayerStack,
    PorousStackExcitation,
};
use fs_material::gas::{GasSpec, GasState};
use fs_material::porous::PorousFluidLaw;
use fs_math::c64::C64;

/// Original thermoviscous material laws for accepted-temperature feedback.
pub mod fluid;
pub use fluid::{IsobaricPorousLiquid, PorousLiquidState, PorousTemperatureFluid};

/// Original molecular calorics refreshed at accepted pore-gas temperatures.
pub mod molecular;
pub use molecular::{IsobaricRelaxingPorousGas, PorousRelaxingState};

/// Partially coherent two-sided stationary pressure spectra.
pub mod covariance;
/// Statistically stationary independent frequencies on one thermal specimen.
pub mod spectrum;

/// One original calorically perfect/Sutherland gas law at prescribed pressure.
/// The explicit temperature support is additionally admitted by the original
/// gas constructor at both ends; no source domain is extended or defaulted.
#[derive(Clone, Copy, Debug)]
pub struct IsobaricPorousGas<'a> {
    spec: &'a GasSpec,
    pressure_pa: f64,
    temperature_support_k: [f64; 2],
}
impl<'a> IsobaricPorousGas<'a> {
    /// Admit a caller-declared pressure and thermal support through the native
    /// gas EOS, caloric and transport checks at both support endpoints.
    pub fn try_new(
        spec: &'a GasSpec,
        pressure_pa: f64,
        temperature_support_k: [f64; 2],
    ) -> Result<Self, PorousHeatingError> {
        if !temperature_support_k
            .into_iter()
            .all(|t| t.is_finite() && t > 0.0)
            || temperature_support_k[0] > temperature_support_k[1]
        {
            return Err(PorousHeatingError::Invalid(
                "isobaric pore gas needs explicit ordered temperature support",
            ));
        }
        for temperature in temperature_support_k {
            GasState::try_new(spec, temperature, pressure_pa)?.validate_coherence()?;
        }
        Ok(Self {
            spec,
            pressure_pa,
            temperature_support_k,
        })
    }
    /// Original fixed composition, caloric and transport specification.
    pub const fn spec(self) -> &'a GasSpec {
        self.spec
    }
    /// Prescribed background pressure [Pa], not an evolving acoustic pressure.
    pub const fn pressure_pa(self) -> f64 {
        self.pressure_pa
    }
    /// Explicit admitted inclusive temperature support [K].
    pub const fn temperature_support_k(self) -> [f64; 2] {
        self.temperature_support_k
    }
    fn evaluate(self, temperature_k: f64) -> Result<GasState, PorousHeatingError> {
        if !temperature_k.is_finite()
            || !(self.temperature_support_k[0]..=self.temperature_support_k[1])
                .contains(&temperature_k)
        {
            return Err(PorousHeatingError::Invalid(
                "thermal pore-gas temperature left its declared support",
            ));
        }
        let state = GasState::try_new(self.spec, temperature_k, self.pressure_pa)?;
        state.validate_coherence()?;
        Ok(state)
    }
}

/// One original thermal cell and its unchanged porous skeleton material.
/// Entries run from the left incident face to the right face. Physical lengths,
/// volumes, masses and caloric laws come only from the original thermal cells.
#[derive(Clone, Copy, Debug)]
pub struct PorousTemperatureCell<'a> {
    /// Stable existing thermal control volume, in incident-face-first order.
    pub cell: CellId,
    /// Original fixed skeleton, pore geometry and characterized drag law.
    pub material: &'a PorousFluidLaw,
}

/// Fixed coherent exterior drive; heating does not silently change its power
/// normalization when the INTERNAL pore-fluid temperature and impedance change.
#[derive(Clone, Copy, Debug)]
pub struct PorousHarmonicDrive {
    /// Shared angular frequency [rad/s], in each skeleton's original support.
    pub angular_frequency_rad_s: f64,
    /// Common full specimen face area [m²].
    pub area_m2: f64,
    /// Fixed positive real exterior wave reference [Pa s/m³].
    pub reference_impedance_pa_s_m3: f64,
    /// Coherent peak pressure phasors at the two EXTERIOR faces [Pa].
    pub incident_pressure_pa: [C64; 2],
}

/// Reusable accepted-state cache plus the original thermal/material bindings.
/// Each changed temperature reevaluates only its own fluid and porous sublayer.
/// A shared stack solve runs once if any cell changed; unchanged/latent plateau
/// temperatures reuse all coefficients and analytic sources. Geometry and
/// material objects remain immutable, including during candidate cloning.
#[derive(Clone)]
pub struct TemperaturePorousHeating<'a, F: PorousTemperatureFluid = IsobaricPorousGas<'a>> {
    cells: Vec<PorousTemperatureCell<'a>>,
    gas: F,
    drive: PorousHarmonicDrive,
    tangential_wavenumber_rad_m: f64,
    termination: Option<PassivePorousTermination>,
    options: ComposeOptions,
    states: Vec<F::State>,
    layers: Vec<PorousLayer<'a>>,
    heat: PreparedHeat,
    boundary: PorousLayerPower,
}

/// Additional evaluations in one advance, excluding initial construction.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct PorousFeedbackWork {
    /// Changed physical cell states evaluated through any original fluid law.
    pub fluid_evaluations: usize,
    /// Gas evaluations, retained separately for existing gas-only consumers.
    pub gas_evaluations: usize,
    /// Shared coherent stack solves, including admitted changed endpoints.
    pub stack_solves: usize,
}
impl PorousFeedbackWork {
    fn record_fluid<F: PorousTemperatureFluid>(&mut self) -> Result<(), PorousHeatingError> {
        self.fluid_evaluations = self
            .fluid_evaluations
            .checked_add(1)
            .ok_or(DepositionError::Budget("porous fluid evaluations"))?;
        if F::IS_GAS {
            self.gas_evaluations = self
                .gas_evaluations
                .checked_add(1)
                .ok_or(DepositionError::Budget("porous gas evaluations"))?;
        }
        Ok(())
    }
}

impl<'a> TemperaturePorousHeating<'a> {
    /// Original gas states, preserving the existing gas-only view.
    pub fn gas_states(&self) -> &[GasState] {
        self.fluid_states()
    }
    /// Prescribed isobaric gas law and its explicit admitted temperature range.
    pub const fn gas(&self) -> IsobaricPorousGas<'a> {
        self.gas
    }
}

/// One jointly accepted cache and native thermal candidate.
pub struct PorousTemperatureStep<'t, 'a, F: PorousTemperatureFluid = IsobaricPorousGas<'a>> {
    /// Reusable cache evaluated at the returned thermal state's temperatures.
    pub heating: TemperaturePorousHeating<'a, F>,
    /// Existing native thermal and acoustic energy accounts, without duplication.
    pub step: PorousHeatingStep<'t>,
    /// Actual changed-state/source work in this advance.
    pub work: PorousFeedbackWork,
}

// Carry the caller's exact thermal cancellation/refusal through a native
// acoustic bool cancellation gate without replacing it by a generic failure.
fn acoustic<T>(
    checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    solve: impl FnOnce(&mut dyn FnMut() -> bool) -> Result<T, DuctError>,
) -> Result<T, PorousHeatingError> {
    let mut refusal = None;
    let result = solve(&mut || match checkpoint() {
        Ok(()) => false,
        Err(error) => {
            refusal = Some(error);
            true
        }
    });
    if let Some(error) = refusal {
        Err(error.into())
    } else {
        result.map_err(Into::into)
    }
}

fn admit_tangential_wavenumber(value: f64) -> Result<(), PorousHeatingError> {
    if value.is_finite() && value >= 0.0 {
        Ok(())
    } else {
        Err(PorousHeatingError::Invalid(
            "porous tangential wavenumber must be finite and nonnegative",
        ))
    }
}

fn prepare_stack_drive<'s, 'l, 'm>(
    stack: &'s PorousLayerStack<'l, 'm>,
    incident: [C64; 2],
    termination: Option<PassivePorousTermination>,
    cancelled: &mut dyn FnMut() -> bool,
) -> Result<PorousStackExcitation<'s, 'l, 'm>, DuctError> {
    match termination {
        Some(load) => stack.prepare_terminated_excitation(incident[0], load, cancelled),
        None => stack.prepare_excitation(incident, cancelled),
    }
}

impl<'a, F: PorousTemperatureFluid> TemperaturePorousHeating<'a, F> {
    /// One-sided normal-incidence illumination with a passive rear load.
    /// The authored right pressure must be zero. The fixed load's counter-wave
    /// is recomputed with every changed-temperature stack; its dissipated work
    /// is exported separately from the original porous thermal cells.
    pub fn try_new_terminated(
        thermal: &DepositionSimulation<'_>,
        cells: &[PorousTemperatureCell<'a>],
        gas: F,
        drive: PorousHarmonicDrive,
        termination: PassivePorousTermination,
        options: ComposeOptions,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, PorousHeatingError> {
        Self::try_new_oblique_terminated(
            thermal,
            cells,
            gas,
            drive,
            0.0,
            termination,
            options,
            checkpoint,
        )
    }

    /// Conserved oblique mode with a passive rear pressure/volume-flow load.
    /// The load impedance must describe this normal mode at the rear face;
    /// neither incidence geometry nor an air gap is inferred from it.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new_oblique_terminated(
        thermal: &DepositionSimulation<'_>,
        cells: &[PorousTemperatureCell<'a>],
        gas: F,
        drive: PorousHarmonicDrive,
        tangential_wavenumber_rad_m: f64,
        termination: PassivePorousTermination,
        options: ComposeOptions,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, PorousHeatingError> {
        Self::prepare(
            thermal,
            cells,
            gas,
            drive,
            tangential_wavenumber_rad_m,
            Some(termination),
            options,
            None,
            &mut checkpoint,
        )
    }

    /// Bind the complete ordered specimen and evaluate its actual initial
    /// temperatures. Every cell becomes one original porous acoustic sublayer.
    /// The native layer/stack heat adapters admit full-area prismatic geometry,
    /// common body/normal, unique cells, original storage identities and power.
    /// The original composition budget and 4096-cell heat ceiling apply before
    /// allocation. No source spectrum, caloric law or conductance is replaced.
    pub fn try_new(
        thermal: &DepositionSimulation<'_>,
        cells: &[PorousTemperatureCell<'a>],
        gas: F,
        drive: PorousHarmonicDrive,
        options: ComposeOptions,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, PorousHeatingError> {
        Self::try_new_oblique(thermal, cells, gas, drive, 0.0, options, checkpoint)
    }

    /// Bind an oblique plane-wave drive to the same native thermal specimen.
    /// `k_t` [rad/m] and the drive's NORMAL exterior wave reference stay fixed
    /// while internal gas temperature changes. The caller supplies a physical
    /// propagating exterior mode; original normal thickness/full area apply.
    /// Each layer retains its normal and transverse viscous absorption.
    pub fn try_new_oblique(
        thermal: &DepositionSimulation<'_>,
        cells: &[PorousTemperatureCell<'a>],
        gas: F,
        drive: PorousHarmonicDrive,
        tangential_wavenumber_rad_m: f64,
        options: ComposeOptions,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, PorousHeatingError> {
        Self::prepare(
            thermal,
            cells,
            gas,
            drive,
            tangential_wavenumber_rad_m,
            None,
            options,
            None,
            &mut checkpoint,
        )
    }

    // Spectral bins share fluid evaluations while retaining the same original
    // single-frequency preparation and complete physical thermal admission.
    #[allow(clippy::too_many_arguments)]
    fn prepare(
        thermal: &DepositionSimulation<'_>,
        cells: &[PorousTemperatureCell<'a>],
        gas: F,
        drive: PorousHarmonicDrive,
        tangential_wavenumber_rad_m: f64,
        termination: Option<PassivePorousTermination>,
        options: ComposeOptions,
        shared_states: Option<&[F::State]>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, PorousHeatingError> {
        checkpoint()?;
        admit_tangential_wavenumber(tangential_wavenumber_rad_m)?;
        if termination.is_some() && drive.incident_pressure_pa[1] != C64::ZERO {
            return Err(PorousHeatingError::Invalid(
                "terminated porous drive requires one-sided left illumination",
            ));
        }
        if cells.is_empty()
            || cells.len() > MAX_POROUS_HEATING_CELLS
            || cells.len() > options.maximum_blocks
        {
            return Err(DepositionError::Budget("temperature-dependent porous cells").into());
        }
        if shared_states.is_some_and(|states| states.len() != cells.len()) {
            return Err(PorousHeatingError::Invalid("shared porous gas state count"));
        }
        let mut states = Vec::with_capacity(cells.len());
        let mut layers = Vec::with_capacity(cells.len());
        let mut partitions = Vec::with_capacity(cells.len());
        for (index, cell) in cells.iter().enumerate() {
            checkpoint()?;
            let original = thermal
                .cells()
                .get(&cell.cell)
                .ok_or(PorousHeatingError::Invalid(
                    "temperature-dependent porous cell is inactive",
                ))?;
            let length = finite(norm(std::array::from_fn(|i| {
                original.geometry.end_m[i] - original.geometry.start_m[i]
            })))?;
            let state = if let Some(states) = shared_states {
                states[index].clone()
            } else {
                fluid::evaluate(&gas, thermal.temperature_k(cell.cell)?)?
            };
            layers.push(gas.prepare_layer(
                cell.material,
                &state,
                drive.angular_frequency_rad_s,
                length,
                drive.area_m2,
                drive.reference_impedance_pa_s_m3,
                tangential_wavenumber_rad_m,
            )?);
            states.push(state);
            partitions.push([PorousHeatCell {
                cell: cell.cell,
                end_depth_m: length,
            }]);
        }
        let (heat, boundary) = {
            let stack = acoustic(checkpoint, |cancelled| {
                PorousLayerStack::try_new(&layers, options, cancelled)
            })?;
            let excitation = acoustic(checkpoint, |cancelled| {
                prepare_stack_drive(&stack, drive.incident_pressure_pa, termination, cancelled)
            })?;
            let partitions: Vec<&[PorousHeatCell]> =
                partitions.iter().map(|p| p.as_slice()).collect();
            let bound = PorousStackHeating::try_new(
                &excitation,
                thermal,
                &partitions,
                cells.len(),
                &mut *checkpoint,
            )?;
            let boundary = *excitation.boundary_power();
            let source = excitation.source_power();
            (
                PreparedHeat {
                    cells: bound.cells().to_vec(),
                    sources_w: bound
                        .cells()
                        .iter()
                        .map(|cell| (cell.cell, cell.dissipation.absorbed_power_w))
                        .collect(),
                    dissipation: bound.dissipation(),
                    power_residual_w: source.residual_w,
                    boundary_power_w: [source.incident_power_w, source.outgoing_power_w],
                    load_absorbed_power_w: source.load_absorbed_power_w,
                    load_receiver: None,
                },
                boundary,
            )
        };
        checkpoint()?;
        Ok(Self {
            cells: cells.to_vec(),
            gas,
            drive,
            termination,
            tangential_wavenumber_rad_m: if tangential_wavenumber_rad_m == 0.0 {
                0.0
            } else {
                tangential_wavenumber_rad_m
            },
            options,
            states,
            layers,
            heat,
            boundary,
        })
    }

    /// Deposit all reported passive-load work in one distinct native backing
    /// cell, in the same accepted steps as pore heat. Its native conduction may
    /// change pore temperatures and hence acoustics. Its own temperature does
    /// not alter the authored fixed load impedance or trigger an acoustic solve.
    /// The backing retains its original mass, caloric law, sources and events.
    /// Bind only work that thermalizes in this body; escaping radiation stays
    /// exported. No heat fraction or dissipation mechanism is inferred.
    pub fn with_load_receiver(
        mut self,
        thermal: &DepositionSimulation<'_>,
        cell: CellId,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, PorousHeatingError> {
        self.heat
            .bind_load_receiver(thermal, cell, &mut checkpoint)?;
        Ok(self)
    }

    /// Explicit backing destination, separate from the porous cell partition.
    pub fn load_receiver(&self) -> Option<CellId> {
        self.heat.load_receiver.map(|receiver| receiver.cell)
    }

    /// Physical thermal destinations and their current analytic heat rates.
    pub fn cells(&self) -> &[PorousCellHeating] {
        &self.heat.cells
    }
    /// Original fluid states evaluated at the current cached temperatures.
    pub fn fluid_states(&self) -> &[F::State] {
        &self.states
    }
    /// Original per-cell porous responses, retaining each skeleton's identity.
    pub fn acoustic_layers(&self) -> &[PorousLayer<'a>] {
        &self.layers
    }
    /// Original fluid source and its prescribed background conditions.
    pub const fn fluid(&self) -> &F {
        &self.gas
    }
    /// Original coherent pressure drive and FIXED exterior normalization.
    pub const fn drive(&self) -> PorousHarmonicDrive {
        self.drive
    }
    /// Fixed conserved tangential wavenumber of this physical frequency [rad/m].
    pub const fn tangential_wavenumber_rad_m(&self) -> f64 {
        self.tangential_wavenumber_rad_m
    }
    /// Current physical work, after coherent illumination from both exteriors.
    pub const fn dissipation(&self) -> PorousLayerDissipation {
        self.heat.dissipation
    }
    /// Raw pressure/flow power at both stack faces. With a rear load this
    /// includes its returning wave; use `power` for prescribed-source energy.
    pub const fn boundary_power(&self) -> &PorousLayerPower {
        &self.boundary
    }
    /// Original fixed rear load, or `None` for two prescribed open faces.
    pub const fn termination(&self) -> Option<PassivePorousTermination> {
        self.termination
    }
    /// Prescribed-source power and separately exported terminal work.
    pub fn power(&self) -> spectrum::PorousSpectrumPower {
        spectrum::PorousSpectrumPower {
            incident_power_w: self.heat.boundary_power_w[0],
            outgoing_power_w: self.heat.boundary_power_w[1],
            load_absorbed_power_w: self.heat.load_absorbed_power_w,
            load_deposited_power_w: self.heat.load_deposited_power_w(),
            load_exported_power_w: self.heat.load_exported_power_w(),
            dissipation: self.heat.dissipation,
            residual_w: self.heat.power_residual_w,
        }
    }

    fn refresh(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        work: &mut PorousFeedbackWork,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), PorousHeatingError> {
        self.heat.validate(thermal, checkpoint)?;
        let mut changed = Vec::new();
        for (i, cell) in self.cells.iter().enumerate() {
            checkpoint()?;
            let temperature = thermal.temperature_k(cell.cell)?;
            if temperature.to_bits() == self.gas.temperature_k(&self.states[i]).to_bits() {
                continue;
            }
            let state = fluid::evaluate(&self.gas, temperature)?;
            changed.push((i, state));
            work.record_fluid::<F>()?;
        }
        self.refresh_changed(&changed, work, checkpoint)
    }

    // Only internally evaluated native states enter here. The single-frequency
    // caller or common spectral binding already admitted the thermal specimen.
    fn refresh_changed(
        &mut self,
        changed: &[(usize, F::State)],
        work: &mut PorousFeedbackWork,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), PorousHeatingError> {
        if changed.is_empty() {
            return Ok(());
        }
        for (i, state) in changed {
            checkpoint()?;
            let layer = self.gas.prepare_layer(
                self.cells[*i].material,
                state,
                self.drive.angular_frequency_rad_s,
                self.layers[*i].thickness_m(),
                self.drive.area_m2,
                self.drive.reference_impedance_pa_s_m3,
                self.tangential_wavenumber_rad_m,
            )?;
            self.states[*i] = state.clone();
            self.layers[*i] = layer;
        }
        let stack = acoustic(checkpoint, |cancelled| {
            PorousLayerStack::try_new(&self.layers, self.options, cancelled)
        })?;
        let excitation = acoustic(checkpoint, |cancelled| {
            prepare_stack_drive(
                &stack,
                self.drive.incident_pressure_pa,
                self.termination,
                cancelled,
            )
        })?;
        for (bound, local) in self
            .heat
            .cells
            .iter_mut()
            .zip(excitation.layer_excitations())
        {
            checkpoint()?;
            bound.dissipation = local.dissipation_between(0.0, local.layer().thickness_m())?;
            *self
                .heat
                .sources_w
                .get_mut(&bound.cell)
                .expect("original admitted cell source") = bound.dissipation.absorbed_power_w;
        }
        self.heat.dissipation = excitation.dissipation();
        self.heat.power_residual_w = excitation.source_power().residual_w;
        self.boundary = *excitation.boundary_power();
        self.heat.boundary_power_w = [
            excitation.source_power().incident_power_w,
            excitation.source_power().outgoing_power_w,
        ];
        self.heat.load_absorbed_power_w = excitation.source_power().load_absorbed_power_w;
        work.stack_solves = work
            .stack_solves
            .checked_add(1)
            .ok_or(DepositionError::Budget("porous stack evaluations"))?;
        checkpoint()?;
        Ok(())
    }

    /// Advance source feedback and thermal state as ONE immutable candidate.
    /// Rates use accepted-start temperatures, then native duration/CFL/event
    /// clipping applies. Endpoint refresh admits the returned fluid state and
    /// caches it for continuation. Conduction, authored heat, latent storage,
    /// original material/source domains and all rollback rules remain native.
    /// Fluid and acoustic refreshes add no separate heat credit or storage term.
    pub fn advance_to<'t>(
        &self,
        initial: &DepositionSimulation<'t>,
        target_s: f64,
        config: StepConfig,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<PorousTemperatureStep<'t, 'a, F>, PorousHeatingError> {
        checkpoint()?;
        let mut observer = FeedbackObserver {
            source: self.clone(),
            account: HeatAccount::new(self.cells.len()),
            work: PorousFeedbackWork::default(),
        };
        observer
            .source
            .refresh(initial, &mut observer.work, &mut checkpoint)?;
        let (thermal, observer, heat_report) =
            initial.advance_observed(&observer, target_s, config, &mut checkpoint)?;
        let step = observer.account.finish(
            &observer.source.heat,
            thermal,
            heat_report,
            &mut checkpoint,
        )?;
        checkpoint()?;
        Ok(PorousTemperatureStep {
            heating: observer.source,
            step,
            work: observer.work,
        })
    }
}

#[derive(Clone)]
struct FeedbackObserver<'a, F: PorousTemperatureFluid> {
    source: TemperaturePorousHeating<'a, F>,
    account: HeatAccount,
    work: PorousFeedbackWork,
}
impl<F: PorousTemperatureFluid> DepositionObserver for FeedbackObserver<'_, F> {
    type Error = PorousHeatingError;
    fn sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.source.refresh(thermal, &mut self.work, checkpoint)?;
        Ok(self.source.heat.sources())
    }
    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        Ok(())
    }
    fn step(
        &mut self,
        _: f64,
        dt: f64,
        _: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        // Charge the heat that actually drove this step BEFORE replacing its
        // accepted-start coefficients with the candidate endpoint response.
        self.account.record(&self.source.heat, dt, checkpoint)?;
        self.source.refresh(thermal, &mut self.work, checkpoint)
    }
}
