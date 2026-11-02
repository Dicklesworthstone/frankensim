//! Quasisteady harmonic absorption refreshed by native cell temperatures.
//!
//! Each thermally lumped specimen cell is one homogeneous acoustic sublayer.
//! Its actual caloric-law temperature reevaluates the ORIGINAL GasSpec/GasState
//! and PorousFluidLaw. Density, sound speed, viscosity, conductivity and pore
//! resistivity change together; no independent temperature multiplier is fitted.
//! Changed layers then participate in ONE coherent two-sided stack solve.
//!
//! Pressure and composition are prescribed isobaric background conditions, not
//! a finite gas inventory. Thermal mass is never replaced by rho(T)*pore volume.
//! The original thermal law must represent the specimen's storage under this
//! background, with a justified common gas/skeleton temperature. Reservoir mass
//! transport, gas pressure evolution and constituent heat fractions are absent.
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
use fs_duct::porous::{PorousLayer, PorousLayerPower, PorousLayerStack};
use fs_material::gas::{GasSpec, GasState};
use fs_material::porous::PorousFluidLaw;
use fs_math::c64::C64;

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
/// normalization when the INTERNAL pore-gas temperature and impedance change.
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
/// Each changed temperature reevaluates only its own gas and porous sublayer.
/// A shared stack solve runs once if any cell changed; unchanged/latent plateau
/// temperatures reuse all coefficients and analytic sources. Geometry and
/// material objects remain immutable, including during candidate cloning.
#[derive(Clone)]
pub struct TemperaturePorousHeating<'a> {
    cells: Vec<PorousTemperatureCell<'a>>,
    gas: IsobaricPorousGas<'a>,
    drive: PorousHarmonicDrive,
    options: ComposeOptions,
    states: Vec<GasState>,
    layers: Vec<PorousLayer<'a>>,
    heat: PreparedHeat,
    boundary: PorousLayerPower,
}

/// Additional evaluations in one advance, excluding initial construction.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct PorousFeedbackWork {
    /// Changed physical cell states evaluated through the original gas law.
    pub gas_evaluations: usize,
    /// Shared coherent stack solves, including admitted changed endpoints.
    pub stack_solves: usize,
}

/// One jointly accepted cache and native thermal candidate.
pub struct PorousTemperatureStep<'t, 'a> {
    /// Reusable cache evaluated at the returned thermal state's temperatures.
    pub heating: TemperaturePorousHeating<'a>,
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

impl<'a> TemperaturePorousHeating<'a> {
    /// Bind the complete ordered specimen and evaluate its actual initial
    /// temperatures. Every cell becomes one original porous acoustic sublayer.
    /// The native layer/stack heat adapters admit full-area prismatic geometry,
    /// common body/normal, unique cells, original storage identities and power.
    /// The original composition budget and 4096-cell heat ceiling apply before
    /// allocation. No source spectrum, caloric law or conductance is replaced.
    pub fn try_new(
        thermal: &DepositionSimulation<'_>,
        cells: &[PorousTemperatureCell<'a>],
        gas: IsobaricPorousGas<'a>,
        drive: PorousHarmonicDrive,
        options: ComposeOptions,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, PorousHeatingError> {
        Self::prepare(thermal, cells, gas, drive, options, None, &mut checkpoint)
    }

    // Spectral bins share gas evaluations while retaining the same original
    // single-frequency preparation and complete physical thermal admission.
    fn prepare(
        thermal: &DepositionSimulation<'_>,
        cells: &[PorousTemperatureCell<'a>],
        gas: IsobaricPorousGas<'a>,
        drive: PorousHarmonicDrive,
        options: ComposeOptions,
        shared_states: Option<&[GasState]>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, PorousHeatingError> {
        checkpoint()?;
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
                states[index]
            } else {
                gas.evaluate(thermal.temperature_k(cell.cell)?)?
            };
            layers.push(PorousLayer::try_new(
                cell.material,
                &state,
                drive.angular_frequency_rad_s,
                length,
                drive.area_m2,
                drive.reference_impedance_pa_s_m3,
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
                stack.prepare_excitation(drive.incident_pressure_pa, cancelled)
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
            (
                PreparedHeat {
                    cells: bound.cells().to_vec(),
                    sources_w: bound
                        .cells()
                        .iter()
                        .map(|cell| (cell.cell, cell.dissipation.absorbed_power_w))
                        .collect(),
                    dissipation: bound.dissipation(),
                    power_residual_w: bound.power_residual_w(),
                    boundary_power_w: [boundary.incident_power_w, boundary.outgoing_power_w],
                },
                boundary,
            )
        };
        checkpoint()?;
        Ok(Self {
            cells: cells.to_vec(),
            gas,
            drive,
            options,
            states,
            layers,
            heat,
            boundary,
        })
    }

    /// Physical thermal destinations and their current analytic heat rates.
    pub fn cells(&self) -> &[PorousCellHeating] {
        &self.heat.cells
    }
    /// Original gas states evaluated at the current cached temperatures.
    pub fn gas_states(&self) -> &[GasState] {
        &self.states
    }
    /// Original per-cell porous responses, retaining each skeleton's identity.
    pub fn acoustic_layers(&self) -> &[PorousLayer<'a>] {
        &self.layers
    }
    /// Prescribed isobaric gas law and its explicit admitted temperature range.
    pub const fn gas(&self) -> IsobaricPorousGas<'a> {
        self.gas
    }
    /// Original coherent pressure drive and FIXED exterior normalization.
    pub const fn drive(&self) -> PorousHarmonicDrive {
        self.drive
    }
    /// Current physical work, after coherent illumination from both exteriors.
    pub const fn dissipation(&self) -> PorousLayerDissipation {
        self.heat.dissipation
    }
    /// Exterior power only; local interface exchanges are never external sources.
    pub const fn boundary_power(&self) -> &PorousLayerPower {
        &self.boundary
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
            if temperature.to_bits() == self.states[i].temperature.to_bits() {
                continue;
            }
            let state = self.gas.evaluate(temperature)?;
            changed.push((i, state));
            work.gas_evaluations = work
                .gas_evaluations
                .checked_add(1)
                .ok_or(DepositionError::Budget("porous gas evaluations"))?;
        }
        self.refresh_changed(&changed, work, checkpoint)
    }

    // Only internally evaluated native states enter here. The single-frequency
    // caller or common spectral binding already admitted the thermal specimen.
    fn refresh_changed(
        &mut self,
        changed: &[(usize, GasState)],
        work: &mut PorousFeedbackWork,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), PorousHeatingError> {
        if changed.is_empty() {
            return Ok(());
        }
        for &(i, state) in changed {
            checkpoint()?;
            let layer = PorousLayer::try_new(
                self.cells[i].material,
                &state,
                self.drive.angular_frequency_rad_s,
                self.layers[i].thickness_m(),
                self.drive.area_m2,
                self.drive.reference_impedance_pa_s_m3,
            )?;
            self.states[i] = state;
            self.layers[i] = layer;
        }
        let stack = acoustic(checkpoint, |cancelled| {
            PorousLayerStack::try_new(&self.layers, self.options, cancelled)
        })?;
        let excitation = acoustic(checkpoint, |cancelled| {
            stack.prepare_excitation(self.drive.incident_pressure_pa, cancelled)
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
        self.heat.power_residual_w = excitation.power_residual_w();
        self.boundary = *excitation.boundary_power();
        self.heat.boundary_power_w = [
            self.boundary.incident_power_w,
            self.boundary.outgoing_power_w,
        ];
        work.stack_solves = work
            .stack_solves
            .checked_add(1)
            .ok_or(DepositionError::Budget("porous stack evaluations"))?;
        checkpoint()?;
        Ok(())
    }

    /// Advance source feedback and thermal state as ONE immutable candidate.
    /// Rates use accepted-start temperatures, then native duration/CFL/event
    /// clipping applies. Endpoint refresh admits the returned gas state and
    /// caches it for continuation. Conduction, authored heat, latent storage,
    /// original material/source domains and all rollback rules remain native.
    /// Gas and acoustic refreshes add no separate heat credit or storage term.
    pub fn advance_to<'t>(
        &self,
        initial: &DepositionSimulation<'t>,
        target_s: f64,
        config: StepConfig,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<PorousTemperatureStep<'t, 'a>, PorousHeatingError> {
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
struct FeedbackObserver<'a> {
    source: TemperaturePorousHeating<'a>,
    account: HeatAccount,
    work: PorousFeedbackWork,
}
impl DepositionObserver for FeedbackObserver<'_> {
    type Error = PorousHeatingError;
    fn sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.source.refresh(thermal, &mut self.work, checkpoint)?;
        Ok(self.source.heat.sources_w.clone())
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
