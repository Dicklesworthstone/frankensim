//! Held coherent Maxwell work on contiguous, well-mixed atomic gas films.
//!
//! Each original finite film supplies its signed electron/heavy friction work
//! exactly once. One native multicomponent column then evolves all charge
//! stages, thermal exchange and spectral recombination recycling together.
//! Original momentum and bound-charge sources refresh only after every layer
//! succeeds. The returned endpoint seeds the next Maxwell solve.
//!
//! The coherent drive must lie below every supplied photoionization threshold;
//! its dielectric law does not include direct bound-free absorption. It supplies
//! collisional work, not an additional ionizing beam. Diffuse photons use the
//! atomic owner's quasi-static escaping
//! boundaries, not the drive's dielectric reflections. Each film is a uniform
//! fixed-volume parcel. Spatial mixing, conduction, flow, continuously updated
//! Maxwell fields and changes of coherent field storage are outside this held
//! interval model. Original background/ion-drag cards must remain physically
//! applicable to the changing mixture; a T/p refresh alone does not prove it.

use super::IonizedGasReactingHeatControls;
use fs_em::plane_wave::layers::{MAX_LAYERS, PreparedIonizedGasLayerIllumination};
use fs_gas::chemistry::collisional_radiative::CollisionalRadiativeThermal;
use fs_gas::chemistry::diffuse_recombination::DiffuseRecombinationSolve;
use fs_gas::chemistry::ionization_network::transport::{
    NetworkRadiationCell, NetworkRadiationError, NetworkRadiationResult, NetworkRadiationSolve,
    diffuse::{NetworkDiffuseSolve, advance_network_diffuse_recombination},
};
use fs_matdb::MaterialCard;
use fs_material::ionized_gas::{
    IonizedGasConductivity, IonizedGasRefreshError, IonizedGasTransportBank,
};
use fs_material::state_point::{MaterialStatePointError, ResolvedMaterialStatePoint};
use fs_math::dd::Dd;
use fs_qty::{Area, Power, Time};
use fs_thermochem::SpeciesId;
use fs_thermochem::equilibrium::composition::FrozenIonizedGasStateV1;
use fs_thermochem::kinetics::atomic_ionization::radiative_recombination::AtomicRadiativeRecombinationV1;
use fs_thermochem::plasma_radiation::ElectronIonBremsstrahlungV1;

mod source;

/// Complete original gas, reaction, momentum and bound-charge support per film.
#[derive(Clone, Copy, Debug)]
pub struct IonizedGasLayerInput<'a> {
    pub initial: &'a FrozenIonizedGasStateV1,
    pub sources: &'a [AtomicRadiativeRecombinationV1],
    pub transport_bank: &'a IonizedGasTransportBank,
    pub original_background_card: &'a MaterialCard,
    /// Optional original frequency-integrated free-free law. Selecting it
    /// explicitly assumes all its emission escapes; no photon count is inferred.
    pub optically_thin_bremsstrahlung: Option<&'a ElectronIonBremsstrahlungV1>,
}
impl<'a> IonizedGasLayerInput<'a> {
    pub const fn new(
        initial: &'a FrozenIonizedGasStateV1,
        sources: &'a [AtomicRadiativeRecombinationV1],
        transport_bank: &'a IonizedGasTransportBank,
        original_background_card: &'a MaterialCard,
    ) -> Self {
        Self {
            initial,
            sources,
            transport_bank,
            original_background_card,
            optically_thin_bremsstrahlung: None,
        }
    }
    pub const fn with_optically_thin_bremsstrahlung(
        mut self,
        source: &'a ElectronIonBremsstrahlungV1,
    ) -> Self {
        self.optically_thin_bremsstrahlung = Some(source);
        self
    }
}

/// Native time/source controls, independent radiation work and endpoint limits.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IonizedGasLayerHeatControls {
    pub reacting: IonizedGasReactingHeatControls,
    pub radiation: DiffuseRecombinationSolve,
    /// Counts accepted and rejected field evaluations across the whole column.
    pub maximum_radiation_evaluations: usize,
    pub maximum_radiation_work: usize,
}
impl IonizedGasLayerHeatControls {
    pub fn new(
        reacting: IonizedGasReactingHeatControls,
        radiation: DiffuseRecombinationSolve,
    ) -> Self {
        let native = NetworkRadiationSolve::new(reacting.network);
        Self {
            reacting,
            radiation,
            maximum_radiation_evaluations: native.maximum_radiation_evaluations,
            maximum_radiation_work: native.maximum_radiation_work,
        }
    }
}

/// Refusal publishes no partially evolved gas, transport or background.
#[derive(Clone, Debug, PartialEq)]
pub enum LayeredIonizedGasHeatingError {
    Invalid {
        layer: Option<usize>,
        reason: &'static str,
    },
    Radiation(NetworkRadiationError),
    Transport {
        layer: usize,
        error: IonizedGasRefreshError,
    },
    Background {
        layer: usize,
        error: MaterialStatePointError,
    },
    TemperatureIncrement {
        layer: usize,
        heavy_k: f64,
        electron_k: f64,
        maximum_k: f64,
    },
    PopulationIncrement {
        layer: usize,
        species: SpeciesId,
        change_mol: f64,
        maximum_mol: f64,
    },
    EnergyBalance {
        residual_j: f64,
        tolerance_j: f64,
    },
    Cancelled,
}
impl core::fmt::Display for LayeredIonizedGasHeatingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "layered ionized gas heating: {self:?}")
    }
}
impl std::error::Error for LayeredIonizedGasHeatingError {}
type Error = LayeredIonizedGasHeatingError;

/// A single prepared coherent wave bound to every original finite gas film.
pub struct LayeredIonizedGasHeating<'w, 'm, 's> {
    wave: &'w PreparedIonizedGasLayerIllumination<'m>,
    inputs: Vec<IonizedGasLayerInput<'s>>,
    layers: Vec<source::BoundLayer<'m>>,
    area: Area,
}
impl<'w, 'm, 's> LayeredIonizedGasHeating<'w, 'm, 's> {
    /// All finite films must be contiguous gas parcels, in Maxwell layer order,
    /// with actual volume = area * thickness. Ordinary absorptive materials and
    /// zero-thickness films are refused. A lossless film is admitted only for an
    /// exactly neutral gas; its complete bank may create carriers subsequently.
    /// Every populated wave carrier must retain the identical original state,
    /// momentum law, Coulomb prescription and transport limits.
    /// The coherent photon energy must be strictly below every supplied atomic
    /// threshold, including initially absent stages. Above-threshold drive
    /// absorption requires a different coupled optical closure.
    pub fn try_new(
        wave: &'w PreparedIonizedGasLayerIllumination<'m>,
        inputs: &[IonizedGasLayerInput<'s>],
        area: Area,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        poll(&mut cancelled)?;
        if !area.value().is_finite()
            || area.value() <= 0.0
            || inputs.is_empty()
            || inputs.len() > MAX_LAYERS
            || inputs.len() != wave.layer_count()
        {
            return Err(invalid(
                None,
                "one original gas per finite film and positive physical area required",
            ));
        }
        let mut layers = Vec::with_capacity(inputs.len());
        for (i, &input) in inputs.iter().enumerate() {
            poll(&mut cancelled)?;
            layers.push(source::prepare(wave, input, area, i, &mut cancelled)?);
        }
        poll(&mut cancelled)?;
        Ok(Self {
            wave,
            inputs: inputs.to_vec(),
            layers,
            area,
        })
    }
    pub const fn wave(&self) -> &'w PreparedIonizedGasLayerIllumination<'m> {
        self.wave
    }
    pub fn inputs(&self) -> &[IonizedGasLayerInput<'s>] {
        &self.inputs
    }
    pub const fn area(&self) -> Area {
        self.area
    }

    /// Deposit each film's two signed Maxwell powers once, with native baseline
    /// electron/heavy exchange, impact ionization, capture and nonlocal spectral
    /// photon recycling. Optional free-free energy escapes independently.
    ///
    /// The source-bound network checks every attempted temperature/population
    /// state. Held-field limits compare accepted endpoints only and are not an
    /// error estimate for unchanged optical coefficients. For an initially empty
    /// electron pool, the source-supported newborn temperature is used; its old
    /// stored Te is not a physical temperature increment. The absolute species
    /// allowance still bounds all carrier birth. A new Maxwell solve must admit
    /// the refreshed transport and full coherent amplitude before another step.
    pub fn advance<'h>(
        &'h self,
        duration: Time,
        controls: IonizedGasLayerHeatControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<LayeredIonizedGasHeatStep<'h, 'w, 'm, 's>, Error> {
        poll(&mut cancelled)?;
        let c = controls.reacting;
        if !duration.value().is_finite()
            || duration.value() <= 0.0
            || !c.maximum_temperature_increment_k.is_finite()
            || c.maximum_temperature_increment_k <= 0.0
            || !c.maximum_absolute_population_increment_mol.is_finite()
            || c.maximum_absolute_population_increment_mol < 0.0
            || !c.maximum_relative_population_increment.is_finite()
            || c.maximum_relative_population_increment <= 0.0
            || c.maximum_relative_population_increment > 1.0
        {
            return Err(invalid(
                None,
                "positive duration and explicit held-field increment limits required",
            ));
        }
        let mut total = Dd::ZERO;
        let mut pool_work_magnitude = Dd::ZERO;
        let mut thermal = Vec::with_capacity(self.layers.len());
        for (i, layer) in self.layers.iter().enumerate() {
            poll(&mut cancelled)?;
            for power in [layer.electron_w, layer.heavy_w, layer.total_w] {
                let energy = power * duration.value();
                if !energy.is_finite() || (power != 0.0 && energy == 0.0) {
                    return Err(invalid(
                        Some(i),
                        "physical Maxwell heating energy is unrepresentable",
                    ));
                }
            }
            total = total + Dd::from_f64(layer.total_w) * Dd::from_f64(duration.value());
            pool_work_magnitude = pool_work_magnitude
                + (Dd::from_f64(layer.electron_w.abs()) + Dd::from_f64(layer.heavy_w.abs()))
                    * Dd::from_f64(duration.value());
            thermal.push(CollisionalRadiativeThermal {
                collisions: Some(layer.collisions.clone()),
                electron_power: Power::new(layer.electron_w),
                heavy_power: Power::new(layer.heavy_w),
                minimum_heavy_temperature_k: c.minimum_heavy_temperature_k,
                maximum_heavy_temperature_k: c.maximum_heavy_temperature_k,
                maximum_iterations: c.maximum_heavy_caloric_iterations,
            });
        }
        let cells: Vec<_> = self
            .inputs
            .iter()
            .zip(&thermal)
            .map(|(input, thermal)| {
                let cell =
                    NetworkRadiationCell::new(input.initial, input.sources).with_thermal(thermal);
                match input.optically_thin_bremsstrahlung {
                    Some(source) => cell.with_optically_thin_bremsstrahlung(source),
                    None => cell,
                }
            })
            .collect();
        let mut network = NetworkRadiationSolve::new(c.network);
        network.maximum_cells = MAX_LAYERS;
        network.maximum_radiation_evaluations = controls.maximum_radiation_evaluations;
        network.maximum_radiation_work = controls.maximum_radiation_work;
        let chemistry = advance_network_diffuse_recombination(
            &cells,
            self.area,
            &[],
            duration,
            NetworkDiffuseSolve::new(network, controls.radiation),
            &mut cancelled,
        )
        .map_err(|error| match error {
            NetworkRadiationError::Cancelled => Error::Cancelled,
            error => Error::Radiation(error),
        })?;
        let radiation = chemistry
            .diffuse_radiation()
            .ok_or_else(|| invalid(None, "native diffuse energy account missing"))?
            .balance();
        let gas = chemistry.balance();
        let work_residual = (Dd::from_f64(gas.external_work_j) - total).to_f64();
        // Independent boundary account: internal recycling is not an external
        // heat source. Free-free escape is separate from bound-free photons.
        let energy_residual_j = (Dd::from_f64(gas.gas_energy_change_j)
            + Dd::from_f64(radiation.escaped_energy_j)
            + Dd::from_f64(gas.free_free_escaped_energy_j)
            - total)
            .to_f64();
        // Oppositely signed Coulomb partitions can be much larger than their
        // total. Include only their floating-point summation allowance here.
        let roundoff = 128.0
            * f64::EPSILON
            * (gas.external_work_j.abs() + total.to_f64().abs() + pool_work_magnitude.to_f64());
        let tolerance = c.network.local.energy_tolerance_j + roundoff;
        for (residual_j, tolerance_j) in [
            (work_residual, tolerance),
            (
                energy_residual_j,
                tolerance + radiation.boundary_energy_roundoff_j,
            ),
        ] {
            if !residual_j.is_finite() || !tolerance_j.is_finite() || residual_j.abs() > tolerance_j
            {
                return Err(Error::EnergyBalance {
                    residual_j,
                    tolerance_j,
                });
            }
        }
        let mut refreshed = Vec::with_capacity(self.inputs.len());
        for (i, cell) in chemistry.cells().iter().enumerate() {
            poll(&mut cancelled)?;
            refreshed.push(source::refresh(
                self.inputs[i],
                &self.layers[i],
                cell.state(),
                c,
                i,
                &mut cancelled,
            )?);
        }
        poll(&mut cancelled)?;
        Ok(LayeredIonizedGasHeatStep {
            heating: self,
            chemistry,
            refreshed,
            absorbed_energy_j: total.to_f64(),
            energy_residual_j,
            controls,
        })
    }
}

/// Original laws evaluated at one accepted layer endpoint.
pub struct RefreshedIonizedGasLayer {
    transport: Option<IonizedGasConductivity>,
    background: ResolvedMaterialStatePoint,
}
impl RefreshedIonizedGasLayer {
    /// None means exact free-carrier extinction; no conductivity floor is used.
    pub const fn transport(&self) -> Option<&IonizedGasConductivity> {
        self.transport.as_ref()
    }
    pub const fn background(&self) -> &ResolvedMaterialStatePoint {
        &self.background
    }
}

/// Atomic publication of the complete column and next Maxwell material inputs.
pub struct LayeredIonizedGasHeatStep<'h, 'w, 'm, 's> {
    heating: &'h LayeredIonizedGasHeating<'w, 'm, 's>,
    chemistry: NetworkRadiationResult,
    refreshed: Vec<RefreshedIonizedGasLayer>,
    absorbed_energy_j: f64,
    energy_residual_j: f64,
    controls: IonizedGasLayerHeatControls,
}
impl<'h, 'w, 'm, 's> LayeredIonizedGasHeatStep<'h, 'w, 'm, 's> {
    pub const fn heating(&self) -> &'h LayeredIonizedGasHeating<'w, 'm, 's> {
        self.heating
    }
    pub const fn chemistry(&self) -> &NetworkRadiationResult {
        &self.chemistry
    }
    pub fn refreshed_layers(&self) -> &[RefreshedIonizedGasLayer] {
        &self.refreshed
    }
    pub const fn absorbed_energy_j(&self) -> f64 {
        self.absorbed_energy_j
    }
    /// Gas storage plus both escaped radiation energies minus exact wave work.
    pub const fn energy_residual_j(&self) -> f64 {
        self.energy_residual_j
    }
    pub const fn controls(&self) -> IonizedGasLayerHeatControls {
        self.controls
    }
}

fn invalid(layer: Option<usize>, reason: &'static str) -> Error {
    Error::Invalid { layer, reason }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
fn transport_error(layer: usize, error: IonizedGasRefreshError) -> Error {
    match error {
        IonizedGasRefreshError::Cancelled => Error::Cancelled,
        error => Error::Transport { layer, error },
    }
}
