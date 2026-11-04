//! Molecular acoustic absorption on the original fixed-volume gas storage.
//!
//! A prepared BEM source bank already owns the coherent pressure field, the
//! molecular compression law and positive physical volume quadrature. This
//! adapter turns its cycle-averaged watts into joules exactly once and sends
//! them through the original parcel's nonlinear internal-energy closure.
//! Fixed volume uses `u(T)` and `Cv`; no solid law or constant-`Cp` surrogate
//! owns the absorbed gas heat.
//!
//! The harmonic field is frozen during each declared heating interval. A
//! returned next relaxation state can initialize a NEW acoustic solve; it
//! cannot update the old source bank in place. This is a uniformly mixed gas
//! storage model, not spatial conduction, mean flow or transient acoustics.
//! Quadrature convergence and a sufficiently small heating interval remain
//! separate physical requirements.

use crate::{
    cavity_phs::piston::thermal::{HeatExchangeSolve, feedback},
    finite_gas_thermal::FiniteGasThermalEnvironment,
};
use fs_bem::observation_transfers::{
    HarmonicVolumeQuadrature, MolecularObservationError, PreparedMolecularAbsorption,
    PreparedMolecularResponse,
};
use fs_conduction::lumped::LumpedThermalEnvironment;
use fs_material::{
    MaterialError,
    relaxation::{RelaxationAbsorption, RelaxationState, RelaxingMedium},
    thermochemical_gas::energy::{
        GasEnergyError, GasEnergySolve, GasHeatStep, IsochoricGasParcel,
    },
};
use fs_math::c64::C64;

/// Source rates in the original bank order, at its one physical frequency.
/// Coherence factors are mode-major: each consecutive source-length row is
/// propagated coherently, then independent rows contribute positive heat.
#[derive(Clone, Copy, Debug)]
pub enum MolecularGasExcitation<'a> {
    /// Peak deterministic source-rate phasors.
    Coherent(&'a [C64]),
    /// Peak second moment C = B B*, with independent coherent modes.
    CoherenceModes { factors: &'a [C64], modes: usize },
    /// Nonnegative peak second moments for mutually independent sources.
    IndependentSources(&'a [f64]),
}

/// Refusal leaves the accepted gas and original acoustic source bank unchanged.
#[derive(Debug)]
pub enum MolecularHeatingError {
    Invalid(&'static str),
    Cancelled,
    Observation(MolecularObservationError),
    Material(MaterialError),
    Gas(GasEnergyError),
    /// Frozen-source heating exceeded the caller's admitted temperature change.
    TemperatureIncrement { actual_k: f64, maximum_k: f64 },
}

impl core::fmt::Display for MolecularHeatingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(reason) => write!(f, "molecular gas heating: {reason}"),
            Self::Cancelled => write!(f, "molecular gas heating cancelled"),
            Self::Observation(error) => error.fmt(f),
            Self::Material(error) => error.fmt(f),
            Self::Gas(error) => error.fmt(f),
            Self::TemperatureIncrement { actual_k, maximum_k } => write!(
                f,
                "molecular gas heating temperature increment {actual_k} K exceeds {maximum_k} K"
            ),
        }
    }
}
impl std::error::Error for MolecularHeatingError {}
impl From<MolecularObservationError> for MolecularHeatingError {
    fn from(error: MolecularObservationError) -> Self { Self::Observation(error) }
}
impl From<MaterialError> for MolecularHeatingError {
    fn from(error: MaterialError) -> Self { Self::Material(error) }
}
impl From<GasEnergyError> for MolecularHeatingError {
    fn from(error: GasEnergyError) -> Self { Self::Gas(error) }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), MolecularHeatingError> {
    if cancelled() { Err(MolecularHeatingError::Cancelled) } else { Ok(()) }
}

/// A finite, uniformly mixed parcel bound to its original acoustic readout.
///
/// The caller supplies the ORIGINAL thermochemical parcel and relaxation law.
/// The acoustic bank retains an evaluated relaxation snapshot, not the NASA-9
/// pack identity, so numerical snapshot equality alone is not a provenance
/// certificate. This binding retains the supplied complete parcel/model and
/// relaxation law immutably for every load and candidate.
///
/// The positive quadrature represents exactly the finite parcel volume.
/// No hidden volume rescaling, extra peak-to-RMS factor or additional
/// attenuation is applied; wall work and outgoing acoustic flux are excluded.
pub struct MolecularGasHeating<'a> {
    absorption: PreparedMolecularAbsorption<'a>,
    gas: &'a IsochoricGasParcel,
    relaxation: &'a RelaxationAbsorption,
}

impl<'a> MolecularGasHeating<'a> {
    /// Bind the original state, all evaluated molecular reservoirs and volume.
    /// Volume agreement admits floating-point roundoff only. The caller still
    /// owns the spatial region represented by its positive quadrature; a lumped
    /// gas parcel has no independent shape from which to infer that region.
    pub fn new(
        absorption: PreparedMolecularAbsorption<'a>,
        gas: &'a IsochoricGasParcel,
        relaxation: &'a RelaxationAbsorption,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, MolecularHeatingError> {
        poll(&mut cancelled)?;
        let source = relaxation.evaluate(&RelaxingMedium::from_thermochemical(gas.state())?)?;
        if &source != absorption.source() {
            return Err(MolecularHeatingError::Invalid(
                "acoustic source must match the original complete gas and relaxation state",
            ));
        }
        let ratio = absorption.volume_m3() / gas.volume_m3();
        if !ratio.is_finite() || (ratio - 1.0).abs() > 64.0 * f64::EPSILON {
            return Err(MolecularHeatingError::Invalid(
                "positive acoustic quadrature must represent the original gas parcel volume",
            ));
        }
        poll(&mut cancelled)?;
        Ok(Self { absorption, gas, relaxation })
    }

    /// Original accepted state and complete thermochemical source model.
    #[must_use]
    pub const fn gas(&self) -> &IsochoricGasParcel { self.gas }

    /// Same molecular law retained for the candidate's updated source state.
    #[must_use]
    pub const fn relaxation_law(&self) -> &RelaxationAbsorption { self.relaxation }

    /// Bound field, frequency, source order, receiver geometry and quadrature.
    #[must_use]
    pub fn absorption(&self) -> &PreparedMolecularAbsorption<'a> { &self.absorption }

    /// Evaluate mean molecular heat without changing the gas.
    ///
    /// This is a cycle-averaged steady heating rate over an explicit physical
    /// duration, not a finite-window instantaneous phasor-energy integral.
    /// Each coherence query uses the prepared owner's bounded workspace; no
    /// boundary solve, dense covariance or per-query field allocation occurs.
    /// Retaining this borrowed load prevents changing its prepared readout.
    pub fn heat_load(
        &mut self,
        excitation: MolecularGasExcitation<'_>,
        duration_s: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MolecularGasHeatLoad<'_>, MolecularHeatingError> {
        poll(&mut cancelled)?;
        if !duration_s.is_finite() || duration_s <= 0.0 {
            return Err(MolecularHeatingError::Invalid(
                "finite positive physical heating duration required",
            ));
        }
        let absorbed_power_w = match excitation {
            MolecularGasExcitation::Coherent(rates) => self.absorption
                .evaluate_coherent(rates, &mut cancelled)?.absorbed_power_w,
            MolecularGasExcitation::CoherenceModes { factors, modes } => self.absorption
                .evaluate_coherence_modes(factors, modes, &mut cancelled)?.absorbed_power_w,
            MolecularGasExcitation::IndependentSources(moments) => self.absorption
                .evaluate_independent_sources(moments, &mut cancelled)?.absorbed_power_w,
        };
        let absorbed_energy_j = absorbed_power_w * duration_s;
        if !absorbed_power_w.is_finite() || absorbed_power_w < 0.0
            || !absorbed_energy_j.is_finite()
            || (absorbed_power_w > 0.0 && absorbed_energy_j <= 0.0)
        {
            return Err(MolecularHeatingError::Invalid(
                "molecular absorbed power or its physical energy is unrepresentable",
            ));
        }
        poll(&mut cancelled)?;
        Ok(MolecularGasHeatLoad {
            gas: self.gas,
            relaxation: self.relaxation,
            response: self.absorption.response(),
            quadrature: self.absorption.quadrature(),
            volume_m3: self.absorption.volume_m3(),
            duration_s,
            absorbed_power_w,
            absorbed_energy_j,
        })
    }
}

/// One immutable, source-bound cycle-averaged load on the original parcel.
///
/// To combine absorption with reciprocal finite-body thermal exchange, pass
/// `gas()` and `absorbed_energy_j()` to the existing
/// `FiniteGasThermalEnvironment` as `gas` and `gas_heat_j`. That owner must
/// start from THIS original parcel. Do not first call `heat_gas` and then
/// inject the same energy into its already-heated candidate.
#[derive(Clone, Copy)]
pub struct MolecularGasHeatLoad<'a> {
    gas: &'a IsochoricGasParcel,
    relaxation: &'a RelaxationAbsorption,
    response: &'a PreparedMolecularResponse,
    quadrature: &'a [HarmonicVolumeQuadrature],
    volume_m3: f64,
    duration_s: f64,
    absorbed_power_w: f64,
    absorbed_energy_j: f64,
}

impl<'a> MolecularGasHeatLoad<'a> {
    #[must_use]
    pub const fn gas(&self) -> &'a IsochoricGasParcel { self.gas }
    #[must_use]
    pub const fn relaxation_law(&self) -> &'a RelaxationAbsorption { self.relaxation }
    /// Original source bank, including source/receiver identity and frequency.
    #[must_use]
    pub const fn response(&self) -> &'a PreparedMolecularResponse { self.response }
    /// Original receiver indices and physical volume weights, without remapping.
    #[must_use]
    pub const fn quadrature(&self) -> &'a [HarmonicVolumeQuadrature] { self.quadrature }
    #[must_use]
    pub const fn volume_m3(&self) -> f64 { self.volume_m3 }
    #[must_use]
    pub const fn duration_s(&self) -> f64 { self.duration_s }
    #[must_use]
    pub const fn absorbed_power_w(&self) -> f64 { self.absorbed_power_w }
    #[must_use]
    pub const fn absorbed_energy_j(&self) -> f64 { self.absorbed_energy_j }

    /// Bind this exact acoustic heat load into the existing reciprocal finite
    /// gas/body thermal solver without exposing a second heat-input seam.
    ///
    /// The environment always starts from this load's ORIGINAL accepted gas and
    /// supplies exactly `absorbed_energy_j()` as the gas-side external heat.
    /// Convection, body enthalpy, radiation and reciprocal solid-to-gas exchange
    /// remain owned by `FiniteGasThermalEnvironment`. This prevents callers
    /// from first heating the gas and then accidentally crediting the same
    /// molecular absorption again during body exchange.
    ///
    /// All numerical controls are admitted before the environment is published.
    /// The returned object is still a proposal owner: neither gas nor body is
    /// mutated until its own `advance` call succeeds and the caller accepts it.
    pub fn finite_body_environment<'b>(
        &'b self,
        body: &'b LumpedThermalEnvironment,
        gas_solve: GasEnergySolve,
        exchange: HeatExchangeSolve,
        maximum_gas_temperature_increment_k: f64,
    ) -> Result<FiniteGasThermalEnvironment<'b>, MolecularHeatingError> {
        gas_solve.validate()?;
        feedback::admit::<fs_conduction::ConductionError>(exchange)
            .map_err(|_| MolecularHeatingError::Invalid(
                "invalid reciprocal finite-body heat-exchange controls",
            ))?;
        if !maximum_gas_temperature_increment_k.is_finite()
            || maximum_gas_temperature_increment_k <= 0.0
        {
            return Err(MolecularHeatingError::Invalid(
                "finite positive maximum gas temperature increment required",
            ));
        }
        Ok(FiniteGasThermalEnvironment {
            gas: self.gas,
            body,
            gas_heat_j: self.absorbed_energy_j,
            gas_solve,
            exchange,
            maximum_gas_temperature_increment_k,
        })
    }

    /// Apply this load once to a candidate cloned from the ORIGINAL parcel.
    ///
    /// The native source-backed inversion owns u(T), Cv, fixed mass/density,
    /// EOS pressure and all transport properties. The next molecular state is
    /// re-evaluated from the SAME relaxation law and must stay inside its
    /// source domain. Any refusal or cancellation discards the whole candidate.
    /// The maximum temperature increment bounds frozen-field heating; it is
    /// not a certificate of acoustic/thermal operator-splitting accuracy.
    pub fn heat_gas(
        &self,
        controls: GasEnergySolve,
        maximum_temperature_increment_k: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MolecularGasHeatingStep<'a>, MolecularHeatingError> {
        poll(&mut cancelled)?;
        controls.validate()?;
        if !maximum_temperature_increment_k.is_finite()
            || maximum_temperature_increment_k <= 0.0
        {
            return Err(MolecularHeatingError::Invalid(
                "finite positive maximum gas temperature increment required",
            ));
        }
        let mut gas = (*self.gas).clone();
        let energy = gas.add_heat(self.absorbed_energy_j, controls, &mut cancelled)?;
        let increment = gas.state().state().temperature - self.gas.state().state().temperature;
        if !increment.is_finite() || increment < 0.0 {
            return Err(MolecularHeatingError::Invalid(
                "passive molecular heat must retain finite nondecreasing gas temperature",
            ));
        }
        if increment > maximum_temperature_increment_k {
            return Err(MolecularHeatingError::TemperatureIncrement {
                actual_k: increment,
                maximum_k: maximum_temperature_increment_k,
            });
        }
        poll(&mut cancelled)?;
        let next_relaxation = self.relaxation
            .evaluate(&RelaxingMedium::from_thermochemical(gas.state())?)?;
        poll(&mut cancelled)?;
        Ok(MolecularGasHeatingStep { load: *self, gas, next_relaxation, energy })
    }
}

/// Complete candidate and independently accounted original-law heat balance.
/// Accept its gas explicitly, then re-solve acoustics using `next_relaxation`
/// before a temperature-feedback step; the old harmonic field remains frozen.
pub struct MolecularGasHeatingStep<'a> {
    load: MolecularGasHeatLoad<'a>,
    gas: IsochoricGasParcel,
    next_relaxation: RelaxationState,
    energy: GasHeatStep,
}
impl<'a> MolecularGasHeatingStep<'a> {
    #[must_use]
    pub const fn load(&self) -> &MolecularGasHeatLoad<'a> { &self.load }
    #[must_use]
    pub const fn gas(&self) -> &IsochoricGasParcel { &self.gas }
    #[must_use]
    pub const fn next_relaxation(&self) -> &RelaxationState { &self.next_relaxation }
    #[must_use]
    pub const fn energy_balance(&self) -> &GasHeatStep { &self.energy }
    /// Accept the candidate without modifying the original parcel implicitly.
    #[must_use]
    pub fn into_gas(self) -> IsochoricGasParcel { self.gas }
}
