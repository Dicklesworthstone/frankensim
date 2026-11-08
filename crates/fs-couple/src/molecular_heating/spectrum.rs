//! Stationary molecular absorption at distinct frequencies on one finite gas.
//!
//! Sources interfere coherently within each frequency; different frequencies
//! contribute mean powers, never summed complex amplitudes. The result is an
//! ensemble or sufficiently long-time average, not a finite-window beat-energy
//! integral. Each excitation uses the original peak convention and already
//! includes any spectral quadrature weight. Spectral/spatial convergence,
//! acoustic startup/storage and finite-time beats are not inferred here.
//!
//! Every bank retains its original molecular law and physical volume rule.
//! Their total heat enters the original nonlinear fixed-volume gas storage
//! once. Neither classical viscous/conductive heat nor outgoing radiation is
//! inferred from the molecular pressure loss.

use super::*;
use fs_math::dd::Dd;

/// Hard ceiling on the number of independently resolved frequencies.
/// Each prepared readout also retains its original source/mode/work limits.
pub const MAX_MOLECULAR_SPECTRUM_BINS: usize = 4096;

/// One original frequency's contribution to the shared thermal load.
/// This is a diagnostic power, not an additional thermal energy source.
#[derive(Clone, Copy)]
pub struct MolecularGasSpectrumBinPower<'a> {
    response: &'a PreparedMolecularResponse,
    quadrature: &'a [HarmonicVolumeQuadrature],
    absorbed_power_w: f64,
}

impl<'a> MolecularGasSpectrumBinPower<'a> {
    /// Original source order, frequency, medium and receiver geometry.
    #[must_use]
    pub const fn response(&self) -> &'a PreparedMolecularResponse { self.response }
    #[must_use]
    pub const fn quadrature(&self) -> &'a [HarmonicVolumeQuadrature] { self.quadrature }
    #[must_use]
    pub const fn absorbed_power_w(&self) -> f64 { self.absorbed_power_w }
}

/// Reusable bounded absorption readouts sharing the ORIGINAL parcel and law.
///
/// Bin and excitation order are the caller's original order. A separately
/// prepared frequency ordering makes the compensated power sum deterministic
/// under a permutation of bins. No dense covariance or acoustic boundary solve
/// is introduced, and heat queries allocate no per-frequency field storage.
pub struct MolecularGasSpectrumHeating<'a> {
    absorptions: Vec<PreparedMolecularAbsorption<'a>>,
    powers: Vec<MolecularGasSpectrumBinPower<'a>>,
    frequency_order: Vec<usize>,
    gas: &'a IsochoricGasParcel,
    relaxation: &'a RelaxationAbsorption,
}

impl<'a> MolecularGasSpectrumHeating<'a> {
    /// Bind each frequency to the same complete gas/relaxation snapshot and
    /// parcel volume. Different positive spatial rules may resolve that same
    /// region; no shape or quadrature accuracy is inferred from volume alone.
    /// As with [`MolecularGasHeating`], the caller supplies the original source
    /// objects; evaluated snapshot equality is not a provenance certificate.
    ///
    /// `maximum_bins` is a positive caller limit at most the hard ceiling.
    /// It bounds the aggregate work multiplier and extra O(bins) storage before
    /// allocation; each input readout already owns its admitted workspace and
    /// per-query source/mode/product limits. Duplicate frequencies refuse:
    /// combine their sources in one coherent bank or coherence-factor query.
    pub fn new(
        absorptions: Vec<PreparedMolecularAbsorption<'a>>,
        gas: &'a IsochoricGasParcel,
        relaxation: &'a RelaxationAbsorption,
        maximum_bins: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, MolecularHeatingError> {
        poll(&mut cancelled)?;
        if maximum_bins == 0 || maximum_bins > MAX_MOLECULAR_SPECTRUM_BINS
            || absorptions.is_empty() || absorptions.len() > maximum_bins
        {
            return Err(MolecularHeatingError::Invalid(
                "molecular spectrum requires a nonempty bounded frequency bank",
            ));
        }
        let source = relaxation.evaluate(&RelaxingMedium::from_thermochemical(gas.state())?)?;
        for absorption in &absorptions {
            poll(&mut cancelled)?;
            admit_binding(absorption, gas, &source)?;
        }
        let mut frequency_order = Vec::new();
        frequency_order.try_reserve_exact(absorptions.len()).map_err(|_| {
            MolecularHeatingError::Invalid("molecular spectrum ordering allocation failed")
        })?;
        frequency_order.extend(0..absorptions.len());
        frequency_order.sort_unstable_by(|&left, &right| {
            absorptions[left].angular_frequency_rad_s()
                .total_cmp(&absorptions[right].angular_frequency_rad_s())
        });
        for pair in frequency_order.windows(2) {
            poll(&mut cancelled)?;
            if absorptions[pair[0]].angular_frequency_rad_s()
                == absorptions[pair[1]].angular_frequency_rad_s()
            {
                return Err(MolecularHeatingError::Invalid(
                    "same-frequency sources must share one coherent or covariance query",
                ));
            }
        }
        let mut powers = Vec::new();
        powers.try_reserve_exact(absorptions.len()).map_err(|_| {
            MolecularHeatingError::Invalid("molecular spectrum power allocation failed")
        })?;
        for absorption in &absorptions {
            poll(&mut cancelled)?;
            powers.push(MolecularGasSpectrumBinPower {
                response: absorption.response(),
                quadrature: absorption.quadrature(),
                absorbed_power_w: 0.0,
            });
        }
        poll(&mut cancelled)?;
        Ok(Self { absorptions, powers, frequency_order, gas, relaxation })
    }

    #[must_use]
    pub const fn gas(&self) -> &'a IsochoricGasParcel { self.gas }
    #[must_use]
    pub const fn relaxation_law(&self) -> &'a RelaxationAbsorption { self.relaxation }
    /// Original readouts in the same order required for excitation queries.
    #[must_use]
    pub fn absorptions(&self) -> &[PreparedMolecularAbsorption<'a>] { &self.absorptions }

    /// Evaluate all bins before publishing a shared immutable heat load.
    ///
    /// Each bin preserves its own coherent/low-rank/independent source model.
    /// Powers add in frequency order using double-double accumulation; the
    /// total converts to joules exactly once. A retained load borrows this
    /// cache, so its original readouts cannot change underneath it. Refusal or
    /// cancellation publishes no load, leaves the gas unchanged, and permits
    /// a complete retry using the same prepared readouts.
    pub fn heat_load(
        &mut self,
        excitations: &[MolecularGasExcitation<'_>],
        duration_s: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MolecularGasSpectrumHeatLoad<'_>, MolecularHeatingError> {
        poll(&mut cancelled)?;
        admit_duration(duration_s)?;
        if excitations.len() != self.absorptions.len() {
            return Err(MolecularHeatingError::Invalid(
                "molecular spectrum requires one excitation per original frequency bank",
            ));
        }
        for ((absorption, power), &excitation) in self.absorptions.iter_mut()
            .zip(&mut self.powers).zip(excitations)
        {
            poll(&mut cancelled)?;
            power.absorbed_power_w = absorbed_power(absorption, excitation, &mut cancelled)?;
        }
        let mut total = Dd::ZERO;
        for &index in &self.frequency_order {
            poll(&mut cancelled)?;
            total = total + Dd::from_f64(self.powers[index].absorbed_power_w);
        }
        if !total.hi.is_finite() || !total.lo.is_finite() {
            return Err(MolecularHeatingError::Invalid(
                "total molecular spectrum power is unrepresentable",
            ));
        }
        let absorbed_power_w = total.to_f64();
        let absorbed_energy_j = absorbed_energy(absorbed_power_w, duration_s)?;
        poll(&mut cancelled)?;
        Ok(MolecularGasSpectrumHeatLoad {
            gas: self.gas,
            relaxation: self.relaxation,
            bins: &self.powers,
            duration_s,
            absorbed_power_w,
            absorbed_energy_j,
        })
    }
}

/// One stationary spectrum's mean heat on its original finite gas inventory.
/// Per-bin powers diagnose this total; they are not separate deposits.
#[derive(Clone, Copy)]
pub struct MolecularGasSpectrumHeatLoad<'a> {
    gas: &'a IsochoricGasParcel,
    relaxation: &'a RelaxationAbsorption,
    bins: &'a [MolecularGasSpectrumBinPower<'a>],
    duration_s: f64,
    absorbed_power_w: f64,
    absorbed_energy_j: f64,
}

impl<'a> MolecularGasSpectrumHeatLoad<'a> {
    #[must_use]
    pub const fn gas(&self) -> &'a IsochoricGasParcel { self.gas }
    #[must_use]
    pub const fn relaxation_law(&self) -> &'a RelaxationAbsorption { self.relaxation }
    /// Original frequency-bank order, unchanged by internal summation order.
    #[must_use]
    pub const fn bins(&self) -> &'a [MolecularGasSpectrumBinPower<'a>] { self.bins }
    #[must_use]
    pub fn volume_m3(&self) -> f64 { self.gas.volume_m3() }
    #[must_use]
    pub const fn duration_s(&self) -> f64 { self.duration_s }
    #[must_use]
    pub const fn absorbed_power_w(&self) -> f64 { self.absorbed_power_w }
    #[must_use]
    pub const fn absorbed_energy_j(&self) -> f64 { self.absorbed_energy_j }

    /// Bind the total once to the existing reciprocal body/gas thermal owner.
    /// It starts from this ORIGINAL gas, not an already-heated candidate.
    /// After acceptance, re-evaluate the original relaxation law before the
    /// next acoustic solve; the returned environment owns thermal exchange.
    pub fn finite_body_environment<'b>(
        &'b self,
        body: &'b LumpedThermalEnvironment,
        gas_solve: GasEnergySolve,
        exchange: HeatExchangeSolve,
        maximum_gas_temperature_increment_k: f64,
    ) -> Result<FiniteGasThermalEnvironment<'b>, MolecularHeatingError> {
        finite_body_environment(
            self.gas, self.absorbed_energy_j, body, gas_solve, exchange,
            maximum_gas_temperature_increment_k,
        )
    }

    /// One original nonlinear u(T) inversion and one source-law refresh for
    /// the entire spectrum. All frequencies remain frozen over this interval;
    /// a feedback step must re-solve them using the returned molecular state.
    pub fn heat_gas(
        &self,
        controls: GasEnergySolve,
        maximum_temperature_increment_k: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MolecularGasSpectrumHeatingStep<'a>, MolecularHeatingError> {
        let (gas, next_relaxation, energy) = heated_gas(
            self.gas, self.relaxation, self.absorbed_energy_j, controls,
            maximum_temperature_increment_k, &mut cancelled,
        )?;
        Ok(MolecularGasSpectrumHeatingStep { load: *self, gas, next_relaxation, energy })
    }
}

/// Shared candidate gas and independently accounted nonlinear storage change.
pub struct MolecularGasSpectrumHeatingStep<'a> {
    load: MolecularGasSpectrumHeatLoad<'a>,
    gas: IsochoricGasParcel,
    next_relaxation: RelaxationState,
    energy: GasHeatStep,
}

impl<'a> MolecularGasSpectrumHeatingStep<'a> {
    #[must_use]
    pub const fn load(&self) -> &MolecularGasSpectrumHeatLoad<'a> { &self.load }
    #[must_use]
    pub const fn gas(&self) -> &IsochoricGasParcel { &self.gas }
    #[must_use]
    pub const fn next_relaxation(&self) -> &RelaxationState { &self.next_relaxation }
    #[must_use]
    pub const fn energy_balance(&self) -> &GasHeatStep { &self.energy }
    /// Accept without implicitly changing the original parcel or source banks.
    #[must_use]
    pub fn into_gas(self) -> IsochoricGasParcel { self.gas }
}
