//! Actual stationary bath entropy on the original SSP boundary transfers.
//! The state account uses the same fixed species reference as the bath account.
//! Chemical/enthalpy terms are combined before source reference offsets return.

use super::{Dd, FlowError, MAX_SPECIES, MixtureGas, R, log_ratio, poll, scalar};
use super::super::evolve::Stage;
use super::super::reservoir::ReservoirEnds;
use fs_qty::Temperature;

#[derive(Clone, Copy)]
pub(super) struct EntropyExchange {
    pub input: Dd,
    pub source_input: [Dd; 2],
    pub species: [Dd; MAX_SPECIES],
    pub throughput: [Dd; MAX_SPECIES],
}

impl Default for EntropyExchange {
    fn default() -> Self {
        Self {
            input: Dd::ZERO,
            source_input: [Dd::ZERO; 2],
            species: [Dd::ZERO; MAX_SPECIES],
            throughput: [Dd::ZERO; MAX_SPECIES],
        }
    }
}

struct Bath {
    temperature: f64,
    positive_species: [bool; MAX_SPECIES],
    // Referenced caloric chemical potential / T, with the state basis removed.
    normalized_mu_over_t: [f64; MAX_SPECIES],
}

pub(super) struct ReservoirEntropyAccount {
    baths: [Option<Bath>; 2],
    molar_mass: [f64; MAX_SPECIES],
    reference_entropy: [f64; MAX_SPECIES],
    pub trial: EntropyExchange,
    pub total: EntropyExchange,
}

impl ReservoirEntropyAccount {
    pub fn new(
        gas: &MixtureGas,
        reference_temperature: f64,
        reference_density: &[f64; MAX_SPECIES],
        reference_entropy: &[f64; MAX_SPECIES],
        reservoirs: ReservoirEnds<'_>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, FlowError> {
        let mut molar_mass = [0.0; MAX_SPECIES];
        for (k, model) in gas.source().models().iter().enumerate() {
            molar_mass[k] = model.molar_mass().value();
            if !molar_mass[k].is_finite() {
                return Err(FlowError::EntropyDomain("reservoir molar mass"));
            }
        }
        let mut baths: [Option<Bath>; 2] = [None, None];
        for (end, trace) in reservoirs.into_iter().enumerate() {
            let Some(trace) = trace else { continue };
            let t = trace.primitive.temperature;
            let mut bath = Bath {
                temperature: t,
                positive_species: [false; MAX_SPECIES],
                normalized_mu_over_t: [0.0; MAX_SPECIES],
            };
            let log_t = log_ratio(t, reference_temperature);
            for (k, model) in gas.source().models().iter().enumerate() {
                poll(cancelled)?;
                let density = trace.cell.partial_densities[k];
                if density == 0.0 { continue; }
                bath.positive_species[k] = true;
                let sensible = model.molar_internal_energy_change(
                    Temperature::new(gas.domain().minimum_temperature),
                    Temperature::new(t),
                ).map_err(super::super::FlowError::from)?.value();
                let ds = model.molar_entropy_change(
                    Temperature::new(reference_temperature),
                    Temperature::new(t),
                ).map_err(super::super::FlowError::from)?.value();
                // s_b-alpha = ds - R ln(Tb/T*) - R ln(rho_b/rho*) + R.
                // href/Tb = sensible/Tb + R; those two R terms cancel exactly.
                bath.normalized_mu_over_t[k] = scalar(
                    Dd::from_f64(sensible) / Dd::from_f64(t) - Dd::from_f64(ds)
                        + Dd::from_f64(R) * Dd::from_f64(log_t)
                        + Dd::from_f64(R)
                            * Dd::from_f64(log_ratio(density, reference_density[k])),
                )?;
            }
            baths[end] = Some(bath);
        }
        Ok(Self {
            baths,
            molar_mass,
            reference_entropy: *reference_entropy,
            trial: EntropyExchange::default(),
            total: EntropyExchange::default(),
        })
    }

    /// Exact native signed left/right transfer arithmetic, with no source query.
    pub fn stage_pair(
        &mut self,
        a: &Stage,
        b: &Stage,
        dt: f64,
        areas: [f64; 2],
        n: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), FlowError> {
        let mut trial = EntropyExchange::default();
        for end in 0..2 {
            poll(cancelled)?;
            let (first, second, sign) = if end == 0 {
                (&a.left, &b.left, 1.0)
            } else {
                (&a.right, &b.right, -1.0)
            };
            let energy = sign * 0.5 * dt * areas[end] * (first[n + 1] + second[n + 1]);
            let mut input = match &self.baths[end] {
                Some(bath) => Dd::from_f64(energy) / Dd::from_f64(bath.temperature),
                None => {
                    if first[n + 1] != 0.0 || second[n + 1] != 0.0 {
                        return Err(FlowError::EntropyDomain("unowned molecular boundary heat"));
                    }
                    Dd::ZERO
                }
            };
            let mut reference = Dd::ZERO;
            for k in 0..n {
                let transfer = sign * 0.5 * dt * areas[end] * (first[k] + second[k]);
                let mass = Dd::from_f64(transfer);
                trial.species[k] = trial.species[k] + mass;
                // Each stage contributes absolute throughput; opposite ends or
                // reversing stage fluxes cannot hide the inventory scale.
                let throughput = 0.5 * dt * areas[end] * (first[k].abs() + second[k].abs());
                trial.throughput[k] = trial.throughput[k] + Dd::from_f64(throughput);
                match &self.baths[end] {
                    Some(bath) => {
                        if !bath.positive_species[k] && (first[k] != 0.0 || second[k] != 0.0) {
                            return Err(FlowError::EntropyDomain(
                                "nonzero molecular transfer at zero bath species activity",
                            ));
                        }
                        let moles = mass / Dd::from_f64(self.molar_mass[k]);
                        input = input - moles * Dd::from_f64(bath.normalized_mu_over_t[k]);
                        reference = reference
                            + moles * Dd::from_f64(self.reference_entropy[k]);
                    }
                    None if first[k] != 0.0 || second[k] != 0.0 => {
                        return Err(FlowError::EntropyDomain("unowned molecular boundary species"));
                    }
                    None => {}
                }
                scalar(trial.species[k])?;
                scalar(trial.throughput[k])?;
            }
            trial.input = trial.input + input;
            trial.source_input[end] = input + reference;
            scalar(trial.source_input[end])?;
        }
        scalar(trial.input)?;
        self.trial = trial;
        Ok(())
    }

    /// Start one native multistage quadrature attempt. Rejected attempts never
    /// contribute to accepted throughput or entropy input.
    pub fn begin_quadrature(&mut self) {
        self.trial = EntropyExchange::default();
    }

    /// Account for each actual RHS contribution with its SIGNED tableau weight.
    /// The final conservative transfer is retained by the native recurrence;
    /// only absolute stage throughput and zero-activity admission are needed here.
    pub fn observe_quadrature_stage(
        &mut self,
        stage: &Stage,
        weighted_dt: f64,
        areas: [f64; 2],
        n: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), FlowError> {
        if !weighted_dt.is_finite() {
            return Err(FlowError::EntropyDomain("non-finite molecular quadrature weight"));
        }
        if weighted_dt == 0.0 {
            return Ok(());
        }
        for end in 0..2 {
            poll(cancelled)?;
            let flux = if end == 0 { &stage.left } else { &stage.right };
            if self.baths[end].is_none() && flux[n + 1] != 0.0 {
                return Err(FlowError::EntropyDomain("unowned molecular boundary heat"));
            }
            let measure = Dd::from_f64(weighted_dt.abs()) * Dd::from_f64(areas[end]);
            for k in 0..n {
                match &self.baths[end] {
                    Some(bath) if !bath.positive_species[k] && flux[k] != 0.0 => {
                        return Err(FlowError::EntropyDomain(
                            "nonzero molecular transfer at zero bath species activity",
                        ));
                    }
                    None if flux[k] != 0.0 => {
                        return Err(FlowError::EntropyDomain("unowned molecular boundary species"));
                    }
                    _ => {}
                }
                self.trial.throughput[k] = self.trial.throughput[k]
                    + measure * Dd::from_f64(flux[k].abs());
                scalar(self.trial.throughput[k])?;
            }
        }
        Ok(())
    }

    /// Use the actual final native signed end transfers, not a separately
    /// rounded reconstruction of them. Stationary bath chemical efforts make
    /// entropy a linear function of these SAME energy/species integrals.
    pub fn finish_quadrature(
        &mut self,
        transfers: &[[f64; super::super::MAX_COMPONENTS]; 2],
        n: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), FlowError> {
        let mut trial = EntropyExchange {
            throughput: self.trial.throughput,
            ..EntropyExchange::default()
        };
        for end in 0..2 {
            poll(cancelled)?;
            let energy = transfers[end][n + 1];
            let mut input = match &self.baths[end] {
                Some(bath) => Dd::from_f64(energy) / Dd::from_f64(bath.temperature),
                None if energy != 0.0 => {
                    return Err(FlowError::EntropyDomain("unowned molecular boundary heat"));
                }
                None => Dd::ZERO,
            };
            let mut reference = Dd::ZERO;
            for k in 0..n {
                let transfer = transfers[end][k];
                let mass = Dd::from_f64(transfer);
                trial.species[k] = trial.species[k] + mass;
                match &self.baths[end] {
                    Some(bath) => {
                        if !bath.positive_species[k] && transfer != 0.0 {
                            return Err(FlowError::EntropyDomain(
                                "nonzero molecular transfer at zero bath species activity",
                            ));
                        }
                        let moles = mass / Dd::from_f64(self.molar_mass[k]);
                        input = input - moles * Dd::from_f64(bath.normalized_mu_over_t[k]);
                        reference = reference + moles * Dd::from_f64(self.reference_entropy[k]);
                    }
                    None if transfer != 0.0 => {
                        return Err(FlowError::EntropyDomain("unowned molecular boundary species"));
                    }
                    None => {}
                }
                scalar(trial.species[k])?;
            }
            trial.input = trial.input + input;
            trial.source_input[end] = input + reference;
            scalar(trial.source_input[end])?;
        }
        scalar(trial.input)?;
        self.trial = trial;
        Ok(())
    }

    pub fn accept(&mut self, n: usize) -> Result<(), FlowError> {
        let mut total = self.total;
        total.input = total.input + self.trial.input;
        scalar(total.input)?;
        for end in 0..2 {
            total.source_input[end] = total.source_input[end] + self.trial.source_input[end];
            scalar(total.source_input[end])?;
        }
        for k in 0..n {
            total.species[k] = total.species[k] + self.trial.species[k];
            total.throughput[k] = total.throughput[k] + self.trial.throughput[k];
            scalar(total.species[k])?;
            scalar(total.throughput[k])?;
        }
        self.total = total;
        Ok(())
    }
}
