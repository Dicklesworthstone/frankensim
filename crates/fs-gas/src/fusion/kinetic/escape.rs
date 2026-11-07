//! Prescribed energy-selective escape from an isotropic alpha population.
//!
//! This is the local loss term `dn(E)/dt = -lambda(E) n(E)`. The caller
//! supplies its physical-energy dependence; no magnetic orbit, pitch-angle,
//! wall-interception or ambipolar-potential calculation is inferred.

use super::{AlphaEnergyGrid, MAX_ALPHA_BINS};
use crate::fusion::solve::Sum;
use crate::fusion::{DtFusionError, poll, positive};
use fs_math::decay::{attenuate, decayed_fraction};
use fs_qty::{Energy, Frequency};
use std::sync::Arc;

const AVOGADRO: f64 = 6.022_140_76e23;
const MAX_KNOTS: usize = 4096;

/// Identified, piecewise-linear escape hazard against physical alpha energy.
///
/// Energies are joules per particle, hazards are nonnegative inverse seconds.
/// The original ordered table is retained without sorting, extrapolation or
/// clamping. Every represented alpha center must remain inside its support,
/// including during prescribed compression/expansion. Zero hazard is allowed.
/// Rates are an authored confinement closure, not a prediction of orbit loss.
#[derive(Clone, Debug, PartialEq)]
pub struct AlphaEscapeTable {
    data_id: String,
    energies: Arc<[Energy]>,
    hazards: Arc<[Frequency]>,
}

impl AlphaEscapeTable {
    /// Admit 2..=4096 strictly increasing positive energies and matching rates.
    pub fn new(
        data_id: String,
        energies: Vec<Energy>,
        hazards: Vec<Frequency>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, DtFusionError> {
        poll(&mut cancelled)?;
        if data_id.trim().is_empty()
            || !(2..=MAX_KNOTS).contains(&energies.len())
            || energies.len() != hazards.len()
        {
            return Err(DtFusionError::Invalid(
                "identified alpha escape table with 2..=4096 matching energy/rate entries required",
            ));
        }
        for (i, (energy, hazard)) in energies.iter().zip(&hazards).enumerate() {
            poll(&mut cancelled)?;
            if !positive(energy.value())
                || (i > 0 && *energy <= energies[i - 1])
                || !hazard.value().is_finite()
                || hazard.value() < 0.0
            {
                return Err(DtFusionError::Invalid(
                    "ordered positive alpha energies and finite nonnegative escape hazards required",
                ));
            }
        }
        poll(&mut cancelled)?;
        Ok(Self {
            data_id,
            energies: energies.into(),
            hazards: hazards.into(),
        })
    }

    pub fn data_id(&self) -> &str {
        &self.data_id
    }
    pub fn energies(&self) -> &[Energy] {
        &self.energies
    }
    pub fn hazards(&self) -> &[Frequency] {
        &self.hazards
    }

    /// Evaluate the declared linear interpolant inside the original support.
    pub fn hazard_at(&self, energy: Energy) -> Result<Frequency, DtFusionError> {
        if !positive(energy.value())
            || energy < self.energies[0]
            || energy > self.energies[self.energies.len() - 1]
        {
            return Err(DtFusionError::OutsideDomain);
        }
        let upper = self.energies.partition_point(|e| *e < energy);
        if self.energies[upper] == energy {
            return Ok(self.hazards[upper]);
        }
        let lower = upper - 1;
        let weight = (energy.value() - self.energies[lower].value())
            / (self.energies[upper].value() - self.energies[lower].value());
        let value =
            (1.0 - weight) * self.hazards[lower].value() + weight * self.hazards[upper].value();
        if !value.is_finite() || value < 0.0 {
            return Err(DtFusionError::Arithmetic("alpha escape interpolation"));
        }
        Ok(Frequency::new(value))
    }

    pub(crate) fn prepare(
        &self,
        grid: &AlphaEnergyGrid,
        scale: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedAlphaEscape, DtFusionError> {
        let mut rates_s_inv = [0.0; MAX_ALPHA_BINS];
        let mut maximum_rate_s_inv = 0.0_f64;
        for (rate, energy) in rates_s_inv.iter_mut().zip(grid.energies()) {
            poll(cancelled)?;
            *rate = self.hazard_at(Energy::new(energy.value() * scale))?.value();
            maximum_rate_s_inv = maximum_rate_s_inv.max(*rate);
        }
        Ok(PreparedAlphaEscape {
            rates_s_inv,
            maximum_rate_s_inv,
        })
    }
}

/// Stack-owned coefficients for one physical energy coordinate.
pub(crate) struct PreparedAlphaEscape {
    rates_s_inv: [f64; MAX_ALPHA_BINS],
    maximum_rate_s_inv: f64,
}

pub(crate) struct AlphaEscapeStep {
    pub(crate) amounts: [f64; MAX_ALPHA_BINS],
    pub(crate) retained_mol: f64,
    pub(crate) removed_mol: f64,
    pub(crate) removed_energy_j: f64,
}

impl PreparedAlphaEscape {
    pub(crate) fn maximum_rate_s_inv(&self) -> f64 {
        self.maximum_rate_s_inv
    }

    /// Separate positive survivors and expm1 losses retain both deep tails and
    /// weak transfers. Energy is removed at each actual physical bin center.
    pub(crate) fn advance(
        &self,
        grid: &AlphaEnergyGrid,
        scale: f64,
        old: &[f64; MAX_ALPHA_BINS],
        dt: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<AlphaEscapeStep, DtFusionError> {
        poll(cancelled)?;
        let mut amounts = *old;
        let mut retained = Sum::default();
        let mut removed = Sum::default();
        let mut energy = Sum::default();
        for (i, center) in grid.energies().iter().enumerate() {
            poll(cancelled)?;
            let exposure = self.rates_s_inv[i] * dt;
            let loss = old[i] * decayed_fraction(exposure);
            amounts[i] = attenuate(old[i], exposure);
            if !loss.is_finite()
                || !amounts[i].is_finite()
                || loss < 0.0
                || amounts[i] < 0.0
                || (old[i] > 0.0 && self.rates_s_inv[i] > 0.0 && dt > 0.0 && exposure == 0.0)
            {
                return Err(DtFusionError::Arithmetic("alpha escape exposure"));
            }
            retained.add(amounts[i]);
            removed.add(loss);
            energy.add((loss * AVOGADRO) * (center.value() * scale));
        }
        if [retained.value, removed.value, energy.value]
            .iter()
            .any(|x| !x.is_finite() || *x < 0.0)
        {
            return Err(DtFusionError::Arithmetic("alpha escape moments"));
        }
        poll(cancelled)?;
        Ok(AlphaEscapeStep {
            amounts,
            retained_mol: retained.value,
            removed_mol: removed.value,
            removed_energy_j: energy.value,
        })
    }
}
