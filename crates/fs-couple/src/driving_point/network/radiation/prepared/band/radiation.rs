//! Each selected covariance must pass its own original angular convergence.
use super::*;
use fs_bem::observation_transfers::RadiatedPowerMoments;

/// Compensated sum over selected bins [W]. Angular convergence estimates stay
/// with `bin_power_moments()`; no fictitious broadband convergence certificate
/// or extra radiation-work debit is attached to this sum.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct RadiatingWelchBandPowerMoments {
    pub total_power_w: f64,
    pub first_directional_moment_w: [f64; 3],
    pub second_directional_moment_w: [[f64; 3]; 3],
}

pub struct PreparedRadiatingWelchBandFarField<'b, 'p, 'a> {
    band: PreparedRadiatingWelchBand<'b, 'p, 'a>,
    directions: &'p [[f64; 3]],
    directional: Vec<f64>,
    sums: Vec<Accumulator>,
    per_bin: Vec<RadiatedPowerMoments>,
    power: RadiatingWelchBandPowerMoments,
    valid: bool,
    workspace_bytes: usize,
}

impl<'b, 'p, 'a> PreparedRadiatingWelchBand<'b, 'p, 'a> {
    /// Require the original angular owner at EVERY selected bin, with identical
    /// requested unit directions and exterior fluid. Each frequency keeps its
    /// own cached quadrature ladder and tests its own measured covariance.
    /// Missing angular coverage is refused, never treated as zero power.
    pub fn far_field(
        &self,
        maximum_workspace_bytes: usize,
        gate: &CancelGate,
    ) -> Result<PreparedRadiatingWelchBandFarField<'b, 'p, 'a>, ExteriorNetworkError> {
        poll(gate)?;
        let owner: &'p PreparedRadiatingNetwork<'a> = self.bins[0].owner;
        let first = owner
            .radiation
            .as_ref()
            .ok_or_else(|| bad("band bin has no original angular preparation"))?;
        let directions = first.directions();
        let medium = first.medium();
        for bin in self.bins {
            poll(gate)?;
            let radiation = bin
                .owner
                .radiation
                .as_ref()
                .ok_or_else(|| bad("band bin has no original angular preparation"))?;
            let other = radiation.medium();
            if radiation.directions().len() != directions.len()
                || medium.density.to_bits() != other.density.to_bits()
                || medium.sound_speed.to_bits() != other.sound_speed.to_bits()
            {
                return Err(bad("band angular directions and physical fluid must match"));
            }
            for (index, (a, b)) in directions.iter().zip(radiation.directions()).enumerate() {
                if index % 256 == 0 {
                    poll(gate)?;
                }
                if a != b {
                    return Err(bad("band angular directions and order must match"));
                }
            }
        }
        let workspace_bytes = add(
            bytes::<PreparedRadiatingWelchBandFarField<'_, '_, '_>>(1)?,
            add(
                bytes::<RadiatedPowerMoments>(self.bins.len())?,
                add(
                    bytes::<Accumulator>(directions.len())?,
                    bytes::<f64>(directions.len())?,
                )?,
            )?,
        )?;
        storage(workspace_bytes, maximum_workspace_bytes)?;
        let result = PreparedRadiatingWelchBandFarField {
            band: *self,
            directions,
            directional: zeros(directions.len())?,
            sums: zeros(directions.len())?,
            per_bin: zeros(self.bins.len())?,
            power: RadiatingWelchBandPowerMoments::default(),
            valid: false,
            workspace_bytes,
        };
        poll(gate)?;
        Ok(result)
    }
}

impl PreparedRadiatingWelchBandFarField<'_, '_, '_> {
    pub fn band(&self) -> &PreparedRadiatingWelchBand<'_, '_, '_> {
        &self.band
    }
    pub fn directions(&self) -> &[[f64; 3]] {
        self.directions
    }
    pub fn workspace_bytes(&self) -> usize {
        self.workspace_bytes
    }
    pub fn is_valid(&self) -> bool {
        self.valid
    }
    pub fn directional_power_w_per_sr(&self) -> Option<&[f64]> {
        self.valid.then_some(&self.directional)
    }
    pub fn power_moments(&self) -> Option<RadiatingWelchBandPowerMoments> {
        self.valid.then_some(self.power)
    }
    /// Original convergence evidence for EACH selected physical frequency.
    pub fn bin_power_moments(&self) -> Option<&[RadiatedPowerMoments]> {
        self.valid.then_some(&self.per_bin)
    }

    /// Allocation-free stationary sum. Whole-band work conservatively counts
    /// bins*(original per-bin angular allowance + requested directions + 13).
    /// Any bin's convergence/budget/finite/cancellation refusal invalidates ALL
    /// band getters, including after a previous successful evaluation.
    pub fn evaluate(
        &mut self,
        per_bin_budget: HarmonicWorkBudget,
        band_budget: HarmonicWorkBudget,
        gate: &CancelGate,
    ) -> Result<RadiatingWelchBandPowerMoments, ExteriorNetworkError> {
        self.valid = false;
        poll(gate)?;
        work(
            self.band.bins.len(),
            per_bin_budget.maximum_products,
            add(self.directions.len(), 13)?,
            band_budget,
        )?;
        self.sums.fill(Accumulator::default());
        let mut power = [Accumulator::default(); 13];
        for (index, bin) in self.band.bins.iter().enumerate() {
            poll(gate)?;
            let far = bin
                .owner
                .far_field()
                .ok_or_else(|| bad("band bin has no original angular preparation"))?;
            let value = far.evaluate_welch_bin(bin, per_bin_budget, &mut self.directional, gate)?;
            self.per_bin[index] = value;
            power[0].add(C64::from_re(value.total_power_w));
            for i in 0..3 {
                power[1 + i].add(C64::from_re(value.first_directional_moment_w[i]));
                for j in 0..3 {
                    power[4 + 3 * i + j].add(C64::from_re(value.second_directional_moment_w[i][j]));
                }
            }
            for (index, (sum, &value)) in self.sums.iter_mut().zip(&self.directional).enumerate() {
                if index % 256 == 0 {
                    poll(gate)?;
                }
                sum.add(C64::from_re(value));
            }
        }
        self.power.total_power_w = power[0].finish()?.re;
        for i in 0..3 {
            self.power.first_directional_moment_w[i] = power[1 + i].finish()?.re;
            for j in 0..3 {
                self.power.second_directional_moment_w[i][j] = power[4 + 3 * i + j].finish()?.re;
            }
        }
        for (index, (out, sum)) in self.directional.iter_mut().zip(&self.sums).enumerate() {
            if index % 256 == 0 {
                poll(gate)?;
            }
            *out = sum.finish()?.re;
        }
        poll(gate)?;
        self.valid = true;
        Ok(self.power)
    }
}
