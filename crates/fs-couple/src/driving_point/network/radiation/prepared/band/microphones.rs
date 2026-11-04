//! Stationary array moments of the original frequency-dependent p/u readout.
use super::*;

pub struct PreparedRadiatingWelchBandMicrophones<'b, 'p, 'a> {
    band: PreparedRadiatingWelchBand<'b, 'p, 'a>,
    microphones: Vec<PreparedRadiatingMicrophones<'p, 'a>>,
    patterns: Vec<FirstOrder>,
    sums: Vec<Accumulator>,
    spectrum: Vec<C64>,
    valid: bool,
    workspace_bytes: usize,
}

impl<'b, 'p, 'a> PreparedRadiatingWelchBand<'b, 'p, 'a> {
    /// Project the same authored first-order pattern at every original receiver
    /// in every frequency bin. Each original near-field p/u operator supplies
    /// its own physical delay, complex response and exterior medium. This does
    /// not replace a directional microphone by pressure or a plane wave.
    pub fn microphones(
        &self,
        patterns: &[FirstOrder],
        maximum_workspace_bytes: usize,
        gate: &CancelGate,
    ) -> Result<PreparedRadiatingWelchBandMicrophones<'b, 'p, 'a>, ExteriorNetworkError> {
        poll(gate)?;
        let r = self.receivers().len();
        if patterns.len() != r {
            return Err(bad(
                "band needs one common microphone pattern per original receiver",
            ));
        }
        let entries = mul(r, r)?;
        let mut workspace_bytes = add(
            bytes::<PreparedRadiatingWelchBandMicrophones<'_, '_, '_>>(1)?,
            add(
                bytes::<FirstOrder>(r)?,
                add(bytes::<Accumulator>(entries)?, bytes::<C64>(entries)?)?,
            )?,
        )?;
        for bin in self.bins {
            poll(gate)?;
            // Original projected transfer, complex Neumaier sums, and CSD.
            let original = bytes::<C64>(add(
                mul(bin.owner.observer.sources().len(), r)?,
                mul(3, entries)?,
            )?)?;
            workspace_bytes = add(
                workspace_bytes,
                add(bytes::<PreparedRadiatingMicrophones<'_, '_>>(1)?, original)?,
            )?;
        }
        storage(workspace_bytes, maximum_workspace_bytes)?;
        let mut microphones = reserve(self.bins.len())?;
        for bin in self.bins {
            poll(gate)?;
            microphones.push(bin.owner.microphones(patterns)?);
        }
        let mut retained_patterns = reserve(r)?;
        retained_patterns.extend_from_slice(patterns);
        let result = PreparedRadiatingWelchBandMicrophones {
            band: *self,
            microphones,
            patterns: retained_patterns,
            sums: zeros(entries)?,
            spectrum: zeros(entries)?,
            valid: false,
            workspace_bytes,
        };
        poll(gate)?;
        Ok(result)
    }
}

impl PreparedRadiatingWelchBandMicrophones<'_, '_, '_> {
    pub fn band(&self) -> &PreparedRadiatingWelchBand<'_, '_, '_> {
        &self.band
    }
    pub fn patterns(&self) -> &[FirstOrder] {
        &self.patterns
    }
    pub fn workspace_bytes(&self) -> usize {
        self.workspace_bytes
    }
    pub fn is_valid(&self) -> bool {
        self.valid
    }
    /// Sum of selected-bin CSD matrices [Pa²], row-major in original receiver
    /// order. Complex off-diagonal phase is retained; this is not a waveform.
    pub fn cross_moments_pa2(&self) -> Option<&[C64]> {
        self.valid.then_some(&self.spectrum)
    }

    /// No allocation or repeated source-to-opening composition. The whole-band
    /// allowance counts bins*(original per-bin allowance + microphones²).
    pub fn evaluate(
        &mut self,
        per_bin_budget: HarmonicWorkBudget,
        band_budget: HarmonicWorkBudget,
        gate: &CancelGate,
    ) -> Result<&[C64], ExteriorNetworkError> {
        self.valid = false;
        poll(gate)?;
        work(
            self.band.bins.len(),
            per_bin_budget.maximum_products,
            self.spectrum.len(),
            band_budget,
        )?;
        self.sums.fill(Accumulator::default());
        for (bin, microphones) in self.band.bins.iter().zip(&mut self.microphones) {
            poll(gate)?;
            let values = microphones.welch_bin(bin, per_bin_budget, gate)?;
            for (index, (sum, &value)) in self.sums.iter_mut().zip(values).enumerate() {
                if index % 256 == 0 {
                    poll(gate)?;
                }
                sum.add(value);
            }
        }
        for (index, (out, sum)) in self.spectrum.iter_mut().zip(&self.sums).enumerate() {
            if index % 256 == 0 {
                poll(gate)?;
            }
            *out = sum.finish()?;
        }
        poll(gate)?;
        self.valid = true;
        Ok(&self.spectrum)
    }
}
