//! Fixed-entropy, fixed-volume chemical equilibrium (SV).
//!
//! On the same element-constrained equilibrium manifold as TV, the first law
//! gives `(dS/dT)_V = Cv_equilibrium/T`. This slope includes reaction and
//! ionization response. A prescribed chamber volume therefore determines its
//! temperature, pressure and composition without a frozen-gamma approximation.
//! The result is a reversible equilibrium endpoint; it does not model finite
//! reaction rates, shocks, irreversible work or entropy generation.

use super::*;

/// Temperature path and entropy residual for a fixed-volume isentropic solve.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EquilibriumVolumeEntropySolveV1 {
    /// Every temperature trial closes volume and chemical equilibrium.
    pub volume: EquilibriumVolumeSolveV1,
    /// Positive source-supported lower temperature [K].
    pub minimum_temperature_k: f64,
    /// Source-supported upper temperature, greater than the lower bound [K].
    pub maximum_temperature_k: f64,
    /// Absolute total entropy tolerance for this elemental inventory [J/K].
    pub entropy_tolerance_j_k: f64,
    /// At most 256 interior temperature trials, excluding bracket endpoints.
    pub maximum_temperature_iterations: usize,
}

impl EquilibriumVolumeEntropySolveV1 {
    fn validate(self) -> Result<(), EquilibriumVolumeErrorV1> {
        self.volume.validate()?;
        if !self.minimum_temperature_k.is_finite()
            || self.minimum_temperature_k <= 0.0
            || !self.maximum_temperature_k.is_finite()
            || self.maximum_temperature_k <= self.minimum_temperature_k
            || !self.entropy_tolerance_j_k.is_finite()
            || self.entropy_tolerance_j_k <= 0.0
            || !(1..=256).contains(&self.maximum_temperature_iterations)
        {
            return Err(EquilibriumVolumeErrorV1::Invalid(
                "invalid SV temperature bracket, entropy tolerance or budget",
            ));
        }
        Ok(())
    }
}

impl IdealGasEquilibriumV1 {
    /// Re-equilibrate at prescribed volume, total entropy and elemental amounts.
    ///
    /// Use an accepted equilibrium state's `entropy_j_k()` to compress or expand
    /// a closed chamber reversibly. All temperature trials use the original TV
    /// owner and source models, while analytic reacting `Cv/T` supplies the
    /// safeguarded Newton slope. Pressure, molecule count and charge are solved
    /// together; the feed specifies conserved elements rather than frozen species.
    ///
    /// Entropy and volume residuals must independently meet their declared
    /// tolerances. Original source validity and caloric continuity apply over the
    /// entire temperature bracket. One bounded workspace is reused; source
    /// receipts are allocated only after convergence, with no partial publication
    /// on cancellation or exhausted work.
    pub fn equilibrate_sv(
        &self,
        volume: Volume,
        feed: &[(SpeciesId, Amount)],
        target_entropy_j_k: f64,
        controls: EquilibriumVolumeEntropySolveV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<IdealGasIsochoricEquilibriumV1, EquilibriumVolumeErrorV1> {
        poll(&mut cancelled)?;
        controls.validate()?;
        if !target_entropy_j_k.is_finite() {
            return Err(EquilibriumVolumeErrorV1::Invalid(
                "SV target entropy must be finite",
            ));
        }
        let mut work = self.volume_workspace(volume, feed)?;
        work.admit_temperature_range(
            self,
            controls.minimum_temperature_k,
            controls.maximum_temperature_k,
        )?;
        let mut counts = Counts::default();
        let mut low = controls.minimum_temperature_k;
        let mut high = controls.maximum_temperature_k;
        let mut endpoints = [0.0; 2];
        for (index, t) in [low, high].into_iter().enumerate() {
            let temperature = Temperature::new(t);
            let volume_residual = work.solve_tv(
                self,
                temperature,
                volume,
                controls.volume,
                &mut counts,
                &mut cancelled,
            )?;
            endpoints[index] = work.entropy(self, temperature, &mut cancelled)?;
            if (endpoints[index] - target_entropy_j_k).abs() <= controls.entropy_tolerance_j_k {
                return work.publish_sv(
                    self,
                    temperature,
                    volume,
                    target_entropy_j_k,
                    controls,
                    counts,
                    0,
                    volume_residual,
                    &mut cancelled,
                );
            }
        }
        if !(endpoints[0] < target_entropy_j_k && target_entropy_j_k < endpoints[1]) {
            return Err(EquilibriumExpansionErrorV1::EntropyNotBracketed {
                lower_entropy_j_k: endpoints[0],
                upper_entropy_j_k: endpoints[1],
            }
            .into());
        }
        let mut t = low + 0.5 * (high - low);
        let mut residual = endpoints[1] - target_entropy_j_k;
        let mut steps = 0;
        for step in 1..=controls.maximum_temperature_iterations {
            steps = step;
            let temperature = Temperature::new(t);
            let volume_residual = work.solve_tv(
                self,
                temperature,
                volume,
                controls.volume,
                &mut counts,
                &mut cancelled,
            )?;
            residual = work.entropy(self, temperature, &mut cancelled)? - target_entropy_j_k;
            if residual.abs() <= controls.entropy_tolerance_j_k {
                return work.publish_sv(
                    self,
                    temperature,
                    volume,
                    target_entropy_j_k,
                    controls,
                    counts,
                    step,
                    volume_residual,
                    &mut cancelled,
                );
            }
            if residual < 0.0 {
                low = t;
            } else {
                high = t;
            }
            let response = work.response(
                self,
                temperature,
                controls.volume.equilibrium,
                &mut cancelled,
            )?;
            let next = t - (residual / response.heat_capacity_v_j_k) * t;
            // Preserve bracket contraction even when chemical response changes
            // rapidly; a small step cannot replace an entropy-residual check.
            let margin = 0.01 * (high - low);
            t = if next.is_finite() && next > low + margin && next < high - margin {
                next
            } else {
                low + 0.5 * (high - low)
            };
            if t == low || t == high {
                break;
            }
        }
        Err(EquilibriumExpansionErrorV1::EntropyNoConvergence {
            iterations: steps,
            residual_j_k: residual,
        }
        .into())
    }
}

impl Workspace {
    #[allow(clippy::too_many_arguments)]
    fn publish_sv(
        &mut self,
        owner: &IdealGasEquilibriumV1,
        temperature: Temperature,
        volume: Volume,
        target: f64,
        controls: EquilibriumVolumeEntropySolveV1,
        counts: Counts,
        temperature_iterations: usize,
        relative_volume_residual: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<IdealGasIsochoricEquilibriumV1, EquilibriumVolumeErrorV1> {
        let internal_energy = Energy::new(self.internal_energy(temperature)?);
        let thermodynamic = self.publish_response(
            owner,
            temperature,
            controls.volume.equilibrium,
            counts.composition,
            temperature_iterations,
            Some(target),
            cancelled,
        )?;
        poll(cancelled)?;
        Ok(IdealGasIsochoricEquilibriumV1 {
            thermodynamic,
            volume,
            internal_energy,
            target_internal_energy: None,
            relative_volume_residual,
            pressure_iterations: counts.pressure,
        })
    }
}
