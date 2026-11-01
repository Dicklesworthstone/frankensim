//! Phase selection from the original two tangent-plane searches and PT flash.
//!
//! A homogeneous result means that both declared searches converged without a
//! negative witness. It is not a global minimum certificate. A negative witness
//! must resolve an admitted two-phase flash; a failed split never falls back to
//! the homogeneous feed. Equal-Gibbs distinct feed roots refuse because PT alone
//! cannot determine their amount (including pure-fluid/azeotropic coexistence).

use super::flash::{
    MAX_PENG_ROBINSON_FLASH_ITERATIONS_V1, MAX_RACHFORD_RICE_ITERATIONS_V1,
    PengRobinsonFlashErrorV1 as E, PengRobinsonFlashOptionsV1, PengRobinsonFlashV1,
};
use super::stability::{PengRobinsonStabilityOptionsV1, PengRobinsonStabilityV1};
use super::{PengRobinsonMixtureV1, PengRobinsonPhaseRootV1};
use fs_qty::{Pressure, Temperature};

/// Complete bounded work and accuracy for one screened PT state.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct PengRobinsonPhaseScreenOptionsV1 {
    /// Both original Wilson-seeded tangent-plane searches must converge.
    pub stability: PengRobinsonStabilityOptionsV1,
    /// Original two-phase solver, used only when instability is witnessed.
    pub flash: PengRobinsonFlashOptionsV1,
}
impl PengRobinsonPhaseScreenOptionsV1 {
    /// Validate both original numerical policies before any EOS evaluation.
    pub fn validate(self) -> Result<(), E> {
        let s = self.stability;
        let f = self.flash;
        if s.maximum_iterations == 0
            || s.maximum_iterations > 1024
            || !s.stationarity_tolerance.is_finite()
            || s.stationarity_tolerance <= 0.0
            || !s.tangent_plane_tolerance.is_finite()
            || s.tangent_plane_tolerance <= 0.0
            || f.maximum_iterations == 0
            || f.maximum_iterations > MAX_PENG_ROBINSON_FLASH_ITERATIONS_V1
            || f.rachford_rice_iterations == 0
            || f.rachford_rice_iterations > MAX_RACHFORD_RICE_ITERATIONS_V1
            || !f.log_k_tolerance.is_finite()
            || f.log_k_tolerance <= 0.0
        {
            return Err(E::InvalidOptions);
        }
        Ok(())
    }

    /// Worst-case original cubic EOS calls for a screen followed by a split.
    pub fn maximum_eos_evaluation_bound(self) -> Result<usize, E> {
        self.validate()?;
        Ok(1 + 2 * self.stability.maximum_iterations + 6 * self.flash.maximum_iterations + 2)
    }
}

/// Original PT admission and its finite, explicitly non-global screening evidence.
#[derive(Clone, Debug, PartialEq)]
pub struct PengRobinsonPhaseScreenV1 {
    temperature: Temperature,
    pressure: Pressure,
    options: PengRobinsonPhaseScreenOptionsV1,
    stability: PengRobinsonStabilityV1,
    split: Option<PengRobinsonFlashV1>,
    eos_evaluation_upper_bound: usize,
}
impl PengRobinsonPhaseScreenV1 {
    /// Common absolute temperature.
    #[must_use]
    pub const fn temperature(&self) -> Temperature {
        self.temperature
    }
    /// Common pressure.
    #[must_use]
    pub const fn pressure(&self) -> Pressure {
        self.pressure
    }
    /// Exact finite screening, split tolerances and iteration caps used here.
    #[must_use]
    pub const fn options(&self) -> PengRobinsonPhaseScreenOptionsV1 {
        self.options
    }
    /// Both converged trials, unchanged feed roots and work evidence.
    #[must_use]
    pub const fn stability(&self) -> &PengRobinsonStabilityV1 {
        &self.stability
    }
    /// Admitted separated phases, when the screen exhibited feed instability.
    #[must_use]
    pub const fn two_phase(&self) -> Option<&PengRobinsonFlashV1> {
        self.split.as_ref()
    }
    /// Screened homogeneous root; no phase amount or liquid/vapor label is inferred.
    #[must_use]
    pub fn homogeneous(&self) -> Option<&PengRobinsonPhaseRootV1> {
        self.split.is_none().then(|| self.stability.feed())
    }
    /// Actual screen EOS work plus a conservative bound for the optional flash.
    /// This is not an observed total EOS count.
    #[must_use]
    pub const fn eos_evaluation_upper_bound(&self) -> usize {
        self.eos_evaluation_upper_bound
    }
}

impl PengRobinsonMixtureV1 {
    /// Select a screened homogeneous state or an original admitted phase split.
    /// Every original iteration polls cancellation; publication polls once more.
    pub fn phase_screen_pt_cancellable(
        &self,
        temperature: Temperature,
        pressure: Pressure,
        options: PengRobinsonPhaseScreenOptionsV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PengRobinsonPhaseScreenV1, E> {
        options.validate()?;
        let stability = self.stability_pt_cancellable(
            temperature,
            pressure,
            options.stability,
            &mut cancelled,
        )?;
        let mut eos_evaluation_upper_bound = stability.eos_evaluations();
        let split = if stability.instability_witness().is_some() {
            let split = self.two_phase_flash_pt_cancellable(
                temperature,
                pressure,
                options.flash,
                &mut cancelled,
            )?;
            for root in [split.liquid(), split.vapor()] {
                let pv = root.pressure_derivatives().1;
                if !pv.is_finite() || pv >= 0.0 {
                    return Err(E::Unrepresentable("screened split mechanical stability"));
                }
            }
            eos_evaluation_upper_bound += 6 * options.flash.maximum_iterations + 2;
            Some(split)
        } else {
            let feed = stability.feed();
            let pv = feed.pressure_derivatives().1;
            if !pv.is_finite() || pv >= 0.0 {
                return Err(E::Unrepresentable(
                    "screened homogeneous mechanical stability",
                ));
            }
            for root in stability.feed_roots() {
                let volume_gap =
                    (root.molar_volume_m3_per_mol() - feed.molar_volume_m3_per_mol()).abs();
                if root.pressure_derivatives().1 < 0.0
                    && volume_gap > 64.0 * f64::EPSILON * feed.molar_volume_m3_per_mol()
                    && (root.residual_gibbs_over_rt() - feed.residual_gibbs_over_rt()).abs()
                        <= options.stability.tangent_plane_tolerance
                {
                    return Err(E::Unrepresentable(
                        "equal-Gibbs phase amount is undetermined",
                    ));
                }
            }
            None
        };
        if cancelled() {
            return Err(E::Cancelled);
        }
        Ok(PengRobinsonPhaseScreenV1 {
            temperature,
            pressure,
            options,
            stability,
            split,
            eos_evaluation_upper_bound,
        })
    }
}
