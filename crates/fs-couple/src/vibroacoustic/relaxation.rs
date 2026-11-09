//! Hereditary structural stiffness in the original cavity/radiation coordinates.
//!
//! This is the exact real-frequency image of the existing fs-phs Maxwell
//! energy H_j=(c_j.q-z_j)^2/2 and resistance 1/tau_j. The shared material
//! spectrum owner supplies every frequency coefficient; no second Maxwell
//! law, fitted modal Q, diagonal approximation or internal-state truncation
//! is introduced.

use super::{C64, VibroError};
use crate::bem_loaded_modes::{BemModesError, memory::RetainedMemory};
use crate::plate_viscosity::relaxation::{ModalPlateRelaxation, PlateBendingRelaxation};
use fs_exec::CancelGate;
use fs_material::visco::PronySpectrumWeights;
use fs_math::dd::Dd;
use fs_phs::RelaxationBranch;

/// Complete positive material strain channels in the ORIGINAL mass-normalized
/// structural displacement coordinates.
///
/// The structural frequencies supplied to `VibroacousticModel` must describe
/// RELAXED equilibrium stiffness. Each branch contributes its excess stiffness
/// only once. Channel signs, equal-time branches and exact nullspaces survive.
/// Source identity, actual basis matching and constitutive applicability belong
/// to the original material/geometry owner; matching dimensions is insufficient.
#[derive(Clone, Debug)]
pub struct ModalStructuralRelaxation {
    modes: usize,
    memory: RetainedMemory,
    band: (f64, f64),
    plate_source: Option<PlateBendingRelaxation>,
}

fn invalid(what: &'static str) -> VibroError {
    VibroError::BadParameter { what }
}

impl ModalStructuralRelaxation {
    /// Admit the SAME declared arms accepted by the native BEM/pHS owner.
    ///
    /// Projections have `modes` displacement entries, not interleaved q/p
    /// entries. The caller supplies the material's angular-frequency support
    /// [rad/s] and a branch ceiling no larger than the shared limit of 256.
    /// All branches retain their original order. No loss is inferred from a
    /// mode frequency or an authored damping ratio.
    ///
    /// # Errors
    /// Invalid source band, shape or material parameters, a branch-cap refusal,
    /// unrepresentable physical channels, or cancellation.
    pub fn try_new(
        modes: usize,
        branches: &[RelaxationBranch],
        angular_frequency_band_rad_s: (f64, f64),
        maximum_branches: usize,
        gate: &CancelGate,
    ) -> Result<Self, VibroError> {
        if gate.is_requested() {
            return Err(VibroError::Cancelled);
        }
        let (low, high) = angular_frequency_band_rad_s;
        if modes == 0 || !low.is_finite() || !high.is_finite() || low < 0.0 || high < low {
            return Err(invalid(
                "structural relaxation requires nonempty modes and a finite ordered source band",
            ));
        }
        let memory = RetainedMemory::from_branches(modes, branches, maximum_branches, &mut || {
            gate.is_requested()
        })
        .map_err(|error| match error {
            BemModesError::Cancelled => VibroError::Cancelled,
            BemModesError::Invalid { what } => invalid(what),
            _ => invalid("original structural relaxation admission refused"),
        })?;
        Ok(Self {
            modes,
            memory,
            band: angular_frequency_band_rad_s,
            plate_source: None,
        })
    }

    /// Reuse the original source-resolved plate projection and every regional
    /// material band. No elastic solve, factorization or channel merging occurs.
    /// The same relaxed plate basis must have supplied the structural modes.
    ///
    /// # Errors
    /// The original source bands have no common support, the complete projected
    /// channels exceed the original cap, or cancellation is observed.
    pub fn from_plate(image: &ModalPlateRelaxation, gate: &CancelGate) -> Result<Self, VibroError> {
        if gate.is_requested() {
            return Err(VibroError::Cancelled);
        }
        let source = image.source();
        let band = source
            .omega_bands_rad_s()
            .iter()
            .fold((0.0_f64, f64::MAX), |(low, high), &(a, b)| {
                (low.max(a), high.min(b))
            });
        let branches: Vec<_> = image
            .coupling_rows_per_s()
            .iter()
            .zip(image.relaxation_times_s())
            .map(|(row, &time)| RelaxationBranch {
                projection: row.clone(),
                stiffness: 1.0,
                relaxation_time_s: time,
            })
            .collect();
        let mut result = Self::try_new(
            image.mode_count(),
            &branches,
            band,
            image.limits().max_memory_states,
            gate,
        )?;
        result.plate_source = Some(source.clone());
        Ok(result)
    }

    #[must_use]
    pub const fn mode_count(&self) -> usize {
        self.modes
    }

    /// Original full signed strain factors c_j [1/s].
    #[must_use]
    pub fn coupling_rows_per_s(&self) -> &[Vec<f64>] {
        &self.memory.rows
    }

    #[must_use]
    pub fn relaxation_times_s(&self) -> &[f64] {
        &self.memory.times
    }

    #[must_use]
    pub const fn angular_frequency_band_rad_s(&self) -> (f64, f64) {
        self.band
    }

    /// Original material/property selection, when constructed from that owner.
    #[must_use]
    pub fn plate_source(&self) -> Option<&PlateBendingRelaxation> {
        self.plate_source.as_ref()
    }

    pub(super) fn at_frequency(&self, omega: f64) -> Result<FrequencyRelaxation<'_>, VibroError> {
        if !omega.is_finite() || omega < self.band.0 || omega > self.band.1 {
            return Err(invalid(
                "structural relaxation frequency lies outside its original material support",
            ));
        }
        let weights = self
            .memory
            .times
            .iter()
            .map(|&tau| {
                PronySpectrumWeights::new(omega, tau).ok_or_else(|| {
                    invalid("original structural relaxation spectrum is unrepresentable")
                })
            })
            .collect::<Result<Vec<_>, _>>()?;
        Ok(FrequencyRelaxation {
            rows: &self.memory.rows,
            weights,
            omega,
        })
    }
}

/// Prepared once per FRF, then shared by the matrix and independent heat/storage
/// contractions. No extra internal coordinates enlarge the frequency solve.
pub(super) struct FrequencyRelaxation<'a> {
    rows: &'a [Vec<f64>],
    weights: Vec<PronySpectrumWeights>,
    omega: f64,
}

impl FrequencyRelaxation<'_> {
    pub(super) fn add_stiffness(
        &self,
        matrix: &mut [C64],
        dimension: usize,
        retained: usize,
    ) -> Result<(), VibroError> {
        if self.rows.is_empty() {
            return Ok(());
        }
        for r in 0..retained {
            for s in r..retained {
                let mut storage = Dd::ZERO;
                let mut loss = Dd::ZERO;
                for (row, weights) in self.rows.iter().zip(&self.weights) {
                    let values = weights.primal(row[r] * row[s]).ok_or_else(|| {
                        invalid("projected structural relaxation stiffness is unrepresentable")
                    })?;
                    storage = storage + Dd::from_f64(values[0]);
                    loss = loss + Dd::from_f64(values[1]);
                }
                // The material owner returns positive E'' for exp(+i omega t).
                // This solver uses exp(-i omega t), hence storage - i loss.
                let value = C64::new(storage.to_f64(), -loss.to_f64());
                matrix[r * dimension + s] = matrix[r * dimension + s] + value;
                if r != s {
                    matrix[s * dimension + r] = matrix[s * dimension + r] + value;
                }
                if !matrix[r * dimension + s].re.is_finite()
                    || !matrix[r * dimension + s].im.is_finite()
                    || !matrix[s * dimension + r].re.is_finite()
                    || !matrix[s * dimension + r].im.is_finite()
                {
                    return Err(invalid(
                        "loaded structural relaxation equation is unrepresentable",
                    ));
                }
            }
        }
        Ok(())
    }

    /// Mean spring storage [J] and nonnegative heat rate [W], from the SAME
    /// original strain channels. Prefix truncation sets omitted displacements
    /// to zero without changing any retained channel or relaxation time.
    pub(super) fn energy_and_heat(&self, displacement: &[C64]) -> Result<(f64, f64), VibroError> {
        let mut energy = Dd::ZERO;
        let mut heat = Dd::ZERO;
        for (row, weights) in self.rows.iter().zip(&self.weights) {
            let mut re = Dd::ZERO;
            let mut im = Dd::ZERO;
            for (&coefficient, &value) in row.iter().zip(displacement) {
                re = re + Dd::from_f64(coefficient) * Dd::from_f64(value.re);
                im = im + Dd::from_f64(coefficient) * Dd::from_f64(value.im);
            }
            let square = (re * re + im * im).to_f64();
            let [stored, loss] = weights
                .primal(square)
                .ok_or_else(|| invalid("structural relaxation channel work is unrepresentable"))?;
            energy = energy + Dd::from_f64(0.25) * Dd::from_f64(stored);
            heat = heat + Dd::from_f64(0.5 * self.omega) * Dd::from_f64(loss);
        }
        let values = (energy.to_f64(), heat.to_f64());
        if !values.0.is_finite() || !values.1.is_finite() || values.0 < 0.0 || values.1 < 0.0 {
            return Err(invalid(
                "structural relaxation energy or heat is unrepresentable",
            ));
        }
        Ok(values)
    }
}

#[cfg(test)]
mod tests;
