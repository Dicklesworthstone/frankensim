//! Actual interval work, heat and recoverable energy in a periodic material field.

use super::super::plane_stress::projection::{Sum, window::HarmonicEnergyWindow};
use super::{TensorPronyHarmonic, ViscoError, invalid, positive, scalar};
use fs_fft::harmonic_window::{HarmonicWindow, HarmonicWindowProfile};
use fs_math::{c64::C64, dd::Dd};

/// Physical, unweighted interval account per original reference volume.
/// The original settled sinusoidal particular response supplies every channel;
/// preparing or observing a window neither starts nor advances material history.
#[derive(Clone, Debug, PartialEq)]
pub struct TensorPronyHarmonicWindowPoint {
    /// Recoverable equilibrium and Maxwell spring energy at the start [J/m3].
    pub start_stored_energy_j_m3: f64,
    /// Recoverable equilibrium and Maxwell spring energy at the end [J/m3].
    pub end_stored_energy_j_m3: f64,
    /// Independently integrated original storage rate [J/m3]. This preserves
    /// short-interval changes which rounded endpoint subtraction can erase.
    pub stored_energy_change_j_m3: f64,
    /// Physical heat in the original Kelvin dashpot [J/m3].
    pub kelvin_heat_j_m3: f64,
    /// Physical heat in the original Maxwell dashpots [J/m3].
    pub maxwell_heat_j_m3: f64,
    /// Sum of independently calculated physical heat channels [J/m3].
    pub material_heat_j_m3: f64,
    /// Independently integrated complete stress : strain_rate [J/m3].
    pub stress_work_j_m3: f64,
    /// Signed work minus stored-energy change minus physical heat [J/m3].
    pub energy_balance_residual_j_m3: f64,
}

/// Allocation-free temporal cache borrowing the same original material spectrum.
/// Time trigonometry is prepared once and reused across points at this frequency.
pub struct TensorPronyHarmonicWindow<'p, 'a> {
    source: &'p TensorPronyHarmonic<'a>,
    window: HarmonicWindow,
    temporal: HarmonicEnergyWindow,
}

impl<'a> TensorPronyHarmonic<'a> {
    /// Prepare an actual physical interval for repeated material observations.
    /// Only a uniform window represents unweighted work, heat and endpoint
    /// storage conservation; shaped detector sensitivities explicitly refuse.
    pub fn prepare_window(
        &self,
        window: HarmonicWindow,
    ) -> Result<TensorPronyHarmonicWindow<'_, 'a>, ViscoError> {
        if window.profile() != HarmonicWindowProfile::Uniform {
            return Err(ViscoError::Parameters {
                what: "physical material energy requires a uniform observation window",
            });
        }
        let temporal = HarmonicEnergyWindow::new(self.omega, window).map_err(|_| invalid())?;
        Ok(TensorPronyHarmonicWindow {
            source: self,
            window,
            temporal,
        })
    }

    /// Convenience observation; use `prepare_window` to reuse temporal weights.
    pub fn observe_window(
        &self,
        strain: &[C64; 6],
        window: HarmonicWindow,
    ) -> Result<TensorPronyHarmonicWindowPoint, ViscoError> {
        self.prepare_window(window)?.evaluate(strain)
    }
}

impl<'p, 'a> TensorPronyHarmonicWindow<'p, 'a> {
    /// The unchanged original material spectrum.
    pub const fn source(&self) -> &'p TensorPronyHarmonic<'a> {
        self.source
    }
    /// Exact admitted interval, retaining start and duration separately.
    pub const fn window(&self) -> HarmonicWindow {
        self.window
    }

    /// Observe original spring storage, physical dashpot heat and stress work.
    /// The complete strain cycle still satisfies the original source limit,
    /// even when this interval samples a small displacement near a phase node.
    /// Original binary Gram/spectrum/rate/tau scales survive until each physical
    /// result, including a representable heat with unrepresentable branch strain.
    pub fn evaluate(
        &self,
        strain: &[C64; 6],
    ) -> Result<TensorPronyHarmonicWindowPoint, ViscoError> {
        let law = self.source.law;
        law.validate_harmonic_mandel_strain(strain)?;
        let equilibrium = self
            .temporal
            .spring(&law.equilibrium_gram, strain, [(1.0, 0), (0.0, 0)], None)
            .map_err(|_| invalid())?;
        let mut storage = [equilibrium[0], equilibrium[1], equilibrium[2]].map(Dd::from_f64);
        let mut maxwell = Dd::ZERO;
        let mut work = Sum::default();
        self.temporal
            .add_stress_work(&law.c_eq, strain, [(1.0, 0), (0.0, 0)], &mut work);
        for (index, branch) in law.branches.iter().enumerate() {
            let gram = branch.resolved_gram_rows(&law.factors[index]);
            let response = self
                .temporal
                .spring(
                    &gram,
                    strain,
                    self.source.weights[index].primal_binary(),
                    Some(branch.tau_s),
                )
                .map_err(|_| invalid())?;
            for j in 0..3 {
                storage[j] = storage[j] + Dd::from_f64(response[j]);
            }
            maxwell = maxwell + Dd::from_f64(response[3]);
            self.temporal.add_stress_work(
                &branch.c_pa,
                strain,
                self.source.weights[index].primal_binary(),
                &mut work,
            );
        }
        let kelvin = if let Some(factor) = &law.kv_factor {
            let gram = core::array::from_fn(|i| core::array::from_fn(|j| factor[j][i]));
            self.temporal
                .kelvin_heat(&gram, strain)
                .map_err(|_| invalid())?
        } else {
            0.0
        };
        if let Some(viscosity) = &law.kv_pa_s {
            self.temporal.add_stress_work(
                viscosity,
                strain,
                [(0.0, 0), (self.source.omega, 0)],
                &mut work,
            );
        }
        let work = scalar(work.value())?;
        let heat = maxwell + Dd::from_f64(kelvin);
        let balance = Dd::from_f64(work) - storage[2] - heat;
        Ok(TensorPronyHarmonicWindowPoint {
            start_stored_energy_j_m3: positive(storage[0].to_f64())?,
            end_stored_energy_j_m3: positive(storage[1].to_f64())?,
            stored_energy_change_j_m3: scalar(storage[2].to_f64())?,
            kelvin_heat_j_m3: kelvin,
            maxwell_heat_j_m3: positive(maxwell.to_f64())?,
            material_heat_j_m3: positive(heat.to_f64())?,
            stress_work_j_m3: work,
            energy_balance_residual_j_m3: scalar(balance.to_f64())?,
        })
    }
}
