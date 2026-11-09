//! Original constitutive stress, spring storage and physical heat in a harmonic cycle.

use super::plane_stress::projection::{Sum, apply, apply_spectrum, weighted_norm, weighted_norms};
use super::{Mandel6, TensorProny, ViscoError};
use crate::visco::prony_spectrum::PronySpectrumWeights;
use fs_math::{c64::C64, dd::Dd};

mod window;
pub use window::{TensorPronyHarmonicWindow, TensorPronyHarmonicWindowPoint};
mod spectrum;
pub use spectrum::{
    MAX_HARMONIC_MATERIAL_LINES, TensorPronyCausalHarmonicPoint, TensorPronyCausalSpectrumWindow,
    TensorPronyHarmonicSpectrum, TensorPronySpectrumWindow,
};

fn invalid() -> ViscoError {
    ViscoError::Parameters {
        what: "original harmonic material observation is not representable",
    }
}
fn finite(vector: &[C64; 6]) -> Result<(), ViscoError> {
    if vector
        .iter()
        .any(|v| !v.re.is_finite() || !v.im.is_finite())
    {
        Err(invalid())
    } else {
        Ok(())
    }
}
fn scalar(value: f64) -> Result<f64, ViscoError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid())
    }
}
fn positive(value: f64) -> Result<f64, ViscoError> {
    if value.is_finite() && value >= 0.0 {
        Ok(value)
    } else {
        Err(invalid())
    }
}

/// Physical cycle means and complete stress at one original material point.
/// Peak phasors use exp(-i omega t), with Mandel order
/// xx, yy, zz, sqrt(2)xy, sqrt(2)yz, sqrt(2)zx.
/// Energies are per original reference volume. The source's strain limit applies
/// to the entire physical cycle before any result is returned.
#[derive(Clone, Debug, PartialEq)]
pub struct TensorPronyHarmonicPoint {
    /// Complete storage/loss stress peak phasor [Pa], including all mutual terms.
    pub stress_mandel_pa: [C64; 6],
    /// Mean energy in the original equilibrium spring [J/m3].
    pub mean_equilibrium_energy_j_m3: f64,
    /// Mean energy in all original Maxwell springs [J/m3].
    pub mean_maxwell_energy_j_m3: f64,
    /// Mean total recoverable material energy [J/m3].
    pub mean_stored_energy_j_m3: f64,
    /// Physical Kelvin dashpot heat rate [W/m3].
    pub mean_kelvin_heat_w_m3: f64,
    /// Physical heat rate in all original Maxwell dashpots [W/m3].
    pub mean_maxwell_heat_w_m3: f64,
    /// Sum of the independently evaluated positive physical heat channels [W/m3].
    pub mean_material_heat_w_m3: f64,
    /// Independent Re(stress^H strain_rate)/2 [W/m3], without forcing closure.
    pub mean_stress_power_w_m3: f64,
    /// Signed stress work minus original physical heat [W/m3].
    pub power_balance_residual_w_m3: f64,
}

/// One original Maxwell branch at the same physical frequency and strain.
#[derive(Clone, Debug, PartialEq)]
pub struct TensorPronyHarmonicBranch {
    /// Peak elastic spring strain q satisfying tau*q_dot+q=tau*epsilon_dot.
    /// Rounded to f64 for observation; stress and heat keep their independent
    /// physical scale even when a tiny reported spring strain rounds to zero.
    pub spring_strain_mandel: [C64; 6],
    /// Original branch stress C_j*q [Pa].
    pub stress_mandel_pa: [C64; 6],
    /// Mean original spring energy q^H C_j q/4 [J/m3].
    pub mean_stored_energy_j_m3: f64,
    /// Mean original dashpot heat q^H C_j q/(2*tau) [W/m3].
    pub mean_heat_w_m3: f64,
}

/// Cached material spectrum, borrowing the immutable original tensor law.
/// Preparation is O(branches*6^2); point observations are allocation-free and
/// use existing Gram factors for storage/heat, independent of stress work.
/// This is the periodic particular response, without startup transients or an
/// inferred thermal temperature evolution. Source conditions remain fixed.
pub struct TensorPronyHarmonic<'a> {
    law: &'a TensorProny,
    omega: f64,
    stress: [[C64; 6]; 6],
    weights: Vec<PronySpectrumWeights>,
}

impl TensorProny {
    /// Prepare original harmonic stress and physical material-energy observations.
    /// A positive finite frequency inside the ORIGINAL source band is required.
    /// Static strain uses the original static/time-domain owner, without applying
    /// sinusoidal averaging factors to DC. The original law is never refitted.
    pub fn prepare_harmonic_material(
        &self,
        omega: f64,
    ) -> Result<TensorPronyHarmonic<'_>, ViscoError> {
        if !omega.is_finite() || omega <= 0.0 {
            return Err(ViscoError::Parameters {
                what: "harmonic material observation requires positive finite frequency",
            });
        }
        let (storage, loss) = self.complex_modulus(omega)?;
        let weights = self
            .branches
            .iter()
            .map(|b| PronySpectrumWeights::new(omega, b.tau_s).ok_or_else(invalid))
            .collect::<Result<Vec<_>, _>>()?;
        Ok(TensorPronyHarmonic {
            law: self,
            omega,
            weights,
            stress: core::array::from_fn(|i| {
                core::array::from_fn(|j| C64::new(storage[i][j], -loss[i][j]))
            }),
        })
    }
}

impl<'a> TensorPronyHarmonic<'a> {
    /// Immutable law which supplies every coefficient, factor and source limit.
    pub const fn law(&self) -> &'a TensorProny {
        self.law
    }
    /// Admitted positive angular frequency [rad/s].
    pub const fn omega_rad_s(&self) -> f64 {
        self.omega
    }
    /// Prepared numeric storage, excluding the separately owned source law and allocator metadata.
    pub fn numerical_payload_bytes(&self) -> usize {
        core::mem::size_of::<Self>()
            + self.weights.len() * core::mem::size_of::<PronySpectrumWeights>()
    }

    /// Observe the complete physical harmonic material response at one point.
    /// The original Gram factors, Maxwell spring strains and relaxation times
    /// determine storage and heat; neither is inferred from stress work or a
    /// modal loss ratio. Signed work imbalance exposes finite-precision closure.
    pub fn evaluate(&self, strain: &[C64; 6]) -> Result<TensorPronyHarmonicPoint, ViscoError> {
        self.law.validate_harmonic_mandel_strain(strain)?;
        let stress = apply(&self.stress, strain);
        finite(&stress)?;
        let equilibrium = self.energy(&self.law.equilibrium_gram, strain, 1.0, 0.25, &[])?;
        let kelvin = if let Some(factor) = &self.law.kv_factor {
            let gram = core::array::from_fn(|i| core::array::from_fn(|j| factor[j][i]));
            self.energy(&gram, strain, self.omega, 0.5, &[])?
        } else {
            0.0
        };
        let (mut storage, mut heat) = (Sum::default(), Sum::default());
        for branch in 0..self.weights.len() {
            let (e, d) = self.branch_energy(branch, strain)?;
            storage.add([e, 1.0, 1.0, 1.0]);
            heat.add([d, 1.0, 1.0, 1.0]);
        }
        let maxwell_energy = positive(storage.value())?;
        let maxwell_heat = positive(heat.value())?;
        let total_energy =
            positive((Dd::from_f64(equilibrium) + Dd::from_f64(maxwell_energy)).to_f64())?;
        let total_heat = positive((Dd::from_f64(kelvin) + Dd::from_f64(maxwell_heat)).to_f64())?;
        let mut work = Sum::default();
        for i in 0..6 {
            work.add([0.5, self.omega, stress[i].re, strain[i].im]);
            work.add([-0.5, self.omega, stress[i].im, strain[i].re]);
        }
        let work = scalar(work.value())?;
        let balance = scalar((Dd::from_f64(work) - Dd::from_f64(total_heat)).to_f64())?;
        Ok(TensorPronyHarmonicPoint {
            stress_mandel_pa: stress,
            mean_equilibrium_energy_j_m3: equilibrium,
            mean_maxwell_energy_j_m3: maxwell_energy,
            mean_stored_energy_j_m3: total_energy,
            mean_kelvin_heat_w_m3: kelvin,
            mean_maxwell_heat_w_m3: maxwell_heat,
            mean_material_heat_w_m3: total_heat,
            mean_stress_power_w_m3: work,
            power_balance_residual_w_m3: balance,
        })
    }

    /// Recover one branch in ORIGINAL source order, with no history mutation.
    pub fn branch_response(
        &self,
        branch: usize,
        strain: &[C64; 6],
    ) -> Result<TensorPronyHarmonicBranch, ViscoError> {
        self.law.validate_harmonic_mandel_strain(strain)?;
        if branch >= self.weights.len() {
            return Err(ViscoError::Parameters {
                what: "harmonic Maxwell branch index is outside the original law",
            });
        }
        let identity: Mandel6 =
            core::array::from_fn(|i| core::array::from_fn(|j| if i == j { 1.0 } else { 0.0 }));
        let spectrum = self.weights[branch].primal_binary();
        let q = apply_spectrum(&identity, strain, spectrum);
        finite(&q)?;
        let (energy, heat) = self.branch_energy(branch, strain)?;
        let stress = apply_spectrum(&self.law.branches[branch].c_pa, strain, spectrum);
        finite(&stress)?;
        Ok(TensorPronyHarmonicBranch {
            spring_strain_mandel: q,
            stress_mandel_pa: stress,
            mean_stored_energy_j_m3: energy,
            mean_heat_w_m3: heat,
        })
    }

    fn branch_energy(&self, index: usize, strain: &[C64; 6]) -> Result<(f64, f64), ViscoError> {
        let branch = &self.law.branches[index];
        let gram = branch.resolved_gram_rows(&self.law.factors[index]);
        // |q|^2 = storage_weight * |epsilon|^2, even when q itself
        // underflows. Retain the weight, Gram and tau scales until output.
        let spectrum = self.weights[index].primal_binary()[0];
        let [energy, heat] = weighted_norms(
            &gram,
            strain,
            1.0,
            spectrum,
            [(0.25, &[]), (0.5, &[branch.tau_s])],
        );
        Ok((positive(energy)?, positive(heat)?))
    }

    fn energy(
        &self,
        gram: &Mandel6,
        strain: &[C64; 6],
        rate: f64,
        weight: f64,
        denominator: &[f64],
    ) -> Result<f64, ViscoError> {
        positive(weighted_norm(
            gram,
            strain,
            rate,
            (1.0, 0),
            weight,
            denominator,
        ))
    }
}
