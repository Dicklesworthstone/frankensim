//! Finite reversible exchange of one original small-strain thermal potential.
//!
//! The same admitted source law resolves both endpoints. Mechanical work uses
//! endpoint-average stress and the actual total Mandel strain increment; this
//! is a trapezoidal path approximation, not an exact arbitrary thermal path.
//! Polarization retains changes of Helmholtz energy, entropy and internal
//! energy without subtracting their rounded absolute endpoint values.
//! Reversible heat is work minus internal-energy change. No residual,
//! numerical loss, caloric reference energy or dissipation is added to it.

use super::{
    AffineThermoelastic, IsotropicThermalModuli, ResolvedThermoelastic,
    ResolvedThermoelasticLaw, ThermoelasticError, ThermoelasticState,
    isotropic_thermal_moduli,
};
use super::tensor::{ThermalCoefficientIncrement, coefficient_increments};
use fs_math::dd::Dd;

/// Extensive endpoint account on the original reference volume.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermoelasticExchange {
    pub before: ThermoelasticState,
    pub after: ThermoelasticState,
    /// V * (sigma0 + sigma1) : (epsilon1 - epsilon0) / 2 [J].
    pub mechanical_work_j: f64,
    /// Signed free-strain contribution -V * mean(sigma) : Delta free [J].
    pub expansion_work_j: f64,
    /// Signed coefficient contribution V * e0 : Delta C : e1 / 2 [J].
    pub stiffness_work_j: f64,
    pub free_energy_change_j: f64,
    pub entropy_change_j_k: f64,
    pub internal_energy_change_j: f64,
    /// Positive into the separate caloric store: W - Delta U [J].
    pub reversible_heat_j: f64,
    /// Rounded endpoint Delta Psi minus the retained increment [J].
    /// This representational diagnostic is never deposited as heat.
    pub free_energy_balance_residual_j: f64,
    retained: [Dd; 7],
}

impl ThermoelasticExchange {
    pub const fn retained_mechanical_work_j(&self) -> Dd {
        self.retained[0]
    }

    pub const fn retained_free_energy_change_j(&self) -> Dd {
        self.retained[1]
    }

    pub const fn retained_entropy_change_j_k(&self) -> Dd {
        self.retained[2]
    }

    pub const fn retained_internal_energy_change_j(&self) -> Dd {
        self.retained[3]
    }

    pub const fn retained_reversible_heat_j(&self) -> Dd {
        self.retained[4]
    }

    pub const fn retained_expansion_work_j(&self) -> Dd {
        self.retained[5]
    }

    pub const fn retained_stiffness_work_j(&self) -> Dd {
        self.retained[6]
    }
}

fn dd(value: f64) -> Dd {
    Dd::from_f64(value)
}

fn finite(value: Dd) -> Result<f64, ThermoelasticError> {
    super::finite(value.to_f64())
}

fn dot(a: &[Dd; 6], b: &[Dd; 6]) -> Dd {
    (0..6).fold(Dd::ZERO, |sum, i| sum + a[i] * b[i])
}

// Original source-declared modes must survive even when the rounded dense
// normal block cannot represent the weak bulk or shear coefficient.
enum Operator {
    Isotropic { bulk: Dd, shear: Dd },
    Tensor([[Dd; 6]; 6]),
}

impl Operator {
    fn bilinear(&self, a: &[Dd; 6], b: &[Dd; 6]) -> Dd {
        match self {
            Self::Isotropic { bulk, shear } => {
                let trace_a = a[0] + a[1] + a[2];
                let trace_b = b[0] + b[1] + b[2];
                // Pair differences preserve an exactly hydrostatic mode.
                let mut dev = ((a[0] - a[1]) * (b[0] - b[1])
                    + (a[0] - a[2]) * (b[0] - b[2])
                    + (a[1] - a[2]) * (b[1] - b[2])) / dd(3.0);
                for i in 3..6 {
                    dev = dev + a[i] * b[i];
                }
                *bulk * trace_a * trace_b + dd(2.0) * *shear * dev
            }
            Self::Tensor(c) => (0..6).fold(Dd::ZERO, |sum, i| {
                sum + a[i] * (0..6).fold(Dd::ZERO, |row, j| row + c[i][j] * b[j])
            }),
        }
    }
}

fn operators(
    source: &ResolvedThermoelastic,
) -> Result<[Operator; 2], ThermoelasticError> {
    let moduli = match &source.response {
        ResolvedThermoelasticLaw::Isotropic(response) => Some(isotropic_thermal_moduli(
            response.material.youngs,
            response.material.poisson,
            response.youngs_slope_pa_k,
            response.poisson_slope_per_k,
        )?),
        ResolvedThermoelasticLaw::Tensor(response) => response.isotropic_moduli,
    };
    if let Some(IsotropicThermalModuli {
        bulk, shear, bulk_first, shear_first, ..
    }) = moduli
    {
        return Ok([
            Operator::Isotropic { bulk: dd(bulk), shear: dd(shear) },
            Operator::Isotropic { bulk: dd(bulk_first), shear: dd(shear_first) },
        ]);
    }
    let c = source.stiffness_mandel_pa();
    let ct = source.temperature_derivatives()?.stiffness_mandel_pa_k;
    Ok([c, ct].map(|matrix| Operator::Tensor(
        core::array::from_fn(|i| core::array::from_fn(|j| dd(matrix[i][j]))),
    )))
}

impl AffineThermoelastic {
    /// Integrate one declared endpoint-work interval on this original source.
    ///
    /// Both temperatures, stiffness/expansion axes, reference strain and source
    /// claims belong to this SAME law. It admits both total and elastic strain
    /// component boxes and uses the original fixed reference volume. The
    /// returned states are the ordinary native potential evaluations, suitable
    /// for comparison with actual mechanical stress recovery.
    ///
    /// The admitted instantaneous expansion is affine, so its finite integral
    /// is mean(alpha) * Delta T, retained before subtracting total strain.
    /// Elastic coefficient increments use the original affine engineering
    /// paths and exact rational/compliance-inverse identities. This avoids
    /// subtracting rounded free-strain, stiffness or derivative endpoints.
    /// The finite work is a trapezoidal approximation and requires refinement
    /// on changing-temperature/load paths. This does not assert entropy
    /// conservation, solve mechanics or advance a caloric history.
    pub fn exchange(
        &self,
        before_temperature_k: f64,
        before_total_strain_mandel: [f64; 6],
        after_temperature_k: f64,
        after_total_strain_mandel: [f64; 6],
        reference_volume_m3: f64,
    ) -> Result<ThermoelasticExchange, ThermoelasticError> {
        let source0 = self.at_temperature(before_temperature_k)?;
        let changed = if after_temperature_k == before_temperature_k {
            None
        } else {
            Some(self.at_temperature(after_temperature_k)?)
        };
        let source1 = changed.as_ref().unwrap_or(&source0);
        let before = source0.state(before_total_strain_mandel, reference_volume_m3)?;
        let after = source1.state(after_total_strain_mandel, reference_volume_m3)?;
        let d0 = source0.temperature_derivatives()?;
        let d1 = source1.temperature_derivatives()?;
        let [c0, ct0] = operators(&source0)?;
        let [c1, ct1] = operators(source1)?;
        let (dc, dct) = match coefficient_increments(self, &source0, source1)? {
            ThermalCoefficientIncrement::Isotropic { moduli, first } => (
                Operator::Isotropic { bulk: moduli[0], shear: moduli[1] },
                Operator::Isotropic { bulk: first[0], shear: first[1] },
            ),
            ThermalCoefficientIncrement::Tensor { stiffness, first } => {
                (Operator::Tensor(stiffness), Operator::Tensor(first))
            }
        };
        let half = dd(0.5);
        let dt = dd(after_temperature_k) - dd(before_temperature_k);
        let total_increment: [Dd; 6] = core::array::from_fn(|i| {
            dd(after_total_strain_mandel[i]) - dd(before_total_strain_mandel[i])
        });
        let alpha: [Dd; 6] = core::array::from_fn(|i| {
            half * (dd(d0.free_strain_mandel_per_k[i])
                + dd(d1.free_strain_mandel_per_k[i]))
        });
        let delta_alpha: [Dd; 6] = core::array::from_fn(|i| {
            half * (dd(d0.free_strain_mandel_per_k2[i])
                + dd(d1.free_strain_mandel_per_k2[i])) * dt
        });
        let delta_free: [Dd; 6] = core::array::from_fn(|i| alpha[i] * dt);
        // Do not recover this from rounded (total - free) endpoints: those
        // may be identical even when the original input increment is nonzero.
        let de: [Dd; 6] = core::array::from_fn(|i| {
            total_increment[i] - delta_free[i]
        });
        let e0 = before.elastic_strain_mandel.map(dd);
        let e1 = after.elastic_strain_mandel.map(dd);
        let mean_e = core::array::from_fn(|i| half * (e0[i] + e1[i]));
        let mean_stress: [Dd; 6] = core::array::from_fn(|i| {
            half * (dd(before.stress_mandel_pa[i]) + dd(after.stress_mandel_pa[i]))
        });
        let volume = dd(reference_volume_m3);
        let work = volume * dot(&mean_stress, &total_increment);
        let expansion_work = -volume * dot(&mean_stress, &delta_free);
        let stiffness_work = volume * half * dc.bilinear(&e0, &e1);
        let free_change = volume * dot(&mean_stress, &de) + stiffness_work;
        // Delta(alpha:sigma) - Delta(e:C_T:e)/2, polarized before summation.
        let entropy_change = volume * (
            dot(&delta_alpha, &mean_stress)
                + half * (c0.bilinear(&alpha, &de) + c1.bilinear(&alpha, &de))
                + dc.bilinear(&alpha, &mean_e)
                - half * (
                    ct0.bilinear(&e0, &de) + ct1.bilinear(&e1, &de)
                        + dct.bilinear(&e0, &e1)
                )
        );
        let temperature = half * (dd(before_temperature_k) + dd(after_temperature_k));
        let entropy = half * (dd(before.entropy_j_k) + dd(after.entropy_j_k));
        let internal_change = free_change + temperature * entropy_change + entropy * dt;
        Ok(ThermoelasticExchange {
            before,
            after,
            mechanical_work_j: finite(work)?,
            expansion_work_j: finite(expansion_work)?,
            stiffness_work_j: finite(stiffness_work)?,
            free_energy_change_j: finite(free_change)?,
            entropy_change_j_k: finite(entropy_change)?,
            internal_energy_change_j: finite(internal_change)?,
            reversible_heat_j: finite(work - internal_change)?,
            free_energy_balance_residual_j: finite(
                dd(after.free_energy_j) - dd(before.free_energy_j) - free_change,
            )?,
            retained: [
                work, free_change, entropy_change, internal_change, work - internal_change,
                expansion_work, stiffness_work,
            ],
        })
    }
}
