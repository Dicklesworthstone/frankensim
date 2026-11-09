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
#[derive(Debug)]
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
    let ct = source.source_side_temperature_derivatives()?.stiffness_mandel_pa_k;
    Ok([c, ct].map(|matrix| Operator::Tensor(
        core::array::from_fn(|i| core::array::from_fn(|j| dd(matrix[i][j]))),
    )))
}

/// One original thermal interval, resolved once and reusable at bulk points.
///
/// Source/card resolution, exact coefficient increments and thermal derivatives
/// are prepared once. Evaluation retains the ordinary source strain admission
/// and potential; it performs no material lookup and allocates no point storage.
/// Reference volume is supplied per evaluation so retained quadrature weights
/// can be streamed without constructing an array of point responses.
#[derive(Debug)]
pub struct PreparedThermoelasticExchange {
    before: ResolvedThermoelastic,
    changed: Option<ResolvedThermoelastic>,
    c0: Operator,
    ct0: Operator,
    c1: Operator,
    ct1: Operator,
    dc: Operator,
    dct: Operator,
    alpha: [Dd; 6],
    delta_alpha: [Dd; 6],
    delta_free: [Dd; 6],
    delta_temperature: Dd,
    mean_temperature: Dd,
}

impl AffineThermoelastic {
    /// Resolve one source-owned interval before integrating many strain points.
    ///
    /// Both endpoints retain this same original law and its admitted range.
    /// This preparation does not change the finite path approximation used by
    /// exchange: endpoint-average mechanical stress still requires temporal
    /// refinement on changing-temperature/load paths.
    pub fn prepare_exchange(
        &self,
        before_temperature_k: f64,
        after_temperature_k: f64,
    ) -> Result<PreparedThermoelasticExchange, ThermoelasticError> {
        let source0 = self.at_temperature(before_temperature_k)?;
        let changed = if after_temperature_k == before_temperature_k {
            None
        } else {
            Some(self.at_temperature(after_temperature_k)?)
        };
        let source1 = changed.as_ref().unwrap_or(&source0);
        let d0 = source0.source_side_temperature_derivatives()?;
        let d1 = source1.source_side_temperature_derivatives()?;
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
        let alpha: [Dd; 6] = core::array::from_fn(|i| {
            half * (dd(d0.free_strain_mandel_per_k[i])
                + dd(d1.free_strain_mandel_per_k[i]))
        });
        // CTE is continuous at its knots, but one endpoint-average alpha
        // generally does not integrate a path crossing different slopes.
        // Retain each original principal segment before the fixed rotation.
        let [delta_free, delta_alpha] = self.expansion_increments(&source0, source1)?;
        Ok(PreparedThermoelasticExchange {
            before: source0,
            changed,
            c0,
            ct0,
            c1,
            ct1,
            dc,
            dct,
            alpha,
            delta_alpha,
            delta_free,
            delta_temperature: dt,
            mean_temperature: half * (dd(before_temperature_k) + dd(after_temperature_k)),
        })
    }

    /// Integrate one declared endpoint-work interval on this original source.
    ///
    /// Both temperatures, stiffness/expansion axes, reference strain and source
    /// claims belong to this SAME law. It admits both total and elastic strain
    /// component boxes and uses the original fixed reference volume.
    ///
    /// Continuous piecewise-affine instantaneous expansion is integrated on
    /// every crossed original source segment before strain subtraction.
    /// Coefficient increments use original engineering paths and exact rational
    /// or compliance-inverse identities, without subtracting rounded free-strain,
    /// stiffness or derivative endpoints. Work is a trapezoidal approximation;
    /// it does not assert entropy conservation or advance a caloric history.
    /// For multiple points at the same temperatures, reuse prepare_exchange.
    pub fn exchange(
        &self,
        before_temperature_k: f64,
        before_total_strain_mandel: [f64; 6],
        after_temperature_k: f64,
        after_total_strain_mandel: [f64; 6],
        reference_volume_m3: f64,
    ) -> Result<ThermoelasticExchange, ThermoelasticError> {
        self.prepare_exchange(before_temperature_k, after_temperature_k)?.evaluate(
            before_total_strain_mandel,
            after_total_strain_mandel,
            reference_volume_m3,
        )
    }
}

impl PreparedThermoelasticExchange {
    /// Original source resolved at the start of this interval.
    #[must_use]
    pub fn before_source(&self) -> &ResolvedThermoelastic {
        &self.before
    }

    /// Original source resolved at the end of this interval.
    #[must_use]
    pub fn after_source(&self) -> &ResolvedThermoelastic {
        self.changed.as_ref().unwrap_or(&self.before)
    }

    /// Evaluate one reference-volume weight with retained source increments.
    ///
    /// Total and elastic component boxes are admitted at BOTH endpoints.
    /// The seven retained extensive increments stay in double-double precision;
    /// callers may sum and scale them before rounding their final report.
    pub fn evaluate(
        &self,
        before_total_strain_mandel: [f64; 6],
        after_total_strain_mandel: [f64; 6],
        reference_volume_m3: f64,
    ) -> Result<ThermoelasticExchange, ThermoelasticError> {
        let before = self.before_source().state(
            before_total_strain_mandel, reference_volume_m3,
        )?;
        let after = self.after_source().state(
            after_total_strain_mandel, reference_volume_m3,
        )?;
        let half = dd(0.5);
        let total_increment: [Dd; 6] = core::array::from_fn(|i| {
            dd(after_total_strain_mandel[i]) - dd(before_total_strain_mandel[i])
        });
        // Do not recover this from rounded (total - free) endpoints: those
        // may be identical even when the original input increment is nonzero.
        let de: [Dd; 6] = core::array::from_fn(|i| {
            total_increment[i] - self.delta_free[i]
        });
        let e0 = before.elastic_strain_mandel.map(dd);
        let e1 = after.elastic_strain_mandel.map(dd);
        let mean_e = core::array::from_fn(|i| half * (e0[i] + e1[i]));
        let mean_stress: [Dd; 6] = core::array::from_fn(|i| {
            half * (dd(before.stress_mandel_pa[i]) + dd(after.stress_mandel_pa[i]))
        });
        let volume = dd(reference_volume_m3);
        let work = volume * dot(&mean_stress, &total_increment);
        let expansion_work = -volume * dot(&mean_stress, &self.delta_free);
        let stiffness_work = volume * half * self.dc.bilinear(&e0, &e1);
        let free_change = volume * dot(&mean_stress, &de) + stiffness_work;
        // Delta(alpha:sigma) - Delta(e:C_T:e)/2, polarized before summation.
        let entropy_change = volume * (
            dot(&self.delta_alpha, &mean_stress)
                + half * (self.c0.bilinear(&self.alpha, &de)
                    + self.c1.bilinear(&self.alpha, &de))
                + self.dc.bilinear(&self.alpha, &mean_e)
                - half * (
                    self.ct0.bilinear(&e0, &de) + self.ct1.bilinear(&e1, &de)
                        + self.dct.bilinear(&e0, &e1)
                )
        );
        let entropy = half * (dd(before.entropy_j_k) + dd(after.entropy_j_k));
        let internal_change = free_change
            + self.mean_temperature * entropy_change + entropy * self.delta_temperature;
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
