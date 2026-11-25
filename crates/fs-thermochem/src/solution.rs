//! Symmetric binary liquid miscibility from one molar mixing free energy.
//!
//! `g_mix/(RT) = x ln x + (1-x) ln(1-x) + chi(T) x(1-x)`, with
//! `chi(T) = A + B/T` by default. An explicitly supplied excess heat-capacity
//! coefficient C adds `C[(T-T0)-T ln(T/T0)] x(1-x)` to the molar Gibbs energy.
//! The same source then supplies excess enthalpy, entropy and heat capacity;
//! no independent mixture-Cp correction is added to the phase free energy.
//! This is an equal-size, non-electrolyte regular-solution model, not a prediction
//! from densities or material labels. The caller owns the pair fit and validity.
//! Pure-component reference chemical potentials cancel between liquid phases.
//! Pressure dependence, unequal molecular sizes, interfacial energy, nucleation
//! and spatial evolution are outside this constitutive model.
//!
//! The molar form follows the regular-solution free-energy implementation in
//! INL MOOSE (its `omega*c*(1-c)` interaction term):
//! <https://mooseframework.inl.gov/docs/doxygen/modules/classRegularSolutionFreeEnergy.html>.
//! Temperature integration uses the constant-Cp enthalpy/entropy primitives
//! documented by Cantera, applied here to the authored excess coefficient:
//! <https://www.cantera.org/stable/reference/thermo/species-thermo.html#constant-heat-capacity>.
//! The binodal is the common tangent; for this symmetric model it is horizontal.
//! A negative thermodynamic factor is a spinodal instability, never a negative
//! molecular diffusivity to feed an ordinary Fick evolution step.

use crate::{SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};
use fs_math::det::{exp, expm1, ln, sqrt};
use fs_qty::{Pressure, Temperature};

/// Pure-component calorics and closed isobaric thermal equilibrium.
pub mod caloric;

/// Homogeneous binary chemical potentials and explicit square-gradient energy.
pub mod gradient;

/// Source-owned multicomponent liquid activities and coupled diffusion forces.
pub mod multicomponent;

/// Asymmetric multicomponent NRTL excess-Gibbs liquid solution model.
pub mod nrtl;

/// Finite solid/liquid partition for an explicitly declared neutral solute.
pub mod solubility;

/// Asymmetric multicomponent UNIQUAC combinatorial-plus-residual excess-Gibbs model.
pub mod uniquac;

/// Maximum absolute sum error admitted for an independently represented pair.
/// Inputs within this floating-point roundoff budget are normalized together.
pub const BINARY_COMPOSITION_SUM_TOLERANCE_V1: f64 = 8.0 * f64::EPSILON;

/// Refusal of an invalid or unrepresentable liquid-mixture state.
#[derive(Clone, Debug, PartialEq)]
pub enum BinarySolutionErrorV1 {
    /// Invalid species, model parameters, composition or numerical controls.
    Invalid(&'static str),
    /// Temperature is outside the caller's declared source domain.
    OutsideTemperatureDomain,
    /// A finite, distinct pair of liquid compositions cannot be represented.
    Unrepresentable(&'static str),
    /// The declared iteration budget did not resolve the coexistence root.
    IterationLimit,
    /// Cancelled before publishing a complete result.
    Cancelled,
}
impl core::fmt::Display for BinarySolutionErrorV1 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "binary liquid solution: {self:?}")
    }
}
impl std::error::Error for BinarySolutionErrorV1 {}

/// Explicit bounded accuracy/work policy for the nontrivial binodal root.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BinaryCoexistenceOptionsV1 {
    /// Relative width of the root bracket in positive log-odds coordinates.
    pub relative_tolerance: f64,
    /// Maximum bisections, from 1 through 256.
    pub maximum_iterations: usize,
}
impl Default for BinaryCoexistenceOptionsV1 {
    fn default() -> Self {
        Self {
            relative_tolerance: 2e-13,
            maximum_iterations: 128,
        }
    }
}
impl BinaryCoexistenceOptionsV1 {
    /// Validate controls even for a state that requires no root solve.
    pub fn validate(self) -> Result<(), BinarySolutionErrorV1> {
        if !self.relative_tolerance.is_finite()
            || self.relative_tolerance < 4.0 * f64::EPSILON
            || self.relative_tolerance > 1e-3
            || !(1..=256).contains(&self.maximum_iterations)
        {
            return Err(BinarySolutionErrorV1::Invalid("coexistence controls"));
        }
        Ok(())
    }
}

/// Stability of the authored homogeneous feed, not the final equilibrated phases.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BinaryMixingStabilityV1 {
    /// Homogeneous composition lies on the stable free-energy envelope.
    Stable,
    /// Homogeneous composition is a coexistence endpoint.
    Binodal,
    /// Locally convex, but a two-liquid split has lower Gibbs energy.
    Metastable,
    /// Zero curvature on the boundary of the spinodal interval.
    SpinodalBoundary,
    /// Negative composition curvature; infinitesimal demixing lowers energy.
    Spinodal,
    /// The symmetric critical point, `chi=2`, `x=1/2`.
    Critical,
}

/// Immutable caller-declared pair law. Species order defines composition `x_A`.
#[derive(Clone, Debug, PartialEq)]
pub struct BinaryRegularSolutionV1 {
    species: [SpeciesId; 2],
    chi_constant: f64,
    chi_inverse_temperature_k: f64,
    excess_heat_capacity_j_per_mol_k: f64,
    excess_heat_capacity_reference_temperature: Temperature,
    minimum_temperature: Temperature,
    maximum_temperature: Temperature,
}
impl BinaryRegularSolutionV1 {
    /// Admit a pair law on an inclusive positive temperature domain.
    /// `chi_constant` is dimensionless; `chi_inverse_temperature_k` is in kelvin.
    /// Neither parameter is inferred from a density, phase label or solubility.
    pub fn new(
        species: [SpeciesId; 2],
        chi_constant: f64,
        chi_inverse_temperature_k: f64,
        minimum_temperature: Temperature,
        maximum_temperature: Temperature,
    ) -> Result<Self, BinarySolutionErrorV1> {
        if species[0] == species[1]
            || !chi_constant.is_finite()
            || !chi_inverse_temperature_k.is_finite()
            || !minimum_temperature.value().is_finite()
            || minimum_temperature.value() <= 0.0
            || !maximum_temperature.value().is_finite()
            || maximum_temperature.value() < minimum_temperature.value()
        {
            return Err(BinarySolutionErrorV1::Invalid(
                "pair law or temperature domain",
            ));
        }
        Ok(Self {
            species,
            chi_constant,
            chi_inverse_temperature_k,
            excess_heat_capacity_j_per_mol_k: 0.0,
            excess_heat_capacity_reference_temperature: minimum_temperature,
            minimum_temperature,
            maximum_temperature,
        })
    }

    /// Retain a signed constant coefficient C [J/(mol K)] for the homogeneous
    /// excess heat capacity `C x_A x_B`, with an explicit source datum T0.
    /// At T0 the original A/B interaction enthalpy and entropy are unchanged.
    /// Away from it, `H_ex/[x_A x_B] = R B + C (T-T0)` and
    /// `S_ex/[x_A x_B] = -R A + C ln(T/T0)` define the SAME Gibbs interaction.
    ///
    /// C may be negative: excess heat capacity need not be positive. A caloric
    /// consumer must separately admit positive total Cp with its original pure
    /// sources. The finite positive datum must lie inside this source domain.
    /// Supplying zero restores the original A+B/T law exactly.
    pub fn with_excess_heat_capacity(
        mut self,
        coefficient_j_per_mol_k: f64,
        reference_temperature: Temperature,
    ) -> Result<Self, BinarySolutionErrorV1> {
        if !coefficient_j_per_mol_k.is_finite() {
            return Err(BinarySolutionErrorV1::Invalid(
                "finite excess heat-capacity coefficient",
            ));
        }
        self.admitted_temperature_value(reference_temperature)?;
        self.excess_heat_capacity_j_per_mol_k = if coefficient_j_per_mol_k == 0.0 {
            0.0
        } else {
            coefficient_j_per_mol_k
        };
        self.excess_heat_capacity_reference_temperature = if coefficient_j_per_mol_k == 0.0 {
            self.minimum_temperature
        } else {
            reference_temperature
        };
        Ok(self)
    }

    /// Authored coefficient C of the homogeneous molar excess Cp `C x_A x_B`.
    pub const fn excess_heat_capacity_coefficient_j_per_mol_k(&self) -> f64 {
        self.excess_heat_capacity_j_per_mol_k
    }

    /// Original dimensionless A coefficient, without reconstructing a fit
    /// from sampled interaction values. Together with B, C, T0 and the source
    /// domain this identifies the complete shared liquid caloric convention.
    pub const fn chi_constant(&self) -> f64 {
        self.chi_constant
    }

    /// Original B coefficient in kelvin, in the declared species order.
    pub const fn chi_inverse_temperature_k(&self) -> f64 {
        self.chi_inverse_temperature_k
    }

    /// Original datum T0 of the authored excess heat capacity. For zero C the
    /// canonical datum is the domain's lower endpoint, since T0 has no effect.
    pub const fn excess_heat_capacity_reference_temperature(&self) -> Temperature {
        self.excess_heat_capacity_reference_temperature
    }

    /// Original component identities in the declared, never silently sorted order.
    pub const fn species(&self) -> &[SpeciesId; 2] {
        &self.species
    }

    /// Inclusive validity limits supplied with this pair law.
    pub const fn temperature_domain(&self) -> [Temperature; 2] {
        [self.minimum_temperature, self.maximum_temperature]
    }

    /// Dimensionless interaction parameter at the actual admitted temperature.
    pub fn chi(&self, temperature: Temperature) -> Result<f64, BinarySolutionErrorV1> {
        let t = self.admitted_temperature_value(temperature)?;
        let original = self.chi_constant + self.chi_inverse_temperature_k / t;
        if self.excess_heat_capacity_j_per_mol_k == 0.0 {
            return finite(original, "chi(T)");
        }
        let reference = self.excess_heat_capacity_reference_temperature.value();
        // This is (H_ex-T*S_ex)/(R*T), rearranged before evaluation. Near
        // the datum its C contribution starts at order (T-T0)^2; subtracting
        // separately rounded H and T*S would erase that entire response.
        let relative = (t - reference) / reference;
        let thermal = if relative.abs() < 0.125 {
            let mut remainder = 0.0;
            // r/(1+r)-ln(1+r) = sum_{n>=2} (-1)^(n-1)(n-1)r^n/n.
            // At |r|<=1/8 the omitted terms are below double roundoff.
            for n in (2..=20).rev() {
                let coefficient =
                    (f64::from(n - 1) / f64::from(n)) * if n % 2 == 0 { -1.0 } else { 1.0 };
                remainder = coefficient + relative * remainder;
            }
            relative * relative * remainder
        } else {
            (t - reference) / t - crate::entropy::log_ratio(t, reference)
        };
        finite(
            original + (self.excess_heat_capacity_j_per_mol_k / R) * thermal,
            "chi(T)",
        )
    }

    /// Coefficient of x_A*x_B in the original molar excess enthalpy [J/mol].
    pub(crate) fn interaction_enthalpy_j_per_mol(
        &self,
        temperature: Temperature,
    ) -> Result<f64, BinarySolutionErrorV1> {
        let t = self.admitted_temperature_value(temperature)?;
        finite(
            R * self.chi_inverse_temperature_k
                + self.excess_heat_capacity_j_per_mol_k
                    * (t - self.excess_heat_capacity_reference_temperature.value()),
            "interaction enthalpy coefficient",
        )
    }

    /// Coefficient of x_A*x_B in the original molar excess entropy [J/(mol K)].
    pub(crate) fn interaction_entropy_j_per_mol_k(
        &self,
        temperature: Temperature,
    ) -> Result<f64, BinarySolutionErrorV1> {
        let t = self.admitted_temperature_value(temperature)?;
        let thermal = if self.excess_heat_capacity_j_per_mol_k == 0.0 {
            0.0
        } else {
            self.excess_heat_capacity_j_per_mol_k
                * crate::entropy::log_ratio(
                    t,
                    self.excess_heat_capacity_reference_temperature.value(),
                )
        };
        finite(
            -R * self.chi_constant + thermal,
            "interaction entropy coefficient",
        )
    }

    /// Direct change of the excess-entropy coefficient, retaining nearby
    /// temperatures and exactly reversing without subtracting reference values.
    pub(crate) fn interaction_entropy_change_j_per_mol_k(
        &self,
        from: Temperature,
        to: Temperature,
    ) -> Result<f64, BinarySolutionErrorV1> {
        let from = self.admitted_temperature_value(from)?;
        let to = self.admitted_temperature_value(to)?;
        finite(
            self.excess_heat_capacity_j_per_mol_k * crate::entropy::log_ratio(to, from),
            "interaction entropy coefficient increment",
        )
    }

    fn admitted_temperature_value(
        &self,
        temperature: Temperature,
    ) -> Result<f64, BinarySolutionErrorV1> {
        let t = temperature.value();
        if !t.is_finite()
            || t < self.minimum_temperature.value()
            || t > self.maximum_temperature.value()
        {
            return Err(BinarySolutionErrorV1::OutsideTemperatureDomain);
        }
        Ok(t)
    }

    /// Evaluate the homogeneous free energy, even on a metastable/unstable branch.
    /// Exact pure states have zero mixing energy; absent-species chemical
    /// potentials remain absent rather than being clipped to a fictitious trace.
    pub fn evaluate(
        &self,
        temperature: Temperature,
        mole_fraction_a: f64,
    ) -> Result<BinarySolutionStateV1, BinarySolutionErrorV1> {
        self.evaluate_composition(temperature, [mole_fraction_a, 1.0 - mole_fraction_a])
    }

    /// Evaluate independently represented component mole fractions in species order.
    /// The sum must be within [`BINARY_COMPOSITION_SUM_TOLERANCE_V1`] of one.
    /// Both entries are normalized together; a real trace is never reconstructed
    /// by subtracting its rounded dominant component from one.
    pub fn evaluate_composition(
        &self,
        temperature: Temperature,
        mole_fractions: [f64; 2],
    ) -> Result<BinarySolutionStateV1, BinarySolutionErrorV1> {
        let fractions = normalized_fractions(mole_fractions)?;
        let chi = self.chi(temperature)?;
        self.state(temperature, fractions, chi)
    }

    fn state(
        &self,
        temperature: Temperature,
        fractions: [f64; 2],
        chi: f64,
    ) -> Result<BinarySolutionStateV1, BinarySolutionErrorV1> {
        let [x, y] = fractions;
        let log = [
            if x > 0.0 {
                Some(if x > 0.5 { ln_complement(y) } else { ln(x) })
            } else {
                None
            },
            if y > 0.0 {
                Some(if y > 0.5 { ln_complement(x) } else { ln(y) })
            } else {
                None
            },
        ];
        let ideal = x * log[0].unwrap_or(0.0) + y * log[1].unwrap_or(0.0);
        let xy = x * y;
        let activity = [
            finite(chi * y * y, "activity")?,
            finite(chi * x * x, "activity")?,
        ];
        let chemical = [
            log[0].map(|v| v + activity[0]),
            log[1].map(|v| v + activity[1]),
        ];
        if chemical.into_iter().flatten().any(|v| !v.is_finite()) {
            return Err(BinarySolutionErrorV1::Unrepresentable("chemical potential"));
        }
        Ok(BinarySolutionStateV1 {
            temperature,
            mole_fractions: fractions,
            chi,
            mixing_gibbs_j_per_mol: finite(
                (ideal + chi * xy) * R * temperature.value(),
                "mixing Gibbs energy",
            )?,
            mixing_enthalpy_j_per_mol: finite(
                if self.excess_heat_capacity_j_per_mol_k == 0.0 {
                    (self.chi_inverse_temperature_k * xy) * R
                } else {
                    self.interaction_enthalpy_j_per_mol(temperature)? * xy
                },
                "mixing enthalpy",
            )?,
            mixing_entropy_j_per_mol_k: finite(
                if self.excess_heat_capacity_j_per_mol_k == 0.0 {
                    -R * (ideal + self.chi_constant * xy)
                } else {
                    -R * ideal + self.interaction_entropy_j_per_mol_k(temperature)? * xy
                },
                "mixing entropy",
            )?,
            log_activity_coefficients: activity,
            reduced_chemical_potentials: chemical,
            // Equivalent for a normalized pair, but preserves the curvature
            // sign near chi=2 without subtracting two nearly equal numbers.
            thermodynamic_factor: finite(
                if chi == 0.0 {
                    1.0
                } else {
                    (x - y) * (x - y) - 2.0 * ((chi - 2.0) * xy)
                },
                "thermodynamic factor",
            )?,
        })
    }

    /// Resolve both coexisting compositions when `chi>2`; otherwise return `None`.
    /// The solve excludes the spurious central stationary root. No allocations
    /// occur, and each bounded iteration polls cancellation. Extreme gaps whose
    /// finite rich-side composition rounds to a pure state are refused.
    pub fn coexistence(
        &self,
        temperature: Temperature,
        options: BinaryCoexistenceOptionsV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Option<BinaryCoexistenceV1>, BinarySolutionErrorV1> {
        options.validate()?;
        if cancelled() {
            return Err(BinarySolutionErrorV1::Cancelled);
        }
        let chi = self.chi(temperature)?;
        if chi <= 2.0 {
            if cancelled() {
                return Err(BinarySolutionErrorV1::Cancelled);
            }
            return Ok(None);
        }
        // For u=ln((1-x)/x)>0, u*coth(u/2)-2 = chi-2 is strictly
        // increasing. Its small-u series avoids cancellation at criticality.
        let target = chi - 2.0;
        let mut bracket = [0.0, chi];
        for iteration in 1..=options.maximum_iterations {
            if cancelled() {
                return Err(BinarySolutionErrorV1::Cancelled);
            }
            let u = 0.5 * bracket[0] + 0.5 * bracket[1];
            if log_odds_curvature(u) < target {
                bracket[0] = u;
            } else {
                bracket[1] = u;
            }
            if bracket[1] - bracket[0] <= options.relative_tolerance * bracket[1] {
                let root = 0.5 * bracket[0] + 0.5 * bracket[1];
                let e = exp(-root);
                let x = e / (1.0 + e);
                let y = 1.0 - x;
                if !(x > 0.0 && x < 0.5 && y > 0.5 && y < 1.0) {
                    return Err(BinarySolutionErrorV1::Unrepresentable(
                        "distinct non-pure binodal compositions",
                    ));
                }
                let lower_a = self.state(temperature, [x, y], chi)?;
                // Retain the original dilute fraction: recomputing 1-y would
                // destroy relative accuracy in the rich phase.
                let upper_a = self.state(temperature, [y, x], chi)?;
                let spinodal_gap = sqrt((chi - 2.0) / chi);
                // Rationalized dilute spinodal endpoint avoids 1-sqrt(...).
                let lower_spinodal = (1.0 / chi) / (1.0 + spinodal_gap);
                if cancelled() {
                    return Err(BinarySolutionErrorV1::Cancelled);
                }
                return Ok(Some(BinaryCoexistenceV1 {
                    lower_a,
                    upper_a,
                    spinodal_mole_fractions: [lower_spinodal, 1.0 - lower_spinodal],
                    log_odds_bracket: bracket,
                    iterations: iteration,
                }));
            }
        }
        Err(BinarySolutionErrorV1::IterationLimit)
    }

    /// Equilibrate a homogeneous feed by the global common tangent and molar
    /// lever rule. The returned phase amounts refer to total MOLES; consumers
    /// must use original molecular masses and density laws for mass or volume.
    pub fn equilibrate(
        &self,
        temperature: Temperature,
        overall_mole_fraction_a: f64,
        options: BinaryCoexistenceOptionsV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<BinaryPhaseEquilibriumV1, BinarySolutionErrorV1> {
        let mut result = self.equilibrate_composition(
            temperature,
            [overall_mole_fraction_a, 1.0 - overall_mole_fraction_a],
            options,
            &mut cancelled,
        )?;
        // Preserve this scalar API's exact first-coordinate binodal endpoint.
        // Paired callers instead retain and compare the actual dilute component.
        if result.coexistence.as_ref().is_some_and(|coexistence| {
            overall_mole_fraction_a == coexistence.upper_a.mole_fractions[0]
        }) {
            result.phases[0] = BinaryEquilibriumPhaseV1 {
                state: result.feed,
                fraction_of_total_moles: 1.0,
            };
            result.phase_count = 1;
            result.feed_stability = BinaryMixingStabilityV1::Binodal;
        }
        if cancelled() {
            return Err(BinarySolutionErrorV1::Cancelled);
        }
        Ok(result)
    }

    /// Equilibrate independently represented, normalized component mole fractions.
    /// Trace inventories survive even when the other component rounds to one.
    /// Binodal comparisons and lever distances use each original minority
    /// fraction, not a subtraction from its rounded dominant counterpart.
    pub fn equilibrate_composition(
        &self,
        temperature: Temperature,
        overall_mole_fractions: [f64; 2],
        options: BinaryCoexistenceOptionsV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<BinaryPhaseEquilibriumV1, BinarySolutionErrorV1> {
        options.validate()?;
        if cancelled() {
            return Err(BinarySolutionErrorV1::Cancelled);
        }
        let feed = self.evaluate_composition(temperature, overall_mole_fractions)?;
        let [x, y] = feed.mole_fractions;
        let phase = BinaryEquilibriumPhaseV1 {
            state: feed,
            fraction_of_total_moles: 1.0,
        };
        let mut result = BinaryPhaseEquilibriumV1 {
            feed,
            phases: [phase; 2],
            phase_count: 1,
            coexistence: None,
            feed_stability: if feed.chi == 2.0 && x == 0.5 && y == 0.5 {
                BinaryMixingStabilityV1::Critical
            } else {
                BinaryMixingStabilityV1::Stable
            },
        };
        if x > 0.0
            && y > 0.0
            && let Some(coexistence) = self.coexistence(temperature, options, &mut cancelled)?
        {
            let dilute = coexistence.lower_a.mole_fractions[0];
            // hi-x equals y-dilute mathematically; only the latter retains an
            // independently represented dilute excess in a nearly pure feed.
            let low_distance = y - dilute;
            let high_distance = x - dilute;
            if low_distance == 0.0 || high_distance == 0.0 {
                result.feed_stability = BinaryMixingStabilityV1::Binodal;
            } else if low_distance > 0.0 && high_distance > 0.0 {
                let distance_sum = low_distance + high_distance;
                result.phases = [
                    BinaryEquilibriumPhaseV1 {
                        state: coexistence.lower_a,
                        fraction_of_total_moles: low_distance / distance_sum,
                    },
                    BinaryEquilibriumPhaseV1 {
                        state: coexistence.upper_a,
                        fraction_of_total_moles: high_distance / distance_sum,
                    },
                ];
                result.phase_count = 2;
                result.feed_stability = if feed.thermodynamic_factor < 0.0 {
                    BinaryMixingStabilityV1::Spinodal
                } else if feed.thermodynamic_factor == 0.0 {
                    BinaryMixingStabilityV1::SpinodalBoundary
                } else {
                    BinaryMixingStabilityV1::Metastable
                };
            }
            result.coexistence = Some(coexistence);
        }
        if cancelled() {
            return Err(BinarySolutionErrorV1::Cancelled);
        }
        Ok(result)
    }
}

/// Homogeneous composition response relative to its two pure-liquid references.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BinarySolutionStateV1 {
    temperature: Temperature,
    mole_fractions: [f64; 2],
    chi: f64,
    mixing_gibbs_j_per_mol: f64,
    mixing_enthalpy_j_per_mol: f64,
    mixing_entropy_j_per_mol_k: f64,
    log_activity_coefficients: [f64; 2],
    reduced_chemical_potentials: [Option<f64>; 2],
    thermodynamic_factor: f64,
}
impl BinarySolutionStateV1 {
    /// Actual absolute temperature [K].
    pub const fn temperature(&self) -> Temperature {
        self.temperature
    }
    /// Both mole fractions in the original species order.
    pub const fn mole_fractions(&self) -> [f64; 2] {
        self.mole_fractions
    }
    /// First declared component's mole fraction.
    pub const fn mole_fraction_a(&self) -> f64 {
        self.mole_fractions[0]
    }
    /// Dimensionless pair interaction evaluated at this temperature.
    pub const fn chi(&self) -> f64 {
        self.chi
    }
    /// Molar Gibbs energy of mixing [J/mol], with `0 ln 0 = 0`.
    pub const fn mixing_gibbs_j_per_mol(&self) -> f64 {
        self.mixing_gibbs_j_per_mol
    }
    /// Molar enthalpy of mixing [J/mol], `R B x(1-x)`.
    pub const fn mixing_enthalpy_j_per_mol(&self) -> f64 {
        self.mixing_enthalpy_j_per_mol
    }
    /// Molar entropy of mixing [J/(mol K)] from the same temperature law.
    pub const fn mixing_entropy_j_per_mol_k(&self) -> f64 {
        self.mixing_entropy_j_per_mol_k
    }
    /// Natural logarithms of activity coefficients in the pure-liquid convention.
    pub const fn log_activity_coefficients(&self) -> [f64; 2] {
        self.log_activity_coefficients
    }
    /// `(mu_i-mu_i_pure)/(RT)`. An absent component has divergent chemical
    /// potential and is represented by `None`, never a clipped logarithm.
    pub const fn reduced_chemical_potentials(&self) -> [Option<f64>; 2] {
        self.reduced_chemical_potentials
    }
    /// `Gamma=1-2 chi x(1-x)`. Its sign is the homogeneous curvature sign.
    pub const fn thermodynamic_factor(&self) -> f64 {
        self.thermodynamic_factor
    }
    /// Isothermal/isobaric Gibbs–Duhem-balanced forces [1/m] for the original
    /// Maxwell–Stefan owner: `[Gamma grad(x_A), -Gamma grad(x_A)]`.
    /// The exact pure endpoint uses the finite one-sided limit; no finite
    /// absent-species chemical potential is asserted. Reorder these forces
    /// with the species if the receiving transport owner uses canonical order.
    pub fn driving_forces(&self, gradient_a_per_m: f64) -> Result<[f64; 2], BinarySolutionErrorV1> {
        if !gradient_a_per_m.is_finite() {
            return Err(BinarySolutionErrorV1::Invalid("composition gradient"));
        }
        let force = finite(
            self.thermodynamic_factor * gradient_a_per_m,
            "driving force",
        )?;
        Ok([force, -force])
    }
}

/// Two finite coexistence compositions and their bounded numerical solve.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BinaryCoexistenceV1 {
    lower_a: BinarySolutionStateV1,
    upper_a: BinarySolutionStateV1,
    spinodal_mole_fractions: [f64; 2],
    log_odds_bracket: [f64; 2],
    iterations: usize,
}
impl BinaryCoexistenceV1 {
    /// Coexisting phase with less of the first component.
    pub const fn lower_a(&self) -> &BinarySolutionStateV1 {
        &self.lower_a
    }
    /// Coexisting phase with more of the first component.
    pub const fn upper_a(&self) -> &BinarySolutionStateV1 {
        &self.upper_a
    }
    /// Lower and upper first-component mole fractions at zero curvature.
    pub const fn spinodal_mole_fractions(&self) -> [f64; 2] {
        self.spinodal_mole_fractions
    }
    /// Final bracket for `ln((1-x_lower)/x_lower)`; numerical, not an interval certificate.
    pub const fn log_odds_bracket(&self) -> [f64; 2] {
        self.log_odds_bracket
    }
    /// Number of bounded root iterations used.
    pub const fn iterations(&self) -> usize {
        self.iterations
    }
}

/// One homogeneous equilibrium phase and its molar share of the whole inventory.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BinaryEquilibriumPhaseV1 {
    state: BinarySolutionStateV1,
    fraction_of_total_moles: f64,
}
impl BinaryEquilibriumPhaseV1 {
    /// Original thermodynamic composition response of this phase.
    pub const fn state(&self) -> &BinarySolutionStateV1 {
        &self.state
    }
    /// First declared component's mole fraction inside this phase.
    pub const fn mole_fraction_a(&self) -> f64 {
        self.state.mole_fractions[0]
    }
    /// Both component mole fractions inside this phase.
    pub const fn mole_fractions(&self) -> [f64; 2] {
        self.state.mole_fractions
    }
    /// Fraction of the entire inventory's total moles, not its mass or volume.
    pub const fn fraction_of_total_moles(&self) -> f64 {
        self.fraction_of_total_moles
    }
}

/// Homogeneous or two-liquid equilibrium of a closed binary molar inventory.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BinaryPhaseEquilibriumV1 {
    feed: BinarySolutionStateV1,
    phases: [BinaryEquilibriumPhaseV1; 2],
    phase_count: usize,
    coexistence: Option<BinaryCoexistenceV1>,
    feed_stability: BinaryMixingStabilityV1,
}
impl BinaryPhaseEquilibriumV1 {
    /// One or two positive-mole phases, ordered by their first-component fraction.
    pub fn phases(&self) -> &[BinaryEquilibriumPhaseV1] {
        &self.phases[..self.phase_count]
    }
    /// Original homogeneous feed, including any metastability/instability.
    pub const fn feed(&self) -> &BinarySolutionStateV1 {
        &self.feed
    }
    /// Free-energy classification of the original homogeneous feed.
    pub const fn feed_stability(&self) -> BinaryMixingStabilityV1 {
        self.feed_stability
    }
    /// Resolved coexistence gap when a non-pure feed required one.
    pub const fn coexistence(&self) -> Option<&BinaryCoexistenceV1> {
        self.coexistence.as_ref()
    }
    /// Equilibrium Gibbs mixing energy per mole of the whole inventory [J/mol].
    pub fn mixing_gibbs_j_per_mol(&self) -> f64 {
        self.phases()
            .iter()
            .map(|p| p.fraction_of_total_moles * p.state.mixing_gibbs_j_per_mol)
            .sum()
    }

    /// Equilibrium mixing enthalpy [J/mol], including the actual phase split.
    /// This is the excess contribution; pure-component enthalpies are separate.
    pub fn mixing_enthalpy_j_per_mol(&self) -> f64 {
        self.phases()
            .iter()
            .map(|p| p.fraction_of_total_moles * p.state.mixing_enthalpy_j_per_mol)
            .sum()
    }

    /// Equilibrium mixing entropy [J/(mol K)] from the same phase amounts.
    pub fn mixing_entropy_j_per_mol_k(&self) -> f64 {
        self.phases()
            .iter()
            .map(|p| p.fraction_of_total_moles * p.state.mixing_entropy_j_per_mol_k)
            .sum()
    }
}

fn finite(value: f64, quantity: &'static str) -> Result<f64, BinarySolutionErrorV1> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(BinarySolutionErrorV1::Unrepresentable(quantity))
    }
}
fn fraction(value: f64) -> Result<(), BinarySolutionErrorV1> {
    if value.is_finite() && (0.0..=1.0).contains(&value) {
        Ok(())
    } else {
        Err(BinarySolutionErrorV1::Invalid(
            "mole fraction outside [0,1]",
        ))
    }
}
fn normalized_fractions(values: [f64; 2]) -> Result<[f64; 2], BinarySolutionErrorV1> {
    fraction(values[0])?;
    fraction(values[1])?;
    let sum = values[0] + values[1];
    if (sum - 1.0).abs() > BINARY_COMPOSITION_SUM_TOLERANCE_V1 {
        return Err(BinarySolutionErrorV1::Invalid(
            "component mole fractions do not sum to one",
        ));
    }
    Ok([values[0] / sum, values[1] / sum])
}
fn ln_complement(x: f64) -> f64 {
    // The other component may round to one while x remains a real trace.
    // Keep its -x contribution to mixing entropy and chemical potential.
    if x < 1e-4 {
        -x * (1.0 + x * (0.5 + x * (1.0 / 3.0 + x * (0.25 + x * 0.2))))
    } else {
        ln(1.0 - x)
    }
}
fn log_odds_curvature(u: f64) -> f64 {
    if u <= 0.25 {
        let u2 = u * u;
        u2 * (1.0 / 6.0
            + u2 * (-1.0 / 360.0
                + u2 * (1.0 / 15120.0
                    + u2 * (-1.0 / 604800.0
                        + u2 * (1.0 / 23950080.0 - u2 * (691.0 / 653837184000.0))))))
    } else {
        // -expm1(-u) is 1-exp(-u), without a close subtraction.
        let one_minus_e = -expm1(-u);
        u * ((2.0 - one_minus_e) / one_minus_e) - 2.0
    }
}

#[cfg(test)]
mod tests;

#[cfg(test)]
mod excess_heat_capacity_tests;

// ============================================================================
// Asymmetric Multicomponent Excess-Gibbs Solutions & Modified Raoult's Law VLE
// ============================================================================

/// Maximum species count admitted by asymmetric excess-Gibbs models (`NRTL`, `UNIQUAC`).
pub const MAX_EXCESS_GIBBS_SPECIES_V1: usize = 64;

/// Maximum admitted deviation of raw mole-fraction sum from unity before normalization.
pub const EXCESS_GIBBS_COMPOSITION_SUM_TOLERANCE_V1: f64 = 64.0 * f64::EPSILON;

/// Maximum outer iterations permitted in bounded VLE and flash solvers.
pub const MAX_VLE_SOLVER_ITERATIONS_V1: usize = 512;

/// Dedicated typed refusal for asymmetric excess-Gibbs models and Modified Raoult VLE.
#[derive(Clone, Debug, PartialEq)]
pub enum ExcessGibbsErrorV1 {
    /// Invalid species list, duplicate species, mismatched dimensions, or invalid parameters.
    InvalidModel(&'static str),
    /// Temperature is non-finite, non-positive, or outside the declared validity domain.
    TemperatureOutOfDomain {
        /// Rejected temperature [K].
        temperature_k: f64,
        /// Inclusive lower validity bound [K].
        min_k: f64,
        /// Inclusive upper validity bound [K].
        max_k: f64,
    },
    /// Pressure is non-finite or non-positive.
    InvalidPressure {
        /// Rejected pressure [Pa].
        pressure_pa: f64,
    },
    /// Composition vector has wrong length, out-of-range/non-finite entries, or sum != 1.
    InvalidComposition(&'static str),
    /// Invalid solver options or vapor-pressure law parameter.
    InvalidOptions(&'static str),
    /// Intermediate or output thermodynamic value overflowed, underflowed, or became non-finite.
    Unrepresentable(&'static str),
    /// Bounded VLE/flash solver did not converge or root was unbracketed in the domain.
    ConvergenceFailed {
        /// Number of iterations executed before refusal.
        iterations: usize,
        /// Final dimensionless or relative pressure residual.
        residual: f64,
    },
    /// Caller cancellation observed before publishing a complete result.
    Cancelled,
}

impl core::fmt::Display for ExcessGibbsErrorV1 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidModel(msg) => write!(f, "invalid excess-Gibbs model: {msg}"),
            Self::TemperatureOutOfDomain {
                temperature_k,
                min_k,
                max_k,
            } => write!(
                f,
                "temperature {temperature_k} K outside excess-Gibbs domain [{min_k}, {max_k}] K"
            ),
            Self::InvalidPressure { pressure_pa } => {
                write!(
                    f,
                    "pressure {pressure_pa} Pa must be finite and strictly positive"
                )
            }
            Self::InvalidComposition(msg) => write!(f, "invalid liquid/vapor composition: {msg}"),
            Self::InvalidOptions(msg) => write!(f, "invalid VLE solver options: {msg}"),
            Self::Unrepresentable(msg) => write!(f, "unrepresentable excess-Gibbs quantity: {msg}"),
            Self::ConvergenceFailed {
                iterations,
                residual,
            } => write!(
                f,
                "VLE solver failed to converge after {iterations} iterations (residual {residual:e})"
            ),
            Self::Cancelled => f.write_str("excess-Gibbs VLE calculation cancelled"),
        }
    }
}

impl core::error::Error for ExcessGibbsErrorV1 {}

/// Exact single-potential evaluation returned by asymmetric excess-Gibbs models (`NRTL`, `UNIQUAC`).
#[derive(Clone, Debug, PartialEq)]
pub struct ExcessGibbsEvaluationV1 {
    pub(crate) temperature: Temperature,
    pub(crate) mole_fractions: Vec<f64>,
    pub(crate) reduced_excess_gibbs: f64,
    pub(crate) excess_gibbs_j_per_mol: f64,
    pub(crate) excess_enthalpy_j_per_mol: f64,
    pub(crate) excess_entropy_j_per_mol_k: f64,
    pub(crate) ln_gamma: Vec<f64>,
    pub(crate) gamma: Vec<f64>,
    pub(crate) d_ln_gamma_dt: Vec<f64>,
    pub(crate) partial_excess_enthalpy_j_per_mol: Vec<f64>,
    pub(crate) d_ln_gamma_dx: Vec<f64>,
}

impl ExcessGibbsEvaluationV1 {
    /// Evaluated temperature [K].
    #[must_use]
    pub const fn temperature(&self) -> Temperature {
        self.temperature
    }

    /// Normalized liquid mole fractions in declared species order.
    #[must_use]
    pub fn mole_fractions(&self) -> &[f64] {
        &self.mole_fractions
    }

    /// Dimensionless reduced molar excess Gibbs energy `g^E / (R T)`.
    #[must_use]
    pub const fn reduced_excess_gibbs(&self) -> f64 {
        self.reduced_excess_gibbs
    }

    /// Molar excess Gibbs energy `g^E` [J/mol].
    #[must_use]
    pub const fn excess_gibbs_j_per_mol(&self) -> f64 {
        self.excess_gibbs_j_per_mol
    }

    /// Molar excess enthalpy `h^E = -R T^2 d(g^E/(RT))/dT` [J/mol].
    #[must_use]
    pub const fn excess_enthalpy_j_per_mol(&self) -> f64 {
        self.excess_enthalpy_j_per_mol
    }

    /// Molar excess entropy `s^E = (h^E - g^E) / T` [J/(mol K)].
    #[must_use]
    pub const fn excess_entropy_j_per_mol_k(&self) -> f64 {
        self.excess_entropy_j_per_mol_k
    }

    /// Natural logarithms of activity coefficients `ln gamma_i`.
    #[must_use]
    pub fn ln_activity_coefficients(&self) -> &[f64] {
        &self.ln_gamma
    }

    /// Activity coefficients `gamma_i = exp(ln gamma_i)`.
    #[must_use]
    pub fn activity_coefficients(&self) -> &[f64] {
        &self.gamma
    }

    /// Analytical temperature derivatives `d(ln gamma_i)/dT` [1/K].
    #[must_use]
    pub fn d_ln_gamma_dt(&self) -> &[f64] {
        &self.d_ln_gamma_dt
    }

    /// Partial molar excess enthalpies `h_i^E = -R T^2 d(ln gamma_i)/dT` [J/mol].
    #[must_use]
    pub fn partial_excess_enthalpy_j_per_mol(&self) -> &[f64] {
        &self.partial_excess_enthalpy_j_per_mol
    }

    /// Row-major `n x n` symmetric homogeneous-degree-0 Jacobian `d(ln gamma_i)/dx_j`.
    #[must_use]
    pub fn d_ln_gamma_dx(&self) -> &[f64] {
        &self.d_ln_gamma_dx
    }

    /// Look up entry `(i, j)` of the Jacobian `d(ln gamma_i)/dx_j`.
    #[must_use]
    pub fn d_ln_gamma_dx_entry(&self, i: usize, j: usize) -> Option<f64> {
        let n = self.mole_fractions.len();
        if i < n && j < n {
            Some(self.d_ln_gamma_dx[i * n + j])
        } else {
            None
        }
    }

    /// Compute each Gibbs-Duhem column residual `sum_i x_i * d(ln gamma_i)/dx_j`.
    ///
    /// For any exact single-potential model this vanishes to floating-point roundoff.
    #[must_use]
    pub fn gibbs_duhem_column_residuals(&self) -> Vec<f64> {
        let n = self.mole_fractions.len();
        let mut residuals = Vec::with_capacity(n);
        for j in 0..n {
            let mut sum = 0.0;
            let mut comp = 0.0;
            for i in 0..n {
                neumaier_add(
                    &mut sum,
                    &mut comp,
                    self.mole_fractions[i] * self.d_ln_gamma_dx[i * n + j],
                );
            }
            residuals.push(sum + comp);
        }
        residuals
    }

    /// Compute the directional Gibbs-Duhem residual `sum_i x_i * (sum_j (d ln gamma_i / dx_j) * delta_x_j)`.
    pub fn gibbs_duhem_residual(&self, delta_x: &[f64]) -> Result<f64, ExcessGibbsErrorV1> {
        let n = self.mole_fractions.len();
        if delta_x.len() != n || delta_x.iter().any(|v| !v.is_finite()) {
            return Err(ExcessGibbsErrorV1::InvalidComposition(
                "Gibbs-Duhem direction must have length n with finite entries",
            ));
        }
        let mut total = 0.0;
        let mut comp = 0.0;
        for i in 0..n {
            let mut d_lng_i = 0.0;
            let mut c_i = 0.0;
            for j in 0..n {
                neumaier_add(
                    &mut d_lng_i,
                    &mut c_i,
                    self.d_ln_gamma_dx[i * n + j] * delta_x[j],
                );
            }
            neumaier_add(
                &mut total,
                &mut comp,
                self.mole_fractions[i] * (d_lng_i + c_i),
            );
        }
        Ok(total + comp)
    }
}

/// Shared contract for asymmetric excess-Gibbs models driving Modified Raoult VLE.
pub trait ExcessGibbsModelV1 {
    /// Declared species in row/column order.
    fn species(&self) -> &[SpeciesId];
    /// Inclusive valid temperature interval `[T_min, T_max]` [K].
    fn temperature_domain(&self) -> [Temperature; 2];
    /// Evaluate excess Gibbs energy, activity coefficients, and exact analytical derivatives.
    fn evaluate(
        &self,
        temperature: Temperature,
        mole_fractions: &[f64],
    ) -> Result<ExcessGibbsEvaluationV1, ExcessGibbsErrorV1>;
}

/// Pure-component saturation vapor pressure evaluation `(P^sat, dP^sat/dT)` in SI units.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct VaporPressurePointV1 {
    /// Saturation vapor pressure `P^sat(T)` [Pa] (`0.0` for non-volatile solutes).
    pub pressure_pa: f64,
    /// Analytical derivative `dP^sat/dT` [Pa/K] (`0.0` for non-volatile solutes).
    pub dp_dt_pa_per_k: f64,
}

/// Pure-component saturation vapor-pressure correlation for Modified Raoult VLE.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum VaporPressureLawV1 {
    /// SI natural-log Antoine equation: `ln(P^sat / Pa) = a - b_k / (T + c_k)`.
    Antoine {
        /// Dimensionless intercept `A` in `ln(P/Pa)`.
        a: f64,
        /// Positive slope coefficient `B` [K].
        b_k: f64,
        /// Temperature shift `C` [K] (`T + C > 0` across validity interval).
        c_k: f64,
        /// Inclusive lower temperature bound [K].
        temperature_min_k: f64,
        /// Inclusive upper temperature bound [K].
        temperature_max_k: f64,
    },
    /// Wagner 2.5-5 equation: `ln(P^sat / P_c) = (a*tau + b*tau^1.5 + c*tau^2.5 + d*tau^5) / T_r`,
    /// where `T_r = T / T_c` and `tau = 1 - T_r`.
    Wagner25_5 {
        /// Critical temperature `T_c` [K].
        critical_temperature_k: f64,
        /// Critical pressure `P_c` [Pa].
        critical_pressure_pa: f64,
        /// Linear `tau` coefficient `a`.
        a: f64,
        /// `tau^1.5` coefficient `b`.
        b: f64,
        /// `tau^2.5` coefficient `c`.
        c: f64,
        /// `tau^5` coefficient `d`.
        d: f64,
        /// Inclusive lower temperature bound [K].
        temperature_min_k: f64,
        /// Inclusive upper temperature bound [K] (`<= T_c`).
        temperature_max_k: f64,
    },
    /// Wagner 3-6 equation: `ln(P^sat / P_c) = (a*tau + b*tau^1.5 + c*tau^3 + d*tau^6) / T_r`.
    Wagner3_6 {
        /// Critical temperature `T_c` [K].
        critical_temperature_k: f64,
        /// Critical pressure `P_c` [Pa].
        critical_pressure_pa: f64,
        /// Linear `tau` coefficient `a`.
        a: f64,
        /// `tau^1.5` coefficient `b`.
        b: f64,
        /// `tau^3` coefficient `c`.
        c: f64,
        /// `tau^6` coefficient `d`.
        d: f64,
        /// Inclusive lower temperature bound [K].
        temperature_min_k: f64,
        /// Inclusive upper temperature bound [K] (`<= T_c`).
        temperature_max_k: f64,
    },
    /// IAPWS SR1-86(1992) ordinary water saturation correlation (`273.16..=647.096 K`).
    IapwsWater,
    /// Non-volatile solute (`P^sat = 0`, `dP^sat/dT = 0`, e.g., sucrose or dissolved salt).
    NonVolatile,
}

impl VaporPressureLawV1 {
    /// Construct a validated SI natural-log Antoine law `ln(P/Pa) = a - b_k / (T + c_k)`.
    pub fn antoine_si_ln_pa(
        a: f64,
        b_k: f64,
        c_k: f64,
        temperature_min: Temperature,
        temperature_max: Temperature,
    ) -> Result<Self, ExcessGibbsErrorV1> {
        validate_temperature_domain(temperature_min, temperature_max)?;
        let t_min = temperature_min.value();
        let t_max = temperature_max.value();
        if !a.is_finite() || !b_k.is_finite() || b_k <= 0.0 || !c_k.is_finite() {
            return Err(ExcessGibbsErrorV1::InvalidOptions(
                "Antoine parameters must be finite with B > 0",
            ));
        }
        if t_min + c_k <= 0.0 {
            return Err(ExcessGibbsErrorV1::InvalidOptions(
                "Antoine denominator T_min + C must be strictly positive",
            ));
        }
        Ok(Self::Antoine {
            a,
            b_k,
            c_k,
            temperature_min_k: t_min,
            temperature_max_k: t_max,
        })
    }

    /// Construct an Antoine law from published DECHEMA base-10 mmHg / Celsius constants:
    /// `log10(P_mmHg) = a10 - b10 / (T_C + c_celsius)`.
    pub fn antoine_base10_mmhg_celsius(
        a10: f64,
        b10: f64,
        c_celsius: f64,
        temperature_min: Temperature,
        temperature_max: Temperature,
    ) -> Result<Self, ExcessGibbsErrorV1> {
        if !a10.is_finite() || !b10.is_finite() || !c_celsius.is_finite() {
            return Err(ExcessGibbsErrorV1::InvalidOptions(
                "DECHEMA Antoine parameters must be finite",
            ));
        }
        let ln10 = core::f64::consts::LN_10;
        let ln_pa_per_mmhg = ln(101_325.0 / 760.0);
        Self::antoine_si_ln_pa(
            ln_pa_per_mmhg + a10 * ln10,
            b10 * ln10,
            c_celsius - 273.15,
            temperature_min,
            temperature_max,
        )
    }

    /// Construct an Antoine law from published NIST WebBook base-10 bar / Kelvin constants:
    /// `log10(P_bar) = a10 - b10 / (T_K + c_kelvin)`.
    pub fn antoine_base10_bar_kelvin(
        a10: f64,
        b10: f64,
        c_kelvin: f64,
        temperature_min: Temperature,
        temperature_max: Temperature,
    ) -> Result<Self, ExcessGibbsErrorV1> {
        if !a10.is_finite() || !b10.is_finite() || !c_kelvin.is_finite() {
            return Err(ExcessGibbsErrorV1::InvalidOptions(
                "NIST Antoine parameters must be finite",
            ));
        }
        let ln10 = core::f64::consts::LN_10;
        Self::antoine_si_ln_pa(
            (a10 + 5.0) * ln10,
            b10 * ln10,
            c_kelvin,
            temperature_min,
            temperature_max,
        )
    }

    /// Construct a validated Wagner 2.5-5 vapor-pressure correlation.
    #[allow(clippy::too_many_arguments)]
    pub fn wagner_25_5(
        critical_temperature: Temperature,
        critical_pressure: Pressure,
        a: f64,
        b: f64,
        c: f64,
        d: f64,
        temperature_min: Temperature,
        temperature_max: Temperature,
    ) -> Result<Self, ExcessGibbsErrorV1> {
        validate_wagner_inputs(
            critical_temperature,
            critical_pressure,
            [a, b, c, d],
            temperature_min,
            temperature_max,
        )?;
        Ok(Self::Wagner25_5 {
            critical_temperature_k: critical_temperature.value(),
            critical_pressure_pa: critical_pressure.value(),
            a,
            b,
            c,
            d,
            temperature_min_k: temperature_min.value(),
            temperature_max_k: temperature_max.value(),
        })
    }

    /// Construct a validated Wagner 3-6 vapor-pressure correlation.
    #[allow(clippy::too_many_arguments)]
    pub fn wagner_3_6(
        critical_temperature: Temperature,
        critical_pressure: Pressure,
        a: f64,
        b: f64,
        c: f64,
        d: f64,
        temperature_min: Temperature,
        temperature_max: Temperature,
    ) -> Result<Self, ExcessGibbsErrorV1> {
        validate_wagner_inputs(
            critical_temperature,
            critical_pressure,
            [a, b, c, d],
            temperature_min,
            temperature_max,
        )?;
        Ok(Self::Wagner3_6 {
            critical_temperature_k: critical_temperature.value(),
            critical_pressure_pa: critical_pressure.value(),
            a,
            b,
            c,
            d,
            temperature_min_k: temperature_min.value(),
            temperature_max_k: temperature_max.value(),
        })
    }

    /// IAPWS SR1-86(1992) ordinary water saturation vapor pressure.
    #[must_use]
    pub const fn iapws_water() -> Self {
        Self::IapwsWater
    }

    /// Non-volatile solute (`P^sat = 0` identically).
    #[must_use]
    pub const fn non_volatile() -> Self {
        Self::NonVolatile
    }

    /// NIST Chemistry WebBook Antoine constants for pure ethanol (valid `273.15..=423.15 K`).
    #[must_use]
    pub fn ethanol_antoine_nist() -> Self {
        Self::antoine_base10_bar_kelvin(
            5.24677,
            1598.673,
            -46.424,
            Temperature::new(273.15),
            Temperature::new(423.15),
        )
        .expect("published NIST ethanol Antoine parameters are valid")
    }

    /// DECHEMA (Gmehling & Onken) Antoine constants for pure ethanol (`273.15..=423.15 K`).
    #[must_use]
    pub fn ethanol_antoine_dechema() -> Self {
        Self::antoine_base10_mmhg_celsius(
            8.11220,
            1592.864,
            226.184,
            Temperature::new(273.15),
            Temperature::new(423.15),
        )
        .expect("published DECHEMA ethanol Antoine parameters are valid")
    }

    /// DECHEMA (Gmehling & Onken) Antoine constants for pure water (`273.15..=423.15 K`).
    #[must_use]
    pub fn water_antoine_dechema() -> Self {
        Self::antoine_base10_mmhg_celsius(
            8.07131,
            1730.63,
            233.426,
            Temperature::new(273.15),
            Temperature::new(423.15),
        )
        .expect("published DECHEMA water Antoine parameters are valid")
    }

    /// Whether this component can enter the vapor phase (`false` only for `NonVolatile`).
    #[must_use]
    pub const fn is_volatile(self) -> bool {
        !matches!(self, Self::NonVolatile)
    }

    /// Inclusive temperature validity domain `[T_min, T_max]` [K], or `None` for `NonVolatile`.
    #[must_use]
    pub const fn temperature_domain(self) -> Option<[f64; 2]> {
        match self {
            Self::Antoine {
                temperature_min_k,
                temperature_max_k,
                ..
            }
            | Self::Wagner25_5 {
                temperature_min_k,
                temperature_max_k,
                ..
            }
            | Self::Wagner3_6 {
                temperature_min_k,
                temperature_max_k,
                ..
            } => Some([temperature_min_k, temperature_max_k]),
            Self::IapwsWater => Some([
                crate::water::TRIPLE_POINT_TEMPERATURE_K,
                crate::water::CRITICAL_TEMPERATURE_K,
            ]),
            Self::NonVolatile => None,
        }
    }

    /// Evaluate `(P^sat(T), dP^sat/dT)` at `temperature`.
    pub fn evaluate(
        self,
        temperature: Temperature,
    ) -> Result<VaporPressurePointV1, ExcessGibbsErrorV1> {
        let t = temperature.value();
        if !t.is_finite() || t <= 0.0 {
            return Err(ExcessGibbsErrorV1::TemperatureOutOfDomain {
                temperature_k: t,
                min_k: 0.0,
                max_k: f64::INFINITY,
            });
        }
        match self {
            Self::NonVolatile => Ok(VaporPressurePointV1 {
                pressure_pa: 0.0,
                dp_dt_pa_per_k: 0.0,
            }),
            Self::IapwsWater => {
                let sat = crate::water::SaturatedWaterState::at_temperature(temperature).map_err(
                    |_| ExcessGibbsErrorV1::TemperatureOutOfDomain {
                        temperature_k: t,
                        min_k: crate::water::TRIPLE_POINT_TEMPERATURE_K,
                        max_k: crate::water::CRITICAL_TEMPERATURE_K,
                    },
                )?;
                Ok(VaporPressurePointV1 {
                    pressure_pa: sat.pressure().value(),
                    dp_dt_pa_per_k: sat.pressure_derivative().value(),
                })
            }
            Self::Antoine {
                a,
                b_k,
                c_k,
                temperature_min_k,
                temperature_max_k,
            } => {
                check_temperature_in_domain(temperature, temperature_min_k, temperature_max_k)?;
                let denom = t + c_k;
                if !denom.is_finite() || denom <= 0.0 {
                    return Err(ExcessGibbsErrorV1::Unrepresentable(
                        "Antoine denominator T + C is non-positive",
                    ));
                }
                let ln_p = check_finite(a - b_k / denom, "Antoine ln(P^sat)")?;
                let p_sat = exp(ln_p);
                let dp_dt = p_sat * (b_k / (denom * denom));
                if !p_sat.is_finite() || p_sat <= 0.0 || !dp_dt.is_finite() {
                    return Err(ExcessGibbsErrorV1::Unrepresentable(
                        "Antoine vapor pressure overflowed or underflowed",
                    ));
                }
                Ok(VaporPressurePointV1 {
                    pressure_pa: p_sat,
                    dp_dt_pa_per_k: dp_dt,
                })
            }
            Self::Wagner25_5 {
                critical_temperature_k,
                critical_pressure_pa,
                a,
                b,
                c,
                d,
                temperature_min_k,
                temperature_max_k,
            } => {
                check_temperature_in_domain(temperature, temperature_min_k, temperature_max_k)?;
                evaluate_wagner(
                    t,
                    critical_temperature_k,
                    critical_pressure_pa,
                    [a, b, c, d],
                    false,
                )
            }
            Self::Wagner3_6 {
                critical_temperature_k,
                critical_pressure_pa,
                a,
                b,
                c,
                d,
                temperature_min_k,
                temperature_max_k,
            } => {
                check_temperature_in_domain(temperature, temperature_min_k, temperature_max_k)?;
                evaluate_wagner(
                    t,
                    critical_temperature_k,
                    critical_pressure_pa,
                    [a, b, c, d],
                    true,
                )
            }
        }
    }
}

fn validate_wagner_inputs(
    critical_temperature: Temperature,
    critical_pressure: Pressure,
    coeffs: [f64; 4],
    temperature_min: Temperature,
    temperature_max: Temperature,
) -> Result<(), ExcessGibbsErrorV1> {
    validate_temperature_domain(temperature_min, temperature_max)?;
    let tc = critical_temperature.value();
    let pc = critical_pressure.value();
    if !tc.is_finite() || tc <= 0.0 || !pc.is_finite() || pc <= 0.0 {
        return Err(ExcessGibbsErrorV1::InvalidOptions(
            "Wagner critical T_c and P_c must be finite and positive",
        ));
    }
    if temperature_max.value() > tc {
        return Err(ExcessGibbsErrorV1::InvalidOptions(
            "Wagner temperature_max must not exceed critical temperature T_c",
        ));
    }
    if coeffs.iter().any(|v| !v.is_finite()) {
        return Err(ExcessGibbsErrorV1::InvalidOptions(
            "Wagner coefficients must be finite",
        ));
    }
    Ok(())
}

fn evaluate_wagner(
    t: f64,
    tc: f64,
    pc: f64,
    coeffs: [f64; 4],
    is_3_6: bool,
) -> Result<VaporPressurePointV1, ExcessGibbsErrorV1> {
    let [a, b, c, d] = coeffs;
    let tr = t / tc;
    let tau = (tc - t) / tc;
    if !tr.is_finite() || tr <= 0.0 || !tau.is_finite() || tau < 0.0 {
        return Err(ExcessGibbsErrorV1::Unrepresentable(
            "Wagner reduced temperature outside (0, 1]",
        ));
    }
    let s = sqrt(tau);
    let tau_15 = tau * s;
    let tau2 = tau * tau;
    let (f, fp) = if is_3_6 {
        let tau3 = tau2 * tau;
        let tau6 = tau3 * tau3;
        (
            a * tau + b * tau_15 + c * tau3 + d * tau6,
            a + 1.5 * b * s + 3.0 * c * tau2 + 6.0 * d * (tau3 * tau2),
        )
    } else {
        let tau25 = tau2 * s;
        let tau5 = tau2 * tau2 * tau;
        (
            a * tau + b * tau_15 + c * tau25 + d * tau5,
            a + 1.5 * b * s + 2.5 * c * tau_15 + 5.0 * d * (tau2 * tau2),
        )
    };
    let ln_pr = check_finite(f / tr, "Wagner ln(P/P_c)")?;
    let p_sat = pc * exp(ln_pr);
    let dlnp_dt = -(fp * tr + f) / (tc * tr * tr);
    let dp_dt = p_sat * dlnp_dt;
    if !p_sat.is_finite() || p_sat <= 0.0 || !dp_dt.is_finite() {
        return Err(ExcessGibbsErrorV1::Unrepresentable(
            "Wagner vapor pressure overflowed or underflowed",
        ));
    }
    Ok(VaporPressurePointV1 {
        pressure_pa: p_sat,
        dp_dt_pa_per_k: dp_dt,
    })
}

/// Explicit tolerances and iteration cap for Modified Raoult VLE solvers.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct VleSolverOptionsV1 {
    /// Absolute temperature bracket convergence tolerance [K].
    pub temperature_tolerance_k: f64,
    /// Relative pressure / equilibrium convergence tolerance (dimensionless).
    pub relative_residual_tolerance: f64,
    /// Maximum solver iterations (`1..=MAX_VLE_SOLVER_ITERATIONS_V1`).
    pub maximum_iterations: usize,
}

impl Default for VleSolverOptionsV1 {
    fn default() -> Self {
        Self {
            temperature_tolerance_k: 1e-11,
            relative_residual_tolerance: 1e-12,
            maximum_iterations: 128,
        }
    }
}

impl VleSolverOptionsV1 {
    /// Validate numerical tolerance and iteration limits before running a solve.
    pub fn validate(self) -> Result<(), ExcessGibbsErrorV1> {
        if !self.temperature_tolerance_k.is_finite()
            || self.temperature_tolerance_k <= 0.0
            || !self.relative_residual_tolerance.is_finite()
            || self.relative_residual_tolerance <= 0.0
            || !(1..=MAX_VLE_SOLVER_ITERATIONS_V1).contains(&self.maximum_iterations)
        {
            return Err(ExcessGibbsErrorV1::InvalidOptions(
                "VLE solver tolerances must be finite positive and iterations in 1..=512",
            ));
        }
        Ok(())
    }
}

/// Modified Raoult's law K-values evaluated at `(T, P, x)`.
#[derive(Clone, Debug, PartialEq)]
pub struct VleKValuesV1 {
    /// Evaluated temperature [K].
    pub temperature: Temperature,
    /// Evaluated total pressure [Pa].
    pub pressure: Pressure,
    /// Normalized liquid mole fractions `x_i`.
    pub liquid_x: Vec<f64>,
    /// Equilibrium partition ratios `K_i = gamma_i(T, x) * P_i^sat(T) / P`.
    pub k_values: Vec<f64>,
    /// Normalized incipient vapor composition `y_i = K_i x_i / sum_j(K_j x_j)`.
    pub vapor_y_normalized: Vec<f64>,
    /// Bubble-point sum `sum_i K_i x_i` (`1.0` exactly at the bubble point).
    pub sum_kx: f64,
    /// Underlying excess-Gibbs liquid evaluation.
    pub excess_gibbs: ExcessGibbsEvaluationV1,
}

/// Converged bubble-point solution `(T_bub, y)` at caller-specified pressure `P` and liquid `x`.
#[derive(Clone, Debug, PartialEq)]
pub struct BubblePointV1 {
    /// Bubble-point temperature `T_bub` [K].
    pub temperature: Temperature,
    /// Specified total pressure `P` [Pa].
    pub pressure: Pressure,
    /// Normalized liquid mole fractions `x`.
    pub liquid_x: Vec<f64>,
    /// Equilibrium vapor mole fractions `y_i = x_i gamma_i P_i^sat / P`.
    pub vapor_y: Vec<f64>,
    /// Partition ratios `K_i = gamma_i P_i^sat / P`.
    pub k_values: Vec<f64>,
    /// Underlying excess-Gibbs evaluation at `T_bub`.
    pub excess_gibbs: ExcessGibbsEvaluationV1,
    /// Final dimensionless pressure residual `(P_calc - P) / P`.
    pub relative_residual: f64,
    /// Iterations executed.
    pub iterations: usize,
}

/// Converged dew-point solution `(T_dew, x)` at caller-specified pressure `P` and vapor `y`.
#[derive(Clone, Debug, PartialEq)]
pub struct DewPointV1 {
    /// Dew-point temperature `T_dew` [K].
    pub temperature: Temperature,
    /// Specified total pressure `P` [Pa].
    pub pressure: Pressure,
    /// Equilibrium liquid mole fractions `x`.
    pub liquid_x: Vec<f64>,
    /// Normalized vapor mole fractions `y`.
    pub vapor_y: Vec<f64>,
    /// Partition ratios `K_i = gamma_i P_i^sat / P`.
    pub k_values: Vec<f64>,
    /// Underlying excess-Gibbs evaluation at `(T_dew, x)`.
    pub excess_gibbs: ExcessGibbsEvaluationV1,
    /// Final dimensionless pressure residual `(P_dew - P) / P`.
    pub relative_residual: f64,
    /// Iterations executed.
    pub iterations: usize,
}

/// Phase regime resolved by an isothermal `(T, P, z)` Modified Raoult flash.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FlashPhaseRegimeV1 {
    /// Total pressure is at or above bubble pressure (`beta = 0`).
    SubcooledLiquid,
    /// Feed splits into coexisting liquid and vapor phases (`0 < beta < 1`).
    TwoPhase,
    /// Total pressure is at or below dew pressure (`beta = 1`).
    SuperheatedVapor,
}

/// Converged isothermal `(T, P, z)` flash result.
#[derive(Clone, Debug, PartialEq)]
pub struct IsothermalFlashV1 {
    /// Specified temperature [K].
    pub temperature: Temperature,
    /// Specified pressure [Pa].
    pub pressure: Pressure,
    /// Normalized overall feed mole fractions `z`.
    pub feed_z: Vec<f64>,
    /// Vapor molar fraction `beta` in `[0, 1]`.
    pub vapor_fraction: f64,
    /// Liquid phase mole fractions `x`.
    pub liquid_x: Vec<f64>,
    /// Vapor phase mole fractions `y`.
    pub vapor_y: Vec<f64>,
    /// Equilibrium K-values `K_i = gamma_i(T, x) P_i^sat(T) / P`.
    pub k_values: Vec<f64>,
    /// Identified phase regime.
    pub phase_regime: FlashPhaseRegimeV1,
    /// Maximum relative fugacity/equilibrium residual across components.
    pub relative_residual: f64,
    /// Outer flash iterations executed.
    pub iterations: usize,
}

/// Modified Raoult's law VLE evaluator (`y_i P = x_i gamma_i(T, x) P_i^sat(T)`)
/// parameterized by an explicit excess-Gibbs model `M` and per-species vapor-pressure laws.
///
/// Never hardcodes sea-level atmospheric pressure: every VLE query takes an explicit
/// caller-supplied `Pressure` (`P_atm(z)` or vessel pressure).
#[derive(Clone, Debug, PartialEq)]
pub struct ModifiedRaoultVleV1<M> {
    model: M,
    vapor_pressures: Vec<VaporPressureLawV1>,
    effective_t_min_k: f64,
    effective_t_max_k: f64,
}

impl<M: ExcessGibbsModelV1> ModifiedRaoultVleV1<M> {
    /// Bind an excess-Gibbs liquid model with one vapor-pressure law per species.
    pub fn new(
        model: M,
        vapor_pressures: Vec<VaporPressureLawV1>,
    ) -> Result<Self, ExcessGibbsErrorV1> {
        let n = model.species().len();
        if vapor_pressures.len() != n {
            return Err(ExcessGibbsErrorV1::InvalidModel(
                "vapor_pressures length must match excess-Gibbs species count",
            ));
        }
        let [t_min, t_max] = model.temperature_domain();
        let mut lo = t_min.value();
        let mut hi = t_max.value();
        for vp in &vapor_pressures {
            if let Some([v_lo, v_hi]) = vp.temperature_domain() {
                lo = lo.max(v_lo);
                hi = hi.min(v_hi);
            }
        }
        if !lo.is_finite() || !hi.is_finite() || lo <= 0.0 || hi < lo {
            return Err(ExcessGibbsErrorV1::InvalidModel(
                "excess-Gibbs and vapor-pressure temperature domains have empty intersection",
            ));
        }
        Ok(Self {
            model,
            vapor_pressures,
            effective_t_min_k: lo,
            effective_t_max_k: hi,
        })
    }

    /// Borrow the underlying excess-Gibbs model.
    #[must_use]
    pub const fn model(&self) -> &M {
        &self.model
    }

    /// Borrow the per-species saturation vapor-pressure laws.
    #[must_use]
    pub fn vapor_pressures(&self) -> &[VaporPressureLawV1] {
        &self.vapor_pressures
    }

    /// Effective inclusive temperature domain `[T_min, T_max]` [K] (intersection of all laws).
    #[must_use]
    pub const fn temperature_domain(&self) -> [Temperature; 2] {
        [
            Temperature::new(self.effective_t_min_k),
            Temperature::new(self.effective_t_max_k),
        ]
    }

    /// Evaluate Modified Raoult's law K-values `K_i = gamma_i(T, x) * P_i^sat(T) / P`.
    pub fn k_values(
        &self,
        temperature: Temperature,
        pressure: Pressure,
        liquid_x: &[f64],
    ) -> Result<VleKValuesV1, ExcessGibbsErrorV1> {
        let p = check_positive_pressure(pressure)?;
        check_temperature_in_domain(temperature, self.effective_t_min_k, self.effective_t_max_k)?;
        let eval = self.model.evaluate(temperature, liquid_x)?;
        let n = self.vapor_pressures.len();
        let mut k_values = Vec::with_capacity(n);
        let mut kx = Vec::with_capacity(n);
        let mut sum_kx = 0.0;
        for i in 0..n {
            let vp = self.vapor_pressures[i].evaluate(temperature)?;
            let ki = check_finite(eval.gamma[i] * vp.pressure_pa / p, "VLE K_i")?;
            let term = ki * eval.mole_fractions[i];
            k_values.push(ki);
            kx.push(term);
            sum_kx += term;
        }
        let sum_kx = check_finite(sum_kx, "VLE sum(K_i x_i)")?;
        let vapor_y_normalized = if sum_kx > 0.0 {
            kx.into_iter().map(|v| v / sum_kx).collect()
        } else {
            vec![0.0; n]
        };
        Ok(VleKValuesV1 {
            temperature,
            pressure,
            liquid_x: eval.mole_fractions.clone(),
            k_values,
            vapor_y_normalized,
            sum_kx,
            excess_gibbs: eval,
        })
    }

    /// Solve bubble-point temperature `T_bub(P, x)` and equilibrium vapor composition `y`
    /// satisfying `sum_i x_i * gamma_i(T, x) * P_i^sat(T) = P` via safeguarded Newton-bisection.
    pub fn bubble_point_t(
        &self,
        pressure: Pressure,
        liquid_x: &[f64],
        options: VleSolverOptionsV1,
    ) -> Result<BubblePointV1, ExcessGibbsErrorV1> {
        options.validate()?;
        let p = check_positive_pressure(pressure)?;
        let n = self.vapor_pressures.len();
        let x = validate_and_normalize_mole_fractions(liquid_x, n)?;
        if !(0..n).any(|i| x[i] > 0.0 && self.vapor_pressures[i].is_volatile()) {
            return Err(ExcessGibbsErrorV1::InvalidComposition(
                "bubble-point solve requires at least one present volatile species",
            ));
        }

        let eval_bubble = |t_k: f64| -> Result<
            (f64, f64, ExcessGibbsEvaluationV1, Vec<f64>),
            ExcessGibbsErrorV1,
        > {
            let temp = Temperature::new(t_k);
            let eg = self.model.evaluate(temp, &x)?;
            let mut p_calc = 0.0;
            let mut dp_dt = 0.0;
            let mut p_sat_vec = Vec::with_capacity(n);
            for i in 0..n {
                let vp = self.vapor_pressures[i].evaluate(temp)?;
                p_sat_vec.push(vp.pressure_pa);
                if x[i] > 0.0 && vp.pressure_pa > 0.0 {
                    let partial_p = x[i] * eg.gamma[i] * vp.pressure_pa;
                    p_calc += partial_p;
                    dp_dt += x[i]
                        * eg.gamma[i]
                        * (eg.d_ln_gamma_dt[i] * vp.pressure_pa + vp.dp_dt_pa_per_k);
                }
            }
            let p_calc = check_finite(p_calc, "bubble pressure")?;
            let dp_dt = check_finite(dp_dt, "bubble pressure derivative")?;
            Ok((p_calc, dp_dt, eg, p_sat_vec))
        };

        let mut lo = self.effective_t_min_k;
        let mut hi = self.effective_t_max_k;
        let (p_lo, _, eg_lo, psat_lo) = eval_bubble(lo)?;
        let (p_hi, _, eg_hi, psat_hi) = eval_bubble(hi)?;
        let f_lo = p_lo - p;
        let f_hi = p_hi - p;

        if f_lo.abs() / p <= options.relative_residual_tolerance {
            return Ok(build_bubble_result(
                Temperature::new(lo),
                pressure,
                &x,
                eg_lo,
                &psat_lo,
                p_lo,
                0,
            ));
        }
        if f_hi.abs() / p <= options.relative_residual_tolerance {
            return Ok(build_bubble_result(
                Temperature::new(hi),
                pressure,
                &x,
                eg_hi,
                &psat_hi,
                p_hi,
                0,
            ));
        }
        if f_lo > 0.0 || f_hi < 0.0 {
            let residual = if f_lo > 0.0 { f_lo / p } else { f_hi / p };
            return Err(ExcessGibbsErrorV1::ConvergenceFailed {
                iterations: 0,
                residual,
            });
        }

        // Clausius-Clapeyron log-linear initial guess in 1/T inside (lo, hi).
        let mut t = if p_lo > 0.0 && p_hi > p_lo {
            let frac = (ln(p) - ln(p_lo)) / (ln(p_hi) - ln(p_lo));
            let inv_t = (1.0 / lo) + frac * ((1.0 / hi) - (1.0 / lo));
            let guess = 1.0 / inv_t;
            if guess.is_finite() && guess > lo && guess < hi {
                guess
            } else {
                0.5 * (lo + hi)
            }
        } else {
            0.5 * (lo + hi)
        };

        let mut last_rel = f64::INFINITY;
        for iter in 1..=options.maximum_iterations {
            let (p_calc, dp_dt, eg, psat_vec) = eval_bubble(t)?;
            let f = p_calc - p;
            last_rel = f / p;
            if last_rel.abs() <= options.relative_residual_tolerance
                || (hi - lo) <= options.temperature_tolerance_k
            {
                return Ok(build_bubble_result(
                    Temperature::new(t),
                    pressure,
                    &x,
                    eg,
                    &psat_vec,
                    p_calc,
                    iter,
                ));
            }
            if f < 0.0 {
                lo = t;
            } else {
                hi = t;
            }
            let newton = if dp_dt > 0.0 && dp_dt.is_finite() {
                t - f / dp_dt
            } else {
                f64::NAN
            };
            let next = if newton.is_finite() && newton > lo && newton < hi && newton != t {
                newton
            } else {
                0.5 * (lo + hi)
            };
            if next <= lo || next >= hi {
                return Ok(build_bubble_result(
                    Temperature::new(t),
                    pressure,
                    &x,
                    eg,
                    &psat_vec,
                    p_calc,
                    iter,
                ));
            }
            t = next;
        }

        Err(ExcessGibbsErrorV1::ConvergenceFailed {
            iterations: options.maximum_iterations,
            residual: last_rel,
        })
    }

    /// Solve dew-point temperature `T_dew(P, y)` and incipient equilibrium liquid composition `x`
    /// satisfying `y_i P = x_i gamma_i(T, x) P_i^sat(T)` and `sum_i x_i = 1`.
    pub fn dew_point_t(
        &self,
        pressure: Pressure,
        vapor_y: &[f64],
        options: VleSolverOptionsV1,
    ) -> Result<DewPointV1, ExcessGibbsErrorV1> {
        options.validate()?;
        let p = check_positive_pressure(pressure)?;
        let n = self.vapor_pressures.len();
        let y = validate_and_normalize_mole_fractions(vapor_y, n)?;
        for i in 0..n {
            if y[i] > 0.0 && !self.vapor_pressures[i].is_volatile() {
                return Err(ExcessGibbsErrorV1::InvalidComposition(
                    "non-volatile species cannot have positive vapor mole fraction in dew-point solve",
                ));
            }
        }

        let mut lo = self.effective_t_min_k;
        let mut hi = self.effective_t_max_k;
        let (p_lo, eg_lo, psat_lo) = self.dew_pressure_at_t(Temperature::new(lo), &y)?;
        let (p_hi, eg_hi, psat_hi) = self.dew_pressure_at_t(Temperature::new(hi), &y)?;
        let f_lo = p_lo - p;
        let f_hi = p_hi - p;

        if f_lo.abs() / p <= options.relative_residual_tolerance {
            return Ok(build_dew_result(
                Temperature::new(lo),
                pressure,
                &y,
                eg_lo,
                &psat_lo,
                p_lo,
                0,
            ));
        }
        if f_hi.abs() / p <= options.relative_residual_tolerance {
            return Ok(build_dew_result(
                Temperature::new(hi),
                pressure,
                &y,
                eg_hi,
                &psat_hi,
                p_hi,
                0,
            ));
        }
        if f_lo > 0.0 || f_hi < 0.0 {
            let residual = if f_lo > 0.0 { f_lo / p } else { f_hi / p };
            return Err(ExcessGibbsErrorV1::ConvergenceFailed {
                iterations: 0,
                residual,
            });
        }

        // Safeguarded secant-bisection in log-pressure coordinates (nearly linear in 1/T).
        let ln_p = ln(p);
        let mut g_lo = ln(p_lo) - ln_p;
        let mut g_hi = ln(p_hi) - ln_p;
        let mut last_rel = f64::INFINITY;
        let mut use_secant = true;
        for iter in 1..=options.maximum_iterations {
            let t = if use_secant
                && g_lo.is_finite()
                && g_hi.is_finite()
                && (g_hi - g_lo).abs() > 0.0
            {
                let sec = lo - g_lo * (hi - lo) / (g_hi - g_lo);
                let margin = 0.05 * (hi - lo);
                if sec.is_finite() && sec > lo + margin && sec < hi - margin {
                    sec
                } else {
                    0.5 * (lo + hi)
                }
            } else {
                0.5 * (lo + hi)
            };
            let (p_calc, eg, psat_vec) = self.dew_pressure_at_t(Temperature::new(t), &y)?;
            let f = p_calc - p;
            last_rel = f / p;
            if last_rel.abs() <= options.relative_residual_tolerance
                || (hi - lo) <= options.temperature_tolerance_k
            {
                return Ok(build_dew_result(
                    Temperature::new(t),
                    pressure,
                    &y,
                    eg,
                    &psat_vec,
                    p_calc,
                    iter,
                ));
            }
            let prev_span = hi - lo;
            let g_mid = ln(p_calc) - ln_p;
            if f < 0.0 {
                lo = t;
                g_lo = g_mid;
            } else {
                hi = t;
                g_hi = g_mid;
            }
            // Alternate to bisection whenever secant shrinks bracket slowly.
            use_secant = (hi - lo) < 0.75 * prev_span;
            if hi <= lo {
                return Ok(build_dew_result(
                    Temperature::new(t),
                    pressure,
                    &y,
                    eg,
                    &psat_vec,
                    p_calc,
                    iter,
                ));
            }
        }

        Err(ExcessGibbsErrorV1::ConvergenceFailed {
            iterations: options.maximum_iterations,
            residual: last_rel,
        })
    }

    /// Compute the mixture boiling range `(bubble, dew, T_dew - T_bub)` at `pressure`
    /// and overall composition `composition` (used by Schluender-Thome mixture boiling
    /// suppression in wine deglazing and sauce reduction).
    pub fn boiling_range_k(
        &self,
        pressure: Pressure,
        composition: &[f64],
        options: VleSolverOptionsV1,
    ) -> Result<(BubblePointV1, DewPointV1, f64), ExcessGibbsErrorV1> {
        let bubble = self.bubble_point_t(pressure, composition, options)?;
        let dew = self.dew_point_t(pressure, composition, options)?;
        let delta_t_bp = dew.temperature.value() - bubble.temperature.value();
        Ok((bubble, dew, delta_t_bp))
    }

    /// Isothermal `(T, P, z)` flash using Modified Raoult's law and Rachford-Rice.
    pub fn flash_pt(
        &self,
        temperature: Temperature,
        pressure: Pressure,
        overall_z: &[f64],
        options: VleSolverOptionsV1,
    ) -> Result<IsothermalFlashV1, ExcessGibbsErrorV1> {
        options.validate()?;
        let p = check_positive_pressure(pressure)?;
        check_temperature_in_domain(temperature, self.effective_t_min_k, self.effective_t_max_k)?;
        let n = self.vapor_pressures.len();
        let z = validate_and_normalize_mole_fractions(overall_z, n)?;

        let mut p_sat = Vec::with_capacity(n);
        for vp in &self.vapor_pressures {
            p_sat.push(vp.evaluate(temperature)?.pressure_pa);
        }

        // 1. Bubble-point pressure test at (T, z).
        let eg_z = self.model.evaluate(temperature, &z)?;
        let mut p_bub = 0.0;
        let mut k_bub = Vec::with_capacity(n);
        for i in 0..n {
            let ki = check_finite(eg_z.gamma[i] * p_sat[i] / p, "flash K_bub")?;
            k_bub.push(ki);
            p_bub += z[i] * eg_z.gamma[i] * p_sat[i];
        }
        let p_bub = check_finite(p_bub, "flash bubble pressure")?;
        if p >= p_bub * (1.0 - options.relative_residual_tolerance) {
            let y_inc = if p_bub > 0.0 {
                (0..n)
                    .map(|i| z[i] * eg_z.gamma[i] * p_sat[i] / p_bub)
                    .collect()
            } else {
                vec![0.0; n]
            };
            return Ok(IsothermalFlashV1 {
                temperature,
                pressure,
                feed_z: z.clone(),
                vapor_fraction: 0.0,
                liquid_x: z,
                vapor_y: y_inc,
                k_values: k_bub,
                phase_regime: FlashPhaseRegimeV1::SubcooledLiquid,
                relative_residual: 0.0,
                iterations: 0,
            });
        }

        // 2. Dew-point pressure test (only possible if every present species is volatile).
        let all_present_volatile =
            (0..n).all(|i| z[i] == 0.0 || self.vapor_pressures[i].is_volatile());
        let mut x_init = z.clone();
        if all_present_volatile {
            let (p_dew, eg_dew, _) = self.dew_pressure_at_t(temperature, &z)?;
            if p <= p_dew * (1.0 + options.relative_residual_tolerance) {
                let k_dew = (0..n).map(|i| eg_dew.gamma[i] * p_sat[i] / p).collect();
                return Ok(IsothermalFlashV1 {
                    temperature,
                    pressure,
                    feed_z: z.clone(),
                    vapor_fraction: 1.0,
                    liquid_x: eg_dew.mole_fractions,
                    vapor_y: z,
                    k_values: k_dew,
                    phase_regime: FlashPhaseRegimeV1::SuperheatedVapor,
                    relative_residual: 0.0,
                    iterations: 0,
                });
            }
            if p_bub > p_dew {
                let theta = ((p_bub - p) / (p_bub - p_dew)).clamp(0.05, 0.95);
                for i in 0..n {
                    x_init[i] = (1.0 - theta) * z[i] + theta * eg_dew.mole_fractions[i];
                }
                normalize_positive_in_place(&mut x_init)?;
            }
        }

        // 3. Two-phase Rachford-Rice + bounded Wegstein-accelerated ln(K) substitution.
        let mut x = x_init;
        let eg0 = self.model.evaluate(temperature, &x)?;
        let mut k_vals: Vec<f64> = (0..n).map(|i| eg0.gamma[i] * p_sat[i] / p).collect();
        let mut y = vec![0.0; n];
        let mut prev_u: Option<Vec<f64>> = None;
        let mut prev_gu: Option<Vec<f64>> = None;
        let mut last_res = f64::INFINITY;

        for iter in 1..=options.maximum_iterations {
            let beta = solve_rachford_rice(&z, &k_vals)?;

            for i in 0..n {
                // (1 - beta) + beta * K_i avoids catastrophic 1 - 1 cancellation for tiny K_i.
                let denom = (1.0 - beta) + beta * k_vals[i];
                if !denom.is_finite() || denom <= 0.0 {
                    return Err(ExcessGibbsErrorV1::Unrepresentable(
                        "Rachford-Rice denominator became non-positive",
                    ));
                }
                x[i] = (z[i] / denom).max(0.0);
                y[i] = (k_vals[i] * x[i]).max(0.0);
            }
            normalize_positive_in_place(&mut x)?;
            normalize_positive_in_place(&mut y)?;

            let eg = self.model.evaluate(temperature, &x)?;
            let mut max_eq_res = 0.0_f64;
            let mut max_mb_res = 0.0_f64;
            let mut k_new = Vec::with_capacity(n);
            for i in 0..n {
                let ki_new = check_finite(eg.gamma[i] * p_sat[i] / p, "flash K_new")?;
                k_new.push(ki_new);
                let eq_res = (y[i] - ki_new * x[i]).abs();
                let mb_res = ((1.0 - beta) * x[i] + beta * y[i] - z[i]).abs();
                max_eq_res = max_eq_res.max(eq_res);
                max_mb_res = max_mb_res.max(mb_res);
            }
            let res = max_eq_res.max(max_mb_res);
            if res <= options.relative_residual_tolerance {
                return Ok(IsothermalFlashV1 {
                    temperature,
                    pressure,
                    feed_z: z,
                    vapor_fraction: beta,
                    liquid_x: x,
                    vapor_y: y,
                    k_values: k_new,
                    phase_regime: FlashPhaseRegimeV1::TwoPhase,
                    relative_residual: res,
                    iterations: iter,
                });
            }
            let u: Vec<f64> = k_vals
                .iter()
                .map(|&ki| if ki > 0.0 { ln(ki) } else { 0.0 })
                .collect();
            let gu: Vec<f64> = k_new
                .iter()
                .map(|&ki| if ki > 0.0 { ln(ki) } else { 0.0 })
                .collect();
            let mut k_next = vec![0.0; n];
            if res > last_res {
                for i in 0..n {
                    k_next[i] = if k_new[i] > 0.0 && k_vals[i] > 0.0 {
                        exp(0.5 * (u[i] + gu[i]))
                    } else {
                        k_new[i]
                    };
                }
            } else if let (Some(pu), Some(pgu)) = (&prev_u, &prev_gu) {
                for i in 0..n {
                    if k_new[i] <= 0.0 || k_vals[i] <= 0.0 {
                        k_next[i] = k_new[i];
                        continue;
                    }
                    let du = u[i] - pu[i];
                    let dgu = gu[i] - pgu[i];
                    if du.abs() > 1e-15 {
                        let slope = dgu / du;
                        let denom = 1.0 - slope;
                        let w = if denom.abs() > 1e-12 {
                            (1.0 / denom).clamp(0.1, 1.5)
                        } else {
                            0.5
                        };
                        let u_next = check_finite(u[i] + w * (gu[i] - u[i]), "Wegstein ln(K)")?;
                        k_next[i] = check_finite(exp(u_next), "Wegstein K")?;
                    } else {
                        k_next[i] = k_new[i];
                    }
                }
            } else {
                for i in 0..n {
                    k_next[i] = if k_new[i] > 0.0 && k_vals[i] > 0.0 {
                        exp(0.5 * (u[i] + gu[i]))
                    } else {
                        k_new[i]
                    };
                }
            }
            prev_u = Some(u);
            prev_gu = Some(gu);
            k_vals = k_next;
            last_res = res;
        }

        Err(ExcessGibbsErrorV1::ConvergenceFailed {
            iterations: options.maximum_iterations,
            residual: last_res,
        })
    }

    fn dew_pressure_at_t(
        &self,
        temperature: Temperature,
        y: &[f64],
    ) -> Result<(f64, ExcessGibbsEvaluationV1, Vec<f64>), ExcessGibbsErrorV1> {
        let n = self.vapor_pressures.len();
        let mut p_sat = Vec::with_capacity(n);
        let mut x = vec![0.0; n];
        let mut w_sum = 0.0;
        for i in 0..n {
            let ps = self.vapor_pressures[i].evaluate(temperature)?.pressure_pa;
            p_sat.push(ps);
            if y[i] > 0.0 {
                if ps <= 0.0 {
                    return Err(ExcessGibbsErrorV1::Unrepresentable(
                        "volatile species has non-positive saturation pressure",
                    ));
                }
                let wi = y[i] / ps;
                x[i] = wi;
                w_sum += wi;
            }
        }
        if !w_sum.is_finite() || w_sum <= 0.0 {
            return Err(ExcessGibbsErrorV1::Unrepresentable(
                "dew-point ideal weight sum is non-positive",
            ));
        }
        for xi in &mut x {
            *xi /= w_sum;
        }

        let mut prev_diff = f64::INFINITY;
        for _ in 0..128 {
            let eg = self.model.evaluate(temperature, &x)?;
            let mut next_x = vec![0.0; n];
            let mut w = 0.0;
            for i in 0..n {
                if y[i] > 0.0 {
                    let denom = eg.gamma[i] * p_sat[i];
                    if !denom.is_finite() || denom <= 0.0 {
                        return Err(ExcessGibbsErrorV1::Unrepresentable(
                            "dew-point denominator gamma_i * P_i^sat is non-positive",
                        ));
                    }
                    let wi = y[i] / denom;
                    next_x[i] = wi;
                    w += wi;
                }
            }
            if !w.is_finite() || w <= 0.0 {
                return Err(ExcessGibbsErrorV1::Unrepresentable(
                    "dew-point liquid normalization sum is non-positive",
                ));
            }
            let mut max_diff = 0.0_f64;
            for i in 0..n {
                next_x[i] /= w;
                max_diff = max_diff.max((next_x[i] - x[i]).abs());
            }
            if max_diff <= 1e-13 {
                let eg_final = self.model.evaluate(temperature, &next_x)?;
                let mut w_final = 0.0;
                for i in 0..n {
                    if y[i] > 0.0 {
                        w_final += y[i] / (eg_final.gamma[i] * p_sat[i]);
                    }
                }
                let p_dew = check_finite(1.0 / w_final, "dew pressure")?;
                return Ok((p_dew, eg_final, p_sat));
            }
            if max_diff > prev_diff {
                for i in 0..n {
                    x[i] = 0.5 * x[i] + 0.5 * next_x[i];
                }
                normalize_positive_in_place(&mut x)?;
            } else {
                x = next_x;
            }
            prev_diff = max_diff;
        }

        Err(ExcessGibbsErrorV1::ConvergenceFailed {
            iterations: 128,
            residual: prev_diff,
        })
    }
}

fn solve_rachford_rice(z: &[f64], k: &[f64]) -> Result<f64, ExcessGibbsErrorV1> {
    let eval_g = |beta: f64| -> (f64, f64) {
        let mut g = 0.0;
        let mut gp = 0.0;
        for (&zi, &ki) in z.iter().zip(k) {
            if zi == 0.0 {
                continue;
            }
            let km1 = ki - 1.0;
            // (1 - beta) + beta * K_i is a sum of nonnegative terms on [0, 1], avoiding
            // cancellation when 0 < K_i < 2^-53 near beta = 1.
            let denom = (1.0 - beta) + beta * ki;
            if denom <= 0.0 {
                return (f64::NEG_INFINITY, f64::NEG_INFINITY);
            }
            let ratio = km1 / denom;
            g += zi * ratio;
            gp -= zi * ratio * ratio;
        }
        (g, gp)
    };

    let (g0, _) = eval_g(0.0);
    if !g0.is_finite() {
        return Err(ExcessGibbsErrorV1::Unrepresentable(
            "Rachford-Rice g(0) is non-finite",
        ));
    }
    if g0 <= 0.0 {
        return Ok(1e-9);
    }
    let has_non_volatile = z.iter().zip(k).any(|(&zi, &ki)| zi > 0.0 && ki == 0.0);
    if !has_non_volatile {
        let (g1, _) = eval_g(1.0);
        if g1 >= 0.0 {
            return Ok(1.0 - 1e-9);
        }
    }

    let mut lo = 0.0;
    let mut hi = 1.0;
    let mut beta = 0.5;
    for _ in 0..100 {
        let (g, gp) = eval_g(beta);
        if !g.is_finite() {
            beta = 0.5 * (lo + hi);
            continue;
        }
        if g.abs() <= 1e-14 || (hi - lo) <= 1e-14 {
            return Ok(beta);
        }
        if g > 0.0 {
            lo = beta;
        } else {
            hi = beta;
        }
        let newton = if gp < 0.0 && gp.is_finite() {
            beta - g / gp
        } else {
            f64::NAN
        };
        beta = if newton.is_finite() && newton > lo && newton < hi {
            newton
        } else {
            0.5 * (lo + hi)
        };
    }
    Ok(0.5 * (lo + hi))
}

fn build_bubble_result(
    temperature: Temperature,
    pressure: Pressure,
    x: &[f64],
    eg: ExcessGibbsEvaluationV1,
    psat_vec: &[f64],
    p_calc: f64,
    iterations: usize,
) -> BubblePointV1 {
    let p = pressure.value();
    let n = x.len();
    let mut k_values = Vec::with_capacity(n);
    let mut vapor_y = Vec::with_capacity(n);
    for i in 0..n {
        let ki = eg.gamma[i] * psat_vec[i] / p;
        k_values.push(ki);
        vapor_y.push(if p_calc > 0.0 {
            x[i] * eg.gamma[i] * psat_vec[i] / p_calc
        } else {
            0.0
        });
    }
    BubblePointV1 {
        temperature,
        pressure,
        liquid_x: x.to_vec(),
        vapor_y,
        k_values,
        excess_gibbs: eg,
        relative_residual: (p_calc - p) / p,
        iterations,
    }
}

fn build_dew_result(
    temperature: Temperature,
    pressure: Pressure,
    y: &[f64],
    eg: ExcessGibbsEvaluationV1,
    psat_vec: &[f64],
    p_dew: f64,
    iterations: usize,
) -> DewPointV1 {
    let p = pressure.value();
    let k_values = (0..y.len())
        .map(|i| eg.gamma[i] * psat_vec[i] / p)
        .collect();
    DewPointV1 {
        temperature,
        pressure,
        liquid_x: eg.mole_fractions.clone(),
        vapor_y: y.to_vec(),
        k_values,
        excess_gibbs: eg,
        relative_residual: (p_dew - p) / p,
        iterations,
    }
}

fn normalize_positive_in_place(values: &mut [f64]) -> Result<(), ExcessGibbsErrorV1> {
    let mut sum = 0.0;
    let mut comp = 0.0;
    for &v in values.iter() {
        neumaier_add(&mut sum, &mut comp, v);
    }
    let total = sum + comp;
    if !total.is_finite() || total <= 0.0 {
        return Err(ExcessGibbsErrorV1::Unrepresentable(
            "phase mole-fraction normalization sum is non-positive",
        ));
    }
    for v in values.iter_mut() {
        *v /= total;
    }
    Ok(())
}

pub(crate) fn check_finite(value: f64, label: &'static str) -> Result<f64, ExcessGibbsErrorV1> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(ExcessGibbsErrorV1::Unrepresentable(label))
    }
}

pub(crate) fn check_positive_pressure(pressure: Pressure) -> Result<f64, ExcessGibbsErrorV1> {
    let p = pressure.value();
    if !p.is_finite() || p <= 0.0 {
        return Err(ExcessGibbsErrorV1::InvalidPressure { pressure_pa: p });
    }
    Ok(p)
}

pub(crate) fn validate_distinct_species(species: &[SpeciesId]) -> Result<(), ExcessGibbsErrorV1> {
    for i in 0..species.len() {
        for j in (i + 1)..species.len() {
            if species[i] == species[j] {
                return Err(ExcessGibbsErrorV1::InvalidModel(
                    "duplicate species identity in excess-Gibbs model",
                ));
            }
        }
    }
    Ok(())
}

pub(crate) fn validate_temperature_domain(
    min_t: Temperature,
    max_t: Temperature,
) -> Result<(), ExcessGibbsErrorV1> {
    let lo = min_t.value();
    let hi = max_t.value();
    if !lo.is_finite() || lo <= 0.0 || !hi.is_finite() || hi < lo {
        return Err(ExcessGibbsErrorV1::InvalidModel(
            "temperature domain must be finite, strictly positive, and ordered",
        ));
    }
    Ok(())
}

pub(crate) fn check_temperature_in_domain(
    temperature: Temperature,
    min_k: f64,
    max_k: f64,
) -> Result<f64, ExcessGibbsErrorV1> {
    let t = temperature.value();
    if !t.is_finite() || t <= 0.0 || t < min_k || t > max_k {
        return Err(ExcessGibbsErrorV1::TemperatureOutOfDomain {
            temperature_k: t,
            min_k,
            max_k,
        });
    }
    Ok(t)
}

pub(crate) fn validate_and_normalize_mole_fractions(
    mole_fractions: &[f64],
    expected_len: usize,
) -> Result<Vec<f64>, ExcessGibbsErrorV1> {
    if mole_fractions.len() != expected_len {
        return Err(ExcessGibbsErrorV1::InvalidComposition(
            "mole fraction length does not match species count",
        ));
    }
    let mut sum = 0.0;
    let mut comp = 0.0;
    for &x in mole_fractions {
        if !x.is_finite() || !(0.0..=1.0).contains(&x) {
            return Err(ExcessGibbsErrorV1::InvalidComposition(
                "mole fraction must be finite and in [0, 1]",
            ));
        }
        neumaier_add(&mut sum, &mut comp, x);
    }
    let total = sum + comp;
    if !total.is_finite() || (total - 1.0).abs() > EXCESS_GIBBS_COMPOSITION_SUM_TOLERANCE_V1 {
        return Err(ExcessGibbsErrorV1::InvalidComposition(
            "mole fractions must sum to unity",
        ));
    }
    let mut normalized = Vec::with_capacity(expected_len);
    for &x in mole_fractions {
        let nx = x / total;
        if x > 0.0 && nx == 0.0 {
            return Err(ExcessGibbsErrorV1::Unrepresentable(
                "positive trace mole fraction underflowed during normalization",
            ));
        }
        normalized.push(nx);
    }
    Ok(normalized)
}

fn neumaier_add(sum: &mut f64, comp: &mut f64, term: f64) {
    let next = *sum + term;
    *comp += if sum.abs() >= term.abs() {
        (*sum - next) + term
    } else {
        (term - next) + *sum
    };
    *sum = next;
}
