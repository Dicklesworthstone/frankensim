//! Symmetric binary liquid miscibility from one molar mixing free energy.
//!
//! `g_mix/(RT) = x ln x + (1-x) ln(1-x) + chi(T) x(1-x)`, with
//! `chi(T) = A + B/T`. This is an equal-size, non-electrolyte regular-solution
//! model with an optional excess-entropy term, not a prediction from densities
//! or material labels. The caller owns the species-pair fit and its validity.
//! Pure-component reference chemical potentials cancel between liquid phases.
//! Pressure dependence, unequal molecular sizes, interfacial energy, nucleation
//! and spatial evolution are outside this constitutive model.
//!
//! The molar form follows the regular-solution free-energy implementation in
//! INL MOOSE (its `omega*c*(1-c)` interaction term):
//! <https://mooseframework.inl.gov/docs/doxygen/modules/classRegularSolutionFreeEnergy.html>.
//! The binodal is the common tangent; for this symmetric model it is horizontal.
//! A negative thermodynamic factor is a spinodal instability, never a negative
//! molecular diffusivity to feed an ordinary Fick evolution step.

use crate::{SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};
use fs_math::det::{exp, expm1, ln, sqrt};
use fs_qty::Temperature;

/// Finite solid/liquid partition for an explicitly declared neutral solute.
pub mod solubility;

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
            minimum_temperature,
            maximum_temperature,
        })
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
        let t = temperature.value();
        if !t.is_finite()
            || t < self.minimum_temperature.value()
            || t > self.maximum_temperature.value()
        {
            return Err(BinarySolutionErrorV1::OutsideTemperatureDomain);
        }
        finite(
            self.chi_constant + self.chi_inverse_temperature_k / t,
            "chi(T)",
        )
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
                (self.chi_inverse_temperature_k * xy) * R,
                "mixing enthalpy",
            )?,
            mixing_entropy_j_per_mol_k: finite(
                -R * (ideal + self.chi_constant * xy),
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
