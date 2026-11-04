//! The original homogeneous liquid law for conservative diffuse interfaces.
//!
//! The free energy is `F = integral[c g_mix(x,T) + kappa |grad x|^2/2] dV`.
//! Here c is a common constant site concentration [mol/m^3], so both components
//! have the same fixed partial molar volume. Kappa [J/m] is independently
//! qualified interfacial data, not a surface tension [J/m^2]. This L1 owner
//! supplies local constitutive terms; an L3 spatial owner supplies gradients,
//! conservative fluxes, boundary conditions and time integration.
//!
//! Local states remain HOMOGENEOUS, including unstable compositions. Replacing
//! them with their bulk equilibrium convex envelope would erase spinodal
//! dynamics and double-count phase separation already represented by the field.

use super::{
    BinaryRegularSolutionV1, BinarySolutionErrorV1, BinarySolutionStateV1, finite, ln,
    ln_complement, normalized_fractions,
};
use crate::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;
use fs_qty::{AmountConcentration, Temperature};

/// Immutable isothermal constitutive owner for a binary square-gradient field.
#[derive(Clone, Debug, PartialEq)]
pub struct BinarySquareGradientV1 {
    thermodynamics: BinaryRegularSolutionV1,
    temperature: Temperature,
    molar_concentration_mol_m3: f64,
    gradient_coefficient_j_m: f64,
    chi: f64,
    rt_j_per_mol: f64,
}

impl BinarySquareGradientV1 {
    /// Retain the exact bulk pair, constant concentration and independent kappa.
    /// A nonconvex spatial continuum requires positive kappa. Zero kappa is
    /// admitted only in the convex regime, including the ordinary Fick limit.
    pub fn new(
        thermodynamics: BinaryRegularSolutionV1,
        temperature: Temperature,
        molar_concentration: AmountConcentration,
        gradient_coefficient_j_m: f64,
    ) -> Result<Self, BinarySolutionErrorV1> {
        let chi = thermodynamics.chi(temperature)?;
        let c = molar_concentration.value();
        if !c.is_finite()
            || c <= 0.0
            || !gradient_coefficient_j_m.is_finite()
            || gradient_coefficient_j_m < 0.0
            || (chi > 2.0 && gradient_coefficient_j_m == 0.0)
        {
            return Err(BinarySolutionErrorV1::Invalid(
                "positive site concentration and nonnegative kappa; nonconvex fields require kappa>0",
            ));
        }
        let rt_j_per_mol = finite(R * temperature.value(), "square-gradient RT")?;
        Ok(Self {
            thermodynamics,
            temperature,
            molar_concentration_mol_m3: c,
            gradient_coefficient_j_m,
            chi,
            rt_j_per_mol,
        })
    }

    pub const fn thermodynamics(&self) -> &BinaryRegularSolutionV1 {
        &self.thermodynamics
    }
    pub const fn temperature(&self) -> Temperature {
        self.temperature
    }
    pub const fn molar_concentration_mol_m3(&self) -> f64 {
        self.molar_concentration_mol_m3
    }
    pub const fn gradient_coefficient_j_m(&self) -> f64 {
        self.gradient_coefficient_j_m
    }
    pub const fn chi(&self) -> f64 {
        self.chi
    }
    pub const fn rt_j_per_mol(&self) -> f64 {
        self.rt_j_per_mol
    }

    /// Evaluate the SAME original homogeneous branch, never a local phase split.
    pub fn evaluate_composition(
        &self,
        fractions: [f64; 2],
    ) -> Result<BinarySolutionStateV1, BinarySolutionErrorV1> {
        self.thermodynamics
            .evaluate_composition(self.temperature, fractions)
    }

    /// Bulk exchange chemical potential mu_A-mu_B [J/mol]. Pure endpoints have
    /// divergent potentials and refuse; their energy remains well-defined.
    pub fn exchange_potential_j_per_mol(
        &self,
        fractions: [f64; 2],
    ) -> Result<f64, BinarySolutionErrorV1> {
        let fractions = interior(fractions)?;
        finite(
            self.rt_j_per_mol * reduced_exchange(fractions, self.chi),
            "binary exchange potential",
        )
    }

    /// Convex-split chemical potential: ideal entropy and negative chi are
    /// implicit; positive-chi concavity is evaluated at the original state.
    /// Old compositions may be pure, provided the new composition is interior.
    pub fn convex_split_exchange_potential_j_per_mol(
        &self,
        new: [f64; 2],
        old: [f64; 2],
    ) -> Result<f64, BinarySolutionErrorV1> {
        let new = interior(new)?;
        let old = normalized_fractions(old)?;
        let convex = reduced_exchange(new, self.chi.min(0.0));
        finite(
            self.rt_j_per_mol * (convex + self.chi.max(0.0) * (old[1] - old[0])),
            "convex-split exchange potential",
        )
    }

    /// Derivative of the implicit convex exchange term with respect to x_A.
    /// Strictly positive for interior fractions; units J/mol since x is unitless.
    pub fn convex_curvature_j_per_mol(
        &self,
        fractions: [f64; 2],
    ) -> Result<f64, BinarySolutionErrorV1> {
        let [x, y] = interior(fractions)?;
        finite(
            self.rt_j_per_mol * (1.0 / x + 1.0 / y - 2.0 * self.chi.min(0.0)),
            "convex mixing curvature",
        )
    }

    /// Direct homogeneous mixing-free-energy increment [J/mol]. Local
    /// increments preserve small changes without subtracting large molar
    /// energies. Reversing endpoints returns the exact opposite operation.
    pub fn mixing_free_energy_change_j_per_mol(
        &self,
        old: [f64; 2],
        new: [f64; 2],
    ) -> Result<f64, BinarySolutionErrorV1> {
        let old = normalized_fractions(old)?;
        let new = normalized_fractions(new)?;
        if old == new {
            return Ok(0.0);
        }
        let (a, b, sign) = if old[0] < new[0] || (old[0] == new[0] && old[1] < new[1]) {
            (old, new, 1.0)
        } else {
            (new, old, -1.0)
        };
        // Use the trace coordinate when both states have the same dominant
        // component. Its real increment survives even if the majority rounds 1.
        let (a, b) = if a[0] >= 0.5 && b[0] >= 0.5 {
            ([a[1], a[0]], [b[1], b[0]])
        } else {
            (a, b)
        };
        let d = b[0] - a[0];
        let change =
            if a[0] > 0.0 && a[1] > 0.0 && (d / a[0]).abs() < 0.125 && (d / a[1]).abs() < 0.125 {
                // Taylor-integrated ideal entropy, with the quadratic interaction
                // combined BEFORE scaling. This retains the fourth-order energy at
                // the symmetric critical point instead of cancelling two O(d^2)s.
                let [x, y] = a;
                let r = [d / x, -d / y];
                let mut powers = r;
                let mut tail = 0.0;
                for k in 3..=20 {
                    let coefficient = if k % 2 == 0 { 1.0 } else { -1.0 } / ((k * (k - 1)) as f64);
                    tail += coefficient * (powers[0] / x + powers[1] / y);
                    powers[0] *= r[0];
                    powers[1] *= r[1];
                }
                let gamma = (x - y) * (x - y) - 2.0 * ((self.chi - 2.0) * x * y);
                d * reduced_exchange(a, self.chi) + d * d * (gamma / (2.0 * x * y) + tail)
            } else {
                // Well-separated states do not need a cancellation-sensitive local
                // expansion; original pure-endpoint 0*ln(0) limits are retained.
                let ideal = ideal_mixing(b) - ideal_mixing(a);
                let product = d * (0.5 * (a[1] - a[0]) + 0.5 * (b[1] - b[0]));
                ideal + self.chi * product
            };
        finite(
            sign * self.rt_j_per_mol * change,
            "mixing free-energy increment",
        )
    }
}

fn interior(fractions: [f64; 2]) -> Result<[f64; 2], BinarySolutionErrorV1> {
    let fractions = normalized_fractions(fractions)?;
    if fractions.iter().any(|&x| x <= 0.0) {
        Err(BinarySolutionErrorV1::Invalid(
            "finite exchange potential requires both species present",
        ))
    } else {
        Ok(fractions)
    }
}

fn log_fraction(pair: [f64; 2], i: usize) -> f64 {
    if pair[i] > 0.5 {
        ln_complement(pair[1 - i])
    } else {
        ln(pair[i])
    }
}
fn ideal_mixing(pair: [f64; 2]) -> f64 {
    (0..2)
        .map(|i| {
            if pair[i] == 0.0 {
                0.0
            } else {
                pair[i] * log_fraction(pair, i)
            }
        })
        .sum()
}
fn reduced_exchange(pair: [f64; 2], chi: f64) -> f64 {
    let z = pair[0] - pair[1];
    if z.abs() < 0.125 {
        let z2 = z * z;
        let mut power = z2;
        let mut remainder = 0.0;
        for k in 1..=12 {
            remainder += power / ((2 * k + 1) as f64);
            power *= z2;
        }
        z * ((2.0 - chi) + 2.0 * remainder)
    } else {
        log_fraction(pair, 0) - log_fraction(pair, 1) - chi * z
    }
}

#[cfg(test)]
mod tests;
