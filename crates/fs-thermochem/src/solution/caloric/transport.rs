//! Source-conjugate binary liquid heat and molecular transport.
//!
//! The reference frame is molar: `j_A = -j_B`, with a common fixed site
//! concentration and equal fixed partial molar volumes. Pressure is fixed.
//! The energy flux includes the original partial-enthalpy exchange ONCE;
//! the reduced heat flux contains Fourier heat and reciprocal Dufour heat.
//! The authored `Q` is an A-minus-B molar heat-of-transport contrast, not
//! a partial enthalpy and not inferred from the excess mixing enthalpy.
//!
//! For beta=1/T, the two-point source identity is
//! `Delta[(mu_A-mu_B)/T] = a_chem + h_star Delta beta`, where
//! `a_chem = R [Delta ln(x_A/x_B) + mean(chi) Delta(x_B-x_A)]`.
//! The divided original Cp and interaction sources supply `h_star`. Thus a
//! temperature difference at uniform composition causes no molecular flux
//! when Q=0, even for different pure Cp or temperature-dependent interactions.
//! With L=c D x_A_face x_B_face/R and a=a_chem-Q Delta beta,
//! `j_A=-L a/d`, `q_reduced=-k Delta T/d+Q j_A`,
//! and `energy=h_star j_A+q_reduced`. The paired entropy rate per area is
//! `L a^2/d + k Delta T^2/(T_left T_right d) >= 0`.
//! This is the Onsager Schur-complement construction; see Van-Brunt et al.,
//! <https://arxiv.org/abs/2109.05082>. The two-point caloric identity is derived
//! here for this original regular-solution/Cp owner. Finite-step entropy and
//! conservation remain the spatial integrator's responsibility.
//!
//! Both species must be present at both endpoints; logarithmic floors are
//! forbidden. Spatial use without square-gradient regularization must admit
//! a convex temperature interval. No pressure gradient, thermal expansion,
//! nonlocal interface-energy localization or equilibrium phase flash is added.

use super::{
    BinaryCaloricErrorV1, BinaryLiquidCaloricModelV1, PureLiquidHeatCapacityV1, R, checked, poll,
    positive,
};
use crate::entropy::log_ratio;
use fs_math::dd::Dd;
use fs_qty::{AmountConcentration, Temperature};

mod state;
pub use state::{
    BinaryLiquidTransportBalanceV1, BinaryLiquidTransportChangeV1, BinaryLiquidTransportSolveV1,
};

/// Actual homogeneous endpoint in the original species order.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BinaryLiquidTransportPointV1 {
    pub temperature: Temperature,
    pub mole_fractions: [f64; 2],
}

/// Independent kinetic coefficients at the actual face state.
/// D is Maxwell--Stefan diffusivity; k is conductivity at zero diffusion flux.
/// Both may be explicit zero. Q is a signed molar heat-of-transport contrast.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BinaryLiquidTransportCoefficientsV1 {
    pub diffusivity_m2_s: f64,
    pub thermal_conductivity_w_m_k: f64,
    pub heat_of_transport_j_per_mol: f64,
}

/// Left-to-right fluxes with one source-conjugate energy/entropy convention.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BinaryLiquidTransportFaceV1 {
    pub molar_flux_mol_m2_s: [f64; 2],
    pub energy_flux_w_m2: f64,
    pub reduced_heat_flux_w_m2: f64,
    /// Divided original caloric exchange enthalpy h_star, excluding Q.
    /// At coincident endpoints this is the local partial-enthalpy difference.
    pub exchange_enthalpy_j_per_mol: f64,
    /// Composition effort a_chem, excluding the authored Soret force.
    pub chemical_effort_j_per_mol_k: f64,
    /// L=c D x_A_face x_B_face/R, before division by distance.
    pub mobility_mol2_k_j_m_s: f64,
    pub entropy_production_w_m2_k: f64,
}

type Result<T> = core::result::Result<T, BinaryCaloricErrorV1>;

impl BinaryLiquidCaloricModelV1 {
    /// Bound chi on an original source-supported temperature interval and
    /// admit homogeneous spatial transport without gradient regularization.
    /// The only possible interior maximum satisfies H_interaction(T)=0:
    /// chi'=-H_interaction/(R T^2), with affine original H_interaction.
    /// Endpoint checks alone would miss this maximum for positive excess Cp.
    /// This is an analytical-model admission, not an interval uncertainty proof.
    pub fn admit_transport_temperature_interval(&self, bracket: [Temperature; 2]) -> Result<f64> {
        let [lo, hi] = bracket.map(|t| t.value());
        let domain = self.temperature_domain().map(|t| t.value());
        if !positive(lo) || !positive(hi) || lo >= hi || lo < domain[0] || hi > domain[1] {
            return Err(BinaryCaloricErrorV1::Invalid(
                "strict transport temperature interval inside original caloric sources",
            ));
        }
        let pair = self.thermodynamics();
        let mut maximum = pair.chi(bracket[0])?.max(pair.chi(bracket[1])?);
        if pair.excess_heat_capacity_coefficient_j_per_mol_k() > 0.0 {
            let h0 = pair.interaction_enthalpy_j_per_mol(bracket[0])?;
            let h1 = pair.interaction_enthalpy_j_per_mol(bracket[1])?;
            if h0 < 0.0 && h1 > 0.0 {
                let scale = (-h0).max(h1);
                let left = -h0 / scale;
                let right = h1 / scale;
                let temperature = lo + (hi - lo) * (left / (left + right));
                maximum = maximum.max(pair.chi(Temperature::new(temperature))?);
            }
        }
        if maximum > 2.0 {
            return Err(BinaryCaloricErrorV1::Invalid(
                "nonconvex nonisothermal spatial liquids require an interface-energy owner",
            ));
        }
        Ok(maximum)
    }

    /// Evaluate one actual two-point face without allocating or mutating state.
    /// Material source selection belongs to the caller; all original caloric
    /// domains and independent kinetic coefficients remain authoritative here.
    /// Zero kinetic coefficients still return both thermodynamic efforts, so
    /// implicit integrators can reuse this exact source identity in their roots.
    pub fn transport_face(
        &self,
        points: [BinaryLiquidTransportPointV1; 2],
        concentration: AmountConcentration,
        distance_m: f64,
        coefficients: BinaryLiquidTransportCoefficientsV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<BinaryLiquidTransportFaceV1> {
        poll(&mut cancelled)?;
        let c = concentration.value();
        let d = coefficients.diffusivity_m2_s;
        let k = coefficients.thermal_conductivity_w_m_k;
        let q = coefficients.heat_of_transport_j_per_mol;
        if !positive(c)
            || !positive(distance_m)
            || !d.is_finite()
            || d < 0.0
            || !k.is_finite()
            || k < 0.0
            || !q.is_finite()
        {
            return Err(BinaryCaloricErrorV1::Invalid(
                "positive concentration/distance, nonnegative D/k and finite molar heat of transport",
            ));
        }
        let domain = self.temperature_domain().map(|t| t.value());
        let mut fractions = [[0.0; 2]; 2];
        let mut chi = [0.0; 2];
        for i in 0..2 {
            let t = points[i].temperature.value();
            if !positive(t) || t < domain[0] || t > domain[1] {
                return Err(BinaryCaloricErrorV1::Invalid(
                    "face temperature outside original caloric sources",
                ));
            }
            fractions[i] = self
                .thermodynamics()
                .evaluate_composition(points[i].temperature, points[i].mole_fractions)?
                .mole_fractions();
            if fractions[i].iter().any(|&x| !positive(x)) {
                return Err(BinaryCaloricErrorV1::Invalid(
                    "entropy-conjugate liquid faces require common positive species support",
                ));
            }
            chi[i] = self.thermodynamics().chi(points[i].temperature)?;
        }
        let [tl, tr] = points.map(|p| p.temperature.value());
        let delta_beta = reciprocal_difference(tl, tr)?;
        let sampled = Temperature::new(if tl == tr { tl } else { 0.5 * tl + 0.5 * tr });
        let mut h = Dd::ZERO;
        for (i, source) in self.heat_capacities().iter().enumerate() {
            poll(&mut cancelled)?;
            let value = divided_enthalpy(
                source,
                points[0].temperature,
                points[1].temperature,
                sampled,
                delta_beta,
            )?;
            h = if i == 0 { h + value } else { h - value };
        }
        // Delta(R chi)/Delta beta is the original affine interaction
        // enthalpy at the reciprocal-logarithmic mean temperature.
        let t_star = if tl == tr {
            tl
        } else {
            let lo = tl.min(tr);
            let hi = tl.max(tr);
            lo * (log_ratio(hi, lo) / ((hi - lo) / hi))
        };
        let z_bar =
            0.5 * (fractions[0][1] - fractions[0][0]) + 0.5 * (fractions[1][1] - fractions[1][0]);
        h = h + Dd::from_f64(z_bar)
            * Dd::from_f64(
                self.thermodynamics()
                    .interaction_enthalpy_j_per_mol(Temperature::new(t_star))?,
            );
        let exchange_enthalpy_j_per_mol = dd_value(h, "divided liquid exchange enthalpy")?;
        let chemical_effort_j_per_mol_k = checked(
            R * composition_effort(fractions, 0.5 * chi[0] + 0.5 * chi[1]),
            "liquid chemical entropy effort",
        )?;
        let effort = dd_value(
            Dd::from_f64(chemical_effort_j_per_mol_k) - Dd::from_f64(q) * Dd::from_f64(delta_beta),
            "coupled liquid entropy effort",
        )?;
        let mean = [
            0.5 * (fractions[0][0] + fractions[1][0]),
            0.5 * (fractions[0][1] + fractions[1][1]),
        ];
        // Pair the smallest/largest factors before multiplying the middle
        // factor, retaining finite trace mobilities without an occupancy floor.
        let mobility = if d == 0.0 {
            0.0
        } else {
            let mut factors = [c, d, mean[0], mean[1], R.recip()];
            factors.sort_by(f64::total_cmp);
            let a = factors[0] * factors[4];
            let b = factors[1] * factors[3];
            checked(
                (a.max(b) * factors[2]) * a.min(b),
                "liquid Onsager mobility",
            )?
        };
        if d > 0.0 && mobility == 0.0 {
            return Err(BinaryCaloricErrorV1::Arithmetic(
                "positive liquid mobility underflows",
            ));
        }
        let j = checked(-mobility * (effort / distance_m), "liquid molecular flux")?;
        let fourier = checked(-k * ((tr - tl) / distance_m), "liquid Fourier heat")?;
        let reduced = dd_value(
            Dd::from_f64(fourier) + Dd::from_f64(q) * Dd::from_f64(j),
            "liquid reduced heat flux",
        )?;
        let energy = dd_value(
            h * Dd::from_f64(j) + Dd::from_f64(reduced),
            "liquid total energy flux",
        )?;
        let entropy = dd_value(
            -Dd::from_f64(j) * Dd::from_f64(effort)
                + Dd::from_f64(fourier) * Dd::from_f64(delta_beta),
            "liquid entropy production",
        )?;
        if entropy < 0.0 {
            return Err(BinaryCaloricErrorV1::Arithmetic(
                "negative liquid face entropy production",
            ));
        }
        poll(&mut cancelled)?;
        Ok(BinaryLiquidTransportFaceV1 {
            molar_flux_mol_m2_s: [j, -j],
            energy_flux_w_m2: energy,
            reduced_heat_flux_w_m2: reduced,
            exchange_enthalpy_j_per_mol,
            chemical_effort_j_per_mol_k,
            mobility_mol2_k_j_m_s: mobility,
            entropy_production_w_m2_k: entropy,
        })
    }
}

fn divided_enthalpy(
    source: &PureLiquidHeatCapacityV1,
    left: Temperature,
    right: Temperature,
    sampled: Temperature,
    delta_beta: f64,
) -> Result<Dd> {
    let sampled_h =
        Dd::from_f64(source.enthalpy_change_j_per_mol(source.reference_temperature(), sampled)?);
    if left == right {
        return Ok(sampled_h);
    }
    let relative = Dd::from_f64(source.enthalpy_change_j_per_mol(sampled, right)?)
        / Dd::from_f64(right.value())
        - Dd::from_f64(source.enthalpy_change_j_per_mol(sampled, left)?)
            / Dd::from_f64(left.value())
        - Dd::from_f64(source.entropy_change_j_per_mol_k(left, right)?);
    let h = sampled_h + relative / Dd::from_f64(delta_beta);
    dd_value(h, "divided pure-liquid enthalpy")?;
    Ok(h)
}

fn reciprocal_difference(left: f64, right: f64) -> Result<f64> {
    if left == right {
        return Ok(0.0);
    }
    let (lo, hi, sign) = if left < right {
        (left, right, -1.0)
    } else {
        (right, left, 1.0)
    };
    let value = ((hi - lo) / hi) / lo;
    if !positive(value) {
        return Err(BinaryCaloricErrorV1::Arithmetic(
            "reciprocal temperature difference unrepresentable",
        ));
    }
    Ok(sign * value)
}

// A divided log-odds expression retains the critical cubic response instead
// of subtracting rounded, nearly equal chemical potentials. At chi=2 and
// x near 1/2 the linear ideal and interaction terms cancel analytically.
fn composition_effort(x: [[f64; 2]; 2], chi: f64) -> f64 {
    let zl = x[0][0] - x[0][1];
    let zr = x[1][0] - x[1][1];
    let dz = (Dd::from_f64(x[1][0]) - Dd::from_f64(x[0][0]))
        - (Dd::from_f64(x[1][1]) - Dd::from_f64(x[0][1]));
    if zl.abs().max(zr.abs()) < 0.125 {
        let product = zl * zr;
        let denominator = 1.0 - product;
        let u = dz.to_f64() / denominator;
        let u2 = u * u;
        let mut tail = 0.0;
        for n in (1..=20).rev() {
            tail = 1.0 / f64::from(2 * n + 1) + u2 * tail;
        }
        dz.to_f64() * ((2.0 - chi) + 2.0 * (product + u2 * tail) / denominator)
    } else {
        (Dd::from_f64(log_ratio(x[1][0], x[0][0]))
            - Dd::from_f64(log_ratio(x[1][1], x[0][1]))
            - Dd::from_f64(chi) * dz)
            .to_f64()
    }
}

fn dd_value(value: Dd, label: &'static str) -> Result<f64> {
    if !value.hi.is_finite() || !value.lo.is_finite() {
        return Err(BinaryCaloricErrorV1::Arithmetic(label));
    }
    checked(value.to_f64(), label)
}

#[cfg(test)]
mod tests;
