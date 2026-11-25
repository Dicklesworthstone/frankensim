//! Asymmetric multicomponent UNIQUAC excess-Gibbs liquid solution model.
//!
//! Implements the Abrams-Prausnitz (1975) Universal Quasi-Chemical model with
//! separate combinatorial (molecular size/shape) and residual (energetic
//! interaction) contributions:
//!
//! - `g^E / (R T) = g^{E,comb} / (R T) + g^{E,res} / (R T)`
//! - `tau_ij(T) = exp(a_ij + b_ij / T)` (`a_ii = 0`, `b_ii = 0`, so `tau_ii = 1`)
//! - Supports both classic `q'_i = q_i` (Abrams-Prausnitz) and modified residual
//!   surface parameters `q'_i != q_i` (Anderson-Prausnitz 1978 for alcohols/water).
//!
//! Both combinatorial and residual terms are homogeneous-degree-0 excess Gibbs
//! potentials and independently satisfy Euler's theorem and the Gibbs-Duhem
//! identity `sum_i x_i d(ln gamma_i) = 0`. Ratio evaluations use `r_i / R_mix`
//! and `q_i / Q_mix` directly, avoiding any `0/0` singularity at `x_i = 0`.
//!
//! All transcendentals route strictly through `fs_math::det`.

use crate::solution::{
    ExcessGibbsErrorV1, ExcessGibbsEvaluationV1, ExcessGibbsModelV1, MAX_EXCESS_GIBBS_SPECIES_V1,
    ModifiedRaoultVleV1, VaporPressureLawV1, check_finite, check_temperature_in_domain,
    validate_and_normalize_mole_fractions, validate_distinct_species, validate_temperature_domain,
};
use crate::{SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy};
use fs_math::det::{exp, ln};
use fs_qty::{Dims, Temperature};

/// Standard Abrams-Prausnitz lattice coordination number (`z = 10`).
pub const UNIQUAC_COORDINATION_NUMBER_V1: f64 = 10.0;
/// Exact `fs-matdb` law id for a UNIQUAC excess-Gibbs card.
pub const UNIQUAC_LAW_ID_V1: &str = "uniquac-excess-gibbs";
/// Supported semantic version for UNIQUAC `fs-matdb` cards.
pub const UNIQUAC_LAW_VERSION_V1: u32 = 1;

/// Pure-component structural parameters `(r_i, q_i, q'_i)` for UNIQUAC.
#[derive(Clone, Debug, PartialEq)]
pub struct UniquacComponentV1 {
    species: SpeciesId,
    volume_r: f64,
    area_q: f64,
    residual_area_q_prime: f64,
}

impl UniquacComponentV1 {
    /// Construct a standard Abrams-Prausnitz component with `q'_i = q_i`.
    pub fn standard(
        species: SpeciesId,
        volume_r: f64,
        area_q: f64,
    ) -> Result<Self, ExcessGibbsErrorV1> {
        Self::modified(species, volume_r, area_q, area_q)
    }

    /// Construct an Anderson-Prausnitz component with distinct residual surface area `q'_i`.
    pub fn modified(
        species: SpeciesId,
        volume_r: f64,
        area_q: f64,
        residual_area_q_prime: f64,
    ) -> Result<Self, ExcessGibbsErrorV1> {
        if !volume_r.is_finite()
            || volume_r <= 0.0
            || !area_q.is_finite()
            || area_q <= 0.0
            || !residual_area_q_prime.is_finite()
            || residual_area_q_prime <= 0.0
        {
            return Err(ExcessGibbsErrorV1::InvalidModel(
                "UNIQUAC structural parameters r, q, q' must be finite and strictly positive",
            ));
        }
        Ok(Self {
            species,
            volume_r,
            area_q,
            residual_area_q_prime,
        })
    }

    /// Component species identity.
    #[must_use]
    pub const fn species(&self) -> &SpeciesId {
        &self.species
    }

    /// van der Waals molecular volume parameter `r_i`.
    #[must_use]
    pub const fn volume_r(&self) -> f64 {
        self.volume_r
    }

    /// van der Waals molecular surface-area parameter `q_i`.
    #[must_use]
    pub const fn area_q(&self) -> f64 {
        self.area_q
    }

    /// Residual surface-area parameter `q'_i` (equals `q_i` for standard UNIQUAC).
    #[must_use]
    pub const fn residual_area_q_prime(&self) -> f64 {
        self.residual_area_q_prime
    }
}

/// Combinatorial and residual breakdown alongside the full combined UNIQUAC evaluation.
#[derive(Clone, Debug, PartialEq)]
pub struct UniquacSplitEvaluationV1 {
    /// Combined excess-Gibbs evaluation (`combinatorial + residual`).
    pub total: ExcessGibbsEvaluationV1,
    /// Combinatorial reduced excess Gibbs energy `g^{E,comb} / (R T)`.
    pub g_ex_combinatorial_over_rt: f64,
    /// Residual reduced excess Gibbs energy `g^{E,res} / (R T)`.
    pub g_ex_residual_over_rt: f64,
    /// Combinatorial log activity coefficients `ln gamma_i^{comb}`.
    pub ln_gamma_combinatorial: Vec<f64>,
    /// Residual log activity coefficients `ln gamma_i^{res}`.
    pub ln_gamma_residual: Vec<f64>,
}

/// Pure validated multicomponent UNIQUAC excess-Gibbs model.
#[derive(Clone, Debug, PartialEq)]
pub struct UniquacModelV1 {
    species: Vec<SpeciesId>,
    components: Vec<UniquacComponentV1>,
    a_matrix: Vec<f64>,
    b_matrix_k: Vec<f64>,
    coordination_number: f64,
    temperature_min: Temperature,
    temperature_max: Temperature,
}

impl UniquacModelV1 {
    /// Construct a validated multicomponent UNIQUAC model with `z = 10`.
    ///
    /// Interaction convention: `tau_ij(T) = exp(a_ij + b_ij / T)` with `a_ii = 0`,
    /// `b_ii = 0`. Note that in terms of interaction energy differences
    /// `Delta u_ij = u_ij - u_jj`, `b_ij = -Delta u_ij / R`.
    pub fn new(
        components: Vec<UniquacComponentV1>,
        a_matrix: Vec<f64>,
        b_matrix_k: Vec<f64>,
        minimum_temperature: Temperature,
        maximum_temperature: Temperature,
    ) -> Result<Self, ExcessGibbsErrorV1> {
        let n = components.len();
        if !(1..=MAX_EXCESS_GIBBS_SPECIES_V1).contains(&n) {
            return Err(ExcessGibbsErrorV1::InvalidModel(
                "UNIQUAC component count must be in 1..=64",
            ));
        }
        let species: Vec<SpeciesId> = components.iter().map(|c| c.species.clone()).collect();
        validate_distinct_species(&species)?;
        validate_temperature_domain(minimum_temperature, maximum_temperature)?;

        let expected_len = n * n;
        if a_matrix.len() != expected_len || b_matrix_k.len() != expected_len {
            return Err(ExcessGibbsErrorV1::InvalidModel(
                "UNIQUAC parameter matrices must have length n * n",
            ));
        }
        if a_matrix.iter().any(|v| !v.is_finite()) || b_matrix_k.iter().any(|v| !v.is_finite()) {
            return Err(ExcessGibbsErrorV1::InvalidModel(
                "UNIQUAC parameter matrices must contain only finite values",
            ));
        }
        for i in 0..n {
            let diag = i * n + i;
            if a_matrix[diag] != 0.0 || b_matrix_k[diag] != 0.0 {
                return Err(ExcessGibbsErrorV1::InvalidModel(
                    "UNIQUAC parameter matrices must have zero diagonal entries",
                ));
            }
        }

        Ok(Self {
            species,
            components,
            a_matrix,
            b_matrix_k,
            coordination_number: UNIQUAC_COORDINATION_NUMBER_V1,
            temperature_min: minimum_temperature,
            temperature_max: maximum_temperature,
        })
    }

    /// Construct a binary UNIQUAC model from scalar off-diagonal interaction parameters.
    #[allow(clippy::too_many_arguments)]
    pub fn binary(
        components: [UniquacComponentV1; 2],
        a12: f64,
        a21: f64,
        b12_k: f64,
        b21_k: f64,
        minimum_temperature: Temperature,
        maximum_temperature: Temperature,
    ) -> Result<Self, ExcessGibbsErrorV1> {
        Self::new(
            components.into(),
            vec![0.0, a12, a21, 0.0],
            vec![0.0, b12_k, b21_k, 0.0],
            minimum_temperature,
            maximum_temperature,
        )
    }

    /// Load a validated UNIQUAC model from an immutable `fs-matdb` constitutive card.
    ///
    /// Expected card schema:
    /// - `law = "uniquac-excess-gibbs"`, `law_version = 1`, stateless
    /// - validity axis `"T"` in kelvin
    /// - for every ordered pair `i != j`: `a_{i}_{j}` (dimensionless), `b_{i}_{j}` [K]
    pub fn from_card(
        card: &ConstitutiveModelCard,
        components: Vec<UniquacComponentV1>,
    ) -> Result<Self, ExcessGibbsErrorV1> {
        card.validate().map_err(|_| {
            ExcessGibbsErrorV1::InvalidModel("UNIQUAC matdb card failed validation")
        })?;
        if card.law.0 != UNIQUAC_LAW_ID_V1 {
            return Err(ExcessGibbsErrorV1::InvalidModel(
                "UNIQUAC matdb card law id mismatch",
            ));
        }
        if card.law_version != UNIQUAC_LAW_VERSION_V1 {
            return Err(ExcessGibbsErrorV1::InvalidModel(
                "UNIQUAC matdb card law version mismatch",
            ));
        }
        if card.state_schema_version != 0
            || card.initial_state != InitialStatePolicy::ZeroInternalState
        {
            return Err(ExcessGibbsErrorV1::InvalidModel(
                "UNIQUAC matdb card must be stateless",
            ));
        }
        let (t_min, t_max) = card
            .validity
            .bound("T")
            .ok_or(ExcessGibbsErrorV1::InvalidModel(
                "UNIQUAC matdb card missing T validity axis",
            ))?;
        let n = components.len();
        if !(1..=MAX_EXCESS_GIBBS_SPECIES_V1).contains(&n) {
            return Err(ExcessGibbsErrorV1::InvalidModel(
                "UNIQUAC component count must be in 1..=64",
            ));
        }
        let mut a_matrix = vec![0.0; n * n];
        let mut b_matrix_k = vec![0.0; n * n];
        for i in 0..n {
            for j in 0..n {
                if i == j {
                    continue;
                }
                let a_key = format!("a_{i}_{j}");
                let b_key = format!("b_{i}_{j}");
                let a_param =
                    card.parameters
                        .get(&a_key)
                        .ok_or(ExcessGibbsErrorV1::InvalidModel(
                            "UNIQUAC matdb card missing a_i_j parameter",
                        ))?;
                if a_param.dims != Dims::NONE {
                    return Err(ExcessGibbsErrorV1::InvalidModel(
                        "UNIQUAC matdb card parameter a_i_j must be dimensionless",
                    ));
                }
                let b_param =
                    card.parameters
                        .get(&b_key)
                        .ok_or(ExcessGibbsErrorV1::InvalidModel(
                            "UNIQUAC matdb card missing b_i_j parameter",
                        ))?;
                if b_param.dims != Temperature::DIMS {
                    return Err(ExcessGibbsErrorV1::InvalidModel(
                        "UNIQUAC matdb card parameter b_i_j must have temperature dimensions [K]",
                    ));
                }
                a_matrix[i * n + j] = a_param.value;
                b_matrix_k[i * n + j] = b_param.value;
            }
        }
        Self::new(
            components,
            a_matrix,
            b_matrix_k,
            Temperature::new(t_min),
            Temperature::new(t_max),
        )
    }

    /// Declared species identities in matrix row/column order.
    #[must_use]
    pub fn species(&self) -> &[SpeciesId] {
        &self.species
    }

    /// Declared structural components `(r_i, q_i, q'_i)`.
    #[must_use]
    pub fn components(&self) -> &[UniquacComponentV1] {
        &self.components
    }

    /// Inclusive valid temperature interval `[T_min, T_max]` [K].
    #[must_use]
    pub const fn temperature_domain(&self) -> [Temperature; 2] {
        [self.temperature_min, self.temperature_max]
    }

    /// Lattice coordination number `z` (default `10.0`).
    #[must_use]
    pub const fn coordination_number(&self) -> f64 {
        self.coordination_number
    }

    /// Evaluate UNIQUAC returning both the combined excess-Gibbs state and the
    /// separate combinatorial/residual contributions.
    pub fn evaluate_split(
        &self,
        temperature: Temperature,
        mole_fractions: &[f64],
    ) -> Result<UniquacSplitEvaluationV1, ExcessGibbsErrorV1> {
        let t = check_temperature_in_domain(
            temperature,
            self.temperature_min.value(),
            self.temperature_max.value(),
        )?;
        let n = self.components.len();
        let x = validate_and_normalize_mole_fractions(mole_fractions, n)?;

        let z_half = 0.5 * self.coordination_number;
        let inv_t = 1.0 / t;
        let inv_t2 = inv_t * inv_t;

        // Mixture structural sums R_mix, Q_mix, Q'_mix, and Staverman-Guggenheim l_i.
        let mut r_mix = 0.0;
        let mut q_mix = 0.0;
        let mut qp_mix = 0.0;
        let mut l_vec = vec![0.0; n];
        let mut l_bar = 0.0;
        for i in 0..n {
            let c = &self.components[i];
            let li = z_half * (c.volume_r - c.area_q) - (c.volume_r - 1.0);
            l_vec[i] = li;
            r_mix += c.volume_r * x[i];
            q_mix += c.area_q * x[i];
            qp_mix += c.residual_area_q_prime * x[i];
            l_bar += x[i] * li;
        }
        if !r_mix.is_finite()
            || r_mix <= 0.0
            || !q_mix.is_finite()
            || q_mix <= 0.0
            || !qp_mix.is_finite()
            || qp_mix <= 0.0
        {
            return Err(ExcessGibbsErrorV1::Unrepresentable(
                "UNIQUAC mixture structural sum is non-positive or non-finite",
            ));
        }

        // Combinatorial ratios rho_i = r_i / R_mix, sigma_i = q_i / Q_mix (safe at x_i = 0).
        let mut rho = vec![0.0; n];
        let mut sigma = vec![0.0; n];
        let mut ln_gamma_comb = vec![0.0; n];
        let mut g_ex_comb_over_rt = 0.0;
        for i in 0..n {
            let c = &self.components[i];
            let rho_i = check_finite(c.volume_r / r_mix, "UNIQUAC rho_i")?;
            let sig_i = check_finite(c.area_q / q_mix, "UNIQUAC sigma_i")?;
            let ratio_sr = check_finite(sig_i / rho_i, "UNIQUAC sigma_i/rho_i")?;
            let ln_rho = ln(rho_i);
            let ln_sr = ln(ratio_sr);
            let lng_c = check_finite(
                ln_rho + z_half * c.area_q * ln_sr + l_vec[i] - rho_i * l_bar,
                "UNIQUAC ln gamma_i^comb",
            )?;
            rho[i] = rho_i;
            sigma[i] = sig_i;
            ln_gamma_comb[i] = lng_c;
            g_ex_comb_over_rt += x[i] * (ln_rho + z_half * c.area_q * ln_sr);
        }
        let g_ex_comb_over_rt = check_finite(g_ex_comb_over_rt, "UNIQUAC g^{E,comb}/(RT)")?;

        // Residual interaction matrix tau_ij = exp(a_ij + b_ij / T).
        let mut tau = vec![0.0; n * n];
        for i in 0..n {
            for j in 0..n {
                let idx = i * n + j;
                if i == j {
                    tau[idx] = 1.0;
                } else {
                    let arg = check_finite(
                        self.a_matrix[idx] + self.b_matrix_k[idx] * inv_t,
                        "UNIQUAC tau exponent",
                    )?;
                    let t_ij = exp(arg);
                    if !t_ij.is_finite() || t_ij <= 0.0 {
                        return Err(ExcessGibbsErrorV1::Unrepresentable(
                            "UNIQUAC tau_ij overflowed or underflowed",
                        ));
                    }
                    tau[idx] = t_ij;
                }
            }
        }

        // Column sums W_k = sum_j q'_j x_j tau_jk and W_tilde_k = sum_j q'_j x_j b_jk tau_jk.
        let mut w_col = vec![0.0; n];
        let mut w_tilde = vec![0.0; n];
        for k in 0..n {
            let mut wk = 0.0;
            let mut wtk = 0.0;
            for j in 0..n {
                let idx = j * n + k;
                let qx_tau = self.components[j].residual_area_q_prime * x[j] * tau[idx];
                wk += qx_tau;
                wtk += qx_tau * self.b_matrix_k[idx];
            }
            if !wk.is_finite() || wk <= 0.0 {
                return Err(ExcessGibbsErrorV1::Unrepresentable(
                    "UNIQUAC residual sum W_k is non-positive or non-finite",
                ));
            }
            w_col[k] = wk;
            w_tilde[k] = check_finite(wtk, "UNIQUAC dW_k/dT")?;
        }

        // Residual excess Gibbs energy, enthalpy, activity coefficients, and T derivatives.
        let mut g_ex_res_over_rt = 0.0;
        let mut h_ex_over_r = 0.0;
        let mut ln_gamma_res = vec![0.0; n];
        let mut ln_gamma = vec![0.0; n];
        let mut gamma = vec![0.0; n];
        let mut d_ln_gamma_dt = vec![0.0; n];
        let mut partial_h_ex = vec![0.0; n];

        for i in 0..n {
            let qp_i = self.components[i].residual_area_q_prime;
            let s_i = check_finite(w_col[i] / qp_mix, "UNIQUAC S_i")?;
            let ln_si = ln(s_i);
            let wt_over_w_i = check_finite(w_tilde[i] / w_col[i], "UNIQUAC W_tilde_i / W_i")?;

            g_ex_res_over_rt -= qp_i * x[i] * ln_si;
            // h^E / R = -sum_k q'_k x_k (W_tilde_k / W_k)
            h_ex_over_r -= qp_i * x[i] * wt_over_w_i;

            let mut res_sum = 0.0;
            let mut res_t_sum = 0.0;
            for k in 0..n {
                let idx = i * n + k;
                let weight =
                    (self.components[k].residual_area_q_prime * x[k] * tau[idx]) / w_col[k];
                res_sum += weight;
                res_t_sum += weight * (self.b_matrix_k[idx] - w_tilde[k] / w_col[k]);
            }
            let lng_r = check_finite(qp_i * (1.0 - ln_si - res_sum), "UNIQUAC ln gamma_i^res")?;
            let lng = check_finite(ln_gamma_comb[i] + lng_r, "UNIQUAC ln gamma_i")?;
            let gam = exp(lng);
            if !gam.is_finite() || gam <= 0.0 {
                return Err(ExcessGibbsErrorV1::Unrepresentable(
                    "UNIQUAC activity coefficient gamma_i overflowed or underflowed",
                ));
            }
            // -T^2 * d(ln gamma_i)/dT = q'_i * [ -W_tilde_i/W_i - sum_k (q'_k x_k tau_ik / W_k)*(b_ik - W_tilde_k/W_k) ]
            let neg_t2_dlng_dt = check_finite(
                qp_i * (-wt_over_w_i - res_t_sum),
                "UNIQUAC partial excess enthalpy / R",
            )?;
            let dlng_dt = check_finite(-neg_t2_dlng_dt * inv_t2, "UNIQUAC d(ln gamma_i)/dT")?;
            let hi_ex = check_finite(R * neg_t2_dlng_dt, "UNIQUAC partial excess enthalpy")?;

            ln_gamma_res[i] = lng_r;
            ln_gamma[i] = lng;
            gamma[i] = gam;
            d_ln_gamma_dt[i] = dlng_dt;
            partial_h_ex[i] = hi_ex;
        }

        let g_ex_res_over_rt = check_finite(g_ex_res_over_rt, "UNIQUAC g^{E,res}/(RT)")?;
        let g_ex_over_rt = check_finite(g_ex_comb_over_rt + g_ex_res_over_rt, "UNIQUAC g^E/(RT)")?;
        let g_ex = check_finite(g_ex_over_rt * R * t, "UNIQUAC g^E")?;
        let h_ex = check_finite(h_ex_over_r * R, "UNIQUAC h^E")?;
        let s_ex = check_finite(R * (h_ex_over_r * inv_t - g_ex_over_rt), "UNIQUAC s^E")?;

        // Exact symmetric homogeneous-degree-0 Jacobian:
        // d(ln gamma_i^comb)/dx_j = (1 - rho_i)*(1 - rho_j)
        //   + (z/2)*(q_i*rho_j + q_j*rho_i - q_i*sigma_j - Q_mix*rho_i*rho_j)
        // d(ln gamma_i^res)/dx_j = q'_i*q'_j * [ 1/Q'_mix - tau_ji/W_i - tau_ij/W_j
        //   + sum_k (q'_k*x_k*tau_ik*tau_jk / W_k^2) ]
        let mut d_ln_gamma_dx = vec![0.0; n * n];
        let inv_qp_mix = 1.0 / qp_mix;
        for i in 0..n {
            let qi = self.components[i].area_q;
            let qpi = self.components[i].residual_area_q_prime;
            for j in i..n {
                let qj = self.components[j].area_q;
                let qpj = self.components[j].residual_area_q_prime;

                let comb_ij = (1.0 - rho[i]) * (1.0 - rho[j])
                    + z_half
                        * (qi * rho[j] + qj * rho[i] - qi * sigma[j] - q_mix * rho[i] * rho[j]);

                let mut res_sum_k = 0.0;
                for k in 0..n {
                    let wk = w_col[k];
                    let qpk = self.components[k].residual_area_q_prime;
                    res_sum_k += qpk * x[k] * (tau[i * n + k] / wk) * (tau[j * n + k] / wk);
                }
                let res_ij = (qpi * qpj)
                    * (inv_qp_mix - tau[j * n + i] / w_col[i] - tau[i * n + j] / w_col[j]
                        + res_sum_k);

                let entry = check_finite(comb_ij + res_ij, "UNIQUAC d(ln gamma_i)/dx_j")?;
                d_ln_gamma_dx[i * n + j] = entry;
                d_ln_gamma_dx[j * n + i] = entry;
            }
        }

        let total = ExcessGibbsEvaluationV1 {
            temperature,
            mole_fractions: x,
            reduced_excess_gibbs: g_ex_over_rt,
            excess_gibbs_j_per_mol: g_ex,
            excess_enthalpy_j_per_mol: h_ex,
            excess_entropy_j_per_mol_k: s_ex,
            ln_gamma,
            gamma,
            d_ln_gamma_dt,
            partial_excess_enthalpy_j_per_mol: partial_h_ex,
            d_ln_gamma_dx,
        };

        Ok(UniquacSplitEvaluationV1 {
            total,
            g_ex_combinatorial_over_rt: g_ex_comb_over_rt,
            g_ex_residual_over_rt: g_ex_res_over_rt,
            ln_gamma_combinatorial: ln_gamma_comb,
            ln_gamma_residual: ln_gamma_res,
        })
    }

    /// Evaluate combined UNIQUAC excess Gibbs properties at `(temperature, mole_fractions)`.
    pub fn evaluate(
        &self,
        temperature: Temperature,
        mole_fractions: &[f64],
    ) -> Result<ExcessGibbsEvaluationV1, ExcessGibbsErrorV1> {
        Ok(self.evaluate_split(temperature, mole_fractions)?.total)
    }

    /// Pre-validated Ethanol(0) - Water(1) UNIQUAC parameter set (Aspen Plus APV88/VLE-IG).
    ///
    /// Structural parameters: Ethanol (`r = 2.1055`, `q = 1.972`), Water (`r = 0.92`, `q = 1.40`);
    /// interaction parameters: `a12 = -2.4936`, `a21 = 2.0046`, `b12 = 756.9477 K`,
    /// `b21 = -728.9705 K`, valid on `273.15..=423.15 K`.
    #[must_use]
    pub fn ethanol_water_aspen_v1() -> Self {
        let components = [
            UniquacComponentV1::standard(
                SpeciesId::new("C2H5OH").expect("valid literal species id"),
                2.1055,
                1.972,
            )
            .expect("valid ethanol structural parameters"),
            UniquacComponentV1::standard(
                SpeciesId::new("H2O").expect("valid literal species id"),
                0.92,
                1.40,
            )
            .expect("valid water structural parameters"),
        ];
        Self::binary(
            components,
            -2.4936,
            2.0046,
            756.9477,
            -728.9705,
            Temperature::new(273.15),
            Temperature::new(423.15),
        )
        .expect("published Aspen Plus ethanol-water UNIQUAC parameters are valid")
    }

    /// Pre-validated Ethanol(0) - Water(1) Modified Raoult VLE model backed by UNIQUAC.
    #[must_use]
    pub fn ethanol_water_vle_v1() -> ModifiedRaoultVleV1<Self> {
        ModifiedRaoultVleV1::new(
            Self::ethanol_water_aspen_v1(),
            vec![
                VaporPressureLawV1::ethanol_antoine_nist(),
                VaporPressureLawV1::water_antoine_dechema(),
            ],
        )
        .expect("published UNIQUAC ethanol-water VLE model is valid")
    }
}

impl ExcessGibbsModelV1 for UniquacModelV1 {
    fn species(&self) -> &[SpeciesId] {
        &self.species
    }

    fn temperature_domain(&self) -> [Temperature; 2] {
        [self.temperature_min, self.temperature_max]
    }

    fn evaluate(
        &self,
        temperature: Temperature,
        mole_fractions: &[f64],
    ) -> Result<ExcessGibbsEvaluationV1, ExcessGibbsErrorV1> {
        Self::evaluate(self, temperature, mole_fractions)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn uniquac_pure_limits_and_infinite_dilution_match_exact_formulas() {
        let model = UniquacModelV1::ethanol_water_aspen_v1();
        let t = Temperature::new(351.30);

        let pure0 = model.evaluate_split(t, &[1.0, 0.0]).unwrap();
        assert!(pure0.total.reduced_excess_gibbs().abs() < 1e-15);
        assert!(pure0.g_ex_combinatorial_over_rt.abs() < 1e-15);
        assert!(pure0.g_ex_residual_over_rt.abs() < 1e-15);
        assert!((pure0.total.activity_coefficients()[0] - 1.0).abs() < 1e-15);

        let tau12 = exp(-2.4936 + 756.9477 / t.value());
        let tau21 = exp(2.0046 - 728.9705 / t.value());
        let q2 = 1.40;
        let expected_res2_inf = q2 * (1.0 - ln(tau12) - tau21);
        assert!((pure0.ln_gamma_residual[1] - expected_res2_inf).abs() < 1e-14);

        let pure1 = model.evaluate_split(t, &[0.0, 1.0]).unwrap();
        assert!((pure1.total.activity_coefficients()[1] - 1.0).abs() < 1e-15);
        let q1 = 1.972;
        let expected_res1_inf = q1 * (1.0 - ln(tau21) - tau12);
        assert!((pure1.ln_gamma_residual[0] - expected_res1_inf).abs() < 1e-14);
    }

    #[test]
    fn uniquac_gibbs_duhem_and_derivatives_match_ternary_with_modified_q_prime() {
        let comps = vec![
            UniquacComponentV1::modified(SpeciesId::new("C2H5OH").unwrap(), 2.1055, 1.972, 0.92)
                .unwrap(),
            UniquacComponentV1::modified(SpeciesId::new("H2O").unwrap(), 0.92, 1.40, 1.00).unwrap(),
            UniquacComponentV1::standard(SpeciesId::new("C3H6O").unwrap(), 2.5735, 2.336).unwrap(),
        ];
        let a = vec![0.0, -0.5, 0.2, 0.7, 0.0, -0.3, -0.1, 0.4, 0.0];
        let b = vec![0.0, 180.0, -90.0, -240.0, 0.0, 150.0, 110.0, -170.0, 0.0];
        let model = UniquacModelV1::new(
            comps,
            a,
            b,
            Temperature::new(273.15),
            Temperature::new(450.0),
        )
        .unwrap();

        let t = Temperature::new(330.0);
        let x = [0.33, 0.45, 0.22];
        let split = model.evaluate_split(t, &x).unwrap();
        let eval = &split.total;

        // Euler identities for combinatorial, residual, and total.
        let euler_comb: f64 = x
            .iter()
            .zip(&split.ln_gamma_combinatorial)
            .map(|(xi, lg)| xi * lg)
            .sum();
        let euler_res: f64 = x
            .iter()
            .zip(&split.ln_gamma_residual)
            .map(|(xi, lg)| xi * lg)
            .sum();
        assert!((euler_comb - split.g_ex_combinatorial_over_rt).abs() < 1e-14);
        assert!((euler_res - split.g_ex_residual_over_rt).abs() < 1e-14);

        // Exact Gibbs-Duhem column residuals < 1e-13.
        for (j, res) in eval.gibbs_duhem_column_residuals().into_iter().enumerate() {
            assert!(res.abs() < 1e-13, "col {j} Gibbs-Duhem residual {res:e}");
        }

        // Finite-difference d(ln gamma_i)/dT and h^E.
        let dt = 1e-4;
        let plus = model
            .evaluate(Temperature::new(t.value() + dt), &x)
            .unwrap();
        let minus = model
            .evaluate(Temperature::new(t.value() - dt), &x)
            .unwrap();
        for i in 0..3 {
            let fd = (plus.ln_activity_coefficients()[i] - minus.ln_activity_coefficients()[i])
                / (2.0 * dt);
            assert!(
                (eval.d_ln_gamma_dt()[i] - fd).abs() < 1e-9,
                "species {i}: analytic={}, fd={fd}",
                eval.d_ln_gamma_dt()[i]
            );
        }
        let h_fd = -R
            * t.value()
            * t.value()
            * (plus.reduced_excess_gibbs() - minus.reduced_excess_gibbs())
            / (2.0 * dt);
        assert!(
            (eval.excess_enthalpy_j_per_mol() - h_fd).abs() < 1e-5,
            "h^E analytic={}, fd={h_fd}",
            eval.excess_enthalpy_j_per_mol()
        );

        // Directional simplex Jacobian check.
        let dir = [0.7, -0.5, -0.2];
        let eps = 1e-6;
        let xp = [
            x[0] + eps * dir[0],
            x[1] + eps * dir[1],
            x[2] + eps * dir[2],
        ];
        let xm = [
            x[0] - eps * dir[0],
            x[1] - eps * dir[1],
            x[2] - eps * dir[2],
        ];
        let ep = model.evaluate(t, &xp).unwrap();
        let em = model.evaluate(t, &xm).unwrap();
        for i in 0..3 {
            let analytic: f64 = (0..3)
                .map(|j| eval.d_ln_gamma_dx_entry(i, j).unwrap() * dir[j])
                .sum();
            let fd =
                (ep.ln_activity_coefficients()[i] - em.ln_activity_coefficients()[i]) / (2.0 * eps);
            assert!(
                (analytic - fd).abs() < 1e-8,
                "directional d(ln gamma_{i}): analytic={analytic}, fd={fd}"
            );
        }
    }
}
