//! Asymmetric multicomponent NRTL excess-Gibbs liquid solution model.
//!
//! Implements the Renon-Prausnitz (1968) Non-Random Two-Liquid model:
//!
//! - `tau_ij(T) = a_ij + b_ij / T` (`a_ii = 0`, `b_ii = 0`)
//! - `G_ij(T) = exp(-alpha_ij * tau_ij(T))` (`alpha_ij = alpha_ji`, `alpha_ii = 0`)
//! - `g^E / (R T) = sum_i x_i * (sum_j tau_ji G_ji x_j) / (sum_k G_ki x_k)`
//!
//! Activity coefficients `ln gamma_i`, excess enthalpy `h^E = -R T^2 d(g^E/RT)/dT`,
//! temperature derivatives `d ln gamma_i / dT`, and the symmetric homogeneous-degree-0
//! composition Jacobian `d ln gamma_i / d x_j` all derive analytically from the same
//! single Gibbs excess potential, satisfying Euler's theorem and the Gibbs-Duhem
//! identity `sum_i x_i d(ln gamma_i) = 0` to floating-point roundoff.
//!
//! All transcendentals route strictly through `fs_math::det` for deterministic
//! evaluation with zero external runtime dependencies.

use crate::solution::{
    ExcessGibbsErrorV1, ExcessGibbsEvaluationV1, ExcessGibbsModelV1, MAX_EXCESS_GIBBS_SPECIES_V1,
    ModifiedRaoultVleV1, VaporPressureLawV1, check_finite, check_temperature_in_domain,
    validate_and_normalize_mole_fractions, validate_distinct_species, validate_temperature_domain,
};
use crate::{SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy};
use fs_math::det::exp;
use fs_qty::{Dims, Temperature};

/// Exact `fs-matdb` law id for an NRTL excess-Gibbs card.
pub const NRTL_LAW_ID_V1: &str = "nrtl-excess-gibbs";
/// Supported semantic version for NRTL `fs-matdb` cards.
pub const NRTL_LAW_VERSION_V1: u32 = 1;

/// Pure validated NRTL parameter set and immutable multicomponent evaluator.
///
/// Row-major `n x n` parameter matrices follow the declared species order:
/// index `i * n + j` stores the `(i, j)` interaction entering `tau_ij = a_ij + b_ij / T`
/// and `G_ij = exp(-alpha_ij * tau_ij)`.
#[derive(Clone, Debug, PartialEq)]
pub struct NrtlModelV1 {
    species: Vec<SpeciesId>,
    a_matrix: Vec<f64>,
    b_matrix_k: Vec<f64>,
    alpha_matrix: Vec<f64>,
    temperature_min: Temperature,
    temperature_max: Temperature,
}

impl NrtlModelV1 {
    /// Construct a validated multicomponent NRTL excess-Gibbs model.
    ///
    /// # Conventions
    /// - `species` must contain `1..=MAX_EXCESS_GIBBS_SPECIES_V1` distinct species.
    /// - `a_matrix` (dimensionless), `b_matrix_k` [K], and `alpha_matrix` (dimensionless)
    ///   must each have length `n * n` with finite entries and zero diagonals.
    /// - `alpha_matrix` must be symmetric (`alpha_ij == alpha_ji`); `a_matrix` and
    ///   `b_matrix_k` are generally asymmetric.
    /// - `minimum_temperature` and `maximum_temperature` define an inclusive positive
    ///   kelvin validity interval.
    pub fn new(
        species: Vec<SpeciesId>,
        a_matrix: Vec<f64>,
        b_matrix_k: Vec<f64>,
        alpha_matrix: Vec<f64>,
        minimum_temperature: Temperature,
        maximum_temperature: Temperature,
    ) -> Result<Self, ExcessGibbsErrorV1> {
        let n = species.len();
        if !(1..=MAX_EXCESS_GIBBS_SPECIES_V1).contains(&n) {
            return Err(ExcessGibbsErrorV1::InvalidModel(
                "NRTL species count must be in 1..=64",
            ));
        }
        validate_distinct_species(&species)?;
        validate_temperature_domain(minimum_temperature, maximum_temperature)?;

        let expected_len = n * n;
        if a_matrix.len() != expected_len
            || b_matrix_k.len() != expected_len
            || alpha_matrix.len() != expected_len
        {
            return Err(ExcessGibbsErrorV1::InvalidModel(
                "NRTL parameter matrices must have length n * n",
            ));
        }
        if a_matrix.iter().any(|v| !v.is_finite())
            || b_matrix_k.iter().any(|v| !v.is_finite())
            || alpha_matrix.iter().any(|v| !v.is_finite())
        {
            return Err(ExcessGibbsErrorV1::InvalidModel(
                "NRTL parameter matrices must contain only finite values",
            ));
        }
        for i in 0..n {
            let diag = i * n + i;
            if a_matrix[diag] != 0.0 || b_matrix_k[diag] != 0.0 || alpha_matrix[diag] != 0.0 {
                return Err(ExcessGibbsErrorV1::InvalidModel(
                    "NRTL parameter matrices must have zero diagonal entries",
                ));
            }
            for j in (i + 1)..n {
                if alpha_matrix[i * n + j].to_bits() != alpha_matrix[j * n + i].to_bits() {
                    return Err(ExcessGibbsErrorV1::InvalidModel(
                        "NRTL non-randomness matrix alpha must be symmetric",
                    ));
                }
            }
        }

        Ok(Self {
            species,
            a_matrix,
            b_matrix_k,
            alpha_matrix,
            temperature_min: minimum_temperature,
            temperature_max: maximum_temperature,
        })
    }

    /// Construct a binary NRTL model from scalar off-diagonal interaction parameters.
    ///
    /// Convention: `tau_12 = a12 + b12_k / T`, `tau_21 = a21 + b21_k / T`,
    /// `alpha_12 = alpha_21 = alpha12`.
    #[allow(clippy::too_many_arguments)]
    pub fn binary(
        species: [SpeciesId; 2],
        a12: f64,
        a21: f64,
        b12_k: f64,
        b21_k: f64,
        alpha12: f64,
        minimum_temperature: Temperature,
        maximum_temperature: Temperature,
    ) -> Result<Self, ExcessGibbsErrorV1> {
        Self::new(
            species.into(),
            vec![0.0, a12, a21, 0.0],
            vec![0.0, b12_k, b21_k, 0.0],
            vec![0.0, alpha12, alpha12, 0.0],
            minimum_temperature,
            maximum_temperature,
        )
    }

    /// Load a validated NRTL model from an immutable `fs-matdb` constitutive card.
    ///
    /// Expected card schema:
    /// - `law = "nrtl-excess-gibbs"`, `law_version = 1`, stateless (`state_schema_version = 0`)
    /// - validity axis `"T"` in kelvin
    /// - for every ordered pair `i != j` (`0..species.len()`):
    ///   - `a_{i}_{j}` (dimensionless)
    ///   - `b_{i}_{j}` (temperature dimensions, kelvin)
    /// - for every unordered pair `i < j`:
    ///   - `alpha_{i}_{j}` (dimensionless)
    pub fn from_card(
        card: &ConstitutiveModelCard,
        species: Vec<SpeciesId>,
    ) -> Result<Self, ExcessGibbsErrorV1> {
        card.validate()
            .map_err(|_| ExcessGibbsErrorV1::InvalidModel("NRTL matdb card failed validation"))?;
        if card.law.0 != NRTL_LAW_ID_V1 {
            return Err(ExcessGibbsErrorV1::InvalidModel(
                "NRTL matdb card law id mismatch",
            ));
        }
        if card.law_version != NRTL_LAW_VERSION_V1 {
            return Err(ExcessGibbsErrorV1::InvalidModel(
                "NRTL matdb card law version mismatch",
            ));
        }
        if card.state_schema_version != 0
            || card.initial_state != InitialStatePolicy::ZeroInternalState
        {
            return Err(ExcessGibbsErrorV1::InvalidModel(
                "NRTL matdb card must be stateless",
            ));
        }
        let (t_min, t_max) = card
            .validity
            .bound("T")
            .ok_or(ExcessGibbsErrorV1::InvalidModel(
                "NRTL matdb card missing T validity axis",
            ))?;
        let n = species.len();
        if !(1..=MAX_EXCESS_GIBBS_SPECIES_V1).contains(&n) {
            return Err(ExcessGibbsErrorV1::InvalidModel(
                "NRTL species count must be in 1..=64",
            ));
        }
        let mut a_matrix = vec![0.0; n * n];
        let mut b_matrix_k = vec![0.0; n * n];
        let mut alpha_matrix = vec![0.0; n * n];
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
                            "NRTL matdb card missing a_i_j parameter",
                        ))?;
                if a_param.dims != Dims::NONE {
                    return Err(ExcessGibbsErrorV1::InvalidModel(
                        "NRTL matdb card parameter a_i_j must be dimensionless",
                    ));
                }
                let b_param =
                    card.parameters
                        .get(&b_key)
                        .ok_or(ExcessGibbsErrorV1::InvalidModel(
                            "NRTL matdb card missing b_i_j parameter",
                        ))?;
                if b_param.dims != Temperature::DIMS {
                    return Err(ExcessGibbsErrorV1::InvalidModel(
                        "NRTL matdb card parameter b_i_j must have temperature dimensions [K]",
                    ));
                }
                a_matrix[i * n + j] = a_param.value;
                b_matrix_k[i * n + j] = b_param.value;
                if i < j {
                    let alpha_key = format!("alpha_{i}_{j}");
                    let alpha_param =
                        card.parameters
                            .get(&alpha_key)
                            .ok_or(ExcessGibbsErrorV1::InvalidModel(
                                "NRTL matdb card missing alpha_i_j parameter",
                            ))?;
                    if alpha_param.dims != Dims::NONE {
                        return Err(ExcessGibbsErrorV1::InvalidModel(
                            "NRTL matdb card parameter alpha_i_j must be dimensionless",
                        ));
                    }
                    alpha_matrix[i * n + j] = alpha_param.value;
                    alpha_matrix[j * n + i] = alpha_param.value;
                }
            }
        }
        Self::new(
            species,
            a_matrix,
            b_matrix_k,
            alpha_matrix,
            Temperature::new(t_min),
            Temperature::new(t_max),
        )
    }

    /// Declared species identities in matrix row/column order.
    #[must_use]
    pub fn species(&self) -> &[SpeciesId] {
        &self.species
    }

    /// Inclusive valid temperature interval `[T_min, T_max]` [K].
    #[must_use]
    pub const fn temperature_domain(&self) -> [Temperature; 2] {
        [self.temperature_min, self.temperature_max]
    }

    /// Row-major dimensionless `a_ij` matrix.
    #[must_use]
    pub fn a_matrix(&self) -> &[f64] {
        &self.a_matrix
    }

    /// Row-major `b_ij` [K] matrix.
    #[must_use]
    pub fn b_matrix_k(&self) -> &[f64] {
        &self.b_matrix_k
    }

    /// Row-major symmetric non-randomness `alpha_ij` matrix.
    #[must_use]
    pub fn alpha_matrix(&self) -> &[f64] {
        &self.alpha_matrix
    }

    /// Evaluate exact NRTL excess Gibbs properties, activity coefficients,
    /// analytical temperature derivatives, excess enthalpy, and analytical
    /// composition Jacobian at `(temperature, mole_fractions)`.
    pub fn evaluate(
        &self,
        temperature: Temperature,
        mole_fractions: &[f64],
    ) -> Result<ExcessGibbsEvaluationV1, ExcessGibbsErrorV1> {
        let t = check_temperature_in_domain(
            temperature,
            self.temperature_min.value(),
            self.temperature_max.value(),
        )?;
        let n = self.species.len();
        let x = validate_and_normalize_mole_fractions(mole_fractions, n)?;

        let inv_t = 1.0 / t;
        let inv_t2 = inv_t * inv_t;

        // Precompute tau_ij, G_ij, and -T^2 * d(G_ij)/dT = -alpha_ij * b_ij * G_ij.
        let mut tau = vec![0.0; n * n];
        let mut g = vec![0.0; n * n];
        let mut g_tilde = vec![0.0; n * n];
        let mut n_tilde_coeff = vec![0.0; n * n];

        for i in 0..n {
            for j in 0..n {
                let idx = i * n + j;
                if i == j {
                    tau[idx] = 0.0;
                    g[idx] = 1.0;
                    g_tilde[idx] = 0.0;
                    n_tilde_coeff[idx] = 0.0;
                } else {
                    let a_ij = self.a_matrix[idx];
                    let b_ij = self.b_matrix_k[idx];
                    let alpha_ij = self.alpha_matrix[idx];
                    let tau_ij = check_finite(a_ij + b_ij * inv_t, "NRTL tau_ij")?;
                    let arg = check_finite(-alpha_ij * tau_ij, "NRTL exp argument")?;
                    let g_ij = exp(arg);
                    if !g_ij.is_finite() || g_ij <= 0.0 {
                        return Err(ExcessGibbsErrorV1::Unrepresentable(
                            "NRTL Boltzmann factor G_ij underflowed or overflowed",
                        ));
                    }
                    tau[idx] = tau_ij;
                    g[idx] = g_ij;
                    // d(tau_ij)/dT = -b_ij / T^2 => -T^2 * d(tau_ij)/dT = b_ij
                    // d(G_ij)/dT = -alpha_ij * G_ij * d(tau_ij)/dT = (alpha_ij * b_ij / T^2) * G_ij
                    // => -T^2 * d(G_ij)/dT = -alpha_ij * b_ij * G_ij
                    g_tilde[idx] = check_finite(-alpha_ij * b_ij * g_ij, "NRTL dG_ij/dT")?;
                    // -T^2 * d(tau_ij * G_ij)/dT = b_ij * G_ij * (1 - alpha_ij * tau_ij)
                    n_tilde_coeff[idx] = check_finite(
                        b_ij * g_ij * (1.0 - alpha_ij * tau_ij),
                        "NRTL d(tau_ij G_ij)/dT",
                    )?;
                }
            }
        }

        // Per-component column sums S_k, R_k = N_k / S_k and temperature-scaled tilde sums.
        let mut s_col = vec![0.0; n];
        let mut r_col = vec![0.0; n];
        let mut s_tilde = vec![0.0; n];
        let mut r_tilde = vec![0.0; n];

        for k in 0..n {
            let mut sk = 0.0;
            let mut nk = 0.0;
            let mut stk = 0.0;
            let mut ntk = 0.0;
            for j in 0..n {
                let idx = j * n + k;
                let xj = x[j];
                sk += g[idx] * xj;
                nk += tau[idx] * g[idx] * xj;
                stk += g_tilde[idx] * xj;
                ntk += n_tilde_coeff[idx] * xj;
            }
            if !sk.is_finite() || sk <= 0.0 {
                return Err(ExcessGibbsErrorV1::Unrepresentable(
                    "NRTL denominator S_k is non-positive or non-finite",
                ));
            }
            let rk = check_finite(nk / sk, "NRTL R_k")?;
            let rtk = check_finite((ntk - rk * stk) / sk, "NRTL -T^2 dR_k/dT")?;
            s_col[k] = sk;
            r_col[k] = rk;
            s_tilde[k] = stk;
            r_tilde[k] = rtk;
        }

        // Excess Gibbs energy, excess enthalpy, and excess entropy.
        let mut g_ex_over_rt = 0.0;
        let mut h_ex_over_r = 0.0;
        for i in 0..n {
            g_ex_over_rt += x[i] * r_col[i];
            h_ex_over_r += x[i] * r_tilde[i];
        }
        let g_ex_over_rt = check_finite(g_ex_over_rt, "NRTL g^E/(RT)")?;
        let g_ex = check_finite(g_ex_over_rt * R * t, "NRTL g^E")?;
        let h_ex = check_finite(h_ex_over_r * R, "NRTL h^E")?;
        let s_ex = check_finite(R * (h_ex_over_r * inv_t - g_ex_over_rt), "NRTL s^E")?;

        // Activity coefficients ln gamma_i, gamma_i, d(ln gamma_i)/dT, and partial h_i^E.
        let mut ln_gamma = vec![0.0; n];
        let mut gamma = vec![0.0; n];
        let mut d_ln_gamma_dt = vec![0.0; n];
        let mut partial_h_ex = vec![0.0; n];

        for i in 0..n {
            let mut sum_ln_g = 0.0;
            let mut sum_tilde_d_ln_g = 0.0;
            for k in 0..n {
                let idx = i * n + k;
                let sk = s_col[k];
                let xk_over_sk = x[k] / sk;
                let diff_tau = tau[idx] - r_col[k];
                sum_ln_g += xk_over_sk * g[idx] * diff_tau;

                // -T^2 * d[ (x_k / S_k) * G_ik * (tau_ik - R_k) ] / dT
                // = (x_k / S_k) * [ (g_tilde_ik - g_ik * s_tilde_k / S_k) * (tau_ik - R_k)
                //                   + g_ik * (b_ik - r_tilde_k) ]
                let term_tilde = (g_tilde[idx] - g[idx] * (s_tilde[k] / sk)) * diff_tau
                    + g[idx] * (self.b_matrix_k[idx] - r_tilde[k]);
                sum_tilde_d_ln_g += xk_over_sk * term_tilde;
            }
            let lng = check_finite(r_col[i] + sum_ln_g, "NRTL ln gamma_i")?;
            let gam = exp(lng);
            if !gam.is_finite() || gam <= 0.0 {
                return Err(ExcessGibbsErrorV1::Unrepresentable(
                    "NRTL activity coefficient gamma_i overflowed or underflowed",
                ));
            }
            // neg_t2_dlng_dt = -T^2 * d(ln gamma_i)/dT
            let neg_t2_dlng_dt = check_finite(r_tilde[i] + sum_tilde_d_ln_g, "NRTL h_i^E / R")?;
            let dlng_dt = check_finite(-neg_t2_dlng_dt * inv_t2, "NRTL d(ln gamma_i)/dT")?;
            let hi_ex = check_finite(R * neg_t2_dlng_dt, "NRTL partial excess enthalpy")?;

            ln_gamma[i] = lng;
            gamma[i] = gam;
            d_ln_gamma_dt[i] = dlng_dt;
            partial_h_ex[i] = hi_ex;
        }

        // Exact symmetric homogeneous-degree-0 Jacobian:
        // d(ln gamma_i)/dx_j = (G_ji / S_i)*(tau_ji - R_i)
        //                    + (G_ij / S_j)*(tau_ij - R_j)
        //                    - sum_k (x_k * G_ik * G_jk / S_k^2) * (tau_ik + tau_jk - 2*R_k)
        let mut d_ln_gamma_dx = vec![0.0; n * n];
        for i in 0..n {
            for j in i..n {
                let idx_ji = j * n + i;
                let idx_ij = i * n + j;
                let term1 = (g[idx_ji] / s_col[i]) * (tau[idx_ji] - r_col[i]);
                let term2 = (g[idx_ij] / s_col[j]) * (tau[idx_ij] - r_col[j]);
                let mut term3 = 0.0;
                for k in 0..n {
                    let idx_ik = i * n + k;
                    let idx_jk = j * n + k;
                    let sk = s_col[k];
                    let weight = x[k] * (g[idx_ik] / sk) * (g[idx_jk] / sk);
                    term3 += weight * (tau[idx_ik] + tau[idx_jk] - 2.0 * r_col[k]);
                }
                let entry = check_finite(term1 + term2 - term3, "NRTL d(ln gamma_i)/dx_j")?;
                d_ln_gamma_dx[i * n + j] = entry;
                d_ln_gamma_dx[j * n + i] = entry;
            }
        }

        Ok(ExcessGibbsEvaluationV1 {
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
        })
    }

    /// Pre-validated Ethanol(0) - Water(1) NRTL parameter set (Aspen Plus APV88/VLE-IG).
    ///
    /// Parameters: `a12 = -0.8009`, `a21 = 3.4578`, `b12 = 246.18 K`,
    /// `b21 = -586.0809 K`, `alpha12 = 0.3`, valid on `273.15..=423.15 K`.
    #[must_use]
    pub fn ethanol_water_aspen_v1() -> Self {
        let species = [
            SpeciesId::new("C2H5OH").expect("valid literal species id"),
            SpeciesId::new("H2O").expect("valid literal species id"),
        ];
        Self::binary(
            species,
            -0.8009,
            3.4578,
            246.18,
            -586.0809,
            0.3,
            Temperature::new(273.15),
            Temperature::new(423.15),
        )
        .expect("published Aspen Plus ethanol-water NRTL parameters are valid")
    }

    /// Pre-validated Ethanol(0) - Water(1) Modified Raoult VLE model using Aspen Plus
    /// NRTL paired with NIST ethanol Antoine and DECHEMA water Antoine vapor pressures.
    ///
    /// Reproduces the 1 atm minimum-boiling azeotrope (`T_az ≈ 351.30 K`, `x_EtOH ≈ 0.897`)
    /// and generalizes to arbitrary barometric pressures (e.g., Denver 83.5 kPa).
    #[must_use]
    pub fn ethanol_water_vle_v1() -> ModifiedRaoultVleV1<Self> {
        ModifiedRaoultVleV1::new(
            Self::ethanol_water_aspen_v1(),
            vec![
                VaporPressureLawV1::ethanol_antoine_nist(),
                VaporPressureLawV1::water_antoine_dechema(),
            ],
        )
        .expect("published ethanol-water VLE model is valid")
    }

    /// Pre-validated Water(0) - Sucrose(1) excess-Gibbs NRTL set reproducing the
    /// Norrish (1966) sucrose water-activity correlation (`ln gamma_w = -6.47 * x_s^2`).
    ///
    /// Realized identically by setting `alpha = 0` and `a12 = a21 = -3.235`
    /// (`b12 = b21 = 0`), where NRTL reduces exactly to the two-suffix Margules form.
    #[must_use]
    pub fn water_sucrose_norrish_v1() -> Self {
        let species = [
            SpeciesId::new("H2O").expect("valid literal species id"),
            SpeciesId::new("C12H22O11").expect("valid literal species id"),
        ];
        Self::binary(
            species,
            -3.235,
            -3.235,
            0.0,
            0.0,
            0.0,
            Temperature::new(273.16),
            Temperature::new(433.15),
        )
        .expect("published Norrish sucrose-water parameters are valid")
    }

    /// Pre-validated Water(0) - Sucrose(1) Modified Raoult boiling-point elevation VLE model
    /// (`IAPWS` water vapor pressure + non-volatile sucrose).
    #[must_use]
    pub fn water_sucrose_vle_v1() -> ModifiedRaoultVleV1<Self> {
        ModifiedRaoultVleV1::new(
            Self::water_sucrose_norrish_v1(),
            vec![
                VaporPressureLawV1::iapws_water(),
                VaporPressureLawV1::non_volatile(),
            ],
        )
        .expect("water-sucrose VLE model is valid")
    }

    /// Pre-validated Water(0) - NaCl(1) effective NRTL model on a dissociated-ion mole
    /// fraction basis (`nu = 2` moles of ions per mole NaCl, effective molar mass
    /// `29.22 g/mol`, `a12 = a21 = -1.65`, `alpha = 0`).
    ///
    /// Note: this is a neutral-potential effective correlation for culinary brine
    /// boiling-point elevation (0-26 wt% NaCl), not an ionic Debye-Huckel/Pitzer
    /// electrolyte model.
    #[must_use]
    pub fn water_nacl_ion_basis_v1() -> Self {
        let species = [
            SpeciesId::new("H2O").expect("valid literal species id"),
            SpeciesId::new("NaCl_ion").expect("valid literal species id"),
        ];
        Self::binary(
            species,
            -1.65,
            -1.65,
            0.0,
            0.0,
            0.0,
            Temperature::new(273.16),
            Temperature::new(423.15),
        )
        .expect("effective brine NRTL parameters are valid")
    }

    /// Pre-validated Water(0) - NaCl(1) Modified Raoult brine boiling-point elevation
    /// model (`IAPWS` water vapor pressure + non-volatile salt on dissociated-ion basis).
    #[must_use]
    pub fn water_nacl_vle_v1() -> ModifiedRaoultVleV1<Self> {
        ModifiedRaoultVleV1::new(
            Self::water_nacl_ion_basis_v1(),
            vec![
                VaporPressureLawV1::iapws_water(),
                VaporPressureLawV1::non_volatile(),
            ],
        )
        .expect("water-NaCl VLE model is valid")
    }
}

impl ExcessGibbsModelV1 for NrtlModelV1 {
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
    fn nrtl_pure_limits_and_infinite_dilution_match_exact_formulas() {
        let model = NrtlModelV1::ethanol_water_aspen_v1();
        let t = Temperature::new(351.30);

        let pure0 = model.evaluate(t, &[1.0, 0.0]).unwrap();
        assert!(pure0.reduced_excess_gibbs().abs() < 1e-15);
        assert!(pure0.excess_enthalpy_j_per_mol().abs() < 1e-12);
        assert!((pure0.activity_coefficients()[0] - 1.0).abs() < 1e-15);

        let tau12 = -0.8009 + 246.18 / t.value();
        let tau21 = 3.4578 - 586.0809 / t.value();
        let g12 = exp(-0.3 * tau12);
        let g21 = exp(-0.3 * tau21);
        let expected_ln_gamma2_inf = tau12 + tau21 * g21;
        assert!((pure0.ln_activity_coefficients()[1] - expected_ln_gamma2_inf).abs() < 1e-14);

        let pure1 = model.evaluate(t, &[0.0, 1.0]).unwrap();
        assert!((pure1.activity_coefficients()[1] - 1.0).abs() < 1e-15);
        let expected_ln_gamma1_inf = tau21 + tau12 * g12;
        assert!((pure1.ln_activity_coefficients()[0] - expected_ln_gamma1_inf).abs() < 1e-14);
    }

    #[test]
    fn nrtl_gibbs_duhem_and_finite_difference_derivatives_match_ternary() {
        let sp = vec![
            SpeciesId::new("C2H5OH").unwrap(),
            SpeciesId::new("H2O").unwrap(),
            SpeciesId::new("CH3COOH").unwrap(),
        ];
        let a = vec![0.0, -0.8009, 0.25, 3.4578, 0.0, -0.40, -0.15, 0.60, 0.0];
        let b = vec![
            0.0, 246.18, -120.0, -586.0809, 0.0, 310.0, 180.0, -210.0, 0.0,
        ];
        let alpha = vec![0.0, 0.3, 0.25, 0.3, 0.0, 0.35, 0.25, 0.35, 0.0];
        let model = NrtlModelV1::new(
            sp,
            a,
            b,
            alpha,
            Temperature::new(273.15),
            Temperature::new(450.0),
        )
        .unwrap();

        let t = Temperature::new(340.0);
        let x = [0.27, 0.51, 0.22];
        let eval = model.evaluate(t, &x).unwrap();

        // Exact Gibbs-Duhem column residuals < 1e-13.
        for (j, res) in eval.gibbs_duhem_column_residuals().into_iter().enumerate() {
            assert!(res.abs() < 1e-13, "col {j} Gibbs-Duhem residual {res:e}");
        }

        // Analytical vs central finite difference d(ln gamma_i)/dT and h^E.
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

        // Directional composition derivative along simplex tangent dx = [1, -0.6, -0.4].
        let dir = [1.0, -0.6, -0.4];
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
