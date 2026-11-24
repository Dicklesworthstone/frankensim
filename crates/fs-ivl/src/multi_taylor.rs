//! Multivariate Taylor models (plan §6.4, bead frankensim-epic-bedrock-6ys.23.2.1).
//!
//! Provides canonical bounded multivariate functional enclosures:
//! `f(x) ∈ P(x − c) + remainder` over a d-dimensional domain box `D = I₁ × … × I_d`.
//!
//! Features:
//! - Graded lexicographical (grlex) multi-index ordering.
//! - Canonical variable identity, physical units, and domain box definitions.
//! - Checked dimension, order, term count, and memory budgets.
//! - Zero-width axis reduction ensuring normalization never divides by zero.
//! - Lossless bijection between 1D TaylorModel1 and d=1 multivariate models.

#![deny(unsafe_code)]

use crate::Interval;
use crate::taylor::{TaylorModel1, TaylorModelError};
use core::cmp::Ordering;

/// Maximum admitted multivariate dimension.
pub const MAX_MULTIVARIATE_DIM: usize = 32;

/// Maximum admitted multivariate order.
pub const MAX_MULTIVARIATE_ORDER: usize = 30;

/// Maximum admitted total polynomial terms across all multi-indices.
pub const MAX_MULTIVARIATE_TERMS: usize = 100_000;

/// Maximum memory allocation per model in bytes (16 MiB).
pub const MAX_MODEL_MEMORY_BYTES: usize = 16 * 1024 * 1024;

/// Specification of a single variable in a multivariate Taylor model domain.
#[derive(Debug, Clone, PartialEq)]
pub struct VariableInfo {
    /// Identifier name of the variable.
    pub name: String,
    /// Domain interval of valid inputs.
    pub domain: Interval,
    /// Optional physical unit (e.g., "m", "s", "rad").
    pub unit: Option<String>,
}

impl VariableInfo {
    /// Create a new variable specification.
    pub fn new(name: impl Into<String>, domain: Interval) -> Result<Self, TaylorModelError> {
        if !domain.lo().is_finite() || !domain.hi().is_finite() {
            return Err(TaylorModelError::NonFiniteDomain);
        }
        if domain.lo() > domain.hi() {
            return Err(TaylorModelError::NonFiniteDomain);
        }
        Ok(Self {
            name: name.into(),
            domain,
            unit: None,
        })
    }

    /// Attach a physical unit.
    #[must_use]
    pub fn with_unit(mut self, unit: impl Into<String>) -> Self {
        self.unit = Some(unit.into());
        self
    }

    /// Center (midpoint) of the domain.
    #[must_use]
    pub fn center(&self) -> f64 {
        self.domain.midpoint()
    }

    /// Half-width (radius) of the domain.
    #[must_use]
    pub fn radius(&self) -> f64 {
        self.domain.width() * 0.5
    }

    /// Whether this axis is degenerate (zero-width).
    #[must_use]
    pub fn is_fixed(&self) -> bool {
        // Degeneracy is exact; adjacent endpoints must retain their axis.
        self.domain.lo().partial_cmp(&self.domain.hi()) == Some(Ordering::Equal)
    }
}

/// Compute binomial coefficient `C(n, k) = n! / (k! (n - k)!)` in 64-bit unsigned arithmetic.
#[must_use]
pub fn binomial(n: usize, k: usize) -> Option<usize> {
    if k > n {
        return Some(0);
    }
    if k == 0 || k == n {
        return Some(1);
    }
    let k = k.min(n - k);
    let mut res: usize = 1;
    for i in 1..=k {
        let num = n - k + i;
        // Check overflow
        res = res.checked_mul(num)? / i;
    }
    Some(res)
}

/// Total number of terms in a d-variable polynomial up to total degree `order`: `binom(d + p, p)`.
pub fn term_count(dim: usize, order: usize) -> Result<usize, TaylorModelError> {
    if dim == 0 {
        return Ok(1); // Single constant term
    }
    if dim > MAX_MULTIVARIATE_DIM {
        return Err(TaylorModelError::OrderTooLarge {
            requested: dim,
            maximum: MAX_MULTIVARIATE_DIM,
        });
    }
    if order > MAX_MULTIVARIATE_ORDER {
        return Err(TaylorModelError::OrderTooLarge {
            requested: order,
            maximum: MAX_MULTIVARIATE_ORDER,
        });
    }
    let terms = binomial(dim + order, order).ok_or(TaylorModelError::AllocationFailed {
        coefficients: usize::MAX,
    })?;
    if terms > MAX_MULTIVARIATE_TERMS {
        return Err(TaylorModelError::AllocationFailed {
            coefficients: terms,
        });
    }
    let bytes = terms.saturating_mul(core::mem::size_of::<f64>());
    if bytes > MAX_MODEL_MEMORY_BYTES {
        return Err(TaylorModelError::AllocationFailed {
            coefficients: terms,
        });
    }
    Ok(terms)
}

/// Graded multi-index generator: generates all multi-indices of dimension `dim` with total degree `<= order`
/// in graded lexicographical (grlex) order.
#[must_use]
pub fn generate_multi_indices(dim: usize, order: usize) -> Vec<Vec<u16>> {
    if dim == 0 {
        return vec![vec![]];
    }
    let mut result = Vec::new();
    for total_deg in 0..=order {
        generate_fixed_degree_indices(dim, total_deg, &mut result);
    }
    result
}

fn generate_fixed_degree_indices(dim: usize, target_deg: usize, out: &mut Vec<Vec<u16>>) {
    fn backtrack(
        dim: usize,
        idx: usize,
        remaining_deg: usize,
        current: &mut Vec<u16>,
        out: &mut Vec<Vec<u16>>,
    ) {
        if idx == dim - 1 {
            current[idx] = remaining_deg as u16;
            out.push(current.clone());
            return;
        }
        for d in (0..=remaining_deg).rev() {
            current[idx] = d as u16;
            backtrack(dim, idx + 1, remaining_deg - d, current, out);
        }
    }
    let mut current = vec![0u16; dim];
    backtrack(dim, 0, target_deg, &mut current, out);
}

/// Canonical bounded multivariate Taylor model.
///
/// Represents an enclosure `P(x - c) + remainder` where `P` is a polynomial in `(x_0 - c_0, ..., x_{d-1} - c_{d-1})`
/// with coefficients ordered by graded lexicographical multi-indices.
#[derive(Debug, Clone, PartialEq)]
pub struct TaylorModel {
    /// Dimension (number of variables).
    dim: usize,
    /// Maximum polynomial total degree.
    order: usize,
    /// Variable metadata (names, physical units, domain intervals).
    variables: Vec<VariableInfo>,
    /// Expansion centers `c_i = midpoint(I_i)`.
    centers: Vec<f64>,
    /// Polynomial coefficients stored in graded multi-index order.
    coefficients: Vec<f64>,
    /// Multi-index mapping table (cached).
    multi_indices: Vec<Vec<u16>>,
    /// Rigorous interval remainder enclosing all truncation, rounding, and composition errors.
    remainder: Interval,
}

impl TaylorModel {
    /// Create an identity variable model for variable `var_idx` among `variables`.
    pub fn variable(
        var_idx: usize,
        variables: Vec<VariableInfo>,
        order: usize,
    ) -> Result<Self, TaylorModelError> {
        let dim = variables.len();
        if dim == 0 || var_idx >= dim {
            return Err(TaylorModelError::IncompatibleModels);
        }
        if order < 1 {
            return Err(TaylorModelError::VariableOrderTooSmall {
                requested: order,
                minimum: 1,
            });
        }
        let terms = term_count(dim, order)?;
        let multi_indices = generate_multi_indices(dim, order);
        let mut coefficients = vec![0.0; terms];

        let centers: Vec<f64> = variables.iter().map(VariableInfo::center).collect();

        // Constant term is the center coordinate c_var
        coefficients[0] = centers[var_idx];

        // Linear term for var_idx (multi-index with 1 at var_idx and 0 elsewhere)
        for (i, mi) in multi_indices.iter().enumerate() {
            if mi.iter().sum::<u16>() == 1 && mi[var_idx] == 1 {
                coefficients[i] = 1.0;
                break;
            }
        }

        Ok(Self {
            dim,
            order,
            variables,
            centers,
            coefficients,
            multi_indices,
            remainder: Interval::point(0.0),
        })
    }

    /// Create a constant Taylor model on `variables`.
    pub fn constant(
        val: f64,
        variables: Vec<VariableInfo>,
        order: usize,
    ) -> Result<Self, TaylorModelError> {
        if !val.is_finite() {
            return Err(TaylorModelError::NonFiniteConstant);
        }
        let dim = variables.len();
        let terms = term_count(dim, order)?;
        let multi_indices = generate_multi_indices(dim, order);
        let mut coefficients = vec![0.0; terms];
        coefficients[0] = val;
        let centers = variables.iter().map(VariableInfo::center).collect();

        Ok(Self {
            dim,
            order,
            variables,
            centers,
            coefficients,
            multi_indices,
            remainder: Interval::point(0.0),
        })
    }

    /// Construct from an existing 1D [`TaylorModel1`].
    pub fn from_tm1(
        tm1: &TaylorModel1,
        var_name: impl Into<String>,
        unit: Option<String>,
    ) -> Result<Self, TaylorModelError> {
        let mut var_info = VariableInfo::new(var_name, tm1.domain())?;
        var_info.unit = unit;
        let variables = vec![var_info];
        let order = tm1.order();
        let dim = 1;
        let terms = term_count(dim, order)?;
        let multi_indices = generate_multi_indices(dim, order);

        let mut coefficients = vec![0.0; terms];
        for (k, &c) in tm1.poly().iter().enumerate().take(terms) {
            coefficients[k] = c;
        }

        Ok(Self {
            dim: 1,
            order,
            variables,
            centers: vec![tm1.center()],
            coefficients,
            multi_indices,
            remainder: tm1.remainder(),
        })
    }

    /// Convert a 1D multivariate model back to [`TaylorModel1`].
    pub fn to_tm1(&self) -> Result<TaylorModel1, TaylorModelError> {
        if self.dim != 1 {
            return Err(TaylorModelError::IncompatibleModels);
        }
        let domain = self.variables[0].domain;
        let mut tm1 = if self.coefficients.len() > 1
            && self.coefficients[1].partial_cmp(&1.0) == Some(Ordering::Equal)
        {
            TaylorModel1::variable(domain, self.order)?
        } else {
            TaylorModel1::constant(self.coefficients[0], domain, self.order)?
        };
        // Copy coefficients and remainder
        for (k, &c) in self.coefficients.iter().enumerate() {
            if k < tm1.poly().len() {
                tm1.set_poly_coeff(k, c);
            }
        }
        tm1.set_remainder(self.remainder);
        Ok(tm1)
    }

    /// Dimension (number of variables).
    #[must_use]
    pub const fn dim(&self) -> usize {
        self.dim
    }

    /// Polynomial order (maximum degree).
    #[must_use]
    pub const fn order(&self) -> usize {
        self.order
    }

    /// Number of polynomial terms.
    #[must_use]
    pub fn term_count(&self) -> usize {
        self.coefficients.len()
    }

    /// Variable metadata slice.
    #[must_use]
    pub fn variables(&self) -> &[VariableInfo] {
        &self.variables
    }

    /// Domain bounding box.
    #[must_use]
    pub fn domain_box(&self) -> Vec<Interval> {
        self.variables.iter().map(|v| v.domain).collect()
    }

    /// Expansion centers.
    #[must_use]
    pub fn centers(&self) -> &[f64] {
        &self.centers
    }

    /// Polynomial coefficients in graded multi-index order.
    #[must_use]
    pub fn coefficients(&self) -> &[f64] {
        &self.coefficients
    }

    /// Multi-index mapping table.
    #[must_use]
    pub fn multi_indices(&self) -> &[Vec<u16>] {
        &self.multi_indices
    }

    /// Rigorous interval remainder.
    #[must_use]
    pub const fn remainder(&self) -> Interval {
        self.remainder
    }

    /// Add another Taylor model with compatible domain and variables.
    pub fn add(&self, other: &Self) -> Result<Self, TaylorModelError> {
        self.check_compatible(other)?;
        let mut res = self.clone();
        res.remainder = res.remainder + other.remainder;
        for (index, (c_res, &c_other)) in res
            .coefficients
            .iter_mut()
            .zip(&other.coefficients)
            .enumerate()
        {
            if *c_res == 0.0 || c_other == 0.0 {
                *c_res += c_other;
            } else {
                *c_res = crate::taylor::split_mid(
                    Interval::point(*c_res) + Interval::point(c_other),
                    self.monomial_range(&self.multi_indices[index]),
                    &mut res.remainder,
                );
            }
        }
        Ok(res)
    }

    /// Subtract another Taylor model.
    pub fn sub(&self, other: &Self) -> Result<Self, TaylorModelError> {
        self.check_compatible(other)?;
        let mut res = self.clone();
        res.remainder = res.remainder - other.remainder;
        for (index, (c_res, &c_other)) in res
            .coefficients
            .iter_mut()
            .zip(&other.coefficients)
            .enumerate()
        {
            if *c_res == 0.0 || c_other == 0.0 {
                *c_res -= c_other;
            } else {
                *c_res = crate::taylor::split_mid(
                    Interval::point(*c_res) - Interval::point(c_other),
                    self.monomial_range(&self.multi_indices[index]),
                    &mut res.remainder,
                );
            }
        }
        Ok(res)
    }

    /// Multiply by a scalar factor.
    pub fn scale(&self, factor: f64) -> Result<Self, TaylorModelError> {
        if !factor.is_finite() {
            return Err(TaylorModelError::NonFiniteScaleFactor);
        }
        let mut res = self.clone();
        res.remainder = res.remainder * Interval::point(factor);
        for (index, c) in res.coefficients.iter_mut().enumerate() {
            // Exact zero and sign changes need no domain power, which may
            // overflow even when this sparse coefficient stays finite.
            if *c == 0.0 || factor == 0.0 || factor.abs() == 1.0 {
                *c *= factor;
            } else {
                *c = crate::taylor::split_mid(
                    Interval::point(*c) * Interval::point(factor),
                    self.monomial_range(&self.multi_indices[index]),
                    &mut res.remainder,
                );
            }
        }
        Ok(res)
    }

    /// Compute the interval range of a multi-index monomial over the centered domain box.
    #[must_use]
    pub fn monomial_range(&self, mi: &[u16]) -> Interval {
        let mut acc = Interval::point(1.0);
        for (i, &power) in mi.iter().enumerate() {
            if power > 0 {
                // The rounded midpoint need not be equidistant from the
                // endpoints (for adjacent floats it is itself an endpoint).
                let centered = self.variables[i].domain - Interval::point(self.centers[i]);
                acc = acc * centered.powi(power as usize);
            }
        }
        acc
    }

    /// Multiply two Taylor models with rigorous remainder truncation.
    pub fn mul(&self, other: &Self) -> Result<Self, TaylorModelError> {
        self.check_compatible(other)?;

        let mut retained = vec![Interval::point(0.0); self.coefficients.len()];
        let mut tail_rem = Interval::point(0.0);

        for (i, &a) in self.coefficients.iter().enumerate() {
            if a == 0.0 {
                continue;
            }
            let mi_a = &self.multi_indices[i];
            let deg_a: usize = mi_a.iter().map(|&d| d as usize).sum();

            for (j, &b) in other.coefficients.iter().enumerate() {
                if b == 0.0 {
                    continue;
                }
                let mi_b = &other.multi_indices[j];
                let deg_b: usize = mi_b.iter().map(|&d| d as usize).sum();

                let mut mi_prod = vec![0u16; self.dim];
                for k in 0..self.dim {
                    mi_prod[k] = mi_a[k] + mi_b[k];
                }
                let prod_ab = Interval::point(a) * Interval::point(b);

                if deg_a + deg_b <= self.order {
                    if let Some(pos) = self.multi_indices.iter().position(|m| m == &mi_prod) {
                        retained[pos] = retained[pos] + prod_ab;
                    }
                } else {
                    let mono = self.monomial_range(&mi_prod);
                    tail_rem = tail_rem + prod_ab * mono;
                }
            }
        }

        // Split each retained interval coefficient back into its polynomial
        // part and the roundoff contribution over the whole centered domain.
        let new_coeffs = retained
            .into_iter()
            .zip(&self.multi_indices)
            .map(|(coefficient, mi)| {
                if coefficient == Interval::point(0.0) {
                    0.0
                } else {
                    crate::taylor::split_mid(coefficient, self.monomial_range(mi), &mut tail_rem)
                }
            })
            .collect();

        // (P + R)(Q + S) = PQ + P*S + Q*R + R*S. Using the full model
        // ranges here would count R*S three times, which can exclude the
        // true product when the remainders have a definite sign.
        let dom_box = self.domain_box();
        let p1_range = self.polynomial_on_box(&dom_box);
        let p2_range = other.polynomial_on_box(&dom_box);
        let total_rem = tail_rem
            + p1_range * other.remainder
            + p2_range * self.remainder
            + self.remainder * other.remainder;

        Ok(Self {
            dim: self.dim,
            order: self.order,
            variables: self.variables.clone(),
            centers: self.centers.clone(),
            coefficients: new_coeffs,
            multi_indices: self.multi_indices.clone(),
            remainder: total_rem,
        })
    }

    /// Compute the reciprocal of the entire model if its range excludes zero.
    pub fn reciprocal(&self) -> Result<Self, TaylorModelError> {
        let r = self.range()?;
        if r.contains_zero() {
            return Err(TaylorModelError::DenominatorContainsZero);
        }
        let c0 = self.coefficients[0];
        let one = Interval::point(1.0);
        let denominator = if c0 == 0.0 { r } else { Interval::point(c0) };
        let mut inverse = Self::constant(0.0, self.variables.clone(), self.order)?;
        inverse.coefficients[0] =
            crate::taylor::split_mid(one / denominator, one, &mut inverse.remainder);
        if c0 == 0.0 {
            // A nonzero remainder can exclude zero even with no polynomial
            // constant. The range reciprocal still supplies an enclosure.
            return Ok(inverse);
        }

        // self = c0 * (1 + u). Retain the rounding of 1/c0 in every use.
        let mut u = self.clone();
        u.coefficients[0] = 0.0;
        u = u.mul(&inverse)?;

        // Expand 1 / (1 + u) = 1 - u + u^2 - u^3 + ... + (-u)^p.
        let mut sum = Self::constant(1.0, self.variables.clone(), self.order)?;
        let mut power = sum.clone();
        let minus_u = u.scale(-1.0)?;
        for _ in 0..self.order {
            power = power.mul(&minus_u)?;
            sum = sum.add(&power)?;
        }

        // The finite geometric identity has tail (-u)^(p+1) / self.
        // Its sign matters for one-sided remainders, including at order zero.
        // Use the admitted denominator range instead of widening 1 + u again.
        let u_box = u.range()?;
        let rem_geom = (-u_box).powi(self.order + 1) / r;

        let mut res = sum.mul(&inverse)?;
        res.remainder = res.remainder + rem_geom;
        Ok(res)
    }

    /// Divide two Taylor models: `self / other`.
    pub fn div(&self, other: &Self) -> Result<Self, TaylorModelError> {
        let recip = other.reciprocal()?;
        self.mul(&recip)
    }

    /// Truncate model polynomial to a smaller degree `new_order <= self.order`.
    pub fn truncate(&self, new_order: usize) -> Result<Self, TaylorModelError> {
        if new_order > self.order {
            return Err(TaylorModelError::TruncationOrderTooLarge {
                requested: new_order,
                current: self.order,
            });
        }
        if new_order == self.order {
            return Ok(self.clone());
        }

        let new_terms = term_count(self.dim, new_order)?;
        let new_multi_indices = generate_multi_indices(self.dim, new_order);
        let mut new_coeffs = vec![0.0; new_terms];
        let mut trunc_rem = Interval::point(0.0);

        for (coeff, mi) in self.coefficients.iter().zip(&self.multi_indices) {
            if *coeff == 0.0 {
                continue;
            }
            let deg: usize = mi.iter().map(|&d| d as usize).sum();
            if deg <= new_order {
                if let Some(pos) = new_multi_indices.iter().position(|m| m == mi) {
                    new_coeffs[pos] = *coeff;
                }
            } else {
                let mono = self.monomial_range(mi);
                trunc_rem = trunc_rem + Interval::point(*coeff) * mono;
            }
        }

        Ok(Self {
            dim: self.dim,
            order: new_order,
            variables: self.variables.clone(),
            centers: self.centers.clone(),
            coefficients: new_coeffs,
            multi_indices: new_multi_indices,
            remainder: self.remainder + trunc_rem,
        })
    }

    /// Range enclosure over the entire domain box.
    pub fn range(&self) -> Result<Interval, TaylorModelError> {
        let dom_box = self.domain_box();
        self.eval_box(&dom_box)
    }

    /// Range enclosure over a sub-box of the domain.
    pub fn range_subdomain(&self, sub_box: &[Interval]) -> Result<Interval, TaylorModelError> {
        self.eval_box(sub_box)
    }

    /// Evaluate the Taylor model over an input interval box `X ⊆ D`.
    pub fn eval_box(&self, box_: &[Interval]) -> Result<Interval, TaylorModelError> {
        if box_.len() != self.dim {
            return Err(TaylorModelError::IncompatibleModels);
        }
        Ok(self.polynomial_on_box(box_) + self.remainder)
    }

    fn polynomial_on_box(&self, box_: &[Interval]) -> Interval {
        // Center deviations (x_i - c_i)
        let mut deltas = Vec::with_capacity(self.dim);
        for (i, &iv) in box_.iter().enumerate() {
            deltas.push(iv - Interval::point(self.centers[i]));
        }

        let mut sum = Interval::point(0.0);
        for (coeff, mi) in self.coefficients.iter().zip(&self.multi_indices) {
            if *coeff == 0.0 {
                continue;
            }
            let mut term = Interval::point(*coeff);
            for (var_i, &power) in mi.iter().enumerate() {
                if power > 0 {
                    term = term * deltas[var_i].powi(power as usize);
                }
            }
            sum = sum + term;
        }

        sum
    }

    /// Recenter the Taylor model to a subdomain box `sub_box ⊆ D`.
    ///
    /// Changes polynomial expansion centers from `c` to `c'` via exact binomial shifts
    /// without any polynomial truncation loss.
    pub fn recenter(&self, sub_box: &[Interval]) -> Result<Self, TaylorModelError> {
        if sub_box.len() != self.dim {
            return Err(TaylorModelError::IncompatibleModels);
        }

        // Validate that sub_box is contained in current domain
        for (&sub_iv, var) in sub_box.iter().zip(&self.variables) {
            if !var.domain.encloses(sub_iv) {
                return Err(TaylorModelError::IncompatibleModels);
            }
            if !sub_iv.lo().is_finite() || !sub_iv.hi().is_finite() {
                return Err(TaylorModelError::NonFiniteDomain);
            }
        }

        let new_centers: Vec<f64> = sub_box.iter().map(|iv| iv.midpoint()).collect();
        let mut shifted_coeffs = self.coefficients.clone();

        // Shift axis by axis
        for k in 0..self.dim {
            let delta = new_centers[k] - self.centers[k];
            if delta == 0.0 {
                continue;
            }

            let mut next_coeffs = vec![0.0; self.coefficients.len()];
            for (idx, &coeff) in shifted_coeffs.iter().enumerate() {
                if coeff == 0.0 {
                    continue;
                }
                let mi = &self.multi_indices[idx];
                let alpha_k = mi[k] as usize;

                for j in 0..=alpha_k {
                    let bin = binomial(alpha_k, j).unwrap_or(1) as f64;
                    let delta_pow = delta.powi((alpha_k - j) as i32);
                    let weight = coeff * bin * delta_pow;

                    let mut beta = mi.clone();
                    beta[k] = j as u16;

                    if let Some(target_pos) = self.multi_indices.iter().position(|m| m == &beta) {
                        next_coeffs[target_pos] += weight;
                    }
                }
            }
            shifted_coeffs = next_coeffs;
        }

        // Update variable domain metadata
        let mut new_vars = self.variables.clone();
        for (var, &new_dom) in new_vars.iter_mut().zip(sub_box) {
            var.domain = new_dom;
        }

        Ok(Self {
            dim: self.dim,
            order: self.order,
            variables: new_vars,
            centers: new_centers,
            coefficients: shifted_coeffs,
            multi_indices: self.multi_indices.clone(),
            remainder: self.remainder,
        })
    }

    /// Subdivide the model along `axis` into two child models: left `[lo, mid]` and right `[mid, hi]`.
    pub fn subdivide_axis(&self, axis: usize) -> Result<(Self, Self), TaylorModelError> {
        if axis >= self.dim {
            return Err(TaylorModelError::IncompatibleModels);
        }
        let dom = self.variables[axis].domain;
        let mid = dom.midpoint();

        let mut left_box = self.domain_box();
        left_box[axis] = Interval::new(dom.lo(), mid);

        let mut right_box = self.domain_box();
        right_box[axis] = Interval::new(mid, dom.hi());

        let left_model = self.recenter(&left_box)?;
        let right_model = self.recenter(&right_box)?;

        Ok((left_model, right_model))
    }

    /// Subdivide all non-fixed axes, producing up to `2^d` child models.
    pub fn subdivide_all_axes(&self) -> Result<Vec<Self>, TaylorModelError> {
        let mut models = vec![self.clone()];
        for axis in 0..self.dim {
            if self.variables[axis].is_fixed() {
                continue;
            }
            let mut next_generation = Vec::with_capacity(models.len() * 2);
            for m in &models {
                let (left, right) = m.subdivide_axis(axis)?;
                next_generation.push(left);
                next_generation.push(right);
            }
            models = next_generation;
        }
        Ok(models)
    }

    /// Multivariate function composition `self(inner[0], ..., inner[d-1])`.
    pub fn compose(&self, inner: &[TaylorModel]) -> Result<Self, TaylorModelError> {
        if inner.len() != self.dim {
            return Err(TaylorModelError::IncompatibleModels);
        }
        let input_vars = inner[0].variables.clone();
        let target_order = inner[0].order;

        // Verify all inner models share compatible domain and order
        for model in inner {
            if model.variables != input_vars || model.order != target_order {
                return Err(TaylorModelError::IncompatibleModels);
            }
        }

        // Deviations (g_i - c_i)
        let mut h = Vec::with_capacity(self.dim);
        for (i, model) in inner.iter().enumerate() {
            let c_i = self.centers[i];
            let c_model = Self::constant(c_i, input_vars.clone(), target_order)?;
            h.push(model.sub(&c_model)?);
        }

        let mut result = Self::constant(self.coefficients[0], input_vars.clone(), target_order)?;

        for (idx, &coeff) in self.coefficients.iter().enumerate() {
            if coeff == 0.0 || idx == 0 {
                continue;
            }
            let mi = &self.multi_indices[idx];
            let mut monomial = Self::constant(1.0, input_vars.clone(), target_order)?;

            for (var_i, &power) in mi.iter().enumerate() {
                for _ in 0..power {
                    monomial = monomial.mul(&h[var_i])?;
                }
            }

            let term = monomial.scale(coeff)?;
            result = result.add(&term)?;
        }

        result.remainder = result.remainder + self.remainder;
        Ok(result)
    }

    /// Partial derivative with respect to `axis`.
    ///
    /// Differentiates the polynomial with outward coefficient rounding. A
    /// nonconstant value remainder provides no derivative bound, so the
    /// derivative remainder is the whole line in that case.
    /// Differentiating along a fixed/degenerate zero-width axis returns exact zero.
    pub fn diff(&self, axis: usize) -> Result<Self, TaylorModelError> {
        if axis >= self.dim {
            return Err(TaylorModelError::AxisIndexOutOfBounds {
                axis,
                dim: self.dim,
            });
        }

        if self.variables[axis].is_fixed() {
            let target_order = self.order.saturating_sub(1).max(1);
            return Self::constant(0.0, self.variables.clone(), target_order);
        }

        let new_order = self.order.saturating_sub(1);
        let new_terms = term_count(self.dim, new_order)?;
        let new_multi_indices = generate_multi_indices(self.dim, new_order);
        let mut new_coeffs = vec![0.0; new_terms];
        // Even a very small bounded value remainder can oscillate arbitrarily
        // quickly. Only a finite singleton establishes a constant offset.
        let mut diff_rem =
            if self.remainder.lo().is_finite() && self.remainder.lo() == self.remainder.hi() {
                Interval::point(0.0)
            } else {
                Interval::WHOLE
            };

        for (coeff, mi) in self.coefficients.iter().zip(&self.multi_indices) {
            if *coeff == 0.0 {
                continue;
            }
            let alpha_k = mi[axis] as usize;
            if alpha_k > 0 {
                let mut beta = mi.clone();
                beta[axis] -= 1;
                if let Some(pos) = new_multi_indices.iter().position(|m| m == &beta) {
                    new_coeffs[pos] = if alpha_k == 1 {
                        *coeff
                    } else {
                        let weight = Interval::point(*coeff) * Interval::point(alpha_k as f64);
                        crate::taylor::split_mid(weight, self.monomial_range(&beta), &mut diff_rem)
                    };
                }
            }
        }

        Ok(Self {
            dim: self.dim,
            order: new_order,
            variables: self.variables.clone(),
            centers: self.centers.clone(),
            coefficients: new_coeffs,
            multi_indices: new_multi_indices,
            remainder: diff_rem,
        })
    }

    /// Compute the gradient vector as a list of partial derivative models `[∂f/∂x₀, ..., ∂f/∂x_{d-1}]`.
    pub fn gradient(&self) -> Result<Vec<Self>, TaylorModelError> {
        (0..self.dim).map(|axis| self.diff(axis)).collect()
    }

    /// Compute certified interval range bounds for each component of the gradient over the domain box.
    pub fn gradient_range(&self) -> Result<Vec<Interval>, TaylorModelError> {
        let grad = self.gradient()?;
        grad.into_iter().map(|g| g.range()).collect()
    }

    /// Certified Lipschitz constant bound over the entire domain box: `L = √(∑ L_k²)`.
    pub fn lipschitz_bound(&self) -> Result<f64, TaylorModelError> {
        let grad_ranges = self.gradient_range()?;
        let scale = grad_ranges
            .iter()
            .map(|iv| iv.lo().abs().max(iv.hi().abs()))
            .fold(0.0_f64, f64::max);
        if scale == 0.0 || !scale.is_finite() {
            return Ok(scale);
        }
        let mut sum_sq = Interval::point(0.0);
        for iv in grad_ranges {
            let max_mag = iv.lo().abs().max(iv.hi().abs());
            let normalized = Interval::point(max_mag) / Interval::point(scale);
            sum_sq = sum_sq + normalized.powi(2);
        }
        Ok((sum_sq.sqrt() * Interval::point(scale)).hi())
    }

    /// Integrate the model along coordinate `axis` over `[a, b] ⊆ D_{axis}`.
    pub fn integrate_axis(
        &self,
        axis: usize,
        a: f64,
        b: f64,
    ) -> Result<Interval, TaylorModelError> {
        if axis >= self.dim {
            return Err(TaylorModelError::AxisIndexOutOfBounds {
                axis,
                dim: self.dim,
            });
        }
        if a > b {
            let rev = self.integrate_axis(axis, b, a)?;
            return Ok(Interval::point(0.0) - rev);
        }

        let var = &self.variables[axis];
        let req_iv = Interval::new(a, b);
        if !var.domain.encloses(req_iv) {
            return Err(TaylorModelError::IncompatibleModels);
        }

        let c_k = self.centers[axis];
        let delta_a = a - c_k;
        let delta_b = b - c_k;

        let mut acc = Interval::point(0.0);

        for (coeff, mi) in self.coefficients.iter().zip(&self.multi_indices) {
            if *coeff == 0.0 {
                continue;
            }
            let alpha_k = mi[axis] as usize;
            let int_k = (delta_b.powi((alpha_k + 1) as i32) - delta_a.powi((alpha_k + 1) as i32))
                / ((alpha_k + 1) as f64);

            let mut term = Interval::point(*coeff * int_k);

            // Bounded monomial across other dimensions
            for (m, &power) in mi.iter().enumerate() {
                if m != axis && power > 0 {
                    let rad = self.variables[m].radius();
                    let centered = Interval::new(-rad, rad);
                    term = term * centered.powi(power as usize);
                }
            }
            acc = acc + term;
        }

        // Remainder integration: R * (b - a)
        let rem_int = self.remainder * Interval::point(b - a);
        Ok(acc + rem_int)
    }

    fn check_compatible(&self, other: &Self) -> Result<(), TaylorModelError> {
        if self.dim != other.dim || self.order != other.order {
            return Err(TaylorModelError::IncompatibleModels);
        }
        for (v1, v2) in self.variables.iter().zip(&other.variables) {
            if v1.name != v2.name || v1.domain != v2.domain {
                return Err(TaylorModelError::IncompatibleModels);
            }
        }
        Ok(())
    }
}
