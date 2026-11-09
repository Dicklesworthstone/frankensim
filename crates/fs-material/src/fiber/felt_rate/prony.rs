//! Cold, source-checked Prony-to-retardation conversion.
//!
//! Equal source times are combined algebraically; no time or weight is fitted.
//! Dense eigenvectors must resolve every active constitutive equation, including
//! weak rows that an absolute eigensolver residual can otherwise hide.

use super::{CreepBranch, MAX_FELT_RATE_BRANCHES, RateFelt, WoolFelt, bad, dd, finite};
use crate::{MaterialError, visco::GeneralizedMaxwell};
use fs_math::dd::Dd;

// Products of positive finite source factors, with the binary scale carried
// separately. A modulus, viscosity or compliance may be representable even
// when an intermediate square, product or quotient is not.
fn product_ratio(numerator: &[f64], denominator: &[f64]) -> f64 {
    fn fraction(value: f64) -> (f64, i32) {
        let shift = if value < f64::MIN_POSITIVE { 52 } else { 0 };
        let bits = (value
            * if shift == 0 {
                1.0
            } else {
                4_503_599_627_370_496.0
            })
        .to_bits();
        let exponent = ((bits >> 52) & 0x7ff) as i32 - 1023 - shift;
        let mantissa = f64::from_bits((bits & 0x000f_ffff_ffff_ffff) | (1023_u64 << 52));
        (mantissa, exponent)
    }
    let mut value = 1.0;
    let mut exponent = 0;
    for &factor in numerator {
        if factor == 0.0 {
            return 0.0;
        }
        let (mantissa, power) = fraction(factor);
        value *= mantissa;
        exponent += power;
    }
    for &factor in denominator {
        let (mantissa, power) = fraction(factor);
        value /= mantissa;
        exponent -= power;
    }
    while exponent != 0 {
        let shift = exponent.clamp(-512, 512);
        value *= f64::from_bits(((1023 + shift) as u64) << 52);
        exponent -= shift;
    }
    value
}

fn admitted_branch(modulus: f64, viscosity: f64) -> Result<CreepBranch, MaterialError> {
    let tau = viscosity / modulus;
    if !modulus.is_finite()
        || modulus <= 0.0
        || !viscosity.is_finite()
        || viscosity <= 0.0
        || !tau.is_finite()
        || tau <= 0.0
    {
        return Err(bad("unrepresentable original Prony compliance arm"));
    }
    Ok(CreepBranch {
        modulus_pa: modulus,
        viscosity_pa_s: viscosity,
    })
}

impl RateFelt {
    /// Algebraic Prony-to-retardation realization of the existing source card.
    /// Its instantaneous linear spring is replaced by the supplied WoolFelt
    /// envelope. This composition does not identify a nonlinear material fit.
    ///
    /// Exactly equal relaxation times combine into one observable arm; their
    /// omitted stress differences have identically zero compliance coupling.
    /// Positive weak weights are retained when numerically resolved. A dense
    /// eigenvector that fails a relative original-source equation check, or an
    /// unresolved coupling, refuses instead of silently changing the spectrum.
    /// The eight-active-arm cap applies before equal-time coalescence.
    pub fn from_prony(
        envelope: WoolFelt,
        card: &GeneralizedMaxwell,
    ) -> Result<Self, MaterialError> {
        Self::new(envelope, realize(card)?)
    }
}

/// Retardation branches of the original source compliance, independent of any
/// chosen nonlinear spring or structural geometry.
pub(super) fn realize(card: &GeneralizedMaxwell) -> Result<Vec<CreepBranch>, MaterialError> {
    GeneralizedMaxwell::new(card.e_inf, card.terms.clone()).map_err(|error| {
        MaterialError::Parameters {
            what: error.to_string(),
        }
    })?;
    let mut active: Vec<_> = card
        .terms
        .iter()
        .copied()
        .filter(|term| term.0 > 0.0)
        .collect();
    if active.len() > MAX_FELT_RATE_BRANCHES {
        return Err(bad("Prony material arm cap"));
    }
    if active.is_empty() {
        return Ok(Vec::new());
    }
    active.sort_by(|a, b| a.1.total_cmp(&b.1));
    let mut terms = Vec::with_capacity(active.len());
    let mut at = 0;
    while at < active.len() {
        let tau = active[at].1;
        let mut modulus = Dd::ZERO;
        while at < active.len() && active[at].1 == tau {
            modulus = modulus + dd(active[at].0);
            at += 1;
        }
        terms.push((
            finite(modulus.to_f64(), "equal-time Prony modulus overflow")?,
            tau,
        ));
    }
    let excess = finite(
        terms
            .iter()
            .fold(Dd::ZERO, |sum, term| sum + dd(term.0))
            .to_f64(),
        "Prony excess modulus overflow",
    )?;
    let instant = finite(
        terms
            .iter()
            .fold(dd(card.e_inf), |sum, term| sum + dd(term.0))
            .to_f64(),
        "instantaneous modulus overflow",
    )?;
    if terms.len() == 1 {
        // E_K = E_inf*E_0/E_1, eta_K = tau_1*E_0^2/E_1.
        // The positive equilibrium modulus must not be reconstructed by
        // subtracting two nearly equal instantaneous stiffnesses.
        let (modulus, tau) = terms[0];
        let branch = admitted_branch(
            product_ratio(&[card.e_inf, instant], &[modulus]),
            product_ratio(&[tau, instant, instant], &[modulus]),
        )?;
        return Ok(vec![branch]);
    }

    let n = terms.len();
    let mut h = vec![0.0; n * n];
    let mut r = vec![0.0; n * n];
    for i in 0..n {
        r[i * n + i] = finite(
            product_ratio(&[terms[i].0, terms[i].1], &[]),
            "Prony viscosity overflow",
        )?;
        if r[i * n + i] <= 0.0 {
            return Err(bad("unresolved positive Prony viscosity"));
        }
        let complement = terms
            .iter()
            .enumerate()
            .filter(|(j, _)| *j != i)
            .fold(dd(card.e_inf), |sum, (_, term)| sum + dd(term.0))
            .to_f64();
        for j in 0..n {
            // H = diag(E_i) - e e^T/E_0, but form each diagonal as
            // E_i*(E_inf + sum_{j!=i} E_j)/E_0, with no cancellation.
            h[i * n + j] = finite(
                if i == j {
                    product_ratio(&[terms[i].0, complement], &[instant])
                } else {
                    -product_ratio(&[terms[i].0, terms[j].0], &[instant])
                },
                "Prony compliance matrix overflow",
            )?;
        }
    }
    let modes = fs_modal::eigh_gen_dense(&h, &r, n).map_err(|error| MaterialError::Parameters {
        what: error.to_string(),
    })?;
    let mut creep = Vec::with_capacity(n);
    let mut sum = Dd::ZERO;
    for mode in modes {
        if !mode.lambda.is_finite()
            || mode.lambda <= 0.0
            || !mode.residual.is_finite()
            || mode.residual > 1e-6 * mode.lambda
            || mode.phi.len() != n
            || mode.phi.iter().any(|x| !x.is_finite())
        {
            return Err(bad("unresolved positive Prony retardation spectrum"));
        }
        let mut coupling = Dd::ZERO;
        let mut magnitude = Dd::ZERO;
        for (&v, term) in mode.phi.iter().zip(&terms) {
            let component = finite(
                product_ratio(&[v.abs(), term.0], &[instant]),
                "unresolved original Prony source coupling",
            )?;
            magnitude = magnitude + dd(component);
            coupling = coupling + dd(component.copysign(v));
        }
        let coupling = finite(
            coupling.to_f64(),
            "unresolved original Prony source coupling",
        )?;
        let magnitude = finite(
            magnitude.to_f64(),
            "unresolved original Prony source coupling",
        )?;
        // Distinct positive-time arms have no exactly unobservable mode.
        // Do not mistake rounded zero or severe cancellation for a dark
        // state; equal-time dark states were removed algebraically above.
        if coupling == 0.0 || 32.0 * f64::EPSILON * magnitude > 1e-9 * coupling.abs() {
            return Err(bad("unresolved relative original Prony source coupling"));
        }

        // Divide H*phi=lambda*R*phi by each POSITIVE original E_i:
        //   (1-lambda*tau_i)*phi_i = sum_j(E_j/E_0)*phi_j.
        // This unweighted row check catches lost weak-row rotations that
        // the mass-norm residual and equilibrium-compliance sum cannot.
        // Near a source pole, a rounded eigenvalue itself contributes
        // delta_lambda*tau_i*phi_i; permit that independently bounded
        // amount, without relaxing a row whose eigenvector entry is zero.
        let lambda_error = mode.residual + 8.0 * f64::EPSILON * mode.lambda;
        for (&v, &(_, tau)) in mode.phi.iter().zip(&terms) {
            let component = (dd(1.0) - dd(mode.lambda) * dd(tau)) * dd(v);
            let value = finite(
                component.to_f64(),
                "unresolved original Prony source equation",
            )?;
            let residual = finite(
                (component - dd(coupling)).to_f64().abs(),
                "unresolved original Prony source equation",
            )?;
            let rounding = finite(
                product_ratio(&[lambda_error, tau, v.abs()], &[]),
                "unresolved original Prony source equation",
            )?;
            let scale = value.abs().max(coupling.abs());
            if residual > rounding + 1e-9 * scale {
                return Err(bad("unresolved relative original Prony source equation"));
            }
        }

        let absolute_coupling = coupling.abs();
        let weight = finite(
            product_ratio(&[absolute_coupling, absolute_coupling], &[mode.lambda]),
            "original Prony compliance overflow",
        )?;
        let branch = admitted_branch(
            product_ratio(&[mode.lambda], &[absolute_coupling, absolute_coupling]),
            product_ratio(&[1.0], &[absolute_coupling, absolute_coupling]),
        )?;
        if weight <= 0.0 {
            return Err(bad("unrepresentable original Prony compliance arm"));
        }
        sum = sum + dd(weight);
        creep.push(branch);
    }
    let target = finite(
        product_ratio(&[excess], &[card.e_inf, instant]),
        "equilibrium creep compliance overflow",
    )?;
    let total = finite(sum.to_f64(), "total Prony creep compliance overflow")?;
    let scale = target.max(total);
    if target <= 0.0 || (sum - dd(target)).to_f64().abs() > 1e-9 * scale {
        return Err(bad(
            "original equilibrium creep compliance was not retained",
        ));
    }
    Ok(creep)
}
