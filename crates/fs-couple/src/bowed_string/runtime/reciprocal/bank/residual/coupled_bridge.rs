//! Two support laws through the SAME omitted transverse modes. No inverse of
//! G is needed, including parallel axes at one station (rank-one compliance).
//! Eliminate the first spring, then the second. If h = k0_eff G01, the elastic
//! factors are e0 and e1-h e0. Their positive coefficients give K_eff exactly;
//! the dashpots give T^T C T, T=(I+G K)^-1, to first order in frequency.

use super::super::polarization::StringBodyTranslationConnection;
use super::*;
use fs_math::dd::Dd;

mod columns;
pub(crate) use columns::condensed_columns;

#[derive(Clone, Copy, Debug)]
pub(crate) struct Pair {
    pub first: CondensedBridge,
    pub second: CondensedBridge,
    /// e1 - mix*e0 is the second elastic/dissipative coordinate.
    pub mix: f64,
    /// e0 - back*(e1-mix*e0) is the first dissipative coordinate.
    pub back: f64,
    pub coupled: bool,
}

impl Pair {
    fn new(g: [f64; 3], k: [f64; 2], c: [f64; 2]) -> Result<Self, ReciprocalStringError> {
        if g.iter().any(|x| !x.is_finite())
            || g[0] < 0.0
            || g[2] < 0.0
            || k.iter().chain(&c).any(|x| !x.is_finite() || *x < 0.0)
            || (k[0] == 0.0 && c[0] == 0.0)
            || (k[1] == 0.0 && c[1] == 0.0)
        {
            return Err(invalid(
                "coupled residual bridges require finite passive support laws and compliance",
            ));
        }
        let mut first = condense(g[0], k[0], c[0])?;
        if g[1] == 0.0 {
            return Ok(Self {
                first,
                second: condense(g[2], k[1], c[1])?,
                mix: 0.0,
                back: 0.0,
                coupled: false,
            });
        }
        if g[0] == 0.0 || g[2] == 0.0 {
            return Err(invalid(
                "coupled residual compliance is not positive semidefinite",
            ));
        }
        let rho2 =
            (Dd::from_f64(g[1]) / Dd::from_f64(g[0])) * (Dd::from_f64(g[1]) / Dd::from_f64(g[2]));
        let independent = Dd::ONE - rho2;
        if !rho2.hi.is_finite()
            || !rho2.lo.is_finite()
            || independent.to_f64() < -128.0 * f64::EPSILON
        {
            return Err(invalid(
                "coupled residual compliance is not positive semidefinite",
            ));
        }
        // g11 - k0_eff*g01^2, without catastrophic subtraction for collinear
        // supports at large stiffness. Only roundoff in a unit correlation is
        // clipped; no compliance/stiffness floor or inverse of a singular G.
        // The shared double-double arithmetic also retains a genuinely weak
        // nonparallel direction; subtracting rounded square roots would not.
        let independent = if independent.to_f64() < 0.0 {
            Dd::ZERO
        } else {
            independent
        };
        let first_factor = 1.0 + k[0] * g[0];
        let remaining = if k[0] == 0.0 {
            g[2]
        } else {
            (Dd::from_f64(g[2])
                * (independent + (Dd::ONE - independent) / Dd::from_f64(first_factor)))
            .to_f64()
        };
        let mut second = condense(remaining, k[1], c[1])?;
        let second_factor = 1.0 + k[1] * remaining;
        // Divide successively so a finite weak loss is not erased by squaring
        // a large denominator. The uncoupled scalar lane keeps its arithmetic.
        first.damping_n_s_m = (c[0] / first_factor) / first_factor;
        second.damping_n_s_m = (c[1] / second_factor) / second_factor;
        let mix = checked(first.stiffness_n_m * g[1])?;
        let back = checked(second.stiffness_n_m * g[1])?;
        checked(1.0 + mix * back)?;
        Ok(Self {
            first,
            second,
            mix,
            back,
            coupled: true,
        })
    }
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct BridgePair {
    pub directed: usize,
    pub law: Pair,
}

/// The existing one-bridge lane plus one independently directed support.
/// Source order remains fixed. Orthogonality uses the original admission
/// tolerance, preserving that already supported lane's arithmetic bit for bit.
pub(crate) fn pairs(
    strings: &[MaterialStringMember],
    cards: &[BowedStringCard],
    directed: &[StringBodyTranslationConnection],
    gate: &CancelGate,
) -> Result<Vec<Option<BridgePair>>, ReciprocalStringError> {
    poll(gate)?;
    if strings.len() != cards.len() || strings.iter().any(|s| s.connections.len() != 1) {
        return Err(invalid(
            "residual flexibility requires one original bridge per material string",
        ));
    }
    let mut pairs = vec![None; cards.len()];
    for (index, link) in directed.iter().enumerate() {
        poll(gate)?;
        let member = strings
            .get(link.string)
            .ok_or_else(|| invalid("residual directed bridge names an unknown string"))?;
        if pairs[link.string].is_some() {
            return Err(invalid(
                "residual flexibility supports at most two bridge connections per material string",
            ));
        }
        let card = &cards[link.string];
        let original = &member.connections[0];
        let axis = directions(&member.source).1;
        let other = link.direction.components();
        let dot = dot2(axis, other);
        let cross = if dot.abs() <= 1e-12 {
            0.0
        } else {
            dot * residual_compliance_m_n(card, original.station_fraction, link.station_fraction)?
        };
        let g = [
            dot2(axis, axis)
                * residual_compliance_m_n(
                    card,
                    original.station_fraction,
                    original.station_fraction,
                )?,
            cross,
            dot2(other, other)
                * residual_compliance_m_n(card, link.station_fraction, link.station_fraction)?,
        ];
        let law = Pair::new(
            g,
            [original.stiffness_n_m, link.stiffness_n_m],
            [original.damping_n_s_m, link.damping_n_s_m],
        )?;
        pairs[link.string] = Some(BridgePair {
            directed: index,
            law,
        });
    }
    Ok(pairs)
}

fn checked(value: f64) -> Result<f64, ReciprocalStringError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid(
            "coupled residual bridge coefficient is not representable",
        ))
    }
}

#[cfg(test)]
mod tests;
