//! Finite axial sections of the original mean Hausen relation.
//!
//! A section uses the difference of cumulative mean transfer, not a fresh
//! thermal entrance at each mesh cell. This is a derived held-property model;
//! both original mean evaluations remain attached and retain their source and
//! domain. Their discrepancy allowance is not a validated local error bound.

use super::{
    CorrelationError, CorrelationId, CorrelationInputs, HeatTransferCoefficient,
    Length, NusseltEvaluation, ThermalConductivity, det, evaluate,
};

const DEVELOPED: f64 = 3.66;
const ENTRANCE: f64 = 0.0668;
const DENOMINATOR: f64 = 0.04;

// The ordinary whole-length card and its derived section share one formula.
pub(super) fn mean_nusselt(graetz: f64) -> f64 {
    DEVELOPED + ENTRANCE * graetz
        / (1.0 + DENOMINATOR * det::pow(graetz, 2.0 / 3.0))
}

/// Mean coefficient over one positive-length section of a uniform circular
/// duct with a hydrodynamically developed inlet and constant wall temperature.
/// All lengths are measured from the SAME heated entrance in diameters.
#[derive(Debug)]
pub struct CircularDuctSection {
    start_over_diameter: f64,
    end_over_diameter: f64,
    upstream: Option<NusseltEvaluation>,
    downstream: NusseltEvaluation,
    nusselt: f64,
}

impl CircularDuctSection {
    pub const fn start_over_diameter(&self) -> f64 { self.start_over_diameter }
    pub const fn end_over_diameter(&self) -> f64 { self.end_over_diameter }
    /// Original cumulative-mean evaluation; absent only at the exact entrance.
    pub const fn upstream(&self) -> Option<&NusseltEvaluation> { self.upstream.as_ref() }
    pub const fn downstream(&self) -> &NusseltEvaluation { &self.downstream }
    /// Derived section-average Nu, based on the original duct diameter.
    pub const fn nusselt(&self) -> f64 { self.nusselt }

    /// Convert the section mean with the SAME fluid conductivity and diameter.
    /// The result has the attached mean-card authority only through this receipt;
    /// neither its numerical error nor its model discrepancy is enclosed.
    pub fn heat_transfer_coefficient(
        &self,
        conductivity: ThermalConductivity,
        diameter: Length,
    ) -> Result<HeatTransferCoefficient, CorrelationError> {
        // Reuse original dimensional admission without manufacturing a new card.
        self.downstream.heat_transfer_coefficient(conductivity, diameter)?;
        let value = self.nusselt * conductivity.value() / diameter.value();
        if !value.is_finite() || value <= 0.0 {
            return Err(CorrelationError::NonFiniteResult {
                stage: "developing section Nu-to-h conversion",
                value_bits: value.to_bits(),
            });
        }
        Ok(HeatTransferCoefficient::new(value))
    }
}

/// Integrate the original Hausen mean from `start_over_diameter` to
/// `end_over_diameter`, retaining its original Re, Pr and heated-length support.
///
/// With `s=x/D` and `P=Re Pr`, cumulative dimensionless transfer is
/// `s Nu_mean(s) = 3.66 s + 0.0668 P s^(2/3)/(s^(2/3)+0.04 P^(2/3))`.
/// Its exact limit is zero at s=0, so the first interval reproduces the original
/// downstream mean without querying an inadmissible zero-length card. A positive
/// upstream distance MUST satisfy the original card; no near-zero floor exists.
///
/// For positive starts, the difference quotient uses
/// `(b^(2/3)-a^(2/3))/(b-a) = (cbrt(a)+cbrt(b)) /
/// (cbrt(a)^2+cbrt(a)cbrt(b)+cbrt(b)^2)`, avoiding cancellation for short cells.
///
/// Subdivision additivity holds for the same Re, Pr, diameter, conductivity and
/// wall state. Re-querying properties or bulk temperatures per section is an
/// explicitly local held-property approximation, not an exact variable-property
/// boundary-layer solution. This does not supply a local model error bound,
/// a new correlation, an adjoint, or an extension of original source support.
pub fn evaluate_circular_duct_section(
    reynolds: f64,
    prandtl: f64,
    start_over_diameter: f64,
    end_over_diameter: f64,
) -> Result<CircularDuctSection, CorrelationError> {
    if !start_over_diameter.is_finite() || start_over_diameter < 0.0 {
        return Err(CorrelationError::InvalidGroup {
            axis: "section_start_over_Dh",
            value_bits: start_over_diameter.to_bits(),
        });
    }
    if !end_over_diameter.is_finite() || end_over_diameter <= start_over_diameter {
        return Err(CorrelationError::InvalidGroup {
            axis: "section_end_over_Dh",
            value_bits: end_over_diameter.to_bits(),
        });
    }
    let original = |length| evaluate(
        CorrelationId::CircularDuctHausen,
        CorrelationInputs::forced(reynolds, prandtl).with_length_ratio(length),
    );
    let downstream = original(end_over_diameter)?;
    let upstream = if start_over_diameter == 0.0 {
        None
    } else {
        Some(original(start_over_diameter)?)
    };
    let nusselt = if upstream.is_none() {
        downstream.evidence().value
    } else {
        let a = det::pow(start_over_diameter, 1.0 / 3.0);
        let b = det::pow(end_over_diameter, 1.0 / 3.0);
        let peclet = reynolds * prandtl;
        let w = DENOMINATOR * det::pow(peclet, 2.0 / 3.0);
        let slope = (a + b) / (a * a + a * b + b * b);
        DEVELOPED + ENTRANCE * peclet * (w / (a * a + w))
            * (slope / (b * b + w))
    };
    if !nusselt.is_finite() || nusselt <= 0.0 {
        return Err(CorrelationError::NonFiniteResult {
            stage: "developing section cumulative-transfer difference",
            value_bits: nusselt.to_bits(),
        });
    }
    Ok(CircularDuctSection {
        start_over_diameter, end_over_diameter, upstream, downstream, nusselt,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn near(actual: f64, expected: f64) {
        assert!((actual - expected).abs() <= 2e-12 * expected.abs().max(1.0),
            "{actual} != {expected}");
    }

    #[test]
    fn g1_g3_section_sums_recover_original_cumulative_mean_without_entrance_resets() {
        let mean = evaluate(CorrelationId::CircularDuctHausen,
            CorrelationInputs::forced(100.0, 8.0).with_length_ratio(30.0)).unwrap();
        let whole = evaluate_circular_duct_section(100.0, 8.0, 0.0, 30.0).unwrap();
        assert_eq!(whole.nusselt(), mean.evidence().value);
        assert!(whole.upstream().is_none());
        for divisions in [2, 8, 32, 128] {
            let length = 30.0 / divisions as f64;
            let sum: f64 = (0..divisions).map(|i| {
                evaluate_circular_duct_section(
                    100.0, 8.0, i as f64 * length, (i + 1) as f64 * length,
                ).unwrap().nusselt() * length
            }).sum();
            near(sum, 30.0 * mean.evidence().value);
            let reset = evaluate_circular_duct_section(100.0, 8.0, 0.0, length).unwrap();
            assert!(30.0 * reset.nusselt() > sum);
        }
        let section = evaluate_circular_duct_section(100.0, 8.0, 10.0, 30.0).unwrap();
        let upper = section.downstream().evidence().value;
        let lower = section.upstream().unwrap().evidence().value;
        near(section.nusselt(), (30.0 * upper - 10.0 * lower) / 20.0);
        assert_eq!(section.downstream().card().id, CorrelationId::CircularDuctHausen);
        near(section.heat_transfer_coefficient(
            ThermalConductivity::new(0.5), Length::new(0.03),
        ).unwrap().value(), section.nusselt() * 0.5 / 0.03);
    }

    #[test]
    fn g1_short_section_retains_the_analytic_local_limit() {
        let start = 10.0_f64;
        let tiny = evaluate_circular_duct_section(100.0, 8.0, start, start.next_up()).unwrap();
        // Differentiate x*Nu_mean(Re*Pr/x) independently. The 3.66 limit stays
        // intact; only the entrance enhancement contributes the derivative.
        let gz = 80.0_f64;
        let b = 0.04 * gz.powf(2.0 / 3.0);
        let expected = 3.66 + (0.0668 * gz / (1.0 + b)) * (2.0 * b / 3.0) / (1.0 + b);
        near(tiny.nusselt(), expected);
        assert!(tiny.nusselt() > 3.66);
        assert!(tiny.nusselt() < tiny.downstream().evidence().value);
    }

    #[test]
    fn g4_original_domains_remain_active_at_both_section_endpoints() {
        assert!(evaluate_circular_duct_section(100.0, 8.0, 0.0, 0.05).is_ok());
        assert!(evaluate_circular_duct_section(100.0, 8.0, 0.049, 1.0).is_err());
        assert!(evaluate_circular_duct_section(100.0, 8.0, 0.0, 0.049).is_err());
        assert!(evaluate_circular_duct_section(100.0, 8.0, 999.0, 1000.0).is_ok());
        for (start, end) in [(-1.0, 1.0), (1.0, 1.0), (2.0, 1.0),
            (0.0, 1001.0), (f64::NAN, 1.0), (0.0, f64::INFINITY)]
        {
            assert!(evaluate_circular_duct_section(100.0, 8.0, start, end).is_err());
        }
        for (re, pr) in [(0.0, 8.0), (2301.0, 8.0), (100.0, 0.5), (100.0, 1001.0)] {
            assert!(evaluate_circular_duct_section(re, pr, 0.0, 1.0).is_err());
        }
        let section = evaluate_circular_duct_section(100.0, 8.0, 0.0, 1.0).unwrap();
        assert!(section.heat_transfer_coefficient(
            ThermalConductivity::new(0.5), Length::new(0.0),
        ).is_err());
    }
}
