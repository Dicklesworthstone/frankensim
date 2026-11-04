//! Finite sections of the original mean duct and flat-plate relations.
//!
//! A section uses the difference of cumulative mean transfer, not a fresh
//! thermal entrance at each mesh cell. This is a derived held-property model;
//! both original mean evaluations remain attached and retain their source and
//! domain. Their discrepancy allowance is not a validated local error bound.

use super::{
    CorrelationError, CorrelationId, CorrelationInputs, CorrelationRobinBoundary,
    HeatTransferCoefficient, Length, NusseltEvaluation, Temperature,
    ThermalConductivity, det, evaluate,
};
use fs_evidence::{Evidence, NumericalCertificate, ProvenanceHash, SensitivitySummary};

mod flat_plate;
pub(super) use flat_plate::{flat_plate_laminar_mean, flat_plate_turbulent_mean};

/// A strip must select one of the original flat-plate mean correlations.
#[derive(Debug)]
pub enum FlatPlateSectionError {
    /// A duct, cylinder or natural-convection relation cannot describe this strip.
    UnsupportedCorrelation(CorrelationId),
    /// Preserve the original group, source-domain or arithmetic refusal.
    Correlation(CorrelationError),
}

impl core::fmt::Display for FlatPlateSectionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::UnsupportedCorrelation(id) => write!(f,
                "flat-plate sections require a laminar or turbulent flat-plate mean, not {}", id.name()),
            Self::Correlation(error) => write!(f, "{error}"),
        }
    }
}
impl std::error::Error for FlatPlateSectionError {}
impl From<CorrelationError> for FlatPlateSectionError {
    fn from(error: CorrelationError) -> Self { Self::Correlation(error) }
}

/// Mean heat transfer over the physical interval `[a,b]` of an isothermal plate.
/// Both positions use the SAME leading edge, fluid state and relative velocity.
/// The Nusselt number is based on the downstream distance `b`, not strip width.
#[derive(Debug)]
pub struct FlatPlateSection {
    start_fraction: f64,
    upstream: Option<NusseltEvaluation>,
    downstream: NusseltEvaluation,
    nusselt: f64,
}

impl FlatPlateSection {
    /// Physical upstream/downstream distance ratio `a/b` in `[0,1)`.
    pub const fn start_fraction(&self) -> f64 { self.start_fraction }
    /// Original mean at `a`; absent only at the exact physical leading edge.
    pub const fn upstream(&self) -> Option<&NusseltEvaluation> { self.upstream.as_ref() }
    /// Original mean at the downstream position `b`.
    pub const fn downstream(&self) -> &NusseltEvaluation { &self.downstream }
    /// Strip-average `h b/k`, retaining `b` as its characteristic length.
    pub const fn nusselt(&self) -> f64 { self.nusselt }

    /// Convert with the SAME fluid conductivity and physical downstream distance.
    ///
    /// A leading-edge strip delegates exactly to the original mean conversion.
    /// For positive starts, original source cards remain attached but their
    /// mean discrepancy cannot bound a difference of cumulative predictions:
    /// the derived strip therefore reports an explicit unbounded model band.
    /// Formula arithmetic remains an estimate, with no local sensitivity claim.
    pub fn heat_transfer_coefficient(
        &self,
        conductivity: ThermalConductivity,
        downstream_distance: Length,
    ) -> Result<Evidence<HeatTransferCoefficient>, CorrelationError> {
        let mut evidence = self.downstream
            .heat_transfer_coefficient(conductivity, downstream_distance)?;
        let Some(upstream) = &self.upstream else { return Ok(evidence); };
        let value = self.nusselt * conductivity.value() / downstream_distance.value();
        if !value.is_finite() || value <= 0.0 {
            return Err(CorrelationError::NonFiniteResult {
                stage: "flat-plate section Nu-to-h conversion", value_bits: value.to_bits(),
            });
        }
        evidence.value = HeatTransferCoefficient::new(value);
        evidence.qoi = value;
        evidence.numerical = NumericalCertificate::estimate(value, value);
        evidence.model.discrepancy_rel = f64::INFINITY;
        evidence.model.assumptions.push(
            "section difference of cumulative flat-plate means at one held fluid and isothermal wall state".into());
        evidence.model.assumptions.push(
            "original mean-card discrepancy does not bound local or section model discrepancy".into());
        evidence.model.assumptions.sort_unstable();
        evidence.model.assumptions.dedup();
        evidence.sensitivity = SensitivitySummary::default();
        let mut values = [0_u8; 24];
        for (bytes, input) in values.chunks_exact_mut(8).zip([
            self.start_fraction, conductivity.value(), downstream_distance.value(),
        ]) {
            bytes.copy_from_slice(&input.to_bits().to_le_bytes());
        }
        evidence.provenance = ProvenanceHash::chain("convection-flat-plate-section-nu-to-h-v1", &[
            upstream.evidence().provenance, self.downstream.evidence().provenance,
            ProvenanceHash::of_bytes(&values),
        ]);
        Ok(evidence)
    }

    /// Lower exactly this strip coefficient into the existing conduction owner,
    /// retaining the source-bearing coefficient beside the private Robin row.
    pub fn robin_boundary(
        &self,
        conductivity: ThermalConductivity,
        downstream_distance: Length,
        reference_temperature: Temperature,
    ) -> Result<CorrelationRobinBoundary, CorrelationError> {
        let temperature = reference_temperature.value();
        if !temperature.is_finite() {
            return Err(CorrelationError::InvalidDimensionalInput {
                field: "Robin reference temperature", value_bits: temperature.to_bits(),
                requirement: "supply a finite coherent-SI absolute temperature",
            });
        }
        let coefficient = self.heat_transfer_coefficient(conductivity, downstream_distance)?;
        let boundary = fs_conduction::ThermalBc::robin(coefficient.value.value(), temperature)
            .map_err(|_| CorrelationError::NonFiniteResult {
                stage: "conduction Robin boundary admission", value_bits: coefficient.value.value().to_bits(),
            })?;
        Ok(CorrelationRobinBoundary { coefficient, reference_temperature, boundary })
    }
}

/// Integrate one ORIGINAL flat-plate mean relation between physical positions
/// `a` and `b`, using `start_fraction=a/b` and Reynolds number at `b`.
///
/// For positive starts, both `Re_b` and `Re_a=Re_b*a/b` must satisfy the SAME
/// selected card. The strip transfer is `h b/k = (Nu_b-Nu_a)/(1-a/b)`.
/// A laminar strip rationalizes the square-root difference. A turbulent strip
/// cancels the shared 871 leading-edge correction analytically and uses an
/// accurate exponential difference for very short strips. At `a=0` the exact
/// original whole-plate mean is retained, including its leading-edge physics.
///
/// No positive unsupported upstream point is replaced by the exact entrance.
/// In particular, no laminar/turbulent card switch or smoothing hides their
/// rounded transition mismatch. The turbulent first interval already includes
/// its original laminar leading-edge correction; a later strip must lie fully
/// within that card's turbulent range.
///
/// Heat additivity assumes uniform relative velocity, fluid properties and
/// isothermal wall temperature. Spatially varying wall/film properties are a
/// separate held-property approximation, not a resolved boundary-layer solve.
pub fn evaluate_flat_plate_section(
    correlation: CorrelationId,
    reynolds_at_end: f64,
    prandtl: f64,
    start_fraction: f64,
) -> Result<FlatPlateSection, FlatPlateSectionError> {
    if !matches!(correlation,
        CorrelationId::FlatPlateLaminarAverage | CorrelationId::FlatPlateTurbulentAverage)
    {
        return Err(FlatPlateSectionError::UnsupportedCorrelation(correlation));
    }
    if !start_fraction.is_finite() || !(0.0..1.0).contains(&start_fraction) {
        return Err(CorrelationError::InvalidGroup {
            axis: "section_start_fraction", value_bits: start_fraction.to_bits(),
        }.into());
    }
    let original = |reynolds| evaluate(correlation, CorrelationInputs::forced(reynolds, prandtl));
    let downstream = original(reynolds_at_end)?;
    let upstream = if start_fraction == 0.0 { None }
        else { Some(original(reynolds_at_end * start_fraction)?) };
    let mean = downstream.evidence().value;
    let nusselt = match correlation {
        CorrelationId::FlatPlateLaminarAverage => flat_plate::laminar_section(mean, start_fraction),
        CorrelationId::FlatPlateTurbulentAverage => flat_plate::turbulent_section(mean, prandtl, start_fraction),
        _ => return Err(FlatPlateSectionError::UnsupportedCorrelation(correlation)),
    };
    if !nusselt.is_finite() || nusselt <= 0.0 {
        return Err(CorrelationError::NonFiniteResult {
            stage: "flat-plate section cumulative-transfer difference", value_bits: nusselt.to_bits(),
        }.into());
    }
    Ok(FlatPlateSection { start_fraction, upstream, downstream, nusselt })
}

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
