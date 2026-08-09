//! Composite fan systems (bead frn2i.1): exact piecewise-linear composition
//! of heterogeneous [`FanBank`]s under an explicit series or parallel
//! topology.
//!
//! Each member bank first materializes its EFFECTIVE curve: the declared
//! base curve with the fan-law speed scaling and bank arrangement factors
//! already applied (`q' = q * flow_factor`, `p' = p * pressure_factor`).
//! Composition is then exact, not numerical:
//!
//! - SERIES banks share one total flow and pressures add: the composite
//!   knots are the merged member knots over the shared flow domain, and
//!   each composite pressure is the exact sum of member interpolants.
//! - PARALLEL banks share one pressure and flows add: each member's
//!   piecewise-linear curve has an exact piecewise-linear inverse, the
//!   composite sums member flows over the merged pressure knots of the
//!   shared pressure domain, and the result is re-inverted exactly.
//!
//! The product is a synthetic [`FanCurve`] whose provenance binds every
//! member bank's identity, wrapped as a one-fan series bank at speed ratio
//! 1 so the certified interval-Newton solve path applies unchanged. The
//! wrapper is a composition adapter, not a collapse: the member structure
//! is fully bound into the composite provenance, and the composite speed
//! domain is pinned to exactly 1 so no further scaling can sneak in.

use crate::{
    AirflowError, FanArrangement, FanBank, FanCurve, FanPoint, SourceProvenance, ToleranceBasis,
};
use fs_qty::{Pressure, VolumetricFlowRate};

/// Compose member banks in SERIES: one total flow passes through every
/// bank and pressure rises add. The composite is a synthetic one-fan bank
/// whose curve is the exact merged piecewise-linear pressure sum.
///
/// Members are canonically ordered by identity before composition:
/// series pressure addition is commutative, so declaration order must not
/// move the composite's bytes or identity.
///
/// # Errors
/// Returns [`AirflowError`] for fewer than two banks or when the members
/// share no common flow domain above every member's stall boundary.
pub fn compose_series(banks: &[FanBank]) -> Result<FanBank, AirflowError> {
    if banks.len() < 2 {
        return Err(AirflowError::EmptyFanComposition { topology: "series" });
    }
    let mut banks = banks.to_vec();
    banks.sort_by(compare_banks);
    let banks = banks.as_slice();
    let members: Vec<Vec<(f64, f64)>> = banks.iter().map(effective_curve).collect();
    // Shared flow domain: above every member's (scaled) stall boundary and
    // below every member's last knot.
    let domain_lo = banks
        .iter()
        .map(|bank| bank.curve().admissible_min_flow().value() * bank.flow_factor())
        .fold(0.0_f64, f64::max);
    let domain_hi = members
        .iter()
        .map(|member| member.last().expect("nonempty curve").0)
        .fold(f64::INFINITY, f64::min);
    if !(domain_lo.is_finite() && domain_hi.is_finite() && domain_lo < domain_hi) {
        return Err(AirflowError::NoCommonSeriesDomain {
            low_bits: domain_lo.to_bits(),
            high_bits: domain_hi.to_bits(),
        });
    }
    // Knots: the domain endpoints plus every member knot strictly inside.
    let mut knots: Vec<f64> = vec![domain_lo, domain_hi];
    for member in &members {
        for &(q, _) in member {
            if q > domain_lo && q < domain_hi {
                knots.push(q);
            }
        }
    }
    knots.sort_by(f64::total_cmp);
    knots.dedup_by(|left, right| left.to_bits() == right.to_bits());
    let mut points = Vec::with_capacity(knots.len());
    for &q in &knots {
        let pressure: f64 = members
            .iter()
            .map(|member| interpolate(member, q))
            .sum();
        points.push(FanPoint::new(
            VolumetricFlowRate::new(q),
            Pressure::new(pressure),
        ));
    }
    let curve = synthetic_curve("series", banks, points, domain_lo)?;
    FanBank::new(curve, 1, FanArrangement::Series, 1.0)
}

/// Compose member banks in PARALLEL: every bank sees one shared pressure
/// and flows add. The composite uses the exact piecewise-linear inverse of
/// each member curve over the shared pressure domain.
///
/// # Errors
/// Returns [`AirflowError`] for fewer than two banks or when the members
/// share no common pressure domain that keeps every bank on its declared
/// curve and above its stall boundary.
pub fn compose_parallel(banks: &[FanBank]) -> Result<FanBank, AirflowError> {
    if banks.len() < 2 {
        return Err(AirflowError::EmptyFanComposition {
            topology: "parallel",
        });
    }
    let mut banks = banks.to_vec();
    banks.sort_by(compare_banks);
    let banks = banks.as_slice();
    let members: Vec<Vec<(f64, f64)>> = banks.iter().map(effective_curve).collect();
    // Shared pressure domain: [max member minimum pressure, min member
    // stall pressure]. The stall pressure is the pressure at which the
    // member's (scaled) flow hits its declared admissible minimum.
    let pressure_lo = members
        .iter()
        .map(|member| member.last().expect("nonempty curve").1)
        .fold(f64::NEG_INFINITY, f64::max);
    let pressure_hi = banks
        .iter()
        .zip(&members)
        .map(|(bank, member)| {
            let stall_flow = bank.curve().admissible_min_flow().value() * bank.flow_factor();
            interpolate_inverse(member, stall_flow)
        })
        .fold(f64::INFINITY, f64::min);
    if !(pressure_lo.is_finite() && pressure_hi.is_finite() && pressure_lo < pressure_hi) {
        return Err(AirflowError::NoCommonParallelDomain {
            low_bits: pressure_lo.to_bits(),
            high_bits: pressure_hi.to_bits(),
        });
    }
    // Knots in descending pressure: endpoints plus interior member knots.
    let mut pressures: Vec<f64> = vec![pressure_lo, pressure_hi];
    for member in &members {
        for &(_, p) in member {
            if p > pressure_lo && p < pressure_hi {
                pressures.push(p);
            }
        }
    }
    pressures.sort_by(|left, right| right.total_cmp(left));
    pressures.dedup_by(|left, right| left.to_bits() == right.to_bits());
    // Flows at each pressure: exact member-inverse sums. Descending
    // pressure gives ascending flow, so reverse into flow-ascending order.
    let mut points: Vec<FanPoint> = Vec::with_capacity(pressures.len());
    for &p in &pressures {
        let flow: f64 = members
            .iter()
            .map(|member| interpolate_inverse(member, p))
            .sum();
        points.push(FanPoint::new(
            VolumetricFlowRate::new(flow),
            Pressure::new(p),
        ));
    }
    points.reverse();
    let admissible_min = points.first().expect("nonempty").flow.value();
    let curve = synthetic_curve("parallel", banks, points, admissible_min)?;
    FanBank::new(curve, 1, FanArrangement::Series, 1.0)
}

/// One member bank's effective curve: the declared base curve with the
/// fan-law speed scaling and bank arrangement factors applied.
fn effective_curve(bank: &FanBank) -> Vec<(f64, f64)> {
    let flow_factor = bank.flow_factor();
    let pressure_factor = bank.pressure_factor();
    bank.curve()
        .points()
        .iter()
        .map(|point| {
            (
                point.flow.value() * flow_factor,
                point.pressure.value() * pressure_factor,
            )
        })
        .collect()
}

/// Exact piecewise-linear pressure at flow `q` on a `(q, p)` curve that is
/// flow-ascending; `q` is always inside the declared domain by caller
/// contract.
fn interpolate(curve: &[(f64, f64)], q: f64) -> f64 {
    if q <= curve[0].0 {
        return curve[0].1;
    }
    let last = curve.len() - 1;
    if q >= curve[last].0 {
        return curve[last].1;
    }
    for pair in curve.windows(2) {
        let (q0, p0) = pair[0];
        let (q1, p1) = pair[1];
        if q >= q0 && q <= q1 {
            return p0 + (p1 - p0) * (q - q0) / (q1 - q0);
        }
    }
    curve[last].1
}

/// Exact piecewise-linear inverse: flow at pressure `p` on a `(q, p)`
/// curve whose pressure is non-increasing in flow.
fn interpolate_inverse(curve: &[(f64, f64)], p: f64) -> f64 {
    let last = curve.len() - 1;
    // Pressure at index 0 is the maximum (shutoff side); the last index is
    // the minimum (free-delivery side).
    if p >= curve[0].1 {
        return curve[0].0;
    }
    if p <= curve[last].1 {
        return curve[last].0;
    }
    for pair in curve.windows(2) {
        let (q0, p0) = pair[0];
        let (q1, p1) = pair[1];
        if p <= p0 && p >= p1 {
            return q0 + (q1 - q0) * (p - p0) / (p1 - p0);
        }
    }
    curve[last].0
}

/// The synthetic composite curve: member-bound provenance, the widest
/// member tolerance, the weakest member authority basis, the computed
/// stall boundary, and a speed domain pinned to exactly 1 (member speeds
/// are already folded into the effective curves).
fn synthetic_curve(
    topology: &'static str,
    banks: &[FanBank],
    points: Vec<FanPoint>,
    admissible_min_flow: f64,
) -> Result<FanCurve, AirflowError> {
    let names: Vec<&str> = banks
        .iter()
        .map(|bank| bank.curve().name())
        .collect();
    let name = format!("{topology}({})", names.join("+"));
    let mut hasher = fs_blake3::DomainHasher::new(COMPOSITE_SOURCE_DOMAIN);
    hasher.update(topology.as_bytes());
    for bank in banks {
        let source = bank.curve().source();
        hasher.update(&(source.citation.len() as u64).to_le_bytes());
        hasher.update(source.citation.as_bytes());
        hasher.update(&(source.identifier.len() as u64).to_le_bytes());
        hasher.update(source.identifier.as_bytes());
        hasher.update(&(bank.curve().points().len() as u64).to_le_bytes());
        for point in bank.curve().points() {
            hasher.update(&point.flow.value().to_le_bytes());
            hasher.update(&point.pressure.value().to_le_bytes());
        }
        hasher.update(&(bank.count() as u64).to_le_bytes());
        hasher.update(&[match bank.arrangement() {
            FanArrangement::Series => 0,
            FanArrangement::Parallel => 1,
        }]);
        hasher.update(&bank.speed_ratio().to_le_bytes());
    }
    let identity = hasher.finalize();
    let source = SourceProvenance::new(
        format!(
            "composite {topology} fan system of {} member banks",
            banks.len()
        ),
        format!("fan-system-composite:v1:{identity}"),
    );
    let tolerance = banks
        .iter()
        .map(|bank| bank.curve().pressure_tolerance_rel())
        .fold(0.0_f64, f64::max);
    let basis = banks
        .iter()
        .map(|bank| bank.curve().tolerance_basis())
        .min_by_key(|basis| basis_strength(*basis))
        .unwrap_or(ToleranceBasis::Analytic);
    FanCurve::new(
        name,
        points,
        source,
        tolerance,
        basis,
        VolumetricFlowRate::new(admissible_min_flow),
        (1.0, 1.0),
    )
}

/// Authority-strength ordering: analytic is strongest, manufacturer next,
/// engineering allowance weakest; the composite takes the weakest member.
const fn basis_strength(basis: ToleranceBasis) -> u8 {
    match basis {
        ToleranceBasis::Analytic => 2,
        ToleranceBasis::ManufacturerDeclared => 1,
        ToleranceBasis::EngineeringAllowance => 0,
    }
}

const COMPOSITE_SOURCE_DOMAIN: &str = "org.frankensim.fs-airflow.fan-system-composite.v1";

/// Canonical member ordering for composition: both topologies are
/// physically commutative, so members sort by their full declared identity
/// and declaration order cannot change the composite.
fn compare_banks(left: &FanBank, right: &FanBank) -> core::cmp::Ordering {
    left.curve()
        .name()
        .cmp(right.curve().name())
        .then(left.count().cmp(&right.count()))
        .then(match (left.arrangement(), right.arrangement()) {
            (FanArrangement::Series, FanArrangement::Series)
            | (FanArrangement::Parallel, FanArrangement::Parallel) => core::cmp::Ordering::Equal,
            (FanArrangement::Series, FanArrangement::Parallel) => core::cmp::Ordering::Less,
            (FanArrangement::Parallel, FanArrangement::Series) => core::cmp::Ordering::Greater,
        })
        .then(left.speed_ratio().total_cmp(&right.speed_ratio()))
}
