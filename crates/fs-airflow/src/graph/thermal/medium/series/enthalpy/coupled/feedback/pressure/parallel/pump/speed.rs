//! Calibrated variable-speed pumps on the original hydraulic/thermal owners.
//!
//! Each speed slice supplies its own head, shaft input and electrical input.
//! Between adjacent speeds, each quantity is linearly interpolated at the SAME
//! inlet volume flow, over the intersection of both slices' source supports.
//! No affinity law, inferred efficiency, extrapolation or loss rescaling is
//! applied. The calibration source must justify this interpolation model.
//!
//! `at_speed` prepares an ordinary `PumpHeadCurve` once per command. It preserves
//! all relevant original flow breakpoints, so repeated operating-point and
//! implicit thermal trials retain the original bounded curve-evaluation cost.
//! The resulting curve works with existing fixed-wall, common-solid, finite
//! reservoir and stored-solid entrypoints. Every trial still checks electrical
//! >= shaft >= hydraulic input/output and the existing heat destinations.
//!
//! A command is held over an admitted quasi-steady interval. This map supplies
//! no rotor inertia, acceleration work, cavitation or motor-control dynamics;
//! changing command does not invent a transition-energy model.
//! An explicit stopped slice can report standby consumption at zero flow;
//! network marching still requires its existing bracketed flowing solution.

use super::power::{PumpDriveCurve, PumpPowerPoint};
use super::{PumpHeadCurve, PumpHeadPoint, PumpOperatingError, SourceProvenance};
use fs_qty::{Length, Power, VolumetricFlowRate};

/// One independently calibrated head and drive-power slice.
#[derive(Clone, Debug, PartialEq)]
pub struct PumpSpeedCurve {
    /// Measured or authored operating speed in revolutions per minute.
    pub speed_rpm: f64,
    /// Original head data with mandatory shaft/electrical input data attached.
    pub curve: PumpHeadCurve,
}

/// Bounded speed family with explicit interpolation applicability.
#[derive(Clone, Debug, PartialEq)]
pub struct PumpSpeedMap {
    curves: Vec<PumpSpeedCurve>,
    source: SourceProvenance,
}

impl PumpSpeedMap {
    /// Admit 2..=64 strictly increasing, finite nonnegative speeds. Every slice
    /// must have its own drive data, the same gravity convention and the same
    /// declared pump/drive heat routing. Adjacent slices must share a nonzero
    /// flow interval and inclusive density/viscosity support. A stopped-speed
    /// slice is used only if explicitly supplied; it is never synthesized.
    /// At zero speed, head and shaft work must be zero; explicitly calibrated
    /// electrical standby consumption is retained.
    pub fn try_new(
        curves: Vec<PumpSpeedCurve>,
        source: SourceProvenance,
    ) -> Result<Self, PumpOperatingError> {
        if !(2..=64).contains(&curves.len())
            || source.citation.trim().is_empty()
            || source.identifier.trim().is_empty()
            || curves.iter().any(|slice| {
                !slice.speed_rpm.is_finite()
                    || slice.speed_rpm < 0.0
                    || slice
                        .curve
                        .speed_rpm()
                        .is_some_and(|speed| speed != slice.speed_rpm)
            })
            || curves.windows(2).any(|p| p[0].speed_rpm >= p[1].speed_rpm)
        {
            return Err(PumpOperatingError::Invalid(
                "speed map needs bounded ascending calibrated speeds and an interpolation source",
            ));
        }
        let first = &curves[0].curve;
        let routing = drive(first)?.routing();
        for slice in &curves {
            let input = drive(&slice.curve)?;
            if slice.curve.gravity_m_per_s2() != first.gravity_m_per_s2()
                || input.routing() != routing
            {
                return Err(PumpOperatingError::Invalid(
                    "speed slices must share gravity convention and explicit loss routing",
                ));
            }
            if slice.speed_rpm == 0.0
                && (slice.curve.points().iter().any(|p| p.head.value() != 0.0)
                    || input.points().iter().any(|p| p.shaft_input.value() != 0.0))
            {
                return Err(PumpOperatingError::Invalid(
                    "a stopped pump cannot supply head or consume mechanical shaft power",
                ));
            }
        }
        for pair in curves.windows(2) {
            common_support(&pair[0].curve, &pair[1].curve)?;
        }
        Ok(Self { curves, source })
    }

    /// Original calibrated slices, in ascending speed order.
    pub fn curves(&self) -> &[PumpSpeedCurve] {
        &self.curves
    }

    /// Source supporting the family and same-flow speed interpolation.
    pub const fn source(&self) -> &SourceProvenance {
        &self.source
    }

    /// Inclusive speed support in revolutions per minute.
    pub fn speed_domain_rpm(&self) -> (f64, f64) {
        (
            self.curves[0].speed_rpm,
            self.curves[self.curves.len() - 1].speed_rpm,
        )
    }

    /// Prepare a held-speed curve for all existing pump/network/thermal APIs.
    /// Exact calibration speeds retain their original samples and full support.
    /// Interior speeds retain only common source support and every original
    /// head/power breakpoint within it, with at most 1024 knots per curve.
    /// Original per-trial power-feasibility and material gates remain active.
    pub fn at_speed(&self, speed_rpm: f64) -> Result<PumpHeadCurve, PumpOperatingError> {
        let (lower, upper) = self.speed_domain_rpm();
        if !speed_rpm.is_finite() || speed_rpm < lower || speed_rpm > upper {
            return Err(PumpOperatingError::OutsideCurve {
                field: "calibrated pump speed (rpm)",
                value: speed_rpm,
                lower,
                upper,
            });
        }
        let right = self.curves.partition_point(|p| p.speed_rpm < speed_rpm);
        let high = &self.curves[right];
        if high.speed_rpm == speed_rpm {
            let mut selected = high.curve.clone();
            selected.speed_rpm = Some(speed_rpm);
            return Ok(selected);
        }
        let low = &self.curves[right - 1];
        let fraction = (speed_rpm - low.speed_rpm) / (high.speed_rpm - low.speed_rpm);
        let blend = |a: f64, b: f64| (1.0 - fraction) * a + fraction * b;
        let support = common_support(&low.curve, &high.curve)?;
        let low_drive = drive(&low.curve)?;
        let high_drive = drive(&high.curve)?;

        let heads = knots(
            low.curve
                .points()
                .iter()
                .chain(high.curve.points())
                .map(|p| p.flow.value()),
            support.flow,
        )?
        .into_iter()
        .map(|q| {
            Ok(PumpHeadPoint {
                flow: VolumetricFlowRate::new(q),
                head: Length::new(blend(
                    low.curve.head_at_flow(q)?,
                    high.curve.head_at_flow(q)?,
                )),
            })
        })
        .collect::<Result<Vec<_>, PumpOperatingError>>()?;
        let inputs = knots(
            low_drive
                .points()
                .iter()
                .chain(high_drive.points())
                .map(|p| p.flow.value()),
            support.flow,
        )?
        .into_iter()
        .map(|q| {
            let (low_shaft, low_electrical) = low_drive.inputs_at_flow(q)?;
            let (high_shaft, high_electrical) = high_drive.inputs_at_flow(q)?;
            Ok(PumpPowerPoint {
                flow: VolumetricFlowRate::new(q),
                shaft_input: Power::new(blend(low_shaft, high_shaft)),
                electrical_input: Power::new(blend(low_electrical, high_electrical)),
            })
        })
        .collect::<Result<Vec<_>, PumpOperatingError>>()?;
        let provenance = |role: &str, a: &SourceProvenance, b: &SourceProvenance| {
            SourceProvenance::new(
                format!(
                    "{}; {role} at {speed_rpm} rpm: linear same-inlet-flow interpolation between {} rpm [{}: {}] and {} rpm [{}: {}]",
                    self.source.citation,
                    low.speed_rpm,
                    a.identifier,
                    a.citation,
                    high.speed_rpm,
                    b.identifier,
                    b.citation,
                ),
                format!("{}:{role}:{speed_rpm}rpm", self.source.identifier),
            )
        };
        let input = PumpDriveCurve::try_new(
            inputs,
            support.density,
            support.viscosity,
            low_drive.routing(),
            provenance("drive", low_drive.source(), high_drive.source()),
        )?;
        let mut selected = PumpHeadCurve::try_new(
            heads,
            low.curve.gravity_m_per_s2(),
            support.viscosity,
            provenance("head", low.curve.source(), high.curve.source()),
        )?
        .with_drive_curve(input);
        selected.speed_rpm = Some(speed_rpm);
        Ok(selected)
    }
}

fn drive(curve: &PumpHeadCurve) -> Result<&PumpDriveCurve, PumpOperatingError> {
    curve.drive_curve().ok_or(PumpOperatingError::Invalid(
        "every calibrated speed needs independent shaft and electrical input data",
    ))
}

struct CommonSupport {
    flow: (f64, f64),
    density: (f64, f64),
    viscosity: (f64, f64),
}

fn common_support(
    a: &PumpHeadCurve,
    b: &PumpHeadCurve,
) -> Result<CommonSupport, PumpOperatingError> {
    let ad = drive(a)?;
    let bd = drive(b)?;
    let head_span = |curve: &PumpHeadCurve| {
        (
            curve.points()[0].flow.value(),
            curve.points()[curve.points().len() - 1].flow.value(),
        )
    };
    let power_span = |curve: &PumpDriveCurve| {
        (
            curve.points()[0].flow.value(),
            curve.points()[curve.points().len() - 1].flow.value(),
        )
    };
    let flow = intersect(
        intersect(head_span(a), head_span(b))?,
        intersect(power_span(ad), power_span(bd))?,
    )?;
    if flow.0 >= flow.1 {
        return Err(PumpOperatingError::Invalid(
            "adjacent speeds need a nonzero common head/power flow interval",
        ));
    }
    Ok(CommonSupport {
        flow,
        density: intersect(ad.density_domain_kg_m3(), bd.density_domain_kg_m3())?,
        viscosity: intersect(
            intersect(
                a.kinematic_viscosity_domain_m2_s(),
                b.kinematic_viscosity_domain_m2_s(),
            )?,
            intersect(ad.viscosity_domain_m2_s(), bd.viscosity_domain_m2_s())?,
        )?,
    })
}

fn intersect(a: (f64, f64), b: (f64, f64)) -> Result<(f64, f64), PumpOperatingError> {
    let range = (a.0.max(b.0), a.1.min(b.1));
    if range.0 > range.1 {
        Err(PumpOperatingError::Invalid(
            "adjacent speed slices have disjoint flow or material support",
        ))
    } else {
        Ok(range)
    }
}

// Include support endpoints even when neither input has a knot there. Each
// original knot is retained; evenly resampling would erase measured changes
// of slope and alter both the operating point and actual electrical power.
fn knots(
    values: impl Iterator<Item = f64>,
    (low, high): (f64, f64),
) -> Result<Vec<f64>, PumpOperatingError> {
    let mut result = vec![low, high];
    result.extend(values.filter(|&q| q > low && q < high));
    result.sort_by(f64::total_cmp);
    result.dedup();
    if result.len() > 1024 {
        return Err(PumpOperatingError::Invalid(
            "calibrated speed interpolation exceeds the 1024-knot curve bound",
        ));
    }
    Ok(result)
}
