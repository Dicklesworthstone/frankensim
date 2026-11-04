//! Supplied coefficient data on the original component force owner.
//!
//! Re = rho*V*L/mu uses the velocity that drives each selected resistance:
//! full local slip for an isotropic component, normal slip for a pressure face.
//! Mach admission always uses the full local slip. Reference definitions:
//! <https://www.grc.nasa.gov/www/k-12/airplane/reynolds.html> and
//! <https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/drag-coefficient/>.
//!
//! Positive piecewise power laws interpolate supplied (Re,Cd) points exactly
//! in log-log coordinates. This is an explicit numerical closure; a citation,
//! interpolation or a successful simulation does not qualify the supplied data
//! for a new geometry, surface state, incidence or Mach range. No default
//! historical coefficient table or extrapolated drag-crisis branch is inferred.

use super::{
    ComponentDragGeometry, ComponentDragShape, DragLedger, DragWrench, Refusal, hash_domain,
};
use fs_math::det;
use std::sync::Arc;

/// Bounded immutable coefficient data per curve.
pub const MAX_REYNOLDS_POINTS: usize = 256;

/// Transport already evaluated by the caller's fluid owner, in SI units.
/// Aircraft adapters obtain all three values from one coherent GasState.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DragAirState {
    /// Mass density [kg/m³].
    pub density_kg_m3: f64,
    /// Dynamic viscosity [Pa·s].
    pub dynamic_viscosity_pa_s: f64,
    /// Local speed of sound [m/s], used for source-domain admission.
    pub sound_speed_mps: f64,
}

impl DragAirState {
    /// Require positive finite transport before evaluating any component flow.
    pub fn admit(self) -> Result<(), Refusal> {
        if [
            self.density_kg_m3,
            self.dynamic_viscosity_pa_s,
            self.sound_speed_mps,
        ]
        .iter()
        .any(|value| !value.is_finite() || *value <= 0.0)
        {
            return Err(refuse(
                "drag-transport-invalid",
                "drag needs positive finite density, dynamic viscosity and sound speed",
            ));
        }
        Ok(())
    }
}

/// One supplied coefficient and its source-declared fractional uncertainty.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReynoldsDragPoint {
    /// Reynolds number based on the bound length and driving velocity.
    pub reynolds: f64,
    /// Drag coefficient for the bound reference-area convention.
    pub cd: f64,
    /// One-sigma fraction in [0,1], as in the existing component ledger.
    pub uncertainty_frac: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct Segment {
    log_re_start: f64,
    inverse_log_span: f64,
    log_cd_start: f64,
    log_cd_change: f64,
}

/// Immutable, source-identified Reynolds coefficient data, prepared once.
/// Evaluation allocates nothing and uses bounded binary search. Exact knots
/// retain the supplied coefficient/uncertainty bits. Interior uncertainty is
/// linear in log(Re); it is not a certified interpolation-error enclosure.
#[derive(Clone, Debug, PartialEq)]
pub struct ReynoldsDragCurve {
    source: Arc<str>,
    maximum_mach: f64,
    points: Arc<[ReynoldsDragPoint]>,
    segments: Arc<[Segment]>,
}

impl ReynoldsDragCurve {
    /// Supply 2..=256 strictly increasing positive Re values, positive Cd,
    /// explicit uncertainty, and the data's admitted subsonic Mach ceiling.
    /// Narrow transition regions require resolved input data; this constructor
    /// cannot discover or certify a missing drag crisis or surface transition.
    pub fn new(
        source: impl Into<String>,
        maximum_mach: f64,
        points: Vec<ReynoldsDragPoint>,
    ) -> Result<Self, Refusal> {
        let source = source.into();
        if source.trim().is_empty()
            || source.len() > 8192
            || !maximum_mach.is_finite()
            || maximum_mach <= 0.0
            || maximum_mach >= 1.0
            || !(2..=MAX_REYNOLDS_POINTS).contains(&points.len())
        {
            return Err(refuse(
                "drag-curve-invalid",
                "coefficient data need a bounded source identity, subsonic Mach ceiling and 2..=256 samples",
            ));
        }
        for (i, point) in points.iter().enumerate() {
            if !point.reynolds.is_finite()
                || point.reynolds <= 0.0
                || !point.cd.is_finite()
                || point.cd <= 0.0
                || !point.uncertainty_frac.is_finite()
                || !(0.0..=1.0).contains(&point.uncertainty_frac)
                || (i > 0 && point.reynolds <= points[i - 1].reynolds)
            {
                return Err(refuse(
                    "drag-curve-invalid",
                    "coefficient samples need ordered positive Re/Cd and uncertainty in [0,1]",
                ));
            }
        }
        let mut segments = Vec::with_capacity(points.len() - 1);
        for pair in points.windows(2) {
            let log_re_start = det::ln(pair[0].reynolds);
            let span = det::ln(pair[1].reynolds) - log_re_start;
            let inverse_log_span = 1.0 / span;
            let log_cd_start = det::ln(pair[0].cd);
            let log_cd_change = det::ln(pair[1].cd) - log_cd_start;
            if !inverse_log_span.is_finite()
                || inverse_log_span <= 0.0
                || !log_re_start.is_finite()
                || !log_cd_start.is_finite()
                || !log_cd_change.is_finite()
            {
                return Err(refuse(
                    "drag-curve-invalid",
                    "coefficient samples are not distinguishable in logarithmic coordinates",
                ));
            }
            segments.push(Segment {
                log_re_start,
                inverse_log_span,
                log_cd_start,
                log_cd_change,
            });
        }
        Ok(Self {
            source: source.into(),
            maximum_mach,
            points: points.into(),
            segments: segments.into(),
        })
    }

    /// Caller-supplied evidence identity for the coefficient data and envelope.
    pub fn source(&self) -> &str {
        &self.source
    }
    /// Immutable, admitted input samples in ascending Reynolds order.
    pub fn points(&self) -> &[ReynoldsDragPoint] {
        &self.points
    }
    /// Largest full-slip Mach number admitted by the supplied source.
    pub const fn maximum_mach(&self) -> f64 {
        self.maximum_mach
    }

    /// Interpolate only inside the supplied closed Reynolds interval.
    /// This coefficient-only query does not admit a physical Mach number;
    /// actual loads additionally enforce it through ReynoldsDragLaw.
    pub fn evaluate(&self, reynolds: f64) -> Result<ReynoldsDragPoint, Refusal> {
        if !reynolds.is_finite()
            || reynolds < self.points[0].reynolds
            || reynolds > self.points[self.points.len() - 1].reynolds
        {
            return Err(refuse(
                "drag-reynolds-outside-domain",
                "local Reynolds number is outside the supplied coefficient interval",
            ));
        }
        let index = self
            .points
            .partition_point(|point| point.reynolds < reynolds);
        if self.points[index].reynolds == reynolds {
            return Ok(self.points[index]);
        }
        let low = self.points[index - 1];
        let high = self.points[index];
        let segment = self.segments[index - 1];
        // Reynolds admission already proved this is interpolation. Bounding
        // the fraction only removes log-coordinate roundoff at a knot.
        let fraction =
            ((det::ln(reynolds) - segment.log_re_start) * segment.inverse_log_span).clamp(0.0, 1.0);
        let cd = if low.cd == high.cd {
            low.cd
        } else {
            det::exp(segment.log_cd_start + fraction * segment.log_cd_change)
        };
        let uncertainty_frac =
            low.uncertainty_frac + fraction * (high.uncertainty_frac - low.uncertainty_frac);
        if !cd.is_finite() || cd <= 0.0 || !uncertainty_frac.is_finite() {
            return Err(refuse(
                "drag-curve-invalid",
                "interpolated coefficient is not representable",
            ));
        }
        Ok(ReynoldsDragPoint {
            reynolds,
            cd,
            uncertainty_frac,
        })
    }

    fn encode(&self, bytes: &mut Vec<u8>) {
        bytes.extend_from_slice(&(self.source.len() as u64).to_le_bytes());
        bytes.extend_from_slice(self.source.as_bytes());
        bytes.extend_from_slice(&self.maximum_mach.to_bits().to_le_bytes());
        bytes.extend_from_slice(&(self.points.len() as u64).to_le_bytes());
        for point in self.points.iter() {
            for value in [point.reynolds, point.cd, point.uncertainty_frac] {
                bytes.extend_from_slice(&value.to_bits().to_le_bytes());
            }
        }
    }
}

/// One physical length definition bound to the coefficient's original area
/// convention. A face law uses normal relative speed for Re; its Mach limit
/// still uses the full relative speed at that material point.
#[derive(Clone, Debug, PartialEq)]
pub struct ReynoldsDragLaw {
    characteristic_length_m: f64,
    curve: ReynoldsDragCurve,
}

impl ReynoldsDragLaw {
    /// Bind prepared data to its positive finite reference length [m].
    pub fn new(characteristic_length_m: f64, curve: ReynoldsDragCurve) -> Result<Self, Refusal> {
        if !characteristic_length_m.is_finite() || characteristic_length_m <= 0.0 {
            return Err(refuse(
                "drag-length-invalid",
                "Reynolds drag needs a positive finite declared characteristic length",
            ));
        }
        Ok(Self {
            characteristic_length_m,
            curve,
        })
    }
    /// Reference length in the coefficient source's Reynolds convention [m].
    pub const fn characteristic_length_m(&self) -> f64 {
        self.characteristic_length_m
    }
    /// The original immutable coefficient source and validity interval.
    pub const fn curve(&self) -> &ReynoldsDragCurve {
        &self.curve
    }

    fn coefficient(
        &self,
        air: DragAirState,
        speed: f64,
        driving_speed: f64,
    ) -> Result<f64, Refusal> {
        let mach = speed / air.sound_speed_mps;
        if !mach.is_finite() || mach > self.curve.maximum_mach {
            return Err(refuse(
                "drag-mach-outside-domain",
                "local full-slip Mach number exceeds the coefficient data's declared ceiling",
            ));
        }
        let reynolds = air.density_kg_m3 * driving_speed * self.characteristic_length_m
            / air.dynamic_viscosity_pa_s;
        Ok(self.curve.evaluate(reynolds)?.cd)
    }

    fn encode(&self, bytes: &mut Vec<u8>) {
        bytes.extend_from_slice(&self.characteristic_length_m.to_bits().to_le_bytes());
        self.curve.encode(bytes);
    }
}

/// Explicit coefficient owner for one component. Each active pressure face
/// may have its own length, source, uncertainty and Reynolds/Mach domain.
#[derive(Clone, Debug, PartialEq)]
pub enum ComponentDragLaw {
    /// Retain the corresponding original ledger coefficient and its source.
    LedgerConstant,
    /// Source data evaluated at the component's full local relative speed.
    Isotropic(ReynoldsDragLaw),
    /// Some exactly for positive-area faces; None exactly for zero-area faces.
    OrthogonalFaces([Option<ReynoldsDragLaw>; 3]),
}

/// A complete identity-bound coefficient choice for one original component.
#[derive(Clone, Debug, PartialEq)]
pub struct ComponentDragCoefficients {
    /// Must match the corresponding ledger and geometry component.
    pub component_id: &'static str,
    /// The explicitly selected constant or Reynolds-dependent coefficient law.
    pub law: ComponentDragLaw,
}

impl DragLedger {
    /// Complete identity/order/shape binding, without fluid sampling or hashes.
    pub fn admit_reynolds_coefficients(
        &self,
        geometry: &[ComponentDragGeometry],
        coefficients: &[ComponentDragCoefficients],
    ) -> Result<(), Refusal> {
        self.admit_geometry(geometry)?;
        if coefficients.len() != self.components.len() {
            return Err(refuse(
                "drag-coefficients-invalid",
                "coefficient selection must cover the complete component ledger",
            ));
        }
        for ((component, geometry), coefficient) in
            self.components.iter().zip(geometry).zip(coefficients)
        {
            let shape_matches = match (&coefficient.law, geometry.shape) {
                (ComponentDragLaw::LedgerConstant, _) => true,
                (ComponentDragLaw::Isotropic(_), ComponentDragShape::Isotropic) => true,
                (
                    ComponentDragLaw::OrthogonalFaces(laws),
                    ComponentDragShape::OrthogonalFaces { areas_m2, .. },
                ) => laws
                    .iter()
                    .zip(areas_m2)
                    .all(|(law, area)| law.is_some() == (area > 0.0)),
                _ => false,
            };
            if component.component_id != coefficient.component_id || !shape_matches {
                return Err(refuse(
                    "drag-coefficients-invalid",
                    "component coefficient identity/order or area model does not match the original geometry",
                ));
            }
        }
        Ok(())
    }

    /// The original spatial force/moment/work assembly with local coefficient
    /// lookup. All selections are admitted before calling the incident field.
    /// Zero driving speed produces zero force without inventing a Cd at Re=0.
    /// No allocation, hashing, curve fitting or material law runs on success.
    pub fn wrench_with_reynolds_geometry(
        &self,
        air: DragAirState,
        geometry: &[ComponentDragGeometry],
        coefficients: &[ComponentDragCoefficients],
        origin_body_m: [f64; 3],
        relative_velocity_at: &dyn Fn([f64; 3]) -> Result<[f64; 3], Refusal>,
    ) -> Result<DragWrench, Refusal> {
        self.admit_reynolds_coefficients(geometry, coefficients)?;
        air.admit()?;
        self.wrench_with_coefficients(
            air.density_kg_m3,
            geometry,
            origin_body_m,
            relative_velocity_at,
            &|index, face, speed, driving_speed| match &coefficients[index].law {
                ComponentDragLaw::LedgerConstant => Ok(self.components[index].cd),
                ComponentDragLaw::Isotropic(law) => law.coefficient(air, speed, driving_speed),
                ComponentDragLaw::OrthogonalFaces(laws) => laws
                    [face.expect("admitted pressure face")]
                .as_ref()
                .expect("admitted positive-area face")
                .coefficient(air, speed, driving_speed),
            },
        )
    }

    /// Identity of the selected geometry, coefficient data, source, uncertainty,
    /// length and Mach envelope. The legacy scalar/spatial digests are unchanged.
    #[must_use]
    pub fn reynolds_digest(
        &self,
        geometry: &[ComponentDragGeometry],
        coefficients: &[ComponentDragCoefficients],
    ) -> String {
        let mut bytes = self.spatial_digest(geometry).into_bytes();
        // The legacy digest predates coefficient-source identity. Include the
        // original source strings in this opt-in identity without changing it.
        bytes.extend_from_slice(&(self.components.len() as u64).to_le_bytes());
        for component in &self.components {
            bytes.extend_from_slice(&(component.cd_source.len() as u64).to_le_bytes());
            bytes.extend_from_slice(component.cd_source.as_bytes());
        }
        bytes.extend_from_slice(&(coefficients.len() as u64).to_le_bytes());
        for coefficient in coefficients {
            bytes.extend_from_slice(&(coefficient.component_id.len() as u64).to_le_bytes());
            bytes.extend_from_slice(coefficient.component_id.as_bytes());
            match &coefficient.law {
                ComponentDragLaw::LedgerConstant => bytes.push(0),
                ComponentDragLaw::Isotropic(law) => {
                    bytes.push(1);
                    law.encode(&mut bytes);
                }
                ComponentDragLaw::OrthogonalFaces(laws) => {
                    bytes.push(2);
                    for law in laws {
                        match law {
                            Some(law) => {
                                bytes.push(1);
                                law.encode(&mut bytes);
                            }
                            None => bytes.push(0),
                        }
                    }
                }
            }
        }
        hash_domain("org.frankensim.fs-flyer.reynolds-drag.v1", &bytes).to_hex()
    }
}

fn refuse(code: &'static str, message: &'static str) -> Refusal {
    Refusal {
        code,
        message: message.into(),
        ranked_repairs: vec![
            "bind matching source data and coherent local gas within their declared domains".into(),
        ],
    }
}
