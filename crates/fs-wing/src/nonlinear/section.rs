//! Admitted normal-section coordinates for the nonlinear lifting-line owner.

use crate::{Panel, Refusal, cross3, dot, norm, sub};
use fs_math::det;

const FRAME_TOLERANCE: f64 = 1e-9;

fn refuse(message: &str) -> Refusal {
    Refusal {
        code: "section-frame-invalid",
        message: message.into(),
        ranked_repairs: vec![
            "supply orthonormal forward-chord/span axes aligned with every owning bound row".into(),
        ],
    }
}

/// Coordinates of a section normal to its directed bound-vortex span.
///
/// The chord axis points forward (toward the leading edge), the span axis
/// follows panel A to B, and their cross product is the positive-incidence
/// normal. For the original FRD horizontal section these are +x, +y and +z.
/// Positive section lift opposes this normal, since the native convention is
/// `Gamma = -0.5 * section_speed * chord * cl`.
///
/// Axes are admitted, not silently orthogonalized. `StripSpec::chord_m` must
/// be the physical chord in this normal section: a swept bound segment does
/// not make a streamwise chord interchangeable with the normal-section chord.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SectionFrame {
    chord: [f64; 3],
    span: [f64; 3],
    normal: [f64; 3],
}

/// Accepted local velocity resolved in one section's coordinates.
///
/// `velocity_mps` is the chord/normal-plane projection. The original complete
/// velocity remains the Kutta--Joukowsky force input; this projection supplies
/// the two-dimensional section speed and profile-drag direction.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SectionFlow {
    pub chordwise_mps: f64,
    pub normal_mps: f64,
    pub spanwise_mps: f64,
    pub speed_mps: f64,
    pub incidence_rad: f64,
    pub velocity_mps: [f64; 3],
    incidence_normal: [f64; 3],
}

impl SectionFrame {
    /// The conventional horizontal FRD section.
    pub const BODY_XZ: Self = Self {
        chord: [1.0, 0.0, 0.0],
        span: [0.0, 1.0, 0.0],
        normal: [0.0, 0.0, 1.0],
    };

    /// Admit finite orthogonal unit forward-chord and directed-span axes.
    ///
    /// # Errors
    /// `section-frame-invalid` if either squared norm or their dot product
    /// differs from the orthonormal value by more than 1e-9.
    pub fn new(chord_unit: [f64; 3], span_unit: [f64; 3]) -> Result<Self, Refusal> {
        if !chord_unit
            .iter()
            .chain(span_unit.iter())
            .all(|x| x.is_finite())
            || (dot(chord_unit, chord_unit) - 1.0).abs() > FRAME_TOLERANCE
            || (dot(span_unit, span_unit) - 1.0).abs() > FRAME_TOLERANCE
            || dot(chord_unit, span_unit).abs() > FRAME_TOLERANCE
        {
            return Err(refuse("section axes must be finite, unit and orthogonal"));
        }
        Ok(Self {
            chord: chord_unit,
            span: span_unit,
            normal: cross3(chord_unit, span_unit),
        })
    }

    /// Derive a normal section from a panel's actual span and chord geometry.
    ///
    /// A to B sets the span. The bound midpoint minus the control point sets
    /// the forward direction after removing its spanwise component. Either
    /// sign of the panel boundary-condition normal is legal: reversing that
    /// normal reverses both sides of the linear boundary condition.
    ///
    /// # Errors
    /// `section-frame-invalid` for degenerate/nonfinite geometry or a panel
    /// normal that is not parallel to the derived incidence normal (and thus
    /// perpendicular to the chord/span lifting plane).
    pub fn from_panel(panel: &Panel) -> Result<Self, Refusal> {
        let segment = sub(panel.b, panel.a);
        let length = norm(segment);
        if !length.is_finite() || length <= 0.0 {
            return Err(refuse(
                "a finite nonzero directed bound segment is required",
            ));
        }
        let span = segment.map(|x| x / length);
        let forward = core::array::from_fn(|i| 0.5 * (panel.a[i] + panel.b[i]) - panel.ctrl[i]);
        let swept = dot(forward, span);
        let chord: [f64; 3] = core::array::from_fn(|i| {
            let along_span = swept * span[i];
            forward[i] - along_span
        });
        let length = norm(chord);
        if !length.is_finite() || length <= 0.0 {
            return Err(refuse(
                "the control point must define a nonzero section chord",
            ));
        }
        let frame = Self::new(chord.map(|x| x / length), span)?;
        frame.admit_panel(panel)?;
        Ok(frame)
    }

    #[must_use]
    pub const fn chord_direction(self) -> [f64; 3] {
        self.chord
    }

    #[must_use]
    pub const fn span_direction(self) -> [f64; 3] {
        self.span
    }

    #[must_use]
    pub const fn normal_direction(self) -> [f64; 3] {
        self.normal
    }

    /// Resolve the complete accepted velocity, including already composed wash.
    ///
    /// The span component contributes neither to the two-dimensional speed nor
    /// to incidence. This is a quasi-steady normal-section model; it does not
    /// add an indicial clock, crossflow separation or reverse-flow polar.
    ///
    /// # Errors
    /// `incident-flow-invalid` for nonfinite/unrepresentable flow or a
    /// nonpositive forward chordwise component. Reversal is never hidden by
    /// an absolute value or a change of section orientation.
    pub fn flow(self, velocity: [f64; 3]) -> Result<SectionFlow, Refusal> {
        let chordwise = dot(velocity, self.chord);
        let normal = dot(velocity, self.normal);
        let spanwise = dot(velocity, self.span);
        let speed = norm([chordwise, 0.0, normal]);
        if !velocity.iter().all(|x| x.is_finite())
            || ![chordwise, normal, spanwise, speed]
                .iter()
                .all(|x| x.is_finite())
            || chordwise <= 0.0
            || speed <= 0.0
        {
            return Err(Refusal {
                code: "incident-flow-invalid",
                message: "the section requires finite positive forward chordwise flow".into(),
                ranked_repairs: vec![
                    "use the physical section frame and a forward-flow polar".into(),
                ],
            });
        }
        let section_velocity = core::array::from_fn(|i| {
            let along_chord = chordwise * self.chord[i];
            let along_normal = normal * self.normal[i];
            along_chord + along_normal
        });
        let cosine = chordwise / speed;
        let sine = normal / speed;
        let incidence_normal = core::array::from_fn(|i| {
            let along_normal = cosine * self.normal[i];
            let along_chord = sine * self.chord[i];
            along_normal - along_chord
        });
        Ok(SectionFlow {
            chordwise_mps: chordwise,
            normal_mps: normal,
            spanwise_mps: spanwise,
            speed_mps: speed,
            incidence_rad: det::atan2(normal, chordwise),
            velocity_mps: section_velocity,
            incidence_normal,
        })
    }

    pub(super) fn admit_panel(self, panel: &Panel) -> Result<(), Refusal> {
        if !panel
            .a
            .iter()
            .chain(panel.b.iter())
            .chain(panel.ctrl.iter())
            .chain(panel.normal.iter())
            .all(|x| x.is_finite())
        {
            return Err(refuse("section panel geometry must be finite"));
        }
        let segment = sub(panel.b, panel.a);
        let length = norm(segment);
        if !length.is_finite() || length <= 0.0 {
            return Err(refuse(
                "section panel bound span must be finite and nonzero",
            ));
        }
        if (0..3).any(|i| (segment[i] / length - self.span[i]).abs() > FRAME_TOLERANCE)
            || (dot(panel.normal, panel.normal) - 1.0).abs() > FRAME_TOLERANCE
            || (dot(panel.normal, self.normal).abs() - 1.0).abs() > FRAME_TOLERANCE
        {
            return Err(refuse(
                "all owning rows must share the directed span and chord/span lifting plane",
            ));
        }
        let forward = core::array::from_fn(|i| 0.5 * (panel.a[i] + panel.b[i]) - panel.ctrl[i]);
        let scale = norm(forward);
        if !scale.is_finite()
            || dot(forward, self.chord) <= 0.0
            || dot(forward, self.normal).abs() > FRAME_TOLERANCE * scale
        {
            return Err(refuse(
                "control points must lie aft in their declared chord/span lifting plane",
            ));
        }
        Ok(())
    }
}

impl SectionFlow {
    /// Numerator of the frozen-incident first-order angle response. Dividing
    /// by speed gives (Uc*Wn - Un*Wc)/(Uc^2 + Un^2), the directional derivative
    /// of atan2(Un,Uc). Cache this projection with the trailing-sheet operator.
    pub(super) fn induction_velocity(self, induced: [f64; 3]) -> f64 {
        dot(induced, self.incidence_normal)
    }
}
