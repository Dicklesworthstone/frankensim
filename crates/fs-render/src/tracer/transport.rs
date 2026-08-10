//! Reciprocal transport bookkeeping shared by camera and light subpaths.
//!
//! The legacy path tracer only needed the density of the direction it had
//! just sampled. Bidirectional techniques additionally need the density of the
//! same edge when the path is generated in reverse, and both densities must be
//! expressed in a common measure before MIS. This module owns those semantics;
//! it deliberately contains no scene- or Euler-specific policy.

/// Direction in which a scattering function transports throughput.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum TransportMode {
    /// Sensor-to-scene paths carry radiance and include the refractive
    /// `eta_i^2 / eta_t^2` Jacobian on transmission.
    Radiance,
    /// Emitter-to-scene paths carry importance. The corresponding refractive
    /// Jacobian is supplied when the light and camera techniques are compared.
    Importance,
}

/// Forward and reverse densities of one sampled scattering event.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct DirectionalPdfPair {
    forward_solid_angle: f64,
    reverse_solid_angle: f64,
    delta: bool,
}

impl DirectionalPdfPair {
    /// Admit one pair of finite solid-angle densities.
    pub(super) fn continuous(
        forward_solid_angle: f64,
        reverse_solid_angle: f64,
    ) -> Option<Self> {
        if valid_pdf(forward_solid_angle) && valid_pdf(reverse_solid_angle) {
            Some(Self {
                forward_solid_angle,
                reverse_solid_angle,
                delta: false,
            })
        } else {
            None
        }
    }

    /// A delta event has no finite solid-angle density. Its discrete branch
    /// probability belongs to path throughput, not an area-measure MIS PDF.
    pub(super) const fn delta() -> Self {
        Self {
            forward_solid_angle: 0.0,
            reverse_solid_angle: 0.0,
            delta: true,
        }
    }

    pub(super) const fn forward_solid_angle(self) -> f64 {
        self.forward_solid_angle
    }

    pub(super) const fn reverse_solid_angle(self) -> f64 {
        self.reverse_solid_angle
    }

    pub(super) const fn is_delta(self) -> bool {
        self.delta
    }

    /// Convert both directions to the area measures at their respective edge
    /// targets. For an edge `x -> y`, the forward density uses the cosine at
    /// `y`; the reverse density uses the cosine at `x`.
    pub(super) fn to_area(
        self,
        distance_squared: f64,
        target_abs_cosine: f64,
        source_abs_cosine: f64,
    ) -> Option<AreaPdfPair> {
        if self.delta {
            return Some(AreaPdfPair::delta());
        }
        if !distance_squared.is_finite()
            || distance_squared <= 0.0
            || !valid_cosine(target_abs_cosine)
            || !valid_cosine(source_abs_cosine)
        {
            return None;
        }
        let forward_area = self.forward_solid_angle * target_abs_cosine / distance_squared;
        let reverse_area = self.reverse_solid_angle * source_abs_cosine / distance_squared;
        if valid_pdf(forward_area) && valid_pdf(reverse_area) {
            Some(AreaPdfPair {
                forward_area,
                reverse_area,
                delta: false,
            })
        } else {
            None
        }
    }
}

/// Forward and reverse densities after solid-angle-to-area conversion.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct AreaPdfPair {
    forward_area: f64,
    reverse_area: f64,
    delta: bool,
}

impl AreaPdfPair {
    const fn delta() -> Self {
        Self {
            forward_area: 0.0,
            reverse_area: 0.0,
            delta: true,
        }
    }

    pub(super) const fn forward_area(self) -> f64 {
        self.forward_area
    }

    pub(super) const fn reverse_area(self) -> f64 {
        self.reverse_area
    }

    pub(super) const fn is_delta(self) -> bool {
        self.delta
    }
}

/// Refractive throughput factor for the declared transport direction.
pub(super) fn refractive_transport_factor(
    mode: TransportMode,
    eta_incident: f64,
    eta_transmitted: f64,
) -> Option<f64> {
    if !eta_incident.is_finite()
        || !eta_transmitted.is_finite()
        || eta_incident <= 0.0
        || eta_transmitted <= 0.0
    {
        return None;
    }
    let factor = match mode {
        TransportMode::Radiance => {
            let ratio = eta_incident / eta_transmitted;
            ratio * ratio
        }
        TransportMode::Importance => 1.0,
    };
    factor.is_finite().then_some(factor)
}

fn valid_pdf(pdf: f64) -> bool {
    pdf.is_finite() && pdf >= 0.0
}

fn valid_cosine(cosine: f64) -> bool {
    cosine.is_finite() && (0.0..=1.0).contains(&cosine)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn g0_area_measure_conversion_is_reciprocal_under_edge_reversal() {
        let pair = DirectionalPdfPair::continuous(0.25, 0.75).unwrap();
        let forward = pair.to_area(4.0, 0.5, 0.8).unwrap();
        let reverse = DirectionalPdfPair::continuous(0.75, 0.25)
            .unwrap()
            .to_area(4.0, 0.8, 0.5)
            .unwrap();
        assert_eq!(
            forward.forward_area().to_bits(),
            reverse.reverse_area().to_bits()
        );
        assert_eq!(
            forward.reverse_area().to_bits(),
            reverse.forward_area().to_bits()
        );
        assert!(!forward.is_delta());
    }

    #[test]
    fn g0_delta_event_never_masquerades_as_finite_area_density() {
        let area = DirectionalPdfPair::delta()
            .to_area(2.0, 0.4, 0.6)
            .unwrap();
        assert!(area.is_delta());
        assert_eq!(area.forward_area().to_bits(), 0.0_f64.to_bits());
        assert_eq!(area.reverse_area().to_bits(), 0.0_f64.to_bits());
    }

    #[test]
    fn g0_radiance_eta_factors_cancel_across_a_lossless_slab() {
        let entry = refractive_transport_factor(TransportMode::Radiance, 1.0, 1.5).unwrap();
        let exit = refractive_transport_factor(TransportMode::Radiance, 1.5, 1.0).unwrap();
        assert!((entry * exit - 1.0).abs() <= 2.0 * f64::EPSILON);
        assert_eq!(
            refractive_transport_factor(TransportMode::Importance, 1.0, 1.5)
                .unwrap()
                .to_bits(),
            1.0_f64.to_bits()
        );
    }
}
