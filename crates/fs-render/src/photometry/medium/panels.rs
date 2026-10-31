//! Metric panel coordinates that retain the complete optical interval.
//!
//! Thin intervals keep the original uniform arithmetic. Thick intervals use
//! t(u) = expm1(u * ln(1 + sigma * L)) / sigma. This resolves the foreground
//! attenuation scale without truncating distant incident-field support.

use fs_math::det;

#[derive(Clone, Copy, Debug)]
pub(super) struct OpticalPanels {
    distance_m: f64,
    extinction_per_m: f64,
    log_span: f64,
    log_extinction: f64,
}

#[derive(Clone, Copy, Debug)]
pub(super) struct Panel {
    pub start_m: f64,
    pub width_m: f64,
    pub midpoint_m: f64,
}

impl OpticalPanels {
    pub(super) fn new(extinction_per_m: f64, distance_m: f64) -> Option<Self> {
        if !extinction_per_m.is_finite()
            || extinction_per_m < 0.0
            || !distance_m.is_finite()
            || distance_m < 0.0
        {
            return None;
        }
        let depth = extinction_per_m * distance_m;
        let log_span = if depth <= 1.0 {
            0.0
        } else if depth.is_finite() {
            det::ln(1.0 + depth)
        } else {
            // Both inputs are finite; only their optical-depth product
            // overflowed. At this scale adding one cannot affect the log.
            det::ln(extinction_per_m) + det::ln(distance_m)
        };
        Some(Self {
            distance_m,
            extinction_per_m,
            log_span,
            log_extinction: if log_span > 0.0 {
                det::ln(extinction_per_m)
            } else {
                0.0
            },
        })
    }

    pub(super) const fn is_uniform(self) -> bool {
        self.log_span == 0.0
    }

    /// The first panel must span at most one optical depth before successive
    /// small discrepancies can qualify as resolved foreground quadratures.
    /// This prerequisite does not bound arbitrary incident-field variation.
    pub(super) fn optical_scale_resolved(self, count: u32) -> bool {
        count > 0
            && (self.is_uniform() || self.log_span <= f64::from(count) * core::f64::consts::LN_2)
    }

    pub(super) fn panel(self, index: u32, count: u32) -> Option<Panel> {
        if count == 0 || index >= count {
            return None;
        }
        if self.is_uniform() {
            // Preserve all three original operations, including midpoint
            // rounding; start + width/2 need not produce the same bits.
            let width_m = self.distance_m / f64::from(count);
            return (width_m > 0.0).then_some(Panel {
                start_m: f64::from(index) * width_m,
                width_m,
                midpoint_m: (f64::from(index) + 0.5) * width_m,
            });
        }
        let start_m = self.boundary(index, count)?;
        let end_m = self.boundary(index + 1, count)?;
        let width_m = end_m - start_m;
        let midpoint_m = start_m + 0.5 * width_m;
        (width_m > 0.0 && midpoint_m > start_m && midpoint_m < end_m).then_some(Panel {
            start_m,
            width_m,
            midpoint_m,
        })
    }

    fn boundary(self, index: u32, count: u32) -> Option<f64> {
        if count == 0 || index > count {
            return None;
        }
        if index == 0 {
            return Some(0.0);
        }
        if index == count {
            return Some(self.distance_m);
        }
        let exponent = self.log_span * (f64::from(index) / f64::from(count));
        let optical_distance = det::expm1(exponent);
        let distance = if optical_distance.is_finite() {
            optical_distance / self.extinction_per_m
        } else {
            // expm1(exponent) may overflow even though division by sigma
            // would restore a finite metric boundary. The subtracted one
            // is negligible here; keep that quotient in the log domain.
            det::exp(exponent - self.log_extinction)
        };
        (distance.is_finite() && distance > 0.0 && distance < self.distance_m).then_some(distance)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn integrate(
        map: OpticalPanels,
        count: u32,
        albedo: f64,
        incident: impl Fn(f64) -> f64,
    ) -> f64 {
        let mut sum = 0.0;
        for index in 0..count {
            let panel = map.panel(index, count).unwrap();
            let weight = albedo * -det::expm1(-map.extinction_per_m * panel.width_m);
            sum += (incident(panel.midpoint_m) * weight)
                * det::exp(-map.extinction_per_m * panel.start_m);
        }
        sum
    }

    #[test]
    fn g0_thin_panels_preserve_original_coordinate_bits() {
        for (sigma, distance, count) in [(0.3, 1.0, 7), (1.0, 1.0, 32), (0.0, 3.0, 19)] {
            let map = OpticalPanels::new(sigma, distance).unwrap();
            assert!(map.is_uniform());
            assert!(map.optical_scale_resolved(count));
            let width = distance / f64::from(count);
            for index in 0..count {
                let panel = map.panel(index, count).unwrap();
                assert_eq!(panel.width_m.to_bits(), width.to_bits());
                assert_eq!(
                    panel.start_m.to_bits(),
                    (f64::from(index) * width).to_bits()
                );
                assert_eq!(
                    panel.midpoint_m.to_bits(),
                    ((f64::from(index) + 0.5) * width).to_bits()
                );
            }
        }
    }

    #[test]
    fn g1_thick_panels_cover_the_interval_and_resolve_its_optical_scale() {
        let map = OpticalPanels::new(100.0, 100.0).unwrap();
        assert!(!map.optical_scale_resolved(8));
        assert!(map.optical_scale_resolved(16));
        assert!(100.0 * map.panel(0, 8).unwrap().width_m > 1.0);
        assert!(100.0 * map.panel(0, 16).unwrap().width_m < 1.0);
        assert_eq!(map.boundary(0, 32), Some(0.0));
        assert_eq!(map.boundary(32, 32), Some(100.0));
        let mut previous = 0.0;
        for index in 1..=32 {
            let boundary = map.boundary(index, 32).unwrap();
            let expected = ((10001.0_f64).ln() * f64::from(index) / 32.0).exp_m1() / 100.0;
            assert!(boundary > previous);
            assert!((boundary / expected - 1.0).abs() < 3e-14);
            previous = boundary;
        }
        // Equal optical depth under an SI length/coefficient rescaling.
        let scaled = OpticalPanels::new(0.1, 100_000.0).unwrap();
        for index in 0..32 {
            let a = map.panel(index, 32).unwrap();
            let b = scaled.panel(index, 32).unwrap();
            assert!((b.midpoint_m / (1000.0 * a.midpoint_m) - 1.0).abs() < 3e-14);
        }
    }

    #[test]
    fn g1_overflowing_optical_depth_retains_finite_metric_panels() {
        let map = OpticalPanels::new(1e308, 1e308).unwrap();
        assert!(map.optical_scale_resolved(4096));
        assert!(!map.optical_scale_resolved(1024));
        for index in [0, 1, 2047, 4094, 4095] {
            let panel = map.panel(index, 4096).unwrap();
            assert!(panel.width_m.is_finite() && panel.width_m > 0.0);
            assert!(panel.midpoint_m > panel.start_m && panel.midpoint_m < 1e308);
        }
        assert_eq!(map.boundary(4096, 4096), Some(1e308));
        // Even the first panel at the largest admitted count is representable.
        assert!(map.panel(0, u32::MAX).unwrap().width_m > 0.0);
    }

    #[test]
    fn g1_nonuniform_beer_weights_preserve_a_constant_incident_field() {
        for (sigma, distance) in [(100.0, 100.0), (0.2, 1000.0), (1e308, 1e308)] {
            let map = OpticalPanels::new(sigma, distance).unwrap();
            let sum = integrate(map, 4096, 1.0, |_| 1.0);
            let expected = -(-sigma * distance).exp_m1();
            assert!((sum - expected).abs() < 3e-13, "{sum} versus {expected}");
        }
    }

    #[test]
    fn g1_thick_source_matches_independent_dimensionless_integral() {
        // I=100 cd, source at z=-0.01 m, camera z=0..100 m, sigma_t=100/m,
        // sigma_s=80/m. Substitute u=t/0.01 in the physical source integral;
        // reference quadrature uses only standard exp and uniform Simpson.
        let count = 8192;
        let step = 40.0 / f64::from(count);
        let mut reference = 0.0;
        for index in 0..=count {
            let u = f64::from(index) * step;
            let weight = if index == 0 || index == count {
                1.0
            } else if index % 2 == 0 {
                2.0
            } else {
                4.0
            };
            reference += weight * (-2.0 * u).exp() / (1.0 + u).powi(2);
        }
        reference *=
            step / 3.0 * (100.0 * 80.0 / (4.0 * core::f64::consts::PI * 0.01)) * (-1.0_f64).exp();
        // The omitted u>40 tail is below 1.3e-34 cd/m².
        let map = OpticalPanels::new(100.0, 100.0).unwrap();
        let incident = |t: f64| {
            let r = 0.01 + t;
            100.0 * det::exp(-100.0 * r) / (4.0 * core::f64::consts::PI * r * r)
        };
        let coarse = integrate(map, 32, 0.8, incident);
        let refined = integrate(map, 64, 0.8, incident);
        assert!(coarse > 0.9 * reference && coarse < refined && refined < reference);
        let fine = integrate(map, 2048, 0.8, incident);
        assert!(
            (fine / reference - 1.0).abs() < 2e-5,
            "{fine} versus {reference}"
        );
    }

    #[test]
    fn g1_rising_incident_field_retains_distant_support() {
        // The incident exp(t) cancels the camera exp(-t), giving exactly 100.
        let map = OpticalPanels::new(1.0, 100.0).unwrap();
        let actual = integrate(map, 4096, 1.0, f64::exp);
        assert!((actual - 100.0).abs() < 0.02, "{actual}");
    }

    #[test]
    fn g0_zero_lengths_and_invalid_coordinates_do_not_create_panels() {
        let zero = OpticalPanels::new(100.0, 0.0).unwrap();
        assert!(zero.is_uniform());
        assert!(zero.panel(0, 1).is_none());
        assert!(!zero.optical_scale_resolved(0));
        assert!(zero.panel(0, 0).is_none());
        assert!(zero.panel(1, 1).is_none());
        for invalid in [-1.0, f64::NAN, f64::INFINITY] {
            assert!(OpticalPanels::new(invalid, 1.0).is_none());
            assert!(OpticalPanels::new(1.0, invalid).is_none());
        }
    }
}
