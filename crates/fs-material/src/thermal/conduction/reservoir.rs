//! One source-free conducting layer and an explicit bulk-side film.
//!
//! The same conductivity chart supplies the exact Kirchhoff integral. A finite
//! film resistance determines the material-face temperature through one bounded
//! monotone solve. No conductivity is frozen at either bulk endpoint.

use super::{ConductivityCurve, LayerConductionConfig};
use fs_thermal::enthalpy::DepositionError;

/// Positive geometric multiplier C in Q = C integral(k(T) dT), in metres.
/// This contains no material coefficient or heat-transfer calibration.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReservoirLayerGeometry {
    multiplier_m: f64,
    surface_area_m2: f64,
}
impl ReservoirLayerGeometry {
    /// Source-free planar layer, with the stated area and actual thickness.
    pub fn planar(area_m2: f64, thickness_m: f64) -> Result<Self, DepositionError> {
        if !positive(area_m2) || !positive(thickness_m) {
            return Err(invalid(
                "conducting reservoir layer needs positive area and thickness",
            ));
        }
        Self::new(area_m2 / thickness_m, area_m2)
    }

    /// Radial annular cylinder with insulated ends; the gas contacts its inner
    /// radius. The outer face is at the declared reservoir temperature.
    pub fn cylindrical(
        inner_radius_m: f64,
        outer_radius_m: f64,
        length_m: f64,
    ) -> Result<Self, DepositionError> {
        if !positive(inner_radius_m)
            || !positive(outer_radius_m)
            || outer_radius_m <= inner_radius_m
            || !positive(length_m)
        {
            return Err(invalid(
                "conducting cylinder needs ordered positive radii and length",
            ));
        }
        let circumference_length = core::f64::consts::TAU * length_m;
        let logarithm = ((outer_radius_m - inner_radius_m) / inner_radius_m).ln_1p();
        Self::new(
            circumference_length / logarithm,
            circumference_length * inner_radius_m,
        )
    }

    fn new(multiplier_m: f64, surface_area_m2: f64) -> Result<Self, DepositionError> {
        if !positive(multiplier_m) || !positive(surface_area_m2) {
            return Err(DepositionError::Arithmetic);
        }
        Ok(Self {
            multiplier_m,
            surface_area_m2,
        })
    }

    /// Actual contacted area; distinct from the geometric resistance multiplier.
    pub const fn surface_area_m2(self) -> f64 {
        self.surface_area_m2
    }

    /// Evaluate Q from the bulk side through a film and the material into a
    /// prescribed reservoir. Film resistance is TOTAL K/W, not areal resistance.
    /// Zero explicitly selects perfect thermal contact. All temperatures must
    /// be covered by the original chart; no extrapolation or hidden film occurs.
    ///
    /// With geometric multiplier C and material face Ts, this solves
    /// Q = C integral from T_reservoir to Ts of k(T)dT and T_bulk - Ts = R_film Q.
    /// This quasi-steady law stores no heat. It does not predict a convective
    /// coefficient, resolve a boundary layer, or account for internal generation.
    pub fn evaluate(
        self,
        curve: &ConductivityCurve,
        temperatures_k: [f64; 2],
        film_resistance_k_w: f64,
        config: LayerConductionConfig,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<ReservoirLayerResponse, DepositionError> {
        checkpoint()?;
        let [bulk, reservoir] = temperatures_k;
        let kb = curve.conductivity_w_m_k(bulk)?;
        let kr = curve.conductivity_w_m_k(reservoir)?;
        if !film_resistance_k_w.is_finite()
            || film_resistance_k_w < 0.0
            || config.max_iterations == 0
            || !config.absolute_heat_tolerance_w.is_finite()
            || config.absolute_heat_tolerance_w < 0.0
            || !config.relative_heat_tolerance.is_finite()
            || config.relative_heat_tolerance < 0.0
        {
            return Err(invalid(
                "invalid conducting reservoir film or solve controls",
            ));
        }
        let delta = bulk - reservoir;
        if delta == 0.0 {
            let g = self.slope(kb, film_resistance_k_w)?;
            checkpoint()?;
            return Ok(ReservoirLayerResponse {
                power_to_reservoir_w: 0.0,
                conductance_w_k: g,
                endpoint_derivative_w_k: [g, -g],
                material_face_temperature_k: bulk,
                heat_residual_w: 0.0,
                iterations: 0,
            });
        }
        let signed_wall_power = self.multiplier_m * curve.integral_w_m(reservoir, bulk)?;
        if !signed_wall_power.is_finite() || signed_wall_power == 0.0 {
            return Err(DepositionError::Arithmetic);
        }
        if film_resistance_k_w == 0.0 {
            let g = signed_wall_power / delta;
            let derivatives = [self.slope(kb, 0.0)?, -self.slope(kr, 0.0)?];
            if !positive(g) {
                return Err(DepositionError::Arithmetic);
            }
            checkpoint()?;
            return Ok(ReservoirLayerResponse {
                power_to_reservoir_w: signed_wall_power,
                conductance_w_k: g,
                endpoint_derivative_w_k: derivatives,
                material_face_temperature_k: bulk,
                heat_residual_w: 0.0,
                iterations: 0,
            });
        }
        let sign = delta.signum();
        let mut low = 0.0;
        // Either the material or film alone gives a rigorous upper flux bound.
        let mut high = signed_wall_power
            .abs()
            .min(delta.abs() / film_resistance_k_w);
        if !positive(high) {
            return Err(DepositionError::Arithmetic);
        }
        // Evaluate the residual across the larger series resistance. When the
        // film dominates, reconstructing Q from the rounded absolute face
        // temperature amplifies its ulp by C*k (large for a thin metal wall).
        // Invert the SAME material integral instead and measure the film flux
        // mismatch. Both residuals increase with Q with derivative >= 1, so
        // their absolute value bounds the remaining flux error in watts.
        let invert_material = film_resistance_k_w >= delta.abs() / signed_wall_power.abs();
        let mut power = 0.5 * high;
        for iteration in 0..config.max_iterations {
            checkpoint()?;
            let (face, kface, other_power, inverse_derivative) = if invert_material {
                let face = curve.temperature_after_integral_w_m(
                    reservoir,
                    sign * (power / self.multiplier_m),
                )?;
                let kface = curve.conductivity_w_m_k(face)?;
                let film_power = sign * (bulk - face) / film_resistance_k_w;
                let material_resistance = (1.0 / kface) / self.multiplier_m;
                (
                    face,
                    kface,
                    film_power,
                    1.0 / (1.0 + material_resistance / film_resistance_k_w),
                )
            } else {
                let face = bulk - sign * (film_resistance_k_w * power);
                let kface = curve.conductivity_w_m_k(face)?;
                let wall_power = sign * self.multiplier_m * curve.integral_w_m(reservoir, face)?;
                (
                    face,
                    kface,
                    wall_power,
                    1.0 / (1.0 + film_resistance_k_w * self.multiplier_m * kface),
                )
            };
            let residual = power - other_power;
            let tolerance = config.absolute_heat_tolerance_w
                + config.relative_heat_tolerance * power.abs().max(other_power.abs());
            if !residual.is_finite() || !tolerance.is_finite() {
                return Err(DepositionError::Arithmetic);
            }
            if residual.abs() <= tolerance {
                let g = power / delta.abs();
                let tangent = self.slope(kface, film_resistance_k_w)?;
                let reservoir_tangent = tangent * (kr / kface);
                if !positive(g) || !positive(reservoir_tangent) {
                    return Err(DepositionError::Arithmetic);
                }
                checkpoint()?;
                return Ok(ReservoirLayerResponse {
                    power_to_reservoir_w: sign * power,
                    conductance_w_k: g,
                    endpoint_derivative_w_k: [tangent, -reservoir_tangent],
                    material_face_temperature_k: face,
                    heat_residual_w: sign * residual,
                    iterations: iteration + 1,
                });
            }
            if residual > 0.0 {
                high = power;
            } else {
                low = power;
            }
            // A safeguarded Newton step; extreme series stiffness falls back to bisection.
            let newton = power - residual * inverse_derivative;
            let next = if newton.is_finite() && low < newton && newton < high {
                newton
            } else {
                low + 0.5 * (high - low)
            };
            if next <= low || next >= high || next == power {
                return Err(DepositionError::Arithmetic);
            }
            power = next;
        }
        Err(DepositionError::Budget(
            "conducting reservoir layer iterations",
        ))
    }

    fn slope(self, conductivity: f64, film: f64) -> Result<f64, DepositionError> {
        let resistance = film + (1.0 / conductivity) / self.multiplier_m;
        let value = 1.0 / resistance;
        if positive(value) {
            Ok(value)
        } else {
            Err(DepositionError::Arithmetic)
        }
    }
}

/// Flux and exact local tangent of the accepted material/film series law.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReservoirLayerResponse {
    /// Signed heat leaving the bulk and entering the prescribed reservoir [W].
    pub power_to_reservoir_w: f64,
    /// Positive Q/(T_bulk - T_reservoir), with its equilibrium limit [W/K].
    pub conductance_w_k: f64,
    /// dQ/d[T_bulk,T_reservoir] [W/K], at fixed material and geometry.
    pub endpoint_derivative_w_k: [f64; 2],
    /// Material temperature after the gas/contact film [K].
    pub material_face_temperature_k: f64,
    /// Signed residual of the selected monotone flux equation [W]. Its absolute
    /// value bounds the flux error; either the film or the original material
    /// integral supplies the inverse law to avoid a stiff temperature residual.
    pub heat_residual_w: f64,
    /// Root iterations; zero for the perfect-contact or equilibrium limit.
    pub iterations: usize,
}

fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
fn invalid(reason: &'static str) -> DepositionError {
    DepositionError::InvalidInput(reason)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn config() -> LayerConductionConfig {
        LayerConductionConfig {
            max_iterations: 96,
            absolute_heat_tolerance_w: 1e-10,
            relative_heat_tolerance: 1e-12,
        }
    }

    fn close(actual: f64, expected: f64, tolerance: f64) {
        assert!(
            (actual - expected).abs() <= tolerance,
            "{actual} differs from {expected} by more than {tolerance}"
        );
    }

    #[test]
    fn constant_source_reproduces_planar_and_cylindrical_series_resistance() {
        let conductivity = ConductivityCurve::new(vec![[200.0, 4.0], [1000.0, 4.0]]).unwrap();
        let planar = ReservoirLayerGeometry::planar(0.02, 0.04).unwrap();
        let cylinder = ReservoirLayerGeometry::cylindrical(0.02, 0.03, 0.1).unwrap();
        for (geometry, material_resistance) in [
            (planar, 0.04 / (0.02 * 4.0)),
            (
                cylinder,
                (1.5_f64).ln() / (core::f64::consts::TAU * 0.1 * 4.0),
            ),
        ] {
            for film in [0.0, 0.75] {
                let expected_g = 1.0 / (material_resistance + film);
                for temperatures in [[700.0, 300.0], [300.0, 700.0], [450.0, 450.0]] {
                    let response = geometry
                        .evaluate(&conductivity, temperatures, film, config(), || Ok(()))
                        .unwrap();
                    let expected_q = expected_g * (temperatures[0] - temperatures[1]);
                    close(response.power_to_reservoir_w, expected_q, 1e-9);
                    close(response.conductance_w_k, expected_g, 1e-11);
                    close(response.endpoint_derivative_w_k[0], expected_g, 1e-11);
                    close(response.endpoint_derivative_w_k[1], -expected_g, 1e-11);
                    close(
                        response.material_face_temperature_k,
                        temperatures[0] - film * expected_q,
                        1e-9,
                    );
                }
            }
        }
        close(
            cylinder.surface_area_m2(),
            core::f64::consts::TAU * 0.02 * 0.1,
            1e-15,
        );
    }

    #[test]
    fn variable_conductivity_integral_and_film_recover_manufactured_faces_and_tangents() {
        // Synthetic k(T)=0.01*T, with a crossed source knot. These exact Fourier
        // solutions are independent of the numerical root finder.
        let conductivity =
            ConductivityCurve::new(vec![[200.0, 2.0], [350.0, 3.5], [1000.0, 10.0]]).unwrap();
        let geometry = ReservoirLayerGeometry::planar(0.02, 0.04).unwrap();
        for (bulk, reservoir, face, expected_q) in
            [(580.0, 300.0, 500.0, 400.0), (220.0, 500.0, 300.0, -400.0)]
        {
            let response = geometry
                .evaluate(&conductivity, [bulk, reservoir], 0.2, config(), || Ok(()))
                .unwrap();
            close(response.power_to_reservoir_w, expected_q, 1e-8);
            close(response.material_face_temperature_k, face, 1e-8);
            close(response.heat_residual_w, 0.0, 1e-8);
            let denominator = 1.0 + 0.2 * 0.5 * 0.01 * face;
            close(
                response.endpoint_derivative_w_k[0],
                0.5 * 0.01 * face / denominator,
                1e-10,
            );
            close(
                response.endpoint_derivative_w_k[1],
                -0.5 * 0.01 * reservoir / denominator,
                1e-10,
            );
            for endpoint in 0..2 {
                let mut plus = [bulk, reservoir];
                let mut minus = plus;
                plus[endpoint] += 1e-3;
                minus[endpoint] -= 1e-3;
                let q_plus = geometry
                    .evaluate(&conductivity, plus, 0.2, config(), || Ok(()))
                    .unwrap()
                    .power_to_reservoir_w;
                let q_minus = geometry
                    .evaluate(&conductivity, minus, 0.2, config(), || Ok(()))
                    .unwrap()
                    .power_to_reservoir_w;
                close(
                    (q_plus - q_minus) / 2e-3,
                    response.endpoint_derivative_w_k[endpoint],
                    2e-7,
                );
            }
        }
    }

    #[test]
    fn film_dominated_metal_wall_retains_small_heat_fluxes_near_equilibrium() {
        let conductivity = ConductivityCurve::new(vec![[200.0, 200.0], [1000.0, 200.0]]).unwrap();
        let geometry = ReservoirLayerGeometry::planar(0.1, 0.0001).unwrap();
        let film = 20.0;
        let expected_g = 1.0 / (film + 1.0 / (geometry.multiplier_m * 200.0));
        let controls = LayerConductionConfig {
            absolute_heat_tolerance_w: 1e-12,
            ..config()
        };
        for bulk in [300.003, 299.997, 320.0] {
            let response = geometry
                .evaluate(&conductivity, [bulk, 300.0], film, controls, || Ok(()))
                .unwrap();
            close(
                response.power_to_reservoir_w,
                expected_g * (bulk - 300.0),
                1e-12,
            );
            assert!(
                response.heat_residual_w.abs()
                    <= controls.absolute_heat_tolerance_w
                        + controls.relative_heat_tolerance * response.power_to_reservoir_w.abs()
            );
        }
    }

    #[test]
    fn aluminum_source_curve_matches_analytic_film_flux_without_temperature_quantization() {
        // Original k(T) knots from aluminum-2024-t3-first-heating-km2024.
        // The real gas-field example exposed stagnation from C*k*ulp(300 K).
        let conductivity =
            ConductivityCurve::new(vec![[298.15, 121.1], [323.15, 125.2], [373.15, 133.0]])
                .unwrap();
        let geometry = ReservoirLayerGeometry::cylindrical(0.0073, 0.0103, 0.285).unwrap();
        let film = 0.01 / geometry.surface_area_m2();
        let controls = LayerConductionConfig {
            max_iterations: 96,
            absolute_heat_tolerance_w: 1e-11,
            relative_heat_tolerance: 1e-12,
        };
        let slope = (125.2 - 121.1) / (323.15 - 298.15);
        let k_reservoir = 121.1 + slope * (300.0 - 298.15);
        // Within the face's first affine source segment, u = T_face - 300
        // satisfies a*u^2 + b*u - delta = 0. This is an independent closed form.
        let a = 0.5 * film * geometry.multiplier_m * slope;
        let b = 1.0 + film * geometry.multiplier_m * k_reservoir;
        for bulk in [300.003, 302.97, 320.0, 340.0, 299.0] {
            let delta = bulk - 300.0;
            let offset = 2.0 * delta / (b + (b * b + 4.0 * a * delta).sqrt());
            let expected_q = (delta - offset) / film;
            let response = geometry
                .evaluate(&conductivity, [bulk, 300.0], film, controls, || Ok(()))
                .unwrap();
            assert!((298.15..323.15).contains(&response.material_face_temperature_k));
            close(response.power_to_reservoir_w, expected_q, 1e-11);
            close(response.material_face_temperature_k, 300.0 + offset, 1e-11);
        }
    }

    #[test]
    fn source_domain_bad_geometry_work_budget_and_cancellation_refuse() {
        let conductivity = ConductivityCurve::new(vec![[200.0, 2.0], [1000.0, 10.0]]).unwrap();
        let geometry = ReservoirLayerGeometry::planar(0.02, 0.04).unwrap();
        assert!(ReservoirLayerGeometry::planar(0.0, 0.04).is_err());
        assert!(ReservoirLayerGeometry::cylindrical(0.02, 0.01, 0.1).is_err());
        assert!(matches!(
            geometry.evaluate(&conductivity, [1001.0, 300.0], 0.2, config(), || Ok(())),
            Err(DepositionError::MaterialDomain)
        ));
        assert!(
            geometry
                .evaluate(&conductivity, [580.0, 300.0], -1.0, config(), || Ok(()))
                .is_err()
        );
        assert!(matches!(
            geometry.evaluate(
                &conductivity,
                [580.0, 300.0],
                0.2,
                LayerConductionConfig {
                    max_iterations: 1,
                    ..config()
                },
                || Ok(())
            ),
            Err(DepositionError::Budget(_))
        ));
        assert!(matches!(
            geometry.evaluate(&conductivity, [580.0, 300.0], 0.2, config(), || Err(
                DepositionError::Cancelled
            )),
            Err(DepositionError::Cancelled)
        ));
    }
}
