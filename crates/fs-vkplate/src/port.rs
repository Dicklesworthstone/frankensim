//! A compact, low-Mach vent with inertance and optional quasi-steady flow separation.
//!
//! The declared loss coefficient K defines the aperture pressure drop
//! `dp = R Q + rho K (Q/A) |Q/A| / 2`. No universal K is inferred from a
//! hole radius: edge shape, flow regime and the pressure measurement locations
//! matter. This is a lumped closure, not resolved vortex shedding or jet noise.
//! See Hersh & Rogers, NASA CR-2682 (1976), sections 2.1.1 and 3:
//! <https://ntrs.nasa.gov/citations/19760016859>.

use std::f64::consts::PI;

mod radiation;
use radiation::BaffledRadiation;

/// A vent's air plug: volume inertance `rho L_eff / A`, driven by cavity pressure.
/// Volume flow is positive outward. Geometry, density and clock must be finite
/// and positive; linear resistance must be finite and nonnegative.
#[derive(Clone)]
pub struct Port {
    pub area: f64,
    /// Throat length plus the end corrections still represented by this plug [m].
    pub length: f64,
    /// Included interior end correction [m]; resolved cavity modes replace it.
    pub interior_end_correction: f64,
    /// Included exterior end correction [m]; the causal radiation load replaces it.
    pub exterior_end_correction: f64,
    /// Linear radiation/flow resistance [Pa s / m^3].
    pub resistance: f64,
    flow_loss_coefficient: f64,
    radiation: Option<BaffledRadiation>,
    /// Integrated work removed by the linear and quadratic pressure drops [J].
    /// This includes linear radiation resistance, so it is not all local heat.
    pub dissipated_j: f64,
    /// Opening position and outward normal sign along z (for microphones).
    pub pos: [f64; 3],
    pub facing: f64,
    /// Displaced volume at the integer step [m^3].
    pub x: f64,
    /// Volume flow at the half step [m^3/s].
    pub u: f64,
    /// Volume-flow acceleration at the integer step [m^3/s^2].
    pub du: f64,
}

impl Port {
    pub fn new(radius: f64, head_thickness: f64, pos: [f64; 3], facing: f64) -> Port {
        let area = PI * radius * radius;
        let interior_end_correction = 0.85 * radius;
        let length = head_thickness + 2.0 * interior_end_correction;
        Port {
            area,
            length,
            interior_end_correction,
            exterior_end_correction: interior_end_correction,
            resistance: 0.0,
            flow_loss_coefficient: 0.0,
            radiation: None,
            dissipated_j: 0.0,
            pos,
            facing,
            x: 0.0,
            u: 0.0,
            du: 0.0,
        }
    }

    /// Select a dimensionless quasi-steady loss coefficient K. Zero disables
    /// the quadratic term. Supply a coefficient for this aperture and regime;
    /// do not add it on top of a resolved model of the same separated flow.
    pub fn set_flow_loss_coefficient(&mut self, coefficient: f64) -> Result<(), String> {
        if !coefficient.is_finite() || coefficient < 0.0 {
            return Err("vent flow-loss coefficient must be finite and nonnegative".into());
        }
        self.flow_loss_coefficient = coefficient;
        Ok(())
    }

    pub fn flow_loss_coefficient(&self) -> f64 {
        self.flow_loss_coefficient
    }

    /// Attach the uniform circular-piston radiation impulse response in an
    /// infinite rigid baffle. Attach at rest, on the eventual integration clock.
    /// It replaces the included exterior end correction; the remaining linear
    /// resistance is retained and must exclude that exterior radiation.
    ///
    /// The finite memory is the aperture transit time rounded up to this clock
    /// (at most 4096 taps). It is passive for every admitted clock; waveform and
    /// impedance accuracy still require time refinement. Aperture flow remains
    /// uniform and low-Mach. Shell diffraction and mutual exterior loading are
    /// outside this isolated baffled-aperture model.
    pub fn enable_baffled_radiation(
        &mut self,
        density: f64,
        sound_speed: f64,
        dt: f64,
    ) -> Result<(), String> {
        if self.radiation.is_some()
            || self.x != 0.0
            || self.u != 0.0
            || self.du != 0.0
            || self.dissipated_j != 0.0
        {
            return Err(
                "attach baffled vent radiation once, with the port initially at rest".into(),
            );
        }
        let length = self.length - self.exterior_end_correction;
        if !self.exterior_end_correction.is_finite()
            || self.exterior_end_correction < 0.0
            || !length.is_finite()
            || length <= 0.0
        {
            return Err("baffled vent radiation needs positive throat/interior inertance after removing the exterior correction".into());
        }
        let radiation = BaffledRadiation::new(self.area, density, sound_speed, dt)?;
        self.radiation = Some(radiation);
        self.length = length;
        self.exterior_end_correction = 0.0;
        Ok(())
    }

    pub fn has_baffled_radiation(&self) -> bool {
        self.radiation.is_some()
    }

    /// Energy stored in the causal exterior load's delay memory [J].
    pub fn radiation_energy_j(&self) -> f64 {
        self.radiation.as_ref().map_or(0.0, |r| r.stored_j)
    }

    /// Energy emitted by the discrete passive radiation load [J], separate
    /// from neck/separation work in `dissipated_j` and from stored exterior energy.
    pub fn radiated_j(&self) -> f64 {
        self.radiation.as_ref().map_or(0.0, |r| r.radiated_j)
    }

    /// Leapfrog inertia with midpoint linear/quadratic drag. All coupled
    /// pressures belong to the old integer state, as do the head/modal forces.
    /// The scalar solve is monotone and closed form: no iteration or allocation.
    pub(crate) fn advance_pressure(&mut self, dt: f64, pressure: f64, density: f64) {
        assert!(
            self.resistance.is_finite() && self.resistance >= 0.0,
            "vent linear resistance must be finite and nonnegative"
        );
        let (radiation_gain, radiation_memory) = self.radiation.as_ref().map_or((0.0, 0.0), |r| {
            assert!(
                r.dt == dt && r.area == self.area && r.density == density,
                "baffled vent radiation must retain its prepared clock, aperture and gas density"
            );
            r.pressure_coefficients()
        });
        // Solve in particle velocity to avoid rho K / (2 A^2) at small apertures.
        // Radiation's current-flow term is implicit in this same scalar solve.
        let h = 0.5 * dt / (density * self.length);
        let a = 1.0 + h * (self.resistance + radiation_gain) * self.area;
        let b = 0.25 * dt * self.flow_loss_coefficient / self.length;
        let s = self.u / self.area + h * (pressure + radiation_memory);
        assert!(
            a.is_finite() && a >= 1.0 && b.is_finite() && b >= 0.0 && s.is_finite(),
            "vent pressure step is outside its finite passive domain"
        );
        let mean_velocity = if b == 0.0 || s == 0.0 {
            s / a
        } else {
            // Positive root of a |v| + b |v|^2 = |s|, scaled before
            // squaring and rationalized to avoid subtracting nearly equal roots.
            let t = b.sqrt() * s.abs().sqrt();
            let scale = a.max(t);
            let d = 0.5 * (a / scale);
            (s / scale) / (d + d.hypot(t / scale))
        };
        let mean_flow = self.area * mean_velocity;
        let next_flow = 2.0 * mean_flow - self.u;
        let loss_pressure = self.resistance * mean_flow
            + 0.5 * density * self.flow_loss_coefficient * mean_velocity * mean_velocity.abs();
        self.dissipated_j += dt * mean_flow * loss_pressure;
        if let Some(radiation) = &mut self.radiation {
            radiation.advance(mean_flow);
        }
        self.du = (next_flow - self.u) / dt;
        self.u = next_flow;
        self.x += dt * self.u;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn port() -> Port {
        Port::new(0.025, 0.0003, [0.0, 0.0, -0.15], -1.0)
    }

    #[test]
    fn coefficient_refusal_preserves_the_admitted_law() {
        let mut p = port();
        assert_eq!(p.flow_loss_coefficient(), 0.0);
        p.set_flow_loss_coefficient(1.7).unwrap();
        for invalid in [-1.0, f64::NAN, f64::INFINITY] {
            assert!(p.set_flow_loss_coefficient(invalid).is_err());
            assert_eq!(p.flow_loss_coefficient(), 1.7);
        }
    }

    #[test]
    fn steady_pressure_flow_law_is_signed_and_independent_of_the_clock() {
        for flow in [-0.003, 0.003] {
            for dt in [1.0 / 192000.0, 0.02] {
                let (density, coefficient) = (1.225, 1.7);
                let mut p = port();
                p.resistance = 1500.0;
                p.set_flow_loss_coefficient(coefficient).unwrap();
                p.u = flow;
                let velocity = flow / p.area;
                let pressure =
                    p.resistance * flow + 0.5 * density * coefficient * velocity * velocity.abs();
                p.advance_pressure(dt, pressure, density);
                assert!((p.u / flow - 1.0).abs() < 2e-15);
                assert!((p.dissipated_j / (dt * pressure * flow) - 1.0).abs() < 3e-15);
            }
        }
    }

    #[test]
    fn pressure_work_closes_even_when_explicit_drag_would_be_unstable() {
        for resistance in [0.0, 1500.0, 1e9] {
            for coefficient in [0.0, 1.7, 1e8] {
                let mut p = port();
                p.resistance = resistance;
                p.set_flow_loss_coefficient(coefficient).unwrap();
                p.u = 0.003;
                let (density, dt) = (1.225, 0.002);
                let mass = density * p.length / p.area;
                for pressure in [100.0, -200.0, 0.0, 15.0, 0.0] {
                    let (before, loss) = (p.u, p.dissipated_j);
                    p.advance_pressure(dt, pressure, density);
                    let input = dt * pressure * 0.5 * (before + p.u);
                    let kinetic = 0.5 * mass * (p.u * p.u - before * before);
                    let removed = p.dissipated_j - loss;
                    let scale = (0.5 * mass * before * before).max(input.abs()).max(1e-20);
                    assert!(removed >= 0.0);
                    assert!(
                        (kinetic + removed - input).abs() < 2e-12 * scale,
                        "R={resistance} K={coefficient}: energy residual {}",
                        kinetic + removed - input
                    );
                    if pressure == 0.0 {
                        assert!(p.u.abs() <= before.abs() * (1.0 + 1e-14));
                    }
                }
            }
        }
    }

    #[test]
    fn nonlinear_flow_converges_at_second_order_to_the_analytic_transient() {
        let (density, coefficient, pressure, duration) = (1.225, 1.7, 100.0, 0.005);
        let template = port();
        let terminal = f64::sqrt(2.0 * pressure / (density * coefficient));
        let inverse_time = f64::sqrt(pressure * coefficient / (2.0 * density)) / template.length;
        let exact = template.area * terminal * f64::tanh(duration * inverse_time);
        let mut previous = None;
        for steps in [16, 32, 64] {
            let mut p = template.clone();
            p.set_flow_loss_coefficient(coefficient).unwrap();
            for _ in 0..steps {
                p.advance_pressure(duration / steps as f64, pressure, density);
            }
            let error = (p.u - exact).abs();
            assert!(error > 1e-12 * exact);
            if let Some(old) = previous {
                let ratio = old / error;
                assert!((3.8..4.2).contains(&ratio), "refinement ratio {ratio}");
            }
            previous = Some(error);
        }
        assert!(previous.unwrap() / exact < 2e-5);
    }

    #[test]
    fn causal_radiation_and_nonlinear_drag_share_one_pressure_work_balance() {
        for dt in [1.0 / 96000.0, 0.002] {
            let (density, c) = (1.225, 343.0);
            let mut p = port();
            p.enable_baffled_radiation(density, c, dt).unwrap();
            p.resistance = 1500.0;
            p.set_flow_loss_coefficient(1.7).unwrap();
            let inertia = density * p.length / p.area;
            let mut work = 0.0;
            let mut scale: f64 = 0.0;
            for i in 0..500 {
                let pressure = if i < 300 {
                    100.0 * (i as f64 * 0.137).cos()
                } else {
                    0.0
                };
                let old_flow = p.u;
                p.advance_pressure(dt, pressure, density);
                let input = dt * pressure * 0.5 * (old_flow + p.u);
                work += input;
                scale += input.abs();
                let accounted = 0.5 * inertia * p.u.powi(2)
                    + p.radiation_energy_j()
                    + p.dissipated_j
                    + p.radiated_j();
                assert!(
                    (accounted - work).abs() < 2e-12 * scale.max(1e-30),
                    "dt={dt}: pressure work {work}, stored+removed {accounted}"
                );
            }
            assert!(p.dissipated_j > 0.0 && p.radiated_j() > 0.0);
        }
    }

    #[test]
    fn radiation_attachment_is_atomic_and_replaces_only_the_exterior_mass() {
        let mut p = port();
        let original = p.length;
        for (rho, c, dt) in [
            (0.0, 343.0, 1e-5),
            (1.225, f64::NAN, 1e-5),
            (1.225, 343.0, 1e-20),
        ] {
            assert!(p.enable_baffled_radiation(rho, c, dt).is_err());
            assert_eq!(p.length, original);
            assert!(!p.has_baffled_radiation());
        }
        p.resistance = 1500.0;
        p.enable_baffled_radiation(1.225, 343.0, 1.0 / 96000.0)
            .unwrap();
        assert_eq!(p.exterior_end_correction, 0.0);
        assert!((p.length - (original - 0.85 * 0.025)).abs() < 1e-15);
        assert_eq!(p.interior_end_correction, 0.85 * 0.025);
        assert_eq!(
            p.resistance, 1500.0,
            "independent neck resistance is retained"
        );
        assert!(
            p.enable_baffled_radiation(1.225, 343.0, 1.0 / 96000.0)
                .is_err()
        );
    }
}
