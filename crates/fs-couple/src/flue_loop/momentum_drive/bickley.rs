//! Advected momentum of the same Bickley profile that supplies the cut flow.
//!
//! For `u(y) = U sech²((y-eta)/b)`, span H and `s = tanh((eta-offset)/b)`,
//! the selected volume flux is `Q = U b H (1+s)`. Its streamwise momentum
//! flux is `J = rho U² b H (2/3+s-s³/3)`, from the integral of `sech⁴`.
//! With `f=Q/Qmax`, `J = rho U Qmax (2f²-4f³/3)` and
//! `dJ/dQ = 4 rho U f(1-f)`, bounded by `rho U`.
//!
//! The series source assumes that the selected streamwise momentum is
//! deposited in the bore. It does not resolve the labium reaction, entrainment,
//! or downstream jet, and is not an experimentally calibrated sound model.

/// Profile-integrated source with the voice's original gas density and bore.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BickleyMomentumDrive {
    density_kg_m3: f64,
    bore_area_m2: f64,
}

impl BickleyMomentumDrive {
    pub fn new(density_kg_m3: f64, bore_area_m2: f64) -> Result<Self, &'static str> {
        if !density_kg_m3.is_finite()
            || density_kg_m3 <= 0.0
            || !bore_area_m2.is_finite()
            || bore_area_m2 <= 0.0
        {
            return Err("Bickley momentum requires positive finite density and bore area");
        }
        Ok(Self {
            density_kg_m3,
            bore_area_m2,
        })
    }

    #[must_use]
    pub const fn bore_area_m2(self) -> f64 {
        self.bore_area_m2
    }

    /// Upper derivative of pressure with respect to cut flow [Pa s/m³].
    /// The source solver separately admits this bound as finite/nonnegative.
    #[must_use]
    pub fn pressure_per_flow_bound(self, centre_speed_m_s: f64) -> f64 {
        self.density_kg_m3 * centre_speed_m_s / self.bore_area_m2
    }

    /// Actual selected streamwise momentum flux [N], including its steady part.
    /// `maximum_flow_m3_s` is the SAME station's `2 U b H`, not exit bulk flow.
    pub fn momentum_flux_n(
        self,
        centre_speed_m_s: f64,
        maximum_flow_m3_s: f64,
        cut_flow_m3_s: f64,
    ) -> Result<f64, &'static str> {
        admit(centre_speed_m_s, maximum_flow_m3_s, &[cut_flow_m3_s])?;
        if maximum_flow_m3_s == 0.0 {
            return Ok(0.0);
        }
        let fraction = cut_flow_m3_s / maximum_flow_m3_s;
        finite(
            self.density_kg_m3
                * centre_speed_m_s
                * cut_flow_m3_s
                * fraction
                * (2.0 - (4.0 / 3.0) * fraction),
        )
    }

    /// Pressure deviation from the quiescent cut of this SAME profile [Pa].
    /// The factored cubic difference retains small changes about a large jet.
    pub fn pressure_pa(
        self,
        centre_speed_m_s: f64,
        maximum_flow_m3_s: f64,
        cut_flow_m3_s: f64,
        rest_flow_m3_s: f64,
    ) -> Result<f64, &'static str> {
        admit(
            centre_speed_m_s,
            maximum_flow_m3_s,
            &[cut_flow_m3_s, rest_flow_m3_s],
        )?;
        if maximum_flow_m3_s == 0.0 || cut_flow_m3_s == rest_flow_m3_s {
            return Ok(0.0);
        }
        let a = cut_flow_m3_s / maximum_flow_m3_s;
        let b = rest_flow_m3_s / maximum_flow_m3_s;
        let ca = (maximum_flow_m3_s - cut_flow_m3_s) / maximum_flow_m3_s;
        let cb = (maximum_flow_m3_s - rest_flow_m3_s) / maximum_flow_m3_s;
        // 2(a+b)-4(a²+ab+b²)/3, expressed as nonnegative products so
        // a cut near full deflection does not subtract nearly equal slopes.
        let secant = (4.0 / 3.0) * (a * ca + b * cb) + (2.0 / 3.0) * (a * cb + b * ca);
        finite(
            self.pressure_per_flow_bound(centre_speed_m_s)
                * (cut_flow_m3_s - rest_flow_m3_s)
                * secant,
        )
    }
}

fn admit(speed: f64, maximum: f64, flows: &[f64]) -> Result<(), &'static str> {
    if !speed.is_finite()
        || speed < 0.0
        || !maximum.is_finite()
        || maximum < 0.0
        || (speed == 0.0) != (maximum == 0.0)
        || flows
            .iter()
            .any(|q| !q.is_finite() || *q < 0.0 || *q > maximum)
    {
        return Err(
            "Bickley momentum needs one finite station speed and cut flows within its total flow",
        );
    }
    Ok(())
}
fn finite(value: f64) -> Result<f64, &'static str> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err("Bickley momentum source is not finitely representable")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn near(a: f64, b: f64) {
        assert!(
            (a - b).abs() <= 2e-13 * b.abs().max(1e-12),
            "{a:e} != {b:e}"
        );
    }

    #[test]
    fn g1_cut_flux_matches_direct_velocity_squared_quadrature() {
        let (rho, u, width, span) = (1.2, 12.0, 0.0007, 0.01);
        let drive = BickleyMomentumDrive::new(rho, 3e-4).unwrap();
        let maximum = 2.0 * u * width * span;
        near(
            drive.momentum_flux_n(u, maximum, maximum).unwrap(),
            (4.0 / 3.0) * rho * u * u * width * span,
        );
        for cut in [-3.0_f64, -0.5, 0.0, 0.75, 4.0] {
            // Integrate rho*u(y)^2 independently through the selected cut.
            let panels = 100_000;
            let h = (cut + 20.0) / panels as f64;
            let mut sum = 0.0;
            for i in 0..=panels {
                let z = -20.0 + i as f64 * h;
                let velocity = u / z.cosh().powi(2);
                let weight = if i == 0 || i == panels {
                    1.0
                } else if i % 2 == 0 {
                    2.0
                } else {
                    4.0
                };
                sum += weight * rho * velocity * velocity * width * span;
            }
            let q = 0.5 * maximum * (1.0 + cut.tanh());
            near(drive.momentum_flux_n(u, maximum, q).unwrap(), sum * h / 3.0);
        }
    }

    #[test]
    fn g0_g3_cut_symmetry_rest_and_derivative_bound() {
        let drive = BickleyMomentumDrive::new(1.2, 3e-4).unwrap();
        let (u, maximum) = (10.0, 2e-4);
        let full = drive.momentum_flux_n(u, maximum, maximum).unwrap();
        for f in [0.0, 0.01, 0.25, 0.5, 0.75, 0.99, 1.0] {
            let q = f * maximum;
            near(
                drive.momentum_flux_n(u, maximum, q).unwrap()
                    + drive.momentum_flux_n(u, maximum, maximum - q).unwrap(),
                full,
            );
            assert_eq!(drive.pressure_pa(u, maximum, q, q).unwrap(), 0.0);
            let p = drive.pressure_pa(u, maximum, q, 0.5 * maximum).unwrap();
            near(
                p,
                -drive
                    .pressure_pa(u, maximum, maximum - q, 0.5 * maximum)
                    .unwrap(),
            );
            let dq = maximum * 1e-5;
            if f > 0.0 && f < 1.0 {
                let slope = drive.pressure_pa(u, maximum, q + dq, q - dq).unwrap() / (2.0 * dq);
                let expected = drive.pressure_per_flow_bound(u) * 4.0 * f * (1.0 - f);
                assert!((slope - expected).abs() < 2e-10 * drive.pressure_per_flow_bound(u));
                assert!(slope >= 0.0 && slope <= drive.pressure_per_flow_bound(u));
            }
        }
    }

    #[test]
    fn g3_spreading_preserves_total_momentum_and_small_cut_changes() {
        let drive = BickleyMomentumDrive::new(1.0, 1.0).unwrap();
        let initial = drive.momentum_flux_n(12.0, 2e-4, 2e-4).unwrap();
        for beta in [1.0_f64, 4.0, 100.0] {
            let u = 12.0 / beta.sqrt();
            let q = 2e-4 * beta.sqrt();
            near(drive.momentum_flux_n(u, q, q).unwrap(), initial);
        }
        let q = 0.75_f64;
        let dq = q.next_up() - q;
        near(drive.pressure_pa(1.0, 1.0, q + dq, q).unwrap() / dq, 0.75);
        assert_eq!(drive.pressure_pa(0.0, 0.0, 0.0, 0.0).unwrap(), 0.0);
        for q in [-0.1, 1.1, f64::NAN] {
            assert!(drive.pressure_pa(1.0, 1.0, q, 0.5).is_err());
        }
        assert!(drive.pressure_pa(0.0, 1.0, 0.5, 0.5).is_err());
    }
}
