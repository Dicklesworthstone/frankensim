//! Quasi-steady constant-section fins with explicit lateral and tip convection.
//!
//! This is the one-dimensional boundary-value problem k A theta'' = h P theta,
//! theta(0) = theta_base, -k A theta'(L) = h_tip A theta(L). Coefficients
//! are uniform along the fin. The original material/film owner must admit any
//! spatial freezing approximation; this model supplies no fin heat storage,
//! radiation, contact resistance, array interference or distributed k(T) solve.
//!
//! Reference: Lienhard & Lienhard, A Heat Transfer Textbook, 6th ed., section 4.5,
//! equations (4.48)--(4.49), https://ahtt.mit.edu/wp-content/uploads/2024/04/AHTTv600.pdf.
//! The insulated-tip reduction and transverse-Biot assumption are also derived
//! in MIT Unified Engineering, Thermodynamics 18.2:
//! https://web.mit.edu/16.unified/www/FALL/thermodynamics/notes/node128.html.

/// Explicit physical or floating-point admission failure.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum FinError {
    Geometry,
    Coefficient,
    Position,
    Temperature,
    Arithmetic,
}
impl core::fmt::Display for FinError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "constant-section fin: {self:?}")
    }
}
impl std::error::Error for FinError {}

/// One axial length, conducting cross section and exposed lateral perimeter.
/// Perimeter excludes insulated lateral faces; no shape or fin count is inferred.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StraightFinGeometry {
    length_m: f64,
    cross_section_m2: f64,
    exposed_perimeter_m: f64,
}
impl StraightFinGeometry {
    pub fn try_new(
        length_m: f64,
        cross_section_m2: f64,
        exposed_perimeter_m: f64,
    ) -> Result<Self, FinError> {
        if [length_m, cross_section_m2, exposed_perimeter_m]
            .into_iter()
            .any(|v| !v.is_finite() || v <= 0.0)
        {
            return Err(FinError::Geometry);
        }
        let lateral_area = length_m * exposed_perimeter_m;
        if !lateral_area.is_finite() || lateral_area <= 0.0 {
            return Err(FinError::Arithmetic);
        }
        Ok(Self {
            length_m,
            cross_section_m2,
            exposed_perimeter_m,
        })
    }
    pub const fn length_m(self) -> f64 {
        self.length_m
    }
    pub const fn cross_section_m2(self) -> f64 {
        self.cross_section_m2
    }
    pub const fn exposed_perimeter_m(self) -> f64 {
        self.exposed_perimeter_m
    }
    pub fn lateral_area_m2(self) -> f64 {
        self.length_m * self.exposed_perimeter_m
    }
}

/// The tip is either insulated or exchanges over the physical cross section.
/// A convecting tip may have a different authored coefficient from the sides.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum FinTip {
    Adiabatic,
    Convection { coefficient_w_m2_k: f64 },
}

/// Exact constant-coefficient BVP response, with no axial mesh or allocation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StraightFinResponse {
    geometry: StraightFinGeometry,
    conductivity_w_m_k: f64,
    side_coefficient_w_m2_k: f64,
    tip_coefficient_w_m2_k: f64,
    m_length: f64,
    tip_to_axial_conductance: f64,
    conductance_w_k: f64,
    ideal_conductance_w_k: f64,
    transverse_biot: f64,
}
impl StraightFinResponse {
    /// Resolve once. The supplied h and k are model inputs, not material evidence.
    ///
    /// The transverse Biot number h A/(k P) is reported for the caller's
    /// explicit one-dimensional-model admission. No empirical cutoff is hidden.
    pub fn try_new(
        geometry: StraightFinGeometry,
        conductivity_w_m_k: f64,
        side_coefficient_w_m2_k: f64,
        tip: FinTip,
    ) -> Result<Self, FinError> {
        let h_tip = match tip {
            FinTip::Adiabatic => 0.0,
            FinTip::Convection { coefficient_w_m2_k } => coefficient_w_m2_k,
        };
        if !conductivity_w_m_k.is_finite()
            || conductivity_w_m_k <= 0.0
            || [side_coefficient_w_m2_k, h_tip]
                .into_iter()
                .any(|h| !h.is_finite() || h < 0.0)
        {
            return Err(FinError::Coefficient);
        }
        let axial = conductivity_w_m_k * (geometry.cross_section_m2 / geometry.length_m);
        let side = side_coefficient_w_m2_k * geometry.lateral_area_m2();
        let tip = h_tip * geometry.cross_section_m2;
        let ideal = side + tip;
        let transverse = (side_coefficient_w_m2_k / conductivity_w_m_k)
            * (geometry.cross_section_m2 / geometry.exposed_perimeter_m);
        if !axial.is_finite()
            || axial <= 0.0
            || !side.is_finite()
            || !tip.is_finite()
            || !ideal.is_finite()
            || !transverse.is_finite()
            || (side_coefficient_w_m2_k > 0.0 && side == 0.0)
            || (h_tip > 0.0 && tip == 0.0)
        {
            return Err(FinError::Arithmetic);
        }
        // sqrt(side)/sqrt(axial) avoids overflowing side/axial. tanh(x)/x
        // keeps the exact zero-side and short-fin limits without cancellation.
        let x = side.sqrt() / axial.sqrt();
        let b = tip / axial;
        if !x.is_finite() || !b.is_finite() {
            return Err(FinError::Arithmetic);
        }
        let t = tanh_over_x(x);
        let denominator = 1.0 + b * t;
        let conductance = (side * t) / denominator + tip / denominator;
        if !conductance.is_finite()
            || conductance < 0.0
            || (ideal > 0.0 && conductance == 0.0)
            || conductance > ideal
        {
            return Err(FinError::Arithmetic);
        }
        Ok(Self {
            geometry,
            conductivity_w_m_k,
            side_coefficient_w_m2_k,
            tip_coefficient_w_m2_k: h_tip,
            m_length: x,
            tip_to_axial_conductance: b,
            conductance_w_k: conductance,
            ideal_conductance_w_k: ideal,
            transverse_biot: transverse,
        })
    }
    pub const fn geometry(self) -> StraightFinGeometry {
        self.geometry
    }
    pub const fn conductivity_w_m_k(self) -> f64 {
        self.conductivity_w_m_k
    }
    pub const fn side_coefficient_w_m2_k(self) -> f64 {
        self.side_coefficient_w_m2_k
    }
    pub const fn tip_coefficient_w_m2_k(self) -> f64 {
        self.tip_coefficient_w_m2_k
    }
    pub const fn m_length(self) -> f64 {
        self.m_length
    }
    pub const fn conductance_w_k(self) -> f64 {
        self.conductance_w_k
    }
    /// Isothermal-at-base heat conductance, with the actual tip coefficient.
    pub const fn ideal_conductance_w_k(self) -> f64 {
        self.ideal_conductance_w_k
    }
    pub const fn transverse_biot(self) -> f64 {
        self.transverse_biot
    }
    /// A fully insulated fin has the continuous zero-exchange value one.
    pub fn efficiency(self) -> f64 {
        if self.ideal_conductance_w_k == 0.0 {
            1.0
        } else {
            self.conductance_w_k / self.ideal_conductance_w_k
        }
    }

    /// theta(x)/theta_base. Uses decaying exponentials rather than cosh(mL),
    /// so long fins cannot overflow. This also resolves the zero-side-convection
    /// limit with a finite tip load: a linear axial conduction profile.
    pub fn temperature_ratio_at(self, position_m: f64) -> Result<f64, FinError> {
        if !position_m.is_finite() || position_m < 0.0 || position_m > self.geometry.length_m {
            return Err(FinError::Position);
        }
        if position_m == 0.0 {
            return Ok(1.0);
        }
        let xi = position_m / self.geometry.length_m;
        let remaining = 1.0 - xi;
        let x = self.m_length;
        let y = x * remaining;
        let cosh_ratio = (-x * xi).exp() * (1.0 + (-2.0 * y).exp()) / (1.0 + (-2.0 * x).exp());
        let numerator = 1.0 + self.tip_to_axial_conductance * remaining * tanh_over_x(y);
        let denominator = 1.0 + self.tip_to_axial_conductance * tanh_over_x(x);
        let ratio = cosh_ratio * (numerator / denominator);
        if !ratio.is_finite() || !(0.0..=1.0).contains(&ratio) {
            return Err(FinError::Arithmetic);
        }
        Ok(ratio)
    }

    pub fn temperature_k_at(
        self,
        position_m: f64,
        base_temperature_k: f64,
        fluid_temperature_k: f64,
    ) -> Result<f64, FinError> {
        admit_temperatures(base_temperature_k, fluid_temperature_k)?;
        let ratio = self.temperature_ratio_at(position_m)?;
        let temperature = ratio.mul_add(
            base_temperature_k - fluid_temperature_k,
            fluid_temperature_k,
        );
        if !temperature.is_finite() || temperature <= 0.0 {
            return Err(FinError::Arithmetic);
        }
        Ok(temperature)
    }

    /// Positive heat leaves the base and enters the fluid; cooling/heating use
    /// the same conductance without clipping the signed temperature difference.
    pub fn base_heat_w(
        self,
        base_temperature_k: f64,
        fluid_temperature_k: f64,
    ) -> Result<f64, FinError> {
        admit_temperatures(base_temperature_k, fluid_temperature_k)?;
        let heat = self.conductance_w_k * (base_temperature_k - fluid_temperature_k);
        if !heat.is_finite() {
            return Err(FinError::Arithmetic);
        }
        Ok(heat)
    }
}

fn admit_temperatures(base: f64, fluid: f64) -> Result<(), FinError> {
    if !base.is_finite() || base <= 0.0 || !fluid.is_finite() || fluid <= 0.0 {
        Err(FinError::Temperature)
    } else {
        Ok(())
    }
}

fn tanh_over_x(x: f64) -> f64 {
    if x < 1e-3 {
        let y = x * x;
        1.0 + y * (-1.0 / 3.0 + y * (2.0 / 15.0 - y * 17.0 / 315.0))
    } else {
        x.tanh() / x
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn near(a: f64, b: f64, tolerance: f64) {
        assert!((a - b).abs() <= tolerance, "{a} != {b}");
    }
    fn geometry() -> StraightFinGeometry {
        StraightFinGeometry::try_new(0.04, 0.0002, 0.042).unwrap()
    }

    #[test]
    fn published_adiabatic_and_convective_tip_relations_and_signed_heat() {
        let g = geometry();
        for tip in [
            FinTip::Adiabatic,
            FinTip::Convection {
                coefficient_w_m2_k: 70.0,
            },
        ] {
            let k = 150.0;
            let h = 25.0;
            let solved = StraightFinResponse::try_new(g, k, h, tip).unwrap();
            let m = (h * g.exposed_perimeter_m() / (k * g.cross_section_m2())).sqrt();
            let ht = solved.tip_coefficient_w_m2_k();
            let r = ht / (k * m);
            let x = m * g.length_m();
            let exact = k * g.cross_section_m2() * m * (x.sinh() + r * x.cosh())
                / (x.cosh() + r * x.sinh());
            near(solved.conductance_w_k(), exact, 1e-15);
            for position in [0.0, 0.01, 0.03, 0.04] {
                let y = m * (g.length_m() - position);
                near(
                    solved.temperature_ratio_at(position).unwrap(),
                    (y.cosh() + r * y.sinh()) / (x.cosh() + r * x.sinh()),
                    1e-15,
                );
            }
            near(
                solved.base_heat_w(350.0, 300.0).unwrap(),
                50.0 * exact,
                5e-14,
            );
            assert_eq!(
                solved.base_heat_w(300.0, 350.0).unwrap(),
                -solved.base_heat_w(350.0, 300.0).unwrap()
            );
        }
    }

    // Independent finite-volume BVP: original dimensional axial conductance,
    // distributed side loss, and a half control volume with a Robin tip.
    fn finite_volume(g: StraightFinGeometry, k: f64, h: f64, ht: f64, n: usize) -> (f64, f64) {
        let dx = g.length_m() / n as f64;
        let c = k * g.cross_section_m2() / dx;
        let loss = h * g.exposed_perimeter_m() * dx;
        let tip = ht * g.cross_section_m2();
        let mut upper = vec![0.0; n];
        let mut rhs = vec![0.0; n];
        for i in 0..n {
            let diagonal = if i + 1 == n {
                c + 0.5 * loss + tip
            } else {
                2.0 * c + loss
            };
            let denominator = diagonal - if i == 0 { 0.0 } else { -c * upper[i - 1] };
            upper[i] = if i + 1 == n { 0.0 } else { -c / denominator };
            rhs[i] = (if i == 0 { c } else { c * rhs[i - 1] }) / denominator;
        }
        for i in (0..n - 1).rev() {
            rhs[i] -= upper[i] * rhs[i + 1];
        }
        let base = c * (1.0 - rhs[0]) + 0.5 * loss;
        let fluid =
            loss * (0.5 + rhs[..n - 1].iter().sum::<f64>() + 0.5 * rhs[n - 1]) + tip * rhs[n - 1];
        near(base, fluid, 5e-11);
        (base, rhs[n - 1])
    }

    #[test]
    fn dimensional_bvp_refinement_agrees_with_conductance_and_tip_temperature() {
        let g = geometry();
        let exact = StraightFinResponse::try_new(
            g,
            7.0,
            80.0,
            FinTip::Convection {
                coefficient_w_m2_k: 120.0,
            },
        )
        .unwrap();
        let mut errors = Vec::new();
        let mut tip_errors = Vec::new();
        for n in [64, 128, 256] {
            let (conductance, tip) = finite_volume(g, 7.0, 80.0, 120.0, n);
            errors.push((conductance - exact.conductance_w_k()).abs());
            tip_errors.push((tip - exact.temperature_ratio_at(g.length_m()).unwrap()).abs());
        }
        assert!(errors[0] / errors[1] > 3.8);
        assert!(errors[1] / errors[2] > 3.8);
        assert!(errors[2] < 1e-6, "conductance errors: {errors:?}");
        assert!(tip_errors[0] / tip_errors[1] > 3.8);
        assert!(tip_errors[1] / tip_errors[2] > 3.8);
        assert!(tip_errors[2] < 6e-6);
    }

    #[test]
    fn zero_side_tip_resistance_short_and_long_limits_remain_finite() {
        let g = geometry();
        let isolated = StraightFinResponse::try_new(g, 100.0, 0.0, FinTip::Adiabatic).unwrap();
        assert_eq!(isolated.conductance_w_k(), 0.0);
        assert_eq!(isolated.efficiency(), 1.0);
        assert_eq!(isolated.temperature_ratio_at(g.length_m()).unwrap(), 1.0);
        let end = StraightFinResponse::try_new(
            g,
            100.0,
            0.0,
            FinTip::Convection {
                coefficient_w_m2_k: 40.0,
            },
        )
        .unwrap();
        let axial = 100.0 * g.cross_section_m2() / g.length_m();
        let tip = 40.0 * g.cross_section_m2();
        near(
            end.conductance_w_k(),
            1.0 / (1.0 / axial + 1.0 / tip),
            1e-18,
        );
        near(
            end.temperature_ratio_at(0.02).unwrap(),
            1.0 - 0.5 * end.conductance_w_k() / axial,
            2e-16,
        );
        let short = StraightFinResponse::try_new(g, 1e20, 10.0, FinTip::Adiabatic).unwrap();
        near(short.efficiency(), 1.0, 1e-15);
        let long = StraightFinResponse::try_new(
            StraightFinGeometry::try_new(1000.0, 0.01, 1.0).unwrap(),
            1.0,
            1.0,
            FinTip::Adiabatic,
        )
        .unwrap();
        assert_eq!(long.m_length(), 10000.0);
        near(long.conductance_w_k(), 0.1, 1e-16);
        assert_eq!(long.temperature_ratio_at(1000.0).unwrap(), 0.0);
    }

    #[test]
    fn geometry_rescaling_and_original_input_refusals_are_explicit() {
        let g = geometry();
        let original = StraightFinResponse::try_new(
            g,
            20.0,
            30.0,
            FinTip::Convection {
                coefficient_w_m2_k: 60.0,
            },
        )
        .unwrap();
        let scale = 1000.0;
        let resized = StraightFinResponse::try_new(
            StraightFinGeometry::try_new(
                scale * g.length_m(),
                scale * scale * g.cross_section_m2(),
                scale * g.exposed_perimeter_m(),
            )
            .unwrap(),
            20.0,
            30.0 / scale,
            FinTip::Convection {
                coefficient_w_m2_k: 60.0 / scale,
            },
        )
        .unwrap();
        near(
            resized.conductance_w_k(),
            scale * original.conductance_w_k(),
            2e-14,
        );
        near(resized.efficiency(), original.efficiency(), 3e-16);
        near(
            resized.temperature_ratio_at(scale * 0.03).unwrap(),
            original.temperature_ratio_at(0.03).unwrap(),
            3e-16,
        );
        assert!(StraightFinGeometry::try_new(0.0, 1.0, 1.0).is_err());
        assert!(StraightFinResponse::try_new(g, 0.0, 1.0, FinTip::Adiabatic).is_err());
        assert!(StraightFinResponse::try_new(g, 1.0, -1.0, FinTip::Adiabatic).is_err());
        assert!(
            original
                .temperature_ratio_at(g.length_m().next_up())
                .is_err()
        );
        assert!(original.base_heat_w(0.0, 300.0).is_err());
    }
}
