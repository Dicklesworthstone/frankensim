//! Air loss laws read from the gas record, with no fitted constants (.372).
//!
//! - [`oscillating_cylinder`]: the drag on a snare-wire coil, an infinite circular
//!   cylinder oscillating in an unbounded viscous gas. This is the exact linear
//!   solution of Stokes (1851) in Stuart's (1963) closed form, as reproduced by
//!   Dolfo, Vigué & Lhuillier, "Stokes force on a cylinder in the presence of
//!   fluid confinement", arXiv:2011.12000, eq. (60):
//!   k - i k' = 1 + 4 K1(alpha) / (alpha K0(alpha)), with alpha = (1 + i) a / delta,
//!   delta = sqrt(2 nu / omega) and time dependence exp(i omega t). The force per
//!   unit length is -(pi a^2 rho k) dU/dt - (omega pi a^2 rho k') U.
//!   Validity: linear unsteady Stokes flow, so Keulegan-Carpenter A/a << 1. The
//!   caller reports the KC it reaches; this function has no amplitude input.
//! - [`plate_viscous_layer`]: the viscous (Stokes) layer on the faces of a plate
//!   vibrating in a subsonic bending wave. Above a free bending wave
//!   w ~ exp(i (k_b x - omega t)) the potential near field phi ~ exp(i k_b x - gamma z),
//!   gamma = sqrt(k_b^2 - k_0^2), has |v_t| / |v_n| = k_b / gamma =
//!   1 / sqrt(1 - f / f_c) at the face; the face itself is tangentially at rest to
//!   O(k_b h / 2). The no-slip layer over an oscillating tangential slip V
//!   dissipates mu V^2 / (2 delta) per unit area (Landau & Lifshitz, Fluid
//!   Mechanics, §24). Dividing by omega times the areal kinetic energy gives
//!   eta = faces * (rho_0 delta / (2 rho_s)) / (1 - f / f_c), which for two faces is
//!   rho_0 delta / rho_s / (1 - f / f_c). Validity: f < f_c and k_b delta << 1.
//!
//! The complex modified Bessel functions come from the integral
//! K_nu(z) = integral_0^inf exp(-z cosh t) cosh(nu t) dt (Re z > 0), scaled by
//! exp(z) so no factor underflows, on an even-integrand trapezoid rule. Its step
//! shrinks as 1/sqrt(|z|) and its cutoff keeps exp(-Re z (cosh t - 1)) above
//! exp(-46), so every call costs about 100 to 250 evaluations. The tests pin it
//! against the independent ascending series (A&S 9.6.10, 9.6.11, 9.6.13), and the
//! drag against the paper's own low- and high-frequency limits.

use crate::GasState;
use std::f64::consts::PI;
use std::ops::{Add, Div, Mul, Sub};

/// A complex number; this module needs only these operations.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Complex {
    pub(crate) re: f64,
    pub(crate) im: f64,
}

impl Complex {
    pub(crate) const fn new(re: f64, im: f64) -> Complex {
        Complex { re, im }
    }
    pub(crate) fn abs(self) -> f64 {
        self.re.hypot(self.im)
    }
    fn scale(self, s: f64) -> Complex {
        Complex::new(self.re * s, self.im * s)
    }
    fn exp(self) -> Complex {
        let m = self.re.exp();
        Complex::new(m * self.im.cos(), m * self.im.sin())
    }
}

impl Add for Complex {
    type Output = Complex;
    fn add(self, o: Complex) -> Complex {
        Complex::new(self.re + o.re, self.im + o.im)
    }
}

impl Sub for Complex {
    type Output = Complex;
    fn sub(self, o: Complex) -> Complex {
        Complex::new(self.re - o.re, self.im - o.im)
    }
}

impl Mul for Complex {
    type Output = Complex;
    fn mul(self, o: Complex) -> Complex {
        Complex::new(
            self.re * o.re - self.im * o.im,
            self.re * o.im + self.im * o.re,
        )
    }
}

impl Div for Complex {
    type Output = Complex;
    fn div(self, o: Complex) -> Complex {
        let d = o.re * o.re + o.im * o.im;
        Complex::new(
            (self.re * o.re + self.im * o.im) / d,
            (self.im * o.re - self.re * o.im) / d,
        )
    }
}

/// exp(z) K0(z) and exp(z) K1(z) for Re z > 0.
pub(crate) fn scaled_bessel_k01(z: Complex) -> Result<(Complex, Complex), String> {
    if !(z.re.is_finite() && z.im.is_finite()) || z.re <= 0.0 {
        return Err(format!(
            "modified Bessel K needs a finite z with Re z > 0, got {z:?}"
        ));
    }
    let t_max = (1.0 + 46.0 / z.re).acosh();
    let h = 0.1 / z.abs().max(1.0).sqrt();
    let steps = (t_max / h).ceil() as usize;
    let mut k0 = Complex::new(0.5, 0.0);
    let mut k1 = Complex::new(0.5, 0.0);
    for j in 1..=steps {
        let t = j as f64 * h;
        let e = (z.scale(-(t.cosh() - 1.0))).exp();
        k0 = k0 + e;
        k1 = k1 + e.scale(t.cosh());
    }
    Ok((k0.scale(h), k1.scale(h)))
}

/// Added-mass and damping coefficients of an oscillating cylinder, and the
/// per-unit-length quantities they give in the declared gas.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CylinderDrag {
    /// Cylinder radius a [m].
    pub radius_m: f64,
    /// Viscous penetration depth delta = sqrt(2 nu / omega) [m].
    pub penetration_depth_m: f64,
    /// Added-mass coefficient k (1 for inviscid potential flow).
    pub added_mass_coefficient: f64,
    /// Damping coefficient k'.
    pub damping_coefficient: f64,
    /// Added mass pi a^2 rho k [kg/m].
    pub added_mass_kg_per_m: f64,
    /// Drag resistance omega pi a^2 rho k' [N s/m^2]: force per unit length per unit speed.
    pub resistance_n_s_per_m2: f64,
}

impl CylinderDrag {
    /// a / delta.
    #[must_use]
    pub fn radius_over_penetration(&self) -> f64 {
        self.radius_m / self.penetration_depth_m
    }
}

fn gas_density_viscosity(gas: &GasState) -> Result<(f64, f64), String> {
    let (rho, mu) = (gas.density, gas.dynamic_viscosity);
    if !(rho.is_finite() && rho > 0.0 && mu.is_finite() && mu > 0.0) {
        return Err(format!(
            "air loss needs the gas record's density and dynamic viscosity, got rho {rho} kg/m^3, mu {mu} Pa s"
        ));
    }
    Ok((rho, mu))
}

/// Stuart's exact linear drag on a cylinder of radius `radius_m` oscillating at
/// `omega_rad_s` in `gas`.
pub fn oscillating_cylinder(
    radius_m: f64,
    omega_rad_s: f64,
    gas: &GasState,
) -> Result<CylinderDrag, String> {
    let (rho, mu) = gas_density_viscosity(gas)?;
    if !(radius_m.is_finite() && radius_m > 0.0 && omega_rad_s.is_finite() && omega_rad_s > 0.0) {
        return Err(format!(
            "oscillating-cylinder drag needs a positive radius and angular frequency, got a {radius_m} m, omega {omega_rad_s} rad/s"
        ));
    }
    let delta = (2.0 * mu / (rho * omega_rad_s)).sqrt();
    let r = radius_m / delta;
    let alpha = Complex::new(r, r);
    let (k0, k1) = scaled_bessel_k01(alpha)?;
    let h = Complex::new(1.0, 0.0) + Complex::new(4.0, 0.0) * k1 / (alpha * k0);
    let (k, k_prime) = (h.re, -h.im);
    let displaced = PI * radius_m * radius_m * rho;
    Ok(CylinderDrag {
        radius_m,
        penetration_depth_m: delta,
        added_mass_coefficient: k,
        damping_coefficient: k_prime,
        added_mass_kg_per_m: displaced * k,
        resistance_n_s_per_m2: omega_rad_s * displaced * k_prime,
    })
}

/// The viscous-layer loss of a plate face pair in a subsonic bending wave.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PlateViscousLayer {
    /// Viscous penetration depth delta = sqrt(2 nu / omega) [m].
    pub penetration_depth_m: f64,
    /// |v_t| / |v_n| at the face, 1 / sqrt(1 - f / f_c).
    pub slip_over_normal: f64,
    /// Modal loss factor eta from the layer on `faces` faces.
    pub loss_factor: f64,
}

/// Viscous-layer loss factor for a plate of areal density `areal_density_kg_m2`
/// at `omega_rad_s`, `frequency_over_coincidence` = f / f_c in this gas, with
/// `faces` (1 or 2) wetted faces.
pub fn plate_viscous_layer(
    areal_density_kg_m2: f64,
    omega_rad_s: f64,
    frequency_over_coincidence: f64,
    faces: u8,
    gas: &GasState,
) -> Result<PlateViscousLayer, String> {
    let (rho, mu) = gas_density_viscosity(gas)?;
    if !(areal_density_kg_m2.is_finite()
        && areal_density_kg_m2 > 0.0
        && omega_rad_s.is_finite()
        && omega_rad_s > 0.0)
    {
        return Err(format!(
            "plate viscous layer needs a positive areal density and angular frequency, got {areal_density_kg_m2} kg/m^2, {omega_rad_s} rad/s"
        ));
    }
    if !(frequency_over_coincidence.is_finite()
        && frequency_over_coincidence >= 0.0
        && frequency_over_coincidence < 1.0)
    {
        return Err(format!(
            "plate viscous layer holds only below coincidence (subsonic bending wave), got f/f_c {frequency_over_coincidence}"
        ));
    }
    if !(faces == 1 || faces == 2) {
        return Err(format!("a plate has one or two wetted faces, got {faces}"));
    }
    let delta = (2.0 * mu / (rho * omega_rad_s)).sqrt();
    let slip2 = 1.0 / (1.0 - frequency_over_coincidence);
    Ok(PlateViscousLayer {
        penetration_depth_m: delta,
        slip_over_normal: slip2.sqrt(),
        loss_factor: f64::from(faces) * 0.5 * rho * delta / areal_density_kg_m2 * slip2,
    })
}
