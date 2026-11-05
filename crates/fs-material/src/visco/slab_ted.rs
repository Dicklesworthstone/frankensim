//! Thermoelastic bending of a homogeneous slab with insulated faces.
//!
//! The odd thermal modes have `tau_n = tau_0/(2n+1)^2` and weights
//! `96/[pi^4 (2n+1)^4]`. Their complete loss sum is evaluated in closed form,
//! without a thermal-mode truncation or a frequency-dependent loop. A positive
//! Taylor quotient avoids cancellation in the isothermal limit; scaled
//! exponentials avoid overflow in the adiabatic limit.
//!
//! Isotropic plane stress couples temperature to
//! `E alpha/(1-nu) (epsilon_x+epsilon_y)`, so only mean curvature drives loss.
//! The anisotropic law retains the full `Q alpha` coupling, including shear.
//! Both use the first-order, weak
//! thermoelastic coupling approximation: heat capacity is `rho c_p`, lateral
//! heat diffusion and surface heat exchange are omitted. Do not add this loss
//! to a model that already resolves the same thermomechanical heat transfer.
//!
//! References: Lifshitz & Roukes, Phys. Rev. B 61, 5600 (2000), eq. (29) and
//! appendix A, <https://doi.org/10.1103/PhysRevB.61.5600>; Norris & Photiadis,
//! Q. J. Mech. Appl. Math. 58, 145 (2005), thermal spectrum and modal
//! participation, <https://arxiv.org/abs/cond-mat/0405323>, eqs. (33)-(45).

use super::{ThermoelasticZener, ViscoError};
use fs_math::det;

/// Complete thermal diffusion loss sum per unit relaxation strength:
/// `sum_n 96/[pi^4(2n+1)^4] omega tau_n/[1+(omega tau_n)^2]`.
///
/// `first_relaxation_time` is `h^2 rho c_p/(pi^2 conductivity)`. Signed finite
/// frequency is supported (the loss part is odd); zero frequency has zero loss.
/// Invalid frequency or nonpositive/nonfinite relaxation time returns NaN.
#[must_use]
pub fn slab_relaxation_shape(omega: f64, first_relaxation_time: f64) -> f64 {
    if !(omega.is_finite() && first_relaxation_time.is_finite() && first_relaxation_time > 0.0) {
        return f64::NAN;
    }
    let pi2 = core::f64::consts::PI * core::f64::consts::PI;
    let wt = omega.abs() * first_relaxation_time;
    let shape = if wt <= 2.0 / pi2 {
        // xi^2 = pi^2 omega tau_0/2. The positive power series of
        // 6 (xi (cosh xi + cos xi) - sinh xi - sin xi) /
        // (xi^3 (cosh xi + cos xi)) has no subtractive cancellation.
        // For xi <= 1 the first omitted terms are below f64 roundoff.
        let xi2 = (0.5 * pi2) * wt;
        let xi4 = xi2 * xi2;
        let numerator = xi2
            * (0.2 + xi4 * (1.0 / 7560.0 + xi4 * (1.0 / 86_486_400.0 + xi4 / 3_705_077_376_000.0)));
        let denominator = 1.0
            + xi4
                * (1.0 / 24.0
                    + xi4
                        * (1.0 / 40_320.0
                            + xi4 * (1.0 / 479_001_600.0 + xi4 / 20_922_789_888_000.0)));
        numerator / denominator
    } else if wt.is_infinite() {
        // omega and tau are finite, but their product need not be. At this
        // scale the relative boundary-layer correction is below roundoff.
        ((12.0 / pi2) / omega.abs()) / first_relaxation_time
    } else {
        let xi = (core::f64::consts::PI / core::f64::consts::SQRT_2) * det::sqrt(wt);
        let ratio = if xi >= 40.0 {
            // The exponential correction is below f64 roundoff here; avoid
            // evaluating trig functions at arbitrarily large arguments.
            1.0
        } else {
            let decay = det::exp(-xi);
            let decay2 = decay * decay;
            (1.0 - decay2 + 2.0 * decay * det::sin(xi))
                / (1.0 + decay2 + 2.0 * decay * det::cos(xi))
        };
        ((12.0 / pi2) / wt) * (1.0 - ratio / xi)
    };
    shape.copysign(omega)
}


/// Homogeneous plane-stress slab with an explicitly supplied thermal expansion
/// tensor. All coefficients are at one frozen specimen state. This retains
/// weak thermoelastic coupling, insulated faces and through-thickness Fourier
/// conduction; lateral transport and a changing temperature field are absent.
///
/// For engineering strains [epsilon_xx, epsilon_yy, 2 epsilon_xy], beta=Q alpha.
/// The bending relaxation tensor is f f^T, with
/// f=sqrt(T h^3/(12 C_p)) beta. Thus every original curvature cross term survives
/// and dissipation vanishes whenever beta dot curvature vanishes.
/// See Norris & Photiadis, equations (3)-(4), (31)-(34), cited above.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SlabThermoelasticBending {
    relaxation_factor: [f64; 3],
    first_relaxation_time_s: f64,
}

impl SlabThermoelasticBending {
    /// Admit the actual isothermal plane-stress Q [Pa], row major, and expansion
    /// [alpha_xx, alpha_yy, 2 alpha_xy] [1/K] in the SAME coordinate frame.
    /// Heat capacity is per volume at constant stress [J/(m^3 K)]; conductivity
    /// is the through-thickness component [W/(m K)], not an in-plane average.
    ///
    /// # Errors
    /// Non-SPD Q, nonfinite expansion, nonpositive thermal/geometry data, or
    /// an unrepresentable relaxation tensor/time.
    pub fn new(
        plane_stress_pa: [f64; 9],
        engineering_expansion_per_k: [f64; 3],
        temperature_k: f64,
        volumetric_heat_capacity_j_m3_k: f64,
        thickness_conductivity_w_m_k: f64,
        thickness_m: f64,
    ) -> Result<Self, ViscoError> {
        let _ = crate::sheet::SheetMatrix::from_coefficients(plane_stress_pa)
            .map_err(|_| ViscoError::Parameters {
                what: "thermoelastic plane-stress stiffness must be symmetric positive definite",
            })?;
        if engineering_expansion_per_k.iter().any(|x| !x.is_finite())
            || [temperature_k, volumetric_heat_capacity_j_m3_k,
                thickness_conductivity_w_m_k, thickness_m]
                .iter().any(|x| !x.is_finite() || *x <= 0.0)
        {
            return Err(ViscoError::Parameters {
                what: "anisotropic slab needs finite expansion and positive thermal and thickness data",
            });
        }
        let scale = det::sqrt(temperature_k / volumetric_heat_capacity_j_m3_k)
            * thickness_m * det::sqrt(thickness_m / 12.0);
        let relaxation_factor = core::array::from_fn(|i| {
            scale * (0..3).map(|j| plane_stress_pa[3*i+j]
                * engineering_expansion_per_k[j]).sum::<f64>()
        });
        let length = thickness_m / core::f64::consts::PI;
        let first_relaxation_time_s = (volumetric_heat_capacity_j_m3_k
            / thickness_conductivity_w_m_k) * length * length;
        if !first_relaxation_time_s.is_finite() || first_relaxation_time_s <= 0.0
            || relaxation_factor.iter().any(|x| !x.is_finite() || !(x*x).is_finite())
        {
            return Err(ViscoError::Parameters {
                what: "anisotropic slab relaxation tensor or time is unrepresentable",
            });
        }
        Ok(Self { relaxation_factor, first_relaxation_time_s })
    }

    /// Original rank-one factor [sqrt(N m)] on engineering curvatures.
    /// Project this factor before squaring to retain weak and null directions.
    #[must_use]
    pub const fn relaxation_factor(&self) -> [f64; 3] { self.relaxation_factor }

    /// First insulated through-thickness thermal relaxation time [s].
    #[must_use]
    pub const fn first_relaxation_time_s(&self) -> f64 { self.first_relaxation_time_s }

    /// Full bending relaxation tensor [N m], row major. Cross terms can be
    /// signed; the quadratic form is nonnegative by construction.
    #[must_use]
    pub fn relaxation_matrix(&self) -> [f64; 9] {
        core::array::from_fn(|ij| self.relaxation_factor[ij/3]
            * self.relaxation_factor[ij%3])
    }

    /// Imaginary bending stiffness [N m] for exp(i omega t), retaining the
    /// complete insulated-slab spectrum. Signed frequency gives conjugacy.
    ///
    /// # Errors
    /// Nonfinite angular frequency.
    pub fn loss_matrix(&self, omega: f64) -> Result<[f64; 9], ViscoError> {
        if !omega.is_finite() {
            return Err(ViscoError::Parameters {
                what: "slab angular frequency must be finite",
            });
        }
        let shape = slab_relaxation_shape(omega, self.first_relaxation_time_s);
        Ok(self.relaxation_matrix().map(|value| value * shape))
    }
}

impl ThermoelasticZener {
    /// Weak-coupling loss factor of uniaxial bending across a slab of full
    /// thickness `thickness`, retaining the complete thermal diffusion spectrum.
    /// Plate consumers additionally need the plane-stress curvature coupling in
    /// [`Self::plate_bending_relaxation_matrix`]; this scalar is a beam law.
    #[must_use]
    pub fn slab_loss_factor(&self, omega: f64, thickness: f64) -> f64 {
        self.relaxation_strength() * slab_relaxation_shape(omega, self.relaxation_time(thickness))
    }

    /// Plane-stress bending relaxation tensor [N m] for engineering curvature
    /// `[kappa_xx, kappa_yy, 2 kappa_xy]`, row major. Its four nonzero entries are
    /// `T (E alpha/(1-nu))^2 h^3/(12 rho c_p)`.
    ///
    /// Multiplication by [`slab_relaxation_shape`] gives the imaginary bending
    /// stiffness in the `exp(i omega t)` convention. Project the full tensor on
    /// the original curvatures; multiplying all elastic bending energy by a
    /// scalar would spuriously damp trace-free bending and miss Poisson coupling.
    ///
    /// # Errors
    /// [`ViscoError::Parameters`] for invalid isotropic/thermal data, thickness,
    /// or unrepresentable relaxation stiffness or time.
    pub fn plate_bending_relaxation_matrix(
        &self,
        thickness: f64,
        poisson: f64,
    ) -> Result<[f64; 9], ViscoError> {
        if !self.alpha_t.is_finite()
            || !poisson.is_finite()
            || poisson <= -1.0
            || poisson >= 0.5
            || [
                self.e,
                self.t0,
                self.rho,
                self.cp,
                self.conductivity,
                thickness,
            ]
            .iter()
            .any(|x| !x.is_finite() || *x <= 0.0)
        {
            return Err(ViscoError::Parameters {
                what: "slab thermoelastic bending requires valid isotropic, thermal and thickness data",
            });
        }
        let poisson_scale = 1.0 - poisson;
        let coefficient = self.relaxation_strength()
            * (self.e / poisson_scale)
            * (thickness * thickness * thickness / (12.0 * poisson_scale));
        let tau = self.relaxation_time(thickness);
        if !(coefficient.is_finite() && coefficient >= 0.0 && tau.is_finite() && tau > 0.0) {
            return Err(ViscoError::Parameters {
                what: "slab thermoelastic bending stiffness or relaxation time is unrepresentable",
            });
        }
        Ok([
            coefficient,
            coefficient,
            0.0,
            coefficient,
            coefficient,
            0.0,
            0.0,
            0.0,
            0.0,
        ])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn thermal() -> ThermoelasticZener {
        ThermoelasticZener {
            e: 70e9,
            alpha_t: 20e-6,
            t0: 300.0,
            rho: 5000.0,
            cp: 600.0,
            conductivity: 100.0,
        }
    }

    #[test]
    fn complete_diffusion_response_matches_independent_thermal_modes() {
        // Direct eigenfunction expansion, not the production closed form.
        // 16,384 terms bound the relative omitted high-frequency loss below
        // 3e-10 over this range (and much more tightly near the peak).
        for wt in [
            1e-12, 1e-6, 0.01, 0.2, 0.203, 0.5, 1.0, 3.0, 10.0, 100.0, 1e4,
        ] {
            let expected: f64 = (0..16_384)
                .map(|n| {
                    let odd = (2 * n + 1) as f64;
                    let odd2 = odd * odd;
                    let x = wt / odd2;
                    let weight = 96.0 / (core::f64::consts::PI.powi(4) * odd2 * odd2);
                    weight * x / (1.0 + x * x)
                })
                .sum();
            let actual = slab_relaxation_shape(wt, 1.0);
            assert!(
                (actual / expected - 1.0).abs() < 3e-10,
                "wt={wt}: {actual} vs {expected}"
            );
            assert!(actual > 0.0 && actual <= 0.5);
        }
    }

    #[test]
    fn diffusion_limits_remain_positive_and_representable() {
        let pi2 = core::f64::consts::PI.powi(2);
        for wt in [1e-300, 1e-100, 1e-8] {
            assert!((slab_relaxation_shape(wt, 1.0) / (pi2 * wt / 10.0) - 1.0).abs() < 1e-14);
        }
        for (omega, tau) in [(1e100, 1.0), (1e250, 1.0), (1e200, 1e110)] {
            let expected = ((12.0 / pi2) / omega) / tau;
            let loss = slab_relaxation_shape(omega, tau);
            assert!(loss.is_finite() && loss > 0.0);
            assert!((loss / expected - 1.0).abs() < 1e-14);
            assert_eq!(slab_relaxation_shape(-omega, tau), -loss);
        }
        assert_eq!(slab_relaxation_shape(0.0, 0.3), 0.0);
        for (omega, tau) in [
            (f64::NAN, 1.0),
            (f64::INFINITY, 1.0),
            (1.0, 0.0),
            (1.0, -1.0),
        ] {
            assert!(slab_relaxation_shape(omega, tau).is_nan());
        }
    }

    #[test]
    fn plane_stress_loss_respects_curvature_and_recovers_uniaxial_bending() {
        let law = thermal();
        let h: f64 = 0.002;
        for nu in [-0.3, 0.0, 0.3, 0.49] {
            let d = law.plate_bending_relaxation_matrix(h, nu).unwrap();
            let energy = |k: [f64; 3]| -> f64 {
                (0..3)
                    .map(|i| (0..3).map(|j| k[i] * d[3 * i + j] * k[j]).sum::<f64>())
                    .sum()
            };
            assert_eq!(energy([1.0, -1.0, 0.0]), 0.0);
            assert_eq!(energy([0.0, 0.0, 2.0]), 0.0);
            // Lateral stress-free anticlastic bending has ky=-nu*kx.
            let beam_stiffness = law.e * h.powi(3) / 12.0;
            assert!(
                (energy([1.0, -nu, 0.0]) / (beam_stiffness * law.relaxation_strength()) - 1.0)
                    .abs()
                    < 1e-14
            );
            let cylindrical_stiffness = beam_stiffness / (1.0 - nu * nu);
            assert!(
                (energy([1.0, 0.0, 0.0])
                    / cylindrical_stiffness
                    / (law.relaxation_strength() * (1.0 + nu) / (1.0 - nu))
                    - 1.0)
                    .abs()
                    < 1e-14
            );
            assert!(energy([1.0, 1.0, 0.0]) > energy([1.0, 0.0, 0.0]));
        }
    }

    #[test]
    fn plate_relaxation_uses_explicit_thermal_data_and_refuses_invalid_inputs() {
        let base = thermal();
        let a = base.plate_bending_relaxation_matrix(0.002, 0.3).unwrap();
        let mut twin = base;
        twin.alpha_t *= 2.0;
        let b = twin.plate_bending_relaxation_matrix(0.002, 0.3).unwrap();
        assert_eq!(b[0], 4.0 * a[0]);
        twin.alpha_t = 0.0;
        assert_eq!(
            twin.plate_bending_relaxation_matrix(0.002, 0.3).unwrap(),
            [0.0; 9]
        );
        for (h, nu) in [
            (0.0, 0.3),
            (-0.1, 0.3),
            (0.002, -1.0),
            (0.002, 0.5),
            (0.002, f64::NAN),
        ] {
            assert!(base.plate_bending_relaxation_matrix(h, nu).is_err());
        }
        for bad in [f64::NAN, f64::INFINITY, -1.0, 0.0] {
            twin = base;
            twin.cp = bad;
            assert!(twin.plate_bending_relaxation_matrix(0.002, 0.3).is_err());
        }
    }
}

#[cfg(test)]
mod anisotropic_tests {
    use super::*;

    #[test]
    fn anisotropic_slab_retains_normal_shear_coupling_and_exact_nulls() {
        let q = [8e9, 0.0, 0.0, 0.0, 4e9, 0.0, 0.0, 0.0, 2e9];
        let law = SlabThermoelasticBending::new(q, [1e-5; 3], 300.0, 3e6, 10.0, 0.002)
            .unwrap();
        let beta = [80_000.0, 40_000.0, 20_000.0];
        let scale = 300.0 * 0.002_f64.powi(3) / (12.0 * 3e6);
        let matrix = law.relaxation_matrix();
        for i in 0..3 { for j in 0..3 {
            assert!((matrix[3*i+j] / (scale*beta[i]*beta[j]) - 1.0).abs() < 1e-14);
        } }
        let f = law.relaxation_factor();
        assert_eq!(f[0] - 2.0*f[1], 0.0);
        assert!(f[2] > 0.0, "engineering thermal shear must drive twist loss");
        let omega = law.first_relaxation_time_s().recip();
        // Independent Fourier eigenfunction expansion of the thermal field.
        let shape: f64 = (0..2048).map(|n| {
            let odd = (2*n+1) as f64;
            let x = 1.0 / (odd*odd);
            96.0 / (core::f64::consts::PI.powi(4)*odd.powi(4)) * x / (1.0+x*x)
        }).sum();
        let loss = law.loss_matrix(omega).unwrap();
        for i in 0..9 {
            assert!((loss[i] / (matrix[i]*shape) - 1.0).abs() < 1e-13);
        }
        assert_eq!(law.loss_matrix(-omega).unwrap(), loss.map(|x| -x));
        assert_eq!(law.loss_matrix(0.0).unwrap(), [0.0; 9]);
        let signed = SlabThermoelasticBending::new(q, [1e-5, 1e-5, -1e-5],
            300.0, 3e6, 10.0, 0.002).unwrap();
        assert!(signed.relaxation_matrix()[2] < 0.0);
    }

    #[test]
    fn general_slab_recovers_the_existing_isotropic_plate_tensor() {
        let original = ThermoelasticZener {
            e: 70e9, alpha_t: 20e-6, t0: 300.0, rho: 5000.0,
            cp: 600.0, conductivity: 100.0,
        };
        for nu in [-0.3, 0.0, 0.3, 0.49] {
            let q = original.e / (1.0-nu*nu);
            let law = SlabThermoelasticBending::new(
                [q, nu*q, 0.0, nu*q, q, 0.0, 0.0, 0.0, original.e/(2.0*(1.0+nu))],
                [original.alpha_t, original.alpha_t, 0.0], original.t0,
                original.rho*original.cp, original.conductivity, 0.002,
            ).unwrap();
            let expected = original.plate_bending_relaxation_matrix(0.002, nu).unwrap();
            for (actual, expected) in law.relaxation_matrix().iter().zip(expected) {
                assert!((*actual-expected).abs() < 1e-14 * expected.abs().max(1e-30));
            }
            assert!((law.first_relaxation_time_s()/original.relaxation_time(0.002)-1.0).abs() < 1e-14);
            let f = law.relaxation_factor();
            assert_eq!(f[0]-f[1], 0.0);
            assert_eq!(f[2], 0.0);
        }
    }

    #[test]
    fn general_slab_requires_actual_positive_thermal_data_and_spd_elasticity() {
        let q = [8.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 2.0];
        for bad in [0.0, -1.0, f64::INFINITY, f64::NAN] {
            assert!(SlabThermoelasticBending::new(q, [1e-5; 3], 300.0, 3e6, bad, 0.002).is_err());
        }
        let zero = SlabThermoelasticBending::new(q, [0.0; 3], 300.0, 3e6, 10.0, 0.002).unwrap();
        assert_eq!(zero.relaxation_matrix(), [0.0; 9]);
        assert!(zero.loss_matrix(f64::NAN).is_err());
        let mut invalid = q;
        invalid[1] = 1.0;
        assert!(SlabThermoelasticBending::new(invalid, [1e-5; 3], 300.0, 3e6, 10.0, 0.002).is_err());
        invalid = q;
        invalid[8] = -2.0;
        assert!(SlabThermoelasticBending::new(invalid, [1e-5; 3], 300.0, 3e6, 10.0, 0.002).is_err());
    }
}
