//! Allocation-free Stokes/Mueller polarization algebra for spectral light transport.
//!
//! Stokes components use the convention S=(I,Q,U,V), with Q positive for
//! horizontal linear polarization, U positive at +45 degrees and V positive
//! for the chosen right-handed circular basis. All operations assume the same
//! propagation direction and transverse reference frame. Rotate frames before
//! composing optical elements; this module does not silently rotate them.
//! Incoherent light is represented by Stokes vectors, not Jones amplitudes.

/// Polarized Fresnel, spectral materials and native smooth-sampler integration.
#[cfg(feature = "tracer")]
pub mod dielectric;

#[cfg(feature = "photometry")]
mod source_coherency;

/// A physical (possibly partially polarized) incoherent light state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Stokes {
    pub i: f64,
    pub q: f64,
    pub u: f64,
    pub v: f64,
}

impl Stokes {
    /// Reject nonfinite, negative-intensity or super-physical polarization.
    pub fn try_new(i: f64, q: f64, u: f64, v: f64) -> Option<Self> {
        if ![i, q, u, v].iter().all(|x| x.is_finite()) || i < 0.0 {
            return None;
        }
        // Division by I avoids squaring large radiances and overflow.
        if i == 0.0 {
            return (q == 0.0 && u == 0.0 && v == 0.0).then_some(Self { i, q, u, v });
        }
        let norm = (q / i).hypot(u / i).hypot(v / i);
        (norm <= 1.0 + 16.0 * f64::EPSILON).then_some(Self { i, q, u, v })
    }

    pub const fn unpolarized(i: f64) -> Self {
        Self { i, q: 0.0, u: 0.0, v: 0.0 }
    }

    /// Rotate the reference axes counterclockwise by angle (radians).
    /// The Stokes Q/U coordinates rotate by minus twice the angle.
    pub fn rotate_axes(self, angle: f64) -> Self {
        let (s, c) = (2.0 * angle).sin_cos();
        Self { i: self.i, q: c * self.q + s * self.u,
            u: -s * self.q + c * self.u, v: self.v }
    }

    /// Degree of polarization; zero for a dark state.
    pub fn degree(self) -> f64 {
        if self.i == 0.0 { 0.0 } else {
            (self.q / self.i).hypot(self.u / self.i).hypot(self.v / self.i)
        }
    }
}

/// Mueller matrix acting on a Stokes column vector. The matrix is row-major.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Mueller(pub [[f64; 4]; 4]);

// Rounding allowance of the weak channel that `apply_diagonal_jones` rebuilds from the cross entries.
//
// In exact arithmetic a nondepolarizing diagonal Jones matrix has |cross|^2 = R_s R_p, so the weak power rebuilt as
// |cross|^2 / large cannot exceed `large`. In binary64 it can, by rounding alone: a total internal reflection at cosine
// 0.6 from index 1.5 to 1.0 published an intensity of 1 + 2^-52 for unit input.
//
// u = 2^-53 and gamma_k = k u / (1 - k u). Higham, Accuracy and Stability of Numerical Algorithms, 2nd ed., Lemma 3.1:
// prod (1 + d_i)^(+-1) = 1 + theta_k with |theta_k| <= gamma_k. For the total internal reflection matrix of
// polarization/dielectric/interface.rs (large = 1 exactly), with (a, b) the admittances handed to `tir_amplitude`:
//   1. d = fl(a^2 + b^2) has relative error gamma_2, both terms being positive. The difference a^2 - b^2 cancels, so its
//      error is gamma_2 (a^2 + b^2): absolute, and relative to |r| = 1, not to the component. Hence sr = fl(n / d) is
//      within gamma_2 (1 + gamma_3) + |x| gamma_3 of x = (a^2 - b^2) / (a^2 + b^2), and si = fl(fl(-2 a b) / d) has
//      relative error gamma_4. With x^2 + y^2 = 1: |sr^2 + si^2 - 1| <= 2 (gamma_2 (1 + gamma_3) + gamma_4), about 12 u.
//   2. c = fl(sr pr + si pi) and s = fl(si pr - sr pi) are each within gamma_2 (|sr pr| + |si pi|) <= gamma_2 |r_s| |r_p|
//      of the exact value (Cauchy-Schwarz), so the pair is within sqrt(2) gamma_2 |r_s| |r_p|. On the stored floats
//      c^2 + s^2 = |r_s|^2 |r_p|^2 (Lagrange), so |(c, s)|^2 <= (1 + 12 u)^2 (1 + sqrt(2) gamma_2)^2, about 1 + 29.7 u.
//   3. The rebuild: the ratio rounds once and its square twice, the fma once and the two products twice, gamma_5 in all,
//      so small <= (c^2 + s^2) (1 + gamma_5), about 1 + 34.7 u: k = 36.
// The matrices that reach the call carry two more roundings: one for the division by `large` when large != 1, and the
// callers' weight scaling v -> fl(fl(v / m) w) of every entry (transport.rs, coated_diffuse.rs), which puts gamma_2 on
// each entry, squares it in |cross|^2 and shows once more in `large`: k = 36 + 1 + 7 = 44.
//
// k is derived for the dielectric and scaled paths only. The conductor thin-film builder (conductor/thin_film) reaches
// the same call (dielectric.rs `apply_passive`) and its k is not derived: measured over the fs-render lib on ts2,
// 2026-11-11, the largest (small - large) / large at each caller site was 10.0 u at transport.rs `scatter`, 4.0 u at the
// contact `apply_power`, 3.1 u at camera/sequential.rs and 0.0 u at `apply_passive` (small never exceeded large there).
const DIAGONAL_JONES_ROUNDING_K: f64 = 44.0;
const DIAGONAL_JONES_UNIT_ROUNDOFF: f64 = f64::EPSILON * 0.5;
const DIAGONAL_JONES_ALLOWANCE: f64 = 1.0
    + DIAGONAL_JONES_ROUNDING_K * DIAGONAL_JONES_UNIT_ROUNDOFF
        / (1.0 - DIAGONAL_JONES_ROUNDING_K * DIAGONAL_JONES_UNIT_ROUNDOFF);

/// A matrix given to [`Mueller::apply_diagonal_jones`] is not a nondepolarizing diagonal Jones response.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum DiagonalJonesError {
    /// The weak power rebuilt from the cross entries exceeds the strong one by more than rounding explains
    /// (`weak > strong * (1 + gamma_44)`). In exact arithmetic that is impossible, so the matrix is malformed.
    WeakChannelExceedsStrong {
        /// `|cross|^2 / strong` as rebuilt.
        weak: f64,
        /// The strong power, `mean + |difference|`.
        strong: f64,
    },
}

impl core::fmt::Display for DiagonalJonesError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::WeakChannelExceedsStrong { weak, strong } => write!(
                f,
                "diagonal Jones matrix: weak channel {weak:e} exceeds strong channel {strong:e} by more than rounding"
            ),
        }
    }
}

impl core::error::Error for DiagonalJonesError {}

impl Mueller {
    pub const IDENTITY: Self = Self([
        [1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0],
    ]);

    /// Passive ideal linear polarizer with intensity throughput 1/2
    /// for unpolarized light, oriented at angle from the Q-positive axis.
    pub fn linear_polarizer(angle: f64) -> Self {
        let (s, c) = (2.0 * angle).sin_cos();
        let a = [1.0, c, s, 0.0];
        Self(core::array::from_fn(|r| core::array::from_fn(|k| 0.5 * a[r] * a[k])))
    }

    /// Ideal lossless linear retarder with phase difference in radians.
    /// At zero orientation it leaves Q unchanged and rotates the U/V plane.
    pub fn linear_retarder(angle: f64, phase: f64) -> Self {
        let (s, c) = (2.0 * angle).sin_cos();
        let (sp, cp) = phase.sin_cos();
        let n = [c, s];
        let a = [[cp + (1.0-cp)*n[0]*n[0], (1.0-cp)*n[0]*n[1]],
                 [(1.0-cp)*n[1]*n[0], cp + (1.0-cp)*n[1]*n[1]]];
        Self([
            [1.0, 0.0, 0.0, 0.0],
            [0.0, a[0][0], a[0][1], -s*sp],
            [0.0, a[1][0], a[1][1], c*sp],
            [0.0, s*sp, -c*sp, cp],
        ])
    }

    pub fn apply(self, input: Stokes) -> Stokes {
        let x = [input.i, input.q, input.u, input.v];
        let y = self.0.map(|row| row.iter().zip(x).map(|(a,b)| a*b).sum::<f64>());
        Stokes { i: y[0], q: y[1], u: y[2], v: y[3] }
    }

    /// Apply a nondepolarizing diagonal Jones response in its s/p frame.
    ///
    /// Native Fresnel matrices retain `r_s conj(r_p)` in the U/V block.
    /// Recover the weaker power from this product instead of subtracting
    /// rounded `(R_s + R_p)/2` and `(R_s - R_p)/2`. That subtraction erases
    /// representable p reflection near Brewster incidence, and can turn a
    /// strongly polarized coherent-film transmission into exact extinction.
    ///
    /// This is only for the native coherent interface matrices (including a
    /// scalar transport weight). Arbitrary or depolarizing Mueller matrices
    /// must use `apply`; their U/V block does not determine s/p powers.
    ///
    /// # Errors
    /// [`DiagonalJonesError::WeakChannelExceedsStrong`] when the rebuilt weak power exceeds the strong one by more
    /// than the rounding allowance derived at `DIAGONAL_JONES_ROUNDING_K`. A smaller excess is that rounding and
    /// the weak power is set to the strong one, so no unpolarized input is published above its own power.
    pub fn apply_diagonal_jones(self, input: Stokes) -> Result<Stokes, DiagonalJonesError> {
        let mean = self.0[0][0];
        let difference = self.0[0][1];
        let large = mean + difference.abs();
        if large == 0.0 {
            return Ok(Stokes::unpolarized(0.0));
        }
        let cross_re = self.0[2][2];
        let cross_im = self.0[3][2];
        let cross_scale = cross_re.abs().max(cross_im.abs());
        let mut small = if cross_scale == 0.0 {
            0.0
        } else {
            let ratio = cross_re.abs().min(cross_im.abs()) / cross_scale;
            // Divide before multiplying: |cross|² can underflow while the
            // weaker power remains representable. Scaling also avoids sqrt
            // and platform-libm hypot, retaining deterministic IEEE ops.
            cross_scale * ((cross_scale / large) * ratio.mul_add(ratio, 1.0))
        };
        if small > large {
            if small > large * DIAGONAL_JONES_ALLOWANCE {
                return Err(DiagonalJonesError::WeakChannelExceedsStrong {
                    weak: small,
                    strong: large,
                });
            }
            small = large;
        }
        let (s_power, p_power) = if difference >= 0.0 {
            (large, small)
        } else {
            (small, large)
        };
        let half_sum = |a: f64, b: f64| {
            let sum = a + b;
            if sum.is_finite() { 0.5 * sum } else { 0.5 * a + 0.5 * b }
        };
        let s = s_power * half_sum(input.i, input.q);
        let p = p_power * half_sum(input.i, -input.q);
        Ok(Stokes {
            i: s + p,
            q: s - p,
            u: cross_re * input.u - cross_im * input.v,
            v: cross_im * input.u + cross_re * input.v,
        })
    }

    /// Composition: self acts after rhs.
    pub fn compose(self, rhs: Self) -> Self {
        Self(core::array::from_fn(|r| core::array::from_fn(|c|
            (0..4).map(|k| self.0[r][k] * rhs.0[k][c]).sum()
        )))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn near(a: f64, b: f64) { assert!((a-b).abs() < 1e-12, "{a} != {b}"); }

    #[test]
    fn malus_law_and_crossed_polarizers() {
        let source = Stokes::unpolarized(10.0);
        let horizontal = Mueller::linear_polarizer(0.0).apply(source);
        near(horizontal.i, 5.0);
        near(horizontal.degree(), 1.0);
        near(Mueller::linear_polarizer(core::f64::consts::FRAC_PI_2).apply(horizontal).i, 0.0);
        near(Mueller::linear_polarizer(core::f64::consts::FRAC_PI_4).apply(horizontal).i, 2.5);
    }

    #[test]
    fn retarder_preserves_intensity_and_polarization_degree() {
        let source = Stokes::try_new(2.0, 1.0, 1.0, 0.0).unwrap();
        let out = Mueller::linear_retarder(0.0, core::f64::consts::FRAC_PI_2).apply(source);
        near(out.i, source.i);
        near(out.q, 1.0);
        near(out.u, 0.0);
        near(out.v, -1.0);
        near(out.degree(), source.degree());
    }

    #[test]
    fn rotations_and_composition() {
        let x = Stokes::try_new(1.0, 1.0, 0.0, 0.0).unwrap();
        let rotated = x.rotate_axes(core::f64::consts::FRAC_PI_4);
        near(rotated.q, 0.0);
        near(rotated.u, -1.0);
        let a = Mueller::linear_polarizer(0.3);
        let b = Mueller::linear_retarder(0.7, 0.8);
        let lhs = a.compose(b).apply(x);
        let rhs = a.apply(b.apply(x));
        near(lhs.i, rhs.i); near(lhs.q, rhs.q); near(lhs.u, rhs.u); near(lhs.v, rhs.v);
    }

    #[test]
    fn invalid_stokes_states_are_rejected() {
        assert!(Stokes::try_new(1.0, 1.1, 0.0, 0.0).is_none());
        assert!(Stokes::try_new(-1.0, 0.0, 0.0, 0.0).is_none());
        assert!(Stokes::try_new(0.0, 0.0, 0.0, 1.0).is_none());
        assert!(Stokes::try_new(f64::NAN, 0.0, 0.0, 0.0).is_none());
    }

    #[test]
    fn g1_diagonal_jones_retains_weak_channels_and_complex_phase() {
        // Independent Jones amplitudes s=0.4+0.2i, p=1e-10+2e-10i.
        // Their powers differ by 4e18, so the ordinary Mueller I/Q block
        // rounds the weaker channel out of both its sum and difference.
        let (rs, rp, re, im) = (0.2, 5e-20, 8e-11, -6e-11);
        for swapped in [false, true] {
            let (s, p, phase_im) = if swapped { (rp, rs, -im) } else { (rs, rp, im) };
            let mean = 0.5 * s + 0.5 * p;
            let difference = 0.5 * s - 0.5 * p;
            let matrix = Mueller([
                [mean, difference, 0.0, 0.0],
                [difference, mean, 0.0, 0.0],
                [0.0, 0.0, re, -phase_im],
                [0.0, 0.0, phase_im, re],
            ]);
            let q = if swapped { 1.0 } else { -1.0 };
            let weak = Stokes::try_new(1.0, q, 0.0, 0.0).unwrap();
            assert_eq!(matrix.apply(weak).i, 0.0, "regression must expose cancellation");
            let actual = matrix.apply_diagonal_jones(weak).unwrap();
            assert!((actual.i / rp - 1.0).abs() < 8.0 * f64::EPSILON);
            assert_eq!(actual.q, q * actual.i);
            let mixed = matrix
                .apply_diagonal_jones(Stokes::try_new(1.0, 0.0, 0.6, 0.8).unwrap())
                .unwrap();
            near(mixed.i, mean);
            near(mixed.q, difference);
            assert!((mixed.u - (0.6 * re - 0.8 * phase_im)).abs() < 1e-25);
            assert!((mixed.v - (0.6 * phase_im + 0.8 * re)).abs() < 1e-25);
        }
    }

    #[test]
    fn g0_diagonal_jones_keeps_subnormal_survival_without_squaring_the_cross_term() {
        let s = 1e-300_f64;
        let p = 1e-320_f64;
        let cross = s.sqrt() * p.sqrt();
        let mean = 0.5 * s + 0.5 * p;
        let difference = 0.5 * s - 0.5 * p;
        let matrix = Mueller([
            [mean, difference, 0.0, 0.0],
            [difference, mean, 0.0, 0.0],
            [0.0, 0.0, cross, 0.0],
            [0.0, 0.0, 0.0, cross],
        ]);
        let out = matrix
            .apply_diagonal_jones(Stokes {
                i: 1.0,
                q: -1.0,
                u: 0.0,
                v: 0.0,
            })
            .unwrap();
        assert!(out.i > 0.0);
        assert!((out.i - p).abs() <= f64::from_bits(2));
        assert_eq!(out.q, -out.i);
        // Each squared cross component alone rounds to zero. Combining the
        // scale factor before the last product retains the least positive power.
        let minimum = f64::from_bits(1);
        let component = minimum.sqrt() * core::f64::consts::FRAC_1_SQRT_2;
        let complex = Mueller([
            [0.5, 0.5, 0.0, 0.0],
            [0.5, 0.5, 0.0, 0.0],
            [0.0, 0.0, component, -component],
            [0.0, 0.0, component, component],
        ]);
        assert_eq!(
            complex
                .apply_diagonal_jones(Stokes {
                    i: 1.0,
                    q: -1.0,
                    u: 0.0,
                    v: 0.0
                })
                .unwrap()
                .i,
            minimum
        );
        let identity = Mueller::IDENTITY;
        for i in [f64::from_bits(1), 1.0, f64::MAX] {
            let state = Stokes { i, q: -i, u: 0.0, v: 0.0 };
            assert_eq!(identity.apply_diagonal_jones(state).unwrap(), state);
        }
        assert_eq!(
            Mueller([[0.0; 4]; 4])
                .apply_diagonal_jones(Stokes::unpolarized(1.0))
                .unwrap(),
            Stokes::unpolarized(0.0)
        );
    }

    // Two unit powers and a cross entry c + 0i: in exact arithmetic the weak channel is c^2 and a lossless
    // interface has c = 1.
    fn unit_power_matrix(cross: f64) -> Mueller {
        Mueller([
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, cross, 0.0],
            [0.0, 0.0, 0.0, cross],
        ])
    }

    #[test]
    fn g0_diagonal_jones_sets_a_rounding_excess_of_the_weak_channel_to_the_strong_one() {
        // |cross|^2 = 1 + 2^-51, as the contact total internal reflection rebuilt it.
        let cross = 1.0 + f64::EPSILON;
        assert!(
            cross * cross > 1.0,
            "the rebuilt weak power exceeds the strong one by rounding alone"
        );
        let out = unit_power_matrix(cross)
            .apply_diagonal_jones(Stokes::unpolarized(1.0))
            .unwrap();
        assert_eq!(out, Stokes::unpolarized(1.0));
        // The allowance is 44 u = 22 * 2^-52 above 1 (DIAGONAL_JONES_ROUNDING_K): cross = 1 + m 2^-52 rebuilds
        // 1 + 2m 2^-52, so m = 10 is inside and m = 12 is outside; m = 11 sits on the edge and is not used.
        let inside = 1.0 + 10.0 * f64::EPSILON;
        let out = unit_power_matrix(inside)
            .apply_diagonal_jones(Stokes::unpolarized(1.0))
            .unwrap();
        assert_eq!(out, Stokes::unpolarized(1.0));
        let polarized = Stokes::try_new(1.0, -0.5, 0.0, 0.0).unwrap();
        let out = unit_power_matrix(inside)
            .apply_diagonal_jones(polarized)
            .unwrap();
        // Strong and weak powers both 1 act as the identity on (I, Q): Q' = P_s (I + Q)/2 - P_p (I - Q)/2 = Q.
        assert_eq!((out.i, out.q), (1.0, -0.5));
    }

    #[test]
    fn g0_diagonal_jones_refuses_a_weak_channel_above_the_strong_one_beyond_rounding() {
        // Planted negative: just above the allowance (m = 12), and a weak channel 1e-3 above the strong one.
        let outside = 1.0 + 12.0 * f64::EPSILON;
        assert_eq!(
            unit_power_matrix(outside).apply_diagonal_jones(Stokes::unpolarized(1.0)),
            Err(DiagonalJonesError::WeakChannelExceedsStrong {
                weak: 1.0 + 24.0 * f64::EPSILON,
                strong: 1.0
            })
        );
        let gross = (1.0_f64 + 1e-3).sqrt();
        let refused = unit_power_matrix(gross).apply_diagonal_jones(Stokes::unpolarized(1.0));
        let Err(DiagonalJonesError::WeakChannelExceedsStrong { weak, strong }) = refused else {
            panic!("a weak channel 1e-3 above the strong one must be refused, got {refused:?}");
        };
        assert_eq!(strong, 1.0);
        near(weak, 1.001);
    }

    #[cfg(feature = "tracer")]
    #[test]
    fn g1_total_internal_reflection_never_publishes_more_than_the_input_power() {
        use super::dielectric::fresnel_dielectric_polarized;
        // The contact case: index 1.5 to 1.0 at cosine 0.6 read 1 + 2^-52 for unit input before the cap.
        let reflection = fresnel_dielectric_polarized(0.6, 1.5, 1.0)
            .unwrap()
            .reflection;
        assert_eq!(
            reflection
                .apply_diagonal_jones(Stokes::unpolarized(1.0))
                .unwrap(),
            Stokes::unpolarized(1.0)
        );
        // A grid over the whole total-internal-reflection range: none is refused and none exceeds 1.
        let mut reflected = 0;
        for transmitted in [1.0_f64, 1.33, 1.5] {
            for ratio in [1.02_f64, 1.1, 1.5, 2.0, 3.0] {
                let incident = transmitted * ratio;
                let critical = (1.0 - (transmitted / incident).powi(2)).sqrt();
                for step in 1..200 {
                    let cosine = critical * 0.999 * f64::from(step) / 200.0;
                    let fresnel =
                        fresnel_dielectric_polarized(cosine, incident, transmitted).unwrap();
                    if !fresnel.scalar.total_internal_reflection {
                        continue;
                    }
                    reflected += 1;
                    let out = fresnel
                        .reflection
                        .apply_diagonal_jones(Stokes::unpolarized(1.0))
                        .unwrap();
                    assert!(
                        out.i <= 1.0,
                        "{incident}/{transmitted} at {cosine}: {}",
                        out.i
                    );
                }
            }
        }
        assert!(
            reflected > 2000,
            "the grid must reach total internal reflection, got {reflected}"
        );
    }
}

/// Polarization-preserving optics connected to calibrated fixture output.
#[cfg(feature = "photometry")]
pub mod photometry {
    //! Polarization optics connected to calibrated fixture spectra and lumens.
    //!
    //! A chain retains positive coherency factors in each existing 380--780 nm
    //! bin and exposes their Stokes state. Independent populations never interfere.
    //! Only its final intensity is projected into a native `PhotometricSource`,
    //! so crossed polarizers extinguish light instead of independently halving it.
    //! The original source spectrum/flux remains the normalization reference.
    //!
    //! This is a uniform, paraxial source-optics model in one transverse frame,
    //! not a spatial ray-optics solve. Optics do not change placement or angular
    //! distribution. The caller places each outgoing Fresnel channel. Birefringent
    //! waveplates model phase delay only, without walk-off or interface reflections.
    //! The native intensity source preserves spectral power, but does NOT make the
    //! existing scalar path tracer propagate Q/U/V after emission.

    use core::f64::consts::{PI, TAU};

    use crate::conductor::{CoherentCoatingStack, ConductorOptics};
    use crate::dielectric::{DielectricCoatingStack, DielectricEvent, DielectricIor};
    use crate::photometry::source::{BIN_WIDTH_NM, PhotometricSource, SPECTRAL_BINS, SourceError, wavelength_nm};
    use crate::polarization::dielectric::spectral_fresnel_dielectric_polarized_borrowed;
    use crate::polarization::source_coherency::{JonesResponse, SourceCoherency};
    use crate::polarization::{Mueller, Stokes};
    use fs_math::c64::C64;

    /// Precomputed passive optics at the existing photometric quadrature bins.
    /// Construction is restricted to physical optical elements, not arbitrary
    /// caller-authored Mueller matrices. No per-ray matrix work is added.
    #[derive(Clone, Debug)]
    pub struct SpectralOptics {
        responses: [JonesResponse; SPECTRAL_BINS],
        #[cfg(test)]
        matrices: [Mueller; SPECTRAL_BINS],
    }

    impl SpectralOptics {
        /// Linear polarizer with distinct power transmissions along its principal
        /// axis and perpendicular axis. Both arrays must lie in `[0,1]`.
        /// This includes finite extinction ratio and wavelength-dependent loss.
        pub fn linear_diattenuator(
            angle_rad: f64,
            parallel: [f64; SPECTRAL_BINS],
            perpendicular: [f64; SPECTRAL_BINS],
        ) -> Result<Self, SourceError> {
            if !angle_rad.is_finite() || parallel.iter().chain(&perpendicular)
                .any(|t| !t.is_finite() || !(0.0..=1.0).contains(t)) {
                return Err(SourceError::InvalidFilter);
            }
            let responses = core::array::from_fn(|bin| {
                let (p, t) = (parallel[bin], perpendicular[bin]);
                if p == t { JonesResponse::neutral(p) } else {
                    JonesResponse::principal(angle_rad, [C64::from_re(p.sqrt()), C64::from_re(t.sqrt())])
                }
            });
            #[cfg(test)]
            let matrices = core::array::from_fn(|bin| {
                let (s, c) = fs_math::det::sin_cos(2.0 * angle_rad.rem_euclid(PI));
                let (p, t) = (parallel[bin], perpendicular[bin]);
                if p == t {
                    return Mueller(Mueller::IDENTITY.0.map(|row| row.map(|value| value * p)));
                }
                let a = 0.5 * p + 0.5 * t;
                let b = 0.5 * p - 0.5 * t;
                let g = p.sqrt() * t.sqrt();
                Mueller([
                    [a, b * c, b * s, 0.0],
                    [b * c, a * c * c + g * s * s, (a - g) * c * s, 0.0],
                    [b * s, (a - g) * c * s, a * s * s + g * c * c, 0.0],
                    [0.0, 0.0, 0.0, g],
                ])
            });
            Ok(Self { responses, #[cfg(test)] matrices })
        }

        /// Ideal achromatic polarizer. Unpolarized light loses half its power.
        pub fn linear_polarizer(angle_rad: f64) -> Result<Self, SourceError> {
            Self::linear_diattenuator(angle_rad, [1.0; SPECTRAL_BINS], [0.0; SPECTRAL_BINS])
        }

        /// A neutral spectral filter attenuates all polarization components equally.
        pub fn neutral(transmission: [f64; SPECTRAL_BINS]) -> Result<Self, SourceError> {
            Self::linear_diattenuator(0.0, transmission, transmission)
        }

        /// Lossless linear retarder. Positive phase delays the perpendicular
        /// component relative to the principal axis, matching `Mueller` conventions.
        pub fn linear_retarder(
            angle_rad: f64,
            phase_rad: [f64; SPECTRAL_BINS],
        ) -> Result<Self, SourceError> {
            if !angle_rad.is_finite() || phase_rad.iter().any(|phase| !phase.is_finite()) {
                return Err(SourceError::InvalidFilter);
            }
            // Reduce before doubling to avoid overflow for finite authored angles.
            let angle = angle_rad.rem_euclid(PI);
            let responses = phase_rad.map(|phase| {
                let (s, c) = fs_math::det::sin_cos(phase.rem_euclid(TAU));
                JonesResponse::principal(angle, [C64::ONE, C64::new(c, s)])
            });
            #[cfg(test)]
            let matrices = phase_rad.map(|phase| Mueller::linear_retarder(angle, phase.rem_euclid(TAU)));
            Ok(Self { responses, #[cfg(test)] matrices })
        }

        /// Phase-only, normal-incidence waveplate using two native spectral index
        /// laws. `thickness_m` is metres; the indices are the two principal modes.
        /// Interface losses and walk-off are deliberately not inferred.
        pub fn waveplate(
            angle_rad: f64,
            principal: impl Into<DielectricIor>,
            perpendicular: impl Into<DielectricIor>,
            thickness_m: f64,
        ) -> Result<Self, SourceError> {
            if !thickness_m.is_finite() || thickness_m < 0.0 {
                return Err(SourceError::InvalidFilter);
            }
            let principal = principal.into();
            let perpendicular = perpendicular.into();
            let mut phase = [0.0; SPECTRAL_BINS];
            for (bin, value) in phase.iter_mut().enumerate() {
                let wavelength = wavelength_nm(bin);
                let a = principal
                    .eval(wavelength)
                    .map_err(|_| SourceError::InvalidFilter)?;
                let b = perpendicular
                    .eval(wavelength)
                    .map_err(|_| SourceError::InvalidFilter)?;
                *value = TAU * ((b - a) * thickness_m / (wavelength * 1e-9));
                if !value.is_finite() {
                    return Err(SourceError::NumericalRange);
                }
            }
            Self::linear_retarder(angle_rad, phase)
        }

        /// One Fresnel power channel at a fixed incidence angle, in s/p frames.
        /// Reflection and transmission can be applied to separate copies of the
        /// original source; the caller owns the corresponding outgoing placement.
        /// Both indices use the native Cauchy, Sellmeier or tabulated phase law.
        pub fn dielectric_channel(
            incident_cosine: f64,
            incident: impl Into<DielectricIor>,
            transmitted: impl Into<DielectricIor>,
            event: DielectricEvent,
        ) -> Result<Self, SourceError> {
            let incident = incident.into();
            let transmitted = transmitted.into();
            let mut matrices = [Mueller::IDENTITY; SPECTRAL_BINS];
            for (bin, matrix) in matrices.iter_mut().enumerate() {
                let f = spectral_fresnel_dielectric_polarized_borrowed(
                    incident_cosine,
                    &incident,
                    &transmitted,
                    wavelength_nm(bin),
                )
                .map_err(|_| SourceError::InvalidFilter)?;
                *matrix = match event {
                    DielectricEvent::Reflection => f.reflection,
                    DielectricEvent::Transmission => f.transmission,
                };
            }
            Ok(Self::from_diagonal_mueller(matrices))
        }

        /// Reflected power from the original metal optics, including any authored
        /// transparent, absorbing or periodic coating. The incident medium retains
        /// its native spectral phase law; no constant-index fit is introduced.
        ///
        /// Incoming and reflected frames use +Q along s and right-handed p axes,
        /// as in `ConductorOptics::fresnel_mueller`. The retained U/V coupling
        /// includes reflection phase. The caller places the outgoing source.
        pub fn conductor_reflection(
            incident_cosine: f64,
            incident: impl Into<DielectricIor>,
            optics: &ConductorOptics,
        ) -> Result<Self, SourceError> {
            let incident = incident.into();
            let mut matrices = [Mueller::IDENTITY; SPECTRAL_BINS];
            for (bin, matrix) in matrices.iter_mut().enumerate() {
                let wavelength = wavelength_nm(bin);
                let eta = incident.eval(wavelength).map_err(|_| SourceError::InvalidFilter)?;
                *matrix = optics.fresnel_mueller(wavelength, eta, incident_cosine)
                    .map_err(|_| SourceError::InvalidFilter)?;
            }
            Ok(Self::from_diagonal_mueller(matrices))
        }

        /// Reflected power from an original ordered spectral stack on a metal.
        /// The borrowed layers are incident-side first and replace any coating
        /// already attached to `substrate`, matching the native stack owner.
        /// Full spectral data and s/p phase are evaluated only during construction;
        /// the resulting optic owns fixed-size responses and borrows no layer data.
        /// Frames and outgoing placement follow `Self::conductor_reflection`.
        pub fn conductor_stack_reflection(
            incident_cosine: f64,
            incident: impl Into<DielectricIor>,
            substrate: &ConductorOptics,
            stack: CoherentCoatingStack<'_>,
        ) -> Result<Self, SourceError> {
            let incident = incident.into();
            let mut matrices = [Mueller::IDENTITY; SPECTRAL_BINS];
            for (bin, matrix) in matrices.iter_mut().enumerate() {
                let wavelength = wavelength_nm(bin);
                let eta = incident.eval(wavelength).map_err(|_| SourceError::InvalidFilter)?;
                *matrix = stack.fresnel_mueller(substrate, wavelength, eta, incident_cosine)
                    .map_err(|_| SourceError::InvalidFilter)?;
            }
            Ok(Self::from_diagonal_mueller(matrices))
        }

        /// One coherent multilayer reflection/transmission channel between the
        /// actual adjacent spectral media. `reverse` traverses the authored layers
        /// inside-to-outside; callers also supply the corresponding incident and
        /// transmitted media and incidence cosine. Layers retain their own support.
        ///
        /// These are native ray-power Mueller matrices, including coating loss
        /// and relative phase. No radiance eta-squared factor or stochastic branch
        /// probability is applied. Outgoing placement follows `dielectric_channel`.
        pub fn coated_dielectric_channel(
            incident_cosine: f64,
            incident: impl Into<DielectricIor>,
            transmitted: impl Into<DielectricIor>,
            stack: DielectricCoatingStack<'_>,
            reverse: bool,
            event: DielectricEvent,
        ) -> Result<Self, SourceError> {
            let incident = incident.into();
            let transmitted = transmitted.into();
            let mut matrices = [Mueller::IDENTITY; SPECTRAL_BINS];
            for (bin, matrix) in matrices.iter_mut().enumerate() {
                let wavelength = wavelength_nm(bin);
                let eta_i = incident.eval(wavelength).map_err(|_| SourceError::InvalidFilter)?;
                let eta_t = transmitted.eval(wavelength).map_err(|_| SourceError::InvalidFilter)?;
                let response = stack.fresnel_polarized(
                    wavelength, reverse, incident_cosine, eta_i, eta_t,
                ).map_err(|_| SourceError::InvalidFilter)?;
                *matrix = match event {
                    DielectricEvent::Reflection => response.reflection,
                    DielectricEvent::Transmission => response.transmission,
                };
            }
            Ok(Self::from_diagonal_mueller(matrices))
        }

        fn from_diagonal_mueller(matrices: [Mueller; SPECTRAL_BINS]) -> Self {
            Self { responses: matrices.map(JonesResponse::from_diagonal_mueller), #[cfg(test)] matrices }
        }
    }

    /// A calibrated source plus its retained wavelength-dependent polarization.
    #[derive(Clone, Debug)]
    pub struct PolarizedSource {
        reference: PhotometricSource,
        intensity: PhotometricSource,
        relative: [Stokes; SPECTRAL_BINS],
        coherency: [SourceCoherency; SPECTRAL_BINS],
    }

    impl PolarizedSource {
        /// Wrap an existing calibrated source without changing its emitted flux.
        #[must_use]
        pub fn unpolarized(source: PhotometricSource) -> Self {
            Self { intensity: source.clone(), reference: source,
                relative: [Stokes::unpolarized(1.0); SPECTRAL_BINS],
                coherency: [SourceCoherency::UNPOLARIZED; SPECTRAL_BINS] }
        }

        /// Set normalized per-bin Stokes states. I must be exactly one so authored
        /// polarization cannot silently replace the source's calibrated brightness.
        pub fn try_new(
            source: PhotometricSource,
            normalized: [Stokes; SPECTRAL_BINS],
        ) -> Result<Self, SourceError> {
            if normalized.iter().any(|s| s.i != 1.0
                || Stokes::try_new(s.i, s.q, s.u, s.v).is_none()) {
                return Err(SourceError::InvalidSpectrum);
            }
            Ok(Self { intensity: source.clone(), reference: source, relative: normalized,
                coherency: normalized.map(SourceCoherency::from_stokes) })
        }

        /// Apply the next optic without discarding Q/U/V or restoring lost flux.
        /// Uses bounded fixed-size work; failures leave the original unchanged.
        pub fn through(&self, optic: &SpectralOptics) -> Result<Self, SourceError> {
            let mut relative = self.relative;
            let mut coherency = self.coherency;
            for (bin, state) in relative.iter_mut().enumerate() {
                coherency[bin] = coherency[bin].through(optic.responses[bin]);
                let transmitted = coherency[bin].stokes();
                *state = passive_roundoff(transmitted, state.i)?;
                if transmitted.i > state.i {
                    // Keep the retained factors consistent with a roundoff-only
                    // passivity correction, without discarding their weak mode.
                    coherency[bin].scale_power(state.i / transmitted.i);
                }
            }
            // Project absolute surviving power only once relative to the ORIGINAL
            // spectrum. Normalizing Stokes after each element would erase extinction.
            let intensity = self.reference.filtered(relative.map(|s| s.i))?;
            Ok(Self { reference: self.reference.clone(), intensity, relative, coherency })
        }

        /// Per-bin Stokes power relative to the original source's spectral power.
        /// At dark bins the state is exactly zero, not an arbitrary normalized axis.
        #[must_use]
        pub const fn relative_stokes(&self) -> &[Stokes; SPECTRAL_BINS] { &self.relative }

        /// Emitted polarization fractions at a vacuum wavelength in nanometres.
        ///
        /// This uses the SAME piecewise-constant 380--780 nm bins as the
        /// calibrated intensity spectrum, including its upper endpoint. It
        /// does not interpolate fractions or restore optical loss. The caller
        /// uses `intensity_source()` for the surviving spectral intensity and
        /// declares the transverse frame in which these coordinates apply.
        /// A dark bin returns an unpolarized unit bookkeeping state; its
        /// calibrated emitted intensity remains zero.
        pub fn normalized_stokes(&self, wavelength_nm: f64) -> Result<Stokes, SourceError> {
            if !wavelength_nm.is_finite()
                || !(crate::spectral::LAMBDA_MIN..=crate::spectral::LAMBDA_MAX)
                    .contains(&wavelength_nm)
            {
                return Err(SourceError::InvalidSpectrum);
            }
            let bin = (((wavelength_nm - crate::spectral::LAMBDA_MIN) / BIN_WIDTH_NM)
                as usize).min(SPECTRAL_BINS - 1);
            let state = self.relative[bin];
            if state.i == 0.0 {
                Ok(Stokes::unpolarized(1.0))
            } else {
                Stokes::try_new(1.0, state.q / state.i, state.u / state.i, state.v / state.i)
                    .ok_or(SourceError::NumericalRange)
            }
        }

        /// Native calibrated intensity view, usable by existing illuminance,
        /// haze and spectral rendering consumers. Q/U/V is retained in this wrapper,
        /// not propagated by those scalar consumers.
        #[must_use]
        pub const fn intensity_source(&self) -> &PhotometricSource { &self.intensity }

        /// Explicitly finish source optics and hand calibrated intensity to a scene.
        #[must_use]
        pub fn into_intensity_source(self) -> PhotometricSource { self.intensity }
    }

    impl PhotometricSource {
        /// Start a polarization-preserving optical chain on an unpolarized fixture.
        #[must_use]
        pub fn with_unpolarized_optics(self) -> PolarizedSource { PolarizedSource::unpolarized(self) }
    }

    // Only private, physically constructed matrices reach this boundary. Project
    // floating-point roundoff at zero intensity/the Stokes cone; never admit a
    // general unphysical matrix or a meaningful power gain by clamping it away.
    fn passive_roundoff(mut s: Stokes, incident: f64) -> Result<Stokes, SourceError> {
        let tolerance = 64.0 * f64::EPSILON * incident;
        let polarized = s.q.hypot(s.u).hypot(s.v);
        if ![s.i, s.q, s.u, s.v].iter().all(|v| v.is_finite())
            || s.i < -tolerance || s.i > incident + tolerance
            || polarized > s.i.max(0.0) + tolerance {
            return Err(SourceError::NumericalRange);
        }
        s.i = s.i.clamp(0.0, incident);
        if s.i == 0.0 { return Ok(Stokes::unpolarized(0.0)); }
        if polarized > s.i {
            let scale = s.i / polarized;
            s.q *= scale;
            s.u *= scale;
            s.v *= scale;
        }
        Stokes::try_new(s.i, s.q, s.u, s.v).ok_or(SourceError::NumericalRange)
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use crate::conductor::{
            CoatingIorSample, CoatingLayer, CoherentAbsorbingCoating, CoherentDielectricCoating,
            ConductorDataStatus, ConductorIorSample, ConductorSource,
        };
        use crate::dielectric::{CauchyIor, DielectricIorSample, SellmeierIor, TabulatedIor};
        use crate::instances::RigidTransform;
        use crate::photometry::source::{AngularDistribution, SourceSpectrum};
        use fs_geom::Vec3;

        fn source() -> PhotometricSource {
            PhotometricSource::try_new(
                RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], [0.0; 3]).unwrap(),
                1000.0, 1.0, AngularDistribution::isotropic(),
                SourceSpectrum::try_from_bins([1.0; SPECTRAL_BINS]).unwrap(),
            ).unwrap()
        }

        fn close(a: f64, b: f64) { assert!((a - b).abs() < 2e-10, "{a} != {b}"); }

        #[test]
        fn source_polarization_lookup_uses_exact_calibrated_bins_without_restoring_loss() {
            let states = core::array::from_fn(|bin| match bin % 3 {
                0 => Stokes::try_new(1.0, 1.0, 0.0, 0.0).unwrap(),
                1 => Stokes::try_new(1.0, 0.0, 0.0, -1.0).unwrap(),
                _ => Stokes::try_new(1.0, 0.3, -0.4, 0.2).unwrap(),
            });
            let transmission = core::array::from_fn(|bin| {
                if bin == 7 { 0.0 } else { 0.2 + 0.6 * bin as f64 / SPECTRAL_BINS as f64 }
            });
            let original = source();
            let polarized = PolarizedSource::try_new(original.clone(), states).unwrap()
                .through(&SpectralOptics::neutral(transmission).unwrap()).unwrap();
            for bin in 0..SPECTRAL_BINS {
                let expected = if transmission[bin] == 0.0 {
                    Stokes::unpolarized(1.0)
                } else { states[bin] };
                for wavelength in [
                    crate::spectral::LAMBDA_MIN + bin as f64 * BIN_WIDTH_NM,
                    wavelength_nm(bin),
                    crate::spectral::LAMBDA_MIN + (bin as f64 + 0.9) * BIN_WIDTH_NM,
                ] {
                    let state = polarized.normalized_stokes(wavelength).unwrap();
                    assert_eq!(state.i, 1.0);
                    close(state.q, expected.q);
                    close(state.u, expected.u);
                    close(state.v, expected.v);
                    let direction = Vec3::new(0.0, 0.0, 1.0);
                    let before = original.spectral_intensity_w_sr_nm(direction, wavelength).unwrap();
                    let after = polarized.intensity_source()
                        .spectral_intensity_w_sr_nm(direction, wavelength).unwrap();
                    close(after, before * transmission[bin]);
                }
            }
            assert_eq!(polarized.normalized_stokes(780.0).unwrap(),
                polarized.normalized_stokes(wavelength_nm(SPECTRAL_BINS - 1)).unwrap());
            for wavelength in [379.0, 780.0_f64.next_up(), f64::NAN, f64::INFINITY] {
                assert_eq!(polarized.normalized_stokes(wavelength), Err(SourceError::InvalidSpectrum));
            }
            assert!(polarized.intensity_source().flux_lm() < original.flux_lm());
        }

        #[test]
        fn crossed_and_rotated_polarizers_change_native_flux_and_candela() {
            let original = source();
            let horizontal = original.clone().with_unpolarized_optics()
                .through(&SpectralOptics::linear_polarizer(0.0).unwrap()).unwrap();
            close(horizontal.intensity_source().flux_lm(), 500.0);
            close(horizontal.intensity_source().intensity_cd(Vec3::new(0.0, 0.0, 1.0)).unwrap(),
                original.intensity_cd(Vec3::new(0.0, 0.0, 1.0)).unwrap() * 0.5);
            let crossed = horizontal.through(&SpectralOptics::linear_polarizer(PI / 2.0).unwrap()).unwrap();
            close(crossed.intensity_source().flux_lm(), 0.0);
            let angled = horizontal.through(&SpectralOptics::linear_polarizer(PI / 4.0).unwrap()).unwrap();
            close(angled.into_intensity_source().flux_lm(), 250.0);
            close(original.flux_lm(), 1000.0);
        }

        #[test]
        fn retarder_retains_circular_state_for_the_next_analyzer() {
            let linear = [Stokes::try_new(1.0, 0.0, 1.0, 0.0).unwrap(); SPECTRAL_BINS];
            let quarter = SpectralOptics::linear_retarder(0.0, [PI / 2.0; SPECTRAL_BINS]).unwrap();
            let circular = PolarizedSource::try_new(source(), linear).unwrap().through(&quarter).unwrap();
            close(circular.intensity_source().flux_lm(), 1000.0);
            close(circular.relative_stokes()[0].v, -1.0);
            for angle in [0.0, 0.37, PI / 4.0] {
                let out = circular.through(&SpectralOptics::linear_polarizer(angle).unwrap()).unwrap();
                close(out.intensity_source().flux_lm(), 500.0);
            }
        }

        #[test]
        fn finite_extinction_is_not_an_ideal_polarizer_or_depolarizer() {
            let optic = SpectralOptics::linear_diattenuator(0.0,
                [0.8; SPECTRAL_BINS], [0.2; SPECTRAL_BINS]).unwrap();
            let first = source().with_unpolarized_optics().through(&optic).unwrap();
            close(first.relative_stokes()[0].i, 0.5);
            close(first.relative_stokes()[0].q, 0.3);
            let crossed = SpectralOptics::linear_diattenuator(PI / 2.0,
                [0.8; SPECTRAL_BINS], [0.2; SPECTRAL_BINS]).unwrap();
            close(first.through(&crossed).unwrap().intensity_source().flux_lm(), 160.0);
            let linear = [Stokes::try_new(1.0, 0.0, 1.0, 0.0).unwrap(); SPECTRAL_BINS];
            let out = PolarizedSource::try_new(source(), linear).unwrap().through(&optic).unwrap();
            close(out.relative_stokes()[0].u, 0.4);
            close(out.relative_stokes()[0].degree(), 1.0);
        }

        #[test]
        fn crossed_finite_extinction_retains_minor_incoherent_population() {
            for leakage in [1e-20, 1e-80, 1e-200] {
                let first = SpectralOptics::linear_diattenuator(0.0,
                    [1.0; SPECTRAL_BINS], [leakage; SPECTRAL_BINS]).unwrap();
                let crossed = SpectralOptics::linear_diattenuator(PI / 2.0,
                    [1.0; SPECTRAL_BINS], [leakage; SPECTRAL_BINS]).unwrap();
                let input = source().with_unpolarized_optics().through(&first).unwrap();
                // Public binary64 Stokes cannot represent this minor population.
                // Its positive coherency factor must survive subsequent optics.
                assert_eq!(input.relative_stokes()[0].i, input.relative_stokes()[0].q);
                let neutral = SpectralOptics::neutral([0.3; SPECTRAL_BINS]).unwrap();
                let retarder = SpectralOptics::linear_retarder(0.0, [1.7; SPECTRAL_BINS]).unwrap();
                let output = input.through(&neutral).unwrap().through(&retarder).unwrap()
                    .through(&crossed).unwrap();
                for state in output.relative_stokes() {
                    assert!((state.i / (0.3 * leakage) - 1.0).abs() < 1e-14);
                    assert!((state.q / state.i).abs() < 1e-14);
                    assert_eq!(state.u, 0.0);
                    assert_eq!(state.v, 0.0);
                }
                assert!((output.intensity_source().flux_lm() / (300.0 * leakage) - 1.0)
                    .abs() < 1e-13);
                let direction = Vec3::new(0.0, 0.0, 1.0);
                assert!((output.intensity_source().intensity_cd(direction).unwrap()
                    / (300.0 * leakage / (4.0 * PI)) - 1.0).abs() < 1e-13);
            }
        }

        #[test]
        fn source_coherency_agrees_with_mueller_for_mixed_and_circular_light() {
            let optics = [
                SpectralOptics::linear_diattenuator(0.37,
                    [0.8; SPECTRAL_BINS], [0.23; SPECTRAL_BINS]).unwrap(),
                SpectralOptics::linear_retarder(-0.2, [1.4; SPECTRAL_BINS]).unwrap(),
                SpectralOptics::conductor_reflection(0.45,
                    CauchyIor::try_constant(1.0).unwrap(), &metal()).unwrap(),
                SpectralOptics::linear_diattenuator(-0.91,
                    [0.17; SPECTRAL_BINS], [0.71; SPECTRAL_BINS]).unwrap(),
            ];
            for (q, u, v) in [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0),
                (-1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0),
                (0.0, 0.0, -1.0), (0.2, -0.4, 0.5)] {
                let mut expected = [Stokes::try_new(1.0, q, u, v).unwrap(); SPECTRAL_BINS];
                let mut output = PolarizedSource::try_new(source(), expected).unwrap();
                for optic in &optics {
                    output = output.through(optic).unwrap();
                    for (bin, state) in expected.iter_mut().enumerate() {
                        *state = optic.matrices[bin].apply(*state);
                        let actual = output.relative_stokes()[bin];
                        for (a, b) in [actual.i, actual.q, actual.u, actual.v]
                            .into_iter().zip([state.i, state.q, state.u, state.v]) {
                            assert!((a - b).abs() < 2e-13);
                        }
                        assert!(actual.degree() <= 1.0 + 16.0 * f64::EPSILON);
                    }
                }
            }
        }

        #[test]
        fn neutral_filter_matches_native_spectral_calibration() {
            let transmission = core::array::from_fn(|bin| (bin + 1) as f64 / SPECTRAL_BINS as f64);
            let original = source();
            let expected = original.filtered(transmission).unwrap();
            let out = original.with_unpolarized_optics()
                .through(&SpectralOptics::neutral(transmission).unwrap()).unwrap();
            assert_eq!(out.intensity_source().flux_lm(), expected.flux_lm());
            for bin in 0..SPECTRAL_BINS {
                let wavelength = wavelength_nm(bin);
                assert_eq!(out.intensity_source().spectrum().value_w_per_lumen_nm(wavelength),
                    expected.spectrum().value_w_per_lumen_nm(wavelength));
            }
        }

        #[test]
        fn dispersive_fresnel_channels_close_calibrated_luminous_power() {
            let air = CauchyIor::try_constant(1.0).unwrap();
            let material = crate::dielectric::DielectricGlass::representative_crown();
            let sampled = DielectricIor::from(
                TabulatedIor::try_new(core::array::from_fn(|bin| {
                    DielectricIorSample::try_new(
                        380.0 + 50.0 * bin as f64,
                        1.6 - 0.005 * bin as f64,
                    )
                    .unwrap()
                }))
                .unwrap(),
            );
            let source = source().with_unpolarized_optics();
            for glass in [material.ior(), &sampled] {
                let r = source
                    .through(
                        &SpectralOptics::dielectric_channel(
                            0.6,
                            air,
                            glass,
                            DielectricEvent::Reflection,
                        )
                        .unwrap(),
                    )
                    .unwrap();
                let t = source
                    .through(
                        &SpectralOptics::dielectric_channel(
                            0.6,
                            air,
                            glass,
                            DielectricEvent::Transmission,
                        )
                        .unwrap(),
                    )
                    .unwrap();
                close(
                    r.intensity_source().flux_lm() + t.intensity_source().flux_lm(),
                    1000.0,
                );
                assert!(r.relative_stokes()[0].q > 0.0);
                assert_ne!(
                    r.relative_stokes()[0].i,
                    r.relative_stokes()[SPECTRAL_BINS - 1].i
                );
            }
        }

        #[test]
        fn g1_weak_brewster_channels_retain_calibrated_source_power() {
            let ni = 1.0_f64;
            let nt = 1.5_f64;
            let cosine = ni / ni.hypot(nt) + 1e-9;
            let transmitted_cosine =
                (1.0 - (ni / nt).powi(2) * (1.0 - cosine * cosine)).sqrt();
            // Independent electric-field Fresnel amplitude. Its power is
            // representable, but below the roundoff of the dominant s channel.
            let rp = (nt * cosine - ni * transmitted_cosine)
                / (nt * cosine + ni * transmitted_cosine);
            let expected = rp * rp;
            assert!(expected > 0.0 && expected < 1e-16);
            let air = CauchyIor::try_constant(ni).unwrap();
            let glass = CauchyIor::try_constant(nt).unwrap();
            let films = DielectricCoatingStack::try_new(&[]).unwrap();
            let input = PolarizedSource::try_new(
                source(), [Stokes::try_new(1.0, -1.0, 0.0, 0.0).unwrap(); SPECTRAL_BINS],
            ).unwrap();
            for optic in [
                SpectralOptics::dielectric_channel(
                    cosine, air, glass, DielectricEvent::Reflection,
                ).unwrap(),
                SpectralOptics::coated_dielectric_channel(
                    cosine, air, glass, films, false, DielectricEvent::Reflection,
                ).unwrap(),
            ] {
                let output = input.through(&optic).unwrap();
                for state in output.relative_stokes() {
                    assert!(state.i > 0.0);
                    assert!((state.i / expected - 1.0).abs() < 2e-6);
                    assert_eq!(state.q, -state.i);
                }
                assert!((output.intensity_source().flux_lm() / (1000.0 * expected) - 1.0)
                    .abs() < 2e-6);
                let analyzer = SpectralOptics::linear_polarizer(PI / 2.0).unwrap();
                let mixed = source().with_unpolarized_optics().through(&optic).unwrap()
                    .through(&analyzer).unwrap();
                assert!((mixed.intensity_source().flux_lm() / (500.0 * expected) - 1.0)
                    .abs() < 2e-6);
                // Reusing the same original source must not normalize this
                // very weak surviving signal back to the source's old flux.
                let twice = output.through(&optic).unwrap();
                assert!((twice.intensity_source().flux_lm() / (1000.0 * expected * expected)
                    - 1.0).abs() < 4e-6);
            }
            assert_eq!(input.intensity_source().flux_lm(), 1000.0);
        }

        #[test]
        fn g1_absorbing_brewster_minimum_survives_bare_and_stack_source_optics() {
            let eta = 1.5_f64;
            let k = 1e-9;
            let cosine = 1.0 / eta.hypot(1.0);
            let optics = ConductorOptics::try_new(core::array::from_fn(|i| {
                ConductorIorSample::try_new(380.0 + 50.0 * i as f64, eta, k).unwrap()
            }), metal().source()).unwrap();
            let stack = CoherentCoatingStack::try_new(&[]).unwrap();
            let air = CauchyIor::try_constant(1.0).unwrap();
            let input = PolarizedSource::try_new(
                source(), [Stokes::try_new(1.0, -1.0, 0.0, 0.0).unwrap(); SPECTRAL_BINS],
            ).unwrap();
            // First-order complex-index derivative of the p Fresnel amplitude
            // at Brewster incidence; the next relative term is O(k²).
            let amplitude = k * (eta * eta - 1.0) / (2.0 * eta.powi(3));
            let expected = amplitude * amplitude;
            for optic in [
                SpectralOptics::conductor_reflection(cosine, air, &optics).unwrap(),
                SpectralOptics::conductor_stack_reflection(cosine, air, &optics, stack).unwrap(),
            ] {
                let output = input.through(&optic).unwrap();
                for state in output.relative_stokes() {
                    assert!(state.i > 0.0);
                    assert!((state.i / expected - 1.0).abs() < 2e-6);
                    assert_eq!(state.q, -state.i);
                }
                assert!((output.intensity_source().flux_lm() / (1000.0 * expected) - 1.0)
                    .abs() < 2e-6);
                let analyzer = SpectralOptics::linear_polarizer(PI / 2.0).unwrap();
                let mixed = source().with_unpolarized_optics().through(&optic).unwrap()
                    .through(&analyzer).unwrap();
                assert!((mixed.intensity_source().flux_lm() / (500.0 * expected) - 1.0)
                    .abs() < 2e-6);
            }
        }

        #[test]
        fn waveplate_uses_principal_spectral_indices_and_metric_thickness() {
            let a = CauchyIor::try_constant(1.5).unwrap();
            let b = CauchyIor::try_constant(1.6).unwrap();
            let bin = 30;
            let thickness = wavelength_nm(bin) * 1e-9 / (4.0 * (1.6 - 1.5));
            let plate = SpectralOptics::waveplate(0.0, a, b, thickness).unwrap();
            let a_sellmeier = SellmeierIor::try_new([1.25, 0.0, 0.0], [0.0; 3]).unwrap();
            let b_sampled = DielectricIor::from(
                TabulatedIor::try_new(core::array::from_fn(|i| {
                    DielectricIorSample::try_new(380.0 + 50.0 * i as f64, 1.6).unwrap()
                }))
                .unwrap(),
            );
            assert_eq!(
                plate.matrices,
                SpectralOptics::waveplate(0.0, a_sellmeier, &b_sampled, thickness)
                    .unwrap()
                    .matrices
            );
            let linear = [Stokes::try_new(1.0, 0.0, 1.0, 0.0).unwrap(); SPECTRAL_BINS];
            let out = PolarizedSource::try_new(source(), linear)
                .unwrap()
                .through(&plate)
                .unwrap();
            close(out.relative_stokes()[bin].v, -1.0);
            close(out.intensity_source().flux_lm(), 1000.0);
            assert!(out.relative_stokes()[0].u.abs() > 0.01);
        }


        fn metal() -> ConductorOptics {
            let source = ConductorSource::try_new(
                fs_blake3::hash_domain("org.frankensim.test.photometric-reflection", b"synthetic"),
                ConductorDataStatus::Representative,
            ).unwrap();
            ConductorOptics::try_new(core::array::from_fn(|i| {
                ConductorIorSample::try_new(380.0 + 50.0 * i as f64, 0.2, 3.0).unwrap()
            }), source).unwrap()
        }

        fn mixed_films(source: ConductorSource) -> (CoherentDielectricCoating, CoherentAbsorbingCoating) {
            let transparent = CoherentDielectricCoating::try_new(
                SellmeierIor::try_new([1.56, 0.0, 0.0], [0.0; 3]).unwrap(), 81e-9, source,
            ).unwrap();
            let absorbing = CoherentAbsorbingCoating::try_new(core::array::from_fn(|i| {
                CoatingIorSample::try_new(
                    380.0 + 50.0 * i as f64, 1.9 - 0.025 * i as f64, 0.05 + 0.01 * i as f64,
                ).unwrap()
            }), 100e-9, source).unwrap();
            (transparent, absorbing)
        }

        #[test]
        fn metal_reflection_phase_changes_calibrated_analyzer_flux() {
            let input = PolarizedSource::try_new(source(),
                [Stokes::try_new(1.0, 0.0, 1.0, 0.0).unwrap(); SPECTRAL_BINS]).unwrap();
            let mirror = SpectralOptics::conductor_reflection(
                0.5, CauchyIor::try_constant(1.0).unwrap(), &metal(),
            ).unwrap();
            let reflected = input.through(&mirror).unwrap();
            // Independent complex Fresnel answer for n=0.2+3i at cos(theta)=0.5.
            let expected = [0.9184110846593692, 0.04400178561613738,
                -0.5749762549092700, 0.7148042176556463];
            for state in reflected.relative_stokes() {
                for (actual, expected) in [state.i, state.q, state.u, state.v].into_iter().zip(expected) {
                    close(actual, expected);
                }
            }
            close(reflected.intensity_source().flux_lm(), 1000.0 * expected[0]);
            let quarter = SpectralOptics::linear_retarder(0.0, [PI / 2.0; SPECTRAL_BINS]).unwrap();
            let analyzed = reflected.through(&quarter).unwrap();
            let plus = analyzed.through(&SpectralOptics::linear_polarizer(PI / 4.0).unwrap()).unwrap();
            let minus = analyzed.through(&SpectralOptics::linear_polarizer(-PI / 4.0).unwrap()).unwrap();
            close(plus.intensity_source().flux_lm(), 500.0 * (expected[0] + expected[3]));
            close(minus.intensity_source().flux_lm(), 500.0 * (expected[0] - expected[3]));
            close(input.intensity_source().flux_lm(), 1000.0);
        }

        #[test]
        fn conductor_stack_retains_spectral_phase_and_source_normalization() {
            let bare = metal();
            let (transparent, absorbing) = mixed_films(bare.source());
            let layers = [CoatingLayer::Dielectric(&transparent), CoatingLayer::Absorbing(&absorbing)];
            let stack = CoherentCoatingStack::try_new(&layers).unwrap();
            let incident = DielectricIor::from(TabulatedIor::try_new(core::array::from_fn(|i| {
                DielectricIorSample::try_new(380.0 + 50.0 * i as f64, 1.2 + 0.04 * i as f64).unwrap()
            })).unwrap());
            let optic = SpectralOptics::conductor_stack_reflection(0.6, &incident, &bare, stack).unwrap();
            assert_eq!(optic.matrices, SpectralOptics::conductor_stack_reflection(
                0.6, &incident, &bare.with_coating(transparent), stack,
            ).unwrap().matrices, "the explicit stack replaces the substrate's attached film");
            let state = Stokes::try_new(1.0, 0.3, 0.4, 0.2).unwrap();
            let original = source();
            let output = PolarizedSource::try_new(original.clone(), [state; SPECTRAL_BINS])
                .unwrap().through(&optic).unwrap();
            let surviving = core::array::from_fn(|bin| {
                let wavelength = wavelength_nm(bin);
                let expected = stack.fresnel_mueller(
                    &bare, wavelength, incident.eval(wavelength).unwrap(), 0.6,
                ).unwrap().apply(state);
                let actual = output.relative_stokes()[bin];
                for (actual, expected) in [actual.i, actual.q, actual.u, actual.v].into_iter()
                    .zip([expected.i, expected.q, expected.u, expected.v]) {
                    close(actual, expected);
                }
                expected.i
            });
            let calibrated = original.filtered(surviving).unwrap();
            close(output.intensity_source().flux_lm(), calibrated.flux_lm());
            for bin in 0..SPECTRAL_BINS {
                let wavelength = wavelength_nm(bin);
                close(output.intensity_source().flux_lm()
                    * output.intensity_source().spectrum().value_w_per_lumen_nm(wavelength),
                    calibrated.flux_lm() * calibrated.spectrum().value_w_per_lumen_nm(wavelength));
            }
            assert!((output.relative_stokes()[0].v - output.relative_stokes()[SPECTRAL_BINS - 1].v).abs() > 1e-4);
            let single = [CoatingLayer::Absorbing(&absorbing)];
            let explicit = SpectralOptics::conductor_stack_reflection(
                0.6, &incident, &bare, CoherentCoatingStack::try_new(&single).unwrap(),
            ).unwrap();
            let attached = SpectralOptics::conductor_reflection(
                0.6, &incident, &bare.with_absorbing_coating(absorbing),
            ).unwrap();
            for (a, b) in explicit.matrices.iter().zip(&attached.matrices) {
                for (a, b) in a.0.iter().flatten().zip(b.0.iter().flatten()) { close(*a, *b); }
            }
        }

        #[test]
        fn antireflection_stack_changes_flux_without_radiance_eta_scaling() {
            let air = CauchyIor::try_constant(1.0).unwrap();
            let glass = CauchyIor::try_constant(1.5).unwrap();
            let design_bin = SPECTRAL_BINS / 2;
            let film_index = 1.5_f64.sqrt();
            let film = CoherentDielectricCoating::try_new(
                CauchyIor::try_constant(film_index).unwrap(),
                wavelength_nm(design_bin) * 1e-9 / (4.0 * film_index), metal().source(),
            ).unwrap();
            let films = [film];
            let stack = DielectricCoatingStack::try_new(&films).unwrap();
            let original = source().with_unpolarized_optics();
            let reflected = original.through(&SpectralOptics::coated_dielectric_channel(
                1.0, air, glass, stack, false, DielectricEvent::Reflection,
            ).unwrap()).unwrap();
            let transmitted = original.through(&SpectralOptics::coated_dielectric_channel(
                1.0, air, glass, stack, false, DielectricEvent::Transmission,
            ).unwrap()).unwrap();
            let bare = original.through(&SpectralOptics::dielectric_channel(
                1.0, air, glass, DielectricEvent::Reflection,
            ).unwrap()).unwrap();
            assert!(reflected.relative_stokes()[design_bin].i < 1e-20);
            close(transmitted.relative_stokes()[design_bin].i, 1.0);
            assert!(reflected.intensity_source().flux_lm() < bare.intensity_source().flux_lm());
            close(reflected.intensity_source().flux_lm() + transmitted.intensity_source().flux_lm(), 1000.0);
            for (r, t) in reflected.relative_stokes().iter().zip(transmitted.relative_stokes()) {
                close(r.i + t.i, 1.0);
            }
        }

        #[test]
        fn absorbing_dielectric_stack_reverses_layers_without_restoring_lost_flux() {
            let (transparent, absorbing) = mixed_films(metal().source());
            let layers = [CoatingLayer::Dielectric(&transparent), CoatingLayer::Absorbing(&absorbing)];
            let reversed = [layers[1], layers[0]];
            let stack = DielectricCoatingStack::try_mixed(&layers).unwrap();
            let explicit_reverse = DielectricCoatingStack::try_mixed(&reversed).unwrap();
            let inside = CauchyIor::try_constant(1.5).unwrap();
            let outside = CauchyIor::try_constant(1.0).unwrap();
            let original = source().with_unpolarized_optics();
            let mut surviving_lumens = 0.0;
            for event in [DielectricEvent::Reflection, DielectricEvent::Transmission] {
                let optic = SpectralOptics::coated_dielectric_channel(
                    0.8, inside, outside, stack, true, event,
                ).unwrap();
                let manual = SpectralOptics::coated_dielectric_channel(
                    0.8, inside, outside, explicit_reverse, false, event,
                ).unwrap();
                assert_eq!(optic.matrices, manual.matrices);
                let output = original.through(&optic).unwrap();
                surviving_lumens += output.intensity_source().flux_lm();
                if event == DielectricEvent::Reflection {
                    let wrong_order = original.through(&SpectralOptics::coated_dielectric_channel(
                        0.8, inside, outside, stack, false, event,
                    ).unwrap()).unwrap();
                    assert!(output.relative_stokes().iter().zip(wrong_order.relative_stokes())
                        .any(|(a, b)| (a.i - b.i).abs() > 1e-4));
                }
            }
            assert!(surviving_lumens > 0.0 && surviving_lumens < 1000.0);
            close(original.intensity_source().flux_lm(), 1000.0);
        }

        #[test]
        fn reflection_channels_refuse_invalid_angles_and_missing_spectral_support() {
            let air = CauchyIor::try_constant(1.0).unwrap();
            let glass = CauchyIor::try_constant(1.5).unwrap();
            let bare = metal();
            let metal_stack = CoherentCoatingStack::try_new(&[]).unwrap();
            let empty = DielectricCoatingStack::try_new(&[]).unwrap();
            for cosine in [f64::NAN, f64::INFINITY, -0.1, 1.1] {
                assert!(SpectralOptics::conductor_reflection(cosine, air, &bare).is_err());
                assert!(SpectralOptics::conductor_stack_reflection(cosine, air, &bare, metal_stack).is_err());
                assert!(SpectralOptics::coated_dielectric_channel(
                    cosine, air, glass, empty, false, DielectricEvent::Reflection,
                ).is_err());
            }
            let narrow = SellmeierIor::try_new_in_band([1.25, 0.0, 0.0], [0.0; 3], (500.0, 700.0)).unwrap();
            assert!(SpectralOptics::conductor_reflection(0.5, narrow, &bare).is_err());
            assert!(SpectralOptics::coated_dielectric_channel(
                0.5, air, narrow, empty, false, DielectricEvent::Transmission,
            ).is_err());
            let narrow_substrate = ConductorOptics::try_new_in_band(core::array::from_fn(|i| {
                ConductorIorSample::try_new_in_band(500.0 + 25.0 * i as f64, 0.2, 3.0).unwrap()
            }), 25.0, bare.source()).unwrap();
            assert!(SpectralOptics::conductor_reflection(0.5, air, &narrow_substrate).is_err());
            assert!(SpectralOptics::conductor_stack_reflection(
                0.5, air, &narrow_substrate, metal_stack,
            ).is_err());
            let zero_film = CoherentDielectricCoating::try_new(narrow, 0.0, bare.source()).unwrap();
            let films = [CoatingLayer::Dielectric(&zero_film)];
            assert!(SpectralOptics::conductor_stack_reflection(
                0.5, air, &bare, CoherentCoatingStack::try_new(&films).unwrap(),
            ).is_err(), "zero thickness must not extend an optical source's support");
            assert!(SpectralOptics::coated_dielectric_channel(
                0.5, air, glass, DielectricCoatingStack::try_new_in_band(&[], (500.0, 700.0)).unwrap(),
                false, DielectricEvent::Reflection,
            ).is_err(), "an empty stack still owns its declared support");
        }

        #[test]
        fn invalid_filters_and_brightness_redefinition_are_refused() {
            assert!(SpectralOptics::linear_polarizer(f64::NAN).is_err());
            assert!(SpectralOptics::neutral([1.01; SPECTRAL_BINS]).is_err());
            assert!(SpectralOptics::linear_retarder(0.0, [f64::INFINITY; SPECTRAL_BINS]).is_err());
            assert!(PolarizedSource::try_new(source(), [Stokes::unpolarized(2.0); SPECTRAL_BINS]).is_err());
            let invalid = Stokes { i: 1.0, q: 2.0, u: 0.0, v: 0.0 };
            assert!(PolarizedSource::try_new(source(), [invalid; SPECTRAL_BINS]).is_err());
            let dark = source().with_unpolarized_optics()
                .through(&SpectralOptics::neutral([0.0; SPECTRAL_BINS]).unwrap()).unwrap()
                .through(&SpectralOptics::linear_polarizer(0.3).unwrap()).unwrap();
            assert_eq!(dark.intensity_source().flux_lm(), 0.0);
            assert!(dark.relative_stokes().iter().all(|s| *s == Stokes::unpolarized(0.0)));
        }
    }
}
