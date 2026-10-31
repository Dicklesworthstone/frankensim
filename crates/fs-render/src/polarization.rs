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
}
