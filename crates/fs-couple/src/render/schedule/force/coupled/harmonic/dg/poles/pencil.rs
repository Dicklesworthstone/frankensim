//! Real polynomial entries; fixed row equilibration never changes coordinates.
use fs_la::nonlinear_eig::beyn::{BeynError, BeynOptions, BeynReport, Ellipse, beyn_with_cancel};
use fs_la::nonlinear_eig::{NonlinearEigError, PoleEvaluation};
use fs_math::{c64::C64, dd::Dd, det};

pub(super) struct Pencil {
    pub constant: Vec<f64>,
    pub linear: Vec<f64>,
    pub quadratic_diagonal: Vec<f64>,
    pub dimension: usize,
    rows: Vec<f64>,
    norms: [f64; 3],
    reference: f64,
}

fn invalid(what: &'static str) -> NonlinearEigError {
    NonlinearEigError::BadInput { what }
}
fn poll(cancel: &mut (impl FnMut() -> bool + ?Sized)) -> Result<(), NonlinearEigError> {
    if cancel() {
        Err(NonlinearEigError::Cancelled)
    } else {
        Ok(())
    }
}
fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}

impl Pencil {
    /// The physical composer admits both allocations before calling this.
    pub fn new(dimension: usize, reference: f64) -> Self {
        Self {
            constant: vec![0.0; dimension * dimension],
            linear: vec![0.0; dimension * dimension],
            quadratic_diagonal: vec![0.0; dimension],
            dimension,
            rows: vec![0.0; dimension],
            norms: [0.0; 3],
            reference,
        }
    }

    pub fn condition(
        &mut self,
        cancel: &mut impl FnMut() -> bool,
    ) -> Result<(), NonlinearEigError> {
        let n = self.dimension;
        if !self.reference.is_finite() || self.reference <= 0.0 {
            return Err(invalid(
                "joint pole reference rate must be positive and finite",
            ));
        }
        self.norms = [0.0; 3];
        for i in 0..n {
            poll(cancel)?;
            if !self.quadratic_diagonal[i].is_finite() {
                return Err(invalid("joint original inertia is not finite"));
            }
            let mut scale = self.reference * self.reference * self.quadratic_diagonal[i].abs();
            for j in 0..n {
                let a = self.constant[i * n + j];
                let b = self.linear[i * n + j];
                if !a.is_finite() || !b.is_finite() {
                    return Err(invalid("joint original pencil is not finite"));
                }
                scale = scale.max(a.abs()).max(self.reference * b.abs());
            }
            if !scale.is_finite() || scale <= 0.0 {
                return Err(invalid("joint pole row has no representable nonzero scale"));
            }
            self.rows[i] = scale;
            for j in 0..n {
                self.norms[0] = det::hypot(self.norms[0], self.constant[i * n + j] / scale);
                self.norms[1] = det::hypot(self.norms[1], self.linear[i * n + j] / scale);
            }
            self.norms[2] = det::hypot(self.norms[2], self.quadratic_diagonal[i] / scale);
        }
        if self.norms.iter().any(|v| !v.is_finite()) {
            return Err(invalid("joint pole coefficient envelope overflows"));
        }
        Ok(())
    }

    pub fn evaluate(
        &self,
        s: C64,
        cancel: &mut dyn FnMut() -> bool,
    ) -> Result<PoleEvaluation, NonlinearEigError> {
        poll(cancel)?;
        let n = self.dimension;
        let rate = s.abs().max(self.reference);
        let residual_scale = self.norms[0] + rate * self.norms[1] + rate * rate * self.norms[2];
        let square = s * s;
        if !finite(s) || !finite(square) || !residual_scale.is_finite() || residual_scale <= 0.0 {
            return Err(NonlinearEigError::FunctionRefused { at: s });
        }
        let mut matrix = vec![C64::ZERO; n * n];
        let mut derivative = vec![C64::ZERO; n * n];
        for i in 0..n {
            poll(cancel)?;
            let scale = self.rows[i];
            for j in 0..n {
                let a = self.constant[i * n + j] / scale;
                let b = self.linear[i * n + j] / scale;
                matrix[i * n + j] = C64::from_re(a) + s.scale(b);
                derivative[i * n + j] = C64::from_re(b);
            }
            let m = self.quadratic_diagonal[i] / scale;
            matrix[i * n + i] = matrix[i * n + i] + square.scale(m);
            derivative[i * n + i] = derivative[i * n + i] + s.scale(2.0 * m);
        }
        Ok(PoleEvaluation {
            matrix,
            derivative,
            residual_scale,
        })
    }

    /// Apply the SAME fixed row equilibration to a complete physical RHS.
    /// Spectral s and state coordinates remain unchanged by this operation.
    pub fn scale_rhs(
        &self,
        rhs: &mut [C64],
        cancel: &mut impl FnMut() -> bool,
    ) -> Result<(), NonlinearEigError> {
        if rhs.len() != self.dimension {
            return Err(invalid("joint pole forcing must retain every original row"));
        }
        for (value, &row) in rhs.iter_mut().zip(&self.rows) {
            poll(cancel)?;
            let scaled = C64::new(value.re / row, value.im / row);
            if !finite(scaled)
                || (value.re != 0.0 && scaled.re == 0.0)
                || (value.im != 0.0 && scaled.im == 0.0)
            {
                return Err(invalid("joint pole forcing row scale is not representable"));
            }
            *value = scaled;
        }
        Ok(())
    }

    /// The same entire pencil, including all memory and algebraic rows. Reserve
    /// callback work and complete physical result vectors before the numerical
    /// owner admits its LU/moment/SVD/QR/refinement resources. Cached real
    /// coefficients were already admitted at preparation and are excluded here.
    pub fn discover(
        &self,
        region: Ellipse,
        options: &BeynOptions,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<BeynReport, BeynError> {
        poll(&mut cancelled)?;
        let n = self.dimension as u128;
        let modes = options.probe_columns as u128;
        if modes == 0 || modes > n {
            return Err(BeynError::Invalid {
                what: "joint spectrum probes must lie between one and the complete dimension",
            });
        }
        let bytes = 128 * n * modes + 256 * (n + modes);
        let callbacks = options.nodes as u128 + modes * options.refinement.max_iterations as u128;
        let work = (256 * modes + 16 * callbacks) * n * n;
        if bytes > options.max_dense_bytes as u128 || work > options.max_work as u128 {
            return Err(BeynError::Invalid {
                what: "joint spectrum callbacks and observations exceed query resources",
            });
        }
        let mut numerical = *options;
        numerical.max_dense_bytes -= bytes as usize;
        numerical.max_work -= work as u64;
        beyn_with_cancel(
            self.dimension,
            region,
            &numerical,
            |s, cancel| self.evaluate(s, cancel),
            &mut cancelled,
        )
    }

    /// Check uncombined ORIGINAL coefficient rows in double-double arithmetic.
    /// Each row's coefficient envelope multiplies the full vector infinity norm.
    /// Thus a nodal/dark coordinate may contain roundoff without an O(1) defect,
    /// while fixed row equilibration preserves small physical equation rows.
    /// Stiffness, rate and inertia remain separate in the envelope at resonance.
    pub fn relative_residual(
        &self,
        s: C64,
        state: &[C64],
        cancel: &mut impl FnMut() -> bool,
    ) -> Result<f64, NonlinearEigError> {
        let n = self.dimension;
        if state.len() != n || state.iter().any(|&v| !finite(v)) || !finite(s) {
            return Err(invalid(
                "joint pole state must have every finite original coordinate",
            ));
        }
        let state_norm = state.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
        if state_norm == 0.0 || !state_norm.is_finite() {
            return Err(invalid(
                "joint pole state must have a representable nonzero norm",
            ));
        }
        let re = Dd::from_f64(s.re);
        let im = Dd::from_f64(s.im);
        let sr = re * re - im * im;
        let si = (re * im) * Dd::from_f64(2.0);
        let mut relative = 0.0_f64;
        for i in 0..n {
            poll(cancel)?;
            let scale = Dd::from_f64(self.rows[i]);
            let mut real = Dd::ZERO;
            let mut imag = Dd::ZERO;
            let mut denominator = 0.0;
            for (j, &x) in state.iter().enumerate() {
                let a = Dd::from_f64(self.constant[i * n + j]) / scale;
                let b = Dd::from_f64(self.linear[i * n + j]) / scale;
                let mut cr = a + re * b;
                let mut ci = im * b;
                let mut magnitude = a.abs().to_f64() + s.abs() * b.abs().to_f64();
                if i == j {
                    let m = Dd::from_f64(self.quadratic_diagonal[i]) / scale;
                    cr = cr + sr * m;
                    ci = ci + si * m;
                    magnitude += s.abs() * s.abs() * m.abs().to_f64();
                }
                let xr = Dd::from_f64(x.re);
                let xi = Dd::from_f64(x.im);
                real = real + cr * xr - ci * xi;
                imag = imag + cr * xi + ci * xr;
                denominator += magnitude;
            }
            denominator *= state_norm;
            let residual = det::hypot(real.to_f64(), imag.to_f64());
            let row = if denominator == 0.0 {
                if residual == 0.0 { 0.0 } else { f64::INFINITY }
            } else {
                residual / denominator
            };
            if !row.is_finite() || !denominator.is_finite() {
                return Err(invalid("joint original pole residual is not representable"));
            }
            relative = relative.max(row);
        }
        Ok(relative)
    }
}

#[cfg(test)]
#[path = "pencil/spectrum_tests.rs"]
mod spectrum_tests;

#[cfg(test)]
mod tests {
    use super::*;
    use fs_la::nonlinear_eig::{PoleOptions, newton_pole_with_cancel};

    fn options() -> PoleOptions {
        PoleOptions {
            spectral_scale: 4.0,
            step_tolerance: 1e-12,
            residual_tolerance: 1e-12,
            max_iterations: 24,
            max_dimension: 8,
            max_dense_bytes: 1 << 20,
            max_work: 1_000_000,
        }
    }
    fn solve(p: &Pencil, s: C64, v: &[C64]) -> fs_la::nonlinear_eig::NonlinearPole {
        let result = newton_pole_with_cancel(
            p.dimension,
            s,
            v,
            &options(),
            |s, c| p.evaluate(s, c),
            &mut || false,
        )
        .unwrap();
        assert!(
            p.relative_residual(result.value, &result.vector, &mut || false)
                .unwrap()
                < 2e-12
        );
        result
    }

    #[test]
    fn original_mass_and_drag_resolve_both_damped_poles_under_row_rescaling() {
        // G1: 2*q''+4*q'+10*q=0 has s=-1+-2i. Coefficient units can
        // differ by 160 orders without changing the original state or roots.
        for scale in [1e-80, 1.0, 1e80] {
            let mut p = Pencil::new(1, 4.0);
            p.constant[0] = 10.0 * scale;
            p.linear[0] = 4.0 * scale;
            p.quadratic_diagonal[0] = 2.0 * scale;
            p.condition(&mut || false).unwrap();
            for sign in [-1.0, 1.0] {
                let r = solve(&p, C64::new(-0.8, sign * 1.9), &[C64::ONE]);
                assert!((r.value - C64::new(-1.0, sign * 2.0)).abs() < 2e-12);
            }
        }
    }

    #[test]
    fn maxwell_coordinates_retain_the_coupled_cubic_and_a_dark_real_relaxation() {
        // G1: q''+6*q'+6*q+(q-z)=0; z'=2*(q-z).
        // Eliminating z gives (s+1)(s+3)(s+4), but the separate dark
        // strain w'=-5w must survive without any displacement observation.
        let mut p = Pencil::new(3, 4.0);
        p.constant = vec![7.0, -1.0, 0.0, -2.0, 2.0, 0.0, 0.0, 0.0, 5.0];
        p.linear = vec![6.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        p.quadratic_diagonal[0] = 1.0;
        p.condition(&mut || false).unwrap();
        for s in [-1.0, -3.0, -4.0] {
            let r = solve(
                &p,
                C64::from_re(s * 0.97),
                &[C64::ONE, C64::from_re(2.0 / (s + 2.0)), C64::ZERO],
            );
            assert!((r.value - C64::from_re(s)).abs() < 1e-11);
        }
        let dark = solve(&p, C64::from_re(-4.8), &[C64::ZERO, C64::ZERO, C64::ONE]);
        assert_eq!(dark.value, C64::from_re(-5.0));
        assert_eq!(dark.vector[0], C64::ZERO);
        assert_eq!(dark.vector[1], C64::ZERO);
        let contaminated = solve(
            &p,
            C64::from_re(-4.8),
            &[C64::new(1e-13, -2e-13), C64::from_re(7e-13), C64::ONE],
        );
        assert!((contaminated.value - C64::from_re(-5.0)).abs() < 1e-12);
        assert!(contaminated.vector[0].abs() < 1e-12);
        assert!(contaminated.vector[1].abs() < 1e-12);
    }

    #[test]
    fn finite_gap_inertia_and_pressure_constraint_give_the_physical_added_mass() {
        // G1: q''+q+dp=0; v'+2v=dp; v=q'. The full singular
        // descriptor has the physical roots of 2*s^2+2*s+1, including
        // the gap velocity and its pressure multiplier at each root.
        let mut p = Pencil::new(3, 1.0);
        p.constant = vec![1.0, 0.0, 1.0, 0.0, 2.0, -1.0, 0.0, 1.0, 0.0];
        p.linear = vec![0.0, 0.0, 0.0, 0.0, 1.0, 0.0, -1.0, 0.0, 0.0];
        p.quadratic_diagonal[0] = 1.0;
        p.condition(&mut || false).unwrap();
        for sign in [-1.0, 1.0] {
            let s = C64::new(-0.5, sign * 0.5);
            let r = solve(
                &p,
                s.scale(1.1),
                &[C64::ONE, s, (s + C64::from_re(2.0)) * s],
            );
            assert!((r.value - s).abs() < 2e-12);
            assert!((r.vector[1] - r.value * r.vector[0]).abs() < 2e-12);
        }
    }

    #[test]
    fn original_row_residual_and_cancellation_do_not_depend_on_global_scales() {
        let mut p = Pencil::new(2, 1.0);
        p.constant = vec![1e80, 0.0, 0.0, 1e-80];
        p.linear = vec![1e80, 0.0, 0.0, 2e-80];
        p.condition(&mut || false).unwrap();
        // The first row is exact at -1; the other equation is wrong by 1/3,
        // despite being 160 decades below it in an unscaled global norm.
        let r = p
            .relative_residual(C64::from_re(-1.0), &[C64::ONE; 2], &mut || false)
            .unwrap();
        assert!((r - 1.0 / 3.0).abs() < 1e-15);
        assert!(matches!(
            p.evaluate(C64::ONE, &mut || true),
            Err(NonlinearEigError::Cancelled)
        ));
        assert!(matches!(
            p.relative_residual(C64::ONE, &[C64::ONE; 2], &mut || true),
            Err(NonlinearEigError::Cancelled)
        ));
        let mut zero = Pencil::new(1, 1.0);
        assert!(zero.condition(&mut || false).is_err());
    }
}
