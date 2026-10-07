//! Exact affine-source lifting in the existing conservative flux space.
//!
//! Let f = sum_i f_i lambda_i and f_bar = sum_i f_i/4 on a tetrahedron.
//! The polynomial b = sum_i f_i lambda_i (x-v_i)/4 has zero normal trace
//! on EVERY face and div b = f-f_bar. Thus q_RT0+b equilibrates the actual
//! P1 source, not an averaged substitute. Shared-face and Robin fluxes remain
//! unchanged. All means, lifting coefficients and degree-four mass moments
//! are evaluated outward; no rounded projection is called the exact source.
use super::{Conductivity, Iv, TetError};
use super::geometry::{Cell, scale, sub};

#[derive(Clone, Copy)]
pub(super) enum Source<'a> {
    Constant(&'a [f64]),
    Affine(&'a [[f64; 4]]),
}
impl Source<'_> {
    pub(super) fn len(self) -> usize {
        match self { Self::Constant(f) => f.len(), Self::Affine(f) => f.len() }
    }
    pub(super) fn is_finite(self) -> bool {
        match self {
            Self::Constant(f) => f.iter().all(|v| v.is_finite()),
            Self::Affine(f) => f.iter().flatten().all(|v| v.is_finite()),
        }
    }
    fn constant(self, e: usize) -> Option<f64> {
        match self {
            Self::Constant(f) => Some(f[e]),
            Self::Affine(f) if f[e].iter().all(|&v| v == f[e][0]) => Some(f[e][0]),
            Self::Affine(_) => None,
        }
    }
    /// Encloses the EXACT integral divided by volume. Scale before adding,
    /// so a finite average need not overflow merely because its sum does.
    pub(super) fn mean(self, e: usize) -> Iv {
        if let Some(f) = self.constant(e) { return Iv::point(f); }
        let Self::Affine(f) = self else { unreachable!("constant source handled above") };
        f[e].iter().fold(Iv::zero(), |s, &v| s.add(Iv::point(v).scale_pos(0.25)))
    }
    /// Exact integral of source times a P1 dual. Multiplying separate averages
    /// would omit the covariance term and move the center of the goal bound.
    pub(super) fn load(self, e: usize, dual: &[Iv; 4], volume: Iv) -> Iv {
        if let Some(f) = self.constant(e) {
            return super::goal::linear_integral(dual, volume).mul(Iv::point(f));
        }
        let Self::Affine(f) = self else { unreachable!("constant source handled above") };
        super::goal::product_integral(&f[e].map(Iv::point), dual, volume)
    }
    /// Integrate |q_RT0+b+K grad v|^2_(K^-1) without quadrature approximation.
    /// The legacy constant-source operation sequence is preserved exactly.
    pub(super) fn defect_integral(self, e: usize, cell: &Cell,
        linear: &[[Iv; 3]; 4], conductivity: Conductivity<'_>) -> Result<Iv, TetError>
    {
        if self.constant(e).is_some() {
            return Ok(conductivity.defect_integral(e, linear, cell.volume));
        }
        let Self::Affine(f) = self else { unreachable!("constant source handled above") };
        let bubble = bubble_coefficients(f[e], &cell.points);
        let mut coefficients = [[Iv::zero(); 3]; 10];
        let mut powers = [[0usize; 4]; 10];
        for i in 0..4 { coefficients[i] = linear[i]; powers[i][i] = 1; }
        let mut term = 4;
        for i in 0..4 {
            for j in (i+1)..4 {
                coefficients[term] = bubble[term-4];
                powers[term][i] = 1; powers[term][j] = 1;
                term += 1;
            }
        }
        let mut integral = Iv::zero();
        for i in 0..10 {
            for j in i..10 {
                let weight = mass_moment(powers[i], powers[j], cell.volume);
                let product = conductivity.inverse_bilinear(e, coefficients[i], coefficients[j]);
                let contribution = product.mul(weight);
                integral = integral.add(if i == j { contribution } else { contribution.scale_pos(2.0) });
            }
        }
        if integral.is_unbounded() || integral.hi < 0.0 { return Err(TetError::Unbounded); }
        // The exact integrand is a squared SPD norm, even when individual
        // mixed coefficients are negative. Never clamp its upper endpoint.
        integral.lo = integral.lo.max(0.0);
        Ok(integral)
    }
}

/// Coefficients of lambda_i*lambda_j in increasing (i,j) order. Pairing terms
/// cancels a common source offset algebraically, before interval arithmetic.
fn bubble_coefficients(source: [f64; 4], points: &[[Iv; 3]; 4]) -> [[Iv; 3]; 6] {
    let mut result = [[Iv::zero(); 3]; 6];
    let mut term = 0;
    for i in 0..4 {
        for j in (i+1)..4 {
            let difference = Iv::point(source[i]).scale_pos(0.25)
                .sub(Iv::point(source[j]).scale_pos(0.25));
            result[term] = scale(sub(points[j], points[i]), difference);
            term += 1;
        }
    }
    result
}

/// integral(lambda^(a+b)) = volume*3!*product((a_i+b_i)!)/(3+|a+b|)!.
/// Our linear/quadratic terms require only degrees 2, 3 and 4. Every integer
/// factor below is represented exactly; the division itself remains outward.
fn mass_moment(a: [usize; 4], b: [usize; 4], volume: Iv) -> Iv {
    let mut degree = 0;
    let mut numerator = 1usize;
    const FACTORIAL: [usize; 5] = [1, 1, 2, 6, 24];
    for i in 0..4 { let power = a[i]+b[i]; degree += power; numerator *= FACTORIAL[power]; }
    let denominator = match degree { 2 => 20.0, 3 => 120.0, 4 => 840.0,
        _ => unreachable!("linear/quadratic polynomial product") };
    volume.scale_pos(numerator as f64).div_pos(Iv::point(denominator))
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::geometry::dot;
    fn contains(x: Iv, truth: f64) { assert!(x.lo <= truth && truth <= x.hi, "{x:?} misses {truth}"); }
    #[test]
    fn exact_source_means_and_affine_dual_products_are_not_products_of_averages() {
        let f = [[0.0, 4.0, 0.0, 0.0]];
        let source = Source::Affine(&f);
        contains(source.mean(0), 1.0);
        let z = [0.0, 1.0, 0.0, 0.0].map(Iv::point);
        // Unit reference tet: int 4*x*x = 1/15; int f * mean(z) = 1/24 is wrong.
        let load = source.load(0, &z, Iv::point(1.0).div_pos(Iv::point(6.0)));
        contains(load, 1.0/15.0);
        assert!(load.lo > 1.0/24.0);
        let huge = [[f64::MAX, f64::MAX, -f64::MAX, -f64::MAX]];
        let mean = Source::Affine(&huge).mean(0);
        assert!(!mean.is_unbounded()); contains(mean, 0.0);
    }
    #[test]
    fn polynomial_moments_include_cross_term_multiplicities() {
        contains(mass_moment([1,0,0,0], [1,0,0,0], Iv::point(1.0)), 0.1);
        contains(mass_moment([1,0,0,0], [0,1,0,0], Iv::point(1.0)), 0.05);
        contains(mass_moment([1,1,0,0], [1,1,0,0], Iv::point(1.0)), 1.0/210.0);
        contains(mass_moment([1,1,0,0], [0,0,1,1], Iv::point(1.0)), 1.0/840.0);
        contains(mass_moment([1,0,0,0], [1,1,0,0], Iv::point(1.0)), 1.0/60.0);
    }
    #[test]
    fn bubble_has_zero_trace_on_all_faces_and_correct_divergence() {
        let vertices: [[f64;3];4] = [[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]];
        let points = vertices.map(|v| v.map(Iv::point));
        let f = [3.0, -2.0, 5.0, 10.0];
        let coefficients = bubble_coefficients(f, &points);
        let gradients = [[-1.,-1.,-1.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]];
        for lambda in [[0.25;4],[0.,0.2,0.3,0.5],[0.2,0.,0.3,0.5],
            [0.2,0.3,0.,0.5],[0.2,0.3,0.5,0.]] {
            let mut flux = [Iv::zero();3]; let mut div = Iv::zero(); let mut term = 0;
            for i in 0..4 { for j in (i+1)..4 {
                for d in 0..3 { flux[d] = flux[d].add(coefficients[term][d]
                    .mul(Iv::point(lambda[i])).mul(Iv::point(lambda[j]))); }
                let grad = std::array::from_fn(|d| Iv::point(gradients[i][d]).mul(Iv::point(lambda[j]))
                    .add(Iv::point(gradients[j][d]).mul(Iv::point(lambda[i]))));
                div = div.add(dot(coefficients[term], grad)); term += 1;
            }}
            let expected: f64 = f.iter().zip(lambda).map(|(f,l)| f*l).sum::<f64>()-4.0;
            contains(div, expected);
            for i in 0..4 { if lambda[i] == 0.0 { contains(dot(flux, gradients[i].map(Iv::point)), 0.0); } }
        }
    }
}
