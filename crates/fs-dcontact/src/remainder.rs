//! Cancellation-stable evaluation of the original Gonzalez contact numerator.
//! A small cubic remainder must not be obtained by subtracting two O(H)
//! endpoint energies. HC's effort derivative amplifies that noise in Newton.
//! No force, potential, loss law, solver request or admission floor is changed.
use super::{ContactStorage, det};
use fs_math::eft::two_sum;

/// `(b^beta-a^beta)/beta - midpoint^(beta-1)*(b-a)` for active a,b.
/// With c=(a+b)/2 and z=(b-a)/(2c), its centered remainder is
/// `2*c^beta/beta * sum_{odd j>=3} binom(beta,j)*z^j`.
/// For |beta*z|<=1/2 and beta>=2, |z|<=1/4 and EVERY successive
/// term contracts by <=1/4. Terms through j=33 leave a tail <2^-63
/// of the leading cubic term; this is roundoff-level evaluation of the
/// power law, not a truncated constitutive model or fitted coefficient.
fn power_remainder(a: f64, b: f64, midpoint: f64, alpha: f64, coefficient: f64) -> f64 {
    let beta = alpha + 1.0;
    let increment = b - a;
    if a > 0.0 && b > 0.0 {
        let center = f64::midpoint(a, b);
        let z = (0.5 * increment) / center;
        if (beta * z).abs() <= 0.5 {
            let mut term = beta * z;
            let mut odd_sum = 0.0;
            for j in 2..=33 {
                term *= ((beta - f64::from(j - 1)) / f64::from(j)) * z;
                if j % 2 == 1 {
                    odd_sum += term;
                }
            }
            // Scale in the original potential's order before the small cubic
            // term, avoiding a new underflow before a large physical K.
            let centered = coefficient / beta * det::pow(center, beta) * (2.0 * odd_sum);
            // An affine collocation midpoint need not round to exactly c.
            // Correct that discrepancy at force*increment scale, not H scale.
            return centered
                + (coefficient * det::pow(center, beta - 1.0)
                    - coefficient * det::pow(midpoint.max(0.0), alpha))
                    * increment;
        }
    }
    // Entry/separation and large increments do not subtract nearly equal H.
    coefficient / beta * (det::pow(b.max(0.0), beta) - det::pow(a.max(0.0), beta))
        - coefficient * det::pow(midpoint.max(0.0), alpha) * increment
}

impl ContactStorage {
    pub(super) fn stable_midpoint_remainder(
        &self,
        a: &[f64],
        b: &[f64],
        midpoint: &[f64],
    ) -> Option<f64> {
        let n = a.len();
        if b.len() != n
            || midpoint.len() != n
            || self.n_modes.checked_mul(2).is_none_or(|prefix| n < prefix)
        {
            return None;
        }
        // A declared exact quadratic has zero mathematical midpoint remainder.
        // Other inner laws must provide their own complete stable contribution.
        // The TOTAL contact storage remains nonlinear and is never advertised
        // as exact quadratic to the Newton-floor/linear-evolution machinery.
        let inner = self
            .inner
            .midpoint_energy_remainder(a, b, midpoint)
            .or_else(|| {
                // Read the CURRENT borrowed declaration; a cold cached flag
                // could outlive a material/history update of an inner law.
                self.inner
                    .exact_quadratic_matrix()
                    .filter(|matrix| n.checked_mul(n) == Some(matrix.len()))
                    .map(|_| 0.0)
            })?;
        let (mut sum, mut low) = (inner, 0.0);
        for ob in &self.obstacles {
            for i in 0..ob.n_points {
                let pa = ob.penetration_at(self.n_modes, a, i);
                let pb = ob.penetration_at(self.n_modes, b, i);
                let pm = ob.penetration_at(self.n_modes, midpoint, i);
                let mut row_increment = 0.0;
                for k in 0..self.n_modes {
                    row_increment += ob.collocation[i * self.n_modes + k] * (b[2 * k] - a[2 * k]);
                }
                // Account for rounding in the affine penetration map as well as
                // the nonlinear power. The gradient contracts Phi with b-a.
                let coefficient = ob.weights[i] * ob.stiffness;
                let term = power_remainder(pa, pb, pm, ob.alpha, coefficient)
                    + coefficient * det::pow(pm.max(0.0), ob.alpha) * ((pb - pa) - row_increment);
                let (next, roundoff) = two_sum(sum, term);
                sum = next;
                low += roundoff;
            }
        }
        Some(sum + low)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_math::dd::Dd;
    use fs_phs::{QuadraticStorage, Storage};

    #[test]
    fn hertz_midpoint_remainder_matches_independent_double_double_at_three_scales() {
        // Independent exact half-integer powers via the existing ~106-bit
        // sqrt/multiply owner, not the centered series being tested.
        for h in [1.0_f64 / 256.0, 1.0 / 4096.0, 1.0 / 65536.0] {
            for sign in [-1.0, 1.0] {
                let (a, b) = (1.0, 1.0 + sign * h);
                let midpoint = f64::midpoint(a, b);
                let (ad, bd, md) = (Dd::from_f64(a), Dd::from_f64(b), Dd::from_f64(midpoint));
                let reference = ((bd * bd * bd.sqrt() - ad * ad * ad.sqrt()) / Dd::from_f64(2.5)
                    - md * md.sqrt() * (bd - ad))
                    .to_f64();
                let actual = power_remainder(a, b, midpoint, 1.5, 1.0);
                assert!(
                    (actual - reference).abs() <= 1e-12 * reference.abs(),
                    "h={h}, sign={sign}: {actual:e} vs {reference:e}"
                );
            }
        }
    }

    #[test]
    fn unilateral_entry_and_separation_keep_the_original_potential() {
        let h = 0.03125_f64;
        for (a, b) in [(-h, h), (h, -h), (-2.0 * h, -h)] {
            let (ad, bd) = (Dd::from_f64(a.max(0.0)), Dd::from_f64(b.max(0.0)));
            let reference =
                ((bd * bd * bd.sqrt() - ad * ad * ad.sqrt()) / Dd::from_f64(2.5)).to_f64();
            let actual = power_remainder(a, b, f64::midpoint(a, b), 1.5, 1.0);
            assert!((actual - reference).abs() <= 1e-13 * reference.abs());
        }
        // Force uses the ACTUAL affine midpoint, not an inferred centered one.
        let exact = power_remainder(1.0, 1.0001, 1.00005, 1.5, 1.0);
        let shifted = power_remainder(1.0, 1.0001, 1.01, 1.5, 1.0);
        let difference = (det::pow(1.00005, 1.5) - det::pow(1.01, 1.5)) * (1.0001 - 1.0);
        assert!((shifted - exact - difference).abs() < 1e-15 * difference.abs());
    }

    #[test]
    fn distributed_remainder_retains_signed_rows_weights_gaps_and_internal_tail() {
        let mut q = vec![0.0; 25];
        for i in 0..5 {
            q[5 * i + i] = (i + 1) as f64;
        }
        let obstacle = super::super::Obstacle::new(
            vec![1.0, -0.5, -1.0, 2.0, 0.0, 1.0],
            3,
            2,
            vec![0.001, 0.0, 0.01],
            vec![0.3, 0.4, 0.8],
            2e5,
            1.5,
            "synthetic signed distributed remainder test".into(),
        )
        .unwrap();
        let storage = ContactStorage::new(
            Box::new(QuadraticStorage::new(q, 5).unwrap()),
            2,
            vec![obstacle],
        )
        .unwrap();
        let a = [0.005, 0.1, 0.003, -0.2, 7.0];
        let b = [0.004, 0.11, 0.0015, -0.19, 7.0];
        let midpoint = std::array::from_fn::<_, 5, _>(|i| f64::midpoint(a[i], b[i]));
        let actual = storage
            .midpoint_energy_remainder(&a, &b, &midpoint)
            .unwrap();
        let mut reference = Dd::ZERO;
        for ob in &storage.obstacles {
            for i in 0..ob.n_points {
                let pa = Dd::from_f64(ob.penetration_at(2, &a, i).max(0.0));
                let pb = Dd::from_f64(ob.penetration_at(2, &b, i).max(0.0));
                let pm = Dd::from_f64(ob.penetration_at(2, &midpoint, i).max(0.0));
                let mut increment = Dd::ZERO;
                for k in 0..2 {
                    increment = increment
                        + Dd::from_f64(ob.collocation[2 * i + k])
                            * (Dd::from_f64(b[2 * k]) - Dd::from_f64(a[2 * k]));
                }
                reference = reference
                    + Dd::from_f64(ob.weights[i] * ob.stiffness)
                        * ((pb * pb * pb.sqrt() - pa * pa * pa.sqrt()) / Dd::from_f64(2.5)
                            - pm * pm.sqrt() * increment);
            }
        }
        let reference = reference.to_f64();
        assert!((actual - reference).abs() <= 1e-12 * reference.abs());
    }
}
