//! Stable work of the original condensed, positive facet stretching energy.
use super::MembraneReduction;
use fs_math::dd::Dd;

impl MembraneReduction {
    /// Original complete mass-normalized tension/bending matrix, row major.
    /// The admitted off-diagonal projection terms remain part of the energy;
    /// replacing them with squared reported frequencies changes the model.
    #[must_use]
    pub fn linear_stiffness_matrix(&self) -> &[f64] {
        &self.linear
    }

    /// `U(b)-U(a)-grad U(mid).(b-a)` for the SAME condensed stretching energy.
    ///
    /// Quadratic strain increments evaluate this quartic work directly, without
    /// subtracting large endpoint energies. Every facet metric and mixed modal
    /// pair survives, including rounding of the represented midpoint. The
    /// original linear tension/bending is excluded so its time owner can retain
    /// that complete quadratic block without subtracting it again.
    ///
    /// Work is O(facets*modal_pairs), with fixed three-component scratch and no
    /// heap allocation or new fourth-order tensor. Two-word sums/products retain
    /// tiny nonlinear work near coincidence. This supplies no integrator, material
    /// history, extra energy or change to the original force.
    ///
    /// Returns None for incomplete/nonfinite states or unrepresentable work.
    #[must_use]
    pub fn stretching_midpoint_energy_remainder(
        &self,
        a: &[f64],
        b: &[f64],
        mid: &[f64],
    ) -> Option<f64> {
        let n = self.mode_count();
        if a.len() != n
            || b.len() != n
            || mid.len() != n
            || a.iter().chain(b).chain(mid).any(|v| !v.is_finite())
        {
            return None;
        }
        let dd = Dd::from_f64;
        let mut work = dd(0.0);
        for facet in &self.facets {
            let mut strain = [dd(0.0); 3];
            let mut change = [dd(0.0); 3];
            let mut cross = [dd(0.0); 3];
            let mut offset = [dd(0.0); 3];
            for (&(i, j), coefficient) in self.pairs.iter().zip(&facet.strains) {
                let (mi, mj) = (dd(mid[i]), dd(mid[j]));
                let (ai, aj) = (dd(a[i]) - mi, dd(a[j]) - mj);
                let (bi, bj) = (dd(b[i]) - mi, dd(b[j]) - mj);
                let (di, dj) = (dd(b[i]) - dd(a[i]), dd(b[j]) - dd(a[j]));
                let (si, sj) = (ai + bi, aj + bj);
                let linear_change = mi * dj + mj * di;
                let curvature_change = dd(0.5) * (si * dj + sj * di);
                let mean_offset = mi * sj + mj * si + ai * aj + bi * bj;
                for channel in 0..3 {
                    let weight = dd(coefficient[channel]);
                    strain[channel] = strain[channel] + weight * mi * mj;
                    change[channel] = change[channel] + weight * (linear_change + curvature_change);
                    cross[channel] = cross[channel] + weight * curvature_change;
                    offset[channel] = offset[channel] + weight * mean_offset;
                }
            }
            // For U=epsilon^T A epsilon/2:
            // remainder = offset^T A change/2 + epsilon(mid)^T A cross.
            // cross vanishes at an exactly represented endpoint midpoint.
            for i in 0..3 {
                for j in 0..3 {
                    work = work
                        + dd(facet.a[3 * i + j])
                            * (dd(0.5) * offset[i] * change[j] + strain[i] * cross[j]);
                }
            }
        }
        let result = work.to_f64();
        result.is_finite().then_some(result)
    }
}

#[cfg(test)]
mod tests {
    use super::super::Facet;
    use super::*;

    fn specimen() -> MembraneReduction {
        MembraneReduction {
            omegas: vec![2.0, 3.0],
            linear: vec![4.0, 0.25, 0.25, 9.0],
            pairs: vec![(0, 0), (0, 1), (1, 1)],
            solve_residual: 0.0,
            facets: vec![
                Facet {
                    a: [4.0, 1.0, 0.25, 1.0, 3.0, -0.125, 0.25, -0.125, 2.0],
                    gradients: vec![[0.0; 2]; 2],
                    strains: vec![[1.0, 0.125, 0.0], [0.25, 0.375, 0.5], [-0.25, 0.5, 0.125]],
                },
                Facet {
                    a: [2.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 3.0],
                    gradients: vec![[0.0; 2]; 2],
                    strains: vec![[0.0, 0.375, 0.25], [0.125, -0.5, 0.25], [0.75, 0.125, 0.5]],
                },
            ],
        }
    }

    // Independent full polynomial and derivative in two-word arithmetic.
    // Only this test subtracts endpoint energies; production uses increments.
    fn reference(law: &MembraneReduction, a: &[f64], b: &[f64], mid: &[f64]) -> f64 {
        let dd = Dd::from_f64;
        let energy = |q: &[f64]| {
            let mut result = dd(0.0);
            for facet in &law.facets {
                let mut strain = [dd(0.0); 3];
                for (&(i, j), coefficient) in law.pairs.iter().zip(&facet.strains) {
                    for channel in 0..3 {
                        strain[channel] =
                            strain[channel] + dd(coefficient[channel]) * dd(q[i]) * dd(q[j]);
                    }
                }
                for i in 0..3 {
                    for j in 0..3 {
                        result = result + dd(0.5 * facet.a[3 * i + j]) * strain[i] * strain[j];
                    }
                }
            }
            result
        };
        let mut gradient = vec![dd(0.0); mid.len()];
        for facet in &law.facets {
            let mut strain = [dd(0.0); 3];
            for (&(i, j), coefficient) in law.pairs.iter().zip(&facet.strains) {
                for channel in 0..3 {
                    strain[channel] =
                        strain[channel] + dd(coefficient[channel]) * dd(mid[i]) * dd(mid[j]);
                }
            }
            let stress: [Dd; 3] = core::array::from_fn(|i| {
                (0..3).fold(dd(0.0), |s, j| s + dd(facet.a[3 * i + j]) * strain[j])
            });
            for (&(i, j), coefficient) in law.pairs.iter().zip(&facet.strains) {
                let force = (0..3).fold(dd(0.0), |s, c| s + dd(coefficient[c]) * stress[c]);
                gradient[i] = gradient[i] + dd(mid[j]) * force;
                gradient[j] = gradient[j] + dd(mid[i]) * force;
            }
        }
        let work = gradient
            .iter()
            .zip(a)
            .zip(b)
            .fold(dd(0.0), |s, ((g, a), b)| s + *g * (dd(*b) - dd(*a)));
        (energy(b) - energy(a) - work).to_f64()
    }

    #[test]
    fn complete_film_work_retains_mixed_strain_and_represented_midpoint() {
        let law = specimen();
        assert_eq!(law.linear_stiffness_matrix(), &[4.0, 0.25, 0.25, 9.0]);
        for (a, b) in [
            ([0.125, -0.375], [0.5, 0.25]),
            ([-0.7, 0.3], [0.2, -0.4]),
            ([0.0, 0.0], [0.3, -0.1]),
        ] {
            for shift in [0.0, 1e-7] {
                let mid = [
                    f64::midpoint(a[0], b[0]) + shift,
                    f64::midpoint(a[1], b[1]) - shift,
                ];
                let expected = reference(&law, &a, &b, &mid);
                let actual = law
                    .stretching_midpoint_energy_remainder(&a, &b, &mid)
                    .unwrap();
                assert!((actual - expected).abs() < 2e-14 * expected.abs().max(1e-5));
                let reverse = law
                    .stretching_midpoint_energy_remainder(&b, &a, &mid)
                    .unwrap();
                assert_eq!(actual, -reverse);
                let mut gradient = [0.0; 2];
                law.gradient(&mid, &mut gradient);
                let linear = law.linear_stiffness_matrix();
                let native_work: f64 = (0..2)
                    .map(|i| {
                        let linear_force: f64 = (0..2).map(|j| linear[2 * i + j] * mid[j]).sum();
                        (gradient[i] - linear_force) * (b[i] - a[i])
                    })
                    .sum();
                let native_change = law.stretching_energy(&b) - law.stretching_energy(&a);
                assert!(
                    (native_change - native_work - actual).abs()
                        < 2e-13 * (1.0 + native_change.abs())
                );
            }
        }
    }

    #[test]
    fn near_coincident_film_work_keeps_tiny_geometric_remainder() {
        let law = specimen();
        let a = [0.31, -0.17];
        for exponent in [4, 12, 20, 24] {
            let step = 2.0_f64.powi(-exponent);
            for sign in [-1.0, 1.0] {
                let b = [a[0] + sign * step, a[1] - 0.375 * sign * step];
                let mid = core::array::from_fn::<_, 2, _>(|i| f64::midpoint(a[i], b[i]));
                let expected = reference(&law, &a, &b, &mid);
                let actual = law
                    .stretching_midpoint_energy_remainder(&a, &b, &mid)
                    .unwrap();
                assert!(actual != 0.0);
                assert!(
                    (actual - expected).abs() <= 2e-7 * expected.abs(),
                    "2^-{exponent}: {actual:e} != {expected:e}"
                );
            }
        }
        assert_eq!(
            law.stretching_midpoint_energy_remainder(&a, &a, &a),
            Some(0.0)
        );
    }

    #[test]
    fn malformed_film_work_refuses_before_indexing_or_publication() {
        let law = specimen();
        assert!(
            law.stretching_midpoint_energy_remainder(&[], &[0.0; 2], &[0.0; 2])
                .is_none()
        );
        assert!(
            law.stretching_midpoint_energy_remainder(&[0.0; 2], &[0.0], &[0.0; 2])
                .is_none()
        );
        assert!(
            law.stretching_midpoint_energy_remainder(&[0.0; 2], &[0.0; 2], &[0.0])
                .is_none()
        );
        assert!(
            law.stretching_midpoint_energy_remainder(&[0.0; 2], &[f64::NAN, 0.0], &[0.0; 2])
                .is_none()
        );
        assert!(
            law.stretching_midpoint_energy_remainder(&[1e300; 2], &[-1e300; 2], &[0.0; 2])
                .is_none()
        );
    }
}
