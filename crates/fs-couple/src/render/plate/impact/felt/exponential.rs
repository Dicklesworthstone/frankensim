//! Original felt potential in an exponential discrete-gradient composition.
//! The spring remains nonlinear and its compaction history is frozen within
//! each solve. Kelvin storage belongs to the quadratic block; the complete
//! original pad still supplies its force, potential and accepted history.
use super::super::Pad;
use super::{FeltPad, WoolFeltState};
use fs_math::{dd::Dd, det};

impl FeltPad {
    // Integrate one ORIGINAL frozen-history segment without subtracting two
    // large endpoint potentials. The shared unilateral power evaluator keeps
    // its cubic remainder; Dd retains the smaller affine/work correction.
    fn segment_increment(
        &self,
        a: f64,
        b: f64,
        anchor: f64,
        origin: f64,
        scale: f64,
        stress: f64,
        power: f64,
    ) -> Dd {
        let normal =
            |value: f64| (Dd::from_f64(value) - Dd::from_f64(origin)) / Dd::from_f64(scale);
        let (exact_x, exact_y) = (normal(a), normal(b));
        let (x, y, m) = (
            exact_x.to_f64(),
            exact_y.to_f64(),
            (anchor - origin) / scale,
        );
        let coefficient = self.area_m2 * self.thickness_m * stress * scale;
        let remainder =
            fs_dcontact::power_potential_midpoint_remainder(x, y, m, power, coefficient);
        // Normalize in two lanes. Otherwise independent endpoint divisions
        // leave O(ulp*H) work that overwhelms a small cubic remainder.
        let endpoint_correction = Dd::from_f64(coefficient)
            * (Dd::from_f64(det::pow(y.max(0.0), power)) * (exact_y - Dd::from_f64(y))
                - Dd::from_f64(det::pow(x.max(0.0), power)) * (exact_x - Dd::from_f64(x)));
        Dd::from_f64(remainder)
            + Dd::from_f64(coefficient * det::pow(m.max(0.0), power)) * Dd::from_f64(y - x)
            + endpoint_correction
    }

    /// H_felt(b)-H_felt(a)-F_felt(mid)*delta_compression, with frozen history.
    /// Splitting only at the original conditioning maximum retains the exact
    /// source law. Open gaps and reversible/virgin crossings use the same
    /// unilateral power primitive; no fitted tangent or loss is introduced.
    pub(super) fn midpoint_remainder(
        &self,
        a: f64,
        b: f64,
        mid: f64,
        compression_increment: Dd,
        history: &WoolFeltState,
    ) -> f64 {
        let (lo, hi, sign) = if a <= b { (a, b, 1.0) } else { (b, a, -1.0) };
        let peak = history.eps_max;
        let mut increment = Dd::ZERO;
        if peak > 0.0 && lo < peak {
            let end = hi.min(peak);
            let residual = self.law.eps_residual(history);
            increment = increment
                + self.segment_increment(
                    lo,
                    end,
                    mid.min(peak),
                    residual,
                    peak - residual,
                    history.sig_max,
                    self.law.q,
                );
        }
        if hi > peak {
            increment = increment
                + self.segment_increment(
                    lo.max(peak),
                    hi,
                    mid.max(peak),
                    0.0,
                    self.law.eps_ref,
                    self.law.f_ref,
                    self.law.p,
                );
        }
        (Dd::from_f64(sign) * increment
            - Dd::from_f64(self.force(mid, history)) * compression_increment)
            .to_f64()
    }
}

impl Pad {
    pub(crate) fn midpoint_remainder(
        &self,
        a: &[f64],
        b: &[f64],
        mid: &[f64],
        modes: usize,
        history: &WoolFeltState,
    ) -> f64 {
        let mut compression_increment = Dd::ZERO;
        for (i, weight) in self.spec.weights.iter().enumerate() {
            compression_increment = compression_increment
                + Dd::from_f64(*weight) * (Dd::from_f64(b[2 * i]) - Dd::from_f64(a[2 * i]));
        }
        for (i, branch) in self.spec.creep.iter().enumerate() {
            let at = 2 * modes + self.creep_start + i;
            compression_increment = compression_increment
                - (Dd::from_f64(b[at]) - Dd::from_f64(a[at]))
                    / Dd::from_f64(branch.stiffness_n_m.sqrt());
        }
        self.spec.midpoint_remainder(
            self.strain(a, modes),
            self.strain(b, modes),
            self.strain(mid, modes),
            compression_increment,
            history,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_material::fiber::felt_rate::{CreepBranch, MAX_FELT_RATE_BRANCHES, RateFelt};
    use fs_material::fiber::{Uniaxial, WoolFelt};

    fn pad() -> FeltPad {
        FeltPad {
            area_m2: 0.125,
            thickness_m: 0.5,
            precompression_m: 0.0,
            weights: vec![1.0],
            law: WoolFelt::new(128.0, 0.25, 2.0, 3.0, 0.25, 0.875).unwrap(),
            prior_maximum_strain: 0.0,
            creep: vec![],
        }
    }

    // Integer powers supply an independent exact original potential in Dd.
    fn energy(p: &FeltPad, e: f64, h: &WoolFeltState) -> Dd {
        let v = Dd::from_f64(p.area_m2 * p.thickness_m);
        let residual = Dd::from_f64(p.law.crush_fraction * h.eps_max);
        let span = Dd::from_f64(h.eps_max) - residual;
        let value = Dd::from_f64(e.min(h.eps_max));
        let x = if e.min(h.eps_max) > p.law.crush_fraction * h.eps_max {
            (value - residual) / span
        } else {
            Dd::ZERO
        };
        let mut u = v * Dd::from_f64(h.sig_max) * span * x * x * x * x / Dd::from_f64(4.0);
        if e > h.eps_max {
            let a = Dd::from_f64(h.eps_max / p.law.eps_ref);
            let b = Dd::from_f64(e / p.law.eps_ref);
            u = u + v * Dd::from_f64(p.law.f_ref * p.law.eps_ref) * (b * b * b - a * a * a)
                / Dd::from_f64(3.0);
        }
        u
    }

    #[test]
    fn g1_original_felt_remainder_keeps_open_conditioned_and_virgin_work() {
        let p = pad();
        let h = p.history_at(0.5);
        for (a, b) in [
            (-0.1, -0.02),
            (-0.1, 0.3),
            (0.2, 0.4),
            (0.49, 0.51),
            (0.6, 0.7),
        ] {
            for (a, b) in [(a, b), (b, a)] {
                let mid = f64::midpoint(a, b);
                let increment = Dd::from_f64(p.thickness_m) * (Dd::from_f64(b) - Dd::from_f64(a));
                let actual = p.midpoint_remainder(a, b, mid, increment, &h);
                let expected = (energy(&p, b, &h)
                    - energy(&p, a, &h)
                    - Dd::from_f64(p.force(mid, &h)) * increment)
                    .to_f64();
                let scale = energy(&p, a, &h)
                    .abs()
                    .to_f64()
                    .max(energy(&p, b, &h).abs().to_f64())
                    .max(1e-100);
                assert!(
                    (actual - expected).abs() < 2e-14 * scale,
                    "{a} -> {b}: {actual:e} vs {expected:e}"
                );
            }
        }
    }

    #[test]
    fn g1_small_felt_steps_retain_the_cubic_original_energy_remainder() {
        let p = pad();
        let h = p.history_at(0.5);
        for (a, exponents) in [(0.25, [8, 12, 16]), (0.625, [8, 12, 16])] {
            for exponent in exponents {
                let b = a + 2.0_f64.powi(-exponent);
                let mid = f64::midpoint(a, b);
                let increment = Dd::from_f64(p.thickness_m) * (Dd::from_f64(b) - Dd::from_f64(a));
                let actual = p.midpoint_remainder(a, b, mid, increment, &h);
                let expected = (energy(&p, b, &h)
                    - energy(&p, a, &h)
                    - Dd::from_f64(p.force(mid, &h)) * increment)
                    .to_f64();
                assert!(expected > 0.0);
                assert!(
                    (actual - expected).abs() < 2e-7 * expected,
                    "{exponent}: {actual:e} vs {expected:e}"
                );
            }
        }
    }

    #[test]
    fn g3_intrinsic_material_binding_preserves_all_eight_arms_and_physical_work() {
        let original = pad();
        let branches: Vec<_> = (0..MAX_FELT_RATE_BRANCHES)
            .map(|i| CreepBranch {
                modulus_pa: 128.0 * (i + 1) as f64,
                viscosity_pa_s: 32.0 * (i + 1) as f64,
            })
            .collect();
        let material = RateFelt::new(original.law.clone(), branches.clone()).unwrap();
        let whole = FeltPad::from_material(0.125, 0.5, 0.01, vec![1.0], 0.25, &material).unwrap();
        let half = FeltPad::from_material(0.0625, 0.5, 0.01, vec![1.0], 0.25, &material).unwrap();
        assert_eq!(whole.creep.len(), MAX_FELT_RATE_BRANCHES);
        assert_eq!(whole.law, *material.envelope());
        for ((source, mechanical), part) in branches.iter().zip(&whole.creep).zip(&half.creep) {
            let strain = 0.03125;
            let strain_rate = 0.125;
            let displacement = whole.thickness_m * strain;
            let velocity = whole.thickness_m * strain_rate;
            assert_eq!(
                mechanical.stiffness_n_m,
                source.modulus_pa * whole.area_m2 / whole.thickness_m
            );
            assert_eq!(
                mechanical.viscosity_n_s_m,
                source.viscosity_pa_s * whole.area_m2 / whole.thickness_m
            );
            assert_eq!(mechanical.stiffness_n_m, 2.0 * part.stiffness_n_m);
            assert_eq!(
                0.5 * mechanical.stiffness_n_m * displacement * displacement,
                whole.area_m2 * whole.thickness_m * 0.5 * source.modulus_pa * strain * strain
            );
            assert_eq!(
                mechanical.viscosity_n_s_m * velocity * velocity,
                whole.area_m2
                    * whole.thickness_m
                    * source.viscosity_pa_s
                    * strain_rate
                    * strain_rate
            );
        }
    }

    #[test]
    fn g0_empty_material_retains_legacy_felt_bits_and_bad_geometry_refuses() {
        let original = pad();
        let material = RateFelt::new(original.law.clone(), vec![]).unwrap();
        let bound = FeltPad::from_material(
            original.area_m2,
            original.thickness_m,
            original.precompression_m,
            original.weights.clone(),
            0.0,
            &material,
        )
        .unwrap();
        assert!(bound.creep.is_empty());
        for peak in [0.0, 0.1, 0.5] {
            let h = original.history_at(peak);
            for e in [-0.1, 0.0, 0.03, 0.2, 0.6] {
                assert_eq!(
                    bound.force(e, &h).to_bits(),
                    original.force(e, &h).to_bits()
                );
                assert_eq!(
                    bound.path_energy(e, &h).to_bits(),
                    original.path_energy(e, &h).to_bits()
                );
                assert_eq!(
                    bound.law.update_state(e, &h),
                    original.law.update_state(e, &h)
                );
            }
        }
        for (area, thickness) in [(0.0, 0.5), (0.1, 0.0), (f64::MAX, f64::MIN_POSITIVE)] {
            assert!(
                FeltPad::from_material(area, thickness, 0.0, vec![1.0], 0.0, &material).is_err()
            );
        }
        assert!(FeltPad::from_material(0.125, 0.5, 0.0, vec![], 0.0, &material).is_err());
    }
}
