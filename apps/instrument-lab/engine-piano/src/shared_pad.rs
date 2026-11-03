//! One physical viscous pad station acting across two modal-state owners.
//! The completed-square modal shapes and bridge lift form one rank-one loss
//! vector. Keeping that vector whole preserves its contact-velocity nullspace.

use fs_math::det;

/// The bank stores other modal state beside velocity; the active prefix uses
/// a contiguous velocity slice. This projection keeps the shared flow generic
/// without packing or allocating on the audio path.
pub(crate) trait ModalVelocity {
    fn velocity(&self) -> f64;
    fn add_velocity(&mut self, increment: f64);
}

impl ModalVelocity for f64 {
    fn velocity(&self) -> f64 { *self }
    fn add_velocity(&mut self, increment: f64) { *self += increment; }
}

#[derive(Clone, Debug)]
pub(crate) struct PadPoint {
    pub(crate) shape: Vec<f64>,
    pub(crate) lift: f64,
    pub(crate) norm_string: f64,
    /// expm1(-drag*norm*duration)/norm for the current engagement and basis.
    pub(crate) coef: f64,
    pub(crate) norm: f64,
}

impl PadPoint {
    /// Refresh when engagement or the loaded board basis changes. Finite
    /// nonnegative drag is admitted by the physical pad specification.
    pub(crate) fn prepare(&mut self, drag: f64, bridge_norm2: f64, duration: f64) {
        self.norm = self.norm_string + self.lift * self.lift * bridge_norm2;
        self.coef = if self.norm == 0.0 || drag == 0.0 { 0.0 }
            else { det::expm1(-drag * self.norm * duration) / self.norm };
    }

    /// Exact rank-one viscous flow. The active prefix and bank tail use the
    /// same modal order regardless of their ownership boundary. Returns the
    /// impulse multiplying the bridge row and the physical dissipated heat.
    /// The caller carries that impulse into the bridge velocity before another
    /// station acts, and may accumulate its vector update once per course.
    pub(crate) fn apply<T: ModalVelocity>(&self, active: &mut [f64], bank: &mut [T], bridge_speed: f64) -> (f64, f64) {
        debug_assert_eq!(self.shape.len(), active.len() + bank.len());
        if self.coef == 0.0 {
            return (0.0, 0.0);
        }
        let speed = self.shape.iter().zip(active.iter().copied().chain(bank.iter().map(ModalVelocity::velocity)))
            .map(|(g, v)| g * v).sum::<f64>() + self.lift * bridge_speed;
        let change = self.coef * speed;
        for (v, g) in active.iter_mut().zip(&self.shape) {
            *v += g * change;
        }
        for (v, g) in bank.iter_mut().zip(&self.shape[active.len()..]) {
            v.add_velocity(g * change);
        }
        // coef*norm lies in [-1,0]. This factorized exact integral stays
        // nonnegative even when the heat is below an energy subtraction's ulp.
        let heat = -change * speed * (1.0 + 0.5 * self.coef * self.norm);
        (self.lift * change, heat)
    }
}

#[cfg(test)]
mod tests {
    use super::PadPoint;

    fn station(exponent: f64) -> PadPoint {
        // One bridge coordinate: the physical station vector is [1,2,2].
        let mut p = PadPoint {
            shape: vec![1.0, 2.0], lift: 2.0, norm_string: 5.0, coef: 0.0, norm: 0.0,
        };
        p.prepare(exponent / 9.0, 1.0, 1.0);
        p
    }

    fn advance(p: &PadPoint, initial: [f64; 3], split: usize) -> ([f64; 3], f64) {
        let mut active = initial[..split].to_vec();
        let mut bank = initial[split..2].to_vec();
        let (impulse, heat) = p.apply(&mut active, &mut bank, initial[2]);
        let mut velocity = initial;
        for (out, value) in velocity[..2].iter_mut().zip(active.iter().chain(&bank)) {
            *out = *value;
        }
        velocity[2] += impulse;
        (velocity, heat)
    }

    #[test]
    fn a_station_preserves_modal_and_bridge_contact_nullspaces() {
        for exponent in [1e-18, 1.0, 1_000.0] {
            let p = station(exponent);
            for initial in [[2.0, -1.0, 0.0], [2.0, 0.0, -1.0]] {
                for split in 0..=2 {
                    let (actual, heat) = advance(&p, initial, split);
                    assert_eq!(actual.map(f64::to_bits), initial.map(f64::to_bits));
                    assert_eq!(heat, 0.0);
                }
            }
        }
    }

    #[test]
    fn recoil_and_heat_match_the_continuous_rank_one_solution_at_every_split() {
        // g.v=3, g.g=9. The contact speed is 3*exp(-exponent*t),
        // so integral(drag*speed^2,0..1) = (1-exp(-2*exponent))/2.
        let initial = [3.0, -1.0, 1.0];
        let initial_energy = 5.5;
        for exponent in [1e-18, 1e-12, 1.0, 1_000.0] {
            let p = station(exponent);
            let delta = (-exponent).exp_m1() / 3.0;
            let expected = [3.0 + delta, -1.0 + 2.0 * delta, 1.0 + 2.0 * delta];
            let expected_heat = -0.5 * (-2.0 * exponent).exp_m1();
            let mut first = None;
            for split in 0..=2 {
                let (actual, heat) = advance(&p, initial, split);
                for (got, want) in actual.iter().zip(expected) {
                    assert!((got - want).abs() < 3e-14, "exponent {exponent}, split {split}: {actual:?}");
                }
                assert!(heat.is_finite() && heat > 0.0 && heat <= 0.5);
                assert!((heat - expected_heat).abs() < 3e-13 * expected_heat,
                    "exponent {exponent}, split {split}: heat {heat:e}, expected {expected_heat:e}");
                let final_energy = 0.5 * actual.iter().map(|v| v * v).sum::<f64>();
                assert!((final_energy + heat - initial_energy).abs() < 3e-14);
                let bits = (actual.map(f64::to_bits), heat.to_bits());
                if let Some(reference) = first {
                    assert_eq!(bits, reference, "ownership must not change the station flow");
                } else {
                    first = Some(bits);
                }
            }
        }
    }

    #[test]
    fn successive_half_flows_give_the_same_recoil_and_integrated_heat() {
        let full = station(0.8);
        let half = station(0.4);
        let initial = [0.75, -0.25, -0.5];
        for split in 0..=2 {
            let (once, heat) = advance(&full, initial, split);
            let (middle, first_heat) = advance(&half, initial, split);
            let (twice, second_heat) = advance(&half, middle, split);
            for (a, b) in once.iter().zip(twice) {
                assert!((a - b).abs() < 2e-15);
            }
            assert!((heat - first_heat - second_heat).abs() < 2e-15);
        }
    }
}
