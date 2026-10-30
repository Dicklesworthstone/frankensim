//! One shared vector-field quadrature for physical work and its full tangent.
use super::shape::{Evaluation, Work};
use super::*;

#[derive(Clone, Copy, Default)]
pub(super) struct Gradient {
    pub value: [f64; 3],
    /// Derivative with respect to each NEW physical coordinate (h,b0,b1).
    pub derivative: [[f64; 3]; 3],
    pub error: [f64; 3],
}
impl SourceAsymmetricCircularContactLineBridge {
    pub(super) fn discrete_gradient(
        &self,
        old: [f64; 3],
        new: [f64; 3],
        seed: Evaluation,
        work: &mut Work,
        gate: &CancelGate,
    ) -> Result<Gradient, ContactLineError> {
        if old == new {
            let e = self.at(old, seed, work, gate)?;
            return Ok(Gradient {
                value: e.gradient,
                derivative: e.hessian.map(|r| r.map(|x| 0.5 * x)),
                error: e.gradient_error,
            });
        }
        let coarse = self.quadrature(old, new, seed, &GAUSS8, work, gate)?;
        let mut fine = self.quadrature(old, new, seed, &GAUSS16, work, gate)?;
        for i in 0..3 {
            fine.error[i] +=
                (fine.value[i] - coarse.value[i]).abs() + 64.0 * f64::EPSILON * fine.value[i].abs();
        }
        Ok(fine)
    }
    fn quadrature(
        &self,
        old: [f64; 3],
        new: [f64; 3],
        mut seed: Evaluation,
        nodes: &[(f64, f64)],
        work: &mut Work,
        gate: &CancelGate,
    ) -> Result<Gradient, ContactLineError> {
        let mut sum = Gradient::default();
        for &(x, w) in nodes {
            for sign in [-1.0, 1.0] {
                let t = 0.5 * (1.0 + sign * x);
                let point = core::array::from_fn(|i| old[i] + t * (new[i] - old[i]));
                seed = self.at(point, seed, work, gate)?;
                for i in 0..3 {
                    sum.value[i] += 0.5 * w * seed.gradient[i];
                    sum.error[i] += 0.5 * w * seed.gradient_error[i];
                    for j in 0..3 {
                        sum.derivative[i][j] += 0.5 * w * t * seed.hessian[i][j];
                    }
                }
            }
        }
        Ok(sum)
    }
}
const GAUSS8: [(f64, f64); 4] = [
    (0.183_434_642_495_649_8, 0.362_683_783_378_362),
    (0.525_532_409_916_329, 0.313_706_645_877_887_3),
    (0.796_666_477_413_626_7, 0.222_381_034_453_374_5),
    (0.960_289_856_497_536_3, 0.101_228_536_290_376_26),
];
const GAUSS16: [(f64, f64); 8] = [
    (0.095_012_509_837_637_44, 0.189_450_610_455_068_5),
    (0.281_603_550_779_258_9, 0.182_603_415_044_923_6),
    (0.458_016_777_657_227_4, 0.169_156_519_395_002_54),
    (0.617_876_244_402_643_8, 0.149_595_988_816_576_73),
    (0.755_404_408_355_003, 0.124_628_971_255_533_87),
    (0.865_631_202_387_831_8, 0.095_158_511_682_492_78),
    (0.944_575_023_073_232_6, 0.062_253_523_938_647_89),
    (0.989_400_934_991_649_9, 0.027_152_459_411_754_095),
];
