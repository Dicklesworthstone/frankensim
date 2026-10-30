//! Symmetric coordinate-increment discrete gradient, without energy secants.
use super::shape::{Evaluation, Work};
use super::*;

#[derive(Clone, Copy, Default)]
pub(super) struct Gradient {
    pub value: f64,
    /// Derivatives with respect to the NEW gap and NEW radius.
    pub derivative: [f64; 2],
    pub error: f64,
}

impl SourceAxisymmetricContactLineBridge {
    pub(super) fn discrete_gradient(
        &self,
        coordinate: usize,
        old: [f64; 2],
        new: [f64; 2],
        seed: Evaluation,
        work: &mut Work,
        gate: &CancelGate,
    ) -> Result<Gradient, ContactLineError> {
        let other = 1 - coordinate;
        let left = self.integral(
            coordinate,
            old[other],
            old[coordinate],
            new[coordinate],
            seed,
            work,
            gate,
        )?;
        let right = if old[other] == new[other] {
            left
        } else {
            self.integral(
                coordinate,
                new[other],
                old[coordinate],
                new[coordinate],
                seed,
                work,
                gate,
            )?
        };
        let mut derivative = [0.0; 2];
        derivative[coordinate] = 0.5 * (left.derivative[coordinate] + right.derivative[coordinate]);
        derivative[other] = 0.5 * right.derivative[other];
        Ok(Gradient {
            value: 0.5 * (left.value + right.value),
            derivative,
            error: 0.5 * (left.error + right.error),
        })
    }

    fn integral(
        &self,
        coordinate: usize,
        fixed: f64,
        start: f64,
        end: f64,
        seed: Evaluation,
        work: &mut Work,
        gate: &CancelGate,
    ) -> Result<Gradient, ContactLineError> {
        if start == end {
            let point = if coordinate == 0 {
                [start, fixed]
            } else {
                [fixed, start]
            };
            let e = self.at(point[0], point[1], seed, work, gate)?;
            let mut derivative = e.hessian[coordinate];
            derivative[coordinate] *= 0.5;
            return Ok(Gradient {
                value: e.gradient[coordinate],
                derivative,
                error: e.gradient_error[coordinate],
            });
        }
        let coarse = self.quadrature(coordinate, fixed, start, end, seed, &GAUSS8, work, gate)?;
        let mut fine =
            self.quadrature(coordinate, fixed, start, end, seed, &GAUSS16, work, gate)?;
        fine.error += (fine.value - coarse.value).abs() + 64.0 * f64::EPSILON * fine.value.abs();
        Ok(fine)
    }

    fn quadrature(
        &self,
        coordinate: usize,
        fixed: f64,
        start: f64,
        end: f64,
        mut seed: Evaluation,
        nodes: &[(f64, f64)],
        work: &mut Work,
        gate: &CancelGate,
    ) -> Result<Gradient, ContactLineError> {
        let mut sum = Gradient::default();
        for &(x, weight) in nodes {
            for sign in [-1.0, 1.0] {
                let t = 0.5 * (1.0 + sign * x);
                let value = start + t * (end - start);
                let point = if coordinate == 0 {
                    [value, fixed]
                } else {
                    [fixed, value]
                };
                seed = self.at(point[0], point[1], seed, work, gate)?;
                sum.value += 0.5 * weight * seed.gradient[coordinate];
                for i in 0..2 {
                    sum.derivative[i] += 0.5
                        * weight
                        * seed.hessian[coordinate][i]
                        * if i == coordinate { t } else { 1.0 };
                }
                sum.error += 0.5 * weight * seed.gradient_error[coordinate];
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
