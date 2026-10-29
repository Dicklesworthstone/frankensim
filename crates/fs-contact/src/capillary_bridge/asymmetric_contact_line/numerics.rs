//! Symmetric gap/footprint-block discrete gradient and analytic Jacobian.
use super::{ContactLineError, SourceAsymmetricContactLineBridge};
use fs_exec::CancelGate;

#[derive(Clone, Copy, Default)]
pub(super) struct FootprintGradient {
    pub value: [f64; 2],
    pub jacobian: [[f64; 2]; 2],
    pub gap_derivative: [f64; 2],
    pub error: [f64; 2],
}
#[derive(Clone, Copy, Default)]
pub(super) struct GapGradient {
    pub value: f64,
    pub gap_derivative: f64,
    pub span_derivatives: [f64; 2],
    pub error: f64,
}

impl SourceAsymmetricContactLineBridge {
    pub(super) fn footprint_gradient(
        &self,
        h0: f64,
        h1: f64,
        b0: [f64; 2],
        b1: [f64; 2],
        evaluations: &mut usize,
        gate: &CancelGate,
    ) -> Result<FootprintGradient, ContactLineError> {
        let left = self.footprint_integral(h0, b0, b1, evaluations, gate)?;
        let right = if h0 == h1 {
            left
        } else {
            self.footprint_integral(h1, b0, b1, evaluations, gate)?
        };
        Ok(FootprintGradient {
            value: core::array::from_fn(|i| 0.5 * (left.value[i] + right.value[i])),
            jacobian: core::array::from_fn(|i| {
                core::array::from_fn(|j| 0.5 * (left.jacobian[i][j] + right.jacobian[i][j]))
            }),
            gap_derivative: right.gap_derivative.map(|v| 0.5 * v),
            error: core::array::from_fn(|i| 0.5 * (left.error[i] + right.error[i])),
        })
    }
    fn footprint_integral(
        &self,
        h: f64,
        b0: [f64; 2],
        b1: [f64; 2],
        evaluations: &mut usize,
        gate: &CancelGate,
    ) -> Result<FootprintGradient, ContactLineError> {
        if b0 == b1 {
            let e = self.at(h, b0, evaluations, gate)?;
            return Ok(FootprintGradient {
                value: [e.gradient[1], e.gradient[2]],
                jacobian: core::array::from_fn(|i| {
                    core::array::from_fn(|j| 0.5 * e.hessian[i + 1][j + 1])
                }),
                gap_derivative: [e.hessian[1][0], e.hessian[2][0]],
                error: [e.gradient_error[1], e.gradient_error[2]],
            });
        }
        let coarse = self.footprint_quadrature(h, b0, b1, &GAUSS8, evaluations, gate)?;
        let mut fine = self.footprint_quadrature(h, b0, b1, &GAUSS16, evaluations, gate)?;
        for i in 0..2 {
            fine.error[i] +=
                (fine.value[i] - coarse.value[i]).abs() + 64.0 * f64::EPSILON * fine.value[i].abs();
        }
        Ok(fine)
    }
    fn footprint_quadrature(
        &self,
        h: f64,
        b0: [f64; 2],
        b1: [f64; 2],
        nodes: &[(f64, f64)],
        evaluations: &mut usize,
        gate: &CancelGate,
    ) -> Result<FootprintGradient, ContactLineError> {
        let mut sum = FootprintGradient::default();
        for &(x, w) in nodes {
            for sign in [-1.0, 1.0] {
                let t = 0.5 * (1.0 + sign * x);
                let e = self.at(
                    h,
                    core::array::from_fn(|i| b0[i] + t * (b1[i] - b0[i])),
                    evaluations,
                    gate,
                )?;
                for i in 0..2 {
                    sum.value[i] += 0.5 * w * e.gradient[i + 1];
                    sum.gap_derivative[i] += 0.5 * w * e.hessian[i + 1][0];
                    sum.error[i] += 0.5 * w * e.gradient_error[i + 1];
                    for j in 0..2 {
                        sum.jacobian[i][j] += 0.5 * w * t * e.hessian[i + 1][j + 1];
                    }
                }
            }
        }
        Ok(sum)
    }

    pub(super) fn gap_gradient(
        &self,
        h0: f64,
        h1: f64,
        b0: [f64; 2],
        b1: [f64; 2],
        evaluations: &mut usize,
        gate: &CancelGate,
    ) -> Result<GapGradient, ContactLineError> {
        let left = self.gap_integral(h0, h1, b0, evaluations, gate)?;
        let right = if b0 == b1 {
            left
        } else {
            self.gap_integral(h0, h1, b1, evaluations, gate)?
        };
        Ok(GapGradient {
            value: 0.5 * (left.value + right.value),
            gap_derivative: 0.5 * (left.gap_derivative + right.gap_derivative),
            span_derivatives: right.span_derivatives.map(|v| 0.5 * v),
            error: 0.5 * (left.error + right.error),
        })
    }
    fn gap_integral(
        &self,
        h0: f64,
        h1: f64,
        b: [f64; 2],
        evaluations: &mut usize,
        gate: &CancelGate,
    ) -> Result<GapGradient, ContactLineError> {
        if h0 == h1 {
            let e = self.at(h0, b, evaluations, gate)?;
            return Ok(GapGradient {
                value: e.gradient[0],
                gap_derivative: 0.5 * e.hessian[0][0],
                span_derivatives: [e.hessian[0][1], e.hessian[0][2]],
                error: e.gradient_error[0],
            });
        }
        let coarse = self.gap_quadrature(h0, h1, b, &GAUSS8, evaluations, gate)?;
        let mut fine = self.gap_quadrature(h0, h1, b, &GAUSS16, evaluations, gate)?;
        fine.error += (fine.value - coarse.value).abs() + 64.0 * f64::EPSILON * fine.value.abs();
        Ok(fine)
    }
    fn gap_quadrature(
        &self,
        h0: f64,
        h1: f64,
        b: [f64; 2],
        nodes: &[(f64, f64)],
        evaluations: &mut usize,
        gate: &CancelGate,
    ) -> Result<GapGradient, ContactLineError> {
        let mut sum = GapGradient::default();
        for &(x, w) in nodes {
            for sign in [-1.0, 1.0] {
                let t = 0.5 * (1.0 + sign * x);
                let e = self.at(h0 + t * (h1 - h0), b, evaluations, gate)?;
                sum.value += 0.5 * w * e.gradient[0];
                sum.gap_derivative += 0.5 * w * t * e.hessian[0][0];
                sum.error += 0.5 * w * e.gradient_error[0];
                for i in 0..2 {
                    sum.span_derivatives[i] += 0.5 * w * e.hessian[0][i + 1];
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
