//! Fixed-order positive Gauss rules on every source segment. Threshold shifting
//! factors out the rare-event Boltzmann exponential before integration.

use super::*;

pub(super) struct Moments {
    pub value: [f64; 3],
    pub relative_error: f64,
    pub relative_tail: f64,
}

const X8: [f64; 4] = [0.183_434_642_495_649_8, 0.525_532_409_916_329,
    0.796_666_477_413_626_7, 0.960_289_856_497_536_3];
const W8: [f64; 4] = [0.362_683_783_378_362, 0.313_706_645_877_887_3,
    0.222_381_034_453_374_5, 0.101_228_536_290_376_3];
const X16: [f64; 8] = [0.095_012_509_837_637_44, 0.281_603_550_779_258_9,
    0.458_016_777_657_227_4, 0.617_876_244_402_643_8, 0.755_404_408_355_003,
    0.865_631_202_387_831_8, 0.944_575_023_073_232_6, 0.989_400_934_991_649_9];
const W16: [f64; 8] = [0.189_450_610_455_068_5, 0.182_603_415_044_923_6,
    0.169_156_519_395_002_5, 0.149_595_988_816_576_7, 0.124_628_971_255_533_9,
    0.095_158_511_682_492_79, 0.062_253_523_938_647_89, 0.027_152_459_411_754_1];

#[derive(Default)]
struct Sum { value: f64, correction: f64 }
impl Sum {
    fn add(&mut self, x: f64) {
        let y = x - self.correction;
        let total = self.value + y;
        self.correction = (total - self.value) - y;
        self.value = total;
    }
}

pub(super) fn moments(
    source: &ElectronCrossSectionV1,
    thermal: f64,
    threshold_ratio: f64,
    power: usize,
    controls: MaxwellianRateControlsV1,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Moments, ElectronCollisionErrorV1> {
    let source_end = (source.maximum_energy.value() - source.threshold.value()) / thermal;
    if !positive(source_end) {
        return Err(ElectronCollisionErrorV1::Arithmetic("covered Maxwellian energy interval"));
    }
    let mut end = source_end.min(32.0);
    let mut panels_used = 0;
    loop {
        let mut refinement = 1_usize;
        let (values, error) = loop {
            let mut values: [Sum; 3] = core::array::from_fn(|_| Sum::default());
            let mut errors: [Sum; 3] = core::array::from_fn(|_| Sum::default());
            let mut segment = |lower: f64, upper: f64| -> Result<(), ElectronCollisionErrorV1> {
                let lower = lower.min(end);
                let upper = upper.min(end);
                if lower >= upper { return Ok(()); }
                let count = (((upper - lower) / 4.0).ceil() as usize).max(1)
                    .checked_mul(refinement).ok_or(ElectronCollisionErrorV1::QuadratureBudget)?;
                panels_used += count;
                if panels_used > controls.maximum_subdivisions {
                    return Err(ElectronCollisionErrorV1::QuadratureBudget);
                }
                let width = (upper - lower) / count as f64;
                for panel in 0..count {
                    if cancelled() { return Err(ElectronCollisionErrorV1::Cancelled); }
                    let a = lower + panel as f64 * width;
                    let b = if panel + 1 == count { upper } else { lower + (panel + 1) as f64 * width };
                    let fine = gauss(source, thermal, threshold_ratio, a, b, &X16, &W16);
                    let coarse = gauss(source, thermal, threshold_ratio, a, b, &X8, &W8);
                    for i in 0..3 {
                        values[i].add(fine[i]);
                        errors[i].add((fine[i] - coarse[i]).abs());
                    }
                }
                Ok(())
            };
            match &source.law {
                CrossSectionLaw::Beb { .. } => segment(0.0, end)?,
                CrossSectionLaw::Linear(points) => {
                    for pair in points.windows(2) {
                        let lower = (pair[0].energy.value() - source.threshold.value()) / thermal;
                        if lower >= end { break; }
                        let upper = (pair[1].energy.value() - source.threshold.value()) / thermal;
                        segment(lower, upper)?;
                    }
                }
            }
            let values = values.map(|sum| sum.value);
            if values.iter().any(|v| !nonnegative(*v)) {
                return Err(ElectronCollisionErrorV1::Arithmetic("Maxwellian integral moment"));
            }
            // Both moments must resolve: a small rate-error estimate alone
            // does not establish the temperature derivative.
            let error = (power - 1..=power).map(|i| {
                if values[i] == 0.0 { 0.0 } else { errors[i].value / values[i] }
            }).fold(0.0, f64::max);
            if error <= controls.relative_tolerance * 0.25 {
                break (values, error);
            }
            refinement = refinement.checked_mul(2).ok_or(ElectronCollisionErrorV1::QuadratureBudget)?;
        };
        let envelope_ratio = if end == source_end {
            source.high_energy_upper_bound.value() / source.scale
        } else { 1.0 };
        let tail = tail_moments(end, threshold_ratio).map(|v| v * envelope_ratio);
        let relative_tail = (power - 1..=power).map(|i| {
            if values[i] == 0.0 { f64::INFINITY } else { tail[i] / values[i] }
        }).fold(0.0, f64::max);
        if relative_tail <= controls.relative_tolerance * 0.25 {
            return Ok(Moments { value: values, relative_error: error, relative_tail });
        }
        if end == source_end || end >= 1024.0 {
            return Err(ElectronCollisionErrorV1::UnresolvedTail { relative_bound: relative_tail });
        }
        end = (end * 2.0).min(source_end);
    }
}

fn gauss<const N: usize>(
    source: &ElectronCrossSectionV1,
    thermal: f64,
    a: f64,
    lower: f64,
    upper: f64,
    nodes: &[f64; N],
    weights: &[f64; N],
) -> [f64; 3] {
    let half = 0.5 * (upper - lower);
    let midpoint = lower + half;
    let mut values: [Sum; 3] = core::array::from_fn(|_| Sum::default());
    for (&node, &weight) in nodes.iter().zip(weights) {
        for sign in [-1.0, 1.0] {
            let y = midpoint + sign * half * node;
            let energy = (source.threshold.value() + thermal * y).min(source.maximum_energy.value());
            let x = (a + y) / (a + 1.0);
            let weighted = weight * (source.value(energy) / source.scale) * exp(-y);
            let mut term = weighted;
            for value in &mut values {
                term *= x;
                value.add(term);
            }
        }
    }
    values.map(|sum| sum.value * half)
}

// Exact incomplete-gamma tails of x^p exp(-x), scaled by exp(a)*(a+1)^-p.
// Their positivity avoids cancellation even when the ionization threshold is
// many thermal energies above zero and the direct physical rate underflows.
fn tail_moments(end: f64, a: f64) -> [f64; 3] {
    let inv = 1.0 / (a + 1.0);
    let b = (a + end) * inv;
    let e = exp(-end);
    [e * (b + inv),
     e * (b * b + 2.0 * inv * b + 2.0 * inv * inv),
     e * (b * b * b + 3.0 * inv * b * b + 6.0 * inv * inv * b + 6.0 * inv * inv * inv)]
}
