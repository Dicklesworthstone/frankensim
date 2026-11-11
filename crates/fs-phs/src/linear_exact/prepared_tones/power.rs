//! Original port work and resistance squares along continuous source motion.
use super::super::{LinearSplit, finite, positive_heat, refuse};
use crate::interconnection::accumulation::Sum as Contraction;
use crate::{PhsError, PortHamiltonian, ResistanceWitness};
use fs_math::det;

struct Observation {
    rows: usize,
    dimension: usize,
    map: Vec<f64>,
    left: Vec<f64>,
    right: Option<Vec<f64>>,
    /// Doublings of the plain f64 chain that integrated this observation.
    doublings: u32,
}

pub(in crate::linear_exact) struct TonePower {
    n: usize,
    m: usize,
    sources: usize,
    stride: usize,
    rate: f64,
    root_rate: f64,
    inverse_root_rate: Contraction,
    scales: Vec<f64>,
    source_scales: Vec<f64>,
    work: Observation,
    heat: Vec<Observation>,
    /// Longest chain among the work and heat observations: the step's energy-balance gate reads it.
    doublings: u32,
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), crate::PreparedStepError> {
    if cancelled() {
        Err(crate::PreparedStepError::Cancelled)
    } else {
        Ok(())
    }
}

impl TonePower {
    #[allow(clippy::too_many_arguments)]
    pub(in crate::linear_exact) fn new(
        system: &PortHamiltonian,
        split: &LinearSplit,
        frequencies: &[f64],
        source_scales: &[f64],
        dt: f64,
        affine: bool,
        channels: &[ResistanceWitness],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, crate::PreparedStepError> {
        let (n, m) = (system.n, system.m);
        let sources = source_scales.len();
        let stride = if affine { 4 } else { 2 };
        let scales = (0..n)
            .map(|i| split.coordinate_scale(i))
            .collect::<Vec<_>>();
        let mut a = vec![0.0; n * n];
        for i in 0..n {
            poll(cancelled)?;
            for j in 0..n {
                let mut value = Contraction::default();
                let refused = || refuse("continuous-tone physical generator overflows");
                for k in 0..n {
                    let q = split.quadratic.matrix[k * n + j];
                    value.product(system.j[i * n + k], q).ok_or_else(refused)?;
                    value.product(-system.r[i * n + k], q).ok_or_else(refused)?;
                }
                a[i * n + j] = value.value().ok_or_else(refused)? * scales[i] / scales[j];
            }
        }
        let mut rate = frequencies.iter().copied().fold(1.0_f64, f64::max);
        for i in 0..n {
            rate = rate
                .max((0..n).map(|j| a[i * n + j].abs()).sum())
                .max((0..n).map(|j| a[j * n + i].abs()).sum());
        }
        if !finite(&a) || !rate.is_finite() {
            return Err(refuse("continuous-tone physical generator/rate is nonfinite").into());
        }
        // Retain the division remainder once, during preparation. Warm effort
        // projections then carry both words through the unchanged rate map.
        let root_rate = det::sqrt(rate);
        let reciprocal = 1.0 / root_rate;
        let mut inverse_root_rate = Contraction::default();
        inverse_root_rate.product(1.0, reciprocal)
            .and_then(|()| inverse_root_rate.product(
                1.0, (-reciprocal).mul_add(root_rate, 1.0) / root_rate,
            ))
            .ok_or_else(|| refuse("continuous-tone observation scale overflows"))?;
        let mut output = vec![0.0; m * n];
        for p in 0..m {
            for i in 0..n {
                output[p * n + i] = system.g[i * m + p];
            }
        }
        let work = prepare_observation(
            system,
            split,
            frequencies,
            source_scales,
            &scales,
            &a,
            rate,
            dt,
            affine,
            output,
            m,
            true,
        )?;
        poll(cancelled)?;
        let original_root = if let Some(witness) = &system.resistance_witness {
            if witness.matrix == system.r {
                witness.root.clone()
            } else {
                positive_heat::factor_resistance(&system.r, n)?
            }
        } else {
            positive_heat::factor_resistance(&system.r, n)?
        };
        let mut heat = Vec::with_capacity(channels.len() + 1);
        for index in 0..=channels.len() {
            poll(cancelled)?;
            let (matrix, root) = if index == 0 {
                (&system.r, original_root.clone())
            } else {
                (
                    &channels[index - 1].matrix,
                    channels[index - 1].root.clone(),
                )
            };
            positive_heat::check_resistance_root(matrix, &root, n)?;
            heat.push(prepare_observation(
                system,
                split,
                frequencies,
                source_scales,
                &scales,
                &a,
                rate,
                dt,
                affine,
                root,
                n,
                false,
            )?);
        }
        poll(cancelled)?;
        let doublings = heat.iter().fold(work.doublings, |s, o| s.max(o.doublings));
        Ok(Self {
            n,
            m,
            sources,
            stride,
            rate,
            root_rate,
            inverse_root_rate,
            scales,
            source_scales: source_scales.to_vec(),
            work,
            heat,
            doublings,
        })
    }

    /// Doublings of the longest plain f64 chain behind the supplied work and dissipated heat.
    pub(in crate::linear_exact) fn doublings(&self) -> u32 {
        self.doublings
    }

    /// Leading and remainder planes for the largest original observation.
    /// This is evaluation scratch, not another physical or propagation state.
    pub(in crate::linear_exact) fn coordinate_len(&self) -> usize {
        2 * (self.n + self.n.max(self.m) + self.sources)
    }

    /// Uses the original physical effort and force imbalance before any
    /// factor contraction. This keeps a stationary resistance nullspace and
    /// tiny forced motion from becoming subtraction of large heat terms.
    pub(in crate::linear_exact) fn evaluate_into(
        &self,
        system: &PortHamiltonian,
        effort: &[f64],
        input: &[f64],
        derivative: &mut [f64],
        coordinates: &mut [f64],
        heat: &mut [f64],
    ) -> Result<f64, PhsError> {
        let tones = self.sources / (self.stride * self.m);
        for i in 0..self.n {
            // Keep the complete physical imbalance in one contraction. A weak
            // drive can survive only after large port/internal forces cancel;
            // rounding those groups first would erase its entire motion.
            let mut value = Contraction::default();
            let refused = || refuse("continuous-tone physical derivative overflows");
            for p in 0..self.m {
                for t in 0..tones {
                    value
                        .product(system.g[i * self.m + p], input[self.stride * (p * tones + t)])
                        .ok_or_else(refused)?;
                }
            }
            for (j, &e) in effort.iter().enumerate() {
                value.product(system.j[i * self.n + j], e).ok_or_else(refused)?;
                value.product(-system.r[i * self.n + j], e).ok_or_else(refused)?;
            }
            derivative[i] = value.value().ok_or_else(refused)? * self.scales[i] / self.rate;
        }
        let supplied = self.evaluate(&self.work, effort, derivative, input, None, coordinates)?;
        for (value, observation) in heat.iter_mut().zip(&self.heat) {
            *value = self.evaluate(observation, effort, derivative, input, None, coordinates)?;
        }
        Ok(supplied)
    }

    /// Integrate a declared source against the output of the COMPLETE physical
    /// trajectory. Only the right bilinear factor changes: source oscillators
    /// are autonomous, so that factor has no initial physical-state columns.
    /// This retains small external work even beside a large interface exchange.
    pub(super) fn source_work_into(
        &self,
        effort: &[f64],
        derivative: &[f64],
        physical_input: &[f64],
        source_input: &[f64],
        coordinates: &mut [f64],
    ) -> Result<f64, PhsError> {
        self.evaluate(
            &self.work,
            effort,
            derivative,
            physical_input,
            Some(source_input),
            coordinates,
        )
    }

    fn evaluate(
        &self,
        observation: &Observation,
        effort: &[f64],
        derivative: &[f64],
        input: &[f64],
        source_input: Option<&[f64]>,
        coordinates: &mut [f64],
    ) -> Result<f64, PhsError> {
        if observation.left.is_empty() {
            return Ok(0.0);
        }
        let Some(right) = observation.right.as_deref() else {
            return self.evaluate_heat(observation, effort, derivative, input, coordinates);
        };
        if source_input.unwrap_or(input).iter().all(|value| *value == 0.0) {
            return Ok(0.0);
        }
        let l = observation.rows;
        let d = observation.dimension;
        let coordinates = coordinates.get_mut(..2 * d)
            .ok_or_else(|| refuse("continuous-tone retained power scratch extent"))?;
        let (high, low) = coordinates.split_at_mut(d);
        let refused = || refuse("continuous-tone physical power contraction overflows");
        for i in 0..l {
            let projection = Contraction::dot(
                &observation.map[i * self.n..(i + 1) * self.n], effort,
            ).ok_or_else(refused)?;
            let mut value = Contraction::default();
            value.add_product(projection, self.inverse_root_rate).ok_or_else(refused)?;
            (high[i], low[i]) = value.components().ok_or_else(refused)?;
        }
        high[l..l + self.n].copy_from_slice(derivative);
        low[l..l + self.n].fill(0.0);
        for i in 0..self.sources {
            let mut value = Contraction::default();
            value.product(self.source_scales[i], input[i]).ok_or_else(refused)?;
            (high[l + self.n + i], low[l + self.n + i]) =
                value.components().ok_or_else(refused)?;
        }
        if !finite(high) || !finite(low) {
            return Err(refuse("continuous-tone power coordinates overflow"));
        }
        let mut integral = Contraction::default();
        for (i, row) in observation.left.chunks_exact(d).enumerate() {
            let mut value = Contraction::default();
            for (j, &a) in right[i * d..(i + 1) * d].iter().enumerate() {
                if a == 0.0 { continue; }
                let (leading, tail) = if let Some(source) = source_input.filter(|_| j >= l + self.n) {
                    let index = j - l - self.n;
                    let mut coordinate = Contraction::default();
                    coordinate.product(self.source_scales[index], source[index])
                        .ok_or_else(refused)?;
                    coordinate.components().ok_or_else(refused)?
                } else {
                    (high[j], low[j])
                };
                value.product(a, leading).ok_or_else(refused)?;
                if tail != 0.0 { value.product_remainder(a, tail).ok_or_else(refused)?; }
            }
            if value.components().ok_or_else(refused)? != (0.0, 0.0) {
                let mut left = Contraction::default();
                for ((&a, &x), &tail) in row.iter().zip(high.iter()).zip(low.iter()) {
                    if a == 0.0 { continue; }
                    left.product(a, x).ok_or_else(refused)?;
                    if tail != 0.0 { left.product_remainder(a, tail).ok_or_else(refused)?; }
                }
                integral.add_product(left, value).ok_or_else(refused)?;
            }
        }
        integral.value().ok_or_else(|| refuse("continuous-tone physical power integral overflow"))
    }

    // Positive heat has no cancellation between its squared factors. Preserve
    // its original scalar evaluation and tiny-loss range; only signed work
    // needs retained conjugate factors through the final cancellation.
    fn evaluate_heat(
        &self,
        observation: &Observation,
        effort: &[f64],
        derivative: &[f64],
        input: &[f64],
        coordinates: &mut [f64],
    ) -> Result<f64, PhsError> {
        let l = observation.rows;
        let d = observation.dimension;
        let z = coordinates.get_mut(..d)
            .ok_or_else(|| refuse("continuous-tone heat scratch extent"))?;
        for i in 0..l {
            let mut value = Sum::default();
            for j in 0..self.n {
                value.add(observation.map[i * self.n + j] * effort[j]);
            }
            z[i] = value.value() / self.root_rate;
        }
        z[l..l + self.n].copy_from_slice(derivative);
        for i in 0..self.sources {
            z[l + self.n + i] = self.source_scales[i] * input[i];
        }
        if !finite(z) {
            return Err(refuse("continuous-tone power coordinates overflow"));
        }
        let mut integral = Sum::default();
        for row in observation.left.chunks_exact(d) {
            let mut value = Sum::default();
            for (&a, &x) in row.iter().zip(z.iter()) {
                value.add(a * x);
            }
            let value = value.value();
            integral.add(value * value);
        }
        let value = integral.value();
        if !value.is_finite() {
            return Err(refuse("continuous-tone physical power integral overflow"));
        }
        Ok(value)
    }
}

// Original observation a=C*effort/sqrt(rate), physical derivative b=D*x'/rate.
// Then a'=sqrt(rate)*C*Q*D^-1*b and b'=D*A*D^-1*b+D*G*u'/rate.
// The SAME rotating sources drive motion, work and every original heat channel.
#[allow(clippy::too_many_arguments)]
fn prepare_observation(
    system: &PortHamiltonian,
    split: &LinearSplit,
    frequencies: &[f64],
    source_scales: &[f64],
    scales: &[f64],
    a: &[f64],
    rate: f64,
    dt: f64,
    affine: bool,
    map: Vec<f64>,
    rows: usize,
    paired: bool,
) -> Result<Observation, PhsError> {
    let (n, m) = (system.n, system.m);
    let s = source_scales.len();
    let d = rows + n + s;
    if map.iter().all(|x| *x == 0.0) {
        return Ok(Observation {
            rows,
            dimension: d,
            map,
            left: Vec::new(),
            right: paired.then(Vec::new),
            doublings: 0,
        });
    }
    let root_rate = det::sqrt(rate);
    let mut generator = vec![0.0; d * d];
    let mut factor = vec![0.0; d * d];
    for i in 0..rows {
        factor[i * d + i] = root_rate;
        for j in 0..n {
            let mut value = Contraction::default();
            let refused = || refuse("continuous-tone physical observation map overflows");
            for k in 0..n {
                value
                    .product(map[i * n + k], split.quadratic.matrix[k * n + j])
                    .ok_or_else(refused)?;
            }
            generator[i * d + rows + j] =
                root_rate * value.value().ok_or_else(refused)? / scales[j];
        }
    }
    for i in 0..n {
        for j in 0..n {
            generator[(rows + i) * d + rows + j] = a[i * n + j];
        }
    }
    for p in 0..m {
        for (tone, &omega) in frequencies.iter().enumerate() {
            let stride = if affine { 4 } else { 2 };
            let c = stride * (p * frequencies.len() + tone);
            generator[(rows + n + c) * d + rows + n + c + 1] = omega;
            generator[(rows + n + c + 1) * d + rows + n + c] = -omega;
            if affine {
                generator[(rows + n + c) * d + rows + n + c + 2] = 1.0 / dt;
                generator[(rows + n + c + 1) * d + rows + n + c + 3] = 1.0 / dt;
                generator[(rows + n + c + 2) * d + rows + n + c + 3] = omega;
                generator[(rows + n + c + 3) * d + rows + n + c + 2] = -omega;
            }
            for i in 0..n {
                generator[(rows + i) * d + rows + n + c + 1] =
                    (scales[i] * system.g[i * m + p] / source_scales[c]) * (omega / rate);
                if affine {
                    // b'=A*b+D*G*u'/rate, including the actual envelope slope.
                    generator[(rows + i) * d + rows + n + c + 2] =
                        (scales[i] * system.g[i * m + p] / source_scales[c]) * ((1.0 / dt) / rate);
                }
            }
        }
    }
    let (left, right) = if paired {
        let mut source = vec![0.0; d * d];
        for p in 0..m {
            for tone in 0..frequencies.len() {
                let c = (if affine { 4 } else { 2 }) * (p * frequencies.len() + tone);
                source[p * d + rows + n + c] = 1.0 / source_scales[c];
            }
        }
        let (left, right) =
            positive_heat::integral_cross_root(&generator, &factor, &source, d, dt)?;
        (left, Some(right))
    } else {
        (
            positive_heat::integral_root(&generator, &factor, d, dt)?,
            None,
        )
    };
    Ok(Observation {
        rows,
        dimension: d,
        map,
        left,
        right,
        doublings: positive_heat::scaling_doublings(&generator, d, dt)?,
    })
}

#[derive(Clone, Copy, Default)]
pub(super) struct Sum {
    sum: f64,
    correction: f64,
}
impl Sum {
    pub(super) fn add(&mut self, value: f64) {
        let next = self.sum + value;
        self.correction += if self.sum.abs() >= value.abs() {
            (self.sum - next) + value
        } else {
            (value - next) + self.sum
        };
        self.sum = next;
    }
    pub(super) fn value(self) -> f64 {
        self.sum + self.correction
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{QuadraticStorage, linear_exact::{LinearSplitBudget, split_quadratic}};

    #[test]
    fn complete_tone_force_keeps_motion_after_internal_and_external_cancellation() {
        let large = 2.0_f64.powi(54);
        let delta = 2.0_f64.powi(-27);
        for (g, cosine, effort, expected) in [
            (vec![1.0, 1.0], vec![large, 1.0], large, 1.0),
            (vec![1.0, 1.0], vec![1.0, large], large, 1.0),
            (vec![1.0 + delta], vec![1.0 - delta], 1.0, -delta * delta),
        ] {
            let m = g.len();
            let system = PortHamiltonian::new(
                1, m, vec![0.0], vec![1.0], g,
                Box::new(QuadraticStorage::new(vec![1.0], 1).unwrap()),
            ).unwrap();
            let split = split_quadratic(&system, LinearSplitBudget {
                max_modes: 1, max_omega_dt: 0.1,
            }).unwrap();
            let power = TonePower::new(
                &system, &split, &[0.0], &vec![1.0; 2 * m], 0.125,
                false, &[], &mut || false,
            ).unwrap();
            let mut input = vec![0.0; 2 * m];
            for (p, value) in cosine.into_iter().enumerate() { input[2 * p] = value; }
            let mut derivative = [f64::NAN];
            let mut coordinates = vec![0.0; power.coordinate_len()];
            let mut heat = [0.0];
            power.evaluate_into(
                &system, &[effort], &input, &mut derivative, &mut coordinates, &mut heat,
            ).unwrap();
            // H=x²/2, J=0, R=1: the independently known physical derivative
            // is sum_p(G_p*u_p)-x, including the product's FMA remainder.
            assert_eq!(derivative[0] * power.rate / power.scales[0], expected);
            assert!(heat[0].is_finite() && heat[0] > 0.0);
        }
    }
}
