//! Original port work and resistance squares along continuous source motion.
use super::super::{LinearSplit, finite, positive_heat, refuse};
use crate::{PhsError, PortHamiltonian, ResistanceWitness};
use fs_math::det;

struct Observation {
    rows: usize,
    dimension: usize,
    map: Vec<f64>,
    left: Vec<f64>,
    right: Option<Vec<f64>>,
}

pub(super) struct TonePower {
    n: usize,
    m: usize,
    sources: usize,
    stride: usize,
    rate: f64,
    root_rate: f64,
    scales: Vec<f64>,
    source_scales: Vec<f64>,
    work: Observation,
    heat: Vec<Observation>,
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
    pub(super) fn new(
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
                a[i * n + j] = (0..n)
                    .map(|k| {
                        (system.j[i * n + k] - system.r[i * n + k])
                            * split.quadratic.matrix[k * n + j]
                    })
                    .sum::<f64>()
                    * scales[i]
                    / scales[j];
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
        Ok(Self {
            n,
            m,
            sources,
            stride,
            rate,
            root_rate: det::sqrt(rate),
            scales,
            source_scales: source_scales.to_vec(),
            work,
            heat,
        })
    }

    /// Uses the original physical effort and force imbalance before any
    /// factor contraction. This keeps a stationary resistance nullspace and
    /// tiny forced motion from becoming subtraction of large heat terms.
    pub(super) fn evaluate_into(
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
            let mut force = Sum::default();
            let mut conservative = Sum::default();
            let mut resistive = Sum::default();
            for p in 0..self.m {
                let mut total = Sum::default();
                for t in 0..tones {
                    total.add(input[self.stride * (p * tones + t)]);
                }
                force.add(system.g[i * self.m + p] * total.value());
            }
            for (j, &e) in effort.iter().enumerate() {
                conservative.add(system.j[i * self.n + j] * e);
                resistive.add(system.r[i * self.n + j] * e);
            }
            derivative[i] = ((force.value() + conservative.value()) - resistive.value())
                * self.scales[i]
                / self.rate;
        }
        let supplied = self.evaluate(&self.work, effort, derivative, input, coordinates)?;
        for (value, observation) in heat.iter_mut().zip(&self.heat) {
            *value = self.evaluate(observation, effort, derivative, input, coordinates)?;
        }
        Ok(supplied)
    }

    fn evaluate(
        &self,
        observation: &Observation,
        effort: &[f64],
        derivative: &[f64],
        input: &[f64],
        coordinates: &mut [f64],
    ) -> Result<f64, PhsError> {
        if observation.left.is_empty() {
            return Ok(0.0);
        }
        let l = observation.rows;
        let d = observation.dimension;
        let z = &mut coordinates[..d];
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
        for (i, row) in observation.left.chunks_exact(d).enumerate() {
            let mut left = Sum::default();
            for (&a, &x) in row.iter().zip(z.iter()) {
                left.add(a * x);
            }
            let left = left.value();
            if let Some(right) = &observation.right {
                let mut value = Sum::default();
                for (&a, &x) in right[i * d..(i + 1) * d].iter().zip(z.iter()) {
                    value.add(a * x);
                }
                integral.add(left * value.value());
            } else {
                integral.add(left * left);
            }
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
        });
    }
    let root_rate = det::sqrt(rate);
    let mut generator = vec![0.0; d * d];
    let mut factor = vec![0.0; d * d];
    for i in 0..rows {
        factor[i * d + i] = root_rate;
        for j in 0..n {
            generator[i * d + rows + j] = root_rate
                * (0..n)
                    .map(|k| map[i * n + k] * split.quadratic.matrix[k * n + j])
                    .sum::<f64>()
                / scales[j];
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
