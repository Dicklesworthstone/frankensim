//! Joint endpoint and interval-mean transition of the advected OU field.
//!
//! For `A = a - i b`, `l = dt/tau`, `theta = (k.U)dt`, `q = l+i theta`,
//! `F = (1-exp(-q))/q` and `h = (1-exp(-2l))/(2l)`, the mean quadrature in
//! the END phase frame is `B = m A + alpha xi + beta chi`, where
//!
//! ```text
//! m     = exp(i theta) F
//! alpha = (2l/sigma) (conj(F)-h)/q
//! beta² = (2l/|q|²) (1-|F|²/h).
//! ```
//!
//! `xi` is the existing endpoint innovation and `chi` is an independent
//! bridge innovation. These follow by integrating the OU stochastic
//! convolution and conditioning its joint Gaussian law on the endpoint.
//! The bridge is essential: endpoint interpolation is not an OU integral.
//! Small-argument divided differences and variances are evaluated directly
//! as analytic series, without subtracting almost equal covariances.

use super::*;

/// Independent bridge stream; the existing endpoint stream is unchanged.
pub const FIXED_STEP_MANN_AVERAGE_KERNEL: u32 = 0x4d46_4156;

type Complex = [f64; 2];

fn multiply(a: Complex, b: Complex) -> Complex {
    [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]]
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct Coefficients {
    mean: Complex,
    endpoint: Complex,
    bridge: f64,
}

impl Coefficients {
    fn prepare(l: f64, theta: f64, sigma: f64) -> Result<Self, Refusal> {
        if !l.is_finite() || l <= 0.0 || !theta.is_finite() {
            return Err(error(
                "mann-step-average-coefficients",
                "interval decay and advective phase must be finite and representable",
            ));
        }
        let h = (0.5 * -det::expm1(-2.0 * l)) / l;
        let scale = l.max(theta.abs());
        let (f, endpoint, bridge) = if scale <= 0.5 {
            // Twenty terms suffice on |q| <= sqrt(1/2). For the divided
            // difference below |b| <= 1; its omitted terms are < 1e-20.
            const TERMS: usize = 20;
            let minus_q = [-l, -theta];
            let mut f = [1.0, 0.0];
            let mut term = [1.0, 0.0];
            for n in 1..=TERMS {
                term = multiply(term, minus_q).map(|x| x / (n + 1) as f64);
                for i in 0..2 {
                    f[i] += term[i];
                }
            }

            // (conj(F)-h)/q = sum_{n>=1} (a^n-b^n)/(a-b)/(n+1)!,
            // with a=-conj(q), b=-2l and a-b=q. Evaluate the polynomial
            // divided difference, including its q -> 0 value 1/2.
            let a = [-l, theta];
            let b = -2.0 * l;
            let mut power_difference = [1.0, 0.0];
            let mut b_power = 1.0;
            let mut inverse_factorial = 0.5;
            let mut divided = [0.0; 2];
            for n in 1..=TERMS {
                for i in 0..2 {
                    divided[i] += inverse_factorial * power_difference[i];
                }
                b_power *= b;
                power_difference = multiply(a, power_difference);
                power_difference[0] += b_power;
                inverse_factorial /= (n + 2) as f64;
            }
            let endpoint = divided.map(|x| (2.0 * l / sigma) * x);

            // V=(h-|F|²)/|q|² = integral |(exp(-qu)-F)/q|² du.
            // Expanding the integral directly gives d_n=(-q)^(n-1)/n!
            // and weight nm/[(n+m+1)(n+1)(m+1)], with leading V=1/12.
            let mut powers = [[0.0; 2]; TERMS];
            powers[0] = [1.0, 0.0];
            for n in 1..TERMS {
                powers[n] = multiply(powers[n - 1], minus_q).map(|x| x / (n + 1) as f64);
            }
            let mut variance = 0.0;
            let mut correction = 0.0;
            for (i, a) in powers.iter().enumerate() {
                for (j, b) in powers.iter().enumerate() {
                    let n = (i + 1) as f64;
                    let m = (j + 1) as f64;
                    let weight = n * m / ((n + m + 1.0) * (n + 1.0) * (m + 1.0));
                    let contribution = weight * (a[0] * b[0] + a[1] * b[1]);
                    let adjusted = contribution - correction;
                    let next = variance + adjusted;
                    correction = (next - variance) - adjusted;
                    variance = next;
                }
            }
            let bridge = core::f64::consts::SQRT_2 * det::sqrt(l) * det::sqrt(variance / h);
            (f, endpoint, bridge)
        } else {
            // Scale before complex division: neither |q|² nor 2l need be
            // representable. expm1 and 2 sin²(theta/2) preserve a small real
            // numerator even near integer periods with slow OU evolution.
            let decay = det::exp(-l);
            let (sin, _) = det::sin_cos(theta);
            let half_sin = det::sin(0.5 * theta);
            let numerator = [
                -det::expm1(-l) + 2.0 * decay * half_sin * half_sin,
                decay * sin,
            ];
            let normalized = [l / scale, theta / scale];
            let norm_squared = normalized[0] * normalized[0] + normalized[1] * normalized[1];
            let conjugate = [normalized[0], -normalized[1]];
            let f = multiply(numerator, conjugate).map(|x| (x / norm_squared) / scale);
            // Preserve the original endpoint sigma, including its rounding.
            // Splitting 2l into sqrt factors avoids an overflowing numerator
            // or prematurely underflowing l/scale at very large theta.
            let root_l = det::sqrt(l);
            let endpoint_scale = (2.0 * root_l / scale) * (root_l / sigma);
            let endpoint =
                multiply([f[0] - h, -f[1]], conjugate).map(|x| endpoint_scale * (x / norm_squared));
            let ratio = det::hypot(f[0], f[1]) / det::sqrt(h);
            let gap = 1.0 - ratio * ratio;
            let bridge = ((core::f64::consts::SQRT_2 * root_l / scale) / det::sqrt(norm_squared))
                * det::sqrt(gap);
            (f, endpoint, bridge)
        };
        let (sin, cos) = det::sin_cos(theta);
        let result = Self {
            mean: multiply([cos, sin], f),
            endpoint,
            bridge,
        };
        if !h.is_finite()
            || h <= 0.0
            || result
                .mean
                .iter()
                .chain(&result.endpoint)
                .any(|x| !x.is_finite())
            || !result.bridge.is_finite()
            || result.bridge <= 0.0
        {
            return Err(error(
                "mann-step-average-coefficients",
                "joint OU interval-mean coefficients are not representable",
            ));
        }
        Ok(result)
    }
}

/// Reusable, exact Gaussian average of one native physical interval.
///
/// Prepare once, then pass this workspace to
/// [`FixedStepMannField::advance_one_with_average`]. Every observation uses
/// the same accepted interval quadratures, including spatial filters, wall
/// blocking and analytic gradients. Sampling neither allocates nor draws.
/// The average is unavailable until the first successful joint advance.
///
/// This averages the modeled **velocity** and its spatial derivatives. It
/// does not average nonlinear momentum/energy fluxes, provide a physical
/// subgrid closure, or make a continuous OU process band-limited. Existing
/// endpoint checkpoints completely determine all future interval averages;
/// this observation workspace does not change their format or stream bits.
#[derive(Clone, Debug)]
pub struct MannStepAverage {
    definition: Arc<Definition>,
    coefficients: Vec<Coefficients>,
    end_step: Option<u64>,
    amplitudes: [Vec<f64>; 4],
    trial: [Vec<f64>; 4],
}

impl PartialEq for MannStepAverage {
    fn eq(&self, other: &Self) -> bool {
        self.definition == other.definition
            && self.coefficients == other.coefficients
            && self.end_step == other.end_step
            && self.amplitudes == other.amplitudes
    }
}

impl MannStepAverage {
    /// Admit coefficients and both amplitude banks before any temporal work.
    /// An ordinary field can admit an infinitely stiff rounded dt/tau; the
    /// exact mean additionally requires a finite representable dt/tau and
    /// advective phase per mode, and positive representable bridge scale.
    pub fn prepare(
        source: &FixedStepMannField,
        maximum_workspace_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Refusal> {
        poll(&mut cancelled)?;
        let n = source.model().mode_count();
        if Self::required_workspace_bytes(n)? > maximum_workspace_bytes {
            return Err(error(
                "mann-step-average-budget",
                "interval-average workspace exceeds its byte budget",
            ));
        }
        let mut coefficients = Vec::new();
        coefficients.try_reserve_exact(n).map_err(|_| {
            error(
                "mann-step-average-allocation",
                "unable to reserve interval-average coefficients",
            )
        })?;
        for (i, mode) in source.model().modes.iter().enumerate() {
            poll(&mut cancelled)?;
            coefficients.push(Coefficients::prepare(
                source.dt_s() / source.definition.times[i],
                dot(mode.k, source.advection_mps()) * source.dt_s(),
                source.definition.transition[i][1],
            )?);
        }
        let result = Self {
            definition: Arc::clone(&source.definition),
            coefficients,
            end_step: None,
            amplitudes: [zeros(n)?, zeros(n)?, zeros(n)?, zeros(n)?],
            trial: [zeros(n)?, zeros(n)?, zeros(n)?, zeros(n)?],
        };
        if result.workspace_bytes() > maximum_workspace_bytes {
            return Err(error(
                "mann-step-average-budget",
                "allocated interval-average capacities exceed the byte budget",
            ));
        }
        poll(&mut cancelled)?;
        Ok(result)
    }

    /// Additional numerical storage: 104 bytes/mode, excluding the shared
    /// source definition, fixed metadata and allocator metadata.
    pub fn required_workspace_bytes(modes: usize) -> Result<usize, Refusal> {
        modes
            .checked_mul(13 * core::mem::size_of::<f64>())
            .ok_or_else(|| {
                error(
                    "mann-step-average-budget",
                    "interval-average workspace size overflows",
                )
            })
    }

    pub fn workspace_bytes(&self) -> usize {
        self.coefficients.capacity() * core::mem::size_of::<Coefficients>()
            + self
                .amplitudes
                .iter()
                .chain(&self.trial)
                .map(Vec::capacity)
                .sum::<usize>()
                * core::mem::size_of::<f64>()
    }

    pub fn model(&self) -> &MannSpectralModel {
        &self.definition.model
    }

    pub fn dt_s(&self) -> f64 {
        self.definition.dt_s
    }

    /// Native endpoint index of the last accepted average.
    pub const fn end_step(&self) -> Option<u64> {
        self.end_step
    }

    /// Physical phase-time endpoints of the last accepted average. The mean
    /// includes the continuously advected OU field throughout this interval.
    pub fn interval(&self) -> Option<[f64; 2]> {
        self.end_step.map(|end| {
            [
                self.definition.origin_time_s + (end - 1) as f64 * self.dt_s(),
                self.definition.origin_time_s + end as f64 * self.dt_s(),
            ]
        })
    }

    /// Check immutable spectrum, clock, phase origin and random-stream
    /// identity. Reusing scratch after a failed consumer trial is allowed;
    /// this workspace need not hold the source's preceding interval.
    pub fn validate_source(&self, source: &FixedStepMannField) -> Result<(), Refusal> {
        if !Arc::ptr_eq(&self.definition, &source.definition)
            && self.definition != source.definition
        {
            return Err(error(
                "mann-step-average-source-mismatch",
                "interval-average workspace belongs to a different model or clock",
            ));
        }
        Ok(())
    }

    fn time_s(&self) -> Result<f64, Refusal> {
        self.interval().map(|interval| interval[1]).ok_or_else(|| {
            error(
                "mann-step-average-unavailable",
                "no interval average has been accepted",
            )
        })
    }

    pub fn sample(&self, point_m: [f64; 3]) -> Result<FlowSample, Refusal> {
        self.sample_with_gains(point_m, None)
    }

    pub fn sample_filtered(
        &self,
        point_m: [f64; 3],
        filter: &MannBoxFilter,
    ) -> Result<FlowSample, Refusal> {
        self.sample_with_gains(point_m, filter.gains(self.model())?)
    }

    fn sample_with_gains(
        &self,
        point_m: [f64; 3],
        gains: Option<&[f64]>,
    ) -> Result<FlowSample, Refusal> {
        sample_modes(
            self.model(),
            self.time_s()?,
            core::array::from_fn(|i| self.amplitudes[i].as_slice()),
            point_m,
            gains,
        )
    }

    pub fn sample_wall_blocked(
        &self,
        point_m: [f64; 3],
        wall: MannPlaneWall,
    ) -> Result<FlowSample, Refusal> {
        wall::sample_modes(
            self.model(),
            self.time_s()?,
            core::array::from_fn(|i| self.amplitudes[i].as_slice()),
            point_m,
            wall,
        )
    }

    pub fn sample_wall_blocked_filtered(
        &self,
        point_m: [f64; 3],
        filter: &MannWallBoxFilter,
    ) -> Result<FlowSample, Refusal> {
        wall::filter::sample_modes(
            self.model(),
            self.time_s()?,
            core::array::from_fn(|i| self.amplitudes[i].as_slice()),
            point_m,
            filter,
        )
    }
}

fn innovation(seed: u64, kernel: u32, path: usize, mode: usize, step: u64) -> f64 {
    let key = StreamKey {
        seed: seed ^ PATH_KEYS[path],
        kernel,
        tile: mode as u32,
    };
    Stream::resume(StreamCheckpoint::current(key, 2 * step))
        .expect("current-version checkpoint and admitted counter always resume")
        .next_normal()
}

impl FixedStepMannField {
    /// Jointly publish the next endpoint and exact interval-mean velocity.
    /// The budget counts eight normal draws/mode: four unchanged endpoint
    /// draws and four independent bridge draws. No numerical allocation is
    /// performed. Cancellation is polled every eight modes (64 draws) and
    /// immediately before publication; either both accepted states advance
    /// or neither does. Endpoint checkpoint bits equal ordinary advance.
    pub fn advance_one_with_average(
        &mut self,
        average: &mut MannStepAverage,
        maximum_normal_draws: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), Refusal> {
        poll(&mut cancelled)?;
        average.validate_source(self)?;
        let step = self
            .step
            .checked_add(1)
            .ok_or_else(|| error("mann-fixed-step-clock", "fixed-step index overflows"))?;
        time_at(self.origin_time_s(), self.dt_s(), step)?;
        let n = self.model().mode_count();
        let work = n
            .checked_mul(8)
            .ok_or_else(|| error("mann-fixed-step-work", "interval-average work overflows"))?;
        if work > maximum_normal_draws {
            return Err(error(
                "mann-fixed-step-work",
                "interval-average work exceeds the declared budget",
            ));
        }
        let seed = self.definition.evolution_seed;
        for mode in 0..n {
            if mode % 8 == 0 {
                poll(&mut cancelled)?;
            }
            let xi: [f64; 4] = core::array::from_fn(|path| {
                innovation(seed, FIXED_STEP_MANN_KERNEL, path, mode, step)
            });
            let chi: [f64; 4] = core::array::from_fn(|path| {
                innovation(seed, FIXED_STEP_MANN_AVERAGE_KERNEL, path, mode, step)
            });
            let [rho, sigma] = self.definition.transition[mode];
            for (path, normal) in xi.iter().enumerate() {
                self.trial[path][mode] = rho * self.amplitudes[path][mode] + sigma * normal;
            }
            let c = average.coefficients[mode];
            for polarization in 0..2 {
                let a = 2 * polarization;
                let b = a + 1;
                let mean = multiply(
                    c.mean,
                    [self.amplitudes[a][mode], -self.amplitudes[b][mode]],
                );
                let endpoint = multiply(c.endpoint, [xi[a], -xi[b]]);
                average.trial[a][mode] = mean[0] + endpoint[0] + c.bridge * chi[a];
                average.trial[b][mode] = -mean[1] - endpoint[1] + c.bridge * chi[b];
            }
            if self
                .trial
                .iter()
                .chain(&average.trial)
                .any(|path| !path[mode].is_finite())
            {
                return Err(error(
                    "mann-step-average-state",
                    "joint endpoint or interval-average amplitude is not finite",
                ));
            }
        }
        poll(&mut cancelled)?;
        core::mem::swap(&mut self.amplitudes, &mut self.trial);
        core::mem::swap(&mut average.amplitudes, &mut average.trial);
        self.step = step;
        average.end_step = Some(step);
        Ok(())
    }
}

#[cfg(test)]
mod tests;
