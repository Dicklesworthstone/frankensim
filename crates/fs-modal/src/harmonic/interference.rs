//! Finite-window mechanical interference through the original coupled pencil.
//!
//! The real motion is the sum of particular harmonic responses, not a sum of
//! their independently cycle-averaged powers. Applying the original M, K and
//! symmetric part of C before observation keeps all spatial and temporal cross
//! terms in kinetic energy, elastic energy and viscous heat. The instantaneous
//! identity is F^T v = v^T M a + v^T K q + v^T D v + flow^T effort,
//! D=(C+C^T)/2; the skew gyroscopic part does no work. A supplied physical
//! impedance contributes signed boundary work. Every term is integrated independently.

use super::{
    HarmonicPortControls, HarmonicResponse, HarmonicResponseControls, HarmonicResponseError,
    PreparedHarmonicResponse, admit_size, finite, invalid, poll, symmetric_entry,
};
use fs_exec::Cx;
use fs_fft::harmonic_window::{HarmonicProduct, HarmonicWindow, WindowError};
use fs_math::{c64::C64, dd::Dd};

mod ensemble;
pub use ensemble::{
    HarmonicEnsembleControls, HarmonicFactorStatistics, HarmonicForceFactor,
    MAX_ENSEMBLE_FORCE_TONES, MAX_INTERFERENCE_FACTORS, PreparedHarmonicInterferenceEnsemble,
};

/// Hard ceiling on authored forcing tones. Exact duplicates are combined;
/// nearby frequencies are never rounded together or discarded.
pub const MAX_INTERFERENCE_TONES: usize = 128;

/// One coherent force tone on the common source clock, under exp(-i omega t).
/// Forces are peak phasors in the original pencil's generalized coordinates.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicForceTone<'a> {
    /// Nonnegative angular frequency [rad/s], including DC.
    pub angular_frequency_rad_s: f64,
    /// One finite peak generalized force per original pencil coordinate.
    pub force: &'a [C64],
}

/// One physical owner's complete impedance at an authored frequency.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicImpedanceSample<'a> {
    /// Nonnegative angular frequency [rad/s], including DC.
    pub angular_frequency_rad_s: f64,
    /// Row-major `ports x ports` effort/flow map under exp(-i omega t).
    pub impedance: &'a [C64],
}

/// Physical load shared by all coherent tones and statistical factors.
///
/// Samples must be strictly increasing and cover exactly the distinct source
/// frequencies. No interpolation, extrapolation or impedance fit is performed.
/// Each Hermitian part passes the existing pointwise PSD admission; the source
/// owner remains responsible for broadband causality and physical port units.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicInterferenceLoad<'a> {
    /// Real row-major `n x ports` map: reaction is `-B effort`, flow is `B^T v`.
    pub force_map: &'a [f64],
    /// Number of physical effort/flow ports.
    pub ports: usize,
    /// One full impedance per distinct authored source frequency.
    pub samples: &'a [HarmonicImpedanceSample<'a>],
    /// Existing physical-port preparation limits, applied at every frequency.
    pub controls: HarmonicPortControls,
}

impl<'a> HarmonicInterferenceLoad<'a> {
    fn validate(
        self,
        n: usize,
        maximum_samples: usize,
        cx: &Cx<'_>,
    ) -> Result<(), HarmonicResponseError> {
        cap("load ports", self.ports, self.controls.maximum_ports)?;
        let np = n
            .checked_mul(self.ports)
            .ok_or(invalid("force map extent overflow"))?;
        let pp = self
            .ports
            .checked_mul(self.ports)
            .ok_or(invalid("impedance extent overflow"))?;
        if self.ports == 0
            || self.force_map.len() != np
            || self.force_map.iter().any(|x| !x.is_finite())
            || self.samples.is_empty()
            || self.samples.len() > maximum_samples
        {
            return Err(invalid(
                "complete physical map and exact source-frequency load samples required",
            ));
        }
        for (i, sample) in self.samples.iter().enumerate() {
            poll(cx)?;
            let omega = sample.angular_frequency_rad_s;
            if !omega.is_finite()
                || omega < 0.0
                || (i > 0 && omega <= self.samples[i - 1].angular_frequency_rad_s)
                || sample.impedance.len() != pp
                || sample.impedance.iter().any(|&z| !finite(z))
            {
                return Err(invalid(
                    "finite, strictly ordered, complete impedance samples required",
                ));
            }
        }
        Ok(())
    }

    fn impedance(self, omega: f64) -> Result<&'a [C64], HarmonicResponseError> {
        self.samples
            .binary_search_by(|sample| {
                // Positive and negative zero represent the same physical DC.
                if sample.angular_frequency_rad_s == omega {
                    core::cmp::Ordering::Equal
                } else {
                    sample.angular_frequency_rad_s.total_cmp(&omega)
                }
            })
            .map(|index| self.samples[index].impedance)
            .map_err(|_| invalid("every source frequency requires its exact impedance sample"))
    }
}

/// Admission for the complete multi-frequency preparation and observation.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicInterferenceControls {
    /// Maximum authored tones, including duplicates, bounded by 128.
    pub maximum_tones: usize,
    /// Conservative peak numerical payload, including one original dense
    /// response preparation/solve, retained responses and coefficient actions.
    /// Caller input, allocator overhead and thread stacks are excluded.
    pub maximum_payload_bytes: usize,
    /// Maximum real products: `(6*n + loaded_ports)*distinct_tones^2`.
    /// Checked before solving. This is a loop count, not a wall-time guarantee.
    pub maximum_product_terms: usize,
}

/// Sensitivity-weighted physical moments of the complete real motion.
///
/// The window is normalized to unit integral. Multiplying powers by its
/// duration gives the weighted work/heat; only a uniform window gives the
/// unweighted physical integral. Finite windows retain sum-frequency terms
/// as well as slow beats. No cycle averaging or source independence is assumed.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HarmonicInterferenceMoments {
    /// <v^T M v>/2 [J].
    pub mean_kinetic_energy_j: f64,
    /// <q^T K q>/2 [J].
    pub mean_elastic_energy_j: f64,
    /// <F^T v> [W]; can be negative while the structure returns stored energy.
    pub mean_input_power_w: f64,
    /// <v^T D v> [W], with full symmetric D=(C+C^T)/2.
    /// Signed floating roundoff is retained rather than clipped.
    pub mean_viscous_power_w: f64,
    /// <flow^T effort> [W] at the supplied physical load, zero when unloaded.
    /// This is signed net boundary work, including reactive energy returned to
    /// the structure; it is not an inferred heat or radiation decomposition.
    pub mean_boundary_power_w: f64,
    /// <v^T M a + v^T K q> [W], computed independently of force and heat.
    /// Equals endpoint stored-energy change divided by duration only for a
    /// uniform window. A shaped detector measures the weighted local rate.
    pub mean_stored_energy_rate_w: f64,
    /// Input minus viscous power, boundary work and structural storage rate.
    pub power_balance_residual_w: f64,
}

#[derive(Clone, Debug)]
struct Tone {
    omega: f64,
    force: Vec<C64>,
    response: HarmonicResponse,
    mass_velocity: Vec<[Dd; 2]>,
    stiffness_displacement: Vec<[Dd; 2]>,
    viscous_velocity: Vec<[Dd; 2]>,
}

/// Prepared coherent interference in one complete real M/C/K pencil.
///
/// Each distinct frequency is solved through [`PreparedHarmonicResponse`];
/// exact duplicate forces are summed in double-double precision BEFORE their
/// solve. Full mutual inertia, stiffness, nonproportional viscous damping and
/// signed gyroscopic terms remain intact. Preparations retain original-coordinate
/// coefficient actions. Observation allocates nothing, performs no linear solve
/// or matrix application, and reuses one analytic kernel per frequency pair.
///
/// This is a finite-dimensional, particular harmonic motion on a shared clock.
/// It does not add startup/free-decay transients, infer stochastic coherence,
/// resolve nonlinear interactions, or infer a broadband boundary impedance.
/// Singular frequencies still refuse. Solver residuals and DD accumulation
/// are numerical evidence, not forward-error or near-null relative bounds.
#[derive(Clone, Debug)]
pub struct PreparedHarmonicInterference {
    n: usize,
    tones: Vec<Tone>,
    maximum_relative_response_residual: f64,
}

impl PreparedHarmonicInterference {
    /// Solve and retain every coherent tone in one original physical pencil.
    /// Frequencies are finite and nonnegative; DC and exact duplicates work.
    /// No force, source matrix or accepted simulation state is mutated.
    /// Numerical and count budgets are admitted before coefficient allocation.
    #[allow(clippy::too_many_arguments)] // Complete pencil plus independent physical budgets.
    pub fn try_new(
        mass: &[f64],
        damping: &[f64],
        stiffness: &[f64],
        n: usize,
        forces: &[HarmonicForceTone<'_>],
        response_controls: HarmonicResponseControls,
        controls: HarmonicInterferenceControls,
        cx: &Cx<'_>,
    ) -> Result<Self, HarmonicResponseError> {
        Self::prepare(
            mass,
            damping,
            stiffness,
            n,
            forces,
            None,
            response_controls,
            controls,
            cx,
        )
    }

    /// Retain loaded coherent responses and their signed finite-window boundary
    /// work. Every exact source frequency uses the original loaded pencil;
    /// mutual and reactive impedance entries survive without cycle averaging.
    #[allow(clippy::too_many_arguments)]
    pub fn try_with_port_impedance(
        mass: &[f64],
        damping: &[f64],
        stiffness: &[f64],
        n: usize,
        forces: &[HarmonicForceTone<'_>],
        load: HarmonicInterferenceLoad<'_>,
        response_controls: HarmonicResponseControls,
        controls: HarmonicInterferenceControls,
        cx: &Cx<'_>,
    ) -> Result<Self, HarmonicResponseError> {
        Self::prepare(
            mass,
            damping,
            stiffness,
            n,
            forces,
            Some(load),
            response_controls,
            controls,
            cx,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn prepare(
        mass: &[f64],
        damping: &[f64],
        stiffness: &[f64],
        n: usize,
        forces: &[HarmonicForceTone<'_>],
        load: Option<HarmonicInterferenceLoad<'_>>,
        response_controls: HarmonicResponseControls,
        controls: HarmonicInterferenceControls,
        cx: &Cx<'_>,
    ) -> Result<Self, HarmonicResponseError> {
        poll(cx)?;
        admit_size(n, response_controls)?;
        let count = forces.len();
        if count == 0 {
            return Err(invalid("at least one coherent force tone is required"));
        }
        cap(
            "interference tones",
            count,
            controls.maximum_tones.min(MAX_INTERFERENCE_TONES),
        )?;
        let extent = n.checked_mul(n).ok_or(invalid("pencil extent overflow"))?;
        if [mass, damping, stiffness].iter().any(|x| x.len() != extent) {
            return Err(invalid("one complete M/C/K pencil is required"));
        }
        for tone in forces {
            poll(cx)?;
            if !tone.angular_frequency_rad_s.is_finite()
                || tone.angular_frequency_rad_s < 0.0
                || tone.force.len() != n
                || tone.force.iter().any(|&z| !finite(z))
            {
                return Err(invalid(
                    "finite nonnegative tone frequencies and complete forces required",
                ));
            }
        }
        // Bounded stack indices: admission itself allocates no numerical arrays.
        let mut order = [0usize; MAX_INTERFERENCE_TONES];
        for (i, slot) in order[..count].iter_mut().enumerate() {
            *slot = i;
        }
        order[..count].sort_unstable_by(|&a, &b| {
            forces[a]
                .angular_frequency_rad_s
                .total_cmp(&forces[b].angular_frequency_rad_s)
                .then_with(|| a.cmp(&b))
        });
        let distinct = 1 + order[..count]
            .windows(2)
            .filter(|pair| {
                forces[pair[0]].angular_frequency_rad_s != forces[pair[1]].angular_frequency_rad_s
            })
            .count();
        if let Some(load) = load {
            load.validate(n, distinct, cx)?;
            if load.samples.len() != distinct {
                return Err(invalid(
                    "exactly one load sample per distinct source frequency required",
                ));
            }
            for tone in forces {
                load.impedance(tone.angular_frequency_rad_s)?;
            }
        }
        admit_payload(n, distinct, load, controls)?;
        let products = distinct
            .checked_mul(distinct)
            .and_then(|x| {
                x.checked_mul(n.checked_mul(6)?.checked_add(load.map_or(0, |l| l.ports))?)
            })
            .ok_or(invalid("interference product extent overflow"))?;
        cap(
            "interference product terms",
            products,
            controls.maximum_product_terms,
        )?;

        let mut tones = Vec::with_capacity(distinct);
        let mut maximum_relative_response_residual = 0.0_f64;
        let mut first = 0;
        while first < count {
            poll(cx)?;
            let omega = forces[order[first]].angular_frequency_rad_s;
            let mut end = first + 1;
            while end < count && forces[order[end]].angular_frequency_rad_s == omega {
                end += 1;
            }
            let mut force = Vec::with_capacity(n);
            for i in 0..n {
                if i % 64 == 0 {
                    poll(cx)?;
                }
                let mut sum = [Dd::ZERO; 2];
                for &source in &order[first..end] {
                    let value = parts(forces[source].force[i]);
                    sum[0] = checked(sum[0] + value[0])?;
                    sum[1] = checked(sum[1] + value[1])?;
                }
                force.push(C64::new(sum[0].to_f64(), sum[1].to_f64()));
            }
            let prepared = prepare_response(
                mass,
                damping,
                stiffness,
                n,
                omega,
                load,
                response_controls,
                cx,
            )?;
            let response = prepared.solve(&force, cx)?;
            maximum_relative_response_residual =
                maximum_relative_response_residual.max(response.relative_residual);
            let mass_velocity = apply(mass, &response.velocity, n, false, cx)?;
            let stiffness_displacement = apply(stiffness, &response.displacement, n, false, cx)?;
            let viscous_velocity = apply(damping, &response.velocity, n, true, cx)?;
            tones.push(Tone {
                omega,
                force,
                response,
                mass_velocity,
                stiffness_displacement,
                viscous_velocity,
            });
            first = end;
        }
        poll(cx)?;
        Ok(Self {
            n,
            tones,
            maximum_relative_response_residual,
        })
    }

    /// Number of distinct frequencies; exact duplicate sources share one solve.
    #[must_use]
    pub fn tone_count(&self) -> usize {
        self.tones.len()
    }

    /// Maximum original dynamic-equation backward residual over all tone solves.
    #[must_use]
    pub const fn maximum_relative_response_residual(&self) -> f64 {
        self.maximum_relative_response_residual
    }

    /// Retained single-tone response and its original angular frequency.
    /// Frequencies are ordered increasingly; this borrows immutable coefficients.
    #[must_use]
    pub fn tone_response(&self, index: usize) -> Option<(f64, &HarmonicResponse)> {
        self.tones
            .get(index)
            .map(|tone| (tone.omega, &tone.response))
    }

    /// Integrate the coherent real motion over the actual detector window.
    /// One shared fs-fft kernel supplies sum/difference frequencies and DD phase
    /// reduction per pair. No tone power or phase is separately averaged away.
    /// Polls cancellation per pair and at 64-coordinate boundaries; no partial
    /// result is published. Work and payload were admitted at preparation.
    pub fn observe(
        &self,
        window: HarmonicWindow,
        cx: &Cx<'_>,
    ) -> Result<HarmonicInterferenceMoments, HarmonicResponseError> {
        let sums = observe_sums(
            self.n,
            &self.tones,
            window,
            HarmonicFactorStatistics::Real,
            cx,
        )?;
        let result = finish_moments(sums)?;
        poll(cx)?;
        Ok(result)
    }
}

// The deterministic owner and mixed-factor owner use exactly the same seven
// physical contractions. Each frequency pair prepares its temporal kernel once;
// dispatch on the explicit law occurs outside the coordinate/product loop.
fn observe_sums(
    n: usize,
    tones: &[Tone],
    window: HarmonicWindow,
    statistics: HarmonicFactorStatistics,
    cx: &Cx<'_>,
) -> Result<[Dd; 7], HarmonicResponseError> {
    poll(cx)?;
    let mut moments = [Dd::ZERO; 7];
    for a in tones {
        for b in tones {
            poll(cx)?;
            match statistics {
                HarmonicFactorStatistics::Real => {
                    let kernel = window
                        .product_kernel(a.omega, b.omega)
                        .map_err(observation)?;
                    accumulate_pair(n, a, b, |x, y| product(kernel, x, y), &mut moments, cx)?;
                }
                HarmonicFactorStatistics::ProperComplex => {
                    let kernel = window
                        .proper_product_kernel(a.omega, b.omega)
                        .map_err(observation)?;
                    accumulate_pair(
                        n,
                        a,
                        b,
                        |x, y| kernel.mean_real_product_parts_dd(x, y).map_err(observation),
                        &mut moments,
                        cx,
                    )?;
                }
            }
        }
    }
    Ok(moments)
}

fn accumulate_pair(
    n: usize,
    a: &Tone,
    b: &Tone,
    product: impl Fn([Dd; 2], [Dd; 2]) -> Result<Dd, HarmonicResponseError>,
    moments: &mut [Dd; 7],
    cx: &Cx<'_>,
) -> Result<(), HarmonicResponseError> {
    for i in 0..n {
        if i % 64 == 0 {
            poll(cx)?;
        }
        let q = parts(a.response.displacement[i]);
        let v = parts(a.response.velocity[i]);
        let mv = b.mass_velocity[i];
        let frequency = Dd::from_f64(b.omega);
        let ma = [checked(mv[1] * frequency)?, checked(-mv[0] * frequency)?];
        let values = [
            product(v, mv)? * Dd::from_f64(0.5),
            product(q, b.stiffness_displacement[i])? * Dd::from_f64(0.5),
            product(parts(a.force[i]), parts(b.response.velocity[i]))?,
            product(v, b.viscous_velocity[i])?,
            product(v, ma)?,
            product(v, b.stiffness_displacement[i])?,
        ];
        for (sum, value) in moments.iter_mut().zip(values) {
            *sum = checked(*sum + value)?;
        }
    }
    if let (Some(a), Some(b)) = (&a.response.port_load, &b.response.port_load) {
        for (i, (&flow, &effort)) in a.flow.iter().zip(&b.effort).enumerate() {
            if i % 64 == 0 {
                poll(cx)?;
            }
            moments[6] = checked(moments[6] + product(parts(flow), parts(effort))?)?;
        }
    }
    Ok(())
}

fn finish_moments(moments: [Dd; 7]) -> Result<HarmonicInterferenceMoments, HarmonicResponseError> {
    let stored = checked(moments[4] + moments[5])?;
    let residual = checked(moments[2] - moments[3] - moments[6] - stored)?;
    Ok(HarmonicInterferenceMoments {
        mean_kinetic_energy_j: moments[0].to_f64(),
        mean_elastic_energy_j: moments[1].to_f64(),
        mean_input_power_w: moments[2].to_f64(),
        mean_viscous_power_w: moments[3].to_f64(),
        mean_boundary_power_w: moments[6].to_f64(),
        mean_stored_energy_rate_w: stored.to_f64(),
        power_balance_residual_w: residual.to_f64(),
    })
}

fn admit_payload(
    n: usize,
    tones: usize,
    load: Option<HarmonicInterferenceLoad<'_>>,
    controls: HarmonicInterferenceControls,
) -> Result<(), HarmonicResponseError> {
    // Original dense owner bound plus 144*n retained scalar bytes per tone
    // (F,q,v and three DD coefficient actions), headers and bounded scratch.
    let additional_load = load_payload(n, tones, load)?;
    let bytes = n
        .checked_mul(n)
        .and_then(|x| x.checked_mul(144))
        .and_then(|x| x.checked_add(n.checked_mul(192)?))
        .and_then(|x| {
            x.checked_add(
                tones.checked_mul(
                    n.checked_mul(144)?
                        .checked_add(core::mem::size_of::<Tone>())?,
                )?,
            )
        })
        .and_then(|x| x.checked_add(additional_load))
        .ok_or(invalid("interference payload extent overflow"))?;
    cap(
        "interference payload bytes",
        bytes,
        controls.maximum_payload_bytes,
    )
}

fn load_payload(
    n: usize,
    tones: usize,
    load: Option<HarmonicInterferenceLoad<'_>>,
) -> Result<usize, HarmonicResponseError> {
    let Some(load) = load else {
        return Ok(0);
    };
    // Retained original flow/effort per tone plus one loaded preparation/solve.
    let preparation = super::impedance::dense_payload(n, load.ports)?;
    tones
        .checked_mul(load.ports)
        .and_then(|x| x.checked_mul(32))
        .and_then(|x| x.checked_add(preparation))
        .ok_or(invalid("interference load payload extent overflow"))
}

#[allow(clippy::too_many_arguments)]
fn prepare_response(
    mass: &[f64],
    damping: &[f64],
    stiffness: &[f64],
    n: usize,
    omega: f64,
    load: Option<HarmonicInterferenceLoad<'_>>,
    controls: HarmonicResponseControls,
    cx: &Cx<'_>,
) -> Result<PreparedHarmonicResponse, HarmonicResponseError> {
    match load {
        Some(load) => PreparedHarmonicResponse::try_with_port_impedance(
            mass,
            damping,
            stiffness,
            n,
            load.force_map,
            load.ports,
            load.impedance(omega)?,
            omega,
            controls,
            load.controls,
            cx,
        ),
        None => PreparedHarmonicResponse::try_new(mass, damping, stiffness, n, omega, controls, cx),
    }
}

fn cap(
    quantity: &'static str,
    requested: usize,
    limit: usize,
) -> Result<(), HarmonicResponseError> {
    if requested > limit {
        Err(HarmonicResponseError::Budget {
            quantity,
            requested,
            cap: limit,
        })
    } else {
        Ok(())
    }
}

fn apply(
    matrix: &[f64],
    vector: &[C64],
    n: usize,
    viscous: bool,
    cx: &Cx<'_>,
) -> Result<Vec<[Dd; 2]>, HarmonicResponseError> {
    let mut result = Vec::with_capacity(n);
    for i in 0..n {
        poll(cx)?;
        let mut sum = [Dd::ZERO; 2];
        for j in 0..n {
            let coefficient = Dd::from_f64(if viscous {
                symmetric_entry(matrix, n, i, j)
            } else {
                matrix[i * n + j]
            });
            let value = parts(vector[j]);
            sum[0] = checked(sum[0] + coefficient * value[0])?;
            sum[1] = checked(sum[1] + coefficient * value[1])?;
        }
        result.push(sum);
    }
    Ok(result)
}

fn parts(z: C64) -> [Dd; 2] {
    [Dd::from_f64(z.re), Dd::from_f64(z.im)]
}
fn observation(error: WindowError) -> HarmonicResponseError {
    HarmonicResponseError::Observation(error)
}
fn product(kernel: HarmonicProduct, a: [Dd; 2], b: [Dd; 2]) -> Result<Dd, HarmonicResponseError> {
    kernel.mean_real_product_parts_dd(a, b).map_err(observation)
}
fn checked(value: Dd) -> Result<Dd, HarmonicResponseError> {
    if value.hi.is_finite() && value.lo.is_finite() {
        Ok(value)
    } else {
        Err(invalid("mechanical interference is not representable"))
    }
}
