//! Replayable conditional OU paths between two original fixed-clock states.
//!
//! A canonical binary tree supplies shared Gaussian bridge innovations, so
//! querying another stage or retrying a smaller step never changes a value
//! already observed. The tree is traversed only along the requested time's
//! ancestry; no complete tree, per-query cache, or interpolated amplitudes
//! are needed. Gaussian conditioning preserves the joint OU law, including
//! its nonzero conditional variance between accepted endpoints.
//!
//! Physical endpoint subtraction can differ from the owner's nominal dt
//! through rounding. The OU clock maps the represented physical interval
//! affinely onto that nominal dt; spatial phase always uses the requested
//! physical time. This preserves the original fixed-step endpoint law.

use super::*;

/// Separate counter namespace for conditional tree innovations. Neither
/// fixed-clock endpoint nor interval-average streams are changed.
pub const MANN_OU_BRIDGE_KERNEL: u32 = 0x4d42_5247;
/// Hard bound on the stack-resident query ancestry and per-mode work. A
/// requested time that needs a deeper tree is refused, never snapped.
pub const MAX_MANN_OU_BRIDGE_DEPTH: usize = 128;

/// Immutable interval retaining the accepted state and its exact original
/// next endpoint. Preparation does not advance or modify the supplied owner.
#[derive(Debug)]
pub struct MannOuBridge {
    definition: Arc<Definition>,
    step: u64,
    bridge_seed: u64,
    maximum_depth: usize,
    start_time_s: f64,
    end_time_s: f64,
    left: [Vec<f64>; 4],
    right: [Vec<f64>; 4],
}

/// Reusable observation storage. A frame binds the original spectrum, fixed
/// clock and realization seeds, and the conditional bridge seed; it can be
/// reused for successive intervals of that source. Accepted amplitudes and
/// time change together only after a successful observation. Sampling the
/// frame and updating it allocate no numerical storage.
#[derive(Debug)]
pub struct MannOuBridgeFrame {
    definition: Arc<Definition>,
    bridge_seed: u64,
    time_s: f64,
    amplitudes: [Vec<f64>; 4],
    trial: [Vec<f64>; 4],
}

/// Actual work for one successful frame observation. A normal draw consumes
/// two uniform counter positions, as specified by `fs-rand`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MannOuBridgeSampleReport {
    pub sampled_time_s: f64,
    pub normal_draws: usize,
    pub tree_depth: usize,
}

#[derive(Clone, Copy, Debug, Default)]
struct Node {
    physical_time_s: f64,
    left_duration_s: f64,
    right_duration_s: f64,
    follow_right: bool,
}

struct QueryPath {
    nodes: [Node; MAX_MANN_OU_BRIDGE_DEPTH],
    depth: usize,
}

fn bridge_poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Refusal> {
    if cancelled() {
        Err(error(
            "mann-ou-bridge-cancelled",
            "OU bridge operation cancelled",
        ))
    } else {
        Ok(())
    }
}

fn numerical_bytes(
    definition: &Definition,
    first: &[Vec<f64>; 4],
    second: &[Vec<f64>; 4],
) -> usize {
    (definition.times.capacity()
        + 2 * definition.transition.capacity()
        + first.iter().chain(second).map(Vec::capacity).sum::<usize>())
        * core::mem::size_of::<f64>()
}

impl MannOuBridge {
    /// Prepare one exact original endpoint (`4 * mode_count` normal draws),
    /// retaining the accepted endpoint unchanged. The explicit depth bound
    /// admits future query ancestry, not a preallocated exponential tree.
    pub fn prepare(
        source: &FixedStepMannField,
        bridge_seed: u64,
        maximum_depth: usize,
        maximum_workspace_bytes: usize,
        maximum_normal_draws: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Refusal> {
        bridge_poll(&mut cancelled)?;
        if !(1..=MAX_MANN_OU_BRIDGE_DEPTH).contains(&maximum_depth) {
            return Err(error(
                "mann-ou-bridge-depth",
                "invalid OU bridge depth bound",
            ));
        }
        let next_step = source.step.checked_add(1).ok_or_else(|| {
            error(
                "mann-ou-bridge-clock",
                "next fixed-clock endpoint overflows",
            )
        })?;
        let start_time_s = time_at(source.origin_time_s(), source.dt_s(), source.step)?;
        let end_time_s = time_at(source.origin_time_s(), source.dt_s(), next_step)?;
        let physical_duration = end_time_s - start_time_s;
        if !physical_duration.is_finite() || physical_duration <= 0.0 {
            return Err(error(
                "mann-ou-bridge-clock",
                "physical endpoint gap is not representable",
            ));
        }
        let n = source.model().mode_count();
        if Self::required_workspace_bytes(n)? > maximum_workspace_bytes {
            return Err(error(
                "mann-ou-bridge-budget",
                "bridge endpoint storage exceeds its byte budget",
            ));
        }
        if 4 * n > maximum_normal_draws {
            return Err(error(
                "mann-ou-bridge-work",
                "next endpoint exceeds its normal-draw budget",
            ));
        }
        let mut result = Self {
            definition: Arc::clone(&source.definition),
            step: source.step,
            bridge_seed,
            maximum_depth,
            start_time_s,
            end_time_s,
            left: [zeros(n)?, zeros(n)?, zeros(n)?, zeros(n)?],
            right: [zeros(n)?, zeros(n)?, zeros(n)?, zeros(n)?],
        };
        if result.workspace_bytes() > maximum_workspace_bytes {
            return Err(error(
                "mann-ou-bridge-budget",
                "allocated bridge capacities exceed the byte budget",
            ));
        }
        for path in 0..4 {
            result.left[path].copy_from_slice(&source.amplitudes[path]);
            for (mode, &[rho, sigma]) in source.definition.transition.iter().enumerate() {
                if mode % 64 == 0 {
                    bridge_poll(&mut cancelled)?;
                }
                // This address and arithmetic are exactly advance_to(s+1).
                let key = StreamKey {
                    seed: source.definition.evolution_seed ^ PATH_KEYS[path],
                    kernel: FIXED_STEP_MANN_KERNEL,
                    tile: mode as u32,
                };
                let mut rng = Stream::resume(StreamCheckpoint::current(key, 2 * next_step))
                    .expect("current-version checkpoint and admitted counter always resume");
                let value = rho * source.amplitudes[path][mode] + sigma * rng.next_normal();
                if !value.is_finite() {
                    return Err(error(
                        "mann-ou-bridge-state",
                        "next endpoint amplitude is not finite",
                    ));
                }
                result.right[path][mode] = value;
            }
        }
        bridge_poll(&mut cancelled)?;
        Ok(result)
    }

    /// Nominal numerical storage is 88 bytes/mode: eight endpoint arrays and
    /// the three original shared temporal arrays. Immutable spectral geometry,
    /// fixed metadata, allocator metadata and bounded stack ancestry are excluded.
    pub fn required_workspace_bytes(modes: usize) -> Result<usize, Refusal> {
        FixedStepMannField::required_workspace_bytes(modes)
    }

    /// Actual capacities, conservatively counting shared temporal arrays once
    /// for each owner. The retained source itself is budgeted independently.
    pub fn workspace_bytes(&self) -> usize {
        numerical_bytes(&self.definition, &self.left, &self.right)
    }

    pub fn model(&self) -> &MannSpectralModel {
        &self.definition.model
    }
    pub fn correlation_times_s(&self) -> &[f64] {
        &self.definition.times
    }
    pub fn advection_mps(&self) -> [f64; 3] {
        self.model().advection_mps
    }
    pub const fn start_time_s(&self) -> f64 {
        self.start_time_s
    }
    pub const fn end_time_s(&self) -> f64 {
        self.end_time_s
    }
    /// Nominal OU duration. The physical endpoint difference can differ by rounding.
    pub fn dt_s(&self) -> f64 {
        self.definition.dt_s
    }
    pub const fn start_step(&self) -> u64 {
        self.step
    }
    pub const fn maximum_depth(&self) -> usize {
        self.maximum_depth
    }

    /// Allocate admitted active/trial storage initially equal to the left
    /// endpoint. No random draws or change of the original owner is required.
    pub fn frame(
        &self,
        maximum_workspace_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MannOuBridgeFrame, Refusal> {
        bridge_poll(&mut cancelled)?;
        let n = self.model().mode_count();
        if Self::required_workspace_bytes(n)? > maximum_workspace_bytes {
            return Err(error(
                "mann-ou-bridge-budget",
                "observation storage exceeds its byte budget",
            ));
        }
        let mut result = MannOuBridgeFrame {
            definition: Arc::clone(&self.definition),
            bridge_seed: self.bridge_seed,
            time_s: self.start_time_s,
            amplitudes: [zeros(n)?, zeros(n)?, zeros(n)?, zeros(n)?],
            trial: [zeros(n)?, zeros(n)?, zeros(n)?, zeros(n)?],
        };
        if result.workspace_bytes() > maximum_workspace_bytes {
            return Err(error(
                "mann-ou-bridge-budget",
                "allocated observation capacities exceed the byte budget",
            ));
        }
        for (target, left) in result.amplitudes.iter_mut().zip(&self.left) {
            bridge_poll(&mut cancelled)?;
            target.copy_from_slice(left);
        }
        bridge_poll(&mut cancelled)?;
        Ok(result)
    }

    /// Exact normal-draw cost of one observation, independent of past query
    /// order. Endpoints cost zero. This also checks time, representable clock
    /// durations and ancestry depth, before any draws or frame mutation.
    pub fn required_normal_draws(&self, time_s: f64) -> Result<usize, Refusal> {
        Ok(self.query_path(time_s)?.depth * 4 * self.model().mode_count())
    }

    fn query_path(&self, time_s: f64) -> Result<QueryPath, Refusal> {
        if !time_s.is_finite() || time_s < self.start_time_s || time_s > self.end_time_s {
            return Err(error(
                "mann-ou-bridge-time",
                "observation is outside the retained interval",
            ));
        }
        let mut result = QueryPath {
            nodes: [Node::default(); MAX_MANN_OU_BRIDGE_DEPTH],
            depth: 0,
        };
        if time_s == self.start_time_s || time_s == self.end_time_s {
            return Ok(result);
        }
        let physical_duration = self.end_time_s - self.start_time_s;
        let mut lo = self.start_time_s;
        let mut hi = self.end_time_s;
        loop {
            if result.depth == self.maximum_depth {
                return Err(error(
                    "mann-ou-bridge-depth",
                    "observation ancestry exceeds the declared depth",
                ));
            }
            let mid = lo + 0.5 * (hi - lo);
            if mid <= lo || mid >= hi {
                return Err(error(
                    "mann-ou-bridge-clock",
                    "interior tree midpoint is not representable",
                ));
            }
            // Local differences avoid loss from subtracting two nearly equal
            // absolute OU offsets near the interval's far endpoint.
            let left_duration_s = ((mid - lo) / physical_duration) * self.dt_s();
            let right_duration_s = ((hi - mid) / physical_duration) * self.dt_s();
            if !left_duration_s.is_finite()
                || left_duration_s <= 0.0
                || !right_duration_s.is_finite()
                || right_duration_s <= 0.0
            {
                return Err(error(
                    "mann-ou-bridge-clock",
                    "interior OU duration is not representable",
                ));
            }
            result.nodes[result.depth] = Node {
                physical_time_s: mid,
                left_duration_s,
                right_duration_s,
                follow_right: time_s > mid,
            };
            result.depth += 1;
            if time_s == mid {
                return Ok(result);
            }
            if time_s < mid {
                hi = mid;
            } else {
                lo = mid;
            }
        }
    }

    /// Atomically observe the same retained Gaussian path at a requested
    /// physical time. A replay consumes the same bounded work and returns the
    /// same bits; neither order, rejected CFD stages, nor retries alter it.
    /// Cancellation is checked at most every 64 normal draws and immediately
    /// before accepted frame publication. No successful call allocates.
    pub fn sample_frame(
        &self,
        time_s: f64,
        frame: &mut MannOuBridgeFrame,
        maximum_normal_draws: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MannOuBridgeSampleReport, Refusal> {
        bridge_poll(&mut cancelled)?;
        if frame.bridge_seed != self.bridge_seed
            || (!Arc::ptr_eq(&frame.definition, &self.definition)
                && frame.definition != self.definition)
        {
            return Err(error(
                "mann-ou-bridge-source-mismatch",
                "observation frame belongs to another source or bridge seed",
            ));
        }
        let path = self.query_path(time_s)?;
        let normal_draws = path.depth * 4 * self.model().mode_count();
        if normal_draws > maximum_normal_draws {
            return Err(error(
                "mann-ou-bridge-work",
                "observation exceeds its normal-draw budget",
            ));
        }
        if path.depth == 0 {
            let endpoint = if time_s == self.start_time_s {
                &self.left
            } else {
                &self.right
            };
            for (target, source) in frame.trial.iter_mut().zip(endpoint) {
                bridge_poll(&mut cancelled)?;
                target.copy_from_slice(source);
            }
        } else {
            let mut completed_draws = 0;
            let mut coefficients_by_node = [[0.0; 3]; MAX_MANN_OU_BRIDGE_DEPTH];
            let mut previous_tau = None;
            for (mode, &tau) in self.definition.times.iter().enumerate() {
                // Equal adjacent lifetimes are common. Reuse their bounded
                // coefficient table rather than repeat transcendental work.
                if previous_tau != Some(tau) {
                    for (depth, node) in path.nodes[..path.depth].iter().enumerate() {
                        if depth % 16 == 0 {
                            bridge_poll(&mut cancelled)?;
                        }
                        coefficients_by_node[depth] =
                            coefficients(node.left_duration_s, node.right_duration_s, tau)?;
                    }
                    previous_tau = Some(tau);
                }
                let mut left = core::array::from_fn::<_, 4, _>(|p| self.left[p][mode]);
                let mut right = core::array::from_fn::<_, 4, _>(|p| self.right[p][mode]);
                for (depth, node) in path.nodes[..path.depth].iter().enumerate() {
                    if completed_draws % 64 == 0 {
                        bridge_poll(&mut cancelled)?;
                    }
                    let coefficients = coefficients_by_node[depth];
                    for p in 0..4 {
                        let value = coefficients[0] * left[p]
                            + coefficients[1] * right[p]
                            + coefficients[2] * self.normal(node.physical_time_s, mode, p);
                        if !value.is_finite() {
                            return Err(error(
                                "mann-ou-bridge-state",
                                "conditional amplitude is not finite",
                            ));
                        }
                        frame.trial[p][mode] = value;
                        if node.follow_right {
                            left[p] = value;
                        } else {
                            right[p] = value;
                        }
                    }
                    completed_draws += 4;
                }
            }
        }
        bridge_poll(&mut cancelled)?;
        core::mem::swap(&mut frame.amplitudes, &mut frame.trial);
        // Preserve the requested physical phase, without clock snapping.
        frame.time_s = time_s;
        Ok(MannOuBridgeSampleReport {
            sampled_time_s: time_s,
            normal_draws,
            tree_depth: path.depth,
        })
    }

    fn normal(&self, time_s: f64, mode: usize, path: usize) -> f64 {
        // For a fixed bridge seed, distinct interval steps give distinct
        // keys. Within one interval, finite midpoint magnitudes address
        // disjoint even counter pairs; sign, mode and quadrature occupy the
        // tile. Finite magnitude bits are < 2^63, so both draws fit in u64.
        // Canonicalize signed zero as the same physical node.
        let bits = if time_s == 0.0 { 0 } else { time_s.to_bits() };
        let key = StreamKey {
            seed: self.bridge_seed ^ self.step,
            kernel: MANN_OU_BRIDGE_KERNEL,
            tile: ((mode * 4 + path) * 2 + (bits >> 63) as usize) as u32,
        };
        let index = 2 * (bits & 0x7fff_ffff_ffff_ffff);
        Stream::resume(StreamCheckpoint::current(key, index))
            .expect("current-version checkpoint and finite midpoint counter always resume")
            .next_normal()
    }
}

/// Gaussian Markov conditioning with stable innovations, including tiny
/// positive durations for which exp(-duration/tau) itself rounds to one.
fn coefficients(left_s: f64, right_s: f64, tau_s: f64) -> Result<[f64; 3], Refusal> {
    let left = left_s / tau_s;
    let right = right_s / tau_s;
    let a = det::exp(-left);
    let b = det::exp(-right);
    let q_left = -det::expm1(-2.0 * left);
    let q_right = -det::expm1(-2.0 * right);
    let q_total = -det::expm1(-2.0 * (left + right));
    // Divide before multiplication to avoid premature variance underflow.
    let variance = (q_left / q_total) * q_right;
    let result = [
        a * (q_right / q_total),
        b * (q_left / q_total),
        det::sqrt(variance),
    ];
    if left <= 0.0
        || right <= 0.0
        || q_left <= 0.0
        || q_right <= 0.0
        || result.iter().any(|x| !x.is_finite())
        || result[2] <= 0.0
    {
        return Err(error(
            "mann-ou-bridge-coefficients",
            "conditional OU transition is not representable",
        ));
    }
    Ok(result)
}

impl MannOuBridgeFrame {
    /// Same nominal 88 bytes/mode policy as bridge preparation.
    pub fn required_workspace_bytes(modes: usize) -> Result<usize, Refusal> {
        MannOuBridge::required_workspace_bytes(modes)
    }
    pub fn workspace_bytes(&self) -> usize {
        numerical_bytes(&self.definition, &self.amplitudes, &self.trial)
    }
    pub fn model(&self) -> &MannSpectralModel {
        &self.definition.model
    }
    pub fn correlation_times_s(&self) -> &[f64] {
        &self.definition.times
    }
    pub fn advection_mps(&self) -> [f64; 3] {
        self.model().advection_mps
    }
    pub const fn time_s(&self) -> f64 {
        self.time_s
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
    /// Retain the original harmonic slip-plane observation at the sampled
    /// physical time, including its analytic spatial derivatives.
    pub fn sample_wall_blocked(
        &self,
        point_m: [f64; 3],
        wall: MannPlaneWall,
    ) -> Result<FlowSample, Refusal> {
        wall::sample_modes(
            self.model(),
            self.time_s,
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
            self.time_s,
            core::array::from_fn(|i| self.amplitudes[i].as_slice()),
            point_m,
            filter,
        )
    }
    fn sample_with_gains(
        &self,
        point_m: [f64; 3],
        gains: Option<&[f64]>,
    ) -> Result<FlowSample, Refusal> {
        sample_modes(
            self.model(),
            self.time_s,
            core::array::from_fn(|i| self.amplitudes[i].as_slice()),
            point_m,
            gains,
        )
    }
}

#[cfg(test)]
mod tests;
