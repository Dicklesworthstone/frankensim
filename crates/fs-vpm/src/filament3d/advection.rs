//! Direct finite-core self-advection of the connected, attached near wake.
//! Circulation remains on material quads; moving their shared vertices
//! transports the vortex lines and their stretching without rescaling Γ.

use super::{FilamentWake, Refusal, SegmentGeometry};
use crate::farfield3d::WakeCoreEvolutionMode;

/// Prepared once per stage, then read by every target's direct sweep.
struct PreparedSource {
    geometry: SegmentGeometry,
    circulation: f64,
    core2_m2: f64,
}

/// Inputs to one self-induced shedding interval of length `1 / shed_hz`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SelfInducedShed {
    /// Uniform imposed fluid velocity [m/s], added to wake induction.
    pub ambient_velocity_mps: [f64; 3],
    /// Positive base radius [m] assigned to the row emitted at the endpoint.
    /// Existing rows retain their own base cores.
    pub emitted_core_radius_m: f64,
    /// Same age-dependent algebraic spreading law used by the far field.
    pub core_mode: WakeCoreEvolutionMode,
    /// Fixed positive number of midpoint RK2 substeps during this interval.
    pub substeps: usize,
    /// Maximum full-substep displacement divided by the smaller of the
    /// smallest active core radius and resolved active segment length, in (0, 1].
    /// Checked at both RK stages; this is a resolution guard, not an error bound.
    /// Edges born with zero length remain subject to the core guard and join
    /// the edge-length guard in the next shedding interval.
    pub max_displacement_fraction: f64,
    /// Upper bound on direct interaction work, including source assembly and
    /// target updates: `2 * substeps * (targets + 1) * (source_bound + 1)`.
    /// The source bound includes zero-strength edges, so admission is conservative.
    pub max_work: u64,
}

impl FilamentWake {
    /// Convect the existing near-wake nodes by their own finite-core velocity
    /// plus uniform ambient flow, keep the lifting line fixed, and shed one row.
    ///
    /// Every substep rebuilds the physical source sheet at the midpoint before
    /// evaluating its second-stage velocity. Cores use the initial and midpoint
    /// elapsed ages. Stored circulations and base cores are retained, existing
    /// ages advance exactly one shedding interval, and the new row has age zero.
    /// The old sheet strength is held during the interval; `gamma` is captured
    /// at its endpoint. RK2 order applies to that smooth held-strength problem,
    /// not to discontinuous boundary-circulation changes.
    ///
    /// The anchored line exchanges impulse with the wake. This method does not
    /// assert exact impulse/energy conservation or resolved viscous tube strain.
    /// It is a direct O(nodes × segments × substeps) solver, without far cells,
    /// image boundaries, spatially varying imposed flow or automatic remeshing.
    ///
    /// `cancelled` is polled before work, at every target, within at most 256
    /// source-target interactions, and immediately before publication. Every
    /// error or cancellation leaves the original wake unchanged.
    ///
    /// # Errors
    /// Invalid configuration, geometry, finite positive source cores, work cap,
    /// excessive stage displacement, non-finite numerics, or cancellation.
    /// Attached far-field interfaces refuse: advancing one side independently
    /// would disconnect its physical filaments.
    pub fn shed_self_induced(
        &mut self,
        gamma: &[f64],
        config: SelfInducedShed,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), Refusal> {
        checkpoint(&mut cancelled)?;
        self.cert.admit()?;
        self.admit_geometry()?;
        config.core_mode.admit()?;
        let dt = 1.0 / self.cert.shed_hz;
        let h = dt / config.substeps as f64;
        let emitted_core2 = config.emitted_core_radius_m.powi(2);
        let growth = config.core_mode.core_growth_m2(dt);
        if config.substeps == 0
            || !h.is_finite()
            || h <= 0.0
            || !dt.is_finite()
            || !growth.is_finite()
            || config.ambient_velocity_mps.iter().any(|v| !v.is_finite())
            || !config.emitted_core_radius_m.is_finite()
            || config.emitted_core_radius_m <= 0.0
            || !emitted_core2.is_finite()
            || emitted_core2 <= 0.0
            || emitted_core2 > f64::MAX / 16.0
            || !config.max_displacement_fraction.is_finite()
            || !(0.0 < config.max_displacement_fraction && config.max_displacement_fraction <= 1.0)
        {
            return Err(refusal(
                "wake-advection-invalid",
                "invalid time, core or step settings",
            ));
        }
        if self.downstream_gamma.is_some() {
            return Err(refusal(
                "wake-advection-interface-attached",
                "self-advection needs the complete connected near wake",
            ));
        }
        // Include linear work so a zero-circulation sheet cannot bypass the
        // admission by requesting arbitrarily many nominally zero-pair steps.
        let targets = self.rows.len() * (self.cert.n_stations + 1);
        let source_bound = self.rows.len() * (2 * self.cert.n_stations + 1) + self.cert.n_stations;
        let work = (config.substeps as u64)
            .checked_mul(2)
            .and_then(|v| v.checked_mul(targets as u64 + 1))
            .and_then(|v| v.checked_mul(source_bound as u64 + 1));
        if work.is_none_or(|v| v > config.max_work) {
            return Err(refusal(
                "wake-advection-work-exhausted",
                "direct sweep exceeds the work cap",
            ));
        }
        for row in &self.rows {
            let final_core = row.core2_m2 + growth * (row.age_ticks + 1.0);
            if !row.core2_m2.is_finite()
                || row.core2_m2 < 0.0
                || !row.age_ticks.is_finite()
                || row.age_ticks < 0.0
                || !final_core.is_finite()
            {
                return Err(refusal(
                    "wake-advection-core-invalid",
                    "invalid retained age or core",
                ));
            }
        }
        // Reuse the existing shed admission and endpoint construction on a
        // private candidate. No published source moves until all stages pass.
        let mut state = self.clone();
        state.shed(gamma, [0.0; 3])?;
        let mut emitted = state.rows.pop().expect("shed appended one row");
        emitted.core2_m2 = emitted_core2;
        if emitted.gamma.windows(2).any(|g| !(g[1] - g[0]).is_finite())
            || state.rows.last().is_some_and(|row| {
                row.gamma
                    .iter()
                    .zip(&emitted.gamma)
                    .any(|(a, b)| !(a - b).is_finite())
            })
        {
            return Err(refusal(
                "wake-advection-circulation-invalid",
                "emission would overflow a shared-edge circulation",
            ));
        }
        for (row, original) in state.rows.iter_mut().zip(&self.rows) {
            row.age_ticks = original.age_ticks;
        }
        let mut midpoint = state.clone();
        let mut sources = Vec::new();
        let mut resolved_edges = Vec::new();
        sources.try_reserve_exact(source_bound).map_err(|_| {
            refusal(
                "wake-advection-allocation",
                "source scratch allocation failed",
            )
        })?;
        for substep in 0..config.substeps {
            checkpoint(&mut cancelled)?;
            let fraction = substep as f64 / config.substeps as f64;
            let first_limit = collect_sources(
                &state,
                growth,
                fraction,
                &mut sources,
                (substep > 0).then_some(resolved_edges.as_slice()),
            )? * config.max_displacement_fraction;
            if substep == 0 {
                resolved_edges = sources
                    .iter()
                    .map(|s| s.geometry.start != s.geometry.end)
                    .collect();
            }
            for (row, midrow) in state.rows.iter().zip(&mut midpoint.rows) {
                for (p, q) in row.nodes.iter().zip(&mut midrow.nodes) {
                    let velocity =
                        velocity_at(*p, &sources, config.ambient_velocity_mps, &mut cancelled)?;
                    let delta = displacement(velocity, h, first_limit)?;
                    *q = translated(*p, delta.map(|v| 0.5 * v))?;
                }
            }
            let fraction = (substep as f64 + 0.5) / config.substeps as f64;
            let middle_limit = collect_sources(
                &midpoint,
                growth,
                fraction,
                &mut sources,
                Some(&resolved_edges),
            )? * config.max_displacement_fraction;
            for (row, midrow) in state.rows.iter_mut().zip(&midpoint.rows) {
                for (p, q) in row.nodes.iter_mut().zip(&midrow.nodes) {
                    let velocity =
                        velocity_at(*q, &sources, config.ambient_velocity_mps, &mut cancelled)?;
                    let delta = displacement(velocity, h, first_limit.min(middle_limit))?;
                    *p = translated(*p, delta)?;
                }
            }
        }
        for (row, original) in state.rows.iter_mut().zip(&self.rows) {
            row.age_ticks = original.age_ticks + 1.0;
        }
        state.rows.push(emitted);
        checkpoint(&mut cancelled)?;
        *self = state;
        Ok(())
    }
}

fn collect_sources(
    wake: &FilamentWake,
    growth: f64,
    extra_age: f64,
    sources: &mut Vec<PreparedSource>,
    resolved_edges: Option<&[bool]>,
) -> Result<f64, Refusal> {
    sources.clear();
    let mut scale = f64::INFINITY;
    let mut valid = true;
    wake.visit_segments_with_core_growth(growth, |mut source| {
        source.core2_m2 += growth * extra_age;
        let geometry = SegmentGeometry::new(source.start, source.end);
        let resolved = resolved_edges.map_or(geometry.length > 0.0, |mask| mask[sources.len()]);
        valid &= source.circulation.is_finite()
            && source.core2_m2.is_finite()
            && source.core2_m2 > 0.0
            && source.core2_m2 <= f64::MAX / 16.0
            && geometry.length.is_finite()
            && (geometry.length > 0.0 || (!resolved && source.start == source.end));
        scale = scale.min(source.core2_m2.sqrt());
        if resolved {
            scale = scale.min(geometry.length);
        }
        sources.push(PreparedSource {
            geometry,
            circulation: source.circulation,
            core2_m2: source.core2_m2,
        });
    });
    if !valid {
        return Err(refusal(
            "wake-advection-core-invalid",
            "every active edge needs a finite positive core and representable length",
        ));
    }
    Ok(scale)
}

fn velocity_at(
    p: [f64; 3],
    sources: &[PreparedSource],
    ambient: [f64; 3],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<[f64; 3], Refusal> {
    checkpoint(cancelled)?;
    let mut velocity = ambient;
    let mut correction = [0.0; 3];
    for (i, s) in sources.iter().enumerate() {
        if i != 0 && i % 256 == 0 {
            checkpoint(cancelled)?;
        }
        // Keep the finite-segment kernel's squared lengths and endpoint
        // products in range; checking only its final velocity could conceal
        // an overflow to an apparent zero far-field contribution.
        let extent = s
            .geometry
            .start
            .iter()
            .zip(s.geometry.end)
            .zip(p)
            .fold(0.0_f64, |v, ((a, b), x)| {
                v.max((x - a).abs()).max((x - b).abs())
            });
        if !extent.is_finite() || extent > f64::MAX.sqrt() / 16.0 {
            return Err(refusal(
                "wake-advection-nonfinite",
                "source-target separation exceeds the numerical domain",
            ));
        }
        let value = s.geometry.velocity(p, s.core2_m2);
        for c in 0..3 {
            let increment = s.circulation * value[c] - correction[c];
            let sum = velocity[c] + increment;
            correction[c] = (sum - velocity[c]) - increment;
            velocity[c] = sum;
        }
    }
    if velocity.iter().any(|v| !v.is_finite()) {
        return Err(refusal(
            "wake-advection-nonfinite",
            "induced velocity exceeds the numerical domain",
        ));
    }
    Ok(velocity)
}

fn displacement(velocity: [f64; 3], h: f64, limit: f64) -> Result<[f64; 3], Refusal> {
    let delta = velocity.map(|v| h * v);
    let distance = delta[0].hypot(delta[1]).hypot(delta[2]);
    if !distance.is_finite() {
        return Err(refusal(
            "wake-advection-nonfinite",
            "node displacement exceeds the numerical domain",
        ));
    }
    if distance > limit {
        return Err(refusal(
            "wake-advection-displacement",
            "substep motion exceeds its core/edge resolution limit; increase substeps",
        ));
    }
    Ok(delta)
}

fn translated(p: [f64; 3], delta: [f64; 3]) -> Result<[f64; 3], Refusal> {
    let next = std::array::from_fn(|i| p[i] + delta[i]);
    if next.iter().any(|v| !v.is_finite()) {
        return Err(refusal(
            "wake-advection-nonfinite",
            "translated node exceeds the numerical domain",
        ));
    }
    Ok(next)
}

fn checkpoint(cancelled: &mut impl FnMut() -> bool) -> Result<(), Refusal> {
    if cancelled() {
        Err(refusal(
            "wake-advection-cancelled",
            "self-induced shedding was cancelled",
        ))
    } else {
        Ok(())
    }
}

fn refusal(code: &'static str, message: &str) -> Refusal {
    Refusal {
        code,
        message: message.into(),
        ranked_repairs: vec![
            "use a complete finite-core wake and admitted integration settings".into(),
        ],
    }
}
