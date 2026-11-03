//! Shared bounded RK2 stages for material vertices. Geometry is prepared once
//! per stage, independently of how physical sources are partitioned/stored.

use super::{Refusal, SelfInducedShed, WakeAdvectionReport, WakeTreeSpec, checkpoint, refusal};

#[path = "tree.rs"]
mod tree;
use crate::filament3d::SegmentGeometry;
use crate::images3d::ImagePlane;
use tree::SourceTree;

struct PreparedSource {
    geometry: SegmentGeometry,
    circulation: f64,
    core2_m2: f64,
}

pub(crate) struct PreparedSources {
    sources: Vec<PreparedSource>,
    bound: usize,
}

impl PreparedSources {
    pub(crate) fn push(
        &mut self,
        start: [f64; 3],
        end: [f64; 3],
        gamma: f64,
        core2: f64,
    ) -> Result<(), Refusal> {
        if !gamma.is_finite() {
            return Err(refusal(
                "wake-advection-circulation-invalid",
                "non-finite source strength",
            ));
        }
        if gamma == 0.0 {
            return Ok(());
        }
        let geometry = SegmentGeometry::new(start, end);
        if start.iter().chain(&end).any(|v| !v.is_finite())
            || !core2.is_finite()
            || core2 <= 0.0
            || core2 > f64::MAX / 16.0
            || !geometry.length.is_finite()
            || (geometry.length == 0.0 && start != end)
        {
            return Err(refusal(
                "wake-advection-core-invalid",
                "every active edge needs a finite positive core and representable length",
            ));
        }
        if self.sources.len() >= self.bound {
            return Err(refusal(
                "wake-advection-work-exhausted",
                "source assembly exceeds its admitted bound",
            ));
        }
        self.sources.push(PreparedSource {
            geometry,
            circulation: gamma,
            core2_m2: core2,
        });
        Ok(())
    }

    fn images(&mut self, plane: Option<ImagePlane>) -> Result<(), Refusal> {
        let Some(plane) = plane else { return Ok(()) };
        let count = self.sources.len();
        for i in 0..count {
            let s = &self.sources[i];
            let (start, end, gamma, core2) = (
                plane.mirror(s.geometry.start),
                plane.mirror(s.geometry.end),
                -s.circulation,
                s.core2_m2,
            );
            self.push(start, end, gamma, core2)?;
        }
        Ok(())
    }

    fn scale(&self, resolved: Option<&[bool]>) -> Result<f64, Refusal> {
        if resolved.is_some_and(|mask| mask.len() != self.sources.len()) {
            return Err(refusal(
                "wake-advection-topology-invalid",
                "stage source identity changed",
            ));
        }
        let mut scale = f64::INFINITY;
        for (i, s) in self.sources.iter().enumerate() {
            let active = resolved.map_or(s.geometry.length > 0.0, |mask| mask[i]);
            if active && s.geometry.length == 0.0 {
                return Err(refusal(
                    "wake-advection-core-invalid",
                    "a resolved source edge collapsed",
                ));
            }
            scale = scale.min(s.core2_m2.sqrt());
            if active {
                scale = scale.min(s.geometry.length);
            }
        }
        Ok(scale)
    }
}

pub(crate) struct StageSettings {
    pub dt_s: f64,
    /// Conservative physical source/traversal bound before optional images.
    pub source_bound: usize,
    pub image_plane: Option<ImagePlane>,
    pub line_step_m: f64,
}

pub(crate) fn integrate_nodes(
    points: &mut [[f64; 3]],
    config: SelfInducedShed,
    settings: StageSettings,
    tree_spec: Option<WakeTreeSpec>,
    mut assemble: impl FnMut(&[[f64; 3]], f64, &mut PreparedSources) -> Result<(), Refusal>,
    variation: &mut impl FnMut([f64; 3], f64) -> Result<[f64; 3], Refusal>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<WakeAdvectionReport, Refusal> {
    if tree_spec.is_some_and(|s| {
        !s.absolute_velocity_error_mps.is_finite()
            || s.absolute_velocity_error_mps < 0.0
            || !(1..=256).contains(&s.leaf_capacity)
    }) {
        return Err(refusal(
            "wake-tree-invalid",
            "finite nonnegative accuracy and leaf size 1..=256 required",
        ));
    }
    let bound = settings
        .source_bound
        .checked_mul(if settings.image_plane.is_some() { 2 } else { 1 });
    let admitted_work = bound.and_then(|bound| {
        let stages = u64::try_from(config.substeps).ok()?.checked_mul(2)?;
        let targets = u64::try_from(points.len()).ok()?;
        let sources = u64::try_from(bound).ok()?.checked_add(1)?;
        if tree_spec.is_some() {
            stages
                .checked_add(1)?
                .checked_mul(sources)?
                .checked_add(stages.checked_mul(targets)?)
        } else {
            stages
                .checked_mul(targets.checked_add(1)?)?
                .checked_add(1)?
                .checked_mul(sources)
        }
    });
    if admitted_work.is_none_or(|v| v > config.max_work) {
        return Err(refusal(
            "wake-advection-work-exhausted",
            "wake evolution exceeds the work cap",
        ));
    }
    let bound = bound.expect("checked work bound");
    let h = settings.dt_s / config.substeps as f64;
    if !settings.line_step_m.is_finite()
        || settings.line_step_m < 0.0
        || settings
            .image_plane
            .is_some_and(|plane| !plane.z0_m.is_finite())
    {
        return Err(refusal(
            "wake-advection-invalid",
            "invalid stage boundary or line displacement",
        ));
    }
    for &point in points.iter() {
        admit_point(point, settings.image_plane)?;
    }
    let mut sources = PreparedSources {
        sources: Vec::new(),
        bound,
    };
    sources
        .sources
        .try_reserve_exact(bound)
        .map_err(|_| allocation())?;
    let mut midpoint = Vec::new();
    midpoint
        .try_reserve_exact(points.len())
        .map_err(|_| allocation())?;
    midpoint.extend_from_slice(points);
    let mut resolved_edges = Vec::new();
    let mut tree = SourceTree::default();
    let mut report = WakeAdvectionReport::default();
    let mut work = WorkMeter {
        limit: config.max_work,
        used: 0,
        metered: tree_spec.is_some(),
    };
    for substep in 0..config.substeps {
        checkpoint(cancelled)?;
        let fraction = substep as f64 / config.substeps as f64;
        work.charge(bound as u64 + 1)?;
        sources.sources.clear();
        assemble(points, fraction, &mut sources)?;
        sources.images(settings.image_plane)?;
        let first_limit = sources.scale((substep > 0).then_some(resolved_edges.as_slice()))?
            * config.max_displacement_fraction;
        if substep == 0 {
            resolved_edges
                .try_reserve_exact(sources.sources.len())
                .map_err(|_| allocation())?;
            resolved_edges.extend(
                sources
                    .sources
                    .iter()
                    .map(|s| s.geometry.start != s.geometry.end),
            );
        }
        admit_line_step(settings.line_step_m, first_limit)?;
        if let Some(spec) = tree_spec.filter(|s| s.absolute_velocity_error_mps > 0.0) {
            tree.rebuild(&sources.sources, spec, &mut work, cancelled)?;
            report.peak_tree_nodes = report.peak_tree_nodes.max(tree.node_count());
        }
        for (p, q) in points.iter().zip(&mut midpoint) {
            let velocity = velocity_at(
                *p,
                fraction * settings.dt_s,
                &sources,
                config.ambient_velocity_mps,
                tree_spec,
                &mut tree,
                &mut work,
                &mut report,
                variation,
                cancelled,
            )?;
            let delta = displacement(velocity, h, first_limit)?;
            *q = translated(*p, delta.map(|v| 0.5 * v), settings.image_plane)?;
        }
        let fraction = (substep as f64 + 0.5) / config.substeps as f64;
        work.charge(bound as u64 + 1)?;
        sources.sources.clear();
        assemble(&midpoint, fraction, &mut sources)?;
        sources.images(settings.image_plane)?;
        let middle_limit = sources.scale(Some(&resolved_edges))? * config.max_displacement_fraction;
        admit_line_step(settings.line_step_m, middle_limit)?;
        if let Some(spec) = tree_spec.filter(|s| s.absolute_velocity_error_mps > 0.0) {
            tree.rebuild(&sources.sources, spec, &mut work, cancelled)?;
            report.peak_tree_nodes = report.peak_tree_nodes.max(tree.node_count());
        }
        for (p, q) in points.iter_mut().zip(&midpoint) {
            let velocity = velocity_at(
                *q,
                fraction * settings.dt_s,
                &sources,
                config.ambient_velocity_mps,
                tree_spec,
                &mut tree,
                &mut work,
                &mut report,
                variation,
                cancelled,
            )?;
            let delta = displacement(velocity, h, first_limit.min(middle_limit))?;
            *p = translated(*p, delta, settings.image_plane)?;
        }
    }
    // The midpoint is not the accepted endpoint. A strong resolved incident
    // strain can leave both sampled stages regular yet collapse an old edge
    // in the final RK update. Admit the final geometry and evolved cores before
    // publication, without another velocity sweep or incident-field sample.
    checkpoint(cancelled)?;
    work.charge(bound as u64 + 1)?;
    sources.sources.clear();
    assemble(points, 1.0, &mut sources)?;
    sources.images(settings.image_plane)?;
    sources.scale(Some(&resolved_edges))?;
    report.work_used = if tree_spec.is_some() {
        work.used
    } else {
        admitted_work.unwrap()
    };
    Ok(report)
}

#[allow(clippy::too_many_arguments)]
fn velocity_at(
    p: [f64; 3],
    elapsed_s: f64,
    sources: &PreparedSources,
    ambient: [f64; 3],
    tree_spec: Option<WakeTreeSpec>,
    tree: &mut SourceTree,
    work: &mut WorkMeter,
    report: &mut WakeAdvectionReport,
    variation: &mut impl FnMut([f64; 3], f64) -> Result<[f64; 3], Refusal>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<[f64; 3], Refusal> {
    checkpoint(cancelled)?;
    work.charge(1)?;
    let extra = variation(p, elapsed_s)?;
    let mut velocity = ambient;
    for c in 0..3 {
        // Avoid changing the default path's signed-zero arithmetic.
        if extra[c] != 0.0 {
            velocity[c] += extra[c];
        }
    }
    if velocity.iter().any(|v| !v.is_finite()) {
        return Err(refusal(
            "wake-advection-nonfinite",
            "incident velocity exceeds the numerical domain",
        ));
    }
    if let Some(spec) = tree_spec.filter(|s| s.absolute_velocity_error_mps > 0.0) {
        return tree.velocity(p, velocity, &sources.sources, spec, work, report, cancelled);
    }
    work.charge(sources.sources.len() as u64)?;
    report.direct_interactions += sources.sources.len() as u64;
    let mut correction = [0.0; 3];
    for (i, s) in sources.sources.iter().enumerate() {
        if i != 0 && i % 256 == 0 {
            checkpoint(cancelled)?;
        }
        let value = s.unit_velocity(p)?;
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

impl PreparedSource {
    #[inline]
    fn unit_velocity(&self, p: [f64; 3]) -> Result<[f64; 3], Refusal> {
        // Bound squared lengths and endpoint products before evaluating the
        // kernel: an overflowing denominator must not become apparent zero.
        let extent = self
            .geometry
            .start
            .iter()
            .zip(self.geometry.end)
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
        Ok(self.geometry.velocity(p, self.core2_m2))
    }
}

struct WorkMeter {
    limit: u64,
    used: u64,
    metered: bool,
}

impl WorkMeter {
    #[inline]
    fn charge(&mut self, amount: u64) -> Result<(), Refusal> {
        if self.metered {
            self.used = self
                .used
                .checked_add(amount)
                .filter(|n| *n <= self.limit)
                .ok_or_else(|| {
                    refusal(
                        "wake-advection-work-exhausted",
                        "metered tree evolution exhausted its work cap",
                    )
                })?;
        }
        Ok(())
    }
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
    admit_line_step(distance, limit)?;
    Ok(delta)
}

fn admit_line_step(distance: f64, limit: f64) -> Result<(), Refusal> {
    if distance > limit {
        Err(refusal(
            "wake-advection-displacement",
            "substep motion exceeds its core/edge resolution limit; increase substeps",
        ))
    } else {
        Ok(())
    }
}

fn translated(
    p: [f64; 3],
    delta: [f64; 3],
    plane: Option<ImagePlane>,
) -> Result<[f64; 3], Refusal> {
    let next = std::array::from_fn(|i| p[i] + delta[i]);
    admit_point(next, plane)?;
    Ok(next)
}

fn admit_point(p: [f64; 3], plane: Option<ImagePlane>) -> Result<(), Refusal> {
    if p.iter().any(|v| !v.is_finite()) {
        return Err(refusal(
            "wake-advection-nonfinite",
            "translated node exceeds the numerical domain",
        ));
    }
    if plane.is_some_and(|plane| p[2] <= plane.z0_m) {
        return Err(refusal(
            "image-plane-invalid",
            "a physical stage node reached or crossed the image plane",
        ));
    }
    Ok(())
}

fn allocation() -> Refusal {
    refusal(
        "wake-advection-allocation",
        "stage scratch allocation failed",
    )
}
