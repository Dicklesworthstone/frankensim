//! Coupled direct evolution of retained near and far physical filaments.
//! Multipoles are evaluation caches; their exact source endpoints, rather
//! than the expansion centers, follow the common material-vertex integrator.

use std::collections::BTreeMap;

use super::{FarCell, FarField, cell_from_segs_with_cores};
use crate::filament3d::advection::stage::{PreparedSources, StageSettings, integrate_nodes};
use crate::filament3d::advection::{LineMotion, append_near_sources, prepare_near};
use crate::filament3d::{FilamentWake, Refusal, SelfInducedShed, WakeBoundary};

impl FarField {
    /// Advance all retained near/far filaments by their common finite-core
    /// induced velocity and uniform ambient flow, then emit one near row.
    ///
    /// Shared free endpoints are evolved once, keeping the near/far interface
    /// and inter-cell vortex lines connected through both midpoint RK2 stages.
    /// All circulations and base source cores survive unchanged; source ages
    /// advance together and every far expansion is rebuilt from moved sources.
    ///
    /// This direct reference path costs O(nodes × sources × substeps). The
    /// exact retained far segments provide the sources, including inside far
    /// cells where a multipole approximation is not admissible. It does not
    /// prune, remesh, approximate the far dynamics, or claim exact impulse
    /// conservation for a wake exchanging impulse with its prescribed line.
    ///
    /// # Errors
    /// Incompatible near/far time or core modes, incomplete vortex lines,
    /// invalid geometry/cores, unrepresentable moments, work/displacement
    /// limits, or cancellation. Both inputs remain unchanged on every error.
    pub fn shed_self_induced(
        &mut self,
        near: &mut FilamentWake,
        gamma: &[f64],
        config: SelfInducedShed,
        cancelled: impl FnMut() -> bool,
    ) -> Result<(), Refusal> {
        self.shed_self_induced_with_environment(
            near,
            gamma,
            config,
            WakeBoundary::default(),
            |_, _| Ok([0.0; 3]),
            cancelled,
        )
    }

    /// Coupled near/far evolution with the same prescribed line motion,
    /// stage-dependent incident velocity and image-plane treatment as the
    /// corresponding near-only method. `variation` adds to the uniform
    /// velocity and receives the physical point and elapsed seconds.
    ///
    /// The far-field shedding interval must equal `1 / near.cert.shed_hz`,
    /// and its core mode must equal `config.core_mode`. Previously pruned or
    /// independently moved inputs that leave unbalanced vortex endpoints
    /// refuse. Rebuild any separately stored image snapshots after success.
    ///
    /// # Errors
    /// The same refusals as [`Self::shed_self_induced`], plus an invalid
    /// boundary, incident-field refusal, or a physical image-plane crossing.
    pub fn shed_self_induced_with_environment(
        &mut self,
        near: &mut FilamentWake,
        gamma: &[f64],
        config: SelfInducedShed,
        boundary: WakeBoundary<'_>,
        mut variation: impl FnMut([f64; 3], f64) -> Result<[f64; 3], Refusal>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), Refusal> {
        checkpoint(&mut cancelled)?;
        let proposal = prepare_near(near, gamma, config)?;
        if self.shed_dt_s != proposal.dt_s || self.mode != config.core_mode {
            return Err(refusal(
                "wake-advection-schedule-mismatch",
                "near and far must share the shedding interval and core mode",
            ));
        }
        if near.downstream_gamma.as_ref().is_some_and(|g| {
            near.rows.is_empty()
                || g.len() != near.cert.n_stations
                || g.iter().any(|v| !v.is_finite())
        }) {
            return Err(refusal(
                "wake-advection-topology-invalid",
                "the attached near boundary has invalid circulation or no retained row",
            ));
        }
        let growth = proposal.growth_per_tick_m2;
        let source_bound = admit_sources(self, near, growth, config.max_work, &mut cancelled)?;
        admit_circulation(near, self, &mut cancelled)?;
        let (mut points, bindings) = shared_nodes(near, self, &mut cancelled)?;
        let mut candidate = proposal.near;
        let line = LineMotion::new(&candidate, boundary)?;
        let settings = StageSettings {
            dt_s: proposal.dt_s,
            source_bound,
            image_plane: boundary.image_plane,
            line_step_m: line.step_distance(config.substeps),
        };
        integrate_nodes(
            &mut points,
            config,
            settings,
            |stage_points, fraction, sources| {
                bindings.scatter_near(&mut candidate, stage_points);
                line.apply(&mut candidate, fraction)?;
                // Preserve the physical source order of the pre-aggregation
                // sheet. This also preserves compensated summation order.
                append_far_sources(self, &bindings, stage_points, growth, fraction, sources)?;
                append_near_sources(&candidate, growth, fraction, sources)
            },
            &mut variation,
            &mut cancelled,
        )?;
        bindings.scatter_near(&mut candidate, &points);
        line.apply(&mut candidate, 1.0)?;
        for (row, original) in candidate.rows.iter_mut().zip(&near.rows) {
            row.age_ticks = original.age_ticks + 1.0;
        }
        let mut emitted = proposal.emitted;
        emitted.nodes.clone_from(&candidate.line_nodes);
        candidate.rows.push(emitted);
        let cells = rebuild_cells(self, &bindings, &points, growth, &mut cancelled)?;
        checkpoint(&mut cancelled)?;
        // No fallible operation remains between these two publications.
        *near = candidate;
        self.cells = cells;
        Ok(())
    }
}

/// Bind coordinate-identical FREE vertices, never prescribed line vertices.
/// A newborn free row can coincide with the line without being pinned to it.
struct Bindings {
    near: Vec<usize>,
    far: Vec<[usize; 2]>,
}

impl Bindings {
    fn scatter_near(&self, near: &mut FilamentWake, points: &[[f64; 3]]) {
        for (p, &index) in near
            .rows
            .iter_mut()
            .flat_map(|r| &mut r.nodes)
            .zip(&self.near)
        {
            *p = points[index];
        }
    }
}

fn node_key(point: [f64; 3]) -> [u64; 3] {
    // Finite geometry is admitted before indexing. Signed zero is a single
    // physical location, while all other coordinates keep their exact bits.
    point.map(|v| if v == 0.0 { 0 } else { v.to_bits() })
}

fn shared_nodes(
    near: &FilamentWake,
    far: &FarField,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Vec<[f64; 3]>, Bindings), Refusal> {
    let mut indices = BTreeMap::new();
    let mut points = Vec::new();
    let near_count = near.rows.len() * (near.cert.n_stations + 1);
    // The earlier checked source count already bounds this sum.
    let far_count: usize = far.cells.iter().map(|cell| cell.segs.len()).sum();
    let point_bound = far_count
        .checked_mul(2)
        .and_then(|n| n.checked_add(near_count))
        .ok_or_else(work_refusal)?;
    points
        .try_reserve_exact(point_bound)
        .map_err(|_| refusal("wake-advection-allocation", "shared-node allocation failed"))?;
    let mut bind = |p| {
        *indices.entry(node_key(p)).or_insert_with(|| {
            let index = points.len();
            points.push(p);
            index
        })
    };
    let mut bindings = Bindings {
        near: Vec::new(),
        far: Vec::new(),
    };
    bindings.near.try_reserve_exact(near_count).map_err(|_| {
        refusal(
            "wake-advection-allocation",
            "near-node bindings allocation failed",
        )
    })?;
    bindings.far.try_reserve_exact(far_count).map_err(|_| {
        refusal(
            "wake-advection-allocation",
            "far-node bindings allocation failed",
        )
    })?;
    for row in &near.rows {
        checkpoint(cancelled)?;
        bindings.near.extend(row.nodes.iter().map(|&p| bind(p)));
    }
    for cell in &far.cells {
        for (index, &(a, b, _)) in cell.segs.iter().enumerate() {
            if index % 256 == 0 {
                checkpoint(cancelled)?;
            }
            bindings.far.push([bind(a), bind(b)]);
        }
    }
    Ok((points, bindings))
}

fn admit_sources(
    far: &FarField,
    near: &FilamentWake,
    growth: f64,
    max_work: u64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<usize, Refusal> {
    let mut count = near.rows.len() * (2 * near.cert.n_stations + 1) + near.cert.n_stations;
    for cell in &far.cells {
        checkpoint(cancelled)?;
        // Include cell traversal even when it contains no active sources.
        count = count
            .checked_add(cell.segs.len())
            .and_then(|n| n.checked_add(1))
            .ok_or_else(work_refusal)?;
        // Bound even the admission/indexing passes before allocating scratch.
        // The common stage engine additionally admits its full interaction cost.
        if u64::try_from(count).ok().is_none_or(|n| n > max_work) {
            return Err(work_refusal());
        }
        let next_offset = cell.diffusion_offset_m2 + growth;
        if cell.segs.len() != cell.source_core2_m2.len()
            || !cell.diffusion_offset_m2.is_finite()
            || cell.diffusion_offset_m2 < 0.0
            || !next_offset.is_finite()
        {
            return Err(refusal(
                "wake-advection-core-invalid",
                "invalid far-source cores or diffusion offset",
            ));
        }
        for (index, ((a, b, gamma), base_core)) in
            cell.segs.iter().zip(&cell.source_core2_m2).enumerate()
        {
            if index % 256 == 0 {
                checkpoint(cancelled)?;
            }
            let core = base_core + cell.diffusion_offset_m2;
            let final_core = base_core + next_offset;
            if !base_core.is_finite()
                || *base_core < 0.0
                || !core.is_finite()
                || (*gamma != 0.0 && core <= 0.0)
                || !final_core.is_finite()
                || final_core > f64::MAX / 16.0
                || !gamma.is_finite()
                || a.iter().chain(b).any(|v| !v.is_finite())
            {
                return Err(refusal(
                    "wake-advection-core-invalid",
                    "far sources need finite geometry, circulation and admissible evolved cores",
                ));
            }
        }
    }
    if u64::try_from(count).ok().is_none_or(|n| n > max_work) {
        return Err(work_refusal());
    }
    Ok(count)
}

#[derive(Default)]
struct Balance {
    sum: f64,
    correction: f64,
    magnitude: f64,
}

impl Balance {
    fn add(&mut self, value: f64) {
        let next = self.sum + value;
        self.correction += if self.sum.abs() >= value.abs() {
            (self.sum - next) + value
        } else {
            (value - next) + self.sum
        };
        self.sum = next;
        self.magnitude += value.abs();
    }

    fn closed(&self) -> bool {
        self.magnitude.is_finite()
            && (self.sum + self.correction).abs() <= 64.0 * f64::EPSILON * self.magnitude
    }
}

fn admit_circulation(
    near: &FilamentWake,
    far: &FarField,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), Refusal> {
    let mut balances = BTreeMap::<[u64; 3], Balance>::new();
    let mut count = 0_usize;
    let mut failure = None;
    let mut add = |a, b, gamma: f64| {
        if failure.is_some() {
            return;
        }
        if count % 256 == 0 {
            if let Err(error) = checkpoint(cancelled) {
                failure = Some(error);
                return;
            }
        }
        count += 1;
        balances.entry(node_key(a)).or_default().add(-gamma);
        balances.entry(node_key(b)).or_default().add(gamma);
    };
    for cell in &far.cells {
        for &(a, b, gamma) in &cell.segs {
            add(a, b, gamma);
        }
    }
    near.visit_segments_with_core_growth(0.0, |source| {
        add(source.start, source.end, source.circulation);
    });
    if let Some(error) = failure {
        return Err(error);
    }
    if balances.values().any(|value| !value.closed()) {
        return Err(refusal(
            "wake-advection-topology-invalid",
            "retained near/far sources have unbalanced vortex endpoints",
        ));
    }
    Ok(())
}

fn append_far_sources(
    far: &FarField,
    bindings: &Bindings,
    points: &[[f64; 3]],
    growth: f64,
    fraction: f64,
    sources: &mut PreparedSources,
) -> Result<(), Refusal> {
    let mut indices = bindings.far.iter();
    for cell in &far.cells {
        for ((_, _, gamma), base_core) in cell.segs.iter().zip(&cell.source_core2_m2) {
            let &[a, b] = indices.next().expect("one binding per admitted source");
            if *gamma == 0.0 {
                continue;
            }
            sources.push(
                points[a],
                points[b],
                *gamma,
                (base_core + cell.diffusion_offset_m2) + growth * fraction,
            )?;
        }
    }
    Ok(())
}

fn rebuild_cells(
    far: &FarField,
    bindings: &Bindings,
    points: &[[f64; 3]],
    growth: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<FarCell>, Refusal> {
    let mut indices = bindings.far.iter();
    let mut cells = Vec::new();
    cells
        .try_reserve_exact(far.cells.len())
        .map_err(|_| refusal("wake-advection-allocation", "far-cell allocation failed"))?;
    for original in &far.cells {
        checkpoint(cancelled)?;
        let mut segs = Vec::new();
        segs.try_reserve_exact(original.segs.len())
            .map_err(|_| refusal("wake-advection-allocation", "far-source allocation failed"))?;
        for (index, &(_, _, gamma)) in original.segs.iter().enumerate() {
            if index % 256 == 0 {
                checkpoint(cancelled)?;
            }
            let &[a, b] = indices.next().expect("one binding per admitted source");
            segs.push((points[a], points[b], gamma));
        }
        let mut cell = cell_from_segs_with_cores(segs, original.source_core2_m2.clone());
        // A common core offset cancels in every core-deviation moment. Keep
        // base sources and elapsed diffusion separate through reconstruction.
        cell.diffusion_offset_m2 = original.diffusion_offset_m2 + growth;
        cell.core2_m2 += cell.diffusion_offset_m2;
        let scalars = [
            cell.radius_m,
            cell.w_abs_sum,
            cell.core2_m2,
            cell.max_abs_core_delta_m2,
            cell.diffusion_offset_m2,
            cell.max_source_core2_m2,
        ];
        if cell
            .centroid
            .iter()
            .chain(&cell.w_sum)
            .chain(&cell.core_moment)
            .chain(&scalars)
            .chain(cell.dip.iter().flatten())
            .chain(cell.quad.iter().flatten().flatten())
            .any(|v| !v.is_finite())
        {
            return Err(refusal(
                "wake-advection-moment-invalid",
                "moved far-source moments are not representable",
            ));
        }
        cells.push(cell);
    }
    Ok(cells)
}

fn checkpoint(cancelled: &mut impl FnMut() -> bool) -> Result<(), Refusal> {
    if cancelled() {
        Err(refusal(
            "wake-advection-cancelled",
            "wake evolution cancelled",
        ))
    } else {
        Ok(())
    }
}

fn work_refusal() -> Refusal {
    refusal(
        "wake-advection-work-exhausted",
        "hybrid sources exceed the work cap",
    )
}

fn refusal(code: &'static str, message: &str) -> Refusal {
    Refusal {
        code,
        message: message.into(),
        ranked_repairs: vec!["retain a complete admitted wake and a resolved bounded step".into()],
    }
}
