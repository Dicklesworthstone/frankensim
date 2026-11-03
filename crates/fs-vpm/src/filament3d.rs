//! 3-D connected filament near-wake (bead wf-root-guzez.5.18.1,
//! E4.7-i). ADDITIVE to fs-vpm: the 2-D particle paths and their
//! goldens are untouched (the battery pins them unchanged).
//!
//! Structure: a lifting line with `n_stations` sheds one ROW of the
//! filament mesh per 120 Hz tick. The mesh is a quad lattice:
//!
//!   - each quad carries the bound circulation captured at shed time,
//!   - SPANWISE segments carry the DIFFERENCE between adjacent quads,
//!   - TRAILING segments carry the spanwise DIFFERENCE of circulation,
//!   - the newest line closes the sheet with the current bound circulation.
//!
//! Adjacent quads share one oriented physical edge and one core. This
//! cancels steady interior spans before regularization, including when
//! adjacent rows have different ages. Vortex lines have no unbalanced
//! endpoints: the signed incident circulations close at every node.
//!
//! Admission: the wake-rate certificate bounds shed_hz × max_rows ×
//! stations against a declared budget — refused at cap AND cap+1,
//! never silently truncated (truncation happens only through the
//! E4.7-ii coarsener or the E4.7-iii audited pruner).

use fs_blake3::hash_domain;

#[path = "filament3d/advection.rs"]
pub(crate) mod advection;
pub use advection::{SelfInducedShed, WakeBoundary};

#[path = "filament3d/remesh.rs"]
mod remesh;
pub use remesh::StreamwiseRefinement;

#[path = "filament3d/impulse.rs"]
mod impulse;
pub(crate) use impulse::ring_angular_impulse;

/// Typed refusal (workspace convention; the 3-D filament lane's own —
/// the 2-D lane keeps its `VpmError`, untouched).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Refusal {
    /// Stable machine-readable code.
    pub code: &'static str,
    /// Human-readable diagnosis.
    pub message: String,
    /// Ranked repairs, most likely fix first.
    pub ranked_repairs: Vec<String>,
}

/// Wake-rate certificate budget: segments the near-wake may hold.
pub const NEAR_WAKE_SEGMENT_BUDGET: usize = 2_000_000;

/// Station cap.
pub const MAX_STATIONS: usize = 256;

/// Row cap.
pub const MAX_ROWS: usize = 4_096;

/// The wake-rate admission certificate (registered; enters identity).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WakeRateCertificate {
    /// Shedding rate [Hz].
    pub shed_hz: f64,
    /// Stations on the lifting line.
    pub n_stations: usize,
    /// Maximum retained near-wake rows.
    pub max_rows: usize,
}

impl WakeRateCertificate {
    /// Admit the certificate.
    ///
    /// # Errors
    /// `wake-rate-uncertified` (caps at cap AND cap+1 on stations,
    /// rows, rate, and the segment budget).
    pub fn admit(&self) -> Result<(), Refusal> {
        let segments = self
            .n_stations
            .saturating_add(1)
            .saturating_mul(self.max_rows)
            .saturating_add(self.n_stations.saturating_mul(self.max_rows))
            .saturating_add(self.n_stations);
        let ok = self.shed_hz.is_finite()
            && self.shed_hz > 0.0
            && self.shed_hz <= 1_000.0
            && (1..=MAX_STATIONS).contains(&self.n_stations)
            && (1..=MAX_ROWS).contains(&self.max_rows)
            && segments <= NEAR_WAKE_SEGMENT_BUDGET;
        if ok {
            Ok(())
        } else {
            Err(Refusal {
                code: "wake-rate-uncertified",
                message: format!("{self:?} (segments {segments})"),
                ranked_repairs: vec![format!(
                    "shed_hz (0, 1000]; stations [1, {MAX_STATIONS}]; rows [1, {MAX_ROWS}]; \
                     segments <= {NEAR_WAKE_SEGMENT_BUDGET}"
                )],
            })
        }
    }
}

/// One shed row of the mesh.
#[derive(Clone, Debug, PartialEq)]
pub struct ShedRow {
    /// Node positions at this row's spanwise line (n_stations + 1).
    pub nodes: Vec<[f64; 3]>,
    /// Bound circulations captured at shed time (n_stations). Each value
    /// is the strength of the quad from this row to its newer successor,
    /// positive along the successor's left-to-right spanwise edge.
    pub gamma: Vec<f64>,
    /// Base core width squared [m²], including the spatial geometric
    /// smoothing added by coarsening. Zero for newly shed rows; viscous
    /// growth is applied separately using retained elapsed age.
    pub core2_m2: f64,
    /// Elapsed shedding intervals since emission; fractional after coarsening.
    /// A newly shed row is age zero. Conversion never resets this age.
    pub age_ticks: f64,
}

/// One oriented physical filament after adjacent panel contributions cancel.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WakeSegment {
    /// Start of the oriented segment [m].
    pub start: [f64; 3],
    /// End of the oriented segment [m].
    pub end: [f64; 3],
    /// Circulation along start → end [m²/s].
    pub circulation: f64,
    /// Squared Rosenhead–Moore core radius of this edge [m²].
    pub core2_m2: f64,
}

/// Velocity and spatial derivatives of the represented physical wake sources.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct WakeKinematics {
    /// Induced velocity [m/s], including signed source circulations.
    pub velocity: [f64; 3],
    /// `gradient[i][j] = ∂velocity[i]/∂position[j]` [1/s].
    pub gradient: [[f64; 3]; 3],
}

impl WakeKinematics {
    /// Physical vorticity `curl(velocity)` [1/s].
    #[must_use]
    pub fn vorticity(&self) -> [f64; 3] {
        let g = &self.gradient;
        [g[2][1] - g[1][2], g[0][2] - g[2][0], g[1][0] - g[0][1]]
    }

    /// Velocity divergence [1/s]; analytically zero for these source fields.
    #[must_use]
    pub fn divergence(&self) -> f64 {
        self.gradient[0][0] + self.gradient[1][1] + self.gradient[2][2]
    }
}

/// The connected near-wake.
#[derive(Clone, Debug, PartialEq)]
pub struct FilamentWake {
    /// Certificate this wake runs under.
    pub cert: WakeRateCertificate,
    /// Lifting-line node positions (n_stations + 1), newest edge.
    pub line_nodes: Vec<[f64; 3]>,
    /// Shed rows, oldest first.
    pub rows: Vec<ShedRow>,
    /// Captured quad strengths immediately downstream of the oldest retained
    /// row, if those quads have been transferred to a far field. Their trailing
    /// edges end at this row; retaining the strengths preserves its shared
    /// spanwise edge. `None` denotes the physical starting edge of a full wake.
    pub downstream_gamma: Option<Vec<f64>>,
    /// Ticks shed so far.
    pub ticks: u64,
}

impl FilamentWake {
    /// New wake for a lifting line (nodes = station edges).
    ///
    /// # Errors
    /// Certificate refusals; `filament-line-invalid` (node count must
    /// be stations + 1 and finite).
    pub fn new(cert: WakeRateCertificate, line_nodes: Vec<[f64; 3]>) -> Result<Self, Refusal> {
        cert.admit()?;
        if line_nodes.len() != cert.n_stations + 1
            || line_nodes.iter().flatten().any(|v| !v.is_finite())
        {
            return Err(Refusal {
                code: "filament-line-invalid",
                message: format!(
                    "{} nodes for {} stations",
                    line_nodes.len(),
                    cert.n_stations
                ),
                ranked_repairs: vec!["nodes = stations + 1, all finite".into()],
            });
        }
        Ok(FilamentWake {
            cert,
            line_nodes,
            rows: Vec::new(),
            downstream_gamma: None,
            ticks: 0,
        })
    }

    /// Shed one row: capture the CURRENT bound circulations and convect
    /// every existing node by `convect` (rigid prescribed step). Use
    /// [`Self::shed_self_induced`] to resolve motion induced by the wake.
    ///
    /// # Errors
    /// `filament-shed-invalid` (gamma count/finiteness),
    /// `filament-convection-invalid` (non-finite geometry/displacement or
    /// translated coordinate overflow); the row cap
    /// DROPS the oldest row only through the caller's coarsener —
    /// here it refuses (`filament-rows-exhausted`), never silently.
    pub fn shed(&mut self, gamma: &[f64], convect: [f64; 3]) -> Result<(), Refusal> {
        if gamma.len() != self.cert.n_stations || gamma.iter().any(|g| !g.is_finite()) {
            return Err(Refusal {
                code: "filament-shed-invalid",
                message: format!(
                    "{} gammas for {} stations",
                    gamma.len(),
                    self.cert.n_stations
                ),
                ranked_repairs: vec!["one finite circulation per station".into()],
            });
        }
        if self.rows.len() >= self.cert.max_rows {
            return Err(Refusal {
                code: "filament-rows-exhausted",
                message: format!("{} rows at the certificate cap", self.rows.len()),
                ranked_repairs: vec![
                    "coarsen (E4.7-ii) or prune with audit (E4.7-iii) — never drop silently".into(),
                ],
            });
        }
        let next_tick = self.ticks.checked_add(1).ok_or_else(|| Refusal {
            code: "filament-age-invalid",
            message: "the shedding tick counter is exhausted".into(),
            ranked_repairs: vec!["start a new explicitly timed wake interval".into()],
        })?;
        if self.rows.iter().any(|row| {
            !row.age_ticks.is_finite()
                || row.age_ticks < 0.0
                || !(row.age_ticks + 1.0).is_finite()
                || row.age_ticks + 1.0 == row.age_ticks
        }) {
            return Err(Refusal {
                code: "filament-age-invalid",
                message: "row age must admit one further finite shedding interval".into(),
                ranked_repairs: vec![
                    "retain finite nonnegative ages with adequate time resolution".into(),
                ],
            });
        }
        if convect.iter().any(|v| !v.is_finite())
            || self.line_nodes.iter().flatten().any(|v| !v.is_finite())
            || self
                .rows
                .iter()
                .flat_map(|r| &r.nodes)
                .any(|p| (0..3).any(|i| !(p[i] + convect[i]).is_finite()))
        {
            return Err(Refusal {
                code: "filament-convection-invalid",
                message: "prescribed convection must retain finite wake coordinates".into(),
                ranked_repairs: vec!["use finite geometry and a representable displacement".into()],
            });
        }
        for row in &mut self.rows {
            row.age_ticks += 1.0;
            for n in &mut row.nodes {
                n[0] += convect[0];
                n[1] += convect[1];
                n[2] += convect[2];
            }
        }
        self.rows.push(ShedRow {
            nodes: self.line_nodes.clone(),
            gamma: gamma.to_vec(),
            core2_m2: 0.0,
            age_ticks: 0.0,
        });
        self.ticks = next_tick;
        Ok(())
    }

    /// Signed incident-circulation residual at the newer, right-hand vertex
    /// of cell `(r, s)`. This audits the same assembled edge strengths used
    /// for induction, rather than a difference that cancels itself.
    /// Use [`Self::node_circulation_residual`] to include outer vertices.
    #[must_use]
    pub fn cell_net_circulation(&self, r: usize, s: usize) -> Option<f64> {
        let newer = r.checked_add(1)?;
        if newer >= self.rows.len() || s >= self.cert.n_stations {
            return None;
        }
        self.node_circulation_residual(newer, s + 1)
    }

    /// Incoming minus outgoing circulation at one lattice node. Rows run
    /// `0..=rows.len()`, with the final index denoting the lifting line;
    /// spanwise node indices run `0..=n_stations`. All four incident edge
    /// strengths use the same differences as velocity evaluation. Roundoff
    /// residuals are possible for non-binary circulation values.
    ///
    /// After aggregation, incoming trailing edges at row zero belong to the
    /// far field; their retained strengths are included in this local balance.
    #[must_use]
    pub fn node_circulation_residual(&self, row: usize, edge: usize) -> Option<f64> {
        if row > self.rows.len() || edge > self.cert.n_stations {
            return None;
        }
        if self.rows.is_empty() {
            return Some(0.0);
        }
        let incoming_span = if edge > 0 {
            self.spanwise_circulation(row, edge - 1)
        } else {
            0.0
        };
        let outgoing_span = if edge < self.cert.n_stations {
            self.spanwise_circulation(row, edge)
        } else {
            0.0
        };
        let incoming_trailing = if row > 0 {
            trailing_circulation(&self.rows[row - 1].gamma, edge)
        } else {
            self.downstream_gamma
                .as_ref()
                .map_or(0.0, |gamma| trailing_circulation(gamma, edge))
        };
        let outgoing_trailing = self
            .rows
            .get(row)
            .map_or(0.0, |r| trailing_circulation(&r.gamma, edge));
        Some((incoming_span - outgoing_span) + (incoming_trailing - outgoing_trailing))
    }

    fn spanwise_circulation(&self, row: usize, station: usize) -> f64 {
        if row == self.rows.len() {
            return self.rows.last().map_or(0.0, |r| r.gamma[station]);
        }
        let previous = if row > 0 {
            self.rows[row - 1].gamma[station]
        } else {
            self.downstream_gamma.as_ref().map_or(0.0, |g| g[station])
        };
        previous - self.rows[row].gamma[station]
    }

    /// Visit each nonzero oriented physical edge once, without allocation.
    /// A spanwise edge uses the newer adjacent row's retained core and age;
    /// trailing edges use their emitting row. The closing line uses the newest
    /// row's core. Steady interior spans cancel before core evaluation.
    ///
    /// This includes the line's closing segment. A coupled lifting-line solver
    /// must not add that same bound segment a second time. After aggregation,
    /// callers add the far field, which owns the removed rows' trailing edges.
    /// As with the velocity helpers, inputs must already be finite/admitted.
    pub fn visit_segments_with_core_growth(
        &self,
        growth_per_tick_m2: f64,
        mut visit: impl FnMut(WakeSegment),
    ) {
        for row in 0..self.rows.len() {
            self.visit_row_segments(row, growth_per_tick_m2, &mut visit);
        }
        if let Some(last) = self.rows.last() {
            let core2_m2 = evolved_core(last, growth_per_tick_m2);
            for s in 0..self.cert.n_stations {
                let circulation = last.gamma[s];
                if circulation != 0.0 {
                    visit(WakeSegment {
                        start: self.line_nodes[s],
                        end: self.line_nodes[s + 1],
                        circulation,
                        core2_m2,
                    });
                }
            }
        }
    }

    pub(crate) fn visit_row_segments(
        &self,
        row_index: usize,
        growth_per_tick_m2: f64,
        visit: &mut impl FnMut(WakeSegment),
    ) {
        let row = &self.rows[row_index];
        let core2_m2 = evolved_core(row, growth_per_tick_m2);
        for s in 0..self.cert.n_stations {
            let circulation = self.spanwise_circulation(row_index, s);
            if circulation != 0.0 {
                visit(WakeSegment {
                    start: row.nodes[s],
                    end: row.nodes[s + 1],
                    circulation,
                    core2_m2,
                });
            }
        }
        let next_nodes = self
            .rows
            .get(row_index + 1)
            .map_or(&self.line_nodes, |r| &r.nodes);
        for edge in 0..=self.cert.n_stations {
            let circulation = trailing_circulation(&row.gamma, edge);
            if circulation != 0.0 {
                visit(WakeSegment {
                    start: row.nodes[edge],
                    end: next_nodes[edge],
                    circulation,
                    core2_m2,
                });
            }
        }
    }

    /// Induced velocity of the complete assembled sheet, including its closing
    /// line, using each physical edge's finite Rosenhead–Moore core.
    #[must_use]
    pub fn induced_velocity(&self, p: [f64; 3]) -> [f64; 3] {
        self.induced_velocity_with_core_growth(p, 0.0)
    }

    /// Evaluate stored row cores plus `growth_per_tick_m2 * age_ticks`.
    /// Growth and retained ages must be finite/nonnegative; the checked
    /// aggregate/image entry points validate their complete input state.
    #[must_use]
    pub fn induced_velocity_with_core_growth(
        &self,
        p: [f64; 3],
        growth_per_tick_m2: f64,
    ) -> [f64; 3] {
        let mut v = [0.0f64; 3];
        self.visit_segments_with_core_growth(growth_per_tick_m2, |segment| {
            let sv = segment_velocity_with_core(p, segment.start, segment.end, segment.core2_m2);
            for i in 0..3 {
                v[i] += segment.circulation * sv[i];
            }
        });
        v
    }

    /// Evaluate physical velocity and its analytic spatial Jacobian using
    /// exactly the same assembled edges and cores as the velocity-only path.
    /// No finite-difference length, source allocation, or extra probe sweeps
    /// are needed. At a finite-core axis the curl remains resolved even when
    /// that segment's velocity is zero.
    ///
    /// `None` denotes an invalid probe/growth, a singular zero-core source at
    /// the probe, or an unrepresentable sampled field. The stored lattice must
    /// already be admitted, as for the velocity-only helpers. After aggregation
    /// this samples retained near sources; it does not differentiate far cells.
    #[must_use]
    pub fn induced_kinematics_with_core_growth(
        &self,
        p: [f64; 3],
        growth_per_tick_m2: f64,
    ) -> Option<WakeKinematics> {
        if !growth_per_tick_m2.is_finite()
            || growth_per_tick_m2 < 0.0
            || p.iter().any(|v| !v.is_finite())
        {
            return None;
        }
        let mut result = WakeKinematics::default();
        let mut regular = true;
        self.visit_segments_with_core_growth(growth_per_tick_m2, |s| {
            if !regular {
                return;
            }
            let Some(value) = segment_kinematics_with_core(p, s.start, s.end, s.core2_m2) else {
                regular = false;
                return;
            };
            for i in 0..3 {
                result.velocity[i] += s.circulation * value.velocity[i];
                for j in 0..3 {
                    result.gradient[i][j] += s.circulation * value.gradient[i][j];
                }
            }
        });
        (regular
            && result
                .velocity
                .iter()
                .chain(result.gradient.iter().flatten())
                .chain(result.vorticity().iter())
                .all(|x| x.is_finite())
            && result.divergence().is_finite())
        .then_some(result)
    }

    /// Geometric hydrodynamic impulse divided by uniform fluid density [m⁴/s].
    /// Each closed ring contributes Γ times its oriented vector area, equal to
    /// `½ Γ ∮ x × dx`. This is invariant under translating the world origin;
    /// coordinates are differenced before cross products. Positive circulation
    /// follows the newer row's left-to-right edge, as in source assembly.
    /// This geometric diagnostic does not represent angular momentum or the
    /// momentum of fluid inside an arbitrary finite control volume.
    ///
    /// # Errors
    /// `impulse-interface-attached` when rows have been transferred to a far
    /// field: only the combined closed near/far sources have an origin-independent
    /// impulse. Invalid geometry or an unrepresentable result also refuses.
    pub fn hydrodynamic_impulse_per_density(&self) -> Result<[f64; 3], Refusal> {
        if self.downstream_gamma.is_some() {
            return Err(Refusal {
                code: "impulse-interface-attached",
                message: "retained near sources alone do not form a complete closed wake".into(),
                ranked_repairs: vec!["measure the complete wake before near/far partition".into()],
            });
        }
        self.admit_geometry()?;
        let mut impulse = [0.0; 3];
        let mut correction = [0.0; 3];
        for (ri, row) in self.rows.iter().enumerate() {
            let next = self
                .rows
                .get(ri + 1)
                .map_or(&self.line_nodes, |row| &row.nodes);
            for s in 0..self.cert.n_stations {
                let area = ring_vector_area(next[s], next[s + 1], row.nodes[s + 1], row.nodes[s]);
                for c in 0..3 {
                    let value = row.gamma[s] * area[c] - correction[c];
                    let sum = impulse[c] + value;
                    correction[c] = (sum - impulse[c]) - value;
                    impulse[c] = sum;
                }
            }
        }
        if impulse.iter().any(|v| !v.is_finite()) {
            return Err(Refusal {
                code: "impulse-unrepresentable",
                message: "circulation times vector area exceeds the numerical range".into(),
                ranked_repairs: vec![
                    "rescale geometry and circulation to representable units".into(),
                ],
            });
        }
        Ok(impulse)
    }

    pub(crate) fn admit_geometry(&self) -> Result<(), Refusal> {
        let n = self.cert.n_stations;
        if !(1..=MAX_STATIONS).contains(&n)
            || self.line_nodes.len() != n + 1
            || self.line_nodes.iter().flatten().any(|v| !v.is_finite())
            || self.rows.iter().any(|row| {
                row.nodes.len() != n + 1
                    || row.gamma.len() != n
                    || row.nodes.iter().flatten().any(|v| !v.is_finite())
                    || row.gamma.iter().any(|v| !v.is_finite())
            })
        {
            return Err(Refusal {
                code: "filament-geometry-invalid",
                message: "wake geometry and circulation must match the finite station lattice"
                    .into(),
                ranked_repairs: vec![
                    "retain stations + 1 finite nodes and one finite gamma per station".into(),
                ],
            });
        }
        Ok(())
    }

    /// Content digest over physical sheet geometry, quad strengths, ages,
    /// cores, and the retained near/far interface circulation.
    #[must_use]
    pub fn digest(&self) -> String {
        let mut b = Vec::new();
        b.extend_from_slice(&(self.cert.n_stations as u64).to_le_bytes());
        for node in &self.line_nodes {
            for c in node {
                b.extend_from_slice(&c.to_bits().to_le_bytes());
            }
        }
        b.push(u8::from(self.downstream_gamma.is_some()));
        if let Some(gamma) = &self.downstream_gamma {
            for g in gamma {
                b.extend_from_slice(&g.to_bits().to_le_bytes());
            }
        }
        for row in &self.rows {
            for n in &row.nodes {
                for c in n {
                    b.extend_from_slice(&c.to_bits().to_le_bytes());
                }
            }
            for g in &row.gamma {
                b.extend_from_slice(&g.to_bits().to_le_bytes());
            }
            b.extend_from_slice(&row.core2_m2.to_bits().to_le_bytes());
            b.extend_from_slice(&row.age_ticks.to_bits().to_le_bytes());
        }
        hash_domain("org.frankensim.fs-vpm.filament3d.v3", &b).to_hex()
    }
}

fn trailing_circulation(gamma: &[f64], edge: usize) -> f64 {
    let left = if edge > 0 { gamma[edge - 1] } else { 0.0 };
    let right = gamma.get(edge).copied().unwrap_or(0.0);
    right - left
}

fn evolved_core(row: &ShedRow, growth_per_tick_m2: f64) -> f64 {
    if growth_per_tick_m2 == 0.0 {
        row.core2_m2
    } else {
        row.core2_m2 + growth_per_tick_m2 * row.age_ticks
    }
}

pub(crate) fn ring_vector_area(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> [f64; 3] {
    // Half the cross product of the oriented diagonals equals the polygon
    // edge integral, including nonplanar quads. It avoids a preferred corner
    // triangulation, preserving mirrored parallel sheets to the last bit.
    let c = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
    let d = [d[0] - b[0], d[1] - b[1], d[2] - b[2]];
    [
        0.5 * (c[1] * d[2] - c[2] * d[1]),
        0.5 * (c[2] * d[0] - c[0] * d[2]),
        0.5 * (c[0] * d[1] - c[1] * d[0]),
    ]
}

/// Biot–Savart velocity of a finite straight segment a→b, unit Γ, at p.
/// Coincident endpoints and exactly axial points return zero.
#[must_use]
pub fn segment_velocity(p: [f64; 3], a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    segment_velocity_with_core(p, a, b, 0.0)
}

#[path = "filament3d/kinematics.rs"]
mod kinematics;
pub use kinematics::{SegmentKinematics, segment_kinematics_with_core};

/// Analytic finite-segment integral of the Rosenhead–Moore kernel
/// `dl × r / (4π (|r|² + core2_m2)^(3/2))`, for unit circulation.
///
/// `core2_m2` must be a finite nonnegative regularization length squared [m²].
/// An infinite straight filament has speed `r / (2π (r² + core2_m2))`.
/// Growing this parameter as `a² = a₀² + 4νt` models viscous core spreading;
/// the algebraic profile is a reduced model, not exact Lamb–Oseen diffusion.
/// No absolute length cutoff is used: scaling all lengths by `s` scales
/// velocity by `1/s` at fixed circulation. Exactly axial or degenerate
/// segments contribute zero, including the zero-core singular convention.
#[must_use]
pub fn segment_velocity_with_core(
    p: [f64; 3],
    a: [f64; 3],
    b: [f64; 3],
    core2_m2: f64,
) -> [f64; 3] {
    SegmentGeometry::new(a, b).velocity(p, core2_m2)
}

/// Immutable source geometry shared by one-off queries and cached RK stages.
/// A stage computes the source length and unit tangent once, not per target.
#[derive(Clone, Copy)]
struct SegmentGeometry {
    start: [f64; 3],
    end: [f64; 3],
    tangent: [f64; 3],
    length: f64,
}

impl SegmentGeometry {
    #[inline]
    fn new(a: [f64; 3], b: [f64; 3]) -> Self {
        let ab = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
        let length = (ab[0] * ab[0] + ab[1] * ab[1] + ab[2] * ab[2]).sqrt();
        let tangent = if length == 0.0 {
            [0.0; 3]
        } else {
            ab.map(|v| v / length)
        };
        Self {
            start: a,
            end: b,
            tangent,
            length,
        }
    }

    #[inline]
    fn velocity(&self, p: [f64; 3], core2_m2: f64) -> [f64; 3] {
        let Self {
            start: a,
            end: b,
            tangent: t,
            length,
        } = *self;
        if length == 0.0 {
            return [0.0; 3];
        }
        let r1 = [p[0] - a[0], p[1] - a[1], p[2] - a[2]];
        let r2 = [p[0] - b[0], p[1] - b[1], p[2] - b[2]];
        let c = [
            t[1] * r1[2] - t[2] * r1[1],
            t[2] * r1[0] - t[0] * r1[2],
            t[0] * r1[1] - t[1] * r1[0],
        ];
        let h2 = c[0] * c[0] + c[1] * c[1] + c[2] * c[2];
        if h2 == 0.0 {
            return [0.0; 3];
        }
        let z1 = t[0] * r1[0] + t[1] * r1[1] + t[2] * r1[2];
        let z2 = t[0] * r2[0] + t[1] * r2[1] + t[2] * r2[2];
        let transverse2 = h2 + core2_m2;
        let l1 = (z1 * z1 + transverse2).sqrt();
        let l2 = (z2 * z2 + transverse2).sqrt();
        // Work with a unit tangent: all squared lengths have units m²,
        // avoiding the old dimensionally inconsistent cross-product cutoff.
        let coefficient = if (z1 > 0.0 && z2 > 0.0) || (z1 < 0.0 && z2 < 0.0) {
            // Rationalize the endpoint difference outside the segment. Directly
            // subtracting z1/l1 - z2/l2 loses the small finite field near the axis
            // or far beyond an endpoint, especially for a short segment.
            length / (l1 * l2) * ((z1 + z2) / (z1 * l2 + z2 * l1))
        } else {
            (z1 / l1 - z2 / l2) / transverse2
        };
        let k = coefficient / (4.0 * core::f64::consts::PI);
        [k * c[0], k * c[1], k * c[2]]
    }
}
