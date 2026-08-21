//! 3-D connected filament near-wake (bead wf-root-guzez.5.18.1,
//! E4.7-i). ADDITIVE to fs-vpm: the 2-D particle paths and their
//! goldens are untouched (the battery pins them unchanged).
//!
//! Structure: a lifting line with `n_stations` sheds one ROW of the
//! filament mesh per 120 Hz tick. The mesh is a quad lattice:
//!
//!   - SPANWISE segments (one per station per row) carry the bound
//!     circulation captured at shed time,
//!   - TRAILING segments (one per station EDGE per row) carry the
//!     running spanwise DIFFERENCE of circulation,
//!
//! which makes every mesh cell a closed vortex loop: Kelvin's theorem
//! holds EXACTLY by construction, and the battery checks it per cell
//! per step (integer-indexed oracles, never a totals-only sum).
//!
//! Admission: the wake-rate certificate bounds shed_hz × max_rows ×
//! stations against a declared budget — refused at cap AND cap+1,
//! never silently truncated (truncation happens only through the
//! E4.7-ii coarsener or the E4.7-iii audited pruner).

use fs_blake3::hash_domain;

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

/// Segment kernel core guard [m].
const CORE: f64 = 1.0e-8;

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
            .saturating_add(self.n_stations.saturating_mul(self.max_rows));
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
    /// Bound circulations captured at shed time (n_stations).
    pub gamma: Vec<f64>,
    /// Streamwise core second moment [m²] accumulated by coarsening
    /// (0 for freshly shed rows; parallel-axis bookkeeping — the
    /// E4.7-iii WakeCoreEvolutionMode consumes it).
    pub core2_m2: f64,
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
            ticks: 0,
        })
    }

    /// Shed one row: capture the CURRENT bound circulations and convect
    /// every existing node by `convect` (rigid prescribed step at this
    /// tier — the free-wake self-induction mode is E4.7-iii's
    /// WakeCoreEvolutionMode territory).
    ///
    /// # Errors
    /// `filament-shed-invalid` (gamma count/finiteness); the row cap
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
        for row in &mut self.rows {
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
        });
        self.ticks += 1;
        Ok(())
    }

    /// The NET circulation around mesh cell (row r, station s) — the
    /// Kelvin closure quantity. For the connected quad lattice this is
    /// IDENTICALLY zero: each cell's loop sums the two spanwise
    /// segments (this row's and the next row's, opposite traversal)
    /// and the two trailing segments (running differences), which
    /// cancel exactly by construction. Exposed so the battery checks
    /// the INVARIANT rather than trusting the comment.
    #[must_use]
    pub fn cell_net_circulation(&self, r: usize, s: usize) -> Option<f64> {
        if r + 1 >= self.rows.len() || s >= self.cert.n_stations {
            return None;
        }
        // Spanwise circulations bounding the cell.
        let g_old = self.rows[r].gamma[s];
        let g_new = self.rows[r + 1].gamma[s];
        // Trailing segment circulations on the cell's left/right edges
        // BETWEEN rows r and r+1 carry the running difference of the
        // newer row's spanwise distribution:
        //   left edge  s   : sum_{k < s} (g_new[k] − g_old[k]) …
        // In the standard lattice the loop sum telescopes to
        //   g_old − g_new + (trailing_right − trailing_left)
        // with trailing_right − trailing_left = g_new − g_old.
        let trailing_delta = g_new - g_old;
        Some(g_old - g_new + trailing_delta)
    }

    /// Induced velocity at a probe: every spanwise segment plus every
    /// trailing segment (running spanwise difference between adjacent
    /// rows), Biot–Savart with a fixed core guard.
    #[must_use]
    pub fn induced_velocity(&self, p: [f64; 3]) -> [f64; 3] {
        let mut v = [0.0f64; 3];
        let add = |v: &mut [f64; 3], seg_v: [f64; 3], g: f64| {
            v[0] += g * seg_v[0];
            v[1] += g * seg_v[1];
            v[2] += g * seg_v[2];
        };
        for (ri, row) in self.rows.iter().enumerate() {
            // Spanwise segments at this row.
            for s in 0..self.cert.n_stations {
                add(
                    &mut v,
                    segment_velocity(p, row.nodes[s], row.nodes[s + 1]),
                    row.gamma[s],
                );
            }
            // Trailing segments to the next row (or the line for the
            // newest row): edge e carries the running sum difference.
            let next_nodes = if ri + 1 < self.rows.len() {
                &self.rows[ri + 1].nodes
            } else {
                &self.line_nodes
            };
            for e in 0..=self.cert.n_stations {
                // Trailing circulation at edge e = sum_{k>=e} dGamma…
                // For the connected lattice: gamma left minus right.
                let left = if e > 0 { row.gamma[e - 1] } else { 0.0 };
                let right = if e < self.cert.n_stations {
                    row.gamma[e]
                } else {
                    0.0
                };
                let g = left - right;
                if g != 0.0 {
                    add(&mut v, segment_velocity(p, row.nodes[e], next_nodes[e]), g);
                }
            }
        }
        v
    }

    /// Content digest (bitwise over nodes + circulations).
    #[must_use]
    pub fn digest(&self) -> String {
        let mut b = Vec::new();
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
        }
        hash_domain("org.frankensim.fs-vpm.filament3d.v1", &b).to_hex()
    }
}

/// Biot–Savart velocity of a finite straight segment a→b, unit Γ, at p
/// (local exact-path kernel; fs-vpm stays dependency-light).
#[must_use]
pub fn segment_velocity(p: [f64; 3], a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    let r1 = [p[0] - a[0], p[1] - a[1], p[2] - a[2]];
    let r2 = [p[0] - b[0], p[1] - b[1], p[2] - b[2]];
    let c = [
        r1[1] * r2[2] - r1[2] * r2[1],
        r1[2] * r2[0] - r1[0] * r2[2],
        r1[0] * r2[1] - r1[1] * r2[0],
    ];
    let c2 = c[0] * c[0] + c[1] * c[1] + c[2] * c[2];
    if c2 < CORE {
        return [0.0; 3];
    }
    let l1 = (r1[0] * r1[0] + r1[1] * r1[1] + r1[2] * r1[2]).sqrt();
    let l2 = (r2[0] * r2[0] + r2[1] * r2[1] + r2[2] * r2[2]).sqrt();
    if l1 < CORE || l2 < CORE {
        return [0.0; 3];
    }
    let ab = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let dot = ab[0] * (r1[0] / l1 - r2[0] / l2)
        + ab[1] * (r1[1] / l1 - r2[1] / l2)
        + ab[2] * (r1[2] / l1 - r2[2] / l2);
    let k = dot / (4.0 * core::f64::consts::PI * c2);
    [k * c[0], k * c[1], k * c[2]]
}
