//! TOPOLOGY CONSTRAINTS VIA PERSISTENCE (plan §9.5/§7.8, bead 7tv.15;
//! [M] — behind `moonshot-topo-persistence` per the ambition-tag
//! discipline): Betti targets enforced through DIAGRAM PENALTIES — the
//! principled way to tell an optimizer "no enclosed voids" or "exactly
//! one tunnel" without hand-crafted heuristics.
//!
//! The two localizable penalty channels of this tier:
//! - COMPONENTS (H₀): excess components penalized by their sublevel
//!   PERSISTENCE (from [`crate::cubical::persistence0`]'s elder-rule
//!   bars); attribution = the voxels of the weakest excess component.
//! - ENCLOSED VOIDS (H₂, the castability constraint): a void of the
//!   solid is a connected component of the EMPTY phase that never
//!   touches the domain boundary (the mold cannot escape); its
//!   persistence is its DEPTH (the margin by which it is empty), and
//!   its attribution is exactly its voxel set — fill there and the
//!   penalty falls. This duality route localizes H₂ with union-find
//!   instead of a full cubical boundary matrix.
//!
//! TUNNELS (H₁) are counted (via [`crate::cubical::betti`]) and
//! penalized by deficit/excess, but NOT localized in this tier — the
//! full cubical PH matrix is the growth path (CONTRACT no-claim).
//!
//! The persistence threshold τ IS the physical feature-size floor:
//! bars/voids shallower than τ are noise (a 2-voxel dimple is not a
//! void), exactly the length-scale tie the bead demands.

use crate::cubical::{Bar, VoxelField, betti, persistence0};

/// The topology target vocabulary.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TopoSpec {
    /// Required connected components of the solid (usually 1).
    pub components: u32,
    /// Required tunnels (H₁ target; e.g. 1 for a routing duct).
    pub tunnels: u32,
    /// Allowed enclosed voids (0 for castability).
    pub enclosed_voids: u32,
    /// The persistence / depth threshold — the feature-size floor in
    /// field units.
    pub tau: f64,
    /// The occupancy level (solid = value < level).
    pub level: f64,
}

/// One localized violation: a penalty contribution plus WHERE the
/// field must change.
#[derive(Debug, Clone)]
pub struct Attribution {
    /// Which constraint channel produced this.
    pub channel: &'static str,
    /// Penalty mass contributed.
    pub amount: f64,
    /// Voxel indices where changing the field reduces the penalty
    /// (fill for voids, empty for excess components).
    pub voxels: Vec<usize>,
    /// The signed direction: +1 = increase density/occupancy (fill),
    /// −1 = decrease (carve).
    pub direction: f64,
}

/// The full penalty report: graded magnitude + attribution maps +
/// the diagram evidence.
#[derive(Debug, Clone)]
pub struct TopoPenalty {
    /// Total penalty (0 iff the diagram matches the target up to τ).
    pub total: f64,
    /// Localized contributions.
    pub attributions: Vec<Attribution>,
    /// The H₀ bars (diagram evidence, ledger-ready).
    pub bars0: Vec<Bar>,
    /// Betti numbers at the level.
    pub betti: (u32, u32, u32),
    /// Enclosed-void depths (each > τ).
    pub void_depths: Vec<f64>,
}

/// Union-find with deterministic tie-breaking.
struct Uf {
    parent: Vec<u32>,
}

impl Uf {
    fn new(n: usize) -> Uf {
        Uf {
            parent: (0..n as u32).collect(),
        }
    }

    fn find(&mut self, mut x: u32) -> u32 {
        while self.parent[x as usize] != x {
            self.parent[x as usize] = self.parent[self.parent[x as usize] as usize];
            x = self.parent[x as usize];
        }
        x
    }

    fn union(&mut self, a: u32, b: u32) {
        let (ra, rb) = (self.find(a), self.find(b));
        // Deterministic: smaller root wins.
        if ra < rb {
            self.parent[rb as usize] = ra;
        } else if rb < ra {
            self.parent[ra as usize] = rb;
        }
    }
}

/// Enclosed voids of the solid at `level`: connected components of the
/// EMPTY phase (value ≥ level) that do not touch the domain boundary.
/// Returns per-void (depth, voxel list), depth = min over the void of
/// (value − level) — the margin that must be erased to fill it.
#[must_use]
pub fn enclosed_voids(field: &VoxelField, level: f64) -> Vec<(f64, Vec<usize>)> {
    let [nx, ny, nz] = field.dims;
    let n = field.values.len();
    let empty: Vec<bool> = field.values.iter().map(|&v| v >= level).collect();
    let mut uf = Uf::new(n);
    let idx = |x: u32, y: u32, z: u32| -> usize { ((z * ny + y) * nx + x) as usize };
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                let i = idx(x, y, z);
                if !empty[i] {
                    continue;
                }
                if x + 1 < nx && empty[idx(x + 1, y, z)] {
                    uf.union(i as u32, idx(x + 1, y, z) as u32);
                }
                if y + 1 < ny && empty[idx(x, y + 1, z)] {
                    uf.union(i as u32, idx(x, y + 1, z) as u32);
                }
                if z + 1 < nz && empty[idx(x, y, z + 1)] {
                    uf.union(i as u32, idx(x, y, z + 1) as u32);
                }
            }
        }
    }
    // Mark roots whose components touch the boundary (mold escape).
    let mut escapes = vec![false; n];
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                if (x == 0 || y == 0 || z == 0 || x == nx - 1 || y == ny - 1 || z == nz - 1)
                    && empty[idx(x, y, z)]
                {
                    let r = uf.find(idx(x, y, z) as u32) as usize;
                    escapes[r] = true;
                }
            }
        }
    }
    // Collect enclosed components (deterministic order by root index).
    let mut groups: std::collections::BTreeMap<u32, Vec<usize>> = std::collections::BTreeMap::new();
    for (i, &is_empty) in empty.iter().enumerate() {
        if is_empty {
            let r = uf.find(i as u32);
            if !escapes[r as usize] {
                groups.entry(r).or_default().push(i);
            }
        }
    }
    groups
        .into_values()
        .map(|voxels| {
            let depth = voxels
                .iter()
                .map(|&i| field.values[i] - level)
                .fold(f64::INFINITY, f64::min);
            (depth, voxels)
        })
        .collect()
}

/// Evaluate the persistence penalty of `field` against `spec`.
///
/// Channels:
/// - `excess-component`: keep the `spec.components` most persistent H₀
///   bars; every further bar with persistence > τ contributes its
///   persistence, attributed to the voxels of that component (carve).
/// - `component-deficit`: `(target − persistent count) · τ` — graded
///   only in count (no localization; documented).
/// - `enclosed-void`: every void deeper than τ beyond the allowance
///   contributes `depth · |void|^{1/3}` (depth-weighted, scale-aware),
///   attributed to its voxels (fill).
/// - `tunnel-mismatch`: `|b₁ − target| · τ` (counted, not localized —
///   the H₁ no-claim of this tier).
#[must_use]
#[allow(clippy::cast_precision_loss)]
pub fn evaluate(field: &VoxelField, spec: &TopoSpec) -> TopoPenalty {
    let bars0 = persistence0(field);
    let b = betti(field, spec.level);
    let mut attributions = Vec::new();
    let mut total = 0.0f64;
    // H0 (bead 84ib): the component count of THE SOLID is the number of
    // sublevel components ALIVE AT the spec level — bars with
    // birth ≤ level < death — which coincides with b₀(level). Counting
    // long-lived bars ANYWHERE in the filtration (the old form) turned
    // internal density basins separated by still-solid saddles into
    // phantom components: a connected solid with two basins carried a
    // false excess-component penalty while the SAME report said
    // betti = (1, ·, ·). b₀ is the count authority; alive bars supply
    // the voxel attribution, ranked so the LEAST persistent (most
    // weakly attached) components are carved first.
    let mut alive: Vec<&Bar> = bars0
        .iter()
        .filter(|bar| bar.birth <= spec.level && bar.death > spec.level)
        .collect();
    alive.sort_by(|x, y| {
        y.persistence()
            .total_cmp(&x.persistence())
            .then(x.birth.total_cmp(&y.birth))
    });
    // The bars and b₀ can DISAGREE legitimately: persistence0's union-
    // find connectivity is not guaranteed to match betti's 26-connected
    // solid (measured on tp_001: 2 alive bars, b₀ = 1). b₀ is the count
    // AUTHORITY (the report can never contradict itself again); bars
    // supply voxel localization only when the two agree.
    let solid_components = b.0 as usize;
    let localizable = alive.len() == solid_components;
    if solid_components > spec.components as usize {
        if localizable {
            for bar in alive.iter().skip(spec.components as usize) {
                let amount = bar.persistence().min(1e12); // essential guard
                total += amount;
                attributions.push(Attribution {
                    channel: "excess-component",
                    amount,
                    voxels: component_voxels(field, spec.level, bar),
                    direction: -1.0,
                });
            }
        } else {
            // Counted, not localized (the honest degradation when the
            // filtration proxy disagrees with the solid's connectivity).
            let amount = f64::from(b.0 - spec.components) * spec.tau;
            total += amount;
            attributions.push(Attribution {
                channel: "excess-component",
                amount,
                voxels: Vec::new(),
                direction: -1.0,
            });
        }
    } else if (solid_components as u32) < spec.components {
        let amount = f64::from(spec.components - solid_components as u32) * spec.tau;
        total += amount;
        attributions.push(Attribution {
            channel: "component-deficit",
            amount,
            voxels: Vec::new(),
            direction: 1.0,
        });
    }
    // H2 via the duality route: enclosed voids of the empty phase.
    let voids = enclosed_voids(field, spec.level);
    let mut deep: Vec<&(f64, Vec<usize>)> = voids.iter().filter(|(d, _)| *d > spec.tau).collect();
    deep.sort_by(|a, b| b.0.total_cmp(&a.0).then(a.1[0].cmp(&b.1[0])));
    let mut void_depths = Vec::new();
    for (k, (depth, voxels)) in deep.iter().enumerate() {
        void_depths.push(*depth);
        if k < spec.enclosed_voids as usize {
            continue; // allowed
        }
        let amount = depth * (voxels.len() as f64).cbrt();
        total += amount;
        attributions.push(Attribution {
            channel: "enclosed-void",
            amount,
            voxels: voxels.clone(),
            direction: 1.0,
        });
    }
    // H1: counted only (the no-claim boundary of this tier).
    if b.1 != spec.tunnels {
        let amount = f64::from(b.1.abs_diff(spec.tunnels)) * spec.tau;
        total += amount;
        attributions.push(Attribution {
            channel: "tunnel-mismatch",
            amount,
            voxels: Vec::new(),
            direction: if b.1 < spec.tunnels { -1.0 } else { 1.0 },
        });
    }
    TopoPenalty {
        total,
        attributions,
        bars0,
        betti: b,
        void_depths,
    }
}

/// The voxels of the sublevel component a bar was born from: flood
/// from the minimum-value voxel among those below the bar's death,
/// restricted to the occupied phase at `level` — the carve-attribution
/// set for an excess component.
fn component_voxels(field: &VoxelField, level: f64, bar: &Bar) -> Vec<usize> {
    let [nx, ny, nz] = field.dims;
    let n = field.values.len();
    // Seed: the voxel whose value equals the birth (deterministic:
    // first index on ties).
    let mut seed = None;
    for (i, &v) in field.values.iter().enumerate() {
        if (v - bar.birth).abs() < 1e-12 {
            seed = Some(i);
            break;
        }
    }
    let Some(seed) = seed else {
        return Vec::new();
    };
    let cap = bar.death.min(level);
    let mut seen = vec![false; n];
    let mut stack = vec![seed];
    let mut out = Vec::new();
    seen[seed] = true;
    let idx = |x: u32, y: u32, z: u32| -> usize { ((z * ny + y) * nx + x) as usize };
    while let Some(i) = stack.pop() {
        if field.values[i] >= cap {
            continue;
        }
        out.push(i);
        let iu = i as u32;
        let x = iu % nx;
        let y = (iu / nx) % ny;
        let z = iu / (nx * ny);
        let mut push = |j: usize| {
            if !seen[j] {
                seen[j] = true;
                stack.push(j);
            }
        };
        if x > 0 {
            push(idx(x - 1, y, z));
        }
        if x + 1 < nx {
            push(idx(x + 1, y, z));
        }
        if y > 0 {
            push(idx(x, y - 1, z));
        }
        if y + 1 < ny {
            push(idx(x, y + 1, z));
        }
        if z > 0 {
            push(idx(x, y, z - 1));
        }
        if z + 1 < nz {
            push(idx(x, y, z + 1));
        }
    }
    out.sort_unstable();
    out
}

/// The HEURISTIC FALLBACK the [M] gate compares against: binary
/// connected-component labeling — penalty = count mismatches, no
/// magnitude, no localization (zero gradient until a violation
/// disappears entirely).
#[must_use]
pub fn heuristic_cc_penalty(field: &VoxelField, spec: &TopoSpec) -> f64 {
    let b = betti(field, spec.level);
    let voids = enclosed_voids(field, spec.level)
        .iter()
        .filter(|(d, _)| *d > spec.tau)
        .count() as u32;
    f64::from(b.0.abs_diff(spec.components))
        + f64::from(b.1.abs_diff(spec.tunnels))
        + f64::from(voids.saturating_sub(spec.enclosed_voids))
}

/// One attribution-guided descent step: move the field toward
/// compliance where the attributions say, by `step` (fill raises
/// occupancy = LOWERS the value; carve raises it). Returns the number
/// of voxels touched.
pub fn apply_attribution_step(field: &mut VoxelField, report: &TopoPenalty, step: f64) -> usize {
    let mut touched = 0usize;
    for att in &report.attributions {
        for &i in &att.voxels {
            // direction +1 = fill (solid = value < level): subtract.
            field.values[i] -= att.direction * step;
            touched += 1;
        }
    }
    touched
}
