//! Declared triangle near rule for the closed-surface Burton-Miller rows
//! (bead frankensim-irp-v5-instrument-realism-pexfp.111).
//!
//! The default [`PanelRule::Centroid`] evaluates every off-diagonal panel as
//! one point kernel times its area. Near the collocation point that weight
//! carries an O(h) error in S, D and the hypersingular N, so the combined rows
//! converge at first order while plain triangle CBIE converges at second.
//!
//! [`PanelRule::TriangleNear`] keeps the centroid rows and replaces every
//! `triangle_is_near` pair (the plain-CBIE near zone, four circumradii) by
//! flat-triangle integrals of constant density:
//! - S and D: the plain CBIE arm's exact static layers plus 8x8 Gauss
//!   remainders;
//! - N: the exact finite part
//!   `N[chi_T](x) = (n_x.n)[k^2 S_T(x) + SUM_e m_e.I_e] - SUM_e (n_x.m_e)(n.I_e)`,
//!   `I_e = INT_e grad_y G dl`, with `m_e` the outward in-plane edge normal.
//!   It follows from `(n.grad)grad G = -k^2 G n - n Lap_t G + grad_t dG/dn`
//!   and the in-plane divergence theorem. The static edge integrals are
//!   closed form; the bounded dynamic remainder uses 8-point Gauss per edge.
//!
//! Self terms are the flat triangle's at its centroid (exact static plus a
//! Duffy fan). The static closures `(D_0+I/2)1 = 0` and `N_0[1] = 0` are
//! recomputed under the same rule, as in the centroid rows. The adjoint double
//! layer D' keeps its centroid kernel and point closure. Exact flat-facet D'
//! is itself first order: a faceted surface has no adjoint solid-angle
//! identity, and `D'_0[1] + 1/2 = O(h)` on the icosphere. Far pairs keep
//! one-point weights, so the FMM route applies the rule as a sparse correction
//! without touching its expansions.

use super::iterative::{NearOperator, NearRows, RADON7};
use super::{
    DENSE_WORKERS, Formulation, GL8_W, GL8_X, HelmholtzError, Medium, RadiationOperator,
    RadiationSolution, alpha_for, assemble_dense, condition_lower_bound, cross, dot, fill_rows,
    kernels, norm, regularized_coplanar_hypersingular, self_terms, sub, triangle_is_near,
    triangle_self_terms, triangle_static_influence, triangle_weak_influence_with_static,
};
use crate::panel3d::SpherePanels;
use fs_la::eigen_complex::lu_complex;
use fs_math::c64::C64;
use std::collections::HashMap;

/// Declared panel quadrature of the closed-surface Burton-Miller rows.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum PanelRule {
    /// One point kernel per off-diagonal panel and equivalent-disc self terms.
    #[default]
    Centroid,
    /// Exact flat-triangle S, D and finite-part N on every near pair (see the
    /// module docs). Requires retained triangles and refuses a surface with
    /// more than `max_pairs` near (target, source) pairs.
    TriangleNear {
        /// Declared ceiling on corrected pairs.
        max_pairs: usize,
    },
    /// [`Self::TriangleNear`] plus the Radon seven-point rule for S, D and N on
    /// every far pair, with both static closures under the same rule. One-point
    /// far weights leave an O(h) remnant that `alpha = i/k` amplifies at low
    /// ka (a dipole's radiation resistance stalls); this rule removes it.
    /// Dense route only: the iterative owner refuses it, since its far field
    /// is one point per panel.
    Triangle {
        /// Declared ceiling on corrected near pairs.
        max_pairs: usize,
    },
}

/// The rule's sparse correction to the centroid rows. Row `i` lists its near
/// sources `j != i` in ascending order with area-integrated
/// `accurate - centroid` values of S, D and N. Diagonal deltas carry the
/// triangle self terms and the static rows re-closed under the rule.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct NearRule {
    starts: Vec<usize>,
    cols: Vec<usize>,
    single: Vec<C64>,
    double: Vec<C64>,
    hyper: Vec<C64>,
    static_double: Vec<f64>,
    static_hyper: Vec<f64>,
    self_single: Vec<C64>,
    self_hyper: Vec<C64>,
}

struct RowValues {
    cols: Vec<usize>,
    single: Vec<C64>,
    double: Vec<C64>,
    hyper: Vec<C64>,
    static_double: f64,
    static_hyper: f64,
    self_single: C64,
    self_hyper: C64,
}

impl NearRule {
    /// Build the correction at real wavenumber `k` on `workers` threads
    /// (bit-identical for every count). `cancelled` is polled between rows.
    pub(crate) fn build(
        surface: &SpherePanels,
        k: f64,
        max_pairs: usize,
        workers: usize,
        cancelled: &(dyn Fn() -> bool + Sync),
    ) -> Result<Self, HelmholtzError> {
        let triangles = surface.triangles().ok_or(HelmholtzError::BadParameter {
            what: "the triangle near rule requires retained panel triangles",
        })?;
        if !(k > 0.0 && k.is_finite()) {
            return Err(HelmholtzError::BadParameter {
                what: "the triangle near rule requires a positive finite wavenumber",
            });
        }
        let sources = near_sources(surface.centroids(), triangles);
        let required = sources.iter().map(Vec::len).sum::<usize>();
        if required > max_pairs {
            return Err(HelmholtzError::NearRuleBudget {
                required_pairs: required,
                maximum_pairs: max_pairs,
            });
        }
        let n = sources.len();
        let workers = workers.clamp(1, n.max(1));
        let row = |i: usize| -> Result<RowValues, HelmholtzError> {
            if cancelled() {
                return Err(HelmholtzError::Cancelled);
            }
            row_values(surface, triangles, k, i, &sources[i])
        };
        let rows: Vec<RowValues> = if workers == 1 {
            (0..n).map(row).collect::<Result<_, _>>()?
        } else {
            let chunks: Vec<Result<Vec<RowValues>, HelmholtzError>> = std::thread::scope(|scope| {
                let handles: Vec<_> = (0..workers)
                    .map(|p| {
                        let row = &row;
                        scope.spawn(move || {
                            (p * n / workers..(p + 1) * n / workers).map(row).collect()
                        })
                    })
                    .collect();
                handles
                    .into_iter()
                    .map(|handle| handle.join().expect("near-rule worker panicked"))
                    .collect()
            });
            let mut rows = Vec::with_capacity(n);
            for chunk in chunks {
                rows.extend(chunk?);
            }
            rows
        };
        let mut out = Self {
            starts: Vec::with_capacity(n + 1),
            cols: Vec::with_capacity(required),
            single: Vec::with_capacity(required),
            double: Vec::with_capacity(required),
            hyper: Vec::with_capacity(required),
            static_double: Vec::with_capacity(n),
            static_hyper: Vec::with_capacity(n),
            self_single: Vec::with_capacity(n),
            self_hyper: Vec::with_capacity(n),
        };
        out.starts.push(0);
        for row in rows {
            out.cols.extend(row.cols);
            out.single.extend(row.single);
            out.double.extend(row.double);
            out.hyper.extend(row.hyper);
            out.starts.push(out.cols.len());
            out.static_double.push(row.static_double);
            out.static_hyper.push(row.static_hyper);
            out.self_single.push(row.self_single);
            out.self_hyper.push(row.self_hyper);
        }
        Ok(out)
    }

    /// Corrected (target, source) pairs.
    pub(crate) fn pairs(&self) -> usize {
        self.cols.len()
    }

    /// Diagonal deltas `(dA_ii, dB_ii)` of the combined rows for `alpha`.
    /// D' keeps its point closure, so B changes only by the triangle S_ii.
    pub(crate) fn diagonal(&self, i: usize, alpha: C64) -> (C64, C64) {
        (
            C64::from_re(-self.static_double[i])
                + alpha * (self.self_hyper[i] - C64::from_re(self.static_hyper[i])),
            self.self_single[i],
        )
    }

    /// Row `i`'s near sources in ascending order, each with its
    /// area-integrated `(dS, dD, dN)`. Off-diagonal BM entries combine as
    /// `A_ij += dD + alpha dN` and `B_ij += dS`; D' is unchanged.
    pub(crate) fn for_row(&self, i: usize) -> impl Iterator<Item = (usize, C64, C64, C64)> + '_ {
        (self.starts[i]..self.starts[i + 1]).map(|at| {
            (
                self.cols[at],
                self.single[at],
                self.double[at],
                self.hyper[at],
            )
        })
    }

    /// Correct a dense row-major centroid pair `(A, B)` in place.
    fn correct_dense(&self, alpha: C64, a: &mut [C64], b: &mut [C64]) {
        let n = self.self_single.len();
        for i in 0..n {
            let (da, db) = self.diagonal(i, alpha);
            a[i * n + i] = a[i * n + i] + da;
            b[i * n + i] = b[i * n + i] + db;
            for (j, single, double, hyper) in self.for_row(i) {
                a[i * n + j] = a[i * n + j] + double + alpha * hyper;
                b[i * n + j] = b[i * n + j] + single;
            }
        }
    }
}

/// Far rows for [`PanelRule::Triangle`]. Every pair outside row `i`'s
/// near zone and off the diagonal replaces its centroid S, D and N by the
/// source panel's Radon seven-point sums of the same point kernels, and
/// both static closures follow the same rule. D' keeps its centroid rows.
/// Rows run on the declared dense workers and read shared data only, so every
/// bit is independent of the worker count.
fn correct_dense_far(
    surface: &SpherePanels,
    k: f64,
    alpha: C64,
    near: &NearRule,
    a: &mut [C64],
    b: &mut [C64],
) {
    let triangles = surface
        .triangles()
        .expect("NearRule::build admitted retained triangles");
    let (centroids, normals, areas) = (surface.centroids(), surface.normals(), surface.areas());
    let n = centroids.len();
    let fill = |i: usize, arow: &mut [C64], brow: &mut [C64]| {
        let (x, nx) = (centroids[i], normals[i]);
        let row = &near.cols[near.starts[i]..near.starts[i + 1]];
        let (mut static_double, mut static_hyper) = (0.0, 0.0);
        let mut next = 0;
        for j in 0..n {
            if next < row.len() && row[next] == j {
                next += 1;
                continue;
            }
            if j == i {
                continue;
            }
            let (ny, area) = (normals[j], areas[j]);
            let ([single, double, hyper], [d0, n0]) =
                triangle_far_layers(k, x, nx, (triangles[j], centroids[j], area), ny);
            let (g, dgdny, _, d2g) = kernels(k, x, nx, centroids[j], ny);
            let (_, point_d0, _, point_n0) = kernels(0.0, x, nx, centroids[j], ny);
            arow[j] = arow[j] + (double - dgdny.scale(area)) + alpha * (hyper - d2g.scale(area));
            brow[j] = brow[j] + (single - g.scale(area));
            static_double += d0 - point_d0.re * area;
            static_hyper += n0 - point_n0.re * area;
        }
        arow[i] = arow[i] + C64::from_re(-static_double) - alpha.scale(static_hyper);
    };
    fill_rows(a, b, n, &fill);
}

/// Radon seven-point sums of the point kernels over a source panel at a
/// well-separated target: `([S, D, N], [D_0, N_0])`. The nodes and weights
/// are the iterative owner's far-field cloud's: the panel centroid first,
/// then barycentric nodes, each weight times the panel area.
pub(crate) fn triangle_far_layers(
    k: f64,
    x: [f64; 3],
    nx: [f64; 3],
    panel: ([[f64; 3]; 3], [f64; 3], f64),
    ny: [f64; 3],
) -> ([C64; 3], [f64; 2]) {
    let (triangle, centroid, area) = panel;
    let mut layers = [C64::ZERO; 3];
    let mut statics = [0.0; 2];
    for (m, &(barycentric, weight)) in RADON7.iter().enumerate() {
        let y = if m == 0 {
            centroid
        } else {
            core::array::from_fn(|axis| {
                barycentric[0] * triangle[0][axis]
                    + barycentric[1] * triangle[1][axis]
                    + barycentric[2] * triangle[2][axis]
            })
        };
        let weight = weight * area;
        let (g, dgdny, _, d2g) = kernels(k, x, nx, y, ny);
        let (_, d0, _, n0) = kernels(0.0, x, nx, y, ny);
        layers[0] = layers[0] + g.scale(weight);
        layers[1] = layers[1] + dgdny.scale(weight);
        layers[2] = layers[2] + d2g.scale(weight);
        statics[0] += d0.re * weight;
        statics[1] += n0.re * weight;
    }
    (layers, statics)
}

/// Sources `j != i` inside the plain-CBIE near zone of each centroid target,
/// ascending. The predicate is exactly [`triangle_is_near`]; a uniform grid
/// whose cell is the largest reach only narrows the candidates.
fn near_sources(centroids: &[[f64; 3]], triangles: &[[[f64; 3]; 3]]) -> Vec<Vec<usize>> {
    let cell = centroids
        .iter()
        .zip(triangles)
        .map(|(&c, triangle)| {
            4.0 * triangle
                .iter()
                .map(|&vertex| norm(sub(vertex, c)))
                .fold(0.0_f64, f64::max)
        })
        .fold(0.0_f64, f64::max);
    // `|x - c_j| <= 4 r_j <= cell` puts c_j in one of the 27 cells around x.
    let key = |x: [f64; 3]| x.map(|v| (v / cell).floor() as i64);
    let mut grid: HashMap<[i64; 3], Vec<usize>> = HashMap::new();
    for (j, &c) in centroids.iter().enumerate() {
        grid.entry(key(c)).or_default().push(j);
    }
    centroids
        .iter()
        .enumerate()
        .map(|(i, &x)| {
            let [a, b, c] = key(x);
            let mut row = Vec::new();
            for da in -1..=1 {
                for db in -1..=1 {
                    for dc in -1..=1 {
                        for &j in grid.get(&[a + da, b + db, c + dc]).into_iter().flatten() {
                            if j != i && triangle_is_near(x, centroids[j], triangles[j]) {
                                row.push(j);
                            }
                        }
                    }
                }
            }
            row.sort_unstable();
            row
        })
        .collect()
}

fn row_values(
    surface: &SpherePanels,
    triangles: &[[[f64; 3]; 3]],
    k: f64,
    i: usize,
    sources: &[usize],
) -> Result<RowValues, HelmholtzError> {
    let centroids = surface.centroids();
    let normals = surface.normals();
    let areas = surface.areas();
    let (x, nx) = (centroids[i], normals[i]);
    let mut out = RowValues {
        cols: sources.to_vec(),
        single: Vec::with_capacity(sources.len()),
        double: Vec::with_capacity(sources.len()),
        hyper: Vec::with_capacity(sources.len()),
        static_double: 0.0,
        static_hyper: 0.0,
        self_single: C64::ZERO,
        self_hyper: C64::ZERO,
    };
    for &j in sources {
        let (triangle, ny, area) = (triangles[j], normals[j], areas[j]);
        let statics = triangle_static_influence(x, triangle, ny);
        let (single, double) =
            triangle_weak_influence_with_static(k, x, triangle, ny, (&GL8_X, &GL8_W), statics);
        let (hyper, static_hyper) =
            triangle_hypersingular_influence(k, x, nx, triangle, ny, single);
        // The centroid rows' own point values, removed exactly.
        let (g, dgdny, _, d2g) = kernels(k, x, nx, centroids[j], ny);
        let (_, d0, _, n0) = kernels(0.0, x, nx, centroids[j], ny);
        let values = [
            single - g.scale(area),
            double - dgdny.scale(area),
            hyper - d2g.scale(area),
        ];
        let static_values = [statics.1 - d0.re * area, static_hyper - n0.re * area];
        if values
            .iter()
            .any(|v| !(v.re.is_finite() && v.im.is_finite()))
            || static_values.iter().any(|v| !v.is_finite())
        {
            return Err(HelmholtzError::BadParameter {
                what: "triangle near rule: non-finite pair integral (target on a source edge?)",
            });
        }
        out.single.push(values[0]);
        out.double.push(values[1]);
        out.hyper.push(values[2]);
        out.static_double += static_values[0];
        out.static_hyper += static_values[1];
    }
    let (single, hyper_regular) = triangle_self_terms(k, x, triangles[i]);
    let (disc_single, disc_hyper) = self_terms(k, areas[i]);
    out.self_single = single - disc_single;
    out.self_hyper = hyper_regular - disc_hyper;
    Ok(out)
}

/// Finite-part hypersingular layer `(N_k, N_0)` of a flat source triangle with
/// unit density, at a target `x` off the closed triangle with normal `nx`.
/// `single` is the same triangle's single layer `S_T(x)`.
fn triangle_hypersingular_influence(
    k: f64,
    x: [f64; 3],
    nx: [f64; 3],
    triangle: [[f64; 3]; 3],
    ny: [f64; 3],
    single: C64,
) -> (C64, f64) {
    let four_pi = 4.0 * core::f64::consts::PI;
    let normal_dot = dot(nx, ny);
    let mut static_part = 0.0;
    let mut tangential = C64::ZERO;
    let mut normal = C64::ZERO;
    for edge in 0..3 {
        let start = triangle[edge];
        let end = triangle[(edge + 1) % 3];
        let direction = sub(end, start);
        let length = norm(direction);
        let tangent = direction.map(|v| v / length);
        let outward = cross(tangent, ny);
        let w = sub(start, x);
        let u0 = dot(w, tangent);
        let u1 = u0 + length;
        let rho = [
            w[0] - u0 * tangent[0],
            w[1] - u0 * tangent[1],
            w[2] - u0 * tangent[2],
        ];
        let (r0, r1) = (norm(w), norm(sub(end, x)));
        // INT_e du/r^3 = [u/(rho^2 r)]. When the foot of x lies off the edge,
        // the equivalent form avoids subtracting nearly equal endpoint terms.
        let inverse_cube = if u0 * u1 > 0.0 {
            length * (u0 + u1) / (r0 * r1 * (u1 * r0 + u0 * r1))
        } else {
            (u1 / r1 - u0 / r0) / dot(rho, rho)
        };
        let axial = 1.0 / r0 - 1.0 / r1;
        // INT_e grad_y G_0 dl = -INT_e (y - x)/(4 pi r^3) dl.
        let gradient = [
            -(rho[0] * inverse_cube + tangent[0] * axial) / four_pi,
            -(rho[1] * inverse_cube + tangent[1] * axial) / four_pi,
            -(rho[2] * inverse_cube + tangent[2] * axial) / four_pi,
        ];
        let outward_dot = dot(nx, outward);
        static_part += normal_dot * dot(outward, gradient) - outward_dot * dot(ny, gradient);
        // INT_e (grad_y G - grad_y G_0) dl: the factor (G' - G_0')/r equals
        // minus the regularized coplanar hypersingular kernel; it is bounded,
        // but its direction turns on the scale of the target distance, so each
        // half edge gets its own 8-point rule.
        let mut remainder = [C64::ZERO; 3];
        let halves = GL8_X
            .iter()
            .zip(&GL8_W)
            .flat_map(|(&s, &w)| [(0.5 * s, 0.5 * w), (0.5 + 0.5 * s, 0.5 * w)]);
        for (s, weight) in halves {
            let y = [
                start[0] + s * direction[0],
                start[1] + s * direction[1],
                start[2] + s * direction[2],
            ];
            let d = sub(y, x);
            let factor = regularized_coplanar_hypersingular(k, norm(d)).scale(-weight * length);
            for (value, &component) in remainder.iter_mut().zip(&d) {
                *value = *value + factor.scale(component);
            }
        }
        let along = |axis: [f64; 3]| {
            remainder[0].scale(axis[0]) + remainder[1].scale(axis[1]) + remainder[2].scale(axis[2])
        };
        tangential = tangential + along(outward);
        normal = normal + along(ny).scale(outward_dot);
    }
    let dynamic = (single.scale(k * k) + tangential).scale(normal_dot) - normal;
    (dynamic + C64::from_re(static_part), static_part)
}

/// [`super::solve_radiation_batch`] under a declared [`PanelRule`]. The
/// centroid rule is that function, bit for bit. The triangle rule corrects
/// the Burton-Miller rows before the same LU and diagnostics. Plain CBIE ignores
/// any rule: it keeps its own triangle quadrature, so routed callers can pass
/// one rule to rows of either formulation.
///
/// # Errors
/// The admissions of [`super::solve_radiation_batch`], plus the rule's
/// refusals: missing triangles or [`HelmholtzError::NearRuleBudget`].
pub fn solve_radiation_batch_with_rule(
    surface: &SpherePanels,
    k: f64,
    medium: Medium,
    velocity_fields: &[&[C64]],
    formulation: Formulation,
    rule: PanelRule,
) -> Result<Vec<RadiationSolution>, HelmholtzError> {
    if formulation == Formulation::PlainCbie {
        return super::solve_radiation_batch(surface, k, medium, velocity_fields, formulation);
    }
    let (max_pairs, gauss_far) = match rule {
        PanelRule::Centroid => {
            return super::solve_radiation_batch(surface, k, medium, velocity_fields, formulation);
        }
        PanelRule::TriangleNear { max_pairs } => (max_pairs, false),
        PanelRule::Triangle { max_pairs } => (max_pairs, true),
    };
    let operator = RadiationOperator::prepare_ruled(
        surface,
        k,
        medium,
        velocity_fields,
        formulation,
        |alpha| {
            let workers = DENSE_WORKERS.with(core::cell::Cell::get);
            let rule = NearRule::build(surface, k, max_pairs, workers, &|| false)?;
            let (mut a, mut b) = assemble_dense(surface, k, alpha);
            rule.correct_dense(alpha, &mut a, &mut b);
            if gauss_far {
                correct_dense_far(surface, k, alpha, &rule, &mut a, &mut b);
            }
            Ok((a, b))
        },
    )?;
    Ok(velocity_fields
        .iter()
        .map(|&velocity| operator.solve(velocity))
        .collect())
}

/// One-field form of [`solve_radiation_batch_with_rule`].
///
/// # Errors
/// As for [`solve_radiation_batch_with_rule`].
pub fn solve_radiation_with_rule(
    surface: &SpherePanels,
    k: f64,
    medium: Medium,
    velocity: &[C64],
    formulation: Formulation,
    rule: PanelRule,
) -> Result<RadiationSolution, HelmholtzError> {
    let mut solutions =
        solve_radiation_batch_with_rule(surface, k, medium, &[velocity], formulation, rule)?;
    Ok(solutions.remove(0))
}

impl<'a> RadiationOperator<'a> {
    /// Admit exactly as [`RadiationOperator::prepare`], then factorize the
    /// pair returned by `assemble(alpha)`.
    fn prepare_ruled(
        surface: &'a SpherePanels,
        k: f64,
        medium: Medium,
        velocity_fields: &[&[C64]],
        formulation: Formulation,
        assemble: impl FnOnce(C64) -> Result<(Vec<C64>, Vec<C64>), HelmholtzError>,
    ) -> Result<Self, HelmholtzError> {
        let ppw = super::admit_operator(surface, k, medium, velocity_fields)?;
        let n = surface.centroids().len();
        let (matrix, bmat) = assemble(alpha_for(formulation, k))?;
        let lu = lu_complex(&matrix, n).map_err(|_| HelmholtzError::Singular)?;
        let condition_lower_bound = condition_lower_bound(&matrix, n, &lu);
        Ok(Self {
            surface,
            k,
            medium,
            panels_per_wavelength: ppw,
            bmat,
            lu,
            condition_lower_bound,
        })
    }
}

/// The iterative owner's near rows are this rule's corrections: pair-major
/// `(dS, dD, dN)` on top of its centroid FMM operator, D' unchanged.
impl NearRows for NearRule {
    fn operators(&self) -> &[NearOperator] {
        &[
            NearOperator::Single,
            NearOperator::Double,
            NearOperator::Hypersingular,
        ]
    }
    fn row_len(&self, i: usize) -> usize {
        self.starts[i + 1] - self.starts[i]
    }
    fn fill_row(&self, i: usize, cols: &mut [u32], values: &mut [C64]) {
        for (p, (j, single, double, hyper)) in self.for_row(i).enumerate() {
            // Surfaces hold at most MAX_SURFACE_PANELS, far below u32::MAX.
            cols[p] = j as u32;
            values[3 * p] = single;
            values[3 * p + 1] = double;
            values[3 * p + 2] = hyper;
        }
    }
    fn diagonal(&self, i: usize, alpha: C64) -> (C64, C64) {
        NearRule::diagonal(self, i, alpha)
    }
    fn stamp(&self) -> u64 {
        // The rule's tag, "tri-near".
        0x7472_692d_6e65_6172
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Brute force: 4^levels congruent subtriangles, each with a 4x4 collapsed
    /// Gauss rule, summing the production point kernels (G, dG/dny, d2G).
    fn subdivided(
        k: f64,
        x: [f64; 3],
        nx: [f64; 3],
        triangle: [[f64; 3]; 3],
        ny: [f64; 3],
        levels: u32,
    ) -> [C64; 3] {
        let mut pieces = vec![triangle];
        for _ in 0..levels {
            let mut next = Vec::with_capacity(4 * pieces.len());
            for [a, b, c] in pieces {
                let mid = |p: [f64; 3], q: [f64; 3]| {
                    [
                        0.5 * (p[0] + q[0]),
                        0.5 * (p[1] + q[1]),
                        0.5 * (p[2] + q[2]),
                    ]
                };
                let (ab, bc, ca) = (mid(a, b), mid(b, c), mid(c, a));
                next.extend([[a, ab, ca], [ab, b, bc], [ca, bc, c], [ab, bc, ca]]);
            }
            pieces = next;
        }
        let mut out = [C64::ZERO; 3];
        for [a, b, c] in pieces {
            let (ab, ac) = (sub(b, a), sub(c, a));
            let double_area = norm(cross(ab, ac));
            for (&u, &wu) in super::super::GL4_X.iter().zip(&super::super::GL4_W) {
                for (&v, &wv) in super::super::GL4_X.iter().zip(&super::super::GL4_W) {
                    let y = super::super::affine3(a, ab, u, ac, (1.0 - u) * v);
                    let weight = wu * wv * (1.0 - u) * double_area;
                    let (g, dgdny, _, d2g) = kernels(k, x, nx, y, ny);
                    out[0] = out[0] + g.scale(weight);
                    out[1] = out[1] + dgdny.scale(weight);
                    out[2] = out[2] + d2g.scale(weight);
                }
            }
        }
        out
    }

    fn unit(v: [f64; 3]) -> [f64; 3] {
        let l = norm(v);
        v.map(|c| c / l)
    }

    #[test]
    fn finite_part_triangle_layers_match_subdivided_point_quadrature() {
        let triangle = [[0.0, 0.0, 0.0], [1.0, 0.1, 0.0], [0.3, 0.9, 0.05]];
        let ny = unit(cross(
            sub(triangle[1], triangle[0]),
            sub(triangle[2], triangle[0]),
        ));
        let c = [
            (triangle[0][0] + triangle[1][0] + triangle[2][0]) / 3.0,
            (triangle[0][1] + triangle[1][1] + triangle[2][1]) / 3.0,
            (triangle[0][2] + triangle[1][2] + triangle[2][2]) / 3.0,
        ];
        // Near-coplanar neighbours across an edge and a vertex, a target above
        // the interior, and oblique well-separated targets with tilted normals.
        let cases = [
            ([1.05, 0.75, 0.02], unit([0.02, -0.01, 1.0])),
            ([0.5, -0.35, -0.01], unit([0.0, 0.03, 1.0])),
            ([1.45, 0.05, 0.03], ny),
            (
                [c[0] + 0.05, c[1] - 0.04, c[2] + 0.3],
                unit([0.3, 0.2, 0.9]),
            ),
            ([c[0] - 0.8, c[1] + 0.6, c[2] - 0.5], unit([-0.5, 0.4, 0.7])),
            ([2.1, 1.7, 0.9], unit([0.1, -0.9, 0.2])),
        ];
        for (x, nx) in cases {
            for k in [0.0, 1.0, 5.0] {
                let statics = triangle_static_influence(x, triangle, ny);
                let (s, d) = triangle_weak_influence_with_static(
                    k,
                    x,
                    triangle,
                    ny,
                    (&GL8_X, &GL8_W),
                    statics,
                );
                let (h, h0) = triangle_hypersingular_influence(k, x, nx, triangle, ny, s);
                let reference = subdivided(k, x, nx, triangle, ny, 6);
                let static_reference = subdivided(0.0, x, nx, triangle, ny, 6);
                for (got, want, what) in [
                    (s, reference[0], "S"),
                    (d, reference[1], "D"),
                    (h, reference[2], "N"),
                ] {
                    let error = (got - want).abs() / want.abs();
                    assert!(
                        error < 1e-6,
                        "{what} k={k} x={x:?}: {got:?} vs {want:?} ({error:e})"
                    );
                }
                let error = (h0 - static_reference[2].re).abs() / static_reference[2].abs();
                assert!(
                    error < 1e-9,
                    "N0 x={x:?}: {h0} vs {:?} ({error:e})",
                    static_reference[2]
                );
            }
        }
    }

    /// Static finite part of a flat triangle's N_0 at its own centroid,
    /// `-(1/4 pi) SUM_e d_e INT_e dl/r^3`, independent of the off-panel path.
    fn static_self_finite_part(triangle: [[f64; 3]; 3]) -> f64 {
        let x = [0, 1, 2].map(|c| (triangle[0][c] + triangle[1][c] + triangle[2][c]) / 3.0);
        let ny = unit(cross(
            sub(triangle[1], triangle[0]),
            sub(triangle[2], triangle[0]),
        ));
        let mut sum = 0.0;
        for edge in 0..3 {
            let (a, b) = (triangle[edge], triangle[(edge + 1) % 3]);
            let tangent = unit(sub(b, a));
            let distance = dot(sub(a, x), cross(tangent, ny));
            let (u0, u1) = (dot(sub(a, x), tangent), dot(sub(b, x), tangent));
            let r = |u: f64| (distance * distance + u * u).sqrt();
            sum += distance * (u1 / r(u1) - u0 / r(u0)) / (distance * distance);
        }
        -sum / (4.0 * core::f64::consts::PI)
    }

    #[test]
    fn closed_polyhedron_static_rows_close_under_the_exact_rule() {
        // Every pair integrated exactly: the double layer of unit density is
        // -1/2 at every face centroid of a closed polyhedron and its normal
        // derivative vanishes, so each N_0 row sum is minus the independent
        // self finite part.
        let surface = SpherePanels::icosphere(1.0, 2).unwrap();
        let triangles = surface.triangles().unwrap();
        let n = triangles.len();
        for i in (0..n).step_by(7) {
            let (x, nx) = (surface.centroids()[i], surface.normals()[i]);
            let (mut d0, mut n0) = (0.0, 0.0);
            for j in (0..n).filter(|&j| j != i) {
                let ny = surface.normals()[j];
                d0 += triangle_static_influence(x, triangles[j], ny).1;
                n0 += triangle_hypersingular_influence(0.0, x, nx, triangles[j], ny, C64::ZERO).1;
            }
            let own = static_self_finite_part(triangles[i]);
            assert!((d0 + 0.5).abs() < 1e-9, "row {i}: D0 sum {d0}");
            assert!(
                (n0 + own).abs() < 1e-9 * own.abs(),
                "row {i}: N0 sum {n0} vs self {own}"
            );
        }
    }

    #[test]
    fn near_rule_refusals_budget_boundary_and_worker_bits() {
        let surface = SpherePanels::icosphere(1.0, 2).unwrap();
        let k = 1.0;
        let full = NearRule::build(&surface, k, usize::MAX, 1, &|| false).unwrap();
        let pairs = full.pairs();
        assert!(pairs > 0);
        assert_eq!(
            NearRule::build(&surface, k, pairs, 1, &|| false).unwrap(),
            full
        );
        assert_eq!(
            NearRule::build(&surface, k, pairs - 1, 1, &|| false),
            Err(HelmholtzError::NearRuleBudget {
                required_pairs: pairs,
                maximum_pairs: pairs - 1
            })
        );
        assert_eq!(
            NearRule::build(&surface, k, pairs, 3, &|| false).unwrap(),
            full
        );
        assert_eq!(
            NearRule::build(&surface, k, pairs, 2, &|| true),
            Err(HelmholtzError::Cancelled)
        );
        // The grid only narrows candidates: every near pair, and no other.
        let triangles = surface.triangles().unwrap();
        let sources = near_sources(surface.centroids(), triangles);
        for (i, row) in sources.iter().enumerate() {
            let want: Vec<usize> = (0..triangles.len())
                .filter(|&j| {
                    j != i
                        && triangle_is_near(
                            surface.centroids()[i],
                            surface.centroids()[j],
                            triangles[j],
                        )
                })
                .collect();
            assert_eq!(row, &want);
        }
        let bare = SpherePanels::new(
            surface.centroids().to_vec(),
            surface.normals().to_vec(),
            surface.areas().to_vec(),
        )
        .unwrap();
        assert!(matches!(
            NearRule::build(&bare, k, usize::MAX, 1, &|| false),
            Err(HelmholtzError::BadParameter { .. })
        ));
        let medium = Medium::air();
        let v = vec![C64::ONE; triangles.len()];
        assert!(matches!(
            solve_radiation_with_rule(
                &surface,
                k,
                medium,
                &v,
                Formulation::BurtonMiller,
                PanelRule::TriangleNear {
                    max_pairs: pairs - 1
                }
            ),
            Err(HelmholtzError::NearRuleBudget { .. })
        ));
    }

    #[test]
    fn centroid_rule_and_plain_cbie_keep_the_existing_owner_bits() {
        let surface = SpherePanels::icosphere(1.0, 2).unwrap();
        let medium = Medium::air();
        let v: Vec<C64> = surface
            .centroids()
            .iter()
            .map(|c| C64::new(c[2], 0.5 * c[0]))
            .collect();
        // Plain CBIE ignores every rule: its own triangle quadrature stands.
        let cases = [
            (Formulation::BurtonMiller, PanelRule::Centroid),
            (Formulation::PlainCbie, PanelRule::Centroid),
            (
                Formulation::PlainCbie,
                PanelRule::TriangleNear { max_pairs: 1 << 20 },
            ),
            (
                Formulation::PlainCbie,
                PanelRule::Triangle { max_pairs: 1 << 20 },
            ),
        ];
        for (formulation, rule) in cases {
            let old =
                super::super::solve_radiation(&surface, 1.3, medium, &v, formulation).unwrap();
            let new =
                solve_radiation_with_rule(&surface, 1.3, medium, &v, formulation, rule).unwrap();
            for (a, b) in old.pressure.iter().zip(&new.pressure) {
                assert_eq!(
                    (a.re.to_bits(), a.im.to_bits()),
                    (b.re.to_bits(), b.im.to_bits())
                );
            }
        }
    }

    #[test]
    fn corrected_rows_hold_the_rule_integrals_and_closures() {
        // Composed independently of NearRule: near pairs carry the accurate
        // triangle integrals, far pairs keep centroid bits, D' stays the
        // centroid kernel, and both static rows close under the same rule.
        let surface = SpherePanels::icosphere(1.0, 1).unwrap();
        let triangles = surface.triangles().unwrap();
        let (k, alpha) = (1.7, C64::new(0.0, 1.0 / 1.7));
        let n = surface.areas().len();
        let rule = NearRule::build(&surface, k, usize::MAX, 1, &|| false).unwrap();
        let (a0, b0) = assemble_dense(&surface, k, alpha);
        let (mut a, mut b) = (a0.clone(), b0.clone());
        rule.correct_dense(alpha, &mut a, &mut b);
        let close = |got: C64, want: C64| (got - want).abs() <= 1e-13 * (1.0 + want.abs());
        for i in 0..n {
            let (x, nx) = (surface.centroids()[i], surface.normals()[i]);
            let (mut d0_rule, mut n0_rule, mut d0_point) = (0.0, 0.0, 0.0);
            for j in (0..n).filter(|&j| j != i) {
                let (c, triangle) = (surface.centroids()[j], triangles[j]);
                let (ny, area) = (surface.normals()[j], surface.areas()[j]);
                let (_, point_d0, _, point_n0) = kernels(0.0, x, nx, c, ny);
                d0_point += point_d0.re * area;
                if triangle_is_near(x, c, triangle) {
                    let statics = triangle_static_influence(x, triangle, ny);
                    let (s, d) = triangle_weak_influence_with_static(
                        k,
                        x,
                        triangle,
                        ny,
                        (&GL8_X, &GL8_W),
                        statics,
                    );
                    let (h, h0) = triangle_hypersingular_influence(k, x, nx, triangle, ny, s);
                    let (_, _, dgdnx, _) = kernels(k, x, nx, c, ny);
                    assert!(close(a[i * n + j], d + alpha * h), "A[{i},{j}]");
                    assert!(
                        close(b[i * n + j], s + alpha * dgdnx.scale(area)),
                        "B[{i},{j}]"
                    );
                    d0_rule += statics.1;
                    n0_rule += h0;
                } else {
                    assert_eq!(a[i * n + j], a0[i * n + j]);
                    assert_eq!(b[i * n + j], b0[i * n + j]);
                    d0_rule += point_d0.re * area;
                    n0_rule += point_n0.re * area;
                }
            }
            let (s_tri, n_tri) = triangle_self_terms(k, x, triangles[i]);
            let want_a = C64::from_re(-1.0 - d0_rule) + alpha * (n_tri - C64::from_re(n0_rule));
            assert!(close(a[i * n + i], want_a), "A[{i},{i}]");
            // B_ii = S_ii + alpha (D'_ii + 1/2), D' closed under its point rule.
            assert!(
                close(b[i * n + i], s_tri - alpha.scale(d0_point)),
                "B[{i},{i}]"
            );
        }
    }

    #[test]
    fn far_layers_match_subdivided_quadrature_beyond_the_near_zone() {
        // Targets at 4, 6 and 8 circumradii from the centroid along two
        // directions, the far side of the triangle_is_near boundary, with a
        // tilted target normal. The seven-point rule must hold every layer
        // within 1e-4 there, a hundredth of the 1e-3 acceptance band it
        // serves. Its static (k = 0) error must also fall with distance; the
        // oscillatory part keeps a kh-set floor that distance does not remove
        // (S at k = 1.5: 2.3e-7 at 4 radii, 6.1e-7 at 6). Measured on ts1:
        // N at 4 radii, k = 0, is 1.05e-5.
        let triangle = [[0.0, 0.0, 0.0], [1.0, 0.1, 0.0], [0.3, 0.9, 0.05]];
        let ny = unit(cross(
            sub(triangle[1], triangle[0]),
            sub(triangle[2], triangle[0]),
        ));
        let c = [0, 1, 2].map(|a| (triangle[0][a] + triangle[1][a] + triangle[2][a]) / 3.0);
        let radius = triangle
            .iter()
            .map(|&v| norm(sub(v, c)))
            .fold(0.0_f64, f64::max);
        let area = 0.5
            * norm(cross(
                sub(triangle[1], triangle[0]),
                sub(triangle[2], triangle[0]),
            ));
        let nx = unit([0.3, 0.2, 0.9]);
        // k = 1.5 puts this panel (sqrt area 0.66) at the solver's resolution
        // floor, six panels per wavelength: the far rule never sees a larger kh.
        for direction in [unit([0.6, 0.7, 0.4]), unit([-0.8, 0.2, -0.5])] {
            for k in [0.0, 1.5] {
                let mut previous = [f64::INFINITY; 5];
                for scale in [4.0, 6.0, 8.0] {
                    let x = [0, 1, 2].map(|a| c[a] + 1.0001 * scale * radius * direction[a]);
                    assert!(!triangle_is_near(x, c, triangle));
                    let (layers, statics) = triangle_far_layers(k, x, nx, (triangle, c, area), ny);
                    let reference = subdivided(k, x, nx, triangle, ny, 5);
                    let static_reference = subdivided(0.0, x, nx, triangle, ny, 5);
                    let cases = [
                        (layers[0], reference[0], "S"),
                        (layers[1], reference[1], "D"),
                        (layers[2], reference[2], "N"),
                        (C64::from_re(statics[0]), static_reference[1], "D0"),
                        (C64::from_re(statics[1]), static_reference[2], "N0"),
                    ];
                    for (slot, (got, want, what)) in cases.into_iter().enumerate() {
                        let error = (got - want).abs() / want.abs();
                        assert!(error < 1e-4, "{what} {scale}r k={k}: {error:e}");
                        assert!(
                            k > 0.0 || error < previous[slot],
                            "{what} {scale}r k={k}: {error:e} after {:e}",
                            previous[slot]
                        );
                        previous[slot] = error;
                    }
                }
            }
        }
    }

    #[test]
    fn far_rows_hold_their_sums_and_close_both_static_rows() {
        let surface = SpherePanels::icosphere(1.0, 1).unwrap();
        let triangles = surface.triangles().unwrap();
        let (k, alpha) = (1.7, C64::new(0.0, 1.0 / 1.7));
        let n = surface.areas().len();
        let rule = NearRule::build(&surface, k, usize::MAX, 1, &|| false).unwrap();
        let (mut a, mut b) = assemble_dense(&surface, k, alpha);
        rule.correct_dense(alpha, &mut a, &mut b);
        let (near_a, near_b) = (a.clone(), b.clone());
        correct_dense_far(&surface, k, alpha, &rule, &mut a, &mut b);
        let close = |got: C64, want: C64| (got - want).abs() <= 1e-13 * (1.0 + want.abs());
        for i in 0..n {
            let (x, nx) = (surface.centroids()[i], surface.normals()[i]);
            let (mut d0_delta, mut n0_delta) = (0.0, 0.0);
            for j in (0..n).filter(|&j| j != i) {
                let (c, ny, area) = (
                    surface.centroids()[j],
                    surface.normals()[j],
                    surface.areas()[j],
                );
                if triangle_is_near(x, c, triangles[j]) {
                    assert_eq!(
                        (a[i * n + j], b[i * n + j]),
                        (near_a[i * n + j], near_b[i * n + j])
                    );
                    continue;
                }
                let ([s, d, h], [d0, n0]) =
                    triangle_far_layers(k, x, nx, (triangles[j], c, area), ny);
                let (_, _, dgdnx, _) = kernels(k, x, nx, c, ny);
                let (_, point_d0, _, point_n0) = kernels(0.0, x, nx, c, ny);
                assert!(close(a[i * n + j], d + alpha * h), "A[{i},{j}]");
                assert!(
                    close(b[i * n + j], s + alpha * dgdnx.scale(area)),
                    "B[{i},{j}]"
                );
                d0_delta += d0 - point_d0.re * area;
                n0_delta += n0 - point_n0.re * area;
            }
            let want = near_a[i * n + i] - C64::from_re(d0_delta) - alpha.scale(n0_delta);
            assert!(close(a[i * n + i], want), "A[{i},{i}]");
            assert_eq!(b[i * n + i], near_b[i * n + i]);
        }
        // Declared workers leave every entry bitwise.
        let (mut wa, mut wb) = assemble_dense(&surface, k, alpha);
        rule.correct_dense(alpha, &mut wa, &mut wb);
        super::super::with_dense_workers(4, || {
            correct_dense_far(&surface, k, alpha, &rule, &mut wa, &mut wb);
        });
        assert_eq!((wa, wb), (a, b));
    }
}
