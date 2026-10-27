//! Multiply-connected free-traction Airy problem, from complementary energy.
//!
//! A free hole permits an independent affine Airy trace (value and two slopes).
//! Fixing that trace to the outer rim's gauge removes compatible stress fields.
//! The load is the weak geometric strain, including its boundary contribution;
//! merely adding three unknowns to the strong determinant load is incorrect.
use super::{Band, Csr, PlateSpec, add};
use std::collections::{BTreeMap, BTreeSet};

type Row = Vec<(usize, f64)>;

/// One symmetric quadratic term in a weak geometric-strain load. Node indices
/// are body-local. For rest height z, evaluate weight*(wa*wb+za*wb+wa*zb).
#[derive(Clone, Debug)]
pub struct MembraneTerm {
    pub a: usize,
    pub b: usize,
    pub weight: f64,
}

#[derive(Clone)]
struct GradientEntry {
    row: usize,
    other: usize,
    weight: f64,
}

/// Transpose only the fixed contribution addresses. Each output node retains
/// the original row/term order and the a-then-b order when a == b. Rest shape,
/// displacement, Airy force, and scale are read afresh for every evaluation.
#[derive(Clone)]
struct OrderedGradient {
    ptr: Vec<usize>,
    entries: Vec<GradientEntry>,
}

impl OrderedGradient {
    fn new(terms: &[Vec<MembraneTerm>], nodes: usize) -> Self {
        let mut rows: Vec<Vec<GradientEntry>> = (0..nodes).map(|_| Vec::new()).collect();
        for (row, terms) in terms.iter().enumerate() {
            for t in terms {
                rows[t.a].push(GradientEntry {
                    row,
                    other: t.b,
                    weight: t.weight,
                });
                rows[t.b].push(GradientEntry {
                    row,
                    other: t.a,
                    weight: t.weight,
                });
            }
        }
        let mut ptr = Vec::with_capacity(nodes + 1);
        let mut entries = Vec::with_capacity(rows.iter().map(Vec::len).sum());
        ptr.push(0);
        for row in rows {
            entries.extend(row);
            ptr.push(entries.len());
        }
        Self { ptr, entries }
    }

    fn add(&self, w: &[f64], z: &[f64], f: &[f64], h2: f64, g: &mut [f64]) {
        for (node, range) in self.ptr.windows(2).enumerate() {
            let mut value = g[node];
            for t in &self.entries[range[0]..range[1]] {
                // Preserve the original terms.iter().zip(f) short-force behavior.
                if let Some(&fc) = f.get(t.row) {
                    let scale = h2 * fc * t.weight;
                    value += scale * (w[t.other] + z[t.other]);
                }
            }
            g[node] = value;
        }
    }
}

#[derive(Clone)]
pub(super) struct AnnularMembrane {
    pub terms: Vec<Vec<MembraneTerm>>,
    gradient_rows: OrderedGradient,
    // The same geometric-strain product contributes to several Airy rows.
    // Keep their original term order, but evaluate each product once per step.
    strain_pairs: Vec<(usize, usize)>,
    strains: Vec<f64>,
    load_terms: Vec<Vec<(usize, f64)>>,
    core: Band,
    cross: [Vec<f64>; 3],
    lift: [Vec<f64>; 3],
    schur: Band,
    nc: usize,
}

fn outer(m: &mut [BTreeMap<usize, f64>], a: &Row, b: &Row, scale: f64) {
    for &(i, x) in a {
        for &(j, y) in b {
            add(m, i, j, scale * x * y);
        }
    }
}

fn quadratic(
    q: &mut [BTreeMap<(usize, usize), f64>],
    stress: &Row,
    dx: &Row,
    dy: &Row,
    scale: f64,
) {
    for &(c, v) in stress {
        for &(a, x) in dx {
            for &(b, y) in dy {
                *q[c].entry((a.min(b), a.max(b))).or_default() += scale * v * x * y;
            }
        }
    }
}

impl AnnularMembrane {
    pub(super) fn new(
        grid: &BTreeMap<(i64, i64), usize>,
        sorted: &[(i64, i64)],
        spec: &PlateSpec<'_>,
        compliance: impl Fn(f64, f64) -> Result<f64, String>,
    ) -> Result<(Self, Csr), String> {
        let (h, hole, radius, nu) = (spec.grid, spec.hole, spec.radius, spec.poisson);
        let ci: BTreeMap<_, _> = sorted.iter().enumerate().map(|(k, &p)| (p, k)).collect();
        let nc = sorted.len();
        let trace = |i: i64, j: i64| -> Row {
            if let Some(&k) = ci.get(&(i, j)) {
                return vec![(k, 1.0)];
            }
            let (x, y) = (i as f64 * h, j as f64 * h);
            if x.hypot(y) < (radius + hole) * 0.5 {
                vec![(nc, 1.0), (nc + 1, x), (nc + 2, y)]
            } else {
                vec![]
            }
        };
        let stencil = |i: i64, j: i64, offsets: &[(i64, i64, f64)]| -> Row {
            let mut row = BTreeMap::<usize, f64>::new();
            for &(di, dj, weight) in offsets {
                for (k, v) in trace(i + di, j + dj) {
                    *row.entry(k).or_default() += weight * v / (h * h);
                }
            }
            row.into_iter().filter(|&(_, v)| v != 0.0).collect()
        };
        let gradient = |i: i64, j: i64, di: i64, dj: i64| -> Result<Row, String> {
            let at = |d| grid.get(&(i + d * di, j + d * dj)).copied();
            let c = *grid
                .get(&(i, j))
                .ok_or("membrane gradient centre missing")?;
            Ok(match (at(-1), at(1)) {
                (Some(a), Some(b)) => vec![(a, -0.5 / h), (b, 0.5 / h)],
                (None, Some(b)) => match at(2) {
                    Some(d) => vec![(c, -1.5 / h), (b, 2.0 / h), (d, -0.5 / h)],
                    None => vec![(c, -1.0 / h), (b, 1.0 / h)],
                },
                (Some(a), None) => match at(-2) {
                    Some(d) => vec![(c, 1.5 / h), (a, -2.0 / h), (d, 0.5 / h)],
                    None => vec![(c, 1.0 / h), (a, -1.0 / h)],
                },
                _ => return Err("annulus grid cannot resolve a surface gradient".into()),
            })
        };
        let mut matrix = vec![BTreeMap::new(); nc + 3];
        let mut load = vec![BTreeMap::new(); nc + 3];
        for &(i, j) in grid.keys() {
            let xx = stencil(i, j, &[(-1, 0, 1.0), (0, 0, -2.0), (1, 0, 1.0)]);
            let yy = stencil(i, j, &[(0, -1, 1.0), (0, 0, -2.0), (0, 1, 1.0)]);
            if xx.is_empty() && yy.is_empty() {
                continue;
            }
            let weight = compliance(i as f64 * h, j as f64 * h)?;
            outer(&mut matrix, &xx, &xx, weight);
            outer(&mut matrix, &yy, &yy, weight);
            outer(&mut matrix, &xx, &yy, -nu * weight);
            outer(&mut matrix, &yy, &xx, -nu * weight);
            for (stress, di, dj) in [(&yy, 1, 0), (&xx, 0, 1)] {
                // Integrate the squared slope on both halves of the nodal
                // control interval. Squaring the centered average instead
                // erases alternating slopes before computing their strain.
                // With piecewise-linear displacement these are the exact
                // half-interval integrals, not an added stiffness or filter.
                let sides = (grid.get(&(i - di, j - dj)), grid.get(&(i + di, j + dj)));
                if let (Some(&a), Some(&b)) = sides {
                    let c = grid[&(i, j)];
                    let backward = vec![(a, -1.0 / h), (c, 1.0 / h)];
                    let forward = vec![(c, -1.0 / h), (b, 1.0 / h)];
                    quadratic(&mut load, stress, &backward, &backward, 0.25);
                    quadratic(&mut load, stress, &forward, &forward, 0.25);
                } else {
                    // Preserve the one-sided boundary reconstruction.
                    let slope = gradient(i, j, di, dj)?;
                    quadratic(&mut load, stress, &slope, &slope, 0.5);
                }
            }
        }
        let cells: BTreeSet<_> = sorted
            .iter()
            .flat_map(|&(i, j)| [(i - 1, j - 1), (i - 1, j), (i, j - 1), (i, j)])
            .collect();
        for (i, j) in cells {
            let xy = stencil(
                i,
                j,
                &[(0, 0, 1.0), (1, 0, -1.0), (0, 1, -1.0), (1, 1, 1.0)],
            );
            let mut nodes = [0; 4];
            for (n, (di, dj)) in nodes.iter_mut().zip([(0, 0), (1, 0), (0, 1), (1, 1)]) {
                *n = *grid
                    .get(&(i + di, j + dj))
                    .ok_or("annulus shear cell is incomplete")?;
            }
            let dx = nodes
                .into_iter()
                .zip([-0.5 / h, 0.5 / h, -0.5 / h, 0.5 / h])
                .collect();
            let dy = nodes
                .into_iter()
                .zip([-0.5 / h, -0.5 / h, 0.5 / h, 0.5 / h])
                .collect();
            let weight = 2.0 * (1.0 + nu) * compliance((i as f64 + 0.5) * h, (j as f64 + 0.5) * h)?;
            outer(&mut matrix, &xy, &xy, weight);
            quadratic(&mut load, &xy, &dx, &dy, -1.0);
        }
        // Keep the narrow band. Only the three boundary modes are dense;
        // eliminate them with the exact Schur complement, factored once.
        let interior: Vec<BTreeMap<usize, f64>> = matrix[..nc]
            .iter()
            .map(|r| r.range(..nc).map(|(&i, &v)| (i, v)).collect())
            .collect();
        let core = Band::factor(&interior)?;
        let cross: [Vec<f64>; 3] = std::array::from_fn(|k| {
            (0..nc)
                .map(|i| *matrix[i].get(&(nc + k)).unwrap_or(&0.0))
                .collect()
        });
        let lift: [Vec<f64>; 3] = std::array::from_fn(|k| {
            let mut x = cross[k].clone();
            core.solve(&mut x);
            x
        });
        let mut border = vec![BTreeMap::new(); 3];
        for (i, row) in border.iter_mut().enumerate() {
            for (j, column) in lift.iter().enumerate() {
                let reduction: f64 = cross[i].iter().zip(column).map(|(a, b)| a * b).sum();
                row.insert(
                    j,
                    matrix[nc + i].get(&(nc + j)).copied().unwrap_or(0.0) - reduction,
                );
            }
        }
        let schur = Band::factor(&border)?;
        let terms: Vec<Vec<MembraneTerm>> = load
            .into_iter()
            .map(|row| {
                row.into_iter()
                    .filter(|&(_, v)| v != 0.0)
                    .map(|((a, b), weight)| MembraneTerm { a, b, weight })
                    .collect()
            })
            .collect();
        let mut indices = BTreeMap::new();
        let mut strain_pairs = Vec::new();
        let load_terms = terms
            .iter()
            .map(|row| {
                row.iter()
                    .map(|term| {
                        let index = *indices.entry((term.a, term.b)).or_insert_with(|| {
                            let index = strain_pairs.len();
                            strain_pairs.push((term.a, term.b));
                            index
                        });
                        (index, term.weight)
                    })
                    .collect()
            })
            .collect();
        let strains = vec![0.0; strain_pairs.len()];
        let gradient_rows = OrderedGradient::new(&terms, grid.len());
        Ok((
            Self {
                terms,
                gradient_rows,
                strain_pairs,
                strains,
                load_terms,
                core,
                cross,
                lift,
                schur,
                nc,
            },
            Csr::from_rows(&matrix),
        ))
    }

    // CUDA plan construction also uses this unfused solve to export A^-1.
    // Keep it available outside tests even though CPU timesteps fuse the load.
    pub(super) fn solve(&self, f: &mut [f64]) {
        let (interior, boundary) = f.split_at_mut(self.nc);
        self.core.solve(interior);
        self.finish_solve(interior, boundary);
    }

    /// Build each load row during forward substitution, as the simply
    /// connected solver does. The Schur correction and summation order stay
    /// unchanged; independent strain work can overlap the triangular chain.
    pub(super) fn solve_with(&self, f: &mut [f64], mut rhs: impl FnMut(usize) -> f64) {
        let (interior, boundary) = f.split_at_mut(self.nc);
        self.core.solve_with(interior, &mut rhs);
        for (k, value) in boundary.iter_mut().enumerate() {
            *value = rhs(self.nc + k);
        }
        self.finish_solve(interior, boundary);
    }

    fn finish_solve(&self, interior: &mut [f64], boundary: &mut [f64]) {
        for (k, cross) in self.cross.iter().enumerate() {
            boundary[k] -= cross
                .iter()
                .zip(interior.iter())
                .map(|(a, b)| a * b)
                .sum::<f64>();
        }
        self.schur.solve(boundary);
        for (k, lift) in self.lift.iter().enumerate() {
            for (x, a) in interior.iter_mut().zip(lift) {
                *x -= a * boundary[k];
            }
        }
    }

    #[cfg(test)]
    fn load(&self, w: &[f64], z: &[f64], eh: f64, q: &mut [f64]) {
        for (c, qc) in q.iter_mut().enumerate() {
            *qc = self.load_row(c, w, z, eh);
        }
    }

    /// Refresh unconditionally: both displacement and rest geometry are public
    /// body state. No timestep/identity shortcut may reuse an older product.
    pub(super) fn prepare_load(&mut self, w: &[f64], z: &[f64]) {
        for (&(a, b), strain) in self.strain_pairs.iter().zip(&mut self.strains) {
            *strain = w[a] * w[b] + z[a] * w[b] + w[a] * z[b];
        }
    }

    #[inline]
    pub(super) fn prepared_load_row(&self, c: usize, eh: f64) -> f64 {
        eh * self.load_terms[c]
            .iter()
            .map(|&(index, weight)| weight * self.strains[index])
            .sum::<f64>()
    }

    #[cfg(test)]
    #[inline]
    pub(super) fn load_row(&self, c: usize, w: &[f64], z: &[f64], eh: f64) -> f64 {
        eh * self.terms[c]
            .iter()
            .map(|t| t.weight * (w[t.a] * w[t.b] + z[t.a] * w[t.b] + w[t.a] * z[t.b]))
            .sum::<f64>()
    }

    pub(super) fn gradient(&self, w: &[f64], z: &[f64], f: &[f64], h2: f64, g: &mut [f64]) {
        self.gradient_rows.add(w, z, f, h2, g);
    }

    /// Exact derivative J of the weak strain load q/(E t) at w=0.
    /// The quadratic displacement term vanishes; the two rest-shape cross
    /// terms remain. Assemble afresh for the supplied geometry, never reuse
    /// a matrix across a rest-shape change.
    pub(super) fn rest_jacobian(&self, z: &[f64]) -> Csr {
        let rows: Vec<_> = self
            .terms
            .iter()
            .map(|terms| {
                let mut row = BTreeMap::<usize, f64>::new();
                for t in terms {
                    *row.entry(t.a).or_default() += t.weight * z[t.b];
                    *row.entry(t.b).or_default() += t.weight * z[t.a];
                }
                row
            })
            .collect();
        Csr::from_rows(&rows)
    }

    /// Hessian of V = h²/(2 E t) q' A^-1 q at rest:
    /// B x = h² E t J' A^-1 J x. No finite displacement or fitted stiffness.
    pub(super) fn rest_hessian_apply(
        &self,
        j: &Csr,
        x: &[f64],
        eh: f64,
        h2: f64,
        airy: &mut [f64],
        out: &mut [f64],
    ) {
        airy.fill(0.0);
        j.mul_add(x, eh, airy);
        self.solve(airy);
        out.fill(0.0);
        for (c, &fc) in airy.iter().enumerate() {
            for k in j.ptr[c]..j.ptr[c + 1] {
                out[j.col[k]] += h2 * j.val[k] * fc;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::super::{Body, Edge, Nonlinear, PlateSpec};
    use super::{MembraneTerm, OrderedGradient};

    // Frozen original scatter, independent of the candidate's transposed rows.
    fn scatter(
        terms: &[Vec<MembraneTerm>],
        w: &[f64],
        z: &[f64],
        f: &[f64],
        h2: f64,
        g: &mut [f64],
    ) {
        for (row, &fc) in terms.iter().zip(f) {
            for t in row {
                let scale = h2 * fc * t.weight;
                g[t.a] += scale * (w[t.b] + z[t.b]);
                g[t.b] += scale * (w[t.a] + z[t.a]);
            }
        }
    }

    #[test]
    fn ordered_gradient_matches_real_annuli_with_live_field_changes() {
        for h in [0.08, 0.04, 0.02] {
            for hole in [0.1, 0.25] {
                for curved in [false, true] {
                    let body = plate(h, hole, curved);
                    let vk = body.vk.as_ref().unwrap();
                    let a = vk.annulus.as_ref().unwrap();
                    assert!(a.terms.iter().flatten().any(|t| t.a == t.b));
                    for (phase, rest_scale, gain) in [
                        (0.0, 1.0, 1.0),
                        (0.7, -2.0, -0.125),
                        (1.9, 0.0, 1e-100),
                        (0.0, 1.0, 1.0),
                    ] {
                        let w: Vec<_> = body
                            .xy
                            .iter()
                            .map(|[x, y]| {
                                0.007 * (19.0 * x + phase).sin() * (13.0 * y - phase).cos()
                            })
                            .collect();
                        let z: Vec<_> = body
                            .xy
                            .iter()
                            .zip(&vk.z0)
                            .map(|([x, y], z)| rest_scale * z + phase * 0.003 * (x - y))
                            .collect();
                        let f: Vec<_> = (0..a.terms.len())
                            .map(|i| gain * ((i % 19) as f64 - 9.0) / 7.0)
                            .collect();
                        let initial: Vec<_> = (0..body.n)
                            .map(|i| match i % 4 {
                                0 => -0.0,
                                1 => 0.0,
                                2 => 1e-6 * (i + 1) as f64,
                                _ => -1e3 / (i + 1) as f64,
                            })
                            .collect();
                        // The short-force control also preserves the original zip.
                        for count in [0, f.len() / 2, f.len()] {
                            let mut old = initial.clone();
                            let mut new = initial.clone();
                            scatter(&a.terms, &w, &z, &f[..count], h * h, &mut old);
                            a.gradient(&w, &z, &f[..count], h * h, &mut new);
                            assert!(old.iter().all(|v| v.is_finite()));
                            for (node, (old, new)) in old.iter().zip(&new).enumerate() {
                                assert_eq!(
                                    old.to_bits(),
                                    new.to_bits(),
                                    "h={h} hole={hole} curved={curved} phase={phase} count={count} node={node}"
                                );
                            }
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn ordered_gradient_keeps_diagonal_updates_and_initial_accumulator() {
        let terms = vec![vec![MembraneTerm {
            a: 0,
            b: 0,
            weight: 1.0,
        }]];
        let rows = OrderedGradient::new(&terms, 2);
        let initial = [9007199254740992.0_f64, -0.0];
        let mut old = initial;
        let mut new = initial;
        scatter(&terms, &[1.0], &[0.0], &[1.0], 1.0, &mut old);
        rows.add(&[1.0], &[0.0], &[1.0], 1.0, &mut new);
        assert_eq!(new.map(f64::to_bits), old.map(f64::to_bits));
        assert_eq!(new[0].to_bits(), initial[0].to_bits());
        assert_eq!(new[1].to_bits(), (-0.0_f64).to_bits());
        // Combining the diagonal's two updates, or summing before adding g,
        // changes rounding. This is a discriminating control, not a tolerance.
        assert_ne!((initial[0] + 2.0).to_bits(), old[0].to_bits());
    }

    #[test]
    fn ordered_gradient_keeps_original_row_order_under_cancellation() {
        let terms: Vec<_> = (0..3)
            .map(|_| {
                vec![MembraneTerm {
                    a: 0,
                    b: 1,
                    weight: 1.0,
                }]
            })
            .collect();
        let rows = OrderedGradient::new(&terms, 2);
        let f = [9007199254740992.0, -9007199254740992.0, 1.0];
        let mut old = [0.0; 2];
        let mut new = old;
        scatter(&terms, &[1.0; 2], &[0.0; 2], &f, 1.0, &mut old);
        rows.add(&[1.0; 2], &[0.0; 2], &f, 1.0, &mut new);
        assert_eq!(new.map(f64::to_bits), old.map(f64::to_bits));
        assert_eq!(new, [1.0; 2]);
        let mut reordered = [0.0; 2];
        scatter(
            &terms,
            &[1.0; 2],
            &[0.0; 2],
            &[f[0], f[2], f[1]],
            1.0,
            &mut reordered,
        );
        assert_ne!(new.map(f64::to_bits), reordered.map(f64::to_bits));
    }

    #[test]
    fn rest_tangent_matches_energy_derivative_and_is_reciprocal() {
        for curved in [false, true] {
            for h in [0.08, 0.04] {
                let mut b = plate(h, 0.1, curved);
                let x: Vec<_> =
                    b.xy.iter()
                        .map(|[x, y]| (11.0 * x + 7.0 * y).sin())
                        .collect();
                let y: Vec<_> =
                    b.xy.iter()
                        .map(|[x, y]| (5.0 * x - 13.0 * y).cos())
                        .collect();
                let vk = b.vk.as_ref().unwrap();
                let a = vk.annulus.as_ref().unwrap();
                let j = a.rest_jacobian(&vk.z0);
                let mut airy = vec![0.0; vk.f.len()];
                let mut bx = vec![0.0; b.n];
                let mut by = bx.clone();
                a.rest_hessian_apply(&j, &x, vk.eh, h * h, &mut airy, &mut bx);
                a.rest_hessian_apply(&j, &y, vk.eh, h * h, &mut airy, &mut by);
                let inner = |v: &[f64], w: &[f64]| v.iter().zip(w).map(|(a, b)| a * b).sum::<f64>();
                let norm = |v: &[f64]| inner(v, v).sqrt();
                assert!(inner(&x, &bx) >= -1e-12 * norm(&x) * norm(&bx));
                assert!(
                    (inner(&x, &by) - inner(&y, &bx)).abs()
                        <= 2e-10 * (norm(&x) * norm(&by) + norm(&y) * norm(&bx))
                );
                if !curved {
                    assert!(bx.iter().chain(&by).all(|&v| v == 0.0));
                }
                // Independent central difference of the actual nonlinear
                // gradient, not of the Jacobian being tested.
                let eps = 1e-7;
                for sign in [1.0, -1.0] {
                    for (w, x) in b.w.iter_mut().zip(&x) {
                        *w = sign * eps * x;
                    }
                    b.nonlinear_gradient();
                    for (value, g) in by.iter_mut().zip(&b.g) {
                        if sign > 0.0 {
                            *value = *g / (2.0 * eps);
                        } else {
                            *value -= *g / (2.0 * eps);
                        }
                    }
                }
                let error = bx
                    .iter()
                    .zip(&by)
                    .map(|(a, b)| (a - b).powi(2))
                    .sum::<f64>()
                    .sqrt();
                assert!(
                    error <= 2e-7 * norm(&bx).max(1.0),
                    "curved={curved} h={h}: error={error:e} norm={:e}",
                    norm(&bx)
                );
            }
        }
    }

    fn plate(h: f64, hole: f64, curved: bool) -> Body {
        Body::plate(
            "annular compatibility",
            &PlateSpec {
                radius: 0.5,
                hole,
                grid: h,
                edge: Edge::Free,
                thickness: &|r| if curved { 0.01 / (1.0 + r * r) } else { 0.01 },
                rest_height: &|r| if curved { 0.03 * r * r } else { 0.0 },
                young: 100.0,
                poisson: 0.33,
                density: 1000.0,
                tension: 0.0,
                bending: true,
                nonlinear: Nonlinear::VonKarman,
                t60: [(100.0, f64::INFINITY), (1000.0, f64::INFINITY)],
                support: None,
                muffle: None,
            },
        )
        .unwrap()
    }

    #[test]
    fn alternating_slopes_with_a_transverse_envelope_retain_strain() {
        // Alternating x slopes have nonzero squared strain. Its y envelope
        // makes it incompatible with an in-plane displacement, so it must
        // load the Airy problem even at an interior point. Averaging the two
        // slopes before squaring incorrectly makes this load vanish.
        let h = 0.02;
        let amplitude = 0.001;
        let transverse = 4.0;
        let body = plate(h, 0.1, false);
        let w: Vec<_> = body
            .xy
            .iter()
            .map(|[x, y]| {
                let sign = if (x / h).round() as i64 % 2 == 0 {
                    1.0
                } else {
                    -1.0
                };
                amplitude * sign * (transverse * y).cos()
            })
            .collect();
        let node = body.node_at(0.2, 0.0).unwrap();
        let vk = body.vk.as_ref().unwrap();
        let c = vk.nb.iter().position(|row| row[0] == Some(node)).unwrap();
        let mut q = vk.q.clone();
        vk.annulus.as_ref().unwrap().load(&w, &vk.z0, vk.eh, &mut q);
        let scale = vk.eh * amplitude * amplitude * transverse * transverse / (h * h);
        assert!(
            q[c].abs() > 0.1 * scale,
            "incompatible alternating-slope strain vanished: {} / {scale}",
            q[c]
        );
    }

    #[test]
    fn fused_load_solve_preserves_separate_load_and_schur_bits() {
        for h in [0.08, 0.04, 0.02] {
            for curved in [false, true] {
                let body = plate(h, 0.1, curved);
                let vk = body.vk.as_ref().unwrap();
                let annulus = vk.annulus.as_ref().unwrap();
                let mut q = vec![0.0; vk.q.len()];
                let mut fused_q = q.clone();
                let mut f = q.clone();
                let mut fused_f = q.clone();
                let mut prepared = annulus.clone();
                let mut prepared_q = q.clone();
                let mut prepared_f = q.clone();
                for phase in [0.0, 0.7, 1.9] {
                    let w: Vec<_> = body
                        .xy
                        .iter()
                        .map(|[x, y]| 0.007 * (19.0 * x + phase).sin() * (13.0 * y - phase).cos())
                        .collect();
                    annulus.load(&w, &vk.z0, vk.eh, &mut q);
                    f.copy_from_slice(&q);
                    annulus.solve(&mut f);
                    annulus.solve_with(&mut fused_f, |c| {
                        fused_q[c] = annulus.load_row(c, &w, &vk.z0, vk.eh);
                        fused_q[c]
                    });
                    prepared.prepare_load(&w, &vk.z0);
                    prepared.solve_with(&mut prepared_f, |c| {
                        prepared_q[c] = prepared.prepared_load_row(c, vk.eh);
                        prepared_q[c]
                    });
                    for (expected, actual) in
                        q.iter().chain(&f).zip(prepared_q.iter().chain(&prepared_f))
                    {
                        assert_eq!(
                            expected.to_bits(),
                            actual.to_bits(),
                            "prepared h={h} curved={curved} phase={phase}"
                        );
                    }
                    for (expected, actual) in q.iter().chain(&f).zip(fused_q.iter().chain(&fused_f))
                    {
                        assert_eq!(
                            expected.to_bits(),
                            actual.to_bits(),
                            "h={h} curved={curved} phase={phase}"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn prepared_load_refreshes_changed_rest_shape_and_displacement() {
        let body = plate(0.04, 0.1, true);
        let vk = body.vk.as_ref().unwrap();
        let original = vk.annulus.as_ref().unwrap();
        let mut prepared = original.clone();
        assert!(prepared.strains.len() < original.terms.iter().map(Vec::len).sum::<usize>());
        for (amplitude, rest_scale) in [
            (0.007, 1.0),
            (-0.003, 0.0),
            (0.0, -2.0),
            (1e-150, 1e-150),
            (0.007, 1.0),
        ] {
            let w: Vec<_> = body
                .xy
                .iter()
                .map(|[x, y]| amplitude * (19.0 * x + 13.0 * y).sin())
                .collect();
            let z: Vec<_> = vk.z0.iter().map(|z| rest_scale * z).collect();
            prepared.prepare_load(&w, &z);
            for c in 0..original.terms.len() {
                assert_eq!(
                    original.load_row(c, &w, &z, vk.eh).to_bits(),
                    prepared.prepared_load_row(c, vk.eh).to_bits(),
                    "row={c} amplitude={amplitude} rest_scale={rest_scale}"
                );
            }
        }
    }

    #[test]
    #[ignore = "paired performance measurement; run explicitly in release mode"]
    fn paired_annular_load_solve_timings() {
        use std::hint::black_box;
        use std::time::Instant;
        let body = plate(0.02, 0.1, true);
        let vk = body.vk.as_ref().unwrap();
        let annulus = vk.annulus.as_ref().unwrap();
        let w: Vec<_> = body
            .xy
            .iter()
            .map(|[x, y]| 0.007 * (19.0 * x + 13.0 * y).sin())
            .collect();
        let mut q = vec![0.0; vk.q.len()];
        let mut f = q.clone();
        // Warm both paths before 20 alternating pairs. Both do the same
        // load, core solve and three-boundary-mode correction on this body.
        for pair in 0..22 {
            let mut seconds = [0.0; 2];
            for path in if pair % 2 == 0 { [0, 1] } else { [1, 0] } {
                let start = Instant::now();
                for _ in 0..100 {
                    let w = black_box(w.as_slice());
                    if path == 0 {
                        annulus.load(w, &vk.z0, vk.eh, &mut q);
                        f.copy_from_slice(&q);
                        annulus.solve(&mut f);
                    } else {
                        annulus.solve_with(&mut f, |c| {
                            q[c] = annulus.load_row(c, w, &vk.z0, vk.eh);
                            q[c]
                        });
                    }
                    black_box((&q, &f));
                }
                seconds[path] = start.elapsed().as_secs_f64();
            }
            if pair >= 2 {
                println!(
                    "pair={} nodes={} separate_s={:.9} fused_s={:.9}",
                    pair - 2,
                    body.n,
                    seconds[0],
                    seconds[1]
                );
            }
        }
    }

    #[test]
    fn g1_free_annulus_energy_converges_to_radial_equilibrium() {
        // Independent plane-stress radial displacement solution for w=A*r^2,
        // with zero radial traction at both a and b. No engine matrix in oracle.
        let (b, amplitude, nu) = (0.5_f64, 0.1_f64, 0.33_f64);
        for a in [0.1_f64, 0.25] {
            let dr = (b - a) / 100_000.0;
            let mut exact = 0.0;
            for i in 0..100_000 {
                let r = a + (i as f64 + 0.5) * dr;
                let rr =
                    amplitude.powi(2) / 4.0 * (a * a + b * b - r * r - a * a * b * b / (r * r));
                let tt = amplitude.powi(2) / 4.0
                    * (a * a + b * b - 3.0 * r * r + a * a * b * b / (r * r));
                exact += (rr * rr + tt * tt - 2.0 * nu * rr * tt) * std::f64::consts::PI * r * dr;
            }
            let mut previous = f64::INFINITY;
            for h in [0.04, 0.02, 0.01] {
                let mut body = plate(h, a, false);
                for (w, &[x, y]) in body.w.iter_mut().zip(&body.xy) {
                    *w = amplitude * (x * x + y * y);
                }
                let computed = body.nonlinear_gradient();
                let error = (computed / exact - 1.0).abs();
                eprintln!(
                    "annulus a={a} h={h} energy={computed:e} exact={exact:e} relative_error={error}"
                );
                assert!(
                    error < 0.65 * previous,
                    "annulus energy must converge under refinement"
                );
                previous = error;
            }
            assert!(
                previous < 0.13,
                "finest staircase grid remains first-order at its rim"
            );
        }
    }

    #[test]
    fn g0_annular_force_is_energy_gradient_with_rest_shape_and_taper() {
        let mut body = plate(0.04, 0.1, true);
        for (w, &[x, y]) in body.w.iter_mut().zip(&body.xy) {
            *w = 0.007 * (3.0 * x + 2.0 * y).sin();
        }
        let energy = body.nonlinear_gradient();
        assert!(energy > 0.0);
        let expected = body.g.clone();
        for i in (0..body.n).step_by(17) {
            let original = body.w[i];
            let eps = 1e-7;
            body.w[i] = original + eps;
            let plus = body.nonlinear_gradient();
            body.w[i] = original - eps;
            let minus = body.nonlinear_gradient();
            body.w[i] = original;
            let actual = (plus - minus) / (2.0 * eps);
            assert!(
                (actual - expected[i]).abs() < 1e-11,
                "node {i}: {actual} != {}",
                expected[i]
            );
        }
        let vk = body.vk.as_ref().unwrap();
        let mut applied = vec![0.0; vk.f.len()];
        vk.matrix.mul_add(&vk.f, 1.0, &mut applied);
        let residual = applied
            .iter()
            .zip(&vk.q)
            .map(|(a, b)| (a - b).powi(2))
            .sum::<f64>()
            .sqrt();
        let norm = vk.q.iter().map(|q| q * q).sum::<f64>().sqrt();
        assert!(
            residual < 1e-8 * norm,
            "Schur solve residual {residual} / {norm}"
        );
    }

    #[test]
    fn g3_free_annular_rigid_tilt_has_no_membrane_energy() {
        for curved in [false, true] {
            let mut body = plate(0.02, 0.1, curved);
            for (w, &[x, y]) in body.w.iter_mut().zip(&body.xy) {
                *w = 0.03 * x - 0.07 * y + 0.02;
            }
            // A tilt's cross strain with the rest shape is compatible with
            // an in-plane displacement too; curvature must not pin rotation.
            assert!(body.nonlinear_gradient().abs() < 1e-20);
            assert!(body.g.iter().all(|g| g.abs() < 1e-16));
        }
    }
}
