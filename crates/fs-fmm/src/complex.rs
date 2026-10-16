//! Complex extension of the uniform-depth Chebyshev FMM. The original real
//! engine is reused only for topology/admission; its arithmetic is unchanged.

use std::collections::BTreeMap;

use fs_math::{c64::C64, det};

use crate::{
    CellKey, Fmm, FmmError, Laplace3d, MAX_COEFFICIENT_SLOTS, MAX_DIRECT_POINTS,
    MAX_INTERPOLATION_ORDER, MAX_TRANSLATION_WORK, MIN_INTERPOLATION_ORDER,
};

mod kernels;
pub use kernels::{
    ComplexKernel, HelmholtzAdjointDouble, HelmholtzDouble, HelmholtzGradient,
    HelmholtzHypersingular, HelmholtzSingle,
};

#[cfg(test)]
mod tests;

/// Non-directional qualification domain, with D the padded root-box side.
/// Higher kD must use the directional successor, not a silently fixed order.
pub const MAX_BB_HELMHOLTZ_KD: f64 = 10.0;
/// Same byte budget as the real path: each C64 occupies two f64 slots.
pub const MAX_COMPLEX_COEFFICIENT_SLOTS: usize = MAX_COEFFICIENT_SLOTS / 2;
const DIRECTIONAL_ROUTE: &str = "M2.5-W5b";

/// Conservative order recommendation within the battery's kD<=10, tol>=1e-6
/// domain. It is not a relative-error certificate for arbitrary cancelling
/// charges. Explicit lower orders are available for convergence measurements.
pub fn helmholtz_order_for(kd: f64, tol: f64) -> Result<usize, FmmError> {
    if !kd.is_finite() || kd < 0.0 || !tol.is_finite() || tol <= 0.0 || tol >= 1.0 {
        return Err(FmmError::ComplexRequest {
            reason: "order selection requires finite kD >= 0 and 0 < tolerance < 1",
        });
    }
    if kd > MAX_BB_HELMHOLTZ_KD || tol < 1.0e-6 {
        return Err(FmmError::HelmholtzRegime {
            kd,
            tolerance: Some(tol),
            route: DIRECTIONAL_ROUTE,
        });
    }
    Ok(MAX_INTERPOLATION_ORDER)
}

enum Operator<'k> {
    Scalar(&'k dyn ComplexKernel),
    Double(HelmholtzDouble),
    Adjoint(HelmholtzAdjointDouble),
    Mixed(HelmholtzHypersingular),
}

impl Operator<'_> {
    fn dimensions(&self) -> (usize, usize) {
        match self {
            Self::Scalar(_) => (1, 1),
            Self::Double(_) => (3, 1),
            Self::Adjoint(_) => (1, 3),
            Self::Mixed(_) => (3, 3),
        }
    }

    fn wavenumber(&self) -> Option<f64> {
        match self {
            Self::Scalar(k) => k.wavenumber(),
            Self::Double(k) => Some(k.k),
            Self::Adjoint(k) => Some(k.k),
            Self::Mixed(k) => Some(k.k),
        }
    }

    // Row-major [target component][source component]. Normals are deliberately
    // absent: source normals weight P2M charges; target normals contract L2P.
    fn block(&self, x: [f64; 3], y: [f64; 3]) -> [C64; 9] {
        let mut block = [C64::ZERO; 9];
        if let Self::Scalar(k) = self {
            block[0] = k.eval(x, y);
            return block;
        }
        let Some(v) = kernels::radial(self.wavenumber().expect("Helmholtz operator"), x, y) else {
            return block;
        };
        match self {
            Self::Double(_) => {
                for (a, slot) in block[..3].iter_mut().enumerate() {
                    *slot = v.gp.scale(-v.unit[a]);
                }
            }
            Self::Adjoint(_) => {
                for (a, slot) in block[..3].iter_mut().enumerate() {
                    *slot = v.gp.scale(v.unit[a]);
                }
            }
            Self::Mixed(_) => {
                let radial = v.gpp - v.gp.scale(v.inv_r);
                for a in 0..3 {
                    for b in 0..3 {
                        block[3 * a + b] = -radial.scale(v.unit[a] * v.unit[b]);
                        if a == b {
                            block[3 * a + b] = block[3 * a + b] - v.gp.scale(v.inv_r);
                        }
                    }
                }
            }
            Self::Scalar(_) => unreachable!("handled above"),
        }
        block
    }
}

/// Complex FMM on a shared source/target cloud, with optional per-point unit
/// normals. Real and complex engines use the same pair partition and tree.
pub struct FmmC64<'k> {
    geometry: Fmm<'static>,
    operator: Operator<'k>,
    normals: Option<Vec<[f64; 3]>>,
    nodes: Vec<f64>,
    barycentric: Vec<f64>,
    interactions: Vec<(CellKey, CellKey)>,
    coefficient_slots: usize,
}

impl<'k> FmmC64<'k> {
    /// Construct a scalar complex FMM with explicit interpolation order.
    /// Lower orders are diagnostic approximations, not tolerance certificates.
    pub fn new(
        kernel: &'k dyn ComplexKernel,
        points: Vec<[f64; 3]>,
        order: usize,
        leaf_cap: usize,
    ) -> Result<Self, FmmError> {
        Self::build(Operator::Scalar(kernel), points, None, order, leaf_cap)
    }

    /// Source-normal double layer, with normals factored into three P2M fields.
    pub fn new_double(
        kernel: HelmholtzDouble,
        points: Vec<[f64; 3]>,
        normals: Vec<[f64; 3]>,
        order: usize,
        leaf_cap: usize,
    ) -> Result<Self, FmmError> {
        Self::build(
            Operator::Double(kernel),
            points,
            Some(normals),
            order,
            leaf_cap,
        )
    }

    /// Target-normal adjoint double layer, contracting three local fields.
    pub fn new_adjoint_double(
        kernel: HelmholtzAdjointDouble,
        points: Vec<[f64; 3]>,
        normals: Vec<[f64; 3]>,
        order: usize,
        leaf_cap: usize,
    ) -> Result<Self, FmmError> {
        Self::build(
            Operator::Adjoint(kernel),
            points,
            Some(normals),
            order,
            leaf_cap,
        )
    }

    /// Mixed derivative. This supplies only the off-diagonal point kernel;
    /// BEM must still supply near-panel quadrature and finite-part self terms.
    pub fn new_hypersingular(
        kernel: HelmholtzHypersingular,
        points: Vec<[f64; 3]>,
        normals: Vec<[f64; 3]>,
        order: usize,
        leaf_cap: usize,
    ) -> Result<Self, FmmError> {
        Self::build(
            Operator::Mixed(kernel),
            points,
            Some(normals),
            order,
            leaf_cap,
        )
    }

    /// Validate a peak complex coefficient count before allocation. The bound
    /// includes both passes and transient cell buffers, in C64 units.
    pub fn validate_coefficient_slots(slots: usize) -> Result<(), FmmError> {
        if slots > MAX_COMPLEX_COEFFICIENT_SLOTS {
            return Err(FmmError::WorkEnvelopeExceeded {
                operation: "complex coefficient-slot (two f64 per C64)",
                requested: slots.saturating_mul(2),
                max: MAX_COEFFICIENT_SLOTS,
            });
        }
        Ok(())
    }

    #[allow(clippy::too_many_lines)] // one construction admission audit
    fn build(
        operator: Operator<'k>,
        points: Vec<[f64; 3]>,
        normals: Option<Vec<[f64; 3]>>,
        order: usize,
        leaf_cap: usize,
    ) -> Result<Self, FmmError> {
        if !(MIN_INTERPOLATION_ORDER..=MAX_INTERPOLATION_ORDER).contains(&order) {
            return Err(FmmError::InvalidOrder {
                order,
                min: MIN_INTERPOLATION_ORDER,
                max: MAX_INTERPOLATION_ORDER,
            });
        }
        if let Some(k) = operator.wavenumber() {
            if !k.is_finite() || k < 0.0 {
                return Err(FmmError::ComplexRequest {
                    reason: "wavenumber must be finite and >= 0",
                });
            }
        }
        if let Some(ns) = &normals {
            if ns.len() != points.len()
                || ns.iter().any(|n| {
                    n.iter().any(|v| !v.is_finite()) || (kernels::dot(*n, *n) - 1.0).abs() > 1.0e-12
                })
            {
                return Err(FmmError::ComplexRequest {
                    reason: "normal field must contain one finite unit vector per point (1e-12 squared-norm tolerance)",
                });
            }
        }
        // Reuse tree construction at its smallest coefficient order; it
        // allocates no coefficients. Actual complex storage and M2L work are
        // validated below before any complex pass. The real implementation's
        // conservative p^6 estimate remains untouched for existing consumers.
        let mut geometry = Fmm::new(&Laplace3d, points, MIN_INTERPOLATION_ORDER, leaf_cap)?;
        geometry.order = order;
        if let Some(k) = operator.wavenumber() {
            let kd = k * geometry.side;
            if !kd.is_finite() || kd > MAX_BB_HELMHOLTZ_KD {
                return Err(FmmError::HelmholtzRegime {
                    kd,
                    tolerance: None,
                    route: DIRECTIONAL_ROUTE,
                });
            }
        }
        let (ns, nt) = operator.dimensions();
        // Two extra cells bound parent clones and per-cell scratch arrays.
        let coefficient_slots = (geometry.cells.len() + 2)
            .saturating_mul(order.pow(3))
            .saturating_mul(ns + nt);
        Self::validate_coefficient_slots(coefficient_slots)?;
        let mut levels: BTreeMap<u32, Vec<CellKey>> = BTreeMap::new();
        for &key in geometry.cells.keys() {
            levels.entry(key.0).or_default().push(key);
        }
        let mut interactions = Vec::new();
        let work_per_pair = order.pow(6).saturating_mul(ns * nt);
        for (&level, keys) in &levels {
            if level == 0 {
                continue;
            }
            for &a in keys {
                for &b in keys {
                    let pa = (level - 1, a.1 >> 1, a.2 >> 1, a.3 >> 1);
                    let pb = (level - 1, b.1 >> 1, b.2 >> 1, b.3 >> 1);
                    if !Fmm::adjacent(a, b) && Fmm::adjacent(pa, pb) {
                        let requested = (interactions.len() + 1).saturating_mul(work_per_pair);
                        if requested > MAX_TRANSLATION_WORK {
                            return Err(FmmError::WorkEnvelopeExceeded {
                                operation: "complex M2L channel multiply-add",
                                requested,
                                max: MAX_TRANSLATION_WORK,
                            });
                        }
                        interactions.push((a, b));
                    }
                }
            }
        }
        let mut nodes = Vec::with_capacity(order);
        let mut barycentric = Vec::with_capacity(order);
        for j in 0..order {
            let angle = (2 * j + 1) as f64 * core::f64::consts::PI / (2 * order) as f64;
            let (s, c) = det::sin_cos(angle);
            nodes.push(-c);
            barycentric.push(if j % 2 == 0 { s } else { -s });
        }
        Ok(Self {
            geometry,
            operator,
            normals,
            nodes,
            barycentric,
            interactions,
            coefficient_slots,
        })
    }

    fn weights(&self, key: CellKey, point: [f64; 3]) -> Vec<f64> {
        let (lo, side) = self.geometry.cell_box(key);
        let basis = |axis: usize| {
            let x = ((point[axis] - lo[axis]) * (2.0 / side) - 1.0).clamp(-1.0, 1.0);
            let mut w = vec![0.0; self.nodes.len()];
            if let Some(j) = self.nodes.iter().position(|&n| (x - n).abs() < 1.0e-14) {
                w[j] = 1.0;
                return w;
            }
            let mut den = 0.0;
            for j in 0..w.len() {
                w[j] = self.barycentric[j] / (x - self.nodes[j]);
                den += w[j];
            }
            for v in &mut w {
                *v /= den;
            }
            w
        };
        let (wx, wy, wz) = (basis(0), basis(1), basis(2));
        let mut w = Vec::with_capacity(self.nodes.len().pow(3));
        for z in wz {
            for &y in &wy {
                for &x in &wx {
                    w.push(x * y * z);
                }
            }
        }
        w
    }

    fn grid(&self, key: CellKey) -> Vec<[f64; 3]> {
        let (lo, side) = self.geometry.cell_box(key);
        let mut out = Vec::with_capacity(self.nodes.len().pow(3));
        for &z in &self.nodes {
            for &y in &self.nodes {
                for &x in &self.nodes {
                    out.push([
                        lo[0] + 0.5 * side * (x + 1.0),
                        lo[1] + 0.5 * side * (y + 1.0),
                        lo[2] + 0.5 * side * (z + 1.0),
                    ]);
                }
            }
        }
        out
    }

    fn source_weight(&self, index: usize, channel: usize) -> f64 {
        if self.operator.dimensions().0 == 1 {
            1.0
        } else {
            self.normals.as_ref().expect("normal operator")[index][channel]
        }
    }

    fn target_weight(&self, index: usize, channel: usize) -> f64 {
        if self.operator.dimensions().1 == 1 {
            1.0
        } else {
            self.normals.as_ref().expect("normal operator")[index][channel]
        }
    }

    fn validate_charges(&self, charges: &[C64]) -> Result<(), FmmError> {
        let bad = charges.iter().position(|q| !finite(*q));
        if charges.len() != self.geometry.points.len() || bad.is_some() {
            return Err(FmmError::InvalidCharges {
                expected: self.geometry.points.len(),
                actual: charges.len(),
                non_finite_index: if charges.len() == self.geometry.points.len() {
                    bad
                } else {
                    None
                },
            });
        }
        Ok(())
    }

    /// P2M/M2M, M2L, L2L/L2P and the adjacent-leaf direct sum. Never invokes
    /// the quadratic oracle. Traversal and accumulation orders are fixed.
    #[allow(clippy::too_many_lines)] // the four standard passes share one channel layout
    pub fn potentials(&self, charges: &[C64]) -> Result<Vec<C64>, FmmError> {
        self.validate_charges(charges)?;
        let (ns, nt) = self.operator.dimensions();
        let np = self.nodes.len().pow(3);
        let mut multipoles: BTreeMap<CellKey, Vec<C64>> = BTreeMap::new();
        for (&key, cell) in self.geometry.cells.iter().rev() {
            let mut m = vec![C64::ZERO; np * ns];
            if cell.leaf {
                for &idx in &cell.points {
                    let weights = self.weights(key, self.geometry.points[idx]);
                    for a in 0..ns {
                        let q = charges[idx].scale(self.source_weight(idx, a));
                        for (j, &w) in weights.iter().enumerate() {
                            m[ns * j + a] = m[ns * j + a] + q.scale(w);
                        }
                    }
                }
            } else {
                for di in 0..2 {
                    for dj in 0..2 {
                        for dk in 0..2 {
                            let child = (key.0 + 1, 2 * key.1 + di, 2 * key.2 + dj, 2 * key.3 + dk);
                            let Some(cm) = multipoles.get(&child) else {
                                continue;
                            };
                            for (j, point) in self.grid(child).into_iter().enumerate() {
                                let weights = self.weights(key, point);
                                for (i, &w) in weights.iter().enumerate() {
                                    for a in 0..ns {
                                        m[ns * i + a] = m[ns * i + a] + cm[ns * j + a].scale(w);
                                    }
                                }
                            }
                        }
                    }
                }
            }
            multipoles.insert(key, m);
        }
        let mut locals: BTreeMap<CellKey, Vec<C64>> = self
            .geometry
            .cells
            .keys()
            .map(|&k| (k, vec![C64::ZERO; np * nt]))
            .collect();
        for &(a, b) in &self.interactions {
            let (ga, gb) = (self.grid(a), self.grid(b));
            let m = &multipoles[&b];
            let local = locals.get_mut(&a).expect("registered target cell");
            for (i, &x) in ga.iter().enumerate() {
                for (j, &y) in gb.iter().enumerate() {
                    let block = self.operator.block(x, y);
                    if block[..ns * nt].iter().any(|&z| !finite(z)) {
                        return Err(FmmError::NonFiniteOutput { index: i });
                    }
                    for ta in 0..nt {
                        for sa in 0..ns {
                            let slot = nt * i + ta;
                            local[slot] = local[slot] + block[ns * ta + sa] * m[ns * j + sa];
                        }
                    }
                }
            }
        }
        let mut out = vec![C64::ZERO; charges.len()];
        for (&key, cell) in &self.geometry.cells {
            if key.0 > 0 {
                let parent = (key.0 - 1, key.1 >> 1, key.2 >> 1, key.3 >> 1);
                let lp = locals[&parent].clone();
                let local = locals.get_mut(&key).expect("registered cell");
                for (i, point) in self.grid(key).into_iter().enumerate() {
                    let weights = self.weights(parent, point);
                    for (j, &w) in weights.iter().enumerate() {
                        for a in 0..nt {
                            local[nt * i + a] = local[nt * i + a] + lp[nt * j + a].scale(w);
                        }
                    }
                }
            }
            if !cell.leaf {
                continue;
            }
            let local = &locals[&key];
            for &idx in &cell.points {
                let weights = self.weights(key, self.geometry.points[idx]);
                for (j, &w) in weights.iter().enumerate() {
                    for a in 0..nt {
                        out[idx] =
                            out[idx] + local[nt * j + a].scale(w * self.target_weight(idx, a));
                    }
                }
            }
            for di in -1i64..=1 {
                for dj in -1i64..=1 {
                    for dk in -1i64..=1 {
                        let n = [
                            i64::from(key.1) + di,
                            i64::from(key.2) + dj,
                            i64::from(key.3) + dk,
                        ];
                        if n.iter().any(|&v| v < 0) {
                            continue;
                        }
                        let neighbor = (key.0, n[0] as u32, n[1] as u32, n[2] as u32);
                        let Some(source) = self.geometry.cells.get(&neighbor) else {
                            continue;
                        };
                        for &t in &cell.points {
                            for &s in &source.points {
                                if t != s {
                                    out[t] = out[t] + self.pair(t, s) * charges[s];
                                }
                            }
                        }
                    }
                }
            }
        }
        validate_output(out)
    }

    fn pair(&self, t: usize, s: usize) -> C64 {
        let block = self
            .operator
            .block(self.geometry.points[t], self.geometry.points[s]);
        let (ns, nt) = self.operator.dimensions();
        let mut value = C64::ZERO;
        for a in 0..nt {
            for b in 0..ns {
                value = value
                    + block[ns * a + b].scale(self.target_weight(t, a) * self.source_weight(s, b));
            }
        }
        value
    }

    /// Bounded O(N²) oracle, excluding equal indices. It is never used as a
    /// substitute for FMM potentials or to conceal an unqualified regime.
    pub fn direct(&self, charges: &[C64]) -> Result<Vec<C64>, FmmError> {
        self.validate_charges(charges)?;
        if charges.len() > MAX_DIRECT_POINTS {
            return Err(FmmError::WorkEnvelopeExceeded {
                operation: "complex direct pair",
                requested: charges.len().saturating_mul(charges.len()),
                max: MAX_DIRECT_POINTS * MAX_DIRECT_POINTS,
            });
        }
        let mut out = vec![C64::ZERO; charges.len()];
        for (t, slot) in out.iter_mut().enumerate() {
            for (s, &q) in charges.iter().enumerate() {
                if t != s {
                    *slot = *slot + self.pair(t, s) * q;
                }
            }
        }
        validate_output(out)
    }

    /// Root-box side in coordinate units, including the existing 1e-9 padding.
    #[must_use]
    pub fn root_box_size(&self) -> f64 {
        self.geometry.side
    }

    /// JSON tree/work statistics; coefficient slots are counted in C64 units.
    #[must_use]
    pub fn stats(&self) -> String {
        format!(
            "{{\"tree\":{},\"m2l_pairs\":{},\"complex_coefficient_slots\":{}}}",
            self.geometry.stats(),
            self.interactions.len(),
            self.coefficient_slots
        )
    }
}

fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}

fn validate_output(out: Vec<C64>) -> Result<Vec<C64>, FmmError> {
    if let Some(index) = out.iter().position(|&z| !finite(z)) {
        return Err(FmmError::NonFiniteOutput { index });
    }
    Ok(out)
}
