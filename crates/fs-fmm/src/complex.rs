//! Complex extension of the uniform-depth Chebyshev FMM. The original real
//! engine is reused only for topology/admission; its arithmetic is unchanged.

use std::collections::BTreeMap;

use fs_math::{c64::C64, det};

use crate::{
    CellKey, Fmm, FmmError, Laplace3d, MAX_COEFFICIENT_SLOTS, MAX_DIRECT_POINTS,
    MAX_INTERPOLATION_ORDER, MAX_TRANSLATION_WORK, MIN_INTERPOLATION_ORDER,
};

pub mod directional;
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

/// Caller-declared W5a resource ceilings. The plain constructors use
/// [`ComplexBudget::V1`]; the `_within` constructors take one explicitly.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ComplexBudget {
    /// Counted M2L channel multiply-adds in one pass. Passes poll no
    /// cancellation, so this also bounds one uncancellable pass.
    pub max_translation_work: usize,
}
impl ComplexBudget {
    /// The v1 ceiling, [`MAX_TRANSLATION_WORK`].
    pub const V1: Self = Self {
        max_translation_work: MAX_TRANSLATION_WORK,
    };
}

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
    // Two diagonal radial channels retain (x_axis-y_axis) exactly, instead
    // of spending one interpolation degree on this vanishing linear factor.
    Gradient(HelmholtzGradient),
    Double(HelmholtzDouble),
    Adjoint(HelmholtzAdjointDouble),
    Mixed(HelmholtzHypersingular),
}

impl Operator<'_> {
    fn dimensions(&self) -> (usize, usize) {
        match self {
            Self::Scalar(_) => (1, 1),
            Self::Gradient(_) => (2, 2),
            Self::Double(_) => (3, 1),
            Self::Adjoint(_) => (1, 3),
            Self::Mixed(_) => (3, 3),
        }
    }

    fn wavenumber(&self) -> Option<f64> {
        match self {
            Self::Scalar(k) => k.wavenumber(),
            Self::Gradient(k) => Some(k.k),
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
            Self::Gradient(_) => {
                let radial = v.gp.scale(v.inv_r);
                block[0] = radial;
                block[3] = radial;
            }
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
/// normals. The non-directional real/complex paths share their partition;
/// directional construction selects a parabolically separated pair tree.
pub struct FmmC64<'k> {
    geometry: Fmm<'static>,
    operator: Operator<'k>,
    normals: Option<Vec<[f64; 3]>>,
    nodes: Vec<f64>,
    barycentric: Vec<f64>,
    interactions: Vec<(CellKey, CellKey)>,
    /// Interaction indices sharing one M2L kernel block, in first-use order.
    /// Translation-invariant operators group by (level, cell offset); other
    /// kernels keep one group per pair.
    translations: Vec<Vec<usize>>,
    /// Optional adjacent-leaf pair values in `potentials` order.
    near: Option<Vec<C64>>,
    coefficient_slots: usize,
    directional: Option<directional::DirectionalPlan>,
    /// Declared worker threads for the non-directional passes; see
    /// [`Self::set_workers`].
    workers: usize,
    /// Optional per-group M2L kernel rows in `potentials` layout; see
    /// [`Self::cache_m2l_rows`].
    m2l_rows: Option<Vec<Vec<C64>>>,
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
        Self::new_within(kernel, points, order, leaf_cap, ComplexBudget::V1)
    }

    /// [`Self::new`] under a caller-declared [`ComplexBudget`].
    pub fn new_within(
        kernel: &'k dyn ComplexKernel,
        points: Vec<[f64; 3]>,
        order: usize,
        leaf_cap: usize,
        budget: ComplexBudget,
    ) -> Result<Self, FmmError> {
        Self::build(
            Operator::Scalar(kernel),
            points,
            None,
            order,
            leaf_cap,
            budget,
        )
    }

    /// Source-normal double layer, with normals factored into three P2M fields.
    pub fn new_double(
        kernel: HelmholtzDouble,
        points: Vec<[f64; 3]>,
        normals: Vec<[f64; 3]>,
        order: usize,
        leaf_cap: usize,
    ) -> Result<Self, FmmError> {
        Self::new_double_within(kernel, points, normals, order, leaf_cap, ComplexBudget::V1)
    }

    /// [`Self::new_double`] under a caller-declared [`ComplexBudget`].
    pub fn new_double_within(
        kernel: HelmholtzDouble,
        points: Vec<[f64; 3]>,
        normals: Vec<[f64; 3]>,
        order: usize,
        leaf_cap: usize,
        budget: ComplexBudget,
    ) -> Result<Self, FmmError> {
        Self::build(
            Operator::Double(kernel),
            points,
            Some(normals),
            order,
            leaf_cap,
            budget,
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
        Self::new_adjoint_double_within(kernel, points, normals, order, leaf_cap, ComplexBudget::V1)
    }

    /// [`Self::new_adjoint_double`] under a caller-declared [`ComplexBudget`].
    pub fn new_adjoint_double_within(
        kernel: HelmholtzAdjointDouble,
        points: Vec<[f64; 3]>,
        normals: Vec<[f64; 3]>,
        order: usize,
        leaf_cap: usize,
        budget: ComplexBudget,
    ) -> Result<Self, FmmError> {
        Self::build(
            Operator::Adjoint(kernel),
            points,
            Some(normals),
            order,
            leaf_cap,
            budget,
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
        Self::new_hypersingular_within(kernel, points, normals, order, leaf_cap, ComplexBudget::V1)
    }

    /// [`Self::new_hypersingular`] under a caller-declared [`ComplexBudget`].
    pub fn new_hypersingular_within(
        kernel: HelmholtzHypersingular,
        points: Vec<[f64; 3]>,
        normals: Vec<[f64; 3]>,
        order: usize,
        leaf_cap: usize,
        budget: ComplexBudget,
    ) -> Result<Self, FmmError> {
        Self::build(
            Operator::Mixed(kernel),
            points,
            Some(normals),
            order,
            leaf_cap,
            budget,
        )
    }

    /// Validate a peak complex coefficient count before allocation. The bound
    /// includes both passes and transient cell buffers, in C64 units.
    pub fn validate_coefficient_slots(slots: usize) -> Result<(), FmmError> {
        Self::validate_coefficient_slots_within(slots, MAX_COMPLEX_COEFFICIENT_SLOTS)
    }

    /// The same check against a caller-declared ceiling, also in C64 units.
    /// The error reports both counts in f64 slots, as the v1 check does.
    pub fn validate_coefficient_slots_within(slots: usize, max: usize) -> Result<(), FmmError> {
        if slots > max {
            return Err(FmmError::WorkEnvelopeExceeded {
                operation: "complex coefficient-slot (two f64 per C64)",
                requested: slots.saturating_mul(2),
                max: max.saturating_mul(2),
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
        budget: ComplexBudget,
    ) -> Result<Self, FmmError> {
        if budget.max_translation_work == 0 {
            return Err(FmmError::ComplexRequest {
                reason: "the declared M2L work budget must be positive",
            });
        }
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
                        if requested > budget.max_translation_work {
                            return Err(FmmError::WorkEnvelopeExceeded {
                                operation: "complex M2L channel multiply-add",
                                requested,
                                max: budget.max_translation_work,
                            });
                        }
                        interactions.push((a, b));
                    }
                }
            }
        }
        // Same-level grids are translates, so b-a fixes every node displacement.
        let invariant = match &operator {
            Operator::Scalar(k) => k.translation_invariant(),
            _ => true,
        };
        let mut translations: Vec<Vec<usize>> = Vec::new();
        let mut offsets: BTreeMap<(u32, i64, i64, i64), usize> = BTreeMap::new();
        for (index, &(a, b)) in interactions.iter().enumerate() {
            let offset = (
                a.0,
                i64::from(b.1) - i64::from(a.1),
                i64::from(b.2) - i64::from(a.2),
                i64::from(b.3) - i64::from(a.3),
            );
            let group = if invariant {
                *offsets.entry(offset).or_insert(translations.len())
            } else {
                translations.len()
            };
            if group == translations.len() {
                translations.push(Vec::new());
            }
            translations[group].push(index);
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
            translations,
            near: None,
            coefficient_slots,
            directional: None,
            workers: 1,
            m2l_rows: None,
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
        if let Operator::Gradient(kernel) = &self.operator {
            return if channel == 0 {
                1.0
            } else {
                self.geometry.points[index][kernel.axis] - self.geometry.lo[kernel.axis]
            };
        }
        if self.operator.dimensions().0 == 1 {
            1.0
        } else {
            self.normals.as_ref().expect("normal operator")[index][channel]
        }
    }

    fn target_weight(&self, index: usize, channel: usize) -> f64 {
        if let Operator::Gradient(kernel) = &self.operator {
            return if channel == 0 {
                self.geometry.points[index][kernel.axis] - self.geometry.lo[kernel.axis]
            } else {
                -1.0
            };
        }
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
        if let Some(plan) = &self.directional {
            return plan.potentials(self, charges);
        }
        if self.workers > 1 {
            return self.potentials_on_workers(charges);
        }
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
        // One target-node row of a group's kernel block is evaluated once and
        // applied to every pair in the group. The row fits the two scratch
        // cells already charged to the coefficient-slot budget. It is stored
        // per target channel in the multipole layout, so each application is
        // one contiguous dot product.
        let span = np * ns;
        let mut buffer = vec![C64::ZERO; nt * span];
        for (g, group) in self.translations.iter().enumerate() {
            let (a0, b0) = self.interactions[group[0]];
            let (ga, gb) = (self.grid(a0), self.grid(b0));
            for (i, &x) in ga.iter().enumerate() {
                let row: &[C64] = if let Some(rows) = &self.m2l_rows {
                    &rows[g][i * nt * span..(i + 1) * nt * span]
                } else {
                    for (j, &y) in gb.iter().enumerate() {
                        let block = self.operator.block(x, y);
                        if block[..ns * nt].iter().any(|&z| !finite(z)) {
                            return Err(FmmError::NonFiniteOutput { index: i });
                        }
                        for ta in 0..nt {
                            buffer[ta * span + ns * j..ta * span + ns * (j + 1)]
                                .copy_from_slice(&block[ns * ta..ns * (ta + 1)]);
                        }
                    }
                    &buffer
                };
                for &index in group {
                    let (a, b) = self.interactions[index];
                    let m = &multipoles[&b];
                    let local = locals.get_mut(&a).expect("registered target cell");
                    for ta in 0..nt {
                        let slot = nt * i + ta;
                        local[slot] = local[slot] + dot(&row[ta * span..(ta + 1) * span], m);
                    }
                }
            }
        }
        let mut out = vec![C64::ZERO; charges.len()];
        let mut next = 0;
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
            for source in self.neighbors(key) {
                for &t in &cell.points {
                    for &s in &source.points {
                        if t != s {
                            let value = match &self.near {
                                Some(cache) => {
                                    next += 1;
                                    cache[next - 1]
                                }
                                None => self.pair(t, s),
                            };
                            out[t] = out[t] + value * charges[s];
                        }
                    }
                }
            }
        }
        validate_output(out)
    }

    /// Declare `workers` threads (at least one) for the non-directional
    /// `potentials` passes. Potentials are bit-identical for every count:
    /// - upward and downward passes split each level's cells;
    /// - M2L splits target cells, and each worker applies the translation
    ///   groups in their fixed order to its own targets only;
    /// - leaves accumulate their own points, in the serial order.
    ///
    /// Each extra worker holds one M2L kernel row, charged to the
    /// coefficient-slot envelope. The directional route stays serial.
    ///
    /// # Errors
    /// `WorkEnvelopeExceeded` when the worker rows would exceed the
    /// coefficient-slot ceiling; the previous declaration is kept.
    pub fn set_workers(&mut self, workers: usize) -> Result<(), FmmError> {
        let workers = workers.max(1);
        if self.directional.is_none() {
            let (ns, nt) = self.operator.dimensions();
            let row = self.nodes.len().pow(3) * ns * nt;
            let base = self.coefficient_slots - (self.workers - 1) * row;
            let slots = base.saturating_add((workers - 1).saturating_mul(row));
            Self::validate_coefficient_slots(slots)?;
            self.coefficient_slots = slots;
        }
        self.workers = workers;
        Ok(())
    }

    /// [`Self::potentials`] on the declared workers. Every coefficient and
    /// output receives the serial path's operations in the serial order.
    #[allow(clippy::too_many_lines)] // the four passes share one channel layout
    fn potentials_on_workers(&self, charges: &[C64]) -> Result<Vec<C64>, FmmError> {
        let (ns, nt) = self.operator.dimensions();
        let np = self.nodes.len().pow(3);
        let workers = self.workers;
        let keys: Vec<CellKey> = self.geometry.cells.keys().copied().collect();
        let cells: Vec<&crate::Cell> = self.geometry.cells.values().collect();
        let index: BTreeMap<CellKey, usize> =
            keys.iter().enumerate().map(|(i, &k)| (k, i)).collect();
        // Keys sort by level first, so each level is one contiguous range.
        let mut levels: Vec<(usize, usize)> = Vec::new();
        for (i, key) in keys.iter().enumerate() {
            if levels
                .last()
                .is_some_and(|r: &(usize, usize)| keys[r.0].0 == key.0)
            {
                levels.last_mut().expect("nonempty").1 = i + 1;
            } else {
                levels.push((i, i + 1));
            }
        }
        // Upward: deepest level first; a cell reads only its children.
        let mut multipoles: Vec<Vec<C64>> = vec![Vec::new(); keys.len()];
        for &(start, end) in levels.iter().rev() {
            let (head, deeper) = multipoles.split_at_mut(end);
            let deeper: &[Vec<C64>] = deeper;
            let upward = |i: usize| -> Vec<C64> {
                let (key, cell) = (keys[i], cells[i]);
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
                                let child =
                                    (key.0 + 1, 2 * key.1 + di, 2 * key.2 + dj, 2 * key.3 + dk);
                                let Some(&c) = index.get(&child) else {
                                    continue;
                                };
                                let cm = &deeper[c - end];
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
                m
            };
            let fill = |i: usize, slot: &mut Vec<C64>| *slot = upward(i);
            split_run(&mut head[start..end], start, workers, &fill);
        }
        // M2L: workers own contiguous target-cell ranges, balanced by pairs.
        let mut incoming = vec![0usize; keys.len()];
        for &(a, _) in &self.interactions {
            incoming[index[&a]] += 1;
        }
        let total: usize = incoming.iter().sum();
        let mut bounds = vec![0usize];
        let mut running = 0usize;
        for (i, &count) in incoming.iter().enumerate() {
            running += count;
            let filled = bounds.len();
            if filled < workers && running * workers >= filled * total.max(1) {
                bounds.push(i + 1);
            }
        }
        while bounds.len() <= workers {
            bounds.push(keys.len());
        }
        let span = np * ns;
        let mut locals: Vec<Vec<C64>> = vec![vec![C64::ZERO; np * nt]; keys.len()];
        let multipoles = &multipoles;
        let m2l = |range: (usize, usize), owned: &mut [Vec<C64>]| -> Result<(), (usize, usize)> {
            let mut buffer = vec![C64::ZERO; nt * span];
            for (g, group) in self.translations.iter().enumerate() {
                let targets: Vec<(usize, usize)> = group
                    .iter()
                    .filter_map(|&pair| {
                        let (a, b) = self.interactions[pair];
                        let ta = index[&a];
                        (range.0..range.1).contains(&ta).then(|| (ta, index[&b]))
                    })
                    .collect();
                if targets.is_empty() {
                    continue;
                }
                let (a0, b0) = self.interactions[group[0]];
                let (ga, gb) = (self.grid(a0), self.grid(b0));
                for (i, &x) in ga.iter().enumerate() {
                    let row: &[C64] = if let Some(rows) = &self.m2l_rows {
                        &rows[g][i * nt * span..(i + 1) * nt * span]
                    } else {
                        for (j, &y) in gb.iter().enumerate() {
                            let block = self.operator.block(x, y);
                            if block[..ns * nt].iter().any(|&z| !finite(z)) {
                                return Err((g, i));
                            }
                            for ta in 0..nt {
                                buffer[ta * span + ns * j..ta * span + ns * (j + 1)]
                                    .copy_from_slice(&block[ns * ta..ns * (ta + 1)]);
                            }
                        }
                        &buffer
                    };
                    for &(ta_cell, b_cell) in &targets {
                        let m = &multipoles[b_cell];
                        let local = &mut owned[ta_cell - range.0];
                        for ta in 0..nt {
                            let slot = nt * i + ta;
                            local[slot] = local[slot] + dot(&row[ta * span..(ta + 1) * span], m);
                        }
                    }
                }
            }
            Ok(())
        };
        let mut first_error: Option<(usize, usize)> = None;
        {
            let mut rest = &mut locals[..];
            let mut parts = Vec::with_capacity(workers);
            for w in 0..workers {
                let (head, tail) =
                    core::mem::take(&mut rest).split_at_mut(bounds[w + 1] - bounds[w]);
                rest = tail;
                parts.push(((bounds[w], bounds[w + 1]), head));
            }
            let m2l = &m2l;
            let results: Vec<Result<(), (usize, usize)>> = std::thread::scope(|scope| {
                let handles: Vec<_> = parts
                    .into_iter()
                    .map(|(range, owned)| scope.spawn(move || m2l(range, owned)))
                    .collect();
                handles
                    .into_iter()
                    .map(|h| h.join().expect("M2L worker"))
                    .collect()
            });
            for result in results {
                if let Err(at) = result {
                    first_error = Some(first_error.map_or(at, |e: (usize, usize)| e.min(at)));
                }
            }
        }
        if let Some((_, i)) = first_error {
            return Err(FmmError::NonFiniteOutput { index: i });
        }
        // Downward: parents first; a cell reads only its finished parent.
        for &(start, end) in levels.iter().skip(1) {
            let (parents, rest) = locals.split_at_mut(start);
            let parents: &[Vec<C64>] = parents;
            split_run(
                &mut rest[..end - start],
                start,
                workers,
                &|c, local: &mut Vec<C64>| {
                    let key = keys[c];
                    let parent = (key.0 - 1, key.1 >> 1, key.2 >> 1, key.3 >> 1);
                    let lp = &parents[index[&parent]];
                    for (i, point) in self.grid(key).into_iter().enumerate() {
                        let weights = self.weights(parent, point);
                        for (j, &w) in weights.iter().enumerate() {
                            for a in 0..nt {
                                local[nt * i + a] = local[nt * i + a] + lp[nt * j + a].scale(w);
                            }
                        }
                    }
                },
            );
        }
        // Leaves: L2P then the adjacent-leaf sums, each point in its own leaf.
        let leaves: Vec<usize> = (0..keys.len()).filter(|&c| cells[c].leaf).collect();
        let mut offsets = Vec::with_capacity(leaves.len());
        let mut next = 0usize;
        for &c in &leaves {
            offsets.push(next);
            for source in self.neighbors(keys[c]) {
                next += cells[c].points.len() * source.points.len();
            }
            next -= cells[c].points.len();
        }
        let locals = &locals;
        let mut values: Vec<Vec<C64>> = vec![Vec::new(); leaves.len()];
        split_run(&mut values, 0, workers, &|l, out: &mut Vec<C64>| {
            let (c, mut next) = (leaves[l], offsets[l]);
            let (key, cell, local) = (keys[c], cells[c], &locals[c]);
            let mut sums = vec![C64::ZERO; cell.points.len()];
            for (p, &idx) in cell.points.iter().enumerate() {
                let weights = self.weights(key, self.geometry.points[idx]);
                for (j, &w) in weights.iter().enumerate() {
                    for a in 0..nt {
                        sums[p] = sums[p] + local[nt * j + a].scale(w * self.target_weight(idx, a));
                    }
                }
            }
            for source in self.neighbors(key) {
                for (p, &t) in cell.points.iter().enumerate() {
                    for &s in &source.points {
                        if t != s {
                            let value = match &self.near {
                                Some(cache) => {
                                    next += 1;
                                    cache[next - 1]
                                }
                                None => self.pair(t, s),
                            };
                            sums[p] = sums[p] + value * charges[s];
                        }
                    }
                }
            }
            *out = sums;
        });
        let mut out = vec![C64::ZERO; charges.len()];
        for (l, sums) in values.into_iter().enumerate() {
            for (&idx, value) in cells[leaves[l]].points.iter().zip(sums) {
                out[idx] = value;
            }
        }
        validate_output(out)
    }

    /// Same-level neighbour cells of `key`, itself included, in the fixed
    /// offset order shared by `potentials` and the near-field cache.
    fn neighbors(&self, key: CellKey) -> impl Iterator<Item = &crate::Cell> + '_ {
        (0..27i64).filter_map(move |o| {
            let n = [
                i64::from(key.1) + o / 9 - 1,
                i64::from(key.2) + (o / 3) % 3 - 1,
                i64::from(key.3) + o % 3 - 1,
            ];
            if n.iter().any(|&v| v < 0) {
                return None;
            }
            self.geometry
                .cells
                .get(&(key.0, n[0] as u32, n[1] as u32, n[2] as u32))
        })
    }

    /// Evaluate every translation group's M2L kernel rows once, on the
    /// declared workers, and reuse them on every later pass. They are the
    /// same `block` evaluations the passes make, so potentials are
    /// bit-identical with and without the cache. Returns its size in bytes.
    ///
    /// # Errors
    /// The engine stays uncached on every error:
    /// - `WorkEnvelopeExceeded` when the cache would exceed `max_bytes` (an
    ///   inclusive cap) or cannot be reserved;
    /// - `NonFiniteOutput`, the passes' own first error, for a non-finite
    ///   kernel row;
    /// - `ComplexRequest` on the directional route.
    pub fn cache_m2l_rows(&mut self, max_bytes: usize) -> Result<usize, FmmError> {
        if self.directional.is_some() {
            return Err(FmmError::ComplexRequest {
                reason: "the M2L row cache applies only to the non-directional route",
            });
        }
        let (ns, nt) = self.operator.dimensions();
        let np = self.nodes.len().pow(3);
        let span = np * ns;
        let per_group = np.saturating_mul(nt).saturating_mul(span);
        let bytes = self
            .translations
            .len()
            .saturating_mul(per_group)
            .saturating_mul(core::mem::size_of::<C64>());
        let refused = FmmError::WorkEnvelopeExceeded {
            operation: "complex M2L row cache bytes",
            requested: bytes,
            max: max_bytes,
        };
        if bytes > max_bytes {
            return Err(refused);
        }
        let mut slots: Vec<(Vec<C64>, Option<usize>)> = Vec::new();
        slots
            .try_reserve_exact(self.translations.len())
            .map_err(|_| refused.clone())?;
        for _ in 0..self.translations.len() {
            let mut rows = Vec::new();
            rows.try_reserve_exact(per_group)
                .map_err(|_| refused.clone())?;
            rows.resize(per_group, C64::ZERO);
            slots.push((rows, None));
        }
        let fill = |g: usize, slot: &mut (Vec<C64>, Option<usize>)| {
            let (a0, b0) = self.interactions[self.translations[g][0]];
            let (ga, gb) = (self.grid(a0), self.grid(b0));
            for (i, &x) in ga.iter().enumerate() {
                let row = &mut slot.0[i * nt * span..(i + 1) * nt * span];
                for (j, &y) in gb.iter().enumerate() {
                    let block = self.operator.block(x, y);
                    if block[..ns * nt].iter().any(|&z| !finite(z)) {
                        slot.1 = Some(i);
                        return;
                    }
                    for ta in 0..nt {
                        row[ta * span + ns * j..ta * span + ns * (j + 1)]
                            .copy_from_slice(&block[ns * ta..ns * (ta + 1)]);
                    }
                }
            }
        };
        split_run(&mut slots, 0, self.workers, &fill);
        if let Some(index) = slots.iter().find_map(|slot| slot.1) {
            return Err(FmmError::NonFiniteOutput { index });
        }
        self.m2l_rows = Some(slots.into_iter().map(|slot| slot.0).collect());
        Ok(bytes)
    }

    /// Evaluate every adjacent-leaf pair once, in the exact order
    /// `potentials` consumes them, and reuse the values on every later pass.
    /// They are the same `pair` evaluations, so potentials are bit-identical
    /// with and without the cache. Returns the cache size in bytes.
    ///
    /// # Errors
    /// `WorkEnvelopeExceeded` when the cache would exceed `max_bytes` or
    /// cannot be reserved; the engine then stays uncached. `ComplexRequest`
    /// on the directional route, which owns a different near partition.
    pub fn cache_near_field(&mut self, max_bytes: usize) -> Result<usize, FmmError> {
        if self.directional.is_some() {
            return Err(FmmError::ComplexRequest {
                reason: "the near-field cache applies only to the non-directional route",
            });
        }
        let mut count = 0usize;
        for (&key, cell) in self.geometry.cells.iter().filter(|(_, c)| c.leaf) {
            for source in self.neighbors(key) {
                count = count.saturating_add(cell.points.len().saturating_mul(source.points.len()));
            }
            // The cell is its own neighbour; its t == s pairs are excluded.
            count = count.saturating_sub(cell.points.len());
        }
        let bytes = count.saturating_mul(core::mem::size_of::<C64>());
        let refused = FmmError::WorkEnvelopeExceeded {
            operation: "complex near-field cache bytes",
            requested: bytes,
            max: max_bytes,
        };
        if bytes > max_bytes {
            return Err(refused);
        }
        let mut values = Vec::new();
        values.try_reserve_exact(count).map_err(|_| refused)?;
        for (&key, cell) in self.geometry.cells.iter().filter(|(_, c)| c.leaf) {
            for source in self.neighbors(key) {
                for &t in &cell.points {
                    for &s in &source.points {
                        if t != s {
                            values.push(self.pair(t, s));
                        }
                    }
                }
            }
        }
        self.near = Some(values);
        Ok(bytes)
    }

    fn pair(&self, t: usize, s: usize) -> C64 {
        // Keep direct/near evaluation on the original analytic point kernel;
        // it must not inherit cancellation in the far-field moment split.
        if let Operator::Gradient(kernel) = &self.operator {
            return kernel.eval(self.geometry.points[t], self.geometry.points[s]);
        }
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

    /// Nonempty direct-interaction leaf groups, in the matvec's fixed tree order.
    /// Every point occurs exactly once. A consumer may subdivide these groups
    /// for a bounded block-Jacobi preconditioner without rebuilding the tree.
    pub fn near_field_blocks(&self) -> impl Iterator<Item = &[usize]> {
        self.geometry
            .cells
            .values()
            .filter(|cell| cell.leaf)
            .map(|cell| cell.points.as_slice())
    }

    /// JSON tree/work statistics; coefficient slots are counted in C64 units.
    #[must_use]
    pub fn stats(&self) -> String {
        if let Some(plan) = &self.directional {
            return plan.stats.to_json();
        }
        format!(
            "{{\"tree\":{},\"m2l_pairs\":{},\"m2l_kernel_blocks\":{},\"complex_coefficient_slots\":{},\"near_cache_bytes\":{},\"m2l_cache_bytes\":{}}}",
            self.geometry.stats(),
            self.interactions.len(),
            self.translations.len(),
            self.coefficient_slots,
            self.near
                .as_ref()
                .map_or(0, |v| v.len() * core::mem::size_of::<C64>()),
            self.m2l_rows.as_ref().map_or(0, |rows| {
                rows.iter().map(Vec::len).sum::<usize>() * core::mem::size_of::<C64>()
            })
        )
    }
}

/// Run `f(offset + i, &mut slots[i])` for every slot, in contiguous chunks
/// on `workers` scoped threads (inline for one). Each call owns its slot.
fn split_run<T: Send>(
    slots: &mut [T],
    offset: usize,
    workers: usize,
    f: &(impl Fn(usize, &mut T) + Sync),
) {
    let parts = workers.min(slots.len());
    if parts <= 1 {
        for (i, slot) in slots.iter_mut().enumerate() {
            f(offset + i, slot);
        }
        return;
    }
    let len = slots.len();
    let mut rest = slots;
    std::thread::scope(|scope| {
        for p in 0..parts {
            let (start, end) = (p * len / parts, (p + 1) * len / parts);
            let (head, tail) = core::mem::take(&mut rest).split_at_mut(end - start);
            rest = tail;
            scope.spawn(move || {
                for (i, slot) in head.iter_mut().enumerate() {
                    f(offset + start + i, slot);
                }
            });
        }
    });
}

fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}

/// Complex dot product in four fixed interleaved partial sums with unfused
/// IEEE products. The order is fixed, so the result is bit-identical on every
/// ISA. The fused `C64` product calls libm `fma` per term when the target
/// lacks FMA, and that call dominated the M2L pass. The interpolation error
/// dwarfs the rounding difference.
fn dot(a: &[C64], b: &[C64]) -> C64 {
    let (mut re, mut im) = ([0.0f64; 4], [0.0f64; 4]);
    let (a4, b4) = (a.chunks_exact(4), b.chunks_exact(4));
    let (ar, br) = (a4.remainder(), b4.remainder());
    for (x, y) in a4.zip(b4) {
        for lane in 0..4 {
            re[lane] += x[lane].re * y[lane].re - x[lane].im * y[lane].im;
            im[lane] += x[lane].re * y[lane].im + x[lane].im * y[lane].re;
        }
    }
    for (x, y) in ar.iter().zip(br) {
        re[0] += x.re * y.re - x.im * y.im;
        im[0] += x.re * y.im + x.im * y.re;
    }
    C64::new(
        (re[0] + re[1]) + (re[2] + re[3]),
        (im[0] + im[1]) + (im[2] + im[3]),
    )
}

fn validate_output(out: Vec<C64>) -> Result<Vec<C64>, FmmError> {
    if let Some(index) = out.iter().position(|&z| !finite(z)) {
        return Err(FmmError::NonFiniteOutput { index });
    }
    Ok(out)
}

#[cfg(test)]
mod near_block_tests {
    use super::*;

    #[test]
    fn near_blocks_partition_indices_in_fixed_order() {
        let points = (0..40)
            .map(|i| [f64::from(i % 5), f64::from((i / 5) % 4), f64::from(i / 20)])
            .collect();
        let kernel = crate::HelmholtzSingle { k: 0.25 };
        let op = FmmC64::new(&kernel, points, 2, 4).unwrap();
        let first: Vec<_> = op.near_field_blocks().map(<[usize]>::to_vec).collect();
        assert_eq!(
            first,
            op.near_field_blocks()
                .map(<[usize]>::to_vec)
                .collect::<Vec<_>>()
        );
        assert!(first.len() > 1 && first.iter().all(|block| !block.is_empty()));
        let mut indices: Vec<_> = first.into_iter().flatten().collect();
        indices.sort_unstable();
        assert_eq!(indices, (0..40).collect::<Vec<_>>());
    }
}
