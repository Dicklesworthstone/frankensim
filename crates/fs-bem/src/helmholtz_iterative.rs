//! Matrix-free version of the existing centroid Burton-Miller discretization.
//!
//! Adjacent-leaf pairs use fs-fmm's exact point kernels. Self integrals and
//! jumps reuse the dense owner's equivalent-disc rule. Static row closures
//! retain the same point quadrature and effective interpolation order. The
//! directional route's zero-k owner uses the low-frequency partition.
//! FMM approximation error is distinct from the recomputed GMRES residual;
//! callers must qualify the explicit order against dense/refinement oracles.
//!
//! Plain CBIE keeps the dense owner's triangle quadrature. Its far field is a
//! far-field-only W5a pass over a seven-point quadrature cloud of every panel.
//! Each target row stores the dense owner's exact `(S, D)` entries for every
//! panel with a node in the FMM's direct-interaction leaves, minus the point
//! sums of those of its nodes the far pass already covers. A panel with all
//! nodes in separated leaves enters through its seven-point sum alone.

use super::triangle_rule::{NearRule, triangle_rule_pair};
use super::{
    Formulation, GL4_W, GL4_X, HelmholtzError, MAX_DENSE_PANELS, MIN_PANELS_PER_WAVELENGTH, Medium,
    PanelRule, RadiationSolution, alpha_for, characteristic_panel_size, evaluate_radiated_power,
    kernels, self_terms, surface_fingerprint, triangle_is_near, triangle_self_terms,
    triangle_weak_influence, triangle_weak_influence_with_rule,
};
use crate::{BemError, panel3d::SpherePanels};
use fs_exec::CancelGate;
use fs_fmm::{
    ComplexBudget, DirectionalOptions, FmmC64, FmmError, HelmholtzAdjointDouble, HelmholtzDouble,
    HelmholtzHypersingular, HelmholtzSingle, MAX_COMPLEX_COEFFICIENT_SLOTS, MAX_TRANSLATION_WORK,
};
use fs_la::eigen_complex::{LuComplex, lu_complex};
use fs_math::c64::C64;
use fs_solver::{
    krylov::{GmresState, SolveReport},
    op::LinearOp,
};
use std::{
    cell::{Cell, RefCell},
    time::Instant,
};

/// Right preconditioning leaves GMRES's residual in the original equation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Preconditioner {
    /// Original equation, without inverse application.
    None,
    /// Exact LU blocks drawn from individual direct-interaction FMM leaves.
    /// Large leaves are split in their stable original point order.
    NearFieldBlockJacobi { max_block_panels: usize },
}

/// Explicit accuracy, restart, FMM and near-field storage budgets.
#[derive(Clone, Debug)]
pub struct IterativeOptions {
    pub tol: f64,
    pub restart: usize,
    pub max_iters: usize,
    /// Explicit interpolation order; never silently increased or reduced.
    pub fmm_order: usize,
    pub leaf_capacity: usize,
    /// Explicit high-kD route. None retains W5a and its named kD refusal.
    /// When present, its order and leaf capacity are the effective FMM values.
    pub directional: Option<DirectionalOptions>,
    pub precond: Preconditioner,
    /// Maximum stored complex entries across all inverse-preconditioner
    /// blocks. Plain CBIE also charges its exact near table here, two complex
    /// entries (S and D) per stored target-panel pair.
    pub max_near_entries: usize,
    /// Byte budget for caching the double and mixed operators' adjacent-leaf
    /// pair values (half each), which are otherwise re-evaluated every matvec.
    /// The cache changes no bits. An operator whose cache exceeds its half,
    /// and the directional route, stay uncached. Zero disables it.
    pub near_field_cache_bytes: usize,
    /// Counted W5a M2L multiply-adds per operator pass; one pass is not
    /// cancellable. Defaults to the fs-fmm v1 `MAX_TRANSLATION_WORK`. The
    /// directional route declares its envelopes in its own options.
    pub max_translation_work: usize,
    /// Peak W5a complex coefficients per operator, in C64 units. Defaults to
    /// the fs-fmm v1 `MAX_COMPLEX_COEFFICIENT_SLOTS`. It changes no bits; the
    /// directional route declares its own.
    pub max_coefficient_slots: usize,
    /// Declared worker threads for every W5a operator pass
    /// ([`FmmC64::set_workers`]) and for the directional route's M2L and
    /// direct near passes. Potentials are bit-identical for every count. The
    /// directional coupling setup declares its own `setup_workers`. Defaults
    /// to one.
    pub workers: usize,
    /// Byte budget for caching the double and mixed operators' M2L kernel
    /// rows (half each), which are otherwise re-evaluated by every pass and
    /// by every declared worker. Plain CBIE's one far-field double layer may
    /// use all of it. The cache changes no bits. An operator whose cache
    /// exceeds its share, and the directional route, stay uncached.
    /// Zero (the default) disables it.
    pub m2l_cache_bytes: usize,
    /// Burton-Miller panel rule. [`PanelRule::TriangleNear`] adds the dense
    /// owner's triangle near rule as this owner's sparse near rows, built on
    /// the declared workers and refused above its declared pair count. The
    /// default [`PanelRule::Centroid`] is the centroid operator, bitwise.
    /// [`PanelRule::Triangle`] runs its far field on the Radon seven-point
    /// cloud (W5a only; directional options are refused) and is the dense
    /// Triangle operator to FMM accuracy. Plain CBIE ignores the rule, so one set of
    /// options can route rows of either formulation: its own near table is
    /// already the dense owner's triangle quadrature.
    pub panel_rule: PanelRule,
}

impl Default for IterativeOptions {
    fn default() -> Self {
        Self {
            tol: 1e-6,
            restart: 60,
            max_iters: 600,
            fmm_order: 4,
            leaf_capacity: 128,
            directional: None,
            precond: Preconditioner::NearFieldBlockJacobi {
                max_block_panels: 32,
            },
            max_near_entries: 4_194_304,
            near_field_cache_bytes: 1 << 29,
            max_translation_work: MAX_TRANSLATION_WORK,
            max_coefficient_slots: MAX_COMPLEX_COEFFICIENT_SLOTS,
            workers: 1,
            m2l_cache_bytes: 0,
            panel_rule: PanelRule::Centroid,
        }
    }
}

enum NormalLayer {
    Double,
    Adjoint,
    Mixed,
}
impl IterativeOptions {
    fn normal_fmm(
        &self,
        layer: NormalLayer,
        k: f64,
        points: Vec<[f64; 3]>,
        normals: Vec<[f64; 3]>,
    ) -> Result<FmmC64<'static>, FmmError> {
        let mut fmm = match (layer, &self.directional) {
            (NormalLayer::Double, Some(o)) => FmmC64::new_double_directional(
                HelmholtzDouble { k },
                points,
                normals,
                self.directional_for(o),
            ),
            (NormalLayer::Adjoint, Some(o)) => FmmC64::new_adjoint_double_directional(
                HelmholtzAdjointDouble { k },
                points,
                normals,
                self.directional_for(o),
            ),
            (NormalLayer::Mixed, Some(o)) => FmmC64::new_hypersingular_directional(
                HelmholtzHypersingular { k },
                points,
                normals,
                self.directional_for(o),
            ),
            (NormalLayer::Double, None) => FmmC64::new_double_within(
                HelmholtzDouble { k },
                points,
                normals,
                self.fmm_order,
                self.leaf_capacity,
                self.budget(),
            ),
            (NormalLayer::Adjoint, None) => FmmC64::new_adjoint_double_within(
                HelmholtzAdjointDouble { k },
                points,
                normals,
                self.fmm_order,
                self.leaf_capacity,
                self.budget(),
            ),
            (NormalLayer::Mixed, None) => FmmC64::new_hypersingular_within(
                HelmholtzHypersingular { k },
                points,
                normals,
                self.fmm_order,
                self.leaf_capacity,
                self.budget(),
            ),
        }?;
        fmm.set_workers(self.workers)?;
        Ok(fmm)
    }
    fn budget(&self) -> ComplexBudget {
        ComplexBudget {
            max_translation_work: self.max_translation_work,
            max_coefficient_slots: self.max_coefficient_slots,
        }
    }
    /// The directional options with this owner's W5a budget, which also
    /// governs a directional request that delegates whole to W5a (a static
    /// k = 0 operator, for one).
    fn directional_for(&self, o: &DirectionalOptions) -> DirectionalOptions {
        DirectionalOptions {
            w5a_budget: self.budget(),
            ..o.clone()
        }
    }
    fn effective_order(&self) -> usize {
        self.directional
            .as_ref()
            .map_or(self.fmm_order, |o| o.order)
    }
}

/// Measured process evidence, separate from numerical convergence authority.
#[derive(Clone, Debug)]
pub struct ResourceReceipt {
    pub host: String,
    /// Preparation, RHS and committed/aborted solve time for this state.
    pub wall_s: f64,
    /// Linux process VmHWM; None where that measurement is unavailable.
    /// This is process-wide high-water RSS, not a per-solve allocation claim.
    pub peak_rss_bytes: Option<u64>,
    pub matvecs: usize,
    pub iterations: usize,
    pub fmm_order: usize,
    pub near_field_bytes: usize,
    /// Bytes of cached adjacent-leaf pair values (zero when uncached).
    pub near_cache_bytes: usize,
    pub fmm_stats: [String; 2],
}

/// Refusals retain diagnostics; an unconverged pressure is never a solution.
#[derive(Clone, Debug)]
pub enum IterativeError {
    Geometry(BemError),
    Helmholtz(HelmholtzError),
    Fmm(FmmError),
    Cancelled {
        matvecs: usize,
    },
    Invalid {
        what: &'static str,
    },
    NearFieldBudget {
        required_entries: usize,
        maximum_entries: usize,
    },
    SingularNearBlock {
        block: usize,
    },
    NotConverged {
        pressure: Vec<C64>,
        report: SolveReport,
        receipt: ResourceReceipt,
    },
}
impl core::fmt::Display for IterativeError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Geometry(e) => e.fmt(f),
            Self::Helmholtz(e) => e.fmt(f),
            Self::Fmm(e) => e.fmt(f),
            Self::Cancelled { matvecs } => {
                write!(f, "FS-BEM-ITERATIVE-CANCELLED: after {matvecs} matvecs")
            }
            Self::Invalid { what } => write!(f, "FS-BEM-ITERATIVE-INVALID: {what}"),
            Self::NearFieldBudget {
                required_entries,
                maximum_entries,
            } => write!(
                f,
                "FS-BEM-ITERATIVE-NEAR-BUDGET: {required_entries} entries exceed {maximum_entries}"
            ),
            Self::SingularNearBlock { block } => {
                write!(f, "FS-BEM-ITERATIVE-SINGULAR-BLOCK: {block}")
            }
            Self::NotConverged { report, .. } => write!(
                f,
                "FS-BEM-ITERATIVE-NOT-CONVERGED: {} iterations, relative residual {}",
                report.iters, report.rel_residual
            ),
        }
    }
}
impl std::error::Error for IterativeError {}
impl From<HelmholtzError> for IterativeError {
    fn from(e: HelmholtzError) -> Self {
        Self::Helmholtz(e)
    }
}
impl From<FmmError> for IterativeError {
    fn from(e: FmmError) -> Self {
        Self::Fmm(e)
    }
}

struct NearBlock {
    indices: Vec<usize>,
    inverse: LuComplex,
}

/// Radon's seven-point degree-5 triangle rule (barycentric node, weight over
/// area). Its first node is the centroid, so the far-field cloud holds every
/// collocation point exactly.
const RADON_A: f64 = 0.101_286_507_323_456_34; // (6 - sqrt 15) / 21
const RADON_B: f64 = 0.797_426_985_353_087_3; // 1 - 2 RADON_A
const RADON_C: f64 = 0.470_142_064_105_115_1; // (6 + sqrt 15) / 21
const RADON_D: f64 = 0.059_715_871_789_769_82; // 1 - 2 RADON_C
const RADON_WA: f64 = 0.125_939_180_544_827_14; // (155 - sqrt 15) / 1200
const RADON_WC: f64 = 0.132_394_152_788_506_2; // (155 + sqrt 15) / 1200
pub(crate) const RADON7: [([f64; 3], f64); 7] = [
    ([1.0 / 3.0; 3], 0.225),
    ([RADON_A, RADON_A, RADON_B], RADON_WA),
    ([RADON_A, RADON_B, RADON_A], RADON_WA),
    ([RADON_B, RADON_A, RADON_A], RADON_WA),
    ([RADON_C, RADON_C, RADON_D], RADON_WC),
    ([RADON_C, RADON_D, RADON_C], RADON_WC),
    ([RADON_D, RADON_C, RADON_C], RADON_WC),
];

/// Smallest W5a leaf side, in panel radii (largest centroid-vertex distance),
/// at which a fully separated panel's seven-point sum stands in for the dense
/// owner's static-exact rule. Against a converged integral of 1/r and its
/// normal derivative over an equilateral triangle, the seven-point error
/// relative to the single layer is 1.1e-6 at a centroid distance of 4 radii,
/// 7.5e-8 at 6 and 1.3e-8 at 8; the dense 4x4 rule's is 2.9e-8 at 4.
pub const CBIE_MIN_LEAF_RADII: f64 = 4.0;

/// One boundary operator a near table carries. The adjoint double layer D'
/// always stays on its centroid rule.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum NearOperator {
    /// Single layer S.
    Single,
    /// Double layer D, with the source normal.
    Double,
    /// Hypersingular N, with both normals.
    Hypersingular,
}

/// A pluggable near-field rule for [`IterativeRadiation::prepare_with_near`]:
/// the sparse rows a Burton-Miller operator gains on top of its centroid
/// operator (P2P near field included), and its diagonal terms. Rows must be
/// finite, strictly ascending in column, and independent of the worker count.
pub(crate) trait NearRows: Sync {
    /// Operators per stored pair, in value order.
    fn operators(&self) -> &[NearOperator];
    /// Stored pairs of target row `i`, known before any allocation.
    fn row_len(&self, i: usize) -> usize;
    /// Ascending columns and pair-major values, `values[p * ops + o]`.
    fn fill_row(&self, i: usize, cols: &mut [u32], values: &mut [C64]);
    /// `(dA_ii, dB_ii)` at coupling `alpha`.
    fn diagonal(&self, i: usize, alpha: C64) -> (C64, C64);
    /// Rule declaration word for the checkpoint stamp.
    fn stamp(&self) -> u64;
}

/// This owner's one sparse near-field mechanism: target rows, ascending
/// columns, one value per declared operator per stored pair, applied
/// matrix-free beside the FMM passes. Plain CBIE stores its exact dense
/// triangle entries over a far-field-only pass; a Burton-Miller [`NearRows`]
/// rule stores corrections over the centroid operator.
struct NearTable {
    operators: Vec<NearOperator>,
    row_start: Vec<usize>,
    cols: Vec<u32>,
    values: Vec<C64>,
}

impl NearTable {
    /// Row offsets and the complex values the table would store, before the
    /// table is allocated.
    fn plan(
        n: usize,
        operators: usize,
        row_len: impl Fn(usize) -> usize,
    ) -> Result<(Vec<usize>, usize), IterativeError> {
        let overflow = IterativeError::Invalid {
            what: "near-field extent overflow",
        };
        if n > u32::MAX as usize || operators == 0 {
            return Err(overflow);
        }
        let mut row_start = Vec::with_capacity(n + 1);
        row_start.push(0usize);
        for i in 0..n {
            let end = row_start[i]
                .checked_add(row_len(i))
                .ok_or(overflow.clone())?;
            row_start.push(end);
        }
        let values = row_start[n].checked_mul(operators).ok_or(overflow)?;
        Ok((row_start, values))
    }

    /// Fill every row on the declared workers. Each row writes only its own
    /// range from shared data, so any worker count gives the same bits.
    fn build(
        operators: &[NearOperator],
        row_start: Vec<usize>,
        workers: usize,
        gate: &CancelGate,
        fill: &(dyn Fn(usize, &mut [u32], &mut [C64]) + Sync),
    ) -> Result<Self, IterativeError> {
        let n = row_start.len() - 1;
        let ops = operators.len();
        let pairs = row_start[n];
        let refused = IterativeError::Invalid {
            what: "near-field table allocation",
        };
        let mut cols = Vec::new();
        cols.try_reserve_exact(pairs).map_err(|_| refused.clone())?;
        cols.resize(pairs, 0u32);
        let mut values = Vec::new();
        values.try_reserve_exact(pairs * ops).map_err(|_| refused)?;
        values.resize(pairs * ops, C64::ZERO);
        let workers = workers.clamp(1, n.max(1));
        let cancelled = core::sync::atomic::AtomicBool::new(false);
        {
            let (mut c_rest, mut v_rest) = (&mut cols[..], &mut values[..]);
            std::thread::scope(|scope| {
                for w in 0..workers {
                    let (start, end) = (w * n / workers, (w + 1) * n / workers);
                    let entries = row_start[end] - row_start[start];
                    let (c_head, c_tail) = core::mem::take(&mut c_rest).split_at_mut(entries);
                    let (v_head, v_tail) = core::mem::take(&mut v_rest).split_at_mut(entries * ops);
                    (c_rest, v_rest) = (c_tail, v_tail);
                    let (row_start, cancelled) = (&row_start, &cancelled);
                    scope.spawn(move || {
                        let base = row_start[start];
                        for i in start..end {
                            if (i - start) % 256 == 0 && gate.is_requested() {
                                cancelled.store(true, core::sync::atomic::Ordering::Relaxed);
                                return;
                            }
                            let range = row_start[i] - base..row_start[i + 1] - base;
                            fill(
                                i,
                                &mut c_head[range.clone()],
                                &mut v_head[range.start * ops..range.end * ops],
                            );
                        }
                    });
                }
            });
        }
        if cancelled.into_inner() {
            return Err(IterativeError::Cancelled { matvecs: 0 });
        }
        for i in 0..n {
            let row = &cols[row_start[i]..row_start[i + 1]];
            if row.windows(2).any(|w| w[0] >= w[1]) || row.iter().any(|&j| j as usize >= n) {
                return Err(IterativeError::Invalid {
                    what: "near-field rows need strictly ascending in-range columns",
                });
            }
        }
        if values.iter().any(|v| !finite(*v)) {
            return Err(IterativeError::Invalid {
                what: "nonfinite near-field entry",
            });
        }
        Ok(Self {
            operators: operators.to_vec(),
            row_start,
            cols,
            values,
        })
    }

    /// `y_i = sum_j (sum over operators c_o T_o,ij) x_j`, rows on the declared
    /// workers. An undeclared operator contributes nothing.
    fn apply(&self, coefficients: &[(NearOperator, C64)], x: &[C64], workers: usize) -> Vec<C64> {
        let ops = self.operators.len();
        let terms: Vec<(usize, C64)> = self
            .operators
            .iter()
            .enumerate()
            .filter_map(|(o, op)| {
                coefficients
                    .iter()
                    .find(|(c, _)| c == op)
                    .map(|&(_, weight)| (o, weight))
            })
            .collect();
        let row = |i: usize| -> C64 {
            let mut sum = C64::ZERO;
            for p in self.row_start[i]..self.row_start[i + 1] {
                let entry = match terms.as_slice() {
                    [(o, weight)] if *weight == C64::ONE => self.values[p * ops + o],
                    _ => terms.iter().fold(C64::ZERO, |s, &(o, weight)| {
                        s + weight * self.values[p * ops + o]
                    }),
                };
                sum = sum + entry * x[self.cols[p] as usize];
            }
            sum
        };
        let n = self.row_start.len() - 1;
        let mut out = vec![C64::ZERO; n];
        let workers = workers.clamp(1, n.max(1));
        let mut rest = &mut out[..];
        std::thread::scope(|scope| {
            for w in 0..workers {
                let (start, end) = (w * n / workers, (w + 1) * n / workers);
                let (head, tail) = core::mem::take(&mut rest).split_at_mut(end - start);
                rest = tail;
                let row = &row;
                scope.spawn(move || {
                    for (slot, i) in head.iter_mut().zip(start..end) {
                        *slot = row(i);
                    }
                });
            }
        });
        out
    }

    /// The stored value of `op` at `(i, j)`, if the pair and operator exist.
    fn get(&self, i: usize, j: usize, op: NearOperator) -> Option<C64> {
        let o = self.operators.iter().position(|&d| d == op)?;
        let row = &self.cols[self.row_start[i]..self.row_start[i + 1]];
        let p = row.binary_search(&u32::try_from(j).ok()?).ok()?;
        Some(self.values[(self.row_start[i] + p) * self.operators.len() + o])
    }

    fn pairs(&self) -> usize {
        self.cols.len()
    }

    fn bytes(&self) -> usize {
        self.values.len() * core::mem::size_of::<C64>()
            + self.cols.len() * core::mem::size_of::<u32>()
            + self.row_start.len() * core::mem::size_of::<usize>()
    }
}

/// Plain CBIE's far-field FMM over the quadrature cloud, and its exact near
/// table (operators S and D). Cloud point j < N is panel j's centroid; node
/// m >= 1 of panel j is point N + 6j + m - 1.
struct CbieLayers {
    points: Vec<[f64; 3]>,
    /// Rule weight times panel area, per cloud point.
    weights: Vec<f64>,
    cloud_capacity: usize,
    double: FmmC64<'static>,
    table: NearTable,
}

/// The prepared layers: the original Burton-Miller FMM pair, with an optional
/// [`NearRows`] correction table, or plain CBIE.
enum Layers {
    BurtonMiller {
        double: FmmC64<'static>,
        mixed: FmmC64<'static>,
        /// The B row's single and adjoint double layers, built once: a Robin
        /// owner applies B on every matvec.
        single: FmmC64<'static>,
        adjoint: FmmC64<'static>,
        diagonal_a: Vec<C64>,
        diagonal_b: Vec<C64>,
        near: Option<NearTable>,
    },
    PlainCbie(Box<CbieLayers>),
    BurtonMillerCloud(Box<BmCloudLayers>),
}

/// Prepared immutable geometry/operator, reusable across boundary velocities.
/// Burton-Miller and its explicit sign-mutation arm use centroid point
/// kernels. Plain CBIE keeps the dense owner's triangle quadrature in its near
/// table and needs retained triangles; it is not a centroid-kernel operator.
pub struct IterativeRadiation<'a> {
    surface: &'a SpherePanels,
    admittance: Option<super::bc::Admittance<'a>>,
    k: f64,
    medium: Medium,
    alpha: C64,
    options: IterativeOptions,
    ppw: f64,
    layers: Layers,
    blocks: Vec<NearBlock>,
    near_entries: usize,
    near_cache_bytes: usize,
    setup_s: f64,
    stamp: u64,
}

/// Clone is a numerical checkpoint at a completed GMRES restart boundary.
/// Wall-time counters are diagnostic and do not participate in replay equality.
#[derive(Clone, Debug)]
pub struct IterativeState {
    gmres: GmresState,
    rhs: Vec<f64>,
    velocity: Vec<C64>,
    stamp: u64,
    matvecs: usize,
    wall_s: f64,
}
impl IterativeState {
    #[must_use]
    pub fn iterations(&self) -> usize {
        self.gmres.iters
    }
    #[must_use]
    pub fn report(&self, tolerance: f64) -> SolveReport {
        SolveReport::from_claim(
            self.gmres.iters,
            self.gmres.residual_claim(),
            tolerance,
            self.gmres.history.clone(),
        )
    }
}

/// Every owner's admission, in this order, before any engine is built:
/// geometry, frequency and medium, GMRES and preconditioner parameters, then
/// resolution. Returns the panels per wavelength.
pub(crate) fn admit(
    surface: &SpherePanels,
    k: f64,
    medium: Medium,
    options: &IterativeOptions,
) -> Result<f64, IterativeError> {
    surface.validate().map_err(IterativeError::Geometry)?;
    if !k.is_finite()
        || k <= 0.0
        || !medium.density.is_finite()
        || medium.density <= 0.0
        || !medium.sound_speed.is_finite()
        || medium.sound_speed <= 0.0
    {
        return Err(IterativeError::Invalid {
            what: "positive finite k, density and sound speed required",
        });
    }
    if !options.tol.is_finite()
        || options.tol <= 0.0
        || options.tol >= 1.0
        || options.restart == 0
        || options.restart > 256
    {
        return Err(IterativeError::Invalid {
            what: "0 < tol < 1 and 1 <= restart <= 256 required",
        });
    }
    if let Preconditioner::NearFieldBlockJacobi { max_block_panels } = options.precond {
        if !(1..=128).contains(&max_block_panels) {
            return Err(IterativeError::Invalid {
                what: "near blocks require 1..=128 panels",
            });
        }
    }
    let ppw = 2.0 * core::f64::consts::PI / (k * characteristic_panel_size(surface));
    if ppw < MIN_PANELS_PER_WAVELENGTH {
        return Err(HelmholtzError::TooCoarse {
            panels_per_wavelength: ppw,
        }
        .into());
    }
    Ok(ppw)
}

fn poll(gate: &CancelGate) -> Result<(), IterativeError> {
    if gate.is_requested() {
        Err(IterativeError::Cancelled { matvecs: 0 })
    } else {
        Ok(())
    }
}

impl<'a> IterativeRadiation<'a> {
    /// Robin extension of this same FMM/GMRES owner. The original A-block
    /// Jacobi inverse remains a right preconditioner; the recomputed residual
    /// is always that of A+ik*B*beta, including both BM rows.
    pub fn prepare_bc(
        surface: &'a SpherePanels,
        k: f64,
        medium: Medium,
        formulation: Formulation,
        options: &IterativeOptions,
        admittance: super::bc::Admittance<'a>,
        gate: &CancelGate,
    ) -> Result<Self, IterativeError> {
        poll(gate)?;
        if formulation == Formulation::PlainCbie {
            return Err(IterativeError::Invalid {
                what: "Robin boundaries keep the Burton-Miller iterative owner",
            });
        }
        admittance.validate(surface.areas().len())?;
        super::bc::validate_closed(surface)?;
        let mut out = Self::prepare(surface, k, medium, formulation, options, gate)?;
        if !admittance.is_zero() {
            out.stamp = (out.stamp ^ admittance.stamp()).wrapping_mul(0x100000001b3);
            out.admittance = Some(admittance);
        }
        Ok(out)
    }
    /// Prepare without allocating any N-by-N matrix. Near blocks are explicitly
    /// bounded, and FMM keeps its original coefficient/translation/P2P budgets.
    pub fn prepare(
        surface: &'a SpherePanels,
        k: f64,
        medium: Medium,
        formulation: Formulation,
        options: &IterativeOptions,
        gate: &CancelGate,
    ) -> Result<Self, IterativeError> {
        // A panel rule is a Burton-Miller rule; plain CBIE's own near table is
        // already the dense owner's triangle quadrature.
        let rule = if formulation == Formulation::PlainCbie {
            PanelRule::Centroid
        } else {
            options.panel_rule
        };
        let max_pairs = match rule {
            PanelRule::Centroid => {
                return Self::prepare_with_near(
                    surface,
                    k,
                    medium,
                    formulation,
                    options,
                    None,
                    gate,
                );
            }
            PanelRule::TriangleNear { max_pairs } => max_pairs,
            PanelRule::Triangle { max_pairs } => {
                return Self::prepare_cloud(
                    surface,
                    k,
                    medium,
                    formulation,
                    options,
                    max_pairs,
                    gate,
                );
            }
        };
        poll(gate)?;
        let rule = NearRule::build(surface, k, max_pairs, options.workers, &|| {
            gate.is_requested()
        })
        .map_err(|e| match e {
            HelmholtzError::Cancelled => IterativeError::Cancelled { matvecs: 0 },
            other => other.into(),
        })?;
        Self::prepare_with_near(surface, k, medium, formulation, options, Some(&rule), gate)
    }

    /// Burton-Miller under [`PanelRule::Triangle`] (see [`BmCloudLayers`]):
    /// the same operator as the dense rule, to FMM accuracy.
    fn prepare_cloud(
        surface: &'a SpherePanels,
        k: f64,
        medium: Medium,
        formulation: Formulation,
        options: &IterativeOptions,
        max_pairs: usize,
        gate: &CancelGate,
    ) -> Result<Self, IterativeError> {
        let started = Instant::now();
        poll(gate)?;
        let ppw = admit(surface, k, medium, options)?;
        let alpha = alpha_for(formulation, k);
        let (cloud, blocks, near_entries) =
            BmCloudLayers::prepare(surface, k, alpha, max_pairs, options, gate)?;
        let mut stamp = surface_fingerprint(surface);
        let preconditioner_key = match options.precond {
            Preconditioner::None => 0,
            Preconditioner::NearFieldBlockJacobi { max_block_panels } => max_block_panels as u64,
        };
        // 0xB3C7 marks the cloud owner; its table size fixes the partition.
        for word in [
            0xB3C7,
            k.to_bits(),
            alpha.im.to_bits(),
            medium.density.to_bits(),
            medium.sound_speed.to_bits(),
            options.tol.to_bits(),
            options.restart as u64,
            options.max_iters as u64,
            options.fmm_order as u64,
            options.leaf_capacity as u64,
            near_entries as u64,
            preconditioner_key,
            cloud.table.pairs() as u64,
        ] {
            stamp = (stamp ^ word).wrapping_mul(0x0000_0100_0000_01b3);
        }
        Ok(Self {
            surface,
            admittance: None,
            k,
            medium,
            alpha,
            options: options.clone(),
            ppw,
            layers: Layers::BurtonMillerCloud(Box::new(cloud)),
            blocks,
            near_entries,
            near_cache_bytes: 0,
            setup_s: started.elapsed().as_secs_f64(),
            stamp,
        })
    }

    /// [`Self::prepare`] with an optional Burton-Miller near-field rule. Its
    /// rows are corrections on top of the centroid operator: A gains
    /// `D + alpha N`, B gains `S`, the rule's diagonal terms join both
    /// diagonals, and the block-Jacobi blocks include all of them.
    /// `max_near_entries` charges the rule's stored values with the blocks.
    /// `None` is [`Self::prepare`], bitwise. Plain CBIE refuses a rule: its
    /// own table already holds the dense owner's triangle entries.
    pub(crate) fn prepare_with_near(
        surface: &'a SpherePanels,
        k: f64,
        medium: Medium,
        formulation: Formulation,
        options: &IterativeOptions,
        near: Option<&dyn NearRows>,
        gate: &CancelGate,
    ) -> Result<Self, IterativeError> {
        let started = Instant::now();
        poll(gate)?;
        let ppw = admit(surface, k, medium, options)?;
        if formulation == Formulation::PlainCbie {
            if near.is_some() {
                return Err(IterativeError::Invalid {
                    what: "plain CBIE already stores the dense triangle near entries",
                });
            }
            let (cbie, blocks, near_entries) = CbieLayers::prepare(surface, k, options, gate)?;
            let mut stamp = surface_fingerprint(surface);
            let preconditioner_key = match options.precond {
                Preconditioner::None => 0,
                Preconditioner::NearFieldBlockJacobi { max_block_panels } => {
                    max_block_panels as u64
                }
            };
            // 0xCB1E marks the formulation; Burton-Miller stamps are unchanged.
            for word in [
                0xCB1E,
                k.to_bits(),
                medium.density.to_bits(),
                medium.sound_speed.to_bits(),
                options.tol.to_bits(),
                options.restart as u64,
                options.max_iters as u64,
                options.fmm_order as u64,
                options.leaf_capacity as u64,
                near_entries as u64,
                preconditioner_key,
            ] {
                stamp = (stamp ^ word).wrapping_mul(0x0000_0100_0000_01b3);
            }
            return Ok(Self {
                surface,
                admittance: None,
                k,
                medium,
                alpha: C64::ZERO,
                options: options.clone(),
                ppw,
                layers: Layers::PlainCbie(Box::new(cbie)),
                blocks,
                near_entries,
                near_cache_bytes: 0,
                setup_s: started.elapsed().as_secs_f64(),
                stamp,
            });
        }
        let points = || surface.centroids().to_vec();
        let normals = || surface.normals().to_vec();
        let mut double = options.normal_fmm(NormalLayer::Double, k, points(), normals())?;
        poll(gate)?;
        let mut mixed = options.normal_fmm(NormalLayer::Mixed, k, points(), normals())?;
        let mut near_cache_bytes = 0;
        if options.directional.is_none() && options.near_field_cache_bytes > 0 {
            for op in [&mut double, &mut mixed] {
                poll(gate)?;
                match op.cache_near_field(options.near_field_cache_bytes / 2) {
                    Ok(bytes) => near_cache_bytes += bytes,
                    Err(FmmError::WorkEnvelopeExceeded { .. }) => {}
                    Err(e) => return Err(e.into()),
                }
            }
        }
        if options.directional.is_none() && options.m2l_cache_bytes > 0 {
            for op in [&mut double, &mut mixed] {
                poll(gate)?;
                match op.cache_m2l_rows(options.m2l_cache_bytes / 2) {
                    Ok(_) | Err(FmmError::WorkEnvelopeExceeded { .. }) => {}
                    Err(e) => return Err(e.into()),
                }
            }
        }
        let weights: Vec<_> = surface.areas().iter().map(|&a| C64::from_re(a)).collect();
        // Both static closures use the same point quadrature and effective
        // order. W5a also shares traversal; the directional zero-k call uses
        // the low-frequency owner's partition. No dense row scans occur.
        let d0 = {
            poll(gate)?;
            let op = options.normal_fmm(NormalLayer::Double, 0.0, points(), normals())?;
            op.potentials(&weights)?
        };
        let n0 = {
            poll(gate)?;
            let op = options.normal_fmm(NormalLayer::Mixed, 0.0, points(), normals())?;
            op.potentials(&weights)?
        };
        poll(gate)?;
        let alpha = alpha_for(formulation, k);
        let mut diagonal_a = Vec::with_capacity(weights.len());
        let mut diagonal_b = Vec::with_capacity(weights.len());
        for i in 0..weights.len() {
            let (s, n_reg) = self_terms(k, surface.areas()[i]);
            let d_self = C64::from_re(-0.5) - d0[i];
            diagonal_a.push(d_self - C64::from_re(0.5) + alpha * (n_reg - n0[i]));
            diagonal_b.push(s + alpha * (d_self + C64::from_re(0.5)));
        }
        let planned = match near {
            Some(rows) => Some(NearTable::plan(
                weights.len(),
                rows.operators().len(),
                |i| rows.row_len(i),
            )?),
            None => None,
        };
        let groups: Vec<Vec<usize>> = match options.precond {
            Preconditioner::None => Vec::new(),
            Preconditioner::NearFieldBlockJacobi { max_block_panels } => double
                .near_field_blocks()
                .flat_map(|leaf| leaf.chunks(max_block_panels).map(<[usize]>::to_vec))
                .collect(),
        };
        let near_entries = groups
            .iter()
            .try_fold(planned.as_ref().map_or(0, |p| p.1), |count, block| {
                block
                    .len()
                    .checked_mul(block.len())
                    .and_then(|entries| count.checked_add(entries))
            })
            .ok_or(IterativeError::Invalid {
                what: "near-field extent overflow",
            })?;
        if near_entries > options.max_near_entries {
            return Err(IterativeError::NearFieldBudget {
                required_entries: near_entries,
                maximum_entries: options.max_near_entries,
            });
        }
        let near_table = match (near, planned) {
            (Some(rows), Some((row_start, _))) => {
                poll(gate)?;
                let table = NearTable::build(
                    rows.operators(),
                    row_start,
                    options.workers,
                    gate,
                    &|i, cols, values| rows.fill_row(i, cols, values),
                )?;
                for i in 0..weights.len() {
                    let (da, db) = rows.diagonal(i, alpha);
                    if !finite(da) || !finite(db) {
                        return Err(IterativeError::Invalid {
                            what: "nonfinite near-field diagonal",
                        });
                    }
                    diagonal_a[i] = diagonal_a[i] + da;
                    diagonal_b[i] = diagonal_b[i] + db;
                }
                Some(table)
            }
            _ => None,
        };
        let correction = |i: usize, j: usize| -> Option<C64> {
            let table = near_table.as_ref()?;
            let d = table.get(i, j, NearOperator::Double);
            let h = table.get(i, j, NearOperator::Hypersingular);
            (d.is_some() || h.is_some())
                .then(|| d.unwrap_or(C64::ZERO) + alpha * h.unwrap_or(C64::ZERO))
        };
        let mut blocks = Vec::with_capacity(groups.len());
        for indices in groups {
            poll(gate)?;
            let m = indices.len();
            let mut matrix = vec![C64::ZERO; m * m];
            for (a, &i) in indices.iter().enumerate() {
                for (b, &j) in indices.iter().enumerate() {
                    matrix[a * m + b] = if i == j {
                        diagonal_a[i]
                    } else {
                        // This is the SAME exact P2P quadrature as dense BM and
                        // FMM's near field, not a second panel integration rule.
                        let (_, d, _, h) = kernels(
                            k,
                            surface.centroids()[i],
                            surface.normals()[i],
                            surface.centroids()[j],
                            surface.normals()[j],
                        );
                        (d + alpha * h).scale(surface.areas()[j])
                    };
                    // A near rule's stored pair enters the block as in apply.
                    if let Some(extra) = correction(i, j) {
                        matrix[a * m + b] = matrix[a * m + b] + extra;
                    }
                }
            }
            let inverse =
                lu_complex(&matrix, m).map_err(|_| IterativeError::SingularNearBlock {
                    block: blocks.len(),
                })?;
            blocks.push(NearBlock { indices, inverse });
        }
        poll(gate)?;
        let mut single = match &options.directional {
            Some(o) => FmmC64::new_single_directional(
                HelmholtzSingle { k },
                surface.centroids(),
                options.directional_for(o),
            )?,
            None => FmmC64::new_single_within(
                HelmholtzSingle { k },
                surface.centroids().to_vec(),
                options.fmm_order,
                options.leaf_capacity,
                options.budget(),
            )?,
        };
        single.set_workers(options.workers)?;
        poll(gate)?;
        let adjoint = options.normal_fmm(
            NormalLayer::Adjoint,
            k,
            surface.centroids().to_vec(),
            surface.normals().to_vec(),
        )?;
        // Cross-wiring guard for numerical checkpoints, not cryptographic proof.
        let mut stamp = surface_fingerprint(surface);
        let preconditioner_key = match options.precond {
            Preconditioner::None => 0,
            Preconditioner::NearFieldBlockJacobi { max_block_panels } => max_block_panels as u64,
        };
        for word in [
            k.to_bits(),
            alpha.im.to_bits(),
            medium.density.to_bits(),
            medium.sound_speed.to_bits(),
            options.tol.to_bits(),
            options.restart as u64,
            options.max_iters as u64,
            options.fmm_order as u64,
            options.leaf_capacity as u64,
            near_entries as u64,
            preconditioner_key,
        ] {
            stamp = (stamp ^ word).wrapping_mul(0x0000_0100_0000_01b3);
        }
        if let Some(o) = &options.directional {
            for word in [
                1,
                o.tolerance.to_bits(),
                o.order as u64,
                o.leaf_capacity as u64,
                o.low_freq_threshold.to_bits(),
                o.cone_resolution as u64,
                o.max_cones_per_level as u64,
                o.compression_tolerance.to_bits(),
                o.max_rank as u64,
            ] {
                stamp = (stamp ^ word).wrapping_mul(0x0000_0100_0000_01b3);
            }
        }
        if let (Some(rows), Some(table)) = (near, &near_table) {
            for word in [2, rows.stamp(), table.pairs() as u64] {
                stamp = (stamp ^ word).wrapping_mul(0x0000_0100_0000_01b3);
            }
        }
        Ok(Self {
            surface,
            admittance: None,
            k,
            medium,
            alpha,
            options: options.clone(),
            ppw,
            layers: Layers::BurtonMiller {
                double,
                mixed,
                single,
                adjoint,
                diagonal_a,
                diagonal_b,
                near: near_table,
            },
            blocks,
            near_entries,
            near_cache_bytes,
            setup_s: started.elapsed().as_secs_f64(),
            stamp,
        })
    }

    /// The shared Bq image; Robin application reuses this exact RHS owner.
    fn boundary_rhs(
        &self,
        velocity: &[C64],
        gate: &CancelGate,
    ) -> Result<Vec<C64>, IterativeError> {
        poll(gate)?;
        if velocity.len() != self.surface.areas().len() || velocity.iter().any(|v| !finite(*v)) {
            return Err(IterativeError::Invalid {
                what: "one finite complex normal velocity per panel required",
            });
        }
        let q: Vec<_> = velocity
            .iter()
            .map(|&v| v * C64::new(0.0, self.k * self.medium.sound_speed * self.medium.density))
            .collect();
        let (diagonal_b, near, single, adjoint) = match &self.layers {
            Layers::BurtonMiller {
                diagonal_b,
                near,
                single,
                adjoint,
                ..
            } => (diagonal_b, near, single, adjoint),
            Layers::PlainCbie(cbie) => {
                let rhs = cbie.single_layer(self.k, &self.options, &q, gate)?;
                if rhs.iter().any(|v| !finite(*v)) {
                    return Err(IterativeError::Invalid {
                        what: "nonfinite radiation right-hand side",
                    });
                }
                return Ok(rhs);
            }
            Layers::BurtonMillerCloud(cloud) => {
                let rhs =
                    cloud.apply_b(self.alpha, &q, self.surface.areas(), self.options.workers)?;
                poll(gate)?;
                if rhs.iter().any(|v| !finite(*v)) {
                    return Err(IterativeError::Invalid {
                        what: "nonfinite radiation right-hand side",
                    });
                }
                return Ok(rhs);
            }
        };
        let weighted: Vec<_> = q
            .iter()
            .zip(self.surface.areas())
            .map(|(&q, &a)| q.scale(a))
            .collect();
        let single = single.potentials(&weighted)?;
        poll(gate)?;
        let adjoint = adjoint.potentials(&weighted)?;
        poll(gate)?;
        let mut rhs_complex: Vec<_> = (0..q.len())
            .map(|i| single[i] + self.alpha * adjoint[i] + diagonal_b[i] * q[i])
            .collect();
        if let Some(table) = near {
            let extra = table.apply(
                &[(NearOperator::Single, C64::ONE)],
                &q,
                self.options.workers,
            );
            for (value, extra) in rhs_complex.iter_mut().zip(extra) {
                *value = *value + extra;
            }
        }
        if rhs_complex.iter().any(|v| !finite(*v)) {
            return Err(IterativeError::Invalid {
                what: "nonfinite radiation right-hand side",
            });
        }
        Ok(rhs_complex)
    }

    /// Compute Bq, with q=i*k*c*rho*v, using the original dense B signs.
    pub fn start(
        &self,
        velocity: &[C64],
        gate: &CancelGate,
    ) -> Result<IterativeState, IterativeError> {
        let started = Instant::now();
        let rhs = stack(&self.boundary_rhs(velocity, gate)?);
        let mut gmres = GmresState::new(&rhs, self.options.restart.min(rhs.len()));
        if rhs.iter().all(|&v| v == 0.0) {
            gmres.run(&ZeroOp(rhs.len()), &rhs, self.options.tol, 1, false);
        }
        Ok(IterativeState {
            gmres,
            rhs,
            velocity: velocity.to_vec(),
            stamp: self.stamp,
            matvecs: 0,
            wall_s: started.elapsed().as_secs_f64(),
        })
    }

    fn inverse(&self, input: &[C64]) -> Vec<C64> {
        if self.blocks.is_empty() {
            return input.to_vec();
        }
        let mut out = vec![C64::ZERO; input.len()];
        for block in &self.blocks {
            let mut values: Vec<_> = block.indices.iter().map(|&i| input[i]).collect();
            block.inverse.solve(&mut values);
            for (&i, value) in block.indices.iter().zip(values) {
                out[i] = value;
            }
        }
        out
    }

    /// Original unpreconditioned complex A application, including self closure.
    pub fn apply(&self, pressure: &[C64], gate: &CancelGate) -> Result<Vec<C64>, IterativeError> {
        poll(gate)?;
        if pressure.len() != self.surface.areas().len() {
            return Err(IterativeError::Invalid {
                what: "pressure/operator dimension mismatch",
            });
        }
        let mut out = match &self.layers {
            Layers::BurtonMiller {
                double,
                mixed,
                diagonal_a,
                near,
                ..
            } => {
                let weighted: Vec<_> = pressure
                    .iter()
                    .zip(self.surface.areas())
                    .map(|(&p, &a)| p.scale(a))
                    .collect();
                let d = double.potentials(&weighted)?;
                poll(gate)?;
                let h = mixed.potentials(&weighted)?;
                poll(gate)?;
                let mut out: Vec<_> = (0..pressure.len())
                    .map(|i| d[i] + self.alpha * h[i] + diagonal_a[i] * pressure[i])
                    .collect();
                if let Some(table) = near {
                    let extra = table.apply(
                        &[
                            (NearOperator::Double, C64::ONE),
                            (NearOperator::Hypersingular, self.alpha),
                        ],
                        pressure,
                        self.options.workers,
                    );
                    for (value, extra) in out.iter_mut().zip(extra) {
                        *value = *value + extra;
                    }
                }
                out
            }
            Layers::PlainCbie(cbie) => {
                let far = cbie.double.far_potentials(&cbie.charges(pressure))?;
                poll(gate)?;
                let mut out = cbie.table.apply(
                    &[(NearOperator::Double, C64::ONE)],
                    pressure,
                    self.options.workers,
                );
                for (i, value) in out.iter_mut().enumerate() {
                    *value = *value + far[i] - pressure[i].scale(0.5);
                }
                out
            }
            Layers::BurtonMillerCloud(cloud) => {
                let out = cloud.apply_a(self.alpha, pressure, self.options.workers)?;
                poll(gate)?;
                out
            }
        };
        if let Some(admittance) = self.admittance {
            let velocity: Vec<_> = admittance
                .image(pressure)?
                .into_iter()
                .map(|v| v.scale(1.0 / (self.medium.density * self.medium.sound_speed)))
                .collect();
            let b = self.boundary_rhs(&velocity, gate)?;
            for (a, b) in out.iter_mut().zip(b) {
                *a = *a + b;
            }
        }
        if out.iter().any(|v| !finite(*v)) {
            return Err(IterativeError::Invalid {
                what: "nonfinite boundary-integral application",
            });
        }
        Ok(out)
    }

    /// Run additional restart cycles; cancellation/refusal discards the current
    /// partial cycle and retains the preceding checkpoint. max_iters counts
    /// actual Arnoldi updates, including a final partial restart cycle.
    pub fn advance(
        &self,
        state: &mut IterativeState,
        cycles: usize,
        gate: &CancelGate,
    ) -> Result<SolveReport, IterativeError> {
        self.advance_controlled(state, cycles, gate, &mut |_| false)
    }

    /// The predicate is polled immediately before each actual matvec, with the
    /// cumulative matvec index. A true result cancels and is latched for this
    /// call. This supports deterministic cancellation without a timing race.
    pub fn advance_controlled(
        &self,
        state: &mut IterativeState,
        cycles: usize,
        gate: &CancelGate,
        cancel: &mut impl FnMut(usize) -> bool,
    ) -> Result<SolveReport, IterativeError> {
        if state.stamp != self.stamp {
            return Err(IterativeError::Invalid {
                what: "checkpoint belongs to a different prepared operator",
            });
        }
        let started = Instant::now();
        let op = RealOp {
            owner: self,
            gate,
            first_error: RefCell::new(None),
            matvecs: Cell::new(0),
            offset: state.matvecs,
            cancel: RefCell::new(cancel),
        };
        for _ in 0..cycles {
            if gate.is_requested() {
                state.matvecs += op.matvecs.get();
                state.wall_s += started.elapsed().as_secs_f64();
                return Err(IterativeError::Cancelled {
                    matvecs: state.matvecs,
                });
            }
            if state.gmres.rel_residual() < self.options.tol
                || state.gmres.iters >= self.options.max_iters
            {
                break;
            }
            let mut trial = state.gmres.clone();
            trial.restart = self
                .options
                .restart
                .min(self.options.max_iters - trial.iters)
                .min(state.rhs.len());
            trial.run(&op, &state.rhs, self.options.tol, 1, false);
            if gate.is_requested() && op.first_error.borrow().is_none() {
                *op.first_error.borrow_mut() = Some(IterativeError::Cancelled { matvecs: 0 });
            }
            if let Some(e) = op.first_error.borrow_mut().take() {
                state.matvecs += op.matvecs.get();
                state.wall_s += started.elapsed().as_secs_f64();
                return Err(match e {
                    IterativeError::Cancelled { .. } => IterativeError::Cancelled {
                        matvecs: state.matvecs,
                    },
                    other => other,
                });
            }
            trial.restart = self.options.restart.min(state.rhs.len());
            state.gmres = trial;
        }
        state.matvecs += op.matvecs.get();
        state.wall_s += started.elapsed().as_secs_f64();
        Ok(state.report(self.options.tol))
    }

    /// Publish only a converged original-equation residual and finite pressure.
    pub fn finish(
        &self,
        state: &IterativeState,
    ) -> Result<(RadiationSolution, SolveReport, ResourceReceipt), IterativeError> {
        if state.stamp != self.stamp {
            return Err(IterativeError::Invalid {
                what: "checkpoint belongs to a different prepared operator",
            });
        }
        let pressure = self.inverse(&unstack(&state.gmres.x));
        let report = state.report(self.options.tol);
        let receipt = self.receipt(state);
        if !report.converged_euclidean() || pressure.iter().any(|v| !finite(*v)) {
            return Err(IterativeError::NotConverged {
                pressure,
                report,
                receipt,
            });
        }
        let velocity = if let Some(admittance) = self.admittance {
            let absorbed = admittance.image(&pressure)?;
            state
                .velocity
                .iter()
                .zip(absorbed)
                .map(|(&v, y)| v - y.scale(1.0 / (self.medium.density * self.medium.sound_speed)))
                .collect()
        } else {
            state.velocity.clone()
        };
        let (radiated_power, radiated_power_roundoff_interval, power_diagnostics) =
            evaluate_radiated_power(&pressure, &velocity, self.surface.areas(), self.medium);
        Ok((
            RadiationSolution {
                pressure,
                velocity,
                k: self.k,
                medium: self.medium,
                surface_fingerprint: surface_fingerprint(self.surface),
                panels_per_wavelength: self.ppw,
                radiated_power,
                radiated_power_roundoff_interval,
                power_diagnostics,
                // Universal condition-number lower bound; no dense probe performed.
                condition_lower_bound: 1.0,
                dense_cap_utilization: self.surface.areas().len() as f64 / MAX_DENSE_PANELS as f64,
            },
            report,
            receipt,
        ))
    }

    fn near_table(&self) -> Option<&NearTable> {
        match &self.layers {
            Layers::BurtonMiller { near, .. } => near.as_ref(),
            Layers::PlainCbie(cbie) => Some(&cbie.table),
            Layers::BurtonMillerCloud(cloud) => Some(&cloud.table),
        }
    }

    #[must_use]
    pub fn receipt(&self, state: &IterativeState) -> ResourceReceipt {
        ResourceReceipt {
            host: host(),
            wall_s: self.setup_s + state.wall_s,
            peak_rss_bytes: peak_rss(),
            matvecs: state.matvecs,
            iterations: state.gmres.iters,
            fmm_order: self.options.effective_order(),
            near_field_bytes: match self.near_table() {
                None => self.near_entries * core::mem::size_of::<C64>(),
                Some(table) => {
                    (self.near_entries - table.values.len()) * core::mem::size_of::<C64>()
                        + table.bytes()
                }
            } + self
                .blocks
                .iter()
                .map(|b| b.indices.len() * core::mem::size_of::<usize>() * 2)
                .sum::<usize>(),
            near_cache_bytes: self.near_cache_bytes,
            fmm_stats: match &self.layers {
                Layers::BurtonMiller { double, mixed, .. } => [double.stats(), mixed.stats()],
                Layers::PlainCbie(cbie) => [
                    cbie.double.stats(),
                    format!(
                        "{{\"cbie_cloud_points\":{},\"cbie_near_pairs\":{}}}",
                        cbie.points.len(),
                        cbie.table.pairs()
                    ),
                ],
                Layers::BurtonMillerCloud(cloud) => [
                    cloud.mixed.stats(),
                    format!(
                        "{{\"bm_cloud_points\":{},\"bm_cloud_near_pairs\":{}}}",
                        cloud.weights.len(),
                        cloud.table.pairs()
                    ),
                ],
            },
        }
    }
}

/// One plain-CBIE dense entry `(S_ij, D_ij)`, by the dense owner's own
/// triangle dispatch (assemble_dense_with_static without a static table).
/// The A row's -1/2 jump is applied by the caller.
fn cbie_entry(
    k: f64,
    surface: &SpherePanels,
    triangles: &[[[f64; 3]; 3]],
    i: usize,
    j: usize,
) -> (C64, C64) {
    let x = surface.centroids()[i];
    if i == j {
        (triangle_self_terms(k, x, triangles[i]).0, C64::ZERO)
    } else if triangle_is_near(x, surface.centroids()[j], triangles[j]) {
        triangle_weak_influence(k, x, triangles[j], surface.normals()[j])
    } else {
        triangle_weak_influence_with_rule(k, x, triangles[j], surface.normals()[j], &GL4_X, &GL4_W)
    }
}

/// The Radon seven-point cloud of every panel, shared by the plain-CBIE and
/// Burton-Miller far fields: point `j < N` is panel j's centroid, node
/// `m >= 1` of panel j is point `N + 6j + m - 1`, and each weight is the rule
/// weight times the panel area. Also returns the largest centroid-vertex
/// distance.
fn radon_cloud(
    surface: &SpherePanels,
    triangles: &[[[f64; 3]; 3]],
) -> (Vec<[f64; 3]>, Vec<[f64; 3]>, Vec<f64>, f64) {
    let n = triangles.len();
    let mut points = vec![[0.0; 3]; n * RADON7.len()];
    let mut normals = vec![[0.0; 3]; points.len()];
    let mut weights = vec![0.0; points.len()];
    let mut radius = 0.0_f64;
    for (j, triangle) in triangles.iter().enumerate() {
        let centroid = surface.centroids()[j];
        for vertex in triangle {
            let d = [
                vertex[0] - centroid[0],
                vertex[1] - centroid[1],
                vertex[2] - centroid[2],
            ];
            radius = radius.max((d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt());
        }
        for (m, &(barycentric, weight)) in RADON7.iter().enumerate() {
            let index = if m == 0 { j } else { n + 6 * j + m - 1 };
            points[index] = if m == 0 {
                centroid
            } else {
                core::array::from_fn(|axis| {
                    barycentric[0] * triangle[0][axis]
                        + barycentric[1] * triangle[1][axis]
                        + barycentric[2] * triangle[2][axis]
                })
            };
            normals[index] = surface.normals()[j];
            weights[index] = weight * surface.areas()[j];
        }
    }
    (points, normals, weights, radius)
}

/// A far-field-only cloud engine's partition by panel: each point's leaf, each
/// leaf's direct-interaction leaves, and the panels with a node in them.
struct CloudPartition {
    leaves: Vec<Vec<usize>>,
    neighbours: Vec<Vec<usize>>,
    point_leaf: Vec<usize>,
    leaf_panels: Vec<Vec<usize>>,
}

impl CloudPartition {
    fn new(engine: &FmmC64<'_>, n: usize, points: usize) -> Self {
        let leaves: Vec<Vec<usize>> = engine.near_field_blocks().map(<[usize]>::to_vec).collect();
        let neighbours = engine.near_field_neighbor_leaves();
        let mut point_leaf = vec![0usize; points];
        for (l, leaf) in leaves.iter().enumerate() {
            for &p in leaf {
                point_leaf[p] = l;
            }
        }
        let panel_of = |p: usize| if p < n { p } else { (p - n) / 6 };
        let leaf_panels = neighbours
            .iter()
            .map(|around| {
                let mut panels: Vec<usize> = around
                    .iter()
                    .flat_map(|&l| leaves[l].iter().map(|&p| panel_of(p)))
                    .collect();
                panels.sort_unstable();
                panels.dedup();
                panels
            })
            .collect();
        Self {
            leaves,
            neighbours,
            point_leaf,
            leaf_panels,
        }
    }

    /// Panels with a node in a direct-interaction leaf of target `i`'s leaf,
    /// ascending.
    fn panels(&self, i: usize) -> &[usize] {
        &self.leaf_panels[self.point_leaf[i]]
    }

    /// Whether the far pass omits cloud point `p` at target `i`.
    fn is_near(&self, i: usize, p: usize) -> bool {
        self.neighbours[self.point_leaf[i]]
            .binary_search(&self.point_leaf[p])
            .is_ok()
    }

    /// Each leaf's centroid targets in stable order, split at `max_block_panels`.
    fn leaf_blocks(&self, n: usize, max_block_panels: usize) -> Vec<Vec<usize>> {
        self.leaves
            .iter()
            .flat_map(|leaf| {
                let targets: Vec<usize> = leaf.iter().copied().filter(|&p| p < n).collect();
                targets
                    .chunks(max_block_panels)
                    .map(<[usize]>::to_vec)
                    .collect::<Vec<_>>()
            })
            .filter(|block| !block.is_empty())
            .collect()
    }
}

impl CbieLayers {
    /// Build the cloud, its far-field double layer, the exact near table and
    /// the leaf-block preconditioner. Returns the stored complex entries.
    #[allow(clippy::too_many_lines)] // one admission and assembly audit
    fn prepare(
        surface: &SpherePanels,
        k: f64,
        options: &IterativeOptions,
        gate: &CancelGate,
    ) -> Result<(Self, Vec<NearBlock>, usize), IterativeError> {
        let triangles = surface.triangles().ok_or(IterativeError::Invalid {
            what: "iterative plain CBIE requires retained triangles",
        })?;
        if options.directional.is_some() {
            return Err(IterativeError::Invalid {
                what: "iterative plain CBIE is W5a only; directional rows are Burton-Miller",
            });
        }
        let n = triangles.len();
        let cloud_capacity =
            options
                .leaf_capacity
                .checked_mul(RADON7.len())
                .ok_or(IterativeError::Invalid {
                    what: "leaf capacity overflow",
                })?;
        let (points, normals, weights, radius) = radon_cloud(surface, triangles);
        poll(gate)?;
        let mut double = FmmC64::new_double_far_field_within(
            HelmholtzDouble { k },
            points.clone(),
            normals.clone(),
            options.fmm_order,
            cloud_capacity,
            options.budget(),
        )?;
        double.set_workers(options.workers)?;
        if double.leaf_box_size() < CBIE_MIN_LEAF_RADII * radius {
            return Err(IterativeError::Invalid {
                what: "iterative plain CBIE leaves must span CBIE_MIN_LEAF_RADII panel radii; raise leaf_capacity",
            });
        }
        let partition = CloudPartition::new(&double, n, points.len());
        let operators = [NearOperator::Single, NearOperator::Double];
        let (row_start, table_values) =
            NearTable::plan(n, operators.len(), |i| partition.panels(i).len())?;
        let groups: Vec<Vec<usize>> = match options.precond {
            Preconditioner::None => Vec::new(),
            Preconditioner::NearFieldBlockJacobi { max_block_panels } => {
                partition.leaf_blocks(n, max_block_panels)
            }
        };
        let near_entries = groups
            .iter()
            .try_fold(table_values, |count, block| {
                block
                    .len()
                    .checked_mul(block.len())
                    .and_then(|entries| count.checked_add(entries))
            })
            .ok_or(IterativeError::Invalid {
                what: "near-field extent overflow",
            })?;
        if near_entries > options.max_near_entries {
            return Err(IterativeError::NearFieldBudget {
                required_entries: near_entries,
                maximum_entries: options.max_near_entries,
            });
        }
        poll(gate)?;
        let fill = |i: usize, cols: &mut [u32], values: &mut [C64]| {
            let (x, nx) = (surface.centroids()[i], surface.normals()[i]);
            for (slot, &j) in partition.panels(i).iter().enumerate() {
                let (mut single, mut double) = cbie_entry(k, surface, triangles, i, j);
                // Remove the nodes the far pass already sums; the stored
                // entry plus that pass is the dense entry to FMM accuracy.
                for m in 0..RADON7.len() {
                    let p = if m == 0 { j } else { n + 6 * j + m - 1 };
                    if !partition.is_near(i, p) {
                        let (g, dgdny, _, _) = kernels(k, x, nx, points[p], normals[p]);
                        single = single - g.scale(weights[p]);
                        double = double - dgdny.scale(weights[p]);
                    }
                }
                cols[slot] = j as u32;
                values[2 * slot] = single;
                values[2 * slot + 1] = double;
            }
        };
        let table = NearTable::build(&operators, row_start, options.workers, gate, &fill)?;
        poll(gate)?;
        let mut blocks = Vec::with_capacity(groups.len());
        for indices in groups {
            poll(gate)?;
            let m = indices.len();
            let mut matrix = vec![C64::ZERO; m * m];
            for (a, &i) in indices.iter().enumerate() {
                for (b, &j) in indices.iter().enumerate() {
                    // j's centroid shares i's leaf, so every node of j is in
                    // a direct-interaction leaf and this is the dense entry.
                    let d =
                        table
                            .get(i, j, NearOperator::Double)
                            .ok_or(IterativeError::Invalid {
                                what: "near table lacks a same-leaf panel",
                            })?;
                    matrix[a * m + b] = d - if a == b { C64::from_re(0.5) } else { C64::ZERO };
                }
            }
            let inverse =
                lu_complex(&matrix, m).map_err(|_| IterativeError::SingularNearBlock {
                    block: blocks.len(),
                })?;
            blocks.push(NearBlock { indices, inverse });
        }
        if options.m2l_cache_bytes > 0 {
            poll(gate)?;
            match double.cache_m2l_rows(options.m2l_cache_bytes) {
                Ok(_) | Err(FmmError::WorkEnvelopeExceeded { .. }) => {}
                Err(e) => return Err(e.into()),
            }
        }
        Ok((
            Self {
                points,
                weights,
                cloud_capacity,
                double,
                table,
            },
            blocks,
            near_entries,
        ))
    }

    /// Panel densities spread to their cloud nodes, with node weights.
    fn charges(&self, density: &[C64]) -> Vec<C64> {
        let n = density.len();
        (0..self.points.len())
            .map(|p| density[if p < n { p } else { (p - n) / 6 }].scale(self.weights[p]))
            .collect()
    }

    /// The dense owner's S q: the near rows plus a far-field single layer on
    /// the same cloud and leaf capacity, hence the same tree and partition.
    fn single_layer(
        &self,
        k: f64,
        options: &IterativeOptions,
        q: &[C64],
        gate: &CancelGate,
    ) -> Result<Vec<C64>, IterativeError> {
        let kernel = HelmholtzSingle { k };
        let mut single = FmmC64::new_far_field_within(
            &kernel,
            self.points.clone(),
            options.fmm_order,
            self.cloud_capacity,
            options.budget(),
        )?;
        single.set_workers(options.workers)?;
        let far = single.far_potentials(&self.charges(q))?;
        poll(gate)?;
        let mut out = self
            .table
            .apply(&[(NearOperator::Single, C64::ONE)], q, options.workers);
        for (value, &f) in out.iter_mut().zip(&far) {
            *value = *value + f;
        }
        Ok(out)
    }
}

/// Burton-Miller under [`PanelRule::Triangle`]: far-field-only S, D and N
/// passes over the Radon cloud, near rows holding the dense rule's entries
/// minus the node sums those passes cover, and the dense rule's diagonals. D'
/// keeps its centroid adjoint pass and point closure.
struct BmCloudLayers {
    weights: Vec<f64>,
    single: FmmC64<'static>,
    double: FmmC64<'static>,
    mixed: FmmC64<'static>,
    adjoint: FmmC64<'static>,
    table: NearTable,
    diagonal_a: Vec<C64>,
    diagonal_b: Vec<C64>,
}

impl BmCloudLayers {
    /// Build the cloud passes, the near table, both diagonals and the
    /// leaf-block preconditioner. Returns the stored complex entries.
    #[allow(clippy::too_many_lines)] // one admission and assembly audit
    fn prepare(
        surface: &SpherePanels,
        k: f64,
        alpha: C64,
        max_pairs: usize,
        options: &IterativeOptions,
        gate: &CancelGate,
    ) -> Result<(Self, Vec<NearBlock>, usize), IterativeError> {
        let triangles = surface.triangles().ok_or(IterativeError::Invalid {
            what: "PanelRule::Triangle requires retained triangles",
        })?;
        if options.directional.is_some() {
            return Err(IterativeError::Invalid {
                what: "PanelRule::Triangle's far cloud is W5a only; directional rows use TriangleNear",
            });
        }
        let n = triangles.len();
        let cloud_capacity =
            options
                .leaf_capacity
                .checked_mul(RADON7.len())
                .ok_or(IterativeError::Invalid {
                    what: "leaf capacity overflow",
                })?;
        let (points, normals, weights, radius) = radon_cloud(surface, triangles);
        let far = |layer: usize, k: f64| -> Result<FmmC64<'static>, IterativeError> {
            let (p, nrm) = (points.clone(), normals.clone());
            let (order, budget) = (options.fmm_order, options.budget());
            let mut engine = match layer {
                0 => FmmC64::new_single_far_field_within(
                    HelmholtzSingle { k },
                    p,
                    order,
                    cloud_capacity,
                    budget,
                )?,
                1 => FmmC64::new_double_far_field_within(
                    HelmholtzDouble { k },
                    p,
                    nrm,
                    order,
                    cloud_capacity,
                    budget,
                )?,
                _ => FmmC64::new_hypersingular_far_field_within(
                    HelmholtzHypersingular { k },
                    p,
                    nrm,
                    order,
                    cloud_capacity,
                    budget,
                )?,
            };
            engine.set_workers(options.workers)?;
            Ok(engine)
        };
        poll(gate)?;
        let double = far(1, k)?;
        // Every far pair is then at least four radii from the target, where
        // the dense rule takes the same seven-point sum this pass makes.
        if double.leaf_box_size() < CBIE_MIN_LEAF_RADII * radius {
            return Err(IterativeError::Invalid {
                what: "PanelRule::Triangle cloud leaves must span CBIE_MIN_LEAF_RADII panel radii; raise leaf_capacity",
            });
        }
        let partition = CloudPartition::new(&double, n, points.len());
        let centroids = surface.centroids();
        let exact = (0..n)
            .map(|i| {
                partition
                    .panels(i)
                    .iter()
                    .filter(|&&j| {
                        j != i && triangle_is_near(centroids[i], centroids[j], triangles[j])
                    })
                    .count()
            })
            .sum::<usize>();
        if exact > max_pairs {
            return Err(HelmholtzError::NearRuleBudget {
                required_pairs: exact,
                maximum_pairs: max_pairs,
            }
            .into());
        }
        let row_len = |i: usize| partition.panels(i).iter().filter(|&&j| j != i).count();
        let operators = [
            NearOperator::Single,
            NearOperator::Double,
            NearOperator::Hypersingular,
        ];
        let (row_start, table_values) = NearTable::plan(n, operators.len(), row_len)?;
        let groups = match options.precond {
            Preconditioner::None => Vec::new(),
            Preconditioner::NearFieldBlockJacobi { max_block_panels } => {
                partition.leaf_blocks(n, max_block_panels)
            }
        };
        let near_entries = groups
            .iter()
            .try_fold(table_values, |count, block| {
                block
                    .len()
                    .checked_mul(block.len())
                    .and_then(|entries| count.checked_add(entries))
            })
            .ok_or(IterativeError::Invalid {
                what: "near-field extent overflow",
            })?;
        if near_entries > options.max_near_entries {
            return Err(IterativeError::NearFieldBudget {
                required_entries: near_entries,
                maximum_entries: options.max_near_entries,
            });
        }
        // Stored entry = dense rule entry minus the nodes the far pass sums,
        // so the row plus that pass is the dense entry to FMM accuracy.
        let covered = |i: usize, j: usize, k: f64| -> [C64; 3] {
            let (x, nx) = (centroids[i], surface.normals()[i]);
            let mut sums = [C64::ZERO; 3];
            for m in 0..RADON7.len() {
                let p = if m == 0 { j } else { n + 6 * j + m - 1 };
                if !partition.is_near(i, p) {
                    let (g, dgdny, _, d2g) = kernels(k, x, nx, points[p], normals[p]);
                    sums[0] = sums[0] + g.scale(weights[p]);
                    sums[1] = sums[1] + dgdny.scale(weights[p]);
                    sums[2] = sums[2] + d2g.scale(weights[p]);
                }
            }
            sums
        };
        poll(gate)?;
        let fill = |i: usize, cols: &mut [u32], values: &mut [C64]| {
            let mut slot = 0;
            for &j in partition.panels(i).iter().filter(|&&j| j != i) {
                let (layers, _) = triangle_rule_pair(surface, triangles, k, i, j);
                let sums = covered(i, j, k);
                cols[slot] = j as u32;
                for o in 0..3 {
                    values[3 * slot + o] = layers[o] - sums[o];
                }
                slot += 1;
            }
        };
        let table = NearTable::build(&operators, row_start, options.workers, gate, &fill)?;
        // Static rows under the same rule: the far static passes over unit
        // densities, plus the near rows' static entries less their covered
        // node sums.
        poll(gate)?;
        let unit: Vec<C64> = weights.iter().map(|&w| C64::from_re(w)).collect();
        let far_d0 = far(1, 0.0)?.far_potentials(&unit)?;
        poll(gate)?;
        let far_n0 = far(2, 0.0)?.far_potentials(&unit)?;
        poll(gate)?;
        let near_statics: Vec<(f64, f64)> = {
            let row = |i: usize| {
                let (mut d0, mut n0) = (0.0, 0.0);
                for &j in partition.panels(i).iter().filter(|&&j| j != i) {
                    let (_, statics) = triangle_rule_pair(surface, triangles, 0.0, i, j);
                    let sums = covered(i, j, 0.0);
                    d0 += statics[0] - sums[1].re;
                    n0 += statics[1] - sums[2].re;
                }
                (d0, n0)
            };
            let workers = options.workers.clamp(1, n.max(1));
            std::thread::scope(|scope| {
                let handles: Vec<_> = (0..workers)
                    .map(|w| {
                        let row = &row;
                        scope.spawn(move || {
                            (w * n / workers..(w + 1) * n / workers)
                                .map(row)
                                .collect::<Vec<_>>()
                        })
                    })
                    .collect();
                handles
                    .into_iter()
                    .flat_map(|h| h.join().expect("near static row worker panicked"))
                    .collect()
            })
        };
        // D' keeps its centroid rows: the point D_0 sum closes it.
        poll(gate)?;
        let areas: Vec<C64> = surface.areas().iter().map(|&a| C64::from_re(a)).collect();
        let point_d0 = options
            .normal_fmm(
                NormalLayer::Double,
                0.0,
                centroids.to_vec(),
                surface.normals().to_vec(),
            )?
            .potentials(&areas)?;
        let mut diagonal_a = Vec::with_capacity(n);
        let mut diagonal_b = Vec::with_capacity(n);
        for i in 0..n {
            let (s, n_reg) = triangle_self_terms(k, centroids[i], triangles[i]);
            let d0 = far_d0[i].re + near_statics[i].0;
            let n0 = far_n0[i].re + near_statics[i].1;
            diagonal_a.push(C64::from_re(-1.0 - d0) + alpha * (n_reg - C64::from_re(n0)));
            diagonal_b.push(s - alpha * point_d0[i]);
        }
        poll(gate)?;
        let mut blocks = Vec::with_capacity(groups.len());
        for indices in groups {
            poll(gate)?;
            let m = indices.len();
            let mut matrix = vec![C64::ZERO; m * m];
            for (a, &i) in indices.iter().enumerate() {
                for (b, &j) in indices.iter().enumerate() {
                    matrix[a * m + b] = if a == b {
                        diagonal_a[i]
                    } else {
                        // j's centroid shares i's leaf, so every node of j is
                        // near and the stored rows are the dense entry.
                        let entry = |op| {
                            table.get(i, j, op).ok_or(IterativeError::Invalid {
                                what: "near table lacks a same-leaf panel",
                            })
                        };
                        entry(NearOperator::Double)? + alpha * entry(NearOperator::Hypersingular)?
                    };
                }
            }
            let inverse =
                lu_complex(&matrix, m).map_err(|_| IterativeError::SingularNearBlock {
                    block: blocks.len(),
                })?;
            blocks.push(NearBlock { indices, inverse });
        }
        poll(gate)?;
        let single = far(0, k)?;
        poll(gate)?;
        let mixed = far(2, k)?;
        poll(gate)?;
        let adjoint = options.normal_fmm(
            NormalLayer::Adjoint,
            k,
            centroids.to_vec(),
            surface.normals().to_vec(),
        )?;
        Ok((
            Self {
                weights,
                single,
                double,
                mixed,
                adjoint,
                table,
                diagonal_a,
                diagonal_b,
            },
            blocks,
            near_entries,
        ))
    }

    /// Panel densities spread to their cloud nodes, with node weights.
    fn charges(&self, density: &[C64]) -> Vec<C64> {
        let n = density.len();
        (0..self.weights.len())
            .map(|p| density[if p < n { p } else { (p - n) / 6 }].scale(self.weights[p]))
            .collect()
    }

    /// The dense rule's A p.
    fn apply_a(&self, alpha: C64, p: &[C64], workers: usize) -> Result<Vec<C64>, IterativeError> {
        let charges = self.charges(p);
        let d = self.double.far_potentials(&charges)?;
        let h = self.mixed.far_potentials(&charges)?;
        let mut out = self.table.apply(
            &[
                (NearOperator::Double, C64::ONE),
                (NearOperator::Hypersingular, alpha),
            ],
            p,
            workers,
        );
        for (i, value) in out.iter_mut().enumerate() {
            *value = *value + d[i] + alpha * h[i] + self.diagonal_a[i] * p[i];
        }
        Ok(out)
    }

    /// The dense rule's B q, with D' on its centroid rows.
    fn apply_b(
        &self,
        alpha: C64,
        q: &[C64],
        areas: &[f64],
        workers: usize,
    ) -> Result<Vec<C64>, IterativeError> {
        let s = self.single.far_potentials(&self.charges(q))?;
        let weighted: Vec<C64> = q.iter().zip(areas).map(|(&q, &a)| q.scale(a)).collect();
        let adjoint = self.adjoint.potentials(&weighted)?;
        let mut out = self
            .table
            .apply(&[(NearOperator::Single, C64::ONE)], q, workers);
        for (i, value) in out.iter_mut().enumerate() {
            *value = *value + s[i] + alpha * adjoint[i] + self.diagonal_b[i] * q[i];
        }
        Ok(out)
    }
}

struct RealOp<'a, 's> {
    owner: &'a IterativeRadiation<'s>,
    gate: &'a CancelGate,
    first_error: RefCell<Option<IterativeError>>,
    matvecs: Cell<usize>,
    offset: usize,
    cancel: RefCell<&'a mut dyn FnMut(usize) -> bool>,
}
impl LinearOp for RealOp<'_, '_> {
    fn n(&self) -> usize {
        2 * self.owner.surface.areas().len()
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        if self.first_error.borrow().is_some() {
            y.fill(f64::NAN);
            return;
        }
        let result = if (self.cancel.borrow_mut())(self.offset + self.matvecs.get()) {
            Err(IterativeError::Cancelled { matvecs: 0 })
        } else {
            poll(self.gate)
        }
        .and_then(|()| {
            self.matvecs.set(self.matvecs.get() + 1);
            self.owner
                .apply(&self.owner.inverse(&unstack(x)), self.gate)
        });
        match result {
            Ok(v) => y.copy_from_slice(&stack(&v)),
            Err(e) => {
                *self.first_error.borrow_mut() = Some(e);
                y.fill(f64::NAN);
            }
        }
    }
    fn apply_transpose(&self, _x: &[f64], y: &mut [f64]) {
        *self.first_error.borrow_mut() = Some(IterativeError::Invalid {
            what: "transpose iterative radiation is not implemented",
        });
        y.fill(f64::NAN);
    }
}
struct ZeroOp(usize);
impl LinearOp for ZeroOp {
    fn n(&self) -> usize {
        self.0
    }
    fn apply(&self, _: &[f64], y: &mut [f64]) {
        y.fill(0.0);
    }
}
fn finite(v: C64) -> bool {
    v.re.is_finite() && v.im.is_finite()
}
fn stack(v: &[C64]) -> Vec<f64> {
    v.iter()
        .map(|z| z.re)
        .chain(v.iter().map(|z| z.im))
        .collect()
}
fn unstack(v: &[f64]) -> Vec<C64> {
    let n = v.len() / 2;
    (0..n).map(|i| C64::new(v[i], v[n + i])).collect()
}
fn host() -> String {
    let name = std::fs::read_to_string("/proc/sys/kernel/hostname").ok();
    let cpu = std::fs::read_to_string("/proc/cpuinfo")
        .ok()
        .and_then(|text| {
            text.lines()
                .find_map(|line| line.strip_prefix("model name").map(str::to_string))
        });
    format!(
        "{} {} {} {}",
        std::env::consts::OS,
        std::env::consts::ARCH,
        name.as_deref().unwrap_or("hostname-unavailable").trim(),
        cpu.as_deref().unwrap_or("cpu-unavailable").trim()
    )
}
fn peak_rss() -> Option<u64> {
    std::fs::read_to_string("/proc/self/status")
        .ok()?
        .lines()
        .find_map(|line| {
            let value = line
                .strip_prefix("VmHWM:")?
                .split_whitespace()
                .next()?
                .parse::<u64>()
                .ok()?;
            value.checked_mul(1024)
        })
}

/// One-shot convenience over the same checkpointable operator/state API.
pub fn solve_radiation_iterative(
    surface: &SpherePanels,
    k: f64,
    medium: Medium,
    velocity: &[C64],
    formulation: Formulation,
    options: &IterativeOptions,
    gate: &CancelGate,
) -> Result<(RadiationSolution, SolveReport, ResourceReceipt), IterativeError> {
    let operator = IterativeRadiation::prepare(surface, k, medium, formulation, options, gate)?;
    let mut state = operator.start(velocity, gate)?;
    operator.advance(&mut state, options.max_iters.saturating_add(1), gate)?;
    operator.finish(&state)
}

#[cfg(test)]
mod cbie_tests {
    use super::*;

    /// G1: every dense plain-CBIE entry is a stored near entry plus the exact
    /// point sum of its far nodes, or, for a fully separated panel, its
    /// seven-point sum within the CBIE_MIN_LEAF_RADII error table.
    #[test]
    fn cbie_near_table_and_far_nodes_rebuild_the_dense_owner_entries() {
        let surface = SpherePanels::icosphere(1.0, 3).unwrap();
        let triangles = surface.triangles().unwrap();
        let n = triangles.len();
        let k = 1.3;
        let options = IterativeOptions {
            leaf_capacity: 32,
            workers: 3,
            ..IterativeOptions::default()
        };
        let gate = CancelGate::new_clock_free();
        let (cbie, blocks, entries) = CbieLayers::prepare(&surface, k, &options, &gate).unwrap();
        assert!(entries > cbie.table.values.len() && !blocks.is_empty());
        let (amat, bmat) = super::super::assemble_dense(&surface, k, C64::ZERO);
        let neighbours = cbie.double.near_field_neighbor_leaves();
        let mut point_leaf = vec![0; cbie.points.len()];
        for (l, leaf) in cbie.double.near_field_blocks().enumerate() {
            for &p in leaf {
                point_leaf[p] = l;
            }
        }
        let (mut exact, mut corrected, mut separated) = (0, 0, 0);
        let mut worst = [0.0_f64; 3];
        for i in 0..n {
            let leaf = point_leaf[i];
            let (x, nx) = (surface.centroids()[i], surface.normals()[i]);
            for j in 0..n {
                let (mut single, mut double) = (C64::ZERO, C64::ZERO);
                let mut far_nodes = 0;
                for m in 0..RADON7.len() {
                    let p = if m == 0 { j } else { n + 6 * j + m - 1 };
                    if neighbours[leaf].binary_search(&point_leaf[p]).is_err() {
                        let (g, dgdny, _, _) =
                            kernels(k, x, nx, cbie.points[p], surface.normals()[j]);
                        single = single + g.scale(cbie.weights[p]);
                        double = double + dgdny.scale(cbie.weights[p]);
                        far_nodes += 1;
                    }
                }
                let dense_d = amat[i * n + j] + if i == j { C64::from_re(0.5) } else { C64::ZERO };
                let dense_s = bmat[i * n + j];
                let stored = (
                    cbie.table.get(i, j, NearOperator::Single),
                    cbie.table.get(i, j, NearOperator::Double),
                );
                assert_eq!(stored.0.is_some(), stored.1.is_some());
                if let (Some(s), Some(d)) = stored {
                    if far_nodes == 0 {
                        exact += 1;
                        assert_eq!(
                            (s.re.to_bits(), s.im.to_bits()),
                            (dense_s.re.to_bits(), dense_s.im.to_bits())
                        );
                        assert_eq!(
                            (d.re.to_bits(), d.im.to_bits()),
                            (dense_d.re.to_bits(), dense_d.im.to_bits())
                        );
                    } else {
                        corrected += 1;
                        worst[0] = worst[0].max((s + single - dense_s).abs() / dense_s.abs());
                        worst[0] = worst[0].max((d + double - dense_d).abs() / dense_s.abs());
                    }
                } else {
                    assert_eq!(far_nodes, RADON7.len());
                    separated += 1;
                    let radius = triangles[j]
                        .iter()
                        .map(|v| {
                            let c = surface.centroids()[j];
                            ((v[0] - c[0]).powi(2) + (v[1] - c[1]).powi(2) + (v[2] - c[2]).powi(2))
                                .sqrt()
                        })
                        .fold(0.0, f64::max);
                    worst[1] = worst[1].max((single - dense_s).abs() / dense_s.abs());
                    worst[2] = worst[2].max((double - dense_d).abs() * radius / dense_s.abs());
                }
            }
        }
        println!(
            "{{\"case\":\"cbie-near-far-split\",\"panels\":{n},\"exact\":{exact},\"corrected\":{corrected},\"separated\":{separated},\"corrected_rel\":{},\"separated_single_rel\":{},\"separated_double_rel\":{}}}",
            worst[0], worst[1], worst[2]
        );
        assert!(exact > 0 && corrected > 0 && separated > 0);
        assert!(worst[0] <= 1e-12, "{worst:?}");
        assert!(worst[1] <= 2e-6 && worst[2] <= 2e-6, "{worst:?}");
    }

    /// G0: leaves narrower than CBIE_MIN_LEAF_RADII panel radii, centroid-only
    /// surfaces, the directional route and Robin boundaries refuse.
    #[test]
    fn cbie_route_refusals() {
        let gate = CancelGate::new_clock_free();
        let surface = SpherePanels::icosphere(1.0, 3).unwrap();
        let deep = IterativeOptions {
            leaf_capacity: 2,
            ..IterativeOptions::default()
        };
        assert!(matches!(
            CbieLayers::prepare(&surface, 1.0, &deep, &gate),
            Err(IterativeError::Invalid { what }) if what.contains("CBIE_MIN_LEAF_RADII")
        ));
        let centroids_only = SpherePanels::new(
            surface.centroids().to_vec(),
            surface.normals().to_vec(),
            surface.areas().to_vec(),
        )
        .unwrap();
        assert!(matches!(
            CbieLayers::prepare(&centroids_only, 1.0, &IterativeOptions::default(), &gate),
            Err(IterativeError::Invalid { what }) if what.contains("retained triangles")
        ));
        let directional = IterativeOptions {
            directional: Some(fs_fmm::DirectionalOptions::for_tolerance(1e-3).unwrap()),
            ..IterativeOptions::default()
        };
        assert!(matches!(
            CbieLayers::prepare(&surface, 1.0, &directional, &gate),
            Err(IterativeError::Invalid { what }) if what.contains("W5a only")
        ));
        let budget = IterativeOptions {
            max_near_entries: 1,
            ..IterativeOptions::default()
        };
        assert!(matches!(
            CbieLayers::prepare(&surface, 1.0, &budget, &gate),
            Err(IterativeError::NearFieldBudget {
                maximum_entries: 1,
                ..
            })
        ));
        let beta = vec![C64::ZERO; surface.areas().len()];
        assert!(matches!(
            IterativeRadiation::prepare_bc(
                &surface,
                1.0,
                Medium::air(),
                Formulation::PlainCbie,
                &IterativeOptions::default(),
                super::super::bc::Admittance::Local {
                    beta: &beta,
                    declared_active: false,
                },
                &gate,
            ),
            Err(IterativeError::Invalid { what }) if what.contains("Robin")
        ));
    }

    /// A planted correction rule: rows at columns i and i + 1 (mod n).
    struct Planted(usize);
    impl NearRows for Planted {
        fn operators(&self) -> &[NearOperator] {
            &[
                NearOperator::Single,
                NearOperator::Double,
                NearOperator::Hypersingular,
            ]
        }
        fn row_len(&self, _: usize) -> usize {
            2
        }
        fn fill_row(&self, i: usize, cols: &mut [u32], values: &mut [C64]) {
            let mut columns = [i, (i + 1) % self.0];
            columns.sort_unstable();
            for (p, &j) in columns.iter().enumerate() {
                cols[p] = j as u32;
                for o in 0..3 {
                    values[3 * p + o] = C64::new(1e-3 * (o + 1) as f64, 1e-4 * (j % 7) as f64);
                }
            }
        }
        fn diagonal(&self, i: usize, alpha: C64) -> (C64, C64) {
            (alpha.scale(1e-3 * (i % 3) as f64), C64::new(2e-3, -1e-3))
        }
        fn stamp(&self) -> u64 {
            0x5EED
        }
    }

    /// G1/G0: a [`NearRows`] rule enters the Burton-Miller operator exactly
    /// as declared: A gains D + alpha N and dA_ii, B gains S and dB_ii.
    /// `None` is prepare; plain CBIE and an over-budget rule refuse.
    #[test]
    fn burton_miller_near_rows_enter_as_declared() {
        let surface = SpherePanels::icosphere(1.0, 2).unwrap();
        let n = surface.areas().len();
        let (k, medium) = (1.0, Medium::air());
        let gate = CancelGate::new_clock_free();
        let options = IterativeOptions {
            workers: 3,
            ..IterativeOptions::default()
        };
        let rule = Planted(n);
        let bm = Formulation::BurtonMiller;
        let prepare = |o: &IterativeOptions, near: Option<&dyn NearRows>| {
            IterativeRadiation::prepare_with_near(&surface, k, medium, bm, o, near, &gate)
        };
        let plain = IterativeRadiation::prepare(&surface, k, medium, bm, &options, &gate).unwrap();
        assert_eq!(prepare(&options, None).unwrap().stamp, plain.stamp);
        let ruled = prepare(&options, Some(&rule)).unwrap();
        assert_ne!(ruled.stamp, plain.stamp);
        let x: Vec<C64> = (0..n)
            .map(|i| C64::new(1.0 + 0.01 * i as f64, 0.5 - 0.003 * i as f64))
            .collect();
        let q: Vec<C64> = x
            .iter()
            .map(|&v| v * C64::new(0.0, k * medium.sound_speed * medium.density))
            .collect();
        let alpha = plain.alpha;
        let pairs = [
            (
                plain.apply(&x, &gate).unwrap(),
                ruled.apply(&x, &gate).unwrap(),
            ),
            (
                plain.boundary_rhs(&x, &gate).unwrap(),
                ruled.boundary_rhs(&x, &gate).unwrap(),
            ),
        ];
        for (side, (base, with)) in pairs.iter().enumerate() {
            let scale = base.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
            let mut worst = 0.0_f64;
            for i in 0..n {
                let (mut cols, mut values) = ([0u32; 2], [C64::ZERO; 6]);
                rule.fill_row(i, &mut cols, &mut values);
                let (da, db) = rule.diagonal(i, alpha);
                let mut extra = if side == 0 { da * x[i] } else { db * q[i] };
                for p in 0..2 {
                    let j = cols[p] as usize;
                    extra = extra
                        + if side == 0 {
                            (values[3 * p + 1] + alpha * values[3 * p + 2]) * x[j]
                        } else {
                            values[3 * p] * q[j]
                        };
                }
                worst = worst.max((with[i] - base[i] - extra).abs() / scale);
            }
            assert!(worst <= 1e-14, "side {side}: {worst}");
        }
        // The corrected operator solves through its own preconditioner.
        let mut state = ruled.start(&x, &gate).unwrap();
        ruled.advance(&mut state, 100, &gate).unwrap();
        assert!(ruled.finish(&state).is_ok());
        // Its values join the near-entry budget; the cap is inclusive.
        let exact = IterativeOptions {
            max_near_entries: plain.near_entries + 6 * n,
            ..options.clone()
        };
        assert!(prepare(&exact, Some(&rule)).is_ok());
        let short = IterativeOptions {
            max_near_entries: plain.near_entries + 6 * n - 1,
            ..options.clone()
        };
        assert!(matches!(
            prepare(&short, Some(&rule)),
            Err(IterativeError::NearFieldBudget { .. })
        ));
        assert!(matches!(
            IterativeRadiation::prepare_with_near(
                &surface,
                k,
                medium,
                Formulation::PlainCbie,
                &options,
                Some(&rule),
                &gate
            ),
            Err(IterativeError::Invalid { what }) if what.contains("already stores")
        ));
    }
}
