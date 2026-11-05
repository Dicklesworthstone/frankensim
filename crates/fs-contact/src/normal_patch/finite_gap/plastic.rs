//! Hardness-capped pressure saturation for finite-gap half-space contact.
//!
//! The pressure-saturation approximation solves the convex box problem
//! `min 0.5 pᵀ C p - (approach - gap - old_plastic_depth - recession)ᵀ p`, with
//! `0 <= p[i] <= hardness[i]`. `C` is the same exact rectangular Boussinesq
//! compliance used by the elastic finite-gap model. Unaccommodated closure at
//! a saturated cell becomes irreversible residual indentation. The original
//! undeformed gap datum is fixed at construction, including on unload/reload.
//!
//! This is the penetration-hardness approximation, not subsurface J2
//! plasticity, hardening, pile-up, a hardness-to-yield-stress conversion, or a
//! thermal model. Plastic work is accounted separately from elastic storage;
//! its conversion to heat is a caller-owned material question. See the model
//! description in <https://doi.org/10.1038/s41467-018-02981-y> and the contrast
//! with full plasticity in <https://doi.org/10.1007/s11249-023-01819-z>.

mod continuation;
pub use continuation::*;
mod load;
pub use load::*;

use core::fmt;

use fs_blake3::{ContentHash, hash_domain};
use fs_exec::Cx;
use fs_la::factor::{FactorError, FactorRunError, cholesky_cancellable};

use super::{
    ComplianceKernel, FiniteGapContactError, FiniteGapContactReceipt, FiniteGapContactRequest,
    FiniteGapGrid, build_receipt, verify_boundary_clearance,
};

/// Identity domain for the pressure-saturation constitutive model.
pub const FINITE_GAP_PLASTIC_MODEL_ID: &str =
    "org.frankensim.fs-contact.finite-gap-pressure-saturation.v1";

/// Maximum simultaneously free cells in the dense reference factorization.
///
/// Bounds its matrix plus factor to about 64 MiB. The full sampled gap can
/// contain up to [`super::MAX_FINITE_GAP_CELLS`]; lower/upper-bound cells do not
/// consume dense factor storage. This is a numerical work limit, not a
/// material parameter or an asserted physical applicability threshold.
pub const MAX_FINITE_GAP_PLASTIC_FREE_CELLS: usize = 2_048;

/// Immutable physical inputs and numerical controls for pressure saturation.
#[derive(Debug, Clone, PartialEq)]
pub struct FiniteGapPlasticConfig {
    /// Cell-centred tangent-plane grid.
    pub grid: FiniteGapGrid,
    /// Original undeformed separation at the grid centres [m].
    ///
    /// Its minimum is subtracted exactly once when the model is compiled.
    pub undeformed_gap_m: Vec<f64>,
    /// Classical combined half-space modulus `E*` [Pa].
    pub reduced_modulus_pa: f64,
    /// Positive indentation pressure cap for each cell [Pa].
    ///
    /// These are physical indentation hardness pressures at the applicable
    /// material state, not raw Rockwell/Mohs values or inferred yield stresses.
    pub hardness_pa: Vec<f64>,
    /// Maximum feasible active-set iterations per attempted advance.
    pub maximum_active_set_iterations: usize,
    /// Absolute box-complementarity tolerance in displacement units [m].
    pub complementarity_tolerance_m: f64,
    /// Required inactive ring between any positive pressure and each edge.
    pub boundary_clearance_cells: usize,
}

/// Compiled immutable half-space kernel and fixed hardness field.
///
/// Construction precomputes the translation-invariant compliance once. Each
/// advance consumes accepted history immutably and returns a new candidate;
/// a refused solve cannot publish partial indentation or energy accounting.
#[derive(Clone, PartialEq)]
pub struct FiniteGapPlasticModel {
    config: FiniteGapPlasticConfig,
    identity: ContentHash,
    gap_datum_m: f64,
    gap_m: Vec<f64>,
    compliance: ComplianceKernel,
}

impl fmt::Debug for FiniteGapPlasticModel {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("FiniteGapPlasticModel")
            .field("config", &self.config)
            .field("identity", &self.identity)
            .field("gap_datum_m", &self.gap_datum_m)
            .finish_non_exhaustive()
    }
}

/// Opaque, model-bound accepted or candidate constitutive history.
#[derive(Debug, Clone, PartialEq)]
pub struct FiniteGapPlasticState {
    owner_identity: ContentHash,
    identity: ContentHash,
    approach_m: f64,
    plastic_depth_m: Vec<f64>,
    recession_m: Vec<f64>,
    pressure_pa: Vec<f64>,
    cumulative_plastic_work_j: f64,
}

impl FiniteGapPlasticState {
    /// Identity of the exact compiled physical inputs and numerical controls.
    pub fn owner_identity(&self) -> ContentHash {
        self.owner_identity
    }

    /// Content identity of this complete constitutive history.
    pub fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Last accepted rigid approach from the original gap minimum [m].
    pub fn approach_m(&self) -> f64 {
        self.approach_m
    }

    /// Cumulative residual indentation at each original grid location [m].
    pub fn plastic_depth_m(&self) -> &[f64] {
        &self.plastic_depth_m
    }

    /// Accumulated relative surface recession, separate from plastic indentation [m].
    ///
    /// Recession increases the original relative gap without moving its datum.
    pub fn recession_m(&self) -> &[f64] {
        &self.recession_m
    }

    /// Last equilibrated compressive pressure field [Pa].
    pub fn pressure_pa(&self) -> &[f64] {
        &self.pressure_pa
    }

    /// Accumulated nonnegative plastic work along accepted advances [J].
    pub fn cumulative_plastic_work_j(&self) -> f64 {
        self.cumulative_plastic_work_j
    }
}

/// Candidate history and the response that produced it.
#[derive(Debug, Clone, PartialEq)]
pub struct FiniteGapPlasticAdvance {
    /// Replace the caller's accepted state with this only when accepting the step.
    pub candidate: FiniteGapPlasticState,
    /// Independently inspectable pressure, storage, and plastic-work response.
    pub receipt: FiniteGapPlasticReceipt,
}

/// Constitutive response for one converged pressure-saturation advance.
#[derive(Debug, Clone, PartialEq)]
pub struct FiniteGapPlasticReceipt {
    /// Identity binding the owner, predecessor, candidate, and response.
    pub identity: ContentHash,
    /// Exact prior state used to compute this response.
    pub predecessor_identity: ContentHash,
    /// Common pressure, force, moments, and elastic storage integrals.
    ///
    /// `reversible_energy_j` is the current half-space energy, not total
    /// indentation work. `complementarity_residual_m` includes the new plastic
    /// depth. Its identity is the pressure-saturation response identity.
    pub normal: FiniteGapContactReceipt,
    /// Nonnegative new residual indentation at each cell [m].
    pub plastic_increment_m: Vec<f64>,
    /// Incremental plastic work `sum(hardness * plastic_increment * cell_area)` [J].
    ///
    /// No fraction is automatically assigned to heat. Endpoint work and stored
    /// energy do not substitute for integrating external work along a load path.
    pub plastic_work_j: f64,
    /// Cumulative plastic work of the candidate history [J].
    pub cumulative_plastic_work_j: f64,
    /// Area of the positive-pressure constant-pressure cells [m²].
    pub contact_area_m2: f64,
    /// Number of cells acquiring strictly positive residual indentation.
    pub yielded_cells: usize,
    /// Maximum violation of the box-QP KKT conditions [m].
    pub kkt_residual_m: f64,
}

/// Refusal without mutation of the caller's accepted plastic history.
#[derive(Debug, Clone, PartialEq)]
pub enum FiniteGapPlasticError {
    /// Invalid physical/numerical input, cancellation, boundary, or convergence failure.
    Contact(FiniteGapContactError),
    /// A state belongs to different physical inputs or numerical controls.
    ForeignState,
    /// Dense free-set factorization would exceed its published work/storage bound.
    FreeSetTooLarge {
        /// Simultaneously unbound cells in the attempted factorization.
        free_cells: usize,
        /// Maximum supported dense dimension.
        maximum: usize,
    },
}

impl From<FiniteGapContactError> for FiniteGapPlasticError {
    fn from(error: FiniteGapContactError) -> Self {
        Self::Contact(error)
    }
}

impl fmt::Display for FiniteGapPlasticError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Contact(error) => write!(f, "{error}"),
            Self::ForeignState => {
                f.write_str("finite-gap plastic history belongs to another model")
            }
            Self::FreeSetTooLarge {
                free_cells,
                maximum,
            } => write!(
                f,
                "finite-gap plastic free set has {free_cells} cells; dense work limit is {maximum}"
            ),
        }
    }
}

impl core::error::Error for FiniteGapPlasticError {}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Bound {
    Lower,
    Free,
    Upper,
}

struct PlasticStepMetrics {
    iterations: usize,
    physical_residual_m: f64,
    kkt_residual_m: f64,
    plastic_increment_m: Vec<f64>,
    plastic_work_j: f64,
    yielded_cells: usize,
}

impl FiniteGapPlasticModel {
    /// Validate physical inputs and compile the shared Boussinesq compliance.
    pub fn new(config: FiniteGapPlasticConfig, cx: &Cx<'_>) -> Result<Self, FiniteGapPlasticError> {
        checkpoint(cx)?;
        let count = config.grid.cell_count()?;
        if config.undeformed_gap_m.len() != count
            || config.undeformed_gap_m.iter().any(|g| !g.is_finite())
        {
            return Err(invalid("undeformed gap shape or values"));
        }
        if config.hardness_pa.len() != count
            || config
                .hardness_pa
                .iter()
                .any(|h| !h.is_finite() || *h <= 0.0)
        {
            return Err(invalid("indentation hardness field"));
        }
        if !config.reduced_modulus_pa.is_finite() || config.reduced_modulus_pa <= 0.0 {
            return Err(invalid("reduced modulus"));
        }
        if config.maximum_active_set_iterations == 0
            || !config.complementarity_tolerance_m.is_finite()
            || config.complementarity_tolerance_m <= 0.0
        {
            return Err(invalid("active-set controls"));
        }
        let grid = config.grid;
        if config.boundary_clearance_cells == 0
            || config.boundary_clearance_cells >= grid.cells_x.div_ceil(2)
            || config.boundary_clearance_cells >= grid.cells_y.div_ceil(2)
        {
            return Err(invalid("boundary clearance"));
        }
        let area = grid.cell_width_m * grid.cell_depth_m;
        if !area.is_finite() || area <= 0.0 || !(area * count as f64).is_finite() {
            return Err(invalid("finite cell area"));
        }
        let gap_datum_m = config
            .undeformed_gap_m
            .iter()
            .copied()
            .fold(f64::INFINITY, f64::min);
        let gap_m: Vec<_> = config
            .undeformed_gap_m
            .iter()
            .map(|g| g - gap_datum_m)
            .collect();
        if gap_m.iter().any(|g| !g.is_finite()) {
            return Err(invalid("relative gap range"));
        }
        let compliance = ComplianceKernel::build(grid, config.reduced_modulus_pa);
        if compliance
            .values
            .iter()
            .any(|c| !c.is_finite() || *c <= 0.0)
        {
            return Err(invalid("finite positive half-space compliance"));
        }
        checkpoint(cx)?;
        let identity = config_identity(&config);
        Ok(Self {
            config,
            identity,
            gap_datum_m,
            gap_m,
            compliance,
        })
    }

    /// Identity of the immutable physical model and numerical controls.
    pub fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Validated physical inputs and numerical controls.
    pub fn config(&self) -> &FiniteGapPlasticConfig {
        &self.config
    }

    /// Original constant height datum subtracted once during construction [m].
    pub fn undeformed_gap_datum_m(&self) -> f64 {
        self.gap_datum_m
    }

    /// Virgin unloaded history at zero approach from the original gap minimum.
    pub fn initial_state(&self) -> FiniteGapPlasticState {
        self.make_state(
            0.0,
            vec![0.0; self.gap_m.len()],
            vec![0.0; self.gap_m.len()],
            vec![0.0; self.gap_m.len()],
            0.0,
        )
    }

    /// Compute equilibrium and irreversible indentation without changing `state`.
    ///
    /// Signed finite approaches support full separation. The fixed gap datum
    /// and residual indentation survive opening and unloading. No contact patch
    /// reaching the requested boundary clearance is accepted.
    pub fn advance(
        &self,
        state: &FiniteGapPlasticState,
        approach_m: f64,
        cx: &Cx<'_>,
    ) -> Result<FiniteGapPlasticAdvance, FiniteGapPlasticError> {
        checkpoint(cx)?;
        self.check_state_owner(state)?;
        self.advance_with_history(
            state,
            approach_m,
            &state.recession_m,
            &state.pressure_pa,
            cx,
        )
    }

    fn check_state_owner(
        &self,
        state: &FiniteGapPlasticState,
    ) -> Result<(), FiniteGapPlasticError> {
        if state.owner_identity != self.identity {
            return Err(FiniteGapPlasticError::ForeignState);
        }
        Ok(())
    }

    // Only checked public entry points may transfer accepted history here.
    // In particular, material refresh must verify the predecessor owner and
    // exact original geometry/controls before rebinding to the successor.
    fn advance_with_history(
        &self,
        state: &FiniteGapPlasticState,
        approach_m: f64,
        recession_m: &[f64],
        pressure_guess: &[f64],
        cx: &Cx<'_>,
    ) -> Result<FiniteGapPlasticAdvance, FiniteGapPlasticError> {
        checkpoint(cx)?;
        if !approach_m.is_finite() {
            return Err(invalid("finite approach"));
        }
        let target: Vec<f64> = self
            .gap_m
            .iter()
            .zip(&state.plastic_depth_m)
            .zip(recession_m)
            .map(|((g, z), r)| approach_m - g - z - r)
            .collect();
        if target.iter().any(|b| !b.is_finite()) {
            return Err(invalid("finite trial closure"));
        }
        let (pressure, opening, iterations, kkt_residual_m) =
            self.solve_box(&target, pressure_guess, cx)?;
        let active: Vec<usize> = pressure
            .iter()
            .enumerate()
            .filter_map(|(i, p)| (*p > 0.0).then_some(i))
            .collect();
        verify_boundary_clearance(
            self.config.grid,
            &active,
            self.config.boundary_clearance_cells,
        )?;
        let area = self.config.grid.cell_width_m * self.config.grid.cell_depth_m;
        let mut plastic_depth = state.plastic_depth_m.clone();
        let mut plastic_increment_m = vec![0.0; pressure.len()];
        let mut plastic_work_per_area = 0.0;
        let mut physical_residual = 0.0_f64;
        let mut yielded_cells = 0;
        for i in 0..pressure.len() {
            let closure = target[i] - opening[i];
            if pressure[i] == self.config.hardness_pa[i] && closure > 0.0 {
                // Use the representable change in the cumulative history for
                // the work account; never record unrepresentable increments.
                plastic_depth[i] += closure;
                plastic_increment_m[i] = plastic_depth[i] - state.plastic_depth_m[i];
                plastic_work_per_area = self.config.hardness_pa[i]
                    .mul_add(plastic_increment_m[i], plastic_work_per_area);
                yielded_cells += usize::from(plastic_increment_m[i] > 0.0);
            }
            let separation = opening[i] - target[i] + plastic_increment_m[i];
            physical_residual = physical_residual.max(if pressure[i] > 0.0 {
                separation.abs()
            } else {
                (-separation).max(0.0)
            });
        }
        let plastic_work_j = plastic_work_per_area * area;
        let cumulative_plastic_work_j = state.cumulative_plastic_work_j + plastic_work_j;
        if !plastic_work_j.is_finite()
            || !cumulative_plastic_work_j.is_finite()
            || plastic_depth.iter().any(|z| !z.is_finite())
        {
            return Err(invalid("finite accumulated plastic history"));
        }
        if physical_residual > self.config.complementarity_tolerance_m {
            return Err(FiniteGapContactError::ActiveSetDidNotConverge {
                iterations,
                residual_m: physical_residual,
            }
            .into());
        }
        let candidate = self.make_state(
            approach_m,
            plastic_depth,
            recession_m.to_vec(),
            pressure,
            cumulative_plastic_work_j,
        );
        let metrics = PlasticStepMetrics {
            iterations,
            physical_residual_m: physical_residual,
            kkt_residual_m,
            plastic_increment_m,
            plastic_work_j,
            yielded_cells,
        };
        let receipt =
            self.receipt_for_candidate(state.identity, &candidate, opening, metrics, cx)?;
        Ok(FiniteGapPlasticAdvance { candidate, receipt })
    }

    fn receipt_for_candidate(
        &self,
        predecessor_identity: ContentHash,
        candidate: &FiniteGapPlasticState,
        opening: Vec<f64>,
        metrics: PlasticStepMetrics,
        cx: &Cx<'_>,
    ) -> Result<FiniteGapPlasticReceipt, FiniteGapPlasticError> {
        let request = FiniteGapContactRequest {
            grid: self.config.grid,
            undeformed_gap_m: self.config.undeformed_gap_m.clone(),
            approach_m: candidate.approach_m,
            reduced_modulus_pa: self.config.reduced_modulus_pa,
            maximum_active_set_iterations: self.config.maximum_active_set_iterations,
            complementarity_tolerance_m: self.config.complementarity_tolerance_m,
            boundary_clearance_cells: self.config.boundary_clearance_cells,
        };
        let mut normal = build_receipt(
            &request,
            self.gap_datum_m,
            candidate.pressure_pa.clone(),
            opening,
            metrics.iterations,
            metrics.physical_residual_m,
        )?;
        if !normal.normal_force_n.is_finite()
            || !normal.reversible_energy_j.is_finite()
            || normal.pressure_centroid_m.iter().any(|x| !x.is_finite())
            || normal
                .equivalent_pressure_semiaxes_m
                .iter()
                .any(|x| !x.is_finite())
        {
            return Err(invalid("finite pressure integrals"));
        }
        let mut bytes = Vec::with_capacity(144);
        bytes.extend_from_slice(&self.identity.0);
        bytes.extend_from_slice(&predecessor_identity.0);
        bytes.extend_from_slice(&candidate.identity.0);
        bytes.extend_from_slice(&normal.identity.0);
        bytes.extend_from_slice(&metrics.plastic_work_j.to_bits().to_le_bytes());
        bytes.extend_from_slice(&metrics.kkt_residual_m.to_bits().to_le_bytes());
        let identity = hash_domain(
            "org.frankensim.fs-contact.finite-gap-plastic-response.v2",
            &bytes,
        );
        normal.identity = identity;
        let area = self.config.grid.cell_width_m * self.config.grid.cell_depth_m;
        let contact_area_m2 = normal.active_cells as f64 * area;
        checkpoint(cx)?;
        Ok(FiniteGapPlasticReceipt {
            identity,
            predecessor_identity,
            normal,
            plastic_increment_m: metrics.plastic_increment_m,
            plastic_work_j: metrics.plastic_work_j,
            cumulative_plastic_work_j: candidate.cumulative_plastic_work_j,
            contact_area_m2,
            yielded_cells: metrics.yielded_cells,
            kkt_residual_m: metrics.kkt_residual_m,
        })
    }

    fn make_state(
        &self,
        approach_m: f64,
        plastic_depth_m: Vec<f64>,
        recession_m: Vec<f64>,
        pressure_pa: Vec<f64>,
        cumulative_plastic_work_j: f64,
    ) -> FiniteGapPlasticState {
        let mut bytes = Vec::with_capacity(48 + 24 * pressure_pa.len());
        bytes.extend_from_slice(&self.identity.0);
        bytes.extend_from_slice(&approach_m.to_bits().to_le_bytes());
        bytes.extend_from_slice(&cumulative_plastic_work_j.to_bits().to_le_bytes());
        for value in plastic_depth_m
            .iter()
            .chain(&recession_m)
            .chain(&pressure_pa)
        {
            bytes.extend_from_slice(&value.to_bits().to_le_bytes());
        }
        let identity = hash_domain(
            "org.frankensim.fs-contact.finite-gap-plastic-state.v2",
            &bytes,
        );
        FiniteGapPlasticState {
            owner_identity: self.identity,
            identity,
            approach_m,
            plastic_depth_m,
            recession_m,
            pressure_pa,
            cumulative_plastic_work_j,
        }
    }

    // Primal feasible active set: solve the free equality subproblem, advance
    // only to the first bound hit, then release the single worst KKT-violating
    // bound once the free minimum is feasible. Thus no infeasible pressure
    // iterates, ad-hoc clipping solve, or simultaneous switching cycle is used.
    fn solve_box(
        &self,
        target: &[f64],
        previous_pressure: &[f64],
        cx: &Cx<'_>,
    ) -> Result<(Vec<f64>, Vec<f64>, usize, f64), FiniteGapPlasticError> {
        let mut pressure = previous_pressure.to_vec();
        let mut opening = self.apply_compliance(&pressure, cx)?;
        let tolerance = self.config.complementarity_tolerance_m;
        let mut bound: Vec<Bound> = pressure
            .iter()
            .enumerate()
            .map(|(i, p)| {
                if *p == 0.0 {
                    Bound::Lower
                } else if *p == self.config.hardness_pa[i] {
                    Bound::Upper
                } else {
                    Bound::Free
                }
            })
            .collect();
        // Seed the first free set with violated bounds. Later releases are
        // strictly one at a time. This makes virgin elastic contact practical
        // while retaining the feasible active-set convergence mechanism.
        for i in 0..pressure.len() {
            let gradient = opening[i] - target[i];
            if (bound[i] == Bound::Lower && gradient < -tolerance)
                || (bound[i] == Bound::Upper && gradient > tolerance)
            {
                bound[i] = Bound::Free;
            }
        }
        let initial_residual = box_residual(&pressure, &opening, target, &self.config.hardness_pa);
        if initial_residual <= tolerance {
            return Ok((pressure, opening, 0, initial_residual));
        }
        for iteration in 1..=self.config.maximum_active_set_iterations {
            checkpoint(cx)?;
            let free: Vec<usize> = bound
                .iter()
                .enumerate()
                .filter_map(|(i, b)| (*b == Bound::Free).then_some(i))
                .collect();
            if free.len() > MAX_FINITE_GAP_PLASTIC_FREE_CELLS {
                return Err(FiniteGapPlasticError::FreeSetTooLarge {
                    free_cells: free.len(),
                    maximum: MAX_FINITE_GAP_PLASTIC_FREE_CELLS,
                });
            }
            if !free.is_empty() {
                let n = free.len();
                let upper: Vec<usize> = bound
                    .iter()
                    .enumerate()
                    .filter_map(|(i, b)| (*b == Bound::Upper).then_some(i))
                    .collect();
                let mut matrix = vec![0.0; n * n];
                let mut solution = Vec::with_capacity(n);
                for (row, &i) in free.iter().enumerate() {
                    checkpoint(cx)?;
                    for (col, &j) in free.iter().enumerate().take(row + 1) {
                        let c = self.compliance.compliance(i, j);
                        matrix[row * n + col] = c;
                        matrix[col * n + row] = c;
                    }
                    let mut rhs = target[i];
                    for &j in &upper {
                        rhs = (-self.compliance.compliance(i, j)).mul_add(pressure[j], rhs);
                    }
                    solution.push(rhs);
                }
                if solution.iter().any(|p| !p.is_finite()) {
                    return Err(invalid("finite constrained pressure right-hand side"));
                }
                let factor = cholesky_cancellable(&matrix, n, || cx.checkpoint().is_err())
                    .map_err(factor_error)?;
                factor
                    .solve_cancellable(&mut solution, || cx.checkpoint().is_err())
                    .map_err(factor_error)?;
                if solution.iter().any(|p| !p.is_finite()) {
                    return Err(invalid("finite pressure solution"));
                }
                let mut fraction = 1.0_f64;
                let mut blocker = None;
                for (slot, &i) in free.iter().enumerate() {
                    let proposed = solution[slot];
                    let hit = if proposed < 0.0 {
                        Some((pressure[i] / (pressure[i] - proposed), Bound::Lower))
                    } else if proposed > self.config.hardness_pa[i] {
                        Some((
                            (self.config.hardness_pa[i] - pressure[i]) / (proposed - pressure[i]),
                            Bound::Upper,
                        ))
                    } else {
                        None
                    };
                    if let Some((alpha, side)) = hit {
                        if alpha < fraction || (alpha == fraction && blocker.is_none()) {
                            fraction = alpha;
                            blocker = Some((i, side));
                        }
                    }
                }
                if !fraction.is_finite() || fraction < 0.0 {
                    return Err(invalid("finite feasible pressure step"));
                }
                for (slot, &i) in free.iter().enumerate() {
                    pressure[i] = if blocker.is_none() {
                        solution[slot]
                    } else {
                        fraction
                            .mul_add(solution[slot] - pressure[i], pressure[i])
                            .clamp(0.0, self.config.hardness_pa[i])
                    };
                }
                if pressure
                    .iter()
                    .enumerate()
                    .any(|(i, p)| !p.is_finite() || *p < 0.0 || *p > self.config.hardness_pa[i])
                {
                    return Err(invalid("finite feasible pressure iterate"));
                }
                if let Some((i, side)) = blocker {
                    bound[i] = side;
                    pressure[i] = if side == Bound::Lower {
                        0.0
                    } else {
                        self.config.hardness_pa[i]
                    };
                    // Degenerate zero-length steps can put several starting
                    // free variables immediately on the same blocking bound.
                    // Fix those exact ties too; only bound releases are singular.
                    if fraction == 0.0 {
                        for (slot, &j) in free.iter().enumerate() {
                            if pressure[j] == 0.0 && solution[slot] < 0.0 {
                                bound[j] = Bound::Lower;
                            } else if pressure[j] == self.config.hardness_pa[j]
                                && solution[slot] > self.config.hardness_pa[j]
                            {
                                bound[j] = Bound::Upper;
                            }
                        }
                    }
                    continue;
                }
            }
            opening = self.apply_compliance(&pressure, cx)?;
            let residual = box_residual(&pressure, &opening, target, &self.config.hardness_pa);
            if residual <= tolerance {
                return Ok((pressure, opening, iteration, residual));
            }
            let mut release = None;
            let mut worst = tolerance;
            for (i, side) in bound.iter().enumerate() {
                let gradient = opening[i] - target[i];
                let violation = match side {
                    Bound::Lower => -gradient,
                    Bound::Upper => gradient,
                    Bound::Free => 0.0,
                };
                if violation > worst {
                    release = Some(i);
                    worst = violation;
                }
            }
            if let Some(i) = release {
                bound[i] = Bound::Free;
            } else {
                // Factor roundoff already exceeds the requested tolerance;
                // repeating exactly the same free solve cannot improve it.
                return Err(FiniteGapContactError::ActiveSetDidNotConverge {
                    iterations: iteration,
                    residual_m: residual,
                }
                .into());
            }
        }
        // Last feasible line search may have changed pressure without a full
        // gradient evaluation. Report the residual of that actual iterate.
        opening = self.apply_compliance(&pressure, cx)?;
        let residual = box_residual(&pressure, &opening, target, &self.config.hardness_pa);
        Err(FiniteGapContactError::ActiveSetDidNotConverge {
            iterations: self.config.maximum_active_set_iterations,
            residual_m: residual,
        }
        .into())
    }

    fn apply_compliance(
        &self,
        pressure: &[f64],
        cx: &Cx<'_>,
    ) -> Result<Vec<f64>, FiniteGapPlasticError> {
        let active: Vec<usize> = pressure
            .iter()
            .enumerate()
            .filter_map(|(i, p)| (*p > 0.0).then_some(i))
            .collect();
        let mut opening = vec![0.0; pressure.len()];
        for (i, displacement) in opening.iter_mut().enumerate() {
            checkpoint(cx)?;
            for &j in &active {
                *displacement = self
                    .compliance
                    .compliance(i, j)
                    .mul_add(pressure[j], *displacement);
            }
            if !displacement.is_finite() {
                return Err(invalid("finite elastic opening"));
            }
        }
        Ok(opening)
    }
}

fn box_residual(pressure: &[f64], opening: &[f64], target: &[f64], hardness: &[f64]) -> f64 {
    let mut residual = 0.0_f64;
    for i in 0..pressure.len() {
        let gradient = opening[i] - target[i];
        residual = residual.max(if pressure[i] == 0.0 {
            (-gradient).max(0.0)
        } else if pressure[i] == hardness[i] {
            gradient.max(0.0)
        } else {
            gradient.abs()
        });
    }
    residual
}

fn config_identity(config: &FiniteGapPlasticConfig) -> ContentHash {
    let mut bytes = Vec::with_capacity(80 + 16 * config.undeformed_gap_m.len());
    bytes.extend_from_slice(&(config.grid.cells_x as u64).to_le_bytes());
    bytes.extend_from_slice(&(config.grid.cells_y as u64).to_le_bytes());
    for value in [
        config.grid.cell_width_m,
        config.grid.cell_depth_m,
        config.reduced_modulus_pa,
        config.complementarity_tolerance_m,
    ] {
        bytes.extend_from_slice(&value.to_bits().to_le_bytes());
    }
    bytes.extend_from_slice(&(config.maximum_active_set_iterations as u64).to_le_bytes());
    bytes.extend_from_slice(&(config.boundary_clearance_cells as u64).to_le_bytes());
    for value in config.undeformed_gap_m.iter().chain(&config.hardness_pa) {
        bytes.extend_from_slice(&value.to_bits().to_le_bytes());
    }
    hash_domain(FINITE_GAP_PLASTIC_MODEL_ID, &bytes)
}

fn factor_error(error: FactorRunError) -> FiniteGapPlasticError {
    match error {
        FactorRunError::Cancelled => FiniteGapContactError::Cancelled.into(),
        FactorRunError::Factor(FactorError::NotSpd { index } | FactorError::Singular { index }) => {
            FiniteGapContactError::ComplianceNotPositiveDefinite {
                active_index: index,
            }
            .into()
        }
    }
}

fn invalid(field: &'static str) -> FiniteGapPlasticError {
    FiniteGapContactError::InvalidInput { field }.into()
}

fn checkpoint(cx: &Cx<'_>) -> Result<(), FiniteGapPlasticError> {
    cx.checkpoint()
        .map_err(|_| FiniteGapContactError::Cancelled.into())
}
