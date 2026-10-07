//! Component means use a region-supported dual on the COMPLETE thermal domain.
//! Restricting the primal domain or slicing a whole-domain mean certificate
//! would change the functional and lose heat transfer through other components.
use super::*;

/// An explicit, nonempty union of complete tetrahedra in the caller's mesh.
/// This is selection admission, not a geometry or material certificate. The
/// caller must preserve the original element numbering when reusing it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RegionMeanSelection {
    cell_count: usize,
    cells: Vec<usize>,
}

impl RegionMeanSelection {
    /// Admit the full cell list before geometry, allocation of dual fields or
    /// physical solves. Repeated cells are errors, not additional weights.
    /// Selection order is canonicalized; disconnected selected regions are valid.
    pub fn new(
        cell_count: usize,
        cells: &[usize],
        budget: FluxBudget,
        mut keep_going: impl FnMut() -> bool,
    ) -> Result<Self, TetError> {
        poll(&mut keep_going)?;
        if cell_count > budget.max_cells || cells.len() > budget.max_cells
            || budget.max_iterations > 1_000_000 {
            return Err(TetError::Budget);
        }
        if cell_count == 0 || cells.is_empty() || cells.len() > cell_count {
            return Err(TetError::Invalid("region must be a nonempty union of distinct domain cells"));
        }
        let mut unique = BTreeSet::new();
        for &cell in cells {
            poll(&mut keep_going)?;
            if cell >= cell_count || !unique.insert(cell) {
                return Err(TetError::Invalid("region cell is out of range or repeated"));
            }
        }
        poll(&mut keep_going)?;
        Ok(Self { cell_count, cells: unique.into_iter().collect() })
    }

    /// Original whole-domain element count, not selected cardinality.
    #[must_use]
    pub const fn domain_cell_count(&self) -> usize { self.cell_count }

    /// Original element indices, sorted and distinct; no vertex averaging.
    #[must_use]
    pub fn cells(&self) -> &[usize] { &self.cells }
}

/// Mean temperature over the explicit cell union, not over the whole assembly.
/// The two energy bounds and their residual correction still use ALL cells.
#[derive(Debug, Clone)]
pub struct RegionMeanBound {
    pub selection: RegionMeanSelection,
    /// Outward volume of only the selected cells. Division happens LAST.
    pub region_volume: Iv,
    /// Enclosure of the continuum integral over the region / region volume.
    pub enclosure: Iv,
    /// Same region average of the supplied P1 field, evaluated outward.
    pub candidate_mean: Iv,
    /// Unnormalized region integral and complete-domain primal/dual evidence.
    /// `integral.domain_volume` is intentionally the whole-domain volume.
    pub integral: GoalBound,
}

impl AffineSourceTetProblem<'_> {
    /// Bound a component/region volume mean with the original tensor, source,
    /// boundary and matching-contact model. The dual source is exactly ONE on
    /// selected cells and ZERO elsewhere; no rounded reciprocal-volume source
    /// and no artificial boundary are introduced at the selection boundary.
    ///
    /// Any finite admissible dual candidate is accepted, including an inexact
    /// solve or a candidate for another load: the majorant bounds its error
    /// against the correct region-supported dual. Such a candidate may be loose.
    /// No whole-domain certificate is merely rescaled or re-labelled.
    ///
    /// # Errors
    /// Selection/count mismatch, original physical/geometry/tensor refusals,
    /// exhausted flux budget, unbounded arithmetic or cancellation yield no
    /// partial bound. This is not a nodal maximum, CAD, nonlinear, nonmatching
    /// interface, material-uncertainty or physical-validation certificate.
    pub fn region_mean_bound(
        &self,
        candidate: &[f64],
        dual_candidate: &[f64],
        selection: &RegionMeanSelection,
        budget: FluxBudget,
        mut keep_going: impl FnMut() -> bool,
    ) -> Result<RegionMeanBound, TetError> {
        poll(&mut keep_going)?;
        if selection.cell_count != self.tets.len() {
            return Err(TetError::Invalid("region selection belongs to a different domain cell count"));
        }
        let tensors = self.prepare(budget, &mut keep_going)?;
        let problem = self.problem(&tensors);
        // Use the same outward geometry and all the original boundary/contact
        // admission. Selected cells are not isolated from the rest of the PDE.
        let (cells, _) = build(&problem, candidate, budget, &mut keep_going)?;
        let mut volume = Iv::zero();
        let mut weights = vec![0.0; cells.len()];
        for &cell in &selection.cells {
            poll(&mut keep_going)?;
            volume = volume.add(cells[cell].volume);
            weights[cell] = 1.0;
        }
        if volume.is_unbounded() || volume.lo <= 0.0 { return Err(TetError::Unbounded); }
        // Geometry is reconstructible; do not retain another full copy while
        // the original goal owner constructs its primal and dual fluxes.
        drop(cells);
        let integral = goal::goal_bound_impl(&problem, candidate, dual_candidate,
            &weights, budget, &mut keep_going)?;
        let enclosure = integral.enclosure.div_pos(volume);
        let candidate_mean = integral.candidate_value.div_pos(volume);
        poll(&mut keep_going)?;
        if enclosure.is_unbounded() || candidate_mean.is_unbounded() {
            return Err(TetError::Unbounded);
        }
        Ok(RegionMeanBound { selection: selection.clone(), region_volume: volume,
            enclosure, candidate_mean, integral })
    }
}
