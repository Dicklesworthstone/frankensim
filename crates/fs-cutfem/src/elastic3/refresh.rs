//! Transactional material reassembly on the original reference geometry.
//!
//! Only constitutive values change. Retained quadrature, support selection,
//! ghost stencils and displacement coordinates are neither rebuilt nor copied.

use super::*;
use std::mem::{size_of, swap};
use std::sync::Arc;

/// Explicit native work/storage bounds before material callbacks or allocation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MaterialRefreshLimits3 {
    /// Maximum original active cells, each queried exactly once.
    pub max_cells: usize,
    /// Sum of retained bulk points, selected embedded-support points, and ghost
    /// points visited by the refresh. No new geometric quadrature is produced.
    pub max_points: usize,
    /// Conservative fixed native staging payload: replacement cell blocks and
    /// kinematics, one fixed-size native tensor cache per cell, and ghost
    /// coefficients. This is NOT a total heap cap. It excludes original dynamic
    /// thermal/source payloads and receipts, caller-owned material resolution/
    /// query allocations, allocator bookkeeping and fixed kernel stack.
    /// Existing geometry is borrowed.
    pub max_scratch_bytes: usize,
}

impl Default for MaterialRefreshLimits3 {
    fn default() -> Self {
        Self {
            max_cells: 32_768,
            max_points: 1_000_000,
            max_scratch_bytes: 256 * 1024 * 1024,
        }
    }
}

/// Outcome of a valid native material trial.
///
/// Rejection returns owned scratch for a coupling iteration while restoring the
/// original material caches and exact preparation token. Any native preparation
/// contained in rejected scratch belongs to the rejected candidate and is stale
/// against the restored operator. Acceptance keeps the candidate owner and its
/// returned preparations only after the final cancellation check.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MaterialRefreshDecision3<R> {
    Accept(R),
    Reject(R),
}

struct CellRefresh3 {
    material: ReferenceMaterial3,
    kinematics: CellKinematics3,
    stiffness: [[f64; 24]; 24],
}

/// Only constitutive caches and their ownership token move between accepted
/// and trial states. After installation these same buffers hold the old state.
struct MaterialRefresh3 {
    cells: Vec<CellRefresh3>,
    ghost_weights: Vec<[f64; 2]>,
    binding: Arc<()>,
}

impl MaterialRefresh3 {
    fn swap_with(&mut self, source: &mut CutElasticity3) {
        // Both vectors were staged from this exact owner. These swaps perform
        // no allocation, indexing, callbacks or destruction of stored values.
        // Geometry, scales and reference inventory retain their exact storage.
        for (cell, replacement) in source.cells.iter_mut().zip(&mut self.cells) {
            swap(&mut cell.material, &mut replacement.material);
            swap(&mut cell.kinematics, &mut replacement.kinematics);
            swap(&mut cell.stiffness, &mut replacement.stiffness);
        }
        for (ghost, weights) in source.ghosts.iter_mut().zip(&mut self.ghost_weights) {
            swap(&mut ghost.weights, weights);
        }
        swap(&mut source.material_body_binding, &mut self.binding);
    }
}

// The guard borrows the COMPLETE original owner. In particular, an adaptive
// callback must retain access to its original master-space solve and maps.
pub(in crate::elastic3) trait MaterialRefreshOwner3 {
    fn material_refresh_raw(&self) -> &CutElasticity3;
    fn material_refresh_raw_mut(&mut self) -> &mut CutElasticity3;
}

impl MaterialRefreshOwner3 for CutElasticity3 {
    fn material_refresh_raw(&self) -> &CutElasticity3 {
        self
    }

    fn material_refresh_raw_mut(&mut self) -> &mut CutElasticity3 {
        self
    }
}

struct MaterialRefreshGuard3<'a, O: MaterialRefreshOwner3> {
    owner: &'a mut O,
    previous: MaterialRefresh3,
    rollback: bool,
}

impl<O: MaterialRefreshOwner3> Drop for MaterialRefreshGuard3<'_, O> {
    fn drop(&mut self) {
        if self.rollback {
            self.previous.swap_with(self.owner.material_refresh_raw_mut());
        }
    }
}

/// Existing accepting callbacks retain their concrete native error type.
pub(in crate::elastic3) fn refresh_material_owner<O: MaterialRefreshOwner3, R>(
    owner: &mut O,
    field: ElasticMaterialField3<'_>,
    limits: MaterialRefreshLimits3,
    attempt: impl FnOnce(&O) -> Result<R, ElasticityError3>,
    checkpoint: impl FnMut() -> ControlFlow<()>,
) -> Result<R, ElasticityError3> {
    let decision = refresh_material_owner_with_decision(
        owner, field, limits,
        |candidate| attempt(candidate).map(MaterialRefreshDecision3::Accept),
        checkpoint,
    )?;
    // This adapter constructs only Accept; extracting its payload introduces
    // no extra cancellation boundary or independent publication decision.
    match decision {
        MaterialRefreshDecision3::Accept(value) | MaterialRefreshDecision3::Reject(value) => Ok(value),
    }
}

/// Shared single-owner decision transaction on either native displacement space.
pub(in crate::elastic3) fn refresh_material_owner_with_decision<O, R, E>(
    owner: &mut O,
    field: ElasticMaterialField3<'_>,
    limits: MaterialRefreshLimits3,
    attempt: impl FnOnce(&O) -> Result<MaterialRefreshDecision3<R>, E>,
    mut checkpoint: impl FnMut() -> ControlFlow<()>,
) -> Result<MaterialRefreshDecision3<R>, E>
where
    O: MaterialRefreshOwner3,
    E: From<ElasticityError3>,
{
    let mut poll = || {
        if checkpoint().is_break() {
            Err(ElasticityError3::Cancelled)
        } else {
            Ok(())
        }
    };
    let replacement = owner
        .material_refresh_raw()
        .stage_material_refresh(field, limits, &mut poll)?;
    // Arm rollback before installation. The swap itself cannot invoke fallible
    // or user code; every callback and the final poll run under this guard.
    let mut trial = MaterialRefreshGuard3 {
        owner,
        previous: replacement,
        rollback: true,
    };
    trial.previous.swap_with(trial.owner.material_refresh_raw_mut());
    let decision = attempt(&*trial.owner)?;
    poll()?;
    if matches!(&decision, MaterialRefreshDecision3::Accept(_)) {
        trial.rollback = false;
    }
    // A rejected payload moves to the caller only after this guard has restored
    // the accepted caches/token. Native and caller errors retain the same guard.
    Ok(decision)
}

fn invalid(reason: &'static str) -> ElasticityError3 {
    ElasticityError3::Invalid(reason)
}

impl CutElasticity3 {
    fn admit_material_refresh(
        &self,
        limits: MaterialRefreshLimits3,
        checkpoint: &mut impl FnMut() -> Result<(), ElasticityError3>,
    ) -> Result<(), ElasticityError3> {
        checkpoint()?;
        if self.cells.len() > limits.max_cells {
            return Err(invalid("material refresh exceeds the cell allowance"));
        }
        if self.reference_cell_masses_kg.is_none() {
            return Err(invalid("material refresh requires original native reference density"));
        }
        let mut points = self.ghosts.len();
        for (index, cell) in self.cells.iter().enumerate() {
            checkpoint()?;
            if cell.material.source_state().is_none() {
                return Err(invalid("material refresh requires original native reference density"));
            }
            let surface_points = self.embedded.as_ref()
                .map_or(0, |support| support.selected_points(index).len());
            points = points.checked_add(cell.rules.bulk().len())
                .and_then(|count| count.checked_add(surface_points))
                .ok_or(invalid("material refresh point count overflow"))?;
            if points > limits.max_points {
                return Err(invalid("material refresh exceeds the retained-point allowance"));
            }
        }
        // Charge a separate fixed native cache even for a uniform provider,
        // which shares one Arc payload. Arc counters and the new owner token
        // are included. Dynamic source payloads and allocator overhead are not.
        let cell_bytes = size_of::<CellRefresh3>()
            .checked_add(size_of::<material::PreparedMaterial3>())
            .and_then(|bytes| bytes.checked_add(2 * size_of::<usize>()))
            .ok_or(invalid("material refresh workspace overflow"))?;
        let bytes = self.cells.len().checked_mul(cell_bytes)
            .and_then(|bytes| self.ghosts.len().checked_mul(size_of::<[f64; 2]>())
                .and_then(|ghosts| bytes.checked_add(ghosts)))
            .and_then(|bytes| bytes.checked_add(2 * size_of::<usize>()))
            .ok_or(invalid("material refresh workspace overflow"))?;
        if bytes > limits.max_scratch_bytes {
            return Err(invalid("material refresh exceeds the staging-byte allowance"));
        }
        checkpoint()
    }

    /// Replace the complete native constitutive field on this SAME reference
    /// geometry. Reassembles bulk, B-bar, selected Nitsche and both ghost sides
    /// from the original retained rules and stencils; no SDF, patch, clamp or
    /// geometric quadrature callback is invoked. Each active-cell material is
    /// resolved once and retained, with no material queries in later applies.
    ///
    /// This is a fixed-reference thermal/material update: every cell's new
    /// native reference-density bits must equal its original density. A current
    /// thermal density cannot change the mass of an undeformed closed specimen.
    /// Legacy states without native inventory refuse. Stiffness-design scales,
    /// physical inventory, geometry, supports and formulation remain unchanged.
    /// B-bar retains its original physical-volume projection and requires the
    /// new source to remain explicitly isotropic; no isotropic fit is inferred.
    ///
    /// All dependent values are staged before publication. Any returned error
    /// or cancellation preserves the old operator and ALL of its preparations.
    /// A successful update rotates ownership: mass, body-load, eigenstrain and
    /// prestress preparations must be recreated even when their dimensions or
    /// reference density match. Source/material domains are never enlarged.
    ///
    /// This commits an OPERATOR update, not a complete heat/mechanical step.
    /// It does not advance temperature, replace a source query, differentiate
    /// C(T), or claim atomicity for a later equilibrium/thermal solve.
    /// Use [Self::refresh_material_field_with] when dependent solves or checks
    /// must succeed before this operator update can be accepted.
    pub fn refresh_material_field(
        &mut self,
        field: ElasticMaterialField3<'_>,
        limits: MaterialRefreshLimits3,
        checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<(), ElasticityError3> {
        self.refresh_material_field_with(field, limits, |_| Ok(()), checkpoint)
    }

    /// Accept a native material update only after a dependent attempt succeeds.
    ///
    /// The immutable candidate is this SAME native displacement owner with its
    /// fully staged material, bulk/B-bar, selected Nitsche and ghost caches.
    /// The callback can use existing thermal/eigenstrain loads, equilibrium
    /// solvers, physical response, adjoints, mass and prestress preparation.
    /// Its own physical/conservation checks and solve budgets remain required.
    /// Staging limits and reference-density rules match [Self::refresh_material_field].
    ///
    /// Returning Err, cancellation after the callback, or an unwinding panic in
    /// the callback/final checkpoint restores the old caches and exact ownership
    /// token. Every previously valid preparation remains valid after rejection.
    /// Successful return keeps the candidate token; owned preparations returned
    /// in R remain valid against the accepted operator, while old ones are stale.
    /// R cannot borrow the temporary candidate view.
    ///
    /// Geometry, scales, supports and reference inventory are never changed.
    /// The same staging buffers hold the old caches during the attempt; no
    /// geometry copy or second old-cache allocation is needed. Final acceptance
    /// and rollback have no callbacks, allocations or fallible work.
    ///
    /// Atomicity covers this operator and the returned owned payload. Mutations
    /// to caller-owned state, source accounts or observations inside the callback
    /// cannot be undone here: return proposed external state and publish it only
    /// after success. Nesting independent transactions does not establish a
    /// multi-operator coupled commit or a whole-system conservation claim.
    pub fn refresh_material_field_with<R>(
        &mut self,
        field: ElasticMaterialField3<'_>,
        limits: MaterialRefreshLimits3,
        attempt: impl FnOnce(&Self) -> Result<R, ElasticityError3>,
        checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<R, ElasticityError3> {
        refresh_material_owner(self, field, limits, attempt, checkpoint)
    }

    /// Evaluate a native material trial and explicitly accept or reject it.
    ///
    /// Both valid decisions pass a final checkpoint. Accept keeps the candidate
    /// material owner and returns its owned payload. Reject restores the prior
    /// caches and exact token before returning its owned scratch; preparations
    /// inside that scratch are stale against the restored operator. This lets
    /// coupling iterations inspect rejected physical proposals without exposing
    /// a half-accepted material field.
    ///
    /// Native admission/staging/cancellation errors convert through E::from.
    /// Callback errors retain E unchanged; errors and unwinding callback or
    /// final-checkpoint panics restore the original owner. R cannot borrow the
    /// candidate view. Reference inventory, geometry, constraints and scales are
    /// fixed. External callback side effects are outside this transaction, as
    /// specified by [Self::refresh_material_field_with].
    pub fn refresh_material_field_with_decision<R, E>(
        &mut self,
        field: ElasticMaterialField3<'_>,
        limits: MaterialRefreshLimits3,
        attempt: impl FnOnce(&Self) -> Result<MaterialRefreshDecision3<R>, E>,
        checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<MaterialRefreshDecision3<R>, E>
    where
        E: From<ElasticityError3>,
    {
        refresh_material_owner_with_decision(self, field, limits, attempt, checkpoint)
    }

    fn stage_material_refresh(
        &self,
        field: ElasticMaterialField3<'_>,
        limits: MaterialRefreshLimits3,
        poll: &mut impl FnMut() -> Result<(), ElasticityError3>,
    ) -> Result<MaterialRefresh3, ElasticityError3> {
        self.admit_material_refresh(limits, poll)?;
        let provider = MaterialProvider3::from_field(field)?;
        poll()?;
        let mut cells = Vec::with_capacity(self.cells.len());
        for (index, cell) in self.cells.iter().enumerate() {
            poll()?;
            let material = provider.at(cell.bounds)?;
            poll()?;
            let original = cell.material.source_state()
                .expect("original native density was admitted");
            let refreshed = material.source_state()
                .ok_or(invalid("material refresh requires native reference density"))?;
            if original.density_kg_m3().to_bits() != refreshed.density_kg_m3().to_bits() {
                return Err(invalid("material refresh cannot change reference density"));
            }
            let kinematics = cell.kinematics.with_material(&material)?;
            let mut stiffness = assemble_cell_stiffness_controlled(
                cell.bounds, &cell.rules, &material, kinematics, poll)?;
            if let Some(support) = &self.embedded {
                support.reassemble_cell_stiffness(index, cell, &material, kinematics,
                    &mut stiffness, poll)?;
            }
            cells.push(CellRefresh3 { material, kinematics, stiffness });
        }
        let mut ghost_weights = Vec::with_capacity(self.ghosts.len());
        for ghost in &self.ghosts {
            poll()?;
            let weights = ghost.cells.map(|index| {
                let cell = &cells[index];
                ghost.geometry.weight(cell.kinematics.ghost_modulus(&cell.material))
            });
            if weights.iter().any(|weight| !weight.is_finite() || *weight <= 0.0) {
                return Err(invalid("material refresh ghost coefficient overflow"));
            }
            ghost_weights.push(weights);
        }
        let binding = Arc::new(());
        poll()?;
        Ok(MaterialRefresh3 { cells, ghost_weights, binding })
    }
}
