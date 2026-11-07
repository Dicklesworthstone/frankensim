//! Independently varying physical shell loads on the original accepted path.
//!
//! Reference fields are scaled by the next declared `load_path` factor, exactly
//! as in the original native problem. Unit factors therefore admit arbitrary
//! absolute displacement/force histories, including phase-shifted membrane,
//! bending and shear. The support DOFs, source laws, geometry and time owner stay
//! fixed. Supplied loads never replace accepted stress or strain observations.

use super::{
    PlasticShellReferenceInertia, ShellFatigueError, ShellFatiguePath, ShellInertialEnergy, poll,
};
use fs_exec::Cx;
use fs_material::fatigue::cyclic::SourceBoundFatigueLaw;

/// Reference boundary fields for one original load-path entry.
///
/// The next native load factor multiplies both fields. Declare unit factors
/// when these are absolute endpoint values. Loads are dead nodal forces/couples
/// in the fixed global frame, with the original shell's native SI convention.
#[derive(Debug, Clone, Copy)]
pub struct ShellFatigueLoading<'a> {
    /// Values in exactly the original `problem.prescribed` order: translations
    /// in metres and rotations in radians. Their DOFs cannot be replaced.
    pub prescribed_values: &'a [f64],
    /// `(6 * node + component, load)` pairs: newtons for components 0..3 and
    /// newton-metres for components 3..6. Native additive duplicates are allowed.
    pub nodal_loads: &'a [(usize, f64)],
}

impl<'p, 'm, L: SourceBoundFatigueLaw> ShellFatiguePath<'p, 'm, L> {
    /// Advance one native static target with independently supplied boundary
    /// fields. All normal and shear observations publish with its accepted
    /// solution. `false` means the original load path has no remaining entry.
    pub fn advance_next_with_loading(
        &mut self,
        loading: ShellFatigueLoading<'_>,
        cx: &Cx<'_>,
    ) -> Result<bool, ShellFatigueError> {
        self.admit_advance(cx)?;
        let prescribed = self.prescribed_for_loading(loading, cx)?;
        self.advance_next_with_fields(&prescribed, loading.nodal_loads, false, cx)
    }

    /// Advance a static endpoint and close every declared representative shear
    /// window in the same transaction. A nonclosed cycle or late calibration
    /// refusal leaves mechanics and every fatigue history at the prior endpoint.
    pub fn advance_next_cycle_with_loading(
        &mut self,
        loading: ShellFatigueLoading<'_>,
        cx: &Cx<'_>,
    ) -> Result<bool, ShellFatigueError> {
        self.require_shear(cx)?;
        let prescribed = self.prescribed_for_loading(loading, cx)?;
        self.advance_next_with_fields(&prescribed, loading.nodal_loads, true, cx)
    }

    /// Advance the next native backward-Euler endpoint with independent boundary
    /// fields and the original reference inertia. Physical velocity, reactions
    /// and energy come from that same native interval, including support motion.
    pub fn advance_next_backward_euler_with_loading(
        &mut self,
        reference: &'p PlasticShellReferenceInertia<'m>,
        dt_s: f64,
        loading: ShellFatigueLoading<'_>,
        cx: &Cx<'_>,
    ) -> Result<Option<ShellInertialEnergy>, ShellFatigueError> {
        self.admit_advance(cx)?;
        let prescribed = self.prescribed_for_loading(loading, cx)?;
        self.advance_next_backward_euler_with_fields(
            reference,
            dt_s,
            &prescribed,
            loading.nodal_loads,
            false,
            cx,
        )
    }

    /// Advance and close the representative shear windows atomically with the
    /// original physical-time solve. Failure preserves accepted displacement,
    /// velocity, constitutive state, elapsed time and every fatigue counter.
    pub fn advance_next_backward_euler_cycle_with_loading(
        &mut self,
        reference: &'p PlasticShellReferenceInertia<'m>,
        dt_s: f64,
        loading: ShellFatigueLoading<'_>,
        cx: &Cx<'_>,
    ) -> Result<Option<ShellInertialEnergy>, ShellFatigueError> {
        self.require_shear(cx)?;
        let prescribed = self.prescribed_for_loading(loading, cx)?;
        self.advance_next_backward_euler_with_fields(
            reference,
            dt_s,
            &prescribed,
            loading.nodal_loads,
            true,
            cx,
        )
    }

    fn prescribed_for_loading(
        &self,
        loading: ShellFatigueLoading<'_>,
        cx: &Cx<'_>,
    ) -> Result<Vec<(usize, f64)>, ShellFatigueError> {
        poll(cx)?;
        if loading.prescribed_values.len() != self.problem.prescribed.len() {
            return Err(ShellFatigueError::Invalid(
                "explicit shell loading must cover exactly the original prescribed DOFs",
            ));
        }
        let dofs = self
            .problem
            .mesh
            .nodes
            .len()
            .checked_mul(6)
            .ok_or(ShellFatigueError::Budget)?;
        // Geometry is already prepared. Admit new physical load fields here
        // before the original solve can index a force or assemble its endpoint.
        for (index, &(dof, value)) in loading.nodal_loads.iter().enumerate() {
            if index % 256 == 0 {
                poll(cx)?;
            }
            if dof >= dofs || !value.is_finite() {
                return Err(ShellFatigueError::Invalid(
                    "explicit shell nodal forces/couples must be finite and inside the original mesh",
                ));
            }
        }
        let mut prescribed = Vec::new();
        prescribed
            .try_reserve_exact(loading.prescribed_values.len())
            .map_err(|_| ShellFatigueError::Budget)?;
        for (index, (&(dof, _), &value)) in self
            .problem
            .prescribed
            .iter()
            .zip(loading.prescribed_values)
            .enumerate()
        {
            if index % 256 == 0 {
                poll(cx)?;
            }
            if !value.is_finite() {
                return Err(ShellFatigueError::Invalid(
                    "explicit shell prescribed coordinates must be finite",
                ));
            }
            prescribed.push((dof, value));
        }
        poll(cx)?;
        Ok(prescribed)
    }
}
