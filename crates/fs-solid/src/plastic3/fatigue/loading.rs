//! Independent physical boundary fields on native cyclic solid histories.
//!
//! Each field is a reference pattern scaled by that endpoint's original load
//! factor. Unit factors admit absolute displacement/force histories, including
//! independently phased normal and shear loads. Support topology, material
//! sources and the native constitutive/time integration owners stay fixed.

use fs_exec::Cx;
use fs_material::SmallStrainLaw;
use fs_material::fatigue::cyclic::SourceBoundFatigueLaw;

use super::{J2FatigueError, J2FatiguePath, SmallStrainTetProblem, poll};

/// Reference boundary fields for one native fatigue endpoint.
///
/// The endpoint's load factor multiplies both fields. Supply a unit factor
/// when these are absolute endpoint values. Forces are dead nodal forces in
/// the original fixed mesh frame; this does not admit follower loads.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct J2FatigueLoading<'a> {
    /// Metres in exactly the original `problem.prescribed_m` order. Only the
    /// values may change; the original constrained DOFs cannot be replaced.
    pub prescribed_values_m: &'a [f64],
    /// `(3 * node + component, force_in_newtons)` pairs. Native additive
    /// duplicates are allowed, including loads on supports for reactions.
    pub nodal_forces_n: &'a [(usize, f64)],
}

impl<'p, 'm, L: SourceBoundFatigueLaw> J2FatiguePath<'p, 'm, L> {
    /// Accept a nonempty block of independently supplied boundary fields. Each
    /// entry consumes the next original `load_path` factor and calls the same
    /// native solve with the previous complete constitutive state. Every actual
    /// stress/strain endpoint enters all declared fatigue banks.
    ///
    /// The whole block is one transaction: invalid loading, nonconvergence,
    /// cancellation or any fatigue refusal preserves the previous accepted
    /// mechanics, increment index and every fatigue history. Calling ordinary
    /// `advance` afterwards returns to the original reference fields.
    pub fn advance_with_loading(
        &mut self,
        loading: &[J2FatigueLoading<'_>],
        cx: &Cx<'_>,
    ) -> Result<(), J2FatigueError> {
        self.advance_block(loading.len(), Some(loading), false, cx)
    }

    /// Accept an independent-loading block and close all representative SWT
    /// and Fatemi-Socie windows in the same transaction. Nonclosed stress or
    /// strain histories refuse the entire block; no cycle is manufactured by
    /// returning just the applied loads to their earlier values.
    pub fn advance_cycle_with_loading(
        &mut self,
        loading: &[J2FatigueLoading<'_>],
        cx: &Cx<'_>,
    ) -> Result<(), J2FatigueError> {
        self.observer.require_critical(cx)?;
        self.advance_block(loading.len(), Some(loading), true, cx)
    }
}

impl J2FatigueLoading<'_> {
    /// Admit reference fields before either native static or inertial work.
    /// The native assembly remains responsible for accumulated-force overflow
    /// and constitutive/geometry admission; this check does not replace it.
    pub(crate) fn validate<L: SmallStrainLaw>(
        self,
        problem: &SmallStrainTetProblem<'_, L>,
        factor: f64,
        cx: &Cx<'_>,
    ) -> Result<(), J2FatigueError> {
        poll(cx)?;
        if self.prescribed_values_m.len() != problem.prescribed_m.len() {
            return Err(J2FatigueError::Invalid(
                "explicit solid loading must cover exactly the original prescribed DOFs",
            ));
        }
        if !factor.is_finite() {
            return Err(J2FatigueError::Invalid(
                "explicit solid loading needs a finite endpoint load factor",
            ));
        }
        let dofs = problem
            .nodes_m
            .len()
            .checked_mul(3)
            .ok_or(J2FatigueError::Budget)?;
        for (index, &value) in self.prescribed_values_m.iter().enumerate() {
            if index % 256 == 0 {
                poll(cx)?;
            }
            if !value.is_finite() || !(factor * value).is_finite() {
                return Err(J2FatigueError::Invalid(
                    "explicit solid prescribed displacement must be finite and representable at its load factor",
                ));
            }
        }
        for (index, &(dof, value)) in self.nodal_forces_n.iter().enumerate() {
            if index % 256 == 0 {
                poll(cx)?;
            }
            if dof >= dofs || !value.is_finite() || !(factor * value).is_finite() {
                return Err(J2FatigueError::Invalid(
                    "explicit solid nodal force must lie inside the original mesh and be finite and representable at its load factor",
                ));
            }
        }
        poll(cx)
    }
}
