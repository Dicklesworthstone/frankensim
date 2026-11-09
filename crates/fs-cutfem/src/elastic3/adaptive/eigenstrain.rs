//! Free-strain mechanics through the original adaptive Q1 constraint map.
//!
//! Preparation retains physical cell quadrature and original source identity.
//! Forces use T^T; observations and duals use the SAME homogeneous-clamp-aware
//! T. Refining or rebuilding an operator requires preparation on its new owner.

use std::ops::ControlFlow;

use fs_solver::op::LinearOp;

use super::AdaptiveElasticity3;
use crate::elastic3::ElasticityError3;
use crate::elastic3::eigenstrain::{EigenstrainResponse3, PreparedEigenstrain3};
use crate::elastic3::stress::BulkStressPullback3;

mod equilibrium;

fn checkpoint(poll: &mut impl FnMut() -> ControlFlow<()>) -> Result<(), ElasticityError3> {
    if poll().is_break() {
        Err(ElasticityError3::Cancelled)
    } else {
        Ok(())
    }
}

impl AdaptiveElasticity3 {
    /// Prepare tensor-shear free strains on original physical cut quadrature.
    /// Cell indices follow [`Self::leaves`]. No cell-average thermal field or
    /// new quadrature is introduced by the hanging-node reduction.
    pub fn prepare_eigenstrain(
        &self,
        free_strain: &dyn Fn(usize, [f64; 3]) -> Result<[f64; 6], ElasticityError3>,
        max_points: usize,
        poll: impl FnMut() -> ControlFlow<()>,
    ) -> Result<PreparedEigenstrain3, ElasticityError3> {
        self.raw.prepare_eigenstrain(free_strain, max_points, poll)
    }

    /// Refuse another geometry/material/boundary owner's preparation.
    /// Changes of the existing positive stiffness scales preserve the binding.
    pub fn validate_eigenstrain(
        &self,
        prepared: &PreparedEigenstrain3,
    ) -> Result<(), ElasticityError3> {
        self.raw.validate_eigenstrain(prepared)
    }

    /// Apply the original hanging-node transpose to the scaled free-strain load.
    /// Strong homogeneous master clamps carry no work. Embedded eigenstress
    /// traction has already been included by the physical preparation owner.
    pub fn eigenstrain_load(
        &self,
        prepared: &PreparedEigenstrain3,
        mut poll: impl FnMut() -> ControlFlow<()>,
    ) -> Result<Vec<f64>, ElasticityError3> {
        checkpoint(&mut poll)?;
        let raw = self.raw.eigenstrain_load(prepared, &mut poll)?;
        let mut rhs = vec![0.0; self.n()];
        self.restrict(&raw, &mut rhs);
        for chunk in rhs.chunks(192) {
            checkpoint(&mut poll)?;
            if chunk.iter().any(|value| !value.is_finite()) {
                return Err(ElasticityError3::Invalid(
                    "adaptive eigenstrain load overflow",
                ));
            }
        }
        checkpoint(&mut poll)?;
        Ok(rhs)
    }

    /// Recover thermal stress and elastic storage from independent masters.
    /// Reconstruct only the current cell's coefficients. Total kinematics and
    /// elastic material domains remain separate in the original response owner;
    /// both compatible and B-bar-projected elastic strains are admitted there.
    pub fn eigenstrain_response(
        &self,
        prepared: &PreparedEigenstrain3,
        u: &[f64],
        maximum_kinematic_strain: f64,
        max_points: usize,
        poll: impl FnMut() -> ControlFlow<()>,
    ) -> Result<Vec<EigenstrainResponse3>, ElasticityError3> {
        self.raw.eigenstrain_response_mapped::<true>(
            prepared,
            u,
            self.n(),
            maximum_kinematic_strain,
            max_points,
            |cell, field| {
                std::array::from_fn(|i| self.physical_coefficient(field, cell.nodes[i / 3], i % 3))
            },
            poll,
        )
    }

    /// Admit the original thermal state at retained points without a global
    /// displacement copy or observation array. The raw owner performs the same
    /// separate total-kinematics and elastic-material checks as response.
    pub fn admit_eigenstrain(
        &self,
        prepared: &PreparedEigenstrain3,
        u: &[f64],
        maximum_kinematic_strain: f64,
        max_points: usize,
        poll: impl FnMut() -> ControlFlow<()>,
    ) -> Result<(), ElasticityError3> {
        self.raw
            .eigenstrain_response_mapped::<false>(
                prepared,
                u,
                self.n(),
                maximum_kinematic_strain,
                max_points,
                |cell, field| {
                    std::array::from_fn(|i| {
                        self.physical_coefficient(field, cell.nodes[i / 3], i % 3)
                    })
                },
                poll,
            )
            .map(|_| ())
    }

    /// Exact scale-load transpose in the independent master space.
    /// Cell-local conforming reconstruction avoids expanding a global dual.
    /// No extra material scale or arbitrary strain-domain check is applied.
    pub fn eigenstrain_load_scale_pullback(
        &self,
        prepared: &PreparedEigenstrain3,
        dual: &[f64],
        poll: impl FnMut() -> ControlFlow<()>,
    ) -> Result<Vec<f64>, ElasticityError3> {
        self.raw.eigenstrain_load_scale_pullback_mapped(
            prepared,
            dual,
            self.n(),
            |cell, field| {
                std::array::from_fn(|i| self.physical_coefficient(field, cell.nodes[i / 3], i % 3))
            },
            poll,
        )
    }

    /// Exact thermal stress VJP with the original hanging-node transpose.
    /// The direct material-scale derivative is unchanged by reconstruction;
    /// displacement cotangents are restricted to independent master coordinates.
    pub fn eigenstrain_stress_pullback(
        &self,
        prepared: &PreparedEigenstrain3,
        u: &[f64],
        maximum_kinematic_strain: f64,
        derivatives: &[[f64; 6]],
        mut poll: impl FnMut() -> ControlFlow<()>,
    ) -> Result<BulkStressPullback3, ElasticityError3> {
        checkpoint(&mut poll)?;
        let full = self.physical_displacements(u)?;
        let mut result = self.raw.eigenstrain_stress_pullback(
            prepared,
            &full,
            maximum_kinematic_strain,
            derivatives,
            &mut poll,
        )?;
        let mut displacement = vec![0.0; self.n()];
        self.restrict(&result.displacement, &mut displacement);
        for chunk in displacement.chunks(192) {
            checkpoint(&mut poll)?;
            if chunk.iter().any(|value| !value.is_finite()) {
                return Err(ElasticityError3::Invalid(
                    "adaptive eigenstrain stress pullback overflow",
                ));
            }
        }
        checkpoint(&mut poll)?;
        result.displacement = displacement;
        Ok(result)
    }
}
