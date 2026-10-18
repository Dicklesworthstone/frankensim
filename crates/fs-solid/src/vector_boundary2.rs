//! Typed vector histories bound to the existing body-fitted 2-D elasticity owner.
//!
//! These are borrowed solver bindings, not a second serialized scenario schema.
//! A payload prescribes a spatially uniform vector on an explicit mesh patch;
//! shape, interpolation, outside policy and moving frames remain fs-scenario's.
//! Coordinates are WORLD x/y, lengths are metres, traction is Pa per unit
//! out-of-plane thickness. Any nonzero WORLD z component refuses: no hidden
//! projection, pressure-normal convention, inertial force or follower load.
//! Each call solves one quasi-static endpoint, not a transient mechanics step.

use crate::{LinearProblem, Mesh2, Patch, SolidError};
use fs_qty::Dims;
use fs_scenario::{BcKind, BcValue, BoundaryCondition, FrameTree, Physics};
use fs_scenario::payload::{OrientationParity, Payload, QuantityContract, VectorPayload};
use fs_scenario::vector::VectorEvaluationError;
use std::collections::{BTreeMap, BTreeSet};
use std::ops::ControlFlow;

const DISPLACEMENT: Dims = Dims([1, 0, 0, 0, 0, 0]);
const TRACTION: Dims = Dims([-1, 1, -2, 0, 0, 0]);

/// Explicit patch and role for an existing, immutable vector payload.
#[derive(Debug, Clone, Copy)]
pub enum VectorBoundary2<'a> {
    /// Both displacement components, in metres.
    Displacement { patch: Patch, value: &'a VectorPayload },
    /// Dead-load traction components, in Pa; not a scalar pressure.
    Traction { patch: Patch, value: &'a VectorPayload },
}

impl<'a> VectorBoundary2<'a> {
    /// Bind the canonical Elasticity/Neumann typed-vector row to an explicitly
    /// resolved mesh patch. Legacy scalar Traction/Dirichlet rows are NOT
    /// reinterpreted as vectors. Prescribed vector displacements use the explicit
    /// `Displacement` binding above; this does not widen the native IR grammar.
    pub fn scenario_traction(patch: Patch, row: &'a BoundaryCondition) -> Result<Self, VectorBoundaryError2> {
        let mut violations = Vec::new();
        row.check(&mut violations);
        if let Some(first) = violations.first() {
            return Err(VectorBoundaryError2::Scenario(first.what.clone()));
        }
        if row.physics != Physics::Elasticity || row.kind != BcKind::Neumann {
            return Err(VectorBoundaryError2::Invalid("only the explicit Elasticity/Neumann vector row lowers through this constructor"));
        }
        match &row.value {
            Some(BcValue::Typed(Payload::Vector(value))) => Ok(Self::Traction { patch, value }),
            _ => Err(VectorBoundaryError2::Invalid("a typed vector payload is required")),
        }
    }
}

/// Errors before publication of a complete equilibrium field.
#[derive(Debug, Clone, PartialEq)]
pub enum VectorBoundaryError2 {
    Invalid(&'static str),
    Scenario(String),
    Vector { patch: Patch, error: VectorEvaluationError },
    Quantity { patch: Patch, expected: Dims, found: QuantityContract },
    OutOfPlane { patch: Patch, value: f64 },
    ConflictingDisplacement { node: usize },
    Solver(SolidError),
    Cancelled,
}
impl std::fmt::Display for VectorBoundaryError2 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "2-D vector boundary solve refused: {self:?}")
    }
}
impl std::error::Error for VectorBoundaryError2 {}

fn poll(checkpoint: &mut impl FnMut() -> ControlFlow<()>) -> Result<(), VectorBoundaryError2> {
    if checkpoint().is_break() { Err(VectorBoundaryError2::Cancelled) } else { Ok(()) }
}
fn patch_nodes(mesh: &Mesh2, patch: Patch) -> Result<Vec<usize>, VectorBoundaryError2> {
    if mesh.patches.iter().filter(|(name, _)| *name == patch).count() != 1 {
        return Err(VectorBoundaryError2::Invalid("boundary patch must resolve exactly once"));
    }
    let edges = mesh.patch_edges(patch).ok_or(VectorBoundaryError2::Invalid("unknown patch"))?;
    if edges.is_empty() { return Err(VectorBoundaryError2::Invalid("boundary patch contains no edges")); }
    for &(a, b) in edges {
        if a >= mesh.node_count() || b >= mesh.node_count() || a == b {
            return Err(VectorBoundaryError2::Invalid("invalid patch-edge node"));
        }
        if !mesh.nodes[a].iter().chain(&mesh.nodes[b]).all(|v| v.is_finite())
            || mesh.nodes[a] == mesh.nodes[b] {
            return Err(VectorBoundaryError2::Invalid("patch edge is nonfinite or has zero length"));
        }
    }
    Ok(mesh.patch_nodes(patch))
}
fn prescribe(fixed: &mut BTreeMap<usize, [f64; 2]>, node: usize, value: [f64; 2])
    -> Result<(), VectorBoundaryError2> {
    if !value.iter().all(|v| v.is_finite()) {
        return Err(VectorBoundaryError2::Invalid("nonfinite prescribed displacement"));
    }
    if let Some(previous) = fixed.insert(node, value) {
        if previous != value { return Err(VectorBoundaryError2::ConflictingDisplacement { node }); }
    }
    Ok(())
}

impl LinearProblem<'_> {
    /// Solve equilibrium with additional spatially uniform vector histories at
    /// `time_s`. Existing body force, static boundary fields, plane reduction,
    /// formulation and symmetry conditions remain in the same physical solve.
    /// Strong prescribed values use the original symmetric RHS lift; traction
    /// uses the original edge Gauss rule. There is no duplicate stiffness code.
    ///
    /// Duplicate patch roles and incompatible shared Dirichlet vertices refuse,
    /// rather than depending on row order. Existing explicit displacement still
    /// takes precedence over symmetry, as in `solve`. Only dimension-contract,
    /// even-parity Cartesian vectors are currently admitted; stronger semantic
    /// contracts refuse rather than being stripped to matching dimensions.
    ///
    /// The supplied mesh/material/body/static fields retain `LinearProblem`'s
    /// requirements (valid P1/Q1 mesh and pure finite callbacks). Added patch
    /// bindings and material scalars are checked here. Neither self nor payloads
    /// are modified. Cancellation is checked while binding and before/after
    /// `solve`; its existing assembly and CG remain indivisible in this adapter.
    /// A late cancellation returns no field but does not pretend work was free.
    pub fn solve_with_vector_boundaries(
        &self, frames: &FrameTree, boundaries: &[VectorBoundary2<'_>], time_s: f64,
        checkpoint: &mut impl FnMut() -> ControlFlow<()>,
    ) -> Result<Vec<[f64; 2]>, VectorBoundaryError2> {
        use VectorBoundaryError2 as E;
        poll(checkpoint)?;
        if !time_s.is_finite() || boundaries.len() > 1024 {
            return Err(E::Invalid("finite endpoint time and at most 1024 vector bindings required"));
        }
        if !self.youngs.is_finite() || self.youngs <= 0.0 || !self.poisson.is_finite()
            || self.poisson <= -1.0 || self.poisson >= 0.5 {
            return Err(E::Invalid("positive finite E and -1 < nu < 0.5 required"));
        }
        let (lambda, mu) = crate::linear::lame(self.youngs, self.poisson, self.plane);
        if !lambda.is_finite() || !mu.is_finite() || mu <= 0.0 {
            return Err(E::Invalid("unrepresentable elastic coefficients"));
        }
        let mut occupied = BTreeSet::new();
        let mut fixed = BTreeMap::new();
        for (patch, field) in &self.dirichlet {
            poll(checkpoint)?;
            if !occupied.insert(*patch) { return Err(E::Invalid("duplicate boundary patch assignment")); }
            for node in patch_nodes(self.mesh, *patch)? {
                poll(checkpoint)?;
                let p = self.mesh.nodes[node];
                prescribe(&mut fixed, node, field(p[0], p[1]))?;
            }
        }
        for (patch, _) in &self.traction {
            poll(checkpoint)?;
            if !occupied.insert(*patch) { return Err(E::Invalid("duplicate or mixed displacement/traction patch assignment")); }
            patch_nodes(self.mesh, *patch)?;
        }
        for (patch, component) in &self.symmetry {
            poll(checkpoint)?;
            if *component >= 2 { return Err(E::Invalid("symmetry component must be x or y")); }
            patch_nodes(self.mesh, *patch)?;
        }
        let mut prescribed = Vec::new();
        let mut tractions = Vec::new();
        for row in boundaries {
            poll(checkpoint)?;
            let (patch, value, expected, displacement) = match *row {
                VectorBoundary2::Displacement { patch, value } => (patch, value, DISPLACEMENT, true),
                VectorBoundary2::Traction { patch, value } => (patch, value, TRACTION, false),
            };
            if !occupied.insert(patch) { return Err(E::Invalid("duplicate or mixed displacement/traction patch assignment")); }
            let nodes = patch_nodes(self.mesh, patch)?;
            if value.meta().contract() != QuantityContract::Dimensions(expected) {
                return Err(E::Quantity { patch, expected, found: value.meta().contract() });
            }
            if value.meta().orientation() != OrientationParity::Even {
                return Err(E::Invalid("displacement and traction require ordinary even-parity vectors"));
            }
            let sample = value.world_vector_at(frames, time_s, checkpoint).map_err(|error| {
                if error == VectorEvaluationError::Cancelled { E::Cancelled }
                else { E::Vector { patch, error } }
            })?;
            if sample.values[2] != 0.0 { return Err(E::OutOfPlane { patch, value: sample.values[2] }); }
            let vector = [sample.values[0], sample.values[1]];
            if displacement {
                for node in nodes { poll(checkpoint)?; prescribe(&mut fixed, node, vector)?; }
                prescribed.push((patch, vector));
            } else { tractions.push((patch, vector)); }
        }
        // Keep the sampled values alive for the complete native invocation.
        let displacements: Vec<_> = prescribed.iter().map(|(_, value)| {
            let value = *value;
            move |_: f64, _: f64| value
        }).collect();
        let loads: Vec<_> = tractions.iter().map(|(_, value)| {
            let value = *value;
            move |_: f64, _: f64| value
        }).collect();
        let mut bound = LinearProblem {
            mesh: self.mesh, youngs: self.youngs, poisson: self.poisson,
            plane: self.plane, formulation: self.formulation, body_force: self.body_force,
            dirichlet: self.dirichlet.clone(), traction: self.traction.clone(),
            symmetry: self.symmetry.clone(),
        };
        for ((patch, _), field) in prescribed.iter().zip(&displacements) { bound.dirichlet.push((*patch, field)); }
        for ((patch, _), field) in tractions.iter().zip(&loads) { bound.traction.push((*patch, field)); }
        poll(checkpoint)?;
        let solution = bound.solve().map_err(E::Solver)?;
        poll(checkpoint)?;
        if !solution.iter().flatten().all(|v| v.is_finite()) { return Err(E::Invalid("nonfinite equilibrium field")); }
        Ok(solution)
    }
}
