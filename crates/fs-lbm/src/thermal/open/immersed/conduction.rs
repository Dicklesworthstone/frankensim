//! Implicit solid conduction and resolved fluid exchange on immersed nodes.
//!
//! A node is an original `ImmersedThermalBody3`, now optionally one spatial
//! solid cell. Zero-interface interior nodes retain their own capacity and
//! source. Authored reciprocal conductances connect these original owners;
//! temperature is never stored in a second authoritative field.
//!
//! With frozen post-collision fluid traces, the simultaneous endpoint satisfies
//! `(C + K_fluid + L_G) delta_T = Q + K_fluid(T_trace-T_old) - L_G T_old`.
//! The native sparse LDLT implementation also used by `WallConduction` factors
//! the diagonally scaled SPD operator once. Subsequent steps reuse that factor
//! and linear workspace. Body storage, actual reflected heat and reciprocal
//! endpoint edge heat are audited by the original thermal transaction.

use super::{
    Body, ImmersedThermalBody3, Sum, ThermalOpenError3, ThermalOpenField3, ThermalOpenLimits3,
    check_boundary_temperature, invalid,
};
use crate::thermal::boundary::conduction::{WallConductionControls, WallConductionError};
use fs_sparse::{Coo, DirectOrdering, LdltError, LdltFactor, LdltOptions, SymbolicLdlt};

/// An authored undirected solid path, using node declaration indices.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ImmersedConductionLink3 {
    pub a: usize,
    pub b: usize,
    /// Nonnegative G_SI*dt/C_f. Zero explicitly insulates. Geometry, material
    /// conductivity and any contact resistance determine this value externally.
    /// Combine parallel paths before installation; repeated pairs are refused.
    pub conductance: f64,
}

pub(super) struct ImmersedConduction3 {
    edges: Vec<ImmersedConductionLink3>,
    scale: Vec<f64>,
    factor: LdltFactor,
    rhs: Vec<f64>,
    tolerance: f64,
}

fn sparse(error: LdltError) -> ThermalOpenError3 {
    ThermalOpenError3::Conduction(WallConductionError::Sparse(error))
}

fn buffer<T: Clone + Default>(len: usize) -> Result<Vec<T>, ThermalOpenError3> {
    let mut data = Vec::new();
    data.try_reserve_exact(len)
        .map_err(|_| ThermalOpenError3::Allocation)?;
    data.resize(len, T::default());
    Ok(data)
}

fn root(parents: &mut [usize], mut node: usize) -> usize {
    while parents[node] != node {
        parents[node] = parents[parents[node]];
        node = parents[node];
    }
    node
}

fn difference(coefficient: f64, first: f64, second: f64) -> f64 {
    let delta = first - second;
    if delta.is_finite() {
        coefficient * delta
    } else {
        coefficient * first - coefficient * second
    }
}

impl ThermalOpenField3 {
    /// Install original immersed nodes and their reciprocal solid graph
    /// atomically. An interior node may have zero fluid links, but its positive
    /// conductance component must reach an exposed node. Interior heat is then
    /// accounted in the same native fluid/solid transaction.
    ///
    /// `maximum_storage_bytes` bounds original node/interface storage and its
    /// preparation scratch as for `set_immersed_thermal_bodies`. The additional
    /// `max_factor_bytes` bounds the retained graph factor and its linear work
    /// arrays. `max_nodes` covers ALL nodes; `max_links` includes zero edges.
    /// The smaller graph/carrier tolerance gates local and complete balances.
    /// Replacing a network briefly retains the old plan until success and
    /// explicitly installs the supplied temperatures and inventories.
    ///
    /// All capacities, resistances and conductances are fixed until replaced.
    /// Held sources may change through `set_immersed_body_heat`. Time accuracy
    /// is first order for the solid solve. No material/geometry law, moving
    /// boundary, nonlinear caloric response or spatial convergence is inferred.
    pub fn set_immersed_thermal_network(
        &mut self,
        declarations: &[ImmersedThermalBody3<'_>],
        links: &[ImmersedConductionLink3],
        maximum_storage_bytes: usize,
        controls: WallConductionControls,
    ) -> Result<(), ThermalOpenError3> {
        if declarations.len() > controls.max_nodes || links.len() > controls.max_links {
            return Err(invalid("immersed conduction node or edge budget exceeded"));
        }
        let mut prepared =
            self.prepare_immersed_bodies(declarations, maximum_storage_bytes, true)?;
        prepared.conduction = Some(ImmersedConduction3::try_new(
            &prepared.bodies,
            links,
            controls,
        )?);
        self.immersed = Some(prepared);
        Ok(())
    }
}

impl ImmersedConduction3 {
    pub(super) fn try_new(
        bodies: &[Body],
        links: &[ImmersedConductionLink3],
        controls: WallConductionControls,
    ) -> Result<Self, ThermalOpenError3> {
        let n = bodies.len();
        if n == 0
            || n > controls.max_nodes
            || links.len() > controls.max_links
            || controls.max_nodes == 0
            || controls.max_links == 0
            || !controls.balance_tolerance.is_finite()
            || controls.balance_tolerance < 0.0
        {
            return Err(invalid("invalid immersed conduction controls"));
        }
        let retained_bytes = n
            .checked_mul(2 * core::mem::size_of::<f64>())
            .and_then(|bytes| {
                links
                    .len()
                    .checked_mul(core::mem::size_of::<ImmersedConductionLink3>())
                    .and_then(|edges| bytes.checked_add(edges))
            })
            .ok_or_else(|| invalid("immersed conduction storage overflow"))?;
        let factor_budget = controls
            .max_factor_bytes
            .checked_sub(retained_bytes)
            .filter(|bytes| *bytes > 0)
            .ok_or_else(|| invalid("immersed conduction factor/workspace budget exceeded"))?;
        let mut edges = Vec::new();
        edges
            .try_reserve_exact(links.len())
            .map_err(|_| ThermalOpenError3::Allocation)?;
        for link in links {
            if link.a >= n
                || link.b >= n
                || link.a == link.b
                || !link.conductance.is_finite()
                || link.conductance < 0.0
            {
                return Err(invalid(
                    "invalid immersed conduction endpoint or conductance",
                ));
            }
            edges.push(ImmersedConductionLink3 {
                a: link.a.min(link.b),
                b: link.a.max(link.b),
                conductance: link.conductance,
            });
        }
        edges.sort_unstable_by_key(|edge| (edge.a, edge.b));
        if edges
            .windows(2)
            .any(|pair| (pair[0].a, pair[0].b) == (pair[1].a, pair[1].b))
        {
            return Err(invalid("duplicate undirected immersed conduction edge"));
        }
        edges.retain(|edge| edge.conductance > 0.0);
        let mut parents: Vec<_> = buffer(n)?;
        for (i, parent) in parents.iter_mut().enumerate() {
            *parent = i;
        }
        for edge in &edges {
            let a = root(&mut parents, edge.a);
            let b = root(&mut parents, edge.b);
            parents[a.max(b)] = a.min(b);
        }
        let mut exposed = buffer::<bool>(n)?;
        for (i, body) in bodies.iter().enumerate() {
            if body.accepted.interface_links != 0 {
                exposed[root(&mut parents, i)] = true;
            }
        }
        for i in 0..n {
            if !exposed[root(&mut parents, i)] {
                return Err(invalid(
                    "each immersed conduction component must reach a fluid interface",
                ));
            }
        }
        // Sum positive diagonals in a local scale, never as C + K + sum G.
        // This admits finite normalized matrices even when that raw sum would
        // overflow. No singular/ill-conditioned factor is silently repaired.
        let mut maximum = buffer::<f64>(n)?;
        for (i, body) in bodies.iter().enumerate() {
            maximum[i] = body.accepted.solid.heat_capacity.max(body.conductance);
        }
        for edge in &edges {
            maximum[edge.a] = maximum[edge.a].max(edge.conductance);
            maximum[edge.b] = maximum[edge.b].max(edge.conductance);
        }
        let mut diagonal = buffer::<Sum>(n)?;
        for (i, body) in bodies.iter().enumerate() {
            diagonal[i].add(body.accepted.solid.heat_capacity / maximum[i]);
            diagonal[i].add(body.conductance / maximum[i]);
        }
        for edge in &edges {
            diagonal[edge.a].add(edge.conductance / maximum[edge.a]);
            diagonal[edge.b].add(edge.conductance / maximum[edge.b]);
        }
        let mut scale = buffer::<f64>(n)?;
        let mut matrix = Coo::new(n, n);
        for i in 0..n {
            scale[i] = (1.0 / maximum[i].sqrt()) / diagonal[i].total().sqrt();
            if !scale[i].is_finite() || scale[i] <= 0.0 {
                return Err(invalid("unrepresentable immersed conduction scaling"));
            }
            matrix.push(i, i, 1.0);
        }
        for edge in &edges {
            let sqrt_g = edge.conductance.sqrt();
            let value = -(sqrt_g * scale[edge.a]) * (sqrt_g * scale[edge.b]);
            matrix.push(edge.a, edge.b, value);
            matrix.push(edge.b, edge.a, value);
        }
        let matrix = matrix.assemble();
        let symbolic = SymbolicLdlt::analyze(&matrix, DirectOrdering::Amd).map_err(sparse)?;
        let factor = symbolic
            .factor(
                &matrix,
                &LdltOptions {
                    max_factor_bytes: Some(factor_budget),
                    ..LdltOptions::default()
                },
            )
            .map_err(sparse)?;
        if factor.inertia().positive != n || factor.inertia().negative != 0 {
            return Err(invalid(
                "immersed conduction operator lost positive definiteness",
            ));
        }
        Ok(Self {
            edges,
            scale,
            factor,
            rhs: buffer(n)?,
            tolerance: controls.balance_tolerance,
        })
    }

    pub(super) fn balance_tolerance(&self) -> f64 {
        self.tolerance
    }

    pub(super) fn solve<E: From<ThermalOpenError3>>(
        &mut self,
        bodies: &mut [Body],
        limits: ThermalOpenLimits3,
        mut poll: impl FnMut() -> Result<(), E>,
    ) -> Result<(), E> {
        // The original body's compensated scratch first accumulates the
        // scaled increment RHS. Endpoint edge heat replaces it after solving.
        for (i, body) in bodies.iter_mut().enumerate() {
            if i % 64 == 0 {
                poll()?;
            }
            body.conducted = Sum::default();
            body.conducted
                .add(self.scale[i] * body.accepted.solid.heat_per_step);
            if body.accepted.interface_links != 0 {
                let k = body.conductance;
                body.conducted.add(difference(
                    k * self.scale[i],
                    body.mean_fluid_trace(),
                    body.accepted.solid.temperature,
                ));
            }
        }
        for (i, edge) in self.edges.iter().enumerate() {
            if i % 64 == 0 {
                poll()?;
            }
            let a = bodies[edge.a].accepted.solid.temperature;
            let b = bodies[edge.b].accepted.solid.temperature;
            bodies[edge.a]
                .conducted
                .add(difference(edge.conductance * self.scale[edge.a], b, a));
            bodies[edge.b]
                .conducted
                .add(difference(edge.conductance * self.scale[edge.b], a, b));
        }
        for (i, body) in bodies.iter().enumerate() {
            if i % 64 == 0 {
                poll()?;
            }
            self.rhs[i] = body.conducted.total();
            if !self.rhs[i].is_finite() {
                return Err(invalid("unrepresentable immersed conduction right-hand side").into());
            }
        }
        poll()?;
        self.factor.solve_in_place(&mut self.rhs);
        poll()?;
        for (i, body) in bodies.iter_mut().enumerate() {
            if i % 64 == 0 {
                poll()?;
            }
            body.trial_temperature = body.accepted.solid.temperature + self.scale[i] * self.rhs[i];
            check_boundary_temperature(6, body.trial_temperature, limits)?;
            body.conducted = Sum::default();
        }
        for (i, edge) in self.edges.iter().enumerate() {
            if i % 64 == 0 {
                poll()?;
            }
            let heat = difference(
                edge.conductance,
                bodies[edge.b].trial_temperature,
                bodies[edge.a].trial_temperature,
            );
            if !heat.is_finite() {
                return Err(invalid("unrepresentable immersed conductive heat").into());
            }
            bodies[edge.a].conducted.add(heat);
            bodies[edge.b].conducted.add(-heat);
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests;
