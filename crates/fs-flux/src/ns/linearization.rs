//! Sparse derivatives of the actual generalized-Newtonian discrete residual.
//!
//! The frozen operator alone omits `(dA/du)u - db/du`. This module adds those
//! terms locally, including viscosity-dependent SIP weights and boundary loads,
//! and (when requested) the velocity transporting the DG convective flux.

use std::cell::RefCell;

use super::{FluxParams, FluxSystem, viscous_normal};
use crate::bdm::{edge_gauss_pub, eval_basis, tri_quad};
use fs_sparse::{Coo, Csr};

/// A constitutive value and its derivative at the same physical shear rate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FluxViscosityTangent {
    /// Positive finite kinematic viscosity [m²/s], at fixed reference density.
    pub nu: f64,
    /// `d nu / d gamma`, holding position and external material fields fixed [m²].
    pub derivative: f64,
}

/// Residual and its sparse velocity/pressure linearization at one exact state.
pub struct FluxRheologyLinearization {
    /// `dR/dx`, with the same strong normal-velocity and pressure-datum rows.
    pub jacobian: Csr,
    /// `R(x) = A(x) x - b(x)`, including the endpoint constitutive field.
    pub residual: Vec<f64>,
    /// Quadrature sites with zero advecting normal flux and unequal upwind
    /// traces. There the returned `>= 0` active-branch derivative is semismooth,
    /// not a unique classical derivative. A smooth adjoint claim must exclude
    /// these switches or establish that admissible variations cannot cross them.
    pub nonsmooth_upwind_points: usize,
}

/// A refused linearization publishes neither a residual nor a Jacobian.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum FluxLinearizationError {
    /// Invalid state, timestep, coefficient, geometry-derived arithmetic or load.
    Invalid(&'static str),
    /// Invalid constitutive response; a smooth radial viscosity needs zero
    /// shear-rate derivative at zero strain.
    Constitutive {
        /// Triangle index.
        cell: usize,
        /// Physical quadrature position [m].
        point: [f64; 2],
        /// Objective strain rate [1/s].
        shear_rate: f64,
        /// Rejected value and derivative.
        tangent: FluxViscosityTangent,
    },
}

impl core::fmt::Display for FluxLinearizationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(reason) => write!(f, "flow linearization: {reason}"),
            Self::Constitutive {
                cell,
                point,
                shear_rate,
                tangent,
            } => write!(
                f,
                "flow linearization: invalid constitutive tangent {tangent:?} in cell {cell} at {point:?}, rate {shear_rate}"
            ),
        }
    }
}
impl std::error::Error for FluxLinearizationError {}

fn sample(
    law: &dyn Fn(usize, [f64; 2], f64) -> FluxViscosityTangent,
    cell: usize,
    point: [f64; 2],
    rate: f64,
) -> Result<FluxViscosityTangent, FluxLinearizationError> {
    let tangent = law(cell, point, rate);
    if !tangent.nu.is_finite()
        || tangent.nu <= 0.0
        || !tangent.derivative.is_finite()
        || (rate == 0.0 && tangent.derivative != 0.0)
    {
        return Err(FluxLinearizationError::Constitutive {
            cell,
            point,
            shear_rate: rate,
            tangent,
        });
    }
    Ok(tangent)
}

fn dot(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}

impl FluxSystem<'_> {
    /// Linearize steady or BDF1 generalized-Newtonian flow with velocity data
    /// on every boundary edge. `state` contains ALL velocity and pressure dofs.
    /// `inverse_timestep == 0` selects steady flow; a positive value requires
    /// exactly `n_u` previous velocity moments. Previous state, force, boundary
    /// data, reference density and captured temperature/composition stay fixed.
    /// `params.nu` is unused; `params.sigma` sets the actual SIP penalty.
    ///
    /// The pure constitutive callback supplies `nu(gamma)` and its analytic
    /// derivative at the same volume/one-sided face quadrature point, where
    /// `gamma = sqrt(2 D:D)`. It may be queried more than once per point. Smooth
    /// regularized laws must supply derivative zero at gamma zero; no numerical
    /// derivative, viscosity floor or constitutive fit is introduced here.
    ///
    /// The Jacobian includes cell stress, harmonic interface weights, penalties,
    /// viscosity-dependent velocity-boundary loads, and both transported and
    /// advecting velocities when `convection` is true. Mass, incompressibility,
    /// normal constraints and pressure pin are inherited from the SAME primal
    /// assembler. Extra derivatives are local sparse updates, O(cells + faces);
    /// no dense Jacobian or perturb-and-reassemble column loop is used.
    ///
    /// Feed `jacobian` to `solve_adjoint_checked` for the discrete implicit
    /// derivative: `dJ/dtheta = J_theta - lambda^T R_theta`. Parameter-dependent
    /// loads and external fields still belong in `R_theta`. At upwind switches
    /// the active `>= 0` branch is a semismooth linearization; inspect
    /// `nonsmooth_upwind_points` before claiming a classical gradient.
    /// This is one fixed-mesh state derivative, not a time-history adjoint.
    ///
    /// # Errors
    /// Refuses invalid dimensions, nonfinite data, nonpositive viscosity or
    /// penalty, invalid constitutive tangents and unrepresentable arithmetic.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn linearize_rheology(
        &self,
        params: FluxParams,
        law: &dyn Fn(usize, [f64; 2], f64) -> FluxViscosityTangent,
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        boundary: &dyn Fn([f64; 2]) -> [f64; 2],
        state: &[f64],
        convection: bool,
        inverse_timestep: f64,
        previous_velocity: Option<&[f64]>,
    ) -> Result<FluxRheologyLinearization, FluxLinearizationError> {
        if self.mesh.tris.is_empty()
            || state.len() != self.n
            || state.iter().any(|v| !v.is_finite())
            || !params.sigma.is_finite()
            || params.sigma <= 0.0
            || !inverse_timestep.is_finite()
            || inverse_timestep < 0.0
            || (inverse_timestep > 0.0 && previous_velocity.is_none())
            || previous_velocity
                .is_some_and(|u| u.len() != self.n_u || u.iter().any(|v| !v.is_finite()))
        {
            return Err(FluxLinearizationError::Invalid(
                "finite full state, positive penalty and valid steady/BDF1 history required",
            ));
        }
        let mesh = self.mesh;
        let rates = self.shear_rates(state);
        if rates.iter().any(|rate| !rate.is_finite()) {
            return Err(FluxLinearizationError::Invalid("nonfinite strain rate"));
        }
        let failure = RefCell::new(None);
        let viscosity = |cell, point| {
            if failure.borrow().is_some() {
                return 1.0;
            }
            match sample(law, cell, point, rates[cell]) {
                Ok(value) => value.nu,
                Err(error) => {
                    *failure.borrow_mut() = Some(error);
                    1.0 // private assembly is discarded before any publication
                }
            }
        };
        let checked_vector = |field: &dyn Fn([f64; 2]) -> [f64; 2], point| {
            if failure.borrow().is_some() {
                return [0.0; 2];
            }
            let value = field(point);
            if value.iter().any(|v| !v.is_finite()) {
                *failure.borrow_mut() = Some(FluxLinearizationError::Invalid(
                    "nonfinite load or velocity data",
                ));
                [0.0; 2]
            } else {
                value
            }
        };
        let (frozen, rhs) = self.assemble_viscous_transient(
            params,
            &viscosity,
            &|point| checked_vector(force, point),
            &|point| checked_vector(boundary, point),
            convection.then_some(&state[..self.n_u]),
            inverse_timestep,
            previous_velocity,
        );
        if let Some(error) = failure.into_inner() {
            return Err(error);
        }
        let mut residual = vec![0.0; self.n];
        frozen.spmv(state, &mut residual);
        for (value, load) in residual.iter_mut().zip(rhs) {
            *value -= load;
        }
        let mut jacobian = Coo::new(self.n, self.n);
        for row in 0..self.n {
            let (columns, values) = frozen.row(row);
            for (&column, &value) in columns.iter().zip(values) {
                jacobian.push(row, column, value);
            }
        }
        let mut constrained = vec![false; self.n_u];
        for (edge, geometry) in mesh.edges.iter().enumerate() {
            if geometry.tris.1 == usize::MAX {
                constrained[2 * edge] = true;
                constrained[2 * edge + 1] = true;
            }
        }
        let gradients: Vec<[[f64; 2]; 2]> = (0..mesh.tris.len())
            .map(|cell| {
                let mut gradient = [[0.0; 2]; 2];
                for i in 0..6 {
                    for (r, row) in gradient.iter_mut().enumerate() {
                        for (c, value) in row.iter_mut().enumerate() {
                            *value += state[self.dof(cell, i)] * self.bases[cell].grad[i][r][c];
                        }
                    }
                }
                gradient
            })
            .collect();
        let rate_derivatives: Vec<[f64; 6]> = gradients
            .iter()
            .enumerate()
            .map(|(cell, g)| {
                core::array::from_fn(|j| {
                    if rates[cell] == 0.0 {
                        return 0.0;
                    }
                    let d = self.bases[cell].grad[j];
                    2.0 * (g[0][0] / rates[cell]) * d[0][0]
                        + 2.0 * (g[1][1] / rates[cell]) * d[1][1]
                        + ((g[0][1] + g[1][0]) / rates[cell]) * (d[0][1] + d[1][0])
                })
            })
            .collect();
        for cell in 0..mesh.tris.len() {
            let triangle = core::array::from_fn(|k| mesh.verts[mesh.tris[cell][k]]);
            let mut integrated_derivative = 0.0;
            for (point, weight) in tri_quad(triangle, mesh.areas[cell]) {
                integrated_derivative += weight * sample(law, cell, point, rates[cell])?.derivative;
            }
            if integrated_derivative == 0.0 {
                continue;
            }
            let g = gradients[cell];
            for i in 0..6 {
                let row = self.dof(cell, i);
                if constrained[row] {
                    continue;
                }
                let v = self.bases[cell].grad[i];
                let stress = 2.0 * g[0][0] * v[0][0]
                    + 2.0 * g[1][1] * v[1][1]
                    + (g[0][1] + g[1][0]) * (v[0][1] + v[1][0]);
                for j in 0..6 {
                    jacobian.push(
                        row,
                        self.dof(cell, j),
                        integrated_derivative * stress * rate_derivatives[cell][j],
                    );
                }
            }
        }
        for edge in &mesh.edges {
            let (a, b) = (mesh.verts[edge.verts.0], mesh.verts[edge.verts.1]);
            let (c0, c1) = edge.tris;
            let h0 = 2.0 * mesh.areas[c0] / edge.len;
            let tau0 = viscous_normal(gradients[c0], edge.normal, true);
            for (point, weight) in edge_gauss_pub(a, b) {
                let first = sample(law, c0, point, rates[c0])?;
                let u0 = self.linearization_velocity(state, c0, point);
                if c1 == usize::MAX {
                    if first.derivative == 0.0 {
                        continue;
                    }
                    let prescribed = boundary(point);
                    if prescribed.iter().any(|v| !v.is_finite()) {
                        return Err(FluxLinearizationError::Invalid(
                            "nonfinite boundary velocity",
                        ));
                    }
                    let jump = [u0[0] - prescribed[0], u0[1] - prescribed[1]];
                    for i in 0..6 {
                        let row = self.dof(c0, i);
                        if constrained[row] {
                            continue;
                        }
                        let v = eval_basis(&self.bases[c0], i, point);
                        let tv = viscous_normal(self.bases[c0].grad[i], edge.normal, true);
                        // Differentiate A(nu)u - b(nu), including BOTH weak
                        // boundary-load terms, with prescribed velocity fixed.
                        let factor =
                            -dot(tau0, v) - dot(tv, jump) + params.sigma / h0 * dot(v, jump);
                        for j in 0..6 {
                            jacobian.push(
                                row,
                                self.dof(c0, j),
                                weight * factor * first.derivative * rate_derivatives[c0][j],
                            );
                        }
                    }
                    continue;
                }
                let second = sample(law, c1, point, rates[c1])?;
                if first.derivative == 0.0 && second.derivative == 0.0 {
                    continue;
                }
                let h1 = 2.0 * mesh.areas[c1] / edge.len;
                let delta = [first.nu / h0, second.nu / h1];
                let scale = delta[0].max(delta[1]);
                let d = delta.map(|x| x / scale);
                let weights = [d[1] / (d[0] + d[1]), d[0] / (d[0] + d[1])];
                let cells = [c0, c1];
                let heights = [h0, h1];
                let derivatives = [first.derivative, second.derivative];
                let tau1 = viscous_normal(gradients[c1], edge.normal, true);
                let traction = [h0 * tau0[0] + h1 * tau1[0], h0 * tau0[1] + h1 * tau1[1]];
                let u1 = self.linearization_velocity(state, c1, point);
                let jump = [u0[0] - u1[0], u0[1] - u1[1]];
                // q=delta0*delta1/(delta0+delta1), traction_s=h_s*q,
                // penalty=2*sigma*q; dq/dnu_s=weight_s^2/h_s. This
                // differentiates the harmonic weighting without products of
                // the possibly very different material viscosities.
                for (side, &cell) in cells.iter().enumerate() {
                    let sign = if side == 0 { 1.0 } else { -1.0 };
                    for i in 0..6 {
                        let row = self.dof(cell, i);
                        if constrained[row] {
                            continue;
                        }
                        let v = eval_basis(&self.bases[cell], i, point);
                        let tv = viscous_normal(self.bases[cell].grad[i], edge.normal, true);
                        let factor = -sign * dot(traction, v) - heights[side] * dot(tv, jump)
                            + 2.0 * params.sigma * sign * dot(v, jump);
                        for (source, &other) in cells.iter().enumerate() {
                            if derivatives[source] == 0.0 {
                                continue;
                            }
                            let dq = weights[source] * weights[source] / heights[source]
                                * derivatives[source];
                            for j in 0..6 {
                                jacobian.push(
                                    row,
                                    self.dof(other, j),
                                    weight * factor * dq * rate_derivatives[other][j],
                                );
                            }
                        }
                    }
                }
            }
        }
        let nonsmooth_upwind_points = if convection {
            self.add_convecting_derivative(&mut jacobian, state, boundary, &constrained)?
        } else {
            0
        };
        let jacobian = jacobian.assemble();
        if residual.iter().any(|v| !v.is_finite())
            || (0..self.n).any(|row| jacobian.row(row).1.iter().any(|v| !v.is_finite()))
        {
            return Err(FluxLinearizationError::Invalid(
                "nonfinite assembled residual or Jacobian",
            ));
        }
        Ok(FluxRheologyLinearization {
            jacobian,
            residual,
            nonsmooth_upwind_points,
        })
    }

    fn linearization_velocity(&self, state: &[f64], cell: usize, point: [f64; 2]) -> [f64; 2] {
        let mut velocity = [0.0; 2];
        for i in 0..6 {
            let basis = eval_basis(&self.bases[cell], i, point);
            velocity[0] += state[self.dof(cell, i)] * basis[0];
            velocity[1] += state[self.dof(cell, i)] * basis[1];
        }
        velocity
    }

    fn add_convecting_derivative(
        &self,
        jacobian: &mut Coo,
        state: &[f64],
        boundary: &dyn Fn([f64; 2]) -> [f64; 2],
        constrained: &[bool],
    ) -> Result<usize, FluxLinearizationError> {
        let mesh = self.mesh;
        for cell in 0..mesh.tris.len() {
            let triangle = core::array::from_fn(|k| mesh.verts[mesh.tris[cell][k]]);
            for (point, weight) in tri_quad(triangle, mesh.areas[cell]) {
                let u = self.linearization_velocity(state, cell, point);
                for i in 0..6 {
                    let row = self.dof(cell, i);
                    if constrained[row] {
                        continue;
                    }
                    let gradient = self.bases[cell].grad[i];
                    for j in 0..6 {
                        let v = eval_basis(&self.bases[cell], j, point);
                        jacobian.push(
                            row,
                            self.dof(cell, j),
                            -weight * (u[0] * dot(gradient[0], v) + u[1] * dot(gradient[1], v)),
                        );
                    }
                }
            }
        }
        let mut switches = 0;
        for edge in &mesh.edges {
            let (a, b) = (mesh.verts[edge.verts.0], mesh.verts[edge.verts.1]);
            let (c0, c1) = edge.tris;
            for (point, weight) in edge_gauss_pub(a, b) {
                let u0 = self.linearization_velocity(state, c0, point);
                let normal = dot(u0, edge.normal);
                let other = if c1 == usize::MAX {
                    boundary(point)
                } else {
                    self.linearization_velocity(state, c1, point)
                };
                if other.iter().any(|v| !v.is_finite()) {
                    return Err(FluxLinearizationError::Invalid(
                        "nonfinite convective trace",
                    ));
                }
                if normal == 0.0 && u0 != other {
                    switches += 1;
                }
                let upwind = if normal >= 0.0 { u0 } else { other };
                let sides = if c1 == usize::MAX { 1 } else { 2 };
                for side in 0..sides {
                    let cell = if side == 0 { c0 } else { c1 };
                    let sign = if side == 0 { 1.0 } else { -1.0 };
                    for i in 0..6 {
                        let row = self.dof(cell, i);
                        if constrained[row] {
                            continue;
                        }
                        let test = eval_basis(&self.bases[cell], i, point);
                        for j in 0..6 {
                            let trial = eval_basis(&self.bases[c0], j, point);
                            jacobian.push(
                                row,
                                self.dof(c0, j),
                                weight * sign * dot(upwind, test) * dot(trial, edge.normal),
                            );
                        }
                    }
                }
            }
        }
        Ok(switches)
    }
}
