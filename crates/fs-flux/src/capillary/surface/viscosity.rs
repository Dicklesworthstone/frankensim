//! Boussinesq--Scriven shear and dilatational viscosity on material triangles.
//!
//! With the tangential projector P and surface rate D, the stress is
//! `tau = 2*mu_s*(D - tr(D)*P/2) + kappa_s*tr(D)*P`. Weak nodal forces are
//! `f_i = -A*tau*grad_s(N_i)`; their power is minus
//! `A*(2*mu_s*D_dev:D_dev + kappa_s*tr(D)^2)`. Both coefficients have units
//! N s/m and must be nonnegative. They are interfacial transport properties,
//! not the bulk dynamic viscosity multiplied by an arbitrary mesh thickness.
//!
//! This is the Newtonian, isotropic surface-fluid constitutive law of Scriven,
//! Chemical Engineering Science 12 (1960), 98--108,
//! <https://doi.org/10.1016/0009-2509(60)87003-0>. Coefficients are prescribed;
//! no surfactant transport, surface elasticity, bending viscosity, bulk fluid,
//! or thermal conversion is inferred. A two-sided film needs both interfaces.
//! Leaves allocate nothing, and the analytic directional derivative includes
//! the changing area, projector and shape gradients. Patch callers own mesh
//! validity, topology, work caps and cancellation.

use super::{CapillarySurfaceError, Sum, Vector, cross, dot, geometry, scaled, sub};

type Matrix = [[f64; 3]; 3];

/// Nonnegative Newtonian interfacial transport coefficients [N s/m].
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct SurfaceViscosity {
    /// Resistance to tangential shear at fixed area [N s/m].
    pub shear_n_s_m: f64,
    /// Resistance to local area dilation [N s/m].
    pub dilatational_n_s_m: f64,
}

impl SurfaceViscosity {
    /// An inviscid interface; recovers the original capillary surface law.
    pub const ZERO: Self = Self {
        shear_n_s_m: 0.0,
        dilatational_n_s_m: 0.0,
    };

    fn validate(self) -> Result<(), CapillarySurfaceError> {
        if [self.shear_n_s_m, self.dilatational_n_s_m]
            .into_iter()
            .any(|x| !x.is_finite() || x < 0.0)
        {
            Err(CapillarySurfaceError::InvalidSurfaceViscosity)
        } else {
            Ok(())
        }
    }
}

/// Constant or explicitly prescribed per-face surface viscosities.
#[derive(Clone, Copy, Debug)]
pub enum SurfaceViscosityField<'a> {
    /// The same pair of coefficients on every face.
    Uniform(SurfaceViscosity),
    /// One pair per triangle in the supplied connectivity order.
    PerFace(&'a [SurfaceViscosity]),
}

impl SurfaceViscosityField<'_> {
    pub(super) fn is_zero(self) -> bool {
        match self {
            Self::Uniform(value) => value == SurfaceViscosity::ZERO,
            Self::PerFace(values) => values.iter().all(|&value| value == SurfaceViscosity::ZERO),
        }
    }

    pub(super) fn validate(self, faces: usize) -> Result<(), CapillarySurfaceError> {
        match self {
            Self::Uniform(value) => value.validate(),
            Self::PerFace(values) if values.len() == faces => {
                values.iter().try_for_each(|value| value.validate())
            }
            Self::PerFace(_) => Err(CapillarySurfaceError::LengthMismatch),
        }
    }

    pub(super) fn at(self, face: usize) -> SurfaceViscosity {
        match self {
            Self::Uniform(value) => value,
            Self::PerFace(values) => values[face],
        }
    }
}

/// One triangle's viscous forces and nonnegative mechanical loss rate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TriangleSurfaceViscosityResponse {
    /// Forces on the vertices, in input order [N].
    pub forces_n: [Vector; 3],
    /// Irreversible mechanical power loss; never stored free energy [W].
    pub dissipation_w: f64,
}

/// Assembled interface viscous forces and total mechanical loss rate.
#[derive(Clone, Debug, PartialEq)]
pub struct SurfaceViscosityResponse {
    /// Nodal viscous forces; unused vertices retain zero [N].
    pub forces_n: Vec<Vector>,
    /// Deterministically accumulated nonnegative mechanical loss rate [W].
    pub dissipation_w: f64,
}

struct Kinematics {
    area: f64,
    normal: Vector,
    a: Vector,
    b: Vector,
    gradients: [Vector; 3],
    projector: Matrix,
    velocity_edges: [Vector; 2],
    velocity_gradient: Matrix,
    rate: Matrix,
    trace: f64,
}

fn kinematics(
    positions_m: [Vector; 3],
    velocities_m_s: [Vector; 3],
) -> Result<Kinematics, CapillarySurfaceError> {
    if velocities_m_s.iter().flatten().any(|x| !x.is_finite()) {
        return Err(CapillarySurfaceError::NonFiniteVelocity);
    }
    let g = geometry(positions_m)?;
    let g1 = scaled(cross(g.b, g.normal), 0.5 / g.area);
    let g2 = scaled(cross(g.normal, g.a), 0.5 / g.area);
    let gradients = [core::array::from_fn(|d| -g1[d] - g2[d]), g1, g2];
    let projector = core::array::from_fn(|i| {
        core::array::from_fn(|j| (if i == j { 1.0 } else { 0.0 }) - g.normal[i] * g.normal[j])
    });
    let velocity_edges = [
        sub(velocities_m_s[1], velocities_m_s[0]),
        sub(velocities_m_s[2], velocities_m_s[0]),
    ];
    let velocity_gradient = core::array::from_fn(|i| {
        core::array::from_fn(|j| velocity_edges[0][i] * g1[j] + velocity_edges[1][i] * g2[j])
    });
    let tangent_gradient = product(projector, velocity_gradient);
    let rate: Matrix = core::array::from_fn(|i| {
        core::array::from_fn(|j| 0.5 * tangent_gradient[i][j] + 0.5 * tangent_gradient[j][i])
    });
    let trace = rate[0][0] + rate[1][1] + rate[2][2];
    if gradients
        .iter()
        .flatten()
        .chain(projector.iter().flatten())
        .chain(velocity_gradient.iter().flatten())
        .chain(rate.iter().flatten())
        .any(|x| !x.is_finite())
        || !trace.is_finite()
    {
        return Err(CapillarySurfaceError::Unrepresentable);
    }
    Ok(Kinematics {
        area: g.area,
        normal: g.normal,
        a: g.a,
        b: g.b,
        gradients,
        projector,
        velocity_edges,
        velocity_gradient,
        rate,
        trace,
    })
}

fn stress(k: &Kinematics, viscosity: SurfaceViscosity) -> Matrix {
    core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            let dilation = k.trace * k.projector[i][j];
            2.0 * viscosity.shear_n_s_m * (k.rate[i][j] - 0.5 * dilation)
                + viscosity.dilatational_n_s_m * dilation
        })
    })
}

fn product(a: Matrix, b: Matrix) -> Matrix {
    core::array::from_fn(|i| core::array::from_fn(|j| (0..3).map(|d| a[i][d] * b[d][j]).sum()))
}

fn apply(a: Matrix, x: Vector) -> Vector {
    a.map(|row| dot(row, x))
}

/// Apply objective surface viscous stress to one current material triangle.
///
/// The returned forces balance linear and angular momentum. Their contraction
/// with the supplied velocities is `-dissipation_w`. Rigid translations and
/// rotations have zero strain rate and do not dissipate. Reversing face winding
/// does not change the physical result. Both endpoints of every edge are
/// retained, so open patches include their proper boundary tractions.
pub fn triangle_surface_viscosity(
    positions_m: [Vector; 3],
    velocities_m_s: [Vector; 3],
    viscosity: SurfaceViscosity,
) -> Result<TriangleSurfaceViscosityResponse, CapillarySurfaceError> {
    viscosity.validate()?;
    let k = kinematics(positions_m, velocities_m_s)?;
    let tau = stress(&k, viscosity);
    let f1 = scaled(apply(tau, k.gradients[1]), -k.area);
    let f2 = scaled(apply(tau, k.gradients[2]), -k.area);
    let forces = [core::array::from_fn(|d| -f1[d] - f2[d]), f1, f2];
    let mut shear = Sum::default();
    for i in 0..3 {
        for j in 0..3 {
            let rate = k.rate[i][j] - 0.5 * k.trace * k.projector[i][j];
            shear.add(rate * rate);
        }
    }
    let dissipation = k.area
        * (2.0 * viscosity.shear_n_s_m * shear.value()
            + viscosity.dilatational_n_s_m * k.trace * k.trace);
    if forces.iter().flatten().any(|x| !x.is_finite()) || !dissipation.is_finite() {
        return Err(CapillarySurfaceError::Unrepresentable);
    }
    Ok(TriangleSurfaceViscosityResponse {
        forces_n: forces,
        dissipation_w: dissipation,
    })
}

/// Exact force directional derivative for simultaneous geometry and velocity
/// perturbations. This differentiates the same physical stress as the force
/// leaf, including surface projection and geometry; coefficients are held.
pub fn triangle_surface_viscosity_tangent(
    positions_m: [Vector; 3],
    velocities_m_s: [Vector; 3],
    viscosity: SurfaceViscosity,
    position_direction: [Vector; 3],
    velocity_direction: [Vector; 3],
) -> Result<[Vector; 3], CapillarySurfaceError> {
    viscosity.validate()?;
    if position_direction
        .iter()
        .flatten()
        .chain(velocity_direction.iter().flatten())
        .any(|x| !x.is_finite())
    {
        return Err(CapillarySurfaceError::NonFiniteDirection);
    }
    let k = kinematics(positions_m, velocities_m_s)?;
    let da = sub(position_direction[1], position_direction[0]);
    let db = sub(position_direction[2], position_direction[0]);
    let c1 = cross(da, k.b);
    let c2 = cross(k.a, db);
    let dc = core::array::from_fn(|d| c1[d] + c2[d]);
    let dc_normal = dot(k.normal, dc);
    let darea = 0.5 * dc_normal;
    let dn = core::array::from_fn(|d| (dc[d] - k.normal[d] * dc_normal) * (0.5 / k.area));
    let dprojector = core::array::from_fn(|i| {
        core::array::from_fn(|j| -dn[i] * k.normal[j] - k.normal[i] * dn[j])
    });
    let g1a = cross(db, k.normal);
    let g1b = cross(k.b, dn);
    let g2a = cross(dn, k.a);
    let g2b = cross(k.normal, da);
    let dg1: Vector =
        core::array::from_fn(|d| (0.5 * (g1a[d] + g1b[d]) - darea * k.gradients[1][d]) / k.area);
    let dg2: Vector =
        core::array::from_fn(|d| (0.5 * (g2a[d] + g2b[d]) - darea * k.gradients[2][d]) / k.area);
    let dvelocity = [
        sub(velocity_direction[1], velocity_direction[0]),
        sub(velocity_direction[2], velocity_direction[0]),
    ];
    let dgradient = core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            dvelocity[0][i] * k.gradients[1][j]
                + dvelocity[1][i] * k.gradients[2][j]
                + k.velocity_edges[0][i] * dg1[j]
                + k.velocity_edges[1][i] * dg2[j]
        })
    });
    let a = product(dprojector, k.velocity_gradient);
    let b = product(k.projector, dgradient);
    let drate: Matrix = core::array::from_fn(|i| {
        core::array::from_fn(|j| 0.5 * (a[i][j] + b[i][j] + a[j][i] + b[j][i]))
    });
    let dtrace = drate[0][0] + drate[1][1] + drate[2][2];
    let dtau: Matrix = core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            let dilation = dtrace * k.projector[i][j] + k.trace * dprojector[i][j];
            2.0 * viscosity.shear_n_s_m * (drate[i][j] - 0.5 * dilation)
                + viscosity.dilatational_n_s_m * dilation
        })
    });
    let tau = stress(&k, viscosity);
    let force_direction = |gradient, dgradient| {
        let a = apply(tau, gradient);
        let b = apply(dtau, gradient);
        let c = apply(tau, dgradient);
        core::array::from_fn(|d| -darea * a[d] - k.area * (b[d] + c[d]))
    };
    let f1: Vector = force_direction(k.gradients[1], dg1);
    let f2: Vector = force_direction(k.gradients[2], dg2);
    let forces = [core::array::from_fn(|d| -f1[d] - f2[d]), f1, f2];
    if forces.iter().flatten().any(|x| !x.is_finite()) {
        return Err(CapillarySurfaceError::Unrepresentable);
    }
    Ok(forces)
}

/// Assemble the surface-fluid stress on caller-owned fixed connectivity.
pub fn assemble_surface_viscosity(
    positions_m: &[Vector],
    velocities_m_s: &[Vector],
    faces: &[[usize; 3]],
    viscosity: SurfaceViscosityField<'_>,
) -> Result<SurfaceViscosityResponse, CapillarySurfaceError> {
    if positions_m.len() != velocities_m_s.len() {
        return Err(CapillarySurfaceError::LengthMismatch);
    }
    viscosity.validate(faces.len())?;
    if positions_m.iter().flatten().any(|x| !x.is_finite()) {
        return Err(CapillarySurfaceError::NonFinitePosition);
    }
    if velocities_m_s.iter().flatten().any(|x| !x.is_finite()) {
        return Err(CapillarySurfaceError::NonFiniteVelocity);
    }
    let mut forces = vec![[0.0; 3]; positions_m.len()];
    let mut dissipation = Sum::default();
    for (face, &indices) in faces.iter().enumerate() {
        for vertex in indices {
            if vertex >= positions_m.len() {
                return Err(CapillarySurfaceError::InvalidVertexIndex { face, vertex });
            }
        }
        let response = triangle_surface_viscosity(
            indices.map(|i| positions_m[i]),
            indices.map(|i| velocities_m_s[i]),
            viscosity.at(face),
        )?;
        for (node, force) in indices.into_iter().zip(response.forces_n) {
            for d in 0..3 {
                forces[node][d] += force[d];
            }
        }
        dissipation.add(response.dissipation_w);
    }
    let dissipation = dissipation.value();
    if forces.iter().flatten().any(|x| !x.is_finite()) || !dissipation.is_finite() {
        return Err(CapillarySurfaceError::Unrepresentable);
    }
    Ok(SurfaceViscosityResponse {
        forces_n: forces,
        dissipation_w: dissipation,
    })
}

#[cfg(test)]
mod tests;
