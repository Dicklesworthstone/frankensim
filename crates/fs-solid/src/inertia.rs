//! Shared consistent reference mass and backward-Euler momentum integration.
//!
//! Extracted from the small-strain inelastic owner without changing its mass
//! kernel, accumulation order, momentum residual or kinetic-work identity.
//! Geometry and constitutive history stay with the calling solid. Reference
//! mass is admitted once and may be retained across finite-deformation steps.

use fs_exec::Cx;
use fs_sparse::{Coo, Csr};

/// Inertia's contribution to the complete mechanical/thermal work balance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct InertialEnergy {
    /// Original closed-matter inventory; it is not updated from endpoint T.
    pub reference_mass_kg: f64,
    /// Kinetic energy of the accepted velocity.
    pub kinetic_before_j: f64,
    /// Kinetic energy of the proposed endpoint velocity.
    pub kinetic_after_j: f64,
    /// Evaluated as `average(v) dot M delta(v)` to avoid subtracting large K.
    pub kinetic_change_j: f64,
    /// Nonnegative backward-Euler kinetic loss, distinct from plastic loss.
    pub numerical_loss_j: f64,
    /// Endpoint inertial force dotted with the actual displacement increment.
    pub inertial_work_j: f64,
    /// Inertial work minus kinetic change and kinetic numerical loss.
    pub balance_error_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) enum InertiaError {
    InvalidInput { what: String },
    Cancelled,
}

impl core::fmt::Display for InertiaError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidInput { what } => write!(f, "invalid reference inertia: {what}"),
            Self::Cancelled => write!(f, "reference inertia cancelled"),
        }
    }
}

impl std::error::Error for InertiaError {}

fn invalid(what: &str) -> InertiaError {
    InertiaError::InvalidInput { what: what.into() }
}

fn checkpoint(cx: &Cx<'_>) -> Result<(), InertiaError> {
    cx.checkpoint().map_err(|_| InertiaError::Cancelled)
}

fn finite(value: f64) -> Result<f64, InertiaError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid("unrepresentable inertial force, tangent or work"))
    }
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).fold(0.0, |s, (x, y)| x.mul_add(*y, s))
}

/// Immutable consistent mass on an already admitted reference mesh.
#[derive(Debug)]
pub(crate) struct ReferenceMass {
    matrix: Csr,
    reference_mass_kg: f64,
}

/// Preserve the original tetrahedral constructor and all existing callers.
pub(crate) type P1ReferenceMass = ReferenceMass;

impl ReferenceMass {
    /// The caller owns reference geometry validation and assembly budgets.
    /// Entries are original connectivity and positive reference volume;
    /// density is original mass per reference volume in the same element order.
    pub(crate) fn new(
        geometry: impl ExactSizeIterator<Item = ([usize; 4], f64)>,
        densities: &[f64],
        ndofs: usize,
        cx: &Cx<'_>,
    ) -> Result<Self, InertiaError> {
        checkpoint(cx)?;
        if densities.len() != geometry.len()
            || densities.iter().any(|rho| !rho.is_finite() || *rho <= 0.0)
            || ndofs == 0
            || !ndofs.is_multiple_of(3)
        {
            return Err(invalid(
                "one positive reference density per tet and complete nodal components required",
            ));
        }
        let mut mass = Coo::new(ndofs, ndofs);
        let mut reference_mass_kg = 0.0;
        for ((nodes, volume), &density) in geometry.zip(densities) {
            checkpoint(cx)?;
            if !volume.is_finite() || volume <= 0.0 || nodes.iter().any(|&node| node >= ndofs / 3) {
                return Err(invalid(
                    "valid reference connectivity and positive volume required",
                ));
            }
            let local = crate::linear3::p1_tet_consistent_mass(volume, density);
            let element_mass = finite(density * volume)?;
            if element_mass <= 0.0 {
                return Err(invalid(
                    "each tetrahedron needs representable positive reference mass",
                ));
            }
            reference_mass_kg = finite(reference_mass_kg + element_mass)?;
            for a in 0..4 {
                for b in 0..4 {
                    for c in 0..3 {
                        let value = finite(local[3 * a + c][3 * b + c])?;
                        if value <= 0.0 {
                            return Err(invalid("consistent element mass underflowed"));
                        }
                        mass.push(3 * nodes[a] + c, 3 * nodes[b] + c, value);
                    }
                }
            }
        }
        Self::from_matrix(mass.assemble(), reference_mass_kg, ndofs, cx)
    }

    /// Consistent Hex8 mass at the original full eight-point quadrature.
    /// The admitted geometry supplies positive reference weights and the native
    /// trilinear shapes; neither current volume nor assumed dilatation enters.
    /// The matrix is exact for affine bricks and uses the original geometry's
    /// quadrature for warped cells. The same nonnegative mass action supports
    /// the existing momentum stopping gate and kinetic-work identity.
    pub(crate) fn from_hex8_quadrature(
        geometry: impl ExactSizeIterator<Item = ([usize; 8], [([f64; 8], f64); 8])>,
        densities: &[f64],
        ndofs: usize,
        cx: &Cx<'_>,
    ) -> Result<Self, InertiaError> {
        checkpoint(cx)?;
        if densities.len() != geometry.len()
            || densities.iter().any(|rho| !rho.is_finite() || *rho <= 0.0)
            || ndofs == 0
            || !ndofs.is_multiple_of(3)
        {
            return Err(invalid(
                "one positive reference density per hex and complete nodal components required",
            ));
        }
        let mut mass = Coo::new(ndofs, ndofs);
        let mut reference_mass_kg = 0.0;
        for ((nodes, points), &density) in geometry.zip(densities) {
            checkpoint(cx)?;
            if nodes.iter().any(|&node| node >= ndofs / 3)
                || (0..8).any(|a| nodes[..a].contains(&nodes[a]))
            {
                return Err(invalid(
                    "distinct original Hex8 nodes required for reference mass",
                ));
            }
            let mut local = [[0.0; 8]; 8];
            for (shapes, weight) in points {
                if !weight.is_finite()
                    || weight <= 0.0
                    || shapes
                        .iter()
                        .any(|n| !n.is_finite() || *n <= 0.0 || *n > 1.0)
                    || (shapes.iter().sum::<f64>() - 1.0).abs() > 64.0 * f64::EPSILON
                {
                    return Err(invalid(
                        "positive full Hex8 quadrature with partition of unity required",
                    ));
                }
                let point_mass = finite(density * weight)?;
                if point_mass <= 0.0 {
                    return Err(invalid("consistent Hex8 mass underflowed"));
                }
                reference_mass_kg = finite(reference_mass_kg + point_mass)?;
                for a in 0..8 {
                    for b in a..8 {
                        local[a][b] = (point_mass * shapes[a]).mul_add(shapes[b], local[a][b]);
                    }
                }
            }
            for a in 0..8 {
                for b in 0..8 {
                    let value = finite(local[a.min(b)][a.max(b)])?;
                    if value <= 0.0 {
                        return Err(invalid("consistent Hex8 mass entry underflowed"));
                    }
                    for component in 0..3 {
                        mass.push(3 * nodes[a] + component, 3 * nodes[b] + component, value);
                    }
                }
            }
        }
        Self::from_matrix(mass.assemble(), reference_mass_kg, ndofs, cx)
    }

    fn from_matrix(
        matrix: Csr,
        reference_mass_kg: f64,
        ndofs: usize,
        cx: &Cx<'_>,
    ) -> Result<Self, InertiaError> {
        for row in 0..ndofs {
            if row % 256 == 0 {
                checkpoint(cx)?;
            }
            if !matrix.get(row, row).is_finite()
                || matrix.get(row, row) <= 0.0
                || matrix.row(row).1.iter().any(|m| !m.is_finite())
            {
                return Err(invalid(
                    "every nodal component needs finite positive reference inertia",
                ));
            }
        }
        Ok(Self {
            matrix,
            reference_mass_kg,
        })
    }

    pub(crate) fn reference_mass_kg(&self) -> f64 {
        self.reference_mass_kg
    }

    pub(crate) fn ndofs(&self) -> usize {
        self.matrix.nrows()
    }
}

/// One endpoint momentum equation borrowing the immutable accepted motion/mass.
pub(crate) struct BackwardEulerInertia<'a> {
    mass: &'a P1ReferenceMass,
    old_u: &'a [f64],
    old_v: &'a [f64],
    inverse_dt: f64,
}

impl<'a> BackwardEulerInertia<'a> {
    pub(crate) fn new(
        mass: &'a P1ReferenceMass,
        old_u: &'a [f64],
        old_v: &'a [f64],
        dt_s: f64,
        cx: &Cx<'_>,
    ) -> Result<Self, InertiaError> {
        checkpoint(cx)?;
        if old_u.len() != mass.ndofs()
            || old_v.len() != old_u.len()
            || old_u.iter().chain(old_v).any(|v| !v.is_finite())
            || !dt_s.is_finite()
            || dt_s <= 0.0
        {
            return Err(invalid(
                "finite full displacement/velocity and positive dt required",
            ));
        }
        let inverse_dt = finite(dt_s.recip())?;
        let coefficient = finite(inverse_dt * inverse_dt)?;
        if coefficient <= 0.0 {
            return Err(invalid(
                "inertial tangent coefficient must be representable and positive",
            ));
        }
        Ok(Self {
            mass,
            old_u,
            old_v,
            inverse_dt,
        })
    }

    fn check_length(&self, x: &[f64]) -> Result<(), InertiaError> {
        if x.len() != self.mass.ndofs() {
            return Err(invalid("inertia requires every nodal component"));
        }
        Ok(())
    }

    // Normalize only physical forces on free rows. Elastic stiffness times
    // rigid travel is unrelated to momentum balance and can accept v_next=0.
    // The positive mass action retains a reference under momentum cancellation;
    // unrelated prescribed-support forces cannot inflate the stopping gate.
    pub(crate) fn force_scale(
        &self,
        u: &[f64],
        free: &[usize],
        elastic_free_norm: f64,
        external_free_norm: f64,
        cx: &Cx<'_>,
    ) -> Result<f64, InertiaError> {
        self.check_length(u)?;
        let mut acceleration_scale = Vec::with_capacity(u.len());
        for (i, ((u, old_u), old_v)) in u.iter().zip(self.old_u).zip(self.old_v).enumerate() {
            if i % 256 == 0 {
                checkpoint(cx)?;
            }
            acceleration_scale.push(finite(
                (finite(u - old_u)?.abs() * self.inverse_dt + old_v.abs()) * self.inverse_dt,
            )?);
        }
        // Both admitted consistent reference mass rules have nonnegative entries,
        // so its existing action is exactly |M| times this positive vector.
        let force = self.action(&acceleration_scale, cx)?;
        let mut inertial_norm = 0.0_f64;
        for &d in free {
            let force = force
                .get(d)
                .ok_or_else(|| invalid("free inertia row is outside the mesh"))?;
            inertial_norm = inertial_norm.hypot(*force);
        }
        Ok(finite(inertial_norm)?
            .max(finite(elastic_free_norm)?)
            .max(finite(external_free_norm)?)
            .max(f64::MIN_POSITIVE))
    }

    pub(crate) fn action(&self, x: &[f64], cx: &Cx<'_>) -> Result<Vec<f64>, InertiaError> {
        self.check_length(x)?;
        let mut output = Vec::with_capacity(x.len());
        for row in 0..x.len() {
            if row % 256 == 0 {
                checkpoint(cx)?;
            }
            let (columns, values) = self.mass.matrix.row(row);
            output.push(finite(
                columns
                    .iter()
                    .zip(values)
                    .fold(0.0, |sum, (&column, &value)| value.mul_add(x[column], sum)),
            )?);
        }
        Ok(output)
    }

    pub(crate) fn add_force(
        &self,
        u: &[f64],
        force: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), InertiaError> {
        self.check_length(u)?;
        self.check_length(force)?;
        let acceleration = u
            .iter()
            .zip(self.old_u)
            .zip(self.old_v)
            .map(|((&u, &old_u), &old_v)| {
                finite(((u - old_u) * self.inverse_dt - old_v) * self.inverse_dt)
            })
            .collect::<Result<Vec<_>, _>>()?;
        let inertial = self.action(&acceleration, cx)?;
        for (force, inertia) in force.iter_mut().zip(inertial) {
            *force = finite(*force + inertia)?;
        }
        Ok(())
    }

    pub(crate) fn add_tangent(&self, matrix: &mut Coo, cx: &Cx<'_>) -> Result<(), InertiaError> {
        let coefficient = self.inverse_dt * self.inverse_dt;
        for row in 0..self.mass.matrix.nrows() {
            if row % 256 == 0 {
                checkpoint(cx)?;
            }
            let (columns, values) = self.mass.matrix.row(row);
            for (&column, &value) in columns.iter().zip(values) {
                let scaled = finite(coefficient * value)?;
                if scaled <= 0.0 {
                    return Err(invalid("scaled consistent inertia underflowed"));
                }
                matrix.push(row, column, scaled);
            }
        }
        Ok(())
    }

    /// Gradient is the endpoint inertial force; Hessian is M/dt².
    /// This reduced kinetic potential is for Newton globalization, and is
    /// distinct from endpoint physical kinetic energy and plastic dissipation.
    pub(crate) fn potential(&self, u: &[f64], cx: &Cx<'_>) -> Result<f64, InertiaError> {
        self.check_length(u)?;
        let delta: Vec<_> = u
            .iter()
            .zip(self.old_u)
            .zip(self.old_v)
            .map(|((&u, &old_u), &old_v)| finite((u - old_u) * self.inverse_dt - old_v))
            .collect::<Result<_, _>>()?;
        let impulse = self.action(&delta, cx)?;
        let potential = finite(0.5 * dot(&delta, &impulse))?;
        if potential < 0.0 {
            return Err(invalid("consistent inertial potential must be nonnegative"));
        }
        Ok(potential)
    }

    pub(crate) fn motion(
        &self,
        u: &[f64],
        cx: &Cx<'_>,
    ) -> Result<(Vec<[f64; 3]>, InertialEnergy), InertiaError> {
        self.check_length(u)?;
        let du: Vec<_> = u
            .iter()
            .zip(self.old_u)
            .map(|(u, old)| finite(u - old))
            .collect::<Result<_, _>>()?;
        let velocity: Vec<_> = du
            .iter()
            .map(|du| finite(du * self.inverse_dt))
            .collect::<Result<_, _>>()?;
        let delta: Vec<_> = velocity
            .iter()
            .zip(self.old_v)
            .map(|(v, old)| finite(v - old))
            .collect::<Result<_, _>>()?;
        let mean: Vec<_> = velocity
            .iter()
            .zip(self.old_v)
            .map(|(v, old)| finite(0.5 * v + 0.5 * old))
            .collect::<Result<_, _>>()?;
        let old_momentum = self.action(self.old_v, cx)?;
        let next_momentum = self.action(&velocity, cx)?;
        let impulse = self.action(&delta, cx)?;
        let kinetic_before_j = finite(0.5 * dot(self.old_v, &old_momentum))?;
        let kinetic_after_j = finite(0.5 * dot(&velocity, &next_momentum))?;
        let kinetic_change_j = finite(dot(&mean, &impulse))?;
        let numerical_loss_j = finite(0.5 * dot(&delta, &impulse))?;
        let mut inertia = vec![0.0; u.len()];
        self.add_force(u, &mut inertia, cx)?;
        let inertial_work_j = finite(dot(&du, &inertia))?;
        let balance_error_j = finite(inertial_work_j - kinetic_change_j - numerical_loss_j)?;
        let scale = kinetic_before_j
            .max(kinetic_after_j)
            .max(inertial_work_j.abs())
            .max(numerical_loss_j);
        if kinetic_before_j < 0.0
            || kinetic_after_j < 0.0
            || numerical_loss_j < 0.0
            || balance_error_j.abs() > 2048.0 * f64::EPSILON * scale
        {
            return Err(invalid(
                "consistent inertia misses its discrete kinetic work identity",
            ));
        }
        Ok((
            velocity.as_chunks::<3>().0.to_vec(),
            InertialEnergy {
                reference_mass_kg: self.mass.reference_mass_kg,
                kinetic_before_j,
                kinetic_after_j,
                kinetic_change_j,
                numerical_loss_j,
                inertial_work_j,
                balance_error_j,
            },
        ))
    }
}
