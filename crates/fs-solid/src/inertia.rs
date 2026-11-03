//! Shared consistent reference mass and native momentum integration.
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
    /// Nonnegative numerical kinetic loss; zero for midpoint kinematics.
    /// Backward-Euler loss is kept separate from physical material loss.
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

/// Exact scalar P1 mass tile, shared by linear and finite-strain mechanics.
pub(crate) struct P1DensityMass {
    pub coefficients: [[f64; 4]; 4],
    pub mean_density_kg_m3: f64,
    pub mass_kg: f64,
}

/// Integrate the cubic polynomial `rho * N_a * N_b` on a reference tetrahedron.
/// Normalize before summing density so a representable mass does not overflow
/// merely because the four positive density samples have a large sum.
pub(crate) fn p1_tet_nodal_mass(
    volume_m3: f64,
    rho: [f64; 4],
) -> Result<P1DensityMass, InertiaError> {
    if !volume_m3.is_finite()
        || volume_m3 <= 0.0
        || rho.iter().any(|value| !value.is_finite() || *value <= 0.0)
    {
        return Err(invalid("positive reference volume and nodal densities required"));
    }
    let uniform = rho.iter().all(|value| value.to_bits() == rho[0].to_bits());
    let scale = rho.into_iter().fold(0.0_f64, f64::max);
    let normalized = rho.map(|value| value / scale);
    let total: f64 = normalized.iter().sum();
    let mean_density_kg_m3 = if uniform {
        rho[0]
    } else {
        scale * (0.25 * total)
    };
    let mass_kg = finite(volume_m3 * mean_density_kg_m3)?;
    if mass_kg <= 0.0 {
        return Err(invalid(
            "each tetrahedron needs representable positive reference mass",
        ));
    }
    let fractions = normalized.map(|value| value / total);
    let uniform_scale = crate::linear3::p1_tet_mass_scale(volume_m3, rho[0]);
    let mut coefficients = [[0.0; 4]; 4];
    for a in 0..4 {
        for b in 0..4 {
            // Diagonal: V*(sum(rho)+2*rho_a)/60.
            // Off diagonal: V*(sum(rho)+rho_a+rho_b)/120.
            let value = finite(if uniform {
                uniform_scale * if a == b { 2.0 } else { 1.0 }
            } else if a == b {
                (mass_kg / 15.0) * (1.0 + 2.0 * fractions[a])
            } else {
                (mass_kg / 30.0) * (1.0 + (fractions[a] + fractions[b]))
            })?;
            if value <= 0.0 {
                return Err(invalid("consistent graded-density mass underflowed"));
            }
            coefficients[a][b] = value;
        }
    }
    Ok(P1DensityMass {
        coefficients,
        mean_density_kg_m3,
        mass_kg,
    })
}

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

    /// Exact consistent tetrahedral mass for a positive P1 reference density.
    /// Nodal density belongs to the original mesh; current volume never enters.
    /// Return cell means separately for inventory reporting, not as a substitute
    /// for the retained nodal field or its consistent operator.
    pub(crate) fn from_nodal_density(
        geometry: impl ExactSizeIterator<Item = ([usize; 4], f64)>,
        densities: &[f64],
        ndofs: usize,
        cx: &Cx<'_>,
    ) -> Result<(Self, Vec<f64>), InertiaError> {
        checkpoint(cx)?;
        if ndofs == 0
            || !ndofs.is_multiple_of(3)
            || densities.len() != ndofs / 3
            || densities.iter().any(|rho| !rho.is_finite() || *rho <= 0.0)
        {
            return Err(invalid(
                "one positive reference density per node and complete nodal components required",
            ));
        }
        if densities
            .iter()
            .all(|rho| rho.to_bits() == densities[0].to_bits())
        {
            let means = vec![densities[0]; geometry.len()];
            return Self::new(geometry, &means, ndofs, cx).map(|mass| (mass, means));
        }
        let mut matrix = Coo::new(ndofs, ndofs);
        let mut means = Vec::with_capacity(geometry.len());
        let mut reference_mass_kg = 0.0;
        for (nodes, volume) in geometry {
            checkpoint(cx)?;
            if !volume.is_finite() || volume <= 0.0 || nodes.iter().any(|&node| node >= ndofs / 3) {
                return Err(invalid(
                    "valid reference connectivity and positive volume required",
                ));
            }
            let local = p1_tet_nodal_mass(volume, nodes.map(|node| densities[node]))?;
            means.push(local.mean_density_kg_m3);
            reference_mass_kg = finite(reference_mass_kg + local.mass_kg)?;
            for a in 0..4 {
                for b in 0..4 {
                    for component in 0..3 {
                        matrix.push(
                            3 * nodes[a] + component,
                            3 * nodes[b] + component,
                            local.coefficients[a][b],
                        );
                    }
                }
            }
        }
        Self::from_matrix(matrix.assemble(), reference_mass_kg, ndofs, cx).map(|mass| (mass, means))
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

    /// Consistent P2 mass from the caller's cached, admitted Tet10 geometry.
    /// Each positive reference weight multiplies the actual signed shape pair;
    /// row-sum lumping is invalid because quadratic corner integrals are negative.
    /// Density is original mass per reference volume in the same element order.
    pub(crate) fn from_tet10_quadrature<Q>(
        geometry: impl ExactSizeIterator<Item = ([usize; 10], Q)>,
        densities: &[f64],
        ndofs: usize,
        maximum_assembly_bytes: usize,
        cx: &Cx<'_>,
    ) -> Result<Self, InertiaError>
    where
        Q: IntoIterator<Item = ([f64; 10], f64)>,
    {
        checkpoint(cx)?;
        if geometry.len() == 0
            || densities.len() != geometry.len()
            || densities.iter().any(|rho| !rho.is_finite() || *rho <= 0.0)
            || ndofs == 0
            || !ndofs.is_multiple_of(3)
        {
            return Err(invalid(
                "one positive reference density per Tet10 and complete nodal components required",
            ));
        }
        Self::from_tet10_density_quadrature(
            geometry.zip(densities).map(|((nodes, points), &density)| {
                (
                    nodes,
                    points
                        .into_iter()
                        .map(move |(shapes, weight)| (shapes, weight, density)),
                )
            }),
            ndofs,
            maximum_assembly_bytes,
            cx,
        )
    }

    /// The same signed P2 mass kernel with an explicit positive reference
    /// density at each original quadrature point. The caller admits its density
    /// interpolation and quadrature degree; current volume never enters.
    /// Scalar callers retain their original point-mass accumulation order.
    pub(crate) fn from_tet10_density_quadrature<Q>(
        geometry: impl ExactSizeIterator<Item = ([usize; 10], Q)>,
        ndofs: usize,
        maximum_assembly_bytes: usize,
        cx: &Cx<'_>,
    ) -> Result<Self, InertiaError>
    where
        Q: IntoIterator<Item = ([f64; 10], f64, f64)>,
    {
        checkpoint(cx)?;
        if geometry.len() == 0 || ndofs == 0 || !ndofs.is_multiple_of(3) {
            return Err(invalid(
                "nonempty Tet10 reference quadrature and complete nodal components required",
            ));
        }
        // Three uncoupled component blocks, with conservative room for COO
        // growth, stable-sort scratch, CSR and its row offsets. Quadrature is
        // borrowed and accumulated into one fixed local upper triangle.
        let estimated = geometry
            .len()
            .checked_mul(10 * 10 * 3 * 128)
            .and_then(|entries| {
                ndofs
                    .checked_add(1)
                    .and_then(|rows| rows.checked_mul(core::mem::size_of::<usize>()))
                    .and_then(|rows| entries.checked_add(rows))
            })
            .and_then(|bytes| bytes.checked_add(core::mem::size_of::<[[f64; 10]; 10]>()))
            .ok_or_else(|| invalid("Tet10 reference mass assembly size overflow"))?;
        if estimated > maximum_assembly_bytes {
            return Err(invalid(
                "Tet10 reference mass exceeds its assembly byte cap",
            ));
        }
        let mut mass = Coo::new(ndofs, ndofs);
        let mut reference_mass_kg = 0.0;
        for (nodes, points) in geometry {
            checkpoint(cx)?;
            if nodes.iter().any(|&node| node >= ndofs / 3)
                || (0..10).any(|a| nodes[..a].contains(&nodes[a]))
            {
                return Err(invalid(
                    "distinct original Tet10 nodes required for reference mass",
                ));
            }
            let mut local = [[0.0; 10]; 10];
            for (shapes, weight, density) in points {
                checkpoint(cx)?;
                let partition = shapes.iter().sum::<f64>();
                if !weight.is_finite()
                    || weight <= 0.0
                    || !density.is_finite()
                    || density <= 0.0
                    || shapes.iter().any(|n| !n.is_finite())
                    || !partition.is_finite()
                    || (partition - 1.0).abs() > 64.0 * f64::EPSILON
                {
                    return Err(invalid(
                        "positive Tet10 reference density and quadrature with finite partition of unity required",
                    ));
                }
                let point_mass = finite(density * weight)?;
                if point_mass <= 0.0 {
                    return Err(invalid("consistent Tet10 mass underflowed"));
                }
                reference_mass_kg = finite(reference_mass_kg + point_mass)?;
                for a in 0..10 {
                    for b in a..10 {
                        local[a][b] = (point_mass * shapes[a]).mul_add(shapes[b], local[a][b]);
                    }
                }
            }
            for a in 0..10 {
                if !local[a][a].is_finite() || local[a][a] <= 0.0 {
                    return Err(invalid(
                        "every Tet10 shape needs positive reference inertia",
                    ));
                }
                for b in 0..10 {
                    let value = finite(local[a.min(b)][a.max(b)])?;
                    if value != 0.0 {
                        for component in 0..3 {
                            mass.push(3 * nodes[a] + component, 3 * nodes[b] + component, value);
                        }
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

    /// The same immutable mass action used by momentum and acceleration loads.
    pub(crate) fn action(&self, x: &[f64], cx: &Cx<'_>) -> Result<Vec<f64>, InertiaError> {
        if x.len() != self.ndofs() {
            return Err(invalid("inertia requires every nodal component"));
        }
        let mut output = Vec::with_capacity(x.len());
        for row in 0..x.len() {
            if row % 256 == 0 {
                checkpoint(cx)?;
            }
            let (columns, values) = self.matrix.row(row);
            output.push(finite(
                columns
                    .iter()
                    .zip(values)
                    .fold(0.0, |sum, (&column, &value)| value.mul_add(x[column], sum)),
            )?);
        }
        Ok(output)
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
    // The absolute mass action retains a reference under momentum cancellation;
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
        // Signed P2 entries must not cancel this bound. For the original
        // positive P1/Hex8 rules, abs retains the identical multiply-add order.
        let mut inertial_norm = 0.0_f64;
        for (i, &d) in free.iter().enumerate() {
            if i % 256 == 0 {
                checkpoint(cx)?;
            }
            if d >= self.mass.ndofs() {
                return Err(invalid("free inertia row is outside the mesh"));
            }
            let (columns, values) = self.mass.matrix.row(d);
            let force = finite(
                columns
                    .iter()
                    .zip(values)
                    .fold(0.0, |sum, (&column, &value)| {
                        value.abs().mul_add(acceleration_scale[column], sum)
                    }),
            )?;
            inertial_norm = inertial_norm.hypot(force);
        }
        Ok(finite(inertial_norm)?
            .max(finite(elastic_free_norm)?)
            .max(finite(external_free_norm)?)
            .max(f64::MIN_POSITIVE))
    }

    pub(crate) fn action(&self, x: &[f64], cx: &Cx<'_>) -> Result<Vec<f64>, InertiaError> {
        self.mass.action(x, cx)
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
                if value == 0.0 {
                    continue;
                }
                let scaled = finite(coefficient * value)?;
                if scaled == 0.0 {
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

/// Centered momentum equation for a work-conjugate discrete material gradient.
///
/// `v_next = 2*(u_next-u_old)/dt-v_old` makes the inertial work equal the
/// physical kinetic change, with no numerical kinetic loss. Conservation of
/// total energy and momenta still requires the caller's actual discrete stress,
/// loads and nonlinear solve; this adapter does not remove frequency warping.
pub(crate) struct MidpointInertia<'a> {
    original: BackwardEulerInertia<'a>,
}

impl<'a> MidpointInertia<'a> {
    pub(crate) fn new(
        mass: &'a ReferenceMass,
        old_u: &'a [f64],
        old_v: &'a [f64],
        dt_s: f64,
        cx: &Cx<'_>,
    ) -> Result<Self, InertiaError> {
        let original = BackwardEulerInertia::new(mass, old_u, old_v, dt_s, cx)?;
        finite(2.0 * (original.inverse_dt * original.inverse_dt))?;
        Ok(Self { original })
    }

    pub(crate) fn action(&self, x: &[f64], cx: &Cx<'_>) -> Result<Vec<f64>, InertiaError> {
        self.original.action(x, cx)
    }

    /// Physical free-row reference bound; no elastic stiffness times rigid travel.
    pub(crate) fn force_scale(
        &self,
        u: &[f64],
        free: &[usize],
        elastic_free_norm: f64,
        external_free_norm: f64,
        cx: &Cx<'_>,
    ) -> Result<f64, InertiaError> {
        let inertial = finite(2.0 * self.original.force_scale(u, free, 0.0, 0.0, cx)?)?;
        Ok(inertial
            .max(finite(elastic_free_norm)?)
            .max(finite(external_free_norm)?))
    }

    pub(crate) fn add_force(
        &self,
        u: &[f64],
        force: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), InertiaError> {
        self.original.check_length(u)?;
        self.original.check_length(force)?;
        let inverse_dt = self.original.inverse_dt;
        let mut acceleration = Vec::with_capacity(u.len());
        for (i, ((&u, &old_u), &old_v)) in u
            .iter()
            .zip(self.original.old_u)
            .zip(self.original.old_v)
            .enumerate()
        {
            if i % 256 == 0 {
                checkpoint(cx)?;
            }
            acceleration.push(finite(
                2.0 * (((u - old_u) * inverse_dt - old_v) * inverse_dt),
            )?);
        }
        let inertial = self.action(&acceleration, cx)?;
        for (force, inertial) in force.iter_mut().zip(inertial) {
            *force = finite(*force + inertial)?;
        }
        Ok(())
    }

    pub(crate) fn add_tangent(&self, matrix: &mut Coo, cx: &Cx<'_>) -> Result<(), InertiaError> {
        let coefficient = 2.0 * (self.original.inverse_dt * self.original.inverse_dt);
        for row in 0..self.original.mass.matrix.nrows() {
            if row % 256 == 0 {
                checkpoint(cx)?;
            }
            let (columns, values) = self.original.mass.matrix.row(row);
            for (&column, &value) in columns.iter().zip(values) {
                if value == 0.0 {
                    continue;
                }
                let scaled = finite(coefficient * value)?;
                if scaled == 0.0 {
                    return Err(invalid("scaled consistent midpoint inertia underflowed"));
                }
                matrix.push(row, column, scaled);
            }
        }
        Ok(())
    }

    pub(crate) fn motion(
        &self,
        u: &[f64],
        cx: &Cx<'_>,
    ) -> Result<(Vec<[f64; 3]>, InertialEnergy), InertiaError> {
        self.original.check_length(u)?;
        let mut du = Vec::with_capacity(u.len());
        let mut velocity = Vec::with_capacity(u.len());
        let mut delta = Vec::with_capacity(u.len());
        let mut mean = Vec::with_capacity(u.len());
        for (i, ((&u, &old_u), &old_v)) in u
            .iter()
            .zip(self.original.old_u)
            .zip(self.original.old_v)
            .enumerate()
        {
            if i % 256 == 0 {
                checkpoint(cx)?;
            }
            let displacement = finite(u - old_u)?;
            let next_v = finite(2.0 * (displacement * self.original.inverse_dt) - old_v)?;
            du.push(displacement);
            velocity.push(next_v);
            delta.push(finite(next_v - old_v)?);
            mean.push(finite(0.5 * next_v + 0.5 * old_v)?);
        }
        let old_momentum = self.action(self.original.old_v, cx)?;
        let next_momentum = self.action(&velocity, cx)?;
        let impulse = self.action(&delta, cx)?;
        let kinetic_before_j = finite(0.5 * dot(self.original.old_v, &old_momentum))?;
        let kinetic_after_j = finite(0.5 * dot(&velocity, &next_momentum))?;
        let kinetic_change_j = finite(dot(&mean, &impulse))?;
        let mut inertial = vec![0.0; u.len()];
        self.add_force(u, &mut inertial, cx)?;
        let inertial_work_j = finite(dot(&du, &inertial))?;
        let balance_error_j = finite(inertial_work_j - kinetic_change_j)?;
        let scale = kinetic_before_j
            .max(kinetic_after_j)
            .max(kinetic_change_j.abs())
            .max(inertial_work_j.abs());
        if kinetic_before_j < 0.0
            || kinetic_after_j < 0.0
            || balance_error_j.abs() > 2048.0 * f64::EPSILON * scale
        {
            return Err(invalid(
                "consistent midpoint inertia misses its discrete kinetic work identity",
            ));
        }
        checkpoint(cx)?;
        Ok((
            velocity.as_chunks::<3>().0.to_vec(),
            InertialEnergy {
                reference_mass_kg: self.original.mass.reference_mass_kg,
                kinetic_before_j,
                kinetic_after_j,
                kinetic_change_j,
                numerical_loss_j: 0.0,
                inertial_work_j,
                balance_error_j,
            },
        ))
    }
}

#[cfg(test)]
mod density_tests {
    use super::*;
    use fs_alloc::{ArenaConfig, ArenaPool};
    use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};

    fn with_cx(run: impl FnOnce(&Cx<'_>)) {
        let gate = CancelGate::new_clock_free();
        ArenaPool::new(ArenaConfig::default()).scope(|arena| {
            run(&Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 31,
                    kernel_id: 7,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            ));
        });
    }
    fn graded(cx: &Cx<'_>) -> ReferenceMass {
        ReferenceMass::from_nodal_density(
            [([0, 1, 2, 3], 1.0 / 6.0)].into_iter(),
            &[3.0, 15.0, 3.0, 3.0],
            12,
            cx,
        )
        .unwrap()
        .0
    }
    fn close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() < 2e-13 * expected.abs().max(1.0),
            "{actual:.17e} != {expected:.17e}"
        );
    }
    // Independent simplex monomial integral: 6V*product(p_i!)/(sum(p_i)+3)!.
    fn exact_entry(a: usize, b: usize) -> f64 {
        let factorial = [1.0, 1.0, 2.0, 6.0];
        [3.0, 15.0, 3.0, 3.0]
            .into_iter()
            .enumerate()
            .map(|(k, rho)| {
                let mut powers = [0; 4];
                powers[a] += 1;
                powers[b] += 1;
                powers[k] += 1;
                rho * powers.into_iter().map(|p| factorial[p]).product::<f64>() / 720.0
            })
            .sum()
    }

    #[test]
    fn graded_mass_matches_cubic_integrals_and_density_weighted_rigid_energy() {
        with_cx(|cx| {
            let mass = graded(cx);
            close(mass.reference_mass_kg(), 1.0);
            for a in 0..4 {
                for b in 0..4 {
                    for i in 0..3 {
                        for j in 0..3 {
                            close(
                                mass.matrix.get(3 * a + i, 3 * b + j),
                                if i == j { exact_entry(a, b) } else { 0.0 },
                            );
                        }
                    }
                    assert_eq!(
                        mass.matrix.get(3 * a, 3 * b).to_bits(),
                        mass.matrix.get(3 * b, 3 * a).to_bits()
                    );
                }
            }
            // v = translation + omega cross x; the complete rigid kinetic
            // energy depends on COM and every inertia cross term, not just M.
            let v = [
                [0.3, -0.4, 0.2],
                [0.3, 0.1, 0.5],
                [-0.2, -0.4, 0.4],
                [0.0, -0.6, 0.2],
            ];
            let momentum = mass.action(v.as_flattened(), cx).unwrap();
            let actual = 0.5 * dot(v.as_flattened(), &momentum);
            let mut expected = 0.0;
            for a in 0..4 {
                for b in 0..4 {
                    expected +=
                        0.5 * exact_entry(a, b) * (0..3).map(|i| v[a][i] * v[b][i]).sum::<f64>();
                }
            }
            close(actual, expected);
            let uniform =
                ReferenceMass::new([([0, 1, 2, 3], 1.0 / 6.0)].into_iter(), &[6.0], 12, cx)
                    .unwrap();
            let averaged = 0.5
                * dot(
                    v.as_flattened(),
                    &uniform.action(v.as_flattened(), cx).unwrap(),
                );
            assert!((actual - averaged).abs() > 1e-3);
        });
    }

    #[test]
    fn graded_gravity_uses_native_momentum_and_kinetic_work_with_conserved_mass() {
        with_cx(|cx| {
            let mass = graded(cx);
            let acceleration = [1.2, -9.81, 0.4];
            let load = mass.action([acceleration; 4].as_flattened(), cx).unwrap();
            for (node, share) in [0.225, 0.325, 0.225, 0.225].into_iter().enumerate() {
                for i in 0..3 {
                    close(load[3 * node + i], share * acceleration[i]);
                }
            }
            let old_v = [[0.2, -0.1, 0.05]; 4];
            let old_u = [0.0; 12];
            let dt = 0.02;
            let v: [f64; 3] = core::array::from_fn(|i| old_v[0][i] + dt * acceleration[i]);
            let displacement = [v.map(|value| dt * value); 4];
            let inertia =
                BackwardEulerInertia::new(&mass, &old_u, old_v.as_flattened(), dt, cx).unwrap();
            let mut force = [0.0; 12];
            inertia
                .add_force(displacement.as_flattened(), &mut force, cx)
                .unwrap();
            for (a, b) in force.into_iter().zip(load) {
                close(a, b);
            }
            let (velocity, energy) = inertia.motion(displacement.as_flattened(), cx).unwrap();
            for actual in velocity {
                for i in 0..3 {
                    close(actual[i], v[i]);
                }
            }
            close(energy.reference_mass_kg, 1.0);
            close(energy.kinetic_before_j, 0.5 * dot(&old_v[0], &old_v[0]));
            close(energy.kinetic_after_j, 0.5 * dot(&v, &v));
            close(
                energy.numerical_loss_j,
                0.5 * dt * dt * dot(&acceleration, &acceleration),
            );
            close(energy.inertial_work_j, dt * dot(&v, &acceleration));
            close(energy.balance_error_j, 0.0);
        });
    }

    #[test]
    fn uniform_nodal_admission_keeps_existing_mass_bits_and_invalid_fields_refuse() {
        with_cx(|cx| {
            let geometry = [([0, 1, 2, 3], 0.17)];
            let uniform = ReferenceMass::new(geometry.into_iter(), &[7.3], 12, cx).unwrap();
            let (nodal, means) =
                ReferenceMass::from_nodal_density(geometry.into_iter(), &[7.3; 4], 12, cx).unwrap();
            assert_eq!(means, [7.3]);
            assert_eq!(
                uniform.reference_mass_kg.to_bits(),
                nodal.reference_mass_kg.to_bits()
            );
            assert_eq!(uniform.matrix, nodal.matrix);
            for rho in [
                &[1.0, 2.0, 3.0][..],
                &[1.0, 0.0, 2.0, 3.0],
                &[1.0, f64::NAN, 2.0, 3.0],
            ] {
                assert!(
                    ReferenceMass::from_nodal_density(geometry.into_iter(), rho, 12, cx).is_err()
                );
            }
            assert!(nodal.action(&[1.0; 11], cx).is_err());
            assert!(nodal.action(&[f64::INFINITY; 12], cx).is_err());
        });
    }
}

#[cfg(test)]
mod tet10_tests {
    use super::*;
    use crate::linear3::tet10::{TET10_EDGES, tet10_shapes};
    use fs_alloc::{ArenaConfig, ArenaPool};
    use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
    use fs_feec::highorder::simplex::duffy_quadrature;

    const NODES: [usize; 10] = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9];
    fn with_cx(cancelled: bool, run: impl FnOnce(&Cx<'_>)) {
        let gate = CancelGate::new_clock_free();
        if cancelled {
            gate.request();
        }
        ArenaPool::new(ArenaConfig::default()).scope(|arena| {
            run(&Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 31,
                    kernel_id: 8,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            ));
        });
    }
    fn quadrature() -> Vec<([f64; 10], f64)> {
        duffy_quadrature(4)
            .into_iter()
            .map(|(lambda, weight)| (tet10_shapes(lambda).0, weight))
            .collect()
    }
    fn mass(cx: &Cx<'_>) -> ReferenceMass {
        ReferenceMass::from_tet10_quadrature(
            [(NODES, quadrature())].into_iter(),
            &[6.0],
            30,
            1 << 20,
            cx,
        )
        .unwrap()
    }
    fn close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() < 2e-13 * expected.abs().max(1.0),
            "{actual:.17e} != {expected:.17e}"
        );
    }
    // Independent exact unit-mass P2 simplex integrals. Expanding each pair of
    // barycentric polynomials and using integral(product(lambda_i^p_i)) =
    // product(p_i!)/(3+sum(p_i))! gives these rational coefficients / 420.
    fn exact_entry(a: usize, b: usize) -> f64 {
        let coefficient = match (a < 4, b < 4) {
            (true, true) => {
                if a == b {
                    6.0
                } else {
                    1.0
                }
            }
            (true, false) => {
                if TET10_EDGES[b - 4].contains(&a) {
                    -4.0
                } else {
                    -6.0
                }
            }
            (false, true) => return exact_entry(b, a),
            (false, false) => {
                if a == b {
                    32.0
                } else if TET10_EDGES[a - 4]
                    .iter()
                    .any(|v| TET10_EDGES[b - 4].contains(v))
                {
                    16.0
                } else {
                    8.0
                }
            }
        };
        coefficient / 420.0
    }

    #[test]
    fn quadratic_mass_matches_rational_integrals_and_linear_velocity_moments() {
        with_cx(false, |cx| {
            let mass = mass(cx);
            close(mass.reference_mass_kg(), 1.0);
            for a in 0..10 {
                for b in 0..10 {
                    for i in 0..3 {
                        for j in 0..3 {
                            close(
                                mass.matrix.get(3 * a + i, 3 * b + j),
                                if i == j { exact_entry(a, b) } else { 0.0 },
                            );
                        }
                    }
                    assert_eq!(
                        mass.matrix.get(3 * a, 3 * b).to_bits(),
                        mass.matrix.get(3 * b, 3 * a).to_bits()
                    );
                }
            }
            let constant = mass.action(&[1.0; 30], cx).unwrap();
            for node in 0..10 {
                close(
                    constant[3 * node],
                    if node < 4 { -1.0 / 20.0 } else { 1.0 / 5.0 },
                );
            }
            // v(X)=X on the reference unit tetrahedron: each mean square
            // coordinate is 1/10. This exercises the complete quadratic mass.
            let positions = [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
                [0.5, 0.0, 0.0],
                [0.5, 0.5, 0.0],
                [0.0, 0.5, 0.0],
                [0.0, 0.0, 0.5],
                [0.5, 0.0, 0.5],
                [0.0, 0.5, 0.5],
            ];
            let v = positions.as_flattened();
            close(0.5 * dot(v, &mass.action(v, cx).unwrap()), 3.0 / 20.0);
        });
    }

    #[test]
    fn signed_mass_keeps_tangent_force_bound_and_kinetic_identity() {
        with_cx(false, |cx| {
            let mass = mass(cx);
            let (old_u, old_v, u) = ([0.0; 30], [1.0; 30], [0.1; 30]);
            let inertia = BackwardEulerInertia::new(&mass, &old_u, &old_v, 0.5, cx).unwrap();
            let mut tangent = Coo::new(30, 30);
            inertia.add_tangent(&mut tangent, cx).unwrap();
            let tangent = tangent.assemble();
            assert!(tangent.get(0, 12) < 0.0);
            for a in 0..10 {
                for b in 0..10 {
                    close(tangent.get(3 * a, 3 * b), 4.0 * exact_entry(a, b));
                }
            }
            let direction: [f64; 30] = core::array::from_fn(|d| (d as f64 - 15.0) * 0.01);
            let (mut plus, mut minus) = ([0.0; 30], [0.0; 30]);
            let h = 1e-6;
            inertia
                .add_force(
                    &core::array::from_fn::<_, 30, _>(|d| u[d] + h * direction[d]),
                    &mut plus,
                    cx,
                )
                .unwrap();
            inertia
                .add_force(
                    &core::array::from_fn::<_, 30, _>(|d| u[d] - h * direction[d]),
                    &mut minus,
                    cx,
                )
                .unwrap();
            for row in 0..30 {
                let expected: f64 = (0..30)
                    .map(|col| tangent.get(row, col) * direction[col])
                    .sum();
                assert!(((plus[row] - minus[row]) / (2.0 * h) - expected).abs() < 1e-9);
            }
            let bound = (0.1_f64 * 2.0 + 1.0) * 2.0;
            let absolute_row: f64 = (0..10).map(|b| exact_entry(0, b).abs()).sum();
            close(
                inertia.force_scale(&u, &[0], 0.0, 0.0, cx).unwrap(),
                bound * absolute_row,
            );
            assert!(absolute_row > 0.05); // |M| bound exceeds the cancelled row sum.
            let (_, energy) = inertia.motion(&u, cx).unwrap();
            close(energy.kinetic_before_j, 1.5);
            close(energy.kinetic_after_j, 0.06);
            close(energy.numerical_loss_j, 0.96);
            close(energy.inertial_work_j, -0.48);
            close(energy.balance_error_j, 0.0);

            // Original positive P1 and full Hex8 rules retain the old
            // mass-action accumulation and precisely the old force-scale bits.
            let p1 = ReferenceMass::new([([0, 1, 2, 3], 1.0 / 6.0)].into_iter(), &[6.0], 12, cx)
                .unwrap();
            let g = 1.0 / 3.0_f64.sqrt();
            let points = core::array::from_fn(|q| {
                let shapes = core::array::from_fn(|a| {
                    (0..3)
                        .map(|i| {
                            let sign = if ((a ^ q) >> i) & 1 == 0 { 1.0 } else { -1.0 };
                            0.5 * (1.0 + sign * g)
                        })
                        .product()
                });
                (shapes, 1.0)
            });
            let hex = ReferenceMass::from_hex8_quadrature(
                [([0, 1, 2, 3, 4, 5, 6, 7], points)].into_iter(),
                &[6.0],
                24,
                cx,
            )
            .unwrap();
            for positive in [p1, hex] {
                let n = positive.ndofs();
                let original = positive
                    .action(&vec![bound; n], cx)
                    .unwrap()
                    .into_iter()
                    .fold(0.0_f64, f64::hypot);
                let (old_u, old_v) = (vec![0.0; n], vec![1.0; n]);
                let positive_inertia =
                    BackwardEulerInertia::new(&positive, &old_u, &old_v, 0.5, cx).unwrap();
                let actual = positive_inertia
                    .force_scale(&vec![0.1; n], &(0..n).collect::<Vec<_>>(), 0.0, 0.0, cx)
                    .unwrap();
                assert_eq!(actual.to_bits(), original.to_bits());
            }
        });
    }

    #[test]
    fn quadratic_mass_refuses_bad_quadrature_budget_and_cancellation() {
        with_cx(false, |cx| {
            let q = quadrature();
            let build = |nodes, points: Vec<([f64; 10], f64)>, rho, cap| {
                ReferenceMass::from_tet10_quadrature(
                    [(nodes, points)].into_iter(),
                    &[rho],
                    30,
                    cap,
                    cx,
                )
            };
            assert!(build(NODES, q.clone(), 6.0, 0).is_err());
            assert!(build(NODES, q.clone(), 0.0, 1 << 20).is_err());
            assert!(build([0; 10], q.clone(), 6.0, 1 << 20).is_err());
            assert!(build(NODES, vec![], 6.0, 1 << 20).is_err());
            for bad in [([0.0; 10], 1.0), ([f64::NAN; 10], 1.0), (q[0].0, -1.0)] {
                assert!(build(NODES, vec![bad], 6.0, 1 << 20).is_err());
            }
        });
        with_cx(true, |cx| {
            assert!(matches!(
                ReferenceMass::from_tet10_quadrature(
                    [(NODES, quadrature())].into_iter(),
                    &[6.0],
                    30,
                    1 << 20,
                    cx,
                ),
                Err(InertiaError::Cancelled)
            ));
        });
    }

    #[test]
    fn midpoint_inertia_matches_constant_acceleration_without_numerical_loss() {
        with_cx(false, |cx| {
            let mass = mass(cx);
            let old_u = [[0.03, -0.02, 0.01]; 10];
            let old_v = [[0.2, -0.1, 0.05]; 10];
            let acceleration = [0.3, -0.2, 0.1];
            let dt = 0.25;
            let next: [[f64; 3]; 10] = core::array::from_fn(|node| {
                core::array::from_fn(|i| {
                    old_u[node][i] + dt * (old_v[node][i] + 0.5 * dt * acceleration[i])
                })
            });
            let inertia =
                MidpointInertia::new(&mass, old_u.as_flattened(), old_v.as_flattened(), dt, cx)
                    .unwrap();
            let mut force = [0.0; 30];
            inertia
                .add_force(next.as_flattened(), &mut force, cx)
                .unwrap();
            let (velocity, energy) = inertia.motion(next.as_flattened(), cx).unwrap();
            for node in 0..10 {
                let share = if node < 4 { -0.05 } else { 0.2 };
                for i in 0..3 {
                    close(force[3 * node + i], share * acceleration[i]);
                    close(velocity[node][i], old_v[node][i] + dt * acceleration[i]);
                }
            }
            let v = core::array::from_fn::<_, 3, _>(|i| old_v[0][i] + dt * acceleration[i]);
            let midpoint =
                core::array::from_fn::<_, 3, _>(|i| old_v[0][i] + 0.5 * dt * acceleration[i]);
            close(energy.kinetic_before_j, 0.5 * dot(&old_v[0], &old_v[0]));
            close(energy.kinetic_after_j, 0.5 * dot(&v, &v));
            close(energy.inertial_work_j, dt * dot(&midpoint, &acceleration));
            close(energy.kinetic_change_j, energy.inertial_work_j);
            assert_eq!(energy.numerical_loss_j.to_bits(), 0.0_f64.to_bits());
            close(energy.balance_error_j, 0.0);
            // Free flight is an exact fixed velocity, with no kinetic decrement.
            let flight = core::array::from_fn::<_, 30, _>(|d| {
                old_u.as_flattened()[d] + dt * old_v.as_flattened()[d]
            });
            let (_, account) = inertia.motion(&flight, cx).unwrap();
            close(account.kinetic_after_j, account.kinetic_before_j);
            close(account.inertial_work_j, 0.0);
        });
    }

    #[test]
    fn midpoint_signed_tangent_and_nonrigid_linear_angular_work_are_consistent() {
        with_cx(false, |cx| {
            let mass = mass(cx);
            let old_u: [f64; 30] = core::array::from_fn(|d| 0.001 * (d + 1) as f64);
            let old_v: [f64; 30] = core::array::from_fn(|d| 0.02 * ((d * 5 % 17) as f64 - 8.0));
            let acceleration: [f64; 30] =
                core::array::from_fn(|d| 0.03 * ((d * 7 % 13) as f64 - 6.0));
            let dt = 0.07;
            let u: [f64; 30] =
                core::array::from_fn(|d| old_u[d] + dt * (old_v[d] + 0.5 * dt * acceleration[d]));
            let inertia = MidpointInertia::new(&mass, &old_u, &old_v, dt, cx).unwrap();
            let mut force = [0.0; 30];
            inertia.add_force(&u, &mut force, cx).unwrap();
            let mut coo = Coo::new(30, 30);
            inertia.add_tangent(&mut coo, cx).unwrap();
            let tangent = coo.assemble();
            assert!(tangent.get(0, 12) < 0.0);
            let direction: [f64; 30] = core::array::from_fn(|d| (d as f64 - 15.0) * 0.01);
            let (mut plus, mut minus) = ([0.0; 30], [0.0; 30]);
            let h = 1e-7;
            inertia
                .add_force(
                    &core::array::from_fn::<_, 30, _>(|d| u[d] + h * direction[d]),
                    &mut plus,
                    cx,
                )
                .unwrap();
            inertia
                .add_force(
                    &core::array::from_fn::<_, 30, _>(|d| u[d] - h * direction[d]),
                    &mut minus,
                    cx,
                )
                .unwrap();
            for row in 0..30 {
                let expected: f64 = (0..30)
                    .map(|col| tangent.get(row, col) * direction[col])
                    .sum();
                assert!(((plus[row] - minus[row]) / (2.0 * h) - expected).abs() < 2e-8);
                for col in 0..30 {
                    close(
                        tangent.get(row, col),
                        2.0 * mass.matrix.get(row, col) / (dt * dt),
                    );
                }
            }
            let (velocity, energy) = inertia.motion(&u, cx).unwrap();
            let next_p = mass.action(velocity.as_flattened(), cx).unwrap();
            let old_p = mass.action(&old_v, cx).unwrap();
            let reference = [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
                [0.5, 0.0, 0.0],
                [0.5, 0.5, 0.0],
                [0.0, 0.5, 0.0],
                [0.0, 0.0, 0.5],
                [0.5, 0.0, 0.5],
                [0.0, 0.5, 0.5],
            ];
            let (mut angular_change, mut midpoint_torque) = ([0.0; 3], [0.0; 3]);
            for axis in 0..3 {
                let change: f64 = (0..10)
                    .map(|a| next_p[3 * a + axis] - old_p[3 * a + axis])
                    .sum();
                let impulse: f64 = (0..10).map(|a| dt * force[3 * a + axis]).sum();
                close(change, impulse);
                let (j, k) = ((axis + 1) % 3, (axis + 2) % 3);
                for (a, x) in reference.iter().enumerate() {
                    let (old_j, old_k) = (x[j] + old_u[3 * a + j], x[k] + old_u[3 * a + k]);
                    let (next_j, next_k) = (x[j] + u[3 * a + j], x[k] + u[3 * a + k]);
                    angular_change[axis] += next_j * next_p[3 * a + k]
                        - next_k * next_p[3 * a + j]
                        - old_j * old_p[3 * a + k]
                        + old_k * old_p[3 * a + j];
                    midpoint_torque[axis] += 0.5 * (old_j + next_j) * force[3 * a + k]
                        - 0.5 * (old_k + next_k) * force[3 * a + j];
                }
                close(angular_change[axis], dt * midpoint_torque[axis]);
            }
            close(energy.kinetic_change_j, energy.inertial_work_j);
            assert_eq!(energy.numerical_loss_j, 0.0);
            let expected_bound: f64 = (0..30)
                .map(|d| {
                    2.0 * mass.matrix.get(0, d).abs()
                        * ((u[d] - old_u[d]).abs() / dt + old_v[d].abs())
                        / dt
                })
                .sum();
            close(
                inertia.force_scale(&u, &[0], 0.0, 0.0, cx).unwrap(),
                expected_bound,
            );
        });
    }

    #[test]
    fn midpoint_inertia_refuses_invalid_intervals_states_and_cancellation() {
        with_cx(false, |cx| {
            let mass = mass(cx);
            let zero = [0.0; 30];
            for dt in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::MIN_POSITIVE] {
                assert!(MidpointInertia::new(&mass, &zero, &zero, dt, cx).is_err());
            }
            assert!(MidpointInertia::new(&mass, &zero[..29], &zero, 0.1, cx).is_err());
            assert!(MidpointInertia::new(&mass, &zero, &[f64::NAN; 30], 0.1, cx).is_err());
            let inertia = MidpointInertia::new(&mass, &zero, &zero, 0.1, cx).unwrap();
            assert!(inertia.motion(&zero[..29], cx).is_err());
            assert!(inertia.motion(&[f64::NAN; 30], cx).is_err());
            assert!(inertia.add_force(&zero, &mut [0.0; 29], cx).is_err());
            with_cx(true, |cancelled| {
                assert!(matches!(
                    MidpointInertia::new(&mass, &zero, &zero, 0.1, cancelled),
                    Err(InertiaError::Cancelled)
                ));
                assert!(matches!(
                    inertia.motion(&zero, cancelled),
                    Err(InertiaError::Cancelled)
                ));
                assert!(matches!(
                    inertia.add_force(&zero, &mut [0.0; 30], cancelled),
                    Err(InertiaError::Cancelled)
                ));
                assert!(matches!(
                    inertia.add_tangent(&mut Coo::new(30, 30), cancelled),
                    Err(InertiaError::Cancelled)
                ));
            });
        });
    }
}
