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
            let rho = nodes.map(|node| densities[node]);
            let uniform = rho.iter().all(|value| value.to_bits() == rho[0].to_bits());
            let scale = rho.into_iter().fold(0.0_f64, f64::max);
            let normalized = rho.map(|value| value / scale);
            let total: f64 = normalized.iter().sum();
            let mean = if uniform {
                rho[0]
            } else {
                scale * (0.25 * total)
            };
            let element_mass = finite(volume * mean)?;
            if element_mass <= 0.0 {
                return Err(invalid(
                    "each tetrahedron needs representable positive reference mass",
                ));
            }
            means.push(mean);
            reference_mass_kg = finite(reference_mass_kg + element_mass)?;
            let fractions = normalized.map(|value| value / total);
            // Integral rho*N_a*N_b has degree three. Its off-diagonal entry
            // is V*(sum(rho)+rho_a+rho_b)/120, and its diagonal entry is
            // V*(sum(rho)+2*rho_a)/60. Normalization avoids overflowing the
            // density sum; every entry remains strictly nonnegative.
            let uniform_mass =
                uniform.then(|| crate::linear3::p1_tet_consistent_mass(volume, rho[0]));
            for a in 0..4 {
                for b in 0..4 {
                    let value = finite(if let Some(local) = &uniform_mass {
                        local[3 * a][3 * b]
                    } else if a == b {
                        (element_mass / 15.0) * (1.0 + 2.0 * fractions[a])
                    } else {
                        (element_mass / 30.0) * (1.0 + (fractions[a] + fractions[b]))
                    })?;
                    if value <= 0.0 {
                        return Err(invalid("consistent graded-density mass underflowed"));
                    }
                    for component in 0..3 {
                        matrix.push(3 * nodes[a] + component, 3 * nodes[b] + component, value);
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
