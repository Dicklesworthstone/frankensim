//! Reference-configuration inertia for the existing frame equations.
//!
//! Linear axial interpolation and cubic Hermite transverse interpolation
//! integrate centerline translational kinetic energy. Nodal rotations carry
//! generalized inertia through that interpolation; cross-section rotary
//! inertia and configuration-dependent mass are not introduced.

use super::{Frame2d, SolidError};
use fs_la::factor::{Cholesky, cholesky};

/// Element inertia selected when the frame is constructed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FrameMassFormulation {
    /// Half the member mass at each endpoint's translational DOFs.
    Lumped,
    /// Consistent Euler–Bernoulli mass in the undeformed reference frame.
    /// Appropriate to linear or moderate-rotation vibration; it does not
    /// provide geometrically exact moving-frame inertia.
    ConsistentReference,
}

fn invalid(what: &str) -> SolidError {
    SolidError::InvalidInput { what: what.into() }
}

/// The diagonal path retains the original nodal/lumped arithmetic and storage.
pub(super) enum MassMatrix {
    Lumped(Vec<f64>),
    Consistent { n: usize, entries: Vec<f64> },
}

/// Retained mass factor for restoring all eigenvectors without refactorization.
pub(super) enum MassBackTransform {
    Diagonal(Vec<f64>),
    Cholesky(Cholesky),
}

impl MassBackTransform {
    pub(super) fn restore(&self, vector: &[f64]) -> Vec<f64> {
        match self {
            Self::Diagonal(roots) => vector.iter().zip(roots).map(|(v, root)| v / root).collect(),
            Self::Cholesky(l) => {
                let mut result = vector.to_vec();
                for i in (0..result.len()).rev() {
                    let mut value = result[i];
                    for j in i + 1..result.len() {
                        value -= l.l(j, i) * result[j];
                    }
                    result[i] = value / l.l(i, i);
                }
                result
            }
        }
    }
}

impl MassMatrix {
    fn n(&self) -> usize {
        match self {
            Self::Lumped(m) => m.len(),
            Self::Consistent { n, .. } => *n,
        }
    }

    pub(super) fn entry(&self, i: usize, j: usize) -> f64 {
        match self {
            Self::Lumped(m) => {
                if i == j {
                    m[i]
                } else {
                    0.0
                }
            }
            Self::Consistent { n, entries } => entries[i * n + j],
        }
    }

    pub(super) fn row_action(&self, i: usize, x: &[f64]) -> f64 {
        match self {
            Self::Lumped(m) => m[i] * x[i],
            Self::Consistent { n, entries } => entries[i * n..(i + 1) * n]
                .iter()
                .zip(x)
                .map(|(m, x)| m * x)
                .sum(),
        }
    }

    pub(super) fn kinetic_energy(&self, v: &[f64]) -> f64 {
        0.5 * match self {
            Self::Lumped(m) => (0..m.len()).map(|i| m[i] * v[i] * v[i]).sum::<f64>(),
            Self::Consistent { n, .. } => (0..*n).map(|i| v[i] * self.row_action(i, v)).sum(),
        }
    }

    pub(super) fn add_scaled(&self, matrix: &mut [f64], coefficient: f64) {
        match self {
            Self::Lumped(m) => {
                for i in 0..m.len() {
                    matrix[i * m.len() + i] += coefficient * m[i];
                }
            }
            Self::Consistent { entries, .. } => {
                for (a, m) in matrix.iter_mut().zip(entries) {
                    *a += coefficient * m;
                }
            }
        }
    }

    fn factor(&self, dofs: &[usize]) -> Result<Cholesky, SolidError> {
        let n = dofs.len();
        let mut matrix = vec![0.0; n * n];
        for i in 0..n {
            for j in 0..n {
                matrix[i * n + j] = self.entry(dofs[i], dofs[j]);
            }
        }
        cholesky(&matrix, n)
            .map_err(|_| invalid("reference inertia is not positive definite on massive free DOFs"))
    }

    pub(super) fn acceleration(
        &self,
        force: &[f64],
        fixed: &[bool],
    ) -> Result<Vec<f64>, SolidError> {
        let mut acceleration = vec![0.0; self.n()];
        match self {
            Self::Lumped(m) => {
                for i in 0..m.len() {
                    if !fixed[i] && m[i] > 0.0 {
                        acceleration[i] = force[i] / m[i];
                    }
                }
            }
            Self::Consistent { .. } => {
                let dofs: Vec<usize> = (0..self.n())
                    .filter(|&i| !fixed[i] && self.entry(i, i) > 0.0)
                    .collect();
                let mut rhs: Vec<f64> = dofs.iter().map(|&i| force[i]).collect();
                self.factor(&dofs)?.solve(&mut rhs);
                for (&i, value) in dofs.iter().zip(rhs) {
                    acceleration[i] = value;
                }
            }
        }
        if acceleration.iter().any(|a| !a.is_finite()) {
            return Err(invalid("nonfinite initial acceleration"));
        }
        Ok(acceleration)
    }

    /// L^-1 K L^-T with its retained inverse-transpose mass transform.
    pub(super) fn modal_matrix(
        &self,
        k: &[f64],
        dofs: &[usize],
    ) -> Result<(Vec<Vec<f64>>, MassBackTransform), SolidError> {
        let n = dofs.len();
        let mut a = vec![vec![0.0; n]; n];
        for i in 0..n {
            for j in 0..n {
                a[i][j] = 0.5 * (k[i * n + j] + k[j * n + i]);
            }
        }
        let transform = match self {
            Self::Lumped(m) => {
                let roots: Vec<f64> = dofs.iter().map(|&d| m[d].sqrt()).collect();
                for i in 0..n {
                    for j in 0..n {
                        a[i][j] /= roots[i] * roots[j];
                    }
                }
                MassBackTransform::Diagonal(roots)
            }
            Self::Consistent { .. } => {
                let l = self.factor(dofs)?;
                for i in 0..n {
                    for j in 0..n {
                        let mut value = a[i][j];
                        for p in 0..i {
                            value -= l.l(i, p) * a[p][j];
                        }
                        a[i][j] = value / l.l(i, i);
                    }
                }
                for row in &mut a {
                    for j in 0..n {
                        let mut value = row[j];
                        for p in 0..j {
                            value -= l.l(j, p) * row[p];
                        }
                        row[j] = value / l.l(j, j);
                    }
                }
                for i in 0..n {
                    for j in 0..i {
                        let symmetric = 0.5 * (a[i][j] + a[j][i]);
                        a[i][j] = symmetric;
                        a[j][i] = symmetric;
                    }
                }
                MassBackTransform::Cholesky(l)
            }
        };
        if a.iter().flatten().any(|x| !x.is_finite()) {
            return Err(invalid("nonfinite mass-reduced vibration pencil"));
        }
        Ok((a, transform))
    }
}

impl Frame2d {
    /// Work-conjugate reference transverse trace at member coordinate xi in [0,1].
    ///
    /// Returns original node-major DOFs, their six displacement weights and
    /// reference member length [m]. The dot product with displacement gives
    /// transverse displacement [m]; with velocity it gives m/s. Multiplying
    /// these SAME weights by a transverse point force returns nodal forces and
    /// moments, so virtual work is unchanged. Rotation weights have metres.
    ///
    /// This is the cubic Hermite interpolation of `ConsistentReference` inertia,
    /// along the member's reference normal (-sin,cos). The selected dry inertia
    /// formulation is not changed. Quadrature weights/areas belong to the caller;
    /// this method supplies neither fluid forces nor a material law. Only linear
    /// reference geometry is admitted; finite-rotation follower loads need their
    /// own current-geometry trace. No allocation or accepted-state change occurs.
    pub fn reference_transverse_trace(
        &self,
        element: usize,
        xi: f64,
    ) -> Result<([usize; 6], [f64; 6], f64), SolidError> {
        if self.geometry != super::Geometry::Linear || !xi.is_finite()
            || !(0.0..=1.0).contains(&xi)
        {
            return Err(invalid("transverse trace requires linear reference geometry and xi in [0,1]"));
        }
        let e = self.elements.get(element)
            .ok_or_else(|| invalid("transverse trace member does not exist"))?;
        if !e.l0.is_finite() || e.l0 <= 0.0 || !e.cos0.is_finite() || !e.sin0.is_finite() {
            return Err(invalid("transverse trace reference geometry is not finite"));
        }
        let one_minus = 1.0 - xi;
        let h0 = one_minus * one_minus * (1.0 + 2.0 * xi);
        let h1 = e.l0 * xi * one_minus * one_minus;
        let h2 = xi * xi * (3.0 - 2.0 * xi);
        let h3 = -e.l0 * xi * xi * one_minus;
        let weights = [-e.sin0*h0, e.cos0*h0, h1, -e.sin0*h2, e.cos0*h2, h3];
        if weights.iter().any(|x| !x.is_finite()) {
            return Err(invalid("transverse trace weights are unrepresentable"));
        }
        Ok(([3*e.i, 3*e.i+1, 3*e.i+2, 3*e.j, 3*e.j+1, 3*e.j+2], weights, e.l0))
    }

    /// Generalized reference inertial force M*a for one acceleration per DOF.
    /// Translational entries are m/s² and angular entries rad/s²; the output
    /// contains forces and moments. The same mass owner drives modal extraction
    /// and Newmark. A uniform gravity acceleration therefore supplies its
    /// consistent nodal body load, including member end moments.
    ///
    /// # Errors
    /// Invalid/nonfinite acceleration or unrepresentable reference inertia.
    pub fn inertial_forces(&self, acceleration: &[f64]) -> Result<Vec<f64>, SolidError> {
        if acceleration.len() != self.ndof() || acceleration.iter().any(|v| !v.is_finite()) {
            return Err(invalid(
                "inertial acceleration needs one finite value per frame DOF",
            ));
        }
        let mass = self.reference_mass()?;
        let force: Vec<f64> = (0..self.ndof())
            .map(|i| mass.row_action(i, acceleration))
            .collect();
        if force.iter().any(|v| !v.is_finite()) {
            return Err(invalid("nonfinite reference inertial force"));
        }
        Ok(force)
    }

    pub(super) fn reference_mass(&self) -> Result<MassMatrix, SolidError> {
        if self
            .nodal_mass
            .iter()
            .chain(&self.point_mass)
            .any(|v| !v.is_finite() || *v < 0.0)
        {
            return Err(invalid("frame masses must be finite and nonnegative"));
        }
        if self.mass_formulation == FrameMassFormulation::Lumped
            || self.elements.iter().all(|e| e.mass_per_length == 0.0)
        {
            return Ok(MassMatrix::Lumped(self.nodal_mass.clone()));
        }
        let n = self.ndof();
        let count = n
            .checked_mul(n)
            .ok_or_else(|| invalid("reference mass matrix size overflow"))?;
        let mut entries = vec![0.0; count];
        for i in 0..n {
            entries[i * n + i] = self.point_mass[i];
        }
        for e in &self.elements {
            if e.mass_per_length == 0.0 {
                continue;
            }
            let l = e.l0;
            let member_mass = e.mass_per_length * l;
            let mut local = [[0.0; 6]; 6];
            for (i, j, factor) in [(0, 0, 2.0), (0, 3, 1.0), (3, 0, 1.0), (3, 3, 2.0)] {
                local[i][j] = member_mass * factor / 6.0;
            }
            let bending = [
                [156.0, 22.0 * l, 54.0, -13.0 * l],
                [22.0 * l, 4.0 * l * l, 13.0 * l, -3.0 * l * l],
                [54.0, 13.0 * l, 156.0, -22.0 * l],
                [-13.0 * l, -3.0 * l * l, -22.0 * l, 4.0 * l * l],
            ];
            let transverse = [1, 2, 4, 5];
            for i in 0..4 {
                for j in 0..4 {
                    local[transverse[i]][transverse[j]] = member_mass * bending[i][j] / 420.0;
                }
            }
            let mut transform = [[0.0; 6]; 6];
            for b in [0, 3] {
                transform[b][b] = e.cos0;
                transform[b][b + 1] = e.sin0;
                transform[b + 1][b] = -e.sin0;
                transform[b + 1][b + 1] = e.cos0;
                transform[b + 2][b + 2] = 1.0;
            }
            let mut product = [[0.0; 6]; 6];
            for i in 0..6 {
                for j in 0..6 {
                    product[i][j] = (0..6).map(|p| local[i][p] * transform[p][j]).sum();
                }
            }
            let dofs = [
                3 * e.i,
                3 * e.i + 1,
                3 * e.i + 2,
                3 * e.j,
                3 * e.j + 1,
                3 * e.j + 2,
            ];
            for i in 0..6 {
                for j in 0..=i {
                    let value: f64 = (0..6).map(|p| transform[p][i] * product[p][j]).sum();
                    entries[dofs[i] * n + dofs[j]] += value;
                    if i != j {
                        entries[dofs[j] * n + dofs[i]] += value;
                    }
                }
            }
        }
        if entries.iter().any(|v| !v.is_finite()) {
            return Err(invalid("nonfinite consistent reference mass"));
        }
        Ok(MassMatrix::Consistent { n, entries })
    }
}

#[cfg(test)]
mod transverse_trace_tests {
    use super::*;
    use crate::frame2d::{Geometry, SectionModel};

    fn frame(length: f64, angle: f64) -> Frame2d {
        let mut frame = Frame2d::with_mass_formulation(Geometry::Linear, FrameMassFormulation::ConsistentReference);
        frame.add_node(0.0, 0.0);
        frame.add_node(length*angle.cos(), length*angle.sin());
        frame.add_element(0, 1, 3, 3.4, &|| SectionModel::Elastic { ea: 200.0, ei: 8.0 }).unwrap();
        frame
    }
    fn near(a: f64, b: f64) {
        assert!((a-b).abs() < 2e-12*(1.0+a.abs()+b.abs()), "{a:e} != {b:e}");
    }

    #[test]
    fn trace_endpoints_and_rotation_weights_use_original_dofs() {
        let frame = frame(2.0, 0.0);
        let (dofs, left, length) = frame.reference_transverse_trace(0, 0.0).unwrap();
        let (_, right, _) = frame.reference_transverse_trace(0, 1.0).unwrap();
        let (_, mid, _) = frame.reference_transverse_trace(0, 0.5).unwrap();
        assert_eq!(dofs, [0,1,2,3,4,5]);
        assert_eq!(left, [0.0,1.0,0.0,0.0,0.0,0.0]);
        assert_eq!(right, [0.0,0.0,0.0,0.0,1.0,0.0]);
        assert_eq!(mid, [0.0,0.5,0.25,0.0,0.5,-0.25]);
        assert_eq!(length, 2.0);
    }

    #[test]
    fn rigid_rotation_is_reproduced_on_rotated_members() {
        for angle in [0.0_f64, 0.37, 1.2] {
            let length = 2.7;
            let frame = frame(length, angle);
            let theta = 0.012;
            let velocity = [0.0, 0.0, theta, -angle.sin()*theta*length, angle.cos()*theta*length, theta];
            for xi in [0.0, 0.13, 0.5, 0.91, 1.0] {
                let (_, weights, _) = frame.reference_transverse_trace(0, xi).unwrap();
                near(weights.iter().zip(velocity).map(|(b,v)| b*v).sum(), theta*length*xi);
            }
        }
    }

    #[test]
    fn four_point_trace_quadrature_recovers_original_consistent_transverse_inertia() {
        let points = [0.06943184420297371, 0.33000947820757187, 0.6699905217924281, 0.9305681557970262];
        let weights = [0.17392742256872693, 0.32607257743127307, 0.32607257743127307, 0.17392742256872693];
        for angle in [0.0_f64, 0.37, 1.2] {
            for length in [0.3, 1.0, 2.7] {
                let frame = frame(length, angle);
                // Pure transverse translations plus rotations: no axial kinetic energy.
                let velocity = [-angle.sin()*0.3, angle.cos()*0.3, -0.2,
                    angle.sin()*0.4, -angle.cos()*0.4, 0.17];
                let mut energy = 0.0;
                for (xi, weight) in points.into_iter().zip(weights) {
                    let (_, b, l) = frame.reference_transverse_trace(0, xi).unwrap();
                    let speed: f64 = b.iter().zip(velocity).map(|(b,v)| b*v).sum();
                    energy += 0.5*3.4*l*weight*speed*speed;
                }
                near(energy, frame.reference_mass().unwrap().kinetic_energy(&velocity));
            }
        }
    }

    #[test]
    fn one_physical_trace_preserves_transverse_force_and_nodal_moment_work() {
        let frame = frame(2.0, 0.7);
        let (_, b, _) = frame.reference_transverse_trace(0, 0.3).unwrap();
        let velocity = [0.1, -0.3, 0.2, 0.7, -0.1, 0.4];
        let force = -2.3;
        let point_power = force*b.iter().zip(velocity).map(|(b,v)| b*v).sum::<f64>();
        let nodal_power: f64 = b.iter().zip(velocity).map(|(b,v)| (force*b)*v).sum();
        near(point_power, nodal_power);
        assert_ne!(b[2], 0.0);
        assert_ne!(b[5], 0.0);
    }

    #[test]
    fn unsupported_trace_requests_refuse_without_changing_the_frame() {
        let mut frame = frame(1.0, 0.0);
        let saved = frame.displacements().to_vec();
        assert!(frame.reference_transverse_trace(1, 0.5).is_err());
        for xi in [-0.1, 1.1, f64::NAN, f64::INFINITY] {
            assert!(frame.reference_transverse_trace(0, xi).is_err());
        }
        assert_eq!(frame.displacements(), saved);
        frame.geometry = Geometry::Corotational;
        assert!(frame.reference_transverse_trace(0, 0.5).is_err());
    }
}
