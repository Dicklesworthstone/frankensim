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

    /// L^-1 K L^-T, using the existing Cholesky owner. The caller retains
    /// its existing symmetric eigensolver and massless stiffness condensation.
    pub(super) fn modal_matrix(
        &self,
        k: &[f64],
        dofs: &[usize],
    ) -> Result<Vec<Vec<f64>>, SolidError> {
        let n = dofs.len();
        let mut a = vec![vec![0.0; n]; n];
        for i in 0..n {
            for j in 0..n {
                a[i][j] = 0.5 * (k[i * n + j] + k[j * n + i]);
            }
        }
        match self {
            Self::Lumped(m) => {
                for i in 0..n {
                    for j in 0..n {
                        a[i][j] /= m[dofs[i]].sqrt() * m[dofs[j]].sqrt();
                    }
                }
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
            }
        }
        if a.iter().flatten().any(|x| !x.is_finite()) {
            return Err(invalid("nonfinite mass-reduced vibration pencil"));
        }
        Ok(a)
    }
}

impl Frame2d {
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
