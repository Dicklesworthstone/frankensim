//! Canonical momentum and kinetic energy for the original tensor-material owner.
use super::*;

pub(super) enum Inertia {
    Diagonal(Vec<f64>),
    Consistent { dimension: usize, mass: Vec<f64>, inverse: Vec<f64> },
}
impl Inertia {
    pub(super) fn diagonal(masses: Vec<f64>) -> Result<Self> {
        if masses.iter().any(|m| !m.is_finite() || *m <= 0.0 || !m.recip().is_finite()) {
            return Err(bad("finite positive diagonal inertia required"));
        }
        Ok(Self::Diagonal(masses))
    }
    pub(super) fn consistent(n: usize, mass: Vec<f64>, gate: &CancelGate) -> Result<Self> {
        if mass.len() != n*n || !finite(&mass) {
            return Err(bad("complete finite consistent mass matrix required"));
        }
        for i in 0..n { for j in 0..i {
            if mass[i*n+j] != mass[j*n+i] {
                return Err(bad("consistent mass matrix must be exactly symmetric"));
            }
        }}
        poll(gate)?;
        let factor = fs_la::factor::cholesky(&mass, n).map_err(TensorMaterialError::Mass)?;
        let mut inverse = vec![0.0; n*n];
        for j in 0..n {
            poll(gate)?;
            let mut column = vec![0.0; n]; column[j] = 1.0;
            factor.solve(&mut column);
            if !finite(&column) { return Err(bad("consistent inverse inertia is unrepresentable")); }
            // An SPD inverse has one independent triangle. Use those original
            // Cholesky solves, reflected exactly; never average or lump mass.
            for i in 0..=j { inverse[i*n+j] = column[i]; inverse[j*n+i] = column[i]; }
        }
        // A dependent basis can leave a tiny positive floating-point pivot.
        // Admit only an inverse that resolves the identity at this dimension;
        // never regularize the physical mass or silently discard coordinates.
        let tolerance = 64.0 * (n as f64) * f64::EPSILON;
        for i in 0..n {
            poll(gate)?;
            for j in 0..n {
                let product: f64 = (0..n)
                    .map(|k| mass[i*n+k] * inverse[k*n+j]).sum();
                let expected = if i == j { 1.0 } else { 0.0 };
                if !product.is_finite() || (product-expected).abs() > tolerance {
                    return Err(bad("consistent mass inverse does not resolve the original inertia"));
                }
            }
        }
        Ok(Self::Consistent { dimension: n, mass, inverse })
    }
    pub(super) fn dimension(&self) -> usize {
        match self { Self::Diagonal(m) => m.len(), Self::Consistent { dimension, .. } => *dimension }
    }
    pub(super) fn inverse_entry(&self, i: usize, j: usize) -> f64 {
        match self {
            Self::Diagonal(m) => if i == j { m[i].recip() } else { 0.0 },
            Self::Consistent { dimension, inverse, .. } => inverse[i*dimension+j],
        }
    }
    pub(super) fn momentum(&self, velocity: &[f64]) -> Vec<f64> {
        match self {
            Self::Diagonal(m) => m.iter().zip(velocity).map(|(m,v)| m*v).collect(),
            Self::Consistent { dimension, mass, .. } => mass.chunks_exact(*dimension)
                .map(|row| row.iter().zip(velocity).map(|(m,v)| m*v).sum()).collect(),
        }
    }
    pub(super) fn velocity(&self, momentum: &[f64]) -> Vec<f64> {
        match self {
            Self::Diagonal(m) => momentum.iter().zip(m).map(|(p,m)| p/m).collect(),
            Self::Consistent { dimension, inverse, .. } => inverse.chunks_exact(*dimension)
                .map(|row| row.iter().zip(momentum).map(|(m,p)| m*p).sum()).collect(),
        }
    }
    pub(super) fn matrix(&self) -> Vec<f64> {
        match self {
            Self::Diagonal(m) => {
                let n=m.len(); let mut matrix=vec![0.0;n*n];
                for i in 0..n {matrix[i*n+i]=m[i];} matrix
            }
            Self::Consistent { mass, .. } => mass.clone(),
        }
    }
}
