//! Resolve zero-stiffness coordinates using ONLY declared bilateral springs.
//!
//! Write C = B sqrt(k), split elastic/free coordinates, and retain the existing
//! A = I + C_e^T D_e C_e factor. The free-coordinate Schur matrix is
//! S = C_f A^-1 C_f^T. Its size cannot exceed the admitted connection count.
//! No full modal stiffness, artificial stiffness, pseudoinverse or fixed pose.
use super::*;

pub(super) struct FreeResponse {
    indices: Vec<usize>,
    rows: Vec<Vec<f64>>,
    scaled_matrix: Vec<f64>,
    scales: Vec<f64>,
    factor: Cholesky,
}

impl FreeResponse {
    // Includes the existing setup screen plus bounded free-block solves,
    // contractions and factor work. All dimensions have already been admitted.
    pub(super) fn admit(network: &CoupledModalSystem, free: usize)
        -> Result<(), ModalCouplingError>
    {
        if free == 0 { return Ok(()); }
        let links = network.columns.len();
        if free > links {
            return Err(invalid("bilateral springs do not constrain every free coordinate"));
        }
        let n = network.mode_count();
        let terms = n.checked_mul(links + 1).and_then(|v| v.checked_mul(links + 1))
            .and_then(|v| v.checked_add(free * links * links))
            .and_then(|v| v.checked_add(free * free * links))
            .and_then(|v| v.checked_add(free * free * free))
            .ok_or_else(|| invalid("supported-free-coordinate preload setup overflow"))?;
        if terms > network.config.max_setup_terms {
            return Err(invalid("supported-free-coordinate preload exceeds max_setup_terms"));
        }
        Ok(())
    }

    pub(super) fn new(network: &CoupledModalSystem, indices: Vec<usize>, roots: &[f64],
        matrix: &[f64], factor: &Cholesky, gate: &CancelGate)
        -> Result<Self, ModalCouplingError>
    {
        let n = indices.len();
        let mut rows = Vec::with_capacity(n);
        let mut responses = Vec::with_capacity(n);
        for &index in &indices {
            poll(Some(gate))?;
            let row: Vec<f64> = network.columns.iter().zip(roots)
                .map(|(b, root)| finite(b[index] * root)).collect::<Result<_, _>>()?;
            let mut response = row.clone();
            factor.solve(&mut response);
            check_solve(matrix, &response, &row, network.config.solve_relative_tolerance)?;
            rows.push(row);
            responses.push(response);
        }
        let mut scales = Vec::with_capacity(n);
        for i in 0..n {
            let diagonal = dot(&rows[i], &responses[i])?;
            if diagonal <= 0.0 {
                return Err(invalid("free coordinate has no positive static spring support"));
            }
            let scale = finite(1.0 / diagonal.sqrt())?;
            if scale <= 0.0 { return Err(invalid("free support equilibration is not representable")); }
            scales.push(scale);
        }
        let mut scaled_matrix = vec![0.0; n*n];
        for i in 0..n {
            poll(Some(gate))?;
            for j in 0..=i {
                let a = finite(finite(dot(&rows[i], &responses[j])? * scales[i])? * scales[j])?;
                let b = finite(finite(dot(&rows[j], &responses[i])? * scales[j])? * scales[i])?;
                if (a-b).abs() > network.config.solve_relative_tolerance * a.abs().max(b.abs()).max(1.0) {
                    return Err(invalid("free support response lost numerical symmetry"));
                }
                let entry = f64::midpoint(a,b);
                scaled_matrix[i*n+j] = entry;
                scaled_matrix[j*n+i] = entry;
            }
        }
        let factor = cholesky(&scaled_matrix, n).map_err(ModalCouplingError::Factor)?;
        // A rounded singular Gram matrix can have a tiny positive pivot. Refuse
        // unresolved support directions rather than selecting a gauge. This is
        // a floating-point degeneracy screen, not a spectral/rank certificate.
        let floor = 256.0 * f64::EPSILON * n as f64;
        if (0..n).any(|i| factor.l(i,i).powi(2) <= floor) {
            return Err(invalid("free support directions are linearly dependent or unresolved at roundoff"));
        }
        poll(Some(gate))?;
        Ok(Self { indices, rows, scaled_matrix, scales, factor })
    }

    // The existing connection solution is A^-1(C_e^T D_e f_e - sqrt(k) rest).
    // Solve S q_f = f_f - C_f solution, then rebuild the FULL connection RHS.
    pub(super) fn complete(&self, external: &[f64], rhs: &mut [f64], solution: &[f64],
        tolerance: f64, gate: &CancelGate) -> Result<Vec<f64>, ModalCouplingError>
    {
        let mut scaled_rhs = Vec::with_capacity(self.indices.len());
        for (i, &index) in self.indices.iter().enumerate() {
            poll(Some(gate))?;
            scaled_rhs.push(finite(finite(external[index] - dot(&self.rows[i], solution)?)? * self.scales[i])?);
        }
        let mut scaled_q = scaled_rhs.clone();
        self.factor.solve(&mut scaled_q);
        check_solve(&self.scaled_matrix, &scaled_q, &scaled_rhs, tolerance)?;
        let q: Vec<f64> = scaled_q.iter().zip(&self.scales)
            .map(|(q,s)| finite(q*s)).collect::<Result<_,_>>()?;
        for (row, &value) in self.rows.iter().zip(&q) {
            poll(Some(gate))?;
            for (r,b) in rhs.iter_mut().zip(row) { *r = finite(*r + b*value)?; }
        }
        Ok(q)
    }

    // This gate is also applied to unit-load compliance queries for contacts;
    // a small Schur residual alone is not a complete free-mode force balance.
    pub(super) fn check_and_insert(&self, external: &[f64], solution: &[f64], free_q: &[f64],
        q: &mut [f64], tolerance: f64, gate: &CancelGate) -> Result<(), ModalCouplingError>
    {
        for (i, &index) in self.indices.iter().enumerate() {
            poll(Some(gate))?;
            let applied = dot(&self.rows[i], solution)?;
            let mut scale = external[index].abs();
            for (a,b) in self.rows[i].iter().zip(solution) { scale = finite(scale + finite(a*b)?.abs())?; }
            let residual = finite(applied - external[index])?.abs();
            let relative = if scale == 0.0 { 0.0 } else { residual / scale };
            if relative > tolerance {
                return Err(ModalCouplingError::SolveResidual { relative, tolerance });
            }
            q[index] = free_q[i];
        }
        Ok(())
    }
}
