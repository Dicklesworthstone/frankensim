//! Resolve zero-stiffness coordinates using declared bilateral springs.
//!
//! Write C = B sqrt(k), split elastic/free coordinates, and retain the existing
//! A = I + C_e^T D_e C_e factor. The free-coordinate Schur matrix is
//! S = C_f A^-1 C_f^T. An explicit admitted null basis N permits compatible
//! free assemblies: (S + sigma N N^T) q_f = b selects N^T q_f = 0.
//! The augmentation is only an algebraic gauge. Original force balance is
//! checked independently; no physical stiffness or source projection is added.
use super::*;

pub(super) struct FreeResponse {
    indices: Vec<usize>,
    rows: Vec<Vec<f64>>,
    scaled_matrix: Vec<f64>,
    scales: Vec<f64>,
    factor: Cholesky,
    nullspace: Vec<Vec<f64>>,
}

impl FreeResponse {
    // Includes the existing setup screen plus bounded free-block solves,
    // contractions, declared-nullspace checks and factor work.
    pub(super) fn admit(network: &StaticSystem<'_>, free: usize, gauges: usize)
        -> Result<(), ModalCouplingError>
    {
        if gauges > free {
            return Err(invalid("more static gauges than free coordinates"));
        }
        if free == 0 { return Ok(()); }
        let links = network.columns.len();
        if free > links + gauges {
            return Err(invalid("bilateral springs and declared gauges do not span every free coordinate"));
        }
        let (n,k,f,g) = (network.mode_count() as u128, links as u128,
            free as u128, gauges as u128);
        let terms = n*(k+1)*(k+1) + f*k*k + f*f*k + f*f*f
            + n*g*(k+g+1) + f*f*g + 2*f*g + if gauges == 0 { 0 } else { f*k };
        if terms > network.config.max_setup_terms as u128 {
            return Err(invalid("supported-free-coordinate preload exceeds max_setup_terms"));
        }
        Ok(())
    }

    pub(super) fn new(network: &StaticSystem<'_>, indices: Vec<usize>, roots: &[f64],
        matrix: &[f64], factor: &Cholesky, nullspace: &[Vec<f64>], gate: &CancelGate)
        -> Result<Self, ModalCouplingError>
    {
        let n = indices.len();
        let nullspace = admit_nullspace(network, &indices, nullspace, gate)?;
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
        let diagonals: Vec<f64> = (0..n).map(|i| dot(&rows[i], &responses[i]))
            .collect::<Result<_,_>>()?;
        if diagonals.iter().any(|d| *d < 0.0) {
            return Err(invalid("free support Schur diagonal is negative"));
        }
        // Matching the physical Schur scale keeps the gauge from drowning a
        // weak or stiff assembly. A completely free system has no such scale;
        // unit augmentation only selects its declared zero pose under zero load.
        let weight = if nullspace.is_empty() { 0.0 } else {
            let maximum = diagonals.iter().copied().fold(0.0_f64, f64::max);
            if maximum == 0.0 { 1.0 } else { maximum }
        };
        let mut scales = Vec::with_capacity(n);
        for (i, &physical) in diagonals.iter().enumerate() {
            let gauge = nullspace.iter().try_fold(0.0, |sum, v| finite(sum + v[i]*v[i]))?;
            let diagonal = finite(physical + weight*gauge)?;
            if diagonal <= 0.0 {
                return Err(invalid("free coordinate has no positive static spring support or declared gauge"));
            }
            let scale = finite(1.0 / diagonal.sqrt())?;
            if scale <= 0.0 { return Err(invalid("free support equilibration is not representable")); }
            scales.push(scale);
        }
        let mut scaled_matrix = vec![0.0; n*n];
        for i in 0..n {
            poll(Some(gate))?;
            for j in 0..=i {
                let gauge = nullspace.iter().try_fold(0.0, |sum, v| finite(sum + v[i]*v[j]))?;
                let a = finite(finite(finite(dot(&rows[i], &responses[j])? + weight*gauge)? * scales[i])? * scales[j])?;
                let b = finite(finite(finite(dot(&rows[j], &responses[i])? + weight*gauge)? * scales[j])? * scales[i])?;
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
        // unresolved support directions or missing gauges instead of selecting
        // a new pose. This is a roundoff screen, not a rank certificate.
        let floor = 256.0 * f64::EPSILON * n as f64;
        if (0..n).any(|i| factor.l(i,i).powi(2) <= floor) {
            return Err(invalid("free support directions are linearly dependent or unresolved at roundoff"));
        }
        poll(Some(gate))?;
        Ok(Self { indices, rows, scaled_matrix, scales, factor, nullspace })
    }

    // The existing connection solution is A^-1(C_e^T D_e f_e - sqrt(k) rest).
    // Solve S q_f = f_f - C_f solution, then rebuild the FULL connection RHS.
    pub(super) fn complete(&self, external: &[f64], rhs: &mut [f64], solution: &[f64],
        tolerance: f64, roundoff: Option<&[f64]>, gate: &CancelGate) -> Result<Vec<f64>, ModalCouplingError>
    {
        if !self.nullspace.is_empty() {
            let forces: Vec<f64> = self.indices.iter().map(|&i| external[i]).collect();
            for basis in &self.nullspace {
                poll(Some(gate))?;
                if let Some(bounds) = roundoff {
                    let bound = self.indices.iter().zip(basis).try_fold(0.0,
                        |sum, (&i, &a)| finite(sum + a.abs()*bounds[i]))?;
                    if dot(basis, &forces)?.abs() > bound {
                        return Err(invalid("static contact load does not balance its free direction to roundoff"));
                    }
                } else {
                    check_zero_dot(basis, &forces, tolerance,
                        "static forces are incompatible with a declared free direction")?;
                }
            }
        }
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
        // A contact-supported caller selects its pose along this nullspace
        // after the solve. Its accepted force roundoff can otherwise create an
        // immaterial nonzero gauge component in the augmented scratch system.
        if roundoff.is_none() {
            for basis in &self.nullspace {
                poll(Some(gate))?;
                check_zero_dot(basis, &q, tolerance,
                    "static displacement exceeds its declared zero-gauge tolerance")?;
            }
        }
        for (row, &value) in self.rows.iter().zip(&q) {
            poll(Some(gate))?;
            for (r,b) in rhs.iter_mut().zip(row) { *r = finite(*r + b*value)?; }
        }
        Ok(q)
    }

    // This gate is also applied to unit-load compliance queries for contacts;
    // a small Schur residual alone is not a complete free-mode force balance.
    pub(super) fn check_and_insert(&self, external: &[f64], solution: &[f64], free_q: &[f64],
        q: &mut [f64], tolerance: f64, roundoff: Option<&[f64]>, gate: &CancelGate) -> Result<(), ModalCouplingError>
    {
        let mut size = solution.iter().fold(0.0_f64, |m, x| m.max(x.abs()));
        if !self.nullspace.is_empty() {
            // Rest-length equilibration may cancel the whole connection RHS:
            // r0 + C_f^T q_f is rounded zero although both terms are nonzero.
            // Preserve that reconstruction scale instead of normalizing the
            // resulting force roundoff by itself. A = I + C_e^T D_e C_e has
            // ||A^-1||_2 <= 1; sqrt(links)*max(column term sums) bounds its
            // propagated Euclidean magnitude without a second inverse solve.
            let mut terms = vec![0.0; solution.len()];
            for (row, &value) in self.rows.iter().zip(free_q) {
                poll(Some(gate))?;
                for (sum, &entry) in terms.iter_mut().zip(row) {
                    *sum = finite(*sum + finite(entry*value)?.abs())?;
                }
            }
            let maximum = terms.iter().copied().fold(0.0_f64, f64::max);
            size = finite(size + finite(maximum*(solution.len() as f64).sqrt())?)?;
        }
        for (i, &index) in self.indices.iter().enumerate() {
            poll(Some(gate))?;
            let applied = dot(&self.rows[i], solution)?;
            // Each solved reaction carries roundoff of the whole solution's
            // size, so scale by |a_k| times its magnitude rather than |a_k x_k|.
            // A small Schur residual never replaces the final original-state
            // spring, material, energy and per-coordinate force admissions.
            let mut scale = external[index].abs();
            for a in &self.rows[i] { scale = finite(scale + finite(a.abs()*size)?)?; }
            let residual = finite(applied - external[index])?.abs();
            let residual = match roundoff {
                Some(bounds) => (residual-bounds[index]).max(0.0),
                None => residual,
            };
            let relative = if scale == 0.0 { 0.0 } else { residual / scale };
            if relative > tolerance {
                return Err(ModalCouplingError::SolveResidual { relative, tolerance });
            }
            q[index] = free_q[i];
        }
        Ok(())
    }
}

// K = diag(omega^2) + sum k b b^T is PSD, so every admitted null vector
// must vanish on elastic modes and annihilate EACH positive-stiffness b.
// Checking columns separately prevents weak springs from being discarded by
// comparison with a stronger spring. No loads or basis vectors are repaired.
fn admit_nullspace(network: &StaticSystem<'_>, indices: &[usize],
    basis: &[Vec<f64>], gate: &CancelGate) -> Result<Vec<Vec<f64>>, ModalCouplingError>
{
    let n = network.mode_count();
    let tolerance = network.config.solve_relative_tolerance.min(64.0*f64::EPSILON*n as f64);
    for (i, vector) in basis.iter().enumerate() {
        poll(Some(gate))?;
        if vector.len() != n || vector.iter().any(|v| !v.is_finite()) {
            return Err(invalid("static null directions must be finite complete modal vectors"));
        }
        for (mode, &value) in network.models.iter().flat_map(|m| m.modes()).zip(vector) {
            if mode.angular_frequency_rad_s != 0.0 && value != 0.0 {
                return Err(invalid("a static gauge cannot remove elastic modal stiffness"));
            }
        }
        for (j, other) in basis[..=i].iter().enumerate() {
            poll(Some(gate))?;
            let target = if i == j { 1.0 } else { 0.0 };
            if finite(dot(vector, other)? - target)?.abs() > tolerance {
                return Err(invalid("static null directions must be orthonormal"));
            }
        }
        for (column, link) in network.columns.iter().zip(network.connections) {
            poll(Some(gate))?;
            if link.stiffness_n_m > 0.0 {
                check_zero_dot(column, vector, tolerance,
                    "a static gauge does not annihilate a physical spring")?;
            }
        }
    }
    Ok(basis.iter().map(|v| indices.iter().map(|&i| v[i]).collect()).collect())
}

fn check_zero_dot(a: &[f64], b: &[f64], tolerance: f64, what: &'static str)
    -> Result<(), ModalCouplingError>
{
    let mut sum = 0.0;
    let mut scale = 0.0;
    for (&a, &b) in a.iter().zip(b) {
        let term = finite(a*b)?;
        if a != 0.0 && b != 0.0 && term == 0.0 {
            return Err(invalid("static gauge projection is below representable range"));
        }
        sum = finite(sum + term)?;
        scale = finite(scale + term.abs())?;
    }
    if scale != 0.0 && sum.abs()/scale > tolerance {
        return Err(invalid(what));
    }
    Ok(())
}
