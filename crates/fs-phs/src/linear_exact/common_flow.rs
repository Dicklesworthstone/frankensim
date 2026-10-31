//! The original quadratic common-flow DAE, with any full-rank retained pairing.
//!
//! x'=(J-R)Qx+B lambda+Gu, B^T Qx=0. Only S=B^T QB must be
//! positive definite: Q may have physical zero-storage coordinates. The
//! differentiated constraint gives lambda=-S^-1 B^T Q((J-R)Qx+Gu).
//! P=I-B S^-1 B^T Q then gives Jc=PJP^T, Rc=PRP^T, Gc=PG.
//! On an admitted state P^T Qx=Qx; thus this is the original constrained
//! equation, not a different passive surrogate. No state is projected.
use super::{MAX_LINEAR_EXPONENTIAL_BLOCK_DIM, finite, product, refuse};
use crate::{PhsError, PortHamiltonian, PreparedStepError};
use fs_la::factor::cholesky;
use fs_math::det;

/// Bounded algebra shared by quadratic fluid, mechanical and circuit junctions.
/// This constructs no integrator and assigns no physical meaning to R or a
/// quadratic observation. Consumers retain the original heat and port work.
pub struct CommonFlowProjection {
    n: usize,
    m: usize,
    count: usize,
    q: Vec<f64>,
    constraint: Vec<f64>,
    root_inverse_inertia: Vec<f64>,
    projector: Vec<f64>,
    j: Vec<f64>,
    r: Vec<f64>,
    g: Vec<f64>,
    reaction_state: Vec<f64>,
    reaction_input: Vec<f64>,
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), PreparedStepError> {
    if cancelled() {
        Err(PreparedStepError::Cancelled)
    } else {
        Ok(())
    }
}

impl CommonFlowProjection {
    /// B is row-major n by count in the ORIGINAL state coordinates. Its
    /// columns are work-conjugate multiplier inputs, not arbitrary flow rows.
    /// The complete state and every independent constraint are retained.
    pub fn new(
        original: &PortHamiltonian,
        b: &[f64],
        count: usize,
        max_states: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, PreparedStepError> {
        poll(cancelled)?;
        let (n, m) = (original.n, original.m);
        if n == 0
            || n > max_states
            || n > MAX_LINEAR_EXPONENTIAL_BLOCK_DIM
            || m > MAX_LINEAR_EXPONENTIAL_BLOCK_DIM
            || count == 0
            || count > n
            || b.len() != n * count
            || !finite(b)
            || original.j.len() != n * n
            || original.r.len() != n * n
            || original.g.len() != n * m
            || !finite(&original.j)
            || !finite(&original.r)
            || !finite(&original.g)
        {
            return Err(refuse("bounded complete common-flow dimensions required").into());
        }
        let part = original
            .storage
            .quadratic_part()
            .ok_or_else(|| refuse("common-flow requires original quadratic storage"))?;
        if !part.exact || part.matrix.len() != n * n || !finite(&part.matrix) {
            return Err(refuse("common-flow requires exact finite quadratic storage").into());
        }
        let q = part.matrix;
        crate::require_skew(&original.j, n, "common-flow original J")?;
        crate::require_symmetric(&original.r, n, "common-flow original R")?;
        crate::require_psd(&original.r, n, "common-flow original R")?;
        crate::require_symmetric(&q, n, "common-flow original Q")?;
        crate::require_psd(&q, n, "common-flow original Q")?;
        let mut constraint = vec![0.0; count * n];
        for a in 0..count {
            poll(cancelled)?;
            for j in 0..n {
                constraint[a * n + j] = (0..n).map(|i| b[i * count + a] * q[i * n + j]).sum();
            }
        }
        let mut inertia = vec![0.0; count * count];
        for a in 0..count {
            for c in 0..count {
                inertia[a * count + c] = (0..n)
                    .map(|i| constraint[a * n + i] * b[i * count + c])
                    .sum();
            }
        }
        let roots: Vec<_> = (0..count)
            .map(|a| det::sqrt(inertia[a * count + a]))
            .collect();
        if !finite(&constraint)
            || !finite(&inertia)
            || roots.iter().any(|s| !s.is_finite() || *s <= 0.0)
        {
            return Err(
                refuse("common-flow inverse inertia must be positive and representable").into(),
            );
        }
        let scaled: Vec<_> = (0..count * count)
            .map(|i| inertia[i] / roots[i / count] / roots[i % count])
            .collect();
        if !finite(&scaled) {
            return Err(refuse("common-flow normalized inverse inertia is unrepresentable").into());
        }
        let factor = cholesky(&scaled, count)
            .map_err(|_| refuse("common-flow constraints must be independent"))?;
        if (0..count).any(|i| factor.l(i, i).powi(2) <= 512.0 * f64::EPSILON * count as f64) {
            return Err(refuse("common-flow constraint rank is unresolved at working precision").into());
        }
        let mut normalized = vec![0.0; count * n];
        for j in 0..n {
            poll(cancelled)?;
            let mut rhs: Vec<_> = (0..count)
                .map(|a| constraint[a * n + j] / roots[a])
                .collect();
            factor.solve(&mut rhs);
            for a in 0..count {
                normalized[a * n + j] = rhs[a] / roots[a];
            }
        }
        // Check against ALL original entries; Cholesky's lower-triangle read
        // cannot conceal an inconsistent upper triangle or unresolved rank.
        for a in 0..count {
            poll(cancelled)?;
            for j in 0..n {
                let mut value = 0.0;
                let mut bound = constraint[a * n + j].abs();
                for c in 0..count {
                    let term = inertia[a * count + c] * normalized[c * n + j];
                    value += term;
                    bound += term.abs();
                }
                if !value.is_finite()
                    || !bound.is_finite()
                    || (value - constraint[a * n + j]).abs()
                        > 512.0 * f64::EPSILON * count as f64 * bound
                {
                    return Err(
                        refuse("original common-flow inertia solve failed its residual").into(),
                    );
                }
            }
        }
        let mut reaction_effort = vec![0.0; count * n];
        let mut reaction_state = vec![0.0; count * n];
        let mut reaction_input = vec![0.0; count * m];
        for a in 0..count {
            poll(cancelled)?;
            for j in 0..n {
                reaction_effort[a * n + j] = -(0..n)
                    .map(|i| {
                        normalized[a * n + i] * (original.j[i * n + j] - original.r[i * n + j])
                    })
                    .sum::<f64>();
            }
            for j in 0..n {
                reaction_state[a * n + j] = (0..n)
                    .map(|i| reaction_effort[a * n + i] * q[i * n + j])
                    .sum();
            }
            for p in 0..m {
                reaction_input[a * m + p] = -(0..n)
                    .map(|i| normalized[a * n + i] * original.g[i * m + p])
                    .sum::<f64>();
            }
        }
        let mut projector = vec![0.0; n * n];
        for i in 0..n {
            poll(cancelled)?;
            for j in 0..n {
                projector[i * n + j] = if i == j { 1.0 } else { 0.0 }
                    - (0..count)
                        .map(|a| b[i * count + a] * normalized[a * n + j])
                        .sum::<f64>();
            }
        }
        let transpose: Vec<_> = (0..n * n).map(|i| projector[(i % n) * n + i / n]).collect();
        let mut j = product(&product(&projector, &original.j, n), &transpose, n);
        poll(cancelled)?;
        let mut r = product(&product(&projector, &original.r, n), &transpose, n);
        let mut g = vec![0.0; n * m];
        for i in 0..n {
            poll(cancelled)?;
            // One contraction stamps each analytically skew/symmetric pair;
            // this does not alter any original constitutive coefficient.
            j[i * n + i] = 0.0;
            for k in i + 1..n {
                j[k * n + i] = -j[i * n + k];
                r[k * n + i] = r[i * n + k];
            }
            for p in 0..m {
                g[i * m + p] = (0..n)
                    .map(|k| projector[i * n + k] * original.g[k * m + p])
                    .sum();
            }
        }
        if [&projector, &j, &r, &g, &reaction_state, &reaction_input]
            .iter()
            .any(|v| !finite(v))
        {
            return Err(refuse("common-flow projection or reaction is unrepresentable").into());
        }
        poll(cancelled)?;
        Ok(Self {
            n,
            m,
            count,
            q,
            constraint,
            root_inverse_inertia: roots,
            projector,
            j,
            r,
            g,
            reaction_state,
            reaction_input,
        })
    }

    /// Row-major P, preserving all original coordinates.
    pub fn projector(&self) -> &[f64] {
        &self.projector
    }
    /// Original B^T Q, count by state dimension.
    pub fn constraint_rows(&self) -> &[f64] {
        &self.constraint
    }
    /// Jc, Rc, Gc in original state/port order.
    pub fn structure(&self) -> (&[f64], &[f64], &[f64]) {
        (&self.j, &self.r, &self.g)
    }
    /// lambda = reaction_state*x + reaction_input*u, both row-major.
    pub fn reaction_maps(&self) -> (&[f64], &[f64]) {
        (&self.reaction_state, &self.reaction_input)
    }

    /// Refuse incompatible initial/endpoint flow without changing its energy.
    /// The scale also covers potential-to-kinetic roundoff at zero crossings,
    /// using |B^T Q dx| <= sqrt(B^T QB)*||Q^(1/2) dx||.
    pub fn admit_state(&self, state: &[f64]) -> Result<(), PhsError> {
        if state.len() != self.n || !finite(state) {
            return Err(refuse(
                "common-flow requires a complete finite physical state",
            ));
        }
        let energy_scale: f64 = (0..self.n)
            .map(|i| {
                (state[i]
                    * (0..self.n)
                        .map(|j| self.q[i * self.n + j] * state[j])
                        .sum::<f64>())
                .abs()
            })
            .sum();
        if !energy_scale.is_finite() {
            return Err(refuse("common-flow energy scale is unrepresentable"));
        }
        for a in 0..self.count {
            let mut residual = 0.0;
            let mut scale = 0.0;
            for (c, x) in self.constraint[a * self.n..(a + 1) * self.n]
                .iter()
                .zip(state)
            {
                residual += c * x;
                scale += (c * x).abs();
            }
            scale = scale.max(self.root_inverse_inertia[a] * det::sqrt(energy_scale));
            if !residual.is_finite()
                || !scale.is_finite()
                || residual.abs() > 64.0 * f64::EPSILON * self.n as f64 * scale
            {
                return Err(refuse(
                    "common-flow original velocity constraint is incompatible",
                ));
            }
        }
        Ok(())
    }

    /// Actual differentiated-constraint multiplier for this held external input.
    pub fn reaction(&self, state: &[f64], input: &[f64]) -> Result<Vec<f64>, PhsError> {
        let mut values = vec![0.0; self.count];
        self.reaction_into(state, input, &mut values)?;
        Ok(values)
    }

    /// Same multiplier map without allocation; caller output is unchanged on refusal.
    pub fn reaction_into(
        &self,
        state: &[f64],
        input: &[f64],
        output: &mut [f64],
    ) -> Result<(), PhsError> {
        self.admit_state(state)?;
        if input.len() != self.m || output.len() != self.count || !finite(input) {
            return Err(refuse("common-flow requires complete finite held inputs"));
        }
        let value = |a: usize| {
            self.reaction_state[a * self.n..(a + 1) * self.n]
                .iter()
                .zip(state)
                .map(|(a, x)| a * x)
                .sum::<f64>()
                + self.reaction_input[a * self.m..(a + 1) * self.m]
                    .iter()
                    .zip(input)
                    .map(|(a, u)| a * u)
                    .sum::<f64>()
        };
        if (0..self.count).any(|a| !value(a).is_finite()) {
            return Err(refuse("common-flow reaction is nonfinite"));
        }
        for (a, out) in output.iter_mut().enumerate() {
            *out = value(a);
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::QuadraticStorage;

    fn masses() -> PortHamiltonian {
        // Three masses and an unrestrained displacement, retaining singular H.
        let mut q = vec![0.0; 16];
        q[0] = 0.5;
        q[5] = 1.0 / 3.0;
        q[10] = 0.2;
        let mut j = vec![0.0; 16];
        j[3] = -1.0;
        j[12] = 1.0;
        let mut r = vec![0.0; 16];
        r[0] = 0.4;
        r[5] = 0.6;
        r[10] = 1.0;
        PortHamiltonian::new(
            4,
            2,
            j,
            r,
            vec![1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
            Box::new(QuadraticStorage::new(q, 4).unwrap()),
        )
        .unwrap()
    }
    fn b() -> Vec<f64> {
        vec![1.0, 0.0, -1.0, 1.0, 0.0, -1.0, 0.0, 0.0]
    }

    #[test]
    fn g1_two_constraints_recover_actual_common_acceleration_reactions_and_heat() {
        let source = masses();
        let projection = CommonFlowProjection::new(&source, &b(), 2, 4, &mut || false).unwrap();
        let state = [1.4, 2.1, 3.5, 99.0];
        let input = [3.0, -1.0];
        let reaction = projection.reaction(&state, &input).unwrap();
        assert!((reaction[0] + 2.60).abs() < 2e-14);
        assert!((reaction[1] + 2.00).abs() < 2e-14);
        let effort = source.effort(&state);
        let (j, r, g) = projection.structure();
        let derivative: Vec<_> = (0..4)
            .map(|i| {
                (0..4)
                    .map(|k| (j[i * 4 + k] - r[i * 4 + k]) * effort[k])
                    .sum::<f64>()
                    + (0..2).map(|p| g[i * 2 + p] * input[p]).sum::<f64>()
            })
            .collect();
        for (actual, expected) in derivative.iter().zip([0.12, 0.18, 0.30, 0.7]) {
            assert!((actual - expected).abs() < 2e-14, "{derivative:?}");
        }
        let (oj, orr, og) = source.structure();
        for i in 0..4 {
            let original = (0..4)
                .map(|k| (oj[i * 4 + k] - orr[i * 4 + k]) * effort[k])
                .sum::<f64>()
                + (0..2)
                    .map(|p| og[i * 2 + p] * input[p] + b()[i * 2 + p] * reaction[p])
                    .sum::<f64>();
            assert!((derivative[i] - original).abs() < 2e-14);
        }
        let rate: f64 = effort.iter().zip(&derivative).map(|(e, d)| e * d).sum();
        assert!((rate + 0.98 - 1.4).abs() < 2e-14);
        assert_eq!(state[3], 99.0);
    }

    #[test]
    fn g3_nonorthogonal_constraint_coordinates_do_not_change_physical_dynamics() {
        let source = masses();
        let original = b();
        let transform = [2.0, 0.5, -0.25, 3.0];
        let mapped: Vec<_> = (0..8)
            .map(|i| {
                (0..2)
                    .map(|k| original[(i / 2) * 2 + k] * transform[k * 2 + i % 2])
                    .sum()
            })
            .collect();
        let a = CommonFlowProjection::new(&source, &original, 2, 4, &mut || false).unwrap();
        let b = CommonFlowProjection::new(&source, &mapped, 2, 4, &mut || false).unwrap();
        for (left, right) in a.projector().iter().zip(b.projector()) {
            assert!((left - right).abs() < 2e-14);
        }
        let x = [1.4, 2.1, 3.5, -70.0];
        let la = a.reaction(&x, &[3.0, -1.0]).unwrap();
        let lb = b.reaction(&x, &[3.0, -1.0]).unwrap();
        for i in 0..2 {
            assert!(
                (la[i] - (0..2).map(|j| transform[i * 2 + j] * lb[j]).sum::<f64>()).abs() < 2e-14
            );
        }
    }

    #[test]
    fn g0_incompatible_rank_overflow_and_malformed_sources_refuse_atomically() {
        let mut source = masses();
        let projection = CommonFlowProjection::new(&source, &b(), 2, 4, &mut || false).unwrap();
        let mut out = [55.0, 66.0];
        assert!(
            projection
                .reaction_into(&[1.4, 2.2, 3.5, 0.0], &[0.0, 0.0], &mut out)
                .is_err()
        );
        assert_eq!(out, [55.0, 66.0]);
        assert!(
            projection
                .admit_state(&[f64::MAX, f64::MAX, f64::MAX, 0.0])
                .is_err()
        );
        let rank_one = [1.0, 1.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0];
        assert!(CommonFlowProjection::new(&source, &rank_one, 2, 4, &mut || false).is_err());
        assert!(CommonFlowProjection::new(&source, &b(), 2, 3, &mut || false).is_err());
        assert!(matches!(
            CommonFlowProjection::new(&source, &b(), 2, 4, &mut || true),
            Err(PreparedStepError::Cancelled)
        ));
        source.j.pop();
        assert!(CommonFlowProjection::new(&source, &b(), 2, 4, &mut || false).is_err());
        let mut source = masses();
        source.g[0] = f64::INFINITY;
        assert!(CommonFlowProjection::new(&source, &b(), 2, 4, &mut || false).is_err());
    }
}
