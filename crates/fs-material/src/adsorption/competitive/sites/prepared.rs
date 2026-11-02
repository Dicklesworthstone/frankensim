//! Fixed-capacity output for the original shared-vacancy kinetic equations.
use super::*;

/// Reusable bounded output. Evaluation stages its result on the stack and
/// publishes only after the original equations and cancellation checks pass.
/// The Jacobian is diagonal minus rank one; no dense species matrix is stored.
#[derive(Clone, Debug, PartialEq)]
pub struct CompetitiveSiteBuffer {
    count: usize,
    duration: f64,
    occupied: [f64; MAX_COMPETITIVE_SPECIES],
    incoming: [f64; MAX_COMPETITIVE_SPECIES],
    outgoing: [f64; MAX_COMPETITIVE_SPECIES],
    kinetic_residual: [f64; MAX_COMPETITIVE_SPECIES],
    retained: [f64; MAX_COMPETITIVE_SPECIES],
    inverse_desorption: [f64; MAX_COMPETITIVE_SPECIES],
    vacancy: [f64; 2],
    inverse_vacancy_denominator: f64,
    maximum_residual: f64,
    scalar_derivative: Option<f64>,
}
impl Default for CompetitiveSiteBuffer {
    fn default() -> Self {
        Self {
            count: 0,
            duration: 0.0,
            occupied: [0.0; MAX_COMPETITIVE_SPECIES],
            incoming: [0.0; MAX_COMPETITIVE_SPECIES],
            outgoing: [0.0; MAX_COMPETITIVE_SPECIES],
            kinetic_residual: [0.0; MAX_COMPETITIVE_SPECIES],
            retained: [0.0; MAX_COMPETITIVE_SPECIES],
            inverse_desorption: [0.0; MAX_COMPETITIVE_SPECIES],
            vacancy: [0.0; 2],
            inverse_vacancy_denominator: 1.0,
            maximum_residual: 0.0,
            scalar_derivative: None,
        }
    }
}
impl CompetitiveSiteBuffer {
    pub fn final_moles(&self) -> &[f64] {
        &self.occupied[..self.count]
    }
    pub fn gross_adsorption_mol(&self) -> &[f64] {
        &self.incoming[..self.count]
    }
    pub fn gross_desorption_mol(&self) -> &[f64] {
        &self.outgoing[..self.count]
    }
    /// Signed stored change minus the independently evaluated kinetic change.
    pub fn kinetic_residuals_mol(&self) -> &[f64] {
        &self.kinetic_residual[..self.count]
    }
    pub const fn vacancy_moles(&self) -> [f64; 2] {
        self.vacancy
    }
    pub const fn maximum_balance_residual_mol(&self) -> f64 {
        self.maximum_residual
    }
    /// Row factor f_i in d(n_i)/d(a_j)=v*dt/(1+dt*d_j)*(delta_ij-f_i).
    pub fn vacancy_coupling_fraction(&self, species: usize) -> Option<f64> {
        (species < self.count).then(|| self.retained[species] * self.inverse_vacancy_denominator)
    }
    /// Positive column factor v*dt/(1+dt*d_j) [mol s]. Consumers multiply by
    /// their pressure/concentration-to-rate derivative before forming a matrix.
    /// A derivative may be unrepresentable even when the original endpoint fits;
    /// this refusal does not narrow the inventory-only entry point.
    pub fn rate_response_mol_s(&self, species: usize) -> Result<f64> {
        if species >= self.count {
            return Err(invalid("competitive derivative species outside roster"));
        }
        finite((self.duration * self.inverse_desorption[species]) * self.vacancy[1])
    }
    pub fn adsorption_rate_derivative_mol_s(&self, row: usize, column: usize) -> Result<f64> {
        if row >= self.count || column >= self.count {
            return Err(invalid("competitive derivative species outside roster"));
        }
        if let Some(derivative) = self.scalar_derivative {
            return Ok(derivative);
        }
        let fraction = self
            .vacancy_coupling_fraction(row)
            .expect("admitted species");
        let factor = if row != column {
            -fraction
        } else if fraction > 0.5 {
            // Avoid 1-f_i cancellation when one strongly adsorbing species
            // dominates. At most one row takes this O(species) branch.
            let other = sum((0..self.count)
                .filter(|&i| i != row)
                .map(|i| self.retained[i]))?;
            finite((1.0 + other) * self.inverse_vacancy_denominator)?
        } else {
            1.0 - fraction
        };
        finite(((self.duration * self.inverse_desorption[column]) * factor) * self.vacancy[1])
    }
}

/// Same physical endpoint as `backward_euler_sites`, without heap allocation.
/// Inputs and the caller's previous output remain unchanged on any refusal.
pub fn backward_euler_sites_into(
    initial: &[f64],
    capacity: f64,
    adsorption: &[f64],
    desorption: &[f64],
    duration: f64,
    output: &mut CompetitiveSiteBuffer,
    mut cancelled: impl FnMut() -> bool,
) -> Result<()> {
    poll(&mut cancelled)?;
    let n = initial.len();
    if n == 0
        || n > MAX_COMPETITIVE_SPECIES
        || adsorption.len() != n
        || desorption.len() != n
        || !positive(capacity)
        || !nonnegative(duration)
        || initial
            .iter()
            .chain(adsorption)
            .chain(desorption)
            .any(|&x| !nonnegative(x))
    {
        return Err(invalid(
            "bounded matched competitive inventories and finite nonnegative rates required",
        ));
    }
    let vacancy0 = finite(capacity - sum(initial.iter().copied())?)?;
    if vacancy0 < 0.0 {
        return Err(invalid(
            "competitive loadings exceed their ONE shared site capacity",
        ));
    }
    let mut result = CompetitiveSiteBuffer {
        count: n,
        duration,
        vacancy: [vacancy0; 2],
        ..Default::default()
    };
    result.occupied[..n].copy_from_slice(initial);
    if n == 1 {
        let step =
            backward_euler_site(initial[0], capacity, adsorption[0], desorption[0], duration)
                .map_err(|error| invalid(error.message()))?;
        result.inverse_desorption[0] = 1.0 / (1.0 + duration * desorption[0]);
        result.retained[0] = (duration * adsorption[0]) / (1.0 + duration * desorption[0]);
        result.inverse_vacancy_denominator = 1.0 / (1.0 + result.retained[0]);
        result.occupied[0] = step.occupied_moles;
        result.incoming[0] = step.gross_adsorption_mol;
        result.outgoing[0] = step.gross_desorption_mol;
        result.kinetic_residual[0] =
            step.occupied_change_mol - (step.gross_adsorption_mol - step.gross_desorption_mol);
        result.vacancy[1] = step.vacancy_moles;
        result.maximum_residual = step.maximum_balance_residual_mol;
        result.scalar_derivative = Some(step.adsorption_rate_derivative_mol_s);
    } else if duration > 0.0 {
        let mut qa = [0.0; MAX_COMPETITIVE_SPECIES];
        let mut qd = [0.0; MAX_COMPETITIVE_SPECIES];
        let mut released = [0.0; MAX_COMPETITIVE_SPECIES];
        for i in 0..n {
            poll(&mut cancelled)?;
            qa[i] = finite(duration * adsorption[i])?;
            qd[i] = finite(duration * desorption[i])?;
            if (adsorption[i] > 0.0 && qa[i] == 0.0) || (desorption[i] > 0.0 && qd[i] == 0.0) {
                return Err(invalid("nonzero competitive kinetic interval underflow"));
            }
            let denominator = finite(1.0 + qd[i])?;
            result.inverse_desorption[i] = 1.0 / denominator;
            result.retained[i] = qa[i] / denominator;
            released[i] = initial[i] * (qd[i] / denominator);
            if (qa[i] > 0.0 && result.retained[i] == 0.0)
                || (initial[i] > 0.0 && qd[i] > 0.0 && released[i] == 0.0)
            {
                return Err(invalid(
                    "competitive vacancy solve lost a nonzero contribution",
                ));
            }
        }
        // Old vacancy plus released old occupancy avoids near-saturation
        // subtraction of two large, nearly equal total capacities.
        let numerator = sum(std::iter::once(vacancy0).chain(released[..n].iter().copied()))?;
        let denominator = finite(1.0 + sum(result.retained[..n].iter().copied())?)?;
        let vacancy1 = numerator / denominator;
        if !nonnegative(vacancy1) || (numerator > 0.0 && vacancy1 == 0.0) {
            return Err(invalid("competitive vacancy inventory is unrepresentable"));
        }
        result.inverse_vacancy_denominator = 1.0 / denominator;
        for i in 0..n {
            poll(&mut cancelled)?;
            let incoming = finite(qa[i] * vacancy1)?;
            let delta = finite((incoming - qd[i] * initial[i]) / (1.0 + qd[i]))?;
            let next = if delta.abs() < 0.5 * initial[i] {
                initial[i] + delta
            } else {
                (initial[i] + incoming) / (1.0 + qd[i])
            };
            let outgoing = finite(qd[i] * next)?;
            if !nonnegative(next)
                || next > capacity
                || ((initial[i] > 0.0 || incoming > 0.0) && next == 0.0)
                || (qa[i] > 0.0 && vacancy1 > 0.0 && incoming == 0.0)
                || (qd[i] > 0.0 && next > 0.0 && outgoing == 0.0)
            {
                return Err(invalid(
                    "competitive occupied sites or gross streams are unrepresentable",
                ));
            }
            result.occupied[i] = next;
            result.incoming[i] = incoming;
            result.outgoing[i] = outgoing;
            let residual = finite((next - initial[i]) - (incoming - outgoing))?;
            result.kinetic_residual[i] = residual;
            result.maximum_residual = result.maximum_residual.max(residual.abs());
        }
        let residual = finite(sum(result.final_moles().iter().copied())? + vacancy1 - capacity)?;
        result.maximum_residual = result.maximum_residual.max(residual.abs());
        if result.maximum_residual > 128.0 * f64::EPSILON * capacity {
            return Err(invalid(
                "original competitive site equations failed their roundoff balance",
            ));
        }
        result.vacancy[1] = vacancy1;
    }
    poll(&mut cancelled)?;
    *output = result;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn prepared_site_jacobian_matches_independent_perturbations_and_keeps_shared_capacity() {
        let initial = [0.2, 0.3, 0.1];
        let rates = [4.0, 2.0, 0.0];
        let desorption = [1.0, 3.0, 0.2];
        let dt = 0.1;
        let mut work = CompetitiveSiteBuffer::default();
        backward_euler_sites_into(&initial, 1.0, &rates, &desorption, dt, &mut work, || false)
            .unwrap();
        for j in 0..3 {
            let mut changed = rates;
            changed[j] += 1e-4;
            let mut other = CompetitiveSiteBuffer::default();
            backward_euler_sites_into(&initial, 1.0, &changed, &desorption, dt, &mut other, || {
                false
            })
            .unwrap();
            changed[j] += 1e-4;
            let mut twice = CompetitiveSiteBuffer::default();
            backward_euler_sites_into(&initial, 1.0, &changed, &desorption, dt, &mut twice, || {
                false
            })
            .unwrap();
            for i in 0..3 {
                let measured = (4.0 * (other.final_moles()[i] - work.final_moles()[i])
                    - (twice.final_moles()[i] - work.final_moles()[i]))
                    / 2e-4;
                let analytic = work.adsorption_rate_derivative_mol_s(i, j).unwrap();
                assert!(
                    (measured - analytic).abs() < 2e-9,
                    "{i},{j}: {measured} {analytic}"
                );
                if i != j {
                    assert!(analytic <= 0.0);
                }
            }
        }
        let old = backward_euler_sites(&initial, 1.0, &rates, &desorption, dt, || false).unwrap();
        assert_eq!(old.final_moles, work.final_moles());
        assert_eq!(old.gross_adsorption_mol, work.gross_adsorption_mol());
        assert_eq!(old.gross_desorption_mol, work.gross_desorption_mol());
        assert!(
            (work.final_moles().iter().sum::<f64>() + work.vacancy_moles()[1] - 1.0).abs() < 1e-15
        );
    }
    #[test]
    fn derivatives_keep_empty_and_dominant_limits_and_refusals_are_atomic() {
        let mut work = CompetitiveSiteBuffer::default();
        backward_euler_sites_into(
            &[0.0, 0.0],
            1.0,
            &[0.0, 0.0],
            &[0.5, 1.0],
            0.125,
            &mut work,
            || false,
        )
        .unwrap();
        assert_eq!(work.final_moles(), &[0.0, 0.0]);
        assert_eq!(
            work.adsorption_rate_derivative_mol_s(0, 0).unwrap(),
            0.125 / 1.0625
        );
        backward_euler_sites_into(
            &[0.0, 0.0],
            1.0,
            &[1e18, 0.0],
            &[0.0, 0.0],
            1.0,
            &mut work,
            || false,
        )
        .unwrap();
        assert!(work.adsorption_rate_derivative_mol_s(0, 0).unwrap() > 0.0);
        let saved = work.clone();
        assert!(
            backward_euler_sites_into(
                &[0.8, 0.8],
                1.0,
                &[1.0, 1.0],
                &[0.0, 0.0],
                1.0,
                &mut work,
                || false
            )
            .is_err()
        );
        assert_eq!(work, saved);
        assert!(matches!(
            backward_euler_sites_into(
                &[0.0, 0.0],
                1.0,
                &[1.0, 1.0],
                &[0.0, 0.0],
                1.0,
                &mut work,
                || true
            ),
            Err(AdsorptionError::Cancelled)
        ));
        assert_eq!(work, saved);
    }
}
