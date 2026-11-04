//! Attached adhesive work on the original native mechanical coordinates.
//!
//! A scalar path effort is projected through its physical contact row. This
//! preserves action/reaction even for unequal masses; a global Gonzalez
//! correction could instead put contact work into unrelated momenta.
use super::{BoundNormalContact, BoundNormalStorage, DContactError, NormalLawSource};
use fs_phs::{ConstitutiveDiscreteGradient, Storage};

pub(super) fn admit(
    inner: &dyn Storage,
    contacts: &[BoundNormalContact],
) -> Result<bool, DContactError> {
    let attached = contacts.iter().any(|s| {
        matches!(
            s.source,
            NormalLawSource::MaugisAttached(_) | NormalLawSource::JkrAttached(_)
        )
    });
    if attached
        && inner.discrete_gradient_rule().is_none()
        && inner.exact_quadratic_matrix().is_none()
    {
        return Err(DContactError::Parameter {
            what: "attached adhesion dynamics needs complete native/quadratic inner work",
        });
    }
    Ok(attached)
}

pub(super) fn coordinate_scale(storage: &BoundNormalStorage, coordinate: usize) -> Option<f64> {
    let mut scale = storage.inner.coordinate_scale(coordinate);
    if !storage.component_rule {
        return scale;
    }
    if coordinate.is_multiple_of(2) && coordinate / 2 < storage.n_modes {
        for contact in &storage.contacts {
            let weight = contact.collocation[coordinate / 2].abs();
            if weight == 0.0 {
                continue;
            }
            let contact_scale = match &contact.source {
                NormalLawSource::MaugisAttached(source) => Some(source.indentation_scale_m()),
                NormalLawSource::JkrAttached(source) => Some(source.indentation_scale_m()),
                NormalLawSource::DmtSphere(source) => Some(source.interaction_range_m()),
                NormalLawSource::HertzPoint(source) => {
                    Some(source.envelope().maximum_indentation_m)
                }
                NormalLawSource::FiniteGap(source) => source.nodes.last().map(|n| n.approach_m),
            };
            if let Some(contact_scale) = contact_scale {
                let value = contact_scale / weight;
                if value.is_finite() && value > 0.0 {
                    scale = Some(scale.map_or(value, |old| old.min(value)));
                }
            }
        }
    }
    scale
}

impl ConstitutiveDiscreteGradient for BoundNormalStorage {
    fn evaluate(&self, a: &[f64], b: &[f64], midpoint: &[f64], out: &mut [f64]) {
        if !self.component_rule
            || a.len() != b.len()
            || a.len() != midpoint.len()
            || a.len() != out.len()
            || self.try_hamiltonian(a).is_err()
            || self.try_hamiltonian(b).is_err()
            || self.try_hamiltonian(midpoint).is_err()
        {
            out.fill(f64::NAN);
            return;
        }
        if let Some(rule) = self.inner.discrete_gradient_rule() {
            rule.evaluate(a, b, midpoint, out);
        } else {
            let Some(q) = self.inner.exact_quadratic_matrix() else {
                out.fill(f64::NAN);
                return;
            };
            let n = a.len();
            if q.len() != n * n {
                out.fill(f64::NAN);
                return;
            }
            // Exact symmetric gradient of the original quadratic Hamiltonian,
            // including its admitted tiny asymmetry and midpoint rounding.
            for i in 0..n {
                out[i] = 0.0;
                for j in 0..n {
                    let mean = 0.5 * a[j] + 0.5 * b[j];
                    out[i] += (0.5 * q[i * n + j] + 0.5 * q[j * n + i]) * mean;
                }
            }
        }
        for contact in &self.contacts {
            let Ok(from) = contact.raw_approach(a) else {
                out.fill(f64::NAN);
                return;
            };
            let Ok(to) = contact.raw_approach(b) else {
                out.fill(f64::NAN);
                return;
            };
            let Ok(Some(force)) = contact.source.mean_force_between(from, to) else {
                out.fill(f64::NAN);
                return;
            };
            for (index, weight) in contact.collocation.iter().enumerate() {
                out[2 * index] += weight * force;
            }
        }
        if out.iter().any(|v| !v.is_finite()) {
            out.fill(f64::NAN);
        }
    }

    fn end_jacobian_action(
        &self,
        a: &[f64],
        b: &[f64],
        direction: &[f64],
        out: &mut [f64],
    ) -> bool {
        // Maugis retains its numerical route until its original path owner
        // supplies a derivative. Other source owners compose their own exact
        // path derivatives, including unilateral and finite-range crossings.
        if !self.component_rule
            || self
                .contacts
                .iter()
                .any(|contact| matches!(contact.source, NormalLawSource::MaugisAttached(_)))
        {
            return false;
        }
        if a.len() != b.len()
            || a.len() != direction.len()
            || a.len() != out.len()
            || direction.iter().any(|v| !v.is_finite())
            || self.try_hamiltonian(a).is_err()
            || self.try_hamiltonian(b).is_err()
        {
            out.fill(f64::NAN);
            return true;
        }
        if let Some(rule) = self.inner.discrete_gradient_rule() {
            // evaluate() selects this same rule before any quadratic view.
            // Its absence of a derivative cannot be replaced with Q/2.
            if !rule.end_jacobian_action(a, b, direction, out) {
                return false;
            }
        } else {
            let Some(q) = self.inner.exact_quadratic_matrix() else {
                return false;
            };
            let n = a.len();
            if n.checked_mul(n) != Some(q.len()) {
                out.fill(f64::NAN);
                return true;
            }
            for i in 0..n {
                out[i] = 0.0;
                for j in 0..n {
                    out[i] += (0.25 * q[i * n + j] + 0.25 * q[j * n + i]) * direction[j];
                }
            }
        }
        for contact in &self.contacts {
            let (Ok(from), Ok(to)) = (contact.raw_approach(a), contact.raw_approach(b)) else {
                out.fill(f64::NAN);
                return true;
            };
            let derivative = match contact.source.mean_force_end_derivative(from, to) {
                Ok(Some(value)) => value,
                Ok(None) => return false,
                Err(_) => {
                    // In particular, coincident JKR fold endpoints have no
                    // finite tangent. Do not probe through that singularity.
                    out.fill(f64::NAN);
                    return true;
                }
            };
            let approach_direction: f64 = contact
                .collocation
                .iter()
                .enumerate()
                .map(|(index, weight)| weight * direction[2 * index])
                .sum();
            let force_direction = derivative * approach_direction;
            for (index, weight) in contact.collocation.iter().enumerate() {
                out[2 * index] += weight * force_direction;
            }
        }
        if out.iter().any(|v| !v.is_finite()) {
            out.fill(f64::NAN);
        }
        true
    }
}
