//! Attached Maugis work on the original native mechanical coordinates.
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
    let attached = contacts
        .iter()
        .any(|s| matches!(s.source, NormalLawSource::MaugisAttached(_)));
    if attached
        && (contacts.iter().any(|s| {
            !matches!(
                s.source,
                NormalLawSource::MaugisAttached(_) | NormalLawSource::DmtSphere(_)
            )
        }) || (inner.discrete_gradient_rule().is_none()
            && inner.exact_quadratic_matrix().is_none()))
    {
        return Err(DContactError::Parameter {
            what: "attached Maugis dynamics needs complete native/quadratic inner work and supported normal paths",
        });
    }
    Ok(attached)
}

pub(super) fn coordinate_scale(storage: &BoundNormalStorage, coordinate: usize) -> Option<f64> {
    let mut scale = storage.inner.coordinate_scale(coordinate);
    if coordinate.is_multiple_of(2) && coordinate / 2 < storage.n_modes {
        for contact in &storage.contacts {
            let weight = contact.collocation[coordinate / 2].abs();
            if weight == 0.0 {
                continue;
            }
            if let NormalLawSource::MaugisAttached(source) = &contact.source {
                let value = source.indentation_scale_m() / weight;
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
            let Ok(from) = contact.approach(a) else {
                out.fill(f64::NAN);
                return;
            };
            let Ok(to) = contact.approach(b) else {
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
}
