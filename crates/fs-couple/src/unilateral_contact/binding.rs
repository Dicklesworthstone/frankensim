//! Source-defined normal potentials on caller-owned mechanical coordinates.
//! Each site is one physical point or one complete finite-gap footprint; the
//! latter already includes its BEM area measures and is never counted twice.
use fs_contact::normal_patch::source::SourcePointContact;
use fs_contact::normal_patch::{NormalLawSource, NormalLawSourceError};
use fs_dcontact::{DContactError, Obstacle};
use fs_phs::Storage;

/// Lower an admitted point owner to the original distributed potential.
/// Point quadrature is exactly one. K/alpha/chi cannot be overridden by a caller.
pub fn bound_point_obstacle(
    source: &SourcePointContact,
    collocation: Vec<f64>,
    gap_m: f64,
    key: &str,
) -> Result<Obstacle, DContactError> {
    let modes = collocation.len();
    if modes == 0 || key.trim().is_empty() {
        return Err(DContactError::Shape {
            what: "bound contact needs an identity and mechanical map",
        });
    }
    Obstacle::new(
        collocation,
        1,
        modes,
        vec![gap_m],
        vec![1.0],
        source.stiffness_n_per_m_3_2(),
        1.5,
        format!("source-point-contact/{key}/{}", source.model().identity()),
    )?
    .with_internal_loss(source.dissipation_s_per_m())
}

/// One immutable work-conjugate map onto a sourced normal potential.
#[derive(Debug, Clone, PartialEq)]
pub struct BoundNormalContact {
    source: NormalLawSource,
    collocation: Vec<f64>,
    gap_m: f64,
}
impl BoundNormalContact {
    /// Positive closure is `collocation dot q - gap_m` [m]. This constructor
    /// requires the complete source object rather than authored K/alpha/chi.
    pub fn new(
        source: NormalLawSource,
        collocation: Vec<f64>,
        gap_m: f64,
    ) -> Result<Self, DContactError> {
        if collocation.is_empty()
            || collocation.iter().any(|v| !v.is_finite())
            || collocation.iter().all(|v| *v == 0.0)
            || !gap_m.is_finite()
        {
            return Err(DContactError::Parameter {
                what: "bound normal contact needs a finite nonzero mechanical map and gap",
            });
        }
        // Unilateral owners need an unloaded origin. Adhesive owners carry
        // their own detached reference and must retain nonzero force at touch.
        let origin = source
            .energy_and_force(0.0)
            .map_err(|_| DContactError::Parameter {
                what: "bound source must cover zero approach",
            })?;
        if !source.supports_signed_approach() && (origin.0 != 0.0 || origin.1 != 0.0) {
            return Err(DContactError::Parameter {
                what: "unilateral source zero approach must have zero energy and force",
            });
        }
        Ok(Self {
            source,
            collocation,
            gap_m,
        })
    }
    /// Original immutable source and its authority/operating envelope.
    pub const fn source(&self) -> &NormalLawSource {
        &self.source
    }

    fn approach(&self, x: &[f64]) -> Result<f64, DContactError> {
        if self
            .collocation
            .len()
            .checked_mul(2)
            .is_none_or(|n| x.len() < n)
            || x.iter().any(|v| !v.is_finite())
        {
            return Err(DContactError::Shape {
                what: "bound contact finite mechanical state prefix",
            });
        }
        let mut delta = -self.gap_m;
        for (k, b) in self.collocation.iter().enumerate() {
            delta += b * x[2 * k];
        }
        if !delta.is_finite() {
            return Err(DContactError::Parameter {
                what: "bound contact approach overflow",
            });
        }
        Ok(if self.source.supports_signed_approach() {
            delta
        } else {
            delta.max(0.0)
        })
    }
}

/// Typed binding or original law refusal; no missing field becomes zero force.
#[derive(Debug, Clone, PartialEq)]
pub enum BoundNormalContactError {
    /// Mechanical shape/map admission failed.
    Contact(DContactError),
    /// Original source owner refused this configuration.
    Source(NormalLawSourceError),
}
impl core::fmt::Display for BoundNormalContactError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Contact(e) => core::fmt::Display::fmt(e, f),
            Self::Source(e) => core::fmt::Display::fmt(e, f),
        }
    }
}
impl std::error::Error for BoundNormalContactError {}

/// Existing pHS storage contract extended with source-defined normal potentials.
/// The caller owns stepping and cancellation. No new contact or stepping trait.
pub struct BoundNormalStorage {
    inner: Box<dyn Storage>,
    n_modes: usize,
    contacts: Vec<BoundNormalContact>,
}
impl BoundNormalStorage {
    /// Compose complete mechanical maps. Every footprint contributes once;
    /// line/area integration belongs to the finite-gap source's geometry solve.
    pub fn new(
        inner: Box<dyn Storage>,
        n_modes: usize,
        contacts: Vec<BoundNormalContact>,
    ) -> Result<Self, DContactError> {
        if n_modes == 0 || contacts.iter().any(|s| s.collocation.len() != n_modes) {
            return Err(DContactError::Shape {
                what: "bound contacts must cover the complete mechanical basis",
            });
        }
        Ok(Self {
            inner,
            n_modes,
            contacts,
        })
    }
    /// Evaluate stored energy while retaining the original typed source refusal.
    pub fn try_hamiltonian(&self, x: &[f64]) -> Result<f64, BoundNormalContactError> {
        if self.n_modes.checked_mul(2).is_none_or(|n| x.len() < n)
            || x.iter().any(|v| !v.is_finite())
        {
            return Err(BoundNormalContactError::Contact(DContactError::Shape {
                what: "bound storage finite state prefix",
            }));
        }
        let mut energy = self.inner.hamiltonian(x);
        for site in &self.contacts {
            let delta = site.approach(x).map_err(BoundNormalContactError::Contact)?;
            energy += site
                .source
                .energy_and_force(delta)
                .map_err(BoundNormalContactError::Source)?
                .0;
        }
        if !energy.is_finite() {
            return Err(BoundNormalContactError::Contact(DContactError::Parameter {
                what: "bound contact energy overflow",
            }));
        }
        Ok(energy)
    }
    /// Fill the exact conservative gradient. On failure the trial buffer must
    /// be discarded. State and source objects remain immutable.
    pub fn try_gradient(&self, x: &[f64], out: &mut [f64]) -> Result<(), BoundNormalContactError> {
        if out.len() != x.len() {
            return Err(BoundNormalContactError::Contact(DContactError::Shape {
                what: "bound gradient output",
            }));
        }
        self.try_hamiltonian(x)?;
        self.inner.gradient(x, out);
        for site in &self.contacts {
            let delta = site.approach(x).map_err(BoundNormalContactError::Contact)?;
            let force = site
                .source
                .energy_and_force(delta)
                .map_err(BoundNormalContactError::Source)?
                .1;
            for (k, b) in site.collocation.iter().enumerate() {
                out[2 * k] += b * force;
            }
        }
        if out.iter().any(|v| !v.is_finite()) {
            return Err(BoundNormalContactError::Contact(DContactError::Parameter {
                what: "bound gradient overflow",
            }));
        }
        Ok(())
    }

    /// Positive resisting state-flow for the same implicit pHS equation.
    /// `effort` is its actual discrete gradient; `step_s` is explicit. The
    /// point owner enforces the complete rate/step envelope before loss is
    /// projected. Finite-gap elastic and conservative DMT owners contribute no
    /// internal loss; adhesive attraction is not dissipative resistance.
    /// On refusal discard the trial output; no physical state is changed.
    pub fn dissipative_flow_into(
        &self,
        x: &[f64],
        effort: &[f64],
        step_s: f64,
        out: &mut [f64],
    ) -> Result<(), BoundNormalContactError> {
        if effort.len() != x.len()
            || out.len() != x.len()
            || effort.iter().any(|e| !e.is_finite())
            || !step_s.is_finite()
            || step_s <= 0.0
        {
            return Err(BoundNormalContactError::Contact(DContactError::Shape {
                what: "bound contact effort, output and positive step",
            }));
        }
        self.try_hamiltonian(x)?;
        out.fill(0.0);
        for site in &self.contacts {
            let delta = site.approach(x).map_err(BoundNormalContactError::Contact)?;
            if let NormalLawSource::HertzPoint(source) = &site.source {
                let rate = site
                    .collocation
                    .iter()
                    .enumerate()
                    .map(|(k, b)| b * effort[2 * k + 1])
                    .sum::<f64>();
                source
                    .validate_kinematics(delta, rate, step_s)
                    .map_err(|e| BoundNormalContactError::Source(NormalLawSourceError::Point(e)))?;
                let elastic = site
                    .source
                    .energy_and_force(delta)
                    .map_err(BoundNormalContactError::Source)?
                    .1;
                let force = elastic * (source.dissipation_s_per_m() * rate);
                for (k, b) in site.collocation.iter().enumerate() {
                    out[2 * k + 1] += b * force;
                }
            }
        }
        if out.iter().any(|v| !v.is_finite()) {
            return Err(BoundNormalContactError::Contact(DContactError::Parameter {
                what: "bound normal loss overflow",
            }));
        }
        Ok(())
    }
}
impl Storage for BoundNormalStorage {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        // Storage is infallible; the pHS owner rejects nonfinite trial values.
        // Direct callers use try_hamiltonian to retain the original source error.
        self.try_hamiltonian(x).unwrap_or(f64::NAN)
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        if self.try_gradient(x, out).is_err() {
            out.fill(f64::NAN);
        }
    }

    fn midpoint_energy_remainder(&self, a: &[f64], b: &[f64], mid: &[f64]) -> Option<f64> {
        if a.len() != b.len()
            || a.len() != mid.len()
            || a.iter().chain(b).chain(mid).any(|v| !v.is_finite())
        {
            return None;
        }
        // Never return only a subset of the stored-energy contributions.
        let mut sum = if let Some(value) = self.inner.midpoint_energy_remainder(a, b, mid) {
            value
        } else {
            let q = self.inner.exact_quadratic_matrix()?;
            let n = a.len();
            if q.len() != n.checked_mul(n)? {
                return None;
            }
            let mut total = 0.0;
            let mut error = 0.0;
            for i in 0..n {
                let mut offset = 0.0;
                let mut offset_error = 0.0;
                for j in 0..n {
                    let (center, low) = fs_math::eft::two_sum(0.5 * a[j], 0.5 * b[j]);
                    compensated_add(
                        &mut offset,
                        &mut offset_error,
                        q[i * n + j] * ((center - mid[j]) + low),
                    );
                }
                compensated_add(
                    &mut total,
                    &mut error,
                    (b[i] - a[i]) * (offset + offset_error),
                );
            }
            total + error
        };
        let mut error = 0.0;
        for site in &self.contacts {
            let da = site.approach(a).ok()?;
            let db = site.approach(b).ok()?;
            let dm = site.approach(mid).ok()?;
            let remainder = site.source.midpoint_energy_remainder(da, db, dm).ok()??;
            let force = site.source.energy_and_force(dm).ok()?.1;
            let mut mapped_dx = 0.0;
            let mut mapped_error = 0.0;
            for (k, weight) in site.collocation.iter().enumerate() {
                compensated_add(
                    &mut mapped_dx,
                    &mut mapped_error,
                    weight * (b[2 * k] - a[2 * k]),
                );
            }
            let (delta, delta_error) = fs_math::eft::two_sum(db, -da);
            let correction = force * (((delta - mapped_dx) + delta_error) - mapped_error);
            compensated_add(&mut sum, &mut error, remainder);
            compensated_add(&mut sum, &mut error, correction);
        }
        let total = sum + error;
        total.is_finite().then_some(total)
    }
}

fn compensated_add(sum: &mut f64, error: &mut f64, value: f64) {
    let (high, low) = fs_math::eft::two_sum(*sum, value);
    *sum = high;
    *error += low;
}

#[cfg(test)]
mod adhesive_tests;

#[cfg(test)]
#[path = "../../../fs-contact/tests/source_point/common.rs"]
mod source_fixture;

#[cfg(test)]
mod tests {
    use super::*;
    struct UnitMass;
    impl Storage for UnitMass {
        fn hamiltonian(&self, x: &[f64]) -> f64 {
            0.5 * x[1] * x[1]
        }
        fn gradient(&self, x: &[f64], out: &mut [f64]) {
            out.fill(0.0);
            out[1] = x[1];
        }
    }

    #[test]
    fn g1_bound_storage_and_point_lowering_share_the_original_potential_and_loss() {
        let source = source_fixture::source(1e7, Some(0.4));
        let obstacle = bound_point_obstacle(&source, vec![1.0], 0.0, "fixture").unwrap();
        let reference =
            fs_dcontact::ContactStorage::new(Box::new(UnitMass), 1, vec![obstacle]).unwrap();
        let site =
            BoundNormalContact::new(NormalLawSource::from(source.clone()), vec![1.0], 0.0).unwrap();
        let storage = BoundNormalStorage::new(Box::new(UnitMass), 1, vec![site]).unwrap();
        for delta in [0.0, 1e-6, 1e-5, 1e-4] {
            let x = [delta, 0.0];
            let mut gradient = [0.0; 2];
            let mut expected = [0.0; 2];
            storage.try_gradient(&x, &mut gradient).unwrap();
            reference.gradient(&x, &mut expected);
            assert!((gradient[0] - expected[0]).abs() <= 1e-12 * expected[0].abs().max(1e-30));
            let energy = reference.hamiltonian(&x);
            assert!(
                (storage.try_hamiltonian(&x).unwrap() - energy).abs()
                    <= 1e-12 * energy.abs().max(1e-30)
            );
            for rate in [-0.2, 0.0, 0.2] {
                let effort = [0.0, rate];
                let mut loss = [0.0; 2];
                let mut other = [0.0; 2];
                storage
                    .dissipative_flow_into(&x, &effort, 1e-4, &mut loss)
                    .unwrap();
                assert!(reference.dissipative_flow_into(&x, &effort, &mut other));
                assert!((loss[1] - other[1]).abs() <= 1e-12 * other[1].abs().max(1e-30));
                assert!(rate * loss[1] >= 0.0);
            }
        }
        let mut out = [0.0; 2];
        assert!(matches!(
            storage.try_hamiltonian(&[source.envelope().maximum_indentation_m.next_up(), 0.0]),
            Err(BoundNormalContactError::Source(
                NormalLawSourceError::Point(_)
            ))
        ));
        assert!(
            storage
                .dissipative_flow_into(&[1e-5, 0.0], &[0.0, 0.201], 1e-4, &mut out)
                .is_err()
        );
        assert!(
            storage
                .dissipative_flow_into(&[1e-5, 0.0], &[0.0, 0.0], 0.002, &mut out)
                .is_err()
        );
    }

    #[test]
    fn g0_bound_contact_refuses_missing_maps_and_mismatched_states() {
        let point = source_fixture::source(1e7, None);
        assert!(bound_point_obstacle(&point, vec![], 0.0, "fixture").is_err());
        assert!(bound_point_obstacle(&point, vec![1.0], 0.0, "").is_err());
        let source = NormalLawSource::from(point);
        assert!(BoundNormalContact::new(source.clone(), vec![], 0.0).is_err());
        assert!(BoundNormalContact::new(source.clone(), vec![0.0], 0.0).is_err());
        assert!(BoundNormalContact::new(source.clone(), vec![f64::NAN], 0.0).is_err());
        let site = BoundNormalContact::new(source, vec![1.0], 0.0).unwrap();
        assert!(BoundNormalStorage::new(Box::new(UnitMass), 2, vec![site.clone()]).is_err());
        let storage = BoundNormalStorage::new(Box::new(UnitMass), 1, vec![site]).unwrap();
        assert!(storage.try_hamiltonian(&[0.0]).is_err());
        assert!(storage.try_gradient(&[0.0, 0.0], &mut [0.0]).is_err());
    }
}
