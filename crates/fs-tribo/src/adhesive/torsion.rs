//! Bonded circular contact under twist about its fixed normal.
//!
//! The elastic half-space solution is K_phi = (16/3) G_bar a^3,
//! 1/G_bar = 1/G_a + 1/G_b, M = -K_phi phi and U = K_phi phi^2/2.
//! This torsional reduced modulus differs from Mindlin's translational modulus.
//! See Wada, *Particle interaction model* (2010), equations 3.6--3.7 and
//! 3.27--3.28: <https://www.perc.it-chiba.ac.jp/~wada/particleinteraction/interactionmodel.zeta.2.0.pdf>;
//! the circular-punch stiffness is also equation 5.3 of Popov, Hess and Willert
//! (2019), <https://doi.org/10.1007/978-3-662-58709-6_5>.
//!
//! An independently supplied critical angle selects an elastic/perfectly plastic
//! constant-moment continuation. There is no universal atomistic threshold.
//! This is a lumped bonded-contact idealization, not the local Coulomb partial-
//! slip solution: that solution develops edge slip before global saturation.
//! The normal owner supplies one attached, frozen circular area and must screen
//! its applicability. Sliding, mixed-mode area change, detachment, frame
//! transport and heat partition are separate physical operations. One signed
//! increment is one monotone path; the exact return/work integration shares
//! the rolling law's allocation-free scalar kernel.

use fs_math::dd::Dd;

/// Explicit physical parameters, with no load cutoff or material-name defaults.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdhesiveTorsionParameters {
    /// Actual attached circular contact radius [m].
    pub contact_radius_m: f64,
    /// Shear moduli of the two original isotropic elastic bodies [Pa].
    pub shear_moduli_pa: [f64; 2],
    /// Independently calibrated onset of irreversible twist [rad].
    pub critical_twist_angle_rad: f64,
    /// Caller small-deformation ceiling for critical rim displacement/a.
    /// This equals the critical angle and must lie strictly between 0 and 1.
    pub maximum_critical_twist_angle_rad: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum AdhesiveTorsionError {
    InvalidInput { field: &'static str },
    OutsideSmallTwist { angle_rad: f64, limit_rad: f64 },
    StateMismatch,
    Unrepresentable { field: &'static str },
}
impl core::fmt::Display for AdhesiveTorsionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adhesive torsion refused: {self:?}")
    }
}
impl std::error::Error for AdhesiveTorsionError {}

/// Complete scalar memory bound to the original geometry, moduli and limits.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdhesiveTorsionState {
    parameters: AdhesiveTorsionParameters,
    elastic_fraction: Dd,
    angle: Dd,
}
impl AdhesiveTorsionState {
    pub fn elastic_twist_rad(self) -> f64 {
        self.parameters.critical_twist_angle_rad * self.elastic_fraction.to_f64()
    }
    /// Total signed twist, rounded once for observation; low bits remain retained.
    pub const fn accumulated_angle_rad(self) -> f64 {
        self.angle.to_f64()
    }
}

/// Resisting moments and exact work on the original monotone twist increment.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdhesiveTorsionStep {
    pub state: AdhesiveTorsionState,
    pub endpoint_moment_n_m: f64,
    /// Path-mean resisting moment; endpoint moment for a zero increment.
    pub mean_moment_n_m: f64,
    pub storage_before_j: f64,
    pub storage_after_j: f64,
    /// Factored change retains weak work beside a large stored energy.
    pub storage_change_j: f64,
    pub irreversible_loss_j: f64,
    /// Work entering the junction = storage change + irreversible loss [J].
    pub applied_work_j: f64,
    pub plastic_twist_increment_rad: f64,
}

/// Fixed circular patch with elastic twist and separately retained yield loss.
#[derive(Clone, Debug, PartialEq)]
pub struct AdhesiveTorsionContact {
    parameters: AdhesiveTorsionParameters,
    reduced_shear: f64,
    stiffness: f64,
    critical_moment: f64,
    yield_energy: f64,
}
impl AdhesiveTorsionContact {
    pub fn new(p: AdhesiveTorsionParameters) -> Result<Self, AdhesiveTorsionError> {
        for (field, value) in [
            ("attached radius", p.contact_radius_m),
            ("first shear modulus", p.shear_moduli_pa[0]),
            ("second shear modulus", p.shear_moduli_pa[1]),
            ("critical twist angle", p.critical_twist_angle_rad),
            ("small-twist limit", p.maximum_critical_twist_angle_rad),
        ] {
            if !value.is_finite() || value <= 0.0 {
                return Err(AdhesiveTorsionError::InvalidInput { field });
            }
        }
        if p.maximum_critical_twist_angle_rad >= 1.0 {
            return Err(AdhesiveTorsionError::InvalidInput {
                field: "small-twist limit",
            });
        }
        if p.critical_twist_angle_rad > p.maximum_critical_twist_angle_rad {
            return Err(AdhesiveTorsionError::OutsideSmallTwist {
                angle_rad: p.critical_twist_angle_rad,
                limit_rad: p.maximum_critical_twist_angle_rad,
            });
        }
        let smaller = p.shear_moduli_pa[0].min(p.shear_moduli_pa[1]);
        let larger = p.shear_moduli_pa[0].max(p.shear_moduli_pa[1]);
        let reduced_shear = positive(smaller / (1.0 + smaller / larger), "reduced shear modulus")?;
        // Avoid first cubing a small radius or multiplying two large moduli.
        let a = p.contact_radius_m;
        let stiffness = positive(
            ((16.0 / 3.0 * a) * (reduced_shear * a)) * a,
            "torsional stiffness",
        )?;
        let critical_moment = positive(stiffness * p.critical_twist_angle_rad, "yield moment")?;
        let yield_energy = positive(
            (0.5 * critical_moment) * p.critical_twist_angle_rad,
            "yield storage",
        )?;
        Ok(Self {
            parameters: p,
            reduced_shear,
            stiffness,
            critical_moment,
            yield_energy,
        })
    }
    pub const fn parameters(&self) -> AdhesiveTorsionParameters {
        self.parameters
    }
    pub const fn reduced_shear_modulus_pa(&self) -> f64 {
        self.reduced_shear
    }
    pub const fn stiffness_nm_per_rad(&self) -> f64 {
        self.stiffness
    }
    pub const fn critical_moment_n_m(&self) -> f64 {
        self.critical_moment
    }
    pub const fn initial_state(&self) -> AdhesiveTorsionState {
        AdhesiveTorsionState {
            parameters: self.parameters,
            elastic_fraction: Dd::ZERO,
            angle: Dd::ZERO,
        }
    }
    pub fn stored_energy_j(
        &self,
        state: &AdhesiveTorsionState,
    ) -> Result<f64, AdhesiveTorsionError> {
        self.admit(state)?;
        finite(
            (Dd::from_f64(self.yield_energy) * state.elastic_fraction * state.elastic_fraction)
                .to_f64(),
            "stored energy",
        )
    }
    fn admit(&self, state: &AdhesiveTorsionState) -> Result<(), AdhesiveTorsionError> {
        if state.parameters != self.parameters {
            return Err(AdhesiveTorsionError::StateMismatch);
        }
        Ok(())
    }
    /// Exact monotone increment, with signed unloading and yield reversal.
    /// No accepted state is mutated, including on arithmetic resolution refusal.
    pub fn advance(
        &self,
        old: &AdhesiveTorsionState,
        angle_increment_rad: f64,
    ) -> Result<AdhesiveTorsionStep, AdhesiveTorsionError> {
        self.admit(old)?;
        if !angle_increment_rad.is_finite() {
            return Err(AdhesiveTorsionError::InvalidInput {
                field: "angle increment",
            });
        }
        let step = super::moment::advance(
            old.elastic_fraction,
            old.angle,
            angle_increment_rad,
            self.parameters.critical_twist_angle_rad,
            self.critical_moment,
            self.yield_energy,
        )
        .map_err(|field| AdhesiveTorsionError::Unrepresentable { field })?;
        Ok(AdhesiveTorsionStep {
            state: AdhesiveTorsionState {
                parameters: self.parameters,
                elastic_fraction: step.fraction,
                angle: step.angle,
            },
            endpoint_moment_n_m: step.endpoint_moment,
            mean_moment_n_m: step.mean_moment,
            storage_before_j: step.storage_before,
            storage_after_j: step.storage_after,
            storage_change_j: step.storage_change,
            irreversible_loss_j: step.loss,
            applied_work_j: step.work,
            plastic_twist_increment_rad: finite(
                self.parameters.critical_twist_angle_rad * step.plastic_fraction.to_f64(),
                "plastic twist increment",
            )?,
        })
    }
}
fn finite(value: f64, field: &'static str) -> Result<f64, AdhesiveTorsionError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(AdhesiveTorsionError::Unrepresentable { field })
    }
}
fn positive(value: f64, field: &'static str) -> Result<f64, AdhesiveTorsionError> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(AdhesiveTorsionError::Unrepresentable { field })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn parameters() -> AdhesiveTorsionParameters {
        AdhesiveTorsionParameters {
            contact_radius_m: 2e-6,
            shear_moduli_pa: [3e9, 6e11],
            critical_twist_angle_rad: 0.002,
            maximum_critical_twist_angle_rad: 0.01,
        }
    }
    fn close(a: f64, b: f64) {
        assert!(
            (a - b).abs() <= 4e-13 * a.abs().max(b.abs()).max(1e-100),
            "{a} != {b}"
        );
    }
    fn balance(step: &AdhesiveTorsionStep, increment: f64) {
        close(
            step.applied_work_j,
            step.storage_change_j + step.irreversible_loss_j,
        );
        close(step.applied_work_j, -step.mean_moment_n_m * increment);
        assert!(step.irreversible_loss_j >= 0.0);
    }
    #[test]
    fn g1_two_half_space_compliances_add_and_elastic_work_returns() {
        let p = parameters();
        let law = AdhesiveTorsionContact::new(p).unwrap();
        let compliance = (3.0 / (16.0 * p.contact_radius_m.powi(3)))
            * (p.shear_moduli_pa[0].recip() + p.shear_moduli_pa[1].recip());
        close(law.stiffness_nm_per_rad(), compliance.recip());
        close(law.reduced_shear_modulus_pa(), 3e9 * (200.0 / 201.0));
        for sign in [-1.0, 1.0] {
            let angle = sign * 0.4 * p.critical_twist_angle_rad;
            let step = law.advance(&law.initial_state(), angle).unwrap();
            close(step.endpoint_moment_n_m, -angle / compliance);
            close(step.storage_after_j, 0.5 * angle * angle / compliance);
            assert_eq!(step.irreversible_loss_j, 0.0);
            balance(&step, angle);
            let unload = law.advance(&step.state, -angle).unwrap();
            close(unload.applied_work_j, -step.applied_work_j);
            assert_eq!(unload.storage_after_j, 0.0);
        }
    }
    #[test]
    fn g1_yield_reversal_and_closed_cycle_keep_storage_separate_from_loss() {
        let law = AdhesiveTorsionContact::new(parameters()).unwrap();
        let c = parameters().critical_twist_angle_rad;
        let moment = law.critical_moment_n_m();
        let loaded = law.advance(&law.initial_state(), 3.0 * c).unwrap();
        close(loaded.irreversible_loss_j, 2.0 * c * moment);
        balance(&loaded, 3.0 * c);
        let reversed = law.advance(&loaded.state, -6.0 * c).unwrap();
        close(reversed.endpoint_moment_n_m, moment);
        close(reversed.irreversible_loss_j, 4.0 * c * moment);
        close(reversed.storage_change_j, 0.0);
        balance(&reversed, -6.0 * c);
        let returned = law.advance(&reversed.state, 6.0 * c).unwrap();
        close(
            returned.state.elastic_twist_rad(),
            loaded.state.elastic_twist_rad(),
        );
        close(
            returned.applied_work_j + reversed.applied_work_j,
            8.0 * c * moment,
        );
    }
    #[test]
    fn g3_radius_cubed_scaling_body_order_and_subdivision() {
        let p = parameters();
        let law = AdhesiveTorsionContact::new(p).unwrap();
        let scaled = AdhesiveTorsionContact::new(AdhesiveTorsionParameters {
            contact_radius_m: 2.0 * p.contact_radius_m,
            ..p
        })
        .unwrap();
        close(
            scaled.stiffness_nm_per_rad(),
            8.0 * law.stiffness_nm_per_rad(),
        );
        let swapped = AdhesiveTorsionContact::new(AdhesiveTorsionParameters {
            shear_moduli_pa: [p.shear_moduli_pa[1], p.shear_moduli_pa[0]],
            ..p
        })
        .unwrap();
        assert_eq!(swapped.stiffness_nm_per_rad(), law.stiffness_nm_per_rad());
        let increment = p.critical_twist_angle_rad / 8.0;
        let one = law.advance(&law.initial_state(), 32.0 * increment).unwrap();
        let (mut state, mut loss, mut work) = (law.initial_state(), 0.0, 0.0);
        for _ in 0..32 {
            let step = law.advance(&state, increment).unwrap();
            state = step.state;
            loss += step.irreversible_loss_j;
            work += step.applied_work_j;
        }
        close(state.elastic_twist_rad(), one.state.elastic_twist_rad());
        close(
            state.accumulated_angle_rad(),
            one.state.accumulated_angle_rad(),
        );
        close(loss, one.irreversible_loss_j);
        close(work, one.applied_work_j);
    }
    #[test]
    fn g3_tiny_twist_is_retained_beyond_rounded_stored_energy() {
        let law = AdhesiveTorsionContact::new(parameters()).unwrap();
        let c = parameters().critical_twist_angle_rad;
        let start = law.advance(&law.initial_state(), 0.5 * c).unwrap().state;
        let tiny = c * 1e-20;
        let step = law.advance(&start, tiny).unwrap();
        assert_ne!(step.state, start);
        assert_eq!(step.storage_after_j, step.storage_before_j);
        assert!(step.storage_change_j > 0.0);
        balance(&step, tiny);
        assert_eq!(law.advance(&start, 0.0).unwrap().state, start);
        assert!(law.advance(&step.state, c * 1e-60).is_err());
    }
    #[test]
    fn g0_foreign_history_small_angle_and_invalid_input_refuse() {
        let p = parameters();
        let law = AdhesiveTorsionContact::new(p).unwrap();
        let state = law.initial_state();
        let other = AdhesiveTorsionContact::new(AdhesiveTorsionParameters {
            critical_twist_angle_rad: p.critical_twist_angle_rad * 2.0,
            ..p
        })
        .unwrap();
        assert_eq!(
            other.advance(&state, 0.0),
            Err(AdhesiveTorsionError::StateMismatch)
        );
        assert!(matches!(
            AdhesiveTorsionContact::new(AdhesiveTorsionParameters {
                critical_twist_angle_rad: 0.02,
                ..p
            }),
            Err(AdhesiveTorsionError::OutsideSmallTwist { .. })
        ));
        for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(
                AdhesiveTorsionContact::new(AdhesiveTorsionParameters {
                    shear_moduli_pa: [bad, 1.0],
                    ..p
                })
                .is_err()
            );
        }
        assert!(law.advance(&state, f64::NAN).is_err());
        assert!(law.advance(&state, f64::MAX).is_err());
        assert_eq!(state, law.initial_state());
    }
}
