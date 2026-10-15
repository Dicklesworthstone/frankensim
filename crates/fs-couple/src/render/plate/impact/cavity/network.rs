//! The existing distributed cavity in the reciprocal modal/friction network.
//!
//! This lowers the SAME volume springs and unit-mass acoustic/neck coordinates
//! used by `CavityCoupling::build`. It adds no constitutive law, assigned pitch,
//! delayed pressure drive or extra gas compliance. One collective spring per
//! pressure mode acts on the solid and air together, including a uniform mode
//! with no fictitious acoustic inertia. The receiving modal network owns the
//! complete solve, force/energy admission, friction, rollback and scheduling.
use super::{CavityCoupling, ImpactError, ModalAcousticState, invalid};
use crate::modal_acoustic_time::{ModalAcousticMode, ModalAcousticTimeBudget, ModalAcousticTimeModel};
use crate::render::schedule::force::coupled::{ModalAttachment, ModalConnection};
use crate::vibroacoustic::AcousticMedium;
use fs_exec::CancelGate;
use fs_math::c64::C64;

/// Original solid bank followed, when needed, by ONE bank of free air/neck
/// inertias. Connection indices are local to this list; the existing string-bank
/// assembler remaps them when it prepends strings. Retain every connection.
/// Passing only `components` to a solver would omit the cavity's physics.
#[derive(Debug)]
pub struct ModalCavityParts {
    /// Original structural component, then the appended inertias when nonempty.
    pub components: Vec<ModalAcousticTimeModel>,
    /// Collective compression springs and explicit acoustic-momentum dashpots.
    pub connections: Vec<ModalConnection>,
}

impl CavityCoupling {
    /// Frozen medium used by the cavity's pressure basis and every neck inertia.
    /// This is not a complete gas card or measured-composition receipt.
    #[must_use]
    pub fn acoustic_medium(&self) -> AcousticMedium { self.medium }

    /// Prepare a distributed body-air receiver for the existing coupled modal
    /// solver, including its implicit bow-friction image.
    ///
    /// `body` MUST use the same mass-normalized basis as this coupling's signed
    /// overlap rows. Its motion and observation coefficients are preserved.
    /// Standing-wave coordinates start at zero; explicit neck volume/flow is
    /// preserved. Nonzero body displacement therefore creates actual initial
    /// compression, not a silently relaxed cavity.
    ///
    /// `reference_area_m2` only converts each volume spring to the network's
    /// displacement/force units: shapes = volume_shapes / area and K = A*area².
    /// It is NOT a fitted material parameter and does not change the physics.
    /// Free acoustic coordinates have zero isolated frequency and damping;
    /// their entire stiffness and momentum loss belong to the returned links.
    /// The original neck compactness guards still apply. Necks requiring an
    /// unconnected exterior radiation load refuse rather than silently losing it.
    ///
    /// # Errors
    /// Mismatched basis/clock, nonrepresentable scaling, acoustic Nyquist or
    /// state-budget refusal, missing exterior load, or cancellation. The final
    /// network must still admit its FULL coupled bandwidth, energy and work.
    pub fn modal_network(
        &self, body: ModalAcousticTimeModel, sample_rate_hz: u32,
        reference_area_m2: f64, acoustic_budget: ModalAcousticTimeBudget,
        gate: &CancelGate,
    ) -> Result<ModalCavityParts, ImpactError> {
        if gate.is_requested() { return Err(ImpactError::Cancelled); }
        if body.modes().len() != self.structural || sample_rate_hz == 0
            || body.sample_period_s().to_bits() != f64::from(sample_rate_hz).recip().to_bits()
            || !reference_area_m2.is_finite() || reference_area_m2 <= 0.0
            || self.necks.iter().any(|neck| neck.exterior_load)
        {
            return Err(invalid("modal cavity requires its exact structural basis/clock, positive area scale and complete neck loads"));
        }
        let dt = body.sample_period_s();
        if self.omegas.iter().any(|w| w*dt >= 0.9*core::f64::consts::PI)
            || self.necks.iter().any(|n| n.fixed_wall_omega*dt >= 0.9*core::f64::consts::PI)
        {
            return Err(invalid("modal cavity exceeds the original acoustic Nyquist guard"));
        }
        let extra = self.total - self.structural;
        let mut components = vec![body];
        let mut drag = vec![0.0; extra];
        if extra > 0 {
            let mut air = ModalAcousticTimeModel::try_new_with_free_coordinates(
                sample_rate_hz, vec![ModalAcousticMode {
                    angular_frequency_rad_s: 0.0, damping_ratio: 0.0,
                    pressure_per_modal_velocity: C64::ZERO,
                }; extra], acoustic_budget,
            ).map_err(|e| ImpactError::Owner(e.to_string()))?;
            let mut initial = vec![ModalAcousticState::default(); extra];
            for (&coordinate, &d) in self.dynamic.iter().zip(&self.damping) {
                if let Some(i) = coordinate { drag[i-self.structural] = d; }
            }
            for neck in &self.necks {
                let i = neck.coordinate-self.structural;
                initial[i] = neck.initial;
                drag[i] = neck.drag_per_s;
            }
            air.restore_states(&initial).map_err(|e| ImpactError::Owner(e.to_string()))?;
            components.push(air);
        }
        let mut connections = Vec::with_capacity(self.springs.len()+extra);
        for spring in &self.springs {
            if gate.is_requested() { return Err(ImpactError::Cancelled); }
            // Factor before squaring to avoid needless intermediate overflow.
            let root_k = (spring.bulk_modulus_pa/spring.volume_m3).sqrt()*reference_area_m2;
            let stiffness = root_k*root_k;
            let shapes: Vec<_> = spring.areas.iter().map(|v| v/reference_area_m2).collect();
            if !stiffness.is_finite() || stiffness <= 0.0
                || shapes.iter().zip(&spring.areas)
                    .any(|(v,original)| !v.is_finite() || (*original != 0.0 && *v == 0.0))
            {
                return Err(invalid("modal cavity spring/port scaling is unrepresentable"));
            }
            // The network measures LEFT MINUS RIGHT. Negating the air row
            // retains the original SUM of outward solid and gas displacements.
            let right = if extra == 0 {
                ModalAttachment { component: 0, shapes: vec![0.0; self.structural] }
            } else {
                ModalAttachment { component: 1,
                    shapes: shapes[self.structural..].iter().map(|v| -v).collect() }
            };
            connections.push(ModalConnection {
                left: ModalAttachment { component: 0, shapes: shapes[..self.structural].to_vec() },
                right, stiffness_n_m: stiffness, damping_n_s_m: 0.0, rest_extension_m: 0.0,
            });
        }
        for (i, &d) in drag.iter().enumerate() {
            if gate.is_requested() { return Err(ImpactError::Cancelled); }
            if d != 0.0 {
                let mut row = vec![0.0; extra]; row[i] = 1.0;
                connections.push(ModalConnection {
                    left: ModalAttachment { component: 1, shapes: row },
                    right: ModalAttachment { component: 1, shapes: vec![0.0; extra] },
                    stiffness_n_m: 0.0, damping_n_s_m: d, rest_extension_m: 0.0,
                });
            }
        }
        if gate.is_requested() { return Err(ImpactError::Cancelled); }
        Ok(ModalCavityParts { components, connections })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vibroacoustic::CavityModes;
    use super::super::neck::CavityNeck;

    fn cavity(dynamic: bool) -> CavityCoupling {
        let count = if dynamic { 2 } else { 1 };
        let basis = CavityModes { omegas: if dynamic { vec![0.0, 900.0] } else { vec![0.0] },
            lambdas: vec![0.01; count], interface: vec![vec![1.0]; count],
            rho0: 1.2, c0: 340.0, loss_factor: 0.0 };
        CavityCoupling::new(&basis, 2, if dynamic { &[0.02,0.005,-0.01,0.003] }
            else { &[0.02,-0.01] }, if dynamic { &[0.0,7.0] } else { &[0.0] }).unwrap()
    }
    fn body() -> ModalAcousticTimeModel {
        let mut body = ModalAcousticTimeModel::try_new(48_000,
            vec![ModalAcousticMode { angular_frequency_rad_s: 300.0,
                damping_ratio: 0.001, pressure_per_modal_velocity: C64::ZERO }; 2],
            ModalAcousticTimeBudget::audible_reference()).unwrap();
        body.restore_states(&[ModalAcousticState { displacement_m_sqrt_kg: 1.0e-5,
            velocity_m_sqrt_kg_per_s: 0.003 }, ModalAcousticState::default()]).unwrap();
        body
    }
    fn parts(cavity: &CavityCoupling, area: f64) -> ModalCavityParts {
        cavity.modal_network(body(), 48_000, area,
            ModalAcousticTimeBudget::audible_reference(), &CancelGate::new_clock_free()).unwrap()
    }
    fn row(link: &ModalConnection, counts: &[usize]) -> Vec<f64> {
        let mut row = vec![0.0; counts.iter().sum()];
        for (side, sign) in [(&link.left,1.0),(&link.right,-1.0)] {
            let offset: usize = counts[..side.component].iter().sum();
            for (i,v) in side.shapes.iter().enumerate() { row[offset+i] += sign*v; }
        }
        row
    }
    fn close(a: f64, b: f64) {
        assert!((a-b).abs() <= 2.0e-12*a.abs().max(b.abs()).max(1.0), "{a} != {b}");
    }

    #[test]
    fn uniform_compression_has_no_fictitious_inertia_and_keeps_initial_motion() {
        let cavity = cavity(false); let image = parts(&cavity, 0.01);
        assert_eq!(image.components.len(), 1);
        assert_eq!(image.components[0].states(), body().states());
        assert_eq!(image.connections.len(), 1);
        let r = row(&image.connections[0], &[2]);
        let q = [1.0e-5,-2.0e-5];
        let extension = r.iter().zip(q).map(|(r,q)|r*q).sum::<f64>();
        let volume = 0.02*q[0]-0.01*q[1];
        close(0.5*image.connections[0].stiffness_n_m*extension*extension,
            0.5*(1.2*340.0*340.0/0.01)*volume*volume);
    }

    #[test]
    fn distributed_springs_and_drag_match_original_storage_in_every_coordinate() {
        let cavity = cavity(true);
        for area in [0.001,0.1,2.0] {
            let image = parts(&cavity, area);
            assert_eq!(image.components.len(),2);
            assert_eq!(image.components[1].modes()[0].angular_frequency_rad_s,0.0);
            let counts = [2,1];
            let mut k = vec![0.0;9]; let mut d = vec![0.0;9];
            for link in &image.connections {
                let r = row(link,&counts);
                for i in 0..3 { for j in 0..3 {
                    k[3*i+j] += link.stiffness_n_m*r[i]*r[j];
                    d[3*i+j] += link.damping_n_s_m*r[i]*r[j];
                } }
            }
            for i in 0..3 { for j in 0..3 {
                let expected = cavity.springs.iter().map(|s|
                    (s.bulk_modulus_pa/s.volume_m3)*s.areas[i]*s.areas[j]).sum();
                close(k[3*i+j],expected);
                close(d[3*i+j],if i==2 && j==2 {7.0} else {0.0});
            } }
        }
    }

    #[test]
    fn neck_volume_flow_and_resistance_are_preserved_not_reset_or_doubled() {
        let gate = CancelGate::new_clock_free();
        let cavity = cavity(true).with_necks(vec![CavityNeck {
            area_m2: 1.0e-4, effective_length_m: 0.01, resistance_pa_s_m3: 2000.0,
            pressure_shape_averages: vec![1.0,0.2], initial_volume_m3: 2.0e-8,
            initial_flow_m3_s: -3.0e-6,
        }],&gate).unwrap();
        let image = parts(&cavity,0.01);
        assert_eq!(image.components[1].states().len(),2);
        let neck = &cavity.necks[0];
        assert_eq!(image.components[1].states()[1],neck.initial);
        let velocity = [0.0,0.0,0.01,neck.initial.velocity_m_sqrt_kg_per_s];
        let loss: f64 = image.connections.iter().map(|link| {
            let speed = row(link,&[2,2]).iter().zip(velocity).map(|(r,v)|r*v).sum::<f64>();
            link.damping_n_s_m*speed*speed
        }).sum();
        close(loss,7.0*0.01*0.01+2000.0*3.0e-6*3.0e-6);
    }

    #[test]
    fn mismatched_clock_area_and_acoustic_aliasing_refuse() {
        let gate = CancelGate::new_clock_free(); let budget = ModalAcousticTimeBudget::audible_reference();
        for area in [0.0,-1.0,f64::NAN,f64::INFINITY,f64::MIN_POSITIVE] {
            assert!(cavity(true).modal_network(body(),48_000,area,budget,&gate).is_err());
        }
        assert!(cavity(true).modal_network(body(),44_100,0.01,budget,&gate).is_err());
        let mut high = cavity(true); high.omegas[1] = 200_000.0;
        assert!(high.modal_network(body(),48_000,0.01,budget,&gate).is_err());
    }
}
