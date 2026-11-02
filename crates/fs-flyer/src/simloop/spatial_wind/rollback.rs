//! Reusable dynamic-state backup for the opt-in homogeneous domain restriction.
//! Immutable flight models and spectral factors are not copied every tick.

use super::{SimLoop, WindField, map_atmo};
use crate::Refusal;
use crate::perception::PerceptionState;
use crate::pilot::PilotState;
use crate::simloop::Phase;
use fs_atmo::mann_field::HomogeneousMannField;

pub(super) struct Rollback {
    tick: u64,
    phase: Phase,
    flight: [f64; 14],
    warm_slip: Option<[f64; 2]>,
    envelope_refusal: Option<Refusal>,
    digest_acc: [u8; 32],
    perception: PerceptionState,
    pilot: Option<PilotState>,
    wind: HomogeneousMannField,
    ground_mps: [f64; 2],
}

impl Rollback {
    pub fn new(sim: &SimLoop) -> Self {
        let spatial = sim.spatial_wind.as_ref().expect("homogeneous wind");
        let WindField::Homogeneous(wind) = &spatial.field else {
            unreachable!("rollback is enabled only for homogeneous wind")
        };
        Self {
            tick: sim.tick,
            phase: sim.phase,
            flight: Self::flight(sim),
            warm_slip: sim.warm_slip,
            envelope_refusal: sim.envelope_refusal.clone(),
            digest_acc: sim.digest_acc,
            perception: sim.perception_state.clone(),
            pilot: sim.pilot.as_ref().map(|(_, state)| state.clone()),
            wind: wind.clone(),
            ground_mps: spatial.ground_mps,
        }
    }

    fn flight(sim: &SimLoop) -> [f64; 14] {
        [
            sim.x_m,
            sim.h_m,
            sim.u,
            sim.w,
            sim.q,
            sim.theta,
            sim.omega,
            sim.warp,
            sim.mech.delta_rad,
            sim.mech.rate_rad_s,
            sim.lateral_state.phi_rad,
            sim.lateral_state.p_rad_s,
            sim.lateral_state.psi_rad,
            sim.lateral_state.r_rad_s,
        ]
    }

    pub fn capture(&mut self, sim: &SimLoop) -> Result<(), Refusal> {
        let spatial = sim.spatial_wind.as_ref().expect("homogeneous wind");
        let WindField::Homogeneous(wind) = &spatial.field else {
            unreachable!("rollback is enabled only for homogeneous wind")
        };
        self.wind.copy_temporal_state_from(wind).map_err(map_atmo)?;
        self.tick = sim.tick;
        self.phase = sim.phase;
        self.flight = Self::flight(sim);
        self.warm_slip = sim.warm_slip;
        self.envelope_refusal.clone_from(&sim.envelope_refusal);
        self.digest_acc = sim.digest_acc;
        self.perception.copy_state_from(&sim.perception_state);
        if let (Some(target), Some((_, source))) = (&mut self.pilot, &sim.pilot) {
            target.copy_state_from(source);
        }
        self.ground_mps = spatial.ground_mps;
        Ok(())
    }

    pub fn restore(&mut self, sim: &mut SimLoop) {
        sim.tick = self.tick;
        sim.phase = self.phase;
        [
            sim.x_m,
            sim.h_m,
            sim.u,
            sim.w,
            sim.q,
            sim.theta,
            sim.omega,
            sim.warp,
            sim.mech.delta_rad,
            sim.mech.rate_rad_s,
            sim.lateral_state.phi_rad,
            sim.lateral_state.p_rad_s,
            sim.lateral_state.psi_rad,
            sim.lateral_state.r_rad_s,
        ] = self.flight;
        sim.warm_slip = self.warm_slip;
        core::mem::swap(&mut sim.envelope_refusal, &mut self.envelope_refusal);
        sim.digest_acc = self.digest_acc;
        core::mem::swap(&mut sim.perception_state, &mut self.perception);
        if let (Some(source), Some((_, target))) = (&mut self.pilot, &mut sim.pilot) {
            core::mem::swap(source, target);
        }
        let spatial = sim.spatial_wind.as_mut().expect("homogeneous wind");
        let WindField::Homogeneous(wind) = &mut spatial.field else {
            unreachable!("rollback is enabled only for homogeneous wind")
        };
        core::mem::swap(wind, &mut self.wind);
        spatial.ground_mps = self.ground_mps;
    }
}

impl SimLoop {
    pub(in crate::simloop) fn step_with_wind_domain(
        &mut self,
        input: Option<crate::simloop::ControlInput>,
    ) -> Result<crate::simloop::SimStateOut, Refusal> {
        let stored = self
            .spatial_wind
            .as_mut()
            .expect("homogeneous wind")
            .rollback
            .take();
        let mut rollback = stored.unwrap_or_else(|| Box::new(Rollback::new(self)));
        if let Err(error) = rollback.capture(self) {
            self.spatial_wind
                .as_mut()
                .expect("homogeneous wind")
                .rollback = Some(rollback);
            return Err(error);
        }
        let result = self.step_inner(input);
        if result.is_err() {
            rollback.restore(self);
        }
        self.spatial_wind
            .as_mut()
            .expect("homogeneous wind")
            .rollback = Some(rollback);
        result
    }
}
