//! Opt-in spatial wind for the longitudinal lifecycle. The ground velocity is
//! inertial state; air-relative velocity is sampled, never time-differentiated.
//! This avoids treating an OU innovation as an impulsive aircraft acceleration.

use super::{SimLoop, map_atmo, motion};
use crate::Refusal;
use crate::atmosphere::MaterialAirProfile;
use fs_atmo::evolving::{EvolvingTurbulenceCheckpoint, EvolvingTurbulenceField};
use fs_atmo::ou::OuCheckpoint;
use fs_atmo::{FlatSiteLogLaw, FlowSample, MAX_MODES, TurbulenceField};
use fs_blake3::hash_domain;
use fs_material::gas::GasState;
use fs_math::det;

pub(super) mod material;
pub use material::MaterialAirSample;

/// Explicit synthetic wind model for [`SimLoop::init_with_spatial_wind`].
/// All lengths are metres, speeds m/s, and correlation times seconds. The
/// correlation-time vector sets the number of spatial modes. Both spatial
/// coefficients and OU innovations use the scenario seed in distinct kernels.
/// The finite-mode spectrum is an analytic construction, not a measured or
/// validated reconstruction of historical turbulence.
#[derive(Clone, Debug, PartialEq)]
pub struct SpatialWindSpec {
    /// Neutral log-law mean. Its reference speed must equal the scenario's
    /// declared headwind; its reference height is explicit.
    pub mean: FlatSiteLogLaw,
    /// Spatial builder's velocity-strength parameter (not a certified RMS).
    pub sigma_mps: f64,
    /// Spatial builder's length-scale parameter.
    pub length_m: f64,
    /// Downwind phase-convection speed, independently specified.
    pub advection_mps: f64,
    /// One positive finite decorrelation time per spatial mode (1..=256).
    pub correlation_times_s: Vec<f64>,
}

pub(super) struct SpatialWindState {
    pub spec: SpatialWindSpec,
    pub field: EvolvingTurbulenceField,
    pub distributed: bool,
    pub material_gas: Option<GasState>,
    pub material_profile: Option<MaterialAirProfile>,
    /// Ground velocity in flight-line-forward / vertical-up coordinates.
    pub ground_mps: [f64; 2],
}

pub(super) struct SampledFlightLoads {
    pub forces: crate::aircraft::ForceBuildup,
    pub airspeed_mps: f64,
    pub alpha_rad: f64,
    /// The same CG density used by wing, propeller and lateral stage loads.
    pub rho_kg_m3: f64,
}

impl SpatialWindState {
    fn new(spec: SpatialWindSpec, seed: u64) -> Result<Self, Refusal> {
        spec.mean.admit().map_err(map_atmo)?;
        let basis = TurbulenceField::build(
            seed,
            spec.correlation_times_s.len(),
            spec.sigma_mps,
            spec.length_m,
            spec.advection_mps,
        )
        .map_err(map_atmo)?;
        let field =
            EvolvingTurbulenceField::at_anchor(basis, seed, spec.correlation_times_s.clone())
                .map_err(map_atmo)?;
        Ok(Self {
            spec,
            field,
            distributed: false,
            material_gas: None,
            material_profile: None,
            ground_mps: [0.0; 2],
        })
    }

    /// Flight-frame X is upwind, Y right, H up: x_atmo = -X. Transform
    /// both velocity and derivative coordinates, including mean shear.
    pub fn sample(&self, [x, y, h]: [f64; 3]) -> Result<FlowSample, Refusal> {
        let mut sample = self.field.sample(-x, y, h).map_err(map_atmo)?;
        sample.u[0] += self.spec.mean.speed(h);
        sample.grad[0][2] += self.spec.mean.dspeed_dh(h);
        let sign = [-1.0, 1.0, 1.0];
        for i in 0..3 {
            sample.u[i] *= sign[i];
            for j in 0..3 {
                sample.grad[i][j] *= sign[i] * sign[j];
            }
        }
        Ok(sample)
    }

    fn model_bytes(&self, bytes: &mut Vec<u8>) {
        let mean = self.spec.mean;
        for value in [
            mean.scenario_effective_z0_m,
            mean.displacement_height_m,
            mean.reference_height_m,
            mean.reference_speed_mps,
            self.spec.sigma_mps,
            self.spec.length_m,
            self.spec.advection_mps,
        ] {
            bytes.extend_from_slice(&value.to_bits().to_le_bytes());
        }
        bytes.extend_from_slice(&(self.spec.correlation_times_s.len() as u32).to_le_bytes());
        for value in &self.spec.correlation_times_s {
            bytes.extend_from_slice(&value.to_bits().to_le_bytes());
        }
    }

    pub fn encode(&self, run_intent_id: &str, bytes: &mut Vec<u8>) {
        self.model_bytes(bytes);
        let cp = self.field.checkpoint();
        bytes.extend_from_slice(&cp.cosine.seed.to_le_bytes());
        bytes.extend_from_slice(&cp.cosine.tick.to_le_bytes());
        for value in self
            .ground_mps
            .iter()
            .chain(&cp.cosine.amplitudes)
            .chain(&cp.sine.amplitudes)
        {
            bytes.extend_from_slice(&value.to_bits().to_le_bytes());
        }
        bytes.extend_from_slice(run_intent_id.as_bytes()); // fixed 64-byte hex
    }

    pub fn decode(bytes: &[u8], seed: u64, tick: u64) -> Result<(Self, String), Refusal> {
        let bad = || invalid("malformed spatial-wind checkpoint");
        let mut pos = 0;
        let mut take = |n: usize| {
            let slice = bytes.get(pos..pos + n).ok_or_else(bad)?;
            pos += n;
            Ok::<_, Refusal>(slice)
        };
        let mut f = [0.0; 7];
        for value in &mut f {
            *value = f64::from_le_bytes(take(8)?.try_into().expect("8"));
        }
        let count = u32::from_le_bytes(take(4)?.try_into().expect("4")) as usize;
        if !(1..=MAX_MODES).contains(&count) {
            return Err(bad());
        }
        let mut times = Vec::with_capacity(count);
        for _ in 0..count {
            times.push(f64::from_le_bytes(take(8)?.try_into().expect("8")));
        }
        let stored_seed = u64::from_le_bytes(take(8)?.try_into().expect("8"));
        let stored_tick = i64::from_le_bytes(take(8)?.try_into().expect("8"));
        if stored_seed != seed || u64::try_from(stored_tick).ok() != Some(tick) {
            return Err(invalid("atmosphere seed or tick disagrees with the run"));
        }
        let mut ground_mps = [0.0; 2];
        for value in &mut ground_mps {
            *value = f64::from_le_bytes(take(8)?.try_into().expect("8"));
        }
        if !ground_mps.iter().all(|v| v.is_finite()) {
            return Err(bad());
        }
        let mut amplitudes = [Vec::with_capacity(count), Vec::with_capacity(count)];
        for quadrature in &mut amplitudes {
            for _ in 0..count {
                quadrature.push(f64::from_le_bytes(take(8)?.try_into().expect("8")));
            }
        }
        let run_intent_id = std::str::from_utf8(take(64)?)
            .map_err(|_| bad())?
            .to_owned();
        if pos != bytes.len() {
            return Err(bad());
        }
        let mut state = Self::new(
            SpatialWindSpec {
                mean: FlatSiteLogLaw {
                    scenario_effective_z0_m: f[0],
                    displacement_height_m: f[1],
                    reference_height_m: f[2],
                    reference_speed_mps: f[3],
                },
                sigma_mps: f[4],
                length_m: f[5],
                advection_mps: f[6],
                correlation_times_s: times,
            },
            seed,
        )?;
        let model = state.field.checkpoint();
        let [cosine, sine] = amplitudes;
        state.field = EvolvingTurbulenceField::resume(EvolvingTurbulenceCheckpoint {
            cosine: OuCheckpoint {
                seed: model.cosine.seed,
                tick: stored_tick,
                amplitudes: cosine,
            },
            sine: OuCheckpoint {
                seed: model.sine.seed,
                tick: stored_tick,
                amplitudes: sine,
            },
            ..model
        })
        .map_err(map_atmo)?;
        state.ground_mps = ground_mps;
        Ok((state, run_intent_id))
    }
}

impl SimLoop {
    /// Opt into local wing/canard strip and separate propeller disk inflow.
    /// Uses the same evolving field, time convention, and inertial ground
    /// state as [`Self::init_with_spatial_wind`], but local incident velocity
    /// and prop wash now feed actual circulation, forces, pitch moment, drag,
    /// and rotor torque. Pitch-rate motion enters each sample exactly once.
    /// This remains quasi-steady strip aerodynamics with a shared straight
    /// wake, axial disk-center BEMT, and reduced longitudinal flight motion.
    /// Accepted distributed force moments also drive the existing reduced
    /// roll/yaw response once airborne. Control and damping remain in the
    /// reduced lateral model; no banked translation, sideslip, rotor reaction
    /// torque, chordwise gust memory or disk harmonics are introduced.
    /// V4 checkpoints retain their layout; the v2 distributed model identity
    /// rejects older distributed trajectories that omitted roll/yaw forcing.
    /// V3 retains its original CG-only sampling and lateral behavior.
    ///
    /// # Errors
    /// As [`Self::init_with_spatial_wind`]. Local-flow/ground-plane refusals
    /// during stepping end the run with the existing envelope receipt.
    pub fn init_with_distributed_spatial_wind(
        scenario: super::ScenarioSpec,
        wind: SpatialWindSpec,
    ) -> Result<Self, Refusal> {
        let mut sim = Self::init_with_spatial_wind(scenario, wind)?;
        sim.spatial_wind
            .as_mut()
            .expect("initialized spatial field")
            .distributed = true;
        sim.run_intent_id = hash_domain(
            "org.frankensim.wf.distributed-wind-intent.v2",
            sim.run_intent_id.as_bytes(),
        )
        .to_hex();
        Ok(sim)
    }

    /// Opt into a spatially coherent, evolving wind sampled at the CG. This
    /// replaces the reduced scalar OU gust; it is never added a second time.
    /// Longitudinal aerodynamic/propeller loads use the sampled relative u,w
    /// and the scenario density. Crosswind remains available to observers but
    /// is not a sideslip load in the existing reduced lateral model. The
    /// mechanics retain the baseline steady-trim initialization; this is not
    /// a gust-equilibrated historical launch or a distributed spanwise solver.
    ///
    /// # Errors
    /// Scenario, mean-law, modal, and correlation-time admission refusals;
    /// `spatial-wind-invalid` if the two declared reference headwinds differ.
    pub fn init_with_spatial_wind(
        scenario: super::ScenarioSpec,
        wind: SpatialWindSpec,
    ) -> Result<Self, Refusal> {
        if wind.mean.reference_speed_mps != scenario.headwind_mps {
            return Err(invalid(
                "mean reference speed must equal the scenario headwind",
            ));
        }
        let mut state = SpatialWindState::new(wind, scenario.seed)?;
        state.field.advance_to(0).map_err(map_atmo)?;
        let mut sim = Self::init(scenario)?;
        let mut initial = Vec::new();
        state.encode("", &mut initial);
        initial.extend_from_slice(sim.tick0.digest.as_bytes());
        sim.tick0.digest =
            hash_domain("org.frankensim.wf.spatial-wind-tick0.v1", &initial).to_hex();
        let cp = state.field.checkpoint();
        sim.tick0.ou_amplitudes = cp.cosine.amplitudes;
        sim.tick0.ou_amplitudes.extend(cp.sine.amplitudes);
        let mut intent = sim.run_intent_id.as_bytes().to_vec();
        intent.extend_from_slice(sim.tick0.digest.as_bytes());
        intent.extend_from_slice(&sim.spec.max_ticks.to_le_bytes());
        intent.extend_from_slice(&sim.spec.rho_kg_m3.to_bits().to_le_bytes());
        if let Some(assist) = sim.spec.assist {
            for value in [
                assist.sas_rate_gain,
                assist.assist_attitude_gain,
                assist.authority_frac,
            ] {
                intent.extend_from_slice(&value.to_bits().to_le_bytes());
            }
        }
        sim.run_intent_id =
            hash_domain("org.frankensim.wf.spatial-wind-intent.v1", &intent).to_hex();
        sim.digest_acc = *hash_domain(
            "org.frankensim.wf.sim-digest.v1",
            sim.tick0.digest.as_bytes(),
        )
        .as_bytes();
        sim.spatial_wind = Some(state);
        sim.sync_spatial_air_velocity()?;
        Ok(sim)
    }

    /// Query the exact current physics wind without advancing its clock or
    /// random streams. Coordinates/components are [flight-line forward,
    /// right, vertical up], in metres and m/s; gradient columns use the same
    /// coordinates. Arbitrary lateral observer positions are supported even
    /// though the reduced aircraft trajectory remains at y = 0.
    ///
    /// # Errors
    /// `spatial-wind-disabled` for the legacy lifecycle; fs-atmo query refusals.
    pub fn sample_spatial_wind(&self, position_m: [f64; 3]) -> Result<FlowSample, Refusal> {
        self.spatial_wind
            .as_ref()
            .ok_or_else(|| Refusal {
                code: "spatial-wind-disabled",
                message: "this run uses the legacy reduced scalar gust".into(),
                ranked_repairs: vec!["initialize with an explicit SpatialWindSpec".into()],
            })?
            .sample(position_m)
    }

    /// Inertial [forward, upward] ground velocity for the opt-in model.
    #[must_use]
    pub fn spatial_ground_velocity_mps(&self) -> Option<[f64; 2]> {
        self.spatial_wind.as_ref().map(|wind| wind.ground_mps)
    }

    pub(super) fn sync_spatial_air_velocity(&mut self) -> Result<f64, Refusal> {
        let wind = self.spatial_wind.as_ref().expect("spatial path");
        let sample = wind.sample([self.x_m, 0.0, self.h_m])?;
        [self.u, self.w] =
            motion::body_air_velocity(wind.ground_mps, self.theta, [sample.u[0], sample.u[2]]);
        Ok(sample.u[2])
    }

    /// Translational relative flow at a body point. The aircraft assembly
    /// owns the additional q×r velocity; observer queries never advance time.
    #[cfg(test)]
    fn incident_at_body_point(&self, point: [f64; 3]) -> Result<[f64; 3], Refusal> {
        self.incident_at_motion(self.current_flight_motion(), point)
    }

    fn incident_at_motion(
        &self,
        state: motion::FlightMotion,
        point: [f64; 3],
    ) -> Result<[f64; 3], Refusal> {
        let wind = self
            .spatial_wind
            .as_ref()
            .expect("distributed spatial path");
        let r = [
            point[0] - self.design.cg_m[0],
            point[1] - self.design.cg_m[1],
            point[2] - self.design.cg_m[2],
        ];
        let theta = state.longitudinal.theta_rad;
        let (s, c) = (det::sin(theta), det::cos(theta));
        let position = [
            state.longitudinal.position_m[0] + c * r[0] + s * r[2],
            r[1],
            state.longitudinal.position_m[1] + s * r[0] - c * r[2],
        ];
        let sample = wind.sample(position)?;
        let vx = state.longitudinal.ground_mps[0] - sample.u[0];
        let vh = state.longitudinal.ground_mps[1] - sample.u[2];
        Ok([c * vx + s * vh, -sample.u[1], s * vx - c * vh])
    }

    #[cfg(test)]
    fn force_buildup_current_wind(
        &self,
        v_air: f64,
        alpha: f64,
        dc: f64,
    ) -> Result<crate::aircraft::ForceBuildup, Refusal> {
        self.sampled_loads_current_wind(v_air, alpha, dc)
            .map(|sample| sample.forces)
    }

    pub(super) fn sampled_loads_current_wind(
        &self,
        v_air: f64,
        alpha: f64,
        dc: f64,
    ) -> Result<SampledFlightLoads, Refusal> {
        self.sampled_loads_at_motion(
            self.current_flight_motion(),
            v_air,
            alpha,
            dc,
            self.warm_slip,
        )
    }

    #[cfg(test)]
    fn force_buildup_at_motion(
        &self,
        state: motion::FlightMotion,
        v_air: f64,
        alpha: f64,
        dc: f64,
        warm_slip: Option<[f64; 2]>,
    ) -> Result<crate::aircraft::ForceBuildup, Refusal> {
        self.sampled_loads_at_motion(state, v_air, alpha, dc, warm_slip)
            .map(|sample| sample.forces)
    }

    fn sampled_loads_at_motion(
        &self,
        state: motion::FlightMotion,
        v_air: f64,
        alpha: f64,
        dc: f64,
        warm_slip: Option<[f64; 2]>,
    ) -> Result<SampledFlightLoads, Refusal> {
        let rho_kg_m3 = self.air_density_at_height(state.longitudinal.position_m[1])?;
        let forces = if self
            .spatial_wind
            .as_ref()
            .is_some_and(|wind| wind.distributed)
        {
            self.design.force_buildup_with_incident(
                v_air,
                alpha,
                dc.clamp(-0.55, 0.55),
                state.rotor_rad_s.clamp(10.5, 119.5),
                state.longitudinal.q_rad_s.clamp(-2.9, 2.9),
                rho_kg_m3,
                warm_slip,
                &|point| self.incident_at_motion(state, point),
            )
        } else {
            self.design.force_buildup_warm(
                v_air,
                alpha,
                dc.clamp(-0.55, 0.55),
                state.rotor_rad_s.clamp(10.5, 119.5),
                state.longitudinal.q_rad_s.clamp(-2.9, 2.9),
                rho_kg_m3,
                warm_slip,
            )
        }?;
        Ok(SampledFlightLoads {
            forces,
            airspeed_mps: v_air,
            alpha_rad: alpha,
            rho_kg_m3,
        })
    }

    fn current_flight_motion(&self) -> motion::FlightMotion {
        motion::FlightMotion {
            longitudinal: motion::LongitudinalMotion {
                position_m: [self.x_m, self.h_m],
                ground_mps: self.spatial_wind.as_ref().map_or_else(
                    || {
                        motion::ground_velocity(
                            [self.u, self.w],
                            self.theta,
                            [-self.spec.headwind_mps, 0.0],
                        )
                    },
                    |wind| wind.ground_mps,
                ),
                theta_rad: self.theta,
                q_rad_s: self.q,
            },
            rotor_rad_s: self.omega,
        }
    }

    fn flight_acceleration(
        &self,
        state: motion::FlightMotion,
        loads: &crate::aircraft::ForceBuildup,
        alpha: f64,
    ) -> motion::FlightAcceleration {
        motion::FlightAcceleration {
            world_mps2: motion::world_acceleration(
                motion::nongravity_acceleration(
                    [loads.force_n[0], loads.force_n[2]],
                    alpha,
                    self.design.gross_mass_kg,
                ),
                state.longitudinal.theta_rad,
                motion::GRAVITY_MPS2,
            ),
            pitch_rad_s2: loads.moment_y_nm / super::IYY_KG_M2,
            rotor_rad_s2: loads.torque_imbalance_nm / 1.6,
        }
    }

    fn publish_flight_motion(&mut self, state: motion::FlightMotion) {
        [self.x_m, self.h_m] = state.longitudinal.position_m;
        self.theta = state.longitudinal.theta_rad;
        self.q = state.longitudinal.q_rad_s;
        self.omega = state.rotor_rad_s;
        if let Some(wind) = &mut self.spatial_wind {
            wind.ground_mps = state.longitudinal.ground_mps;
            // The sole endpoint clock advance also updates observed airspeed.
        } else {
            [self.u, self.w] = motion::body_air_velocity(
                state.longitudinal.ground_mps,
                self.theta,
                [-self.spec.headwind_mps, 0.0],
            );
        }
    }

    /// One extra native load solve at a predicted midpoint couples translation,
    /// pitch and rotor inertia. Control and temporal wind remain held during the
    /// tick, while spatial wind is sampled at the predicted material positions.
    /// A refused stage publishes neither predicted motion nor its warm start.
    /// `None` is the terminal fallback when the force point itself is at/below
    /// ground. A positive midpoint always gets its native force solve, and its
    /// accepted endpoint determines ground contact. Neither terminal branch
    /// claims smooth order or localized contact.
    pub(super) fn advance_airborne_motion(
        &mut self,
        start_loads: &crate::aircraft::ForceBuildup,
        start_alpha: f64,
        dc: f64,
        gust_w: f64,
        dt: f64,
    ) -> Result<Option<SampledFlightLoads>, Refusal> {
        let start = self.current_flight_motion();
        let initial_acceleration = self.flight_acceleration(start, start_loads, start_alpha);
        let stage = motion::midpoint_stage(start, initial_acceleration, dt);
        if !stage.is_finite() {
            return Err(invalid_flight_stage());
        }
        if stage.longitudinal.position_m[1] <= 0.0 {
            // Only this inadmissible force point needs the bounded held-load
            // terminal fallback. Do not preempt a valid midpoint because an
            // obsolete held-load endpoint would cross the ground.
            let mut terminal = motion::held_load_step(
                start.longitudinal,
                motion::nongravity_acceleration(
                    [start_loads.force_n[0], start_loads.force_n[2]],
                    start_alpha,
                    self.design.gross_mass_kg,
                ),
                initial_acceleration.pitch_rad_s2,
                motion::GRAVITY_MPS2,
                dt,
            );
            let terminal_rotor = start.rotor_rad_s + dt * initial_acceleration.rotor_rad_s2;
            if !(motion::FlightMotion {
                longitudinal: terminal,
                rotor_rad_s: terminal_rotor,
            })
            .is_finite()
            {
                return Err(invalid_flight_stage());
            }
            terminal.position_m[1] = 0.0;
            self.publish_flight_motion(motion::FlightMotion {
                longitudinal: terminal,
                rotor_rad_s: terminal_rotor,
            });
            self.phase = super::Phase::Ended(super::TerminalEvent::GroundContact);
            self.warm_slip = Some(start_loads.coupled.w_slip);
            return Ok(None);
        }
        let wind = if let Some(wind) = &self.spatial_wind {
            let [x, h] = stage.longitudinal.position_m;
            let sample = wind.sample([x, 0.0, h])?;
            [sample.u[0], sample.u[2]]
        } else {
            [-self.spec.headwind_mps, 0.0]
        };
        let [u, w] = motion::body_air_velocity(
            stage.longitudinal.ground_mps,
            stage.longitudinal.theta_rad,
            wind,
        );
        let airspeed_mps = det::sqrt(u * u + w * w).max(0.1);
        let alpha_rad = (det::atan2(w, u)
            + if self.spatial_wind.is_some() {
                0.0
            } else {
                gust_w / airspeed_mps
            })
        .clamp(-0.55, 0.55);
        let sample = self.sampled_loads_at_motion(
            stage,
            airspeed_mps,
            alpha_rad,
            dc,
            Some(start_loads.coupled.w_slip),
        )?;
        let mut next = motion::finish_midpoint(
            start,
            stage,
            self.flight_acceleration(stage, &sample.forces, alpha_rad),
            dt,
        );
        if !next.is_finite() {
            return Err(invalid_flight_stage());
        }
        let ground_contact = next.longitudinal.position_m[1] <= 0.0;
        if ground_contact {
            next.longitudinal.position_m[1] = 0.0;
        }
        // Admit the profile at the accepted endpoint before publishing motion.
        // A valid midpoint must not commit an out-of-domain atmospheric state.
        self.air_density_at_height(next.longitudinal.position_m[1])?;
        if ground_contact {
            self.phase = super::Phase::Ended(super::TerminalEvent::GroundContact);
        }
        self.publish_flight_motion(next);
        self.warm_slip = Some(sample.forces.coupled.w_slip);
        Ok(Some(sample))
    }

    /// Use the already accepted force solve, not a second wind sample after
    /// longitudinal integration. This callback is reached only on an airborne
    /// step; the rail still carries the lateral reaction. The distributed
    /// aerodynamic path has no warp or p/r damping, so adding its resolved
    /// torque to the existing reduced owners does not duplicate either law.
    pub(super) fn step_lateral_current_wind(
        &mut self,
        loaded_twist_rad: f64,
        v_air: f64,
        rho_kg_m3: f64,
        dc_eff: f64,
        loads: &crate::aircraft::ForceBuildup,
        dt: f64,
    ) -> Result<crate::lateral::LateralTick, Refusal> {
        if self.spatial_wind.as_ref().is_some_and(|wind| wind.distributed) {
            let incident = self.design.distributed_lateral_loads(
                loads, dc_eff.clamp(-0.55, 0.55), rho_kg_m3,
            )?.net();
            self.lateral_model.step_with_moments(
                &mut self.lateral_state, loaded_twist_rad, v_air,
                rho_kg_m3, dt, incident,
            ).map(|step| step.reduced)
        } else {
            self.lateral_model.step(
                &mut self.lateral_state, loaded_twist_rad, v_air,
                rho_kg_m3, dt,
            )
        }
    }

    pub(super) fn finish_spatial_tick(&mut self) -> Result<Option<f64>, Refusal> {
        if let Some(wind) = &mut self.spatial_wind {
            // A lateral envelope exit can coincide with crossing the ground.
            // Preserve that primary terminal event, but publish its endpoint
            // on the wall rather than refusing an atmosphere query below it.
            if matches!(self.phase, super::Phase::Ended(_)) && self.h_m < 0.0 {
                self.h_m = 0.0;
            }
            wind.field.advance_to(self.tick as i64).map_err(map_atmo)?;
            Ok(Some(self.sync_spatial_air_velocity()?))
        } else {
            Ok(None)
        }
    }

    /// Shared physical force projection for rail admission and constrained
    /// motion. Gravity is removed using the force assembly's incidence frame.
    pub(super) fn rail_world_forces(
        &self,
        loads: &crate::aircraft::ForceBuildup,
        alpha: f64,
    ) -> [f64; 2] {
        let [fx, fz] = crate::aircraft::nongravity_force_xz(
            [loads.force_n[0], loads.force_n[2]],
            alpha,
            self.design.gross_mass_kg,
        );
        let (s, c) = (det::sin(self.theta), det::cos(self.theta));
        [c * fx + s * fz, s * fx - c * fz]
    }

    /// Advance a rail that still carries positive load at the tick start.
    /// Already-unloaded states take the ordinary airborne path before this
    /// method is reached. Loaded ticks retain discrete constraint/tow timing;
    /// no within-tick support-loss localization is claimed.
    pub(super) fn advance_rail_motion(
        &mut self,
        loads: &crate::aircraft::ForceBuildup,
        [forward, upward]: [f64; 2],
        dt: f64,
    ) {
        let mass = self.design.gross_mass_kg;
        let weight = mass * motion::GRAVITY_MPS2;
        let normal = (weight - upward).max(0.0);
        let tow = match self.spec.catapult {
            Some(tow) if self.x_m < tow.pull_length_m => tow.pull_force_n,
            _ => 0.0,
        };
        let mut next = self.current_flight_motion();
        next.longitudinal.ground_mps[0] = (next.longitudinal.ground_mps[0]
            + dt * (forward + tow - super::RAIL_MU * normal) / mass)
            .max(0.0);
        next.longitudinal.ground_mps[1] = 0.0;
        next.longitudinal.position_m[0] += dt * next.longitudinal.ground_mps[0];
        next.rotor_rad_s += dt * loads.torque_imbalance_nm / 1.6;
        self.publish_flight_motion(next);
        if self.x_m > self.spec.rail_length_m {
            self.phase = super::Phase::Ended(super::TerminalEvent::RailEndWithoutLift);
        }
    }
}

fn invalid(message: &str) -> Refusal {
    Refusal {
        code: "spatial-wind-invalid",
        message: message.into(),
        ranked_repairs: vec!["use the original admitted wind model and matching checkpoint".into()],
    }
}

fn invalid_flight_stage() -> Refusal {
    Refusal {
        code: "flight-stage-invalid",
        message: "airborne load integration produced a nonfinite stage".into(),
        ranked_repairs: vec![
            "keep the state and held loads within the admitted flight envelope".into(),
        ],
    }
}

#[cfg(test)]
mod tests;

#[cfg(test)]
mod lateral_tests;
