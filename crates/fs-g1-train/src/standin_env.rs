//! Kinematic stand-in environment for the G1 walker — exact port of the
//! default cmaes_explainer `app/lib/g1StepwiseEnv.ts` contract (f64 state,
//! identical reward decomposition and observation layout), so policies
//! trained here face the same default dynamics the browser evaluates them on.
//!
//! Disclosed scope: this is the stepwise-explainer stand-in, NOT the v068
//! kernel SE(3) dynamics. Forward displacement is an action-causal kinematic
//! proxy derived from lower-body joint motion and bilateral hip opposition;
//! zero action produces exactly zero locomotion distance.

use crate::ppo::G1Env;

pub const STANDIN_DT: f64 = 1.0 / 60.0;
pub const STANDIN_GAIT_FREQ_HZ: f64 = 1.5;
pub const STANDIN_TARGET_SPEED_MPS: f64 = 0.65;
pub const STANDIN_FALL_HEIGHT_M: f64 = 0.40;
pub const STANDIN_FALL_TILT_RAD: f64 = 0.85;
pub const STANDIN_INITIAL_HEIGHT_M: f64 = 0.75;
/// Identity shared with the browser-side action-causal stand-in.
pub const STANDIN_CONTRACT_ID: &str = "action-causal-standin-v2";

/// 42-D observation, same layout as the TS env's `rawVector`:
/// jointPos(15), jointVel(15), roll, pitch, yaw, ω(3), accel(3), sinφ, cosφ, target.
pub const STANDIN_OBS_DIM: usize = 42;

pub struct StandinEnv {
    max_steps: usize,
    step_count: usize,
    joint_pos: [f64; 15],
    joint_vel: [f64; 15],
    roll: f64,
    pitch: f64,
    yaw: f64,
    height: f64,
    phase: f64,
    current_forward_speed: f64,
    cumulative_dist: f64,
    cumulative_reward: f64,
    cumulative_work: f64,
}

impl StandinEnv {
    pub fn new(max_steps: usize) -> Self {
        Self {
            max_steps,
            step_count: 0,
            joint_pos: [0.0; 15],
            joint_vel: [0.0; 15],
            roll: 0.0,
            pitch: 0.0,
            yaw: 0.0,
            height: STANDIN_INITIAL_HEIGHT_M,
            phase: 0.0,
            current_forward_speed: 0.0,
            cumulative_dist: 0.0,
            cumulative_reward: 0.0,
            cumulative_work: 0.0,
        }
    }

    pub fn cumulative_distance(&self) -> f64 {
        self.cumulative_dist
    }

    pub fn max_steps(&self) -> usize {
        self.max_steps
    }

    pub fn cumulative_reward(&self) -> f64 {
        self.cumulative_reward
    }

    pub fn cumulative_work_joules(&self) -> f64 {
        self.cumulative_work
    }

    pub fn current_forward_speed(&self) -> f64 {
        self.current_forward_speed
    }

    fn observation_f64(&self) -> [f64; STANDIN_OBS_DIM] {
        let omega_dot = 2.0 * std::f64::consts::PI * STANDIN_GAIT_FREQ_HZ;
        let omega_roll = self.phase.cos() * 0.04 * omega_dot;
        let omega_pitch = -(self.phase * 2.0).sin() * 0.06 * omega_dot;
        let mut obs = [0.0f64; STANDIN_OBS_DIM];
        let mut i = 0;
        for v in self.joint_pos.iter() {
            obs[i] = *v;
            i += 1;
        }
        for v in self.joint_vel.iter() {
            obs[i] = *v;
            i += 1;
        }
        obs[30] = self.roll;
        obs[31] = self.pitch;
        obs[32] = self.yaw;
        obs[33] = omega_roll;
        obs[34] = omega_pitch;
        obs[35] = 0.0; // omega_yaw
        obs[36] = 0.0;
        obs[37] = 0.0;
        obs[38] = -9.81; // linear accel
        obs[39] = self.phase.sin();
        obs[40] = self.phase.cos();
        obs[41] = STANDIN_TARGET_SPEED_MPS;
        obs
    }
}

impl G1Env for StandinEnv {
    fn reset(&mut self, _seed: u64) -> Vec<f32> {
        self.step_count = 0;
        self.joint_pos = [0.0; 15];
        self.joint_vel = [0.0; 15];
        self.roll = 0.0;
        self.pitch = 0.0;
        self.yaw = 0.0;
        self.height = STANDIN_INITIAL_HEIGHT_M;
        self.phase = 0.0;
        self.current_forward_speed = 0.0;
        self.cumulative_dist = 0.0;
        self.cumulative_reward = 0.0;
        self.cumulative_work = 0.0;
        self.observation_f64().iter().map(|v| *v as f32).collect()
    }

    fn step(&mut self, action: &[f32]) -> (Vec<f32>, f32, bool) {
        self.step_count += 1;
        let dt = STANDIN_DT;
        self.phase = (self.phase + 2.0 * std::f64::consts::PI * STANDIN_GAIT_FREQ_HZ * dt)
            % (2.0 * std::f64::consts::PI);

        // PD-tracked joints; work and effort accumulate exactly as in TS.
        let mut step_work = 0.0f64;
        let mut action_effort = 0.0f64;
        for j in 0..15 {
            let candidate = f64::from(action.get(j).copied().unwrap_or(0.0));
            let act = if candidate.is_finite() {
                candidate.clamp(-1.0, 1.0)
            } else {
                0.0
            };
            let target_pos = act * 0.5;
            let torque = 120.0 * (target_pos - self.joint_pos[j]) - 8.0 * self.joint_vel[j];
            self.joint_vel[j] += (torque / 1.5) * dt;
            self.joint_pos[j] += self.joint_vel[j] * dt;
            step_work += (torque * self.joint_vel[j]).abs() * dt;
            action_effort += act.abs();
        }
        self.cumulative_work += step_work;

        // Damped-pendulum pelvis: mean effort torque tilts, 0.6 damping
        // recovers. Normalizing by actuator count keeps an ordinary bounded
        // gait from behaving like fifteen simultaneous full-body kicks.
        let mean_action_effort = action_effort / 15.0;
        let effort_torque = 0.12 * mean_action_effort * mean_action_effort;
        let angular_vel_roll = effort_torque - 0.6 * self.roll;
        let angular_vel_pitch = effort_torque - 0.6 * self.pitch;
        self.roll += angular_vel_roll * dt;
        self.pitch += angular_vel_pitch * dt;

        // Height sinks as tilt grows. Forward displacement is action-causal:
        // a gait-drive proxy uses actual lower-body joint velocity plus
        // bilateral hip opposition. This is still only a disclosed kinematic
        // stand-in, but an inert policy can no longer earn locomotion.
        let tilt_sq = self.roll * self.roll + self.pitch * self.pitch;
        self.height -= 0.5 * tilt_sq * dt;
        let tilt = self.roll.hypot(self.pitch);
        let upright_factor = (1.0 - tilt).max(0.0);
        let leg_velocity_rms = (self.joint_vel[..12]
            .iter()
            .map(|velocity| velocity * velocity)
            .sum::<f64>()
            / 12.0)
            .sqrt();
        let hip_opposition = ((self.joint_vel[0] - self.joint_vel[6]).abs() / 2.0).min(1.0);
        let motion_drive = (leg_velocity_rms / 2.0).min(1.0);
        let gait_drive = motion_drive * (0.25 + 0.75 * hip_opposition);
        let requested_forward_speed = STANDIN_TARGET_SPEED_MPS * gait_drive;
        let speed_response = (6.0 * dt).min(1.0);
        self.current_forward_speed +=
            (requested_forward_speed - self.current_forward_speed) * speed_response;
        let delta_x = self.current_forward_speed.max(0.0) * dt * upright_factor;
        self.cumulative_dist += delta_x;

        let fall = self.height < STANDIN_FALL_HEIGHT_M || tilt > STANDIN_FALL_TILT_RAD;
        let timeout = self.step_count >= self.max_steps;
        let done = fall || timeout;

        // Dense reward decomposition (identical weights to TS).
        let r_progress = 15.0 * delta_x;
        let r_upright = 0.5 * self.pitch.cos() * self.roll.cos();
        let r_energy = -0.002 * step_work;
        let r_fall = if fall { -50.0 } else { 0.0 };
        let step_reward = r_progress + r_upright + r_energy + r_fall;
        self.cumulative_reward += step_reward;

        let obs = self.observation_f64().iter().map(|v| *v as f32).collect();
        (obs, step_reward as f32, done)
    }

    fn obs_dim(&self) -> usize {
        STANDIN_OBS_DIM
    }

    fn action_dim(&self) -> usize {
        15
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Planted negative: an inert policy cannot receive commanded progress.
    #[test]
    fn zero_action_cannot_earn_locomotion_distance() {
        let mut env = StandinEnv::new(720);
        let _obs = env.reset(0);
        let zero = [0.0f32; 15];
        for step in 1..=720 {
            let (_obs, _reward, done) = env.step(&zero);
            assert!(step == 720 || !done, "early done at step {step}");
        }
        assert_eq!(env.cumulative_distance(), 0.0);
        assert_eq!(env.current_forward_speed(), 0.0);
        assert_eq!(env.cumulative_work_joules(), 0.0);
    }

    /// Planted positive: phase-opposed leg motion must diverge from the
    /// no-action control through the same causal path used by the browser.
    #[test]
    fn phase_opposed_leg_actions_cause_forward_progress() {
        let mut active = StandinEnv::new(180);
        let mut idle = StandinEnv::new(180);
        let mut active_obs = active.reset(7);
        let _idle_obs = idle.reset(7);
        for _ in 0..180 {
            let phase_sin = active_obs[39];
            let mut action = [0.0f32; 15];
            action[0] = 0.55 * phase_sin;
            action[3] = -0.4 * phase_sin;
            action[6] = -0.55 * phase_sin;
            action[9] = 0.4 * phase_sin;
            let (next, _reward, _done) = active.step(&action);
            let _ = idle.step(&[0.0; 15]);
            active_obs = next;
        }
        assert!(
            active.cumulative_distance() > 0.05,
            "active distance {}",
            active.cumulative_distance()
        );
        assert_eq!(idle.cumulative_distance(), 0.0);
        assert!(active.cumulative_distance() > idle.cumulative_distance());
    }

    /// A large constant action must drive tilt past the fall threshold.
    #[test]
    fn violent_action_falls() {
        let mut env = StandinEnv::new(720);
        let _obs = env.reset(0);
        let violent = [1.0f32; 15];
        let mut fell = false;
        for _ in 0..720 {
            let (_obs, _r, done) = env.step(&violent);
            if done {
                fell = true;
                break;
            }
        }
        assert!(fell, "constant full-scale action must trip the fall gate");
    }

    /// Observation layout spot-checks against the TS rawVector contract.
    #[test]
    fn observation_layout_matches_ts() {
        let mut env = StandinEnv::new(10);
        let obs = env.reset(0);
        assert_eq!(obs.len(), STANDIN_OBS_DIM);
        assert_eq!(obs[30], 0.0); // roll
        assert!((obs[38] - (-9.81)).abs() < 1e-6); // accel z
        assert!((obs[39] - 0.0).abs() < 1e-6); // sin(0)
        assert!((obs[40] - 1.0).abs() < 1e-6); // cos(0)
        assert!((obs[41] - STANDIN_TARGET_SPEED_MPS as f32).abs() < 1e-6); // target
        let (obs2, _, _) = env.step(&[0.0f32; 15]);
        // phase advanced by 2π·1.5/60 ⇒ sin positive
        assert!(obs2[39] > 0.0);
    }
}
