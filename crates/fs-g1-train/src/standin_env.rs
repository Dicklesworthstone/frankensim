//! Kinematic stand-in environment for the G1 walker — EXACT port of
//! cmaes_explainer `app/lib/g1StepwiseEnv.ts` (f64 state, identical reward
//! decomposition and observation layout), so policies trained here face the
//! same dynamics the browser evaluates them on.
//!
//! Disclosed scope: this is the stepwise-explainer stand-in, NOT the v068
//! kernel SE(3) dynamics. Its optimum is a near-zero action (any effort
//! torque induces tilt), which makes it a clean sample-efficiency benchmark.

use crate::ppo::G1Env;

pub const STANDIN_DT: f64 = 1.0 / 60.0;
pub const STANDIN_GAIT_FREQ_HZ: f64 = 1.5;
pub const STANDIN_TARGET_SPEED_MPS: f64 = 0.65;
pub const STANDIN_FALL_HEIGHT_M: f64 = 0.40;
pub const STANDIN_FALL_TILT_RAD: f64 = 0.85;
pub const STANDIN_INITIAL_HEIGHT_M: f64 = 0.75;

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
            let act = action.get(j).copied().unwrap_or(0.0) as f64;
            let target_pos = act * 0.5;
            let torque = 120.0 * (target_pos - self.joint_pos[j]) - 8.0 * self.joint_vel[j];
            self.joint_vel[j] += (torque / 1.5) * dt;
            self.joint_pos[j] += self.joint_vel[j] * dt;
            step_work += (torque * self.joint_vel[j]).abs() * dt;
            action_effort += act.abs();
        }
        self.cumulative_work += step_work;

        // Damped-pendulum pelvis: effort torque tilts, 0.6 damping recovers.
        let effort_torque = 0.04 * action_effort * action_effort;
        let angular_vel_roll = effort_torque - 0.6 * self.roll;
        let angular_vel_pitch = effort_torque - 0.6 * self.pitch;
        self.roll += angular_vel_roll * dt;
        self.pitch += angular_vel_pitch * dt;

        // Height sinks as tilt grows; tilt also gates forward progress.
        let tilt_sq = self.roll * self.roll + self.pitch * self.pitch;
        self.height -= 0.5 * tilt_sq * dt;
        let tilt = self.roll.hypot(self.pitch);
        let upright_factor = (1.0 - tilt).max(0.0);
        let delta_x = STANDIN_TARGET_SPEED_MPS * dt * upright_factor;
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

    /// Zero action is the disclosed near-optimum: no effort torque, no tilt,
    /// upright factor 1 ⇒ exact analytic distance and reward over 240 steps.
    #[test]
    fn zero_action_walks_at_target_speed() {
        let mut env = StandinEnv::new(240);
        let _obs = env.reset(0);
        let zero = [0.0f32; 15];
        let mut total = 0.0f32;
        for step in 1..=240 {
            let (_obs, r, done) = env.step(&zero);
            total += r;
            // done must fire ONLY on the horizon step (timeout), never a fall.
            assert!(step == 240 || !done, "early done at step {step}");
        }
        let expected_dx_per_step = STANDIN_TARGET_SPEED_MPS * STANDIN_DT;
        assert!(
            (env.cumulative_distance() - 240.0 * expected_dx_per_step).abs() < 1e-9,
            "distance {}",
            env.cumulative_distance()
        );
        // reward/step = 15·Δx + 0.5·cos0·cos0 − 0 = 0.1625 + 0.5
        let expected_total = 240.0 * (15.0 * expected_dx_per_step + 0.5);
        // f32 accumulation of 240 terms drifts ~1e-4 — that is the env's
        // output precision, not a port bug (state/distance stay f64-exact).
        assert!(
            ((total as f64) - expected_total).abs() < 1e-4,
            "total reward {total} vs {expected_total}"
        );
    }

    /// A large constant action must drive tilt past the fall threshold.
    #[test]
    fn violent_action_falls() {
        let mut env = StandinEnv::new(240);
        let _obs = env.reset(0);
        let violent = [1.0f32; 15];
        let mut fell = false;
        for _ in 0..240 {
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
