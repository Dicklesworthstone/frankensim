//! ONNX export metadata: the browser needs a JSON sidecar alongside the
//! model file so the inference runtime knows the observation layout, action
//! scaling, and joint limits — the same honesty pattern as the WASM kernel's
//! schema envelope.

use serde::Serialize;

#[derive(Serialize)]
pub struct PolicyMetadata {
    pub model_name: String,
    pub architecture: String,
    pub param_count: usize,
    pub obs_layout: ObsLayout,
    pub action_scaling: ActionScaling,
    pub training_provenance: TrainingProvenance,
}

#[derive(Serialize)]
pub struct ObsLayout {
    pub total_dims: usize,
    pub joint_positions: (usize, usize),
    pub joint_velocities: (usize, usize),
    pub base_quaternion: (usize, usize),
    pub phase_sin_cos: (usize, usize),
    pub terrain_heights: Option<(usize, usize)>,
}

#[derive(Serialize)]
pub struct ActionScaling {
    pub action_dim: usize,
    pub lower_limits: Vec<f32>,
    pub upper_limits: Vec<f32>,
    pub tanh_scaled: bool,
}

#[derive(Serialize)]
pub struct TrainingProvenance {
    pub algorithm: String,
    pub optimizer: String,
    pub outer_hpo: String,
    pub total_training_steps: u64,
    pub final_mean_reward: f32,
    pub environments: usize,
    pub seed: u64,
    pub checkpoint_hash: String,
}

/// Generate the metadata JSON that accompanies the ONNX file.
pub fn generate_metadata(
    model_name: &str,
    param_count: usize,
    total_steps: u64,
    final_reward: f32,
    checkpoint_hash: &str,
) -> PolicyMetadata {
    PolicyMetadata {
        model_name: model_name.to_string(),
        architecture: "transformer_4L_256d_gqa_swiglu_rope".to_string(),
        param_count,
        obs_layout: ObsLayout {
            total_dims: 42,
            joint_positions: (0, 15),
            joint_velocities: (15, 30),
            base_quaternion: (30, 34),
            phase_sin_cos: (34, 36),
            terrain_heights: Some((36, 42)),
        },
        action_scaling: ActionScaling {
            action_dim: 29,
            lower_limits: vec![-2.0; 29],
            upper_limits: vec![2.0; 29],
            tanh_scaled: true,
        },
        training_provenance: TrainingProvenance {
            final_mean_reward: final_reward,
            algorithm: "PPO + GAE (lambda 0.95, gamma 0.99)".to_string(),
            optimizer: "Muon (hidden) + Adam (embed/head/norm)".to_string(),
            outer_hpo: "CMA-ES (1+lambda) over 8 hyperparameters".to_string(),
            total_training_steps: total_steps,
            environments: 4096,
            seed: 0x4731_5050,
            checkpoint_hash: checkpoint_hash.to_string(),
        },
    }
}
