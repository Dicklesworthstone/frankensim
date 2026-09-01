//! ONNX export metadata: the browser needs a JSON sidecar alongside the
//! model file so the inference runtime knows the observation layout, action
//! scaling, and joint limits — the same honesty pattern as the WASM kernel's
//! schema envelope.

pub struct PolicyMetadata {
    pub model_name: String,
    pub architecture: String,
    pub param_count: usize,
    pub obs_layout: ObsLayout,
    pub action_scaling: ActionScaling,
    pub training_provenance: TrainingProvenance,
}

pub struct ObsLayout {
    pub total_dims: usize,
    pub joint_positions: (usize, usize),
    pub joint_velocities: (usize, usize),
    pub base_quaternion: (usize, usize),
    pub phase_sin_cos: (usize, usize),
    pub terrain_heights: Option<(usize, usize)>,
}

pub struct ActionScaling {
    pub action_dim: usize,
    pub lower_limits: Vec<f32>,
    pub upper_limits: Vec<f32>,
    pub tanh_scaled: bool,
}

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

impl PolicyMetadata {
    /// Emit deterministic JSON for the ONNX sidecar without adding a runtime
    /// serialization dependency. Non-finite numeric metadata is refused so
    /// this always returns valid JSON on success.
    pub fn to_json(&self) -> Result<String, &'static str> {
        use std::fmt::Write as _;

        let mut json = String::from("{\"model_name\":");
        push_json_string(&mut json, &self.model_name);
        json.push_str(",\"architecture\":");
        push_json_string(&mut json, &self.architecture);
        let _ = write!(json, ",\"param_count\":{}", self.param_count);
        let _ = write!(
            json,
            ",\"obs_layout\":{{\"total_dims\":{},\"joint_positions\":[{},{}],\
             \"joint_velocities\":[{},{}],\"base_quaternion\":[{},{}],\"phase_sin_cos\":[{},{}],\
             \"terrain_heights\":",
            self.obs_layout.total_dims,
            self.obs_layout.joint_positions.0,
            self.obs_layout.joint_positions.1,
            self.obs_layout.joint_velocities.0,
            self.obs_layout.joint_velocities.1,
            self.obs_layout.base_quaternion.0,
            self.obs_layout.base_quaternion.1,
            self.obs_layout.phase_sin_cos.0,
            self.obs_layout.phase_sin_cos.1,
        );
        match self.obs_layout.terrain_heights {
            Some((start, end)) => {
                let _ = write!(json, "[{start},{end}]");
            }
            None => json.push_str("null"),
        }
        let _ = write!(
            json,
            "}},\"action_scaling\":{{\"action_dim\":{}",
            self.action_scaling.action_dim
        );
        json.push_str(",\"lower_limits\":");
        push_json_f32_array(&mut json, &self.action_scaling.lower_limits, "lower_limits")?;
        json.push_str(",\"upper_limits\":");
        push_json_f32_array(&mut json, &self.action_scaling.upper_limits, "upper_limits")?;
        let _ = write!(
            json,
            ",\"tanh_scaled\":{}}}",
            self.action_scaling.tanh_scaled
        );
        json.push_str(",\"training_provenance\":{\"algorithm\":");
        push_json_string(&mut json, &self.training_provenance.algorithm);
        json.push_str(",\"optimizer\":");
        push_json_string(&mut json, &self.training_provenance.optimizer);
        json.push_str(",\"outer_hpo\":");
        push_json_string(&mut json, &self.training_provenance.outer_hpo);
        let _ = write!(
            json,
            ",\"total_training_steps\":{}",
            self.training_provenance.total_training_steps
        );
        json.push_str(",\"final_mean_reward\":");
        push_json_f32(
            &mut json,
            self.training_provenance.final_mean_reward,
            "final_mean_reward",
        )?;
        let _ = write!(
            json,
            ",\"environments\":{},\"seed\":{}",
            self.training_provenance.environments, self.training_provenance.seed
        );
        json.push_str(",\"checkpoint_hash\":");
        push_json_string(&mut json, &self.training_provenance.checkpoint_hash);
        json.push_str("}}");
        Ok(json)
    }
}

fn push_json_string(json: &mut String, value: &str) {
    use std::fmt::Write as _;

    json.push('\"');
    for ch in value.chars() {
        match ch {
            '\"' => json.push_str("\\\""),
            '\\' => json.push_str("\\\\"),
            '\u{08}' => json.push_str("\\b"),
            '\u{0C}' => json.push_str("\\f"),
            '\n' => json.push_str("\\n"),
            '\r' => json.push_str("\\r"),
            '\t' => json.push_str("\\t"),
            ch if ch <= '\u{1F}' => {
                let _ = write!(json, "\\u{:04x}", ch as u32);
            }
            ch => json.push(ch),
        }
    }
    json.push('\"');
}

fn push_json_f32_array(
    json: &mut String,
    values: &[f32],
    field: &'static str,
) -> Result<(), &'static str> {
    json.push('[');
    for (index, value) in values.iter().enumerate() {
        if index != 0 {
            json.push(',');
        }
        push_json_f32(json, *value, field)?;
    }
    json.push(']');
    Ok(())
}

fn push_json_f32(json: &mut String, value: f32, field: &'static str) -> Result<(), &'static str> {
    use std::fmt::Write as _;

    if !value.is_finite() {
        return Err(field);
    }
    let _ = write!(json, "{value}");
    Ok(())
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

#[cfg(test)]
mod tests {
    use super::generate_metadata;

    #[test]
    fn metadata_json_escapes_strings_and_preserves_named_fields() {
        let metadata = generate_metadata(
            "g1\"sidecar\\newline\n",
            17,
            42,
            1.25,
            "checkpoint\t\r\u{7}",
        );

        let json = metadata.to_json().expect("finite metadata emits JSON");
        assert_eq!(json, metadata.to_json().expect("deterministic JSON"));
        assert!(json.contains(r#""model_name":"g1\"sidecar\\newline\n""#));
        assert!(json.contains(r#""checkpoint_hash":"checkpoint\t\r\u0007""#));
        assert!(json.contains(r#""param_count":17"#));
        assert!(json.contains(r#""joint_positions":[0,15]"#));
        assert!(json.contains(r#""terrain_heights":[36,42]"#));
        assert!(json.contains(r#""tanh_scaled":true"#));
        assert!(json.contains(r#""total_training_steps":42"#));
    }
}
