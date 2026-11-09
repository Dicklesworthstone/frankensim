//! Fixed regional tangent of the original elastic/Coulomb return map.
//!
//! For a yielding region, e = K (u-z), z = lambda C^-2 e and |C^-1 e| = 1.
//! Differentiation gives D = H - (H b)(H b)^T / (b^T H b), where
//! H = (K^-1 + lambda C^-2)^-1 and b = C^-2 e. The original published
//! plastic increment recovers lambda, without another constitutive update.
//! Elastic regions use D=K. Consumers own their coordinate maps and whether a
//! fixed generalized selection or a classical interior derivative is needed.

use fs_tribo::partial_slip::{PartialSlipLaw, PartialSlipStep};

pub(super) type Block = [[f64; 3]; 3];

/// Row/column-scaled positive effort tangent. The body's friction force has
/// the opposite sign. Scaling roots avoids separately forming dt^2 or ell^2.
pub(super) fn scaled_return_tangent(
    law: &PartialSlipLaw,
    step: &PartialSlipStep,
    row_scale: [f64; 3],
    column_scale: [f64; 3],
) -> Result<Block, &'static str> {
    let p = law.parameters();
    let hard = p.partial_slip_hardening_fraction;
    let fractions = [
        p.partial_slip_onset_fraction * (1.0 - hard),
        1.0 - p.partial_slip_onset_fraction * (1.0 - hard),
    ];
    let mut block = [[0.0; 3]; 3];
    for (region, weight) in [1.0 - hard, hard].into_iter().enumerate() {
        let stiffness = [
            p.tangential_stiffness_n_per_m * weight,
            p.tangential_stiffness_n_per_m * weight,
            p.torsional_stiffness_nm_per_rad * weight,
        ];
        let plastic = step.region_plastic_increment[region];
        let mut diagonal = stiffness;
        let mut normal = [0.0; 3];
        if plastic.iter().any(|&v| v != 0.0) {
            let capacity = [
                step.capacity.kinetic_force_n,
                step.capacity.kinetic_force_n,
                step.capacity.kinetic_torque_nm,
            ]
            .map(|v| v * fractions[region]);
            let elastic = [
                step.next_state.region_elastic_displacement_m()[region][0],
                step.next_state.region_elastic_displacement_m()[region][1],
                step.next_state.region_elastic_twist_rad()[region],
            ];
            let mut q = [0.0; 3];
            let mut numerator = 0.0;
            let mut denominator = 0.0;
            for j in 0..3 {
                q[j] = stiffness[j] * elastic[j] / capacity[j];
                numerator += q[j] * (plastic[j] * capacity[j]);
                denominator += q[j] * q[j];
            }
            let lambda = numerator / denominator;
            if !lambda.is_finite() || lambda < 0.0 {
                return Err("unrepresentable regional return multiplier");
            }
            for j in 0..3 {
                diagonal[j] =
                    stiffness[j] / (1.0 + (lambda / capacity[j]) * (stiffness[j] / capacity[j]));
                normal[j] = diagonal[j].sqrt() * q[j] / capacity[j];
            }
            let norm = normal[0].hypot(normal[1]).hypot(normal[2]);
            if !norm.is_finite() || norm <= 0.0 {
                return Err("unrepresentable regional return tangent");
            }
            normal = normal.map(|v| v / norm);
        }
        let row_root: [f64; 3] = std::array::from_fn(|j| row_scale[j] * diagonal[j].sqrt());
        let column_root: [f64; 3] =
            std::array::from_fn(|j| column_scale[j] * diagonal[j].sqrt());
        for i in 0..3 {
            for j in 0..3 {
                let identity = if i == j { 1.0 } else { 0.0 };
                block[i][j] += row_root[i] * column_root[j] * (identity - normal[i] * normal[j]);
            }
        }
    }
    if block.iter().flatten().any(|v| !v.is_finite()) {
        return Err("unrepresentable contact return tangent");
    }
    Ok(block)
}
