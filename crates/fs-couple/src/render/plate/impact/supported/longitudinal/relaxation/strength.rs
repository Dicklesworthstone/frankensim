//! Source layer identity and actual pHS histories in tensile admission.
use super::*;

#[derive(Clone, Debug)]
struct LayerArm {
    first_memory: usize,
    delta_modulus_pa: f64,
    normalization: [f64; 2],
}

#[derive(Clone, Debug)]
pub(super) struct MemoryTensile {
    pub(super) guard: Arc<StringTensileGuard>,
    layers: Vec<Vec<LayerArm>>,
}

impl SupportedStrings {
    // Only the material adapter knows which physical layer supplied an arm.
    // Public aggregate-EA arms retain their explicit refusal for selected
    // strengths; guessing equal or equilibrium load sharing is not admissible.
    pub(crate) fn with_longitudinal_layer_relaxation(
        mut self,
        arms: Vec<SupportedRelaxationArm>,
        layers: &[(usize, f64)],
        initial: InitialMemory,
        maximum_states: usize,
    ) -> Result<Self, ImpactError> {
        if layers.len() != arms.len() {
            return Err(invalid(
                "every hereditary tensile arm needs its original layer and modulus",
            ));
        }
        let mut tensile: Vec<_> = self
            .0
            .potentials
            .iter()
            .enumerate()
            .map(|(wire, p)| {
                let SupportedPotential::Longitudinal(p) = p else {
                    return None;
                };
                if !arms.iter().any(|arm| arm.string == wire) {
                    return None;
                }
                p.tensile.as_ref().map(|guard| MemoryTensile {
                    guard: Arc::clone(guard),
                    layers: vec![Vec::new(); guard.capacity.section().layers().len()],
                })
            })
            .collect();
        let mut first_memory = 0;
        for (arm, &(layer, delta)) in arms.iter().zip(layers) {
            let Some(SupportedPotential::Longitudinal(p)) = self.0.potentials.get(arm.string)
            else {
                return Err(invalid(
                    "hereditary tensile arm lost its longitudinal string",
                ));
            };
            if !delta.is_finite() || delta <= 0.0 {
                return Err(invalid(
                    "hereditary tensile arm needs its positive source modulus",
                ));
            }
            if let Some(binding) = &mut tensile[arm.string] {
                let geometry = binding
                    .guard
                    .capacity
                    .section()
                    .layer_integrals()
                    .get(layer)
                    .ok_or_else(|| {
                        invalid("hereditary tensile arm names an absent source layer")
                    })?;
                let ea = delta * geometry.area_m2();
                if !ea.is_finite()
                    || ea <= 0.0
                    || !arm.axial_rigidity_n.is_finite()
                    || (ea - arm.axial_rigidity_n).abs() > 32.0 * f64::EPSILON * ea
                {
                    return Err(invalid(
                        "hereditary tensile arm must retain its original layer area and modulus",
                    ));
                }
                let stiffness = arm.axial_rigidity_n * self.0.wires[arm.string].length_m();
                let normalization = [det::sqrt(stiffness), det::sqrt(0.5 * stiffness)];
                if normalization.iter().any(|s| !s.is_finite() || *s <= 0.0) {
                    return Err(invalid(
                        "hereditary tensile memory normalization is unrepresentable",
                    ));
                }
                binding.layers[layer].push(LayerArm {
                    first_memory,
                    delta_modulus_pa: delta,
                    normalization,
                });
            }
            first_memory += 1 + p.model.axial_modes().len();
        }
        // The displacement-only body gate still checks its original geometric
        // domain. The one actual memory gate below now owns total tensile
        // admission, so a relieving arm cannot be rejected as elastic overload.
        for (p, binding) in Rc::make_mut(&mut self.0)
            .potentials
            .iter_mut()
            .zip(&tensile)
        {
            if binding.is_some() {
                let SupportedPotential::Longitudinal(p) = p else {
                    unreachable!()
                };
                p.tensile = None;
            }
        }
        let mut result = self.with_longitudinal_relaxation(arms, initial, maximum_states)?;
        Rc::make_mut(&mut result.0)
            .axial_relaxation
            .as_mut()
            .ok_or_else(|| invalid("hereditary tensile binding lost its material owner"))?
            .tensile = tensile;
        Ok(result)
    }

    pub(in crate::render::plate::impact::supported) fn has_memory_tensile(
        &self,
        wire: usize,
    ) -> bool {
        self.0
            .axial_relaxation
            .as_ref()
            .and_then(|memory| memory.tensile.get(wire))
            .is_some_and(Option::is_some)
    }
}

impl MemoryTensile {
    pub(super) fn assess(
        &self,
        p: &LongitudinalPotential,
        local: &[f64],
        strain: &[f64],
        maximum_tension_n: f64,
        memory: &[f64],
    ) -> Result<StringTensileAssessment, ImpactError> {
        self.guard
            .assess_normal_stress(local, 2, maximum_tension_n, |layer, bending_strain| {
                let young =
                    self.guard.capacity.section().layer_integrals()[layer].young_modulus_pa();
                let axial = layer_stress_bound(
                    young,
                    p.tension / p.ea,
                    strain,
                    &self.layers[layer],
                    memory,
                )?;
                let bending = young * bending_strain;
                if bending_strain != 0.0 && bending == 0.0 {
                    return Err(invalid(
                        "hereditary tensile bending stress is unrepresentable",
                    ));
                }
                Ok(axial + bending)
            })
    }
}

// Sum SIGNED equilibrium/arm coefficients before the continuous-span absolute
// bound. Arm stresses can reinforce or cancel one another, and a layer's mean
// stress can be compressive even while the complete string remains tensioned.
fn layer_stress_bound(
    young: f64,
    prestress_strain: f64,
    strain: &[f64],
    arms: &[LayerArm],
    memory: &[f64],
) -> Result<f64, ImpactError> {
    let mut bound = 0.0;
    for (channel, &epsilon) in strain.iter().enumerate() {
        let elastic_strain = epsilon + if channel == 0 { prestress_strain } else { 0.0 };
        let mut coefficient = young * elastic_strain;
        if elastic_strain != 0.0 && coefficient == 0.0 {
            return Err(invalid(
                "hereditary tensile elastic stress is unrepresentable",
            ));
        }
        for arm in arms {
            let history = memory
                .get(arm.first_memory + channel)
                .ok_or_else(|| invalid("hereditary tensile layer history is incomplete"))?
                / arm.normalization[usize::from(channel != 0)];
            let residual = epsilon - history;
            let stress = arm.delta_modulus_pa * residual;
            if !history.is_finite() || !stress.is_finite() || (residual != 0.0 && stress == 0.0) {
                return Err(invalid("hereditary tensile arm stress is unrepresentable"));
            }
            coefficient += stress;
        }
        if !coefficient.is_finite() {
            return Err(invalid("hereditary tensile layer coefficient overflows"));
        }
        bound += if channel == 0 {
            coefficient
        } else {
            coefficient.abs()
        };
    }
    if !bound.is_finite() {
        return Err(invalid("hereditary tensile layer envelope overflows"));
    }
    Ok(bound)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn signed_layer_histories_cancel_before_the_spatial_bound() {
        let root_two = det::sqrt(2.0);
        let arms = [
            LayerArm {
                first_memory: 0,
                delta_modulus_pa: 2.0,
                normalization: [root_two, 1.0],
            },
            LayerArm {
                first_memory: 2,
                delta_modulus_pa: 3.0,
                normalization: [2.0 * root_two, 2.0],
            },
        ];
        // E*epsilon_1 = 8, arm 1 = -2 and arm 2 = -6: exact cancellation.
        // The independent mean is 4*(prestrain 2 + epsilon 1) + 2 - 3 = 11.
        let stress = layer_stress_bound(4.0, 2.0, &[1.0, 2.0], &arms,
            &[0.0, 3.0, 4.0 * root_two, 8.0]).unwrap();
        assert_eq!(stress, 11.0);
        let reinforced = layer_stress_bound(4.0, 2.0, &[1.0, 2.0], &arms, &[0.0; 4]).unwrap();
        assert_eq!(reinforced, 35.0);
        // Relaxed memories remove every dynamic arm without removing preload.
        let relaxed = layer_stress_bound(4.0, 2.0, &[1.0, 2.0], &arms,
            &[root_two, 2.0, 2.0 * root_two, 4.0]).unwrap();
        assert_eq!(relaxed, 20.0);
        assert!(layer_stress_bound(4.0, 2.0, &[1.0, 2.0], &arms, &[0.0; 3]).is_err());
    }
}
