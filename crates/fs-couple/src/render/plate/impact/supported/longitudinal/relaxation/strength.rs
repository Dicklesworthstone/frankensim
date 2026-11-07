//! Source layer identity, actual velocity and pHS histories in tensile admission.
use super::*;

#[derive(Clone, Debug)]
struct LayerArm {
    first_memory: usize,
    delta_modulus_pa: f64,
    normalization: [f64; 2],
}

#[derive(Clone, Debug)]
pub(in crate::render::plate::impact::supported) struct LayerTensile {
    pub(in crate::render::plate::impact::supported) guard: Arc<StringTensileGuard>,
    layers: Vec<Vec<LayerArm>>,
    viscosity_pa_s: Vec<f64>,
}

impl SupportedStrings {
    pub(in crate::render::plate::impact::supported) fn bind_viscous_tensile(mut self)
        -> Result<Self, ImpactError>
    {
        let Some(sources) = &self.0.axial_viscosity_sources else {
            return Err(invalid("viscous tensile binding needs its original layer sources"));
        };
        let mut bindings = self.0.axial_tensile.clone();
        for (wire, source) in sources.iter().enumerate() {
            let SupportedPotential::Longitudinal(p) = &self.0.potentials[wire] else { continue; };
            if source.viscous_rigidity_n_s() == 0.0 { continue; }
            if let Some(guard) = &p.tensile {
                bindings[wire] = Some(LayerTensile::new(Arc::clone(guard), Some(source))?);
            }
        }
        let data = Rc::make_mut(&mut self.0);
        for (p, binding) in data.potentials.iter_mut().zip(&bindings) {
            if binding.is_some() {
                let SupportedPotential::Longitudinal(p) = p else { unreachable!() };
                p.tensile = None;
            }
        }
        data.axial_tensile = bindings;
        Ok(self)
    }

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
                    return Ok(None);
                };
                if !arms.iter().any(|arm| arm.string == wire) {
                    return Ok(self.0.axial_tensile[wire].clone());
                }
                let guard = p.tensile.as_ref().or_else(||
                    self.0.axial_tensile[wire].as_ref().map(|binding| &binding.guard));
                guard.map(|guard| LayerTensile::new(
                    Arc::clone(guard),
                    self.0.axial_viscosity_sources.as_ref().map(|sources| &sources[wire]),
                )).transpose()
            })
            .collect::<Result<_, ImpactError>>()?;
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
        // Raw aggregate arms cannot bypass an already selected viscous layer
        // guard. Only this source path may temporarily remove that metadata
        // while binding the original Maxwell mechanics, then restore it with
        // the exact layer-to-memory projection above.
        Rc::make_mut(&mut self.0).axial_tensile.fill(None);
        let mut result = self.with_longitudinal_relaxation(arms, initial, maximum_states)?;
        Rc::make_mut(&mut result.0).axial_tensile = tensile;
        Ok(result)
    }

    pub(in crate::render::plate::impact) fn has_state_tensile(
        &self,
        wire: usize,
    ) -> bool {
        self.0.axial_tensile.get(wire).is_some_and(Option::is_some)
            || (self.0.axial_viscosity_n_s.as_ref().is_some_and(|loss|
                loss.get(wire).is_some_and(|rigidity| *rigidity > 0.0))
                && matches!(self.0.potentials.get(wire),
                    Some(SupportedPotential::Longitudinal(p)) if p.tensile.is_some()))
    }
}

impl LayerTensile {
    fn new(
        guard: Arc<StringTensileGuard>,
        viscosity: Option<&fs_material::section::axial::ConcentricAxialViscosity>,
    ) -> Result<Self, ImpactError> {
        let count = guard.capacity.section().layers().len();
        let mut viscosity_pa_s = vec![0.0; count];
        if let Some(source) = viscosity {
            if source.elastic() != guard.capacity.section() {
                return Err(invalid("viscous tensile stress needs the same original section and material query"));
            }
            for (eta, layer) in viscosity_pa_s.iter_mut().zip(source.elastic().layers()) {
                *eta = layer.material.property(source.viscosity_property())
                    .ok_or_else(|| invalid("tensile layer lost its selected axial viscosity"))?
                    .value_si();
            }
        }
        Ok(Self { guard, layers: vec![Vec::new(); count], viscosity_pa_s })
    }

    pub(super) fn assess(
        &self,
        p: &LongitudinalPotential,
        local: &[f64],
        strain: &[f64],
        strain_rate: &[f64],
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
                    self.viscosity_pa_s[layer],
                    strain_rate,
                    &self.layers[layer],
                    memory,
                )?;
                let bending = young * bending_strain;
                if bending_strain != 0.0 && bending == 0.0 {
                    return Err(invalid(
                        "supported tensile bending stress is unrepresentable",
                    ));
                }
                Ok(axial + bending)
            })
    }
}

// Sum SIGNED equilibrium/arm/dashpot coefficients before the continuous-span
// absolute bound. All three can reinforce or cancel, and a layer's mean stress
// can be compressive even while the complete string remains tensioned.
fn layer_stress_bound(
    young: f64,
    prestress_strain: f64,
    strain: &[f64],
    viscosity_pa_s: f64,
    strain_rate: &[f64],
    arms: &[LayerArm],
    memory: &[f64],
) -> Result<f64, ImpactError> {
    if strain_rate.len() != strain.len() {
        return Err(invalid("tensile layer needs every original strain-rate channel"));
    }
    let mut bound = 0.0;
    for (channel, &epsilon) in strain.iter().enumerate() {
        let elastic_strain = epsilon + if channel == 0 { prestress_strain } else { 0.0 };
        let mut coefficient = young * elastic_strain;
        if elastic_strain != 0.0 && coefficient == 0.0 {
            return Err(invalid(
                "supported tensile elastic stress is unrepresentable",
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
        let viscous = viscosity_pa_s * strain_rate[channel];
        if !viscous.is_finite()
            || (viscosity_pa_s > 0.0 && strain_rate[channel] != 0.0 && viscous == 0.0)
        {
            return Err(invalid("tensile layer viscous stress is unrepresentable"));
        }
        coefficient += viscous;
        if !coefficient.is_finite() {
            return Err(invalid("supported tensile layer coefficient overflows"));
        }
        bound += if channel == 0 {
            coefficient
        } else {
            coefficient.abs()
        };
    }
    if !bound.is_finite() {
        return Err(invalid("supported tensile layer envelope overflows"));
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
        let stress = layer_stress_bound(4.0, 2.0, &[1.0, 2.0], 0.0, &[0.0; 2], &arms,
            &[0.0, 3.0, 4.0 * root_two, 8.0]).unwrap();
        assert_eq!(stress, 11.0);
        let reinforced = layer_stress_bound(4.0, 2.0, &[1.0, 2.0], 0.0, &[0.0; 2], &arms, &[0.0; 4]).unwrap();
        assert_eq!(reinforced, 35.0);
        // Relaxed memories remove every dynamic arm without removing preload.
        let relaxed = layer_stress_bound(4.0, 2.0, &[1.0, 2.0], 0.0, &[0.0; 2], &arms,
            &[root_two, 2.0, 2.0 * root_two, 4.0]).unwrap();
        assert_eq!(relaxed, 20.0);
        assert!(layer_stress_bound(4.0, 2.0, &[1.0, 2.0], 0.0, &[0.0; 2], &arms, &[0.0; 3]).is_err());
        // The same unrelaxed elastic + Maxwell harmonic is 18 Pa. A signed
        // dashpot can cancel it exactly or reinforce it; |sigma_e|+|sigma_v|
        // would incorrectly reject the cancellation. Its mean also stays signed.
        let cancelled = layer_stress_bound(4.0, 2.0, &[1.0, 2.0], 3.0,
            &[-1.0, -6.0], &arms, &[0.0; 4]).unwrap();
        assert_eq!(cancelled, 14.0);
        let viscous_reinforced = layer_stress_bound(4.0, 2.0, &[1.0, 2.0], 3.0,
            &[1.0, 6.0], &arms, &[0.0; 4]).unwrap();
        assert_eq!(viscous_reinforced, 56.0);
    }
}
