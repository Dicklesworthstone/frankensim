//! Full matching-P1 jump contractions on the producer's retained face mapping.
use super::{BoundFacePair, ConductionError, Cx, ThermalInterfaces};
use super::super::{ContactResistanceGradient, sensitivity::{error, finite, poll}};

impl ThermalInterfaces {
    pub(crate) fn resistance_scale_pullback(
        &self, cx: &Cx<'_>, temperature: &[f64], lambda: &[f64], max_face_pairs: usize,
    ) -> Result<Vec<ContactResistanceGradient>, ConductionError> {
        poll(cx, 0)?;
        let mut count = 0_usize;
        for surface in &self.surfaces {
            poll(cx, count)?;
            count = count.checked_add(surface.faces.len())
                .ok_or_else(|| error("contact face count overflow"))?;
            if count > max_face_pairs { return Err(error("contact controls exceed the paired-face allowance")); }
        }
        if temperature.len() != lambda.len() {
            return Err(error("contact controls require primal and dual fields with equal lengths"));
        }
        let mut result = Vec::new();
        result.try_reserve_exact(self.surfaces.len())
            .map_err(|_| error("contact gradient allocation refused"))?;
        let mut at = 0_usize;
        for surface in &self.surfaces {
            poll(cx, at)?;
            let mut derivative = 0.0;
            for face in &surface.faces {
                if at.is_multiple_of(crate::assemble::ASSEMBLY_TILE) { poll(cx, at)?; }
                derivative = finite(derivative + face.resistance_contraction(temperature, lambda)?)?;
                at += 1;
            }
            result.push(ContactResistanceGradient { interface: surface.name.clone(), derivative,
                card_identity: surface.card_identity, mapped: surface.mapped, face_pairs: surface.faces.len() });
        }
        poll(cx, at)?;
        Ok(result)
    }
}

impl BoundFacePair {
    fn resistance_contraction(&self, temperature: &[f64], lambda: &[f64]) -> Result<f64, ConductionError> {
        let sample = |field: &[f64], vertex: usize| -> Result<f64, ConductionError> {
            finite(*field.get(vertex).ok_or_else(|| error("contact field omits a retained trace vertex"))?)
        };
        let mut jump = [0.0; 3];
        let mut dual_jump = [0.0; 3];
        for i in 0..3 {
            let a = self.side_a_vertices[i];
            let b = self.side_b_vertices[i];
            jump[i] = finite(sample(temperature, a)? - sample(temperature, b)?)?;
            dual_jump[i] = finite(sample(lambda, a)? - sample(lambda, b)?)?;
        }
        let g = finite(1.0 / self.resistance.value_m2_k_per_w())?;
        let mut result = 0.0;
        for (a, &dual) in dual_jump.iter().enumerate() {
            for (b, &primal) in jump.iter().enumerate() {
                // Same conductance and consistent mass as assemble_into.
                let mass = if a == b { self.area_m2 / 6.0 } else { self.area_m2 / 12.0 };
                let contribution = finite(finite(finite(g * mass)? * dual)? * primal)?;
                result = finite(result + contribution)?;
            }
        }
        Ok(result)
    }
}
