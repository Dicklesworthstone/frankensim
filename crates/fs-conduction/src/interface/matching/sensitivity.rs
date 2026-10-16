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

    /// dK_contact/dT * T, using the SAME matched trace map and consistent face
    /// mass as the primal. The caller supplies admitted source dR/dT_interface;
    /// T_interface is the mean of both P1 face means (six weights of 1/6).
    /// Absent pairs have fixed resistance. No conductivity or storage is added.
    pub(crate) fn temperature_feedback(
        &self, cx: &Cx<'_>, temperature: &[f64],
        slopes: &std::collections::BTreeMap<super::InterfaceFacePair, f64>,
        max_face_pairs: usize,
    ) -> Result<fs_sparse::Csr, ConductionError> {
        poll(cx, 0)?;
        let count = self.surfaces.iter().try_fold(0usize, |n, surface|
            n.checked_add(surface.faces.len())).ok_or_else(|| error("contact face count overflow"))?;
        if count > max_face_pairs || count.checked_mul(36 * 24)
            .is_none_or(|bytes| bytes > isize::MAX as usize)
        {
            return Err(error("contact temperature feedback exceeds its paired-face allowance"));
        }
        let mut coo = fs_sparse::Coo::new(temperature.len(), temperature.len());
        let mut found = 0usize;
        for surface in &self.surfaces {
            for face in &surface.faces {
                poll(cx, found)?;
                let Some(&slope) = slopes.get(&face.pair) else { continue; };
                found += 1;
                let g = finite(1.0 / face.resistance.value_m2_k_per_w())?;
                let factor = finite(-finite(slope)? * g)?;
                let mut jump = [0.0; 3];
                for (i, value) in jump.iter_mut().enumerate() {
                    let sample = |v| temperature.get(v).copied()
                        .ok_or_else(|| error("temperature field omits a contact trace vertex"));
                    *value = finite(sample(face.side_a_vertices[i])? - sample(face.side_b_vertices[i])?)?;
                }
                for i in 0..3 {
                    let mut load = 0.0;
                    for (j, &jump) in jump.iter().enumerate() {
                        let mass = if i == j { face.area_m2 / 6.0 } else { face.area_m2 / 12.0 };
                        load = finite(load + finite(finite(g * mass)? * jump)?)?;
                    }
                    let entry = finite(finite(factor * load)? / 6.0)?;
                    for &column in face.side_a_vertices.iter().chain(&face.side_b_vertices) {
                        coo.push(face.side_a_vertices[i], column, entry);
                        coo.push(face.side_b_vertices[i], column, -entry);
                    }
                }
            }
        }
        if found != slopes.len() { return Err(error("contact slope names a pair outside the bound operator")); }
        poll(cx, found)?;
        Ok(coo.assemble())
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
