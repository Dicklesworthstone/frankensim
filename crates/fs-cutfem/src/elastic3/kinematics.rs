//! Cell-constant volumetric projection on the retained physical integration rule.
//!
//! For isotropy, W = mu*dev(eps):dev(eps) + K*mean(div u)^2/2.
//! This is the condensed constant-pressure/mean-dilatation B-bar formulation;
//! see Hughes' projection as described in MOOSE's Volumetric Locking Correction.
//! Deviatoric quadrature is unchanged. No pressure unknown, quadrature rebuild,
//! per-apply material query or additional global solve is introduced.

use super::*;

#[derive(Clone, Copy)]
pub(super) enum CellKinematics3 {
    Displacement,
    BBar {
        mean_gradients: [[f64; 3]; 8],
        bulk: f64,
        shear: f64,
    },
}

impl CellKinematics3 {
    fn isotropic_moduli(material: &ReferenceMaterial3) -> Result<(f64, f64), ElasticityError3> {
        let (bulk, shear) = material.isotropic_bulk_shear().ok_or(ElasticityError3::Invalid(
            "B-bar requires an explicitly isotropic material",
        ))?;
        if !bulk.is_finite() || !shear.is_finite() || bulk <= 0.0 || shear <= 0.0 {
            return Err(ElasticityError3::Invalid(
                "B-bar requires representable positive bulk and shear moduli",
            ));
        }
        Ok((bulk, shear))
    }

    /// Retain the original physical-volume projection; only material changes.
    pub(super) fn with_material(self, material: &ReferenceMaterial3) -> Result<Self, ElasticityError3> {
        match self {
            Self::Displacement => Ok(Self::Displacement),
            Self::BBar { mean_gradients, .. } => {
                let (bulk, shear) = Self::isotropic_moduli(material)?;
                Ok(Self::BBar { mean_gradients, bulk, shear })
            }
        }
    }

    pub(super) fn prepare(
        formulation: ElasticityFormulation3,
        bounds: HexCell,
        rules: &CutRules3,
        material: &ReferenceMaterial3,
        control: &mut QuadratureControl3<'_>,
    ) -> Result<Self, ElasticityError3> {
        if formulation == ElasticityFormulation3::Displacement {
            return Ok(Self::Displacement);
        }
        let (bulk, shear) = Self::isotropic_moduli(material)?;
        let volume = rules.volume();
        if !volume.is_finite() || volume <= 0.0 {
            return Err(ElasticityError3::Invalid(
                "B-bar physical volume is invalid",
            ));
        }
        let mut mean_gradients = [[0.0; 3]; 8];
        for &(point, weight) in rules.bulk() {
            control.poll()?;
            let (_, gradients) = q1(bounds, point);
            let fraction = weight / volume;
            for (mean, gradient) in mean_gradients.iter_mut().zip(gradients) {
                for axis in 0..3 {
                    mean[axis] += fraction * gradient[axis];
                }
            }
        }
        // Retain the partition-of-unity null mode explicitly; volumetric
        // stiffness must not amplify a rounded gradient sum into rigid stress.
        for axis in 0..3 {
            mean_gradients[0][axis] = -(1..8).map(|a| mean_gradients[a][axis]).sum::<f64>();
        }
        if mean_gradients.iter().flatten().any(|v| !v.is_finite()) {
            return Err(ElasticityError3::Invalid("B-bar mean gradient overflow"));
        }
        control.poll()?;
        Ok(Self::BBar {
            mean_gradients,
            bulk,
            shear,
        })
    }

    pub(super) const fn is_projected(self) -> bool {
        matches!(self, Self::BBar { .. })
    }

    pub(super) const fn projected_moduli(self) -> Option<(f64, f64, [[f64; 3]; 8])> {
        match self {
            Self::Displacement => None,
            Self::BBar {
                mean_gradients,
                bulk,
                shear,
            } => Some((bulk, shear, mean_gradients)),
        }
    }
    /// Work-conjugate constitutive basis used by both recovery and surface terms.
    pub(super) fn basis_stress(
        self,
        material: &ReferenceMaterial3,
        gradient: [f64; 3],
        corner: usize,
        component: usize,
    ) -> [f64; 6] {
        let CellKinematics3::BBar {
            mean_gradients,
            bulk,
            shear,
        } = self
        else {
            return material.basis_stress(gradient, component);
        };
        let strain = engineering_strain(gradient, component);
        let pressure = bulk * mean_gradients[corner][component];
        let third_trace = gradient[component] / 3.0;
        std::array::from_fn(|i| {
            if i < 3 {
                (2.0 * shear) * (strain[i] - third_trace) + pressure
            } else {
                shear * strain[i]
            }
        })
    }

    pub(super) fn basis_traction(
        self,
        material: &ReferenceMaterial3,
        gradient: [f64; 3],
        corner: usize,
        component: usize,
        normal: [f64; 3],
    ) -> [f64; 3] {
        if !self.is_projected() {
            return material.traction(gradient, component, normal);
        }
        let stress = self.basis_stress(material, gradient, corner, component);
        [
            stress[0] * normal[0] + stress[3] * normal[1] + stress[5] * normal[2],
            stress[3] * normal[0] + stress[1] * normal[1] + stress[4] * normal[2],
            stress[5] * normal[0] + stress[4] * normal[1] + stress[2] * normal[2],
        ]
    }

    // Stabilization must not reintroduce the projected-out bulk constraint.
    pub(super) fn ghost_modulus(self, material: &ReferenceMaterial3) -> f64 {
        match self {
            Self::Displacement => material.ghost_modulus(),
            Self::BBar { shear, .. } => shear,
        }
    }
}

impl Cell3 {
    pub(super) fn basis_stress(&self, gradient: [f64; 3], corner: usize, component: usize) -> [f64; 6] {
        self.kinematics.basis_stress(&self.material, gradient, corner, component)
    }

    pub(super) fn basis_traction(&self, gradient: [f64; 3], corner: usize,
        component: usize, normal: [f64; 3]) -> [f64; 3] {
        self.kinematics.basis_traction(&self.material, gradient, corner, component, normal)
    }

    pub(super) fn projected_mandel_strain(
        &self,
        gradients: &[[f64; 3]; 8],
        local: &[f64; 24],
    ) -> [f64; 6] {
        self.project_mandel_strain(super::stress::mandel_strain(gradients, local), local)
    }

    pub(super) fn project_mandel_strain(
        &self,
        mut strain: [f64; 6],
        local: &[f64; 24],
    ) -> [f64; 6] {
        if let CellKinematics3::BBar { mean_gradients, .. } = self.kinematics {
            let mut mean_divergence = 0.0;
            for a in 1..8 {
                for c in 0..3 {
                    mean_divergence += (local[3 * a + c] - local[c]) * mean_gradients[a][c];
                }
            }
            let correction = (mean_divergence - strain[..3].iter().sum::<f64>()) / 3.0;
            for value in &mut strain[..3] {
                *value += correction;
            }
        }
        strain
    }

    // Stabilization must not reintroduce the projected-out bulk constraint.
    // A shear-scale derivative jump controls the vector extension. Adequacy for
    // an arbitrary cut/support configuration still requires a convergence study.
    pub(super) fn ghost_modulus(&self) -> f64 {
        self.kinematics.ghost_modulus(&self.material)
    }
}
