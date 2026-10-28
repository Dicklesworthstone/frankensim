//! Move a native trajectory across caller-owned material-state lifetimes.
//!
//! This is an in-process continuation, not a serialization or material rebind.
//! It retains the exact assembled operators, numerical controls, forces, thermal
//! field, motion, accepted step count and support state. Rebinding only admits
//! the identical mesh, constraints and material values/source identities.
use super::*;
use crate::linear3::{TetElasticMaterial, TetMaterialField};

#[derive(Clone)]
enum Materials {
    Uniform(TetElasticMaterial),
    PerElement(Vec<TetElasticMaterial>),
}

/// Owned dynamic history with shared matrices and retained original geometry.
/// No public mutable state can bypass the native acceptance gates.
#[derive(Clone)]
pub struct TetTransientCheckpoint3<'g, 'm> {
    geometry: &'g TetLoadGeometry3<'m>,
    materials: Materials,
    fixed_dofs: Vec<usize>,
    operators: Arc<TetElasticAssembly>,
    state: DynamicState,
    thermal: Vec<f64>,
    thermal_strains: Option<response::FrozenThermal3>,
    force: Vec<f64>,
    full_force: Vec<(usize, f64)>,
    config: TransientConfig3,
    prescribed: Option<motion::PrescribedMotion3>,
}

impl<'g, 'm> TetTransientRun3<'g, 'm, '_, '_> {
    /// Detach the accepted motion from the borrowed problem without solving or
    /// reassembling anything. Consuming the run moves its state buffers; only
    /// the material binding and support indices need new owned storage.
    pub fn into_checkpoint(
        self,
        cx: &Cx<'_>,
    ) -> Result<TetTransientCheckpoint3<'g, 'm>, TransientError3> {
        check(cx)?;
        let materials = match self.problem.materials {
            TetMaterialField::Uniform(material) => Materials::Uniform(material.clone()),
            TetMaterialField::PerElement(values) => {
                let mut retained = Vec::with_capacity(values.len());
                for material in values {
                    check(cx)?;
                    retained.push(material.clone());
                }
                Materials::PerElement(retained)
            }
        };
        let fixed_dofs = self.problem.fixed_dofs.to_vec();
        check(cx)?;
        Ok(TetTransientCheckpoint3 {
            geometry: self.geometry,
            materials,
            fixed_dofs,
            operators: self.operators,
            state: self.state,
            thermal: self.thermal,
            thermal_strains: self.thermal_strains,
            force: self.force,
            full_force: self.full_force,
            config: self.config,
            prescribed: self.prescribed,
        })
    }
}

impl<'g, 'm> TetTransientCheckpoint3<'g, 'm> {
    pub fn state(&self) -> &DynamicState {
        &self.state
    }

    /// Bind the EXACT original material problem to an immutable fork of the
    /// retained history. Neither physical time nor the native step budget resets.
    /// Changed coefficients, source identities, density, mesh or supports require
    /// the original material/thermal update APIs and are refused here.
    pub fn bind<'p, 's>(
        &self,
        problem: &'p TetLinearElasticProblem<'s>,
        cx: &Cx<'_>,
    ) -> Result<TetTransientRun3<'g, 'm, 'p, 's>, TransientError3> {
        check(cx)?;
        if problem.nodes_m != self.geometry.nodes
            || problem.tetrahedra != self.geometry.tets
            || problem.fixed_dofs != self.fixed_dofs
        {
            return Err(TransientError3::Invalid(
                "checkpoint requires its exact original mesh and supports",
            ));
        }
        let count = problem.tetrahedra.len();
        let material = |field: TetMaterialField<'s>, i| match field {
            TetMaterialField::Uniform(value) => Some(value),
            TetMaterialField::PerElement(values) => values.get(i),
        };
        if let TetMaterialField::PerElement(values) = problem.materials {
            if values.len() != count {
                return Err(TransientError3::Invalid(
                    "checkpoint requires one material per tetrahedron",
                ));
            }
        }
        for i in 0..count {
            check(cx)?;
            let expected = match &self.materials {
                Materials::Uniform(value) => Some(value),
                Materials::PerElement(values) => values.get(i),
            };
            if expected != material(problem.materials, i) {
                return Err(TransientError3::Invalid(
                    "checkpoint material or source identity differs",
                ));
            }
        }
        check(cx)?;
        Ok(TetTransientRun3 {
            geometry: self.geometry,
            problem,
            operators: self.operators.clone(),
            state: self.state.clone(),
            thermal: self.thermal.clone(),
            thermal_strains: self.thermal_strains.clone(),
            force: self.force.clone(),
            full_force: self.full_force.clone(),
            config: self.config,
            prescribed: self.prescribed.clone(),
        })
    }
}
