//! Retained native material actions for repeated midpoint Newton/Krylov work.
//!
//! Opt-in local storage holds one original constitutive Hessian per tetrahedron.
//! Follower pressure and gravity retain their own sparse derivatives. Global
//! sparse storage remains the default, and the energy-momentum endpoint Jacobian
//! keeps its existing assembled owner for either ordinary-storage selection.

use fs_solid::hyper3::prepared::PreparedHyperTetTangent;

use super::{
    Coo, Csr, Cx, DynamicOperator, FiniteSolidTangentStorage, HyperTetError, checkpoint, inflate,
    reduce_matrix,
};

pub(super) enum CachedTangent<'a> {
    Material {
        operator: PreparedHyperTetTangent<'a>,
        external: Option<Csr>,
        diagonal: Vec<f64>,
        direction: Vec<f64>,
        force: Vec<f64>,
        scratch: Vec<f64>,
        external_force: Vec<f64>,
    },
    Assembled(Csr),
}

impl CachedTangent<'_> {
    pub(super) fn diagonal(&self, row: usize) -> f64 {
        match self {
            Self::Material { diagonal, .. } => diagonal[row],
            Self::Assembled(matrix) => matrix.get(row, row),
        }
    }

    pub(super) fn apply(
        &mut self,
        input: &[f64],
        output: &mut [f64],
        free: &[usize],
        cx: &Cx<'_>,
    ) -> Result<(), HyperTetError> {
        checkpoint(cx)?;
        if input.len() != free.len()
            || output.len() != free.len()
            || input.iter().any(|value| !value.is_finite())
        {
            return Err(invalid("finite free-DOF tangent direction required"));
        }
        match self {
            Self::Assembled(matrix) => matrix.spmv(input, output),
            Self::Material {
                operator,
                external,
                direction,
                force,
                scratch,
                external_force,
                ..
            } => {
                // A variation of a prescribed coordinate is zero, including
                // when its held displacement itself is nonzero.
                direction.fill(0.0);
                for (&dof, &value) in free.iter().zip(input) {
                    direction[dof] = value;
                }
                operator.apply(direction, force, scratch, cx)?;
                if let Some(matrix) = external {
                    matrix.spmv(input, external_force);
                    for (row, (&dof, value)) in free.iter().zip(output).enumerate() {
                        *value = force[dof] + external_force[row];
                    }
                } else {
                    for (&dof, value) in free.iter().zip(output) {
                        *value = force[dof];
                    }
                }
            }
        }
        checkpoint(cx)
    }
}

impl<'a, 'b, 'cx> DynamicOperator<'a, 'b, 'cx> {
    pub(super) fn ensure_tangent(&self, q: &[f64]) -> Result<(), HyperTetError> {
        let rebuild = self
            .tangent
            .borrow()
            .as_ref()
            .is_none_or(|(at, _)| at.as_slice() != q);
        if !rebuild {
            return Ok(());
        }
        let displacement = inflate(q, self.free, self.fixed);
        let tangent = if let Some(before) = self.energy_momentum_start {
            let after: Vec<f64> = displacement
                .iter()
                .zip(before)
                .map(|(mid, old)| 2.0 * mid - old)
                .collect();
            let full = match self.pressure {
                Some((surface, pressure_pa)) => {
                    self.problem.energy_momentum_tangent_with_pressure(
                        before,
                        &after,
                        surface,
                        pressure_pa,
                        self.cx,
                    )?
                }
                None => self
                    .problem
                    .energy_momentum_tangent(before, &after, self.cx)?,
            };
            // This owner retains the endpoint/midpoint factor for gravity.
            // The callback applies the existing factor when gravity is absent.
            let full = self.combine_gravity_tangent(full, &displacement)?;
            CachedTangent::Assembled(reduce_matrix(&full, self.free, self.cx)?)
        } else {
            match self.tangent_storage {
                FiniteSolidTangentStorage::Assembled => {
                    self.assembled_material_tangent(&displacement)?
                }
                FiniteSolidTangentStorage::LocalHessians {
                    maximum_tangent_bytes,
                } => self.material_tangent(&displacement, maximum_tangent_bytes)?,
            }
        };
        checkpoint(self.cx)?;
        self.tangent.replace(Some((q.to_vec(), tangent)));
        Ok(())
    }

    fn assembled_material_tangent(
        &self,
        displacement: &[f64],
    ) -> Result<CachedTangent<'b>, HyperTetError> {
        let mut full = self.problem.tangent_matrix(displacement, self.cx)?;
        if let Some((surface, pressure_pa)) = self.pressure {
            let external = surface.force_tangent(displacement, pressure_pa, self.cx)?;
            let mut combined = Coo::new(full.nrows(), full.nrows());
            for row in 0..full.nrows() {
                checkpoint(self.cx)?;
                for (matrix, scale) in [(&full, 1.0), (&external, -1.0)] {
                    let (columns, values) = matrix.row(row);
                    for (&column, &value) in columns.iter().zip(values) {
                        combined.push(row, column, scale * value);
                    }
                }
            }
            full = combined.assemble();
            for row in 0..full.nrows() {
                checkpoint(self.cx)?;
                if full.row(row).1.iter().any(|value| !value.is_finite()) {
                    return Err(invalid(
                        "native solid minus pressure tangent is not representable",
                    ));
                }
            }
        }
        let full = self.combine_gravity_tangent(full, displacement)?;
        Ok(CachedTangent::Assembled(reduce_matrix(
            &full, self.free, self.cx,
        )?))
    }

    fn material_tangent(
        &self,
        displacement: &[f64],
        maximum_tangent_bytes: usize,
    ) -> Result<CachedTangent<'b>, HyperTetError> {
        // The immutable reference problem already admitted its element budget.
        // The retained payload is 81 f64s per cell, with no 12x12 cell matrices
        // or global material triplets. Krylov actions reuse all vector storage.
        let problem = self.problem;
        let operator = problem.linearized_tangent(displacement, maximum_tangent_bytes, self.cx)?;
        let mut force = vec![0.0; displacement.len()];
        let mut scratch = vec![0.0; displacement.len()];
        operator.diagonal(&mut force, &mut scratch, self.cx)?;
        let mut diagonal: Vec<_> = self.free.iter().map(|&dof| force[dof]).collect();
        let external = if self.pressure.is_none() && self.gravity.is_none() {
            None
        } else {
            let mut correction = Coo::new(displacement.len(), displacement.len());
            if let Some((surface, pressure_pa)) = self.pressure {
                let pressure = surface.force_tangent(displacement, pressure_pa, self.cx)?;
                for row in 0..pressure.nrows() {
                    checkpoint(self.cx)?;
                    let (columns, values) = pressure.row(row);
                    for (&column, &value) in columns.iter().zip(values) {
                        correction.push(row, column, -value);
                    }
                }
            }
            let correction = self.combine_gravity_tangent(correction.assemble(), displacement)?;
            let correction = reduce_matrix(&correction, self.free, self.cx)?;
            for (row, value) in diagonal.iter_mut().enumerate() {
                *value += correction.get(row, row);
            }
            Some(correction)
        };
        if diagonal.iter().any(|value| !value.is_finite()) {
            return Err(invalid(
                "native dynamic tangent diagonal is not representable",
            ));
        }
        let external_force = if external.is_some() {
            vec![0.0; self.free.len()]
        } else {
            Vec::new()
        };
        Ok(CachedTangent::Material {
            operator,
            external,
            diagonal,
            direction: vec![0.0; displacement.len()],
            force,
            scratch,
            external_force,
        })
    }
}

fn invalid(what: &str) -> HyperTetError {
    HyperTetError::InvalidInput { what: what.into() }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_alloc::{ArenaConfig, ArenaPool};
    use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
    use fs_material::hyper::cards::HyperelasticMaterial;
    use fs_material::hyper::{Hyperelastic, HyperelasticModel};
    use fs_solid::{HyperTetProblem, HyperTetSettings, HyperelasticLawField, TetAssemblyBudget};
    use std::cell::Cell;

    struct CountingLaw {
        law: Hyperelastic,
        tangent_calls: Cell<usize>,
    }

    impl HyperelasticMaterial for CountingLaw {
        fn hyperelastic(&self) -> &Hyperelastic {
            &self.law
        }
        fn consistent_tangent(
            &self,
            f: &[f64; 9],
        ) -> Result<[[f64; 9]; 9], fs_material::MaterialError> {
            self.tangent_calls.set(self.tangent_calls.get() + 1);
            self.law.tangent(f)
        }
    }

    #[test]
    fn cached_free_actions_and_diagonals_retain_external_derivatives_without_new_material_work() {
        let gate = CancelGate::new_clock_free();
        ArenaPool::new(ArenaConfig::default()).scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 21,
                    kernel_id: 109,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            let nodes = [
                [0.0, 0.0, 0.0],
                [1.0, 0.1, 0.0],
                [0.0, 1.2, 0.0],
                [0.2, 0.0, 0.9],
            ];
            let tets = [[0, 1, 2, 3]];
            let law = CountingLaw {
                law: Hyperelastic::new(
                    HyperelasticModel::NeoHookean {
                        mu: 3.0,
                        lambda: 5.0,
                    },
                    4.0,
                )
                .unwrap(),
                tangent_calls: Cell::new(0),
            };
            let definition = HyperTetProblem {
                nodes_m: &nodes,
                tetrahedra: &tets,
                material: HyperelasticLawField::Uniform(&law),
                prescribed_m: &[(0, 0.3)],
                nodal_forces_n: &[],
                body_force_n_m3: [0.0; 3],
                budget: TetAssemblyBudget::standard(),
                settings: HyperTetSettings::default(),
            };
            let prepared = definition.prepare(&cx).unwrap();
            let u = [
                0.3, 0.0, 0.0, 0.4, 0.02, 0.0, 0.3, -0.04, 0.01, 0.32, 0.0, 0.08,
            ];
            let free = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11];
            let original = prepared.tangent_matrix(&u, &cx).unwrap();
            let reduced = reduce_matrix(&original, &free, &cx).unwrap();
            for with_external in [false, true] {
                let before = law.tangent_calls.get();
                let operator = prepared.linearized_tangent(&u, 81 * 8, &cx).unwrap();
                assert_eq!(law.tangent_calls.get(), before + 1);
                let mut force = vec![0.0; 12];
                let mut scratch = vec![0.0; 12];
                operator.diagonal(&mut force, &mut scratch, &cx).unwrap();
                let mut diagonal: Vec<_> = free.iter().map(|&dof| force[dof]).collect();
                let external = if with_external {
                    let mut matrix = Coo::new(free.len(), free.len());
                    matrix.push(0, 0, -0.7);
                    matrix.push(0, 3, 0.25);
                    matrix.push(3, 0, -0.5);
                    matrix.push(3, 3, 1.1);
                    Some(matrix.assemble())
                } else {
                    None
                };
                for row in 0..free.len() {
                    diagonal[row] += external.as_ref().map_or(0.0, |m| m.get(row, row));
                }
                let mut cache = CachedTangent::Material {
                    operator,
                    external: external.clone(),
                    diagonal,
                    direction: vec![0.0; 12],
                    force,
                    scratch,
                    external_force: vec![0.0; free.len()],
                };
                for scale in [1e-25, 1.0, -3.0] {
                    let input: Vec<_> = (0..free.len())
                        .map(|i| scale * (0.13 * i as f64 - 0.2))
                        .collect();
                    let mut actual = vec![7.0; free.len()];
                    let mut expected = vec![0.0; free.len()];
                    reduced.spmv(&input, &mut expected);
                    if let Some(matrix) = &external {
                        let mut extra = vec![0.0; free.len()];
                        matrix.spmv(&input, &mut extra);
                        for (value, extra) in expected.iter_mut().zip(extra) {
                            *value += extra;
                        }
                    }
                    cache.apply(&input, &mut actual, &free, &cx).unwrap();
                    let norm = expected.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
                    for (a, b) in actual.iter().zip(&expected) {
                        assert!((a - b).abs() <= 5e-13 * norm);
                    }
                    assert_eq!(law.tangent_calls.get(), before + 1);
                }
                for row in 0..free.len() {
                    let expected =
                        reduced.get(row, row) + external.as_ref().map_or(0.0, |m| m.get(row, row));
                    assert!(
                        (cache.diagonal(row) - expected).abs() <= 5e-13 * (1.0 + expected.abs())
                    );
                }
                let mut untouched = vec![123.0; free.len()];
                assert!(cache.apply(&[0.0; 1], &mut untouched, &free, &cx).is_err());
                assert!(untouched.iter().all(|&v| v == 123.0));
            }
        });
    }
}
