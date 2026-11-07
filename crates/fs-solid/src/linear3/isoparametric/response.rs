//! Full-row P1 response using the canonical element geometry, Mandel B and K/M.
//! The caller is the accepted transient owner, not a second equilibrium solver.
use super::super::*;

pub(crate) struct TetKinematicResponse3 {
    pub stress: TetThermalStressReport,
    pub inertia_n: Vec<[f64; 3]>,
    pub damping_n: Vec<[f64; 3]>,
}

/// Full-row actions for prescribed-support dynamics. These retain free/fixed
/// coupling using the same element matrices as the reduced acceleration solve.
pub(crate) struct TetFullActions3 {
    pub ku: Vec<f64>,
    pub mv: Vec<f64>,
    pub kv: Vec<f64>,
    pub ma: Vec<f64>,
}

impl TetLinearElasticProblem<'_> {
    /// Fully prescribed dynamics still has physical inertia and support work.
    /// Keep the public modal assembly's NoFreeDofs contract unchanged.
    pub(crate) fn assemble_prescribed_dynamics(
        &self,
        nodal_reference_density: Option<&[f64]>,
        cx: &Cx<'_>,
    ) -> Result<TetElasticAssembly, TetElasticError> {
        match self.assemble_reference_mass(nodal_reference_density, cx) {
            Ok(assembly) => return Ok(assembly),
            Err(TetElasticError::NoFreeDofs) => {}
            Err(error) => return Err(error),
        }
        // assemble already admitted budgets/materials and every fixed DOF,
        // but its empty-free-set return precedes connectivity/geometry work.
        let mut incident = vec![false; self.nodes_m.len()];
        let mut volumes = Vec::with_capacity(self.tetrahedra.len());
        let mut minimum = f64::INFINITY;
        let mut total_mass = 0.0;
        for (element, tet) in self.tetrahedra.iter().enumerate() {
            cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
            validate_connectivity(element, tet, self.nodes_m.len())?;
            for &node in tet {
                incident[node] = true;
            }
            let geometry = TetGeometry::try_new(
                element,
                tet.map(|n| self.nodes_m[n]),
                self.budget.minimum_scaled_jacobian,
            )?;
            minimum = minimum.min(geometry.scaled_jacobian);
            total_mass += match nodal_reference_density {
                Some(density) => crate::inertia::p1_tet_nodal_mass(
                    geometry.volume, tet.map(|node| density[node]),
                ).map_err(|_| TetElasticError::InvalidMaterial {
                    what: "unrepresentable prescribed-body reference mass",
                })?.mass_kg,
                None => self.materials.at(element).density_kg_m3 * geometry.volume,
            };
            if !total_mass.is_finite() {
                return Err(TetElasticError::InvalidMaterial {
                    what: "nonfinite prescribed-body mass",
                });
            }
            volumes.push(geometry.volume);
        }
        if let Some(node) = incident.iter().position(|seen| !seen) {
            return Err(TetElasticError::UnreferencedNode { node });
        }
        cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
        Ok(TetElasticAssembly {
            stiffness: Coo::new(0, 0).assemble(),
            mass: Coo::new(0, 0).assemble(),
            free_dofs: Vec::new(),
            total_mass_kg: total_mass,
            minimum_scaled_jacobian: minimum,
            element_volumes_m3: volumes,
        })
    }

    pub(crate) fn apply_full_kinematics(
        &self,
        motion: &[Vec<f64>; 3],
        nodal_reference_density: Option<&[f64]>,
        cx: &Cx<'_>,
    ) -> Result<TetFullActions3, TetElasticError> {
        cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
        self.validate_budget()?;
        self.materials.validate(self.tetrahedra.len())?;
        let full = self
            .nodes_m
            .len()
            .checked_mul(3)
            .ok_or(TetElasticError::SizeOverflow)?;
        let finite = |value: f64| {
            if value.is_finite() {
                Ok(value)
            } else {
                Err(TetElasticError::InvalidStressObservation {
                    what: "nonfinite full dynamic operator action",
                })
            }
        };
        for field in motion {
            if field.len() != full {
                return Err(TetElasticError::InvalidStressObservation {
                    what: "full dynamic action needs node-major xyz fields",
                });
            }
            for chunk in field.chunks(256) {
                cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
                for &value in chunk {
                    finite(value)?;
                }
            }
        }
        let mut out = TetFullActions3 {
            ku: vec![0.0; full],
            mv: vec![0.0; full],
            kv: vec![0.0; full],
            ma: vec![0.0; full],
        };
        for (element, tet) in self.tetrahedra.iter().enumerate() {
            cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
            validate_connectivity(element, tet, self.nodes_m.len())?;
            let geometry = TetGeometry::try_new(
                element,
                tet.map(|n| self.nodes_m[n]),
                self.budget.minimum_scaled_jacobian,
            )?;
            let material = self.materials.at(element);
            let b = strain_displacement_matrix(&geometry);
            let me = p1_tet_reference_mass(geometry.volume, material.density_kg_m3,
                *tet, nodal_reference_density)?;
            let stiffness_action = |field: &[f64]| {
                let strain = relative_mandel_strain(
                    tet.len(),
                    |d| core::array::from_fn(|i| b[i][d]),
                    |d| field[3 * tet[d / 3] + d % 3],
                );
                let stress = core::array::from_fn(|i| {
                    (0..6).fold(0.0, |sum, j| {
                        material.stiffness_mandel_pa[i][j].mul_add(strain[j], sum)
                    })
                });
                relative_mandel_force(&b, &stress, geometry.volume)
            };
            let ku = stiffness_action(&motion[0]);
            let kv = stiffness_action(&motion[1]);
            for i in 0..12 {
                let mut mv = 0.0;
                let mut ma = 0.0;
                for j in 0..12 {
                    let column = 3 * tet[j / 3] + j % 3;
                    mv = me[i][j].mul_add(motion[1][column], mv);
                    ma = me[i][j].mul_add(motion[2][column], ma);
                }
                let row = 3 * tet[i / 3] + i % 3;
                out.ku[row] = finite(out.ku[row] + finite(ku[i])?)?;
                out.mv[row] = finite(out.mv[row] + finite(mv)?)?;
                out.kv[row] = finite(out.kv[row] + finite(kv[i])?)?;
                out.ma[row] = finite(out.ma[row] + finite(ma)?)?;
            }
        }
        cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
        Ok(out)
    }

    /// One original full storage action into caller-owned workspace. This keeps
    /// support cross rows without assembling a full matrix or evaluating the
    /// three unrelated kinematic actions needed by endpoint stress recovery.
    pub(crate) fn apply_full_storage(
        &self, field: &[f64], out: &mut [f64], mass: bool,
        nodal_reference_density: Option<&[f64]>, cx: &Cx<'_>,
    ) -> Result<(), TetElasticError> {
        cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
        self.validate_budget()?;
        self.materials.validate(self.tetrahedra.len())?;
        let full = self.nodes_m.len().checked_mul(3).ok_or(TetElasticError::SizeOverflow)?;
        if field.len() != full || out.len() != full {
            return Err(TetElasticError::InvalidStressObservation {
                what: "full storage action needs node-major xyz fields",
            });
        }
        for (input, output) in field.chunks(256).zip(out.chunks_mut(256)) {
            cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
            if input.iter().any(|v| !v.is_finite()) {
                return Err(TetElasticError::InvalidStressObservation {
                    what: "nonfinite full storage input",
                });
            }
            output.fill(0.0);
        }
        for (element, tet) in self.tetrahedra.iter().enumerate() {
            cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
            validate_connectivity(element, tet, self.nodes_m.len())?;
            let geometry = TetGeometry::try_new(element, tet.map(|n| self.nodes_m[n]),
                self.budget.minimum_scaled_jacobian)?;
            let material = self.materials.at(element);
            let action = if mass {
                let matrix = p1_tet_reference_mass(geometry.volume, material.density_kg_m3,
                    *tet, nodal_reference_density)?;
                core::array::from_fn(|i| (0..12).fold(0.0, |value, j|
                    matrix[i][j].mul_add(field[3 * tet[j / 3] + j % 3], value)))
            } else {
                let b = strain_displacement_matrix(&geometry);
                let strain = relative_mandel_strain(tet.len(),
                    |d| core::array::from_fn(|i| b[i][d]),
                    |d| field[3 * tet[d / 3] + d % 3]);
                let stress = core::array::from_fn(|i| (0..6).fold(0.0, |value, j|
                    material.stiffness_mandel_pa[i][j].mul_add(strain[j], value)));
                relative_mandel_force(&b, &stress, geometry.volume)
            };
            for i in 0..12 {
                let row = 3 * tet[i / 3] + i % 3;
                out[row] += action[i];
                if !out[row].is_finite() {
                    return Err(TetElasticError::InvalidStressObservation {
                        what: "nonfinite full storage action",
                    });
                }
            }
        }
        cx.checkpoint().map_err(|_| TetElasticError::Cancelled)
    }

    /// Actual transposed element actions, including all free/fixed cross rows.
    /// The seed is zero on fixed rows when used by the support-motion adjoint.
    /// Transpose the relative strain map and canonical material entries, rather
    /// than assuming a separately assembled stiffness tile is bitwise symmetric.
    pub(crate) fn transpose_full_kinematic_actions(
        &self,
        seed: &[f64],
        nodal_reference_density: Option<&[f64]>,
        cx: &Cx<'_>,
    ) -> Result<(Vec<f64>, Vec<f64>), TetElasticError> {
        cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
        self.validate_budget()?;
        self.materials.validate(self.tetrahedra.len())?;
        let full = self
            .nodes_m
            .len()
            .checked_mul(3)
            .ok_or(TetElasticError::SizeOverflow)?;
        let finite = |value: f64| {
            if value.is_finite() {
                Ok(value)
            } else {
                Err(TetElasticError::InvalidStressObservation {
                    what: "nonfinite transposed dynamic action",
                })
            }
        };
        if seed.len() != full {
            return Err(TetElasticError::InvalidStressObservation {
                what: "transposed dynamic action needs a full node-major cotangent",
            });
        }
        for chunk in seed.chunks(256) {
            cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
            for &value in chunk {
                finite(value)?;
            }
        }
        let mut mass = vec![0.0; full];
        let mut stiffness = vec![0.0; full];
        for (element, tet) in self.tetrahedra.iter().enumerate() {
            cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
            validate_connectivity(element, tet, self.nodes_m.len())?;
            let geometry = TetGeometry::try_new(
                element,
                tet.map(|n| self.nodes_m[n]),
                self.budget.minimum_scaled_jacobian,
            )?;
            let material = self.materials.at(element);
            let b = strain_displacement_matrix(&geometry);
            let me = p1_tet_reference_mass(geometry.volume, material.density_kg_m3,
                *tet, nodal_reference_density)?;
            let strain_seed = relative_mandel_strain(
                tet.len(),
                |d| core::array::from_fn(|i| b[i][d]),
                |d| seed[3 * tet[d / 3] + d % 3],
            );
            let stress_seed = core::array::from_fn(|j| {
                (0..6).fold(0.0, |sum, i| {
                    material.stiffness_mandel_pa[i][j].mul_add(strain_seed[i], sum)
                })
            });
            let stiffness_seed = relative_mandel_force(&b, &stress_seed, geometry.volume);
            for column in 0..12 {
                let mut mm = 0.0;
                for row in 0..12 {
                    let d = 3 * tet[row / 3] + row % 3;
                    mm = me[row][column].mul_add(seed[d], mm);
                }
                let d = 3 * tet[column / 3] + column % 3;
                mass[d] = finite(mass[d] + finite(mm)?)?;
                stiffness[d] = finite(stiffness[d] + finite(stiffness_seed[column])?)?;
            }
        }
        cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
        Ok((mass, stiffness))
    }
}

#[cfg(test)]
mod strain_tests {
    use super::*;
    use fs_alloc::{ArenaConfig, ArenaPool};
    use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};

    const NODES: [[f64; 3]; 4] = [
        [0.25, -0.3, 0.7],
        [1.6, 0.1, 0.2],
        [-0.2, 1.7, 0.3],
        [0.1, 0.2, 1.9],
    ];
    const TETS: [[usize; 4]; 1] = [[0, 1, 2, 3]];
    const FRAME: ContentHash = ContentHash([203; 32]);

    fn with_cx(run: impl FnOnce(&Cx<'_>)) {
        let gate = CancelGate::new_clock_free();
        ArenaPool::new(ArenaConfig::default()).scope(|arena| {
            run(&Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 203,
                    kernel_id: 203,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            ));
        });
    }

    fn material() -> TetElasticMaterial {
        let a = [1.0, 0.5, -0.25, 0.75, -0.5, 0.25];
        let c = core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                (if i == j { 40.0 + 8.0 * i as f64 } else { 0.0 }) + 2.0 * a[i] * a[j]
            })
        });
        TetElasticMaterial::try_new_mandel(3.25, c, [203; 32]).unwrap()
    }

    fn problem(material: &TetElasticMaterial) -> TetLinearElasticProblem<'_> {
        TetLinearElasticProblem {
            nodes_m: &NODES,
            tetrahedra: &TETS,
            materials: TetMaterialField::Uniform(material),
            fixed_dofs: &[],
            budget: TetAssemblyBudget::standard(),
        }
    }

    fn close(a: f64, b: f64) {
        assert!(
            (a - b).abs() <= 2e-11 * a.abs().max(b.abs()).max(1.0),
            "{a:e} != {b:e}"
        );
    }

    fn flatten(motion: &[Vec<[f64; 3]>; 3]) -> [Vec<f64>; 3] {
        core::array::from_fn(|i| motion[i].iter().flatten().copied().collect())
    }

    #[test]
    fn g0_rigid_translation_retains_inertia_without_elastic_strain_or_damping() {
        with_cx(|cx| {
            let law = material();
            let p = problem(&law);
            let u = [2.0_f64.powi(35), -2.0_f64.powi(36), 2.0_f64.powi(34)];
            let v = [2.0_f64.powi(29), -2.0_f64.powi(30), 2.0_f64.powi(28)];
            let a = [0.25, -0.5, 0.75];
            let motion = [vec![u; 4], vec![v; 4], vec![a; 4]];
            let response = p
                .recover_kinematic_response(&motion, None, 0.0, 0, 0.25, 1e6, FRAME, 1e-20, None, cx)
                .unwrap();
            assert_eq!(response.stress.elements[0].elastic_strain_mandel, [0.0; 6]);
            assert_eq!(response.stress.elements[0].stress_mandel_pa, [0.0; 6]);
            assert_eq!(response.stress.elastic_energy_j, 0.0);
            assert_eq!(response.stress.nodal_internal_force_n, vec![[0.0; 3]; 4]);
            let mass = law.density_kg_m3
                * TetGeometry::try_new(0, NODES, p.budget.minimum_scaled_jacobian)
                    .unwrap()
                    .volume;
            let actions = p.apply_full_kinematics(&flatten(&motion), None, cx).unwrap();
            assert_eq!(actions.ku, vec![0.0; 12]);
            assert_eq!(actions.kv, vec![0.0; 12]);
            let (mt, kt) = p
                .transpose_full_kinematic_actions(&flatten(&motion)[1], None, cx)
                .unwrap();
            assert_eq!(kt, vec![0.0; 12]);
            for node in 0..4 {
                for c in 0..3 {
                    close(response.inertia_n[node][c], 0.25 * mass * a[c]);
                    close(response.damping_n[node][c], 0.25 * 0.25 * mass * v[c]);
                    close(actions.mv[3 * node + c], 0.25 * mass * v[c]);
                    close(actions.ma[3 * node + c], 0.25 * mass * a[c]);
                    close(mt[3 * node + c], 0.25 * mass * v[c]);
                }
            }
            let (scales, damping) = p
                .transient_material_contraction(&motion, &vec![v; 4], None, 0.25, 1e6, None, cx)
                .unwrap();
            assert_eq!(scales[0][0], 0.0);
            assert_eq!(damping[1], 0.0);
        });
    }

    #[test]
    fn g3_translation_preserves_thermal_strain_force_energy_and_pullbacks() {
        with_cx(|cx| {
            let law = material();
            let p = problem(&law);
            let thermal = TetThermalStrainState::try_new(
                301.0,
                300.0,
                [1.0, -2.0, 1.0, 0.5, -0.5, 0.25].map(|x| x / 8192.0),
                law.material_state_identity,
                FRAME,
                0.05,
            )
            .unwrap();
            let field = Some(TetThermalStrainField::Uniform(&thermal));
            let motion: [Vec<[f64; 3]>; 3] = core::array::from_fn(|s| {
                (0..4)
                    .map(|n| {
                        core::array::from_fn(|c| ((s + 1) * (n + 1) * (c + 2)) as f64 / 4096.0)
                    })
                    .collect()
            });
            let shift = [2.0_f64.powi(35), -2.0_f64.powi(36), 2.0_f64.powi(34)];
            let mut translated = motion.clone();
            for state in &mut translated[..2] {
                for node in state {
                    for c in 0..3 {
                        node[c] += shift[c];
                    }
                }
            }
            let recover = |m: &[Vec<[f64; 3]>; 3]| {
                p.recover_kinematic_response(m, field, 0.0, 0, 0.0, 0.2, FRAME, 0.05, None, cx)
                    .unwrap()
            };
            let a = recover(&motion);
            let b = recover(&translated);
            assert_eq!(a.stress.elements[0], b.stress.elements[0]);
            assert_eq!(
                a.stress.nodal_internal_force_n,
                b.stress.nodal_internal_force_n
            );
            assert_eq!(a.damping_n, b.damping_n);
            assert_eq!(a.inertia_n, b.inertia_n);
            let plain = p.apply_full_kinematics(&flatten(&motion), None, cx).unwrap();
            let shifted = p.apply_full_kinematics(&flatten(&translated), None, cx).unwrap();
            assert_eq!(plain.ku, shifted.ku);
            assert_eq!(plain.kv, shifted.kv);

            let weights = [0.75, -0.5, 0.25, 1.0, -0.25, 0.5];
            let pullback = p
                .stress_observation_contraction(&motion[0], &[(0, weights)], field, 0.05, cx)
                .unwrap();
            let translated_pullback = p
                .stress_observation_contraction(&translated[0], &[(0, weights)], field, 0.05, cx)
                .unwrap();
            assert_eq!(pullback, translated_pullback);
            let h = 1.0 / 1024.0;
            let observation = |step: f64| {
                let mut trial = motion.clone();
                for n in 0..4 {
                    for c in 0..3 {
                        trial[0][n][c] += step * motion[1][n][c];
                    }
                }
                recover(&trial).stress.elements[0]
                    .stress_mandel_pa
                    .iter()
                    .zip(weights)
                    .map(|(s, w)| s * w)
                    .sum::<f64>()
            };
            let directional: f64 = pullback
                .0
                .iter()
                .flatten()
                .zip(motion[1].iter().flatten())
                .map(|(x, y)| x * y)
                .sum();
            close(directional, (observation(h) - observation(-h)) / (2.0 * h));

            let seed: Vec<_> = motion[2].iter().flatten().copied().collect();
            let (_, kt) = p.transpose_full_kinematic_actions(&seed, None, cx).unwrap();
            close(
                seed.iter().zip(&plain.ku).map(|(x, y)| x * y).sum(),
                flatten(&motion)[0].iter().zip(kt).map(|(x, y)| x * y).sum(),
            );
            let gradient = p
                .transient_material_contraction(&motion, &motion[2], field, 0.0, 0.2, None, cx)
                .unwrap();
            let translated_gradient = p
                .transient_material_contraction(&translated, &translated[1], field, 0.0, 0.2, None, cx)
                .unwrap();
            // Use the same relative cotangent for the translated comparison.
            let base_gradient = p
                .transient_material_contraction(&motion, &motion[1], field, 0.0, 0.2, None, cx)
                .unwrap();
            assert_eq!(base_gradient.0[0][0], translated_gradient.0[0][0]);
            assert_eq!(base_gradient.1[1], translated_gradient.1[1]);
            let force_contraction: f64 = motion[2]
                .iter()
                .flatten()
                .zip(
                    a.stress
                        .nodal_internal_force_n
                        .iter()
                        .flatten()
                        .zip(a.damping_n.iter().flatten()),
                )
                .map(|(seed, (internal, damping))| -seed * (internal + damping))
                .sum();
            close(gradient.0[0][0], force_contraction);
        });
    }
}

impl TetLinearElasticProblem<'_> {
    /// Internal response kernel for an already accepted, immutable P1 trajectory.
    /// Full rows retain inertia/damping on constrained nodes. No reduced matrix
    /// can supply those missing reaction rows; contract the original elements.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn recover_kinematic_response(
        &self,
        motion: &[Vec<[f64; 3]>; 3],
        thermal: Option<TetThermalStrainField<'_>>,
        time_s: f64,
        step: usize,
        eta: f64,
        zeta: f64,
        mesh_frame: ContentHash,
        maximum_strain: f64,
        nodal_reference_density: Option<&[f64]>,
        cx: &Cx<'_>,
    ) -> Result<TetKinematicResponse3, TetElasticError> {
        let invalid = |what| TetElasticError::InvalidStressObservation { what };
        let finite = |x: f64| {
            if x.is_finite() {
                Ok(x)
            } else {
                Err(invalid("nonfinite dynamic response arithmetic"))
            }
        };
        cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
        self.validate_budget()?;
        self.materials.validate(self.tetrahedra.len())?;
        if let Some(field) = thermal {
            field.validate(self.tetrahedra.len())?;
        }
        if mesh_frame == ContentHash([0; 32])
            || !time_s.is_finite()
            || !maximum_strain.is_finite()
            || maximum_strain <= 0.0
            || [eta, zeta].iter().any(|x| !x.is_finite() || *x < 0.0)
            || motion.iter().any(|x| x.len() != self.nodes_m.len())
        {
            return Err(invalid(
                "invalid dynamic response frame, motion shape or observation controls",
            ));
        }
        let mut source = DomainHasher::new("org.frankensim.fs-solid.dynamic-motion-response.v1");
        source.update(tet_geometry_identity(self.nodes_m, self.tetrahedra)?.as_bytes());
        source.update(&time_s.to_bits().to_le_bytes());
        if let Some(density) = nodal_reference_density {
            source.update(b"nodal-reference-density");
            for value in density {
                cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
                source.update(&value.to_bits().to_le_bytes());
            }
        }
        hash_usize(&mut source, step)?;
        for value in [eta, zeta] {
            source.update(&value.to_bits().to_le_bytes());
        }
        for field in motion {
            for p in field {
                cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
                for &v in p {
                    finite(v)?;
                    source.update(&v.to_bits().to_le_bytes());
                }
            }
        }
        let mut elements = Vec::with_capacity(self.tetrahedra.len());
        let mut internal = vec![[0.0; 3]; self.nodes_m.len()];
        let mut inertia = internal.clone();
        let mut damping = internal.clone();
        let mut elastic_energy_j = 0.0;
        for (element, tet) in self.tetrahedra.iter().enumerate() {
            cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
            validate_connectivity(element, tet, self.nodes_m.len())?;
            let geometry = TetGeometry::try_new(
                element,
                tet.map(|n| self.nodes_m[n]),
                self.budget.minimum_scaled_jacobian,
            )?;
            let material = self.materials.at(element);
            let free_strain = match thermal {
                Some(field) => {
                    let state = field.at(element);
                    if state.elastic_material_state_identity != material.material_state_identity {
                        return Err(TetElasticError::ThermalElasticStateMismatch { element });
                    }
                    source.update(&[1]);
                    source.update(state.identity.as_bytes());
                    state.free_strain_mandel
                }
                None => {
                    source.update(&[0]);
                    [0.0; 6]
                }
            };
            source.update(&material.material_state_identity);
            source.update(&material.density_kg_m3.to_bits().to_le_bytes());
            for v in material.stiffness_mandel_pa.iter().flatten() {
                source.update(&v.to_bits().to_le_bytes());
            }
            let b = strain_displacement_matrix(&geometry);
            let total = relative_mandel_strain(
                tet.len(),
                |d| core::array::from_fn(|i| b[i][d]),
                |d| motion[0][tet[d / 3]][d % 3],
            );
            let elastic_strain_mandel =
                core::array::from_fn::<_, 6, _>(|i| total[i] - free_strain[i]);
            if total
                .iter()
                .chain(&elastic_strain_mandel)
                .any(|v| !v.is_finite() || v.abs() > maximum_strain)
            {
                return Err(invalid(
                    "dynamic total or elastic strain exceeds the finite small-strain ceiling",
                ));
            }
            let stress_mandel_pa: [f64; 6] = core::array::from_fn(|i| {
                (0..6).fold(0.0, |sum, j| {
                    material.stiffness_mandel_pa[i][j].mul_add(elastic_strain_mandel[j], sum)
                })
            });
            let energy = finite(
                (0..6).fold(0.0, |sum, i| {
                    elastic_strain_mandel[i].mul_add(stress_mandel_pa[i], sum)
                }) * (0.5 * geometry.volume),
            )?;
            if energy < 0.0 {
                return Err(invalid("negative dynamic elastic strain energy"));
            }
            elastic_energy_j = finite(elastic_energy_j + energy)?;
            // The same constitutive map supplies stiffness damping; only mass
            // acts on absolute motion. Uniform velocity still carries momentum
            // and mass-proportional damping, but no elastic strain rate.
            let rate = relative_mandel_strain(
                tet.len(),
                |d| core::array::from_fn(|i| b[i][d]),
                |d| motion[1][tet[d / 3]][d % 3],
            );
            let rate_stress = core::array::from_fn(|i| {
                (0..6).fold(0.0, |sum, j| {
                    material.stiffness_mandel_pa[i][j].mul_add(rate[j], sum)
                })
            });
            let force = relative_mandel_force(&b, &stress_mandel_pa, geometry.volume);
            let kv = relative_mandel_force(&b, &rate_stress, geometry.volume);
            let me = p1_tet_reference_mass(geometry.volume, material.density_kg_m3,
                *tet, nodal_reference_density)?;
            for d in 0..12 {
                let fi = finite(force[d])?;
                let mut ma = 0.0;
                let mut mv = 0.0;
                for j in 0..12 {
                    let node = tet[j / 3];
                    let c = j % 3;
                    ma = me[d][j].mul_add(motion[2][node][c], ma);
                    mv = me[d][j].mul_add(motion[1][node][c], mv);
                }
                let node = tet[d / 3];
                let c = d % 3;
                internal[node][c] = finite(internal[node][c] + fi)?;
                inertia[node][c] = finite(inertia[node][c] + finite(ma)?)?;
                damping[node][c] = finite(damping[node][c] + finite(eta * mv + zeta * kv[d])?)?;
            }
            elements.push(TetElementThermalStress {
                elastic_strain_mandel,
                stress_mandel_pa,
                elastic_energy_j: energy,
                material_state_identity: ContentHash(material.material_state_identity),
            });
        }
        let solution_identity = source.finalize();
        let mut identity = DomainHasher::new("org.frankensim.fs-solid.dynamic-stress-field.v1");
        identity.update(solution_identity.as_bytes());
        identity.update(mesh_frame.as_bytes());
        identity.update(&maximum_strain.to_bits().to_le_bytes());
        for element in &elements {
            cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
            for v in element
                .elastic_strain_mandel
                .iter()
                .chain(&element.stress_mandel_pa)
                .chain(core::iter::once(&element.elastic_energy_j))
            {
                finite(*v)?;
                identity.update(&v.to_bits().to_le_bytes());
            }
        }
        cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
        Ok(TetKinematicResponse3 {
            stress: TetThermalStressReport {
                elements,
                nodal_internal_force_n: internal,
                elastic_energy_j,
                mesh_frame,
                solution_identity,
                identity: identity.finalize(),
            },
            inertia_n: inertia,
            damping_n: damping,
        })
    }
}
