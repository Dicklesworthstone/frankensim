//! Full-row P1 response using the canonical element geometry, Mandel B and K/M.
//! The caller is the accepted transient owner, not a second equilibrium solver.
use super::super::*;

pub(crate) struct TetKinematicResponse3 {
    pub stress: TetThermalStressReport,
    pub inertia_n: Vec<[f64; 3]>,
    pub damping_n: Vec<[f64; 3]>,
}

impl TetLinearElasticProblem<'_> {
    /// Internal response kernel for an already accepted, immutable P1 trajectory.
    /// Full rows retain inertia/damping on constrained nodes. No reduced matrix
    /// can supply those missing reaction rows; contract the original elements.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn recover_kinematic_response(
        &self, motion: &[Vec<[f64; 3]>; 3], thermal: Option<TetThermalStrainField<'_>>,
        time_s: f64, step: usize, eta: f64, zeta: f64, mesh_frame: ContentHash,
        maximum_strain: f64, cx: &Cx<'_>,
    ) -> Result<TetKinematicResponse3, TetElasticError> {
        let invalid = |what| TetElasticError::InvalidStressObservation { what };
        let finite = |x: f64| if x.is_finite() { Ok(x) } else {
            Err(invalid("nonfinite dynamic response arithmetic"))
        };
        cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
        self.validate_budget()?;
        self.materials.validate(self.tetrahedra.len())?;
        if let Some(field) = thermal { field.validate(self.tetrahedra.len())?; }
        if mesh_frame == ContentHash([0; 32]) || !time_s.is_finite()
            || !maximum_strain.is_finite() || maximum_strain <= 0.0
            || [eta, zeta].iter().any(|x| !x.is_finite() || *x < 0.0)
            || motion.iter().any(|x| x.len() != self.nodes_m.len()) {
            return Err(invalid("invalid dynamic response frame, motion shape or observation controls"));
        }
        let mut source = DomainHasher::new("org.frankensim.fs-solid.dynamic-motion-response.v1");
        source.update(tet_geometry_identity(self.nodes_m, self.tetrahedra)?.as_bytes());
        source.update(&time_s.to_bits().to_le_bytes());
        hash_usize(&mut source, step)?;
        for value in [eta, zeta] { source.update(&value.to_bits().to_le_bytes()); }
        for field in motion {
            for p in field {
                cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
                for &v in p { finite(v)?; source.update(&v.to_bits().to_le_bytes()); }
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
            let geometry = TetGeometry::try_new(element, tet.map(|n| self.nodes_m[n]),
                self.budget.minimum_scaled_jacobian)?;
            let material = self.materials.at(element);
            let free_strain = match thermal {
                Some(field) => {
                    let state = field.at(element);
                    if state.elastic_material_state_identity != material.material_state_identity {
                        return Err(TetElasticError::ThermalElasticStateMismatch { element });
                    }
                    source.update(&[1]); source.update(state.identity.as_bytes());
                    state.free_strain_mandel
                }
                None => { source.update(&[0]); [0.0; 6] }
            };
            source.update(&material.material_state_identity);
            source.update(&material.density_kg_m3.to_bits().to_le_bytes());
            for v in material.stiffness_mandel_pa.iter().flatten() {
                source.update(&v.to_bits().to_le_bytes());
            }
            let b = strain_displacement_matrix(&geometry);
            let total: [f64; 6] = core::array::from_fn(|i| (0..12).fold(0.0, |sum, d|
                b[i][d].mul_add(motion[0][tet[d/3]][d%3], sum)));
            let elastic_strain_mandel = core::array::from_fn::<_, 6, _>(|i| total[i]-free_strain[i]);
            if total.iter().chain(&elastic_strain_mandel)
                .any(|v| !v.is_finite() || v.abs() > maximum_strain) {
                return Err(invalid("dynamic total or elastic strain exceeds the finite small-strain ceiling"));
            }
            let stress_mandel_pa: [f64; 6] = core::array::from_fn(|i| (0..6).fold(0.0, |sum, j|
                material.stiffness_mandel_pa[i][j].mul_add(elastic_strain_mandel[j], sum)));
            let energy = finite((0..6).fold(0.0, |sum, i|
                elastic_strain_mandel[i].mul_add(stress_mandel_pa[i], sum)) * (0.5*geometry.volume))?;
            if energy < 0.0 { return Err(invalid("negative dynamic elastic strain energy")); }
            elastic_energy_j = finite(elastic_energy_j + energy)?;
            // These are the exact same consistent-mass and stiffness tiles
            // assembled by linear3. No lumping and no second constitutive law.
            let (ke, me) = element_matrices(&geometry, material);
            for d in 0..12 {
                let fi = finite((0..6).fold(0.0, |sum, i|
                    b[i][d].mul_add(stress_mandel_pa[i], sum)) * geometry.volume)?;
                let mut ma = 0.0; let mut mv = 0.0; let mut kv = 0.0;
                for j in 0..12 {
                    let node = tet[j/3]; let c = j%3;
                    ma = me[d][j].mul_add(motion[2][node][c], ma);
                    mv = me[d][j].mul_add(motion[1][node][c], mv);
                    kv = ke[d][j].mul_add(motion[1][node][c], kv);
                }
                let node = tet[d/3]; let c = d%3;
                internal[node][c] = finite(internal[node][c] + fi)?;
                inertia[node][c] = finite(inertia[node][c] + finite(ma)?)?;
                damping[node][c] = finite(damping[node][c] + finite(eta*mv + zeta*kv)?)?;
            }
            elements.push(TetElementThermalStress { elastic_strain_mandel,
                stress_mandel_pa, elastic_energy_j: energy });
        }
        let solution_identity = source.finalize();
        let mut identity = DomainHasher::new("org.frankensim.fs-solid.dynamic-stress-field.v1");
        identity.update(solution_identity.as_bytes()); identity.update(mesh_frame.as_bytes());
        identity.update(&maximum_strain.to_bits().to_le_bytes());
        for element in &elements {
            cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
            for v in element.elastic_strain_mandel.iter().chain(&element.stress_mandel_pa)
                .chain(core::iter::once(&element.elastic_energy_j)) {
                finite(*v)?; identity.update(&v.to_bits().to_le_bytes());
            }
        }
        cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
        Ok(TetKinematicResponse3 {
            stress: TetThermalStressReport { elements, nodal_internal_force_n: internal,
                elastic_energy_j, mesh_frame, solution_identity, identity: identity.finalize() },
            inertia_n: inertia, damping_n: damping,
        })
    }
}
