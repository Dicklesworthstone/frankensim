//! Original three-dimensional material response and shared Gauss iteration seam.
use super::*;

impl gauss::Problem for SpatialNonlinearElectrostatics {
    type Evaluation = Evaluation;

    fn new_evaluation(&self) -> Evaluation {
        let cells = self.materials.len();
        let nodes = self.mesh().vertex_count;
        Evaluation {
            electric: vec![[0.0; 3]; cells],
            displacement: vec![[0.0; 3]; cells],
            tangent: vec![[[0.0; 3]; 3]; cells],
            cell_energy: vec![0.0; cells],
            cell_coenergy: vec![0.0; cells],
            nodal: vec![0.0; nodes],
            correction: vec![0.0; nodes],
            energy: 0.0,
            coenergy: 0.0,
        }
    }

    fn evaluate(
        &self,
        potential: &[f64],
        out: &mut Evaluation,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), Error> {
        out.nodal.fill(0.0);
        out.correction.fill(0.0);
        let (mut energy, mut ce, mut coenergy, mut cc) = (0.0, 0.0, 0.0, 0.0);
        for cell in 0..self.materials.len() {
            poll(cancelled)?;
            let element = self.element(cell);
            let e = element.field(potential)?;
            let magnitude = finite(e[0].hypot(e[1]).hypot(e[2]))?;
            let (d, tangent, u, w) = if let Some(linear) = self.linear[cell] {
                let d = linear.displacement_c_m2(e);
                let d = [finite(d[0])?, finite(d[1])?, finite(d[2])?];
                let u = finite(0.5 * dot(e, d))?;
                (d, linear.matrix_f_m(), u, u)
            } else {
                let SpatialDielectricMaterial::Reversible(law) = &self.materials[cell] else {
                    unreachable!("prepared original linear law")
                };
                if magnitude > law.maximum_field_strength_v_m() {
                    return Err(Error::FieldRange {
                        cell,
                        observed_v_m: magnitude,
                        maximum_v_m: law.maximum_field_strength_v_m(),
                    });
                }
                let response = law
                    .evaluate(magnitude)
                    .map_err(ElectrostaticError::Material)?;
                let p = response.differential_permittivity_f_m;
                let s = response.secant_permittivity_f_m;
                let mut tangent = [[0.0; 3]; 3];
                let d = if magnitude == 0.0 {
                    for (i, row) in tangent.iter_mut().enumerate() {
                        row[i] = p;
                    }
                    [0.0; 3]
                } else {
                    let n = e.map(|v| v / magnitude);
                    for i in 0..3 {
                        // Positive diagonal terms avoid subtracting nearly
                        // equal secant and longitudinal contributions.
                        let transverse =
                            n[(i + 1) % 3].mul_add(n[(i + 1) % 3], n[(i + 2) % 3] * n[(i + 2) % 3]);
                        tangent[i][i] = (p * n[i]).mul_add(n[i], s * transverse);
                        for j in 0..i {
                            let value = (p - s) * n[i] * n[j];
                            tangent[i][j] = value;
                            tangent[j][i] = value;
                        }
                    }
                    [
                        finite(response.displacement_c_m2 * n[0])?,
                        finite(response.displacement_c_m2 * n[1])?,
                        finite(response.displacement_c_m2 * n[2])?,
                    ]
                };
                (
                    d,
                    tangent,
                    response.energy_density_j_m3,
                    response.coenergy_density_j_m3,
                )
            };
            let u = finite(element.volume_m3 * u)?;
            let w = finite(element.volume_m3 * w)?;
            if u < 0.0 || w < 0.0 || (magnitude > 0.0 && (u == 0.0 || w == 0.0)) {
                return Err(invalid(
                    "tetrahedral dielectric storage is not representable",
                ));
            }
            // Native FEEC interval admission also applies to an accepted
            // zero-Newton field, whose tangent has not passed through CG.
            CellWeight::Tensor(tangent)
                .validate(cell, 1)
                .map_err(weighted)?;
            for (vertex, charge) in element.vertices.iter().copied().zip(element.reactions(d)?) {
                let v = vertex as usize;
                add(&mut out.nodal[v], &mut out.correction[v], charge)?;
            }
            out.electric[cell] = e;
            out.displacement[cell] = d;
            out.tangent[cell] = tangent;
            out.cell_energy[cell] = u;
            out.cell_coenergy[cell] = w;
            add(&mut energy, &mut ce, u)?;
            add(&mut coenergy, &mut cc, w)?;
        }
        for (charge, correction) in out.nodal.iter_mut().zip(&out.correction) {
            *charge = finite(*charge + *correction)?;
        }
        out.energy = finite(energy + ce)?;
        out.coenergy = finite(coenergy + cc)?;
        Ok(())
    }

    fn nodal_reactions(evaluation: &Evaluation) -> &[f64] {
        &evaluation.nodal
    }

    fn coenergy_change(
        &self,
        before: &Evaluation,
        after: &Evaluation,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, Error> {
        let (mut sum, mut correction) = (0.0, 0.0);
        for cell in 0..self.materials.len() {
            poll(cancelled)?;
            let a = before.electric[cell];
            let b = after.electric[cell];
            let change = if let Some(linear) = self.linear[cell] {
                let delta: [f64; 3] = std::array::from_fn(|i| b[i] - a[i]);
                let midpoint = std::array::from_fn(|i| f64::midpoint(a[i], b[i]));
                let offset = std::array::from_fn(|i| 0.5 * delta[i] - (midpoint[i] - a[i]));
                let d = linear.displacement_c_m2(midpoint);
                let rounding = linear.displacement_c_m2(offset);
                let (mut value, mut cv) = (0.0, 0.0);
                for i in 0..3 {
                    add(&mut value, &mut cv, delta[i] * d[i])?;
                    add(&mut value, &mut cv, delta[i] * rounding[i])?;
                }
                finite(value + cv)?
            } else {
                let SpatialDielectricMaterial::Reversible(law) = &self.materials[cell] else {
                    unreachable!("prepared original linear law")
                };
                law.coenergy_change_vector_density_j_m3(a, b)
                    .map_err(ElectrostaticError::Material)?
            };
            add(
                &mut sum,
                &mut correction,
                self.element(cell).volume_m3 * change,
            )?;
        }
        Ok(finite(sum + correction)?)
    }

    fn direction(
        &self,
        evaluation: &Evaluation,
        dofs: &Dofs,
        rhs: &[f64],
        linear: CheckedCgConfig,
        maximum_applications: usize,
        applications: &mut usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Vec<f64>, SolveReport), Error> {
        tangent::direction(
            self,
            &evaluation.tangent,
            dofs,
            rhs,
            linear,
            maximum_applications,
            applications,
            cancelled,
        )
    }
}
