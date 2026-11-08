//! Fixed-temperature affine loads from the original thermal constitutive owner.
//!
//! sigma = C(S-a)-e^T E, D = eS+epsilon E+P. Here P is the CLAMPED
//! polarization increment; free expansion contributes e*a separately. Native
//! tangents are unchanged. Initial nodal force/charge and reference storage are
//! integrated once, then reused by every static load case and field recovery.

use super::*;
use fs_material::piezoelectric::thermal::ResolvedThermalTensorPiezoelectric;

pub(super) struct ThermalReference {
    pub(super) force_n: Vec<[f64; 3]>,
    pub(super) charge_c: Vec<f64>,
    pub(super) energy_j: f64,
}

impl<'a> TetPiezoelectricModel<'a> {
    /// Assemble the original reciprocal operators with source-bound thermal
    /// strain and clamped polarization at fixed cell temperatures. Samples are
    /// uniform or in connectivity order and already in `mesh_frame`.
    ///
    /// The samples retain their temperature, thermal reference and bias-source
    /// identity. Native tensor storage is reused for constitutive evaluation;
    /// there are no source queries or factorizations during field recovery.
    /// A thermal evolution owner must account for changing these parameters.
    pub fn assemble_thermal_materials(
        nodes_m: &'a [[f64; 3]],
        tetrahedra: &'a [[usize; 4]],
        materials: Vec<ResolvedThermalTensorPiezoelectric>,
        mesh_frame: ContentHash,
        budget: TetAssemblyBudget,
        cx: &Cx<'_>,
    ) -> Result<Self> {
        poll(cx)?;
        let tangents = materials
            .iter()
            .map(|sample| sample.material().clone())
            .collect();
        let mut model =
            Self::assemble_owned_materials(nodes_m, tetrahedra, tangents, mesh_frame, budget, cx)?;
        let active = materials.iter().any(|sample| {
            sample
                .thermal_strain_mandel()
                .iter()
                .any(|&value| value != 0.0)
                || sample
                    .clamped_polarization_c_m2()
                    .iter()
                    .any(|&value| value != 0.0)
        });
        if active {
            let mut reference = ThermalReference {
                force_n: vec![[0.0; 3]; nodes_m.len()],
                charge_c: vec![0.0; nodes_m.len()],
                energy_j: 0.0,
            };
            for (cell, tet) in tetrahedra.iter().enumerate() {
                poll(cx)?;
                let sample = &materials[if materials.len() == 1 { 0 } else { cell }];
                let response = sample.field_response([0.0; 6], [0.0; 3])?;
                let gradients = &model.dielectric.gradients()[cell];
                let b = strain_matrix(gradients);
                let volume = model.dielectric.volumes_m3()[cell];
                reference.energy_j += volume * response.stored_energy_j_m3;
                for local in 0..12 {
                    let value = (0..6)
                        .map(|a| b[a][local] * response.stress_mandel_pa[a])
                        .sum::<f64>();
                    reference.force_n[tet[local / 3]][local % 3] += volume * value;
                }
                for local in 0..4 {
                    reference.charge_c[tet[local]] -=
                        volume * dot3(gradients[local], response.displacement_c_m2);
                }
            }
            if !reference.energy_j.is_finite()
                || reference
                    .force_n
                    .iter()
                    .flatten()
                    .chain(&reference.charge_c)
                    .any(|value| !value.is_finite())
            {
                return Err(invalid(
                    "thermal reference force, charge or storage is not representable",
                ));
            }
            model.thermal_reference = Some(reference);
        }
        model.thermal_materials = Some(materials);
        poll(cx)?;
        Ok(model)
    }

    /// Original resolved thermal samples, uniform or in connectivity order.
    /// Their bias-source identity and reference remain available even when the
    /// current offsets vanish. Ordinary tangent-only assembly returns `None`.
    pub fn thermal_materials(&self) -> Option<&[ResolvedThermalTensorPiezoelectric]> {
        self.thermal_materials.as_deref()
    }

    pub(super) fn has_thermal_bias(&self) -> bool {
        self.thermal_reference.is_some()
    }

    pub(super) fn stamp_thermal_bias(&self, boundary: &mut BoundaryMap, cx: &Cx<'_>) -> Result<()> {
        if let Some(reference) = &self.thermal_reference {
            let nu = 3 * self.nodes_m.len();
            for node in 0..self.nodes_m.len() {
                poll(cx)?;
                for axis in 0..3 {
                    if let Some(row) = boundary.index[3 * node + axis] {
                        boundary.rhs[row] -= reference.force_n[node][axis];
                    }
                }
                if let Some(row) = boundary.index[nu + node] {
                    // Electrical saddle rows represent MINUS physical charge.
                    boundary.rhs[row] += reference.charge_c[node];
                }
            }
        }
        Ok(())
    }

    pub(super) fn elastic_strain(&self, cell: usize, strain: [f64; 6]) -> [f64; 6] {
        if self.has_thermal_bias() {
            let thermal = self.thermal_strain(cell);
            core::array::from_fn(|i| strain[i] - thermal[i])
        } else {
            strain
        }
    }

    pub(super) fn thermal_strain(&self, cell: usize) -> [f64; 6] {
        self.thermal_materials
            .as_ref()
            .map_or([0.0; 6], |materials| {
                materials[if materials.len() == 1 { 0 } else { cell }].thermal_strain_mandel()
            })
    }

    pub(super) fn cell_response(
        &self,
        cell: usize,
        strain: [f64; 6],
        field: [f64; 3],
    ) -> Result<TensorPiezoelectricResponse> {
        if self.has_thermal_bias()
            && let Some(materials) = &self.thermal_materials
        {
            Ok(materials[if materials.len() == 1 { 0 } else { cell }]
                .field_response(strain, field)?)
        } else {
            // Preserve the original response exactly at a zero-bias reference.
            Ok(
                self.materials[if self.materials.len() == 1 { 0 } else { cell }]
                    .field_response(strain, field)?,
            )
        }
    }
}
