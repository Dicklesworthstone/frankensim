//! Material refresh in a permanently retained physical coordinate basis.
//!
//! K(T) = K_ref + T0^T [K_material(T)-K_material(ref)] T0. The original
//! numerical gauge and physical mass remain fixed, so this is a constitutive
//! update, not a time-dependent coordinate map or a hidden momentum transfer.

use super::*;

impl EnrichedMitcShell {
    /// Resolve new material stiffness while retaining EXACT original physical
    /// coordinates, reference mass, support frames and numerical drilling
    /// stiffness. Already accepted states are used directly and are not copied
    /// into the current material's static-condensation basis.
    ///
    /// The current model supplies thermal loading and constitutive recovery;
    /// its independently condensed mass and numerical gauge are not adopted.
    /// Every refresh is computed from the original reference operators, never
    /// by accumulating stiffness differences across accepted temperature steps.
    pub fn with_material_model(
        &self,
        current: Arc<MitcShellModel>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, PlateError> {
        check(&mut cancelled)?;
        if current.mass_model != ShellMass::Consistent
            || current.mesh.nodes != self.model.mesh.nodes
            || current.mesh.tris != self.model.mesh.tris
            || current.directors.vectors != self.model.directors.vectors
            || current.frames != self.model.frames
            || current.pencil.dof_map != self.model.pencil.dof_map
            || current.pencil.free != self.model.pencil.free
            || current.sections.len() != self.model.sections.len()
        {
            return Err(invalid(
                "material refresh must retain exact reference geometry, directors and support coordinates",
            ));
        }
        let mut same_stiffness = true;
        for (reference, next) in self.model.sections.iter().zip(&current.sections) {
            check(&mut cancelled)?;
            if reference.areal_mass_kg_m2 != next.areal_mass_kg_m2
                || reference.mass_first_moment_kg_m != next.mass_first_moment_kg_m
                || reference.rotary_inertia_kg != next.rotary_inertia_kg
            {
                return Err(invalid(
                    "material refresh cannot change original reference mass or inertia moments",
                ));
            }
            same_stiffness &= reference.coupled_matrices() == next.coupled_matrices();
        }
        let stiffness = if same_stiffness {
            Arc::clone(&self.reference_stiffness)
        } else {
            let mut out = Coo::new(self.dimension(), self.dimension());
            for row in 0..self.dimension() {
                row_check(row, &mut cancelled)?;
                let (columns, values) = self.reference_stiffness.row(row);
                for (&column, &value) in columns.iter().zip(values) {
                    out.push(row, column, value);
                }
            }
            for (element, field) in self.reference_fields.iter().enumerate() {
                check(&mut cancelled)?;
                let delta = field.stiffness_change(
                    &self.model.sections[element],
                    &current.sections[element],
                    &mut cancelled,
                )?;
                scatter(&mut out, &self.model, element, &delta, 1.0, &mut cancelled)?;
            }
            Arc::new(out.assemble())
        };
        check(&mut cancelled)?;
        Ok(Self {
            model: Arc::clone(&self.model),
            material_model: current,
            reference_fields: Arc::clone(&self.reference_fields),
            reference_stiffness: Arc::clone(&self.reference_stiffness),
            gauge_stiffness: Arc::clone(&self.gauge_stiffness),
            stiffness,
            mass: Arc::clone(&self.mass),
            gauges: Arc::clone(&self.gauges),
        })
    }
}

/// Scatter one symmetric enriched matrix through the ORIGINAL physical support
/// frames. The two extra columns are the unchanged element director axes.
/// Mirrored entries share exactly the same rounded value and reduction order.
pub(super) fn scatter(
    out: &mut Coo,
    model: &MitcShellModel,
    element: usize,
    matrix: &[[f64; 20]; 20],
    scale: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), PlateError> {
    let triangle = model.mesh.tris[element];
    let map: [Option<usize>; 20] = core::array::from_fn(|local| {
        if local < 18 {
            model.pencil.dof_map[6 * triangle[local / 6] + local % 6]
        } else {
            Some(model.pencil.free + 2 * element + local - 18)
        }
    });
    let columns: [[(usize, f64); 3]; 20] = core::array::from_fn(|local| {
        if local < 18 {
            let node = triangle[local / 6];
            let block = 6 * (local / 6) + 3 * ((local % 6) / 3);
            core::array::from_fn(|axis| (block + axis, model.frames[node][local % 3][axis]))
        } else {
            [(local, 1.0), (local, 0.0), (local, 0.0)]
        }
    });
    for i in 0..20 {
        check(cancelled)?;
        let Some(row) = map[i] else { continue };
        for j in i..20 {
            let Some(column) = map[j] else { continue };
            let mut value = 0.0;
            for &(a, x) in &columns[i] {
                for &(b, y) in &columns[j] {
                    value += (x * y) * matrix[a][b];
                }
            }
            let value = finite(scale * value)?;
            if value != 0.0 {
                out.push(row, column, value);
                if row != column {
                    out.push(column, row, value);
                }
            }
        }
    }
    Ok(())
}
