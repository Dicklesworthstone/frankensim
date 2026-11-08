//! Uniform body acceleration in the original retained director coordinates.

use super::{
    EnrichedMitcShell, MitcShellModel, PlateError, PreparedThermalShell, check, dd_add, finite,
    invalid, row_check, two_prod,
};

impl EnrichedMitcShell {
    /// Generalized force from a uniform world acceleration [m/s²].
    ///
    /// The original physical reference inertia supplies all nodal forces,
    /// offset-mass moments and independent interior-director forces. Every
    /// translational mass column contributes before homogeneous supports are
    /// applied: multiplying the restricted mass by a restricted acceleration
    /// would incorrectly remove force from constrained columns.
    ///
    /// Three response columns are retained during the original assembly;
    /// evaluation is linear in the generalized dimension, with no quadrature,
    /// reassembly or sparse solve. Material refreshes retain this same cache.
    /// The load uses the fixed reference geometry and excludes numerical
    /// drilling inertia. It is not a finite-rotation follower or centrifugal
    /// load. The returned vector is directly accepted by the native midpoint
    /// integrator, including its ordinary load-work accounting.
    pub fn body_acceleration_load(
        &self,
        acceleration_m_s2: [f64; 3],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<f64>, PlateError> {
        check(&mut cancelled)?;
        if acceleration_m_s2.iter().any(|v| !v.is_finite()) {
            return Err(invalid("finite world body acceleration required"));
        }
        let mut load = Vec::with_capacity(self.dimension());
        for (row, response) in self.body_acceleration.iter().enumerate() {
            row_check(row, &mut cancelled)?;
            let mut sum = (0.0, 0.0);
            for axis in 0..3 {
                sum = dd_add(sum, two_prod(response[axis], acceleration_m_s2[axis]));
            }
            load.push(finite(sum.0 + sum.1)?);
        }
        self.admit_gauge_load(&load)?;
        check(&mut cancelled)?;
        Ok(load)
    }

    /// Combine thermal, explicit nodal and uniform body loads in the SAME
    /// original work-conjugate coordinates. Interior body moments remain
    /// separate generalized forces; they are never squeezed into nodal loads.
    pub fn generalized_load_with_body_acceleration(
        &self,
        thermal: &PreparedThermalShell<'_>,
        external_load_6n: &[f64],
        acceleration_m_s2: [f64; 3],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<f64>, PlateError> {
        let body = self.body_acceleration_load(acceleration_m_s2, &mut cancelled)?;
        let mut load = self.generalized_load(thermal, external_load_6n, &mut cancelled)?;
        for (row, (value, body)) in load.iter_mut().zip(body).enumerate() {
            row_check(row, &mut cancelled)?;
            *value = finite(*value + body)?;
        }
        self.admit_gauge_load(&load)?;
        check(&mut cancelled)?;
        Ok(load)
    }
}

/// Pull back rows through the original support frames AFTER taking the complete
/// rigid-translation response. The two residual-director rows are unchanged.
pub(super) fn accumulate(
    response: &mut [[f64; 3]],
    model: &MitcShellModel,
    element: usize,
    inertia: &[[f64; 20]; 20],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), PlateError> {
    let triangle = model.mesh.tris[element];
    for local in 0..20 {
        check(cancelled)?;
        let row = if local < 18 {
            let Some(row) = model.pencil.dof_map[6 * triangle[local / 6] + local % 6] else {
                continue;
            };
            row
        } else {
            model.pencil.free + 2 * element + local - 18
        };
        let terms: [(usize, f64); 3] = if local < 18 {
            let node = triangle[local / 6];
            let block = 6 * (local / 6) + 3 * ((local % 6) / 3);
            core::array::from_fn(|axis| (block + axis, model.frames[node][local % 3][axis]))
        } else {
            [(local, 1.0), (local, 0.0), (local, 0.0)]
        };
        for axis in 0..3 {
            let mut sum = (0.0, 0.0);
            for &(physical_row, coefficient) in &terms {
                if coefficient == 0.0 {
                    continue;
                }
                for node in 0..3 {
                    sum = dd_add(
                        sum,
                        two_prod(coefficient, inertia[physical_row][6 * node + axis]),
                    );
                }
            }
            response[row][axis] = finite(response[row][axis] + finite(sum.0 + sum.1)?)?;
        }
    }
    Ok(())
}
