//! Uniform reference-frame acceleration through the original full shell mass.
//!
//! Preparation accompanies native panel, bonded-beam and MITC assembly. Each of the
//! three cached columns is the full mass acting on a unit rigid translation,
//! before support restriction. Constrained columns therefore still contribute
//! to free loads, and constrained rows retain their applied force and moment.
//! Eccentric reinforcement uses the original reciprocal rigid-arm transform.
//!
//! This is a dead acceleration load about the declared reference geometry.
//! It does not update offsets, density or directors with finite rotation and
//! does not infer gravity, centrifugal acceleration or support motion.

use crate::PlateError;

/// Immutable response of one assembled reference inventory to uniform acceleration.
///
/// Returned with the graded reinforced-shell pencil or retained by the native
/// MITC shell. The private coefficients include every original coordinate;
/// no reduced-mass reconstruction or density resampling occurs during evaluation.
#[derive(Debug, Clone)]
pub struct ShellBodyAcceleration {
    response: Vec<[f64; 3]>,
}

impl ShellBodyAcceleration {
    /// Complete reference-node loads: three Cartesian forces [N], then three
    /// axial moments [N m] per node, including supported coordinates.
    ///
    /// A positive acceleration produces force in that same direction. Supply
    /// gravity explicitly in global Cartesian axes [m/s²]. Evaluation costs
    /// one three-term scalar product per full coordinate and does not reassemble mass.
    /// These loads may be added to other dead loads before the original support
    /// reduction; multiplying the reduced mass by free translations omits the
    /// consistent-mass contribution from supported columns.
    ///
    /// # Errors
    /// Nonfinite acceleration or unrepresentable force/moment. No partial load
    /// vector is returned, and the prepared inventory remains reusable.
    pub fn full_load(&self, acceleration_m_s2: [f64; 3]) -> Result<Vec<f64>, PlateError> {
        self.full_load_with_cancellation(acceleration_m_s2, || false)
    }

    /// The same complete load, with cancellation before allocation, every
    /// bounded group of coordinates, and before publication.
    ///
    /// # Errors
    /// Cancellation or the same physical/numerical refusals as [`Self::full_load`].
    /// Neither a partial vector nor a changed inventory is published.
    pub fn full_load_with_cancellation(
        &self,
        acceleration_m_s2: [f64; 3],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<f64>, PlateError> {
        if cancelled() {
            return Err(PlateError::Cancelled);
        }
        if acceleration_m_s2.iter().any(|x| !x.is_finite()) {
            return Err(bad("shell body acceleration must be finite in m/s^2"));
        }
        let mut load = Vec::with_capacity(self.response.len());
        for (index, row) in self.response.iter().enumerate() {
            if index % 256 == 0 && cancelled() {
                return Err(PlateError::Cancelled);
            }
            let value = row
                .iter()
                .zip(acceleration_m_s2)
                .map(|(m, a)| m * a)
                .sum::<f64>();
            if !value.is_finite() {
                return Err(bad("shell body acceleration load overflow"));
            }
            load.push(value);
        }
        if cancelled() {
            return Err(PlateError::Cancelled);
        }
        Ok(load)
    }

    pub(super) fn empty() -> Self {
        Self {
            response: Vec::new(),
        }
    }

    // Called only after the native shell geometry, sections, reference inventory and
    // support map have been admitted, before element assembly.
    pub(super) fn initialize(&mut self, coordinates: usize) {
        self.response = vec![[0.0; 3]; coordinates];
    }

    // The native assembler supplies each already transformed symmetric mass
    // entry once. Only translational columns act on a rigid acceleration;
    // their force may lie on either a translational or rotational row.
    pub(super) fn accumulate_symmetric(
        &mut self,
        row: usize,
        column: usize,
        value: f64,
    ) -> Result<(), PlateError> {
        if column % 6 < 3 {
            let entry = &mut self.response[row][column % 6];
            *entry += value;
            if !entry.is_finite() {
                return Err(bad("shell body acceleration response overflow"));
            }
        }
        if row != column && row % 6 < 3 {
            let entry = &mut self.response[column][row % 6];
            *entry += value;
            if !entry.is_finite() {
                return Err(bad("shell body acceleration response overflow"));
            }
        }
        Ok(())
    }
}

fn bad(what: &'static str) -> PlateError {
    PlateError::BadSection { what }
}

#[cfg(test)]
mod tests;
