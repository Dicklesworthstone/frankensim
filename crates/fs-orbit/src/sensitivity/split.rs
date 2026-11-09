//! Shared projection of local source and frequency-port parameter directions.

use super::{
    ForcedOrbitSensitivity, ForcedParameterPoint, ForcedSensitivityDirection,
    ForcedSensitivityError, PreparedForcedSensitivity, Result, admit_direction, invalid, poll,
    project_direction,
};

impl PreparedForcedSensitivity {
    /// Differentiate a complete source whose residual is a linear frequency
    /// port minus a time-local island. Both callbacks hold Fourier coefficients
    /// FIXED. `island_direction` includes the total fixed-PHASE parameter
    /// partial: for an explicit f(t,p), its clock term is
    /// `-(t/Omega)*Omega_p*partial_t f`. `port_direction` writes EVERY packed
    /// row of the remaining direct frequency-port derivative, including clock,
    /// descriptor and coefficient changes. It excludes the island partial.
    ///
    /// This combines the original cached waveform and centered AFT projection
    /// with the authored port column, without evaluating the baseline source
    /// or tangent again. The prepared factor solves the full combined direction and
    /// checks both grids. Unlike the island-only entry, a changing clock is
    /// admitted because both pieces are explicitly supplied.
    pub fn solve_split_direction(
        &mut self,
        direction: ForcedSensitivityDirection,
        mut island_direction: impl FnMut(f64, &[f64], &mut [f64]) -> bool,
        mut port_direction: impl FnMut(ForcedParameterPoint<'_>, &mut [f64]) -> bool,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ForcedOrbitSensitivity> {
        admit_direction(direction)?;
        let [coarse, fine] = self.project_split_direction(
            direction,
            &mut island_direction,
            &mut port_direction,
            &mut cancelled,
        )?;
        self.solve_columns(direction, 0.0, &coarse, &fine, &mut cancelled)
    }

    pub(super) fn project_split_direction(
        &mut self,
        direction: ForcedSensitivityDirection,
        island_direction: &mut impl FnMut(f64, &[f64], &mut [f64]) -> bool,
        port_direction: &mut impl FnMut(ForcedParameterPoint<'_>, &mut [f64]) -> bool,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<[Vec<f64>; 2]> {
        poll(cancelled)?;
        // Three linear columns fit the existing 16*n numerical scratch bound.
        // Source callback storage/work retain the original explicit exclusion.
        let mut linear = vec![f64::NAN; self.pack.coeff_len()];
        let mut columns = [Vec::new(), Vec::new()];
        for (column, aft) in columns.iter_mut().zip([&mut self.coarse, &mut self.fine]) {
            *column = project_direction(&self.pack, self.omega, aft, island_direction, cancelled)?;
            poll(cancelled)?;
            linear.fill(f64::NAN);
            let point = ForcedParameterPoint {
                parameter: direction.parameter,
                angular_frequency: self.omega,
                harmonics: self.pack.n_harm,
                aft_samples: aft.samples,
                coefficients: &self.coefficients[..self.pack.coeff_len()],
            };
            if !port_direction(point, &mut linear) || linear.iter().any(|x| !x.is_finite()) {
                return Err(ForcedSensitivityError::DirectionRefused {
                    samples: aft.samples,
                });
            }
            for (value, direct) in column.iter_mut().zip(&linear) {
                poll(cancelled)?;
                *value += direct;
                if !value.is_finite() {
                    return Err(invalid(
                        "combined physical residual direction must be finite",
                    ));
                }
            }
        }
        poll(cancelled)?;
        Ok(columns)
    }
}
