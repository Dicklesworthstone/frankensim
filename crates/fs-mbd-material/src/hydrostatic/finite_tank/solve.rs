//! One bounded scalar waterline owner for single and shared displacement.
use super::{FiniteHydrostaticTank, HydrostaticError, finite, finite_vec, poll};
use fs_ga::Vec3;
use fs_material::fluid_transport::FluidTransport;

/// One common free surface and one liquid gravitational-storage account.
/// Its displacement is the sum over all admitted, non-overlapping envelopes.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteTankSurface {
    /// Height above the fixed tank floor [m].
    pub surface_height_m: f64,
    /// Corresponding point in the world frame [m].
    pub surface_point_world_m: Vec3,
    /// Total envelope volume below this surface [m3].
    pub displaced_volume_m3: f64,
    /// Signed A*h - total displaced volume - liquid volume residual [m3].
    pub liquid_volume_defect_m3: f64,
    /// Bounded evaluations of the complete displacement configuration.
    pub volume_evaluations: usize,
    /// Liquid gravitational energy relative to the same empty-tank inventory
    /// [J]. This is one shared storage term, not an independent body potential.
    pub liquid_gravitational_energy_change_j: f64,
}

impl<F: FluidTransport + ?Sized> FiniteHydrostaticTank<'_, F> {
    pub(super) fn solve_surface<C: FnMut() -> bool>(
        &self,
        whole_volume_m3: f64,
        cancelled: &mut C,
        mut displacement: impl FnMut(f64, &mut C) -> Result<(f64, f64), HydrostaticError>,
    ) -> Result<FiniteTankSurface, HydrostaticError> {
        let mut evaluations = 0;
        let mut evaluate = |height: f64, cancel: &mut C| {
            poll(cancel)?;
            evaluations += 1;
            let (volume, depth_integral) = displacement(height, cancel)?;
            let residual = finite(self.area_m2 * (height - self.empty_height_m) - volume)?;
            Ok::<_, HydrostaticError>((volume, depth_integral, residual))
        };
        let mut low = self.empty_height_m;
        let (low_volume, low_depth, mut f_low) = evaluate(low, cancelled)?;
        if f_low.abs() <= self.controls.maximum_volume_defect_m3 {
            drop(evaluate);
            return self.finish_surface(low, low_volume, low_depth, f_low, evaluations, cancelled);
        }
        let unconstrained_high = finite(self.empty_height_m + whole_volume_m3 / self.area_m2)?;
        let mut high = unconstrained_high.min(self.geometry.inside_dimensions_m.z);
        let (high_volume, high_depth, mut f_high) = evaluate(high, cancelled)?;
        if f_high.abs() <= self.controls.maximum_volume_defect_m3 {
            drop(evaluate);
            return self.finish_surface(
                high,
                high_volume,
                high_depth,
                f_high,
                evaluations,
                cancelled,
            );
        }
        if f_high < 0.0 {
            drop(evaluate);
            if high == self.geometry.inside_dimensions_m.z {
                return Err(HydrostaticError::TankOverflow {
                    excess_volume_m3: -f_high,
                });
            }
            return Err(HydrostaticError::TankVolumeNotConverged {
                defect_m3: f_high,
                maximum_m3: self.controls.maximum_volume_defect_m3,
                evaluations,
            });
        }
        let mut best_defect = if f_low.abs() < f_high.abs() {
            f_low
        } else {
            f_high
        };
        for evaluation in 3..=self.controls.maximum_volume_evaluations {
            let span = high - low;
            let fraction = -f_low / (f_high - f_low);
            // Unchanged safeguarded secant; periodic bisection bounds progress
            // through small caps and near an endpoint.
            let trial = if evaluation % 4 == 0 || !(0.01..0.99).contains(&fraction) {
                low + 0.5 * span
            } else {
                low + fraction * span
            };
            if trial <= low || trial >= high {
                break;
            }
            let (volume, depth, residual) = evaluate(trial, cancelled)?;
            if residual.abs() <= self.controls.maximum_volume_defect_m3 {
                drop(evaluate);
                return self.finish_surface(trial, volume, depth, residual, evaluations, cancelled);
            }
            if residual.abs() < best_defect.abs() {
                best_defect = residual;
            }
            if residual < 0.0 {
                low = trial;
                f_low = residual;
            } else {
                high = trial;
                f_high = residual;
            }
        }
        Err(HydrostaticError::TankVolumeNotConverged {
            defect_m3: best_defect,
            maximum_m3: self.controls.maximum_volume_defect_m3,
            evaluations,
        })
    }

    fn finish_surface(
        &self,
        height: f64,
        volume: f64,
        depth_integral: f64,
        residual: f64,
        evaluations: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<FiniteTankSurface, HydrostaticError> {
        let geometric_potential = finite(depth_integral - 0.5 * volume * (volume / self.area_m2))?;
        if geometric_potential < 0.0 {
            return Err(HydrostaticError::Invalid(
                "finite tank pressure potential is unrepresentable",
            ));
        }
        let potential =
            finite(self.transport.density_kg_m3() * self.gravity_magnitude * geometric_potential)?;
        let point = finite_vec(
            self.geometry
                .world_from_tank
                .transform_point(Vec3::new(0.0, 0.0, height))?,
        )?;
        poll(cancelled)?;
        Ok(FiniteTankSurface {
            surface_height_m: height,
            surface_point_world_m: point,
            displaced_volume_m3: volume,
            liquid_volume_defect_m3: residual,
            volume_evaluations: evaluations,
            liquid_gravitational_energy_change_j: potential,
        })
    }
}
