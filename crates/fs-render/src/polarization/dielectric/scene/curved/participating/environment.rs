//! Exact environment-tail admission for original homogeneous density grids.
//! The native camera path keeps stochastic flights; explicit source rays use
//! deterministic scattering survival. The two estimators share angular MIS.

use super::*;
use crate::charts::Ray;

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct EnvironmentScatteringSurvival {
    pub log_survival: f64,
    pub cell_visits: u32,
}

enum Tail<'a> {
    Vacuum,
    Infinite,
    Grid(&'a VolumeGrid<'a>, f64),
}

impl<'a> CurvedRayleighMedium<'a> {
    pub(super) fn finite_grid(&self) -> Option<&'a VolumeGrid<'a>> {
        match self.density {
            CurvedDensityProfile::Field { grid, .. } => grid,
            CurvedDensityProfile::Uniform(_) => None,
        }
    }
}

impl CurvedRayleighScene<'_, '_, '_> {
    fn environment_tail(
        &self,
        medium: usize,
        wavelength_nm: f64,
    ) -> Result<Option<Tail<'_>>, CurvedSceneError> {
        let host = self.scene.media.get(medium).ok_or(SceneError::InvalidInitialMedium)?;
        // Only the native straight homogeneous owner can prove this infinite
        // source leg. An equal local index never supplies a GRIN escape proof.
        if !matches!(host, SceneMedium::Homogeneous(_)) {
            return Ok(None);
        }
        let Some(law) = self.media.get(medium) else {
            return Ok(Some(Tail::Vacuum));
        };
        if law.scattering_per_m(wavelength_nm)? == 0.0 {
            return Ok(Some(Tail::Vacuum));
        }
        let CurvedRayleighBinding::Density(law) = law else {
            return Ok(None);
        };
        Ok(match law.density {
            CurvedDensityProfile::Uniform(_) => Some(Tail::Infinite),
            CurvedDensityProfile::Field { grid: Some(grid), .. } => {
                Some(Tail::Grid(grid, law.scattering.at(wavelength_nm).map_err(volume_error)?))
            }
            CurvedDensityProfile::Field { grid: None, .. } => None,
        })
    }

    /// Same admission as the explicit source leg, without traversing cells.
    /// The actual outgoing medium owns complementary camera MIS.
    pub(crate) fn supports_environment_tail(
        &self,
        medium: usize,
        wavelength_nm: f64,
    ) -> Result<bool, CurvedSceneError> {
        Ok(self.environment_tail(medium, wavelength_nm)?.is_some())
    }

    /// None denotes an unsupported active field, never unit survival.
    /// Native finite grids use exact original density cells and face ownership.
    pub(crate) fn environment_scattering_survival(
        &self,
        cx: &Cx<'_>,
        medium: usize,
        ray: Ray,
        wavelength_nm: f64,
        max_cells: u32,
    ) -> Result<Option<EnvironmentScatteringSurvival>, CurvedSceneError> {
        cx.checkpoint()?;
        let result = match self.environment_tail(medium, wavelength_nm)? {
            None => return Ok(None),
            Some(Tail::Vacuum) => EnvironmentScatteringSurvival {
                log_survival: 0.0,
                cell_visits: 0,
            },
            Some(Tail::Infinite) => EnvironmentScatteringSurvival {
                log_survival: f64::NEG_INFINITY,
                cell_visits: 0,
            },
            Some(Tail::Grid(grid, coefficient)) => {
                let (depth, cell_visits) = crate::volumes::transport::grid_scattering_optical_depth(
                    cx, grid, ray.origin, ray.dir, coefficient, max_cells,
                ).map_err(volume_error)?;
                EnvironmentScatteringSurvival { log_survival: -depth, cell_visits }
            }
        };
        cx.checkpoint()?;
        Ok(Some(result))
    }
}

/// Place an original majorant flight inside its exact finite support. No
/// collision after the far face means a surviving camera sample, not an
/// additional Beer weight. The caller's native geometry can interrupt it.
pub(super) fn supported_distance(
    (entry, length): (f64, f64),
    local_distance: f64,
) -> Result<Option<f64>, CurvedSceneError> {
    if !entry.is_finite() || entry < 0.0 || !length.is_finite() || length <= 0.0
        || local_distance.is_nan() || local_distance < 0.0
    {
        return Err(volume_error(VolumeTransportError::Unrepresentable));
    }
    if local_distance >= length {
        return Ok(None);
    }
    let absolute = entry + local_distance;
    if !absolute.is_finite() || (local_distance > 0.0 && absolute <= entry) {
        return Err(volume_error(VolumeTransportError::Unrepresentable));
    }
    Ok(Some(absolute))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn g1_finite_support_retains_the_vacuum_prefix_and_truncates_only_after_exit() {
        assert_eq!(supported_distance((7.0, 3.0), 0.0).unwrap(), Some(7.0));
        assert_eq!(supported_distance((7.0, 3.0), 2.0).unwrap(), Some(9.0));
        assert_eq!(supported_distance((7.0, 3.0), 3.0).unwrap(), None);
        assert_eq!(supported_distance((7.0, 3.0), f64::INFINITY).unwrap(), None);
        assert_eq!(supported_distance((0.0, 3.0), 0.25).unwrap(), Some(0.25));
        assert!(supported_distance((1e30, 3.0), 1.0).is_err());
        assert!(supported_distance((0.0, 3.0), f64::NAN).is_err());
    }
}
