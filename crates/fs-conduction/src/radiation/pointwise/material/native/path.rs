//! Whole-chart support and feedback preparation for the original gray source.
use super::{ConductionError, MaterialGraySurfaceLaw, feedback, invalid};
use crate::radiation::STEFAN_BOLTZMANN_W_M2_K4;
use fs_blake3::{ContentHash, DomainHasher};
use fs_matdb::PropertyPathSlice;

const MAX_GRAY_SURFACE_PATH_KNOTS: usize = 65_536;

#[derive(Debug)]
struct Interval {
    temperature_bounds_k: [f64; 2],
    emissivity_enthalpy_slope: f64,
}

/// Original emissivity paths and sufficient positive inward-power feedback.
/// This is immutable admission data, not a replacement emissivity evaluator.
/// The original `MaterialGraySurfaceLaw::evaluate` remains the runtime owner.
#[derive(Debug)]
pub struct MaterialGraySurfacePath {
    paths: Vec<PropertyPathSlice>,
    intervals: Vec<Interval>,
    maximum_emissivity: f64,
    identity: ContentHash,
}

impl MaterialGraySurfaceLaw<'_> {
    /// Admit uninterrupted support across the COMPLETE original caloric chart.
    /// Every source knot/cell crossing survives; no box support or phase mixture
    /// is inferred. Work is capped and polls cancellation at source boundaries.
    pub fn prepare_caloric_path(
        &self,
        mut checkpoint: impl FnMut() -> Result<(), ConductionError>,
    ) -> Result<MaterialGraySurfacePath, ConductionError> {
        checkpoint()?;
        let curve = self.caloric_curve();
        let knots = curve.knots();
        if self.bounds_j_kg
            != [
                knots[0].specific_enthalpy_j_kg,
                knots[knots.len() - 1].specific_enthalpy_j_kg,
            ]
        {
            return Err(invalid(
                "gray path bounds must cover the complete original caloric chart",
            ));
        }
        let material = self.material_point();
        let claim = material
            .material
            .claims()
            .claim(self.source.claim())
            .ok_or_else(|| invalid("original gray surface claim is absent"))?;
        let mut paths = Vec::new();
        let mut intervals = Vec::new();
        let mut maximum_emissivity: f64 = 0.0;
        let mut count = 0usize;
        let mut identity =
            DomainHasher::new("org.frankensim.fs-conduction.material-gray-surface-path.v1");
        identity.update(curve.identity().as_bytes());
        identity.update(
            self.reference_emissivity()
                .receipt()
                .content_hash()
                .as_bytes(),
        );
        for (index, pair) in knots.windows(2).enumerate() {
            checkpoint()?;
            if pair[0].liquid_mass_fraction != pair[1].liquid_mass_fraction
                && material.liquid_fraction_axis.is_none()
            {
                return Err(invalid(
                    "phase radiation needs its explicitly evolving mass fraction coordinate",
                ));
            }
            let point = |h| {
                curve
                    .state_at_specific_enthalpy(h)
                    .map_err(|e| invalid(&e.to_string()))
                    .and_then(|state| material.at(state).map_err(|e| invalid(&e.to_string())))
            };
            let start = point(pair[0].specific_enthalpy_j_kg)?;
            let end = point(pair[1].specific_enthalpy_j_kg)?;
            let path = material
                .material
                .claims()
                .query_path_pinned_typed(&claim.key, &start, &end, self.source.claim())
                .map_err(|e| invalid(&e.to_string()))?;
            count = count
                .checked_add(path.knots.len())
                .filter(|count| *count <= MAX_GRAY_SURFACE_PATH_KNOTS)
                .ok_or_else(|| invalid("gray source path exceeds its retained knot budget"))?;
            identity.update(&(path.knots.len() as u64).to_le_bytes());
            for (offset, (parameter, answer)) in path.knots.iter().enumerate() {
                if offset % 256 == 0 {
                    checkpoint()?;
                }
                let epsilon = answer.evidence.value.value;
                if !epsilon.is_finite() || epsilon <= 0.0 || epsilon > 1.0 {
                    return Err(invalid(
                        "original gray path emissivity must remain in (0,1]",
                    ));
                }
                maximum_emissivity = maximum_emissivity.max(epsilon);
                identity.update(&parameter.to_bits().to_le_bytes());
                identity.update(answer.receipt.content_hash().as_bytes());
            }
            let sensible = curve
                .sensible_heat_capacities_on_interval(index)
                .map(|_| {
                    curve
                        .sensible_heat_capacity_bounds_j_kg_k(index)
                        .map(|cp| (pair[1].temperature_k - pair[0].temperature_k, cp[0]))
                        .map_err(|e| invalid(&e.to_string()))
                })
                .transpose()?;
            for (offset, pair_of_sources) in path.knots.windows(2).enumerate() {
                if offset % 256 == 0 {
                    checkpoint()?;
                }
                let slope = feedback::emissivity_enthalpy_slope(
                    [
                        pair_of_sources[0].1.evidence.value.value,
                        pair_of_sources[1].1.evidence.value.value,
                    ],
                    pair_of_sources[1].0 - pair_of_sources[0].0,
                    pair[1].specific_enthalpy_j_kg - pair[0].specific_enthalpy_j_kg,
                    sensible,
                )
                .ok_or_else(|| invalid("gray emissivity feedback slope is unrepresentable"))?;
                let temperature = |side: usize| {
                    pair_of_sources[side]
                        .1
                        .receipt
                        .query_point
                        .iter()
                        .find(|(axis, _)| axis == material.temperature_axis)
                        .map(|(_, value)| *value)
                        .ok_or_else(|| invalid("gray source path lost its temperature coordinate"))
                };
                intervals.push(Interval {
                    temperature_bounds_k: [temperature(0)?, temperature(1)?],
                    emissivity_enthalpy_slope: slope,
                });
            }
            paths.push(path);
        }
        checkpoint()?;
        Ok(MaterialGraySurfacePath {
            paths,
            intervals,
            maximum_emissivity,
            identity: identity.finalize(),
        })
    }
}

impl MaterialGraySurfacePath {
    pub fn source_paths(&self) -> &[PropertyPathSlice] {
        &self.paths
    }
    pub const fn maximum_emissivity(&self) -> f64 {
        self.maximum_emissivity
    }
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Maximum positive inward radiative dP/dh [W/(J/kg)]. Constant convection
    /// and the -4*epsilon*T^3*dT/dh emission term are nonpositive and need no
    /// additional bound. Latent epsilon changes remain active when dT/dh=0.
    /// This bounds the declared model, not source uncertainty or time error.
    pub fn maximum_positive_radiative_feedback_w_per_j_kg(
        &self,
        area_m2: f64,
        reservoir_k: f64,
        mut checkpoint: impl FnMut() -> Result<(), ConductionError>,
    ) -> Result<f64, ConductionError> {
        checkpoint()?;
        let mut bound: f64 = 0.0;
        for (index, interval) in self.intervals.iter().enumerate() {
            if index % 256 == 0 {
                checkpoint()?;
            }
            bound = bound.max(
                feedback::positive_radiative_feedback(
                    interval.emissivity_enthalpy_slope,
                    interval.temperature_bounds_k,
                    area_m2,
                    reservoir_k,
                    STEFAN_BOLTZMANN_W_M2_K4,
                )
                .ok_or_else(|| invalid("gray radiative feedback bound or geometry is invalid"))?,
            );
        }
        checkpoint()?;
        Ok(bound)
    }
}
