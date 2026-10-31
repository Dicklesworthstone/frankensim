//! A uniform reversible core in series with the common moving gap.
//! Source knots supply both the exact scalar inversion and the energy defect;
//! there is no second fitted saturation curve or iterative magnetic solver.

use std::sync::Arc;

use fs_material::field_storage::nonlinear::{
    ReversibleMagneticMaterial, ReversibleMagneticResponse,
};

use super::{MaterialCircuitError, invalid, positive};

/// One source-resolved, uniform reversible core carrying the common gap flux.
///
/// Area and path length are fixed physical geometry. All material context stays
/// pinned; remanence, hysteresis, eddy loss, fringing and thermal evolution are
/// outside this owner. Its stored energy is `area * length * integral H dB`,
/// never `B * H / 2` after a saturation knee.
#[derive(Clone, Debug, PartialEq)]
pub struct ReversibleMagneticCore {
    material: Arc<ReversibleMagneticMaterial>,
    area_m2: f64,
    path_length_m: f64,
    volume_m3: f64,
    maximum_flux_wb: f64,
}

impl ReversibleMagneticCore {
    /// Bind the original reversible source path and positive uniform geometry.
    /// A nonzero explicit thermal scale derivative needs a field-entropy owner
    /// and cannot be silently frozen into this purely electromechanical store.
    pub fn new(
        material: ReversibleMagneticMaterial,
        area_m2: f64,
        path_length_m: f64,
    ) -> Result<Self, MaterialCircuitError> {
        let volume_m3 = area_m2 * path_length_m;
        let maximum_flux_wb = area_m2 * material.maximum_flux_density_t();
        if [area_m2, path_length_m, volume_m3, maximum_flux_wb]
            .iter()
            .any(|&value| !positive(value) || !positive(1.0 / value))
            || !positive(path_length_m / area_m2)
            || !positive(path_length_m * material.field_scale())
            || material
                .field_scale_slope_per_k()
                .is_some_and(|slope| slope != 0.0)
        {
            return Err(invalid(
                "shared-gap core requires representable positive geometry and no unowned field entropy",
            ));
        }
        let core = Self {
            material: Arc::new(material),
            area_m2,
            path_length_m,
            volume_m3,
            maximum_flux_wb,
        };
        let response = core
            .material
            .evaluate(core.material.maximum_flux_density_t())?;
        if !positive(core.volume_m3 * response.energy_density_j_m3)
            || !positive(core.path_length_m * response.field_strength_a_m)
            || core.material.source_path().knots.windows(2).any(|pair| {
                let slope = (pair[1].1.evidence.value.value - pair[0].1.evidence.value.value)
                    / (pair[1].0 - pair[0].0)
                    * core.material.field_scale();
                !positive(core.path_length_m * slope)
                    || !positive((core.path_length_m / core.area_m2) * slope)
            })
        {
            return Err(invalid(
                "shared-gap core energy or reluctance is unrepresentable",
            ));
        }
        Ok(core)
    }

    /// Original material, complete source path, fixed context and receipts.
    pub fn material(&self) -> &ReversibleMagneticMaterial {
        &self.material
    }
    /// Fixed core cross section [m²]; it may differ from the gap pole area.
    pub const fn area_m2(&self) -> f64 {
        self.area_m2
    }
    /// Fixed magnetic path length through the core [m].
    pub const fn path_length_m(&self) -> f64 {
        self.path_length_m
    }
    /// Largest admitted magnitude of the common physical flux [Wb].
    pub const fn maximum_flux_wb(&self) -> f64 {
        self.maximum_flux_wb
    }

    pub(super) fn admit_reluctance(&self, reluctance: f64) -> Result<(), MaterialCircuitError> {
        if !positive(reluctance * self.area_m2)
            || !positive(
                reluctance * self.maximum_flux_wb
                    + self.path_length_m * self.material.maximum_field_strength_a_m(),
            )
            || self.material.source_path().knots.windows(2).any(|pair| {
                let slope = (pair[1].1.evidence.value.value - pair[0].1.evidence.value.value)
                    / (pair[1].0 - pair[0].0)
                    * self.material.field_scale();
                let derivative = reluctance * self.area_m2 + self.path_length_m * slope;
                !positive(derivative) || !positive(1.0 / derivative)
            })
        {
            return Err(invalid("shared-gap core inversion is unrepresentable"));
        }
        Ok(())
    }

    // q*Phi + length*H(Phi/area) = drive. Both terms are strictly increasing;
    // locate their source segment in O(log knots), then invert that affine span.
    // Out-of-source trials return NaN to the original native nonlinear solver.
    pub(super) fn flux(&self, reluctance: f64, drive: f64) -> f64 {
        if !positive(reluctance) || !drive.is_finite() {
            return f64::NAN;
        }
        let knots = &self.material.source_path().knots;
        let linear_b = reluctance * self.area_m2;
        let core_scale = self.path_length_m * self.material.field_scale();
        let at = |index: usize| {
            linear_b * knots[index].0 + core_scale * knots[index].1.evidence.value.value
        };
        let magnitude = drive.abs();
        let maximum = at(knots.len() - 1);
        if !maximum.is_finite() || magnitude > maximum {
            return f64::NAN;
        }
        if magnitude == maximum {
            return self.maximum_flux_wb.copysign(drive);
        }
        let segment = knots
            .partition_point(|(b, answer)| {
                linear_b * b + core_scale * answer.evidence.value.value <= magnitude
            })
            .saturating_sub(1)
            .min(knots.len() - 2);
        let (b0, b1) = (knots[segment].0, knots[segment + 1].0);
        let slope = (knots[segment + 1].1.evidence.value.value
            - knots[segment].1.evidence.value.value)
            / (b1 - b0);
        let b = b0 + (magnitude - at(segment)) / (linear_b + core_scale * slope);
        // Only roundoff in an already admitted segment is bounded here. No
        // drive beyond the original source range is projected onto its boundary.
        let flux = self.area_m2 * b.clamp(b0, b1);
        if !flux.is_finite() || magnitude > 0.0 && flux <= 0.0 {
            return f64::NAN;
        }
        flux.copysign(drive)
    }

    fn density(&self, flux: f64) -> Option<f64> {
        if !flux.is_finite() || flux.abs() > self.maximum_flux_wb {
            return None;
        }
        // A supported endpoint Phi = A*Bmax may round up on division by A.
        Some((flux.abs() / self.area_m2).min(self.material.maximum_flux_density_t()))
    }

    fn response(&self, flux: f64) -> Option<ReversibleMagneticResponse> {
        self.material.evaluate(self.density(flux)?).ok()
    }

    pub(super) fn energy(&self, flux: f64) -> f64 {
        self.response(flux)
            .map_or(f64::NAN, |r| self.volume_m3 * r.energy_density_j_m3)
    }

    pub(super) fn secant_reluctance(&self, flux: f64) -> f64 {
        self.response(flux).map_or(f64::NAN, |r| {
            (self.path_length_m / self.area_m2) * r.secant_reluctivity_m_h
        })
    }

    pub(super) fn differential_reluctance(&self, flux: f64) -> f64 {
        self.response(flux).map_or(f64::NAN, |r| {
            (self.path_length_m / self.area_m2) * r.differential_reluctivity_m_h
        })
    }

    // W(b)-W(a)-avg(W'(a),W'(b))*(b-a), evaluated from source slope jumps
    // instead of subtracting large endpoint energies. A single affine H span
    // has exactly zero defect; only crossed saturation knots contribute.
    pub(super) fn trapezoid_remainder(&self, before_flux: f64, after_flux: f64) -> f64 {
        let (Some(before), Some(after)) = (self.density(before_flux), self.density(after_flux))
        else {
            return f64::NAN;
        };
        let a = before.copysign(before_flux);
        let b = after.copysign(after_flux);
        if a == b || a == -b {
            return 0.0;
        }
        let (lower, upper) = (a.min(b), a.max(b));
        let knots = &self.material.source_path().knots;
        let mut defect = 0.0;
        for sign in [-1.0, 1.0] {
            let (start, end) = if sign > 0.0 {
                (lower.max(0.0), upper)
            } else {
                ((-upper).max(0.0), -lower)
            };
            if end <= start {
                continue;
            }
            let first = knots.partition_point(|(k, _)| *k <= start).max(1);
            let last = knots
                .partition_point(|(k, _)| *k < end)
                .min(knots.len() - 1);
            for index in first..last {
                let (k0, k, k1) = (knots[index - 1].0, knots[index].0, knots[index + 1].0);
                let (h0, h, h1) = (
                    knots[index - 1].1.evidence.value.value,
                    knots[index].1.evidence.value.value,
                    knots[index + 1].1.evidence.value.value,
                );
                let jump = sign
                    * ((h1 - h) / (k1 - k) - (h - h0) / (k - k0))
                    * self.material.field_scale();
                let location = sign * k;
                defect -= 0.5 * jump * (location - lower) * (upper - location);
            }
        }
        (b - a).signum() * self.volume_m3 * defect
    }
}
