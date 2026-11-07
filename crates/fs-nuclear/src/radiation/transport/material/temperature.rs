//! Temperature-dependent original photon mass-interaction data.
//!
//! A source supplies the ENTIRE energy kernel at every temperature: total
//! collisions, locally thermalized heat and outgoing photon energy. Convex
//! coefficient interpolation preserves their physical separation and closure.
//! No opacity law, electron stopping power, phase mixture or material datum is
//! inferred from temperature or from a mass attenuation coefficient alone.
//! Density and geometry remain the responsibility of the original specimen.

use super::*;

/// One supplied complete photon mass kernel at a declared temperature.
#[derive(Clone, Debug)]
pub struct PhotonMassTemperatureKnot {
    /// Absolute temperature [K], finite and strictly positive.
    pub temperature_k: f64,
    /// Original immutable source, including its own data and material identities.
    pub material: PhotonMassMaterial,
}

/// Immutable piecewise-affine mass coefficients over a CLOSED kelvin domain.
/// Knots share exactly one material, photon kind and energy grid. Source IDs at
/// individual knots remain available; the interpolation is an authored model,
/// not a new measured dataset or a claim that temperature is a complete state.
#[derive(Clone, Debug)]
pub struct TemperaturePhotonMassMaterial(Arc<TemperatureData>);

#[derive(Debug)]
struct TemperatureData {
    data_id: String,
    knots: Vec<PhotonMassTemperatureKnot>,
    limits: PhotonSlabLimits,
    values: usize,
}

impl TemperaturePhotonMassMaterial {
    /// Retain at least two strictly ordered samples without extrapolation or
    /// coefficient renormalization. Limits cover all retained coefficient
    /// scalars, grids and temperature values, not merely one knot.
    pub fn new(
        data_id: String,
        knots: Vec<PhotonMassTemperatureKnot>,
        limits: PhotonSlabLimits,
        mut checkpoint: impl FnMut() -> Result<(), RadiationError>,
    ) -> Result<Self, RadiationError> {
        checkpoint()?;
        if data_id.trim().is_empty() || knots.len() < 2 {
            return Err(RadiationError::InvalidInput(
                "identified photon temperature data and at least two knots required",
            ));
        }
        let first = &knots[0].material;
        let groups = first.groups().len();
        let values = groups
            .checked_add(2)
            .and_then(|width| groups.checked_mul(width))
            .and_then(|n| n.checked_add(first.group_edges_j().len()))
            .and_then(|n| n.checked_add(1))
            .and_then(|n| n.checked_mul(knots.len()))
            .filter(|&n| {
                groups <= limits.max_groups && n <= limits.max_values && n <= limits.max_work
            })
            .ok_or(RadiationError::Budget(
                "photon temperature table storage/work",
            ))?;
        let mut previous = 0.0;
        for knot in &knots {
            checkpoint()?;
            if !knot.temperature_k.is_finite()
                || knot.temperature_k <= previous
                || knot.material.material_id() != first.material_id()
                || knot.material.radiation() != first.radiation()
                || knot.material.group_edges_j() != first.group_edges_j()
            {
                return Err(RadiationError::InvalidInput(
                    "photon knots need ordered positive kelvin and one original material/grid",
                ));
            }
            previous = knot.temperature_k;
        }
        checkpoint()?;
        Ok(Self(Arc::new(TemperatureData {
            data_id,
            knots,
            limits,
            values,
        })))
    }

    /// Authored identity of this interpolation model and its selected sources.
    pub fn data_id(&self) -> &str {
        &self.0.data_id
    }
    /// Original material shared by every temperature sample.
    pub fn material_id(&self) -> u64 {
        self.0.knots[0].material.material_id()
    }
    /// Original gamma or X-ray kind; it never changes with temperature.
    pub fn radiation(&self) -> DecayRadiation {
        self.0.knots[0].material.radiation()
    }
    /// Original shared photon energy grid [J].
    pub fn group_edges_j(&self) -> &[f64] {
        self.0.knots[0].material.group_edges_j()
    }
    /// Original immutable samples, including their individual source identities.
    pub fn knots(&self) -> &[PhotonMassTemperatureKnot] {
        &self.0.knots
    }
    /// Closed temperature support [K]; interpolation never extrapolates.
    pub fn temperature_bounds_k(&self) -> (f64, f64) {
        (
            self.0.knots[0].temperature_k,
            self.0.knots[self.0.knots.len() - 1].temperature_k,
        )
    }
    /// Retained scalar count used by aggregate consumer preparation budgets.
    pub fn value_count(&self) -> usize {
        self.0.values
    }

    /// Evaluate the original mass-kernel owner at one supported temperature.
    /// Exact knots return that original immutable material unchanged. Interior
    /// samples carry this table's identity and retain explicit interpolation.
    pub fn sample(
        &self,
        temperature_k: f64,
        mut checkpoint: impl FnMut() -> Result<(), RadiationError>,
    ) -> Result<PhotonMassMaterial, RadiationError> {
        checkpoint()?;
        let (minimum, maximum) = self.temperature_bounds_k();
        if !temperature_k.is_finite() || temperature_k < minimum || temperature_k > maximum {
            return Err(RadiationError::InvalidInput(
                "temperature is outside the original photon mass-data domain",
            ));
        }
        let right = self
            .0
            .knots
            .partition_point(|knot| knot.temperature_k < temperature_k);
        let upper = &self.0.knots[right];
        if upper.temperature_k == temperature_k {
            return Ok(upper.material.clone());
        }
        let lower = &self.0.knots[right - 1];
        let coefficient = |a, b| {
            interpolate(
                a,
                b,
                lower.temperature_k,
                upper.temperature_k,
                temperature_k,
            )
            .ok_or(RadiationError::Arithmetic)
        };
        let mut groups = Vec::with_capacity(lower.material.groups().len());
        for (a, b) in lower.material.groups().iter().zip(upper.material.groups()) {
            checkpoint()?;
            let mut redistribution = Vec::with_capacity(a.redistribution_m2_kg.len());
            for (index, (&a, &b)) in a
                .redistribution_m2_kg
                .iter()
                .zip(&b.redistribution_m2_kg)
                .enumerate()
            {
                if index % 64 == 0 {
                    checkpoint()?;
                }
                redistribution.push(coefficient(a, b)?);
            }
            groups.push(PhotonMassCollision {
                total_m2_kg: coefficient(a.total_m2_kg, b.total_m2_kg)?,
                heat_m2_kg: coefficient(a.heat_m2_kg, b.heat_m2_kg)?,
                redistribution_m2_kg: redistribution,
            });
        }
        // Reuse the original closure/downscatter admission after rounding.
        PhotonMassMaterial::new(
            self.0.data_id.clone(),
            self.material_id(),
            self.radiation(),
            self.group_edges_j().to_vec(),
            groups,
            self.0.limits,
            checkpoint,
        )
    }
}

// All arguments have been admitted: nonnegative finite coefficients and a
// positive ordered kelvin interval containing t. Evaluation from the nearer
// endpoint retains a weak positive branch near a zero-valued endpoint.
fn interpolate(a: f64, b: f64, left: f64, right: f64, t: f64) -> Option<f64> {
    if t == left || a == b {
        return Some(a);
    }
    if t == right {
        return Some(b);
    }
    let value = if t - left <= right - t {
        (b - a).mul_add((t - left) / (right - left), a)
    } else {
        (a - b).mul_add((right - t) / (right - left), b)
    };
    if value.is_finite() && value >= 0.0 && (value > 0.0 || (a == 0.0 && b == 0.0)) {
        Some(value)
    } else {
        None
    }
}

#[cfg(test)]
mod tests;
