//! Temperature-dependent original photon mass-interaction data.
//!
//! A source supplies the ENTIRE energy kernel at every temperature: total
//! collisions, locally thermalized heat and outgoing photon energy. Convex
//! coefficient interpolation preserves their physical separation and closure.
//! No opacity law, electron stopping power, phase mixture or material datum is
//! inferred from temperature or from a mass attenuation coefficient alone.
//! Density and geometry remain the responsibility of the original specimen.

use std::cmp::Reverse;
use std::collections::BinaryHeap;

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

/// One original temperature law at a constant fraction of mixture MASS.
#[derive(Clone, Debug)]
pub struct TemperaturePhotonMassConstituent {
    /// Fraction of total mixture mass, finite in [0, 1]; all fractions close one.
    pub mass_fraction: f64,
    /// Original immutable temperature law, source identities and kelvin support.
    pub material: TemperaturePhotonMassMaterial,
}

/// Original temperature laws and the physical basis of their homogenization.
#[derive(Clone, Debug)]
pub struct TemperaturePhotonMassMixture {
    /// Applicability of each complete constituent kernel IN the target mixture
    /// throughout the common temperature domain, including local deposition
    /// and a common within-group spectral weighting.
    pub kernel_applicability: String,
    /// Original authored order, constant mass fractions and immutable sources.
    pub constituents: Vec<TemperaturePhotonMassConstituent>,
}

#[derive(Debug)]
struct TemperatureData {
    data_id: String,
    knots: Vec<PhotonMassTemperatureKnot>,
    limits: PhotonSlabLimits,
    values: usize,
    mixture: Option<TemperaturePhotonMassMixture>,
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
            mixture: None,
        })))
    }

    /// Prepare a homogeneous, constant-composition mixture of temperature laws.
    ///
    /// Every complete constituent kernel is sampled at the SAME temperature.
    /// The union of source knots on their common closed domain preserves every
    /// piecewise-affine change of slope; mixing only the endpoints would erase
    /// intermediate opacity structure. Preparation is done once, so subsequent
    /// sampling has the same cost as a directly supplied temperature table.
    /// All sources, including zero-weight ones, must share the photon kind and
    /// exact energy grid and have a common domain of positive width.
    ///
    /// Fractions and physical applicability have the same meaning as
    /// [`PhotonMassMaterial::from_mass_fractions`]. The caller must establish
    /// complete heat/redistribution kernels applicable IN this target mixture
    /// throughout the supported domain. Attenuation additivity alone does not
    /// supply them. No extrapolation, changing phase/composition, density law
    /// or measured-data qualification is inferred.
    ///
    /// Limits include original source table values, the conservatively counted
    /// sampled constituent kernels retained at every mixed knot, and preparation
    /// work. Shared source objects are never mutated or renormalized.
    pub fn from_mass_fractions(
        data_id: String,
        material_id: u64,
        constituents: Vec<TemperaturePhotonMassConstituent>,
        kernel_applicability: String,
        limits: PhotonSlabLimits,
        mut checkpoint: impl FnMut() -> Result<(), RadiationError>,
    ) -> Result<Self, RadiationError> {
        checkpoint()?;
        if data_id.trim().is_empty()
            || kernel_applicability.trim().is_empty()
            || constituents.is_empty()
        {
            return Err(RadiationError::InvalidInput(
                "identified photon temperature mixture, mass fractions and kernel applicability required",
            ));
        }
        let budget = || RadiationError::Budget("photon temperature mixture storage/work");
        let first = &constituents[0].material;
        let groups = first.group_edges_j().len() - 1;
        let knot_values = groups
            .checked_add(2)
            .and_then(|width| groups.checked_mul(width))
            .and_then(|n| n.checked_add(first.group_edges_j().len()))
            .and_then(|n| n.checked_mul(constituents.len().checked_add(1)?))
            .and_then(|n| n.checked_add(constituents.len()))
            // Kelvin value and the original mass-mixture reduction allowance.
            .and_then(|n| n.checked_add(3))
            .filter(|_| groups <= limits.max_groups)
            .ok_or_else(budget)?;
        let mut source_values = 0_usize;
        let mut source_knots = 0_usize;
        let mut minimum = 0.0_f64;
        let mut maximum = f64::MAX;
        let mut fractions = Sum::default();
        for constituent in &constituents {
            checkpoint()?;
            if !nonnegative(constituent.mass_fraction)
                || constituent.mass_fraction > 1.0
                || constituent.material.radiation() != first.radiation()
                || constituent.material.group_edges_j() != first.group_edges_j()
            {
                return Err(RadiationError::InvalidInput(
                    "photon temperature mixture needs fractions in [0,1] and one photon kind/group grid",
                ));
            }
            fractions.add(constituent.mass_fraction)?;
            source_values = source_values
                .checked_add(constituent.material.value_count())
                .and_then(|n| n.checked_add(1))
                .filter(|&n| n <= limits.max_values && n <= limits.max_work)
                .ok_or_else(budget)?;
            source_knots = source_knots
                .checked_add(constituent.material.knots().len())
                .ok_or_else(budget)?;
            let (lower, upper) = constituent.material.temperature_bounds_k();
            minimum = minimum.max(lower);
            maximum = maximum.min(upper);
        }
        fractions.add(-1.0)?;
        if fractions.hi.abs() + fractions.lo.abs() > ROW_CLOSURE_TOLERANCE {
            return Err(RadiationError::InvalidInput(
                "photon constituent mass fractions must close one without normalization",
            ));
        }
        if minimum >= maximum {
            return Err(RadiationError::InvalidInput(
                "photon temperature constituents need a common domain of positive width",
            ));
        }
        // Merge already ordered sources in O(N log C), with a poll at each
        // source knot. Charge the heap comparison bound and three kernel passes
        // (sampling, mixing, closure admission) before coefficient allocation.
        let merge_work = source_knots
            .checked_mul(2 * constituents.len().ilog2() as usize + 3)
            .and_then(|n| n.checked_add(source_values))
            .ok_or_else(budget)?;
        let admit_knots = |count: usize| {
            let prepared = count.checked_mul(knot_values).ok_or_else(budget)?;
            let values = source_values.checked_add(prepared).ok_or_else(budget)?;
            let work = prepared
                .checked_mul(3)
                .and_then(|n| n.checked_add(merge_work))
                .ok_or_else(budget)?;
            if values > limits.max_values || work > limits.max_work {
                return Err(budget());
            }
            Ok(values)
        };
        // The two required endpoints also cover the temporary merge heap.
        admit_knots(2)?;
        let mut pending = BinaryHeap::with_capacity(constituents.len());
        for (source, constituent) in constituents.iter().enumerate() {
            checkpoint()?;
            let knots = constituent.material.knots();
            let next = knots.partition_point(|knot| knot.temperature_k <= minimum);
            if next < knots.len() && knots[next].temperature_k < maximum {
                // Positive finite kelvin has the same order as its IEEE bits.
                pending.push(Reverse((knots[next].temperature_k.to_bits(), source, next)));
            }
        }
        let mut temperatures = vec![minimum];
        while let Some(Reverse((bits, source, index))) = pending.pop() {
            checkpoint()?;
            let temperature = f64::from_bits(bits);
            if temperatures.last().copied() != Some(temperature) {
                admit_knots(temperatures.len().checked_add(2).ok_or_else(budget)?)?;
                temperatures.push(temperature);
            }
            let knots = constituents[source].material.knots();
            let next = index + 1;
            if next < knots.len() && knots[next].temperature_k < maximum {
                pending.push(Reverse((knots[next].temperature_k.to_bits(), source, next)));
            }
        }
        drop(pending);
        temperatures.push(maximum);
        let values = admit_knots(temperatures.len())?;
        let mut knots = Vec::with_capacity(temperatures.len());
        for temperature_k in temperatures {
            checkpoint()?;
            let mut samples = Vec::with_capacity(constituents.len());
            for constituent in &constituents {
                samples.push(PhotonMassConstituent {
                    mass_fraction: constituent.mass_fraction,
                    material: constituent
                        .material
                        .sample(temperature_k, &mut checkpoint)?,
                });
            }
            knots.push(PhotonMassTemperatureKnot {
                temperature_k,
                material: PhotonMassMaterial::from_mass_fractions(
                    data_id.clone(),
                    material_id,
                    samples,
                    kernel_applicability.clone(),
                    limits,
                    &mut checkpoint,
                )?,
            });
        }
        let mut result = Self::new(data_id, knots, limits, &mut checkpoint)?;
        let data =
            Arc::get_mut(&mut result.0).expect("newly admitted unshared temperature mixture");
        data.values = values;
        data.mixture = Some(TemperaturePhotonMassMixture {
            kernel_applicability,
            constituents,
        });
        Ok(result)
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

    /// Original temperature sources and constant composition, when prepared by
    /// [`Self::from_mass_fractions`]. Directly supplied tables have no mixture.
    pub fn mixture(&self) -> Option<&TemperaturePhotonMassMixture> {
        self.0.mixture.as_ref()
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
