//! Depth-conditioned escape on the original coherent finite-film absorption.
//!
//! Each source bin supplies an escape probability conditional on absorption at
//! that depth, rather than the aggregate yield of an entire illuminated film.
//! Its physical interval uses the native positive loss operator, preserving
//! standing-wave interference and longitudinal TM absorption without depth
//! quadrature. Bins are ordered by distance from the declared escape face;
//! exposed physical intervals are always reported in native entrance depth.
//!
//! Per-bin residual energy is attributed to its absorption bin for accounting.
//! Depositing that heat locally requires an explicit thermalization assumption
//! in the thermal consumer; no spatial electron-energy transport is solved.
//! Distinct wavelengths require independent threshold evaluations.

use super::{EscapeFace, LayerPhotoelectricError};
use crate::plane_wave::layers::{
    IlluminationBudget, PreparedLayerDepthInterval, PreparedLayerIllumination,
};
use fs_material::photoelectric::source::depth::{
    DepthPhotoelectricMaterial, MAX_PHOTOELECTRIC_DEPTH_BINS,
};
use fs_material::photoelectric::{PhotoelectricLaw, PhotoelectricResponse};
use fs_math::{c64::C64, dd::Dd};

struct PreparedBin<'p, 's, 'm> {
    absorption: PreparedLayerDepthInterval<'p, 's, 'm>,
    law: PhotoelectricLaw,
}

/// Fixed-frequency source profile and original positive interval operators.
/// Preparation retains O(B) bounded storage. Coherent queries cost O(B), and
/// factor-bank queries O(B*M), with no query allocation or optical re-solve.
pub struct PreparedLayerDepthPhotoelectric<'p, 's, 'm, 'd> {
    illumination: &'p PreparedLayerIllumination<'s, 'm>,
    material: &'d DepthPhotoelectricMaterial,
    bins: Vec<PreparedBin<'p, 's, 'm>>,
    layer_index: usize,
    face: EscapeFace,
    area_m2: f64,
    vacuum_wavelength_m: f64,
}

/// Additive physical rates and powers. No arithmetic average of bin electron
/// energies is reported: individual observations retain their calibrated means.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DepthPhotoelectricTotals {
    pub absorbed_power_w: f64,
    pub absorbed_photon_rate_s: f64,
    pub emitted_electron_rate_s: f64,
    pub surface_charging_current_a: f64,
    pub emitted_electron_current_a: f64,
    pub kinetic_power_w: f64,
    pub removal_power_w: f64,
    pub heat_power_w: f64,
}

/// One original source bin's observation and native entrance-depth support.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DepthPhotoelectricBinResponse {
    depth_interval_m: [f64; 2],
    absorbed_power_w_m2: f64,
    emission: PhotoelectricResponse,
}

impl DepthPhotoelectricBinResponse {
    /// Ordered depths from the physical layer entrance, including exit emission.
    pub const fn depth_interval_m(&self) -> [f64; 2] {
        self.depth_interval_m
    }
    pub const fn absorbed_power_w_m2(&self) -> f64 {
        self.absorbed_power_w_m2
    }
    pub const fn photoelectric(&self) -> &PhotoelectricResponse {
        &self.emission
    }
}

/// Immutable, complete observation with original source and optical ownership.
/// Bin order follows the source's distance from the declared escape face.
pub struct LayerDepthPhotoelectricResponse<'p, 's, 'm, 'd> {
    illumination: &'p PreparedLayerIllumination<'s, 'm>,
    material: &'d DepthPhotoelectricMaterial,
    layer_index: usize,
    face: EscapeFace,
    area_m2: f64,
    vacuum_wavelength_m: f64,
    bins: [Option<DepthPhotoelectricBinResponse>; MAX_PHOTOELECTRIC_DEPTH_BINS],
    bin_count: usize,
    totals: DepthPhotoelectricTotals,
}

impl<'p, 's, 'm, 'd> PreparedLayerDepthPhotoelectric<'p, 's, 'm, 'd> {
    /// Bind the exact complete source, calibrated thickness and exposed vacuum
    /// face. Normal Poynting density already contains incidence projection;
    /// `area_m2` is physical film area and adds no second cosine factor.
    pub fn try_new(
        illumination: &'p PreparedLayerIllumination<'s, 'm>,
        layer_index: usize,
        material: &'d DepthPhotoelectricMaterial,
        face: EscapeFace,
        area_m2: f64,
    ) -> Result<Self, LayerPhotoelectricError> {
        let invalid = LayerPhotoelectricError::Invalid;
        if !area_m2.is_finite() || area_m2 <= 0.0 {
            return Err(invalid("positive finite physical film area required"));
        }
        let stack = illumination.stack();
        let layer = stack
            .layer(layer_index)
            .ok_or_else(|| invalid("emitting film index is out of range"))?;
        let exposed = match face {
            EscapeFace::EntranceVacuum => layer_index == 0,
            EscapeFace::ExitVacuum => layer_index + 1 == stack.layer_count(),
        };
        if !exposed {
            return Err(invalid(
                "depth-profile emission requires an exposed vacuum face",
            ));
        }
        if material.resolved() != layer.response.material()
            || layer.thickness_m != material.thickness_m()
        {
            return Err(invalid(
                "optical film and depth profile require identical source and thickness",
            ));
        }
        let vacuum_wavelength_m =
            (core::f64::consts::TAU * 299_792_458.0) / stack.angular_frequency_rad_s();
        if !vacuum_wavelength_m.is_finite() || vacuum_wavelength_m <= 0.0 {
            return Err(invalid("vacuum photon wavelength is unrepresentable"));
        }
        let count = material.bins().len();
        if count == 0 || count > MAX_PHOTOELECTRIC_DEPTH_BINS {
            return Err(invalid(
                "depth-profile bin count is outside its bounded source domain",
            ));
        }
        let mut bins = Vec::new();
        bins.try_reserve_exact(count)
            .map_err(|_| invalid("depth-profile operator allocation refused"))?;
        for index in 0..count {
            let [start, end] = material.bin_bounds_fraction(index)?;
            let [start, end] = match face {
                EscapeFace::EntranceVacuum => [start, end],
                EscapeFace::ExitVacuum => [1.0 - end, 1.0 - start],
            }
            .map(|fraction| fraction * layer.thickness_m);
            if !start.is_finite() || !end.is_finite() || start >= end {
                return Err(invalid(
                    "a positive source depth bin has unrepresentable physical support",
                ));
            }
            bins.push(PreparedBin {
                absorption: illumination.prepare_depth_interval(layer_index, start, end)?,
                law: material.at_wavelength(index, vacuum_wavelength_m * 1e9)?,
            });
        }
        Ok(Self {
            illumination,
            material,
            bins,
            layer_index,
            face,
            area_m2,
            vacuum_wavelength_m,
        })
    }

    pub const fn illumination(&self) -> &'p PreparedLayerIllumination<'s, 'm> {
        self.illumination
    }
    pub const fn material(&self) -> &'d DepthPhotoelectricMaterial {
        self.material
    }
    pub const fn layer_index(&self) -> usize {
        self.layer_index
    }
    pub const fn escape_face(&self) -> EscapeFace {
        self.face
    }
    pub const fn area_m2(&self) -> f64 {
        self.area_m2
    }
    pub const fn vacuum_wavelength_m(&self) -> f64 {
        self.vacuum_wavelength_m
    }
    pub fn bin_count(&self) -> usize {
        self.bins.len()
    }
    /// Native depth support for preparing intersections with thermal cells.
    pub fn bin_depth_interval_m(&self, index: usize) -> Option<[f64; 2]> {
        self.bins
            .get(index)
            .map(|bin| [bin.absorption.start_depth_m(), bin.absorption.end_depth_m()])
    }
    /// Original conditional escape law at this preparation's vacuum wavelength.
    pub fn bin_law(&self, index: usize) -> Option<&PhotoelectricLaw> {
        self.bins.get(index).map(|bin| &bin.law)
    }

    /// One coherent pair of incoming peak tangential E amplitudes [V/m].
    pub fn evaluate(
        &self,
        input: [C64; 2],
        cancelled: impl FnMut() -> bool,
    ) -> Result<LayerDepthPhotoelectricResponse<'p, 's, 'm, 'd>, LayerPhotoelectricError> {
        self.evaluate_coherence_modes(
            &[input],
            IlluminationBudget {
                maximum_modes: 1,
                maximum_layer_mode_pairs: self.bins.len(),
            },
            cancelled,
        )
    }

    /// Independent face sources add powers after each native positive square.
    pub fn evaluate_independent_sources(
        &self,
        input: [C64; 2],
        cancelled: impl FnMut() -> bool,
    ) -> Result<LayerDepthPhotoelectricResponse<'p, 's, 'm, 'd>, LayerPhotoelectricError> {
        self.evaluate_coherence_modes(
            &[[input[0], C64::ZERO], [C64::ZERO, input[1]]],
            IlluminationBudget {
                maximum_modes: 2,
                maximum_layer_mode_pairs: 2 * self.bins.len(),
            },
            cancelled,
        )
    }

    /// Same-frequency covariance factors. The complete B*M work product is
    /// admitted first. All bins and totals remain private until every native
    /// interval, physical response and final cancellation checkpoint succeeds.
    pub fn evaluate_coherence_modes(
        &self,
        factors: &[[C64; 2]],
        budget: IlluminationBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<LayerDepthPhotoelectricResponse<'p, 's, 'm, 'd>, LayerPhotoelectricError> {
        let invalid = LayerPhotoelectricError::Invalid;
        if factors.len() > 4096
            || factors.len() > budget.maximum_modes
            || factors
                .len()
                .checked_mul(self.bins.len())
                .is_none_or(|work| work > budget.maximum_layer_mode_pairs)
        {
            return Err(invalid(
                "depth-profile source factors exceed the complete bin-mode work budget",
            ));
        }
        if cancelled() {
            return Err(invalid("depth-profile observation cancelled"));
        }
        let mut output = [None; MAX_PHOTOELECTRIC_DEPTH_BINS];
        let mut sums = [Dd::ZERO; 7];
        for (index, bin) in self.bins.iter().enumerate() {
            let density =
                bin.absorption
                    .evaluate_coherence_modes(factors, budget, &mut cancelled)?;
            let power = density * self.area_m2;
            if !power.is_finite() || (density > 0.0 && power == 0.0) {
                return Err(invalid(
                    "physical absorbed depth-bin power is unrepresentable",
                ));
            }
            let emission = bin.law.evaluate(self.vacuum_wavelength_m, power)?;
            for (sum, value) in sums.iter_mut().zip([
                emission.absorbed_power_w,
                emission.absorbed_photon_rate_s,
                emission.emitted_electron_rate_s,
                emission.surface_charging_current_a,
                emission.kinetic_power_w,
                emission.removal_power_w,
                emission.heat_power_w,
            ]) {
                *sum = *sum + Dd::from_f64(value);
            }
            output[index] = Some(DepthPhotoelectricBinResponse {
                depth_interval_m: [bin.absorption.start_depth_m(), bin.absorption.end_depth_m()],
                absorbed_power_w_m2: density,
                emission,
            });
        }
        let values = sums.map(Dd::to_f64);
        if values.iter().any(|v| !v.is_finite() || *v < 0.0) {
            return Err(invalid(
                "depth-profile aggregate power or charge is unrepresentable",
            ));
        }
        if cancelled() {
            return Err(invalid("depth-profile observation cancelled"));
        }
        Ok(LayerDepthPhotoelectricResponse {
            illumination: self.illumination,
            material: self.material,
            layer_index: self.layer_index,
            face: self.face,
            area_m2: self.area_m2,
            vacuum_wavelength_m: self.vacuum_wavelength_m,
            bins: output,
            bin_count: self.bins.len(),
            totals: DepthPhotoelectricTotals {
                absorbed_power_w: values[0],
                absorbed_photon_rate_s: values[1],
                emitted_electron_rate_s: values[2],
                surface_charging_current_a: values[3],
                emitted_electron_current_a: -values[3],
                kinetic_power_w: values[4],
                removal_power_w: values[5],
                heat_power_w: values[6],
            },
        })
    }
}

impl<'p, 's, 'm, 'd> LayerDepthPhotoelectricResponse<'p, 's, 'm, 'd> {
    pub const fn illumination(&self) -> &'p PreparedLayerIllumination<'s, 'm> {
        self.illumination
    }
    pub const fn material(&self) -> &'d DepthPhotoelectricMaterial {
        self.material
    }
    pub const fn layer_index(&self) -> usize {
        self.layer_index
    }
    pub const fn escape_face(&self) -> EscapeFace {
        self.face
    }
    pub const fn area_m2(&self) -> f64 {
        self.area_m2
    }
    pub const fn vacuum_wavelength_m(&self) -> f64 {
        self.vacuum_wavelength_m
    }
    pub const fn bin_count(&self) -> usize {
        self.bin_count
    }
    pub fn bin(&self, index: usize) -> Option<&DepthPhotoelectricBinResponse> {
        self.bins.get(index).and_then(Option::as_ref)
    }
    pub fn bins(&self) -> impl Iterator<Item = &DepthPhotoelectricBinResponse> {
        self.bins[..self.bin_count].iter().flatten()
    }
    pub const fn totals(&self) -> &DepthPhotoelectricTotals {
        &self.totals
    }
}
