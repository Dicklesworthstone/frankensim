//! Polarized thermal emission with the original temperature of every crystal.
//!
//! The complete coherent Maxwell stack determines reciprocal absorption in
//! each original layer. Independent local thermal sources then contribute
//! `B_lambda(T_layer) * A_layer / 2` to a complex outgoing analyzer. A layer's
//! positive dielectric work includes reflections and polarization conversion
//! throughout the whole stack; it is not the absorption of an isolated film.
//! This is the planar-layer form of the local Kirchhoff law in Greffet et al.,
//! Phys. Rev. X 8, 021008 (2018), <https://doi.org/10.1103/PhysRevX.8.021008>.
//!
//! Each layer is internally isothermal and in local thermal equilibrium.
//! Original optical coefficients, global axes and geometry remain fixed.
//! Temperatures are never averaged, and thermal source amplitudes never
//! interfere. This far-field entrance-vacuum observation does not refresh
//! material laws, debit thermal energy or introduce near-field exchange.

use super::{
    LIGHT_M_S, ThermalSurfaceError,
    crystal_emission::{
        CrystalEmissionError, admit_temperature, admit_vacuum, material_temperature,
        reciprocal_analyzer_electric, unit_plane,
    },
    planck_radiance_w_m2_sr_nm,
};
use crate::polarization::Stokes;
use fs_em::plane_wave::{
    angular_spectrum::AngularSpectrumError,
    uniaxial::oriented::{
        absorbing::stack::PreparedAbsorbingStackPolarization, stack::OrientedUniaxialStack,
    },
};
use fs_math::{c64::C64, dd::Dd, det};

type Result<T> = core::result::Result<T, CrystalStackEmissionError>;

/// One uniform temperature per original finite crystal, in authored order.
/// The source retains each material, principal spectrum, permeability, optic
/// axis, thickness and exterior law. Preparation never copies an index fit.
pub struct UniaxialStackThermalEmission<'f, 's, 'm> {
    source: &'f OrientedUniaxialStack<'s, 'm>,
    temperatures_k: Vec<f64>,
}

impl<'f, 's, 'm> UniaxialStackThermalEmission<'f, 's, 'm> {
    /// Read an existing typed absolute-temperature coordinate from EACH
    /// original layer's complete material snapshot. Axis names are in the
    /// same order as `source.layers()`; distinct materials may name them
    /// differently. Every layer needs exactly one coordinate.
    pub fn try_from_material_temperatures(
        source: &'f OrientedUniaxialStack<'s, 'm>,
        temperature_axes: &[&str],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(&mut cancelled)?;
        admit_temperature_count(source.layers().len(), temperature_axes.len())?;
        let mut temperatures_k = reserved(source.layers().len())?;
        for (index, (layer, &axis)) in source.layers().iter().zip(temperature_axes).enumerate() {
            poll(&mut cancelled)?;
            let temperature = material_temperature(layer.crystal_response().material(), axis)
                .map_err(|source| CrystalStackEmissionError::LayerSource {
                    layer: index,
                    source,
                })?;
            temperatures_k.push(temperature);
        }
        poll(&mut cancelled)?;
        Ok(Self {
            source,
            temperatures_k,
        })
    }

    /// Declare one uniform temperature per layer while holding every original
    /// optical snapshot fixed. This does NOT resolve material claims at the
    /// supplied temperatures or imply thermo-optical feedback. The temperature
    /// slice is copied once so future observations retain this declaration.
    pub fn try_layer_temperatures(
        source: &'f OrientedUniaxialStack<'s, 'm>,
        temperatures_k: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(&mut cancelled)?;
        admit_temperature_count(source.layers().len(), temperatures_k.len())?;
        let mut temperatures = reserved(temperatures_k.len())?;
        for (index, &temperature) in temperatures_k.iter().enumerate() {
            poll(&mut cancelled)?;
            admit_temperature(temperature).map_err(|source| {
                CrystalStackEmissionError::LayerSource {
                    layer: index,
                    source,
                }
            })?;
            temperatures.push(temperature);
        }
        poll(&mut cancelled)?;
        Ok(Self {
            source,
            temperatures_k: temperatures,
        })
    }

    pub const fn source(&self) -> &'f OrientedUniaxialStack<'s, 'm> {
        self.source
    }

    /// Uniform temperatures [K], in original authored layer order.
    pub fn layer_temperatures_k(&self) -> &[f64] {
        &self.temperatures_k
    }

    /// Prepare an outgoing mode in the entrance vacuum half-space.
    ///
    /// The ray is `(sin(theta)*u.x, sin(theta)*u.y, -cos(theta))`, with
    /// `0 < outgoing_cosine <= 1`. The nonzero transverse u is normalized;
    /// it also declares the analyzer frame at normal incidence. The native
    /// reciprocal query keeps cosine and reverses u, retaining every crystal
    /// axis in the original global frame.
    ///
    /// Two native polarization responses and one Planck radiance per original
    /// layer are prepared once. Subsequent analyzer queries allocate nothing,
    /// reevaluate no materials and perform no Maxwell boundary solve.
    pub fn prepare(
        &self,
        wavelength_nm: f64,
        outgoing_cosine: f64,
        outgoing_plane_axis: [f64; 2],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedUniaxialStackEmission<'f, 's, 'm>> {
        poll(&mut cancelled)?;
        if !wavelength_nm.is_finite() || wavelength_nm <= 0.0 {
            return Err(ThermalSurfaceError::Wavelength.into());
        }
        if !outgoing_cosine.is_finite() || outgoing_cosine <= 0.0 || outgoing_cosine > 1.0 {
            return Err(CrystalEmissionError::EmissionDirection.into());
        }
        let plane = unit_plane(outgoing_plane_axis)?;
        let omega = (core::f64::consts::TAU * LIGHT_M_S * 1.0e9) / wavelength_nm;
        let polarization = self.source.prepare_absorbing_polarization(
            omega,
            outgoing_cosine,
            [-plane[0], -plane[1]],
            &mut cancelled,
        )?;
        admit_vacuum(
            polarization.incident_phase_index(),
            polarization.incident_permeability_h_m(),
        )?;
        let mut layers = reserved(self.temperatures_k.len())?;
        for (index, &temperature_k) in self.temperatures_k.iter().enumerate() {
            poll(&mut cancelled)?;
            let blackbody_radiance = planck_radiance_w_m2_sr_nm(temperature_k, wavelength_nm)
                .map_err(|source| CrystalStackEmissionError::LayerSource {
                    layer: index,
                    source: source.into(),
                })?;
            layers.push(ThermalLayer {
                temperature_k,
                blackbody_radiance,
            });
        }
        poll(&mut cancelled)?;
        Ok(PreparedUniaxialStackEmission {
            polarization,
            layers,
            wavelength_nm,
            cosine: outgoing_cosine,
            plane,
        })
    }
}

struct ThermalLayer {
    temperature_k: f64,
    blackbody_radiance: f64,
}

/// Polarized spectral radiance from the whole coherent stack, with independent
/// thermal contributions from each original layer. The original temperatures
/// and source remain retained even after the preparation owner is dropped.
pub struct PreparedUniaxialStackEmission<'f, 's, 'm> {
    polarization: PreparedAbsorbingStackPolarization<'f, 's, 'm>,
    layers: Vec<ThermalLayer>,
    wavelength_nm: f64,
    cosine: f64,
    plane: [f64; 2],
}

impl<'f, 's, 'm> PreparedUniaxialStackEmission<'f, 's, 'm> {
    pub const fn source(&self) -> &'f OrientedUniaxialStack<'s, 'm> {
        self.polarization.source()
    }

    pub const fn wavelength_nm(&self) -> f64 {
        self.wavelength_nm
    }

    pub fn layer_count(&self) -> usize {
        self.layers.len()
    }

    pub fn layer_temperature_k(&self, index: usize) -> Option<f64> {
        self.layers.get(index).map(|layer| layer.temperature_k)
    }

    pub fn outgoing_direction(&self) -> [f64; 3] {
        let sine = det::sqrt((1.0 - self.cosine) * (1.0 + self.cosine));
        [sine * self.plane[0], sine * self.plane[1], -self.cosine]
    }

    /// Real global `[s,p]` axes, with `s cross p` along the outgoing ray.
    /// The same frame is used for every layer and every returned Stokes value.
    pub fn polarization_axes(&self) -> [[f64; 3]; 2] {
        let sine = det::sqrt((1.0 - self.cosine) * (1.0 + self.cosine));
        [
            [-self.plane[1], self.plane[0], 0.0],
            [
                self.cosine * self.plane[0],
                self.cosine * self.plane[1],
                sine,
            ],
        ]
    }

    /// One original layer's radiance through a complex `[a_s,a_p]` analyzer,
    /// W/(m² sr nm). The analyzer is normalized; its conjugated global field
    /// probes the reciprocal stack. Work is fixed-size for this one layer.
    /// A transparent layer contributes exactly zero at any temperature.
    ///
    /// This direct per-layer result retains a weak thermal contribution that
    /// a rounded total or Stokes summary may be unable to represent separately.
    pub fn layer_analyzer_radiance_w_m2_sr_nm(
        &self,
        layer_index: usize,
        analyzer: [C64; 2],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<f64> {
        poll(&mut cancelled)?;
        let (electric, incident) = self.incident_for_analyzer(analyzer)?;
        let (_, radiance) = self.layer_contribution(layer_index, electric, incident)?;
        poll(&mut cancelled)?;
        Ok(radiance)
    }

    /// Whole-stack radiance through one normalized complex analyzer.
    /// Independent layer radiances are accumulated with double-double sums.
    /// Each layer uses its own Planck temperature and its positive reciprocal
    /// volume work in the COMPLETE stack, including all internal reflections.
    pub fn analyzer_radiance_w_m2_sr_nm(
        &self,
        analyzer: [C64; 2],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<f64> {
        poll(&mut cancelled)?;
        let (electric, incident) = self.incident_for_analyzer(analyzer)?;
        let mut radiance = Dd::ZERO;
        let mut absorbed = Dd::ZERO;
        for index in 0..self.layers.len() {
            poll(&mut cancelled)?;
            let (fraction, contribution) = self.layer_contribution(index, electric, incident)?;
            absorbed = absorbed + Dd::from_f64(fraction);
            radiance = radiance + Dd::from_f64(contribution);
        }
        poll(&mut cancelled)?;
        let radiance = radiance.to_f64();
        let absorbed = absorbed.to_f64();
        if !radiance.is_finite()
            || radiance < 0.0
            || !absorbed.is_finite()
            || absorbed > 1.0 + 1.0e-9
        {
            return Err(ThermalSurfaceError::Unrepresentable.into());
        }
        Ok(radiance)
    }

    /// Total unanalysed radiance; independent orthogonal polarizations add.
    pub fn radiance_w_m2_sr_nm(&self, mut cancelled: impl FnMut() -> bool) -> Result<f64> {
        let s = self.analyzer_radiance_w_m2_sr_nm([C64::ONE, C64::ZERO], &mut cancelled)?;
        let p = self.analyzer_radiance_w_m2_sr_nm([C64::ZERO, C64::ONE], &mut cancelled)?;
        let radiance = (Dd::from_f64(s) + Dd::from_f64(p)).to_f64();
        if !radiance.is_finite() {
            return Err(ThermalSurfaceError::Unrepresentable.into());
        }
        Ok(radiance)
    }

    /// Conventional Stokes radiance in the common outgoing `[s,p]` frame.
    /// +Q is s and `V = 2 Im(<E_s conj(E_p)>)`. Six analyzer probes reuse
    /// the same native cache. Differences may lose a weak component; direct
    /// analyzer and per-layer queries preserve their own positive work.
    pub fn stokes(&self, mut cancelled: impl FnMut() -> bool) -> Result<Stokes> {
        let s = self.analyzer_radiance_w_m2_sr_nm([C64::ONE, C64::ZERO], &mut cancelled)?;
        let p = self.analyzer_radiance_w_m2_sr_nm([C64::ZERO, C64::ONE], &mut cancelled)?;
        let plus = self.analyzer_radiance_w_m2_sr_nm([C64::ONE, C64::ONE], &mut cancelled)?;
        let minus = self.analyzer_radiance_w_m2_sr_nm([C64::ONE, -C64::ONE], &mut cancelled)?;
        let plus_i =
            self.analyzer_radiance_w_m2_sr_nm([C64::ONE, C64::new(0.0, 1.0)], &mut cancelled)?;
        let minus_i =
            self.analyzer_radiance_w_m2_sr_nm([C64::ONE, C64::new(0.0, -1.0)], &mut cancelled)?;
        Stokes::try_new(s + p, s - p, plus - minus, minus_i - plus_i)
            .ok_or_else(|| ThermalSurfaceError::Unrepresentable.into())
    }

    fn incident_for_analyzer(&self, analyzer: [C64; 2]) -> Result<([C64; 2], f64)> {
        let electric = reciprocal_analyzer_electric(self.cosine, self.plane, analyzer)?;
        let incident = self.polarization.incident_normal_power_w_m2(electric)?;
        if !incident.is_finite() || incident <= 0.0 {
            return Err(ThermalSurfaceError::Unrepresentable.into());
        }
        Ok((electric, incident))
    }

    fn layer_contribution(
        &self,
        index: usize,
        electric: [C64; 2],
        incident: f64,
    ) -> Result<(f64, f64)> {
        let layer = self
            .layers
            .get(index)
            .ok_or(CrystalStackEmissionError::LayerIndex {
                index,
                count: self.layers.len(),
            })?;
        let depth = self.source().layers()[index].thickness_m();
        let power = self
            .polarization
            .layer_absorbed_power_w(index, electric, 1.0, 0.0, depth)?;
        let fraction = power / incident;
        if !fraction.is_finite() || fraction < 0.0 || fraction > 1.0 + 1.0e-9 {
            return Err(ThermalSurfaceError::Unrepresentable.into());
        }
        // B contains both thermal polarizations. Multiply the weak fraction
        // before halving, preserving representable small contributions.
        let radiance = 0.5 * (fraction * layer.blackbody_radiance);
        if !radiance.is_finite() || radiance < 0.0 {
            return Err(ThermalSurfaceError::Unrepresentable.into());
        }
        Ok((fraction, radiance))
    }
}

/// Admission or cancellation returns no partial emitted observation.
#[derive(Debug)]
pub enum CrystalStackEmissionError {
    /// Every original layer requires exactly one temperature or axis name.
    LayerTemperatures { expected: usize, observed: usize },
    /// Original temperature coordinate or Planck source for this layer refused.
    LayerSource {
        layer: usize,
        source: CrystalEmissionError,
    },
    /// The requested authored layer does not exist.
    LayerIndex { index: usize, count: usize },
    /// Shared direction, analyzer, vacuum-port or radiance admission refused.
    Emission(CrystalEmissionError),
    /// Original Maxwell source, wavelength support or field integration refused.
    Electromagnetic(AngularSpectrumError),
    /// Bounded preparation storage could not be reserved.
    Allocation,
    /// The caller cancelled before an observation was published.
    Cancelled,
}

impl From<CrystalEmissionError> for CrystalStackEmissionError {
    fn from(value: CrystalEmissionError) -> Self {
        Self::Emission(value)
    }
}
impl From<ThermalSurfaceError> for CrystalStackEmissionError {
    fn from(value: ThermalSurfaceError) -> Self {
        Self::Emission(value.into())
    }
}
impl From<AngularSpectrumError> for CrystalStackEmissionError {
    fn from(value: AngularSpectrumError) -> Self {
        match value {
            AngularSpectrumError::Cancelled => Self::Cancelled,
            other => Self::Electromagnetic(other),
        }
    }
}
impl core::fmt::Display for CrystalStackEmissionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "crystal-stack thermal emission refused: {self:?}")
    }
}
impl core::error::Error for CrystalStackEmissionError {}

fn admit_temperature_count(expected: usize, observed: usize) -> Result<()> {
    if expected != observed {
        Err(CrystalStackEmissionError::LayerTemperatures { expected, observed })
    } else {
        Ok(())
    }
}
fn reserved<T>(count: usize) -> Result<Vec<T>> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(count)
        .map_err(|_| CrystalStackEmissionError::Allocation)?;
    Ok(values)
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(CrystalStackEmissionError::Cancelled)
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests;
