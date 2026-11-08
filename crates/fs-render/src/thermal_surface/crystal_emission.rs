//! Polarized far-field emission from one original absorbing uniaxial film.
//!
//! For a reciprocal film in local thermal equilibrium, emissivity into a
//! complex outgoing analyzer equals absorptivity of its phase-conjugated,
//! reversed incoming wave. See Miller, Zhu and Fan, PNAS 114, 4336 (2017),
//! <https://doi.org/10.1073/pnas.1701606114>. Native Maxwell preparation retains
//! the original two principal spectra, scalar permeability and optic axis.
//! The positive volume integral of dielectric work supplies absorption, never
//! the rounded difference of incident, reflected and transmitted powers.
//!
//! This owner observes the entrance vacuum half-space of a homogeneous,
//! isothermal film. The opposite lossless medium retains its original law but
//! contributes no invented thermal source. This is a radiance observation;
//! it neither advances/debits heat nor refreshes temperature-dependent optics.

use super::{LIGHT_M_S, ThermalSurfaceError, planck_radiance_w_m2_sr_nm};
use crate::polarization::Stokes;
use fs_em::plane_wave::{
    ReflectionError,
    uniaxial::oriented::{
        absorbing::slab::PreparedAbsorbingSlabPolarization, slab::OrientedUniaxialSlab,
    },
};
use fs_math::{c64::C64, det};
use fs_qty::{
    QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

// Admit customary rounded SI vacuum constants, including pre-2019 mu0 and
// modern measured values. This relative 1e-9 tolerance is an explicit port
// convention, not an exterior-air or arbitrary-material approximation.
const VACUUM_MU_H_M: f64 = 4.0e-7 * core::f64::consts::PI;
const VACUUM_RELATIVE_TOLERANCE: f64 = 1.0e-9;
type Result<T> = core::result::Result<T, CrystalEmissionError>;

/// Original finite film at one declared uniform thermal temperature.
///
/// All optical coefficients and geometry remain with the borrowed native
/// source. Use [`Self::try_from_material_temperature`] to read temperature
/// from that exact source snapshot. Changed optical material state requires
/// preparing a new native slab, without a fitted emissivity or index copy.
pub struct UniaxialSlabThermalEmission<'f, 's, 'm> {
    source: &'f OrientedUniaxialSlab<'s, 'm>,
    temperature_k: f64,
}

impl<'f, 's, 'm> UniaxialSlabThermalEmission<'f, 's, 'm> {
    /// Bind the crystal's existing typed static absolute-temperature axis.
    /// Both principal optical sources already share this complete material
    /// snapshot. A missing, untyped or differently typed coordinate refuses.
    pub fn try_from_material_temperature(
        source: &'f OrientedUniaxialSlab<'s, 'm>,
        temperature_axis: &str,
    ) -> Result<Self> {
        let point = source
            .entrance_interface()
            .crystal_response()
            .material()
            .query();
        let absolute = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature,
            ValueForm::Static,
        ));
        if point.axis_quantities().get(temperature_axis) != Some(&absolute) {
            return Err(CrystalEmissionError::TemperatureAxis);
        }
        let temperature_k = point
            .axes()
            .get(temperature_axis)
            .copied()
            .ok_or(CrystalEmissionError::TemperatureAxis)?;
        Self::try_isothermal(source, temperature_k)
    }

    /// Explicitly declare one uniform temperature while holding the original
    /// optical snapshot fixed. This constructor does NOT substitute this
    /// temperature into material claims or imply thermo-optical feedback.
    pub fn try_isothermal(
        source: &'f OrientedUniaxialSlab<'s, 'm>,
        temperature_k: f64,
    ) -> Result<Self> {
        if !temperature_k.is_finite() || temperature_k <= 0.0 {
            return Err(ThermalSurfaceError::Temperature.into());
        }
        Ok(Self {
            source,
            temperature_k,
        })
    }

    pub const fn source(&self) -> &'f OrientedUniaxialSlab<'s, 'm> {
        self.source
    }

    pub const fn temperature_k(&self) -> f64 {
        self.temperature_k
    }

    /// Prepare one outgoing entrance-side direction and vacuum wavelength.
    ///
    /// The original crystal occupies `0 <= z <= d`. The outgoing ray is
    /// `(sin(theta)*u.x, sin(theta)*u.y, -cos(theta))`, with strictly positive
    /// `outgoing_cosine`. `outgoing_plane_axis` supplies nonzero transverse u;
    /// it is normalized and also fixes the polarization frame at normal
    /// incidence. No absolute cosine, side swap or azimuth inference occurs.
    ///
    /// Exactly two native polarization responses are solved. Subsequent
    /// analyzers combine these original coherent fields before integrating
    /// loss; they need no material queries, boundary solves or allocations.
    /// Both original exterior laws and principal spectra admit wavelength.
    pub fn prepare(
        &self,
        wavelength_nm: f64,
        outgoing_cosine: f64,
        outgoing_plane_axis: [f64; 2],
    ) -> Result<PreparedUniaxialSlabEmission<'f, 's, 'm>> {
        if !wavelength_nm.is_finite() || wavelength_nm <= 0.0 {
            return Err(ThermalSurfaceError::Wavelength.into());
        }
        if !outgoing_cosine.is_finite() || outgoing_cosine <= 0.0 || outgoing_cosine > 1.0 {
            return Err(CrystalEmissionError::EmissionDirection);
        }
        let plane = unit_plane(outgoing_plane_axis)?;
        let omega = (core::f64::consts::TAU * LIGHT_M_S * 1.0e9) / wavelength_nm;
        let polarization = self.source.prepare_absorbing_polarization(
            omega,
            outgoing_cosine,
            [-plane[0], -plane[1]],
        )?;
        let index = polarization.incident_phase_index();
        let mu = polarization.incident_permeability_h_m();
        if !index.is_finite()
            || !mu.is_finite()
            || (index - 1.0).abs() > VACUUM_RELATIVE_TOLERANCE
            || (mu / VACUUM_MU_H_M - 1.0).abs() > VACUUM_RELATIVE_TOLERANCE
        {
            return Err(CrystalEmissionError::ObservationMedium);
        }
        Ok(PreparedUniaxialSlabEmission {
            polarization,
            temperature_k: self.temperature_k,
            wavelength_nm,
            cosine: outgoing_cosine,
            plane,
            blackbody_radiance: planck_radiance_w_m2_sr_nm(self.temperature_k, wavelength_nm)?,
        })
    }
}

/// Source-bound spectral radiance in one fixed outgoing polarization frame.
///
/// For transverse unit u and positive cosine c, `s=(-u.y,u.x,0)` and
/// `p=(c*u.x,c*u.y,sin(theta))`; `s cross p` is the outgoing direction.
/// Analyzer components are complex electric amplitudes along these real
/// axes. No artificial coherence is assigned to independent thermal sources.
pub struct PreparedUniaxialSlabEmission<'f, 's, 'm> {
    polarization: PreparedAbsorbingSlabPolarization<'f, 's, 'm>,
    temperature_k: f64,
    wavelength_nm: f64,
    cosine: f64,
    plane: [f64; 2],
    blackbody_radiance: f64,
}

impl<'f, 's, 'm> PreparedUniaxialSlabEmission<'f, 's, 'm> {
    pub const fn source(&self) -> &'f OrientedUniaxialSlab<'s, 'm> {
        self.polarization.source()
    }

    pub const fn temperature_k(&self) -> f64 {
        self.temperature_k
    }

    pub const fn wavelength_nm(&self) -> f64 {
        self.wavelength_nm
    }

    /// Global outgoing unit wavevector into the original entrance half-space.
    pub fn outgoing_direction(&self) -> [f64; 3] {
        let sine = det::sqrt((1.0 - self.cosine) * (1.0 + self.cosine));
        [sine * self.plane[0], sine * self.plane[1], -self.cosine]
    }

    /// Real global transverse `[s,p]` unit axes for every returned observable.
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

    /// Radiance through one ideal complex analyzer, W/(m² sr nm).
    ///
    /// Nonzero finite `[a_s,a_p]` is normalized, so overall amplitude and
    /// phase do not change the result. Its conjugated global electric field
    /// illuminates the reversed ray; original positive dielectric work gives
    /// the matched absorptance. Planck B already contains two independent
    /// polarizations, hence the single-analyzer factor of one half.
    ///
    /// This direct method retains weak channels that cannot be recovered
    /// from rounded I and Q. It is authoritative for subsequent analyzers;
    /// do not reconstruct a weak radiance by subtracting Stokes components.
    pub fn analyzer_radiance_w_m2_sr_nm(&self, analyzer: [C64; 2]) -> Result<f64> {
        let [s, p] = unit_analyzer(analyzer)?;
        let [ux, uy] = self.plane;
        let electric = [
            (s.scale(-uy) + p.scale(self.cosine * ux)).conj(),
            (s.scale(ux) + p.scale(self.cosine * uy)).conj(),
        ];
        let field = self.polarization.prepare_field(electric)?;
        let power = field
            .slab_absorption(1.0, 0.0, self.source().thickness_m())?
            .absorbed_power_w();
        let incident = field.response().incident_normal_power_w_m2;
        let absorptance = power / incident;
        if !incident.is_finite()
            || incident <= 0.0
            || !absorptance.is_finite()
            || absorptance < 0.0
            || absorptance > 1.0 + 1.0e-9
        {
            return Err(ThermalSurfaceError::Unrepresentable.into());
        }
        // Multiply weak absorptance before halving; do not underflow a
        // representable weak-channel radiance by prematurely halving it.
        let radiance = 0.5 * (absorptance * self.blackbody_radiance);
        if !radiance.is_finite() || radiance < 0.0 {
            return Err(ThermalSurfaceError::Unrepresentable.into());
        }
        Ok(radiance)
    }

    /// Total unanalysed radiance; two independent orthogonal modes add.
    pub fn radiance_w_m2_sr_nm(&self) -> Result<f64> {
        let s = self.analyzer_radiance_w_m2_sr_nm([C64::ONE, C64::ZERO])?;
        let p = self.analyzer_radiance_w_m2_sr_nm([C64::ZERO, C64::ONE])?;
        let radiance = s + p;
        if !radiance.is_finite() {
            return Err(ThermalSurfaceError::Unrepresentable.into());
        }
        Ok(radiance)
    }

    /// Conventional incoherent Stokes radiance in W/(m² sr nm).
    ///
    /// +Q is s, +U is the equal in-phase s/p analyzer, and
    /// `V = 2 Im(<E_s conj(E_p)>)`. Six physical analyzer queries reuse the
    /// same prepared modes. I/Q and U/V differences can lose a much weaker
    /// contribution through rounding; direct analyzer queries retain their
    /// own positive work. A numerically unphysical summary refuses instead
    /// of clipping polarization or changing any analyzer's native result.
    pub fn stokes(&self) -> Result<Stokes> {
        let s = self.analyzer_radiance_w_m2_sr_nm([C64::ONE, C64::ZERO])?;
        let p = self.analyzer_radiance_w_m2_sr_nm([C64::ZERO, C64::ONE])?;
        let plus = self.analyzer_radiance_w_m2_sr_nm([C64::ONE, C64::ONE])?;
        let minus = self.analyzer_radiance_w_m2_sr_nm([C64::ONE, -C64::ONE])?;
        let plus_i = self.analyzer_radiance_w_m2_sr_nm([C64::ONE, C64::new(0.0, 1.0)])?;
        let minus_i = self.analyzer_radiance_w_m2_sr_nm([C64::ONE, C64::new(0.0, -1.0)])?;
        Stokes::try_new(s + p, s - p, plus - minus, minus_i - plus_i)
            .ok_or(ThermalSurfaceError::Unrepresentable.into())
    }
}

/// Unsupported source, observation or arithmetic leaves no emitted snapshot.
#[derive(Debug)]
pub enum CrystalEmissionError {
    /// The original material needs an existing typed absolute-temperature axis.
    TemperatureAxis,
    /// The outgoing entrance-side direction must have positive normal cosine.
    EmissionDirection,
    /// The native incident port must match vacuum index and permeability.
    ObservationMedium,
    /// A finite nonzero pair of complex transverse analyzer amplitudes is required.
    Analyzer,
    /// Original Maxwell source, spectral support, field or work refused.
    Electromagnetic(ReflectionError),
    /// Original Planck source or final radiance is unsupported/unrepresentable.
    Thermal(ThermalSurfaceError),
}

impl From<ReflectionError> for CrystalEmissionError {
    fn from(value: ReflectionError) -> Self {
        Self::Electromagnetic(value)
    }
}

impl From<ThermalSurfaceError> for CrystalEmissionError {
    fn from(value: ThermalSurfaceError) -> Self {
        Self::Thermal(value)
    }
}

impl core::fmt::Display for CrystalEmissionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "crystal thermal emission refused: {self:?}")
    }
}

impl core::error::Error for CrystalEmissionError {}

fn unit_plane(axis: [f64; 2]) -> Result<[f64; 2]> {
    let scale = axis[0].abs().max(axis[1].abs());
    if axis.iter().any(|x| !x.is_finite()) || scale == 0.0 {
        return Err(CrystalEmissionError::EmissionDirection);
    }
    let scaled = axis.map(|x| x / scale);
    let norm = det::hypot(scaled[0], scaled[1]);
    Ok(scaled.map(|x| x / norm))
}

fn unit_analyzer(analyzer: [C64; 2]) -> Result<[C64; 2]> {
    let values = [
        analyzer[0].re,
        analyzer[0].im,
        analyzer[1].re,
        analyzer[1].im,
    ];
    let scale = values.iter().fold(0.0_f64, |scale, x| scale.max(x.abs()));
    if values.iter().any(|x| !x.is_finite()) || scale == 0.0 {
        return Err(CrystalEmissionError::Analyzer);
    }
    let a = values.map(|x| x / scale);
    let norm = det::hypot(det::hypot(a[0], a[1]), det::hypot(a[2], a[3]));
    Ok([
        C64::new(a[0] / norm, a[1] / norm),
        C64::new(a[2] / norm, a[3] / norm),
    ])
}

#[cfg(test)]
mod tests;
