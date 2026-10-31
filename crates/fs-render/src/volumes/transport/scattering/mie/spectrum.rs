//! Visible-band preparation into the native absorption/scattering owners.

use super::super::{MassAbsorptionSpectrum, MassScatteringSpectrum, PhaseFunction};
use super::{
    MieAngularSolution, MieBudget, MieError, MiePhaseTable, MiePopulation, MiePopulationBudget,
    MieSphere, poll,
};
use fs_exec::Cx;

/// Ready-to-bind spectral optical laws with explicit angular and wavelength
/// interpolation diagnostics. Grid density must be dispersed particle mass
/// in kg/m³; the carrier's own optical response is a separate medium law.
///
/// The native ambient scattering grid has refractive index one. A nonunit
/// host therefore requires the existing dielectric-interior binding with the
/// matching host index; these coefficients do not update an interface index,
/// bend rays or create a refractive medium by themselves.
#[derive(Clone, Debug, PartialEq)]
pub struct MieOpticalSpectrum {
    /// Real host index used for the scattering solve; ray/interface ownership
    /// must supply this same host independently when it is not one.
    pub host_index_n: f64,
    /// Particle-mass absorption, m²/kg, linear between declared wavelengths.
    pub absorption: MassAbsorptionSpectrum,
    /// Particle-mass scattering, m²/kg, on the same wavelength grid.
    pub scattering: MassScatteringSpectrum,
    /// Wavelength-dependent normalized angular law and matched sampler.
    pub phase: PhaseFunction,
    /// One report per scattering wavelength in the same order as the supplied
    /// knots. Zero-scattering knots have no angular reports.
    pub angular_reports: Vec<MiePhaseTable>,
    /// Maximum relative absorption/scattering discrepancy at spectral interval
    /// midpoints, divided by the larger direct or interpolated coefficient.
    /// Zero if both are zero; not a continuous spectral error bound.
    pub max_relative_spectral_midpoint_error: f64,
    /// Maximum phase discrepancy at spectral midpoints on a separate cosine
    /// grid, scaled by the larger density or `1/(4 pi)`. Unprobed resonances
    /// are not bounded by this observation.
    pub max_sampled_phase_midpoint_error: f64,
}

impl MieSphere {
    /// Prepare native visible transport from 2–16 explicit vacuum wavelength
    /// knots, starting exactly at 380 nm and ending exactly at 780 nm.
    ///
    /// Every knot is a full Mie solve with its own resolved angular table.
    /// Between wavelengths, coefficients interpolate linearly and normalized
    /// phases interpolate convexly. Comparisons against new full solves at
    /// interval midpoints expose the spectral interpolation error. They are
    /// not a continuous error bound: callers must inspect these diagnostics
    /// and choose wavelengths appropriate to their particles' resonances.
    /// No HG fit or monochromatic angular law replaces the spectral phase.
    pub fn visible_spectrum(
        self,
        cx: &Cx<'_>,
        wavelengths_nm: &[f64],
        budget: MieBudget,
        angular_tolerance: f64,
        max_phase_knots: usize,
    ) -> Result<MieOpticalSpectrum, MieError> {
        prepare_visible_spectrum(
            cx,
            self.host_index_n,
            wavelengths_nm,
            angular_tolerance,
            max_phase_knots,
            |wavelength| self.at_wavelength(cx, wavelength, budget),
        )
    }
}

impl MiePopulation {
    /// Prepare wavelength-dependent mass coefficients and the full
    /// scattering-weighted population phase on the native 380--780 nm grid.
    ///
    /// Uses the same angular refinement and spectral midpoint diagnostics as
    /// [`MieSphere::visible_spectrum`]. The supplied weights represent particle
    /// numbers; grid density remains total dispersed particle mass in kg/m³.
    /// Work limits apply at every spectral knot and diagnostic midpoint.
    pub fn visible_spectrum(
        &self,
        cx: &Cx<'_>,
        wavelengths_nm: &[f64],
        budget: MiePopulationBudget,
        angular_tolerance: f64,
        max_phase_knots: usize,
    ) -> Result<MieOpticalSpectrum, MieError> {
        prepare_visible_spectrum(
            cx,
            self.host_index_n(),
            wavelengths_nm,
            angular_tolerance,
            max_phase_knots,
            |wavelength| self.at_wavelength(cx, wavelength, budget),
        )
    }
}

/// Native spectral preparation shared with wavelength-dependent material
/// queries. Their errors retain their original type instead of becoming an
/// apparent Mie arithmetic failure. The caller owns a constant real host.
pub(crate) fn prepare_visible_spectrum<S, E>(
    cx: &Cx<'_>,
    host_index_n: f64,
    wavelengths_nm: &[f64],
    angular_tolerance: f64,
    max_phase_knots: usize,
    mut solve: impl FnMut(f64) -> Result<S, E>,
) -> Result<MieOpticalSpectrum, E>
where
    S: MieAngularSolution,
    E: From<MieError>,
{
    poll(cx)?;
    if !host_index_n.is_finite()
        || host_index_n <= 0.0
        || !(2..=16).contains(&wavelengths_nm.len())
        || wavelengths_nm.first() != Some(&380.0)
        || wavelengths_nm.last() != Some(&780.0)
        || wavelengths_nm.iter().any(|v| !v.is_finite())
        || wavelengths_nm.windows(2).any(|pair| pair[0] >= pair[1])
        || !angular_tolerance.is_finite()
        || angular_tolerance <= 0.0
        || angular_tolerance > 0.1
        || !(2..=super::MAX_PHASE_KNOTS).contains(&max_phase_knots)
    {
        return Err(MieError::Parameters.into());
    }
    let mut absorption = Vec::with_capacity(wavelengths_nm.len());
    let mut scattering = Vec::with_capacity(wavelengths_nm.len());
    let mut phases = Vec::with_capacity(wavelengths_nm.len());
    let mut angular_reports = Vec::with_capacity(wavelengths_nm.len());
    for &wavelength in wavelengths_nm {
        let solution = solve(wavelength)?;
        let sections = solution.cross_sections();
        absorption.push((wavelength, sections.mass_absorption_m2_kg));
        scattering.push((wavelength, sections.mass_scattering_m2_kg));
        let phase = if sections.scattering_m2 == 0.0 {
            PhaseFunction::Isotropic
        } else {
            let report = solution.phase_table(cx, angular_tolerance, max_phase_knots)?;
            let phase = PhaseFunction::Tabulated(report.phase.clone());
            angular_reports.push(report);
            phase
        };
        phases.push((wavelength, phase));
    }
    let absorption = MassAbsorptionSpectrum::table(&absorption).map_err(|_| MieError::Numerical)?;
    let scattering = MassScatteringSpectrum::table(&scattering).map_err(|_| MieError::Numerical)?;
    let phase = PhaseFunction::spectral_table(&phases).map_err(|_| MieError::Numerical)?;
    let mut max_relative_spectral_midpoint_error = 0.0_f64;
    let mut max_sampled_phase_midpoint_error = 0.0_f64;
    for pair in wavelengths_nm.windows(2) {
        let wavelength = 0.5 * (pair[0] + pair[1]);
        let solution = solve(wavelength)?;
        let sections = solution.cross_sections();
        for (direct, retained) in [
            (
                sections.mass_absorption_m2_kg,
                absorption.at(wavelength).map_err(|_| MieError::Numerical)?,
            ),
            (
                sections.mass_scattering_m2_kg,
                scattering.at(wavelength).map_err(|_| MieError::Numerical)?,
            ),
        ] {
            if direct.max(retained) > 0.0 {
                max_relative_spectral_midpoint_error = max_relative_spectral_midpoint_error
                    .max((direct - retained).abs() / direct.max(retained));
            }
        }
        if sections.scattering_m2 > 0.0 {
            let intervals = (8 * solution.multipoles()).max(32);
            for index in 0..=intervals {
                if index % 32 == 0 {
                    poll(cx)?;
                }
                let cosine =
                    -fs_math::det::cos(core::f64::consts::PI * index as f64 / intervals as f64);
                let direct = solution.eval_cos(cosine)?;
                let retained = phase
                    .eval_at(wavelength, cosine)
                    .map_err(|_| MieError::Numerical)?;
                let scale = direct
                    .max(retained)
                    .max(1.0 / (4.0 * core::f64::consts::PI));
                max_sampled_phase_midpoint_error =
                    max_sampled_phase_midpoint_error.max((direct - retained).abs() / scale);
            }
        }
    }
    poll(cx)?;
    Ok(MieOpticalSpectrum {
        host_index_n,
        absorption,
        scattering,
        phase,
        angular_reports,
        max_relative_spectral_midpoint_error,
        max_sampled_phase_midpoint_error,
    })
}
