//! Source-bounded preparation into the native absorption/scattering owners.

use super::super::{MassAbsorptionSpectrum, MassScatteringSpectrum, PhaseFunction};
use super::{
    MieAngularSolution, MieBudget, MieError, MiePhaseTable, MiePopulation, MiePopulationBudget,
    MieSphere, poll,
};
use fs_blake3::ContentHash;
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

impl MieOpticalSpectrum {
    /// Resolve dispersive particles from their original optical source at every
    /// knot and diagnostic midpoint, then prepare the native transport laws.
    ///
    /// `population_at` returns the same physical radius/number distribution,
    /// particle densities and real host index at each vacuum wavelength in nm.
    /// Only the bins' complex particle indices may vary with wavelength. This
    /// lets existing optical/material owners provide their own supported n/k
    /// laws and preserve their original domain errors; no extrapolation, index
    /// fit, replacement material or gray angular law is introduced here.
    ///
    /// The 2--16 knots define the output support exactly. A nonzero source
    /// identity binds the supplied optical declaration without asserting that
    /// it is measured. Every knot and midpoint retains native Mie work limits;
    /// the returned interpolation discrepancies do not bound unsampled spectral
    /// resonances. Source evaluation and preparation finish before publication.
    #[allow(clippy::too_many_arguments)]
    pub fn from_population_source<E>(
        cx: &Cx<'_>,
        wavelengths_nm: &[f64],
        source_identity: ContentHash,
        budget: MiePopulationBudget,
        angular_tolerance: f64,
        max_phase_knots: usize,
        mut population_at: impl FnMut(f64) -> Result<MiePopulation, E>,
    ) -> Result<Self, E>
    where
        E: From<MieError>,
    {
        poll(cx)?;
        validate_preparation(
            wavelengths_nm,
            Some(source_identity),
            angular_tolerance,
            max_phase_knots,
        )?;
        let first = population_at(wavelengths_nm[0])?;
        if first.bins().len() > budget.max_bins {
            return Err(MieError::Budget.into());
        }
        let host = first.host_index_n();
        let geometry: Vec<_> = first
            .bins()
            .iter()
            .map(|(sphere, fraction)| (sphere.radius_m, sphere.particle_density_kg_m3, *fraction))
            .collect();
        let mut first = Some(first);
        prepare_spectrum(
            cx,
            host,
            wavelengths_nm,
            Some(source_identity),
            angular_tolerance,
            max_phase_knots,
            |wavelength| {
                poll(cx)?;
                let population = match first.take() {
                    Some(value) => value,
                    None => population_at(wavelength)?,
                };
                if population.host_index_n() != host
                    || population.bins().len() != geometry.len()
                    || population.bins().iter().zip(&geometry).any(
                        |((sphere, fraction), &(radius, density, weight))| {
                            sphere.radius_m != radius
                                || sphere.particle_density_kg_m3 != density
                                || *fraction != weight
                        },
                    )
                {
                    return Err(MieError::Parameters.into());
                }
                Ok(population.at_wavelength(cx, wavelength, budget)?)
            },
        )
    }
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

    /// Prepare native absorption, scattering and angular laws over the exact
    /// positive vacuum-wavelength support of an explicitly identified source.
    ///
    /// The 2--16 increasing knots may cover UV, visible or infrared wavelengths.
    /// The declared sphere index is constant over this band; this operation
    /// neither supplies missing dispersion nor confers measurement authority.
    /// The original size/recurrence limits apply at every knot and midpoint.
    /// Spectral/angular discrepancy reports have the same sampled semantics as
    /// [`Self::visible_spectrum`]. No response is extrapolated outside the band.
    #[allow(clippy::too_many_arguments)]
    pub fn source_spectrum(
        self,
        cx: &Cx<'_>,
        wavelengths_nm: &[f64],
        source_identity: ContentHash,
        budget: MieBudget,
        angular_tolerance: f64,
        max_phase_knots: usize,
    ) -> Result<MieOpticalSpectrum, MieError> {
        prepare_spectrum(
            cx,
            self.host_index_n,
            wavelengths_nm,
            Some(source_identity),
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

    /// Prepare an explicitly source-bound UV/visible/IR population spectrum.
    /// Each bin retains its declared constant index, density, radius and number
    /// fraction. Native mass and scattered-power weights remain distinct; a
    /// longer wavelength does not replace the population by an averaged sphere.
    /// See [`MieSphere::source_spectrum`] for support and error semantics.
    #[allow(clippy::too_many_arguments)]
    pub fn source_spectrum(
        &self,
        cx: &Cx<'_>,
        wavelengths_nm: &[f64],
        source_identity: ContentHash,
        budget: MiePopulationBudget,
        angular_tolerance: f64,
        max_phase_knots: usize,
    ) -> Result<MieOpticalSpectrum, MieError> {
        prepare_spectrum(
            cx,
            self.host_index_n(),
            wavelengths_nm,
            Some(source_identity),
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
    solve: impl FnMut(f64) -> Result<S, E>,
) -> Result<MieOpticalSpectrum, E>
where
    S: MieAngularSolution,
    E: From<MieError>,
{
    poll(cx)?;
    if wavelengths_nm.first() != Some(&380.0) || wavelengths_nm.last() != Some(&780.0) {
        return Err(MieError::Parameters.into());
    }
    prepare_spectrum(
        cx,
        host_index_n,
        wavelengths_nm,
        None,
        angular_tolerance,
        max_phase_knots,
        solve,
    )
}

// One preparation and diagnostic owner for visible and source-bounded spectra.
// None preserves the legacy visible table identities and arithmetic.
#[allow(clippy::too_many_arguments)]
pub(crate) fn prepare_spectrum<S, E>(
    cx: &Cx<'_>,
    host_index_n: f64,
    wavelengths_nm: &[f64],
    source_identity: Option<ContentHash>,
    angular_tolerance: f64,
    max_phase_knots: usize,
    mut solve: impl FnMut(f64) -> Result<S, E>,
) -> Result<MieOpticalSpectrum, E>
where
    S: MieAngularSolution,
    E: From<MieError>,
{
    poll(cx)?;
    if !host_index_n.is_finite() || host_index_n <= 0.0 {
        return Err(MieError::Parameters.into());
    }
    validate_preparation(
        wavelengths_nm,
        source_identity,
        angular_tolerance,
        max_phase_knots,
    )?;
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
    let (absorption, scattering, phase) = if let Some(source) = source_identity {
        (
            MassAbsorptionSpectrum::source_table(&absorption, source),
            MassScatteringSpectrum::source_table(&scattering, source),
            PhaseFunction::source_spectral_table(&phases, source),
        )
    } else {
        (
            MassAbsorptionSpectrum::table(&absorption),
            MassScatteringSpectrum::table(&scattering),
            PhaseFunction::spectral_table(&phases),
        )
    };
    let absorption = absorption.map_err(|_| MieError::Numerical)?;
    let scattering = scattering.map_err(|_| MieError::Numerical)?;
    let phase = phase.map_err(|_| MieError::Numerical)?;
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

fn validate_preparation(
    wavelengths_nm: &[f64],
    source_identity: Option<ContentHash>,
    angular_tolerance: f64,
    max_phase_knots: usize,
) -> Result<(), MieError> {
    if !(2..=16).contains(&wavelengths_nm.len())
        || wavelengths_nm.iter().any(|v| !v.is_finite() || *v <= 0.0)
        || wavelengths_nm.windows(2).any(|pair| pair[0] >= pair[1])
        || source_identity == Some(ContentHash([0; 32]))
        || !angular_tolerance.is_finite()
        || angular_tolerance <= 0.0
        || angular_tolerance > 0.1
        || !(2..=super::MAX_PHASE_KNOTS).contains(&max_phase_knots)
    {
        Err(MieError::Parameters)
    } else {
        Ok(())
    }
}
