//! Bounded SI film-width studies for browser and native clients.
//!
//! This adapter calls the native continuous-band film objective and CVaR
//! owners. It does not map dimensionless demo controls to invented SI values,
//! generate source identities, or attach the independent smoke pour animation
//! to a physical study. Results are local Newtonian long-wave model estimates.

use fs_lbm::film::{FilmError, FilmMaterial, LongWaveLimits, NusseltFilm};
use fs_vessel::{
    LaminarFilmBandStudy, LaminarFilmSection, LaminarPourConditions, LaminarRobustError,
    robustify_laminar,
};

/// Maximum sections per browser study.
pub const LAMINAR_MAX_SECTIONS: usize = 64;
/// Maximum authored width candidates.
pub const LAMINAR_MAX_CANDIDATES: usize = 128;
/// Maximum fluid-band points, in addition to one nominal point.
pub const LAMINAR_MAX_FLUIDS: usize = 64;
/// Maximum section/operating-point evaluations, including final profile output.
pub const LAMINAR_MAX_SECTION_EVALUATIONS: usize = 65_536;
/// Number of values in the version-1 output header.
pub const LAMINAR_HEADER_LEN: usize = 18;

/// Typed refusal; no empty vector, clipped input, partial result or NaN verdict.
#[derive(Debug, Clone, PartialEq)]
pub enum LaminarVesselError {
    /// The named packed input has the wrong row shape or is empty.
    Shape(&'static str),
    /// A declared browser size or work limit would be exceeded.
    Budget,
    /// Native candidate, fluid, physical or risk admission failed.
    Study(LaminarRobustError),
    /// Native profile output was not representable.
    Profile { section: usize, source: FilmError },
}

impl core::fmt::Display for LaminarVesselError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Shape(name) => write!(f, "laminar-vessel input shape: {name}"),
            Self::Budget => {
                f.write_str("laminar-vessel section, candidate, fluid or work budget exceeded")
            }
            Self::Study(source) => write!(f, "laminar-vessel study refused: {source}"),
            Self::Profile { section, source } => {
                write!(f, "laminar-vessel section {section}: {source}")
            }
        }
    }
}
impl std::error::Error for LaminarVesselError {}

/// Optimize actual section widths for a physical liquid/flow family.
///
/// Inputs are explicit SI packed rows, copied only after size admission:
/// - `section_rows`: `[width_m, tangential_acceleration_m_s2,
///   normal_acceleration_m_s2]` for each locally planar section;
/// - `width_scales`: positive dimensionless candidates, in tie-breaking order;
/// - `fluid_rows`: `[volume_flow_m3_s, density_kg_m3,
///   kinematic_viscosity_m2_s, surface_tension_n_m]`. The FIRST row is nominal;
///   every subsequent row is one equally weighted CVaR-band operating point;
/// - `settings`: `[minimum_k_rad_m, maximum_k_rad_m, maximum_kh,
///   maximum_Re_kh, beta]`. Every candidate must admit the entire interval.
///
/// Output version 1 has an 18-value header:
/// `[version, sections, candidates, band_points, nominal_index, robust_index,
/// beta, minimum_k, maximum_k, maximum_kh, maximum_Re_kh, nominal_width_scale,
/// robust_width_scale, robust_nominal_growth, nominal_band_cvar,
/// robust_band_cvar, nominal_band_worst, robust_band_worst]`.
/// It is followed by four values per candidate:
/// `[width_scale, nominal_growth_per_s, band_cvar_per_s, band_worst_per_s]`.
/// Finally, each section at the robust width and NOMINAL fluid contributes:
/// `[width_m, thickness_m, mean_velocity_m_s, surface_velocity_m_s, Reynolds,
/// wall_shear_pa, worst_k_rad_m, growth_per_s, max_kh, max_Re_kh]`.
///
/// No stochastic jitter or dense eigensolve is used. Negative growth applies
/// to the declared local long-wave band, not global stability or dribble.
/// Explicit zero surface tension omits capillarity. Values supplied here are
/// authored coefficients; source-backed clients resolve them through the
/// native `MaterialFilmLiquid` owner before crossing this numeric boundary.
pub fn run_laminar_vessel(
    section_rows: &[f64],
    width_scales: &[f64],
    fluid_rows: &[f64],
    settings: &[f64],
) -> Result<Vec<f64>, LaminarVesselError> {
    if section_rows.is_empty() || !section_rows.len().is_multiple_of(3) {
        return Err(LaminarVesselError::Shape(
            "sections require width/acceleration triples",
        ));
    }
    if width_scales.is_empty() {
        return Err(LaminarVesselError::Shape(
            "at least one width candidate is required",
        ));
    }
    if fluid_rows.len() < 8 || !fluid_rows.len().is_multiple_of(4) {
        return Err(LaminarVesselError::Shape(
            "fluids require a nominal quartet and at least one band quartet",
        ));
    }
    if settings.len() != 5 {
        return Err(LaminarVesselError::Shape(
            "settings require k_min/k_max/kh/Re_kh/beta",
        ));
    }
    let sections_count = section_rows.len() / 3;
    let fluid_count = fluid_rows.len() / 4;
    if sections_count > LAMINAR_MAX_SECTIONS
        || width_scales.len() > LAMINAR_MAX_CANDIDATES
        || fluid_count - 1 > LAMINAR_MAX_FLUIDS
        || sections_count * (width_scales.len() * fluid_count + 1) > LAMINAR_MAX_SECTION_EVALUATIONS
    {
        return Err(LaminarVesselError::Budget);
    }
    let sections: Vec<_> = section_rows
        .chunks_exact(3)
        .map(|row| LaminarFilmSection {
            width_m: row[0],
            tangential_acceleration_m_s2: row[1],
            normal_acceleration_m_s2: row[2],
        })
        .collect();
    let fluids: Vec<_> = fluid_rows
        .chunks_exact(4)
        .map(|row| LaminarPourConditions {
            volume_flow_m3_s: row[0],
            material: FilmMaterial {
                density_kg_m3: row[1],
                kinematic_viscosity_m2_s: row[2],
                surface_tension_n_m: row[3],
            },
        })
        .collect();
    let study = LaminarFilmBandStudy {
        sections: &sections,
        wavenumbers_rad_m: [settings[0], settings[1]],
        limits: LongWaveLimits {
            max_kh: settings[2],
            max_re_kh: settings[3],
        },
    };
    let nominal = fluids[0];
    let report = robustify_laminar(&study, width_scales, nominal, &fluids[1..], settings[4])
        .map_err(LaminarVesselError::Study)?;
    let mut output =
        Vec::with_capacity(LAMINAR_HEADER_LEN + 4 * width_scales.len() + 10 * sections_count);
    output.extend_from_slice(&[
        1.0,
        sections_count as f64,
        width_scales.len() as f64,
        (fluid_count - 1) as f64,
        report.nominal_candidate as f64,
        report.robust_candidate as f64,
        settings[4],
        settings[0],
        settings[1],
        settings[2],
        settings[3],
        report.nominal_width_scale,
        report.robust_width_scale,
        report.candidates[report.robust_candidate].nominal_growth_per_s,
        report.nominal_band_cvar_per_s,
        report.robust_band_cvar_per_s,
        report.nominal_worst_growth_per_s,
        report.robust_worst_growth_per_s,
    ]);
    for candidate in &report.candidates {
        output.extend_from_slice(&[
            candidate.width_scale,
            candidate.nominal_growth_per_s,
            candidate.band_cvar_per_s,
            candidate.band_worst_growth_per_s,
        ]);
    }
    for (section, geometry) in sections.iter().enumerate() {
        let width = geometry.width_m * report.robust_width_scale;
        let film = NusseltFilm::new(
            nominal.material,
            nominal.volume_flow_m3_s / width,
            geometry.tangential_acceleration_m_s2,
            geometry.normal_acceleration_m_s2,
        )
        .map_err(|source| LaminarVesselError::Profile { section, source })?;
        let peak = film
            .long_wave_band(study.wavenumbers_rad_m, study.limits)
            .map_err(|source| LaminarVesselError::Profile { section, source })?;
        output.extend_from_slice(&[
            width,
            film.thickness_m(),
            film.mean_velocity_m_s(),
            film.surface_velocity_m_s(),
            film.reynolds(),
            film.wall_shear_pa(),
            peak.wavenumber_rad_m,
            peak.mode.growth_per_s,
            peak.max_kh,
            peak.max_re_kh,
        ]);
    }
    Ok(output)
}

#[cfg(target_arch = "wasm32")]
mod wasm {
    use wasm_bindgen::prelude::*;

    /// Physical SI film study; JavaScript receives a typed-array result or an error.
    #[wasm_bindgen]
    pub fn run_laminar_vessel(
        section_rows: &[f64],
        width_scales: &[f64],
        fluid_rows: &[f64],
        settings: &[f64],
    ) -> Result<Vec<f64>, JsValue> {
        super::run_laminar_vessel(section_rows, width_scales, fluid_rows, settings)
            .map_err(|error| JsValue::from_str(&error.to_string()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const SECTIONS: [f64; 6] = [0.04, 4.0, 3.0, 0.08, 4.0, 3.0];
    const WIDTHS: [f64; 5] = [0.7, 1.0, 1.5, 2.0, 3.0];
    const FLUIDS: [f64; 16] = [
        4e-7, 1000.0, 1e-5, 0.072, 3.2e-7, 1000.0, 1.6e-5, 0.072, 5.2e-7, 1000.0, 0.7e-5, 0.072,
        4.4e-7, 1000.0, 1.2e-5, 0.072,
    ];
    const SETTINGS: [f64; 5] = [2.0, 200.0, 0.3, 0.5, 0.5];

    #[test]
    fn physical_packet_keeps_units_selection_and_native_profile() {
        let result = run_laminar_vessel(&SECTIONS, &WIDTHS, &FLUIDS, &SETTINGS).unwrap();
        assert_eq!(result.len(), LAMINAR_HEADER_LEN + 4 * WIDTHS.len() + 20);
        assert!(result.iter().all(|x| x.is_finite()));
        assert_eq!(&result[..4], &[1.0, 2.0, 5.0, 3.0]);
        let nominal = result[4] as usize;
        let robust = result[5] as usize;
        let rows = &result[LAMINAR_HEADER_LEN..LAMINAR_HEADER_LEN + 4 * WIDTHS.len()];
        assert_eq!(result[14], rows[4 * nominal + 2]);
        assert_eq!(result[15], rows[4 * robust + 2]);
        assert!(result[15] <= result[14]);
        assert_eq!(result[11], WIDTHS[nominal]);
        assert_eq!(result[12], WIDTHS[robust]);
        let first = &result[LAMINAR_HEADER_LEN + 4 * WIDTHS.len()..];
        let q = FLUIDS[0] / first[0];
        assert!((first[1].powi(3) / (3.0 * FLUIDS[2] * q / 4.0) - 1.0).abs() < 3e-13);
        assert!((first[2] * first[1] / q - 1.0).abs() < 3e-13);
        assert_eq!(first[4], q / FLUIDS[2]);
        assert_eq!(
            run_laminar_vessel(&SECTIONS, &WIDTHS, &FLUIDS, &SETTINGS).unwrap(),
            result
        );
    }

    #[test]
    fn unsupported_physics_has_no_packet_or_silent_clamping() {
        let mut settings = SETTINGS;
        settings[1] = 1e9;
        assert!(matches!(
            run_laminar_vessel(&SECTIONS, &WIDTHS, &FLUIDS, &settings),
            Err(LaminarVesselError::Study(LaminarRobustError::Film {
                candidate: 0,
                fluid: None,
                ..
            }))
        ));
        let mut fluid = FLUIDS;
        fluid[6] = f64::NAN;
        assert!(run_laminar_vessel(&SECTIONS, &WIDTHS, &fluid, &SETTINGS).is_err());
        assert!(run_laminar_vessel(&SECTIONS, &[0.0], &FLUIDS, &SETTINGS).is_err());
        settings = SETTINGS;
        settings[4] = 1.0;
        assert!(matches!(
            run_laminar_vessel(&SECTIONS, &WIDTHS, &FLUIDS, &settings),
            Err(LaminarVesselError::Study(LaminarRobustError::Risk(_)))
        ));
    }

    #[test]
    fn shape_and_work_are_admitted_before_numeric_execution() {
        assert!(matches!(
            run_laminar_vessel(&SECTIONS[..5], &WIDTHS, &FLUIDS, &SETTINGS),
            Err(LaminarVesselError::Shape(_))
        ));
        assert!(matches!(
            run_laminar_vessel(&SECTIONS, &WIDTHS, &FLUIDS[..4], &SETTINGS),
            Err(LaminarVesselError::Shape(_))
        ));
        assert_eq!(
            run_laminar_vessel(
                &SECTIONS,
                &[1.0; LAMINAR_MAX_CANDIDATES + 1],
                &FLUIDS,
                &SETTINGS
            ),
            Err(LaminarVesselError::Budget)
        );
        assert_eq!(
            run_laminar_vessel(
                &[1.0; 3 * LAMINAR_MAX_SECTIONS],
                &[1.0; LAMINAR_MAX_CANDIDATES],
                &[1.0; 4 * (LAMINAR_MAX_FLUIDS + 1)],
                &SETTINGS
            ),
            Err(LaminarVesselError::Budget)
        );
    }
}
