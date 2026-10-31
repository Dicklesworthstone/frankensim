//! Material-derived longitudinal and transverse motion of a tensioned ribbon.
//!
//! A cell is a quadrature strip of parallel, noninteracting straight fibers,
//! NOT one fictitious thick fiber. Its mass, EA, EI, loss and pre-tension all
//! scale with its share of the original fiber count. Thus width refinement
//! changes neither total mass nor compliance under a uniform line traction.
//! Each strip retains the sine Galerkin modes of the fixed-end axial wave
//! equation and the pinned, pre-tensioned Euler--Bernoulli equation. There is
//! no selected musical frequency, damping ratio, output gain or sound target.
//!
//! The physical reason for retaining hair stretch and finite width in bowing
//! is documented by Pitteroff & Woodhouse, Acta Acustica 84 (1998), 744--757,
//! https://www.repository.cam.ac.uk/handle/1810/244810 (review, section 1.1).
//! The base rung is the straight-fiber, rigid-end, small-motion limit. Explicit
//! `with_axial_transverse_strain` adds the existing moderate-slope extensional
//! law to the SAME contact solver, so normal deflection changes fiber tension
//! and axial waves. Both rungs omit stick camber/flexure, fiber-to-fiber contact,
//! large rotations, air drag and thermal rosin. The separate explicit
//! `with_material_relaxation` binding retains each source Young-modulus arm
//! in both axial and bending response, through the same coupled contact owner.
//! With nonlinear strain selected, axial memory follows that same complete
//! strain and retains the geometric mean; bending memory remains separate.
//! Fixed-state Kelvin--Voigt coefficients require explicit source data;
//! missing loss is never replaced with an invented Q. Refine BOTH the span
//! modes and ribbon cells; the sample-rate owner refuses unresolved modes.

use crate::modal_acoustic_time::{
    MAX_TIME_DOMAIN_ACOUSTIC_MODES, ModalAcousticMode, ModalAcousticTimeBudget,
    ModalAcousticTimeError, ModalAcousticTimeModel,
};
use fs_exec::CancelGate;
use fs_material::section::{ConcentricKelvinVoigtSection, axial::ConcentricAxialViscosity};
use fs_material::state_point::MaterialStatePointError;
use fs_math::{c64::C64, det};

/// Conservative coupling of fiber extension and normal deflection.
pub mod strain;
/// Source-bound hereditary Young response in both retained fiber directions.
pub mod relaxation;

/// Dimensions and fixed pre-tension of the original physical fiber bundle.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiberRibbonGeometry {
    /// Unsupported length of every straight fiber [m].
    pub length_m: f64,
    /// Total ribbon width [m], along the contacted object's crossing direction.
    pub width_m: f64,
    /// Actual number of equivalent fibers, independent of numerical cells.
    pub fiber_count: usize,
    /// Total tensile preload across all original fibers [N].
    pub tension_n: f64,
}

/// Explicit independent discretizations of width and fiber span.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct FiberRibbonResolution {
    /// Odd midpoint quadrature count in 1..=63, matching the contact owner.
    pub cells: usize,
    /// Consecutive longitudinal sine modes starting at one.
    pub axial_modes: usize,
    /// Consecutive normal bending/tension sine modes starting at one.
    pub transverse_modes: usize,
}

/// Projection onto the two independently retained displacement directions.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FiberDirection {
    /// Along the fiber (bowing tangent).
    Axial,
    /// Across the fiber (contact normal).
    Normal,
}

/// Admission retains the physical or numerical owner's error.
#[derive(Debug)]
pub enum FiberRibbonError {
    Invalid(&'static str),
    Material(MaterialStatePointError),
    Modal(ModalAcousticTimeError),
    Strain(fs_nlmodal::NlModalError),
    Cancelled,
}
impl core::fmt::Display for FiberRibbonError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(what) => write!(f, "fiber ribbon: {what}"),
            Self::Material(error) => write!(f, "fiber ribbon material: {error:?}"),
            Self::Modal(error) => write!(f, "fiber ribbon modes: {error}"),
            Self::Strain(error) => write!(f, "fiber ribbon strain: {error}"),
            Self::Cancelled => write!(f, "fiber ribbon preparation cancelled"),
        }
    }
}
impl std::error::Error for FiberRibbonError {}

/// Original section sources and a convergent basis, reusable by any contact.
#[derive(Clone, Debug)]
pub struct TensionedFiberRibbon {
    geometry: FiberRibbonGeometry,
    resolution: FiberRibbonResolution,
    axial: ConcentricAxialViscosity,
    bending: ConcentricKelvinVoigtSection,
    cell: ModalAcousticTimeModel,
    cell_line_density_kg_m: f64,
    strain: Option<strain::FiberStrain>,
    relaxation: Option<relaxation::FiberRibbonRelaxation>,
}

impl TensionedFiberRibbon {
    /// Derive all modes from a SINGLE fiber section and its explicit loss laws.
    /// The axial and bending laws must retain exactly the same section/source.
    /// All selected modes must fit both material support and the physical clock.
    pub fn try_new(
        geometry: FiberRibbonGeometry,
        resolution: FiberRibbonResolution,
        axial: ConcentricAxialViscosity,
        bending: ConcentricKelvinVoigtSection,
        sample_rate_hz: u32,
        budget: ModalAcousticTimeBudget,
        gate: &CancelGate,
    ) -> Result<Self, FiberRibbonError> {
        let invalid = FiberRibbonError::Invalid;
        if gate.is_requested() {
            return Err(FiberRibbonError::Cancelled);
        }
        if !geometry.length_m.is_finite()
            || geometry.length_m <= 0.0
            || !geometry.width_m.is_finite()
            || geometry.width_m <= 0.0
            || !geometry.tension_n.is_finite()
            || geometry.tension_n <= 0.0
            || geometry.fiber_count == 0
            || geometry.fiber_count > (1_usize << 24)
            || !(1..=63).contains(&resolution.cells)
            || resolution.cells % 2 == 0
            || resolution.axial_modes == 0
            || resolution.transverse_modes == 0
            || resolution
                .axial_modes
                .checked_add(resolution.transverse_modes)
                .is_none_or(|n| n > MAX_TIME_DOMAIN_ACOUSTIC_MODES)
        {
            return Err(invalid(
                "positive finite geometry/tension and bounded nonempty resolutions required",
            ));
        }
        if axial.elastic() != bending.elastic() {
            return Err(invalid(
                "axial and bending laws must use the same original fiber section",
            ));
        }
        let section = axial.elastic();
        let per_cell = geometry.fiber_count as f64 / resolution.cells as f64;
        let line_density = section.line_density_kg_m();
        let cell_line_density_kg_m = line_density * per_cell;
        let normalization = (0.5 * cell_line_density_kg_m * geometry.length_m).sqrt();
        let total_mass = line_density * geometry.length_m * geometry.fiber_count as f64;
        if !cell_line_density_kg_m.is_finite()
            || cell_line_density_kg_m <= 0.0
            || !normalization.is_finite()
            || normalization <= 0.0
            || !total_mass.is_finite()
            || total_mass <= 0.0
        {
            return Err(invalid("fiber-cell mass is unrepresentable"));
        }
        let mut modes = Vec::with_capacity(resolution.axial_modes + resolution.transverse_modes);
        for mode in 1..=resolution.axial_modes {
            if gate.is_requested() {
                return Err(FiberRibbonError::Cancelled);
            }
            let k = mode as f64 * core::f64::consts::PI / geometry.length_m;
            let omega = (section.axial_stiffness_n() / line_density).sqrt() * k;
            let damping = axial
                .longitudinal_damping_ratio(omega)
                .map_err(FiberRibbonError::Material)?;
            modes.push(ModalAcousticMode {
                angular_frequency_rad_s: omega,
                damping_ratio: damping,
                pressure_per_modal_velocity: C64::ZERO,
            });
        }
        let (lo, hi) = bending.omega_band_rad_s();
        for mode in 1..=resolution.transverse_modes {
            if gate.is_requested() {
                return Err(FiberRibbonError::Cancelled);
            }
            let k = mode as f64 * core::f64::consts::PI / geometry.length_m;
            let omega = fs_nlmodal::prestressed_beam_omega(
                geometry.length_m,
                geometry.tension_n / geometry.fiber_count as f64,
                line_density,
                section.bending_stiffness_n_m2(),
                mode,
            );
            if !omega.is_finite() || omega <= 0.0 || omega < lo || omega > hi {
                return Err(invalid(
                    "transverse mode exceeds fiber material/loss frequency support",
                ));
            }
            let damping =
                0.5 * bending.viscous_stiffness_n_m2_s() * k.powi(4) / line_density / omega;
            if bending.viscous_stiffness_n_m2_s() > 0.0 && damping == 0.0 {
                return Err(invalid(
                    "positive bending loss is below representable modal damping",
                ));
            }
            modes.push(ModalAcousticMode {
                angular_frequency_rad_s: omega,
                damping_ratio: damping,
                pressure_per_modal_velocity: C64::ZERO,
            });
        }
        let cell = ModalAcousticTimeModel::try_new(sample_rate_hz, modes, budget)
            .map_err(FiberRibbonError::Modal)?;
        Ok(Self {
            geometry,
            resolution,
            axial,
            bending,
            cell,
            cell_line_density_kg_m,
            strain: None,
            relaxation: None,
        })
    }

    pub const fn geometry(&self) -> FiberRibbonGeometry {
        self.geometry
    }
    pub const fn resolution(&self) -> FiberRibbonResolution {
        self.resolution
    }
    pub const fn axial_material(&self) -> &ConcentricAxialViscosity {
        &self.axial
    }
    pub const fn bending_material(&self) -> &ConcentricKelvinVoigtSection {
        &self.bending
    }
    /// Prototype in the SAME mass-normalized basis as every returned contact row.
    pub const fn cell_model(&self) -> &ModalAcousticTimeModel {
        &self.cell
    }
    pub const fn cell_line_density_kg_m(&self) -> f64 {
        self.cell_line_density_kg_m
    }
    /// Mass of the original fibers, never multiplied by numerical cell count twice.
    pub fn mass_kg(&self) -> f64 {
        self.axial.elastic().line_density_kg_m()
            * self.geometry.length_m
            * self.geometry.fiber_count as f64
    }
    /// Independent strip states; contact coupling occurs in the shared network.
    pub fn cell_models(&self) -> Vec<ModalAcousticTimeModel> {
        vec![self.cell.clone(); self.resolution.cells]
    }
    /// Mass projection of rigid axial support translation, in this cell's
    /// mass-normalized basis. Relative motion receives `-Gamma_n X''`, where
    /// `Gamma_n = integral(mu phi_n ds) = sqrt(2 mu L) (1-(-1)^n)/(n pi)`.
    /// Normal modes have zero axial participation; cell refinement partitions
    /// the original fiber mass instead of adding any support inertia.
    pub fn axial_translation_participation(&self) -> Vec<f64> {
        self.translation_participation(FiberDirection::Axial)
    }
    /// Mass projection of rigid support translation in either retained direction.
    /// The normal projection uses the SAME original fiber mass as the axial one;
    /// it is not a lumped bow mass. Only odd fixed-end sine modes participate.
    pub fn translation_participation(&self, direction: FiberDirection) -> Vec<f64> {
        let scale = 2.0 * (0.5 * self.cell_line_density_kg_m * self.geometry.length_m).sqrt();
        let mut participation = vec![0.0; self.cell.modes().len()];
        let (start, count) = match direction {
            FiberDirection::Axial => (0, self.resolution.axial_modes),
            FiberDirection::Normal => (self.resolution.axial_modes, self.resolution.transverse_modes),
        };
        for mode in (1..=count).step_by(2) {
            participation[start + mode - 1] = 2.0 * scale / (mode as f64 * core::f64::consts::PI);
        }
        participation
    }
    /// Material-point displacement/force row at a fixed interior fiber station.
    /// Both directions refer to this SAME component: axial modes precede normal.
    pub fn displacement_map(
        &self,
        station_fraction: f64,
        direction: FiberDirection,
    ) -> Result<Vec<f64>, FiberRibbonError> {
        if !station_fraction.is_finite() || !(0.0 < station_fraction && station_fraction < 1.0) {
            return Err(FiberRibbonError::Invalid(
                "fiber contact must lie strictly inside its supported span",
            ));
        }
        let scale = (0.5 * self.cell_line_density_kg_m * self.geometry.length_m)
            .sqrt()
            .recip();
        let mut shape = vec![0.0; self.cell.modes().len()];
        let (start, count) = match direction {
            FiberDirection::Axial => (0, self.resolution.axial_modes),
            FiberDirection::Normal => (
                self.resolution.axial_modes,
                self.resolution.transverse_modes,
            ),
        };
        for mode in 1..=count {
            shape[start + mode - 1] =
                det::sin(mode as f64 * core::f64::consts::PI * station_fraction) * scale;
        }
        if shape.iter().any(|x| !x.is_finite()) {
            return Err(FiberRibbonError::Invalid(
                "fiber displacement map is unrepresentable",
            ));
        }
        Ok(shape)
    }
}
