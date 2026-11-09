//! Conservative long-wave evolution of an original laminar film.
//!
//! `h_t + q_x = 0`, `p = rho*g_normal*h - sigma*h_xx`, and
//! `tau(y) = tau_surface + (rho*g_tangent-p_x)*(h-y)`. Each face uses
//! the existing Newtonian film or original prepared Carreau slit; no fitted
//! mobility or superposition of nonlinear pressure/shear profiles is used.
//! The Newtonian reduction agrees with Oron, Davis & Bankoff, Rev. Mod. Phys.
//! 69 (1997), sections II.B-D, <https://doi.org/10.1103/RevModPhys.69.931>.
//!
//! Backward Euler advances actual cell volumes. Its quadratic gravitational
//! and capillary energy defect remains numerical, separate from viscous heat.
//! The same accepted face profile supplies flux, traction work and heat.
//! This is a fully wet, quasi-steady, small-slope reduction on a stationary
//! planar wall. It omits momentum storage, rupture/contact lines, evaporation,
//! thermal evolution and transverse structure. Frozen temperatures may vary
//! by face; density and surface tension remain uniform. Prescribed traction
//! can carry separately authored gas/Marangoni stress, without inferring it.

use std::cell::{Cell, RefCell};

use fs_material::MaterialError;
use fs_material::fluid_transport::rheology::ThermalCarreauYasuda;
use fs_material::fluid_transport::rheology::slit::{CarreauSlitControls, LaminarCarreauSlit};
use fs_math::dd::Dd;
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::{NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport};

use super::{FilmError, FilmMaterial, ShearedFilm};

/// Periodic cells, or a closed film with natural zero-slope ends and no flux.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ThinFilmBoundary {
    Periodic,
    Closed,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThinFilmGeometry {
    pub spacing_m: f64,
    /// Nonnegative acceleration into the supporting wall.
    pub normal_acceleration_m_s2: f64,
    pub boundary: ThinFilmBoundary,
}

/// Authored physical admission; these screens are not calibrated error bounds.
#[derive(Clone, Copy, Debug)]
pub struct ThinFilmLimits {
    /// A strictly wet domain. Drying refuses rather than clipping mass.
    pub minimum_thickness_m: f64,
    pub maximum_thickness_m: f64,
    /// Physical streamwise variation scale, independent of grid spacing.
    pub characteristic_length_m: f64,
    pub maximum_height_to_length: f64,
    pub maximum_slope: f64,
    /// rho*max|u|*h^2/(minimum local mu * characteristic_length).
    /// Uses profile extrema, so zero-net-flux counterflow is not treated as rest.
    /// Must be below one; a smaller authored budget controls model admission.
    pub maximum_reduced_reynolds: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct ThinFilmControls {
    pub maximum_newton_iterations: usize,
    pub maximum_endpoint_evaluations: usize,
    /// Whole-call original-profile budget, including discarded Godunov states
    /// and every finite-difference/Krylov probe. Native slit work is also bounded.
    pub maximum_profile_evaluations: u64,
    pub height_scale_m: f64,
    pub finite_difference_step: f64,
    pub cell_tolerance_m: f64,
    pub volume_tolerance_m2: f64,
    pub energy_tolerance_j_per_m: f64,
    pub solver: NewtonKrylovConfig,
    pub slit: CarreauSlitControls,
}

impl Default for ThinFilmControls {
    fn default() -> Self {
        Self {
            maximum_newton_iterations: 32,
            maximum_endpoint_evaluations: 2048,
            maximum_profile_evaluations: 16_777_216,
            height_scale_m: 1e-3,
            finite_difference_step: 1e-5,
            cell_tolerance_m: 1e-10,
            volume_tolerance_m2: 1e-10,
            energy_tolerance_j_per_m: 1e-8,
            solver: NewtonKrylovConfig {
                forcing_minimum: 1e-4,
                ..NewtonKrylovConfig::default()
            },
            slit: CarreauSlitControls::default(),
        }
    }
}

#[derive(Debug)]
pub enum ThinFilmError {
    Invalid(&'static str),
    Film(FilmError),
    Material(MaterialError),
    Solver(NewtonError),
    NotConverged(NewtonReport),
    OutsideRegime(&'static str),
    Cancelled,
    EndpointBudget,
    ProfileBudget,
    CellResidual {
        residual_m: f64,
        tolerance_m: f64,
    },
    VolumeResidual {
        residual_m2: f64,
        tolerance_m2: f64,
    },
    EnergyResidual {
        residual_j_per_m: f64,
        tolerance_j_per_m: f64,
    },
}
impl From<FilmError> for ThinFilmError {
    fn from(value: FilmError) -> Self {
        Self::Film(value)
    }
}
impl From<MaterialError> for ThinFilmError {
    fn from(value: MaterialError) -> Self {
        Self::Material(value)
    }
}
impl core::fmt::Display for ThinFilmError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "thin-film evolution refused: {self:?}")
    }
}
impl std::error::Error for ThinFilmError {}

#[derive(Clone, Debug)]
enum Source {
    Newtonian(FilmMaterial),
    Carreau(Vec<LaminarCarreauSlit>),
}

/// Original steady profile at the accepted finite-volume face.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThinFilmFace {
    pub from: usize,
    pub to: usize,
    /// Godunov state selected from the two cells or the physical turning point.
    pub thickness_m: f64,
    pub flow_per_width_m2_s: f64,
    pub surface_velocity_m_s: f64,
    pub wall_shear_pa: f64,
    pub force_density_pa_per_m: f64,
    pub surface_shear_pa: f64,
    pub viscous_dissipation_w_m2: f64,
    pub peak_speed_m_s: f64,
    pub minimum_viscosity_pa_s: f64,
    pub maximum_shear_rate_per_s: f64,
}

/// Physical and numerical energy accounts per unit out-of-plane width.
#[derive(Debug)]
pub struct ThinFilmStep {
    pub time_s: f64,
    pub duration_s: f64,
    pub faces: Vec<ThinFilmFace>,
    pub volume_change_m2: f64,
    pub maximum_cell_residual_m: f64,
    pub stored_energy_change_j_per_m: f64,
    /// E(h_new-h_old), nonnegative for the admitted gravity/capillarity.
    /// This is numerical damping and must never be deposited as physical heat.
    pub backward_euler_energy_defect_j_per_m: f64,
    pub viscous_heat_j_per_m: f64,
    pub body_work_j_per_m: f64,
    pub surface_work_j_per_m: f64,
    /// Delta E + BE defect + physical heat - body work - surface work.
    pub energy_residual_j_per_m: f64,
    pub maximum_slope: f64,
    pub maximum_height_to_length: f64,
    pub maximum_reduced_reynolds: f64,
    pub endpoint_evaluations: usize,
    pub profile_evaluations: u64,
    pub nonlinear: NewtonReport,
}

/// Uniform one-dimensional native film volumes with frozen original sources.
#[derive(Clone, Debug)]
pub struct ThinFilmLine {
    geometry: ThinFilmGeometry,
    density: f64,
    tension: f64,
    source: Source,
    heights: Vec<f64>,
    time: f64,
}

struct Evaluation {
    heights: Vec<f64>,
    residual: Vec<f64>,
    faces: Vec<ThinFilmFace>,
    slope: f64,
    aspect: f64,
    reynolds: f64,
}

impl ThinFilmLine {
    pub fn newtonian(
        material: FilmMaterial,
        geometry: ThinFilmGeometry,
        heights: Vec<f64>,
        time_s: f64,
    ) -> Result<Self, ThinFilmError> {
        Self::admit_geometry(geometry, &heights, time_s)?;
        // Native admission owns the actual Newtonian material coefficients.
        ShearedFilm::new(material, heights[0], 0.0, 0.0)?;
        Ok(Self {
            geometry,
            density: material.density_kg_m3,
            tension: material.surface_tension_n_m,
            source: Source::Newtonian(material),
            heights,
            time: time_s,
        })
    }

    /// Prepare each original face temperature once; all later thickness trials
    /// reuse its native constitutive preparation through with_gap.
    #[allow(clippy::too_many_arguments)]
    pub fn carreau(
        law: ThermalCarreauYasuda,
        face_temperatures_k: &[f64],
        density_kg_m3: f64,
        surface_tension_n_m: f64,
        maximum_shear_rate_per_s: f64,
        geometry: ThinFilmGeometry,
        heights: Vec<f64>,
        time_s: f64,
    ) -> Result<Self, ThinFilmError> {
        Self::admit_geometry(geometry, &heights, time_s)?;
        let count = faces(geometry.boundary, heights.len());
        if !positive(density_kg_m3)
            || !nonnegative(surface_tension_n_m)
            || face_temperatures_k.len() != count
        {
            return Err(ThinFilmError::Invalid(
                "original density, surface tension and one temperature per face required",
            ));
        }
        let source = face_temperatures_k
            .iter()
            .map(|&temperature| {
                LaminarCarreauSlit::new(
                    law,
                    temperature,
                    1.0,
                    heights[0],
                    1.0,
                    maximum_shear_rate_per_s,
                )
            })
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Self {
            geometry,
            density: density_kg_m3,
            tension: surface_tension_n_m,
            source: Source::Carreau(source),
            heights,
            time: time_s,
        })
    }

    fn admit_geometry(g: ThinFilmGeometry, h: &[f64], time: f64) -> Result<(), ThinFilmError> {
        if !(3..=65_536).contains(&h.len())
            || h.iter().any(|&x| !positive(x))
            || !positive(g.spacing_m)
            || !nonnegative(g.normal_acceleration_m_s2)
            || !time.is_finite()
            || time < 0.0
            || !positive(g.spacing_m * h.len() as f64)
        {
            return Err(ThinFilmError::Invalid(
                "bounded positive uniform film cells, geometry and finite time required",
            ));
        }
        Ok(())
    }

    pub fn heights_m(&self) -> &[f64] {
        &self.heights
    }
    pub const fn accepted_time_s(&self) -> f64 {
        self.time
    }
    pub const fn geometry(&self) -> ThinFilmGeometry {
        self.geometry
    }
    pub const fn density_kg_m3(&self) -> f64 {
        self.density
    }
    pub const fn surface_tension_n_m(&self) -> f64 {
        self.tension
    }
    pub fn face_count(&self) -> usize {
        faces(self.geometry.boundary, self.heights.len())
    }

    /// Conservative backward Euler, with original nonlinear face profiles.
    ///
    /// A face numbered i joins cell i to i+1; the final periodic face wraps to
    /// zero. Closed ends have no additional face or imposed flux. The caller
    /// supplies one held surface traction per face and a complete height guess.
    /// Every candidate must remain strictly wet and source-supported; no clip,
    /// mass repair, apparent-viscosity substitution or post-solve rescaling is
    /// performed. Failed proposals preserve heights and accepted time.
    ///
    /// First-order time and Godunov transport require refinement for accuracy.
    /// The long-wave scales must be physical, not chosen to relax admission.
    /// The static native profile omits fast forcing/startup momentum storage.
    #[allow(clippy::too_many_arguments)]
    pub fn advance(
        &mut self,
        target_s: f64,
        acceleration_m_s2: f64,
        surface_shear_pa: &[f64],
        limits: ThinFilmLimits,
        controls: ThinFilmControls,
        height_guess_m: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<ThinFilmStep, ThinFilmError> {
        let checkpoint = RefCell::new(cancelled);
        let stopped = Cell::new(false);
        let check = || {
            if stopped.get() || (checkpoint.borrow_mut())() {
                stopped.set(true);
                Err(ThinFilmError::Cancelled)
            } else {
                Ok(())
            }
        };
        check()?;
        let dt = target_s - self.time;
        admit(limits, controls)?;
        if !positive(dt)
            || !target_s.is_finite()
            || !acceleration_m_s2.is_finite()
            || surface_shear_pa.len() != self.face_count()
            || surface_shear_pa.iter().any(|x| !x.is_finite())
            || height_guess_m.len() != self.heights.len()
        {
            return Err(ThinFilmError::Invalid(
                "finite interval, complete face traction and height guess required",
            ));
        }
        self.admit_heights(&self.heights, limits, &check)?;
        let body = finite(self.density * acceleration_m_s2)?;
        let gravity = finite(self.density * self.geometry.normal_acceleration_m_s2)?;
        let count = self.heights.len();
        let evaluations = Cell::new(0usize);
        let profiles = Cell::new(0u64);
        let evaluate = |x: &[f64]| -> Result<Evaluation, ThinFilmError> {
            check()?;
            if evaluations.get() >= controls.maximum_endpoint_evaluations {
                return Err(ThinFilmError::EndpointBudget);
            }
            evaluations.set(evaluations.get() + 1);
            if x.len() != count {
                return Err(ThinFilmError::Invalid(
                    "complete normalized film coordinates required",
                ));
            }
            let h = x
                .iter()
                .zip(&self.heights)
                .map(|(&x, &old)| finite(old + x * controls.height_scale_m))
                .collect::<Result<Vec<_>, _>>()?;
            let (slope, aspect) = self.admit_heights(&h, limits, &check)?;
            let dx = self.geometry.spacing_m;
            // Remove a constant pressure datum before differencing. It has
            // no effect on any flux; the final energy uses actual heights.
            let mut pressure = h
                .iter()
                .map(|&height| finite(gravity * (height - self.heights[0])))
                .collect::<Result<Vec<_>, _>>()?;
            for i in 0..self.face_count() {
                let j = (i + 1) % count;
                let correction = finite(self.tension * ((h[j] - h[i]) / dx) / dx)?;
                pressure[i] = finite(pressure[i] - correction)?;
                pressure[j] = finite(pressure[j] + correction)?;
            }
            let mut residual: Vec<_> = h.iter().zip(&self.heights).map(|(a, b)| a - b).collect();
            let mut selected = Vec::with_capacity(self.face_count());
            let mut reynolds = 0.0_f64;
            for (i, &traction) in surface_shear_pa.iter().enumerate() {
                check()?;
                let j = (i + 1) % count;
                let force = finite(body - (pressure[j] - pressure[i]) / dx)?;
                let profile = |height| {
                    if profiles.get() >= controls.maximum_profile_evaluations {
                        return Err(ThinFilmError::ProfileBudget);
                    }
                    profiles.set(profiles.get() + 1);
                    self.profile(i, j, height, force, traction, controls.slit, &check)
                };
                let mut chosen = profile(h[i])?;
                let minimum = h[i] <= h[j];
                let prefer = |a: &ThinFilmFace, b: &ThinFilmFace| {
                    if minimum {
                        a.flow_per_width_m2_s < b.flow_per_width_m2_s
                    } else {
                        a.flow_per_width_m2_s > b.flow_per_width_m2_s
                    }
                };
                if h[i] != h[j] {
                    let other = profile(h[j])?;
                    if prefer(&other, &chosen) {
                        chosen = other;
                    }
                    // dq/dh = h*gamma(tau_surface + force*h). The original
                    // monotone stress inverse has exactly this possible turn.
                    if force != 0.0 {
                        let turning = -traction / force;
                        if turning > h[i].min(h[j]) && turning < h[i].max(h[j]) {
                            let other = profile(turning)?;
                            if prefer(&other, &chosen) {
                                chosen = other;
                            }
                        }
                    }
                }
                let re = finite(
                    self.density * chosen.peak_speed_m_s * chosen.thickness_m
                        / chosen.minimum_viscosity_pa_s
                        * (chosen.thickness_m / limits.characteristic_length_m),
                )?;
                if re > limits.maximum_reduced_reynolds {
                    return Err(ThinFilmError::OutsideRegime(
                        "reduced inertia exceeds long-wave budget",
                    ));
                }
                reynolds = reynolds.max(re);
                let transferred = finite(dt * chosen.flow_per_width_m2_s / dx)?;
                residual[i] = finite(residual[i] + transferred)?;
                residual[j] = finite(residual[j] - transferred)?;
                selected.push(chosen);
            }
            Ok(Evaluation {
                heights: h,
                residual,
                faces: selected,
                slope,
                aspect,
                reynolds,
            })
        };
        let difference = FallibleFiniteDifference::new(
            count,
            controls.finite_difference_step,
            |x: &[f64]| -> Result<Vec<f64>, ThinFilmError> {
                Ok(evaluate(x)?
                    .residual
                    .into_iter()
                    .map(|v| v / controls.height_scale_m)
                    .collect())
            },
        )
        .map_err(difference_error)?;
        let initial = height_guess_m
            .iter()
            .zip(&self.heights)
            .map(|(&h, &old)| (h - old) / controls.height_scale_m)
            .collect();
        let initialized = NewtonKrylovState::new(&difference, initial, controls.solver);
        if let Some(error) = difference.take_failure() {
            return Err(difference_error(error));
        }
        let mut solver = initialized.map_err(ThinFilmError::Solver)?;
        let result =
            solver.run_controlled(&difference, controls.maximum_newton_iterations, &mut || {
                difference.has_failed() || check().is_err()
            });
        if let Some(error) = difference.take_failure() {
            return Err(difference_error(error));
        }
        check()?;
        if !result.report.converged {
            return Err(ThinFilmError::NotConverged(result.report));
        }
        let accepted = evaluate(&solver.x)?;
        let maximum_cell_residual_m = accepted
            .residual
            .iter()
            .fold(0.0_f64, |a, b| a.max(b.abs()));
        if maximum_cell_residual_m > controls.cell_tolerance_m {
            return Err(ThinFilmError::CellResidual {
                residual_m: maximum_cell_residual_m,
                tolerance_m: controls.cell_tolerance_m,
            });
        }
        let dx = self.geometry.spacing_m;
        let delta: Vec<_> = accepted
            .heights
            .iter()
            .zip(&self.heights)
            .map(|(a, b)| a - b)
            .collect();
        let volume_change_m2 = total(delta.iter().map(|&h| h * dx))?;
        if volume_change_m2.abs() > controls.volume_tolerance_m2 {
            return Err(ThinFilmError::VolumeResidual {
                residual_m2: volume_change_m2,
                tolerance_m2: controls.volume_tolerance_m2,
            });
        }
        let mut change = Dd::from_f64(0.0);
        let mut defect = Dd::from_f64(0.0);
        for i in 0..count {
            if i % 128 == 0 {
                check()?;
            }
            change = change
                + Dd::from_f64(finite(
                    0.5 * gravity * dx * delta[i] * (accepted.heights[i] + self.heights[i]),
                )?);
            defect = defect + Dd::from_f64(finite(0.5 * gravity * dx * delta[i] * delta[i])?);
        }
        for i in 0..self.face_count() {
            let j = (i + 1) % count;
            let old = self.heights[j] - self.heights[i];
            let new = accepted.heights[j] - accepted.heights[i];
            let d = delta[j] - delta[i];
            change = change + Dd::from_f64(finite(0.5 * self.tension / dx * d * (new + old))?);
            defect = defect + Dd::from_f64(finite(0.5 * self.tension / dx * d * d)?);
        }
        let stored_energy_change_j_per_m = finite(change.to_f64())?;
        let backward_euler_energy_defect_j_per_m = finite(defect.to_f64())?;
        let viscous_heat_j_per_m = total(
            accepted
                .faces
                .iter()
                .map(|f| dt * dx * f.viscous_dissipation_w_m2),
        )?;
        let body_work_j_per_m = total(
            accepted
                .faces
                .iter()
                .map(|f| dt * dx * body * f.flow_per_width_m2_s),
        )?;
        let surface_work_j_per_m = total(
            accepted
                .faces
                .iter()
                .map(|f| dt * dx * f.surface_shear_pa * f.surface_velocity_m_s),
        )?;
        let energy_residual_j_per_m = total([
            stored_energy_change_j_per_m,
            backward_euler_energy_defect_j_per_m,
            viscous_heat_j_per_m,
            -body_work_j_per_m,
            -surface_work_j_per_m,
        ])?;
        if energy_residual_j_per_m.abs() > controls.energy_tolerance_j_per_m {
            return Err(ThinFilmError::EnergyResidual {
                residual_j_per_m: energy_residual_j_per_m,
                tolerance_j_per_m: controls.energy_tolerance_j_per_m,
            });
        }
        check()?;
        let report = ThinFilmStep {
            time_s: target_s,
            duration_s: dt,
            faces: accepted.faces,
            volume_change_m2,
            maximum_cell_residual_m,
            stored_energy_change_j_per_m,
            backward_euler_energy_defect_j_per_m,
            viscous_heat_j_per_m,
            body_work_j_per_m,
            surface_work_j_per_m,
            energy_residual_j_per_m,
            maximum_slope: accepted.slope,
            maximum_height_to_length: accepted.aspect,
            maximum_reduced_reynolds: accepted.reynolds,
            endpoint_evaluations: evaluations.get(),
            profile_evaluations: profiles.get(),
            nonlinear: result.report,
        };
        self.heights = accepted.heights;
        self.time = target_s;
        Ok(report)
    }

    fn admit_heights(
        &self,
        h: &[f64],
        limits: ThinFilmLimits,
        check: &dyn Fn() -> Result<(), ThinFilmError>,
    ) -> Result<(f64, f64), ThinFilmError> {
        let mut aspect = 0.0_f64;
        for (i, &height) in h.iter().enumerate() {
            if i % 128 == 0 {
                check()?;
            }
            if !height.is_finite()
                || height < limits.minimum_thickness_m
                || height > limits.maximum_thickness_m
            {
                return Err(ThinFilmError::OutsideRegime(
                    "film must remain inside the strictly wet thickness interval",
                ));
            }
            aspect = aspect.max(finite(height / limits.characteristic_length_m)?);
        }
        let mut slope = 0.0_f64;
        for i in 0..self.face_count() {
            slope = slope.max(finite(
                (h[(i + 1) % h.len()] - h[i]).abs() / self.geometry.spacing_m,
            )?);
        }
        if aspect > limits.maximum_height_to_length || slope > limits.maximum_slope {
            return Err(ThinFilmError::OutsideRegime(
                "thin-film height or slope exceeds its long-wave budget",
            ));
        }
        Ok((slope, aspect))
    }

    #[allow(clippy::too_many_arguments)]
    fn profile(
        &self,
        from: usize,
        to: usize,
        height: f64,
        force: f64,
        traction: f64,
        controls: CarreauSlitControls,
        check: &dyn Fn() -> Result<(), ThinFilmError>,
    ) -> Result<ThinFilmFace, ThinFilmError> {
        check()?;
        let turning_y = if force == 0.0 {
            f64::NAN
        } else {
            height + traction / force
        };
        let (q, surface, wall, dissipation, mut peak, viscosity, shear) = match &self.source {
            Source::Newtonian(material) => {
                let film = ShearedFilm::new(*material, height, force / self.density, traction)?;
                let viscosity = finite(material.density_kg_m3 * material.kinematic_viscosity_m2_s)?;
                let mut peak = film.surface_velocity_m_s().abs();
                if turning_y > 0.0 && turning_y < height {
                    peak = peak.max(film.velocity_m_s(turning_y)?.abs());
                }
                (
                    film.flow_per_width_m2_s(),
                    film.surface_velocity_m_s(),
                    film.wall_shear_pa(),
                    film.power().viscous_dissipation_w_m2,
                    peak,
                    viscosity,
                    film.wall_shear_pa().abs().max(traction.abs()) / viscosity,
                )
            }
            Source::Carreau(sources) => {
                let slit = sources[from].with_gap(height)?;
                let interrupted = Cell::new(false);
                let native_check = || {
                    if check().is_err() {
                        interrupted.set(true);
                        Err(MaterialError::State {
                            what: "thin-film evolution cancelled".into(),
                        })
                    } else {
                        Ok(())
                    }
                };
                let response = (|| -> Result<_, MaterialError> {
                    let profile = slit.profile_with_upper_traction_check(
                        force,
                        0.0,
                        traction,
                        controls,
                        native_check,
                    )?;
                    let flow = profile.flow();
                    let wall = profile.sample_with_check(0.0, native_check)?;
                    let surface = profile.sample_with_check(height, native_check)?;
                    let mut peak = surface.axial_velocity_m_s.abs();
                    if turning_y > 0.0 && turning_y < height {
                        peak = peak.max(
                            profile
                                .sample_with_check(turning_y, native_check)?
                                .axial_velocity_m_s
                                .abs(),
                        );
                    }
                    Ok((
                        flow.volumetric_flow_m3_s,
                        surface.axial_velocity_m_s,
                        wall.shear_stress_pa,
                        flow.dissipated_power_w,
                        peak,
                        slit.zero_shear_viscosity_pa_s()
                            .min(wall.dynamic_viscosity_pa_s)
                            .min(surface.dynamic_viscosity_pa_s),
                        wall.velocity_gradient_per_s
                            .abs()
                            .max(surface.velocity_gradient_per_s.abs()),
                    ))
                })();
                if interrupted.get() {
                    return Err(ThinFilmError::Cancelled);
                }
                response?
            }
        };
        peak = peak.max((q / height).abs());
        if ![q, surface, wall, dissipation, peak, shear]
            .iter()
            .all(|x| x.is_finite())
            || dissipation < 0.0
            || !positive(viscosity)
        {
            return Err(ThinFilmError::Invalid(
                "unrepresentable original film profile",
            ));
        }
        check()?;
        Ok(ThinFilmFace {
            from,
            to,
            thickness_m: height,
            flow_per_width_m2_s: q,
            surface_velocity_m_s: surface,
            wall_shear_pa: wall,
            force_density_pa_per_m: force,
            surface_shear_pa: traction,
            viscous_dissipation_w_m2: dissipation,
            peak_speed_m_s: peak,
            minimum_viscosity_pa_s: viscosity,
            maximum_shear_rate_per_s: shear,
        })
    }
}

fn faces(boundary: ThinFilmBoundary, cells: usize) -> usize {
    match boundary {
        ThinFilmBoundary::Periodic => cells,
        ThinFilmBoundary::Closed => cells - 1,
    }
}
fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
fn nonnegative(x: f64) -> bool {
    x.is_finite() && x >= 0.0
}
fn finite(x: f64) -> Result<f64, ThinFilmError> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(ThinFilmError::Invalid("unrepresentable film arithmetic"))
    }
}
fn total(values: impl IntoIterator<Item = f64>) -> Result<f64, ThinFilmError> {
    let mut sum = Dd::from_f64(0.0);
    for value in values {
        sum = sum + Dd::from_f64(finite(value)?);
    }
    finite(sum.to_f64())
}
fn admit(limits: ThinFilmLimits, controls: ThinFilmControls) -> Result<(), ThinFilmError> {
    controls.solver.validate().map_err(ThinFilmError::Solver)?;
    if ![
        limits.minimum_thickness_m,
        limits.maximum_thickness_m,
        limits.characteristic_length_m,
        limits.maximum_height_to_length,
        limits.maximum_slope,
        limits.maximum_reduced_reynolds,
        controls.height_scale_m,
        controls.finite_difference_step,
    ]
    .into_iter()
    .all(positive)
        || limits.minimum_thickness_m >= limits.maximum_thickness_m
        || limits.maximum_height_to_length >= 1.0
        || limits.maximum_slope >= 1.0
        || limits.maximum_reduced_reynolds >= 1.0
        || !(1..=1024).contains(&controls.maximum_newton_iterations)
        || !(1..=1_000_000).contains(&controls.maximum_endpoint_evaluations)
        || controls.maximum_profile_evaluations == 0
        || ![
            controls.cell_tolerance_m,
            controls.volume_tolerance_m2,
            controls.energy_tolerance_j_per_m,
        ]
        .into_iter()
        .all(nonnegative)
    {
        return Err(ThinFilmError::Invalid(
            "positive physical scales and bounded film solve controls required",
        ));
    }
    Ok(())
}
fn difference_error(error: FiniteDifferenceError<ThinFilmError>) -> ThinFilmError {
    match error {
        FiniteDifferenceError::Source(error) => error,
        FiniteDifferenceError::Invalid(message) => ThinFilmError::Invalid(message),
    }
}

#[cfg(test)]
mod tests;
