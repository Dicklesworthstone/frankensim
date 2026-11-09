//! Transient planar films with the original thermal Carreau constitutive law.
//!
//! Each face uses tau = mu(T, |du/dy|) du/dy and its positive stress tangent.
//! Implicit midpoint resolves momentum diffusion through startup and reversal;
//! it neither freezes an apparent viscosity nor adds steady flow profiles.
//! Newton steps use the original tridiagonal solver and a bounded line search.
//! The Newtonian owner supplies the same geometry, reusable storage and exact
//! discrete work account. Temperature may vary across the original shear faces;
//! the prescribed field, density and thickness are held during each interval.
//! An explicitly prescribed Navier length b gives u_wall_fluid-U_wall=b*gamma.
//! The ORIGINAL wall viscosity still depends on temperature and shear: a
//! constant b is NOT a constant interfacial friction coefficient mu/b.
//! This is a locally parallel generalized-Newtonian model: no elastic memory,
//! free-surface evolution, heat deposition or turbulent-transition prediction.

use fs_material::MaterialError;
use fs_material::fluid_transport::rheology::ThermalCarreauYasuda;
use fs_material::fluid_transport::rheology::pipe::{CarreauPipeControls, PreparedCarreau};

/// Native caloric evolution with this same resolved momentum history.
#[cfg(feature = "implicit-enthalpy-coupling")]
pub mod thermal;

use super::{
    FilmDrive, FilmMaterial, FilmStep, FilmStepLimits, MAX_FILM_CELLS, TransientFilm,
    TransientFilmError, filled, finite, poll, positive, stage_cell_heat,
};

/// Fixed SI geometry/inertia and a physical source-supported shear ceiling.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TransientCarreauFilmSpec {
    pub density_kg_m3: f64,
    pub thickness_m: f64,
    /// Includes wall, interior and free-surface rates at initial, midpoint and
    /// final profiles. This authored envelope is not a transition certificate.
    pub maximum_laminar_shear_rate_per_s: f64,
}

/// Finite nonlinear and constitutive work bounds, independent of physical speed
/// and absolute momentum-residual limits supplied through [`FilmStepLimits`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TransientCarreauControls {
    pub maximum_newton_iterations: u32,
    pub maximum_backtracks: u32,
    pub maximum_material_evaluations: u64,
    /// Relative stress residual AND Newton rate correction in the canonical
    /// free-surface stress inverse, strictly between zero and one.
    pub relative_surface_stress_tolerance: f64,
    pub maximum_surface_iterations: u32,
}

impl Default for TransientCarreauControls {
    fn default() -> Self {
        Self {
            maximum_newton_iterations: 64,
            maximum_backtracks: 24,
            maximum_material_evaluations: 1_000_000,
            relative_surface_stress_tolerance: 1e-12,
            maximum_surface_iterations: 96,
        }
    }
}

/// Refusal changes only reusable scratch, never the accepted profile or clock.
#[derive(Debug, Clone, PartialEq)]
pub enum TransientCarreauError {
    Flow(TransientFilmError),
    Material(MaterialError),
    ShearRateLimit { rate_per_s: f64, maximum_per_s: f64 },
    MaterialBudget,
    NewtonBudget,
    BacktrackBudget,
}
impl core::fmt::Display for TransientCarreauError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "transient Carreau film: {self:?}")
    }
}
impl std::error::Error for TransientCarreauError {}
impl From<TransientFilmError> for TransientCarreauError {
    fn from(value: TransientFilmError) -> Self {
        Self::Flow(value)
    }
}
impl From<MaterialError> for TransientCarreauError {
    fn from(value: MaterialError) -> Self {
        Self::Material(value)
    }
}
impl From<fs_la::tridiagonal::TridiagonalError> for TransientCarreauError {
    fn from(value: fs_la::tridiagonal::TridiagonalError) -> Self {
        Self::Flow(value.into())
    }
}

/// The common native flow/work account plus nonlinear constitutive diagnostics.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TransientCarreauStep {
    pub flow: FilmStep,
    pub maximum_shear_rate_per_s: f64,
    pub newton_iterations: u32,
    pub material_evaluations: u64,
    /// mu(g_surface)*g_surface minus the imposed signed traction [Pa].
    /// The common work account uses the imposed traction, with this separate
    /// inverse defect retained instead of silently changing the boundary load.
    pub surface_stress_residual_pa: f64,
}

/// One retained velocity history and original immutable constitutive source.
/// All numerical arrays are prepared once; warm steps allocate no workspace.
#[derive(Debug)]
pub struct TransientCarreauFilm {
    // The wall-face preparation also retains the common original law and rate
    // ceiling. Empty face_materials selects the allocation-free uniform lane.
    material: PreparedCarreau,
    face_materials: Vec<PreparedCarreau>,
    candidate_materials: Vec<PreparedCarreau>,
    flow: TransientFilm,
    trial: Vec<f64>,
    stress: Vec<f64>,
    tangent: Vec<f64>,
}

fn admit_rate(rate: f64, maximum: f64) -> Result<(), TransientCarreauError> {
    if !rate.is_finite() || rate.abs() > maximum {
        Err(TransientCarreauError::ShearRateLimit {
            rate_per_s: rate.abs(),
            maximum_per_s: maximum,
        })
    } else {
        Ok(())
    }
}

fn query(
    material: &PreparedCarreau,
    rate: f64,
    evaluations: &mut u64,
    budget: u64,
) -> Result<[f64; 2], TransientCarreauError> {
    if *evaluations >= budget {
        return Err(TransientCarreauError::MaterialBudget);
    }
    *evaluations += 1;
    Ok(material.viscosity_and_tangent(rate.abs())?)
}

impl TransientCarreauFilm {
    /// Prepare an initially resting film using the original shifted material.
    pub fn new(
        law: ThermalCarreauYasuda,
        temperature_k: f64,
        spec: TransientCarreauFilmSpec,
        cells: usize,
    ) -> Result<Self, TransientCarreauError> {
        if !(2..=MAX_FILM_CELLS).contains(&cells) {
            return Err(TransientFilmError::Invalid("film cell count").into());
        }
        Self::from_profile(law, temperature_k, spec, filled(cells)?, 0.0)
    }

    /// Import actual cell-center velocities. The first drive supplies the wall
    /// velocity, so its wall-rate admission is deferred to the first step.
    pub fn from_profile(
        law: ThermalCarreauYasuda,
        temperature_k: f64,
        spec: TransientCarreauFilmSpec,
        velocities_m_s: Vec<f64>,
        time_s: f64,
    ) -> Result<Self, TransientCarreauError> {
        let material =
            PreparedCarreau::new(law, temperature_k, spec.maximum_laminar_shear_rate_per_s)?;
        let viscosity = positive(material.zero_shear_viscosity_pa_s() / spec.density_kg_m3)?;
        // Capillarity does not enter a fixed flat film's tangential momentum.
        // This reference viscosity supplies only scales/storage admission;
        // every evolving face uses the original nonlinear material below.
        let flow = TransientFilm::from_profile(
            FilmMaterial {
                density_kg_m3: spec.density_kg_m3,
                kinematic_viscosity_m2_s: viscosity,
                surface_tension_n_m: 0.0,
            },
            spec.thickness_m,
            velocities_m_s,
            time_s,
        )?;
        for pair in flow.velocity.windows(2) {
            admit_rate(
                (pair[1] - pair[0]) / flow.spacing,
                material.maximum_shear_rate_per_s(),
            )?;
        }
        let n = flow.velocity.len();
        Ok(Self {
            material,
            face_materials: Vec::new(),
            candidate_materials: Vec::new(),
            flow,
            trial: filled(n)?,
            stress: filled(n + 1)?,
            tangent: filled(n + 1)?,
        })
    }

    /// Original wall-face material. For a spatial temperature field, use
    /// [`Self::material_at_face`] to query each face's actual preparation.
    pub const fn material(&self) -> &PreparedCarreau {
        &self.material
    }

    /// Original source at shear face 0..=N, located at y = face * spacing.
    /// Face zero is the wall and face N is the free surface. No cell-averaged
    /// apparent viscosity replaces these face-local constitutive evaluations.
    pub fn material_at_face(&self, face: usize) -> Option<&PreparedCarreau> {
        if face > self.flow.velocity.len() {
            None
        } else if self.face_materials.is_empty() {
            Some(&self.material)
        } else {
            self.face_materials.get(face)
        }
    }

    /// Whether one temperature covers every face. Spatial updates whose values
    /// are all equal retain the same compact uniform representation.
    pub fn is_isothermal(&self) -> bool {
        self.face_materials.is_empty()
    }
    pub const fn density_kg_m3(&self) -> f64 {
        self.flow.material.density_kg_m3
    }
    pub const fn thickness_m(&self) -> f64 {
        self.flow.thickness
    }
    pub const fn spacing_m(&self) -> f64 {
        self.flow.spacing
    }
    /// Prescribed finite Navier length; zero is the original no-slip boundary.
    pub const fn wall_slip_length_m(&self) -> f64 {
        self.flow.wall_slip_length_m()
    }
    /// Hold u_fluid_wall-U_wall = b*gamma_wall with the original nonlinear
    /// constitutive stress tau=mu(T,|gamma|)*gamma. This is a prescribed length,
    /// not a fitted/material-sourced constant friction coefficient. Changes
    /// neither accepted momentum nor temperature preparations; refusal is atomic.
    pub fn set_wall_slip_length_m(&mut self, length_m: f64) -> Result<(), TransientCarreauError> {
        self.flow.set_wall_slip_length_m(length_m)?;
        Ok(())
    }
    pub const fn accepted_time_s(&self) -> f64 {
        self.flow.time
    }
    pub fn velocities_m_s(&self) -> &[f64] {
        &self.flow.velocity
    }
    /// Shortest cell-scale time over the faces' zero-shear viscosities. Actual
    /// nonlinear diffusion and forcing still require time/space refinement.
    pub const fn zero_shear_diffusion_time_s(&self) -> f64 {
        self.flow.diffusion_time
    }

    /// Reprepare the SAME material at a prescribed temperature without changing
    /// velocity, physical time, density, thickness or the physical rate ceiling.
    /// This is not thermal evolution or heat deposition. Refusal is atomic.
    pub fn set_temperature(
        &mut self,
        temperature_k: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), TransientCarreauError> {
        poll(&mut cancelled)?;
        let material = PreparedCarreau::new(
            *self.material.law(),
            temperature_k,
            self.material.maximum_shear_rate_per_s(),
        )?;
        let nu = positive(material.zero_shear_viscosity_pa_s() / self.density_kg_m3())?;
        let conductance = positive(material.zero_shear_viscosity_pa_s() / self.flow.spacing)?;
        positive(2.0 * conductance)?;
        let diffusion_time = positive(self.flow.mass_per_cell_area / conductance)?;
        poll(&mut cancelled)?;
        self.material = material;
        self.face_materials.clear();
        self.flow.material.kinematic_viscosity_m2_s = nu;
        self.flow.shear_conductance = conductance;
        self.flow.diffusion_time = diffusion_time;
        self.flow.factored_duration = None;
        Ok(())
    }

    /// Bind N+1 prescribed face temperatures to the SAME original material.
    ///
    /// The caller supplies actual wall/surface temperatures and interior face
    /// traces from its thermal owner; no conductivity, interpolation rule,
    /// caloric inventory or boundary temperature is inferred. The entire field
    /// is admitted before publication. Velocity, clock, density, thickness and
    /// the physical shear-rate ceiling remain unchanged. Spatial preparation
    /// and its staging buffer are allocated on first use and reused thereafter;
    /// ordinary momentum steps do not re-evaluate the temperature shift.
    pub fn set_face_temperatures(
        &mut self,
        temperatures_k: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), TransientCarreauError> {
        poll(&mut cancelled)?;
        let faces = self.flow.velocity.len() + 1;
        if temperatures_k.len() != faces {
            return Err(TransientFilmError::Invalid("film needs N+1 face temperatures").into());
        }
        let mut uniform = true;
        let mut unchanged = true;
        for (i, &temperature) in temperatures_k.iter().enumerate() {
            if i % 128 == 0 {
                poll(&mut cancelled)?;
            }
            uniform &= temperature == temperatures_k[0];
            unchanged &= temperature == self.material_at_face(i).unwrap().temperature_k();
        }
        if unchanged {
            poll(&mut cancelled)?;
            return Ok(());
        }
        if uniform {
            return self.set_temperature(temperatures_k[0], cancelled);
        }
        self.candidate_materials.clear();
        self.candidate_materials
            .try_reserve_exact(faces)
            .map_err(|_| TransientFilmError::Capacity)?;
        // Reserve both sides now, so a later successful rebind can reuse its
        // swapped staging buffer without allocating.
        if self.face_materials.capacity() < faces {
            self.face_materials
                .try_reserve_exact(faces - self.face_materials.len())
                .map_err(|_| TransientFilmError::Capacity)?;
        }
        let mut maximum_viscosity = 0.0_f64;
        for (i, &temperature) in temperatures_k.iter().enumerate() {
            if i % 128 == 0 {
                poll(&mut cancelled)?;
            }
            let material = PreparedCarreau::new(
                *self.material.law(),
                temperature,
                self.material.maximum_shear_rate_per_s(),
            )?;
            maximum_viscosity = maximum_viscosity.max(material.zero_shear_viscosity_pa_s());
            self.candidate_materials.push(material);
        }
        let nu = positive(maximum_viscosity / self.density_kg_m3())?;
        let conductance = positive(maximum_viscosity / self.flow.spacing)?;
        positive(2.0 * conductance)?;
        let diffusion_time = positive(self.flow.mass_per_cell_area / conductance)?;
        poll(&mut cancelled)?;
        core::mem::swap(&mut self.face_materials, &mut self.candidate_materials);
        self.material = self.face_materials[0];
        self.flow.material.kinematic_viscosity_m2_s = nu;
        self.flow.shear_conductance = conductance;
        self.flow.diffusion_time = diffusion_time;
        self.flow.factored_duration = None;
        Ok(())
    }

    fn admit_profile(
        &self,
        profile: &[f64],
        drive: FilmDrive,
        limits: FilmStepLimits,
        surface_offset: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, TransientCarreauError> {
        let mut maximum = 0.0_f64;
        for (i, &velocity) in profile.iter().enumerate() {
            if i % 128 == 0 {
                poll(cancelled)?;
            }
            if !velocity.is_finite() || velocity.abs() > limits.maximum_speed_m_s {
                return Err(TransientFilmError::SpeedLimit.into());
            }
            let rate = if i == 0 {
                (velocity - drive.wall_velocity_m_s) / self.flow.wall_distance()
            } else {
                (velocity - profile[i - 1]) / self.flow.spacing
            };
            admit_rate(rate, self.material.maximum_shear_rate_per_s())?;
            maximum = maximum.max(rate.abs());
        }
        let surface = finite(profile[profile.len() - 1] + surface_offset)?;
        if surface.abs() > limits.maximum_speed_m_s {
            return Err(TransientFilmError::SpeedLimit.into());
        }
        Ok(maximum)
    }

    // Evaluate an endpoint's exact midpoint residual and positive tridiagonal
    // Jacobian together. Accepted line-search trials retain these same arrays,
    // avoiding a second material traversal at the next Newton iteration.
    fn evaluate(
        &mut self,
        trial: bool,
        duration_s: f64,
        drive: FilmDrive,
        evaluations: &mut u64,
        budget: u64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, TransientCarreauError> {
        let flow = &mut self.flow;
        let profile = if trial { &self.trial } else { &flow.candidate };
        let n = profile.len();
        let weight = positive(duration_s / flow.mass_per_cell_area)?;
        let body_increment = finite(duration_s * drive.acceleration_m_s2)?;
        for (i, &velocity) in profile.iter().enumerate() {
            if i % 128 == 0 {
                poll(cancelled)?;
            }
            let delta = finite(velocity - flow.velocity[i])?;
            flow.midpoint[i] = finite(0.5f64.mul_add(delta, flow.velocity[i]))?;
            let distance = if i == 0 {
                flow.wall_distance()
            } else {
                flow.spacing
            };
            let previous = if i == 0 {
                drive.wall_velocity_m_s
            } else {
                flow.midpoint[i - 1]
            };
            let rate = finite((flow.midpoint[i] - previous) / distance)?;
            admit_rate(rate, self.material.maximum_shear_rate_per_s())?;
            let material = if self.face_materials.is_empty() {
                &self.material
            } else {
                &self.face_materials[i]
            };
            let [mu, tangent] = query(material, rate, evaluations, budget)?;
            self.stress[i] = finite(mu * rate)?;
            if rate != 0.0 && self.stress[i] == 0.0 {
                return Err(TransientFilmError::NonRepresentable.into());
            }
            self.tangent[i] = positive(0.5 * weight * (tangent / distance))?;
        }
        self.stress[n] = drive.surface_shear_pa;
        self.tangent[n] = 0.0;
        let mut maximum = 0.0_f64;
        for i in 0..n {
            if i % 128 == 0 {
                poll(cancelled)?;
            }
            flow.rhs[i] = finite(
                (profile[i] - flow.velocity[i])
                    - weight.mul_add(self.stress[i + 1] - self.stress[i], body_increment),
            )?;
            maximum = maximum.max(flow.rhs[i].abs());
            flow.diagonal[i] = finite(1.0 + self.tangent[i] + self.tangent[i + 1])?;
            if i + 1 < n {
                flow.adjacent[i] = -self.tangent[i + 1];
            }
        }
        Ok(maximum)
    }

    /// Advance one held-temperature/drive interval. Source preparation and all
    /// Newton work retain the original local stress law. Every numerical sweep
    /// polls at least every 128 cells; surface inversion is bounded by 128
    /// canonical iterations. Failed trials cannot publish a partial endpoint.
    pub fn step(
        &mut self,
        duration_s: f64,
        drive: FilmDrive,
        limits: FilmStepLimits,
        controls: TransientCarreauControls,
        cancelled: impl FnMut() -> bool,
    ) -> Result<TransientCarreauStep, TransientCarreauError> {
        self.step_inner(duration_s, drive, limits, controls, None, cancelled)
    }

    /// Advance and export the SAME interval's resolved viscous heat into N
    /// caller-owned cells, in J per square metre of wetted wall area.
    ///
    /// Each interior shear segment contributes half its nonnegative work to
    /// each adjacent cell; boundary half-segments heat the adjacent cell only.
    /// Their sum is the flow report's viscous heat up to accumulation roundoff.
    /// The separate wall-slip loss is not deposited in these fluid cells; its
    /// receiver must be supplied explicitly by a thermal/interface owner.
    /// Multiply by the explicitly represented wall area before depositing into
    /// a thermal owner's joule inventory. Never deposit body/wall/traction work
    /// again: kinetic storage has already been removed from these viscous losses.
    /// Heat, velocity and time publish together, after the final cancellation
    /// check. Any refusal leaves the caller's output and accepted history intact.
    pub fn step_with_heat(
        &mut self,
        duration_s: f64,
        drive: FilmDrive,
        limits: FilmStepLimits,
        controls: TransientCarreauControls,
        cell_heat_j_m2: &mut [f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<TransientCarreauStep, TransientCarreauError> {
        self.step_inner(
            duration_s,
            drive,
            limits,
            controls,
            Some(cell_heat_j_m2),
            cancelled,
        )
    }

    fn step_inner(
        &mut self,
        duration_s: f64,
        drive: FilmDrive,
        limits: FilmStepLimits,
        controls: TransientCarreauControls,
        cell_heat_j_m2: Option<&mut [f64]>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<TransientCarreauStep, TransientCarreauError> {
        poll(&mut cancelled)?;
        if cell_heat_j_m2
            .as_ref()
            .is_some_and(|heat| heat.len() != self.flow.velocity.len())
        {
            return Err(TransientFilmError::Invalid("film needs N heat output cells").into());
        }
        if !duration_s.is_finite()
            || duration_s <= 0.0
            || !limits.maximum_speed_m_s.is_finite()
            || limits.maximum_speed_m_s <= 0.0
            || !limits.maximum_velocity_residual_m_s.is_finite()
            || limits.maximum_velocity_residual_m_s <= 0.0
            || ![
                drive.acceleration_m_s2,
                drive.wall_velocity_m_s,
                drive.surface_shear_pa,
            ]
            .into_iter()
            .all(f64::is_finite)
            || !(1..=128).contains(&controls.maximum_newton_iterations)
            || controls.maximum_backtracks > 64
            || controls.maximum_material_evaluations == 0
            || !controls.relative_surface_stress_tolerance.is_finite()
            || controls.relative_surface_stress_tolerance <= 0.0
            || controls.relative_surface_stress_tolerance >= 1.0
            || !(1..=128).contains(&controls.maximum_surface_iterations)
        {
            return Err(TransientFilmError::Invalid("Carreau film step inputs").into());
        }
        let end_time = finite(self.flow.time + duration_s)?;
        if end_time <= self.flow.time {
            return Err(TransientFilmError::NonRepresentable.into());
        }
        if drive.wall_velocity_m_s.abs() > limits.maximum_speed_m_s {
            return Err(TransientFilmError::SpeedLimit.into());
        }
        let surface_material = *self.material_at_face(self.flow.velocity.len()).unwrap();
        let (surface_rate, inverse_iterations) = surface_material.signed_shear_rate_at_stress(
            drive.surface_shear_pa,
            CarreauPipeControls {
                relative_constitutive_tolerance: controls.relative_surface_stress_tolerance,
                maximum_constitutive_iterations: controls.maximum_surface_iterations.min(
                    controls
                        .maximum_material_evaluations
                        .min(u64::from(u32::MAX)) as u32,
                ),
                ..CarreauPipeControls::default()
            },
        )?;
        let mut evaluations = u64::from(inverse_iterations);
        poll(&mut cancelled)?;
        let [surface_mu, _] = query(
            &surface_material,
            surface_rate,
            &mut evaluations,
            controls.maximum_material_evaluations,
        )?;
        let surface_residual = finite(surface_mu * surface_rate - drive.surface_shear_pa)?;
        let surface_offset = finite(0.5 * self.flow.spacing * surface_rate)?;
        let initial_rate = self.admit_profile(
            &self.flow.velocity,
            drive,
            limits,
            surface_offset,
            &mut cancelled,
        )?;
        let n = self.flow.velocity.len();
        for i in 0..n {
            if i % 128 == 0 {
                poll(&mut cancelled)?;
            }
            self.flow.candidate[i] = self.flow.velocity[i];
        }
        let mut residual = self.evaluate(
            false,
            duration_s,
            drive,
            &mut evaluations,
            controls.maximum_material_evaluations,
            &mut cancelled,
        )?;
        let mut iterations = 0;
        while residual > limits.maximum_velocity_residual_m_s {
            if iterations == controls.maximum_newton_iterations {
                return Err(TransientCarreauError::NewtonBudget);
            }
            let flow = &mut self.flow;
            flow.linear.factor(
                &flow.adjacent,
                &flow.diagonal,
                &flow.adjacent,
                &mut cancelled,
            )?;
            for (i, value) in flow.rhs.iter_mut().enumerate() {
                if i % 128 == 0 {
                    poll(&mut cancelled)?;
                }
                *value = -*value;
            }
            flow.linear
                .solve_factored_into(&flow.rhs, &mut flow.delta, &mut cancelled)?;
            let mut accepted = false;
            let mut scale = 1.0_f64;
            for _ in 0..=controls.maximum_backtracks {
                for i in 0..n {
                    if i % 128 == 0 {
                        poll(&mut cancelled)?;
                    }
                    self.trial[i] = scale.mul_add(self.flow.delta[i], self.flow.candidate[i]);
                }
                match self.admit_profile(&self.trial, drive, limits, surface_offset, &mut cancelled)
                {
                    Ok(_) => {
                        let next = self.evaluate(
                            true,
                            duration_s,
                            drive,
                            &mut evaluations,
                            controls.maximum_material_evaluations,
                            &mut cancelled,
                        )?;
                        if next < residual || next <= limits.maximum_velocity_residual_m_s {
                            core::mem::swap(&mut self.flow.candidate, &mut self.trial);
                            residual = next;
                            accepted = true;
                            break;
                        }
                    }
                    Err(
                        TransientCarreauError::ShearRateLimit { .. }
                        | TransientCarreauError::Flow(TransientFilmError::SpeedLimit),
                    ) => {}
                    Err(error) => return Err(error),
                }
                scale *= 0.5;
            }
            if !accepted {
                return Err(TransientCarreauError::BacktrackBudget);
            }
            iterations += 1;
        }
        let final_rate = self.admit_profile(
            &self.flow.candidate,
            drive,
            limits,
            surface_offset,
            &mut cancelled,
        )?;
        for i in 0..n {
            if i % 128 == 0 {
                poll(&mut cancelled)?;
            }
            self.flow.delta[i] = finite(self.flow.candidate[i] - self.flow.velocity[i])?;
        }
        let report = self.flow.report_candidate(
            duration_s,
            drive,
            limits,
            surface_offset,
            |face| self.stress[face],
            &mut cancelled,
        )?;
        let step = TransientCarreauStep {
            flow: report,
            maximum_shear_rate_per_s: initial_rate.max(final_rate).max(surface_rate.abs()),
            newton_iterations: iterations,
            material_evaluations: evaluations,
            surface_stress_residual_pa: surface_residual,
        };
        if cell_heat_j_m2.is_some() {
            let wall_bulk_difference = self
                .flow
                .wall_differences(self.flow.midpoint[0] - drive.wall_velocity_m_s)?[0];
            stage_cell_heat(
                &self.flow.midpoint,
                duration_s,
                drive,
                wall_bulk_difference,
                surface_offset,
                |face| self.stress[face],
                &mut self.flow.heat,
                &mut cancelled,
            )?;
        }
        poll(&mut cancelled)?;
        if let Some(heat) = cell_heat_j_m2 {
            heat.copy_from_slice(&self.flow.heat);
        }
        core::mem::swap(&mut self.flow.velocity, &mut self.flow.candidate);
        self.flow.time = end_time;
        Ok(step)
    }
}

#[cfg(test)]
mod tests;
