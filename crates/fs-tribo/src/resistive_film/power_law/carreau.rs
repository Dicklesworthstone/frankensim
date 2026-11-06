//! Original thermal Carreau passages in the live Reynolds pressure graph.
//!
//! Squeeze, pressure supplies and both tangential walls participate in one
//! continuity solve. The original slit owner evaluates the combined nonlinear
//! pressure/sliding profile. This module shares the existing bounded pressure
//! iteration, factorization and work-conjugate directional tangent; it neither
//! fits a power law nor substitutes a frozen apparent viscosity.

use fs_material::MaterialError;
use fs_material::fluid_transport::rheology::slit::{CarreauSlitControls, LaminarCarreauSlit};

pub use super::super::shear::WallVelocityPorts;
pub use super::PowerLawControls as ReynoldsControls;
use super::nonlinear::ReynoldsLaw;
use super::shear::{Coefficients, Motion, WallResponse, compile_walls};
use super::{
    Channel, ChannelWork, FilmError, FilmReport, MAX_CHANNELS, ReservoirPressures, ResistiveFilm,
    dot, fail,
};

/// Independent bounded pressure-island and original constitutive controls.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CarreauFilmControls {
    pub pressure: ReynoldsControls,
    pub passage: CarreauSlitControls,
}

impl Default for CarreauFilmControls {
    fn default() -> Self {
        Self {
            pressure: ReynoldsControls {
                maximum_iterations: 64,
                maximum_line_search_steps: 48,
                absolute_flow_tolerance_m3_s: 1e-18,
                relative_flow_tolerance: 1e-8,
                absolute_power_tolerance_w: 1e-12,
                relative_power_tolerance: 1e-8,
            },
            passage: CarreauSlitControls::default(),
        }
    }
}

/// Actual passage flow, resisting traction and physical viscous power.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct CarreauFilmChannelReport {
    pub volume_flow_m3_s: f64,
    pub pressure_flow_m3_s: f64,
    pub entrained_flow_m3_s: f64,
    /// Lower/upper resisting tractions along the original channel [Pa].
    pub wall_tractions_pa: [f64; 2],
    /// Signed work at the requested wall velocities [W].
    pub wall_power_w: f64,
    /// Original accepted-profile dissipation, not a numerical work remainder.
    pub dissipated_power_w: f64,
    /// Requested pressure/sliding work minus accepted-profile heat [W].
    pub constitutive_work_defect_w: f64,
}

#[derive(Debug, Clone, Copy)]
pub struct CarreauFilmReport {
    pub film: FilmReport,
    pub mechanical_power_w: f64,
    pub reservoir_power_w: f64,
    pub maximum_flow_defect_m3_s: f64,
    /// Total mechanical plus reservoir work minus physical heat [W]. Includes
    /// numerical pressure-continuity and local constitutive work residuals.
    pub work_defect_w: f64,
    /// Sum of the original profiles' local constitutive work residuals [W].
    pub constitutive_work_defect_w: f64,
    pub iterations: usize,
    pub backtracks: usize,
}

/// Incompressible full-film Carreau resistance with explicit pressure supplies.
///
/// Every passage retains its original law, temperature and supported rate
/// ceiling. The supplied slit width/length must match the original graph;
/// its prepared aperture is rebound to the graph's actual affine gap. Thus a
/// reference gap need not be the current gap. A closed fixed-wall passage has
/// exactly zero flow; any authored moving-wall passage must remain open.
///
/// Configuration, source states and supplies are frozen during evaluation.
/// Successful evaluation and tangents use bounded stack scratch and no per-query heap
/// allocation. All outputs publish together after source, continuity and work
/// admission. This is the thin-gap, quasi-steady Reynolds approximation: no
/// cavitation, Barus pressure dependence, fluid inertia, thermal clock, yield
/// memory or independent slip model is inferred.
#[derive(Debug, Clone)]
pub struct CarreauFilm {
    film: ResistiveFilm,
    passages: Vec<LaminarCarreauSlit>,
    reservoirs: ReservoirPressures,
    walls: Vec<WallVelocityPorts>,
    coefficients: Vec<Coefficients>,
    constant: Vec<bool>,
    all_newtonian: bool,
    controls: CarreauFilmControls,
}

impl CarreauFilm {
    /// Bind original-order material passages, reservoir pressures and wall rows.
    /// Each exterior channel needs Some(pressure), and each internal one None.
    pub fn new(
        geometry: &ResistiveFilm,
        passages: Vec<LaminarCarreauSlit>,
        reservoir_pressures_pa: &[Option<f64>],
        walls: Vec<WallVelocityPorts>,
        controls: CarreauFilmControls,
    ) -> Result<Self, FilmError> {
        controls.pressure.validate()?;
        if geometry.pressure_viscosity_coefficient_per_pa() != 0.0 {
            return fail("Carreau film cannot discard an existing Barus pressure law");
        }
        if passages.len() != geometry.channel_count() || walls.len() != geometry.channel_count() {
            return fail(
                "Carreau film needs one original material passage and wall row per channel",
            );
        }
        let mut viscosities = Vec::with_capacity(passages.len());
        let mut constant = Vec::with_capacity(passages.len());
        for (passage, edge) in passages.iter().zip(geometry.channels.iter()) {
            if passage.width_m() != edge.width_m || passage.length_m() != edge.length_m {
                return fail("Carreau passage width/length differs from the original film graph");
            }
            // Admit the original constitutive controls, without inventing a
            // material envelope or rebuilding its frozen thermal preparation.
            passage
                .profile_at_pressure_and_walls(0.0, [0.0; 2], controls.passage)
                .map_err(|_| {
                    FilmError("Carreau passage refused its numerical controls or rest state")
                })?;
            viscosities.push(passage.zero_shear_viscosity_pa_s());
            let [zero, infinity, time, index, _] = passage.law().curve().parameters();
            constant.push(time == 0.0 || index == 1.0 || zero == infinity);
        }
        let film = geometry.with_channel_viscosities(&viscosities)?;
        let reservoirs = ReservoirPressures::new(&film, reservoir_pressures_pa)?;
        let (coefficients, _) = compile_walls(&walls, geometry.port_count())?;
        let all_newtonian = constant.iter().all(|value| *value);
        Ok(Self {
            film,
            passages,
            reservoirs,
            walls,
            coefficients,
            constant,
            all_newtonian,
            controls,
        })
    }

    /// Original material preparations; actual evaluation gaps come from q.
    pub fn passages(&self) -> &[LaminarCarreauSlit] {
        &self.passages
    }
    pub const fn controls(&self) -> CarreauFilmControls {
        self.controls
    }
    pub fn wall_ports(&self) -> &[WallVelocityPorts] {
        &self.walls
    }
    pub fn reservoir_pressures_pa(&self) -> &[Option<f64>] {
        &self.reservoirs.reservoirs_pa
    }
    /// Admitted geometry with original source zero-shear coefficients. This
    /// view alone does not represent the nonlinear Carreau mechanics.
    pub const fn pressure_operator(&self) -> &ResistiveFilm {
        &self.film
    }

    // The conservative cavitation owner needs the SAME combined pressure and
    // sliding profile, independently of mean-wall liquid-content transport.
    // Even constant original laws pass through their source-rate admission.
    pub(in crate::resistive_film) fn passage_pressure_response(
        &self,
        index: usize,
        gap: f64,
        core_conductance: f64,
        pressure_drop: f64,
        relative_wall_speed: f64,
        check: &mut impl FnMut() -> Result<(), FilmError>,
    ) -> Result<(f64, f64), FilmError> {
        if gap <= 0.0 && relative_wall_speed != 0.0 {
            return fail("moving Carreau passage needs a positive gap");
        }
        let channel = self.channel_response(
            index,
            gap,
            core_conductance,
            pressure_drop,
            Motion {
                mean: 0.0,
                relative: relative_wall_speed,
                entrainment: 0.0,
            },
            relative_wall_speed != 0.0,
            check,
        )?;
        Ok((channel.flow, channel.differential))
    }

    /// Refresh actual source temperatures/laws without inventing a thermal
    /// state or elapsed time. Geometry, wall mappings and supplies are retained.
    pub fn with_passages(&self, passages: Vec<LaminarCarreauSlit>) -> Result<Self, FilmError> {
        Self::new(
            &self.film,
            passages,
            self.reservoir_pressures_pa(),
            self.walls.clone(),
            self.controls,
        )
    }

    pub fn with_reservoir_pressures(
        &self,
        pressures_pa: &[Option<f64>],
    ) -> Result<Self, FilmError> {
        Self::new(
            &self.film,
            self.passages.clone(),
            pressures_pa,
            self.walls.clone(),
            self.controls,
        )
    }

    pub fn evaluate_into(
        &self,
        q: &[f64],
        velocity: &[f64],
        force: &mut [f64],
        pressure: &mut [f64],
        channels: &mut [CarreauFilmChannelReport],
    ) -> Result<CarreauFilmReport, FilmError> {
        self.evaluate_with_check_into(q, velocity, force, pressure, channels, || Ok(()))
    }

    /// Same source-owned solve with bounded checks inside each original slit
    /// profile, pressure correction and final publication. Refusal changes none
    /// of the caller's output slices or the immutable source owner.
    pub fn evaluate_with_check_into(
        &self,
        q: &[f64],
        velocity: &[f64],
        force: &mut [f64],
        pressure: &mut [f64],
        channels: &mut [CarreauFilmChannelReport],
        mut check: impl FnMut() -> Result<(), FilmError>,
    ) -> Result<CarreauFilmReport, FilmError> {
        if force.len() != self.film.port_count()
            || pressure.len() != self.film.cell_count()
            || channels.len() != self.film.channel_count()
        {
            return fail("invalid Carreau film output shape");
        }
        let state = self.solve_with_motion(
            q,
            velocity,
            Some(&self.reservoirs),
            Some(&self.coefficients),
            &mut check,
        )?;
        let mut candidate = [CarreauFilmChannelReport::default(); MAX_CHANNELS];
        let mut constitutive_work_defect_w = 0.0;
        for (k, edge) in self.film.channels.iter().enumerate() {
            if k.is_multiple_of(32) {
                check()?;
            }
            let local = state.channels[k];
            let Some(work) = local.work else {
                continue;
            };
            let dp = state.system.potential[edge.from]
                - edge
                    .to
                    .map_or(self.reservoirs.potential[k], |j| state.system.potential[j]);
            let h = edge.gap.reference_m - dot(&edge.gap.closure, q);
            let motion = self.coefficients[k].motion(h, edge.width_m, velocity)?;
            let half_stress = dp * (0.5 * h / edge.length_m);
            let wall = local.wall;
            candidate[k] = CarreauFilmChannelReport {
                volume_flow_m3_s: local.flow,
                pressure_flow_m3_s: wall.pressure_flow,
                entrained_flow_m3_s: motion.entrainment,
                wall_tractions_pa: [
                    -wall.mid_stress - half_stress,
                    wall.mid_stress - half_stress,
                ],
                wall_power_w: wall.mean_effort * motion.mean
                    + wall.relative_effort * motion.relative,
                dissipated_power_w: work.dissipated_power_w,
                constitutive_work_defect_w: work.constitutive_work_defect_w,
            };
            constitutive_work_defect_w += work.constitutive_work_defect_w;
            if candidate[k]
                .wall_tractions_pa
                .iter()
                .any(|v| !v.is_finite())
            {
                return fail("Carreau wall traction is unrepresentable");
            }
        }
        if !constitutive_work_defect_w.is_finite() {
            return fail("Carreau work sum is unrepresentable");
        }
        let report = CarreauFilmReport {
            film: FilmReport {
                minimum_gap_m: state.system.minimum_gap,
                peak_pressure_pa: state.system.pressure[..self.film.cell_count()]
                    .iter()
                    .fold(0.0_f64, |m, p| m.max(p.abs())),
                dissipated_power_w: state.power,
            },
            mechanical_power_w: state.port_work,
            reservoir_power_w: state.reservoir_power,
            maximum_flow_defect_m3_s: state.maximum_residual,
            work_defect_w: state.port_work + state.reservoir_power - state.power,
            constitutive_work_defect_w,
            iterations: state.iterations,
            backtracks: state.backtracks,
        };
        check()?;
        force.copy_from_slice(&state.force[..self.film.port_count()]);
        pressure.copy_from_slice(&state.system.pressure[..self.film.cell_count()]);
        channels.copy_from_slice(&candidate[..self.film.channel_count()]);
        Ok(report)
    }

    /// Differentiate q, normal/tangential velocities and supplied pressures.
    /// The accepted profile's original positive mobility and cross response
    /// give the complete tangent through the same factored pressure Jacobian.
    #[allow(clippy::too_many_arguments)]
    pub fn tangent_into(
        &self,
        q: &[f64],
        velocity: &[f64],
        dq: &[f64],
        dv: &[f64],
        reservoir_directions_pa: &[f64],
        out: &mut [f64],
    ) -> Result<(), FilmError> {
        self.tangent_with_check_into(q, velocity, dq, dv, reservoir_directions_pa, out, || Ok(()))
    }

    #[allow(clippy::too_many_arguments)]
    pub fn tangent_with_check_into(
        &self,
        q: &[f64],
        velocity: &[f64],
        dq: &[f64],
        dv: &[f64],
        reservoir_directions_pa: &[f64],
        out: &mut [f64],
        mut check: impl FnMut() -> Result<(), FilmError>,
    ) -> Result<(), FilmError> {
        self.tangent_with_walls(
            q,
            velocity,
            dq,
            dv,
            reservoir_directions_pa,
            out,
            Some(&self.reservoirs),
            &self.coefficients,
            &mut check,
        )
    }
}

impl ReynoldsLaw for CarreauFilm {
    fn pressure_operator(&self) -> &ResistiveFilm {
        &self.film
    }
    fn pressure_controls(&self) -> ReynoldsControls {
        self.controls.pressure
    }
    fn native_newtonian(&self) -> bool {
        self.all_newtonian
    }
    fn core_requires_positive_conductance(&self, index: usize) -> bool {
        !self.constant[index]
    }

    fn channel_response(
        &self,
        index: usize,
        gap: f64,
        _core_conductance: f64,
        dp: f64,
        motion: Motion,
        _moving: bool,
        check: &mut impl FnMut() -> Result<(), FilmError>,
    ) -> Result<Channel, FilmError> {
        if gap <= 0.0 {
            return Ok(Channel::default());
        }
        let slit = self.passages[index]
            .with_gap(gap)
            .map_err(|_| FilmError("Carreau source refused the actual film aperture"))?;
        // Galilean translation is exact kinematics, not another constitutive
        // solve. Centering the walls preserves tiny sliding under huge mean
        // entrainment; the original source still resolves pressure and sliding
        // together, and supplies their actual accepted heat and cross tangent.
        let half_relative = 0.5 * motion.relative;
        if motion.relative != 0.0 && half_relative == 0.0 {
            return fail("Carreau relative wall velocity is unrepresentable");
        }
        let mut callback_failure = None;
        let profile = slit.profile_with_check(
            dp,
            [-half_relative, half_relative],
            self.controls.passage,
            || {
                check().map_err(|error| {
                    callback_failure = Some(error);
                    MaterialError::State {
                        what: "Carreau film caller check refused".into(),
                    }
                })
            },
        );
        if let Some(error) = callback_failure {
            return Err(error);
        }
        let profile = profile
            .map_err(|_| FilmError("Carreau passage source support or solve budget refused"))?;
        let flow = *profile.flow();
        let c = profile.stress_offset_pa();
        let pressure_flow = flow.pressure_flow_m3_s;
        let cross = flow.flow_relative_wall_velocity_derivative_m2;
        // Exact geometric derivative for a source independent of geometry:
        // Q_pressure = w*h^2*Phi(dp*h/L, relative/h).
        let aperture_derivative = (2.0 * pressure_flow + dp * flow.pressure_conductance_m3_s_pa
            - motion.relative * cross)
            / gap;
        let mean_effort = -slit.width_m() * gap * dp;
        let relative_effort = (slit.width_m() * slit.length_m()) * c;
        let total_flow = pressure_flow + motion.entrainment;
        if [
            aperture_derivative,
            mean_effort,
            relative_effort,
            total_flow,
        ]
        .iter()
        .any(|v| !v.is_finite())
        {
            return fail("Carreau film response or gap derivative is unrepresentable");
        }
        Ok(Channel {
            flow: total_flow,
            differential: flow.pressure_conductance_m3_s_pa,
            aperture_derivative,
            wall: WallResponse {
                mean_effort,
                relative_effort,
                relative_flow_derivative: cross,
                relative_effort_derivative: flow.relative_wall_resistance_n_s_m,
                pressure_flow,
                mid_stress: c,
                tangent_ambiguous: false,
            },
            work: Some(ChannelWork {
                dissipated_power_w: flow.dissipated_power_w,
                constitutive_work_defect_w: flow.work_residual_w,
            }),
        })
    }
}

#[cfg(test)]
mod tests;
