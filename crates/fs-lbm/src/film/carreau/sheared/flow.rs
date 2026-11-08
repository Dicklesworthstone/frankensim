//! Signed-flux thickness inversion on a caller-selected physical branch.
//!
//! With F=rho*g and z=h-y, q(h)=integral_0^h z*gamma(tau_s+F*z) dz.
//! Therefore dq/dh=h*gamma(tau_s+F*h): the only positive turning point
//! is a zero of WALL shear, independent of the particular monotone fluid law.
//! Opposed surface traction and body force can give two positive thicknesses
//! for the same flux. Neither a globally monotone inverse nor selection by
//! the requested flux sign alone is physically justified.

use super::super::CarreauFilmControls;
use super::{
    CarreauFilmError, CarreauShearedFilm, CarreauShearedFilmSpec, material_call, poll, positive,
};
use crate::film::FilmError;
use fs_material::fluid_transport::rheology::ThermalCarreauYasuda;
use fs_material::fluid_transport::rheology::slit::LaminarCarreauSlit;

/// Signed flow and fixed forcing, with an explicit positive thickness branch.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CarreauShearedFilmFlowSpec {
    pub density_kg_m3: f64,
    /// Signed requested volume flux per width [m²/s]; zero is permitted.
    pub flow_per_width_m2_s: f64,
    pub tangential_acceleration_m_s2: f64,
    pub surface_shear_pa: f64,
    /// Strictly ordered positive thickness bounds [m]. Both endpoints must
    /// satisfy the original material support. The interval must not cross
    /// `tau_s + rho*g*h = 0`; a turning point may be an endpoint.
    pub thickness_interval_m: [f64; 2],
}

/// One accepted film on the selected branch. Work and heat always use its
/// actual flux; the requested flow is retained separately with its residual.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CarreauShearedFilmFlowSolution {
    spec: CarreauShearedFilmFlowSpec,
    film: CarreauShearedFilm,
    profile_evaluations: u32,
    flow_thickness_derivative_m_s: f64,
}

impl CarreauShearedFilmFlowSolution {
    pub const fn spec(&self) -> CarreauShearedFilmFlowSpec {
        self.spec
    }
    pub const fn film(&self) -> &CarreauShearedFilm {
        &self.film
    }
    pub fn flow_residual_m2_s(&self) -> f64 {
        self.film.flow_per_width_m2_s() - self.spec.flow_per_width_m2_s
    }
    /// Includes both admitted endpoint profiles, not only Newton corrections.
    pub const fn profile_evaluations(&self) -> u32 {
        self.profile_evaluations
    }
    /// dq/dh = h*gamma_wall [m/s], retaining the selected branch's sign.
    pub const fn flow_thickness_derivative_m_s(&self) -> f64 {
        self.flow_thickness_derivative_m_s
    }
    /// dh/dq [s/m]. At a turning point this inverse derivative is singular;
    /// it is also absent when the nonzero reciprocal cannot be represented.
    pub fn thickness_flow_derivative_s_m(&self) -> Option<f64> {
        let derivative = 1.0 / self.flow_thickness_derivative_m_s;
        (derivative.is_finite() && derivative != 0.0).then_some(derivative)
    }
}

#[derive(Clone, Copy)]
struct Evaluated {
    film: CarreauShearedFilm,
    derivative: f64,
}

impl CarreauShearedFilm {
    pub fn at_flow(
        law: ThermalCarreauYasuda,
        temperature_k: f64,
        spec: CarreauShearedFilmFlowSpec,
        maximum_laminar_shear_rate_per_s: f64,
        controls: CarreauFilmControls,
    ) -> Result<CarreauShearedFilmFlowSolution, CarreauFilmError> {
        Self::at_flow_with_cancel(
            law,
            temperature_k,
            spec,
            maximum_laminar_shear_rate_per_s,
            controls,
            || false,
        )
    }

    /// Resolve thickness on the declared monotone interval. A single prepared
    /// source supplies every trial gap; the exact wall-rate derivative drives
    /// safeguarded Newton steps, with bisection at a turn or an unsafe proposal.
    /// No apparent-viscosity fit or implicit choice between two roots is made.
    ///
    /// Both endpoints count against `maximum_thickness_iterations`, which must
    /// be at least two. Acceptance requires the relative flow residual AND a
    /// relative thickness correction (or remaining bracket width) to satisfy
    /// `relative_flow_tolerance`. For counterflow, the flow scale is the larger
    /// of the requested flux and the native pressure/entrainment flux terms;
    /// this remains meaningful at exactly zero requested net flux.
    ///
    /// Source refusal, unsupported endpoints, cancellation and exhausted work
    /// never return a partially accepted thickness. A quiescent unforced film
    /// has no unique flux-to-thickness inverse and is refused.
    pub fn at_flow_with_cancel(
        law: ThermalCarreauYasuda,
        temperature_k: f64,
        spec: CarreauShearedFilmFlowSpec,
        maximum_laminar_shear_rate_per_s: f64,
        controls: CarreauFilmControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CarreauShearedFilmFlowSolution, CarreauFilmError> {
        poll(&mut cancelled)?;
        let [lower_h, upper_h] = spec.thickness_interval_m;
        if !positive(spec.density_kg_m3)
            || !spec.flow_per_width_m2_s.is_finite()
            || !spec.tangential_acceleration_m_s2.is_finite()
            || !spec.surface_shear_pa.is_finite()
            || !positive(lower_h)
            || !positive(upper_h)
            || lower_h >= upper_h
            || !positive(controls.relative_flow_tolerance)
            || controls.relative_flow_tolerance >= 1.0
            || !(2..=128).contains(&controls.maximum_thickness_iterations)
        {
            return Err(FilmError::InvalidInput(
                "Carreau film flux/thickness interval or controls",
            )
            .into());
        }
        let force = spec.density_kg_m3 * spec.tangential_acceleration_m_s2;
        if !force.is_finite() || (force == 0.0 && spec.tangential_acceleration_m_s2 != 0.0) {
            return Err(FilmError::NonRepresentable.into());
        }
        if force == 0.0 && spec.surface_shear_pa == 0.0 {
            return Err(FilmError::InvalidInput(
                "flux-to-thickness inversion requires tangential forcing",
            )
            .into());
        }
        let slit = LaminarCarreauSlit::new(
            law,
            temperature_k,
            1.0,
            lower_h,
            1.0,
            maximum_laminar_shear_rate_per_s,
        )?;
        let mut lower = evaluate(&slit, spec, lower_h, force, controls, &mut cancelled)?;
        let mut upper = evaluate(&slit, spec, upper_h, force, controls, &mut cancelled)?;
        let wall = [lower.film.wall_shear_pa(), upper.film.wall_shear_pa()];
        if (wall[0] < 0.0 && wall[1] > 0.0)
            || (wall[0] > 0.0 && wall[1] < 0.0)
            || (wall[0] == 0.0 && wall[1] == 0.0)
        {
            return Err(CarreauFilmError::NonMonotoneThicknessInterval {
                wall_shear_at_bounds_pa: wall,
            });
        }
        let increasing = wall[0] > 0.0 || wall[1] > 0.0;
        let flows = [
            lower.film.flow_per_width_m2_s(),
            upper.film.flow_per_width_m2_s(),
        ];
        if (increasing && flows[1] <= flows[0]) || (!increasing && flows[1] >= flows[0]) {
            return Err(FilmError::NonRepresentable.into());
        }
        let target = spec.flow_per_width_m2_s;
        if target < flows[0].min(flows[1]) || target > flows[0].max(flows[1]) {
            return Err(CarreauFilmError::FlowOutsideThicknessInterval {
                requested_m2_s: target,
                flow_at_bounds_m2_s: flows,
            });
        }
        let mut current = if (flows[0] - target).abs() <= (flows[1] - target).abs() {
            lower
        } else {
            upper
        };
        let mut profiles = 2;
        let mut previous_step = upper_h - lower_h;
        loop {
            poll(&mut cancelled)?;
            let actual = current.film.flow_per_width_m2_s();
            let residual = actual - target;
            let flow = current.film.profile.flow();
            let scale = target
                .abs()
                .max(flow.pressure_flow_m3_s.abs())
                .max(flow.entrained_flow_m3_s.abs());
            if !positive(scale) {
                return Err(FilmError::NonRepresentable.into());
            }
            let relative = if residual.is_finite() {
                residual / scale
            } else {
                actual / scale - target / scale
            };
            let correction = if residual == 0.0 {
                0.0
            } else {
                residual / current.derivative
            };
            let h = current.film.thickness_m();
            let width = upper.film.thickness_m() - lower.film.thickness_m();
            if residual == 0.0
                || (relative.abs() <= controls.relative_flow_tolerance
                    && (correction.abs() / h <= controls.relative_flow_tolerance
                        || width / h <= controls.relative_flow_tolerance))
            {
                poll(&mut cancelled)?;
                return Ok(CarreauShearedFilmFlowSolution {
                    spec,
                    film: current.film,
                    profile_evaluations: profiles,
                    flow_thickness_derivative_m_s: current.derivative,
                });
            }
            if profiles == controls.maximum_thickness_iterations {
                return Err(CarreauFilmError::UnresolvedFlow {
                    iterations: profiles,
                    residual_m2_s: residual,
                });
            }
            let proposal = h - correction;
            let lower_h = lower.film.thickness_m();
            let upper_h = upper.film.thickness_m();
            let next = if proposal.is_finite()
                && proposal > lower_h
                && proposal < upper_h
                && correction.abs() <= 0.5 * previous_step
            {
                proposal
            } else {
                f64::midpoint(lower_h, upper_h)
            };
            if next <= lower_h || next >= upper_h {
                return Err(CarreauFilmError::UnresolvedFlow {
                    iterations: profiles,
                    residual_m2_s: residual,
                });
            }
            previous_step = (next - h).abs();
            current = evaluate(&slit, spec, next, force, controls, &mut cancelled)?;
            profiles += 1;
            if (current.film.flow_per_width_m2_s() < target) == increasing {
                lower = current;
            } else {
                upper = current;
            }
        }
    }
}

fn evaluate(
    slit: &LaminarCarreauSlit,
    spec: CarreauShearedFilmFlowSpec,
    height: f64,
    force: f64,
    controls: CarreauFilmControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Evaluated, CarreauFilmError> {
    poll(cancelled)?;
    let passage = slit.with_gap(height)?;
    let profile = material_call(cancelled, |check| {
        passage.profile_with_upper_traction_check(
            force,
            0.0,
            spec.surface_shear_pa,
            controls.slit,
            check,
        )
    })?;
    let wall = material_call(cancelled, |check| profile.sample_with_check(0.0, check))?;
    let derivative = height * wall.velocity_gradient_per_s;
    if !derivative.is_finite() || (derivative == 0.0 && wall.velocity_gradient_per_s != 0.0) {
        return Err(FilmError::NonRepresentable.into());
    }
    let film = CarreauShearedFilm::from_profile(
        CarreauShearedFilmSpec {
            density_kg_m3: spec.density_kg_m3,
            thickness_m: height,
            tangential_acceleration_m_s2: spec.tangential_acceleration_m_s2,
            surface_shear_pa: spec.surface_shear_pa,
        },
        profile,
    )?;
    poll(cancelled)?;
    Ok(Evaluated { film, derivative })
}
