//! Distributed mechanical forcing of the original resolved fluid.
//!
//! A force density f [N/m^3] enters the vector momentum equation directly.
//! The weak load is integral(N_i f dV), with no density or i*omega multiplier.
//! Its work is Re(integral(conj(v) . f dV))/2 for exp(+i omega t).
//! These are the volume-force terms of the linearized Navier--Stokes model:
//! <https://www.comsol.com/support/learning-center/course/modeling-speaker-drivers-in-comsol-multiphysics-202/supplement-a-equations-76561>.
//!
//! The original P2 velocity map includes its full moving-wall lifting. Both
//! free fluid equations and constrained structural equations receive the SAME
//! weak load, preserving work and the fluid traction on the solid. The fixed
//! equilibrium and force mechanism remain caller-owned; this adds no mean
//! flow, turbulence model, finite-amplitude momentum source or feedback law.

use super::super::helmholtz::heterogeneous::MAX_QUADRATURE_ORDER;
use super::*;

/// Largest number of simultaneously driven original volume-force profiles.
pub const MAX_BODY_FORCE_SOURCES: usize = 256;

/// Vector source shape [m^-3] driven by a scalar force amplitude [N].
///
/// Components use the original global Cartesian frame. Signed, nonunit and
/// zero-resultant distributions retain their authored magnitude and shape.
/// Separate profiles can carry independently phased force components. The
/// source borrows the exact fixed fluid/structure owner on which it was sampled.
pub struct ThermoviscousBodyForceSource<'a, 'c> {
    pub(super) model: &'a ThermoviscousFe<'c>,
    // Complete constrained momentum rows: Cartesian fluid then wall coordinates.
    pub(super) load: Vec<f64>,
    integral: [f64; 3],
}

impl ThermoviscousBodyForceSource<'_, '_> {
    /// Original dimensionless volume integral; multiply by the amplitude [N]
    /// for resultant force [N]. A zero resultant does not imply zero excitation.
    #[must_use]
    pub const fn integrated_shape(&self) -> [f64; 3] {
        self.integral
    }
}

impl<'c> ThermoviscousFe<'c> {
    /// Integrate an original vector force profile against the P2 velocity tests.
    ///
    /// The callback receives tetrahedron index and physical position [m]. It
    /// returns a vector density [m^-3]; its error is preserved. Quadrature is
    /// explicit and independently refinable from the original volume order
    /// through MAX_QUADRATURE_ORDER. Original curved positions and Jacobians
    /// are used. Each point is sampled once and warm solves never resample it.
    /// Stationary no-slip tests are eliminated; moving-wall tests receive the
    /// exact transpose of their original vector lifting, including tangential
    /// motion. No surface normal, clipping or normalization replaces the force.
    pub fn volume_force_source<'a, E: From<ThermoviscousError>>(
        &'a self,
        quadrature_order: usize,
        density_per_m3: &mut impl FnMut(usize, [f64; 3]) -> Result<[f64; 3], E>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ThermoviscousBodyForceSource<'a, 'c>, E> {
        poll(cancelled)?;
        if !(self.volume_quadrature_order()..=MAX_QUADRATURE_ORDER).contains(&quadrature_order) {
            return Err(
                ThermoviscousError::InvalidInput("volume force-source quadrature order").into(),
            );
        }
        let mut full = vec![[Dd::ZERO; 3]; self.space.ndof];
        let mut integral = [Dd::ZERO; 3];
        let quad = duffy_quadrature(quadrature_order);
        for (element, vertices) in self.space.complex.tets.iter().enumerate() {
            poll(cancelled)?;
            let dofs = self.space.element_dofs(element);
            for &(bary, weight) in &quad {
                poll(cancelled)?;
                let (position, jacobian) = if self.fe.is_curved() {
                    let point = self
                        .fe
                        .geometry_point(element, bary)
                        .map_err(ThermoviscousError::from)?;
                    (point.position_m, point.determinant_m3)
                } else {
                    let mut position = [0.0; 3];
                    for (&vertex, &value) in vertices.iter().zip(&bary) {
                        for (axis, coordinate) in position.iter_mut().enumerate() {
                            *coordinate += value * self.fe.positions()[vertex as usize][axis];
                        }
                    }
                    (position, 6.0 * self.geometry.vol_signed[element].abs())
                };
                let force = density_per_m3(element, position)?;
                let measure = weight * jacobian;
                if force.iter().any(|value| !value.is_finite())
                    || !measure.is_finite()
                    || measure <= 0.0
                {
                    return Err(ThermoviscousError::NonFinite.into());
                }
                let weighted = force.map(|value| Dd::from_f64(measure) * Dd::from_f64(value));
                for axis in 0..3 {
                    integral[axis] = integral[axis] + weighted[axis];
                }
                for &(global, shape) in &dofs {
                    let value = Dd::from_f64(shape.eval(bary, 2));
                    for axis in 0..3 {
                        full[global][axis] = full[global][axis] + weighted[axis] * value;
                    }
                }
            }
        }
        let mut load = vec![Dd::ZERO; self.pressure_start];
        for (global, force) in full.iter().enumerate() {
            poll(cancelled)?;
            if let Some(index) = self.velocity_map[global] {
                for axis in 0..3 {
                    load[axis * self.velocity_count + index] = force[axis];
                }
            }
            if let Some(wall) = &self.moving_wall {
                for (mode, shape) in wall.velocity_shapes.iter().enumerate() {
                    let row = 3 * self.velocity_count + mode;
                    for axis in 0..3 {
                        load[row] = load[row] + force[axis] * Dd::from_f64(shape[global][axis]);
                    }
                }
            }
        }
        let integral = integral.map(Dd::to_f64);
        let load: Vec<_> = load.into_iter().map(Dd::to_f64).collect();
        if integral.iter().chain(&load).any(|value| !value.is_finite()) {
            return Err(ThermoviscousError::NonFinite.into());
        }
        poll(cancelled)?;
        Ok(ThermoviscousBodyForceSource {
            model: self,
            load,
            integral,
        })
    }
}

pub(super) fn add_sources(
    model: &ThermoviscousFe<'_>,
    rhs: &mut [C64],
    sources: &[(&ThermoviscousBodyForceSource<'_, '_>, C64)],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), ThermoviscousError> {
    poll(cancelled)?;
    if sources.len() > MAX_BODY_FORCE_SOURCES {
        return Err(ThermoviscousError::Budget {
            what: "simultaneous volume force sources",
            requested: sources.len(),
            limit: MAX_BODY_FORCE_SOURCES,
        });
    }
    for &(source, amplitude) in sources {
        poll(cancelled)?;
        if !std::ptr::eq(source.model, model)
            || source.load.len() != model.pressure_start
            || !finite(amplitude)
        {
            return Err(ThermoviscousError::InvalidInput(
                "original force-source owner and finite force drive required",
            ));
        }
        for (row, &coefficient) in source.load.iter().enumerate() {
            rhs[row] = rhs[row] + amplitude.scale(coefficient);
        }
    }
    Ok(())
}

pub(super) fn add_work_and_reaction(
    model: &ThermoviscousFe<'_>,
    x: &[C64],
    sources: &[(&ThermoviscousBodyForceSource<'_, '_>, C64)],
    field: &mut ThermoviscousField,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), ThermoviscousError> {
    if sources.is_empty() {
        return Ok(());
    }
    let mut work = Dd::ZERO;
    let wall_start = 3 * model.velocity_count;
    let mut projected = vec![(Dd::ZERO, Dd::ZERO); model.pressure_start - wall_start];
    for &(source, amplitude) in sources {
        poll(cancelled)?;
        for (row, &coefficient) in source.load.iter().enumerate() {
            let force_re = Dd::from_f64(coefficient) * Dd::from_f64(amplitude.re);
            let force_im = Dd::from_f64(coefficient) * Dd::from_f64(amplitude.im);
            work = work + Dd::from_f64(x[row].re) * force_re + Dd::from_f64(x[row].im) * force_im;
            if row >= wall_start {
                projected[row - wall_start].0 = projected[row - wall_start].0 + force_re;
                projected[row - wall_start].1 = projected[row - wall_start].1 + force_im;
            }
        }
    }
    field.power.body_force_work_w = 0.5 * work.to_f64();
    field.power.input_w += field.power.body_force_work_w;
    // The original wall reaction is -(A_fluid*x - B_pressure*P). A volume
    // force adds its original constrained load inside that SAME residual.
    for (mode, force) in field
        .moving_wall
        .fluid_force_n_per_sqrt_kg
        .iter_mut()
        .enumerate()
    {
        let (re, im) = projected[mode];
        let added = C64::new(re.to_f64(), im.to_f64());
        *force = *force + added;
        field.moving_wall.power_into_fluid_w -=
            0.5 * (added.conj() * field.moving_wall.modal_velocity_m_sqrt_kg_per_s[mode]).re;
    }
    if !field.power.body_force_work_w.is_finite()
        || !field.power.input_w.is_finite()
        || !field.moving_wall.power_into_fluid_w.is_finite()
        || field
            .moving_wall
            .fluid_force_n_per_sqrt_kg
            .iter()
            .any(|&force| !finite(force))
    {
        return Err(ThermoviscousError::NonFinite);
    }
    poll(cancelled)
}
