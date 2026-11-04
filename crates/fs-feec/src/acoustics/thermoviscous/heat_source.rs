//! Heat-release forcing of the ORIGINAL resolved fluid equations.
//!
//! A prescribed volumetric heat rate q' [W/m^3] enters the temperature equation
//! before eliminating either pressure or temperature. Thus expansion, thermal
//! diffusion, molecular relaxation, viscous losses and moving-wall reaction all
//! follow the existing mixed operator. For exp(+i omega t), acoustic thermal
//! work is Re(integral(conj(theta) q'))/2, theta=T'/T0. It is distinct from the
//! first-order heat supplied to the equilibrium fluid. Geometry and equilibrium
//! remain fixed; heat deposition and chemistry are inputs, not inferred models.

use super::super::helmholtz::heterogeneous::MAX_QUADRATURE_ORDER;
use super::*;

/// An original-volume spatial heat distribution, evaluated once during setup.
///
/// Source density is in m^-3, so a drive amplitude in W produces W/m^3.
/// The integral need not be one: no normalization, clipping or replacement of
/// a signed heating/cooling distribution is performed. The retained model
/// identity prevents reusing a source on a different mesh or thermal trace.
pub struct ThermoviscousHeatSource<'a, 'c> {
    model: &'a ThermoviscousFe<'c>,
    load: Vec<f64>,
    integral: f64,
}

impl ThermoviscousHeatSource<'_, '_> {
    /// Dimensionless integral of the original source shape over fluid volume.
    /// Multiply by the drive amplitude [W] for the heat-release phasor [W].
    #[must_use]
    pub fn integrated_shape(&self) -> f64 {
        self.integral
    }
}

impl<'c> ThermoviscousFe<'c> {
    /// Integrate q's shape against every original P2 temperature test function.
    ///
    /// The callback receives the original tetrahedron index and SI position;
    /// it returns a finite density [m^-3]. Quadrature is explicit, 4 through
    /// [`MAX_QUADRATURE_ORDER`], and is independently refinable for localized
    /// sources. Source errors and cancellation return without publishing a
    /// partial source. Isothermal constraints remain exact; their eliminated
    /// temperature tests receive no acoustic excitation.
    pub fn volume_heat_source<'a>(
        &'a self,
        quadrature_order: usize,
        density_per_m3: &mut impl FnMut(usize, [f64; 3]) -> Result<f64, ThermoviscousError>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ThermoviscousHeatSource<'a, 'c>, ThermoviscousError> {
        poll(cancelled)?;
        if !(4..=MAX_QUADRATURE_ORDER).contains(&quadrature_order) {
            return Err(ThermoviscousError::InvalidInput(
                "volume heat-source quadrature order",
            ));
        }
        let mut load = vec![0.0; self.n - self.temperature_start];
        let mut integral = Dd::from_f64(0.0);
        let quad = duffy_quadrature(quadrature_order);
        for (t, vertices) in self.space.complex.tets.iter().enumerate() {
            poll(cancelled)?;
            let dofs = self.space.element_dofs(t);
            let jacobian = 6.0 * self.geometry.vol_signed[t].abs();
            for &(bary, w) in &quad {
                poll(cancelled)?;
                let mut point = [0.0; 3];
                for (&vertex, &weight) in vertices.iter().zip(&bary) {
                    for (axis, coordinate) in point.iter_mut().enumerate() {
                        *coordinate += weight * self.fe.positions()[vertex as usize][axis];
                    }
                }
                let density = density_per_m3(t, point)?;
                let weight = density * w * jacobian;
                if !density.is_finite() || !weight.is_finite() {
                    return Err(ThermoviscousError::NonFinite);
                }
                integral = integral + Dd::from_f64(weight);
                for &(global, shape) in &dofs {
                    if let Some(thermal) = self.temperature_map[global] {
                        load[thermal] += weight * shape.eval(bary, 2);
                    }
                }
            }
        }
        let integral = integral.to_f64();
        if !integral.is_finite() || load.iter().any(|x| !x.is_finite()) {
            return Err(ThermoviscousError::NonFinite);
        }
        poll(cancelled)?;
        Ok(ThermoviscousHeatSource {
            model: self,
            load,
            integral,
        })
    }
}

impl PreparedThermoviscous<'_, '_> {
    /// Simultaneous traction, structural-force and distributed heat excitation.
    ///
    /// Heat amplitudes are phasors [W]. Sources must belong to this exact
    /// original FE owner; multiple independently phased sources superpose on
    /// the same full mixed LU. The fluid/solid response, losses and thermal
    /// source work use that SAME solution. No source sampling or factorization
    /// occurs here, so frequency sweeps reuse the original spatial source.
    pub fn response_with_heat_sources(
        &self,
        drive: &[C64],
        wall_forces: &[C64],
        heat_sources: &[(&ThermoviscousHeatSource<'_, '_>, C64)],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ThermoviscousField, ThermoviscousError> {
        poll(cancelled)?;
        let model = self.model;
        let wall_count = model.pressure_start - 3 * model.velocity_count;
        if drive.len() != model.loads.len()
            || wall_forces.len() != wall_count
            || drive.iter().chain(wall_forces).any(|z| !finite(*z))
            || heat_sources.iter().any(|(source, amplitude)| {
                !std::ptr::eq(source.model, model) || !finite(*amplitude)
            })
        {
            return Err(ThermoviscousError::InvalidInput(
                "complete finite source drives and original heat-source owner required",
            ));
        }
        let mut rhs = vec![C64::ZERO; model.n];
        for (load, &pressure) in model.loads.iter().zip(drive) {
            for (entry, &coefficient) in rhs.iter_mut().zip(load) {
                *entry = *entry + pressure.scale(coefficient);
            }
        }
        for (i, &force) in wall_forces.iter().enumerate() {
            let row = 3 * model.velocity_count + i;
            rhs[row] = rhs[row] + force;
        }
        for &(source, amplitude) in heat_sources {
            poll(cancelled)?;
            for (entry, &coefficient) in rhs[model.temperature_start..].iter_mut().zip(&source.load)
            {
                *entry = *entry + amplitude.scale(coefficient);
            }
        }
        let x = model.solve_system(
            &self.factor,
            &self.rows,
            &self.columns,
            std::f64::consts::TAU * self.frequency_hz,
            self.cp,
            &self.wall_admittances,
            None,
            &rhs,
            cancelled,
        )?;
        let mut field = self.field(drive, wall_forces, &x, cancelled)?;
        // Only heat loads occupy temperature rows. Keep the original assembled
        // source integral, not a loss remainder, as the acoustic work input.
        let mut heat_work = Dd::from_f64(0.0);
        for (&theta, &heat) in x[model.temperature_start..]
            .iter()
            .zip(&rhs[model.temperature_start..])
        {
            heat_work = heat_work
                + Dd::from_f64(theta.re) * Dd::from_f64(heat.re)
                + Dd::from_f64(theta.im) * Dd::from_f64(heat.im);
        }
        field.power.heat_source_work_w = 0.5 * heat_work.to_f64();
        field.power.input_w += field.power.heat_source_work_w;
        if !field.power.heat_source_work_w.is_finite() || !field.power.input_w.is_finite() {
            return Err(ThermoviscousError::NonFinite);
        }
        poll(cancelled)?;
        Ok(field)
    }
}
