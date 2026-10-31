//! A source-bound isentrope satisfies h(T)-h(T_ref)=g(x-x_ref).
//! Pressure follows the original entropy difference, not a frozen gamma law.

use super::*;

/// Adiabatically stratified reference with uniform frozen-mixture entropy.
///
/// The original heat capacity determines `dT/dx = g/Cp(T)`. This is a
/// convectively neutral reference for this frozen ideal gas, not a prescribed
/// constant lapse rate. Subsequent disturbances evolve through ordinary Euler
/// energy conservation; their entropy is not constrained to this reference.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IsentropicAtmosphere {
    /// Signed acceleration along the increasing coordinate [m/s²].
    pub acceleration: Acceleration,
    /// Temperature at the pressure datum, within the source interval [K].
    pub reference_temperature: Temperature,
    /// Absolute pressure at the datum [Pa].
    pub reference_pressure: Pressure,
    /// Datum position and zero of gravitational potential [m].
    pub reference_position: Length,
}

/// Bounded preparation of conserved cell averages along the source isentrope.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AtmosphereQuadrature {
    /// Relative change required between successive composite Gauss rules.
    /// This is a convergence estimate, not a rigorous integration error bound.
    pub relative_tolerance: f64,
    /// Maximum four-point Gauss panels per cell, in 2..=128.
    pub maximum_panels: usize,
}

impl Default for AtmosphereQuadrature {
    fn default() -> Self {
        Self {
            relative_tolerance: 1e-10,
            maximum_panels: 32,
        }
    }
}

impl Euler1d {
    /// Initialize and preserve an isentropic hydrostatic atmosphere.
    ///
    /// Temperature is found from the original source enthalpy difference and
    /// pressure from its entropy difference. Formation/entropy reference
    /// constants cancel before arithmetic. No constant-Cp or gamma fit is used.
    /// The source's bounded caloric inversion controls apply to every point.
    /// Conserved density and sensible energy are integrated, rather than formed
    /// from spatially averaged primitive fields. All preparation is cancellable.
    ///
    /// Constant area, consistent volumes and physical end boundaries are
    /// required, as for [`Self::new_isothermal_atmosphere`]. The retained face
    /// pressures and converged cell averages define the discrete reference;
    /// the original gravity-balanced reconstruction and conservative work then
    /// preserve it without any quadrature or allocation during time stepping.
    /// Nonrepresentable profiles, source-domain escape, or exhausted inversion
    /// or quadrature budgets refuse before an owner is returned.
    pub fn new_isentropic_atmosphere(
        gas: FrozenGas,
        mesh: Mesh1d,
        reference: IsentropicAtmosphere,
        quadrature: AtmosphereQuadrature,
        boundaries: Boundaries,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        poll(&mut cancelled)?;
        if boundaries == Boundaries::Periodic
            || !positive(quadrature.relative_tolerance)
            || quadrature.relative_tolerance > 1e-2
            || !(2..=128).contains(&quadrature.maximum_panels)
        {
            return Err(Error::Invalid(
                "isentropic atmosphere boundaries or quadrature",
            ));
        }
        admit_geometry(&mesh, &mut cancelled)?;
        let isentrope = SourceIsentrope::new(&gas, reference, &mut cancelled)?;
        let n = mesh.cell_count();
        let empty = (Conserved::default(), Primitive::default());
        let mut states = allocated(n, Conserved::default())?;
        let mut cells = allocated(n, empty)?;
        let mut faces = allocated(n + 1, empty)?;
        for (face, &x) in faces.iter_mut().zip(mesh.positions()) {
            *face = isentrope.at(x, &mut cancelled)?;
        }
        for i in 0..n {
            let q = isentrope.average(
                mesh.positions()[i],
                mesh.positions()[i + 1],
                quadrature,
                &mut cancelled,
            )?;
            let guess = 0.5 * (faces[i].1.temperature + faces[i + 1].1.temperature);
            let primitive = gas.decode(q, Some(guess), &mut cancelled)?;
            states[i] = q;
            cells[i] = (q, primitive);
        }
        let mut flow = Self::new_on_mesh(gas, mesh, states, boundaries, &mut cancelled)?;
        for (primitive, cell) in flow.work.primitive.iter_mut().zip(&cells) {
            *primitive = cell.1;
        }
        poll(&mut cancelled)?;
        flow.hydrostatic = Some(HydrostaticBalance {
            reference: HydrostaticReference::Isentropic(reference),
            cells,
            faces,
        });
        Ok(flow)
    }
}

struct SourceIsentrope<'a> {
    gas: &'a FrozenGas,
    reference: IsentropicAtmosphere,
    datum: (Conserved, Primitive),
    cp_reference: f64,
    minimum_enthalpy_change: f64,
    maximum_enthalpy_change: f64,
}

impl<'a> SourceIsentrope<'a> {
    fn new(
        gas: &'a FrozenGas,
        reference: IsentropicAtmosphere,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        poll(cancelled)?;
        let t = reference.reference_temperature.value();
        let p = reference.reference_pressure.value();
        if !reference.acceleration.value().is_finite()
            || !reference.reference_position.value().is_finite()
            || !positive(p)
            || !(gas.domain().minimum_temperature..=gas.domain().maximum_temperature).contains(&t)
        {
            return Err(Error::Invalid("isentropic atmosphere reference"));
        }
        let datum = gas.state_at(p / (gas.gas_constant() * t), 0.0, p)?;
        let cp_reference = gas
            .source()
            .evaluate_properties(
                reference.reference_temperature,
                reference.reference_pressure,
            )?
            .mass_specific_cp()
            .value();
        let mut result = Self {
            gas,
            reference,
            datum,
            cp_reference,
            minimum_enthalpy_change: 0.0,
            maximum_enthalpy_change: 0.0,
        };
        result.minimum_enthalpy_change =
            result.enthalpy_change(gas.domain().minimum_temperature)?;
        result.maximum_enthalpy_change =
            result.enthalpy_change(gas.domain().maximum_temperature)?;
        if !positive(cp_reference)
            || !result.minimum_enthalpy_change.is_finite()
            || !result.maximum_enthalpy_change.is_finite()
            || result.minimum_enthalpy_change >= result.maximum_enthalpy_change
        {
            return Err(Error::Invalid("unrepresentable isentropic source enthalpy"));
        }
        Ok(result)
    }

    fn enthalpy_change(&self, t: f64) -> Result<f64, Error> {
        let initial = self.reference.reference_temperature;
        Ok(self
            .gas
            .source()
            .specific_internal_energy_change(initial, Temperature::new(t))?
            .value()
            + self.gas.gas_constant() * (t - initial.value()))
    }

    fn temperature(&self, target: f64, cancelled: &mut impl FnMut() -> bool) -> Result<f64, Error> {
        if !target.is_finite()
            || target < self.minimum_enthalpy_change
            || target > self.maximum_enthalpy_change
        {
            return Err(Error::CaloricPath(
                "isentropic atmosphere leaves the source interval",
            ));
        }
        let t0 = self.reference.reference_temperature.value();
        if target == 0.0 {
            return Ok(t0);
        }
        let domain = self.gas.domain();
        let mut lower = domain.minimum_temperature;
        let mut upper = domain.maximum_temperature;
        let mut t = (t0 + target / self.cp_reference).clamp(lower, upper);
        // Prepare the reference near floating-point accuracy where the owner's
        // requested tolerance permits it. This avoids introducing root noise
        // into the independently refined cell quadrature.
        let tolerance = (domain.absolute_energy_tolerance
            + domain.relative_energy_tolerance * target.abs())
        .min(32.0 * f64::EPSILON * (self.cp_reference * t0 + target.abs()));
        for _ in 0..domain.maximum_iterations {
            poll(cancelled)?;
            let residual = self.enthalpy_change(t)? - target;
            if residual.abs() <= tolerance {
                return Ok(t);
            }
            if residual < 0.0 {
                lower = t;
            } else {
                upper = t;
            }
            let cp = self
                .gas
                .source()
                .evaluate_properties(Temperature::new(t), self.reference.reference_pressure)?
                .mass_specific_cp()
                .value();
            if !positive(cp) {
                return Err(Error::CaloricPath(
                    "isentropic atmosphere needs positive Cp",
                ));
            }
            let newton = t - residual / cp;
            t = if newton.is_finite() && newton > lower && newton < upper {
                newton
            } else {
                lower + 0.5 * (upper - lower)
            };
        }
        Err(Error::TemperatureIterations)
    }

    fn at(
        &self,
        x: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Conserved, Primitive), Error> {
        poll(cancelled)?;
        let target =
            self.reference.acceleration.value() * (x - self.reference.reference_position.value());
        let t = self.temperature(target, cancelled)?;
        if target == 0.0 {
            return Ok(self.datum);
        }
        let entropy = self
            .gas
            .source()
            .specific_entropy_change(
                self.reference.reference_temperature,
                self.reference.reference_pressure,
                Temperature::new(t),
                self.reference.reference_pressure,
            )?
            .value();
        let p = self.reference.reference_pressure.value()
            * fs_math::det::exp(entropy / self.gas.gas_constant());
        self.gas.state_at(p / (self.gas.gas_constant() * t), 0.0, p)
    }

    fn average(
        &self,
        left: f64,
        right: f64,
        control: AtmosphereQuadrature,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Conserved, Error> {
        if self.reference.acceleration.value() == 0.0 {
            poll(cancelled)?;
            return Ok(self.datum.0);
        }
        let mut panels = 1;
        let mut coarse = self.integrate(left, right, panels, cancelled)?;
        while panels < control.maximum_panels {
            panels = (2 * panels).min(control.maximum_panels);
            let fine = self.integrate(left, right, panels, cancelled)?;
            let density_tolerance = control.relative_tolerance * fine.density.abs();
            let energy_tolerance = control.relative_tolerance * fine.energy.abs()
                + fine.density * self.gas.domain().absolute_energy_tolerance;
            if (fine.density - coarse.density).abs() <= density_tolerance
                && (fine.energy - coarse.energy).abs() <= energy_tolerance
            {
                return Ok(fine);
            }
            coarse = fine;
        }
        Err(Error::Invalid(
            "isentropic atmosphere quadrature budget exhausted",
        ))
    }

    fn integrate(
        &self,
        left: f64,
        right: f64,
        panels: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Conserved, Error> {
        const NODES: [f64; 4] = [
            -0.861_136_311_594_052_6,
            -0.339_981_043_584_856_3,
            0.339_981_043_584_856_3,
            0.861_136_311_594_052_6,
        ];
        const WEIGHTS: [f64; 4] = [
            0.347_854_845_137_453_9,
            0.652_145_154_862_546_1,
            0.652_145_154_862_546_1,
            0.347_854_845_137_453_9,
        ];
        let mut sums = [0.0; 2];
        let mut corrections = [0.0; 2];
        let width = (right - left) / panels as f64;
        for panel in 0..panels {
            let middle = left + (panel as f64 + 0.5) * width;
            for (node, weight) in NODES.into_iter().zip(WEIGHTS) {
                let (q, _) = self.at(middle + 0.5 * width * node, cancelled)?;
                let factor = weight / (2 * panels) as f64;
                for (field, value) in [q.density, q.energy].into_iter().enumerate() {
                    let y = factor * value - corrections[field];
                    let next = sums[field] + y;
                    corrections[field] = (next - sums[field]) - y;
                    sums[field] = next;
                }
            }
        }
        Ok(Conserved {
            density: sums[0],
            momentum: 0.0,
            energy: sums[1],
        })
    }
}

#[cfg(test)]
mod tests;
