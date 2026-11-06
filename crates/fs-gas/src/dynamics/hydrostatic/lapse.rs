//! Prescribed thermal stratification in hydrostatic mechanical equilibrium.
//!
//! For T(x)=T0+L(x-x0), integrating dp/dx=rho*g with p=rho*R*T gives
//! ln(p/p0)=g*(x-x0)/(R*T0) * log1p(z)/z, z=L*(x-x0)/T0.
//! The regular small-z form includes the isothermal limit without dividing by
//! a small lapse rate. Cp(T) remains the original source law in conserved
//! energy, acoustic response and convective stability. This is one dry, frozen
//! ideal-gas layer, not moist convection or chemical/diffusive equilibrium.

use super::*;

/// A linear temperature profile with an independently specified gravity.
///
/// Unlike the source isentrope, this reference can be stably or unstably
/// stratified. It is mechanical equilibrium of the inviscid adiabatic Euler
/// equations; a nonzero lapse is not thermal equilibrium under conduction.
/// Temperature is an initial reference, never a thermostat during evolution.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LapseRateAtmosphere {
    /// Signed acceleration along the increasing coordinate [m/s²].
    pub acceleration: Acceleration,
    /// Temperature at the pressure datum [K].
    pub reference_temperature: Temperature,
    /// Absolute pressure at the datum [Pa].
    pub reference_pressure: Pressure,
    /// Datum position and zero of gravitational potential [m].
    pub reference_position: Length,
    /// Signed dT/dx [K/m], with either sign allowed.
    pub temperature_gradient_k_m: f64,
}

impl LapseRateAtmosphere {
    fn at(
        self,
        gas: &FrozenGas,
        position: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Conserved, Primitive), Error> {
        poll(cancelled)?;
        let t0 = self.reference_temperature.value();
        let p0 = self.reference_pressure.value();
        let x0 = self.reference_position.value();
        let g = self.acceleration.value();
        let lapse = self.temperature_gradient_k_m;
        let domain = gas.domain();
        if !g.is_finite()
            || !x0.is_finite()
            || !position.is_finite()
            || !lapse.is_finite()
            || !positive(p0)
            || !(domain.minimum_temperature..=domain.maximum_temperature).contains(&t0)
        {
            return Err(Error::Invalid("lapse-rate atmosphere reference"));
        }
        let dx = position - x0;
        let dt = lapse * dx;
        let temperature = t0 + dt;
        if !dx.is_finite()
            || !dt.is_finite()
            || !(domain.minimum_temperature..=domain.maximum_temperature).contains(&temperature)
        {
            return Err(Error::CaloricPath(
                "lapse-rate atmosphere leaves the original source interval",
            ));
        }
        let potential = -g * dx;
        let z = dt / t0;
        // Alternating Taylor series through z^8: at |z|<0.01 the omitted
        // term is <1.01e-19. This avoids a rounded 1+z destroying small lapses.
        let log_ratio = if z.abs() < 0.01 {
            let mut ratio = 1.0 / 9.0;
            for coefficient in [
                -0.125,
                1.0 / 7.0,
                -1.0 / 6.0,
                0.2,
                -0.25,
                1.0 / 3.0,
                -0.5,
                1.0,
            ] {
                ratio = coefficient + z * ratio;
            }
            ratio
        } else {
            fs_math::det::ln(temperature / t0) / z
        };
        let exponent = -potential / (gas.gas_constant() * t0) * log_ratio;
        let pressure = p0 * fs_math::det::exp(exponent);
        if !potential.is_finite() || !exponent.is_finite() || !positive(pressure) {
            return Err(Error::Invalid(
                "unrepresentable lapse-rate pressure or gravitational potential",
            ));
        }
        gas.state_at(pressure / (gas.gas_constant() * temperature), 0.0, pressure)
    }

    /// Local squared buoyancy frequency [s^-2] from the original source Cp(T).
    ///
    /// N² = (g/T) [g/Cp(T) - dT/dx]. Positive is locally restoring and
    /// negative is convectively unstable for a dry frozen adiabatic parcel.
    /// A negative result is physical, not clipped or refused. This local
    /// diagnostic is not a stability certificate for a discrete flow or a
    /// rotating, moist, reacting or thermally conducting atmosphere.
    pub fn buoyancy_frequency_squared(
        self,
        gas: &FrozenGas,
        position: Length,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<f64, Error> {
        let (_, point) = self.at(gas, position.value(), &mut cancelled)?;
        let cp = gas
            .source()
            .evaluate_properties(
                Temperature::new(point.temperature),
                Pressure::new(point.pressure),
            )?
            .mass_specific_cp()
            .value();
        let g = self.acceleration.value();
        let result = (g / point.temperature) * (g / cp - self.temperature_gradient_k_m);
        if !positive(cp) || !result.is_finite() {
            return Err(Error::CaloricPath(
                "unrepresentable lapse-rate buoyancy response",
            ));
        }
        poll(&mut cancelled)?;
        Ok(result)
    }
}

impl Euler1d {
    /// Prepare a prescribed thermal lapse and retain its hydrostatic gravity.
    ///
    /// Pressure is integrated analytically. Conserved density and source
    /// sensible energy use the same bounded quadrature as the source isentrope;
    /// cell averages are never replaced by primitive samples at cell centres.
    /// Zero lapse reuses the exact exponential isothermal integrals. Reference
    /// profiles are cached once; ordinary stepping adds no quadrature or
    /// profile evaluation. Nonuniform constant-area meshes are supported.
    ///
    /// Stable and unstable references both start at mechanical rest. Actual
    /// perturbations evolve with the original gas law, all pressure forces and
    /// conservative gravitational work; they are not reset to this profile.
    pub fn new_lapse_rate_atmosphere(
        gas: FrozenGas,
        mesh: Mesh1d,
        reference: LapseRateAtmosphere,
        quadrature: AtmosphereQuadrature,
        boundaries: Boundaries,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        poll(&mut cancelled)?;
        quadrature.admit()?;
        if boundaries == Boundaries::Periodic {
            return Err(Error::Invalid(
                "lapse-rate atmosphere requires physical ends",
            ));
        }
        admit_geometry(&mesh, &mut cancelled)?;
        reference.at(&gas, reference.reference_position.value(), &mut cancelled)?;
        if reference.temperature_gradient_k_m == 0.0 {
            let mut flow = Self::new_isothermal_atmosphere(
                gas,
                mesh,
                IsothermalAtmosphere {
                    acceleration: reference.acceleration,
                    temperature: reference.reference_temperature,
                    reference_pressure: reference.reference_pressure,
                    reference_position: reference.reference_position,
                },
                boundaries,
                &mut cancelled,
            )?;
            flow.hydrostatic
                .as_mut()
                .expect("prepared atmosphere")
                .reference = HydrostaticReference::LapseRate(reference);
            return Ok(flow);
        }
        let n = mesh.cell_count();
        let empty = (Conserved::default(), Primitive::default());
        let mut states = allocated(n, Conserved::default())?;
        let mut cells = allocated(n, empty)?;
        let mut faces = allocated(n + 1, empty)?;
        // Linear temperature extrema occur at the faces. Checking every face
        // before quadrature prevents an unadmitted endpoint being skipped.
        for (face, &x) in faces.iter_mut().zip(mesh.positions()) {
            *face = reference.at(&gas, x, &mut cancelled)?;
        }
        for i in 0..n {
            let q = quadrature::average(
                &gas,
                mesh.positions()[i],
                mesh.positions()[i + 1],
                quadrature,
                &mut |x, cancelled| reference.at(&gas, x, cancelled),
                &mut cancelled,
            )?;
            let guess = 0.5 * faces[i].1.temperature + 0.5 * faces[i + 1].1.temperature;
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
            reference: HydrostaticReference::LapseRate(reference),
            cells,
            faces,
        });
        Ok(flow)
    }
}

#[cfg(test)]
mod tests;
