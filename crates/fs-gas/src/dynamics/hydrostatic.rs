//! Isothermal hydrostatic reconstruction for a frozen, thermally perfect gas.
//!
//! This is a reference-state well-balanced finite-volume method: reconstruct
//! perturbations about an analytically integrated atmosphere, pair its pressure
//! traction with gravity, and transport gravitational work with the mass flux.
//! The evolving gas remains adiabatic and uses the original NASA calorics.
//! See Käppeli and Mishra, JCP 259 (2014), doi:10.1016/j.jcp.2013.11.028,
//! for hydrostatic reconstruction; this bounded implementation preserves only
//! the declared isothermal, constant-acceleration, constant-area reference.

use super::{
    Boundaries, Conserved, Error, Euler1d, FaceFlux, FrozenGas, Mesh1d, Primitive, allocated, poll,
    positive,
};
use fs_qty::{Acceleration, Length, Pressure, Temperature};

/// Physical reference atmosphere used for gravity-balanced reconstruction.
///
/// `p(x) = p_ref exp[g (x-x_ref)/(R T_ref)]`, with signed axial acceleration
/// `g` and the original frozen mixture's gas constant `R`. Temperature is a
/// reference profile, not a thermostat: disturbed gas evolves its total energy
/// and original variable heat capacity. No constant-gamma law is introduced.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IsothermalAtmosphere {
    /// Signed acceleration along the increasing mesh coordinate [m/s²].
    pub acceleration: Acceleration,
    /// Uniform reference temperature inside the original source support [K].
    pub temperature: Temperature,
    /// Absolute pressure at `reference_position` [Pa].
    pub reference_pressure: Pressure,
    /// Physical position of the pressure datum and zero of potential [m].
    pub reference_position: Length,
}

pub(super) struct HydrostaticBalance {
    pub(super) reference: IsothermalAtmosphere,
    pub(super) cells: Vec<(Conserved, Primitive)>,
    pub(super) faces: Vec<(Conserved, Primitive)>,
}

impl HydrostaticBalance {
    pub(super) fn perturbation_force(&self, i: usize, q: Conserved, mesh: &Mesh1d) -> f64 {
        self.reference.acceleration.value()
            * (q.density - self.cells[i].0.density)
            * mesh.volumes()[i]
    }

    pub(super) fn force(&self, i: usize, q: Conserved, mesh: &Mesh1d) -> f64 {
        mesh.face_areas()[i] * (self.faces[i + 1].1.pressure - self.faces[i].1.pressure)
            + self.perturbation_force(i, q, mesh)
    }

    pub(super) fn power(&self, i: usize, mesh: &Mesh1d, faces: &[FaceFlux]) -> f64 {
        // With phi=-g*(x-x_ref), this is -sum_faces (phi_face-phi_i)*mass_out.
        // Therefore sum_i V_i (E_i + phi_i*rho_i) changes only by boundary
        // transport, for both RK stages and their convex combination.
        let center = mesh.center(i);
        self.reference.acceleration.value()
            * (mesh.face_areas()[i] * (center - mesh.positions()[i]) * faces[i].mass
                + mesh.face_areas()[i + 1] * (mesh.positions()[i + 1] - center) * faces[i + 1].mass)
    }

    fn potential(&self, x: f64) -> f64 {
        -self.reference.acceleration.value() * (x - self.reference.reference_position.value())
    }
}

impl Euler1d {
    /// Initialize an analytically integrated atmosphere and retain its gravity.
    ///
    /// The mesh may be nonuniform, but must have a constant area and volumes
    /// consistent with that area. Each initial density and energy is the exact
    /// exponential cell average up to floating-point/source arithmetic. The
    /// reference faces are resolved once and reused during subsequent steps.
    /// MUSCL reconstructs deviations from these averages; even first-order
    /// reconstruction preserves the resting atmosphere to roundoff.
    ///
    /// Physical end boundaries are required. A transmissive end extrapolates
    /// the perturbation while continuing the hydrostatic profile to the face.
    /// Reservoir conditions refer to the physical end face. Source updates can
    /// perturb the initial state before ordinary [`Self::step`] evolution.
    /// End ports and gas-network assembly currently refuse this owner because
    /// their shared gravitational-potential datum has not yet been declared.
    pub fn new_isothermal_atmosphere(
        gas: FrozenGas,
        mesh: Mesh1d,
        reference: IsothermalAtmosphere,
        boundaries: Boundaries,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        poll(&mut cancelled)?;
        let g = reference.acceleration.value();
        let temperature = reference.temperature.value();
        let datum = reference.reference_position.value();
        let pressure = reference.reference_pressure.value();
        let rt = gas.gas_constant() * temperature;
        if !g.is_finite()
            || !datum.is_finite()
            || !positive(pressure)
            || !positive(rt)
            || !(gas.domain().minimum_temperature..=gas.domain().maximum_temperature)
                .contains(&temperature)
            || boundaries == Boundaries::Periodic
        {
            return Err(Error::Invalid(
                "hydrostatic reference or physical boundaries",
            ));
        }
        let area = mesh.uniform_cross_section().ok_or(Error::Invalid(
            "isothermal atmosphere requires a constant-area passage",
        ))?;
        let n = mesh.cell_count();
        let mut states = allocated(n, Conserved::default())?;
        let empty = (Conserved::default(), Primitive::default());
        let mut cells = allocated(n, empty)?;
        let mut faces = allocated(n + 1, empty)?;
        let beta = g / rt;
        for (i, &position) in mesh.positions().iter().enumerate() {
            poll(&mut cancelled)?;
            let exponent = beta * (position - datum);
            let potential = -g * (position - datum);
            let p = pressure * fs_math::det::exp(exponent);
            if !exponent.is_finite() || !potential.is_finite() || !positive(p) {
                return Err(Error::Invalid(
                    "unrepresentable hydrostatic pressure or potential",
                ));
            }
            faces[i] = gas.state_at(p / rt, 0.0, p)?;
        }
        for i in 0..n {
            poll(&mut cancelled)?;
            let volume = area * mesh.width(i);
            if !positive(volume)
                || (volume - mesh.volumes()[i]).abs()
                    > 64.0 * f64::EPSILON * volume.max(mesh.volumes()[i])
            {
                return Err(Error::Invalid(
                    "atmosphere volume disagrees with constant-area geometry",
                ));
            }
            // The average of exp(beta*x) is evaluated without subtracting
            // nearly equal face pressures. The series is even and regular at g=0.
            let z = 0.5 * beta * mesh.width(i);
            let average_factor = if z.abs() < 1e-3 {
                let z2 = z * z;
                1.0 + z2 * (1.0 / 6.0 + z2 * (1.0 / 120.0 + z2 / 5040.0))
            } else {
                (fs_math::det::exp(z) - fs_math::det::exp(-z)) / (2.0 * z)
            };
            let average_pressure =
                pressure * fs_math::det::exp(beta * (mesh.center(i) - datum)) * average_factor;
            let state = gas.conserved(average_pressure / rt, 0.0, average_pressure)?;
            let primitive = gas.decode(state, Some(temperature), &mut cancelled)?;
            states[i] = state;
            cells[i] = (state, primitive);
        }
        let mut flow = Self::new_on_mesh(gas, mesh, states, boundaries, &mut cancelled)?;
        for (primitive, cell) in flow.work.primitive.iter_mut().zip(&cells) {
            *primitive = cell.1;
        }
        poll(&mut cancelled)?;
        flow.hydrostatic = Some(HydrostaticBalance {
            reference,
            cells,
            faces,
        });
        Ok(flow)
    }

    /// Retained physical hydrostatic reference, when gravity is enabled.
    pub fn atmosphere(&self) -> Option<IsothermalAtmosphere> {
        self.hydrostatic.as_ref().map(|balance| balance.reference)
    }

    /// Discrete gravitational potential energy, `sum rho_i V_i phi(x_i)` [J].
    ///
    /// Adding this to [`Self::totals`]' referenced energy gives the conserved
    /// gas-plus-potential coordinate. Physical energy additionally restores the
    /// original gas formation reference. This midpoint potential quadrature is
    /// second order, not the exact integral of an arbitrary within-cell field.
    pub fn gravitational_potential_energy(&self) -> f64 {
        let Some(balance) = &self.hydrostatic else {
            return 0.0;
        };
        let mut sum = 0.0;
        let mut correction = 0.0;
        for (i, q) in self.states.iter().enumerate() {
            let term = q.density * self.mesh.volumes()[i] * balance.potential(self.mesh.center(i));
            let y = term - correction;
            let next = sum + y;
            correction = (next - sum) - y;
            sum = next;
        }
        sum
    }
}

#[cfg(test)]
mod tests;
