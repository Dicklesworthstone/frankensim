//! Hydrostatic reconstruction for a frozen, thermally perfect gas.
//!
//! This is a reference-state well-balanced finite-volume method: reconstruct
//! perturbations about a source-integrated atmosphere, pair its pressure
//! traction with gravity, and transport gravitational work with the mass flux.
//! The evolving gas remains adiabatic and uses the original NASA calorics.
//! See Käppeli and Mishra, JCP 259 (2014), doi:10.1016/j.jcp.2013.11.028,
//! for hydrostatic reconstruction; this bounded implementation preserves only
//! the declared isothermal or isentropic, constant-acceleration, constant-area
//! reference. Isentropes retain the original variable heat capacity.

mod isentropic;
pub use isentropic::{AtmosphereQuadrature, IsentropicAtmosphere};

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

/// Thermodynamic family of the retained hydrostatic reference.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum HydrostaticReference {
    /// Uniform temperature, with exponential pressure and density.
    Isothermal(IsothermalAtmosphere),
    /// Uniform specific entropy, with source-dependent temperature lapse.
    Isentropic(IsentropicAtmosphere),
}

impl HydrostaticReference {
    fn acceleration(self) -> f64 {
        match self {
            Self::Isothermal(reference) => reference.acceleration.value(),
            Self::Isentropic(reference) => reference.acceleration.value(),
        }
    }

    fn reference_position(self) -> f64 {
        match self {
            Self::Isothermal(reference) => reference.reference_position.value(),
            Self::Isentropic(reference) => reference.reference_position.value(),
        }
    }
}

pub(super) struct HydrostaticBalance {
    pub(super) reference: HydrostaticReference,
    pub(super) cells: Vec<(Conserved, Primitive)>,
    pub(super) faces: Vec<(Conserved, Primitive)>,
}

impl HydrostaticBalance {
    pub(super) fn perturbation_force(&self, i: usize, q: Conserved, mesh: &Mesh1d) -> f64 {
        self.reference.acceleration() * (q.density - self.cells[i].0.density) * mesh.volumes()[i]
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
        self.reference.acceleration()
            * (mesh.face_areas()[i] * (center - mesh.positions()[i]) * faces[i].mass
                + mesh.face_areas()[i + 1] * (mesh.positions()[i + 1] - center) * faces[i + 1].mass)
    }

    pub(super) fn potential(&self, x: f64) -> f64 {
        -self.reference.acceleration() * (x - self.reference.reference_position())
    }

    /// The same zero-slope physical-face reconstruction used by a passage end.
    /// A compact-volume port must not substitute the end cell's volume average.
    pub(super) fn end_trace(
        &self,
        gas: &FrozenGas,
        state: Conserved,
        primitive: Primitive,
        cell: usize,
        face: usize,
    ) -> Result<(Conserved, Primitive), Error> {
        if state == self.cells[cell].0 {
            return Ok(self.faces[face]);
        }
        let reference = self.cells[cell].1;
        let boundary = self.faces[face].1;
        gas.state_at(
            boundary.density + (primitive.density - reference.density),
            boundary.velocity + (primitive.velocity - reference.velocity),
            boundary.pressure + (primitive.pressure - reference.pressure),
        )
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
    /// Gas networks require their explicit shared gravitational-potential
    /// declaration. Arbitrary caller-owned end ports still require a separate
    /// potential convention and remain unsupported by `step_with_ports`.
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
        admit_geometry(&mesh, &mut cancelled)?;
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
            reference: HydrostaticReference::Isothermal(reference),
            cells,
            faces,
        });
        Ok(flow)
    }

    /// Retained isothermal reference, if this is an isothermal atmosphere.
    /// Use [`Self::hydrostatic_reference`] to inspect either supported family.
    pub fn atmosphere(&self) -> Option<IsothermalAtmosphere> {
        match self.hydrostatic_reference() {
            Some(HydrostaticReference::Isothermal(reference)) => Some(reference),
            _ => None,
        }
    }

    /// Retained physical reference whenever hydrostatic gravity is enabled.
    pub fn hydrostatic_reference(&self) -> Option<HydrostaticReference> {
        self.hydrostatic.as_ref().map(|balance| balance.reference)
    }

    /// Original reference atmosphere at a physical face, with its cached source
    /// evaluation. This initializes a matching compact plenum without replacing
    /// the original caloric law or mistaking a cell average for a face state.
    /// Returns `None` without an atmosphere or for an out-of-range face index.
    pub fn hydrostatic_reference_face(&self, face: usize) -> Option<(Conserved, Primitive)> {
        self.hydrostatic.as_ref()?.faces.get(face).copied()
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

fn admit_geometry(mesh: &Mesh1d, cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    let area = mesh.uniform_cross_section().ok_or(Error::Invalid(
        "hydrostatic atmosphere requires a constant-area passage",
    ))?;
    for i in 0..mesh.cell_count() {
        poll(cancelled)?;
        let volume = area * mesh.width(i);
        if !positive(volume)
            || (volume - mesh.volumes()[i]).abs()
                > 64.0 * f64::EPSILON * volume.max(mesh.volumes()[i])
        {
            return Err(Error::Invalid(
                "atmosphere volume disagrees with constant-area geometry",
            ));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests;
