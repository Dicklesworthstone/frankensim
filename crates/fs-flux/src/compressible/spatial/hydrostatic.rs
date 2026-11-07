//! Axis-aligned pressure/gravity balance in the original volumetric gas solver.
//!
//! The thermodynamic owner integrates one axial reference, including variable
//! heat capacity; all three momenta and total energy still evolve on the native
//! Cartesian field. Reconstruction is relative to that reference, and its
//! physical face traction is removed before taking the momentum divergence.
//! Gravitational work uses the same mass flux as continuity, conserving gas
//! plus midpoint potential energy. A reference is neither a thermostat nor a
//! projection of the evolving state. No arbitrary/rotating equilibrium is claimed.
//!
//! This extends the existing planar implementation. Hydrostatic reconstruction
//! and a matching source are essential to resolving small pressure disturbances
//! (Käppeli & Mishra, JCP 259 (2014), doi:10.1016/j.jcp.2013.11.028).

use super::*;
use crate::compressible::hydrostatic::{acceleration, reference_position, reference_profile};
pub use fs_gas::dynamics::{
    AtmosphereQuadrature, HydrostaticReference, IsentropicAtmosphere, IsothermalAtmosphere,
    LapseRateAtmosphere,
};

pub(super) struct HydrostaticBalance {
    pub axis: usize,
    pub reference: HydrostaticReference,
    pub quadrature: AtmosphereQuadrature,
    cells: Vec<(State, Primitive)>,
    faces: Vec<(Primitive, f64)>,
}

impl HydrostaticBalance {
    fn new(
        gas: &FrozenGas,
        grid: &CartesianGrid,
        boundaries: Boundaries3d,
        axis: usize,
        reference: HydrostaticReference,
        quadrature: AtmosphereQuadrature,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        poll(cancelled)?;
        if axis >= 3 || grid.shape[axis] < 2 || boundaries.axis(axis) == AxisBoundaries::Periodic {
            return Err(Error::Invalid(
                "atmosphere axis needs at least two cells and physical ends",
            ));
        }
        let n = grid.shape[axis];
        // Share source-owned cell integration and physical face values with
        // 1-D and 2-D. Only O(n_axis) reference storage survives preparation.
        let source =
            reference_profile(gas, n, grid.spacing[axis], reference, quadrature, cancelled)?;
        let mut cells = allocated(n, (State::default(), Primitive::default()))?;
        let mut faces = allocated(n + 1, (Primitive::default(), 0.0))?;
        let primitive = |p: fs_gas::dynamics::Primitive| Primitive {
            density: p.density,
            velocity: [0.0; 3],
            pressure: p.pressure,
            temperature: p.temperature,
            sound_speed: p.sound_speed,
        };
        for (i, cell) in cells.iter_mut().enumerate() {
            poll(cancelled)?;
            let (q, p) = source.hydrostatic_reference_cell(i).expect("admitted cell");
            *cell = (
                State {
                    density: q.density,
                    momentum: [0.0; 3],
                    energy: q.energy,
                },
                primitive(p),
            );
        }
        for (i, face) in faces.iter_mut().enumerate() {
            poll(cancelled)?;
            let (_, p) = source.hydrostatic_reference_face(i).expect("admitted face");
            let trace = SpatialPrimitive {
                density: p.density,
                normal_velocity: 0.0,
                transverse_velocity: [0.0; 2],
                pressure: p.pressure,
            };
            // Preserve the exact floating-point path used by the actual face
            // solver rather than decoding pressure again during cancellation.
            let traction = gas.hllc_spatial_flux_from_primitive(trace, trace, &mut *cancelled)?;
            *face = (primitive(p), traction.normal_momentum);
        }
        Ok(Self {
            axis,
            reference,
            quadrature,
            cells,
            faces,
        })
    }

    pub fn cell(&self, grid: &CartesianGrid, cell: usize) -> (State, Primitive) {
        self.cells[grid.coordinates(cell)[self.axis]]
    }

    pub fn primitive(
        &self,
        gas: &FrozenGas,
        grid: &CartesianGrid,
        cell: usize,
        state: State,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Primitive, Error> {
        poll(cancelled)?;
        let (q, p) = self.cell(grid, cell);
        // Only cache thermodynamic recovery. Equilibrium cells still traverse
        // reconstruction, Riemann solves, conservative sources and both stages.
        let p = if state == q {
            p
        } else {
            state.primitive(gas, cancelled)?
        };
        self.variables(grid, cell, p)?;
        Ok(p)
    }

    pub fn variables(
        &self,
        grid: &CartesianGrid,
        cell: usize,
        p: Primitive,
    ) -> Result<[f64; 5], Error> {
        let reference = self.cell(grid, cell).1;
        let density = p.density / reference.density;
        let pressure = p.pressure / reference.pressure;
        if !positive(density) || !positive(pressure) {
            return Err(Error::Invalid(
                "unrepresentable relative atmosphere perturbation",
            ));
        }
        Ok([
            density,
            p.velocity[0],
            p.velocity[1],
            p.velocity[2],
            pressure,
        ])
    }

    pub fn trace(&self, grid: &CartesianGrid, cell: usize, axis: usize, sign: f64) -> [f64; 5] {
        let i = grid.coordinates(cell)[self.axis];
        if axis == self.axis {
            self.faces[i + usize::from(sign > 0.0)].0.array()
        } else {
            // A transverse face spans the same complete stratified cell.
            self.cells[i].1.array()
        }
    }

    pub fn traction(&self, grid: &CartesianGrid, face: usize) -> f64 {
        let mut shape = grid.shape;
        shape[self.axis] += 1;
        let stride = match self.axis {
            0 => 1,
            1 => shape[0],
            _ => shape[0] * shape[1],
        };
        self.faces[(face / stride) % shape[self.axis]].1
    }

    /// Perturbation force density, total force density, gravitational power density.
    pub fn sources(
        &self,
        grid: &CartesianGrid,
        cell: usize,
        state: State,
        fluxes: &[Vec<SpatialFaceFlux>; 3],
    ) -> [f64; 3] {
        let coordinates = grid.coordinates(cell);
        let i = coordinates[self.axis];
        let g = acceleration(self.reference);
        let perturbation = g * (state.density - self.cells[i].0.density);
        let force =
            (self.faces[i + 1].1 - self.faces[i].1) / grid.spacing[self.axis] + perturbation;
        let left = grid.cell_face(self.axis, coordinates);
        let right = left + grid.face_strides()[self.axis];
        // phi=-g*(x-x_ref); continuity and this work conserve sum(E+rho*phi)V.
        let power = g * (0.5 * fluxes[self.axis][left].mass + 0.5 * fluxes[self.axis][right].mass);
        [perturbation, force, power]
    }

    fn potential_energy(
        &self,
        grid: &CartesianGrid,
        states: &[State],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, Error> {
        let g = acceleration(self.reference);
        let datum = reference_position(self.reference);
        Ok(sum_arrays(
            states.iter().enumerate().map(|(cell, q)| {
                let x = grid.cell_center(cell).expect("admitted cell")[self.axis];
                [q.density * grid.cell_volume() * (-g * (x - datum))]
            }),
            cancelled,
        )?[0])
    }
}

impl Euler3d {
    /// Initialize integrated hydrostatic cell averages under one constant axial
    /// acceleration. The original gas source owns isothermal, isentropic and
    /// prescribed thermal-lapse profiles, including temperature-dependent Cp.
    ///
    /// `axis` is x=0, y=1 or z=2, with at least two cells and physical ends.
    /// Other directions may be periodic. Reflecting ends remain impermeable;
    /// transmissive ends extrapolate perturbations to the physical reference
    /// face. Every subsequent native `step` includes gravity and its work.
    /// This inviscid balance does not prevent separately modeled heat flow or
    /// externally imposed reservoirs from changing the reference equilibrium.
    #[allow(clippy::too_many_arguments)]
    pub fn new_hydrostatic_atmosphere(
        gas: FrozenGas,
        grid: CartesianGrid,
        axis: usize,
        reference: HydrostaticReference,
        quadrature: AtmosphereQuadrature,
        boundaries: Boundaries3d,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        let balance = HydrostaticBalance::new(
            &gas,
            &grid,
            boundaries,
            axis,
            reference,
            quadrature,
            &mut cancelled,
        )?;
        let mut states = allocated(grid.len(), State::default())?;
        for (cell, state) in states.iter_mut().enumerate() {
            poll(&mut cancelled)?;
            *state = balance.cell(&grid, cell).0;
        }
        if !balance
            .potential_energy(&grid, &states, &mut cancelled)?
            .is_finite()
        {
            return Err(Error::Invalid(
                "unrepresentable gravitational potential energy",
            ));
        }
        let mut flow = Self::new(gas, grid, states, boundaries, &mut cancelled)?;
        poll(&mut cancelled)?;
        flow.hydrostatic = Some(balance);
        Ok(flow)
    }

    /// Bind the same gravity/reconstruction to an authored initial disturbance
    /// without changing any cell average. Binding is atomic and occurs once,
    /// before evolution. Repeating the identical binding is idempotent.
    /// Existing complete-state reservoirs and pressure outlets retain their
    /// physical meaning; their boundary solve uses the reconstructed face.
    pub fn bind_hydrostatic_reference(
        &mut self,
        axis: usize,
        reference: HydrostaticReference,
        quadrature: AtmosphereQuadrature,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), Error> {
        poll(&mut cancelled)?;
        if self.hydrostatic.as_ref().is_some_and(|balance| {
            balance.axis == axis
                && balance.reference == reference
                && balance.quadrature == quadrature
        }) {
            return Ok(());
        }
        if self.hydrostatic.is_some() || self.time != 0.0 {
            return Err(Error::Invalid(
                "gravity reference must be bound once before evolution",
            ));
        }
        let balance = HydrostaticBalance::new(
            &self.gas,
            &self.grid,
            self.boundaries,
            axis,
            reference,
            quadrature,
            &mut cancelled,
        )?;
        for (cell, state) in self.states.iter().enumerate() {
            balance.primitive(&self.gas, &self.grid, cell, *state, &mut cancelled)?;
        }
        if !balance
            .potential_energy(&self.grid, &self.states, &mut cancelled)?
            .is_finite()
        {
            return Err(Error::Invalid(
                "unrepresentable gravitational potential energy",
            ));
        }
        poll(&mut cancelled)?;
        self.hydrostatic = Some(balance);
        Ok(())
    }

    /// Retained gravity axis and original thermodynamic reference.
    pub fn hydrostatic_reference(&self) -> Option<(usize, HydrostaticReference)> {
        self.hydrostatic
            .as_ref()
            .map(|balance| (balance.axis, balance.reference))
    }

    /// `sum rho_i V_i phi(x_i)` [J], zero without a retained reference.
    /// Gas plus this midpoint potential is conserved in a closed inviscid box.
    /// An open face also transports its potential times its outward mass flux.
    /// Midpoint potential quadrature is second order for a disturbed field.
    pub fn gravitational_potential_energy(&self) -> f64 {
        self.hydrostatic.as_ref().map_or(0.0, |balance| {
            balance
                .potential_energy(&self.grid, &self.states, &mut || false)
                .expect("infallible cancellation callback")
        })
    }
}

#[cfg(test)]
mod tests;
