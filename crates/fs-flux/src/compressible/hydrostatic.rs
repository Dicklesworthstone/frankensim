//! Axis-aligned hydrostatic reconstruction over the original gas reference.
//!
//! Isothermal, isentropic and thermal-lapse references share the original gas integration.
//! Perturbations evolve adiabatically with both momentum components and the
//! retained variable-Cp caloric source; the reference is not a thermostat.
//! Gravity work uses the shared mass flux, conserving gas plus midpoint
//! potential energy. This does not claim arbitrary equilibria or rotating flow.

use super::*;
pub use fs_gas::dynamics::{
    AtmosphereQuadrature, HydrostaticReference, IsentropicAtmosphere, IsothermalAtmosphere,
    LapseRateAtmosphere,
};
use fs_gas::dynamics::{Boundaries, Boundary, Euler1d, Mesh1d};
use fs_qty::Length;

/// Prepare one compact axial profile in the thermodynamic owner. Both Cartesian
/// dimensions reuse its integrated conserved averages and physical face states.
pub(super) fn reference_profile(
    gas: &FrozenGas,
    cells: usize,
    spacing: f64,
    reference: HydrostaticReference,
    quadrature: AtmosphereQuadrature,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Euler1d, Error> {
    let mesh = Mesh1d::uniform(cells, Length::new(spacing), 1.0)?;
    let ends = Boundaries::Ends {
        left: Boundary::Reflecting,
        right: Boundary::Reflecting,
    };
    Ok(match reference {
        HydrostaticReference::Isothermal(value) => {
            Euler1d::new_isothermal_atmosphere(gas.clone(), mesh, value, ends, cancelled)?
        }
        HydrostaticReference::Isentropic(value) => Euler1d::new_isentropic_atmosphere(
            gas.clone(),
            mesh,
            value,
            quadrature,
            ends,
            cancelled,
        )?,
        HydrostaticReference::LapseRate(value) => Euler1d::new_lapse_rate_atmosphere(
            gas.clone(),
            mesh,
            value,
            quadrature,
            ends,
            cancelled,
        )?,
    })
}

pub(super) fn acceleration(reference: HydrostaticReference) -> f64 {
    match reference {
        HydrostaticReference::Isothermal(value) => value.acceleration.value(),
        HydrostaticReference::Isentropic(value) => value.acceleration.value(),
        HydrostaticReference::LapseRate(value) => value.acceleration.value(),
    }
}

pub(super) fn reference_position(reference: HydrostaticReference) -> f64 {
    match reference {
        HydrostaticReference::Isothermal(value) => value.reference_position.value(),
        HydrostaticReference::Isentropic(value) => value.reference_position.value(),
        HydrostaticReference::LapseRate(value) => value.reference_position.value(),
    }
}

pub(super) struct HydrostaticBalance {
    pub axis: usize,
    pub reference: HydrostaticReference,
    pub quadrature: AtmosphereQuadrature,
    cells: Vec<(State, Primitive)>,
    // Physical-face primitive and its source-resolved normal pressure flux.
    // Retaining that actual traction also avoids a second EOS roundoff in
    // subtracting the reference from the evolving numerical face flux.
    faces: Vec<(Primitive, f64)>,
}

impl HydrostaticBalance {
    fn new(
        gas: &FrozenGas,
        grid: &CartesianGrid,
        boundaries: Boundaries2d,
        axis: usize,
        reference: HydrostaticReference,
        quadrature: AtmosphereQuadrature,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        poll(cancelled)?;
        if axis >= 2 || grid.shape[axis] < 2 || boundaries.axis(axis) == AxisBoundaries::Periodic {
            return Err(Error::Invalid(
                "atmosphere axis needs at least two cells and physical ends",
            ));
        }
        let n = grid.shape[axis];
        // Integrate the profile once in the existing thermodynamic owner.
        // Its temporary 1-D evolution workspace is discarded, not retained or
        // stepped alongside the planar solver. Only O(n_axis) data survive.
        let source = reference_profile(
            gas,
            n,
            grid.spacing[axis],
            reference,
            quadrature,
            &mut *cancelled,
        )?;
        let mut cells = allocated(n, (State::default(), Primitive::default()))?;
        let mut faces = allocated(n + 1, (Primitive::default(), 0.0))?;
        let primitive = |p: fs_gas::dynamics::Primitive| Primitive {
            density: p.density,
            velocity: [0.0; 2],
            pressure: p.pressure,
            temperature: p.temperature,
            sound_speed: p.sound_speed,
        };
        for (i, cell) in cells.iter_mut().enumerate() {
            poll(cancelled)?;
            let (q, p) = source
                .hydrostatic_reference_cell(i)
                .expect("admitted reference cell");
            *cell = (
                State {
                    density: q.density,
                    momentum: [0.0; 2],
                    energy: q.energy,
                },
                primitive(p),
            );
        }
        for (i, face) in faces.iter_mut().enumerate() {
            poll(cancelled)?;
            let (_, p) = source
                .hydrostatic_reference_face(i)
                .expect("admitted reference face");
            let trace = NormalPrimitive {
                density: p.density,
                normal_velocity: 0.0,
                transverse_velocity: 0.0,
                pressure: p.pressure,
            };
            let traction = gas.hllc_normal_flux_from_primitive(trace, trace, &mut *cancelled)?;
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

    fn position(&self, grid: &CartesianGrid, cell: usize) -> usize {
        if self.axis == 0 {
            cell % grid.shape[0]
        } else {
            cell / grid.shape[0]
        }
    }

    pub fn cell(&self, grid: &CartesianGrid, cell: usize) -> (State, Primitive) {
        self.cells[self.position(grid, cell)]
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
        // Cache only thermodynamic decoding; every cell still runs the full
        // reconstruction, Riemann, conservative source and RK update.
        let primitive = if state == q {
            p
        } else {
            state.primitive(gas, cancelled)?
        };
        self.variables(grid, cell, primitive)?;
        Ok(primitive)
    }

    /// Positive relative density/pressure, with ordinary physical velocities.
    pub fn variables(
        &self,
        grid: &CartesianGrid,
        cell: usize,
        p: Primitive,
    ) -> Result<[f64; 4], Error> {
        let reference = self.cell(grid, cell).1;
        let rho = p.density / reference.density;
        let pressure = p.pressure / reference.pressure;
        if !positive(rho) || !positive(pressure) {
            return Err(Error::Invalid(
                "unrepresentable relative atmosphere perturbation",
            ));
        }
        Ok([rho, p.velocity[0], p.velocity[1], pressure])
    }

    pub fn trace(&self, grid: &CartesianGrid, cell: usize, axis: usize, sign: f64) -> [f64; 4] {
        let i = self.position(grid, cell);
        if axis == self.axis {
            self.faces[i + usize::from(sign > 0.0)].0.array()
        } else {
            // A transverse face spans the same complete stratified cell.
            self.cells[i].1.array()
        }
    }

    pub fn traction(&self, grid: &CartesianGrid, face: usize) -> f64 {
        let i = if self.axis == 0 {
            face % (grid.shape[0] + 1)
        } else {
            face / grid.shape[0]
        };
        self.faces[i].1
    }

    /// Perturbation force density, total force density, gravitational power density.
    pub fn sources(
        &self,
        grid: &CartesianGrid,
        cell: usize,
        state: State,
        x_flux: &[NormalFaceFlux],
        y_flux: &[NormalFaceFlux],
    ) -> [f64; 3] {
        let i = self.position(grid, cell);
        let g = acceleration(self.reference);
        let perturbation = g * (state.density - self.cells[i].0.density);
        let force =
            (self.faces[i + 1].1 - self.faces[i].1) / grid.spacing[self.axis] + perturbation;
        let (left, right) = if self.axis == 0 {
            let f = cell + cell / grid.shape[0];
            (x_flux[f].mass, x_flux[f + 1].mass)
        } else {
            (y_flux[cell].mass, y_flux[cell + grid.shape[0]].mass)
        };
        // -sum_f (phi_f - phi_cell)*mass_out / V, with phi=-g*(x-x_ref).
        // This exactly pairs the energy source with conservative mass transport.
        [perturbation, force, g * (0.5 * left + 0.5 * right)]
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

impl Euler2d {
    /// Initialize exact exponential cell averages and retain uniform gravity.
    ///
    /// `axis` is 0 for x or 1 for y, with at least two cells and physical ends.
    /// The other axis may be periodic. Reflecting ends remain impermeable;
    /// transmissive ends extrapolate perturbations to the physical hydrostatic
    /// face. Temperature describes the reference, not a thermostat.
    pub fn new_isothermal_atmosphere(
        gas: FrozenGas,
        grid: CartesianGrid,
        axis: usize,
        reference: IsothermalAtmosphere,
        boundaries: Boundaries2d,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        Self::new_hydrostatic_atmosphere(
            gas,
            grid,
            axis,
            HydrostaticReference::Isothermal(reference),
            AtmosphereQuadrature::default(),
            boundaries,
            cancelled,
        )
    }

    /// Initialize a source-integrated adiabatic atmosphere and retain gravity.
    ///
    /// The original gas enthalpy and entropy determine the temperature lapse
    /// and pressure; no frozen heat-capacity ratio is used. Conserved density
    /// and energy are integrated with the explicit preparation budget. The
    /// resulting axial profile is cached, so stepping performs no quadrature.
    /// Perturbations retain the full Euler energy equation, not fixed entropy.
    pub fn new_isentropic_atmosphere(
        gas: FrozenGas,
        grid: CartesianGrid,
        axis: usize,
        reference: IsentropicAtmosphere,
        quadrature: AtmosphereQuadrature,
        boundaries: Boundaries2d,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        Self::new_hydrostatic_atmosphere(
            gas,
            grid,
            axis,
            HydrostaticReference::Isentropic(reference),
            quadrature,
            boundaries,
            cancelled,
        )
    }

    /// Initialize a mechanically balanced prescribed thermal lapse.
    ///
    /// Stable and unstable profiles use the original source Cp(T), integrated
    /// conserved averages and shared gravity work. Only one O(n_axis) profile
    /// is prepared; stepping retains the full planar Euler dynamics.
    pub fn new_lapse_rate_atmosphere(
        gas: FrozenGas,
        grid: CartesianGrid,
        axis: usize,
        reference: LapseRateAtmosphere,
        quadrature: AtmosphereQuadrature,
        boundaries: Boundaries2d,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        Self::new_hydrostatic_atmosphere(
            gas,
            grid,
            axis,
            HydrostaticReference::LapseRate(reference),
            quadrature,
            boundaries,
            cancelled,
        )
    }

    fn new_hydrostatic_atmosphere(
        gas: FrozenGas,
        grid: CartesianGrid,
        axis: usize,
        reference: HydrostaticReference,
        quadrature: AtmosphereQuadrature,
        boundaries: Boundaries2d,
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
        let mut flow = Self::new(gas, grid, states, boundaries, &mut cancelled)?;
        if !balance
            .potential_energy(&grid, &flow.states, &mut cancelled)?
            .is_finite()
        {
            return Err(Error::Invalid(
                "unrepresentable gravitational potential energy",
            ));
        }
        poll(&mut cancelled)?;
        flow.hydrostatic = Some(balance);
        Ok(flow)
    }

    /// Bind the same gravity/reference to caller-authored perturbed averages.
    ///
    /// Only an initial field may acquire gravity; an existing reference cannot
    /// be replaced. States and their caloric energy gauge are never changed.
    /// Relative density and pressure must be representable in the reference.
    /// Refusal or cancellation leaves the owner unchanged.
    pub fn bind_isothermal_atmosphere(
        &mut self,
        axis: usize,
        reference: IsothermalAtmosphere,
        cancelled: impl FnMut() -> bool,
    ) -> Result<(), Error> {
        self.bind_hydrostatic_reference(
            axis,
            HydrostaticReference::Isothermal(reference),
            AtmosphereQuadrature::default(),
            cancelled,
        )
    }

    /// Bind the source-owned isentrope to caller-authored initial averages.
    /// The accepted state and its caloric gauge are unchanged. A different
    /// reference or preparation budget cannot replace an existing binding.
    pub fn bind_isentropic_atmosphere(
        &mut self,
        axis: usize,
        reference: IsentropicAtmosphere,
        quadrature: AtmosphereQuadrature,
        cancelled: impl FnMut() -> bool,
    ) -> Result<(), Error> {
        self.bind_hydrostatic_reference(
            axis,
            HydrostaticReference::Isentropic(reference),
            quadrature,
            cancelled,
        )
    }

    /// Retain a prescribed thermal lapse around caller-authored perturbations.
    /// This only binds an initial reference and never changes the gas averages.
    /// It shares the same source-domain, cancellation and single-binding rules
    /// as the isentropic and isothermal paths.
    pub fn bind_lapse_rate_atmosphere(
        &mut self,
        axis: usize,
        reference: LapseRateAtmosphere,
        quadrature: AtmosphereQuadrature,
        cancelled: impl FnMut() -> bool,
    ) -> Result<(), Error> {
        self.bind_hydrostatic_reference(
            axis,
            HydrostaticReference::LapseRate(reference),
            quadrature,
            cancelled,
        )
    }

    pub(super) fn bind_hydrostatic_reference(
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
        if self.reservoirs.any() {
            return Err(Error::Invalid(
                "finite tanks need declared hydrostatic port potentials",
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

    /// Retained isothermal reference; use [`Self::hydrostatic_reference`] for
    /// every supported physical family.
    pub fn atmosphere(&self) -> Option<(usize, IsothermalAtmosphere)> {
        match self.hydrostatic_reference() {
            Some((axis, HydrostaticReference::Isothermal(reference))) => Some((axis, reference)),
            _ => None,
        }
    }

    /// Retained gravity axis and original thermodynamic reference.
    pub fn hydrostatic_reference(&self) -> Option<(usize, HydrostaticReference)> {
        self.hydrostatic.as_ref().map(|b| (b.axis, b.reference))
    }

    /// `sum rho_i V_i phi(x_i)` [J], zero without a gravity reference.
    /// Adding this to gas energy gives the conserved coordinate in a closed
    /// box. Open boundaries also transport `phi_face * mass_flux`.
    /// The midpoint potential quadrature is second order for a disturbed field.
    pub fn gravitational_potential_energy(&self) -> f64 {
        self.hydrostatic.as_ref().map_or(0.0, |b| {
            b.potential_energy(&self.grid, &self.states, &mut || false)
                .expect("infallible cancellation callback")
        })
    }
}

#[cfg(test)]
mod tests;
