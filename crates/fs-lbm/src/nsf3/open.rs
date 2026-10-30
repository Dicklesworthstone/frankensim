//! Prescribed ghost states for the hybrid NSF lattice and energy equation.
//!
//! A ghost centre lies one cell spacing outside the adjacent fluid centre;
//! their common face is halfway between them. Incoming kinetic populations use
//! non-equilibrium extrapolation (Guo, Zheng & Shi, Chinese Physics 11 (2002),
//! 366–374, doi:10.1088/1009-1963/11/4/310). The thermal-pressure force and total
//! energy flux use the SAME ghost state. This is a boundary primitive, not an
//! acoustic non-reflection condition. In particular a prescribed state does not
//! preserve the outgoing characteristic of an impinging sound wave.

use super::{NsfError, NsfGrid3, THETA0, check_state};
use crate::lattice3::Lattice3;

/// Physical state at an exterior ghost centre, held until explicitly changed.
///
/// This specifies incoming populations, not an exact Dirichlet condition on
/// the adjacent fluid centre. Values must come from the caller's physical
/// boundary or a coupled solver; no pressure, temperature or velocity defaults.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NsfOpenState {
    /// Mass density [kg/m³].
    pub density_kg_m3: f64,
    /// Velocity in the grid's Cartesian frame [m/s].
    pub velocity_m_s: [f64; 3],
    /// Absolute temperature [K].
    pub temperature_k: f64,
}

/// A finite set of boundary-face tiles connected to one prescribed reservoir.
///
/// `cells` names the adjacent fluid cells, all on the selected domain face.
/// Each tile has area `dx^2`. The reservoir shares the grid's material law;
/// it does not introduce an independent gas composition or an acoustic load.
#[derive(Clone, Debug, PartialEq)]
pub struct NsfReservoirAperture {
    pub axis: usize,
    pub high: bool,
    pub cells: Vec<[usize; 3]>,
    pub state: NsfOpenState,
}

/// Original continuum pressure/flow projection on one localized reservoir tile.
/// Face pressure and velocity use the midpoint of the adjacent fluid and
/// prescribed ghost states, integrated with the original endpoint trapezoid.
/// The pressure-work fields contain the integral of that pointwise product,
/// not the product of interval means. They omit enthalpy, heat and viscous
/// stress work and never replace the original conservative total-energy flux.
/// Compressible streamed mass remains independent of continuum volume flow.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NsfReservoirFlow {
    /// Adjacent original fluid cell, in canonical z/y/x order.
    pub cell: [usize; 3],
    /// Interval-mean outward volume flow [m³/s], retaining its sign.
    pub outward_volume_flow_m3_s: f64,
    /// Interval-mean absolute face pressure [Pa].
    pub mean_pressure_pa: f64,
    /// Interval-mean face pressure relative to the grid's original source
    /// reference pressure [Pa], not to the mutable prescribed ghost.
    pub mean_acoustic_pressure_pa: f64,
    /// Integral of absolute face pressure times outward volume flow [J].
    pub pressure_work_out_j: f64,
    /// Integral of reference-relative pressure times outward volume flow [J].
    /// Positive work leaves the resolved fluid. Spatial/time covariance stays
    /// in this record rather than being converted into dissipation.
    pub acoustic_pressure_work_out_j: f64,
}
impl NsfReservoirFlow {
    pub(super) fn is_finite(&self) -> bool {
        [self.outward_volume_flow_m3_s, self.mean_pressure_pa,
            self.mean_acoustic_pressure_pa, self.pressure_work_out_j,
            self.acoustic_pressure_work_out_j].iter().all(|x| x.is_finite())
    }
}

#[derive(Clone)]
pub(super) struct BoundReservoirAperture {
    pub source: NsfReservoirAperture,
    indices: Vec<usize>,
    ghost: OpenGhost,
}

#[derive(Clone, Copy)]
pub(super) struct OpenGhost {
    pub rho: f64,
    pub u: [f64; 3],
    pub theta: f64,
    pub energy: f64,
}

impl<L: Lattice3<27>> NsfGrid3<L> {
    /// Open one localized reservoir aperture on a stationary transverse wall
    /// of an already attached characteristic axis, before the first step.
    ///
    /// Tiles must be unique fluid cells, away from other domain edges and the
    /// two-cell characteristic collars. Their one-cell interior/tangential
    /// halo must also be fluid. Thus a lattice link never intersects this
    /// aperture and a characteristic face or internal solid simultaneously.
    /// A diagonal crossing exactly on an aperture/wall tile edge belongs to
    /// the wall unless every incident tile is open. Streaming uses this rule
    /// for both directions of the same link; no axis-order corner choice is
    /// introduced. Normal finite-volume fluxes retain their original tile area.
    ///
    /// Incoming populations, thermal-pressure force and enthalpy/conduction
    /// use the original prescribed-ghost law. Its signed energy is reported in
    /// `EnergyLedger::open_energy_j`, separately from characteristic exchange.
    /// The full Guo force chart is rebased without changing physical cell
    /// density, velocity or energy. Only one aperture is admitted, with at most
    /// the grid's declared cell count; no dense extra grid is allocated.
    ///
    /// This is a prescribed reservoir, not a transparent acoustic boundary,
    /// finite gas inventory, qualified aperture-edge scheme or flute model.
    pub fn with_reservoir_aperture(
        mut self,
        source: NsfReservoirAperture,
    ) -> Result<Self, NsfError> {
        let open = self.open.ok_or(NsfError::Input(
            "a localized reservoir aperture requires an attached characteristic axis",
        ))?;
        if self.steps != 0 || self.reservoir_aperture.is_some() || self.sources.is_some()
            || self.moving_walls
            || self.reservoir_faces.iter().flatten().any(Option::is_some)
        {
            return Err(NsfError::Input(
                "one reservoir aperture requires a cold grid with stationary boundaries",
            ));
        }
        if source.axis >= 3 || source.axis == open.axis || self.spec.periodic[source.axis]
            || self.spec.dims[source.axis] < 2
            || source.cells.is_empty() || source.cells.len() > self.n
        {
            return Err(NsfError::Input(
                "reservoir aperture needs bounded nonempty tiles on a transverse non-periodic face",
            ));
        }
        let plane = if source.high { self.spec.dims[source.axis] - 1 } else { 0 };
        let mut indices = Vec::with_capacity(source.cells.len());
        for &x in &source.cells {
            if x.iter().zip(self.spec.dims).any(|(&v, n)| v >= n)
                || x[source.axis] != plane
                || x[open.axis] < 2 || x[open.axis] >= self.spec.dims[open.axis] - 2
                || (0..3).any(|a| a != source.axis && (x[a] == 0 || x[a] + 1 >= self.spec.dims[a]))
            {
                return Err(NsfError::Input(
                    "reservoir aperture tiles must avoid domain corners and characteristic collars",
                ));
            }
            // The ghost extrapolation has no solid/open junction rule. Admit
            // the local two-layer fluid collar explicitly instead of allowing
            // a diagonal population to bypass an interior obstacle.
            for da in -1..=1 {
                for db in -1..=1 {
                    for depth in 0..=1 {
                        let mut y = x;
                        let a = (source.axis + 1) % 3;
                        let b = (source.axis + 2) % 3;
                        y[a] = (y[a] as isize + da) as usize;
                        y[b] = (y[b] as isize + db) as usize;
                        y[source.axis] = if source.high { plane - depth } else { depth };
                        if self.solid[self.index(y)] {
                            return Err(NsfError::Input(
                                "reservoir aperture requires a fluid interior and tangential collar",
                            ));
                        }
                    }
                }
            }
            indices.push(self.index(x));
        }
        indices.sort_unstable();
        if indices.windows(2).any(|pair| pair[0] == pair[1]) {
            return Err(NsfError::Input("reservoir aperture tiles must be unique"));
        }
        let ghost = self.admit_open_state(source.state)?;
        self.reservoir_aperture = Some(BoundReservoirAperture { source, indices, ghost });
        self.last_reservoir_mass_out_kg = Some(0.0);
        self.last_reservoir_flows = self.reservoir_aperture.as_ref().expect("bound aperture")
            .indices.iter().map(|&index| NsfReservoirFlow {
                cell: self.at(index), outward_volume_flow_m3_s: 0.0,
                mean_pressure_pa: 0.0, mean_acoustic_pressure_pa: 0.0,
                pressure_work_out_j: 0.0, acoustic_pressure_work_out_j: 0.0,
            }).collect();
        let drives = self.boundary_drives().expect("admitted characteristic axis");
        self.set_boundary_drives(drives)?;
        Ok(self)
    }

    /// Original aperture geometry and physical reservoir state, if attached.
    #[must_use]
    pub fn reservoir_aperture(&self) -> Option<&NsfReservoirAperture> {
        self.reservoir_aperture.as_ref().map(|bound| &bound.source)
    }

    /// Resolve an attached reservoir at absolute pressure [Pa] and temperature
    /// [K] using this grid's original gas source. Its declared Cartesian ghost
    /// velocity is retained; pressure is not converted into a prescribed jet
    /// speed. Composition, caloric and transport support remain source-owned.
    /// This query changes no boundary or fluid state.
    pub fn reservoir_aperture_state_at(
        &self, pressure_pa: f64, temperature_k: f64,
    ) -> Result<NsfOpenState, NsfError> {
        let aperture = self.reservoir_aperture.as_ref().ok_or(NsfError::Input(
            "reservoir conditions require an attached aperture",
        ))?;
        let gas = self.spec.gas.at(temperature_k, pressure_pa)?;
        let state = NsfOpenState {
            density_kg_m3: gas.density,
            velocity_m_s: aperture.source.state.velocity_m_s,
            temperature_k: gas.temperature,
        };
        self.admit_open_state(state)?;
        Ok(state)
    }

    /// Change an existing aperture's held ghost state without changing its
    /// geometry or the fluid's physical density, velocity, stress or energy.
    /// The original force chart and energy-flux history are rebound together.
    /// Material/boundary refusal leaves every accepted field and control intact.
    /// Topology remains cold-only; an accepted warm grid may change conditions.
    pub fn set_reservoir_aperture_state(&mut self, state: NsfOpenState) -> Result<(), NsfError> {
        let aperture = self.reservoir_aperture.as_ref().ok_or(NsfError::Input(
            "reservoir state control requires an attached aperture",
        ))?;
        if aperture.source.state == state { return Ok(()); }
        let previous = (aperture.source.state, aperture.ghost);
        let ghost = self.admit_open_state(state)?;
        let drives = self.boundary_drives().ok_or(NsfError::Input(
            "reservoir state control requires its original characteristic axis",
        ))?;
        let aperture = self.reservoir_aperture.as_mut().expect("checked aperture");
        aperture.source.state = state;
        aperture.ghost = ghost;
        // The existing rebind stages its populations, force and fluxes before
        // publication. Restore only our staged ghost if that owner refuses.
        if let Err(error) = self.set_boundary_drives(drives) {
            let aperture = self.reservoir_aperture.as_mut().expect("retained aperture");
            aperture.source.state = previous.0;
            aperture.ghost = previous.1;
            return Err(error);
        }
        Ok(())
    }

    /// Actual streamed mass leaving the localized reservoir aperture in the
    /// last accepted step [kg]. Positive is outward; `None` means no aperture.
    /// Zero before the first step. Counted from original population links,
    /// independently of whole-grid storage and the characteristic-face accounts.
    #[must_use]
    pub const fn last_reservoir_mass_out_kg(&self) -> Option<f64> {
        self.last_reservoir_mass_out_kg
    }

    /// Immutable source-reference pressure used by all reservoir acoustic
    /// pressure and work accounts [Pa]. Warm ghost controls do not change it.
    #[must_use]
    pub const fn reservoir_reference_pressure_pa(&self) -> f64 {
        self.spec.pressure_pa
    }

    /// Accepted per-tile continuum pressure/flow accounts in canonical grid
    /// order. Cold records are zero because no interval has been integrated;
    /// an absent aperture gives an empty slice. Failed steps, warm control
    /// rebinding and dropped trials preserve this complete accepted record.
    #[must_use]
    pub fn last_reservoir_flows(&self) -> &[NsfReservoirFlow] {
        &self.last_reservoir_flows
    }

    pub(super) fn reservoir_volume_flows(
        &self, rho0: &[f64], before: &[[f64; 3]], theta0: &[f64],
        rho1: &[f64], after: &[[f64; 3]], theta1: &[f64],
    ) -> Vec<NsfReservoirFlow> {
        let Some(aperture) = &self.reservoir_aperture else { return Vec::new(); };
        let axis = aperture.source.axis;
        let sign = if aperture.source.high { 1.0 } else { -1.0 };
        let unit = self.law.dx_m.powi(3) / self.law.dt_s;
        let ghost_pressure = self.law.pressure_pa(aperture.ghost.rho, aperture.ghost.theta);
        let reference = self.spec.pressure_pa;
        aperture.indices.iter().map(|&index| {
            let pressure0 = 0.5 * (self.law.pressure_pa(rho0[index], theta0[index]) + ghost_pressure);
            let pressure1 = 0.5 * (self.law.pressure_pa(rho1[index], theta1[index]) + ghost_pressure);
            let acoustic0 = pressure0 - reference;
            let acoustic1 = pressure1 - reference;
            let q0 = sign * unit * 0.5 * (before[index][axis] + aperture.ghost.u[axis]);
            let q1 = sign * unit * 0.5 * (after[index][axis] + aperture.ghost.u[axis]);
            NsfReservoirFlow {
                cell: self.at(index),
                // Preserve the original flow arithmetic and observer replay.
                outward_volume_flow_m3_s: sign * unit * 0.25
                    * (before[index][axis] + after[index][axis] + 2.0 * aperture.ghost.u[axis]),
                mean_pressure_pa: 0.5 * (pressure0 + pressure1),
                mean_acoustic_pressure_pa: 0.5 * (acoustic0 + acoustic1),
                pressure_work_out_j: 0.5 * self.law.dt_s * (pressure0 * q0 + pressure1 * q1),
                acoustic_pressure_work_out_j: 0.5 * self.law.dt_s * (acoustic0 * q0 + acoustic1 * q1),
            }
        }).collect()
    }

    pub(super) fn reservoir_link_mass_out_kg(
        &self, post: &[[f64; 27]], streamed: &[[f64; 27]],
    ) -> Option<f64> {
        let aperture = self.reservoir_aperture.as_ref()?;
        let axis = aperture.source.axis;
        let side = usize::from(aperture.source.high);
        let normal = if aperture.source.high { 1 } else { -1 };
        let mut out = 0.0;
        // Sparse indices are sorted at admission: tile input order and worker
        // scheduling cannot change the reduction. Use the same half-link
        // ownership as streaming, and its actual incoming population value.
        for &i in &aperture.indices {
            let x = self.at(i);
            for (q, c) in L::E.iter().copied().enumerate() {
                if c[axis] == normal && self.reservoir_link_is_open(x, c, 1, axis, side) {
                    out += post[i][q] - streamed[i][L::OPP[q]];
                }
            }
        }
        Some(out * self.law.rho_ref_kg_m3 * self.law.dx_m.powi(3))
    }

    /// Replace one non-periodic face with a prescribed exterior state.
    /// `None` restores its stationary wall. Solid cells remain solid.
    ///
    /// This first boundary rung admits one open axis with periodic transverse
    /// domain edges. Internal solids may describe channels/nozzles, but may not
    /// intersect the open plane: solid/open diagonal-link corners are unqualified.
    /// It deliberately refuses intersecting open/wall domain edges rather than
    /// choosing a corner rule by axis order. It does not claim transparency,
    /// characteristic coupling, jet validation or a converged open-boundary order.
    ///
    /// Updating the ghost state restarts force/energy-flux history while
    /// preserving the current fluid's density, physical velocity and energy.
    /// A refusal leaves the boundary and all published state unchanged.
    ///
    /// # Errors
    /// Invalid axis/topology, a moving face, manufactured sources, or an
    /// exterior state outside the same material/Mach/temperature bounds as the grid.
    pub fn set_open_face(
        &mut self,
        axis: usize,
        high: bool,
        state: Option<NsfOpenState>,
    ) -> Result<(), NsfError> {
        if self.sgs.is_some() {
            return Err(NsfError::Input(
                "NSF SGS does not support prescribed-reservoir faces",
            ));
        }
        if self.open.is_some() {
            return Err(NsfError::Input(
                "prescribed reservoir and characteristic faces are exclusive",
            ));
        }
        if axis >= 3 || self.spec.periodic[axis] || self.spec.dims[axis] < 2 {
            return Err(NsfError::Input(
                "open face requires a non-periodic axis with at least two cells",
            ));
        }
        if (0..3).any(|a| a != axis && !self.spec.periodic[a]) {
            return Err(NsfError::Input(
                "open-face transverse domain edges must be periodic; use interior solids for walls",
            ));
        }
        let plane = if high { self.spec.dims[axis] - 1 } else { 0 };
        if state.is_some()
            && self
                .solid
                .iter()
                .enumerate()
                .any(|(i, solid)| *solid && self.at(i)[axis] == plane)
        {
            return Err(NsfError::Input(
                "open plane must be free of solids; solid/open junctions are unqualified",
            ));
        }
        let side = usize::from(high);
        if self.edge_velocity[axis][side] != [0.0; 3] || self.sources.is_some() {
            return Err(NsfError::Input(
                "open face cannot combine with a moving face or manufactured sources",
            ));
        }
        let ghost = state.map(|s| self.admit_open_state(s)).transpose()?;
        let previous = self.reservoir_faces[axis][side];
        let (rho, mom) = self.macro_fields(&self.f);
        let u = self.velocity_field(&rho, &mom, &self.force);
        self.reservoir_faces[axis][side] = ghost;
        let force = self.force_field(&rho, &self.theta, None);
        let mut populations = self.f.clone();
        for (i, f) in populations.iter_mut().enumerate() {
            if self.solid[i] {
                continue;
            }
            // Preserve the whole physical distribution, including its stress,
            // when changing the thermal-pressure force of a reservoir face.
            super::rebase_force::<L>(
                f,
                u[i],
                std::array::from_fn(|a| self.force[i][a] - force[i][a]),
            );
        }
        let sig = self.stress_field(&populations, &force, &rho, &self.theta, &self.material, &[]);
        let history = self.flux_divergence(&self.e, &rho, &u, &self.theta, &sig, &[], None);
        let (div, wall, open) = match history {
            Ok(value) => value,
            Err(error) => {
                self.reservoir_faces[axis][side] = previous;
                return Err(error);
            }
        };
        self.f = populations;
        self.force = force;
        self.div = div;
        self.wall_out = wall;
        self.open_out = open;
        Ok(())
    }

    fn admit_open_state(&self, s: NsfOpenState) -> Result<OpenGhost, NsfError> {
        let rho = s.density_kg_m3 / self.law.rho_ref_kg_m3;
        let u = s.velocity_m_s.map(|v| v * self.law.dt_s / self.law.dx_m);
        let theta = self.law.theta_of_kelvin(s.temperature_k);
        if !(rho.is_finite()
            && rho > 0.0
            && u.iter().all(|v| v.is_finite())
            && theta.is_finite()
            && theta > 0.0)
        {
            return Err(NsfError::Input(
                "open-face density, velocity and temperature must be finite and physical",
            ));
        }
        let [low_k, high_k] = self.spec.temperature_window_k;
        if s.temperature_k < low_k || s.temperature_k > high_k {
            return Err(NsfError::Input(
                "open-face temperature outside the declared run window",
            ));
        }
        let material = self.spec.gas.resolve(&self.law, rho, theta)?;
        let gas = self
            .spec
            .gas
            .at(s.temperature_k, self.law.pressure_pa(rho, theta))?;
        check_state(&self.spec, gas.gamma, self.n, rho, u, theta)?;
        let internal = material
            .as_ref()
            .map_or(theta / (self.law.gamma - 1.0), |m| {
                m.sensible_energy_lattice
            });
        let energy = rho * (internal + 0.5 * u.iter().map(|v| v * v).sum::<f64>());
        if !energy.is_finite() {
            return Err(NsfError::Input("open-face energy is unrepresentable"));
        }
        Ok(OpenGhost {
            rho,
            u,
            theta,
            energy,
        })
    }

    pub(super) fn open_at(&self, i: usize, axis: usize, direction: i64) -> Option<OpenGhost> {
        if self.shifted(self.at(i), axis, direction).is_none() {
            self.reservoir_ghost(i, axis, usize::from(direction > 0))
        } else {
            None
        }
    }

    pub(super) fn reservoir_ghost(&self, i: usize, axis: usize, side: usize) -> Option<OpenGhost> {
        if let Some(ghost) = self.reservoir_faces[axis][side] {
            return Some(ghost);
        }
        let aperture = self.reservoir_aperture.as_ref()?;
        (aperture.source.axis == axis && usize::from(aperture.source.high) == side
            && aperture.indices.binary_search(&i).is_ok()).then_some(aperture.ghost)
    }

    /// Wall precedence at the geometric half-link intersection. A tangential
    /// half-cell coordinate touches two tiles (four for a corner diagonal);
    /// each must be part of the same aperture. Reversing the link and its
    /// direction leaves this intersection and its ownership unchanged.
    pub(super) fn reservoir_link_is_open(
        &self, x: [usize; 3], c: [i32; 3], direction: i32, axis: usize, side: usize,
    ) -> bool {
        if self.reservoir_faces[axis][side].is_some() {
            return true;
        }
        let Some(aperture) = self.reservoir_aperture.as_ref() else { return false; };
        if aperture.source.axis != axis || usize::from(aperture.source.high) != side {
            return false;
        }
        let axes = [(axis + 1) % 3, (axis + 2) % 3];
        for mask in 0..4 {
            let mut y = x;
            for (bit, &a) in axes.iter().enumerate() {
                if mask & (1 << bit) != 0 {
                    let next = x[a] as i64 + i64::from(direction) * i64::from(c[a]);
                    if next < 0 || next >= self.spec.dims[a] as i64 { return false; }
                    y[a] = next as usize;
                }
            }
            if aperture.indices.binary_search(&self.index(y)).is_err() { return false; }
        }
        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::d3q27::D3q27;
    use crate::nsf3::{MAX_QUALIFIED_MACH, NsfCharacteristicAxis, NsfCharacteristicBoundary,
        NsfSpec, WallThermal};
    use fs_material::gas::{GasSpec, GasState};

    fn spec(axis: usize, workers: usize) -> NsfSpec {
        let mut dims = [3; 3];
        dims[axis] = 24;
        let mut periodic = [true; 3];
        periodic[axis] = false;
        NsfSpec {
            dims,
            periodic,
            gas: GasSpec::dry_air_ussa1976(),
            temperature_k: 293.15,
            pressure_pa: 101325.0,
            dx_m: 5e-5,
            temperature_window_k: [273.15, 333.15],
            maximum_mach: MAX_QUALIFIED_MACH,
            maximum_temperature_ratio: 0.2,
            bulk_to_shear_ratio: 0.73,
            bulk_ratio_source: "Cramer 2012, Phys. Fluids 24 066102, N2 at 293 K",
            maximum_cells: 216,
            wall: WallThermal::Adiabatic,
            workers,
            transport_points: 33,
        }
    }

    fn state(velocity_m_s: [f64; 3], temperature_k: f64) -> NsfOpenState {
        let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), temperature_k, 101325.0).unwrap();
        NsfOpenState {
            density_kg_m3: gas.density,
            velocity_m_s,
            temperature_k,
        }
    }

    fn uniform(axis: usize, workers: usize, s: NsfOpenState) -> NsfGrid3<D3q27> {
        let spec = spec(axis, workers);
        let reference = GasState::try_new(&spec.gas, spec.temperature_k, spec.pressure_pa).unwrap();
        let speed_unit = (3.0 * reference.specific_gas_constant * reference.temperature).sqrt();
        NsfGrid3::new(
            spec,
            |_| false,
            |_| {
                (
                    s.density_kg_m3 / reference.density,
                    s.velocity_m_s.map(|u| u / speed_unit),
                    s.temperature_k,
                )
            },
        )
        .unwrap()
    }

    fn aperture_grid(axis: usize, workers: usize) -> NsfGrid3<D3q27> {
        let mut s = spec(axis, workers);
        s.dims = [9; 3];
        s.maximum_cells = 9usize.pow(3);
        s.periodic = [false; 3];
        NsfGrid3::new(s, |_| false, |_| (1.0, [0.0; 3], 293.15))
            .unwrap()
            .with_characteristic_axis(NsfCharacteristicAxis {
                axis,
                faces: [NsfCharacteristicBoundary {
                    reference_pressure_pa: 101325.0,
                    reference_temperature_k: 293.15,
                    reference_velocity_m_s: [0.0; 3],
                    incoming_characteristic_m_s: 0.0,
                }; 2],
                maximum_relative_perturbation: 0.05,
            })
            .unwrap()
    }

    fn aperture(axis: usize, high: bool, temperature_k: f64) -> NsfReservoirAperture {
        let mut cells = Vec::new();
        for a in 3..=5 {
            for b in 3..=5 {
                let mut x = [0; 3];
                x[axis] = if high { 8 } else { 0 };
                x[(axis + 1) % 3] = a;
                x[(axis + 2) % 3] = b;
                cells.push(x);
            }
        }
        NsfReservoirAperture { axis, high, cells, state: state([0.0; 3], temperature_k) }
    }

    #[test]
    fn g0_reservoir_aperture_preserves_equilibrium_all_transverse_faces() {
        for characteristic_axis in 0..3 {
            for axis in 0..3 {
                if axis == characteristic_axis { continue; }
                for high in [false, true] {
                    let source = aperture(axis, high, 293.15);
                    let mut grid = aperture_grid(characteristic_axis, 1)
                        .with_reservoir_aperture(source.clone()).unwrap();
                    assert_eq!(grid.reservoir_aperture(), Some(&source));
                    let before = grid.total_energy_j();
                    for _ in 0..4 {
                        let ledger = grid.step().unwrap();
                        assert!(ledger.residual_rel < 2e-12, "{ledger:?}");
                    }
                    for i in 0..grid.n {
                        let (rho, u, temperature) = grid.cell(grid.at(i));
                        assert!((rho - 1.0).abs() < 2e-13, "axis={axis} high={high} rho={rho}");
                        assert!(u.iter().all(|v| v.abs() < 2e-13));
                        assert!((temperature / 293.15 - 1.0).abs() < 2e-13);
                    }
                    assert!((grid.total_energy_j() / before - 1.0).abs() < 2e-12);
                }
            }
        }
    }

    #[test]
    fn g0_reservoir_aperture_diagonal_edges_keep_wall_ownership() {
        for axis in 0..3 {
            for high in [false, true] {
                let source = aperture(axis, high, 293.15);
                let grid = aperture_grid((axis + 1) % 3, 1)
                    .with_reservoir_aperture(source).unwrap();
                let mut centre = [4; 3];
                centre[axis] = if high { 8 } else { 0 };
                let mut edge = centre;
                edge[(axis + 1) % 3] = 3;
                let mut corner = edge;
                corner[(axis + 2) % 3] = 3;
                let side = usize::from(high);
                for x in [centre, edge, corner] {
                    for c in D3q27::E {
                        if c[axis] != if high { -1 } else { 1 } { continue; }
                        let expected = (0..3).filter(|&a| a != axis)
                            .all(|a| (3..=5).contains(&(x[a] as i32 - c[a])));
                        let owned = grid.reservoir_link_is_open(x, c, -1, axis, side);
                        assert_eq!(owned, expected, "axis={axis} high={high} x={x:?} c={c:?}");
                        assert_eq!(owned, grid.reservoir_link_is_open(
                            x, c.map(|v| -v), 1, axis, side));
                    }
                }
                // Exercise the real streamer, with distinguishable opposite
                // populations. A diagonal crossing the aperture edge must
                // bounce back; a wholly open diagonal uses the original ghost.
                // A third-order xyz mode has zero mass and momentum but
                // opposite signs on opposite corner populations.
                let post: Vec<_> = (0..grid.n).map(|_| std::array::from_fn(|q| {
                    let [x, y, z] = D3q27::E[q];
                    D3q27::W[q] * (1.0 + 1e-3 * f64::from(x * y * z))
                })).collect();
                let (rho, mom) = grid.macro_fields(&grid.f);
                let u = grid.velocity_field(&rho, &mom, &grid.force);
                let incoming = grid.open_incoming_populations(&post, &rho, &u);
                let streamed = grid.stream(&post, incoming.as_ref());
                let i = grid.index(edge);
                for (q, c) in D3q27::E.iter().copied().enumerate() {
                    if c[axis] != if high { -1 } else { 1 } { continue; }
                    if grid.reservoir_link_is_open(edge, c, -1, axis, side) {
                        // At a matched rest ghost, the extrapolated equilibrium
                        // difference vanishes and the local post population stays.
                        assert!((streamed[i][q] - post[i][q]).abs() < 1e-15);
                    } else {
                        assert_eq!(streamed[i][q], post[i][D3q27::OPP[q]]);
                    }
                }
            }
        }
    }

    #[test]
    fn g0_reservoir_aperture_heat_uses_original_flux_and_separate_ledger() {
        let mut grid = aperture_grid(0, 1);
        let cells: Vec<_> = (0..grid.n).map(|i| grid.cell(grid.at(i))).collect();
        let before = grid.total_energy_j();
        let mass_before = grid.total_mass_kg();
        grid = grid.with_reservoir_aperture(aperture(1, false, 293.35)).unwrap();
        assert_eq!(grid.total_energy_j(), before);
        for (i, expected) in cells.iter().enumerate() {
            let actual = grid.cell(grid.at(i));
            assert!((actual.0 - expected.0).abs() < 2e-14);
            assert_eq!(actual.2, expected.2);
            for a in 0..3 { assert!((actual.1[a] - expected.1[a]).abs() < 2e-14); }
        }
        let mut mass_out = 0.0;
        let mut reservoir_mass_out = 0.0;
        for _ in 0..12 {
            let ledger = grid.step().unwrap();
            assert_eq!(ledger.wall_heat_j, 0.0);
            assert_eq!(ledger.viscous_heating_j, 0.0);
            assert!(ledger.residual_rel < 2e-12, "{ledger:?}");
            let aperture_mass = grid.last_reservoir_mass_out_kg().unwrap();
            reservoir_mass_out += aperture_mass;
            mass_out += aperture_mass + grid.last_boundary_exchange().unwrap()
                .iter().map(|face| face.mass_out_kg).sum::<f64>();
        }
        let ledger = grid.cumulative_ledger();
        assert!(ledger.open_energy_j > 0.0, "{ledger:?}");
        assert!(grid.cell([4, 0, 4]).2 > 293.15);
        let error = grid.total_energy_j() - before
            - ledger.open_energy_j - ledger.open_boundary_energy_j - ledger.wall_heat_j;
        assert!(error.abs() / before < 2e-12, "{ledger:?} error={error:e}");
        let mass_error = grid.total_mass_kg() - mass_before + mass_out;
        assert!(mass_error.abs() / mass_before < 2e-12, "mass_error={mass_error:e}");
        assert!(reservoir_mass_out.abs() > 1e-12 * mass_before);
        let mut sealed = aperture_grid(0, 1);
        for _ in 0..12 { sealed.step().unwrap(); }
        assert!(grid.cell([4, 0, 4]).2 > sealed.cell([4, 0, 4]).2 + 1e-4);
    }

    #[test]
    fn g0_reservoir_aperture_trials_cancel_repeat_and_accept_exactly() {
        let source = aperture(1, false, 293.35);
        let mut grid = aperture_grid(0, 1).with_reservoir_aperture(source.clone()).unwrap();
        let before = (grid.state_hash(), grid.steps(), grid.cumulative_ledger(),
            grid.last_reservoir_mass_out_kg());
        let drives = grid.boundary_drives().unwrap();
        let mut candidate_hash = 0;
        let mut candidate_ledger = None;
        let mut candidate_mass = None;
        for _ in 0..2 {
            let trial = grid.preview_boundary_interval(drives, 3, 3 * 729, || false).unwrap();
            assert_eq!(trial.grid().reservoir_aperture(), Some(&source));
            if let Some(expected) = candidate_ledger {
                assert_eq!(*trial.energy_ledger(), expected);
                assert_eq!(trial.grid().state_hash(), candidate_hash);
                assert_eq!(Some(trial.reservoir_mass_out_kg()), candidate_mass);
            }
            candidate_hash = trial.grid().state_hash();
            candidate_ledger = Some(*trial.energy_ledger());
            candidate_mass = Some(trial.reservoir_mass_out_kg());
            drop(trial);
            assert_eq!((grid.state_hash(), grid.steps(), grid.cumulative_ledger(),
                grid.last_reservoir_mass_out_kg()), before);
            assert_eq!(grid.reservoir_aperture(), Some(&source));
        }
        let calls = std::cell::Cell::new(0);
        assert!(grid.preview_boundary_interval(drives, 3, 3 * 729, || {
            calls.set(calls.get() + 1);
            calls.get() == 3
        }).is_err());
        assert_eq!((grid.state_hash(), grid.steps(), grid.cumulative_ledger(),
            grid.last_reservoir_mass_out_kg()), before);
        grid.preview_boundary_interval(drives, 3, 3 * 729, || false).unwrap().accept();
        assert_eq!(grid.state_hash(), candidate_hash);
        assert_eq!(grid.steps(), 3);
        assert_eq!(grid.reservoir_aperture(), Some(&source));
        assert!(candidate_ledger.unwrap().open_energy_j > 0.0);
    }

    #[test]
    fn g0_warm_reservoir_conditions_preserve_fluid_and_refuse_invalid_sources() {
        let mut grid = aperture_grid(0, 1)
            .with_reservoir_aperture(aperture(1, false, 293.35)).unwrap();
        for _ in 0..3 { grid.step().unwrap(); }
        let cells: Vec<_> = (0..grid.n).map(|i| grid.cell(grid.at(i))).collect();
        let stored = (grid.total_energy_j(), grid.steps(), grid.cumulative_ledger(),
            grid.last_reservoir_mass_out_kg());
        let hash = grid.state_hash();
        let state = grid.reservoir_aperture_state_at(101345.0, 293.45).unwrap();
        let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.45, 101345.0).unwrap();
        assert_eq!(state.density_kg_m3, gas.density);
        assert_eq!(state.velocity_m_s, [0.0; 3]);
        assert_eq!(grid.state_hash(), hash, "resolving conditions must be read-only");
        grid.set_reservoir_aperture_state(state).unwrap();
        assert_eq!(grid.reservoir_aperture().unwrap().state, state);
        assert_ne!(grid.state_hash(), hash, "accepted reservoir conditions are part of replay");
        assert_eq!((grid.total_energy_j(), grid.steps(), grid.cumulative_ledger(),
            grid.last_reservoir_mass_out_kg()), stored);
        for (i, expected) in cells.iter().enumerate() {
            let actual = grid.cell(grid.at(i));
            assert!((actual.0 - expected.0).abs() < 2e-14);
            assert_eq!(actual.2, expected.2);
            for a in 0..3 { assert!((actual.1[a] - expected.1[a]).abs() < 2e-14); }
        }
        let before = (grid.state_hash(), grid.reservoir_aperture().cloned(),
            grid.cumulative_ledger(), grid.boundary_drives());
        for (pressure, temperature) in [(0.0, 293.15), (101325.0, 400.0),
            (f64::NAN, 293.15), (101325.0, f64::NAN)]
        {
            assert!(grid.reservoir_aperture_state_at(pressure, temperature).is_err());
        }
        for invalid in [NsfOpenState { density_kg_m3: 0.0, ..state },
            NsfOpenState { velocity_m_s: [1000.0, 0.0, 0.0], ..state },
            NsfOpenState { temperature_k: 400.0, ..state }]
        {
            assert!(grid.set_reservoir_aperture_state(invalid).is_err());
            assert_eq!((grid.state_hash(), grid.reservoir_aperture().cloned(),
                grid.cumulative_ledger(), grid.boundary_drives()), before);
        }
        let mut closed = aperture_grid(0, 1);
        assert!(closed.reservoir_aperture_state_at(101325.0, 293.15).is_err());
        assert!(closed.set_reservoir_aperture_state(state).is_err());
    }

    #[test]
    fn g3_reservoir_control_trials_replay_mass_energy_and_roll_back_conditions() {
        reservoir_control_replay(false);
        reservoir_control_replay(true);
    }

    fn reservoir_control_replay(modeled_transport: bool) {
        let build = || {
            let grid = aperture_grid(0, 1)
                .with_reservoir_aperture(aperture(1, false, 293.15)).unwrap();
            if modeled_transport {
                grid.with_sgs(crate::nsf3::NsfSgsConfig {
                    closure: crate::sgs3::SgsConfig {
                        model: crate::sgs3::SgsModel::Smagorinsky { cs: 0.17 },
                        filter_width_lu: 1.0,
                        maximum_eddy_to_molecular: 100.0,
                        usage: crate::sgs3::SgsUse::Exploration,
                    },
                    turbulent_prandtl: 0.9,
                }).unwrap()
            } else { grid }
        };
        let mut grid = build();
        let mut original = build();
        for _ in 0..2 { grid.step().unwrap(); original.step().unwrap(); }
        let before = (grid.state_hash(), grid.steps(), grid.reservoir_aperture().cloned(),
            grid.cumulative_ledger(), grid.last_reservoir_mass_out_kg(), grid.sgs_transport.clone());
        let state = grid.reservoir_aperture_state_at(101345.0, 293.35).unwrap();
        let drives = grid.boundary_drives().unwrap();
        let first = grid.preview_boundary_interval_with_reservoir(
            drives, None, state, 3, 3 * 729, || false).unwrap();
        let expected = (first.grid().state_hash(), *first.energy_ledger(), first.reservoir_mass_out_kg());
        let mass_error = first.grid().total_mass_kg() - first.previous_mass_kg()
            + first.reservoir_mass_out_kg()
            + first.boundary_exchange().iter().map(|face| face.mass_out_kg).sum::<f64>();
        assert!(mass_error.abs() / first.previous_mass_kg() < 2e-12);
        assert!(first.energy_ledger().residual_rel < 2e-12);
        assert!(first.energy_ledger().open_energy_j.abs() > 0.0);
        drop(first);
        assert_eq!((grid.state_hash(), grid.steps(), grid.reservoir_aperture().cloned(),
            grid.cumulative_ledger(), grid.last_reservoir_mass_out_kg(), grid.sgs_transport.clone()), before);
        let calls = std::cell::Cell::new(0);
        assert!(grid.preview_boundary_interval_with_reservoir(
            drives, None, state, 3, 3 * 729, || {
                calls.set(calls.get() + 1); calls.get() == 3
            }).is_err());
        assert_eq!((grid.state_hash(), grid.steps(), grid.reservoir_aperture().cloned(),
            grid.cumulative_ledger(), grid.last_reservoir_mass_out_kg(), grid.sgs_transport.clone()), before);
        let retry = grid.preview_boundary_interval_with_reservoir(
            drives, None, state, 3, 3 * 729, || false).unwrap();
        assert_eq!((retry.grid().state_hash(), *retry.energy_ledger(), retry.reservoir_mass_out_kg()), expected);
        retry.accept();
        original.set_reservoir_aperture_state(state).unwrap();
        original.set_boundary_drives(drives).unwrap();
        for _ in 0..3 { original.step().unwrap(); }
        assert_eq!(grid.state_hash(), original.state_hash());
        assert_eq!(grid.cumulative_ledger(), original.cumulative_ledger());
        assert_eq!(grid.reservoir_aperture(), original.reservoir_aperture());
        assert_eq!(grid.sgs_transport, original.sgs_transport);
        if modeled_transport {
            assert!(grid.sgs_transport.iter().any(|t| t.eddy_viscosity_lu > 0.0));
        }
    }

    #[test]
    fn g0_reservoir_aperture_refuses_ambiguous_topology_and_warm_attachment() {
        let valid = aperture(1, false, 293.15);
        let mut bad_sources = vec![aperture(0, false, 293.15)];
        let mut duplicate = valid.clone();
        duplicate.cells.push(duplicate.cells[0]);
        bad_sources.push(duplicate);
        for cell in [[1, 0, 4], [4, 0, 0], [4, 1, 4], [9, 0, 4]] {
            let mut invalid = valid.clone(); invalid.cells = vec![cell]; bad_sources.push(invalid);
        }
        let mut empty = valid.clone(); empty.cells.clear(); bad_sources.push(empty);
        let mut hot = valid.clone(); hot.state.temperature_k = 400.0; bad_sources.push(hot);
        for source in bad_sources {
            assert!(aperture_grid(0, 1).with_reservoir_aperture(source).is_err());
        }
        assert!(uniform(0, 1, state([0.0; 3], 293.15))
            .with_reservoir_aperture(valid.clone()).is_err());
        let mut grid = aperture_grid(0, 1); grid.step().unwrap();
        assert!(grid.with_reservoir_aperture(valid.clone()).is_err());
        let mut grid = aperture_grid(0, 1).with_reservoir_aperture(valid.clone()).unwrap();
        let before = grid.state_hash();
        assert!(grid.set_edge_velocity(2, false, [0.001, 0.0, 0.0]).is_err());
        assert!(grid.set_open_face(1, true, Some(valid.state)).is_err());
        assert_eq!(grid.state_hash(), before);
        assert!(grid.with_reservoir_aperture(valid).is_err());
    }

    #[test]
    fn g0_reservoir_aperture_is_worker_count_invariant() {
        let mut results = Vec::new();
        for workers in [1, 3] {
            let mut grid = aperture_grid(0, workers)
                .with_reservoir_aperture(aperture(1, false, 293.35)).unwrap();
            for _ in 0..8 { grid.step().unwrap(); }
            results.push((grid.state_hash(), grid.cumulative_ledger()));
        }
        assert_eq!(results[0], results[1]);
    }

    #[test]
    fn g0_open_faces_preserve_uniform_throughflow_all_axes() {
        for axis in 0..3 {
            for speed in [-10.0, 10.0] {
                let mut velocity = [0.0; 3];
                velocity[axis] = speed;
                let s = state(velocity, 303.15);
                let mut grid = uniform(axis, 1, s);
                let expected = grid.cell([1; 3]);
                for high in [false, true] {
                    grid.set_open_face(axis, high, Some(s)).unwrap();
                }
                let initial_energy = grid.total_energy_j();
                for _ in 0..32 {
                    let ledger = grid.step().unwrap();
                    assert!(ledger.residual_rel < 1e-12);
                }
                let mut max_error: f64 = 0.0;
                for i in 0..grid.n {
                    let (rho, u, temperature) = grid.cell(grid.at(i));
                    max_error = max_error
                        .max((rho - expected.0).abs())
                        .max((temperature / expected.2 - 1.0).abs());
                    for a in 0..3 {
                        max_error = max_error.max((u[a] - expected.1[a]).abs());
                    }
                }
                assert!(
                    max_error < 1e-11,
                    "axis={axis} speed={speed} error={max_error:e}"
                );
                assert!((grid.total_energy_j() / initial_energy - 1.0).abs() < 1e-12);
                // A sealed domain cannot satisfy this throughflow fixed point.
                let mut sealed = uniform(axis, 1, s);
                sealed.step().unwrap();
                assert!((sealed.cell([0; 3]).1[axis] - expected.1[axis]).abs() > 1e-4);
                println!("NSF_OPEN uniform axis={axis} speed_m_s={speed} error={max_error:e}");
            }
        }
    }

    #[test]
    fn g0_open_energy_flux_is_separate_from_wall_heat() {
        let cold = state([8.0, 0.0, 0.0], 293.15);
        let hot = state([8.0, 0.0, 0.0], 294.15);
        let mut grid = uniform(0, 1, cold);
        grid.set_open_face(0, false, Some(hot)).unwrap();
        grid.set_open_face(0, true, Some(cold)).unwrap();
        let before = grid.total_energy_j();
        for _ in 0..64 {
            let ledger = grid.step().unwrap();
            assert_eq!(ledger.wall_heat_j, 0.0);
            assert!(ledger.residual_rel < 1e-12, "{:?}", ledger);
        }
        let change = grid.total_energy_j() - before;
        assert!(change > 0.0);
        assert!(grid.cell([0, 1, 1]).2 > cold.temperature_k);
        let incoming = grid.cumulative_ledger().open_energy_j;
        assert!((change - incoming).abs() / before < 1e-12);
        println!(
            "NSF_OPEN heat change_j={change:e} open_energy_j={incoming:e} inlet_k={}",
            grid.cell([0, 1, 1]).2
        );
    }

    #[test]
    fn g0_open_refusals_and_updates_preserve_published_state() {
        let s = state([0.0; 3], 293.15);
        let mut grid = uniform(0, 1, s);
        let before = grid.state_hash();
        for axis in [1, 2, 3] {
            assert!(grid.set_open_face(axis, false, Some(s)).is_err());
        }
        for bad in [
            NsfOpenState {
                density_kg_m3: 0.0,
                ..s
            },
            NsfOpenState {
                density_kg_m3: f64::NAN,
                ..s
            },
            NsfOpenState {
                temperature_k: 400.0,
                ..s
            },
            NsfOpenState {
                temperature_k: f64::NAN,
                ..s
            },
            NsfOpenState {
                velocity_m_s: [1000.0, 0.0, 0.0],
                ..s
            },
            NsfOpenState {
                velocity_m_s: [f64::NAN, 0.0, 0.0],
                ..s
            },
        ] {
            assert!(grid.set_open_face(0, false, Some(bad)).is_err());
            assert_eq!(grid.state_hash(), before);
            assert!(grid.reservoir_faces[0][0].is_none());
        }
        let cells: Vec<_> = (0..grid.n).map(|i| grid.cell(grid.at(i))).collect();
        grid.set_open_face(0, false, Some(state([0.0; 3], 303.15)))
            .unwrap();
        for (i, expected) in cells.iter().enumerate() {
            let actual = grid.cell(grid.at(i));
            assert!((actual.0 - expected.0).abs() < 1e-14);
            assert!((actual.2 - expected.2).abs() < 1e-12);
            for a in 0..3 {
                assert!((actual.1[a] - expected.1[a]).abs() < 1e-14);
            }
        }
        assert!(grid.set_edge_velocity(0, false, [0.0, 1e-3, 0.0]).is_err());
        grid.set_open_face(0, false, None).unwrap();
        grid.set_edge_velocity(0, false, [0.0, 1e-3, 0.0]).unwrap();
        assert!(grid.set_open_face(0, false, Some(s)).is_err());
        let mut corner =
            NsfGrid3::<D3q27>::new(spec(0, 1), |x| x == [0; 3], |_| (1.0, [0.0; 3], 293.15))
                .unwrap();
        let before = corner.state_hash();
        assert!(corner.set_open_face(0, false, Some(s)).is_err());
        assert_eq!(corner.state_hash(), before);
        assert!(corner.reservoir_faces[0][0].is_none());
    }

    #[test]
    fn g0_reservoir_rebinding_preserves_physical_stress_and_energy() {
        use crate::nsf3::{Q, neq_second};

        for axis in 0..3 {
            let s = state([8.0, -4.0, 3.0], 293.15);
            let mut grid = uniform(axis, 1, s);
            // Nonzero shear and bulk stresses must survive the rebase too.
            // These Hermite coefficients have zero mass and first moment.
            let physical = [2e-4, -1e-4, 3e-4, 1.5e-4, -2.5e-4, 5e-5];
            let trace = physical[0] + physical[1] + physical[2];
            for f in &mut grid.f {
                for q in 0..Q {
                    let [x, y, z] = D3q27::E[q].map(f64::from);
                    let cpc = x * x * physical[0]
                        + y * y * physical[1]
                        + z * z * physical[2]
                        + 2.0 * (x * y * physical[3] + x * z * physical[4] + y * z * physical[5]);
                    f[q] += D3q27::W[q] * (cpc - THETA0 * trace) / (2.0 * THETA0.powi(2));
                }
            }
            let moments = |grid: &NsfGrid3<D3q27>| -> Vec<[f64; 6]> {
                (0..grid.n)
                    .map(|i| {
                        let (rho, u, _) = grid.cell(grid.at(i));
                        let eq = super::super::equilibrium::<D3q27>(rho, u, THETA0);
                        let mut second = neq_second::<D3q27>(&grid.f[i], &eq);
                        let f = grid.force[i];
                        second[0] += u[0] * f[0];
                        second[1] += u[1] * f[1];
                        second[2] += u[2] * f[2];
                        second[3] += 0.5 * (u[0] * f[1] + u[1] * f[0]);
                        second[4] += 0.5 * (u[0] * f[2] + u[2] * f[0]);
                        second[5] += 0.5 * (u[1] * f[2] + u[2] * f[1]);
                        second
                    })
                    .collect()
            };
            let before = moments(&grid);
            let cells: Vec<_> = (0..grid.n).map(|i| grid.cell(grid.at(i))).collect();
            let energy = grid.total_energy_j();
            for ghost in [
                Some(state(s.velocity_m_s, 303.15)),
                Some(state(s.velocity_m_s, 283.15)),
                None,
            ] {
                grid.set_open_face(axis, false, ghost).unwrap();
                let after = moments(&grid);
                for i in 0..grid.n {
                    for component in 0..6 {
                        assert!(
                            (after[i][component] - before[i][component]).abs() < 2e-14,
                            "axis={axis} cell={i} component={component}: {:?} -> {:?}",
                            before[i],
                            after[i]
                        );
                    }
                    let cell = grid.cell(grid.at(i));
                    assert!((cell.0 - cells[i].0).abs() < 2e-14);
                    assert_eq!(cell.2, cells[i].2);
                    for a in 0..3 {
                        assert!((cell.1[a] - cells[i].1[a]).abs() < 2e-14);
                    }
                }
                assert_eq!(grid.total_energy_j(), energy);
                if ghost.is_some() {
                    assert!(grid.force.iter().any(|f| f[axis].abs() > 1e-5));
                }
            }
        }
    }

    #[test]
    fn g0_local_density_controls_viscosity_admission_and_step_rollback() {
        use crate::nsf3::MAX_VERIFIED_LATTICE_VISCOSITY;

        let s = state([0.0; 3], 293.15);
        let mut grid = uniform(0, 1, s);
        let minimum_density = grid.law.mu_lattice(THETA0) / MAX_VERIFIED_LATTICE_VISCOSITY;
        let admitted = 1.01 * minimum_density;
        let refused = 0.99 * minimum_density;
        assert!(matches!(
            NsfGrid3::<D3q27>::new(spec(0, 1), |_| false, |_| (refused, [0.0; 3], 293.15)),
            Err(NsfError::Viscosity { .. })
        ));
        NsfGrid3::<D3q27>::new(spec(0, 1), |_| false, |_| (admitted, [0.0; 3], 293.15)).unwrap();

        // The same source law must reject a dilute prescribed reservoir
        // before publishing its boundary or any half-force rebase.
        let before = (grid.state_hash(), grid.total_energy_j());
        assert!(matches!(
            grid.set_open_face(
                0,
                false,
                Some(NsfOpenState {
                    density_kg_m3: refused * grid.law.rho_ref_kg_m3,
                    ..s
                })
            ),
            Err(NsfError::Viscosity { .. })
        ));
        assert_eq!((grid.state_hash(), grid.total_energy_j()), before);
        assert!(grid.reservoir_faces[0][0].is_none());

        // A constant-temperature manufactured rarefaction crosses the local
        // envelope after streaming. Match removed energy to removed mass so
        // that a temperature-domain check cannot mask the viscosity defect.
        let removed = refused - admitted;
        let cv_theta = THETA0 / (grid.law.gamma - 1.0);
        let mut periodic = spec(0, 1);
        periodic.periodic = [true; 3];
        let mut evolving = NsfGrid3::<D3q27>::new_with_sources(
            periodic,
            |_| false,
            |_| (admitted, [0.0; 3], 293.15),
            Some(Box::new(move |_, _| {
                (removed, [0.0; 3], removed * cv_theta)
            })),
        )
        .unwrap();
        let before = (
            evolving.state_hash(),
            evolving.steps(),
            evolving.cumulative_ledger(),
        );
        assert!(matches!(evolving.step(), Err(NsfError::Viscosity { .. })));
        assert_eq!(
            (
                evolving.state_hash(),
                evolving.steps(),
                evolving.cumulative_ledger()
            ),
            before
        );
    }

    #[test]
    fn g0_prescribed_and_characteristic_faces_are_mutually_exclusive() {
        use crate::nsf3::{NsfCharacteristicAxis, NsfCharacteristicBoundary};
        let s = state([0.0; 3], 293.15);
        let face = NsfCharacteristicBoundary {
            reference_pressure_pa: 101325.0,
            reference_temperature_k: s.temperature_k,
            reference_velocity_m_s: s.velocity_m_s,
            incoming_characteristic_m_s: 0.0,
        };
        let axis = NsfCharacteristicAxis {
            axis: 0,
            faces: [face; 2],
            maximum_relative_perturbation: 0.05,
        };
        let mut prescribed = uniform(0, 1, s);
        prescribed.set_open_face(0, false, Some(s)).unwrap();
        assert!(prescribed.with_characteristic_axis(axis).is_err());
        let mut characteristic = uniform(0, 1, s).with_characteristic_axis(axis).unwrap();
        let before = characteristic.state_hash();
        let traces = characteristic.open_boundary_traces(false);
        assert!(characteristic.set_open_face(0, false, Some(s)).is_err());
        assert!(characteristic.set_open_face(0, true, None).is_err());
        assert_eq!(characteristic.state_hash(), before);
        assert_eq!(characteristic.open_boundary_traces(false), traces);
    }

    #[test]
    fn g5_open_face_evolution_is_worker_count_invariant() {
        let s = state([8.0, 0.0, 0.0], 293.15);
        let hot = state([8.0, 0.0, 0.0], 294.15);
        let mut hashes = Vec::new();
        for workers in [1, 3] {
            let mut grid = uniform(0, workers, s);
            grid.set_open_face(0, false, Some(hot)).unwrap();
            grid.set_open_face(0, true, Some(s)).unwrap();
            for _ in 0..8 {
                grid.step().unwrap();
            }
            hashes.push(grid.state_hash());
        }
        assert_eq!(hashes[0], hashes[1]);
    }

    #[test]
    fn g2_prescribed_open_face_reflection_diagnostic() {
        // Fixed physical pulse and domain; only dx changes. Prediction: a
        // prescribed exterior state is not a characteristic-matched outlet.
        // This measures its reflection; the test does NOT assert transparency.
        for n in [64, 128, 256] {
            let mut s = spec(0, 1);
            s.dims = [n, 1, 1];
            s.maximum_cells = n;
            s.dx_m = 0.02 / n as f64;
            let gas = GasState::try_new(&s.gas, s.temperature_k, s.pressure_pa).unwrap();
            let c = (gas.gamma * THETA0).sqrt();
            let amplitude = 1e-6;
            let mut grid = NsfGrid3::<D3q27>::new(
                s,
                |_| false,
                |x| {
                    let position = (x[0] as f64 + 0.5) / n as f64;
                    let a = amplitude * fs_math::det::exp(-0.5 * ((position - 0.3) / 0.05).powi(2));
                    (
                        1.0 + a,
                        [c * a, 0.0, 0.0],
                        293.15 * (1.0 + (gas.gamma - 1.0) * a),
                    )
                },
            )
            .unwrap();
            for high in [false, true] {
                grid.set_open_face(0, high, Some(state([0.0; 3], 293.15)))
                    .unwrap();
            }
            let mut incident: f64 = 0.0;
            let mut reflected: f64 = 0.0;
            for step in 0..(1.5 * n as f64 / c).ceil() as usize {
                grid.step().unwrap();
                let (rho, u, temperature) = grid.cell([n / 2, 0, 0]);
                let pressure = rho * grid.law.theta_of_kelvin(temperature) - THETA0;
                let travel = step as f64 * c / n as f64;
                if travel < 0.5 {
                    incident = incident.max((pressure + c * u[0]).abs());
                }
                if travel > 0.8 {
                    reflected = reflected.max((pressure - c * u[0]).abs());
                }
            }
            assert!(incident > amplitude * 0.5 * c * c);
            assert!(reflected.is_finite());
            println!(
                "NSF_OPEN prescribed_reflection cells={n} dx_m={} amplitude_ratio={:e} transparency_qualified=false",
                grid.law.dx_m,
                reflected / incident
            );
        }
    }
}
