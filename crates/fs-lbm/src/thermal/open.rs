//! Conservative D3Q7 temperature transport for native open 3-D grids.
//!
//! Each reservoir supplies a temperature even on a nominal outlet: the
//! accepted, force-corrected normal velocity chooses inflow or outflow at
//! every boundary cell. Inflow uses the same Danckwerts total-flux condition
//! as `nsf3::open_boundary::enthalpy`; outflow has zero normal conduction.
//! Receipts measure the populations actually admitted at each open link.
//!
//! This is the fixed reference volumetric heat-capacity, low-Mach scalar
//! approximation. Summed temperature times the common cell-capacity and
//! temperature scales gives energy. It is not compressible mass-specific
//! enthalpy transport, and its receipt is not the momentum reservoir's mass
//! receipt times temperature. Velocity belongs to the beginning of the heat
//! interval, as in the existing D3Q7 collision/history kernel.

use super::advection::{D3Q7, FluxHistory};
use super::boundary::conduction::{WallConduction, WallConductionError};
use super::boundary::{
    FiniteThermalWallExchange, ThermalWall, ThermalWallError, ThermalWallExchange, WallReflection,
};

/// Implicit conduction between existing finite walls of the native open grid.
pub mod conduction;

/// Heat ownership on one native face, in XMin/XMax/YMin/YMax/ZMin/ZMax order.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ThermalOpenFace3 {
    /// Match the native periodic axis exactly.
    Periodic,
    /// Match a stationary native wall, retaining the existing scalar law.
    Wall(ThermalWall),
    /// An explicit reservoir temperature, including every possible backflow.
    Reservoir { temperature: f64 },
}

/// No thermal inlet, outlet or obstacle condition is inferred from momentum.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ThermalOpenBoundaries3 {
    pub faces: [ThermalOpenFace3; 6],
    /// Law on fluid-to-solid links. Finite solid cells require independent
    /// link identity and are supported on external planes only by this owner.
    pub internal_wall: ThermalWall,
}

/// Fixed resource and numerical envelope for one thermal carrier.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ThermalOpenLimits3 {
    pub maximum_cells: usize,
    /// Bounds all persistent thermal arrays, including retained source/flux
    /// scratch and both copies of every finite wall state. Native momentum
    /// storage and allocator metadata are outside this thermal-only budget.
    /// Optional wall graphs have additional explicit per-pair limits through
    /// `WallConductionControls`; their factors/workspace are not counted here.
    pub maximum_storage_bytes: usize,
    /// Admitted fluid temperatures, reservoir data and resolved wall traces.
    pub temperature_range: [f64; 2],
    /// Absolute tolerance on summed scalar conservation and solid-link defects.
    pub balance_tolerance: f64,
}

/// Measured incoming-minus-outgoing scalar populations on an open plane.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct ThermalReservoirExchange3 {
    pub heat_in: f64,
    /// Signed contributions selected by actual local velocity, not face labels.
    pub inflow_heat: f64,
    pub outflow_heat: f64,
    /// Positive magnitudes of inward/outward velocity summed over face links.
    /// These are reference-volume fluxes, not native compressible mass fluxes.
    pub inward_volume: f64,
    pub outward_volume: f64,
    pub links: usize,
}

/// Independent first-law account of one accepted thermal interval.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct ThermalOpenBalance3 {
    pub source: f64,
    pub stored_change: f64,
    pub reservoirs: [ThermalReservoirExchange3; 6],
    pub walls: [Option<ThermalWallExchange>; 6],
    pub internal_wall: Option<ThermalWallExchange>,
    /// Fluid storage minus source, actual reservoir heat and actual wall heat.
    pub defect: f64,
    /// Fluid defect plus every finite-solid defect. Reciprocal fluid/solid
    /// exchange cancels, leaving only authored external sources and reservoirs.
    pub combined_defect: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub enum ThermalOpenError3 {
    Invalid(&'static str),
    Allocation,
    Wall(ThermalWallError),
    Conduction(WallConductionError),
    Temperature {
        cell: usize,
        value: f64,
    },
    BoundaryTemperature {
        face: usize,
        value: f64,
    },
    Velocity {
        cell: usize,
        axis: usize,
        value: f64,
    },
    Balance {
        defect: f64,
    },
    Cancelled,
}

impl core::fmt::Display for ThermalOpenError3 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native open thermal transport: {self:?}")
    }
}
impl std::error::Error for ThermalOpenError3 {}

/// Accepted heat populations and reusable, bounded transactional workspace.
/// Construct through `OpenBoundaryGrid3::prepare_thermal_field` so geometry,
/// solid links, initial physical velocity and completed time agree.
pub struct ThermalOpenField3 {
    pub(crate) dims: [usize; 3],
    pub(crate) solid: Vec<bool>,
    boundaries: ThermalOpenBoundaries3,
    limits: ThermalOpenLimits3,
    tau: f64,
    populations: Vec<[f64; 7]>,
    trial: Vec<[f64; 7]>,
    post: Vec<f64>,
    pub(crate) velocities: Vec<[f64; 3]>,
    history: FluxHistory<3>,
    trial_history: FluxHistory<3>,
    walls: [Vec<ThermalWall>; 6],
    trial_walls: [Vec<ThermalWall>; 6],
    conduction: [Option<WallConduction>; 3],
    completed_steps: usize,
    last_balance: Option<ThermalOpenBalance3>,
}

fn invalid(reason: &'static str) -> ThermalOpenError3 {
    ThermalOpenError3::Invalid(reason)
}

fn filled<T: Clone>(len: usize, value: T) -> Result<Vec<T>, ThermalOpenError3> {
    let mut data = Vec::new();
    data.try_reserve_exact(len)
        .map_err(|_| ThermalOpenError3::Allocation)?;
    data.resize(len, value);
    Ok(data)
}

#[derive(Default, Clone, Copy)]
struct Sum {
    value: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) {
        let next = self.value + value;
        self.correction += if self.value.abs() >= value.abs() {
            (self.value - next) + value
        } else {
            (value - next) + self.value
        };
        self.value = next;
    }
    fn total(self) -> f64 {
        self.value + self.correction
    }
}

#[derive(Default, Clone, Copy)]
struct WallSum {
    heat: Sum,
    temperature: Sum,
    minimum: f64,
    maximum: f64,
    links: usize,
    solid_heat: Sum,
    solid_storage: Sum,
    solid_source: Sum,
    solid_conducted: Sum,
    maximum_defect: f64,
    has_solid: bool,
}
impl WallSum {
    fn add(&mut self, outgoing: f64, response: &WallReflection) {
        self.heat.add(response.incoming - outgoing);
        self.temperature.add(response.temperature);
        if self.links == 0 {
            self.minimum = response.temperature;
            self.maximum = response.temperature;
        } else {
            self.minimum = self.minimum.min(response.temperature);
            self.maximum = self.maximum.max(response.temperature);
        }
        self.links += 1;
        if let Some(solid) = response.solid {
            self.has_solid = true;
            self.solid_heat.add(solid.heat_in);
            self.solid_storage.add(solid.stored_change);
            self.solid_source.add(solid.source);
            self.solid_conducted.add(solid.conducted_in);
            self.maximum_defect = self.maximum_defect.max(solid.maximum_link_defect);
        }
    }
    fn finish(self) -> Option<ThermalWallExchange> {
        (self.links != 0).then(|| ThermalWallExchange {
            heat_in: self.heat.total(),
            mean_temperature: self.temperature.total() / self.links as f64,
            min_temperature: self.minimum,
            max_temperature: self.maximum,
            solid: self.has_solid.then(|| FiniteThermalWallExchange {
                heat_in: self.solid_heat.total(),
                stored_change: self.solid_storage.total(),
                source: self.solid_source.total(),
                conducted_in: self.solid_conducted.total(),
                maximum_link_defect: self.maximum_defect,
            }),
        })
    }
}

/// Existing Danckwerts convention: total inward flux u*T_res for inflow,
/// zero normal conductive flux for outflow. The latter solves
/// gin-gout = u*(gin+gout)/(2w). No prescribed outflow temperature is inserted.
fn reservoir_incoming(outgoing: f64, inward_velocity: f64, temperature: f64) -> f64 {
    if inward_velocity > 0.0 {
        inward_velocity.mul_add(temperature, outgoing)
    } else if inward_velocity < 0.0 {
        // The admitted |u| < 2w makes this ratio positive and <= 1. Form it
        // before multiplying: an intermediate gout*(2w+u) need not be finite.
        let pair_weight = 2.0 * D3Q7.weights[1];
        outgoing * ((pair_weight + inward_velocity) / (pair_weight - inward_velocity))
    } else {
        outgoing
    }
}

impl ThermalOpenField3 {
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn try_new<E: From<ThermalOpenError3>>(
        dims: [usize; 3],
        solid: &[bool],
        initial: &[f64],
        tau: f64,
        boundaries: ThermalOpenBoundaries3,
        limits: ThermalOpenLimits3,
        completed_steps: usize,
        mut velocity: impl FnMut(usize) -> Result<[f64; 3], E>,
        mut poll: impl FnMut() -> Result<(), E>,
    ) -> Result<Self, E> {
        poll()?;
        let cells = dims
            .into_iter()
            .try_fold(1usize, |n, d| if d == 0 { None } else { n.checked_mul(d) })
            .ok_or_else(|| invalid("empty or overflowing thermal dimensions"))?;
        if cells > limits.maximum_cells || solid.len() != cells || initial.len() != cells {
            return Err(invalid("thermal cell cap or initial field shape mismatch").into());
        }
        if !tau.is_finite()
            || tau <= 0.5
            || !limits.temperature_range.iter().all(|v| v.is_finite())
            || limits.temperature_range[0] > limits.temperature_range[1]
            || !limits.balance_tolerance.is_finite()
            || limits.balance_tolerance < 0.0
        {
            return Err(
                invalid("invalid relaxation, temperature interval or balance tolerance").into(),
            );
        }
        boundaries
            .internal_wall
            .validate()
            .map_err(ThermalOpenError3::Wall)?;
        if let ThermalWall::Temperature(value) = boundaries.internal_wall {
            check_boundary_temperature(6, value, limits)?;
        }
        if matches!(boundaries.internal_wall, ThermalWall::FiniteSolid(_)) {
            return Err(invalid("internal finite solids need individually owned links").into());
        }
        let mut wall_entries = 0usize;
        for (face, law) in boundaries.faces.iter().enumerate() {
            match *law {
                ThermalOpenFace3::Wall(wall) => {
                    wall.validate().map_err(ThermalOpenError3::Wall)?;
                    if let ThermalWall::Temperature(value) = wall {
                        check_boundary_temperature(face, value, limits)?;
                    }
                    wall_entries = wall_entries
                        .checked_add(cells / dims[face / 2])
                        .ok_or_else(|| invalid("thermal wall storage overflow"))?;
                }
                ThermalOpenFace3::Reservoir { temperature } => {
                    check_boundary_temperature(face, temperature, limits)?;
                }
                ThermalOpenFace3::Periodic => {}
            }
            if matches!(law, ThermalOpenFace3::Periodic)
                != matches!(boundaries.faces[face ^ 1], ThermalOpenFace3::Periodic)
            {
                return Err(invalid("periodic thermal faces must be paired").into());
            }
        }
        let per_cell = 2 * core::mem::size_of::<[f64; 7]>()
            + 7 * core::mem::size_of::<f64>()
            + 3 * core::mem::size_of::<[f64; 3]>()
            + core::mem::size_of::<bool>();
        let bytes = cells
            .checked_mul(per_cell)
            .and_then(|n| {
                wall_entries
                    .checked_mul(2 * core::mem::size_of::<ThermalWall>())
                    .and_then(|walls| n.checked_add(walls))
            })
            .ok_or_else(|| invalid("thermal storage overflow"))?;
        if bytes > limits.maximum_storage_bytes {
            return Err(invalid("thermal storage budget exceeded").into());
        }
        let mut field = Self {
            dims,
            solid: filled(cells, false)?,
            boundaries,
            limits,
            tau,
            populations: filled(cells, [0.0; 7])?,
            trial: filled(cells, [0.0; 7])?,
            post: filled(
                cells
                    .checked_mul(7)
                    .ok_or_else(|| invalid("population storage overflow"))?,
                0.0,
            )?,
            velocities: filled(cells, [0.0; 3])?,
            history: FluxHistory::default(),
            trial_history: FluxHistory::default(),
            walls: core::array::from_fn(|_| Vec::new()),
            trial_walls: core::array::from_fn(|_| Vec::new()),
            completed_steps,
            last_balance: None,
            conduction: core::array::from_fn(|_| None),
        };
        field.solid.copy_from_slice(solid);
        field
            .history
            .try_prepare(cells)
            .map_err(|_| ThermalOpenError3::Allocation)?;
        field
            .trial_history
            .try_prepare(cells)
            .map_err(|_| ThermalOpenError3::Allocation)?;
        for (face, law) in boundaries.faces.iter().enumerate() {
            if let ThermalOpenFace3::Wall(wall) = law {
                let links = cells / dims[face / 2];
                field.walls[face] = filled(links, *wall)?;
                field.trial_walls[face] = filled(links, *wall)?;
                check_solid_temperature(face, *wall, limits)?;
            }
        }
        for cell in 0..cells {
            if cell % 64 == 0 {
                poll()?;
            }
            if !solid[cell] {
                check_temperature(cell, initial[cell], limits)?;
                let u = velocity(cell)?;
                check_velocity(cell, u)?;
                field.velocities[cell] = u;
                field.populations[cell] = D3Q7.equilibrium(initial[cell], u);
            }
        }
        poll()?;
        Ok(field)
    }

    pub fn populations(&self) -> &[[f64; 7]] {
        &self.populations
    }
    pub fn dimensions(&self) -> [usize; 3] {
        self.dims
    }
    pub fn boundaries(&self) -> ThermalOpenBoundaries3 {
        self.boundaries
    }
    pub fn completed_steps(&self) -> usize {
        self.completed_steps
    }
    pub fn last_balance(&self) -> Option<ThermalOpenBalance3> {
        self.last_balance
    }
    /// Current external wall cells in x-fast surface order. The authored
    /// declarations retain their initial values; evolving solids live here.
    /// Slots adjacent to solid grid cells are inactive and do not consume heat.
    pub fn wall_states(&self, face: usize) -> Option<&[ThermalWall]> {
        self.walls.get(face).map(Vec::as_slice)
    }
    pub fn temperature(&self, cell: usize) -> Option<f64> {
        self.populations
            .get(cell)
            .filter(|_| !self.solid[cell])
            .map(|g| g.iter().sum())
    }

    /// Replace a reservoir temperature without stepping either physical field.
    /// Face kind is preserved; the value is admitted even if currently outflow.
    pub fn set_reservoir_temperature(
        &mut self,
        face: usize,
        temperature: f64,
    ) -> Result<(), ThermalOpenError3> {
        if !matches!(
            self.boundaries.faces.get(face),
            Some(ThermalOpenFace3::Reservoir { .. })
        ) {
            return Err(invalid("thermal target is not a reservoir face"));
        }
        check_boundary_temperature(face, temperature, self.limits)?;
        self.boundaries.faces[face] = ThermalOpenFace3::Reservoir { temperature };
        Ok(())
    }

    pub(crate) fn begin_step<E: From<ThermalOpenError3>>(
        &mut self,
        heat: &[f64],
        thermal_tau: Option<&[f64]>,
        mut poll: impl FnMut() -> Result<(), E>,
    ) -> Result<ThermalOpenBalance3, E> {
        let cells = self.populations.len();
        if heat.len() != cells || thermal_tau.is_some_and(|tau| tau.len() != cells) {
            return Err(invalid("thermal source or relaxation field shape mismatch").into());
        }
        let mut source = Sum::default();
        for cell in 0..cells {
            if cell % 64 == 0 {
                poll()?;
            }
            if !heat[cell].is_finite() || (self.solid[cell] && heat[cell] != 0.0) {
                return Err(
                    invalid("nonfinite source or heat assigned to a solid grid cell").into(),
                );
            }
            if thermal_tau.is_some_and(|tau| !tau[cell].is_finite() || tau[cell] <= 0.5) {
                return Err(invalid(
                    "local thermal relaxation must be finite and greater than one half",
                )
                .into());
            }
            if !self.solid[cell] {
                check_velocity(cell, self.velocities[cell])?;
                source.add(heat[cell]);
            } else {
                self.velocities[cell] = [0.0; 3];
            }
        }
        self.trial.copy_from_slice(&self.populations);
        self.trial_history.copy_from(&self.history);
        for face in 0..6 {
            self.trial_walls[face].copy_from_slice(&self.walls[face]);
        }
        let velocities = &self.velocities;
        D3Q7.source_half_kick(
            &mut self.trial,
            heat,
            0..cells,
            &mut self.trial_history,
            |cell| velocities[cell],
        );
        self.check_candidate(&mut poll)?;
        poll()?;
        D3Q7.collide(
            &self.trial,
            &mut self.post,
            0..cells,
            self.tau,
            thermal_tau,
            &mut self.trial_history,
            |cell| velocities[cell],
        );
        poll()?;
        self.prepare_wall_conduction()?;
        poll()?;
        let mut balance = self.stream(&mut poll)?;
        balance.source = source.total();
        self.check_candidate(&mut poll)?;
        Ok(balance)
    }

    fn stream<E: From<ThermalOpenError3>>(
        &mut self,
        mut poll: impl FnMut() -> Result<(), E>,
    ) -> Result<ThermalOpenBalance3, E> {
        let mut balance = ThermalOpenBalance3::default();
        let mut walls = [WallSum::default(); 6];
        let mut internal = WallSum::default();
        let mut reservoir_heat = [Sum::default(); 6];
        let mut inflow_heat = [Sum::default(); 6];
        let mut outflow_heat = [Sum::default(); 6];
        let mut inward_volume = [Sum::default(); 6];
        let mut outward_volume = [Sum::default(); 6];
        let strides = [1, self.dims[0], self.dims[0] * self.dims[1]];
        for cell in 0..self.populations.len() {
            if cell % 64 == 0 {
                poll()?;
            }
            if self.solid[cell] {
                continue;
            }
            self.trial[cell][0] = self.post[7 * cell];
            for axis in 0..3 {
                let coordinate = cell / strides[axis] % self.dims[axis];
                let (positive, negative) = D3Q7.direction_pair(axis);
                for (q, opposite, upper) in
                    [(positive, negative, false), (negative, positive, true)]
                {
                    let outside = if upper {
                        coordinate + 1 == self.dims[axis]
                    } else {
                        coordinate == 0
                    };
                    let face = 2 * axis + usize::from(upper);
                    let outgoing = self.post[7 * cell + opposite];
                    let weight = D3Q7.weights[q];
                    let incoming = if outside
                        && !matches!(self.boundaries.faces[face], ThermalOpenFace3::Periodic)
                    {
                        match self.boundaries.faces[face] {
                            ThermalOpenFace3::Reservoir { temperature } => {
                                let u = if upper {
                                    -self.velocities[cell][axis]
                                } else {
                                    self.velocities[cell][axis]
                                };
                                let incoming = reservoir_incoming(outgoing, u, temperature);
                                let trace = (0.5 * incoming + 0.5 * outgoing) / weight;
                                check_boundary_temperature(face, trace, self.limits)?;
                                let heat = incoming - outgoing;
                                reservoir_heat[face].add(heat);
                                if u > 0.0 {
                                    inflow_heat[face].add(heat);
                                    inward_volume[face].add(u);
                                } else if u < 0.0 {
                                    outflow_heat[face].add(heat);
                                    outward_volume[face].add(-u);
                                }
                                balance.reservoirs[face].links += 1;
                                incoming
                            }
                            ThermalOpenFace3::Wall(_) => {
                                let plane = cell / (strides[axis] * self.dims[axis])
                                    * strides[axis]
                                    + cell % strides[axis];
                                let law = &mut self.trial_walls[face][plane];
                                check_solid_temperature(face, *law, self.limits)?;
                                let response = if let Some(conduction) = &self.conduction[axis] {
                                    conduction.reflect(
                                        usize::from(upper),
                                        plane,
                                        law,
                                        outgoing,
                                        weight,
                                    )
                                } else {
                                    law.reflect_and_advance(outgoing, weight)
                                };
                                check_solid_temperature(face, *law, self.limits)?;
                                check_wall_response(face, &response, self.limits)?;
                                walls[face].add(outgoing, &response);
                                response.incoming
                            }
                            ThermalOpenFace3::Periodic => unreachable!(),
                        }
                    } else {
                        let from = if outside {
                            if upper {
                                cell - (self.dims[axis] - 1) * strides[axis]
                            } else {
                                cell + (self.dims[axis] - 1) * strides[axis]
                            }
                        } else if upper {
                            cell + strides[axis]
                        } else {
                            cell - strides[axis]
                        };
                        if self.solid[from] {
                            let mut law = self.boundaries.internal_wall;
                            let response = law.reflect_and_advance(outgoing, weight);
                            check_wall_response(6, &response, self.limits)?;
                            internal.add(outgoing, &response);
                            response.incoming
                        } else {
                            self.post[7 * from + q]
                        }
                    };
                    self.trial[cell][q] = incoming;
                }
            }
        }
        for face in 0..6 {
            balance.walls[face] = walls[face].finish();
            let exchange = &mut balance.reservoirs[face];
            exchange.heat_in = reservoir_heat[face].total();
            exchange.inflow_heat = inflow_heat[face].total();
            exchange.outflow_heat = outflow_heat[face].total();
            exchange.inward_volume = inward_volume[face].total();
            exchange.outward_volume = outward_volume[face].total();
        }
        balance.internal_wall = internal.finish();
        Ok(balance)
    }

    pub(crate) fn finish_step<E: From<ThermalOpenError3>>(
        &mut self,
        heat: &[f64],
        mut balance: ThermalOpenBalance3,
        mut poll: impl FnMut() -> Result<(), E>,
    ) -> Result<ThermalOpenBalance3, E> {
        let cells = self.populations.len();
        for cell in 0..cells {
            if cell % 64 == 0 {
                poll()?;
            }
            if !self.solid[cell] {
                check_velocity(cell, self.velocities[cell])?;
            }
        }
        let velocities = &self.velocities;
        D3Q7.source_half_kick(
            &mut self.trial,
            heat,
            0..cells,
            &mut self.trial_history,
            |cell| velocities[cell],
        );
        self.check_candidate(&mut poll)?;
        let mut stored = Sum::default();
        for cell in 0..cells {
            if cell % 64 == 0 {
                poll()?;
            }
            if !self.solid[cell] {
                for q in 0..7 {
                    stored.add(self.trial[cell][q] - self.populations[cell][q]);
                }
            }
        }
        balance.stored_change = stored.total();
        for reservoir in balance.reservoirs {
            if [
                reservoir.heat_in,
                reservoir.inflow_heat,
                reservoir.outflow_heat,
                reservoir.inward_volume,
                reservoir.outward_volume,
            ]
            .iter()
            .any(|v| !v.is_finite())
            {
                return Err(invalid("nonfinite reservoir heat account").into());
            }
        }
        for wall in balance
            .walls
            .into_iter()
            .chain([balance.internal_wall])
            .flatten()
        {
            if [
                wall.heat_in,
                wall.mean_temperature,
                wall.min_temperature,
                wall.max_temperature,
            ]
            .iter()
            .any(|v| !v.is_finite())
            {
                return Err(invalid("nonfinite wall heat account").into());
            }
            if let Some(solid) = wall.solid {
                let defect = solid.defect();
                if !defect.is_finite() || defect.abs() > self.limits.balance_tolerance {
                    return Err(ThermalOpenError3::Balance { defect }.into());
                }
            }
        }
        let mut defect = Sum::default();
        defect.add(balance.stored_change);
        defect.add(-balance.source);
        for reservoir in balance.reservoirs {
            defect.add(-reservoir.heat_in);
        }
        for wall in balance
            .walls
            .into_iter()
            .chain([balance.internal_wall])
            .flatten()
        {
            defect.add(-wall.heat_in);
        }
        balance.defect = defect.total();
        let mut combined = Sum::default();
        combined.add(balance.defect);
        for wall in balance
            .walls
            .into_iter()
            .chain([balance.internal_wall])
            .flatten()
        {
            if let Some(solid) = wall.solid {
                combined.add(solid.defect());
            }
        }
        balance.combined_defect = combined.total();
        for value in [balance.defect, balance.combined_defect] {
            if !value.is_finite() || value.abs() > self.limits.balance_tolerance {
                return Err(ThermalOpenError3::Balance { defect: value }.into());
            }
        }
        poll()?;
        Ok(balance)
    }

    fn check_candidate<E: From<ThermalOpenError3>>(
        &self,
        mut poll: impl FnMut() -> Result<(), E>,
    ) -> Result<(), E> {
        for (cell, populations) in self.trial.iter().enumerate() {
            if cell % 64 == 0 {
                poll()?;
            }
            if !self.solid[cell] {
                if populations.iter().any(|v| !v.is_finite()) {
                    return Err(invalid("nonfinite thermal population").into());
                }
                check_temperature(cell, populations.iter().sum(), self.limits)?;
            }
        }
        Ok(())
    }

    pub(crate) fn publish(&mut self, step: usize, balance: ThermalOpenBalance3) {
        core::mem::swap(&mut self.populations, &mut self.trial);
        core::mem::swap(&mut self.history, &mut self.trial_history);
        core::mem::swap(&mut self.walls, &mut self.trial_walls);
        self.completed_steps = step;
        self.last_balance = Some(balance);
    }
}

fn check_velocity(cell: usize, velocity: [f64; 3]) -> Result<(), ThermalOpenError3> {
    for (axis, value) in velocity.into_iter().enumerate() {
        if !value.is_finite() || value.abs() >= D3Q7.cs2 {
            return Err(ThermalOpenError3::Velocity { cell, axis, value });
        }
    }
    Ok(())
}
fn check_temperature(
    cell: usize,
    value: f64,
    limits: ThermalOpenLimits3,
) -> Result<(), ThermalOpenError3> {
    if !value.is_finite()
        || value < limits.temperature_range[0]
        || value > limits.temperature_range[1]
    {
        return Err(ThermalOpenError3::Temperature { cell, value });
    }
    Ok(())
}
fn check_boundary_temperature(
    face: usize,
    value: f64,
    limits: ThermalOpenLimits3,
) -> Result<(), ThermalOpenError3> {
    if !value.is_finite()
        || value < limits.temperature_range[0]
        || value > limits.temperature_range[1]
    {
        return Err(ThermalOpenError3::BoundaryTemperature { face, value });
    }
    Ok(())
}
fn check_solid_temperature(
    face: usize,
    law: ThermalWall,
    limits: ThermalOpenLimits3,
) -> Result<(), ThermalOpenError3> {
    if let ThermalWall::FiniteSolid(solid) = law {
        check_boundary_temperature(face, solid.temperature, limits)?;
    }
    Ok(())
}
fn check_wall_response(
    face: usize,
    response: &WallReflection,
    limits: ThermalOpenLimits3,
) -> Result<(), ThermalOpenError3> {
    check_boundary_temperature(face, response.temperature, limits)?;
    if let Some(solid) = response.solid {
        for value in [
            solid.heat_in,
            solid.stored_change,
            solid.source,
            solid.conducted_in,
            solid.maximum_link_defect,
            solid.defect(),
        ] {
            if !value.is_finite() {
                return Err(invalid("nonfinite finite-wall heat account"));
            }
        }
        if solid.maximum_link_defect > limits.balance_tolerance {
            return Err(ThermalOpenError3::Balance {
                defect: solid.maximum_link_defect,
            });
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests;
