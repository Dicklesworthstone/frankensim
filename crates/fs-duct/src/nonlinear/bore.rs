//! A sounding bore on the shared nonlinear gas owner.
//!
//! [`NonlinearBore`] keeps one [`Euler1d`] field of the original NASA-9 source
//! and adds what a played wind bore needs that the inviscid duct lacks:
//!
//! - Radius steps inside cells. Cell volumes are exact; the step's wall force
//!   enters through the same quasi-one-dimensional `p dA` source as a taper.
//!   Like the plane-wave transfer-matrix chain, this resolves no evanescent
//!   step modes.
//! - Zwikker–Kosten wall memory per cell: the `fs_phs::ZkDiffusive` drag and
//!   diffusive arms of the local radius, applied between accepted steps by
//!   their exact propagator (a symmetric negative semidefinite generator, so
//!   the update is contractive). Wall friction turns kinetic into internal
//!   energy; wall heat exchange follows temperature relative to the room wall,
//!   changing internal energy at the cell's current density.
//! - A radiating end: a rigid circular piston in an infinite baffle. Its
//!   radiation impedance has the exact impulse response
//!   `z(t) = (rho c / S) [delta(t) - (4/(pi T)) sqrt(1 - (t/T)^2)]`, `0 <= t <= T = 2a/c`
//!   (Stepanishen, JASA 49 (1971) 841), the transform of
//!   `rho c / S [1 - J1(2ka)/ka + i H1(2ka)/ka]`. It is product-integrated
//!   against the piecewise-linear end flow on the fixed step clock and solved
//!   implicitly with the interior's own end wave at every RK stage.
//!   Optionally the end takes instead the plane-mode equivalent of the complete
//!   flanged multimodal mouth ([`MouthKernel`], [`NonlinearBore::with_mouth_kernel`]):
//!   the same pressure law `p = Z_inf (U_t - sum_j w_j U_{t-j})` with the kernel's
//!   causal real weights, so the plane gas still sees one impedance, now the one a
//!   plane wave at the mouth meets once the higher mouth modes are counted.
//!
//! - Optionally, a breathing wall ([`NonlinearBore::bind_shell`]): one radial
//!   oscillator per cell, driven by the cell pressure, that changes the volume
//!   the gas fills (the wall of a flaring bell moves under the bore pressure;
//!   see the `shell` module for the law, the gas update and the energy account).
//!
//! The driven end is a caller-owned port. Walls are rigid unless a shell is
//! bound, the room is at rest, and the bell radiates as a baffled piston: no
//! flare-mode or unflanged-lip radiation is claimed.

use super::{NonlinearDuctError, area, positive, profile_volume, radius, reserved, unsupported};
use crate::modal::MouthKernel;
use crate::{Duct, DuctError, Segment};
use fs_gas::dynamics::{
    Boundaries, Boundary, Conserved, Error as GasError, Euler1d, FaceFlux, FrozenGas, Mesh1d,
    PortStage, Primitive, Reconstruction, StepControl, StepReport,
};
use fs_material::fluid_transport::FluidTransportState;
use fs_math::det;
use fs_phs::ZkDiffusiveBasis;
use fs_qty::{Pressure, Temperature, Time};

mod shell;
use shell::{DrivenWall, Shell};
pub use shell::{ShellLaw, ShellLedger, ShellReceipt, ShellStrip, WallDriver, WallStepInput};

/// Zwikker–Kosten wall-memory admission for every cell.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WallMemory {
    /// Viscous and thermal diffusive poles per cell (`fs_phs` cap applies).
    pub poles: usize,
    /// Admitted acoustic band [Hz]; the binding checks both wall operators on it.
    pub band_hz: [f64; 2],
    /// Largest sampled relative error of either wall operator on the band.
    pub maximum_relative_error: f64,
}

/// Transport at the stationary room/wall reference state.
///
/// Density and gamma must match the field's original caloric source at the
/// room temperature and pressure. [`NonlinearBore::new`] checks this before
/// binding the wall memory. Viscosity and Prandtl number remain caller-owned
/// transport data at that same state; [`Self::from_source`] derives Prandtl
/// from the original Cp(T) and a supplied thermal conductivity.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WallGas {
    /// Density [kg/m³].
    pub density_kg_m3: f64,
    /// Dynamic viscosity [Pa s].
    pub dynamic_viscosity_pa_s: f64,
    /// Heat-capacity ratio.
    pub gamma: f64,
    /// Prandtl number.
    pub prandtl: f64,
}

impl WallGas {
    /// Combine the field's original thermodynamics at `(T [K], p [Pa])` with
    /// viscosity [Pa s] and thermal conductivity [W/(m K)] at that same state.
    ///
    /// The source supplies density, gamma = Cp(T)/Cv(T), and the Cp(T) used in
    /// Pr = mu Cp(T)/k. The reference must lie in the field's declared caloric
    /// domain. Transport coefficients carry their caller's source authority;
    /// this constructor does not supply transport fits or change their domain.
    /// The resulting wall coefficients stay fixed during bore evolution.
    pub fn from_source(
        gas: &FrozenGas,
        temperature_k: f64,
        pressure_pa: f64,
        dynamic_viscosity_pa_s: f64,
        thermal_conductivity_w_m_k: f64,
    ) -> Result<Self, NonlinearDuctError> {
        let (density, cp, gamma) = Self::reference(gas, temperature_k, pressure_pa)?;
        let transport = FluidTransportState::try_new(
            density,
            cp,
            thermal_conductivity_w_m_k,
            dynamic_viscosity_pa_s,
        )
        .map_err(|_| unsupported("wall viscosity, conductivity or derived transport is invalid"))?;
        Ok(Self {
            density_kg_m3: density,
            dynamic_viscosity_pa_s,
            gamma,
            prandtl: transport.prandtl(),
        })
    }

    fn reference(
        gas: &FrozenGas,
        temperature_k: f64,
        pressure_pa: f64,
    ) -> Result<(f64, f64, f64), NonlinearDuctError> {
        let domain = gas.domain();
        if !(domain.minimum_temperature..=domain.maximum_temperature).contains(&temperature_k) {
            return Err(GasError::CaloricPath(
                "wall reference temperature outside declared interval",
            )
            .into());
        }
        // Evaluate the explicit reference directly, without perturbing its
        // temperature through a conservative-energy inversion.
        let properties = gas
            .source()
            .evaluate_properties(Temperature::new(temperature_k), Pressure::new(pressure_pa))
            .map_err(GasError::from)?;
        let density = pressure_pa / (gas.gas_constant() * temperature_k);
        let cp = properties.mass_specific_cp().value();
        let gamma = cp / properties.mass_specific_cv().value();
        if ![density, cp, gamma].into_iter().all(positive) {
            return Err(unsupported(
                "wall reference thermodynamics is unrepresentable",
            ));
        }
        Ok((density, cp, gamma))
    }

    fn validate_at(
        self,
        gas: &FrozenGas,
        temperature_k: f64,
        pressure_pa: f64,
    ) -> Result<(), NonlinearDuctError> {
        let (density, _, gamma) = Self::reference(gas, temperature_k, pressure_pa)?;
        let close = |actual: f64, expected: f64| {
            positive(actual) && (actual - expected).abs() <= 2.0e-12 * actual.max(expected)
        };
        if !close(self.density_kg_m3, density) {
            return Err(unsupported(
                "wall density disagrees with the source at the room state",
            ));
        }
        if !close(self.gamma, gamma) {
            return Err(unsupported(
                "wall gamma disagrees with source Cp/Cv at the room state",
            ));
        }
        Ok(())
    }
}

/// Discretization of a sounding bore.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BoreSpec {
    /// Largest axial cell [m]. The bore is cut into equal cells; a face that
    /// would fall within a quarter cell of a radius step is moved off it.
    pub maximum_cell_m: f64,
    /// Compact extra volume on the first cell (a lumped cup) [m³], or zero.
    pub inlet_volume_m3: f64,
    /// Fixed accepted step [s]. The radiating end's memory runs on this clock;
    /// a field that cannot take it refuses rather than shortening it.
    pub step_s: f64,
    /// Admitted stage CFL, in (0, 0.45].
    pub cfl: f64,
    /// Wall memory, or `None` for rigid inviscid adiabatic walls.
    pub wall: Option<WallMemory>,
}

/// One accepted bore step.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BoreStep {
    /// The shared gas owner's accepted-step report, before the wall update.
    pub report: StepReport,
    /// Mean outward mass flow through the radiating end over the step [kg/s].
    pub outlet_mass_flow_kg_s: f64,
    /// That mass flow divided by the room density: the piston's volume
    /// velocity, the source strength the exterior pressure is referenced to [m³/s].
    pub outlet_volume_flow_m3_s: f64,
    /// Radiating-end pressure solved at the step's start [Pa].
    pub outlet_pressure_pa: f64,
}

/// Wall-memory receipt: what the bound arms admit.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WallReceipt {
    /// Distinct wall radii bound (cells of one radius share one propagator).
    pub distinct_radii: usize,
    /// Largest viscous or thermal arm count in any cell.
    pub maximum_arms: usize,
    /// Worst sampled complex error of the excess wall operators.
    pub maximum_complex_error: f64,
    /// Worst sampled error of their dissipative parts.
    pub maximum_loss_error: f64,
}

/// A bore whose gas evolves nonlinearly between a caller port and a radiating end.
pub struct NonlinearBore {
    solver: Euler1d,
    gas: FrozenGas,
    room: Primitive,
    radiation: PistonRadiation,
    walls: Option<Walls>,
    shell: Option<Shell>,
    driven: Option<DrivenWall>,
    inlet_volume_m3: f64,
    step_s: f64,
    cfl: f64,
    centers_m: Vec<f64>,
    fallbacks: u64,
}

impl NonlinearBore {
    /// Mesh `duct`, fill it with room gas at rest and bind its wall memory.
    ///
    /// `room` is the stationary exterior (and initial bore) state: its density
    /// and sound speed set the piston's radiation, and its temperature is the
    /// reference of the thermal wall arms. `wall_gas` is required with
    /// `spec.wall` and is ignored without it. Its density and gamma must agree
    /// with the original source at the declared room temperature and pressure
    /// to relative roundoff (2e-12); stale snapshots refuse before wall binding.
    pub fn new(
        duct: &Duct,
        gas: FrozenGas,
        room_temperature_k: f64,
        room_pressure_pa: f64,
        wall_gas: Option<WallGas>,
        spec: BoreSpec,
    ) -> Result<Self, NonlinearDuctError> {
        Self::with_mouth_kernel(duct, gas, room_temperature_k, room_pressure_pa, wall_gas, spec, None)
    }

    /// As [`Self::new`], with the radiating end optionally the plane-mode equivalent of the flanged multimodal mouth.
    /// `None` is [`Self::new`] bit for bit. A kernel must have been built with this bore's own gas as its bore gas (`z_inf x pi a^2 = rho c` of `room`, the gas the bore
    /// is filled with, to 1e-6; the exterior of the flanged load may be another gas), on the bore's own fixed step, and be at least two taps long.
    ///
    /// # Errors
    /// As [`Self::new`], or a kernel that does not fit this bore.
    pub fn with_mouth_kernel(
        duct: &Duct,
        gas: FrozenGas,
        room_temperature_k: f64,
        room_pressure_pa: f64,
        wall_gas: Option<WallGas>,
        spec: BoreSpec,
        mouth: Option<&MouthKernel>,
    ) -> Result<Self, NonlinearDuctError> {
        if !positive(spec.step_s) || !(spec.cfl > 0.0 && spec.cfl <= 0.45) {
            return Err(unsupported(
                "bore step must be positive and CFL in (0, 0.45]",
            ));
        }
        let geometry = bore_mesh(duct, spec.maximum_cell_m, spec.inlet_volume_m3)?;
        let room_density = room_pressure_pa / (gas.gas_constant() * room_temperature_k);
        let rest = gas.conserved(room_density, 0.0, room_pressure_pa)?;
        let room = gas.primitive(rest, || false)?;
        let n = geometry.wall_radii.len();
        let mut cells = reserved(n)?;
        cells.resize(n, rest);
        let mouth_radius = duct
            .segments
            .last()
            .map(Segment::outlet_radius)
            .ok_or(DuctError::EmptyDuct)?;
        let radiation = match mouth {
            None => PistonRadiation::new(mouth_radius, room, spec.step_s)?,
            Some(kernel) => PistonRadiation::from_kernel(mouth_radius, room, spec.step_s, kernel)?,
        };
        let walls = match (spec.wall, wall_gas) {
            (None, _) => None,
            (Some(_), None) => return Err(unsupported("wall memory needs the wall gas transport")),
            (Some(memory), Some(transport)) => {
                transport.validate_at(&gas, room_temperature_k, room_pressure_pa)?;
                Some(Walls::new(
                    &geometry.wall_radii,
                    memory,
                    transport,
                    room_temperature_k,
                    spec.step_s,
                )?)
            }
        };
        let solver = Euler1d::new_on_stepped_mesh(
            gas.clone(),
            geometry.mesh,
            &geometry.steps,
            cells,
            Boundaries::Ends {
                left: Boundary::Reflecting,
                right: Boundary::Reflecting,
            },
            || false,
        )?;
        Ok(Self {
            solver,
            gas,
            room,
            radiation,
            walls,
            shell: None,
            driven: None,
            inlet_volume_m3: spec.inlet_volume_m3,
            step_s: spec.step_s,
            cfl: spec.cfl,
            centers_m: geometry.centers_m,
            fallbacks: 0,
        })
    }

    /// Accept one fixed step. `inlet` supplies the driven end's flux at each
    /// RK stage (`None`: a closed end); the radiating end is this bore's own.
    /// The wall memory is applied after the gas step, over the same duration.
    /// A refused gas step (including one whose CFL bound is shorter than the
    /// fixed step, which the radiating end refuses) leaves the field and both
    /// memories unchanged. A wall refusal after an accepted gas step leaves
    /// the bore mid-step and must end the run.
    ///
    /// A gas step whose second-order (MUSCL) update leaves the source domain
    /// at either RK stage is redone once, from the unchanged field, at the same
    /// fixed step with first-order reconstruction: the a posteriori fallback of
    /// the MOOD scheme (Clain, Diot and Loubère, J. Comput. Phys. 230 (2011)
    /// 4028-4050). The same `admissible_state` test decides the retried step, no
    /// bound or domain moves, and a step the retry also fails is refused as
    /// [`NonlinearDuctError::GasAfterFirstOrderRetry`]. Any other refusal (a
    /// CFL bound shorter than the fixed step, a port or radiating-end
    /// refusal) is not retried. [`Self::first_order_fallbacks`] counts the
    /// steps taken this way.
    pub fn step(
        &mut self,
        mut inlet: impl FnMut(&PortStage) -> Result<Option<FaceFlux>, GasError>,
    ) -> Result<BoreStep, NonlinearDuctError> {
        let gas = &self.gas;
        let (step_s, cfl) = (self.step_s, self.cfl);
        let mut attempt = |solver: &mut Euler1d,
                           radiation: &mut PistonRadiation,
                           reconstruction: Reconstruction| {
            solver.step_with_ports(
                Time::new(step_s),
                StepControl {
                    cfl,
                    reconstruction,
                    maximum_halvings: 0,
                },
                |stage| Ok([inlet(stage)?, Some(radiation.flux(gas, stage)?)]),
                || false,
            )
        };
        let guess = self.radiation.guess;
        let report = match attempt(&mut self.solver, &mut self.radiation, Reconstruction::Muscl) {
            Ok(report) => report,
            Err(first) if source_domain_refusal(&first) => {
                // The refused attempt left the field and the end memory unchanged; only the end's
                // uncommitted predictor flow and its Newton starting pressure belong to the attempt.
                self.radiation.pending = None;
                self.radiation.guess = guess;
                let mut report = attempt(
                    &mut self.solver,
                    &mut self.radiation,
                    Reconstruction::FirstOrder,
                )
                .map_err(NonlinearDuctError::GasAfterFirstOrderRetry)?;
                // The solver sets this flag only when it retried a MUSCL control itself; this retry
                // passed first-order control, so the flag is set here with the same meaning.
                report.first_order_retry = true;
                self.fallbacks += 1;
                report
            }
            Err(error) => return Err(error.into()),
        };
        let outlet_pressure_pa = self.radiation.commit()?;
        if let Some(walls) = self.walls.as_mut() {
            walls.apply(&mut self.solver, &self.gas)?;
        }
        if let Some(shell) = self.shell.as_mut() {
            shell.apply(&mut self.solver, self.room.pressure)?;
        }
        if let Some(driven) = self.driven.as_mut() {
            driven.apply(&mut self.solver, &self.gas, self.room.pressure)?;
        }
        let outlet_mass_flow_kg_s = -report.right_boundary_change.mass / self.step_s;
        Ok(BoreStep {
            report,
            outlet_mass_flow_kg_s,
            outlet_volume_flow_m3_s: outlet_mass_flow_kg_s / self.room.density,
            outlet_pressure_pa,
        })
    }

    /// Gas steps accepted so far only after the first-order retry of
    /// [`Self::step`] (the MUSCL step failed the source-domain admission).
    #[must_use]
    pub const fn first_order_fallbacks(&self) -> u64 {
        self.fallbacks
    }

    /// The shared gas field, its mesh, source and time.
    #[must_use]
    pub const fn solver(&self) -> &Euler1d {
        &self.solver
    }

    /// The source-owned gas law (the same one the field evolves).
    #[must_use]
    pub const fn gas(&self) -> &FrozenGas {
        &self.gas
    }

    /// The stationary room state outside the radiating end.
    #[must_use]
    pub const fn room(&self) -> Primitive {
        self.room
    }

    /// Cell midpoints [m].
    #[must_use]
    pub fn centers_m(&self) -> &[f64] {
        &self.centers_m
    }

    /// Driven-end face area [m²].
    #[must_use]
    pub fn inlet_area_m2(&self) -> f64 {
        self.solver.mesh().face_areas()[0]
    }

    /// Number of past end flows the radiation memory convolves.
    #[must_use]
    pub fn radiation_taps(&self) -> usize {
        self.radiation.weights.len()
    }

    /// Product-integration weights of the regular piston kernel (or of the mouth kernel the bore was built with), newest
    /// first (index 0 multiplies the current end flow). They sum to one (the mouth kernel's within its dropped tail): the
    /// radiation impedance vanishes at zero frequency.
    #[must_use]
    pub fn radiation_weights(&self) -> &[f64] {
        &self.radiation.weights
    }

    /// What the wall memory admitted, or `None` for rigid inviscid walls.
    #[must_use]
    pub fn wall_receipt(&self) -> Option<WallReceipt> {
        self.walls.as_ref().map(|walls| walls.receipt)
    }

    /// Bind a breathing wall to every cell, before the first step. `duct` must be the duct this bore was built
    /// from (its length and every cell volume are checked), `law` supplies each cell's strip mass and stiffness
    /// at the cell's mean radius, and `quality` is the quality factor of every cell's breathing mode. The walls
    /// start at rest at the room pressure. Not binding a shell leaves every byte of the bore as it was.
    ///
    /// # Errors
    /// A bore that already stepped or already has a shell, a duct that is not this bore's own, a non-positive
    /// quality factor, or a strip the law refuses or returns non-positive.
    pub fn bind_shell(
        &mut self,
        duct: &Duct,
        law: &dyn ShellLaw,
        quality: f64,
    ) -> Result<ShellReceipt, NonlinearDuctError> {
        if self.shell.is_some() || self.driven.is_some() {
            return Err(unsupported("a bore binds its shell once"));
        }
        if self.solver.time().value() != 0.0 {
            return Err(unsupported("a shell binds before the first step"));
        }
        let mesh = self.solver.mesh();
        let (shell, receipt) = Shell::new(
            duct,
            mesh.positions(),
            mesh.volumes(),
            law,
            quality,
            self.step_s,
            self.inlet_volume_m3,
        )?;
        self.shell = Some(shell);
        Ok(receipt)
    }

    /// What binding the breathing wall admitted, or `None` without a shell.
    #[must_use]
    pub fn shell_receipt(&self) -> Option<ShellReceipt> {
        self.shell.as_ref().map(|shell| shell.receipt)
    }

    /// The breathing wall's energy account so far, or `None` without a shell or a wall driver.
    #[must_use]
    pub fn shell_ledger(&self) -> Option<ShellLedger> {
        self.shell
            .as_ref()
            .map(Shell::ledger)
            .or_else(|| self.driven.as_ref().map(DrivenWall::ledger))
    }

    /// Bind a wall that couples the cells (a modal shell, a [`WallDriver`]), before the first step. It covers every
    /// cell of this bore; the bore keeps the gas side (the rescale for the volume swept, the work, the account).
    /// A bore binds one wall, a ring shell ([`Self::bind_shell`]) or a driver.
    ///
    /// # Errors
    /// A bore that already stepped or already has a wall, or a driver that covers a different number of cells.
    pub fn bind_wall_driver(&mut self, driver: Box<dyn WallDriver>) -> Result<(), NonlinearDuctError> {
        if self.shell.is_some() || self.driven.is_some() {
            return Err(unsupported("a bore binds one breathing wall"));
        }
        if self.solver.time().value() != 0.0 {
            return Err(unsupported("a wall driver binds before the first step"));
        }
        let nominal = self.solver.mesh().volumes().to_vec();
        self.driven = Some(DrivenWall::new(driver, nominal)?);
        Ok(())
    }
}

// The gas owner's two refusals of an RK stage that leaves the source domain (fs-gas dynamics/evolve.rs), as
// it reports them once its retries are spent. A test provokes the real refusal, so a changed text fails there.
fn source_domain_refusal(error: &GasError) -> bool {
    matches!(
        error,
        GasError::StepRetries { last, .. }
            if matches!(
                **last,
                GasError::Invalid(
                    "first RK stage outside source domain" | "second RK stage outside source domain"
                )
            )
    )
}

struct BoreMesh {
    mesh: Mesh1d,
    wall_radii: Vec<f64>,
    centers_m: Vec<f64>,
    /// `(face, area seen by the cell on the right)` of every abrupt radius step.
    steps: Vec<(usize, f64)>,
}

// Equal cells over each stretch of bore between radius steps (one stretch, and the
// same equal cells as always, when there is no step), a junction face exactly at every
// radius step so no cell contains one, exact piecewise volumes, face areas from the
// original geometry (at a step: the outlet area of the segment before it; the inlet
// area of the segment after it travels in `steps`).
fn bore_mesh(
    duct: &Duct,
    maximum_cell_m: f64,
    inlet_volume_m3: f64,
) -> Result<BoreMesh, NonlinearDuctError> {
    if duct.segments.is_empty() {
        return Err(DuctError::EmptyDuct.into());
    }
    if !positive(maximum_cell_m) || !(inlet_volume_m3 >= 0.0 && inlet_volume_m3.is_finite()) {
        return Err(unsupported("bore cell size and inlet volume"));
    }
    let mut starts = reserved(duct.segments.len() + 1)?;
    starts.push(0.0);
    for segment in &duct.segments {
        segment.validate()?;
        if matches!(segment, Segment::ToneHole { .. }) {
            return Err(unsupported(
                "tone holes require a conservative gas junction",
            ));
        }
        let end = starts.last().copied().unwrap_or(0.0) + segment.axial_length();
        if !positive(end) {
            return Err(unsupported("bore length is not representable"));
        }
        starts.push(end);
    }
    let length = starts[duct.segments.len()];
    // Radius steps: the segments whose inlet radius differs from the outlet radius before them.
    let mut step_segments: Vec<usize> = reserved(duct.segments.len())?;
    for k in 1..duct.segments.len() {
        let before = duct.segments[k - 1].outlet_radius();
        let after = duct.segments[k].inlet_radius();
        if (before - after).abs() > 64.0 * f64::EPSILON * before.max(after) {
            step_segments.push(k);
        }
    }
    // Equal cells over each stretch between steps, none longer than `maximum_cell_m`.
    let stretches = step_segments.len() + 1;
    let bounds = |j: usize| -> (f64, f64) {
        (
            if j == 0 { 0.0 } else { starts[step_segments[j - 1]] },
            if j + 1 == stretches { length } else { starts[step_segments[j]] },
        )
    };
    let mut counts: Vec<f64> = reserved(stretches)?;
    let mut cells = 0.0_f64;
    for j in 0..stretches {
        let (a, b) = bounds(j);
        let count = ((b - a) / maximum_cell_m).ceil().max(1.0);
        counts.push(count);
        cells += count;
    }
    if !(1.0..=super::MAX_NONLINEAR_DUCT_CELLS as f64).contains(&cells) {
        return Err(unsupported("cell count exceeds MAX_NONLINEAR_DUCT_CELLS"));
    }
    let n = cells as usize;
    let mut faces = reserved(n + 1)?;
    let mut steps: Vec<(usize, f64)> = reserved(step_segments.len())?;
    for j in 0..stretches {
        let (a, b) = bounds(j);
        let m = counts[j] as usize;
        if j > 0 {
            let k = step_segments[j - 1];
            steps.push((faces.len(), area(duct.segments[k].inlet_radius())?));
        }
        for i in 0..(if j + 1 == stretches { m + 1 } else { m }) {
            faces.push(if i == m {
                b
            } else {
                a + (b - a) * (i as f64 / counts[j])
            });
        }
    }
    if stretches > 1
        && faces
            .windows(2)
            .any(|w| !(w[1] - w[0] >= 0.4 * maximum_cell_m))
    {
        return Err(unsupported("a stretch between radius steps is shorter than 0.4 cell"));
    }
    let locate = |x: f64| -> usize {
        let k = starts.partition_point(|&s| s <= x);
        k.saturating_sub(1).min(duct.segments.len() - 1)
    };
    let mut areas = reserved(n + 1)?;
    for (i, &x) in faces.iter().enumerate() {
        if let Some(j) = steps.iter().position(|&(face, _)| face == i) {
            // The cell on the left of a step sees the outlet of the segment before it.
            areas.push(area(duct.segments[step_segments[j] - 1].outlet_radius())?);
            continue;
        }
        let k = locate(x);
        areas.push(area(radius(
            &duct.segments[k],
            (x - starts[k]).clamp(0.0, duct.segments[k].axial_length()),
        )?)?);
    }
    let mut volumes = reserved(n)?;
    let mut wall_radii = reserved(n)?;
    let mut centers_m = reserved(n)?;
    for i in 0..n {
        let (left, right) = (faces[i], faces[i + 1]);
        let mut volume = 0.0;
        for (k, segment) in duct.segments.iter().enumerate() {
            let a = left.max(starts[k]) - starts[k];
            let b = right.min(starts[k + 1]) - starts[k];
            if b <= a {
                continue;
            }
            let b = b.min(segment.axial_length());
            volume += match segment {
                Segment::Cylinder { radius, .. } => area(*radius)? * (b - a),
                Segment::Cone { .. } => {
                    let (ra, rb) = (radius(segment, a)?, radius(segment, b)?);
                    core::f64::consts::PI * ((b - a) / 3.0) * (ra * ra + ra * rb + rb * rb)
                }
                Segment::Profile { profile, .. } => profile_volume(profile, a, b, &mut || false)?,
                Segment::ToneHole { .. } => unreachable!("topology admitted above"),
            };
        }
        if !positive(volume) {
            return Err(unsupported("cell volume is not positive and finite"));
        }
        // The wall memory sees the cylinder of the same volume and width.
        wall_radii.push(det::sqrt(volume / (core::f64::consts::PI * (right - left))));
        centers_m.push(left.midpoint(right));
        volumes.push(if i == 0 {
            volume + inlet_volume_m3
        } else {
            volume
        });
    }
    Ok(BoreMesh {
        mesh: Mesh1d::new(faces, areas, volumes)?,
        wall_radii,
        centers_m,
        steps,
    })
}

// Share the normalized impulse kernel with the other pressure/flow owners.
fn piston_weights(step_s: f64, transit_s: f64) -> Result<Vec<f64>, NonlinearDuctError> {
    fs_math::special::semicircle_hat_weights(step_s / transit_s, 1_048_576).map_err(unsupported)
}

struct PistonRadiation {
    weights: Vec<f64>,
    history: Vec<f64>,
    newest: usize,
    pending: Option<(f64, f64)>,
    // Read by the piston's own impedance test only; `flux` takes `flow_gain`.
    #[allow(dead_code)]
    rho_c: f64,
    impedance: f64,
    // Pressure per unit end velocity at the newest flow: `impedance x area`. For the piston it is `rho_c` itself (stored, not
    // recomputed, so the piston's bits do not move); for a mouth kernel it is `z_inf x area` at the bore's own outlet area.
    flow_gain: f64,
    area: f64,
    room: Primitive,
    step_s: f64,
    guess: f64,
}

impl PistonRadiation {
    fn new(radius_m: f64, room: Primitive, step_s: f64) -> Result<Self, NonlinearDuctError> {
        let area = area(radius_m)?;
        let transit = 2.0 * radius_m / room.sound_speed;
        let weights = piston_weights(step_s, transit)?;
        let mut history = reserved(weights.len())?;
        history.resize(weights.len(), 0.0);
        let rho_c = room.density * room.sound_speed;
        Ok(Self {
            weights,
            history,
            newest: 0,
            pending: None,
            rho_c,
            impedance: rho_c / area,
            flow_gain: rho_c,
            area,
            room,
            step_s,
            guess: room.pressure,
        })
    }

    // The causal real weights of the plane-mode equivalent of the flanged multimodal mouth, in the piston's own law.
    fn from_kernel(
        radius_m: f64,
        room: Primitive,
        step_s: f64,
        kernel: &MouthKernel,
    ) -> Result<Self, NonlinearDuctError> {
        let area = area(radius_m)?;
        if !(kernel.step_s.is_finite() && (kernel.step_s - step_s).abs() <= 1e-12 * step_s) {
            return Err(unsupported("mouth kernel step differs from the bore's fixed step"));
        }
        if kernel.weights.len() < 2 || !kernel.weights.iter().all(|w| w.is_finite()) {
            return Err(unsupported("mouth kernel needs at least two finite weights"));
        }
        let rho_c = room.density * room.sound_speed;
        let kernel_rho_c =
            kernel.z_inf * core::f64::consts::PI * kernel.mouth_radius_m * kernel.mouth_radius_m;
        if !(positive(kernel.z_inf) && (kernel_rho_c - rho_c).abs() <= 1e-6 * rho_c) {
            return Err(unsupported("mouth kernel was not built for this bore's gas"));
        }
        let weights = kernel.weights.clone();
        let mut history = reserved(weights.len())?;
        history.resize(weights.len(), 0.0);
        Ok(Self {
            weights,
            history,
            newest: 0,
            pending: None,
            rho_c,
            impedance: kernel.z_inf,
            flow_gain: kernel.z_inf * area,
            area,
            room,
            step_s,
            guess: room.pressure,
        })
    }

    // Past end flows U_{t-j}, j >= 1, convolved with the kernel weights. The
    // corrector's newest past flow is the predictor's uncommitted one.
    fn memory(&self, stage: usize) -> f64 {
        let m = self.history.len();
        let mut sum = 0.0;
        let mut first = 1;
        if stage == 1 {
            sum += self.weights[1] * self.pending.map_or(0.0, |(flow, _)| flow);
            first = 2;
        }
        for j in first..self.weights.len() {
            let back = j - first;
            sum += self.weights[j] * self.history[(self.newest + m - back) % m];
        }
        sum
    }

    fn flux(&mut self, gas: &FrozenGas, stage: &PortStage) -> Result<FaceFlux, GasError> {
        if stage.dt != self.step_s
            || stage.stage > 1
            || (stage.stage == 1 && self.pending.is_none())
        {
            return Err(GasError::Invalid(
                "radiating end runs on its declared fixed step",
            ));
        }
        let interior = stage.right.1;
        let held = self.impedance * self.memory(stage.stage);
        let gain = self.flow_gain * (1.0 - self.weights[0]);
        // G(p) = p - p_room - gain u(p) + held is strictly increasing in p.
        // A previous step's pressure is only a seed: a changed interior trace
        // may put that seed outside the source's caloric path. The unchanged
        // interior state always supplies an admissible starting point.
        let mut pressure = self.guess;
        let mut wave = match gas.end_wave(interior, 1.0, pressure, || false) {
            Ok(wave) => wave,
            Err(GasError::CaloricPath(_)) => {
                pressure = interior.pressure;
                gas.end_wave(interior, 1.0, pressure, || false)?
            }
            Err(error) => return Err(error),
        };
        for _ in 0..64 {
            let residual = pressure - self.room.pressure - gain * wave.primitive.velocity + held;
            if residual.abs() <= 1e-12 * self.room.pressure {
                break;
            }
            let slope = 1.0 - gain * wave.velocity_per_pressure;
            let mut increment = -residual / slope;
            let mut accepted = None;
            let mut source_refusal = None;
            // Keep the normal Newton path to one source evaluation. Large
            // pressure changes can otherwise propose negative absolute pressure,
            // or a wave beyond the source interval despite an admissible root.
            // Backtracking changes the numerical iterate, never the physical
            // pressure law or the declared source temperature bounds.
            for _ in 0..32 {
                let next = pressure + increment;
                if next == pressure {
                    break;
                }
                if positive(next) {
                    match gas.end_wave(interior, 1.0, next, || false) {
                        Ok(candidate) => {
                            let next_residual =
                                next - self.room.pressure - gain * candidate.primitive.velocity
                                    + held;
                            if next_residual.abs() < residual.abs() {
                                accepted = Some((next, candidate));
                                break;
                            }
                        }
                        Err(error @ GasError::CaloricPath(_)) => source_refusal = Some(error),
                        Err(error) => return Err(error),
                    }
                }
                increment *= 0.5;
            }
            let Some((next, candidate)) = accepted else {
                // Retain the original source refusal when the physical root
                // cannot be reached inside its supported caloric path.
                return Err(source_refusal.unwrap_or(GasError::Invalid(
                    "radiating end characteristic did not converge",
                )));
            };
            let converged = (next - pressure).abs() <= 1e-13 * pressure;
            pressure = next;
            wave = candidate;
            if converged {
                break;
            }
        }
        let velocity = wave.primitive.velocity;
        let residual = pressure - self.room.pressure - gain * velocity + held;
        if !(residual.abs() <= 1e-9 * self.room.pressure) {
            return Err(GasError::Invalid(
                "radiating end characteristic did not converge",
            ));
        }
        // Outflow carries the interior's entropy. Inflow carries the room's: the
        // room brought to the end pressure along this source's wave curve
        // (isentropic in expansion, its Hugoniot in compression).
        let face = if velocity >= 0.0 {
            wave.primitive
        } else {
            let room_side = gas.end_wave(self.room, 1.0, pressure, || false)?;
            Primitive {
                velocity,
                ..room_side.primitive
            }
        };
        if stage.stage == 0 {
            self.pending = Some((self.area * velocity, pressure));
        }
        self.guess = pressure;
        gas.face_flux(face, wave.signal_speed)
    }

    fn commit(&mut self) -> Result<f64, NonlinearDuctError> {
        let (flow, pressure) = self
            .pending
            .take()
            .ok_or_else(|| unsupported("radiating end has no accepted predictor flow"))?;
        let m = self.history.len();
        self.newest = (self.newest + 1) % m;
        self.history[self.newest] = flow;
        self.guess = pressure;
        Ok(pressure)
    }
}

// One cell's exactly propagated wall arms, held in the eigen-coordinates of
// their energy-scaled generator, so a step costs a few passes over the arms.
struct Propagator {
    size: usize,
    // G = Q diag(lambda) Q^T: each eigen-coordinate decays by exp(lambda dt);
    // `head` is Q's first row, the eigen-components of the mean coordinate.
    decay: Vec<f64>,
    head: Vec<f64>,
    #[cfg(test)]
    basis: Vec<f64>,
}

impl Propagator {
    // G = [[-(d + sum g), sqrt(g p)^T], [sqrt(g p), -diag(p)]] is symmetric and
    // negative semidefinite (x^T G x = -d x0^2 - sum (sqrt(g) x0 - sqrt(p) xk)^2),
    // so exp(G dt) = Q diag(exp(lambda dt)) Q^T is a contraction.
    fn new(drag: f64, arms: &[(f64, f64)], dt: f64) -> Result<Self, NonlinearDuctError> {
        let size = arms.len() + 1;
        let mut g = vec![0.0; size * size];
        g[0] = -(drag + arms.iter().map(|&(gain, _)| gain).sum::<f64>());
        for (k, &(gain, rate)) in arms.iter().enumerate() {
            let coupling = det::sqrt(gain * rate);
            g[k + 1] = coupling;
            g[(k + 1) * size] = coupling;
            g[(k + 1) * size + k + 1] = -rate;
        }
        let scale = g.iter().map(|x| x.abs()).fold(0.0, f64::max);
        if !scale.is_finite() || !(dt.is_finite() && dt > 0.0) {
            return Err(unsupported("wall arm generator is not finite"));
        }
        let (values, vectors) = fs_la::eigen::jacobi_eigh(&g, size);
        for (k, &lambda) in values.iter().enumerate() {
            let mut residual = 0.0_f64;
            let mut norm = 0.0;
            for i in 0..size {
                let product: f64 = (0..size)
                    .map(|j| g[i * size + j] * vectors[j * size + k])
                    .sum();
                residual = residual.max((product - lambda * vectors[i * size + k]).abs());
                norm += vectors[i * size + k] * vectors[i * size + k];
            }
            if !lambda.is_finite()
                || lambda > 1e-12 * scale
                || residual > 1e-12 * scale
                || (norm - 1.0).abs() > 1e-12
            {
                return Err(unsupported(
                    "wall arm generator is not resolved as negative semidefinite",
                ));
            }
        }
        Ok(Self {
            size,
            // The exact generator is semidefinite; rounding may leave its
            // (thermal) zero mode a hair above zero, whose decay is one.
            decay: values.iter().map(|&l| det::exp(l.min(0.0) * dt)).collect(),
            head: vectors[..size].to_vec(),
            #[cfg(test)]
            basis: vectors,
        })
    }

    // Set the mean coordinate to `mean` (the arms keep their state), advance
    // one step in place, and return the new mean with |x|^2 before and after.
    fn advance(&self, state: &mut [f64], mean: f64) -> (f64, f64, f64) {
        let current: f64 = self.head.iter().zip(&*state).map(|(h, y)| h * y).sum();
        let shift = mean - current;
        let (mut before, mut after, mut next) = (0.0, 0.0, 0.0);
        for ((y, &h), &d) in state.iter_mut().zip(&self.head).zip(&self.decay) {
            let shifted = *y + shift * h;
            before += shifted * shifted;
            *y = shifted * d;
            after += *y * *y;
            next += h * *y;
        }
        (next, before, after)
    }

    // The dense exp(G dt), for the tests.
    #[cfg(test)]
    fn matrix(&self) -> Vec<f64> {
        let n = self.size;
        let mut out = vec![0.0; n * n];
        for r in 0..n {
            for c in 0..n {
                out[r * n + c] = (0..n)
                    .map(|k| self.basis[r * n + k] * self.decay[k] * self.basis[c * n + k])
                    .sum();
            }
        }
        out
    }
}

struct Walls {
    // Per cell: indices of its viscous and thermal propagators, and arm states.
    viscous_of: Vec<usize>,
    thermal_of: Vec<usize>,
    propagators: Vec<Propagator>,
    viscous_state: Vec<f64>,
    thermal_state: Vec<f64>,
    spare_viscous: Vec<f64>,
    spare_thermal: Vec<f64>,
    state_stride: usize,
    wall_temperature: f64,
    receipt: WallReceipt,
}

impl Walls {
    fn new(
        radii: &[f64],
        memory: WallMemory,
        gas: WallGas,
        wall_temperature: f64,
        dt: f64,
    ) -> Result<Self, NonlinearDuctError> {
        let basis = ZkDiffusiveBasis::new(memory.poles)
            .map_err(|_| unsupported("wall memory pole count"))?;
        let band = [
            core::f64::consts::TAU * memory.band_hz[0],
            core::f64::consts::TAU * memory.band_hz[1],
        ];
        let mut bound: Vec<(u64, usize, usize)> = Vec::new();
        let mut propagators = Vec::new();
        let mut viscous_of = reserved(radii.len())?;
        let mut thermal_of = reserved(radii.len())?;
        let mut receipt = WallReceipt {
            distinct_radii: 0,
            maximum_arms: 0,
            maximum_complex_error: 0.0,
            maximum_loss_error: 0.0,
        };
        for &radius in radii {
            let key = radius.to_bits();
            let (v, t) = if let Some(&(_, v, t)) = bound.iter().find(|(k, _, _)| *k == key) {
                (v, t)
            } else {
                let zk = basis
                    .bind(
                        radius,
                        gas.density_kg_m3,
                        gas.dynamic_viscosity_pa_s,
                        gas.gamma,
                        gas.prandtl,
                        band,
                        memory.maximum_relative_error,
                    )
                    .map_err(|_| {
                        unsupported(
                            "wall memory does not resolve this radius on the band: add poles, narrow the band or relax the error",
                        )
                    })?;
                let arms = |terms: &[fs_phs::ZkDiffusiveTerm]| {
                    terms
                        .iter()
                        .map(|t| (t.gain_per_s, t.rate_per_s))
                        .collect::<Vec<_>>()
                };
                let viscous = arms(zk.series_terms());
                // A nonlinear caloric cell owns both density and temperature.
                // Its heat source must vanish at T_wall, for every density.
                // With theta = p - rho R T_wall, the circular-wall kernel is
                // H(s) = gamma s F_t/(1-F_t), rather than the pressure-only
                // acoustic shunt (gamma-1) s F_t. Indeed,
                // s p + H theta = gamma R T_wall s rho gives the same ZK
                // compliance [1+(gamma-1) F_t]/(gamma p0).
                // Reuse the admitted positive series spectrum at thermal
                // diffusivity nu/Pr, with gamma multiplying its drag/gains.
                let thermal_model = basis
                    .bind(
                        radius,
                        gas.density_kg_m3,
                        gas.dynamic_viscosity_pa_s / gas.prandtl,
                        1.0,
                        1.0,
                        band,
                        memory.maximum_relative_error,
                    )
                    .map_err(|_| {
                        unsupported(
                            "thermal wall memory does not resolve this radius on the band: add poles, narrow the band or relax the error",
                        )
                    })?;
                let thermal: Vec<_> = thermal_model
                    .series_terms()
                    .iter()
                    .map(|term| (gas.gamma * term.gain_per_s, term.rate_per_s))
                    .collect();
                receipt.maximum_arms = receipt.maximum_arms.max(viscous.len()).max(thermal.len());
                receipt.maximum_complex_error = receipt
                    .maximum_complex_error
                    .max(zk.max_complex_error())
                    .max(thermal_model.max_complex_error());
                receipt.maximum_loss_error = receipt
                    .maximum_loss_error
                    .max(zk.max_real_loss_error())
                    .max(thermal_model.max_real_loss_error());
                propagators.push(Propagator::new(zk.series_drag_per_s(), &viscous, dt)?);
                propagators.push(Propagator::new(
                    gas.gamma * thermal_model.series_drag_per_s(),
                    &thermal,
                    dt,
                )?);
                let entry = (key, propagators.len() - 2, propagators.len() - 1);
                bound.push(entry);
                (entry.1, entry.2)
            };
            viscous_of.push(v);
            thermal_of.push(t);
        }
        receipt.distinct_radii = bound.len();
        let stride = receipt.maximum_arms + 1;
        Ok(Self {
            viscous_of,
            thermal_of,
            propagators,
            viscous_state: vec![0.0; radii.len() * stride],
            thermal_state: vec![0.0; radii.len() * stride],
            spare_viscous: vec![0.0; radii.len() * stride],
            spare_thermal: vec![0.0; radii.len() * stride],
            state_stride: stride,
            wall_temperature,
            receipt,
        })
    }

    fn apply(&mut self, solver: &mut Euler1d, gas: &FrozenGas) -> Result<(), NonlinearDuctError> {
        let stride = self.state_stride;
        let reference = self.wall_temperature;
        let Self {
            propagators,
            viscous_of,
            thermal_of,
            viscous_state,
            thermal_state,
            spare_viscous,
            spare_thermal,
            ..
        } = self;
        // Every cell writes its whole next arm state into the spare buffers,
        // which are published only with the accepted cells.
        solver.apply_cell_sources(
            |i, q, p| {
                let at = i * stride;
                let viscous = &propagators[viscous_of[i]];
                let n = viscous.size;
                let arms = &mut spare_viscous[at..at + n];
                arms.copy_from_slice(&viscous_state[at..at + n]);
                let (velocity, before, after) = viscous.advance(arms, p.velocity);
                // The arms hold the boundary layer's extra kinetic energy,
                // rho |z|^2 / 2; only what the propagator removes from the
                // mean flow plus arms is dissipated, and it heats the gas.
                let dissipated = 0.5 * q.density * (before - after);
                let kinetic_change =
                    0.5 * q.density * (velocity * velocity - p.velocity * p.velocity);
                let moved = Conserved {
                    density: q.density,
                    momentum: q.density * velocity,
                    energy: q.energy + (kinetic_change + dissipated),
                };
                // Finite heating follows the original integral of Cv(T),
                // not a locally frozen gamma. At zero dissipation the old
                // temperature/pressure remain exact and need no inversion.
                let heated = if dissipated == 0.0 {
                    p.temperature
                } else {
                    gas.primitive(moved, || false)?.temperature
                };
                let thermal = &propagators[thermal_of[i]];
                let n = thermal.size;
                let arms = &mut spare_thermal[at..at + n];
                arms.copy_from_slice(&thermal_state[at..at + n]);
                let (exchanged, _, _) = thermal.advance(arms, heated - reference);
                // No wall heat exchange must retain the conservative
                // viscous work exactly, without a caloric round trip.
                if exchanged == heated - reference {
                    return Ok(moved);
                }
                // Store the thermal history in temperature coordinates so
                // transported density cannot reinterpret it as a heat pulse.
                // Recover the energy at the CURRENT density through the same
                // NASA-9 owner, including its original temperature limits.
                let pressure = q.density * gas.gas_constant() * (exchanged + reference);
                gas.conserved(q.density, velocity, pressure)
            },
            || false,
        )?;
        core::mem::swap(viscous_state, spare_viscous);
        core::mem::swap(thermal_state, spare_thermal);
        Ok(())
    }
}

#[cfg(test)]
#[path = "bore_tests.rs"]
mod tests;
