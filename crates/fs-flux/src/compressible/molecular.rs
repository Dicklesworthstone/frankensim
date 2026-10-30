//! Molecular transport for the same frozen, source-bearing Cartesian gas.
//!
//! The original `fs-gas::multicomponent::dissipation::TransportLaw` evaluates
//! viscosity, bulk viscosity and conductivity at each actual face state. The
//! complete planar Newtonian stress uses the three-dimensional `2/3 div(u)`
//! trace convention with zero out-of-plane velocity. Its work and Fourier heat
//! share one conservative total-energy flux; viscous heat is not added again.
//! Composition is spatially uniform and frozen: this is not a species-diffusion
//! or reacting-flow solver. A composition change requires a new physical state.

use super::*;
use fs_gas::multicomponent::dissipation::{DissipationError, TransportLaw};
use fs_qty::{Pressure, Temperature};

mod coupled;
pub use coupled::{MolecularEuler2d, MolecularStepControl, MolecularStepReport};

/// Explicit molecular substep and source-admissibility budgets.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MolecularControl {
    /// Maximum conservative parabolic estimate, in `(0, 0.45]`.
    pub parabolic_cfl: f64,
    /// Largest internal SSP RK2 substep [s].
    pub maximum_step_s: f64,
    /// Minimum retry step [s]; an exact final remainder may be smaller.
    pub minimum_step_s: f64,
    /// Maximum accepted plus rejected attempts for the complete duration.
    pub maximum_steps: usize,
    /// Maximum halvings of one substep, at most 32.
    pub maximum_halvings: usize,
}
impl Default for MolecularControl {
    fn default() -> Self {
        Self {
            parabolic_cfl: 0.4,
            maximum_step_s: 1.0,
            minimum_step_s: 1e-16,
            maximum_steps: 10000,
            maximum_halvings: 20,
        }
    }
}
impl MolecularControl {
    fn validate(self) -> Result<(), MolecularError> {
        if !positive(self.parabolic_cfl)
            || self.parabolic_cfl > 0.45
            || !positive(self.minimum_step_s)
            || !positive(self.maximum_step_s)
            || self.minimum_step_s > self.maximum_step_s
            || !(1..=1_000_000).contains(&self.maximum_steps)
            || self.maximum_halvings > 32
        {
            return Err(MolecularError::Invalid("planar molecular controls"));
        }
        Ok(())
    }
}

/// Complete requested-duration molecular update; the advection clock is unchanged.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MolecularReport {
    /// Complete represented duration.
    pub elapsed: Time,
    /// Accepted SSP RK2 substeps.
    pub accepted_steps: usize,
    /// Rejected internal attempts.
    pub rejected_steps: usize,
    /// Largest estimated parabolic number over accepted stages.
    pub maximum_parabolic_number: f64,
    /// Original extensive state.
    pub before: Totals,
    /// Accepted extensive state.
    pub after: Totals,
    /// Net molecular transfer into the domain.
    pub boundary_change: Totals,
    /// `after - before - boundary_change`.
    pub conservation_residual: Totals,
}

/// A refusal never publishes a partial field, clock, or boundary impulse.
#[derive(Clone, Debug, PartialEq)]
pub enum MolecularError {
    /// Source, geometry, control, or arithmetic mismatch.
    Invalid(&'static str),
    /// Original gas/caloric or Euler transport refusal.
    Flow(Error),
    /// Original molecular material-law refusal.
    Law(DissipationError),
    /// An SSP stage exceeds the declared parabolic estimate.
    ParabolicNumber { observed: f64, maximum: f64 },
    /// The complete duration exceeds the caller's attempt budget.
    WorkLimit {
        elapsed_s: f64,
        accepted: usize,
        rejected: usize,
    },
    /// No source-admissible internal step remains.
    StepFailure {
        elapsed_s: f64,
        step_s: f64,
        last: Box<MolecularError>,
    },
    /// A complete split step could not retain a common physical duration.
    CouplingRetries {
        attempts: usize,
        last: Box<MolecularError>,
    },
    /// Cancellation before publication.
    Cancelled,
}
impl core::fmt::Display for MolecularError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "planar molecular gas transport: {self:?}")
    }
}
impl std::error::Error for MolecularError {}
impl From<Error> for MolecularError {
    fn from(error: Error) -> Self {
        if error == Error::Cancelled {
            Self::Cancelled
        } else {
            Self::Flow(error)
        }
    }
}
impl From<DissipationError> for MolecularError {
    fn from(error: DissipationError) -> Self {
        if error == DissipationError::Cancelled {
            Self::Cancelled
        } else {
            Self::Law(error)
        }
    }
}

struct MolecularWorkspace {
    current: Vec<State>,
    stage: Vec<State>,
    candidate: Vec<State>,
    primitive: Vec<Primitive>,
    cv_density: Vec<f64>,
    x_flux: Vec<NormalFaceFlux>,
    y_flux: Vec<NormalFaceFlux>,
    binary: Vec<f64>,
    boundary_step: Vec<BoundaryImpulse>,
    boundary_trial: Vec<BoundaryImpulse>,
}

/// Prepared Fourier/Newtonian operator on the exact original gas and geometry.
///
/// Periodic faces share one flux. Reflecting boundaries are stationary adiabatic
/// slip walls: their normal viscous traction loads the support but transfers no
/// heat or work. Transmissive boundaries have zero molecular flux. They do not
/// represent a prescribed-temperature or no-slip wall.
pub struct MolecularTransport2d<L> {
    gas: FrozenGas,
    grid: CartesianGrid,
    boundaries: Boundaries2d,
    law: L,
    boundary_impulses: Vec<BoundaryImpulse>,
    work: MolecularWorkspace,
}
impl<L: TransportLaw> MolecularTransport2d<L> {
    /// Bind the same original species models; labels alone never establish identity.
    pub fn new(
        flow: &Euler2d,
        law: L,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, MolecularError> {
        poll(&mut cancelled)?;
        if law.source().models() != flow.gas.source().components() {
            return Err(MolecularError::Invalid(
                "transport and Euler caloric sources differ",
            ));
        }
        let n = flow.states.len();
        let [nx, ny] = flow.grid.shape;
        let species = law.source().models().len();
        let result = Self {
            gas: flow.gas.clone(),
            grid: flow.grid,
            boundaries: flow.boundaries,
            law,
            boundary_impulses: allocated(2 * (nx + ny), Totals::default())?,
            work: MolecularWorkspace {
                current: allocated(n, State::default())?,
                stage: allocated(n, State::default())?,
                candidate: allocated(n, State::default())?,
                primitive: allocated(n, Primitive::default())?,
                cv_density: allocated(n, 0.0)?,
                x_flux: allocated((nx + 1) * ny, NormalFaceFlux::default())?,
                y_flux: allocated(nx * (ny + 1), NormalFaceFlux::default())?,
                binary: allocated(species * species, 0.0)?,
                boundary_step: allocated(2 * (nx + ny), Totals::default())?,
                boundary_trial: allocated(2 * (nx + ny), Totals::default())?,
            },
        };
        poll(&mut cancelled)?;
        Ok(result)
    }
    /// Original immutable law identity and retained collision/source data.
    pub const fn law(&self) -> &L {
        &self.law
    }
    /// Outward molecular transfer during the last complete accepted update.
    pub fn boundary_impulses(&self, side: BoundarySide) -> &[BoundaryImpulse] {
        &self.boundary_impulses[self.grid.boundary_range(side)]
    }

    /// Advance the entire molecular interval atomically, without advancing time.
    ///
    /// The caller that composes operators owns physical time. The original
    /// `Euler2d::last_step` and its advective boundary report are left intact.
    pub fn advance(
        &mut self,
        flow: &mut Euler2d,
        duration: Time,
        controls: MolecularControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MolecularReport, MolecularError> {
        poll(&mut cancelled)?;
        controls.validate()?;
        if self.gas != flow.gas || self.grid != flow.grid || self.boundaries != flow.boundaries {
            return Err(MolecularError::Invalid(
                "molecular operator and flow differ",
            ));
        }
        if !positive(duration.value()) {
            return Err(MolecularError::Invalid("molecular duration"));
        }
        let before = checked_totals(&flow.states, self.grid.cell_volume(), &mut cancelled)?;
        self.work.current.copy_from_slice(&flow.states);
        self.work.boundary_trial.fill(Totals::default());
        let mut elapsed = 0.0;
        let mut accepted = 0;
        let mut rejected = 0;
        let mut maximum_parabolic: f64 = 0.0;
        while elapsed < duration.value() {
            poll(&mut cancelled)?;
            let first_rate = spatial(
                &self.gas,
                &self.grid,
                self.boundaries,
                &self.work.current,
                &mut self.law,
                &mut self.work.primitive,
                &mut self.work.cv_density,
                &mut self.work.x_flux,
                &mut self.work.y_flux,
                &mut self.work.binary,
                &mut cancelled,
            )?;
            let remaining = duration.value() - elapsed;
            let mut dt = remaining.min(controls.maximum_step_s);
            if first_rate > 0.0 {
                dt = dt.min(controls.parabolic_cfl / first_rate);
            }
            for halving in 0..=controls.maximum_halvings {
                poll(&mut cancelled)?;
                if accepted + rejected >= controls.maximum_steps {
                    return Err(MolecularError::WorkLimit {
                        elapsed_s: elapsed,
                        accepted,
                        rejected,
                    });
                }
                if !positive(dt)
                    || elapsed + dt <= elapsed
                    || (dt < controls.minimum_step_s && dt != remaining)
                {
                    return Err(MolecularError::Invalid("unrepresentable molecular substep"));
                }
                // A rejected second stage overwrites face work; reconstruct the
                // unchanged first stage before every retry.
                let trial = (|| {
                    let r1 = spatial(
                        &self.gas,
                        &self.grid,
                        self.boundaries,
                        &self.work.current,
                        &mut self.law,
                        &mut self.work.primitive,
                        &mut self.work.cv_density,
                        &mut self.work.x_flux,
                        &mut self.work.y_flux,
                        &mut self.work.binary,
                        &mut cancelled,
                    )?;
                    parabolic(r1 * dt, controls.parabolic_cfl)?;
                    update(
                        &self.gas,
                        &self.grid,
                        &self.work.current,
                        &self.work.current,
                        &self.work.x_flux,
                        &self.work.y_flux,
                        &mut self.work.stage,
                        dt,
                        false,
                        &mut cancelled,
                    )?;
                    capture_boundaries(
                        &self.grid,
                        &self.work.x_flux,
                        &self.work.y_flux,
                        0.5 * dt,
                        &mut self.work.boundary_step,
                        false,
                        &mut cancelled,
                    )?;
                    let r2 = spatial(
                        &self.gas,
                        &self.grid,
                        self.boundaries,
                        &self.work.stage,
                        &mut self.law,
                        &mut self.work.primitive,
                        &mut self.work.cv_density,
                        &mut self.work.x_flux,
                        &mut self.work.y_flux,
                        &mut self.work.binary,
                        &mut cancelled,
                    )?;
                    parabolic(r2 * dt, controls.parabolic_cfl)?;
                    update(
                        &self.gas,
                        &self.grid,
                        &self.work.current,
                        &self.work.stage,
                        &self.work.x_flux,
                        &self.work.y_flux,
                        &mut self.work.candidate,
                        dt,
                        true,
                        &mut cancelled,
                    )?;
                    capture_boundaries(
                        &self.grid,
                        &self.work.x_flux,
                        &self.work.y_flux,
                        0.5 * dt,
                        &mut self.work.boundary_step,
                        true,
                        &mut cancelled,
                    )?;
                    Ok::<f64, MolecularError>(r1.max(r2) * dt)
                })();
                match trial {
                    Ok(number) => {
                        maximum_parabolic = maximum_parabolic.max(number);
                        for (total, step) in self
                            .work
                            .boundary_trial
                            .iter_mut()
                            .zip(&self.work.boundary_step)
                        {
                            *total = add_totals(*total, *step);
                        }
                        std::mem::swap(&mut self.work.current, &mut self.work.candidate);
                        elapsed = if dt == remaining {
                            duration.value()
                        } else {
                            elapsed + dt
                        };
                        accepted += 1;
                        break;
                    }
                    Err(MolecularError::Cancelled) => return Err(MolecularError::Cancelled),
                    Err(error) => {
                        rejected += 1;
                        if halving == controls.maximum_halvings
                            || dt * 0.5 < controls.minimum_step_s
                        {
                            return Err(MolecularError::StepFailure {
                                elapsed_s: elapsed,
                                step_s: dt,
                                last: Box::new(error),
                            });
                        }
                        dt *= 0.5;
                    }
                }
            }
        }
        let after = checked_totals(&self.work.current, self.grid.cell_volume(), &mut cancelled)?;
        let outward = sum_arrays(
            self.work.boundary_trial.iter().map(|q| q.array()),
            &mut cancelled,
        )?;
        let boundary_change = Totals::from_array(outward.map(|v| -v));
        let conservation_residual = balance(before, after, boundary_change)?;
        poll(&mut cancelled)?;
        flow.states.copy_from_slice(&self.work.current);
        self.boundary_impulses
            .copy_from_slice(&self.work.boundary_trial);
        Ok(MolecularReport {
            elapsed: duration,
            accepted_steps: accepted,
            rejected_steps: rejected,
            maximum_parabolic_number: maximum_parabolic,
            before,
            after,
            boundary_change,
            conservation_residual,
        })
    }
}

fn parabolic(observed: f64, maximum: f64) -> Result<(), MolecularError> {
    if !observed.is_finite() || observed > maximum * (1.0 + 32.0 * f64::EPSILON) {
        Err(MolecularError::ParabolicNumber { observed, maximum })
    } else {
        Ok(())
    }
}
fn add_totals(a: Totals, b: Totals) -> Totals {
    let a = a.array();
    let b = b.array();
    Totals::from_array(std::array::from_fn(|k| a[k] + b[k]))
}
fn balance(before: Totals, after: Totals, boundary: Totals) -> Result<Totals, MolecularError> {
    let a = after.array();
    let b = before.array();
    let c = boundary.array();
    let result = Totals::from_array(std::array::from_fn(|k| a[k] - b[k] - c[k]));
    if !a
        .into_iter()
        .chain(b)
        .chain(c)
        .chain(result.array())
        .all(f64::is_finite)
    {
        return Err(MolecularError::Invalid(
            "unrepresentable molecular extensive balance",
        ));
    }
    Ok(result)
}

#[allow(clippy::too_many_arguments)]
fn spatial<L: TransportLaw>(
    gas: &FrozenGas,
    grid: &CartesianGrid,
    boundaries: Boundaries2d,
    states: &[State],
    law: &mut L,
    primitive: &mut [Primitive],
    cv_density: &mut [f64],
    x_flux: &mut [NormalFaceFlux],
    y_flux: &mut [NormalFaceFlux],
    binary: &mut [f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, MolecularError> {
    let mut minimum_density = f64::INFINITY;
    let mut minimum_cv = f64::INFINITY;
    for ((q, p), cv) in states.iter().zip(primitive.iter_mut()).zip(cv_density) {
        *p = q.primitive(gas, &mut *cancelled)?;
        let props = gas
            .source()
            .evaluate_properties(Temperature::new(p.temperature), Pressure::new(p.pressure))
            .map_err(Error::from)?;
        *cv = p.density * props.mass_specific_cv().value();
        if !positive(*cv) {
            return Err(MolecularError::Invalid(
                "molecular volumetric heat capacity",
            ));
        }
        minimum_density = minimum_density.min(p.density);
        minimum_cv = minimum_cv.min(*cv);
    }
    let [nx, _] = grid.shape;
    let mut maximum_diffusivity: f64 = 0.0;
    for axis in 0..2 {
        let count = grid.shape[axis];
        let tangent = 1 - axis;
        let pair = boundaries.axis(axis);
        for line in 0..grid.shape[tangent] {
            let index = |position: usize| {
                if axis == 0 {
                    line * nx + position
                } else {
                    position * nx + line
                }
            };
            for face in 0..=count {
                poll(cancelled)?;
                let fluxes = if axis == 0 {
                    &mut *x_flux
                } else {
                    &mut *y_flux
                };
                let f = face_index(grid, axis, face, line);
                if face == count && pair == AxisBoundaries::Periodic {
                    fluxes[f] = fluxes[face_index(grid, axis, 0, line)];
                    continue;
                }
                let exterior = face == 0 || face == count;
                let wall = if let AxisBoundaries::Ends { lower, upper } = pair {
                    exterior
                        && (if face == 0 { lower } else { upper }) == BoundaryCondition::Reflecting
                } else {
                    false
                };
                if exterior && pair != AxisBoundaries::Periodic && !wall {
                    fluxes[f] = NormalFaceFlux::default();
                    continue;
                }
                let left = index(if face == 0 {
                    if wall { 0 } else { count - 1 }
                } else {
                    face - 1
                });
                let right = index(if face == count { count - 1 } else { face });
                let pl = primitive[left];
                let pr = primitive[right];
                let t = 0.5 * pl.temperature + 0.5 * pr.temperature;
                let p = 0.5 * pl.pressure + 0.5 * pr.pressure;
                let coefficients = law.evaluate(
                    Temperature::new(t),
                    Pressure::new(p),
                    gas.source().mole_composition().fractions(),
                    binary,
                    cancelled,
                )?;
                let mu = coefficients.dynamic_viscosity.value();
                let bulk = coefficients.bulk_viscosity.value();
                let conductivity = coefficients.thermal_conductivity.value();
                if ![mu, bulk, conductivity]
                    .into_iter()
                    .all(|v| v.is_finite() && v >= 0.0)
                {
                    return Err(MolecularError::Invalid(
                        "nonnegative finite molecular coefficients required",
                    ));
                }
                let longitudinal = (4.0 / 3.0) * mu + bulk;
                let cross = bulk - (2.0 / 3.0) * mu;
                let diffusivity = ((longitudinal + cross.abs() + 2.0 * mu) / minimum_density)
                    .max(conductivity / minimum_cv);
                if !diffusivity.is_finite() {
                    return Err(MolecularError::Invalid(
                        "molecular diffusion scale overflow",
                    ));
                }
                maximum_diffusivity = maximum_diffusivity.max(diffusivity);
                let tangent_derivative = |cell: usize, component: usize| {
                    let i = cell % nx;
                    let j = cell / nx;
                    let lower = neighbour(primitive, grid, boundaries, i, j, tangent, false);
                    let upper = neighbour(primitive, grid, boundaries, i, j, tangent, true);
                    (upper[component + 1] - lower[component + 1]) / (2.0 * grid.spacing[tangent])
                };
                let (normal_gradient, transverse_gradient, tangent_normal, tangent_transverse) =
                    if wall {
                        let cell = if face == 0 { right } else { left };
                        let sign = if face == 0 { 1.0 } else { -1.0 };
                        (
                            sign * 2.0 * primitive[cell].velocity[axis] / grid.spacing[axis],
                            0.0,
                            0.0,
                            tangent_derivative(cell, tangent),
                        )
                    } else {
                        (
                            (pr.velocity[axis] - pl.velocity[axis]) / grid.spacing[axis],
                            (pr.velocity[tangent] - pl.velocity[tangent]) / grid.spacing[axis],
                            0.5 * tangent_derivative(left, axis)
                                + 0.5 * tangent_derivative(right, axis),
                            0.5 * tangent_derivative(left, tangent)
                                + 0.5 * tangent_derivative(right, tangent),
                        )
                    };
                let normal_stress = longitudinal * normal_gradient + cross * tangent_transverse;
                let shear_stress = if wall {
                    0.0
                } else {
                    mu * (transverse_gradient + tangent_normal)
                };
                let energy = if wall {
                    0.0
                } else {
                    -normal_stress * (0.5 * pl.velocity[axis] + 0.5 * pr.velocity[axis])
                        - shear_stress * (0.5 * pl.velocity[tangent] + 0.5 * pr.velocity[tangent])
                        - conductivity * (pr.temperature - pl.temperature) / grid.spacing[axis]
                };
                if ![normal_stress, shear_stress, energy]
                    .into_iter()
                    .all(f64::is_finite)
                {
                    return Err(MolecularError::Invalid("nonfinite molecular face flux"));
                }
                fluxes[f] = NormalFaceFlux {
                    mass: 0.0,
                    normal_momentum: -normal_stress,
                    transverse_momentum: -shear_stress,
                    energy,
                    signal_speed: 0.0,
                    used_hll_fallback: false,
                };
            }
        }
    }
    // A deliberately conservative global estimate covers both normal and cross
    // derivative stencils and the doubled half-cell wall gradient. Positivity
    // is additionally checked through the original EOS after every stage.
    let rate = 4.0
        * maximum_diffusivity
        * (grid.spacing[0].recip().powi(2) + grid.spacing[1].recip().powi(2));
    if !rate.is_finite() {
        return Err(MolecularError::Invalid("parabolic scale overflow"));
    }
    Ok(rate)
}
