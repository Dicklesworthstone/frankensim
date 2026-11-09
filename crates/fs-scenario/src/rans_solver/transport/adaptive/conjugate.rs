//! Whole-window time control for the original coupled solid/fluid step.
//! The common adaptive controller sees both temperature histories; accepted
//! interface exchanges and all external/generated heat are integrated once.

use super::super::{
    RansConjugateHeatLoad, RansConjugateReport, RansConjugateWall, RansConjugateWorkspace,
};
use super::*;

/// Reusable complete fluid trials, one wall material image and three wall
/// temperature images. Numerical payload is 24 f64 fields per fluid node,
/// plus ten f64 values per wall cell (including retained wall material).
/// Stepping makes no numerical allocation.
pub struct RansConjugateAdaptiveWorkspace {
    kernel: RansConjugateWorkspace,
    accepted: RansFieldState,
    coarse: RansFieldState,
    fine: RansFieldState,
    wall: RansConjugateWall,
    accepted_wall: Vec<f64>,
    coarse_wall: Vec<f64>,
}

impl RansConjugateAdaptiveWorkspace {
    /// Reserve under separate fluid-node and wall-layer ceilings. The
    /// workspace binds the exact wall material, but accepts new temperatures.
    pub fn try_new(
        grid: &RansChannelGrid,
        wall: &RansConjugateWall,
        maximum_fluid_nodes: usize,
        maximum_wall_layers: usize,
    ) -> Result<Self, RansTransportError> {
        let n = grid
            .num_cells
            .checked_add(1)
            .ok_or(RansTransportError::Allocation)?;
        let m = wall.layers().len();
        if n > maximum_fluid_nodes || m > maximum_wall_layers {
            return Err(RansTransportError::NodeBudget {
                needed: if n > maximum_fluid_nodes { n } else { m },
                maximum: if n > maximum_fluid_nodes {
                    maximum_fluid_nodes
                } else {
                    maximum_wall_layers
                },
            });
        }
        n.checked_mul(24)
            .and_then(|fluid| m.checked_mul(10).and_then(|solid| fluid.checked_add(solid)))
            .and_then(|count| count.checked_mul(core::mem::size_of::<f64>()))
            .ok_or(RansTransportError::Allocation)?;
        Ok(Self {
            kernel: RansConjugateWorkspace::try_new(
                grid,
                wall,
                maximum_fluid_nodes,
                maximum_wall_layers,
            )?,
            accepted: state_image(n)?,
            coarse: state_image(n)?,
            fine: state_image(n)?,
            wall: RansConjugateWall::try_new(
                wall.layers(),
                wall.temperatures_k(),
                maximum_wall_layers,
            )?,
            accepted_wall: field_image(m)?,
            coarse_wall: field_image(m)?,
        })
    }

    /// Retained numerical capacities, including the private wall material.
    /// Fixed headers and caller-owned state are excluded.
    #[must_use]
    pub fn numerical_bytes(&self) -> usize {
        let states = [&self.accepted, &self.coarse, &self.fine];
        let fields = states
            .into_iter()
            .flat_map(|s| [&s.u, &s.k, &s.eps, &s.nu_t, &s.temp_k]);
        self.kernel.numerical_bytes()
            + self.wall.numerical_bytes()
            + (fields.map(Vec::capacity).sum::<usize>()
                + self.accepted_wall.capacity()
                + self.coarse_wall.capacity())
                * core::mem::size_of::<f64>()
    }
}

/// Wall and interface heat integrated only along accepted fine half steps,
/// per unit wall area. These are diagnostics, not enforced error bounds.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RansAdaptiveWallBalance {
    /// Original solid sensible heat [J/m²].
    pub energy_before_j_m2: f64,
    /// Final solid sensible heat [J/m²].
    pub energy_after_j_m2: f64,
    /// Compensated represented solid inventory increments [J/m²].
    pub energy_change_j_m2: f64,
    /// Integrated imposed exterior heat, positive into the solid [J/m²].
    pub exterior_heat_j_m2: f64,
    /// Integrated solid volumetric heating [J/m²].
    pub volume_heat_j_m2: f64,
    /// Integrated actual solid-to-fluid interface heat [J/m²].
    pub interface_heat_j_m2: f64,
    /// Actual interface heat minus the fluid boundary exchange [J/m²].
    pub interface_equation_defect_j_m2: f64,
    /// Solid storage minus exterior/volume heat plus interface heat [J/m²].
    pub balance_defect_j_m2: f64,
    /// Combined solid/fluid thermal storage minus external/generated heat [J/m²].
    pub combined_thermal_balance_defect_j_m2: f64,
}

/// A completed joint physical window. The fluid heat account's wall input
/// is the eliminated interface exchange, not the solid's exterior input.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RansConjugateAdaptiveReport {
    /// Shared adaptive work/error diagnostics, mechanics and fluid heat.
    pub fluid: RansAdaptiveReport,
    /// Independent solid and joint thermal accounts.
    pub wall: RansAdaptiveWallBalance,
}

impl RansSolver {
    /// Advance a finite conducting wall and transported channel over a complete
    /// constant-load window. The same full-versus-two-half-step controller
    /// controls all original fluid fields AND every solid temperature. Wall
    /// temperatures use `controls.absolute.temperature_k` and the common
    /// relative allowance. The fractional-change bound is checked at both fine
    /// endpoints; local error compares the full and refined final endpoints.
    ///
    /// Every trial invokes the original monolithic wall/fluid thermal kernel;
    /// turbulence and optional dissipative heating refresh at each half step.
    /// Only accepted fine steps enter the energy accounts. Rejected trials,
    /// cancellation, work exhaustion and late arithmetic cannot publish either
    /// caller state. No extrapolation, thermal clipping or heat repair is used.
    ///
    /// This controls first-order local time error, not global error, forcing
    /// events, empirical model error or unresolved thermal modes below the
    /// authored tolerances. The original constant-property channel scope holds.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_launder_sharma_conjugate_adaptive(
        &self,
        grid: &RansChannelGrid,
        state: &mut RansFieldState,
        wall: &mut RansConjugateWall,
        pressure_gradient_pa_m: f64,
        duration_s: f64,
        heat: RansConjugateHeatLoad<'_>,
        controls: RansAdaptiveControls,
        work: &mut RansConjugateAdaptiveWorkspace,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RansConjugateAdaptiveReport, RansAdaptiveError> {
        checkpoint(0, &mut cancelled)?;
        validate(grid, state, duration_s, controls, work.accepted.u.len())?;
        if wall.layers().len() != work.wall.layers().len() {
            return Err(RansAdaptiveError::InvalidInput(
                "adaptive wall material binding",
            ));
        }
        for (i, (actual, bound)) in wall.layers().iter().zip(work.wall.layers()).enumerate() {
            checkpoint(i, &mut cancelled)?;
            if actual != bound {
                return Err(RansAdaptiveError::InvalidInput(
                    "adaptive wall material binding",
                ));
            }
        }
        copy_temperatures(
            &mut work.accepted_wall,
            wall.temperatures_k(),
            &mut cancelled,
        )?;
        let mut stepper = ConjugateStepper {
            solver: self,
            grid,
            pressure: pressure_gradient_pa_m,
            heat,
            kernel: &mut work.kernel,
            wall: &mut work.wall,
            accepted: &mut work.accepted_wall,
            coarse: &mut work.coarse_wall,
            account: WallAccount::default(),
        };
        let fluid = advance_window(
            state,
            duration_s,
            controls,
            &mut work.accepted,
            &mut work.coarse,
            &mut work.fine,
            &mut stepper,
            &mut cancelled,
        )?;
        let wall_report = stepper.account.finish()?;
        checkpoint(0, &mut cancelled)?;
        // Exclusive borrows and no fallible work beyond this joint publication.
        publish_state(state, &work.accepted);
        wall.temperatures_mut().copy_from_slice(&work.accepted_wall);
        Ok(RansConjugateAdaptiveReport {
            fluid,
            wall: wall_report,
        })
    }
}

fn copy_temperatures(
    to: &mut [f64],
    from: &[f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), RansTransportError> {
    for i in (0..from.len()).step_by(64) {
        checkpoint(i, cancelled)?;
        let end = (i + 64).min(from.len());
        to[i..end].copy_from_slice(&from[i..end]);
    }
    Ok(())
}

struct ConjugateStepper<'a> {
    solver: &'a RansSolver,
    grid: &'a RansChannelGrid,
    pressure: f64,
    heat: RansConjugateHeatLoad<'a>,
    kernel: &'a mut RansConjugateWorkspace,
    wall: &'a mut RansConjugateWall,
    accepted: &'a mut [f64],
    coarse: &'a mut [f64],
    account: WallAccount,
}

impl WindowStepper for ConjugateStepper<'_> {
    fn heat(&self) -> bool {
        true
    }

    fn coarse(&mut self, cancelled: &mut impl FnMut() -> bool) -> Result<(), RansTransportError> {
        copy_temperatures(self.wall.temperatures_mut(), self.accepted, cancelled)
    }

    fn fine(&mut self, cancelled: &mut impl FnMut() -> bool) -> Result<(), RansTransportError> {
        copy_temperatures(self.coarse, self.wall.temperatures_k(), cancelled)?;
        copy_temperatures(self.wall.temperatures_mut(), self.accepted, cancelled)
    }

    fn step(
        &mut self,
        state: &mut RansFieldState,
        dt: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Step, RansTransportError> {
        let report = self.solver.advance_launder_sharma_with_conjugate_wall(
            self.grid,
            state,
            self.wall,
            self.pressure,
            dt,
            self.heat,
            self.kernel,
            cancelled,
        )?;
        Ok(Step {
            mechanical: report.fluid.mechanical,
            heat: Some(report.fluid),
            wall: Some(report),
        })
    }

    fn change(
        &self,
        c: RansAdaptiveControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, RansTransportError> {
        self.ratio(c, false, cancelled)
    }

    fn error(
        &self,
        c: RansAdaptiveControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, RansTransportError> {
        self.ratio(c, true, cancelled)
    }

    fn accept(
        &mut self,
        first: Step,
        second: Step,
        dt: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), RansTransportError> {
        self.account.add(
            first
                .wall
                .ok_or(RansTransportError::InvalidInput("missing wall step"))?,
            dt,
        )?;
        self.account.add(
            second
                .wall
                .ok_or(RansTransportError::InvalidInput("missing wall step"))?,
            dt,
        )?;
        copy_temperatures(self.accepted, self.wall.temperatures_k(), cancelled)
    }
}

impl ConjugateStepper<'_> {
    fn ratio(
        &self,
        c: RansAdaptiveControls,
        error: bool,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, RansTransportError> {
        let mut result = 0.0_f64;
        for (i, &fine) in self.wall.temperatures_k().iter().enumerate() {
            checkpoint(i, cancelled)?;
            let magnitude = self.accepted[i].abs().max(fine.abs());
            let scale = c.absolute.temperature_k
                + if error {
                    c.relative_tolerance * magnitude
                } else {
                    magnitude
                };
            if !scale.is_finite() || scale <= 0.0 {
                return Err(RansTransportError::InvalidInput(
                    "adaptive wall temperature scale",
                ));
            }
            let other = if error {
                self.coarse[i]
            } else {
                self.accepted[i]
            };
            result = result.max((fine - other).abs() / scale);
        }
        Ok(result)
    }
}

#[derive(Default)]
struct WallAccount {
    first: Option<RansConjugateReport>,
    last: Option<RansConjugateReport>,
    sums: [Sum; 7],
}

impl WallAccount {
    fn add(&mut self, r: RansConjugateReport, dt: f64) -> Result<(), RansTransportError> {
        let values = [
            r.wall_energy_change_j_m2,
            dt * r.exterior_heat_flux_w_m2,
            dt * r.wall_volumetric_heat_input_w_m2,
            dt * r.interface_heat_flux_w_m2,
            dt * r.interface_equation_defect_w_m2,
            dt * r.wall_balance_defect_w_m2,
            dt * r.combined_thermal_balance_defect_w_m2,
        ];
        for (sum, value) in self.sums.iter_mut().zip(values) {
            sum.add(value);
        }
        if self.sums.iter().any(|sum| !sum.value.is_finite()) {
            return Err(RansTransportError::Arithmetic {
                equation: "adaptive wall balance",
                node: 0,
            });
        }
        self.first.get_or_insert(r);
        self.last = Some(r);
        Ok(())
    }

    fn finish(self) -> Result<RansAdaptiveWallBalance, RansTransportError> {
        let first = self.first.ok_or(RansTransportError::InvalidInput(
            "no accepted wall interval",
        ))?;
        let last = self.last.ok_or(RansTransportError::InvalidInput(
            "no accepted wall interval",
        ))?;
        Ok(RansAdaptiveWallBalance {
            energy_before_j_m2: first.wall_energy_before_j_m2,
            energy_after_j_m2: last.wall_energy_after_j_m2,
            energy_change_j_m2: self.sums[0].value,
            exterior_heat_j_m2: self.sums[1].value,
            volume_heat_j_m2: self.sums[2].value,
            interface_heat_j_m2: self.sums[3].value,
            interface_equation_defect_j_m2: self.sums[4].value,
            balance_defect_j_m2: self.sums[5].value,
            combined_thermal_balance_defect_j_m2: self.sums[6].value,
        })
    }
}
