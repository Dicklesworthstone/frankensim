//! Two finite, material-resolved gas volumes driving ONE translating wall.
//!
//! This is the missing mechanical connection between closed chambers, not two
//! independently stepped pistons. The original gas midpoint law and piston
//! impulse solver are reused. Each gas keeps its own composition, source cards,
//! caloric state, transport, pressure and acoustic properties.

/// Reciprocal pressure/flow work with additional fixed acoustic domains.
pub mod acoustic;
/// Reciprocal circuit force/back-emf on the original shared mechanical state.
pub mod electrical;
/// Geometry-driven internal gas exchange on the same two-sided motion.
pub mod valve;

use super::{
    ClosedGasParcel, GasEnergySolve, GasHeatVolumeCandidate, GasHeatVolumeStep,
    GasPiston, PistonError, PistonGas, PistonGasState, PistonInput, PistonSolve,
    ResolvedThermochemicalGasState,
    Trial, impulse, invalid, poll, positive,
};

/// Two opposed constant-area faces on the SAME rigid translating body.
///
/// Positive velocity expands the primary gas and compresses the opposing gas.
/// Only the primary owner holds mass and velocity. The opposing gas has no
/// fictitious independent inertia or motion coordinate. Different positive face
/// areas are allowed (e.g. a piston rod); any other load is explicit input.
/// Both parcels use the original ideal or Peng--Robinson material owner; each
/// retains its own source model, composition, conserved mass and density.
#[derive(Clone, Debug)]
pub struct DoubleActingGasPiston<G: PistonGas = ClosedGasParcel> {
    primary: GasPiston<G>,
    opposing: G,
    opposing_area_m2: f64,
    ports: crate::port_adoption::ModulePorts,
}
impl<G: PistonGas> PartialEq for DoubleActingGasPiston<G> {
    fn eq(&self, other: &Self) -> bool {
        self.primary == other.primary
            && self.opposing == other.opposing
            && self.opposing_area_m2 == other.opposing_area_m2
    }
}

/// External energy and force over one common physical interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DoublePistonInput {
    /// Positive elapsed time [s].
    pub duration_s: f64,
    /// Integrated heat into [primary, opposing] gas [J], not watts or reactions.
    pub heat_j: [f64; 2],
    /// Additional constant force [N], positive on primary expansion.
    /// Must exclude BOTH gas pressure loads, which this owner supplies.
    pub external_force_n: f64,
}

/// Original nonlinear controls plus the second material's own support bracket.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DoublePistonSolve {
    /// Primary volume bracket and gas controls; the energy tolerance now bounds
    /// BOTH gas energies plus the SINGLE wall kinetic energy and external work.
    pub piston: PistonSolve,
    /// Opposing gas temperature/energy/source work controls, never borrowed from
    /// the primary material. Both gases must admit every requested proposal.
    pub opposing_gas: GasEnergySolve,
    /// Absolute allowed roundoff in dV_opposing + A_opposing*dx [m3].
    /// This does not relax the independently checked joule/impulse balances.
    pub volume_tolerance_m3: f64,
}

/// Joint accepted state accounting. There is NO exterior-pressure supply term:
/// pressure work between these two finite gases and their wall is internal.
#[derive(Clone, Debug, PartialEq)]
pub struct DoublePistonStep<S = ResolvedThermochemicalGasState> {
    /// Original gas heat/work/entropy reports, [primary, opposing].
    pub gas: [GasHeatVolumeStep; 2],
    /// Source-bearing pressure/thermal/acoustic midpoint states in the same order.
    pub midpoint: [S; 2],
    /// Common signed wall displacement [m].
    pub displacement_m: f64,
    /// Initial and accepted final wall velocity [m/s].
    pub velocity_m_s: [f64; 2],
    /// Independently measured change of the ONE wall's kinetic energy [J].
    pub kinetic_energy_change_j: f64,
    /// Actual accepted additional-force work [J]. Acoustic candidates include
    /// their reciprocal pressure reaction; standalone steps use prescribed force.
    pub external_force_work_j: f64,
    /// M*dv - dt*(A_primary*p_primary - A_opposing*p_opposing + F) [N s].
    pub momentum_residual_n_s: f64,
    /// dK - W_primary - W_opposing - W_external [J].
    pub mechanical_energy_residual_j: f64,
    /// dU_primary + dU_opposing + dK - Q_primary - Q_opposing - W_external [J].
    pub total_energy_residual_j: f64,
    /// Actual dV_opposing + A_opposing*dx [m3].
    pub volume_exchange_residual_m3: f64,
    /// Interior joint momentum iterations, excluding initial/bracket proposals.
    pub iterations: usize,
    /// Total gas proposals, including discarded trials: TWO per impulse trial.
    /// Each retains its own bounded caloric inversion and source resolution.
    pub material_trials: usize,
}

/// Evaluate a held affine midpoint load on the original displacement coordinate.
/// Zero resistance preserves the original prescribed-force value bit for bit.
pub(in crate::cavity_phs::piston) fn midpoint_force_n(
    input: DoublePistonInput, resistance_n_s_m: f64, displacement_m: f64,
) -> Result<f64, PistonError> {
    if !positive(input.duration_s) || !input.external_force_n.is_finite()
        || !resistance_n_s_m.is_finite() || resistance_n_s_m < 0.0
        || !displacement_m.is_finite()
    {
        return Err(invalid("invalid affine midpoint force, resistance or displacement"));
    }
    if resistance_n_s_m == 0.0 { return Ok(input.external_force_n); }
    let velocity = displacement_m / input.duration_s;
    let force = input.external_force_n - resistance_n_s_m * velocity;
    if !velocity.is_finite() || !force.is_finite() {
        return Err(invalid("affine midpoint force is unrepresentable"));
    }
    Ok(force)
}

struct PairedTrial<G: PistonGas> {
    primary: Trial<G>,
    opposing: GasHeatVolumeCandidate<G, G::State>,
    primary_input: PistonInput,
    volume_residual: f64,
}

impl<G: PistonGas> DoubleActingGasPiston<G> {
    /// Bind an existing accepted piston to a second closed material volume.
    /// No source is resampled, no gas state is reset, and no second body mass
    /// can be supplied. This constructs a new coupling, not an impulse-free
    /// connection of two formerly moving mechanical bodies.
    pub fn new(
        primary: GasPiston<G>,
        opposing: G,
        opposing_area_m2: f64,
    ) -> Result<Self, PistonError> {
        if primary.support() == super::PistonSupport::Fixed {
            return Err(invalid("double-acting shared motion requires a free primary wall; use independently fixed GasPiston chambers for stationary volumes"));
        }
        if !positive(opposing_area_m2)
            || !positive(opposing.volume_m3() / opposing_area_m2)
            || !positive(opposing_area_m2 / primary.area_m2())
        {
            return Err(invalid("opposed gas needs a positive represented face area and cavity length"));
        }
        use fs_ports::PortKind::{MechanicalForceVelocity as M, ThermalTemperatureEntropy as T};
        let mut ports = crate::port_adoption::ModulePorts::new(
            "double-piston",
            vec![
                ("face/primary", M),
                ("face/opposing", M),
                ("wall", M),
                ("heat/primary", T),
                ("heat/opposing", T),
            ],
        );
        // Both face efforts act on the SAME translating wall coordinate.
        ports.junctions.push(crate::port_adoption::junction(
            "double-piston",
            "wall",
            3,
            M,
            true,
        ));
        Ok(Self {
            primary,
            opposing,
            opposing_area_m2,
            ports,
        })
    }

    /// Original primary dynamics owner. Immutable: a caller cannot advance only
    /// one side and leave the opposing volume or shared momentum behind.
    #[must_use]
    pub const fn primary(&self) -> &GasPiston<G> { &self.primary }
    /// Original opposing material state, conserved gas mass and source context.
    #[must_use]
    pub const fn opposing(&self) -> &G { &self.opposing }
    /// Positive pressure-exposed area of the opposing face [m2].
    #[must_use]
    pub const fn opposing_area_m2(&self) -> f64 { self.opposing_area_m2 }

    fn admit_linear_clock(
        &self, states: [&G::State; 2], volumes_m3: [f64; 2], duration_s: f64,
    ) -> Result<(), PistonError> {
        // Both closed gases supply restoring stiffness to the SAME inertial
        // coordinate. Actual rho*c^2 includes the source EOS departures;
        // either individual face alone can understate the shared frequency.
        let stiffness: f64 = states.into_iter()
            .zip(volumes_m3)
            .zip([self.primary.area_m2(), self.opposing_area_m2])
            .map(|((state, volume), area)| {
                let speed = state.sound_speed_m_s();
                state.density_kg_m3() * speed * speed * area * area / volume
            }).sum();
        let omega = (stiffness / self.primary.mass_kg()).sqrt();
        fs_phs::linear_exact::admit_midpoint_band(omega, duration_s)
            .map_err(PistonError::TimeResolution)
    }

    fn opposing_volume(&self, primary_volume: f64) -> Result<f64, PistonError> {
        let displacement = (primary_volume - self.primary.gas().volume_m3()) / self.primary.area_m2();
        let volume = self.opposing.volume_m3() - self.opposing_area_m2 * displacement;
        if !displacement.is_finite() || !positive(volume)
            || !positive(volume / self.opposing_area_m2)
        {
            return Err(invalid("primary volume proposal collapses or overflows the opposing cavity"));
        }
        Ok(volume)
    }

    fn trial(
        &self, volume: f64, input: DoublePistonInput, controls: DoublePistonSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<impulse::Trial<PairedTrial<G>>, PistonError> {
        self.trial_with_midpoint_resistance(volume, input, controls, 0.0, cancelled)
    }

    fn trial_with_midpoint_resistance(
        &self, volume: f64, input: DoublePistonInput, controls: DoublePistonSolve,
        resistance_n_s_m: f64, cancelled: &mut impl FnMut() -> bool,
    ) -> Result<impulse::Trial<PairedTrial<G>>, PistonError> {
        poll(cancelled)?;
        let opposing_volume = self.opposing_volume(volume)?;
        let opposing = self.opposing.midpoint_heat_volume(
            opposing_volume, input.heat_j[1], controls.opposing_gas, &mut *cancelled,
        )?;
        let ratio = self.opposing_area_m2 / self.primary.area_m2();
        // The original one-sided trial already owns its gas/momentum equations.
        // Its exterior load is now a material proposal, NOT a fixed reservoir.
        let exterior = opposing.midpoint().pressure_pa() * ratio;
        if !positive(exterior) {
            return Err(invalid("opposing pressure load is unrepresentable"));
        }
        let displacement = (volume - self.primary.gas().volume_m3()) / self.primary.area_m2();
        let primary_input = PistonInput {
            duration_s: input.duration_s, heat_j: input.heat_j[0],
            exterior_pressure_pa: exterior,
            external_force_n: midpoint_force_n(input, resistance_n_s_m, displacement)?,
        };
        let primary = self.primary.trial(volume, primary_input, controls.piston, cancelled)?;
        // dV_opposing/dV_primary=-A_opposing/A_primary. Both material tangents
        // are of the SAME discrete heat-volume update, not -gamma*p/V.
        let tangent = primary.tangent - input.duration_s * self.opposing_area_m2 * ratio
            * opposing.step().pressure_final_volume_derivative_pa_m3;
        // d[-dt*(-R*v_mid)]/dV_primary = R/A_primary, with
        // v_mid=(V_primary-V_initial)/(A_primary*dt). This belongs to the
        // original impulse Jacobian; no outer pressure iteration is needed.
        let tangent = if resistance_n_s_m == 0.0 { tangent } else {
            tangent + resistance_n_s_m / self.primary.area_m2()
        };
        let volume_residual = (opposing.endpoint().volume_m3() - self.opposing.volume_m3())
            + self.opposing_area_m2 * primary.displacement;
        if !positive(tangent) || !volume_residual.is_finite()
            || volume_residual.abs() > controls.volume_tolerance_m3
        {
            return Err(invalid("shared-piston discrete tangent or volume exchange failed admission"));
        }
        Ok(impulse::Trial {
            residual: primary.residual, tangent,
            value: PairedTrial { primary, opposing, primary_input, volume_residual },
        })
    }

    /// Jointly solve both gas energies and the single wall momentum.
    ///
    /// dV_primary=A_primary*dx, dV_opposing=-A_opposing*dx,
    /// dx=dt*(v0+v1)/2, M*dv=dt*(A_primary*pm_primary-A_opposing*pm_opposing+F).
    /// Every proposal starts from the same two accepted material states.
    /// Heat, pressure, transport and acoustic response are recomputed by each
    /// original gas owner; different species or caloric data need no solver edit.
    ///
    /// Both geometric endpoint volumes must be positive BEFORE querying sources.
    /// At most 2*(maximum_iterations+3) gas proposals are made. Sources, geometry,
    /// balances, exhausted work and cancellation cannot publish either side alone.
    /// This closed step has no mass/species exchange; the ideal-parcel
    /// [`DoubleActingGasPiston::valve_candidate`] adds compressible transport
    /// through a moving metering edge. Initial, midpoint and endpoint acoustic
    /// admission include both source-backed gas stiffnesses on the single wall.
    /// Contact stops, crank linkages, combustion, wall heat capacity and acoustic
    /// propagation remain separate owners.
    pub fn step(
        &mut self,
        input: DoublePistonInput,
        controls: DoublePistonSolve,
        cancelled: impl FnMut() -> bool,
    ) -> Result<DoublePistonStep<G::State>, PistonError> {
        self.step_with_midpoint_resistance(input, controls, 0.0, cancelled)
    }

    /// Internal reciprocal-load seam. The resistance is a reaction coefficient,
    /// not automatically heat; its work remains in the native external-force
    /// account for the composing acoustic owner to cancel.
    pub(in crate::cavity_phs::piston) fn step_with_midpoint_resistance(
        &mut self,
        input: DoublePistonInput,
        controls: DoublePistonSolve,
        resistance_n_s_m: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DoublePistonStep<G::State>, PistonError> {
        poll(&mut cancelled)?;
        let [lo, hi] = controls.piston.volume_bracket_m3;
        if !positive(input.duration_s) || !input.heat_j.iter().all(|q| q.is_finite())
            || !input.external_force_n.is_finite()
            || !resistance_n_s_m.is_finite() || resistance_n_s_m < 0.0
            || !positive(controls.volume_tolerance_m3)
            || !positive(controls.piston.energy_tolerance_j)
            || !positive(lo) || !positive(hi) || lo >= hi
            || self.primary.gas().volume_m3() < lo || self.primary.gas().volume_m3() > hi
        {
            return Err(invalid("invalid shared-piston input, volume bracket or balance tolerances"));
        }
        self.opposing_volume(lo)?;
        self.opposing_volume(hi)?;
        self.admit_linear_clock(
            [self.primary.gas().state(), self.opposing.state()],
            [self.primary.gas().volume_m3(), self.opposing.volume_m3()], input.duration_s,
        )?;
        let solved = impulse::solve(self.primary.gas().volume_m3(), controls.piston, |volume| {
            if resistance_n_s_m == 0.0 {
                self.trial(volume, input, controls, &mut cancelled)
            } else {
                self.trial_with_midpoint_resistance(
                    volume, input, controls, resistance_n_s_m, &mut cancelled)
            }
        })?;
        let PairedTrial { primary, opposing, primary_input, volume_residual } = solved.value;
        self.admit_linear_clock(
            [primary.gas.midpoint(), opposing.midpoint()],
            [
                f64::midpoint(self.primary.gas().volume_m3(), primary.gas.endpoint().volume_m3()),
                f64::midpoint(self.opposing.volume_m3(), opposing.endpoint().volume_m3()),
            ], input.duration_s,
        )?;
        self.admit_linear_clock(
            [primary.gas.endpoint().state(), opposing.endpoint().state()],
            [primary.gas.endpoint().volume_m3(), opposing.endpoint().volume_m3()], input.duration_s,
        )?;
        // The accepted force, not the affine intercept, supplies all original
        // momentum, energy and port-account audits at the actual displacement.
        let input = DoublePistonInput {
            external_force_n: primary_input.external_force_n, ..input
        };
        // Finish the original owner only on a private candidate. The accepted
        // pair remains untouched even if the final all-owner audit is refused.
        let opposing_midpoint = opposing.midpoint().clone();
        let mut next_primary = self.primary.clone();
        let primary_step = next_primary.finish(primary, primary_input, controls.piston,
            solved.iterations, solved.proposals, &mut cancelled)?;
        let (next_opposing, opposing_step) = opposing.into_parts();
        let force = self.primary.area_m2() * primary_step.midpoint.pressure_pa()
            - self.opposing_area_m2 * opposing_midpoint.pressure_pa() + input.external_force_n;
        let momentum = self.primary.mass_kg()
            * (primary_step.velocity_m_s[1] - primary_step.velocity_m_s[0]) - input.duration_s * force;
        let mechanical = (primary_step.kinetic_energy_change_j - primary_step.gas.work_by_gas_j)
            - opposing_step.work_by_gas_j - primary_step.external_force_work_j;
        let total = ((primary_step.gas.stored_energy_change_j + opposing_step.stored_energy_change_j)
            + primary_step.kinetic_energy_change_j) - input.heat_j[0] - input.heat_j[1]
            - primary_step.external_force_work_j;
        if ![force, momentum, mechanical, total].iter().all(|x| x.is_finite())
            || momentum.abs() > controls.piston.momentum_tolerance_n_s
            || mechanical.abs() > controls.piston.energy_tolerance_j
            || total.abs() > controls.piston.energy_tolerance_j
        {
            return Err(invalid("joint gas/piston momentum or independently measured total energy failed admission"));
        }
        let report = DoublePistonStep {
            gas: [primary_step.gas, opposing_step],
            midpoint: [primary_step.midpoint, opposing_midpoint],
            displacement_m: primary_step.displacement_m,
            velocity_m_s: primary_step.velocity_m_s,
            kinetic_energy_change_j: primary_step.kinetic_energy_change_j,
            external_force_work_j: primary_step.external_force_work_j,
            momentum_residual_n_s: momentum, mechanical_energy_residual_j: mechanical,
            total_energy_residual_j: total, volume_exchange_residual_m3: volume_residual,
            iterations: solved.iterations, material_trials: 2 * solved.proposals,
        };
        poll(&mut cancelled)?;
        self.primary = next_primary;
        self.opposing = next_opposing;
        let force = [
            self.primary.area_m2() * report.midpoint[0].pressure_pa(),
            -self.opposing_area_m2 * report.midpoint[1].pressure_pa(),
        ];
        let face = [
            force[0] * report.displacement_m,
            force[1] * report.displacement_m,
        ];
        let wall_work = -(force[0] + force[1]) * report.displacement_m;
        self.ports.accept(
            crate::port_adoption::StepAccount {
                tick: self.ports.next_tick(),
                duration_s: input.duration_s,
                storage_change_j: Some(
                    report.gas[0].stored_energy_change_j
                        + report.gas[1].stored_energy_change_j
                        + report.kinetic_energy_change_j,
                ),
                loss_j: Some(0.0),
                included_work_j: report.external_force_work_j,
                entropy_change_j_k: report.gas[0].entropy_change_j_k
                    + report.gas[1].entropy_change_j_k,
                tolerance_j: Some(controls.piston.energy_tolerance_j),
                unavailable: None,
            },
            [
                crate::port_adoption::PortTransfer::work(face[0], false),
                crate::port_adoption::PortTransfer::work(face[1], false),
                crate::port_adoption::PortTransfer::work(wall_work, false),
                crate::port_adoption::PortTransfer::heat(
                    input.heat_j[0],
                    report.midpoint[0].temperature_k(),
                ),
                crate::port_adoption::PortTransfer::heat(
                    input.heat_j[1],
                    report.midpoint[1].temperature_k(),
                ),
            ],
        );
        Ok(report)
    }
}

#[cfg(test)]
mod tests;

#[cfg(test)]
pub(super) mod real_gas;

impl<G: PistonGas> fs_ports::PortedModule for DoubleActingGasPiston<G> {
    fn port_schemas(&self) -> &[fs_ports::PortSchema] { &self.ports.schemas }
    fn port_coordinate_charts(&self) -> &[fs_ports::PortCoordinateChart] { &self.ports.charts }
    fn junction_schemas(&self) -> &[fs_ports::JunctionSchema] { &self.ports.junctions }
    fn last_window_duration_s(&self) -> Option<f64> { self.ports.duration_s() }
    fn last_window_transfers(&self) -> Option<&[fs_ports::IntegratedWindowTransfer]> {
        self.ports.transfers()
    }
    fn last_window_transfer_error(&self) -> Option<&fs_ports::CoupleError> {
        self.ports.transfer_error()
    }
    fn audit_last_window(&self) -> Result<fs_ports::ModuleWindowAudit, fs_ports::CoupleError> {
        self.ports.audit()
    }
}
