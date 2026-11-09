//! Moving-conductor Joule work on the original material/enthalpy endpoint.
//!
//! The native field supplies `J = Sigma(T,f) (E + u cross B)`. Only its
//! conductor-frame loss heats the original receiving cells. Signed Lorentz
//! work leaves through prescribed motion, so the thermal driver receives
//! terminal work minus mechanical output, including mechanical braking input.

use super::*;
use fs_em::spatial::conduction::motion::{
    SpatialConductionMotion, SpatialMovingConductionSnapshot,
};

/// Original moving-field solution and separately integrated physical work.
#[derive(Clone, Debug)]
pub struct SpatialMovingConductionThermalStep {
    pub start_s: f64,
    pub end_s: f64,
    pub field: SpatialMovingConductionSnapshot,
    /// Direct sparse electrode observations retain physical voltage differences
    /// and signed authored source work independently of absolute voltage bias.
    pub terminals: Option<ConductionTerminalReport>,
    /// Stationary interface losses and their explicit receiving-cell partition.
    /// This surface heat is included in `heat_j` and deposited exactly once.
    pub contacts: Vec<SpatialConductionContactHeat>,
    /// Conductor-frame Joule heat already deposited into the finite cells [J].
    pub heat_j: f64,
    /// Signed electrical work INTO the conductor at its native boundary [J].
    pub terminal_work_j: f64,
    /// Signed Lorentz work OUT to the prescribed motion owner [J]. Negative
    /// work is mechanical input, as in short-circuit magnetic braking.
    pub mechanical_work_j: f64,
    /// Electrical input minus mechanical output [J]. This is the driver's
    /// supply to finite heat; it is not an additional thermal source.
    pub net_supplied_work_j: f64,
}

struct Candidate<C> {
    problem: SpatialConduction<C>,
    step: SpatialMovingConductionThermalStep,
}

impl<C: SpatialConductorMaterial> ThermalDriver for Candidate<C> {
    type Model = SpatialConduction<C>;
    type Step = SpatialMovingConductionThermalStep;

    fn account(&self) -> DriverAccount<'_> {
        DriverAccount {
            resistors: &[],
            contacts: &[],
            heat_j: self.step.heat_j,
            numerical_dissipation_j: 0.0,
            stored_change_j: 0.0,
            supplied_j: self.step.net_supplied_work_j,
        }
    }

    fn into_parts(self) -> (Self::Model, Self::Step) {
        (self.problem, self.step)
    }
}

/// One jointly admitted thermal history, source operator and motion-work report.
pub struct SpatialMovingConductionEnthalpyCandidate<
    't,
    C = ResolvedConductorState,
    S = MaterialNetworkSample,
> {
    thermal: DepositionSimulation<'t>,
    problem: SpatialConduction<C>,
    report: CircuitEnthalpyReport<S, SpatialMovingConductionThermalStep>,
}

impl<'t, C: SpatialConductorMaterial, S> SpatialMovingConductionEnthalpyCandidate<'t, C, S> {
    pub fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }

    pub fn problem(&self) -> &SpatialConduction<C> {
        &self.problem
    }

    pub fn report(&self) -> &CircuitEnthalpyReport<S, SpatialMovingConductionThermalStep> {
        &self.report
    }

    /// Continue with both returned owners. Joule heat has already entered this
    /// history; mechanical work belongs to the external prescribed-motion owner.
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        SpatialConduction<C>,
        CircuitEnthalpyReport<S, SpatialMovingConductionThermalStep>,
    ) {
        (self.thermal, self.problem, self.report)
    }
}

impl SpatialConductionThermalNetwork<'_> {
    /// Close prescribed moving-current transport, original source conductivity
    /// and finite-cell heat at the SAME backward-Euler endpoint.
    ///
    /// Every trial refreshes the original scalar or fixed-orientation tensor
    /// source from its actual receiving temperature/phase, then solves charge
    /// continuity with the held laboratory velocity and magnetic flux density.
    /// Each original tetrahedron's `integral((E + u cross B).J)` enters its
    /// original receiving cell once. Laboratory `E.J` and signed Lorentz work
    /// are not heat. Terminal work minus mechanical output must independently
    /// close the deposited Joule energy and native thermal storage.
    ///
    /// Voltage and distributed-current boundary traces retain their native
    /// meaning; current-only open components retain their charge-neutral gauge.
    /// [`SpatialConductionDrive`] also admits total-current and mixed-terminal
    /// controls. Direct `Terminals` drives identify the original electrodes and
    /// their motional load in one sparse solve at each actual endpoint. Matrix
    /// drives rebuild their conductance response and add the motional current
    /// before solving the authored controls.
    /// An open-circuit electrode therefore develops its motional voltage; zero
    /// authored current does not mean zero electric field. Held terminal power
    /// is not admitted in this moving path.
    /// For finite-resistance stationary interfaces use
    /// [`Self::step_moving_contact_candidate`] with explicit heat bindings.
    /// All original material identities, volume/mass guards, source domains,
    /// thermal links, events, solve budgets and cancellation gates apply.
    /// Rejected trials publish neither the thermal nor electrical endpoint.
    ///
    /// This fixed-domain quasistatic approximation requires the prescribed
    /// motion to preserve the represented material/source occupancy and finite
    /// thermal inventory during the interval. It does not transport material,
    /// evolve induction or momentum, or account for a prescribed magnet's field
    /// storage. The caller owns that physical validity and the reported motion
    /// work. Geometry, velocity and B are held; temporal/spatial refinement is
    /// still required independently of the conservation gates.
    #[allow(clippy::too_many_arguments)]
    pub fn step_moving_candidate<'t, 'd, C: SpatialConductorMaterial>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        problem: &SpatialConduction<C>,
        drive: impl Into<SpatialConductionDrive<'d>>,
        motion: SpatialConductionMotion<'_>,
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<SpatialMovingConductionEnthalpyCandidate<'t, C>, Error> {
        self.step_moving_contact_candidate(
            thermal,
            target_s,
            problem,
            drive,
            motion,
            &[],
            assembly,
            field_controls,
            controls,
            initial_energy_increments_j,
            cancelled,
        )
    }

    /// Advance moving bulk conductors and stationary resistive interfaces in
    /// the same finite-enthalpy endpoint.
    ///
    /// Each interface must have zero prescribed velocity in both adjacent
    /// tetrahedra. The native field refuses moving or sliding interfaces, which
    /// require a separate surface-work law. Every accepted trial refreshes the
    /// original interface source at its two receiving temperature/phase states.
    /// Its declared heat partition deposits the native contact loss alongside
    /// bulk conductor-frame loss exactly once. Terminal work minus bulk Lorentz
    /// output closes both losses and the original finite thermal storage.
    ///
    /// Geometry, field and prescribed-motion validity, work budgets, source
    /// domains and transactional refusals are those of
    /// [`Self::step_moving_candidate`].
    #[allow(clippy::too_many_arguments)]
    pub fn step_moving_contact_candidate<'t, 'd, C: SpatialConductorMaterial>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        problem: &SpatialConduction<C>,
        drive: impl Into<SpatialConductionDrive<'d>>,
        motion: SpatialConductionMotion<'_>,
        contact_bindings: &[ConductionContactThermal<'_>],
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<SpatialMovingConductionEnthalpyCandidate<'t, C>, Error> {
        let cancel = RefCell::new(cancelled);
        check(&cancel)?;
        let network = MaterialNetwork {
            cells: self.cells,
            resistors: &[],
            links: self.links,
        };
        let prepared = network
            .prepare_with_bases(thermal, &[], || (*cancel.borrow_mut())())
            .map_err(material)?;
        self.step_moving_resolved_candidate(
            thermal,
            target_s,
            problem,
            drive,
            motion,
            contact_bindings,
            assembly,
            field_controls,
            controls,
            initial_energy_increments_j,
            |trial| {
                prepared
                    .resolve(trial, || (*cancel.borrow_mut())())
                    .map_err(material)
            },
            |sample| sample,
            &cancel,
        )
    }

    // One moving-field/contact driver for the ordinary material network and
    // complete sourced heat assemblies; only their existing heat resolver varies.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub(crate) fn step_moving_resolved_candidate<
        't,
        'd,
        C: SpatialConductorMaterial,
        S,
        E: From<ConductionThermalError>,
    >(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        problem: &SpatialConduction<C>,
        drive: impl Into<SpatialConductionDrive<'d>>,
        motion: SpatialConductionMotion<'_>,
        contact_bindings: &[ConductionContactThermal<'_>],
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        resolve: impl Fn(
            &DepositionSimulation<'t>,
        ) -> Result<CircuitThermalCoefficients<S>, CircuitEnthalpyError<E>>,
        sample: impl Fn(&S) -> &MaterialNetworkSample,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<SpatialMovingConductionEnthalpyCandidate<'t, C, S>, CircuitEnthalpyError<E>> {
        let drive = drive.into();
        // Direct electrodes close their motional load in one actual-state solve.
        // Prepared responses need N-1 basis solves and one motional-source solve.
        // Charge every solve, including a zero-current open-circuit field.
        let field_solves = match drive {
            SpatialConductionDrive::Terminals(_) | SpatialConductionDrive::Boundary(_) => 1,
            SpatialConductionDrive::TotalCurrent { .. } => 2,
            SpatialConductionDrive::TerminalCurrents { terminals, .. }
            | SpatialConductionDrive::MixedTerminals { terminals, .. } => {
                if terminals.len() < 2 {
                    return Err(map_conduction_error(invalid(
                        "moving terminal drive needs at least two electrodes",
                    )));
                }
                terminals.len()
            }
            SpatialConductionDrive::Power { .. } => {
                return Err(map_conduction_error(invalid(
                    "held terminal power is not admitted with prescribed motional EMF",
                )));
            }
        };
        let dt = self
            .admit_step(
                thermal,
                target_s,
                problem,
                assembly,
                controls,
                initial_energy_increments_j,
            )
            .map_err(map_conduction_error)?;
        let (thermal, problem, report) = self.step_driver_candidate_with_geometry(
            thermal,
            target_s,
            problem,
            contact_bindings,
            assembly,
            controls,
            initial_energy_increments_j,
            resolve,
            sample,
            None,
            |_| Ok(None),
            &[],
            |_| Ok(Vec::new()),
            &[],
            |next, _resistors, contact_owner, remaining| {
                check(cancel).map_err(map_conduction_error)?;
                if field_solves > remaining {
                    return Err(CircuitEnthalpyError::CircuitBudget);
                }
                let mut terminal_report = None;
                let field = match drive {
                    SpatialConductionDrive::Terminals(drive) => next
                        .solve_terminals_with_motion(
                            drive.terminals,
                            drive.control,
                            drive.limits,
                            motion,
                            field_controls,
                            &mut || (*cancel.borrow_mut())(),
                        )
                        .map(|(field, report)| {
                            terminal_report = Some(report);
                            field
                        }),
                    SpatialConductionDrive::Boundary(boundary) => {
                        next.solve_with_motion(boundary, motion, field_controls, &mut || {
                            (*cancel.borrow_mut())()
                        })
                    }
                    SpatialConductionDrive::TotalCurrent {
                        positive_faces,
                        negative_faces,
                        current_a,
                        reference_voltage_v,
                    } => {
                        let response = next
                            .resistance_at_reference_voltage(
                                positive_faces,
                                negative_faces,
                                reference_voltage_v,
                                field_controls,
                                &mut || (*cancel.borrow_mut())(),
                            )
                            .map_err(field_error)
                            .map_err(map_conduction_error)?;
                        response.field_at_current_with_motion(
                            current_a,
                            motion,
                            field_controls,
                            &mut || (*cancel.borrow_mut())(),
                        )
                    }
                    SpatialConductionDrive::TerminalCurrents {
                        terminals,
                        currents_a,
                        response,
                    } => {
                        let matrix = next
                            .conductance_matrix(terminals, response, field_controls, &mut || {
                                (*cancel.borrow_mut())()
                            })
                            .map_err(field_error)
                            .map_err(map_conduction_error)?;
                        matrix.field_at_currents_with_motion(
                            currents_a,
                            motion,
                            field_controls,
                            &mut || (*cancel.borrow_mut())(),
                        )
                    }
                    SpatialConductionDrive::MixedTerminals {
                        terminals,
                        current_terminals,
                        voltages_v,
                        currents_a,
                        response,
                        mixed,
                    } => {
                        let matrix = next
                            .conductance_matrix(terminals, response, field_controls, &mut || {
                                (*cancel.borrow_mut())()
                            })
                            .map_err(field_error)
                            .map_err(map_conduction_error)?;
                        let mixed = matrix
                            .prepare_current_control(current_terminals, mixed, &mut || {
                                (*cancel.borrow_mut())()
                            })
                            .map_err(field_error)
                            .map_err(map_conduction_error)?;
                        mixed.field_at_controls_with_motion(
                            voltages_v,
                            currents_a,
                            motion,
                            field_controls,
                            &mut || (*cancel.borrow_mut())(),
                        )
                    }
                    SpatialConductionDrive::Power { .. } => {
                        return Err(map_conduction_error(invalid(
                            "held terminal power is not admitted with prescribed motional EMF",
                        )));
                    }
                }
                .map_err(field_error)
                .map_err(map_conduction_error)?;
                let heat_j = finite(dt * field.joule_power_w())?;
                let terminal_power_w = match &terminal_report {
                    Some(report) => report.supplied_power_w,
                    None => drive
                        .supplied_power_from_terminals_w(
                            next.mesh(),
                            field.voltage_v(),
                            field.nodal_boundary_current_a(),
                            field.terminal_power_w(),
                            cancel,
                        )
                        .map_err(map_conduction_error)?,
                };
                let terminal_work_j = finite(dt * terminal_power_w)?;
                let mechanical_work_j = finite(dt * field.mechanical_power_w())?;
                let net_supplied_work_j = finite(terminal_work_j - mechanical_work_j)?;
                if (field.joule_power_w() > 0.0 && heat_j == 0.0)
                    || finite(terminal_work_j - finite(dt * field.terminal_power_w())?)?.abs()
                        > controls.total_energy_tolerance_j
                    || finite(heat_j - net_supplied_work_j)?.abs()
                        > controls.total_energy_tolerance_j
                {
                    return Err(map_conduction_error(invalid(
                        "conductor-frame heat, terminal work and Lorentz work do not close",
                    )));
                }
                let contacts = contact_owner
                    .heat(field.contact_response())
                    .map_err(map_conduction_error)?;
                let field = field.into_snapshot();
                check(cancel).map_err(map_conduction_error)?;
                Ok((
                    Candidate {
                        problem: next,
                        step: SpatialMovingConductionThermalStep {
                            start_s: thermal.time_s(),
                            end_s: target_s,
                            field,
                            terminals: terminal_report,
                            contacts,
                            heat_j,
                            terminal_work_j,
                            mechanical_work_j,
                            net_supplied_work_j,
                        },
                    },
                    field_solves,
                ))
            },
            |candidate| {
                (
                    candidate.step.field.cell_joule_power_w(),
                    &candidate.step.contacts,
                )
            },
            cancel,
        )?;
        Ok(SpatialMovingConductionEnthalpyCandidate {
            thermal,
            problem,
            report,
        })
    }
}
