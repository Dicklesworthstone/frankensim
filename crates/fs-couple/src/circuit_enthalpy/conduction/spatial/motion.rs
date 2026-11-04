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
    pub fn step_moving_candidate<'t, C: SpatialConductorMaterial>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        problem: &SpatialConduction<C>,
        boundary: &SpatialConductionBoundary,
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
            boundary,
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
    pub fn step_moving_contact_candidate<'t, C: SpatialConductorMaterial>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        problem: &SpatialConduction<C>,
        boundary: &SpatialConductionBoundary,
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
            boundary,
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
        C: SpatialConductorMaterial,
        S,
        E: From<ConductionThermalError>,
    >(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        problem: &SpatialConduction<C>,
        boundary: &SpatialConductionBoundary,
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
                if remaining == 0 {
                    return Err(CircuitEnthalpyError::CircuitBudget);
                }
                let field = next
                    .solve_with_motion(boundary, motion, field_controls, &mut || {
                        (*cancel.borrow_mut())()
                    })
                    .map_err(field_error)
                    .map_err(map_conduction_error)?;
                let heat_j = finite(dt * field.joule_power_w())?;
                let terminal_work_j = finite(dt * field.terminal_power_w())?;
                let mechanical_work_j = finite(dt * field.mechanical_power_w())?;
                let net_supplied_work_j = finite(terminal_work_j - mechanical_work_j)?;
                if (field.joule_power_w() > 0.0 && heat_j == 0.0)
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
                            contacts,
                            heat_j,
                            terminal_work_j,
                            mechanical_work_j,
                            net_supplied_work_j,
                        },
                    },
                    1,
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
