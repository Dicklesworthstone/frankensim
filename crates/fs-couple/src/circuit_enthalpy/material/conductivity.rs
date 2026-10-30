//! Data-bound electrical operating response on the original thermal network.
//!
//! Cards, source selections, conductor geometry, thermal laws and heat
//! destinations come from the existing MaterialNetwork. These bindings only
//! declare which electrical property and operating coordinates the original
//! conductivity solver must close; no additional material or state owner exists.

use super::*;
use crate::circuit_enthalpy::conductivity::{
    self as coupled_conductivity, ConductivityEnthalpyCandidate,
};
use crate::material_circuit::conductivity::{
    ConductionAxis, ConductionCoordinate, ConductivitySample, ConductivitySolve, ConductorResponse,
};
use crate::circuit_enthalpy::motion::{FieldMotionEnthalpyInput, conductivity as motion_conductivity};
use crate::material_circuit::motion::MaterialFieldMotion;

/// Electrical response of an existing whole-cell resistor.
/// Every original resistor needs exactly one entry, including constant laws.
#[derive(Clone, Copy, Debug)]
pub struct CellConductionResponse<'a> {
    /// Original graph branch index from CellResistor::destination.
    pub branch: usize,
    /// Explicit source property; no conductivity/resistivity fallback is used.
    pub basis: ElectricalPropertyBasis,
    /// Existing typed field/current coordinates. At most one electric-field
    /// and one current-density role; an empty slice declares neither role.
    /// Signed versus magnitude conventions stay with each ConductionAxis.
    pub axes: &'a [ConductionAxis<'a>],
}

pub(super) type ElectricalAxes = BTreeMap<usize, Vec<(String, QuantitySpec)>>;

pub(super) fn admitted_axes(
    axes: Option<&ElectricalAxes>,
    branch: usize,
) -> Vec<(&str, QuantitySpec)> {
    axes.and_then(|bindings| bindings.get(&branch))
        .map(|axes| {
            axes.iter()
                .map(|(name, quantity)| (name.as_str(), *quantity))
                .collect()
        })
        .unwrap_or_default()
}

struct ResponsePlan<'a> {
    bases: BTreeMap<usize, ElectricalPropertyBasis>,
    axes: ElectricalAxes,
    sources: Vec<ConductorResponse<'a>>,
}

fn check_cancelled(cancelled: &mut impl FnMut() -> bool) -> Result<(), MaterialNetworkError> {
    if cancelled() {
        Err(DepositionError::Cancelled.into())
    } else {
        Ok(())
    }
}

fn response_plan<'a>(
    network: &'a MaterialNetwork<'_>,
    responses: &'a [CellConductionResponse<'_>],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<ResponsePlan<'a>, MaterialNetworkError> {
    check_cancelled(cancelled)?;
    let invalid = MaterialNetworkError::Invalid;
    if network.resistors.is_empty()
        || network.resistors.len() > MAX_MATERIAL_CIRCUIT_BRANCHES
        || responses.len() != network.resistors.len()
    {
        return Err(invalid(
            "one response for every bounded original resistor is required",
        ));
    }
    let choices: Vec<_> = responses
        .iter()
        .map(|response| (response.branch, response.basis))
        .collect();
    let bases = network.basis_plan(&choices)?;
    let mut cells = BTreeMap::new();
    for cell in network.cells {
        check_cancelled(cancelled)?;
        if cells.insert(cell.cell, cell).is_some() {
            return Err(invalid("duplicate material cell binding"));
        }
    }
    let mut resistors = BTreeMap::new();
    for resistor in network.resistors {
        check_cancelled(cancelled)?;
        if resistors
            .insert(resistor.destination.branch, resistor)
            .is_some()
        {
            return Err(invalid("duplicate whole-cell resistor branch binding"));
        }
    }
    let mut axes = BTreeMap::new();
    let mut sources = Vec::with_capacity(responses.len());
    for response in responses {
        check_cancelled(cancelled)?;
        let resistor = resistors
            .get(&response.branch)
            .ok_or(invalid("response branch lacks its original resistor"))?;
        let cell = cells
            .get(&resistor.destination.cell)
            .ok_or(invalid("resistor cell lacks material binding"))?;
        let context = cell.point.context;
        if response.axes.len() > 2 {
            return Err(invalid(
                "electrical response permits at most one field and one current-density axis",
            ));
        }
        for (index, axis) in response.axes.iter().enumerate() {
            if !context.axes().contains_key(axis.name)
                || context.axis_quantities().get(axis.name)
                    != Some(&QuantitySpec::dimensional(axis.coordinate.dimensions()))
                || response.axes[..index].iter().any(|other| {
                    other.name == axis.name
                        || other.coordinate.dimensions() == axis.coordinate.dimensions()
                })
            {
                return Err(invalid(
                    "conductivity axes require unique existing SI field/current-density coordinates",
                ));
            }
        }
        for (name, quantity) in context.axis_quantities() {
            if [
                ConductionCoordinate::ElectricField.dimensions(),
                ConductionCoordinate::CurrentDensity.dimensions(),
            ]
            .contains(&quantity.dims())
                && !response.axes.iter().any(|axis| axis.name == name.as_str())
            {
                return Err(invalid(
                    "an electrical operating coordinate cannot remain frozen",
                ));
            }
        }
        axes.insert(
            response.branch,
            response
                .axes
                .iter()
                .map(|axis| {
                    (
                        axis.name.to_owned(),
                        QuantitySpec::dimensional(axis.coordinate.dimensions()),
                    )
                })
                .collect(),
        );
        sources.push(ConductorResponse {
            branch: response.branch,
            card: cell.point.material,
            axes: response.axes,
        });
    }
    check_cancelled(cancelled)?;
    Ok(ResponsePlan {
        bases,
        axes,
        sources,
    })
}

impl<'a> MaterialNetwork<'a> {
    fn resolve_conductivity_transport(
        &self,
        thermal: &DepositionSimulation<'_>,
        electrical: &[ConductivitySample],
        coefficients: StepCoefficients<MaterialNetworkSample>,
        directional: Option<&directional::DirectionalPlan>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<StepCoefficients<MaterialNetworkSample>, MaterialNetworkError> {
        check_cancelled(&mut cancelled)?;
        let invalid = MaterialNetworkError::Invalid;
        let evaluated: BTreeMap<_, _> = electrical
            .iter()
            .map(|sample| (sample.point.branch, sample))
            .collect();
        if evaluated.len() != electrical.len() || evaluated.len() != self.resistors.len()
        {
            return Err(invalid(
                "exactly one solved-field sample per resistor is required",
            ));
        }
        let bindings: BTreeMap<_, _> = self.cells.iter().map(|cell| (cell.cell, cell)).collect();
        let mut samples: BTreeMap<_, _> = coefficients
            .payload
            .cells
            .into_iter()
            .map(|sample| (sample.cell, sample))
            .collect();
        for resistor in self.resistors {
            check_cancelled(&mut cancelled)?;
            let source = evaluated
                .get(&resistor.destination.branch)
                .ok_or(invalid("resistor has no solved-field material sample"))?;
            let sample = samples
                .get_mut(&resistor.destination.cell)
                .ok_or(invalid("resistor has no thermal material sample"))?;
            if let Some(material) = &mut sample.conductivity {
                // The electrical owner refreshed the complete pinned bundle at
                // E=V/L and J=I/A, including its thermal conductivity. Reuse it;
                // a second query could select a different source or time level.
                *material = source.conductor.material().clone();
            }
        }
        let (conductances_w_k, links) =
            self.resolve_links(thermal, &bindings, &samples, directional, &mut || {
                if cancelled() {
                    Err(DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            })?;
        Ok(StepCoefficients {
            sources_w: coefficients.sources_w,
            conductances_w_k,
            payload: MaterialNetworkSample {
                cells: samples.into_values().collect(),
                links,
            },
        })
    }

    /// Pin the original complete thermal/electrical source bundle while admitting
    /// declared electrical operating axes throughout each resistor's bundle.
    ///
    /// Use resolve_conductivity_seeds to construct the original MaterialCircuit;
    /// its resistor bundles must include these same thermal properties and pins.
    /// The ordinary prepare/resolve routes retain their static thermal admission.
    /// Thermal conductivity may depend on E/J only where an original resistor
    /// supplies those coordinates. Its solved-field value drives the same heat
    /// residual; passive thermal cells cannot invent an electrical field owner.
    pub fn prepare_conductivity<'n>(
        &'n self,
        thermal: &DepositionSimulation<'_>,
        responses: &[CellConductionResponse<'_>],
        cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedMaterialNetwork<'n, 'a>, MaterialNetworkError> {
        self.prepare_conductivity_with_directional_conduction(thermal, responses, &[], cancelled)
    }

    /// Bind principal conductivity roles to the same live E/J material bundles
    /// used by the original electrical response solver. All three principal
    /// claims are re-evaluated at solved field/current and trial thermal state;
    /// only the declared one-dimensional gradient projection is added.
    pub fn prepare_conductivity_with_directional_conduction<'n>(
        &'n self,
        thermal: &DepositionSimulation<'_>,
        responses: &[CellConductionResponse<'_>],
        links: &[DirectionalConductionLink],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedMaterialNetwork<'n, 'a>, MaterialNetworkError> {
        let plan = response_plan(self, responses, &mut cancelled)?;
        let directional = if links.is_empty() { None } else {
            Some(directional::DirectionalPlan::new(self, thermal, links, &mut cancelled)?)
        };
        self.prepare_bound_with_links(thermal, plan.bases, Some(&plan.axes), None, directional, cancelled)
    }
}

impl PreparedMaterialNetwork<'_, '_> {
    /// Resolve complete pinned bundles at the supplied thermal state and the
    /// original context's seed E/J coordinates, without advancing either history.
    ///
    /// Build the initial MaterialCircuit from these resistors. Seed coordinates
    /// are starting values, not a claim that they satisfy Kirchhoff/constitutive
    /// equations. conductivity_candidate closes them at every thermal trial;
    /// both its electrical and thermal reports retain the actual solved-field
    /// source evaluations and the resulting heat conductances.
    /// The requested property bases must match this prepared network exactly.
    pub fn resolve_conductivity_seeds(
        &self,
        thermal: &DepositionSimulation<'_>,
        responses: &[CellConductionResponse<'_>],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CircuitThermalCoefficients<MaterialNetworkSample>, MaterialNetworkError> {
        let plan = response_plan(self.network, responses, &mut cancelled)?;
        if plan.bases != self.bases {
            return Err(MaterialNetworkError::Invalid(
                "conductivity response bases must retain the prepared source choices",
            ));
        }
        self.network.resolve_bound(
            thermal,
            &self.bases,
            Some(&self.pins),
            Some(&plan.axes),
            self.directional.as_ref(),
            cancelled,
        )
    }

    /// Joint conductivity, Joule heating and enthalpy advance using only the
    /// existing cell/link/resistor bindings plus declared operating axes.
    ///
    /// Every trial reuses the original source pins, whole-cell geometry and heat
    /// destinations. The existing electrical owner closes E=V/L and J=I/A at
    /// the proposed thermal state and the shared enthalpy owner receives its
    /// native Joule heat once. Trial histories remain private until acceptance.
    /// Thermal links and their samples use the complete original bundles at that
    /// solved field/current, so electrical operating response also drives heat
    /// transport. Seed E/J coordinates never stand in for the accepted fields.
    ///
    /// Source/domain refusal, work-budget exhaustion and cancellation publish
    /// neither circuit nor thermal state. This inherits the scalar, fixed-geometry
    /// and no-contact scope of the existing conductivity solver.
    #[allow(clippy::too_many_arguments)]
    pub fn conductivity_candidate<'m>(
        &self,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        circuit: &MaterialCircuit,
        ports: &[f64],
        responses: &[CellConductionResponse<'_>],
        controls: CircuitEnthalpySolve,
        conductivity: ConductivitySolve,
        guess_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<
        ConductivityEnthalpyCandidate<'m, MaterialNetworkSample>,
        CircuitEnthalpyError<MaterialNetworkError>,
    > {
        let cancel = RefCell::new(cancelled);
        let plan = match response_plan(self.network, responses, &mut || (*cancel.borrow_mut())()) {
            Err(MaterialNetworkError::Thermal(DepositionError::Cancelled)) => {
                return Err(CircuitEnthalpyError::Cancelled);
            }
            other => other.map_err(CircuitEnthalpyError::Source)?,
        };
        if plan.bases != self.bases {
            return Err(CircuitEnthalpyError::Source(MaterialNetworkError::Invalid(
                "conductivity response bases must retain the prepared source choices",
            )));
        }
        let destinations = self.resistor_destinations();
        let result = coupled_conductivity::step_candidate_with_transport(
            thermal,
            target_s,
            CircuitEnthalpyInput {
                circuit,
                ports,
                destinations: &destinations,
            },
            &plan.sources,
            controls,
            conductivity,
            guess_j,
            |state| {
                self.network.resolve_bound(
                    state,
                    &self.bases,
                    Some(&self.pins),
                    Some(&plan.axes),
                    self.directional.as_ref(),
                    || (*cancel.borrow_mut())(),
                )
            },
            |state, electrical, coefficients| {
                self.network
                    .resolve_conductivity_transport(state, &electrical.evaluated, coefficients,
                        self.directional.as_ref(), || (*cancel.borrow_mut())())
            },
            || (*cancel.borrow_mut())(),
        );
        match result {
            Err(CircuitEnthalpyError::Source(MaterialNetworkError::Thermal(
                DepositionError::Cancelled,
            ))) => Err(CircuitEnthalpyError::Cancelled),
            other => other,
        }
    }

    /// The same data-bound E/J and thermal transport closure for an actual
    /// electric/magnetic gap actuator. All source cards, pins, branch geometry,
    /// caloric laws and heat destinations come from this prepared network.
    /// Electrical operating axes come from `responses`; actual field motion
    /// comes from the original reciprocal descriptor, including optional support.
    /// The SAME solved-field bundle supplies resistor response and heat transport.
    ///
    /// Every trial restarts accepted electrical, body and thermal history. Only
    /// original native Joule work is deposited; finite nonlinear residuals remain
    /// separate. Stationary conductor geometry, fixed body/field material state
    /// and bulk-only no-contact restrictions remain explicit.
    #[allow(clippy::too_many_arguments)]
    pub fn field_motion_conductivity_candidate<'m>(
        &self, thermal: &DepositionSimulation<'m>, target_s: f64,
        motion: &MaterialFieldMotion, ports: &[f64], external_force_n: f64,
        responses: &[CellConductionResponse<'_>], controls: CircuitEnthalpySolve,
        conductivity: ConductivitySolve, guess_j: &[f64], cancelled: impl FnMut() -> bool,
    ) -> Result<motion_conductivity::FieldConductivityEnthalpyCandidate<'m, MaterialNetworkSample>,
        CircuitEnthalpyError<MaterialNetworkError>> {
        let cancel = RefCell::new(cancelled);
        let plan = match response_plan(self.network, responses, &mut || (*cancel.borrow_mut())()) {
            Err(MaterialNetworkError::Thermal(DepositionError::Cancelled)) =>
                return Err(CircuitEnthalpyError::Cancelled),
            other => other.map_err(CircuitEnthalpyError::Source)?,
        };
        if plan.bases != self.bases {
            return Err(CircuitEnthalpyError::Source(MaterialNetworkError::Invalid(
                "conductivity response bases must retain the prepared source choices")));
        }
        let destinations = self.resistor_destinations();
        let result = motion_conductivity::step_candidate_with_transport(thermal, target_s,
            FieldMotionEnthalpyInput { motion, ports, external_force_n, destinations: &destinations },
            &plan.sources, controls, conductivity, guess_j,
            |state| self.network.resolve_bound(state, &self.bases, Some(&self.pins),
                Some(&plan.axes), self.directional.as_ref(), || (*cancel.borrow_mut())()),
            |state, electrical, coefficients| self.network.resolve_conductivity_transport(
                state, &electrical.evaluated, coefficients, self.directional.as_ref(),
                || (*cancel.borrow_mut())()),
            || (*cancel.borrow_mut())());
        match result {
            Err(CircuitEnthalpyError::Source(MaterialNetworkError::Thermal(
                DepositionError::Cancelled))) => Err(CircuitEnthalpyError::Cancelled),
            other => other,
        }
    }
}
