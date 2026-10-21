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
    ConductionAxis, ConductionCoordinate, ConductivitySolve, ConductorResponse,
};

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
    /// Pin the original complete thermal/electrical source bundle while admitting
    /// the declared electrical operating axes on its electrical property only.
    ///
    /// Use resolve_conductivity_seeds to construct the original MaterialCircuit;
    /// its resistor bundles must include these same thermal properties and pins.
    /// The ordinary prepare/resolve routes retain their static thermal admission.
    /// Thermal conductivity depending on E/J is refused: heat conductances are
    /// resolved before the inner electrical iteration and cannot freeze its fields.
    pub fn prepare_conductivity<'n>(
        &'n self,
        thermal: &DepositionSimulation<'_>,
        responses: &[CellConductionResponse<'_>],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedMaterialNetwork<'n, 'a>, MaterialNetworkError> {
        let plan = response_plan(self, responses, &mut cancelled)?;
        self.prepare_bound(thermal, plan.bases, Some(&plan.axes), cancelled)
    }
}

impl PreparedMaterialNetwork<'_, '_> {
    /// Resolve complete pinned bundles at the supplied thermal state and the
    /// original context's seed E/J coordinates, without advancing either history.
    ///
    /// Build the initial MaterialCircuit from these resistors. Seed coordinates
    /// are starting values, not a claim that they satisfy Kirchhoff/constitutive
    /// equations. conductivity_candidate closes them at every thermal trial;
    /// its electrical report retains the actual solved-field source evaluations.
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
    /// Thermal network samples retain the seed E/J context; actual applied and
    /// evaluated electrical bundles are in the candidate's electrical report.
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
        let result = coupled_conductivity::step_candidate(
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
                    || (*cancel.borrow_mut())(),
                )
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
}
