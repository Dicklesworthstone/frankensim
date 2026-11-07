//! Original Cp(T) through native hydraulic graphs and their mixing junctions.
//!
//! One immutable, fixed-composition isobar owns every supply, branch and node.
//! Junctions mix mass-weighted sensible enthalpy, then recover temperature with
//! that same source's inverse. Branch exchange uses `CaloricAirPath`; the
//! native graph traversal, row ownership, cancellation and global watt gate
//! remain shared with `TransportNetwork`.
//!
//! Density, hydraulic flow and film coefficients are prescribed. The source
//! pressure is the common thermodynamic approximation, not a reinterpretation
//! of the graph's hydraulic pressure drops. No pressure work, phase change,
//! composition mixing or density/flow feedback is inferred. Imposed return-air
//! links close a bounded nonlinear enthalpy balance using the original source
//! response. Closed native circulations use simultaneous source-enthalpy balances;
//! hydraulic feedback and return-loop adjoints are not supplied.

use std::sync::Arc;

use fs_material::thermal::transport::ThermalTransportError;
use fs_material::thermal::transport::isobaric::IsobaricHeatCapacityLaw;

use super::super::coupled_transport::{
    CoupledTransportModel, CoupledTransportSolution, solve_driver,
};
use super::*;
use crate::conjugate::caloric::CaloricAirPath;
use crate::conjugate::{ConjugateConfig, IqnIlsConfig, RegionBalance, SolidRegionState};

mod recirculation;
pub use recirculation::CaloricRecirculationConfig;
mod circulation;
pub use circulation::CaloricCirculationConfig;
mod recovery;
pub use recovery::{CaloricHeatRecoveryConfig, CaloricHeatRecoveryLink,
    CaloricHeatRecoveryReport, CaloricRecoveredStream};

/// Source-caloric specialization of the existing admitted hydraulic network.
/// The source is shared, never selected independently at a junction or branch.
#[derive(Debug)]
pub struct CaloricTransportNetwork<'a> {
    topology: TransportNetwork<'a>,
    caloric: Arc<IsobaricHeatCapacityLaw>,
    reference_temperature: Temperature,
    external_mass: Vec<f64>,
    mass_residual: Vec<f64>,
    feedback: Option<recirculation::CaloricFeedback>,
    circulation: Option<circulation::CaloricCirculation>,
    recovery: Option<recovery::CaloricRecovery>,
}

impl<'a> CaloricTransportNetwork<'a> {
    /// Bind the original source to a solved, directed-acyclic hydraulic graph.
    ///
    /// The physical reference temperature determines the Cp accompanying the
    /// fixed density at native admission (including gravity-snapshot checks).
    /// It is independent of the caloric datum, which may change without
    /// changing mixing or heat transfer. Actual source support is checked at
    /// every external inlet and every subsequent fluid temperature.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        cx: &Cx<'_>,
        flow: &'a GraphSolution,
        reference_density: Density,
        reference_temperature: Temperature,
        caloric: Arc<IsobaricHeatCapacityLaw>,
        models: Vec<BranchThermalModel>,
        inlets: &[TransportInlet],
        config: TransportConfig,
    ) -> Result<Self, TransportError> {
        Self::new_inner(cx, flow, reference_density, reference_temperature, caloric,
            models, inlets, config, None)
    }

    #[allow(clippy::too_many_arguments)]
    fn new_inner(
        cx: &Cx<'_>,
        flow: &'a GraphSolution,
        reference_density: Density,
        reference_temperature: Temperature,
        caloric: Arc<IsobaricHeatCapacityLaw>,
        models: Vec<BranchThermalModel>,
        inlets: &[TransportInlet],
        config: TransportConfig,
        circulation: Option<CaloricCirculationConfig>,
    ) -> Result<Self, TransportError> {
        poll(cx)?;
        let cp = caloric
            .specific_heat_cp_j_kg_k(reference_temperature.value())
            .map_err(source_error)?;
        let air = TransportAir {
            density: reference_density,
            specific_heat_j_kg_k: cp,
        };
        let mut topology = if circulation.is_some() {
            TransportNetwork::new_circulating(cx, flow, air, models, reference_temperature, config)?
        } else {
            TransportNetwork::new(cx, flow, air, models, inlets, config)?
        };
        for inlet in inlets {
            poll(cx)?;
            caloric
                .specific_heat_cp_j_kg_k(inlet.temperature.value())
                .map_err(source_error)?;
        }
        // A wholly stagnant adiabatic graph has no inlet-derived reference.
        if topology.reference_k == 0.0 {
            topology.reference_k = reference_temperature.value();
        }
        // Recompute MASS balances from the admitted original branch rates.
        // Dividing rounded capacity differences by a frozen Cp would make
        // these balances depend on the choice of reference temperature.
        let count = flow.pressures.len();
        let mut incoming = vec![0.0; count];
        let mut outgoing = vec![0.0; count];
        for (node, branches) in topology.outgoing.iter().enumerate() {
            for &branch in branches {
                poll(cx)?;
                let mass = topology.mass_flows[branch];
                outgoing[node] = finite(outgoing[node] + mass, "outgoing source mass flow")?;
                let down = topology.downstream[branch];
                incoming[down] = finite(incoming[down] + mass, "incoming source mass flow")?;
            }
        }
        let fixed: BTreeSet<_> = flow
            .boundaries
            .iter()
            .map(|boundary| boundary.node)
            .collect();
        let mut external_mass = vec![0.0; count];
        let mut mass_residual = vec![0.0; count];
        for node in 0..count {
            poll(cx)?;
            let residual = finite(outgoing[node] - incoming[node], "source mass residual")?;
            if fixed.contains(&node) && circulation.is_none() {
                external_mass[node] = residual;
            } else {
                mass_residual[node] = residual;
            }
            if (external_mass[node] > 0.0) != topology.inlet_temperatures[node].is_some() {
                return Err(TransportError::Node {
                    node,
                    reason: "source mass balance and external inlet ownership disagree",
                });
            }
        }
        let mut network = Self {
            topology,
            caloric,
            reference_temperature,
            external_mass,
            mass_residual,
            feedback: None,
            circulation: None,
            recovery: None,
        };
        if let Some(config) = circulation {
            network.circulation = Some(circulation::CaloricCirculation::new(&network, cx, config)?);
        }
        Ok(network)
    }

    /// Exact original Cp source, including its pressure, identity and datum.
    #[must_use]
    pub fn caloric_law(&self) -> &IsobaricHeatCapacityLaw {
        &self.caloric
    }

    /// Fixed density converting native volume flow to transported mass flow.
    #[must_use]
    pub const fn reference_density(&self) -> Density {
        self.topology.air.density
    }

    /// Physical reference of the native density/Cp snapshot, not an energy datum.
    #[must_use]
    pub const fn reference_temperature(&self) -> Temperature {
        self.reference_temperature
    }

    /// The unchanged hydraulic solution and its original branch ordering.
    #[must_use]
    pub const fn hydraulics(&self) -> &GraphSolution {
        self.topology.hydraulics()
    }

    /// Solid rows in branch-major, actual flow-direction order.
    #[must_use]
    pub fn regions(&self) -> Vec<&str> {
        self.topology.regions()
    }

    /// Zero-wall-heat initialization with actual enthalpy mixing of all supplies.
    /// A closed circulation instead uses its declared reference as an initial
    /// solid guess; there is no unique adiabatic closed-loop temperature.
    pub fn initial_references(&self, cx: &Cx<'_>) -> Result<Vec<f64>, TransportError> {
        if self.circulation.is_some() {
            return self.topology.initial_references(cx);
        }
        Ok(self.march_inner(cx, None)?.reference_temperatures_k)
    }

    /// March all branches with source enthalpy and the native global watt gate.
    /// Failed source admission, inversion or balance returns no partial field.
    pub fn march(&self, cx: &Cx<'_>, walls_k: &[f64]) -> Result<TransportMarch, TransportError> {
        if walls_k.len() != self.topology.offsets.last().copied().unwrap_or(0) {
            return Err(TransportError::InvalidInput(
                "one wall temperature per thermal region is required",
            ));
        }
        for &wall in walls_k {
            poll(cx)?;
            positive(wall, "wall temperature")?;
        }
        self.march_inner(cx, Some(walls_k))
    }

    fn march_inner(
        &self,
        cx: &Cx<'_>,
        walls: Option<&[f64]>,
    ) -> Result<TransportMarch, TransportError> {
        if let Some(circulation) = &self.circulation {
            return circulation.march(self, cx, walls.ok_or(TransportError::InvalidInput(
                "closed source-caloric circulation needs actual wall temperatures",
            ))?);
        }
        if let Some(recovery) = &self.recovery {
            return recovery.march(self, cx, walls);
        }
        if let Some(feedback) = &self.feedback {
            return feedback.march(self, cx, walls);
        }
        self.topology.march_once_with_source(
            cx,
            walls,
            &self.topology.inlet_temperatures,
            None,
            Some(self.context()),
            None,
        )
    }

    /// Eliminate the native DAG for one external-supply enthalpy direction.
    /// Shared by mass returns and finite thermal recovery; arithmetic and
    /// actual source branch gains retain the original return owner's order.
    fn inlet_enthalpy_response(
        &self, cx: &Cx<'_>, source: usize, gains: &[f64],
    ) -> Result<Vec<f64>, TransportError> {
        let topology = &self.topology;
        if source >= self.external_mass.len() || gains.len() != topology.models.len() {
            return Err(TransportError::InvalidInput("source inlet response shape"));
        }
        let mut mass: Vec<_> = self.external_mass.iter().map(|flow| flow.max(0.0)).collect();
        let mut values = vec![0.0; mass.len()];
        values[source] = 1.0;
        for &node in &topology.order {
            poll(cx)?;
            for &branch in &topology.outgoing[node] {
                poll(cx)?;
                let down = topology.downstream[branch];
                let total = positive(mass[down]+topology.mass_flows[branch], "response mixing mass")?;
                let arriving = values[node]*gains[branch];
                values[down] = finite(values[down]
                    +(arriving-values[down])*(topology.mass_flows[branch]/total),
                    "return enthalpy response")?;
                mass[down] = total;
            }
        }
        Ok(values)
    }

    fn context(&self) -> CaloricContext<'_> {
        CaloricContext {
            law: &self.caloric,
            external_mass: &self.external_mass,
            mass_residual: &self.mass_residual,
            reference_k: None,
            inlet_enthalpies: None,
        }
    }

    fn path(&self, branch: usize, inlet: f64) -> Result<CaloricAirPath, TransportError> {
        self.context().path(&self.topology, branch, inlet)
    }

    /// Run the existing shared-solid driver with source-caloric branches and
    /// mixing. Each iteration calls the solid once for every native region.
    pub fn solve<F>(
        &self,
        cx: &Cx<'_>,
        config: &ConjugateConfig,
        solid: F,
    ) -> Result<CoupledTransportSolution, TransportError>
    where
        F: FnMut(&Cx<'_>, &[f64]) -> Result<Vec<SolidRegionState>, AirflowError>,
    {
        let initial = self.initial_references(cx)?;
        self.solve_from(cx, config, &initial, solid)
    }

    /// Continue from explicit native Robin references. Fixed relaxation has
    /// the existing replay semantics; scalar acceleration history restarts.
    pub fn solve_from<F>(
        &self,
        cx: &Cx<'_>,
        config: &ConjugateConfig,
        initial_references_k: &[f64],
        solid: F,
    ) -> Result<CoupledTransportSolution, TransportError>
    where
        F: FnMut(&Cx<'_>, &[f64]) -> Result<Vec<SolidRegionState>, AirflowError>,
    {
        solve_driver(cx, self, config, initial_references_k, None, solid)
    }

    /// Use the native bounded IQN-ILS history across the complete source-caloric
    /// network. Temperature and each branch's original watt gate still apply.
    pub fn solve_iqn<F>(
        &self,
        cx: &Cx<'_>,
        config: &ConjugateConfig,
        acceleration: IqnIlsConfig,
        solid: F,
    ) -> Result<CoupledTransportSolution, TransportError>
    where
        F: FnMut(&Cx<'_>, &[f64]) -> Result<Vec<SolidRegionState>, AirflowError>,
    {
        let initial = self.initial_references(cx)?;
        self.solve_iqn_from(cx, config, &initial, acceleration, solid)
    }

    /// Warm-start the same native IQN driver with empty secant history.
    pub fn solve_iqn_from<F>(
        &self,
        cx: &Cx<'_>,
        config: &ConjugateConfig,
        initial_references_k: &[f64],
        acceleration: IqnIlsConfig,
        solid: F,
    ) -> Result<CoupledTransportSolution, TransportError>
    where
        F: FnMut(&Cx<'_>, &[f64]) -> Result<Vec<SolidRegionState>, AirflowError>,
    {
        solve_driver(
            cx,
            self,
            config,
            initial_references_k,
            Some(acceleration),
            solid,
        )
    }
}

fn source_error(error: ThermalTransportError) -> TransportError {
    AirflowError::CaloricSource(error).into()
}

/// Optional primal output retained by the nonlinear return owner. Branch gains
/// come from the same original channel pass as its heat and node enthalpies.
pub(super) struct CaloricCapture<'a> {
    pub(super) node_enthalpies: &'a mut [Option<f64>],
    pub(super) branch_enthalpy_gains: &'a mut [f64],
    // Private nonlinear closed-loop trials retain their physical mixer watts
    // before final native admission. Published marches always run both gates.
    pub(super) circulation_residuals_w: Option<&'a mut [f64]>,
}

/// Borrowed constitutive data for the native graph walk; topology stays native.
#[derive(Clone, Copy)]
pub(super) struct CaloricContext<'a> {
    law: &'a Arc<IsobaricHeatCapacityLaw>,
    pub(super) external_mass: &'a [f64],
    pub(super) mass_residual: &'a [f64],
    // Only the nonlinear return owner supplies primary energy coordinates.
    // Ordinary source marches retain the original allocation-free context.
    pub(super) reference_k: Option<f64>,
    pub(super) inlet_enthalpies: Option<&'a [Option<f64>]>,
}

impl CaloricContext<'_> {
    pub(super) fn inlet_enthalpy(
        self,
        node: usize,
        reference_k: f64,
        temperature: f64,
    ) -> Result<f64, TransportError> {
        match self.inlet_enthalpies {
            Some(values) => values[node].ok_or(TransportError::InvalidInput(
                "source supply is missing its retained enthalpy",
            )),
            None => self.enthalpy(reference_k, temperature),
        }
    }

    pub(super) fn enthalpy(
        self,
        reference_k: f64,
        temperature: f64,
    ) -> Result<f64, TransportError> {
        self.law
            .specific_enthalpy_change_j_kg(reference_k, temperature)
            .map_err(source_error)
    }

    pub(super) fn temperature(
        self,
        reference_k: f64,
        specific_h: f64,
    ) -> Result<f64, TransportError> {
        self.law
            .temperature_after_specific_heat(reference_k, specific_h)
            .map_err(source_error)
    }

    pub(super) fn path(
        self,
        topology: &TransportNetwork<'_>,
        branch: usize,
        inlet: f64,
    ) -> Result<CaloricAirPath, TransportError> {
        let BranchThermalModel::Exchange(segments) = &topology.models[branch] else {
            return Err(TransportError::InvalidInput(
                "an adiabatic branch has no caloric air path",
            ));
        };
        Ok(CaloricAirPath::new(
            inlet,
            topology.mass_flows[branch],
            self.law.clone(),
            segments.clone(),
        )?)
    }
}

impl CoupledTransportModel for CaloricTransportNetwork<'_> {
    fn topology(&self) -> &TransportNetwork<'_> {
        &self.topology
    }

    fn initial_references(&self, cx: &Cx<'_>) -> Result<Vec<f64>, TransportError> {
        CaloricTransportNetwork::initial_references(self, cx)
    }

    fn march(&self, cx: &Cx<'_>, walls: &[f64]) -> Result<TransportMarch, TransportError> {
        CaloricTransportNetwork::march(self, cx, walls)
    }

    fn admission_path(&self, branch: usize, inlet: f64) -> Result<AirPath, TransportError> {
        let cp = self
            .caloric
            .specific_heat_cp_j_kg_k(inlet)
            .map_err(source_error)?;
        let BranchThermalModel::Exchange(segments) = &self.topology.models[branch] else {
            return Err(TransportError::InvalidInput(
                "an adiabatic branch has no caloric air path",
            ));
        };
        Ok(AirPath::new(
            inlet,
            self.topology.mass_flows[branch],
            cp,
            segments.clone(),
        )?)
    }

    fn probe(
        &self,
        cx: &Cx<'_>,
        branch: usize,
        path: &AirPath,
        config: &ConjugateConfig,
        references: &[f64],
        states: &[SolidRegionState],
    ) -> Result<Vec<RegionBalance>, TransportError> {
        Ok(self
            .path(branch, path.inlet_temperature_k())?
            .solve_from(cx, config, references, |_, _| Ok(states.to_vec()))?
            .balance
            .regions)
    }
}

#[cfg(test)]
mod tests;
