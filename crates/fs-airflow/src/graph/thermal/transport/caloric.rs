//! Original Cp(T) through the native hydraulic DAG and its mixing junctions.
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
//! composition mixing or density/flow feedback is inferred. Native return-air
//! and closed-circulation solvers assume affine temperature response and are
//! unavailable on this distinct source-caloric type, as are their adjoints.

use std::sync::Arc;

use fs_material::thermal::transport::ThermalTransportError;
use fs_material::thermal::transport::isobaric::IsobaricHeatCapacityLaw;

use super::super::coupled_transport::{
    CoupledTransportModel, CoupledTransportSolution, solve_driver,
};
use super::*;
use crate::conjugate::caloric::CaloricAirPath;
use crate::conjugate::{ConjugateConfig, IqnIlsConfig, RegionBalance, SolidRegionState};

/// Source-caloric specialization of the existing admitted once-through network.
/// The source is shared, never selected independently at a junction or branch.
#[derive(Debug)]
pub struct CaloricTransportNetwork<'a> {
    topology: TransportNetwork<'a>,
    caloric: Arc<IsobaricHeatCapacityLaw>,
    reference_temperature: Temperature,
    external_mass: Vec<f64>,
    mass_residual: Vec<f64>,
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
        poll(cx)?;
        let cp = caloric
            .specific_heat_cp_j_kg_k(reference_temperature.value())
            .map_err(source_error)?;
        let mut topology = TransportNetwork::new(
            cx,
            flow,
            TransportAir {
                density: reference_density,
                specific_heat_j_kg_k: cp,
            },
            models,
            inlets,
            config,
        )?;
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
            if fixed.contains(&node) {
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
        Ok(Self {
            topology,
            caloric,
            reference_temperature,
            external_mass,
            mass_residual,
        })
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
    pub fn initial_references(&self, cx: &Cx<'_>) -> Result<Vec<f64>, TransportError> {
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
        self.topology.march_once_with_source(
            cx,
            walls,
            &self.topology.inlet_temperatures,
            None,
            Some(self.context()),
        )
    }

    fn context(&self) -> CaloricContext<'_> {
        CaloricContext {
            law: &self.caloric,
            external_mass: &self.external_mass,
            mass_residual: &self.mass_residual,
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

/// Borrowed constitutive data for the native graph walk; topology stays native.
#[derive(Clone, Copy)]
pub(super) struct CaloricContext<'a> {
    law: &'a Arc<IsobaricHeatCapacityLaw>,
    pub(super) external_mass: &'a [f64],
    pub(super) mass_residual: &'a [f64],
}

impl CaloricContext<'_> {
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
