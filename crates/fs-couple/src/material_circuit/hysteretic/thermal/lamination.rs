//! Physical thin-lamination loss on the original core port and heat inventory.
use super::*;
use crate::{
    circuit_enthalpy::{ResistorCell, material::CellResistor},
    material_circuit::MaterialResistor,
};
use fs_material::{
    conductor::ConductorError, field_storage::hysteresis::lamination::LaminatedCoreEddyCurrent,
    state_point::MaterialPropertySelection,
};
use fs_phs::circuit::Branch;

/// One explicit native resistor in parallel with the ORIGINAL magnetic port.
/// Its source and equivalent transport coordinates come only from the physical
/// lamination reduction. It must receive heat in that same whole core cell.
#[derive(Clone, Debug, PartialEq)]
pub struct HystereticLamination {
    pub branch: usize,
    pub loss: LaminatedCoreEddyCurrent,
}

fn bad(quantity: &'static str) -> ConductorError {
    ConductorError::InvalidInput { quantity }
}

impl HystereticLamination {
    /// Build the original material-network binding for the physical core.
    /// L_e and A_e are derived transport coordinates with L_e*A_e=core volume;
    /// they are neither a second inventory nor actual wire/contact geometry.
    /// Keep the corresponding basis from `loss.equivalent_conductor()` in the
    /// material network's electrical basis plan.
    pub fn resistor_binding(&self, core_cell: CellId) -> Result<CellResistor, ConductorError> {
        let conductor = self.loss.equivalent_conductor();
        let basis = conductor.source_basis();
        if conductor.source_key().name() != basis.property_name()
            || conductor.source_key().quantity()
                != fs_qty::QuantitySpec::dimensional(basis.dimensions())
        {
            return Err(bad(
                "native caloric lamination requires its exact canonical electrical basis",
            ));
        }
        Ok(CellResistor {
            destination: ResistorCell {
                branch: self.branch,
                cell: core_cell,
            },
            selection: MaterialPropertySelection::PinnedByProperty(vec![(
                basis.property_name().to_owned(),
                conductor.source_property().answer().receipt.selected,
            )]),
            length_m: conductor.length_m(),
            area_m2: conductor.area_m2(),
        })
    }

    pub(super) fn refresh(
        &self,
        core: &UniformHystereticCore,
        resistors: &[MaterialResistor],
    ) -> Result<Self, ConductorError> {
        if !self.loss.matches_core(core) {
            return Err(bad(
                "lamination must retain its original magnetic core and geometry",
            ));
        }
        let mut matching = resistors.iter().filter(|r| r.branch == self.branch);
        let conductor = &matching
            .next()
            .ok_or_else(|| bad("lamination shunt is missing its native sourced resistor"))?
            .conductor;
        if matching.next().is_some() {
            return Err(bad("one native shunt must own classical lamination loss"));
        }
        if let Some(source) = core.source() {
            if source.material().card_identity() != conductor.material().card_identity()
                || source.material().query() != conductor.material().query()
            {
                return Err(bad(
                    "lamination electrical and magnetic sources require the same actual core state",
                ));
            }
        }
        Ok(Self {
            branch: self.branch,
            loss: self.loss.with_equivalent_conductor(conductor)?,
        })
    }

    fn admit_binding(
        &self,
        circuit: &HystereticCircuit,
        network: &MaterialNetwork<'_>,
        bases: &[(usize, ElectricalPropertyBasis)],
        core_cell: CellId,
    ) -> Result<(), ConductorError> {
        let graph = circuit.native_circuit().graph();
        let &(a, b, _) = graph
            .branches
            .get(circuit.core_branch)
            .ok_or_else(|| bad("original hysteretic core port is missing"))?;
        let &(x, y, branch) = graph
            .branches
            .get(self.branch)
            .ok_or_else(|| bad("native lamination shunt branch is missing"))?;
        if !matches!(branch, Branch::Resistor { .. }) || !((a == x && b == y) || (a == y && b == x))
        {
            return Err(bad(
                "lamination shunt must be parallel to the original single-core port",
            ));
        }
        let expected = self.resistor_binding(core_cell)?;
        let actual = network
            .resistors
            .iter()
            .find(|r| r.destination.branch == self.branch)
            .ok_or_else(|| bad("lamination requires its original core heat/source binding"))?;
        let core = network
            .cells
            .iter()
            .find(|c| c.cell == core_cell)
            .ok_or_else(|| bad("lamination core caloric/material binding is absent"))?;
        if actual.destination.cell != core_cell
            || actual.length_m.to_bits() != expected.length_m.to_bits()
            || actual.area_m2.to_bits() != expected.area_m2.to_bits()
            || core.point.material.content_hash()
                != self.loss.equivalent_conductor().material().card_identity()
            || bases
                .iter()
                .filter(|(branch, _)| *branch == self.branch)
                .count()
                != 1
            || !bases.contains(&(self.branch, self.loss.equivalent_conductor().source_basis()))
        {
            return Err(bad(
                "lamination must retain its reduced transport, core receiver and source basis",
            ));
        }
        self.refresh(circuit.core(), circuit.native_circuit().resistors())?;
        Ok(())
    }
}

impl<'n, 'a> ThermalHystereticCircuit<'n, 'a> {
    /// Add classical eddy reaction and heat from an explicit insulated
    /// lamination to the original single-core magnetic/thermal solve.
    ///
    /// The graph must already contain the shunt resistor in parallel with the
    /// magnetic port. Build its `CellResistor` with `resistor_binding`; do not
    /// author a fitted core-loss resistance. The existing Kirchhoff owner then
    /// solves the added reaction current, and the existing enthalpy owner
    /// refreshes its original rho(T) or sigma(T) source on every trial.
    ///
    /// The shunt's native I²R work and the original irreversible pinning work
    /// accumulate ONCE in the same physical core. Their receipts remain
    /// separate in the existing circuit-resistor and core reports. Ordinary
    /// winding losses retain distinct whole-cell material inventories.
    ///
    /// This opt-in uniform-flux, low-frequency, fixed insulated-sheet reduction
    /// does not resolve skin/edge currents, inter-sheet leakage, excess losses,
    /// variable magnetic-fit coefficients or a second electromagnetic history.
    /// Do not add another eddy/core-loss law for the same volume. Source failure,
    /// cancellation and rejected coupled trials publish no partial history.
    #[allow(clippy::too_many_arguments)]
    pub fn bind_with_lamination(
        circuit: HystereticCircuit,
        thermal: &DepositionSimulation<'_>,
        network: &'n MaterialNetwork<'a>,
        electrical_bases: &[(usize, ElectricalPropertyBasis)],
        core_cell: CellId,
        core_temperature_interval_k: [f64; 2],
        lamination: HystereticLamination,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        Self::bind_with_lamination_and_link_conductances(
            circuit,
            thermal,
            network,
            electrical_bases,
            core_cell,
            core_temperature_interval_k,
            lamination,
            &BTreeMap::new(),
            cancelled,
        )
    }

    /// Same physical lamination with additional original finite coolant links.
    /// Continue through `step_candidate_with_heat_assembly` so their sources
    /// are freshly resolved at every core/winding/coolant trial endpoint.
    #[allow(clippy::too_many_arguments)]
    pub fn bind_with_lamination_and_link_conductances(
        circuit: HystereticCircuit,
        thermal: &DepositionSimulation<'_>,
        network: &'n MaterialNetwork<'a>,
        electrical_bases: &[(usize, ElectricalPropertyBasis)],
        core_cell: CellId,
        core_temperature_interval_k: [f64; 2],
        lamination: HystereticLamination,
        initial_link_conductances_w_k: &BTreeMap<u64, f64>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        if cancelled() {
            return Err(CircuitEnthalpyError::Cancelled);
        }
        lamination
            .admit_binding(&circuit, network, electrical_bases, core_cell)
            .map_err(lamination_error)?;
        Self::bind_with_lamination_option(
            circuit,
            thermal,
            network,
            electrical_bases,
            core_cell,
            core_temperature_interval_k,
            initial_link_conductances_w_k,
            Some(lamination),
            cancelled,
        )
    }
}
