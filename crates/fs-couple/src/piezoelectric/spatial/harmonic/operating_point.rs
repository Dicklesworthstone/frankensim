//! Steady electromechanical bias and harmonic perturbations on one original mesh.
//!
//! The affine, small-strain Helmholtz law has the same tangent at every admitted
//! bias. DC stress/polarization are retained in the equilibrium, never inserted
//! as phasors. This is not geometric stress stiffening, field-dependent material
//! tuning or temperature modulation. Temperature and reference geometry stay
//! fixed during each periodic solve. See the stationary/perturbation distinction:
//! <https://doc.comsol.com/6.4/doc/com.comsol.help.comsol/comsol_ref_solver.36.059.html>

use super::super::{
    PiezoelectricElectrode, Result, TetPiezoelectricControls, TetPiezoelectricModel,
    TetPiezoelectricSolution, invalid, poll,
};
use super::acoustic::HarmonicAcousticPort;
use super::envelope::biased_cycle_norm;
use super::{
    HarmonicElectrode, HarmonicElectrodeLoad, HarmonicPiezoelectricControls,
    HarmonicPiezoelectricResponse, ReactiveElectrodeLoad,
};
use fs_exec::Cx;
use fs_math::c64::C64;

/// Solved, immutable DC equilibrium tied to the exact original model.
/// Preparation retains relative potentials and boundary topology once, then
/// frequency sweeps reuse the original stiffness, mass, coupling and geometry.
pub struct HarmonicOperatingPoint<'model, 'mesh> {
    model: &'model TetPiezoelectricModel<'mesh>,
    state: TetPiezoelectricSolution,
    relative_voltage_v: Vec<f64>,
    support_dofs: Vec<usize>,
    electrode_nodes: Vec<Vec<usize>>,
}

/// One branch's steady state, in the same order as `harmonic.loads`.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicBiasLoadResponse {
    /// Gauge-independent DC voltage from the first terminal to the second [V].
    pub voltage_v: f64,
    /// DC current from the first terminal to the second [A].
    pub current_a: f64,
    /// Additional external capacitor charge on its first terminal [C].
    /// Zero for resistive and inductive branches; ceramic charge is separate.
    pub charge_c: f64,
    /// Steady magnetic or external dielectric storage [J].
    pub stored_energy_j: f64,
    /// Original conductor's DC Joule power [W].
    pub joule_power_w: f64,
}

/// Incremental response plus the retained DC state and complete cycle admission.
/// The equilibrium itself remains available from `HarmonicOperatingPoint::state`.
#[derive(Clone, Debug)]
pub struct BiasedHarmonicResponse {
    /// Zero-mean phasors, incremental cycle storage and incremental source power.
    pub harmonic: HarmonicPiezoelectricResponse,
    /// Steady body free energy [J], with no sinusoidal factor of one half.
    pub bias_stored_energy_j: f64,
    /// DC body free energy plus mean incremental body storage [J]. Does not
    /// include kinetic energy, external branch storage or caloric energy.
    pub average_total_stored_energy_j: f64,
    /// Upper bound on the largest total-strain norm anywhere in the cycle.
    pub maximum_total_strain_norm: f64,
    /// Upper bound on the largest thermally shifted elastic-strain norm.
    pub maximum_elastic_strain_norm: f64,
    /// Upper bound on the largest DC-plus-AC electric-field norm [V/m].
    pub maximum_electric_field_v_m: f64,
    /// Original branch order; these are steady offsets, not AC phasors.
    pub dc_loads: Vec<HarmonicBiasLoadResponse>,
    /// Required ideal DC source current INTO each original electrode [A].
    /// These currents hold the supplied equilibrium against steady shunt flow.
    pub dc_electrode_source_current_a: Vec<f64>,
    /// Required external DC electrical supply, evaluated from branch voltage
    /// differences to preserve gauge invariance [W].
    pub dc_electrical_supply_w: f64,
    /// DC conductor heat [W], distinct from incremental mean load heat.
    pub dc_load_heat_w: f64,
    /// DC plus incremental mean conductor heat [W]. Maintaining temperature
    /// requires this heat to leave; this fixed-temperature API does not evolve it.
    pub average_total_load_heat_w: f64,
    /// DC external capacitor/inductor storage, separate from the ceramic [J].
    pub dc_external_stored_energy_j: f64,
    /// DC plus cycle-average external reactive storage [J].
    pub average_total_external_stored_energy_j: f64,
}

impl<'mesh> TetPiezoelectricModel<'mesh> {
    /// Solve the actual DC mechanical/electrical boundary problem, including
    /// sourced thermal strain and clamped pyroelectric polarization. This
    /// immutable equilibrium can then supply many harmonic frequencies.
    ///
    /// Static electrode charge is the ceramic's free charge, not external
    /// capacitor charge. Attached shunts are accounted when solving a frequency:
    /// required ideal DC source currents/heat are reported explicitly. The
    /// caller supplies the DC bias; this is not an autonomous bias-circuit solve.
    #[allow(clippy::too_many_arguments)]
    pub fn harmonic_operating_point(
        &self,
        prescribed_displacements_m: &[(usize, f64)],
        electrodes: &[PiezoelectricElectrode<'_>],
        nodal_forces_n: &[[f64; 3]],
        nodal_free_charge_c: &[f64],
        controls: TetPiezoelectricControls,
        cx: &Cx<'_>,
    ) -> Result<HarmonicOperatingPoint<'_, 'mesh>> {
        let (state, full) = self.solve_with_state(
            prescribed_displacements_m,
            electrodes,
            nodal_forces_n,
            nodal_free_charge_c,
            controls,
            cx,
        )?;
        HarmonicOperatingPoint::from_static(
            self,
            state,
            full,
            prescribed_displacements_m,
            electrodes,
            cx,
        )
    }
}

impl<'model, 'mesh> HarmonicOperatingPoint<'model, 'mesh> {
    fn from_static(
        model: &'model TetPiezoelectricModel<'mesh>,
        state: TetPiezoelectricSolution,
        full: Vec<f64>,
        supports: &[(usize, f64)],
        electrodes: &[PiezoelectricElectrode<'_>],
        cx: &Cx<'_>,
    ) -> Result<Self> {
        poll(cx)?;
        let mut support_dofs = supports.iter().map(|&(dof, _)| dof).collect::<Vec<_>>();
        support_dofs.sort_unstable();
        let mut electrode_nodes = Vec::with_capacity(electrodes.len());
        for electrode in electrodes {
            poll(cx)?;
            let mut nodes = electrode.nodes.to_vec();
            nodes.sort_unstable();
            electrode_nodes.push(nodes);
        }
        poll(cx)?;
        Ok(Self {
            model,
            state,
            relative_voltage_v: full[3 * model.nodes_m.len()..].to_vec(),
            support_dofs,
            electrode_nodes,
        })
    }

    /// Exact original mesh/operator/material owner; no tangent copy is created.
    pub const fn model(&self) -> &'model TetPiezoelectricModel<'mesh> {
        self.model
    }

    /// Actual steady fields, forces, electrode charges and free-energy audit.
    pub const fn state(&self) -> &TetPiezoelectricSolution {
        &self.state
    }

    /// Solve zero-mean perturbations about this prepared equilibrium. Mechanical
    /// support DOFs and electrode node sets/order must match the DC problem;
    /// electrode source kinds may differ between DC and AC. Port maps must
    /// belong to the same original model. Their waves have zero DC pressure.
    ///
    /// Resistive/RL DC flow and capacitor DC charge follow the actual bias
    /// voltage and sourced branches. Ideal lossless inductors need an explicit
    /// persistent current through `solve_harmonic_with_inductor_bias`.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_harmonic(
        &self,
        omega_rad_s: f64,
        acoustic_ports: &[HarmonicAcousticPort<'_>],
        prescribed_displacements: &[(usize, C64)],
        electrodes: &[HarmonicElectrode<'_>],
        loads: &[HarmonicElectrodeLoad],
        nodal_forces_n: &[[C64; 3]],
        nodal_free_charge_c: &[C64],
        controls: HarmonicPiezoelectricControls,
        cx: &Cx<'_>,
    ) -> Result<BiasedHarmonicResponse> {
        self.solve_harmonic_with_inductor_bias(
            omega_rad_s,
            acoustic_ports,
            &[],
            prescribed_displacements,
            electrodes,
            loads,
            nodal_forces_n,
            nodal_free_charge_c,
            controls,
            cx,
        )
    }

    /// As `solve_harmonic`, with `(load_index, persistent_current_a)` for every
    /// ideal lossless inductor. Explicit zero is allowed. Such a branch must
    /// have exactly zero DC voltage; otherwise its flux cannot be stationary.
    /// Extra, duplicate, nonfinite or non-inductor entries refuse. Finite-R RL
    /// currents are determined by their DC voltage and original resistance.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_harmonic_with_inductor_bias(
        &self,
        omega_rad_s: f64,
        acoustic_ports: &[HarmonicAcousticPort<'_>],
        inductor_bias_current_a: &[(usize, f64)],
        prescribed_displacements: &[(usize, C64)],
        electrodes: &[HarmonicElectrode<'_>],
        loads: &[HarmonicElectrodeLoad],
        nodal_forces_n: &[[C64; 3]],
        nodal_free_charge_c: &[C64],
        controls: HarmonicPiezoelectricControls,
        cx: &Cx<'_>,
    ) -> Result<BiasedHarmonicResponse> {
        poll(cx)?;
        self.admit_topology(prescribed_displacements, electrodes, cx)?;
        let harmonic = self.model.solve_harmonic_internal(
            Some(self),
            omega_rad_s,
            None,
            acoustic_ports,
            prescribed_displacements,
            electrodes,
            loads,
            nodal_forces_n,
            nodal_free_charge_c,
            controls,
            cx,
        )?;
        let mut maximum_total_strain_norm = 0.0_f64;
        let mut maximum_elastic_strain_norm = 0.0_f64;
        let mut maximum_electric_field_v_m = 0.0_f64;
        for (cell, (bias, ac)) in self.state.cells.iter().zip(&harmonic.cells).enumerate() {
            poll(cx)?;
            let strain = biased_cycle_norm(bias.strain_mandel, ac.strain_mandel);
            let elastic = biased_cycle_norm(
                self.model.elastic_strain(cell, bias.strain_mandel),
                ac.strain_mandel,
            );
            let electric =
                biased_cycle_norm(bias.response.electric_field_v_m, ac.electric_field_v_m);
            if !strain.is_finite()
                || !elastic.is_finite()
                || !electric.is_finite()
                || strain > controls.field.maximum_strain_norm
                || elastic > controls.field.maximum_strain_norm
                || electric > controls.field.maximum_electric_field_v_m
            {
                return Err(invalid(
                    "DC plus harmonic cycle leaves the sourced field envelope",
                ));
            }
            maximum_total_strain_norm = maximum_total_strain_norm.max(strain);
            maximum_elastic_strain_norm = maximum_elastic_strain_norm.max(elastic);
            maximum_electric_field_v_m = maximum_electric_field_v_m.max(electric);
        }
        // The native solve has already admitted every load, its endpoints,
        // fixed source conditions and common voltage-reference component.
        let mut explicit_current = vec![None; loads.len()];
        for &(id, current) in inductor_bias_current_a {
            poll(cx)?;
            if !current.is_finite()
                || id >= loads.len()
                || explicit_current[id].is_some()
                || !matches!(
                    &loads[id],
                    HarmonicElectrodeLoad::Reactive(ReactiveElectrodeLoad::SeriesRl {
                        conductor: None,
                        ..
                    })
                )
            {
                return Err(invalid(
                    "one explicit finite DC current per lossless inductor required",
                ));
            }
            explicit_current[id] = Some(current);
        }
        let mut dc_loads = Vec::with_capacity(loads.len());
        let mut dc_electrode_source_current_a = vec![0.0; electrodes.len()];
        let mut dc_electrical_supply_w = 0.0;
        let mut dc_load_heat_w = 0.0;
        let mut dc_external_stored_energy_j = 0.0;
        let mut ac_external_stored_energy_j = 0.0;
        for (id, (load, ac)) in loads.iter().zip(&harmonic.loads).enumerate() {
            poll(cx)?;
            let [a, b] = load.electrodes();
            let voltage_v = self.relative_voltage_v[electrodes[a].nodes[0]]
                - self.relative_voltage_v[electrodes[b].nodes[0]];
            let mut response = HarmonicBiasLoadResponse {
                voltage_v,
                current_a: 0.0,
                charge_c: 0.0,
                stored_energy_j: 0.0,
                joule_power_w: 0.0,
            };
            match load {
                HarmonicElectrodeLoad::Resistive(source) => {
                    response.current_a = voltage_v / source.conductor.resistance_ohm();
                    response.joule_power_w = source
                        .conductor
                        .joule_power_w(response.current_a)
                        .map_err(|_| invalid("DC resistor power is not representable"))?;
                }
                HarmonicElectrodeLoad::Reactive(ReactiveElectrodeLoad::SeriesRl {
                    inductor,
                    conductor,
                    ..
                }) => {
                    response.current_a = if let Some(conductor) = conductor {
                        let current = voltage_v / conductor.resistance_ohm();
                        response.joule_power_w = conductor
                            .joule_power_w(current)
                            .map_err(|_| invalid("DC winding power is not representable"))?;
                        current
                    } else {
                        if voltage_v != 0.0 {
                            return Err(invalid(
                                "lossless inductor cannot sustain nonzero DC voltage",
                            ));
                        }
                        explicit_current[id].ok_or_else(|| {
                            invalid("lossless inductor requires an explicit persistent DC current")
                        })?
                    };
                    response.stored_energy_j = inductor
                        .stored_energy_j(inductor.inductance_h() * response.current_a)
                        .map_err(|_| invalid("DC winding storage is not representable"))?;
                    ac_external_stored_energy_j += 0.5
                        * inductor
                            .stored_energy_j(inductor.inductance_h() * ac.current_a.abs())
                            .map_err(|_| invalid("AC winding storage is not representable"))?;
                }
                HarmonicElectrodeLoad::Reactive(ReactiveElectrodeLoad::Capacitor {
                    capacitor,
                    ..
                }) => {
                    response.charge_c = capacitor.capacitance_f() * voltage_v;
                    response.stored_energy_j =
                        capacitor.stored_energy_j(response.charge_c).map_err(|_| {
                            invalid("DC external capacitor storage is not representable")
                        })?;
                    ac_external_stored_energy_j += 0.5
                        * capacitor
                            .stored_energy_j(capacitor.capacitance_f() * ac.voltage_v.abs())
                            .map_err(|_| {
                                invalid("AC external capacitor storage is not representable")
                            })?;
                }
            }
            if [
                response.voltage_v,
                response.current_a,
                response.charge_c,
                response.stored_energy_j,
                response.joule_power_w,
            ]
            .iter()
            .any(|value| !value.is_finite())
            {
                return Err(invalid("DC load state is not representable"));
            }
            dc_electrode_source_current_a[a] += response.current_a;
            dc_electrode_source_current_a[b] -= response.current_a;
            dc_electrical_supply_w += response.voltage_v * response.current_a;
            dc_load_heat_w += response.joule_power_w;
            dc_external_stored_energy_j += response.stored_energy_j;
            dc_loads.push(response);
        }
        let average_total_load_heat_w = dc_load_heat_w + harmonic.average_load_heat_w;
        let average_total_stored_energy_j =
            self.state.stored_energy_j + harmonic.average_stored_energy_j;
        let average_total_external_stored_energy_j =
            dc_external_stored_energy_j + ac_external_stored_energy_j;
        if [
            dc_electrical_supply_w,
            dc_load_heat_w,
            dc_external_stored_energy_j,
            average_total_load_heat_w,
            average_total_stored_energy_j,
            average_total_external_stored_energy_j,
        ]
        .iter()
        .chain(&dc_electrode_source_current_a)
        .any(|value| !value.is_finite())
        {
            return Err(invalid("DC plus AC work or storage is not representable"));
        }
        if (dc_electrical_supply_w - dc_load_heat_w).abs() > controls.power_tolerance_w {
            return Err(invalid(
                "steady electrical supply does not match conductor heat",
            ));
        }
        poll(cx)?;
        Ok(BiasedHarmonicResponse {
            harmonic,
            bias_stored_energy_j: self.state.stored_energy_j,
            average_total_stored_energy_j,
            maximum_total_strain_norm,
            maximum_elastic_strain_norm,
            maximum_electric_field_v_m,
            dc_loads,
            dc_electrode_source_current_a,
            dc_electrical_supply_w,
            dc_load_heat_w,
            average_total_load_heat_w,
            dc_external_stored_energy_j,
            average_total_external_stored_energy_j,
        })
    }

    fn admit_topology(
        &self,
        supports: &[(usize, C64)],
        electrodes: &[HarmonicElectrode<'_>],
        cx: &Cx<'_>,
    ) -> Result<()> {
        if supports.len() != self.support_dofs.len()
            || electrodes.len() != self.electrode_nodes.len()
        {
            return Err(invalid(
                "DC and AC must retain the same supports and electrode topology",
            ));
        }
        for &(dof, _) in supports {
            poll(cx)?;
            if self.support_dofs.binary_search(&dof).is_err() {
                return Err(invalid(
                    "harmonic support is absent from the DC operating point",
                ));
            }
        }
        for (electrode, expected) in electrodes.iter().zip(&self.electrode_nodes) {
            poll(cx)?;
            if electrode.nodes.len() != expected.len() {
                return Err(invalid(
                    "harmonic electrode differs from the DC operating point",
                ));
            }
            for node in electrode.nodes {
                poll(cx)?;
                if expected.binary_search(node).is_err() {
                    return Err(invalid(
                        "harmonic electrode differs from the DC operating point",
                    ));
                }
            }
        }
        // Duplicate indices/nodes and every remaining boundary invariant are
        // checked by the native harmonic boundary owner before solving.
        Ok(())
    }
}
