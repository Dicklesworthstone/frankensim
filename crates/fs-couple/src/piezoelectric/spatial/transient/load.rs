//! Source-bound resistors on the original spatial electromechanical midpoint solve.
//!
//! A load changes the floating electrode's charge equation, not the solid's
//! damping. Its positive Joule heat belongs to the retained conductor; no heat
//! is assigned to the ceramic. Any attached intrinsic material law retains its
//! own separate ceramic heat account; no temperature evolution is inferred.

use super::*;
use fs_material::conductor::{OhmicDrive, ResolvedConductor};

/// Retained material capacitance and inductance on the same spatial solve.
pub mod reactive;
/// Finite load heating and material feedback through the original enthalpy solve.
pub mod thermal;

/// External electrode source with the existing voltage/charge topology fixed.
#[derive(Clone, Copy, Debug)]
pub enum ElectrodeDrive {
    /// Prescribed endpoint potential [V], linearly interpolated over the step.
    Voltage(f64),
    /// Held current supplied INTO an existing charge electrode [A]. Zero is open.
    Current(f64),
}

/// Uniform Ohmic conductor joining two original electrode indices.
#[derive(Clone, Debug)]
pub struct ResistiveElectrodeLoad {
    /// Positive current flows from the first electrode to the second.
    pub electrodes: [usize; 2],
    /// Original immutable source, query conditions and actual resistor geometry.
    pub conductor: Arc<ResolvedConductor>,
}

/// Dissipation in one named physical conductor, never numerical/ceramic damping.
#[derive(Clone, Debug)]
pub struct ElectrodeLoadHeat {
    /// Original endpoints and source shared without cloning the material bundle.
    pub load: ResistiveElectrodeLoad,
    /// Signed midpoint current from the first electrode to the second [A].
    pub midpoint_current_a: f64,
    /// Integrated signed charge transferred through this conductor [C].
    pub transferred_charge_c: f64,
    /// Nonnegative native Ohmic Joule power integrated over the window [J].
    pub joule_heat_j: f64,
}

/// Independent account enclosing the body, external sources and load losses.
#[derive(Clone, Debug)]
pub struct LoadedPiezoelectricStep {
    /// Work of prescribed currents/voltages and supplied nodal charge [J].
    pub electrical_source_work_j: f64,
    /// Sum of actual conductor losses [J].
    pub joule_heat_j: f64,
    /// Delta(all body storage) + conductor and intrinsic material heat + signed
    /// acoustic outflow - mechanical/source work [J]. Native reports keep each
    /// physical heat receiver and outgoing wave energy distinct.
    pub energy_residual_j: f64,
    /// One source-bound record per input resistor, in original order.
    pub loads: Vec<ElectrodeLoadHeat>,
    source_work: Dd,
}

impl LoadedPiezoelectricStep {
    /// Independent external current/voltage work before diagnostic rounding [J].
    pub const fn retained_electrical_source_work_j(&self) -> Dd {
        self.source_work
    }
}

impl<'model, 'mesh> TetPiezoelectricTransient<'model, 'mesh> {
    /// Solve the SAME coupled midpoint system with source-bound shunts. A held
    /// external current controls charge supply, while resistor current is solved
    /// simultaneously from the actual midpoint potential. Source types remain
    /// fixed within this step; explicit electrode release occurs between steps.
    /// Resistors may join only the same retained dielectric/gauge component.
    /// The candidate's ordinary step report accounts for net electrode work into
    /// the body. The returned load report instead encloses body AND resistors,
    /// separating external electrical source work from their dissipated heat.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_loads_candidate(
        &self,
        dt_s: f64,
        drives: &[ElectrodeDrive],
        loads: &[ResistiveElectrodeLoad],
        midpoint_nodal_forces_n: &[[f64; 3]],
        endpoint_nodal_free_charge_c: &[f64],
        controls: TetPiezoelectricControls,
        cx: &Cx<'_>,
    ) -> Result<(
        TetPiezoelectricTransientCandidate<'model, 'mesh>,
        LoadedPiezoelectricStep,
    )> {
        self.advance_loaded_candidate(
            dt_s,
            None,
            drives,
            loads,
            midpoint_nodal_forces_n,
            endpoint_nodal_free_charge_c,
            controls,
            cx,
        )
    }

    /// Source-bound shunts with explicit motion of every original support DOF.
    /// Support velocities use the same kinematics and reaction-work convention
    /// as [`Self::step_with_support_motion_candidate`]. Resistor heat and direct
    /// sensing back-reaction are solved simultaneously with that imposed motion.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_loaded_support_motion_candidate(
        &self,
        dt_s: f64,
        endpoint_support_velocities_m_s: &[(usize, f64)],
        drives: &[ElectrodeDrive],
        loads: &[ResistiveElectrodeLoad],
        midpoint_nodal_forces_n: &[[f64; 3]],
        endpoint_nodal_free_charge_c: &[f64],
        controls: TetPiezoelectricControls,
        cx: &Cx<'_>,
    ) -> Result<(
        TetPiezoelectricTransientCandidate<'model, 'mesh>,
        LoadedPiezoelectricStep,
    )> {
        self.advance_loaded_candidate(
            dt_s,
            Some(endpoint_support_velocities_m_s),
            drives,
            loads,
            midpoint_nodal_forces_n,
            endpoint_nodal_free_charge_c,
            controls,
            cx,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn advance_loaded_candidate(
        &self,
        dt_s: f64,
        endpoint_support_velocities_m_s: Option<&[(usize, f64)]>,
        drives: &[ElectrodeDrive],
        loads: &[ResistiveElectrodeLoad],
        midpoint_nodal_forces_n: &[[f64; 3]],
        endpoint_nodal_free_charge_c: &[f64],
        controls: TetPiezoelectricControls,
        cx: &Cx<'_>,
    ) -> Result<(
        TetPiezoelectricTransientCandidate<'model, 'mesh>,
        LoadedPiezoelectricStep,
    )> {
        self.advance_loaded_with_acoustic_candidate(
            dt_s,
            endpoint_support_velocities_m_s,
            &[],
            drives,
            loads,
            midpoint_nodal_forces_n,
            endpoint_nodal_free_charge_c,
            controls,
            cx,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn advance_loaded_with_acoustic_candidate(
        &self,
        dt_s: f64,
        endpoint_support_velocities_m_s: Option<&[(usize, f64)]>,
        acoustic_ports: &[acoustic::TransientAcousticPort<'_>],
        drives: &[ElectrodeDrive],
        loads: &[ResistiveElectrodeLoad],
        midpoint_nodal_forces_n: &[[f64; 3]],
        endpoint_nodal_free_charge_c: &[f64],
        controls: TetPiezoelectricControls,
        cx: &Cx<'_>,
    ) -> Result<(
        TetPiezoelectricTransientCandidate<'model, 'mesh>,
        LoadedPiezoelectricStep,
    )> {
        poll(cx)?;
        let endpoints = endpoint_conditions(self, dt_s, drives, cx)?;
        let (candidate, heat, _) = self.advance_with_acoustic_loads(
            dt_s,
            endpoint_support_velocities_m_s,
            &endpoints,
            loads,
            &[],
            &[],
            acoustic_ports,
            midpoint_nodal_forces_n,
            endpoint_nodal_free_charge_c,
            controls,
            cx,
        )?;
        let next = &candidate.body.state;
        let mut outgoing = vec![Dd::ZERO; drives.len()];
        let mut joule_heat = Dd::ZERO;
        for flow in &heat {
            poll(cx)?;
            outgoing[flow.load.electrodes[0]] =
                outgoing[flow.load.electrodes[0]] + Dd::from_f64(flow.transferred_charge_c);
            outgoing[flow.load.electrodes[1]] =
                outgoing[flow.load.electrodes[1]] - Dd::from_f64(flow.transferred_charge_c);
            joule_heat = joule_heat + Dd::from_f64(flow.joule_heat_j);
        }
        let source_work = source_work(
            self,
            next,
            dt_s,
            drives,
            endpoint_nodal_free_charge_c,
            &outgoing,
            cx,
        )?;
        let electrical_source_work_j = source_work.to_f64();
        let joule_heat_j = joule_heat.to_f64();
        let energy_residual_j = (candidate.step.retained_energy_change_j()
            + joule_heat
            + Dd::from_f64(candidate.step.material_heat_j)
            + Dd::from_f64(candidate.step.acoustic_outward_work_j)
            - candidate.step.retained_mechanical_work_j()
            - source_work)
            .to_f64();
        if [electrical_source_work_j, joule_heat_j, energy_residual_j]
            .iter()
            .any(|v| !v.is_finite())
        {
            return Err(invalid(
                "loaded piezoelectric energy account is not representable",
            ));
        }
        if energy_residual_j.abs() > controls.energy_tolerance_j {
            return Err(TetPiezoelectricError::Balance(
                "body energy, physical heat, acoustic outflow and external source work disagree",
            ));
        }
        poll(cx)?;
        Ok((
            candidate,
            LoadedPiezoelectricStep {
                electrical_source_work_j,
                joule_heat_j,
                energy_residual_j,
                loads: heat,
                source_work,
            },
        ))
    }
}

fn endpoint_conditions(
    body: &TetPiezoelectricTransient<'_, '_>,
    dt_s: f64,
    drives: &[ElectrodeDrive],
    cx: &Cx<'_>,
) -> Result<Vec<ElectrodeCondition>> {
    if !positive(dt_s) || drives.len() != body.state.electrode_conditions.len() {
        return Err(invalid(
            "positive timestep and complete external electrode drives required",
        ));
    }
    drives
        .iter()
        .zip(&body.state.electrode_conditions)
        .map(|(&drive, &old)| {
            poll(cx)?;
            match (drive, old) {
                (ElectrodeDrive::Voltage(v), ElectrodeCondition::Voltage(_)) if v.is_finite() => {
                    Ok(ElectrodeCondition::Voltage(v))
                }
                (ElectrodeDrive::Current(i), ElectrodeCondition::Charge(q)) if i.is_finite() => {
                    Ok(ElectrodeCondition::Charge(dt_s.mul_add(i, q)))
                }
                _ => Err(invalid(
                    "finite drives must preserve voltage/current electrode topology",
                )),
            }
        })
        .collect()
}

#[allow(clippy::too_many_arguments)]
fn source_work(
    body: &TetPiezoelectricTransient<'_, '_>,
    next: &TetPiezoelectricTransientState,
    dt_s: f64,
    drives: &[ElectrodeDrive],
    endpoint_nodal_free_charge_c: &[f64],
    outgoing: &[Dd],
    cx: &Cx<'_>,
) -> Result<Dd> {
    let midpoint = |node: usize| {
        Dd::from_f64(0.5 * body.state.relative_voltage_v[node])
            + Dd::from_f64(0.5 * next.relative_voltage_v[node])
    };
    let mut work = Products::default();
    for (node, (&a, &b)) in body
        .state
        .nodal_free_charge_c
        .iter()
        .zip(endpoint_nodal_free_charge_c)
        .enumerate()
    {
        poll(cx)?;
        work.add(midpoint(node), Dd::from_f64(b) - Dd::from_f64(a));
    }
    for (id, drive) in drives.iter().enumerate() {
        poll(cx)?;
        let supplied = match *drive {
            ElectrodeDrive::Current(i) => Dd::from_f64(dt_s) * Dd::from_f64(i),
            ElectrodeDrive::Voltage(_) => {
                Dd::from_f64(next.field.electrodes[id].charge_c)
                    - Dd::from_f64(body.state.field.electrodes[id].charge_c)
                    + outgoing[id]
            }
        };
        work.add(midpoint(body.boundary.electrode_nodes[id][0]), supplied);
    }
    let work = work.value();
    if !work.hi.is_finite() || !work.lo.is_finite() {
        return Err(invalid("electrode source work is not representable"));
    }
    Ok(work)
}

/// Stamp in the original COO assembly once; never rebuild an assembled matrix.
pub(super) fn stamp(
    body: &TetPiezoelectricTransient<'_, '_>,
    loads: &[ResistiveElectrodeLoad],
    dt: f64,
    boundary: &mut BoundaryMap,
    cx: &Cx<'_>,
) -> Result<Vec<(usize, usize, f64)>> {
    let mut entries = Vec::with_capacity(
        loads
            .len()
            .checked_mul(4)
            .ok_or_else(|| invalid("electrode load count overflows"))?,
    );
    let nu = 3 * body.model.nodes_m.len();
    for load in loads {
        poll(cx)?;
        let [a, b] = load.electrodes;
        if a == b
            || a >= body.boundary.electrode_nodes.len()
            || b >= body.boundary.electrode_nodes.len()
            || !load
                .conductor
                .source_property()
                .is_static_at_fixed_external_conditions()
        {
            return Err(invalid(
                "distinct in-range electrodes and a fixed scalar Ohmic source required",
            ));
        }
        let [a, b] = [
            body.boundary.electrode_nodes[a][0],
            body.boundary.electrode_nodes[b][0],
        ];
        if boundary.voltage_reference_electrode[a] != boundary.voltage_reference_electrode[b] {
            return Err(invalid(
                "resistor endpoints must share the same dielectric voltage-reference component",
            ));
        }
        let half_conductance_time = 0.5 * dt / load.conductor.resistance_ohm();
        if !positive(half_conductance_time) {
            return Err(invalid(
                "timestep times electrode conductance is not representable",
            ));
        }
        let charge = half_conductance_time
            * (body.state.relative_voltage_v[a] - body.state.relative_voltage_v[b]);
        for (node, sign) in [(a, 1.0), (b, -1.0)] {
            if let Some(row) = boundary.index[nu + node] {
                boundary.rhs[row] += sign * charge;
            }
        }
        for (row, col, sign) in [(a, a, -1.0), (a, b, 1.0), (b, a, 1.0), (b, b, -1.0)] {
            entries.push((nu + row, nu + col, sign * half_conductance_time));
        }
    }
    Ok(entries)
}

/// Actual endpoint charge includes solved resistor flow. Native recovery then
/// verifies both Gauss constraints with those same midpoint/endpoint charges.
#[allow(clippy::too_many_arguments)]
pub(super) fn recover_load_charge(
    body: &TetPiezoelectricTransient<'_, '_>,
    loads: &[ResistiveElectrodeLoad],
    dt: f64,
    midpoint_voltage: &[f64],
    endpoint: &mut [ElectrodeCondition],
    normalized_endpoint: &mut [ElectrodeCondition],
    midpoint: &mut [ElectrodeCondition],
    cx: &Cx<'_>,
) -> Result<Vec<ElectrodeLoadHeat>> {
    let mut heat = Vec::with_capacity(loads.len());
    for load in loads {
        poll(cx)?;
        let [a, b] = load.electrodes;
        let voltage = midpoint_voltage[body.boundary.electrode_nodes[a][0]]
            - midpoint_voltage[body.boundary.electrode_nodes[b][0]];
        let current = load
            .conductor
            .current_a_for_drive(OhmicDrive::Voltage(voltage))
            .map_err(|_| invalid("loaded electrode current is not representable"))?;
        let charge = dt * current;
        let joule_heat_j = dt
            * load
                .conductor
                .joule_power_w(current)
                .map_err(|_| invalid("loaded electrode Joule power is not representable"))?;
        if !charge.is_finite()
            || !joule_heat_j.is_finite()
            || current != 0.0 && (charge == 0.0 || joule_heat_j == 0.0)
        {
            return Err(invalid(
                "electrode load charge or heat is not representable",
            ));
        }
        for (id, flow) in [(a, charge), (b, -charge)] {
            if let ElectrodeCondition::Charge(value) = &mut endpoint[id] {
                *value -= flow;
            }
        }
        heat.push(ElectrodeLoadHeat {
            load: load.clone(),
            midpoint_current_a: current,
            transferred_charge_c: charge,
            joule_heat_j,
        });
    }
    update_charge_conditions(body, endpoint, normalized_endpoint, midpoint, cx)?;
    Ok(heat)
}

fn update_charge_conditions(
    body: &TetPiezoelectricTransient<'_, '_>,
    endpoint: &[ElectrodeCondition],
    normalized_endpoint: &mut [ElectrodeCondition],
    midpoint: &mut [ElectrodeCondition],
    cx: &Cx<'_>,
) -> Result<()> {
    for (id, condition) in endpoint.iter().enumerate() {
        poll(cx)?;
        if let ElectrodeCondition::Charge(next) = *condition {
            if !next.is_finite() {
                return Err(invalid("loaded electrode charge overflows"));
            }
            let ElectrodeCondition::Charge(old) = body.state.electrode_conditions[id] else {
                unreachable!()
            };
            normalized_endpoint[id] = ElectrodeCondition::Charge(next);
            midpoint[id] = ElectrodeCondition::Charge(0.5 * old + 0.5 * next);
        }
    }
    Ok(())
}
