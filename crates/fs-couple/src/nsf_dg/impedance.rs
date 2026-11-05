//! Causal passive acoustic loads on original resolved NSF faces or apertures.
//!
//! The existing fs-vfit series R/L/C and RL/RC relaxation owner supplies the
//! frequency-dependent reaction and its physical storage. An original NSF trial
//! supplies mean acoustic pressure and outward volume flow over the SAME
//! interval. Their characteristic waves are matched in one bounded solve.
//! Rejected roots never advance either owner; final candidates publish together.
//!
//! This is a zero-reference-flow, uniform acoustic-port projection. It is not a
//! pointwise reactive wall, an exterior mass inventory, or full thermodynamic
//! transport into the load. Original fluid mass/enthalpy and acoustic work stay
//! separate. Signed pressure/flow projection covariance is reported and bounded,
//! never relabeled as material heat. The original boundary and interval require
//! refinement; passivity of the load alone does not qualify interface accuracy.

use super::{
    NsfBoundaryDrive, NsfBoundaryExchange, NsfCharacteristicPatch, NsfDgError, NsfDgLimits,
    NsfGrid3, coupled_exchange, poll, preview_interval, upstream,
};
use fs_exec::CancelGate;
use fs_lbm::lattice3::Lattice3;
use fs_lbm::nsf3::{EnergyLedger, NsfReservoirExchange, NsfReservoirFlow};
use fs_material::gas::GasState;
use fs_vfit::relaxation::{RelaxationFrame, RelaxationImpedance, RelaxationImpedanceSpec};
use fs_vfit::waveguide::WaveguideError;

/// Original gas and clock bound to an explicitly declared circular exterior.
pub mod radiation;

/// One accepted original fluid interval and its complete causal load response.
#[derive(Debug, Clone)]
pub struct NsfImpedanceFrame {
    pub step: u64,
    /// Includes the final original-fluid replay; all evaluations share one budget.
    pub evaluations: usize,
    /// Accepted incoming characteristic on the coupled face or aperture [m/s].
    pub incoming_characteristic_m_s: f64,
    /// Held control on the opposite source face or its selected aperture.
    pub drive: NsfBoundaryDrive,
    /// Ideal series effort: p_zone - p_load [Pa].
    pub series_pressure_pa: f64,
    pub faces: [NsfBoundaryExchange; 2],
    pub coupled_patch: Option<NsfBoundaryExchange>,
    pub source_patch: Option<NsfBoundaryExchange>,
    /// Actual original exchange of the selected loaded face or aperture.
    pub load_exchange: NsfBoundaryExchange,
    /// Complete original RLC and relaxation states, storage, loss and port work.
    pub load: RelaxationFrame,
    pub flow: EnergyLedger,
    pub flow_mass_change_kg: f64,
    pub flow_energy_change_j: f64,
    pub reservoir_mass_out_kg: f64,
    pub reservoir_flows: Vec<NsfReservoirFlow>,
    pub reservoir_exchanges: Vec<NsfReservoirExchange>,
    /// p_zone - p_load - series_pressure [Pa].
    pub pressure_residual_pa: f64,
    /// Original outward volume flow minus load flow [m^3/s].
    pub flow_residual_m3_s: f64,
    /// Integral p' u_n minus mean(p') mean(Q) dt [J].
    pub projection_covariance_j: f64,
    /// Work delivered by the ideal series source: -series_pressure Q_out dt [J].
    pub source_work_j: f64,
    /// Resolved acoustic work + source work - complete load work [J].
    pub acoustic_work_residual_j: f64,
    /// Actual load storage change + nonnegative loss - supplied work [J].
    pub load_energy_residual_j: f64,
    /// Native streamed mass minus rho_ref times acoustic displacement [kg].
    pub mass_projection_difference_kg: f64,
}

/// An original NSF owner coupled to one causal passive acoustic impedance.
///
/// The load's volume-flow units are Pa s/m^3, not specific impedance Pa s/m.
/// Its characteristic reference rho*c/area comes from the selected original
/// face material and quadrature. Reactive states start at zero only on cold
/// construction and are retained across every accepted step.
pub struct NsfImpedance<L: Lattice3<27>> {
    zone: NsfGrid3<L>,
    load: RelaxationImpedance,
    limits: NsfDgLimits,
    high: bool,
    reference: GasState,
    area_m2: f64,
    port_impedance_pa_s_m3: f64,
    time_step_s: f64,
    patch_background: [Option<NsfBoundaryDrive>; 2],
    steps: u64,
}

impl<L: Lattice3<27>> NsfImpedance<L> {
    /// Couple one complete original face; the other face remains a typed source.
    pub fn new(
        zone: NsfGrid3<L>,
        load: RelaxationImpedanceSpec,
        high: bool,
        limits: NsfDgLimits,
        gate: &CancelGate,
    ) -> Result<Self, NsfDgError> {
        Self::bind(zone, load, high, limits, gate, [false; 2])
    }

    /// Load one declared aperture, retaining its complementary face drive.
    pub fn new_on_patch(
        zone: NsfGrid3<L>,
        load: RelaxationImpedanceSpec,
        high: bool,
        limits: NsfDgLimits,
        gate: &CancelGate,
    ) -> Result<Self, NsfDgError> {
        let mut selected = [false; 2];
        selected[usize::from(high)] = true;
        Self::bind(zone, load, high, limits, gate, selected)
    }

    /// Load a complete face and drive the explicit aperture on the opposite face.
    pub fn new_with_source_patch(
        zone: NsfGrid3<L>,
        load: RelaxationImpedanceSpec,
        high: bool,
        limits: NsfDgLimits,
        gate: &CancelGate,
    ) -> Result<Self, NsfDgError> {
        let mut selected = [false; 2];
        selected[usize::from(!high)] = true;
        Self::bind(zone, load, high, limits, gate, selected)
    }

    /// Select independent loaded and driven apertures on the two original faces.
    pub fn new_on_patches(
        zone: NsfGrid3<L>,
        load: RelaxationImpedanceSpec,
        high: bool,
        limits: NsfDgLimits,
        gate: &CancelGate,
    ) -> Result<Self, NsfDgError> {
        Self::bind(zone, load, high, limits, gate, [true; 2])
    }

    fn bind(
        zone: NsfGrid3<L>,
        spec: RelaxationImpedanceSpec,
        high: bool,
        limits: NsfDgLimits,
        gate: &CancelGate,
        selected: [bool; 2],
    ) -> Result<Self, NsfDgError> {
        poll(gate)?;
        limits.validate_for_cells(zone.cell_count())?;
        let open = zone.characteristic_axis().ok_or(NsfDgError::Input(
            "causal impedance requires an original characteristic axis",
        ))?;
        if zone.steps() != 0
            || open
                .faces
                .iter()
                .any(|face| face.reference_velocity_m_s != [0.0; 3])
        {
            return Err(NsfDgError::Input(
                "causal impedance requires a cold original NSF grid and zero reference flow",
            ));
        }
        let drives = zone.boundary_drives().ok_or(NsfDgError::Input(
            "causal impedance requires the original surrounding drives",
        ))?;
        let mut patch_background = [None; 2];
        for side in 0..2 {
            if selected[side] != zone.characteristic_patch_on(side != 0).is_some() {
                return Err(NsfDgError::Input(
                    "select every original loaded or source aperture explicitly",
                ));
            }
            if selected[side] {
                patch_background[side] = Some(drives[side]);
            }
        }
        let traces = if selected[usize::from(high)] {
            zone.open_patch_traces_on(high)
        } else {
            zone.open_boundary_traces(high)
        }
        .ok_or(NsfDgError::Input(
            "causal impedance requires original face quadrature",
        ))?;
        let area = traces
            .first()
            .ok_or(NsfDgError::Input(
                "causal impedance requires a nonempty original face",
            ))?
            .area_m2
            * traces.len() as f64;
        let reference = zone
            .characteristic_reference_gas(high)
            .map_err(NsfDgError::Flow)?;
        let port = reference.density * reference.sound_speed / area;
        let dt = zone.law().dt_s * limits.substeps as f64;
        if !area.is_finite() || area <= 0.0 || !port.is_finite() || port <= 0.0 {
            return Err(NsfDgError::Input(
                "causal impedance needs finite original area and material wave impedance",
            ));
        }
        let load = RelaxationImpedance::new(spec, port, dt).map_err(load_error)?;
        Ok(Self {
            zone,
            load,
            limits,
            high,
            reference,
            area_m2: area,
            port_impedance_pa_s_m3: port,
            time_step_s: dt,
            patch_background,
            steps: 0,
        })
    }

    #[must_use]
    pub const fn zone(&self) -> &NsfGrid3<L> {
        &self.zone
    }
    #[must_use]
    pub const fn load(&self) -> &RelaxationImpedance {
        &self.load
    }
    #[must_use]
    pub const fn reference_gas(&self) -> &GasState {
        &self.reference
    }
    #[must_use]
    pub const fn area_m2(&self) -> f64 {
        self.area_m2
    }
    #[must_use]
    pub const fn port_impedance_pa_s_m3(&self) -> f64 {
        self.port_impedance_pa_s_m3
    }
    #[must_use]
    pub const fn time_step_s(&self) -> f64 {
        self.time_step_s
    }
    #[must_use]
    pub const fn coupled_face_is_high(&self) -> bool {
        self.high
    }
    #[must_use]
    pub const fn limits(&self) -> &NsfDgLimits {
        &self.limits
    }
    #[must_use]
    pub const fn accepted_steps(&self) -> u64 {
        self.steps
    }
    #[must_use]
    pub fn coupled_patch(&self) -> Option<NsfCharacteristicPatch> {
        self.zone.characteristic_patch_on(self.high)
    }
    #[must_use]
    pub fn source_patch(&self) -> Option<NsfCharacteristicPatch> {
        self.zone.characteristic_patch_on(!self.high)
    }
    #[must_use]
    pub const fn patch_background_drive_on(&self, high: bool) -> Option<NsfBoundaryDrive> {
        self.patch_background[high as usize]
    }

    /// Transfer both original owners, including all accepted reactive memory.
    #[must_use]
    pub fn into_parts(self) -> (NsfGrid3<L>, RelaxationImpedance) {
        (self.zone, self.load)
    }

    /// Solve one complete original interval with a live causal reaction.
    ///
    /// The source drive acts on the other face/aperture. The optional ideal
    /// series effort obeys p_zone = p_load + series_pressure_pa at the loaded
    /// face, whose outward volume flow is INTO the load on either grid side.
    /// Its independent work is -series_pressure_pa Q_out dt.
    ///
    /// The existing root budget includes one final original-fluid replay.
    /// Refusal, cancellation or projection excess publishes neither owner.
    pub fn step(
        &mut self,
        external: NsfBoundaryDrive,
        series_pressure_pa: f64,
        gate: &CancelGate,
    ) -> Result<NsfImpedanceFrame, NsfDgError> {
        poll(gate)?;
        if !external.is_finite()
            || !series_pressure_pa.is_finite()
            || self.steps >= self.limits.maximum_steps
        {
            return Err(NsfDgError::Input(
                "causal impedance input is nonfinite or its interval horizon is exhausted",
            ));
        }
        let limits = self.limits;
        let side = usize::from(self.high);
        let patches = self.patch_background;
        let port = self.port_impedance_pa_s_m3;
        let dt = self.time_step_s;
        let reference = self.reference;
        let load = self.load;
        let accepted = if patches[side].is_some() {
            self.zone.characteristic_patch_drive_on(self.high)
        } else {
            self.zone.boundary_drives().map(|drives| drives[side])
        };
        let seed = match accepted {
            Some(NsfBoundaryDrive::IncomingCharacteristic(value)) => Some(value),
            _ => None,
        };
        let mut evaluations = 0;
        let mut last_residual = 0.0;
        let zone = &mut self.zone;
        let (incoming, ()) = upstream::solve_scalar(limits, seed, |incoming| {
            poll(gate)?;
            // Reserve a complete replay so publication never reuses a stale
            // candidate or exceeds the caller's total original-cell work.
            if evaluations + 1 >= limits.maximum_evaluations {
                return Err(NsfDgError::Convergence {
                    residual_pa: last_residual,
                    evaluations,
                });
            }
            evaluations += 1;
            let mut drives = [external; 2];
            drives[side] = NsfBoundaryDrive::IncomingCharacteristic(incoming);
            let trial = preview_interval(zone, drives, patches, &[], limits, gate)?;
            let face = coupled_exchange(&trial, side)?;
            let (_, acoustic) = candidate_load(
                load,
                port,
                face.mean_acoustic_pressure_pa,
                face.normal_volume_flow_m3_s,
                series_pressure_pa,
            )
            .map_err(load_error)?;
            let (residual, scale) = matching_residual(face, acoustic, port, series_pressure_pa)?;
            last_residual = residual;
            Ok((residual, scale, ()))
        })?;
        poll(gate)?;
        evaluations += 1;
        let mut drives = [external; 2];
        drives[side] = NsfBoundaryDrive::IncomingCharacteristic(incoming);
        let fluid = preview_interval(zone, drives, patches, &[], limits, gate)?;
        let face = coupled_exchange(&fluid, side)?;
        let (candidate, acoustic) = candidate_load(
            load,
            port,
            face.mean_acoustic_pressure_pa,
            face.normal_volume_flow_m3_s,
            series_pressure_pa,
        )
        .map_err(load_error)?;
        let (wave_residual, scale) = matching_residual(face, acoustic, port, series_pressure_pa)?;
        let q = face.normal_volume_flow_m3_s;
        let pressure_residual =
            face.mean_acoustic_pressure_pa - series_pressure_pa - acoustic.port.pressure_pa;
        let covariance =
            face.acoustic_pressure_work_out_j - face.mean_acoustic_pressure_pa * (q * dt);
        let source_work = -series_pressure_pa * (q * dt);
        let work = face.acoustic_pressure_work_out_j + source_work - acoustic.port.supplied_work_j;
        let work_scale = face.acoustic_pressure_work_out_j.abs()
            + source_work.abs()
            + acoustic.port.supplied_work_j.abs();
        let allowance = limits.projection_absolute_tolerance_j
            + limits.projection_relative_tolerance * work_scale;
        let load_balance = acoustic.port.balance_residual_j();
        let load_scale = load.stored_energy_j()
            + acoustic.port.stored_energy_j
            + acoustic.port.dissipated_energy_j
            + acoustic.port.supplied_work_j.abs();
        if ![
            pressure_residual,
            covariance,
            source_work,
            work,
            work_scale,
            allowance,
            load_balance,
            load_scale,
        ]
        .iter()
        .all(|value| value.is_finite())
            || wave_residual.abs()
                > limits.pressure_absolute_tolerance_pa + limits.pressure_relative_tolerance * scale
            || pressure_residual.abs()
                > limits.pressure_absolute_tolerance_pa + limits.pressure_relative_tolerance * scale
            || acoustic.port.dissipated_energy_j < 0.0
            || load_balance.abs() > 1e-10 * load_scale.max(1e-300)
        {
            return Err(NsfDgError::Input(
                "original causal load pressure, flow or stored-energy work disagrees",
            ));
        }
        if work.abs() > allowance {
            return Err(NsfDgError::Projection {
                residual_j: work,
                allowance_j: allowance,
            });
        }
        let frame = NsfImpedanceFrame {
            step: self.steps + 1,
            evaluations,
            incoming_characteristic_m_s: incoming,
            drive: external,
            series_pressure_pa,
            faces: *fluid.boundary_exchange(),
            coupled_patch: fluid.patch_exchange_on(self.high).copied(),
            source_patch: fluid.patch_exchange_on(!self.high).copied(),
            load_exchange: face,
            load: acoustic,
            flow: *fluid.energy_ledger(),
            flow_mass_change_kg: fluid.grid().total_mass_kg() - fluid.previous_mass_kg(),
            flow_energy_change_j: fluid.grid().total_energy_j() - fluid.previous_energy_j(),
            reservoir_mass_out_kg: fluid.reservoir_mass_out_kg(),
            reservoir_flows: fluid.reservoir_flows().to_vec(),
            reservoir_exchanges: fluid.reservoir_exchanges().to_vec(),
            pressure_residual_pa: pressure_residual,
            flow_residual_m3_s: q - acoustic.port.flow_m3_s,
            projection_covariance_j: covariance,
            source_work_j: source_work,
            acoustic_work_residual_j: work,
            load_energy_residual_j: load_balance,
            mass_projection_difference_kg: face.mass_out_kg - reference.density * (q * dt),
        };
        if ![
            frame.flow_mass_change_kg,
            frame.flow_energy_change_j,
            frame.mass_projection_difference_kg,
            frame.reservoir_mass_out_kg,
        ]
        .iter()
        .all(|value| value.is_finite())
        {
            return Err(NsfDgError::Input(
                "causal impedance interval account overflowed",
            ));
        }
        poll(gate)?;
        fluid.accept();
        self.load = candidate;
        self.steps += 1;
        Ok(frame)
    }
}

fn load_error(error: WaveguideError) -> NsfDgError {
    NsfDgError::Input(error.0)
}

// Reconstruct the ACTUAL incident wave in the original outward volume-flow
// coordinates. Step a private Copy of the existing causal owner; no rewritten
// impedance integrator and no accepted memory change occurs during root trials.
fn candidate_load(
    mut load: RelaxationImpedance,
    port: f64,
    pressure: f64,
    flow: f64,
    series_pressure: f64,
) -> Result<(RelaxationImpedance, RelaxationFrame), WaveguideError> {
    if !port.is_finite()
        || port <= 0.0
        || ![pressure, flow, series_pressure]
            .iter()
            .all(|value| value.is_finite())
    {
        return Err(WaveguideError(
            "causal impedance requires finite physical port coordinates",
        ));
    }
    let incident = 0.5 * pressure - 0.5 * series_pressure + 0.5 * (port * flow);
    let frame = load.step(incident)?;
    Ok((load, frame))
}

fn matching_residual(
    face: NsfBoundaryExchange,
    load: RelaxationFrame,
    port: f64,
    series_pressure: f64,
) -> Result<(f64, f64), NsfDgError> {
    let residual = port * (face.normal_volume_flow_m3_s - load.port.flow_m3_s);
    let scale = face.mean_acoustic_pressure_pa.abs()
        + series_pressure.abs()
        + load.port.pressure_pa.abs()
        + port * (face.normal_volume_flow_m3_s.abs() + load.port.flow_m3_s.abs());
    if !residual.is_finite() || !scale.is_finite() {
        return Err(NsfDgError::Input(
            "causal impedance pressure equation overflowed",
        ));
    }
    Ok((residual, scale))
}

#[cfg(test)]
mod tests;
