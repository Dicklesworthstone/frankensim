//! Spatial supplied heat closes the original thermal-volume pressure work.
//!
//! Every trial uses the SAME original M/P/h and native pressure/material routes.
//! Only the accepted net caloric increments, transport and work are returned.

use super::{
    ChtError, FlowFace, FlowField, ThermalVolumeMomentumControls, ThermalVolumeMomentumInput,
    ThermalVolumeMomentumStep, VoxelDomain, poll,
};
use fs_exec::CancelGate;
use fs_material::phase::EquilibriumEnthalpyPhaseCurve;
use fs_math::dd::Dd;

/// Original implicit conduction and finite-solid storage inside heat/work closure.
#[cfg(feature = "implicit-enthalpy-coupling")]
pub mod implicit;

/// Native transport/pressure gates plus a bounded heat/work source iteration.
#[derive(Clone, Copy, Debug)]
pub struct SuppliedHeatThermalVolumeControls {
    pub native: ThermalVolumeMomentumControls,
    /// Maximum complete native thermal-volume/momentum proposals. Each retains
    /// its own original source, pressure, Krylov and transport work ceilings.
    pub maximum_heat_work_iterations: usize,
    /// Relaxation in (0,1] of NET caloric increments, not of a physical state.
    pub work_relaxation: f64,
    /// Maximum sum_i |represented Delta H_i + pressure work_i - supplied Q_i|.
    /// Opposing local defects cannot hide inside a zero global residual [J].
    pub maximum_absolute_local_heat_defect_j: f64,
    /// Maximum |actual inventory Delta H - carried boundary H + work - Q| [J].
    pub maximum_total_heat_defect_j: f64,
}

/// Accepted physical heat/work readback; discarded trials contribute no energy.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SuppliedHeatThermalVolumeReport {
    pub supplied_heat_j: f64,
    pub represented_net_enthalpy_increment_j: f64,
    /// Sum of the actual native correction-pressure work partition [J].
    pub source_pressure_work_j: f64,
    pub maximum_local_heat_defect_j: f64,
    pub absolute_local_heat_defect_j: f64,
    /// Partitioned work minus the ORIGINAL projection's work report [J].
    pub pressure_work_partition_defect_j: f64,
    /// Actual returned inventory change minus native carried boundary enthalpy,
    /// plus original source-pressure work, minus supplied physical heat [J].
    pub total_heat_defect_j: f64,
    pub heat_work_iterations: usize,
    /// Includes complete mixing-volume trials from every discarded heat guess.
    pub volume_source_iterations: usize,
}

/// One original conservative spatial endpoint with its supplied-heat closure.
#[derive(Clone, Debug, PartialEq)]
pub struct SuppliedHeatThermalVolumeStep {
    pub evolution: ThermalVolumeMomentumStep,
    /// Accepted requested NET cell increments [J], zero in stationary solids.
    pub net_enthalpy_increment_j: Vec<f64>,
    /// dt*p_c*S_c, signed from caloric storage into native mechanics [J].
    pub cell_pressure_work_j: Vec<f64>,
    pub report: SuppliedHeatThermalVolumeReport,
}

impl FlowField {
    /// Apply spatial physical heat while paying its native expansion work.
    ///
    /// The existing `input.volume.net_enthalpy_increment_j` is explicitly an
    /// INITIAL NET GUESS for this entry, not an additional heat source. The
    /// separate `supplied_heat_j` is the physical heat deposited into each
    /// original cell mass during the native interval. All trials restart the
    /// same original M/P/h and once-applied body force. The existing owner alone
    /// computes phase-density and enthalpy-mixing volume sources, pressure,
    /// impulse and conservative M/Mh/P transport.
    ///
    /// For each native trial, this entry enforces
    /// M_i*(h_star_i-h_old_i) + dt*p_i*S_i = Q_i. The actual represented caloric
    /// increment is measured through the original chart. S includes the
    /// original owner's nonlinear mixing-volume closure; its signed work is
    /// not omitted or deposited again as heat. Bounded relaxed source iteration
    /// changes only the guessed net increments. It is not another time march.
    ///
    /// Every fluid component MUST have an original zero-correction-pressure
    /// Free opening. The native component report enforces this, without a
    /// second topology map. The caller declares that this zero is the original
    /// chart's maintained reference pressure. Only correction-pressure work
    /// beyond the chart's reference p/rho is charged: ambient p0*dV is already
    /// included in h. A sealed or unanchored component refuses even when its
    /// total source is compatible; its arbitrary pressure gauge cannot set
    /// local caloric work. No thermodynamic pressure EOS is inferred.
    ///
    /// Held-flow/body-force/boundary work and projection, scatter and donor
    /// losses remain in their original mechanical accounts. None is recycled
    /// as heat. The returned native carried enthalpy is removed exactly once
    /// from the independent total heat audit; open-stream pressure-flow work
    /// is not charged again. All original source, phase-mobility, inventory,
    /// volume and numerical gates remain active.
    ///
    /// This is the original first-order fixed-grid thermal/source/transport
    /// split. Heat is held for this interval; spatial conduction, radiation or
    /// reservoir feedback must supply their own physical heat owner. It adds
    /// no free surface, solid stress, viscosity, arbitrary ALE, new phase
    /// kinetics or coupled entropy claim. A source-domain trial, stalled or
    /// exhausted iteration, failed physical budget or cancellation returns no
    /// candidate. Smaller timesteps or a supported initial guess may be needed;
    /// unsupported source trials are never clipped.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn advance_supplied_heat_thermal_volume_momentum(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        curve: &EquilibriumEnthalpyPhaseCurve,
        input: ThermalVolumeMomentumInput<'_>,
        supplied_heat_j: &[f64],
        controls: SuppliedHeatThermalVolumeControls,
        gate: &CancelGate,
    ) -> Result<SuppliedHeatThermalVolumeStep, ChtError> {
        let accepted = close_native_heat_work(
            domain,
            curve,
            input,
            supplied_heat_j,
            controls,
            gate,
            |net| {
                self.advance_piecewise_thermal_volume_momentum(
                    domain,
                    faces,
                    curve,
                    ThermalVolumeMomentumInput {
                        volume: super::CoexistenceVolumeInput {
                            net_enthalpy_increment_j: net,
                            ..input.volume
                        },
                        ..input
                    },
                    controls.native,
                    gate,
                )
                .map(|evolution| NativeHeatWorkTrial {
                    evolution,
                    conductive_heat_j: None,
                    receipt: (),
                })
            },
        )?;
        let NativeHeatWorkTrial {
            evolution,
            receipt: (),
            ..
        } = accepted.trial;
        Ok(SuppliedHeatThermalVolumeStep {
            evolution,
            net_enthalpy_increment_j: accepted.authored_net_enthalpy_increment_j,
            cell_pressure_work_j: accepted.cell_pressure_work_j,
            report: accepted.report,
        })
    }
}

// One concrete native trial. Optional Fourier heat belongs to the original
// implicit graph; its finite-solid/source receipt stays with this exact trial.
struct NativeHeatWorkTrial<T> {
    evolution: ThermalVolumeMomentumStep,
    conductive_heat_j: Option<Vec<f64>>,
    receipt: T,
}

struct NativeHeatWorkAcceptance<T> {
    trial: NativeHeatWorkTrial<T>,
    #[cfg(feature = "implicit-enthalpy-coupling")]
    absolute_local_heat_defect: Dd,
    authored_net_enthalpy_increment_j: Vec<f64>,
    cell_pressure_work_j: Vec<f64>,
    report: SuppliedHeatThermalVolumeReport,
}

// Both public entries use this SAME source-work closure and acceptance gates.
// A native thermal producer may refresh Fourier heat, never mass routes,
// pressure, source laws or a separate time trajectory.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn close_native_heat_work<T>(
    domain: &VoxelDomain,
    curve: &EquilibriumEnthalpyPhaseCurve,
    input: ThermalVolumeMomentumInput<'_>,
    supplied_heat_j: &[f64],
    controls: SuppliedHeatThermalVolumeControls,
    gate: &CancelGate,
    mut advance: impl FnMut(&[f64]) -> Result<NativeHeatWorkTrial<T>, ChtError>,
) -> Result<NativeHeatWorkAcceptance<T>, ChtError> {
    poll(gate)?;
    let n = domain.cell_count();
    if n > controls.native.volume.thermal.volume.transport.max_cells
        || supplied_heat_j.len() != n
        || input.volume.net_enthalpy_increment_j.len() != n
        || controls.maximum_heat_work_iterations == 0
        || !controls.work_relaxation.is_finite()
        || controls.work_relaxation <= 0.0
        || controls.work_relaxation > 1.0
        || [
            controls.maximum_absolute_local_heat_defect_j,
            controls.maximum_total_heat_defect_j,
        ]
        .iter()
        .any(|v| !v.is_finite() || *v < 0.0)
    {
        return Err(refuse(
            "invalid spatial heat/work input or iteration budget",
        ));
    }
    let mut supplied = Dd::ZERO;
    for (c, &heat) in supplied_heat_j.iter().enumerate() {
        poll_cell(gate, c)?;
        if !heat.is_finite() || (!domain.is_fluid(c) && heat != 0.0) {
            return Err(refuse(
                "physical heat must be finite and zero in solid grid cells",
            ));
        }
        supplied = supplied + Dd::from_f64(heat);
    }
    if !number(supplied).is_finite() {
        return Err(refuse("supplied spatial heat is unrepresentable"));
    }
    let dt = controls
        .native
        .volume
        .thermal
        .volume
        .transport
        .projection
        .duration_s;
    let mut net = input.volume.net_enthalpy_increment_j.to_vec();
    let mut next_net = vec![0.0; n];
    let mut work = vec![0.0; n];
    let mut volume_source_iterations = 0usize;
    for iteration in 1..=controls.maximum_heat_work_iterations {
        poll(gate)?;
        let trial = advance(&net)?;
        if trial
            .conductive_heat_j
            .as_ref()
            .is_some_and(|heat| heat.len() != n)
        {
            return Err(refuse(
                "native conduction returned an invalid fluid heat field",
            ));
        }
        let volume = &trial.evolution.thermal.volume;
        let projection = &volume.transport.projection;
        if projection.report.projection.pinned_components != 0 {
            return Err(refuse(
                "heat/work closure requires a reference-pressure Free opening in every fluid component",
            ));
        }
        volume_source_iterations = volume_source_iterations
            .checked_add(trial.evolution.thermal.source_iterations)
            .ok_or_else(|| refuse("spatial heat/work source count overflow"))?;
        let mut trial_supplied = supplied;
        let mut represented = Dd::ZERO;
        let mut partitioned_work = Dd::ZERO;
        let mut absolute_local = Dd::ZERO;
        let mut maximum_local = 0.0_f64;
        let mut changed = false;
        for c in 0..n {
            poll_cell(gate, c)?;
            let conductive = trial.conductive_heat_j.as_ref().map(|heat| heat[c]);
            if conductive
                .is_some_and(|heat| !heat.is_finite() || (!domain.is_fluid(c) && heat != 0.0))
            {
                return Err(refuse("native conduction returned invalid physical heat"));
            }
            if !domain.is_fluid(c) {
                continue;
            }
            let original = input.volume.phase_states[c]
                .ok_or_else(|| refuse("native heat/work trial lost its original phase state"))?;
            // Same original fixed-mass caloric advance as the native
            // volume owner; no subtraction of extensive absolute H.
            // Native implicit conduction assembles this same rounded sum
            // before the original caloric/volume advance. The physical heat
            // sum below retains its low part for the independent local gate.
            let (combined, physical_heat) = match conductive {
                Some(heat) => {
                    trial_supplied = trial_supplied + Dd::from_f64(heat);
                    (
                        net[c] + heat,
                        Dd::from_f64(supplied_heat_j[c]) + Dd::from_f64(heat),
                    )
                }
                None => (net[c], Dd::from_f64(supplied_heat_j[c])),
            };
            let heated = curve
                .advance_specific_energy(original, combined / input.volume.mass_kg[c])
                .map_err(|error| refuse(error.to_string()))?;
            let increment = Dd::from_f64(input.volume.mass_kg[c])
                * (Dd::from_f64(heated.specific_enthalpy_j_kg())
                    - Dd::from_f64(original.specific_enthalpy_j_kg()));
            let cell_work = Dd::from_f64(dt)
                * Dd::from_f64(projection.pressure_correction_pa[c])
                * Dd::from_f64(volume.volume_source_m3_s[c]);
            let residual = increment + cell_work - physical_heat;
            let residual_j = number(residual);
            let proposed =
                number(Dd::from_f64(net[c]) - Dd::from_f64(controls.work_relaxation) * residual);
            work[c] = number(cell_work);
            if !residual_j.is_finite() || !work[c].is_finite() || !proposed.is_finite() {
                return Err(refuse("spatial heat/work arithmetic is unrepresentable"));
            }
            represented = represented + increment;
            partitioned_work = partitioned_work + cell_work;
            absolute_local = absolute_local + magnitude(residual);
            maximum_local = maximum_local.max(residual_j.abs());
            next_net[c] = proposed;
            changed |= proposed.to_bits() != net[c].to_bits();
        }
        let work_defect = partitioned_work - Dd::from_f64(projection.report.pressure_source_work_j);
        let total_defect = Dd::from_f64(volume.report.material_enthalpy_change_j)
            - Dd::from_f64(volume.transport.report.carried_enthalpy_j)
            + Dd::from_f64(projection.report.pressure_source_work_j)
            - trial_supplied;
        if !number(represented).is_finite()
            || !number(partitioned_work).is_finite()
            || !number(work_defect).is_finite()
            || !number(absolute_local).is_finite()
            || !number(total_defect).is_finite()
            || !number(trial_supplied).is_finite()
        {
            return Err(refuse("spatial heat/work readback is unrepresentable"));
        }
        if exceeds(
            work_defect,
            controls.native.mechanics.projection.maximum_energy_defect_j,
        ) {
            return Err(refuse(
                "cell pressure work disagrees with the original native account",
            ));
        }
        if !exceeds(
            absolute_local,
            controls.maximum_absolute_local_heat_defect_j,
        ) {
            if exceeds(total_defect, controls.maximum_total_heat_defect_j) {
                return Err(refuse(
                    "actual spatial inventory, carried enthalpy and pressure work exceed the heat budget",
                ));
            }
            let report = SuppliedHeatThermalVolumeReport {
                supplied_heat_j: number(trial_supplied),
                represented_net_enthalpy_increment_j: number(represented),
                source_pressure_work_j: number(partitioned_work),
                maximum_local_heat_defect_j: maximum_local,
                absolute_local_heat_defect_j: number(absolute_local),
                pressure_work_partition_defect_j: number(work_defect),
                total_heat_defect_j: number(total_defect),
                heat_work_iterations: iteration,
                volume_source_iterations,
            };
            poll(gate)?;
            return Ok(NativeHeatWorkAcceptance {
                trial,
                #[cfg(feature = "implicit-enthalpy-coupling")]
                absolute_local_heat_defect: absolute_local,
                authored_net_enthalpy_increment_j: net,
                cell_pressure_work_j: work,
                report,
            });
        }
        if !changed || iteration == controls.maximum_heat_work_iterations {
            return Err(refuse("spatial heat/work iteration stalled or exhausted"));
        }
        core::mem::swap(&mut net, &mut next_net);
    }
    Err(refuse("spatial heat/work iteration exhausted"))
}

fn number(value: Dd) -> f64 {
    value.hi + value.lo
}

fn magnitude(value: Dd) -> Dd {
    if value.hi < 0.0 || (value.hi == 0.0 && value.lo < 0.0) {
        Dd::ZERO - value
    } else {
        value
    }
}

fn exceeds(value: Dd, budget: f64) -> bool {
    let excess = magnitude(value) - Dd::from_f64(budget);
    excess.hi > 0.0 || (excess.hi == 0.0 && excess.lo > 0.0)
}

fn poll_cell(gate: &CancelGate, cell: usize) -> Result<(), ChtError> {
    if cell % 256 == 0 {
        poll(gate)?;
    }
    Ok(())
}

fn refuse(reason: impl Into<String>) -> ChtError {
    ChtError::InvalidInput {
        field: "flow.supplied_heat_thermal_volume",
        reason: reason.into(),
    }
}

#[cfg(test)]
mod tests;
