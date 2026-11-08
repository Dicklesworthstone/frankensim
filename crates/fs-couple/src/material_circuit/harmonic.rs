//! Original material-resistor heat at a settled harmonic operating point.
//!
//! Native branch efforts are peak phasors with passive first-to-second current.
//! Their actual complex power is checked against the same original conductor
//! law used to lower the circuit resistance. No circuit solve, source refresh,
//! reference-resistor loss or second Joule law is introduced here.

use super::{
    MaterialCircuitError, MaterialCircuitSolve, MaterialResistor, admit_graph_shape, poll,
    positive, resistors::admit_resistor_sources,
};
use fs_em::planar::eddy::harmonic::circuit::HarmonicCircuitBranch;
use fs_material::conductor::OhmicDrive;
use fs_math::{c64::C64, dd::Dd};
use fs_phs::circuit::{Branch, CircuitGraph};

/// One original source-bound resistor's actual settled AC work and heat.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialHarmonicResistorHeat {
    /// Complete source, geometry and index in the original graph's branches.
    pub material: MaterialResistor,
    /// Actual first-minus-second terminal voltage, as a peak phasor [V].
    pub voltage_v: C64,
    /// Actual passive first-to-second current, as a peak phasor [A].
    pub current_a: C64,
    /// Cycle-mean nonnegative Joule power, including the peak factor 1/2 [W].
    pub mean_power_w: f64,
    /// Mean power times the declared settled thermal interval [J].
    /// This is not a phase-resolved finite-window integral.
    pub heat_j: f64,
}

/// Reconstruct each material resistor's cycle-mean heat from native AC efforts.
///
/// `graph` must be the actual source-lowered graph used for the native response.
/// Every original resistor needs one complete material binding; `branches`
/// contains every original graph branch once, without numerical reference
/// terminations. Inputs may be unordered; receipts return in original branch
/// order. Capacitors, inductors, transformers and sources retain their native
/// owners and produce no material-resistor heat here.
///
/// At the held material state, the selected source must be an instantaneous
/// scalar resistance. Field-, rate- and frequency-dependent sources refuse.
/// Temperature-dependent sources are allowed at a fixed resolved temperature;
/// their refresh and caloric destination remain the thermal owner's task.
/// Electrical series contacts require a separate work and heat owner.
///
/// Mean heat uses `(P(I.re) + P(I.im))/2` through the original conductor kernel.
/// Actual current, complex branch work and total native resistor loss are
/// independently checked. These absolute tolerances concern accounting, not
/// carrier settling, thermal integration or electromagnetic discretization.
#[allow(clippy::too_many_arguments)]
pub fn account_harmonic_resistor_heat(
    graph: &CircuitGraph,
    resistors: &[MaterialResistor],
    branches: &[HarmonicCircuitBranch],
    mean_resistor_power_w: f64,
    duration_s: f64,
    controls: MaterialCircuitSolve,
    current_tolerance_a: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Vec<MaterialHarmonicResistorHeat>, f64), MaterialCircuitError> {
    use MaterialCircuitError::{Balance, Invalid};
    poll(cancelled)?;
    if !positive(duration_s)
        || [
            controls.energy_tolerance_j,
            controls.power_tolerance_w,
            current_tolerance_a,
            mean_resistor_power_w,
        ]
        .iter()
        .any(|value| !value.is_finite() || *value < 0.0)
    {
        return Err(Invalid(
            "finite AC heat controls and positive interval required",
        ));
    }
    let count = admit_graph_shape(graph)?;
    if resistors.len() != count || branches.len() != graph.branches.len() {
        return Err(Invalid(
            "complete original AC material and branch receipts required",
        ));
    }
    // Reuse source geometry, resistance-scale, unique-binding and sort admission.
    // Each actual resistor must already contain this exact source resistance;
    // unrelated storage and field-placeholder numbers are not interpreted here.
    let sources = admit_resistor_sources(graph, resistors.to_vec())?;
    let mut by_branch = vec![None; graph.branches.len()];
    for branch in branches {
        poll(cancelled)?;
        let slot = by_branch
            .get_mut(branch.branch)
            .ok_or(Invalid("AC branch receipt is outside the original graph"))?;
        if slot.replace(branch).is_some() {
            return Err(Invalid("duplicate original AC branch receipt"));
        }
        if !finite(branch.voltage_v)
            || !finite(branch.current_a)
            || !finite(branch.absorbed_power_va)
        {
            return Err(Balance("nonfinite native AC branch receipt"));
        }
    }
    let mut receipts = Vec::with_capacity(count);
    let mut power = Dd::ZERO;
    let mut heat = Dd::ZERO;
    for source in sources {
        poll(cancelled)?;
        let Branch::Resistor { ohms } = graph.branches[source.branch].2 else {
            unreachable!("the original material owner admits resistor bindings")
        };
        if ohms.to_bits() != source.conductor.resistance_ohm().to_bits() {
            return Err(Invalid(
                "native AC graph resistance differs from its material source",
            ));
        }
        if !source
            .conductor
            .source_property()
            .is_static_at_fixed_external_conditions()
        {
            return Err(Invalid(
                "AC resistor heat requires an instantaneous scalar source",
            ));
        }
        let branch =
            by_branch[source.branch].ok_or(Invalid("missing original AC resistor receipt"))?;
        let material_current = C64::new(
            source
                .conductor
                .current_a_for_drive(OhmicDrive::Voltage(branch.voltage_v.re))?,
            source
                .conductor
                .current_a_for_drive(OhmicDrive::Voltage(branch.voltage_v.im))?,
        );
        let current_defect = (material_current - branch.current_a).abs();
        if !current_defect.is_finite() || current_defect > current_tolerance_a {
            return Err(Balance(
                "native AC resistor current and material law disagree",
            ));
        }
        // midpoint avoids overflow of the sum and loss of two subnormal halves.
        let mean_power_w = f64::midpoint(
            source.conductor.joule_power_w(branch.current_a.re)?,
            source.conductor.joule_power_w(branch.current_a.im)?,
        );
        let native_power = (branch.voltage_v * branch.current_a.conj()).scale(0.5);
        let defects = [
            native_power.re - mean_power_w,
            native_power.im,
            branch.absorbed_power_va.re - mean_power_w,
            branch.absorbed_power_va.re - native_power.re,
            branch.absorbed_power_va.im - native_power.im,
            branch.absorbed_power_va.im,
        ];
        let heat_j = duration_s * mean_power_w;
        if !mean_power_w.is_finite()
            || mean_power_w < 0.0
            || (branch.current_a != C64::ZERO && mean_power_w == 0.0)
            || !heat_j.is_finite()
            || (mean_power_w != 0.0 && heat_j == 0.0)
            || branch.absorbed_power_va.re < 0.0
            || defects
                .iter()
                .any(|defect| !defect.is_finite() || defect.abs() > controls.power_tolerance_w)
        {
            return Err(Balance(
                "AC resistor complex work and material heat disagree",
            ));
        }
        power = power + Dd::from_f64(mean_power_w);
        heat = heat + Dd::from_f64(heat_j);
        receipts.push(MaterialHarmonicResistorHeat {
            material: source,
            voltage_v: branch.voltage_v,
            current_a: branch.current_a,
            mean_power_w,
            heat_j,
        });
    }
    let native_heat_j = duration_s * mean_resistor_power_w;
    let power_defect = (power - Dd::from_f64(mean_resistor_power_w)).to_f64();
    let heat_defect = (heat - Dd::from_f64(native_heat_j)).to_f64();
    let heat_j = heat.to_f64();
    if !power.to_f64().is_finite()
        || !heat_j.is_finite()
        || !native_heat_j.is_finite()
        || (mean_resistor_power_w != 0.0 && native_heat_j == 0.0)
        || !power_defect.is_finite()
        || !heat_defect.is_finite()
        || power_defect.abs() > controls.power_tolerance_w
        || heat_defect.abs() > controls.energy_tolerance_j
    {
        return Err(Balance(
            "summed AC material heat and native resistor loss disagree",
        ));
    }
    poll(cancelled)?;
    Ok((receipts, heat_j))
}

fn finite(value: C64) -> bool {
    value.re.is_finite() && value.im.is_finite()
}
