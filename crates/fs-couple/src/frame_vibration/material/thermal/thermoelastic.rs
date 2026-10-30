//! Source E(T), expansion and finite heat on the original complete frame owner.
//!
//! The fixed, complete mass basis is only a coordinate chart. Its kinetic
//! coordinates use fs-modal's exact free-coordinate owner; each original layer
//! supplies three positive physical springs to the SAME coupled propagator.
//! Axial rest extension is that layer's integrated alpha dT times L. Separate
//! bending factors retain the full off-diagonal stiffness and mismatch energy
//! of layers at unequal temperatures. Original material dashpots stay separate.
//!
//! An interval first follows the exact accepted-temperature mechanics. Its
//! original spring-port work W and dashpot losses D then enter the existing
//! finite-cell residual, with Q_rev = W - delta(Psi + T S). Endpoint stiffness,
//! rest extension and viscosity are prepared at unchanged endpoint coordinates
//! and installed only with the matching thermal state. This documented
//! first-order split needs time refinement; conservation alone is not a
//! temporal-accuracy or reversible-entropy certificate.
//!
//! Geometry/reference density stay fixed. Scalar expansion and scalar or
//! globally affine E(T) use the original ConcentricThermoelasticLaw. Every
//! layer has a whole solid caloric cell, positive combined capacity and an
//! explicit small-strain envelope. Radial constraint, buckling/geometric
//! stiffness, changing mass, topology and phase transitions require other
//! owners. No residual-energy correction or second time integrator is used.
//!
//! [`FrameVibrationComponent::into_thermoelastic_parts`] retains this physical
//! source when the frame joins other native mechanical components, reciprocal
//! acoustic loads and constitutive memories. Multiple bindings join through
//! [`FrameThermoelasticBinding::into_shared_system`]: ONE native mechanical step,
//! ONE complete thermal residual and ONE endpoint coefficient replacement.
//! Each frame's original dashpots heat only its own cells. Cross-frame thermal
//! links use the same source network; unbound losses remain separate sinks.
//! [`FrameThermoelasticSystem::with_material_circuit`] adds source-bound R(T)
//! and native finite electrical storage to that SAME thermal residual. Joule
//! heat, material dissipation and reversible exchange remain distinct accounts.

mod assembly;
mod connection_heat;
mod electrical;
pub use assembly::{FrameThermoelasticBinding, FrameThermoelasticFrame};
pub use connection_heat::{FrameConnectionHeat, FrameConnectionHeatDestination};
pub use electrical::FrameElectricalInput;

use super::*;
use crate::modal_acoustic_time::{ModalAcousticMode, ModalAcousticState, ModalAcousticTimeModel};
use crate::render::schedule::force::coupled::{
    ModalConnection, coefficients::ConnectionCoefficients,
};
use fs_material::section_temperature::thermoelastic::{
    ConcentricThermoelasticLaw, ConcentricThermoelasticState,
};
use fs_math::c64::C64;

struct Member<'a> {
    region: Region<'a>,
    law: ConcentricThermoelasticLaw<'a>,
    columns: [Vec<f64>; 3],
    spring_start: usize,
    component: usize,
    damping_start: usize,
    linear_strain_limit: f64,
    axial_property: String,
    bending_property: String,
}

/// One physical layer's mechanical contribution, additional to its caloric U.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FrameThermoelasticLayerState {
    pub destination: FrameLayerCell,
    pub temperature_k: f64,
    pub free_energy_j: f64,
    pub entropy_j_k: f64,
    pub internal_energy_j: f64,
    pub maximum_absolute_elastic_strain: f64,
}

/// Complete pinned source evaluation at the actual accepted layer temperatures.
#[derive(Clone, Debug)]
pub struct FrameThermoelasticMemberState {
    /// Original native component; element numbers are local to this frame.
    pub component: usize,
    pub element: usize,
    pub source: ConcentricThermoelasticState,
    pub viscosity: FrameViscousSection,
    pub layers: Vec<FrameThermoelasticLayerState>,
}

/// Physical work and reversible exchange of one whole original layer.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FrameThermoelasticLayerExchange {
    pub states: [FrameThermoelasticLayerState; 2],
    /// Original spring-port work into elastic storage during frozen-T motion.
    pub stress_work_j: f64,
    /// Original eta*A / eta*I partition of independent dashpot work.
    pub material_heat_j: f64,
    /// Psi(T_new, q_end) - Psi(T_old, q_end), at unchanged physical coordinates.
    pub coefficient_free_energy_change_j: f64,
    pub mechanical_entropy_energy_change_j: f64,
    /// Signed caloric exchange W_stress - delta(U_material), not dissipation.
    pub reversible_heat_j: f64,
}

/// Original caloric entropy plus the mechanical source entropy, in J/K.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FrameThermoelasticEntropy {
    pub thermal_storage_change_j_k: f64,
    pub mechanical_storage_change_j_k: f64,
    pub reservoir_entropy_input_j_k: f64,
    pub conduction_generation_j_k: f64,
    /// All deposited material/connector loss and, when bound, Joule loss over
    /// the actual receiving endpoint Kelvin. Reversible heat is separate.
    pub material_loss_generation_j_k: f64,
    pub unclassified_source_at_cell_j_k: f64,
    pub has_unclassified_sources: bool,
    /// Sum Q_rev/T_endpoint; the mechanical entropy balances its continuum limit.
    pub reversible_transfer_entropy_j_k: f64,
    pub thermal_temporal_defect_j_k: f64,
    /// Delta S_mechanical + sum Q_rev/T_endpoint, a signed splitting error.
    pub thermoelastic_temporal_defect_j_k: f64,
    pub algebraic_residual_j_k: f64,
    pub accounting_residual_j_k: f64,
}

/// One jointly accepted sample; native mechanics describes its frozen-T interval.
pub struct FrameThermoelasticReport {
    pub mechanical: CoupledModalFrame,
    /// Original electrical work, finite storage and branch heat, when bound.
    pub electrical: Option<crate::material_circuit::MaterialCircuitStep>,
    pub thermal: EnergyAudit,
    pub source: MaterialNetworkSample,
    pub members: Vec<FrameThermoelasticMemberState>,
    pub layers: Vec<FrameThermoelasticLayerExchange>,
    /// Original non-frame connection losses routed into explicit caloric cells.
    pub connection_heat: Vec<FrameConnectionHeat>,
    pub entropy: FrameThermoelasticEntropy,
    pub parameter_work_j: f64,
    pub reversible_heat_j: f64,
    /// Source-frame damping plus explicitly routed connector losses, counted once.
    pub deposited_material_heat_j: f64,
    /// Native circuit heat, separate from mechanical losses and reversible heat.
    pub deposited_joule_heat_j: f64,
    pub external_thermal_heat_j: f64,
    pub external_mechanical_loss_j: f64,
    /// Includes kinetics and source Psi+TS, excluding caloric and circuit storage.
    pub mechanical_internal_energy_before_j: f64,
    pub mechanical_internal_energy_after_j: f64,
    /// Complete mechanical, electrical and finite-thermal balance when bound.
    pub total_energy_residual_j: f64,
    pub nonlinear: NewtonReport,
    pub material_trials: usize,
    pub maximum_cell_residual_j: f64,
}

/// Original mechanical history, fixed physical chart and finite-cell sources.
pub struct FrameThermoelasticSystem<'n, 'a> {
    frames: Vec<FrameThermoelasticFrame>,
    /// Sorted disjoint source-owned connection ranges, prepared once.
    owned_connections: Vec<core::ops::Range<usize>>,
    connection_heat_destinations: Vec<FrameConnectionHeatDestination>,
    mechanical: CoupledModalSystem,
    prepared: PreparedMaterialNetwork<'n, 'a>,
    source: &'n MaterialNetwork<'a>,
    electrical: Option<electrical::Binding<'n, 'a>>,
    members: Vec<Member<'a>>,
    states: Vec<FrameThermoelasticMemberState>,
    origin_s: f64,
}

fn modal(error: crate::modal_acoustic_time::ModalAcousticTimeError) -> FrameThermalError {
    CircuitEnthalpyError::Source(FrameThermalSourceError::Frame(error.into()))
}
fn solid(error: fs_solid::SolidError) -> FrameThermalError {
    CircuitEnthalpyError::Source(FrameThermalSourceError::Frame(error.into()))
}
fn sum(values: impl IntoIterator<Item = f64>) -> Result<f64> {
    values.into_iter().try_fold(0.0, |s, x| finite(s + x))
}
fn positive(value: f64) -> Result<f64> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(invalid(
            "positive representable frame thermoelastic coefficient required",
        ))
    }
}
fn temperatures(region: &Region<'_>, thermal: &DepositionSimulation<'_>) -> Result<Vec<f64>> {
    region
        .receivers
        .iter()
        .map(|receiver| {
            let cell = thermal
                .cells()
                .get(&receiver.binding.cell)
                .ok_or_else(|| invalid("frame thermoelastic receiver disappeared"))?;
            let original = &receiver.original;
            if cell.id != original.id
                || cell.material_id != original.material_id
                || cell.segment_id != original.segment_id
                || cell.body_id != original.body_id
                || cell.mass_kg.to_bits() != original.mass_kg.to_bits()
                || cell.geometry != original.geometry
            {
                return Err(invalid(
                    "frame thermoelastic receiver changed its original inventory or geometry",
                ));
            }
            thermal.temperature_k(cell.id).map_err(Into::into)
        })
        .collect()
}
fn resolve_member(
    member: &Member<'_>,
    temperatures: &[f64],
    q: &[ModalAcousticState],
    dt: f64,
) -> Result<FrameThermoelasticMemberState> {
    let strain_limit = member.linear_strain_limit;
    let mut x = [0.0; 3];
    for (value, column) in x.iter_mut().zip(&member.columns) {
        *value = sum(column
            .iter()
            .zip(q)
            .map(|(b, q)| b * q.displacement_m_sqrt_kg))?;
    }
    let length = member.region.length_m;
    let strain = finite(x[0] / length)?;
    let source = member
        .law
        .at_temperatures_k(temperatures, strain, strain_limit)
        .map_err(material)?;
    let curvature =
        finite((3.0 * x[1] + x[2]).abs().max((3.0 * x[1] - x[2]).abs()) / length / length)?;
    let curvature_integral = finite((3.0 * x[1] * x[1] + x[2] * x[2]) / length / length / length)?;
    let mut layers = Vec::with_capacity(source.layers.len());
    for (i, layer) in source.layers.iter().enumerate() {
        let radius = source.section.section().layers()[i].outer_radius_m;
        let maximum_absolute_elastic_strain =
            finite(layer.elastic_strain.abs() + radius * curvature)?;
        if maximum_absolute_elastic_strain > strain_limit
            || finite(strain.abs() + radius * curvature)? > strain_limit
        {
            return Err(invalid(
                "frame endpoint exceeds its total or elastic fibre-strain envelope",
            ));
        }
        let bending_geometry = finite(layer.second_moment_m4 * curvature_integral)?;
        let free_energy_j = finite(
            length * layer.axial_free_energy_per_reference_length_j_m
                + 0.5 * layer.young_modulus_pa * bending_geometry,
        )?;
        let entropy_j_k = finite(
            length * layer.entropy_per_reference_length_j_m_k
                - 0.5 * layer.young_modulus_derivative_pa_k * bending_geometry,
        )?;
        layers.push(FrameThermoelasticLayerState {
            destination: member.region.receivers[i].binding,
            temperature_k: layer.temperature_k,
            free_energy_j,
            entropy_j_k,
            internal_energy_j: finite(free_energy_j + layer.temperature_k * entropy_j_k)?,
            maximum_absolute_elastic_strain,
        });
    }
    let viscosity = FrameViscousSection::try_new(
        source.section.section().clone(),
        &member.axial_property,
        &member.bending_property,
    )
    .map_err(solid)?;
    let band = viscosity.omega_band_rad_s();
    if band.0 != 0.0 || band.1 < core::f64::consts::PI / dt {
        return Err(invalid(
            "thermoelastic frame sources must support DC through sampled Nyquist",
        ));
    }
    Ok(FrameThermoelasticMemberState {
        component: member.component,
        element: member.region.element,
        source,
        viscosity,
        layers,
    })
}
fn layer_springs(
    state: &FrameThermoelasticMemberState,
    layer: usize,
    length: f64,
) -> Result<[(f64, f64); 3]> {
    let value = &state.source.layers[layer];
    let axial = positive(value.young_modulus_pa * value.area_m2 / length)?;
    let bending =
        positive(value.young_modulus_pa * value.second_moment_m4 / length / length / length)?;
    let free_strain = state.source.section.expansion_paths()[layer].free_linear_strain();
    let rest = finite(length * free_strain)?;
    if free_strain != 0.0 && rest == 0.0 {
        return Err(invalid("thermal layer rest extension underflow"));
    }
    Ok([
        (axial, rest),
        (positive(3.0 * bending)?, 0.0),
        (bending, 0.0),
    ])
}
fn admit_capacity(
    member: &Member<'_>,
    state: &FrameThermoelasticMemberState,
    source: &MaterialNetwork<'_>,
    strain_limit: f64,
) -> Result<()> {
    for (receiver, layer) in member.region.receivers.iter().zip(&state.source.layers) {
        let caloric = source
            .cells
            .iter()
            .find(|s| s.cell == receiver.binding.cell)
            .ok_or_else(|| invalid("missing original layer caloric chart"))?;
        let knots = caloric.curve.knots();
        let mut capacity = f64::INFINITY;
        for pair in knots.windows(2) {
            capacity = capacity.min(positive(
                receiver.original.mass_kg
                    * (pair[1].specific_enthalpy_j_kg - pair[0].specific_enthalpy_j_kg)
                    / positive(pair[1].temperature_k - pair[0].temperature_k)?,
            )?);
        }
        let (lo, hi) = (knots[0].temperature_k, knots[knots.len() - 1].temperature_k);
        // The original support thermal owner uses this same conservative bound.
        // Affine arithmetic bounds Psi_TT; actual source queries never extrapolate.
        let modulus = |t: f64| {
            finite(
                layer
                    .young_modulus_derivative_pa_k
                    .mul_add(t - layer.temperature_k, layer.young_modulus_pa),
            )
        };
        let maximum_modulus = modulus(lo)?.max(modulus(hi)?).max(0.0);
        let curvature = finite(
            (1.0 + 2048.0 * f64::EPSILON)
                * receiver.original.geometry.volume_m3
                * (maximum_modulus * layer.alpha_per_k.powi(2)
                    + 2.0
                        * (layer.young_modulus_derivative_pa_k * layer.alpha_per_k).abs()
                        * strain_limit),
        )?;
        positive(finite(capacity - hi * curvature)?)?;
    }
    Ok(())
}

impl FrameVibrationComponent {
    /// Bind complete source-defined elasticity, expansion and viscosity to the
    /// same whole-layer caloric cells. Restored q/v are retained exactly.
    ///
    /// Every original member must have its explicit material viscosity binding
    /// (zero eta is valid), and ALL original free modes must be retained. The
    /// fixed basis therefore spans the complete unchanged mass space: changed
    /// E(T) may couple arbitrary modes without a truncation/reprojection claim.
    /// Frozen preload and additional authored modal damping currently refuse.
    /// Each layer requires explicit scalar alpha, fixed scalar rho, scalar or
    /// globally affine E(T), a solid caloric chart and positive combined capacity.
    /// `linear_strain_limit` bounds total, thermal and extreme fibre strains.
    #[allow(clippy::too_many_arguments)]
    pub fn into_thermoelastic_system<'n, 'a>(
        self,
        thermal: &DepositionSimulation<'_>,
        source: &'n MaterialNetwork<'a>,
        destinations: &[FrameLayerCell],
        linear_strain_limit: f64,
        config: ModalCouplingConfig,
        prepared: Option<PreparedModalBudget>,
        gate: &CancelGate,
    ) -> Result<FrameThermoelasticSystem<'n, 'a>> {
        let (binding, model, connections) = self.into_thermoelastic_parts(
            thermal,
            source,
            destinations,
            linear_strain_limit,
            0,
            0,
            config,
            gate,
        )?;
        let mechanical = match prepared {
            Some(budget) => {
                CoupledModalSystem::new_prepared(vec![model], connections, config, budget, gate)?
            }
            None => CoupledModalSystem::new(vec![model], connections, config, gate)?,
        };
        binding.into_system(mechanical, thermal, gate)
    }

    /// Prepare the SAME source-bound frame for a larger original coupled owner.
    /// Append its returned model at `component` and EVERY returned connection
    /// contiguously at `connection_start`, in unchanged order. Other components
    /// and reciprocal connectors may precede or follow those addresses.
    /// Configure native radiation, material memories and propagation normally,
    /// then call [`FrameThermoelasticBinding::into_system`] or collect bindings
    /// in [`FrameThermoelasticBinding::into_shared_system`]. No extra integrator
    /// or acoustic pressure transfer is supplied by these mechanical parts.
    /// The explicit config bounds this source projection; the final native
    /// constructor independently admits the complete assembly and its budgets.
    #[allow(clippy::too_many_arguments)]
    pub fn into_thermoelastic_parts<'n, 'a>(
        self,
        thermal: &DepositionSimulation<'_>,
        source: &'n MaterialNetwork<'a>,
        destinations: &[FrameLayerCell],
        linear_strain_limit: f64,
        component: usize,
        connection_start: usize,
        config: ModalCouplingConfig,
        gate: &CancelGate,
    ) -> Result<(
        FrameThermoelasticBinding<'n, 'a>,
        ModalAcousticTimeModel,
        Vec<ModalConnection>,
    )> {
        check(gate)?;
        if !linear_strain_limit.is_finite()
            || linear_strain_limit <= 0.0
            || self.basis.modes().len() != self.basis.available_mode_count()
            || self
                .basis
                .reference_displacement()
                .iter()
                .any(|&u| u != 0.0)
            || self.model.modes().iter().any(|m| m.damping_ratio != 0.0)
            || self
                .basis
                .material_damping()
                .is_none_or(|m| !m.covers_all_members())
        {
            return Err(invalid(
                "thermoelastic frame requires complete source-bound mass basis, zero frozen preload, explicit strain limit and no authored modal damping",
            ));
        }
        let (binding, regions, viscosity) =
            bind_regions(&self.basis, thermal, source, destinations, false, gate)?;
        let FrameVibrationComponent {
            basis,
            model,
            mut connections,
        } = self;
        let count = model.modes().len();
        let dt = model.sample_period_s();
        let total_ports = destinations
            .len()
            .checked_mul(3)
            .and_then(|n| n.checked_add(connections.len()))
            .ok_or_else(|| invalid("thermoelastic layer port count overflow"))?;
        if count > config.max_modes
            || component >= config.max_modes
            || connection_start
                .checked_add(total_ports)
                .is_none_or(|n| n > config.max_connections)
            || total_ports
                .checked_mul(count)
                .and_then(|n| n.checked_mul(count))
                .is_none_or(|n| n > config.max_setup_terms)
        {
            return Err(invalid(
                "thermoelastic complete layer projection exceeds native setup budget",
            ));
        }
        let rate = (1.0 / dt).round();
        if rate < 1.0 || rate > f64::from(u32::MAX) {
            return Err(invalid("unrepresentable native frame clock"));
        }
        let mut free = ModalAcousticTimeModel::try_new_with_free_coordinates(
            rate as u32,
            vec![
                ModalAcousticMode {
                    angular_frequency_rad_s: 0.0,
                    damping_ratio: 0.0,
                    pressure_per_modal_velocity: C64::ZERO
                };
                count
            ],
            model.budget(),
        )
        .map_err(modal)?;
        if free.sample_period_s().to_bits() != dt.to_bits() {
            return Err(invalid("native frame clock changed"));
        }
        free.restore_states(model.states()).map_err(modal)?;
        let mut members = Vec::with_capacity(regions.len());
        let mut states = Vec::with_capacity(regions.len());
        let mut projected = vec![0.0; count * count];
        for (index, (region, (_, eta))) in regions.into_iter().zip(viscosity).enumerate() {
            check(gate)?;
            let member = Member {
                law: ConcentricThermoelasticLaw::try_new(region.law.clone()).map_err(material)?,
                columns: core::array::from_fn(|k| connections[3 * index + k].left.shapes.clone()),
                spring_start: connections.len(),
                component,
                damping_start: connection_start + 3 * index,
                linear_strain_limit,
                axial_property: eta.axial().viscosity_property().into(),
                bending_property: eta.bending().viscosity_property().into(),
                region,
            };
            let state = resolve_member(
                &member,
                &temperatures(&member.region, thermal)?,
                model.states(),
                dt,
            )?;
            admit_capacity(&member, &state, source, linear_strain_limit)?;
            for layer in 0..state.layers.len() {
                for (k, (stiffness, rest)) in layer_springs(&state, layer, member.region.length_m)?
                    .into_iter()
                    .enumerate()
                {
                    let mut connection: ModalConnection = connections[3 * index + k].clone();
                    connection.stiffness_n_m = stiffness;
                    connection.damping_n_s_m = 0.0;
                    connection.rest_extension_m = rest;
                    for row in 0..count {
                        for col in 0..count {
                            projected[row * count + col] = finite(
                                projected[row * count + col]
                                    + stiffness * member.columns[k][row] * member.columns[k][col],
                            )?;
                        }
                    }
                    connections.push(connection);
                }
            }
            members.push(member);
            states.push(state);
        }
        // Independently reconstruct the original native K in its COMPLETE Phi.
        // This also prevents counting the source elasticity twice.
        let scale = model
            .modes()
            .iter()
            .map(|m| m.angular_frequency_rad_s.powi(2))
            .fold(0.0, f64::max);
        let tolerance = finite(
            4096.0 * f64::EPSILON * count as f64 * scale
                + basis
                    .modes()
                    .iter()
                    .map(|m| m.eigenvalue_resolution_per_s2)
                    .sum::<f64>(),
        )?;
        for row in 0..count {
            for col in 0..count {
                let expected = if row == col {
                    model.modes()[row].angular_frequency_rad_s.powi(2)
                } else {
                    0.0
                };
                if finite(projected[row * count + col] - expected)?.abs() > tolerance {
                    return Err(invalid(
                        "layer elastic ports do not reconstruct the complete original frame stiffness",
                    ));
                }
            }
        }
        for link in &mut connections {
            link.left.component = component;
            link.right.component = component;
        }
        for member in &mut members {
            member.spring_start += connection_start;
        }
        check(gate)?;
        let binding = FrameThermoelasticBinding {
            basis,
            prepared: binding,
            source,
            members,
            states,
            origin_s: thermal.time_s(),
            component,
            connection_start,
            initial_model: free.clone(),
            definitions: connections.clone(),
        };
        Ok((binding, free, connections))
    }
}

impl FrameThermoelasticSystem<'_, '_> {
    /// All source frames in canonical native-component order. Each range indexes
    /// the same flattened member states returned by `members()` and step reports.
    pub fn frames(&self) -> &[FrameThermoelasticFrame] {
        &self.frames
    }
    /// Immutable mass chart of the lowest-index bound frame. Existing single-
    /// frame calls retain their meaning; use `frames()` for a shared assembly.
    /// Reference eigenfrequencies are not the evolving elastic operator.
    pub const fn basis(&self) -> &FrameVibrationBasis {
        &self.frames.as_slice()[0].basis
    }
    pub const fn mechanical(&self) -> &CoupledModalSystem {
        &self.mechanical
    }
    /// Lowest-index bound component; see `frames()` for the complete assembly.
    pub const fn component_index(&self) -> usize {
        self.frames.as_slice()[0].component
    }
    /// First material connection of the lowest-index bound component.
    pub const fn connection_start(&self) -> usize {
        self.frames.as_slice()[0].connections.start
    }
    pub fn members(&self) -> &[FrameThermoelasticMemberState] {
        &self.states
    }
    pub fn mechanical_internal_energy_j(&self) -> Result<f64> {
        finite(
            self.mechanical.total_energy_j()?
                + sum(self
                    .states
                    .iter()
                    .flat_map(|m| m.layers.iter().map(|s| s.temperature_k * s.entropy_j_k)))?,
        )
    }

    /// Exact frozen-start-temperature motion, finite caloric solve and endpoint
    /// source replacement on one atomic boundary. The new physical K, rest and
    /// eta act NEXT interval; q/v, port IDs and sample history never reset.
    /// Independent caloric and Psi+TS accounts close all supplied energy.
    /// Failed source/strain/capacity/thermal/native update or cancellation
    /// publishes neither history and permits deterministic retry.
    /// An electrically bound assembly requires `step_electrical` instead, so
    /// ordinary stepping cannot silently freeze its charge and flux history.
    pub fn step(
        &mut self,
        thermal: &mut DepositionSimulation<'_>,
        external_modal_forces: &[f64],
        controls: FrameThermalSolve,
        energy_guess_j: &[f64],
        gate: &CancelGate,
    ) -> Result<FrameThermoelasticReport> {
        self.step_inner(thermal, external_modal_forces, None, controls, energy_guess_j, gate)
    }

    #[allow(clippy::too_many_arguments)]
    fn step_inner(
        &mut self,
        thermal: &mut DepositionSimulation<'_>,
        external_modal_forces: &[f64],
        electrical_input: Option<FrameElectricalInput<'_>>,
        controls: FrameThermalSolve,
        energy_guess_j: &[f64],
        gate: &CancelGate,
    ) -> Result<FrameThermoelasticReport> {
        check(gate)?;
        if self.electrical.is_some() != electrical_input.is_some() {
            return Err(invalid("electrically bound frame histories require explicit held electrical inputs"));
        }
        let sample = self.mechanical.samples_rendered();
        let next = sample
            .checked_add(1)
            .ok_or_else(|| invalid("thermoelastic frame clock overflow"))?;
        let dt = self.mechanical.sample_period_s();
        let start = finite(self.origin_s + sample as f64 * dt)?;
        let target = finite(self.origin_s + next as f64 * dt)?;
        if thermal.time_s().to_bits() != start.to_bits() || target <= start {
            return Err(invalid(
                "thermoelastic frame and caloric histories require their same accepted clock",
            ));
        }
        let dt_heat = target - start;
        let before = self.mechanical.total_energy_j()?;
        let internal_before = self.mechanical_internal_energy_j()?;
        let members = &self.members;
        let states = &self.states;
        let prepared = &self.prepared;
        let electrical = &self.electrical;
        let owned_connections = &self.owned_connections;
        let connection_heat_destinations = &self.connection_heat_destinations;
        for (member, state) in members.iter().zip(states) {
            for (temperature, old) in temperatures(&member.region, thermal)?
                .iter()
                .zip(&state.layers)
            {
                if temperature.to_bits() != old.temperature_k.to_bits() {
                    return Err(invalid(
                        "thermoelastic source and caloric temperatures must remain jointly accepted",
                    ));
                }
            }
        }
        let (next_thermal, report, update, next_states, next_circuit) = self.mechanical.step_with_admission(
            external_modal_forces, gate, |mechanical, frame, transfers| {
                let mut frozen = Vec::with_capacity(members.len());
                let mut works = Vec::new();
                let mut losses = Vec::new();
                let mut powers = BTreeMap::new();
                let mut material_heat_j = 0.0;
                for (index, (member, old)) in members.iter().zip(states).enumerate() {
                    check(gate)?;
                    let endpoint = mechanical.staged_components()[member.component].states();
                    frozen.push(resolve_member(member, &old.layers.iter().map(|l| l.temperature_k).collect::<Vec<_>>(),
                        endpoint, dt)?);
                    let port_loss = |k: usize| -> Result<f64> {
                        let loss = finite(-transfers[member.damping_start+k].energy_into().value())?;
                        if loss < 0.0 { return Err(invalid("frame dashpot supplied negative material heat")); }
                        Ok(loss)
                    };
                    let axial_loss = port_loss(0)?;
                    let bending_loss = finite(port_loss(1)?+port_loss(2)?)?;
                    material_heat_j = finite(material_heat_j+axial_loss+bending_loss)?;
                    for (layer, receiver) in member.region.receivers.iter().enumerate() {
                        let eta = &old.viscosity;
                        let partition = |fraction: f64, loss: f64| -> Result<f64> {
                            let part = finite(fraction*loss)?;
                            if fraction > 0.0 && loss > 0.0 && part == 0.0 {
                                return Err(invalid("physical frame layer heat underflow"));
                            }
                            Ok(part)
                        };
                        let loss = finite(partition(share(eta.axial().layer_viscous_rigidity_n_s()[layer],
                            eta.axial().viscous_rigidity_n_s())?, axial_loss)?
                            + partition(share(eta.bending().layer_viscous_stiffness_n_m2_s()[layer],
                                eta.bending().viscous_stiffness_n_m2_s())?, bending_loss)?)?;
                        let first = member.spring_start+3*layer;
                        let work = sum(transfers[first..first+3].iter().map(|p| -p.energy_into().value()))?;
                        let delta = frozen[index].layers[layer].free_energy_j-old.layers[layer].free_energy_j;
                        if finite(work-delta)?.abs() > controls.total_energy_tolerance_j {
                            return Err(invalid("native layer spring work disagrees with its source potential"));
                        }
                        let power = finite(loss/dt_heat)?;
                        if loss > 0.0 && power == 0.0 { return Err(invalid("frame loss power underflow")); }
                        powers.insert(receiver.binding.cell, power);
                        works.push(work); losses.push(loss);
                    }
                }
                let loss_account = assembly::outside_loss(mechanical, frame, transfers,
                    owned_connections, connection_heat_destinations, material_heat_j,
                    controls.total_energy_tolerance_j, gate)?;
                // These are actual connection losses, not spring work, an
                // energy residual, or a second copy of source-frame damping.
                for receipt in &loss_account.receivers {
                    check(gate)?;
                    let power = finite(receipt.heat_j / dt_heat)?;
                    if receipt.heat_j > 0.0 && power == 0.0 {
                        return Err(invalid("routed connection heat power underflow"));
                    }
                    let source = powers.entry(receipt.cell).or_insert(0.0);
                    *source = finite(*source + power)?;
                }
                let material_heat_j = finite(material_heat_j + loss_account.routed_j)?;
                let outside = loss_account.outside_j;
                let (candidate, next_circuit, audit, entropy) = cell_driver_candidate_entropy_option(
                    thermal, target, controls.driver(), energy_guess_j,
                    |trial, _remaining| {
                        check(gate)?;
                        let (mut network_source, circuit) = electrical::resolve(
                            electrical.as_ref(), electrical_input, prepared, trial, dt_heat, gate,
                        )?;
                        // Capture independent sources before adding signed reversible exchange.
                        let mut independent_entropy = 0.0;
                        let mut has_independent = false;
                        for (&id, cell) in trial.cells() {
                            for power in [cell.boundary.source_w, network_source.sources_w.get(&id).copied().unwrap_or(0.0)] {
                                if power != 0.0 { has_independent = true; }
                                independent_entropy = finite(independent_entropy+dt_heat*power/trial.temperature_k(id)?)?;
                            }
                        }
                        let mut next_states = Vec::with_capacity(members.len());
                        let mut exchanges = Vec::with_capacity(works.len());
                        let mut flat = 0;
                        for (index, (member, old)) in members.iter().zip(states).enumerate() {
                            check(gate)?;
                            let endpoint = mechanical.staged_components()[member.component].states();
                            let state = resolve_member(member, &temperatures(&member.region, trial)?, endpoint, dt)?;
                            for (layer, (a, b)) in old.layers.iter().zip(&state.layers).enumerate() {
                                let parameter = finite(b.free_energy_j-frozen[index].layers[layer].free_energy_j)?;
                                let entropy_energy = finite(b.temperature_k*b.entropy_j_k-a.temperature_k*a.entropy_j_k)?;
                                // Derive reversible heat from the source potential, keeping
                                // native spring-work roundoff OUT of physical heat.
                                let reversible = finite(-parameter-entropy_energy)?;
                                if finite(works[flat]-(b.internal_energy_j-a.internal_energy_j)-reversible)?.abs() > controls.total_energy_tolerance_j {
                                    return Err(invalid("thermoelastic spring work and source internal energy disagree"));
                                }
                                let power = finite(reversible/dt_heat)?;
                                if reversible != 0.0 && power == 0.0 { return Err(invalid("reversible frame heat power underflow")); }
                                let source = network_source.sources_w.entry(b.destination.cell).or_insert(0.0);
                                *source = finite(*source+power)?;
                                exchanges.push(FrameThermoelasticLayerExchange { states: [*a, *b],
                                    stress_work_j: works[flat], material_heat_j: losses[flat],
                                    coefficient_free_energy_change_j: parameter,
                                    mechanical_entropy_energy_change_j: entropy_energy,
                                    reversible_heat_j: reversible });
                                flat += 1;
                            }
                            next_states.push(state);
                        }
                        Ok(DrivenCoefficients {
                            circuit: electrical::Driver::new(MechanicalCandidate { frame: frame.clone(), stored_before_j: before,
                                material_heat_j, outside_loss_j: outside }, circuit)?,
                            sources_w: network_source.sources_w, conductances_w_k: network_source.conductances_w_k,
                            source: (network_source.source, next_states, exchanges, independent_entropy, has_independent),
                            circuit_trials: 1,
                        })
                    },
                    |driver, dt, checkpoint| {
                        checkpoint()?;
                        let mut combined = powers.clone();
                        match (&driver.electrical, electrical) {
                            (Some(candidate), Some(binding)) => binding.deposit(candidate, dt, &mut combined, checkpoint)?,
                            (None, None) => {},
                            _ => return Err(invalid("electrical heat and its frame-network binding disagree")),
                        }
                        Ok(combined)
                    },
                    true, || gate.is_requested(),
                )?;
                if next_circuit.is_some() != electrical.is_some() {
                    return Err(invalid("complete frame electrical candidate is missing"));
                }
                let (mechanical_report, electrical_report) = audit.electrical;
                let electrical_storage_change = electrical_report.as_ref().map_or(0.0,
                    |step| step.stored_after_j-step.stored_before_j);
                let electrical_work = electrical_report.as_ref().map_or(0.0, |step| step.native.supplied);
                let electrical_heat = electrical_report.as_ref().map_or(0.0, |step| step.heat_j);
                let (source, next_states, layers, independent_entropy, has_independent) = audit.source;
                let mut updates = Vec::with_capacity(transfers.len());
                for (member, state) in members.iter().zip(&next_states) {
                    let l = member.region.length_m;
                    let a = finite(state.viscosity.axial().viscous_rigidity_n_s()/l)?;
                    let b = finite(state.viscosity.bending().viscous_stiffness_n_m2_s()/l/l/l)?;
                    if (state.viscosity.axial().viscous_rigidity_n_s() > 0.0 && a == 0.0)
                        || (state.viscosity.bending().viscous_stiffness_n_m2_s() > 0.0 && b == 0.0) {
                        return Err(invalid("refreshed frame viscosity underflows its original geometry"));
                    }
                    for (k, damping_n_s_m) in [a, finite(3.0*b)?, b].into_iter().enumerate() {
                        updates.push(ConnectionCoefficients { connection: member.damping_start+k,
                            stiffness_n_m: 0.0, damping_n_s_m, rest_extension_m: 0.0 });
                    }
                }
                for (member, state) in members.iter().zip(&next_states) {
                    for layer in 0..state.layers.len() {
                        for (k, (stiffness_n_m, rest_extension_m)) in layer_springs(state, layer, member.region.length_m)?.into_iter().enumerate() {
                            updates.push(ConnectionCoefficients { connection: member.spring_start+3*layer+k,
                                stiffness_n_m, damping_n_s_m: 0.0, rest_extension_m });
                        }
                    }
                }
                // Frame/member order is physical reporting order, not native
                // connection order. The original coefficient owner requires a
                // strictly increasing batch, including noncontiguous frames.
                check(gate)?;
                updates.sort_unstable_by_key(|update| update.connection);
                let update = mechanical.prepare_staged_connection_coefficients(&updates, gate)?;
                let mut flat = 0;
                for (index, member) in members.iter().enumerate() {
                    for layer in 0..member.region.receivers.len() {
                        let first = member.spring_start+3*layer;
                        // Ordered native receipts locate one layer in O(log P)
                        // rather than repeatedly scanning every assembly port.
                        let all = update.parameter_changes();
                        let entries = if all.is_empty() { all } else {
                            let start = all.partition_point(|p| p.connection < first);
                            let entries = all.get(start..start+3)
                                .ok_or_else(|| invalid("native layer coefficient receipts are incomplete"))?;
                            if entries.iter().enumerate().any(|(k, p)| p.connection != first+k) {
                                return Err(invalid("native layer coefficient receipt addresses changed"));
                            }
                            entries
                        };
                        let native = sum(entries.iter().map(|p| p.parameter_work_j))?;
                        if finite(native-layers[flat].coefficient_free_energy_change_j)?.abs()
                            > controls.total_energy_tolerance_j {
                            return Err(invalid("native spring parameter work reached the wrong material layer"));
                        }
                        if !update.parameter_changes().is_empty() {
                            let before = sum(entries.iter().map(|p| p.free_energy_before_j))?;
                            let after = sum(entries.iter().map(|p| p.free_energy_after_j))?;
                            if finite(before-frozen[index].layers[layer].free_energy_j)?.abs()
                                > controls.total_energy_tolerance_j
                                || finite(after-next_states[index].layers[layer].free_energy_j)?.abs()
                                    > controls.total_energy_tolerance_j {
                                return Err(invalid("native layer spring storage differs from its source potential"));
                            }
                        }
                        flat += 1;
                    }
                }
                let source_parameter_work = sum(layers.iter().map(|l| l.coefficient_free_energy_change_j))?;
                let parameter_work_j = update.parameter_work_j();
                if finite(parameter_work_j-source_parameter_work)?.abs() > controls.total_energy_tolerance_j {
                    return Err(invalid("native endpoint stiffness work and original thermoelastic potential disagree"));
                }
                let reversible_heat_j = sum(layers.iter().map(|l| l.reversible_heat_j))?;
                let external_heat = finite(audit.non_electrical_heat_j-reversible_heat_j)?;
                let internal_after = finite(frame.total_energy_j()+parameter_work_j
                    +sum(next_states.iter().flat_map(|s| s.layers.iter().map(|l| l.temperature_k*l.entropy_j_k)))?)?;
                let residual = finite(internal_after-internal_before+audit.thermal.stored_energy_change_j
                    +outside-frame.external_work_j-external_heat
                    +electrical_storage_change-electrical_work)?;
                if residual.abs() > controls.total_energy_tolerance_j {
                    return Err(invalid("frame kinetics, electrical storage, material internal energy and finite caloric heat do not close"));
                }
                let e = entropy.ok_or_else(|| invalid("thermoelastic frame requires original entropy audit"))?;
                let mechanical_entropy = sum(layers.iter().map(|l| l.states[1].entropy_j_k-l.states[0].entropy_j_k))?;
                let reversible_entropy = sum(layers.iter().map(|l| l.reversible_heat_j/l.states[1].temperature_k))?;
                let split = finite(mechanical_entropy+reversible_entropy)?;
                let entropy = FrameThermoelasticEntropy {
                    thermal_storage_change_j_k: e.thermal_storage_change_j_k,
                    mechanical_storage_change_j_k: mechanical_entropy,
                    reservoir_entropy_input_j_k: e.reservoir_entropy_input_j_k,
                    conduction_generation_j_k: e.conduction_generation_j_k,
                    material_loss_generation_j_k: e.electrical_loss_generation_j_k,
                    unclassified_source_at_cell_j_k: independent_entropy, has_unclassified_sources: has_independent,
                    reversible_transfer_entropy_j_k: reversible_entropy,
                    thermal_temporal_defect_j_k: e.temporal_defect_j_k,
                    thermoelastic_temporal_defect_j_k: split,
                    algebraic_residual_j_k: e.algebraic_residual_j_k,
                    accounting_residual_j_k: sum([e.thermal_storage_change_j_k, mechanical_entropy,
                        -e.reservoir_entropy_input_j_k, -e.conduction_generation_j_k,
                        -e.electrical_loss_generation_j_k, -independent_entropy, -e.temporal_defect_j_k,
                        -split, -e.algebraic_residual_j_k])?,
                };
                let report = FrameThermoelasticReport {
                    mechanical: mechanical_report, electrical: electrical_report,
                    thermal: audit.thermal, source, members: next_states.clone(),
                    layers, connection_heat: loss_account.receivers, entropy, parameter_work_j, reversible_heat_j,
                    deposited_material_heat_j: material_heat_j, deposited_joule_heat_j: electrical_heat,
                    external_thermal_heat_j: external_heat, external_mechanical_loss_j: outside,
                    mechanical_internal_energy_before_j: internal_before,
                    mechanical_internal_energy_after_j: internal_after, total_energy_residual_j: residual,
                    nonlinear: audit.nonlinear, material_trials: audit.circuit_trials,
                    maximum_cell_residual_j: audit.maximum_cell_residual_j,
                };
                Ok((candidate, report, update, next_states, next_circuit))
            },
        )?;
        // Native admission checked cancellation immediately before publication.
        self.mechanical.install_connection_coefficients(update);
        self.states = next_states;
        if let (Some(binding), Some(circuit)) = (&mut self.electrical, next_circuit) {
            binding.circuit = circuit;
        }
        *thermal = next_thermal;
        Ok(report)
    }
}

/// Complete exterior-fluid memory on this same thermoelastic/electrical timeline.
pub mod fluid;
