//! Source-derived bilateral bridge laws in the existing reciprocal mechanics.
//! The original physical body footprint determines the receiving material;
//! modal coordinates, rigid tie representatives and names do not select it.

use super::plate::SharedPlateBridge;
use super::{
    MaterialStringBankSource, MaterialStringMember, MaterialStringSource, ReciprocalStringBank,
    ReciprocalStringError, interior, invalid, poll,
};
use crate::render::schedule::force::coupled::ConstantLossConnection;
use crate::render::schedule::force::coupled::contact::friction::FrictionModalSystem;
use fs_exec::CancelGate;
use fs_matdb::MaterialStateId;
use fs_material::state_point::{ResolvedLinearInterface, ResolvedMaterialStatePoint};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{QuantitySpec, Temperature};
use std::collections::BTreeMap;

/// The original ordered surface law used by one existing physical bridge.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialBridgeInterface {
    /// Original string declaration index; its material body has one bridge.
    pub string: usize,
    /// Complete source, selected properties, physical area and integrated k/c.
    pub interface: ResolvedLinearInterface,
}

impl SharedPlateBridge {
    /// Lower a sourced scalar law into the existing work-conjugate connection.
    /// The span station and stress-free extension remain explicit geometry.
    /// Bind the same law with `with_material_bridge_interfaces` after the body
    /// constructor has established its actual receiving footprint and modes.
    pub fn from_interface(
        station_fraction: f64,
        rest_extension_m: f64,
        law: &ResolvedLinearInterface,
    ) -> Result<Self, ReciprocalStringError> {
        if !interior(station_fraction) || !rest_extension_m.is_finite() {
            return Err(invalid(
                "material bridge requires an interior station and finite rest extension",
            ));
        }
        Ok(Self {
            station_fraction,
            rest_extension_m,
            stiffness_n_m: law.stiffness_n_m(),
            damping_n_s_m: law.damping_n_s_m(),
        })
    }
}

impl MaterialStringBankSource {
    /// Original bindings, retained through gas/torque preparation, nonlinear
    /// storage, multiple bows and scheduled rendering without replacing sources.
    #[must_use]
    pub fn material_bridge_interfaces(&self) -> &[MaterialBridgeInterface] {
        &self.material_bridge_interfaces
    }

    /// Actual outer-string and receiving-body material states, in that order.
    /// All nonzero original physical footprint nodes, including shell moments,
    /// must meet one identical complete resolved body state. A footprint across
    /// different states refuses; no expected interface material selects a side.
    /// Raw modal bodies have no source geometry and cannot answer this query.
    pub fn bridge_surface_materials(
        &self,
        string: usize,
    ) -> Result<[MaterialStateId; 2], ReciprocalStringError> {
        let [first, second] = self.bridge_surface_states(string, None)?;
        Ok([first.material().clone(), second.material().clone()])
    }

    fn bridge_surface_states(
        &self,
        string: usize,
        gate: Option<&CancelGate>,
    ) -> Result<[&ResolvedMaterialStatePoint; 2], ReciprocalStringError> {
        let member = self
            .strings
            .get(string)
            .ok_or_else(|| invalid("material bridge names an unknown string"))?;
        let first = outer_material(&member.source)?;
        let second =
            if let Some(body) = &self.material_plate {
                body.specimen.material()
            } else {
                let mut selected = None;
                if let Some(body) = &self.material_chart {
                    let mesh = &body.specimen.chart().mesh;
                    let row =
                        match &self.material_chart_bridges {
                            Some(ports) => ports
                                .weights
                                .get(string)
                                .map(Vec::as_slice)
                                .ok_or_else(|| {
                                    invalid("missing original regional plate bridge footprint")
                                })?,
                            None => &body.reduction.unit_force_weights,
                        };
                    footprint(
                        mesh.nodes.len(),
                        &mesh.tris,
                        row,
                        1,
                        |e| body.specimen.material_at_element(e),
                        &mut selected,
                        gate,
                    )?;
                } else if let Some(body) = &self.material_shell {
                    let mesh = body.specimen.mesh();
                    let row = body
                        .force_ports
                        .get(string)
                        .ok_or_else(|| invalid("missing original shell bridge footprint"))?;
                    footprint(
                        mesh.nodes.len(),
                        &mesh.tris,
                        row,
                        6,
                        |e| body.specimen.material_at_element(e),
                        &mut selected,
                        gate,
                    )?;
                } else if let Some(body) = &self.material_shell_assembly {
                    let port = body.force_ports.get(string).ok_or_else(|| {
                        invalid("missing original assembled body bridge footprint")
                    })?;
                    if port.components.len() != body.specimen.components.len() {
                        return Err(invalid(
                            "assembled bridge needs every original component row",
                        ));
                    }
                    for (part, row) in body.specimen.components.iter().zip(&port.components) {
                        let mesh = part.specimen.mesh();
                        footprint(
                            mesh.nodes.len(),
                            &mesh.tris,
                            row,
                            6,
                            |e| part.specimen.material_at_element(e),
                            &mut selected,
                            gate,
                        )?;
                    }
                }
                selected.ok_or_else(|| {
                    invalid("material bridge requires an actual sourced body footprint")
                })?
            };
        Ok([first, second])
    }

    fn original_bridge(&self, string: usize) -> Result<&SharedPlateBridge, ReciprocalStringError> {
        self.material_plate
            .as_ref()
            .and_then(|b| b.bridges.get(string))
            .or_else(|| {
                self.material_chart
                    .as_ref()
                    .and_then(|b| b.bridges.get(string))
            })
            .or_else(|| {
                self.material_shell
                    .as_ref()
                    .and_then(|b| b.bridges.get(string))
            })
            .or_else(|| {
                self.material_shell_assembly
                    .as_ref()
                    .and_then(|b| b.bridges.get(string))
            })
            .ok_or_else(|| {
                invalid("material interface needs the original sole material-body bridge")
            })
    }

    /// Every bound bridge whose interface states a constant-loss anchor, named
    /// as its connection in the assembled network.
    pub(super) fn constant_loss_bridges(&self) -> Vec<ConstantLossConnection> {
        bridge_anchors(&self.strings, &self.material_bridge_interfaces)
    }

    pub(super) fn validate_material_bridge_frequencies(
        &self,
        system: &FrictionModalSystem,
        gate: &CancelGate,
    ) -> Result<(), ReciprocalStringError> {
        validate_frequencies(
            self.strings.len() + 1 + usize::from(self.material_cavity.is_some())
                + usize::from(self.spatial_cavity.as_ref().is_some_and(|cavity|
                    cavity.coupling().total_modes() > cavity.coupling().structural_modes())),
            &self.material_bridge_interfaces,
            system,
            gate,
        )
    }
}

impl ReciprocalStringBank {
    /// Retain source laws only after validating the already constructed bridges.
    /// Coefficients must match exactly; this never edits a connection, rebuilds
    /// modes, changes vibration or advances the clock. Every supplied binding
    /// uses the actual outer string and physical body footprint, never a name
    /// inferred from the interface. Supply all desired bindings in one call.
    ///
    /// A scalar law requires one resolved state over the entire footprint.
    /// Mixed-material seats need explicitly separated physical interfaces.
    /// Where frozen absolute temperatures are supplied on both sides, they must
    /// equal the interface temperature. All retained physical string/body/cavity modes, plus
    /// the existing owner's conservative COMPLETE linear stiffness upper bound,
    /// must lie inside both original coefficient domains. This is neither a
    /// nonlinear spectrum qualification nor continuum convergence certification.
    /// A physical cavity neck whose restoring spring is a body connection gives
    /// a conservative lower bound of zero; appended viscous gas memory does not
    /// masquerade as an independent physical vibration mode.
    ///
    /// The one exception to "never edits": an interface that states a
    /// constant-loss anchor (`ResolvedLinearInterface::constant_loss_reference_omega_rad_s`)
    /// binds its bridge's dashpot at a constant loss factor, which prepares the
    /// SAME network again before playback (the prepared owner is required).
    /// Later fluid, cavity and rotational rebuilds keep the anchor.
    pub fn with_material_bridge_interfaces(
        mut self,
        bindings: Vec<MaterialBridgeInterface>,
        gate: &CancelGate,
    ) -> Result<Self, ReciprocalStringError> {
        poll(gate)?;
        if bindings.is_empty()
            || bindings.len() > self.source.strings.len()
            || !self.source.material_bridge_interfaces.is_empty()
        {
            return Err(invalid(
                "material bridges require one nonempty bounded set of source bindings",
            ));
        }
        if self.source.residual_flexibility {
            // Its links are already softened; an interface law would re-anchor them.
            return Err(invalid("material bridge interfaces are not wired with residual flexibility"));
        }
        let mut seen = vec![false; self.source.strings.len()];
        for binding in &bindings {
            poll(gate)?;
            let member = self
                .source
                .strings
                .get(binding.string)
                .ok_or_else(|| invalid("material bridge names an unknown string"))?;
            if seen[binding.string] || member.connections.len() != 1 {
                return Err(invalid(
                    "material bridge bindings must be unique sole string/body connections",
                ));
            }
            seen[binding.string] = true;
            let bridge = self.source.original_bridge(binding.string)?;
            let connection = &member.connections[0];
            let law = &binding.interface;
            if connection.body.component != 0
                || connection.station_fraction.to_bits() != bridge.station_fraction.to_bits()
                || connection.rest_extension_m.to_bits() != bridge.rest_extension_m.to_bits()
                || [connection.stiffness_n_m, bridge.stiffness_n_m]
                    .iter()
                    .any(|v| v.to_bits() != law.stiffness_n_m().to_bits())
                || [connection.damping_n_s_m, bridge.damping_n_s_m]
                    .iter()
                    .any(|v| v.to_bits() != law.damping_n_s_m().to_bits())
            {
                return Err(invalid(
                    "material interface coefficients differ from the actual original bridge",
                ));
            }
            let [first, second] = self
                .source
                .bridge_surface_states(binding.string, Some(gate))?;
            law.require_surface_materials(first.material(), second.material())
                .map_err(|_| invalid("material bridge interface does not match its actual ordered surface states"))?;
            let source = law.resolved();
            if let Some(t) = temperature(source.query_point(), source.axis_quantities())? {
                for material in [first, second] {
                    if temperature(material.query_point(), material.axis_quantities())?
                        .is_some_and(|actual| actual.to_bits() != t.to_bits())
                    {
                        return Err(invalid(
                            "material bridge interface and receiving surfaces have different frozen temperatures",
                        ));
                    }
                }
            }
        }
        validate_frequencies(
            self.source.strings.len() + 1 + usize::from(self.source.material_cavity.is_some())
                + usize::from(self.source.spatial_cavity.as_ref().is_some_and(|cavity|
                    cavity.coupling().total_modes() > cavity.coupling().structural_modes())),
            &bindings,
            &self.system,
            gate,
        )?;
        let anchors = bridge_anchors(&self.source.strings, &bindings);
        if !anchors.is_empty() {
            self.system = self.system.with_constant_loss_connections(&anchors, gate)
                .map_err(ReciprocalStringError::Friction)?;
        }
        poll(gate)?;
        self.source.material_bridge_interfaces = bindings;
        Ok(self)
    }
}

/// A bridge's connection index: every string's own connections come first in
/// the assembled network, in string order. Bindings name existing strings.
fn bridge_anchors(
    strings: &[MaterialStringMember],
    bindings: &[MaterialBridgeInterface],
) -> Vec<ConstantLossConnection> {
    bindings
        .iter()
        .filter_map(|binding| {
            let reference_omega_rad_s = binding.interface.constant_loss_reference_omega_rad_s()?;
            let connection = strings[..binding.string]
                .iter()
                .map(|member| member.connections.len())
                .sum();
            Some(ConstantLossConnection { connection, reference_omega_rad_s })
        })
        .collect()
}

fn outer_material(
    source: &MaterialStringSource,
) -> Result<&ResolvedMaterialStatePoint, ReciprocalStringError> {
    match source {
        MaterialStringSource::WithPolarizations(p) => outer_material(p.source()),
        MaterialStringSource::WithTorsion(t) => outer_material(t.bending()),
        MaterialStringSource::Homogeneous(s) => Ok(s.material()),
        MaterialStringSource::Wound(s) => Ok(s.section().winding_material()),
        MaterialStringSource::Concentric { section, .. } => section
            .elastic()
            .layers()
            .last()
            .map(|layer| &layer.material)
            .ok_or_else(|| invalid("material string has no outer surface")),
    }
}

#[allow(clippy::too_many_arguments)]
fn footprint<'a>(
    nodes: usize,
    triangles: &[[usize; 3]],
    row: &[f64],
    dofs: usize,
    mut material: impl FnMut(usize) -> Option<&'a ResolvedMaterialStatePoint>,
    selected: &mut Option<&'a ResolvedMaterialStatePoint>,
    gate: Option<&CancelGate>,
) -> Result<(), ReciprocalStringError> {
    if nodes.checked_mul(dofs) != Some(row.len()) || row.iter().any(|v| !v.is_finite()) {
        return Err(invalid(
            "material bridge requires its complete finite original physical row",
        ));
    }
    let active: Vec<_> = row
        .chunks_exact(dofs)
        .map(|d| d.iter().any(|v| *v != 0.0))
        .collect();
    let mut incident = vec![false; nodes];
    for (e, triangle) in triangles.iter().enumerate() {
        if let Some(gate) = gate {
            poll(gate)?;
        }
        if triangle.iter().any(|&node| node >= nodes) {
            return Err(invalid(
                "material bridge footprint has an invalid incident triangle",
            ));
        }
        if triangle.iter().any(|&node| active[node]) {
            let state =
                material(e).ok_or_else(|| invalid("bridge triangle has no material source"))?;
            if selected.is_some_and(|previous| previous.identity() != state.identity()) {
                return Err(invalid(
                    "scalar material bridge footprint crosses different resolved material states",
                ));
            }
            *selected = Some(state);
            for &node in triangle {
                incident[node] |= active[node];
            }
        }
    }
    if active
        .iter()
        .zip(incident)
        .any(|(&active, incident)| active && !incident)
    {
        return Err(invalid(
            "material bridge acts on a node without an incident material surface",
        ));
    }
    Ok(())
}

fn validate_frequencies(
    physical_components: usize,
    bindings: &[MaterialBridgeInterface],
    system: &FrictionModalSystem,
    gate: &CancelGate,
) -> Result<(), ReciprocalStringError> {
    if bindings.is_empty() {
        return Ok(());
    }
    let physical = system
        .components()
        .get(..physical_components)
        .ok_or_else(|| invalid("material bridge lost its original physical component prefix"))?;
    for binding in bindings {
        poll(gate)?;
        for component in physical {
            for mode in component.modes() {
                // Physical zero-stiffness coordinates, including cavity neck
                // inertia whose spring is a body link, make the conservative
                // lower bound zero. Appended cylinder memory is outside this
                // physical prefix, not an independent vibration frequency.
                binding.interface.require_angular_frequency_rad_s(mode.angular_frequency_rad_s)
                    .map_err(|_| invalid("retained string/body mode is outside material bridge frequency applicability"))?;
            }
        }
        binding.interface.require_angular_frequency_rad_s(system.linear_frequency_upper_bound_rad_s())
            .map_err(|_| invalid("coupled linear stiffness bound is outside material bridge frequency applicability"))?;
    }
    Ok(())
}

fn temperature(
    point: &[(String, f64)],
    quantities: &BTreeMap<String, QuantitySpec>,
) -> Result<Option<f64>, ReciprocalStringError> {
    let absolute = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    let mut found = None;
    for (axis, value) in point {
        let quantity = quantities.get(axis).copied();
        let is_temperature = quantity == Some(absolute)
            || axis == "T"
                && quantity.is_none_or(|q| q == QuantitySpec::dimensional(Temperature::DIMS));
        if is_temperature {
            if !value.is_finite() || *value <= 0.0 || found.replace(*value).is_some() {
                return Err(invalid(
                    "material bridge requires an unambiguous positive absolute temperature",
                ));
            }
        } else if axis == "T" {
            return Err(invalid(
                "material bridge temperature has incompatible quantity semantics",
            ));
        }
    }
    Ok(found)
}

#[cfg(test)]
mod tests;
