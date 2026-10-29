//! Source ownership through the original complete mechanical assembly.
use super::*;
use core::ops::Range;

/// Immutable placement of one source frame in the complete coupled assembly.
/// Member ranges index the accepted system states and its step reports.
pub struct FrameThermoelasticFrame {
    pub(super) basis: FrameVibrationBasis,
    pub(super) component: usize,
    pub(super) connections: Range<usize>,
    members: Range<usize>,
}
impl FrameThermoelasticFrame {
    pub const fn basis(&self) -> &FrameVibrationBasis { &self.basis }
    pub const fn component_index(&self) -> usize { self.component }
    pub fn connection_range(&self) -> Range<usize> { self.connections.clone() }
    pub fn member_range(&self) -> Range<usize> { self.members.clone() }
}

/// Source/caloric ownership for the returned free-coordinate frame and ports.
/// Keep every definition unchanged when assembling the original native owner.
/// Single and multiple frame bindings use the same admission and step owner.
/// A shared assembly requires one original MaterialNetwork and disjoint
/// component, material-port and whole-layer caloric inventories.
pub struct FrameThermoelasticBinding<'n, 'a> {
    pub(super) basis: FrameVibrationBasis,
    pub(super) prepared: PreparedMaterialNetwork<'n, 'a>,
    pub(super) source: &'n MaterialNetwork<'a>,
    pub(super) members: Vec<Member<'a>>,
    pub(super) states: Vec<FrameThermoelasticMemberState>,
    pub(super) origin_s: f64,
    pub(super) component: usize,
    pub(super) connection_start: usize,
    pub(super) initial_model: ModalAcousticTimeModel,
    pub(super) definitions: Vec<ModalConnection>,
}

impl<'n, 'a> FrameThermoelasticBinding<'n, 'a> {
    /// Original fixed mass basis for physical forcing and observations.
    pub const fn basis(&self) -> &FrameVibrationBasis {
        &self.basis
    }

    /// Admit the source frame into the SAME fully configured native owner.
    ///
    /// The network must be cold, with the unchanged explicit initial frame
    /// q/v, modal limits, free coordinates and every original material port.
    /// Signed columns are checked over ALL components, including foreign zeros.
    /// A source dashpot cannot be reinterpreted as a constant-loss fit. Native
    /// preparation, other components, initialized Maxwell memory, acoustic
    /// state, port identities and all complete-system budgets remain untouched.
    ///
    /// The reciprocal acoustic producer still owns its geometry/gas/basis and
    /// fitted-band validity. Native passivity does not establish that fidelity.
    /// This is not a warm rebind: unobserved prior frame work/heat cannot be
    /// inferred from endpoint motion. Subsequent accepted history stays on the
    /// same mechanical/finite-thermal transaction, including refused retries.
    pub fn into_system(
        self,
        mechanical: CoupledModalSystem,
        thermal: &DepositionSimulation<'_>,
        gate: &CancelGate,
    ) -> Result<FrameThermoelasticSystem<'n, 'a>> {
        Self::into_shared_system(vec![self], mechanical, thermal, gate)
    }

    /// Bind several independently heated source frames to one original native
    /// assembly and one finite thermal timeline. Mechanical components and
    /// protected port blocks may appear in any order, with ordinary connectors
    /// between them. All bindings must borrow the SAME complete MaterialNetwork.
    /// Cross-frame conduction is then resolved once, at the same thermal endpoint.
    ///
    /// Every original frame state, material port, source, whole-cell inventory
    /// and clock is admitted before ownership is consumed. Duplicate components,
    /// overlapping port ranges and reused caloric cells refuse. No coefficients
    /// or coordinates are changed here. Binding order cannot change reductions:
    /// accepted frames and member reports use native component order.
    ///
    /// Evolution retains one native step/admit/publish boundary. No frame can
    /// advance ahead of another, and one failed source or heat solve rolls back
    /// the complete assembly. Fixed mass/geometry and first-order splitting
    /// limits are the same as for a single frame.
    pub fn into_shared_system(
        mut bindings: Vec<Self>,
        mechanical: CoupledModalSystem,
        thermal: &DepositionSimulation<'_>,
        gate: &CancelGate,
    ) -> Result<FrameThermoelasticSystem<'n, 'a>> {
        check(gate)?;
        if bindings.is_empty() || bindings.len() > mechanical.components().len() {
            return Err(invalid("nonempty bounded thermoelastic frame bindings required"));
        }
        bindings.sort_unstable_by_key(|binding| binding.component);
        let source = bindings[0].source;
        let origin_s = bindings[0].origin_s;
        let mut used_cells = BTreeSet::new();
        let mut owned_connections = Vec::with_capacity(bindings.len());
        let mut previous_component = None;
        for binding in &bindings {
            check(gate)?;
            if previous_component == Some(binding.component)
                || !core::ptr::eq(source, binding.source)
                || binding.origin_s.to_bits() != origin_s.to_bits()
            {
                return Err(invalid(
                    "shared thermoelastic frames need distinct components, one original material network and one clock",
                ));
            }
            binding.admit(&mechanical, thermal, gate)?;
            previous_component = Some(binding.component);
            let end = binding.connection_start.checked_add(binding.definitions.len())
                .ok_or_else(|| invalid("thermoelastic connection range overflows"))?;
            owned_connections.push(binding.connection_start..end);
            for receiver in binding.members.iter().flat_map(|member| &member.region.receivers) {
                check(gate)?;
                if !used_cells.insert(receiver.binding.cell) {
                    return Err(invalid("a whole caloric cell cannot belong to multiple frame layers"));
                }
            }
        }
        owned_connections.sort_unstable_by_key(|range| range.start);
        if owned_connections.windows(2).any(|pair| pair[0].end > pair[1].start) {
            return Err(invalid("thermoelastic material connection ranges overlap"));
        }
        let mut frames = Vec::with_capacity(bindings.len());
        let mut members = Vec::new();
        let mut states = Vec::new();
        let mut prepared = None;
        for mut binding in bindings {
            check(gate)?;
            let first = members.len();
            members.append(&mut binding.members);
            states.append(&mut binding.states);
            frames.push(FrameThermoelasticFrame {
                basis: binding.basis,
                component: binding.component,
                connections: binding.connection_start..binding.connection_start + binding.definitions.len(),
                members: first..members.len(),
            });
            // All preparations select from the same immutable source network
            // at the same admitted initial state. Retain one transport owner,
            // not one repeated full-network resolution per frame per trial.
            if prepared.is_none() {
                prepared = Some(binding.prepared);
            }
        }
        check(gate)?;
        Ok(FrameThermoelasticSystem {
            frames, owned_connections, mechanical,
            prepared: prepared.ok_or_else(|| invalid("missing shared material preparation"))?,
            members, states, origin_s,
        })
    }

    fn admit(
        &self,
        mechanical: &CoupledModalSystem,
        thermal: &DepositionSimulation<'_>,
        gate: &CancelGate,
    ) -> Result<()> {
        check(gate)?;
        let model = mechanical
            .components()
            .get(self.component)
            .ok_or_else(|| invalid("thermoelastic frame component is absent from assembly"))?;
        if mechanical.samples_rendered() != 0
            || thermal.time_s().to_bits() != self.origin_s.to_bits()
            || model.sample_period_s().to_bits() != self.initial_model.sample_period_s().to_bits()
            || model.modes() != self.initial_model.modes()
            || model.states() != self.initial_model.states()
            || model.budget() != self.initial_model.budget()
        {
            return Err(invalid(
                "thermoelastic assembly requires its cold original frame state, clock, coordinates and limits",
            ));
        }
        let offset: usize = mechanical.components()[..self.component]
            .iter()
            .map(|m| m.modes().len())
            .sum();
        let end = self.connection_start + self.definitions.len();
        if mechanical
            .constant_loss_connections()
            .iter()
            .any(|c| (self.connection_start..end).contains(&c.connection))
        {
            return Err(invalid(
                "source thermoelastic ports cannot become constant-loss fits",
            ));
        }
        for (i, expected) in self.definitions.iter().enumerate() {
            check(gate)?;
            let index = self.connection_start + i;
            let law = mechanical.connection_coefficients(index).ok_or_else(|| {
                invalid("thermoelastic source connection is absent from assembly")
            })?;
            let column = mechanical
                .connection_column(index)
                .ok_or_else(|| invalid("thermoelastic source column is absent from assembly"))?;
            if law.stiffness_n_m.to_bits() != expected.stiffness_n_m.to_bits()
                || law.damping_n_s_m.to_bits() != expected.damping_n_s_m.to_bits()
                || law.rest_extension_m.to_bits() != expected.rest_extension_m.to_bits()
                || column.len() != mechanical.mode_count()
                || column.iter().enumerate().any(|(j, &b)| {
                    let value = if (offset..offset + model.modes().len()).contains(&j) {
                        expected.left.shapes[j - offset] - expected.right.shapes[j - offset]
                    } else {
                        0.0
                    };
                    b != value
                })
            {
                return Err(invalid(
                    "assembled frame material columns or coefficients differ from their original sources",
                ));
            }
        }
        for (member, state) in self.members.iter().zip(&self.states) {
            check(gate)?;
            if temperatures(&member.region, thermal)?
                .iter()
                .zip(&state.layers)
                .any(|(&t, s)| t.to_bits() != s.temperature_k.to_bits())
            {
                return Err(invalid(
                    "thermoelastic caloric state changed before assembly binding",
                ));
            }
        }
        check(gate)?;
        Ok(())
    }
}

// A foreign spring/dashpot's integrated work is -delta(Psi)-D. Recover its
// OWN loss from that original ledger and original endpoint coordinates. Never
// call conservative support work heat, and never route a foreign damper into
// this frame's caloric inventory. The native independent total must agree.
pub(super) fn outside_loss(
    mechanical: &CoupledModalSystem,
    frame: &CoupledModalFrame,
    transfers: &[fs_ports::IntegratedWindowTransfer],
    owned_connections: &[Range<usize>],
    material_heat_j: f64,
    tolerance: f64,
    gate: &CancelGate,
) -> Result<f64> {
    let mut owned = 0;
    let mut foreign = 0.0;
    for (i, transfer) in transfers.iter().enumerate() {
        check(gate)?;
        // Both lists are canonical; this is linear in connections + frames,
        // including unbound connectors between noncontiguous source blocks.
        while owned < owned_connections.len() && owned_connections[owned].end <= i {
            owned += 1;
        }
        if owned < owned_connections.len() && owned_connections[owned].contains(&i) {
            continue;
        }
        let law = mechanical
            .connection_coefficients(i)
            .ok_or_else(|| invalid("foreign connection work has no original law"))?;
        if law.damping_n_s_m == 0.0 {
            continue;
        }
        let column = mechanical
            .connection_column(i)
            .ok_or_else(|| invalid("foreign connection work has no original column"))?;
        let extension = |models: &[ModalAcousticTimeModel]| -> Result<f64> {
            finite(
                sum(column
                    .iter()
                    .zip(models.iter().flat_map(|m| m.states()))
                    .map(|(&b, s)| b * s.displacement_m_sqrt_kg))?
                    - law.rest_extension_m,
            )
        };
        let before = extension(mechanical.components())?;
        let after = extension(mechanical.staged_components())?;
        let stored_change = finite(0.5 * law.stiffness_n_m * (after - before) * (after + before))?;
        let loss = finite(-transfer.energy_into().value() - stored_change)?;
        if loss < -tolerance {
            return Err(invalid(
                "foreign connection work implies negative dissipation",
            ));
        }
        foreign = finite(foreign + loss.max(0.0))?;
    }
    if finite(foreign + material_heat_j - frame.connection_dissipation_j)?.abs() > tolerance {
        return Err(invalid(
            "frame and foreign connector losses disagree with the complete native account",
        ));
    }
    sum([
        foreign,
        frame.component_dissipation_j,
        frame.material_dissipation_j,
        frame.radiation_dissipation_j,
    ])
}
