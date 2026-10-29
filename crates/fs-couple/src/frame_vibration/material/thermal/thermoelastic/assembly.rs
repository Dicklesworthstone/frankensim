//! Source ownership through the original complete mechanical assembly.
use super::*;

/// Source/caloric ownership for the returned free-coordinate frame and ports.
/// Keep every definition unchanged when assembling the original native owner.
/// This binds ONE live thermoelastic frame; companions retain their own laws
/// and their losses are not silently deposited into this frame's finite cells.
pub struct FrameThermoelasticBinding<'n, 'a> {
    pub(super) basis: FrameVibrationBasis,
    pub(super) prepared: PreparedMaterialNetwork<'n, 'a>,
    pub(super) members: Vec<Member<'a>>,
    pub(super) states: Vec<FrameThermoelasticMemberState>,
    pub(super) origin_s: f64,
    pub(super) linear_strain_limit: f64,
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
        Ok(FrameThermoelasticSystem {
            basis: self.basis,
            mechanical,
            prepared: self.prepared,
            members: self.members,
            states: self.states,
            origin_s: self.origin_s,
            linear_strain_limit: self.linear_strain_limit,
            component: self.component,
            connection_start: self.connection_start,
        })
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
    connection_start: usize,
    members: &[Member<'_>],
    material_heat_j: f64,
    tolerance: f64,
    gate: &CancelGate,
) -> Result<f64> {
    let end = members.last().map_or(connection_start, |m| {
        m.spring_start + 3 * m.region.receivers.len()
    });
    let mut foreign = 0.0;
    for (i, transfer) in transfers.iter().enumerate() {
        check(gate)?;
        if (connection_start..end).contains(&i) {
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
