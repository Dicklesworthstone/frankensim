//! Source-bound unsteady gas loading in the SAME string/body/friction network.
//!
//! Every transverse mode gets the cylinder continuum's potential-flow inertia
//! and positive viscous memory branches. The branches are real free coordinates
//! coupled through dashpots in the existing exact linear owner; their kinetic
//! energy, reactions and dissipation enter the original work ledger. Nothing is
//! applied to the acoustic samples. Both physical polarizations share one gas
//! and cross-section. Selected twist modes also carry the rotational Stokes
//! continuum, using the original solid polar inertia and the actual radius.
//!
//! Local two-dimensional flow is a slender-cylinder approximation along a
//! vibrating string. Spectral cutoffs, quadrature and string modes must be
//! refined. No nonlinear wake, confinement or radiation claim.

use super::*;
use crate::air_path::cylinder::{CylinderFlow, CylinderFlowResolution};
use crate::air_path::cylinder::rotary::RotaryCylinderFlow;
use fs_math::det;

/// One original solid realization and its separately retained fluid continuum.
#[derive(Clone, Debug)]
pub struct StringCylinderLoading {
    solid_card: BowedStringCard,
    flow: CylinderFlow,
    rotary: Option<RotaryCylinderFlow>,
}
impl StringCylinderLoading {
    /// Original solid mass, bending, tension and old sampled loss realization.
    /// The old card is provenance only; its sampled gas resistance is replaced.
    #[must_use]
    pub const fn solid_card(&self) -> &BowedStringCard { &self.solid_card }
    /// Source-derived added inertia, memory masses, resistance and resolution.
    #[must_use]
    pub const fn flow(&self) -> &CylinderFlow { &self.flow }
    /// Rotational fluid memory for the same section, when solid torsion exists.
    #[must_use]
    pub const fn rotary_flow(&self) -> Option<&RotaryCylinderFlow> { self.rotary.as_ref() }
}

/// Frozen fluid and structural provenance retained through every performance.
#[derive(Clone, Debug)]
pub struct StringBankFluidLoading {
    resolution: CylinderFlowResolution,
    strings: Vec<StringCylinderLoading>,
    pub(super) body_components: usize,
}
impl StringBankFluidLoading {
    #[must_use]
    pub const fn resolution(&self) -> CylinderFlowResolution { self.resolution }
    /// Same original material-string order. Fluid is never merged into a solid
    /// density property or hidden inside per-mode authored damping.
    #[must_use]
    pub fn strings(&self) -> &[StringCylinderLoading] { &self.strings }

    // The exact steady rotary resistance and collapsed low tail act against
    // stationary gas. Every remaining torque is a reciprocal fluid branch.
    // Recompute from the original solid modes on each assembly, never add to
    // an already loaded image or change the material's polar inertia.
    pub(super) fn add_rotary_damping(&self, index: usize, source: &MaterialStringSource,
        card: &BowedStringCard, modes: &mut [ModalAcousticMode])
        -> Result<(), ReciprocalStringError>
    {
        let loading = self.strings.get(index)
            .ok_or_else(|| invalid("rotary loading lost its original string"))?;
        let Some(flow) = &loading.rotary else {
            if source.torsion_source().is_some() {
                return Err(invalid("cylinder loading omitted the original torsion"));
            }
            return Ok(());
        };
        let twist = source.torsion_source()
            .ok_or_else(|| invalid("rotary fluid has no original torsional material"))?;
        let start = transverse_count(source, card);
        if modes.len() != start + twist.mode_count() {
            return Err(invalid("rotary fluid lost the original torsional basis"));
        }
        let rate = (flow.steady_resistance_per_length_n_s()
            + flow.low_tail_resistance_per_length_n_s()) / twist.section().rotary_inertia_kg_m();
        for mode in &mut modes[start..] {
            let added = 0.5 * rate / mode.angular_frequency_rad_s;
            let loss = mode.damping_ratio + added;
            if !added.is_finite() || added <= 0.0 || !loss.is_finite() {
                return Err(invalid("rotary gas resistance is not representable in the solid basis"));
            }
            mode.damping_ratio = loss;
        }
        Ok(())
    }

    pub(super) fn append_components(
        &self,
        cards: &[BowedStringCard],
        sources: &[MaterialStringMember],
        components: &mut Vec<ModalAcousticTimeModel>,
        links: &mut Vec<ModalConnection>,
        options: &StringBankOptions,
        gate: &CancelGate,
    ) -> Result<(), ReciprocalStringError> {
        if self.strings.len() != cards.len() || sources.len() != cards.len()
            || components.len() != cards.len() + self.body_components {
            return Err(invalid("cylinder memory must retain original strings and physical body layout"));
        }
        let fluid_start = components.len();
        for (index, (card, loading)) in cards.iter().zip(&self.strings).enumerate() {
            poll(gate)?;
            let transverse_modes = transverse_count(&sources[index].source, card);
            let fluid_modes = sources[index].source.modes().checked_mul(loading.flow.branches().len())
                .ok_or_else(|| invalid("cylinder memory dimension overflow"))?;
            let mut modes = Vec::with_capacity(fluid_modes);
            let solid_shape = det::sqrt(2.0 / card.linear_density_kg_m / card.length_m);
            for k in 0..transverse_modes {
                for (branch_index, branch) in loading.flow.branches().iter().enumerate() {
                    poll(gate)?;
                    let branch_shape = det::sqrt(2.0 / branch.mass_per_length_kg_m / card.length_m);
                    let resistance = branch.resistance_per_length_n_s_m2 * (0.5 * card.length_m);
                    if [solid_shape, branch_shape, resistance].iter()
                        .any(|v| !v.is_finite() || *v <= 0.0) {
                        return Err(invalid("fluid/string work-conjugate port is unrepresentable"));
                    }
                    let mut left = vec![0.0; sources[index].source.modes()];
                    let mut right = vec![0.0; fluid_modes];
                    left[k] = solid_shape;
                    right[k * loading.flow.branches().len() + branch_index] = branch_shape;
                    links.push(ModalConnection {
                        left: ModalAttachment { component: index, shapes: left },
                        right: ModalAttachment { component: fluid_start + index, shapes: right },
                        stiffness_n_m: 0.0,
                        damping_n_s_m: resistance,
                        rest_extension_m: 0.0,
                    });
                    modes.push(ModalAcousticMode {
                        angular_frequency_rad_s: 0.0,
                        damping_ratio: 0.0,
                        pressure_per_modal_velocity: fs_math::c64::C64::ZERO,
                    });
                }
            }
            if let Some(flow) = &loading.rotary {
                let twist = sources[index].source.torsion_source()
                    .ok_or_else(|| invalid("rotary fluid has no original torsion"))?;
                let radius = flow.radius_m();
                let solid_surface = radius * det::sqrt(2.0
                    / twist.section().rotary_inertia_kg_m() / card.length_m);
                let offset = transverse_modes * loading.flow.branches().len();
                for k in 0..twist.mode_count() {
                    for (branch_index, branch) in flow.branches().iter().enumerate() {
                        poll(gate)?;
                        // Convert angle/torque to surface displacement/force on
                        // BOTH sides of one dashpot; r cancels from its work.
                        let fluid_surface = radius * det::sqrt(2.0
                            / branch.inertia_per_length_kg_m / card.length_m);
                        let resistance = branch.resistance_per_length_n_s
                            * (0.5 * card.length_m) / radius / radius;
                        if [solid_surface, fluid_surface, resistance].iter()
                            .any(|value| !value.is_finite() || *value <= 0.0) {
                            return Err(invalid("rotary fluid/string conjugate port is unrepresentable"));
                        }
                        let mut left = vec![0.0; sources[index].source.modes()];
                        let mut right = vec![0.0; fluid_modes];
                        left[transverse_modes + k] = solid_surface;
                        right[offset + k * flow.branches().len() + branch_index] = fluid_surface;
                        links.push(ModalConnection {
                            left: ModalAttachment { component: index, shapes: left },
                            right: ModalAttachment { component: fluid_start + index, shapes: right },
                            stiffness_n_m: 0.0, damping_n_s_m: resistance, rest_extension_m: 0.0,
                        });
                        modes.push(ModalAcousticMode {
                            angular_frequency_rad_s: 0.0, damping_ratio: 0.0,
                            pressure_per_modal_velocity: fs_math::c64::C64::ZERO,
                        });
                    }
                }
            }
            if modes.len() != fluid_modes {
                return Err(invalid("cylinder memory did not retain every transverse and rotary coordinate"));
            }
            components.push(ModalAcousticTimeModel::try_new_with_free_coordinates(
                options.sample_rate_hz, modes, options.string_budget,
            ).map_err(|source| ReciprocalStringError::Network(
                ModalCouplingError::Component { component: fluid_start + index, source },
            ))?);
        }
        Ok(())
    }
}

impl MaterialStringBankSource {
    /// The actual causal gas law and original solid mass survive scheduling.
    /// Absent loading retains the existing sampled-resistance model.
    #[must_use]
    pub const fn cylinder_gas_loading(&self) -> Option<&StringBankFluidLoading> {
        self.fluid_loading.as_ref()
    }
}

fn transverse_count(source: &MaterialStringSource, card: &BowedStringCard) -> usize {
    card.mode_count * if matches!(source, MaterialStringSource::WithPolarizations(_)) { 2 } else { 1 }
}

fn geometry_and_band(source: &MaterialStringSource)
    -> Result<(f64, (f64, f64)), ReciprocalStringError>
{
    match source {
        MaterialStringSource::WithPolarizations(p) => geometry_and_band(p.source()),
        MaterialStringSource::WithTorsion(t) => geometry_and_band(t.bending()),
        MaterialStringSource::Wound(s) => Ok((s.section().outer_radius_m(), s.section().omega_band_rad_s())),
        // No loss-property band is manufactured for a measured elastic section.
        // The finite sentinel is numerical admission only; source metadata keeps
        // an unstated frequency range as None throughout the bank/performance.
        MaterialStringSource::Measured(s) =>
            Ok((s.radius_m(), s.omega_band_rad_s().unwrap_or((0.0, f64::MAX)))),
        MaterialStringSource::Homogeneous(s) => {
            let string = s.string();
            let band = string.kelvin_voigt_bending.map(|law| law.omega_band_rad_s)
                .or_else(|| string.relaxation_bending.as_ref().map(|law| law.omega_band_rad_s))
                .ok_or_else(|| invalid("cylinder loading needs the original material bending-loss law"))?;
            Ok((s.radius_m(), band))
        }
        MaterialStringSource::Concentric { section, .. } =>
            Ok((section.elastic().outer_radius_m(), section.omega_band_rad_s())),
        MaterialStringSource::ConcentricProny { section, .. } =>
            Ok((section.elastic().outer_radius_m(), section.omega_band_rad_s())),
        MaterialStringSource::ConcentricModalLoss { section, .. } =>
            Ok((section.elastic().outer_radius_m(), section.omega_band_rad_s())),
    }
}

impl ReciprocalStringBank {
    /// Replace sampled air damping with a causal unsteady-cylinder continuum.
    ///
    /// Source solid density is untouched. The returned cards carry the actual
    /// transverse inertia (solid plus instantaneous fluid added mass); their
    /// damping re-evaluates the original material bending law and keeps the independent
    /// tension loss and the collapsed slow fluid tail. Every remaining fluid
    /// contribution is a coupled mass, so it cannot be double-counted as old
    /// per-mode gas resistance. Selected torsion receives its distinct rotary
    /// diffusion spectrum and steady viscous torque in the same solve. The
    /// material's polar inertia and shear storage remain unchanged; gas memory
    /// exchanges torque through the bow's existing surface-velocity port.
    /// Original Prony specimens retain their actual internal strains and use
    /// the loaded transverse mass in their bending-memory projection. Bind an
    /// explicit layered-string spectrum only after selecting gas loading.
    ///
    /// Gas density and viscosity come from source().gas(), radius from each
    /// original circular section. Additional free components follow all physical
    /// bodies, preserving every existing component index and body observation.
    /// One binding is allowed before any string motion/playback. Body initial
    /// states, complete tensor damping, cavity connections, torque bridges and
    /// source interface are retained. Subsequent torque attachment retains fluid.
    ///
    /// The explicit resolution is numerical data, not an accuracy guarantee.
    /// Refine both cutoffs, branch count, string modes and coupled response.
    /// The underlying continuum is small-motion, incompressible, local, circular
    /// and unconfined; it is not the general finite-amplitude atmosphere solver.
    pub fn with_cylinder_gas_loading(
        mut self,
        resolution: CylinderFlowResolution,
        gate: &CancelGate,
    ) -> Result<Self, ReciprocalStringError> {
        poll(gate)?;
        resolution.validate().map_err(ReciprocalStringError::Material)?;
        if self.system.has_nonlinear_potentials() {
            return Err(invalid("bind cylinder gas loading before geometric stretching; a rebuild must not discard its storage"));
        }
        if self.source.fluid_loading.is_some() || self.system.samples_rendered() != 0
            || self.system.components()[..self.cards.len()].iter().flat_map(|c| c.states())
                .any(|s| s.displacement_m_sqrt_kg != 0.0 || s.velocity_m_sqrt_kg_per_s != 0.0) {
            return Err(invalid("cylinder loading requires one cold binding before string motion or playback"));
        }
        // Bound complete dimensions BEFORE allocating quadrature or port arrays.
        let moving_modes = self.source.strings.iter()
            .try_fold(0usize, |n, member| n.checked_add(member.source.modes()))
            .ok_or_else(|| invalid("complete transverse and torsional mode count overflow"))?;
        let added_modes = moving_modes.checked_mul(resolution.intervals)
            .ok_or_else(|| invalid("complete fluid memory dimension overflow"))?;
        let total_modes = self.system.mode_count().checked_add(added_modes)
            .ok_or_else(|| invalid("complete fluid/string/body dimension overflow"))?;
        let old_links = self.source.strings.iter().try_fold(
            self.source.body_connections.len() + self.source.translational_body_connections.len()
                + self.source.rotational_body_connections.len(),
            |n, s| n.checked_add(s.connections.len()),
        ).ok_or_else(|| invalid("fluid/string/body connection count overflow"))?;
        if total_modes > self.source.options.network.max_modes
            || old_links.checked_add(added_modes)
                .is_none_or(|n| n > self.source.options.network.max_connections) {
            return Err(invalid("cylinder memory exceeds the declared COMPLETE network mode or connection budget"));
        }
        let mut cards = Vec::with_capacity(self.cards.len());
        let mut loading = Vec::with_capacity(self.cards.len());
        for (solid, member) in self.cards.iter().zip(&self.source.strings) {
            poll(gate)?;
            let (radius, (lo, hi)) = geometry_and_band(&member.source)?;
            let flow = CylinderFlow::from_gas(radius, self.source.gas, resolution, gate)
                .map_err(ReciprocalStringError::Material)?;
            let mut card = solid.clone();
            card.linear_density_kg_m = solid.linear_density_kg_m + flow.added_mass_per_length_kg_m();
            if !card.linear_density_kg_m.is_finite()
                || card.linear_density_kg_m <= solid.linear_density_kg_m {
                return Err(invalid("added gas inertia is not representable beside the actual solid line mass"));
            }
            let mut source_string = member.source.resolved_transverse_string()?;
            source_string.lin_density_kg_m = card.linear_density_kg_m;
            for k in 0..card.mode_count {
                let wave_number = (k + 1) as f64 * core::f64::consts::PI / card.length_m;
                let omega = card.mode_omega_rad_s(k + 1);
                // Re-derive intrinsic loss from its source; subtracting two close
                // sampled damping values could silently discard real viscosity.
                let material = crate::acoustic_realize::intrinsic_string_mode_zeta(
                    &source_string, omega, wave_number,
                ).map_err(ReciprocalStringError::Material)?;
                let fluid_tail = 0.5 * flow.low_tail_resistance_per_length_n_s_m2()
                    / card.linear_density_kg_m / omega;
                // Loss-angle/TED is evaluated again at the loaded mode, using
                // its actual inertia. Keeping the dry sampled zeta would erase
                // thermal dispersion and give the wrong tension dilution.
                let modal_rate = member.source.modal_bending_loss().map(|section|
                    section.damping_rate_per_s(omega, wave_number, card.linear_density_kg_m))
                    .transpose().map_err(|error| ReciprocalStringError::Material(
                        crate::acoustic_realize::AcousticRealizeError::Nonlinear(error.to_string())))?
                    .unwrap_or(0.0);
                let zeta = material + 0.5 * modal_rate / omega + fluid_tail;
                if !omega.is_finite() || omega <= 0.0 || omega < lo || omega > hi
                    || !zeta.is_finite() || zeta <= 0.0 {
                    return Err(invalid("fluid-loaded string mode exceeds its original material band or representability"));
                }
                card.zetas[k] = zeta;
            }
            card.validate().map_err(|_| invalid("fluid-loaded physical string card is not admissible"))?;
            cards.push(card);
            let rotary = member.source.torsion_source().map(|_| {
                RotaryCylinderFlow::from_gas(radius, self.source.gas, resolution, gate)
                    .map_err(ReciprocalStringError::Material)
            }).transpose()?;
            loading.push(StringCylinderLoading { solid_card: solid.clone(), flow, rotary });
        }
        let bodies = self.system.components()[self.cards.len()..].to_vec();
        let fluid = StringBankFluidLoading { resolution, strings: loading,
            body_components: bodies.len() };
        let bridges: Vec<_> = self.source.strings.iter().map(|s| s.connections.as_slice()).collect();
        let mut system = assemble_sources(
            &cards, Some(&self.source.strings), bodies, &bridges,
            &self.source.body_connections, &self.source.translational_body_connections,
            &self.source.rotational_body_connections,
            Some(&fluid), self.source.prepared_linear, &self.source.constant_loss_bridges(),
            self.source.residual_flexibility, self.system.port().law.clone(), &self.source.options, gate,
        )?;
        if let Some(original) = self.system.acceleration_pressure_weights() {
            let mut weights = original.to_vec();
            let mut offset = 0;
            for (i, member) in self.source.strings.iter().enumerate() {
                let scale = det::sqrt(self.cards[i].linear_density_kg_m / cards[i].linear_density_kg_m);
                for weight in &mut weights[offset..offset + transverse_count(&member.source, &cards[i])] {
                    *weight *= scale;
                }
                offset += member.source.modes();
            }
            weights.resize(system.mode_count(), 0.0);
            system = system.with_acceleration_pressure(weights).map_err(ReciprocalStringError::Friction)?;
        }
        if let Some(binding) = self.system.interface() {
            system = system.with_interface(binding.clone()).map_err(ReciprocalStringError::Friction)?;
        }
        system = self.source.restore_material_memory(system, &cards,
            Some(self.system.material_memory().to_vec()), gate)?;
        self.source.validate_material_bridge_frequencies(&system, gate)?;
        poll(gate)?;
        self.cards = cards;
        self.source.fluid_loading = Some(fluid);
        self.system = system;
        Ok(self)
    }
}
