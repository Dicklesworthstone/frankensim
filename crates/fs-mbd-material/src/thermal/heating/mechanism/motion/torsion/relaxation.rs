//! Per-layer Maxwell memory, coupled through the original shaft and RK stages.
use super::*;
use fs_material::section::relaxation::{LayerRelaxationBranch, LayerRelaxationProperties};

/// Original equilibrium/arm source roles and explicit initial spring history.
///
/// A spring twist is the elastic part of that arm's total end rotation, in
/// layer/declaration order. It is NOT inferred from the current joint angle:
/// zero describes relaxed arms; an instantaneous initial twist is supplied by
/// the caller. Its initial stored energy belongs to the initial mechanical state.
#[derive(Clone, Debug)]
pub struct ShaftRelaxation {
    pub properties: Vec<LayerRelaxationProperties>,
    pub initial_spring_twists_rad: Vec<f64>,
    /// Absolute RK accuracy scale for each spring twist [rad].
    pub absolute_spring_twist_tolerance_rad: f64,
}

/// Original integrated source arms and actual accepted elastic spring twists.
/// Original source answers remain in the containing shaft's section layers.
#[derive(Clone, Debug, PartialEq)]
pub struct ShaftRelaxationState {
    pub branches: Vec<LayerRelaxationBranch>,
    pub spring_twists_rad: Vec<f64>,
}

impl BoundTorsion<'_> {
    pub fn install_relaxation(
        &mut self,
        definition: ShaftRelaxation,
        q: f64,
        rate: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Error> {
        if !positive(definition.absolute_spring_twist_tolerance_rad) {
            return Err(Error::Invalid(
                "Maxwell spring twists need a positive RK accuracy scale",
            ));
        }
        let state = self.resolve_relaxation(
            &definition,
            &self.state.section,
            &definition.initial_spring_twists_rad,
            checkpoint,
        )?;
        self.state.relaxation = Some(state);
        self.relaxation = Some(definition);
        self.update(q, rate)?;
        checkpoint()?;
        Ok(())
    }

    pub fn spring_twists(&self) -> &[f64] {
        self.state
            .relaxation
            .as_ref()
            .map_or(&[], |memory| memory.spring_twists_rad.as_slice())
    }

    pub fn spring_tolerance(&self) -> f64 {
        self.relaxation.as_ref().map_or(1.0, |definition| {
            definition.absolute_spring_twist_tolerance_rad
        })
    }

    /// The section owns Delta_G * J; only the actual reference length converts
    /// it to a generalized rotational spring. The constitutive owner supplies
    /// the differential strain and heat law, without an extra time integrator.
    pub fn memory_rate(&self, index: usize, twist: f64, rate: f64) -> Result<(f64, f64), Error> {
        let arm = self
            .state
            .relaxation
            .as_ref()
            .and_then(|s| s.branches.get(index))
            .ok_or(Error::Invalid("Maxwell shaft arm is absent"))?;
        let (derivative, normalized_heat) =
            fs_material::visco::prony_rate(twist, rate, arm.relaxation_time_s).map_err(|_| {
                Error::Invalid("Maxwell shaft spring rate or heat is outside its domain")
            })?;
        let stiffness = finite(arm.relaxing_stiffness_n_m2 / self.definition.reference_length_m)?;
        let heat = finite(stiffness * normalized_heat)?;
        if arm.relaxing_stiffness_n_m2 != 0.0
            && (stiffness == 0.0 || normalized_heat != 0.0 && heat == 0.0)
        {
            return Err(Error::Invalid(
                "Maxwell shaft stiffness or physical power underflow",
            ));
        }
        Ok((derivative, heat))
    }

    /// Effort, nonnegative power and spring storage from the SAME trial memory.
    pub fn memory_response(&self, twists: &[f64]) -> Result<(f64, f64, f64), Error> {
        let Some(memory) = &self.state.relaxation else {
            return if twists.is_empty() {
                Ok((0.0, 0.0, 0.0))
            } else {
                Err(Error::Invalid("shaft has no Maxwell memory"))
            };
        };
        if twists.len() != memory.branches.len() {
            return Err(Error::Invalid(
                "one elastic twist per original Maxwell shaft arm required",
            ));
        }
        let length = self.definition.reference_length_m;
        let mut torque = 0.0;
        let mut heat = 0.0;
        let mut energy = 0.0;
        for (i, (arm, &twist)) in memory.branches.iter().zip(twists).enumerate() {
            let radius =
                self.state.section.shear().elastic().layers()[arm.layer_index].outer_radius_m;
            let shear = finite((twist / length) * radius)?;
            if shear.abs() > self.definition.maximum_absolute_shear_strain
                || twist != 0.0 && shear == 0.0
            {
                return Err(Error::Invalid(
                    "Maxwell shaft spring left its small-shear domain",
                ));
            }
            let stiffness = finite(arm.relaxing_stiffness_n_m2 / length)?;
            let effort = finite(stiffness * twist)?;
            let stored = finite(0.5 * effort * twist)?;
            if arm.relaxing_stiffness_n_m2 != 0.0
                && (stiffness == 0.0 || twist != 0.0 && stored == 0.0)
            {
                return Err(Error::Invalid("Maxwell shaft storage underflow"));
            }
            torque = finite(torque - effort)?;
            heat = finite(heat + self.memory_rate(i, twist, 0.0)?.1)?;
            energy = finite(energy + stored)?;
        }
        Ok((torque, heat, energy))
    }

    pub fn accept_spring_twists(&mut self, twists: Vec<f64>) -> Result<(), Error> {
        self.memory_response(&twists)?;
        if let Some(memory) = &mut self.state.relaxation {
            memory.spring_twists_rad = twists;
        }
        Ok(())
    }

    /// Every branch deposits its own integrated dashpot work into its original
    /// finite layer cell. Different time constants are never averaged.
    pub fn deposit_memory_losses(
        &self,
        losses: &[f64],
        dt_s: f64,
        heat_w: &mut BTreeMap<CellId, f64>,
    ) -> Result<f64, Error> {
        if losses.len() != self.spring_twists().len() || !positive(dt_s) {
            return Err(Error::Invalid(
                "Maxwell heat needs the original arm count and positive interval",
            ));
        }
        let mut total = 0.0;
        if let Some(memory) = &self.state.relaxation {
            for (arm, &loss) in memory.branches.iter().zip(losses) {
                if !nonnegative(loss) || arm.relaxing_stiffness_n_m2 == 0.0 && loss != 0.0 {
                    return Err(Error::Budget(
                        "Maxwell integrated physical heat is invalid; refine RK tolerances",
                    ));
                }
                let watts = heat_w
                    .entry(self.definition.cells[arm.layer_index])
                    .or_insert(0.0);
                *watts = finite(*watts + loss / dt_s)?;
                total = finite(total + loss)?;
            }
        }
        Ok(total)
    }

    pub(super) fn resolve_relaxation(
        &self,
        definition: &ShaftRelaxation,
        section: &ConcentricTorsionalSection,
        twists: &[f64],
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<ShaftRelaxationState, Error> {
        checkpoint()?;
        let source = section.clone().with_prony_torsion(&definition.properties)?;
        if twists.len() != source.branches().len()
            || source.omega_band_rad_s().0 != 0.0
            || source.omega_band_rad_s().1 < self.definition.maximum_angular_frequency_rad_s
        {
            return Err(Error::Invalid(
                "shaft memory needs every original arm and complete DC-to-band source support",
            ));
        }
        for (((current, reference), region), properties) in section
            .shear()
            .elastic()
            .layers()
            .iter()
            .zip(self.definition.law.reference().layers())
            .zip(&self.inventory.regions)
            .zip(&definition.properties)
        {
            checkpoint()?;
            let original = reference.material.query();
            let actual = current.material.query();
            let (lower, upper) = if actual.axes()["T"] < original.axes()["T"] {
                (&actual, &original)
            } else {
                (&original, &actual)
            };
            // Keep the exact original source pin AND its supported path, not
            // just a positive endpoint or an independently filled state box.
            let roles = std::iter::once((&properties.equilibrium_modulus, false)).chain(
                properties
                    .branches
                    .iter()
                    .flat_map(|(delta, tau)| [(delta, true), (tau, false)]),
            );
            for (name, permits_zero) in roles {
                let property = current
                    .material
                    .property(name)
                    .ok_or(Error::Invalid("Maxwell shaft lost an original source role"))?;
                let claim = region
                    .card
                    .claims()
                    .claim(property.answer().receipt.selected)
                    .ok_or(Error::Invalid("Maxwell shaft source is absent"))?;
                let path = region
                    .card
                    .claims()
                    .query_line_typed(
                        &claim.key,
                        "T",
                        lower,
                        upper,
                        ClaimSelection::Pinned(property.answer().receipt.selected),
                    )
                    .map_err(|_| {
                        Error::Invalid(
                            "Maxwell shaft source does not support its complete thermal path",
                        )
                    })?;
                if path.knots.iter().any(|(_, answer)| {
                    let value = answer.evidence.value.value;
                    if permits_zero {
                        !nonnegative(value)
                    } else {
                        !positive(value)
                    }
                }) {
                    return Err(Error::Invalid(
                        "Maxwell shaft coefficients leave their physical domain along the thermal path",
                    ));
                }
            }
        }
        Ok(ShaftRelaxationState {
            branches: source.branches().to_vec(),
            spring_twists_rad: twists.to_vec(),
        })
    }
}
