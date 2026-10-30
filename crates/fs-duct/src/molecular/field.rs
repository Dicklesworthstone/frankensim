//! Bind the existing primitive-field quadrature to the same molecular duct.
//! Geometry never gets a second mesh, thermodynamic coefficients are never
//! interpolated, and the scalar/exterior response owners remain unchanged.

use super::*;
use crate::profile::RiccatiSteps;
use crate::profile::gas_field::{GasFieldDuct, GasFieldDuctError, MolecularFieldStep};
use fs_material::gas::GasState;
use fs_material::gas_field::AxialGasField;
use fs_material::relaxation::{RelaxationAbsorption, RelaxingMedium};

#[derive(Debug, Clone)]
pub(super) struct FieldBinding {
    field: AxialGasField,
    law: RelaxationAbsorption,
    outlet: RelaxationState,
    /// Index zero is the exact inlet state, reused by frequency preparation.
    pub(super) states: Vec<RelaxationState>,
}

// Exact bit keys, never rounded bins. The common declared source/law survives
// in FieldBinding. Reservoir identities are retained even when two responses
// happen to coincide numerically, so source distinctions cannot disappear.
type GasKey = [u64; 12];
type StateKey = (GasKey, u64, Vec<(String, usize, u64, u64)>);

fn gas_key(gas: &GasState) -> GasKey {
    [
        gas.temperature,
        gas.pressure,
        gas.density,
        gas.sound_speed,
        gas.dynamic_viscosity,
        gas.thermal_conductivity,
        gas.gamma,
        gas.specific_gas_constant,
        gas.specific_heat_cp,
        gas.prandtl,
        gas.characteristic_impedance,
        gas.water_mole_fraction,
    ]
    .map(f64::to_bits)
}

fn state_key(state: &RelaxationState) -> StateKey {
    (
        gas_key(state.gas()),
        state.frozen_cv_j_kg_k().to_bits(),
        state
            .processes()
            .iter()
            .map(|p| {
                (
                    p.species.to_string(),
                    p.index,
                    p.heat_capacity_j_kg_k.to_bits(),
                    p.relaxation_time_s.to_bits(),
                )
            })
            .collect(),
    )
}

fn evaluate(
    field: &AxialGasField,
    law: &RelaxationAbsorption,
    x: f64,
) -> Result<RelaxationState, MolecularDuctError> {
    Ok(law.evaluate(&RelaxingMedium::from_thermochemical(&field.state_at(x)?)?)?)
}

impl From<GasFieldDuctError> for MolecularDuctError {
    fn from(error: GasFieldDuctError) -> Self {
        match error {
            GasFieldDuctError::Material(error) => Self::Material(error),
            GasFieldDuctError::Duct(error) => Self::Duct(error),
        }
    }
}

impl MolecularDuct {
    /// Resolve temperature, composition and molecular relaxation at every
    /// original gas-field quadrature sample and tone-hole station.
    ///
    /// The primitive field and declared relaxation law retain their sources.
    /// The complete field must cover the exact duct length at common static
    /// pressure. Nonuniform cylinders, cones and profiles use the original
    /// fourth-order Webster/Magnus geometry/field quadrature, split at every
    /// geometry boundary, field station and profile knot. `steps_per_span`
    /// refines that operator; it is not a three-dimensional acoustic ladder.
    /// Uniform fields delegate `new` verbatim, including `uniform_cones` and
    /// each original profile's own step axis.
    ///
    /// Every distinct admitted material state evaluates its caloric response
    /// once per frequency. Diffusion spectra are prepared once per exact gas
    /// state/radius pair. No interpolation of density, transport, cp or rates
    /// occurs, and no molecular or wall loss is added twice. Wall pole count
    /// and axial steps remain independently refinable.
    ///
    /// Tone holes retain the existing homogeneous local chimney and compact
    /// junction law. A varying primitive field between multiple open holes
    /// refuses because their compact inner mutual correction has no such law.
    /// All response APIs use these same local states.
    /// This is a prescribed frozen background: mean flow, axial heat/species
    /// transport, condensation and reservoir remapping remain separate owners.
    pub fn from_gas_field(
        duct: &Duct,
        field: &AxialGasField,
        law: &RelaxationAbsorption,
        wall: Option<MolecularWallResolution>,
        steps_per_span: RiccatiSteps,
        uniform_cones: ConeSlices,
        maximum_blocks: usize,
    ) -> Result<Self, MolecularDuctError> {
        Self::from_gas_field_controlled(
            duct,
            field,
            law,
            wall,
            steps_per_span,
            uniform_cones,
            maximum_blocks,
            &mut || false,
        )
    }

    /// Cancellable cold binding; cancellation/refusal returns no partial duct.
    #[allow(clippy::too_many_arguments)]
    pub fn from_gas_field_controlled(
        duct: &Duct,
        field: &AxialGasField,
        law: &RelaxationAbsorption,
        wall: Option<MolecularWallResolution>,
        steps_per_span: RiccatiSteps,
        uniform_cones: ConeSlices,
        maximum_blocks: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, MolecularDuctError> {
        Self::from_gas_field_bound_controlled(
            duct, field, law, wall, None, steps_per_span, uniform_cones,
            maximum_blocks, cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn from_gas_field_bound_controlled(
        duct: &Duct,
        field: &AxialGasField,
        law: &RelaxationAbsorption,
        wall: Option<MolecularWallResolution>,
        compliant_wall: Option<MolecularCompliantWall>,
        steps_per_span: RiccatiSteps,
        uniform_cones: ConeSlices,
        maximum_blocks: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, MolecularDuctError> {
        poll(cancelled, 0)?;
        if maximum_blocks == 0 || maximum_blocks > compose::MAX_COMPOSITION_BLOCKS {
            return Err(bad("molecular duct composition block cap is invalid").into());
        }
        // The existing owner admits coverage, primitive material evaluations,
        // geometric breaks, exact hole stations, and the inner mutual law.
        let geometry = GasFieldDuct::new_controlled(duct, field, steps_per_span, cancelled)?;
        let inlet = evaluate(field, law, 0.0)?;
        let outlet = evaluate(field, law, field.length_m())?;
        for &x in field.stations_m() {
            poll(cancelled, 0)?;
            evaluate(field, law, x)?;
        }
        if field.uniform_state().is_some() {
            let mut result = Self::new_bound_controlled(
                duct, inlet, wall, compliant_wall, uniform_cones, maximum_blocks, cancelled,
            )?;
            result.field_binding = Some(FieldBinding {
                field: field.clone(),
                law: law.clone(),
                outlet,
                states: Vec::new(),
            });
            return Ok(result);
        }
        let count = geometry.molecular_steps().try_fold(0usize, |total, step| {
            let blocks = match step {
                MolecularFieldStep::Magnus { .. } => 1,
                MolecularFieldStep::Hole { .. } => 3,
            };
            total
                .checked_add(blocks)
                .filter(|&n| n <= maximum_blocks)
                .ok_or_else(|| bad("molecular field composition block work cap exceeded"))
        })?;
        let reference = geometry.inlet_impedance();
        if !positive(reference) {
            return Err(bad("molecular field reference impedance is unrepresentable").into());
        }
        let mut state_cache = BTreeMap::new();
        state_cache.insert(state_key(&inlet), 0);
        let mut result = Self {
            state: inlet.clone(),
            field_binding: Some(FieldBinding {
                field: field.clone(),
                law: law.clone(),
                outlet,
                states: vec![inlet],
            }),
            points: Vec::new(),
            pieces: Vec::with_capacity(count),
            reference_impedance: reference,
            compliant_wall,
        };
        let mut point_cache = BTreeMap::new();
        let mut wall_cache: BTreeMap<(u64, GasKey), ZkDiffusive> = BTreeMap::new();
        poll(cancelled, 0)?;
        // All physical circular diffusion operators share one dimensionless
        // radial eigensystem; only local scaling and Bessel admission repeat.
        let basis = wall
            .map(|r| {
                fs_phs::ZkDiffusiveBasis::new(r.poles_per_field).map_err(|reason| {
                    DuctError::WallDiffusionUnresolved {
                        reason,
                        radius_m: duct.segments[0].inlet_radius(),
                        x_m: 0.0,
                    }
                })
            })
            .transpose()?;
        poll(cancelled, 0)?;
        let open_sites = crate::collect_open_hole_sites(duct);
        for step in geometry.molecular_steps() {
            poll(cancelled, result.pieces.len())?;
            let mut sample = |radius: f64, x: f64| -> Result<usize, MolecularDuctError> {
                poll(cancelled, result.pieces.len())?;
                let state = evaluate(field, law, x)?;
                let key = state_key(&state);
                let binding = result.field_binding.as_mut().expect("field binding");
                let state_index = if let Some(&index) = state_cache.get(&key) {
                    index
                } else {
                    let index = binding.states.len();
                    binding.states.push(state);
                    state_cache.insert(key, index);
                    index
                };
                let key = (radius.to_bits(), state_index);
                if let Some(&index) = point_cache.get(&key) {
                    return Ok(index);
                }
                let gas = binding.states[state_index].gas();
                let wall = if let Some(r) = wall {
                    let key = (radius.to_bits(), gas_key(gas));
                    if let Some(model) = wall_cache.get(&key) {
                        Some(model.clone())
                    } else {
                        let model = basis
                            .as_ref()
                            .expect("wall basis")
                            .bind(
                                radius,
                                gas.density,
                                gas.dynamic_viscosity,
                                gas.gamma,
                                gas.prandtl,
                                r.band_rad_s,
                                r.maximum_relative_error,
                            )
                            .map_err(|reason| DuctError::WallDiffusionUnresolved {
                                reason,
                                radius_m: radius,
                                x_m: x,
                            })?;
                        wall_cache.insert(key, model.clone());
                        Some(model)
                    }
                } else {
                    None
                };
                let index = result.points.len();
                result.points.push(MaterialPoint {
                    radius,
                    x_m: x,
                    wall,
                    state: Some(state_index),
                });
                point_cache.insert(key, index);
                Ok(index)
            };
            let piece = match step {
                MolecularFieldStep::Magnus {
                    h,
                    radii,
                    positions_m,
                    slants,
                } => Piece::Magnus {
                    h,
                    points: [
                        sample(radii[0], positions_m[0])?,
                        sample(radii[1], positions_m[1])?,
                    ],
                    slants,
                },
                MolecularFieldStep::Hole {
                    geometry,
                    axial_position_m,
                } => {
                    let Segment::ToneHole { hole_radius, .. } = geometry else {
                        unreachable!()
                    };
                    let point = sample(*hole_radius, axial_position_m)?;
                    Piece::ToneHole(holes::ToneHole::new(
                        geometry,
                        point,
                        axial_position_m,
                        &open_sites,
                    )?)
                }
            };
            result.pieces.push(piece);
        }
        poll(cancelled, count)?;
        Ok(result)
    }

    /// Original primitive material field, when this duct was field-bound.
    pub fn gas_field(&self) -> Option<&AxialGasField> {
        self.field_binding.as_ref().map(|binding| &binding.field)
    }

    /// Original source-backed molecular law for a field-bound duct.
    pub fn field_relaxation_law(&self) -> Option<&RelaxationAbsorption> {
        self.field_binding.as_ref().map(|binding| &binding.law)
    }

    /// Exact outlet molecular material; this does not specify exterior gas.
    pub fn outlet_state(&self) -> &RelaxationState {
        self.field_binding
            .as_ref()
            .map(|binding| &binding.outlet)
            .unwrap_or(&self.state)
    }

    /// Exact admitted material samples, counting the inlet once. Multiple
    /// geometry radii in one state share its frequency-domain caloric sum.
    pub fn material_states(&self) -> usize {
        self.field_binding
            .as_ref()
            .map_or(1, |binding| binding.states.len().max(1))
    }
}

#[cfg(test)]
#[path = "field_tests.rs"]
mod tests;
