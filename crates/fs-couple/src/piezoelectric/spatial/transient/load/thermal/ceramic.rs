//! Complete original-cell heat routes into the existing finite caloric owner.
//!
//! Several tetrahedra may share a thermal cell only when their combined
//! reference volume and mass match that cell's actual inventory. Every source
//! retains the same caloric card, density claim and non-temperature conditions.
//! The material owner admits temperature changes only for coefficients that
//! remain independent of temperature within their original source validity.
//! Source-bound thermal samples additionally retain their exact cell temperature:
//! this route cannot freeze evolving thermal strain or pyroelectric polarization.

use super::{LoadThermalError, Result};
use crate::circuit_enthalpy::{
    CircuitEnthalpyError,
    material::{MaterialNetworkError, MaterialNetworkSample, PreparedMaterialNetwork},
};
use crate::piezoelectric::spatial::{TetPiezoelectricError, transient::TetPiezoelectricTransient};
use fs_blake3::ContentHash;
use fs_exec::Cx;
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation, ThermalLaw};
use std::collections::BTreeMap;

/// Whole-cell intrinsic heat from one original piezoelectric tetrahedron.
/// Each original tetrahedron must appear exactly once. Shared receiving cells
/// retain the sum of their original volumes and masses; no heat fraction or
/// inferred thermal mixing is introduced.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ThermalCeramicCell {
    /// Index in the original piezoelectric tetrahedral connectivity.
    pub piezoelectric_cell: usize,
    /// Stable active thermal control volume receiving both intrinsic losses.
    pub thermal_cell: CellId,
}

// This only accommodates rounding of the SAME authored geometric inventory.
// It is not an allowance for adding a heat sink or changing ceramic density.
const INVENTORY_RELATIVE_TOLERANCE: f64 = 1e-10;

struct Inventory {
    cell: CellId,
    volume_m3: f64,
    mass_kg: f64,
}

struct MaterialReceiver {
    original_cell: usize,
    material_index: usize,
    thermal_cell: CellId,
    curve_identity: ContentHash,
}

pub(super) struct CeramicPlan {
    destinations: Vec<ThermalCeramicCell>,
    temperature_axis: String,
    inventory: Vec<Inventory>,
    groups: Vec<MaterialReceiver>,
}

fn invalid(message: &'static str) -> CircuitEnthalpyError<LoadThermalError> {
    CircuitEnthalpyError::Invalid(message)
}

fn check(cx: &Cx<'_>) -> Result<()> {
    cx.checkpoint().map_err(|_| CircuitEnthalpyError::Cancelled)
}

fn material_error(error: MaterialNetworkError) -> CircuitEnthalpyError<LoadThermalError> {
    CircuitEnthalpyError::Source(LoadThermalError::Material(error))
}

fn same_inventory(actual: f64, expected: f64) -> bool {
    actual.is_finite()
        && actual > 0.0
        && expected.is_finite()
        && expected > 0.0
        && (actual / expected - 1.0).abs() <= INVENTORY_RELATIVE_TOLERANCE
}

impl CeramicPlan {
    pub(super) fn prepare(
        body: &TetPiezoelectricTransient<'_, '_>,
        thermal: &DepositionSimulation<'_>,
        material: &PreparedMaterialNetwork<'_, '_>,
        bindings: &[ThermalCeramicCell],
        temperature_axis: &str,
        cx: &Cx<'_>,
    ) -> Result<Self> {
        check(cx)?;
        let Some(intrinsic) = &body.intrinsic else {
            if !bindings.is_empty() {
                return Err(invalid(
                    "ceramic heat routes require retained intrinsic material memories",
                ));
            }
            return Ok(Self {
                destinations: Vec::new(),
                temperature_axis: String::new(),
                inventory: Vec::new(),
                groups: Vec::new(),
            });
        };
        let count = body.model.tetrahedra.len();
        if bindings.len() != count {
            return Err(invalid(
                "each original piezoelectric cell requires exactly one thermal destination",
            ));
        }
        let mut destinations = vec![None; count];
        for binding in bindings {
            check(cx)?;
            if binding.piezoelectric_cell >= count
                || destinations[binding.piezoelectric_cell]
                    .replace(*binding)
                    .is_some()
                || !thermal.cells().contains_key(&binding.thermal_cell)
            {
                return Err(invalid(
                    "ceramic heat routes must be complete, unique and active",
                ));
            }
        }
        let destinations = destinations
            .into_iter()
            .map(|route| {
                route.ok_or_else(|| {
                    invalid("original piezoelectric cell has no thermal destination")
                })
            })
            .collect::<Result<Vec<_>>>()?;
        let mut by_cell = BTreeMap::<CellId, Inventory>::new();
        let mut representatives = BTreeMap::new();
        for route in &destinations {
            check(cx)?;
            let (material_index, source) = intrinsic.material(route.piezoelectric_cell);
            representatives
                .entry((route.thermal_cell, material_index))
                .or_insert(route.piezoelectric_cell);
            let volume = body.model.dielectric.volumes_m3()[route.piezoelectric_cell];
            let mass = volume * source.base().density_kg_m3();
            let inventory = by_cell.entry(route.thermal_cell).or_insert(Inventory {
                cell: route.thermal_cell,
                volume_m3: 0.0,
                mass_kg: 0.0,
            });
            inventory.volume_m3 += volume;
            inventory.mass_kg += mass;
            if !volume.is_finite()
                || volume <= 0.0
                || !mass.is_finite()
                || mass <= 0.0
                || !inventory.volume_m3.is_finite()
                || !inventory.mass_kg.is_finite()
            {
                return Err(invalid("ceramic thermal inventory is not representable"));
            }
        }
        let mut caloric = BTreeMap::new();
        for binding in material.cell_bindings() {
            check(cx)?;
            if caloric.insert(binding.cell, binding).is_some() {
                return Err(invalid("duplicate original thermal material cell"));
            }
        }
        let mut groups = Vec::with_capacity(representatives.len());
        for ((thermal_cell, material_index), original_cell) in representatives {
            check(cx)?;
            let (_, source) = intrinsic.material(original_cell);
            let base = source.base().material();
            let density = source.base().density_kg_m3();
            let binding = caloric.get(&thermal_cell).ok_or_else(|| {
                invalid("ceramic receiver lacks its original caloric material binding")
            })?;
            let cell = &thermal.cells()[&thermal_cell];
            let curve = binding.curve;
            if binding.point.temperature_axis != temperature_axis
                || cell.material_id != binding.material_id
                || !thermal.uses_material_law(cell.material_id, curve as &dyn ThermalLaw)
                || curve.material_card_identity() != base.card_identity()
                || binding.point.material.content_hash() != base.card_identity()
            {
                return Err(invalid(
                    "ceramic and caloric owners must retain the same card, temperature axis and actual law",
                ));
            }
            let phase = curve
                .state_at_specific_enthalpy(cell.specific_h_j_kg)
                .map_err(|error| material_error(MaterialNetworkError::Phase(error)))?;
            let query = binding
                .point
                .at(phase)
                .map_err(|error| material_error(MaterialNetworkError::Material(error)))?;
            let expected = base
                .query()
                .with(temperature_axis, phase.temperature_k())
                .map_err(|_| invalid("ceramic thermal source query is not representable"))?;
            if query != expected {
                return Err(invalid(
                    "ceramic and caloric owners must retain identical non-temperature conditions",
                ));
            }
            let density_claim = base
                .property("density")
                .ok_or_else(|| invalid("original ceramic density source is absent"))?
                .answer()
                .receipt
                .selected;
            let receipts = curve.source_receipts();
            if receipts.is_empty()
                || receipts.len() % 2 != 0
                || receipts
                    .chunks_exact(2)
                    .any(|pair| pair[1].property != "density" || pair[1].selected != density_claim)
                || receipts.iter().any(|receipt| {
                    receipt.axis_quantities != *base.axis_quantities()
                        || !receipt
                            .query_point
                            .iter()
                            .filter(|(axis, _)| axis != temperature_axis)
                            .eq(base
                                .query_point()
                                .iter()
                                .filter(|(axis, _)| axis != temperature_axis))
                        || receipt.property == "density" && receipt.selected != density_claim
                })
            {
                return Err(invalid(
                    "ceramic caloric receipts must retain original density and complete fixed conditions",
                ));
            }
            if curve.knots().iter().any(|knot| {
                knot.liquid_mass_fraction != 0.0
                    || knot.bulk_density_kg_m3.to_bits() != density.to_bits()
            }) {
                return Err(invalid(
                    "fixed ceramic inventory requires a solid caloric chart with the original constant density",
                ));
            }
            groups.push(MaterialReceiver {
                original_cell,
                material_index,
                thermal_cell,
                curve_identity: curve.identity(),
            });
        }
        let plan = Self {
            destinations,
            temperature_axis: temperature_axis.to_owned(),
            inventory: by_cell.into_values().collect(),
            groups,
        };
        plan.admit_inventory(thermal, cx)?;
        check(cx)?;
        Ok(plan)
    }

    fn admit_inventory(&self, thermal: &DepositionSimulation<'_>, cx: &Cx<'_>) -> Result<()> {
        for inventory in &self.inventory {
            check(cx)?;
            let cell = thermal
                .cells()
                .get(&inventory.cell)
                .ok_or_else(|| invalid("ceramic thermal destination is inactive"))?;
            if !same_inventory(cell.geometry.volume_m3, inventory.volume_m3)
                || !same_inventory(cell.mass_kg, inventory.mass_kg)
            {
                return Err(invalid(
                    "ceramic thermal cell must retain its aggregate original volume and mass",
                ));
            }
        }
        Ok(())
    }

    pub(super) fn admit(
        &self,
        body: &TetPiezoelectricTransient<'_, '_>,
        thermal: &DepositionSimulation<'_>,
        sample: &MaterialNetworkSample,
        cx: &Cx<'_>,
    ) -> Result<()> {
        check(cx)?;
        let Some(intrinsic) = &body.intrinsic else {
            return if self.destinations.is_empty() {
                Ok(())
            } else {
                Err(invalid(
                    "ceramic heat plan lost its original material memory",
                ))
            };
        };
        if self.destinations.len() != body.model.tetrahedra.len() {
            return Err(invalid(
                "ceramic heat plan must cover the complete original body",
            ));
        }
        self.admit_inventory(thermal, cx)?;
        let mut phases = BTreeMap::new();
        for cell in &sample.cells {
            check(cx)?;
            if phases.insert(cell.cell, cell.phase).is_some() {
                return Err(invalid("duplicate thermal material phase sample"));
            }
        }
        // A retained thermal sample owns its actual temperature even when
        // its current expansion and pyroelectric offsets vanish at T_ref.
        // This heat route advances only fixed constitutive references; changing
        // their temperature requires reversible thermal work and memory-aware
        // tensor evolution. Check every original cell, not only source groups.
        if let Some(samples) = body.model.thermal_materials() {
            for route in &self.destinations {
                check(cx)?;
                let original = &samples[if samples.len() == 1 {
                    0
                } else {
                    route.piezoelectric_cell
                }];
                let phase = phases.get(&route.thermal_cell).ok_or_else(|| {
                    invalid("ceramic receiver lacks its actual thermal phase sample")
                })?;
                if phase.temperature_k().to_bits() != original.temperature_k().to_bits() {
                    return Err(invalid(
                        "retained thermal strain and polarization require their original cell temperature",
                    ));
                }
            }
        }
        for group in &self.groups {
            check(cx)?;
            let (index, material) = intrinsic.material(group.original_cell);
            let phase = phases
                .get(&group.thermal_cell)
                .ok_or_else(|| invalid("ceramic receiver lacks its actual thermal phase sample"))?;
            let cell = &thermal.cells()[&group.thermal_cell];
            if index != group.material_index
                || phase.material_card_identity() != material.base().material().card_identity()
                || phase.phase_curve_identity() != group.curve_identity
                || phase.specific_enthalpy_j_kg().to_bits() != cell.specific_h_j_kg.to_bits()
                || !phase.is_fully_solid()
                || phase.bulk_density_kg_m3().to_bits() != material.base().density_kg_m3().to_bits()
            {
                return Err(invalid(
                    "ceramic heat requires the same current solid caloric source and reference density",
                ));
            }
            material
                .admit_held_temperature(&self.temperature_axis, phase.temperature_k())
                .map_err(|error| {
                    CircuitEnthalpyError::Source(LoadThermalError::Piezoelectric(
                        TetPiezoelectricError::Material(error),
                    ))
                })?;
        }
        check(cx)
    }

    pub(super) fn deposit(
        &self,
        body: &TetPiezoelectricTransient<'_, '_>,
        dt: f64,
        powers: &mut BTreeMap<CellId, f64>,
        checkpoint: &mut dyn FnMut() -> core::result::Result<(), DepositionError>,
    ) -> Result<()> {
        checkpoint()?;
        if !dt.is_finite() || dt <= 0.0 {
            return Err(invalid(
                "ceramic heat deposition requires a positive finite interval",
            ));
        }
        let Some(steps) = body.material_relaxation_steps() else {
            return if self.destinations.is_empty() {
                Ok(())
            } else {
                Err(invalid(
                    "ceramic heat deposition lost its retained physical step",
                ))
            };
        };
        if steps.len() != self.destinations.len() {
            return Err(invalid(
                "every intrinsic cell loss requires its original receiving cell",
            ));
        }
        for route in &self.destinations {
            checkpoint()?;
            let step = &steps[route.piezoelectric_cell];
            let heat = step.mechanical_heat_j + step.dielectric_heat_j;
            let power = heat / dt;
            let sum = powers.entry(route.thermal_cell).or_insert(0.0);
            if [
                step.mechanical_heat_j,
                step.dielectric_heat_j,
                heat,
                power,
                *sum,
            ]
            .iter()
            .any(|x| !x.is_finite() || *x < 0.0)
                || heat != 0.0 && power == 0.0
            {
                return Err(DepositionError::Arithmetic.into());
            }
            *sum += power;
            if !sum.is_finite() {
                return Err(DepositionError::Arithmetic.into());
            }
        }
        checkpoint()?;
        Ok(())
    }

    pub(super) fn destinations(&self) -> &[ThermalCeramicCell] {
        &self.destinations
    }
}
