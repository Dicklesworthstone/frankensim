//! Original surface-temperature contact laws inside the native tensor solve.
//!
//! The scalar energy owner and this tensor owner share SourceFaceContact.
//! Only the spatial equations differ: tensor faces already retain both physical
//! traces, so their nonlinear contact response is assembled directly.

use std::collections::BTreeMap;

use fs_exec::CancelGate;

use super::{ChtError, Joint, Prepared, ThermalSetup, VoxelDomain, invalid, poll};
use crate::conjugate::energy::{SourceContactState, SourceContactTangent, SourceFaceContact};

mod local;
pub(super) use local::ContactHeat;

#[derive(Clone, Copy)]
pub(super) struct SourceJoint {
    pub index: usize,
    pub reversed: bool,
}

struct Binding<'a> {
    law: &'a SourceFaceContact,
    tangent: SourceContactTangent<'a>,
    face_from_a: usize,
    area_m2: f64,
    reference_resistance_m2_k_w: f64,
}

/// Immutable source preparation, absent for the frozen tensor entry.
/// A source is prepared once, while every physical trial refreshes the same
/// selected original claim and caches its complete two-trace tangent.
pub(super) struct PreparedSources<'a> {
    by_face: BTreeMap<(usize, usize), usize>,
    bindings: Vec<Binding<'a>>,
}

impl<'a> PreparedSources<'a> {
    pub fn new(
        domain: &VoxelDomain,
        setup: &'a ThermalSetup,
        gate: &CancelGate,
    ) -> Result<Self, ChtError> {
        poll(gate)?;
        let count = setup.source_face_contacts.len();
        if count > domain.cell_count().saturating_mul(3) {
            return Err(invalid(
                "source contact count exceeds the Cartesian face count",
            ));
        }
        let mut by_face = BTreeMap::new();
        let mut bindings = Vec::with_capacity(count);
        for (index, law) in setup.source_face_contacts.iter().enumerate() {
            if index.is_multiple_of(256) {
                poll(gate)?;
            }
            let [a, b] = law.cells();
            if a >= domain.cell_count() || b >= domain.cell_count() || a == b {
                return Err(invalid("source contact needs two distinct in-range cells"));
            }
            let face_from_a = (0..6)
                .find(|&face| domain.neighbor(a, face) == Some(b))
                .ok_or_else(|| invalid("source contact cells must be face-adjacent"))?;
            if by_face.insert((a.min(b), a.max(b)), index).is_some() {
                return Err(invalid(
                    "duplicate source contact, including reversed surface order",
                ));
            }
            let area_m2 = domain.face_area(a, face_from_a / 2);
            let reference_resistance_m2_k_w = law.reference_resistance()?;
            let tangent = law.prepare_temperature_tangent()?;
            poll(gate)?;
            bindings.push(Binding {
                law,
                tangent,
                face_from_a,
                area_m2,
                reference_resistance_m2_k_w,
            });
        }
        Ok(Self { by_face, bindings })
    }

    pub fn cell_pairs(&self) -> impl Iterator<Item = (usize, usize)> + '_ {
        self.by_face.keys().copied()
    }

    pub fn face(&self, first: usize, second: usize) -> Option<(SourceJoint, f64)> {
        self.by_face
            .get(&(first.min(second), first.max(second)))
            .map(|&index| {
                let binding = &self.bindings[index];
                (
                    SourceJoint {
                        index,
                        reversed: binding.law.cells()[0] != first,
                    },
                    binding.reference_resistance_m2_k_w,
                )
            })
    }

    /// Use the declared supported interface state, not an inferred bulk mean.
    /// The original bulk conductivity must still admit every cell-face interval.
    pub fn seed_traces(
        &self,
        data: &Prepared,
        trace: &mut [f64],
        gate: &CancelGate,
    ) -> Result<(), ChtError> {
        for (index, binding) in self.bindings.iter().enumerate() {
            if index.is_multiple_of(256) {
                poll(gate)?;
            }
            let [a, b] = binding.law.cells();
            let f = binding.face_from_a;
            let ids = [data.face_ids[a][f], data.face_ids[b][f ^ 1]];
            if ids[0] == ids[1] {
                return Err(invalid(
                    "a finite source contact must retain both physical traces",
                ));
            }
            let temperatures = binding.law.reference_surfaces()?;
            for side in 0..2 {
                trace[ids[side]] = temperatures[side];
            }
        }
        poll(gate)
    }

    /// One original source query per face and physical Newton trial. Numerical
    /// heat/tangent values are cached for the shared residual, audit and solve;
    /// receipts report actual bulk face heat in the original A/B source order.
    pub fn evaluate(
        &self,
        data: &Prepared,
        trace: &[f64],
        outward_heat_w: &[[f64; 6]],
        gate: &CancelGate,
    ) -> Result<(Vec<ContactHeat>, Vec<SourceContactState>), ChtError> {
        let mut heat = Vec::with_capacity(self.bindings.len());
        let mut states = Vec::with_capacity(self.bindings.len());
        for (index, binding) in self.bindings.iter().enumerate() {
            if index.is_multiple_of(256) {
                poll(gate)?;
            }
            let cells = binding.law.cells();
            let f = binding.face_from_a;
            let temperatures = [
                trace[data.face_ids[cells[0]][f]],
                trace[data.face_ids[cells[1]][f ^ 1]],
            ];
            let (resistance, partials, source) = binding
                .law
                .at_with_temperature_tangent(temperatures, &binding.tangent)?;
            poll(gate)?;
            let response =
                local::original_resistance(binding.area_m2, temperatures, resistance, partials)
                    .ok_or_else(|| {
                        invalid("original contact heat or full trace tangent is unrepresentable")
                    })?;
            let physical_heat = outward_heat_w[cells[0]][f];
            let flux = physical_heat / binding.area_m2;
            let residual = (-resistance).mul_add(flux, temperatures[0] - temperatures[1]);
            if !physical_heat.is_finite() || !flux.is_finite() || !residual.is_finite() {
                return Err(invalid(
                    "original contact surface balance is unrepresentable",
                ));
            }
            heat.push(response);
            states.push(SourceContactState {
                cells,
                surface_temperature_k: temperatures,
                heat_rate_first_to_second_w: physical_heat,
                resistance_m2_k_w: resistance,
                applied_resistance_m2_k_w: resistance,
                resistance_relative_residual: 0.0,
                contact_temperature_residual_k: residual,
                source,
            });
        }
        poll(gate)?;
        Ok((heat, states))
    }
}

/// Convert original A/B orientation to the prepared spatial face orientation.
/// Fixed contacts keep their original arithmetic and require no source query.
pub(super) fn joint_heat(
    joint: &Joint,
    temperatures: [f64; 2],
    source_heat: &[ContactHeat],
) -> Result<ContactHeat, ChtError> {
    if let Some(source) = joint.source {
        let response = *source_heat
            .get(source.index)
            .ok_or_else(|| invalid("source contact lost its original trial response"))?;
        Ok(if source.reversed {
            response.reversed()
        } else {
            response
        })
    } else {
        local::linear(joint.conductance, temperatures)
            .ok_or_else(|| invalid("fixed contact heat is unrepresentable"))
    }
}
