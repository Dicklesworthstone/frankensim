//! Two-component conservative rezoning in retained-mass coordinate.
//!
//! Reconstruct char fraction and sensible enthalpy together. The limiter keeps
//! both component masses nonnegative and the reconstructed mixture inside its
//! common caloric domain. Formation references travel with component masses,
//! so changing a component's energy datum cannot create remap heating.

use super::*;
use crate::enthalpy::pyrolysis::{Arrhenius, PyrolysisCell, solve_cell};

impl RecedingSlab<'_> {
    pub(super) fn rezone_char(
        &mut self,
        controls: SlabSolve,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Remap, SlabError> {
        let n = self.ids.len();
        let char = self.char.expect("admitted charring geometry");
        let model = self.solid.materials[&self.solid.cells()[&self.ids[0]].material_id];
        let h_reference = model.enthalpies(model.temperature_range_k[0])?;
        let h_max = model.enthalpies(model.temperature_range_k[1])?;
        let relative_max = [h_max[0] - h_reference[0], h_max[1] - h_reference[1]];
        let char_datum = finite(h_reference[1] - h_reference[0])?;
        let mut edges = Vec::with_capacity(n + 1);
        let mut fractions = Vec::with_capacity(n);
        let mut relative_h = Vec::with_capacity(n);
        let mut before = [
            Sum::default(),
            Sum::default(),
            Sum::default(),
            Sum::default(),
        ];
        edges.push(0.0);
        for &id in &self.ids {
            checkpoint()?;
            let cell = &self.solid.cells()[&id];
            let inventory = self.solid.inventories[&id];
            let fraction = finite(inventory.char_mass_kg / cell.mass_kg)?;
            let h = finite((cell.specific_h_j_kg - h_reference[0]) - fraction * char_datum)?;
            before[0].add(cell.mass_kg)?;
            before[1].add(inventory.virgin_mass_kg)?;
            before[2].add(inventory.char_mass_kg)?;
            before[3].add(finite(cell.mass_kg * h)?)?;
            if before[0].value() <= *edges.last().unwrap() || !(0.0..=1.0).contains(&fraction) {
                return Err(DepositionError::Arithmetic.into());
            }
            edges.push(before[0].value());
            fractions.push(fraction);
            relative_h.push(h);
        }
        let mass_before = before[0].value();
        let mut fraction_slopes = limited_slopes(&edges, &fractions)?;
        let mut energy_slopes = limited_slopes(&edges, &relative_h)?;
        // The admissible h-Y set has affine lower and upper caloric bounds.
        // Scaling the two slopes together preserves each donor mean exactly.
        for i in 0..n {
            checkpoint()?;
            let mut factor = 1.0_f64;
            let width = edges[i + 1] - edges[i];
            let upper = relative_max[0] + fractions[i] * (relative_max[1] - relative_max[0]);
            for offset in [-0.5 * width, 0.5 * width] {
                let dy = fraction_slopes[i] * offset;
                let dh = energy_slopes[i] * offset;
                limit_nonnegative(fractions[i], dy, &mut factor);
                limit_nonnegative(1.0 - fractions[i], -dy, &mut factor);
                limit_nonnegative(relative_h[i], dh, &mut factor);
                limit_nonnegative(
                    upper - relative_h[i],
                    dy * (relative_max[1] - relative_max[0]) - dh,
                    &mut factor,
                );
            }
            if factor < 1.0 {
                factor *= 1.0 - 32.0 * f64::EPSILON;
            }
            fraction_slopes[i] *= factor;
            energy_slopes[i] *= factor;
        }

        let mut new_mass = Vec::with_capacity(n);
        let mut new_h = Vec::with_capacity(n);
        let mut inventories = Vec::with_capacity(n);
        let mut volumes = Vec::with_capacity(n);
        let mut after = [
            Sum::default(),
            Sum::default(),
            Sum::default(),
            Sum::default(),
        ];
        let mut volume_total = Sum::default();
        let mut redistribution = Sum::default();
        let mut temperature_evaluations = 0usize;
        let mut source = 0usize;
        let mut caloric_model = model;
        caloric_model.kinetics = Arrhenius {
            prefactor_s_inv: 0.0,
            activation_energy_j_mol: 0.0,
        };
        let caloric_controls = PyrolysisControls {
            cell_energy_tolerance_j: controls.energy_tolerance_j / (8.0 * n as f64),
            ..Default::default()
        };
        for j in 0..n {
            checkpoint()?;
            let left = mass_before * (j as f64 / n as f64);
            let right = if j + 1 == n {
                mass_before
            } else {
                mass_before * ((j + 1) as f64 / n as f64)
            };
            if right <= left {
                return Err(DepositionError::Arithmetic.into());
            }
            while source + 1 < n && edges[source + 1] <= left {
                source += 1;
            }
            let first = source;
            let mut target = [Sum::default(), Sum::default(), Sum::default()];
            let mut last = source;
            while source < n && edges[source] < right {
                let a = left.max(edges[source]);
                let b = right.min(edges[source + 1]);
                if b > a {
                    let overlap = b - a;
                    let offset = 0.5 * ((a - edges[source]) + (b - edges[source]))
                        - 0.5 * (edges[source + 1] - edges[source]);
                    target[0].add(overlap)?;
                    target[1].add(finite(
                        overlap * (fractions[source] + fraction_slopes[source] * offset),
                    )?)?;
                    target[2].add(finite(
                        overlap * (relative_h[source] + energy_slopes[source] * offset),
                    )?)?;
                    last = source;
                }
                if edges[source + 1] >= right {
                    break;
                }
                source += 1;
            }
            let char_mass = target[1].value();
            let virgin_mass = finite(target[0].value() - char_mass)?;
            let mass = finite(virgin_mass + char_mass)?;
            if !positive(mass) || !nonnegative(virgin_mass) || !nonnegative(char_mass) {
                return Err(DepositionError::Arithmetic.into());
            }
            // The ORIGINAL scalar energy solver inverts the mixed calorics.
            // Reaction is disabled; no volatile transfer is generated here.
            let recovered = solve_cell(
                PyrolysisCell {
                    virgin_mass_kg: virgin_mass,
                    char_mass_kg: char_mass,
                    temperature_k: model.temperature_range_k[0],
                },
                caloric_model,
                1.0,
                target[2].value(),
                caloric_controls,
                checkpoint,
            )?;
            temperature_evaluations += recovered.evaluations;
            let h = recovered.specific_h_j_kg;
            let fraction = char_mass / mass;
            let relative = finite((h - h_reference[0]) - fraction * char_datum)?;
            let volume = component_volume(recovered.inventory, self.geometry, Some(char))?;
            after[0].add(mass)?;
            after[1].add(virgin_mass)?;
            after[2].add(char_mass)?;
            after[3].add(finite(mass * relative)?)?;
            volume_total.add(volume)?;
            for i in first..=last {
                let overlap = right.min(edges[i + 1]) - left.max(edges[i]);
                if overlap > 0.0 {
                    redistribution.add(finite(overlap * (relative_h[i] - relative).abs())?)?;
                }
            }
            new_mass.push(mass);
            new_h.push(h);
            inventories.push(recovered.inventory);
            volumes.push(volume);
        }
        let mass_residual = finite(after[0].value() - before[0].value())?;
        let virgin_residual = finite(after[1].value() - before[1].value())?;
        let char_residual = finite(after[2].value() - before[2].value())?;
        let energy_residual = finite(
            after[3].value() - before[3].value()
                + mass_residual * h_reference[0]
                + char_residual * char_datum,
        )?;
        if mass_residual.abs() > controls.mass_tolerance_kg
            || virgin_residual.abs() + char_residual.abs() > controls.mass_tolerance_kg
            || energy_residual.abs() > controls.energy_tolerance_j
        {
            return Err(DepositionError::Budget("conservative char slab remap discrepancy").into());
        }

        let area = self.geometry.area_m2;
        let density_bound = self.geometry.density_kg_m3.max(char.density_kg_m3);
        let mut prefix = Sum::default();
        let mut coordinate = finite(self.back_m - volume_total.value() / area)?;
        let mut geometry_mass = 0.0_f64;
        for (i, &id) in self.ids.iter().enumerate() {
            checkpoint()?;
            prefix.add(volumes[i])?;
            let end = if i + 1 == n {
                self.back_m
            } else {
                finite(self.back_m - (volume_total.value() - prefix.value()) / area)?
            };
            if end <= coordinate {
                return Err(DepositionError::Arithmetic.into());
            }
            geometry_mass = geometry_mass
                .max(finite((area * (end - coordinate) - volumes[i]) * density_bound)?.abs());
            let cell = self.solid.thermal.cells.get_mut(&id).unwrap();
            cell.mass_kg = new_mass[i];
            cell.specific_h_j_kg = new_h[i];
            cell.geometry.start_m[0] = coordinate;
            cell.geometry.end_m[0] = end;
            cell.geometry.volume_m3 = volumes[i];
            self.solid.inventories.insert(id, inventories[i]);
            coordinate = end;
        }
        if geometry_mass > controls.mass_tolerance_kg {
            return Err(DepositionError::Budget(
                "rounded char slab geometry exceeds mass-equivalent budget",
            )
            .into());
        }
        let mut conductances = BTreeMap::new();
        for (i, &link_id) in self.link_ids.iter().enumerate() {
            checkpoint()?;
            let pair = [self.ids[i], self.ids[i + 1]];
            let halves = pair.map(|id| {
                let g = self.solid.thermal.cells[&id].geometry;
                0.5 * (g.end_m[0] - g.start_m[0])
            });
            conductances.insert(
                link_id,
                series_conductance_w_k(
                    area,
                    halves,
                    [
                        self.cell_conductivity_w_m_k(pair[0])?,
                        self.cell_conductivity_w_m_k(pair[1])?,
                    ],
                )?,
            );
        }
        self.solid
            .thermal
            .set_step_conductances(&conductances, checkpoint)?;
        checkpoint()?;
        Ok(Remap {
            mass_kg: mass_residual,
            virgin_mass_kg: virgin_residual,
            char_mass_kg: char_residual,
            energy_j: energy_residual,
            redistribution_j: redistribution.value(),
            geometry_mass_kg: geometry_mass,
            temperature_evaluations,
        })
    }
}

fn limit_nonnegative(mean: f64, deviation: f64, factor: &mut f64) {
    if mean + deviation < 0.0 && deviation < 0.0 {
        // A source on the exact caloric boundary must reconstruct constantly
        // in that constrained direction; do not alter its admitted mean.
        *factor = factor.min((mean / -deviation).max(0.0));
    }
}
