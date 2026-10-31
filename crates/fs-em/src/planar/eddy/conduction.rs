//! Shared isolated-conductor P1 projection and actual spatial Ohmic heating.
//!
//! The field owners admit material sources, topology and geometry before
//! constructing these coefficients. Both midpoint and nonlinear induction use
//! this exact algebra; no averaged-current or fitted magnetic loss is added.

pub(super) struct ConductionCell {
    pub vertices: [usize; 3],
    pub group: usize,
    pub sigma_volume: f64,
    pub conductivity_s_m: f64,
    pub volume_m3: f64,
}
pub(super) struct Conduction {
    pub cells: Vec<ConductionCell>,
    pub regions: Vec<u32>,
    pub weights: Vec<f64>,
    pub anchors: Vec<usize>,
}
pub(super) struct CellHeat {
    pub electric_field_v_m: [f64; 3],
    pub current_density_a_m2: [f64; 3],
    pub volume_m3: f64,
    pub power_w: f64,
    pub heat_j: f64,
}
pub(super) struct JouleHeat {
    pub cells: Vec<CellHeat>,
    pub heat_j: f64,
    pub region_currents_a: Vec<(u32, f64)>,
}
#[derive(Debug, PartialEq)]
pub(super) enum ConductionError {
    Cancelled,
    Invalid(&'static str),
    Balance(&'static str),
}

fn finite(value: f64) -> Result<f64, ConductionError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(ConductionError::Invalid(
            "unrepresentable induction arithmetic",
        ))
    }
}

impl Conduction {
    // An additional cell-local P1 field need not be continuous across magnetic
    // material interfaces. Subtract its own regional anchor before taking the
    // mean, just as for the original nodal potential. This also annihilates a
    // spatially constant motional EMF without cancellation of large values.
    fn offset_anchors(&self, offsets: &[[f64; 3]], check: &dyn Fn() -> bool) -> Option<Vec<f64>> {
        let mut anchors = vec![None; self.regions.len()];
        for (cell, values) in self.cells.iter().zip(offsets) {
            if check() {
                return None;
            }
            anchors[cell.group].get_or_insert(values[0]);
        }
        Some(anchors.into_iter().map(|v| v.unwrap_or(0.0)).collect())
    }

    /// Exact conductivity-weighted means of a cell-local P1 field. Retain a
    /// regional anchor so a large constant motional field does not lose its
    /// small spatial correction through cancellation in a raw weighted sum.
    pub fn mean_offsets(
        &self,
        offsets: &[[f64; 3]],
        check: &dyn Fn() -> bool,
    ) -> Result<Vec<f64>, ConductionError> {
        if offsets.len() != self.cells.len() {
            return Err(ConductionError::Invalid(
                "complete regional offset field required",
            ));
        }
        let anchors = self
            .offset_anchors(offsets, check)
            .ok_or(ConductionError::Cancelled)?;
        let mut means = vec![0.0; self.regions.len()];
        for (cell, offset) in self.cells.iter().zip(offsets) {
            if check() {
                return Err(ConductionError::Cancelled);
            }
            let delta = offset
                .iter()
                .map(|&value| (value - anchors[cell.group]) / 3.0)
                .sum::<f64>();
            means[cell.group] = finite(cell.sigma_volume.mul_add(delta, means[cell.group]))?;
        }
        for ((mean, &weight), &anchor) in means.iter_mut().zip(&self.weights).zip(&anchors) {
            *mean = finite(*mean / weight + anchor)?;
        }
        Ok(means)
    }

    fn means_with_offsets(
        &self,
        x: &[f64],
        offsets: &[[f64; 3]],
        offset_anchors: &[f64],
        check: &dyn Fn() -> bool,
    ) -> Option<Vec<f64>> {
        let mut means = vec![0.0; self.regions.len()];
        for (cell, offsets) in self.cells.iter().zip(offsets) {
            if check() {
                return None;
            }
            let anchor = x[self.anchors[cell.group]];
            means[cell.group] += cell.sigma_volume
                * (0..3)
                    .map(|i| {
                        ((x[cell.vertices[i]] - anchor) + (offsets[i] - offset_anchors[cell.group]))
                            / 3.0
                    })
                    .sum::<f64>();
        }
        for (mean, weight) in means.iter_mut().zip(&self.weights) {
            *mean /= weight;
        }
        Some(means)
    }

    /// Apply the same isolated-current projection to a nodal P1 field plus
    /// cell-local P1 offsets. Moving induction uses the latter for v.grad(A).
    pub fn apply_with_offsets(
        &self,
        x: &[f64],
        offsets: &[[f64; 3]],
        y: &mut [f64],
        check: &dyn Fn() -> bool,
    ) -> bool {
        assert_eq!(offsets.len(), self.cells.len());
        y.fill(0.0);
        let Some(offset_anchors) = self.offset_anchors(offsets, check) else {
            return false;
        };
        let Some(means) = self.means_with_offsets(x, offsets, &offset_anchors, check) else {
            return false;
        };
        for (cell, offsets) in self.cells.iter().zip(offsets) {
            if check() {
                return false;
            }
            let anchor = x[self.anchors[cell.group]];
            let d: [f64; 3] = std::array::from_fn(|i| {
                ((x[cell.vertices[i]] - anchor) + (offsets[i] - offset_anchors[cell.group]))
                    - means[cell.group]
            });
            let sum = d.iter().sum::<f64>();
            for i in 0..3 {
                y[cell.vertices[i]] += cell.sigma_volume / 12.0 * (d[i] + sum);
            }
        }
        true
    }

    /// Exactly integrate the conductor-frame sigma E² for a nodal rate plus
    /// cell-local v.grad(A), with the same projected field as the weak solve.
    pub fn heat_with_offsets(
        &self,
        rate: &[f64],
        offsets: &[[f64; 3]],
        dt: f64,
        thickness_m: f64,
        current_tolerance_a: f64,
        drives: Option<(&[f64], &[f64])>,
        check: &dyn Fn() -> bool,
    ) -> Result<JouleHeat, ConductionError> {
        assert_eq!(offsets.len(), self.cells.len());
        if let Some((electric, targets)) = drives {
            if electric.len() != self.regions.len()
                || targets.len() != self.regions.len()
                || electric.iter().chain(targets).any(|v| !v.is_finite())
            {
                return Err(ConductionError::Invalid(
                    "complete finite regional electric drives and current targets required",
                ));
            }
        }
        let offset_anchors = self
            .offset_anchors(offsets, check)
            .ok_or(ConductionError::Cancelled)?;
        let means = self
            .means_with_offsets(rate, offsets, &offset_anchors, check)
            .ok_or(ConductionError::Cancelled)?;
        self.integrate_heat(
            dt,
            thickness_m,
            current_tolerance_a,
            drives.map(|(_, targets)| targets),
            check,
            |index, cell| {
                let anchor = rate[self.anchors[cell.group]];
                std::array::from_fn(|i| {
                    drives.map_or(0.0, |(electric, _)| electric[cell.group])
                        - (((rate[cell.vertices[i]] - anchor)
                            + (offsets[index][i] - offset_anchors[cell.group]))
                            - means[cell.group])
                })
            },
        )
    }

    // Subtract a regional coordinate BEFORE summation; constant gauge values
    // cannot manufacture axial field through cancellation of large offsets.
    pub fn means(&self, x: &[f64], check: &dyn Fn() -> bool) -> Option<Vec<f64>> {
        let mut means = vec![0.0; self.regions.len()];
        for cell in &self.cells {
            if check() {
                return None;
            }
            let anchor = x[self.anchors[cell.group]];
            means[cell.group] += cell.sigma_volume
                * cell
                    .vertices
                    .iter()
                    .map(|&v| (x[v] - anchor) / 3.0)
                    .sum::<f64>();
        }
        for (v, &weight) in means.iter_mut().zip(&self.weights) {
            *v /= weight;
        }
        Some(means)
    }

    /// Consistent P1 mass minus the regional conductivity-weighted mean.
    pub fn apply(&self, x: &[f64], y: &mut [f64], check: &dyn Fn() -> bool) -> bool {
        y.fill(0.0);
        let Some(means) = self.means(x, check) else {
            return false;
        };
        for cell in &self.cells {
            if check() {
                return false;
            }
            let anchor = x[self.anchors[cell.group]];
            let d = cell.vertices.map(|v| (x[v] - anchor) - means[cell.group]);
            let sum = d.iter().sum::<f64>();
            for i in 0..3 {
                y[cell.vertices[i]] += cell.sigma_volume / 12.0 * (d[i] + sum);
            }
        }
        true
    }

    /// Exactly integrate sigma E^2, with E from the SAME projected increment.
    pub fn heat(
        &self,
        increments: &[f64],
        dt: f64,
        thickness_m: f64,
        current_tolerance_a: f64,
        check: &dyn Fn() -> bool,
    ) -> Result<JouleHeat, ConductionError> {
        let means = self
            .means(increments, check)
            .ok_or(ConductionError::Cancelled)?;
        self.integrate_heat(
            dt,
            thickness_m,
            current_tolerance_a,
            None,
            check,
            |_, cell| {
                let anchor = increments[self.anchors[cell.group]];
                cell.vertices
                    .map(|v| -((increments[v] - anchor) - means[cell.group]) / dt)
            },
        )
    }

    /// Add a physical uniform axial electric drive AFTER the zero-current
    /// projection. Solid winding sections prescribe their integrated current;
    /// unbound regions pass zero drive/target and retain floating behavior.
    /// The same exact P1 quadrature includes all local drive/eddy cross terms.
    #[allow(clippy::too_many_arguments)]
    pub fn heat_with_region_drives(
        &self,
        increments: &[f64],
        dt: f64,
        thickness_m: f64,
        current_tolerance_a: f64,
        electric_drives_v_m: &[f64],
        target_currents_a: &[f64],
        check: &dyn Fn() -> bool,
    ) -> Result<JouleHeat, ConductionError> {
        if electric_drives_v_m.len() != self.regions.len()
            || target_currents_a.len() != self.regions.len()
            || electric_drives_v_m
                .iter()
                .chain(target_currents_a)
                .any(|v| !v.is_finite())
        {
            return Err(ConductionError::Invalid(
                "complete finite regional electric drives and current targets required",
            ));
        }
        let means = self
            .means(increments, check)
            .ok_or(ConductionError::Cancelled)?;
        self.integrate_heat(
            dt,
            thickness_m,
            current_tolerance_a,
            Some(target_currents_a),
            check,
            |_, cell| {
                let anchor = increments[self.anchors[cell.group]];
                cell.vertices.map(|v| {
                    electric_drives_v_m[cell.group]
                        - ((increments[v] - anchor) - means[cell.group]) / dt
                })
            },
        )
    }

    // One quadrature/accounting owner serves stationary, nonlinear and moving
    // induction. Only the constitutive-frame electric field differs.
    fn integrate_heat(
        &self,
        dt: f64,
        thickness_m: f64,
        current_tolerance_a: f64,
        target_currents_a: Option<&[f64]>,
        check: &dyn Fn() -> bool,
        electric: impl Fn(usize, &ConductionCell) -> [f64; 3],
    ) -> Result<JouleHeat, ConductionError> {
        let mut cells = Vec::with_capacity(self.cells.len());
        let mut net_current = vec![0.0; self.regions.len()];
        let mut heat_j = 0.0;
        for (index, cell) in self.cells.iter().enumerate() {
            if check() {
                return Err(ConductionError::Cancelled);
            }
            let electric_field_v_m = electric(index, cell);
            let current_density_a_m2 = electric_field_v_m.map(|e| cell.conductivity_s_m * e);
            let e_sum = electric_field_v_m.iter().sum::<f64>();
            let power_w = finite(
                cell.sigma_volume / 12.0
                    * (e_sum * e_sum + electric_field_v_m.iter().map(|e| e * e).sum::<f64>()),
            )?;
            let heat = finite(power_w * dt)?;
            if current_density_a_m2.iter().any(|v| !v.is_finite())
                || power_w < 0.0
                || electric_field_v_m != [0.0; 3] && power_w == 0.0
                || power_w != 0.0 && heat == 0.0
            {
                return Err(ConductionError::Balance(
                    "represented eddy heat is not finite and nonnegative",
                ));
            }
            net_current[cell.group] = finite(
                net_current[cell.group]
                    + current_density_a_m2.iter().sum::<f64>() / 3.0 * cell.volume_m3 / thickness_m,
            )?;
            heat_j = finite(heat_j + heat)?;
            cells.push(CellHeat {
                electric_field_v_m,
                current_density_a_m2,
                volume_m3: cell.volume_m3,
                power_w,
                heat_j: heat,
            });
        }
        if net_current.iter().enumerate().any(|(group, &i)| {
            (i - target_currents_a.map_or(0.0, |targets| targets[group])).abs()
                > current_tolerance_a
        }) {
            return Err(ConductionError::Balance(
                "conductor has unresolved prescribed axial-current balance",
            ));
        }
        Ok(JouleHeat {
            cells,
            heat_j,
            region_currents_a: self.regions.iter().copied().zip(net_current).collect(),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn conductor() -> Conduction {
        let sigma = 18.0;
        let volume = 0.3 / 4.0;
        Conduction {
            cells: [[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]]
                .map(|vertices| ConductionCell {
                    vertices,
                    group: 0,
                    sigma_volume: sigma * volume,
                    conductivity_s_m: sigma,
                    volume_m3: volume,
                })
                .into(),
            regions: vec![7],
            weights: vec![sigma * 0.3],
            anchors: vec![0],
        }
    }
    fn near(a: f64, b: f64) {
        assert!((a - b).abs() < 2e-12 * (1.0 + b.abs()), "{a} != {b}");
    }

    #[test]
    fn consistent_mass_projection_heat_and_current_share_the_same_p1_mode() {
        let conductor = conductor();
        let increment = [0.0, 0.0, 0.0, 0.0, 3.0];
        let mut applied = [0.0; 5];
        assert!(conductor.apply(&increment, &mut applied, &|| false));
        // Integral N4=1/3, integral N4^2=1/6 on the unit square.
        // C44 = sigma*depth*(1/6 - 1/9) = 0.3, not the unprojected 0.9.
        near(applied[4], 0.3 * 3.0);
        near(applied.iter().sum(), 0.0);
        let heat = conductor
            .heat(&increment, 0.2, 0.3, 1e-10, &|| false)
            .unwrap();
        near(heat.heat_j, 0.3 * 3.0_f64.powi(2) / 0.2);
        near(heat.region_currents_a[0].1, 0.0);
        for (cell, values) in conductor.cells.iter().zip(&heat.cells) {
            for (local, &vertex) in cell.vertices.iter().enumerate() {
                near(
                    values.electric_field_v_m[local],
                    if vertex == 4 { -10.0 } else { 5.0 },
                );
                near(
                    values.current_density_a_m2[local],
                    18.0 * values.electric_field_v_m[local],
                );
            }
            near(values.heat_j, heat.heat_j / 4.0);
            near(values.power_w, values.heat_j / 0.2);
            near(values.volume_m3, 0.075);
        }
        let bilinear = increment
            .iter()
            .zip(applied)
            .map(|(x, y)| x * y)
            .sum::<f64>();
        near(bilinear / 0.2, heat.heat_j);
    }

    #[test]
    fn gauge_current_reversal_and_cancellation_do_not_change_the_projection() {
        let conductor = conductor();
        let zero = conductor
            .heat(&[1e8; 5], 0.2, 0.3, 1e-10, &|| false)
            .unwrap();
        near(zero.heat_j, 0.0);
        let a = conductor
            .heat(&[0.0, 0.0, 0.0, 0.0, 3.0], 0.2, 0.3, 1e-10, &|| false)
            .unwrap();
        let b = conductor
            .heat(&[1e8, 1e8, 1e8, 1e8, 1e8 - 3.0], 0.2, 0.3, 1e-10, &|| false)
            .unwrap();
        near(a.heat_j, b.heat_j);
        for (a, b) in a.cells.iter().zip(b.cells) {
            for axis in 0..3 {
                near(a.current_density_a_m2[axis], -b.current_density_a_m2[axis]);
            }
        }
        assert!(matches!(
            conductor.heat(&[0.0; 5], 0.2, 0.3, 1e-10, &|| true),
            Err(ConductionError::Cancelled)
        ));
    }
}
