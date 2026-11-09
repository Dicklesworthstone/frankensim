//! Component energy rows for an already admitted closed thermal tangent.
//!
//! This numerical core does not invent a bath or a reference temperature. With
//! exp(-i omega t), summing the original equations gives -i omega sum(E)=sum(P).
//! At DC the caller instead supplies the conserved energy increment, normally
//! zero, and each component's heat drive must sum to zero. The native harmonic
//! owner retains its original LU and checks every original power equation.

use fs_math::{c64::C64, dd::Dd};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IsolatedThermalImpedanceControls {
    /// Independent absolute residual of every original cell equation, W.
    pub cell_power_tolerance_w: f64,
    /// Absolute complex net-power compatibility/conservation tolerance, W.
    pub component_power_tolerance_w: f64,
    /// Absolute complex extensive stored-energy residual tolerance, J.
    /// This remains absolute even near a physical 1/omega storage pole.
    pub component_energy_tolerance_j: f64,
}
impl Default for IsolatedThermalImpedanceControls {
    fn default() -> Self {
        Self {
            cell_power_tolerance_w: 1e-9,
            component_power_tolerance_w: 1e-9,
            component_energy_tolerance_j: 1e-9,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum ThermalEnergyConstraintError<E> {
    Invalid(&'static str),
    Checkpoint(E),
    Nonconservative { column: usize, residual_per_s: f64 },
    Underdetermined { component: usize, latent_cells: usize },
    IncompatibleHeat { component: usize, net_power_w: C64 },
    EnergyResidual { component: usize, residual_j: C64 },
    PowerResidual { component: usize, residual_w: C64 },
}
type Result<T, E> = core::result::Result<T, ThermalEnergyConstraintError<E>>;

#[derive(Clone, Debug, PartialEq)]
pub struct ThermalComponentResponse {
    /// Actual compensated sum of the original heat phasors, W.
    pub net_heat_w: C64,
    /// Actual compensated sum of extensive original stored-energy changes, J.
    pub energy_j: C64,
    /// sum(E)-prescribed energy at DC, or sum(E)-i*sum(P)/omega otherwise.
    pub energy_residual_j: C64,
    /// -i*omega*sum(E)-sum(P), independently checked in watts.
    pub power_residual_w: C64,
}

/// Numerical inventory constraints; physical closed-graph admission belongs to
/// the immutable original isolated-equilibrium owner, not this helper.
pub struct ConservedThermalEnergy {
    components: Vec<Vec<usize>>,
    latent_cells: Vec<usize>,
    cell_count: usize,
    controls: IsolatedThermalImpedanceControls,
}

pub struct ConservedThermalMatrix {
    pub values: Vec<C64>,
    pub row_scales: Vec<f64>,
}

pub struct ConservedThermalDrive {
    pub rhs: Vec<C64>,
    omega_rad_s: f64,
    // Keep the original-coordinate sums through the final admission. Rounding
    // a large component energy into its RHS must not hide a small residual.
    net_heat: Vec<[Dd; 2]>,
    target_energy: Vec<[Dd; 2]>,
}

fn real<E>(value: Dd) -> Result<f64, E> {
    if !value.hi.is_finite()
        || !value.lo.is_finite()
        || (value.hi == 0.0 && value.lo != 0.0)
    {
        return Err(ThermalEnergyConstraintError::Invalid(
            "component energy arithmetic is not representable",
        ));
    }
    Ok(value.to_f64())
}
fn complex<E>(value: [Dd; 2]) -> Result<C64, E> {
    Ok(C64::new(real(value[0])?, real(value[1])?))
}
fn input<E>(value: C64) -> Result<[Dd; 2], E> {
    if !value.re.is_finite() || !value.im.is_finite() {
        return Err(ThermalEnergyConstraintError::Invalid(
            "finite thermal phasors required",
        ));
    }
    Ok([Dd::from_f64(value.re), Dd::from_f64(value.im)])
}
fn magnitude<E>(value: [Dd; 2]) -> Result<f64, E> {
    let result = complex(value)?.abs();
    if !result.is_finite() {
        return Err(ThermalEnergyConstraintError::Invalid(
            "component residual magnitude is not representable",
        ));
    }
    Ok(result)
}
fn sum<E>(values: &[C64], cells: &[usize]) -> Result<[Dd; 2], E> {
    let mut total = [Dd::ZERO; 2];
    for &cell in cells {
        let value = input(values[cell])?;
        for k in 0..2 {
            total[k] = total[k] + value[k];
            real(total[k])?;
        }
    }
    Ok(total)
}

impl ConservedThermalEnergy {
    pub fn new<E>(
        components: Vec<Vec<usize>>,
        temperature_per_energy: &[f64],
        tangent_per_s: &[f64],
        controls: IsolatedThermalImpedanceControls,
        checkpoint: &mut impl FnMut() -> core::result::Result<(), E>,
    ) -> Result<Self, E> {
        checkpoint().map_err(ThermalEnergyConstraintError::Checkpoint)?;
        let n = temperature_per_energy.len();
        if n == 0
            || n > 1024
            || tangent_per_s.len() != n * n
            || components.is_empty()
            || !controls.cell_power_tolerance_w.is_finite()
            || controls.cell_power_tolerance_w < 0.0
            || !controls.component_power_tolerance_w.is_finite()
            || controls.component_power_tolerance_w < 0.0
            || !controls.component_energy_tolerance_j.is_finite()
            || controls.component_energy_tolerance_j < 0.0
        {
            return Err(ThermalEnergyConstraintError::Invalid(
                "invalid closed thermal dimensions or tolerances",
            ));
        }
        let mut membership = vec![usize::MAX; n];
        let mut latent_cells = Vec::with_capacity(components.len());
        for (component, cells) in components.iter().enumerate() {
            checkpoint().map_err(ThermalEnergyConstraintError::Checkpoint)?;
            if cells.is_empty() {
                return Err(ThermalEnergyConstraintError::Invalid(
                    "empty thermal component",
                ));
            }
            let mut latent = 0;
            for &cell in cells {
                if cell >= n || membership[cell] != usize::MAX {
                    return Err(ThermalEnergyConstraintError::Invalid(
                        "thermal components must partition the original cells",
                    ));
                }
                let slope = temperature_per_energy[cell];
                if !slope.is_finite() || slope < 0.0 {
                    return Err(ThermalEnergyConstraintError::Invalid(
                        "invalid caloric tangent",
                    ));
                }
                latent += usize::from(slope == 0.0);
                membership[cell] = component;
            }
            latent_cells.push(latent);
        }
        if membership.contains(&usize::MAX) {
            return Err(ThermalEnergyConstraintError::Invalid(
                "thermal components omit original cells",
            ));
        }
        // The original isolated owner proves conducting topology. Here audit
        // its represented tangent, without projecting or repairing any entry.
        for column in 0..n {
            checkpoint().map_err(ThermalEnergyConstraintError::Checkpoint)?;
            let (mut total, mut scale) = (Dd::ZERO, Dd::ZERO);
            for row in 0..n {
                let entry = tangent_per_s[row * n + column];
                if !entry.is_finite()
                    || (membership[row] != membership[column] && entry != 0.0)
                {
                    return Err(ThermalEnergyConstraintError::Invalid(
                        "closed tangent crosses an original inventory or is nonfinite",
                    ));
                }
                total = total + Dd::from_f64(entry);
                scale = scale + Dd::from_f64(entry.abs());
            }
            let residual = real(total)?;
            // Covers original diagonal assembly roundoff; final unmodified
            // cell-power and component-energy equations remain authoritative.
            if residual.abs() > 8.0 * n as f64 * f64::EPSILON * real(scale)? {
                return Err(ThermalEnergyConstraintError::Nonconservative {
                    column,
                    residual_per_s: residual,
                });
            }
        }
        Ok(Self {
            components,
            latent_cells,
            cell_count: n,
            controls,
        })
    }

    /// Scale original rows and replace exactly one dependent power row per
    /// component by its extensive-energy sum. The caller uses the original LU.
    pub fn matrix<E>(
        &self,
        tangent_per_s: &[f64],
        omega_rad_s: f64,
        checkpoint: &mut impl FnMut() -> core::result::Result<(), E>,
    ) -> Result<ConservedThermalMatrix, E> {
        let n = self.cell_count;
        if !omega_rad_s.is_finite()
            || omega_rad_s < 0.0
            || tangent_per_s.len() != n * n
        {
            return Err(ThermalEnergyConstraintError::Invalid(
                "invalid thermal frequency or tangent",
            ));
        }
        let mut replaced = vec![false; n];
        let mut values = vec![C64::ZERO; n * n];
        for (component, cells) in self.components.iter().enumerate() {
            checkpoint().map_err(ThermalEnergyConstraintError::Checkpoint)?;
            if omega_rad_s == 0.0 && self.latent_cells[component] > 1 {
                return Err(ThermalEnergyConstraintError::Underdetermined {
                    component,
                    latent_cells: self.latent_cells[component],
                });
            }
            let row = cells[0];
            replaced[row] = true;
            for &column in cells {
                values[row * n + column] = C64::ONE;
            }
        }
        let mut row_scales = vec![1.0; n];
        for row in 0..n {
            checkpoint().map_err(ThermalEnergyConstraintError::Checkpoint)?;
            if replaced[row] {
                continue;
            }
            let original = &tangent_per_s[row * n..(row + 1) * n];
            let scale = original
                .iter()
                .fold(omega_rad_s, |scale, value| scale.max(value.abs()));
            if !scale.is_finite() || scale == 0.0 {
                return Err(ThermalEnergyConstraintError::Invalid(
                    "singular closed thermal tangent",
                ));
            }
            row_scales[row] = scale;
            for column in 0..n {
                values[row * n + column] = C64::new(
                    original[column] / scale,
                    if row == column { -omega_rad_s / scale } else { 0.0 },
                );
            }
        }
        Ok(ConservedThermalMatrix { values, row_scales })
    }

    /// None means unchanged component energy at DC. Explicit energies are
    /// admitted only at DC; harmonic energy follows the actual summed heat.
    pub fn drive<E>(
        &self,
        omega_rad_s: f64,
        heat_w: &[C64],
        dc_energy_j: Option<&[C64]>,
        checkpoint: &mut impl FnMut() -> core::result::Result<(), E>,
    ) -> Result<ConservedThermalDrive, E> {
        if !omega_rad_s.is_finite()
            || omega_rad_s < 0.0
            || heat_w.len() != self.cell_count
            || dc_energy_j.is_some_and(|energy| {
                omega_rad_s != 0.0 || energy.len() != self.components.len()
            })
        {
            return Err(ThermalEnergyConstraintError::Invalid(
                "invalid closed thermal drive dimensions",
            ));
        }
        let mut drive = ConservedThermalDrive {
            rhs: heat_w.to_vec(),
            omega_rad_s,
            net_heat: Vec::with_capacity(self.components.len()),
            target_energy: Vec::with_capacity(self.components.len()),
        };
        for (component, cells) in self.components.iter().enumerate() {
            checkpoint().map_err(ThermalEnergyConstraintError::Checkpoint)?;
            let heat = sum(heat_w, cells)?;
            let target = if omega_rad_s == 0.0 {
                if magnitude(heat)? > self.controls.component_power_tolerance_w {
                    return Err(ThermalEnergyConstraintError::IncompatibleHeat {
                        component,
                        net_power_w: complex(heat)?,
                    });
                }
                input(dc_energy_j.map_or(C64::ZERO, |energy| energy[component]))?
            } else {
                let omega = Dd::from_f64(omega_rad_s);
                let target = [-heat[1] / omega, heat[0] / omega];
                for k in 0..2 {
                    if real(target[k])? == 0.0 && heat[1 - k] != Dd::ZERO {
                        return Err(ThermalEnergyConstraintError::Invalid(
                            "harmonic component energy underflows",
                        ));
                    }
                }
                target
            };
            drive.rhs[cells[0]] = complex(target)?;
            drive.net_heat.push(heat);
            drive.target_energy.push(target);
        }
        Ok(drive)
    }

    /// Independent absolute joule and watt gates, after the native LU has
    /// already checked every original row, including those replaced above.
    pub fn admit<E>(
        &self,
        drive: &ConservedThermalDrive,
        energy_j: &[C64],
        checkpoint: &mut impl FnMut() -> core::result::Result<(), E>,
    ) -> Result<Vec<ThermalComponentResponse>, E> {
        if energy_j.len() != self.cell_count
            || drive.net_heat.len() != self.components.len()
        {
            return Err(ThermalEnergyConstraintError::Invalid(
                "invalid thermal response dimensions",
            ));
        }
        let mut reports = Vec::with_capacity(self.components.len());
        for (component, cells) in self.components.iter().enumerate() {
            checkpoint().map_err(ThermalEnergyConstraintError::Checkpoint)?;
            let energy = sum(energy_j, cells)?;
            let residual = [
                energy[0] - drive.target_energy[component][0],
                energy[1] - drive.target_energy[component][1],
            ];
            if magnitude(residual)? > self.controls.component_energy_tolerance_j {
                return Err(ThermalEnergyConstraintError::EnergyResidual {
                    component,
                    residual_j: complex(residual)?,
                });
            }
            let omega = Dd::from_f64(drive.omega_rad_s);
            let storage = [omega * energy[1], -omega * energy[0]];
            for k in 0..2 {
                if drive.omega_rad_s != 0.0
                    && energy[1 - k] != Dd::ZERO
                    && real(storage[k])? == 0.0
                {
                    return Err(ThermalEnergyConstraintError::Invalid(
                        "harmonic component power underflows",
                    ));
                }
            }
            let power = [
                storage[0] - drive.net_heat[component][0],
                storage[1] - drive.net_heat[component][1],
            ];
            if magnitude(power)? > self.controls.component_power_tolerance_w {
                return Err(ThermalEnergyConstraintError::PowerResidual {
                    component,
                    residual_w: complex(power)?,
                });
            }
            reports.push(ThermalComponentResponse {
                net_heat_w: complex(drive.net_heat[component])?,
                energy_j: complex(energy)?,
                energy_residual_j: complex(residual)?,
                power_residual_w: complex(power)?,
            });
        }
        checkpoint().map_err(ThermalEnergyConstraintError::Checkpoint)?;
        Ok(reports)
    }
}

#[cfg(test)]
#[path = "conservation_tests.rs"]
mod tests;
