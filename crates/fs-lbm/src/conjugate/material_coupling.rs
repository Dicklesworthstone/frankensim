//! Material-law feedback through the existing electrical and thermal owners.
//!
//! Every electrical evaluation and thermal coefficient update sees an explicit
//! temperature field. No material-name defaults, constant-property inference,
//! source extrapolation, or new current/heat discretization is introduced here.

use fs_exec::CancelGate;

use super::{
    ChtError, ElectricSetup, ElectricSolution, EnergySolution, ThermalSetup, VoxelDomain,
    finite, poll,
};

/// Explicit numerical budget for a frozen-geometry electrothermal fixed point.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MaterialCouplingConfig {
    /// Maximum number of thermal solves; must be positive.
    pub max_passes: usize,
    /// Absolute temperature-update tolerance, K.
    pub temperature_tolerance_k: f64,
    /// Absolute per-cell Joule-source mismatch tolerance, W.
    pub power_tolerance_w: f64,
    /// Relative tolerance, added against temperature span and peak cell power.
    /// Must lie in `[0, 1)`; zero requests only the absolute tolerances.
    pub relative_tolerance: f64,
}

/// Measured defects at the returned fixed point, not a physical certificate.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MaterialCouplingReport {
    /// Thermal solves performed, including the accepted pass.
    pub passes: usize,
    /// Maximum temperature change across cells, K.
    pub temperature_change_k: f64,
    /// Maximum difference between the Joule heat used by the thermal solve and
    /// the electrical solve re-evaluated at its returned temperatures, W/cell.
    pub joule_source_change_w: f64,
}

/// Coupled fields, with the electrical field evaluated at the returned thermal
/// temperatures. The heat-source consistency defect remains explicit.
#[derive(Debug, Clone, PartialEq)]
pub struct MaterialCoupledSolution {
    /// Temperature field and the thermal owner's accounting report.
    pub thermal: EnergySolution,
    /// Potential/current/Joule heat at `thermal.temperature`.
    pub electric: ElectricSolution,
    /// Actual coupling defects and work performed.
    pub report: MaterialCouplingReport,
}

impl MaterialCouplingConfig {
    fn validate(self) -> Result<(), ChtError> {
        if self.max_passes == 0 {
            return Err(invalid("material_coupling.max_passes", "must be positive"));
        }
        for (field, value) in [
            (
                "material_coupling.temperature_tolerance_k",
                self.temperature_tolerance_k,
            ),
            ("material_coupling.power_tolerance_w", self.power_tolerance_w),
            ("material_coupling.relative_tolerance", self.relative_tolerance),
        ] {
            if !value.is_finite() || value < 0.0 {
                return Err(invalid(field, "must be finite and non-negative"));
            }
        }
        if self.relative_tolerance >= 1.0 {
            return Err(invalid(
                "material_coupling.relative_tolerance",
                "must be less than one",
            ));
        }
        Ok(())
    }
}

/// Couple a material-resolved electrical law to a state-aware thermal owner.
///
/// `resistivity(cell, material_index, temperature_k)` has the same semantics as
/// [`ElectricSetup::solve_with_resistivity`]. `energy(heated, temperature_k)`
/// receives the current coefficient-evaluation state as well as the combined
/// authored and Joule sources. It may resolve thermal/flow properties from the
/// same cards and state before calling the existing energy or conjugate owner.
/// Neither closure is permitted to silently invent unsupported material states;
/// source-owner errors should be propagated as `ChtError` with source context.
///
/// The driver re-evaluates the electrical law at *every proposed new thermal
/// state*, including the final one. A small temperature update alone cannot
/// accept a large change in electrical heating. The accepted electrical field
/// is the re-evaluated field, not one cached at the previous temperature.
///
/// One electrical solve is cached between passes at exactly the same state;
/// no stale-temperature material cache or duplicate discretization is used.
/// Source receipts captured by either callback remain caller-owned. The caller
/// must keep external coordinates (pressure, composition, etc.) consistent.
///
/// Fixed geometry, DC ohmic conduction, steady Picard coupling. No melting,
/// chemical evolution, topology changes, stability theorem or experimental
/// qualification is implied. Source tables are not extended by convergence.
///
/// # Errors
/// Invalid budgets or field shapes, non-finite/non-positive temperatures,
/// non-finite sources, either owner's refusal, cancellation, or exhaustion of
/// the coupling budget. Inputs are not mutated and no partial field is returned.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)] // explicit owners, state and budgets
pub fn solve_material_electrothermal(
    domain: &VoxelDomain,
    electric: &ElectricSetup,
    thermal: &ThermalSetup,
    initial_temperature_k: &[f64],
    config: MaterialCouplingConfig,
    mut resistivity: impl FnMut(usize, u16, f64) -> Result<f64, ChtError>,
    mut energy: impl FnMut(&ThermalSetup, &[f64]) -> Result<EnergySolution, ChtError>,
    gate: &CancelGate,
) -> Result<MaterialCoupledSolution, ChtError> {
    poll(gate)?;
    config.validate()?;
    let cells = domain.cell_count();
    validate_temperature(initial_temperature_k, cells)?;
    if !thermal.power_w.is_empty() && thermal.power_w.len() != cells {
        return Err(invalid(
            "material_coupling.power_w",
            "needs one source per cell or an empty field",
        ));
    }
    for &value in &thermal.power_w {
        finite("material_coupling.power_w", value)?;
    }
    let mut temperature = initial_temperature_k.to_vec();
    let mut field = electric.solve_with_resistivity(domain, &temperature, &mut resistivity, gate)?;
    let mut heated = thermal.clone();
    let mut residual = f64::INFINITY;
    for pass in 1..=config.max_passes {
        poll(gate)?;
        validate_joule(&field, cells)?;
        heated.power_w = field.joule_w.clone();
        if !thermal.power_w.is_empty() {
            for (sum, &authored) in heated.power_w.iter_mut().zip(&thermal.power_w) {
                *sum += authored;
                finite("material_coupling.combined_power_w", *sum)?;
            }
        }
        let solution = energy(&heated, &temperature)?;
        poll(gate)?;
        validate_temperature(&solution.temperature, cells)?;
        let next = electric.solve_with_resistivity(
            domain,
            &solution.temperature,
            &mut resistivity,
            gate,
        )?;
        validate_joule(&next, cells)?;
        let temperature_change = maximum_difference(&temperature, &solution.temperature);
        let power_change = maximum_difference(&field.joule_w, &next.joule_w);
        let (lo, hi) = solution.temperature.iter().fold(
            (f64::INFINITY, f64::NEG_INFINITY),
            |(lo, hi), &t| (lo.min(t), hi.max(t)),
        );
        let peak_power = field
            .joule_w
            .iter()
            .chain(&next.joule_w)
            .copied()
            .fold(0.0, f64::max);
        let temperature_limit =
            config.temperature_tolerance_k + config.relative_tolerance * (hi - lo);
        let power_limit = config.power_tolerance_w + config.relative_tolerance * peak_power;
        finite("material_coupling.temperature_limit_k", temperature_limit)?;
        finite("material_coupling.power_limit_w", power_limit)?;
        residual = scaled_defect(temperature_change, temperature_limit)
            .max(scaled_defect(power_change, power_limit));
        if temperature_change <= temperature_limit && power_change <= power_limit {
            poll(gate)?;
            return Ok(MaterialCoupledSolution {
                thermal: solution,
                electric: next,
                report: MaterialCouplingReport {
                    passes: pass,
                    temperature_change_k: temperature_change,
                    joule_source_change_w: power_change,
                },
            });
        }
        temperature = solution.temperature;
        field = next;
    }
    Err(ChtError::SolverNotConverged {
        system: "material electrothermal coupling",
        iterations: config.max_passes,
        relative_residual: residual,
        tolerance: 1.0,
    })
}

fn invalid(field: &'static str, reason: &str) -> ChtError {
    ChtError::InvalidInput {
        field,
        reason: reason.into(),
    }
}

fn validate_temperature(values: &[f64], cells: usize) -> Result<(), ChtError> {
    if values.len() != cells {
        return Err(invalid(
            "material_coupling.temperature",
            "needs one temperature per cell",
        ));
    }
    if values.iter().any(|t| !t.is_finite() || *t <= 0.0) {
        return Err(invalid(
            "material_coupling.temperature",
            "all temperatures must be finite and positive kelvin",
        ));
    }
    Ok(())
}

fn validate_joule(field: &ElectricSolution, cells: usize) -> Result<(), ChtError> {
    if field.joule_w.len() != cells
        || field.joule_w.iter().any(|p| !p.is_finite() || *p < 0.0)
    {
        return Err(invalid(
            "material_coupling.joule_w",
            "needs one finite, non-negative Joule source per cell",
        ));
    }
    Ok(())
}

fn maximum_difference(left: &[f64], right: &[f64]) -> f64 {
    debug_assert_eq!(left.len(), right.len());
    left.iter()
        .zip(right)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0, f64::max)
}

fn scaled_defect(change: f64, limit: f64) -> f64 {
    if change == 0.0 { 0.0 } else { change / limit }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::conjugate::{
        Conductor, Electrode, ElectrodeDrive, EnergyConfig, FlowField, FluidProperties,
        SolidMaterial, ThermalFace, Voxel, solve_energy,
    };

    fn fixture() -> (VoxelDomain, ElectricSetup, ThermalSetup) {
        let domain = VoxelDomain::from_fn(3, 1, 1, 1.0, |_| Voxel::Solid(0)).unwrap();
        let electric = ElectricSetup {
            conductors: vec![Some(Conductor::new(2.0))],
            electrodes: vec![
                Electrode {
                    lo: [0, 0, 0],
                    hi: [1, 1, 1],
                    drive: ElectrodeDrive::Voltage(0.0),
                },
                Electrode {
                    lo: [2, 0, 0],
                    hi: [3, 1, 1],
                    drive: ElectrodeDrive::Voltage(3.0),
                },
            ],
        };
        let mut faces = [ThermalFace::Adiabatic; 6];
        faces[0] = ThermalFace::Temperature(300.0);
        faces[1] = ThermalFace::Temperature(300.0);
        (domain, electric, ThermalSetup::new(faces))
    }

    fn config() -> MaterialCouplingConfig {
        MaterialCouplingConfig {
            max_passes: 100,
            temperature_tolerance_k: 1e-9,
            power_tolerance_w: 1e-9,
            relative_tolerance: 1e-10,
        }
    }

    fn heat(
        domain: &VoxelDomain,
        thermal: &ThermalSetup,
        gate: &CancelGate,
    ) -> Result<EnergySolution, ChtError> {
        solve_energy(
            domain,
            &FluidProperties {
                density_kg_m3: 1.0,
                specific_heat_j_kg_k: 1.0,
                conductivity_w_m_k: 1.0,
                kinematic_viscosity_m2_s: 1.0,
            },
            &[SolidMaterial::new("analytic unit conductor", 1.0)],
            &FlowField::quiescent(domain),
            thermal,
            &EnergyConfig::default(),
            gate,
        )
    }

    #[test]
    fn g1_g3_nonlinear_material_swap_reaches_analytic_electrothermal_fixed_point() {
        let (domain, electric, thermal) = fixture();
        let gate = CancelGate::new();
        let original = thermal.clone();
        for scale in [1.0, 2.0] {
            let law = |t: f64| {
                let rise = t - 300.0;
                scale * 2.0 * (1.0 + 0.1 * rise + 0.02 * rise * rise)
            };
            let result = solve_material_electrothermal(
                &domain,
                &electric,
                &thermal,
                &[300.0; 3],
                config(),
                |_, _, t| Ok(law(t)),
                |heated, state| {
                    assert_eq!(state.len(), 3);
                    heat(&domain, heated, &gate)
                },
                &gate,
            )
            .unwrap();
            // Centre-to-ambient thermal resistance is 0.75 K/W (the two
            // 1.5-m paths are in parallel). Therefore the rise d satisfies
            // d*(1+.1*d+.02*d^2)=3.375/scale.
            let (mut lo, mut hi) = (0.0, 4.0);
            for _ in 0..80 {
                let mid = 0.5 * (lo + hi);
                if mid * (1.0 + 0.1 * mid + 0.02 * mid * mid) < 3.375 / scale {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            assert!((result.thermal.temperature[1] - 300.0 - 0.5 * (lo + hi)).abs() < 1e-8);
            let power = 9.0 / law(result.thermal.temperature[1]);
            assert!((result.electric.total_joule_w - power).abs() < 1e-10);
            assert!((result.electric.delivered_w - power).abs() < 1e-10);
            assert!(result.report.passes > 1);
            assert!(result.report.temperature_change_k < 2e-9);
            assert!(result.report.joule_source_change_w < 2e-9);
            assert!(result.thermal.report.balance.relative_residual < 1e-9);
        }
        assert_eq!(thermal, original);
    }

    #[test]
    fn g0_final_candidate_must_still_be_inside_material_domain() {
        let (domain, electric, thermal) = fixture();
        let gate = CancelGate::new();
        let refusal = invalid("material.temperature_domain", "source supports at most 301 K");
        let result = solve_material_electrothermal(
            &domain,
            &electric,
            &thermal,
            &[300.0; 3],
            config(),
            |_, _, t| {
                if t > 301.0 { Err(refusal.clone()) } else { Ok(2.0) }
            },
            |heated, _| heat(&domain, heated, &gate),
            &gate,
        );
        assert_eq!(result, Err(refusal));
    }

    #[test]
    fn g0_bad_thermal_state_cannot_be_accepted_as_converged() {
        let (domain, electric, thermal) = fixture();
        let gate = CancelGate::new();
        for bad_shape in [false, true] {
            let result = solve_material_electrothermal(
                &domain,
                &electric,
                &thermal,
                &[300.0; 3],
                config(),
                |_, _, _| Ok(2.0),
                |heated, _| {
                    let mut solution = heat(&domain, heated, &gate)?;
                    if bad_shape {
                        solution.temperature.pop();
                    } else {
                        solution.temperature[1] = f64::NAN;
                    }
                    Ok(solution)
                },
                &gate,
            );
            assert!(matches!(
                result,
                Err(ChtError::InvalidInput {
                    field: "material_coupling.temperature",
                    ..
                })
            ));
        }
    }

    #[test]
    fn g0_small_temperature_update_cannot_hide_a_large_joule_source_change() {
        let (domain, electric, thermal) = fixture();
        let gate = CancelGate::new();
        let result = solve_material_electrothermal(
            &domain,
            &electric,
            &thermal,
            &[300.0; 3],
            MaterialCouplingConfig {
                max_passes: 1,
                temperature_tolerance_k: 1e-6,
                power_tolerance_w: 1e-12,
                relative_tolerance: 0.0,
            },
            |_, _, t| Ok(if t > 300.0 { 4.0 } else { 2.0 }),
            |heated, _| {
                // Deliberately synthetic proposal isolates the coupling gate;
                // the other test exercises the real thermal fixed point.
                let mut solution = heat(&domain, heated, &gate)?;
                solution.temperature.fill(300.0 + 1e-8);
                Ok(solution)
            },
            &gate,
        );
        assert!(matches!(result, Err(ChtError::SolverNotConverged { .. })));
    }

    #[test]
    fn g0_source_shape_and_iteration_budget_refuse_before_evaluation() {
        let (domain, electric, mut thermal) = fixture();
        let gate = CancelGate::new();
        thermal.power_w = vec![1.0];
        let result = solve_material_electrothermal(
            &domain,
            &electric,
            &thermal,
            &[300.0; 3],
            config(),
            |_, _, _| panic!("malformed source must refuse before material evaluation"),
            |_, _| panic!("malformed source must refuse before a thermal solve"),
            &gate,
        );
        assert!(matches!(
            result,
            Err(ChtError::InvalidInput {
                field: "material_coupling.power_w",
                ..
            })
        ));
        assert!(MaterialCouplingConfig { max_passes: 0, ..config() }.validate().is_err());
        assert_eq!(scaled_defect(0.0, 0.0), 0.0);
        assert_eq!(scaled_defect(1.0, 0.0), f64::INFINITY);
    }
}
