//! Circular compact radiation attached to the original resolved NSF port.
//!
//! The caller declares a circular reduced exterior and its flange condition.
//! The selected Cartesian face or aperture supplies its equal-area uniform
//! acoustic port; area agreement does not establish a circular resolved mesh.
//! Density, sound speed and the physical interval come from the original grid.
//! The existing radiation owner supplies all coefficients and model checks:
//! ka <= 0.5, dt*fmax <= 0.1 and sampled analogue/discrete response limits.
//! The returned model retains those limits and their measured discrepancies.
//! No transient bandwidth, exterior field, mass or enthalpy inventory is inferred.

use super::NsfImpedance;
use super::super::{NsfDgError, NsfDgLimits, poll};
use crate::bernoulli_aperture::radiation::{BaffledRadiationLoad, RadiationFlange};
use fs_exec::CancelGate;
use fs_lbm::lattice3::Lattice3;
use fs_lbm::nsf3::NsfGrid3;

/// An explicit circular radiation model, with no inferred flange or use band.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CircularRadiationSpec {
    /// Physical radius of the declared circular reduced exterior [m].
    /// Its area must match the selected original uniform acoustic port.
    pub radius_m: f64,
    /// Positive declared use band [Hz]; this does not filter a transient.
    pub maximum_frequency_hz: f64,
    /// Baffled piston or thin-wall unflanged circular pipe end.
    pub flange: RadiationFlange,
}

impl<L: Lattice3<27>> NsfImpedance<L> {
    /// Couple a complete original face to its declared circular radiation load.
    ///
    /// The gas and clock are resolved from the same NSF owner. The returned
    /// immutable model reports the original coefficient source and admitted
    /// band; near-field storage and radiated energy remain in the causal load.
    /// This is a uniform acoustic projection, not a resolved circular aperture.
    pub fn new_with_circular_radiation(
        zone: NsfGrid3<L>,
        radiation: CircularRadiationSpec,
        high: bool,
        limits: NsfDgLimits,
        gate: &CancelGate,
    ) -> Result<(Self, BaffledRadiationLoad), NsfDgError> {
        let model = prepare(&zone, radiation, high, false, limits, gate)?;
        let system = Self::new(zone, model.load(), high, limits, gate)?;
        Ok((system, model))
    }

    /// Couple only the selected original aperture to circular radiation.
    ///
    /// The complementary opening retains its original declared boundary drive
    /// and contributes no flow or work to this radiation owner.
    pub fn new_on_patch_with_circular_radiation(
        zone: NsfGrid3<L>,
        radiation: CircularRadiationSpec,
        high: bool,
        limits: NsfDgLimits,
        gate: &CancelGate,
    ) -> Result<(Self, BaffledRadiationLoad), NsfDgError> {
        let model = prepare(&zone, radiation, high, true, limits, gate)?;
        let system = Self::new_on_patch(zone, model.load(), high, limits, gate)?;
        Ok((system, model))
    }
}

fn prepare<L: Lattice3<27>>(
    zone: &NsfGrid3<L>,
    radiation: CircularRadiationSpec,
    high: bool,
    patch: bool,
    limits: NsfDgLimits,
    gate: &CancelGate,
) -> Result<BaffledRadiationLoad, NsfDgError> {
    poll(gate)?;
    limits.validate_for_cells(zone.cell_count())?;
    let reference = zone
        .characteristic_reference_gas(high)
        .map_err(NsfDgError::Flow)?;
    let traces = if patch {
        zone.open_patch_traces_on(high)
    } else {
        zone.open_boundary_traces(high)
    }
    .ok_or(NsfDgError::Input(
        "circular radiation requires the original selected characteristic quadrature",
    ))?;
    let area = traces
        .first()
        .ok_or(NsfDgError::Input("circular radiation port is empty"))?
        .area_m2
        * traces.len() as f64;
    let declared_area = core::f64::consts::PI * radiation.radius_m * radiation.radius_m;
    if !radiation.radius_m.is_finite()
        || radiation.radius_m <= 0.0
        || !area.is_finite()
        || area <= 0.0
        || !declared_area.is_finite()
        || declared_area <= 0.0
        || (area - declared_area).abs()
            > 32.0 * f64::EPSILON * area.max(declared_area)
    {
        return Err(NsfDgError::Input(
            "declared circular radiation area must match the original selected acoustic port",
        ));
    }
    let dt = zone.law().dt_s * limits.substeps as f64;
    if !dt.is_finite() || dt <= 0.0 {
        return Err(NsfDgError::Input(
            "circular radiation requires a representable original physical interval",
        ));
    }
    let model = BaffledRadiationLoad::with_flange(
        radiation.radius_m,
        reference.density,
        reference.sound_speed,
        dt,
        radiation.maximum_frequency_hz,
        radiation.flange,
        gate,
    )
    .map_err(NsfDgError::Acoustic)?;
    poll(gate)?;
    Ok(model)
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_lbm::d3q27::D3q27;
    use fs_lbm::nsf3::{
        NsfBoundaryDrive, NsfCharacteristicAxis, NsfCharacteristicBoundary,
        NsfCharacteristicPatch, NsfSpec, WallThermal,
    };
    use fs_material::gas::{GasSpec, GasState};

    fn zone(pressure: f64, high: bool, patch: bool) -> NsfGrid3<D3q27> {
        let gas = GasSpec::dry_air_ussa1976();
        let reference = GasState::try_new(&gas, 300.0, pressure).unwrap();
        let velocity_unit = fs_math::det::sqrt(
            3.0 * reference.specific_gas_constant * reference.temperature,
        );
        let mut grid = NsfGrid3::new(
            NsfSpec {
                dims: [8, 2, 2],
                periodic: [false, true, true],
                gas,
                temperature_k: reference.temperature,
                pressure_pa: reference.pressure,
                dx_m: velocity_unit / 96000.0,
                temperature_window_k: [250.0, 600.0],
                maximum_mach: 0.2,
                maximum_temperature_ratio: 1.0,
                bulk_to_shear_ratio: 0.73,
                bulk_ratio_source: "authored circular radiation coupling fixture",
                maximum_cells: 32,
                wall: WallThermal::Adiabatic,
                workers: 1,
                transport_points: 17,
            },
            |_| false,
            |_| (1.0, [0.0; 3], reference.temperature),
        )
        .unwrap()
        .with_characteristic_axis(NsfCharacteristicAxis {
            axis: 0,
            faces: [NsfCharacteristicBoundary {
                reference_pressure_pa: reference.pressure,
                reference_temperature_k: reference.temperature,
                reference_velocity_m_s: [0.0; 3],
                incoming_characteristic_m_s: 0.0,
            }; 2],
            maximum_relative_perturbation: 0.05,
        })
        .unwrap();
        if patch {
            grid = grid
                .with_characteristic_patch(
                    high,
                    NsfCharacteristicPatch { start: [0, 0], end: [1, 1] },
                    NsfBoundaryDrive::PressureOffset(0.0),
                )
                .unwrap();
        }
        grid
    }

    fn limits() -> NsfDgLimits {
        NsfDgLimits {
            substeps: 2,
            maximum_cells: 32,
            maximum_cell_steps: 32 * 2 * 64,
            maximum_evaluations: 64,
            maximum_steps: 32,
            incoming_bracket_m_s: [-0.05, 0.05],
            pressure_absolute_tolerance_pa: 1e-7,
            pressure_relative_tolerance: 1e-9,
            projection_absolute_tolerance_j: 0.0,
            projection_relative_tolerance: 1.0,
        }
    }

    fn description(
        grid: &NsfGrid3<D3q27>,
        patch: bool,
        flange: RadiationFlange,
    ) -> CircularRadiationSpec {
        let count = if patch { 1.0 } else { 4.0 };
        CircularRadiationSpec {
            radius_m: fs_math::det::sqrt(
                count * grid.law().dx_m * grid.law().dx_m / core::f64::consts::PI,
            ),
            maximum_frequency_hz: 1000.0,
            flange,
        }
    }

    /// G1/G3: original medium, selected area and physical clock determine the load.
    #[test]
    fn original_gas_and_aperture_bind_both_circular_exterior_models() {
        let gate = CancelGate::new_clock_free();
        for high in [false, true] {
            for patch in [false, true] {
                for flange in [RadiationFlange::Baffled, RadiationFlange::Unflanged] {
                    let mut resistances = [0.0; 2];
                    for (index, pressure) in [100000.0, 200000.0].into_iter().enumerate() {
                        let grid = zone(pressure, high, patch);
                        let dt = grid.law().dt_s * limits().substeps as f64;
                        let gas = grid.characteristic_reference_gas(high).unwrap();
                        let spec = description(&grid, patch, flange);
                        let (system, model) = if patch {
                            NsfImpedance::new_on_patch_with_circular_radiation(
                                grid, spec, high, limits(), &gate,
                            )
                        } else {
                            NsfImpedance::new_with_circular_radiation(
                                grid, spec, high, limits(), &gate,
                            )
                        }
                        .unwrap();
                        assert_eq!(model.flange(), flange);
                        assert_eq!(model.density_kg_m3(), gas.density);
                        assert_eq!(model.sound_speed_m_s(), gas.sound_speed);
                        assert_eq!(model.time_step_s(), dt);
                        assert_eq!(system.time_step_s(), dt);
                        assert_eq!(system.load().spec(), &model.load());
                        let area = core::f64::consts::PI * spec.radius_m * spec.radius_m;
                        assert!((system.area_m2() / area - 1.0).abs() < 2e-14);
                        resistances[index] = model.load().terms()[0].resistance_pa_s_m3;
                    }
                    assert!((resistances[1] / resistances[0] - 2.0).abs() < 2e-13);
                }
            }
        }
    }

    /// G4: undeclared geometry, unsupported bands and cancellation refuse.
    #[test]
    fn circular_radiation_refuses_area_mismatch_missing_patch_and_unresolved_band() {
        let gate = CancelGate::new_clock_free();
        let grid = zone(100000.0, true, false);
        let mut wrong_area = description(&grid, false, RadiationFlange::Baffled);
        wrong_area.radius_m *= 1.01;
        assert!(NsfImpedance::new_with_circular_radiation(
            grid, wrong_area, true, limits(), &gate,
        ).is_err());

        let grid = zone(100000.0, true, false);
        let spec = description(&grid, false, RadiationFlange::Baffled);
        assert!(NsfImpedance::new_on_patch_with_circular_radiation(
            grid, spec, true, limits(), &gate,
        ).is_err());

        let grid = zone(100000.0, true, true);
        let mut too_wide = description(&grid, true, RadiationFlange::Unflanged);
        too_wide.maximum_frequency_hz = 100000.0;
        assert!(NsfImpedance::new_on_patch_with_circular_radiation(
            grid, too_wide, true, limits(), &gate,
        ).is_err());

        let grid = zone(100000.0, false, false);
        let spec = description(&grid, false, RadiationFlange::Baffled);
        gate.request();
        assert!(matches!(
            NsfImpedance::new_with_circular_radiation(grid, spec, false, limits(), &gate),
            Err(NsfDgError::Cancelled)
        ));
    }
}
