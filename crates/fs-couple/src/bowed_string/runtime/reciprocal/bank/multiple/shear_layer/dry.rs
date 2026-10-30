//! Moving material-coordinate contact heat on the existing bow/body owner.

use super::*;
use crate::render::schedule::force::coupled::contact::friction::multiple::layer::thermal::dry::DryContactHeat;

/// Coordinate of the actual opposite surface. A stationary body has a held
/// material footprint. A source fiber ribbon follows its actual material
/// crossing, including acceleration and reversal, using the SAME held map.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum DryReceiverHeatPosition {
    Fixed,
    Fiber { coordinate_origin_m: f64 },
}

/// A source bow's thermal strips. The first footprint follows the actual string
/// contact in material metres. Its original interval supplies its finite length;
/// `string_coordinate_origin_m` converts the string's source origin to the strip
/// coordinates. Width and nominal contact area stay unchanged. Receiver strips
/// are explicit, never inferred from a mechanical spring or modal mass.
#[derive(Clone, Debug)]
pub struct BowDryContactHeat<'a> {
    pub heat: DryContactHeat<'a>,
    pub string_coordinate_origin_m: f64,
    pub receiver: DryReceiverHeatPosition,
}
impl<'a> BowDryContactHeat<'a> {
    fn at<E>(
        &self,
        string_centres_m: &[f64],
        fiber_centre_m: Option<f64>,
    ) -> Result<DryContactHeat<'a>, LayerThermalError<E>> {
        let bad = LayerThermalError::Invalid;
        let string = *string_centres_m
            .get(self.heat.port)
            .ok_or(bad("dry thermal bow port is absent"))?;
        let mut heat = self.heat.clone();
        let mut place =
            |side: usize, origin: f64, position: f64| -> Result<(), LayerThermalError<E>> {
                let original = heat.footprint.intervals_m[side];
                let half = 0.5 * (original[1] - original[0]);
                let centre = origin + position;
                let mut interval = [centre - half, centre + half];
                if !centre.is_finite()
                    || !half.is_finite()
                    || half <= 0.0
                    || interval.iter().any(|x| !x.is_finite())
                    || interval[0] >= interval[1]
                {
                    return Err(bad("moving dry footprint is unrepresentable"));
                }
                // Cancellation near a local zero retains rounding error from the
                // ORIGINAL translation operands. Snap only a coincident cell edge,
                // translating the whole footprint so its physical area is retained.
                let roundoff = 8.0 * f64::EPSILON * origin.abs().max(position.abs()).max(half);
                'snap: for edge in 0..2 {
                    for boundary in heat.surfaces[side].iter().flat_map(|cell| cell.interval_m) {
                        if boundary.is_finite() && (interval[edge] - boundary).abs() <= roundoff {
                            interval = if edge == 0 {
                                [boundary, boundary + 2.0 * half]
                            } else {
                                [boundary - 2.0 * half, boundary]
                            };
                            break 'snap;
                        }
                    }
                }
                heat.footprint.intervals_m[side] = interval;
                Ok(())
            };
        place(0, self.string_coordinate_origin_m, string)?;
        match (self.receiver, fiber_centre_m) {
            (DryReceiverHeatPosition::Fixed, None) => {}
            (
                DryReceiverHeatPosition::Fiber {
                    coordinate_origin_m,
                },
                Some(centre),
            ) => place(1, coordinate_origin_m, centre)?,
            _ => {
                return Err(bad(
                    "fiber heat must follow its actual source material crossing; a stationary receiver cannot substitute for it",
                ));
            }
        }
        Ok(heat)
    }
}

impl MultiBowStringBank {
    /// Advance moving dry bow/string/fiber footprints, contact work, original
    /// finite heat cells and pressure together. The usual station projection,
    /// normal co-location, support inertia and endpoint source checks are shared
    /// with ordinary rendering. A failed thermal/source solve advances none of
    /// the mechanical state, station cursor, fiber support, heat or observer.
    ///
    /// This is dry heating with a massless isothermal contact plane, not rosin
    /// transport. Moving material shear layers remain refused. Source dry laws
    /// keep their original coefficients and are re-admitted at the solved T.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_dry_contact_thermal_motion<'t, S, E: core::fmt::Display>(
        &mut self,
        thermal: &mut DepositionSimulation<'t>,
        dry: &[BowDryContactHeat<'_>],
        series_heat: &[SeriesDriveHeatDestination],
        controls: LayerThermalSolve,
        energy_guess_j: &[f64],
        external: &[f64],
        drives: &[FrictionDrive],
        bow_stations: Option<&[f64]>,
        stop_stations: Option<&[f64]>,
        offsets_m: Option<&[f64]>,
        gate: &CancelGate,
        resolve: impl Fn(&DepositionSimulation<'t>) -> Result<StepCoefficients<S>, E>,
    ) -> Result<LayerThermalFrame<S>, RenderError> {
        let normal = stop_stations
            .map(|stations| {
                self.prepare_stop_station_geometry_with_bows(stations, bow_stations, gate)
            })
            .transpose()?;
        let failure = |error: LayerThermalError<E>| {
            RenderError::Voice(AcousticRealizeError::Nonlinear(error.to_string()))
        };
        let accepted = if self.bow_station_motion_limits().is_some() {
            self.with_bow_geometry(
                external,
                drives,
                bow_stations,
                normal.as_deref(),
                offsets_m,
                gate,
                |system, step| {
                    let routes = dry
                        .iter()
                        .map(|route| route.at(step.string_centres_m, step.fiber_centre_m))
                        .collect::<Result<Vec<_>, LayerThermalError<E>>>()?;
                    system.step_with_dry_contact_thermal_geometry(
                        thermal,
                        &routes,
                        series_heat,
                        controls,
                        energy_guess_j,
                        step.external,
                        step.drives,
                        step.geometry,
                        step.endpoints,
                        step.normal_geometry,
                        step.stop_offsets_m,
                        step.config,
                        gate,
                        resolve,
                    )
                },
            )
            .map_err(failure)?
        } else {
            if bow_stations.is_some() {
                return Err(RenderError::Control {
                    what: "bow migration needs admitted station motion",
                });
            }
            let centres = self
                .source
                .bows
                .iter()
                .map(|bow| bow.location.station_fraction * self.cards[bow.location.string].length_m)
                .collect::<Vec<_>>();
            let fiber_centre = self
                .source
                .fiber_ribbon()
                .map(|fiber| fiber.station_fraction() * fiber.model().geometry().length_m);
            let routes = dry
                .iter()
                .map(|route| route.at(&centres, fiber_centre))
                .collect::<Result<Vec<_>, LayerThermalError<E>>>()
                .map_err(failure)?;
            self.system
                .step_with_dry_contact_thermal_network(
                    thermal,
                    &routes,
                    series_heat,
                    controls,
                    energy_guess_j,
                    external,
                    drives,
                    normal.as_deref(),
                    offsets_m,
                    gate,
                    resolve,
                )
                .map_err(failure)?
        };
        if let Some(stations) = stop_stations {
            self.stop_stations
                .as_mut()
                .expect("admitted stop geometry")
                .accepted
                .copy_from_slice(stations);
        }
        Ok(accepted)
    }
}
