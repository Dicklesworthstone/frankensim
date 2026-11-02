//! Centre-regular generated heat with the original material-integral owner.
//! Symmetry fixes Q(0)=0. Enclosed generation determines every heat rate;
//! the outer film fixes its face temperature, then exact signed Kirchhoff
//! inversion marches inward. No artificial inner radius or reservoir is used.

use super::super::geometry::{LayerShape, RadialKind, SolidCore};
use super::super::nonlinear::ConductivityIntegral;
use super::super::{
    LayeredWall, NonlinearGenerationOptions, checkpoint, finite, invalid, positive,
};
use super::{admit_shell_shape, reference_resistance, scaled};
use crate::ConductionError;
use fs_exec::Cx;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SolidRadialGeometry {
    /// Full circular cylinder with insulated endcaps.
    Cylindrical { length_m: f64 },
    /// Full solid sphere.
    Spherical,
}

/// A solid core followed by zero or more concentric material shells.
/// The first original layer thickness is the core radius. Its inner surface
/// resistance must be explicit zero because the centre is a symmetry point,
/// not a heat reservoir. All remaining contacts and the outer film retain
/// their original physical areal resistance and source context.
#[derive(Debug, Clone, PartialEq)]
pub struct SolidRadialLayeredBody {
    wall: LayeredWall,
    geometry: SolidRadialGeometry,
    core: SolidCore,
    shells: Vec<LayerShape>,
    radii_m: Vec<f64>,
    reference_area_m2: f64,
    contacts_reference_m2_k_w: Vec<f64>,
    outer_resistance_reference_m2_k_w: f64,
}

impl SolidRadialLayeredBody {
    pub fn cylindrical(wall: LayeredWall, length_m: f64) -> Result<Self, ConductionError> {
        positive("solid cylinder length", length_m)?;
        Self::new(wall, SolidRadialGeometry::Cylindrical { length_m })
    }

    pub fn spherical(wall: LayeredWall) -> Result<Self, ConductionError> {
        Self::new(wall, SolidRadialGeometry::Spherical)
    }

    fn new(wall: LayeredWall, geometry: SolidRadialGeometry) -> Result<Self, ConductionError> {
        if wall.surface_resistances_m2_k_per_w[0] != 0.0 {
            return Err(invalid(
                "a solid radial centre cannot carry an inner surface film",
            ));
        }
        let radius = wall.layers[0].thickness_m;
        let (kind, area) = match geometry {
            SolidRadialGeometry::Cylindrical { length_m } => (
                RadialKind::Cylindrical,
                (2.0 * core::f64::consts::PI * radius) * length_m,
            ),
            SolidRadialGeometry::Spherical => (
                RadialKind::Spherical,
                (4.0 * core::f64::consts::PI * radius) * radius,
            ),
        };
        let core = SolidCore::new(kind, radius)
            .ok_or_else(|| invalid("solid core geometry is not representable"))?;
        positive("solid core outer area", area)?;
        positive("solid core physical volume", core.volume() * area)?;
        let mut radii_m = vec![0.0, radius];
        let mut shells = Vec::with_capacity(wall.layers.len() - 1);
        let mut contacts = Vec::with_capacity(wall.layers.len() - 1);
        if let Some(&resistance) = wall.contact_resistances_m2_k_per_w.first() {
            contacts.push(resistance);
        }
        let mut outer = radius;
        let mut outer_ratio = 1.0;
        for (index, layer) in wall.layers.iter().enumerate().skip(1) {
            let shape = admit_shell_shape(kind, outer, radius, layer.thickness_m, area)?;
            outer = finite("outer solid layer radius", outer + layer.thickness_m)?;
            outer_ratio = shape.area_ratio(layer.thickness_m);
            if let Some(&resistance) = wall.contact_resistances_m2_k_per_w.get(index) {
                contacts.push(reference_resistance(resistance, outer_ratio)?);
            }
            shells.push(shape);
            radii_m.push(outer);
        }
        let outer_resistance_reference_m2_k_w =
            reference_resistance(wall.surface_resistances_m2_k_per_w[1], outer_ratio)?;
        Ok(Self {
            wall,
            geometry,
            core,
            shells,
            radii_m,
            reference_area_m2: area,
            contacts_reference_m2_k_w: contacts,
            outer_resistance_reference_m2_k_w,
        })
    }

    pub fn wall(&self) -> &LayeredWall {
        &self.wall
    }
    pub const fn geometry(&self) -> SolidRadialGeometry {
        self.geometry
    }
    pub fn radii_m(&self) -> &[f64] {
        &self.radii_m
    }
    /// Outer area of the solid core, used only for numerical normalization.
    pub const fn reference_area_m2(&self) -> f64 {
        self.reference_area_m2
    }

    fn volume_length(&self, layer: usize) -> f64 {
        if layer == 0 {
            self.core.volume()
        } else {
            self.shells[layer - 1].volume()
        }
    }

    fn drop(&self, layer: usize, inlet: f64, generated: f64, offset: f64) -> f64 {
        if layer == 0 {
            self.core.potential_drop(generated, offset)
        } else {
            self.shells[layer - 1].potential_drop(inlet, generated, offset)
        }
    }

    fn flux(&self, layer: usize, inlet: f64, generated: f64, offset: f64) -> f64 {
        if layer == 0 {
            self.core.local_flux(generated, offset)
        } else {
            let shape = self.shells[layer - 1];
            (inlet + generated * (shape.volume_length(offset) / shape.volume()))
                / shape.area_ratio(offset)
        }
    }

    /// Steady solid-cylinder/sphere conduction with uniform signed generation
    /// per original material layer and one physical outer reservoir.
    ///
    /// Positive conductivity makes each retained Kirchhoff inverse unique.
    /// This direct inward pass needs no shooting iterations; the supplied
    /// options retain the original source interval and residual tolerances.
    /// Their absolute flux tolerance is referenced to the CORE OUTER area.
    /// Work/cancellation remain bounded by the original layer/knot caps.
    /// The model assumes fixed radial symmetry and generation independent of
    /// temperature. It supplies no transient storage or phase-front evolution.
    pub fn solve_temperature_dependent_with_generation(
        &self,
        cx: &Cx<'_>,
        outer_reservoir_temperature_k: f64,
        volumetric_generation_w_m3: &[f64],
        options: NonlinearGenerationOptions,
    ) -> Result<GeneratedSolidRadialSolution, ConductionError> {
        checkpoint(cx, 0)?;
        options.solver.validate()?;
        positive(
            "outer reservoir absolute temperature",
            outer_reservoir_temperature_k,
        )?;
        let [low, high] = options.temperature_interval_k;
        positive("minimum admitted material temperature", low)?;
        positive("maximum admitted material temperature", high)?;
        if low >= high {
            return Err(invalid(
                "material temperature interval must have positive width",
            ));
        }
        if volumetric_generation_w_m3.len() != self.wall.layers.len() {
            return Err(invalid(
                "supply exactly one volumetric source per solid radial layer",
            ));
        }
        let mut laws = Vec::with_capacity(self.wall.layers.len());
        let mut generated = Vec::with_capacity(self.wall.layers.len());
        let mut physical_generated = Vec::with_capacity(self.wall.layers.len());
        let (mut generated_heat, mut heat_correction) = (0.0, 0.0);
        let mut faces = Vec::with_capacity(self.wall.layers.len() + 1);
        faces.push(0.0);
        for (index, (layer, &source)) in self
            .wall
            .layers
            .iter()
            .zip(volumetric_generation_w_m3)
            .enumerate()
        {
            checkpoint(cx, index)?;
            finite("solid radial volumetric generation", source)?;
            let amount = finite(
                "solid radial generation per reference area",
                source * self.volume_length(index),
            )?;
            if source != 0.0 && amount == 0.0 {
                return Err(invalid("nonzero solid radial generation underflowed"));
            }
            faces.push(finite("enclosed radial generation", faces[index] + amount)?);
            generated.push(amount);
            // Independent dimensional source inventory; do not report a
            // tautological heat-minus-itself energy residual.
            let power = scaled(
                source,
                self.volume_length(index) * self.reference_area_m2,
                "volume-integrated solid generation",
            )?;
            let adjusted = finite("compensated solid generation", power - heat_correction)?;
            let next = finite("total physical solid generation", generated_heat + adjusted)?;
            heat_correction = finite(
                "solid generation summation correction",
                (next - generated_heat) - adjusted,
            )?;
            generated_heat = next;
            physical_generated.push(power);
            laws.push(ConductivityIntegral::new(
                cx,
                &layer.conductivity,
                low,
                high,
            )?);
        }
        let total = faces[self.wall.layers.len()];
        let film = self.outer_resistance_reference_m2_k_w;
        let outer_face = finite(
            "solid outer face temperature",
            total.mul_add(film, outer_reservoir_temperature_k),
        )?;
        let mut outer = outer_face;
        let mut rows = Vec::with_capacity(self.wall.layers.len());
        let mut maximum_residual = 0.0_f64;
        for index in (0..self.wall.layers.len()).rev() {
            checkpoint(cx, index)?;
            admitted(outer, [low, high])?;
            let layer = &self.wall.layers[index];
            let law = &laws[index];
            let inlet = faces[index];
            let outlet = faces[index + 1];
            let expected_drop = finite(
                "solid layer Kirchhoff drop",
                self.drop(index, inlet, generated[index], layer.thickness_m),
            )?;
            let expected_flux = if index == 0 {
                generated[index]
            } else {
                self.shells[index - 1].weighted_flux(inlet, generated[index])
            };
            if expected_flux != 0.0 && expected_drop == 0.0 {
                return Err(invalid("nonzero solid layer Kirchhoff drop underflowed"));
            }
            let inner = law.shift(cx, outer, -expected_drop)?.ok_or_else(|| {
                invalid("solid radial profile leaves the admitted material interval")
            })?;
            let (mut minimum, mut min_x, mut maximum, mut max_x) = if inner <= outer {
                (inner, 0.0, outer, layer.thickness_m)
            } else {
                (outer, layer.thickness_m, inner, 0.0)
            };
            if index > 0 && ((inlet < 0.0 && outlet > 0.0) || (inlet > 0.0 && outlet < 0.0)) {
                let fraction = finite(
                    "solid shell stationary volume fraction",
                    -inlet / generated[index],
                )?;
                let x = finite(
                    "solid shell stationary position",
                    self.shells[index - 1].stationary_position(fraction),
                )?;
                if fraction <= 0.0 || fraction >= 1.0 || x <= 0.0 || x >= layer.thickness_m {
                    return Err(invalid(
                        "solid shell interior extremum is not representable",
                    ));
                }
                let drop = finite(
                    "solid shell interior Kirchhoff drop",
                    self.drop(index, inlet, generated[index], x),
                )?;
                let temperature = law.shift(cx, inner, drop)?.ok_or_else(|| {
                    invalid("solid radial interior extremum leaves material support")
                })?;
                if temperature < minimum {
                    minimum = temperature;
                    min_x = x;
                }
                if temperature > maximum {
                    maximum = temperature;
                    max_x = x;
                }
            }
            let integral = law.integral(cx, inner.min(outer), inner.max(outer))?;
            let signed = if inner >= outer { integral } else { -integral };
            let recovered = if index == 0 {
                2.0 * (signed / self.core.radius())
            } else {
                signed / self.shells[index - 1].resistance()
            };
            let residual = finite(
                "solid radial constitutive flux residual",
                recovered - expected_flux,
            )?;
            maximum_residual = maximum_residual.max(residual.abs());
            let ratios = if index == 0 {
                [
                    self.core.area_ratio(0.0),
                    self.core.area_ratio(self.core.radius()),
                ]
            } else {
                let shape = self.shells[index - 1];
                [shape.area_ratio(0.0), shape.area_ratio(layer.thickness_m)]
            };
            rows.push(SolidRadialLayerTemperatures {
                face_radii_m: [self.radii_m[index], self.radii_m[index + 1]],
                face_areas_m2: ratios.map(|ratio| ratio * self.reference_area_m2),
                face_temperatures_k: [inner, outer],
                radial_heat_rates_w: [
                    scaled(inlet, self.reference_area_m2, "inner solid face heat rate")?,
                    scaled(outlet, self.reference_area_m2, "outer solid face heat rate")?,
                ],
                radial_heat_fluxes_w_per_m2: [
                    if index == 0 {
                        0.0
                    } else {
                        finite("inner solid local flux", inlet / ratios[0])?
                    },
                    finite("outer solid local flux", outlet / ratios[1])?,
                ],
                generated_heat_rate_w: physical_generated[index],
                minimum_temperature_k: minimum,
                minimum_radius_m: finite("solid minimum radius", self.radii_m[index] + min_x)?,
                maximum_temperature_k: maximum,
                maximum_radius_m: finite("solid maximum radius", self.radii_m[index] + max_x)?,
                constitutive_heat_rate_residual_w: finite(
                    "solid layer heat-rate defect",
                    residual * self.reference_area_m2,
                )?,
            });
            outer = inner;
            if index > 0 {
                let contact = self.contacts_reference_m2_k_w[index - 1];
                outer = finite(
                    "inner-side solid contact temperature",
                    inlet.mul_add(contact, inner),
                )?;
                if contact > 0.0 {
                    maximum_residual = maximum_residual.max(
                        finite(
                            "solid contact flux residual",
                            (outer - inner) / contact - inlet,
                        )?
                        .abs(),
                    );
                }
            }
        }
        rows.reverse();
        if film > 0.0 {
            maximum_residual = maximum_residual.max(
                finite(
                    "solid outer film flux residual",
                    (outer_face - outer_reservoir_temperature_k) / film - total,
                )?
                .abs(),
            );
        }
        let boundary_residual = finite(
            "solid outer reservoir residual",
            (-total).mul_add(film, outer_face) - outer_reservoir_temperature_k,
        )?;
        let scale = faces.iter().copied().map(f64::abs).fold(0.0_f64, f64::max);
        if maximum_residual > options.solver.flux_tolerance(scale)?
            || boundary_residual.abs() > options.solver.temperature_tolerance_k
        {
            return Err(invalid(
                "solid radial residual exceeds the declared numerical tolerance",
            ));
        }
        let heat = scaled(total, self.reference_area_m2, "solid outer reservoir heat")?;
        let energy_residual = finite("solid radial energy residual", heat - generated_heat)?;
        let energy_tolerance = finite(
            "solid radial energy tolerance",
            options.solver.flux_tolerance(scale)? * self.reference_area_m2,
        )?;
        if energy_residual.abs() > energy_tolerance {
            return Err(invalid(
                "solid radial source and reservoir heat fail energy closure",
            ));
        }
        checkpoint(cx, self.wall.layers.len())?;
        Ok(GeneratedSolidRadialSolution {
            body: self.clone(),
            outer_reservoir_temperature_k,
            volumetric_generation_w_m3: volumetric_generation_w_m3.to_vec(),
            reservoir_heat_rate_w: heat,
            generated_heat_rate_w: generated_heat,
            energy_residual_w: energy_residual,
            layer_temperatures: rows,
            maximum_heat_rate_residual_w: finite(
                "solid radial maximum heat-rate defect",
                maximum_residual * self.reference_area_m2,
            )?,
            outer_boundary_residual_k: boundary_residual,
            options,
            laws,
            generated,
            faces,
        })
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct SolidRadialLayerTemperatures {
    pub face_radii_m: [f64; 2],
    /// The centre has zero physical area and the regular zero-flux limit.
    pub face_areas_m2: [f64; 2],
    pub face_temperatures_k: [f64; 2],
    /// Signed toward increasing radius.
    pub radial_heat_rates_w: [f64; 2],
    pub radial_heat_fluxes_w_per_m2: [f64; 2],
    pub generated_heat_rate_w: f64,
    pub minimum_temperature_k: f64,
    pub minimum_radius_m: f64,
    pub maximum_temperature_k: f64,
    pub maximum_radius_m: f64,
    /// Recovered minus solved constitutive heat rate: outer heat for the core,
    /// reciprocal-area-weighted heat for a hollow shell. A numerical diagnostic.
    pub constitutive_heat_rate_residual_w: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct GeneratedSolidRadialSolution {
    pub body: SolidRadialLayeredBody,
    pub outer_reservoir_temperature_k: f64,
    pub volumetric_generation_w_m3: Vec<f64>,
    /// Signed heat delivered to the only physical reservoir.
    pub reservoir_heat_rate_w: f64,
    pub generated_heat_rate_w: f64,
    pub energy_residual_w: f64,
    pub layer_temperatures: Vec<SolidRadialLayerTemperatures>,
    pub maximum_heat_rate_residual_w: f64,
    pub outer_boundary_residual_k: f64,
    pub options: NonlinearGenerationOptions,
    laws: Vec<ConductivityIntegral>,
    generated: Vec<f64>,
    faces: Vec<f64>,
}

impl GeneratedSolidRadialSolution {
    pub fn temperature_at(
        &self,
        cx: &Cx<'_>,
        layer: usize,
        offset_m: f64,
    ) -> Result<f64, ConductionError> {
        self.check_query(cx, layer, offset_m)?;
        let row = &self.layer_temperatures[layer];
        if offset_m == 0.0 {
            return Ok(row.face_temperatures_k[0]);
        }
        if offset_m == self.body.wall.layers[layer].thickness_m {
            return Ok(row.face_temperatures_k[1]);
        }
        let drop = finite(
            "solid profile Kirchhoff drop",
            self.body
                .drop(layer, self.faces[layer], self.generated[layer], offset_m),
        )?;
        self.laws[layer]
            .shift(cx, row.face_temperatures_k[0], drop)?
            .ok_or_else(|| invalid("solid profile query left its solved material interval"))
    }

    pub fn heat_flux_at(
        &self,
        cx: &Cx<'_>,
        layer: usize,
        offset_m: f64,
    ) -> Result<f64, ConductionError> {
        self.check_query(cx, layer, offset_m)?;
        finite(
            "solid local radial heat flux",
            self.body
                .flux(layer, self.faces[layer], self.generated[layer], offset_m),
        )
    }

    fn check_query(&self, cx: &Cx<'_>, index: usize, offset: f64) -> Result<(), ConductionError> {
        checkpoint(cx, index)?;
        let layer = self
            .body
            .wall
            .layers
            .get(index)
            .ok_or_else(|| invalid("solid profile layer index is out of range"))?;
        if !offset.is_finite() || offset < 0.0 || offset > layer.thickness_m {
            return Err(invalid(
                "solid profile query must lie inside its physical layer",
            ));
        }
        Ok(())
    }
}

fn admitted(temperature: f64, interval: [f64; 2]) -> Result<(), ConductionError> {
    if !temperature.is_finite() || temperature < interval[0] || temperature > interval[1] {
        return Err(invalid(
            "solid radial face temperature leaves the admitted material interval",
        ));
    }
    Ok(())
}
