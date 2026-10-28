//! The solved heat interface supplies the SAME planar pressure and traction.
//!
//! Each tank column owns one surface state. Its sigma extends normally onto
//! native wet/zero-mass ghost lanes; it never replaces their bulk material T.

use super::*;
use crate::d3q19::TILE;
use crate::material::free_surface::{MaterialFreeSurfaceScaling, capillary_scaling};
use fs_material::state_point::{MaterialPropertySelection, resolve_material_state_point};

fn column(p: [usize; 3], normal: usize) -> [usize; 2] {
    match normal {
        0 => [p[1], p[2]],
        1 => [p[0], p[2]],
        _ => [p[0], p[1]],
    }
}

/// Query only the original surface-tension claim at interface T. Keeping the
/// bulk state's other coordinates does not re-evaluate or overwrite viscosity,
/// density, caloric enthalpy, mass, or the held rheological bulk state.
pub(super) fn resolve_surface_tension(
    law: &LiquidTemperatureLaw<'_>,
    liquid: &ResolvedLiquidState,
    temperature_k: f64,
) -> Result<ResolvedScalarProperty, Error> {
    let original = law
        .reference_state()
        .material()
        .property(SURFACE_TENSION_PROPERTY)
        .ok_or(Error::Invalid(
            "solved surface stress requires the original sigma claim",
        ))?;
    let point = liquid
        .material()
        .query()
        .with(law.temperature_axis(), temperature_k)
        .map_err(|_| Error::Invalid("invalid surface-tension temperature coordinate"))?;
    let state = resolve_material_state_point(
        law.card(),
        &point,
        std::slice::from_ref(original.requirement()),
        MaterialPropertySelection::PinnedByProperty(vec![(
            SURFACE_TENSION_PROPERTY.into(),
            original.answer().receipt.selected,
        )]),
    )?;
    let sigma = state
        .property(SURFACE_TENSION_PROPERTY)
        .ok_or(Error::Invalid(
            "solved surface stress lost its sigma source",
        ))?;
    if !sigma.value_si().is_finite() || sigma.value_si() <= 0.0 {
        return Err(Error::Invalid(
            "solved surface stress needs positive finite sigma",
        ));
    }
    Ok(sigma.clone())
}

impl MaterialFreeSurfaceEnthalpy3<'_, '_> {
    pub(super) fn admit_surface_stress_geometry(
        &self,
        surfaces: &[FreeSurfaceInterfaceContact<'_>],
        dx: f64,
        gate: &CancelGate,
    ) -> Result<(), Error> {
        let plane = self.thermocapillary_surface().ok_or(Error::Invalid(
            "solved surface stress requires the original planar geometry",
        ))?;
        let axis = plane.normal_axis;
        if axis >= 3 || ![-1, 1].contains(&plane.outward_sign) {
            return Err(Error::Invalid("invalid solved-surface normal"));
        }
        let (nx, ny, nz) = self.flow.lattice.grid().tile_dims();
        let size = [nx * TILE, ny * TILE, nz * TILE];
        let tangent = column(size, axis);
        let count = tangent[0]
            .checked_mul(tangent[1])
            .ok_or(Error::Invalid("surface column count overflows"))?;
        if surfaces.len() != count {
            return Err(Error::Invalid(
                "surface-temperature stress needs one complete outward heat patch per tank column",
            ));
        }
        let position = plane.surface_position_m / dx;
        let layer = if plane.outward_sign > 0 {
            position.ceil() as usize - 1
        } else {
            position.floor() as usize
        };
        let mut columns = BTreeSet::new();
        for contact in surfaces {
            poll(gate)?;
            let face = contact.surface.face;
            let mut outward = face.liquid;
            outward[axis] = if plane.outward_sign > 0 {
                layer.checked_add(1)
            } else {
                layer.checked_sub(1)
            }
            .ok_or(Error::Invalid("surface gas neighbor is outside the tank"))?;
            if face.liquid[axis] != layer
                || face.gas != outward
                || face.area_m2 != dx * dx
                || !columns.insert(column(face.liquid, axis))
            {
                return Err(Error::Invalid(
                    "surface-temperature stress needs one full normal patch on each surface cell",
                ));
            }
        }
        Ok(())
    }
}

impl Conduction<'_, '_> {
    /// Same material/scalar closure at the actual native epoch. The first
    /// invocation consumes the accepted thermal endpoint; the second resolves
    /// advected h/rheology with the held boundary geometry and finite gas state.
    /// Neither invocation integrates another thermal interval.
    #[allow(clippy::too_many_arguments)]
    pub(in super::super) fn surface_force_field(
        &mut self,
        _lattice: &FreeSurface3,
        states: &BTreeMap<MassLocation3, ResolvedLiquidState>,
        normal_axis: usize,
        scaling: MaterialFreeSurfaceScaling<'_>,
        prepared: bool,
        gate: &CancelGate,
    ) -> Result<BTreeMap<(u64, usize), f64>, Error> {
        let samples = if prepared {
            self.surface_endpoint.clone()
        } else {
            let mut samples = Vec::with_capacity(self.surfaces.len());
            for (index, contact) in self.surfaces.iter().enumerate() {
                poll(gate)?;
                let location = FreeSurface3::transport_location(contact.surface.face.liquid)
                    .ok_or(Error::Invalid("unrepresentable solved-surface cell"))?;
                let liquid = states.get(&location).ok_or(Error::Invalid(
                    "solved surface left its admitted liquid cell",
                ))?;
                let point = liquid.material().query();
                let temperature = *point
                    .axes()
                    .get(self.law.temperature_axis())
                    .ok_or(Error::Invalid("solved surface lost its bulk temperature"))?;
                samples.push(self.solve_surface(
                    index,
                    liquid.clone(),
                    temperature,
                    &mut || checkpoint(gate),
                )?);
            }
            samples
        };
        if samples.len() != self.surfaces.len() {
            return Err(Error::Invalid(
                "solved surface temperature field is incomplete",
            ));
        }
        let mut columns = BTreeMap::new();
        for sample in &samples {
            poll(gate)?;
            let sigma = sample.surface_tension.as_ref().ok_or(Error::Invalid(
                "surface temperature did not resolve its original sigma",
            ))?;
            let native = capillary_scaling(
                Some(sigma.value_si()),
                scaling.pressure_unit_pa,
                scaling.cell_size_m(),
            )?;
            if columns
                .insert(column(sample.surface.face.liquid, normal_axis), native)
                .is_some()
            {
                return Err(Error::Invalid("duplicate solved-surface stress column"));
            }
        }
        let mut field = BTreeMap::new();
        for &location in states.keys() {
            poll(gate)?;
            let MassLocation3::Cell(key, lane) = location else {
                continue;
            };
            let p = FreeSurface3::transport_coordinates(location).ok_or(Error::Invalid(
                "solved surface stress has an unlocated wet cell",
            ))?;
            let sigma = *columns.get(&column(p, normal_axis)).ok_or(Error::Invalid(
                "solved surface stress does not cover every wet column",
            ))?;
            field.insert((key, lane), sigma);
        }
        if !prepared {
            self.surface_native_endpoint = samples;
        }
        Ok(field)
    }
}
