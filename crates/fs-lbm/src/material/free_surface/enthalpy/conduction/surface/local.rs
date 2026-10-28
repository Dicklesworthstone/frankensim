//! Solved interface temperature on the native detached-interface PLIC patches.

use super::*;
use crate::material::free_surface::enthalpy::thermocapillary::LocalThermocapillarity;
use crate::material::free_surface::{MaterialFreeSurfaceScaling, capillary_scaling};

/// How prescribed boundary data continue when detached PLIC owners change.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LocalSurfaceBoundaryEvolution {
    /// Each current partial-cell owner retains its supplied boundary. New or
    /// vanished patch owners refuse; no spatial interpolation is invented.
    FixedOwners,
    /// One explicitly spatially uniform, prescribed reservoir covers every
    /// reconstructed patch, including newly exposed owners. All model/context,
    /// distance and radiation inputs must be identical. Finite gas refuses.
    UniformReservoir,
}

/// Current local PLIC patch and its actual occupied-liquid boundary area.
///
/// The interface separates liquid and gas subvolumes of the SAME partial cell;
/// it does not require an adjacent wholly empty Cartesian cell. Positive normal
/// conduction distances remain physical inputs, not inferred from voxel fill.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FreeSurfaceLocalPatch {
    pub liquid: [usize; 3],
    pub centroid_m: [f64; 3],
    pub outward_normal: [f64; 3],
    pub area_m2: f64,
    /// Cut surface plus all clipped wet Cartesian faces, for the original
    /// thermal owner's total-area budget. It does not declare thermal links.
    pub liquid_boundary_area_m2: f64,
}

impl FreeSurfaceLocalPatch {
    /// Attach a bath to this resolved partial cell. Only the Local solved-heat
    /// entry admits equal liquid/gas coordinates for its distinct subvolumes.
    #[must_use]
    pub fn with_exchange(
        self,
        distances_m: [f64; 2],
        emissivity: f64,
        surroundings_k: f64,
    ) -> FreeSurfaceInterfaceFace {
        FreeSurfaceInterfaceFace {
            face: FreeSurfaceGasFace {
                liquid: self.liquid,
                gas: self.liquid,
                area_m2: self.area_m2,
                distances_m,
            },
            emissivity,
            surroundings_k,
        }
    }
}

impl MaterialFreeSurfaceEnthalpy3<'_, '_> {
    /// Actual partial-cell patches from the same native reconstruction used for
    /// local Marangoni stress. Empty/full closure ghosts own no patch or heat.
    pub fn local_interface_patches(
        &self,
        gate: &CancelGate,
    ) -> Result<Vec<FreeSurfaceLocalPatch>, Error> {
        let model = self.local_thermocapillary_surface().ok_or(Error::Invalid(
            "local heat patches require the original detached-interface geometry",
        ))?;
        let dx = self.flow.scaling.cell_size_m();
        self.flow
            .lattice
            .local_interface_geometry(model.native(), gate)?
            .into_iter()
            .map(|patch| {
                poll(gate)?;
                let area_m2 = environment::finite(patch.area * dx * dx)?;
                let liquid_boundary_area_m2 =
                    environment::finite(patch.liquid_boundary_area() * dx * dx)?;
                let centroid_m = patch.centroid.map(|v| v * dx);
                if area_m2 <= 0.0
                    || liquid_boundary_area_m2 < area_m2
                    || centroid_m.iter().any(|v| !v.is_finite())
                {
                    return Err(Error::Invalid(
                        "local heat geometry overflows physical units",
                    ));
                }
                Ok(FreeSurfaceLocalPatch {
                    liquid: patch.liquid,
                    centroid_m,
                    outward_normal: patch.outward_normal,
                    area_m2,
                    liquid_boundary_area_m2,
                })
            })
            .collect()
    }

    pub(super) fn admit_local_surface_geometry(
        &self,
        surfaces: &[FreeSurfaceInterfaceContact<'_>],
        gate: &CancelGate,
    ) -> Result<(), Error> {
        let patches: BTreeMap<_, _> = self
            .local_interface_patches(gate)?
            .into_iter()
            .map(|p| (p.liquid, p))
            .collect();
        if patches.len() != surfaces.len() {
            return Err(Error::Invalid(
                "local surface temperature needs every resolved partial-cell patch",
            ));
        }
        let mut covered = BTreeSet::new();
        for contact in surfaces {
            poll(gate)?;
            let face = contact.surface.face;
            let patch = patches
                .get(&face.liquid)
                .ok_or(Error::Invalid("local heat patch has no live native owner"))?;
            if face.gas != face.liquid
                || !stress::same_area(face.area_m2, patch.area_m2)
                || !covered.insert(face.liquid)
            {
                return Err(Error::Invalid(
                    "local heat patch must uniquely match its cut-cell liquid/gas owner and actual area",
                ));
            }
        }
        Ok(())
    }

    pub(super) fn admit_local_boundary_evolution(
        &self,
        surfaces: &[FreeSurfaceInterfaceContact<'_>],
        evolution: LocalSurfaceBoundaryEvolution,
    ) -> Result<(), Error> {
        if evolution == LocalSurfaceBoundaryEvolution::FixedOwners {
            return Ok(());
        }
        let Some(first) = surfaces.first() else {
            return Err(Error::Invalid(
                "uniform local reservoir needs a complete nonempty patch field",
            ));
        };
        if self.local_thermocapillary_surface().is_none()
            || surfaces.iter().any(|contact| {
                !core::ptr::eq(contact.model, first.model)
                    || contact.context != first.context
                    || contact.surface.face.distances_m != first.surface.face.distances_m
                    || contact.surface.emissivity != first.surface.emissivity
                    || contact.surface.surroundings_k != first.surface.surroundings_k
            })
        {
            return Err(Error::Invalid(
                "uniform local reservoir requires one identical source binding, complete context and physical boundary everywhere",
            ));
        }
        Ok(())
    }
}

impl Conduction<'_, '_> {
    pub(super) fn local_surface_force_field(
        &mut self,
        lattice: &FreeSurface3,
        states: &BTreeMap<MassLocation3, ResolvedLiquidState>,
        model: LocalThermocapillarity,
        scaling: MaterialFreeSurfaceScaling<'_>,
        prepared: bool,
        gate: &CancelGate,
    ) -> Result<BTreeMap<(u64, usize), f64>, Error> {
        let patches: BTreeMap<_, _> = lattice
            .local_interface_geometry(model.native(), gate)?
            .into_iter()
            .map(|p| (p.liquid, p))
            .collect();
        if patches.len() > self.maximum_surface_patches {
            return Err(Error::Invalid(
                "evolving local heat geometry exceeds the declared thermal face budget",
            ));
        }
        let uniform = self.surface_solve.local_boundary_evolution
            == LocalSurfaceBoundaryEvolution::UniformReservoir;
        if uniform && self.finite_atmosphere.is_some() {
            return Err(Error::Invalid(
                "uniform-reservoir evolution cannot remap a finite gas boundary",
            ));
        }
        let dx2 = scaling.cell_size_m().powi(2);
        let samples = if prepared {
            if self.surface_endpoint.len() != patches.len() {
                return Err(Error::Invalid(
                    "prepared local surface-temperature field is incomplete",
                ));
            }
            for sample in &self.surface_endpoint {
                let face = sample.surface.face;
                let patch = patches.get(&face.liquid).ok_or(Error::Invalid(
                    "prepared local surface temperature has a stale owner",
                ))?;
                if face.gas != face.liquid || !stress::same_area(face.area_m2, patch.area * dx2) {
                    return Err(Error::Invalid(
                        "prepared local surface temperature has stale geometry",
                    ));
                }
            }
            self.surface_endpoint.clone()
        } else {
            let contacts: BTreeMap<_, _> = self
                .surfaces
                .iter()
                .enumerate()
                .map(|(i, contact)| (contact.surface.face.liquid, i))
                .collect();
            if !uniform && contacts.len() != patches.len() {
                return Err(Error::Invalid(
                    "fixed local heat owners changed during native flow",
                ));
            }
            let mut samples = Vec::with_capacity(patches.len());
            for (&owner, patch) in &patches {
                poll(gate)?;
                let index = if uniform {
                    0
                } else {
                    *contacts.get(&owner).ok_or(Error::Invalid(
                        "new local patch has no declared heat boundary",
                    ))?
                };
                let mut live = self.surfaces[index].surface;
                live.face.liquid = owner;
                live.face.gas = owner;
                live.face.area_m2 = environment::finite(patch.area * dx2)?;
                let location = MassLocation3::Cell(patch.location.0, patch.location.1);
                let liquid = states.get(&location).ok_or(Error::Invalid(
                    "live local heat patch lost its original liquid material state",
                ))?;
                let temperature = *liquid
                    .material()
                    .query()
                    .axes()
                    .get(self.law.temperature_axis())
                    .ok_or(Error::Invalid("local surface lost its bulk temperature"))?;
                samples.push(self.solve_surface_on(
                    index,
                    live,
                    liquid.clone(),
                    temperature,
                    &mut || checkpoint(gate),
                )?);
            }
            samples
        };
        // Complete native fields retain actual bulk values away from patches.
        // The dedicated local surface fit below samples ONLY patch centroids.
        let mut field = BTreeMap::new();
        for (&location, state) in states {
            poll(gate)?;
            let MassLocation3::Cell(key, lane) = location else {
                continue;
            };
            let sigma =
                state
                    .material()
                    .property(SURFACE_TENSION_PROPERTY)
                    .ok_or(Error::Invalid(
                        "local surface field lost its original bulk sigma",
                    ))?;
            field.insert(
                (key, lane),
                capillary_scaling(
                    Some(sigma.value_si()),
                    scaling.pressure_unit_pa,
                    scaling.cell_size_m(),
                )?,
            );
        }
        for sample in &samples {
            poll(gate)?;
            let sigma = sample.surface_tension.as_ref().ok_or(Error::Invalid(
                "local interface temperature did not resolve the original sigma source",
            ))?;
            let patch = &patches[&sample.surface.face.liquid];
            field.insert(
                patch.location,
                capillary_scaling(
                    Some(sigma.value_si()),
                    scaling.pressure_unit_pa,
                    scaling.cell_size_m(),
                )?,
            );
        }
        if !prepared {
            self.surface_native_endpoint = samples;
        }
        Ok(field)
    }
}
