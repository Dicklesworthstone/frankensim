//! A zero-storage surface shared by liquid conduction, gas film and radiation.
//!
//! The finite-volume liquid half-cell uses its current bulk conductivity. The
//! gas conductivity is resolved at the solved mean surface/bulk-gas temperature.
//! Eliminating that surface gives three positive pairwise conductances, so the
//! existing thermal owner can conserve heat and bound BOTH finite capacities.

use super::*;
use crate::conjugate::radiation::STEFAN_BOLTZMANN;
use fs_matdb::{PropertyValue, QueryPoint};
use fs_material::state_point::ResolvedScalarProperty;
use fs_material::thermochemical_gas::energy::{GasEnergySolve, IsochoricGasParcel};
use fs_material::thermochemical_gas::{ResolvedThermochemicalGasState, ThermochemicalGasModel};

mod stress;

/// Explicit geometry of one liquid/gas interface facing black surroundings.
/// Both gas-film and liquid-centroid distances are positive physical inputs.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FreeSurfaceInterfaceFace {
    pub face: FreeSurfaceGasFace,
    /// Diffuse-gray opaque-surface emissivity in [0, 1]; zero omits radiation.
    pub emissivity: f64,
    /// Black isothermal surroundings [K], with prescribed view factor one.
    /// This reservoir is independent of the gas and its source temperature range.
    pub surroundings_k: f64,
}

/// One physical patch of the same height graph used by native surface stress.
///
/// The centroid and normal belong to the reconstructed graph, not to a
/// Cartesian voxel face. The positive liquid cell owns this patch's thermal
/// storage; its occupied-liquid centroid is NOT inferred from a fill fraction.
/// Supply its normal conduction distance explicitly when attaching exchange.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FreeSurfaceHeightGraphPatch {
    pub liquid: [usize; 3],
    pub gas: [usize; 3],
    pub centroid_m: [f64; 3],
    pub outward_normal: [f64; 3],
    pub area_m2: f64,
}

impl FreeSurfaceHeightGraphPatch {
    /// Attach the original positive liquid-normal and gas-film distances [m].
    /// Emissivity, temperatures, distances and live geometry are admitted by
    /// the transactional surface-exchange step before any heat is published.
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
                gas: self.gas,
                area_m2: self.area_m2,
                distances_m,
            },
            emissivity,
            surroundings_k,
        }
    }
}

/// Prescribed gas reservoir on a surface with a solved interface temperature.
#[derive(Clone, Copy, Debug)]
pub struct FreeSurfaceInterfaceContact<'a> {
    pub surface: FreeSurfaceInterfaceFace,
    pub model: &'a ThermochemicalGasModel,
    pub context: &'a QueryPoint,
}

/// Bounded nonlinear work and residual controls, per face and sampling epoch.
#[derive(Clone, Copy, Debug)]
pub struct FreeSurfaceInterfaceSolve {
    /// Includes unsupported/backtracked trials; at most 256 attempts are admitted.
    pub maximum_iterations: usize,
    pub temperature_tolerance_k: f64,
    pub absolute_heat_tolerance_w: f64,
    pub relative_heat_tolerance: f64,
}

impl Default for FreeSurfaceInterfaceSolve {
    fn default() -> Self {
        Self {
            maximum_iterations: 64,
            temperature_tolerance_k: 1e-9,
            absolute_heat_tolerance_w: 1e-12,
            relative_heat_tolerance: 1e-10,
        }
    }
}

/// Sources, surface residual and powers actually used by one thermal substep.
#[derive(Clone, Debug, PartialEq)]
pub struct FreeSurfaceInterfaceSample {
    pub surface: FreeSurfaceInterfaceFace,
    /// Original liquid source at the cell temperature, used for its half-cell k.
    pub liquid: ResolvedLiquidState,
    /// Original gas sources at (surface_temperature_k + gas_temperature_k)/2.
    pub gas: ResolvedThermochemicalGasState,
    pub liquid_temperature_k: f64,
    pub surface_temperature_k: f64,
    /// Original sigma source at the solved interface T when surface stress is
    /// active. The bulk liquid bundle retains its actual bulk temperature.
    pub surface_tension: Option<ResolvedScalarProperty>,
    pub gas_temperature_k: f64,
    pub liquid_conductance_w_k: f64,
    pub gas_conductance_w_k: f64,
    /// Exact fourth-power secant, eps sigma A (Ts+Tr)(Ts²+Tr²).
    pub radiation_conductance_w_k: f64,
    /// Sum of the two positive eliminated conductances incident on liquid.
    pub liquid_decay_conductance_w_k: f64,
    /// Includes the gas-to-radiative-reservoir path through the shared surface.
    pub gas_decay_conductance_w_k: f64,
    /// Sum of the following two signed supplies; applied to liquid once.
    pub heat_into_liquid_w: f64,
    /// Positive out of gas. A finite atmosphere stores its exact opposite.
    pub heat_from_gas_w: f64,
    /// Positive out of the radiative reservoir, external to liquid + finite gas.
    pub heat_from_radiation_w: f64,
    /// Raw zero-storage surface residual before conservative elimination [W].
    pub balance_residual_w: f64,
    pub iterations: usize,
}

#[derive(Clone, Debug)]
pub struct FreeSurfaceInterfaceExchange {
    /// Already counted once in the liquid's external account. Its opposite is
    /// finite gas storage when present, making it internal to their union.
    pub heat_from_gas_j: f64,
    /// Already counted once; external also to liquid + finite gas + finite walls.
    pub heat_from_radiation_j: f64,
    pub last_applied: Vec<FreeSurfaceInterfaceSample>,
    pub thermal_endpoint: Vec<FreeSurfaceInterfaceSample>,
    /// Present for temperature-dependent surface stress: re-solved after native
    /// advection using actual bulk state and the same finite gas endpoint.
    /// These diagnostics update pressure/traction only and supply no heat.
    pub native_endpoint: Vec<FreeSurfaceInterfaceSample>,
}

impl MaterialFreeSurfaceEnthalpy3<'_, '_> {
    /// Current complete curved patches from the native height reconstruction.
    ///
    /// Call again after an accepted flow interval. Normals, centroids and areas
    /// use exactly the geometry that supplies Marangoni traction, including
    /// both slope components and either liquid orientation. A patch may span
    /// several cut cells vertically; heat is attached to its outermost positive
    /// liquid cell as a declared column-patch finite-volume reduction. Its
    /// liquid-normal distance remains a physical input, not a guessed half dx.
    /// Native conversion hysteresis can temporarily admit over/underfilled
    /// cells that this monotone graph cannot represent. Such a candidate
    /// refuses transactionally; its mass or geometry is never clipped here.
    pub fn height_graph_interface_patches(
        &self,
        gate: &CancelGate,
    ) -> Result<Vec<FreeSurfaceHeightGraphPatch>, Error> {
        let graph = self
            .height_graph_thermocapillary_surface()
            .ok_or(Error::Invalid(
                "curved heat patches require the original height-graph surface",
            ))?;
        let dx = self.flow.scaling.cell_size_m();
        self.flow
            .lattice
            .height_graph_geometry(graph.native(), gate)?
            .into_iter()
            .map(|patch| {
                poll(gate)?;
                let area_m2 = environment::finite(patch.area * dx * dx)?;
                let centroid_m = patch.centroid.map(|coordinate| coordinate * dx);
                if area_m2 <= 0.0 || centroid_m.iter().any(|coordinate| !coordinate.is_finite()) {
                    return Err(Error::Invalid(
                        "height-graph heat geometry overflows physical units",
                    ));
                }
                Ok(FreeSurfaceHeightGraphPatch {
                    liquid: patch.liquid,
                    gas: patch.gas,
                    centroid_m,
                    outward_normal: patch.outward_normal,
                    area_m2,
                })
            })
            .collect()
    }

    /// Replace the projected outer Cartesian face by its actual graph area in
    /// the original thermal owner's exposed-area budget. Values are in m².
    pub(super) fn height_graph_area_corrections(
        &self,
        gate: &CancelGate,
    ) -> Result<BTreeMap<[usize; 3], f64>, Error> {
        if self.height_graph_thermocapillary_surface().is_none() {
            return Ok(BTreeMap::new());
        }
        let dx = self.flow.scaling.cell_size_m();
        self.height_graph_interface_patches(gate)?
            .into_iter()
            .map(|patch| {
                let extra = environment::finite(patch.area_m2 - dx * dx)?;
                if extra < 0.0 {
                    return Err(Error::Invalid(
                        "graph area is smaller than its projected area",
                    ));
                }
                Ok((patch.liquid, extra))
            })
            .collect()
    }

    pub(super) fn admit_surfaces(
        &self,
        surfaces: &[FreeSurfaceInterfaceContact<'_>],
        solve: FreeSurfaceInterfaceSolve,
        cell_id: &impl Fn([usize; 3]) -> Result<CellId, Error>,
        dx: f64,
        gate: &CancelGate,
    ) -> Result<Vec<CellId>, Error> {
        if solve.maximum_iterations == 0
            || solve.maximum_iterations > 256
            || !solve.temperature_tolerance_k.is_finite()
            || solve.temperature_tolerance_k <= 0.0
            || !solve.absolute_heat_tolerance_w.is_finite()
            || solve.absolute_heat_tolerance_w < 0.0
            || !solve.relative_heat_tolerance.is_finite()
            || solve.relative_heat_tolerance < 0.0
        {
            return Err(Error::Invalid("invalid surface-temperature solve budget"));
        }
        if !surfaces.is_empty()
            && (self.height_graph_thermocapillary_surface().is_some()
                || self.temperature_dependent_thermocapillarity())
        {
            self.admit_surface_stress_geometry(surfaces, dx, gate)?;
        }
        if !surfaces.is_empty() && self.height_graph_thermocapillary_surface().is_some() {
            // The shared graph reconstruction has already checked each live
            // liquid/gas owner, unique column and actual curved patch area.
            // Cartesian face-area limits remain on the separate voxel APIs.
            let mut endpoints = Vec::with_capacity(surfaces.len());
            for contact in surfaces {
                poll(gate)?;
                let surface = contact.surface;
                if surface
                    .face
                    .distances_m
                    .iter()
                    .any(|d| !d.is_finite() || *d <= 0.0)
                    || !surface.emissivity.is_finite()
                    || !(0.0..=1.0).contains(&surface.emissivity)
                    || !surface.surroundings_k.is_finite()
                    || surface.surroundings_k <= 0.0
                {
                    return Err(Error::Invalid(
                        "curved surface exchange needs positive normal distances and kelvin surroundings with supported emissivity",
                    ));
                }
                contact
                    .model
                    .evaluate(contact.context)
                    .map_err(MaterialLbmError::from)?;
                endpoints.push(cell_id(surface.face.liquid)?);
            }
            return Ok(endpoints);
        }
        let mut gas = Vec::with_capacity(surfaces.len());
        let mut radiation = Vec::with_capacity(surfaces.len());
        for contact in surfaces {
            poll(gate)?;
            let surface = contact.surface;
            gas.push(FreeSurfaceGasContact {
                face: surface.face,
                model: contact.model,
                context: contact.context,
            });
            radiation.push(FreeSurfaceRadiationFace {
                liquid: surface.face.liquid,
                gas: surface.face.gas,
                area_m2: surface.face.area_m2,
                emissivity: surface.emissivity,
                surroundings_k: surface.surroundings_k,
            });
        }
        self.admit_radiation(&radiation, cell_id, dx, gate)?;
        self.admit_environment(&gas, cell_id, dx, gate)
    }

    /// Solve a shared liquid/gas/radiative surface on each currently exposed face,
    /// then advance the original thermal, wall, rheology and native VOF owners.
    ///
    /// Bulk liquid k and declared centroid distance form a finite-volume
    /// half-cell resistance. Gas k is reevaluated at the solved film temperature;
    /// radiation uses that same surface temperature. The surface has no heat
    /// capacity, mass transfer or phase change. It must stay inside the current
    /// pure-liquid caloric interval. Gas claims are pinned to the admitted bath
    /// sources; no source switch or extrapolation supplies an interface root.
    /// A safeguarded scalar iteration must meet BOTH declared tolerances.
    /// For temperature-dependent planar or height-graph capillarity, supply
    /// exactly one full outward patch per tank column. Curved patch owners and
    /// actual areas come from [`Self::height_graph_interface_patches`]; the
    /// liquid-normal conduction distance remains explicit. Its original sigma(Ts) source drives BOTH
    /// Marangoni traction and native capillary pressure. The thermal endpoint
    /// starts the flow split; actual advected bulk states re-solve the surface
    /// before acceptance, without supplying extra heat. Graph owners, positions,
    /// normals and areas are reconstructed at that native endpoint. Positive
    /// film/centroid normal distances are held over this interval and must be
    /// supplied again on continuation with current complete geometry.
    /// Once active, steps omitting this complete boundary refuse transactionally.
    /// Bulk viscosity, density and caloric storage retain their actual bulk state.
    ///
    /// The eliminated positive conductances share the original explicit CFL
    /// and timeline. Time/mesh refinement remains required for this first-order,
    /// frozen-coefficient split. Geometry is supplied anew after native motion.
    /// The older resolved-cell radiation API retains its explicit cell-T model.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_surface_exchange(
        &mut self,
        faces: &[FreeSurfaceConductionFace],
        heat: &[FreeSurfaceCellHeat],
        controls: FreeSurfaceConductionControls,
        bodies: &mut [FreeSurfaceThermalBody<'_>],
        contacts: &[FreeSurfaceBodyContact],
        interfaces: &[FreeSurfaceBodyInterface<'_>],
        surfaces: &[FreeSurfaceInterfaceContact<'_>],
        solve: FreeSurfaceInterfaceSolve,
        with_viscous_heating: bool,
        gate: &CancelGate,
    ) -> Result<viscous::FreeSurfaceViscousConductionReport, Error> {
        let applied = if with_viscous_heating {
            self.viscous_sources(gate)?
        } else {
            Vec::new()
        };
        let conduction = self.advance_with_thermal_boundaries(
            faces,
            heat,
            controls,
            &applied,
            bodies,
            contacts,
            interfaces,
            &[],
            &[],
            surfaces,
            solve,
            None,
            gate,
        )?;
        Ok(viscous::FreeSurfaceViscousConductionReport {
            conduction,
            applied,
        })
    }

    /// The same solved surface with a well-mixed, fixed-mass/volume atmosphere.
    /// Gas stores the exact opposite gas supply on EACH thermal substep. Updated
    /// NASA-9 internal energy, T and EOS pressure feed the next surface solve.
    /// Radiation may heat gas through the surface even when gas and liquid bulk
    /// temperatures coincide. Its contribution to gas decay bounds the step.
    ///
    /// No gas momentum, pressure work, evaporation, participating-media radiation
    /// or moving-solid/phase-density handoff is added. Liquid, wall and gas
    /// histories publish together only after every endpoint/native gate passes.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_finite_gas_surface_exchange(
        &mut self,
        faces: &[FreeSurfaceConductionFace],
        heat: &[FreeSurfaceCellHeat],
        controls: FreeSurfaceConductionControls,
        bodies: &mut [FreeSurfaceThermalBody<'_>],
        contacts: &[FreeSurfaceBodyContact],
        interfaces: &[FreeSurfaceBodyInterface<'_>],
        gas: &mut IsochoricGasParcel,
        surfaces: &[FreeSurfaceInterfaceFace],
        gas_controls: GasEnergySolve,
        solve: FreeSurfaceInterfaceSolve,
        with_viscous_heating: bool,
        gate: &CancelGate,
    ) -> Result<viscous::FreeSurfaceViscousConductionReport, Error> {
        poll(gate)?;
        if surfaces.is_empty() || surfaces.len() > controls.maximum_faces {
            return Err(Error::Invalid(
                "finite atmosphere requires bounded exposed liquid patches",
            ));
        }
        let mut candidate = atmosphere::FiniteAtmosphere::try_new(gas, gas_controls)?;
        let boundaries: Vec<_> = surfaces
            .iter()
            .map(|&surface| FreeSurfaceInterfaceContact {
                surface,
                model: gas.model(),
                context: gas.point(),
            })
            .collect();
        let applied = if with_viscous_heating {
            self.viscous_sources(gate)?
        } else {
            Vec::new()
        };
        let conduction = self.advance_with_thermal_boundaries(
            faces,
            heat,
            controls,
            &applied,
            bodies,
            contacts,
            interfaces,
            &[],
            &[],
            &boundaries,
            solve,
            Some(&mut candidate),
            gate,
        )?;
        *gas = candidate.parcel;
        Ok(viscous::FreeSurfaceViscousConductionReport {
            conduction,
            applied,
        })
    }
}

// A connected source interval containing the seed; touching NASA regions and
// liquid caloric segments join, but a gap never becomes invented support.
fn containing_interval(ranges: impl IntoIterator<Item = [f64; 2]>, seed: f64) -> Option<[f64; 2]> {
    let mut run: Option<[f64; 2]> = None;
    for [lo, hi] in ranges {
        if let Some([start, end]) = run {
            if lo <= end {
                run = Some([start, hi.max(end)]);
            } else if seed >= start && seed <= end {
                return run;
            } else {
                run = Some([lo, hi]);
            }
        } else {
            run = Some([lo, hi]);
        }
    }
    run.filter(|[lo, hi]| seed >= *lo && seed <= *hi)
}

// This bounds trials, not source evaluation: every trial still goes through
// the original complete evaluator. In particular a hot radiative bath is NEVER
// queried as gas; its temperature need not be inside any gas material domain.
fn film_interval(
    model: &ThermochemicalGasModel,
    gas: &ResolvedThermochemicalGasState,
) -> Result<[f64; 2], Error> {
    let mut range = [0.0_f64, f64::MAX];
    for species in model.thermodynamics().components() {
        let [lo, hi] = containing_interval(
            species
                .regions()
                .iter()
                .map(|r| [r.temperature_min().value(), r.temperature_max().value()]),
            gas.state().temperature,
        )
        .ok_or(Error::Invalid(
            "surface film lacks connected gas caloric support",
        ))?;
        range[0] = range[0].max(lo);
        range[1] = range[1].min(hi);
    }
    for (material, component) in model.transport_materials().iter().zip(gas.components()) {
        for property in component.properties() {
            let answer = property.answer();
            if let Some((lo, hi)) = answer.evidence.model.validity.bound("temperature") {
                range[0] = range[0].max(lo);
                range[1] = range[1].min(hi);
            }
            let claim = material
                .material
                .claims()
                .claim(answer.receipt.selected)
                .ok_or(Error::Invalid("surface film source claim is absent"))?;
            if let PropertyValue::Curve {
                abscissa, knots, ..
            } = &claim.value
            {
                if abscissa == "temperature" {
                    range[0] = range[0].max(knots[0].0);
                    range[1] = range[1].min(knots[knots.len() - 1].0);
                }
            }
        }
    }
    Ok(range)
}

impl Conduction<'_, '_> {
    pub(super) fn sample_surfaces(
        &self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Vec<FreeSurfaceInterfaceSample>, Error> {
        let mut samples = Vec::with_capacity(self.surfaces.len());
        let mut covered = BTreeMap::new();
        for (index, (contact, &id)) in self
            .surfaces
            .iter()
            .zip(&self.surface_endpoints)
            .enumerate()
        {
            check()?;
            let surface = contact.surface;
            let area = covered.entry(id).or_insert(0.0);
            *area = environment::finite(*area + surface.face.area_m2)?;
            // Curved area replaces a projected outer face in the cell budget.
            // Subtracting internal faces can lose a few ulps; permit only that
            // bounded arithmetic error, never additional physical coverage.
            let available = thermal.exposed_area_m2(id)?;
            if *area - available > 64.0 * f64::EPSILON * area.abs().max(available.abs()) {
                return Err(Error::Invalid("solved surfaces exceed exposed liquid area"));
            }
            let liquid = self.sample_liquid(id, thermal.cells()[&id].specific_h_j_kg, check)?;
            let tl = thermal.temperature_k(id)?;
            samples.push(self.solve_surface(index, liquid, tl, check)?);
        }
        Ok(samples)
    }
    // The same scalar/source owner serves thermal substeps and native endpoint
    // diagnostics. Endpoint calls do not integrate either side's energy.
    fn solve_surface(
        &self,
        index: usize,
        liquid: ResolvedLiquidState,
        tl: f64,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<FreeSurfaceInterfaceSample, Error> {
        self.solve_surface_on(index, self.surfaces[index].surface, liquid, tl, check)
    }

    // Native graph motion can change the exposed cell and area. Reuse exactly
    // the original scalar/source solve on that live patch, with the declared
    // normal distances held for this split. This diagnostic supplies no heat.
    fn solve_surface_on(
        &self,
        index: usize,
        surface: FreeSurfaceInterfaceFace,
        liquid: ResolvedLiquidState,
        tl: f64,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<FreeSurfaceInterfaceSample, Error> {
        let contact = &self.surfaces[index];
        let (model, context) = self
            .finite_atmosphere
            .as_ref()
            .map_or((contact.model, contact.context), |gas| {
                (gas.parcel.model(), gas.parcel.point())
            });
        let bulk = if let Some(atmosphere) = &self.finite_atmosphere {
            // The original parcel has already resolved this exact T/p;
            // reuse its accepted source bundle for every exposed patch.
            atmosphere.parcel.state().clone()
        } else {
            model.evaluate(context).map_err(MaterialLbmError::from)?
        };
        if self
            .surface_reference
            .get(index)
            .is_some_and(|old| !environment::same_sources(&old.gas, &bulk))
        {
            return Err(Error::Invalid(
                "surface film left its original gas source claims",
            ));
        }
        let tg = bulk.state().temperature;
        let tr = surface.surroundings_k;
        let [film_lo, film_hi] = film_interval(model, &bulk)?;
        let [liquid_lo, liquid_hi] = containing_interval(
            self.caloric
                .knots()
                .windows(2)
                .filter(|pair| {
                    pair.iter().all(|k| {
                        k.liquid_mass_fraction == 1.0
                            && k.bulk_density_kg_m3.to_bits() == self.density.to_bits()
                    })
                })
                .map(|pair| [pair[0].temperature_k, pair[1].temperature_k]),
            tl,
        )
        .ok_or(Error::Invalid(
            "surface temperature needs a supported pure-liquid caloric interval",
        ))?;
        let mut lo = tl.min(tg).min(tr).max(liquid_lo).max(2.0 * film_lo - tg);
        let mut hi = tl.max(tg).max(tr).min(liquid_hi).min(2.0 * film_hi - tg);
        if !lo.is_finite() || !hi.is_finite() || lo > hi {
            return Err(Error::Invalid(
                "surface and film temperature domains do not intersect",
            ));
        }
        let gl = environment::finite(
            liquid.thermal_conductivity_w_m_k() * surface.face.area_m2
                / surface.face.distances_m[0],
        )?;
        let radiation = surface.emissivity * STEFAN_BOLTZMANN * surface.face.area_m2;
        let mut ts = tl.clamp(lo, hi);
        let mut supported = None;
        let mut accepted = None;
        for iteration in 1..=self.surface_solve.maximum_iterations {
            check()?;
            let point = context
                .clone()
                .with("temperature", 0.5 * ts + 0.5 * tg)
                .map_err(|_| Error::Invalid("invalid solved gas-film temperature"))?;
            let gas = match model.evaluate(&point) {
                Ok(gas) if environment::same_sources(&bulk, &gas) => gas,
                _ => {
                    // A triangulated claim can have less support than its
                    // validity box. Backtrack a trial toward a supported
                    // sample; rejected trials consume the SAME solve budget.
                    // Bulk gas need not lie in the pure-liquid surface
                    // interval. Move toward its admissible endpoint, not
                    // an out-of-phase temperature that rejects a valid root.
                    let anchor = supported.unwrap_or(tg.clamp(lo, hi));
                    if anchor < lo || anchor > hi || anchor == ts {
                        return Err(Error::Invalid(
                            "surface solve has no supported gas-film trial",
                        ));
                    }
                    if ts > anchor {
                        hi = ts;
                    } else {
                        lo = ts;
                    }
                    ts = 0.5 * ts + 0.5 * anchor;
                    continue;
                }
            };
            supported = Some(ts);
            let gg = environment::finite(
                gas.state().thermal_conductivity * surface.face.area_m2
                    / surface.face.distances_m[1],
            )?;
            let gr = if radiation == 0.0 {
                0.0
            } else {
                environment::finite(radiation * (ts + tr) * (ts * ts + tr * tr))?
            };
            let total = environment::finite(gl + gg + gr)?;
            let ql = environment::finite(gl * (tl - ts))?;
            let qg = environment::finite(gg * (tg - ts))?;
            let qr = environment::finite(gr * (tr - ts))?;
            let residual = environment::finite(ql + qg + qr)?;
            let correction = if total == 0.0 { 0.0 } else { residual / total };
            let tolerance = environment::finite(
                self.surface_solve.absolute_heat_tolerance_w
                    + self.surface_solve.relative_heat_tolerance * (ql.abs() + qg.abs() + qr.abs()),
            )?;
            if residual.abs() <= tolerance
                && correction.abs() <= self.surface_solve.temperature_tolerance_k
            {
                // Schur elimination conserves every pair exactly even when
                // the scalar residual is finite. Its effect on interface T
                // and each original flux is bounded by the admitted solve.
                let (lg, lr, gr_pair) = if total == 0.0 {
                    (0.0, 0.0, 0.0)
                } else {
                    (gl * (gg / total), gl * (gr / total), gg * (gr / total))
                };
                let gas_heat = environment::finite(lg * (tg - tl) + gr_pair * (tg - tr))?;
                let radiant_heat = environment::finite(lr * (tr - tl) + gr_pair * (tr - tg))?;
                let surface_tension = if self.surface_stress {
                    Some(stress::resolve_surface_tension(self.law, &liquid, ts)?)
                } else {
                    None
                };
                accepted = Some(FreeSurfaceInterfaceSample {
                    surface,
                    liquid,
                    gas,
                    liquid_temperature_k: tl,
                    surface_temperature_k: ts,
                    surface_tension,
                    gas_temperature_k: tg,
                    liquid_conductance_w_k: gl,
                    gas_conductance_w_k: gg,
                    radiation_conductance_w_k: gr,
                    liquid_decay_conductance_w_k: environment::finite(lg + lr)?,
                    gas_decay_conductance_w_k: environment::finite(lg + gr_pair)?,
                    heat_into_liquid_w: environment::finite(gas_heat + radiant_heat)?,
                    heat_from_gas_w: gas_heat,
                    heat_from_radiation_w: radiant_heat,
                    balance_residual_w: residual,
                    iterations: iteration,
                });
                break;
            }
            if residual > 0.0 {
                lo = ts;
            } else {
                hi = ts;
            }
            let slope = environment::finite(gl + gg + 4.0 * radiation * ts * ts * ts)?;
            let newton = ts + residual / slope;
            // Periodic bisection guarantees bounded bracket progress even
            // when temperature-dependent gas k makes this slope inexact.
            let next = if iteration % 4 != 0 && newton > lo && newton < hi {
                newton
            } else {
                0.5 * lo + 0.5 * hi
            };
            if next == ts {
                break;
            }
            ts = next;
        }
        accepted.ok_or(Error::Invalid(
            "surface temperature exhausted its supported solve budget",
        ))
    }
}
