//! An admitted planar thermal expansion step into the original native liquid.
//!
//! The pressure owner determines the required normal impulse. Exact affine
//! material overlaps move mass, uniform enthalpy and momentum into actual VOF
//! cells. Ordinary native collision/advection resumes from the returned state.

use super::super::{
    MaterialFreeSurfaceDomain3, MaterialFreeSurfaceInitialCell3,
    MaterialFreeSurfaceVolumeControls3, MaterialFreeSurfaceVolumeInput3,
    MaterialFreeSurfaceVolumeProjection3, plan_free_surface_from_material,
};
use super::*;
use crate::conjugate::{FlowField, Voxel, VoxelDomain};
use crate::d3q19::Face3;
use crate::d3q19::freesurface3::transport::VolumeSlabProfile3;
use fs_exec::Cx;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_math::dd::Dd;

mod slab;

/// Only this source split may create a continuation binding. Exact native
/// moments are checked before its affine velocity is used as a predictor.
#[derive(Clone)]
pub(super) struct MappedSlabState {
    anchored_face: Face3,
    specific_enthalpy_j_kg: f64,
    affine_rate_s_inv: f64,
    profile: VolumeSlabProfile3,
}

/// Explicit discrepancy between the original caloric/source density and the
/// retained constant physical reference. This is numerical accuracy, not a
/// material uncertainty or permission to continue unmodeled thermal expansion.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FreeSurfaceDensityAgreement {
    pub maximum_allowed_kg_m3: f64,
    pub maximum_defect_kg_m3: f64,
}

/// Physical accuracy limits for one uniform planar expansion/contraction.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FreeSurfaceSlabVolumeControls {
    pub volume: MaterialFreeSurfaceVolumeControls3,
    /// Maximum original projected-face flux vs exact affine slab flux [m3/s].
    pub maximum_face_flux_defect_m3_s: f64,
    /// Maximum global/native and noncancelling destination mass defects [kg].
    pub maximum_mass_defect_kg: f64,
    /// Maximum noncancelling momentum discrepancy [kg m/s].
    pub maximum_momentum_defect_kg_m_s: f64,
    /// Maximum actual native volume vs original caloric volume demand [m3].
    pub maximum_geometry_volume_defect_m3: f64,
    /// Maximum offset between a terminal liquid-volume centroid and its
    /// native Cartesian sample point [m], at either endpoint. The explicit
    /// cut metric represents the actual liquid; native velocities store means.
    pub maximum_native_centroid_offset_m: f64,
    /// Sum of ABSOLUTE incoming native-to-face, outgoing face-to-affine and
    /// affine-to-native kinetic changes [J]. Unresolved loss and gain consume
    /// this physical accuracy allowance; neither discrepancy becomes heat.
    pub maximum_kinetic_mapping_change_j: f64,
}

/// Independent physical inventories and pressure/remap discrepancies.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FreeSurfaceSlabVolumeReport {
    pub anchored_face: Face3,
    /// Number of initially occupied layers, including a partial terminal one.
    pub initial_depth_cells: usize,
    pub initial_height_m: f64,
    pub initial_terminal_fill_fraction: f64,
    pub final_height_m: f64,
    pub newly_wet_cells: usize,
    pub vacated_cells: usize,
    pub mass_before_kg: f64,
    pub mass_after_kg: f64,
    pub mass_defect_kg: f64,
    pub enthalpy_before_j: f64,
    pub enthalpy_after_j: f64,
    /// NET extensive enthalpy; pressure-work convention is unchanged.
    pub requested_enthalpy_increment_j: f64,
    /// Stable extensive increment defect including the explicit mass-datum term.
    pub enthalpy_defect_j: f64,
    /// Actual remap mass defect times original h0 [J]. Retained separately
    /// because an enthalpy datum must not hide a real inventory discrepancy.
    pub mass_defect_enthalpy_datum_j: f64,
    pub pressure_impulse_kg_m_s: [f64; 3],
    pub momentum_before_kg_m_s: [f64; 3],
    pub momentum_after_kg_m_s: [f64; 3],
    pub momentum_defect_kg_m_s: f64,
    pub native_kinetic_energy_before_j: f64,
    pub affine_kinetic_energy_j: f64,
    pub native_kinetic_energy_j: f64,
    /// Pressure/source mechanical work required by this volume-changing step,
    /// separate from the supplied NET caloric increment. The caller owns its
    /// physical source; fixed-pressure Cp already includes ambient p*dV.
    pub required_mechanical_work_j: f64,
    /// Incoming face kinetic storage minus actual incoming native storage.
    /// Reconstructing a retained affine field is explicitly bounded; missing
    /// subcell motion is never silently restored as physical energy.
    pub native_to_face_kinetic_change_j: f64,
    /// Affine material kinetic storage minus pressure-face kinetic storage.
    /// This includes the discrete face quadrature difference, not physical heat.
    pub face_to_affine_kinetic_change_j: f64,
    /// Native cell-mean kinetic storage minus exact carried affine storage.
    /// Negative values are unresolved subcell motion, not physical heat.
    pub affine_to_native_kinetic_change_j: f64,
    /// Original pressure work minus the native kinetic INCREMENT and pressure
    /// projection loss, plus the signed incoming/outgoing mapping changes [J].
    pub mechanical_energy_defect_j: f64,
    pub maximum_face_flux_defect_m3_s: f64,
    pub geometry_volume_defect_m3: f64,
    pub maximum_native_centroid_offset_m: f64,
    pub density_agreement: FreeSurfaceDensityAgreement,
}

/// The moved native interface, continuing original enthalpy owner and its
/// pressure/volume account. The input trajectory remains unchanged.
#[derive(Clone)]
pub struct FreeSurfaceSlabVolumeStep<'s, 'm> {
    pub liquid: MaterialFreeSurfaceEnthalpy3<'s, 'm>,
    pub projection: MaterialFreeSurfaceVolumeProjection3,
    pub report: FreeSurfaceSlabVolumeReport,
}

impl<'s, 'm> MaterialFreeSurfaceEnthalpy3<'s, 'm> {
    /// Current accepted source/chart/reference discrepancy for the explicit
    /// accuracy binding. All ordinary constructors retain None and their exact
    /// original density admission. A continued heat step must remain inside
    /// this physical allowance; it cannot silently expand at constant volume.
    pub fn density_agreement(&self) -> Result<Option<FreeSurfaceDensityAgreement>, Error> {
        let Some(maximum_allowed_kg_m3) = self.density_accuracy_kg_m3 else {
            return Ok(None);
        };
        let reference = self.flow.scaling.material().density_kg_m3();
        let mut maximum_defect_kg_m3 = 0.0_f64;
        for (&location, state) in &self.states {
            let h = *self
                .h
                .get(&location)
                .ok_or(Error::Invalid("density readback lost its caloric state"))?;
            let chart = self
                .caloric
                .state_at_specific_enthalpy(h)?
                .bulk_density_kg_m3();
            maximum_defect_kg_m3 =
                maximum_defect_kg_m3.max(density_defect(reference, chart, state.density_kg_m3())?);
        }
        Ok(Some(FreeSurfaceDensityAgreement {
            maximum_allowed_kg_m3,
            maximum_defect_kg_m3,
        }))
    }

    /// Advance one uniform planar thermal source step on the native liquid.
    ///
    /// The liquid fills the entire transverse cross section from the selected
    /// stationary wall, with one optional partial terminal layer touching an
    /// unoccupied native gas layer at both endpoints.
    /// An initial state must be at canonical rest. A later moving state is
    /// admitted only with this method's retained exact per-layer native
    /// population/mass binding, unchanged enthalpy and anchored orientation.
    /// Ordinary native ticks invalidate that binding. Holes, arbitrary moving
    /// liquid, internal solids, capillary/force and coupled thermal owners refuse.
    ///
    /// The pressure solve uses the actual partial-cell width and gas face,
    /// corrects the admitted incoming affine predictor, then carries its new
    /// velocity on the exact material map x -> (rho_old/rho_new)*x. Actual
    /// incoming and outgoing native mass, momentum and kinetic storage are
    /// audited. Gas supplies no mass or enthalpy during either direction of
    /// interface motion. Terminal velocity is the liquid-volume mean, with its
    /// offset from the native sample point explicitly reported and bounded.
    ///
    /// Repeated calls compose this THERMAL SOURCE SPLIT. They do not perform
    /// native viscous collision, wall relaxation or mass-route advection, and
    /// therefore do not replace a general time integrator. Known affine
    /// molecular stress is reseeded for the new source state; unrelated stored
    /// nonequilibrium cannot enter through the retained-state admission.
    ///
    /// Endpoint source and law are explicitly borrowed at the new caloric
    /// temperature. Their original card, claims and external conditions must
    /// match. Replanning their physical units changes rho*dx^3 consistently
    /// with the remapped inventory; native population density remains the
    /// pressure coordinate one, never the thermodynamic density ratio.
    ///
    /// The original enthalpy chart/datum and clock are retained. The supplied
    /// increment is NET specific enthalpy [J/kg], including the caller's
    /// pressure-work convention; neither projection loss nor unresolved affine
    /// motion heats the liquid. Only this complete accepted result advances the
    /// clock by one original dt. Refusal/cancellation changes no input state.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn advance_uniform_slab_volume<'next>(
        &self,
        anchored_face: Face3,
        net_specific_enthalpy_increment_j_kg: f64,
        endpoint_source: &'next ResolvedLiquidSurfaceState,
        endpoint_law: &'next LiquidTemperatureLaw<'next>,
        controls: FreeSurfaceSlabVolumeControls,
        gate: &CancelGate,
        cx: &Cx<'_>,
    ) -> Result<FreeSurfaceSlabVolumeStep<'next, 'next>, Error>
    where
        'm: 'next,
    {
        checkpoint(gate, cx)?;
        for value in [
            controls.maximum_face_flux_defect_m3_s,
            controls.maximum_mass_defect_kg,
            controls.maximum_momentum_defect_kg_m_s,
            controls.maximum_geometry_volume_defect_m3,
            controls.maximum_native_centroid_offset_m,
            controls.maximum_kinetic_mapping_change_j,
        ] {
            if !value.is_finite() || value < 0.0 {
                return Err(Error::Invalid(
                    "invalid physical slab remap accuracy budget",
                ));
            }
        }
        if !net_specific_enthalpy_increment_j_kg.is_finite()
            || self.phase_conduction.is_some()
            || self.buoyancy.is_some()
            || self.rheology.is_some()
            || self.wetting.is_some()
            || self.thermocapillary.is_some()
            || self.surface_temperature_forces
            || self.flow.scaling.bulk_viscosity.is_some()
        {
            return Err(Error::Invalid(
                "uniform thermal slab cannot discard a coupled source or nonfinite heat",
            ));
        }
        #[cfg(feature = "enthalpy-coupling")]
        if self.electrical.is_some() {
            return Err(Error::Invalid(
                "uniform thermal slab cannot discard its electrical owner",
            ));
        }
        let dt = self.flow.scaling.timestep_s();
        let dx = self.flow.scaling.cell_size_m();
        if controls.volume.projection.duration_s != dt {
            return Err(Error::Invalid(
                "thermal slab must retain the original native timestep",
            ));
        }
        let next_step = self
            .accepted_steps
            .checked_add(1)
            .filter(|v| *v <= 1_u64 << 53)
            .ok_or(Error::Invalid("native liquid clock exhausted"))?;
        let next_time = (next_step as f64).mul_add(dt, self.time_origin_s);
        if !next_time.is_finite() || next_time <= self.physical_time_s() {
            return Err(Error::Invalid(
                "thermal slab time must advance representably",
            ));
        }
        let axis = anchored_face as usize / 2;
        let upper = anchored_face as usize % 2 == 1;
        let direction = if upper { -1.0 } else { 1.0 };
        let size = self.flow.lattice.domain_cells();
        size[0]
            .checked_mul(size[1])
            .and_then(|n| n.checked_mul(size[2]))
            .filter(|n| *n > 0 && *n <= controls.volume.max_cells)
            .ok_or(Error::Invalid(
                "thermal slab exceeds its complete grid budget",
            ))?;
        if self
            .volume_slab
            .as_ref()
            .is_some_and(|state| state.anchored_face != anchored_face)
        {
            return Err(Error::Invalid(
                "thermal slab continuation changed its anchored face",
            ));
        }
        let redistribution = self.flow.lattice.admit_volume_slab(
            axis,
            upper,
            self.volume_slab.as_ref().map(|state| &state.profile),
            gate,
        )?;
        let incoming_rate = self
            .volume_slab
            .as_ref()
            .map_or(0.0, |state| state.affine_rate_s_inv);
        let transverse = [0, 1, 2]
            .into_iter()
            .filter(|&a| a != axis)
            .collect::<Vec<_>>();
        let columns = size[transverse[0]]
            .checked_mul(size[transverse[1]])
            .ok_or(Error::Invalid("thermal slab cross section overflows"))?;
        let masses = self.flow.lattice.mass_snapshot();
        let m0 = self.flow.scaling.mass_unit_kg();
        let mut depth = 0;
        let mut occupied = 0usize;
        let mut h0 = None;
        let mut source_mass = Dd::ZERO;
        let mut source_momentum = [Dd::ZERO; 3];
        let mut source_kinetic = Dd::ZERO;
        let mut layer_fills = BTreeMap::new();
        let mut input = Vec::new();
        for (&location, &mass) in &masses {
            checkpoint(gate, cx)?;
            if !mass.is_finite() || mass < 0.0 {
                return Err(Error::Invalid("thermal slab has invalid native mass"));
            }
            if mass == 0.0 {
                continue;
            }
            let coordinates = FreeSurface3::transport_coordinates(location)
                .ok_or(Error::Invalid("thermal slab has unlocated inventory"))?;
            let specific = self.h[&location];
            if h0.is_some_and(|h: f64| h.to_bits() != specific.to_bits()) {
                return Err(Error::Invalid(
                    "thermal slab needs uniform original specific enthalpy",
                ));
            }
            h0 = Some(specific);
            let index = if upper {
                size[axis] - 1 - coordinates[axis]
            } else {
                coordinates[axis]
            };
            depth = depth.max(index + 1);
            occupied += 1;
            let sample = self
                .flow
                .sample_at_si(coordinates[0], coordinates[1], coordinates[2])?
                .ok_or(Error::Invalid("thermal slab lost occupied native sample"))?;
            let fill = sample.fill_fraction;
            if !fill.is_finite()
                || fill <= 0.0
                || fill > 1.0
                || layer_fills.get(&index).is_some_and(|f| *f != fill)
            {
                return Err(Error::Invalid(
                    "thermal slab requires one common planar fill per layer",
                ));
            }
            layer_fills.insert(index, fill);
            let physical_mass = Dd::from_f64(mass) * Dd::from_f64(m0);
            source_mass = source_mass + physical_mass;
            for a in 0..3 {
                let p = physical_mass * Dd::from_f64(sample.velocity_m_s[a]);
                source_momentum[a] = source_momentum[a] + p;
                source_kinetic =
                    source_kinetic + Dd::from_f64(0.5) * p * Dd::from_f64(sample.velocity_m_s[a]);
            }
            input.push(MaterialFreeSurfaceVolumeInput3 {
                coordinates,
                state: self.caloric.state_at_specific_enthalpy(specific)?,
                net_enthalpy_increment_j: value(
                    physical_mass * Dd::from_f64(net_specific_enthalpy_increment_j_kg),
                ),
            });
        }
        if depth == 0 || depth >= size[axis] || depth.checked_mul(columns) != Some(occupied) {
            return Err(Error::Invalid(
                "thermal slab must be one complete planar section with gas above it",
            ));
        }
        let terminal_fill = layer_fills[&(depth - 1)];
        if layer_fills
            .iter()
            .any(|(&index, &fill)| index + 1 < depth && fill != 1.0)
        {
            return Err(Error::Invalid(
                "only the terminal thermal slab layer may be partial",
            ));
        }
        let original_depth = (depth - 1) as f64 + terminal_fill;
        let initial_height = original_depth * dx;
        let h0 = h0.ok_or(Error::Invalid("thermal slab has no liquid"))?;
        if self
            .volume_slab
            .as_ref()
            .is_some_and(|state| state.specific_enthalpy_j_kg.to_bits() != h0.to_bits())
        {
            return Err(Error::Invalid(
                "thermal slab continuation changed its original enthalpy",
            ));
        }
        let initial = self.caloric.state_at_specific_enthalpy(h0)?;
        let proposed = self
            .caloric
            .advance_specific_energy(initial, net_specific_enthalpy_increment_j_kg)?;
        let h1 = proposed.specific_enthalpy_j_kg();
        let expected_source = self
            .law
            .at_temperature(Temperature::new(proposed.temperature_k()))?;
        let original = expected_source.material();
        let supplied = endpoint_source.material();
        if original.card_identity() != supplied.card_identity()
            || original.query_point() != supplied.query_point()
            || original.axis_quantities() != supplied.axis_quantities()
            || endpoint_source.viscosity_basis() != self.law.viscosity_basis()
            || endpoint_law.temperature_axis() != self.law.temperature_axis()
            || endpoint_source.surface_tension_n_m().is_some()
            || endpoint_source.rheology_axes().is_some()
        {
            return Err(Error::Invalid(
                "thermal slab endpoint must retain the original source at its new temperature",
            ));
        }
        let viscosity = match self.law.viscosity_basis() {
            fs_material::liquid_temperature::LiquidViscosityBasis::Dynamic => "dynamic_viscosity",
            fs_material::liquid_temperature::LiquidViscosityBasis::Kinematic => {
                "kinematic_viscosity"
            }
        };
        for name in [
            "density",
            "specific_heat_capacity",
            "thermal_conductivity",
            viscosity,
        ] {
            if original.property(name).map(|p| p.answer().receipt.selected)
                != supplied.property(name).map(|p| p.answer().receipt.selected)
                || original.property(name).map(|p| p.answer().receipt.selected)
                    != endpoint_law
                        .reference_state()
                        .material()
                        .property(name)
                        .map(|p| p.answer().receipt.selected)
            {
                return Err(Error::Invalid(
                    "thermal slab endpoint changed an original source claim",
                ));
            }
        }
        let mut scaling = plan_free_surface_from_material(
            endpoint_source,
            dx,
            dt,
            self.flow.scaling.flow.u_lattice * self.flow.scaling.velocity_scale,
        )?;
        if let Some(reference) = self.flow.scaling.reference_pressure_pa {
            scaling = scaling.with_pressure_reference(reference)?;
        }
        let m1 = scaling.mass_unit_kg();
        let column = slab::remap(original_depth, size[axis], dx, dt, m0, m1, || {
            gate.is_requested() || cx.is_cancel_requested()
        });
        checkpoint(gate, cx)?;
        let column = column.ok_or(Error::Invalid(
            "thermal slab overlap exceeds its domain or arithmetic",
        ))?;
        // A gas layer remains: the next ordinary native step must still have
        // its actual pressure boundary, not an accidentally sealed domain.
        if column.cells.len() >= size[axis] {
            return Err(Error::Invalid(
                "thermal slab requires an unoccupied gas-side native layer",
            ));
        }
        let domain = planar_domain(size, axis, upper, depth, terminal_fill, dx, gate, cx)?;
        let predictor = affine_predictor(&domain, axis, upper, incoming_rate, gate, cx)?;
        let projection = self.flow.project_enthalpy_volume_demand_on_domain(
            self.law,
            self.caloric,
            &input,
            |_, _| 0.0,
            controls.volume,
            gate,
            domain,
            Some(predictor.clone()),
        )?;
        let rate = ((m0 - m1) / m1) / dt;
        let mut flux_defect = 0.0_f64;
        let mut pressure_impulse = [Dd::ZERO; 3];
        let mut incoming_momentum_defect = 0.0;
        let mut centroid_offset = 0.5 * (1.0 - terminal_fill) * dx;
        let extent = projection.domain.face_coord(axis, size[axis]);
        for c in 0..projection.domain.cell_count() {
            checkpoint(gate, cx)?;
            if !projection.domain.is_fluid(c) {
                continue;
            }
            let position = projection.domain.coords(c);
            let center = projection
                .domain
                .center(position[0], position[1], position[2]);
            let sample = self
                .flow
                .sample_at_si(position[0], position[1], position[2])?
                .ok_or(Error::Invalid("thermal slab lost incoming velocity"))?;
            let mass = Dd::from_f64(projection.cells[&position].mass_kg);
            for f in 0..6 {
                let a = f / 2;
                let expected = if a == axis {
                    let coordinate = projection.domain.face_coord(a, position[a] + f % 2);
                    let distance = if upper {
                        extent - coordinate
                    } else {
                        coordinate
                    };
                    direction
                        * rate
                        * distance
                        * if f % 2 == 0 { -1.0 } else { 1.0 }
                        * projection.domain.face_area(c, a)
                } else {
                    0.0
                };
                let actual = projection
                    .projection
                    .field
                    .outward(&projection.domain, c, f);
                flux_defect = flux_defect.max((actual - expected).abs());
            }
            for a in 0..3 {
                centroid_offset =
                    centroid_offset.max((center[a] - (position[a] as f64 + 0.5) * dx).abs());
                let area = projection.domain.face_area(c, a);
                let velocity = 0.5
                    * (projection
                        .projection
                        .field
                        .outward(&projection.domain, c, 2 * a + 1)
                        - projection
                            .projection
                            .field
                            .outward(&projection.domain, c, 2 * a))
                    / area;
                let previous = 0.5
                    * (predictor.outward(&projection.domain, c, 2 * a + 1)
                        - predictor.outward(&projection.domain, c, 2 * a))
                    / area;
                incoming_momentum_defect +=
                    value(mass * (Dd::from_f64(previous) - Dd::from_f64(sample.velocity_m_s[a])))
                        .abs();
                pressure_impulse[a] =
                    pressure_impulse[a] + mass * (Dd::from_f64(velocity) - Dd::from_f64(previous));
            }
        }
        for target in &column.cells {
            centroid_offset = centroid_offset.max(0.5 * (1.0 - target.fill_fraction) * dx);
        }
        if !rate.is_finite()
            || !flux_defect.is_finite()
            || flux_defect > controls.maximum_face_flux_defect_m3_s
        {
            return Err(Error::Invalid(
                "projected thermal volume flux is incompatible with the planar material map",
            ));
        }
        let mut cells = Vec::new();
        let target_count = column
            .cells
            .len()
            .checked_mul(columns)
            .filter(|n| *n <= controls.volume.max_cells && *n <= self.controls.max_cells)
            .ok_or(Error::Invalid(
                "thermal slab endpoint exceeds its cell budget",
            ))?;
        cells
            .try_reserve_exact(target_count)
            .map_err(|_| Error::Invalid("thermal slab endpoint allocation refused"))?;
        for a in 0..size[transverse[0]] {
            for b in 0..size[transverse[1]] {
                checkpoint(gate, cx)?;
                for target in &column.cells {
                    let mut coordinates = [0; 3];
                    coordinates[transverse[0]] = a;
                    coordinates[transverse[1]] = b;
                    coordinates[axis] = if upper {
                        size[axis] - 1 - target.index
                    } else {
                        target.index
                    };
                    let mut velocity_m_s = [0.0; 3];
                    velocity_m_s[axis] = direction * target.momentum_kg_m_s / target.mass_kg;
                    cells.push(MaterialFreeSurfaceInitialCell3 {
                        coordinates,
                        fill_fraction: target.fill_fraction,
                        velocity_m_s,
                    });
                }
            }
        }
        let mut flow = scaling
            .free_surface3_from_cells(
                MaterialFreeSurfaceDomain3 {
                    cells: size,
                    gravity_m_s2: [0.0; 3],
                    contact: crate::d3q19::freesurface3::ContactModel3::Neutral,
                    solid_tiles: &[],
                    max_cells: controls.volume.max_cells,
                },
                &cells,
                cx,
            )?
            .with_conversion_redistribution(redistribution);
        // The carried physical velocity is affine in the NEW material
        // coordinate. Seed its actual molecular stress through the existing
        // native strain owner, not equilibrium populations at that velocity.
        let mut strain = [0.0; 6];
        strain[axis] = (m0 - m1) / m0;
        for target in &cells {
            checkpoint(gate, cx)?;
            flow.lattice.initialize_strain_at_with_bulk(
                target.coordinates,
                strain,
                flow.scaling.flow.tau,
                None,
            )?;
        }
        let mut liquid = flow.with_enthalpy_model_seeded_density(
            endpoint_law,
            self.caloric,
            self.controls,
            |_| h1,
            None,
            None,
            None,
            None,
            Some(controls.volume.maximum_source_density_defect_kg_m3),
            gate,
        )?;
        liquid.time_origin_s = self.time_origin_s;
        liquid.accepted_steps = next_step;
        let mut native_momentum = [Dd::ZERO; 3];
        let mut native_mass = Dd::ZERO;
        let mut native_kinetic = Dd::ZERO;
        let mut mass_mapping_defect = 0.0;
        let mut momentum_mapping_defect = 0.0;
        let mut actual_volume = Dd::ZERO;
        for target in &cells {
            checkpoint(gate, cx)?;
            let p = target.coordinates;
            let index = if upper {
                size[axis] - 1 - p[axis]
            } else {
                p[axis]
            };
            let expected = &column.cells[index];
            let location = FreeSurface3::transport_location(p)
                .ok_or(Error::Invalid("thermal slab lost target address"))?;
            let mass =
                Dd::from_f64(liquid.flow.lattice.transport_cell_mass(location)?) * Dd::from_f64(m1);
            native_mass = native_mass + mass;
            let sample = liquid
                .flow
                .sample_at_si(p[0], p[1], p[2])?
                .ok_or(Error::Invalid("thermal slab lost target liquid"))?;
            mass_mapping_defect += value(mass - Dd::from_f64(expected.mass_kg)).abs();
            actual_volume =
                actual_volume + Dd::from_f64(sample.fill_fraction) * Dd::from_f64(dx * dx * dx);
            for a in 0..3 {
                let momentum = mass * Dd::from_f64(sample.velocity_m_s[a]);
                native_momentum[a] = native_momentum[a] + momentum;
                let expected_p = if a == axis {
                    direction * expected.momentum_kg_m_s
                } else {
                    0.0
                };
                momentum_mapping_defect += (value(momentum) - expected_p).abs();
                native_kinetic = native_kinetic
                    + Dd::from_f64(0.5) * momentum * Dd::from_f64(sample.velocity_m_s[a]);
            }
        }
        let mass_before = self.flow.ledger_mass_kg()?;
        let mass_after = liquid.flow.ledger_mass_kg()?;
        let mass_defect = (mass_after - mass_before)
            .abs()
            .max(value(native_mass - source_mass).abs())
            .max(mass_mapping_defect);
        let pressure_impulse = pressure_impulse.map(value);
        let source_momentum = source_momentum.map(value);
        let native_momentum = native_momentum.map(value);
        let analytic_p = direction * column.source_momentum_kg_m_s * columns as f64;
        let momentum_defect = (0..3)
            .map(|a| {
                let expected = if a == axis { analytic_p } else { 0.0 };
                (source_momentum[a] + pressure_impulse[a] - expected).abs()
                    + (native_momentum[a] - expected).abs()
            })
            .sum::<f64>()
            .max(momentum_mapping_defect + incoming_momentum_defect);
        let demanded_volume = projection.cells.values().fold(Dd::ZERO, |sum, cell| {
            sum + Dd::from_f64(cell.proposed_volume_m3)
        });
        let volume_defect = value(actual_volume - demanded_volume).abs();
        let affine_kinetic = column.source_kinetic_energy_j * columns as f64;
        let pr = &projection.projection.report;
        // Retain the pressure owner's compensated kinetic increment. Rounded
        // absolute endpoint diagnostics are not subtracted to recover work.
        let face_before = Dd::from_f64(pr.kinetic_energy_before_j);
        let face_after = face_before + Dd::from_f64(pr.kinetic_energy_change_j);
        let native_to_face_dd = face_before - source_kinetic;
        let face_to_affine_dd = Dd::from_f64(affine_kinetic) - face_after;
        let affine_to_native_dd = native_kinetic - Dd::from_f64(affine_kinetic);
        let native_to_face = value(native_to_face_dd);
        let face_to_affine = value(face_to_affine_dd);
        let affine_to_native = value(affine_to_native_dd);
        let mapping_change = native_to_face.abs() + face_to_affine.abs() + affine_to_native.abs();
        let enthalpy_before = self.total_enthalpy_j()?;
        let enthalpy_after = liquid.total_enthalpy_j()?;
        // Absolute H values are diagnostics. Form the original-chart
        // increment before adding its datum, retaining product and summation
        // roundoff in physical mass. A real mass defect still carries h0 and
        // is exposed separately; no convenient reference erases it.
        let requested_dd = source_mass * Dd::from_f64(net_specific_enthalpy_increment_j_kg);
        let mass_datum = (native_mass - source_mass) * Dd::from_f64(h0);
        let enthalpy_defect =
            value(native_mass * (Dd::from_f64(h1) - Dd::from_f64(h0)) + mass_datum - requested_dd);
        let requested = value(requested_dd);
        let mass_datum = value(mass_datum);
        let required_work_dd = Dd::from_f64(pr.pressure_source_work_j)
            - Dd::from_f64(pr.fixed_boundary_pressure_work_j)
            - Dd::from_f64(pr.open_boundary_pressure_work_j);
        let mechanical_defect = value(
            required_work_dd
                - (native_kinetic - source_kinetic)
                - Dd::from_f64(pr.projection_loss_j)
                - Dd::from_f64(pr.resistance_correction_work_j)
                + native_to_face_dd
                + face_to_affine_dd
                + affine_to_native_dd,
        );
        let required_work = value(required_work_dd);
        let source_kinetic = value(source_kinetic);
        let native_kinetic = value(native_kinetic);
        for measured in [
            mass_defect,
            momentum_defect,
            volume_defect,
            centroid_offset,
            source_kinetic,
            affine_kinetic,
            native_kinetic,
            mapping_change,
            mass_datum,
            enthalpy_defect,
            mechanical_defect,
        ] {
            if !measured.is_finite() {
                return Err(Error::Invalid(
                    "thermal slab physical inventory is unrepresentable",
                ));
            }
        }
        if mass_defect > controls.maximum_mass_defect_kg
            || momentum_defect > controls.maximum_momentum_defect_kg_m_s
            || volume_defect > controls.maximum_geometry_volume_defect_m3
            || centroid_offset > controls.maximum_native_centroid_offset_m
            || mapping_change > controls.maximum_kinetic_mapping_change_j
            || enthalpy_defect.abs() > controls.volume.maximum_enthalpy_defect_j
            || mechanical_defect.abs() > controls.volume.projection.maximum_energy_defect_j
        {
            return Err(Error::Invalid(
                "thermal slab native mass/enthalpy/momentum/work readback exceeds its budget",
            ));
        }
        let density_agreement = liquid
            .density_agreement()?
            .ok_or(Error::Invalid("thermal slab lost density accuracy binding"))?;
        liquid.volume_slab = Some(MappedSlabState {
            anchored_face,
            specific_enthalpy_j_kg: h1,
            affine_rate_s_inv: ((m0 - m1) / m0) / dt,
            profile: liquid.flow.lattice.volume_slab_profile(axis, upper, gate)?,
        });
        checkpoint(gate, cx)?;
        Ok(FreeSurfaceSlabVolumeStep {
            report: FreeSurfaceSlabVolumeReport {
                anchored_face,
                initial_depth_cells: depth,
                initial_height_m: initial_height,
                initial_terminal_fill_fraction: terminal_fill,
                final_height_m: column.height_m,
                newly_wet_cells: column.cells.len().saturating_sub(depth) * columns,
                vacated_cells: depth.saturating_sub(column.cells.len()) * columns,
                mass_before_kg: mass_before,
                mass_after_kg: mass_after,
                mass_defect_kg: mass_defect,
                enthalpy_before_j: enthalpy_before,
                enthalpy_after_j: enthalpy_after,
                requested_enthalpy_increment_j: requested,
                enthalpy_defect_j: enthalpy_defect,
                mass_defect_enthalpy_datum_j: mass_datum,
                pressure_impulse_kg_m_s: pressure_impulse,
                momentum_before_kg_m_s: source_momentum,
                momentum_after_kg_m_s: native_momentum,
                momentum_defect_kg_m_s: momentum_defect,
                native_kinetic_energy_before_j: source_kinetic,
                affine_kinetic_energy_j: affine_kinetic,
                native_kinetic_energy_j: native_kinetic,
                required_mechanical_work_j: required_work,
                native_to_face_kinetic_change_j: native_to_face,
                face_to_affine_kinetic_change_j: face_to_affine,
                affine_to_native_kinetic_change_j: affine_to_native,
                mechanical_energy_defect_j: mechanical_defect,
                maximum_face_flux_defect_m3_s: flux_defect,
                geometry_volume_defect_m3: volume_defect,
                maximum_native_centroid_offset_m: centroid_offset,
                density_agreement,
            },
            liquid,
            projection,
        })
    }
}

#[allow(clippy::too_many_arguments)]
fn planar_domain(
    size: [usize; 3],
    axis: usize,
    upper: bool,
    depth: usize,
    terminal_fill: f64,
    dx: f64,
    gate: &CancelGate,
    cx: &Cx<'_>,
) -> Result<VoxelDomain, Error> {
    let mut widths = size.map(|n| vec![dx; n]);
    let terminal = if upper { size[axis] - depth } else { depth - 1 };
    let gas = if upper { terminal - 1 } else { terminal + 1 };
    widths[axis][terminal] = terminal_fill * dx;
    widths[axis][gas] = (2.0 - terminal_fill) * dx;
    let mut domain = VoxelDomain::graded(widths)?;
    for c in 0..domain.cell_count() {
        checkpoint(gate, cx)?;
        let p = domain.coords(c);
        let index = if upper {
            size[axis] - 1 - p[axis]
        } else {
            p[axis]
        };
        if index >= depth {
            domain.set(p[0], p[1], p[2], Voxel::Solid(0));
        }
    }
    Ok(domain)
}

/// This known affine predictor includes the actual internal gas-face velocity.
/// Adjacent-center copying at a partial opening would lose the retained motion.
fn affine_predictor(
    domain: &VoxelDomain,
    axis: usize,
    upper: bool,
    rate: f64,
    gate: &CancelGate,
    cx: &Cx<'_>,
) -> Result<FlowField, Error> {
    let field = FlowField::quiescent(domain);
    let mut arrays = field.face_arrays().map(<[f64]>::to_vec);
    let extent = domain.face_coord(axis, domain.dims()[axis]);
    let direction = if upper { -1.0 } else { 1.0 };
    for c in 0..domain.cell_count() {
        checkpoint(gate, cx)?;
        if !domain.is_fluid(c) {
            continue;
        }
        let p = domain.coords(c);
        for f in [2 * axis, 2 * axis + 1] {
            let coordinate = domain.face_coord(axis, p[axis] + f % 2);
            let distance = if upper {
                extent - coordinate
            } else {
                coordinate
            };
            let flux = direction * rate * distance * domain.face_area(c, axis);
            if !flux.is_finite() {
                return Err(Error::Invalid(
                    "retained thermal slab face motion overflows",
                ));
            }
            let (slot, _) = field.slot(domain, c, f);
            arrays[axis][slot] = flux;
        }
    }
    let [x, y, z] = arrays;
    Ok(FlowField::from_face_arrays(domain, x, y, z))
}

fn checkpoint(gate: &CancelGate, cx: &Cx<'_>) -> Result<(), Error> {
    poll(gate)?;
    cx.checkpoint()
        .map_err(|_| FreeSurfaceError3::Cancelled.into())
}

fn value(number: Dd) -> f64 {
    number.hi + number.lo
}

fn density_defect(reference: f64, chart: f64, source: f64) -> Result<f64, Error> {
    let defect = (chart - reference)
        .abs()
        .max((source - reference).abs())
        .max((chart - source).abs());
    if !defect.is_finite() {
        Err(Error::Invalid(
            "source/chart/reference density difference overflows",
        ))
    } else {
        Ok(defect)
    }
}

/// Retain the original source sample and fixed reference storage. Only the
/// explicit physical arithmetic allowance relaxes equality; source claims and
/// chart interpolation are unchanged, and any larger thermal density change
/// still refuses the ordinary constant-volume native step.
pub(super) fn resolve_density_agreement(
    law: &LiquidTemperatureLaw<'_>,
    caloric: &EquilibriumEnthalpyPhaseCurve,
    h: f64,
    reference_density: f64,
    maximum: f64,
) -> Result<(ResolvedLiquidState, FluidTransportState, f64), Error> {
    let phase = caloric.state_at_specific_enthalpy(h)?;
    if phase.phase() != SolidLiquidPhase::Liquid {
        return Err(Error::Invalid(
            "thermal volume continuation requires pure liquid",
        ));
    }
    let state = law.at_temperature(Temperature::new(phase.temperature_k()))?;
    let defect = density_defect(
        reference_density,
        phase.bulk_density_kg_m3(),
        state.density_kg_m3(),
    )?;
    if defect > maximum {
        return Err(Error::Invalid(
            "ordinary native step exceeds its fixed-reference density accuracy budget",
        ));
    }
    let actual = state.transport_state().map_err(MaterialLbmError::from)?;
    let transport = FluidTransportState::try_new(
        reference_density,
        actual.specific_heat_j_kg_k(),
        actual.thermal_conductivity_w_m_k(),
        actual.dynamic_viscosity_pa_s(),
    )
    .map_err(MaterialLbmError::from)?;
    Ok((state, transport, defect))
}
