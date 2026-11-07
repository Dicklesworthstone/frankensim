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
use crate::d3q19::Face3;
use fs_exec::Cx;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_math::dd::Dd;

mod slab;

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
    /// Sum of ABSOLUTE face-to-affine and affine-to-native kinetic changes
    /// [J]. Both unresolved loss and gain consume this physical accuracy
    /// allowance; neither discrepancy becomes material heat.
    pub maximum_kinetic_mapping_change_j: f64,
}

/// Independent physical inventories and pressure/remap discrepancies.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FreeSurfaceSlabVolumeReport {
    pub anchored_face: Face3,
    pub initial_depth_cells: usize,
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
    pub momentum_after_kg_m_s: [f64; 3],
    pub momentum_defect_kg_m_s: f64,
    pub affine_kinetic_energy_j: f64,
    pub native_kinetic_energy_j: f64,
    /// Pressure/source mechanical work required by this volume-changing step,
    /// separate from the supplied NET caloric increment. The caller owns its
    /// physical source; fixed-pressure Cp already includes ambient p*dV.
    pub required_mechanical_work_j: f64,
    /// Affine material kinetic storage minus pressure-face kinetic storage.
    /// This includes the discrete face quadrature difference, not physical heat.
    pub face_to_affine_kinetic_change_j: f64,
    /// Native cell-mean kinetic storage minus exact carried affine storage.
    /// Negative values are unresolved subcell motion, not physical heat.
    pub affine_to_native_kinetic_change_j: f64,
    /// Original pressure work minus native kinetic storage, pressure projection
    /// loss and the separately reported mapping losses [J].
    pub mechanical_energy_defect_j: f64,
    pub maximum_face_flux_defect_m3_s: f64,
    pub geometry_volume_defect_m3: f64,
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

    /// Advance one uniformly heated/cooled planar liquid slab from exact rest.
    ///
    /// The initial positive inventory must fill one complete Cartesian slab
    /// anchored at the selected domain face. Its opposite face touches gas;
    /// other domain faces are stationary walls. Partial initial cells, holes,
    /// disconnected bodies, internal solids, numerical pressure, stored
    /// non-equilibrium stress, force/capillary and coupled thermal owners refuse.
    ///
    /// This first-order split applies the original physical pressure impulse,
    /// then carries its affine normal velocity on the exact material map
    /// x -> (rho_old/rho_new)*x. The pressure flux and native momentum/kinetic
    /// readback must agree within explicit physical budgets. The new surface
    /// actually expands into gas or retreats within old liquid; gas supplies no
    /// mass or enthalpy. This is not a general cut-cell or repeated thermal VOF
    /// integrator. Continued ordinary steps use the existing native mass routes.
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
        let redistribution = self.flow.lattice.admit_quiescent_volume_slab(gate)?;
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
        let mut input = Vec::new();
        for (&location, &mass) in &masses {
            checkpoint(gate, cx)?;
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
            source_mass = source_mass + Dd::from_f64(mass) * Dd::from_f64(m0);
            input.push(MaterialFreeSurfaceVolumeInput3 {
                coordinates,
                state: self.caloric.state_at_specific_enthalpy(specific)?,
                net_enthalpy_increment_j: m0 * net_specific_enthalpy_increment_j_kg,
            });
        }
        if depth == 0 || depth >= size[axis] || depth.checked_mul(columns) != Some(occupied) {
            return Err(Error::Invalid(
                "thermal slab must be one complete full-cell section with gas above it",
            ));
        }
        let h0 = h0.ok_or(Error::Invalid("thermal slab has no liquid"))?;
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
        let column = slab::remap(depth, size[axis], dx, dt, m0, m1, || {
            gate.is_requested() || cx.is_cancel_requested()
        });
        checkpoint(gate, cx)?;
        let column = column.ok_or(Error::Invalid(
            "thermal slab overlap exceeds its domain or arithmetic",
        ))?;
        // A gas layer remains: the next ordinary native step must still have
        // its actual pressure boundary, not an accidentally sealed domain.
        if column.height_m >= size[axis] as f64 * dx {
            return Err(Error::Invalid(
                "thermal slab reaches the opposite closed wall",
            ));
        }
        let projection = self.flow.project_enthalpy_volume_demand(
            self.law,
            self.caloric,
            &input,
            |_, _| 0.0,
            controls.volume,
            gate,
        )?;
        let rate = ((m0 - m1) / m1) / dt;
        let mut flux_defect = 0.0_f64;
        let mut pressure_impulse = [Dd::ZERO; 3];
        for c in 0..projection.domain.cell_count() {
            checkpoint(gate, cx)?;
            if !projection.domain.is_fluid(c) {
                continue;
            }
            let position = projection.domain.coords(c);
            for f in 0..6 {
                let a = f / 2;
                let expected = if a == axis {
                    let face_index = position[a] + f % 2;
                    let distance = if upper {
                        size[axis] - face_index
                    } else {
                        face_index
                    };
                    direction
                        * rate
                        * (distance as f64 * dx)
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
                let velocity = 0.5
                    * (projection
                        .projection
                        .field
                        .outward(&projection.domain, c, 2 * a + 1)
                        - projection
                            .projection
                            .field
                            .outward(&projection.domain, c, 2 * a))
                    / projection.domain.face_area(c, a);
                pressure_impulse[a] =
                    pressure_impulse[a] + Dd::from_f64(m0) * Dd::from_f64(velocity);
            }
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
            let mass = Dd::from_f64(liquid.flow.lattice.transport_cell_mass(location)?)
                * Dd::from_f64(m1);
            native_mass = native_mass + mass;
            let sample = liquid
                .flow
                .sample_at_si(p[0], p[1], p[2])?
                .ok_or(Error::Invalid("thermal slab lost target liquid"))?;
            mass_mapping_defect += value(mass - Dd::from_f64(expected.mass_kg)).abs();
            actual_volume = actual_volume
                + Dd::from_f64(sample.fill_fraction) * Dd::from_f64(projection.domain.volume(0));
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
        let mass_defect = (mass_after - mass_before).abs()
            .max(value(native_mass - source_mass).abs()).max(mass_mapping_defect);
        let pressure_impulse = pressure_impulse.map(value);
        let native_momentum = native_momentum.map(value);
        let analytic_p = direction * column.source_momentum_kg_m_s * columns as f64;
        let momentum_defect = (0..3)
            .map(|a| {
                let expected = if a == axis { analytic_p } else { 0.0 };
                (pressure_impulse[a] - expected).abs() + (native_momentum[a] - expected).abs()
            })
            .sum::<f64>()
            .max(momentum_mapping_defect);
        let demanded_volume = projection.cells.values().fold(Dd::ZERO, |sum, cell| {
            sum + Dd::from_f64(cell.proposed_volume_m3)
        });
        let volume_defect = value(actual_volume - demanded_volume).abs();
        let affine_kinetic = column.source_kinetic_energy_j * columns as f64;
        let native_kinetic = value(native_kinetic);
        let face_kinetic = projection.projection.report.kinetic_energy_after_j;
        let face_to_affine = affine_kinetic - face_kinetic;
        let affine_to_native = native_kinetic - affine_kinetic;
        let mapping_change = face_to_affine.abs() + affine_to_native.abs();
        let enthalpy_before = self.total_enthalpy_j()?;
        let enthalpy_after = liquid.total_enthalpy_j()?;
        // Absolute H values are diagnostics. Form the original-chart
        // increment before adding its datum, retaining product and summation
        // roundoff in physical mass. A real mass defect still carries h0 and
        // is exposed separately; no convenient reference erases it.
        let requested_dd = source_mass * Dd::from_f64(net_specific_enthalpy_increment_j_kg);
        let mass_datum = (native_mass - source_mass) * Dd::from_f64(h0);
        let enthalpy_defect = value(
            native_mass * (Dd::from_f64(h1) - Dd::from_f64(h0))
                + mass_datum - requested_dd,
        );
        let requested = value(requested_dd);
        let mass_datum = value(mass_datum);
        let pr = &projection.projection.report;
        let required_work = pr.pressure_source_work_j
            - pr.fixed_boundary_pressure_work_j
            - pr.open_boundary_pressure_work_j;
        let mechanical_defect = required_work - native_kinetic - pr.projection_loss_j
            + face_to_affine
            + affine_to_native;
        for measured in [
            mass_defect,
            momentum_defect,
            volume_defect,
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
        checkpoint(gate, cx)?;
        Ok(FreeSurfaceSlabVolumeStep {
            report: FreeSurfaceSlabVolumeReport {
                anchored_face,
                initial_depth_cells: depth,
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
                momentum_after_kg_m_s: native_momentum,
                momentum_defect_kg_m_s: momentum_defect,
                affine_kinetic_energy_j: affine_kinetic,
                native_kinetic_energy_j: native_kinetic,
                required_mechanical_work_j: required_work,
                face_to_affine_kinetic_change_j: face_to_affine,
                affine_to_native_kinetic_change_j: affine_to_native,
                mechanical_energy_defect_j: mechanical_defect,
                maximum_face_flux_defect_m3_s: flux_defect,
                geometry_volume_defect_m3: volume_defect,
                density_agreement,
            },
            liquid,
            projection,
        })
    }
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
