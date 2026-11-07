//! Fixed-mass liquid calorics drive the existing physical pressure projection.
//! This returns a volume-compatible flux proposal, not a VOF evolution step.

use super::*;
use std::collections::BTreeMap;
use std::sync::Arc;

use super::enthalpy::FreeSurfaceEnthalpyError as Error;
use crate::conjugate::{
    DensityProjection, DensityProjectionControls, FlowFace, FlowField, InternalPressureFace, Voxel,
    VoxelDomain,
};
use crate::d3q19::Face3;
use crate::d3q19::freesurface3::transport::MassLocation3;
use fs_material::liquid::ResolvedLiquidState;
use fs_material::liquid_temperature::LiquidTemperatureLaw;
use fs_material::phase::{EquilibriumEnthalpyPhaseCurve, EquilibriumPhaseState, SolidLiquidPhase};
use fs_qty::Temperature;

/// Original caloric state and NET extensive enthalpy change at one occupied
/// native cell. The existing native inventory supplies the fixed mass.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialFreeSurfaceVolumeInput3 {
    pub coordinates: [usize; 3],
    pub state: EquilibriumPhaseState,
    /// J, including the external owner's already determined heat/work split.
    /// This is not an extra latent heat or automatically an external heat input.
    pub net_enthalpy_increment_j: f64,
}

/// Frozen-geometry pressure and source-accuracy budgets.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialFreeSurfaceVolumeControls3 {
    pub projection: DensityProjectionControls,
    /// Maximum complete Cartesian cell count, including gas and solid cells.
    pub max_cells: usize,
    /// Maximum |M/rho_initial - occupied voxel volume| [m3].
    /// An admitted initial mismatch is never turned into an expansion source.
    pub maximum_initial_volume_defect_m3: f64,
    /// Maximum |rho_chart - rho_source(T)| at either endpoint [kg/m3].
    /// This explicit arithmetic/consistency allowance never replaces either
    /// original density law or describes empirical source uncertainty.
    pub maximum_source_density_defect_kg_m3: f64,
    /// Sum of |requested - represented M*delta_h| [J].
    pub maximum_enthalpy_defect_j: f64,
}

/// Fixed inventory, original source snapshots and thermal volume demand.
/// Identical temperatures share source snapshots across cells.
#[derive(Clone, Debug)]
pub struct MaterialFreeSurfaceVolumeCell3 {
    pub mass_kg: f64,
    pub initial: EquilibriumPhaseState,
    pub proposed: EquilibriumPhaseState,
    pub initial_material: Arc<ResolvedLiquidState>,
    pub proposed_material: Arc<ResolvedLiquidState>,
    pub requested_enthalpy_increment_j: f64,
    pub represented_enthalpy_increment_j: f64,
    pub initial_volume_m3: f64,
    pub proposed_volume_m3: f64,
    /// Stable M*(1/rho_proposed - 1/rho_initial) [m3].
    pub volume_change_m3: f64,
}

/// Original-chart volume demand and its native pressure/flux proposal.
#[derive(Clone, Debug)]
pub struct MaterialFreeSurfaceVolumeProjection3 {
    /// Exact binary occupied-voxel domain for the returned staggered fluxes.
    /// Inactive entries represent gas or walls, not an inferred gas material.
    pub domain: VoxelDomain,
    pub projection: DensityProjection,
    pub cells: BTreeMap<[usize; 3], MaterialFreeSurfaceVolumeCell3>,
    pub volume_source_m3_s: Vec<f64>,
    pub gas_pressure_faces: Vec<InternalPressureFace>,
    pub maximum_initial_volume_defect_m3: f64,
    /// Largest actually measured chart/source density discrepancy [kg/m3].
    pub maximum_source_density_defect_kg_m3: f64,
    pub absolute_enthalpy_defect_j: f64,
}

impl MaterialFreeSurface3<'_> {
    /// Compute expansion/contraction-compatible physical fluxes from the
    /// original liquid h/rho chart, with every native cell's mass held fixed.
    ///
    /// For each occupied cell, S=M*(v(h_new)-v(h_old))/dt. The original source
    /// temperature law must agree with BOTH chart densities within the explicit
    /// physical density budget and retain the frozen material's density claim
    /// and external conditions. The maximum discrepancy is returned; independently
    /// rounded T(h) and rho(T) evaluations are never hidden by a guessed
    /// epsilon. Heating can expand or contract according to that source; no
    /// expansion coefficient, material-name rule or pressure-density rescaling
    /// is introduced.
    ///
    /// The pressure operator uses the actual initial M/voxel-volume inertia,
    /// actual Guo-corrected velocity predictor, impermeable voxel walls and
    /// explicit pressure faces towards gas. Its source, boundary work and
    /// numerical projection loss remain in projection.report. Those mechanical
    /// correction pressures are not the fixed thermodynamic chart pressure;
    /// their required work is not silently inserted into the supplied enthalpy.
    ///
    /// This first adapter admits full occupied voxels and exactly empty
    /// interface ghosts. Partial fill, conversion carry and piston-remapped
    /// geometry require an explicit cut-cell/remap owner and refuse. In
    /// particular the historical binary classifier seeds half-filled
    /// interfaces; use explicitly supplied full-cell initial states here.
    /// The phase must remain liquid. The source/core projection does not invent
    /// a mobility law for newly solid or mixed material.
    ///
    /// `gas_pressure` returns the mechanical CORRECTION pressure [Pa] on each
    /// liquid-to-gas voxel face, in one shared gauge; zero is a pressure-release
    /// correction. Existing velocity already carries its accepted force history.
    /// No additional gravity, capillary or pressure impulse is applied here.
    ///
    /// The result is a frozen-geometry pressure proposal. Accepted lattice
    /// populations, VOF fill, mass, enthalpy, force ownership and time are
    /// unchanged, including on success. A conservative moving-interface owner
    /// must consume these fluxes with its actual material/energy transport;
    /// ordinary free-surface stepping does not automatically use this proposal.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn project_enthalpy_volume_demand(
        &self,
        law: &LiquidTemperatureLaw<'_>,
        caloric: &EquilibriumEnthalpyPhaseCurve,
        input: &[MaterialFreeSurfaceVolumeInput3],
        mut gas_pressure: impl FnMut([usize; 3], Face3) -> f64,
        controls: MaterialFreeSurfaceVolumeControls3,
        gate: &CancelGate,
    ) -> Result<MaterialFreeSurfaceVolumeProjection3, Error> {
        poll(gate)?;
        for value in [
            controls.maximum_initial_volume_defect_m3,
            controls.maximum_source_density_defect_kg_m3,
            controls.maximum_enthalpy_defect_j,
        ] {
            if !value.is_finite() || value < 0.0 {
                return Err(Error::Invalid("invalid thermal-volume accuracy budget"));
            }
        }
        let dt = positive(
            controls.projection.duration_s,
            "invalid volume projection duration",
        )?;
        if self.lattice.carry() != 0.0 || self.lattice.geometry_epoch() != 0 {
            return Err(Error::Invalid(
                "volume projection needs located mass and original stationary voxel geometry",
            ));
        }
        let source = self.scaling.material().material();
        let reference = law.reference_state().material();
        if caloric.material_card_identity() != source.card_identity()
            || reference.card_identity() != source.card_identity()
            || reference.query_point() != source.query_point()
            || reference.axis_quantities() != source.axis_quantities()
            || reference
                .property("density")
                .map(|p| p.answer().receipt.selected)
                != source
                    .property("density")
                    .map(|p| p.answer().receipt.selected)
        {
            return Err(Error::Invalid(
                "volume projection requires the original material, density claim and query context",
            ));
        }
        let context = source.query();
        for receipt in caloric.source_receipts() {
            poll(gate)?;
            if matches!(
                receipt.property.as_str(),
                "density" | "specific_heat_capacity"
            ) && reference
                .property(&receipt.property)
                .map(|property| property.answer().receipt.selected)
                != Some(receipt.selected)
            {
                return Err(Error::Invalid(
                    "caloric volume law changed an original liquid source claim",
                ));
            }
            for (axis, value) in &receipt.query_point {
                if axis != law.temperature_axis()
                    && (context.axes().get(axis) != Some(value)
                        || context.axis_quantities().get(axis) != receipt.axis_quantities.get(axis))
                {
                    return Err(Error::Invalid(
                        "caloric volume law changed an external source condition",
                    ));
                }
            }
        }
        let dims = self.lattice.grid().tile_dims();
        let [Some(nx), Some(ny), Some(nz)] = [dims.0, dims.1, dims.2].map(|n| n.checked_mul(TILE))
        else {
            return Err(Error::Invalid("thermal-volume cell dimensions overflow"));
        };
        let size = [nx, ny, nz];
        let count = size[0]
            .checked_mul(size[1])
            .and_then(|n| n.checked_mul(size[2]))
            .filter(|&n| n > 0 && n <= controls.max_cells)
            .ok_or(Error::Invalid(
                "thermal-volume grid exceeds its cell budget",
            ))?;
        if input.len() > count {
            return Err(Error::Invalid(
                "thermal-volume input exceeds the original grid",
            ));
        }
        let mut domain = VoxelDomain::new(size[0], size[1], size[2], self.scaling.cell_size_m())?;
        for c in 0..count {
            poll(gate)?;
            let [x, y, z] = domain.coords(c);
            domain.set(x, y, z, Voxel::Solid(0));
        }
        let mut masses = vec![0.0; count];
        let mut velocities = vec![[0.0; 3]; count];
        let mut densities = vec![0.0; count];
        let mut sources = vec![0.0; count];
        let native_mass = self.lattice.mass_snapshot();
        let mut occupied = 0usize;
        for (&location, &mass) in &native_mass {
            poll(gate)?;
            if !mass.is_finite() || mass < 0.0 {
                return Err(Error::Invalid(
                    "thermal-volume input has invalid native inventory",
                ));
            }
            if location == MassLocation3::Carry || mass == 0.0 {
                continue;
            }
            let coordinates = FreeSurface3::transport_coordinates(location).ok_or(
                Error::Invalid("thermal-volume mass has no Cartesian location"),
            )?;
            let [x, y, z] = coordinates;
            let sample = self.lattice.sample_at(x, y, z).ok_or(Error::Invalid(
                "thermal-volume inventory is not an occupied native cell",
            ))?;
            // Check the actual tracked inventory, not its clamped fill readback.
            if sample.fill_fraction != 1.0 || mass != sample.density {
                return Err(Error::Invalid(
                    "partial or overfilled interfaces need an explicit cut-cell volume metric",
                ));
            }
            let c = domain.index(x, y, z);
            domain.set(x, y, z, Voxel::Fluid);
            masses[c] = positive(
                mass * self.scaling.mass_unit_kg(),
                "native physical mass is unrepresentable",
            )?;
            densities[c] = positive(
                masses[c] / domain.volume(c),
                "physical voxel inertia is unrepresentable",
            )?;
            for axis in 0..3 {
                velocities[c][axis] = self.scaling.velocity_si(sample.velocity[axis])?;
            }
            occupied += 1;
        }
        if occupied == 0 || input.len() != occupied {
            return Err(Error::Invalid(
                "thermal-volume states must cover every occupied native cell exactly once",
            ));
        }
        let mut cache = BTreeMap::new();
        let mut cells = BTreeMap::new();
        let mut maximum_initial_volume_defect_m3 = 0.0_f64;
        let mut maximum_source_density_defect_kg_m3 = 0.0_f64;
        let mut absolute_enthalpy_defect_j = 0.0_f64;
        for entry in input {
            poll(gate)?;
            let [x, y, z] = entry.coordinates;
            if x >= size[0]
                || y >= size[1]
                || z >= size[2]
                || cells.contains_key(&entry.coordinates)
            {
                return Err(Error::Invalid(
                    "duplicate or out-of-domain thermal-volume cell",
                ));
            }
            let c = domain.index(x, y, z);
            if !domain.is_fluid(c) || entry.state.phase_curve_identity() != caloric.identity() {
                return Err(Error::Invalid(
                    "thermal-volume state needs occupied liquid on the original caloric chart",
                ));
            }
            let mass = masses[c];
            let requested = finite(
                entry.net_enthalpy_increment_j,
                "nonfinite enthalpy volume increment",
            )?;
            let dh = finite(
                requested / mass,
                "specific enthalpy increment is unrepresentable",
            )?;
            if requested != 0.0 && dh == 0.0 {
                return Err(Error::Invalid(
                    "specific enthalpy volume increment underflows",
                ));
            }
            let next = caloric.advance_specific_energy(entry.state, dh)?;
            let (initial_material, initial_density_defect) =
                source_state(law, entry.state, &mut cache)?;
            let (proposed_material, proposed_density_defect) = source_state(law, next, &mut cache)?;
            maximum_source_density_defect_kg_m3 = maximum_source_density_defect_kg_m3
                .max(initial_density_defect)
                .max(proposed_density_defect);
            if maximum_source_density_defect_kg_m3 > controls.maximum_source_density_defect_kg_m3 {
                return Err(Error::Invalid(
                    "caloric/source density discrepancy exceeds its physical accuracy budget",
                ));
            }
            let represented = finite(
                mass * (next.specific_enthalpy_j_kg() - entry.state.specific_enthalpy_j_kg()),
                "represented enthalpy increment is unrepresentable",
            )?;
            if requested != 0.0 && represented == 0.0 {
                return Err(Error::Invalid(
                    "represented extensive enthalpy increment underflows",
                ));
            }
            absolute_enthalpy_defect_j = finite(
                absolute_enthalpy_defect_j + (requested - represented).abs(),
                "enthalpy volume residual overflows",
            )?;
            if absolute_enthalpy_defect_j > controls.maximum_enthalpy_defect_j {
                return Err(Error::Invalid(
                    "thermal-volume enthalpy increment exceeds its accuracy budget",
                ));
            }
            let rho0 = entry.state.bulk_density_kg_m3();
            let rho1 = next.bulk_density_kg_m3();
            let v0 = positive(mass / rho0, "initial material volume is unrepresentable")?;
            let v1 = positive(mass / rho1, "proposed material volume is unrepresentable")?;
            let defect = (v0 - domain.volume(c)).abs();
            maximum_initial_volume_defect_m3 = maximum_initial_volume_defect_m3.max(defect);
            if defect > controls.maximum_initial_volume_defect_m3 {
                return Err(Error::Invalid(
                    "original material mass and density do not fill the occupied voxel",
                ));
            }
            let change = finite(
                v0.max(v1) * ((rho0 - rho1) / rho0.max(rho1)),
                "thermal material volume change is unrepresentable",
            )?;
            if rho0 != rho1 && change == 0.0 {
                return Err(Error::Invalid("thermal material volume change underflows"));
            }
            sources[c] = finite(change / dt, "thermal volume rate is unrepresentable")?;
            if change != 0.0 && sources[c] == 0.0 {
                return Err(Error::Invalid("thermal volume rate underflows"));
            }
            cells.insert(
                entry.coordinates,
                MaterialFreeSurfaceVolumeCell3 {
                    mass_kg: mass,
                    initial: entry.state,
                    proposed: next,
                    initial_material,
                    proposed_material,
                    requested_enthalpy_increment_j: requested,
                    represented_enthalpy_increment_j: represented,
                    initial_volume_m3: v0,
                    proposed_volume_m3: v1,
                    volume_change_m3: change,
                },
            );
        }
        let mut openings = Vec::new();
        for c in 0..count {
            poll(gate)?;
            if !domain.is_fluid(c) {
                continue;
            }
            for face in Face3::ALL {
                let Some(n) = domain.neighbor(c, face as usize) else {
                    continue;
                };
                if !domain.is_fluid(n)
                    && !self
                        .lattice
                        .transport_wall_at(domain.coords(n).map(|v| v as i64))
                {
                    openings.push(InternalPressureFace {
                        cell: c,
                        face,
                        pressure_pa: gas_pressure(domain.coords(c), face),
                    });
                }
            }
        }
        let faces = [FlowFace::Wall; 6];
        let predictor = FlowField::from_cell_velocity_predictor_with_internal_pressure_faces(
            &domain,
            &velocities,
            faces,
            &openings,
            gate,
        )?;
        let projection = predictor.project_density_volume_source_with_internal_pressure_faces(
            &domain,
            faces,
            &openings,
            &densities,
            &sources,
            controls.projection,
            gate,
        )?;
        poll(gate)?;
        Ok(MaterialFreeSurfaceVolumeProjection3 {
            domain,
            projection,
            cells,
            volume_source_m3_s: sources,
            gas_pressure_faces: openings,
            maximum_initial_volume_defect_m3,
            maximum_source_density_defect_kg_m3,
            absolute_enthalpy_defect_j,
        })
    }
}

fn poll(gate: &CancelGate) -> Result<(), Error> {
    if gate.is_requested() {
        Err(crate::d3q19::freesurface3::FreeSurfaceError3::Cancelled.into())
    } else {
        Ok(())
    }
}

fn source_state(
    law: &LiquidTemperatureLaw<'_>,
    phase: EquilibriumPhaseState,
    cache: &mut BTreeMap<u64, Arc<ResolvedLiquidState>>,
) -> Result<(Arc<ResolvedLiquidState>, f64), Error> {
    if phase.phase() != SolidLiquidPhase::Liquid {
        return Err(Error::Invalid(
            "free-surface volume projection needs an explicit liquid mobility state",
        ));
    }
    let key = phase.temperature_k().to_bits();
    if let std::collections::btree_map::Entry::Vacant(entry) = cache.entry(key) {
        entry.insert(Arc::new(
            law.at_temperature(Temperature::new(phase.temperature_k()))?,
        ));
    }
    let source = &cache[&key];
    let defect = finite(
        (source.density_kg_m3() - phase.bulk_density_kg_m3()).abs(),
        "caloric/source density discrepancy is unrepresentable",
    )?;
    Ok((source.clone(), defect))
}
