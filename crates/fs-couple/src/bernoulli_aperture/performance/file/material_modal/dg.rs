//! Authored material plates and an original causal DG bore on one source clock.
//!
//! The material/mesh, loss, wetted-face, spatial-pressure, closure and force
//! records are the same owners used by the two-sided material-modal source.
//! Only the acoustic host changes. Each requested signed mode and material
//! history remains in the original ExactZoh ModalAperture.
//!
//! Record order: schema; audio RATE RATE SAMPLES FULL_SCALE_PA;
//! propagator exact-zoh; ambient GAS_RECORD; shared dg_tube/dg_bore and its
//! geometry/gas/radiation records; slit_length none|LENGTH; observation
//! inlet-pressure|baffled-outlet X Y Z RINGS ANGLES FMAX_HZ; valve REST_GAP;
//! limits; contact; loss; blades and original plate records; optional
//! wetted_surfaces, pressure_shapes, closure, source_forces, source_pressures;
//! compile_limits; schedule and the canonical blowing-pressure gesture suffix.
//!
//! Audio frames are mechanical samples. Explicit downstream decimation retains
//! the pressure owner and its reported delay. A baffled observer requires the
//! actual computed DG terminal and its independently declared exterior gas.

use super::{
    BladePropagator, BladeState, LossSelection, MaterialCache,
    MaterialModalValveInfo, ModalAperture, ModalValveMaterialRecord,
    NormalizedMaterialCardPack, Obstacle, PairedApertureLimits,
    PlateValveInputError, Reader, MAX_CONTROLS, MAX_PLATE_VALVE_PERFORMANCE_BYTES,
    bad, checkpoint, closure, forces, gas, optional_positive, plate, pressure,
    read_propagator, read_receiver, wetted,
};
use crate::bernoulli_aperture::dg::LoadedApertureDg;
use crate::bernoulli_aperture::dg::loaded::performance::LoadedDgPerformance;
use crate::bernoulli_aperture::performance::{ApertureObservation, AperturePerformanceConfig};
use crate::dg_port::input::DgSourceInput;
use crate::dg_port::input::wall::DgMaterialWallRecord;
use crate::dg_port::input::wall::regions::{DgMaterialWallLayout, DgMaterialWallLayoutInput};
use fs_blake3::hash_domain;
use fs_exec::CancelGate;
use fs_material::gas::GasState;
use fs_material::gas_species::input::{GasRecord, ThermochemicalSources};

/// Exact selector for material modal mechanics coupled to the original DG bore.
pub const MATERIAL_MODAL_DG_VALVE_SCHEMA: &str = "frankensim-material-modal-dg-valve-v1";

/// Original source records and the finite renderer; no duplicate mechanical
/// states, radiation fits, schedules or material coefficients are stored here.
pub struct MaterialModalDgValvePerformance {
    info: MaterialModalValveInfo,
    renderer: LoadedDgPerformance,
    materials: Vec<ModalValveMaterialRecord>,
    gases: Vec<GasRecord>,
    ambient: GasState,
    wall: DgMaterialWallLayout,
}

impl MaterialModalDgValvePerformance {
    /// Bind all authored material modes to the actual DG inlet gas and mean
    /// pressure/flow impedance, then compile the original mechanical controls.
    ///
    /// The caller resolves normalized material artifacts; this library performs
    /// no file I/O. Every additional Maxwell, spatial traction, contact and
    /// viscous-slit choice is admitted by its original material/physics owner.
    pub fn from_bytes_with_sources(
        bytes: &[u8],
        max_block: usize,
        gate: &CancelGate,
        gas_sources: Option<&ThermochemicalSources<'_>>,
        mut load_material: impl FnMut(&str) -> Result<NormalizedMaterialCardPack, String>,
    ) -> Result<Self, PlateValveInputError> {
        checkpoint(gate)?;
        if bytes.len() > MAX_PLATE_VALVE_PERFORMANCE_BYTES
            || !(1..=65_536).contains(&max_block)
        {
            return Err(bad(0, "material-modal DG source exceeds its byte or callback budget"));
        }
        let text = core::str::from_utf8(bytes).map_err(|_| bad(0, "input must be UTF-8"))?;
        let (header, gestures) = text.split_once("\nschedule\n")
            .ok_or_else(|| bad(0, "missing LF-delimited schedule and canonical gesture suffix"))?;
        let mut rows = Reader { lines: header.lines(), line: 0 };
        rows.row(MATERIAL_MODAL_DG_VALVE_SCHEMA)?.finish()?;
        let mut row = rows.row("audio")?;
        let mechanics: u32 = row.parse()?;
        let output: u32 = row.parse()?;
        let samples: u64 = row.parse()?;
        let full_scale_pa = row.scalar()?;
        row.finish()?;
        if !(8_000..=3_072_000).contains(&mechanics) || output != mechanics
            || samples == 0 || samples > 28_800_000
            || samples > 600 * u64::from(mechanics) || full_scale_pa <= 0.0
        {
            return Err(bad(rows.line,
                "DG audio requires one 8..3072 kHz source/output clock, positive PCM pascals and a bounded finite window"));
        }
        let propagator = read_propagator(&mut rows)?;
        if propagator != BladePropagator::ExactZoh {
            return Err(bad(rows.line, "material-modal DG mechanics require explicit propagator exact-zoh"));
        }
        let dt = 1.0 / f64::from(mechanics);
        let mut gases = Vec::new();
        let mut row = rows.row("ambient")?;
        let ambient = gas::read_with_sources(&mut row, None, gas_sources, &mut gases)?;
        row.finish()?;
        let kind = if rows.lines.clone().next().and_then(|line|
            line.split_ascii_whitespace().next()) == Some("dg_tube")
        {
            "dg_tube"
        } else { "dg_bore" };
        let mut row = rows.row(kind)?;
        let duct = DgSourceInput::read(
            kind, row.line, &mut row.fields,
            &mut || {
                let row = rows.next()?;
                Ok::<_, PlateValveInputError>((row.line, row.fields))
            },
            &ambient, gas_sources, &mut gases, gate,
        )?;
        row.finish()?;
        let mut cache = MaterialCache::default();
        let wall = if rows.lines.clone().next().and_then(|line|
            line.split_ascii_whitespace().next()) == Some("dg_wall")
        {
            let mut row = rows.row("dg_wall")?;
            let declaration = DgMaterialWallLayoutInput::read(
                row.line, &mut row.fields, &mut || {
                    let row = rows.next()?;
                    Ok::<_, PlateValveInputError>((row.line, row.fields))
                },
            )?;
            row.finish()?;
            declaration.resolve(gate, &mut |path, line|
                cache.load(path, line, &mut load_material, gate))?
        } else { DgMaterialWallLayout::default() };
        let mut row = rows.row("slit_length")?;
        let slit = optional_positive(&mut row)?;
        row.finish()?;
        let mut row = rows.row("observation")?;
        let observation = match row.word()? {
            "inlet-pressure" => ApertureObservation::Inlet,
            "baffled-outlet" => ApertureObservation::TubeBaffled(read_receiver(&mut row)?),
            _ => return Err(bad(row.line,
                "DG observation requires inlet-pressure or baffled-outlet X Y Z RINGS ANGLES FMAX_HZ")),
        };
        row.finish()?;
        let mut row = rows.row("valve")?;
        let rest_opening = row.scalar()?;
        row.finish()?;
        let mut row = rows.row("limits")?;
        let limits = PairedApertureLimits {
            maximum_abs_displacement_m: row.scalar()?,
            maximum_abs_velocity_m_s: row.scalar()?,
            maximum_energy_j: row.scalar()?,
            max_evaluations: row.count(65_536)?,
            velocity_tolerance_m_s: row.scalar()?,
            relative_tolerance: row.scalar()?,
        };
        row.finish()?;
        let mut row = rows.row("contact")?;
        let stiffness = row.scalar()?;
        let alpha = row.scalar()?;
        let contact_loss = row.scalar()?;
        row.finish()?;
        let contact = Obstacle::new(
            vec![-1.0], 1, 1, vec![0.0], vec![1.0], stiffness, alpha,
            "material-modal DG file: declared effective uniform slit contact".into(),
        ).and_then(|contact| contact.with_internal_loss(contact_loss))
            .map_err(|error| PlateValveInputError::Physics(
                crate::acoustic_realize::AcousticRealizeError::Nonlinear(error.to_string()),
            ))?;
        let loss = LossSelection::read(&mut rows)?;
        let blades = rows.count("blades", 2)?;
        if blades == 0 {
            return Err(bad(rows.line, "one or two independently supplied material blades are required"));
        }
        let mut plates = Vec::with_capacity(blades);
        let mut states = Vec::<BladeState>::new();
        let mut materials = Vec::new();
        let (mut nodes, mut triangles) = (0, 0);
        for blade in 0..blades {
            let (source, options, initial, records) = plate::read(
                &mut rows, blade, states.len(), &loss, &mut cache, &mut load_material, gate,
            )?;
            nodes += source.chart().mesh.nodes.len();
            triangles += source.chart().mesh.tris.len();
            states.extend(initial);
            materials.extend(records);
            plates.push((source, options));
        }
        let modes = states.len();
        let wetted_surfaces = wetted::read(&mut rows, blades, gate)?;
        let pressure_shapes = pressure::shapes(&mut rows, blades, gate)?;
        let closure = closure::read(&mut rows, gate)?;
        let forces = forces::read(&mut rows, blades, samples, gate)?;
        let pressure_events = pressure::events(&mut rows, blades, samples, gate)?;
        let mut row = rows.row("compile_limits")?;
        let max_compile_work: u64 = row.parse()?;
        let max_controls = row.count(MAX_CONTROLS)?;
        row.finish()?;
        if max_compile_work > 16_777_216 {
            return Err(bad(rows.line, "pressure compilation exceeds 16777216 visits"));
        }
        if rows.lines.next().is_some() {
            return Err(bad(rows.line + 1, "unexpected trailing material-modal DG record"));
        }
        let schedule = super::super::decode_schedule(gestures)?;
        checkpoint(gate)?;
        let load = duct.prepare_with_material_wall_layout(&ambient, dt, gate, &wall)?;
        let inlet = *load.duct().inlet_gas();
        let impedance = load.impedance_pa_s_m3();
        let mut prepared_plates = Vec::with_capacity(plates.len());
        for (source, options) in plates {
            let tensor = if let LossSelection::Tensor { properties, limits: projection } = &loss {
                let law = crate::plate_viscosity::PlateViscosity::from_chart(
                    &source, properties.each_ref().map(String::as_str), gate,
                ).map_err(PlateValveInputError::Physics)?;
                Some((law, *projection))
            } else { None };
            prepared_plates.push((source, options, tensor));
        }
        let mut valve = ModalAperture::from_material_images_with_wetted(
            prepared_plates, wetted_surfaces, rest_opening, inlet, impedance, dt,
            samples, limits, states, contact, propagator, gate,
        ).map_err(PlateValveInputError::Physics)?;
        if let LossSelection::Maxwell { properties, initial, regional } = loss {
            valve = if regional {
                valve.with_regional_material_relaxation(&properties, initial, gate)
            } else {
                valve.with_material_relaxation(&properties, initial, gate)
            }.map_err(PlateValveInputError::Physics)?;
        }
        for pressure in pressure_shapes {
            valve = pressure.bind(valve, gate)?;
        }
        if let Some(closure) = closure {
            valve = valve.with_plate_closure(closure).map_err(PlateValveInputError::Physics)?;
        }
        if let Some(length) = slit {
            valve.set_material_viscous_slit(length).map_err(PlateValveInputError::Physics)?;
        }
        let system = LoadedApertureDg::new_loaded(valve, load)
            .map_err(PlateValveInputError::Physics)?;
        let mut renderer = LoadedDgPerformance::new(system, observation, schedule,
            AperturePerformanceConfig {
                sample_rate_hz: mechanics, samples, max_block, max_compile_work, max_controls,
            },
        ).map_err(PlateValveInputError::Gesture)?;
        for (blade, events) in forces.into_iter().enumerate() {
            if !events.is_empty() {
                renderer = renderer.with_material_plate_force_controls(blade, events)
                    .map_err(PlateValveInputError::Gesture)?;
            }
        }
        for (blade, events) in pressure_events.into_iter().enumerate() {
            if !events.is_empty() {
                renderer = renderer.with_material_pressure_controls(blade, events)
                    .map_err(PlateValveInputError::Gesture)?;
            }
        }
        checkpoint(gate)?;
        Ok(Self {
            info: MaterialModalValveInfo {
                input_hash: hash_domain("org.frankensim.fs-couple.material-modal-dg-valve.v1", bytes),
                mechanics_rate_hz: mechanics, propagator, output_rate_hz: output,
                samples, full_scale_pa, blades, modes, nodes, triangles,
            },
            renderer, materials, gases, ambient, wall,
        })
    }

    /// Exact original source dimensions, mode counts, physical clock and PCM scale.
    pub const fn info(&self) -> MaterialModalValveInfo { self.info }
    /// Original modal mechanics, material memory, DG bore and pressure observer.
    pub const fn renderer(&self) -> &LoadedDgPerformance { &self.renderer }
    /// Every original normalized material artifact and regional address.
    pub fn material_records(&self) -> &[ModalValveMaterialRecord] { &self.materials }
    /// Ambient, prescribed local and computed exterior material receipts.
    pub fn gas_records(&self) -> &[GasRecord] { &self.gases }
    /// Original material wall source, shared with any plate using the same pack.
    pub const fn wall_record(&self) -> Option<&DgMaterialWallRecord> { self.wall.uniform_record() }
    /// Original physical wall regions, retaining the shared source artifacts.
    pub const fn wall_layout(&self) -> &DgMaterialWallLayout { &self.wall }
    /// Declared default gas; the actual DG inlet/field is retained by the renderer.
    pub const fn ambient_gas(&self) -> &GasState { &self.ambient }
    /// Move the complete accepted physical and scheduling state without a reset.
    pub fn into_renderer(self) -> LoadedDgPerformance { self.renderer }
}
