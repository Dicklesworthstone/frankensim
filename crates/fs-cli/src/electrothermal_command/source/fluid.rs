//! Scene-owned fluid role and prescribed conservative transport. Material
//! properties still resolve through the common state-point and transport owners.

use fs_exec::CancelGate;
use fs_lbm::conjugate::{ChtError, FlowFace, FlowField, FluidProperties, ThermalFace, ThermalSetup, Voxel, VoxelDomain};
use fs_lbm::material::conjugate::MaterialChtFluid;
use fs_matdb::{ClaimSelection, PropertyKey, SelectionPolicy};
use fs_material::fluid_transport::FluidTransportState;
use fs_material::liquid::{DYNAMIC_VISCOSITY_DIMS, DYNAMIC_VISCOSITY_PROPERTY};
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, ResolvedMaterialStatePoint,
    SPECIFIC_HEAT_CAPACITY_DIMS, SPECIFIC_HEAT_CAPACITY_PROPERTY,
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_qty::Density;
use super::{Source, invalid, poll};
use super::super::{J, array, fields, integer, number, positive, quote, required, values};

pub(crate) const NO_CLAIM: &str = "steady fixed-geometry scalar ohmic and isotropic or grid-aligned orthotropic solid thermal conduction with prescribed volume-conservative isotropic fluid transport; density, heat capacity and viscosity frozen at an explicit sourced reference state, local conductivity and resistivity re-resolved; no momentum/pressure solve, variable-density enthalpy transport, phase change, electrochemistry, Lorentz force, AC effects, propagated uncertainty or experimental qualification";

pub(super) const TRANSIENT_NO_CLAIM: &str = "fixed-geometry backward-Euler electrothermal evolution with prescribed volume-conservative fluid velocity; fluid density verified constant over the complete source interval, one sourced Cp/enthalpy chart for storage and first-order upwind face transport, viscosity frozen at the sourced reference, local conductivity and resistivity re-resolved; single-phase fixed cell masses, sourced fluid and solid Cp and exterior gray radiation; no momentum/pressure solve, variable-density fluid transport, phase change, deformation, electrochemistry, Lorentz force, AC effects, propagated uncertainty or experimental qualification; temporal accuracy requires time-step refinement";

/// One fluid material, its exact reference receipt and the actual energy fluxes.
/// No material identity is inferred from a mechanical occupancy or phase name.
pub(crate) struct Fluid {
    pub material_index: usize,
    reference_temperature_k: f64,
    reference: ResolvedMaterialStatePoint,
    pub properties: FluidProperties,
    pub flow: FlowField,
    velocity_m_s: [f64; 3],
    divergence_m3_s: f64,
}

/// Omission preserves the all-solid scene. The optional role assigns exactly
/// one existing material to the existing single-fluid heat-advection model.
pub(crate) fn load_fluid(
    root: &J,
    domain: &mut VoxelDomain,
    material_indices: &[usize],
    sources: &[Source],
) -> Result<Option<Fluid>, String> {
    let Some(input) = root.get("fluid") else { return Ok(None); };
    fields(input, &["material_index", "reference_temperature_k", "velocity_m_s", "flow_faces"])?;
    let material_index = integer(required(input, "material_index")?)?;
    let source = sources.get(material_index).ok_or("fluid names an unknown material")?;
    if material_indices.len() != domain.cell_count() || !material_indices.contains(&material_index) {
        return Err("fluid material must be assigned to at least one cell".into());
    }
    let reference_temperature_k = positive(input, "reference_temperature_k")?;
    let velocity = array(input, "velocity_m_s", 3)?;
    if velocity.len() != 3 { return Err("fluid velocity needs three components in m/s".into()); }
    let velocity_m_s = [number(&velocity[0])?, number(&velocity[1])?, number(&velocity[2])?];
    let rules = array(input, "flow_faces", 6)?;
    if rules.len() != 6 { return Err("fluid needs six flow faces in x-,x+,y-,y+,z-,z+ order".into()); }
    let mut faces = [FlowFace::Wall; 6];
    for (face, rule) in faces.iter_mut().zip(rules) {
        *face = match rule.as_str() {
            Some("wall") => FlowFace::Wall,
            Some("open") => FlowFace::Fixed,
            _ => return Err("fluid flow faces must be wall or open (prescribed flux)".into()),
        };
    }
    // The native reference transport is isotropic. Orthotropic solids can
    // coexist with it, but a directional fluid must not lose its tensor here.
    let mut requirements = vec![source.conductivity.isotropic_requirement()
        .map_err(|e| e.to_string())?.clone()];
    for (name, dims) in [
        (DENSITY_PROPERTY, Density::DIMS),
        (SPECIFIC_HEAT_CAPACITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_DIMS),
        (DYNAMIC_VISCOSITY_PROPERTY, DYNAMIC_VISCOSITY_DIMS),
    ] {
        requirements.push(ScalarPropertyRequirement::try_new(
            name, dims, ScalarAdmissibility::StrictlyPositive,
        ).map_err(|e| e.to_string())?);
    }
    // Reuse the source's exact conductivity quantity and state coordinates,
    // including its pressure/composition, instead of stripping semantic kinds.
    let point = source.point_at(reference_temperature_k).map_err(|e| e.to_string())?;
    let reference = resolve_material_state_point(source.pack.card(), &point, &requirements,
        MaterialPropertySelection::SingleClaimOnly).map_err(|e| format!("fluid reference: {e}"))?;
    let value = |key| reference.property(key).expect("required fluid property was resolved").value_si();
    let transport = FluidTransportState::try_new(
        value(DENSITY_PROPERTY), value(SPECIFIC_HEAT_CAPACITY_PROPERTY),
        value(super::THERMAL_CONDUCTIVITY_PROPERTY), value(DYNAMIC_VISCOSITY_PROPERTY),
    ).map_err(|e| e.to_string())?;
    let properties = MaterialChtFluid::try_new(&transport).map_err(|e| e.to_string())?
        .properties().clone();

    // Geometry stays fixed and all other material roles stay solid. Electrical
    // electrodes do not implicitly create mechanical walls or block flow.
    let mut candidate = domain.clone();
    for (cell, &index) in material_indices.iter().enumerate() {
        if index == material_index {
            let [x, y, z] = candidate.coords(cell);
            candidate.set(x, y, z, Voxel::Fluid);
        }
    }
    let area = candidate.dx() * candidate.dx();
    if !area.is_finite() || area <= 0.0 || !candidate.volume(0).is_finite() || candidate.volume(0) <= 0.0 {
        return Err("fluid face area or cell volume is unrepresentable".into());
    }
    let mut flux_scale = 0.0;
    for speed in velocity_m_s {
        let flux = speed.abs() * area;
        if !flux.is_finite() || (speed != 0.0 && flux == 0.0) {
            return Err("fluid prescribed volume flux is unrepresentable".into());
        }
        flux_scale += flux;
    }
    if !flux_scale.is_finite() { return Err("fluid flux scale is unrepresentable".into()); }
    let flow = FlowField::from_face_velocity(&candidate, faces, |axis, _| velocity_m_s[axis])
        .map_err(|e| e.to_string())?;
    let divergence_m3_s = flow.max_divergence(&candidate);
    if !divergence_m3_s.is_finite() || divergence_m3_s > 128.0 * f64::EPSILON * flux_scale {
        return Err("prescribed fluid velocity is not volume-conservative around walls/obstacles; supply a compatible straight-through field or zero velocity".into());
    }
    *domain = candidate;
    Ok(Some(Fluid { material_index, reference_temperature_k, reference,
        properties, flow, velocity_m_s, divergence_m3_s }))
}

impl Fluid {
    /// Constant density preserves reference control-volume mass. Cp may vary:
    /// the transient owner builds ONE sourced caloric chart for both storage
    /// and conservative face enthalpy. Full source-line queries retain every
    /// crossed breakpoint and refuse holes, switching or ambiguous claims.
    pub(super) fn admit_transient(
        &self, source: &Source, lower_k: f64, upper_k: f64,
        thermal: &ThermalSetup, gate: &CancelGate,
    ) -> Result<(), ChtError> {
        poll(gate)?;
        if self.reference.card_identity() != source.pack.card().content_hash() {
            return Err(invalid("fluid.transient", "transport and storage must use the same material card"));
        }
        let lower = source.point_at(lower_k)?;
        let upper = source.point_at(upper_k)?;
        for (name, dims) in [
            (DENSITY_PROPERTY, Density::DIMS),
            (SPECIFIC_HEAT_CAPACITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_DIMS),
        ] {
            poll(gate)?;
            let reference = self.reference.property(name).expect("fluid reference required rho and Cp");
            let line = source.pack.card().claims().query_line_typed(
                &PropertyKey::new(name, dims), &source.temperature_axis,
                &lower, &upper, ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
            ).map_err(|e| invalid("fluid.transient", format!("{name} requires continuous source support: {e:?}")))?;
            for (_, answer) in &line.knots {
                poll(gate)?;
                if answer.receipt.selected != reference.answer().receipt.selected {
                    return Err(invalid("fluid.transient", format!(
                        "{name} must retain its original reference source claim throughout the interval",
                    )));
                }
                let value = answer.evidence.value.value;
                if !value.is_finite() || value <= 0.0 {
                    return Err(invalid("fluid.transient", format!("{name} must remain finite and positive")));
                }
                if name == DENSITY_PROPERTY && value != reference.value_si() {
                    return Err(invalid("fluid.transient", "density must remain constant throughout the complete source interval"));
                }
            }
        }
        for face in thermal.faces {
            let temperature = match face {
                ThermalFace::Inflow { temperature } => temperature,
                ThermalFace::Outflow { backflow_temperature } => backflow_temperature,
                _ => continue,
            };
            if !temperature.is_finite() || temperature < lower_k || temperature > upper_k {
                return Err(invalid("fluid.transient", "inflow and backflow temperatures must lie inside the sourced storage interval"));
            }
        }
        poll(gate)
    }

    /// This is reference-state evidence; transient calorics retain their own chart.
    pub fn render(&self) -> String {
        format!(concat!(
            "{{\"material_index\":{},\"temperature_k\":{},\"material_state\":{},",
            "\"density_kg_m3\":{},\"specific_heat_j_kg_k\":{},\"kinematic_viscosity_m2_s\":{},",
            "\"velocity_m_s\":{},\"max_divergence_m3_s\":{}}}"),
            self.material_index, self.reference_temperature_k, quote(&self.reference.identity().to_hex()),
            self.properties.density_kg_m3, self.properties.specific_heat_j_kg_k,
            self.properties.kinematic_viscosity_m2_s, values(self.velocity_m_s), self.divergence_m3_s)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::electrothermal_command::source::{Problem, Source};
    use fs_blake3::hash_bytes;
    use fs_evidence::ValidityDomain;
    use fs_exec::CancelGate;
    use fs_lbm::conjugate::{ElectricSetup, Electrode, ElectrodeDrive, EnergyConfig, ThermalFace, ThermalSetup};
    use fs_lbm::conjugate::material_coupling::MaterialCouplingConfig;
    use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialStateId, NormalizedPack,
        NormalizedMaterialCardPack, PropertyClaim, PropertyKey, PropertyValue, Provenance,
        QueryPoint, UncertaintyModel};
    use fs_material::conductor::{ELECTRICAL_RESISTIVITY_DIMS, ELECTRICAL_RESISTIVITY_PROPERTY};
    use fs_material::state_point::{THERMAL_CONDUCTIVITY_DIMS, THERMAL_CONDUCTIVITY_PROPERTY};
    use fs_qty::QuantitySpec;
    use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};

    fn source(cp: f64, resistivity: f64, conductivity: f64, typed: bool) -> Source {
        material_source(cp, resistivity, conductivity, typed, "liquid", |_| {})
    }
    fn material_source(
        cp: f64, resistivity: f64, conductivity: f64, typed: bool, phase: &str,
        mut edit: impl FnMut(&mut PropertyClaim),
    ) -> Source {
        let mut claims = ClaimSet::new();
        for (key, value, dims) in [
            (DENSITY_PROPERTY, 1.0, Density::DIMS),
            (SPECIFIC_HEAT_CAPACITY_PROPERTY, cp, SPECIFIC_HEAT_CAPACITY_DIMS),
            (DYNAMIC_VISCOSITY_PROPERTY, 1.0, DYNAMIC_VISCOSITY_DIMS),
            (THERMAL_CONDUCTIVITY_PROPERTY, conductivity, THERMAL_CONDUCTIVITY_DIMS),
            (ELECTRICAL_RESISTIVITY_PROPERTY, resistivity, ELECTRICAL_RESISTIVITY_DIMS),
        ] {
            let quantity = if typed && key == THERMAL_CONDUCTIVITY_PROPERTY {
                QuantitySpec::semantic(SemanticType::new(QuantityKind::ThermalConductivity, ValueForm::Static))
            } else { QuantitySpec::dimensional(dims) };
            let mut claim = PropertyClaim {
                key: PropertyKey::with_quantity(key, quantity), value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained().with("T", 290.0, 350.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity, observations: vec![],
                provenance: Provenance { source: "synthetic conducting fluid".into(),
                    license: "test-only".into(), artifact: None },
            };
            edit(&mut claim);
            claims.insert_claim(claim).unwrap();
        }
        crate::electrothermal_command::tests::authored_observation(&mut claims, "synthetic conducting fluid");
        let pack = NormalizedPack::new("fluid-fixture", "test-v1", hash_bytes(b"synthetic fluid"),
            "test-only", claims, vec![], vec![]).unwrap();
        let pack = NormalizedMaterialCardPack::new(MaterialStateId {
            chemistry: "synthetic".into(), phase: phase.into(), process: "authored test".into(), revision: 0,
        }, pack).unwrap();
        // Exercise the ordinary pack decoder, not an alternate data path.
        Source::new(NormalizedMaterialCardPack::from_bytes(&pack.to_bytes()).unwrap(),
            QueryPoint::new(), "T".into(), None, true).unwrap()
    }
    fn input(speed: f64) -> J {
        J::parse(&format!(r#"{{"fluid":{{"material_index":0,"reference_temperature_k":300,
            "velocity_m_s":[{speed},0,0],"flow_faces":["open","open","wall","wall","wall","wall"]}}}}"#)).unwrap()
    }
    fn problem(cp: f64, rho: f64, k: f64, speed: f64, typed: bool) -> Problem {
        let mut domain = VoxelDomain::from_fn(3, 1, 1, 1.0, |_| Voxel::Solid(0)).unwrap();
        let materials = vec![source(cp, rho, k, typed)];
        let indices = vec![0; 3];
        let fluid = load_fluid(&input(speed), &mut domain, &indices, &materials).unwrap();
        let mut faces = [ThermalFace::Adiabatic; 6];
        faces[0] = ThermalFace::Temperature(300.0);
        faces[1] = ThermalFace::Temperature(300.0);
        if speed != 0.0 {
            faces[0] = ThermalFace::Inflow { temperature: 300.0 };
            faces[1] = ThermalFace::Outflow { backflow_temperature: 300.0 };
        }
        Problem { domain, material_indices: indices, materials, fluid,
            electrical_interfaces: Vec::new(), thermal_interfaces: Vec::new(),
            thermal_contact_tolerance: 1e-10, radiation: None,
            electric: ElectricSetup { conductors: vec![], electrodes: vec![
                Electrode { lo: [0,0,0], hi: [1,1,1], drive: ElectrodeDrive::Voltage(0.0) },
                Electrode { lo: [2,0,0], hi: [3,1,1], drive: ElectrodeDrive::Voltage(3.0) },
            ] }, thermal: ThermalSetup::new(faces), initial_temperature: vec![300.0; 3],
            coupling: MaterialCouplingConfig { max_passes: 100, temperature_tolerance_k: 1e-9,
                power_tolerance_w: 1e-9, relative_tolerance: 1e-10 },
            energy: EnergyConfig::default(), max_thermal_passes: 100, conductivity_tolerance: 1e-10,
        }
    }

    #[test]
    fn source_only_electrical_and_thermal_swap_changes_fluid_current_and_temperature() {
        for typed in [false, true] {
            for (rho, k, t, power) in [(2.0,1.0,303.375,4.5), (4.0,2.0,300.84375,2.25)] {
                let mut p = problem(1.0, rho, k, 0.0, typed);
                let result = p.solve(&CancelGate::new()).unwrap();
                assert!((result.coupled.thermal.temperature[1]-t).abs() < 1e-7);
                assert!((result.coupled.electric.total_joule_w-power).abs() < 1e-8);
                assert!((result.coupled.electric.delivered_w-power).abs() < 1e-8);
                assert!((0..3).all(|cell| p.domain.is_fluid(cell)));
                let rendered = J::parse(&crate::electrothermal_command::render(&p, &result, true)).unwrap();
                assert!(rendered.get("fluid_reference").unwrap().str_field("material_state").is_some());
            }
        }
    }

    #[test]
    fn sourced_heat_capacity_changes_advective_cooling_without_changing_prescribed_flow() {
        let mut a = problem(1.0, 2.0, 1.0, 1.0, false);
        let mut b = problem(4.0, 2.0, 1.0, 1.0, false);
        let before = b.fluid.as_ref().unwrap().reference.clone();
        let ra = a.solve(&CancelGate::new()).unwrap();
        let rb = b.solve(&CancelGate::new()).unwrap();
        let peak = |r: &super::super::ResultField| r.coupled.thermal.temperature.iter().copied().fold(0.0, f64::max);
        assert!(peak(&rb) < peak(&ra));
        assert_eq!(a.fluid.as_ref().unwrap().flow, b.fluid.as_ref().unwrap().flow);
        assert_eq!(ra.coupled.electric.total_joule_w, rb.coupled.electric.total_joule_w);
        assert_eq!(b.fluid.as_ref().unwrap().reference, before);
        assert!(rb.coupled.thermal.report.balance.relative_residual < 1e-7);
    }

    #[test]
    fn incompatible_flow_or_unsupported_reference_refuses_without_changing_occupancy() {
        let mut domain = VoxelDomain::from_fn(3, 1, 1, 1.0, |_| Voxel::Solid(0)).unwrap();
        let original = domain.clone();
        let sources = vec![source(1.0, 2.0, 1.0, false)];
        let base = r#"{"fluid":{"material_index":0,"reference_temperature_k":300,
            "velocity_m_s":[1,0,0],"flow_faces":["open","open","wall","wall","wall","wall"]}}"#;
        for bad in [base.replace("\"material_index\":0", "\"material_index\":9"),
            base.replace("\"reference_temperature_k\":300", "\"reference_temperature_k\":400"),
            base.replace("[\"open\",\"open\"", "[\"wall\",\"open\""),
            base.replace("[1,0,0]", "[1,0]"),
        ] {
            assert!(load_fluid(&J::parse(&bad).unwrap(), &mut domain, &[0;3], &sources).is_err());
            assert_eq!(domain, original);
        }
        assert!(load_fluid(&J::parse("{}").unwrap(), &mut domain, &[0;3], &sources).unwrap().is_none());
        assert_eq!(domain, original);
    }

    #[test]
    fn scene_file_reaches_both_steady_and_transient_fluid_solvers() {
        use std::sync::atomic::{AtomicU64, Ordering};
        static NEXT: AtomicU64 = AtomicU64::new(0);
        let stamp = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos();
        let dir = std::env::temp_dir().join(format!("frankensim-fluid-cli-{}-{stamp}-{}",
            std::process::id(), NEXT.fetch_add(1, Ordering::Relaxed)));
        std::fs::create_dir(&dir).unwrap();
        let bytes = source(1.0, 2.0, 1.0, false).pack.to_bytes();
        std::fs::write(dir.join("fluid.fsmcdpk"), &bytes).unwrap();
        let scene = r#"{"schema":"frankensim.electrothermal.v1","cells":[3,1,1],"cell_size_m":1,
            "materials":[{"pack":"fluid.fsmcdpk","conducting":true,"temperature_axis":"T",
                "temperature_kind":"legacy","state":[]}],
            "material_indices":[0,0,0],"initial_temperature_k":300,
            "fluid":{"material_index":0,"reference_temperature_k":300,"velocity_m_s":[1,0,0],
                "flow_faces":["open","open","wall","wall","wall","wall"]},
            "thermal_faces":[{"inflow_temperature_k":300},{"backflow_temperature_k":300},"adiabatic","adiabatic","adiabatic","adiabatic"],
            "electrodes":[{"lo":[0,0,0],"hi":[1,1,1],"voltage_v":0},
                {"lo":[2,0,0],"hi":[3,1,1],"voltage_v":3}],
            "numerics":{"max_passes":100,"max_thermal_passes":100,"temperature_tolerance_k":1e-9,
                "power_tolerance_w":1e-9,"relative_tolerance":1e-10,"conductivity_relative_tolerance":1e-10,
                "energy_tolerance":1e-12,"energy_max_iterations":2000}}"#;
        let path = dir.join("scene.json");
        std::fs::write(&path, scene).unwrap();
        let run = crate::electrothermal_command::run;
        let result = run(&[path.as_os_str().to_owned()], true);
        assert_eq!(result.exit_code, crate::exit::SUCCESS, "{}", result.stderr);
        let value = J::parse(&result.stdout).unwrap();
        assert!((value.f64_field("total_joule_w").unwrap() - 4.5).abs() < 1e-8);
        assert!(value.get("fluid_reference").unwrap().str_field("material_state").is_some());
        assert!(value.str_field("no_claim").unwrap().contains("prescribed"));
        let temperatures = value.get("temperature_k").unwrap().as_array().unwrap();
        assert!(temperatures[2].as_f64().unwrap() > temperatures[0].as_f64().unwrap());
        assert_eq!(result, run(&[path.as_os_str().to_owned()], true));
        assert_eq!(std::fs::read_to_string(&path).unwrap(), scene);
        assert_eq!(std::fs::read(dir.join("fluid.fsmcdpk")).unwrap(), bytes);
        let time_path = dir.join("time.json");
        std::fs::write(&time_path, transient::request_text()).unwrap();
        let args = [path.as_os_str().to_owned(), "--transient".into(), time_path.into_os_string()];
        let transient = run(&args, true);
        assert_eq!(transient.exit_code, crate::exit::SUCCESS, "{}", transient.stderr);
        let value = J::parse(&transient.stdout).unwrap();
        assert_eq!(value.f64_field("steps"), Some(8.0));
        assert_eq!(value.f64_field("total_joule_w"), Some(0.0));
        assert!(value.get("fluid_reference").unwrap().str_field("material_state").is_some());
        assert!(value.str_field("no_claim").unwrap().contains("verified constant"));
        assert_eq!(transient, run(&args, true));
        assert_eq!(std::fs::read_to_string(&path).unwrap(), scene);
        assert_eq!(std::fs::read(dir.join("fluid.fsmcdpk")).unwrap(), bytes);
        let bad_path = dir.join("bad-boundary.json");
        for bad in [
            scene.replace("\"inflow_temperature_k\":300", "\"inflow_temperature_k\":0"),
            scene.replace("\"backflow_temperature_k\":300", "\"backflow_temperature_k\":-1"),
        ] {
            std::fs::write(&bad_path, bad).unwrap();
            let refused = run(&[bad_path.as_os_str().to_owned()], true);
            assert_eq!(refused.exit_code, crate::exit::INPUT);
            assert!(refused.stdout.is_empty());
        }
    }

    mod transient;
}
