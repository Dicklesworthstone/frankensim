//! Data-only electrical/material change -> heat history -> actual pressure samples.
//! Cold and hot sounds are frozen-state observations of the SAME specimen and
//! excitation. This does not simulate sound while heating or close mechanical
//! work against the thermal account. Demo materials are explicitly synthetic.
use std::io::{Read, Write};
use fs_conduction::transient::deposition::{Boundary, DepositionBudget, StepConfig};
use fs_couple::material_assembly::{AcousticMaterialBindings,
    concentric::{compile_concentric_material_assembly, thermal::heating::
        compile_heated_concentric_material_assembly}};
use fs_couple::string_specimen::StringPrestress;
use fs_couple::pcm_wav::encode_pcm16_wav_interleaved;
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId,
    NormalizedMaterialCardPack, PropertyClaim, PropertyKey, PropertyValue,
    Provenance, QueryPoint, SelectionPolicy, UncertaintyModel};
use fs_material::conductor::{ELECTRICAL_RESISTIVITY_DIMS, ELECTRICAL_RESISTIVITY_PROPERTY, OhmicDrive};
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer, thermal::AxialConstraint};
use fs_material::section_temperature::{ConcentricSectionTemperatureLaw,
    heating::{HeatedConcentricSection, SectionThermalInterface}};
use fs_material::state_point::{DENSITY_PROPERTY, YOUNG_MODULUS_PROPERTY,
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    INVERSE_TEMPERATURE_DIMS, LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
    SPECIFIC_HEAT_CAPACITY_DIMS, SPECIFIC_HEAT_CAPACITY_PROPERTY,
    THERMAL_CONDUCTIVITY_DIMS, THERMAL_CONDUCTIVITY_PROPERTY, resolve_material_state_point};
use fs_qty::{Density, Pressure, QuantitySpec, Temperature};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_scenario::{AcousticAssembly, AmbientGas, Listener, Pluck, PrestressedString};

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
const LENGTH: f64 = 0.65;
const RADII: [f64; 2] = [0.00025, 0.0005];
const SAMPLE_RATE: u32 = 16000;

fn point(t: f64) -> Result<QueryPoint> { Ok(QueryPoint::new().with("T", t)?) }
fn load(path: &str) -> Result<MaterialCard> {
    const MAX_BYTES: u64 = 16 * 1024 * 1024;
    let mut bytes = Vec::new();
    std::fs::File::open(path)?.take(MAX_BYTES + 1).read_to_end(&mut bytes)?;
    if bytes.len() as u64 > MAX_BYTES { return Err("material pack exceeds 16 MiB".into()); }
    Ok(NormalizedMaterialCardPack::from_bytes(&bytes)?.card().clone())
}
fn caloric(card: &MaterialCard) -> Result<ResolvedHeatCapacityCurve> {
    Ok(ResolvedHeatCapacityCurve::try_new(card, HeatCapacitySourceQuery {
        heat_capacity: &PropertyKey::new(SPECIFIC_HEAT_CAPACITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_DIMS),
        density: &PropertyKey::new(DENSITY_PROPERTY, Density::DIMS),
        temperature_axis: "T", lower: &point(300.0)?, upper: &point(400.0)?,
        selection: SelectionPolicy::SingleClaimOnly,
    }, SolidLiquidPhase::Solid, 0.0, 0.001)?)
}
fn section(cards: [&MaterialCard; 2], expansion: bool) -> Result<ConcentricElasticSection> {
    let mut layers = Vec::new();
    for (card, radius) in cards.into_iter().zip(RADII) {
        let k_quantity = card.claims().registered_quantity(THERMAL_CONDUCTIVITY_PROPERTY)
            .ok_or("missing thermal conductivity")?;
        if k_quantity != QuantitySpec::dimensional(THERMAL_CONDUCTIVITY_DIMS)
            && k_quantity != QuantitySpec::semantic(SemanticType::new(
                QuantityKind::ThermalConductivity, ValueForm::Static)) {
            return Err("thermal conductivity must use its static quantity convention".into());
        }
        let mut requirements = Vec::new();
        for (name, dims) in [(DENSITY_PROPERTY, Density::DIMS), (YOUNG_MODULUS_PROPERTY, Pressure::DIMS),
            (SPECIFIC_HEAT_CAPACITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_DIMS),
            (ELECTRICAL_RESISTIVITY_PROPERTY, ELECTRICAL_RESISTIVITY_DIMS)] {
            requirements.push(ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive)?);
        }
        requirements.push(ScalarPropertyRequirement::try_with_quantity(THERMAL_CONDUCTIVITY_PROPERTY,
            k_quantity, ScalarAdmissibility::StrictlyPositive)?);
        if expansion {
            requirements.push(ScalarPropertyRequirement::try_new(
                LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY, INVERSE_TEMPERATURE_DIMS,
                ScalarAdmissibility::Finite)?);
        }
        layers.push(ConcentricMaterialLayer { outer_radius_m: radius,
            material: resolve_material_state_point(card, &point(300.0)?, &requirements,
                MaterialPropertySelection::SingleClaimOnly)? });
    }
    Ok(ConcentricElasticSection::try_new(0.0, layers)?)
}
fn template() -> AcousticAssembly {
    AcousticAssembly { ambient: AmbientGas::sea_level(), string: Some(PrestressedString {
        length_m: LENGTH, tension_n: 10.0, lin_density_kg_m: 1.0, axial_stiffness_n: 0.0,
        width_m: 0.01, n_modes: 4, damping_ratio: 0.0, rayleigh: None,
        bending_stiffness_n_m2: 0.0, kelvin_voigt_bending: None, relaxation_bending: None,
        polarization_detune: 0.0, moving_end: false, tension_loss_factor: None, bending_loss_factor: None }), duct: None,
        pluck: Some(Pluck { station_frac: 0.25, height_m: 1e-5 }), bow: None, blow: None,
        reed: None, soundboard: None, body_modes: vec![], plate: None, cavity: None,
        obstacles: vec![], contact_texture: None, listener: Listener { distance_m: 1.0 },
        sample_rate_hz: SAMPLE_RATE, duration_s: 0.05 }
}
struct SoundComparison {
    cold_pa: Vec<f64>, hot_pa: Vec<f64>,
    temperatures_k: [f64; 2],
    cold_axial_n: f64, hot_axial_n: f64,
    cold_tension_n: f64, hot_tension_n: f64,
    // Only present when thermal expansion and the actual supports were resolved.
    hot_layer_stresses_pa: Option<Vec<f64>>,
    electrical_energy_j: f64, stored_heat_j: f64, thermal_residual_j: f64,
}
fn synthesize(cards: [&MaterialCard; 2], voltage: f64, duration: f64,
    support: Option<FixedSupport>) -> Result<SoundComparison> {
    if !voltage.is_finite() || !duration.is_finite() || duration <= 0.0 || duration > 10.0 {
        return Err("voltage must be finite and heating duration must lie in (0,10] seconds".into());
    }
    let charts = [caloric(cards[0])?, caloric(cards[1])?];
    let reference = section(cards, support.is_some())?;
    let law = ConcentricSectionTemperatureLaw::try_new(reference.clone(), &cards, &["T", "T"])?;
    // Reduced annular network: uniform layer temperatures are represented at
    // their area-midpoint radii. Logarithmic half paths are explicit geometry,
    // not a claim of a resolved radial PDE or a measured thermal interface.
    let two_pi_l = 2.0 * std::f64::consts::PI * LENGTH;
    let core_node = RADII[0] / 2.0_f64.sqrt();
    let sleeve_node = ((RADII[0] * RADII[0] + RADII[1] * RADII[1]) / 2.0).sqrt();
    let interface = SectionThermalInterface {
        area_m2: two_pi_l * RADII[0],
        resistance_geometry_per_m: [(RADII[0] / core_node).ln() / two_pi_l,
            (sleeve_node / RADII[0]).ln() / two_pi_l],
        contact_resistance_k_w: 0.0,
    };
    let boundary = Boundary { ambient_k: 300.0, htc_w_m2_k: 0.0, reservoir: None, source_w: 0.0 };
    let body = HeatedConcentricSection::try_new(law, &[&charts[0], &charts[1]], LENGTH,
        &[boundary; 2], vec![interface], DepositionBudget {
            max_cells: 2, max_links: 1, max_events: 3, max_substeps: 200000,
        })?;
    let body = match support {
        Some(support) => body.with_axial_constraint(
            AxialConstraint::FixedStrain(support.strain), support.limit)?,
        None => body,
    };
    let (heated, report) = body.advance_to(duration, OhmicDrive::Voltage(voltage),
        StepConfig { max_dt_s: 0.001, energy_tolerance_j: 1e-7 }, || Ok(()))?;
    let template = template();
    let bindings = AcousticMaterialBindings::default();
    // Keep the actual support prescription, not a target pitch. The original
    // mode retains 10 N; fixed strain uses the accepted thermal support force.
    let (cold, hot, cold_tension_n, hot_tension_n, hot_layer_stresses_pa) = if support.is_some() {
        let cold = compile_heated_concentric_material_assembly(&template, &body, &bindings)?;
        let hot = compile_heated_concentric_material_assembly(&template, &heated, &bindings)?;
        let response = heated.axial_state().expect("support was bound").response();
        (cold.realize()?, hot.realize()?,
            body.axial_state().expect("support was bound").response().force_n(), response.force_n(),
            Some(response.layers().iter().map(|layer| layer.stress_pa()).collect()))
    } else {
        (compile_concentric_material_assembly(&template, reference.clone(),
            StringPrestress::FixedTension(10.0), &bindings)?.realize()?,
        compile_concentric_material_assembly(&template, heated.section().clone(),
            StringPrestress::FixedTension(10.0), &bindings)?.realize()?, 10.0, 10.0, None)
    };
    if cold.gas != hot.gas || cold.pressure_pa.len() != hot.pressure_pa.len()
        || cold.pressure_pa.iter().chain(&hot.pressure_pa).any(|p| !p.is_finite()) {
        return Err("cold/hot acoustic observations must share the medium and finite sample layout".into());
    }
    Ok(SoundComparison { cold_pa: cold.pressure_pa, hot_pa: hot.pressure_pa,
        cold_tension_n, hot_tension_n, hot_layer_stresses_pa,
        temperatures_k: [heated.thermal().temperature_k(fs_conduction::transient::deposition::CellId(0))?,
            heated.thermal().temperature_k(fs_conduction::transient::deposition::CellId(1))?],
        cold_axial_n: reference.axial_stiffness_n(), hot_axial_n: heated.section().axial_stiffness_n(),
        electrical_energy_j: report.electrical_energy_j,
        stored_heat_j: report.thermal.stored_energy_change_j, thermal_residual_j: report.thermal.residual_j })
}

#[cfg(test)]
fn demo(resistivity_scale: f64) -> Result<MaterialCard> {
    demo_with_expansion(resistivity_scale, None)
}
fn demo_with_expansion(resistivity_scale: f64, alpha: Option<f64>) -> Result<MaterialCard> {
    if !resistivity_scale.is_finite() || resistivity_scale <= 0.0 {
        return Err("synthetic resistivity scale must be finite and positive".into());
    }
    let mut claims = ClaimSet::new();
    let mut properties = vec![(DENSITY_PROPERTY, Density::DIMS, 8000.0, 8000.0),
        (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, 2e11, 1e11),
        (SPECIFIC_HEAT_CAPACITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_DIMS, 400.0, 800.0),
        (THERMAL_CONDUCTIVITY_PROPERTY, THERMAL_CONDUCTIVITY_DIMS, 20.0, 40.0),
        (ELECTRICAL_RESISTIVITY_PROPERTY, ELECTRICAL_RESISTIVITY_DIMS,
            1e-7 * resistivity_scale, 2e-7 * resistivity_scale)];
    if let Some(alpha) = alpha {
        if !alpha.is_finite() { return Err("synthetic expansion coefficient must be finite".into()); }
        properties.push((LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
            INVERSE_TEMPERATURE_DIMS, alpha, 2.0*alpha));
    }
    for (name, dims, a, b) in properties {
        claims.insert_claim(PropertyClaim { key: PropertyKey::new(name, dims),
            value: PropertyValue::Curve { abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
                knots: vec![(300.0, a), (400.0, b)], dims },
            validity: ValidityDomain::unconstrained().with("T", 300.0, 400.0),
            uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::LinearInside,
            observations: vec![], provenance: Provenance {
                source: "synthetic heated-section demonstration, not measured material data".into(),
                license: "CC0-1.0".into(), artifact: None },
        })?;
    }
    Ok(MaterialCard::assemble(MaterialStateId { chemistry: "synthetic".into(), phase: "solid".into(),
        process: "illustrative data only".into(), revision: 0 }, claims, vec![])?)
}
#[derive(Clone, Copy, Debug, PartialEq)]
struct FixedSupport { strain: f64, limit: f64 }
impl FixedSupport {
    fn try_new(strain: f64, limit: f64) -> Result<Self> {
        if !strain.is_finite() || !limit.is_finite() || strain <= 0.0 || limit >= 1.0 || strain > limit {
            return Err("fixed-end sound requires 0 < STRAIN <= STRAIN_LIMIT < 1".into());
        }
        Ok(Self { strain, limit })
    }
}

// An exact optional suffix; malformed support inputs refuse before loading any
// source pack or creating output. Positional arguments retain their old meaning.
fn split_support(args: &[String]) -> Result<(&[String], Option<FixedSupport>)> {
    let count = match args.first().map(String::as_str) {
        Some("--demo") => 6, Some("--packs") => 7,
        _ => return Err("choose --demo or --packs, optionally followed by --fixed-strain STRAIN STRAIN_LIMIT".into()),
    };
    if args.len() == count { return Ok((args, None)); }
    if args.len() != count + 3 || args[count] != "--fixed-strain" {
        return Err("support option must be the suffix --fixed-strain STRAIN STRAIN_LIMIT".into());
    }
    let support = FixedSupport::try_new(args[count+1].parse()?, args[count+2].parse()?)?;
    Ok((&args[..count], Some(support)))
}

fn run(args: &[String]) -> Result<()> {
    let (args, support) = split_support(args)?;
    let (cards, voltage, duration, full_scale, output) = match args.first().map(String::as_str) {
        Some("--demo") if args.len() == 6 => ([
            demo_with_expansion(1.0, support.map(|_| 1e-5))?,
            demo_with_expansion(args[1].parse()?, support.map(|_| 2e-5))?],
            args[2].parse()?, args[3].parse()?, args[4].parse::<f64>()?, &args[5]),
        Some("--packs") if args.len() == 7 => ([load(&args[1])?, load(&args[2])?],
            args[3].parse()?, args[4].parse()?, args[5].parse::<f64>()?, &args[6]),
        _ => return Err("usage: heated_section_sound --demo SLEEVE_RESISTIVITY_SCALE VOLTAGE_V HEAT_SECONDS FULL_SCALE_PA OUTPUT_PREFIX | --packs CORE.pack SLEEVE.pack VOLTAGE_V HEAT_SECONDS FULL_SCALE_PA OUTPUT_PREFIX".into()),
    };
    if !full_scale.is_finite() || full_scale <= 0.0 { return Err("PCM full scale must be positive finite pascals".into()); }
    let result = synthesize([&cards[0], &cards[1]], voltage, duration, support)?;
    let stereo: Vec<f64> = result.cold_pa.iter().zip(&result.hot_pa)
        .flat_map(|(&cold, &hot)| [cold, hot]).collect();
    let (wav, clips) = encode_pcm16_wav_interleaved(&stereo, SAMPLE_RATE, 2, full_scale)?;
    // No artifact is created unless both thermal and acoustic owners succeed.
    // Existing files are never overwritten. Samples retain physical pascals.
    let file = std::fs::OpenOptions::new().write(true).create_new(true).open(format!("{output}.csv"))?;
    let mut writer = std::io::BufWriter::new(file);
    writeln!(writer, "time_s,cold_pressure_pa,heated_pressure_pa")?;
    for (i, (cold, hot)) in result.cold_pa.iter().zip(&result.hot_pa).enumerate() {
        writeln!(writer, "{},{cold},{hot}", i as f64 / f64::from(SAMPLE_RATE))?;
    }
    writer.flush()?;
    let mut audio = std::fs::OpenOptions::new().write(true).create_new(true).open(format!("{output}.wav"))?;
    audio.write_all(&wav)?;
    eprintln!("stereo comparison: left=cold, right=heated; full_scale_Pa={full_scale}; clips={clips}; no peak normalization");
    eprintln!("support={support:?}; tension_cold_N={}; tension_heated_N={}; heated_layer_stresses_Pa={:?}",
        result.cold_tension_n, result.hot_tension_n, result.hot_layer_stresses_pa);
    eprintln!("T_layers_K={:?}; EA_cold_N={}; EA_heated_N={}; electrical_J={}; stored_heat_J={}; thermal_residual_J={}",
        result.temperatures_k, result.cold_axial_n, result.hot_axial_n,
        result.electrical_energy_j, result.stored_heat_j, result.thermal_residual_j);
    Ok(())
}
fn main() -> Result<()> { run(&std::env::args().skip(1).collect::<Vec<_>>()) }

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn actual_pressure_changes_after_material_driven_heating() {
        let a = demo(1.0).unwrap();
        let b = demo(2.0).unwrap();
        let heated = synthesize([&a, &b], 0.2, 0.1, None).unwrap();
        assert!(heated.temperatures_k.iter().all(|t| *t > 300.0));
        assert!(heated.hot_axial_n < heated.cold_axial_n);
        assert_eq!(heated.cold_pa.len(), 800);
        assert!(heated.cold_pa.iter().any(|p| *p != 0.0));
        assert_ne!(heated.cold_pa, heated.hot_pa);
        assert!((heated.electrical_energy_j - heated.stored_heat_j).abs() < 1e-7);
        let idle = synthesize([&a, &b], 0.0, 0.1, None).unwrap();
        assert_eq!(idle.temperatures_k, [300.0; 2]);
        assert_eq!(idle.cold_pa, idle.hot_pa);
    }
    #[test]
    fn changing_only_sleeve_data_reaches_the_emitted_pressure() {
        let a = demo(1.0).unwrap();
        let b = demo(2.0).unwrap();
        let x = synthesize([&a, &a], 0.2, 0.1, None).unwrap();
        let y = synthesize([&a, &b], 0.2, 0.1, None).unwrap();
        assert_eq!(x.cold_pa, y.cold_pa);
        assert_ne!(x.temperatures_k, y.temperatures_k);
        assert_ne!(x.hot_pa, y.hot_pa);
    }

    #[test]
    fn fixed_support_suffix_is_explicit_and_bounded_for_both_source_modes() {
        for base in ["--demo 2 0.2 1 0.01 out", "--packs a.pack b.pack 0.2 1 0.01 out"] {
            let args: Vec<String> = base.split_whitespace().map(str::to_owned).collect();
            assert_eq!(split_support(&args).unwrap().1, None);
            for (suffix, accepted) in [("--fixed-strain 0.0001 0.01", true),
                ("--fixed-strain 0.01 0.01", true), ("--fixed-strain 0 0.01", false),
                ("--fixed-strain NaN 0.01", false), ("--fixed-strain 0.01 inf", false),
                ("--fixed-strain 0.02 0.01", false), ("--fixed-strain 0.1 1", false),
                ("--fixed-strain 0.01", false), ("--unknown 0.0001 0.01", false)] {
                let mut complete = args.clone();
                complete.extend(suffix.split_whitespace().map(str::to_owned));
                assert_eq!(split_support(&complete).is_ok(), accepted, "{suffix}");
            }
        }
        assert!(split_support(&[]).is_err());
    }

    #[test]
    fn sourced_expansion_reaches_sound_and_missing_data_or_slack_refuse() {
        let plain = demo(1.0).unwrap();
        let a = demo_with_expansion(1.0, Some(1e-5)).unwrap();
        let b = demo_with_expansion(1.0, Some(2e-5)).unwrap();
        let support = Some(FixedSupport::try_new(1e-4, 0.01).unwrap());
        assert!(synthesize([&plain, &a], 0.2, 0.1, support).is_err());
        let x = synthesize([&a, &a], 0.2, 0.1, support).unwrap();
        let y = synthesize([&a, &b], 0.2, 0.1, support).unwrap();
        assert_eq!(x.cold_pa, y.cold_pa);
        assert_eq!(x.temperatures_k, y.temperatures_k);
        assert_eq!(x.electrical_energy_j, y.electrical_energy_j);
        assert!(y.hot_tension_n < x.hot_tension_n && x.hot_tension_n < x.cold_tension_n);
        assert_ne!(x.hot_pa, y.hot_pa);
        assert_eq!(y.hot_layer_stresses_pa.as_ref().unwrap().len(), 2);
        let idle = synthesize([&a, &b], 0.0, 0.1, support).unwrap();
        assert_eq!(idle.hot_tension_n, idle.cold_tension_n);
        assert_eq!(idle.cold_pa, idle.hot_pa);
        let slack = Some(FixedSupport::try_new(1e-9, 0.01).unwrap());
        assert!(synthesize([&a, &b], 0.2, 0.1, slack).is_err());
    }
}
