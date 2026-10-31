//! Data-driven bonded concentric speaking-string materials.
//!
//! One or more layer records per course, ordered from centre outward:
//! `layer,key,outer_radius_m,density_kg_m3,young_pa,source,license`.
//! Version 2 also supplies solid bending viscosity and its frequency support;
//! the shared section owner integrates the resulting Kelvin-Voigt bending law.
//! The source and license are explicit input claims, NOT measurements verified
//! by this parser. No preset or representative metal is silently selected.
//! Helically wound strings must use their own geometry/slip model instead.
use std::collections::{BTreeMap, BTreeSet};
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer};
use fs_material::state_point::{DENSITY_PROPERTY, YOUNG_MODULUS_PROPERTY,
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point};
use fs_qty::{Density, DynViscosity, Pressure};

use super::geometry::Course;

pub const HEADER: &str = "frankensim-string-sections-v1";
pub const ROW: &str = "layer,key,outer_radius_m,density_kg_m3,young_pa,source,license";
pub const LOSS_HEADER: &str = "frankensim-string-sections-v2";
pub const LOSS_ROW: &str = "layer,key,outer_radius_m,density_kg_m3,young_pa,bending_viscosity_pa_s,omega_min_rad_s,omega_max_rad_s,source,license";
const BENDING_VISCOSITY_PROPERTY: &str = "solid_bending_viscosity";

#[derive(Clone, Copy)]
struct BendingData {
    viscosity_pa_s: f64,
    omega_band_rad_s: (f64, f64),
}

#[allow(clippy::too_many_arguments)]
fn layer(key: u8, ordinal: usize, radius: f64, density: f64, young: f64,
    bending: Option<BendingData>, source: &str, license: &str)
    -> Result<ConcentricMaterialLayer, String> {
    if source.is_empty() || license.is_empty() || source.contains(',')
        || license.contains(',') || !density.is_finite() || density <= 0.0
        || !young.is_finite() || young <= 0.0 {
        return Err("layer needs positive finite properties, source and license".into());
    }
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        (DENSITY_PROPERTY, Density::DIMS, density),
        (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, young),
    ] {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained(),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: Provenance {
                source: source.into(), license: license.into(), artifact: None,
            },
        }).map_err(|e| format!("invalid property claim: {e}"))?;
    }
    let mut query = QueryPoint::new();
    let mut requirements = vec![
        ScalarPropertyRequirement::try_new(DENSITY_PROPERTY, Density::DIMS,
            ScalarAdmissibility::StrictlyPositive).map_err(|e| e.to_string())?,
        ScalarPropertyRequirement::try_new(YOUNG_MODULUS_PROPERTY, Pressure::DIMS,
            ScalarAdmissibility::StrictlyPositive).map_err(|e| e.to_string())?,
    ];
    if let Some(loss) = bending {
        let (lo, hi) = loss.omega_band_rad_s;
        if !loss.viscosity_pa_s.is_finite() || loss.viscosity_pa_s < 0.0
            || !lo.is_finite() || !hi.is_finite() || lo < 0.0 || hi < lo {
            return Err("bending viscosity needs finite nonnegative eta and ordered nonnegative rad/s support".into());
        }
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(BENDING_VISCOSITY_PROPERTY, DynViscosity::DIMS),
            value: PropertyValue::Scalar { value: loss.viscosity_pa_s, dims: DynViscosity::DIMS },
            validity: ValidityDomain::unconstrained().with("omega", lo, hi),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: Provenance {
                source: source.into(), license: license.into(), artifact: None,
            },
        }).map_err(|e| format!("invalid bending claim: {e}"))?;
        requirements.push(ScalarPropertyRequirement::try_new(
            BENDING_VISCOSITY_PROPERTY, DynViscosity::DIMS, ScalarAdmissibility::Finite,
        ).map_err(|e| e.to_string())?);
        query = query.with("omega", lo + 0.5 * (hi - lo)).map_err(|e| e.to_string())?;
    }
    let card = MaterialCard::assemble(MaterialStateId {
        chemistry: format!("declared-string-{key}-layer-{ordinal}"),
        phase: "solid".into(), process: "authored-source-claim".into(), revision: 0,
    }, claims, vec![]).map_err(|e| format!("invalid material card: {e}"))?;
    let material = resolve_material_state_point(&card, &query, &requirements,
        MaterialPropertySelection::SingleClaimOnly).map_err(|e| e.to_string())?;
    Ok(ConcentricMaterialLayer { outer_radius_m: radius, material })
}

/// Replace mass, bending and axial coefficients for every key, all-or-nothing.
/// Does not change tension or the equations in the string/bridge/contact solver.
/// The retained section binding also supplies the stretching file's `section` EA.
/// V2 additionally binds C_b and the intersection of all layer frequency bands;
/// its material loss replaces the authored demonstration bending spectrum.
/// The source model is perfectly bonded concentric layers; winding is NOT bonded.
///
/// # Errors
/// Missing keys, repeated or out-of-order layers, invalid material claims,
/// hollow/non-solid geometry, and unrepresentable derived dynamics are refused.
pub fn apply(text: &str, courses: &[Course]) -> Result<Vec<Course>, String> {
    if courses.is_empty() { return Err("string material override needs a nonempty scale".into()); }
    let keys: BTreeSet<u8> = courses.iter().map(|c| c.midi).collect();
    if keys.len() != courses.len() { return Err("string scale has duplicate keys".into()); }
    let mut rows: BTreeMap<u8, Vec<ConcentricMaterialLayer>> = BTreeMap::new();
    let mut header_seen = false;
    let mut column_seen = false;
    let mut with_loss = false;
    for (index, raw) in text.lines().enumerate() {
        let row = raw.split('#').next().unwrap_or("").trim();
        if row.is_empty() { continue; }
        let error = |why: &str| format!("string sections line {}: {why}", index + 1);
        if !header_seen {
            with_loss = match row {
                HEADER => false,
                LOSS_HEADER => true,
                _ => return Err(error("expected string-sections v1 or v2 header")),
            };
            header_seen = true;
            continue;
        }
        if !column_seen {
            let expected = if with_loss { LOSS_ROW } else { ROW };
            if row != expected {
                return Err(error("expected explicit SI column header for the declared version"));
            }
            column_seen = true;
            continue;
        }
        let fields: Vec<_> = row.split(',').map(str::trim).collect();
        let expected_fields = if with_loss { 10 } else { 7 };
        if fields.len() != expected_fields || fields[0] != "layer" {
            return Err(error("wrong layer field count for the declared version"));
        }
        let key: u8 = fields[1].parse().map_err(|_| error("invalid key"))?;
        if !keys.contains(&key) { return Err(error("key absent from the selected scale")); }
        let radius: f64 = fields[2].parse().map_err(|_| error("invalid outer radius"))?;
        let rho: f64 = fields[3].parse().map_err(|_| error("invalid density"))?;
        let young: f64 = fields[4].parse().map_err(|_| error("invalid Young modulus"))?;
        let bending = if with_loss {
            Some(BendingData {
                viscosity_pa_s: fields[5].parse().map_err(|_| error("invalid solid bending viscosity"))?,
                omega_band_rad_s: (
                    fields[6].parse().map_err(|_| error("invalid angular-frequency lower bound"))?,
                    fields[7].parse().map_err(|_| error("invalid angular-frequency upper bound"))?,
                ),
            })
        } else { None };
        let source_index = if with_loss { 8 } else { 5 };
        let layers = rows.entry(key).or_default();
        if layers.len() >= fs_material::section::MAX_CONCENTRIC_LAYERS {
            return Err(error("too many layers"));
        }
        if !radius.is_finite() || radius <= 0.0 ||
            layers.last().is_some_and(|last| radius <= last.outer_radius_m) {
            return Err(error("solid concentric layers need strictly increasing positive radii"));
        }
        layers.push(layer(key, layers.len(), radius, rho, young, bending,
            fields[source_index], fields[source_index + 1]).map_err(|e| error(&e))?);
    }
    if !header_seen || !column_seen { return Err("missing string section header or columns".into()); }
    courses.iter().map(|course| {
        let layers = rows.remove(&course.midi)
            .ok_or_else(|| format!("missing bonded section for key {}; no fallback", course.midi))?;
        let section = ConcentricElasticSection::try_new(0.0, layers)
            .map_err(|e| format!("key {}: {e}", course.midi))?;
        let changed = if with_loss {
            let law = section.with_kelvin_voigt_bending(BENDING_VISCOSITY_PROPERTY)
                .map_err(|e| format!("key {}: {e}", course.midi))?;
            course.with_bonded_kelvin_voigt_section(&law)
        } else { course.with_bonded_material_section(&section) };
        changed.map_err(|e| format!("key {}: {e}", course.midi))
    }).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::linear::{Bank, BoardMode, string_stretching::{Specification, MAX_TRIALS}};

    fn material_course(density: f64, young: f64) -> Course {
        let template = Course { length_m: 1.0, tension_n: 400.0, unison: 1,
            duplex_length_m: 0.0, detune_cents: 0.0,
            ..super::super::geometry::demonstration_scale().unwrap()[48] };
        apply(&format!("{HEADER}\n{ROW}\nlayer,69,0.001,{density},{young},synthetic,CC0-1.0\n"),
            &[template]).unwrap()[0]
    }

    fn bank(course: Course) -> Bank {
        Bank::new(&[course], &[BoardMode { frequency_hz: 100.0, damping_ratio: 0.0,
            bridge: [0.0; 88], volume: 0.0 }], 48_000, 3_000.0, 2, false).unwrap()
    }

    fn specification(row: &str, course: Course) -> Result<Specification, String> {
        Specification::read(&format!("{}\n{row}\n", super::super::linear::string_stretching::HEADER), &[course])
    }

    fn near(actual: f64, expected: f64) {
        assert!((actual - expected).abs() <= 1e-11 * expected.abs().max(1e-12),
            "{actual} != {expected}");
    }

    // G1/G3: one data substitution drives every section coefficient and the
    // actual nonlinear owner. The independent reference is sinusoidal strain.
    #[test]
    fn section_data_drives_mass_bending_and_nonlinear_tension_together() {
        let mut responses = Vec::new();
        for (rho, young) in [(8000.0, 80e9), (10000.0, 200e9)] {
            let course = material_course(rho, young);
            let area = core::f64::consts::PI * 1e-6;
            let ea = young * area;
            near(course.linear_density_kg_m, rho * area);
            near(course.flexural_rigidity_nm2, young * core::f64::consts::PI * 1e-12 / 4.0);
            near(course.bonded_section.unwrap().axial_rigidity_n(), ea);
            let selection = specification("stretch,69,section,0.2", course).unwrap();
            near(selection.courses[&69].unwrap().axial_rigidity_n, ea);
            let mut mechanics = bank(course);
            mechanics.configure_string_stretching(&[course], &selection).unwrap();
            let amplitude = 1e-3;
            mechanics.q[0] = amplitude * course.modal_mass_kg().sqrt();
            let strain = (core::f64::consts::PI * amplitude / course.length_m).powi(2) / 4.0;
            let response = mechanics.string_stretching_observation(0).unwrap();
            near(response.additional_strain, strain);
            near(response.tension_n - course.tension_n, ea * strain);
            near(response.stretching_energy_j, 0.5 * ea * course.length_m * strain * strain);
            responses.push((course.partial_hz(1, course.tension_n), response));
        }
        assert_ne!(responses[0].0, responses[1].0);
        assert!(responses[1].1.tension_n > responses[0].1.tension_n);
        near(responses[1].1.stretching_energy_j / responses[0].1.stretching_energy_j, 2.5);
    }

    // The same waves row is re-resolved without copying EA. Changes in the
    // actual material wave speed alter retained modes and the solved trajectory.
    #[test]
    fn section_data_reaches_longitudinal_wave_retention_and_motion() {
        let mut counts = Vec::new();
        let mut states = Vec::new();
        for (rho, young) in [(8000.0_f64, 80e9_f64), (10000.0, 200e9)] {
            let course = material_course(rho, young);
            // This bank's board has an all-zero vertical bridge row, so its fixture axial row is zero too.
            let selection = specification("waves,69,section,0.2,4000,4,0,0.05", course).unwrap()
                .with_fixture_axial_rows(&[BoardMode { frequency_hz: 100.0, damping_ratio: 0.0,
                    bridge: [0.0; 88], volume: 0.0 }]);
            let mut mechanics = bank(course);
            mechanics.configure_string_stretching(&[course], &selection).unwrap();
            let f_axial = 0.5 * (young / rho).sqrt() / course.length_m;
            let expected_count = (4000.0 / f_axial).floor() as usize;
            assert_eq!(mechanics.longitudinal_mode_count(), expected_count);
            let scale = course.modal_mass_kg().sqrt();
            mechanics.q[0] = 1e-4 * scale;
            mechanics.q[1] = 0.5e-4 * scale;
            let initial_energy = mechanics.energy();
            let mut loss = 0.0;
            let forces = vec![0.0; mechanics.contact_strings.len()];
            for _ in 0..32 {
                mechanics.begin_string_stretching_step();
                let mut converged = false;
                for _ in 0..MAX_TRIALS {
                    mechanics.predict(); mechanics.finish(&forces);
                    if mechanics.correct_string_stretching_step().unwrap() {
                        converged = true; break;
                    }
                }
                assert!(converged);
                loss += mechanics.last_modal_loss_j;
                mechanics.commit();
                assert!((mechanics.energy() + loss - initial_energy).abs() < 1e-10);
            }
            let (q, v) = mechanics.string_longitudinal_state(0).unwrap();
            assert!(q.iter().chain(v).any(|value| value.abs() > 1e-16));
            counts.push(expected_count);
            states.push((mechanics.q.clone(), v.to_vec()));
        }
        assert_eq!(counts, vec![2, 1]);
        assert_ne!(states[0], states[1]);
    }

    #[test]
    fn missing_duplicate_and_stale_axial_material_bindings_refuse() {
        let ordinary = super::super::geometry::demonstration_scale().unwrap()[48];
        assert!(specification("stretch,69,section,0.2", ordinary).is_err());
        let original = material_course(8000.0, 80e9);
        assert!(specification("stretch,69,250000,0.2", original).is_err());
        assert!(specification("waves,69,250000,0.2,4000,4,0,0.05", original).is_err());
        assert!(Course { linear_density_kg_m: 2.0 * original.linear_density_kg_m, ..original }.validate().is_err());
        assert!(Course { flexural_rigidity_nm2: 2.0 * original.flexural_rigidity_nm2, ..original }.validate().is_err());
        let stale = specification("stretch,69,section,0.2", original).unwrap();
        let replaced = material_course(10000.0, 200e9);
        let mut mechanics = bank(replaced);
        let before = (mechanics.q.clone(), mechanics.v.clone());
        assert!(mechanics.configure_string_stretching(&[replaced], &stale).is_err());
        assert!(!mechanics.has_string_stretching());
        assert_eq!((mechanics.q.clone(), mechanics.v.clone()), before);
        mechanics.configure_string_stretching(&[replaced],
            &specification("stretch,69,section,0.2", replaced).unwrap()).unwrap();
    }

    #[test]
    fn material_change_reaches_string_spectrum_without_tension_retuning() {
        let source = [super::super::geometry::demonstration_scale().unwrap()[48]];
        let row = |rho: f64, young: f64| format!(
            "{HEADER}\n{ROW}\nlayer,69,0.0005,{rho},{young},synthetic,CC0-1.0\n");
        let light = apply(&row(7800.0, 200e9), &source).unwrap();
        let heavy = apply(&row(10500.0, 80e9), &source).unwrap();
        assert_eq!(light[0].tension_n, heavy[0].tension_n);
        assert!(heavy[0].linear_density_kg_m > light[0].linear_density_kg_m);
        assert!(heavy[0].flexural_rigidity_nm2 < light[0].flexural_rigidity_nm2);
        assert!(heavy[0].partial_hz(1, heavy[0].tension_n) <
            light[0].partial_hz(1, light[0].tension_n));
        assert_eq!(light[0].midi, source[0].midi);
        assert!(apply(&row(7800.0, 200e9), &[]).is_err());
        assert!(apply(&row(-1.0, 200e9), &source).is_err());
        assert!(apply(&row(7800.0, f64::NAN), &source).is_err());
        assert!(apply(&row(7800.0, 200e9).replace("synthetic", ""), &source).is_err());
        assert!(apply(&format!("{}layer,69,0.0004,7800,200000000000,synthetic,CC0\n",
            row(7800.0, 200e9)), &source).is_err());
    }
}

#[cfg(test)]
#[path = "string_section_loss_tests.rs"]
mod loss_tests;
