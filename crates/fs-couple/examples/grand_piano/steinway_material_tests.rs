//! The board's materials are DATA. Pack values reach the board rows exactly; a
//! pack that lacks a required property refuses by name; scaling every pack's
//! E_L by s and density by r moves each mode by sqrt(s/r) and each bridge
//! amplitude by 1/sqrt(r), as K -> sK, M -> rM require; another spruce panel
//! moves the modes inside the Courant-Fischer bounds of its plate law; and the
//! environment record sets the wood's moisture content, so drier air raises
//! the modes as the FPL moisture laws predict.
use super::*;
use fs_matdb::{ClaimSet, PropertyValue};

/// The FPL tables' own state, exactly 12 % moisture content: the moisture
/// laws are the identity there, so the packs' 12 % values reach the board.
const AT_12_PERCENT: BoardEnvironment = BoardEnvironment {
    temperature_k: 293.15,
    relative_humidity: 0.65,
    moisture_content: 0.12,
};

fn humid(temperature_k: f64, relative_humidity: f64) -> BoardEnvironment {
    BoardEnvironment::from_gas(
        &GasState::try_new_moist_air(temperature_k, 101_325.0, relative_humidity).unwrap(),
    )
    .unwrap()
}

const ENGELMANN_SPRUCE: &[u8] = include_bytes!("packs/spruce-engelmann-fpl-gtr282.fsmatpk");
const CFRP_IM7_8552: &[u8] =
    include_bytes!("packs/cfrp-im7-8552-tape-ncamp-cam-rp-2009-015.fsmatpk");

fn pack(bytes: &[u8]) -> NormalizedPack {
    load_pack(bytes).unwrap()
}

fn scalar(pack: &NormalizedPack, name: &str) -> f64 {
    let claims = pack.claims().claims_for(name);
    assert_eq!(claims.len(), 1, "{} {name}", pack.pack_id());
    match &claims[0].1.value {
        PropertyValue::Scalar { value, .. } => *value,
        other => panic!("{} {name}: {other:?}", pack.pack_id()),
    }
}

/// A copy of `source` under `id`, with its observations, and with every scalar
/// claim passed through `edit`; `None` drops the claim.
fn edited(
    source: &NormalizedPack,
    id: &str,
    edit: impl Fn(&str, f64) -> Option<f64>,
) -> NormalizedPack {
    let mut claims = ClaimSet::new();
    for id in source.claims().observation_ids() {
        let dataset = source.claims().observation(id).unwrap().clone();
        assert_eq!(claims.register_observation(dataset).unwrap(), id);
    }
    for (_, claim) in source.claims().claims_ordered() {
        let mut claim = claim.clone();
        if let PropertyValue::Scalar { value, .. } = &mut claim.value {
            match edit(claim.key.name(), *value) {
                Some(edited) => *value = edited,
                None => continue,
            }
        }
        claims.insert_claim(claim).unwrap();
    }
    NormalizedPack::new(
        id,
        "steinway-material-test-v1",
        source.source_artifact(),
        source.redistribution_terms(),
        claims,
        vec![],
        vec![],
    )
    .unwrap()
}

fn close(a: f64, b: f64, relative: f64) -> bool {
    (a - b).abs() <= relative * a.abs().max(b.abs())
}

#[test]
fn board_materials_are_the_source_species_pack_values() {
    let m = model_d_materials(&AT_12_PERCENT).unwrap();
    let (sitka, pine, maple) = (pack(SITKA_SPRUCE), pack(SUGAR_PINE), pack(HARD_MAPLE));
    let e_l = scalar(&sitka, "young_modulus_longitudinal");
    assert!(close(m.panel.e_l, e_l, 1e-15));
    assert!(close(
        m.panel.e_r,
        scalar(&sitka, "er_over_el") * e_l,
        1e-15
    ));
    assert!(close(
        m.panel.g_lr,
        scalar(&sitka, "glr_over_el") * e_l,
        1e-15
    ));
    assert!(close(m.panel.nu_lr, scalar(&sitka, "nu_lr"), 1e-15));
    assert!(close(m.panel.density, scalar(&sitka, "density"), 1e-15));
    assert_eq!(m.panel.pack.pack_id, "spruce-sitka-fpl-gtr282");
    for (beam, source) in [
        (&m.ribs, &pine),
        (&m.cutoff_bar, &pine),
        (&m.bridges, &maple),
    ] {
        let e_l = scalar(source, "young_modulus_longitudinal");
        let g = (scalar(source, "glr_over_el") * e_l * scalar(source, "glt_over_el") * e_l).sqrt();
        assert!(close(beam.e_l, e_l, 1e-15));
        assert!(close(beam.g, g, 1e-15));
        assert!(close(beam.density, scalar(source, "density"), 1e-15));
        assert_eq!(beam.pack.pack_id, source.pack_id());
        assert_eq!(beam.pack.content_hash, source.content_hash().to_hex());
    }
    // FPL-GTR-282 Table 5-3a as the seeds transcribe it: E_L is the bending
    // MOE +10 % (shear deflection removed), density 1120 kg/m3 x SG at 12 % MC.
    for (material_e_l, density, moe_gpa, sg) in [
        (m.panel.e_l, m.panel.density, 10.8, 0.40),
        (m.ribs.e_l, m.ribs.density, 8.2, 0.36),
        (m.bridges.e_l, m.bridges.density, 12.6, 0.63),
    ] {
        assert!(
            close(material_e_l, 1.1 * moe_gpa * 1e9, 1e-12),
            "{material_e_l}"
        );
        assert!(close(density, 1120.0 * sg, 1e-12), "{density}");
    }

    // The board rows carry exactly these values, and the source row names
    // every pack by id and content hash.
    let board = build_with_materials(8, true, &m).unwrap();
    let panel = [
        m.panel.density,
        m.panel.e_l,
        m.panel.e_r,
        m.panel.nu_lr,
        m.panel.g_lr,
    ];
    let triangles: Vec<&str> = board
        .geometry
        .lines()
        .filter(|row| row.starts_with("triangle,"))
        .collect();
    assert!(!triangles.is_empty());
    for row in triangles {
        let fields: Vec<f64> = row
            .split(',')
            .skip(6)
            .take(5)
            .map(|x| x.parse().unwrap())
            .collect();
        assert_eq!(fields, panel, "{row}");
    }
    let beams = [&m.ribs, &m.bridges];
    let mut used = [false; 2];
    for row in board
        .geometry
        .lines()
        .filter(|row| row.starts_with("stiffener,"))
    {
        let fields: Vec<f64> = row
            .split(',')
            .skip(1)
            .take(7)
            .map(|x| x.parse().unwrap())
            .collect();
        let found = beams
            .iter()
            .position(|b| [b.e_l, b.g, b.density] == [fields[0], fields[1], fields[6]]);
        used[found.unwrap_or_else(|| panic!("stiffener material from no pack: {row}"))] = true;
    }
    assert_eq!(used, [true, true]);
    let source = board
        .geometry
        .lines()
        .find(|row| row.starts_with("source,"))
        .unwrap();
    for identity in [
        &m.panel.pack,
        &m.ribs.pack,
        &m.bridges.pack,
        &m.cutoff_bar.pack,
    ] {
        assert!(
            source.contains(&format!("{} {}", identity.pack_id, identity.content_hash)),
            "{source}"
        );
    }
}

#[test]
fn a_missing_or_incomplete_pack_refuses_and_never_falls_back() {
    assert!(load_pack(&[]).is_err());
    assert!(load_pack(&SITKA_SPRUCE[..SITKA_SPRUCE.len() / 2]).is_err());
    let sitka = pack(SITKA_SPRUCE);
    for missing in [
        "young_modulus_longitudinal",
        "density",
        "er_over_el",
        "glr_over_el",
        "nu_lr",
        "moisture_law_green_bending_modulus",
        "volumetric_shrinkage_green_to_ovendry",
    ] {
        let partial = edited(&sitka, "sitka-without-property", |key, value| {
            (key != missing).then_some(value)
        });
        let error = PanelMaterial::from_pack(&partial, &AT_12_PERCENT).unwrap_err();
        assert!(
            error.contains("sitka-without-property") && error.contains(missing),
            "{missing}: {error}"
        );
    }
    // The beam's torsion modulus also needs the L-T plane.
    let pine = edited(&pack(SUGAR_PINE), "pine-without-glt", |key, value| {
        (key != "glt_over_el").then_some(value)
    });
    let error = BeamMaterial::from_pack(&pine, &AT_12_PERCENT).unwrap_err();
    assert!(
        error.contains("pine-without-glt") && error.contains("glt_over_el"),
        "{error}"
    );
}

const KEYS: [u8; 88] = {
    let mut keys = [0u8; 88];
    let mut i = 0;
    while i < 88 {
        keys[i] = 21 + i as u8;
        i += 1;
    }
    keys
};

fn modes(materials: &BoardMaterials, upper_hz: f64) -> Vec<super::super::linear::BoardMode> {
    let preset = build_with_materials(8, true, materials).unwrap();
    super::super::board_geometry::BoardGeometry::read(&preset.geometry)
        .unwrap()
        .prepare(&KEYS, upper_hz)
        .unwrap()
        .modes
}

fn materials(
    panel: &NormalizedPack,
    ribs: &NormalizedPack,
    bridges: &NormalizedPack,
    environment: &BoardEnvironment,
) -> BoardMaterials {
    let ribs = BeamMaterial::from_pack(ribs, environment).unwrap();
    BoardMaterials {
        environment: *environment,
        panel: PanelMaterial::from_pack(panel, environment).unwrap(),
        cutoff_bar: ribs.clone(),
        ribs,
        bridges: BeamMaterial::from_pack(bridges, environment).unwrap(),
    }
}

/// Every elastic constant is a ratio times E_L, so E_L x s scales the whole
/// stiffness by s; density x r scales the whole inertia by r. Then each
/// eigenvalue scales by s/r and each mass-normalized shape by 1/sqrt(r): the
/// drive-point mobility Y(sqrt(s/r) w) = Y(w)/sqrt(s r), the infinite-plate
/// 1/(8 sqrt(D rho h)) law, exactly.
#[test]
fn scaling_every_pack_moves_modes_and_mobility_exactly() {
    let (s, r): (f64, f64) = (1.3, 0.8);
    let k = (s / r).sqrt();
    let scale = |bytes: &[u8], id: &str| {
        edited(&pack(bytes), id, |key, value| {
            Some(match key {
                "young_modulus_longitudinal" => value * s,
                "density" => value * r,
                _ => value,
            })
        })
    };
    let base = reference_materials().unwrap();
    let scaled = materials(
        &scale(SITKA_SPRUCE, "sitka-scaled"),
        &scale(SUGAR_PINE, "pine-scaled"),
        &scale(HARD_MAPLE, "maple-scaled"),
        &base.environment,
    );
    assert!(close(scaled.panel.e_r, s * base.panel.e_r, 1e-15));
    assert!(close(scaled.panel.density, r * base.panel.density, 1e-15));
    assert!(close(scaled.bridges.g, s * base.bridges.g, 1e-15));

    let upper_hz = 450.0;
    let before = modes(&base, upper_hz);
    let after = modes(&scaled, upper_hz * k);
    assert!(
        before.len() >= 5,
        "{} modes below {upper_hz} Hz",
        before.len()
    );
    assert_eq!(before.len(), after.len());
    let (mut frequency_error, mut shape_error) = (0.0f64, 0.0f64);
    for (a, b) in before.iter().zip(&after) {
        frequency_error = frequency_error.max((b.frequency_hz / (k * a.frequency_hz) - 1.0).abs());
        let dot: f64 = a.bridge.iter().zip(&b.bridge).map(|(x, y)| x * y).sum();
        let sign = dot.signum();
        let peak = a.bridge.iter().fold(0.0f64, |m, x| m.max(x.abs()));
        for (x, y) in a.bridge.iter().zip(&b.bridge) {
            shape_error = shape_error.max((sign * y * r.sqrt() - x).abs() / peak);
        }
        assert!(
            close(sign * b.volume * r.sqrt(), a.volume, 1e-7),
            "{} {}",
            a.volume,
            b.volume
        );
    }
    eprintln!(
        "scaled packs s={s} r={r}: {} modes to {upper_hz} Hz, max |f'/(k f) - 1| = {frequency_error:.3e}, max bridge-shape error = {shape_error:.3e}",
        before.len()
    );
    // Two independent slices, each converged to the default 1e-10 relative
    // Ritz residual: measured 2.0e-9 in frequency and 1.7e-9 in bridge shape
    // (ts1, 2026-10-23). A term the packs did not scale shows at 1e-2 and up.
    assert!(frequency_error < 1e-8, "{frequency_error}");
    assert!(shape_error < 1e-8, "{shape_error}");
}

/// Generalized eigenvalue range of the swapped panel's plane-stress law Q'
/// against the default Q. Bending is D = h^3/12 Q at every element and grain
/// angle (a congruence), so K'_panel lies between min and max times K_panel.
fn law_ratio_range(from: &PanelMaterial, to: &PanelMaterial) -> (f64, f64) {
    let q = |m: &PanelMaterial| {
        let nu_rl = m.nu_lr * m.e_r / m.e_l;
        let d = 1.0 - m.nu_lr * nu_rl;
        ([m.e_l / d, m.nu_lr * m.e_r / d, m.e_r / d], m.g_lr)
    };
    let ((a, ga), (b, gb)) = (q(from), q(to));
    // det(B - mu A) = 0 over the 2x2 normal block.
    let qa = a[0] * a[2] - a[1] * a[1];
    let qb = -(b[0] * a[2] + b[2] * a[0] - 2.0 * b[1] * a[1]);
    let qc = b[0] * b[2] - b[1] * b[1];
    let root = (qb * qb - 4.0 * qa * qc).sqrt();
    let mu = [
        (-qb - root) / (2.0 * qa),
        (-qb + root) / (2.0 * qa),
        gb / ga,
    ];
    (
        mu.iter().copied().fold(f64::INFINITY, f64::min),
        mu.iter().copied().fold(0.0, f64::max),
    )
}

/// Another spruce panel (FPL Engelmann), ribs and bridges unchanged. With
/// K'_panel in [lo, hi] K_panel and the panel's inertia scaled by its density
/// ratio, Courant-Fischer bounds every sorted eigenvalue ratio. The swap must
/// move the board, and stay inside what its plate law allows.
#[test]
fn another_spruce_panel_moves_modes_inside_its_plate_law_bounds() {
    let base = reference_materials().unwrap();
    let swapped = BoardMaterials {
        panel: PanelMaterial::from_pack(&pack(ENGELMANN_SPRUCE), &base.environment).unwrap(),
        ..base.clone()
    };
    let (lo, hi) = law_ratio_range(&base.panel, &swapped.panel);
    let rho = swapped.panel.density / base.panel.density;
    let lower = (lo.min(1.0) / rho.max(1.0)).sqrt();
    let upper = (hi.max(1.0) / rho.min(1.0)).sqrt();

    let upper_hz = 450.0;
    let before = modes(&base, upper_hz);
    let after = modes(&swapped, upper_hz * upper * (1.0 + 1e-9));
    assert!(
        before.len() >= 5,
        "{} modes below {upper_hz} Hz",
        before.len()
    );
    assert!(
        after.len() >= before.len(),
        "{} < {}",
        after.len(),
        before.len()
    );
    let ratios: Vec<f64> = before
        .iter()
        .zip(&after)
        .map(|(a, b)| b.frequency_hz / a.frequency_hz)
        .collect();
    eprintln!(
        "Engelmann panel: Q ratio range [{lo:.4}, {hi:.4}], density ratio {rho:.4}, frequency bounds [{lower:.4}, {upper:.4}], measured {ratios:.4?}"
    );
    for (i, ratio) in ratios.iter().enumerate() {
        assert!(
            (lower..=upper).contains(ratio),
            "mode {i}: {ratio} outside [{lower}, {upper}]"
        );
    }
    assert!(
        ratios.iter().any(|ratio| (ratio - 1.0).abs() > 1e-3),
        "the swap did not move the board"
    );
}

/// FPL eq. (4-5) through the environment record: the gas owner's moist air at
/// (T, RH) gives back its RH, and the isotherm pack reproduces Table 4-2 to its
/// printed 0.1 % of moisture content. Dry air is outside the isotherm and refuses.
#[test]
fn the_environment_record_sets_the_wood_moisture_content() {
    for (celsius, rh, table_percent) in [
        (-1.1, 0.40, 7.9),
        (-1.1, 0.65, 12.4),
        (4.4, 0.65, 12.3),
        (15.6, 0.65, 12.1),
        (21.1, 0.40, 7.7),
        (21.1, 0.65, 12.0),
        (21.1, 0.80, 16.0),
    ] {
        let environment = humid(273.15 + celsius, rh);
        assert!(
            (environment.relative_humidity - rh).abs() < 1e-12,
            "{environment:?}"
        );
        assert!(
            (100.0 * environment.moisture_content - table_percent).abs() <= 0.05 + 1e-9,
            "{celsius} C, RH {rh}: {} % against Table 4-2's {table_percent} %",
            100.0 * environment.moisture_content
        );
    }
    let dry =
        BoardEnvironment::from_gas(&GasState::try_new_moist_air(293.15, 101_325.0, 0.0).unwrap());
    assert!(dry.unwrap_err().contains("wood-emc-hailwood-horrobin-fpl-gtr282"));
    let reference = BoardEnvironment::reference().unwrap();
    assert_eq!(reference.temperature_k, 293.15);
    assert!((reference.relative_humidity - 0.65).abs() < 1e-12);
    // At exactly 12 % the moisture laws are the identity, bit for bit.
    let tabled = model_d_materials(&AT_12_PERCENT).unwrap();
    assert_eq!(tabled.panel.density, scalar(&pack(SITKA_SPRUCE), "density"));
    assert!(reference_materials().unwrap().panel != tabled.panel);
}

/// RH 65 -> 40 % at 20 C: the wood dries from about 12.0 to 7.7 % MC, its
/// moduli rise by eq. (5-3) and its density falls, so every mode rises. Each
/// species scales K and M by its own factors, so Courant-Fischer bounds every
/// sorted eigenvalue ratio by the extreme stiffness and mass factors.
#[test]
fn drier_air_raises_the_board_modes_as_the_moisture_laws_predict() {
    let (humid65, dry40) = (humid(293.15, 0.65), humid(293.15, 0.40));
    let (wet, dry) = (
        model_d_materials(&humid65).unwrap(),
        model_d_materials(&dry40).unwrap(),
    );
    let stiffness = [
        dry.panel.e_l / wet.panel.e_l,
        dry.ribs.e_l / wet.ribs.e_l,
        dry.bridges.e_l / wet.bridges.e_l,
    ];
    let mass = [
        dry.panel.density / wet.panel.density,
        dry.ribs.density / wet.ribs.density,
        dry.bridges.density / wet.bridges.density,
    ];
    let min = |v: [f64; 3]| v.iter().copied().fold(f64::INFINITY, f64::min);
    let max = |v: [f64; 3]| v.iter().copied().fold(0.0, f64::max);
    let lower = (min(stiffness) / max(mass)).sqrt();
    let upper = (max(stiffness) / min(mass)).sqrt();
    let upper_hz = 450.0;
    let before = modes(&wet, upper_hz);
    let after = modes(&dry, upper_hz * upper * (1.0 + 1e-9));
    assert!(before.len() >= 5 && after.len() >= before.len());
    let ratios: Vec<f64> = before
        .iter()
        .zip(&after)
        .map(|(a, b)| b.frequency_hz / a.frequency_hz)
        .collect();
    eprintln!(
        "RH 0.65 -> 0.40 at 293.15 K: MC {:.4} -> {:.4}; stiffness factors {stiffness:.4?}, density factors {mass:.4?}; predicted mode ratio bounds [{lower:.4}, {upper:.4}], measured {ratios:.4?}",
        humid65.moisture_content, dry40.moisture_content
    );
    for (i, ratio) in ratios.iter().enumerate() {
        assert!(
            *ratio > 1.0 && (lower..=upper).contains(ratio),
            "mode {i}: {ratio} outside ({lower}, {upper})"
        );
    }
}

/// The shipped environment (20 C, 65 % RH) gives EMC 11.996 %, not the
/// tables' 12.000 %: the default board moves by the measured delta below,
/// reported and never closed by choosing an RH.
#[test]
fn the_default_environment_moves_the_board_by_a_measured_delta() {
    let reference = reference_materials().unwrap();
    let tabled = model_d_materials(&AT_12_PERCENT).unwrap();
    let before = modes(&tabled, 450.0);
    let after = modes(&reference, 450.0 * 1.01);
    let worst = before
        .iter()
        .zip(&after)
        .map(|(a, b)| (b.frequency_hz / a.frequency_hz - 1.0).abs())
        .fold(0.0f64, f64::max);
    eprintln!(
        "default environment {:?}: panel E_L x{:.6}, density x{:.6}; max |mode ratio - 1| to 450 Hz {worst:.3e}",
        reference.environment,
        reference.panel.e_l / tabled.panel.e_l,
        reference.panel.density / tabled.panel.density
    );
    assert!(worst > 0.0 && worst < 1e-3, "{worst}");
}

/// A non-wood panel: NCAMP's Hexcel IM7/8552 unidirectional tape, fiber along
/// the board's grain, ribs and bridges unchanged. Its claims hold over the
/// report's own RTD band, 70 +- 10 F (288.71-299.82 K), which the default
/// 20 C environment lies in; it has no moisture law, so the wood EMC does not
/// touch it. The swap moves every mode inside the Courant-Fischer bounds of its
/// plate law and density. Outside the band the pack refuses by name.
#[test]
fn a_carbon_fibre_panel_moves_modes_inside_its_plate_law_bounds_and_refuses_outside_its_temperature()
 {
    let base = reference_materials().unwrap();
    let cfrp = PanelMaterial::from_pack(&pack(CFRP_IM7_8552), &base.environment).unwrap();
    let msi = 6.894757e9;
    assert!(close(cfrp.e_l, 23.51 * msi, 1e-6), "{}", cfrp.e_l);
    assert!(close(cfrp.e_r, 1.30 * msi, 1e-6), "{}", cfrp.e_r);
    assert!(close(cfrp.g_lr, 0.68 * msi, 1e-6), "{}", cfrp.g_lr);
    assert_eq!((cfrp.nu_lr, cfrp.density), (0.316, 1580.0));
    let swapped = BoardMaterials {
        panel: cfrp,
        ..base.clone()
    };
    let (lo, hi) = law_ratio_range(&base.panel, &swapped.panel);
    let rho = swapped.panel.density / base.panel.density;
    let lower = (lo.min(1.0) / rho.max(1.0)).sqrt();
    let upper = (hi.max(1.0) / rho.min(1.0)).sqrt();
    let upper_hz = 450.0;
    let before = modes(&base, upper_hz);
    let after = modes(&swapped, upper_hz * upper * (1.0 + 1e-9));
    assert!(before.len() >= 5 && after.len() >= before.len());
    let ratios: Vec<f64> = before
        .iter()
        .zip(&after)
        .map(|(a, b)| b.frequency_hz / a.frequency_hz)
        .collect();
    eprintln!(
        "CFRP panel at {} K: Q ratio range [{lo:.4}, {hi:.4}], density ratio {rho:.4}, frequency bounds [{lower:.4}, {upper:.4}], measured {ratios:.4?}",
        base.environment.temperature_k
    );
    for (i, ratio) in ratios.iter().enumerate() {
        assert!(
            (lower..=upper).contains(ratio),
            "mode {i}: {ratio} outside [{lower}, {upper}]"
        );
    }
    assert!(
        ratios.iter().any(|ratio| (ratio - 1.0).abs() > 1e-2),
        "the carbon panel did not move the board"
    );
    // 30 C is outside the RTD band the claims hold for: refused, never extrapolated.
    let error = PanelMaterial::from_pack(&pack(CFRP_IM7_8552), &humid(303.15, 0.65)).unwrap_err();
    eprintln!("CFRP panel at 303.15 K: {error}");
    assert!(
        error.contains("cfrp-im7-8552-tape-ncamp-cam-rp-2009-015"),
        "{error}"
    );
    assert!(
        error.contains("outside the validity of") && error.contains("density [288.706, 299.817] K"),
        "{error}"
    );
}
