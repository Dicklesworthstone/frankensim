//! Source -> store -> consumer for the materials of the fs-couple literature
//! tests (bead `frankensim-irp-v5-instrument-realism-pexfp.32`, DONE-WHEN 2).
//!
//! Each case resolves the compiled packs those tests use, which xtask guards
//! byte for byte against their seeds. It asserts that every resolved property
//! equals its source claim bitwise (the claims are the truth), and drives one
//! real consumer: the string specimen and its stiff-string law, the guitar and
//! casebook plate sections, the free-free bar law, and the reed plate. The
//! constants the tests copied before pexfp.32 are kept here only as a diff
//! oracle for the consumer output: bitwise equal where the copy was exact,
//! otherwise the change is printed with its cause. They are never the asserted
//! truth. (The bead placed these cases in xtask; with guarded compiled packs
//! they live beside the tests, which also avoids an xtask -> fs-plate edge.)

#[path = "support/material_packs.rs"]
mod material_packs;

use fs_couple::string_specimen::{StringPrestress, with_uniform_circular_material_and_prestress};
use fs_matdb::{ClaimSet, MaterialCard, MaterialStateId, PropertyValue, QueryPoint};
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    YOUNG_MODULUS_PROPERTY, resolve_material_state_point,
};
use fs_plate::PlateSection;
use fs_qty::{Density, Dims, Pressure};
use fs_scenario::PrestressedString;
use material_packs::{
    ALUMINUM_2024_T3, CANE_REED_BLANK, CANE_REED_FIXTURE, CANE_REED_HEEL, FPL_MOISTURE_CONTENT,
    MUSIC_WIRE_DENSITY, MUSIC_WIRE_MODULI, MaterialPack, ROSEWOOD_BACK_FIXTURE,
    ROSEWOOD_BAR_FIXTURE, ROSEWOOD_BRAZILIAN, ROSEWOOD_INDIAN, SPRUCE_SITKA, at_moisture_content,
};

/// The values the five tests copied into constants before pexfp.32: a diff
/// oracle for consumer output, never a material source.
mod copied {
    pub const MUSIC_WIRE_E_PA: f64 = 200.0e9;
    pub const MUSIC_WIRE_DENSITY: f64 = 7850.0;
    pub const SPRUCE_E_L: f64 = 11.88e9;
    pub const SPRUCE_DENSITY: f64 = 448.0;
    pub const ROSEWOOD_BRAZILIAN_E_L: f64 = 14.3e9;
    pub const ROSEWOOD_BRAZILIAN_DENSITY: f64 = 800.0;
    pub const ROSEWOOD_INDIAN_E_L: f64 = 13.53e9;
    pub const ROSEWOOD_INDIAN_DENSITY: f64 = 750.0;
    pub const ALUMINUM_E: f64 = 72.4e9;
    pub const ALUMINUM_NU: f64 = 0.33;
    pub const ALUMINUM_DENSITY: f64 = 2705.0;
    pub const CANE_E_HEEL_DRY: f64 = 5.0e9;
    pub const CANE_DENSITY: f64 = 778.9;
    pub const CANE_NU: f64 = 0.30;
}

/// The value of the one scalar claim of `property` [SI].
fn claim_value(pack: &MaterialPack, property: &str) -> f64 {
    let claims = pack.card().claims().claims_for(property);
    let [(_, claim)] = claims.as_slice() else {
        panic!("{}: {} claims of '{property}'", pack.id, claims.len());
    };
    match claim.value {
        PropertyValue::Scalar { value, .. } => value,
        ref other => panic!("{}: '{property}' is not a scalar: {other:?}", pack.id),
    }
}

/// Resolve `property` at its measured state and require it to be its claim, bit for bit.
fn resolved_claim(pack: &MaterialPack, property: &str, dims: Dims) -> f64 {
    let resolved = pack.scalar_as_measured(property, dims).expect("resolves");
    assert_eq!(
        resolved.to_bits(),
        claim_value(pack, property).to_bits(),
        "{} '{property}': resolution must return the claim",
        pack.id
    );
    resolved
}

fn bits(d: [f64; 9]) -> [u64; 9] {
    d.map(f64::to_bits)
}

/// G1/G3: the music-wire density (NBS C447) and Young's modulus (Fuchs 1968)
/// reach the string specimen, and through it the stiff-string law
/// `B = pi^2 E I / (T L^2)`. The bake-off copied E 200 GPa and rho 7850, so B
/// and the fundamental move, by exactly the modulus and density ratios.
#[test]
fn g1_g3_sourced_music_wire_reaches_stiff_string() {
    let (moduli, density) = (
        MaterialPack::seed(MUSIC_WIRE_MODULI),
        MaterialPack::seed(MUSIC_WIRE_DENSITY),
    );
    let e = resolved_claim(&moduli, YOUNG_MODULUS_PROPERTY, Pressure::DIMS);
    let rho = resolved_claim(&density, DENSITY_PROPERTY, Density::DIMS);
    // One card for the specimen: both claims, each keeping its own provenance.
    let mut claims = ClaimSet::new();
    for (pack, property) in [
        (&moduli, YOUNG_MODULUS_PROPERTY),
        (&density, DENSITY_PROPERTY),
    ] {
        let mut claim = pack.card().claims().claims_for(property)[0].1.clone();
        claim.observations.clear();
        claims.insert_claim(claim).expect("claim");
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "music-wire".into(),
            phase: "solid".into(),
            process: "composed: young_modulus from fuchs-1968, density from nbs-c447".into(),
            revision: 0,
        },
        claims,
        Vec::new(),
    )
    .expect("music-wire card");
    let requirement = |name: &str, dims| {
        ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive)
            .expect("requirement")
    };
    let state = resolve_material_state_point(
        &card,
        &QueryPoint::new().with("temperature", 293.0).expect("point"),
        &[
            requirement(DENSITY_PROPERTY, Density::DIMS),
            requirement(YOUNG_MODULUS_PROPERTY, Pressure::DIMS),
        ],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .expect("music-wire state");
    // The bake-off string: d = 1 mm, L = 0.65 m, T = 700 N.
    let (diameter, length, tension) = (1.0e-3, 0.65, 700.0);
    let template = PrestressedString {
        length_m: length,
        tension_n: tension,
        lin_density_kg_m: 1.0,
        axial_stiffness_n: 0.0,
        width_m: diameter,
        n_modes: 8,
        damping_ratio: 0.0,
        rayleigh: None,
        bending_stiffness_n_m2: 0.0,
        kelvin_voigt_bending: None,
        relaxation_bending: None,
        polarization_detune: 0.0,
        moving_end: false,
        tension_loss_factor: None,
        bending_loss_factor: None,
    };
    let string = with_uniform_circular_material_and_prestress(
        template,
        0.5 * diameter,
        &state,
        StringPrestress::FixedTension(tension),
    )
    .expect("music-wire specimen")
    .string();
    let pi = core::f64::consts::PI;
    let (area, inertia) = (pi * diameter * diameter / 4.0, pi * diameter.powi(4) / 64.0);
    // The specimen carries the sourced values.
    assert!(
        (string.bending_stiffness_n_m2 / (e * inertia) - 1.0).abs() < 1e-14,
        "EI"
    );
    assert!(
        (string.lin_density_kg_m / (rho * area) - 1.0).abs() < 1e-14,
        "mu"
    );
    let law = |ei: f64, mu: f64| {
        (
            pi * pi * ei / (tension * length * length),
            (tension / mu).sqrt() / (2.0 * length),
        )
    };
    let (b, f0) = law(string.bending_stiffness_n_m2, string.lin_density_kg_m);
    let (b_copied, f0_copied) = law(
        copied::MUSIC_WIRE_E_PA * inertia,
        copied::MUSIC_WIRE_DENSITY * area,
    );
    // The copies were not the claims, so the law moves, by exactly their ratios.
    assert!(
        ((b / b_copied) / (e / copied::MUSIC_WIRE_E_PA) - 1.0).abs() < 1e-13,
        "B moves with E only"
    );
    assert!(
        ((f0 / f0_copied) / (copied::MUSIC_WIRE_DENSITY / rho).sqrt() - 1.0).abs() < 1e-13,
        "f0 with rho only"
    );
    println!(
        "music wire: B {b:.6e} vs copied {b_copied:.6e} ({:+.3} %, E {e:.0} Pa from {} vs copied 200 GPa); \
         f0 {f0:.4} vs {f0_copied:.4} Hz ({:+.3} %, rho {rho} kg/m3 from {} vs copied 7850)",
        100.0 * (b / b_copied - 1.0),
        moduli.id,
        100.0 * (f0 / f0_copied - 1.0),
        density.id
    );
}

/// G1/G3: the Sitka spruce reference ratios reach the guitar top through the
/// M0.2-W1 adapter, and the rosewood E_L plus the labelled back fixture reach
/// the guitar back. Both copies were exact in the four plane constants the DKT
/// section reads, so the sections are bitwise the copied ones.
#[test]
fn g1_g3_sourced_tonewoods_reach_guitar_plate() {
    let at = at_moisture_content(FPL_MOISTURE_CONTENT);
    let spruce = MaterialPack::seed(SPRUCE_SITKA);
    let e_l = resolved_claim(&spruce, "young_modulus_longitudinal", Pressure::DIMS);
    let er = resolved_claim(&spruce, "er_over_el", Dims::NONE);
    let nu = resolved_claim(&spruce, "nu_lr", Dims::NONE);
    let glr = resolved_claim(&spruce, "glr_over_el", Dims::NONE);
    let rho = resolved_claim(&spruce, "density", Density::DIMS);
    let (plane_rho, law) = spruce.plane_lr(&at).expect("spruce L-R plane");
    assert_eq!(
        [law.e[0], law.e[1], law.nu, law.g, plane_rho].map(f64::to_bits),
        [e_l, e_l * er, nu, e_l * glr, rho].map(f64::to_bits),
        "the adapter derives the plane from the claims"
    );
    let top = PlateSection::orthotropic_plane_stress(
        law.e[0], law.e[1], law.nu, law.g, 0.0029, plane_rho,
    )
    .expect("top");
    let copied_law = fs_material::elastic::OrthotropicElastic::new(
        [
            copied::SPRUCE_E_L,
            0.078 * copied::SPRUCE_E_L,
            0.078 * copied::SPRUCE_E_L,
        ],
        [0.372, 0.372, 0.435],
        [
            0.064 * copied::SPRUCE_E_L,
            0.003 * copied::SPRUCE_E_L,
            0.061 * copied::SPRUCE_E_L,
        ],
        1.0,
    )
    .expect("copied spruce law");
    let copied_top =
        PlateSection::orthotropic(&copied_law, 0.0029, copied::SPRUCE_DENSITY).expect("copied top");
    assert_eq!(
        bits(top.d),
        bits(copied_top.d),
        "spruce top: exact copy, bitwise section"
    );
    assert_eq!(top.density.to_bits(), copied_top.density.to_bits());

    let rosewood = MaterialPack::seed(ROSEWOOD_BRAZILIAN);
    let rw_e_l = rosewood
        .scalar("young_modulus_longitudinal", Pressure::DIMS, &at)
        .expect("rosewood E_L");
    assert_eq!(
        rw_e_l.to_bits(),
        claim_value(&rosewood, "young_modulus_longitudinal").to_bits()
    );
    let back = MaterialPack::fixture(ROSEWOOD_BACK_FIXTURE);
    let value = |p, dims| back.scalar(p, dims, &at).expect("back fixture");
    let back_section = PlateSection::orthotropic_plane_stress(
        rw_e_l,
        value("er_over_el", Dims::NONE) * rw_e_l,
        value("nu_lr", Dims::NONE),
        value("glr_over_el", Dims::NONE) * rw_e_l,
        0.0024,
        value("density", Density::DIMS),
    )
    .expect("back");
    let copied_back_law = fs_material::elastic::OrthotropicElastic::new(
        [
            copied::ROSEWOOD_BRAZILIAN_E_L,
            0.10 * copied::ROSEWOOD_BRAZILIAN_E_L,
            0.10 * copied::ROSEWOOD_BRAZILIAN_E_L,
        ],
        [0.372, 0.372, 0.435],
        [
            0.07 * copied::ROSEWOOD_BRAZILIAN_E_L,
            0.004 * copied::ROSEWOOD_BRAZILIAN_E_L,
            0.065 * copied::ROSEWOOD_BRAZILIAN_E_L,
        ],
        1.0,
    )
    .expect("copied rosewood law");
    let copied_back =
        PlateSection::orthotropic(&copied_back_law, 0.0024, copied::ROSEWOOD_BRAZILIAN_DENSITY)
            .expect("copied back");
    assert_eq!(
        bits(back_section.d),
        bits(copied_back.d),
        "rosewood back: bitwise section"
    );
    assert_eq!(
        back_section.density.to_bits(),
        copied_back.density.to_bits()
    );
    println!(
        "guitar plate: spruce top and rosewood back sections bitwise equal to the copied ones"
    );
}

/// G1/G3: the 2024-T3 property set reaches the casebook's simply supported
/// plate; the copy was exact, so the section is bitwise the copied one.
#[test]
fn g1_g3_sourced_aluminium_reaches_casebook_plate() {
    let pack = MaterialPack::seed(ALUMINUM_2024_T3);
    let e = resolved_claim(&pack, YOUNG_MODULUS_PROPERTY, Pressure::DIMS);
    let nu = resolved_claim(&pack, "poisson_ratio", Dims::NONE);
    let rho = resolved_claim(&pack, DENSITY_PROPERTY, Density::DIMS);
    let plate = PlateSection::isotropic(e, nu, 0.002, rho).expect("plate");
    let copied_plate = PlateSection::isotropic(
        copied::ALUMINUM_E,
        copied::ALUMINUM_NU,
        0.002,
        copied::ALUMINUM_DENSITY,
    )
    .expect("copied plate");
    assert_eq!(
        bits(plate.d),
        bits(copied_plate.d),
        "aluminium: bitwise section"
    );
    assert_eq!(plate.density.to_bits(), copied_plate.density.to_bits());
}

/// G1/G3: the Indian rosewood E_L reaches the free-free Euler-Bernoulli bar
/// `f1 = (beta1 L)^2 / (2 pi L^2) sqrt(E I / (rho A))`; the pack states no
/// 12 % density and refuses one by name, so the density is the labelled bar
/// fixture. Both values were exact copies: f1 is bitwise the copied one.
#[test]
fn g1_g3_sourced_rosewood_reaches_free_free_bar() {
    let at = at_moisture_content(FPL_MOISTURE_CONTENT);
    let pack = MaterialPack::seed(ROSEWOOD_INDIAN);
    let e_l = pack
        .scalar("young_modulus_longitudinal", Pressure::DIMS, &at)
        .expect("E_L");
    assert_eq!(
        e_l.to_bits(),
        claim_value(&pack, "young_modulus_longitudinal").to_bits()
    );
    let refused = pack.scalar("density", Density::DIMS, &at).unwrap_err();
    assert!(refused.contains("'density'"), "{refused}");
    let rho = MaterialPack::fixture(ROSEWOOD_BAR_FIXTURE)
        .scalar("density", Density::DIMS, &at)
        .expect("bar density");
    // The xylophone fixture's bar: 0.38 x 0.038 x 0.018 m; beta1 L of cosh cos = 1.
    let (length, width, thickness, beta_l): (f64, f64, f64, f64) =
        (0.38, 0.038, 0.018, 4.730_040_744_862_704);
    let f1 = |e: f64, rho: f64| {
        let (area, inertia) = (width * thickness, width * thickness.powi(3) / 12.0);
        beta_l * beta_l / (2.0 * core::f64::consts::PI * length * length)
            * (e * inertia / (rho * area)).sqrt()
    };
    assert_eq!(
        f1(e_l, rho).to_bits(),
        f1(copied::ROSEWOOD_INDIAN_E_L, copied::ROSEWOOD_INDIAN_DENSITY).to_bits(),
        "bar f1 bitwise"
    );
    println!(
        "free-free rosewood bar: f1 {:.3} Hz, bitwise equal to the copied one",
        f1(e_l, rho)
    );
}

/// G1/G3: the cane heel E_L reaches the reed plate in both of the pack's
/// states, dry and water-soaked (the played state). Density (the blank's low
/// value) and the fixture Poisson ratio are held at their dry state: no wet
/// density or Poisson ratio is sourced. The dry plate is bitwise the copied
/// one; the wet plate is new, its bending rigidity scaled by the modulus ratio.
#[test]
fn g1_g3_sourced_cane_reaches_reed_plate() {
    let heel = MaterialPack::seed(CANE_REED_HEEL);
    let dry = resolved_claim(&heel, "young_modulus_longitudinal", Pressure::DIMS);
    let wet = resolved_claim(
        &heel,
        "young_modulus_longitudinal_water_soaked",
        Pressure::DIMS,
    );
    let rho = resolved_claim(
        &MaterialPack::seed(CANE_REED_BLANK),
        "density_printed_range_low",
        Density::DIMS,
    );
    let nu = MaterialPack::fixture(CANE_REED_FIXTURE)
        .scalar_as_measured("poisson_ratio", Dims::NONE)
        .expect("fixture Poisson ratio");
    // The vamp's mid thickness, 0.5 mm.
    let plate = |e| PlateSection::isotropic(e, nu, 0.5e-3, rho).expect("reed plate");
    let (dry_plate, wet_plate) = (plate(dry), plate(wet));
    let copied_plate = PlateSection::isotropic(
        copied::CANE_E_HEEL_DRY,
        copied::CANE_NU,
        0.5e-3,
        copied::CANE_DENSITY,
    )
    .expect("copied reed plate");
    assert_eq!(
        bits(dry_plate.d),
        bits(copied_plate.d),
        "dry reed plate bitwise"
    );
    assert!(wet < dry, "water-soaked cane is softer");
    for (w, d) in wet_plate.d.iter().zip(&dry_plate.d) {
        if *d != 0.0 {
            assert!((w / d / (wet / dry) - 1.0).abs() < 1e-14, "D scales with E");
        }
    }
    println!(
        "reed plate: water-soaked (played) E {wet:.4e} Pa vs dry {dry:.4e} Pa; D x{:.4}, plate frequencies x{:.4} \
         (a state the pre-pack test never used)",
        wet / dry,
        (wet / dry).sqrt()
    );
}
