//! Primary dilute-gas references and source-to-acoustics checks. These formulas
//! are test oracles only; production resolves the existing material curves.

use super::common_material_acquisition::{compile, compile_thermochemical_species, sample};
use fs_matdb::{
    MaterialStateId, NormalizedMaterialCardPack, NormalizedModelPack, NormalizedSpeciesPack,
    QueryPoint,
};
use fs_material::{
    gas::{GasSpec, GasState, R_USSA_1976},
    state_point::MaterialPropertySelection,
    thermochemical_gas::{
        ThermochemicalGasModel,
        packs::{SpeciesPackRef, from_packs},
    },
};
use fs_math::det;
use fs_qty::{
    QuantitySpec,
    semantic::{Composition, CompositionBasis, QuantityKind, SemanticType, ValueForm},
};

pub(super) const SPECIES: [(&str, &str, &str); 5] = [
    ("N2", "nitrogen", "nitrogen-transport-lemmon2004"),
    ("O2", "oxygen", "oxygen-transport-lemmon2004"),
    ("Ar", "argon", "argon-transport-lemmon2004"),
    (
        "CO2",
        "carbon-dioxide",
        "carbon-dioxide-transport-laesecke2017-huber2016",
    ),
    ("H2O", "water-vapor", "water-vapor-transport-iapws2008-2011"),
];

fn power(x: f64, exponent: f64) -> f64 {
    det::exp(exponent * det::ln(x))
}
fn inverse_polynomial(coefficients: &[f64], x: f64) -> f64 {
    let inverse = 1.0 / x;
    coefficients
        .iter()
        .rev()
        .fold(0.0, |value, &coefficient| value * inverse + coefficient)
}

/// SI dilute-gas viscosity and conductivity; constants are the source's facts.
/// N2/O2/Ar: Lemmon/Jacobsen2004 Eqs2/5, TablesI/II/IV (PDF pp3,4,6).
/// CO2: Laesecke/Muzny2017 Eq4/Table2 and Huber2016 Eq3/Table3.
/// H2O: IAPWS R12-08 Eq11/Table1 and R15-11 Eq16/Table1, rho->0.
pub(super) fn reference_transport(species: &str, t: f64) -> [f64; 2] {
    match species {
        "N2" | "O2" | "Ar" => {
            let (mass, epsilon, sigma, tc, n1, n2, t2, n3, t3) = match species {
                "N2" => (
                    28.01348, 98.94, 0.3656, 126.192, 1.511, 2.117, -1.0, -3.332, -0.7,
                ),
                "O2" => (
                    31.9988, 118.5, 0.3428, 154.581, 1.036, 6.283, -0.9, -4.262, -0.6,
                ),
                _ => (
                    39.948, 143.2, 0.335, 150.687, 0.8158, -0.4320, -0.77, 0.0, -1.0,
                ),
            };
            let logarithm = det::ln(t / epsilon);
            let b = [0.431, -0.4623, 0.08406, 0.005341, -0.00331];
            let collision = det::exp(
                b.iter()
                    .rev()
                    .fold(0.0, |value, &coefficient| value * logarithm + coefficient),
            );
            let viscosity_micro = 0.0266958 * det::sqrt(mass * t) / (sigma * sigma * collision);
            [
                viscosity_micro * 1e-6,
                (n1 * viscosity_micro + n2 * power(tc / t, t2) + n3 * power(tc / t, t3)) * 1e-3,
            ]
        }
        "CO2" => {
            let a = [
                1749.354893188350,
                -369.069300007128,
                5423856.34887691,
                -2.21283852168356,
                -269503.247933569,
                73145.021531826,
                5.34368649509278,
            ];
            let third = power(t, 1.0 / 3.0);
            let denominator = a[0]
                + a[1] * power(t, 1.0 / 6.0)
                + a[2] * det::exp(a[3] * third)
                + (a[4] + a[5] * third) / det::exp(third)
                + a[6] * det::sqrt(t);
            let viscosity = 1.0055 * det::sqrt(t) / denominator * 1e-3; // source mPa s
            let tr = t / 304.1282;
            let conductivity = det::sqrt(tr)
                / inverse_polynomial(
                    &[0.0151874307, 0.0280674040, 0.0228564190, -0.00741624210],
                    tr,
                )
                * 1e-3;
            [viscosity, conductivity]
        }
        "H2O" => {
            let tr = t / 647.096;
            [
                100.0 * det::sqrt(tr)
                    / inverse_polynomial(&[1.67752, 2.20462, 0.6366564, -0.241605], tr)
                    * 1e-6,
                det::sqrt(tr)
                    / inverse_polynomial(
                        &[
                            0.002443221,
                            0.01323095,
                            0.006770357,
                            -0.003454586,
                            0.0004096266,
                        ],
                        tr,
                    )
                    * 1e-3,
            ]
        }
        _ => panic!("no primary transport reference for {species}"),
    }
}

pub(super) struct Source {
    pub species: NormalizedSpeciesPack,
    pub nasa9: NormalizedModelPack,
    pub transport: NormalizedMaterialCardPack,
}

pub(super) fn sources() -> Vec<Source> {
    SPECIES
        .into_iter()
        .map(|(id, slug, transport_slug)| {
            let (species, nasa9) = compile_thermochemical_species(slug);
            assert_eq!(species.association().species().as_str(), id);
            let (claims, _) = compile(transport_slug);
            let transport = NormalizedMaterialCardPack::new(
                MaterialStateId {
                    chemistry: id.into(),
                    phase: "dilute gas".into(),
                    process: "explicit primary correlation source association".into(),
                    revision: 0,
                },
                claims,
            )
            .unwrap();
            let bytes = transport.to_bytes();
            let transport =
                NormalizedMaterialCardPack::from_bytes_verified(transport.content_hash(), &bytes)
                    .unwrap();
            Source {
                species,
                nasa9,
                transport,
            }
        })
        .collect()
}

pub(super) fn point(t: f64, p: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "temperature",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            t,
        )
        .unwrap()
        .with_quantity(
            "pressure",
            QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static)),
            p,
        )
        .unwrap()
}

/// Explicit caller-authored dry carrier; no source or species is inferred.
pub(super) fn mixture(
    sources: &[Source],
    carrier: [f64; 3],
    water: f64,
    co2: f64,
) -> ThermochemicalGasModel {
    let balance = 1.0 - water - co2;
    let oxygen = carrier[1] * balance;
    let argon = carrier[2] * balance;
    // Nitrogen is the caller-declared balance gas. Form its residual in the
    // owner's canonical Ar/CO2/H2O/N2/O2 summation order to retain exact unity;
    // H2O and CO2 stay exactly the authored fractions (not renormalized).
    let mut nitrogen = (1.0 - oxygen) - ((argon + co2) + water);
    for _ in 0..4 {
        let sum = (((argon + co2) + water) + nitrogen) + oxygen;
        if sum == 1.0 {
            break;
        }
        nitrogen = if sum < 1.0 {
            nitrogen.next_up()
        } else {
            nitrogen.next_down()
        };
    }
    assert_eq!((((argon + co2) + water) + nitrogen) + oxygen, 1.0);
    assert!((nitrogen - carrier[0] * balance).abs() <= 4.0 * f64::EPSILON);
    let fractions = [nitrogen, oxygen, argon, co2, water];
    let mut bindings = Vec::new();
    let mut active = Vec::new();
    for (source, &fraction) in sources.iter().zip(&fractions) {
        if fraction == 0.0 {
            continue;
        }
        bindings.push(SpeciesPackRef {
            species: &source.species,
            nasa9: Some(&source.nasa9),
            transport: Some((
                &source.transport,
                MaterialPropertySelection::SingleClaimOnly,
            )),
        });
        active.push(fraction);
    }
    let water_id = if water > 0.0 {
        Some(sources[4].species.association().species())
    } else {
        None
    };
    from_packs(
        &bindings,
        Composition::new(CompositionBasis::MoleFraction, active).unwrap(),
        water_id,
    )
    .unwrap()
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a / b - 1.0).abs() <= tolerance,
        "{a} vs {b}, relative gate {tolerance}"
    );
}

/// Three measured interpolation rungs plus verification of every retained knot.
/// Published uncertainty is separate from this representation error.
fn transport_ladder(sources: &[Source]) {
    let pins = [
        ("N2", 17.8771e-6, 25.9361e-3),
        ("Ar", 22.7241e-6, 17.8042e-3),
        ("O2", 20.6307e-6, 26.4403e-3),
    ];
    for (species, mu, k) in pins {
        let reference = reference_transport(species, 300.0);
        assert!((reference[0] - mu).abs() <= 0.00005e-6);
        assert!((reference[1] - k).abs() <= 0.00005e-3);
    }
    for source in sources {
        let species = source.species.association().species().as_str();
        let low = if species == "H2O" { 250.0 } else { 200.0 };
        for (column, name) in ["dynamic_viscosity", "thermal_conductivity"]
            .into_iter()
            .enumerate()
        {
            let claim = source.transport.card().claims().claims_for(name)[0].1;
            let fs_matdb::PropertyValue::Curve { knots, .. } = &claim.value else {
                panic!("transport must be an explicit curve")
            };
            for &(t, value) in knots {
                close(value, reference_transport(species, t)[column], 1e-12);
            }
            let mut errors = Vec::new();
            for intervals in [64, 128, 256] {
                let step = (1000.0 - low) / intervals as f64;
                let mut error: f64 = 0.0;
                for cell in 0..intervals {
                    let t0 = low + cell as f64 * step;
                    let t1 = t0 + step;
                    let a = reference_transport(species, t0)[column];
                    let b = reference_transport(species, t1)[column];
                    for fraction in [0.25, 0.5, 0.75] {
                        let t = t0 + fraction * step;
                        error = error.max(
                            ((a * (1.0 - fraction) + b * fraction)
                                / reference_transport(species, t)[column]
                                - 1.0)
                                .abs(),
                        );
                    }
                }
                errors.push(error);
            }
            assert!(
                errors[2] < errors[1] && errors[1] < errors[0],
                "{species}/{name}: {errors:?}"
            );
            let order = det::ln(errors[1] / errors[2]) / det::ln(2.0);
            assert!(
                order >= 1.8,
                "{species}/{name}: interpolation order {order}"
            );
            assert!(
                errors[2] <= 1e-4,
                "{species}/{name}: final interpolation error {}",
                errors[2]
            );
            // The shipped curve uses the final uniform rung plus additional
            // query knots. Check the actual resolver between ALL stored knots.
            for pair in knots.windows(2) {
                let t = (pair[0].0 + pair[1].0) * 0.5;
                let bundle = fs_material::state_point::resolve_material_state_point(
                    source.transport.card(),
                    &point(t, 101325.0),
                    &[
                        fs_material::state_point::ScalarPropertyRequirement::try_new(
                            name,
                            claim.key.dims(),
                            fs_material::state_point::ScalarAdmissibility::StrictlyPositive,
                        )
                        .unwrap(),
                    ],
                    MaterialPropertySelection::SingleClaimOnly,
                )
                .unwrap();
                close(
                    bundle.property(name).unwrap().value_si(),
                    reference_transport(species, t)[column],
                    1e-4,
                );
            }
            println!(
                "{{\"suite\":\"matdb-pack\",\"bead\":\"M0.1-W1\",\"test\":\"transport_interpolation_ladder\",\"species\":\"{species}\",\"property\":\"{name}\",\"rungs\":[64,128,256],\"observed\":{{\"relative_errors\":{:?},\"order\":{order}}},\"tolerance\":0.0001,\"verdict\":\"matched\"}}",
                errors
            );
        }
    }
}

fn lossless_resonance(gas: &GasState, mode: usize) -> f64 {
    let expected = (2 * mode - 1) as f64 * gas.sound_speed / (4.0 * 0.57);
    let admittance_imaginary = |f| {
        let z = impedance(gas, f, fs_duct::LossModel::Lossless);
        -z.im / z.norm_sq()
    };
    let (mut lo, mut hi) = (0.9 * expected, 1.1 * expected);
    let mut left = admittance_imaginary(lo);
    assert!(left.is_sign_positive() != admittance_imaginary(hi).is_sign_positive());
    for _ in 0..60 {
        let mid = 0.5 * (lo + hi);
        if mid == lo || mid == hi {
            break;
        }
        let value = admittance_imaginary(mid);
        if value == 0.0 {
            return mid;
        }
        if value.is_sign_positive() == left.is_sign_positive() {
            lo = mid;
            left = value;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

fn impedance(gas: &GasState, f: f64, loss: fs_duct::LossModel) -> fs_math::c64::C64 {
    fs_duct::input_impedance(
        &fs_duct::Duct {
            segments: vec![fs_duct::Segment::Cylinder {
                radius: 0.0073,
                length: 0.57,
            }],
        },
        gas,
        core::f64::consts::TAU * f,
        loss,
        fs_duct::Termination::IdealOpen,
    )
    .unwrap()
    .impedance
}

fn bessel_peak(gas: &GasState, mode: usize) -> f64 {
    let f0 = (2 * mode - 1) as f64 * gas.sound_speed / (4.0 * 0.57);
    let (mut lo, mut hi) = (0.85 * f0, 1.05 * f0);
    for _ in 0..90 {
        let a = lo + (hi - lo) / 3.0;
        let b = hi - (hi - lo) / 3.0;
        if impedance(gas, a, fs_duct::LossModel::Bessel).norm_sq()
            < impedance(gas, b, fs_duct::LossModel::Bessel).norm_sq()
        {
            lo = a;
        } else {
            hi = b;
        }
    }
    (lo + hi) * 0.5
}

fn cents(a: f64, b: f64) -> f64 {
    1200.0 * det::ln(a / b) / det::ln(2.0)
}

pub(super) fn mixture_reaches_consumers() {
    let started = std::time::Instant::now();
    let sources = sources();
    transport_ladder(&sources);
    // This older information-only pack has legacy dimensional flag axes;
    // compile() above deliberately requires typed v3 for the NEW transport data.
    let manifest = super::workspace_path(
        "data/matdb/seed-v1/nist-srm-1720-northern-continental-air/manifest.tsv",
    );
    let output = super::fixture_dir().join("srm1720.fsmatpk");
    let run = super::run_compiler(&manifest, &output);
    assert!(
        run.status.success(),
        "{}",
        String::from_utf8_lossy(&run.stderr)
    );
    let bytes = std::fs::read(output).unwrap();
    let srm = fs_matdb::NormalizedPack::from_bytes(&bytes).unwrap();
    let srm = fs_matdb::NormalizedPack::from_bytes_verified(srm.content_hash(), &bytes).unwrap();
    let oxygen = sample(&srm, "information_oxygen_amount_fraction", &[]);
    let argon = sample(&srm, "information_argon_amount_fraction", &[]);
    let carrier = [1.0 - oxygen - argon, oxygen, argon];
    let saturation = fs_material::gas::saturation_pressure_water_pa(308.15).unwrap() / 101325.0;
    for (fixture, t, water, co2) in [
        (
            "SRM1720 information components; CO2 absent",
            293.15,
            0.0,
            0.0,
        ),
        (
            "breath; explicit4percent CO2; saturated liquid-water Buck condition",
            308.15,
            saturation,
            0.04,
        ),
        (
            "same SRM1720 information carrier at600K; CO2 absent",
            600.0,
            0.0,
            0.0,
        ),
    ] {
        let model = mixture(&sources, carrier, water, co2);
        let resolved = if water > 0.0 {
            model.evaluate_with_liquid_water_limit(&point(t, 101325.0))
        } else {
            model.evaluate(&point(t, 101325.0))
        }
        .unwrap();
        let gas = resolved.state();
        gas.validate_coherence().unwrap();
        let properties = resolved.thermodynamics().properties();
        let cp = properties.mass_specific_cp().value();
        let gamma = cp / properties.mass_specific_cv().value();
        let r = resolved.equation_of_state().specific_gas_constant().value();
        let rho = resolved.equation_of_state().mass_density().value();
        let c = det::sqrt(gamma * r * t);
        let independent = GasState {
            temperature: t,
            pressure: 101325.0,
            density: rho,
            sound_speed: c,
            dynamic_viscosity: gas.dynamic_viscosity,
            thermal_conductivity: gas.thermal_conductivity,
            gamma,
            specific_gas_constant: r,
            specific_heat_cp: cp,
            prandtl: gas.prandtl,
            characteristic_impedance: rho * c,
            water_mole_fraction: gas.water_mole_fraction,
        };
        assert_eq!(*gas, independent);
        for f in [100.0, 300.0, 900.0] {
            let a = impedance(gas, f, fs_duct::LossModel::Bessel);
            let b = impedance(&independent, f, fs_duct::LossModel::Bessel);
            assert_eq!(
                (a.re.to_bits(), a.im.to_bits()),
                (b.re.to_bits(), b.im.to_bits())
            );
        }
        let bem = fs_bem::helmholtz::Medium::from(gas);
        let rayleigh = fs_couple::pcm_wav::baffled::RayleighMedium::from(gas);
        for (density, speed) in [
            (bem.density, bem.sound_speed),
            (rayleigh.density, rayleigh.sound_speed),
        ] {
            assert_eq!(density.to_bits(), gas.density.to_bits());
            assert_eq!(speed.to_bits(), gas.sound_speed.to_bits());
        }
        let legacy = GasState::try_new(&GasSpec::dry_air_ussa1976(), t, 101325.0).unwrap();
        use fs_couple::cavity_modal::{CavityNeckLoss, HelmholtzModalCavity};
        let geometry = fs_scenario::HelmholtzCavity {
            volume_m3: 0.003,
            neck_radius_m: 0.0073,
            neck_length_m: 0.057,
            inner_end: fs_scenario::NeckEnd::Flanged,
            outer_end: fs_scenario::NeckEnd::Flanged,
        };
        let cavity =
            HelmholtzModalCavity::new(geometry, *gas, CavityNeckLoss::Resistance(0.0)).unwrap();
        let original =
            HelmholtzModalCavity::new(geometry, legacy, CavityNeckLoss::Resistance(0.0)).unwrap();
        close(
            cavity.isolated_omega_rad_s() / original.isolated_omega_rad_s(),
            gas.sound_speed / legacy.sound_speed,
            1e-12,
        );
        let mut original_transport = *gas;
        original_transport.dynamic_viscosity = legacy.dynamic_viscosity;
        original_transport.thermal_conductivity = legacy.thermal_conductivity;
        original_transport.prandtl = legacy.dynamic_viscosity * cp / legacy.thermal_conductivity;
        original_transport.validate_coherence().unwrap();
        for mode in 1..=3 {
            let analytic = (2 * mode - 1) as f64 * gas.sound_speed / (4.0 * 0.57);
            let f = lossless_resonance(gas, mode);
            let f_legacy = lossless_resonance(&legacy, mode);
            close(f, analytic, 1e-12);
            close(f / f_legacy, gas.sound_speed / legacy.sound_speed, 1e-12);
            let bessel = bessel_peak(gas, mode);
            let bessel_legacy = bessel_peak(&legacy, mode);
            let bessel_original_transport = bessel_peak(&original_transport, mode);
            println!(
                "{{\"suite\":\"matdb-pack\",\"bead\":\"M0.1-W1\",\"test\":\"g1_g3_sourced_nasa9_mixture_reaches_duct_radiation_cavity\",\"fixture\":\"{fixture}\",\"inputs\":{{\"T_K\":{t},\"p_Pa\":101325,\"x_water\":{water},\"x_CO2\":{co2}}},\"observed\":{{\"mode\":{mode},\"lossless_Hz\":{f},\"bessel_Hz\":{bessel},\"legacy_bessel_Hz\":{bessel_legacy},\"shift_cents\":{},\"transport_only_cents\":{},\"R_legacy\":{R_USSA_1976},\"R_CODATA\":{},\"R_only_cents\":{},\"NASA_molar_mass_kg_mol\":{},\"legacy_molar_mass_kg_mol\":0.0289644}},\"tolerance\":{{\"frequency_ratio\":1e-12,\"medium_bits\":0}},\"verdict\":\"matched\",\"elapsed_s\":{}}}",
                cents(bessel, bessel_legacy),
                cents(bessel, bessel_original_transport),
                fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K,
                cents(
                    det::sqrt(fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / R_USSA_1976),
                    1.0
                ),
                properties.molar_mass().value(),
                started.elapsed().as_secs_f64()
            );
        }
        if t == 600.0 {
            assert!(
                (gas.gamma / 1.4 - 1.0).abs() >= 0.01,
                "600K must expose frozen-gamma error"
            );
        }
        let receipt = resolved.thermodynamics().receipt();
        let composition = receipt
            .species()
            .iter()
            .zip(receipt.mole_fraction_bits())
            .map(|(id, bits)| {
                format!(
                    "{{\"species\":\"{id}\",\"mole_fraction\":{}}}",
                    f64::from_bits(*bits)
                )
            })
            .collect::<Vec<_>>()
            .join(",");
        let packs=resolved.source_packs().iter().map(|pack|format!(
            "{{\"species\":\"{}\",\"species_pack\":\"{}\",\"nasa9_pack\":\"{}\",\"transport_pack\":\"{}\"}}",
            pack.species,pack.species_pack,pack.nasa9_pack,pack.transport_pack)).collect::<Vec<_>>().join(",");
        println!(
            "{{\"suite\":\"matdb-pack\",\"bead\":\"M0.1-W1\",\"test\":\"g1_g3_sourced_nasa9_mixture_reaches_duct_radiation_cavity\",\"fixture\":\"{fixture}\",\"inputs\":{{\"T_K\":{t},\"p_Pa\":101325,\"composition\":[{composition}],\"packs\":[{packs}]}},\"provenance\":\"NASA/TP-2002-211556 AppendixD; primary dilute-transport citations in original packs; SRM1720 information values, not certified complete air\",\"rungs\":[1,2,3],\"verdict\":\"matched\",\"elapsed_s\":{}}}",
            started.elapsed().as_secs_f64()
        );
    }
}

/// Exercise the real source-backed humid states. The comparison deliberately
/// remains red until Cramer's primary coefficients, convention and uncertainty
/// are available: secondary transcriptions cannot qualify a literature gate.
pub(super) fn humid_sound_speed_grid() -> Result<(), String> {
    let started = std::time::Instant::now();
    let sources = sources();
    let carrier = [1.0 - 0.2093 - 0.00935, 0.2093, 0.00935];
    let mut count = 0;
    for t in [273.15, 283.15, 293.15, 303.15] {
        for p in [75000.0, 90000.0, 101325.0] {
            for rh in [0.0, 0.5, 1.0] {
                let water = rh * fs_material::gas::saturation_pressure_water_pa(t).unwrap() / p;
                for co2 in [0.0, 0.0004, 0.01] {
                    assert!(water <= 0.06);
                    let model = mixture(&sources, carrier, water, co2);
                    let state = if water > 0.0 {
                        model.evaluate_with_liquid_water_limit(&point(t, p))
                    } else {
                        model.evaluate(&point(t, p))
                    }
                    .unwrap();
                    state.state().validate_coherence().unwrap();
                    println!(
                        "{{\"suite\":\"matdb-pack\",\"bead\":\"M0.1-W1\",\"test\":\"g2_sourced_humid_air_sound_speed_vs_cramer\",\"inputs\":{{\"T_K\":{t},\"p_Pa\":{p},\"RH\":{rh},\"x_water\":{water},\"x_CO2\":{co2}}},\"observed\":{{\"NASA9_c_m_s\":{}}},\"reference\":null,\"relative_error\":null,\"tolerance\":0.0005,\"verdict\":\"NO-DATA\",\"reason\":\"source-unverified: primary Cramer1993 coefficients, dry-carrier convention and uncertainty unavailable; no secondary transcription used\",\"carrier\":\"SRM1720 information carrier proportionally diluted by explicit H2O and CO2; Cramer convention unverified\"}}",
                        state.state().sound_speed
                    );
                    count += 1;
                }
            }
        }
    }
    assert_eq!(count, 108);
    let t = 308.15;
    let p = 101325.0;
    let water = fs_material::gas::saturation_pressure_water_pa(t).unwrap() / p;
    let breath = mixture(&sources, carrier, water, 0.04)
        .evaluate_with_liquid_water_limit(&point(t, p))
        .unwrap();
    println!(
        "{{\"suite\":\"matdb-pack\",\"bead\":\"M0.1-W1\",\"test\":\"g2_sourced_humid_air_sound_speed_vs_cramer\",\"fixture\":\"breath\",\"inputs\":{{\"T_K\":{t},\"p_Pa\":{p},\"x_water\":{water},\"x_CO2\":0.04}},\"observed\":{{\"NASA9_c_m_s\":{}}},\"verdict\":\"outside-reference-domain\",\"gate_applied\":false,\"elapsed_s\":{}}}",
        breath.state().sound_speed,
        started.elapsed().as_secs_f64()
    );
    Err(format!(
        "NO-DATA/source-unverified: {count} NASA-9 humid states computed; primary Cramer (1993, JASA 93:2510) coefficients, composition convention and uncertainty are required before the 5e-4 comparison or residual-times-two tightening can pass"
    ))
}
