//! Synthetic source cards exercise the actual material -> bow -> body path.
//! These are not measured alloys, wrapped strings or calibrated instruments.
use fs_couple::acoustic_realize::AcousticRealizeError;
use fs_couple::air_path::oscillating_cylinder_air_resistance_per_length;
use fs_couple::bowed_string::{
    BowGesture, BowedRunConfig, BowedStringCard, FrictionIsland, Termination,
};
use fs_couple::bowed_string::runtime::BowedStringState;
use fs_couple::string_specimen::{
    KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY as ETA, ResolvedStringSpecimen,
    StringPrestress, with_uniform_circular_material_state,
};
use fs_couple::thin_plate::CompactBody;
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
use fs_material::gas::{ConductivityModel, GasSpec, GasState};
use fs_material::section::{ConcentricElasticSection, ConcentricKelvinVoigtSection, ConcentricMaterialLayer};
use fs_material::state_point::{MaterialPropertySelection, ResolvedMaterialStatePoint,
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point};
use fs_qty::{Density, Dims, DynViscosity, Pressure};
use fs_scenario::{AmbientGas, PrestressedString, RadiatingPlate};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>, &CancelGate) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| f(&Cx::new(&gate, arena,
        StreamKey { seed: 0, kernel_id: 81, tile: 0, iteration: 0 },
        Budget::INFINITE, ExecMode::Deterministic), &gate))
}
fn template() -> PrestressedString {
    PrestressedString {
        length_m: 0.65, tension_n: 60.0, lin_density_kg_m: 1.0,
        axial_stiffness_n: 0.0, width_m: 0.01, n_modes: 4,
        damping_ratio: 0.0, rayleigh: None, bending_stiffness_n_m2: 0.0,
        kelvin_voigt_bending: None, relaxation_bending: None,
        polarization_detune: 0.0, moving_end: false, tension_loss_factor: None, bending_loss_factor: None,
    }
}
fn source(rho: f64, young: f64, eta: f64, band: (f64, f64)) -> ResolvedMaterialStatePoint {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (name, value, dims) in [("density", rho, Density::DIMS),
        ("young_modulus", young, Pressure::DIMS), (ETA, eta, DynViscosity::DIMS)] {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims), value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained().with("T", 290.0, 300.0)
                .with("omega", band.0, band.1),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            provenance: Provenance { source: "synthetic bowed-material regression".into(),
                license: "CC0-1.0".into(), artifact: None }, observations: vec![],
        }).unwrap();
        requirements.push(ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::Finite).unwrap());
    }
    let card = MaterialCard::assemble(MaterialStateId {
        chemistry: "synthetic".into(), phase: "solid".into(), process: "test".into(), revision: 0,
    }, claims, vec![]).unwrap();
    let point = QueryPoint::new().with("T", 293.15).unwrap().with("omega", band.0).unwrap();
    resolve_material_state_point(&card, &point, &requirements,
        MaterialPropertySelection::SingleClaimOnly).unwrap()
}
fn specimen(state: &ResolvedMaterialStatePoint) -> ResolvedStringSpecimen {
    with_uniform_circular_material_state(template(), 0.0005, state).unwrap()
        .with_kelvin_voigt_bending_loss().unwrap()
}
fn section(core: ResolvedMaterialStatePoint, skin: ResolvedMaterialStatePoint)
    -> ConcentricKelvinVoigtSection {
    ConcentricElasticSection::try_new(0.0, vec![
        ConcentricMaterialLayer { outer_radius_m: 0.0003, material: core },
        ConcentricMaterialLayer { outer_radius_m: 0.0005, material: skin },
    ]).unwrap().with_kelvin_voigt_bending(ETA).unwrap()
}
fn gas() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap()
}
fn close(a: f64, b: f64) {
    assert!((a-b).abs() < 2e-12 * a.abs().max(b.abs()).max(1e-30), "{a:e} != {b:e}");
}
fn compare(a: &BowedStringCard, b: &BowedStringCard) {
    close(a.length_m, b.length_m);
    close(a.tension_n, b.tension_n);
    close(a.linear_density_kg_m, b.linear_density_kg_m);
    close(a.bending_stiffness_n_m2, b.bending_stiffness_n_m2);
    close(a.viscous_bending_n_m2_s, b.viscous_bending_n_m2_s);
    assert_eq!(a.sample_rate_hz, b.sample_rate_hz);
    assert_eq!(a.mode_count, b.mode_count);
    for k in 0..a.mode_count {
        close(a.zetas[k], b.zetas[k]);
        close(a.mode_omega_rad_s(k+1), b.mode_omega_rad_s(k+1));
    }
}
fn response(card: BowedStringCard, gas: GasState) -> (Vec<f64>, Vec<f64>) {
    let config = BowedRunConfig {
        card, island: FrictionIsland::ViscousOnly { viscous_n_s_per_m: 0.0001 },
        gesture: BowGesture::admit(0.2, 1.0, 0.2).unwrap(), steps: 128, subsamples: 4,
        termination: Termination::PlateOnePort {
            body: Box::new(CompactBody::from_radiator(RadiatingPlate {
                area_m2: 0.08, mass_kg: 0.2, frequency_hz: 220.0, damping_ratio: 0.02,
            }).unwrap()), ambient: gas,
        // A 0.1 m listener hears this 128-sample window: the time of flight
        // is about 14 samples (1 m would be about 140, past the window).
        }, listener_m: 0.1,
    };
    let mut runtime = BowedStringState::new(&config, 128).unwrap();
    let (mut shear, mut pressure) = (Vec::new(), Vec::new());
    for index in 0..128 {
        if index == 64 { runtime.set_bow(0.0, 0.0, 0.2).unwrap(); }
        let sample = runtime.step().unwrap();
        shear.push(sample.bridge_force_n);
        pressure.push(sample.radiated_pressure_pa.unwrap());
    }
    assert!(shear.iter().chain(&pressure).all(|x| x.is_finite()));
    assert!(pressure.iter().any(|p| *p != 0.0));
    (shear, pressure)
}

#[test]
fn ambient_and_explicit_gas_share_identical_material_lowering() {
    with_cx(|cx, _| {
        let s = specimen(&source(1000.0, 2e9, 3e7, (1.0, 20_000.0)));
        for rh in [0.0, 0.5] {
            let ambient = AmbientGas { temperature_k: 293.15, pressure_pa: 101_325.0, relative_humidity: rh };
            let gas = GasState::try_new_moist_air(293.15, 101_325.0, rh).unwrap();
            let a = BowedStringCard::from_material_specimen(cx, &s, ambient, 48_000).unwrap();
            let b = BowedStringCard::from_material_specimen_with_gas(cx, &s, &gas, 48_000).unwrap();
            compare(&a, &b);
            assert_eq!(a.zetas, b.zetas);
            assert_eq!(response(a, gas), response(b, gas));
        }
    });
}

#[test]
fn homogeneous_and_identically_layered_sources_reduce_to_the_same_bow() {
    with_cx(|cx, _| {
        let material = source(1000.0, 2e9, 3e7, (1.0, 20_000.0));
        let a = BowedStringCard::from_material_specimen_with_gas(cx, &specimen(&material), &gas(), 48_000).unwrap();
        let layers = section(material.clone(), material);
        let b = BowedStringCard::from_concentric_material_section(cx, template(), &layers,
            StringPrestress::FixedTension(60.0), &gas(), 48_000).unwrap();
        compare(&a, &b);
        let (a_shear, a_pressure) = response(a, gas());
        let (b_shear, b_pressure) = response(b, gas());
        for (a, b) in a_shear.iter().zip(&b_shear).chain(a_pressure.iter().zip(&b_pressure)) {
            assert!((a-b).abs() < 1e-10 * a.abs().max(b.abs()).max(1e-12));
        }
    });
}

#[test]
fn independent_layer_density_elasticity_and_loss_changes_reach_bridge_and_sound() {
    with_cx(|cx, _| {
        let mut cards = Vec::new();
        let mut signals = Vec::new();
        for (rho, young, eta) in [(4000.0, 4e9, 0.0), (8000.0, 4e9, 0.0),
            (4000.0, 8e9, 0.0), (4000.0, 4e9, 3e7)] {
            let core = source(1000.0, 2e9, 0.0, (1.0, 20_000.0));
            let skin = source(rho, young, eta, (1.0, 20_000.0));
            let layers = section(core, skin);
            let before = layers.clone();
            let input = template(); let snapshot = input.clone();
            let card = BowedStringCard::from_concentric_material_section(cx, input, &layers,
                StringPrestress::FixedTension(60.0), &gas(), 48_000).unwrap();
            assert_eq!(layers, before);
            assert_eq!(template(), snapshot);
            close(card.linear_density_kg_m, layers.elastic().line_density_kg_m());
            close(card.bending_stiffness_n_m2, layers.elastic().bending_stiffness_n_m2());
            assert_eq!(card.viscous_bending_n_m2_s, layers.viscous_stiffness_n_m2_s());
            for k in 1..=4 {
                let wave = k as f64 * core::f64::consts::PI / card.length_m;
                let omega = card.mode_omega_rad_s(k);
                let drag = oscillating_cylinder_air_resistance_per_length(0.0005, omega, &gas()).unwrap();
                close(card.zetas[k-1], (drag + layers.viscous_stiffness_n_m2_s() * wave.powi(4))
                    / (2.0 * card.linear_density_kg_m * omega));
            }
            signals.push(response(card.clone(), gas()));
            cards.push(card);
        }
        assert!(cards[1].fundamental_hz() < cards[0].fundamental_hz());
        assert!(cards[2].mode_omega_rad_s(4) > cards[0].mode_omega_rad_s(4));
        assert_eq!(cards[3].fundamental_hz(), cards[0].fundamental_hz());
        assert!(cards[3].zetas.iter().zip(&cards[0].zetas).all(|(a,b)| a > b));
        for signal in &signals[1..] {
            assert_ne!(signal.0, signals[0].0);
            assert_ne!(signal.1, signals[0].1);
        }
    });
}

#[test]
fn arbitrary_gas_and_pressure_change_drag_and_rendering_without_retuning_solid_modes() {
    with_cx(|cx, _| {
        let s = specimen(&source(1000.0, 2e9, 0.0, (1.0, 20_000.0)));
        // An explicitly synthetic gas, not an unsupported claim about helium.
        let synthetic = GasSpec { molar_mass: 0.01, gamma: 5.0/3.0,
            sutherland_beta: 9e-7, sutherland_s: 100.0, conductivity: ConductivityModel::Eucken };
        let gases = [gas(), GasState::try_new(&GasSpec::dry_air_ussa1976(),293.15,50_000.0).unwrap(),
            GasState::try_new(&synthetic,293.15,101_325.0).unwrap()];
        let mut cards = Vec::new(); let mut signals = Vec::new();
        for gas in gases {
            let card = BowedStringCard::from_material_specimen_with_gas(cx, &s, &gas, 48_000).unwrap();
            for k in 1..=4 {
                let omega = card.mode_omega_rad_s(k);
                let drag = oscillating_cylinder_air_resistance_per_length(s.radius_m(), omega, &gas).unwrap();
                close(card.zetas[k-1], drag / (2.0 * card.linear_density_kg_m * omega));
            }
            signals.push(response(card.clone(), gas)); cards.push(card);
        }
        for k in 1..=4 {
            assert_eq!(cards[0].mode_omega_rad_s(k), cards[1].mode_omega_rad_s(k));
            assert_eq!(cards[0].mode_omega_rad_s(k), cards[2].mode_omega_rad_s(k));
            assert!(cards[1].zetas[k-1] < cards[0].zetas[k-1]);
        }
        for i in 1..3 {
            assert_ne!(cards[0].zetas, cards[i].zetas);
            assert_ne!(signals[0].0, signals[i].0);
            assert_ne!(signals[0].1, signals[i].1);
        }
    });
}

#[test]
fn source_band_support_loss_and_gas_refusals_remain_atomic() {
    with_cx(|cx, gate| {
        let material = source(1000.0, 2e9, 0.0, (1.0,20_000.0));
        let s = specimen(&material);
        let layers = section(material.clone(), material);
        let compile = |input, gas: &GasState, rate| BowedStringCard::from_concentric_material_section(
            cx, input, &layers, StringPrestress::FixedTension(60.0), gas, rate);
        let mut bad_gas = gas(); bad_gas.sound_speed *= 2.0;
        assert!(compile(template(), &bad_gas, 48_000).is_err());
        assert!(BowedStringCard::from_material_specimen_with_gas(cx, &s, &bad_gas, 48_000).is_err());
        assert!(compile(template(), &gas(), 100).is_err());
        for kind in 0..4 {
            let mut input = template();
            match kind { 0 => input.damping_ratio = 0.01, 1 => input.moving_end = true,
                2 => input.polarization_detune = 0.01, _ => input.n_modes = 0 }
            assert!(compile(input, &gas(), 48_000).is_err());
        }
        let narrow = section(source(1000.0,2e9,0.0,(1.0,2000.0)),
            source(1000.0,2e9,0.0,(1.0,2000.0)));
        assert!(BowedStringCard::from_concentric_material_section(cx, template(), &narrow,
            StringPrestress::FixedTension(60.0), &gas(), 48_000).is_err());
        let before = s.clone();
        compile(template(), &gas(), 48_000).unwrap();
        gate.request();
        assert!(matches!(compile(template(), &gas(), 48_000), Err(AcousticRealizeError::Cancelled)));
        assert!(matches!(BowedStringCard::from_material_specimen_with_gas(cx, &s, &gas(),48_000),
            Err(AcousticRealizeError::Cancelled)));
        assert_eq!(s, before);
    });
}
