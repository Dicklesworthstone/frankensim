//! Manufactured local-factor data exercise the original public face and parcel owners.
use super::*;
use fs_gas::material::transport::{
    BinaryThermalDiffusionFactorMaterial, MaterialBinarySoretDufourTransport,
};
use fs_gas::multicomponent::NormalPrimitive;
use fs_gas::multicomponent::dissipation::{TransportFaceReport, TransportLaw};
use fs_qty::Dims;

type FactorLaw = MaterialBinarySoretDufourTransport<ConstantTransport>;
type FactorOperator = ParcelMolecularTransport<FactorLaw>;

fn mole_quantity() -> QuantitySpec {
    quantity(QuantityKind::Composition(CoordinateBasis::MoleFraction))
}
fn reference(t: f64) -> QueryPoint {
    point(t, 1e5)
        .with_quantity("x_A", mole_quantity(), 0.5).unwrap()
        .with_quantity("x_B", mole_quantity(), 0.5).unwrap()
}
fn axes() -> Vec<(SpeciesId, String)> {
    vec![(id("B"), "x_B".into()), (id("A"), "x_A".into())]
}
fn factor_claim(alpha: f64, maximum: f64) -> PropertyClaim {
    let mut claim = scalar("thermal_diffusion_factor", Dims::NONE, alpha);
    claim.validity = domain(maximum)
        .with_quantity("x_A", mole_quantity(), 0.0, 1.0)
        .with_quantity("x_B", mole_quantity(), 0.0, 1.0);
    claim
}
fn factor_source(alpha: f64, reversed: bool) -> BinaryThermalDiffusionFactorMaterial {
    BinaryThermalDiffusionFactorMaterial {
        species: if reversed { [id("B"), id("A")] } else { [id("A"), id("B")] },
        material: card("binary-local-thermal-factor", [factor_claim(alpha, 800.0)]),
        property: PropertyKey::new("thermal_diffusion_factor", Dims::NONE),
        selection: MaterialPropertySelection::SingleClaimOnly,
    }
}
fn binary(mass_b: f64) -> (ThermochemicalGasModel, ConstantTransport) {
    // Reuse the original source inventory; the old heat law is not evaluated.
    let (bank, original) = fixture(mass_b, 0.0, 0.0);
    let law = ConstantTransport::new(
        original.face().gas().source().clone(),
        "manufactured-local-factor-pair",
        1,
        TransportDomain {
            minimum_temperature_k: 200.0,
            maximum_temperature_k: 1000.0,
            minimum_pressure_pa: 100.0,
            maximum_pressure_pa: 1e8,
        },
        TransportCoefficients {
            dynamic_viscosity: DynViscosity::new(1e-5),
            bulk_viscosity: DynViscosity::new(0.0),
            thermal_conductivity: ThermalConductivityQuantityV1::new(CONDUCTIVITY),
        },
        vec![0.0, DIFFUSIVITY, DIFFUSIVITY, 0.0],
    ).unwrap();
    (bank, law)
}
fn sourced(mass_b: f64, source: BinaryThermalDiffusionFactorMaterial) -> (ThermochemicalGasModel, FactorLaw) {
    let (bank, law) = binary(mass_b);
    let law = MaterialBinarySoretDufourTransport::new(
        law, ThermalConductivityConvention::ZeroDiffusiveSpeciesFlux,
        source, axes(), reference(400.0), || false,
    ).unwrap();
    (bank, law)
}
fn face(law: FactorLaw) -> PlanarTransport<FactorLaw> {
    let gas = MixtureGas::new(
        law.source().clone(),
        CaloricDomain::new(Temperature::new(200.0), Temperature::new(900.0)),
        || false,
    ).unwrap();
    PlanarTransport::new(gas, law, || false).unwrap()
}
fn operator(mass_b: f64, alpha: f64, reversed: bool) -> (ThermochemicalGasModel, FactorOperator) {
    let (bank, law) = sourced(mass_b, factor_source(alpha, reversed));
    (bank, ParcelMolecularTransport::new(face(law)))
}
fn flux(
    face: &mut PlanarTransport<FactorLaw>,
    t: f64, p: f64, x: f64, gx: f64, gt: f64,
) -> ([f64; 5], TransportFaceReport) {
    let m = face.gas().source().models().iter()
        .map(|s| s.molar_mass().value()).collect::<Vec<_>>();
    let distance = 0.1;
    let density = |x: f64, t: f64| [
        p * x * m[0] / (R * t), p * (1.0 - x) * m[1] / (R * t),
    ];
    let left = density(x - 0.5 * distance * gx, t - 0.5 * distance * gt);
    let right = density(x + 0.5 * distance * gx, t + 0.5 * distance * gt);
    let primitive = |partial| NormalPrimitive {
        partial_densities: partial, normal_velocity: 0.0,
        transverse_velocity: 0.0, pressure: p,
    };
    let mut output = [0.0; 5];
    let report = face.face_from_primitives(
        primitive(&left), primitive(&right), distance, [0.0; 2],
        MaxwellStefanSolveV1::default(), &mut output, || false,
    ).unwrap();
    (output, report)
}
fn relative(actual: f64, expected: f64) {
    close(actual, expected, 2e-10 * expected.abs().max(1e-14));
}

#[test]
fn local_factor_matches_binary_mass_frame_and_reciprocal_heat_oracles() {
    let mass_b = 2.0 * MASS;
    for alpha in [-0.8, 0.0, 0.75] {
        let (_, law) = sourced(mass_b, factor_source(alpha, false));
        let mut face = face(law);
        for (t, p, x, gx, gt) in [
            (400.0, 1e5, 0.4, 0.0, 200.0),
            (650.0, 1.7e5, 0.7, 0.6, 0.0),
            (500.0, 1.3e5, 0.3, -0.5, 150.0),
        ] {
            let (out, report) = flux(&mut face, t, p, x, gx, gt);
            let mean_mass = x * MASS + (1.0 - x) * mass_b;
            let force = gx + x * (1.0 - x) * alpha * gt / t;
            let concentration = p / (R * t);
            let j = -concentration * DIFFUSIVITY * MASS * mass_b / mean_mass * force;
            let delta_heat = alpha * R * t * mean_mass / (MASS * mass_b);
            let reduced_heat = -CONDUCTIVITY * gt + delta_heat * j;
            close(out[0], j, 2e-10 * j.abs() + 1e-17);
            close(out[1], -j, 2e-10 * j.abs() + 1e-17);
            close(report.reduced_heat_flux_w_m2, reduced_heat,
                2e-10 * reduced_heat.abs() + 1e-10);
            close(report.mass_frame_residual_kg_m2_s, 0.0, 1e-16);
            // The independent conjugate-force entropy equals friction + Fourier.
            let entropy = -R * out[0] * gx * mean_mass / (MASS * mass_b * x * (1.0 - x))
                - report.reduced_heat_flux_w_m2 * gt / (t * t);
            let pair = R * out[0].powi(2) * mean_mass.powi(2)
                / (concentration * DIFFUSIVITY * MASS.powi(2) * mass_b.powi(2) * x * (1.0 - x));
            relative(entropy, pair + CONDUCTIVITY * (gt / t).powi(2));
            assert!(entropy > 0.0);
            let resolved = face.law().last_resolved();
            relative(resolved.heats()[0].value(), delta_heat);
            assert_eq!(resolved.heats()[1].value(), 0.0);
            assert_eq!(resolved.material().card_identity(), face.law().factor_source().material.content_hash());
        }
    }
}

#[test]
fn measured_null_force_and_pair_orientation_use_one_sign_convention() {
    let (t, p, x, gt, alpha) = (500.0, 1e5, 0.4, 250.0, 0.75);
    let gx = -x * (1.0 - x) * alpha * gt / t;
    let (_, forward) = sourced(2.0 * MASS, factor_source(alpha, false));
    let (_, reverse) = sourced(2.0 * MASS, factor_source(-alpha, true));
    let mut forward = face(forward);
    let mut reverse = face(reverse);
    let (out, report) = flux(&mut forward, t, p, x, gx, gt);
    close(out[0], 0.0, 1e-16);
    close(out[1], 0.0, 1e-16);
    relative(report.reduced_heat_flux_w_m2, -CONDUCTIVITY * gt);
    // Off equilibrium, reversing the species pair and alpha leaves every flux.
    for (gx, gt) in [(0.0, 250.0), (0.7, 0.0), (-0.4, -150.0)] {
        let (a, _) = flux(&mut forward, t, p, x, gx, gt);
        let (b, _) = flux(&mut reverse, t, p, x, gx, gt);
        for (a, b) in a.into_iter().zip(b) {
            close(a, b, 2e-10 * a.abs().max(1e-14));
        }
    }
}

#[test]
fn source_surface_refreshes_actual_temperature_pressure_and_composition() {
    let mut claim = factor_claim(0.0, 800.0);
    let alpha = |t: f64, x: f64| 0.2 + 0.001 * t + 0.4 * x;
    claim.value = PropertyValue::Surface2D {
        axes: ["temperature".into(), "x_A".into()],
        axis_dims: [Temperature::DIMS, Dims::NONE],
        samples: [(200.0, 0.0), (800.0, 0.0), (200.0, 1.0), (800.0, 1.0)]
            .map(|(t, x)| ([t, x], alpha(t, x))).to_vec(),
        triangles: vec![[0, 1, 2], [1, 3, 2]],
        dims: Dims::NONE,
    };
    claim.interpolation = InterpolationPolicy::LinearInside;
    let mut source = factor_source(0.0, false);
    source.material = card("binary-factor-surface", [claim]);
    let (_, mut law) = sourced(2.0 * MASS, source);
    for (t, p, x) in [(300.0, 8e4, 0.3), (700.0, 2e5, 0.7)] {
        let mut heat = [0.0; 2];
        law.heats_of_transport(
            Temperature::new(t), Pressure::new(p), &[x, 1.0 - x],
            &mut heat, &mut || false,
        ).unwrap();
        let state = law.last_resolved();
        relative(state.factor(), alpha(t, x));
        relative(heat[0], alpha(t, x) * R * t * (x / (2.0 * MASS) + (1.0 - x) / MASS));
        assert_eq!(state.material().query().axes()["pressure"], p);
        relative(state.material().query().axes()["x_A"], x);
        assert_eq!(state.material().query().axis_quantities()["x_A"], mole_quantity());
    }
}

#[test]
fn source_pin_domain_and_cancellation_refusals_do_not_publish_coefficients() {
    let low = factor_claim(0.75, 450.0);
    let mut high = factor_claim(-1.0, 800.0);
    high.validity = high.validity.with_quantity(
        "temperature", quantity(QuantityKind::AbsoluteTemperature), 600.0, 800.0,
    );
    let mut source = factor_source(0.75, false);
    source.material = card("disjoint-factor-claims", [low, high]);
    let (_, base) = binary(2.0 * MASS);
    assert!(MaterialBinarySoretDufourTransport::new(
        base, ThermalConductivityConvention::ZeroDiffusiveSpeciesFlux,
        source.clone(), axes(), reference(700.0), || false,
    ).is_ok());
    let (_, mut law) = sourced(2.0 * MASS, source);
    let original = law.last_resolved().clone();
    let mut heat = [17.0; 2];
    let mut diffusion = [19.0; 4];
    for (t, p, x) in [
        (700.0, 1e5, [0.5, 0.5]),
        (400.0, 1e9, [0.5, 0.5]),
        (400.0, 1e5, [0.4, 0.5]),
    ] {
        assert!(law.heats_of_transport(
            Temperature::new(t), Pressure::new(p), &x, &mut heat, &mut || false,
        ).is_err());
        assert!(law.evaluate(
            Temperature::new(t), Pressure::new(p), &x, &mut diffusion, &mut || false,
        ).is_err());
        assert_eq!(heat, [17.0; 2]);
        assert_eq!(diffusion, [19.0; 4]);
        assert_eq!(law.last_resolved(), &original);
    }
    assert!(matches!(law.heats_of_transport(
        Temperature::new(410.0), Pressure::new(1e5), &[0.4, 0.6],
        &mut heat, &mut || true,
    ), Err(DissipationError::Cancelled)));
    assert_eq!(heat, [17.0; 2]);
    assert_eq!(law.last_resolved(), &original);
}

#[test]
fn exact_absent_species_are_unseeded_and_bad_factor_assignments_refuse() {
    let (_, law) = sourced(2.0 * MASS, factor_source(0.75, false));
    let mut face = face(law);
    for x in [0.0, 1.0] {
        let (out, report) = flux(&mut face, 400.0, 1e5, x, 0.0, 200.0);
        assert_eq!(out[0], 0.0);
        assert_eq!(out[1], 0.0);
        relative(report.reduced_heat_flux_w_m2, -CONDUCTIVITY * 200.0);
        assert_eq!(face.law().last_resolved().mole_fractions(), &[x, 1.0 - x]);
    }
    for bad in 0..3 {
        let (_, base) = binary(2.0 * MASS);
        let mut source = factor_source(0.75, false);
        match bad {
            0 => source.species[1] = id("A"),
            1 => source.species[1] = id("missing"),
            _ => source.property = PropertyKey::new("soret_per_kelvin", Dims([0, 0, 0, -1, 0, 0])),
        }
        assert!(MaterialBinarySoretDufourTransport::new(
            base, ThermalConductivityConvention::ZeroDiffusiveSpeciesFlux,
            source, axes(), reference(400.0), || false,
        ).is_err());
    }
}


#[test]
fn binary_source_refuses_a_third_original_caloric_species() {
    let (_, binary) = binary(2.0 * MASS);
    let mut models = binary.source().models().to_vec();
    models.push(species("C", 3.0 * MASS, 0.0));
    let ids = models.iter().map(|s| s.species().clone()).collect::<Vec<_>>();
    let source = IdealGasEquilibriumV1::new(
        ElementalMatrix::new(
            vec![ElementId::new("H").unwrap()], ids.clone(), vec![vec![1, 2, 3]],
        ).unwrap(),
        ChargeVector::new(ids, vec![0, 0, 0]).unwrap(), models,
    ).unwrap();
    let law = ConstantTransport::new(
        source, "manufactured-three-species", 1,
        TransportDomain {
            minimum_temperature_k: 200.0, maximum_temperature_k: 1000.0,
            minimum_pressure_pa: 100.0, maximum_pressure_pa: 1e8,
        },
        TransportCoefficients {
            dynamic_viscosity: DynViscosity::new(1e-5),
            bulk_viscosity: DynViscosity::new(0.0),
            thermal_conductivity: ThermalConductivityQuantityV1::new(CONDUCTIVITY),
        },
        vec![0.0, DIFFUSIVITY, DIFFUSIVITY, DIFFUSIVITY, 0.0,
             DIFFUSIVITY, DIFFUSIVITY, DIFFUSIVITY, 0.0],
    ).unwrap();
    assert!(MaterialBinarySoretDufourTransport::new(
        law, ThermalConductivityConvention::ZeroDiffusiveSpeciesFlux,
        factor_source(0.75, false), axes(), reference(400.0), || false,
    ).is_err());
}

#[test]
fn local_factor_drives_closed_parcels_with_reciprocal_energy_and_source_receipts() {
    let alpha = 0.75;
    let (bank, mut op) = operator(MASS, alpha, false);
    let original = BTreeMap::from([
        (1, parcel(&bank, 0.5, 390.0, 1e5)),
        (2, parcel(&bank, 0.5, 410.0, 1e5)),
    ]);
    let endpoints = original.keys().map(|&id| (id, endpoint())).collect();
    let next = op.step_candidate(
        &original, &[link()], GasTransportPaths::default(), 0.1,
        &endpoints, controls(), || false,
    ).unwrap();
    let j = -1e5 / (R * 400.0) * DIFFUSIVITY * MASS * 0.25 * alpha * 200.0 / 400.0;
    let heat = -CONDUCTIVITY * 200.0 + alpha * R * 400.0 / MASS * j;
    close(next.step().links[0].mass_left_to_right_kg[0], j * 0.001, 1e-18);
    relative(next.step().links[0].face.reduced_heat_flux_w_m2, heat);
    for (id, sign) in [(1, -1.0), (2, 1.0)] {
        let old = &original[&id];
        let new = &next.parcels()[&id];
        close(mass(new, "A") - mass(old, "A"), sign * j * 0.001, 1e-18);
        close(new.mass_kg(), old.mass_kg(), 1e-18);
        close(new.state().state().temperature,
            old.state().state().temperature + sign * heat * 0.001 / (old.mass_kg() * 3.0 * R / MASS),
            3e-10);
    }
    let total = |parcels: &BTreeMap<u64, ClosedGasParcel>| parcels.values()
        .map(|p| p.mass_kg() * p.state().specific_internal_energy_j_kg()).sum::<f64>();
    close(total(next.parcels()), total(&original), 1e-9);
    assert!(next.step().transport.entropy_change_j_k > 0.0);
    let resolved = op.face().law().last_resolved();
    close(resolved.temperature().value(), next.parcels()[&2].state().state().temperature, 1e-12);
    relative(resolved.heats()[0].value(), alpha * R * resolved.temperature().value() / MASS);
    assert_eq!(resolved.material().card_identity(), op.face().law().factor_source().material.content_hash());
}
