use super::*;
use crate::multicomponent::{NormalPrimitive, SpatialPrimitive};
use fs_evidence::ValidityDomain;
use fs_thermochem::{
    MassSpecificEnergyQuantityV1, Nasa9RegionV1, Nasa9StandardStateModelV1, SpeciesId,
};

const MASSES: [f64; 2] = [0.002, 0.028];
const HEATS: [f64; 2] = [125000.0, -75000.0];
const DIFFUSIVITY: f64 = 1e-4;
const CONDUCTIVITY: f64 = 0.04;
const PRESSURE: f64 = 1e5;

fn source_gas() -> MixtureGas {
    let models = [("A", MASSES[0], 3.5), ("B", MASSES[1], 4.0)]
        .into_iter()
        .map(|(name, mass, cp)| {
            let original = source(name, mass, cp, 0.0, 0.0);
            let mut card = original.regions()[0].card().clone();
            card.validity = ValidityDomain::unconstrained().with("T", 200.0, 1500.0);
            Nasa9StandardStateModelV1::new(
                original.species().clone(),
                original.molar_mass(),
                original.convention().clone(),
                vec![Nasa9RegionV1::from_card(card).unwrap()],
            )
            .unwrap()
        })
        .collect();
    MixtureGas::new(
        bank(models),
        CaloricDomain::new(Temperature::new(200.0), Temperature::new(1500.0)),
        || false,
    )
    .unwrap()
}
fn domain() -> TransportDomain {
    TransportDomain {
        minimum_temperature_k: 200.0,
        maximum_temperature_k: 1500.0,
        minimum_pressure_pa: 1e4,
        maximum_pressure_pa: 1e6,
    }
}
fn molecular(g: &MixtureGas) -> ConstantTransport {
    ConstantTransport::new(
        g.source().clone(),
        "manufactured-molecular",
        1,
        domain(),
        TransportCoefficients {
            dynamic_viscosity: DynViscosity::new(0.0),
            bulk_viscosity: DynViscosity::new(0.0),
            thermal_conductivity: ThermalConductivityQuantityV1::new(CONDUCTIVITY),
        },
        vec![0.0, DIFFUSIVITY, DIFFUSIVITY, 0.0],
    )
    .unwrap()
}
fn heats(values: [f64; 2]) -> Vec<(SpeciesId, MassSpecificEnergyQuantityV1)> {
    // Reversed authored order exercises canonical source binding.
    vec![
        (
            SpeciesId::new("B").unwrap(),
            MassSpecificEnergyQuantityV1::new(values[1]),
        ),
        (
            SpeciesId::new("A").unwrap(),
            MassSpecificEnergyQuantityV1::new(values[0]),
        ),
    ]
}
fn thermal(g: &MixtureGas, values: [f64; 2]) -> SoretDufourTransport<ConstantTransport> {
    SoretDufourTransport::new(
        molecular(g),
        "manufactured-cross-coefficients",
        1,
        domain(),
        heats(values),
    )
    .unwrap()
}
fn densities(x: f64, temperature: f64) -> [f64; 2] {
    [
        PRESSURE * x * MASSES[0] / (R * temperature),
        PRESSURE * (1.0 - x) * MASSES[1] / (R * temperature),
    ]
}
fn primitive(densities: &[f64]) -> NormalPrimitive<'_> {
    NormalPrimitive {
        partial_densities: densities,
        normal_velocity: 0.0,
        transverse_velocity: 0.0,
        pressure: PRESSURE,
    }
}
fn near(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 2e-10 * b.abs().max(1e-18),
        "actual={a:e}, expected={b:e}"
    );
}

#[test]
fn g1_soret_dufour_binary_flux_and_onsager_entropy_match_independent_equations() {
    let g = source_gas();
    let mut face = PlanarTransport::new(g.clone(), thermal(&g, HEATS), || false).unwrap();
    let capacity = face.workspace_bytes();
    for (xl, xr, tl, tr) in [
        (0.4, 0.4, 450.0, 550.0),
        (0.3, 0.5, 500.0, 500.0),
        (0.3, 0.5, 450.0, 550.0),
        (0.5, 0.3, 550.0, 450.0),
        (0.4018, 0.3982, 450.0, 550.0),
    ] {
        let left = densities(xl, tl);
        let right = densities(xr, tr);
        let mut flux = [0.0; 5];
        let report = face
            .face_from_primitives(
                primitive(&left),
                primitive(&right),
                0.1,
                [0.0; 2],
                MaxwellStefanSolveV1::default(),
                &mut flux,
                || false,
            )
            .unwrap();
        let x = 0.5 * (xl + xr);
        let t = 0.5 * (tl + tr);
        let gx = (xr - xl) / 0.1;
        let gt = (tr - tl) / 0.1;
        let mean_mass = x * MASSES[0] + (1.0 - x) * MASSES[1];
        let y = x * MASSES[0] / mean_mass;
        let mean_heat = y * HEATS[0] + (1.0 - y) * HEATS[1];
        let driving = gx + x * MASSES[0] * (HEATS[0] - mean_heat) / (R * t) * (gt / t);
        let c = PRESSURE / (R * t);
        let j = -c * DIFFUSIVITY * driving * MASSES[0] * MASSES[1] / mean_mass;
        near(flux[0], j);
        near(flux[1], -j);
        if gx < 0.0 && gt > 0.0 {
            // An admitted thermal force can overcome ordinary mixing and drive
            // the light species up its concentration gradient without negative
            // total entropy production.
            assert!(j < 0.0);
        }
        let enthalpies = [
            (2.5 * R * (t - 200.0) + R * t) / MASSES[0],
            (3.0 * R * (t - 200.0) + R * t) / MASSES[1],
        ];
        let reduced_heat = -CONDUCTIVITY * gt + (HEATS[0] - HEATS[1]) * j;
        near(flux[4], (enthalpies[0] - enthalpies[1]) * j + reduced_heat);
        let entropy = -R * (j / MASSES[0] * gx / x + (-j / MASSES[1]) * (-gx) / (1.0 - x))
            - reduced_heat * gt / (t * t);
        let relative_velocity = j / (c * x * MASSES[0]) + j / (c * (1.0 - x) * MASSES[1]);
        let pair_entropy =
            R * c * x * (1.0 - x) / DIFFUSIVITY * relative_velocity * relative_velocity;
        near(entropy, pair_entropy + CONDUCTIVITY * (gt / t).powi(2));
        assert!(entropy > 0.0);
        assert!(report.maximum_gradient_residual_per_m < 1e-12);
        assert!(report.mass_frame_residual_kg_m2_s.abs() < 1e-16);
    }
    assert_eq!(face.workspace_bytes(), capacity);
}

#[test]
fn g3_thermal_heat_gauge_zero_coupling_and_spatial_reduction_preserve_physics() {
    let g = source_gas();
    let mut plain = PlanarTransport::new(g.clone(), molecular(&g), || false).unwrap();
    let mut zero = PlanarTransport::new(g.clone(), thermal(&g, [0.0, 0.0]), || false).unwrap();
    let mut coupled = PlanarTransport::new(g.clone(), thermal(&g, HEATS), || false).unwrap();
    let mut shifted = PlanarTransport::new(
        g.clone(),
        thermal(&g, [HEATS[0] + 1e12, HEATS[1] + 1e12]),
        || false,
    )
    .unwrap();
    let mut spatial = SpatialTransport::new(g.clone(), thermal(&g, HEATS), || false).unwrap();
    let left = densities(0.3, 450.0);
    let right = densities(0.5, 550.0);
    let mut a = [0.0; 5];
    let mut b = [0.0; 5];
    let mut c = [0.0; 5];
    let mut d = [0.0; 5];
    let controls = MaxwellStefanSolveV1::default();
    let plain_report = plain
        .face_from_primitives(
            primitive(&left),
            primitive(&right),
            0.1,
            [0.0; 2],
            controls,
            &mut a,
            || false,
        )
        .unwrap();
    zero.face_from_primitives(
        primitive(&left),
        primitive(&right),
        0.1,
        [0.0; 2],
        controls,
        &mut b,
        || false,
    )
    .unwrap();
    let coupled_report = coupled
        .face_from_primitives(
            primitive(&left),
            primitive(&right),
            0.1,
            [0.0; 2],
            controls,
            &mut c,
            || false,
        )
        .unwrap();
    shifted
        .face_from_primitives(
            primitive(&left),
            primitive(&right),
            0.1,
            [0.0; 2],
            controls,
            &mut d,
            || false,
        )
        .unwrap();
    for k in 0..5 {
        near(a[k], b[k]);
        near(c[k], d[k]);
    }
    assert!(coupled_report.diffusivity_m2_s > plain_report.diffusivity_m2_s);
    let mut three = [0.0; 6];
    spatial
        .face_from_primitives(
            SpatialPrimitive {
                partial_densities: &left,
                normal_velocity: 0.0,
                transverse_velocity: [0.0; 2],
                pressure: PRESSURE,
            },
            SpatialPrimitive {
                partial_densities: &right,
                normal_velocity: 0.0,
                transverse_velocity: [0.0; 2],
                pressure: PRESSURE,
            },
            0.1,
            [[0.0; 3]; 2],
            controls,
            &mut three,
            || false,
        )
        .unwrap();
    near(three[0], c[0]);
    near(three[1], c[1]);
    near(three[5], c[4]);
    let left = densities(0.0, 450.0);
    let right = densities(0.0, 550.0);
    coupled
        .face_from_primitives(
            primitive(&left),
            primitive(&right),
            0.1,
            [0.0; 2],
            controls,
            &mut c,
            || false,
        )
        .unwrap();
    assert_eq!(&c[..2], &[0.0, 0.0]);
    near(c[4], -CONDUCTIVITY * 1000.0);
}

#[test]
fn g0_soret_source_axis_domain_and_cancellation_are_enforced() {
    let g = source_gas();
    let mut wrong = heats(HEATS);
    wrong[0].0 = SpeciesId::new("foreign").unwrap();
    assert!(SoretDufourTransport::new(molecular(&g), "wrong-species", 1, domain(), wrong).is_err());
    assert!(
        SoretDufourTransport::new(
            molecular(&g),
            "bad-coefficient",
            1,
            domain(),
            heats([f64::NAN, 0.0])
        )
        .is_err()
    );
    let mut restricted = domain();
    restricted.maximum_temperature_k = 499.0;
    let law =
        SoretDufourTransport::new(molecular(&g), "bounded", 1, restricted, heats(HEATS)).unwrap();
    let mut face = PlanarTransport::new(g.clone(), law, || false).unwrap();
    let left = densities(0.4, 450.0);
    let right = densities(0.4, 550.0);
    let mut out = [0.0; 5];
    assert!(
        face.face_from_primitives(
            primitive(&left),
            primitive(&right),
            0.1,
            [0.0; 2],
            MaxwellStefanSolveV1::default(),
            &mut out,
            || false
        )
        .is_err()
    );
    let mut law = thermal(&g, HEATS);
    let mut output = [7.0; 2];
    assert_eq!(
        law.heats_of_transport(
            Temperature::new(500.0),
            fs_qty::Pressure::new(PRESSURE),
            &[0.4, 0.6],
            &mut output,
            &mut || true
        )
        .unwrap_err(),
        DissipationError::Cancelled
    );
    assert_eq!(output, [7.0; 2]);
    assert!(
        law.heats_of_transport(
            Temperature::new(500.0),
            fs_qty::Pressure::new(PRESSURE),
            &[0.4, 0.6],
            &mut output,
            &mut || false
        )
        .unwrap()
    );
    assert_eq!(output, HEATS);
}

#[test]
fn g0_soret_field_evolves_composition_and_conserves_all_species_and_energy() {
    let g = source_gas();
    let mut flow = field(g.clone(), 24, Boundaries::Periodic, |x| {
        let temperature = 500.0 + 30.0 * (2.0 * core::f64::consts::PI * x).sin();
        (densities(0.4, temperature).to_vec(), 0.0, PRESSURE)
    });
    let before = flow.totals();
    let first = flow.cell(0).unwrap().partial_densities.to_vec();
    let mut operator = DissipationOperator::new(&flow, thermal(&g, HEATS), || false).unwrap();
    let report = operator
        .advance_for(
            &mut flow,
            Time::new(0.1),
            DissipationControl::default(),
            || false,
        )
        .unwrap();
    let after = flow.totals();
    for k in 0..2 {
        near(after.species_masses[k], before.species_masses[k]);
    }
    near(after.energy, before.energy);
    assert!(report.absolute_energy_residual_j.abs() < 1e-8);
    assert!((flow.cell(0).unwrap().partial_densities[0] - first[0]).abs() > 1e-9);
    assert_eq!(report.elapsed.value(), 0.1);
    assert_eq!(flow.time().value(), 0.0);
}
