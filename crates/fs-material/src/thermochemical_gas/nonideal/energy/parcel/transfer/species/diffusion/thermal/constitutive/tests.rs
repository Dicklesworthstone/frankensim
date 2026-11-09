use super::*;
use fs_evidence::ValidityDomain;
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance};
use fs_qty::{Density, MolarMass, Pressure};
use fs_thermochem::{
    ElementalReferenceIdV1, NASA9_LAW_ID_V1, NASA9_LAW_VERSION_V1, NASA9_STATE_SCHEMA_VERSION_V1,
    Nasa9RegionV1, ReferenceEquationOfStateV1, SpeciesId, StandardStateConventionV1,
    StandardStatePhaseV1,
    eos::peng_robinson::{PengRobinsonComponentV1, PengRobinsonMixtureV1, PengRobinsonPhaseV1},
    required_nasa9_parameter_dims,
};
use std::collections::BTreeMap;

fn calorics(references: [f64; 2]) -> Vec<Nasa9StandardStateModelV1> {
    ["A", "B"]
        .into_iter()
        .enumerate()
        .map(|(i, name)| {
            let coefficients = [
                0.,
                0.,
                [3.5, 2.5][i],
                [0.001, 0.][i],
                0.,
                0.,
                0.,
                references[i] / R,
                0.,
            ];
            let mut parameters: BTreeMap<_, _> = required_nasa9_parameter_dims()
                .into_iter()
                .take(9)
                .zip(coefficients)
                .map(|((name, dims), value)| (name.into(), LawParameter { value, dims }))
                .collect();
            parameters.insert(
                "reference_pressure".into(),
                LawParameter {
                    value: 100000.,
                    dims: Pressure::DIMS,
                },
            );
            let region = Nasa9RegionV1::from_card(ConstitutiveModelCard {
                law: LawId(NASA9_LAW_ID_V1.into()),
                law_version: NASA9_LAW_VERSION_V1,
                parameters,
                state_schema_version: NASA9_STATE_SCHEMA_VERSION_V1,
                initial_state: InitialStatePolicy::ZeroInternalState,
                validity: ValidityDomain::unconstrained().with("T", 200., 1000.),
                sources: vec![],
                provenance: Provenance {
                    source: "manufactured original NASA thermal source".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
            Nasa9StandardStateModelV1::new(
                SpeciesId::new(name).unwrap(),
                MolarMass::new([0.028, 0.004][i]),
                StandardStateConventionV1::new(
                    StandardStatePhaseV1::Gas,
                    ReferenceEquationOfStateV1::IdealGas,
                    Pressure::new(100000.),
                    ElementalReferenceIdV1::new("thermal-reference").unwrap(),
                )
                .unwrap(),
                vec![region],
            )
            .unwrap()
        })
        .collect()
}
fn ms() -> MaxwellStefanV1 {
    MaxwellStefanV1::new(vec![
        (SpeciesId::new("A").unwrap(), MolarMass::new(0.028)),
        (SpeciesId::new("B").unwrap(), MolarMass::new(0.004)),
    ])
    .unwrap()
}
struct Fixture {
    x: [f64; 2],
    t: [f64; 2],
    logf: [f64; 2],
    hr: [[f64; 2]; 2],
    c: f64,
}
fn fixture(xa: [f64; 2], t: [f64; 2]) -> Fixture {
    let component = |m, tc, pc, w| PengRobinsonComponentV1 {
        molar_mass: MolarMass::new(m),
        critical_temperature: Temperature::new(tc),
        critical_pressure: Pressure::new(pc),
        acentric_factor: w,
    };
    let components = vec![
        component(0.028, 190.564, 4.5992e6, 0.01142),
        component(0.004, 150.687, 4.863e6, -0.00219),
    ];
    let mixtures = xa.map(|a| {
        PengRobinsonMixtureV1::new(
            components.clone(),
            vec![a, 1. - a],
            vec![0., 0.03, 0.03, 0.],
        )
        .unwrap()
    });
    let roots = std::array::from_fn::<_, 2, _>(|i| {
        mixtures[i]
            .evaluate_pt(Temperature::new(t[i]), Pressure::new(4e6))
            .unwrap()
    });
    let roots = roots
        .each_ref()
        .map(|s| s.select(PengRobinsonPhaseV1::Vapor));
    let partial = std::array::from_fn::<_, 2, _>(|i| {
        mixtures[i]
            .evaluate_density_partial_calorics(
                Temperature::new(t[i]),
                Density::new(mixtures[i].molar_mass().value() / roots[i].molar_volume_m3_per_mol()),
            )
            .unwrap()
    });
    Fixture {
        x: [0.5 * xa[0] + 0.5 * xa[1], 1. - (0.5 * xa[0] + 0.5 * xa[1])],
        t,
        logf: std::array::from_fn(|i| {
            fs_math::det::ln(mixtures[1].mole_fractions()[i] / mixtures[0].mole_fractions()[i])
                + roots[1].log_fugacity_coefficients()[i]
                - roots[0].log_fugacity_coefficients()[i]
        }),
        hr: std::array::from_fn(|s| {
            std::array::from_fn(|i| partial[s].components()[i].residual_enthalpy_j_mol())
        }),
        c: 0.5 / roots[0].molar_volume_m3_per_mol() + 0.5 / roots[1].molar_volume_m3_per_mol(),
    }
}
fn call<'a>(
    f: &Fixture,
    sources: &[Nasa9StandardStateModelV1],
    h: [f64; 2],
    k: f64,
    work: &'a mut MaxwellStefanWorkspaceV1,
) -> Evaluation<'a> {
    evaluate(
        &ms(),
        sources,
        &f.x,
        f.t,
        &f.logf,
        [&f.hr[0], &f.hr[1]],
        f.c,
        &[0., 3e-8, 3e-8, 0.],
        &h,
        k,
        0.02,
        MaxwellStefanSolveV1::default(),
        1e-10,
        work,
        || false,
    )
    .unwrap()
}
fn close(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol, "{a:e} != {b:e} (tol {tol:e})");
}

#[test]
fn g0_equal_temperature_matches_original_pr_diffusion_with_centered_energy() {
    let f = fixture([0.75, 0.35], [330., 330.]);
    let sources = calorics([0., 0.]);
    let m = ms();
    let mut a = m.workspace().unwrap();
    let mut b = m.workspace().unwrap();
    let thermal = call(&f, &sources, [0., 0.], 0., &mut a);
    let g = f.logf.map(|v| v / 0.02);
    let ordinary = force::evaluate(
        &m,
        &f.x,
        &g,
        f.c,
        &[0., 3e-8, 3e-8, 0.],
        MaxwellStefanSolveV1::default(),
        1e-10,
        &mut b,
        || false,
    )
    .unwrap();
    for i in 0..2 {
        close(
            thermal.kinetic.diffusion.molar_fluxes()[i],
            ordinary.diffusion.molar_fluxes()[i],
            1e-17,
        );
    }
    close(
        thermal.entropy_from_sources,
        ordinary.force_entropy * 0.02,
        1e-13,
    );
    assert_eq!(thermal.fourier, 0.);
    assert_eq!(thermal.reciprocal, 0.);
    // At equal T, source-relative centered PR enthalpy can still differ from
    // existing upwind enthalpy when composition differs. It is not hidden heat.
    let expected: f64 = (0..2)
        .map(|i| {
            let n = thermal.kinetic.diffusion.molar_fluxes()[i];
            n * (0.5 * f.hr[0][i] + 0.5 * f.hr[1][i] - f.hr[usize::from(n < 0.)][i])
        })
        .sum();
    close(thermal.correction, expected, 1e-14);
    assert!(thermal.correction.abs() > 1e-8);
}

#[test]
fn g0_actual_nonideal_chemical_and_donor_energy_power_closes_positive_entropy() {
    let f = fixture([0.75, 0.35], [310., 350.]);
    let sources = calorics([0., 0.]);
    let m = ms();
    let mut work = m.workspace().unwrap();
    let value = call(&f, &sources, [6000., -9000.], 0.08, &mut work);
    let beta = 1. / f.t[1] - 1. / f.t[0];
    let mut energy = value.additional;
    let mut chemical = 0.;
    for i in 0..2 {
        let n = value.kinetic.diffusion.molar_fluxes()[i];
        let side = usize::from(n < 0.);
        let donor = sources[i]
            .evaluate_properties(Temperature::new(f.t[side]))
            .unwrap()
            .h()
            .value()
            + f.hr[side][i];
        let a = sources[i]
            .evaluate_properties(Temperature::new(f.t[0]))
            .unwrap();
        let b = sources[i]
            .evaluate_properties(Temperature::new(f.t[1]))
            .unwrap();
        let delta = b.h().value() / f.t[1] - b.s().value() - a.h().value() / f.t[0]
            + a.s().value()
            + R * f.logf[i];
        energy += n * donor;
        chemical += n * delta;
    }
    let independently_reconstructed = energy * beta - chemical;
    close(
        independently_reconstructed,
        value.entropy_from_sources,
        2e-12,
    );
    close(
        value.entropy_from_sources,
        value.entropy_from_dissipation,
        2e-12,
    );
    assert!(value.entropy_from_sources > 0.);
    assert!((f.c / (4e6 / (R * 330.)) - 1.).abs() > 0.03);
}

#[test]
fn g0_soret_dufour_cross_response_is_reciprocal_on_the_mass_subspace() {
    let sources = calorics([0., 0.]);
    let m = ms();
    let mut a = m.workspace().unwrap();
    let mut b = m.workspace().unwrap();
    let mut chemical = fixture([0.6, 0.6], [330., 330.]);
    // Independent constitutive coordinate: a species chemical effort. Source
    // binding is tested separately; no altered fixture claims a measured gas.
    chemical.logf = [0.01, 0.];
    let c = call(&chemical, &sources, [6000., -9000.], 0., &mut a);
    let thermal = fixture([0.6, 0.6], [329.999, 330.001]);
    let t = call(&thermal, &sources, [6000., -9000.], 0., &mut b);
    let l_qi = c.reciprocal / (-R * chemical.logf[0] / 0.02);
    let beta = (thermal.t[0] - thermal.t[1]) / (thermal.t[0] * thermal.t[1]);
    let l_iq = t.kinetic.diffusion.molar_fluxes()[0] / (beta / 0.02);
    close(l_qi, l_iq, 2e-7 * l_qi.abs());
    assert!(t.entropy_from_sources >= 0. && c.entropy_from_sources >= 0.);
    close(t.kinetic.diffusion.mass_fluxes().iter().sum(), 0., 1e-18);
}

#[test]
fn g1_zero_transport_heat_temperature_only_force_refines_to_zero() {
    let sources = calorics([0., 0.]);
    let m = ms();
    let mut work = m.workspace().unwrap();
    let mut currents = Vec::new();
    for delta in [8., 4., 2.] {
        let f = fixture([0.6, 0.6], [330. - delta, 330. + delta]);
        let v = call(&f, &sources, [0., 0.], 0., &mut work);
        currents.push(v.kinetic.diffusion.molar_fluxes()[0].abs());
    }
    // Symmetric residual h* is a finite-face quadrature: it removes the linear
    // thermal derivative, with a cubic remaining force at constant p and x.
    assert!(currents[1] < 0.2 * currents[0]);
    assert!(currents[2] < 0.2 * currents[1]);
}

#[test]
fn g3_formation_and_common_heat_datums_leave_physical_fluxes_unchanged() {
    let f = fixture([0.75, 0.35], [310., 350.]);
    let m = ms();
    let mut a = m.workspace().unwrap();
    let mut b = m.workspace().unwrap();
    let base = call(&f, &calorics([0., 0.]), [8192., -4096.], 0.08, &mut a);
    let common = 2.0_f64.powi(40);
    let shifted = call(
        &f,
        &calorics([1e19, -2e19]),
        [common + 8192., common - 4096.],
        0.08,
        &mut b,
    );
    assert_eq!(
        base.kinetic.diffusion.molar_fluxes(),
        shifted.kinetic.diffusion.molar_fluxes()
    );
    assert_eq!(base.centered_heats, shifted.centered_heats);
    assert_eq!(base.correction, shifted.correction);
    assert_eq!(base.reciprocal, shifted.reciprocal);
    assert_eq!(base.entropy_from_sources, shifted.entropy_from_sources);
}

#[test]
fn g4_source_and_cancellation_refusals_return_no_constitutive_result() {
    let f = fixture([0.75, 0.35], [310., 350.]);
    let m = ms();
    let sources = calorics([0., 0.]);
    let mut work = m.workspace().unwrap();
    assert!(matches!(
        evaluate(
            &m,
            &sources,
            &f.x,
            f.t,
            &f.logf,
            [&f.hr[0], &f.hr[1]],
            f.c,
            &[0., 3e-8, 3e-8, 0.],
            &[0., 0.],
            -0.1,
            0.02,
            MaxwellStefanSolveV1::default(),
            1e-10,
            &mut work,
            || false
        ),
        Err(Error::Invalid(_))
    ));
    assert!(matches!(
        evaluate(
            &m,
            &sources,
            &f.x,
            f.t,
            &f.logf,
            [&f.hr[0], &f.hr[1]],
            f.c,
            &[0., 3e-8, 3e-8, 0.],
            &[0., 0.],
            0.1,
            0.02,
            MaxwellStefanSolveV1::default(),
            1e-10,
            &mut work,
            || true
        ),
        Err(Error::Cancelled)
    ));
}
