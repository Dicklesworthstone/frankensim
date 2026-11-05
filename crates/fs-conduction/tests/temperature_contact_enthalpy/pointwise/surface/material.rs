//! Source-surface contact with evolving nodal emissivity and exact P1 moments.
use super::*;
use fs_conduction::{
    MaterialAmbientRadiationPatch, MaterialEnclosureRadiationPatch,
    MaterialPointwiseAmbientRadiation, MaterialPointwiseEnclosureRadiation,
};
use fs_matdb::MaterialCard;
use fs_material::state_point::MaterialPropertySelection;
use fs_material::thermal::source::ThermalMaterialPoint;

const STEP_S: f64 = 0.025;
const FRACTIONS: [f64; 8] = [0.2, 1., 0., 0., 0., 0.5, 1., 0.];
const EMISSIVITY: [[f64; 2]; 2] = [[0.2, 0.8], [0.25, 0.7]];

fn material_card(body: usize) -> MaterialCard {
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS),
        value: PropertyValue::Curve {
            abscissa: "melt".into(), abscissa_dims: fs_qty::Dims::NONE,
            knots: vec![(0., EMISSIVITY[body][0]), (1., EMISSIVITY[body][1])],
            dims: EMISSIVITY_DIMS,
        },
        validity: ValidityDomain::unconstrained().with("T", 280., 380.).with("melt", 0., 1.),
        uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::LinearInside,
        observations: vec![], provenance: Provenance {
            source: "synthetic source-surface contact and phase-emissivity reference".into(),
            license: "internal-test-use".into(), artifact: None,
        },
    }).unwrap();
    MaterialCard::assemble(MaterialStateId {
        chemistry: ["A", "B"][body].into(), phase: "solid".into(),
        process: "fixture".into(), revision: 0,
    }, claims, vec![]).unwrap()
}

fn epsilon(body: usize, fraction: f64) -> f64 {
    EMISSIVITY[body][0] + fraction * (EMISSIVITY[body][1] - EMISSIVITY[body][0])
}

// Exact simplex integral of N_i epsilon_h T_h^4 on an area-1/2 triangle.
// The multinomial expansion and factorial moments are independent of the
// production degree-six quadrature, including correlated T/epsilon variation.
fn moments(face: [usize; 3], body: usize) -> ([f64; 3], [f64; 3], f64) {
    let temperature = face.map(|v| TARGET_T[v]);
    let emissivity = face.map(|v| epsilon(body, FRACTIONS[v]));
    let mut emission = [0.; 3];
    for a in 0..=4 {
        for b in 0..=4 - a {
            let powers = [a, b, 4 - a - b];
            let monomial = (0..3).map(|i| temperature[i].powi(powers[i])).product::<f64>();
            for i in 0..3 {
                for j in 0..3 {
                    let factor = f64::from((powers[i] + 1)
                        * (powers[j] + 1 + if i == j { 1 } else { 0 }));
                    emission[i] += 0.5 / 840. * factor * emissivity[j] * monomial;
                }
            }
        }
    }
    let sum = emissivity.iter().sum::<f64>();
    let absorption = emissivity.map(|value| 0.5 / 12. * (sum + value));
    (emission.map(|value| value * STEFAN_BOLTZMANN_W_M2_K4), absorption, sum / 3.)
}

// Independent two-equation gray radiosity, with an explicit black environment
// only for mode 2. Mode 0 is ambient; mode 1 is closed internal exchange.
fn radiation_load(mode: u8, bath_k: f64) -> [f64; 8] {
    let cold = moments(COLD_FACE, 0);
    let hot = moments(HOT_FACE, 1);
    let black = STEFAN_BOLTZMANN_W_M2_K4 * bath_k.powi(4);
    let irradiation = if mode == 0 {
        [black, 0.]
    } else {
        let exchange = if mode == 1 { 1. } else { 0.4 };
        let escape = 1. - exchange;
        let emitted = [cold.0.iter().sum::<f64>() / 0.5, hot.0.iter().sum::<f64>() / 0.5];
        let reflected = [1. - cold.2, 1. - hot.2];
        let rhs = [emitted[0] + reflected[0] * escape * black,
            emitted[1] + reflected[1] * escape * black];
        let j0 = (rhs[0] + reflected[0] * exchange * rhs[1])
            / (1. - reflected[0] * reflected[1] * exchange * exchange);
        let j1 = rhs[1] + reflected[1] * exchange * j0;
        [exchange * j1 + escape * black, exchange * j0 + escape * black]
    };
    let mut load = [0.; 8];
    for (body, face, (emission, absorption, _)) in [(0, COLD_FACE, cold), (1, HOT_FACE, hot)] {
        if body == 1 && mode == 0 { continue; }
        for i in 0..3 { load[face[i]] = emission[i] - irradiation[body] * absorption[i]; }
    }
    load
}

fn checked_advance(
    cx: &Cx<'_>, old: &[f64],
    advance: impl Fn(&Cx<'_>, ContactIterationConfig) -> Result<ContactEnthalpyStep, EnthalpyError>,
) -> ContactEnthalpyStep {
    let before = old.to_vec();
    assert!(advance(cx, ContactIterationConfig { max_iterations: 1, ..coupling() }).is_err());
    let gate = CancelGate::new_clock_free();
    gate.request();
    context(&gate, |cancelled| assert!(matches!(advance(cancelled, coupling()),
        Err(EnthalpyError::Conduction(ConductionError::Cancelled { .. })))));
    assert_eq!(old, before.as_slice());
    let next = advance(cx, coupling()).unwrap();
    let repeat = advance(cx, coupling()).unwrap();
    assert_eq!(old, before.as_slice());
    assert_eq!(next.conduction.specific_enthalpy_j_kg, repeat.conduction.specific_enthalpy_j_kg);
    assert_eq!(next.physical_contacts.states(), repeat.physical_contacts.states());
    assert_eq!(next.conduction.pointwise_radiation, repeat.conduction.pointwise_radiation);
    assert_eq!(next.conduction.pointwise_enclosure_radiation, repeat.conduction.pointwise_enclosure_radiation);
    next
}

fn case(cx: &Cx<'_>, mode: u8, bath_k: f64) {
    let mesh = mesh();
    let boundary = if mode == 0 { ambient_boundary(&mesh) } else { enclosure_boundary(&mesh) };
    let contact_card = surface_card(280., 380.);
    let contacts = TemperatureContacts::new(cx, &mesh, &boundary, vec![],
        vec![law(&mesh, &contact_card)], 1).unwrap();
    let cards = [material_card(0), material_card(1)];
    let curves = [0, 1].map(|body| EquilibriumEnthalpyPhaseCurve::try_new(
        cards[body].content_hash(), chart(body == 1).knots().to_vec(),
    ).unwrap());
    let engine = engine(cx, &mesh, &curves[0], &curves[1]);
    let material = ConductivityModel::isotropic_declared(2.).unwrap();
    let source = ScalarField::Uniform(0.);
    let problem = ConductionProblem { mesh: &mesh, boundary: &boundary,
        material: &material, element_materials: None, source: &source };
    let point = QueryPoint::new().with("T", 320.).unwrap().with("melt", 0.5).unwrap();
    let binding = |body: usize| ThermalMaterialPoint { material: &cards[body], context: &point,
        temperature_axis: "T", liquid_fraction_axis: Some("melt") };
    let radiative = radiation_load(mode, bath_k);
    let convection = if mode == 0 { convection_load() } else { [0.; 8] };
    let exterior = std::array::from_fn(|i| radiative[i] + convection[i]);
    // Reuse the independent full stiffness/contact moments, changing only
    // the manufactured step length. No production solve supplies this history.
    let old: Vec<_> = manufactured_history(exterior).iter().zip(TARGET_H)
        .map(|(&h, target)| target + (STEP_S / DT) * (h - target)).collect();
    let before = old.clone();
    let accepted = if mode == 0 {
        let patches = [MaterialAmbientRadiationPatch {
            region: "exterior", ambient_temperature_k: bath_k, curve: &curves[0],
            material: binding(0), selection: MaterialPropertySelection::SingleClaimOnly,
        }];
        let radiation = MaterialPointwiseAmbientRadiation::new(
            cx, &mesh, &boundary, &patches,
        ).unwrap();
        let copied_curve = curves[0].clone();
        let foreign_patches = [MaterialAmbientRadiationPatch {
            region: "exterior", ambient_temperature_k: bath_k, curve: &copied_curve,
            material: binding(0), selection: MaterialPropertySelection::SingleClaimOnly,
        }];
        let foreign = MaterialPointwiseAmbientRadiation::new(
            cx, &mesh, &boundary, &foreign_patches,
        ).unwrap();
        assert!(contacts.advance_enthalpy_with_material_ambient_radiation(
            cx, &engine, problem, &old, STEP_S, config(), coupling(), &foreign,
        ).is_err(), "an equal but different caloric object cannot replace the storage chart");
        checked_advance(cx, &old, |cx, coupling|
            contacts.advance_enthalpy_with_material_ambient_radiation(
                cx, &engine, problem, &old, STEP_S, config(), coupling, &radiation,
            ))
    } else {
        let surfaces = [0_usize, 1].map(|body| {
            let name = ["left", "right"][body];
            RadiationSurface::new(&mesh, name, |face| face.vertices.iter().all(|v|
                [4 * body, 4 * body + 1, 4 * body + 3].contains(&(*v as usize))),
                SurfaceEmissivity::from_condition(name, binding(body),
                    MaterialPropertySelection::SingleClaimOnly).unwrap().0,
            ).unwrap()
        }).to_vec();
        let closed = GrayDiffuseEnclosure::new(surfaces.clone(),
            ViewFactorMatrix::infinite_parallel_plates(0.5).unwrap()).unwrap();
        let open = OpenDiffuseEnclosure::new(cx, surfaces, vec![vec![0., 0.4], vec![0.4, 0.]],
            ViewFactorEvidence::Analytic { geometry: "manufactured reciprocal material exchange".into() },
            ViewFactorTolerance::default()).unwrap();
        let base = if mode == 1 {
            PointwiseEnclosureRadiation::new(cx, &mesh, &boundary, &closed, PointwiseEnclosureTolerance::default())
        } else {
            PointwiseEnclosureRadiation::new_open(cx, &mesh, &boundary, &open, bath_k,
                PointwiseEnclosureTolerance::default())
        }.unwrap();
        let radiation = MaterialPointwiseEnclosureRadiation::new(cx, &base, &[0_usize, 1].map(|body|
            MaterialEnclosureRadiationPatch { region: ["left", "right"][body],
                curve: &curves[body], material: binding(body),
                selection: MaterialPropertySelection::SingleClaimOnly },
        )).unwrap();
        checked_advance(cx, &old, |cx, coupling|
            contacts.advance_enthalpy_with_material_enclosure_radiation(
                cx, &engine, problem, &old, STEP_S, config(), coupling, &radiation,
            ))
    };
    assert_eq!(old, before);
    audit(&accepted, &TARGET_H, &TARGET_T);
    let heat = radiative.iter().sum::<f64>();
    close(accepted.physical.radiation_out_w, heat, 2e-6);
    close(accepted.physical.convective_out_w, convection.iter().sum(), 2e-6);
    close(accepted.physical.stored_energy_change_j, -STEP_S * exterior.iter().sum::<f64>(), 2e-8);
    let stored: f64 = accepted.conduction.specific_enthalpy_j_kg.iter().zip(&old).enumerate()
        .map(|(i, (&h, &previous))| (h - previous) / if i < 4 { 24. } else { 12. }).sum();
    close(stored, -STEP_S * exterior.iter().sum::<f64>(), 2e-8);
    if mode == 1 { close(heat, 0., 2e-12); }
    else { assert!(heat.abs() > 0.1, "actual external radiation must exchange stored energy"); }
    let reports = if mode == 0 {
        accepted.conduction.pointwise_radiation.as_ref().unwrap().patches.iter()
            .map(|row| &row.material_states).collect::<Vec<_>>()
    } else {
        let report = accepted.conduction.pointwise_enclosure_radiation.as_ref().unwrap();
        if mode == 2 {
            close(report.ambient.as_ref().unwrap().outward_heat_w, -heat, 2e-6);
        } else { assert!(report.ambient.is_none()); }
        report.patches.iter().map(|row| &row.material_states).collect()
    };
    for (body, states) in reports.iter().enumerate() {
        assert_eq!(states.len(), 3);
        let snapshot = epsilon(body, 0.5);
        assert!(states.iter().any(|state| (state.emissivity.value() - snapshot).abs() > 0.05));
        for state in states.iter() {
            let vertex = state.vertex;
            assert_eq!(state.emissivity.card_identity(), cards[body].content_hash());
            assert_eq!(state.phase.material_card_identity(), cards[body].content_hash());
            assert_eq!(state.phase.specific_enthalpy_j_kg(), accepted.conduction.specific_enthalpy_j_kg[vertex]);
            close(state.phase.liquid_mass_fraction(), FRACTIONS[vertex], 4e-8);
            close(state.emissivity.value(), epsilon(body, FRACTIONS[vertex]), 4e-8);
            cards[body].claims().verify_receipt(state.emissivity.receipt()).unwrap();
        }
    }
}

#[test]
fn g1_surface_contact_material_ambient_recovers_latent_states_and_exact_local_emission() {
    context(&CancelGate::new_clock_free(), |cx| {
        for bath_k in [300., 340.] { case(cx, 0, bath_k); }
    });
}

#[test]
fn g1_surface_contact_material_enclosure_recovers_correlated_emission_and_reflected_heat() {
    context(&CancelGate::new_clock_free(), |cx| {
        for mode in [1, 2] { case(cx, mode, 320.); }
    });
}
