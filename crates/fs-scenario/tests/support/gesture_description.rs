//! A declared address fixture; no geometry, material or realism assertion.
use fs_qty::{Dims, QtyAny};
use fs_scenario::assembly::*;
use fs_scenario::gesture::{GestureTarget as T, GestureValue as V};
use fs_scenario::{EntityCatalog, EntityDeclaration, Environment, Scenario};

pub fn fixture() -> (AssemblyDescription, Vec<(T, ParameterPath, V)>) {
    let s = Scenario::new(
        "test-fixture:typed-gesture-addresses",
        1,
        Environment::earth_lab(),
    );
    let mut catalog = EntityCatalog::new();
    let root = catalog
        .declare(EntityDeclaration::assembly("test-fixture:root"))
        .unwrap();
    let mut part = |name| {
        catalog
            .declare(EntityDeclaration::part(root, name))
            .unwrap()
    };
    let holes = [part("hole-zero"), part("hole-one")];
    let string = part("string-zero");
    let valve = part("valve-zero");
    let length = Dims([1, 0, 0, 0, 0, 0]);
    let pressure = Dims([-1, 1, -2, 0, 0, 0]);
    let force = Dims([1, 1, -2, 0, 0, 0]);
    let speed = Dims([1, 0, -1, 0, 0, 0]);
    let specs = vec![
        (
            T::ValveInsertedLength { valve: 0 },
            valve,
            "inserted_length",
            length,
            V::LengthM(0.1),
        ),
        (
            T::ToneHoleSigma { hole: 0 },
            holes[0],
            "open_fraction",
            Dims::NONE,
            V::Fraction(0.5),
        ),
        (
            T::ToneHoleSigma { hole: 1 },
            holes[1],
            "open_fraction",
            Dims::NONE,
            V::Fraction(0.75),
        ),
        (
            T::SlideLength,
            root,
            "slide_length",
            length,
            V::LengthM(0.2),
        ),
        (
            T::StringTension { string: 0 },
            string,
            "tension",
            force,
            V::TensionN(20.0),
        ),
        (
            T::Fret { string: 0 },
            string,
            "fret",
            Dims::NONE,
            V::Fret {
                engaged: true,
                position_m: 0.1,
                height_m: 0.001,
                velocity_m_per_s: 0.2,
            },
        ),
        (
            T::BowStroke { string: 0 },
            string,
            "bow",
            Dims::NONE,
            V::Bow {
                velocity_m_per_s: 0.1,
                normal_force_n: 0.01,
                station: 0.2,
            },
        ),
        (
            T::HammerStrike { string: 0 },
            string,
            "strike",
            speed,
            V::StrikeVelocity {
                velocity_m_per_s: 0.2,
            },
        ),
        (
            T::SustainPedal,
            root,
            "sustain",
            Dims::NONE,
            V::Fraction(0.5),
        ),
        (
            T::UnaCordaPedal,
            root,
            "una_corda",
            Dims::NONE,
            V::Fraction(0.5),
        ),
        (
            T::SostenutoHold { string: 0 },
            string,
            "sostenuto_hold",
            Dims::NONE,
            V::Fraction(1.0),
        ),
        (
            T::TerminationSwap,
            root,
            "termination",
            Dims::NONE,
            V::TerminationSwap {
                digest_hex: "ab".repeat(32),
                fade_s: 0.01,
            },
        ),
        (
            T::BlowingPressure,
            root,
            "blowing_pressure",
            pressure,
            V::PressurePa(2000.0),
        ),
        (
            T::PreStress,
            root,
            "prestress",
            Dims::NONE,
            V::Multiplier(1.1),
        ),
        (
            T::RestAperture,
            root,
            "rest_aperture",
            length,
            V::LengthM(0.001),
        ),
        (T::JetSpeed, root, "jet_speed", speed, V::VelocityMPerS(3.0)),
        (
            T::JetAngle,
            root,
            "jet_angle",
            Dims::NONE,
            V::AngleRad(-0.2),
        ),
    ];
    let cases: Vec<_> = specs
        .into_iter()
        .map(|(t, entity, name, dims, v)| {
            (
                t,
                ParameterPath {
                    entity,
                    parameter: ParameterKey::new(name, dims),
                },
                v,
            )
        })
        .collect();
    let assignments = cases
        .iter()
        .map(|(_, p, _)| (p.clone(), QtyAny::new(0.0, p.parameter.dims)))
        .collect();
    let parts = [holes[0], holes[1], string, valve]
        .into_iter()
        .zip([
            PartKind::Duct,
            PartKind::Duct,
            PartKind::String,
            PartKind::ValvePattern,
        ])
        .map(|(entity, kind)| PartDecl {
            entity,
            kind,
            geometry: GeometrySource::Parametric {
                generator: "test-fixture:address-only".into(),
                params: vec![],
            },
            materials: vec![],
            uncertainty: vec![],
        })
        .collect();
    let description = AssemblyDescription::new(
        &s,
        AssemblyInput {
            root,
            catalog,
            parts,
            interfaces: vec![],
            configurations: vec![ConfigurationState {
                name: "declared-inputs".into(),
                assignments,
                provenance: StateProvenance::Declared {
                    rationale: "address fixture only, no physical specimen".into(),
                },
            }],
            observers: vec![],
            reference_state: Some(ReferenceState::Assembled),
        },
        DEFAULT_ASSEMBLY_BUDGET,
    )
    .unwrap();
    (description, cases)
}
