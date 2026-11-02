//! Manufactured temperature-dependent original cards shared by source/native tests.
#![allow(dead_code)]

use fs_contact::interface_binding::adhesive_shear::{
    KINETIC_JUNCTION_SHEAR_PROPERTY, STATIC_JUNCTION_SHEAR_PROPERTY,
};
use fs_contact::interface_binding::*;
use fs_contact::normal_patch::adhesive::{JkrEnvelope, SourceJkrSphereContact};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialCard, MaterialStateId,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, SurfaceSpec, SystemContext,
    UncertaintyModel,
};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_interface_state_point, resolve_isotropic_solid_state_point,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, Pressure, QuantitySpec, Temperature};
use fs_tribo::InputAuthority;

pub fn point() -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "T",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            300.0,
        )
        .unwrap()
        .with_quantity(
            "p",
            QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static)),
            2e5,
        )
        .unwrap()
}

fn scalar(name: &str, value: f64, dims: Dims) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Scalar { value, dims },
        validity: ValidityDomain::unconstrained()
            .with("T", 280.0, 320.0)
            .with("p", 1e5, 3e5),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: vec![],
        provenance: Provenance {
            source: format!("manufactured thermal adhesion {name}; not measured material data"),
            license: "CC0-1.0".into(),
            artifact: None,
        },
    }
}

fn curve(name: &str, dims: Dims, values: [f64; 2], maximum_t: f64) -> PropertyClaim {
    let mut claim = scalar(name, values[0], dims);
    claim.value = PropertyValue::Curve {
        abscissa: "T".into(),
        abscissa_dims: Temperature::DIMS,
        knots: vec![(280.0, values[0]), (320.0, values[1])],
        dims,
    };
    claim.validity = claim.validity.with("T", 280.0, maximum_t);
    claim.interpolation = InterpolationPolicy::LinearInside;
    claim
}

fn alternative(name: &str, dims: Dims, value: f64) -> PropertyClaim {
    let mut claim = scalar(name, value, dims);
    claim.validity = claim.validity.with("T", 310.0, 320.0);
    claim
        .provenance
        .source
        .push_str("; competing hot-only source");
    claim
}

fn bulk(name: &str, young_300: f64, competing: bool, maximum_t: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for property in [
        scalar("density", 1000.0, Density::DIMS),
        curve(
            "young_modulus",
            Pressure::DIMS,
            [1.2 * young_300, 0.8 * young_300],
            maximum_t,
        ),
        scalar("poisson_ratio", 0.25, Dims::NONE),
        scalar("yield_stress", 1e10, Pressure::DIMS),
    ] {
        claims.insert_claim(property).unwrap();
    }
    if competing {
        claims
            .insert_claim(alternative(
                "young_modulus",
                Pressure::DIMS,
                8.0 * young_300,
            ))
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: name.into(),
            phase: "solid".into(),
            process: "manufactured thermal contact fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

pub struct Fixture {
    pub cards: [MaterialCard; 2],
    pub normal: BoundNormalInterfaceState,
}

impl Fixture {
    pub fn new() -> Self {
        Self::with_options(1.0, 1.0, false, 320.0)
    }

    /// E_a(300)=1e6*scale, E_b(300)=2e6*scale; each changes by -1%/K.
    /// W(300)=work_scale and changes by -3%/K. Optional hot-only alternatives
    /// test claim pinning; maximum_t bounds only the originally selected E/W.
    pub fn with_options(
        modulus_scale: f64,
        work_scale: f64,
        competing: bool,
        maximum_t: f64,
    ) -> Self {
        let cards = [
            bulk(
                "thermal-contact-A",
                1e6 * modulus_scale,
                competing,
                maximum_t,
            ),
            bulk(
                "thermal-contact-B",
                2e6 * modulus_scale,
                competing,
                maximum_t,
            ),
        ];
        let mut claims = ClaimSet::new();
        for property in [
            curve(
                ADHESION_ENERGY_PROPERTY,
                ADHESION_ENERGY_DIMS,
                [1.6 * work_scale, 0.4 * work_scale],
                maximum_t,
            ),
            curve(
                STATIC_JUNCTION_SHEAR_PROPERTY,
                Pressure::DIMS,
                [3600.0, 2400.0],
                320.0,
            ),
            curve(
                KINETIC_JUNCTION_SHEAR_PROPERTY,
                Pressure::DIMS,
                [2400.0, 1600.0],
                320.0,
            ),
        ] {
            claims.insert_claim(property).unwrap();
        }
        if competing {
            claims
                .insert_claim(alternative(
                    ADHESION_ENERGY_PROPERTY,
                    ADHESION_ENERGY_DIMS,
                    8.0 * work_scale,
                ))
                .unwrap();
        }
        let interface = InterfaceSystemCard::assemble(
            SurfaceSpec {
                material: cards[0].id().clone(),
                texture_frame: "A-frame".into(),
            },
            SurfaceSpec {
                material: cards[1].id().clone(),
                texture_frame: "B-frame".into(),
            },
            SystemContext {
                medium: "dry".into(),
                third_body: None,
                environment: "manufactured temperature test".into(),
                history: "initial".into(),
            },
            claims,
            vec![],
        )
        .unwrap();
        let query = point();
        let selection = MaterialPropertySelection::SingleClaimOnly;
        let a = resolve_isotropic_solid_state_point(&cards[0], &query, selection.clone()).unwrap();
        let b = resolve_isotropic_solid_state_point(&cards[1], &query, selection.clone()).unwrap();
        let dry =
            bind_dry_interface_system_card(&interface, InputAuthority::SyntheticFixture).unwrap();
        let elastic = bind_isotropic_elastic_interface(&dry, &a, &b).unwrap();
        let resolved = resolve_interface_state_point(
            &interface,
            &query,
            &[ScalarPropertyRequirement::try_new(
                ADHESION_ENERGY_PROPERTY,
                ADHESION_ENERGY_DIMS,
                ScalarAdmissibility::NonNegative,
            )
            .unwrap()],
            selection,
        )
        .unwrap();
        Self {
            cards,
            normal: bind_normal_interface_state(&elastic, &resolved).unwrap(),
        }
    }

    pub fn jkr(&self) -> SourceJkrSphereContact {
        SourceJkrSphereContact::new(self.normal.clone(), 1.0, 1e-7, envelope()).unwrap()
    }

    pub fn cards(&self) -> [&MaterialCard; 2] {
        [&self.cards[0], &self.cards[1]]
    }
}

pub fn envelope() -> JkrEnvelope {
    JkrEnvelope {
        maximum_indentation_m: 0.001,
        half_space_depth_m: 1.0,
        layer_thickness_m: 1.0,
        maximum_patch_to_radius: 0.1,
        maximum_patch_to_depth: 0.1,
        maximum_patch_to_layer: 0.1,
        maximum_nominal_strain: 0.1,
        maximum_repulsive_pressure_to_yield: 0.1,
        maximum_process_zone_to_patch: 0.01,
        minimum_tabor_parameter: 5.0,
    }
}
