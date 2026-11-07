//! Manufactured ordered-interface spectra shared by native consumer tests.
//! These values verify mechanics and source continuity; they are not measured data.
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, SurfaceSpec, SystemContext, UncertaintyModel,
};
use fs_material::cohesive::relaxation::{ViscoelasticBondLimits, source::*};
use fs_qty::{Dims, Length, Pressure, Temperature, Time};

pub fn limits() -> ViscoelasticBondLimits {
    ViscoelasticBondLimits {
        maximum_jump_m: 1e-4,
        maximum_rate_m_s: 1.0,
        maximum_elapsed_s: 2.0,
    }
}

pub fn card(missing_shear: bool, wrong_units: bool) -> InterfaceSystemCard {
    build(missing_shear, wrong_units, false, 0.001)
}

/// The same indexed springs soften and relax faster between 280 K and 320 K.
pub fn temperature_card() -> InterfaceSystemCard {
    build(false, false, true, 0.001)
}

/// Deliberately inadmissible geometry change for held-jump source-continuity tests.
pub fn expanding_temperature_card() -> InterfaceSystemCard {
    build(false, false, true, 0.0011)
}

fn build(
    missing_shear: bool,
    wrong_units: bool,
    thermal: bool,
    hot_thickness: f64,
) -> InterfaceSystemCard {
    let mut claims = ClaimSet::new();
    let mut parameters = vec![(
        BONDLINE_THICKNESS_PROPERTY.to_owned(),
        Length::DIMS,
        0.001,
        hot_thickness,
    )];
    for (prefix, e0, e1, tau, hot0, hot1, hot_tau) in [
        (NORMAL_RELAXATION_PREFIX, 1e6, 2e6, 0.01, 5e5, 1.2e6, 0.002),
        (SHEAR_RELAXATION_PREFIX, 4e5, 8e5, 0.025, 3e5, 4e5, 0.005),
    ] {
        if missing_shear && prefix == SHEAR_RELAXATION_PREFIX {
            continue;
        }
        parameters.extend([
            (format!("{prefix}_branch_count"), Dims::NONE, 1.0, 1.0),
            (
                format!("{prefix}_equilibrium_modulus"),
                Pressure::DIMS,
                e0,
                hot0,
            ),
            (format!("{prefix}_modulus_1"), Pressure::DIMS, e1, hot1),
            (
                format!("{prefix}_time_1"),
                if wrong_units {
                    Length::DIMS
                } else {
                    Time::DIMS
                },
                tau,
                hot_tau,
            ),
        ]);
    }
    for (name, dims, cold, hot) in parameters {
        let changing = thermal && cold != hot;
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: if changing {
                    PropertyValue::Curve {
                        abscissa: "T".into(),
                        abscissa_dims: Temperature::DIMS,
                        knots: vec![(280.0, cold), (320.0, hot)],
                        dims,
                    }
                } else {
                    PropertyValue::Scalar { value: cold, dims }
                },
                validity: ValidityDomain::unconstrained().with("T", 280.0, 320.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: if changing {
                    InterpolationPolicy::LinearInside
                } else {
                    InterpolationPolicy::ConstantWithinValidity
                },
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured intact joint relaxation; not measured adhesive data"
                        .into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let surface = |name: &str| SurfaceSpec {
        material: MaterialStateId {
            chemistry: name.into(),
            phase: "solid".into(),
            process: "specified prepared substrate".into(),
            revision: 1,
        },
        texture_frame: name.into(),
    };
    InterfaceSystemCard::assemble(
        surface("A"),
        surface("B"),
        SystemContext {
            medium: "fixed intact thin adhesive".into(),
            third_body: None,
            environment: if thermal {
                "bounded temperature"
            } else {
                "fixed temperature"
            }
            .into(),
            history: "declared cure and normal restraint".into(),
        },
        claims,
        vec![],
    )
    .unwrap()
}
