//! Manufactured positive table with a unique, stiff nonlinear operating point.
//! This is an analytical solver fixture, not measured material data.
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    conductor::{
        ElectricalPropertyBasis, ResolvedConductorState,
        response::{ConductorMagnitude, ScalarConductorResponse},
    },
    state_point::{MaterialPropertySelection, resolve_material_state_point},
};

pub fn law(current_control: bool) -> ScalarConductorResponse {
    let (basis, coordinate) = if current_control {
        (
            ElectricalPropertyBasis::Conductivity,
            ConductorMagnitude::ElectricField,
        )
    } else {
        (
            ElectricalPropertyBasis::Resistivity,
            ConductorMagnitude::CurrentDensity,
        )
    };
    let key = PropertyKey::new(basis.property_name(), basis.dimensions());
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: key.clone(),
            value: PropertyValue::Curve {
                abscissa: "operating".into(),
                abscissa_dims: coordinate.quantity().dims(),
                knots: vec![(0.0, 1.0), (0.9, 1.0), (1.1, 9.0), (10.0, 9.0)],
                dims: basis.dimensions(),
            },
            validity: ValidityDomain::unconstrained().with_quantity(
                "operating",
                coordinate.quantity(),
                0.0,
                10.0,
            ),
            interpolation: InterpolationPolicy::LinearInside,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured stiff passive conductor; not measured material".into(),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured".into(),
            phase: "solid".into(),
            process: "piecewise positive source response".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let query = QueryPoint::new()
        .with_quantity("operating", coordinate.quantity(), 0.0)
        .unwrap();
    let material = resolve_material_state_point(
        &card,
        &query,
        &[basis.requirement_for(&key).unwrap()],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let original = ResolvedConductorState::try_from_material_with_basis(material, basis).unwrap();
    ScalarConductorResponse::new(&card, original, &[("operating", coordinate)]).unwrap()
}
