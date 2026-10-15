//! The solver-facing and composition-facing APIs share one nominal owner.
use fs_couple::{CoordinateBinding, PortKind, PortOrientation, PortTimestamp};

#[test]
fn electrical_schema_crosses_the_boundary_without_conversion() {
    let stable = |s| fs_ports::StableId::new(s).unwrap();
    let schema: fs_ports::PortSchema = PortKind::Electrical.scalar_seed_schema(
        stable("terminal/a"),
        CoordinateBinding::new(stable("basis/scalar"), stable("frame/local"),
            PortOrientation::OutwardFromOwner),
        PortTimestamp::new(stable("clock/main"), 0),
    ).unwrap();
    let composed: fs_couple::PortSchema = schema.clone();
    let direct: fs_ports::PortSchema = composed;
    assert_eq!(direct, schema);
    assert_eq!(fs_couple::PORT_SCHEMA_VERSION, fs_ports::PORT_SCHEMA_VERSION);
}

#[test]
fn errors_and_all_retained_core_types_have_one_identity() {
    fn same<T: 'static, U: 'static>() {
        assert_eq!(std::any::TypeId::of::<T>(), std::any::TypeId::of::<U>());
    }
    same::<fs_ports::CoupleError, fs_couple::CoupleError>();
    same::<fs_ports::StreamPort, fs_couple::StreamPort>();
    same::<fs_ports::ConservativeJunction, fs_couple::ConservativeJunction>();
    same::<fs_ports::WindowAuditReport, fs_couple::WindowAuditReport>();
    let error: fs_couple::CoupleError = fs_ports::StableId::new("").unwrap_err();
    assert_eq!(error, fs_couple::StableId::new("").unwrap_err());
}
