//! Assembly-owned catalog transport. Identity operations are replayed through
//! the existing catalog; every resulting receipt must equal the retained one.
use super::*;
use crate::assembly::wire::{Reader, Writer};
use crate::assembly::{AssemblyBudget, AssemblyError, invalid};

pub(crate) fn write_assembly_catalog(c: &EntityCatalog, w: &mut Writer, b: AssemblyBudget) {
    let q = c.budget;
    for n in [
        q.max_entities,
        q.max_receipts,
        q.max_name_bytes,
        q.max_supersession_hops,
        q.max_hierarchy_depth,
        q.max_datums,
        q.max_tolerances,
        q.max_placements,
        q.max_bindings,
    ] {
        w.integer(n as u64);
    }
    w.count(c.entities.len(), b.max_catalog_entities, "catalog entities");
    for e in &c.entities {
        let d = &e.declaration;
        w.entity(e.id);
        w.maybe_entity(d.parent);
        w.text(&d.declared_name);
        w.text(&d.display_name);
        w.tag(u8::from(d.fingerprint.is_some()));
        if let Some(f) = d.fingerprint {
            w.hash(f.hash());
        }
        w.tag(u8::from(d.pair.is_some()));
        if let Some(p) = d.pair {
            w.tag(p.pairing.tag());
            w.entity(p.from);
            w.entity(p.to);
        }
        w.tag(u8::from(d.legacy));
    }
    w.count(c.receipts.len(), b.max_catalog_receipts, "catalog receipts");
    for receipt in &c.receipts {
        w.integer(receipt.sequence);
        w.tag(receipt.event.tag());
        w.entity(receipt.subject);
        w.maybe_entity(receipt.predecessor);
        w.tag(receipt.basis.tag());
        w.maybe_text(receipt.display_before.as_deref());
        w.maybe_text(receipt.display_after.as_deref());
        w.hash(receipt.previous);
        w.hash(receipt.digest);
    }
    w.hash(c.chain);
    w.count(c.datums.len(), b.max_nested_items, "catalog records");
    for d in &c.datums {
        w.hash(d.id.0);
        w.entity(d.owner);
        w.text(&d.declared_name);
        w.tag(d.feature.tag());
        w.count(d.references.len(), MAX_DATUM_FRAME_LEN, "datum references");
        for id in &d.references {
            w.hash(id.0);
        }
    }
    w.count(c.tolerances.len(), b.max_nested_items, "catalog records");
    for t in &c.tolerances {
        w.entity(t.subject);
        w.tag(match t.kind {
            ToleranceKind::Position => 0,
            ToleranceKind::Flatness => 1,
            ToleranceKind::Parallelism => 2,
            ToleranceKind::Perpendicularity => 3,
            ToleranceKind::Profile => 4,
            ToleranceKind::InterfaceGap => 5,
        });
        w.qty(t.magnitude);
        w.count(t.datum_frame.len(), MAX_DATUM_FRAME_LEN, "datum frame");
        for id in &t.datum_frame {
            w.hash(id.0);
        }
        w.source(&t.source);
    }
    w.count(c.placements.len(), b.max_nested_items, "catalog records");
    for p in &c.placements {
        w.entity(p.occurrence);
        w.integer(u64::from(p.frame.0));
        match p.basis {
            PlacementBasis::Nominal => w.tag(0),
            PlacementBasis::AsBuilt { registration_ref } => {
                w.tag(1);
                w.hash(registration_ref);
            }
        }
    }
}

fn size(r: &mut Reader<'_>) -> Result<usize, AssemblyError> {
    usize::try_from(r.integer()?).map_err(|_| r.fail("catalog limit does not fit host"))
}
fn declaration(r: &mut Reader<'_>) -> Result<(EntityId, EntityDeclaration), AssemblyError> {
    let id = r.entity()?;
    let parent = r.maybe_entity()?;
    let declared_name = r.text()?;
    let display_name = r.text()?;
    let fingerprint = if r.boolean()? {
        Some(GeometryFingerprint::from_hash(r.hash()?))
    } else {
        None
    };
    let pair = if r.boolean()? {
        let pairing = match r.tag()? {
            1 => InterfacePairing::Ordered,
            2 => InterfacePairing::Unordered,
            _ => return Err(r.fail("interface pairing")),
        };
        let from = r.entity()?;
        let to = r.entity()?;
        if pairing == InterfacePairing::Unordered && from > to {
            return Err(r.fail("noncanonical unordered interface"));
        }
        Some(InterfacePair { from, to, pairing })
    } else {
        None
    };
    let legacy = r.boolean()?;
    let d = EntityDeclaration {
        kind: id.kind(),
        parent,
        declared_name,
        display_name,
        fingerprint,
        pair,
        legacy,
    };
    if d.identity() != id {
        return Err(r.fail("entity declaration digest mismatch"));
    }
    Ok((id, d))
}
fn receipt(r: &mut Reader<'_>) -> Result<IdentityReceipt, AssemblyError> {
    let sequence = r.integer()?;
    let event = match r.tag()? {
        1 => RebindEvent::Declaration,
        2 => RebindEvent::Rename,
        3 => RebindEvent::Import,
        4 => RebindEvent::Remesh,
        5 => RebindEvent::RevisionMigration,
        6 => RebindEvent::LegacyMigration,
        7 => RebindEvent::Retirement,
        _ => return Err(r.fail("receipt event")),
    };
    let subject = r.entity()?;
    let predecessor = r.maybe_entity()?;
    let basis = match r.tag()? {
        1 => MatchBasis::NewIdentity,
        2 => MatchBasis::Unchanged,
        3 => MatchBasis::GeometryFingerprint,
        4 => MatchBasis::DeclaredPath,
        5 => MatchBasis::LegacyName,
        6 => MatchBasis::Asserted,
        7 => MatchBasis::Absent,
        _ => return Err(r.fail("receipt basis")),
    };
    Ok(IdentityReceipt {
        sequence,
        event,
        subject,
        predecessor,
        basis,
        display_before: r.maybe_text()?,
        display_after: r.maybe_text()?,
        previous: r.hash()?,
        digest: r.hash()?,
    })
}

pub(crate) fn read_assembly_catalog(r: &mut Reader<'_>) -> Result<EntityCatalog, AssemblyError> {
    let budget = EntityBudget {
        max_entities: size(r)?,
        max_receipts: size(r)?,
        max_name_bytes: size(r)?,
        max_supersession_hops: size(r)?,
        max_hierarchy_depth: size(r)?,
        max_datums: size(r)?,
        max_tolerances: size(r)?,
        max_placements: size(r)?,
        max_bindings: size(r)?,
    };
    let mut catalog = EntityCatalog::with_budget(budget);
    let n = r.count(
        r.budget.max_catalog_entities.min(budget.max_entities),
        "catalog entities",
    )?;
    let mut declarations = Vec::new();
    for _ in 0..n {
        declarations.push(declaration(r)?);
    }
    let n = r.count(
        r.budget.max_catalog_receipts.min(budget.max_receipts),
        "catalog receipts",
    )?;
    let mut next = 0;
    for _ in 0..n {
        let expected = receipt(r)?;
        if !expected.verifies() {
            return Err(r.fail("receipt digest mismatch"));
        }
        match (expected.event, expected.basis, expected.predecessor) {
            (event, MatchBasis::NewIdentity | MatchBasis::LegacyName, None) => {
                let (id, d) = declarations
                    .get(next)
                    .cloned()
                    .ok_or_else(|| r.fail("missing entity declaration"))?;
                if id != expected.subject {
                    return Err(r.fail("receipt declaration order"));
                }
                match (event, expected.basis, d.legacy) {
                    (RebindEvent::Declaration, MatchBasis::NewIdentity, false)
                    | (RebindEvent::LegacyMigration, MatchBasis::LegacyName, true) => {
                        catalog.declare(d)?;
                    }
                    (e, MatchBasis::NewIdentity, _) if e.is_revision_event() => {
                        catalog.declare_imported(d, e)?;
                    }
                    _ => return Err(r.fail("invalid declaration receipt")),
                }
                next += 1;
            }
            (RebindEvent::Rename, MatchBasis::Unchanged, None) => {
                catalog.rename(
                    expected.subject,
                    expected
                        .display_after
                        .as_deref()
                        .ok_or_else(|| r.fail("rename has no display name"))?,
                )?;
            }
            (event, MatchBasis::Unchanged, None) if event.is_revision_event() => {
                if !catalog.is_active(expected.subject) {
                    return Err(r.fail("unchanged receipt references inactive entity"));
                }
                catalog
                    .receipts
                    .try_reserve(1)
                    .map_err(|_| invalid("allocation", "catalog receipt"))?;
                catalog.append_receipt(
                    event,
                    expected.subject,
                    None,
                    MatchBasis::Unchanged,
                    None,
                    None,
                );
            }
            (event, MatchBasis::Absent, None)
                if event.is_revision_event() || event == RebindEvent::Retirement =>
            {
                catalog.retire(expected.subject, event)?;
            }
            (event, basis, Some(predecessor))
                if event.is_revision_event()
                    && matches!(
                        basis,
                        MatchBasis::GeometryFingerprint
                            | MatchBasis::DeclaredPath
                            | MatchBasis::Asserted
                    ) =>
            {
                if !catalog.is_active(expected.subject) {
                    return Err(r.fail("supersession target is not active"));
                }
                let old = catalog
                    .get(predecessor)
                    .ok_or_else(|| r.fail("unknown supersession predecessor"))?;
                let new = catalog
                    .get(expected.subject)
                    .ok_or_else(|| r.fail("unknown supersession target"))?;
                if basis == MatchBasis::GeometryFingerprint
                    && (old.kind() != new.kind()
                        || old.fingerprint().is_none()
                        || old.fingerprint() != new.fingerprint())
                {
                    return Err(
                        r.fail("content-matched receipt has no equal geometry fingerprints")
                    );
                }
                if basis == MatchBasis::DeclaredPath
                    && (old.kind() != new.kind()
                        || old.declared_name() != new.declared_name()
                        || !catalog.parent_corresponds(old.parent(), new.parent()))
                {
                    return Err(r.fail("path-matched receipt has no corresponding declared path"));
                }
                catalog.supersede(predecessor, expected.subject, event, basis)?;
            }
            _ => return Err(r.fail("invalid identity receipt operation")),
        }
        if catalog.receipts.last() != Some(&expected) {
            return Err(r.fail("receipt does not match replayed identity operation"));
        }
    }
    if next != declarations.len() || catalog.chain != r.hash()? {
        return Err(r.fail("catalog receipt root/declaration count mismatch"));
    }
    let n = r.count(
        r.budget.max_nested_items.min(budget.max_datums),
        "catalog records",
    )?;
    for row in 0..n {
        let id = DatumId(r.hash()?);
        let owner = r.entity()?;
        let declared_name = r.text()?;
        crate::assembly::text(&declared_name, "datum name", r.budget)?;
        if catalog.get(owner).is_none() || owner.kind() == EntityKind::Interface {
            return Err(r.fail("invalid datum owner"));
        }
        let feature = match r.tag()? {
            1 => DatumFeature::Plane,
            2 => DatumFeature::Axis,
            3 => DatumFeature::Point,
            _ => return Err(r.fail("datum feature")),
        };
        let count = r.count(MAX_DATUM_FRAME_LEN, "datum references")?;
        let mut references = Vec::new();
        for _ in 0..count {
            let id = DatumId(r.hash()?);
            if catalog.datum(id).is_none() {
                return Err(r.fail("datum reference is not declared"));
            }
            references.push(id);
        }
        let mut hasher = DomainHasher::new(DATUM_ID_DOMAIN);
        absorb_entity_id(&mut hasher, Some(owner));
        absorb_bytes(&mut hasher, declared_name.as_bytes());
        hasher.update(&[feature.tag()]);
        hasher.update(&(references.len() as u64).to_le_bytes());
        for id in &references {
            hasher.update(id.0.as_bytes());
        }
        if DatumId(hasher.finalize()) != id || catalog.datum(id).is_some() {
            return Err(r.fail("datum digest or uniqueness mismatch"));
        }
        catalog.datums.push(Datum {
            id,
            owner,
            declared_name,
            feature,
            references,
            row,
        });
        catalog.datum_index.push((id, row));
        catalog.datum_index.sort_by_key(|(id, _)| *id);
    }
    let n = r.count(
        r.budget.max_nested_items.min(budget.max_tolerances),
        "catalog records",
    )?;
    for _ in 0..n {
        let subject = r.entity()?;
        let kind = match r.tag()? {
            0 => ToleranceKind::Position,
            1 => ToleranceKind::Flatness,
            2 => ToleranceKind::Parallelism,
            3 => ToleranceKind::Perpendicularity,
            4 => ToleranceKind::Profile,
            5 => ToleranceKind::InterfaceGap,
            _ => return Err(r.fail("tolerance kind")),
        };
        let magnitude = r.qty()?;
        let count = r.count(MAX_DATUM_FRAME_LEN, "datum frame")?;
        let mut datum_frame = Vec::new();
        for _ in 0..count {
            datum_frame.push(DatumId(r.hash()?));
        }
        catalog.declare_tolerance(Tolerance::new(
            subject,
            kind,
            magnitude,
            datum_frame,
            r.source()?,
        ))?;
    }
    let n = r.count(
        r.budget.max_nested_items.min(budget.max_placements),
        "catalog records",
    )?;
    for row in 0..n {
        let occurrence = r.entity()?;
        let frame = FrameId(u32::try_from(r.integer()?).map_err(|_| r.fail("frame overflow"))?);
        let basis = match r.tag()? {
            0 => PlacementBasis::Nominal,
            1 => {
                let registration_ref = r.hash()?;
                if registration_ref.0 == [0; 32] {
                    return Err(r.fail("zero placement registration"));
                }
                PlacementBasis::AsBuilt { registration_ref }
            }
            _ => return Err(r.fail("placement basis")),
        };
        if catalog.get(occurrence).is_none() || occurrence.kind() != EntityKind::Part {
            return Err(r.fail("invalid placement occurrence"));
        }
        catalog.placements.push(Placement {
            occurrence,
            frame,
            basis,
            row,
        });
    }
    Ok(catalog)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn catalog_wire_refuses_a_forged_content_matched_receipt() {
        let mut c = EntityCatalog::new();
        let root = c
            .declare(EntityDeclaration::assembly("test-fixture:root"))
            .unwrap();
        let a = c
            .declare(
                EntityDeclaration::part(root, "test-fixture:a")
                    .with_fingerprint(GeometryFingerprint::of_bytes(b"a")),
            )
            .unwrap();
        let b = c
            .declare(
                EntityDeclaration::part(root, "test-fixture:b")
                    .with_fingerprint(GeometryFingerprint::of_bytes(b"b")),
            )
            .unwrap();
        // The private primitive normally trusts auto_match. Forge a locally
        // hash-consistent higher-tier receipt to falsify the wire boundary.
        c.supersede(a, b, RebindEvent::Import, MatchBasis::GeometryFingerprint)
            .unwrap();
        assert!(c.verify_receipts());
        let budget = crate::assembly::DEFAULT_ASSEMBLY_BUDGET;
        let mut w = Writer::new(budget);
        write_assembly_catalog(&c, &mut w, budget);
        let bytes = w.finish().unwrap();
        let mut reader = Reader::new(&bytes, budget).unwrap();
        assert!(
            matches!(read_assembly_catalog(&mut reader),Err(AssemblyError::Wire {detail,..}) if detail.contains("equal geometry fingerprints"))
        );
    }
}
