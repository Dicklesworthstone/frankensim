//! sample_field over a snapshot LEASE (bead wf-root-guzez.8.1.2,
//! E7.1-ii, plan §5.5 Round-4 S-02). Callers supply a LEASED
//! immutable E5.0 ring snapshot — never a raw tick as a state
//! substitute. A stale snapshot yields a valid HISTORICAL result
//! under its own `FieldSourceSnapshotId`, never mislabeled current.
//!
//! The field state derives from the leased payload alone: an
//! airborne payload publishes the lift-carrying bound system
//! (Kutta–Joukowski at the published gross mass); the atmosphere is
//! the canonical v1 ambient built from the run scenario. Everything
//! else — masks, duals, meta, forbidden claims — is fs-flyer's
//! fieldsvc, unchanged: one physics, two entry points.

use crate::Refusal;
use crate::ring::{ABI_VERSION, SnapshotRing};
use fs_blake3::hash_domain;
use fs_flyer::fieldsvc::{
    FieldSampleSet, FieldSourceStateV1, GridSpec, ambient_atmosphere_v1, sample_field,
};
use fs_flyer::simloop::SNAPSHOT_LEN;

fn refuse(code: &'static str, message: String, repair: &str) -> Refusal {
    Refusal {
        code,
        message,
        ranked_repairs: vec![repair.into()],
    }
}

/// Staleness threshold [ticks] beyond which a result is labeled
/// HISTORICAL (2 s at 120 Hz).
pub const STALENESS_THRESHOLD_TICKS: u64 = 240;

/// Point cap for the JSON-exporting wasm path (the full native cap
/// stays available to array consumers).
pub const MAX_JSON_POINTS: usize = 4_096;

/// A lease-bound field sample: the §5.5 set plus the lease receipt.
#[derive(Clone, Debug)]
pub struct LeasedFieldSample {
    /// The §5.5 sample (meta carries the snapshot's OWN id).
    pub sample: FieldSampleSet,
    /// Labeled historical? (staleness beyond the threshold).
    pub historical: bool,
    /// How stale the snapshot was when sampled [ticks].
    pub staleness_ticks: u64,
    /// Ring slot the lease held.
    pub slot: u32,
    /// Run epoch the lease validated against.
    pub run_epoch: u64,
}

/// Sample the field under a fresh ring lease.
///
/// # Errors
/// Ring refusals (torn/stale-epoch/ABI) propagate; payload, grid,
/// mask, and atmosphere refusals propagate from fieldsvc.
pub fn sample_field_leased(
    ring: &mut SnapshotRing,
    current_tick: u64,
    scenario_seed: u64,
    headwind_mps: f64,
    rho_kg_m3: f64,
    grid: &GridSpec,
    component_mask: u32,
) -> Result<LeasedFieldSample, Refusal> {
    let lease = ring.acquire(ABI_VERSION)?;
    // The v1 payload layout is versioned: exactly SNAPSHOT_LEN floats
    // (the ring refuses any other reader buffer size).
    let mut payload = [0.0f64; SNAPSHOT_LEN];
    let read = ring.read(&lease, &mut payload);
    // The lease is released on EVERY path — torn reads included.
    let header = match read {
        Ok(h) => h,
        Err(e) => {
            ring.release(&lease);
            return Err(e);
        }
    };
    ring.release(&lease);
    let n = header.payload_len as usize;
    // The snapshot's OWN identity: header identity fields + payload
    // bits. A stale snapshot keeps this id — it is never re-labeled.
    let mut b = Vec::with_capacity(8 * (n + 4));
    for v in [
        header.run_epoch,
        header.run_anchor_digest_prefix,
        header.model_id_prefix,
        header.tick,
    ] {
        b.extend_from_slice(&v.to_le_bytes());
    }
    for v in &payload[..n] {
        b.extend_from_slice(&v.to_bits().to_le_bytes());
    }
    let source_digest = hash_domain("org.frankensim.wf.leased-snapshot.v1", &b).to_hex();
    let atmosphere = ambient_atmosphere_v1(scenario_seed, headwind_mps).map_err(Refusal::from)?;
    let state = FieldSourceStateV1::from_ring_payload(
        header.tick,
        source_digest,
        &payload[..n],
        atmosphere,
        rho_kg_m3,
    )
    .map_err(Refusal::from)?;
    let sample = sample_field(&state, grid, component_mask).map_err(Refusal::from)?;
    let staleness = current_tick.saturating_sub(header.tick);
    Ok(LeasedFieldSample {
        sample,
        historical: staleness > STALENESS_THRESHOLD_TICKS,
        staleness_ticks: staleness,
        slot: lease.slot,
        run_epoch: lease.run_epoch,
    })
}

/// Round-4 S-02: a raw tick is FORBIDDEN as a state substitute. The
/// internal replay convenience (reconstruct → publish → lease) is not
/// wired in v1, so this entry point refuses TYPED — it never invents
/// a field from a tick number.
///
/// # Errors
/// Always `tick-as-state-forbidden`.
pub fn sample_field_by_tick(_tick: u64) -> Result<LeasedFieldSample, Refusal> {
    Err(refuse(
        "tick-as-state-forbidden",
        "a raw tick is not a state (Round-4 S-02)".into(),
        "acquire a snapshot lease; replay reconstruction publishes a snapshot first",
    ))
}

/// Bounded JSON summary of a leased sample (the wasm export payload):
/// meta + masks summary + payload digest — never the raw arrays.
///
/// # Errors
/// `json-points-exceeded` (AT the JSON cap admits, one more refuses).
pub fn leased_sample_json(s: &LeasedFieldSample) -> Result<String, Refusal> {
    let n = s.sample.u.len();
    if n > MAX_JSON_POINTS {
        return Err(refuse(
            "json-points-exceeded",
            format!("{n} points > {MAX_JSON_POINTS}"),
            "sample a coarser grid for the JSON path; array consumers have the full cap",
        ));
    }
    let valid = s.sample.validity_mask.iter().filter(|v| **v).count();
    let cores = s
        .sample
        .singularity_core_mask
        .iter()
        .filter(|v| **v)
        .count();
    let omitted: Vec<String> = s
        .sample
        .meta
        .omitted_components
        .iter()
        .map(|c| format!("\"{c}\""))
        .collect();
    Ok(format!(
        "{{\"schema\":\"wf-field-sample-v1\",\"field_source_snapshot_id\":\"{}\",\
         \"source_tick\":{},\"historical\":{},\"staleness_ticks\":{},\
         \"component_mask\":{},\"omitted_components\":[{}],\"n_points\":{n},\
         \"valid_points\":{valid},\"core_points\":{cores},\
         \"payload_digest\":\"{}\",\"slot\":{},\"run_epoch\":{}}}",
        s.sample.meta.field_source_snapshot_id,
        s.sample.meta.source_tick,
        s.historical,
        s.staleness_ticks,
        s.sample.component_mask,
        omitted.join(","),
        s.sample.digest(),
        s.slot,
        s.run_epoch,
    ))
}

/// Whole-path self-test (also the wasm-lane harness surface): build
/// a 3-slot ring, publish one canonical airborne payload, sample the
/// full supported sum under a lease, and return the bounded JSON.
/// Any refusal on the way returns its envelope instead.
#[must_use]
pub fn field_selftest_json() -> String {
    let inner = || -> Result<String, Refusal> {
        let mut ring = SnapshotRing::new(3, 12, 0x1903, 0xd17)?;
        let payload = [
            50.0, 3.0, 12.0, 0.8, 0.01, 0.05, 0.02, 0.0, 35.0, 0.0, 0.0, 1.0,
        ];
        ring.publish(100, &payload)?;
        let grid = GridSpec {
            origin_m: [45.0, -4.0, 0.5],
            dx_m: 2.0,
            nx: 3,
            ny: 3,
            nz: 2,
        };
        use fs_flyer::fieldsvc::{C_BOUND_CIRCULATION, C_GROUND_IMAGES, C_MEAN_ATMO, C_TURB_ATMO};
        let leased = sample_field_leased(
            &mut ring,
            110,
            1903,
            8.0,
            1.294,
            &grid,
            C_MEAN_ATMO | C_TURB_ATMO | C_BOUND_CIRCULATION | C_GROUND_IMAGES,
        )?;
        leased_sample_json(&leased)
    };
    match inner() {
        Ok(j) => j,
        Err(e) => format!(
            "{{\"schema\":\"wf-field-sample-v1\",\"refusal\":\"{}\",\"message\":\"{}\"}}",
            e.code, e.message
        ),
    }
}

impl From<fs_flyer::Refusal> for Refusal {
    fn from(e: fs_flyer::Refusal) -> Self {
        Refusal {
            code: e.code,
            message: e.message,
            ranked_repairs: e.ranked_repairs,
        }
    }
}
