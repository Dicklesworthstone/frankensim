//! Exact, bounded snapshots for the projected executable. Every snapshot is an
//! immutable file; a torn new file cannot overwrite an earlier accepted state.
//! The checksum detects damaged bytes, not physical/model validity. Recovery
//! still re-admits the constraints and independently solves the retained field.
use fs_topols::projected::{ProjectedOptimizer, ProjectedSettings, ProjectedSetupStage};
use fs_topols::volume::VolumeProjectionSettings;
use fs_topols::{Cantilever, GridSdf, OptimizeCheckpoint, OptimizeSettings};
use std::error::Error;
use std::fs::{File, OpenOptions};
use std::io::{Read, Write};
use std::ops::ControlFlow;
use std::path::Path;

const MAGIC: &[u8] = b"fs-marquee-projected-checkpoint-v1\0";
const MAX_BYTES: u64 = 8 * 1024 * 1024;

pub(super) fn filename(iteration: usize) -> String { format!("checkpoint-{iteration:06}.bin") }

fn integer(out: &mut Vec<u8>, value: u64) { out.extend_from_slice(&value.to_le_bytes()); }
fn number(out: &mut Vec<u8>, value: f64) { integer(out, value.to_bits()); }

fn encode(
    checkpoint: &OptimizeCheckpoint,
    fixed: &[(usize, f64)],
    projection: VolumeProjectionSettings,
    controls: ProjectedSettings,
) -> Vec<u8> {
    let mut body = Vec::new();
    let settings = checkpoint.settings();
    integer(&mut body, u64::from(settings.level));
    integer(&mut body, settings.iterations as u64);
    number(&mut body, checkpoint.fixture().load);
    number(&mut body, checkpoint.fixture().band);
    for value in [settings.volfrac, settings.band_cells, settings.move_cells,
        settings.ell0, settings.mu_al, settings.sobolev_alpha] { number(&mut body, value); }
    integer(&mut body, settings.nucleation_period as u64);
    for value in [settings.hole_radius_cells, settings.youngs, settings.poisson] {
        number(&mut body, value);
    }
    integer(&mut body, checkpoint.next_iteration() as u64);
    number(&mut body, checkpoint.ell());
    for value in [projection.target, projection.tolerance, projection.max_shift] {
        number(&mut body, value);
    }
    integer(&mut body, projection.max_evaluations as u64);
    integer(&mut body, controls.max_candidates as u64);
    number(&mut body, controls.contraction);
    number(&mut body, controls.min_relative_improvement);
    integer(&mut body, controls.poll_iters as u64);
    integer(&mut body, fixed.len() as u64);
    for &(index, value) in fixed {
        integer(&mut body, index as u64);
        number(&mut body, value);
    }
    for &value in checkpoint.geometry().nodes() { number(&mut body, value); }
    seal(&body)
}

fn seal(body: &[u8]) -> Vec<u8> {
    let mut bytes = MAGIC.to_vec();
    bytes.extend_from_slice(fs_ledger::hash_bytes(body).to_hex().as_bytes());
    bytes.push(b'\n');
    bytes.extend_from_slice(body);
    bytes
}

struct Reader<'a>(&'a [u8]);
impl Reader<'_> {
    fn integer(&mut self) -> Result<u64, Box<dyn Error>> {
        let raw = self.0.get(..8).ok_or("truncated projected checkpoint")?;
        let value = u64::from_le_bytes(raw.try_into()?);
        self.0 = &self.0[8..];
        Ok(value)
    }
    fn count(&mut self) -> Result<usize, Box<dyn Error>> { Ok(self.integer()?.try_into()?) }
    fn number(&mut self) -> Result<f64, Box<dyn Error>> { Ok(f64::from_bits(self.integer()?)) }
}

type Decoded = (OptimizeCheckpoint, Vec<(usize, f64)>, VolumeProjectionSettings, ProjectedSettings);

fn decode(bytes: &[u8]) -> Result<Decoded, Box<dyn Error>> {
    let header = MAGIC.len() + 64 + 1;
    if bytes.len() as u64 > MAX_BYTES || bytes.len() < header || !bytes.starts_with(MAGIC)
        || bytes[header - 1] != b'\n'
    {
        return Err("unsupported, oversized or truncated projected checkpoint".into());
    }
    let body = &bytes[header..];
    if &bytes[MAGIC.len()..header - 1] != fs_ledger::hash_bytes(body).to_hex().as_bytes() {
        return Err("projected checkpoint checksum mismatch".into());
    }
    let mut reader = Reader(body);
    let level: u32 = reader.integer()?.try_into()?;
    let iterations = reader.count()?;
    // Bound the allocation/work envelope BEFORE allocating a level-set lattice.
    if !(2..=7).contains(&level) || !(1..=200).contains(&iterations) {
        return Err("checkpoint exceeds projected executable level/update limits".into());
    }
    let fixture = Cantilever { load: reader.number()?, band: reader.number()? };
    let settings = OptimizeSettings {
        level, iterations, volfrac: reader.number()?, band_cells: reader.number()?,
        move_cells: reader.number()?, ell0: reader.number()?, mu_al: reader.number()?,
        sobolev_alpha: reader.number()?, nucleation_period: reader.count()?,
        hole_radius_cells: reader.number()?, youngs: reader.number()?, poisson: reader.number()?,
    };
    let next_iteration = reader.count()?;
    let ell = reader.number()?;
    let projection = VolumeProjectionSettings {
        target: reader.number()?, tolerance: reader.number()?, max_shift: reader.number()?,
        max_evaluations: reader.count()?,
    };
    let controls = ProjectedSettings {
        max_candidates: reader.count()?, contraction: reader.number()?,
        min_relative_improvement: reader.number()?, poll_iters: reader.count()?,
    };
    if !(1..=16).contains(&controls.max_candidates) || !(1..=1024).contains(&controls.poll_iters)
        || !(1..=256).contains(&projection.max_evaluations)
    {
        return Err("checkpoint exceeds projected candidate/projection/polling limits".into());
    }
    let n = 1usize << level;
    let node_count = (n + 1) * (n + 1);
    let fixed_count = reader.count()?;
    if fixed_count != 2 * (n + 1) {
        return Err("projected checkpoint must retain the complete left/right boundary policy".into());
    }
    let mut fixed = Vec::with_capacity(fixed_count);
    for ordinal in 0..fixed_count {
        let index = reader.count()?;
        let value = reader.number()?;
        let expected = (ordinal / 2) * (n + 1) + if ordinal % 2 == 0 { 0 } else { n };
        if index != expected || !value.is_finite() {
            return Err("checkpoint contains an invalid fixed-node assignment".into());
        }
        fixed.push((index, value));
    }
    if reader.0.len() != node_count * 8 {
        return Err("checkpoint must contain exactly one complete nodal lattice".into());
    }
    let mut geometry = GridSdf::from_fn(n, &|_, _| 0.0);
    for j in 0..=n {
        for i in 0..=n { *geometry.node_mut(i, j) = reader.number()?; }
    }
    Ok((OptimizeCheckpoint::restore(geometry, fixture, settings, next_iteration, ell)?,
        fixed, projection, controls))
}

pub(super) fn save(directory: &Path, optimizer: &ProjectedOptimizer) -> Result<(), Box<dyn Error>> {
    let bytes = encode(optimizer.checkpoint(), optimizer.fixed_nodes(),
        optimizer.projection_settings(), optimizer.controls());
    let path = directory.join(filename(optimizer.checkpoint().next_iteration()));
    let mut output = OpenOptions::new().write(true).create_new(true).open(path)?;
    output.write_all(&bytes)?;
    output.sync_all()?;
    Ok(())
}

pub(super) fn load_controlled<B>(
    path: &Path,
    mut control: impl FnMut(ProjectedSetupStage) -> ControlFlow<B>,
) -> Result<ControlFlow<B, ProjectedOptimizer>, Box<dyn Error>> {
    if let ControlFlow::Break(reason) = control(ProjectedSetupStage::Prepare) {
        return Ok(ControlFlow::Break(reason));
    }
    let mut bytes = Vec::new();
    File::open(path)?.take(MAX_BYTES + 1).read_to_end(&mut bytes)?;
    let (checkpoint, fixed, projection, controls) = decode(&bytes)?;
    // Never rebuild fixed nodes from the imported field or re-project on resume.
    // The canonical owner verifies exact fixed bits and area and re-solves it.
    Ok(ProjectedOptimizer::from_checkpoint_controlled(
        &checkpoint, fixed, projection, controls, control,
    )?)
}

#[cfg(test)]
fn load(path: &Path) -> Result<ProjectedOptimizer, Box<dyn Error>> {
    match load_controlled(path, |_| ControlFlow::<std::convert::Infallible>::Continue(()))? {
        ControlFlow::Continue(optimizer) => Ok(optimizer),
        ControlFlow::Break(never) => match never {},
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn encoded() -> Vec<u8> {
        let settings = OptimizeSettings { level: 3, iterations: 5, youngs: 2.0,
            poisson: -0.2, nucleation_period: 2, ..OptimizeSettings::default() };
        let geometry = GridSdf::from_fn(8, &|x, y| ((y - 0.5).abs() - 0.35) * (1.0 + x));
        let fixed: Vec<_> = geometry.nodes().iter().copied().enumerate()
            .filter(|(index, _)| index % 9 == 0 || index % 9 == 8).collect();
        let checkpoint = OptimizeCheckpoint::restore(geometry,
            Cantilever { load: 1.5, band: 0.125 }, settings, 2, 0.375).unwrap();
        encode(&checkpoint, &fixed, VolumeProjectionSettings {
            target: 0.5, tolerance: 1e-4, max_shift: 2.0, max_evaluations: 64,
        }, ProjectedSettings::default())
    }

    #[test]
    fn exact_state_and_policy_round_trip_without_replaying_updates() {
        let bytes = encoded();
        let (checkpoint, fixed, projection, controls) = decode(&bytes).unwrap();
        assert_eq!(checkpoint.next_iteration(), 2);
        assert_eq!(checkpoint.ell().to_bits(), 0.375_f64.to_bits());
        assert_eq!(checkpoint.settings().youngs, 2.0);
        assert_eq!(checkpoint.settings().poisson, -0.2);
        assert_eq!(checkpoint.settings().nucleation_period, 2);
        assert_eq!(checkpoint.fixture().load, 1.5);
        assert_eq!(encode(&checkpoint, &fixed, projection, controls), bytes);
    }

    #[test]
    fn torn_corrupt_and_unknown_snapshots_are_refused() {
        let bytes = encoded();
        for length in [0, MAGIC.len(), MAGIC.len() + 64, bytes.len() - 1] {
            assert!(decode(&bytes[..length]).is_err());
        }
        let mut changed = bytes.clone();
        *changed.last_mut().unwrap() ^= 1;
        assert!(decode(&changed).is_err());
        let mut extra = bytes.clone();
        extra.push(0);
        assert!(decode(&extra).is_err());
        let mut wrong_version = bytes;
        wrong_version[0] ^= 1;
        assert!(decode(&wrong_version).is_err());
    }

    #[test]
    fn correctly_checksummed_oversized_lattice_is_refused_before_allocation() {
        let mut body = Vec::new();
        integer(&mut body, 63);
        integer(&mut body, 5);
        assert!(decode(&seal(&body)).is_err());
        let bytes = encoded();
        let mut body = bytes[MAGIC.len() + 65..].to_vec();
        body.push(0);
        assert!(decode(&seal(&body)).is_err(), "valid checksum is not shape admission");
    }

    #[test]
    fn a_valid_checksum_cannot_remove_prescribed_boundaries() {
        let bytes = encoded();
        let (checkpoint, mut fixed, projection, controls) = decode(&bytes).unwrap();
        fixed.pop();
        assert!(decode(&encode(&checkpoint, &fixed, projection, controls)).is_err());
        let (_, mut fixed, _, _) = decode(&bytes).unwrap();
        fixed.swap(0, 1);
        assert!(decode(&encode(&checkpoint, &fixed, projection, controls)).is_err());
    }

    #[test]
    fn accepted_state_round_trip_continues_the_same_numerical_trajectory() {
        let geometry = GridSdf::from_fn(8, &|_, y| (y - 0.5).abs() - 0.35);
        let fixed = geometry.nodes().iter().copied().enumerate()
            .filter(|(index, _)| index % 9 == 0 || index % 9 == 8).collect();
        let settings = OptimizeSettings { level: 3, iterations: 2, volfrac: 0.6,
            move_cells: 0.1, nucleation_period: 0, ..OptimizeSettings::default() };
        let mut original = ProjectedOptimizer::new(geometry, Cantilever { load: 1.0, band: 0.125 },
            settings, fixed, VolumeProjectionSettings {
                target: 0.6, tolerance: 1e-4, max_shift: 2.0, max_evaluations: 64,
            }, ProjectedSettings { max_candidates: 8, ..ProjectedSettings::default() }).unwrap();
        assert!(matches!(original.advance_one().unwrap(),
            fs_topols::projected::ProjectedProgress::Accepted(_)),
            "use the backend's non-vacuous accepted-step fixture");
        let root = std::env::temp_dir().join(format!("frankensim-projected-accepted-checkpoint-{}", std::process::id()));
        std::fs::create_dir(&root).unwrap();
        save(&root, &original).unwrap();
        let mut resumed = load(&root.join(filename(1))).unwrap();
        assert_eq!(resumed.checkpoint().next_iteration(), 1);
        assert_eq!(resumed.current().compliance.to_bits(), original.current().compliance.to_bits());
        let a = original.advance_one().unwrap();
        let b = resumed.advance_one().unwrap();
        assert_eq!(std::mem::discriminant(&a), std::mem::discriminant(&b));
        assert_eq!(resumed.current().compliance.to_bits(), original.current().compliance.to_bits());
        assert_eq!(encode(resumed.checkpoint(), resumed.fixed_nodes(),
            resumed.projection_settings(), resumed.controls()),
            encode(original.checkpoint(), original.fixed_nodes(),
                original.projection_settings(), original.controls()));
    }

    #[test]
    fn cancelled_recovery_stops_before_opening_or_decoding_a_checkpoint() {
        let mut polls = 0;
        let result = load_controlled(Path::new("unused-cancelled-checkpoint.bin"), |stage| {
            polls += 1;
            assert_eq!(stage, ProjectedSetupStage::Prepare);
            ControlFlow::Break("budget")
        }).unwrap();
        assert!(matches!(result, ControlFlow::Break("budget")));
        assert_eq!(polls, 1);
    }
}
