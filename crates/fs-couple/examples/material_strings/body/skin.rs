//! A closed acoustic skin of a shell body, for the one-way observer. Each skin
//! vertex moves as a declared weighted sum of points on the profile's own parts
//! (component, triangle, barycentric), and the observer checks every anchored
//! point against the vertex before any boundary work. Aperture caps carry
//! declared shares of the lumped cavity's neck flow. Nothing here is inferred:
//! the generator that writes the profile states how the skin was made.
use super::*;
use fs_blake3::ContentHash;
use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::performance::observation::one_way::{
    MAX_ONE_WAY_SETUP_BYTES, OneWaySetupStore, ShellAcousticSkin, SkinAnchor, SkinAperture,
};
use fs_bem::modal_loading::{MAX_SWEEP_FREQUENCY_BYTES, SweepFrequencyStore};
use fs_couple::bowed_string::runtime::reciprocal::bank::radiation::bem::{
    MAX_SKIN_LOAD_BYTES, SkinLoadStore,
};
use std::fs::File;
use std::io::{ErrorKind, Read as _, Write as _};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

const MAX_VERTICES: usize = fs_bem::helmholtz::MAX_DENSE_PANELS;
const MAX_ANCHORS: usize = 8;
const MAX_APERTURES: usize = 8;

/// What a [`SetupDir`] keeps: its option, its entries' suffix and their byte ceiling.
#[derive(Clone, Copy, Debug)]
pub(super) struct Kept {
    option: &'static str,
    suffix: &'static str,
    maximum_bytes: usize,
}
/// `--observer-store`: one-way observer setups, under their setup keys.
pub(super) const ONE_WAY_SETUPS: Kept = Kept {
    option: "--observer-store",
    suffix: "one-way-setup",
    maximum_bytes: MAX_ONE_WAY_SETUP_BYTES,
};
/// `--skin-load-store`: two-way skin loads, under their loading identities.
pub(super) const SKIN_LOADS: Kept = Kept {
    option: "--skin-load-store",
    suffix: "skin-load",
    maximum_bytes: MAX_SKIN_LOAD_BYTES,
};
/// `--skin-load-store`: each frequency's solve of a two-way skin load's sweep, under its frequency key, so that a bake
/// stopped part way resumes from the frequencies it solved.
pub(super) const SKIN_FREQUENCIES: Kept = Kept {
    option: "--skin-load-store",
    suffix: "skin-frequency",
    maximum_bytes: MAX_SWEEP_FREQUENCY_BYTES,
};
/// `--skin-load-store`: an assembled body's modal solve, under its solve key.
pub(super) const SHELL_MODES: Kept = Kept {
    option: "--skin-load-store",
    suffix: "shell-modes",
    maximum_bytes: fs_couple::shear_shell::assembly::MAX_SHELL_MODES_ENTRY_BYTES,
};

/// `--observer-store DIR` or `--skin-load-store DIR`: this executable's
/// one-way observer setups, or two-way skin loads, their sweeps' frequencies
/// and assembled bodies' modal solves, kept in an existing directory. Each entry's name holds its key and
/// the start of the executable's own hash, because the key cannot see a solver
/// change that keeps its hash domains, so a rebuilt binary bakes afresh. An
/// entry is written once and never replaced.
#[derive(Clone, Debug)]
pub(super) struct SetupDir {
    dir: PathBuf,
    build: String,
    kept: Kept,
}

impl SetupDir {
    pub(super) fn open(dir: &Path, kept: Kept) -> Result<Self> {
        if !dir.is_dir() {
            return Err(format!("{} needs an existing directory", kept.option).into());
        }
        let executable = std::fs::read(std::env::current_exe()?)?;
        let build = fs_blake3::hash_domain(
            "frankensim:material-strings:observer-store-build:v1",
            &executable,
        );
        Ok(Self::with_build(dir, build, kept))
    }

    fn with_build(dir: &Path, build: ContentHash, kept: Kept) -> Self {
        Self {
            dir: dir.to_path_buf(),
            build: build.to_hex()[..32].into(),
            kept,
        }
    }

    /// The same directory and build, keeping `kept`.
    pub(super) fn keeping(self, kept: Kept) -> Self {
        Self { kept, ..self }
    }

    pub(super) fn entry(&self, key: &ContentHash) -> PathBuf {
        self.dir
            .join(format!("{}-{}.{}", key.to_hex(), self.build, self.kept.suffix))
    }

    fn read(&self, key: &ContentHash) -> std::result::Result<Option<Vec<u8>>, String> {
        let path = self.entry(key);
        let failed = |e: std::io::Error| format!("{}: {e}", path.display());
        match File::open(&path) {
            Ok(file) => {
                let mut bytes = Vec::new();
                file.take(self.kept.maximum_bytes as u64 + 1)
                    .read_to_end(&mut bytes)
                    .map_err(failed)?;
                Ok(Some(bytes))
            }
            Err(e) if e.kind() == ErrorKind::NotFound => Ok(None),
            Err(e) => Err(failed(e)),
        }
    }

    /// Writes and syncs `entry` under a temporary name of its own beside it (the entry's name and `.PID-N.tmp`),
    /// then hard-links it to the entry's name, which an existing entry keeps: a write stopped part way leaves only
    /// the temporary, which no read opens. A file system without hard links refuses every entry.
    fn write(&self, key: &ContentHash, entry: &[u8]) -> std::result::Result<(), String> {
        static WRITES: AtomicU64 = AtomicU64::new(0);
        let path = self.entry(key);
        let mut temporary = path.clone().into_os_string();
        let n = WRITES.fetch_add(1, Ordering::Relaxed);
        temporary.push(format!(".{}-{n}.tmp", std::process::id()));
        let temporary = PathBuf::from(temporary);
        let written = File::create(&temporary).and_then(|mut file| {
            file.write_all(entry)?;
            file.sync_all()
        });
        let linked = written.and_then(|()| match std::fs::hard_link(&temporary, &path) {
            Err(e) if e.kind() == ErrorKind::AlreadyExists => Ok(()),
            linked => linked,
        });
        let removed = std::fs::remove_file(&temporary);
        linked
            .and(removed)
            .map_err(|e| format!("{}: {e}", path.display()))
    }
}

impl OneWaySetupStore for SetupDir {
    fn load(&self, key: &ContentHash) -> std::result::Result<Option<Vec<u8>>, String> {
        self.read(key)
    }

    fn keep(&self, key: &ContentHash, entry: &[u8]) -> std::result::Result<(), String> {
        self.write(key, entry)
    }
}

impl SkinLoadStore for SetupDir {
    fn load(&self, key: &ContentHash) -> std::result::Result<Option<Vec<u8>>, String> {
        self.read(key)
    }

    fn keep(&self, key: &ContentHash, entry: &[u8]) -> std::result::Result<(), String> {
        self.write(key, entry)
    }

    fn frequencies(&self) -> Option<&dyn SweepFrequencyStore> {
        Some(self)
    }
}

/// A skin load's frequencies, in the same directory and build, as [`SKIN_FREQUENCIES`].
impl SweepFrequencyStore for SetupDir {
    fn load(&self, key: &ContentHash) -> std::result::Result<Option<Vec<u8>>, String> {
        self.clone().keeping(SKIN_FREQUENCIES).read(key)
    }

    fn keep(&self, key: &ContentHash, entry: &[u8]) -> std::result::Result<(), String> {
        self.clone().keeping(SKIN_FREQUENCIES).write(key, entry)
    }
}

#[derive(Clone, Debug)]
pub(super) struct Spec {
    // x/y scale with the original performance dimensions; z is in metres, as the shell nodes.
    vertices: Vec<[f64; 3]>,
    anchors: Vec<Vec<SkinAnchor>>,
    triangles: Vec<[usize; 3]>,
    apertures: Vec<SkinAperture>,
}

/// The optional skin follows a complete v3/v4 body and precedes any interior
/// volume (which [`super::cavity::split`] has already taken off the tail).
pub(super) fn split(text: &str) -> Result<(&str, Option<Spec>)> {
    let mut offset = 0;
    for line in text.split_inclusive('\n') {
        if line.split_ascii_whitespace().next() == Some("acoustic-skin-v1") {
            return Ok((&text[..offset], Some(Spec::parse(&text[offset..])?)));
        }
        offset += line.len();
    }
    Ok((text, None))
}

fn count(value: &str, minimum: usize, maximum: usize) -> Result<usize> {
    let value: usize = value.parse()?;
    if !(minimum..=maximum).contains(&value) {
        return Err("acoustic skin count exceeds its explicit range".into());
    }
    Ok(value)
}

impl Spec {
    fn parse(text: &str) -> Result<Self> {
        let mut r = Records {
            lines: text.lines(),
            line: 0,
        };
        r.row("acoustic-skin-v1", 0)?;
        let n = count(r.row("skin-vertices", 1)?[0], 3, MAX_VERTICES)?;
        let (mut vertices, mut anchors) = (Vec::with_capacity(n), Vec::with_capacity(n));
        for _ in 0..n {
            let v = r.row("skin-vertex", 4)?;
            vertices.push([scalar(v[0])?, scalar(v[1])?, scalar(v[2])?]);
            let k = count(v[3], 0, MAX_ANCHORS)?;
            let mut points = Vec::with_capacity(k);
            for _ in 0..k {
                let a = r.row("skin-anchor", 6)?;
                points.push(SkinAnchor {
                    weight: scalar(a[0])?,
                    component: a[1].parse()?,
                    triangle: a[2].parse()?,
                    barycentric: [scalar(a[3])?, scalar(a[4])?, scalar(a[5])?],
                });
            }
            anchors.push(points);
        }
        let m = count(
            r.row("skin-triangles", 1)?[0],
            4,
            fs_bem::helmholtz::MAX_DENSE_PANELS,
        )?;
        let mut triangles = Vec::with_capacity(m);
        for _ in 0..m {
            let v = r.row("skin-triangle", 3)?;
            triangles.push([
                count(v[0], 0, n - 1)?,
                count(v[1], 0, n - 1)?,
                count(v[2], 0, n - 1)?,
            ]);
        }
        let p = count(r.row("skin-apertures", 1)?[0], 0, MAX_APERTURES)?;
        let mut apertures = Vec::<SkinAperture>::with_capacity(p);
        for _ in 0..p {
            let v = r.row("skin-aperture", 3)?;
            if v[0].len() > 128 || apertures.iter().any(|a| a.id == v[0]) {
                return Err("skin aperture ids must be unique bounded tokens".into());
            }
            let caps = count(v[2], 1, m)?;
            let mut cap = Vec::with_capacity(caps);
            for _ in 0..caps {
                cap.push(count(r.row("skin-cap", 1)?[0], 0, m - 1)?);
            }
            apertures.push(SkinAperture {
                id: v[0].into(),
                triangles: cap,
                flow_share: scalar(v[1])?,
            });
        }
        if r.lines.next().is_some() {
            return Err("the acoustic skin must end the body section".into());
        }
        // Closure, orientation, anchors and shares are the observer's own admission.
        Ok(Self {
            vertices,
            anchors,
            triangles,
            apertures,
        })
    }

    /// Source admission before material I/O; the observer repeats every check on the built parts.
    pub(super) fn preflight(&self, spec: &input::Spec) -> Result<()> {
        let load = match spec.body_radiation {
            Some(input::BodyRadiation::Skin(load)) => Some(load),
            _ => None,
        };
        if load.is_some_and(|load| self.triangles.len() > load.options.maximum_surface_panels) {
            return Err(
                "acoustic skin exceeds the body-radiation skin load's declared panels".into(),
            );
        }
        if let Some(observer) = &spec.one_way_observer {
            if self.triangles.len() > observer.maximum_skin_panels {
                return Err("acoustic skin exceeds the one_way_observer's declared panels".into());
            }
        } else if load.is_none() {
            return Ok(());
        }
        let compact =
            |c: &fs_couple::bowed_string::runtime::reciprocal::plate::MaterialBodyCavity| {
                c.cavity.loss() != fs_couple::cavity_modal::CavityNeckLoss::ExteriorMouth
            };
        if load.is_some() && spec.cavity.as_ref().is_some_and(compact) {
            return Err(
                "a compact outer neck end already carries the mouth's exterior load; declare the cavity exterior-neck"
                    .into(),
            );
        }
        if spec.cavity.is_some() == self.apertures.is_empty() {
            return Err(
                "a lumped cavity radiates through declared skin apertures, and only then".into(),
            );
        }
        Ok(())
    }

    /// The skin in the performance body's frame [m].
    pub(super) fn resolve(&self, spec: &input::Spec) -> ShellAcousticSkin {
        ShellAcousticSkin {
            vertices: self
                .vertices
                .iter()
                .map(|p| [p[0] * spec.body.length, p[1] * spec.body.width, p[2]])
                .collect(),
            anchors: self.anchors.clone(),
            triangles: self.triangles.clone(),
            apertures: self.apertures.clone(),
        }
    }
}
