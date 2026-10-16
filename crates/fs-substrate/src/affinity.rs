//! Shard-by-CCD tile placement (plan §5.1 consequence 3): tiles are
//! partitioned into contiguous Z-ORDER ranges, one shard per CCD/cluster,
//! so a shard's working set stays inside its L3 island and its pages are
//! first-touched by the owning worker. This crate produces the MAP; thread
//! pinning and actual first-touch scheduling are fs-exec's contract.

use crate::CapabilityProbe;
use crate::tile::TileGrid;
use core::ops::Range;
use std::fmt::Write as _;

/// Core-complex topology: how many shards (CCDs / clusters) and how many
/// cores each holds. Fixture constants encode the reference machines;
/// [`CcdTopology::from_probe`] is an honest heuristic, recorded not trusted.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CcdTopology {
    /// Shard count (CCDs on x86, core clusters on Apple).
    pub ccds: u32,
    /// Physical cores per shard.
    pub cores_per_ccd: u32,
}

impl CcdTopology {
    /// Threadripper PRO 7995WX fixture: 12 CCDs × 8 cores (plan §5.1).
    pub const TR_7995WX: CcdTopology = CcdTopology {
        ccds: 12,
        cores_per_ccd: 8,
    };

    /// 128-core EPYC-class fixture: 16 CCDs × 8 cores.
    pub const EPYC_128C: CcdTopology = CcdTopology {
        ccds: 16,
        cores_per_ccd: 8,
    };

    /// Apple M4 Pro/Max-class fixture: one P cluster + one E cluster
    /// modeled as two shards (UMA makes sharding a locality nicety, not a
    /// bandwidth necessity).
    pub const APPLE_M_CLASS: CcdTopology = CcdTopology {
        ccds: 2,
        cores_per_ccd: 8,
    };

    /// Heuristic topology from a capability probe. Apple: two shards
    /// (P/E clusters) when both core classes exist, else one. x86/other:
    /// NUMA node count when exposed, else one shard per 16 logical CPUs
    /// (SMT-2 × 8-core CCD), at least one. The result is a SCHEDULING HINT
    /// recorded in logs — never a hardware claim (see CONTRACT no-claims).
    #[must_use]
    pub fn from_probe(p: &CapabilityProbe) -> CcdTopology {
        if p.perf_cores.is_some() && p.eff_cores.is_some() {
            let perf = p.perf_cores.unwrap_or(0).max(1);
            return CcdTopology {
                ccds: 2,
                cores_per_ccd: perf,
            };
        }
        if let Some(nodes) = p.numa_nodes
            && nodes > 0
        {
            return CcdTopology {
                ccds: nodes,
                cores_per_ccd: (p.logical_cpus / nodes.max(1)).max(1),
            };
        }
        let ccds = (p.logical_cpus / 16).max(1);
        CcdTopology {
            ccds,
            cores_per_ccd: (p.logical_cpus / (2 * ccds)).max(1),
        }
    }

    /// The logical-core id range owned by `shard` under this topology
    /// (physical cores, SMT siblings excluded — pinning is fs-exec's job).
    #[must_use]
    pub fn cores_of(&self, shard: u16) -> Range<u32> {
        let s = u32::from(shard).min(self.ccds - 1);
        s * self.cores_per_ccd..(s + 1) * self.cores_per_ccd
    }

    /// Topology from MEASURED L3 groups (fz2.2): ground truth where
    /// [`measured_l3_groups`] finds it, `None` otherwise (caller falls
    /// back to [`CcdTopology::from_probe`] and records which it used).
    #[must_use]
    pub fn from_l3_groups(groups: &[Vec<u32>]) -> Option<CcdTopology> {
        if groups.is_empty() {
            return None;
        }
        let per = groups.iter().map(Vec::len).min().unwrap_or(1).max(1);
        Some(CcdTopology {
            ccds: groups.len() as u32,
            cores_per_ccd: per as u32,
        })
    }
}

/// MEASURED L3 cache groups (fz2.2): the logical CPUs sharing each L3
/// instance — the CCD islands on Zen, cluster caches elsewhere. Ground
/// truth from Linux sysfs (`cache/index3/shared_cpu_list`); empty when
/// the target exposes no such tree (macOS), so callers can distinguish
/// "measured" from "heuristic" and ledger which one they scheduled by.
/// Groups are deduplicated and sorted; SMT siblings are INCLUDED (the
/// caller decides whether to use one or both per physical core).
#[must_use]
pub fn measured_l3_groups() -> Vec<Vec<u32>> {
    let mut groups: Vec<Vec<u32>> = Vec::new();
    let Ok(cpus) = std::fs::read_dir("/sys/devices/system/cpu") else {
        return groups;
    };
    for entry in cpus.flatten() {
        let name = entry.file_name();
        let name = name.to_string_lossy();
        let Some(idx) = name.strip_prefix("cpu") else {
            continue;
        };
        if idx.parse::<u32>().is_err() {
            continue;
        }
        let path = entry.path().join("cache/index3/shared_cpu_list");
        let Ok(list) = std::fs::read_to_string(path) else {
            continue;
        };
        let cpus = parse_cpu_list(list.trim());
        if !cpus.is_empty() && !groups.contains(&cpus) {
            groups.push(cpus);
        }
    }
    groups.sort();
    groups
}

/// The kernel's transparent-hugepage mode ("always" / "madvise" /
/// "never"), read from sysfs — `None` where the tree does not exist
/// (macOS). Ledger this beside any hugepage measurement: advice only
/// does anything under `madvise` (or is redundant under `always`).
#[must_use]
pub fn thp_mode() -> Option<String> {
    let s = std::fs::read_to_string("/sys/kernel/mm/transparent_hugepage/enabled").ok()?;
    let start = s.find('[')?;
    let end = s[start + 1..].find(']')? + start + 1;
    Some(s[start + 1..end].to_string())
}

/// Parse a sysfs cpu list ("0-7,64-71") into sorted ids. Malformed
/// pieces are SKIPPED, not guessed at.
fn parse_cpu_list(s: &str) -> Vec<u32> {
    let mut out = Vec::new();
    for piece in s.split(',') {
        let piece = piece.trim();
        if let Some((a, b)) = piece.split_once('-') {
            if let (Ok(a), Ok(b)) = (a.parse::<u32>(), b.parse::<u32>())
                && a <= b
            {
                out.extend(a..=b);
            }
        } else if let Ok(v) = piece.parse::<u32>() {
            out.push(v);
        }
    }
    out.sort_unstable();
    out.dedup();
    out
}

/// Tile→shard assignment: contiguous, balanced z-order slot ranges (the
/// space-filling-curve partition — each shard's tiles are spatially
/// compact, deterministically).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AffinityMap {
    tile_count: u32,
    /// `bounds[s]..bounds[s+1]` is shard s's z-order slot range.
    bounds: Vec<u32>,
}

impl AffinityMap {
    /// Partition `grid`'s tiles across `topo.ccds` shards: `n / shards`
    /// tiles each, the first `n % shards` shards taking one extra
    /// (deterministic balance within ±1).
    #[must_use]
    pub fn assign(grid: &TileGrid, topo: &CcdTopology) -> AffinityMap {
        let n = grid.tile_count() as u32;
        let shards = topo.ccds.max(1).min(n.max(1));
        let base = n / shards;
        let extra = n % shards;
        let mut bounds = Vec::with_capacity(shards as usize + 1);
        let mut at = 0u32;
        bounds.push(0);
        for s in 0..shards {
            at += base + u32::from(s < extra);
            bounds.push(at);
        }
        AffinityMap {
            tile_count: n,
            bounds,
        }
    }

    /// Number of shards in the map.
    #[must_use]
    pub fn shard_count(&self) -> u16 {
        (self.bounds.len() - 1) as u16
    }

    /// Number of tiles in the map.
    #[must_use]
    pub fn tile_count(&self) -> u32 {
        self.tile_count
    }

    /// The z-order slot range owned by `shard`.
    #[must_use]
    pub fn slots_of(&self, shard: u16) -> Range<u32> {
        let s = shard as usize;
        self.bounds[s]..self.bounds[s + 1]
    }

    /// The shard owning z-order slot `slot` (binary search over bounds).
    #[must_use]
    pub fn shard_of_slot(&self, slot: u32) -> u16 {
        debug_assert!(slot < self.tile_count);
        (self.bounds.partition_point(|&b| b <= slot) - 1) as u16
    }

    /// Canonical JSON affinity table (deterministic; the "affinity tables"
    /// log artifact of the tile-layout bead).
    #[must_use]
    pub fn to_json(&self) -> String {
        let mut s = String::from("{\"shards\":[");
        for shard in 0..self.shard_count() {
            if shard > 0 {
                s.push(',');
            }
            let r = self.slots_of(shard);
            let _ = write!(
                s,
                "{{\"shard\":{shard},\"slot_start\":{},\"slot_end\":{}}}",
                r.start, r.end
            );
        }
        let _ = write!(s, "],\"tiles\":{}}}", self.tile_count);
        s
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tile::TileEdge;

    #[test]
    fn fixture_partitions_are_balanced_and_respect_ccd_boundaries() {
        // 40³ tiles on the Threadripper fixture.
        let g = TileGrid::new([320, 320, 320], TileEdge::E8).expect("grid");
        let topo = CcdTopology::TR_7995WX;
        let map = AffinityMap::assign(&g, &topo);
        assert_eq!(map.shard_count(), 12);
        let n = map.tile_count();
        let sizes: Vec<u32> = (0..map.shard_count())
            .map(|s| map.slots_of(s).len() as u32)
            .collect();
        assert_eq!(
            sizes.iter().sum::<u32>(),
            n,
            "every tile owned exactly once"
        );
        let (min, max) = (sizes.iter().min().unwrap(), sizes.iter().max().unwrap());
        assert!(max - min <= 1, "balance within one tile: {sizes:?}");
        // Worker ranges never straddle a CCD boundary.
        for s in 0..map.shard_count() {
            let cores = topo.cores_of(s);
            assert_eq!(cores.len() as u32, topo.cores_per_ccd);
            assert_eq!(cores.start % topo.cores_per_ccd, 0);
        }
    }

    #[test]
    fn shard_lookup_inverts_the_ranges_and_more_shards_than_tiles_clamps() {
        let g = TileGrid::new([16, 16, 16], TileEdge::E8).expect("grid"); // 8 tiles
        let map = AffinityMap::assign(&g, &CcdTopology::EPYC_128C);
        assert_eq!(map.shard_count(), 8, "shards clamp to tile count");
        for s in 0..map.shard_count() {
            for slot in map.slots_of(s) {
                assert_eq!(map.shard_of_slot(slot), s);
            }
        }
    }

    #[test]
    fn probe_heuristic_yields_a_sane_hint_on_this_machine() {
        let topo = CcdTopology::from_probe(&CapabilityProbe::topology_only());
        assert!(topo.ccds >= 1);
        assert!(topo.cores_per_ccd >= 1);
    }

    #[test]
    fn cpu_list_parsing_and_measured_topology_agree() {
        assert_eq!(parse_cpu_list("0-7,64-71"), {
            let mut v: Vec<u32> = (0..8).collect();
            v.extend(64..72);
            v
        });
        assert_eq!(parse_cpu_list("3"), vec![3]);
        assert_eq!(parse_cpu_list(""), Vec::<u32>::new());
        assert_eq!(
            parse_cpu_list("7-3,junk,5"),
            vec![5],
            "malformed pieces skipped"
        );
        // Measured groups (empty on targets without the sysfs tree) are
        // disjoint, sorted, and convert to a sane topology.
        let groups = measured_l3_groups();
        for w in groups.windows(2) {
            assert!(w[0] < w[1], "groups sorted and deduplicated");
        }
        match CcdTopology::from_l3_groups(&groups) {
            Some(t) => {
                assert_eq!(t.ccds as usize, groups.len());
                assert!(t.cores_per_ccd >= 1);
            }
            None => assert!(groups.is_empty(), "None only when nothing was measured"),
        }
    }

    #[test]
    fn affinity_json_is_deterministic_and_balanced_bookkeeping() {
        let g = TileGrid::new([24, 24, 24], TileEdge::E8).expect("grid"); // 27 tiles
        let map = AffinityMap::assign(&g, &CcdTopology::APPLE_M_CLASS);
        let j = map.to_json();
        assert_eq!(j, map.to_json());
        assert_eq!(
            j,
            "{\"shards\":[{\"shard\":0,\"slot_start\":0,\"slot_end\":14},\
             {\"shard\":1,\"slot_start\":14,\"slot_end\":27}],\"tiles\":27}"
        );
    }
}
