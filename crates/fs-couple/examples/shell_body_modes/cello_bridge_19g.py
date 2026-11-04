"""The one bridge change to the cello_specimen.py directory: the bridge volume is set so that the uniform bridge shell weighs
the published typical cello bridge mass. Source: Zhang, Woodhouse and Stoppani, 'Motion of the cello bridge', J. Acoust.
Soc. Am. 140(4) 2016, doi 10.1121/1.4964609, Sec. II p. 2: 'a typical mass of a cello bridge around 19 g'. Density = the
maple pack's 626 kg/m3 (Sonderegger 2013, normal climate). The silhouette area, the old volume and the thickness are read
from the PROVENANCE bridge row (the specimen moves the top's skin under the feet, so they are not the Cannone ratios).
usage: python3 -I cello_bridge_19g.py SPECIMEN_DIR NEW_DIR   (NEW_DIR must not exist)"""
import re
import shutil
import sys
from pathlib import Path

src, dst = Path(sys.argv[1]), Path(sys.argv[2])
MASS_G, RHO = 19.0, 626.0
vol = MASS_G * 1e-3 / RHO * 1e9  # mm3
# read and check everything before writing anything
pv = (src / "part-volumes.tsv").read_text().split("\n")
kv = [i for i, l in enumerate(pv) if l.startswith("bridge\t")]
assert len(kv) == 1
f = pv[kv[0]].split("\t")
old = float(f[1])
pr = (src / "PROVENANCE.tsv").read_text().split("\n")
k = [i for i, l in enumerate(pr) if l.startswith("bridge\t")]
assert len(k) == 1
fields = pr[k[0]].split("\t")  # item, value, source, note: only the value field is parsed and changed
value = fields[1]
m_t = re.findall(r"uniform shell thickness ([0-9.]+) mm", value)
m_a = re.findall(r"area ([0-9]+) mm2", value)
m_v = re.findall(r"volume ([0-9]+) mm3", value)
assert len(m_t) == 1 and len(m_a) == 1 and len(m_v) == 1, (m_t, m_a, m_v)
area, v_row = float(m_a[0]), m_v[0]
assert abs(float(v_row) - old) < 1.0, (v_row, old)  # the row and part-volumes.tsv agree before the change
shutil.copytree(src, dst)  # refuses if dst exists
f[1] = "%.1f" % vol
pv[kv[0]] = "\t".join(f)
(dst / "part-volumes.tsv").write_text("\n".join(pv))
fields[1] = value.replace("uniform shell thickness %s mm" % m_t[0], "uniform shell thickness %.2f mm" % (vol / area)).replace("volume %s mm3" % v_row, "volume %.0f mm3" % vol)
row = "\t".join(fields)
row += ("; SPEC8 (MB 136, as SPEC7 / MB 16000): the volume is NOT the Cannone-ratio estimate (%s mm3, %.1f g at 626 kg/m3) but the published "
        "typical cello bridge mass 19 g (Zhang, Woodhouse and Stoppani, JASA 140(4) 2016, Sec. II: 'a typical mass of a cello bridge around 19 g') "
        "at the maple pack density 626 kg/m3 (Sonderegger 2013): %.0f mm3, uniform %.2f mm; a typical value, not this bridge's; "
        "the real bridge tapers (JOH foot thickness 11) and a uniform shell keeps only its mass" % (v_row, old * RHO * 1e-6, vol, vol / area))
pr[k[0]] = row
(dst / "PROVENANCE.tsv").write_text("\n".join(pr))
print("bridge volume %.1f -> %.1f mm3, uniform thickness %.3f -> %.3f mm (area %.0f mm2), mass at 626: %.2f -> %.2f g" % (
    old, vol, old / area, vol / area, area, old * RHO * 1e-6, vol * RHO * 1e-6))
