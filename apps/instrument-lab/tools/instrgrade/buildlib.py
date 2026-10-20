import sys, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from instrgrade import *
S = Path(sys.argv[1])
jobs = int(sys.argv[3]) if len(sys.argv) > 3 else 1
lib = build_library(S/'real'/'wav48', S/'grade'/(sys.argv[2] if len(sys.argv)>2 else 'library.json'), jobs=jobs)
from collections import Counter
c = Counter((e['inst'], e['art'], e['dyn']) for e in lib)
for k, v in sorted(c.items()): print(k, v)
bad = [e for e in lib if abs(e['cents_vs_et']) > 40]
print('entries', len(lib), 'with |cents|>40:', len(bad))
for e in bad[:20]: print('  ', e['file'], e['name'], round(e['cents_vs_et'],1), round(e['desc']['f0'],1))
