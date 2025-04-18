# copy_qt_lgpl.py  (drop in project root, run once)
import pathlib, shutil, sys, importlib.util

try:
    import PySide6
except ImportError:
    sys.exit("PySide6 not importable – activate MicroSeq env first.")

lic_dir = pathlib.Path("LICENSES")
lic_dir.mkdir(exist_ok=True)

root = pathlib.Path(PySide6.__file__).parent          # …/site-packages/PySide6

candidates = []

# 1) Wheel or pip layout
candidates += list(root.rglob("*LGPL*"))

# 2) Older conda‑forge layout
candidates += list(
    (root.parent.parent / "share" / "licenses" / "PySide6").glob("*LGPL*")
)

# 3) New conda‑forge layout: into $PREFIX/info/licenses/{PKG}/
info_lic_dir = root.parent.parent.parent / "info" / "licenses"
if info_lic_dir.exists():
    candidates += list(info_lic_dir.glob("*PySide6*LGPL*"))

if not candidates:
    sys.exit(
        "Still couldn't locate the Qt LGPL text. "
        "Download it manually:\n"
        "  curl -L https://www.gnu.org/licenses/lgpl-3.0.txt"
        " -o LICENSES/Qt_LGPLv3.txt"
    )

dst = lic_dir / "Qt_LGPLv3.txt"
shutil.copy(candidates[0], dst)
print("Copied", candidates[0], "→", dst)

