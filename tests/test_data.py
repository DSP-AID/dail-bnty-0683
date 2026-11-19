"""Clean-machine data path: from tarball+manifest only, ensure_data() rebuilds every member."""
import shutil
from pathlib import Path

from seestar_eff import pipeline

DATA = Path(__file__).resolve().parents[1] / "data"


def test_extract_from_tarball_only(tmp_path):
    man = pipeline.json.load(open(DATA / "manifest.json"))
    for f in (man["tarball"], "manifest.json", "vsx_mwcam.json"):
        shutil.copy(DATA / f, tmp_path / f)
    out = pipeline.ensure_data(tmp_path)
    names = set(pipeline.tarfile.open(out / man["tarball"], "r:gz").getnames())
    assert set(p.name for p in out.iterdir()) == names | {man["tarball"], "manifest.json",
                                                          "vsx_mwcam.json"}
    for name, want in man["members_sha256"].items():
        assert pipeline.sha256(out / name) == want
