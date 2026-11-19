"""EFF v1.2 writer/validator: clean roundtrip + every documented failure mode caught."""
from seestar_eff import aavso_eff

ROW = dict(starid="000-BFQ-649", jd=2461163.4522, mag=9.230, merr=0.0203, filt="TG",
           trans="NO", mtype="STD", cname="ENSEMBLE", cmag=None,
           kname="G1719707967631113856", kmag=9.94, amass="na", group=None,
           chart="S50_98ac1c17", notes="cube; sidecar-wcs")


def _doc():
    return aavso_eff.write_eff([ROW], "TEST", "seestar_eff test", delim=",")


def test_roundtrip_valid():
    text = _doc()
    assert aavso_eff.validate_eff(text) == []
    assert text.splitlines()[0] == "#TYPE=EXTENDED"
    row = text.splitlines()[7].split(",")
    assert len(row) == 15 and row[2] == "9.230" and row[3] == "0.0203"


def _mutate(fn):
    return aavso_eff.validate_eff(fn(_doc()))


def test_missing_obscode():
    errs = _mutate(lambda t: "\n".join(l for l in t.splitlines() if not l.startswith("#OBSCODE")))
    assert any("#OBSCODE" in e for e in errs)


def test_mag_without_dot_rejected():
    errs = _mutate(lambda t: t.replace(",9.230,", ",9,"))
    assert any("decimal point" in e for e in errs)


def test_bad_filter_trans_mtype_merr_cname_limits():
    assert any("FILTER" in e for e in _mutate(lambda t: t.replace(",TG,", ",XX,")))
    assert any("TRANS" in e for e in _mutate(lambda t: t.replace(",NO,STD", ",MAYBE,STD")))
    assert any("MTYPE" in e for e in _mutate(lambda t: t.replace(",NO,STD", ",NO,DIFX")))
    assert any("MERR" in e for e in _mutate(lambda t: t.replace("0.0203", "-2.0")))
    long = "X" * 25
    assert any("CNAME" in e and "limit" in e for e in _mutate(lambda t: t.replace("ENSEMBLE", long)))


def test_field_count_and_blank_line():
    errs = _mutate(lambda t: t.replace(",S50_98ac1c17,cube; sidecar-wcs", ",S50_98ac1c17"))
    assert any("14 fields" in e for e in errs)
    errs = _mutate(lambda t: t + "\n")
    assert any("blank line" in e for e in errs)


def test_delim_rules():
    tab = aavso_eff.write_eff([ROW], "TEST", "sw", delim="tab")
    assert aavso_eff.validate_eff(tab) == []
    assert any("DELIM" in e for e in _mutate(lambda t: t.replace("#DELIM=,", "#DELIM=|")))
