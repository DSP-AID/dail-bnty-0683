"""AAVSO Extended File Format (EFF) v1.2 writer + validator (spec: archive.aavso.org/
aavso-extended-file-format, v1.2, latest update 2021-12-02, retrieved 2026-11-19).

Field order (15): STARID,DATE,MAG,MERR,FILT,TRANS,MTYPE,CNAME,CMAG,KNAME,KMAG,AMASS,
GROUP,CHART,NOTES. Missing values are the literal "na". MAG must contain a dot. Limits:
STARID 30, DATE 16, MAG 8, MERR 6, CNAME 20, CMAG 8, KNAME 20, KMAG 8, AMASS 7, GROUP 5,
CHART 20. DELIM may be any ASCII 32-126 except | # and space (tab/comma spelled out).
Blank lines must be commented with '#'. FILTER codes accepted here: a subset of the
published list (TG/TB/TR/V/B/R/I/CV/CR/SG/SR/SI) — enough for smart-telescope work.
"""
from __future__ import annotations

FIELDS = ["STARID", "DATE", "MAG", "MERR", "FILT", "TRANS", "MTYPE", "CNAME", "CMAG",
          "KNAME", "KMAG", "AMASS", "GROUP", "CHART", "NOTES"]
LIMITS = {"STARID": 30, "DATE": 16, "MAG": 8, "MERR": 6, "CNAME": 20, "CMAG": 8,
          "KNAME": 20, "KMAG": 8, "AMASS": 7, "GROUP": 5, "CHART": 20, "NOTES": 4000}
FILTERS = {"U", "B", "V", "R", "I", "J", "H", "K", "TG", "TB", "TR", "CV", "CR",
           "SZ", "SU", "SG", "SR", "SI", "STU", "STV", "STB", "STY", "STHBW", "STHBN",
           "MA", "MB", "MI", "ZS", "Y", "HA", "HAC"}
REQUIRED_PARAMS = ["TYPE", "OBSCODE", "SOFTWARE", "DELIM", "DATE", "OBSTYPE"]

DELIM_NAME = {"tab": "\t", "comma": ","}


def _delim(spec: str) -> str:
    s = DELIM_NAME.get(str(spec).lower(), str(spec))
    return s if len(s) == 1 else str(spec)


def write_eff(obslist, obscode, software, obstype="CCD", delim=",") -> str:
    """obslist: dicts with keys starid, jd (float), mag, merr, filt, trans, mtype,
    cname, cmag, kname, kmag, amass, group, chart, notes (values pre-formatted or None
    -> 'na'; mag/merr/jd numbers are formatted here)."""
    d = _delim(delim)
    lines = [f"#TYPE=EXTENDED", f"#OBSCODE={obscode}", f"#SOFTWARE={software}",
             f"#DELIM={delim}", "#DATE=JD", f"#OBSTYPE={obstype}",
             "#" + d.join(FIELDS)]
    for o in obslist:
        mag = o["mag"] if isinstance(o["mag"], str) else f"{o['mag']:.3f}"
        merr = o["merr"] if isinstance(o["merr"], str) or o["merr"] is None \
            else f"{o['merr']:.4f}"
        jd = o["jd"] if isinstance(o["jd"], str) else f"{o['jd']:.5f}"
        row = [str(o.get("starid") or "na"), jd, mag, merr or "na",
               o["filt"], o["trans"], o["mtype"],
               str(o.get("cname") or "na"),
               f"{o['cmag']:.3f}" if isinstance(o.get("cmag"), float) else (o.get("cmag") or "na"),
               str(o.get("kname") or "na"),
               f"{o['kmag']:.3f}" if isinstance(o.get("kmag"), float) else (o.get("kmag") or "na"),
               o.get("amass") or "na", o.get("group") or "na", o["chart"], o["notes"]]
        lines.append(d.join(row))
    return "\n".join(lines) + "\n"


def validate_eff(text: str) -> list[str]:
    """Return a list of human-readable spec violations (empty list == valid)."""
    errs: list[str] = []
    lines = text.splitlines()
    params: dict[str, str] = {}
    delim = None
    for ln in lines:
        if ln.startswith("#"):
            body = ln[1:].strip()
            if "=" in body and "," not in body.split("=")[0]:
                k, v = body.split("=", 1)
                params[k.strip().upper()] = v.strip()
    for p in REQUIRED_PARAMS:
        if p not in params or not params[p]:
            errs.append(f"missing/empty required parameter #{p}")
    if "DELIM" in params and params["DELIM"]:
        spec = params["DELIM"]
        d = _delim(spec)
        spelled = spec.lower() in DELIM_NAME       # spec allows "tab"/"comma" spelled out
        if len(d) != 1 or d in "|# " or (not spelled and not (32 <= ord(d) <= 126)):
            errs.append(f"#DELIM invalid: {spec!r}")
        else:
            delim = d
    if params.get("TYPE", "").upper() != "EXTENDED":
        errs.append("#TYPE must be Extended")
    if params.get("DATE", "").upper() not in ("JD", "HJD", "EXCEL"):
        errs.append("#DATE must be JD, HJD or EXCEL")
    if params.get("OBSTYPE", "CCD").upper() not in ("CCD", "DSLR", "PEP"):
        errs.append("#OBSTYPE must be CCD, DSLR or PEP")
    if len(params.get("SOFTWARE", "")) > 255:
        errs.append("#SOFTWARE longer than 255 characters")

    n_data = 0
    for i, ln in enumerate(lines, start=1):
        if not ln.strip():
            errs.append(f"line {i}: uncommented blank line")
            continue
        if ln.startswith("#"):
            continue
        if delim is None:
            errs.append("no usable #DELIM; cannot parse data lines")
            return errs
        f = ln.split(delim)
        if len(f) != len(FIELDS):
            errs.append(f"line {i}: {len(f)} fields, expected {len(FIELDS)}")
            continue
        n_data += 1
        for name, val in zip(FIELDS, f):
            if val.startswith(" ") or val.endswith(" "):
                errs.append(f"line {i}: leading/trailing space in {name}")
            lim = LIMITS.get(name)
            if lim and len(val) > lim:
                errs.append(f"line {i}: {name} exceeds limit {lim} ({val!r})")
        starid, date, mag, merr, filt, trans, mtype = f[:7]
        if starid != "na" and not mag:
            errs.append(f"line {i}: empty STARID")
        try:
            float(date)
        except ValueError:
            errs.append(f"line {i}: DATE {date!r} not numeric for #DATE={params.get('DATE')}")
        m = mag.lstrip("<")
        if "." not in m:
            errs.append(f"line {i}: MAG {mag!r} lacks a decimal point")
        try:
            mv = float(m)
            if not 0 < mv < 30:
                errs.append(f"line {i}: MAG {mag!r} implausible")
        except ValueError:
            errs.append(f"line {i}: MAG {mag!r} not numeric")
        if merr != "na":
            try:
                e = float(merr)
                if not 0 <= e < 10:
                    errs.append(f"line {i}: MERR {merr!r} implausible")
            except ValueError:
                errs.append(f"line {i}: MERR {merr!r} not numeric and not 'na'")
        if filt.upper() not in FILTERS:
            errs.append(f"line {i}: FILTER {filt!r} not an accepted code")
        if trans.upper() not in ("YES", "NO"):
            errs.append(f"line {i}: TRANS {trans!r} must be YES or NO")
        if mtype.upper() not in ("STD", "DIF", "ABS"):
            errs.append(f"line {i}: MTYPE {mtype!r} must be STD, DIF (or legacy ABS)")
    if n_data == 0:
        errs.append("no data rows")
    return errs
